module dist_solver_module
  use perf_mod, only: t_startf, t_stopf

  use prec,only:rp
  use iso_c_binding
  use cam_logfile, only: iulog
  use dist_spmv_mod

  include 'netcdf.inc'

  !max nonzeros per row (this does not incl the dense row at the pole)
  integer, parameter :: MAX_NNZ=12
  integer :: output_matrix_count=0

  contains
!-----------------------------------------------------------------------
  subroutine dist_linear_system(mlatd0,mlatd1,mlond0,mlond1, &
    bij,pot_hl,fac_hl,coef_ns,pot)
    ! construct linear system based on bij and coef and solve in pot

    ! note: this is only called by the active processes
    
    ! if FAC is read in, pot_hl is not used, only fac_hl is used
    ! if potential is read in, pot_hl is used, fac_hl is output

    
    use params_module, only:nmlat_h,nmlat_T1,nmlon
    use cons_module, only:read_fac
    use mpi_module, only: mpi_rank,dynamo_world,lat_rank,lon_rank,&
                          nmlat_task,nmlon_task,mygrid_size,&
                          sync_mlat_5d, sync_mlon_5d,& 
                          sync_mlat_3d, sync_mlon_3d, &
                          mlat0, mlat1, mlon0, mlon1, &
                          task_csr_rowstarts, mpi_size, &
                          task_csr_mapping
    
    ! the processor grid only covers one hemisphere ((nmlat_h, nmlon)
    ! nmlat_h => # mag latitudes in one hemisphere
    ! nmlon => num of mag longitudes
    ! nmlat_T1 = 2*nmlat_h-1 # num mag latitudes globally
    
    ! input data here includes halo cells
    integer,intent(in) :: mlatd0,mlatd1,mlond0,mlond1
    real(kind=rp),dimension(mlatd0:mlatd1,mlond0:mlond1),intent(in) :: bij
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: pot_hl
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1),intent(inout) :: fac_hl
    real(kind=rp),dimension(10,2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: coef_ns
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1),intent(out) :: pot
    
    integer,parameter :: root = 0
    integer :: nlonlat,i,j,ic,nnz, nnz_est

    integer,dimension(:), allocatable :: rowptr
    integer,dimension(:), allocatable :: colind
    real(kind=rp),dimension(:), allocatable :: values_csr
    real(kind=rp),dimension(10,mlatd0:mlatd1,mlond0:mlond1) :: coef_s
    real(kind=rp),dimension(10,mlatd0:mlatd1,mlond0:mlond1) :: coef_n

    real(kind=rp),dimension(:), allocatable :: rhs,z,pot_hl_f,sol

    integer :: ierr, fst_row

    type(halo_t) :: halo
    
    !for optional output
    logical :: output_matrix = .false.
    integer :: fileid, dd(4), &
         mygridsize_id, nmlaatt1_idd, nmlon_id, nmlath_id, size_id, sizep1_id, &
         nnz_id, mygridsizep1_id, nprocs_id, nprocsp1_id,  &
         potid, zid, rhsBid, Xid, valuesid, colsid, rowptrid, procrowstartsid, &
         procmappingid, countB(1), startB(1) 
    character(len=200) :: fbuf
    character(kind=c_char,len=:), allocatable :: newfile


    if (output_matrix_count == 0) then
       output_matrix = .true.
       output_matrix_count = output_matrix_count + 1
       write(*,*) 'AB: Output matrix this time through ...'
    else
       output_matrix = .false.
    endif
    
    call t_startf('dist_linear_system')

    ! global matrix size
    nlonlat = nmlat_T1*nmlon 

    !proc domain extents without halo pts (these are global: mlat0, mlon0, mlat1, mlon1)

    !number of grid points I will own (after hemisphere exchange) is mygrid_size (global var)

    !allocate space for rowptr,colind,values_csr
    ! AB why does MAX_NNZ=12? seems like 10 is max?
    allocate(rowptr(mygrid_size+1))
    if (mpi_rank == 0) then ! make room for dense row
       nnz_est = (mygrid_size-1)*MAX_NNZ + (nmlon + 2)
    else
       nnz_est = mygrid_size*MAX_NNZ
    endif
    allocate(colind(nnz_est))
    allocate(values_csr(nnz_est))

    allocate(rhs(mygrid_size))
    allocate(z(mygrid_size))
    allocate(pot_hl_f(mygrid_size))
    allocate(sol(mygrid_size))

    rhs = 0.0
    sol = 0.0
    pot_hl_f = 0.0
    colind = 0
    row_ptr = 0
    valuses_csr = 0.0
    z = 0.0
    
    ! for now, split two hemispheres (keep halo pts)
    do concurrent (i = mlond0:mlond1, j = mlatd0:mlatd1, ic = 1:10)
       coef_s(ic,j,i) = coef_ns(ic,1,j,i) !south
       coef_n(ic,j,i) = coef_ns(ic,2,j,i) !north (note: row = nmlat_T1-j+1) 
    enddo

    ! construct LHS matrix in Block CSR format
    ! Note: rows are contiguous on each proc (for LHS matrix
    ! and RHS) 
    call dist_construct_lhs(nnz_est,bij,coef_s(1:9,:,:),coef_n(1:9,:,:),rowptr,colind,values_csr)
    !need nnz for solver
    nnz = rowptr(mygrid_size+1)-1
    write(*,*) "AB: (after lhs) nnz = ", nnz
    
    ! RHS in Block format to match LHS
    rhs = dist_construct_rhs(coef_s(10,:,:), coef_n(10,:,:))

    !0-based indexing 
    !need matrix to be 0-based index for superlu and the matvec
    colind=colind-1
    rowptr=rowptr-1


    ! determine FAC forcing (dense)
    if (read_fac) then ! input is corrected fac_hl, pot_hl is not used

       z = dist_flatten(fac_hl)

    else ! input is pot_hl, fac_hl is to be calculated (output)

       ! A. Maute 2023/11/21: put the high latitude potential in X
       ! and then use LHS to calculate the RHS FAC
       pot_hl_f = dist_flatten(pot_hl)

       !DEBUG
       do i = 1,mygrid_size
          if (isnan(pot_hl_f(i))) write(*,*) 'AB: pot_hl_f(i) is NaN, i = ', i
       enddo
       
       !Parallel matmult
       ! z = matmul(lhs, pot_hl_f)
       fst_row = task_csr_rowstarts(mpi_rank) ! 0-index
       !write(*,*) "AB: fst_row = ", fst_row
       !write(*,*) "AB: POT_HL_F= ", pot_hl_f

       !TO DO - this init should be called just the first timestep because the nonzero
       !matrix pattern does not change
       call dist_spmv_init(mygrid_size, fst_row, nlonlat, mpi_size, mpi_rank, rowptr, colind, task_csr_rowstarts, dynamo_world, halo, ierr)

       !DEBUG
       !call write_halo_to_file(halo, dynamo_world)

       call dist_spmv(rowptr, colind, values_csr, pot_hl_f, z, halo, dynamo_world, ierr)

       call dist_spmv_free(halo)
       
       ! reconstruct 2D distribution of FAC based on z
       fac_hl(:,mlat0:mlat1,mlon0:mlon1) = dist_unravel(z)

       !get ghost/halo points
       call sync_mlat_3d(fac_hl(:,:,mlon0:mlon1), 2)
       call sync_mlon_3d(fac_hl, 2)
       
     endif !FAC

     ! add FAC forcing to RHS
     !(these are both hemisphere swapped for contiguous rows already)
     do i = 1,mygrid_size
        if (isnan(rhs(i))) write(*,*) 'AB: rhs(i) is NaN, i = ', i
        if (isnan(z(i))) write(*,*) 'AB: z(i) is NaN, i = ', i 
        rhs(i) = rhs(i)+z(i)
     enddo

     !0-based indexing 
     !moved above

     !do we want to output matrix and rhs for debugging
     if (output_matrix == .true. ) then
        write(*,*) 'AB: writing netcdf file ...'

        !newfile = '/glade/derecho/scratch/abaker/dynamo_test.nc'
         write(fbuf,'(A,".",I4.4,".nc")') &
             '/glade/derecho/scratch/abaker/dynamo_test', mpi_rank
         !trim it
         newfile = trim(fbuf)//c_null_char

         write(*,*) 'creating nc-file:', newfile
         ierr = nf_create(newfile, NF_CLOBBER, fileid)

         !Define the dimensions
         ierr = nf_def_dim(fileid,'mygrid_size', mygrid_size, mygridsize_id)
         ierr = nf_def_dim(fileid, 'mygrid_sizep1', mygrid_size+1, mygridsizep1_id)
         ierr = nf_def_dim(fileid, 'mynnz', nnz, nnz_id)
         ierr = nf_def_dim(fileid,'nmlat_T1', nmlat_T1, nmlatt1_id)
         ierr = nf_def_dim(fileid,'nmlat_h', nmlat_h, nmlath_id)
         ierr = nf_def_dim(fileid,'nmlon', nmlon, nmlon_id)
         ierr = nf_def_dim(fileid,'nprocs', mpi_size, nprocs_id)
         ierr = nf_def_dim(fileid,'nprocsp1', mpi_size+1, nprocsp1_id)
         ierr = nf_def_dim(fileid,'global_size', nlonlat, size_id)
         ierr = nf_def_dim(fileid,'global_sizep1', nlonlat + 1, sizep1_id)

         ! define vars and their associated dimensions         
         dd(1) = mygridsize_id
         ierr = nf_def_var(fileid, 'z', NF_DOUBLE, 1, dd, zid)
         dd(1) = mygridsize_id
         ierr = nf_def_var(fileid, 'pot', NF_DOUBLE, 1, dd, potid)
         dd(1) = mygridsize_id
         ierr = nf_def_var(fileid, 'rhsB', NF_DOUBLE, 1, dd, rhsBid)
         dd(1) = mygridsize_id
         ierr = nf_def_var(fileid, 'X', NF_DOUBLE, 1, dd, Xid)
         dd(1) = nnz_id
         ierr = nf_def_var(fileid, 'valuesA', NF_DOUBLE, 1, dd, valuesid)
         dd(1) = nnz_id
         ierr = nf_def_var(fileid, 'colindxA',NF_INT, 1, dd, colsid)
         dd(1) = mygridsizep1_id
         ierr = nf_def_var(fileid, 'rowptrA',NF_INT, 1, dd, rowptrid)
         dd(1) = nprocsp1_id
         ierr = nf_def_var(fileid, 'procRowStarts',NF_INT, 1, dd, procrowstartsid)
         dd(1) = nprocs_id
         ierr = nf_def_var(fileid, 'procMapping',NF_INT, 1, dd, procmappingid)

         !now fill fields 
         !change mode
         ierr=nf_enddef(fileid)
         !Output the  potential
         startB(1)=1
         countB(1)=mygrid_size
         ierr=NF_PUT_VARA_DOUBLE(fileid,potid,startB,countB,pot_hl_f)
         !Output the  z
         startB(1)=1
         countB(1)=mygrid_size
         ierr=NF_PUT_VARA_DOUBLE(fileid,zid,startB,countB,z)
         !Output the  b
         startB(1)=1
         countB(1)=mygrid_size
         ierr=NF_PUT_VARA_DOUBLE(fileid,rhsBid,startB,countB,rhs)
         !Output A in CSR format
         startB(1)=1
         countB(1)=mygrid_size +1
         ierr=NF_PUT_VARA_INT(fileid,rowptrid,startB,countB,rowptr)
         startB(1)=1
         countB(1)=nnz
         ierr=NF_PUT_VARA_INT(fileid,colsid,startB,countB,colind)
         startB(1)=1
         countB(1)=nnz
         ierr=NF_PUT_VARA_DOUBLE(fileid,valuesid,startB,countB,values_csr)

         startB(1)=1
         countB(1)=mpi_size+1
         ierr=NF_PUT_VARA_INT(fileid,procrowstartsid,startB,countB,task_csr_rowstarts)
         startB(1)=1
         countB(1)=mpi_size
         ierr=NF_PUT_VARA_INT(fileid,procmappingid,startB,countB,task_csr_mapping)


         !Close the file up (later after solve)
         !ierr=NF_CLOSE(fileid)

      endif
     
     call t_startf('linear_system->solve_superlu')
     sol = dist_solve_superlu(nlonlat, mygrid_size, nnz, rowptr, colind(1:nnz), values_csr(1:nnz), rhs)
     call t_stopf('linear_system->solve_superlu')

     !if output turned on for debugging 
     if (output_matrix == .true.) then
        startB(1)=1
        countB(1)=mygrid_size
        ierr=NF_PUT_VARA_DOUBLE(fileid,Xid,startB,countB,rhs)
        !close file
        ierr=NF_CLOSE(fileid)
     endif
#if 0     
     ! reconstruct 2D distribution of potential based on the solution
     pot(1:2,mlat0:mlat1,mlon0:mlon1) = dist_unravel(sol)
     !get ghost/halo points
     call sync_mlat_3d(pot(:,:,mlon0:mlon1), 2)
     call sync_mlon_3d(pot, 2)
#endif

     !pot and fac_hl are ready to return (updated grid + halo)

     ! Deallocate arrays to free memory
     if (allocated(rowptr)) deallocate(rowptr)
     if (allocated(colind)) deallocate(colind) 
     if (allocated(values_csr)) deallocate(values_csr)
     if (allocated(rhs)) deallocate(rhs)
     if (allocated(z)) deallocate(z)
     if (allocated(pot_hl_f)) deallocate(pot_hl_f)
     if (allocated(sol)) deallocate(sol)
      
     call t_stopf('dist_linear_system')

  endsubroutine dist_linear_system
!-----------------------------------------------------------------------
  subroutine dist_construct_lhs(nnz_est,bij,coef_s, coef_n,my_rowptr,my_colind, my_values)
! construct LHS matrix (CSR format - block row format for each task)

! need to set where the two hemispheres are connected
! This is not in the 9-point stencil but needs to be done manually

    
    use params_module,only:nmlat_h,nmlat_T1,nmlon
    use cons_module,only:jlatm_JT
    use mpi_module,only:mpi_rank,mpi_size,dynamo_world,lat_rank,lon_rank, &
         nmlat_task,nmlon_task,mlatd0,mlatd1,mlat0,mlat1, &
         mlond0,mlond1,mlon0,mlon1, &
         lat_size,lon_size,task_lat_offset,ij_start_n,ij_stop_n, &
         ij_start_s,ij_stop_s,partner_hgridsize,my_hgridsize, &
         calc_grid_ij, partner_exchange_hemisphere_mat, &
         gather_lon_1d, mygrid_size, mpi_partner, &
         mygrid_size_n, mygrid_size_s
    
    integer,intent(in) :: nnz_est
    real(kind=rp),dimension(mlatd0:mlatd1,mlond0:mlond1),intent(in) :: bij
    real(kind=rp),dimension(9,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: coef_s
    real(kind=rp),dimension(9,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: coef_n
    integer,dimension(mygrid_size+1),intent(out) :: my_rowptr
    integer,dimension(nnz_est),intent(out) :: my_colind
    real(kind=rp),dimension(nnz_est),intent(out) :: my_values

!mlatd0 and mlond0 can be 0    
!mlatd1 and mlond1 can be nmlat+1 and nmlon+1

    
! if two hemispheres are uncoupled at high latitudes, set bijSum to zero
    real(kind=rp),parameter :: bijSum = 0
    integer :: nlonlat,i,j,jS,jN,ij,isub,im,ip, k, cnt
    integer :: loop_start_i, loop_stop_i, loop_start_j, loop_stop_j
    integer :: nnz_s, nnz_n, row_counter_s, row_counter_n
    integer :: nnz_south, nnz_north, num_row_n, num_row_s
    
    ! the first row has most elements (nmlon+2)
    !special case (only 1 processor does this)
    integer,dimension(nmlon+2) :: jcol1
    real(kind=rp),dimension(nmlon+2) :: nzval1
    integer :: counter


    ! other rows have at most MAX_NNZ elements
    !we want the indexes to be the global ranges of rows....
    integer,dimension(MAX_NNZ,ij_start_s:ij_stop_s) :: jcol_s
    real(kind=rp),dimension(MAX_NNZ,ij_start_s:ij_stop_s) :: nzval_s
    
    integer,dimension(MAX_NNZ,ij_start_n:ij_stop_n) :: jcol_n
    real(kind=rp),dimension(MAX_NNZ,ij_start_n:ij_stop_n) :: nzval_n 

    integer,dimension(ij_start_s:ij_stop_s) :: rowcnt_s
    integer,dimension(ij_start_n:ij_stop_n) :: rowcnt_n
    
    real(kind=rp),dimension(nmlon) :: coef3_j1_buf
    real(kind=rp),dimension(nmlon) :: coef9_j1_buf

    !csr for each hemisphere - south
    integer, dimension(mygrid_size_s+1):: rowptr_s
    integer, dimension((mygrid_size_s)*MAX_NNZ + nmlon) ::  colind_s
    real(kind=rp), dimension((mygrid_size_s)*MAX_NNZ + nmlon) :: values_s

    !csr for each hemisphere - north
    integer, dimension(mygrid_size_n+1):: rowptr_n
    integer, dimension((mygrid_size_n)*MAX_NNZ) ::  colind_n
    real(kind=rp), dimension((mygrid_size_n)*MAX_NNZ) :: values_n
    
    !get hemisphere partner info
    integer, dimension(partner_hgridsize+1) :: partner_rowptr
    real(kind=rp),dimension(partner_hgridsize*MAX_NNZ) :: partner_values
    integer, dimension(partner_hgridsize*MAX_NNZ) :: partner_cols

    !global size 
    nlonlat = nmlat_T1*nmlon

    !initialize
    rowcnt_s = 0
    rowcnt_n = 0

    rowptr_s= 0
    colind_s = 0
    values_s = 0.0

    rowptr_n = 0
    colind_n = 0
    values_n = 0.0
    
    
    jcol1 = 0
    nzval1 = 0.0
    coef3_j1_buf = 0.0

    !   some debugging info
    !write(*, *) 'AB: LHS: mpi_rank, ij_start_s, ij_stop_s, s_grid_pts = ', mpi_rank, ij_start_s, ij_stop_s, ij_stop_s - ij_start_s + 1
    !write(*, *) 'AB: LHS: mpi_rank, ij_start_n, ij_stop_n, n_grid_pts = ',mpi_rank, ij_start_n, ij_stop_n,  ij_stop_n - ij_start_n + 1
    !write(*, *) 'AB: LHS: mpi_rank, my_hgridsize, mpi_partner, partner_hgridsize, mygrid_size = ',mpi_rank, my_hgridsize, mpi_partner, partner_hgridsize, mygrid_size
    !write(*, *) 'AB: LHS: mpi_rank, lon_rank, lat_rank',mpi_rank, lon_rank, lat_rank

    
    
!_______________    

    ! i, j, ij are all global indices - each task has a subset of the grid
    ! i loops thru longitude, j thru latitude, and ij is the global row in the matrix

    ! recall the proc grid looks like this:
    ! stack along latitudes first then longitudes
    ! (lat_size=3)
    !  8  9 10 11
    !  4  5  6  7
    !  0  1  2  3 (lon_size=4)

    ! within each proc grid, the row numbering looks like this  (assume each has a 2x2 grid):
    ! (each grid point is a row in the marix)
    !task 0:
    ! 2 4
    ! 1 3
    !task 1:
    ! 6 8
    ! 5 7
    !task 4:
    ! 18 20
    ! 17 19
    ! etc.

    ! so now calculating ij is bit more involved (compared to serial) and we need info
    ! about the processor grid

    ! procs will set up rows for their domain in the southern hemisphere and northern hemisphere
    ! we'll keep in two different arrays to combine when we assemble the CSR matricx
    ! each proc has to have a continguous block od rows, so we'll have to do a bit of
    ! rearraging
    
    ! Set up matrix with loops through the following regions:
    ! poles, lower/upper lats (through latm_JT), then latm_Jt, through equator

    !coefficients:
    ! ^   equatorward
! coef(4) (i-1,j+1)      coef(3) (i,j+1)    coef(2) (i+1,j+1)
! coef(5) (i-1,j)        coef(9) (i,j)      coef(1) (i+1,j)
! coef(6) (i-1,j-1)      coef(7) (i,j-1)    coef(8) (i+1,j-1)
! v   poleward

    
    !Outer if statements to see if a proc owns rows in that region

    !POLES------------------------------------------_!
    if (lat_rank == 0) then ! I own the pole regions (j=1)
       j = 1
       ! there are no c6,c7,c8 values at the south pole

       !do communication for coef(3,1,i) here so that root proc can complete
       !it's row below
       !gather_lon_1d()

       !the proc that owns j=1,i=1 (proc 0)  needs  coef(3,j=1,isub)
       !for isub=2:nmlon, but only owns 2:mlon1
       !so we need to communicate within proc lat_rank = 0 and gather to
       ! proc 0 (check scalability here at large proc counts)
       coef3_j1_buf = gather_lon_1d(coef_s(3,j,mlon0:mlon1))
       !also needs coef9  for the whole row corres to j=1
       coef9_j1_buf = gather_lon_1d(coef_s(9,j,mlon0:mlon1))

       counter = 0
       if (lon_rank == 0) then !I also own i=1 (special case - applies to 1 processor)
          i=1
          counter = counter + 1
          jcol1(counter) = calc_grid_ij(i,j,0) !this will be 1
          !this requires the whole row (needs to get the rest)
          !nzval1(counter) = sum(coef_s(9,j,:))-bijSum
          nzval1(counter) = sum(coef9_j1_buf(:))-bijSum
          
          counter = counter + 1
          jcol1(counter) = calc_grid_ij(i,j+1,0) 
          nzval1(counter) = coef_s(3,j,i)

          counter = counter + 1
          jcol1(counter) = calc_grid_ij(i,nmlat_T1,0)    
          nzval1(counter) = bijSum

          ! this processor only owns the coef(3,1,i) for i <= mlond1!
          ! (needs to get the rest from procs in this row! done above)
          do isub = 2,nmlon
             ! Sum_i=1^nmlon C3(i,1) Phi(i,2)
             ! lhs(ij,(isub-1)*nmlat_T1+j+1) = coef(3,j,isub)
             !jcol1(isub+2) = (isub-1)*nmlat_T1+j+1
             !nzval1(isub+2) = coef(3,j,isub)
             counter = counter +1
             jcol1(counter) =  calc_grid_ij(isub,j+1,0)
             nzval1(counter) = coef3_j1_buf(isub)
          enddo

          rowcnt_s(1) = counter !should be = nmlon+2
          
          !jcol1 will be sorted already
          
          !get ready for next loop
          loop_start_i = 2
       else !I don't own i=1
          loop_start_i = mlon0
       endif ! end i=1

       !loop_start_i is set in above if statement (based on owning i=1 or not)
       loop_stop_i = mlon1

       !South pole - other longitudes (i/=1)
       do i = loop_start_i, loop_stop_i   
          !calculate matrix row for South
          jS=j
          ij = calc_grid_ij(i,jS,lat_rank)
          
          !connection to pole (jS=1)
          jcol_s(rowcnt_s(ij)+1,ij) = jS  !this will be 1
          nzval_s(rowcnt_s(ij)+1,ij) = -1
          !center (ij)
          jcol_s(rowcnt_s(ij)+2,ij) = ij
          nzval_s(rowcnt_s(ij)+2,ij) = 1
          !increase row count for ij
          rowcnt_s(ij) = rowcnt_s(ij)+2

          !sort by col indices
          call insert_sort(jcol_s(:,ij),nzval_s(:,ij),rowcnt_s(ij))
          
       enddo
       
       ! North pole: for each i, set Phi(i,nmlat_T1) = phi_pol
       loop_start_i = mlon0
       do i = loop_start_i, loop_stop_i   
          !calculate matrix row for North
          jN = nmlat_T1-j+1 !j=1, so jN= nmlat_T1
          ij = calc_grid_ij(i,jN,lat_rank)

          rowcnt_n(ij) = rowcnt_n(ij)+1
          ! TO DO: verify this (doesn't make intuive sense to me)
          jcol_n(rowcnt_n(ij),ij) = calc_grid_ij(i,jN,0) ! j=jN, lat_rank=0 owns
          nzval_n(rowcnt_n(ij),ij) = 1

          !sort by col indices
          call insert_sort(jcol_n(:,ij),nzval_n(:,ij),rowcnt_n(ij))
        enddo
       
    endif !end of loop for procs owning j=1

    
    !REGION BETWEEN POLES AND jlatm_JT (i.e., 2:jlatm_JT-1)
    if (mlat0<jlatm_JT) then !I own latitudes in this region

       !loop through the longitudes in my grid
       do i = mlon0, mlon1

          !loop through latitudes
          if (lat_rank == 0) then !skip the pole
             loop_start_j = 2
          else
             loop_start_j = mlat0
          endif
          loop_stop_j = min(mlat1, jlatm_JT-1)
          
          do j = loop_start_j, loop_stop_j            
             !!!!!!!South Hemishere
             jS=j
             jN = nmlat_T1-j+1
             ij = calc_grid_ij(i,jS,lat_rank)

             if (ij > ij_stop_s .or. ij < ij_start_s) then
                write(*,*) 'AB: Error ij south index 1 for lhs', mpi_rank, ij
             endif

             !the edges of the poc domain are handles in calc_grid_ij
             
             !coef 6 (i-1, j-1)
             rowcnt_s(ij) = rowcnt_s(ij)+1
             jcol_s(rowcnt_s(ij),ij)= calc_grid_ij(i-1, jS-1, lat_rank)
             nzval_s(rowcnt_s(ij),ij)= coef_s(6,j,i)
             !coef 5 (i-1, j)
             rowcnt_s(ij) = rowcnt_s(ij)+1
             jcol_s(rowcnt_s(ij),ij)= calc_grid_ij(i-1, jS, lat_rank)
             nzval_s(rowcnt_s(ij),ij)= coef_s(5,j,i)
             !coef 4 (i-1, j+1)
             rowcnt_s(ij) = rowcnt_s(ij)+1
             jcol_s(rowcnt_s(ij),ij)= calc_grid_ij(i-1, jS+1, lat_rank)
             nzval_s(rowcnt_s(ij),ij)= coef_s(4,j,i)

             !coef 7 (i, j-1)
             rowcnt_s(ij) = rowcnt_s(ij)+1
             jcol_s(rowcnt_s(ij),ij)= calc_grid_ij(i, jS-1, lat_rank)
             nzval_s(rowcnt_s(ij),ij)= coef_s(7,j,i)
             !coef 9 (i, j) ! should be 1
             rowcnt_s(ij) = rowcnt_s(ij)+1
             jcol_s(rowcnt_s(ij),ij)= calc_grid_ij(i, jS, lat_rank)
             nzval_s(rowcnt_s(ij),ij)= coef_s(9,j,i)-bij(j,i)
             !coef 3 (i, j+1)
             rowcnt_s(ij) = rowcnt_s(ij)+1
             jcol_s(rowcnt_s(ij),ij)= calc_grid_ij(i, jS+1, lat_rank)
             nzval_s(rowcnt_s(ij),ij)= coef_s(3,j,i)

             !bij -> connects to jN (conjugate point)
             rowcnt_s(ij) = rowcnt_s(ij)+1
             jcol_s(rowcnt_s(ij),ij)= calc_grid_ij(i, jN, lat_rank)
             nzval_s(rowcnt_s(ij),ij)= bij(j,i)
                
             !coef 8 (i+1, j-1)
             rowcnt_s(ij) = rowcnt_s(ij)+1
             jcol_s(rowcnt_s(ij),ij)= calc_grid_ij(i+1, jS-1, lat_rank)
             nzval_s(rowcnt_s(ij),ij)= coef_s(8,j,i)
             !coef 1 (i+1, j)
             rowcnt_s(ij) = rowcnt_s(ij)+1
             jcol_s(rowcnt_s(ij),ij)= calc_grid_ij(i+1, jS, lat_rank)
             nzval_s(rowcnt_s(ij),ij)= coef_s(1,j,i)
             !coef 2 (i+1, j+1)
             rowcnt_s(ij) = rowcnt_s(ij)+1
             jcol_s(rowcnt_s(ij),ij)= calc_grid_ij(i+1, jS+1, lat_rank)
             nzval_s(rowcnt_s(ij),ij)= coef_s(2,j,i)

             !sort by col indices
             call insert_sort(jcol_s(:,ij),nzval_s(:,ij),rowcnt_s(ij))
          
             
             !!!!!!!!!!North Hemisphere
             ! note that coef has the direction switched in the NH
             ij = calc_grid_ij(i,jN,lat_rank)
             if (ij > ij_stop_n .or. ij < ij_start_n) then
                write(iulog,*) 'AB: Error ij north index 1 for lhs', mpi_rank, ij
             endif
             
             !coef 4 (i-1, j-1)
             rowcnt_n(ij) = rowcnt_n(ij)+1
             jcol_n(rowcnt_n(ij),ij)= calc_grid_ij(i-1, jN-1, lat_rank)
             nzval_n(rowcnt_n(ij),ij)= coef_n(4,j,i)
             !coef 5 (i-1, j)
             rowcnt_n(ij) = rowcnt_n(ij)+1
             jcol_n(rowcnt_n(ij),ij)= calc_grid_ij(i-1, jN, lat_rank)
             nzval_n(rowcnt_n(ij),ij)= coef_n(5,j,i)
             !coef 6 (i-1, j+1)
             rowcnt_n(ij) = rowcnt_n(ij)+1
             jcol_n(rowcnt_n(ij),ij)= calc_grid_ij(i-1, jN+1, lat_rank)
             nzval_n(rowcnt_n(ij),ij)= coef_n(6,j,i)

             !coef 3 (i, j-1)
             rowcnt_n(ij) = rowcnt_n(ij)+1
             jcol_n(rowcnt_n(ij),ij)= calc_grid_ij(i, jN-1, lat_rank)
             nzval_n(rowcnt_n(ij),ij)= coef_n(3,j,i)
             !coef 9 (i, j)
             rowcnt_n(ij) = rowcnt_n(ij)+1
             jcol_n(rowcnt_n(ij),ij)= calc_grid_ij(i, jN, lat_rank)
             nzval_n(rowcnt_n(ij),ij)= coef_n(9,j,i) - bij(j,i)
             !coef 7 (i, j+1)
             rowcnt_n(ij) = rowcnt_n(ij)+1
             jcol_n(rowcnt_n(ij),ij)= calc_grid_ij(i, jN+1, lat_rank)
             nzval_n(rowcnt_n(ij),ij)= coef_n(7,j,i)
             !bij -> connects to jS (conjugate point)
             rowcnt_n(ij) = rowcnt_n(ij)+1
             jcol_n(rowcnt_n(ij),ij)= calc_grid_ij(i, jS, lat_rank)
             nzval_n(rowcnt_n(ij),ij)= bij(j,i)             
             
             !coef 2 (i+1, j-1)
             rowcnt_n(ij) = rowcnt_n(ij)+1
             jcol_n(rowcnt_n(ij),ij)= calc_grid_ij(i+1, jN-1, lat_rank)
             nzval_n(rowcnt_n(ij),ij)= coef_n(2,j,i)
             !coef 1 (i+1, j)
             rowcnt_n(ij) = rowcnt_n(ij)+1
             jcol_n(rowcnt_n(ij),ij)= calc_grid_ij(i+1, jN, lat_rank)
             nzval_n(rowcnt_n(ij),ij)= coef_n(1,j,i)
             !coef 8 (i+1, j+1)
             rowcnt_n(ij) = rowcnt_n(ij)+1
             jcol_n(rowcnt_n(ij),ij)= calc_grid_ij(i+1, jN+1, lat_rank)
             nzval_n(rowcnt_n(ij),ij)= coef_n(8,j,i)

             !sort by col indices
             call insert_sort(jcol_n(:,ij),nzval_n(:,ij),rowcnt_n(ij))
             
          enddo ! end j loop through latitudes
       enddo !end i loop through longitudes
    endif ! end of REGION BETWEEN POLES AND jlatm_JT

    !---------------------------------------!   
          
    !jlatm_JT REGION- two hemispheres are coupled at j-1
    if (mlat0 <= jlatm_JT .and. mlat1 >= jlatm_JT) then ! I own grid points at latm_JT

       !loop through the longitudes in my grid
       do i = mlon0, mlon1

          !just one latitude (so no loop)
          j = jlatm_JT
          jS = j
          jN = nmlat_T1-j+1

          !Southern hemisphere
          ij = calc_grid_ij(i,jS,lat_rank)
          if (ij > ij_stop_s .or. ij < ij_start_s) then
             write(*,*) 'AB: Error ij south index 2 for lhs', mpi_rank, ij
          endif
          
          !coef 6 (i-1, j-1)
          rowcnt_s(ij) = rowcnt_s(ij)+1
          jcol_s(rowcnt_s(ij),ij)= calc_grid_ij(i-1, jS-1, lat_rank)
          nzval_s(rowcnt_s(ij),ij)= coef_s(6,j,i)
          !coef 5 (i-1, j)
          rowcnt_s(ij) = rowcnt_s(ij)+1
          jcol_s(rowcnt_s(ij),ij)= calc_grid_ij(i-1, jS, lat_rank)
          nzval_s(rowcnt_s(ij),ij)= coef_s(5,j,i)
          !coef 4 (i-1, j+1)
          rowcnt_s(ij) = rowcnt_s(ij)+1
          jcol_s(rowcnt_s(ij),ij)= calc_grid_ij(i-1, jS+1, lat_rank)
          nzval_s(rowcnt_s(ij),ij)= coef_s(4,j,i)
          !extra connection to north at 6
          !      lhs(ij,(im-1)*nmlat_T1+jN+1) = coef(6,jN,i)
          rowcnt_s(ij) = rowcnt_s(ij)+1
          !!!!AB - double check
          jcol_s(rowcnt_s(ij),ij)= calc_grid_ij(i-1, jN+1, lat_rank)
          nzval_s(rowcnt_s(ij),ij)= coef_n(6,j,i)

          !coef 7 (i, j-1)
          rowcnt_s(ij) = rowcnt_s(ij)+1
          jcol_s(rowcnt_s(ij),ij)= calc_grid_ij(i, jS-1, lat_rank)
          nzval_s(rowcnt_s(ij),ij)= coef_s(7,j,i)
          !coef 9 (i, j) ! should be 1
          rowcnt_s(ij) = rowcnt_s(ij)+1
          jcol_s(rowcnt_s(ij),ij)= calc_grid_ij(i, jS, lat_rank)
          nzval_s(rowcnt_s(ij),ij)= coef_s(9,j,i)
          !coef 3 (i, j+1)
          rowcnt_s(ij) = rowcnt_s(ij)+1
          jcol_s(rowcnt_s(ij),ij)= calc_grid_ij(i, jS+1, lat_rank)
          nzval_s(rowcnt_s(ij),ij)= coef_s(3,j,i)
          !extra connection to north at 7
          !      lhs(ij, (i-1)*nmlat_T1+jN+1) = coef(7,jN,i)
          rowcnt_s(ij) = rowcnt_s(ij)+1
          !!!!!AB- double check
          jcol_s(rowcnt_s(ij),ij)= calc_grid_ij(i, jN+1, lat_rank)
          nzval_s(rowcnt_s(ij),ij)= coef_n(7,j,i)

          !coef 8 (i+1, j-1)
          rowcnt_s(ij) = rowcnt_s(ij)+1
          jcol_s(rowcnt_s(ij),ij)= calc_grid_ij(i+1, jS-1, lat_rank)
          nzval_s(rowcnt_s(ij),ij)= coef_s(8,j,i)
          !coef 1 (i+1, j)
          rowcnt_s(ij) = rowcnt_s(ij)+1
          jcol_s(rowcnt_s(ij),ij)= calc_grid_ij(i+1, jS, lat_rank)
          nzval_s(rowcnt_s(ij),ij)= coef_s(1,j,i)
          !coef 2 (i+1, j+1)
          rowcnt_s(ij) = rowcnt_s(ij)+1
          jcol_s(rowcnt_s(ij),ij)= calc_grid_ij(i+1, jS+1, lat_rank)
          nzval_s(rowcnt_s(ij),ij)= coef_s(2,j,i)
          !extra connection to north at 8
          !      lhs(ij,(ip-1)*nmlat_T1+jN+1) = coef(8,jN,i)
          rowcnt_s(ij) = rowcnt_s(ij)+1
          !!!!!AB double -check
          jcol_s(rowcnt_s(ij),ij)= calc_grid_ij(i+1, jN+1, lat_rank)
          nzval_s(rowcnt_s(ij),ij)= coef_n(8,j,i)

          !sort by col indices
          call insert_sort(jcol_s(:,ij),nzval_s(:,ij),rowcnt_s(ij))
          
          
          !!Northern Hemi

          !AB-calc_ij for North hemisphere
          ij = calc_grid_ij(i,jN,lat_rank)
          if (ij > ij_stop_n .or. ij < ij_start_n) then
             write(*,*) 'AB: Error ij north index 2 for lhs', mpi_rank, ij
          endif
             
          !extra connection to South at 6
          !      lhs(ij,(im-1)*nmlat_T1+jS-1) = coef(6,jS,i)
          rowcnt_n(ij) = rowcnt_n(ij)+1
          jcol_n(rowcnt_n(ij),ij)= calc_grid_ij(i-1, jS-1, lat_rank)
          nzval_n(rowcnt_n(ij),ij)= coef_s(6,j,i)
          !coef 4 (i-1, j-1)
          rowcnt_n(ij) = rowcnt_n(ij)+1
          jcol_n(rowcnt_n(ij),ij)= calc_grid_ij(i-1, jN-1, lat_rank)
          nzval_n(rowcnt_n(ij),ij)= coef_n(4,j,i)
          !coef 5 (i-1, j)
          rowcnt_n(ij) = rowcnt_n(ij)+1
          jcol_n(rowcnt_n(ij),ij)= calc_grid_ij(i-1, jN, lat_rank)
          nzval_n(rowcnt_n(ij),ij)= coef_n(5,j,i)
          !coef 6 (i-1, j+1)
          rowcnt_n(ij) = rowcnt_n(ij)+1
          jcol_n(rowcnt_n(ij),ij)= calc_grid_ij(i-1, jN+1, lat_rank)
          nzval_n(rowcnt_n(ij),ij)= coef_n(6,j,i)

          !extra connection to South at 7
          !      lhs(ij, (i-1)*nmlat_T1+jS-1) = coef(7,jS,i)
          rowcnt_n(ij) = rowcnt_n(ij)+1
          jcol_n(rowcnt_n(ij),ij)= calc_grid_ij(i, jS-1, lat_rank)
          nzval_n(rowcnt_n(ij),ij)= coef_s(7,j,i)
          !coef 3 (i, j-1)
          rowcnt_n(ij) = rowcnt_n(ij)+1
          jcol_n(rowcnt_n(ij),ij)= calc_grid_ij(i, jN-1, lat_rank)
          nzval_n(rowcnt_n(ij),ij)= coef_n(3,j,i)
          !coef 9 (i, j)
          rowcnt_n(ij) = rowcnt_n(ij)+1
          jcol_n(rowcnt_n(ij),ij)= calc_grid_ij(i, jN, lat_rank)
          nzval_n(rowcnt_n(ij),ij)= coef_n(9,j,i) 
          !coef 7 (i, j+1)
          rowcnt_n(ij) = rowcnt_n(ij)+1
          jcol_n(rowcnt_n(ij),ij)= calc_grid_ij(i, jN+1, lat_rank)
          nzval_n(rowcnt_n(ij),ij)= coef_n(7,j,i)

          !extra connection to South at 8
          !      lhs(ij,(ip-1)*nmlat_T1+jS-1) = coef(8,jS,i)
          rowcnt_n(ij) = rowcnt_n(ij)+1
          jcol_n(rowcnt_n(ij),ij)= calc_grid_ij(i+1, jS-1, lat_rank)
          nzval_n(rowcnt_n(ij),ij)= coef_s(8,j,i)
          !coef 2 (i+1, j-1)
          rowcnt_n(ij) = rowcnt_n(ij)+1
          jcol_n(rowcnt_n(ij),ij)= calc_grid_ij(i+1, jN-1, lat_rank)
          nzval_n(rowcnt_n(ij),ij)= coef_n(2,j,i)
          !coef 1 (i+1, j)
          rowcnt_n(ij) = rowcnt_n(ij)+1
          jcol_n(rowcnt_n(ij),ij)= calc_grid_ij(i+1, jN, lat_rank)
          nzval_n(rowcnt_n(ij),ij)= coef_n(1,j,i)
          !coef 8 (i+1, j+1)
          rowcnt_n(ij) = rowcnt_n(ij)+1
          jcol_n(rowcnt_n(ij),ij)= calc_grid_ij(i+1, jN+1, lat_rank)
          nzval_n(rowcnt_n(ij),ij)= coef_n(8,j,i)

          !sort by col indices
          call insert_sort(jcol_n(:,ij),nzval_n(:,ij),rowcnt_n(ij))
          
          
       enddo !end of loop through longitude
       !only 1 latitude
       
    endif
    
    ! between latm_JT and equator (nmlat_h) REGION (symmetric solution)
    ! don't include equator
    if ((mlat0 > jlatm_JT .and. mlat0 < nmlat_h ) .or. (mlat1 > jlatm_JT .and. mlat0 < nmlat_h)) then

        !loop through the longitudes in my grid
       do i = mlon0, mlon1

          !loop through relevant latitudes (don't do equator)
          loop_start_j = max(mlat0, jlatm_JT+1)
          loop_stop_j = min(mlat1, nmlat_h-1)
          
          do j = loop_start_j, loop_stop_j            
             !!!!!!!South Hemishere
             jS=j
             jN = nmlat_T1-j+1
             ij = calc_grid_ij(i,jS,lat_rank)

             if (ij > ij_stop_s .or. ij < ij_start_s) then
                write(*,*) 'AB: Error ij south index 3 for lhs', mpi_rank, ij
             endif
             
             !coef 6 (i-1, j-1)
             rowcnt_s(ij) = rowcnt_s(ij)+1
             jcol_s(rowcnt_s(ij),ij)= calc_grid_ij(i-1, jS-1, lat_rank)
             nzval_s(rowcnt_s(ij),ij)= coef_s(6,j,i)
             !coef 5 (i-1, j)
             rowcnt_s(ij) = rowcnt_s(ij)+1
             jcol_s(rowcnt_s(ij),ij)= calc_grid_ij(i-1, jS, lat_rank)
             nzval_s(rowcnt_s(ij),ij)= coef_s(5,j,i)
             !coef 4 (i-1, j+1)
             rowcnt_s(ij) = rowcnt_s(ij)+1
             jcol_s(rowcnt_s(ij),ij)= calc_grid_ij(i-1, jS+1, lat_rank)
             nzval_s(rowcnt_s(ij),ij)= coef_s(4,j,i)

             !coef 7 (i, j-1)
             rowcnt_s(ij) = rowcnt_s(ij)+1
             jcol_s(rowcnt_s(ij),ij)= calc_grid_ij(i, jS-1, lat_rank)
             nzval_s(rowcnt_s(ij),ij)= coef_s(7,j,i)
             !coef 9 (i, j) ! should be 1
             rowcnt_s(ij) = rowcnt_s(ij)+1
             jcol_s(rowcnt_s(ij),ij)= calc_grid_ij(i, jS, lat_rank)
             nzval_s(rowcnt_s(ij),ij)= coef_s(9,j,i)
             !coef 3 (i, j+1)
             rowcnt_s(ij) = rowcnt_s(ij)+1
             jcol_s(rowcnt_s(ij),ij)= calc_grid_ij(i, jS+1, lat_rank)
             nzval_s(rowcnt_s(ij),ij)= coef_s(3,j,i)

             !coef 8 (i+1, j-1)
             rowcnt_s(ij) = rowcnt_s(ij)+1
             jcol_s(rowcnt_s(ij),ij)= calc_grid_ij(i+1, jS-1, lat_rank)
             nzval_s(rowcnt_s(ij),ij)= coef_s(8,j,i)
             !coef 1 (i+1, j)
             rowcnt_s(ij) = rowcnt_s(ij)+1
             jcol_s(rowcnt_s(ij),ij)= calc_grid_ij(i+1, jS, lat_rank)
             nzval_s(rowcnt_s(ij),ij)= coef_s(1,j,i)
             !coef 2 (i+1, j+1)
             rowcnt_s(ij) = rowcnt_s(ij)+1
             jcol_s(rowcnt_s(ij),ij)= calc_grid_ij(i+1, jS+1, lat_rank)
             nzval_s(rowcnt_s(ij),ij)= coef_s(2,j,i)


             !sort by col indices
             call insert_sort(jcol_s(:,ij),nzval_s(:,ij),rowcnt_s(ij))
          
             
             !!!!!!!!!!!!North hemisphere (coef direction switched)
             ij = calc_grid_ij(i,jN,lat_rank)

             if (ij > ij_stop_n .or. ij < ij_start_n) then
                write(iulog,*) 'AB: Error ij north index 3 for lhs', mpi_rank, ij
             endif
             
             !coef 4 (i-1, j-1)
             rowcnt_n(ij) = rowcnt_n(ij)+1
             jcol_n(rowcnt_n(ij),ij)= calc_grid_ij(i-1, jN-1, lat_rank)
             nzval_n(rowcnt_n(ij),ij)= coef_n(4,j,i)
             !coef 5 (i-1, j)
             rowcnt_n(ij) = rowcnt_n(ij)+1
             jcol_n(rowcnt_n(ij),ij)= calc_grid_ij(i-1, jN, lat_rank)
             nzval_n(rowcnt_n(ij),ij)= coef_n(5,j,i)
             !coef 6 (i-1, j+1)
             rowcnt_n(ij) = rowcnt_n(ij)+1
             jcol_n(rowcnt_n(ij),ij)= calc_grid_ij(i-1, jN+1, lat_rank)
             nzval_n(rowcnt_n(ij),ij)= coef_n(6,j,i)

             !coef 3 (i, j-1)
             rowcnt_n(ij) = rowcnt_n(ij)+1
             jcol_n(rowcnt_n(ij),ij)= calc_grid_ij(i, jN-1, lat_rank)
             nzval_n(rowcnt_n(ij),ij)= coef_n(3,j,i)
             !coef 9 (i, j)
             rowcnt_n(ij) = rowcnt_n(ij)+1
             jcol_n(rowcnt_n(ij),ij)= calc_grid_ij(i, jN, lat_rank)
             nzval_n(rowcnt_n(ij),ij)= coef_n(9,j,i) 
             !coef 7 (i, j+1)
             rowcnt_n(ij) = rowcnt_n(ij)+1
             jcol_n(rowcnt_n(ij),ij)= calc_grid_ij(i, jN+1, lat_rank)
             nzval_n(rowcnt_n(ij),ij)= coef_n(7,j,i)
             
             !coef 2 (i+1, j-1)
             rowcnt_n(ij) = rowcnt_n(ij)+1
             jcol_n(rowcnt_n(ij),ij)= calc_grid_ij(i+1, jN-1, lat_rank)
             nzval_n(rowcnt_n(ij),ij)= coef_n(2,j,i)
             !coef 1 (i+1, j)
             rowcnt_n(ij) = rowcnt_n(ij)+1
             jcol_n(rowcnt_n(ij),ij)= calc_grid_ij(i+1, jN, lat_rank)
             nzval_n(rowcnt_n(ij),ij)= coef_n(1,j,i)
             !coef 8 (i+1, j+1)
             rowcnt_n(ij) = rowcnt_n(ij)+1
             jcol_n(rowcnt_n(ij),ij)= calc_grid_ij(i+1, jN+1, lat_rank)
             nzval_n(rowcnt_n(ij),ij)= coef_n(8,j,i)

             !sort by col indices
             call insert_sort(jcol_n(:,ij),nzval_n(:,ij),rowcnt_n(ij))
             
          enddo !latitudes (j)

       enddo !longitudes (i)
    
    endif !  between latm_JT and equator

    !equator REGION
    !we could check in the process row, but if we own to the equator then it will be mlat1
    if (mlat1 == nmlat_h) then !I own the equator

       !add this to the south data: jcol_s and nzval_s
       j = nmlat_h       

       !loop through the longitudes in my grid
       do i = mlon0, mlon1
          ij = calc_grid_ij(i,j,lat_rank)
          
          if (ij > ij_stop_s .or. ij < ij_start_s) then
             write(*,*) 'AB: Error ij south index 4 for lhs', mpi_rank, ij
          endif

          
          !coef 6 (i-1, j-1)
          rowcnt_s(ij) = rowcnt_s(ij)+1
          jcol_s(rowcnt_s(ij),ij)= calc_grid_ij(i-1, j-1, lat_rank)
          nzval_s(rowcnt_s(ij),ij)= coef_s(6,j,i)
          !coef 5 (i-1, j)
          rowcnt_s(ij) = rowcnt_s(ij)+1
          jcol_s(rowcnt_s(ij),ij)= calc_grid_ij(i-1, j, lat_rank)
          nzval_s(rowcnt_s(ij),ij)= coef_s(5,j,i)
          !coef 4 (i-1, j+1)
          rowcnt_s(ij) = rowcnt_s(ij)+1
          jcol_s(rowcnt_s(ij),ij)= calc_grid_ij(i-1, j+1, lat_rank)
          nzval_s(rowcnt_s(ij),ij)= coef_s(4,j,i)

          !coef 7 (i, j-1)
          rowcnt_s(ij) = rowcnt_s(ij)+1
          jcol_s(rowcnt_s(ij),ij)= calc_grid_ij(i, j-1, lat_rank)
          nzval_s(rowcnt_s(ij),ij)= coef_s(7,j,i)
          !coef 9 (i, j) ! should be 1
          rowcnt_s(ij) = rowcnt_s(ij)+1
          jcol_s(rowcnt_s(ij),ij)= calc_grid_ij(i, j, lat_rank)
          nzval_s(rowcnt_s(ij),ij)= coef_s(9,j,i)
          !coef 3 (i, j+1)
          rowcnt_s(ij) = rowcnt_s(ij)+1
          jcol_s(rowcnt_s(ij),ij)= calc_grid_ij(i, j+1, lat_rank)
          nzval_s(rowcnt_s(ij),ij)= coef_s(3,j,i)
          
          !coef 8 (i+1, j-1)
          rowcnt_s(ij) = rowcnt_s(ij)+1
          jcol_s(rowcnt_s(ij),ij)= calc_grid_ij(i+1, j-1, lat_rank)
          nzval_s(rowcnt_s(ij),ij)= coef_s(8,j,i)
          !coef 1 (i+1, j)
          rowcnt_s(ij) = rowcnt_s(ij)+1
          jcol_s(rowcnt_s(ij),ij)= calc_grid_ij(i+1, j, lat_rank)
          nzval_s(rowcnt_s(ij),ij)= coef_s(1,j,i)
          !coef 2 (i+1, j+1)
          rowcnt_s(ij) = rowcnt_s(ij)+1
          jcol_s(rowcnt_s(ij),ij)= calc_grid_ij(i+1, j+1, lat_rank)
          nzval_s(rowcnt_s(ij),ij)= coef_s(2,j,i)

          !sort by col indices
          call insert_sort(jcol_s(:,ij),nzval_s(:,ij),rowcnt_s(ij))
          
       enddo ! longitudes
       
    endif !equator

    !------------
    !now each proc has rowcnt_s, jcol_s, nzval_s
    !              rowcnt_n, jcol_n, nzval_n
    ! popular rowptr_s, colind_s and values_s (csr format), then same for north also
    
    ! jcol and nzval contain up to MAX_NNZ entries per grid point, except for grid point 1
    ! which is in jcol1 and nzval1
    !my range of rows (grid points) is ij_start_s: ij_stop_s and ij_start_n: ij_stop_n
    
    !might need to change to 0-based indexing (yes, but happens later for superlu)
    !first SOUTH hemisphere rows
    if (mpi_rank > 0) then !don't own row 1
       rowptr_s(1) = 1
       row_counter_s = 1
       nnz_s = 0

       !loop through grid pts/ matrix rows in southern hemisphere
       do ij = ij_start_s, ij_stop_s
          cnt = rowcnt_s(ij) !entries in row
           if (nnz_s + cnt > size(colind_s)) then
              write(*,*) 'AB: Error: sparse matrix arrays too small nnz_s, cnt = ', nnz_s, cnt
           endif
          do k = 1, cnt
             nnz_s = nnz_s + 1
             colind_s(nnz_s) = jcol_s(k,ij)
             values_s(nnz_s) = nzval_s(k,ij)
          enddo
          row_counter_s = row_counter_s + 1 !add a row
          rowptr_s(row_counter_s) = nnz_s + 1 !pointer to the next row
       enddo
        
    elseif (mpi_rank== 0) then ! I own row 1 (grid point i=1, j=1)
       rowptr_s(1) = 1
       row_counter_s = 1
       nnz_s = 0
       ! first row seperately
       cnt = rowcnt_s(1)
       !write(iulog,*) 'AB: mpi_rank 0, row1 count = ', cnt
       if (cnt + nnz_s > size(colind_s)) then
          write(*,*) 'AB: Error: 1st row, sparse matrix arrays too small, cnt = ', cnt
       endif
       do k = 1, cnt
          nnz_s = nnz_s + 1
          colind_s(nnz_s) = jcol1(k)
          values_s(nnz_s) = nzval1(k)
       enddo
       row_counter_s = row_counter_s + 1
       rowptr_s(row_counter_s) = nnz_s + 1
       !now remaining rows
       !write(iulog,*) "AB: ij_start_s + 1 = (2)", ij_start_s+1
       !write(iulog,*) "AB: ij_stop_s ", ij_stop_s
       !start on row 2
       do ij = ij_start_s+1, ij_stop_s
          cnt = rowcnt_s(ij)
          if (nnz_s + cnt > size(colind_s)) then
             write(*,*) 'AB: Error: sparse matrix arrays too small nnz_s, cnt = ', nnz_s, cnt
          endif
          do k = 1, cnt
             nnz_s = nnz_s + 1
             colind_s(nnz_s) = jcol_s(k,ij)
             values_s(nnz_s) = nzval_s(k,ij)
          enddo
          row_counter_s = row_counter_s + 1
          rowptr_s(row_counter_s) = nnz_s + 1
       enddo
    endif !end of SOUTH
    !now row_counter_s needs to be decremented by 1 so it = #rows of s
    num_row_s = row_counter_s -1
    nnz_south = nnz_s
    !write(*,*) 'AB: mpi_rank, num_rows_s, nnz_south =  ', mpi_rank, num_row_s, nnz_south
    
    !loop through NORTH hemisphere
    rowptr_n(1) = 1
    row_counter_n = 1
    nnz_n = 0

    !loop through grid pts/ matrix rows in north hemisphere
    do ij = ij_start_n, ij_stop_n
       cnt = rowcnt_n(ij)
       if (nnz_n + cnt > size(colind_n)) then
          write(*,*) 'AB: Error: north sparse matrix arrays too small nnz_n, cnt = ', nnz_n, cnt
       endif
       do k = 1, cnt
          nnz_n = nnz_n + 1
          colind_n(nnz_n) = jcol_n(k,ij)
          values_n(nnz_n) = nzval_n(k,ij)
       enddo
       row_counter_n = row_counter_n + 1
       rowptr_n(row_counter_n) = nnz_n + 1
    enddo
    num_row_n = row_counter_n - 1
    nnz_north = nnz_n
    !write(*,*) 'AB: mpi_rank, num_row_n, nnz_north =  ', mpi_rank, num_row_n, nnz_north


    !SET UP BLOCK CSR (exchange with partner)
    if (mpi_size > 1) then
       if (mpi_rank >= 0) then !active proc
          !now set up my block of the csr matrix (my_rowptr, my_values, my_colind)
          !still using 1-based indices
          my_rowptr(1) = 1

          ! now swap hemisphere data with partner
          !for N and S hemi, the even proc rows go first to maintain grid order
          !(so the south owning process)
          if (mod(mpi_rank,2) == 0) then !EVEN, own south, **send north**
             call partner_exchange_hemisphere_mat(MAX_NNZ, rowptr_n, values_n, colind_n, &
                  partner_rowptr, partner_values, partner_cols)
             
             !my south (goes first into csr block)
             if (num_row_s + 1 > size(my_rowptr)) then
                write(*,*) 'AB: Error: my_rowptr array too small for south data, num_row_s +1, size = ', num_row_s+1, size(my_rowptr)
             endif
             do concurrent (i = 2:num_row_s + 1)
                my_rowptr(i) = rowptr_s(i)
             enddo
             !nnz_south set above - but double check
             if (nnz_south /= my_rowptr(num_row_s + 1) - 1) then
                write(*,*) 'AB: Error mpi_rank, nnz_south = ', mpi_rank, nnz_south
             endif
          
             !write(iulog,*) 'AB: again: mpi_rank, num_row_s, nnz_south = ', mpi_rank, num_row_s, nnz_south

             if (nnz_south > size(my_colind) .or. nnz_south > size(my_values)) then
                write(*,*) 'AB: Error: my_colind or my_values array too small for south data'
             endif
          
             do concurrent (i = 1:nnz_south)
                my_colind(i) = colind_s(i)
                my_values(i) = values_s(i)
             enddo
          
             !partner has north
             do i = 1,  partner_hgridsize
                cnt = partner_rowptr(i+1) - partner_rowptr(i)
                !write(iulog,*) 'AB: extra: mpi_rank, i, cnt= ', mpi_rank, i, cnt 
                my_rowptr(num_row_s + 1 + i) = my_rowptr(num_row_s + i) + cnt
             enddo
             
             !north proc's data goes second
             !CHECK mygrid_size = num_row_s + partner_hgridsize
             if (mygrid_size /= num_row_s + partner_hgridsize) then
                write(*, *) 'AB: Error: mygrid_size & partner_hgridsize check', mygrid_size, partner_hgridsize
             endif
          
             !nnz_north needs to be set from partner's info
             nnz_north = my_rowptr(mygrid_size + 1) -1  - nnz_south

             !write(iulog,*) 'AB: again: mpi_rank, nnz_north = ', mpi_rank, nnz_north

             
             if (nnz_south + nnz_north > size(my_colind)) then
                write(*,*) 'AB: Error: my_colind array too small for combined data: mpi_rank, nnz_south, nnz_north, size', mpi_rank, nnz_south, nnz_north, size(my_colind)
             endif
          
             do concurrent (i = 1:nnz_north)
                my_colind(nnz_south + i) = partner_cols(i)
                my_values(nnz_south + i) = partner_values(i)
             enddo
          
          else !ODD, own north, **send south**
             call partner_exchange_hemisphere_mat(MAX_NNZ, rowptr_s, values_s, colind_s, &
                  partner_rowptr, partner_values, partner_cols)
          
             !partner has south - partner's data goes first into csr block
             do concurrent (i = 2:partner_hgridsize + 1)
                my_rowptr(i) = partner_rowptr(i)
             enddo

             !calc based on partner info
             nnz_south = my_rowptr(partner_hgridsize + 1) -1

             do concurrent (i = 1:nnz_south)
                my_colind(i) = partner_cols(i)
                my_values(i) = partner_values(i)
             enddo
          
             !my north - goes second
             do i = 1, num_row_n
                cnt = rowptr_n(i+1) - rowptr_n(i)
                my_rowptr( partner_hgridsize + i + 1) = my_rowptr(partner_hgridsize + i) + cnt
             enddo
          
             !nnz_north already calculated for me above

             do concurrent (i = 1:nnz_north)
                my_colind(nnz_south + i) = colind_n(i)
                my_values(nnz_south + i) = values_n(i)
             enddo
          endif
       endif! active tasks
    else !one task
       !south 
       do concurrent (i = 1:num_row_s + 1)
          my_rowptr(i) = rowptr_s(i)
       enddo
       !nnz_south is calculated above    
       do concurrent (i=1:nnz_south)
          my_colind(i) = colind_s(i)
          my_values(i) = values_s(i)
       enddo
       !north
       do i = 1, num_row_n
          cnt = rowptr_n(i+1) - rowptr_n(i)
          my_rowptr(num_row_s +i + 1) = my_rowptr(num_row_s + i) + cnt
       enddo
       !nnz_north is calc above
       do concurrent (i = 1:nnz_north)
          my_colind(nnz_south + i) = colind_n(i)
          my_values(nnz_south + i) = values_n(i)
       enddo
       
    endif
    
  endsubroutine dist_construct_lhs
!-----------------------------------------------------------------------
  function dist_construct_rhs(coef_10_s, coef_10_n) result(rhs)
! construct vector RHS

    use params_module,only:nmlat_h,nmlat_T1,nmlon
    use cons_module,only:phi_pol
    use mpi_module, only:mlatd0, mlatd1, mlat0, mlat1, mpi_rank, &
         mlon0, mlon1, mlond0, mlond1, &
         lat_rank, lon_rank, partner_hgridsize, my_hgridsize, &
         ij_start_s, ij_stop_s, ij_start_n, ij_stop_n, &
         calc_grid_ij, partner_exchange_hemisphere_vec, mygrid_size, &
         gather_lon_1d, mpi_size, my_sendgrid_size

    real(kind=rp),dimension(mlatd0:mlatd1,mlond0:mlond1),intent(in) :: coef_10_s, coef_10_n
    real(kind=rp),dimension(mygrid_size) :: rhs
    real(kind=rp),dimension(ij_start_s:ij_stop_s) :: rhs_s
    real(kind=rp),dimension(ij_start_n:ij_stop_n) :: rhs_n

    integer :: i,j,ij, jN, j_start, cnt, istart

    real(kind=rp),dimension(nmlon) :: coef10_j1_buf
    
    rhs = 0.0
    rhs_n = 0.0
    rhs_s = 0.0

    !   some debugging info
    !write(*, *) 'AB: RHS mpi_rank, ij_start_s, ij_stop_s, s_grid_pts = ', mpi_rank, ij_start_s, ij_stop_s, ij_stop_s - ij_start_s + 1
    !write(*, *) 'AB: RHS mpi_rank, ij_start_n, ij_stop_n, n_grid_pts = ',mpi_rank, ij_start_n, ij_stop_n,  ij_stop_n - ij_start_n + 1
    !write(*, *) 'AB: RHS mpi_rank, mlat0, mlat1, mlon0,mlon1',mpi_rank, mlat0, mlat1, mlon0,mlon1
    !write(*, *) 'AB: RHS mpi_rank, lon_rank, lat_rank',mpi_rank, lon_rank, lat_rank
    !write(*, *) 'AB: RHS nmlat_h', nmlat_h

    !first do j=1
    if (lat_rank == 0) then ! I own the pole regions (j=1) -this is the first row of procs
       j=1
       !SOUTH POLE
       !proc in lat_rows 0 have to share info with proc 0
       coef10_j1_buf = gather_lon_1d(coef_10_s(j,mlon0:mlon1))
       
       if (lon_rank == 0) then !I also own i=1 (special case for 1 processor)
          ! this is mpi_rank =  0 
          ! for longitude i=1 at the south pole
          i = 1
          ij = calc_grid_ij(i,j, lat_rank) !this should be 1
          
          if (ij > ij_stop_s .or. ij < ij_start_s) then
             write(*,*) 'AB: Error ij south index for rhs', mpi_rank, ij
          endif
          !has to get rest of row from other tasks
          rhs_s(ij) = sum(coef10_j1_buf(:))
       endif

       !nothing else for south pole j=1
       
       !NORTH POLE (still j=1), all i
       do i = mlon0, mlon1
          jN = nmlat_T1-j+1 !j=1, so jN = nmlat_T1
          ij = calc_grid_ij(i, jN, lat_rank)

          if (ij > ij_stop_n .or. ij < ij_start_n) then
             write(*,*) 'AB: Error ij north index for rhs', mpi_rank, ij
          endif
          rhs_n(ij) = phi_pol
       enddo
       !now we've done j=1, so increment
       j_start = 2
    else !done with j=1 for procs owning j=1 (i.e., lat_rank = 0)
       j_start = mlat0
    endif !treatment for j=1

    !populate rhs_s and rhs_n
    !everyone loop through remaining grid points (lat_rank = 0 procs did j=1 already)
    !here we need to not do the equator twice :) (onlysave with the south)
    do i = mlon0,mlon1
       do j = j_start,mlat1
          !SOUTH
          ij = calc_grid_ij(i,j,lat_rank)
          if (ij > ij_stop_s .or. ij < ij_start_s) then
             write(*,*) 'AB: Error ij south index 2 for rhs: rank, ij = ', mpi_rank, ij
          endif
          rhs_s(ij) = coef_10_s(j,i)
          
          !NORTH
          ! don't do the equator in the N
          if (j == nmlat_h) then
             cycle
          endif
          
          jN = nmlat_T1-j+1 !needed to calc ij
          ij = calc_grid_ij(i, jN, lat_rank)
          if (ij > ij_stop_n .or. ij < ij_start_n) then
             write(*,*) 'AB: Error ij north index 2 for rhs: rank, i,j,jN, ij =', mpi_rank, i,j,jN,ij
          endif
          rhs_n(ij) = coef_10_n(j,i)
       enddo
    enddo


    !redistribution with partners
    if (mpi_size > 1) then
       !now do partner hemisphere exchange for continguous rows
       
       !for N and S hemi, the even proc rows go first to maintain grid order
       !(so the south owning process)
       if (mod(mpi_rank,2) == 0) then !even, own south, **send north**
          !copy rhs_s into rhs
          cnt = 0
          do ij = ij_start_s, ij_stop_s
             cnt = cnt + 1
             rhs(cnt) = rhs_s(ij)
             if (isnan(rhs(cnt))) write(*,*) 'AB: ERROR rhs(cnt) is NaN, cnt = ', cnt, ' rank = ', mpi_rank

          enddo
          !send rhs_n and get the rest from my partner
          istart = my_hgridsize + 1
          call partner_exchange_hemisphere_vec(rhs_n, rhs(istart:istart+partner_hgridsize))
          
       else !odd, own north, **send south**
          cnt = partner_hgridsize
          !copy rhs_n into rhs
          do ij = ij_start_n, ij_stop_n
             cnt = cnt + 1
             rhs(cnt) = rhs_n(ij)
             if (isnan(rhs(cnt))) write(*,*) 'AB: ERROR rhs(cnt) is NaN, cnt = ', cnt, ' rank = ', mpi_rank
          enddo
          istart = partner_hgridsize + 1
          !get the rest of the rhs info from partner
          call partner_exchange_hemisphere_vec(rhs_s, rhs(1:partner_hgridsize))
          
       endif
     else !mpi_size = 1
        cnt = 0
        do ij = ij_start_s, ij_stop_s
           cnt = cnt + 1
           rhs(cnt) = rhs_s(ij)
        enddo
        do ij = ij_start_n, ij_stop_n
           cnt = cnt + 1
           rhs(cnt) = rhs_n(ij)
        enddo
     endif
     
  endfunction dist_construct_rhs
 !-----------------------------------------------------------------------

  function dist_solve_superlu(n_global, n_loc, nnz_loc, rowptr, colind, values, rhs) result(sol)

    #include "superlu_dist_config.fh"
    use superlu_mod    
    use mpi_module,only: lat_size,lon_size,dynamo_world,&
         task_csr_rowstarts, mpi_rank, task_csr_mapping, &
         mpi_size
    
    integer,intent(in) :: n_loc,nnz_loc, n_global
    integer(kind=c_int),dimension(n_loc+1),intent(in) :: rowptr
    integer(kind=c_int),dimension(nnz_loc),intent(in) :: colind
    real(kind=c_double),dimension(nnz_loc),intent(in) :: values

    !warning! we assume that kind=rp is same as c_double for values
    !and rhs input

    real(kind=c_double),dimension(n_loc),intent(in) :: rhs
    real(kind=c_double),dimension(n_loc) :: sol
    
    ! for SuperLU sparse matrix solver
    integer,parameter :: nrhs = 1
    integer :: i, iopt, first_row, nprow, npcol

    !superlu structures
    integer(superlu_ptr) :: grid
    integer(superlu_ptr) :: options
    integer(superlu_ptr) :: ScalePermstruct
    integer(superlu_ptr) :: LUstruct
    integer(superlu_ptr) :: SOLVEstruct
    integer(superlu_ptr) :: A
    integer(superlu_ptr) :: stat
      
    ! Other variables
    integer(kind=c_int) :: info, ierr
    real(kind=c_double), target :: berr_array(nrhs)
 
    ! Create Fortran handles for the C structures used in SuperLU_DIST
    call f_create_gridinfo_handle(grid)
    call f_create_options_handle(options)
    call f_dcreate_ScalePerm_handle(ScalePermstruct)
    call f_dcreate_LUstruct_handle(LUstruct)
    call f_dcreate_SOLVEstruct_handle(SOLVEstruct)
    call f_create_SuperMatrix_handle(A)
    call f_create_SuperLUStat_handle(stat)

    write(*, *) 'AB: SUPERLU mpi_rank = ', mpi_rank, 'n_global = ', n_global, 'n_loc = ', n_loc, 'nnz_loc = ', nnz_loc, 'lat_size = ', lat_size, 'lon_size = ', lon_size
           
    ! Initialize the SuperLU_DIST process grid
    !i'll use the same layout as the dynamo
    nprow = lat_size
    npcol = lon_size

    call f_superlu_gridinit(dynamo_world, nprow, npcol, grid)
    if (nprow * npcol /= mpi_size) then
       write(*,*) "ERROR: nprow*npcol != nprocs"
    endif
    
    !get my first row in distributed matrix
    first_row = task_csr_rowstarts(task_csr_mapping(mpi_rank))     !these are 0-based already 
    
    !create the distributed compressed row matrix A (O-index)
    !some debugging
    if (rowptr(n_loc+1) /= nnz_loc) then
       write(*,*)  "AB: ERROR rowptr(n_loc+1) /= nnz_loc ", rowptr(n_loc+1), nnz_loc
    endif
    if (minval(colind) < 0 .or. maxval(colind) >= n_global) then
       write(*,*)  "AB: Error colind out of bounds: min(col) max(col), n_global" , minval(colind), maxval(colind), n_global
    endif
    write(*,*) 'AB: SUPERLU first_row (0-index) = ', first_row
    if (first_row < 0 .or. first_row >= n_global) then
       write(*,*)  "first_row out of range", first_row, n_global
    endif

    call f_dCreate_CompRowLoc_Mat_dist(A, n_global, n_global, nnz_loc, n_loc, first_row, &
         values, colind, rowptr, SLU_NR_loc, SLU_D, SLU_GE) 

    ! Setup the right hand side (rhs contains local data)
    sol=rhs ! Copy RHS to solution vector

    ! Set the default input options
    call f_set_default_options(options)

    ! Change one or more options
    !these below are the defaults
    !call set_superlu_options(options,ColPerm=MMD_AT_PLUS_A)
    !call set_superlu_options(options,RowPerm=LargeDiag_MC64)
    !refinement: (or none = 0 or  single = 1, double = 2)
    !call set_superlu_options(options, IterRefine = 2)
    ! Optionally, enable equilibration/scaling for stability
    ! 1 - on, 0 = off
    !call set_superlu_options(options, Equil=1)

    ! Initialize ScalePermstruct and LUstruct
    call get_SuperMatrix(A, nrow=n_global, ncol=n_global)
    call f_dScalePermstructInit(n_global, n_global, ScalePermstruct)
    call f_dLUstructInit(n_global, n_global, LUstruct)
 
    ! Initialize the statistics variables
    call f_PStatInit(stat)
 
    ! Call the linear equation solver (writes over rhs (sol))
    call f_pdgssvx(options, A, ScalePermstruct, sol, n_loc, nrhs, &
         grid, LUstruct, SOLVEstruct, berr_array, stat, info)
    
    if (info /= 0) then
       write(iulog,*) 'ERROR: pdgssvx failed with mpi_rank, INFO = ', mpi_rank, info
    endif
    if (info == 0 .and. mpi_rank == 0) then
       write(iulog,*) 'Backward error: ', berr_array(1)
    endif

    write(*,*) 'DONE WITH PDGSSSVX'
    
    ! result is sol (already assigned by reference in pdgssvx)

    !TO DO Can I save some of the LU structures after the first time since
    !the sparsity pattern of A will not change

    !  deallocate the storage allocated by SuperLU_DIST
    call f_PStatFree(stat)
    !do not call f_Destroy_CompRowLoc_Mat_dist(A)
    ! - tries to free the fortran-allocated rowptr,colind and nzval array

    !TEMP: temp because of tree issues
    !call f_dDestroy_LU_SOLVE_struct(options, n, grid, LUstruct, SOLVEstruct)

    call f_dScalePermstructFree(ScalePermstruct)
    ! Release the SuperLU process grid
    call f_superlu_gridexit(grid)

    ! Deallocate the C structures pointed to by the Fortran handles
    call f_destroy_SuperLUStat_handle(stat)    
    call f_destroy_SOLVEstruct_handle(SOLVEstruct)
    call f_destroy_LUstruct_handle(LUstruct)
    call f_destroy_ScalePerm_handle(ScalePermstruct)
    call f_destroy_options_handle(options)
    call f_destroy_SuperMatrix_handle(A)
    call f_destroy_gridinfo_handle(grid)

    
  endfunction dist_solve_superlu

!-----------------------------------------------------------------------
  function dist_flatten(fin) result(fout)
! reorder 2D fields (lat-lon) into 1D vector (RHS)
! northern/southern hemispheres are either separate or averaged
! based on their latitude ranges (high-lat, transition, low-lat, equator)

    use params_module,only:nmlat_h,nmlat_T1,nmlon
    use cons_module,only:jlatm_JT
    use mpi_module, only:lat_rank, mlat0, mlat1, &
         mlon0, mlon1, mlatd0, mlatd1, mlond0, mlond1, &
         ij_start_n, ij_stop_n, mpi_size, mpi_rank, &
         ij_start_s, ij_stop_s, partner_hgridsize, mygrid_size, &
         calc_grid_ij, my_hgridsize, partner_exchange_hemisphere_vec, &
         my_sendgrid_size

    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: fin
    real(kind=rp),dimension(mygrid_size) :: fout
    
    real(kind=rp),dimension(ij_start_s:ij_stop_s) :: fout_s
    real(kind=rp),dimension(ij_start_n:ij_stop_n) :: fout_n
    integer :: i,j,ij, loop_start_j, loop_stop_j, jS, jN, cnt
    integer :: istart
    real(kind=rp) :: avg

    real(kind=rp),dimension(my_sendgrid_size) :: sendbuf

    ! Initialize arrays to avoid uninitialized values
    fout_s = 0.0
    fout_n = 0.0
    sendbuf = 0.0
    fout = 0.0
    
    if (mlat0<jlatm_JT) then !includes jlatm_JT
       ! from pole to jlatm_JT, two hemispheres are uncoupled
       ! my longitudes
       loop_start_j = mlat0
       loop_stop_j = min(mlat1, jlatm_JT)

       do concurrent (i = mlon0:mlon1, j=loop_start_j:loop_stop_j)
          !south
          jS=j
          ij = calc_grid_ij(i,jS,lat_rank)
          !checkbounds
          if (ij >= ij_start_s .and. ij <= ij_stop_s) then
             fout_s(ij) = fin(1,j,i)
             if (isnan(fout_s(ij))) write(*,*) 'AB: fout_s(ij) is NaN, ij, i, jS = ', ij, i, jS
           else
             write(iulog,*) "Error in flatten south, mpirank, ij = ", mpi_rank, ij
          endif
          
          !north
          jN = nmlat_T1-j+1
          ij = calc_grid_ij(i,jN,lat_rank)
          !checkbounds
          if (ij >= ij_start_n .and. ij <= ij_stop_n) then
             fout_n(ij) = fin(2,j,i)
             if (isnan(fout_n(ij))) write(*,*) 'AB: fout_n(ij) is NaN, ij, i, jN = ', ij, i, jN
          else
             write(iulog,*) "Error in flatten north, mpirank, ij = ", mpi_rank, ij
          endif
       enddo
    endif
    
    if ((mlat0 > jlatm_JT) .or. (mlat1 > jlatm_JT)) then
       ! from jlatm_JT+1 to equator, symmetric solution
       loop_start_j = max(mlat0, jlatm_JT+1)
       loop_stop_j = mlat1

       do concurrent (i = mlon0:mlon1, j = loop_start_j:loop_stop_j)
          avg = (fin(1,j,i)+fin(2,j,i))/2

          !south
          jS=j
          ij = calc_grid_ij(i,jS,lat_rank)
          !checkbounds
          if (ij >= ij_start_s .and. ij <= ij_stop_s) then
             fout_s(ij) = avg
             if (isnan(fout_s(ij))) write(*,*) 'AB: fout_s(ij) is NaN, ij, i, jS = ', ij, i, jS
          else
             write(iulog,*) "Error in flatten 2 south, mpirank, ij = ", mpi_rank, ij
          endif
             
          !north
          ! don't do the equator in the N 
          if (j == nmlat_h) then
             cycle
          endif
          jN = nmlat_T1-j+1
          ij = calc_grid_ij(i,jN,lat_rank)
          !checkbounds
          if (ij >= ij_start_n .and. ij <= ij_stop_n) then
             fout_n(ij) = avg
             if (isnan(fout_n(ij))) write(*,*) 'AB: fout_n(ij) is NaN, ij, i, jN = ', ij, i, jN

          else
             write(iulog,*) "Error in flatten 2 north, mpirank, ij = ", mpi_rank, ij
          endif

       enddo
    endif

    if (mpi_size == 1) then
       !now combine the rows so they are the same order as in matrix
       cnt = 0
       do ij = ij_start_s, ij_stop_s
          cnt = cnt + 1
          fout(cnt) = fout_s(ij)
       enddo
       do ij = ij_start_n, ij_stop_n
          cnt = cnt + 1
          fout(cnt) = fout_n(ij)
       enddo
    else !multiple procs

       write(*,*) 'AB: flatten: mygrid_size,  my_sendgrid_size, my_hgridsize ', mygrid_size,  my_sendgrid_size, my_hgridsize 
       
       !now do partner hemisphere exchange for continguous rows
       !for N and S hemi, the even proc rows go first to maintain grid order
       !(so the south owning process)
       if (mod(mpi_rank,2) == 0) then !even, own south, **send north**
          !copy south into fout
          cnt = 0
          do ij = ij_start_s, ij_stop_s
             cnt = cnt + 1
             fout(cnt) = fout_s(ij)
             if (isnan(fout(cnt))) write(*,*) 'AB: #1 fout(cnt) is NaN, cnt, ij = ', cnt, ij

          enddo
          !copy north to send
          cnt = 0
          do ij = ij_start_n, ij_stop_n
             cnt = cnt + 1
             sendbuf(cnt) = fout_n(ij)
          enddo
          
          istart = my_hgridsize + 1
          !DEBUG
          do i = 1,my_sendgrid_size
             if (isnan(sendbuf(i))) write(*,*) 'AB: SEND 1 sendbuf(i) is NaN, i = ', i
          enddo
          !get partner's north
          call partner_exchange_hemisphere_vec(sendbuf, fout(istart:istart+partner_hgridsize))
          
       else !odd, own north, **send south**
          !copy north to fout
          cnt = partner_hgridsize
          do ij = ij_start_n, ij_stop_n
             cnt = cnt + 1
             fout(cnt) = fout_n(ij)
             if (isnan(fout(cnt))) write(*,*) 'AB: #2 fout(cnt) is NaN, cnt, ij = ', cnt, ij
          enddo

          !copy south to sendbuf
          cnt = 0
          do ij = ij_start_s, ij_stop_s
             cnt = cnt + 1
             sendbuf(cnt) = fout_s(ij)
          enddo
          istart = partner_hgridsize + 1
          do i = 1, my_sendgrid_size
             if (isnan(sendbuf(i))) write(*,*) 'AB: SEND 2 sendbuf(i) is NaN, i = ', i
          enddo
          call partner_exchange_hemisphere_vec(sendbuf, fout(1:partner_hgridsize))
       endif
    endif ! multiple procs
       
  endfunction dist_flatten

!-----------------------------------------------------------------------
  function dist_unravel(fin) result(fout)
    ! reorder 1D vector (RHS) into 2D fields (lat-lon)
    ! 1d vector has the reordering for process rows, so we have to undo that also
    ! this is the reverse of the dist_flatten
    
    use params_module,only:nmlat_h,nmlat_T1,nmlon
    use mpi_module, only:lat_rank, mlat0, mlat1, &
         mlon0, mlon1, &
         ij_start_n, ij_stop_n, &
         ij_start_s, ij_stop_s, lat_rank, partner_hgridsize, &
         mygrid_size, mpi_rank, calc_grid_ij, &
         partner_exchange_hemisphere_vec, my_sendgrid_size


    real(kind=rp),dimension(mygrid_size),intent(in) :: fin
    real(kind=rp),dimension(2,mlat0:mlat1,mlon0:mlon1) :: fout

    integer :: i,j,ij, cnt, jS, jN
    real(kind=rp),dimension(ij_start_s:ij_stop_s) :: fin_s
    real(kind=rp),dimension(ij_start_n:ij_stop_n) :: fin_n
    real(kind=rp),dimension(partner_hgridsize) :: buffer

    !output 3D array for n=2 & s=1
    fout = 0.0

    !first reverse the hemisphere swap
    if (mod(mpi_rank,2) == 0) then !EVEN, own then south, send other
                                   !half back to partner
       !south
       cnt = 0
       !this is mine for the south, so copy
       do ij = ij_start_s, ij_stop_s
          cnt = cnt + 1
          fin_s(ij) = fin(cnt)
       enddo
       !now put south into 3D array
       do concurrent (i = mlon0:mlon1, j = mlat0:mlat1)
          jS = j
          ij =  calc_grid_ij(i,jS,lat_rank)
          if (ij >= ij_start_s .and. ij <= ij_stop_s) then
             fout(1,j,i) = fin_s(ij)
             if (j == nmlat_h) then !equator, set in north also
                fout(2,j,i) = fin_s(ij)
             endif
          else
             write(iulog,*) "Error in unravel south 1, mpirank, ij, ij_start, ij_stop = ", mpi_rank, ij, ij_start_s, ij_stop_s
          endif
       enddo
       
       !send other half info to partner and receive my north data
       do i=1, partner_hgridsize
          if ((cnt +i) <= mygrid_size) then 
             buffer(i) = fin(cnt + i)
          else
             write(iulog,*) "Error in unravel buffer, mpirank, cnt+i = ", mpi_rank, cnt+i
          endif
       enddo

       call partner_exchange_hemisphere_vec(buffer, fin_n)

       !now fin_n contains my north data, so put back into 3D array   
       do concurrent (i = mlon0:mlon1, j = mlat0:mlat1)
          jN = nmlat_T1-j+1 !north j
          if (j == nmlat_h) then !no equator data from north in fin_n
             !(already populated)
             cycle
          endif
          ij =  calc_grid_ij(i,jN,lat_rank)
          if (ij >= ij_start_n .and. ij <= ij_stop_n) then
             fout(2,j,i) = fin_n(ij)
           else
             write(iulog,*) "Error in unravel north 1, mpirank, ij = ", mpi_rank, ij, ij_start_n, ij_stop_n
          endif
       enddo
       
    else !ODD, so own north, get south back from  partner
       !first bit to send buffer
       do i=1, partner_hgridsize
          buffer(i) = fin(i)
       enddo
       
       call partner_exchange_hemisphere_vec(buffer, fin_s)

       do concurrent (i = mlon0:mlon1, j = mlat0:mlat1)
          jS = j
          ij =  calc_grid_ij(i,jS,lat_rank)
          if (ij >= ij_start_s .and. ij <= ij_stop_s) then
             fout(1,j,i) = fin_s(ij)
             if (j == nmlat_h) then !equator, set in north also
                fout(2,j,i) = fin_s(ij)
             endif
          else
             write(iulog,*) "Error in unravel south 2, mpirank, ij = ", mpi_rank, ij, ij_start_s, ij_stop_s
          endif
       enddo
       
       !north (I own, so copy)
       cnt = partner_hgridsize
       do ij = ij_start_n, ij_stop_n
          cnt = cnt + 1
          fin_n(ij) = fin(cnt)
       enddo
       !put in 3D array
       do concurrent (i = mlon0:mlon1, j = mlat0:mlat1)
          jN = nmlat_T1-j+1 !north index
          if (j == nmlat_h) then !no equator data from north in fin_n
             !(already populated)
             cycle
          endif
          ij =  calc_grid_ij(i,jN,lat_rank)
          if (ij >= ij_start_n .and. ij <= ij_stop_n) then
             fout(2,j,i) = fin_n(ij)
          else
             write(iulog,*) "Error in unravel north 2, mpirank, ij, ij_start, ij_stop = ", mpi_rank, ij, ij_start_n, ij_stop_n
          endif
          
       enddo
    endif !end north
       
  endfunction dist_unravel
  
   
 
  !-----------------------------------------------------------------------

  subroutine insert_sort(array_i, array_r, len)
   ! this is only an ok sorting approach for small arrays
   ! since its O(n^2)
   !LIMITED to MAX_NNZ in length
   ! sorting by the int array but moving the reals 
   ! then removes zeros
     integer, dimension(MAX_NNZ), intent(inout) :: array_i
     real(kind=rp), dimension(MAX_NNZ), intent(inout) :: array_r
     integer, intent(inout) :: len

     integer :: i, j, temp_i, z
     real(kind=rp) :: temp_r
     integer :: keep(MAX_NNZ)

     do i=2,len
       temp_i = array_i(i)
       temp_r = array_r(i)
       do j=i-1,1,-1
          if (array_i(j).le.temp_i) exit
          array_i(j+1) = array_i(j)
          array_r(j+1) = array_r(j)
       enddo
       array_i(j+1) = temp_i
       array_r(j+1) = temp_r
    enddo

    !now remove zeros
    keep = 1
    do i=1,len
       if (array_r(i) ==  0.0) then
          keep(i) = 0
       endif
    enddo
    z = 0
    do i=1,len
       if (keep(i) == 1) then
          z = z+1
          array_i(z) = array_i(i)
          array_r(z) = array_r(i)
       endif
    enddo
    len = z

  endsubroutine insert_sort

!-----------------------------------------------------------------------------
  
endmodule dist_solver_module
