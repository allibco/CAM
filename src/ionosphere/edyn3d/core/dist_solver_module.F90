module dist_solver_module
  use perf_mod, only: t_startf, t_stopf
  use superlu_mod    

  use prec,only:rp
  use iso_c_binding
  use cam_logfile, only: iulog
  use dist_spmv_mod

  include 'netcdf.inc'
  #include "superlu_dist_config.fh"
 

  !max nonzeros per row (this does not incl the dense row at the pole)
  integer, parameter :: MAX_NNZ=12
  !if you want output the first interation, set to 0, otherwise set to 1
  integer :: output_matrix_count=1

  !only want to setup the halo once for the spmv
  type(halo_t) :: halo

  !superlu (here to faciltate reuse across iterations)
  logical, save :: superlu_initialized = .false.

  integer(superlu_ptr), save :: grid
  integer(superlu_ptr), save :: options
  integer(superlu_ptr), save :: ScalePermstruct
  integer(superlu_ptr), save :: LUstruct
  integer(superlu_ptr), save :: SOLVEstruct
  integer(superlu_ptr), save :: A
  integer(superlu_ptr), save :: stat

  ! to reuse matrix structure for superlu, we need to only allocate block csr_matix
  ! arrays once
  integer,dimension(:), allocatable, save :: g_rowptr
  integer,dimension(:), allocatable, save :: g_colind
  real(kind=rp),dimension(:), allocatable, save :: g_values_csr
  integer, save :: g_n_global, g_n_loc, g_nnz_loc, g_first_row

  integer(kind=c_int64_t), save :: g_pattern_hash = 0_c_int64_t
  logical, save :: pattern_initialized = .false.
  
  contains
!-----------------------------------------------------------------------
  subroutine dist_linear_system(mlatd0,mlatd1,mlond0,mlond1, &
    bij,pot_hl,fac_hl,coef_ns,pot)
    ! construct linear system based on bij and coef and solve in pot

    ! note: this is only called by the active processes (in union_world)
    
    ! if FAC is read in, pot_hl is not used, only fac_hl is used
    ! if potential is read in, pot_hl is used, fac_hl is output

    
    use params_module, only:nmlat_h,nmlat_T1,nmlon
    use cons_module, only:read_fac
    use mpi_module, only: mpi_rank,dynamo_world,union_world, lat_rank,lon_rank,&
                          nmlat_task,nmlon_task,mygrid_size,&
                          sync_mlat_5d, sync_mlon_5d,& 
                          sync_mlat_3d, sync_mlon_3d, &
                          mlat0, mlat1, mlon0, mlon1, &
                          task_csr_rowstarts, mpi_size, &
                          un_mpi_size, un_mpi_rank

    
    ! the processor grid only covers one hemisphere (nmlat_h, nmlon)
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

    !integer,dimension(:), allocatable :: rowptr
    !integer,dimension(:), allocatable :: colind
    !real(kind=rp),dimension(:), allocatable :: values_csr
    real(kind=rp),dimension(10,mlatd0:mlatd1,mlond0:mlond1) :: coef_s
    real(kind=rp),dimension(10,mlatd0:mlatd1,mlond0:mlond1) :: coef_n

    real(kind=rp),dimension(:), allocatable :: rhs,z,pot_hl_f,sol

    integer :: ierr, fst_row

    !type(halo_t) :: halo
    
    !for optional output
    logical :: output_matrix = .false.
    integer :: fileid, Xid, potid
    character(len=200) :: fbuf
    character(kind=c_char,len=:), allocatable :: newfile

    if (output_matrix_count == 0) then
       output_matrix = .true.
       output_matrix_count = output_matrix_count + 1
       !write(*,*) 'Dist_ls: Output matrix this time through ...'
    else
       output_matrix = .false.
    endif
    
    call t_startf('dist_linear_system')

    ! global matrix size
    nlonlat = nmlat_T1*nmlon 

    !proc domain extents without halo pts (these are global: mlat0, mlon0, mlat1, mlon1)

    !number of grid points I will own (after hemisphere exchange) is mygrid_size (global var)

    !allocate space for g_rowptr,g_colind,g_values_csr
    ! (only first time through)
    ! why does MAX_NNZ=12? seems like 10 is max?
    
    if (mpi_rank == 0) then ! make room for dense row
       nnz_est = (mygrid_size-1)*MAX_NNZ + (nmlon + 2)
    else
       nnz_est = mygrid_size*MAX_NNZ
    endif

    !allocate only the first time through
    if (.not. superlu_initialized) then
       allocate(g_rowptr(mygrid_size+1))
       allocate(g_colind(nnz_est))
       allocate(g_values_csr(nnz_est))
    endif
       
    allocate(rhs(mygrid_size))
    allocate(z(mygrid_size))
    allocate(pot_hl_f(mygrid_size))
    allocate(sol(mygrid_size))

    rhs = 0.0
    sol = 0.0
    pot_hl_f = 0.0
    g_colind = 0
    g_row_ptr = 0
    g_values_csr = 0.0
    z = 0.0
    
    ! for now, split two hemispheres (keep halo pts)
    do concurrent (i = mlond0:mlond1, j = mlatd0:mlatd1, ic = 1:10)
       coef_s(ic,j,i) = coef_ns(ic,1,j,i) !south
       coef_n(ic,j,i) = coef_ns(ic,2,j,i) !north (note: row = nmlat_T1-j+1) 
    enddo

    ! construct LHS matrix in Block CSR format
    ! Note: rows are contiguous on each proc (for LHS matrix
    ! and RHS) 
    call dist_construct_lhs(nnz_est,bij,coef_s(1:9,:,:),coef_n(1:9,:,:),g_rowptr,g_colind,g_values_csr)
    !need nnz for solver
    nnz = g_rowptr(mygrid_size+1)-1
    
    ! RHS in Block format to match LHS
    rhs = dist_construct_rhs(coef_s(10,:,:), coef_n(10,:,:))

    !change matrix to 0-based indexing 
    !need matrix to be 0-based index for superlu and the matvec
    g_colind=g_colind-1
    g_rowptr=g_rowptr-1

    ! determine FAC forcing (dense)
    if (read_fac) then ! input is corrected fac_hl, pot_hl is not used

       z = dist_flatten(fac_hl)

    else ! input is pot_hl, fac_hl is to be calculated (output)

       ! A. Maute 2023/11/21: put the high latitude potential in X
       ! and then use LHS to calculate the RHS FAC
       pot_hl_f = dist_flatten(pot_hl)
       
       !Parallel matmult (uses union group)
       ! z = matmul(lhs, pot_hl_f)
       fst_row = task_csr_rowstarts(un_mpi_rank) ! 0-index

       !This init should be called just the first timestep because the nonzero
       !matrix pattern does not change
       if (halo%nprocs_local == 0) then
          call dist_spmv_init(mygrid_size, fst_row, nlonlat, un_mpi_size, un_mpi_rank, g_rowptr, g_colind, task_csr_rowstarts, union_world, halo, ierr)
       endif
          
       !optional DEBUG
       !call write_halo_to_file(halo, union_world)

       call dist_spmv(g_rowptr, g_colind, g_values_csr, pot_hl_f, z, halo, union_world, ierr)

       !free later (only after done because we are just initializing once)
       !call dist_spmv_free(halo)

       !write(*,*) 'Dist_ls: Finished spmv'
          
       ! reconstruct 2D distribution of FAC based on z
       fac_hl(:,mlat0:mlat1,mlon0:mlon1) = dist_unravel(z)

       !write(*,*) 'Dist_ls: Finished unravel'
       
       !get ghost/halo points (only dynamo procs)
       if (mpi_rank >=0 ) then
          if (mpi_size > 0) then  
             call sync_mlat_3d(fac_hl(:,:,mlon0:mlon1), 2)
             call sync_mlon_3d(fac_hl, 2)
          else ! add periodic points                                                     
             do j = 1,nmlat_h
                do isn = 1,2
                   fac_hl(isn,j,0) = fac_hl(isn,j,nmlon)
                   fac_hl(isn,j,nmlon+1) = fac_hl(isn,j,1)
                enddo
             enddo
          endif
       endif !dynamo procs
     endif !FAC

     !write(*,*) 'Dist_ls: Finished fac section'
     
     ! add FAC forcing to RHS
     !(these are both set for contiguous rows already - all union procs own)
     do i = 1,mygrid_size
        !if (isnan(rhs(i))) write(*,*) 'Dist_ls warning: rhs(i) is NaN, i = ', i
        !if (isnan(z(i))) write(*,*) 'Dist_ls warning: z(i) is NaN, i = ', i 
        rhs(i) = rhs(i)+z(i)
     enddo

     !do we want to output matrix and rhs for debugging
     if (output_matrix == .true. ) then
        write(*,*) 'Dist_ls: writing netcdf file ...'
        write(fbuf,'(A,".",I4.4,".nc")') &
             '/glade/derecho/scratch/abaker/dynamo_test', un_mpi_rank
        !trim it
        newfile = trim(fbuf)//c_null_char

        !call write fcun (1 or 2) - phase 0
        call write_dist_to_file(newfile, 0, mlatd0,mlatd1,mlond0,mlond1, nnz, g_colind, g_rowptr, g_values_csr, pot_hl_f, z, rhs,fac_hl, pot, sol, fileid, Xid, potid)

      endif
     
     call t_startf('linear_system->solve_superlu')
     sol = dist_solve_superlu(nlonlat, mygrid_size, nnz, g_rowptr, g_colind(1:nnz), g_values_csr(1:nnz), rhs)
     call t_stopf('linear_system->solve_superlu')

     !print *, 'Dist_ls: done with solve'

     
     ! reconstruct 2D distribution of potential based on the solution
     pot(:,mlat0:mlat1,mlon0:mlon1) = dist_unravel(sol)

     !get ghost/halo points (only dynamo procs)
     if (mpi_rank >=0 ) then
        if (mpi_size > 0) then
           call sync_mlat_3d(pot(:,:,mlon0:mlon1), 2)
           call sync_mlon_3d(pot, 2)
        else
           !periodic points
           do j = 1,nmlat_h
              do isn = 1,2
                 pot(isn,j,0) = pot(isn,j,nmlon)
                 pot(isn,j,nmlon+1) = pot(isn,j,1)
              enddo
           enddo
        endif
     endif
     
     !print *, 'Dist_ls: done with unravel'

     !if output turned on for debugging 
     if (output_matrix == .true.) then
        call write_dist_to_file(newfile, 1, mlatd0,mlatd1,mlond0,mlond1, nnz, g_colind, g_rowptr, g_values_csr, pot_hl_f, z, rhs, fac_hl, pot,sol, fileid, Xid, potid)
        !print *, 'Dist_ls: done with netcdf file'
     endif

     !pot and fac_hl are ready to return (updated grid + halo)

     ! Deallocate arrays to free memory
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
         un_mpi_size, un_mpi_rank, union_world, ex_mpi_rank, &
         nmlat_task,nmlon_task,mlatd0,mlatd1,mlat0,mlat1, &
         mlond0,mlond1,mlon0,mlon1, &
         lat_size,lon_size,task_lat_offset,ij_start_n,ij_stop_n, &
         ij_start_s,ij_stop_s,my_recvgrid_size, &
         calc_grid_ij, mpi_sendnorth_mat, &
         gather_lon_1d, mygrid_size, mpi_partner, &
         mygrid_size_n, mygrid_size_s, my_sendgrid_size
    
    integer,intent(in) :: nnz_est
    real(kind=rp),dimension(mlatd0:mlatd1,mlond0:mlond1),intent(in) :: bij
    real(kind=rp),dimension(9,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: coef_s
    real(kind=rp),dimension(9,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: coef_n
    integer,dimension(mygrid_size+1),intent(out) :: my_rowptr
    integer,dimension(nnz_est),intent(out) :: my_colind
    real(kind=rp),dimension(nnz_est),intent(out) :: my_values

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
    integer, dimension(my_recvgrid_size+1) :: partner_rowptr
    real(kind=rp),dimension(my_recvgrid_size*MAX_NNZ) :: partner_values
    integer, dimension(my_recvgrid_size*MAX_NNZ) :: partner_cols

    !global size 
    nlonlat = nmlat_T1*nmlon

    !---------------------------------
    if (mpi_rank >= 0) then !just dynamo procs for getting the coefficients

       !-------------------------------------

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
       ! each proc has to have a continguous block of rows, so we'll have to send the
       ! north  hemispher data to the extra procs later
    
       ! Set up matrix with loops through the following regions:
       ! poles, lower/upper lats (through latm_JT), then latm_Jt, through equator

       !coefficients:
       ! ^   equatorward
       ! coef(4) (i-1,j+1)      coef(3) (i,j+1)    coef(2) (i+1,j+1)
       ! coef(5) (i-1,j)        coef(9) (i,j)      coef(1) (i+1,j)
       ! coef(6) (i-1,j-1)      coef(7) (i,j-1)    coef(8) (i+1,j-1)
       ! v   poleward


       ! some initialization
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
       
    
       !The outer if statements to see if a proc owns rows in that region

       !POLES------------------------------------------_!
       if (lat_rank == 0) then ! I own the pole regions (j=1)
          j = 1
          ! there are no c6,c7,c8 values at the south pole

          !do communication for coef(3,1,i) here so that root proc can complete
          !its row below
          !use gather_lon_1d()

          !the proc that owns j=1,i=1 (proc 0)  needs  coef(3,j=1,isub)
          !for isub=2:nmlon, but only owns 2:mlon1
          !so we need to communicate within proc lat_rank = 0 and gather to
          ! proc 0 
          coef3_j1_buf = gather_lon_1d(coef_s(3,j,mlon0:mlon1))
          !also needs coef9  for the whole row corres to j=1
          coef9_j1_buf = gather_lon_1d(coef_s(9,j,mlon0:mlon1))

          counter = 0
          if (lon_rank == 0) then !I also own i=1 (special case - applies to 1 processor)
             i=1
             counter = counter + 1
             jcol1(counter) = calc_grid_ij(i,j,0) !this will be 1
             !this requires the whole row (needs to get the rest)
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
                counter = counter +1
                jcol1(counter) =  calc_grid_ij(isub,j+1,0)
                nzval1(counter) = coef3_j1_buf(isub)
             enddo

             rowcnt_s(1) = counter !should be = nmlon+2
             !jcol1 will be sorted already
             !but *not* for > 1 process
             !sort by col indices
             call insert_sort(jcol1, nzval1, counter)
             
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
             ! this doesn't make intuive sense to me .. but matched original
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
                   write(*,*) 'Dist_ls: Error ij south index 1 for lhs', mpi_rank, ij
                endif
                !note: the edges of the proc domain are handled in calc_grid_ij
             
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
                   write(iulog,*) 'Dist_ls: Error ij north index 1 for lhs', mpi_rank, ij
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
                write(*,*) 'Dist_ls: Error ij south index 2 for lhs', mpi_rank, ij
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
             jcol_s(rowcnt_s(ij),ij)= calc_grid_ij(i+1, jN+1, lat_rank)
             nzval_s(rowcnt_s(ij),ij)= coef_n(8,j,i)

             !sort by col indices
             call insert_sort(jcol_s(:,ij),nzval_s(:,ij),rowcnt_s(ij))
          
             !!Northern Hemi

             !calc_ij for North hemisphere
             ij = calc_grid_ij(i,jN,lat_rank)
             if (ij > ij_stop_n .or. ij < ij_start_n) then
                write(*,*) 'Dist_ls: Error ij north index 2 for lhs', mpi_rank, ij
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
                   write(*,*) 'Dist_construct_lhs: Error ij south index 3 for lhs', mpi_rank, ij
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
                   write(iulog,*) 'Dist_construct_lhs: Error ij north index 3 for lhs', mpi_rank, ij
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
                write(*,*) 'DISH_LHS: Error ij south index 4 for lhs', mpi_rank, ij
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
       ! now each proc has rowcnt_s, jcol_s, nzval_s
       !              rowcnt_n, jcol_n, nzval_n
       ! populate rowptr_s, colind_s and values_s (csr format), then same for north also
       ! jcol and nzval contain up to MAX_NNZ entries per grid point, except for grid point 1
       ! which is in jcol1 and nzval1
       ! my range of rows (grid points) is ij_start_s: ij_stop_s and ij_start_n: ij_stop_n
    
       !first SOUTH hemisphere rows
       if (mpi_rank > 0) then !don't own row 1
          rowptr_s(1) = 1
          row_counter_s = 1
          nnz_s = 0

          !loop through grid pts/ matrix rows in southern hemisphere
          do ij = ij_start_s, ij_stop_s
             cnt = rowcnt_s(ij) !entries in row
             if (nnz_s + cnt > size(colind_s)) then
                write(*,*) 'Dist_construct_lhs: Error: sparse matrix arrays too small nnz_s, cnt = ', nnz_s, cnt
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
          ! do first row seperately
          cnt = rowcnt_s(1)
          if (cnt + nnz_s > size(colind_s)) then
             write(*,*) 'Dist_construct_lhs: Error: 1st row, sparse matrix arrays too small, cnt = ', cnt
          endif
          do k = 1, cnt
             nnz_s = nnz_s + 1
             colind_s(nnz_s) = jcol1(k)
             values_s(nnz_s) = nzval1(k)
          enddo
          row_counter_s = row_counter_s + 1
          rowptr_s(row_counter_s) = nnz_s + 1
          !now remaining rows
          !start on row 2
          do ij = ij_start_s+1, ij_stop_s
             cnt = rowcnt_s(ij)
             if (nnz_s + cnt > size(colind_s)) then
                write(*,*) 'Dist_construct_lhs: Error: sparse matrix arrays too small nnz_s, cnt = ', nnz_s, cnt
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
    
       !loop through NORTH hemisphere
       rowptr_n(1) = 1
       row_counter_n = 1
       nnz_n = 0

       !loop through grid pts/ matrix rows in north hemisphere
       do ij = ij_start_n, ij_stop_n
          cnt = rowcnt_n(ij)
          if (nnz_n + cnt > size(colind_n)) then
             write(*,*) 'Dist_construct_lhs: Error: north sparse matrix arrays too small nnz_n, cnt = ', nnz_n, cnt
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

    endif !this is the end of just the dynamo procs

    ! SET UP BLOCK CSR (and send north hemisphere to partner)
    ! We need all the procs
    if (un_mpi_size > 1) then
       
       call mpi_sendnorth_mat(MAX_NNZ, rowptr_n, values_n, colind_n, &
            partner_rowptr, partner_values, partner_cols)

       !note: still using 1-based indices
       my_rowptr(1) = 1
       
       if (mpi_rank >= 0) then !dynamo proc
          !now set up my block of the csr matrix (my_rowptr, my_values, my_colind)
          !from the south data
          if (num_row_s /= mygrid_size) then
             write(*,*) 'Dist_construct_lhs Error: my_rowptr array too small for south data, num_row_s +1, size = ', num_row_s+1, size(my_rowptr)
          endif
          do concurrent (i = 1:mygrid_size + 1)
             my_rowptr(i) = rowptr_s(i)
          enddo
          !nnz_south set above - but double check
          if (nnz_south /= my_rowptr(mygrid_size + 1) - 1) then
             write(*,*) 'Dist_construct_lhs Error: mpi_rank, nnz_south = ', mpi_rank, nnz_south
          endif          

          if (nnz_south > size(my_colind) .or. nnz_south > size(my_values)) then
             write(*,*) 'Dist_construct_lhs Error: my_colind or my_values array too small for south data'
          endif
          
          do concurrent (i = 1:nnz_south)
             my_colind(i) = colind_s(i)
             my_values(i) = values_s(i)
          enddo
          
       elseif (ex_mpi_rank >=0) then !extra procs
          ! we received data from partner
          if (my_recvgrid_size /= mygrid_size) then
             write(*,*) 'Dist_cons=truct_lhs error: un_mpi_rank, my_recvgrid_size, mygrid_size = ', un_mpi_rank, my_recvgrid_size, mygrid_size
          endif
             
          do concurrent (i = 1: mygrid_size + 1)
             my_rowptr(i)  = partner_rowptr(i)
          enddo

          nnz_north = partner_rowptr(mygrid_size + 1) - 1
             
          do concurrent (i = 1:nnz_north)
             my_colind(i) = partner_cols(i)
             my_values(i) = partner_values(i)
          enddo
       endif !end extra procs
    else !one task (un_mpi_task = 1)
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
         un_mpi_rank, un_mpi_size, ex_mpi_rank, ex_mpi_size, &
         mlon0, mlon1, mlond0, mlond1, &
         lat_rank, lon_rank, &
         ij_start_s, ij_stop_s, ij_start_n, ij_stop_n, &
         calc_grid_ij, mpi_sendnorth_vec, mygrid_size, &
         gather_lon_1d, mpi_size, my_sendgrid_size, my_recvgrid_size

    real(kind=rp),dimension(mlatd0:mlatd1,mlond0:mlond1),intent(in) :: coef_10_s, coef_10_n
    real(kind=rp),dimension(mygrid_size) :: rhs
    real(kind=rp),dimension(ij_start_s:ij_stop_s) :: rhs_s
    real(kind=rp),dimension(ij_start_n:ij_stop_n) :: rhs_n

    integer :: i,j,ij, jN, j_start, cnt, istart

    real(kind=rp),dimension(nmlon) :: coef10_j1_buf

    real(kind=rp),dimension(my_sendgrid_size) :: sendbuf
    real(kind=rp),dimension(my_recvgrid_size) :: recvbuf

    
    rhs = 0.0

    if (mpi_rank >= 0) then !dynamo procs
    
       rhs_n = 0.0
       rhs_s = 0.0

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
                write(*,*) 'Dist_construct_rhs Error: ij south index for rhs', mpi_rank, ij
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
                write(*,*) 'Dist_construct_rhs Error: ij north index for rhs', mpi_rank, ij
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
                write(*,*) ' Dist_construct_rhs Error: ij south index 2 for rhs: rank, ij = ', mpi_rank, ij
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
                write(*,*) 'Dist_construct_rhs Error: ij north index 2 for rhs: rank, i,j,jN, ij =', mpi_rank, i,j,jN,ij
             endif
             rhs_n(ij) = coef_10_n(j,i)
          enddo
       enddo
    endif !dynamo procs only
    
    !redistribution with partner (i.e., send north hemisphere to extra procs)    
    if (un_mpi_size > 1) then
       !now do partner hemisphere exchange for continguous rows on incr. procs

       if (mpi_rank >= 0) then !dynamo procs
          !copy rhs_n to sendbuf
          cnt = 0
          do ij = ij_start_n, ij_stop_n
             cnt = cnt + 1
             sendbuf(cnt) = rhs_n(ij) 
          enddo
       endif !dynamo procs only
       
       call mpi_sendnorth_vec(sendbuf, recvbuf)

       if (mpi_rank >=0) then !dynamo procs
          !copy rhs_s into rhs
          cnt = 0
          do ij = ij_start_s, ij_stop_s
             cnt = cnt + 1
             rhs(cnt) = rhs_s(ij)
             if (isnan(rhs(cnt))) write(*,*) 'Dist_construct_rhs ERROR: rhs(cnt) is NaN, cnt = ', cnt, ' rank = ', mpi_rank
          enddo
       elseif (ex_mpi_rank >=0) then
          do i = 1, mygrid_size
             rhs(i) = recvbuf(i)
             if (isnan(rhs(i))) write(*,*) 'Dist_construct_rhs error: #2 rhs(i) is NaN, i = ', i
          enddo
       endif
    else !un_mpi_size = 1
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

    use mpi_module,only: lat_size,lon_size,union_world,&
         task_csr_rowstarts, un_mpi_rank, &
         un_mpi_size
    
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
    logical :: A_changed
    !most superlu structures are modulue-level variable
         
    ! Other variables
    integer(kind=c_int) :: info, ierr
    real(kind=c_double), target :: berr_array(nrhs)

    !get my first row in distributed matrix
    first_row = task_csr_rowstarts(un_mpi_rank)     !these are 0-based already 
    if (first_row < 0 .or. first_row >= n_global) then
       write(*,*)  "Superlu Error: first_row out of range", first_row, n_global
    endif

    !This is to detect if rowptr if colind changed
    if (.not. pattern_initialized) then
       g_pattern_hash = compute_pattern_hash(rowptr, colind)
       pattern_initialized = .true.
    endif
    
    if (.not. superlu_initialized) then

       !print *, 'initializing superlu ...'

       !save for error checking on further iteratons
       g_n_loc = n_loc
       g_nnz_loc = nnz_loc
       g_n_global = n_global
       g_first_row = first_row
       
       ! Create Fortran handles for the C structures used in SuperLU_DIST
       call f_create_gridinfo_handle(grid)
       call f_create_options_handle(options)
       call f_dcreate_ScalePerm_handle(ScalePermstruct)
       call f_dcreate_LUstruct_handle(LUstruct)
       call f_dcreate_SOLVEstruct_handle(SOLVEstruct)
       call f_create_SuperMatrix_handle(A)
       call f_create_SuperLUStat_handle(stat)

       ! Initialize the SuperLU_DIST process grid
       !i'll use the same layout as the dynamo, but double the lat for the extra procs
       nprow = 2*lat_size
       npcol = lon_size

       call f_superlu_gridinit(union_world, nprow, npcol, grid)
       if (nprow * npcol /= un_mpi_size) then
          write(*,*) "SuperLU ERROR: nprow*npcol != nprocs"
       endif
    
       !some debugging
       if (rowptr(n_loc+1) /= nnz_loc) then
          write(*,*)  "SuperLU ERROR rowptr(n_loc+1) /= nnz_loc ", rowptr(n_loc+1), nnz_loc
       endif

       if (minval(colind) < 0 .or. maxval(colind) >= n_global) then
          write(*,*)  "SuperLU Error colind out of bounds: min(col) max(col), n_global" , minval(colind), maxval(colind), n_global
       endif
       
       ! Set the default input options
       call f_set_default_options(options)
       !do the factorization from scratch (do not assume values are similar to previous)
       call set_superlu_options(options, Fact = DOFACT)

       ! Change one or more options
       !these below are the defaults
       call set_superlu_options(options,ColPerm=MMD_AT_PLUS_A)
       call set_superlu_options(options,RowPerm=LargeDiag_MC64)
       
       !for IterRefine SLU_DOUBLE=2 9want double precision for refinement)
       call set_superlu_options(options, IterRefine = 2)
       call set_superlu_options(options, Equil=1)
       call set_superlu_options(options, ReplaceTinyPivot=1) !for more stability
       !you can turn off print stat (with 0) for less solver output
       call set_superlu_options(options, PrintStat=1)
       !
       
       ! Initialize ScalePermstruct and LUstruct
       call f_dScalePermstructInit(n_global, n_global, ScalePermstruct)
       call f_dLUstructInit(n_global, n_global, LUstruct)
       
       !create the distributed compressed row matrix A (O-index)
       call f_dCreate_CompRowLoc_Mat_dist(A, n_global, n_global, nnz_loc, n_loc, first_row, &
            values, colind, rowptr, SLU_NR_loc, SLU_D, SLU_GE) 

       ! Initialize the statistics variables
       call f_PStatInit(stat)
    
       
       superlu_initialized = .true.
       
    else
       ! reuse symbolic structure
       ! only values should have changed since previous iteration (otherwise need DOFACT)
       ! pointer should be the same (so let's do a check here)

       !error check 
       A_changed = .false.
       
       if (g_n_loc /= n_loc .or. g_nnz_loc /= nnz_loc .or. &
            g_n_global /= n_global .or. g_first_row /= first_row) then
          write (*, *) 'SUPERLU ERROR: matrix structure has changed: forcing DOFACT!'
          A_changed = .true.
       elseif (compute_pattern_hash(rowptr, colind) /= g_pattern_hash) then
          write(*,*) 'SUPERLU: sparsity pattern changed, forcing DOFACT'
          A_changed = .true.
       endif

       if (A_changed ) then
          ! don't believe this will ever be triggered in the current code,
          ! but just to be safe :)
          g_pattern_hash = compute_pattern_hash(rowptr, colind)
          g_n_loc = n_loc
          g_nnz_loc = nnz_loc
          g_n_global = n_global
          g_first_row = first_row
          
          ! Destroy old symbolic data
          call f_Destroy_CompRowLoc_Mat_dist(A)
          call f_dScalePermstructFree(ScalePermstruct)
          call f_dLUstructFree(LUstruct)
          
          ! Reinitialize symbolic containers
          call f_dScalePermstructInit(n_global, n_global, ScalePermstruct)
          call f_dLUstructInit(n_global, n_global, LUstruct)

          ! Recreate A with new sparsity pattern
          call f_dCreate_CompRowLoc_Mat_dist(A, n_global, n_global, nnz_loc, n_loc, first_row, &
               values, colind, rowptr, SLU_NR_loc, SLU_D, SLU_GE)

          call set_superlu_options(options, Fact = DOFACT)

       else
          !we do not need to refactor (COlPerm and ROwPerm stays the same)
          call set_superlu_options(options, Fact = SamePattern_SameRowPerm)
       endif
    endif
    
    ! Setup the right hand side (rhs contains local data)
    sol=rhs ! Copy RHS to solution vector

   
    ! Call the linear equation solver (writes over rhs (sol))
    call f_pdgssvx(options, A, ScalePermstruct, sol, n_loc, nrhs, &
         grid, LUstruct, SOLVEstruct, berr_array, stat, info)

    if (info /= 0) then
       write(*,*) 'SuperLU ERROR: pdgssvx failed with mpi_rank, INFO = ', un_mpi_rank, info
    endif
    if (info == 0 .and. un_mpi_rank == 0) then
       write(*,*) 'Success: SuperLU Backward error: ', berr_array(1)
    endif

    ! result is sol (already assigned by reference in pdgssvx)

    
  endfunction dist_solve_superlu
  !-----------------------------------------------------------------------

  subroutine finalize_superlu()

    if (superlu_initialized) then

       !release storage allocated by superlu
       call f_PStatFree(stat)  !freed after each solve)
       call f_Destroy_CompRowLoc_Mat_dist(A)
       call f_dScalePermstructFree(ScalePermstruct)
       call f_dDestroy_LU_SOLVE_struct(options, g_n_global, grid, LUstruct, SOLVEstruct)

       ! Release the SuperLU process grid
       call f_superlu_gridexit(grid)

       !free the fortran handle wrappers
       call f_destroy_gridinfo_handle(grid)
       call f_destroy_options_handle(options)
       call f_destroy_ScalePerm_handle(ScalePermstruct)
       call f_destroy_LUstruct_handle(LUstruct)
       call f_destroy_SOLVEstruct_handle(SOLVEstruct)
       call f_destroy_SuperMatrix_handle(A)
       call f_destroy_SuperLUStat_handle(stat)
       
       superlu_initialized = .false.
    endif

  end subroutine finalize_superlu

!-----------------------------------------------------------------------

function compute_pattern_hash(rowptr, colind) result(h)
  use iso_c_binding, only: c_int, c_int64_t
  implicit none

  integer(kind=c_int), intent(in) :: rowptr(:)
  integer(kind=c_int), intent(in) :: colind(:)
  integer(kind=c_int64_t) :: h

  integer :: i
  integer(kind=c_int64_t), parameter :: FNV_OFFSET = &
       1469598103934665603_c_int64_t
  integer(kind=c_int64_t), parameter :: FNV_PRIME  = &
       1099511628211_c_int64_t

  ! Initialize hash
  h = FNV_OFFSET

  ! Hash rowptr (row structure)
  do i = 1, size(rowptr)
     h = ieor(h, int(rowptr(i), c_int64_t))
     h = h * FNV_PRIME
  end do

  ! Hash colind (column pattern + ordering)
  do i = 1, size(colind)
     h = ieor(h, int(colind(i), c_int64_t))
     h = h * FNV_PRIME
  end do

end function compute_pattern_hash
  
!-----------------------------------------------------------------------
  function dist_flatten(fin) result(fout)
! reorder 2D fields (lat-lon) into 1D vector (RHS)
! northern/southern hemispheres are either separate or averaged
! based on their latitude ranges (high-lat, transition, low-lat, equator)
    ! give north hemispher to extra processes
    
    use params_module,only:nmlat_h,nmlat_T1,nmlon
    use cons_module,only:jlatm_JT
    use mpi_module, only:lat_rank, mlat0, mlat1, &
         mlon0, mlon1, mlatd0, mlatd1, mlond0, mlond1, &
         ij_start_n, ij_stop_n, mpi_size, mpi_rank, &
         ij_start_s, ij_stop_s, mygrid_size, &
         calc_grid_ij, my_recvgrid_size, mpi_sendnorth_vec, &
         my_sendgrid_size, ex_mpi_rank, un_mpi_rank, un_mpi_size

    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: fin
    real(kind=rp),dimension(mygrid_size) :: fout
    
    real(kind=rp),dimension(ij_start_s:ij_stop_s) :: fout_s
    real(kind=rp),dimension(ij_start_n:ij_stop_n) :: fout_n
    integer :: i,j,ij, loop_start_j, loop_stop_j, jS, jN, cnt
    integer :: istart
    real(kind=rp) :: avg

    real(kind=rp),dimension(my_sendgrid_size) :: sendbuf
    real(kind=rp),dimension(my_recvgrid_size) :: recvbuf

    fout = 0.0

    if (mpi_rank >=0) then !dynamo_procs
       ! Initialize arrays to avoid uninitialized values
       fout_s = 0.0
       fout_n = 0.0
       sendbuf = 0.0
    
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
                if (isnan(fout_s(ij))) write(*,*) 'Flatten warning: fout_s(ij) is NaN, ij, i, jS = ', ij, i, jS
             else
                write(iulog,*) "Error in flatten south, mpirank, ij = ", mpi_rank, ij
             endif
             
             !north
             jN = nmlat_T1-j+1
             ij = calc_grid_ij(i,jN,lat_rank)
             !checkbounds
             if (ij >= ij_start_n .and. ij <= ij_stop_n) then
                fout_n(ij) = fin(2,j,i)
                if (isnan(fout_n(ij))) write(*,*) 'Flatten warning: fout_n(ij) is NaN, ij, i, jN = ', ij, i, jN
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
                if (isnan(fout_s(ij))) write(*,*) 'Flatten warning: fout_s(ij) is NaN, ij, i, jS = ', ij, i, jS
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
                if (isnan(fout_n(ij))) write(*,*) 'Flatten warning: fout_n(ij) is NaN, ij, i, jN = ', ij, i, jN

             else
                write(iulog,*) "Error in flatten 2 north, mpirank, ij = ", mpi_rank, ij
             endif

          enddo
       endif
    endif !dynamo procs

    !!now give north to extra procs
    if (un_mpi_size == 1) then !no extra procs
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
       !write(*,*) 'AB: flatten: mygrid_size,  my_sendgrid_size, my_recvgrid_size ', mygrid_size,  my_sendgrid_size, my_recvgrid_size 
       if (mpi_rank >=0) then !dynamo procs
          !copy north to sendbuf
          cnt = 0
          do ij = ij_start_n, ij_stop_n
             cnt = cnt + 1
             sendbuf(cnt) = fout_n(ij)
          enddo
       endif !dynam procs only
       
       call mpi_sendnorth_vec(sendbuf, recvbuf)

       !extra procs 
       if (mpi_rank >=0) then !dynamo procs
          !copy south into fout
          cnt = 0
          do ij = ij_start_s, ij_stop_s
             cnt = cnt + 1
             fout(cnt) = fout_s(ij)
             if (isnan(fout(cnt))) write(*,*) 'Flatten warning: #1 fout(cnt) is NaN, cnt, ij = ', cnt, ij
          enddo
       elseif (ex_mpi_rank >=0) then
          do i = 1, mygrid_size
             fout(i) = recvbuf(i)
             if (isnan(fout(i))) write(*,*) 'Flatten warning: #2 fout(i) is NaN, i = ', i
          enddo
       endif
    endif ! multiple procs
  endfunction dist_flatten

!-----------------------------------------------------------------------
  function dist_unravel(fin) result(fout)
    ! reorder 1D vector (RHS) into 2D fields (lat-lon)
    ! 1d vector has the distributed north process rows, so we have to undo that also
    ! this is the reverse of the dist_flatten
    ! this needs all union procs
    
    use params_module,only:nmlat_h,nmlat_T1,nmlon
    use mpi_module, only:lat_rank, mlat0, mlat1, &
         mlon0, mlon1, &
         un_mpi_rank, un_mpi_size, &
         ij_start_n, ij_stop_n, &
         ij_start_s, ij_stop_s, lat_rank,  &
         mygrid_size, mpi_rank, calc_grid_ij, mpi_size, &
         mpi_recvnorth_vec, my_sendgrid_size, my_recvgrid_size

    real(kind=rp),dimension(mygrid_size),intent(in) :: fin
    real(kind=rp),dimension(2,mlat0:mlat1,mlon0:mlon1) :: fout

    integer :: i,j,ij, cnt, jS, jN
    real(kind=rp),dimension(ij_start_s:ij_stop_s) :: fin_s
    real(kind=rp),dimension(ij_start_n:ij_stop_n) :: fin_n

    !purposely reversed for the recv cal
    real(kind=rp),dimension(my_recvgrid_size) :: sendbuf
    real(kind=rp),dimension(my_sendgrid_size) :: recvbuf
    
    !output 3D array for north=2 & south=1
    fout = 0.0

    if (un_mpi_size == 1) then
       cnt = 0
       !south
       do ij = ij_start_s, ij_stop_s
          cnt = cnt + 1
          fin_s(ij) = fin(cnt)
       enddo
       !north
       do ij = ij_start_n, ij_stop_n
          cnt = cnt + 1
          fin_n(ij) = fin(cnt)
       enddo
       !fill 3D array
       do concurrent (i = mlon0:mlon1, j = mlat0:mlat1)
          jS = j
          ij =  calc_grid_ij(i,jS,lat_rank)
          if (ij >= ij_start_s .and. ij <= ij_stop_s) then
             fout(1,j,i) = fin_s(ij)
             if (j == nmlat_h) then !equator, set in north also
                fout(2,j,i) = fin_s(ij)
             endif
          else
             write(*,*) "Error in unravel south (1 proc case), ij, ij_start_s, ij_stop_s = ", j, ij_start_s, ij_stop_s
          endif
          if (j == nmlat_h) then
             cycle
          endif
          jN = nmlat_T1-j+1 !north j
          ij =  calc_grid_ij(i,jN,lat_rank)
          if (ij >= ij_start_n .and. ij <= ij_stop_n) then
             fout(2,j,i) = fin_n(ij)
          else
             write(*,*) "Error in unravel north (1 proc case), ij = ",ij, ij_start_n, ij_stop_n
          endif
      enddo    
    else ! multiple procs
       if (mpi_rank >=0 ) then !dynamo procs - own south and get north from extra
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
                write(*,*) "Error in unravel south 1, mpirank, ij, ij_start, ij_stop = ", mpi_rank, ij, ij_start_s, ij_stop_s
             endif
          enddo
       else !now extra procs: copy north to buffer to send
          do i=1, my_recvgrid_size
             sendbuf(i) = fin(i)
          enddo
       endif !extra procs

       !get north from extra procs (all call)
       call mpi_recvnorth_vec(sendbuf, recvbuf)

       if (mpi_rank >= 0) then
          !dynamo procs copy recvbuf into fin_n
          cnt = 0
          do ij = ij_start_n, ij_stop_n
             cnt = cnt + 1
             fin_n(ij) = recvbuf(cnt)
          enddo

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
                write(*,*) "Error in unravel north 1, mpirank, ij = ", mpi_rank, ij, ij_start_n, ij_stop_n
             endif
          enddo
       endif
    endif !multiprocs


  endfunction dist_unravel

  !-----------------------------------------------------------------------

  subroutine dist_solver_final()

    !clean up
    call dist_spmv_free(halo)

    if (allocated(g_rowptr)) deallocate(g_rowptr)
    if (allocated(g_colind)) deallocate(g_colind) 
    if (allocated(g_values_csr)) deallocate(g_values_csr)
     
    call finalize_superlu()

  endsubroutine dist_solver_final
  
  !-----------------------------------------------------------------------

  subroutine insert_sort(array_i, array_r, len)
   ! this is only an ok sorting approach for small arrays
   ! since its O(n^2)
    ! (only row one of matrix sort will be a little longer than MAX_NNZ
    !, but is already mostly sorted, so it's fine)
   ! sorting by the int array but moving the reals 
   ! then removes zeros
   
     integer, dimension(:), intent(inout) :: array_i     ! Assumed-shape
     real(kind=rp), dimension(:), intent(inout) :: array_r  ! Assumed-shape
     integer, intent(inout) :: len
     integer :: i, j, temp_i, z
     real(kind=rp) :: temp_r
     integer :: keep(size(array_i))  ! Allocate based on actual size
     
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


subroutine write_dist_to_file( filename, phase, mlatd0,mlatd1,mlond0,mlond1, nnz, colind, rowptr, values_csr, pot_hl_f, z, rhs, fac_hl, pot, sol, fileid, Xid, potid)

! phase 0 is initialization, and before solver
! phase 1 is after solve - must call both in order!  
    use iso_c_binding
    use mpi
    use params_module, only:nmlat_h,nmlat_T1,nmlon
    use mpi_module, only: mpi_rank,dynamo_world,union_world, lat_rank,lon_rank,&
         nmlat_task,nmlon_task,mygrid_size,&
         mlat0, mlat1, mlon0, mlon1, &
         task_csr_rowstarts, mpi_size, &
         un_mpi_size, un_mpi_rank
    
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: phase, nnz
    integer,intent(in) :: mlatd0,mlatd1,mlond0,mlond1

    real(kind=rp), intent(in) :: z(:), rhs(:), pot_hl_f(:), &
         fac_hl(:,:,:), pot(:,:,:), values_csr(:), sol(:)
    integer, intent(in) :: colind(:), rowptr(:)
    integer, intent (inout) :: fileid, Xid, potid
    
    integer :: ierr, unitno
    integer :: nlonlat
    integer :: dd(4)
    integer :: mygridsize_id, mygridsizep1_id, nmlatt1_id, nmlon_id, nmlath_id, &
         size_id, sizep1_id, &
         nnz_id, nprocsp1_id, hemsize_id, &
         mylatsize_id, mylonsize_id

    integer :: fachlid, potflatid, zid, rhsBid, valuesid, colsid, rowptrid, &
         procrowstartsid, countB(3), startB(3)

    ! global matrix size
    nlonlat = nmlat_T1*nmlon 

    if (phase == 0) then
       write(*,*) 'creating nc-file:', filename
       
       ierr = nf_create(filename, NF_CLOBBER, fileid)
       
       !Define the dimensions
       ierr = nf_def_dim(fileid,'mygrid_size', mygrid_size, mygridsize_id)
       ierr = nf_def_dim(fileid, 'mygrid_sizep1', mygrid_size+1, mygridsizep1_id)
       ierr = nf_def_dim(fileid, 'mynnz', nnz, nnz_id)
       ierr = nf_def_dim(fileid,'nmlat_T1', nmlat_T1, nmlatt1_id)
       ierr = nf_def_dim(fileid,'nmlat_h', nmlat_h, nmlath_id)
       ierr = nf_def_dim(fileid,'nmlon', nmlon, nmlon_id)
       ierr = nf_def_dim(fileid,'allprocsp1', un_mpi_size+1, nprocsp1_id)
       ierr = nf_def_dim(fileid,'global_size', nlonlat, size_id)
       ierr = nf_def_dim(fileid,'global_sizep1', nlonlat + 1, sizep1_id)
       ierr = nf_def_dim(fileid,'hem', 2, hemsize_id)
       ierr = nf_def_dim(fileid,'mylat_size', mlatd1-mlatd0+1, mylatsize_id)
       ierr = nf_def_dim(fileid,'mylon_size',mlond1-mlond0+1, mylonsize_id)

       ! define vars and their associated dimensions         
       dd(1) = mygridsize_id
       ierr = nf_def_var(fileid, 'z', NF_DOUBLE, 1, dd, zid)
       
       dd(1) = mygridsize_id
       ierr = nf_def_var(fileid, 'pot_in_flat', NF_DOUBLE, 1, dd, potflatid)
       
       dd(1) = mygridsize_id
       ierr = nf_def_var(fileid, 'rhsB', NF_DOUBLE, 1, dd, rhsBid)

       dd(1) = nnz_id
       ierr = nf_def_var(fileid, 'valuesA', NF_DOUBLE, 1, dd, valuesid)

       dd(1) = nnz_id
       ierr = nf_def_var(fileid, 'colindxA',NF_INT, 1, dd, colsid)

       dd(1) = mygridsizep1_id
       ierr = nf_def_var(fileid, 'rowptrA',NF_INT, 1, dd, rowptrid)

       dd(1) = nprocsp1_id
       ierr = nf_def_var(fileid, 'procRowStarts',NF_INT, 1, dd, procrowstartsid)

       dd(1) = hemsize_id
       dd(2) = mylatsize_id
       dd(3) = mylonsize_id
       ierr = nf_def_var(fileid, 'fac_hl',NF_DOUBLE, 3, dd, fachlid)
       
       dd(1) = mygridsize_id
       ierr = nf_def_var(fileid, 'X', NF_DOUBLE, 1, dd, Xid)
         
       dd(1) = hemsize_id
       dd(2) = mylatsize_id
       dd(3) = mylonsize_id
       ierr = nf_def_var(fileid, 'pot',NF_DOUBLE, 3, dd, potid)

       !now fill pre-solve fields (everything but sol and pot)
       !change mode
       ierr=nf_enddef(fileid)

       !Output the flat potential
       startB(1)=1
       countB(1)=mygrid_size
       ierr=NF_PUT_VARA_DOUBLE(fileid,potflatid,startB,countB,pot_hl_f)
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
       !row starts
       startB(1)=1
       countB(1)=un_mpi_size+1
       ierr=NF_PUT_VARA_INT(fileid,procrowstartsid,startB,countB,task_csr_rowstarts)
       !fac_hl
       startB(1)=1
       countB(1)=2
       startB(2)=1
       countB(2)=mlatd1-mlatd0+1
       startB(3)=1
       countB(3)=mlond1-mlond0+1
       ierr=NF_PUT_VARA_DOUBLE(fileid,fachlid,startB,countB,fac_hl)
       
       !don't close
         
    elseif (phase ==1) then

       startB(1)=1
       countB(1)=mygrid_size
        ierr=NF_PUT_VARA_DOUBLE(fileid,Xid,startB,countB,sol)

        startB(1)=1
        countB(1)=2
        startB(2)=1
        countB(2)=mlatd1-mlatd0+1
        startB(3)=1
        countB(3)=mlond1-mlond0+1
        ierr=NF_PUT_VARA_DOUBLE(fileid,potid,startB,countB,pot)

        !close file
        ierr=NF_CLOSE(fileid)
         
    endif
    
  end subroutine write_dist_to_file


  
endmodule dist_solver_module
