module dist_solver_module
  use perf_mod, only: t_startf, t_stopf

  use prec,only:rp

  implicit none

  contains
!-----------------------------------------------------------------------
  subroutine dist_linear_system(mlatd0,mlatd1,mlond0,mlond1, &
    bij,pot_hl,fac_hl,coef_ns,pot)
! construct linear system based on bij and coef and solve in pot
    
! if FAC is read in, pot_hl is not used, only fac_hl is used
! if potential is read in, pot_hl is used, fac_hl is output

    
    use params_module,only:nmlat_h,nmlat_T1,nmlon
    use cons_module,only:read_fac
    use mpi_module,only:mpi_rank,dynamo_world,lat_rank,lon_rank,nmlat_task,nmlon_task


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
    integer :: nlonlat,mlat0,mlat1,mlon0,mlon1,i,j,isn,ic,nnz

    !integer,dimension(nmlat_T1*nmlon+1) :: rowptr,colptr
    !integer,dimension(12*nmlat_T1*nmlon) :: colind,rowind
    !real(kind=rp),dimension(12*nmlat_T1*nmlon) :: values_csr,values_csc
    !real(kind=rp),dimension(nmlat_h,nmlon) :: bij_full
    !real(kind=rp),dimension(2,nmlat_h,nmlon) :: pot_hl_full,fac_hl_    !real(kind=rp),dimension(10,2,nmlat_h,nmlon) :: coef_ns_full
    !real(kind=rp),dimension(10,nmlat_T1,nmlon) :: coef_full
    !real(kind=rp),dimension(2,nmlat_h,0:nmlon+1) :: fac_hl_2,pot_2
    !real(kind=rp),dimension(nmlat_T1*nmlon) :: rhs,z,pot_hl_f,sol

    integer,dimension(:), allocatable :: rowptr
    integer,dimension(:), allocatable :: colind
    real(kind=rp),dimension(::), allocatable :: values_csr
    real(kind=rp),dimension(10,mlatd0:mlatd1,mlond0:mlond1) :: coef_s
    real(kind=rp),dimension(10,mlatd0:mlatd1,mlond0:mlond1) :: coef_n

    integer :: ier
    integer :: mygrid_size
    
    call t_startf('dist_linear_system')

    ! global matrix size
    nlonlat = nmlat_T1*nmlon 

    !proc domain extents without halo pts
    mlat0 = mlatd0+1
    mlat1 = mlatd1-1
    mlon0 = mlond0+1
    mlon1 = mlond1-1

    !number of grid points I own (no halo)
    !these should be the same
    !mygrid_size = (mlat1-mlat0+1)*(mlon1-mlon0+1)
    mygrid_size = nmlat_task(lat_rank)*nmlon_task(lon_rank)

    !allocate space for rowptr,colind,values_csr
    !double for both hemispheres
    allocate(rowptr(2*mygrid_size+1))
    allocate(colind(12*2*mygrid_size))
    allocate(values_csr(12*2*mygrid_size))
    
! for now, split two hemispheres (keep halo pts)
    do concurrent (i = mlond0:mlond1, j = mlatd0:mlatd1, ic = 1:10)
       coef_s(ic,j,i) = coef_ns(ic,1,j,i) !south
       coef_n(ic,j,i) = coef_ns(ic,2,j,i) ! north (note: row = nmlat_T1-j+1) 
    enddo

! construct LHS matrix in CSR format
    call dist_construct_lhs(mygrid_size,bij,coef_s(1:9,:,:),coef_n(1:9,:,:),rowptr,colind,values_csr)
    nnz = rowptr(nlonlat+1)-1

! RHS is dense
    rhs = construct_rhs(coef_full(10,:,:))

! determine FAC forcing (dense)
    if (read_fac) then ! input is corrected fac_hl, pot_hl is not used
       z = flatten(fac_hl_full)

    else ! input is pot_hl, fac_hl is to be calculated (output)

! A. Maute 2023/11/21: put the high latitude potential in X
! and then use LHS to calculate the RHS FAC
       pot_hl_f = flatten(pot_hl_full)

! z = matmul(lhs, pot_hl)
       z = 0
       do i = 1,nlonlat
          do j = rowptr(i),rowptr(i+1)-1
             z(i) = z(i)+values_csr(j)*pot_hl_f(colind(j))
          enddo
       enddo

! no need for correction since it is from the divergence of horizontal current

! reconstruct 2D distribution of FAC based on z
       fac_hl_2(:,:,1:nmlon) = unravel(z)

! add periodic points
        do j = 1,nmlat_h
          do isn = 1,2
            fac_hl_2(isn,j,0) = fac_hl_2(isn,j,nmlon)
            fac_hl_2(isn,j,nmlon+1) = fac_hl_2(isn,j,1)
          enddo
        enddo
     endif !FAC

! add FAC forcing to RHS
     do i = 1,nlonlat
        rhs(i) = rhs(i)+z(i)
     enddo

!superlu
     call csr_to_csc(nlonlat,nlonlat,nnz, &
           rowptr,colind(1:nnz),values_csr(1:nnz), &
           colptr,rowind(1:nnz),values_csc(1:nnz))
      
     call t_startf('linear_system->solve_superlu')
     sol = solve_superlu(nlonlat,nnz,colptr,rowind(1:nnz),values_csc(1:nnz),rhs)
     call t_stopf('linear_system->solve_superlu')

! reconstruct 2D distribution of potential based on the solution
     pot_2(:,:,1:nmlon) = unravel(sol)

! periodic points
     do j = 1,nmlat_h
        do isn = 1,2
          pot_2(isn,j,0) = pot_2(isn,j,nmlon)
          pot_2(isn,j,nmlon+1) = pot_2(isn,j,1)
        enddo
     enddo
  

     call mpi_barrier (dynamo_world, ier)

    call t_startf('linear_system->bcast_3d')
    if (.not. read_fac) then
       call bcast_3d(fac_hl_2,2,nmlat_h,nmlon+2,root)

       do concurrent (i = mlond0:mlond1, j = mlatd0:mlatd1, isn = 1:2, j>=1 .and. j<=nmlat_h)
        fac_hl(isn,j,i) = fac_hl_2(isn,j,i)
      enddo
    endif

    call bcast_3d(pot_2,2,nmlat_h,nmlon+2,root)
    do concurrent (i = mlond0:mlond1, j = mlatd0:mlatd1, isn = 1:2, j>=1 .and. j<=nmlat_h)
      pot(isn,j,i) = pot_2(isn,j,i)
    enddo
    call t_stopf('linear_system->bcast_3d')

    call t_stopf('dist_linear_system')

  endsubroutine dist_linear_system
!-----------------------------------------------------------------------
  pure subroutine dist_construct_lhs(mygrid_size,bij,coef_s, coef_n,rowptr,colind,values)
! construct LHS matrix (CSR format - block row format for each task)

! need to set where the two hemispheres are connected
! This is not in the 9-point stencil but needs to be done manually

    use params_module,only:nmlat_h,nmlat_T1,nmlon
    use cons_module,only:jlatm_JT
    use mpi_module,only:mpi_rank,dynamo_world,lat_rank,lon_rank, &
         nmlat_task,nmlon_task, mlatd0, mlatd1, mlat0, mlat1, &
         lat_size, lon_size, task_lat_offset, ij_start_n, ij_stop_n, &
         ij_start_s, ij_stop_s
    
    integer,intent(in) :: mygrid_size
    real(kind=rp),dimension(mlatd0:mlatd1,mlond0:mlond1),intent(in) :: bij
    real(kind=rp),dimension(10,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: coef_s
    real(kind=rp),dimension(10,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: coef_n
    integer,dimension(2*mygrid_size),intent(out) :: rowptr
    integer,dimension(12*2*mygrid_size),intent(out) :: colind
    real(kind=rp),dimension(12*2*mygrid_size),intent(out) :: values

!mlatd0 and mlond0 can be 0    
!mlatd1 and mlond1 can be nmlat+1 and nmlon+1

    
! if two hemispheres are uncoupled at high latitudes, set bijSum to zero
    real(kind=rp),parameter :: bijSum = 0
    integer :: nlonlat,i,j,jS,jN,ij,isub,im,ip
    integer :: loop_start_i, loop_stop_i, loop_start_j, loop_stop_j
    integer,dimension(mygrid_size) :: rowcnt_s, rowcnt_n

    
! the first row has most elements (nmlon+2)
!TO DO - handle this special case (only 1 processor does this_
    integer,dimension(nmlon+2) :: jcol1
    real(kind=rp),dimension(nmlon+2) :: nzval1

    ! other rows have at most 12 elements
    !we want the indexes to be the global ranges of rows....
    integer,dimension(12,ij_start_s:ij_stop_s) :: jcol_s
    real(kind=rp),dimension(12,ij_start_s:ij_stop_s) :: nzval_s
    
    integer,dimension(12,ij_start_n:ij_stop_n) :: jcol_n
    real(kind=rp),dimension(12,ij_start_n:ij_stop_n) :: nzval_n 
    
    !global size (might not need this?)
    nlonlat = nmlat_T1*nmlon

    !initialize
    rowcnt_s = 0
    rowcnt_n = 0
    
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

    ! procs will set up rows for their domain in the southern hemisphere and nothern hemisphere
    ! we'll keep in two different arrays to combine when we assemble the CSR matricx
    ! each proc has to have a continguous block od rows, so we'll have to apply a permutation
    
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

       if (lon_rank == 0) then !I also own i=1 (special case)


          !TODO


          
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

          !old code for reference
          !jcol_s(rowcnt_s(ij)+1:rowcnt_s(ij)+2,ij) = (/(1-1)*nmlat_T1+j,(i-1)*nmlat_T1+j/)
          !nzval_s(rowcnt_s(ij)+1:rowcnt_s(ij)+2,ij) = (/-1,1/)
          !rowcnt_s(ij) = rowcnt_s(ij)+2
       enddo
       
       ! North pole: for each i, set Phi(i,nmlat_T1) = phi_pol
       loop_start_i = mlon0
       do i = loop_start_i, loop_stop_i   
          !calculate matrix row for North
          jN = nmlat_T1-j+1 !j=1, so jN= nmlat_T1
          ij = calc_grid_ij(i,jN,lat_rank)

          rowcnt_n(ij) = rowcnt_n(ij)+1
          !the column will be the same i position but at j=1 (instead of j=15)
          ! TO DO: verify this (doesn't make intuive sense to me)
          jcol_n(rowcnt(ij),ij) = calc_grid_ij(i,1,0) ! j=1, lat_rank=1
          nzval_n(rowcnt(ij),ij) = 1
          
          !old for reference
          !rowcnt(ij) = rowcnt(ij)+1
          !jcol(rowcnt(ij),ij) = (i-1)*nmlat_T1+nmlat_T1-j+1
          !nzval(rowcnt(ij),ij) = 1

          
       enddo
       
    endif !end of loop for procs owning j=1


    
    !REGION BETWEEN POLES AND latm_JT (i.e., 2:latm_JT-1)
    if (mlat0<latm_JT) then !I own latitudes in this region

       m = task_lat_offset[lat_rank] ! this is needed to calculate the row

       !loop through the longitudes in my grid
       do i = mlon0, mlon1
          !we have the halo regions
          im = i-1
          ip = 1+1

          !loop through latitudes
          if (lat_rank == 0) then !skip the pole
             loop_start_j = 2
          else
             loop_start_j = mlat0
          endif
          loop_end_j = min(mlat1, latm_JT-1)
          
          do j = loop_start_j, loop_end_j
             
             !South Hemishere
             jS=j
             jN = nmlat_T1-j+1
             ij = calc_grid_ij(i,jS,lat_rank)

             if (j==mlat0) then !bottom edge of proc domain
                
             elseif (i=mlat1) then !top edge of proc domain

                
             else !interior proc domain
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
                nzval_s(rowcnt_s(ij),ij)= bij(,j,i)
                
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
             endif ! interior

             !North Hemisphere
             ! note that coef has the direction switched in the NH
             jN = nmlat_T1-j+1
             ij = (nmlon*m) + (jN-m) + ((i-1)*nmlat_task(lat_rank))

             if (i==1) then !would be better for perf to have these if statements
                            ! outside the j loop (but for now, following old structure)

                
             elseif (i=nmlon) then

                
             else  

                
             endif   


             
          enddo ! end j loop through latitudes
       

       enddo !end i loop through longitudes

    endif ! end of REGION BETWEEN POLES AND latm_JT


    !---------------------------------------!   
          
    !latm_JT
    if (mlat0 <= latm_JT .and. mlat1 >= latm_JT) then ! I own latm_JT

       !TO DO (j=latmJT)
       
    endif
    ! between latm_JT and equator (nmlat_h)
    if (mlat0 > latm_JT .and. mlat0 /= nmlat_h ) .or. (mlat1 > latm_JT .and. mlat0 \= nmlat_h) then


    !TO DO
    
    endif
    !equator
    !we could check in the process row, but if we own to the equator then it will be mlat1
    if (mlat_1 == nmlat_h) then !I own equator

       ! TO DO
       
    endif
       

!!!OLD CODE
    ! are my cols in ascending order (may not be at the proc boundaries)
    ! permute so that H& S rows are contiguous
    
    rowptr(1) = 1
    do i = 2,nlonlat+1
      rowptr(i) = rowptr(i-1)+rowcnt(i-1)
    enddo

    do j = 1,nmlon+2
      colind(j) = jcol1(j)
      values(j) = nzval1(j)
    enddo

    do i = 2,nlonlat
      do j = 1,rowcnt(i)
        colind(rowptr(i)+j-1) = jcol(j,i)
        values(rowptr(i)+j-1) = nzval(j,i)
      enddo
    enddo

  endsubroutine dist_construct_lhs
!-----------------------------------------------------------------------
  pure function dist_construct_rhs(coef_10) result(rhs)
! construct vector RHS
! this is different from flatten

    use params_module,only:nmlat_h,nmlat_T1,nmlon
    use cons_module,only:phi_pol

    real(kind=rp),dimension(nmlat_T1,nmlon),intent(in) :: coef_10
    real(kind=rp),dimension(nmlat_T1*nmlon) :: rhs

    integer :: i,j,ij

    rhs = 0

! set up poles, there are no c6,c7,c8 values
    j = 1

! for longitude i=1 at the south pole
    i = 1
    ij = (i-1)*nmlat_T1+j
    rhs(ij) = sum(coef_10(j,:))

! for each i, set Phi(i,nmlat_T1) = phi_pol at the north pole
    do concurrent (i = 1:nmlon)
      ij = (i-1)*nmlat_T1+nmlat_T1-j+1
      rhs(ij) = phi_pol
    enddo

    do concurrent (i = 1:nmlon, j = 2:nmlat_h)
      ij = (i-1)*nmlat_T1+j
      rhs(ij) = coef_10(j,i)

      ij = (i-1)*nmlat_T1+nmlat_T1-j+1
      rhs(ij) = coef_10(nmlat_T1-j+1,i)
    enddo

  endfunction dist_construct_rhs

!-----------------------------------------------------------------------
  function solve_superlu_dist(n,nnz,colptr,rowind,values,rhs) result(sol)
    use iso_c_binding,only:c_int,c_long_long,c_double

    integer,intent(in) :: n,nnz
    integer(kind=c_int),dimension(n+1),intent(in) :: colptr
    integer(kind=c_int),dimension(nnz),intent(in) :: rowind
    real(kind=c_double),dimension(nnz),intent(in) :: values
    real(kind=rp),dimension(n),intent(in) :: rhs
    real(kind=rp),dimension(n) :: sol

! for SuperLU sparse matrix solver
    integer,parameter :: nrhs = 1
    integer :: i,iopt,info
    integer(kind=c_long_long) :: f_factors

    interface
      subroutine c_fortran_dgssv(iopt,n,nnz,nrhs, &
        values,rowind,colptr,b,ldb,f_factors,info) &
        bind(c,name='c_fortran_dgssv_')
        use iso_c_binding,only:c_int,c_long_long,c_double
        integer(kind=c_int) :: iopt,n,nnz,nrhs,ldb,info
        real(kind=c_double),dimension(nnz) :: values
        integer(kind=c_int),dimension(nnz) :: rowind
        integer(kind=c_int),dimension(n+1) :: colptr
        real(kind=c_double),dimension(ldb) :: b
        integer(kind=c_long_long) :: f_factors
      endsubroutine c_fortran_dgssv
    endinterface

    do concurrent (i = 1:n)
      sol(i) = rhs(i)
    enddo

! first, factorize the matrix, the factors are stored in *f_factors* handle
    iopt = 1
    call c_fortran_dgssv(iopt, n, nnz, nrhs, &
      values, rowind, colptr, sol, n, f_factors, info)
    write(6,"('INFO from LU decomposition = ',i4)") info

! second, solve the system using the existing factors
    iopt = 2
    call c_fortran_dgssv(iopt, n, nnz, nrhs, &
      values, rowind, colptr, sol, n, f_factors, info)
    write(6,"('INFO from triangular solve = ',i4)") info

! last, free the storage allocated inside SuperLU
    iopt = 3
    call c_fortran_dgssv(iopt, n, nnz, nrhs, &
      values, rowind, colptr, sol, n, f_factors, info)

  endfunction solve_superlu_dist

!-----------------------------------------------------------------------
  pure function flatten(fin) result(fout)
! reorder 2D fields (lat-lon) into 1D vector (RHS)
! northern/southern hemispheres are either separate or averaged
! based on their latitude ranges (high-lat, transition, low-lat, equator)

    use params_module,only:nmlat_h,nmlat_T1,nmlon
    use cons_module,only:jlatm_JT

    real(kind=rp),dimension(2,nmlat_h,nmlon),intent(in) :: fin
    real(kind=rp),dimension(nmlat_T1*nmlon) :: fout

    integer :: i,j,ij
    real(kind=rp) :: avg

! from pole to latm_JT, two hemispheres are uncoupled
    do concurrent (i = 1:nmlon, j = 1:jlatm_JT)
      ij = (i-1)*nmlat_T1+j
      fout(ij) = fin(1,j,i)

      ij = (i-1)*nmlat_T1+nmlat_T1-j+1
      fout(ij) = fin(2,j,i)
    enddo

! from latm_JT to equator, symmetric solution
    do concurrent (i = 1:nmlon, j = jlatm_JT+1:nmlat_h)
      avg = (fin(1,j,i)+fin(2,j,i))/2

      ij = (i-1)*nmlat_T1+j
      fout(ij) = avg

      ij = (i-1)*nmlat_T1+nmlat_T1-j+1
      fout(ij) = avg
    enddo

  endfunction flatten
!-----------------------------------------------------------------------
  pure function unravel(fin) result(fout)
! reorder 1D vector (RHS) into 2D fields (lat-lon)

    use params_module,only:nmlat_h,nmlat_T1,nmlon

    real(kind=rp),dimension(nmlat_T1*nmlon),intent(in) :: fin
    real(kind=rp),dimension(2,nmlat_h,nmlon) :: fout

    integer :: i,j,isn,ij

    do concurrent (i = 1:nmlon, j = 1:nmlat_h, isn = 1:2)
      if (isn == 1) then
        ij = (i-1)*nmlat_T1+j
      else
        ij = (i-1)*nmlat_T1+nmlat_T1-j+1
      endif
      fout(isn,j,i) = fin(ij)
    enddo

  endfunction unravel
  
   
  !-----------------------------------------------------------------------

  ! find the coef matrix row number for the grid location  i,j
  ! this works for north and south
  ! j can be in north or south
  ! note that my_latrank could be calculated from j, but it would be a
  ! lookup in mlat0_task() and the calling processor should know this info
  pure function calc_grid_ij(i,j,my_latrank) return(ij)
    
    use params_module,only:nmlon,nmlat_T1, nmlat_h
    use mpi_module, only:nmlat_task,task_lat_offset,lat_size,mlat0,mlat1


    integer, intent(in) :: i,j,my_latrank
    integer, intent(out):: ij
    
    integer:: m, my_numlat


     
  
    if (j <= nmlat_h) then !south hemi

       !if j is in ghost layer, then adjust my_latrank
       if (j == mlat1 + 1) then
          my_latrank  = my_latrank + 1
       elseif (j == mlat0 -1) then
          my_latrank  = my_latrank + 1
       endif
        
       !num of latitude points in proc parition    
       my_numlat = nmlat_task(my_latrank)
       
       !adjust if i is on edge of global domain
       if (i == 0) then
          i = nmlon
       elseif (i == nmlon+1)
          i = 1
       endif
       
       m = task_lat_offset(my_latrank)

    else !north hemi

       !if j is in ghost layer, then adjust my_latrank
       !TO DO - fix j
       if (j == mlat1 + 1) then
          my_latrank  = my_latrank + 1
       elseif (j == mlat0 -1) then
          my_latrank  = my_latrank + 1
       endif
         
       !num of latitude points in proc parition    
       my_numlat = nmlat_task(my_latrank)
       
       !adjust if i is on edge of global domain
       if (i == 0) then
          i = nmlon
       elseif (i == nmlon+1)
          i = 1
       endif
       
       !if i am in the task row that includes equator, then
       ! subtract 1 from num_lat
       if (my_latrank == lat_size - 1 ) then
          my_numlat = my_numlat - 1
       endif
       
       !m is different in the in the north hemisphere
       m = nmlat_T1 - task_lat_offset(my_latrank) - my_numlat
    endif
    
    ij = (nmlon*m) + (j-m) + ((i-1)*my_numlat)
    
  endfunction calc_grid_ij
  !-----------------------------------------------------------------------
  
endmodule dist_solver_module
