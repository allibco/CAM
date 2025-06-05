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
    !real(kind=rp),dimension(2,nmlat_h,nmlon) :: pot_hl_full,fac_hl_full
    !real(kind=rp),dimension(10,2,nmlat_h,nmlon) :: coef_ns_full
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
    use mpi_module,only:mpi_rank,dynamo_world,lat_rank,lon_rank,nmlat_task,nmlon_task

    integer,intent(in) :: mygrid_size
    real(kind=rp),dimension(mlatd0:mlatd1,mlond0:mlond1),intent(in) :: bij
    !real(kind=rp),dimension(nmlat_h,nmlon),intent(in) :: bij
    real(kind=rp),dimension(10,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: coef_s
    real(kind=rp),dimension(10,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: coef_n
!    real(kind=rp),dimension(9,nmlat_T1,nmlon),intent(in) :: coef
    !integer,dimension(nmlat_T1*nmlon+1),intent(out) :: rowptr
    !integer,dimension(12*nmlat_T1*nmlon),intent(out) :: colind
    !real(kind=rp),dimension(12*nmlat_T1*nmlon),intent(out) :: values
    integer,dimension(2*mygrid_size),intent(out) :: rowptr
    integer,dimension(12*2*mygrid_size),intent(out) :: colind
    real(kind=rp),dimension(12*2*mygrid_size),intent(out) :: values

    
! if two hemispheres are uncoupled at high latitudes, set bijSum to zero
    real(kind=rp),parameter :: bijSum = 0
    integer :: nlonlat,i,j,jS,jN,ij,isub,im,ip
    !integer,dimension(nmlat_T1*nmlon) :: rowcnt
    !TO DO: double?
    integer,dimension(mygrid_size) :: rowcnt

    
! the first row has most elements (nmlon+2)
!TO DO - handle this special case
    integer,dimension(nmlon+2) :: jcol1
    real(kind=rp),dimension(nmlon+2) :: nzval1

! other rows have at most 12 elements
    integer,dimension(12,mygrid_size) :: jcol
    real(kind=rp),dimension(12,mygrid_size) :: nzval

    !gloabl size (might not need this?)
    nlonlat = nmlat_T1*nmlon
    rowcnt = 0

!_______________    

!    I think one big loop over lat lines and then have to examin which
!    case we are in (pole or < latm_JT or otherwise)
! though that would mean lots of if statements, so ...figure out which lat_ranks have what?
    

 !
!Do the poles first (only affects subset of tasks)
    if (lat_rank == 0) then

       !POLES

    endif
! set up poles
    j = 1

! there are no c6,c7,c8 values at the south pole

! for longitude i=1 at the south pole
    i = 1
!    ij = (i-1)*nmlat_T1+j
!    lhs(ij,(i-1)*nmlat_T1+j           ) = sum(coef(9,j,:))-bijSum ! A. Richmond 2023/06/20: add bij effects
!    lhs(ij,(i-1)*nmlat_T1+j+1         ) = coef(3,j,i)             ! Sum_i=1^nmlon C3(i,1) Phi(i,2)
!    lhs(ij,(i-1)*nmlat_T1+nmlat_T1-j+1) = bijSum                  ! conjugate point (north pole at i=1)
    jcol1(1:3) = (/(i-1)*nmlat_T1+j,(i-1)*nmlat_T1+j+1,(i-1)*nmlat_T1+nmlat_T1-j+1/)
    nzval1(1:3) = (/sum(coef(9,j,:))-bijSum,coef(3,j,i),bijSum/)
    do isub = 2,nmlon ! Sum_i=1^nmlon C3(i,1) Phi(i,2)
!      lhs(ij,(isub-1)*nmlat_T1+j+1) = coef(3,j,isub)
      jcol1(isub+2) = (isub-1)*nmlat_T1+j+1
      nzval1(isub+2) = coef(3,j,isub)
    enddo
    rowcnt(1) = nmlon+2

! for other longitudes (i/=1) at the south pole
    do i = 2,nmlon
      ij = (i-1)*nmlat_T1+j
!      lhs(ij,(1-1)*nmlat_T1+j) = -1 ! Phi(i,1) = Phi(1,1)
!      lhs(ij,(i-1)*nmlat_T1+j) = 1
      jcol(rowcnt(ij)+1:rowcnt(ij)+2,ij) = (/(1-1)*nmlat_T1+j,(i-1)*nmlat_T1+j/)
      nzval(rowcnt(ij)+1:rowcnt(ij)+2,ij) = (/-1,1/)
      rowcnt(ij) = rowcnt(ij)+2
    enddo

! for each i, set Phi(i,nmlat_T1) = phi_pol at the north pole
    do i = 1,nmlon
      ij = (i-1)*nmlat_T1+nmlat_T1-j+1
!      lhs(ij,(i-1)*nmlat_T1+nmlat_T1-j+1) = 1
      rowcnt(ij) = rowcnt(ij)+1
      jcol(rowcnt(ij),ij) = (i-1)*nmlat_T1+nmlat_T1-j+1
      nzval(rowcnt(ij),ij) = 1
    enddo

    do i = 1,nmlon
      if (i == 1) then
        im = nmlon
      else
        im = i-1
      endif
      if (i == nmlon) then
        ip = 1
      else
        ip = i+1
      endif
!_______________    

! from pole to latm_JT, two hemispheres are uncoupled
      do j = 2,jlatm_JT-1
        jS = j
        jN = nmlat_T1-j+1

        ij = (i-1)*nmlat_T1+jS
!        lhs(ij,(im-1)*nmlat_T1+jS-1) = coef(6,jS,i)
!        lhs(ij,(im-1)*nmlat_T1+jS  ) = coef(5,jS,i)
!        lhs(ij,(im-1)*nmlat_T1+jS+1) = coef(4,jS,i)
!        lhs(ij, (i-1)*nmlat_T1+jS-1) = coef(7,jS,i)
!        lhs(ij, (i-1)*nmlat_T1+jS  ) = coef(9,jS,i)-bij(j,i) ! should be 1
!        lhs(ij, (i-1)*nmlat_T1+jS+1) = coef(3,jS,i)
!        lhs(ij, (i-1)*nmlat_T1+jN  ) = bij(j,i)              ! conjugate point, b(i,j) Phi*(i,j)
!        lhs(ij,(ip-1)*nmlat_T1+jS-1) = coef(8,jS,i)
!        lhs(ij,(ip-1)*nmlat_T1+jS  ) = coef(1,jS,i)
!        lhs(ij,(ip-1)*nmlat_T1+jS+1) = coef(2,jS,i)
        if (i == 1) then
          jcol(rowcnt(ij)+1:rowcnt(ij)+10,ij) = (/ &
            (i-1)*nmlat_T1+jS-1,(i-1)*nmlat_T1+jS,(i-1)*nmlat_T1+jS+1,(i-1)*nmlat_T1+jN, &
            (ip-1)*nmlat_T1+jS-1,(ip-1)*nmlat_T1+jS,(ip-1)*nmlat_T1+jS+1, &
            (im-1)*nmlat_T1+jS-1,(im-1)*nmlat_T1+jS,(im-1)*nmlat_T1+jS+1/)
          nzval(rowcnt(ij)+1:rowcnt(ij)+10,ij) = (/ &
            coef(7,jS,i),coef(9,jS,i)-bij(j,i),coef(3,jS,i),bij(j,i), &
            coef(8,jS,i),coef(1,jS,i),coef(2,jS,i), &
            coef(6,jS,i),coef(5,jS,i),coef(4,jS,i)/)
        elseif (i == nmlon) then
          jcol(rowcnt(ij)+1:rowcnt(ij)+10,ij) = (/ &
            (ip-1)*nmlat_T1+jS-1,(ip-1)*nmlat_T1+jS,(ip-1)*nmlat_T1+jS+1, &
            (im-1)*nmlat_T1+jS-1,(im-1)*nmlat_T1+jS,(im-1)*nmlat_T1+jS+1, &
            (i-1)*nmlat_T1+jS-1,(i-1)*nmlat_T1+jS,(i-1)*nmlat_T1+jS+1,(i-1)*nmlat_T1+jN/)
          nzval(rowcnt(ij)+1:rowcnt(ij)+10,ij) = (/ &
            coef(8,jS,i),coef(1,jS,i),coef(2,jS,i), &
            coef(6,jS,i),coef(5,jS,i),coef(4,jS,i), &
            coef(7,jS,i),coef(9,jS,i)-bij(j,i),coef(3,jS,i),bij(j,i)/)
        else
          jcol(rowcnt(ij)+1:rowcnt(ij)+10,ij) = (/ &
            (im-1)*nmlat_T1+jS-1,(im-1)*nmlat_T1+jS,(im-1)*nmlat_T1+jS+1, &
            (i-1)*nmlat_T1+jS-1,(i-1)*nmlat_T1+jS,(i-1)*nmlat_T1+jS+1,(i-1)*nmlat_T1+jN, &
            (ip-1)*nmlat_T1+jS-1,(ip-1)*nmlat_T1+jS,(ip-1)*nmlat_T1+jS+1/)
          nzval(rowcnt(ij)+1:rowcnt(ij)+10,ij) = (/ &
            coef(6,jS,i),coef(5,jS,i),coef(4,jS,i), &
            coef(7,jS,i),coef(9,jS,i)-bij(j,i),coef(3,jS,i),bij(j,i), &
            coef(8,jS,i),coef(1,jS,i),coef(2,jS,i)/)
        endif
        rowcnt(ij) = rowcnt(ij)+10

! note that coef has the direction switched in the NH
        ij = (i-1)*nmlat_T1+jN
!        lhs(ij,(im-1)*nmlat_T1+jN-1) = coef(4,jN,i)
!        lhs(ij,(im-1)*nmlat_T1+jN  ) = coef(5,jN,i)
!        lhs(ij,(im-1)*nmlat_T1+jN+1) = coef(6,jN,i)
!        lhs(ij, (i-1)*nmlat_T1+jS  ) = bij(j,i)              ! conjugate point, b(i,j) Phi*(i,j)
!        lhs(ij, (i-1)*nmlat_T1+jN-1) = coef(3,jN,i)
!        lhs(ij, (i-1)*nmlat_T1+jN  ) = coef(9,jN,i)-bij(j,i)
!        lhs(ij, (i-1)*nmlat_T1+jN+1) = coef(7,jN,i)
!        lhs(ij,(ip-1)*nmlat_T1+jN-1) = coef(2,jN,i)
!        lhs(ij,(ip-1)*nmlat_T1+jN  ) = coef(1,jN,i)
!        lhs(ij,(ip-1)*nmlat_T1+jN+1) = coef(8,jN,i)
        if (i == 1) then
          jcol(rowcnt(ij)+1:rowcnt(ij)+10,ij) = (/ &
            (i-1)*nmlat_T1+jS,(i-1)*nmlat_T1+jN-1,(i-1)*nmlat_T1+jN,(i-1)*nmlat_T1+jN+1, &
            (ip-1)*nmlat_T1+jN-1,(ip-1)*nmlat_T1+jN,(ip-1)*nmlat_T1+jN+1, &
            (im-1)*nmlat_T1+jN-1,(im-1)*nmlat_T1+jN,(im-1)*nmlat_T1+jN+1/)
          nzval(rowcnt(ij)+1:rowcnt(ij)+10,ij) = (/ &
            bij(j,i),coef(3,jN,i),coef(9,jN,i)-bij(j,i),coef(7,jN,i), &
            coef(2,jN,i),coef(1,jN,i),coef(8,jN,i), &
            coef(4,jN,i),coef(5,jN,i),coef(6,jN,i)/)
        elseif (i == nmlon) then
          jcol(rowcnt(ij)+1:rowcnt(ij)+10,ij) = (/ &
            (ip-1)*nmlat_T1+jN-1,(ip-1)*nmlat_T1+jN,(ip-1)*nmlat_T1+jN+1, &
            (im-1)*nmlat_T1+jN-1,(im-1)*nmlat_T1+jN,(im-1)*nmlat_T1+jN+1, &
            (i-1)*nmlat_T1+jS,(i-1)*nmlat_T1+jN-1,(i-1)*nmlat_T1+jN,(i-1)*nmlat_T1+jN+1/)
          nzval(rowcnt(ij)+1:rowcnt(ij)+10,ij) = (/ &
            coef(2,jN,i),coef(1,jN,i),coef(8,jN,i), &
            coef(4,jN,i),coef(5,jN,i),coef(6,jN,i), &
            bij(j,i),coef(3,jN,i),coef(9,jN,i)-bij(j,i),coef(7,jN,i)/)
        else
          jcol(rowcnt(ij)+1:rowcnt(ij)+10,ij) = (/ &
            (im-1)*nmlat_T1+jN-1,(im-1)*nmlat_T1+jN,(im-1)*nmlat_T1+jN+1, &
            (i-1)*nmlat_T1+jS,(i-1)*nmlat_T1+jN-1,(i-1)*nmlat_T1+jN,(i-1)*nmlat_T1+jN+1, &
            (ip-1)*nmlat_T1+jN-1,(ip-1)*nmlat_T1+jN,(ip-1)*nmlat_T1+jN+1/)
          nzval(rowcnt(ij)+1:rowcnt(ij)+10,ij) = (/ &
            coef(4,jN,i),coef(5,jN,i),coef(6,jN,i), &
            bij(j,i),coef(3,jN,i),coef(9,jN,i)-bij(j,i),coef(7,jN,i), &
            coef(2,jN,i),coef(1,jN,i),coef(8,jN,i)/)
        endif
        rowcnt(ij) = rowcnt(ij)+10
      enddo

! at latm_JT, two hemispheres are coupled at j-1
      j = jlatm_JT
      jS = j
      jN = nmlat_T1-j+1

! note that coef has the direction switched in the NH
! therefore the original coef_ns2 c6,c7,c8 at j-1 poleward
! become coef_ns2 c4,c3,c2 at j'+1 poleward
      ij = (i-1)*nmlat_T1+jS
!      lhs(ij,(im-1)*nmlat_T1+jS-1) = coef(6,jS,i)
!      lhs(ij,(im-1)*nmlat_T1+jS  ) = coef(5,jS,i)
!      lhs(ij,(im-1)*nmlat_T1+jS+1) = coef(4,jS,i)
!      lhs(ij,(im-1)*nmlat_T1+jN+1) = coef(6,jN,i)
!      lhs(ij, (i-1)*nmlat_T1+jS-1) = coef(7,jS,i)
!      lhs(ij, (i-1)*nmlat_T1+jS  ) = coef(9,jS,i)
!      lhs(ij, (i-1)*nmlat_T1+jS+1) = coef(3,jS,i)
!      lhs(ij, (i-1)*nmlat_T1+jN+1) = coef(7,jN,i)
!      lhs(ij,(ip-1)*nmlat_T1+jS-1) = coef(8,jS,i)
!      lhs(ij,(ip-1)*nmlat_T1+jS  ) = coef(1,jS,i)
!      lhs(ij,(ip-1)*nmlat_T1+jS+1) = coef(2,jS,i)
!      lhs(ij,(ip-1)*nmlat_T1+jN+1) = coef(8,jN,i)
      if (i == 1) then
        jcol(rowcnt(ij)+1:rowcnt(ij)+12,ij) = (/ &
          (i-1)*nmlat_T1+jS-1,(i-1)*nmlat_T1+jS,(i-1)*nmlat_T1+jS+1,(i-1)*nmlat_T1+jN+1, &
          (ip-1)*nmlat_T1+jS-1,(ip-1)*nmlat_T1+jS,(ip-1)*nmlat_T1+jS+1,(ip-1)*nmlat_T1+jN+1, &
          (im-1)*nmlat_T1+jS-1,(im-1)*nmlat_T1+jS,(im-1)*nmlat_T1+jS+1,(im-1)*nmlat_T1+jN+1/)
        nzval(rowcnt(ij)+1:rowcnt(ij)+12,ij) = (/ &
          coef(7,jS,i),coef(9,jS,i),coef(3,jS,i),coef(7,jN,i), &
          coef(8,jS,i),coef(1,jS,i),coef(2,jS,i),coef(8,jN,i), &
          coef(6,jS,i),coef(5,jS,i),coef(4,jS,i),coef(6,jN,i)/)
      elseif (i == nmlon) then
        jcol(rowcnt(ij)+1:rowcnt(ij)+12,ij) = (/ &
          (ip-1)*nmlat_T1+jS-1,(ip-1)*nmlat_T1+jS,(ip-1)*nmlat_T1+jS+1,(ip-1)*nmlat_T1+jN+1, &
          (im-1)*nmlat_T1+jS-1,(im-1)*nmlat_T1+jS,(im-1)*nmlat_T1+jS+1,(im-1)*nmlat_T1+jN+1, &
          (i-1)*nmlat_T1+jS-1,(i-1)*nmlat_T1+jS,(i-1)*nmlat_T1+jS+1,(i-1)*nmlat_T1+jN+1/)
        nzval(rowcnt(ij)+1:rowcnt(ij)+12,ij) = (/ &
          coef(8,jS,i),coef(1,jS,i),coef(2,jS,i),coef(8,jN,i), &
          coef(6,jS,i),coef(5,jS,i),coef(4,jS,i),coef(6,jN,i), &
          coef(7,jS,i),coef(9,jS,i),coef(3,jS,i),coef(7,jN,i)/)
      else
        jcol(rowcnt(ij)+1:rowcnt(ij)+12,ij) = (/ &
          (im-1)*nmlat_T1+jS-1,(im-1)*nmlat_T1+jS,(im-1)*nmlat_T1+jS+1,(im-1)*nmlat_T1+jN+1, &
          (i-1)*nmlat_T1+jS-1,(i-1)*nmlat_T1+jS,(i-1)*nmlat_T1+jS+1,(i-1)*nmlat_T1+jN+1, &
          (ip-1)*nmlat_T1+jS-1,(ip-1)*nmlat_T1+jS,(ip-1)*nmlat_T1+jS+1,(ip-1)*nmlat_T1+jN+1/)
        nzval(rowcnt(ij)+1:rowcnt(ij)+12,ij) = (/ &
          coef(6,jS,i),coef(5,jS,i),coef(4,jS,i),coef(6,jN,i), &
          coef(7,jS,i),coef(9,jS,i),coef(3,jS,i),coef(7,jN,i), &
          coef(8,jS,i),coef(1,jS,i),coef(2,jS,i),coef(8,jN,i)/)
      endif
      rowcnt(ij) = rowcnt(ij)+12

! note that coef has the direction switched in the NH
! SH coef does not have a direction switch
! therefore coef at j-1 is still coef at j'-1
      ij = (i-1)*nmlat_T1+jN
!      lhs(ij,(im-1)*nmlat_T1+jS-1) = coef(6,jS,i)
!      lhs(ij,(im-1)*nmlat_T1+jN-1) = coef(4,jN,i)
!      lhs(ij,(im-1)*nmlat_T1+jN  ) = coef(5,jN,i)
!      lhs(ij,(im-1)*nmlat_T1+jN+1) = coef(6,jN,i)
!      lhs(ij, (i-1)*nmlat_T1+jS-1) = coef(7,jS,i)
!      lhs(ij, (i-1)*nmlat_T1+jN-1) = coef(3,jN,i)
!      lhs(ij, (i-1)*nmlat_T1+jN  ) = coef(9,jN,i)
!      lhs(ij, (i-1)*nmlat_T1+jN+1) = coef(7,jN,i)
!      lhs(ij,(ip-1)*nmlat_T1+jS-1) = coef(8,jS,i)
!      lhs(ij,(ip-1)*nmlat_T1+jN-1) = coef(2,jN,i)
!      lhs(ij,(ip-1)*nmlat_T1+jN  ) = coef(1,jN,i)
!      lhs(ij,(ip-1)*nmlat_T1+jN+1) = coef(8,jN,i)
      if (i == 1) then
        jcol(rowcnt(ij)+1:rowcnt(ij)+12,ij) = (/ &
          (i-1)*nmlat_T1+jS-1,(i-1)*nmlat_T1+jN-1,(i-1)*nmlat_T1+jN,(i-1)*nmlat_T1+jN+1, &
          (ip-1)*nmlat_T1+jS-1,(ip-1)*nmlat_T1+jN-1,(ip-1)*nmlat_T1+jN,(ip-1)*nmlat_T1+jN+1, &
          (im-1)*nmlat_T1+jS-1,(im-1)*nmlat_T1+jN-1,(im-1)*nmlat_T1+jN,(im-1)*nmlat_T1+jN+1/)
        nzval(rowcnt(ij)+1:rowcnt(ij)+12,ij) = (/ &
          coef(7,jS,i),coef(3,jN,i),coef(9,jN,i),coef(7,jN,i), &
          coef(8,jS,i),coef(2,jN,i),coef(1,jN,i),coef(8,jN,i), &
          coef(6,jS,i),coef(4,jN,i),coef(5,jN,i),coef(6,jN,i)/)
      elseif (i == nmlon) then
        jcol(rowcnt(ij)+1:rowcnt(ij)+12,ij) = (/ &
          (ip-1)*nmlat_T1+jS-1,(ip-1)*nmlat_T1+jN-1,(ip-1)*nmlat_T1+jN,(ip-1)*nmlat_T1+jN+1, &
          (im-1)*nmlat_T1+jS-1,(im-1)*nmlat_T1+jN-1,(im-1)*nmlat_T1+jN,(im-1)*nmlat_T1+jN+1, &
          (i-1)*nmlat_T1+jS-1,(i-1)*nmlat_T1+jN-1,(i-1)*nmlat_T1+jN,(i-1)*nmlat_T1+jN+1/)
        nzval(rowcnt(ij)+1:rowcnt(ij)+12,ij) = (/ &
          coef(8,jS,i),coef(2,jN,i),coef(1,jN,i),coef(8,jN,i), &
          coef(6,jS,i),coef(4,jN,i),coef(5,jN,i),coef(6,jN,i), &
          coef(7,jS,i),coef(3,jN,i),coef(9,jN,i),coef(7,jN,i)/)
      else
        jcol(rowcnt(ij)+1:rowcnt(ij)+12,ij) = (/ &
          (im-1)*nmlat_T1+jS-1,(im-1)*nmlat_T1+jN-1,(im-1)*nmlat_T1+jN,(im-1)*nmlat_T1+jN+1, &
          (i-1)*nmlat_T1+jS-1,(i-1)*nmlat_T1+jN-1,(i-1)*nmlat_T1+jN,(i-1)*nmlat_T1+jN+1, &
          (ip-1)*nmlat_T1+jS-1,(ip-1)*nmlat_T1+jN-1,(ip-1)*nmlat_T1+jN,(ip-1)*nmlat_T1+jN+1/)
        nzval(rowcnt(ij)+1:rowcnt(ij)+12,ij) = (/ &
          coef(6,jS,i),coef(4,jN,i),coef(5,jN,i),coef(6,jN,i), &
          coef(7,jS,i),coef(3,jN,i),coef(9,jN,i),coef(7,jN,i), &
          coef(8,jS,i),coef(2,jN,i),coef(1,jN,i),coef(8,jN,i)/)
      endif
      rowcnt(ij) = rowcnt(ij)+12

! from latm_JT to equator, symmetric solution
      do j = jlatm_JT+1,nmlat_h-1
        jS = j
        jN = nmlat_T1-j+1

        ij = (i-1)*nmlat_T1+jS
!        lhs(ij,(im-1)*nmlat_T1+jS-1) = coef(6,jS,i)
!        lhs(ij,(im-1)*nmlat_T1+jS  ) = coef(5,jS,i)
!        lhs(ij,(im-1)*nmlat_T1+jS+1) = coef(4,jS,i)
!        lhs(ij, (i-1)*nmlat_T1+jS-1) = coef(7,jS,i)
!        lhs(ij, (i-1)*nmlat_T1+jS  ) = coef(9,jS,i)
!        lhs(ij, (i-1)*nmlat_T1+jS+1) = coef(3,jS,i)
!        lhs(ij,(ip-1)*nmlat_T1+jS-1) = coef(8,jS,i)
!        lhs(ij,(ip-1)*nmlat_T1+jS  ) = coef(1,jS,i)
!        lhs(ij,(ip-1)*nmlat_T1+jS+1) = coef(2,jS,i)
        if (i == 1) then
          jcol(rowcnt(ij)+1:rowcnt(ij)+9,ij) = (/ &
            (i-1)*nmlat_T1+jS-1,(i-1)*nmlat_T1+jS,(i-1)*nmlat_T1+jS+1, &
            (ip-1)*nmlat_T1+jS-1,(ip-1)*nmlat_T1+jS,(ip-1)*nmlat_T1+jS+1, &
            (im-1)*nmlat_T1+jS-1,(im-1)*nmlat_T1+jS,(im-1)*nmlat_T1+jS+1/)
          nzval(rowcnt(ij)+1:rowcnt(ij)+9,ij) = (/ &
            coef(7,jS,i),coef(9,jS,i),coef(3,jS,i), &
            coef(8,jS,i),coef(1,jS,i),coef(2,jS,i), &
            coef(6,jS,i),coef(5,jS,i),coef(4,jS,i)/)
        elseif (i == nmlon) then
          jcol(rowcnt(ij)+1:rowcnt(ij)+9,ij) = (/ &
            (ip-1)*nmlat_T1+jS-1,(ip-1)*nmlat_T1+jS,(ip-1)*nmlat_T1+jS+1, &
            (im-1)*nmlat_T1+jS-1,(im-1)*nmlat_T1+jS,(im-1)*nmlat_T1+jS+1, &
            (i-1)*nmlat_T1+jS-1,(i-1)*nmlat_T1+jS,(i-1)*nmlat_T1+jS+1/)
          nzval(rowcnt(ij)+1:rowcnt(ij)+9,ij) = (/ &
            coef(8,jS,i),coef(1,jS,i),coef(2,jS,i), &
            coef(6,jS,i),coef(5,jS,i),coef(4,jS,i), &
            coef(7,jS,i),coef(9,jS,i),coef(3,jS,i)/)
        else
          jcol(rowcnt(ij)+1:rowcnt(ij)+9,ij) = (/ &
            (im-1)*nmlat_T1+jS-1,(im-1)*nmlat_T1+jS,(im-1)*nmlat_T1+jS+1, &
            (i-1)*nmlat_T1+jS-1,(i-1)*nmlat_T1+jS,(i-1)*nmlat_T1+jS+1, &
            (ip-1)*nmlat_T1+jS-1,(ip-1)*nmlat_T1+jS,(ip-1)*nmlat_T1+jS+1/)
          nzval(rowcnt(ij)+1:rowcnt(ij)+9,ij) = (/ &
            coef(6,jS,i),coef(5,jS,i),coef(4,jS,i), &
            coef(7,jS,i),coef(9,jS,i),coef(3,jS,i), &
            coef(8,jS,i),coef(1,jS,i),coef(2,jS,i)/)
        endif
        rowcnt(ij) = rowcnt(ij)+9

! note that coef has the direction switched in the NH
        ij = (i-1)*nmlat_T1+jN
!        lhs(ij,(im-1)*nmlat_T1+jN-1) = coef(4,jN,i)
!        lhs(ij,(im-1)*nmlat_T1+jN  ) = coef(5,jN,i)
!        lhs(ij,(im-1)*nmlat_T1+jN+1) = coef(6,jN,i)
!        lhs(ij, (i-1)*nmlat_T1+jN-1) = coef(3,jN,i)
!        lhs(ij, (i-1)*nmlat_T1+jN  ) = coef(9,jN,i)
!        lhs(ij, (i-1)*nmlat_T1+jN+1) = coef(7,jN,i)
!        lhs(ij,(ip-1)*nmlat_T1+jN-1) = coef(2,jN,i)
!        lhs(ij,(ip-1)*nmlat_T1+jN  ) = coef(1,jN,i)
!        lhs(ij,(ip-1)*nmlat_T1+jN+1) = coef(8,jN,i)
        if (i == 1) then
          jcol(rowcnt(ij)+1:rowcnt(ij)+9,ij) = (/ &
            (i-1)*nmlat_T1+jN-1,(i-1)*nmlat_T1+jN,(i-1)*nmlat_T1+jN+1, &
            (ip-1)*nmlat_T1+jN-1,(ip-1)*nmlat_T1+jN,(ip-1)*nmlat_T1+jN+1, &
            (im-1)*nmlat_T1+jN-1,(im-1)*nmlat_T1+jN,(im-1)*nmlat_T1+jN+1/)
          nzval(rowcnt(ij)+1:rowcnt(ij)+9,ij) = (/ &
            coef(3,jN,i),coef(9,jN,i),coef(7,jN,i), &
            coef(2,jN,i),coef(1,jN,i),coef(8,jN,i), &
            coef(4,jN,i),coef(5,jN,i),coef(6,jN,i)/)
        elseif (i == nmlon) then
          jcol(rowcnt(ij)+1:rowcnt(ij)+9,ij) = (/ &
            (ip-1)*nmlat_T1+jN-1,(ip-1)*nmlat_T1+jN,(ip-1)*nmlat_T1+jN+1, &
            (im-1)*nmlat_T1+jN-1,(im-1)*nmlat_T1+jN,(im-1)*nmlat_T1+jN+1, &
            (i-1)*nmlat_T1+jN-1,(i-1)*nmlat_T1+jN,(i-1)*nmlat_T1+jN+1/)
          nzval(rowcnt(ij)+1:rowcnt(ij)+9,ij) = (/ &
            coef(2,jN,i),coef(1,jN,i),coef(8,jN,i), &
            coef(4,jN,i),coef(5,jN,i),coef(6,jN,i), &
            coef(3,jN,i),coef(9,jN,i),coef(7,jN,i)/)
        else
          jcol(rowcnt(ij)+1:rowcnt(ij)+9,ij) = (/ &
            (im-1)*nmlat_T1+jN-1,(im-1)*nmlat_T1+jN,(im-1)*nmlat_T1+jN+1, &
            (i-1)*nmlat_T1+jN-1,(i-1)*nmlat_T1+jN,(i-1)*nmlat_T1+jN+1, &
            (ip-1)*nmlat_T1+jN-1,(ip-1)*nmlat_T1+jN,(ip-1)*nmlat_T1+jN+1/)
          nzval(rowcnt(ij)+1:rowcnt(ij)+9,ij) = (/ &
            coef(4,jN,i),coef(5,jN,i),coef(6,jN,i), &
            coef(3,jN,i),coef(9,jN,i),coef(7,jN,i), &
            coef(2,jN,i),coef(1,jN,i),coef(8,jN,i)/)
        endif
        rowcnt(ij) = rowcnt(ij)+9
      enddo

! equator
      j = nmlat_h

      ij = (i-1)*nmlat_T1+j
!      lhs(ij,(im-1)*nmlat_T1+j-1) = coef(6,j,i)
!      lhs(ij,(im-1)*nmlat_T1+j  ) = coef(5,j,i)
!      lhs(ij,(im-1)*nmlat_T1+j+1) = coef(4,j,i)
!      lhs(ij, (i-1)*nmlat_T1+j-1) = coef(7,j,i)
!      lhs(ij, (i-1)*nmlat_T1+j  ) = coef(9,j,i)
!      lhs(ij, (i-1)*nmlat_T1+j+1) = coef(3,j,i)
!      lhs(ij,(ip-1)*nmlat_T1+j-1) = coef(8,j,i)
!      lhs(ij,(ip-1)*nmlat_T1+j  ) = coef(1,j,i)
!      lhs(ij,(ip-1)*nmlat_T1+j+1) = coef(2,j,i)
      if (i == 1) then
        jcol(rowcnt(ij)+1:rowcnt(ij)+9,ij) = (/ &
          (i-1)*nmlat_T1+j-1,(i-1)*nmlat_T1+j,(i-1)*nmlat_T1+j+1, &
          (ip-1)*nmlat_T1+j-1,(ip-1)*nmlat_T1+j,(ip-1)*nmlat_T1+j+1, &
          (im-1)*nmlat_T1+j-1,(im-1)*nmlat_T1+j,(im-1)*nmlat_T1+j+1/)
        nzval(rowcnt(ij)+1:rowcnt(ij)+9,ij) = (/ &
          coef(7,j,i),coef(9,j,i),coef(3,j,i), &
          coef(8,j,i),coef(1,j,i),coef(2,j,i), &
          coef(6,j,i),coef(5,j,i),coef(4,j,i)/)
      elseif (i == nmlon) then
        jcol(rowcnt(ij)+1:rowcnt(ij)+9,ij) = (/ &
          (ip-1)*nmlat_T1+j-1,(ip-1)*nmlat_T1+j,(ip-1)*nmlat_T1+j+1, &
          (im-1)*nmlat_T1+j-1,(im-1)*nmlat_T1+j,(im-1)*nmlat_T1+j+1, &
          (i-1)*nmlat_T1+j-1,(i-1)*nmlat_T1+j,(i-1)*nmlat_T1+j+1/)
        nzval(rowcnt(ij)+1:rowcnt(ij)+9,ij) = (/ &
          coef(8,j,i),coef(1,j,i),coef(2,j,i), &
          coef(6,j,i),coef(5,j,i),coef(4,j,i), &
          coef(7,j,i),coef(9,j,i),coef(3,j,i)/)
      else
        jcol(rowcnt(ij)+1:rowcnt(ij)+9,ij) = (/ &
          (im-1)*nmlat_T1+j-1,(im-1)*nmlat_T1+j,(im-1)*nmlat_T1+j+1, &
          (i-1)*nmlat_T1+j-1,(i-1)*nmlat_T1+j,(i-1)*nmlat_T1+j+1, &
          (ip-1)*nmlat_T1+j-1,(ip-1)*nmlat_T1+j,(ip-1)*nmlat_T1+j+1/)
        nzval(rowcnt(ij)+1:rowcnt(ij)+9,ij) = (/ &
          coef(6,j,i),coef(5,j,i),coef(4,j,i), &
          coef(7,j,i),coef(9,j,i),coef(3,j,i), &
          coef(8,j,i),coef(1,j,i),coef(2,j,i)/)
      endif
      rowcnt(ij) = rowcnt(ij)+9
    enddo

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
endmodule dist_solver_module
