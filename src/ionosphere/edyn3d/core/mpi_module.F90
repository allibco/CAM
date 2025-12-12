#define PARALLEL
module mpi_module

  use prec, only: rp
  use cam_logfile, only: iulog
  use cam_abortutils, only: endrun
#ifdef PARALLEL
  use MPI
  use iso_fortran_env, only: real32,real64
#endif

  implicit none

!nmlon is also in params module, which can cause confusion
  
  integer :: dynamo_world=-huge(1), extra_world=-huge(1),&
       union_world=-huge(1),&
       mpi_rp=-huge(1), mpi_size=0, mpi_rank=-1, &
       ex_mpi_rank=-1,ex_mpi_size =0,& 
       un_mpi_rank=-1, un_mpi_size=0,&
       lat_size=0, lon_size=0, lat_rank=-1, lon_rank=-1, &
       nmlat=0, maxmlat=-1, mlat0=1, mlat1=0, mlatd0=1, mlatd1=0, &
       nmlon=0, maxmlon=-1, mlon0=1, mlon1=0, mlond0=1, mlond1=0, &
       ij_start_s=1, ij_stop_s=0, ij_start_n=1, ij_stop_n=0, &
       mpi_partner=-1, my_recvgrid_size=0,  &
       mygrid_size =0, mygrid_size_s=0, mygrid_size_n=0, &
       mpi_comm_host_rank=-1, my_sendgrid_size=0

  integer, dimension(:), allocatable :: &
       nmlat_task, mlat0_task, mlat1_task, &
       nmlon_task, mlon0_task, mlon1_task, &
       task_lat_offset, task_csr_rowstarts

#ifdef PARALLEL
   integer :: host_group, dynamo_group, extra_group, union_group
#endif
  
  interface gather_mag ! gather magnetic fields
    module procedure gather_mag_2d, gather_mag_3d, gather_mag_4d, gather_mag_5d
  endinterface

  contains
!-----------------------------------------------------------------------
  subroutine init(mpi_comm_host, npes_edyn3d)

    integer, intent(in) :: mpi_comm_host
    integer, intent(in) :: npes_edyn3d

#ifdef PARALLEL
    integer :: ierror, mid
    integer :: npes_host
    integer :: ranges(3, 2)
    
    if (rp == real32) then
      mpi_rp = MPI_REAL4
    elseif (rp == real64) then
      mpi_rp = MPI_REAL8
    else
      stop 'unknown real precision'
    endif

    ! host comm rank and size
    call mpi_comm_size(mpi_comm_host, npes_host, ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Comm_size', ierror)

    call MPI_Comm_rank(mpi_comm_host, mpi_comm_host_rank, ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Comm_rank', ierror)
    
    !create dynamo communicator and addition one for the solve
    ! MUST be even (or 1) and at most half the number of procs in mpi_comm_host (unless host is size 1)
    mid = npes_host / 2
    if (npes_host == 1 .and. npes_edyn3d == 1) then
       write(iulog,*) 'MPI init check: size for edyn3d and cam are both size 1'
    else
       if (npes_edyn3d > mid) then
          write(iulog,*) 'MPI init ERROR: npes_edyn3d can be at most half the size of the number of procs used for cam, so must be <=', mid
          call endrun('MPI init ERROR: npes_edyn3d can be at most half the size of the number of procs used for cam')
       elseif ((npes_edyn3d /= 1 .and. mod(npes_edyn3d, 2) /= 0)) then
          write(iulog,*) 'MPI init ERROR: npes_edyn3d must be even'
          call endrun('MPI init ERROR: npes_edyn3d must be even')
       endif
    endif
    
    
    if (npes_edyn3d < npes_host) then ! create two subgroups
       !one for the majority of the dynamo work and a second for the linear solve
       !double check
       if (npes_edyn3d*2 > npes_host) then
          write(iulog,*) 'MPI init ERROR: npes_edyn3d size'
          call endrun('MPI init ERROR: npes_edyn3d size')
       end if
    
       !first get original group
       call MPI_Comm_group(mpi_comm_host, host_group, ierror)
       if (ierror /= MPI_SUCCESS) call handle_error('MPI_Comm_group', ierror)

       ! Define the rank range for the first subgroup (first, last, stride)
       ranges(1, 1) = 0               ! First rank
       ranges(2, 1) = npes_edyn3d - 1 ! Last rank                                                   
       ranges(3, 1) = 1               ! Stride of 1

       ! Define the rank range for the second subgroup (first, last, stride)
       ranges(1, 2) = npes_edyn3d       ! First rank
       ranges(2, 2) = 2*npes_edyn3d - 1 ! Last rank                                                   
       ranges(3, 2) = 1                 ! Stride of 1
  
       ! Step 3: Create the new groups using the rank range                                      
       call MPI_Group_range_incl(host_group, 1, ranges(:,1:1), dynamo_group, ierror)
       if (ierror /= MPI_SUCCESS) call handle_error('MPI_Group_range_incl', ierror)

       call MPI_Group_range_incl(host_group, 1, ranges(:,2:2), extra_group, ierror)
       if (ierror /= MPI_SUCCESS) call handle_error('MPI_Group_range_incl group2', ierror)
       
       ! Step 4: Create new communicators (dynamo_world, extra_world) from the new groups
        call MPI_Comm_create(mpi_comm_host, dynamo_group, dynamo_world, ierror)
       if (ierror /= MPI_SUCCESS) call handle_error('MPI_Comm_create dynamo_world', ierror)

       call MPI_Comm_create(mpi_comm_host, extra_group, extra_world, ierror)
       if (ierror /= MPI_SUCCESS) call handle_error('MPI_Comm_create extra_world', ierror)

       !Now we want the union group which compines dynamo_world and extra_world
       call MPI_Group_union(dynamo_group, extra_group, union_group, ierror)
       call MPI_Comm_create(mpi_comm_host, union_group, union_world, ierror)
       if (ierror /= MPI_SUCCESS) call handle_error('MPI_Comm_create union_world', ierror)
       
       ! Check if the process is part of the new communicators                
       if (dynamo_world /= MPI_COMM_NULL) then
          call MPI_Comm_rank(dynamo_world, mpi_rank, ierror)
          call MPI_Comm_size(dynamo_world, mpi_size, ierror)
       else
          mpi_size = npes_edyn3d !let non-participating procs know the sizes (for                                                    !other functions, all need to know lat_size and lon_size)    
          mpi_rank = -1 !not participating                                                    
       endif

       if (extra_world /= MPI_COMM_NULL) then
          call MPI_Comm_rank(extra_world, ex_mpi_rank, ierror)
          call MPI_Comm_size(extra_world, ex_mpi_size, ierror)
       else
          ex_mpi_size = npes_edyn3d
          ex_mpi_rank = -1
       endif

       if (union_world /= MPI_COMM_NULL) then
          call MPI_Comm_rank(union_world, un_mpi_rank, ierror)
          call MPI_Comm_size(union_world, un_mpi_size, ierror)
       else
          un_mpi_size = 2*npes_edyn3d
          un_mpi_rank = -1
       endif
       
    else   !just duplicate      (both # of cam procs and number of dyno are size 1)                                                                     
       call MPI_Comm_dup(mpi_comm_host, dynamo_world, ierror)
       if (ierror /= MPI_SUCCESS) call handle_error('MPI_Comm_dup', ierror)

       call MPI_Comm_dup(mpi_comm_host, extra_world, ierror)
       if (ierror /= MPI_SUCCESS) call handle_error('MPI_Comm_dup', ierror)

       call MPI_Comm_dup(mpi_comm_host, union_world, ierror)
       if (ierror /= MPI_SUCCESS) call handle_error('MPI_Comm_dup', ierror)
       
       ! Get the rank in the new communicators                                                 
       call MPI_Comm_rank(dynamo_world, mpi_rank, ierror)
       call MPI_Comm_size(dynamo_world, mpi_size, ierror)

       call MPI_Comm_rank(extra_world, ex_mpi_rank, ierror)
       call MPI_Comm_size(extra_world, ex_mpi_size, ierror)

       call MPI_Comm_rank(union_world, un_mpi_rank, ierror)
       call MPI_Comm_size(union_world, un_mpi_size, ierror)
       
    endif
    
#else
    mpi_rp = rp
    mpi_size = 1
    un_mpi_size = 1
    ex_mpi_size = 1
    ex_mpi_rank = 0
    un_mpi_rank = 0
    mpi_rank = 0
#endif

    !figure out the process grid for the dynamo procs
    if (mpi_size == 1) then
       lat_size = 1
       lon_size = 1
    elseif (mpi_size == 2) then
       lat_size = 1
       lon_size = 2
    else
       ! Original factorization loop for other cases
       do lat_size = int(sqrt(real(mpi_size, kind=rp))), 1, -1
          lon_size = mpi_size / lat_size
          if (lon_size*lat_size == mpi_size) exit ! lon_size >= lat_size
       enddo
    endif

! stack along latitudes first then longitudes
! (lat_size=3)
!  8  9 10 11
!  4  5  6  7
!  0  1  2  3 (lon_size=4)
    if (mpi_rank >= 0) then ! main group
       lat_rank = mpi_rank / lon_size
       lon_rank = modulo(mpi_rank, lon_size)
    elseif (ex_mpi_rank >= 0) then !extra group (use this later to find partner)
       lat_rank = ex_mpi_rank / lon_size
       lon_rank = modulo(ex_mpi_rank, lon_size)
       !swap bc north hemi is mirror
       lat_rank = lat_size - lat_rank - 1
    endif

  endsubroutine init
!-----------------------------------------------------------------------
  subroutine setup_topology(nmlat_in, nmlon_in)
    ! setup MPI decompositions in geo and mag coordinates and the connectivity matrix
    ! note: non-active procs have mpi_rank < 0
    
    use params_module, only:nmlat_h,nmlat_T1
    
    integer, intent(in) :: nmlat_in, nmlon_in

    integer :: i, j, rnk, rnki, rnkj, mlat0_n, mlat1_n, &
               mysize_n, mysize_s, cnt, out_int
    integer, dimension(:), allocatable :: task_mygrid_size
    
    allocate(nmlat_task(0:lat_size-1))
    allocate(nmlon_task(0:lon_size-1))

    nmlat_task = 0
    nmlon_task = 0

    allocate(mlat0_task(0:mpi_size-1))
    allocate(mlat1_task(0:mpi_size-1))
    allocate(mlon0_task(0:mpi_size-1))
    allocate(mlon1_task(0:mpi_size-1))
    mlat0_task = 1
    mlon0_task = 1
    mlat1_task = -1
    mlon1_task = -1
 
    !for dist
    allocate(task_csr_rowstarts(0:un_mpi_size))
    allocate(task_mygrid_size(0:un_mpi_size-1))
    allocate(task_lat_offset(0:lat_size-1))
    task_csr_rowstarts = 0
    task_mygrid_size = 0
    task_lat_offset = 0

    !set global vars in this module
    nmlat = nmlat_in
    nmlon = nmlon_in

! setup magnetic grid decomposition
! each process can have unequal number of latitudes or longitudes

    nmlat_task = generate_minvar_list(nmlat, lat_size)
    maxmlat = maxval(nmlat_task)
    if (mpi_rank >= 0) then
       mlat0 = 1
       do j = 0, lat_rank-1
          mlat0 = mlat0 + nmlat_task(j)
       enddo
       mlat1 = mlat0 + nmlat_task(lat_rank) - 1
    endif

    nmlon_task = generate_minvar_list(nmlon, lon_size)
    maxmlon = maxval(nmlon_task)
    if (mpi_rank >= 0) then
       mlon0 = 1
       do i = 0, lon_rank-1
          mlon0 = mlon0 + nmlon_task(i)
       enddo
       mlon1 = mlon0 + nmlon_task(lon_rank) - 1
    endif

    ! each process keeps a record of the lat-lon decomposition
    do concurrent (rnk = 0:mpi_size-1)
      !position on proc grid 
      rnkj = rnk / lon_size
      rnki = modulo(rnk, lon_size)

      mlat0_task(rnk) = 1
      do j = 0, rnkj-1
        mlat0_task(rnk) = mlat0_task(rnk) + nmlat_task(j)
      enddo
      mlat1_task(rnk) = mlat0_task(rnk) + nmlat_task(rnkj) - 1

      mlon0_task(rnk) = 1
      do i = 0, rnki-1
        mlon0_task(rnk) = mlon0_task(rnk) + nmlon_task(i)
      enddo
      mlon1_task(rnk) = mlon0_task(rnk) + nmlon_task(rnki) - 1
    enddo
        
    !calc offset for lat to faciliate matrix creation
    do i = 1, lat_size-1
       task_lat_offset(i) = task_lat_offset(i-1) + nmlat_task(i-1)
    enddo


    !Find partner for north hemispher data for the LU
    !! here is south hemisphere:
    ! stack along latitudes first then longitudes
    ! (lat_size=3)
    !  8  9 10 11
    !  4  5  6  7
    !  0  1  2  3 (lon_size=4)

    ! Then north is mirror:
    !  0  1  2  3 (lon_size=4)
    !  4  5  6  7
    !  8  9 10 11

    ! and we map north to
    ! 20 21 22 23
    ! 16 17 18 19
    ! 12 13 14 15
    ! Why? LU solver needs contiguous rows on increasing proc #s
    ! recall in extra group for north, proc numbering is also
    !  8  9 10 11  (lat rank 0)
    !  4  5  6  7
    !  0  1  2  3  (lat rank 2)

    ! number of procs is even (partners are with the global numbering in larger group)
    if (un_mpi_size > 1) then
       if (mpi_rank >=0 ) then !active in main group
          mpi_partner = 2*mpi_size - lon_size*(lat_rank +1) + lon_rank 
       elseif (ex_mpi_rank >= 0) then ! active extra group
          mpi_partner = lon_size*(lat_rank) + lon_rank 
       else !not active in either group
          mpi_partner = -1
       endif
    else !one proc
       mpi_partner = 0
    endif

    !initial matrix row start and stops
    if (un_mpi_rank >=0) then !active for union
       if (mpi_rank >=0 ) then !active for dyno_group
          !s hemi
          ij_start_s = calc_grid_ij(mlon0,mlat0,lat_rank)
          ij_stop_s = calc_grid_ij(mlon1,mlat1,lat_rank)
          !n hemi
          !corresponding j for n hemisphere
          mlat0_n = nmlat_T1 - mlat0 + 1
          mlat1_n =  nmlat_T1 - mlat1 + 1
          if (mlat1_n == nmlat_h) then !equator (don't double count)
             mlat1_n  = mlat1_n + 1
          endif
          !now mlat0_n will be bigger than mlat1_n in north hemisphere
          ij_start_n = calc_grid_ij(mlon0,mlat1_n,lat_rank)
          ij_stop_n = calc_grid_ij(mlon1,mlat0_n,lat_rank)
       
          !my sizes in each hemisphere
          mysize_n =  (ij_stop_n -ij_start_n + 1)
          mysize_s = (ij_stop_s -ij_start_s + 1)
          !write(*,*) 'AB: SETUP TOPO mysize_s, mysize_n', mysize_s, mysize_n
          !write(*, *) 'AB: TOPO ij_start_s, ij_stop_s, s_grid_pts = ', ij_start_s, ij_stop_s, mysize_s
          !write(*, *) 'AB: TOPO ij_start_n, ij_stop_n, n_grid_pts = ', ij_start_n, ij_stop_n, mysize_n
       
          !set global vars (n & s row counts will be diff for procs on equator)
          mygrid_size_n = mysize_n
          mygrid_size_s = mysize_s
       !elseif (ex_mpi_rank >= 0) then
       !   write(*,*) 'SETUP TOPO ex_mpi_rank', ex_mpi_rank
       endif 
       
       !dynamo group will send north hemi to extra group
       !so north needs to know the size
       if (un_mpi_size > 1) then
           out_int = mpi_partner_size(mysize_n)
           if (mpi_rank >= 0) then !south
              my_sendgrid_size = mysize_n
              my_recvgrid_size = 0
              mygrid_size = mysize_s
           elseif (ex_mpi_rank >= 0) then !north
              my_sendgrid_size = 0
              my_recvgrid_size = out_int
              mygrid_size = out_int
           endif
           
          !now we need to calculate the rowstarts for the global block
          !csr martix - this will be 0-based indeing for superlu
          !do an allgather to get each procs grid size
          task_mygrid_size = all_gather_int(mygrid_size, union_world, un_mpi_size)
          !write(*,*) 'TASK grid size = ', task_mygrid_size
          do i=0, un_mpi_size-1
             task_csr_rowstarts(i+1) = task_csr_rowstarts(i) &
                  + task_mygrid_size(i)
          enddo
       else !one proc
          mygrid_size = mysize_n + mysize_s
          my_sendgrid_size = 0
          my_recvgrid_size = 0
          task_csr_rowstarts(1) = mygrid_size
       endif

       ! set halo global variables
       mlatd0 = mlat0 - 1
       mlatd1 = mlat1 + 1
       mlond0 = mlon0 - 1
       mlond1 = mlon1 + 1

       !write(*, *) 'AB: TOPO un_mpi_rank, mpi_rank, ex_mpi_rank, mlat0, mlat1, mlon0,mlon1',un_mpi_rank, mpi_rank, ex_mpi_rank, mlat0, mlat1, mlon0,mlon1
       
    else !completely non-active
       mygrid_size = 0
       my_sendgrid_size = 0
       my_recvgrid_size = 0
    endif

    deallocate(task_mygrid_size)

    
  endsubroutine setup_topology
!-----------------------------------------------------------------------
              
  subroutine dynamo_mpi_finalize

#ifdef PARALLEL
    use MPI

    integer :: ierror

    ! Clean up                                                                                                    
    if (mpi_rank >=0) then
       call MPI_Comm_free(dynamo_world, ierror)
    endif
    if (ex_mpi_rank >=0) then
       call MPI_Comm_free(extra_world, ierror)
    endif
    if (un_mpi_rank >=0) then
       call MPI_Comm_free(union_world, ierror)
    endif
    
    call MPI_Group_free(host_group, ierror)
    Call MPI_Group_free(dynamo_group, ierror)
    call MPI_Group_free(extra_group, ierror)
    call MPI_Group_free(union_group, ierror)

#endif

    deallocate(mlat0_task)
    deallocate(mlon0_task)
    deallocate(mlat1_task)
    deallocate(mlon1_task)

    deallocate(nmlat_task)
    deallocate(nmlon_task)

    deallocate(task_csr_rowstarts)
    deallocate(task_lat_offset)
    

    
  endsubroutine dynamo_mpi_finalize

 !-----------------------------------------------------------------------

!-----------------------------------------------------------------------

function mpi_partner_size(intin) result(intout)

#ifdef PARALLEL
    use MPI
#endif

    integer,  intent(in) :: intin
    integer :: intout
   
#ifdef PARALLEL
    integer :: ierr
    integer :: tag = 88
    integer :: send_request, recv_request

    intout = 0

    if (mpi_rank >= 0) then
       !post send
       call MPI_Isend(intin, 1, MPI_INTEGER, mpi_partner, tag, &
            union_world, send_request, ierr)
       if (ierr /= MPI_SUCCESS) call handle_error('MPI_Isend', ierr)

       ! Wait for send to complete
       call MPI_Wait(send_request, MPI_STATUS_IGNORE, ierr)
       if (ierr /= MPI_SUCCESS) call handle_error('MPI_Wait', ierr)

       intout = 0

    elseif (ex_mpi_rank >=0) then 
       !post receive
       call MPI_Irecv(intout, 1, MPI_INTEGER, mpi_partner, tag, &
            union_world, recv_request, ierr)
       if (ierr /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierr)
    
       ! Wait for receive to complete
       call MPI_Wait(recv_request, MPI_STATUS_IGNORE, ierr)
       if (ierr /= MPI_SUCCESS) call handle_error('MPI_Wait', ierr)

    endif
#else

    intout = intin

#endif
    
endfunction mpi_partner_size

!-----------------------------------------------------------------------

subroutine mpi_sendnorth_mat(nnz_per_row, my_rowptr, my_values, my_cols, &
     partner_rowptr, partner_values, partner_cols)

  !send my csr matrix info  or recv my partner's matrix info
  
  !note: the south pole's denser row does not get sent (only north data)  

#ifdef PARALLEL
    use MPI
#endif
    integer, intent(in) :: nnz_per_row
    
    integer, dimension(my_sendgrid_size+1), intent(in) :: my_rowptr
    integer, dimension(my_sendgrid_size*nnz_per_row), intent(in) :: my_cols
    real(kind=rp), dimension(my_sendgrid_size*nnz_per_row), intent(in) :: my_values

    integer, dimension(my_recvgrid_size+1), intent(out) :: partner_rowptr
    integer, dimension(my_recvgrid_size*nnz_per_row), intent(out) :: partner_cols
    real(kind=rp), dimension(my_recvgrid_size*nnz_per_row), intent(out) :: partner_values

    integer :: i, nnz_count
    
#ifdef PARALLEL

    integer ierr
    integer, dimension(3) :: send_request, recv_request

   
    !send north to my partner
    if (mpi_rank >= 0) then 
       call MPI_Isend(my_rowptr, my_sendgrid_size + 1, MPI_INTEGER, mpi_partner, 200, &
            union_world, send_request(1), ierr)
       if (ierr /= MPI_SUCCESS) call handle_error('MPI_Isend', ierr)

       call MPI_Isend(my_cols, my_sendgrid_size*nnz_per_row, MPI_INTEGER, mpi_partner, 201, union_world, &
            send_request(2), ierr)
       if (ierr /= MPI_SUCCESS) call handle_error('MPI_Isend', ierr)

       call MPI_Isend(my_values, my_sendgrid_size*nnz_per_row, mpi_rp, mpi_partner, 202, union_world, &
            send_request(3), ierr)
       if (ierr /= MPI_SUCCESS) call handle_error('MPI_Isend', ierr)

         ! Wait for all sends to complete
       call MPI_Waitall(3, send_request, MPI_STATUSES_IGNORE, ierr)
       if (ierr /= MPI_SUCCESS) call handle_error('MPI_Waitall', ierr)


    elseif (ex_mpi_rank >=0) then
       
       partner_rowptr = 0
       partner_cols = 0
       partner_values = 0.0
       
       !recv from my partner
       call MPI_Irecv(partner_rowptr, my_recvgrid_size + 1, MPI_INTEGER, mpi_partner, &
            200, union_world, recv_request(1), ierr)
       if (ierr /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierr)
       
       call MPI_Irecv(partner_cols, my_recvgrid_size*nnz_per_row, MPI_INTEGER, mpi_partner, &
            201, union_world, recv_request(2), ierr)
       if (ierr /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierr)
       
       call MPI_Irecv(partner_values, my_recvgrid_size*nnz_per_row, mpi_rp, mpi_partner, &
            202, union_world, recv_request(3), ierr)
       if (ierr /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierr)
      
       !Wait for all recvs
       call MPI_Waitall(3, recv_request, MPI_STATUSES_IGNORE, ierr)
       if (ierr /= MPI_SUCCESS) call handle_error('MPI_Waitall', ierr)

    endif
    
#else
!serial
    do concurrent (i = 1: mygrid_size + 1)
       partner_rowptr(i) = my_rowptr(i)
    enddo
    !note for nnz calc that we are 1-based?
    nnz_count = my_rowptr(mygrid_size + 1) -1
    do concurrent (i = 1: nnz_count
       partner_cols(i) = my_cols(i)
       partner_values(i) = my_values(i)
    enddo
    
#endif
    
    
  endsubroutine mpi_sendnorth_mat
  
!-----------------------------------------------------------------------
subroutine mpi_sendnorth_vec(send_values, recv_values)
!send my_values to partner or recv from partner_values  

#ifdef PARALLEL
    use MPI
#endif

    real(kind=rp), dimension(my_sendgrid_size), intent(in) :: send_values
    real(kind=rp), dimension(my_recvgrid_size), intent(out) :: recv_values
#ifdef PARALLEL

    integer :: ierr
    integer :: send_request, recv_request

    !send north to my partner
    if (mpi_rank >= 0) then !send
   
       call MPI_Isend(send_values, my_sendgrid_size, mpi_rp, mpi_partner, 400, union_world, &
            send_request, ierr)
       if (ierr /= MPI_SUCCESS) call handle_error('MPI_Isend', ierr)

       ! Wait for my send to complete
       call MPI_Wait(send_request, MPI_STATUSES_IGNORE, ierr)
       if (ierr /= MPI_SUCCESS) call handle_error('MPI_Wait', ierr)

    elseif (ex_mpi_rank >= 0) then !recv north from my partner
       
       !recv from my partner
       call MPI_Irecv(recv_values, my_recvgrid_size, mpi_rp, mpi_partner, &
            400, union_world, recv_request, ierr)
       if (ierr /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierr)
       
    
       !Wait for my recv
       call MPI_Wait(recv_request, MPI_STATUSES_IGNORE, ierr)
       if (ierr /= MPI_SUCCESS) call handle_error('MPI_Wait', ierr)

    endif

#else
    !serial
    !do nothing
    
#endif
    
endsubroutine mpi_sendnorth_vec

!-----------------------------------------------------------------------
subroutine mpi_recvnorth_vec(send_values, recv_values)
!send my_values to partner or recv from partner_values  

#ifdef PARALLEL
    use MPI
#endif
    !note these are reversed for the recvnorth
    real(kind=rp), dimension(my_recvgrid_size), intent(in) :: send_values
    real(kind=rp), dimension(my_sendgrid_size), intent(out) :: recv_values
#ifdef PARALLEL

    integer :: ierr
    integer :: send_request, recv_request

    !send north back to my partner
    if (ex_mpi_rank >= 0) then !send
   
       call MPI_Isend(send_values, my_recvgrid_size, mpi_rp, mpi_partner, 400, union_world, &
            send_request, ierr)
       if (ierr /= MPI_SUCCESS) call handle_error('MPI_Isend', ierr)

       ! Wait for my send to complete
       call MPI_Wait(send_request, MPI_STATUSES_IGNORE, ierr)
       if (ierr /= MPI_SUCCESS) call handle_error('MPI_Wait', ierr)

    elseif (mpi_rank >= 0) then !recv north from partner
       
       !recv from my partner
       call MPI_Irecv(recv_values, my_sendgrid_size, mpi_rp, mpi_partner, &
            400, union_world, recv_request, ierr)
       if (ierr /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierr)
    
       !Wait for my recv
       call MPI_Wait(recv_request, MPI_STATUSES_IGNORE, ierr)
       if (ierr /= MPI_SUCCESS) call handle_error('MPI_Wait', ierr)

    endif

#else
    !serial
    !do nothing
    
#endif
    
endsubroutine mpi_recvnorth_vec
!-----------------------------------------------------------------------

  subroutine sync_mlat_5d(var, l, m, n)
! longitude halo points are not included

#ifdef PARALLEL
    use MPI
#endif

    integer, intent(in) :: l, m, n
    real(kind=rp), dimension(l, m, n, mlatd0:mlatd1, mlon0:mlon1), intent(inout) :: var

#ifdef PARALLEL
    integer :: below, above, cnt, i, lc, mc, nc, ierr
    integer, dimension(4) :: request
    real(kind=rp), dimension(l, m, n, maxmlon) :: &
      send_to_below, send_to_above, recv_from_below, recv_from_above

! find the rank of adjacent processes
    if (lat_rank == 0) then
      below = MPI_PROC_NULL
    else
      below = mpi_rank - lon_size
    endif

    if (lat_rank == lat_size-1) then
      above = MPI_PROC_NULL
    else
      above = mpi_rank + lon_size
    endif

    cnt = l * m * n * maxmlon

! load to work array
    do concurrent (i = 1:mlon1-mlon0+1, nc = 1:n, mc = 1:m, lc = 1:l)
      send_to_below(lc, mc, nc, i) = var(lc, mc, nc, mlat0, i+mlon0-1)
      send_to_above(lc, mc, nc, i) = var(lc, mc, nc, mlat1, i+mlon0-1)
    enddo

! sync in latitude
    call MPI_Isend(send_to_below, cnt, mpi_rp, &
      below, 0, dynamo_world, request(1), ierr)
    if (ierr /= MPI_SUCCESS) call handle_error('MPI_Isend', ierr)

    call MPI_Isend(send_to_above, cnt, mpi_rp, &
      above, 1, dynamo_world, request(2), ierr)
    if (ierr /= MPI_SUCCESS) call handle_error('MPI_Isend', ierr)

    call MPI_Irecv(recv_from_above, cnt, mpi_rp, &
      above, 0, dynamo_world, request(3), ierr)
    if (ierr /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierr)

    call MPI_Irecv(recv_from_below, cnt, mpi_rp, &
      below, 1, dynamo_world, request(4), ierr)
    if (ierr /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierr)

! wait for sync to complete
    call MPI_Waitall(4, request, MPI_STATUSES_IGNORE, ierr)
    if (ierr /= MPI_SUCCESS) call handle_error('MPI_Waitall', ierr)

! unpack to model fields
    if (lat_rank /= 0) then
      do concurrent (i = mlon0:mlon1, nc = 1:n, mc = 1:m, lc = 1:l)
        var(lc, mc, nc, mlatd0, i) = recv_from_below(lc, mc, nc, i-mlon0+1)
      enddo
    endif

    if (lat_rank /= lat_size-1) then
      do concurrent (i = mlon0:mlon1, nc = 1:n, mc = 1:m, lc = 1:l)
        var(lc, mc, nc, mlatd1, i) = recv_from_above(lc, mc, nc, i-mlon0+1)
      enddo
    endif
#endif

  endsubroutine sync_mlat_5d
!-----------------------------------------------------------------------
  subroutine sync_mlon_5d(var, l, m, n)

#ifdef PARALLEL
    use MPI
#endif

    integer, intent(in) :: l, m, n
    real(kind=rp), dimension(l, m, n, mlatd0:mlatd1, mlond0:mlond1), intent(inout) :: var

    integer :: j, lc, mc, nc

#ifdef PARALLEL
    integer :: left, right, cnt, ierr
    integer, dimension(4) :: request
    real(kind=rp), dimension(l, m, n, maxmlat+4) :: &
      send_to_left, send_to_right, recv_from_left, recv_from_right

! find the rank of adjacent processes
    if (lon_rank == 0) then
      left = mpi_rank - 1 + lon_size
    else
      left = mpi_rank - 1
    endif

    if (lon_rank == lon_size-1) then
      right = mpi_rank + 1 - lon_size
    else
      right = mpi_rank + 1
    endif

    cnt = l * m * n * (maxmlat + 4)

! load to work array
    do concurrent (j = 1:mlatd1-mlatd0+1, nc = 1:n, mc = 1:m, lc = 1:l)
      send_to_left(lc, mc, nc, j) = var(lc, mc, nc, j+mlatd0-1, mlon0)
      send_to_right(lc, mc, nc, j) = var(lc, mc, nc, j+mlatd0-1, mlon1)
    enddo

! sync in longitude
    call MPI_Isend(send_to_left, cnt, mpi_rp, &
      left, 0, dynamo_world, request(1), ierr)
    if (ierr /= MPI_SUCCESS) call handle_error('MPI_Isend', ierr)

    call MPI_Isend(send_to_right, cnt, mpi_rp, &
      right, 1, dynamo_world, request(2), ierr)
    if (ierr /= MPI_SUCCESS) call handle_error('MPI_Isend', ierr)

    call MPI_Irecv(recv_from_right, cnt, mpi_rp, &
      right, 0, dynamo_world, request(3), ierr)
    if (ierr /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierr)

    call MPI_Irecv(recv_from_left, cnt, mpi_rp, &
      left, 1, dynamo_world, request(4), ierr)
    if (ierr /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierr)

! wait for sync to complete
    call MPI_Waitall(4, request, MPI_STATUSES_IGNORE, ierr)
    if (ierr /= MPI_SUCCESS) call handle_error('MPI_Waitall', ierr)

! unpack to model fields
    do concurrent (j = mlatd0:mlatd1, nc = 1:n, mc = 1:m, lc = 1:l)
      var(lc, mc, nc, j, mlond0) = recv_from_left(lc, mc, nc, j-mlatd0+1)
      var(lc, mc, nc, j, mlond1) = recv_from_right(lc, mc, nc, j-mlatd0+1)
    enddo
#else
    do concurrent (j = mlatd0:mlatd1, nc = 1:n, mc = 1:m, lc = 1:l)
      var(lc, mc, nc, j, mlond0) = var(lc, mc, nc, j, mlon1)
      var(lc, mc, nc, j, mlond1) = var(lc, mc, nc, j, mlon0)
    enddo
#endif

  endsubroutine sync_mlon_5d
!-----------------------------------------------------------------------
  subroutine sync_mlat_3d(var, n)
! longitude halo points are not included

#ifdef PARALLEL
    use MPI
#endif

    integer, intent(in) ::  n
    real(kind=rp), dimension(n, mlatd0:mlatd1, mlon0:mlon1), intent(inout) :: var

#ifdef PARALLEL
    integer :: below, above, cnt, i, nc, ierr
    integer, dimension(4) :: request
    real(kind=rp), dimension(n, maxmlon) :: &
      send_to_below, send_to_above, recv_from_below, recv_from_above

! find the rank of adjacent processes
    if (lat_rank == 0) then
      below = MPI_PROC_NULL
    else
      below = mpi_rank - lon_size
    endif

    if (lat_rank == lat_size-1) then
      above = MPI_PROC_NULL
    else
      above = mpi_rank + lon_size
    endif

    cnt = n * maxmlon

! load to work array
    do concurrent (i = 1:mlon1-mlon0+1, nc = 1:n)
      send_to_below(nc, i) = var(nc, mlat0, i+mlon0-1)
      send_to_above(nc, i) = var(nc, mlat1, i+mlon0-1)
    enddo

! sync in latitude
    call MPI_Isend(send_to_below, cnt, mpi_rp, &
      below, 0, dynamo_world, request(1), ierr)
    if (ierr /= MPI_SUCCESS) call handle_error('MPI_Isend', ierr)

    call MPI_Isend(send_to_above, cnt, mpi_rp, &
      above, 1, dynamo_world, request(2), ierr)
    if (ierr /= MPI_SUCCESS) call handle_error('MPI_Isend', ierr)

    call MPI_Irecv(recv_from_above, cnt, mpi_rp, &
      above, 0, dynamo_world, request(3), ierr)
    if (ierr /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierr)

    call MPI_Irecv(recv_from_below, cnt, mpi_rp, &
      below, 1, dynamo_world, request(4), ierr)
    if (ierr /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierr)

! wait for sync to complete
    call MPI_Waitall(4, request, MPI_STATUSES_IGNORE, ierr)
    if (ierr /= MPI_SUCCESS) call handle_error('MPI_Waitall', ierr)

! unpack to model fields
    if (lat_rank /= 0) then
      do concurrent (i = mlon0:mlon1, nc = 1:n)
        var(nc, mlatd0, i) = recv_from_below(nc, i-mlon0+1)
      enddo
    endif

    if (lat_rank /= lat_size-1) then
      do concurrent (i = mlon0:mlon1, nc = 1:n)
        var(nc, mlatd1, i) = recv_from_above(nc, i-mlon0+1)
      enddo
    endif
#endif

  endsubroutine sync_mlat_3d
  !-----------------------------------------------------------------------

  subroutine sync_mlon_3d(var, n)

#ifdef PARALLEL
    use MPI
#endif

    integer, intent(in) ::  n
    real(kind=rp), dimension( n, mlatd0:mlatd1, mlond0:mlond1), intent(inout) :: var

    integer :: j, nc

#ifdef PARALLEL
    integer :: left, right, cnt, ierr
    integer, dimension(4) :: request
    real(kind=rp), dimension( n, maxmlat+4) :: &
      send_to_left, send_to_right, recv_from_left, recv_from_right

! find the rank of adjacent processes
    if (lon_rank == 0) then
      left = mpi_rank - 1 + lon_size
    else
      left = mpi_rank - 1
    endif

    if (lon_rank == lon_size-1) then
      right = mpi_rank + 1 - lon_size
    else
      right = mpi_rank + 1
    endif

    cnt = n * (maxmlat + 4)

! load to work array
    do concurrent (j = 1:mlatd1-mlatd0+1, nc = 1:n)
      send_to_left(nc, j) = var(nc, j+mlatd0-1, mlon0)
      send_to_right(nc, j) = var(nc, j+mlatd0-1, mlon1)
    enddo

! sync in longitude
    call MPI_Isend(send_to_left, cnt, mpi_rp, &
      left, 0, dynamo_world, request(1), ierr)
    if (ierr /= MPI_SUCCESS) call handle_error('MPI_Isend', ierr)

    call MPI_Isend(send_to_right, cnt, mpi_rp, &
      right, 1, dynamo_world, request(2), ierr)
    if (ierr /= MPI_SUCCESS) call handle_error('MPI_Isend', ierr)

    call MPI_Irecv(recv_from_right, cnt, mpi_rp, &
      right, 0, dynamo_world, request(3), ierr)
    if (ierr /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierr)

    call MPI_Irecv(recv_from_left, cnt, mpi_rp, &
      left, 1, dynamo_world, request(4), ierr)
    if (ierr /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierr)

! wait for sync to complete
    call MPI_Waitall(4, request, MPI_STATUSES_IGNORE, ierr)
    if (ierr /= MPI_SUCCESS) call handle_error('MPI_Waitall', ierr)

! unpack to model fields
    do concurrent (j = mlatd0:mlatd1, nc = 1:n)
      var(nc, j, mlond0) = recv_from_left(nc, j-mlatd0+1)
      var(nc, j, mlond1) = recv_from_right(nc, j-mlatd0+1)
    enddo
#else
    do concurrent (j = mlatd0:mlatd1, nc = 1:n)
      var(nc, j, mlond0) = var(nc, j, mlon1)
      var(nc, j, mlond1) = var(nc, j, mlon0)
    enddo
#endif

  endsubroutine sync_mlon_3d
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

function all_gather_int(intin, comm, commsize) result(intarrayout)

#ifdef PARALLEL
    use MPI
#endif

    integer, intent(in) :: intin, comm, commsize

    integer, dimension(0:commsize-1) :: intarrayout


#ifdef PARALLEL
    integer :: ierr, cnt

    cnt = 1
    call MPI_Allgather(intin, cnt, MPI_INTEGER, &
        intarrayout, cnt, MPI_INTEGER, comm, ierr)
    if (ierr /= MPI_SUCCESS) call handle_error('MPI_Allgather', ierr)

#else

    intarrayout(0) = intin

#endif
    
endfunction all_gather_int

    
  
!-----------------------------------------------------------------------
  function gather_lon_1d(varin) result(varout)
    ! collect a 1d array from other procs w/lat_rank 0 to root proc 0
    ! this could be genearlized to have any root proc and any proc row (Or column)
    !DYNAMO WORLD
#ifdef PARALLEL

    use MPI
#endif

    real(kind=rp), dimension(mlon0:mlon1), intent(in) :: varin
    real(kind=rp), dimension(nmlon) :: varout

    integer :: i,j

#ifdef PARALLEL
    integer :: ierr, cnt, myrequest, tag, rs, re
    integer, dimension(1:lon_size-1) :: requests  
    real(kind=rp), dimension(maxmlon) :: sendbuf
    real(kind=rp), dimension(nmlon) :: recvbuf


    ! Add this debug output at the start
!    write(*,*) 'AB0: DEBUG - nmlon =', nmlon
!    write(*,*) 'AB: DEBUG - mlon0, mlon1 =', mlon0, mlon1
!    write(*,*) 'AB: DEBUG - size(varin) =', size(varin)
!#ifdef PARALLEL
!    write(*,*) 'AB: DEBUG - lon_size =', lon_size
!    write(*,*) 'AB: DEBUG - mpi_rank =', mpi_rank
!#endif
    
    tag = 34
    
    ! load to work array
    sendbuf=0.0
    recvbuf=0.0
    varout=0.0
   
    ! gather data to 0 from proc_row 0
    if (mpi_rank == 0) then !root proc receives data

       !copy my data
       do concurrent (i = mlon0:mlon1)
          varout(i) = varin(i)
       enddo
       if (lon_size > 1) then !get data from other procs 1:lon_size -1
          do i = 1,lon_size -1
             rs = mlon0_task(i)
             re = mlon1_task(i)
             cnt = re-rs+1
             
             !write(*,*) 'AB: Gather, mpi_rank, i, rs, re', mpi_rank, i, rs, re
             
             call MPI_Irecv(recvbuf(rs:re), cnt, mpi_rp, &
                  i, tag, dynamo_world, &
                  requests(i), ierr)
             if (ierr /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierr)
          enddo

          !now wait to receive all data
          call MPI_Waitall(lon_size-1, requests, MPI_STATUSES_IGNORE, ierr)
          if (ierr /= MPI_SUCCESS) call handle_error('MPI_Waitall', ierr)
          
          !write(*,*) 'AB: size(varout), size(recvbuf)', size(varout), size(recvbuf)

          
          !copy from all other procs
          do i = 1, lon_size - 1
             rs = mlon0_task(i)
             re = mlon1_task(i)
             do concurrent (j = rs:re)
                varout(j) = recvbuf(j)
             end do
          end do
       endif
       
    elseif (mpi_rank > 0 .and. lat_rank == 0) then ! send info to root (0)

       cnt = mlon1-mlon0+1
       sendbuf(1:cnt) = varin(mlon0:mlon1)
       
       call MPI_Isend(sendbuf, cnt, mpi_rp, &
            0, tag, dynamo_world, &
            myrequest, ierr)
       if (ierr /= MPI_SUCCESS) call handle_error('MPI_Isend', ierr)


       call MPI_WAIT(myrequest, MPI_STATUS_IGNORE, ierr)
       if (ierr /= MPI_SUCCESS) call handle_error('MPI_Wait', ierr)
       
    endif

   
#else
    do concurrent (i = mlon0:mlon1)
       varout(i) = varin(i)
    enddo
#endif
    
  endfunction gather_lon_1d

!-----------------------------------------------------------------------
  function gather_mlon_3d(varin, m, n) result(varout)

#ifdef PARALLEL
    use MPI
#endif

    integer, intent(in) :: m, n
    real(kind=rp), dimension(m, n, mlon0:mlon1), intent(in) :: varin
    real(kind=rp), dimension(m, n, nmlon) :: varout

    integer :: i, mc, nc

#ifdef PARALLEL
    integer :: cnt, rnki, i0, i1, ierr
    integer, dimension(0:lon_size*2-1) :: request
    real(kind=rp), dimension(m, n, maxmlon) :: sendbuf
    real(kind=rp), dimension(m, n, maxmlon, 0:lon_size-1) :: recvbuf

    cnt = m * n * maxmlon

! load to work array
    do concurrent (i = 1:mlon1-mlon0+1, nc = 1:n, mc = 1:m)
      sendbuf(mc, nc, i) = varin(mc, nc, i+mlon0-1)
    enddo

! gather longitudes
    do rnki = 0, lon_size-1
      call MPI_Isend(sendbuf, cnt, mpi_rp, &
        lat_rank*lon_size + rnki, 0, dynamo_world, &
        request(rnki), ierr)
      if (ierr /= MPI_SUCCESS) call handle_error('MPI_Isend', ierr)

      call MPI_Irecv(recvbuf(:, :, :, rnki), cnt, mpi_rp, &
        lat_rank*lon_size + rnki, 0, dynamo_world, &
        request(lon_size+rnki), ierr)
      if (ierr /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierr)
    enddo

! wait for gather to complete
    call MPI_Waitall(lon_size*2, request, MPI_STATUSES_IGNORE, ierr)
    if (ierr /= MPI_SUCCESS) call handle_error('MPI_Waitall', ierr)

! reconstruct longitude rings
    do concurrent (rnki = 0:lon_size-1)
      i0 = mlon0_task(lat_rank*lon_size + rnki)
      i1 = mlon1_task(lat_rank*lon_size + rnki)
      do concurrent (i = i0:i1, nc = 1:n, mc = 1:m)
        varout(mc, nc, i) = recvbuf(mc, nc, i-i0+1, rnki)
      enddo
    enddo
#else
    do concurrent (i = mlon0:mlon1, nc = 1:n, mc = 1:m)
      varout(mc, nc, i) = varin(mc, nc, i)
    enddo
#endif

  endfunction gather_mlon_3d
!-----------------------------------------------------------------------
  function gather_mag_2d(varin, root) result(varout)

#ifdef PARALLEL
    use MPI
#endif

    integer, intent(in) :: root
    real(kind=rp), dimension(mlat0:mlat1, mlon0:mlon1), intent(in) :: varin
    real(kind=rp), dimension(nmlat, nmlon) :: varout

    integer :: i, j

#ifdef PARALLEL
    integer :: cnt, rnk, i0, i1, j0, j1, ierr
    real(kind=rp), dimension(maxmlat, maxmlon) :: sendbuf
    real(kind=rp), dimension(maxmlat, maxmlon, 0:mpi_size-1) :: recvbuf

    cnt = maxmlat * maxmlon

! load to work array
    do concurrent (i = 1:mlon1-mlon0+1, j = 1:mlat1-mlat0+1)
      sendbuf(j, i) = varin(j+mlat0-1, i+mlon0-1)
    enddo

    if (root < 0) then
      call MPI_Allgather(sendbuf, cnt, mpi_rp, &
        recvbuf, cnt, mpi_rp, dynamo_world, ierr)
      if (ierr /= MPI_SUCCESS) call handle_error('MPI_Allgather', ierr)
    else
      call MPI_Gather(sendbuf, cnt, mpi_rp, &
        recvbuf, cnt, mpi_rp, root, dynamo_world, ierr)
      if (ierr /= MPI_SUCCESS) call handle_error('MPI_Gather', ierr)
    endif

! reconstruct based on lat-lon decomposition
    if (root<0 .or. mpi_rank==root) then
      do concurrent (rnk = 0:mpi_size-1)
        i0 = mlon0_task(rnk)
        i1 = mlon1_task(rnk)
        j0 = mlat0_task(rnk)
        j1 = mlat1_task(rnk)
        do concurrent (i = i0:i1, j = j0:j1)
          varout(j, i) = recvbuf(j-j0+1, i-i0+1, rnk)
        enddo
      enddo
    endif
#else
    do concurrent (i = mlon0:mlon1, j = mlat0:mlat1)
      varout(j, i) = varin(j, i)
    enddo
#endif

  endfunction gather_mag_2d
!-----------------------------------------------------------------------
  function gather_mag_3d(varin, n, root) result(varout)

#ifdef PARALLEL
    use MPI
#endif

    integer, intent(in) :: n, root
    real(kind=rp), dimension(n, mlat0:mlat1, mlon0:mlon1), intent(in) :: varin
    real(kind=rp), dimension(n, nmlat, nmlon) :: varout

    integer :: i, j, nc

#ifdef PARALLEL
    integer :: cnt, rnk, i0, i1, j0, j1, ierr
    real(kind=rp), dimension(n, maxmlat, maxmlon) :: sendbuf
    real(kind=rp), dimension(n, maxmlat, maxmlon, 0:mpi_size-1) :: recvbuf

    cnt = n * maxmlat * maxmlon

! load to work array
    do concurrent (i = 1:mlon1-mlon0+1, j = 1:mlat1-mlat0+1, nc = 1:n)
      sendbuf(nc, j, i) = varin(nc, j+mlat0-1, i+mlon0-1)
    enddo

    if (root < 0) then
      call MPI_Allgather(sendbuf, cnt, mpi_rp, &
        recvbuf, cnt, mpi_rp, dynamo_world, ierr)
      if (ierr /= MPI_SUCCESS) call handle_error('MPI_Allgather', ierr)
    else
      call MPI_Gather(sendbuf, cnt, mpi_rp, &
        recvbuf, cnt, mpi_rp, root, dynamo_world, ierr)
      if (ierr /= MPI_SUCCESS) call handle_error('MPI_Gather', ierr)
    endif

! reconstruct based on lat-lon decomposition
    if (root<0 .or. mpi_rank==root) then
      do concurrent (rnk = 0:mpi_size-1)
        i0 = mlon0_task(rnk)
        i1 = mlon1_task(rnk)
        j0 = mlat0_task(rnk)
        j1 = mlat1_task(rnk)
        do concurrent (i = i0:i1, j = j0:j1, nc = 1:n)
          varout(nc, j, i) = recvbuf(nc, j-j0+1, i-i0+1, rnk)
        enddo
      enddo
    endif
#else
    do concurrent (i = mlon0:mlon1, j = mlat0:mlat1, nc = 1:n)
      varout(nc, j, i) = varin(nc, j, i)
    enddo
#endif

  endfunction gather_mag_3d
!-----------------------------------------------------------------------
  function gather_mag_4d(varin, m, n, root) result(varout)

#ifdef PARALLEL
    use MPI
#endif

    integer, intent(in) :: m, n, root
    real(kind=rp), dimension(m, n, mlat0:mlat1, mlon0:mlon1), intent(in) :: varin
    real(kind=rp), dimension(m, n, nmlat, nmlon) :: varout

    integer :: i, j, mc, nc

#ifdef PARALLEL
    integer :: cnt, rnk, i0, i1, j0, j1, ierr
    real(kind=rp), dimension(m, n, maxmlat, maxmlon) :: sendbuf
    real(kind=rp), dimension(m, n, maxmlat, maxmlon, 0:mpi_size-1) :: recvbuf

    cnt = m * n * maxmlat * maxmlon

! load to work array
    do concurrent (i = 1:mlon1-mlon0+1, j = 1:mlat1-mlat0+1, nc = 1:n, mc = 1:m)
      sendbuf(mc, nc, j, i) = varin(mc, nc, j+mlat0-1, i+mlon0-1)
    enddo

    if (root < 0) then
      call MPI_Allgather(sendbuf, cnt, mpi_rp, &
        recvbuf, cnt, mpi_rp, dynamo_world, ierr)
      if (ierr /= MPI_SUCCESS) call handle_error('MPI_Allgather', ierr)
    else
      call MPI_Gather(sendbuf, cnt, mpi_rp, &
        recvbuf, cnt, mpi_rp, root, dynamo_world, ierr)
      if (ierr /= MPI_SUCCESS) call handle_error('MPI_Gather', ierr)
    endif

! reconstruct based on lat-lon decomposition
    if (root<0 .or. mpi_rank==root) then
      do concurrent (rnk = 0:mpi_size-1)
        i0 = mlon0_task(rnk)
        i1 = mlon1_task(rnk)
        j0 = mlat0_task(rnk)
        j1 = mlat1_task(rnk)
        do concurrent (i = i0:i1, j = j0:j1, nc = 1:n, mc = 1:m)
          varout(mc, nc, j, i) = recvbuf(mc, nc, j-j0+1, i-i0+1, rnk)
        enddo
      enddo
    endif
#else
    do concurrent (i = mlon0:mlon1, j = mlat0:mlat1, nc = 1:n, mc = 1:m)
      varout(mc, nc, j, i) = varin(mc, nc, j, i)
    enddo
#endif

  endfunction gather_mag_4d
!-----------------------------------------------------------------------
  function gather_mag_5d(varin, l, m, n, root) result(varout)

#ifdef PARALLEL
    use MPI
#endif

    integer, intent(in) :: l, m, n, root
    real(kind=rp), dimension(l, m, n, mlat0:mlat1, mlon0:mlon1), intent(in) :: varin
    real(kind=rp), dimension(l, m, n, nmlat, nmlon) :: varout

    integer :: i, j, lc, mc, nc

#ifdef PARALLEL
    integer :: cnt, rnk, i0, i1, j0, j1, ierr
    real(kind=rp), dimension(l, m, n, maxmlat, maxmlon) :: sendbuf
    real(kind=rp), dimension(l, m, n, maxmlat, maxmlon, 0:mpi_size-1) :: recvbuf

    cnt = l * m * n * maxmlat * maxmlon

! load to work array
    do concurrent (i = 1:mlon1-mlon0+1, j = 1:mlat1-mlat0+1, nc = 1:n, mc = 1:m, lc = 1:l)
      sendbuf(lc, mc, nc, j, i) = varin(lc, mc, nc, j+mlat0-1, i+mlon0-1)
    enddo

    if (root < 0) then
      call MPI_Allgather(sendbuf, cnt, mpi_rp, &
        recvbuf, cnt, mpi_rp, dynamo_world, ierr)
      if (ierr /= MPI_SUCCESS) call handle_error('MPI_Allgather', ierr)
    else
      call MPI_Gather(sendbuf, cnt, mpi_rp, &
        recvbuf, cnt, mpi_rp, root, dynamo_world, ierr)
      if (ierr /= MPI_SUCCESS) call handle_error('MPI_Gather', ierr)
    endif

! reconstruct based on lat-lon decomposition
    if (root<0 .or. mpi_rank==root) then
      do concurrent (rnk = 0:mpi_size-1)
        i0 = mlon0_task(rnk)
        i1 = mlon1_task(rnk)
        j0 = mlat0_task(rnk)
        j1 = mlat1_task(rnk)
        do concurrent (i = i0:i1, j = j0:j1, nc = 1:n, mc = 1:m, lc = 1:l)
          varout(lc, mc, nc, j, i) = recvbuf(lc, mc, nc, j-j0+1, i-i0+1, rnk)
        enddo
      enddo
    endif
#else
    do concurrent (i = mlon0:mlon1, j = mlat0:mlat1, nc = 1:n, mc = 1:m, lc = 1:l)
      varout(lc, mc, nc, j, i) = varin(lc, mc, nc, j, i)
    enddo
#endif

  endfunction gather_mag_5d
!-----------------------------------------------------------------------
  subroutine bcast_3d(var, l, m, n, root)

#ifdef PARALLEL
    use MPI
#endif

    integer, intent(in) :: l, m, n, root
    real(kind=rp), dimension(l, m, n), intent(inout) :: var

#ifdef PARALLEL
    integer :: cnt, lc, mc, nc, ierr
    real(kind=rp), dimension(l, m, n) :: buffer

    cnt = l * m * n

! load to work array
    do concurrent (nc = 1:n, mc = 1:m, lc = 1:l)
      buffer(lc, mc, nc) = var(lc, mc, nc)
    enddo

    call MPI_Bcast(buffer, cnt, mpi_rp, root, dynamo_world, ierr)
    if (ierr /= MPI_SUCCESS) call handle_error('MPI_Bcast', ierr)

! unpack to model fields
    do concurrent (nc = 1:n, mc = 1:m, lc = 1:l)
      var(lc, mc, nc) = buffer(lc, mc, nc)
    enddo
#endif

  endsubroutine bcast_3d
!-----------------------------------------------------------------------
  function reduce_sum_1d(varin, n, root) result(varout)

#ifdef PARALLEL
    use MPI
#endif

    integer, intent(in) :: n, root
    real(kind=rp), dimension(n), intent(in) :: varin
    real(kind=rp), dimension(n) :: varout

#ifdef PARALLEL
    integer :: ierr

    if (root < 0) then
      call MPI_Allreduce(varin, varout, n, mpi_rp, &
        MPI_SUM, dynamo_world, ierr)
      if (ierr /= MPI_SUCCESS) call handle_error('MPI_Allreduce', ierr)
    else
      call MPI_Reduce(varin, varout, n, mpi_rp, &
        MPI_SUM, root, dynamo_world, ierr)
      if (ierr /= MPI_SUCCESS) call handle_error('MPI_Reduce', ierr)
    endif
#else
    integer :: nc

    do concurrent (nc = 1:n)
      varout(nc) = varin(nc)
    enddo
#endif

  endfunction reduce_sum_1d
!-----------------------------------------------------------------------
  subroutine handle_error(funcname, errorcode)

#ifdef PARALLEL
    use MPI
#endif

    character(len=*), intent(in) :: funcname
    integer, intent(in) :: errorcode

#ifdef PARALLEL
    character(len=MPI_MAX_ERROR_STRING) :: string
    integer :: resultlen, ierr

    call MPI_Error_string(errorcode, string, resultlen, ierr)
    write(6, "('MPI error encountered: ', a, ', when calling ', a, '. Finalizing...')") &
      trim(string), trim(funcname)
    call MPI_Finalize(ierr)
#endif

  endsubroutine handle_error
!-----------------------------------------------------------------------
  pure function generate_minvar_list(summation, nfactor) result(factor_list)
! generate a list of the given length summing up to the given number
! which has the minimum variances (all factors differ at most by 1)
! e.g., 5 out of 26 would be 5,5,5,5,6

    integer, intent(in) :: summation, nfactor
    integer, dimension(nfactor) :: factor_list

    integer :: minfactor, maxfactor, mincnt
    real(kind=rp) :: factor

    factor = real(summation, kind=rp) / nfactor
    minfactor = floor(factor)
    maxfactor = ceiling(factor)

    mincnt = maxfactor*nfactor - summation

    factor_list(1:mincnt) = minfactor
    factor_list(mincnt+1:nfactor) = maxfactor

  endfunction generate_minvar_list
  !-----------------------------------------------------------------------
   !-----------------------------------------------------------------------

  ! find the coef matrix row number for the grid location  i,j
  ! j_in (latitude) can be in north or south hemisphere
  ! my_latrank is the position of the *calling* processor
  ! we do not calculate my_latrank  from j, because this is how we determine
  ! if j lives on the calling processor
  pure function calc_grid_ij(i_in,j_in,my_latrank) result(ij)
    
    use params_module,only:nmlat_T1, nmlat_h

    integer, intent(in) :: i_in, j_in, my_latrank
    integer:: ij
    
    integer:: m, my_numlat, jS, latrank, i, j
    
    if (mpi_size == 1) then
       ! this is not the same as the original ordering
       !south and north hemispheres will be a numbered partition
       !old: ij = (i_in-1)*nmlat_T1+j_in
       if (j_in <= nmlat_h) then !south hemi or equator
          ij = (i_in-1)*nmlat_h+j_in
       else
          ij = nmlat_h*nmlon + (i_in-1)*(nmlat_h-1) + (j_in-nmlat_h)
       endif
    endif

    if (j_in <= nmlat_h) then !south hemi or equator

       !calc latrank
       !if j_in is in ghost layer for the calling proc, then adjust latrank
       if (j_in == mlat1 + 1) then
          latrank  = my_latrank + 1
       elseif (j_in == mlat0 -1) then
          latrank  = my_latrank - 1
       else
          latrank = my_latrank
       endif
       
       if (latrank < 0 .OR. latrank >= lat_size) then
          !write(iulog,*) "Error with latrank in calc_grid_ij: mpi_rank, latrank", mpi_rank, latrank
          ij = -1
          return
       endif
       
       !num of latitude points in proc parition    
       my_numlat = nmlat_task(latrank)

       !calc i
       !adjust if i is on edge of global domain
       if (i_in == 0) then
          i = nmlon
       elseif (i_in == nmlon+1) then
          i = 1
       else
          i = i_in
       endif

       !no mods to j for south
       j = j_in

       !offset
       m = task_lat_offset(latrank)
       
    else !north hemi (j_in is in north hemisphere)

       !calc latrank
       !need to know if j_in is in the ghost point for the
       !calling proc (and adjust latrank )
       jS = nmlat_T1 -j_in +1
       if (jS == mlat1 + 1) then
          latrank  = my_latrank + 1
       elseif (jS == mlat0 -1) then
          latrank  = my_latrank - 1
       else
          latrank = my_latrank
       endif

        if (latrank < 0 .OR. latrank >= lat_size) then
           !write(iulog,*) "Error with latrank in calc_grid_ij: mpi_rank, latrank", mpi_rank, latrank
           ij = -2
           return
       endif
              
       !num of latitude points in proc parition    
       my_numlat = nmlat_task(latrank)

       !calc i
       !adjust if i is on edge of global domain
       if (i_in == 0) then
          
          i = nmlon
       elseif (i_in == nmlon+1) then
          i = 1
       else
          i = i_in
       endif
       
       !if i am in the task row that includes equator, then
       ! subtract 1 from num_lat (becuz don't count equator
       ! twice)
       if (latrank == lat_size - 1 ) then
          my_numlat = my_numlat - 1
       endif

       !no change to j
       j = j_in
       
       !m is different in the in the north hemisphere
       m = nmlat_T1 - task_lat_offset(latrank) - my_numlat
    endif

    !calculate ij
    ij = (nmlon*m) + (j-m) + ((i-1)*my_numlat)

  end function calc_grid_ij

  ! -----------------------------------------------------------------------


endmodule mpi_module
