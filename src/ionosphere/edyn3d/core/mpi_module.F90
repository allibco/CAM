#define PARALLEL
module mpi_module

  use prec, only: rp
  use cam_logfile, only: iulog
#ifdef PARALLEL
  use MPI
  use iso_fortran_env, only: real32,real64
#endif

  implicit none

!nmlon is also in params module, which can cause confusion
  
  integer :: dynamo_world=-huge(1), &
    mpi_rp=-huge(1), mpi_size=0, mpi_rank=-1, &
    lat_size=0, lon_size=0, lat_rank=-1, lon_rank=-1, &
    nmlat=0, maxmlat=-1, mlat0=1, mlat1=0, mlatd0=1, mlatd1=0, &
    nmlon=0, maxmlon=-1, mlon0=1, mlon1=0, mlond0=1, mlond1=0, &
    ij_start_s=1, ij_stop_s=0, ij_start_n=1, ij_stop_n=0, &
    mpi_partner=-1, partner_hgridsize=0, my_hgridsize=0, &
    mygrid_size =0, mygrid_size_s=0, mygrid_size_n=0, &
    mpi_comm_host_rank=-1, my_sendgrid_size=0

  integer, dimension(:), allocatable :: &
    nmlat_task, mlat0_task, mlat1_task, &
    nmlon_task, mlon0_task, mlon1_task, &
    task_lat_offset, task_csr_rowstarts, &
    task_csr_mapping

#ifdef PARALLEL
   integer :: host_group, dynamo_group
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
    integer :: ierror
    integer :: npes_host
    integer :: ranges(3, 1)
    
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
    
    !create dynamo communicator (could be a subset of mpi_comm_host or the same size)                                                                       
    if (npes_edyn3d < npes_host) then ! create a subgroup                                                                                                                                  
       !first get original group
       call MPI_Comm_group(mpi_comm_host, host_group, ierror)
       if (ierror /= MPI_SUCCESS) call handle_error('MPI_Comm_group', ierror)

       ! Define the rank range for the subgroup (first, last, stride)
       ranges(1, 1) = 0               ! First rank
       ranges(2, 1) = npes_edyn3d - 1 ! Last rank                                                   
       ranges(3, 1) = 1               ! Stride of 1                                                                                                         
       ! Step 3: Create the new group using the rank range                                          
       call MPI_Group_range_incl(host_group, 1, ranges, dynamo_group, ierror)
       if (ierror /= MPI_SUCCESS) call handle_error('MPI_Group_range_incl', ierror)

       ! Step 4: Create a new communicator (dynamo_world) from the new group                        
       call MPI_Comm_create(mpi_comm_host, dynamo_group, dynamo_world, ierror)
       if (ierror /= MPI_SUCCESS) call handle_error('MPI_Comm_create', ierror)
       
       ! Check if the process is part of the new communicator                                       
       if (dynamo_world /= MPI_COMM_NULL) then
          call MPI_Comm_rank(dynamo_world, mpi_rank, ierror)
          call MPI_Comm_size(dynamo_world, mpi_size, ierror)
       else
          mpi_size = npes_edyn3d !let non-participating procs know the sizes for                   
                                 !other function all need to know lat_size and lon_size            
          mpi_rank = -1 !not participating                                                    
       endif
    else
    !just duplicate                                                                           
       call MPI_Comm_dup(mpi_comm_host, dynamo_world, ierror)
       if (ierror /= MPI_SUCCESS) call handle_error('MPI_Comm_dup', ierror)

       ! Get the rank in the new communicator                                                 
       call MPI_Comm_rank(dynamo_world, mpi_rank, ierror)
       if (ierror /= MPI_SUCCESS) call handle_error('MPI_Comm_rank', ierror)

       ! Get the size
       call MPI_Comm_size(dynamo_world, mpi_size, ierror)
       if (ierror /= MPI_SUCCESS) call handle_error('MPI_Comm_size', ierror)
    endif
    
#else
    mpi_rp = rp
    mpi_size = 1
    mpi_rank = 0
#endif

    !Modified for distributed version, lon_size must be an even number
    !BEST RESULTS when mpi_size is divisible by 4. 
    ! Check if mpi_size is valid for even lon_size constraint
    if (mpi_size /= 1 .and. mpi_size /= 2 .and. mod(mpi_size, 4) /= 0) then
       write(iulog,*) 'MPI WARNING: mpi_size should be divisible by 4, or equal to 1 or 2 (OR THERE WILL BE PROBLEMS)'
       write(iulog,*) 'Current mpi_size =', mpi_size
    endif

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
          if (lon_size*lat_size == mpi_size .and. mod(lon_size,2) == 0) exit ! lon_size >= lat_size
       enddo
    endif
   
! stack along latitudes first then longitudes
! (lat_size=3)
!  8  9 10 11
!  4  5  6  7
!  0  1  2  3 (lon_size=4)
    lat_rank = mpi_rank / lon_size
    lon_rank = modulo(mpi_rank, lon_size)

  endsubroutine init
!-----------------------------------------------------------------------
  subroutine setup_topology(nmlat_in, nmlon_in)
    ! setup MPI decompositions in geo and mag coordinates and the connectivity matrix
    ! note: non-active procs have mpi_rank < 0
    
    use params_module, only:nmlat_h,nmlat_T1
    
    integer, intent(in) :: nmlat_in, nmlon_in

    integer :: i, j, rnk, rnki, rnkj, mlat0_n, mlat1_n, &
               mysize_n, mysize_s, cnt
    integer, dimension(:), allocatable :: task_mygrid_size
    
    allocate(nmlat_task(0:lat_size-1))
    allocate(nmlon_task(0:lon_size-1))

    nmlat_task = 0
    nmlon_task = 0

    ! AB: TO DO: the first 2 aren't needed for the dist version?
    ! TO DO: deallocate
    allocate(mlat0_task(0:mpi_size-1))
    allocate(mlat1_task(0:mpi_size-1))
    allocate(mlon0_task(0:mpi_size-1))
    allocate(mlon1_task(0:mpi_size-1))
    mlat0_task = 1
    mlon0_task = 1
    mlat1_task = -1
    mlon1_task = -1
 
    !for dist
    allocate(task_csr_rowstarts(0:mpi_size))
    allocate(task_csr_mapping(0:mpi_size-1))
    allocate(task_mygrid_size(0:mpi_size-1))
    allocate(task_lat_offset(0:lat_size-1))
    task_csr_rowstarts = 0
    task_csr_mapping = 0
    task_mygrid_size = 0
    task_lat_offset = 0

    !set globale vars in this module
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
    ! AB: i don't think we need these 4 arrays for the dist version 
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


    !Find partner for distributed grid (0,1), (2,3), (3,4) etc.
    ! number of procs is even
    if (mpi_size > 1) then
       if (mpi_rank >=0 ) then !active
          if (mod(mpi_rank, 2) == 0) then
             !Number is even or zero'
             mpi_partner = mpi_rank + 1
          else
             !Number is odd'
             mpi_partner = mpi_rank - 1
          endif
       else !not active
          mpi_partner = -1
       endif
    else !one proc
       mpi_partner = 0
    endif

    
    !initial matrix row start and stops
    if (mpi_rank >=0 ) then !active
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
       
       !sizes in each hemisphere
       mysize_n =  (ij_stop_n -ij_start_n + 1)
       mysize_s = (ij_stop_s -ij_start_s + 1)
       write(*,*) 'AB: SETUP TOPO mysize_s, mysize_n', mysize_s, mysize_n
       write(*, *) 'AB: TOPO ij_start_s, ij_stop_s, s_grid_pts = ', ij_start_s, ij_stop_s, mysize_s
       write(*, *) 'AB: TOPO ij_start_n, ij_stop_n, n_grid_pts = ', ij_start_n, ij_stop_n, mysize_n
       
       !set global vars (n & s row counts will be diff for procs on equator)
       mygrid_size_n = mysize_n
       mygrid_size_s = mysize_s
       
       !get partner sizes and then my grid size for block csr matrix
       !for each partner pair, even owns s hemi and odd owns north hemi
       if (mpi_size > 1) then
          if (mod(mpi_rank,2) == 0) then !EVEN, own south, send north
             partner_hgridsize = partner_exchange_int(mysize_n)
             my_hgridsize = mysize_s
             mygrid_size =  mysize_s + partner_hgridsize
             my_sendgrid_size = mysize_n !i send to my partner
             
          else !ODD, own north, send south
             partner_hgridsize = partner_exchange_int(mysize_s)
             my_hgridsize = mysize_n
             mygrid_size =  mysize_n + partner_hgridsize
             my_sendgrid_size = mysize_s !i send to my partner
             
          endif
      
          !now we need to calculate the rowstarts for the global block
          !csr martix - this will be 0-based indeing for superlu
          !do an allgather to get each procs grid size
          task_mygrid_size = all_gather_int(mygrid_size)
          
          !now task_csr_rowstarts - set all to zero
          !south hemisphere is even, then northern is odd, so for 6 tasks
          ! the order of block rows:
          !1
          !3
          !5
          !4
          !2
          !0
          ! the task_csr_mapping will tell each task it's position in the rowstarts
          ! according to above, so rank 1 needs to know that its accesses rowstarts(5)
          ! so task_csr_mapping(1)= 5
          !south (even)
          cnt = 0
          do i=0, mpi_size-1, 2
             cnt = cnt + 1
             task_csr_rowstarts(cnt) = task_csr_rowstarts(cnt-1) &
                  + task_mygrid_size(i)
             task_csr_mapping(i) = cnt - 1
          enddo
          !north 
          if (mpi_size > 1) then !mpi_size is even
             do i= mpi_size-1, 1, -2
                cnt = cnt + 1
                task_csr_rowstarts(cnt) = task_csr_rowstarts(cnt-1) &
                     + task_mygrid_size(i)
                task_csr_mapping(i) = cnt -1 
             enddo
          endif
       else !one proc
          mygrid_size = mysize_n + mysize_s
          partner_hgridsize = 0
          my_sendgrid_size = 0
          my_hgridsize = mysize_s
          task_csr_rowstarts(1) = mygrid_size
       endif
       
       ! halos
       mlatd0 = mlat0 - 1
       mlatd1 = mlat1 + 1
       mlond0 = mlon0 - 1
       mlond1 = mlon1 + 1

       write(*, *) 'AB: TOPO2 mpi_rank, mlat0, mlat1, mlon0,mlon1',mpi_rank, mlat0, mlat1, mlon0,mlon1
       
    else !non-active
       mygrid_size = 0
       partner_hgridsize = 0
    endif

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
    call MPI_Group_free(host_group, ierror)
    Call MPI_Group_free(dynamo_group, ierror)

#endif

  endsubroutine dynamo_mpi_finalize

 !-----------------------------------------------------------------------

function partner_exchange_int(intin) result(intout)

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

    
    !post receive
    call MPI_Irecv(intout, 1, MPI_INTEGER, mpi_partner, tag, &
         dynamo_world, recv_request, ierr)
    if (ierr /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierr)

    !post send
    call MPI_Isend(intin, 1, MPI_INTEGER, mpi_partner, tag, &
         dynamo_world, send_request, ierr)
    if (ierr /= MPI_SUCCESS) call handle_error('MPI_Isend', ierr)

    ! Wait for send to complete
    call MPI_Wait(send_request, MPI_STATUS_IGNORE, ierr)
    if (ierr /= MPI_SUCCESS) call handle_error('MPI_Wait', ierr)

    ! Wait for receive to complete
    call MPI_Wait(recv_request, MPI_STATUS_IGNORE, ierr)
    if (ierr /= MPI_SUCCESS) call handle_error('MPI_Wait', ierr)

    !write(*,*) 'AB: intin = ', intin, '  intout = ' , intout

#else

    intout = intin

#endif
    
endfunction partner_exchange_int

!-----------------------------------------------------------------------

subroutine partner_exchange_hemisphere_mat(nnz_per_row, my_rowptr, my_values, my_cols, &
     partner_rowptr, partner_values, partner_cols)

  !send my csr matrix info and recv my partner's matrix info
  
  !note: the south pole's denser row does not get sent (it's owned by rank 0, which is even & sends north)  

#ifdef PARALLEL
    use MPI
#endif
    integer, intent(in) :: nnz_per_row
    
    integer, dimension(my_sendgrid_size+1), intent(in) :: my_rowptr
    integer, dimension(my_sendgrid_size*nnz_per_row), intent(in) :: my_cols
    real(kind=rp), dimension(my_sendgrid_size*nnz_per_row), intent(in) :: my_values

    integer, dimension(partner_hgridsize+1), intent(out) :: partner_rowptr
    integer, dimension(partner_hgridsize*nnz_per_row), intent(out) :: partner_cols
    real(kind=rp), dimension(partner_hgridsize*nnz_per_row), intent(out) :: partner_values

    
    integer :: i, nnz_count
    
#ifdef PARALLEL

    integer ierr
    integer, dimension(3) :: send_request, recv_request

    partner_rowptr = 0
    partner_cols = 0
    partner_values = 0.0

    !send to my partner
    call MPI_Isend(my_rowptr, my_sendgrid_size + 1, MPI_INTEGER, mpi_partner, 200, &
         dynamo_world, send_request(1), ierr)
    if (ierr /= MPI_SUCCESS) call handle_error('MPI_Isend', ierr)

    call MPI_Isend(my_cols, my_sendgrid_size*nnz_per_row, MPI_INTEGER, mpi_partner, 201, dynamo_world, &
         send_request(2), ierr)
    if (ierr /= MPI_SUCCESS) call handle_error('MPI_Isend', ierr)

    call MPI_Isend(my_values, my_sendgrid_size*nnz_per_row, mpi_rp, mpi_partner, 202, dynamo_world, &
         send_request(3), ierr)
    if (ierr /= MPI_SUCCESS) call handle_error('MPI_Isend', ierr)

    !recv from my partner
    call MPI_Irecv(partner_rowptr, partner_hgridsize + 1, MPI_INTEGER, mpi_partner, &
         200, dynamo_world, recv_request(1), ierr)
    if (ierr /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierr)

    call MPI_Irecv(partner_cols,partner_hgridsize*nnz_per_row, MPI_INTEGER, mpi_partner, &
         201, dynamo_world, recv_request(2), ierr)
    if (ierr /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierr)

    call MPI_Irecv(partner_values, partner_hgridsize*nnz_per_row, mpi_rp, mpi_partner, &
         202, dynamo_world, recv_request(3), ierr)
    if (ierr /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierr)
       
    ! Wait for all sends to complete
    call MPI_Waitall(3, send_request, MPI_STATUSES_IGNORE, ierr)
    if (ierr /= MPI_SUCCESS) call handle_error('MPI_Waitall', ierr)

    !Wait for all recvs
    call MPI_Waitall(3, recv_request, MPI_STATUSES_IGNORE, ierr)
    if (ierr /= MPI_SUCCESS) call handle_error('MPI_Waitall', ierr)

    
#else
!serial
    do concurrent (i = 1: my_hgridsize + 1)
       partner_rowptr(i) = my_rowptr(i)
    enddo
    !note for nnz calc that we are 1-based?
    nnz_count = my_rowptr(my_hgridsize + 1) -1
    do concurrent (i = 1: nnz_count
       partner_cols(i) = my_cols(i)
       partner_values(i) = my_values(i)
    enddo
    
#endif
    
    
endsubroutine partner_exchange_hemisphere_mat
!-----------------------------------------------------------------------

subroutine partner_exchange_hemisphere_vec(my_values, partner_values)
!send my_values to partner and recv partner_values  

#ifdef PARALLEL
    use MPI
#endif

    real(kind=rp), dimension(my_sendgrid_size), intent(in) :: my_values
    real(kind=rp), dimension(partner_hgridsize), intent(out) :: partner_values

    
#ifdef PARALLEL

    integer :: ierr
    integer :: send_request, recv_request

    !send to my partner
   
    call MPI_Isend(my_values, my_sendgrid_size, mpi_rp, mpi_partner, 400, dynamo_world, &
         send_request, ierr)
    if (ierr /= MPI_SUCCESS) call handle_error('MPI_Isend', ierr)

    !recv from my partner
    call MPI_Irecv(partner_values, partner_hgridsize, mpi_rp, mpi_partner, &
         400, dynamo_world, recv_request, ierr)
    if (ierr /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierr)
       
    ! Wait for my send to complete
    call MPI_Wait(send_request, MPI_STATUSES_IGNORE, ierr)
    if (ierr /= MPI_SUCCESS) call handle_error('MPI_Wait', ierr)

    !Wait for my recv
    call MPI_Wait(recv_request, MPI_STATUSES_IGNORE, ierr)
    if (ierr /= MPI_SUCCESS) call handle_error('MPI_Wait', ierr)
    
#else
    !serial
    integer :: i
    
    do concurrent i =1:my_hgridsize 
       partner_values(i) = my_values(i)
    enddo
    
#endif
    
    
endsubroutine partner_exchange_hemisphere_vec
        
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

function all_gather_int(intin) result(intarrayout)

#ifdef PARALLEL
    use MPI
#endif

    integer, intent(in) :: intin

    integer, dimension(0:mpi_size-1) :: intarrayout


#ifdef PARALLEL
    integer :: ierr, cnt

    cnt = 1
    call MPI_Allgather(intin, cnt, MPI_INTEGER, &
        intarrayout, cnt, MPI_INTEGER, dynamo_world, ierr)
    if (ierr /= MPI_SUCCESS) call handle_error('MPI_Allgather', ierr)

#else

    intarrayout(0) = intin

#endif
    
endfunction all_gather_int

    
  
!-----------------------------------------------------------------------
  function gather_lon_1d(varin) result(varout)
    ! collect a 1d array from other procs w/lat_rank 0 to root proc 0
    ! this could be genearlized to have any root proc and any proc row (Or column)
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
!    write(*,*) 'AB: DEBUG - nmlon =', nmlon
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
            !    write(iulog,*) 'AB: j=', j
                varout(j) = recvbuf(j)
             end do
          end do
       endif
       
    elseif (mpi_rank > 0 .and. lat_rank == 0) then ! send info to root (0)

       cnt = mlon1-mlon0+1
       !do concurrent (i = 1:cnt)
       !   sendbuf(i) = varin(i+mlon0-1)
       !enddo
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
