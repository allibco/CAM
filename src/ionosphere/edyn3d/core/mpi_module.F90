#define PARALLEL
module mpi_module

  use prec, only: rp
  use dist_solver_module, only:  calc_grid_ij
#ifdef PARALLEL
  use MPI
  use iso_fortran_env, only: real32,real64
#endif

  implicit none

  integer :: dynamo_world=-huge(1), &
    mpi_rp=-huge(1), mpi_size=0, mpi_rank=-1, &
    lat_size=0, lon_size=0, lat_rank=-1, lon_rank=-1, &
    nmlat=0, maxmlat=-1, mlat0=1, mlat1=0, mlatd0=1, mlatd1=0, &
    nmlon=0, maxmlon=-1, mlon0=1, mlon1=0, mlond0=1, mlond1=0, &
    ij_start_s=0, ij_stop_s=0, ij_start_n=0, ij_stop_n=0
  integer, dimension(:), allocatable :: &
    nmlat_task, mlat0_task, mlat1_task, &
    nmlon_task, mlon0_task, mlon1_task, &
    task_lat_offset, task_lon_offset

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
    integer :: color, npes_host

    if (rp == real32) then
      mpi_rp = MPI_REAL4
    elseif (rp == real64) then
      mpi_rp = MPI_REAL8
    else
      stop 'unknown real precision'
    endif

    call mpi_comm_size(mpi_comm_host, npes_host, ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Comm_size', ierror)

    mpi_size = min(npes_host,npes_edyn3d)

    call MPI_Comm_rank(mpi_comm_host, mpi_rank, ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Comm_rank', ierror)

    color = mpi_rank/mpi_size
    call mpi_comm_split(mpi_comm_host, color, mpi_rank, dynamo_world, ierror)

#else
    mpi_rp = rp
    mpi_size = 1
    mpi_rank = 0
#endif

! factorize MPI process number to the nearest two numbers
    do lat_size = int(sqrt(real(mpi_size, kind=rp))), 1, -1
      lon_size = mpi_size / lat_size
      if (lon_size*lat_size == mpi_size) exit ! lon_size >= lat_size
    enddo
!!$    do lon_size = int(sqrt(real(mpi_size, kind=rp))), 1, -1
!!$      lat_size = mpi_size / lon_size
!!$      if (lat_size*lon_size == mpi_size) exit ! lon_size >= lat_size
!!$    enddo

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

    integer, intent(in) :: nmlat_in, nmlon_in

    integer :: i, j, rnk, rnki, rnkj, mlat0_n, mlat1_n

    allocate(nmlat_task(0:lat_size-1))
    allocate(nmlon_task(0:lon_size-1))

    nmlat_task = 0
    nmlon_task = 0

    allocate(mlat0_task(0:mpi_size-1))
    allocate(mlat1_task(0:mpi_size-1))
    allocate(mlon0_task(0:mpi_size-1))
    allocate(mlon1_task(0:mpi_size-1))

    allocate(task_lat_offset(0:lat_size-1))
    allocate(task_lon_offset(0:lon_size-1))

    mlat0_task = 1
    mlon0_task = 1
    mlat1_task = -1
    mlon1_task = -1

    task_lat_offset = 0
    task_lon_offset = 0
    
    nmlat = nmlat_in
    nmlon = nmlon_in

! setup magnetic grid decomposition
! each process can have unequal number of latitudes or longitudes

    nmlat_task = generate_minvar_list(nmlat, lat_size)
    maxmlat = maxval(nmlat_task)
    if (lat_rank<lat_size) then
       mlat0 = 1
       do j = 0, lat_rank-1
          mlat0 = mlat0 + nmlat_task(j)
       enddo
       mlat1 = mlat0 + nmlat_task(lat_rank) - 1
    endif

    nmlon_task = generate_minvar_list(nmlon, lon_size)
    maxmlon = maxval(nmlon_task)
    if (lat_rank<lat_size) then
       mlon0 = 1
       do i = 0, lon_rank-1
          mlon0 = mlon0 + nmlon_task(i)
       enddo
       mlon1 = mlon0 + nmlon_task(lon_rank) - 1
    endif

! each process keeps a record of the lat-lon decomposition
    do concurrent (rnk = 0:mpi_size-1)
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
    !calc offset for lon to faciliate matrix creation
    !TO DO: we might not need this
    do j = 1, lon_size-1
       task_lon_offset(j) = task_lon_offset(j-1) + nmlon_task(j-1)       
    enddo
    
    !matrix row start and stops
    !s hemi
    ij_start_s = calc_grid_ij(mlon0,mlat0,lat_rank)
    ij_stop_s = calc_grid_ij(mlon1,mlat1,lat_rank)
    !n hemi
    !adjust j for n hemisphere
    mlat0_n = nmlat_T1 - mlat0 + 1
    mlat1_n =  nmlat_T1 - mlat1 + 1
    if (mlat1_n == nmlat-h) then !equator
       mlat1_n  = mlat1_n + 1
    endif
    !now mlat0_n will be bigger than mlat1_n in north hemisphere
    ij_start_n = calc_grid_ij(mlon0,mlat1_n,lat_rank)
    ij_stop_n = calc_grid_ij(mlon1,mlat0_n,lat_rank)

       
! halos
    mlatd0 = mlat0 - 1
    mlatd1 = mlat1 + 1
    mlond0 = mlon0 - 1
    mlond1 = mlon1 + 1

  endsubroutine setup_topology
!-----------------------------------------------------------------------
  subroutine finalize

#ifdef PARALLEL
    use MPI

    integer :: ierror

    call MPI_Finalize(ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Finalize', ierror)
#endif

  endsubroutine finalize
!-----------------------------------------------------------------------
  subroutine sync_mlat_5d(var, l, m, n)
! longitude halo points are not included

#ifdef PARALLEL
    use MPI
#endif

    integer, intent(in) :: l, m, n
    real(kind=rp), dimension(l, m, n, mlatd0:mlatd1, mlon0:mlon1), intent(inout) :: var

#ifdef PARALLEL
    integer :: below, above, cnt, i, lc, mc, nc, ierror
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
      below, 0, dynamo_world, request(1), ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Isend', ierror)

    call MPI_Isend(send_to_above, cnt, mpi_rp, &
      above, 1, dynamo_world, request(2), ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Isend', ierror)

    call MPI_Irecv(recv_from_above, cnt, mpi_rp, &
      above, 0, dynamo_world, request(3), ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierror)

    call MPI_Irecv(recv_from_below, cnt, mpi_rp, &
      below, 1, dynamo_world, request(4), ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierror)

! wait for sync to complete
    call MPI_Waitall(4, request, MPI_STATUSES_IGNORE, ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Waitall', ierror)

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
    integer :: left, right, cnt, ierror
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
      left, 0, dynamo_world, request(1), ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Isend', ierror)

    call MPI_Isend(send_to_right, cnt, mpi_rp, &
      right, 1, dynamo_world, request(2), ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Isend', ierror)

    call MPI_Irecv(recv_from_right, cnt, mpi_rp, &
      right, 0, dynamo_world, request(3), ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierror)

    call MPI_Irecv(recv_from_left, cnt, mpi_rp, &
      left, 1, dynamo_world, request(4), ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierror)

! wait for sync to complete
    call MPI_Waitall(4, request, MPI_STATUSES_IGNORE, ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Waitall', ierror)

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
  function gather_lon_1d(varin) result(varout)
    ! collect a 1d array from other procs w/lat_rank 0 to root proc 0
    
#ifdef PARALLEL
    use MPI
#endif

    real(kind=rp), dimension(mlon0:mlon1), intent(in) :: varin
    real(kind=rp), dimension(nmlon) :: varout

    integer :: i

#ifdef PARALLEL
    integer :: error, cnt, myrequest, tag, rs, re
    integer, dimension(1:lon_size-1) :: requests  
    real(kind=rp), dimension(maxmlon) :: sendbuf
    real(kind=rp), dimension(nmlon) :: recvbuf
 
    tag = 34
    
    ! load to work array
    sendbuf=0.0
   
    ! gather data to 0 from proc_row 0
    if (mpi_rank == 0 .and. lon_size > 1 ) then ! receive from other procs in row

       do i = 1:lon_size -1
          rs = mlon0_task(i)
          re = mlon1_task(i)
          cnt = re-rs+1
          
          call MPI_Irecv(recvbuf(rs:re), cnt, mpi_rp, &
               i, tag, dynamo_world, &
               requests(i), ierror)
          if (ierror /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierror)
       enddo


       !now wait to receive all data
       call MPI_Waitall(lon_size-1, requests, MPI_STATUSES_IGNORE, ierror)
       if (ierror /= MPI_SUCCESS) call handle_error('MPI_Waitall', ierror)

       
    elseif (mpi_rank > 0 .and. lat_rank == 0) then ! send info to root

       cnt = mlon1-mlon0+1
       do concurrent (i = 1:cnt)
          sendbuf(i) = varin(i+mlon0-1)
       enddo
       
       call MPI_Isend(sendbuf, cnt, mpi_rp, &
            root, tag, dynamo_world, &
            myrequest, ierror)
       if (ierror /= MPI_SUCCESS) call handle_error('MPI_Isend', ierror)


       call MPI_WAIT(myrequest, MPI_STATUS_IGNORE, ierror)
       if (ierror /= MPI_SUCCESS) call handle_error('MPI_Wait', ierror)

    endif

    
    !now copy to output
    !for root 0
    do concurrent (i = mlon0:mlon1)
       varout(i) = varin(i)
    enddo
    if (lon_size > 1) then
       is = mlon0_task(1)
       do concurrent (i = is:nmlon)
          varout(i) = recvbuf(i)
       enddo
    endif
#else
    do concurrent (i = mlon0:mlon1)
       varout(i) = varin(i)
    enddo
#endif
    

  
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
    integer :: cnt, rnki, i0, i1, ierror
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
        request(rnki), ierror)
      if (ierror /= MPI_SUCCESS) call handle_error('MPI_Isend', ierror)

      call MPI_Irecv(recvbuf(:, :, :, rnki), cnt, mpi_rp, &
        lat_rank*lon_size + rnki, 0, dynamo_world, &
        request(lon_size+rnki), ierror)
      if (ierror /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierror)
    enddo

! wait for gather to complete
    call MPI_Waitall(lon_size*2, request, MPI_STATUSES_IGNORE, ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Waitall', ierror)

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
    integer :: cnt, rnk, i0, i1, j0, j1, ierror
    real(kind=rp), dimension(maxmlat, maxmlon) :: sendbuf
    real(kind=rp), dimension(maxmlat, maxmlon, 0:mpi_size-1) :: recvbuf

    cnt = maxmlat * maxmlon

! load to work array
    do concurrent (i = 1:mlon1-mlon0+1, j = 1:mlat1-mlat0+1)
      sendbuf(j, i) = varin(j+mlat0-1, i+mlon0-1)
    enddo

    if (root < 0) then
      call MPI_Allgather(sendbuf, cnt, mpi_rp, &
        recvbuf, cnt, mpi_rp, dynamo_world, ierror)
      if (ierror /= MPI_SUCCESS) call handle_error('MPI_Allgather', ierror)
    else
      call MPI_Gather(sendbuf, cnt, mpi_rp, &
        recvbuf, cnt, mpi_rp, root, dynamo_world, ierror)
      if (ierror /= MPI_SUCCESS) call handle_error('MPI_Gather', ierror)
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
    integer :: cnt, rnk, i0, i1, j0, j1, ierror
    real(kind=rp), dimension(n, maxmlat, maxmlon) :: sendbuf
    real(kind=rp), dimension(n, maxmlat, maxmlon, 0:mpi_size-1) :: recvbuf

    cnt = n * maxmlat * maxmlon

! load to work array
    do concurrent (i = 1:mlon1-mlon0+1, j = 1:mlat1-mlat0+1, nc = 1:n)
      sendbuf(nc, j, i) = varin(nc, j+mlat0-1, i+mlon0-1)
    enddo

    if (root < 0) then
      call MPI_Allgather(sendbuf, cnt, mpi_rp, &
        recvbuf, cnt, mpi_rp, dynamo_world, ierror)
      if (ierror /= MPI_SUCCESS) call handle_error('MPI_Allgather', ierror)
    else
      call MPI_Gather(sendbuf, cnt, mpi_rp, &
        recvbuf, cnt, mpi_rp, root, dynamo_world, ierror)
      if (ierror /= MPI_SUCCESS) call handle_error('MPI_Gather', ierror)
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
    integer :: cnt, rnk, i0, i1, j0, j1, ierror
    real(kind=rp), dimension(m, n, maxmlat, maxmlon) :: sendbuf
    real(kind=rp), dimension(m, n, maxmlat, maxmlon, 0:mpi_size-1) :: recvbuf

    cnt = m * n * maxmlat * maxmlon

! load to work array
    do concurrent (i = 1:mlon1-mlon0+1, j = 1:mlat1-mlat0+1, nc = 1:n, mc = 1:m)
      sendbuf(mc, nc, j, i) = varin(mc, nc, j+mlat0-1, i+mlon0-1)
    enddo

    if (root < 0) then
      call MPI_Allgather(sendbuf, cnt, mpi_rp, &
        recvbuf, cnt, mpi_rp, dynamo_world, ierror)
      if (ierror /= MPI_SUCCESS) call handle_error('MPI_Allgather', ierror)
    else
      call MPI_Gather(sendbuf, cnt, mpi_rp, &
        recvbuf, cnt, mpi_rp, root, dynamo_world, ierror)
      if (ierror /= MPI_SUCCESS) call handle_error('MPI_Gather', ierror)
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
    integer :: cnt, rnk, i0, i1, j0, j1, ierror
    real(kind=rp), dimension(l, m, n, maxmlat, maxmlon) :: sendbuf
    real(kind=rp), dimension(l, m, n, maxmlat, maxmlon, 0:mpi_size-1) :: recvbuf

    cnt = l * m * n * maxmlat * maxmlon

! load to work array
    do concurrent (i = 1:mlon1-mlon0+1, j = 1:mlat1-mlat0+1, nc = 1:n, mc = 1:m, lc = 1:l)
      sendbuf(lc, mc, nc, j, i) = varin(lc, mc, nc, j+mlat0-1, i+mlon0-1)
    enddo

    if (root < 0) then
      call MPI_Allgather(sendbuf, cnt, mpi_rp, &
        recvbuf, cnt, mpi_rp, dynamo_world, ierror)
      if (ierror /= MPI_SUCCESS) call handle_error('MPI_Allgather', ierror)
    else
      call MPI_Gather(sendbuf, cnt, mpi_rp, &
        recvbuf, cnt, mpi_rp, root, dynamo_world, ierror)
      if (ierror /= MPI_SUCCESS) call handle_error('MPI_Gather', ierror)
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
    integer :: cnt, lc, mc, nc, ierror
    real(kind=rp), dimension(l, m, n) :: buffer

    cnt = l * m * n

! load to work array
    do concurrent (nc = 1:n, mc = 1:m, lc = 1:l)
      buffer(lc, mc, nc) = var(lc, mc, nc)
    enddo

    call MPI_Bcast(buffer, cnt, mpi_rp, root, dynamo_world, ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Bcast', ierror)

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
    integer :: ierror

    if (root < 0) then
      call MPI_Allreduce(varin, varout, n, mpi_rp, &
        MPI_SUM, dynamo_world, ierror)
      if (ierror /= MPI_SUCCESS) call handle_error('MPI_Allreduce', ierror)
    else
      call MPI_Reduce(varin, varout, n, mpi_rp, &
        MPI_SUM, root, dynamo_world, ierror)
      if (ierror /= MPI_SUCCESS) call handle_error('MPI_Reduce', ierror)
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
    integer :: resultlen, ierror

    call MPI_Error_string(errorcode, string, resultlen, ierror)
    write(6, "('MPI error encountered: ', a, ', when calling ', a, '. Finalizing...')") &
      trim(string), trim(funcname)
    call MPI_Finalize(ierror)
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
endmodule mpi_module
