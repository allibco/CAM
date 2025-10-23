module dist_spmv_mod
  use iso_c_binding, only: c_int, c_double
  use mpi
  implicit none
   !public :: dist_spmv_init
  !public :: dist_spmv, dist_spmv_free

  type halo_t
     integer :: nprocs_local      ! communicator size
     integer :: rank
     integer :: n_global !global size
     integer :: m_loc !local size
     integer :: fst_row, last_row
     integer, allocatable :: col_owners(:)! owners for each halo col (size nhalo)
     integer, allocatable :: halo_cols(:) ! global column indices of halo entries to recv (sorted, size nhalo)
     integer :: nhalo !number of halo entries that we need to recv
     integer :: nhalo_send !how many entries we need to send
     integer, allocatable :: sendcounts(:), sdispls(:)  !length nprocs
     integer, allocatable :: recvcounts(:), rdispls(:)  !length nprocs
     integer, allocatable :: recv_from(:)  ! ranks we expect data from (list)
     integer, allocatable :: send_to(:)  ! ranks we send data to (list)
     integer, allocatable :: send_cols(:) ! global indices that we need to send to others
     real(c_double), allocatable :: sendbuf(:), recvbuf(:) !for the matvecs
  end type halo_t

contains

  !----------------------------------------------------------------------
  ! dist_spmv_init
  !  Precompute halo pattern and communication groups.
  !  Input:
  !    m_loc, fst_row, n_global, ntasks, myrank,rowptr, colind, task_row_starts
  !    (task_row_starts are 0-index)
  !  Output:
  !    halo (type halo_t) - must be passed to dist_spmv and freed with dist_spmv_free
  !----------------------------------------------------------------------

  subroutine dist_spmv_init(m_loc, fst_row, n_global, ntasks, myrank, rowptr, colind, task_row_starts, comm, halo, ierr)
    integer(c_int), intent(in) :: m_loc, fst_row, n_global, ntasks, myrank
    integer(c_int), intent(in) :: rowptr(:), colind(:)
    integer(c_int), intent(in) :: task_row_starts(ntasks+1)
    integer, intent(in) :: comm
    type(halo_t), intent(out) :: halo
    integer, intent(out) :: ierr

    integer :: i, j, k, nn, rank, nprocs, jp, tag
    integer, allocatable :: tmp(:)
    integer :: owner
    integer :: send_to_size, recv_from_size, cnt, indx
    integer :: last_row
    integer :: nowners
    integer, allocatable :: requests(:)  
    integer, allocatable :: stats(:,:)

    ierr = 0
    tag = 1314

    halo%rank = myrank
    halo%nprocs_local = ntasks
    halo%n_global = n_global
    halo%fst_row = fst_row
    halo%m_loc = m_loc
    last_row = fst_row + m_loc - 1
    halo%last_row = last_row
    ! 1) Collect remote column indices (may have duplicates)
    allocate(tmp(0)) !empty array
    nn = 0
    do i = 1, m_loc
       do j = rowptr(i), rowptr(i+1)-1
          !add one to j since row_ptr is 0-based
          jp = j+1
          if (colind(jp) < fst_row .or. colind(jp) > last_row ) then !i don't own
             nn = nn + 1 !increase halo
             if (size(tmp) == 0) then
                deallocate(tmp)
                allocate(tmp(1))
                tmp(1) = colind(jp)
             else
                tmp = [tmp, colind(jp)]
             end if
          end if
       end do
    end do

    ! 2) Make unique and sort 
    if (nn > 0) then
       call unique_sort_int(tmp, nn)
    endif
    
    halo%nhalo = nn
    allocate(halo%halo_cols(nn))
    halo%halo_cols = tmp(1:nn)

    ! 3) Determine owner of each halo column i need (use task_row_starts)
    !   halo%halo_cols and task_row_starts are sorted, so we can go in order
    allocate(halo%col_owners(halo%nhalo))
    allocate(halo%recvcounts(halo%nprocs_local))
    halo%recvcounts = 0
    halo%col_owners = 0

    rank = 0 ! first rank
    ! Advance to the correct rank if halo_cols(k) has crossed next boundary
    do k = 1, halo%nhalo
       do while (rank < halo%nprocs_local -1 .and. halo%halo_cols(k) >= task_row_starts(rank+2))
          rank = rank + 1
       end do
       halo%col_owners(k) = rank !0-based owner
       !add to recvcounts for that proc
       halo%recvcounts(rank + 1) = halo%recvcounts(rank + 1) + 1  ! +1 since Fortran arrays are 1-based
    end do

    ! Determine sendcounts by sending recvcounts to each rank 
    allocate(halo%sendcounts(halo%nprocs_local))
    halo%sendcounts = 0

    !print *, 'D1: iam = ', myrank,'sendcounts:', halo%sendcounts
    !print *, 'D1: iam = ', myrank,'recvcounts:', halo%recvcounts

    
    call MPI_Alltoall(halo%recvcounts, 1, MPI_INTEGER, halo%sendcounts, 1, MPI_INTEGER, comm, ierr)
    !recvcounts is how many i need to recv from each proc
    !sendcounts is how many i need to send to each

    !print *, 'D2: iam = ', myrank,'sendcounts:', halo%sendcounts
    !print *, 'D2: iam = ', myrank,'recvcounts:', halo%recvcounts

    
    ! 4) Compute displacements for data to recv (indexes into halo_cols)
    allocate(halo%rdispls(halo%nprocs_local))
    halo%rdispls = 0
    do rank = 2, ntasks
        halo%rdispls(rank) = halo%rdispls(rank - 1) + halo%recvcounts(rank - 1)
    end do

    !calc displacements for data that we need to send 
    allocate(halo%sdispls(halo%nprocs_local))
    halo%sdispls = 0
    cnt = halo%sendcounts(1)
    do rank = 2, halo%nprocs_local
       halo%sdispls(rank) = halo%sdispls(rank-1) + halo%sendcounts(rank-1)
       cnt = cnt + halo%sendcounts(rank)
    end do
    !how many entries to send
    halo%nhalo_send = cnt
    
    
    ! Build list of ranks we will receive from (nonzero recvcounts)
    nowners = 0
    do rank = 1, halo%nprocs_local
       if (halo%recvcounts(rank) > 0) then
          nowners = nowners + 1
       end if
    end do
    allocate(halo%recv_from(nowners))
    nowners = 0
    do rank = 1, halo%nprocs_local
       if (halo%recvcounts(rank) > 0) then
          nowners = nowners + 1
          halo%recv_from(nowners) = rank-1 !0-based
       end if
    end do
    recv_from_size = nowners
    
    ! Build list of ranks we will send to (nonzero sendcounts)
    nowners = 0
    do rank = 1, halo%nprocs_local
       if (halo%sendcounts(rank) > 0) then
          nowners = nowners + 1
       end if
    end do
    allocate(halo%send_to(nowners))
    nowners = 0
    do rank = 1, halo%nprocs_local
       if (halo%sendcounts(rank) > 0) then
          nowners = nowners + 1
          halo%send_to(nowners) = rank-1 !0-based
       end if
    end do
    send_to_size = nowners

    
    !print *, 'D3: iam = ', myrank,'send_to_size =', send_to_size
    !print *, 'D3: iam = ', myrank,'recv_from_size =', recv_from_size
    !print *, 'D3: iam = ', myrank,'recv_from =', halo%recv_from
    !print *, 'D3: iam = ', myrank,'send_to =', halo%send_to
    !print *, 'D3: iam = ', myrank,'rdispls =', halo%rdispls
    !print *, 'D3: iam = ', myrank,'sdispls =', halo%sdispls
    !print *, 'D3: iam = ', myrank,'halo_cols =', halo%halo_cols
    !print *, 'D3: iam = ', myrank,'nhalo_send =', halo%nhalo_send
    !print *, 'D3: iam = ', myrank,'nhalo =', halo%nhalo

    
    !need to do a communication to get the indices to send (so i send what i need to recv in halo_cols)
    allocate(requests(recv_from_size + send_to_size))
    allocate(stats(MPI_STATUS_SIZE, recv_from_size + send_to_size))

    requests = MPI_REQUEST_NULL
    stats = 0

    !allocate send_cols
    cnt = 0
    do k = 1, send_to_size
       rank = halo%send_to(k) ! 0-indexed
       cnt = cnt + halo%sendcounts(rank+1)
    enddo
    allocate(halo%send_cols(cnt))

    do k = 1, recv_from_size
       rank = halo%recv_from(k) ! this is 0-based rank (so add 1 when indexing into arrays)
       if (rank < 0 .or. rank >= ntasks) then
          print*, 'BAD RANK in recv_from(',k,') = ', rank
       endif
       cnt = halo%recvcounts(rank+1)
       indx = halo%rdispls(rank+1) + 1
       !print*,'D4 SEND: iam = ', myrank,'rank =', rank, 'cnt =', cnt, 'indx = ',indx

       call MPI_Isend(halo%halo_cols(indx), cnt, MPI_INTEGER, &
            rank, tag, comm, &
            requests(k), ierr)
       if (ierr /= 0) then
          print*, 'Isend returned nonzero ierr=',ierr, 'myrank=', myrank
       endif
    enddo

    do k = 1, send_to_size
       rank = halo%send_to(k)
       if (rank < 0 .or. rank >= ntasks) then
          print*, 'BAD RANK in recv_from(',k,') = ', rank
       endif
       cnt = halo%sendcounts(rank+1)
       indx = halo%sdispls(rank+1) + 1
       !print*,'D4 RECV: iam = ', myrank,'rank =', rank, 'cnt =', cnt, 'indx = ',indx
       call MPI_Irecv(halo%send_cols(indx), cnt, MPI_INTEGER, &
            rank, tag, comm, &
            requests(recv_from_size + k), ierr)
       if (ierr /= 0) then
          print*, 'Irecv returned nonzero ierr=',ierr, 'myrank=', myrank
       endif
    enddo

    !print*,'DID SEND & RECV: iam = ', myrank
    
    call MPI_Waitall(recv_from_size+send_to_size, requests, stats, ierr)

    !print *, 'D5: iam = ', myrank,'send_cols =', halo%send_cols

    
    !print*,'DID WAITALL: iam = ', myrank

    !allocate the sendbuf and recvbuf here so they can be reusued
    ! Build send buffer:for the data we will send to other procs
    allocate(halo%sendbuf(halo%nhalo_send))
    !allocate recv buf - for getting our halo data from other procs
    allocate(halo%recvbuf(halo%nhalo))
    
    
    ! Clean up
    deallocate(tmp, requests, stats)

  end subroutine dist_spmv_init

  !----------------------------------------------------------------------
  ! dist_spmv - perform y_local = A_local * x_global
  !   Inputs:
  !      rowptr, colind, nzval  - local CSR (rowptr(1:m_loc+1), colind(1:nnz), nzval(1:nnz))
  !      x_local(1:m_loc)       - local portion of x (corresponds to fst_row..last_row)
  !      halo (precomputed in the init)
  !   Output:
  !      y_local(1:m_loc)
  ! NOTE: assumes that col_ind within each row are sorted!!!!
  !----------------------------------------------------------------------

  subroutine dist_spmv(rowptr, colind, nzval, x_local, y_local, halo, comm, ierr)
    integer(c_int), intent(in) :: rowptr(:), colind(:)
    real(c_double), intent(in) :: nzval(:)
    real(c_double), intent(in) :: x_local(:)
    real(c_double), intent(out) :: y_local(:)
    type(halo_t), intent(inout) :: halo
    integer, intent(in) :: comm
    integer, intent(out) :: ierr

    integer :: i, j, jp, p, k, id, starter_k
    integer :: m_loc, first_row, loc
    integer :: nprocs
    integer, parameter :: dp_kind = c_double
    integer :: reqs_count, total_reqs
    integer, allocatable :: reqs(:)
    integer, allocatable :: stats(:,:)
    integer :: sdis, rdis
    integer :: num_recvs, num_sends
    integer :: tag, owner, myrank

    ierr = 0
    tag = 1315
    
    m_loc = halo%m_loc
    nprocs = halo%nprocs_local
    first_row = halo%fst_row
    call MPI_Comm_rank(comm, myrank, ierr)

    ! Zero output
    y_local = 0.0d0

    
    ! Build send buffer: extract x values for halo_cols, packed in send_order
    halo%sendbuf = 0.0
    do p = 1, halo%nhalo_send
       id = halo%send_cols(p) ! this is the global col id to send (0-based)
       halo%sendbuf(p) = x_local(id - first_row + 1) 
    end do

    ! Post Irecv from each src rank where recvcounts>0
    num_recvs = size(halo%recv_from)
    num_sends = size(halo%send_to)

    !allocate stats and requests
    total_reqs = num_recvs + num_sends
    
    allocate(reqs(total_reqs))
    allocate(stats(MPI_STATUS_SIZE, total_reqs))
    reqs_count = 0

    halo%recvbuf =0.0
    
    ! Irecv
    do k = 1, num_recvs
       owner = halo%recv_from(k) !0-based rank
       rdis = halo%rdispls(owner+1) ! 0-based displacement
       call MPI_Irecv(halo%recvbuf(rdis+1), halo%recvcounts(owner+1), MPI_DOUBLE_PRECISION, owner, &
                      tag, comm, reqs(reqs_count+1), ierr)
       reqs_count = reqs_count + 1
    end do

    ! Isend
    do k = 1, num_sends
       owner = halo%send_to(k) !0-based rank
       sdis = halo%sdispls(owner+1) !0-based
       call MPI_Isend(halo%sendbuf(sdis+1), halo%sendcounts(owner+1), MPI_DOUBLE_PRECISION, owner, &
            tag, comm, reqs(reqs_count+1), ierr)
       reqs_count = reqs_count + 1
    end do


    !we should do our diag block (local) multiply part here
   

    ! Waitall
    if (reqs_count > 0) then
       call MPI_Waitall(reqs_count, reqs, stats, ierr)
    end if

    !print *, 'MATVEC: iam = ', myrank,'SENDbuf =', halo%sendbuf
    !print *, 'MATVEC: iam = ', myrank,'RECVbuf =', halo%recvbuf

   ! for each owner O, the entries I receive from O correspond to
   ! those halo_cols i have whose owner == O, and in the same order as they appear in halo%send_cols for owner O.
   ! which matches my halo_cols by construction

    !print*,'IN spmv: iam = ', myrank, 'm_loc = ', m_loc

    ! Finally do local SpMV using halo_values when needed
    do i = 1, m_loc ! for each row
       starter_k = 1
       do j = rowptr(i), rowptr(i+1)-1 !for each col ind
          jp = j+1 !becuz rowptr is 0-based
          if (colind(jp) >= halo%fst_row .and. colind(jp) <= halo%last_row) then
             loc =  colind(jp) - halo%fst_row + 1 !colind is 0-bases 
             y_local(i) = y_local(i) + nzval(jp) * x_local(loc)
             !print*, 'S:iam = ', myrank, 'local i, jp', i,  jp
          else
             ! find index into halo_cols
             ! linear search; can be replaced with hash if halo large
             ! but since halo_cols and col_ind are sorted, then linear
             ! for each row, so shouldn't be too bad
             do k = starter_k, halo%nhalo
                if (halo%halo_cols(k) == colind(jp)) then
                   y_local(i) = y_local(i) + nzval(jp) * halo%recvbuf(k)
                   !print*, 'S:iam = ', myrank, 'NOT local i, jp', i,  jp
                   starter_k = k+1
                   exit
                end if
             end do
          end if
       end do
    end do

    !print*,'IN spmv: iam = ', myrank, 'y_local = ', y_local

    ! cleanup temporaries !FINISH
    deallocate(reqs, stats)

   end subroutine dist_spmv

  
  !----------------------------------------------------------------------
  ! dist_spmv_free - free arrays inside halo
  !----------------------------------------------------------------------

  subroutine dist_spmv_free(halo)
    type(halo_t), intent(inout) :: halo
    if (allocated(halo%halo_cols)) deallocate(halo%halo_cols)
    if (allocated(halo%col_owners)) deallocate(halo%col_owners)
    if (allocated(halo%sendcounts)) deallocate(halo%sendcounts)
    if (allocated(halo%sdispls)) deallocate(halo%sdispls)
    if (allocated(halo%recvcounts)) deallocate(halo%recvcounts)
    if (allocated(halo%rdispls)) deallocate(halo%rdispls)
    if (allocated(halo%send_cols)) deallocate(halo%send_cols)
    if (allocated(halo%send_to)) deallocate(halo%send_to)
    if (allocated(halo%recv_from)) deallocate(halo%recv_from)
    if (allocated(halo%recvbuf)) deallocate(halo%recvbuf)
    if (allocated(halo%sendbuf)) deallocate(halo%sendbuf)
  end subroutine dist_spmv_free

  !----------------------------------------------------------------------
  ! Helper: return x value for a global index using local x_local if owned
  !----------------------------------------------------------------------

  !----------------------------------------------------------------------
  ! small utility: unique_sort_int 
  !----------------------------------------------------------------------
  subroutine unique_sort_int(arr, n)
    !! Sorts arr(1:n) in ascending order and removes duplicates in-place.
    !! On return, n = number of unique elements (<= original n).

    implicit none
    integer, intent(inout) :: n
    integer, intent(inout) :: arr(:)

    if (n <= 1) return
    call quicksort_int(arr, 1, n)
    
    ! Remove duplicates in-place
    call unique_inplace_int(arr, n)
  end subroutine unique_sort_int


  recursive subroutine quicksort_int(a, left, right)
    implicit none
    integer, intent(inout) :: a(:)
    integer, intent(in) :: left, right
    integer :: i, j, pivot, tmp

    if (left >= right) return
    pivot = a((left + right) / 2)
    i = left
    j = right
    do
       do while (a(i) < pivot)
          i = i + 1
       end do
       do while (a(j) > pivot)
          j = j - 1
       end do
       if (i <= j) then
          tmp = a(i)
          a(i) = a(j)
          a(j) = tmp
          i = i + 1
          j = j - 1
       end if
       if (i > j) exit
    end do

    if (left < j) call quicksort_int(a, left, j)
    if (i < right) call quicksort_int(a, i, right)
  end subroutine quicksort_int
  
  
  subroutine unique_inplace_int(a, n)
    !! Removes duplicates in sorted array a(1:n), in-place.
    !! Returns new n = number of unique values.
    implicit none
    integer, intent(inout) :: n
    integer, intent(inout) :: a(:)
    integer :: i, k
    
    if (n <= 1) return
    
    k = 1
    do i = 2, n
       if (a(i) /= a(k)) then
          k = k + 1
          a(k) = a(i)
       end if
    end do
    n = k
  end subroutine unique_inplace_int

  
end module dist_spmv_mod
