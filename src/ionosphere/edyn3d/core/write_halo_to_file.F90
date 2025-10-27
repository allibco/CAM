subroutine write_halo_to_file(halo, comm)
  use iso_c_binding
  use mpi
  implicit none
  type(halo_t), intent(in) :: halo
  integer, intent(in) :: comm
  integer :: rank, ierr, unitno, nprint
  character(len=128) :: filename

  call MPI_Comm_rank(comm, rank, ierr)

  write(filename, '(A,I0,A)') 'halo_rank', rank, '.txt'
  open(newunit=unitno, file=filename, status='replace', action='write', iostat=ierr)
  if (ierr /= 0) then
     print *, 'Rank', rank, 'could not open', trim(filename)
     return
  end if

  write(unitno, '(A,I0)') '=== HALO INFO (rank ', rank, ') ==='
  write(unitno, '(A,I0)') 'nprocs_local = ', halo%nprocs_local
  write(unitno, '(A,I0)') 'n_global     = ', halo%n_global
  write(unitno, '(A,I0)') 'm_loc        = ', halo%m_loc
  write(unitno, '(A,2I10)') 'fst_row, last_row = ', halo%fst_row, halo%last_row
  write(unitno, '(A,2I10)') 'nhalo, nhalo_send = ', halo%nhalo, halo%nhalo_send

  !---- Helper internal subroutine for printing arrays
contains
  subroutine print_array_int(name, arr)
    character(len=*), intent(in) :: name
    integer, intent(in) :: arr(:)
    integer :: n, i
    n = size(arr)
    write(unitno, '(A,I0)') trim(name)//' size = ', n
    if (n > 0) then
       do i = 1, n
          write(unitno, '(I8,1X)', advance='no') arr(i)
          if (mod(i,10) == 0) write(unitno,*)
       end do
       write(unitno,*)
    end if
  end subroutine print_array_int

  subroutine print_array_real(name, arr)
    character(len=*), intent(in) :: name
    real(c_double), intent(in) :: arr(:)
    integer :: n, i
    n = size(arr)
    write(unitno, '(A,I0)') trim(name)//' size = ', n
    if (n > 0) then
       do i = 1, n
          write(unitno, '(ES14.6,1X)', advance='no') arr(i)
          if (mod(i,5) == 0) write(unitno,*)
       end do
       write(unitno,*)
    end if
  end subroutine print_array_real

  !---- Print all allocated arrays
  if (allocated(halo%halo_cols))   call print_array_int('halo_cols', halo%halo_cols)
  if (allocated(halo%col_owners))  call print_array_int('col_owners', halo%col_owners)
  if (allocated(halo%recvcounts))  call print_array_int('recvcounts', halo%recvcounts)
  if (allocated(halo%rdispls))     call print_array_int('rdispls', halo%rdispls)
  if (allocated(halo%sendcounts))  call print_array_int('sendcounts', halo%sendcounts)
  if (allocated(halo%sdispls))     call print_array_int('sdispls', halo%sdispls)
  if (allocated(halo%recv_from))   call print_array_int('recv_from', halo%recv_from)
  if (allocated(halo%send_to))     call print_array_int('send_to', halo%send_to)
  if (allocated(halo%send_cols))   call print_array_int('send_cols', halo%send_cols)
  if (allocated(halo%sendbuf))     call print_array_real('sendbuf', halo%sendbuf)
  if (allocated(halo%recvbuf))     call print_array_real('recvbuf', halo%recvbuf)

  close(unitno)
end subroutine write_halo_to_file
