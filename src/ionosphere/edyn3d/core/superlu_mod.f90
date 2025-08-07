
module superlu_mod

  use iso_c_binding

!----------------------------------------------------
  ! This module contains c bindings for superlu functions
  ! that I am using
!----------------------------------------------------

  implicit none


  type, bind(c) :: superlu_options_t
     integer(kind=c_int) :: Fact
     integer(kind=c_int) :: Equil
     integer(kind=c_int) :: ColPerm
     integer(kind=c_int) :: RowPerm 
     integer(kind=c_int) :: ReplaceTinyPivot
     integer(kind=c_int) :: IterRefine
     integer(kind=c_int) :: Trans
     integer(kind=c_int) :: SymmetricMode
     integer(kind=c_int) :: PrintStat
     ! ... more fields exist but these are the most common
  end type superlu_options_t

  
 ! Interface declarations for SuperLU_DIST functions
    interface
         ! Initialize SuperLU process grid
        subroutine superlu_gridinit(comm, nprow, npcol, grid) &
            bind(c, name='superlu_gridinit')
            use iso_c_binding
            integer(c_int), value :: comm    ! MPI communicator (converted to C)
            integer(c_int), value :: nprow   ! Number of process rows
            integer(c_int), value :: npcol   ! Number of process columns
            type(c_ptr) :: grid              ! Output: grid handle
        end subroutine


        ! Create distributed matrix A
        subroutine dCreate_CompRowLoc_Matrix_dist(A, m, n, nnz_loc, m_loc, &
                                                   fst_row, nzval, colind, rowptr, &
                                                   stype, dtype, mtype) &
            bind(c, name='dCreate_CompRowLoc_Matrix_dist')
            use iso_c_binding
            type(c_ptr) :: A                    ! Output: matrix handle
            integer(c_int), value :: m, n       ! Global matrix dimensions
            integer(c_int), value :: nnz_loc    ! Local non-zeros
            integer(c_int), value :: m_loc      ! Local rows
            integer(c_int), value :: fst_row    ! First row (1-based)
            type(c_ptr), value :: nzval         ! Pointer to values
            type(c_ptr), value :: colind        ! Pointer to column indices  
            type(c_ptr), value :: rowptr        ! Pointer to row pointers
            integer(c_int), value :: stype      ! Storage type
            integer(c_int), value :: dtype      ! Data type
            integer(c_int), value :: mtype      ! Matrix type
        end subroutine

       
        ! Set default options
        subroutine set_default_options_dist(options) &
            bind(c, name='set_default_options_dist')
            use iso_c_binding
            type(c_ptr) :: options
        end subroutine

        ! Initialize scale/permutation structure
        subroutine dScalePermstructInit(m, n, ScalePermstruct) &
            bind(c, name='dScalePermstructInit')
            use iso_c_binding
            integer(c_int), value :: m, n
            type(c_ptr) :: ScalePermstruct
        end subroutine

        ! Initialize LU structure
        subroutine dLUstructInit(n, LUstruct) &
            bind(c, name='dLUstructInit')
            use iso_c_binding
            integer(c_int), value :: n
            type(c_ptr) :: LUstruct
        end subroutine

        ! Initialize statistics
        subroutine PStatInit(stat) &
            bind(c, name='PStatInit')
            use iso_c_binding
            type(c_ptr) :: stat
          end subroutine PStatInit
          
        ! Main solver routine
        subroutine pdgssvx(options, A, ScalePermstruct, X, ldx, nrhs, grid, &
                          LUstruct, berr, stat, info) &
            bind(c, name='pdgssvx')
            use iso_c_binding
            type(c_ptr), value :: options, A, ScalePermstruct, X, grid, LUstruct, stat
            integer(c_int), value :: ldx, nrhs
            type(c_ptr), value :: berr
            integer(c_int) :: info
        end subroutine

        subroutine PStatPrint(options, stat, grid) &
            bind(c, name="PStatPrint")
            use iso_c_binding
            type(superlu_options_t) :: options
            type(c_ptr)     :: stat
            type(c_ptr)        :: grid
        end subroutine PStatPrint

        ! Cleanup functions
        subroutine superlu_gridexit(grid) &
            bind(c, name='superlu_gridexit')
            use iso_c_binding
            type(c_ptr), value :: grid
        end subroutine

        ! Destroy SuperLU distributed matrix
        subroutine Destroy_SuperMatrix_Store_dist(A) &
            bind(c, name='Destroy_SuperMatrix_Store_dist')
            use iso_c_binding
            type(c_ptr), value :: A
        end subroutine

        subroutine dDestroy_LU(n, grid, LUstruct) &
            bind(c, name="dDestroy_LU")
            use iso_c_binding
            integer(c_int), value :: n
            type(c_ptr)      :: grid
            type(c_ptr)     :: LUstruct
        end subroutine dDestroy_LU
        
        subroutine dScalePermstructFree(ScalePermstruct) &
            bind(c, name='dScalePermstructFree')
            use iso_c_binding
            type(c_ptr), value :: ScalePermstruct
        end subroutine

        subroutine dLUstructFree(LUstruct) &
            bind(c, name='dLUstructFree')
            use iso_c_binding
            type(c_ptr), value :: LUstruct
        end subroutine

        subroutine PStatFree(stat) &
            bind(c, name='PStatFree')
            use iso_c_binding
            type(c_ptr), value :: stat
        end subroutine

     end interface
  

end module superlu_mod



