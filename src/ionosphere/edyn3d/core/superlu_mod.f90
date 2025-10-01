
module superlu_mod

  use iso_c_binding
  implicit none

!----------------------------------------------------
  ! This module contains c bindings for superlu_dist functions
  ! that I am using - version 9.1.0, use superlu_defs.h and superlu_ddefs.h
!----------------------------------------------------

 
  type, bind(C) :: superlu_dist_options_t
    integer(C_INT) :: Fact              ! fact_t
    integer(C_INT) :: Equil             ! yes_no_t
    integer(C_INT) :: DiagInv           ! yes_no_t
    integer(C_INT) :: ColPerm           ! colperm_t
    integer(C_INT) :: Trans             ! trans_t
    integer(C_INT) :: IterRefine        ! IterRefine_t
    real(C_DOUBLE) :: DiagPivotThresh
    integer(C_INT) :: SymmetricMode     ! yes_no_t
    integer(C_INT) :: PivotGrowth       ! yes_no_t
    integer(C_INT) :: ConditionNumber   ! yes_no_t
    integer(C_INT) :: RowPerm           ! rowperm_t
    integer(C_INT) :: SolveOnly         ! yes_no_t
    integer(C_INT) :: ILU_level
    integer(C_INT) :: ILU_DropRule
    real(C_DOUBLE) :: ILU_DropTol
    real(C_DOUBLE) :: ILU_FillFactor
    integer(C_INT) :: ILU_Norm          ! norm_t
    real(C_DOUBLE) :: ILU_FillTol
    integer(C_INT) :: ILU_MILU          ! milu_t
    real(C_DOUBLE) :: ILU_MILU_Dim
    integer(C_INT) :: UserDefineSupernode ! yes_no_t
    integer(C_INT) :: ParSymbFact         ! yes_no_t
    integer(C_INT) :: ReplaceTinyPivot    ! yes_no_t
    integer(C_INT) :: SolveInitialized    ! yes_no_t
    integer(C_INT) :: RefineInitialized   ! yes_no_t
    integer(C_INT) :: PrintStat           ! yes_no_t
    integer(C_INT) :: lookahead_etree     ! yes_no_t
    integer(C_INT) :: num_lookaheads
    integer(C_INT) :: superlu_relax
    integer(C_INT) :: superlu_maxsup
    character(kind=C_CHAR) :: superlu_rankorder(4)
    character(kind=C_CHAR) :: superlu_lbs(4)
    integer(C_INT) :: superlu_n_gemm
    integer(C_INT) :: superlu_max_buffer_size
    integer(C_INT) :: superlu_num_gpu_streams
    integer(C_INT) :: superlu_acc_offload
    integer(C_INT) :: batchCount
    integer(C_INT) :: SymPattern          ! yes_no_t
    integer(C_INT) :: Use_TensorCore      ! yes_no_t
    integer(C_INT) :: Algo3d              ! yes_no_t
  end type superlu_dist_options_t


  type, bind(C) :: ScalePermstruct_t
     integer(C_INT) :: DiagScale   ! DiagScale_t enum
     type(C_PTR)    :: R           ! double*
     type(C_PTR)    :: C           ! double*
     type(C_PTR)    :: perm_r      ! int_t* (maps to C int*)
     type(C_PTR)    :: perm_c      ! int_t*
  end type ScalePermstruct_t
 
  type, bind(C) :: dLUstruct_t
    type(C_PTR) :: etree         ! int_t*
    type(C_PTR) :: Glu_persist   ! Glu_persist_t*
    type(C_PTR) :: Llu           ! dLocalLU_t*
    type(C_PTR) :: trf3Dpart     ! dtrf3Dpartition_t*
    character(kind=C_CHAR) :: dt ! single char (datatype tag: 'd','s','c','z')
 end type dLUstruct_t
 
 type, bind(C) :: SuperLUStat_t
    type(C_PTR)    :: panel_histo    ! int*
    type(C_PTR)    :: utime          ! double*
    type(C_PTR)    :: ops            ! flops_t* (double*)
    integer(C_INT) :: TinyPivots
    integer(C_INT) :: RefineSteps
    integer(C_INT) :: num_look_aheads
    real(C_FLOAT)  :: current_buffer
    real(C_FLOAT)  :: peak_buffer
    real(C_FLOAT)  :: gpu_buffer
    integer(C_LONG) :: MaxActiveBTrees  ! int_t maps to long on most builds
    integer(C_LONG) :: MaxActiveRTrees
    ! GPU_ACC fields would go here if you enabled GPU support
  end type SuperLUStat_t
  
  
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
            type(c_ptr):: A                    ! Output: matrix handle
            integer(c_int), value :: m, n       ! Global matrix dimensions
            integer(c_int), value :: nnz_loc    ! Local non-zeros
            integer(c_int), value :: m_loc      ! Local rows
            integer(c_int), value :: fst_row    ! First row (o-based)
            type(c_ptr), value :: nzval         ! Pointer to values
            type(c_ptr), value :: colind        ! Pointer to column indices  
            type(c_ptr), value :: rowptr        ! Pointer to row pointers
            integer(c_int), value :: stype      ! Storage type
            integer(c_int), value :: dtype      ! Data type
            integer(c_int), value :: mtype      ! Matrix type
        end subroutine
        
        ! Set default options
        subroutine set_default_options_dist(opt) bind(C, name="set_default_options_dist")
          import :: superlu_dist_options_t
          type(superlu_dist_options_t), intent(out) :: opt
        end subroutine set_default_options_dist

        ! Initialize scale/permutation structure
        subroutine dScalePermstructInit(m, n, ScalePermstruct) &
             bind(c, name='dScalePermstructInit')
          use iso_c_binding
          import :: ScalePermstruct_t
          integer(c_int), value :: m, n
          type(ScalePermstruct_t) :: ScalePermstruct
        end subroutine

        ! Initialize LU structure
        subroutine dLUstructInit(n, LUstruct) &
             bind(c, name='dLUstructInit')
          use iso_c_binding
          import :: dLUstruct_t
          integer(c_int), value :: n
          type(dLUstruct_t) :: LUstruct
        end subroutine dLUstructInit
        
        ! Initialize statistics
        subroutine PStatInit(stat) &
            bind(c, name='PStatInit')
          use iso_c_binding
          import :: SuperLUStat_t
          type(SuperLUStat_t) :: stat
        end subroutine PStatInit

        
        ! Main solver routine
        subroutine pdgssvx(options, A, ScalePermstruct, X, ldx, nrhs, grid, LUstruct, &
                   SolveStruct, berr, stat, info) bind(c, name='pdgssvx')
          use iso_c_binding
          import :: ScalePermstruct_t, dLUstruct_t, SuperLUStat_t, superlu_dist_options_t
          type(superlu_dist_options_t) :: options         ! by reference
          type(c_ptr), value :: A, grid                   ! c_ptr to SuperMatrix and grid
          type(ScalePermstruct_t) :: ScalePermstruct      ! by reference
          type(c_ptr), value :: X                          ! c_loc(sol)
          integer(c_int), value :: ldx, nrhs
          type(dLUstruct_t) :: LUstruct                    ! by reference
          type(c_ptr), value :: SolveStruct               ! pass c_null_ptr if unused
          type(c_ptr), value :: berr                       ! c_loc(berr_array)
          type(SuperLUStat_t) :: stat                      ! by reference
          integer(c_int) :: info                            ! by reference
        end subroutine pdgssvx

        subroutine PStatPrint(options, stat, grid) bind(C, name="PStatPrint")
          use iso_c_binding
          import :: superlu_dist_options_t, SuperLUStat_t
          type(superlu_dist_options_t) :: options
          type(SuperLUStat_t)          :: stat
          type(c_ptr), value           :: grid   ! gridinfo_t*
        end subroutine PStatPrint
        
        !super lu matvec routine (internal - not typically called by user, would need the init also)
        function pdgsmv(n, A, x, ax) &
             bind(C, name="pdgsmv")
          use iso_c_binding
          integer(c_int), value :: n
          type(c_ptr), value :: A
          real(c_double) :: x(*)
          real(c_double) :: ax(*)
          integer(c_int) :: pdgsmv
        end function pdgsmv
       
          
        ! Cleanup functions
        subroutine superlu_gridexit(grid) &
            bind(c, name='superlu_gridexit')
            use iso_c_binding
            type(c_ptr) :: grid
        end subroutine

        ! Destroy SuperLU distributed matrix
        subroutine Destroy_SuperMatrix_Store_dist(A) &
            bind(c, name='Destroy_SuperMatrix_Store_dist')
            use iso_c_binding
            type(c_ptr) :: A
        end subroutine

        subroutine dDestroy_LU(n, grid, LUstruct) &
            bind(c, name="dDestroy_LU")
            use iso_c_binding
            import :: dLUstruct_t
            integer(c_int), value :: n  ! Problem size (number of columns in A)
            type(c_ptr), value    :: grid ! we treat gridinfo_t* as opaque handle
            type(dLUstruct_t)      :: LUstruct
          end subroutine dDestroy_LU

        subroutine dScalePermstructFree(ScalePermstruct) &
            bind(c, name='dScalePermstructFree')
            use iso_c_binding
            import :: ScalePermstruct_t
            type(ScalePermstruct_t) :: ScalePermstruct
        end subroutine dScalePermstructFree

        subroutine dLUstructFree(LUstruct) &
            bind(c, name='dLUstructFree')
          use iso_c_binding
          import :: dLUstruct_t
          type(dLUstruct_t) :: LUstruct
        end subroutine dLUstructFree

        subroutine PStatFree(stat) &
            bind(c, name='PStatFree')
            use iso_c_binding
            import 
            type(SuperLUStat_t) :: stat
        end subroutine PStatFree

     end interface
  

end module superlu_mod



