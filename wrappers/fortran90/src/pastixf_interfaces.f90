!>
!> @file pastixf_interfaces.f90
!>
!> PaStiX Fortran 90 wrapper
!>
!> @copyright 2017-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
!>                      Univ. Bordeaux. All rights reserved.
!>
!> @version 6.2.0
!> @author Mathieu Faverge
!> @author Tony Delarue
!> @author Selmane Lebdaoui
!> @date 2022-01-13
!>
!> This file has been automatically generated with gen_wrappers.py
!>
!> @ingroup wrap_fortran
!>
module pastixf_interfaces
  interface pastixOrderInit
     subroutine pastixOrderInit_f08(ordeptr, baseval, vertnbr, cblknbr, perm, &
          invp, rang, tree, info)
       use :: iso_c_binding, only : c_int
       use :: pastixf_enums, only : pastix_int_t, pastix_order_t
       implicit none
       type(pastix_order_t),       intent(inout), target   :: ordeptr
       integer(kind=pastix_int_t), intent(in)              :: baseval
       integer(kind=pastix_int_t), intent(in)              :: vertnbr
       integer(kind=pastix_int_t), intent(in)              :: cblknbr
       integer(kind=pastix_int_t), intent(inout), target   :: perm(:)
       integer(kind=pastix_int_t), intent(inout), target   :: invp(:)
       integer(kind=pastix_int_t), intent(inout), target   :: rang(:)
       integer(kind=pastix_int_t), intent(inout), target   :: tree(:)
       integer(kind=c_int),        intent(out),   optional :: info
     end subroutine pastixOrderInit_f08
  end interface pastixOrderInit

  interface pastixOrderAlloc
     subroutine pastixOrderAlloc_f08(ordeptr, vertnbr, cblknbr, info)
       use :: iso_c_binding, only : c_int
       use :: pastixf_enums, only : pastix_int_t, pastix_order_t
       implicit none
       type(pastix_order_t),       intent(inout), target   :: ordeptr
       integer(kind=pastix_int_t), intent(in)              :: vertnbr
       integer(kind=pastix_int_t), intent(in)              :: cblknbr
       integer(kind=c_int),        intent(out),   optional :: info
     end subroutine pastixOrderAlloc_f08
  end interface pastixOrderAlloc

  interface pastixOrderAllocId
     subroutine pastixOrderAllocId_f08(ordeptr, vertnbr, info)
       use :: iso_c_binding, only : c_int
       use :: pastixf_enums, only : pastix_int_t, pastix_order_t
       implicit none
       type(pastix_order_t),       intent(inout), target   :: ordeptr
       integer(kind=pastix_int_t), intent(in)              :: vertnbr
       integer(kind=c_int),        intent(out),   optional :: info
     end subroutine pastixOrderAllocId_f08
  end interface pastixOrderAllocId

  interface pastixOrderExit
     subroutine pastixOrderExit_f08(ordeptr)
       use :: pastixf_enums, only : pastix_order_t
       implicit none
       type(pastix_order_t), intent(inout), target :: ordeptr
     end subroutine pastixOrderExit_f08
  end interface pastixOrderExit

  interface pastixOrderBase
     subroutine pastixOrderBase_f08(ordeptr, baseval)
       use :: pastixf_enums, only : pastix_int_t, pastix_order_t
       implicit none
       type(pastix_order_t),       intent(inout), target :: ordeptr
       integer(kind=pastix_int_t), intent(in)            :: baseval
     end subroutine pastixOrderBase_f08
  end interface pastixOrderBase

  interface pastixOrderCheck
     subroutine pastixOrderCheck_f08(ordeptr, info)
       use :: iso_c_binding, only : c_int
       use :: pastixf_enums, only : pastix_order_t
       implicit none
       type(pastix_order_t), intent(in),  target   :: ordeptr
       integer(kind=c_int),  intent(out), optional :: info
     end subroutine pastixOrderCheck_f08
  end interface pastixOrderCheck

  interface pastixOrderCopy
     subroutine pastixOrderCopy_f08(ordedst, ordesrc, info)
       use :: iso_c_binding, only : c_int
       use :: pastixf_enums, only : pastix_order_t
       implicit none
       type(pastix_order_t), intent(inout), target   :: ordedst
       type(pastix_order_t), intent(in),    target   :: ordesrc
       integer(kind=c_int),  intent(out),   optional :: info
     end subroutine pastixOrderCopy_f08
  end interface pastixOrderCopy

  interface pastixOrderGet
     subroutine pastixOrderGet_f08(pastix_data, order)
       use :: pastixf_enums, only : pastix_data_t, pastix_order_t
       implicit none
       type(pastix_data_t),  intent(in),  target  :: pastix_data
       type(pastix_order_t), intent(out), pointer :: order
     end subroutine pastixOrderGet_f08
  end interface pastixOrderGet

  interface pastixOrderBcast
     subroutine pastixOrderBcast_f08(ordemesh, root, pastix_comm)
       use :: iso_c_binding, only : c_int
       use :: pastixf_enums, only : pastix_order_t
       use :: spmf_enums,    only : MPI_Comm
       implicit none
       type(pastix_order_t), intent(inout), target :: ordemesh
       integer(kind=c_int),  intent(in)            :: root
       type(MPI_Comm),       intent(in)            :: pastix_comm
     end subroutine pastixOrderBcast_f08
  end interface pastixOrderBcast

  interface pastixOrderGrid
     subroutine pastixOrderGrid_f08(myorder, nx, ny, nz, info)
       use :: iso_c_binding, only : c_int, c_ptr
       use :: pastixf_enums, only : pastix_int_t, pastix_order_t
       implicit none
       type(pastix_order_t),       intent(inout), pointer  :: myorder
       integer(kind=pastix_int_t), intent(in)              :: nx
       integer(kind=pastix_int_t), intent(in)              :: ny
       integer(kind=pastix_int_t), intent(in)              :: nz
       integer(kind=c_int),        intent(out),   optional :: info
     end subroutine pastixOrderGrid_f08
  end interface pastixOrderGrid

  interface pastixOrderLoad
     subroutine pastixOrderLoad_f08(pastix_data, ordeptr, info)
       use :: iso_c_binding, only : c_int
       use :: pastixf_enums, only : pastix_data_t, pastix_order_t
       implicit none
       type(pastix_data_t),  intent(in),    target   :: pastix_data
       type(pastix_order_t), intent(inout), target   :: ordeptr
       integer(kind=c_int),  intent(out),   optional :: info
     end subroutine pastixOrderLoad_f08
  end interface pastixOrderLoad

  interface pastixOrderSave
     subroutine pastixOrderSave_f08(pastix_data, ordeptr, info)
       use :: iso_c_binding, only : c_int
       use :: pastixf_enums, only : pastix_data_t, pastix_order_t
       implicit none
       type(pastix_data_t),  intent(inout), target   :: pastix_data
       type(pastix_order_t), intent(in),    target   :: ordeptr
       integer(kind=c_int),  intent(out),   optional :: info
     end subroutine pastixOrderSave_f08
  end interface pastixOrderSave

  interface pastix
     subroutine pastix_f08(pastix_data, pastix_comm, n, colptr, rowptr, &
          values, perm, invp, B, nrhs, iparm, dparm, info)
       use :: iso_c_binding,    only : c_double, c_int, c_ptr
       use :: pastixf_bindings, only : pastixGetCptrFrom1dArray, pastixGetCptrFrom2dArray
       use :: pastixf_enums,    only : pastix_data_t, pastix_int_t
       use :: spmf_enums,       only : MPI_Comm
       implicit none
       type(pastix_data_t),        intent(inout), pointer  :: pastix_data
       type(MPI_Comm),             intent(in)              :: pastix_comm
       integer(kind=pastix_int_t), intent(in)              :: n
       integer(kind=pastix_int_t), intent(inout), target   :: colptr(:)
       integer(kind=pastix_int_t), intent(inout), target   :: rowptr(:)
       class(*),                   intent(inout), target   :: values(:)
       integer(kind=pastix_int_t), intent(inout), target   :: perm(:)
       integer(kind=pastix_int_t), intent(inout), target   :: invp(:)
       class(*),                   intent(inout), target   :: B(:,:)
       integer(kind=pastix_int_t), intent(in)              :: nrhs
       integer(kind=pastix_int_t), intent(inout), target   :: iparm(:)
       real(kind=c_double),        intent(inout), target   :: dparm(:)
       integer(kind=c_int),        intent(out),   optional :: info
     end subroutine pastix_f08
  end interface pastix

  interface pastixInitParam
     subroutine pastixInitParam_f08(iparm, dparm)
       use :: iso_c_binding, only : c_double
       use :: pastixf_enums, only : pastix_int_t
       implicit none
       integer(kind=pastix_int_t), intent(inout), target :: iparm(:)
       real(kind=c_double),        intent(inout), target :: dparm(:)
     end subroutine pastixInitParam_f08
  end interface pastixInitParam

  interface pastixInit
     subroutine pastixInit_f08(pastix_data, pastix_comm, iparm, dparm)
       use :: iso_c_binding, only : c_double, c_ptr
       use :: pastixf_enums, only : pastix_data_t, pastix_int_t
       use :: spmf_enums,    only : MPI_Comm
       implicit none
       type(pastix_data_t),        intent(inout), pointer :: pastix_data
       type(MPI_Comm),             intent(in)             :: pastix_comm
       integer(kind=pastix_int_t), intent(inout), target  :: iparm(:)
       real(kind=c_double),        intent(inout), target  :: dparm(:)
     end subroutine pastixInit_f08
  end interface pastixInit

  interface pastixInitWithAffinity
     subroutine pastixInitWithAffinity_f08(pastix_data, pastix_comm, iparm, &
          dparm, bindtab)
       use :: iso_c_binding, only : c_double, c_int, c_ptr
       use :: pastixf_enums, only : pastix_data_t, pastix_int_t
       use :: spmf_enums,    only : MPI_Comm
       implicit none
       type(pastix_data_t),        intent(inout), pointer :: pastix_data
       type(MPI_Comm),             intent(in)             :: pastix_comm
       integer(kind=pastix_int_t), intent(inout), target  :: iparm(:)
       real(kind=c_double),        intent(inout), target  :: dparm(:)
       integer(kind=c_int),        intent(in),    target  :: bindtab(:)
     end subroutine pastixInitWithAffinity_f08
  end interface pastixInitWithAffinity

  interface pastixFinalize
     subroutine pastixFinalize_f08(pastix_data)
       use :: iso_c_binding, only : c_ptr
       use :: pastixf_enums, only : pastix_data_t
       implicit none
       type(pastix_data_t), intent(inout), pointer :: pastix_data
     end subroutine pastixFinalize_f08
  end interface pastixFinalize

  interface pastix_task_analyze
     subroutine pastix_task_analyze_f08(pastix_data, spm, info)
       use :: iso_c_binding, only : c_int
       use :: pastixf_enums, only : pastix_data_t
       use :: spmf_enums,    only : spmatrix_t
       implicit none
       type(pastix_data_t), intent(inout), target   :: pastix_data
       type(spmatrix_t),    intent(in),    target   :: spm
       integer(kind=c_int), intent(out),   optional :: info
     end subroutine pastix_task_analyze_f08
  end interface pastix_task_analyze

  interface pastix_task_numfact
     subroutine pastix_task_numfact_f08(pastix_data, spm, info)
       use :: iso_c_binding, only : c_int
       use :: pastixf_enums, only : pastix_data_t
       use :: spmf_enums,    only : spmatrix_t
       implicit none
       type(pastix_data_t), intent(inout), target   :: pastix_data
       type(spmatrix_t),    intent(inout), target   :: spm
       integer(kind=c_int), intent(out),   optional :: info
     end subroutine pastix_task_numfact_f08
  end interface pastix_task_numfact

  interface pastix_task_solve
     subroutine pastix_task_solve_f08(pastix_data, nrhs, B, ldb, info)
       use :: iso_c_binding,    only : c_int, c_ptr
       use :: pastixf_bindings, only : pastixGetCptrFrom2dArray
       use :: pastixf_enums,    only : pastix_data_t, pastix_int_t
       implicit none
       type(pastix_data_t),        intent(inout), target   :: pastix_data
       integer(kind=pastix_int_t), intent(in)              :: nrhs
       class(*),                   intent(inout), target   :: B(:,:)
       integer(kind=pastix_int_t), intent(in)              :: ldb
       integer(kind=c_int),        intent(out),   optional :: info
     end subroutine pastix_task_solve_f08
  end interface pastix_task_solve

  interface pastix_task_refine
     subroutine pastix_task_refine_f08(pastix_data, n, nrhs, B, ldb, X, ldx, &
          info)
       use :: iso_c_binding,    only : c_int, c_ptr
       use :: pastixf_bindings, only : pastixGetCptrFrom2dArray
       use :: pastixf_enums,    only : pastix_data_t, pastix_int_t
       implicit none
       type(pastix_data_t),        intent(inout), target   :: pastix_data
       integer(kind=pastix_int_t), intent(in)              :: n
       integer(kind=pastix_int_t), intent(in)              :: nrhs
       class(*),                   intent(inout), target   :: B(:,:)
       integer(kind=pastix_int_t), intent(in)              :: ldb
       class(*),                   intent(inout), target   :: X(:,:)
       integer(kind=pastix_int_t), intent(in)              :: ldx
       integer(kind=c_int),        intent(out),   optional :: info
     end subroutine pastix_task_refine_f08
  end interface pastix_task_refine

  interface pastix_subtask_order
     subroutine pastix_subtask_order_f08(pastix_data, spm, myorder, info)
       use :: iso_c_binding, only : c_int
       use :: pastixf_enums, only : pastix_data_t, pastix_order_t
       use :: spmf_enums,    only : spmatrix_t
       implicit none
       type(pastix_data_t),  intent(inout), target   :: pastix_data
       type(spmatrix_t),     intent(in),    target   :: spm
       type(pastix_order_t), intent(inout), target   :: myorder
       integer(kind=c_int),  intent(out),   optional :: info
     end subroutine pastix_subtask_order_f08
  end interface pastix_subtask_order

  interface pastix_subtask_symbfact
     subroutine pastix_subtask_symbfact_f08(pastix_data, info)
       use :: iso_c_binding, only : c_int
       use :: pastixf_enums, only : pastix_data_t
       implicit none
       type(pastix_data_t), intent(inout), target   :: pastix_data
       integer(kind=c_int), intent(out),   optional :: info
     end subroutine pastix_subtask_symbfact_f08
  end interface pastix_subtask_symbfact

  interface pastix_subtask_reordering
     subroutine pastix_subtask_reordering_f08(pastix_data, info)
       use :: iso_c_binding, only : c_int
       use :: pastixf_enums, only : pastix_data_t
       implicit none
       type(pastix_data_t), intent(inout), target   :: pastix_data
       integer(kind=c_int), intent(out),   optional :: info
     end subroutine pastix_subtask_reordering_f08
  end interface pastix_subtask_reordering

  interface pastix_subtask_blend
     subroutine pastix_subtask_blend_f08(pastix_data, info)
       use :: iso_c_binding, only : c_int
       use :: pastixf_enums, only : pastix_data_t
       implicit none
       type(pastix_data_t), intent(inout), target   :: pastix_data
       integer(kind=c_int), intent(out),   optional :: info
     end subroutine pastix_subtask_blend_f08
  end interface pastix_subtask_blend

  interface pastix_subtask_spm2bcsc
     subroutine pastix_subtask_spm2bcsc_f08(pastix_data, spm, info)
       use :: iso_c_binding, only : c_int
       use :: pastixf_enums, only : pastix_data_t
       use :: spmf_enums,    only : spmatrix_t
       implicit none
       type(pastix_data_t), intent(inout), target   :: pastix_data
       type(spmatrix_t),    intent(inout), target   :: spm
       integer(kind=c_int), intent(out),   optional :: info
     end subroutine pastix_subtask_spm2bcsc_f08
  end interface pastix_subtask_spm2bcsc

  interface pastix_subtask_bcsc2ctab
     subroutine pastix_subtask_bcsc2ctab_f08(pastix_data, info)
       use :: iso_c_binding, only : c_int
       use :: pastixf_enums, only : pastix_data_t
       implicit none
       type(pastix_data_t), intent(inout), target   :: pastix_data
       integer(kind=c_int), intent(out),   optional :: info
     end subroutine pastix_subtask_bcsc2ctab_f08
  end interface pastix_subtask_bcsc2ctab

  interface pastix_subtask_sopalin
     subroutine pastix_subtask_sopalin_f08(pastix_data, info)
       use :: iso_c_binding, only : c_int
       use :: pastixf_enums, only : pastix_data_t
       implicit none
       type(pastix_data_t), intent(inout), target   :: pastix_data
       integer(kind=c_int), intent(out),   optional :: info
     end subroutine pastix_subtask_sopalin_f08
  end interface pastix_subtask_sopalin

  interface pastix_subtask_applyorder
     subroutine pastix_subtask_applyorder_f08(pastix_data, flttype, dir, m, n, &
          B, ldb, info)
       use :: iso_c_binding,    only : c_int, c_ptr
       use :: pastixf_bindings, only : pastixGetCptrFrom2dArray
       use :: pastixf_enums,    only : pastix_data_t, pastix_int_t
       implicit none
       type(pastix_data_t),        intent(inout), target   :: pastix_data
       integer(c_int),             intent(in)              :: flttype
       integer(c_int),             intent(in)              :: dir
       integer(kind=pastix_int_t), intent(in)              :: m
       integer(kind=pastix_int_t), intent(in)              :: n
       class(*),                   intent(inout), target   :: B(:,:)
       integer(kind=pastix_int_t), intent(in)              :: ldb
       integer(kind=c_int),        intent(out),   optional :: info
     end subroutine pastix_subtask_applyorder_f08
  end interface pastix_subtask_applyorder

  interface pastix_subtask_trsm
     subroutine pastix_subtask_trsm_f08(pastix_data, flttype, side, uplo, &
          trans, diag, nrhs, B, ldb, info)
       use :: iso_c_binding,    only : c_int, c_ptr
       use :: pastixf_bindings, only : pastixGetCptrFrom2dArray
       use :: pastixf_enums,    only : pastix_data_t, pastix_int_t
       implicit none
       type(pastix_data_t),        intent(inout), target   :: pastix_data
       integer(c_int),             intent(in)              :: flttype
       integer(c_int),             intent(in)              :: side
       integer(c_int),             intent(in)              :: uplo
       integer(c_int),             intent(in)              :: trans
       integer(c_int),             intent(in)              :: diag
       integer(kind=pastix_int_t), intent(in)              :: nrhs
       class(*),                   intent(inout), target   :: B(:,:)
       integer(kind=pastix_int_t), intent(in)              :: ldb
       integer(kind=c_int),        intent(out),   optional :: info
     end subroutine pastix_subtask_trsm_f08
  end interface pastix_subtask_trsm

  interface pastix_subtask_diag
     subroutine pastix_subtask_diag_f08(pastix_data, flttype, nrhs, B, ldb, &
          info)
       use :: iso_c_binding,    only : c_int, c_ptr
       use :: pastixf_bindings, only : pastixGetCptrFrom2dArray
       use :: pastixf_enums,    only : pastix_data_t, pastix_int_t
       implicit none
       type(pastix_data_t),        intent(inout), target   :: pastix_data
       integer(c_int),             intent(in)              :: flttype
       integer(kind=pastix_int_t), intent(in)              :: nrhs
       class(*),                   intent(inout), target   :: B(:,:)
       integer(kind=pastix_int_t), intent(in)              :: ldb
       integer(kind=c_int),        intent(out),   optional :: info
     end subroutine pastix_subtask_diag_f08
  end interface pastix_subtask_diag

  interface pastix_subtask_solve
     subroutine pastix_subtask_solve_f08(pastix_data, nrhs, B, ldb, info)
       use :: iso_c_binding,    only : c_int, c_ptr
       use :: pastixf_bindings, only : pastixGetCptrFrom2dArray
       use :: pastixf_enums,    only : pastix_data_t, pastix_int_t
       implicit none
       type(pastix_data_t),        intent(inout), target   :: pastix_data
       integer(kind=pastix_int_t), intent(in)              :: nrhs
       class(*),                   intent(inout), target   :: B(:,:)
       integer(kind=pastix_int_t), intent(in)              :: ldb
       integer(kind=c_int),        intent(out),   optional :: info
     end subroutine pastix_subtask_solve_f08
  end interface pastix_subtask_solve

  interface pastix_subtask_refine
     subroutine pastix_subtask_refine_f08(pastix_data, n, nrhs, B, ldb, X, &
          ldx, info)
       use :: iso_c_binding,    only : c_int, c_ptr
       use :: pastixf_bindings, only : pastixGetCptrFrom2dArray
       use :: pastixf_enums,    only : pastix_data_t, pastix_int_t
       implicit none
       type(pastix_data_t),        intent(inout), target   :: pastix_data
       integer(kind=pastix_int_t), intent(in)              :: n
       integer(kind=pastix_int_t), intent(in)              :: nrhs
       class(*),                   intent(in),    target   :: B(:,:)
       integer(kind=pastix_int_t), intent(in)              :: ldb
       class(*),                   intent(inout), target   :: X(:,:)
       integer(kind=pastix_int_t), intent(in)              :: ldx
       integer(kind=c_int),        intent(out),   optional :: info
     end subroutine pastix_subtask_refine_f08
  end interface pastix_subtask_refine

  interface pastix_subtask_solve_adv
     subroutine pastix_subtask_solve_adv_f08(pastix_data, transA, nrhs, B, &
          ldb, info)
       use :: iso_c_binding,    only : c_int, c_ptr
       use :: pastixf_bindings, only : pastixGetCptrFrom2dArray
       use :: pastixf_enums,    only : pastix_data_t, pastix_int_t
       implicit none
       type(pastix_data_t),        intent(inout), target   :: pastix_data
       integer(c_int),             intent(in)              :: transA
       integer(kind=pastix_int_t), intent(in)              :: nrhs
       class(*),                   intent(inout), target   :: B(:,:)
       integer(kind=pastix_int_t), intent(in)              :: ldb
       integer(kind=c_int),        intent(out),   optional :: info
     end subroutine pastix_subtask_solve_adv_f08
  end interface pastix_subtask_solve_adv

  interface pastixSetSchurUnknownList
     subroutine pastixSetSchurUnknownList_f08(pastix_data, n, list)
       use :: pastixf_enums, only : pastix_data_t, pastix_int_t
       implicit none
       type(pastix_data_t),        intent(inout), target :: pastix_data
       integer(kind=pastix_int_t), intent(in)            :: n
       integer(kind=pastix_int_t), intent(in),    target :: list
     end subroutine pastixSetSchurUnknownList_f08
  end interface pastixSetSchurUnknownList

  interface pastixGetSchur
     subroutine pastixGetSchur_f08(pastix_data, S, lds, info)
       use :: iso_c_binding,    only : c_int, c_ptr
       use :: pastixf_bindings, only : pastixGetCptrFrom2dArray
       use :: pastixf_enums,    only : pastix_data_t, pastix_int_t
       implicit none
       type(pastix_data_t),        intent(in),    target   :: pastix_data
       class(*),                   intent(inout), target   :: S(:,:)
       integer(kind=pastix_int_t), intent(in)              :: lds
       integer(kind=c_int),        intent(out),   optional :: info
     end subroutine pastixGetSchur_f08
  end interface pastixGetSchur

  interface pastixExpand
     subroutine pastixExpand_f08(pastix_data, spm)
       use :: pastixf_enums, only : pastix_data_t
       use :: spmf_enums,    only : spmatrix_t
       implicit none
       type(pastix_data_t), intent(in),    target :: pastix_data
       type(spmatrix_t),    intent(inout), target :: spm
     end subroutine pastixExpand_f08
  end interface pastixExpand

  interface pastixGetDiag
     subroutine pastixGetDiag_f08(pastix_data, x, incx, info)
       use :: iso_c_binding,    only : c_int, c_ptr
       use :: pastixf_bindings, only : pastixGetCptrFrom1dArray
       use :: pastixf_enums,    only : pastix_data_t, pastix_int_t
       implicit none
       type(pastix_data_t),        intent(in),    target   :: pastix_data
       class(*),                   intent(inout), target   :: x(:)
       integer(kind=pastix_int_t), intent(in)              :: incx
       integer(kind=c_int),        intent(out),   optional :: info
     end subroutine pastixGetDiag_f08
  end interface pastixGetDiag

  interface pastixGetOptions
     subroutine pastixGetOptions_f08(argc, argv, iparm, dparm, check, driver, &
          filename)
       use :: iso_c_binding, only : c_char, c_double, c_int, c_ptr
       use :: pastixf_enums, only : pastix_int_t
       implicit none
       integer(kind=c_int),        intent(in)             :: argc
       character(kind=c_char),     intent(inout), pointer :: argv
       integer(kind=pastix_int_t), intent(inout), target  :: iparm(:)
       real(kind=c_double),        intent(inout), target  :: dparm(:)
       integer(kind=c_int),        intent(inout), target  :: check
       integer(c_int),             intent(inout), target  :: driver
       character(kind=c_char),     intent(inout), pointer :: filename
     end subroutine pastixGetOptions_f08
  end interface pastixGetOptions

  interface pastixDumpParam
     subroutine pastixDumpParam_f08(pastix_data)
       use :: pastixf_enums, only : pastix_data_t
       implicit none
       type(pastix_data_t), intent(in), target :: pastix_data
     end subroutine pastixDumpParam_f08
  end interface pastixDumpParam

  interface pastixCheckParam
     subroutine pastixCheckParam_f08(iparm, dparm, info)
       use :: iso_c_binding, only : c_double, c_int
       use :: pastixf_enums, only : pastix_int_t
       implicit none
       integer(kind=pastix_int_t), intent(in),  target   :: iparm(:)
       real(kind=c_double),        intent(in),  target   :: dparm(:)
       integer(kind=c_int),        intent(out), optional :: info
     end subroutine pastixCheckParam_f08
  end interface pastixCheckParam


  interface pastixOrderGetArray
     subroutine pastixOrderGetArray_f08( order, permtab, peritab, rangtab, treetab, sndetab )
       use :: pastixf_enums, only : pastix_order_t, pastix_int_t
       implicit none
       type(pastix_order_t),                intent(in),            target  :: order
       integer(pastix_int_t), dimension(:), intent(out), optional, pointer :: permtab
       integer(pastix_int_t), dimension(:), intent(out), optional, pointer :: peritab
       integer(pastix_int_t), dimension(:), intent(out), optional, pointer :: rangtab
       integer(pastix_int_t), dimension(:), intent(out), optional, pointer :: treetab
       integer(pastix_int_t), dimension(:), intent(out), optional, pointer :: sndetab
     end subroutine pastixOrderGetArray_f08
  end interface pastixOrderGetArray

end module pastixf_interfaces
