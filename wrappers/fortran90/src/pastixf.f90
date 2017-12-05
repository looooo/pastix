!
! @file pastixf.f90
!
! PaStiX routine wrappers for Fortan 90
!
! @copyright 2017-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
!                      Univ. Bordeaux. All rights reserved.
!
! @version 6.0.0
! @author Mathieu Faverge
! @date 2017-01-01
!
module pastixf
  use iso_c_binding
  use pastix_enums
  use spmf
  implicit none

  ! C structs converted to derived types.
  type, bind(c) :: pastix_data_t
     type(c_ptr) :: ptr
  end type pastix_data_t

  type, bind(c) :: pastix_order_t
    integer(kind=pastix_int_t) :: baseval
    integer(kind=pastix_int_t) :: vertnbr
    integer(kind=pastix_int_t) :: cblknbr
    type(c_ptr)                :: permtab
    type(c_ptr)                :: peritab
    type(c_ptr)                :: rangtab
    type(c_ptr)                :: treetab
  end type pastix_order_t

  ! Interfaces of the C functions of the Order module.
  interface
     function pastixOrderInit_c( order, baseval, vertnbr, cblknbr, permtab, peritab, rangtab, treetab ) &
          bind(c, name='pastixOrderInit')
       use iso_c_binding
       import pastix_int_t
       implicit none
       integer(kind=c_int) :: pastixOrderInit_c
       type(c_ptr),                value :: order
       integer(kind=pastix_int_t), value :: baseval
       integer(kind=pastix_int_t), value :: vertnbr
       integer(kind=pastix_int_t), value :: cblknbr
       type(c_ptr),                value :: permtab
       type(c_ptr),                value :: peritab
       type(c_ptr),                value :: rangtab
       type(c_ptr),                value :: treetab
     end function pastixOrderInit_c
  end interface

  interface
     function pastixOrderAlloc_c( order, vertnbr, cblknbr ) &
          bind(c, name='pastixOrderAlloc')
       use iso_c_binding
       import pastix_int_t
       implicit none
       integer(kind=c_int) :: pastixOrderAlloc_c
       type(c_ptr),                value :: order
       integer(kind=pastix_int_t), value :: vertnbr
       integer(kind=pastix_int_t), value :: cblknbr
     end function pastixOrderAlloc_c
  end interface

  interface
     subroutine pastixOrderExit_c( order ) &
          bind(c, name='pastixOrderExit')
       use iso_c_binding
       import pastix_int_t
       implicit none
       type(c_ptr), value :: order
     end subroutine pastixOrderExit_c
  end interface

  interface
     function pastixOrderGet_c( pastix_data ) &
          bind(c, name='pastixOrderGet')
       use iso_c_binding
       import pastix_int_t
       implicit none
       type(c_ptr)        :: pastixOrderGet_c
       type(c_ptr), value :: pastix_data
     end function pastixOrderGet_c
  end interface

  ! Interfaces of the C functions of the PaStiX module.
  interface
     function pastix_c(pastix_data, pastix_comm, n, colptr, row, avals, perm, invp, &
          b, nrhs, iparm, dparm) &
          bind(c, name='pastix')
       use iso_c_binding
       import pastix_int_t
       implicit none
       integer(kind=c_int) :: pastix_c
       type(c_ptr) :: pastix_data
       integer(kind=c_int), value :: pastix_comm
       integer(kind=pastix_int_t), value :: n
       type(c_ptr), value :: colptr
       type(c_ptr), value :: row
       type(c_ptr), value :: avals
       type(c_ptr), value :: perm
       type(c_ptr), value :: invp
       type(c_ptr), value :: b
       integer(kind=pastix_int_t), value :: nrhs
       type(c_ptr), value :: iparm
       type(c_ptr), value :: dparm
     end function pastix_c
  end interface

  interface
     subroutine pastixInitParam_c(iparm, dparm) &
          bind(c, name='pastixInitParam')
       use iso_c_binding
       implicit none
       type(c_ptr), value :: iparm
       type(c_ptr), value :: dparm
     end subroutine pastixInitParam_c
  end interface

  interface
     subroutine pastixInit_c(pastix_data, pastix_comm, iparm, dparm) &
          bind(c, name='pastixInit')
       use iso_c_binding
       implicit none
       type(c_ptr) :: pastix_data
       integer(kind=c_int), value :: pastix_comm
       type(c_ptr), value :: iparm
       type(c_ptr), value :: dparm
     end subroutine pastixInit_c
  end interface

  interface
     subroutine pastixInitWithAffinity_c(pastix_data, pastix_comm, iparm, dparm, bindtab) &
          bind(c, name='pastixInitWithAffinity')
       use iso_c_binding
       implicit none
       type(c_ptr) :: pastix_data
       integer(kind=c_int), value :: pastix_comm
       type(c_ptr), value :: iparm
       type(c_ptr), value :: dparm
       type(c_ptr), value :: bindtab
     end subroutine pastixInitWithAffinity_c
  end interface

  interface
     subroutine pastixFinalize_c(pastix_data) &
          bind(c, name='pastixFinalize')
       use iso_c_binding
       implicit none
       type(c_ptr) :: pastix_data
     end subroutine pastixFinalize_c
  end interface

  interface
     function pastix_task_analyze_c(pastix_data, spm) &
          bind(c, name='pastix_task_analyze')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: pastix_task_analyze_c
       type(c_ptr), value :: pastix_data
       type(c_ptr), value :: spm
     end function pastix_task_analyze_c
  end interface

  interface
     function pastix_task_numfact_c(pastix_data, spm) &
          bind(c, name='pastix_task_numfact')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: pastix_task_numfact_c
       type(c_ptr), value  :: pastix_data
       type(c_ptr), value  :: spm
     end function pastix_task_numfact_c
  end interface

  interface
     function pastix_task_solve_c(pastix_data, nrhs, b, ldb) &
          bind(c, name='pastix_task_solve')
       use iso_c_binding
       import pastix_int_t
       implicit none
       integer(kind=c_int)               :: pastix_task_solve_c
       type(c_ptr),                value :: pastix_data
       integer(kind=pastix_int_t), value :: nrhs
       type(c_ptr),                value :: b
       integer(kind=pastix_int_t), value :: ldb
     end function pastix_task_solve_c
  end interface

  interface
     function pastix_task_refine_c(pastix_data, x, nrhs, b) &
          bind(c, name='pastix_task_refine')
       use iso_c_binding
       import pastix_int_t
       implicit none
       integer(kind=c_int)               :: pastix_task_refine_c
       type(c_ptr),                value :: pastix_data
       type(c_ptr),                value :: x
       integer(kind=pastix_int_t), value :: nrhs
       type(c_ptr),                value :: b
     end function pastix_task_refine_c
  end interface

  interface
     function pastix_subtask_order_c(pastix_data, spm, myorder) &
          bind(c, name='pastix_subtask_order')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: pastix_subtask_order_c
       type(c_ptr), value :: pastix_data
       type(c_ptr), value :: spm
       type(c_ptr), value :: myorder
     end function pastix_subtask_order_c
  end interface

  interface
     function pastix_subtask_symbfact_c(pastix_data) &
          bind(c, name='pastix_subtask_symbfact')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: pastix_subtask_symbfact_c
       type(c_ptr), value :: pastix_data
     end function pastix_subtask_symbfact_c
  end interface

  interface
     function pastix_subtask_reordering_c(pastix_data) &
          bind(c, name='pastix_subtask_reordering')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: pastix_subtask_reordering_c
       type(c_ptr), value :: pastix_data
     end function pastix_subtask_reordering_c
  end interface

  interface
     function pastix_subtask_blend_c(pastix_data) &
          bind(c, name='pastix_subtask_blend')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: pastix_subtask_blend_c
       type(c_ptr), value :: pastix_data
     end function pastix_subtask_blend_c
  end interface

  interface
     function pastix_subtask_spm2bcsc_c(pastix_data, spm) &
          bind(c, name='pastix_subtask_spm2bcsc')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: pastix_subtask_spm2bcsc_c
       type(c_ptr), value :: pastix_data
       type(c_ptr), value :: spm
     end function pastix_subtask_spm2bcsc_c
  end interface

  interface
     function pastix_subtask_bcsc2ctab_c(pastix_data) &
          bind(c, name='pastix_subtask_bcsc2ctab')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: pastix_subtask_bcsc2ctab_c
       type(c_ptr), value :: pastix_data
     end function pastix_subtask_bcsc2ctab_c
  end interface

  interface
     function pastix_subtask_sopalin_c(pastix_data) &
          bind(c, name='pastix_subtask_sopalin')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: pastix_subtask_sopalin_c
       type(c_ptr), value :: pastix_data
     end function pastix_subtask_sopalin_c
  end interface

  interface
     function pastix_subtask_applyorder_c(pastix_data, flttype, dir, m, n, b, ldb) &
          bind(c, name='pastix_subtask_applyorder')
       use iso_c_binding
       import pastix_int_t
       implicit none
       integer(kind=c_int) :: pastix_subtask_applyorder_c
       type(c_ptr), value :: pastix_data
       integer(c_int), value :: flttype
       integer(c_int), value :: dir
       integer(kind=pastix_int_t), value :: m
       integer(kind=pastix_int_t), value :: n
       type(c_ptr), value :: b
       integer(kind=pastix_int_t), value :: ldb
     end function pastix_subtask_applyorder_c
  end interface

  interface
     function pastix_subtask_trsm_c(pastix_data, flttype, side, uplo, trans, diag, nrhs, b, &
          ldb) &
          bind(c, name='pastix_subtask_trsm')
       use iso_c_binding
       import pastix_int_t
       implicit none
       integer(kind=c_int) :: pastix_subtask_trsm_c
       type(c_ptr), value :: pastix_data
       integer(c_int), value :: flttype
       integer(c_int), value :: side
       integer(c_int), value :: uplo
       integer(c_int), value :: trans
       integer(c_int), value :: diag
       integer(kind=pastix_int_t), value :: nrhs
       type(c_ptr), value :: b
       integer(kind=pastix_int_t), value :: ldb
     end function pastix_subtask_trsm_c
  end interface

  interface
     function pastix_subtask_diag_c(pastix_data, flttype, nrhs, b, ldb) &
          bind(c, name='pastix_subtask_diag')
       use iso_c_binding
       import pastix_int_t
       implicit none
       integer(kind=c_int) :: pastix_subtask_diag_c
       type(c_ptr), value :: pastix_data
       integer(c_int), value :: flttype
       integer(kind=pastix_int_t), value :: nrhs
       type(c_ptr), value :: b
       integer(kind=pastix_int_t), value :: ldb
     end function pastix_subtask_diag_c
  end interface

  interface
     subroutine pastixSetSchurUnknownList_c(pastix_data, n, list) &
          bind(c, name='pastixSetSchurUnknownList')
       use iso_c_binding
       import pastix_int_t
       implicit none
       type(c_ptr), value :: pastix_data
       integer(kind=pastix_int_t), value :: n
       type(c_ptr), value :: list
     end subroutine pastixSetSchurUnknownList_c
  end interface

  interface
     function pastixGetSchur_c(pastix_data, S, lds) &
          bind(c, name='pastixGetSchur')
       use iso_c_binding
       import pastix_int_t
       implicit none
       integer(kind=c_int) :: pastixGetSchur_c
       type(c_ptr), value :: pastix_data
       type(c_ptr), value :: S
       integer(kind=pastix_int_t), value :: lds
     end function pastixGetSchur_c
  end interface

  interface
     subroutine pastixGetOptions_c(argc, argv, iparam, dparam, check, driver, filename) &
          bind(c, name='pastixGetOptions')
       use iso_c_binding
       implicit none
       integer(kind=c_int), value :: argc
       type(c_ptr) :: argv
       type(c_ptr), value :: iparam
       type(c_ptr), value :: dparam
       type(c_ptr), value :: check
       type(c_ptr), value :: driver
       type(c_ptr) :: filename
     end subroutine pastixGetOptions_c
  end interface

contains


  ! Wrappers of the C functions of the order module.
  subroutine pastixOrderInit( order, baseval, vertnbr, cblknbr, permtab, peritab, rangtab, treetab, info )
    use iso_c_binding
    implicit none
    type(pastix_order_t),       intent(in), target               :: order
    integer(kind=pastix_int_t), intent(in)                       :: baseval
    integer(kind=pastix_int_t), intent(in)                       :: vertnbr
    integer(kind=pastix_int_t), intent(in)                       :: cblknbr
    integer(kind=pastix_int_t), dimension(:), intent(in), target :: permtab
    integer(kind=pastix_int_t), dimension(:), intent(in), target :: peritab
    integer(kind=pastix_int_t), dimension(:), intent(in), target :: rangtab
    integer(kind=pastix_int_t), dimension(:), intent(in), target :: treetab
    integer(kind=c_int),        intent(out)                      :: info

    info = pastixOrderInit_c( c_loc(order), baseval, vertnbr, cblknbr, c_loc(permtab), c_loc(peritab), &
         c_loc(rangtab), c_loc(treetab) )
  end subroutine pastixOrderInit

  subroutine pastixOrderAlloc( order, vertnbr, cblknbr, info )
    use iso_c_binding
    implicit none
    type(pastix_order_t),       intent(in), target :: order
    integer(kind=pastix_int_t), intent(in)         :: vertnbr
    integer(kind=pastix_int_t), intent(in)         :: cblknbr
    integer(kind=c_int),        intent(out)        :: info

    info = pastixOrderAlloc_c( c_loc(order), vertnbr, cblknbr )
  end subroutine pastixOrderAlloc

  subroutine pastixOrderExit( order )
    use iso_c_binding
    implicit none
    type(pastix_order_t), intent(in), target :: order

    call pastixOrderExit_c( c_loc(order) )

  end subroutine pastixOrderExit

  subroutine pastixOrderGet( pastix_data, order )
    use iso_c_binding
    implicit none
    type(pastix_order_t), intent(inout), pointer :: order
    type(pastix_data_t),  intent(in),    target  :: pastix_data
    type(c_ptr) :: order_aux

    order_aux = pastixOrderGet_c( c_loc(pastix_data) )
    call c_f_pointer( order_aux, order )

  end subroutine pastixOrderGet

  ! Wrappers of the C functions of the PaStiX module.
  subroutine pastix(pastix_data, pastix_comm, n, colptr, row, avals, perm, invp, &
       b, nrhs, iparm, dparm, info)
    use iso_c_binding
    implicit none
    type(pastix_data_t),        intent(inout), pointer                       :: pastix_data
    integer(kind=c_int),        intent(in)                                   :: pastix_comm
    integer(kind=pastix_int_t), intent(in)                                   :: n
    integer(kind=pastix_int_t), intent(inout), target                        :: colptr
    integer(kind=pastix_int_t), intent(inout), target                        :: row
    type(c_ptr),                intent(inout), target                        :: avals
    integer(kind=pastix_int_t), intent(inout), target                        :: perm
    integer(kind=pastix_int_t), intent(inout), target                        :: invp
    type(c_ptr),                intent(inout), target                        :: b
    integer(kind=pastix_int_t), intent(in)                                   :: nrhs
    integer(kind=pastix_int_t), intent(inout), dimension(iparm_size), target :: iparm
    real(kind=c_double),        intent(inout), dimension(dparm_size), target :: dparm
    integer(kind=c_int),        intent(out)                                  :: info

    type(c_ptr) :: pastix_data_aux

    pastix_data_aux = c_loc(pastix_data)

    info = pastix_c(pastix_data_aux, pastix_comm, n, c_loc(colptr), c_loc(row), avals, c_loc(perm), c_loc(invp), &
         b, nrhs, c_loc(iparm), c_loc(dparm))
    call c_f_pointer(pastix_data_aux, pastix_data)
  end subroutine pastix

  subroutine pastixInitParam(iparm, dparm)
    use iso_c_binding
    implicit none
    integer(kind=pastix_int_t), intent(inout), dimension(iparm_size), target :: iparm
    real(kind=c_double),        intent(inout), dimension(dparm_size), target :: dparm

    call pastixInitParam_c(c_loc(iparm), c_loc(dparm))
  end subroutine pastixInitParam

  subroutine pastixInit(pastix_data, pastix_comm, iparm, dparm)
    use iso_c_binding
    implicit none
    type(pastix_data_t),        intent(inout), pointer                       :: pastix_data
    integer(kind=c_int),        intent(in)                                   :: pastix_comm
    integer(kind=pastix_int_t), intent(inout), dimension(iparm_size), target :: iparm
    real(kind=c_double),        intent(inout), dimension(dparm_size), target :: dparm

    type(c_ptr) :: pastix_data_aux

    call pastixInit_c(pastix_data_aux, pastix_comm, c_loc(iparm), c_loc(dparm))
    call c_f_pointer(pastix_data_aux, pastix_data)
  end subroutine pastixInit

  subroutine pastixInitWithAffinity(pastix_data, pastix_comm, iparm, dparm, bindtab)
    use iso_c_binding
    implicit none
    type(pastix_data_t),        intent(inout), pointer                       :: pastix_data
    integer(kind=c_int),        intent(in)                                   :: pastix_comm
    integer(kind=pastix_int_t), intent(inout), dimension(iparm_size), target :: iparm
    real(kind=c_double),        intent(inout), dimension(dparm_size), target :: dparm
    integer(kind=c_int),        intent(in),    dimension(:), target          :: bindtab

    type(c_ptr) :: pastix_data_aux

    call pastixInitWithAffinity_c(pastix_data_aux, pastix_comm, c_loc(iparm), c_loc(dparm), c_loc(bindtab))
    call c_f_pointer(pastix_data_aux, pastix_data)
  end subroutine pastixInitWithAffinity

  subroutine pastixFinalize(pastix_data)
    use iso_c_binding
    implicit none
    type(pastix_data_t), intent(inout), target :: pastix_data

    call pastixFinalize_c(c_loc(pastix_data))
  end subroutine pastixFinalize

  subroutine pastix_task_analyze(pastix_data, spm, info)
    use iso_c_binding
    implicit none
    type(pastix_data_t), intent(inout), target :: pastix_data
    type(pastix_spm_t),  intent(inout), target :: spm
    integer(kind=c_int), intent(out)           :: info

    info = pastix_task_analyze_c(c_loc(pastix_data), c_loc(spm))
  end subroutine pastix_task_analyze

  subroutine pastix_task_numfact(pastix_data, spm, info)
    use iso_c_binding
    implicit none
    type(pastix_data_t), intent(inout), target :: pastix_data
    type(pastix_spm_t),  intent(inout), target :: spm
    integer(kind=c_int), intent(out)           :: info

    info = pastix_task_numfact_c(c_loc(pastix_data), c_loc(spm))
  end subroutine pastix_task_numfact

  subroutine pastix_task_solve(pastix_data, nrhs, b, ldb, info)
    use iso_c_binding
    implicit none
    type(pastix_data_t),        intent(inout), target :: pastix_data
    integer(kind=pastix_int_t), intent(in)            :: nrhs
    type(c_ptr),                intent(inout)         :: b
    integer(kind=pastix_int_t), intent(in)            :: ldb
    integer(kind=c_int),        intent(out)           :: info

    info = pastix_task_solve_c(c_loc(pastix_data), nrhs, b, ldb)
  end subroutine pastix_task_solve

  subroutine pastix_task_refine(pastix_data, x, nrhs, b, info)
    use iso_c_binding
    implicit none
    type(pastix_data_t),        intent(inout), target :: pastix_data
    type(c_ptr),                intent(inout)         :: x
    integer(kind=pastix_int_t), intent(in)            :: nrhs
    type(c_ptr),                intent(inout)         :: b
    integer(kind=c_int),        intent(out)           :: info

    info = pastix_task_refine_c(c_loc(pastix_data), x, nrhs, b)
  end subroutine pastix_task_refine

  subroutine pastix_subtask_order(pastix_data, spm, myorder, info)
    use iso_c_binding
    implicit none
    type(pastix_data_t),  intent(inout), target  :: pastix_data
    type(pastix_spm_t),   intent(in),    target  :: spm
    type(pastix_order_t), intent(in),    pointer :: myorder
    integer(kind=c_int),  intent(out)            :: info

    info = pastix_subtask_order_c(c_loc(pastix_data), c_loc(spm), c_loc(myorder))
  end subroutine pastix_subtask_order

  subroutine pastix_subtask_symbfact(pastix_data, info)
    use iso_c_binding
    implicit none
    type(pastix_data_t), intent(inout), target :: pastix_data
    integer(kind=c_int), intent(out)           :: info

    info = pastix_subtask_symbfact_c(c_loc(pastix_data))
  end subroutine pastix_subtask_symbfact

  subroutine pastix_subtask_reordering(pastix_data, info)
    use iso_c_binding
    implicit none
    type(pastix_data_t), intent(inout), target :: pastix_data
    integer(kind=c_int), intent(out)           :: info

    info = pastix_subtask_reordering_c(c_loc(pastix_data))
  end subroutine pastix_subtask_reordering

  subroutine pastix_subtask_blend(pastix_data, info)
    use iso_c_binding
    implicit none
    type(pastix_data_t), intent(inout), target :: pastix_data
    integer(kind=c_int), intent(out)           :: info

    info = pastix_subtask_blend_c(c_loc(pastix_data))
  end subroutine pastix_subtask_blend

  subroutine pastix_subtask_spm2bcsc(pastix_data, spm, info)
    use iso_c_binding
    implicit none
    type(pastix_data_t), intent(in),    target :: pastix_data
    type(pastix_spm_t),  intent(inout), target :: spm
    integer(kind=c_int), intent(out)           :: info

    info = pastix_subtask_spm2bcsc_c(c_loc(pastix_data), c_loc(spm))
  end subroutine pastix_subtask_spm2bcsc

  subroutine pastix_subtask_bcsc2ctab(pastix_data, info)
    use iso_c_binding
    implicit none
    type(pastix_data_t), intent(in), target :: pastix_data
    integer(kind=c_int), intent(out)        :: info

    info = pastix_subtask_bcsc2ctab_c(c_loc(pastix_data))
  end subroutine pastix_subtask_bcsc2ctab

  subroutine pastix_subtask_sopalin(pastix_data, info)
    use iso_c_binding
    implicit none
    type(pastix_data_t), intent(in), target :: pastix_data
    integer(kind=c_int), intent(out)        :: info

    info = pastix_subtask_sopalin_c(c_loc(pastix_data))
  end subroutine pastix_subtask_sopalin

  subroutine pastix_subtask_applyorder(pastix_data, flttype, dir, m, n, b, ldb, info)
    use iso_c_binding
    implicit none
    type(pastix_data_t),        intent(inout), target :: pastix_data
    integer(c_int),             intent(in)            :: flttype
    integer(c_int),             intent(in)            :: dir
    integer(kind=pastix_int_t), intent(in)            :: m
    integer(kind=pastix_int_t), intent(in)            :: n
    type(c_ptr),                intent(inout), target :: b
    integer(kind=pastix_int_t), intent(in)            :: ldb
    integer(kind=c_int),        intent(out)           :: info

    info = pastix_subtask_applyorder_c(c_loc(pastix_data), flttype, dir, m, n, b, ldb)
  end subroutine pastix_subtask_applyorder

  subroutine pastix_subtask_trsm(pastix_data, flttype, side, uplo, trans, diag, nrhs, b, &
       ldb, info)
    use iso_c_binding
    implicit none
    type(pastix_data_t),        intent(inout), target :: pastix_data
    integer(c_int),             intent(in)            :: flttype
    integer(c_int),             intent(in)            :: side
    integer(c_int),             intent(in)            :: uplo
    integer(c_int),             intent(in)            :: trans
    integer(c_int),             intent(in)            :: diag
    integer(kind=pastix_int_t), intent(in)            :: nrhs
    type(c_ptr),                intent(inout), target :: b
    integer(kind=pastix_int_t), intent(in)            :: ldb
    integer(kind=c_int),        intent(out)           :: info

    info = pastix_subtask_trsm_c(c_loc(pastix_data), flttype, side, uplo, trans, diag, &
         nrhs, b, ldb)
  end subroutine pastix_subtask_trsm

  subroutine pastix_subtask_diag(pastix_data, flttype, nrhs, b, ldb, info)
    use iso_c_binding
    implicit none
    type(pastix_data_t),        intent(inout), target :: pastix_data
    integer(c_int),             intent(in)            :: flttype
    integer(kind=pastix_int_t), intent(in)            :: nrhs
    type(c_ptr),                intent(inout), target :: b
    integer(kind=pastix_int_t), intent(in)            :: ldb
    integer(kind=c_int),        intent(out)           :: info

    info = pastix_subtask_diag_c(c_loc(pastix_data), flttype, nrhs, b, ldb)
  end subroutine pastix_subtask_diag

  subroutine pastixSetSchurUnknownList(pastix_data, n, list)
    use iso_c_binding
    implicit none
    type(pastix_data_t),        intent(inout), target :: pastix_data
    integer(kind=pastix_int_t), intent(in)            :: n
    integer(kind=pastix_int_t), intent(in),    target :: list

    call pastixSetSchurUnknownList_c(c_loc(pastix_data), n, c_loc(list))
  end subroutine pastixSetSchurUnknownList

  subroutine pastixGetSchur(pastix_data, S, lds, info)
    use iso_c_binding
    implicit none
    type(pastix_data_t),        intent(in),    target :: pastix_data
    type(c_ptr),                intent(inout), target :: S
    integer(kind=pastix_int_t), intent(in)            :: lds
    integer(kind=c_int),        intent(out)           :: info

    info = pastixGetSchur_c(c_loc(pastix_data), S, lds)
  end subroutine pastixGetSchur

  subroutine pastixGetOptions(argc, argv, iparm, dparm, check, driver, filename)
    use iso_c_binding
    implicit none
    integer(kind=c_int),        intent(in)                                   :: argc
    character(kind=c_char),     intent(inout), pointer                       :: argv
    integer(kind=pastix_int_t), intent(inout), dimension(iparm_size), target :: iparm
    real(kind=c_double),        intent(inout), dimension(dparm_size), target :: dparm
    integer(kind=c_int),        intent(inout), target                        :: check
    integer(c_int),             intent(inout), target                        :: driver
    character(kind=c_char),     intent(inout), pointer                       :: filename

    type(c_ptr) :: argv_aux
    type(c_ptr) :: filename_aux

    call pastixGetOptions_c(argc, argv_aux, c_loc(iparm), c_loc(dparm), c_loc(check), c_loc(driver), filename_aux)
    call c_f_pointer(argv_aux, argv)
    call c_f_pointer(filename_aux, filename)
  end subroutine pastixGetOptions

end module pastixf
