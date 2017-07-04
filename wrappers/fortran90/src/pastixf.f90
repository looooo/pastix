module pastixf
  use iso_c_binding
  use pastix_enums
  use spmf
  implicit none

  ! C structs converted to derived types.
  type, bind(c) :: pastix_data_t
     type(c_ptr) :: ptr
  end type pastix_data_t

  type, bind(c) :: Order
     type(c_ptr) :: ptr
  end type Order

  ! Interfaces of the C functions.
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
     function pastix_subtask_bcsc2ctab_c(pastix_data, spm) &
          bind(c, name='pastix_subtask_bcsc2ctab')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: pastix_subtask_bcsc2ctab_c
       type(c_ptr), value :: pastix_data
       type(c_ptr), value :: spm
     end function pastix_subtask_bcsc2ctab_c
  end interface

  interface
     function pastix_subtask_sopalin_c(pastix_data, spm) &
          bind(c, name='pastix_subtask_sopalin')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: pastix_subtask_sopalin_c
       type(c_ptr), value :: pastix_data
       type(c_ptr), value :: spm
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
     subroutine pastix_setSchurUnknownList_c(pastix_data, n, list) &
          bind(c, name='pastix_setSchurUnknownList')
       use iso_c_binding
       import pastix_int_t
       implicit none
       type(c_ptr), value :: pastix_data
       integer(kind=pastix_int_t), value :: n
       type(c_ptr), value :: list
     end subroutine pastix_setSchurUnknownList_c
  end interface

  interface
     function pastix_getSchur_c(pastix_data, S, lds) &
          bind(c, name='pastix_getSchur')
       use iso_c_binding
       import pastix_int_t
       implicit none
       integer(kind=c_int) :: pastix_getSchur_c
       type(c_ptr), value :: pastix_data
       type(c_ptr), value :: S
       integer(kind=pastix_int_t), value :: lds
     end function pastix_getSchur_c
  end interface

  interface
     subroutine pastix_getOptions_c(argc, argv, iparam, dparam, check, driver, filename) &
          bind(c, name='pastix_getOptions')
       use iso_c_binding
       implicit none
       integer(kind=c_int), value :: argc
       type(c_ptr) :: argv
       type(c_ptr), value :: iparam
       type(c_ptr), value :: dparam
       type(c_ptr), value :: check
       type(c_ptr), value :: driver
       type(c_ptr) :: filename
     end subroutine pastix_getOptions_c
  end interface

contains

  ! Wrappers of the C functions.
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

    info = pastix_c(pastix_data_aux, pastix_comm, n, c_loc(colptr), c_loc(row), c_loc(avals), c_loc(perm), c_loc(invp), &
         c_loc(b), nrhs, c_loc(iparm), c_loc(dparm))
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

  subroutine pastixFinalize(pastix_data)
    use iso_c_binding
    implicit none
    type(pastix_data_t), intent(inout), pointer :: pastix_data

    type(c_ptr) :: pastix_data_aux

    call pastixFinalize_c(pastix_data_aux)
    call c_f_pointer(pastix_data_aux, pastix_data)
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
    type(pastix_data_t), intent(inout), target :: pastix_data
    type(pastix_spm_t),  intent(in),    target :: spm
    type(Order),         intent(inout), target :: myorder
    integer(kind=c_int), intent(out)           :: info

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
    type(pastix_data_t), intent(inout), target :: pastix_data
    type(pastix_spm_t),  intent(inout), target :: spm
    integer(kind=c_int), intent(out)           :: info

    info = pastix_subtask_spm2bcsc_c(c_loc(pastix_data), c_loc(spm))
  end subroutine pastix_subtask_spm2bcsc

  subroutine pastix_subtask_bcsc2ctab(pastix_data, spm, info)
    use iso_c_binding
    implicit none
    type(pastix_data_t), intent(inout), target :: pastix_data
    type(pastix_spm_t),  intent(inout), target :: spm
    integer(kind=c_int), intent(out)           :: info

    info = pastix_subtask_bcsc2ctab_c(c_loc(pastix_data), c_loc(spm))
  end subroutine pastix_subtask_bcsc2ctab

  subroutine pastix_subtask_sopalin(pastix_data, spm, info)
    use iso_c_binding
    implicit none
    type(pastix_data_t), intent(inout), target :: pastix_data
    type(pastix_spm_t),  intent(inout), target :: spm
    integer(kind=c_int), intent(out)           :: info

    info = pastix_subtask_sopalin_c(c_loc(pastix_data), c_loc(spm))
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

    info = pastix_subtask_applyorder_c(c_loc(pastix_data), flttype, dir, m, n, c_loc(b), ldb)
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
         nrhs, c_loc(b), ldb)
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

    info = pastix_subtask_diag_c(c_loc(pastix_data), flttype, nrhs, c_loc(b), ldb)
  end subroutine pastix_subtask_diag

  subroutine pastix_setSchurUnknownList(pastix_data, n, list)
    use iso_c_binding
    implicit none
    type(pastix_data_t),        intent(inout), target :: pastix_data
    integer(kind=pastix_int_t), intent(in)            :: n
    integer(kind=pastix_int_t), intent(in),    target :: list

    call pastix_setSchurUnknownList_c(c_loc(pastix_data), n, c_loc(list))
  end subroutine pastix_setSchurUnknownList

  subroutine pastix_getSchur(pastix_data, S, lds, info)
    use iso_c_binding
    implicit none
    type(pastix_data_t),        intent(in),    target :: pastix_data
    type(c_ptr),                intent(inout), target :: S
    integer(kind=pastix_int_t), intent(in)            :: lds
    integer(kind=c_int),        intent(out)           :: info

    info = pastix_getSchur_c(c_loc(pastix_data), c_loc(S), lds)
  end subroutine pastix_getSchur

  subroutine pastix_getOptions(argc, argv, iparm, dparm, check, driver, filename)
    use iso_c_binding
    implicit none
    integer(kind=c_int),        intent(in)                                   :: argc
    character(kind=c_char),     intent(inout), pointer                       :: argv
    integer(kind=pastix_int_t), intent(inout), dimension(iparm_size), target :: iparm
    real(kind=c_double),        intent(inout), dimension(dparm_size), target :: dparm
    integer(kind=c_int),        intent(inout), target                        :: check
    integer(c_int),             intent(inout), target                        :: driver
    character(kind=c_char),     intent(in),    pointer                       :: filename

    type(c_ptr) :: argv_aux
    type(c_ptr) :: filename_aux

    call pastix_getOptions_c(argc, argv_aux, c_loc(iparam), c_loc(dparam), c_loc(check), c_loc(driver), filename_aux)
    call c_f_pointer(argv_aux, argv)
    call c_f_pointer(filename_aux, filename)
  end subroutine pastix_getOptions

end module pastixf
