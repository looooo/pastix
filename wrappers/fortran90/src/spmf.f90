
!
! @file spmf.f90
!
! SPM Fortran 90 wrapper
!
! @copyright 2017      Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
!                      Univ. Bordeaux. All rights reserved.
!
! @version 6.0.0
! @author Mathieu Faverge
! @date 2017-01-01
!
! This file has been automatically generated with gen_wrappers.py
!
module spmf
  use iso_c_binding
  use pastix_enums
  implicit none

  type, bind(c) :: pastix_spm_t
     integer(c_int)             :: mtxtype
     integer(c_int)             :: flttype
     integer(c_int)             :: fmttype
     integer(kind=pastix_int_t) :: gN
     integer(kind=pastix_int_t) :: n
     integer(kind=pastix_int_t) :: gnnz
     integer(kind=pastix_int_t) :: nnz
     integer(kind=pastix_int_t) :: gNexp
     integer(kind=pastix_int_t) :: nexp
     integer(kind=pastix_int_t) :: gnnzexp
     integer(kind=pastix_int_t) :: nnzexp
     integer(kind=pastix_int_t) :: dof
     type(c_ptr)                :: dofs
     integer(c_int)             :: layout
     type(c_ptr)                :: colptr
     type(c_ptr)                :: rowptr
     type(c_ptr)                :: loc2glob
     type(c_ptr)                :: values
  end type pastix_spm_t

  interface
     function spmNew_c(mtxtype, flttype, fmttype, n, nnz, colptr, rowptr, &
          values, loc2glob, dof, layout, dofs) &
          bind(c, name='spmNew')
       use iso_c_binding
       import pastix_int_t
       import pastix_spm_t
       implicit none
       type(c_ptr)                       :: spmNew_c
       integer(c_int),             value :: mtxtype
       integer(c_int),             value :: flttype
       integer(c_int),             value :: fmttype
       integer(kind=pastix_int_t), value :: n
       integer(kind=pastix_int_t), value :: nnz
       type(c_ptr),                value :: colptr
       type(c_ptr),                value :: rowptr
       type(c_ptr),                value :: values
       type(c_ptr),                value :: loc2glob
       integer(kind=pastix_int_t), value :: dof
       integer(c_int),             value :: layout
       type(c_ptr),                value :: dofs
     end function spmNew_c
  end interface

  interface
     subroutine spmInit_c(spm) &
          bind(c, name='spmInit')
       use iso_c_binding
       import pastix_spm_t
       implicit none
       type(c_ptr), value :: spm
     end subroutine spmInit_c
  end interface

  interface
     subroutine spmExit_c(spm) &
          bind(c, name='spmExit')
       use iso_c_binding
       import pastix_spm_t
       implicit none
       type(c_ptr), value :: spm
     end subroutine spmExit_c
  end interface

  interface
     subroutine spmFree_c(spm) &
          bind(c, name='spmFree')
       use iso_c_binding
       import pastix_spm_t
       implicit none
       type(c_ptr), value :: spm
     end subroutine spmFree_c
  end interface

  interface
     function spmCopy_c(spm) &
          bind(c, name='spmCopy')
       use iso_c_binding
       import pastix_spm_t
       implicit none
       type(c_ptr)        :: spmCopy_c
       type(c_ptr), value :: spm
     end function spmCopy_c
  end interface

  interface
     subroutine spmBase_c(spm, baseval) &
          bind(c, name='spmBase')
       use iso_c_binding
       import pastix_spm_t
       implicit none
       type(c_ptr),         value :: spm
       integer(kind=c_int), value :: baseval
     end subroutine spmBase_c
  end interface

  interface
     function spmFindBase_c(spm) &
          bind(c, name='spmFindBase')
       use iso_c_binding
       import pastix_int_t
       import pastix_spm_t
       implicit none
       integer(kind=pastix_int_t)   :: spmFindBase_c
       type(c_ptr),           value :: spm
     end function spmFindBase_c
  end interface

  interface
     function spmConvert_c(ofmttype, ospm) &
          bind(c, name='spmConvert')
       use iso_c_binding
       import pastix_spm_t
       implicit none
       integer(kind=c_int)        :: spmConvert_c
       integer(kind=c_int), value :: ofmttype
       type(c_ptr),         value :: ospm
     end function spmConvert_c
  end interface

  interface
     subroutine spmUpdateComputedFields_c(spm) &
          bind(c, name='spmUpdateComputedFields')
       use iso_c_binding
       import pastix_spm_t
       implicit none
       type(c_ptr), value :: spm
     end subroutine spmUpdateComputedFields_c
  end interface

  interface
     subroutine spmGenFakeValues_c(spm) &
          bind(c, name='spmGenFakeValues')
       use iso_c_binding
       import pastix_spm_t
       implicit none
       type(c_ptr), value :: spm
     end subroutine spmGenFakeValues_c
  end interface

  interface
     function spmNorm_c(ntype, spm) &
          bind(c, name='spmNorm')
       use iso_c_binding
       import pastix_spm_t
       implicit none
       real(kind=c_double)   :: spmNorm_c
       integer(c_int), value :: ntype
       type(c_ptr),    value :: spm
     end function spmNorm_c
  end interface

  interface
     function spmMatVec_c(trans, alpha, spm, x, beta, y) &
          bind(c, name='spmMatVec')
       use iso_c_binding
       import pastix_spm_t
       implicit none
       integer(kind=c_int)   :: spmMatVec_c
       integer(c_int), value :: trans
       type(c_ptr),    value :: alpha
       type(c_ptr),    value :: spm
       type(c_ptr),    value :: x
       type(c_ptr),    value :: beta
       type(c_ptr),    value :: y
     end function spmMatVec_c
  end interface

  interface
     subroutine spmScalMatrix_c(alpha, spm) &
          bind(c, name='spmScalMatrix')
       use iso_c_binding
       import pastix_spm_t
       implicit none
       complex(kind=c_double_complex), value :: alpha
       type(c_ptr),                    value :: spm
     end subroutine spmScalMatrix_c
  end interface

  interface
     subroutine spmScalVector_c(alpha, spm, x) &
          bind(c, name='spmScalVector')
       use iso_c_binding
       import pastix_spm_t
       implicit none
       real(kind=c_double), value :: alpha
       type(c_ptr),         value :: spm
       type(c_ptr),         value :: x
     end subroutine spmScalVector_c
  end interface

  interface
     function spmSort_c(spm) &
          bind(c, name='spmSort')
       use iso_c_binding
       import pastix_spm_t
       implicit none
       integer(kind=c_int)   :: spmSort_c
       type(c_ptr),    value :: spm
     end function spmSort_c
  end interface

  interface
     function spmMergeDuplicate_c(spm) &
          bind(c, name='spmMergeDuplicate')
       use iso_c_binding
       import pastix_int_t
       import pastix_spm_t
       implicit none
       integer(kind=pastix_int_t)   :: spmMergeDuplicate_c
       type(c_ptr),           value :: spm
     end function spmMergeDuplicate_c
  end interface

  interface
     function spmSymmetrize_c(spm) &
          bind(c, name='spmSymmetrize')
       use iso_c_binding
       import pastix_int_t
       import pastix_spm_t
       implicit none
       integer(kind=pastix_int_t)   :: spmSymmetrize_c
       type(c_ptr),           value :: spm
     end function spmSymmetrize_c
  end interface

  interface
     function spmCheckAndCorrect_c(spm) &
          bind(c, name='spmCheckAndCorrect')
       use iso_c_binding
       import pastix_spm_t
       implicit none
       type(c_ptr)        :: spmCheckAndCorrect_c
       type(c_ptr), value :: spm
     end function spmCheckAndCorrect_c
  end interface

  interface
     function spmGenRHS_c(type, nrhs, spm, x, ldx, b, ldb) &
          bind(c, name='spmGenRHS')
       use iso_c_binding
       import pastix_int_t
       import pastix_spm_t
       implicit none
       integer(kind=c_int)               :: spmGenRHS_c
       integer(c_int),             value :: type
       integer(kind=pastix_int_t), value :: nrhs
       type(c_ptr),                value :: spm
       type(c_ptr),                value :: x
       integer(kind=pastix_int_t), value :: ldx
       type(c_ptr),                value :: b
       integer(kind=pastix_int_t), value :: ldb
     end function spmGenRHS_c
  end interface

  interface
     function spmCheckAxb_c(nrhs, spm, x0, ldx0, b, ldb, x, ldx) &
          bind(c, name='spmCheckAxb')
       use iso_c_binding
       import pastix_int_t
       import pastix_spm_t
       implicit none
       integer(kind=c_int)               :: spmCheckAxb_c
       integer(kind=pastix_int_t), value :: nrhs
       type(c_ptr),                value :: spm
       type(c_ptr),                value :: x0
       integer(kind=pastix_int_t), value :: ldx0
       type(c_ptr),                value :: b
       integer(kind=pastix_int_t), value :: ldb
       type(c_ptr),                value :: x
       integer(kind=pastix_int_t), value :: ldx
     end function spmCheckAxb_c
  end interface

  interface
     function spmIntConvert_c(n, input) &
          bind(c, name='spmIntConvert')
       use iso_c_binding
       import pastix_int_t
       implicit none
       type(c_ptr)                       :: spmIntConvert_c
       integer(kind=pastix_int_t), value :: n
       type(c_ptr),                value :: input
     end function spmIntConvert_c
  end interface

  interface
     function spmLoad_c(spm, infile) &
          bind(c, name='spmLoad')
       use iso_c_binding
       import pastix_spm_t
       implicit none
       integer(kind=c_int)   :: spmLoad_c
       type(c_ptr),    value :: spm
       type(c_ptr),    value :: infile
     end function spmLoad_c
  end interface

  interface
     function spmSave_c(spm, outfile) &
          bind(c, name='spmSave')
       use iso_c_binding
       import pastix_spm_t
       implicit none
       integer(kind=c_int)   :: spmSave_c
       type(c_ptr),    value :: spm
       type(c_ptr),    value :: outfile
     end function spmSave_c
  end interface

  interface
     function spmReadDriver_c(driver, filename, spm, pastix_comm) &
          bind(c, name='spmReadDriver')
       use iso_c_binding
       import pastix_spm_t
       implicit none
       integer(kind=c_int)        :: spmReadDriver_c
       integer(c_int),      value :: driver
       type(c_ptr),         value :: filename
       type(c_ptr),         value :: spm
       integer(kind=c_int), value :: pastix_comm
     end function spmReadDriver_c
  end interface

  interface
     subroutine spm2Dense_c(spm) &
          bind(c, name='spm2Dense')
       use iso_c_binding
       import pastix_spm_t
       implicit none
       type(c_ptr), value :: spm
     end subroutine spm2Dense_c
  end interface

  interface
     subroutine spmPrint_c(spm, f) &
          bind(c, name='spmPrint')
       use iso_c_binding
       import pastix_spm_t
       implicit none
       type(c_ptr), value :: spm
       type(c_ptr), value :: f
     end subroutine spmPrint_c
  end interface

  interface
     subroutine spmPrintInfo_c(spm, f) &
          bind(c, name='spmPrintInfo')
       use iso_c_binding
       import pastix_spm_t
       implicit none
       type(c_ptr), value :: spm
       type(c_ptr), value :: f
     end subroutine spmPrintInfo_c
  end interface

  interface
     function spmExpand_c(spm) &
          bind(c, name='spmExpand')
       use iso_c_binding
       import pastix_spm_t
       implicit none
       type(c_ptr)        :: spmExpand_c
       type(c_ptr), value :: spm
     end function spmExpand_c
  end interface

  interface
     function spmDofExtend_c(spm, type, dof) &
          bind(c, name='spmDofExtend')
       use iso_c_binding
       import pastix_spm_t
       implicit none
       type(c_ptr)                :: spmDofExtend_c
       type(c_ptr),         value :: spm
       integer(kind=c_int), value :: type
       integer(kind=c_int), value :: dof
     end function spmDofExtend_c
  end interface

contains

  ! Wrappers of the C functions.
  subroutine spmNew(mtxtype, flttype, fmttype, n, nnz, colptr, rowptr, values, &
       loc2glob, dof, layout, dofs, spmo)
    use iso_c_binding
    implicit none
    integer(c_int),             intent(in)             :: mtxtype
    integer(c_int),             intent(in)             :: flttype
    integer(c_int),             intent(in)             :: fmttype
    integer(kind=pastix_int_t), intent(in)             :: n
    integer(kind=pastix_int_t), intent(in)             :: nnz
    integer(kind=pastix_int_t), intent(inout), target  :: colptr(*)
    integer(kind=pastix_int_t), intent(inout), target  :: rowptr(*)
    type(c_ptr),                intent(inout), target  :: values(*)
    integer(kind=pastix_int_t), intent(inout), target  :: loc2glob(*)
    integer(kind=pastix_int_t), intent(in)             :: dof
    integer(c_int),             intent(in)             :: layout
    integer(kind=pastix_int_t), intent(inout), target  :: dofs(*)
    type(pastix_spm_t),         intent(out),   pointer :: spmo

    call c_f_pointer(spmNew_c(mtxtype, flttype, fmttype, n, nnz, c_loc(colptr), &
         c_loc(rowptr), c_loc(values), c_loc(loc2glob), dof, layout, &
         c_loc(dofs)), spmo)
  end subroutine spmNew

  subroutine spmInit(spm)
    use iso_c_binding
    implicit none
    type(pastix_spm_t), intent(inout), target :: spm

    call spmInit_c(c_loc(spm))
  end subroutine spmInit

  subroutine spmExit(spm)
    use iso_c_binding
    implicit none
    type(pastix_spm_t), intent(inout), target :: spm

    call spmExit_c(c_loc(spm))
  end subroutine spmExit

  subroutine spmFree(spm)
    use iso_c_binding
    implicit none
    type(pastix_spm_t), intent(inout), target :: spm

    call spmFree_c(c_loc(spm))
  end subroutine spmFree

  subroutine spmCopy(spm, spmo)
    use iso_c_binding
    implicit none
    type(pastix_spm_t), intent(in),  target  :: spm
    type(pastix_spm_t), intent(out), pointer :: spmo

    call c_f_pointer(spmCopy_c(c_loc(spm)), spmo)
  end subroutine spmCopy

  subroutine spmBase(spm, baseval)
    use iso_c_binding
    implicit none
    type(pastix_spm_t),  intent(inout), target :: spm
    integer(kind=c_int), intent(in)            :: baseval

    call spmBase_c(c_loc(spm), baseval)
  end subroutine spmBase

  subroutine spmFindBase(spm, value)
    use iso_c_binding
    implicit none
    type(pastix_spm_t),         intent(in), target :: spm
    integer(kind=pastix_int_t), intent(out)        :: value

    value = spmFindBase_c(c_loc(spm))
  end subroutine spmFindBase

  subroutine spmConvert(ofmttype, ospm, info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(in)            :: ofmttype
    type(pastix_spm_t),  intent(inout), target :: ospm
    integer(kind=c_int), intent(out)           :: info

    info = spmConvert_c(ofmttype, c_loc(ospm))
  end subroutine spmConvert

  subroutine spmUpdateComputedFields(spm)
    use iso_c_binding
    implicit none
    type(pastix_spm_t), intent(inout), target :: spm

    call spmUpdateComputedFields_c(c_loc(spm))
  end subroutine spmUpdateComputedFields

  subroutine spmGenFakeValues(spm)
    use iso_c_binding
    implicit none
    type(pastix_spm_t), intent(inout), target :: spm

    call spmGenFakeValues_c(c_loc(spm))
  end subroutine spmGenFakeValues

  subroutine spmNorm(ntype, spm, value)
    use iso_c_binding
    implicit none
    integer(c_int),      intent(in)         :: ntype
    type(pastix_spm_t),  intent(in), target :: spm
    real(kind=c_double), intent(out)        :: value

    value = spmNorm_c(ntype, c_loc(spm))
  end subroutine spmNorm

  subroutine spmMatVec(trans, alpha, spm, x, beta, y, info)
    use iso_c_binding
    implicit none
    integer(c_int),      intent(in)            :: trans
    type(c_ptr),         intent(in),    target :: alpha
    type(pastix_spm_t),  intent(in),    target :: spm
    type(c_ptr),         intent(in),    target :: x
    type(c_ptr),         intent(in),    target :: beta
    type(c_ptr),         intent(inout), target :: y
    integer(kind=c_int), intent(out)           :: info

    info = spmMatVec_c(trans, c_loc(alpha), c_loc(spm), c_loc(x), c_loc(beta), &
         c_loc(y))
  end subroutine spmMatVec

  subroutine spmScalMatrix(alpha, spm)
    use iso_c_binding
    implicit none
    complex(kind=c_double_complex), intent(in)            :: alpha
    type(pastix_spm_t),             intent(inout), target :: spm

    call spmScalMatrix_c(alpha, c_loc(spm))
  end subroutine spmScalMatrix

  subroutine spmScalVector(alpha, spm, x)
    use iso_c_binding
    implicit none
    real(kind=c_double), intent(in)            :: alpha
    type(pastix_spm_t),  intent(inout), target :: spm
    type(c_ptr),         intent(inout), target :: x

    call spmScalVector_c(alpha, c_loc(spm), c_loc(x))
  end subroutine spmScalVector

  subroutine spmSort(spm, info)
    use iso_c_binding
    implicit none
    type(pastix_spm_t),  intent(inout), target :: spm
    integer(kind=c_int), intent(out)           :: info

    info = spmSort_c(c_loc(spm))
  end subroutine spmSort

  subroutine spmMergeDuplicate(spm, value)
    use iso_c_binding
    implicit none
    type(pastix_spm_t),         intent(inout), target :: spm
    integer(kind=pastix_int_t), intent(out)           :: value

    value = spmMergeDuplicate_c(c_loc(spm))
  end subroutine spmMergeDuplicate

  subroutine spmSymmetrize(spm, value)
    use iso_c_binding
    implicit none
    type(pastix_spm_t),         intent(inout), target :: spm
    integer(kind=pastix_int_t), intent(out)           :: value

    value = spmSymmetrize_c(c_loc(spm))
  end subroutine spmSymmetrize

  subroutine spmCheckAndCorrect(spm, spmo)
    use iso_c_binding
    implicit none
    type(pastix_spm_t), intent(inout), target  :: spm
    type(pastix_spm_t), intent(out),   pointer :: spmo

    call c_f_pointer(spmCheckAndCorrect_c(c_loc(spm)), spmo)
  end subroutine spmCheckAndCorrect

  subroutine spmGenRHS(type, nrhs, spm, x, ldx, b, ldb, info)
    use iso_c_binding
    implicit none
    integer(c_int),             intent(in)            :: type
    integer(kind=pastix_int_t), intent(in)            :: nrhs
    type(pastix_spm_t),         intent(in),    target :: spm
    type(c_ptr),                intent(inout), target :: x
    integer(kind=pastix_int_t), intent(in)            :: ldx
    type(c_ptr),                intent(inout), target :: b
    integer(kind=pastix_int_t), intent(in)            :: ldb
    integer(kind=c_int),        intent(out)           :: info

    info = spmGenRHS_c(type, nrhs, c_loc(spm), c_loc(x), ldx, c_loc(b), ldb)
  end subroutine spmGenRHS

  subroutine spmCheckAxb(nrhs, spm, x0, ldx0, b, ldb, x, ldx, info)
    use iso_c_binding
    implicit none
    integer(kind=pastix_int_t), intent(in)            :: nrhs
    type(pastix_spm_t),         intent(in),    target :: spm
    type(c_ptr),                intent(inout), target :: x0
    integer(kind=pastix_int_t), intent(in)            :: ldx0
    type(c_ptr),                intent(inout), target :: b
    integer(kind=pastix_int_t), intent(in)            :: ldb
    type(c_ptr),                intent(in),    target :: x
    integer(kind=pastix_int_t), intent(in)            :: ldx
    integer(kind=c_int),        intent(out)           :: info

    info = spmCheckAxb_c(nrhs, c_loc(spm), c_loc(x0), ldx0, c_loc(b), ldb, &
         c_loc(x), ldx)
  end subroutine spmCheckAxb

  subroutine spmIntConvert(n, input, value)
    use iso_c_binding
    implicit none
    integer(kind=pastix_int_t), intent(in)             :: n
    integer(kind=c_int),        intent(inout), target  :: input
    integer(kind=pastix_int_t), intent(out),   pointer :: value

    call c_f_pointer(spmIntConvert_c(n, c_loc(input)), value)
  end subroutine spmIntConvert

  subroutine spmLoad(spm, info)
    use iso_c_binding
    implicit none
    type(pastix_spm_t),  intent(inout), target :: spm
    integer(kind=c_int), intent(out)           :: info

    info = spmLoad_c(c_loc(spm), c_null_ptr)
  end subroutine spmLoad

  subroutine spmSave(spm, info)
    use iso_c_binding
    implicit none
    type(pastix_spm_t),  intent(in), target :: spm
    integer(kind=c_int), intent(out)        :: info

    info = spmSave_c(c_loc(spm), c_null_ptr)
  end subroutine spmSave

  subroutine spmReadDriver(driver, filename, spm, pastix_comm, info)
    use iso_c_binding
    implicit none
    integer(c_int),         intent(in)            :: driver
    character(kind=c_char), intent(in),    target :: filename
    type(pastix_spm_t),     intent(inout), target :: spm
    integer(kind=c_int),    intent(in)            :: pastix_comm
    integer(kind=c_int),    intent(out)           :: info

    info = spmReadDriver_c(driver, c_loc(filename), c_loc(spm), pastix_comm)
  end subroutine spmReadDriver

  subroutine spm2Dense(spm)
    use iso_c_binding
    implicit none
    type(pastix_spm_t), intent(in), target :: spm

    call spm2Dense_c(c_loc(spm))
  end subroutine spm2Dense

  subroutine spmPrint(spm)
    use iso_c_binding
    implicit none
    type(pastix_spm_t), intent(in), target :: spm

    call spmPrint_c(c_loc(spm), c_null_ptr)
  end subroutine spmPrint

  subroutine spmPrintInfo(spm)
    use iso_c_binding
    implicit none
    type(pastix_spm_t), intent(in), target :: spm

    call spmPrintInfo_c(c_loc(spm), c_null_ptr)
  end subroutine spmPrintInfo

  subroutine spmExpand(spm, spmo)
    use iso_c_binding
    implicit none
    type(pastix_spm_t), intent(in),  target  :: spm
    type(pastix_spm_t), intent(out), pointer :: spmo

    call c_f_pointer(spmExpand_c(c_loc(spm)), spmo)
  end subroutine spmExpand

  subroutine spmDofExtend(spm, type, dof, spmo)
    use iso_c_binding
    implicit none
    type(pastix_spm_t),  intent(in),  target  :: spm
    integer(kind=c_int), intent(in)           :: type
    integer(kind=c_int), intent(in)           :: dof
    type(pastix_spm_t),  intent(out), pointer :: spmo

    call c_f_pointer(spmDofExtend_c(c_loc(spm), type, dof), spmo)
  end subroutine spmDofExtend


end module spmf
