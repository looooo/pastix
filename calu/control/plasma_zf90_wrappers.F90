!
!     Copyright © 2011 The Numerical Algorithms Group Ltd. All rights reserved.
!   
!     Redistribution and use in source and binary forms, with or without
!     modification, are permitted provided that the following conditions are
!     met:
!     - Redistributions of source code must retain the above copyright notice,
!       this list of conditions, and the following disclaimer.
!     - Redistributions in binary form must reproduce the above copyright
!       notice, this list of conditions and the following disclaimer listed in
!       this license in the documentation and/or other materials provided with
!       the distribution.
!     - Neither the name of the copyright holders nor the names of its
!       contributors may be used to endorse or promote products derived from
!       this software without specific prior written permission.
!     
!     This software is provided by the copyright holders and contributors "as
!     is" and any express or implied warranties, including, but not limited
!     to, the implied warranties of merchantability and fitness for a
!     particular purpose are disclaimed. in no event shall the copyright owner
!     or contributors be liable for any direct, indirect, incidental, special,
!     exemplary, or consequential damages (including, but not limited to,
!     procurement of substitute goods or services; loss of use, data, or
!     profits; or business interruption) however caused and on any theory of
!     liability, whether in contract, strict liability, or tort (including
!     negligence or otherwise) arising in any way out of the use of this
!     software, even if advised of the possibility of such damage.
!
!
!
! @file plasma_zf90_wrappers.F90
!
!  PLASMA fortran wrapper for BLAS and LAPACK subroutines.
!  PLASMA is a software package provided by Univ. of Tennessee,
!  Univ. of California Berkeley and Univ. of Colorado Denver
!
! @version 2.4.2
! @author Numerical Algorithm Group
! @author Mathieu Faverge
! @date 2011-09-15
! @precisions normal z -> c d s
!
!
! Wrappers to PLASMA functions are provided for the following BLAS
! subroutines since the PLASMA and BLAS interfaces match exactly:
! ZGEMM  PLASMA_zgemm
! ZHEMM  PLASMA_zhemm
! ZHER2K PLASMA_zher2k
! ZHERK  PLASMA_zherk
! ZSYMM  PLASMA_zsymm
! ZSYR2K PLASMA_zsyr2k
! ZSYRK  PLASMA_zsyrk
! ZTRMM  PLASMA_ztrmm
! ZTRSM  PLASMA_ztrsm
!
! Wrappers to PLASMA functions are provided for the following LAPACK
! subroutines since the PLASMA and LAPACK interfaces match exactly:
! ZGESV  PLASMA_zgesv
! ZGETRF PLASMA_zgetrf
! ZGETRS PLASMA_zgetrs
! ZHEGST PLASMA_zhegst
! ZLASWP PLASMA_zlaswp
! ZLAUUM PLASMA_zlauum
! ZPOSV  PLASMA_zposv
! ZPOTRF PLASMA_zpotrf
! ZPOTRI PLASMA_zpotri
! ZPOTRS PLASMA_zpotrs
! ZTRTRI PLASMA_ztrtri
! ZLACPY PLASMA_zlacpy
! ZLASET PLASMA_zlaset
#define PRECISION_z

      subroutine plasma_wrap_ZGESV(N,NRHS,A,LDA,IPIV,B,LDB,INFO)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: N
            integer, intent(in) :: NRHS
            integer, intent(out) :: INFO
            integer, intent(out), target :: IPIV(*)
            complex(kind=wp), intent(inout), target :: A(LDA,*)
            complex(kind=wp), intent(inout), target :: B(LDB,*)
            if (.not. plasma_initialized) call plasma_init(24,INFO)
            ! write(*,*) " Calling PLASMA_ZGESV"
            call PLASMA_ZGESV(N,NRHS,A,LDA,IPIV,B,LDB,INFO)
      end subroutine plasma_wrap_ZGESV

      subroutine plasma_wrap_ZGETRF(M,N,A,LDA,IPIV,INFO)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: M
            integer, intent(in) :: N
            integer, intent(out) :: INFO
            integer, intent(out), target :: IPIV(*)
            complex(kind=wp), intent(inout), target :: A(LDA,*)
            if (.not. plasma_initialized) call plasma_init(24,INFO)
            ! write(*,*) " Calling PLASMA_ZGETRF"
            call PLASMA_ZGETRF(M,N,A,LDA,IPIV,INFO)
      end subroutine plasma_wrap_ZGETRF

      subroutine plasma_wrap_ZGETRS(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: N
            integer, intent(in) :: NRHS
            integer, intent(in), target :: IPIV(*)
            integer, intent(out) :: INFO
            character, intent(in) :: TRANS
            complex(kind=wp), intent(in), target :: A(LDA,*)
            complex(kind=wp), intent(inout), target :: B(LDB,*)
            integer :: local_TRANS
            if(TRANS=='N' .or. TRANS=='n')then
               local_TRANS = PlasmaNoTrans
            else if(TRANS=='T' .or. TRANS=='t')then
               local_TRANS = PlasmaTrans
            else if(TRANS=='C' .or. TRANS=='c')then
               local_TRANS = PlasmaConjTrans
            else
               local_TRANS = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,INFO)
            ! write(*,*) " Calling PLASMA_ZGETRS"
            call PLASMA_ZGETRS(local_TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
      end subroutine plasma_wrap_ZGETRS

      subroutine plasma_wrap_ZHEGST(ITYPE,UPLO,N,A,LDA,B,LDB,INFO)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: ITYPE
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: N
            integer, intent(out) :: INFO
            character, intent(in) :: UPLO
            complex(kind=wp), intent(in), target :: B(LDB,*)
            complex(kind=wp), intent(inout), target :: A(LDA,*)
            integer :: local_UPLO
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,INFO)
            ! write(*,*) " Calling PLASMA_ZHEGST"
            call PLASMA_ZHEGST(ITYPE,local_UPLO,N,A,LDA,B,LDB,INFO)
      end subroutine plasma_wrap_ZHEGST

      subroutine plasma_wrap_ZLASWP(N,A,LDA,K1,K2,IPIV,INCX)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: INCX
            integer, intent(in) :: K1
            integer, intent(in) :: K2
            integer, intent(in) :: LDA
            integer, intent(in) :: N
            integer, intent(in), target :: IPIV(*)
            complex(kind=wp), intent(inout), target :: A(LDA,*)
            integer :: local_ret
            if (.not. plasma_initialized) call plasma_init(24,local_ret)
            ! write(*,*) " Calling PLASMA_ZLASWP"
            call PLASMA_ZLASWP(N,A,LDA,K1,K2,IPIV,INCX,local_ret)
      end subroutine plasma_wrap_ZLASWP

      subroutine plasma_wrap_ZLAUUM(UPLO,N,A,LDA,INFO)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: N
            integer, intent(out) :: INFO
            character, intent(in) :: UPLO
            complex(kind=wp), intent(inout), target :: A(LDA,*)
            integer :: local_UPLO
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,INFO)
            ! write(*,*) " Calling PLASMA_ZLAUUM"
            call PLASMA_ZLAUUM(local_UPLO,N,A,LDA,INFO)
      end subroutine plasma_wrap_ZLAUUM

      subroutine plasma_wrap_ZPOSV(UPLO,N,NRHS,A,LDA,B,LDB,INFO)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: N
            integer, intent(in) :: NRHS
            integer, intent(out) :: INFO
            character, intent(in) :: UPLO
            complex(kind=wp), intent(inout), target :: A(LDA,*)
            complex(kind=wp), intent(inout), target :: B(LDB,*)
            integer :: local_UPLO
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,INFO)
            ! write(*,*) " Calling PLASMA_ZPOSV"
            call PLASMA_ZPOSV(local_UPLO,N,NRHS,A,LDA,B,LDB,INFO)
      end subroutine plasma_wrap_ZPOSV

      subroutine plasma_wrap_ZPOTRF(UPLO,N,A,LDA,INFO)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: N
            integer, intent(out) :: INFO
            character, intent(in) :: UPLO
            complex(kind=wp), intent(inout), target :: A(LDA,*)
            integer :: local_UPLO
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,INFO)
            ! write(*,*) " Calling PLASMA_ZPOTRF"
            call PLASMA_ZPOTRF(local_UPLO,N,A,LDA,INFO)
      end subroutine plasma_wrap_ZPOTRF

      subroutine plasma_wrap_ZPOTRI(UPLO,N,A,LDA,INFO)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: N
            integer, intent(out) :: INFO
            character, intent(in) :: UPLO
            complex(kind=wp), intent(inout), target :: A(LDA,*)
            integer :: local_UPLO
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,INFO)
            ! write(*,*) " Calling PLASMA_ZPOTRI"
            call PLASMA_ZPOTRI(local_UPLO,N,A,LDA,INFO)
      end subroutine plasma_wrap_ZPOTRI

      subroutine plasma_wrap_ZPOTRS(UPLO,N,NRHS,A,LDA,B,LDB,INFO)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: N
            integer, intent(in) :: NRHS
            integer, intent(out) :: INFO
            character, intent(in) :: UPLO
            complex(kind=wp), intent(in), target :: A(LDA,*)
            complex(kind=wp), intent(inout), target :: B(LDB,*)
            integer :: local_UPLO
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,INFO)
            ! write(*,*) " Calling PLASMA_ZPOTRS"
            call PLASMA_ZPOTRS(local_UPLO,N,NRHS,A,LDA,B,LDB,INFO)
      end subroutine plasma_wrap_ZPOTRS

      subroutine plasma_wrap_ZTRTRI(UPLO,DIAG,N,A,LDA,INFO)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: N
            integer, intent(out) :: INFO
            character, intent(in) :: DIAG
            character, intent(in) :: UPLO
            complex(kind=wp), intent(inout), target :: A(LDA,*)
            integer :: local_DIAG
            integer :: local_UPLO
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if(DIAG=='U' .or. DIAG=='u')then
               local_DIAG = PlasmaUnit
            else if(DIAG=='N' .or. DIAG=='n')then
               local_DIAG = PlasmaNonUnit
            else
               local_DIAG = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,INFO)
            ! write(*,*) " Calling PLASMA_ZTRTRI"
            call PLASMA_ZTRTRI(local_UPLO,local_DIAG,N,A,LDA,INFO)
      end subroutine plasma_wrap_ZTRTRI

      subroutine plasma_wrap_ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: K
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: LDC
            integer, intent(in) :: M
            integer, intent(in) :: N
            character, intent(in) :: TRANSA
            character, intent(in) :: TRANSB
            complex(kind=wp), intent(in) :: ALPHA
            complex(kind=wp), intent(in) :: BETA
            complex(kind=wp), intent(in), target :: A(LDA,*)
            complex(kind=wp), intent(in), target :: B(LDB,*)
            complex(kind=wp), intent(inout), target :: C(LDC,*)
            integer :: local_TRANSA
            integer :: local_TRANSB
            integer :: local_ret
            if(TRANSA=='N' .or. TRANSA=='n')then
               local_TRANSA = PlasmaNoTrans
            else if(TRANSA=='T' .or. TRANSA=='t')then
               local_TRANSA = PlasmaTrans
            else if(TRANSA=='C' .or. TRANSA=='c')then
               local_TRANSA = PlasmaConjTrans
            else
               local_TRANSA = -1
            end if
            if(TRANSB=='N' .or. TRANSB=='n')then
               local_TRANSB = PlasmaNoTrans
            else if(TRANSB=='T' .or. TRANSB=='t')then
               local_TRANSB = PlasmaTrans
            else if(TRANSB=='C' .or. TRANSB=='c')then
               local_TRANSB = PlasmaConjTrans
            else
               local_TRANSB = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,local_ret)
            ! write(*,*) " Calling PLASMA_ZGEMM"
            call PLASMA_ZGEMM(local_TRANSA,local_TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC,local_ret)
      end subroutine plasma_wrap_ZGEMM

#if defined(PRECISION_z) || defined(PRECISION_c)
      subroutine plasma_wrap_ZHEMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: LDC
            integer, intent(in) :: M
            integer, intent(in) :: N
            character, intent(in) :: SIDE
            character, intent(in) :: UPLO
            complex(kind=wp), intent(in) :: ALPHA
            complex(kind=wp), intent(in) :: BETA
            complex(kind=wp), intent(in), target :: A(LDA,*)
            complex(kind=wp), intent(in), target :: B(LDB,*)
            complex(kind=wp), intent(inout), target :: C(LDC,*)
            integer :: local_SIDE
            integer :: local_UPLO
            integer :: local_ret
            if(SIDE=='L' .or. SIDE=='l')then
               local_SIDE = PlasmaLeft
            else if(SIDE=='R' .or. SIDE=='r')then
               local_SIDE = PlasmaRight
            else
               local_SIDE = -1
            end if
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,local_ret)
            ! write(*,*) " Calling PLASMA_ZHEMM"
            call PLASMA_ZHEMM(local_SIDE,local_UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC,local_ret)
      end subroutine plasma_wrap_ZHEMM

      subroutine plasma_wrap_ZHER2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: K
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: LDC
            integer, intent(in) :: N
            character, intent(in) :: TRANS
            character, intent(in) :: UPLO
            complex(kind=wp), intent(in) :: ALPHA
            complex(kind=wp), intent(in), target :: A(LDA,*)
            complex(kind=wp), intent(in), target :: B(LDB,*)
            complex(kind=wp), intent(inout), target :: C(LDC,*)
            double precision, intent(in) :: BETA
            integer :: local_TRANS
            integer :: local_UPLO
            integer :: local_ret
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if(TRANS=='N' .or. TRANS=='n')then
               local_TRANS = PlasmaNoTrans
            else if(TRANS=='T' .or. TRANS=='t')then
               local_TRANS = PlasmaTrans
            else if(TRANS=='C' .or. TRANS=='c')then
               local_TRANS = PlasmaConjTrans
            else
               local_TRANS = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,local_ret)
            ! write(*,*) " Calling PLASMA_ZHER2K"
            call PLASMA_ZHER2K(local_UPLO,local_TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC,local_ret)
      end subroutine plasma_wrap_ZHER2K

      subroutine plasma_wrap_ZHERK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: K
            integer, intent(in) :: LDA
            integer, intent(in) :: LDC
            integer, intent(in) :: N
            character, intent(in) :: TRANS
            character, intent(in) :: UPLO
            complex(kind=wp), intent(in), target :: A(LDA,*)
            complex(kind=wp), intent(inout), target :: C(LDC,*)
            double precision, intent(in) :: ALPHA
            double precision, intent(in) :: BETA
            integer :: local_TRANS
            integer :: local_UPLO
            integer :: local_ret
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if(TRANS=='N' .or. TRANS=='n')then
               local_TRANS = PlasmaNoTrans
            else if(TRANS=='T' .or. TRANS=='t')then
               local_TRANS = PlasmaTrans
            else if(TRANS=='C' .or. TRANS=='c')then
               local_TRANS = PlasmaConjTrans
            else
               local_TRANS = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,local_ret)
            ! write(*,*) " Calling PLASMA_ZHERK"
            call PLASMA_ZHERK(local_UPLO,local_TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC,local_ret)
      end subroutine plasma_wrap_ZHERK
#endif

      subroutine plasma_wrap_ZSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: LDC
            integer, intent(in) :: M
            integer, intent(in) :: N
            character, intent(in) :: SIDE
            character, intent(in) :: UPLO
            complex(kind=wp), intent(in) :: ALPHA
            complex(kind=wp), intent(in) :: BETA
            complex(kind=wp), intent(in), target :: A(LDA,*)
            complex(kind=wp), intent(in), target :: B(LDB,*)
            complex(kind=wp), intent(inout), target :: C(LDC,*)
            integer :: local_SIDE
            integer :: local_UPLO
            integer :: local_ret
            if(SIDE=='L' .or. SIDE=='l')then
               local_SIDE = PlasmaLeft
            else if(SIDE=='R' .or. SIDE=='r')then
               local_SIDE = PlasmaRight
            else
               local_SIDE = -1
            end if
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,local_ret)
            ! write(*,*) " Calling PLASMA_ZSYMM"
            call PLASMA_ZSYMM(local_SIDE,local_UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC,local_ret)
      end subroutine plasma_wrap_ZSYMM

      subroutine plasma_wrap_ZSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: K
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: LDC
            integer, intent(in) :: N
            character, intent(in) :: TRANS
            character, intent(in) :: UPLO
            complex(kind=wp), intent(in) :: ALPHA
            complex(kind=wp), intent(in) :: BETA
            complex(kind=wp), intent(in), target :: A(LDA,*)
            complex(kind=wp), intent(in), target :: B(LDB,*)
            complex(kind=wp), intent(inout), target :: C(LDC,*)
            integer :: local_TRANS
            integer :: local_UPLO
            integer :: local_ret
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if(TRANS=='N' .or. TRANS=='n')then
               local_TRANS = PlasmaNoTrans
            else if(TRANS=='T' .or. TRANS=='t')then
               local_TRANS = PlasmaTrans
            else if(TRANS=='C' .or. TRANS=='c')then
               local_TRANS = PlasmaConjTrans
            else
               local_TRANS = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,local_ret)
            ! write(*,*) " Calling PLASMA_ZSYR2K"
            call PLASMA_ZSYR2K(local_UPLO,local_TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC,local_ret)
      end subroutine plasma_wrap_ZSYR2K

      subroutine plasma_wrap_ZSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: K
            integer, intent(in) :: LDA
            integer, intent(in) :: LDC
            integer, intent(in) :: N
            character, intent(in) :: TRANS
            character, intent(in) :: UPLO
            complex(kind=wp), intent(in) :: ALPHA
            complex(kind=wp), intent(in) :: BETA
            complex(kind=wp), intent(in), target :: A(LDA,*)
            complex(kind=wp), intent(inout), target :: C(LDC,*)
            integer :: local_TRANS
            integer :: local_UPLO
            integer :: local_ret
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if(TRANS=='N' .or. TRANS=='n')then
               local_TRANS = PlasmaNoTrans
            else if(TRANS=='T' .or. TRANS=='t')then
               local_TRANS = PlasmaTrans
            else if(TRANS=='C' .or. TRANS=='c')then
               local_TRANS = PlasmaConjTrans
            else
               local_TRANS = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,local_ret)
            ! write(*,*) " Calling PLASMA_ZSYRK"
            call PLASMA_ZSYRK(local_UPLO,local_TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC,local_ret)
      end subroutine plasma_wrap_ZSYRK

      subroutine plasma_wrap_ZTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: M
            integer, intent(in) :: N
            character, intent(in) :: DIAG
            character, intent(in) :: SIDE
            character, intent(in) :: TRANSA
            character, intent(in) :: UPLO
            complex(kind=wp), intent(in) :: ALPHA
            complex(kind=wp), intent(in), target :: A(LDA,*)
            complex(kind=wp), intent(inout), target :: B(LDB,*)
            integer :: local_DIAG
            integer :: local_SIDE
            integer :: local_TRANSA
            integer :: local_UPLO
            integer :: local_ret
            if(SIDE=='L' .or. SIDE=='l')then
               local_SIDE = PlasmaLeft
            else if(SIDE=='R' .or. SIDE=='r')then
               local_SIDE = PlasmaRight
            else
               local_SIDE = -1
            end if
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if(TRANSA=='N' .or. TRANSA=='n')then
               local_TRANSA = PlasmaNoTrans
            else if(TRANSA=='T' .or. TRANSA=='t')then
               local_TRANSA = PlasmaTrans
            else if(TRANSA=='C' .or. TRANSA=='c')then
               local_TRANSA = PlasmaConjTrans
            else
               local_TRANSA = -1
            end if
            if(DIAG=='U' .or. DIAG=='u')then
               local_DIAG = PlasmaUnit
            else if(DIAG=='N' .or. DIAG=='n')then
               local_DIAG = PlasmaNonUnit
            else
               local_DIAG = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,local_ret)
            ! write(*,*) " Calling PLASMA_ZTRMM"
            call PLASMA_ZTRMM(local_SIDE,local_UPLO,local_TRANSA,local_DIAG,M,N,ALPHA,A,LDA,B,LDB,local_ret)
      end subroutine plasma_wrap_ZTRMM

      subroutine plasma_wrap_ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: M
            integer, intent(in) :: N
            character, intent(in) :: DIAG
            character, intent(in) :: SIDE
            character, intent(in) :: TRANSA
            character, intent(in) :: UPLO
            complex(kind=wp), intent(in) :: ALPHA
            complex(kind=wp), intent(in), target :: A(LDA,*)
            complex(kind=wp), intent(inout), target :: B(LDB,*)
            integer :: local_DIAG
            integer :: local_SIDE
            integer :: local_TRANSA
            integer :: local_UPLO
            integer :: local_ret
            if(SIDE=='L' .or. SIDE=='l')then
               local_SIDE = PlasmaLeft
            else if(SIDE=='R' .or. SIDE=='r')then
               local_SIDE = PlasmaRight
            else
               local_SIDE = -1
            end if
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if(TRANSA=='N' .or. TRANSA=='n')then
               local_TRANSA = PlasmaNoTrans
            else if(TRANSA=='T' .or. TRANSA=='t')then
               local_TRANSA = PlasmaTrans
            else if(TRANSA=='C' .or. TRANSA=='c')then
               local_TRANSA = PlasmaConjTrans
            else
               local_TRANSA = -1
            end if
            if(DIAG=='U' .or. DIAG=='u')then
               local_DIAG = PlasmaUnit
            else if(DIAG=='N' .or. DIAG=='n')then
               local_DIAG = PlasmaNonUnit
            else
               local_DIAG = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,local_ret)
            ! write(*,*) " Calling PLASMA_ZTRSM"
            call PLASMA_ZTRSM(local_SIDE,local_UPLO,local_TRANSA,local_DIAG,M,N,ALPHA,A,LDA,B,LDB,local_ret)
      end subroutine plasma_wrap_ZTRSM

      subroutine plasma_wrap_ZLACPY(UPLO,M,N,A,LDA,B,LDB)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: M
            integer, intent(in) :: N
            character, intent(in) :: UPLO
            complex(kind=wp), intent(in), target :: A(LDA,*)
            complex(kind=wp), intent(out), target :: B(LDB,*)
            integer :: local_UPLO
            integer :: local_ret
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,local_ret)
            ! write(*,*) " Calling PLASMA_ZLACPY"
            call PLASMA_ZLACPY(local_UPLO,M,N,A,LDA,B,LDB,local_ret)
      end subroutine plasma_wrap_ZLACPY

      subroutine plasma_wrap_ZLASET(UPLO,M,N,ALPHA,BETA,A,LDA)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: M
            integer, intent(in) :: N
            character, intent(in) :: UPLO
            complex(kind=wp), intent(in) :: ALPHA
            complex(kind=wp), intent(in) :: BETA
            complex(kind=wp), intent(inout), target :: A(LDA,*)
            integer :: local_UPLO
            integer :: local_ret
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,local_ret)
            ! write(*,*) " Calling PLASMA_ZLASET"
            call PLASMA_ZLASET(local_UPLO,M,N,ALPHA,BETA,A,LDA,local_ret)
      end subroutine plasma_wrap_ZLASET
