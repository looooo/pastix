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
! @file plasma_zf90.F90
!
!  PLASMA Fortran 90 interfaces using Fortran 2003 ISO C bindings
!  PLASMA is a software package provided by Univ. of Tennessee,
!  Univ. of California Berkeley and Univ. of Colorado Denver
!
! @version 2.4.6
! @author Numerical Algorithms Group
! @author Mathieu Faverge
! @date 2011-12-15
! @precisions normal z -> c d s
!
#define PRECISION_z

module plasma_z
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  FORTRAN API - math functions (simple interface)
    !
      interface
         function PLASMA_zLapack_to_Tile_c(Af77,LDA,A) &
          & bind(c, name='PLASMA_zLapack_to_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zLapack_to_Tile_c
            type(c_ptr), value :: Af77
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: A
          end function PLASMA_zLapack_to_Tile_c
      end interface

      interface
         function PLASMA_zLapack_to_Tile_Async_c(Af77,LDA,A,sequence,request) &
          & bind(c, name='PLASMA_zLapack_to_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zLapack_to_Tile_Async_c
            type(c_ptr), value :: Af77
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: A
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zLapack_to_Tile_Async_c
      end interface

      interface
         function PLASMA_zTile_to_Lapack_c(A,Af77,LDA) &
          & bind(c, name='PLASMA_zTile_to_Lapack')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zTile_to_Lapack_c
            type(c_ptr), value :: A
            type(c_ptr), value :: Af77
            integer(kind=c_int), value :: LDA
          end function PLASMA_zTile_to_Lapack_c
      end interface

      interface
         function PLASMA_zTile_to_Lapack_Async_c(A,Af77,LDA,sequence,request) &
          & bind(c, name='PLASMA_zTile_to_Lapack_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zTile_to_Lapack_Async_c
            type(c_ptr), value :: A
            type(c_ptr), value :: Af77
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zTile_to_Lapack_Async_c
      end interface

      interface
         function PLASMA_zgebrd_c(M,N,A,LDA,D,E,descT) &
          & bind(c, name='PLASMA_zgebrd')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgebrd_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: D
            type(c_ptr), value :: E
            type(c_ptr), value :: descT
          end function PLASMA_zgebrd_c
      end interface

      interface
         function PLASMA_zgecfi_c(m,n,A,fin,imb,inb,fout,omb,onb) &
          & bind(c, name='PLASMA_zgecfi')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgecfi_c
            integer(kind=c_int), value :: m
            integer(kind=c_int), value :: n
            type(c_ptr), value :: A
            integer(kind=c_int), value :: fin
            integer(kind=c_int), value :: imb
            integer(kind=c_int), value :: inb
            integer(kind=c_int), value :: fout
            integer(kind=c_int), value :: omb
            integer(kind=c_int), value :: onb
          end function PLASMA_zgecfi_c
      end interface

      interface
         function PLASMA_zgecfi_Async_c(m,n,A,f_in,imb,inb,f_out,omb,onb,sequence,request) &
          & bind(c, name='PLASMA_zgecfi_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgecfi_Async_c
            integer(kind=c_int), value :: m
            integer(kind=c_int), value :: n
            type(c_ptr), value :: A
            integer(kind=c_int), value :: f_in
            integer(kind=c_int), value :: imb
            integer(kind=c_int), value :: inb
            integer(kind=c_int), value :: f_out
            integer(kind=c_int), value :: omb
            integer(kind=c_int), value :: onb
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zgecfi_Async_c
      end interface

      interface
         function PLASMA_zgelqf_c(M,N,A,LDA,descT) &
          & bind(c, name='PLASMA_zgelqf')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgelqf_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: descT
          end function PLASMA_zgelqf_c
      end interface

      interface
         function PLASMA_zgelqs_c(M,N,NRHS,A,LDA,descT,B,LDB) &
          & bind(c, name='PLASMA_zgelqs')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgelqs_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: descT
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
          end function PLASMA_zgelqs_c
      end interface

      interface
         function PLASMA_zgels_c(trans,M,N,NRHS,A,LDA,descT,B,LDB) &
          & bind(c, name='PLASMA_zgels')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgels_c
            integer(kind=c_int), value :: trans
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: descT
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
          end function PLASMA_zgels_c
      end interface

      interface
         function PLASMA_zgemm_c(transA,transB,M,N,K,alpha,A,LDA,B,LDB,beta,C,LDC) &
          & bind(c, name='PLASMA_zgemm')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgemm_c
            integer(kind=c_int), value :: transA
            integer(kind=c_int), value :: transB
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: K
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
            complex(kind=c_double_complex), value :: beta
            type(c_ptr), value :: C
            integer(kind=c_int), value :: LDC
          end function PLASMA_zgemm_c
      end interface

      interface
         function PLASMA_zgeqrf_c(M,N,A,LDA,descT) &
          & bind(c, name='PLASMA_zgeqrf')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgeqrf_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: descT
          end function PLASMA_zgeqrf_c
      end interface

      interface
         function PLASMA_zgeqrs_c(M,N,NRHS,A,LDA,descT,B,LDB) &
          & bind(c, name='PLASMA_zgeqrs')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgeqrs_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: descT
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
          end function PLASMA_zgeqrs_c
      end interface

      interface
         function PLASMA_zgesv_c(N,NRHS,A,LDA,IPIV,B,LDB) &
          & bind(c, name='PLASMA_zgesv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgesv_c
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
          end function PLASMA_zgesv_c
      end interface

      interface
         function PLASMA_zgesv_incpiv_c(N,NRHS,A,LDA,descL,IPIV,B,LDB) &
          & bind(c, name='PLASMA_zgesv_incpiv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgesv_incpiv_c
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: descL
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
          end function PLASMA_zgesv_incpiv_c
      end interface

      interface
         function PLASMA_zgesvd_c(jobu,jobvt,M,N,A,LDA,S,U,LDU,VT,LDVT,descT) &
          & bind(c, name='PLASMA_zgesvd')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgesvd_c
            integer(kind=c_int), value :: jobu
            integer(kind=c_int), value :: jobvt
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: S
            type(c_ptr), value :: U
            integer(kind=c_int), value :: LDU
            type(c_ptr), value :: VT
            integer(kind=c_int), value :: LDVT
            type(c_ptr), value :: descT
          end function PLASMA_zgesvd_c
      end interface

      interface
         function PLASMA_zgetmi_c(m,n,A,fin,mb,nb) &
          & bind(c, name='PLASMA_zgetmi')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgetmi_c
            integer(kind=c_int), value :: m
            integer(kind=c_int), value :: n
            type(c_ptr), value :: A
            integer(kind=c_int), value :: fin
            integer(kind=c_int), value :: mb
            integer(kind=c_int), value :: nb
          end function PLASMA_zgetmi_c
      end interface

      interface
         function PLASMA_zgetmi_Async_c(m,n,A,f_in,mb,inb,sequence,request) &
          & bind(c, name='PLASMA_zgetmi_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgetmi_Async_c
            integer(kind=c_int), value :: m
            integer(kind=c_int), value :: n
            type(c_ptr), value :: A
            integer(kind=c_int), value :: f_in
            integer(kind=c_int), value :: mb
            integer(kind=c_int), value :: inb
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zgetmi_Async_c
      end interface

      interface
         function PLASMA_zgetrf_c(M,N,A,LDA,IPIV) &
          & bind(c, name='PLASMA_zgetrf')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgetrf_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: IPIV
          end function PLASMA_zgetrf_c
      end interface

      interface
         function PLASMA_zgetrf_incpiv_c(M,N,A,LDA,descL,IPIV) &
          & bind(c, name='PLASMA_zgetrf_incpiv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgetrf_incpiv_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: descL
            type(c_ptr), value :: IPIV
          end function PLASMA_zgetrf_incpiv_c
      end interface

      interface
         function PLASMA_zgetri_c(N,A,LDA,IPIV) &
          & bind(c, name='PLASMA_zgetri')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgetri_c
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: IPIV
          end function PLASMA_zgetri_c
      end interface

      interface
         function PLASMA_zgetrs_c(trans,N,NRHS,A,LDA,IPIV,B,LDB) &
          & bind(c, name='PLASMA_zgetrs')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgetrs_c
            integer(kind=c_int), value :: trans
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
          end function PLASMA_zgetrs_c
      end interface

      interface
         function PLASMA_zgetrs_incpiv_c(trans,N,NRHS,A,LDA,descL,IPIV,B,LDB) &
          & bind(c, name='PLASMA_zgetrs_incpiv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgetrs_incpiv_c
            integer(kind=c_int), value :: trans
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: descL
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
          end function PLASMA_zgetrs_incpiv_c
      end interface

      interface
         function PLASMA_zheev_c(jobz,uplo,N,A,LDA,W,descT,Q,LDQ) &
          & bind(c, name='PLASMA_zheev')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zheev_c
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: W
            type(c_ptr), value :: descT
            type(c_ptr), value :: Q
            integer(kind=c_int), value :: LDQ
          end function PLASMA_zheev_c
      end interface

      interface
         function PLASMA_zheevd_c(jobz,uplo,N,A,LDA,W,descT,Q,LDQ) &
          & bind(c, name='PLASMA_zheevd')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zheevd_c
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: W
            type(c_ptr), value :: descT
            type(c_ptr), value :: Q
            integer(kind=c_int), value :: LDQ
          end function PLASMA_zheevd_c
      end interface

      interface
         function PLASMA_zhegst_c(itype,uplo,N,A,LDA,B,LDB) &
          & bind(c, name='PLASMA_zhegst')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zhegst_c
            integer(kind=c_int), value :: itype
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
          end function PLASMA_zhegst_c
      end interface

      interface
         function PLASMA_zhegv_c(itype,jobz,uplo,N,A,LDA,B,LDB,W,descT,Q,LDQ) &
          & bind(c, name='PLASMA_zhegv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zhegv_c
            integer(kind=c_int), value :: itype
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
            type(c_ptr), value :: W
            type(c_ptr), value :: descT
            type(c_ptr), value :: Q
            integer(kind=c_int), value :: LDQ
          end function PLASMA_zhegv_c
      end interface

      interface
         function PLASMA_zhegvd_c(itype,jobz,uplo,N,A,LDA,B,LDB,W,descT,Q,LDQ) &
          & bind(c, name='PLASMA_zhegvd')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zhegvd_c
            integer(kind=c_int), value :: itype
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
            type(c_ptr), value :: W
            type(c_ptr), value :: descT
            type(c_ptr), value :: Q
            integer(kind=c_int), value :: LDQ
          end function PLASMA_zhegvd_c
      end interface

#if defined(PRECISION_z) || defined(PRECISION_c)
      interface
         function PLASMA_zhemm_c(side,uplo,M,N,alpha,A,LDA,B,LDB,beta,C,LDC) &
          & bind(c, name='PLASMA_zhemm')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zhemm_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
            complex(kind=c_double_complex), value :: beta
            type(c_ptr), value :: C
            integer(kind=c_int), value :: LDC
          end function PLASMA_zhemm_c
      end interface
#endif

#if defined(PRECISION_z) || defined(PRECISION_c)
      interface
         function PLASMA_zher2k_c(uplo,trans,N,K,alpha,A,LDA,B,LDB,beta,C,LDC) &
          & bind(c, name='PLASMA_zher2k')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zher2k_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: K
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
            real(kind=c_double), value :: beta
            type(c_ptr), value :: C
            integer(kind=c_int), value :: LDC
          end function PLASMA_zher2k_c
      end interface
#endif

#if defined(PRECISION_z) || defined(PRECISION_c)
      interface
         function PLASMA_zherk_c(uplo,trans,N,K,alpha,A,LDA,beta,C,LDC) &
          & bind(c, name='PLASMA_zherk')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zherk_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: K
            real(kind=c_double), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            real(kind=c_double), value :: beta
            type(c_ptr), value :: C
            integer(kind=c_int), value :: LDC
          end function PLASMA_zherk_c
      end interface
#endif

      interface
         function PLASMA_zhetrd_c(jobz,uplo,N,A,LDA,D,E,descT,Q,LDQ) &
          & bind(c, name='PLASMA_zhetrd')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zhetrd_c
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: D
            type(c_ptr), value :: E
            type(c_ptr), value :: descT
            type(c_ptr), value :: Q
            integer(kind=c_int), value :: LDQ
          end function PLASMA_zhetrd_c
      end interface

      interface
         function PLASMA_zlacpy_c(uplo,M,N,A,LDA,B,LDB) &
          & bind(c, name='PLASMA_zlacpy')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zlacpy_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
          end function PLASMA_zlacpy_c
      end interface

      interface
         function PLASMA_zlange_c(norm,M,N,A,LDA) &
          & bind(c, name='PLASMA_zlange')
            use iso_c_binding
            implicit none
            real(kind=c_double) :: PLASMA_zlange_c
            integer(kind=c_int), value :: norm
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
          end function PLASMA_zlange_c
      end interface

#if defined(PRECISION_z) || defined(PRECISION_c)
      interface
         function PLASMA_zlanhe_c(norm,uplo,N,A,LDA) &
          & bind(c, name='PLASMA_zlanhe')
            use iso_c_binding
            implicit none
            real(kind=c_double) :: PLASMA_zlanhe_c
            integer(kind=c_int), value :: norm
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
          end function PLASMA_zlanhe_c
      end interface
#endif

      interface
         function PLASMA_zlansy_c(norm,uplo,N,A,LDA) &
          & bind(c, name='PLASMA_zlansy')
            use iso_c_binding
            implicit none
            real(kind=c_double) :: PLASMA_zlansy_c
            integer(kind=c_int), value :: norm
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
          end function PLASMA_zlansy_c
      end interface

      interface
         function PLASMA_zlaset_c(uplo,M,N,alpha,beta,A,LDA) &
          & bind(c, name='PLASMA_zlaset')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zlaset_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            complex(kind=c_double_complex), value :: alpha
            complex(kind=c_double_complex), value :: beta
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
          end function PLASMA_zlaset_c
      end interface

      interface
         function PLASMA_zlaswp_c(N,A,LDA,K1,K2,IPIV,INCX) &
          & bind(c, name='PLASMA_zlaswp')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zlaswp_c
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            integer(kind=c_int), value :: K1
            integer(kind=c_int), value :: K2
            type(c_ptr), value :: IPIV
            integer(kind=c_int), value :: INCX
          end function PLASMA_zlaswp_c
      end interface

      interface
         function PLASMA_zlaswpc_c(N,A,LDA,K1,K2,IPIV,INCX) &
          & bind(c, name='PLASMA_zlaswpc')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zlaswpc_c
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            integer(kind=c_int), value :: K1
            integer(kind=c_int), value :: K2
            type(c_ptr), value :: IPIV
            integer(kind=c_int), value :: INCX
          end function PLASMA_zlaswpc_c
      end interface

      interface
         function PLASMA_zlauum_c(uplo,N,A,LDA) &
          & bind(c, name='PLASMA_zlauum')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zlauum_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
          end function PLASMA_zlauum_c
      end interface

#if defined(PRECISION_z) || defined(PRECISION_c)
      interface
         function PLASMA_zplghe_c(bump,N,A,LDA,seed) &
          & bind(c, name='PLASMA_zplghe')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zplghe_c
            real(kind=c_double), value :: bump
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            integer(kind=c_long_long), value :: seed
          end function PLASMA_zplghe_c
      end interface
#endif

      interface
         function PLASMA_zplgsy_c(bump,N,A,LDA,seed) &
          & bind(c, name='PLASMA_zplgsy')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zplgsy_c
            complex(kind=c_double_complex), value :: bump
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            integer(kind=c_long_long), value :: seed
          end function PLASMA_zplgsy_c
      end interface

      interface
         function PLASMA_zplrnt_c(M,N,A,LDA,seed) &
          & bind(c, name='PLASMA_zplrnt')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zplrnt_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            integer(kind=c_long_long), value :: seed
          end function PLASMA_zplrnt_c
      end interface

      interface
         function PLASMA_zposv_c(uplo,N,NRHS,A,LDA,B,LDB) &
          & bind(c, name='PLASMA_zposv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zposv_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
          end function PLASMA_zposv_c
      end interface

      interface
         function PLASMA_zpotrf_c(uplo,N,A,LDA) &
          & bind(c, name='PLASMA_zpotrf')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zpotrf_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
          end function PLASMA_zpotrf_c
      end interface

      interface
         function PLASMA_zpotri_c(uplo,N,A,LDA) &
          & bind(c, name='PLASMA_zpotri')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zpotri_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
          end function PLASMA_zpotri_c
      end interface

      interface
         function PLASMA_zpotrs_c(uplo,N,NRHS,A,LDA,B,LDB) &
          & bind(c, name='PLASMA_zpotrs')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zpotrs_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
          end function PLASMA_zpotrs_c
      end interface

      interface
         function PLASMA_zsymm_c(side,uplo,M,N,alpha,A,LDA,B,LDB,beta,C,LDC) &
          & bind(c, name='PLASMA_zsymm')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zsymm_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
            complex(kind=c_double_complex), value :: beta
            type(c_ptr), value :: C
            integer(kind=c_int), value :: LDC
          end function PLASMA_zsymm_c
      end interface

      interface
         function PLASMA_zsyr2k_c(uplo,trans,N,K,alpha,A,LDA,B,LDB,beta,C,LDC) &
          & bind(c, name='PLASMA_zsyr2k')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zsyr2k_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: K
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
            complex(kind=c_double_complex), value :: beta
            type(c_ptr), value :: C
            integer(kind=c_int), value :: LDC
          end function PLASMA_zsyr2k_c
      end interface

      interface
         function PLASMA_zsyrk_c(uplo,trans,N,K,alpha,A,LDA,beta,C,LDC) &
          & bind(c, name='PLASMA_zsyrk')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zsyrk_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: K
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            complex(kind=c_double_complex), value :: beta
            type(c_ptr), value :: C
            integer(kind=c_int), value :: LDC
          end function PLASMA_zsyrk_c
      end interface

      interface
         function PLASMA_ztrmm_c(side,uplo,transA,diag,N,NRHS,alpha,A,LDA,B,LDB) &
          & bind(c, name='PLASMA_ztrmm')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_ztrmm_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: transA
            integer(kind=c_int), value :: diag
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
          end function PLASMA_ztrmm_c
      end interface

      interface
         function PLASMA_ztrsm_c(side,uplo,transA,diag,N,NRHS,alpha,A,LDA,B,LDB) &
          & bind(c, name='PLASMA_ztrsm')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_ztrsm_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: transA
            integer(kind=c_int), value :: diag
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
          end function PLASMA_ztrsm_c
      end interface

      interface
         function PLASMA_ztrsmpl_c(N,NRHS,A,LDA,descL,IPIV,B,LDB) &
          & bind(c, name='PLASMA_ztrsmpl')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_ztrsmpl_c
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: descL
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
          end function PLASMA_ztrsmpl_c
      end interface

      interface
         function PLASMA_ztrsmrv_c(side,uplo,transA,diag,N,NRHS,alpha,A,LDA,B,LDB) &
          & bind(c, name='PLASMA_ztrsmrv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_ztrsmrv_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: transA
            integer(kind=c_int), value :: diag
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
          end function PLASMA_ztrsmrv_c
      end interface

      interface
         function PLASMA_ztrtri_c(uplo,diag,N,A,LDA) &
          & bind(c, name='PLASMA_ztrtri')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_ztrtri_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: diag
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
          end function PLASMA_ztrtri_c
      end interface

      interface
         function PLASMA_zunglq_c(M,N,K,A,LDA,descT,B,LDB) &
          & bind(c, name='PLASMA_zunglq')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zunglq_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: K
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: descT
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
          end function PLASMA_zunglq_c
      end interface

      interface
         function PLASMA_zungqr_c(M,N,K,A,LDA,descT,B,LDB) &
          & bind(c, name='PLASMA_zungqr')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zungqr_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: K
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: descT
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
          end function PLASMA_zungqr_c
      end interface

      interface
         function PLASMA_zunmlq_c(side,trans,M,N,K,A,LDA,descT,B,LDB) &
          & bind(c, name='PLASMA_zunmlq')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zunmlq_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: trans
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: K
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: descT
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
          end function PLASMA_zunmlq_c
      end interface

      interface
         function PLASMA_zunmqr_c(side,trans,M,N,K,A,LDA,descT,B,LDB) &
          & bind(c, name='PLASMA_zunmqr')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zunmqr_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: trans
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: K
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: descT
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
          end function PLASMA_zunmqr_c
      end interface

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  FORTRAN API - math functions (native interface)
    !
      interface
         function PLASMA_zgebrd_Tile_c(A,D,E,T) &
          & bind(c, name='PLASMA_zgebrd_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgebrd_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: D
            type(c_ptr), value :: E
            type(c_ptr), value :: T
          end function PLASMA_zgebrd_Tile_c
      end interface

      interface
         function PLASMA_zgelqf_Tile_c(A,T) &
          & bind(c, name='PLASMA_zgelqf_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgelqf_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: T
          end function PLASMA_zgelqf_Tile_c
      end interface

      interface
         function PLASMA_zgelqs_Tile_c(A,T,B) &
          & bind(c, name='PLASMA_zgelqs_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgelqs_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
          end function PLASMA_zgelqs_Tile_c
      end interface

      interface
         function PLASMA_zgels_Tile_c(trans,A,T,B) &
          & bind(c, name='PLASMA_zgels_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgels_Tile_c
            integer(kind=c_int), value :: trans
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
          end function PLASMA_zgels_Tile_c
      end interface

      interface
         function PLASMA_zgemm_Tile_c(transA,transB,alpha,A,B,beta,C) &
          & bind(c, name='PLASMA_zgemm_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgemm_Tile_c
            integer(kind=c_int), value :: transA
            integer(kind=c_int), value :: transB
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            complex(kind=c_double_complex), value :: beta
            type(c_ptr), value :: C
          end function PLASMA_zgemm_Tile_c
      end interface

      interface
         function PLASMA_zgeqrf_Tile_c(A,T) &
          & bind(c, name='PLASMA_zgeqrf_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgeqrf_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: T
          end function PLASMA_zgeqrf_Tile_c
      end interface

      interface
         function PLASMA_zgeqrs_Tile_c(A,T,B) &
          & bind(c, name='PLASMA_zgeqrs_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgeqrs_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
          end function PLASMA_zgeqrs_Tile_c
      end interface

      interface
         function PLASMA_zgesv_Tile_c(A,IPIV,B) &
          & bind(c, name='PLASMA_zgesv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgesv_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
          end function PLASMA_zgesv_Tile_c
      end interface

      interface
         function PLASMA_zgesv_incpiv_Tile_c(A,L,IPIV,B) &
          & bind(c, name='PLASMA_zgesv_incpiv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgesv_incpiv_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: L
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
          end function PLASMA_zgesv_incpiv_Tile_c
      end interface

      interface
         function PLASMA_zgesvd_Tile_c(jobu,jobvt,A,S,U,VT,T) &
          & bind(c, name='PLASMA_zgesvd_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgesvd_Tile_c
            integer(kind=c_int), value :: jobu
            integer(kind=c_int), value :: jobvt
            type(c_ptr), value :: A
            type(c_ptr), value :: S
            type(c_ptr), value :: U
            type(c_ptr), value :: VT
            type(c_ptr), value :: T
          end function PLASMA_zgesvd_Tile_c
      end interface

      interface
         function PLASMA_zgetrf_Tile_c(A,IPIV) &
          & bind(c, name='PLASMA_zgetrf_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgetrf_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: IPIV
          end function PLASMA_zgetrf_Tile_c
      end interface

      interface
         function PLASMA_zgetrf_incpiv_Tile_c(A,L,IPIV) &
          & bind(c, name='PLASMA_zgetrf_incpiv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgetrf_incpiv_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: L
            type(c_ptr), value :: IPIV
          end function PLASMA_zgetrf_incpiv_Tile_c
      end interface

      interface
         function PLASMA_zgetri_Tile_c(A,IPIV) &
          & bind(c, name='PLASMA_zgetri_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgetri_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: IPIV
          end function PLASMA_zgetri_Tile_c
      end interface

      interface
         function PLASMA_zgetrs_Tile_c(trans,A,IPIV,B) &
          & bind(c, name='PLASMA_zgetrs_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgetrs_Tile_c
            integer(kind=c_int), value :: trans
            type(c_ptr), value :: A
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
          end function PLASMA_zgetrs_Tile_c
      end interface

      interface
         function PLASMA_zgetrs_incpiv_Tile_c(A,L,IPIV,B) &
          & bind(c, name='PLASMA_zgetrs_incpiv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgetrs_incpiv_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: L
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
          end function PLASMA_zgetrs_incpiv_Tile_c
      end interface

      interface
         function PLASMA_zheev_Tile_c(jobz,uplo,A,W,T,Q,LDQ) &
          & bind(c, name='PLASMA_zheev_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zheev_Tile_c
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: W
            type(c_ptr), value :: T
            type(c_ptr), value :: Q
            integer(kind=c_int), value :: LDQ
          end function PLASMA_zheev_Tile_c
      end interface

      interface
         function PLASMA_zheevd_Tile_c(jobz,uplo,A,W,T,Q,LDQ) &
          & bind(c, name='PLASMA_zheevd_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zheevd_Tile_c
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: W
            type(c_ptr), value :: T
            type(c_ptr), value :: Q
            integer(kind=c_int), value :: LDQ
          end function PLASMA_zheevd_Tile_c
      end interface

      interface
         function PLASMA_zhegst_Tile_c(itype,uplo,A,B) &
          & bind(c, name='PLASMA_zhegst_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zhegst_Tile_c
            integer(kind=c_int), value :: itype
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
          end function PLASMA_zhegst_Tile_c
      end interface

      interface
         function PLASMA_zhegv_Tile_c(itype,jobz,uplo,A,B,W,T,Q) &
          & bind(c, name='PLASMA_zhegv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zhegv_Tile_c
            integer(kind=c_int), value :: itype
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            type(c_ptr), value :: W
            type(c_ptr), value :: T
            type(c_ptr), value :: Q
          end function PLASMA_zhegv_Tile_c
      end interface

      interface
         function PLASMA_zhegvd_Tile_c(itype,jobz,uplo,A,B,W,T,Q) &
          & bind(c, name='PLASMA_zhegvd_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zhegvd_Tile_c
            integer(kind=c_int), value :: itype
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            type(c_ptr), value :: W
            type(c_ptr), value :: T
            type(c_ptr), value :: Q
          end function PLASMA_zhegvd_Tile_c
      end interface

#if defined(PRECISION_z) || defined(PRECISION_c)
      interface
         function PLASMA_zhemm_Tile_c(side,uplo,alpha,A,B,beta,C) &
          & bind(c, name='PLASMA_zhemm_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zhemm_Tile_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            complex(kind=c_double_complex), value :: beta
            type(c_ptr), value :: C
          end function PLASMA_zhemm_Tile_c
      end interface
#endif

#if defined(PRECISION_z) || defined(PRECISION_c)
      interface
         function PLASMA_zher2k_Tile_c(uplo,trans,alpha,A,B,beta,C) &
          & bind(c, name='PLASMA_zher2k_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zher2k_Tile_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            real(kind=c_double), value :: beta
            type(c_ptr), value :: C
          end function PLASMA_zher2k_Tile_c
      end interface
#endif

#if defined(PRECISION_z) || defined(PRECISION_c)
      interface
         function PLASMA_zherk_Tile_c(uplo,trans,alpha,A,beta,C) &
          & bind(c, name='PLASMA_zherk_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zherk_Tile_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            real(kind=c_double), value :: alpha
            type(c_ptr), value :: A
            real(kind=c_double), value :: beta
            type(c_ptr), value :: C
          end function PLASMA_zherk_Tile_c
      end interface
#endif

      interface
         function PLASMA_zhetrd_Tile_c(jobz,uplo,A,D,E,T,Q,LDQ) &
          & bind(c, name='PLASMA_zhetrd_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zhetrd_Tile_c
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: D
            type(c_ptr), value :: E
            type(c_ptr), value :: T
            type(c_ptr), value :: Q
            integer(kind=c_int), value :: LDQ
          end function PLASMA_zhetrd_Tile_c
      end interface

      interface
         function PLASMA_zlacpy_Tile_c(uplo,A,B) &
          & bind(c, name='PLASMA_zlacpy_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zlacpy_Tile_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
          end function PLASMA_zlacpy_Tile_c
      end interface

      interface
         function PLASMA_zlange_Tile_c(norm,A) &
          & bind(c, name='PLASMA_zlange_Tile')
            use iso_c_binding
            implicit none
            real(kind=c_double) :: PLASMA_zlange_Tile_c
            integer(kind=c_int), value :: norm
            type(c_ptr), value :: A
          end function PLASMA_zlange_Tile_c
      end interface

#if defined(PRECISION_z) || defined(PRECISION_c)
      interface
         function PLASMA_zlanhe_Tile_c(norm,uplo,A) &
          & bind(c, name='PLASMA_zlanhe_Tile')
            use iso_c_binding
            implicit none
            real(kind=c_double) :: PLASMA_zlanhe_Tile_c
            integer(kind=c_int), value :: norm
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
          end function PLASMA_zlanhe_Tile_c
      end interface
#endif

      interface
         function PLASMA_zlansy_Tile_c(norm,uplo,A) &
          & bind(c, name='PLASMA_zlansy_Tile')
            use iso_c_binding
            implicit none
            real(kind=c_double) :: PLASMA_zlansy_Tile_c
            integer(kind=c_int), value :: norm
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
          end function PLASMA_zlansy_Tile_c
      end interface

      interface
         function PLASMA_zlaset_Tile_c(uplo,alpha,beta,A) &
          & bind(c, name='PLASMA_zlaset_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zlaset_Tile_c
            integer(kind=c_int), value :: uplo
            complex(kind=c_double_complex), value :: alpha
            complex(kind=c_double_complex), value :: beta
            type(c_ptr), value :: A
          end function PLASMA_zlaset_Tile_c
      end interface

      interface
         function PLASMA_zlaswp_Tile_c(A,K1,K2,IPIV,INCX) &
          & bind(c, name='PLASMA_zlaswp_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zlaswp_Tile_c
            type(c_ptr), value :: A
            integer(kind=c_int), value :: K1
            integer(kind=c_int), value :: K2
            type(c_ptr), value :: IPIV
            integer(kind=c_int), value :: INCX
          end function PLASMA_zlaswp_Tile_c
      end interface

      interface
         function PLASMA_zlaswpc_Tile_c(A,K1,K2,IPIV,INCX) &
          & bind(c, name='PLASMA_zlaswpc_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zlaswpc_Tile_c
            type(c_ptr), value :: A
            integer(kind=c_int), value :: K1
            integer(kind=c_int), value :: K2
            type(c_ptr), value :: IPIV
            integer(kind=c_int), value :: INCX
          end function PLASMA_zlaswpc_Tile_c
      end interface

      interface
         function PLASMA_zlauum_Tile_c(uplo,A) &
          & bind(c, name='PLASMA_zlauum_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zlauum_Tile_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
          end function PLASMA_zlauum_Tile_c
      end interface

#if defined(PRECISION_z) || defined(PRECISION_c)
      interface
         function PLASMA_zplghe_Tile_c(bump,A,seed) &
          & bind(c, name='PLASMA_zplghe_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zplghe_Tile_c
            real(kind=c_double), value :: bump
            type(c_ptr), value :: A
            integer(kind=c_long_long), value :: seed
          end function PLASMA_zplghe_Tile_c
      end interface
#endif

      interface
         function PLASMA_zplgsy_Tile_c(bump,A,seed) &
          & bind(c, name='PLASMA_zplgsy_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zplgsy_Tile_c
            complex(kind=c_double_complex), value :: bump
            type(c_ptr), value :: A
            integer(kind=c_long_long), value :: seed
          end function PLASMA_zplgsy_Tile_c
      end interface

      interface
         function PLASMA_zplrnt_Tile_c(A,seed) &
          & bind(c, name='PLASMA_zplrnt_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zplrnt_Tile_c
            type(c_ptr), value :: A
            integer(kind=c_long_long), value :: seed
          end function PLASMA_zplrnt_Tile_c
      end interface

      interface
         function PLASMA_zposv_Tile_c(uplo,A,B) &
          & bind(c, name='PLASMA_zposv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zposv_Tile_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
          end function PLASMA_zposv_Tile_c
      end interface

      interface
         function PLASMA_zpotrf_Tile_c(uplo,A) &
          & bind(c, name='PLASMA_zpotrf_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zpotrf_Tile_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
          end function PLASMA_zpotrf_Tile_c
      end interface

      interface
         function PLASMA_zpotri_Tile_c(uplo,A) &
          & bind(c, name='PLASMA_zpotri_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zpotri_Tile_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
          end function PLASMA_zpotri_Tile_c
      end interface

      interface
         function PLASMA_zpotrs_Tile_c(uplo,A,B) &
          & bind(c, name='PLASMA_zpotrs_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zpotrs_Tile_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
          end function PLASMA_zpotrs_Tile_c
      end interface

      interface
         function PLASMA_zsymm_Tile_c(side,uplo,alpha,A,B,beta,C) &
          & bind(c, name='PLASMA_zsymm_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zsymm_Tile_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            complex(kind=c_double_complex), value :: beta
            type(c_ptr), value :: C
          end function PLASMA_zsymm_Tile_c
      end interface

      interface
         function PLASMA_zsyr2k_Tile_c(uplo,trans,alpha,A,B,beta,C) &
          & bind(c, name='PLASMA_zsyr2k_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zsyr2k_Tile_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            complex(kind=c_double_complex), value :: beta
            type(c_ptr), value :: C
          end function PLASMA_zsyr2k_Tile_c
      end interface

      interface
         function PLASMA_zsyrk_Tile_c(uplo,trans,alpha,A,beta,C) &
          & bind(c, name='PLASMA_zsyrk_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zsyrk_Tile_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            complex(kind=c_double_complex), value :: beta
            type(c_ptr), value :: C
          end function PLASMA_zsyrk_Tile_c
      end interface

      interface
         function PLASMA_ztrmm_Tile_c(side,uplo,transA,diag,alpha,A,B) &
          & bind(c, name='PLASMA_ztrmm_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_ztrmm_Tile_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: transA
            integer(kind=c_int), value :: diag
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
          end function PLASMA_ztrmm_Tile_c
      end interface

      interface
         function PLASMA_ztrsm_Tile_c(side,uplo,transA,diag,alpha,A,B) &
          & bind(c, name='PLASMA_ztrsm_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_ztrsm_Tile_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: transA
            integer(kind=c_int), value :: diag
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
          end function PLASMA_ztrsm_Tile_c
      end interface

      interface
         function PLASMA_ztrsmpl_Tile_c(A,L,IPIV,B) &
          & bind(c, name='PLASMA_ztrsmpl_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_ztrsmpl_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: L
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
          end function PLASMA_ztrsmpl_Tile_c
      end interface

      interface
         function PLASMA_ztrsmrv_Tile_c(side,uplo,transA,diag,alpha,A,B) &
          & bind(c, name='PLASMA_ztrsmrv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_ztrsmrv_Tile_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: transA
            integer(kind=c_int), value :: diag
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
          end function PLASMA_ztrsmrv_Tile_c
      end interface

      interface
         function PLASMA_ztrtri_Tile_c(uplo,diag,A) &
          & bind(c, name='PLASMA_ztrtri_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_ztrtri_Tile_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: diag
            type(c_ptr), value :: A
          end function PLASMA_ztrtri_Tile_c
      end interface

      interface
         function PLASMA_zunglq_Tile_c(A,T,B) &
          & bind(c, name='PLASMA_zunglq_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zunglq_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
          end function PLASMA_zunglq_Tile_c
      end interface

      interface
         function PLASMA_zungqr_Tile_c(A,T,B) &
          & bind(c, name='PLASMA_zungqr_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zungqr_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
          end function PLASMA_zungqr_Tile_c
      end interface

      interface
         function PLASMA_zunmlq_Tile_c(side,trans,A,T,B) &
          & bind(c, name='PLASMA_zunmlq_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zunmlq_Tile_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: trans
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
          end function PLASMA_zunmlq_Tile_c
      end interface

      interface
         function PLASMA_zunmqr_Tile_c(side,trans,A,T,B) &
          & bind(c, name='PLASMA_zunmqr_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zunmqr_Tile_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: trans
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
          end function PLASMA_zunmqr_Tile_c
      end interface

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  FORTRAN API - math functions (asynchronous interface)
    !
      interface
         function PLASMA_zgebrd_Tile_Async_c(A,D,E,T,sequence,request) &
          & bind(c, name='PLASMA_zgebrd_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgebrd_Tile_Async_c
            type(c_ptr), value :: A
            type(c_ptr), value :: D
            type(c_ptr), value :: E
            type(c_ptr), value :: T
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zgebrd_Tile_Async_c
      end interface

      interface
         function PLASMA_zgelqf_Tile_Async_c(A,T,sequence,request) &
          & bind(c, name='PLASMA_zgelqf_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgelqf_Tile_Async_c
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zgelqf_Tile_Async_c
      end interface

      interface
         function PLASMA_zgelqs_Tile_Async_c(A,T,B,sequence,request) &
          & bind(c, name='PLASMA_zgelqs_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgelqs_Tile_Async_c
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zgelqs_Tile_Async_c
      end interface

      interface
         function PLASMA_zgels_Tile_Async_c(trans,A,T,B,sequence,request) &
          & bind(c, name='PLASMA_zgels_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgels_Tile_Async_c
            integer(kind=c_int), value :: trans
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zgels_Tile_Async_c
      end interface

      interface
         function PLASMA_zgemm_Tile_Async_c(transA,transB,alpha,A,B,beta,C,sequence,request) &
          & bind(c, name='PLASMA_zgemm_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgemm_Tile_Async_c
            integer(kind=c_int), value :: transA
            integer(kind=c_int), value :: transB
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            complex(kind=c_double_complex), value :: beta
            type(c_ptr), value :: C
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zgemm_Tile_Async_c
      end interface

      interface
         function PLASMA_zgeqrf_Tile_Async_c(A,T,sequence,request) &
          & bind(c, name='PLASMA_zgeqrf_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgeqrf_Tile_Async_c
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zgeqrf_Tile_Async_c
      end interface

      interface
         function PLASMA_zgeqrs_Tile_Async_c(A,T,B,sequence,request) &
          & bind(c, name='PLASMA_zgeqrs_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgeqrs_Tile_Async_c
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zgeqrs_Tile_Async_c
      end interface

      interface
         function PLASMA_zgesv_Tile_Async_c(A,IPIV,B,sequence,request) &
          & bind(c, name='PLASMA_zgesv_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgesv_Tile_Async_c
            type(c_ptr), value :: A
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zgesv_Tile_Async_c
      end interface

      interface
         function PLASMA_zgesv_incpiv_Tile_Async_c(A,L,IPIV,B,sequence,request) &
          & bind(c, name='PLASMA_zgesv_incpiv_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgesv_incpiv_Tile_Async_c
            type(c_ptr), value :: A
            type(c_ptr), value :: L
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zgesv_incpiv_Tile_Async_c
      end interface

      interface
         function PLASMA_zgesvd_Tile_Async_c(jobu,jobvt,A,S,U,VT,T,sequence,request) &
          & bind(c, name='PLASMA_zgesvd_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgesvd_Tile_Async_c
            integer(kind=c_int), value :: jobu
            integer(kind=c_int), value :: jobvt
            type(c_ptr), value :: A
            type(c_ptr), value :: S
            type(c_ptr), value :: U
            type(c_ptr), value :: VT
            type(c_ptr), value :: T
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zgesvd_Tile_Async_c
      end interface

      interface
         function PLASMA_zgetrf_Tile_Async_c(A,IPIV,sequence,request) &
          & bind(c, name='PLASMA_zgetrf_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgetrf_Tile_Async_c
            type(c_ptr), value :: A
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zgetrf_Tile_Async_c
      end interface

      interface
         function PLASMA_zgetrf_incpiv_Tile_Async_c(A,L,IPIV,sequence,request) &
          & bind(c, name='PLASMA_zgetrf_incpiv_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgetrf_incpiv_Tile_Async_c
            type(c_ptr), value :: A
            type(c_ptr), value :: L
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zgetrf_incpiv_Tile_Async_c
      end interface

      interface
         function PLASMA_zgetri_Tile_Async_c(A,IPIV,W,sequence,request) &
          & bind(c, name='PLASMA_zgetri_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgetri_Tile_Async_c
            type(c_ptr), value :: A
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: W
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zgetri_Tile_Async_c
      end interface

      interface
         function PLASMA_zgetrs_Tile_Async_c(trans,A,IPIV,B,sequence,request) &
          & bind(c, name='PLASMA_zgetrs_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgetrs_Tile_Async_c
            integer(kind=c_int), value :: trans
            type(c_ptr), value :: A
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zgetrs_Tile_Async_c
      end interface

      interface
         function PLASMA_zgetrs_incpiv_Tile_Async_c(A,L,IPIV,B,sequence,request) &
          & bind(c, name='PLASMA_zgetrs_incpiv_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zgetrs_incpiv_Tile_Async_c
            type(c_ptr), value :: A
            type(c_ptr), value :: L
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zgetrs_incpiv_Tile_Async_c
      end interface

      interface
         function PLASMA_zheev_Tile_Async_c(jobz,uplo,A,W,T,Q,LDQ,sequence,request) &
          & bind(c, name='PLASMA_zheev_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zheev_Tile_Async_c
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: W
            type(c_ptr), value :: T
            type(c_ptr), value :: Q
            integer(kind=c_int), value :: LDQ
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zheev_Tile_Async_c
      end interface

      interface
         function PLASMA_zheevd_Tile_Async_c(jobz,uplo,A,W,T,Q,LDQ,sequence,request) &
          & bind(c, name='PLASMA_zheevd_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zheevd_Tile_Async_c
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: W
            type(c_ptr), value :: T
            type(c_ptr), value :: Q
            integer(kind=c_int), value :: LDQ
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zheevd_Tile_Async_c
      end interface

      interface
         function PLASMA_zhegst_Tile_Async_c(itype,uplo,A,B,sequence,request) &
          & bind(c, name='PLASMA_zhegst_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zhegst_Tile_Async_c
            integer(kind=c_int), value :: itype
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zhegst_Tile_Async_c
      end interface

      interface
         function PLASMA_zhegv_Tile_Async_c(itype,jobz,uplo,A,B,W,T,Q,sequence,request) &
          & bind(c, name='PLASMA_zhegv_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zhegv_Tile_Async_c
            integer(kind=c_int), value :: itype
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            type(c_ptr), value :: W
            type(c_ptr), value :: T
            type(c_ptr), value :: Q
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zhegv_Tile_Async_c
      end interface

      interface
         function PLASMA_zhegvd_Tile_Async_c(itype,jobz,uplo,A,B,W,T,Q,sequence,request) &
          & bind(c, name='PLASMA_zhegvd_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zhegvd_Tile_Async_c
            integer(kind=c_int), value :: itype
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            type(c_ptr), value :: W
            type(c_ptr), value :: T
            type(c_ptr), value :: Q
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zhegvd_Tile_Async_c
      end interface

#if defined(PRECISION_z) || defined(PRECISION_c)
      interface
         function PLASMA_zhemm_Tile_Async_c(side,uplo,alpha,A,B,beta,C,sequence,request) &
          & bind(c, name='PLASMA_zhemm_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zhemm_Tile_Async_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            complex(kind=c_double_complex), value :: beta
            type(c_ptr), value :: C
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zhemm_Tile_Async_c
      end interface
#endif

#if defined(PRECISION_z) || defined(PRECISION_c)
      interface
         function PLASMA_zher2k_Tile_Async_c(uplo,trans,alpha,A,B,beta,C,sequence,request) &
          & bind(c, name='PLASMA_zher2k_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zher2k_Tile_Async_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            real(kind=c_double), value :: beta
            type(c_ptr), value :: C
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zher2k_Tile_Async_c
      end interface
#endif

#if defined(PRECISION_z) || defined(PRECISION_c)
      interface
         function PLASMA_zherk_Tile_Async_c(uplo,trans,alpha,A,beta,C,sequence,request) &
          & bind(c, name='PLASMA_zherk_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zherk_Tile_Async_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            real(kind=c_double), value :: alpha
            type(c_ptr), value :: A
            real(kind=c_double), value :: beta
            type(c_ptr), value :: C
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zherk_Tile_Async_c
      end interface
#endif

      interface
         function PLASMA_zhetrd_Tile_Async_c(jobz,uplo,A,D,E,T,Q,LDQ,sequence,request) &
          & bind(c, name='PLASMA_zhetrd_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zhetrd_Tile_Async_c
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: D
            type(c_ptr), value :: E
            type(c_ptr), value :: T
            type(c_ptr), value :: Q
            integer(kind=c_int), value :: LDQ
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zhetrd_Tile_Async_c
      end interface

      interface
         function PLASMA_zlacpy_Tile_Async_c(uplo,A,B,sequence,request) &
          & bind(c, name='PLASMA_zlacpy_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zlacpy_Tile_Async_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zlacpy_Tile_Async_c
      end interface

      interface
         function PLASMA_zlange_Tile_Async_c(norm,A,value,sequence,request) &
          & bind(c, name='PLASMA_zlange_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zlange_Tile_Async_c
            integer(kind=c_int), value :: norm
            type(c_ptr), value :: A
            type(c_ptr), value :: value
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zlange_Tile_Async_c
      end interface

#if defined(PRECISION_z) || defined(PRECISION_c)
      interface
         function PLASMA_zlanhe_Tile_Async_c(norm,uplo,A,value,sequence,request) &
          & bind(c, name='PLASMA_zlanhe_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zlanhe_Tile_Async_c
            integer(kind=c_int), value :: norm
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: value
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zlanhe_Tile_Async_c
      end interface
#endif

      interface
         function PLASMA_zlansy_Tile_Async_c(norm,uplo,A,value,sequence,request) &
          & bind(c, name='PLASMA_zlansy_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zlansy_Tile_Async_c
            integer(kind=c_int), value :: norm
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: value
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zlansy_Tile_Async_c
      end interface

      interface
         function PLASMA_zlaset_Tile_Async_c(uplo,alpha,beta,A,sequence,request) &
          & bind(c, name='PLASMA_zlaset_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zlaset_Tile_Async_c
            integer(kind=c_int), value :: uplo
            complex(kind=c_double_complex), value :: alpha
            complex(kind=c_double_complex), value :: beta
            type(c_ptr), value :: A
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zlaset_Tile_Async_c
      end interface

      interface
         function PLASMA_zlaswp_Tile_Async_c(A,K1,K2,IPIV,INCX,sequence,request) &
          & bind(c, name='PLASMA_zlaswp_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zlaswp_Tile_Async_c
            type(c_ptr), value :: A
            integer(kind=c_int), value :: K1
            integer(kind=c_int), value :: K2
            type(c_ptr), value :: IPIV
            integer(kind=c_int), value :: INCX
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zlaswp_Tile_Async_c
      end interface

      interface
         function PLASMA_zlaswpc_Tile_Async_c(A,K1,K2,IPIV,INCX,sequence,request) &
          & bind(c, name='PLASMA_zlaswpc_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zlaswpc_Tile_Async_c
            type(c_ptr), value :: A
            integer(kind=c_int), value :: K1
            integer(kind=c_int), value :: K2
            type(c_ptr), value :: IPIV
            integer(kind=c_int), value :: INCX
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zlaswpc_Tile_Async_c
      end interface

      interface
         function PLASMA_zlauum_Tile_Async_c(uplo,A,sequence,request) &
          & bind(c, name='PLASMA_zlauum_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zlauum_Tile_Async_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zlauum_Tile_Async_c
      end interface

#if defined(PRECISION_z) || defined(PRECISION_c)
      interface
         function PLASMA_zplghe_Tile_Async_c(bump,A,seed,sequence,request) &
          & bind(c, name='PLASMA_zplghe_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zplghe_Tile_Async_c
            real(kind=c_double), value :: bump
            type(c_ptr), value :: A
            integer(kind=c_long_long), value :: seed
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zplghe_Tile_Async_c
      end interface
#endif

      interface
         function PLASMA_zplgsy_Tile_Async_c(bump,A,seed,sequence,request) &
          & bind(c, name='PLASMA_zplgsy_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zplgsy_Tile_Async_c
            complex(kind=c_double_complex), value :: bump
            type(c_ptr), value :: A
            integer(kind=c_long_long), value :: seed
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zplgsy_Tile_Async_c
      end interface

      interface
         function PLASMA_zplrnt_Tile_Async_c(A,seed,sequence,request) &
          & bind(c, name='PLASMA_zplrnt_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zplrnt_Tile_Async_c
            type(c_ptr), value :: A
            integer(kind=c_long_long), value :: seed
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zplrnt_Tile_Async_c
      end interface

      interface
         function PLASMA_zposv_Tile_Async_c(uplo,A,B,sequence,request) &
          & bind(c, name='PLASMA_zposv_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zposv_Tile_Async_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zposv_Tile_Async_c
      end interface

      interface
         function PLASMA_zpotrf_Tile_Async_c(uplo,A,sequence,request) &
          & bind(c, name='PLASMA_zpotrf_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zpotrf_Tile_Async_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zpotrf_Tile_Async_c
      end interface

      interface
         function PLASMA_zpotri_Tile_Async_c(uplo,A,sequence,request) &
          & bind(c, name='PLASMA_zpotri_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zpotri_Tile_Async_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zpotri_Tile_Async_c
      end interface

      interface
         function PLASMA_zpotrs_Tile_Async_c(uplo,A,B,sequence,request) &
          & bind(c, name='PLASMA_zpotrs_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zpotrs_Tile_Async_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zpotrs_Tile_Async_c
      end interface

      interface
         function PLASMA_zsymm_Tile_Async_c(side,uplo,alpha,A,B,beta,C,sequence,request) &
          & bind(c, name='PLASMA_zsymm_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zsymm_Tile_Async_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            complex(kind=c_double_complex), value :: beta
            type(c_ptr), value :: C
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zsymm_Tile_Async_c
      end interface

      interface
         function PLASMA_zsyr2k_Tile_Async_c(uplo,trans,alpha,A,B,beta,C,sequence,request) &
          & bind(c, name='PLASMA_zsyr2k_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zsyr2k_Tile_Async_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            complex(kind=c_double_complex), value :: beta
            type(c_ptr), value :: C
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zsyr2k_Tile_Async_c
      end interface

      interface
         function PLASMA_zsyrk_Tile_Async_c(uplo,trans,alpha,A,beta,C,sequence,request) &
          & bind(c, name='PLASMA_zsyrk_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zsyrk_Tile_Async_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            complex(kind=c_double_complex), value :: beta
            type(c_ptr), value :: C
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zsyrk_Tile_Async_c
      end interface

      interface
         function PLASMA_ztrmm_Tile_Async_c(side,uplo,transA,diag,alpha,A,B,sequence,request) &
          & bind(c, name='PLASMA_ztrmm_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_ztrmm_Tile_Async_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: transA
            integer(kind=c_int), value :: diag
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_ztrmm_Tile_Async_c
      end interface

      interface
         function PLASMA_ztrsm_Tile_Async_c(side,uplo,transA,diag,alpha,A,B,sequence,request) &
          & bind(c, name='PLASMA_ztrsm_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_ztrsm_Tile_Async_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: transA
            integer(kind=c_int), value :: diag
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_ztrsm_Tile_Async_c
      end interface

      interface
         function PLASMA_ztrsmpl_Tile_Async_c(A,L,IPIV,B,sequence,request) &
          & bind(c, name='PLASMA_ztrsmpl_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_ztrsmpl_Tile_Async_c
            type(c_ptr), value :: A
            type(c_ptr), value :: L
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_ztrsmpl_Tile_Async_c
      end interface

      interface
         function PLASMA_ztrsmrv_Tile_Async_c(side,uplo,transA,diag,alpha,A,B,sequence,request) &
          & bind(c, name='PLASMA_ztrsmrv_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_ztrsmrv_Tile_Async_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: transA
            integer(kind=c_int), value :: diag
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_ztrsmrv_Tile_Async_c
      end interface

      interface
         function PLASMA_ztrtri_Tile_Async_c(uplo,diag,A,sequence,request) &
          & bind(c, name='PLASMA_ztrtri_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_ztrtri_Tile_Async_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: diag
            type(c_ptr), value :: A
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_ztrtri_Tile_Async_c
      end interface

      interface
         function PLASMA_zunglq_Tile_Async_c(A,T,B,sequence,request) &
          & bind(c, name='PLASMA_zunglq_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zunglq_Tile_Async_c
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zunglq_Tile_Async_c
      end interface

      interface
         function PLASMA_zungqr_Tile_Async_c(A,T,B,sequence,request) &
          & bind(c, name='PLASMA_zungqr_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zungqr_Tile_Async_c
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zungqr_Tile_Async_c
      end interface

      interface
         function PLASMA_zunmlq_Tile_Async_c(side,trans,A,T,B,sequence,request) &
          & bind(c, name='PLASMA_zunmlq_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zunmlq_Tile_Async_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: trans
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zunmlq_Tile_Async_c
      end interface

      interface
         function PLASMA_zunmqr_Tile_Async_c(side,trans,A,T,B,sequence,request) &
          & bind(c, name='PLASMA_zunmqr_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_zunmqr_Tile_Async_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: trans
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function PLASMA_zunmqr_Tile_Async_c
      end interface

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  FORTRAN API - workspace allocation
    !
      interface
         function PLASMA_Alloc_Workspace_zgebrd_c(M,N,descT) &
          & bind(c, name='PLASMA_Alloc_Workspace_zgebrd')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_zgebrd_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
          end function PLASMA_Alloc_Workspace_zgebrd_c
      end interface

      interface
         function PLASMA_Alloc_Workspace_zgeev_c(N,descT) &
          & bind(c, name='PLASMA_Alloc_Workspace_zgeev')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_zgeev_c
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
          end function PLASMA_Alloc_Workspace_zgeev_c
      end interface

      interface
         function PLASMA_Alloc_Workspace_zgehrd_c(N,descT) &
          & bind(c, name='PLASMA_Alloc_Workspace_zgehrd')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_zgehrd_c
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
          end function PLASMA_Alloc_Workspace_zgehrd_c
      end interface

      interface
         function PLASMA_Alloc_Workspace_zgelqf_c(M,N,T) &
          & bind(c, name='PLASMA_Alloc_Workspace_zgelqf')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_zgelqf_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: T ! T is **, so pass by reference
          end function PLASMA_Alloc_Workspace_zgelqf_c
      end interface

      interface
         function PLASMA_Alloc_Workspace_zgelqf_Tile_c(M,N,descT) &
          & bind(c, name='PLASMA_Alloc_Workspace_zgelqf_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_zgelqf_Tile_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
          end function PLASMA_Alloc_Workspace_zgelqf_Tile_c
      end interface

      interface
         function PLASMA_Alloc_Workspace_zgels_c(M,N,T) &
          & bind(c, name='PLASMA_Alloc_Workspace_zgels')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_zgels_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: T ! T is **, so pass by reference
          end function PLASMA_Alloc_Workspace_zgels_c
      end interface

      interface
         function PLASMA_Alloc_Workspace_zgels_Tile_c(M,N,descT) &
          & bind(c, name='PLASMA_Alloc_Workspace_zgels_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_zgels_Tile_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
          end function PLASMA_Alloc_Workspace_zgels_Tile_c
      end interface

      interface
         function PLASMA_Alloc_Workspace_zgeqrf_c(M,N,T) &
          & bind(c, name='PLASMA_Alloc_Workspace_zgeqrf')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_zgeqrf_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: T ! T is **, so pass by reference
          end function PLASMA_Alloc_Workspace_zgeqrf_c
      end interface

      interface
         function PLASMA_Alloc_Workspace_zgeqrf_Tile_c(M,N,descT) &
          & bind(c, name='PLASMA_Alloc_Workspace_zgeqrf_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_zgeqrf_Tile_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
          end function PLASMA_Alloc_Workspace_zgeqrf_Tile_c
      end interface

      interface
         function PLASMA_Alloc_Workspace_zgesv_incpiv_c(N,descL,IPIV) &
          & bind(c, name='PLASMA_Alloc_Workspace_zgesv_incpiv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_zgesv_incpiv_c
            integer(kind=c_int), value :: N
            type(c_ptr) :: descL ! descL is **, so pass by reference
            type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
          end function PLASMA_Alloc_Workspace_zgesv_incpiv_c
      end interface

      interface
         function PLASMA_Alloc_Workspace_zgesv_incpiv_Tile_c(N,descL,IPIV) &
          & bind(c, name='PLASMA_Alloc_Workspace_zgesv_incpiv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_zgesv_incpiv_Tile_c
            integer(kind=c_int), value :: N
            type(c_ptr) :: descL ! descL is **, so pass by reference
            type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
          end function PLASMA_Alloc_Workspace_zgesv_incpiv_Tile_c
      end interface

      interface
         function PLASMA_Alloc_Workspace_zgesvd_c(M,N,descT) &
          & bind(c, name='PLASMA_Alloc_Workspace_zgesvd')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_zgesvd_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
          end function PLASMA_Alloc_Workspace_zgesvd_c
      end interface

      interface
         function PLASMA_Alloc_Workspace_zgetrf_incpiv_c(M,N,descL,IPIV) &
          & bind(c, name='PLASMA_Alloc_Workspace_zgetrf_incpiv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_zgetrf_incpiv_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descL ! descL is **, so pass by reference
            type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
          end function PLASMA_Alloc_Workspace_zgetrf_incpiv_c
      end interface

      interface
         function PLASMA_Alloc_Workspace_zgetrf_incpiv_Tile_c(N,descL,IPIV) &
          & bind(c, name='PLASMA_Alloc_Workspace_zgetrf_incpiv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_zgetrf_incpiv_Tile_c
            integer(kind=c_int), value :: N
            type(c_ptr) :: descL ! descL is **, so pass by reference
            type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
          end function PLASMA_Alloc_Workspace_zgetrf_incpiv_Tile_c
      end interface

      interface
         function PLASMA_Alloc_Workspace_zgetri_Tile_Async_c(A,W) &
          & bind(c, name='PLASMA_Alloc_Workspace_zgetri_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_zgetri_Tile_Async_c
            type(c_ptr), value :: A
            type(c_ptr), value :: W
          end function PLASMA_Alloc_Workspace_zgetri_Tile_Async_c
      end interface

      interface
         function PLASMA_Alloc_Workspace_zheev_c(M,N,descT) &
          & bind(c, name='PLASMA_Alloc_Workspace_zheev')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_zheev_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
          end function PLASMA_Alloc_Workspace_zheev_c
      end interface

      interface
         function PLASMA_Alloc_Workspace_zheevd_c(M,N,descT) &
          & bind(c, name='PLASMA_Alloc_Workspace_zheevd')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_zheevd_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
          end function PLASMA_Alloc_Workspace_zheevd_c
      end interface

      interface
         function PLASMA_Alloc_Workspace_zhegv_c(M,N,descT) &
          & bind(c, name='PLASMA_Alloc_Workspace_zhegv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_zhegv_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
          end function PLASMA_Alloc_Workspace_zhegv_c
      end interface

      interface
         function PLASMA_Alloc_Workspace_zhegvd_c(M,N,descT) &
          & bind(c, name='PLASMA_Alloc_Workspace_zhegvd')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_zhegvd_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
          end function PLASMA_Alloc_Workspace_zhegvd_c
      end interface

      interface
         function PLASMA_Alloc_Workspace_zhetrd_c(M,N,descT) &
          & bind(c, name='PLASMA_Alloc_Workspace_zhetrd')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_zhetrd_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
          end function PLASMA_Alloc_Workspace_zhetrd_c
      end interface

  contains

       subroutine PLASMA_zgebrd(M,N,A,LDA,D,E,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         real(kind=c_double), intent(out), target :: D(*)
         real(kind=c_double), intent(out), target :: E(*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zgebrd_c(M,N,c_loc(A),LDA,c_loc(D),c_loc(E),T)
      end subroutine PLASMA_zgebrd

      subroutine PLASMA_zgelqf(M,N,A,LDA,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zgelqf_c(M,N,c_loc(A),LDA,T)
      end subroutine PLASMA_zgelqf

      subroutine PLASMA_zgelqs(M,N,NRHS,A,LDA,T,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zgelqs_c(M,N,NRHS,c_loc(A),LDA,T,c_loc(B),LDB)
      end subroutine PLASMA_zgelqs

      subroutine PLASMA_zgels(trans,M,N,NRHS,A,LDA,T,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         integer(kind=c_int), intent(in) :: trans
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zgels_c(trans,M,N,NRHS,c_loc(A),LDA,T,c_loc(B),LDB)
      end subroutine PLASMA_zgels

      subroutine PLASMA_zgemm(transA,transB,M,N,K,alpha,A,LDA,B,LDB,beta,C,LDC,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: K
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: LDC
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: transA
         integer(kind=c_int), intent(in) :: transB
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in) :: beta
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(in), target :: B(LDB,*)
         complex(kind=c_double_complex), intent(inout), target :: C(LDC,*)
         info = PLASMA_zgemm_c(transA,transB,M,N,K,alpha,c_loc(A),LDA,c_loc(B),LDB,beta,c_loc(C),LDC)
      end subroutine PLASMA_zgemm

      subroutine PLASMA_zgeqrf(M,N,A,LDA,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zgeqrf_c(M,N,c_loc(A),LDA,T)
      end subroutine PLASMA_zgeqrf

      subroutine PLASMA_zgeqrs(M,N,NRHS,A,LDA,T,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zgeqrs_c(M,N,NRHS,c_loc(A),LDA,T,c_loc(B),LDB)
      end subroutine PLASMA_zgeqrs

      subroutine PLASMA_zgesv(N,NRHS,A,LDA,IPIV,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         integer(kind=c_int), intent(out), target :: IPIV(*)
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
         info = PLASMA_zgesv_c(N,NRHS,c_loc(A),LDA,c_loc(IPIV),c_loc(B),LDB)
      end subroutine PLASMA_zgesv

      subroutine PLASMA_zgesv_incpiv(N,NRHS,A,LDA,L,IPIV,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
         type(c_ptr), value :: IPIV ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: L ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zgesv_incpiv_c(N,NRHS,c_loc(A),LDA,L,IPIV,c_loc(B),LDB)
      end subroutine PLASMA_zgesv_incpiv

      subroutine PLASMA_zgesvd(jobu,jobvt,M,N,A,LDA,S,U,LDU,VT,LDVT,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDU
         integer(kind=c_int), intent(in) :: LDVT
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: jobu
         integer(kind=c_int), intent(in) :: jobvt
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(out), target :: U(LDU,*)
         complex(kind=c_double_complex), intent(out), target :: VT(LDVT,*)
         real(kind=c_double), intent(out), target :: S(*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zgesvd_c(jobu,jobvt,M,N,c_loc(A),LDA,c_loc(S),c_loc(U),LDU,c_loc(VT),LDVT,T)
      end subroutine PLASMA_zgesvd

      subroutine PLASMA_zgetrf(M,N,A,LDA,IPIV,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(out), target :: IPIV(*)
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         info = PLASMA_zgetrf_c(M,N,c_loc(A),LDA,c_loc(IPIV))
      end subroutine PLASMA_zgetrf

      subroutine PLASMA_zgetrf_incpiv(M,N,A,LDA,L,IPIV,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         type(c_ptr), value :: IPIV ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: L    ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zgetrf_incpiv_c(M,N,c_loc(A),LDA,L,IPIV)
      end subroutine PLASMA_zgetrf_incpiv

      subroutine PLASMA_zgetrs(trans,N,NRHS,A,LDA,IPIV,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in), target :: IPIV(*)
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
         info = PLASMA_zgetrs_c(trans,N,NRHS,c_loc(A),LDA,c_loc(IPIV),c_loc(B),LDB)
      end subroutine PLASMA_zgetrs

      subroutine PLASMA_zgetrs_incpiv(trans,N,NRHS,A,LDA,L,IPIV,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         integer(kind=c_int), intent(in) :: trans
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
         type(c_ptr), value :: IPIV ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: L    ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zgetrs_incpiv_c(trans,N,NRHS,c_loc(A),LDA,L,IPIV,c_loc(B),LDB)
      end subroutine PLASMA_zgetrs_incpiv

#if defined(PRECISION_z) || defined(PRECISION_c)
      subroutine PLASMA_zhemm(side,uplo,M,N,alpha,A,LDA,B,LDB,beta,C,LDC,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: LDC
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in) :: beta
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(in), target :: B(LDB,*)
         complex(kind=c_double_complex), intent(inout), target :: C(LDC,*)
         info = PLASMA_zhemm_c(side,uplo,M,N,alpha,c_loc(A),LDA,c_loc(B),LDB,beta,c_loc(C),LDC)
      end subroutine PLASMA_zhemm

      subroutine PLASMA_zherk(uplo,trans,N,K,alpha,A,LDA,beta,C,LDC,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: K
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDC
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(inout), target :: C(LDC,*)
         real(kind=c_double), intent(in) :: alpha
         real(kind=c_double), intent(in) :: beta
         info = PLASMA_zherk_c(uplo,trans,N,K,alpha,c_loc(A),LDA,beta,c_loc(C),LDC)
      end subroutine PLASMA_zherk

      subroutine PLASMA_zher2k(uplo,trans,N,K,alpha,A,LDA,B,LDB,beta,C,LDC,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: K
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: LDC
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(in), target :: B(LDB,*)
         complex(kind=c_double_complex), intent(inout), target :: C(LDC,*)
         real(kind=c_double), intent(in) :: beta
         info = PLASMA_zher2k_c(uplo,trans,N,K,alpha,c_loc(A),LDA,c_loc(B),LDB,beta,c_loc(C),LDC)
      end subroutine PLASMA_zher2k
#endif

      subroutine PLASMA_zheev(jobz,uplo,N,A,LDA,W,T,Q,LDQ,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDQ
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         real(kind=c_double), intent(out), target :: W(*)
         complex(kind=c_double_complex), intent(out), target :: Q(LDQ,*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zheev_c(jobz,uplo,N,c_loc(A),LDA,c_loc(W),T,c_loc(Q),LDQ)
      end subroutine PLASMA_zheev

      subroutine PLASMA_zheevd(jobz,uplo,N,A,LDA,W,T,Q,LDQ,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDQ
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         real(kind=c_double), intent(out), target :: W(*)
         complex(kind=c_double_complex), intent(out), target :: Q(LDQ,*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zheevd_c(jobz,uplo,N,c_loc(A),LDA,c_loc(W),T,c_loc(Q),LDQ)
      end subroutine PLASMA_zheevd

      subroutine PLASMA_zhegv(itype,jobz,uplo,N,A,LDA,B,LDB,W,T,Q,LDQ,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: LDQ
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: itype
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
         real(kind=c_double), intent(out), target :: W(*)
         complex(kind=c_double_complex), intent(out), target :: Q(LDQ,*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zhegv_c(itype,jobz,uplo,N,c_loc(A),LDA,c_loc(B),LDB,c_loc(W),T,c_loc(Q),LDQ)
      end subroutine PLASMA_zhegv

      subroutine PLASMA_zhegvd(itype,jobz,uplo,N,A,LDA,B,LDB,W,T,Q,LDQ,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: LDQ
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: itype
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
         real(kind=c_double), intent(out), target :: W(*)
         complex(kind=c_double_complex), intent(out), target :: Q(LDQ,*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zhegvd_c(itype,jobz,uplo,N,c_loc(A),LDA,c_loc(B),LDB,c_loc(W),T,c_loc(Q),LDQ)
      end subroutine PLASMA_zhegvd

      subroutine PLASMA_zhegst(itype,uplo,N,A,LDA,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: itype
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in), target :: B(LDB,*)
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         info = PLASMA_zhegst_c(itype,uplo,N,c_loc(A),LDA,c_loc(B),LDB)
      end subroutine PLASMA_zhegst

      subroutine PLASMA_zhetrd(jobz,uplo,N,A,LDA,D,E,descT,Q,LDQ,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDQ
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         real(kind=c_double), intent(out), target :: D(*)
         real(kind=c_double), intent(out), target :: E(*)
         type(c_ptr), value :: descT ! Arg managed by PLASMA: opaque to Fortran
         complex(kind=c_double_complex), intent(inout), target :: Q(LDQ,*)
         info = PLASMA_zhetrd_c(jobz,uplo,N,c_loc(A),LDA,c_loc(D),c_loc(E),descT,c_loc(Q),LDQ)
      end subroutine PLASMA_zhetrd

      function PLASMA_zlange(norm,M,N,A,LDA)
         use iso_c_binding
         implicit none
         real(kind=c_double) :: PLASMA_zlange
         integer(kind=c_int), intent(in) :: norm
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         PLASMA_zlange = PLASMA_zlange_c(norm,M,N,c_loc(A),LDA)
      end function PLASMA_zlange

#if defined(PRECISION_z) || defined(PRECISION_c)
      function PLASMA_zlanhe(norm,uplo,N,A,LDA)
         use iso_c_binding
         implicit none
         real(kind=c_double) :: PLASMA_zlanhe
         integer(kind=c_int), intent(in) :: norm
         integer(kind=c_int), intent(in) :: uplo
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: N
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         PLASMA_zlanhe = PLASMA_zlanhe_c(norm,uplo,N,c_loc(A),LDA)
      end function PLASMA_zlanhe
#endif

      function PLASMA_zlansy(norm,uplo,N,A,LDA)
         use iso_c_binding
         implicit none
         real(kind=c_double) :: PLASMA_zlansy
         integer(kind=c_int), intent(in) :: norm
         integer(kind=c_int), intent(in) :: uplo
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: N
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         PLASMA_zlansy = PLASMA_zlansy_c(norm,uplo,N,c_loc(A),LDA)
      end function PLASMA_zlansy

      subroutine PLASMA_zlaswp(N,A,LDA,K1,K2,IPIV,INCX,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: INCX
         integer(kind=c_int), intent(in) :: K1
         integer(kind=c_int), intent(in) :: K2
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in), target :: IPIV(*)
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         info = PLASMA_zlaswp_c(N,c_loc(A),LDA,K1,K2,c_loc(IPIV),INCX)
      end subroutine PLASMA_zlaswp

      subroutine PLASMA_zlauum(uplo,N,A,LDA,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         info = PLASMA_zlauum_c(uplo,N,c_loc(A),LDA)
      end subroutine PLASMA_zlauum

#if defined(PRECISION_z) || defined(PRECISION_c)
      subroutine PLASMA_zplghe(bump,N,A,LDA,seed,info) 
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         real(kind=c_double), intent(in) :: bump
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_long_long), intent(in) :: seed
         complex(kind=c_double_complex), intent(out), target :: A(LDA,*)
         info = PLASMA_zplghe_c(bump,N,c_loc(A),LDA,seed)
       end subroutine PLASMA_zplghe
#endif

      subroutine PLASMA_zplgsy(bump,N,A,LDA,seed,info) 
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         complex(kind=c_double_complex), intent(in) :: bump
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_long_long), intent(in) :: seed
         complex(kind=c_double_complex), intent(out), target :: A(LDA,*)
         info = PLASMA_zplgsy_c(bump,N,c_loc(A),LDA,seed)
      end subroutine PLASMA_zplgsy
       
      subroutine PLASMA_zplrnt(M,N,A,LDA,seed,info) 
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_long_long), intent(in) :: seed
         complex(kind=c_double_complex), intent(out), target :: A(LDA,*)
         info = PLASMA_zplrnt_c(M,N,c_loc(A),LDA,seed)
      end subroutine PLASMA_zplrnt
       
      subroutine PLASMA_zposv(uplo,N,NRHS,A,LDA,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
         info = PLASMA_zposv_c(uplo,N,NRHS,c_loc(A),LDA,c_loc(B),LDB)
      end subroutine PLASMA_zposv

      subroutine PLASMA_zpotrf(uplo,N,A,LDA,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         info = PLASMA_zpotrf_c(uplo,N,c_loc(A),LDA)
      end subroutine PLASMA_zpotrf

      subroutine PLASMA_zpotri(uplo,N,A,LDA,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         info = PLASMA_zpotri_c(uplo,N,c_loc(A),LDA)
      end subroutine PLASMA_zpotri

      subroutine PLASMA_zpotrs(uplo,N,NRHS,A,LDA,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
         info = PLASMA_zpotrs_c(uplo,N,NRHS,c_loc(A),LDA,c_loc(B),LDB)
      end subroutine PLASMA_zpotrs

      subroutine PLASMA_zsymm(side,uplo,M,N,alpha,A,LDA,B,LDB,beta,C,LDC,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: LDC
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in) :: beta
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(in), target :: B(LDB,*)
         complex(kind=c_double_complex), intent(inout), target :: C(LDC,*)
         info = PLASMA_zsymm_c(side,uplo,M,N,alpha,c_loc(A),LDA,c_loc(B),LDB,beta,c_loc(C),LDC)
      end subroutine PLASMA_zsymm

      subroutine PLASMA_zsyrk(uplo,trans,N,K,alpha,A,LDA,beta,C,LDC,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: K
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDC
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in) :: beta
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(inout), target :: C(LDC,*)
         info = PLASMA_zsyrk_c(uplo,trans,N,K,alpha,c_loc(A),LDA,beta,c_loc(C),LDC)
      end subroutine PLASMA_zsyrk

      subroutine PLASMA_zsyr2k(uplo,trans,N,K,alpha,A,LDA,B,LDB,beta,C,LDC,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: K
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: LDC
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in) :: beta
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(in), target :: B(LDB,*)
         complex(kind=c_double_complex), intent(inout), target :: C(LDC,*)
         info = PLASMA_zsyr2k_c(uplo,trans,N,K,alpha,c_loc(A),LDA,c_loc(B),LDB,beta,c_loc(C),LDC)
      end subroutine PLASMA_zsyr2k

      subroutine PLASMA_ztrmm(side,uplo,transA,diag,N,NRHS,alpha,A,LDA,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         integer(kind=c_int), intent(in) :: diag
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: transA
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
         info = PLASMA_ztrmm_c(side,uplo,transA,diag,N,NRHS,alpha,c_loc(A),LDA,c_loc(B),LDB)
      end subroutine PLASMA_ztrmm

      subroutine PLASMA_ztrsm(side,uplo,transA,diag,N,NRHS,alpha,A,LDA,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         integer(kind=c_int), intent(in) :: diag
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: transA
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
         info = PLASMA_ztrsm_c(side,uplo,transA,diag,N,NRHS,alpha,c_loc(A),LDA,c_loc(B),LDB)
      end subroutine PLASMA_ztrsm

      subroutine PLASMA_ztrsmpl(N,NRHS,A,LDA,L,IPIV,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         type(c_ptr), value :: L    ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: IPIV ! Arg managed by PLASMA: opaque to Fortran
         complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
         info = PLASMA_ztrsmpl_c(N,NRHS,c_loc(A),LDA,L,IPIV,c_loc(B),LDB)
      end subroutine PLASMA_ztrsmpl

      subroutine PLASMA_ztrtri(uplo,diag,N,A,LDA,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: diag
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         info = PLASMA_ztrtri_c(uplo,diag,N,c_loc(A),LDA)
      end subroutine PLASMA_ztrtri

      subroutine PLASMA_zunglq(M,N,K,A,LDA,T,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         integer(kind=c_int), intent(in) :: K
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(out), target :: B(LDB,*)
         info = PLASMA_zunglq_c(M,N,K,c_loc(A),LDA,T,c_loc(B),LDB)
      end subroutine PLASMA_zunglq

      subroutine PLASMA_zungqr(M,N,K,A,LDA,T,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         integer(kind=c_int), intent(in) :: K
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(out), target :: B(LDB,*)
         info = PLASMA_zungqr_c(M,N,K,c_loc(A),LDA,T,c_loc(B),LDB)
      end subroutine PLASMA_zungqr

      subroutine PLASMA_zunmlq(side,trans,M,N,K,A,LDA,T,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         integer(kind=c_int), intent(in) :: K
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: trans
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
         info = PLASMA_zunmlq_c(side,trans,M,N,K,c_loc(A),LDA,T,c_loc(B),LDB)
      end subroutine PLASMA_zunmlq

      subroutine PLASMA_zunmqr(side,trans,M,N,K,A,LDA,T,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         integer(kind=c_int), intent(in) :: K
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: trans
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
         info = PLASMA_zunmqr_c(side,trans,M,N,K,c_loc(A),LDA,T,c_loc(B),LDB)
      end subroutine PLASMA_zunmqr

      subroutine PLASMA_zgecfi(m,n,A,fin,imb,inb,fout,omb,onb,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         complex(kind=c_double_complex), intent(inout), target :: A(*)
         integer(kind=c_int), intent(in) :: fin
         integer(kind=c_int), intent(in) :: fout
         integer(kind=c_int), intent(in) :: imb
         integer(kind=c_int), intent(in) :: inb
         integer(kind=c_int), intent(in) :: omb
         integer(kind=c_int), intent(in) :: onb
         integer(kind=c_int), intent(in) :: m
         integer(kind=c_int), intent(in) :: n
         info = PLASMA_zgecfi_c(m,n,c_loc(A),fin,imb,inb,fout,omb,onb)
      end subroutine PLASMA_zgecfi

      subroutine PLASMA_zgetmi(m,n,A,fin,mb,nb,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         complex(kind=c_double_complex), intent(inout), target :: A(*)
         integer(kind=c_int), intent(in) :: fin
         integer(kind=c_int), intent(in) :: mb
         integer(kind=c_int), intent(in) :: nb
         integer(kind=c_int), intent(in) :: m
         integer(kind=c_int), intent(in) :: n
         info = PLASMA_zgetmi_c(m,n,c_loc(A),fin,mb,nb)
      end subroutine PLASMA_zgetmi

      subroutine PLASMA_zgetri(N,A,LDA,IPIV,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in), target :: IPIV(*)
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         info = PLASMA_zgetri_c(N,c_loc(A),LDA,c_loc(IPIV))
      end subroutine PLASMA_zgetri

      subroutine PLASMA_zlacpy(uplo,M,N,A,LDA,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(out), target :: B(LDB,*)
         info = PLASMA_zlacpy_c(uplo,M,N,c_loc(A),LDA,c_loc(B),LDB)
      end subroutine PLASMA_zlacpy

      subroutine PLASMA_zlaset(uplo,M,N,alpha,beta,A,LDA,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in) :: beta
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         info = PLASMA_zlaset_c(uplo,M,N,alpha,beta,c_loc(A),LDA)
      end subroutine PLASMA_zlaset

      subroutine PLASMA_zlaswpc(N,A,LDA,K1,K2,IPIV,INCX,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in), target :: IPIV(*)
         integer(kind=c_int), intent(in) :: INCX
         integer(kind=c_int), intent(in) :: K1
         integer(kind=c_int), intent(in) :: K2
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: N
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         info = PLASMA_zlaswpc_c(N,c_loc(A),LDA,K1,K2,c_loc(IPIV),INCX)
      end subroutine PLASMA_zlaswpc

      subroutine PLASMA_ztrsmrv(side,uplo,transA,diag,N,NRHS,alpha,A,LDA,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: diag
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: transA
         integer(kind=c_int), intent(in) :: uplo
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
         info = PLASMA_ztrsmrv_c(side,uplo,transA,diag,N,NRHS,alpha,c_loc(A),LDA,c_loc(B),LDB)
      end subroutine PLASMA_ztrsmrv

      subroutine PLASMA_zgebrd_Tile(A,D,E,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         real(kind=c_double), intent(out), target :: D(*)
         real(kind=c_double), intent(out), target :: E(*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zgebrd_Tile_c(A,c_loc(D),c_loc(E),T)
      end subroutine PLASMA_zgebrd_Tile

      subroutine PLASMA_zgelqf_Tile(A,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zgelqf_Tile_c(A,T)
      end subroutine PLASMA_zgelqf_Tile

      subroutine PLASMA_zgelqs_Tile(A,T,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zgelqs_Tile_c(A,T,B)
      end subroutine PLASMA_zgelqs_Tile

      subroutine PLASMA_zgels_Tile(trans,A,T,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: trans
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zgels_Tile_c(trans,A,T,B)
      end subroutine PLASMA_zgels_Tile

      subroutine PLASMA_zgemm_Tile(transA,transB,alpha,A,B,beta,C,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: transA
         integer(kind=c_int), intent(in) :: transB
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zgemm_Tile_c(transA,transB,alpha,A,B,beta,C)
      end subroutine PLASMA_zgemm_Tile

      subroutine PLASMA_zgeqrf_Tile(A,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zgeqrf_Tile_c(A,T)
      end subroutine PLASMA_zgeqrf_Tile

      subroutine PLASMA_zgeqrs_Tile(A,T,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zgeqrs_Tile_c(A,T,B)
      end subroutine PLASMA_zgeqrs_Tile

      subroutine PLASMA_zgesv_Tile(A,IPIV,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(out), target :: IPIV(*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zgesv_Tile_c(A,c_loc(IPIV),B)
      end subroutine PLASMA_zgesv_Tile

      subroutine PLASMA_zgesv_incpiv_Tile(A,L,IPIV,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: IPIV ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: L ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zgesv_incpiv_Tile_c(A,L,IPIV,B)
      end subroutine PLASMA_zgesv_incpiv_Tile

      subroutine PLASMA_zgesvd_Tile(jobu,jobvt,A,S,U,VT,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: jobu
         integer(kind=c_int), intent(in) :: jobvt
         real(kind=c_double), intent(out), target :: S(*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: U ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: VT ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zgesvd_Tile_c(jobu,jobvt,A,c_loc(S),U,VT,T)
      end subroutine PLASMA_zgesvd_Tile

      subroutine PLASMA_zgetrf_Tile(A,IPIV,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(out), target :: IPIV(*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zgetrf_Tile_c(A,c_loc(IPIV))
      end subroutine PLASMA_zgetrf_Tile

      subroutine PLASMA_zgetrf_incpiv_Tile(A,L,IPIV,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: IPIV ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: L ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zgetrf_incpiv_Tile_c(A,L,IPIV)
      end subroutine PLASMA_zgetrf_incpiv_Tile

      subroutine PLASMA_zgetrs_Tile(trans,A,IPIV,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: trans
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         integer(kind=c_int), intent(in), target :: IPIV(*)
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zgetrs_Tile_c(trans,A,c_loc(IPIV),B)
      end subroutine PLASMA_zgetrs_Tile

      subroutine PLASMA_zgetrs_incpiv_Tile(A,L,IPIV,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: L ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: IPIV ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zgetrs_incpiv_Tile_c(A,L,IPIV,B)
      end subroutine PLASMA_zgetrs_incpiv_Tile

#if defined(PRECISION_z) || defined(PRECISION_c)
      subroutine PLASMA_zhemm_Tile(side,uplo,alpha,A,B,beta,C,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zhemm_Tile_c(side,uplo,alpha,A,B,beta,C)
      end subroutine PLASMA_zhemm_Tile

      subroutine PLASMA_zherk_Tile(uplo,trans,alpha,A,beta,C,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_double), intent(in) :: alpha
         real(kind=c_double), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zherk_Tile_c(uplo,trans,alpha,A,beta,C)
      end subroutine PLASMA_zherk_Tile

      subroutine PLASMA_zher2k_Tile(uplo,trans,alpha,A,B,beta,C,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         real(kind=c_double), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zher2k_Tile_c(uplo,trans,alpha,A,B,beta,C)
      end subroutine PLASMA_zher2k_Tile
#endif

      subroutine PLASMA_zheev_Tile(jobz,uplo,A,W,T,Q,LDQ,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         integer(kind=c_int), intent(in) :: LDQ
         real(kind=c_double), intent(out), target :: W(*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         complex(kind=c_double_complex), intent(out), target :: Q(LDQ,*)
         info = PLASMA_zheev_Tile_c(jobz,uplo,A,c_loc(W),T,c_loc(Q),LDQ)
      end subroutine PLASMA_zheev_Tile

      subroutine PLASMA_zheevd_Tile(jobz,uplo,A,W,T,Q,LDQ,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         integer(kind=c_int), intent(in) :: LDQ
         real(kind=c_double), intent(out), target :: W(*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         complex(kind=c_double_complex), intent(out), target :: Q(LDQ,*)
         info = PLASMA_zheevd_Tile_c(jobz,uplo,A,c_loc(W),T,c_loc(Q),LDQ)
      end subroutine PLASMA_zheevd_Tile

      subroutine PLASMA_zhegv_Tile(itype,jobz,uplo,A,B,W,T,Q,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: itype
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_double), intent(out), target :: W(*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: Q ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zhegv_Tile_c(itype,jobz,uplo,A,B,c_loc(W),T,Q)
      end subroutine PLASMA_zhegv_Tile

      subroutine PLASMA_zhegvd_Tile(itype,jobz,uplo,A,B,W,T,Q,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: itype
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_double), intent(out), target :: W(*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: Q ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zhegvd_Tile_c(itype,jobz,uplo,A,B,c_loc(W),T,Q)
      end subroutine PLASMA_zhegvd_Tile

      subroutine PLASMA_zhegst_Tile(itype,uplo,A,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: itype
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zhegst_Tile_c(itype,uplo,A,B)
      end subroutine PLASMA_zhegst_Tile

      subroutine PLASMA_zhetrd_Tile(jobz,uplo,A,D,E,T,Q,LDQ,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         integer(kind=c_int), intent(in) :: LDQ
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         real(kind=c_double), intent(out), target :: D(*)
         real(kind=c_double), intent(out), target :: E(*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         complex(kind=c_double_complex), intent(out), target :: Q(LDQ,*)
         info = PLASMA_zhetrd_Tile_c(jobz,uplo,A,c_loc(D),c_loc(E),T,c_loc(Q),LDQ)
      end subroutine PLASMA_zhetrd_Tile

      function PLASMA_zlange_Tile(norm,A)
         use iso_c_binding
         implicit none
         real(kind=c_double) :: PLASMA_zlange_Tile
         integer(kind=c_int), intent(in) :: norm
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         PLASMA_zlange_Tile = PLASMA_zlange_Tile_c(norm,A)
       end function PLASMA_zlange_Tile

#if defined(PRECISION_z) || defined(PRECISION_c)
      function PLASMA_zlanhe_Tile(norm,uplo,A)
         use iso_c_binding
         implicit none
         real(kind=c_double) :: PLASMA_zlanhe_Tile
         integer(kind=c_int), intent(in) :: norm
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         PLASMA_zlanhe_Tile = PLASMA_zlanhe_Tile_c(norm,uplo,A)
      end function PLASMA_zlanhe_Tile
#endif

      function PLASMA_zlansy_Tile(norm,uplo,A)
         use iso_c_binding
         implicit none
         real(kind=c_double) :: PLASMA_zlansy_Tile
         integer(kind=c_int), intent(in) :: norm
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         PLASMA_zlansy_Tile = PLASMA_zlansy_Tile_c(norm,uplo,A)
      end function PLASMA_zlansy_Tile

      subroutine PLASMA_zlaswp_Tile(A,K1,K2,IPIV,INCX,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: INCX
         integer(kind=c_int), intent(in) :: K1
         integer(kind=c_int), intent(in) :: K2
         integer(kind=c_int), intent(in), target :: IPIV(*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zlaswp_Tile_c(A,K1,K2,c_loc(IPIV),INCX)
      end subroutine PLASMA_zlaswp_Tile

      subroutine PLASMA_zlauum_Tile(uplo,A,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zlauum_Tile_c(uplo,A)
      end subroutine PLASMA_zlauum_Tile

#if defined(PRECISION_z) || defined(PRECISION_c)
      subroutine PLASMA_zplghe_Tile(bump,A,seed,info) 
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         real(kind=c_double), intent(in) :: bump
         integer(kind=c_long_long), intent(in) :: seed
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zplghe_Tile_c(bump,A,seed)
      end subroutine PLASMA_zplghe_Tile
#endif

      subroutine PLASMA_zplgsy_Tile(bump,A,seed,info) 
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         complex(kind=c_double_complex), intent(in) :: bump
         integer(kind=c_long_long), intent(in) :: seed
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zplgsy_Tile_c(bump,A,seed)
      end subroutine PLASMA_zplgsy_Tile
       
      subroutine PLASMA_zplrnt_Tile(A,seed,info) 
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_long_long), intent(in) :: seed
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zplrnt_Tile_c(A,seed)
      end subroutine PLASMA_zplrnt_Tile
       
      subroutine PLASMA_zposv_Tile(uplo,A,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zposv_Tile_c(uplo,A,B)
      end subroutine PLASMA_zposv_Tile

      subroutine PLASMA_zpotrf_Tile(uplo,A,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zpotrf_Tile_c(uplo,A)
      end subroutine PLASMA_zpotrf_Tile

      subroutine PLASMA_zpotri_Tile(uplo,A,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zpotri_Tile_c(uplo,A)
      end subroutine PLASMA_zpotri_Tile

      subroutine PLASMA_zpotrs_Tile(uplo,A,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zpotrs_Tile_c(uplo,A,B)
      end subroutine PLASMA_zpotrs_Tile

      subroutine PLASMA_zsymm_Tile(side,uplo,alpha,A,B,beta,C,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zsymm_Tile_c(side,uplo,alpha,A,B,beta,C)
      end subroutine PLASMA_zsymm_Tile

      subroutine PLASMA_zsyrk_Tile(uplo,trans,alpha,A,beta,C,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zsyrk_Tile_c(uplo,trans,alpha,A,beta,C)
      end subroutine PLASMA_zsyrk_Tile

      subroutine PLASMA_zsyr2k_Tile(uplo,trans,alpha,A,B,beta,C,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zsyr2k_Tile_c(uplo,trans,alpha,A,B,beta,C)
      end subroutine PLASMA_zsyr2k_Tile

      subroutine PLASMA_ztrmm_Tile(side,uplo,transA,diag,alpha,A,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: diag
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: transA
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_ztrmm_Tile_c(side,uplo,transA,diag,alpha,A,B)
      end subroutine PLASMA_ztrmm_Tile

      subroutine PLASMA_ztrsm_Tile(side,uplo,transA,diag,alpha,A,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: diag
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: transA
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_ztrsm_Tile_c(side,uplo,transA,diag,alpha,A,B)
      end subroutine PLASMA_ztrsm_Tile

      subroutine PLASMA_ztrsmpl_Tile(A,L,IPIV,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: L ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: IPIV ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_ztrsmpl_Tile_c(A,L,IPIV,B)
      end subroutine PLASMA_ztrsmpl_Tile

      subroutine PLASMA_ztrtri_Tile(uplo,diag,A,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: diag
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_ztrtri_Tile_c(uplo,diag,A)
      end subroutine PLASMA_ztrtri_Tile

      subroutine PLASMA_zunglq_Tile(A,T,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zunglq_Tile_c(A,T,B)
      end subroutine PLASMA_zunglq_Tile

      subroutine PLASMA_zungqr_Tile(A,T,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zungqr_Tile_c(A,T,B)
      end subroutine PLASMA_zungqr_Tile

      subroutine PLASMA_zunmlq_Tile(side,trans,A,T,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: trans
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zunmlq_Tile_c(side,trans,A,T,B)
      end subroutine PLASMA_zunmlq_Tile

      subroutine PLASMA_zunmqr_Tile(side,trans,A,T,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: trans
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zunmqr_Tile_c(side,trans,A,T,B)
      end subroutine PLASMA_zunmqr_Tile

      subroutine PLASMA_zgetri_Tile(A,IPIV,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in), target :: IPIV(*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zgetri_Tile_c(A,c_loc(IPIV))
      end subroutine PLASMA_zgetri_Tile

      subroutine PLASMA_zlacpy_Tile(uplo,A,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zlacpy_Tile_c(uplo,A,B)
      end subroutine PLASMA_zlacpy_Tile

      subroutine PLASMA_zlaset_Tile(uplo,alpha,beta,A,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zlaset_Tile_c(uplo,alpha,beta,A)
      end subroutine PLASMA_zlaset_Tile

      subroutine PLASMA_zlaswpc_Tile(A,K1,K2,IPIV,INCX,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in), target :: IPIV(*)
         integer(kind=c_int), intent(in) :: INCX
         integer(kind=c_int), intent(in) :: K1
         integer(kind=c_int), intent(in) :: K2
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zlaswpc_Tile_c(A,K1,K2,c_loc(IPIV),INCX)
      end subroutine PLASMA_zlaswpc_Tile

      subroutine PLASMA_ztrsmrv_Tile(side,uplo,transA,diag,alpha,A,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: diag
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: transA
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_ztrsmrv_Tile_c(side,uplo,transA,diag,alpha,A,B)
      end subroutine PLASMA_ztrsmrv_Tile

      subroutine PLASMA_zgetri_Tile_Async(A,IPIV,W,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in), target :: IPIV(*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: W ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zgetri_Tile_Async_c(A,c_loc(IPIV),W,sequence,request)
      end subroutine PLASMA_zgetri_Tile_Async

      subroutine PLASMA_zlange_Tile_Async(norm,A,value,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int) :: info
         real(kind=c_double), intent(out), target :: value
         integer(kind=c_int), intent(in) :: norm
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zlange_Tile_Async_c(norm,A,c_loc(value),sequence,request)
      end subroutine PLASMA_zlange_Tile_Async

      subroutine PLASMA_zlansy_Tile_Async(norm,uplo,A,value,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int) :: info
         real(kind=c_double), intent(out), target :: value
         integer(kind=c_int), intent(in) :: norm
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zlansy_Tile_Async_c(norm,uplo,A,c_loc(value),sequence,request)
      end subroutine PLASMA_zlansy_Tile_Async

      subroutine PLASMA_zgebrd_Tile_Async(A,D,E,T,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         real(kind=c_double), intent(out), target :: D(*)
         real(kind=c_double), intent(out), target :: E(*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zgebrd_Tile_Async_c(A,c_loc(D),c_loc(E),T,sequence,request)
      end subroutine PLASMA_zgebrd_Tile_Async

      subroutine PLASMA_zgecfi_Async(m,n,A,fin,imb,inb,fout,omb,onb,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         complex(kind=c_double_complex), intent(inout), target :: A(*)
         integer(kind=c_int), intent(in) :: fin
         integer(kind=c_int), intent(in) :: fout
         integer(kind=c_int), intent(in) :: imb
         integer(kind=c_int), intent(in) :: inb
         integer(kind=c_int), intent(in) :: omb
         integer(kind=c_int), intent(in) :: onb
         integer(kind=c_int), intent(in) :: m
         integer(kind=c_int), intent(in) :: n
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zgecfi_Async_c(m,n,c_loc(A),fin,imb,inb,fout,omb,onb,sequence,request)
      end subroutine PLASMA_zgecfi_Async

      subroutine PLASMA_zgelqf_Tile_Async(A,T,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zgelqf_Tile_Async_c(A,T,sequence,request)
      end subroutine PLASMA_zgelqf_Tile_Async

      subroutine PLASMA_zgelqs_Tile_Async(A,T,B,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zgelqs_Tile_Async_c(A,T,B,sequence,request)
      end subroutine PLASMA_zgelqs_Tile_Async

      subroutine PLASMA_zgels_Tile_Async(trans,A,T,B,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: trans
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zgels_Tile_Async_c(trans,A,T,B,sequence,request)
      end subroutine PLASMA_zgels_Tile_Async

      subroutine PLASMA_zgemm_Tile_Async(transA,transB,alpha,A,B,beta,C,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: transA
         integer(kind=c_int), intent(in) :: transB
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zgemm_Tile_Async_c(transA,transB,alpha,A,B,beta,C,sequence,request)
      end subroutine PLASMA_zgemm_Tile_Async

      subroutine PLASMA_zgeqrf_Tile_Async(A,T,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zgeqrf_Tile_Async_c(A,T,sequence,request)
      end subroutine PLASMA_zgeqrf_Tile_Async

      subroutine PLASMA_zgeqrs_Tile_Async(A,T,B,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zgeqrs_Tile_Async_c(A,T,B,sequence,request)
      end subroutine PLASMA_zgeqrs_Tile_Async

      subroutine PLASMA_zgesv_Tile_Async(A,IPIV,B,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(out), target :: IPIV(*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zgesv_Tile_Async_c(A,c_loc(IPIV),B,sequence,request)
      end subroutine PLASMA_zgesv_Tile_Async

      subroutine PLASMA_zgesv_incpiv_Tile_Async(A,L,IPIV,B,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: IPIV ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: L ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zgesv_incpiv_Tile_Async_c(A,L,IPIV,B,sequence,request)
      end subroutine PLASMA_zgesv_incpiv_Tile_Async

      subroutine PLASMA_zgesvd_Tile_Async(jobu,jobvt,A,S,U,VT,T,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: jobu
         integer(kind=c_int), intent(in) :: jobvt
         real(kind=c_double), intent(out), target :: S(*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: U ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: VT ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zgesvd_Tile_Async_c(jobu,jobvt,A,c_loc(S),U,VT,T,sequence,request)
      end subroutine PLASMA_zgesvd_Tile_Async

      subroutine PLASMA_zgetmi_Async(m,n,A,fin,mb,nb,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         complex(kind=c_double_complex), intent(inout), target :: A(*)
         integer(kind=c_int), intent(in) :: fin
         integer(kind=c_int), intent(in) :: mb
         integer(kind=c_int), intent(in) :: nb
         integer(kind=c_int), intent(in) :: m
         integer(kind=c_int), intent(in) :: n
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zgetmi_Async_c(m,n,c_loc(A),fin,mb,nb,sequence,request)
      end subroutine PLASMA_zgetmi_Async

      subroutine PLASMA_zgetrf_Tile_Async(A,IPIV,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(out), target :: IPIV(*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zgetrf_Tile_Async_c(A,c_loc(IPIV),sequence,request)
      end subroutine PLASMA_zgetrf_Tile_Async

      subroutine PLASMA_zgetrf_incpiv_Tile_Async(A,L,IPIV,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: IPIV ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: L ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zgetrf_incpiv_Tile_Async_c(A,L,IPIV,sequence,request)
      end subroutine PLASMA_zgetrf_incpiv_Tile_Async

      subroutine PLASMA_zgetrs_Tile_Async(trans,A,IPIV,B,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: trans
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         integer(kind=c_int), intent(in), target :: IPIV(*)
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zgetrs_Tile_Async_c(trans,A,c_loc(IPIV),B,sequence,request)
      end subroutine PLASMA_zgetrs_Tile_Async

      subroutine PLASMA_zgetrs_incpiv_Tile_Async(A,L,IPIV,B,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: L ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: IPIV ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zgetrs_incpiv_Tile_Async_c(A,L,IPIV,B,sequence,request)
      end subroutine PLASMA_zgetrs_incpiv_Tile_Async

#if defined(PRECISION_z) || defined(PRECISION_c)
      subroutine PLASMA_zhemm_Tile_Async(side,uplo,alpha,A,B,beta,C,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zhemm_Tile_Async_c(side,uplo,alpha,A,B,beta,C,sequence,request)
      end subroutine PLASMA_zhemm_Tile_Async

      subroutine PLASMA_zherk_Tile_Async(uplo,trans,alpha,A,beta,C,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_double), intent(in) :: alpha
         real(kind=c_double), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zherk_Tile_Async_c(uplo,trans,alpha,A,beta,C,sequence,request)
      end subroutine PLASMA_zherk_Tile_Async

      subroutine PLASMA_zher2k_Tile_Async(uplo,trans,alpha,A,B,beta,C,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         real(kind=c_double), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zher2k_Tile_Async_c(uplo,trans,alpha,A,B,beta,C,sequence,request)
      end subroutine PLASMA_zher2k_Tile_Async
#endif

      subroutine PLASMA_zheev_Tile_Async(jobz,uplo,A,W,T,Q,LDQ,sequence,request,info) 
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         integer(kind=c_int), intent(in) :: LDQ
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         real(kind=c_double), intent(out), target :: W(*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         complex(kind=c_double_complex), intent(out), target :: Q(LDQ, *)
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zheev_Tile_Async_c(jobz,uplo,A,c_loc(W),T,c_loc(Q),LDQ,sequence,request)
      end subroutine PLASMA_zheev_Tile_Async

      subroutine PLASMA_zheevd_Tile_Async(jobz,uplo,A,W,T,Q,LDQ,sequence,request,info) 
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         integer(kind=c_int), intent(in) :: LDQ
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         real(kind=c_double), intent(out), target :: W(*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         complex(kind=c_double_complex), intent(out), target :: Q(LDQ, *)
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zheevd_Tile_Async_c(jobz,uplo,A,c_loc(W),T,c_loc(Q),LDQ,sequence,request)
      end subroutine PLASMA_zheevd_Tile_Async

      subroutine PLASMA_zhegv_Tile_Async(itype,jobz,uplo,A,B,W,T,Q,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: itype
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_double), intent(out), target :: W(*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: Q ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zhegv_Tile_Async_c(itype,jobz,uplo,A,B,c_loc(W),T,Q,sequence,request)
      end subroutine PLASMA_zhegv_Tile_Async

      subroutine PLASMA_zhegvd_Tile_Async(itype,jobz,uplo,A,B,W,T,Q,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: itype
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_double), intent(out), target :: W(*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: Q ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zhegvd_Tile_Async_c(itype,jobz,uplo,A,B,c_loc(W),T,Q,sequence,request)
      end subroutine PLASMA_zhegvd_Tile_Async

      subroutine PLASMA_zhegst_Tile_Async(itype,uplo,A,B,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: itype
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zhegst_Tile_Async_c(itype,uplo,A,B,sequence,request)
      end subroutine PLASMA_zhegst_Tile_Async

      subroutine PLASMA_zhetrd_Tile_Async(jobz,uplo,A,D,E,T,Q,LDQ,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         integer(kind=c_int), intent(in) :: LDQ
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         real(kind=c_double), intent(out), target :: D(*)
         real(kind=c_double), intent(out), target :: E(*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         complex(kind=c_double_complex), intent(out), target :: Q(LDQ, *)
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zhetrd_Tile_Async_c(jobz,uplo,A,c_loc(D),c_loc(E),T,c_loc(Q),LDQ,sequence,request)
      end subroutine PLASMA_zhetrd_Tile_Async

#if defined(PRECISION_z) || defined(PRECISION_c)
      subroutine PLASMA_zlanhe_Tile_Async(norm,uplo,A,value,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: norm
         integer(kind=c_int), intent(in) :: uplo
         integer(kind=c_int), intent(out), target :: value
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zlanhe_Tile_Async_c(norm,uplo,A,c_loc(value),sequence,request)
      end subroutine PLASMA_zlanhe_Tile_Async
#endif

      subroutine PLASMA_zlaswp_Tile_Async(A,K1,K2,IPIV,INCX,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: INCX
         integer(kind=c_int), intent(in) :: K1
         integer(kind=c_int), intent(in) :: K2
         integer(kind=c_int), intent(in), target :: IPIV(*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zlaswp_Tile_Async_c(A,K1,K2,c_loc(IPIV),INCX,sequence,request)
      end subroutine PLASMA_zlaswp_Tile_Async

      subroutine PLASMA_zlauum_Tile_Async(uplo,A,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zlauum_Tile_Async_c(uplo,A,sequence,request)
      end subroutine PLASMA_zlauum_Tile_Async

#if defined(PRECISION_z) || defined(PRECISION_c)
      subroutine PLASMA_zplghe_Tile_Async(bump,A,seed,sequence,request,info) 
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         real(kind=c_double), intent(in) :: bump
         integer(kind=c_long_long), intent(in) :: seed
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zplghe_Tile_Async_c(bump,A,seed,sequence,request)
      end subroutine PLASMA_zplghe_Tile_Async
#endif

      subroutine PLASMA_zplgsy_Tile_Async(bump,A,seed,sequence,request,info) 
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         complex(kind=c_double_complex), intent(in) :: bump
         integer(kind=c_long_long), intent(in) :: seed
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zplgsy_Tile_Async_c(bump,A,seed,sequence,request)
      end subroutine PLASMA_zplgsy_Tile_Async
       
      subroutine PLASMA_zplrnt_Tile_Async(A,seed,sequence,request,info) 
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_long_long), intent(in) :: seed
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zplrnt_Tile_Async_c(A,seed,sequence,request)
      end subroutine PLASMA_zplrnt_Tile_Async

      subroutine PLASMA_zposv_Tile_Async(uplo,A,B,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zposv_Tile_Async_c(uplo,A,B,sequence,request)
      end subroutine PLASMA_zposv_Tile_Async

      subroutine PLASMA_zpotrf_Tile_Async(uplo,A,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zpotrf_Tile_Async_c(uplo,A,sequence,request)
      end subroutine PLASMA_zpotrf_Tile_Async

      subroutine PLASMA_zpotri_Tile_Async(uplo,A,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zpotri_Tile_Async_c(uplo,A,sequence,request)
      end subroutine PLASMA_zpotri_Tile_Async

      subroutine PLASMA_zpotrs_Tile_Async(uplo,A,B,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zpotrs_Tile_Async_c(uplo,A,B,sequence,request)
      end subroutine PLASMA_zpotrs_Tile_Async

      subroutine PLASMA_zsymm_Tile_Async(side,uplo,alpha,A,B,beta,C,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zsymm_Tile_Async_c(side,uplo,alpha,A,B,beta,C,sequence,request)
      end subroutine PLASMA_zsymm_Tile_Async

      subroutine PLASMA_zsyrk_Tile_Async(uplo,trans,alpha,A,beta,C,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zsyrk_Tile_Async_c(uplo,trans,alpha,A,beta,C,sequence,request)
      end subroutine PLASMA_zsyrk_Tile_Async

      subroutine PLASMA_zsyr2k_Tile_Async(uplo,trans,alpha,A,B,beta,C,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zsyr2k_Tile_Async_c(uplo,trans,alpha,A,B,beta,C,sequence,request)
      end subroutine PLASMA_zsyr2k_Tile_Async

      subroutine PLASMA_ztrmm_Tile_Async(side,uplo,transA,diag,alpha,A,B,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: diag
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: transA
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_ztrmm_Tile_Async_c(side,uplo,transA,diag,alpha,A,B,sequence,request)
      end subroutine PLASMA_ztrmm_Tile_Async

      subroutine PLASMA_ztrsm_Tile_Async(side,uplo,transA,diag,alpha,A,B,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: diag
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: transA
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_ztrsm_Tile_Async_c(side,uplo,transA,diag,alpha,A,B,sequence,request)
      end subroutine PLASMA_ztrsm_Tile_Async

      subroutine PLASMA_ztrsmpl_Tile_Async(A,L,IPIV,B,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: L ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: IPIV ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_ztrsmpl_Tile_Async_c(A,L,IPIV,B,sequence,request)
      end subroutine PLASMA_ztrsmpl_Tile_Async

      subroutine PLASMA_ztrtri_Tile_Async(uplo,diag,A,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: diag
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_ztrtri_Tile_Async_c(uplo,diag,A,sequence,request)
      end subroutine PLASMA_ztrtri_Tile_Async

      subroutine PLASMA_zunglq_Tile_Async(A,T,B,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zunglq_Tile_Async_c(A,T,B,sequence,request)
      end subroutine PLASMA_zunglq_Tile_Async

      subroutine PLASMA_zungqr_Tile_Async(A,T,B,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zungqr_Tile_Async_c(A,T,B,sequence,request)
      end subroutine PLASMA_zungqr_Tile_Async

      subroutine PLASMA_zunmlq_Tile_Async(side,trans,A,T,B,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: trans
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zunmlq_Tile_Async_c(side,trans,A,T,B,sequence,request)
      end subroutine PLASMA_zunmlq_Tile_Async

      subroutine PLASMA_zunmqr_Tile_Async(side,trans,A,T,B,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: trans
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zunmqr_Tile_Async_c(side,trans,A,T,B,sequence,request)
      end subroutine PLASMA_zunmqr_Tile_Async

      subroutine PLASMA_zlacpy_Tile_Async(uplo,A,B,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zlacpy_Tile_Async_c(uplo,A,B,sequence,request)
      end subroutine PLASMA_zlacpy_Tile_Async

      subroutine PLASMA_zlaset_Tile_Async(uplo,alpha,beta,A,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zlaset_Tile_Async_c(uplo,alpha,beta,A,sequence,request)
      end subroutine PLASMA_zlaset_Tile_Async

      subroutine PLASMA_zlaswpc_Tile_Async(A,K1,K2,IPIV,INCX,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in), target :: IPIV(*)
         integer(kind=c_int), intent(in) :: INCX
         integer(kind=c_int), intent(in) :: K1
         integer(kind=c_int), intent(in) :: K2
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zlaswpc_Tile_Async_c(A,K1,K2,c_loc(IPIV),INCX,sequence,request)
      end subroutine PLASMA_zlaswpc_Tile_Async

      subroutine PLASMA_ztrsmrv_Tile_Async(side,uplo,transA,diag,alpha,A,B,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: diag
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: transA
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_ztrsmrv_Tile_Async_c(side,uplo,transA,diag,alpha,A,B,sequence,request)
      end subroutine PLASMA_ztrsmrv_Tile_Async

      subroutine PLASMA_Alloc_Workspace_zgelqf(M,N,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: T ! T is **, so pass by reference
         info = PLASMA_Alloc_Workspace_zgelqf_c(M,N,T)
      end subroutine PLASMA_Alloc_Workspace_zgelqf

      subroutine PLASMA_Alloc_Workspace_zgels(M,N,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: T ! T is **, so pass by reference
         info = PLASMA_Alloc_Workspace_zgels_c(M,N,T)
      end subroutine PLASMA_Alloc_Workspace_zgels

      subroutine PLASMA_Alloc_Workspace_zgeqrf(M,N,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: T ! T is **, so pass by reference
         info = PLASMA_Alloc_Workspace_zgeqrf_c(M,N,T)
      end subroutine PLASMA_Alloc_Workspace_zgeqrf

      subroutine PLASMA_Alloc_Workspace_zgesv_incpiv(N,L,IPIV,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
         type(c_ptr) :: L ! L is **, so pass by reference
         info = PLASMA_Alloc_Workspace_zgesv_incpiv_c(N,L,IPIV)
      end subroutine PLASMA_Alloc_Workspace_zgesv_incpiv

      subroutine PLASMA_Alloc_Workspace_zgetrf_incpiv(M,N,L,IPIV,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
         type(c_ptr) :: L ! L is **, so pass by reference
         info = PLASMA_Alloc_Workspace_zgetrf_incpiv_c(M,N,L,IPIV)
      end subroutine PLASMA_Alloc_Workspace_zgetrf_incpiv

      subroutine PLASMA_Alloc_Workspace_zgeev(N,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: T ! T is **, so pass by reference
         info = PLASMA_Alloc_Workspace_zgeev_c(N,T)
      end subroutine PLASMA_Alloc_Workspace_zgeev

      subroutine PLASMA_Alloc_Workspace_zgebrd(M,N,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: T ! T is **, so pass by reference
         info = PLASMA_Alloc_Workspace_zgebrd_c(M,N,T)
      end subroutine PLASMA_Alloc_Workspace_zgebrd

      subroutine PLASMA_Alloc_Workspace_zgehrd(N,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: T ! T is **, so pass by reference
         info = PLASMA_Alloc_Workspace_zgehrd_c(N,T)
      end subroutine PLASMA_Alloc_Workspace_zgehrd

      subroutine PLASMA_Alloc_Workspace_zgesvd(M,N,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: T ! T is **, so pass by reference
         info = PLASMA_Alloc_Workspace_zgesvd_c(M,N,T)
      end subroutine PLASMA_Alloc_Workspace_zgesvd

      subroutine PLASMA_Alloc_Workspace_zheev(M,N,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: T ! T is **, so pass by reference
         info = PLASMA_Alloc_Workspace_zheev_c(M,N,T)
      end subroutine PLASMA_Alloc_Workspace_zheev

      subroutine PLASMA_Alloc_Workspace_zheevd(M,N,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: T ! T is **, so pass by reference
         info = PLASMA_Alloc_Workspace_zheevd_c(M,N,T)
      end subroutine PLASMA_Alloc_Workspace_zheevd

      subroutine PLASMA_Alloc_Workspace_zhegv(M,N,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: T ! T is **, so pass by reference
         info = PLASMA_Alloc_Workspace_zhegv_c(M,N,T)
      end subroutine PLASMA_Alloc_Workspace_zhegv

      subroutine PLASMA_Alloc_Workspace_zhegvd(M,N,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: T ! T is **, so pass by reference
         info = PLASMA_Alloc_Workspace_zhegvd_c(M,N,T)
      end subroutine PLASMA_Alloc_Workspace_zhegvd

      subroutine PLASMA_Alloc_Workspace_zhetrd(M,N,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: T ! T is **, so pass by reference
         info = PLASMA_Alloc_Workspace_zhetrd_c(M,N,T)
      end subroutine PLASMA_Alloc_Workspace_zhetrd

      subroutine PLASMA_Alloc_Workspace_zgelqf_Tile(M,N,descT,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: descT ! descT is **, so pass by reference
         info = PLASMA_Alloc_Workspace_zgelqf_Tile_c(M,N,descT)
      end subroutine PLASMA_Alloc_Workspace_zgelqf_Tile

      subroutine PLASMA_Alloc_Workspace_zgels_Tile(M,N,descT,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: descT ! descT is **, so pass by reference
         info = PLASMA_Alloc_Workspace_zgels_Tile_c(M,N,descT)
      end subroutine PLASMA_Alloc_Workspace_zgels_Tile

      subroutine PLASMA_Alloc_Workspace_zgeqrf_Tile(M,N,descT,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: descT ! descT is **, so pass by reference
         info = PLASMA_Alloc_Workspace_zgeqrf_Tile_c(M,N,descT)
      end subroutine PLASMA_Alloc_Workspace_zgeqrf_Tile

      subroutine PLASMA_Alloc_Workspace_zgesv_incpiv_Tile(N,descL,IPIV,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
         type(c_ptr) :: descL ! descL is **, so pass by reference
         info = PLASMA_Alloc_Workspace_zgesv_incpiv_Tile_c(N,descL,IPIV)
      end subroutine PLASMA_Alloc_Workspace_zgesv_incpiv_Tile

      subroutine PLASMA_Alloc_Workspace_zgetrf_incpiv_Tile(N,descL,IPIV,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
         type(c_ptr) :: descL ! descL is **, so pass by reference
         info = PLASMA_Alloc_Workspace_zgetrf_incpiv_Tile_c(N,descL,IPIV)
      end subroutine PLASMA_Alloc_Workspace_zgetrf_incpiv_Tile

      subroutine PLASMA_Alloc_Workspace_zgetri_Tile_Async(A,W,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: W ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_Alloc_Workspace_zgetri_Tile_Async_c(A,W)
      end subroutine PLASMA_Alloc_Workspace_zgetri_Tile_Async

      subroutine PLASMA_zLapack_to_Tile(Af77,LDA,A,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         complex(kind=c_double_complex), intent(in), target :: Af77(LDA,*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zLapack_to_Tile_c(c_loc(Af77),LDA,A)
      end subroutine PLASMA_zLapack_to_Tile

      subroutine PLASMA_zTile_to_Lapack(A,Af77,LDA,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         complex(kind=c_double_complex), intent(out), target :: Af77(LDA,*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zTile_to_Lapack_c(A,c_loc(Af77),LDA)
      end subroutine PLASMA_zTile_to_Lapack

      subroutine PLASMA_zLapack_to_Tile_Async(Af77,LDA,A,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         complex(kind=c_double_complex), intent(in), target :: Af77(LDA,*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zLapack_to_Tile_Async_c(c_loc(Af77),LDA,A,sequence,request)
      end subroutine PLASMA_zLapack_to_Tile_Async

      subroutine PLASMA_zTile_to_Lapack_Async(A,Af77,LDA,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         complex(kind=c_double_complex), intent(out), target :: Af77(LDA,*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_zTile_to_Lapack_Async_c(A,c_loc(Af77),LDA,sequence,request)
      end subroutine PLASMA_zTile_to_Lapack_Async

end module plasma_z
