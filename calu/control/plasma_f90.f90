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
! @file plasma_f90.f90
!
!  PLASMA fortran 90 interface
!  PLASMA is a software package provided by Univ. of Tennessee,
!  Univ. of California Berkeley and Univ. of Colorado Denver
!
! @version 2.4.2
! @author Numerical Algorithm Group
! @date 2011-09-15
! @precisions normal z -> c d s
!
module plasma

      use plasma_s
      use plasma_d
      use plasma_ds
      use plasma_c
      use plasma_z
      use plasma_zc
      include 'plasmaf.h'

      logical :: plasma_initialized = .false.
      integer, parameter :: sp = kind(0.0)
      integer, parameter :: dp = kind(0.0d0)

      interface
         function PLASMA_Init_c(cores) &
          & bind(c, name='PLASMA_Init')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Init_c
            integer(kind=c_int), value :: cores
         end function PLASMA_Init_c
      end interface

      interface
         function PLASMA_Finalize_c() &
          & bind(c, name='PLASMA_Finalize')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Finalize_c
         end function PLASMA_Finalize_c
      end interface

      interface
         function PLASMA_Set_c(param, pval) &
          & bind(c, name='PLASMA_Set')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Set_c
            integer(kind=c_int), value :: param
            integer(kind=c_int), value :: pval
         end function PLASMA_Set_c
      end interface

      interface
         function PLASMA_Get_c(param, pval) &
          & bind(c, name='PLASMA_Get')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Get_c
            integer(kind=c_int), value :: param
            type(c_ptr), value :: pval
         end function PLASMA_Get_c
      end interface

      interface
         function PLASMA_Enable_c(param) &
          & bind(c, name='PLASMA_Enable')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Enable_c
            integer(kind=c_int), value :: param
         end function PLASMA_Enable_c
      end interface

      interface
         function PLASMA_Disable_c(param) &
          & bind(c, name='PLASMA_Disable')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Disable_c
            integer(kind=c_int), value :: param
         end function PLASMA_Disable_c
      end interface

      interface
         function PLASMA_Lapack_to_Tile_c(a_lpk,lda,a_pma) &
          & bind(c, name='PLASMA_Lapack_to_Tile')
            use iso_c_binding
            integer(kind=c_int) :: PLASMA_Lapack_to_Tile_c
            type(c_ptr), value :: a_lpk, a_pma
            integer(kind=c_int), value :: lda
         end function PLASMA_Lapack_to_Tile_c
      end interface

      interface
         function PLASMA_Tile_to_Lapack_c(a_pma,a_lpk,lda) &
          & bind(c, name='PLASMA_Tile_to_Lapack')
            use iso_c_binding
            integer(kind=c_int) :: PLASMA_Tile_to_Lapack_c
            type(c_ptr), value :: a_lpk, a_pma
            integer(kind=c_int), value :: lda
         end function PLASMA_Tile_to_Lapack_c
      end interface

      interface
         function PLASMA_Desc_Create_c(desc, mat, dtyp, mb, nb, bsiz, lm, ln, i, j, m, n) &
          & bind(c, name='PLASMA_Desc_Create')
            use iso_c_binding
            integer(kind=c_int) :: PLASMA_Desc_Create_c
            type(c_ptr) :: desc
            type(c_ptr), value :: mat
            integer(kind=c_int), value :: dtyp
            integer(kind=c_int), value :: mb, nb, bsiz, lm, ln, i, j, m, n
         end function PLASMA_Desc_Create_c
      end interface

      interface
         function PLASMA_Desc_Destroy_c(desc) &
          & bind(c, name='PLASMA_Desc_Destroy')
            use iso_c_binding
            integer(kind=c_int) :: PLASMA_Desc_Destroy_c
            type(c_ptr) :: desc
         end function PLASMA_Desc_Destroy_c
      end interface

      interface
         subroutine free_c(ptr) bind(c, name='free')
            use iso_c_binding
            implicit none
            type(c_ptr), value :: ptr
         end subroutine free_c
      end interface

      interface
         function PLASMA_Version_c(maj,min,mic) &
          & bind(c, name='PLASMA_Version')
            use iso_c_binding
            integer(kind=c_int) :: PLASMA_Version_c
            type(c_ptr), value ::  maj,min,mic
         end function PLASMA_Version_c
      end interface

      interface
         function PLASMA_Init_Affinity_c(cores,bindtab) &
          & bind(c, name='PLASMA_Init_Affinity')
            use iso_c_binding
            integer(kind=c_int) :: PLASMA_Init_Affinity_c
            integer(kind=c_int), value ::  cores
            type(c_ptr), value :: bindtab
         end function PLASMA_Init_Affinity_c
      end interface

      interface
         function PLASMA_Dealloc_Handle_c(handle) &
          & bind(c, name='PLASMA_Dealloc_Handle')
            use iso_c_binding
            integer(kind=c_int) :: PLASMA_Dealloc_Handle_c
            type(c_ptr) :: handle
         end function PLASMA_Dealloc_Handle_c
      end interface

      interface
         function PLASMA_Dealloc_Handle_Tile_c(desc) &
          & bind(c, name='PLASMA_Dealloc_Handle_Tile')
            use iso_c_binding
            integer(kind=c_int) :: PLASMA_Dealloc_Handle_Tile_c
            type(c_ptr) :: desc
         end function PLASMA_Dealloc_Handle_Tile_c
      end interface

      interface
         function PLASMA_Sequence_Create_c(seq) &
          & bind(c, name='PLASMA_Sequence_Create')
            use iso_c_binding
            integer(kind=c_int) :: PLASMA_Sequence_Create_c
            type(c_ptr) :: seq
         end function PLASMA_Sequence_Create_c
      end interface

      interface
         function PLASMA_Sequence_Destroy_c(seq) &
          & bind(c, name='PLASMA_Sequence_Destroy')
            use iso_c_binding
            integer(kind=c_int) :: PLASMA_Sequence_Destroy_c
            type(c_ptr), value :: seq
         end function PLASMA_Sequence_Destroy_c
      end interface

      interface
         function PLASMA_Sequence_Wait_c(seq) &
          & bind(c, name='PLASMA_Sequence_Wait')
            use iso_c_binding
            integer(kind=c_int) :: PLASMA_Sequence_Wait_c
            type(c_ptr), value :: seq
         end function PLASMA_Sequence_Wait_c
      end interface

      interface
         function PLASMA_Sequence_Flush_c(seq,req) &
          & bind(c, name='PLASMA_Sequence_Flush')
            use iso_c_binding
            integer(kind=c_int) :: PLASMA_Sequence_Flush_c
            type(c_ptr), value :: seq
            type(c_ptr), value :: req
         end function PLASMA_Sequence_Flush_c
      end interface

      interface plasma_lapack_to_tile
         module procedure plasma_lapack_to_tile_s
         module procedure plasma_lapack_to_tile_d
         module procedure plasma_lapack_to_tile_cpx
         module procedure plasma_lapack_to_tile_z
      end interface plasma_lapack_to_tile

      interface plasma_tile_to_lapack
         module procedure plasma_tile_to_lapack_s
         module procedure plasma_tile_to_lapack_d
         module procedure plasma_tile_to_lapack_cpx
         module procedure plasma_tile_to_lapack_z
      end interface plasma_tile_to_lapack

      interface plasma_desc_create
         module procedure plasma_desc_create_s
         module procedure plasma_desc_create_d
         module procedure plasma_desc_create_cpx
         module procedure plasma_desc_create_z
      end interface plasma_desc_create

   contains

   subroutine plasma_init(ncores,info)
      use iso_c_binding
      implicit none
      integer(kind=c_int), intent(in) :: ncores
      integer(kind=c_int), intent(out) :: info
      info = plasma_init_c(ncores)
      plasma_initialized = .true.
   end subroutine plasma_init

   subroutine plasma_finalize(info)
      use iso_c_binding
      implicit none
      integer(kind=c_int), intent(out) :: info
      info = plasma_finalize_c()
      plasma_initialized = .false.
   end subroutine plasma_finalize

   subroutine plasma_set(param,pval,info)
      use iso_c_binding
      implicit none
      integer(kind=c_int), intent(in) :: param
      integer(kind=c_int), intent(in) :: pval
      integer(kind=c_int), intent(out) :: info
      info = plasma_set_c(param,pval)
   end subroutine plasma_set

   subroutine plasma_get(param,pval,info)
      use iso_c_binding
      implicit none
      integer(kind=c_int), intent(in) :: param
      integer(kind=c_int), intent(out), target :: pval
      integer(kind=c_int), intent(out) :: info
      info = plasma_get_c(param,c_loc(pval))
   end subroutine plasma_get

   subroutine plasma_enable(param,info)
      use iso_c_binding
      implicit none
      integer(kind=c_int), intent(in) :: param
      integer(kind=c_int), intent(out) :: info
      info = plasma_enable_c(param)
   end subroutine plasma_enable

   subroutine plasma_disable(param,info)
      use iso_c_binding
      implicit none
      integer(kind=c_int), intent(in) :: param
      integer(kind=c_int), intent(out) :: info
      info = plasma_disable_c(param)
   end subroutine plasma_disable

! overloaded: single precision
   subroutine plasma_lapack_to_tile_s(a_lpk,lda,a_pma,info)
      use iso_c_binding
      implicit none
      integer(kind=c_int), intent(in) :: lda
      real(kind=sp), intent(out), target :: a_lpk(lda,*)
      type(c_ptr), intent(out) ::  a_pma
      integer(kind=c_int), intent(out) :: info
      info = plasma_lapack_to_tile_c(c_loc(a_lpk),lda,a_pma)
   end subroutine plasma_lapack_to_tile_s
! overloaded: double precision
   subroutine plasma_lapack_to_tile_d(a_lpk,lda,a_pma,info)
      use iso_c_binding
      implicit none
      integer(kind=c_int), intent(in) :: lda
      real(kind=dp), intent(out), target :: a_lpk(lda,*)
      type(c_ptr), intent(out) ::  a_pma
      integer(kind=c_int), intent(out) :: info
      info = plasma_lapack_to_tile_c(c_loc(a_lpk),lda,a_pma)
   end subroutine plasma_lapack_to_tile_d
! overloaded: single precision complex
   subroutine plasma_lapack_to_tile_cpx(a_lpk,lda,a_pma,info)
      use iso_c_binding
      implicit none
      integer(kind=c_int), intent(in) :: lda
      complex(kind=sp), intent(out), target :: a_lpk(lda,*)
      type(c_ptr), intent(out) ::  a_pma
      integer(kind=c_int), intent(out) :: info
      info = plasma_lapack_to_tile_c(c_loc(a_lpk),lda,a_pma)
   end subroutine plasma_lapack_to_tile_cpx
! overloaded: double precision complex
   subroutine plasma_lapack_to_tile_z(a_lpk,lda,a_pma,info)
      use iso_c_binding
      implicit none
      integer(kind=c_int), intent(in) :: lda
      complex(kind=dp), intent(out), target :: a_lpk(lda,*)
      type(c_ptr), intent(out) ::  a_pma
      integer(kind=c_int), intent(out) :: info
      info = plasma_lapack_to_tile_c(c_loc(a_lpk),lda,a_pma)
   end subroutine plasma_lapack_to_tile_z

! overloaded: single precision
   subroutine plasma_tile_to_lapack_s(a_pma,a_lpk,lda,info)
      use iso_c_binding
      implicit none
      integer(kind=c_int), intent(in) :: lda
      real(kind=sp), intent(out), target :: a_lpk(lda,*)
      type(c_ptr), intent(in) ::  a_pma
      integer(kind=c_int), intent(out) :: info
      info = plasma_tile_to_lapack_c(a_pma,c_loc(a_lpk),lda)
   end subroutine plasma_tile_to_lapack_s
! overloaded: double precision
   subroutine plasma_tile_to_lapack_d(a_pma,a_lpk,lda,info)
      use iso_c_binding
      implicit none
      integer(kind=c_int), intent(in) :: lda
      real(kind=dp), intent(out), target :: a_lpk(lda,*)
      type(c_ptr), intent(in) ::  a_pma
      integer(kind=c_int), intent(out) :: info
      info = plasma_tile_to_lapack_c(a_pma,c_loc(a_lpk),lda)
   end subroutine plasma_tile_to_lapack_d
! overloaded: single precision complex
   subroutine plasma_tile_to_lapack_cpx(a_pma,a_lpk,lda,info)
      use iso_c_binding
      implicit none
      integer(kind=c_int), intent(in) :: lda
      complex(kind=sp), intent(out), target :: a_lpk(lda,*)
      type(c_ptr), intent(in) ::  a_pma
      integer(kind=c_int), intent(out) :: info
      info = plasma_tile_to_lapack_c(a_pma,c_loc(a_lpk),lda)
   end subroutine plasma_tile_to_lapack_cpx
! overloaded: double precision complex
   subroutine plasma_tile_to_lapack_z(a_pma,a_lpk,lda,info)
      use iso_c_binding
      implicit none
      integer(kind=c_int), intent(in) :: lda
      complex(kind=dp), intent(out), target :: a_lpk(lda,*)
      type(c_ptr), intent(in) ::  a_pma
      integer(kind=c_int), intent(out) :: info
      info = plasma_tile_to_lapack_c(a_pma,c_loc(a_lpk),lda)
   end subroutine plasma_tile_to_lapack_z

! overloaded: single precision
   subroutine plasma_desc_create_s(desc,mat,dtyp,mb,nb,bsiz,lm,ln,i,j,m,n,info)
      use iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: desc
      integer(kind=c_int), intent(in) :: mb, nb, bsiz, lm, ln, i, j, m, n
      real(kind=sp), intent(in), target :: mat(lm,*)
      integer(kind=c_int), intent(in) :: dtyp
      integer(kind=c_int), intent(out) :: info
      info = plasma_desc_create_c(desc,c_loc(mat),dtyp,mb,nb,bsiz,lm,ln,i,j,m,n)
   end subroutine plasma_desc_create_s
! overloaded: double precision
   subroutine plasma_desc_create_d(desc,mat,dtyp,mb,nb,bsiz,lm,ln,i,j,m,n,info)
      use iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: desc
      integer(kind=c_int), intent(in) :: mb, nb, bsiz, lm, ln, i, j, m, n
      real(kind=dp), intent(in), target :: mat(lm,*)
      integer(kind=c_int), intent(in) :: dtyp
      integer(kind=c_int), intent(out) :: info
      info = plasma_desc_create_c(desc,c_loc(mat),dtyp,mb,nb,bsiz,lm,ln,i,j,m,n)
   end subroutine plasma_desc_create_d
! overloaded: single precision complex
   subroutine plasma_desc_create_cpx(desc,mat,dtyp,mb,nb,bsiz,lm,ln,i,j,m,n,info)
      use iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: desc
      integer(kind=c_int), intent(in) :: mb, nb, bsiz, lm, ln, i, j, m, n
      complex(kind=sp), intent(in), target :: mat(lm,*)
      integer(kind=c_int), intent(in) :: dtyp
      integer(kind=c_int), intent(out) :: info
      info = plasma_desc_create_c(desc,c_loc(mat),dtyp,mb,nb,bsiz,lm,ln,i,j,m,n)
   end subroutine plasma_desc_create_cpx
! overloaded: double precision complex
   subroutine plasma_desc_create_z(desc,mat,dtyp,mb,nb,bsiz,lm,ln,i,j,m,n,info)
      use iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: desc
      integer(kind=c_int), intent(in) :: mb, nb, bsiz, lm, ln, i, j, m, n
      complex(kind=dp), intent(in), target :: mat(lm,*)
      integer(kind=c_int), intent(in) :: dtyp
      integer(kind=c_int), intent(out) :: info
      info = plasma_desc_create_c(desc,c_loc(mat),dtyp,mb,nb,bsiz,lm,ln,i,j,m,n)
   end subroutine plasma_desc_create_z

   subroutine plasma_desc_destroy(desc,info)
      use iso_c_binding
      implicit none
      type(c_ptr), intent(inout) :: desc
      integer(kind=c_int), intent(out) :: info
      info = plasma_desc_destroy_c(desc)
   end subroutine plasma_desc_destroy

   subroutine plasma_free(ptr)
      use iso_c_binding
      implicit none
      type(c_ptr), intent(in) :: ptr
      call free_c(ptr)
   end subroutine plasma_free

   subroutine plasma_version(ver_major,ver_minor,ver_micro,info)
      use iso_c_binding
      implicit none
      integer(kind=c_int), intent(out), target :: ver_major,ver_minor,ver_micro
      integer(kind=c_int), intent(out) :: info
      info = plasma_version_c(c_loc(ver_major),c_loc(ver_minor),c_loc(ver_micro))
   end subroutine plasma_version

   subroutine plasma_init_affinity(cores,bindtab,info)
      use iso_c_binding
      implicit none
      integer(kind=c_int), intent(in) :: cores
      integer(kind=c_int), intent(out), target :: bindtab 
      integer(kind=c_int), intent(out) :: info
      info = plasma_init_affinity_c(cores,c_loc(bindtab))
   end subroutine plasma_init_affinity

   subroutine plasma_dealloc_handle(handle,info)
      use iso_c_binding
      implicit none
      type(c_ptr), intent(inout) :: handle
      integer(kind=c_int), intent(out) :: info
      info = plasma_dealloc_handle_c(handle)
   end subroutine plasma_dealloc_handle

   subroutine plasma_dealloc_handle_tile(desc,info)
      use iso_c_binding
      implicit none
      type(c_ptr), intent(inout) :: desc
      integer(kind=c_int), intent(out) :: info
      info = plasma_dealloc_handle_tile_c(desc)
   end subroutine plasma_dealloc_handle_tile

   subroutine plasma_sequence_create(sequence,info)
      use iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: sequence
      integer(kind=c_int), intent(out) :: info
      info = plasma_sequence_create_c(sequence)
   end subroutine plasma_sequence_create

   subroutine plasma_sequence_destroy(sequence,info)
      use iso_c_binding
      implicit none
      type(c_ptr), intent(in) :: sequence
      integer(kind=c_int), intent(out) :: info
      info = plasma_sequence_destroy_c(sequence)
   end subroutine plasma_sequence_destroy

   subroutine plasma_sequence_wait(sequence,info)
      use iso_c_binding
      implicit none
      type(c_ptr), intent(in) :: sequence
      integer(kind=c_int), intent(out) :: info
      info = plasma_sequence_wait_c(sequence)
   end subroutine plasma_sequence_wait

   subroutine plasma_sequence_flush(sequence,request,info)
      use iso_c_binding
      implicit none
      type(c_ptr), intent(in) :: sequence, request
      integer(kind=c_int), intent(out) :: info
      info = plasma_sequence_flush_c(sequence,request)
   end subroutine plasma_sequence_flush

end module plasma
