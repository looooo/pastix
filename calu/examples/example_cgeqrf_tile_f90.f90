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
! @file example_cgeqrf_tile_f90.f90
!
!  PLASMA Fortran 90 example using Fortran 2003 ISO C bindings
!  PLASMA is a software package provided by Univ. of Tennessee,
!  Univ. of California Berkeley and Univ. of Colorado Denver
!
! @version 2.4.6
! @author Numerical Algorithms Group
! @date 2012-08-14
!
! ./example_cgeqrf_tile_f90 < geqrf_example.d
!
    program example_cgeqrf_tile_f90
       use plasma
       use iso_c_binding
       implicit none
       integer, parameter              :: wp = kind(0.0)
       integer, parameter              :: nin = 5, nout = 6
       integer                         :: info, lda, m, n
       complex(kind=wp), allocatable   :: a_tile(:,:), a_lpk(:,:), a(:,:)
       real(kind=wp), allocatable      :: a_r(:,:), a_i(:,:)
       character                       :: sched
       integer                         :: npcores, nb, ib
       logical                         :: atun
       intrinsic                       :: min
       intrinsic                       :: random_number
       type(c_ptr)                     :: desc_t, desc_a

       ! Read and initialise problem
       read (nin,*) m, n
       read (nin,*) atun
       read (nin,*) nb, ib
       read (nin,*) npcores
       read (nin,*) sched
       lda = m
       allocate (a_tile(lda,n),a(lda,n),a_r(lda,n),a_i(lda,n))
       allocate(a_lpk(lda,n))
       call random_number(a_r)
       call random_number(a_i)
       a(:,:) = cmplx(a_r(:,:),a_i(:,:))
       deallocate(a_r,a_i)

       ! Initialize PLASMA
       call plasma_init(npcores,info)

       ! Choose whether to perform autotuning
       if (.not. atun) then
          call plasma_disable(plasma_autotuning,info)
          call plasma_set(plasma_tile_size,nb,info)
          call plasma_set(plasma_inner_block_size,ib,info)
       else
          call plasma_enable(plasma_autotuning,info)
       end if
       ! Choose whether to use dynamic scheduling
       if (sched=='d') then
          call plasma_set(plasma_scheduling_mode,plasma_static_scheduling,info)
       end if

       ! Create workspace point to be desc_t, which is type(c_ptr)
       call plasma_alloc_workspace_cgeqrf_tile(m,n,desc_t,info)

       ! Create a descriptor for the array A. Nte that desc_a is type(c_ptr)
       call plasma_desc_create(desc_a,a_tile,PlasmaComplexFloat,nb,nb,nb*nb,lda,n,0,0,m,n,info)
       if (info/=plasma_success) stop 'plasma_desc_create failed'

       ! Convert array into tiled layout
       call plasma_lapack_to_tile(a,lda,desc_a,info)
       if (info/=plasma_success) stop 'plasma_lapack_to_tile failed'

       write(nout,*) "Sanity check in: ",a(1,1:min(5,n))
         
       ! Call tiled routine passing in descriptor for A
       call plasma_cgeqrf_tile(desc_a,desc_t,info)
       if (info/=0) write (nout,*) 'plasma call failed'

       ! Convert back into column layout
       call plasma_tile_to_lapack(desc_a,a_lpk,lda,info)
       if (info/=plasma_success) stop 'plasma_tile_to_lapack failed'
       write(nout,*) "Sanity check out: ",a_lpk(1,1:min(5,n))

       ! Close down
       call plasma_finalize(info)

    end program example_cgeqrf_tile_f90
