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
! @file example_dgesv_tile_async_f90.f90
!
!  PLASMA Fortran 90 example using Fortran 2003 ISO C bindings
!  PLASMA is a software package provided by Univ. of Tennessee,
!  Univ. of California Berkeley and Univ. of Colorado Denver
!
! @version 2.4.6
! @author Numerical Algorithms Group
! @date 2012-08-14
!
! ./example_dgesv_tile_async_f90 < gesv_example.d
!

    program example_dgesv_tile_async_f90
       use plasma
       use iso_c_binding
       implicit none
       integer, parameter              :: wp = kind(0.0d0)
       integer, parameter              :: nin = 5, nout = 6
       integer                         :: info, lda, m, n, nrhs, ldb
       real(kind=wp), allocatable      :: a_lpk(:,:), a_tile(:,:), a(:,:), b(:,:), b_tile(:,:), b_lpk(:,:)
       integer, allocatable            :: ipiv(:)
       character                       :: sched
       integer                         :: npcores, nb, ib
       integer, target                 :: request_tmp
       logical                         :: atun
       intrinsic                          min
       intrinsic                       :: random_number
       type(c_ptr)                     :: a_desc, sequence, request, b_desc

       ! Read and initialise problem
       read (nin,*) m, nrhs
       n = m
       read (nin,*) atun
       read (nin,*) nb, ib
       read (nin,*) npcores
       read (nin,*) sched
       lda = m
       ldb = m
       allocate (a_tile(lda,n),ipiv(min(m,n)),a(lda,n),b_tile(ldb,nrhs),b(ldb,nrhs))
       allocate(a_lpk(lda,n),b_lpk(ldb,nrhs))
       call random_number(a)
       call random_number(b)

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

       ! Create a descriptor for the array A. Note that a_desc is type(c_ptr)
       call plasma_desc_create(a_desc,a_tile,PlasmaRealDouble,nb,nb,nb*nb,lda,n,0,0,m,n,info)
       if (info/=plasma_success) stop 'plasma_desc_create failed'
       ! Convert array into tiled layout
       call plasma_lapack_to_tile(a,lda,a_desc,info)
       if (info/=plasma_success) stop 'plasma_lapack_to_tile failed'
       write(nout,*) "Sanity check in: ",a(1:min(5,m),1)

       ! Create a descriptor for the array B. Note that b_desc is type(c_ptr)
       call plasma_desc_create(b_desc,b_tile,PlasmaRealDouble,nb,nb,nb*nb,ldb,nrhs,0,0,m,nrhs,info)
       if (info/=plasma_success) stop 'plasma_desc_create failed'
       ! Convert array into tiled layout
       call plasma_lapack_to_tile(b,ldb,b_desc,info)
       if (info/=plasma_success) stop 'plasma_lapack_to_tile failed'

       ! Create a sequence handle. Sequence is type(c_ptr).
       call plasma_sequence_create(sequence,info)
       ! Create a place-holder for the intent-out request.
       request_tmp = 0
       request = c_loc(request_tmp)

       ! Overlap asynchronous factorization and solve calls.
       call plasma_dgetrf_tile_async(a_desc,ipiv,sequence,request,info)
       call plasma_dgetrs_tile_async(PLASMANOTRANS,a_desc,ipiv,b_desc,sequence,request,info)
       if (info/=0) write (nout,*) 'plasma call failed'

       ! Wait for the async calls to finish.
       call plasma_sequence_wait(sequence,info)
       if (info/=0) write (nout,*) 'plasma wait call failed'

       ! Convert back to column layout
       call plasma_tile_to_lapack(a_desc,a_lpk,lda,info)
       call plasma_tile_to_lapack(b_desc,b_lpk,ldb,info)
       if (info/=plasma_success) stop 'plasma_tile_to_lapack failed'
       write(nout,*) "Sanity check out: ",a_lpk(1:min(5,m),1), ipiv(1:min(5,m))
       write(nout,*) "Sanity check out: ",b_lpk(1:min(5,m),1)

       ! Close down
       call plasma_finalize(info)

    end program example_dgesv_tile_async_f90
