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
! @file example_zgetrf_tile_f90.f90
!
!  PLASMA Fortran 90 example using Fortran 2003 ISO C bindings
!  PLASMA is a software package provided by Univ. of Tennessee,
!  Univ. of California Berkeley and Univ. of Colorado Denver
!
! @version 2.4.6
! @author Numerical Algorithms Group
! @date 2012-08-14
!
! ./example_zgetrf_tile_f90 < getrf_example.d
!

    program example_zgetrf_tile_f90
       use plasma
       use iso_c_binding
       implicit none
       integer, parameter              :: wp = kind(0.0d0)
       integer, parameter              :: nin = 5, nout = 6
       integer                         :: info, lda, m, n
       complex(kind=wp), allocatable   :: a_tile(:,:), a(:,:), a_lpk(:,:)
       real(kind=wp), allocatable      :: a_r(:,:), a_i(:,:)
       integer, allocatable            :: ipiv(:)
       character                       :: sched
       integer                         :: npcores, nb, ib
       logical                         :: atun
       intrinsic                       :: min
       intrinsic                       :: random_number
       type(c_ptr)                     :: a_desc

       ! Read and initialise problem
       write(*,*) 'Enter the dimensions M N:'
       read (nin,*) m, n
       write(*,*) 'Do you want autotuning ? [.false./.true.]:'
       read (nin,*) atun
       write(*,*) 'Enter NB IB:'
       read (nin,*) nb, ib
       write(*,*) 'Enter the number of cores:'
       read (nin,*) npcores
       write(*,*) 'Static/dynamic scheduler? [s/d]:'
       read (nin,*) sched
       lda = m
       allocate (a_tile(lda,n),ipiv(min(m,n)),a(lda,n),a_lpk(lda,n))

       allocate (a_r(lda,n),a_i(lda,n))
       call random_number(a_r)
       call random_number(a_i)
       a(:,:) = cmplx(a_r(:,:),a_i(:,:),kind=wp)
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

       ! Create a descriptor for the array A. Note that a_desc is type(c_ptr)
       call plasma_desc_create(a_desc, & ! Descriptor
            &                  a_tile, & ! Array to store the matrix
            &                  PlasmaComplexDouble, & ! Precision
            &                  nb, nb, nb*nb, & ! Tile dimensions
            &                  lda, n, & ! Entire matrix dimensions
            &                  0, 0,   & ! Row and column indexes of the beginning of the submatrix
            &                  m, n,   & ! Submatrix dimensions
            &                  info)
       if (info/=plasma_success) stop 'plasma_desc_create failed'

       ! Convert array into tiled layout
       call plasma_lapack_to_tile(a,lda,a_desc,info)
       if (info/=plasma_success) stop 'plasma_lapack_to_tile failed'

       write(nout,*) "Sanity check in: ",a(1:min(5,m),1)

       ! Call computational routine passing in descriptor for A
       call plasma_zgetrf_tile(a_desc,ipiv,info)
       if (info/=0) write (nout,*) 'plasma call failed'

       ! Convert back into column layout
       call plasma_tile_to_lapack(a_desc,a_lpk,lda,info)
       if (info/=plasma_success) stop 'plasma_tile_to_lapack failed'
       write(nout,*) "Sanity check out: ",a_lpk(1:min(5,m),1), ipiv(1:min(5,m))

       ! Close down
       call plasma_finalize(info)

    end program example_zgetrf_tile_f90
