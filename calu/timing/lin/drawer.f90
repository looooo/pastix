      SUBROUTINE DRAWER( NAMEIO, TYPE, LAB1, LAB2, LAB3, NK, KVAL, NN,  &
     &                    NVAL, NLDA, RESLTS, LDR1, LDR2 ) 
      USE GNUFOR2
      IMPLICIT NONE
!*
!*     .. Scalar Arguments ..
      CHARACTER*( * )  :: NAMEIO, TYPE, LAB1, LAB2, LAB3
      INTEGER          :: LDR1, LDR2, NK, NLDA, NN
!*     ..
!*     .. Array Arguments ..
      INTEGER          :: KVAL( NK ), NVAL( NN )
      DOUBLE PRECISION :: RESLTS( LDR1, LDR2, * )
      INTEGER          :: i,j,k, N0
      CHARACTER( len=50 )  :: NEWNAMEIO, NEWTITLE
      CHARACTER( len=50 )  :: NUM

      Write (Unit=NUM, FMT="(I5.1)") KVAL(NK)
      NEWNAMEIO = NAMEIO //'.'// TYPE
      NEWTITLE = NAMEIO //' ('//LAB1//' = '//trim(NUM)//')'

      CALL PLOT( dble( NVAL ), RESLTS( NK, 1:NN, 1 ), '22-',            &
     &           terminal = TYPE, filename = NEWNAMEIO,                 &
     &           title = NEWTITLE, xlabel = LAB2, ylabel = LAB3 )

      END SUBROUTINE DRAWER
