      SUBROUTINE ZCHKGE( DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NNS,
     $                   IBVAL, NSVAL, THRESH, TSTERR, NMAX, A, AFAC,
     $                   AINV, B, X, XACT, WORK, RWORK, IWORK, 
     $                   NOUT )
*
      INCLUDE 'plasmaf.h'
*
*  -- LAPACK test routine (version 3.1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     January 2007
*
*     .. Scalar Arguments ..
      LOGICAL            TSTERR
      INTEGER            NM, NMAX, NN, NNB, NNS, NOUT
      DOUBLE PRECISION   THRESH
*     ..
*     .. Array Arguments ..
      LOGICAL            DOTYPE( * )
      INTEGER            IWORK( * ), MVAL( * ), NBVAL( * ), NSVAL( * ),
     $                   IBVAL( * ), NVAL( * )
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         A( * ), AFAC( * ), AINV( * ), B( * ),
     $                   WORK( * ), X( * ), XACT( * )
*     ..
*
*  Purpose
*  =======
*
*  ZCHKGE tests ZGETRF, -TRI, -TRS, -RFS, and -CON.
*
*  Arguments
*  =========
*
*  DOTYPE  (input) LOGICAL array, dimension (NTYPES)
*          The matrix types to be used for testing.  Matrices of type j
*          (for 1 <= j <= NTYPES) are used for testing if DOTYPE(j) =
*          .TRUE.; if DOTYPE(j) = .FALSE., then type j is not used.
*
*  NM      (input) INTEGER
*          The number of values of M contained in the vector MVAL.
*
*  MVAL    (input) INTEGER array, dimension (NM)
*          The values of the matrix row dimension M.
*
*  NN      (input) INTEGER
*          The number of values of N contained in the vector NVAL.
*
*  NVAL    (input) INTEGER array, dimension (NN)
*          The values of the matrix column dimension N.
*
*  NNB     (input) INTEGER
*          The number of values of NB contained in the vector NBVAL.
*
*  NBVAL   (input) INTEGER array, dimension (NBVAL)
*          The values of the blocksize NB.
*
*  IBVAL   (input) INTEGER array, dimension (NBVAL)
*          The values of the inner block size IB.
*
*  NNS     (input) INTEGER
*          The number of values of NRHS contained in the vector NSVAL.
*
*  NSVAL   (input) INTEGER array, dimension (NNS)
*          The values of the number of right hand sides NRHS.
*
*  NRHS    (input) INTEGER
*          The number of right hand side vectors to be generated for
*          each linear system.
*
*  THRESH  (input) DOUBLE PRECISION
*          The threshold value for the test ratios.  A result is
*          included in the output file if RESULT >= THRESH.  To have
*          every test ratio printed, use THRESH = 0.
*
*  TSTERR  (input) LOGICAL
*          Flag that indicates whether error exits are to be tested.
*
*  NMAX    (input) INTEGER
*          The maximum value permitted for M or N, used in dimensioning
*          the work arrays.
*
*  A       (workspace) COMPLEX*16 array, dimension (NMAX*NMAX)
*
*  AFAC    (workspace) COMPLEX*16 array, dimension (NMAX*NMAX)
*
*  AINV    (workspace) COMPLEX*16 array, dimension (NMAX*NMAX)
*
*  B       (workspace) COMPLEX*16 array, dimension (NMAX*NSMAX)
*          where NSMAX is the largest entry in NSVAL.
*
*  X       (workspace) COMPLEX*16 array, dimension (NMAX*NSMAX)
*
*  XACT    (workspace) COMPLEX*16 array, dimension (NMAX*NSMAX)
*
*  WORK    (workspace) COMPLEX*16 array, dimension
*                      (NMAX*max(3,NSMAX))
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension
*                      (max(2*NMAX,2*NSMAX+NWORK))
*
*  IWORK   (workspace) INTEGER array, dimension (NMAX)
*
*  NOUT    (input) INTEGER
*          The unit number for output.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      INTEGER            NTYPES
      PARAMETER          ( NTYPES = 11 )
      INTEGER            NTESTS
      PARAMETER          ( NTESTS = 8 )
      INTEGER            NTRAN
      PARAMETER          ( NTRAN = 1 )
*     ..
*     .. Local Scalars ..
      LOGICAL            TRFCON, ZEROT
      CHARACTER          DIST, NORM, TRANS, TYPE, XTYPE
      CHARACTER*3        PATH
      INTEGER            I, IM, IMAT, IN, INB, INFO, IOFF, IRHS, ITRAN,
     $                   IZERO, K, KL, KU, LDA, LWORK, M, MODE, N, NB,
     $                   NERRS, NFAIL, NIMAT, NRHS, NRUN, NT, IB,
     $                   PLASMA_TRANS
      DOUBLE PRECISION   AINVNM, ANORM, ANORMI, ANORMO, CNDNUM, DUMMY,
     $                   RCOND, RCONDC, RCONDI, RCONDO
      INTEGER            HL( 2 ), HPIV( 2 )
*     ..
*     .. Local Arrays ..
      CHARACTER          TRANSS( NTRAN )
      INTEGER            ISEED( 4 ), ISEEDY( 4 ), PLASMA_TRANSS( NTRAN )
      DOUBLE PRECISION   RESULT( NTESTS )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DGET06, ZLANGE
      EXTERNAL           DGET06, ZLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAERH, ALAHD, ALASUM, XLAENV, ZERRGE, ZGECON,
     $                   ZGERFS, ZGET02, ZGET04,
     $                   ZGETRF, ZGETRI, ZGETRS, ZLACPY, ZLARHS, ZLASET,
     $                   ZLATB4, ZLATMS
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DCMPLX, MAX, MIN
*     ..
*     .. Scalars in Common ..
      LOGICAL            LERR, OK
      CHARACTER*32       SRNAMT
      INTEGER            INFOT, NUNIT
*     ..
*     .. Common blocks ..
      COMMON             / INFOC / INFOT, NUNIT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Data statements ..
      DATA               ISEEDY / 1988, 1989, 1990, 1991 / ,
*     $                   TRANSS / 'N', 'T', 'C' /
     $                   TRANSS / 'N' /
     $                   PLASMA_TRANSS / PLASMANOTRANS /
*     ..
*     .. Executable Statements ..
*
*     Initialize constants and the random number seed.
*
      PATH( 1: 1 ) = 'Zomplex precision'
      PATH( 2: 3 ) = 'GE'
      RCONDO = ZERO
      RCONDI = ZERO
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE
*
*     Test the error exits
*
      CALL XLAENV( 1, 1 )
      IF( TSTERR )
     $   CALL ZERRGE( PATH, NOUT )
      INFOT = 0
      CALL XLAENV( 2, 2 )
*
*     Do for each value of M in MVAL
*
      DO 120 IM = 1, NM
         M = MVAL( IM )
         LDA = MAX( 1, M )
*
*        Do for each value of N in NVAL
*
         DO 110 IN = 1, NN
            N = NVAL( IN )
            XTYPE = 'N'
            NIMAT = NTYPES
            IF( M.LE.0 .OR. N.LE.0 )
     $         NIMAT = 1
*
            DO 100 IMAT = 1, NIMAT
*
*              Do the tests only if DOTYPE( IMAT ) is true.
*
               IF( .NOT.DOTYPE( IMAT ) )
     $            GO TO 100
*
*              Skip types 5, 6, or 7 if the matrix size is too small.
*
               ZEROT = IMAT.GE.5 .AND. IMAT.LE.7
               IF( ZEROT .AND. N.LT.IMAT-4 )
     $            GO TO 100
*
*              Set up parameters with ZLATB4 and generate a test matrix
*              with ZLATMS.
*
               CALL ZLATB4( PATH, IMAT, M, N, TYPE, KL, KU, ANORM, MODE,
     $                      CNDNUM, DIST )
*
               SRNAMT = 'ZLATMS'
               CALL ZLATMS( M, N, DIST, ISEED, TYPE, RWORK, MODE,
     $                      CNDNUM, ANORM, KL, KU, 'No packing', A, LDA,
     $                      WORK, INFO )
*
*              Check error code from ZLATMS.
*
               IF( INFO.NE.0 ) THEN
                  CALL ALAERH( PATH, 'ZLATMS', INFO, 0, ' ', M, N, -1,
     $                         -1, -1, IMAT, NFAIL, NERRS, NOUT )
                  GO TO 100
               END IF
*
*              For types 5-7, zero one or more columns of the matrix to
*              test that INFO is returned correctly.
*
               IF( ZEROT ) THEN
                  IF( IMAT.EQ.5 ) THEN
                     IZERO = 1
                  ELSE IF( IMAT.EQ.6 ) THEN
                     IZERO = MIN( M, N )
                  ELSE
                     IZERO = MIN( M, N ) / 2 + 1
                  END IF
                  IOFF = ( IZERO-1 )*LDA
                  IF( IMAT.LT.7 ) THEN
                     DO 20 I = 1, M
                        A( IOFF+I ) = ZERO
   20                CONTINUE
                  ELSE
                     CALL ZLASET( 'Full', M, N-IZERO+1, DCMPLX( ZERO ),
     $                            DCMPLX( ZERO ), A( IOFF+1 ), LDA )
                  END IF
               ELSE
                  IZERO = 0
               END IF
*
*              These lines, if used in place of the calls in the DO 60
*              loop, cause the code to bomb on a Sun SPARCstation.
*
*               ANORMO = ZLANGE( 'O', M, N, A, LDA, RWORK )
*               ANORMI = ZLANGE( 'I', M, N, A, LDA, RWORK )
*
*              Do for each blocksize in NBVAL
*
               DO 90 INB = 1, NNB
                  NB = NBVAL( INB )
                  CALL XLAENV( 1, NB )
                  IB = IBVAL( INB )
                  IF ( (MAX(M, N) / 25) .GT. NB ) THEN
                     GOTO 90
                  END IF
                  CALL PLASMA_SET( PLASMA_TILE_SIZE, NB, INFO )
                  CALL PLASMA_SET( PLASMA_INNER_BLOCK_SIZE, IB, INFO )
*
*                 ALLOCATE HL and HPIV
*
c$$$                  CALL PLASMA_ALLOC_WORKSPACE_ZGETRF_INCPIV(
c$$$     $                 M, N, HL, HPIV, INFO ) 
*
*                 Compute the LU factorization of the matrix.
*
                  CALL ZLACPY( 'Full', M, N, A, LDA, AFAC, LDA )
                  SRNAMT = 'ZGETRF'
c$$$                  CALL PLASMA_ZGETRF_INCPIV( M, N, AFAC, LDA, HL, HPIV, 
c$$$     $                 INFO )
                  CALL PLASMA_ZGETRF( M, N, AFAC, LDA, IWORK, 
     $                 INFO )
*
*                 Check error code from ZGETRF.
*
                  IF( INFO.NE.IZERO )
     $               CALL ALAERH( PATH, 'ZGETRF', INFO, IZERO, ' ', M,
     $                            N, -1, -1, NB, IMAT, NFAIL, NERRS,
     $                            NOUT )
                  TRFCON = .FALSE.
                  NT = 0
*
                  IF( M.NE.N .OR. INFO.GT.0 ) THEN
*
*                    Do only the condition estimate if INFO > 0.
*
                     TRFCON = .TRUE.
                     ANORMO = ZLANGE( 'O', M, N, A, LDA, RWORK )
                     ANORMI = ZLANGE( 'I', M, N, A, LDA, RWORK )
                     RCONDO = ZERO
                     RCONDI = ZERO
                  END IF
*
*                 Print information about the tests so far that did not
*                 pass the threshold.
*
                  DO 30 K = 1, NT
                     IF( RESULT( K ).GE.THRESH ) THEN
                        IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $                     CALL ALAHD( NOUT, PATH )
                        WRITE( NOUT, FMT = 9999 )M, N, NB, IMAT, K,
     $                     RESULT( K )
                        NFAIL = NFAIL + 1
                     END IF
   30             CONTINUE
                  NRUN = NRUN + NT
*
*                 Skip the remaining tests if this is not the first
*                 block size or if M .ne. N.  Skip the solve tests if
*                 the matrix is singular.
*
*                  IF( INB.GT.1 .OR. M.NE.N )
*     $               GO TO 90
                  IF( TRFCON )
     $               GO TO 70
*
                  DO 60 IRHS = 1, NNS
                     NRHS = NSVAL( IRHS )
                     XTYPE = 'N'
*
                     DO 50 ITRAN = 1, NTRAN
                        TRANS = TRANSS( ITRAN )
                        PLASMA_TRANS = PLASMA_TRANSS( ITRAN )
                        IF( ITRAN.EQ.1 ) THEN
                           RCONDC = RCONDO
                        ELSE
                           RCONDC = RCONDI
                        END IF
*
*+    TEST 3
*                       Solve and compute residual for A * X = B.
*
                        SRNAMT = 'ZLARHS'
                        CALL ZLARHS( PATH, XTYPE, ' ', TRANS, N, N, KL,
     $                               KU, NRHS, A, LDA, XACT, LDA, B,
     $                               LDA, ISEED, INFO )
                        XTYPE = 'C'
*
                        CALL ZLACPY( 'Full', N, NRHS, B, LDA, X, LDA )
                        SRNAMT = 'ZGETRS'
c$$$                        CALL PLASMA_ZGETRS_INCPIV( PLASMA_TRANS, N, 
c$$$     $                       NRHS, AFAC, LDA, HL, HPIV,
c$$$     $                       X, LDA, INFO )
                        CALL PLASMA_ZGETRS( PLASMA_TRANS, N, 
     $                       NRHS, AFAC, LDA, IWORK,
     $                       X, LDA, INFO )
*
*                       Check error code from ZGETRS.
*
                        IF( INFO.NE.0 )
     $                     CALL ALAERH( PATH, 'ZGETRS', INFO, 0, TRANS,
     $                                  N, N, -1, -1, NRHS, IMAT, NFAIL,
     $                                  NERRS, NOUT )
*
                        CALL ZLACPY( 'Full', N, NRHS, B, LDA, WORK,
     $                               LDA )
                        CALL ZGET02( TRANS, N, N, NRHS, A, LDA, X, LDA,
     $                               WORK, LDA, RWORK, RESULT( 3 ) )
*
*+    TEST 4
*                       Check solution from generated exact solution.
*
                        CALL ZGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC,
     $                               RESULT( 4 ) )
*
*                       Print information about the tests that did not
*                       pass the threshold.
*
                        DO 40 K = 3, 4
                           IF( RESULT( K ).GE.THRESH ) THEN
                              IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $                           CALL ALAHD( NOUT, PATH )
                              WRITE( NOUT, FMT = 9998 )TRANS, N, NRHS,
     $                           IMAT, K, RESULT( K )
                              NFAIL = NFAIL + 1
                           END IF
   40                   CONTINUE
                        NRUN = NRUN + 2
   50                CONTINUE
   60             CONTINUE
*
   70             CONTINUE
*
*                 DEALLOCATE HL and HPIV
*
c$$$                  CALL PLASMA_DEALLOC_HANDLE( HL, INFO )
c$$$                  CALL PLASMA_DEALLOC_HANDLE( HPIV, INFO )
   90          CONTINUE
  100       CONTINUE
*
  110    CONTINUE
  120 CONTINUE
*
*     Print a summary of the results.
*
      CALL ALASUM( PATH, NOUT, NFAIL, NRUN, NERRS )
*
 9999 FORMAT( ' M = ', I5, ', N =', I5, ', NB =', I4, ', type ', I2,
     $      ', test(', I2, ') =', G12.5 )
 9998 FORMAT( ' TRANS=''', A1, ''', N =', I5, ', NRHS=', I3, ', type ',
     $      I2, ', test(', I2, ') =', G12.5 )
 9997 FORMAT( ' NORM =''', A1, ''', N =', I5, ',', 10X, ' type ', I2,
     $      ', test(', I2, ') =', G12.5 )
      RETURN
*
*     End of ZCHKGE
*
      END
