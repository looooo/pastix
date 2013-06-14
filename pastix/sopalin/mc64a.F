CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC
CC COPYRIGHT (c) 1998-2002  I.S. Duff and J. Koster
CC
CC PACKAGE MC64A/AD/AC/AZ - v1.2b 
CC Generated on Tue Nov 26 03:20:26 MET 2002 at Parallab
CC
CC AUTHORS Iain S. Duff   (1)   i.duff@rl.ac.uk
CC         Jacko Koster   (2)   jak@ii.uib.no
CC
CC (1) Rutherford Appleton Laboratory, Chilton OX11 0QX, ENGLAND.
CC     Also CERFACS, 42 Avenue G. Coriolis, 31057 Toulouse, FRANCE.
CC
CC (2) Parallab/UNIFOB, Dept. Informatics, University of Bergen, 
CC     N-5020 Bergen, NORWAY. 
CC
CC Problems, queries, and requests for up-to-date versions of the
CC MC64 package should be sent to Jacko Koster.
CC
CC THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
CC EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.
CC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC
CC *** Conditions on external use ***
CC
CC The user shall acknowledge the contribution of this package in any 
CC publication of material dependent upon the use of the package. The 
CC user shall use reasonable endeavours to notify the authors of the 
CC package of this publication.
CC
CC The user can modify this code but, at no time shall the right or 
CC title to all or any part of this package pass to the user. The user 
CC shall make available free of charge to the authors for any purpose 
CC all information relating to any alteration or addition made to this 
CC package for the purposes of extending the capabilities or enhancing 
CC the performance of this package.
CC
CC This package and modified versions of this package may not be sold.
CC Users wanting to licence their own copy of these routines must first 
CC contact the authors.
CC
CC One may make copies of the package or modify it for their use
CC provided that the copies, modified or otherwise, are not sold, and 
CC are used under the same terms and conditions as listed here.
CC
CC None of the comments from the Copyright notice up to and including 
CC this one shall be removed or altered in any way.
CC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCC LAST UPDATE Tue Nov 26 03:20:26 MET 2002
      SUBROUTINE MC64ID(ICNTL,CNTL)
      IMPLICIT NONE
C
C *** Copyright (c) 2002  I.S. Duff and J. Koster                   ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C
C Purpose
C =======
C
C The components of the arrays ICNTL/CNTL control the action of 
C MC64A. Default values for these are set in this subroutine.
C See MC64A for a description of the components of ICNTL/CNTL.
C
C Parameters
C ==========
C
      INTEGER NICNTL, NCNTL
      PARAMETER (NICNTL=10, NCNTL=10)

      INTEGER ICNTL(NICNTL)
      DOUBLE PRECISION CNTL(NCNTL)
C
C ICNTL(1) is set to default value 6.
C ICNTL(2) is set to default value 6.
C ICNTL(3) is set to default value -1.
C ICNTL(4) is set to default value -1.
C ICNTL(5) is set to default value 0.
C
C ICNTL(6) to ICNTL(NICNTL) are not used by MC64A but are set to 
C    zero in this routine.
C
C CNTL(1) is set to default value 0.0.
C CNTL(2) is set to default value 0.0. 
C CNTL(3) to CNTL(NCNTL) are not used by MC64A but are set to zero
C    in this routine.
C
C Local variables
      INTEGER I
C   
C Initialization of the ICNTL array.
      ICNTL(1) =  6
      ICNTL(2) =  6
      ICNTL(3) = -1
      ICNTL(4) = -1
      ICNTL(5) =  0
      DO 10 I = 6,NICNTL
        ICNTL(I) = 0
   10 CONTINUE
C Initialization of the CNTL array.
      CNTL(1) = 0.0D0
      CNTL(2) = 0.0D0
      DO 20 I = 3,NCNTL
        CNTL(I) = 0.0D0
   20 CONTINUE

      RETURN
      END

C**********************************************************************
CCCCC LAST UPDATE Tue Nov 26 03:20:26 MET 2002
      SUBROUTINE MC64AD(JOB,M,N,NE,IP,IRN,A,NUM,PERM,LIW,IW,LDW,DW,
     &           ICNTL,CNTL,INFO)
      IMPLICIT NONE
C
C *** Copyright (c) 2002  I.S. Duff and J. Koster                   ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C
C Purpose
C =======
C
C This subroutine attempts to find a permutation for an MxN, M>=N, 
C sparse matrix A = {a_ij} that makes the permuted matrix have N 
C entries on its diagonal.
C If the matrix is structurally nonsingular, the subroutine optionally
C returns a permutation that maximizes the smallest element on the 
C diagonal, maximizes the sum of the diagonal entries, or maximizes 
C the product of the diagonal entries of the permuted matrix. 
C For the latter option, the subroutine also finds scaling factors that 
C may be used to scale the matrix so that the nonzero diagonal entries 
C of the permuted matrix are one in absolute value and all the 
C off-diagonal entries are less than or equal to one in absolute value.
C The natural logarithms of the scaling factors u(i), i=1..M, for the 
C rows and v(j), j=1..N, for the columns are returned so that the 
C scaled matrix B = {b_ij} has entries b_ij = a_ij * EXP(u_i + v_j).
C The scaling factors are returned by this subroutine, but the actual
C scaling of the matrix has to be performed by the calling program.
C
C Parameters
C ==========
C
      INTEGER NICNTL, NCNTL, NINFO
      PARAMETER (NICNTL=10, NCNTL=10, NINFO=10)

      INTEGER JOB,M,N,NE,NUM,LIW,LDW
      INTEGER IP(N+1),IRN(NE),PERM(M),IW(LIW)
      INTEGER ICNTL(NICNTL),INFO(NINFO)
      DOUBLE PRECISION A(NE)
      DOUBLE PRECISION DW(LDW),CNTL(NCNTL)
C
C JOB is an INTEGER variable which must be set by the user to control
C the action. It is not altered by the subroutine.
C Possible values for JOB are:
C   1 Compute a column permutation of the matrix so that the
C     permuted matrix has as many entries on its diagonal as possible. 
C     The values on the diagonal are of arbitrary size. HSL subroutine 
C     MC21A/MC64Z is used for this. See [1].
C   2 Compute a column permutation of the matrix so that the smallest 
C     value on the diagonal of the permuted matrix is maximized.
C     See [3].
C   3 Compute a column permutation of the matrix so that the smallest
C     value on the diagonal of the permuted matrix is maximized.
C     The algorithm differs from the one used for JOB = 2 and may
C     have quite a different performance. See [2].
C   4 Compute a column permutation of the matrix so that the sum
C     of the diagonal entries of the permuted matrix is maximized.
C     See [3].
C   5 Compute a column permutation of the matrix so that the product
C     of the diagonal entries of the permuted matrix is maximized
C     and vectors to scale the matrix so that the nonzero diagonal 
C     entries of the permuted matrix are one in absolute value and 
C     all the off-diagonal entries are less than or equal to one in 
C     absolute value. See [3].
C   6 Compute a column permutation of the matrix so that the product
C     of the diagonal entries of the permuted matrix is maximized
C     and vectors to scale the matrix so that the nonzero diagonal
C     entries of the permuted matrix are one in absolute value and
C     all the off-diagonal entries are less than or equal to one in
C     absolute value. 
C     The algorithm differs from the one used for JOB = 5 and may
C     have quite a different performance.
C  Restriction: 1 <= JOB <= 6.
C 
C M is an INTEGER variable which must be set by the user to the
C   number of rows of the matrix A. It is not altered by the 
C   subroutine. Restriction: M >= N.
C
C N is an INTEGER variable which must be set by the user to the
C   number of columns of the matrix A. It is not altered by the 
C   subroutine. Restriction: N >= 1.
C
C NE is an INTEGER variable which must be set by the user to the
C   number of entries in the matrix. It is not altered by the 
C   subroutine. Restriction: NE >= 1.
C
C IP is an INTEGER array of length N+1.
C   IP(J), J=1..N, must be set by the user to the position in array IRN 
C   of the first row index of an entry in column J. IP(N+1) must be set
C   to NE+1. It is not altered by the subroutine.
C
C IRN is an INTEGER array of length NE. 
C   IRN(K), K=1..NE, must be set by the user to hold the row indices of
C   the entries of the matrix. Those belonging to column J must be 
C   stored contiguously in the positions IP(J)..IP(J+1)-1. The ordering
C   of the row indices within each column is unimportant. Repeated 
C   entries are not allowed. The array IRN is not altered by the 
C   subroutine.
C
C A is a DOUBLE PRECISION array of length NE. 
C   The user must set A(K), K=1..NE, to the numerical value of the 
C   entry that corresponds to IRN(K). 
C   It is not used by the subroutine when JOB = 1. 
C   It is not altered by the subroutine.
C
C NUM is an INTEGER variable that need not be set by the user.
C   On successful exit, NUM will be the number of entries on the 
C   diagonal of the permuted matrix.
C   If NUM < N, the matrix is structurally singular.
C
C PERM is an INTEGER array of length M that need not be set by the 
C   user. On successful exit, PERM can be interpreted in any of the
C   following ways: 
C
C   1. If M=N, PERM contains the column permutation.
C      Column PERM(I) of the original matrix is column I in the 
C      permuted matrix, I=1..N. 
C      (This was the definition of parameter CPERM in versions of 
C      MC64A before version 1.2b)
C
C   2. If M>=N, PERM contains the row permutation.
C      Row I of the original matrix is row PERM(I) in the 
C      permuted matrix, I=1..M.
C
C LIW is an INTEGER variable that must be set by the user to
C   the dimension of array IW. It is not altered by the subroutine.
C   Restriction:
C     JOB = 1 :  LIW >=  4N +  M
C     JOB = 2 :  LIW >=  2N + 2M
C     JOB = 3 :  LIW >=  8N + 2M + NE
C     JOB = 4 :  LIW >=  3N + 2M
C     JOB = 5 :  LIW >=  3N + 2M
C     JOB = 6 :  LIW >=  3N + 2M + NE
C 
C IW is an INTEGER array of length LIW that is used for workspace.
C 
C LDW is an INTEGER variable that must be set by the user to the
C   dimension of array DW. It is not altered by the subroutine.
C   Restriction:
C     JOB = 1 :  LDW not used
C     JOB = 2 :  LDW >=      M
C     JOB = 3 :  LDW >=          NE
C     JOB = 4 :  LDW >=     2M + NE
C     JOB = 5 :  LDW >= N + 2M + NE
C     JOB = 6 :  LDW >= N + 3M + NE
C 
C DW is a REAL array of length LDW used for workspace. 
C   If JOB = 5,6, on return, DW(i) contains u_i, i=1..M, and 
C   DW(M+j) contains v_j, j=1..N.
C 
C ICNTL is an INTEGER array of length NICNTL.
C   Its components control the output of MC64A and must be set by the 
C   user before calling MC64A. They are not altered by the subroutine.
C
C   ICNTL(1) must be set to specify the output stream for
C   error messages. If ICNTL(1) < 0, messages are suppressed.
C   The default value set by MC64I is 6 (standard output).
C
C   ICNTL(2) must be set by the user to specify the output stream for
C   warning messages. If ICNTL(2) < 0, messages are suppressed.
C   The default value set by MC64I is 6 (standard output).
C
C   ICNTL(3) must be set by the user to specify the output stream for
C   diagnostic messages. If ICNTL(3) < 0, messages are suppressed.
C   The default value set by MC64I is -1 (no diagnostics).
C
C   ICNTL(4) must be set by the user to control the printing of
C   diagnostic messages.
C   If ICNTL(4) = 0, terse information on the input and output
C   parameters is printed (only the first ten entries of arrays).
C   If ICNTL(4) = 1, all information on the input and output
C   parameters is printed. (This can create large output if the input
C   data is large.)
C   If ICNTL(4) is set to any other value, messages are suppressed.
C   The default value set by MC64I is -1 (no diagnostics).
C
C   ICNTL(5) must be set by the user to a value other than 0 to avoid
C   checking of the input data.
C   The default value set by MC64I is 0 (checking of data).
C   If left at the default value, the incoming data is checked for
C   out-of-range indices and duplicates.  Setting ICNTL(5) to any
C   other will avoid the checks but is likely to cause problems
C   later if out-of-range indices or duplicates are present.
C   The user should set ICNTL(5) non-zero, if the data is known not
C   to contain such problems. The code will exhibit undefined
C   behaviour in case data checking is not done and the
C   input data does not satisfy the restrictions as listed
C   elsewhere.
C    
C CNTL is a DOUBLE PRECISION array of length NCNTL. 
C   Its components control the output of MC64A and must be set by the 
C   user before calling MC64A. They are not altered by the subroutine.
C
C   CNTL(1) must be set to specify the relaxation parameter.
C   It is used by MC64 only if JOB = 3,4,5,6.
C   It must be set to a non-negative value (usually close to zero).
C   If CNTL(1) < 0.0, it is treated as 0.0.
C   The default value set by MC64I is 0.0.
C
C   CNTL(1) is a relaxation parameter. A positive value will lead to 
C   matchings computed by MC64A that are not optimal/maximal in some 
C   sense but only nearly so. However, these non-optimal matchings are 
C   often computed more quickly. Appropriate values for CNTL(1) are 
C   problem dependent but usually slightly larger than 0.0.
C
C   CNTL(2) must be set to the largest positive real value that 
C   can be stored. CNTL(2) can be initialized with the 
C   LAPACK function call SLAMCH('O') for the single precision versions
C   and with DLAMCH('O') for the double precision versions.
C   If a Fortran 90/95 compiler is used, CNTL(2) can also be 
C   initialized with the Fortran 90/95 intrinsic function HUGE().
C
C INFO is an INTEGER array of length NINFO which need not be set by the 
C   user. INFO(1) is set non-negative to indicate success. A negative 
C   value is returned if an error occurred, a positive value if a 
C   warning occurred. INFO(2) holds further information on the error.
C   On exit from the subroutine, INFO(1) will take one of the 
C   following values:
C    0 : successful entry (for structurally nonsingular matrix).
C   +1 : successful entry (for structurally singular matrix).
C   +2 : the returned scaling factors are large and may cause
C        overflow when used to scale the matrix. 
C        (For JOB = 4,5,6 entries only.)
C   +4 : CNTL(1) is negative and treated as zero.
C   -1 : JOB < 1 or JOB > 6.  Value of JOB held in INFO(2).
C   -2 : M < 1 or N < 1.  Value of invalid M or N held in INFO(2).
C   -3 : NE < 1.  Value of NE held in INFO(2).
C   -4 : the defined length LIW violates the restriction on LIW.
C        Value of LIW required given by INFO(2).
C   -5 : the defined length LDW violates the restriction on LDW.
C        Value of LDW required given by INFO(2).
C   -6 : entries are found whose row indices are out of range. INFO(2)
C        contains the index of a column in which such an entry is found.
C        (This value can be returned only if ICNTL(5) was set to zero.)
C   -7 : repeated entries are found. INFO(2) contains the index of a 
C        column in which such entries are found.
C        (This value can be returned only if ICNTL(5) was set to zero.)
C   -8 : CNTL(2) <= 0.0.
C
C   A return with one of the values INFO(1)=+3,+5,+6,+7 is also possible. 
C   These values are combinations of the above warnings (+1,+2,+4) and 
C   correspond to the sum of the constituent warnings.
C
C   INFO(3) to INFO(NINFO) are not currently used and are set to zero
C        by the routine.
C
C References:
C  [1] I. S. Duff, (1981),
C      "Algorithm 575. Permutations for a zero-free diagonal",
C      ACM Trans. Math. Software 7(3), 387-390.
C  [2] I. S. Duff and J. Koster, (1998),
C      "The design and use of algorithms for permuting large
C      entries to the diagonal of sparse matrices",
C      SIAM J. Matrix Anal. Appl., vol. 20, no. 4, pp. 889-901.
C  [3] I. S. Duff and J. Koster, (2001),
C      "On algorithms for permuting large entries to the diagonal 
C      of sparse matrices",
C      SIAM J. Matrix Anal. Appl., vol. 22, no. 4, pp. 973-996.

C Local variables and parameters
      INTEGER I,J,K,wARN1,WARN2,WARN4
      DOUBLE PRECISION FACT,ZERO,ONE,RINF
      PARAMETER (ZERO=0.0D+00,ONE=1.0D+0)
C External routines and functions
      EXTERNAL MC64ZD,MC64BD,MC64RD,MC64SD,MC64WD
C Intrinsic functions
      INTRINSIC ABS,LOG

C Set RINF to largest positive real number (infinity)
C# if defined(MC64F90)
C     INTRINSIC HUGE
C     RINF = HUGE(RINF)
C# else
       RINF = CNTL(2)
C# endif
      WARN1 = 0
      WARN2 = 0
      WARN4 = 0
C Check value of JOB
      IF (JOB.LT.1 .OR. JOB.GT.6) THEN
        INFO(1) = -1
        INFO(2) = JOB
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'JOB',JOB
        GO TO 99
      ENDIF
C Check value of M 
      IF (M.LT.1 .OR. M.LT.N) THEN
        INFO(1) = -2
        INFO(2) = M
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'M',M
        GO TO 99
      ENDIF
C Check value of N
      IF (N.LT.1) THEN
        INFO(1) = -2
        INFO(2) = N
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'N',N
        GO TO 99
      ENDIF
C Check value of NE
      IF (NE.LT.1) THEN
        INFO(1) = -3
        INFO(2) = NE
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'NE',NE
        GO TO 99
      ENDIF
C Check LIW
      IF (JOB.EQ.1) K = 4*N +   M
      IF (JOB.EQ.2) K = 2*N + 2*M
      IF (JOB.EQ.3) K = 8*N + 2*M + NE
      IF (JOB.EQ.4) K = 3*N + 2*M
      IF (JOB.EQ.5) K = 3*N + 2*M
      IF (JOB.EQ.6) K = 3*N + 2*M + NE
      IF (LIW.LT.K) THEN
        INFO(1) = -4
        INFO(2) = K
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9004) INFO(1),K
        GO TO 99
      ENDIF
C Check LDW; If JOB = 1, do not check
      IF (JOB.GT.1) THEN
        IF (JOB.EQ.2) K =       M
        IF (JOB.EQ.3) K =           NE
        IF (JOB.EQ.4) K =     2*M + NE
        IF (JOB.EQ.5) K = N + 2*M + NE
        IF (JOB.EQ.6) K = N + 3*M + NE
        IF (LDW.LT.K) THEN
          INFO(1) = -5
          INFO(2) = K
          IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9005) INFO(1),K
          GO TO 99
        ENDIF
      ENDIF
      IF (ICNTL(5).EQ.0) THEN
C Check row indices. Use IW(1:M) as workspace
        DO 3 I = 1,M
          IW(I) = 0
    3   CONTINUE
        DO 6 J = 1,N
          DO 4 K = IP(J),IP(J+1)-1
            I = IRN(K)
C Check for row indices that are out of range
            IF (I.LT.1 .OR. I.GT.M) THEN
              INFO(1) = -6
              INFO(2) = J
              IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9006) INFO(1),J,I
              GO TO 99
            ENDIF
C Check for repeated row indices within a column
            IF (IW(I).EQ.J) THEN
              INFO(1) = -7
              INFO(2) = J
              IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9007) INFO(1),J,I 
              GO TO 99
            ELSE
              IW(I) = J
            ENDIF
    4     CONTINUE
    6   CONTINUE
      ENDIF
C Check that CNTL(2) was set by the user to a positive value
      IF (CNTL(2).LE.ZERO) THEN
        INFO(1) = -8
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9008) INFO(1) 
        GO TO 99
      ENDIF

C Print diagnostics on input
      IF (ICNTL(3).GE.0) THEN
        IF (ICNTL(4).EQ.0 .OR. ICNTL(4).EQ.1) THEN
          WRITE(ICNTL(3),9020) JOB,M,N,NE
          IF (ICNTL(4).EQ.0) THEN
            WRITE(ICNTL(3),9021) (IP(J),J=1,MIN(10,N+1))
            WRITE(ICNTL(3),9022) (IRN(J),J=1,MIN(10,NE))
            IF (JOB.GT.1) WRITE(ICNTL(3),9023) (A(J),J=1,MIN(10,NE))
          ELSEIF (ICNTL(4).EQ.1) THEN
            WRITE(ICNTL(3),9021) (IP(J),J=1,N+1)
            WRITE(ICNTL(3),9022) (IRN(J),J=1,NE)
            IF (JOB.GT.1) WRITE(ICNTL(3),9023) (A(J),J=1,NE)
          ENDIF
          WRITE(ICNTL(3),9024) (ICNTL(J),J=1,NICNTL)
          WRITE(ICNTL(3),9025) (CNTL(J),J=1,NCNTL)
        ENDIF
      ENDIF

C Set components of INFO to zero
      DO 8 I=1,NINFO
        INFO(I) = 0
    8 CONTINUE

C Compute maximum matching
      IF (JOB.EQ.1) THEN
C Put length of column J in IW(J)
        DO 10 J = 1,N
          IW(J) = IP(J+1) - IP(J)
   10   CONTINUE
C IW(N+1:3N+M+N) is workspace 
        CALL MC64ZD(M,N,IRN,NE,IP,IW(1),PERM,NUM,
     &       IW(N+1),IW(2*N+1),IW(3*N+1),IW(3*N+M+1))
        GO TO 90
      ENDIF

C Compute bottleneck matching
      IF (JOB.EQ.2) THEN
C Pass CNTL(1) to MC64B through DW(1)
        DW(1) = MAX(ZERO,CNTL(1))
C IW(1:2N+2M), DW(1:M) are workspaces 
        CALL MC64BD(M,N,NE,IP,IRN,A,PERM,NUM,
     &       IW(1),IW(N+1),IW(2*N+1),IW(2*N+M+1),DW,RINF)
        GO TO 90
      ENDIF

C Compute bottleneck matching
      IF (JOB.EQ.3) THEN
C Copy IRN(K) into IW(K), ABS(A(K)) into DW(K), K=1..NE
        DO 20 K = 1,NE
          IW(K) = IRN(K)
          DW(K) = ABS(A(K))
   20   CONTINUE
C Sort entries in each column by decreasing value. 
        CALL MC64RD(N,NE,IP,IW,DW)
C Pass CNTL(1) to MC64S through FACT
        FACT = MAX(ZERO,CNTL(1))
C IW(NE+1:NE+5N+M+(3N+M)) is workspace 
        CALL MC64SD(M,N,NE,IP,IW(1),DW,PERM,NUM,IW(NE+1),
     &       IW(NE+N+1),IW(NE+2*N+1),IW(NE+3*N+1),IW(NE+4*N+1),
     &       IW(NE+5*N+1),IW(NE+5*N+M+1),FACT,RINF)
        GO TO 90
      ENDIF

      IF (JOB.EQ.4) THEN
        DO 50 J = 1,N
          FACT = ZERO
          DO 30 K = IP(J),IP(J+1)-1
            IF (ABS(A(K)).GT.FACT) FACT = ABS(A(K))
   30     CONTINUE
          DO 40 K = IP(J),IP(J+1)-1
            DW(2*M+K) = FACT - ABS(A(K))
   40     CONTINUE
   50   CONTINUE
C B = DW(2M+1:2M+NE); IW(1:3N+2M) and DW(1:2M) are workspaces
C Pass CNTL(1) to MC64W through DW(1)
C Pass JOB to MC64W through IW(1)
        DW(1) = MAX(ZERO,CNTL(1))
        IW(1) = JOB
C Call MC64W
        CALL MC64WD(M,N,NE,IP,IRN,DW(2*M+1),PERM,NUM,
     &       IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),IW(3*N+M+1),
     &       DW(1),DW(M+1),RINF)
        GO TO 90
      ENDIF

      IF (JOB.EQ.5 .or. JOB.EQ.6) THEN
        IF (JOB.EQ.5) THEN
          DO 75 J = 1,N
            FACT = ZERO
            DO 60 K = IP(J),IP(J+1)-1
              DW(2*M+N+K) = ABS(A(K))
              IF (DW(2*M+N+K).GT.FACT) FACT = DW(2*M+N+K)
   60       CONTINUE
            DW(2*M+J) = FACT
            IF (FACT.NE.ZERO) THEN
              FACT = LOG(FACT)
              DO 70 K = IP(J),IP(J+1)-1
                IF (DW(2*M+N+K).NE.ZERO) THEN
                  DW(2*M+N+K) = FACT - LOG(DW(2*M+N+K))
                ELSE
                  DW(2*M+N+K) = RINF/N
                ENDIF
   70         CONTINUE
            ELSE
              DO 71 K = IP(J),IP(J+1)-1
                DW(2*M+N+K) = ONE
   71         CONTINUE
            ENDIF
   75     CONTINUE
        ENDIF

        IF (JOB.EQ.6) THEN
          DO 175 K = 1,NE
            IW(3*N+2*M+K) = IRN(K)
            DW(2*M+N+K) = ABS(A(K))
  175     CONTINUE
          DO 61 I = 1,M
            DW(2*M+N+NE+I) = ZERO
   61     CONTINUE
          DO 63 J = 1,N
            DO 62 K = IP(J),IP(J+1)-1
              I = IRN(K)
              IF (DW(2*M+N+K).GT.DW(2*M+N+NE+I)) THEN
                DW(2*M+N+NE+I) = DW(2*M+N+K)
              ENDIF
   62       CONTINUE
   63     CONTINUE
          DO 64 I = 1,M
            IF (DW(2*M+N+NE+I).NE.ZERO) THEN
              DW(2*M+N+NE+I) = 1/DW(2*M+N+NE+I)
            ENDIF
   64     CONTINUE
          DO 66 J = 1,N
            DO 65 K = IP(J),IP(J+1)-1
              I = IRN(K)
              DW(2*M+N+K) = DW(2*M+N+NE+I) * DW(2*M+N+K)
   65       CONTINUE
   66     CONTINUE
          CALL MC64RD(N,NE,IP,IW(3*N+2*M+1),DW(2*M+N+1))
          DO 176 J = 1,N
            IF (IP(J).NE.IP(J+1)) THEN
              FACT = DW(2*M+N+IP(J))
            ELSE
              FACT = ZERO
            ENDIF
            DW(2*M+J) = FACT
            IF (FACT.NE.ZERO) THEN
              FACT = LOG(FACT)
              DO 170 K = IP(J),IP(J+1)-1
                IF (DW(2*M+N+K).NE.ZERO) THEN
                  DW(2*M+N+K) = FACT - LOG(DW(2*M+N+K))
                ELSE
                  DW(2*M+N+K) = RINF/N
                ENDIF
  170         CONTINUE
            ELSE
              DO 171 K = IP(J),IP(J+1)-1
                DW(2*M+N+K) = ONE
  171         CONTINUE
            ENDIF
  176     CONTINUE
        ENDIF
C Pass CNTL(1) to MC64W through DW(1)
C Pass JOB to MC64W through IW(1)
        DW(1) = MAX(ZERO,CNTL(1))
        IW(1) = JOB
C Call MC64W
        IF (JOB.EQ.5) THEN
          CALL MC64WD(M,N,NE,IP,IRN,DW(2*M+N+1),PERM,NUM,
     &         IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),IW(3*N+M+1),
     &         DW(1),DW(M+1),RINF)
        ENDIF

        IF (JOB.EQ.6) THEN
          CALL MC64WD(M,N,NE,IP,IW(3*N+2*M+1),DW(2*M+N+1),PERM,NUM,
     &         IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),IW(3*N+M+1),
     &         DW(1),DW(M+1),RINF)
        ENDIF
        IF (JOB.EQ.6) THEN
          DO 79 I = 1,M
            IF (DW(2*M+N+NE+I).NE.0) THEN
              DW(I) = DW(I) + LOG(DW(2*M+N+NE+I))
            ENDIF
   79     CONTINUE
        ENDIF
        IF (NUM.EQ.N) THEN
          DO 80 J = 1,N
            IF (DW(2*M+J).NE.ZERO) THEN
              DW(M+J) = DW(M+J) - LOG(DW(2*M+J))
            ELSE
              DW(M+J) = ZERO
            ENDIF
   80     CONTINUE
        ENDIF
C Check size of row and column scaling factors
        FACT = 0.5*LOG(RINF)
        DO 86 I = 1,M
          IF (DW(I).LT.FACT) GO TO 86
          WARN2 = 2
C Scaling factor is large, return with warning
          GO TO 90
   86   CONTINUE 
        DO 87 J = 1,N
          IF (DW(M+J).LT.FACT) GO TO 87
          WARN2 = 2
C Scaling factor is large, return with warning
          GO TO 90
   87   CONTINUE 
C       GO TO 90
      ENDIF

C If matrix is structurally singular, return with warning
   90 IF (NUM.LT.N) WARN1 = 1

C If CNTL(1) is negative and treated as zero, return with warning
      IF (JOB.EQ.4 .OR. JOB.EQ.5 .OR. JOB.EQ.6) THEN 
        IF (CNTL(1).LT.ZERO) WARN4 = 4
      ENDIF

C Set warning flag and print warnings (only if no errors were found)
      IF (INFO(1).EQ.0) THEN
        INFO(1) = WARN1 + WARN2 + WARN4
        IF (INFO(1).GT.0 .AND. ICNTL(2).GE.0) THEN
          WRITE(ICNTL(2),9010) INFO(1)
          IF (WARN1.EQ.1) WRITE(ICNTL(2),9011)
          IF (WARN2.EQ.2) WRITE(ICNTL(2),9012)
          IF (WARN4.EQ.4) WRITE(ICNTL(2),9014)
        ENDIF
      ENDIF
C Print diagnostics on output
      IF (ICNTL(3).GE.0) THEN
        IF (ICNTL(4).EQ.0 .OR. ICNTL(4).EQ.1) THEN
          WRITE(ICNTL(3),9030) (INFO(J),J=1,2)
          WRITE(ICNTL(3),9031) NUM
          IF (ICNTL(4).EQ.0) THEN
            WRITE(ICNTL(3),9032) (PERM(J),J=1,MIN(10,M))
            IF (JOB.EQ.5 .OR. JOB.EQ.6) THEN
              WRITE(ICNTL(3),9033) (DW(J),J=1,MIN(10,M))
              WRITE(ICNTL(3),9034) (DW(M+J),J=1,MIN(10,N))
            ENDIF
          ELSEIF (ICNTL(4).EQ.1) THEN
            WRITE(ICNTL(3),9032) (PERM(J),J=1,M)
            IF (JOB.EQ.5 .OR. JOB.EQ.6) THEN
              WRITE(ICNTL(3),9033) (DW(J),J=1,M)
              WRITE(ICNTL(3),9034) (DW(M+J),J=1,N)
            ENDIF
          ENDIF
        ENDIF
      ENDIF

C Return from subroutine.
   99 RETURN

 9001 FORMAT (' ****** Error in MC64AD. INFO(1) = ',I2,
     &        ' because ',(A),' = ',I10)
 9004 FORMAT (' ****** Error in MC64AD. INFO(1) = ',I2/
     &        '        LIW too small, must be at least ',I8)
 9005 FORMAT (' ****** Error in MC64AD. INFO(1) = ',I2/
     &        '        LDW too small, must be at least ',I8)
 9006 FORMAT (' ****** Error in MC64AD. INFO(1) = ',I2/
     &        '        Column ',I8,
     &        ' contains an entry with invalid row index ',I8)
 9007 FORMAT (' ****** Error in MC64AD. INFO(1) = ',I2/
     &        '        Column ',I8,
     &        ' contains two or more entries with row index ',I8)
 9008 FORMAT (' ****** Error in MC64AD. INFO(1) = ',I2/
     &        '        CNTL(2) is negative or zero.')

 9010 FORMAT (' ****** Warning from MC64AD. INFO(1) = ',I2)
 9011 FORMAT ('        - The matrix is structurally singular.')
 9012 FORMAT ('        - Some scaling factors may be too large.')
 9014 FORMAT ('        - CNTL(1) is negative and was treated as zero.')

 9020 FORMAT (' ****** Input parameters for MC64AD:'/
     &        ' JOB =',I10/' M   =',I10/' N   =',I10/' NE  =',I10)
 9021 FORMAT (' IP(1:N+1)   = ',8I8/(15X,8I8))
 9022 FORMAT (' IRN(1:NE)   = ',8I8/(15X,8I8))
 9023 FORMAT (' A(1:NE)     = ',4(1PD14.4)/(15X,4(1PD14.4)))
 9024 FORMAT (' ICNTL(1:10) = ',8I8/(15X,2I8))
 9025 FORMAT (' CNTL(1:10)  = ',4(1PD14.4)/(15X,4(1PD14.4)))
 9030 FORMAT (' ****** Output parameters for MC64AD:'/
     &        ' INFO(1:2)   = ',2I8)
 9031 FORMAT (' NUM         = ',I8)
 9032 FORMAT (' PERM(1:M)   = ',8I8/(15X,8I8))
 9033 FORMAT (' DW(1:M)     = ',5(F11.3)/(15X,5(F11.3)))
 9034 FORMAT (' DW(M+1:M+N) = ',5(F11.3)/(15X,5(F11.3)))
      END

C**********************************************************************
CCCCC LAST UPDATE Tue Nov 26 03:20:26 MET 2002
      SUBROUTINE MC64BD
     &(M,N,NE,IP,IRN,A,IPERM,NUM,JPERM,PR,Q,L,D,RINF) 
      IMPLICIT NONE
C
C *** Copyright (c) 2002  I.S. Duff and J. Koster                   ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C
      INTEGER M,N,NE,NUM
      INTEGER IP(N+1),IRN(NE),IPERM(M),JPERM(N),PR(N),Q(M),L(M)
      DOUBLE PRECISION A(NE)
      DOUBLE PRECISION D(M), RINF

C N, NE, IP, IRN are described in MC64A.
C A is a DOUBLE PRECISION array of length NE.
C   A(K), K=1..NE, must be set to the value of the entry
C   that corresponds to IRN(K). It is not altered.
C IPERM is an INTEGER array of length M. On exit, it contains the 
C    matching: IPERM(I) = 0 or row I is matched to column IPERM(I).
C NUM is INTEGER variable. On exit, it contains the cardinality of the
C    matching stored in IPERM.
C D is a REAL (or DOUBLE PRECISION) work array of length M.
C    On entry, D(1) contains the relaxation parameter RLX.
C RINF is the largest positive real number

C Local variables
      INTEGER I,II,J,JJ,JORD,Q0,QLEN,IDUM,JDUM,ISP,JSP,
     &        K,KK,KK1,KK2,I0,UP,LOW
      DOUBLE PRECISION CSP,DI,DNEW,DQ0,AI,A0,BV,TBV,RLX
C Local parameters 
      DOUBLE PRECISION ZERO,MINONE
      PARAMETER (ZERO=0.0D+0,MINONE=-1.0D+0)
C Intrinsic functions
      INTRINSIC ABS,MIN
C External subroutines and/or functions
      EXTERNAL MC64DD,MC64ED,MC64FD,MC64XD

      RLX = D(1)
C Initialization
      NUM = 0
      BV = RINF
      DO 10 K = 1,N
        JPERM(K) = 0
        PR(K) = IP(K)
   10 CONTINUE
      DO 12 K = 1,M
        IPERM(K) = 0
        D(K) = ZERO
   12 CONTINUE

      DO 30 J = 1,N
        A0 = MINONE
        DO 20 K = IP(J),IP(J+1)-1
          I = IRN(K)
          AI = ABS(A(K))
          IF (AI.GT.D(I)) D(I) = AI
          IF (JPERM(J).NE.0) GO TO 20
          IF (AI.GE.BV) THEN
            A0 = BV
            IF (IPERM(I).NE.0) GO TO 20
            JPERM(J) = I 
            IPERM(I) = J
            NUM = NUM + 1
          ELSE
            IF (AI.LE.A0) GO TO 20
            A0 = AI
            I0 = I
          ENDIF
   20   CONTINUE
        IF (A0.NE.MINONE .AND. A0.LT.BV) THEN
          BV = A0
          IF (IPERM(I0).NE.0) GO TO 30
          IPERM(I0) = J
          JPERM(J) = I0
          NUM = NUM + 1
        ENDIF
   30 CONTINUE

      IF (M.EQ.N) THEN
C Update BV with smallest of all the largest maximum absolute values 
C of the rows. D(I) contains the largest absolute value in row I.
        DO 35 I = 1,M
          BV = MIN(BV,D(I))
   35   CONTINUE
      ENDIF

C Shortcut if all columns are matched at this stage.
      IF (NUM.EQ.N) GO TO 1000

C Rescan unassigned columns; improve initial assignment
      DO 95 J = 1,N
        IF (JPERM(J).NE.0) GO TO 95
        DO 50 K = IP(J),IP(J+1)-1
          I = IRN(K)
          AI = ABS(A(K))
          IF (AI.LT.BV) GO TO 50
          IF (IPERM(I).EQ.0) GO TO 90
          JJ = IPERM(I)
          KK1 = PR(JJ)
          KK2 = IP(JJ+1) - 1
          IF (KK1.GT.KK2) GO TO 50
          DO 70 KK = KK1,KK2
            II = IRN(KK)
            IF (IPERM(II).NE.0) GO TO 70
            IF (ABS(A(KK)).GE.BV) GO TO 80
   70     CONTINUE
          PR(JJ) = KK2 + 1
   50   CONTINUE
        GO TO 95
   80   JPERM(JJ) = II
        IPERM(II) = JJ
        PR(JJ) = KK + 1
   90   NUM = NUM + 1
        JPERM(J) = I
        IPERM(I) = J
        PR(J) = K + 1
   95 CONTINUE

C Shortcut if all columns are matched at this stage.
      IF (NUM.EQ.N) GO TO 1000

C Prepare for main loop
      DO 99 I = 1,M
        D(I) = MINONE
        L(I) = 0
   99 CONTINUE
C TBV is a relaxed value of BV (ie TBV is slightly smaller than BV).
      TBV = BV * (1-RLX)

C Main loop ... each pass round this loop is similar to Dijkstra's 
C algorithm for solving the single source shortest path problem 

      DO 100 JORD = 1,N

        IF (JPERM(JORD).NE.0) GO TO 100
        QLEN = 0
        LOW = M + 1
        UP = M + 1
C CSP is cost of shortest path to any unassigned row
C ISP is matrix position of unassigned row element in shortest path
C JSP is column index of unassigned row element in shortest path
        CSP = MINONE
C Build shortest path tree starting from unassigned column JORD
        J = JORD
        PR(J) = -1

C Scan column J
        DO 115 K = IP(J),IP(J+1)-1
          I = IRN(K)
          DNEW = ABS(A(K))
          IF (CSP.GE.DNEW) GO TO 115
          IF (IPERM(I).EQ.0) THEN
C Row I is unassigned; update shortest path info
            CSP = DNEW
            ISP = I
            JSP = J
            IF (CSP.GE.TBV) GO TO 160
          ELSE
            D(I) = DNEW
            IF (DNEW.GE.TBV) THEN
C Add row I to Q2
              LOW = LOW - 1
              Q(LOW) = I
            ELSE
C Add row I to Q, and push it 
              QLEN = QLEN + 1
              L(I) = QLEN
              CALL MC64DD(I,M,Q,D,L,1)
            ENDIF
            JJ = IPERM(I)
            PR(JJ) = J
          ENDIF
  115   CONTINUE

        DO 150 JDUM = 1,NUM
C If Q2 is empty, extract new rows from Q
          IF (LOW.EQ.UP) THEN
            IF (QLEN.EQ.0) GO TO 160
            I = Q(1)
            IF (CSP.GE.D(I)) GO TO 160
            BV = D(I)
            TBV = BV * (1-RLX)
            DO 152 IDUM = 1,M
              CALL MC64ED(QLEN,M,Q,D,L,1)
              L(I) = 0
              LOW = LOW - 1
              Q(LOW) = I
              IF (QLEN.EQ.0) GO TO 153
              I = Q(1)
              IF (D(I).LT.TBV) GO TO 153
  152       CONTINUE
C End of dummy loop; this point is never reached
          ENDIF
C Move row Q0 
  153     UP = UP - 1
          Q0 = Q(UP)
          DQ0 = D(Q0)
          L(Q0) = UP
C Scan column that matches with row Q0
          J = IPERM(Q0)
          DO 155 K = IP(J),IP(J+1)-1
            I = IRN(K)
C Update D(I); only if row I is not marked
            IF (L(I).GE.UP) GO TO 155
            DNEW = MIN(DQ0,ABS(A(K)))
            IF (CSP.GE.DNEW) GO TO 155
            IF (IPERM(I).EQ.0) THEN
C Row I is unassigned; update shortest path info
              CSP = DNEW
              ISP = I
              JSP = J
              IF (CSP.GE.TBV) GO TO 160
            ELSE
              DI = D(I)
              IF (DI.GE.TBV .OR. DI.GE.DNEW) GO TO 155
              D(I) = DNEW
              IF (DNEW.GE.TBV) THEN
C Delete row I from Q (if necessary); add row I to Q2
                IF (DI.NE.MINONE) THEN
                  CALL MC64FD(L(I),QLEN,M,Q,D,L,1)
                ENDIF
                L(I) = 0
                LOW = LOW - 1
                Q(LOW) = I
              ELSE
C Add row I to Q (if necessary); push row I up Q
                IF (DI.EQ.MINONE) THEN
                  QLEN = QLEN + 1
                  L(I) = QLEN
                ENDIF
                CALL MC64DD(I,M,Q,D,L,1)
              ENDIF
C Update tree
              JJ = IPERM(I)
              PR(JJ) = J
            ENDIF
  155     CONTINUE
  150   CONTINUE
       
C If CSP = MINONE, no augmenting path is found
  160   IF (CSP.EQ.MINONE) GO TO 190
C Update bottleneck value
        BV = MIN(BV,CSP)
        TBV = BV * (1-RLX)
C Find augmenting path by tracing backward in PR; update IPERM,JPERM
        NUM = NUM + 1
        I = ISP
        J = JSP
        DO 170 JDUM = 1,NUM+1
          I0 = JPERM(J)
          JPERM(J) = I
          IPERM(I) = J
          J = PR(J)
          IF (J.EQ.-1) GO TO 190
          I = I0
  170   CONTINUE
C End of dummy loop; this point is never reached
  190   DO 191 KK = UP,M
          I = Q(KK)
          D(I) = MINONE
          L(I) = 0
  191   CONTINUE 
        DO 192 KK = LOW,UP-1
          I = Q(KK)
          D(I) = MINONE
  192   CONTINUE
        DO 193 KK = 1,QLEN
          I = Q(KK)
          D(I) = MINONE
          L(I) = 0
  193   CONTINUE

  100 CONTINUE
C End of main loop

C BV is now bottleneck value of final matching

C IPERM is complete if M = N and NUM = N
 1000 IF (M.EQ.N .and. NUM.EQ.N) GO TO 2000

C Complete IPERM; L, JPERM are work arrays
      CALL MC64XD(M,N,IPERM,L,JPERM)

 2000 RETURN
      END

C**********************************************************************
CCCCC LAST UPDATE Tue Nov 26 03:20:26 MET 2002
      SUBROUTINE MC64DD(I,N,Q,D,L,IWAY)
      IMPLICIT NONE
C
C *** Copyright (c) 2002  I.S. Duff and J. Koster                   ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C
      INTEGER I,N,IWAY
      INTEGER Q(N),L(N)
      DOUBLE PRECISION D(N)

C Variables N,Q,D,L are described in MC64B
C IF IWAY = 1, node I is pushed from its current position upwards
C IF IWAY /= 1, node I is pushed from its current position downwards

C Local variables and parameters
      INTEGER IDUM,K,POS,POSK,QK
      PARAMETER (K=2)
      DOUBLE PRECISION DI

      POS = L(I)
      IF (POS.LE.1) GO TO 20
      DI = D(I)
C POS is index of current position of I in the tree
      IF (IWAY.EQ.1) THEN
        DO 10 IDUM = 1,N
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.LE.D(QK)) GO TO 20 
          Q(POS) = QK
          L(QK) = POS 
          POS = POSK
          IF (POS.LE.1) GO TO 20
   10   CONTINUE
C End of dummy loop; this point is never reached
      ELSE
        DO 15 IDUM = 1,N
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.GE.D(QK)) GO TO 20
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
          IF (POS.LE.1) GO TO 20
   15   CONTINUE
C End of dummy loop; this point is never reached
      ENDIF
C End of dummy if; this point is never reached
   20 Q(POS) = I
      L(I) = POS

      RETURN
      END
 
C**********************************************************************
CCCCC LAST UPDATE Tue Nov 26 03:20:26 MET 2002
      SUBROUTINE MC64ED(QLEN,N,Q,D,L,IWAY)
      IMPLICIT NONE
C
C *** Copyright (c) 2002  I.S. Duff and J. Koster                   ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C
      INTEGER QLEN,N,IWAY
      INTEGER Q(N),L(N)
      DOUBLE PRECISION D(N)

C Variables QLEN,N,Q,D,L are described in MC64B (IWAY = 1) or
C     MC64W (IWAY = 2)
C The root node is deleted from the binary heap.

C Local variables and parameters
      INTEGER I,IDUM,K,POS,POSK
      PARAMETER (K=2)
      DOUBLE PRECISION DK,DR,DI

C Move last element to begin of Q
      I = Q(QLEN)
      DI = D(I)
      QLEN = QLEN - 1
      POS = 1
      IF (IWAY.EQ.1) THEN
        DO 10 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 20
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.LT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.GE.DK) GO TO 20
C Exchange old last element with larger priority child
          Q(POS) = Q(POSK)
          L(Q(POS)) = POS
          POS = POSK
   10   CONTINUE
C End of dummy loop; this point is never reached
      ELSE
        DO 15 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 20
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.GT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.LE.DK) GO TO 20
C Exchange old last element with smaller child
          Q(POS) = Q(POSK)
          L(Q(POS)) = POS
          POS = POSK
   15   CONTINUE
C End of dummy loop; this point is never reached
      ENDIF
C End of dummy if; this point is never reached
   20 Q(POS) = I
      L(I) = POS

      RETURN
      END

C**********************************************************************
CCCCC LAST UPDATE Tue Nov 26 03:20:26 MET 2002
      SUBROUTINE MC64FD(POS0,QLEN,N,Q,D,L,IWAY)
      IMPLICIT NONE
C
C *** Copyright (c) 2002  I.S. Duff and J. Koster                   ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C
      INTEGER POS0,QLEN,N,IWAY
      INTEGER Q(N),L(N)
      DOUBLE PRECISION D(N)

C Variables QLEN,N,Q,D,L are described in MC64B (IWAY = 1) or
C     MC64W (IWAY = 2).
C Move last element in the heap 

      INTEGER I,IDUM,K,POS,POSK,QK
      PARAMETER (K=2)
      DOUBLE PRECISION DK,DR,DI
 
C Quick return, if possible
      IF (QLEN.EQ.POS0) THEN
        QLEN = QLEN - 1
        RETURN
      ENDIF

C Move last element from queue Q to position POS0
C POS is current position of node I in the tree
      I = Q(QLEN)
      DI = D(I)
      QLEN = QLEN - 1
      POS = POS0
      IF (IWAY.EQ.1) THEN
        IF (POS.LE.1) GO TO 20
        DO 10 IDUM = 1,N
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.LE.D(QK)) GO TO 20 
          Q(POS) = QK
          L(QK) = POS 
          POS = POSK
          IF (POS.LE.1) GO TO 20
   10   CONTINUE
C End of dummy loop; this point is never reached
   20   Q(POS) = I
        L(I) = POS
        IF (POS.NE.POS0) RETURN
        DO 30 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 40
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.LT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.GE.DK) GO TO 40
          QK = Q(POSK)
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
   30   CONTINUE
C End of dummy loop; this point is never reached
      ELSE
        IF (POS.LE.1) GO TO 34
        DO 32 IDUM = 1,N
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.GE.D(QK)) GO TO 34 
          Q(POS) = QK
          L(QK) = POS 
          POS = POSK
          IF (POS.LE.1) GO TO 34
   32   CONTINUE
C End of dummy loop; this point is never reached
   34   Q(POS) = I
        L(I) = POS
        IF (POS.NE.POS0) RETURN
        DO 36 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 40
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.GT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.LE.DK) GO TO 40
          QK = Q(POSK)
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
   36   CONTINUE
C End of dummy loop; this point is never reached
      ENDIF
C End of dummy if; this point is never reached
   40 Q(POS) = I
      L(I) = POS

      RETURN
      END

C**********************************************************************
CCCCC LAST UPDATE Tue Nov 26 03:20:26 MET 2002
      SUBROUTINE MC64PD(M,N,NE,IP,IRN,A,MODE,PERM,IW,INFO)
      IMPLICIT NONE
C
C *** Copyright (c) 2002  I.S. Duff and J. Koster                   ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C
C Purpose 
C =======
C
C Permute the matrix defined by M,N,NE,IP,IRN,A according to the
C permutation defined by PERM.
C
C Parameters
C ==========
C
      INTEGER M,N,NE,MODE,INFO
      DOUBLE PRECISION A(NE)
      INTEGER IP(N+1),IRN(NE),PERM(*),IW(*)
C
C M is an INTEGER variable. 
C     On entry, it must be set by the user.
C     See MC64A for its meaning.
C     Restriction: M > 0.
C
C N is an INTEGER variable. 
C     On entry, it must be set by the user.
C     See MC64A for its meaning.
C     Restriction: N > 0.
C
C NE is an INTEGER variable. 
C     On entry, it must be set by the user.
C     See MC64A for its meaning.
C     Restriction: NE > 0.
C
C IP is an INTEGER array of length N+1. 
C     On entry, it must be set by the user.
C     See MC64A for its meaning.
C
C IRN is an INTEGER array of length NE. 
C     On entry, it must be set by the user.
C     See MC64A for its meaning.
C
C A is a DOUBLE PRECISION array of length NE.
C     On entry, it must be set by the user.
C     See MC64A for its meaning.
C
C MODE is an INTEGER variable.
C     On entry, it must be set by the user. It defines whether the
C     rows or columns of the matrix shall be permuted.
C     Possible values for MODE are:
C       1 : perform a row permutation of the matrix
C       2 : perform a column permutation of the matrix
C     Restriction: 1 <= MODE <= 2.
C
C PERM is an INTEGER array.
C     If MODE = 1, PERM is an array of length M.
C     On entry, PERM must define a row permutation of the matrix.
C     It must be set by the user such that row I of the original 
C     matrix is row PERM(I) in the permuted matrix, I=1..M. 
C
C     If MODE = 2, PERM is an array of length N.
C     On entry, PERM must define a column permutation of the matrix. 
C     It must be set by the user such that column PERM(J) of the
C     original matrix is column J in the permuted matrix, J=1..N.
C
C IW is an INTEGER array.
C     If MODE = 1, PERM is an array of length M.
C     If MODE = 2, PERM is an array of length N+NE.
C     On entry, it need not be set by the user.
C     It is used as workspace.
C     On exit, its contents are undefined.
C
C INFO is an INTEGER variable. 
C     On entry, it need not be set by the user.
C     On exit, INFO will take one of the following values:
C        0 : successful return
C       -2 : Invalid value for M or N  : M < 1 or N < 1
C       -3 : Invalid value for NE : NE < 1
C      -10 : Invalid value for MODE : MODE < 1 or MODE > 2
C      -11 : Permutation PERM is invalid.
C
C The subroutine does not write output to any message stream.

C Local variables and parameters
      INTEGER I, K

      EXTERNAL MC64YD,MC64ND

      INFO = 0
      IF (M.LE.0 .OR. N.LE.0) THEN
        INFO = -2
        GO TO 299
      ENDIF
      IF (NE.LE.0) THEN
        INFO = -3
        GO TO 299
      ENDIF
      IF (MODE.LT.1 .or. MODE.GT.2) THEN
        INFO = -10
        GO TO 299
      ENDIF
      IF (MODE.EQ.1) THEN
        CALL MC64YD(M,PERM,IW,INFO)
      ELSE
        CALL MC64YD(N,PERM,IW,INFO)
      ENDIF
      IF (INFO.NE.0) THEN
        INFO = -11
        GO TO 299
      ENDIF

C Input is OK, permute.

      IF (MODE.EQ.1) THEN
C perform a row permutation of the matrix
C 1. Check if PERM is an identity permutation
        K = 0
        DO 100 I = 1,M
          IF (PERM(I).EQ.I) GO TO 100
          K = 1
          GO TO 120
  100   CONTINUE
  120   IF (K.EQ.0) GO TO 299
C 2. Change the row indices
        DO 150 K = 1,NE
          I = IRN(K)
          IRN(K) = PERM(I) 
  150   CONTINUE
        GO TO 299
      ENDIF

C MODE == 2; perform a column permutation of the matrix
      CALL MC64ND(M,N,NE,IP,IRN,A,PERM,IW(1),IW(N+1))

  299 RETURN
      END

C**********************************************************************
CCCCC LAST UPDATE Tue Nov 26 03:20:26 MET 2002
      SUBROUTINE MC64ND(M,N,NE,IP,IRN,A,PERM,IW,DEST)
      IMPLICIT NONE
C
C *** Copyright (c) 2002  I.S. Duff and J. Koster                   ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C
      INTEGER M,N,NE
      DOUBLE PRECISION A(NE)
      INTEGER IP(N+1),IRN(NE),PERM(N),IW(N),DEST(NE)
C
C Local variables and parameters
      DOUBLE PRECISION AOLD,ANEW
      INTEGER I,J,K,IDUM,IOLD,INEW,I0,I1,JOLD,JNEW,LEN

C perform a column permutation of the matrix
C 1. Check if PERM is an identity permutation
      K = 0
      DO 200 J = 1,N
        IF (PERM(J).EQ.J) GO TO 200
        K = 1
        GO TO 220
  200 CONTINUE
  220 IF (K.EQ.0) GO TO 299
C 2. Compute for each index in IRN its index after permutation.
      K = 1
      DO 240 JNEW = 1,N
        JOLD = PERM(JNEW)
        LEN = IP(JOLD+1) - IP(JOLD)
        IW(JNEW) = LEN
        J = IP(JOLD)
        DO 230 I = 0,LEN-1
          DEST(J+I) = K+I
  230   CONTINUE
        K = K + LEN
  240 CONTINUE
C Compute IP for the permuted matrix from IW
      IP(1) = 1
      DO 250 JNEW = 2,N+1
        IP(JNEW) = IP(JNEW-1) + IW(JNEW-1)
  250 CONTINUE
C 3. Permute A and IRN in place.
C Each pass through this loop places a closed chain of row indices
C in their new positions ... this is recorded by setting the DEST 
C entry to zero so that any which can be encountered again can be
C bypassed.
      DO 290 I = 1,NE
        I0 = I
        I1 = DEST(I)
        IF (I1.EQ.0) GO TO 290
        IF (I1.EQ.I) GO TO 290
        IOLD = IRN(I0)
        AOLD = A(I0)
        DO 270 IDUM = 1,NE
          IF (I1.EQ.I) THEN
            IRN(I1) = IOLD
            A(I1) = AOLD
            GO TO 290
          ELSE
            INEW = IRN(I1)
            ANEW = A(I1)
            IRN(I1) = IOLD
            A(I1) = AOLD
            IOLD = INEW
            AOLD = ANEW
            I0 = I1
            I1 = DEST(I1)
            DEST(I0) = 0
          ENDIF
  270   CONTINUE
  290 CONTINUE
  
  299 RETURN
      END

C**********************************************************************
CCCCC LAST UPDATE Tue Nov 26 03:20:26 MET 2002
      SUBROUTINE MC64QD(IP,LENL,LENH,W,WLEN,A,NVAL,VAL)
      IMPLICIT NONE
C
C *** Copyright (c) 2002  I.S. Duff and J. Koster                   ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C
      INTEGER WLEN,NVAL
      INTEGER IP(*),LENL(*),LENH(*),W(*)
      DOUBLE PRECISION A(*),VAL

C This routine searches for at most XX different numerical values 
C in the columns W(1:WLEN). XX>=2. 
C Each column J is scanned between IP(J)+LENL(J) and IP(J)+LENH(J)-1 
C until XX values are found or all columns have been considered.
C On output, NVAL is the number of different values that is found 
C and SPLIT(1:NVAL) contains the values in decreasing order.
C If NVAL > 0, the routine returns VAL = SPLIT((NVAL+1)/2). 
C
      INTEGER XX,J,K,II,S,POS
      PARAMETER (XX=10)
      DOUBLE PRECISION SPLIT(XX),HA

C Scan columns in W(1:WLEN). For each encountered value, if value not
C already present in SPLIT(1:NVAL), insert value such that SPLIT 
C remains sorted by decreasing value. 
C The sorting is done by straightforward insertion; therefore the use
C of this routine should be avoided for large XX (XX < 20).
      NVAL = 0 
      DO 10 K = 1,WLEN
        J = W(K)
        DO 15 II = IP(J)+LENL(J),IP(J)+LENH(J)-1
          HA = A(II)
          IF (NVAL.EQ.0) THEN
            SPLIT(1) = HA
            NVAL = 1
          ELSE
C Check presence of HA in SPLIT
            DO 20 S = NVAL,1,-1
              IF (SPLIT(S).EQ.HA) GO TO 15
              IF (SPLIT(S).GT.HA) THEN
                POS = S + 1
                GO TO 21
              ENDIF
  20        CONTINUE
            POS = 1
C The insertion 
  21        DO 22 S = NVAL,POS,-1
              SPLIT(S+1) = SPLIT(S)
  22        CONTINUE
            SPLIT(POS) = HA
            NVAL = NVAL + 1
          ENDIF
C Exit loop if XX values are found
          IF (NVAL.EQ.XX) GO TO 11
  15    CONTINUE
  10  CONTINUE
C Determine VAL
  11  IF (NVAL.GT.0) VAL = SPLIT((NVAL+1)/2)

      RETURN
      END

C**********************************************************************
CCCCC LAST UPDATE Tue Nov 26 03:20:26 MET 2002
      SUBROUTINE MC64RD(N,NE,IP,IRN,A)
      IMPLICIT NONE
C
C *** Copyright (c) 2002  I.S. Duff and J. Koster                   ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C
      INTEGER N,NE
      INTEGER IP(N+1),IRN(NE)
      DOUBLE PRECISION A(NE)

C This subroutine sorts the entries in each column of the 
C sparse matrix (defined by N,NE,IP,IRN,A) by decreasing
C numerical value.

C Local constants
      INTEGER THRESH,TDLEN
      PARAMETER (THRESH=15,TDLEN=50)
C Local variables
      INTEGER J,IPJ,K,LEN,R,S,HI,FIRST,MID,LAST,TD
      DOUBLE PRECISION HA,KEY
C Local arrays
      INTEGER TODO(TDLEN)
      
      DO 100 J = 1,N
        LEN = IP(J+1) - IP(J)
        IF (LEN.LE.1) GO TO 100
        IPJ = IP(J)

C Sort array roughly with partial quicksort
        IF (LEN.LT.THRESH) GO TO 400
        TODO(1) = IPJ
        TODO(2) = IPJ + LEN
        TD = 2
  500   CONTINUE
        FIRST = TODO(TD-1)
        LAST = TODO(TD)
C KEY is the smallest of two values present in interval [FIRST,LAST)
        KEY = A((FIRST+LAST)/2)
        DO 475 K = FIRST,LAST-1
          HA = A(K)
          IF (HA.EQ.KEY) GO TO 475
          IF (HA.GT.KEY) GO TO 470
          KEY = HA
          GO TO 470
  475   CONTINUE
C Only one value found in interval, so it is already sorted
        TD = TD - 2
        GO TO 425

C Reorder interval [FIRST,LAST) such that entries before MID are > KEY
  470   MID = FIRST
        DO 450 K = FIRST,LAST-1
          IF (A(K).LE.KEY) GO TO 450
          HA = A(MID)
          A(MID) = A(K)
          A(K) = HA
          HI = IRN(MID)
          IRN(MID) = IRN(K)
          IRN(K) = HI
          MID = MID + 1
  450   CONTINUE
C Both subintervals [FIRST,MID), [MID,LAST) are nonempty
C Stack the longest of the two subintervals first
        IF (MID-FIRST.GE.LAST-MID) THEN
          TODO(TD+2) = LAST
          TODO(TD+1) = MID
          TODO(TD) = MID
C          TODO(TD-1) = FIRST
        ELSE
          TODO(TD+2) = MID
          TODO(TD+1) = FIRST
          TODO(TD) = LAST
          TODO(TD-1) = MID
        ENDIF
        TD = TD + 2

  425   CONTINUE
        IF (TD.EQ.0) GO TO 400 
C There is still work to be done
        IF (TODO(TD)-TODO(TD-1).GE.THRESH) GO TO 500
C Next interval is already short enough for straightforward insertion
        TD = TD - 2
        GO TO 425

C Complete sorting with straightforward insertion
  400   DO 200 R = IPJ+1,IPJ+LEN-1
          IF (A(R-1) .LT. A(R)) THEN
            HA = A(R)
            HI = IRN(R)
            A(R) = A(R-1)
            IRN(R) = IRN(R-1)
            DO 300 S = R-1,IPJ+1,-1
              IF (A(S-1) .LT. HA) THEN
                A(S) = A(S-1)
                IRN(S) = IRN(S-1)
              ELSE
                A(S) = HA
                IRN(S) = HI
                GO TO 200 
              END IF
  300       CONTINUE
            A(IPJ) = HA
            IRN(IPJ) = HI
          END IF
  200   CONTINUE

  100 CONTINUE

      RETURN
      END

C**********************************************************************
CCCCC LAST UPDATE Tue Nov 26 03:20:26 MET 2002
      SUBROUTINE MC64SD(M,N,NE,IP,IRN,A,IPERM,NUMX,
     &           W,LEN,LENL,LENH,FC,IW,IW4,RLX,RINF)
      IMPLICIT NONE
C
C *** Copyright (c) 2002  I.S. Duff and J. Koster                   ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C
      INTEGER M,N,NE,NUMX
      INTEGER IP(N+1),IRN(NE),IPERM(N), 
     &        W(N),LEN(N),LENL(N),LENH(N),FC(N),IW(M),IW4(3*N+M)
      DOUBLE PRECISION A(NE),RLX,RINF

C M, N, NE, IP, IRN, are described in MC64A.
C A is a REAL (or DOUBLE PRECISION) array of length NE.
C   A(K), K=1..NE, must be set to the value of the entry that 
C   corresponds to IRN(k). The entries in each column must be 
C   non-negative and ordered by decreasing value. 
C IPERM is an INTEGER array of length M. On exit, it contains the 
C   bottleneck matching: IPERM(I) - 0 or row I is matched to column 
C   IPERM(I).
C NUMX is an INTEGER variable. On exit, it contains the cardinality 
C   of the matching stored in IPERM.

C FC is an integer array of length N that contains the list of 
C   unmatched columns. 
C LEN(J), LENL(J), LENH(J) are integer arrays of length N that point 
C   to entries in matrix column J.
C   In the matrix defined by the column parts IP(J)+LENL(J) we know 
C   a matching does not exist; in the matrix defined by the column 
C   parts IP(J)+LENH(J) we know one exists.
C   LEN(J) lies between LENL(J) and LENH(J) and determines the matrix
C   that is tested for a maximum matching.
C W is an integer array of length N and contains the indices of the 
C   columns for which LENL /= LENH.
C WLEN is number of indices stored in array W.
C IW is integer work array of length M.
C IW4 is integer work array of length 3N+M used by MC64U.
C
C RLX is a REAL (or DOUBLE PRECISION) variable. It is a relaxation 
C   parameter for finding the optimal matching.
C
C RINF is the largest positive real number

      INTEGER NUM,NVAL,WLEN,II,I,J,K,L,CNT,MOD,IDUM1,IDUM2,IDUM3
      DOUBLE PRECISION BVAL,BMIN,BMAX
C External subroutines and/or functions
      EXTERNAL MC64QD,MC64UD,MC64XD
C Intrinsic functions

C BMIN and BMAX are such that a maximum matching exists for the input
C   matrix in which all entries smaller than BMIN are dropped. 
C   For BMAX, a maximum matching does not exist.
C BVAL is a value between BMIN and BMAX.
C CNT is the number of calls made to MC64U so far.
C NUM is the cardinality of last matching found.

C Compute a first maximum matching from scratch on whole matrix.
      DO 20 J = 1,N
        FC(J) = J
        LEN(J) = IP(J+1) - IP(J)
   20 CONTINUE
      DO 21 I = 1,M
        IW(I) = 0
   21 CONTINUE
C The first call to MC64U
      CNT = 1
      MOD = 1
      NUMX = 0
      CALL MC64UD(CNT,MOD,M,N,IRN,NE,IP,LEN,FC,IW,NUMX,N,
     &            IW4(1),IW4(N+1),IW4(2*N+1),IW4(2*N+M+1))

C IW contains a maximum matching of length NUMX.
      NUM = NUMX

      IF (NUM.NE.N) THEN
C Matrix is structurally singular
        BMAX = RINF
      ELSE
C Matrix is structurally nonsingular, NUM=NUMX=N;
C Set BMAX just above the smallest of all the maximum absolute 
C values of the columns
        BMAX = RINF
        DO 30 J = 1,N
          BVAL = 0.0
          DO 25 K = IP(J),IP(J+1)-1
            IF (A(K).GT.BVAL) BVAL = A(K)
   25     CONTINUE
          IF (BVAL.LT.BMAX) BMAX = BVAL
   30   CONTINUE
C ... should print warning if BMAX == RINF
        BMAX = 1.001 * BMAX
      ENDIF

C Initialize BVAL,BMIN
      BVAL = 0.0
      BMIN = 0.0
C Initialize LENL,LEN,LENH,W,WLEN according to BMAX.
C Set LEN(J), LENH(J) just after last entry in column J. 
C Set LENL(J) just after last entry in column J with value >= BMAX.
      WLEN = 0
      DO 48 J = 1,N
        L = IP(J+1) - IP(J)
        LENH(J) = L
        LEN(J) = L
        DO 45 K = IP(J),IP(J+1)-1
          IF (A(K).LT.BMAX) GO TO 46
   45   CONTINUE
C Column J is empty or all entries are >= BMAX
        K = IP(J+1)
   46   LENL(J) = K - IP(J)
C Add J to W if LENL(J) /= LENH(J)
        IF (LENL(J).EQ.L) GO TO 48
        WLEN = WLEN + 1
        W(WLEN) = J
   48 CONTINUE

C Main loop
      DO 90 IDUM1 = 1,NE
        IF (NUM.EQ.NUMX) THEN
C We have a maximum matching in IW; store IW in IPERM
          DO 50 I = 1,M
            IPERM(I) = IW(I)
   50     CONTINUE
C Keep going round this loop until matching IW is no longer maximum.
          DO 80 IDUM2 = 1,NE
            BMIN = BVAL
            IF (BMAX-BMIN .LE. RLX) GO TO 1000
C Find splitting value BVAL
            CALL MC64QD(IP,LENL,LEN,W,WLEN,A,NVAL,BVAL)
            IF (NVAL.LE.1) GO TO 1000
C Set LEN such that all matrix entries with value < BVAL are 
C discarded. Store old LEN in LENH. Do this for all columns W(K).
C Each step, either K is incremented or WLEN is decremented.
            K = 1
            DO 70 IDUM3 = 1,N
              IF (K.GT.WLEN) GO TO 71
              J = W(K)
              DO 55 II = IP(J)+LEN(J)-1,IP(J)+LENL(J),-1
                IF (A(II).GE.BVAL) GO TO 60 
                I = IRN(II)
                IF (IW(I).NE.J) GO TO 55
C Remove entry from matching
                IW(I) = 0
                NUM = NUM - 1
                FC(N-NUM) = J
   55         CONTINUE
   60         LENH(J) = LEN(J)
C IP(J)+LEN(J)-1 is last entry in column >= BVAL
              LEN(J) = II - IP(J) + 1
C If LENH(J) = LENL(J), remove J from W
              IF (LENL(J).EQ.LENH(J)) THEN
                W(K) = W(WLEN)
                WLEN = WLEN - 1
              ELSE
                K = K + 1
              ENDIF
   70       CONTINUE
   71       IF (NUM.LT.NUMX) GO TO 81
   80     CONTINUE
C End of dummy loop; this point is never reached
C Set mode for next call to MC64U
   81     MOD = 1
        ELSE
C We do not have a maximum matching in IW. 
          BMAX = BVAL
C BMIN is the bottleneck value of a maximum matching;
C for BMAX the matching is not maximum, so BMAX>BMIN
C and following condition is always false if RLX = 0.0
          IF (BMAX-BMIN .LE. RLX) GO TO 1000
C Find splitting value BVAL
          CALL MC64QD(IP,LEN,LENH,W,WLEN,A,NVAL,BVAL)
          IF (NVAL.EQ.0. OR. BVAL.EQ.BMIN) GO TO 1000
C Set LEN such that all matrix entries with value >= BVAL are
C inside matrix. Store old LEN in LENL. Do this for all columns W(K).
C Each step, either K is incremented or WLEN is decremented.
          K = 1
          DO 87 IDUM3 = 1,N
            IF (K.GT.WLEN) GO TO 88
            J = W(K)
            DO 85 II = IP(J)+LEN(J),IP(J)+LENH(J)-1
              IF (A(II).LT.BVAL) GO TO 86
   85       CONTINUE
   86       LENL(J) = LEN(J)
            LEN(J) = II - IP(J)
            IF (LENL(J).EQ.LENH(J)) THEN
              W(K) = W(WLEN)
              WLEN = WLEN - 1
            ELSE
              K = K + 1
            ENDIF
   87     CONTINUE
C End of dummy loop; this point is never reached
C Set mode for next call to MC64U
   88     MOD = 0
        ENDIF
        CNT = CNT + 1
        CALL MC64UD(CNT,MOD,M,N,IRN,NE,IP,LEN,FC,IW,NUM,NUMX,
     &              IW4(1),IW4(N+1),IW4(2*N+1),IW4(2*N+M+1))

C IW contains maximum matching of length NUM
   90 CONTINUE 
C End of dummy loop; this point is never reached

C BMIN is now bottleneck value of final matching

C IPERM is complete if M = N and NUMX = N
 1000 IF (M.EQ.N .and. NUMX.EQ.N) GO TO 2000

C Complete IPERM; IW, W are work arrays
      CALL MC64XD(M,N,IPERM,IW,W)

 2000 RETURN
      END
      
C**********************************************************************
CCCCC LAST UPDATE Tue Nov 26 03:20:26 MET 2002
      SUBROUTINE MC64UD
     &           (ID,MOD,M,N,IRN,LIRN,IP,LENC,FC,IPERM,NUM,NUMX,
     &           PR,ARP,CV,OUT)
      IMPLICIT NONE
C
C *** Copyright (c) 2002  I.S. Duff and J. Koster                   ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C
      INTEGER ID,MOD,M,N,LIRN,NUM,NUMX
      INTEGER ARP(N),CV(M),IRN(LIRN),IP(N),
     &        FC(N),IPERM(M),LENC(N),OUT(N),PR(N)

C PR(J) is the previous column to J in the depth first search.
C   Array PR is used as workspace in the sorting algorithm.
C Elements (I,IPERM(I)) I=1,..,M are entries at the end of the
C   algorithm unless N assignments have not been made in which case
C   N-NUM pairs (I,IPERM(I)) will not be entries in the matrix.
C CV(I) is the most recent loop number (ID+JORD) at which row I
C   was visited.
C ARP(J) is the number of entries in column J which have been scanned 
C   when looking for a cheap assignment.
C OUT(J) is one less than the number of entries in column J which have 
C   not been scanned during one pass through the main loop.
C NUMX is maximum possible size of matching.

      INTEGER I,II,IN1,IN2,J,J1,JORD,K,KK,LAST,NFC,
     &        NUM0,NUM1,NUM2,ID0,ID1

      IF (ID.EQ.1) THEN
C The first call to MC64U.
C Initialize CV and ARP; parameters MOD, NUMX are not accessed
        DO 5 I = 1,M
          CV(I) = 0
    5   CONTINUE
        DO 6 J = 1,N
          ARP(J) = 0
    6   CONTINUE
        NUM1 = N
        NUM2 = N
      ELSE
C Not the first call to MC64U.
C Re-initialize ARP if entries were deleted since last call to MC64U
        IF (MOD.EQ.1) THEN
          DO 8 J = 1,N
            ARP(J) = 0
    8     CONTINUE
        ENDIF
        NUM1 = NUMX
        NUM2 = N - NUMX
      ENDIF
      NUM0 = NUM

C NUM0 is size of input matching
C NUM1 is maximum possible size of matching
C NUM2 is maximum allowed number of unassigned rows/columns
C NUM is size of current matching

C Quick return if possible
C      IF (NUM.EQ.N) GO TO 199
C NFC is number of rows/columns that could not be assigned
      NFC = 0
C Integers ID0+1 to ID0+N are unique numbers for call ID to MC64U,
C so 1st call uses 1..N, 2nd call uses N+1..2N, etc
      ID0 = (ID-1)*N 

C Main loop. Each pass round this loop either results in a new
C assignment or gives a column with no assignment

      DO 100 JORD = NUM0+1,N

C Each pass uses unique number ID1
        ID1 = ID0 + JORD
C J is unmatched column
        J = FC(JORD-NUM0)
        PR(J) = -1
        DO 70 K = 1,JORD
C Look for a cheap assignment
          IF (ARP(J).GE.LENC(J)) GO TO 30
          IN1 = IP(J) + ARP(J)
          IN2 = IP(J) + LENC(J) - 1
          DO 20 II = IN1,IN2
            I = IRN(II)
            IF (IPERM(I).EQ.0) GO TO 80
   20     CONTINUE
C No cheap assignment in row
          ARP(J) = LENC(J)
C Begin looking for assignment chain starting with row J
   30     OUT(J) = LENC(J) - 1
C Inner loop.  Extends chain by one or backtracks
          DO 60 KK = 1,JORD
            IN1 = OUT(J)
            IF (IN1.LT.0) GO TO 50
            IN2 = IP(J) + LENC(J) - 1
            IN1 = IN2 - IN1
C Forward scan
            DO 40 II = IN1,IN2
              I = IRN(II)
              IF (CV(I).EQ.ID1) GO TO 40
C Column J has not yet been accessed during this pass
              J1 = J
              J = IPERM(I)
              CV(I) = ID1
              PR(J) = J1
              OUT(J1) = IN2 - II - 1
              GO TO 70
   40       CONTINUE
C Backtracking step.
   50       J1 = PR(J)
            IF (J1.EQ.-1) THEN
C No augmenting path exists for column J.
              NFC = NFC + 1
              FC(NFC) = J
              IF (NFC.GT.NUM2) THEN
C A matching of maximum size NUM1 is not possible
                LAST = JORD
                GO TO 101
              ENDIF
              GO TO 100
            ENDIF
            J = J1
   60     CONTINUE
C End of dummy loop; this point is never reached
   70   CONTINUE
C End of dummy loop; this point is never reached

C New assignment is made.
   80   IPERM(I) = J
        ARP(J) = II - IP(J) + 1
        NUM = NUM + 1
        DO 90 K = 1,JORD
          J = PR(J)
          IF (J.EQ.-1) GO TO 95
          II = IP(J) + LENC(J) - OUT(J) - 2
          I = IRN(II)
          IPERM(I) = J
   90   CONTINUE
C End of dummy loop; this point is never reached

   95   IF (NUM.EQ.NUM1) THEN
C A matching of maximum size NUM1 is found
          LAST = JORD
          GO TO 101
        ENDIF
C
  100 CONTINUE

C All unassigned columns have been considered
      LAST = N

C Now, a transversal is computed or is not possible.
C Complete FC before returning.
  101 DO 110 JORD = LAST+1,N
        NFC = NFC + 1
        FC(NFC) = FC(JORD-NUM0)
  110 CONTINUE

C  199 RETURN
      RETURN
      END

C**********************************************************************
CCCCC LAST UPDATE Tue Nov 26 03:20:26 MET 2002
      SUBROUTINE MC64WD(M,N,NE,IP,IRN,A,IPERM,NUM,
     &           JPERM,OUT,PR,Q,L,U,D,RINF) 
      IMPLICIT NONE
C
C *** Copyright (c) 2002  I.S. Duff and J. Koster                   ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C
      INTEGER M,N,NE,NUM
      INTEGER IP(N+1),IRN(NE),IPERM(M),JPERM(N),OUT(N),PR(N),Q(M),L(M)
      DOUBLE PRECISION A(NE),U(M),D(M),RINF,RINF3

C M, N, NE, IP, IRN are described in MC64A.
C A is a REAL (or DOUBLE PRECISION) array of length NE.
C   A(K), K=1..NE, must be set to the value of the entry that
C   corresponds to IRN(K). It is not altered.
C   All values A(K) must be non-negative.
C IPERM is an INTEGER array of length M. On exit, it contains the 
C   weighted matching: IPERM(I) = 0 or row I is matched to column 
C   IPERM(I). 
C NUM is an INTEGER variable. On exit, it contains the cardinality of 
C   the matching stored in IPERM.
C D is a REAL (or DOUBLE PRECISION) array of length M.
C   On exit, V = D(1:N) contains the dual column variable.
C   If U(1:M) denotes the dual row variable and if the matrix 
C   is structurally nonsingular (NUM = N), the following holds:
C      U(I)+V(J) <= A(I,J)  if IPERM(I) /= J
C      U(I)+V(J)  = A(I,J)  if IPERM(I)  = J
C      U(I) = 0  if IPERM(I) = 0
C      V(J) = 0  if there is no I for which IPERM(I) = J
C   On entry, U(1) contains the relaxation parameter RLX.
C RINF is the largest positive real number

C Local variables
      INTEGER I,I0,II,J,JJ,JORD,Q0,QLEN,JDUM,ISP,JSP,
     &        K,K0,K1,K2,KK,KK1,KK2,UP,LOW
      DOUBLE PRECISION CSP,DI,DMIN,DNEW,DQ0,VJ,RLX
      LOGICAL LORD
C Local parameters
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
C External subroutines and/or functions
      EXTERNAL MC64DD,MC64ED,MC64FD,MC64XD

      RLX = U(1)
      RINF3 = U(2)
      LORD = (JPERM(1).EQ.6)
C Initialization
      NUM = 0
      DO 10 K = 1,N
        JPERM(K) = 0
        PR(K) = IP(K)
        D(K) = RINF
   10 CONTINUE
      DO 15 K = 1,M
        U(K) = RINF3
        IPERM(K) = 0
        L(K) = 0
   15 CONTINUE
      DO 30 J = 1,N
CSP QD
C     Don't choose cheap assignment from dense columns
         IF (IP(J+1)-IP(J) .GT. N/10 .AND. N.GT.50) GO TO 30
CSPEND QD
        DO 20 K = IP(J),IP(J+1)-1
          I = IRN(K)
          IF (A(K).GT.U(I)) GO TO 20
          U(I) = A(K)
          IPERM(I) = J
          L(I) = K
   20   CONTINUE
   30 CONTINUE
      DO 40 I = 1,M
        J = IPERM(I)
        IF (J.EQ.0) GO TO 40
        IF (JPERM(J).EQ.0) THEN
          JPERM(J) = L(I)
          D(J) = U(I)
          NUM = NUM + 1
        ELSEIF (D(J).GT.U(I)) THEN
          K = JPERM(J)
          II = IRN(K)
          IPERM(II) = 0
          JPERM(J) = L(I)
          D(J) = U(I)
        ELSE
          IPERM(I) = 0
        ENDIF
   40 CONTINUE

C      WRITE(*,*) 'lazy assignement 1 done NUM = ',NUM
#if 0
      DO I=1,N
         IF(JPERM(I) .NE. 0) THEN
            J = IRN(JPERM(I))
            IF(D(I) .GT. LOG(HUGE(RINF))) THEN
               write(*,*) 'warning I,D(I),U(I)',I,D(I),U(I)
            ENDIF
         ENDIF
         IF(J .EQ. I) THEN
            write(*,*) 'diag in lazy I = ',I
         ENDIF
      ENDDO
#endif
      IF (NUM.EQ.N) GO TO 1000

      DO 45 K = 1,M
        D(K) = ZERO
   45 CONTINUE
C Scan unassigned columns; improve assignment
      DO 95 J = 1,N
C JPERM(J) /= 0 iff column J is already assigned
        IF (JPERM(J).NE.0) GO TO 95
        K1 = IP(J)
        K2 = IP(J+1) - 1
C Continue only if column J is not empty
        IF (K1.GT.K2) GO TO 95
        VJ = RINF
        DO 50 K = K1,K2
          I = IRN(K)
          DI = A(K) - U(I)
          IF (DI.GT.VJ) GO TO 50
          IF (DI.LT.VJ .OR. DI.EQ.RINF) GO TO 55
          IF (IPERM(I).NE.0 .OR. IPERM(I0).EQ.0) GO TO 50
   55     VJ = DI
          I0 = I
          K0 = K
   50   CONTINUE
        D(J) = VJ
        K = K0
        I = I0
        IF (IPERM(I).EQ.0) GO TO 90
        DO 60 K = K0,K2
          I = IRN(K)
          IF (A(K)-U(I).GT.VJ) GO TO 60 
          JJ = IPERM(I)
C Scan remaining part of assigned column JJ 
          KK1 = PR(JJ)
          KK2 = IP(JJ+1) - 1
          IF (KK1.GT.KK2) GO TO 60
          DO 70 KK = KK1,KK2
            II = IRN(KK)
            IF (IPERM(II).GT.0) GO TO 70
            IF (A(KK)-U(II).LE.D(JJ)) GO TO 80
   70     CONTINUE
          PR(JJ) = KK2 + 1
   60   CONTINUE
        GO TO 95
   80   JPERM(JJ) = KK
        IPERM(II) = JJ
        PR(JJ) = KK + 1
   90   NUM = NUM + 1
        JPERM(J) = K
        IPERM(I) = J
        PR(J) = K + 1
   95 CONTINUE

C      WRITE(*,*) 'lazy assignement 2 done NUM = ',NUM
      IF (NUM.EQ.N) GO TO 1000
  
C Prepare for main loop
      DO 99 I = 1,M
        D(I) = RINF
        L(I) = 0
   99 CONTINUE

C Main loop ... each pass round this loop is similar to Dijkstra's 
C algorithm for solving the single source shortest path problem 

      DO 100 JORD = 1,N
        IF (JPERM(JORD).NE.0) GO TO 100
C JORD is next unmatched column
C DMIN is the length of shortest path in the tree
        DMIN = RINF
        QLEN = 0
        LOW = M + 1
        UP = M + 1
C CSP is the cost of the shortest augmenting path to unassigned row
C IRN(ISP). The corresponding column index is JSP.
        CSP = RINF
C Build shortest path tree starting from unassigned column (root) JORD
        J = JORD
        PR(J) = -1

C Scan column J
        DO 115 K = IP(J),IP(J+1)-1
          I = IRN(K)
          DNEW = A(K) - U(I)
          IF (DNEW.GE.CSP) GO TO 115
          IF (IPERM(I).EQ.0) THEN
            CSP = DNEW
            ISP = K
            JSP = J
          ELSE
            IF (DNEW.LT.DMIN) DMIN = DNEW
            D(I) = DNEW
            QLEN = QLEN + 1
            Q(QLEN) = K
          ENDIF
  115   CONTINUE
C Initialize heap and Q2 with rows held in Q(1:QLEN)
        Q0 = QLEN
        QLEN = 0
        DO 120 KK = 1,Q0
          K = Q(KK)
          I = IRN(K)
          IF (CSP.LE.D(I)) THEN
            D(I) = RINF
            GO TO 120
          ENDIF
          IF (D(I).LE.DMIN) THEN
            LOW = LOW - 1
            Q(LOW) = I
            L(I) = LOW
          ELSE
            QLEN = QLEN + 1
            L(I) = QLEN
            CALL MC64DD(I,M,Q,D,L,2)
          ENDIF
C Update tree
          JJ = IPERM(I)
          OUT(JJ) = K
          PR(JJ) = J
  120   CONTINUE

        DO 150 JDUM = 1,NUM

C If Q2 is empty, extract rows from Q 
          IF (LOW.EQ.UP) THEN
            IF (QLEN.EQ.0) GO TO 160
            I = Q(1)
            IF (D(I).LT.RINF) DMIN = D(I)*(1+RLX)
            IF (DMIN.GE.CSP) GO TO 160
  152       CALL MC64ED(QLEN,M,Q,D,L,2)
            LOW = LOW - 1
            Q(LOW) = I
            L(I) = LOW
            IF (QLEN.EQ.0) GO TO 153
            I = Q(1)
            IF (D(I).GT.DMIN) GO TO 153
            GO TO 152
          ENDIF
C Q0 is a row whose distance D(Q0) to the root is smallest
  153     Q0 = Q(UP-1)
          DQ0 = D(Q0)
C Exit loop if path to Q0 is longer than the shortest augmenting path 
          IF (DQ0.GE.CSP) GO TO 160
          IF (DMIN.GE.CSP) GO TO 160
C Move row Q0 from Q2 to Q3 
          UP = UP - 1

          J = IPERM(Q0)
          VJ = DQ0 - A(JPERM(J)) + U(Q0)
          K1 = IP(J+1)-1
          IF (LORD) THEN
            IF (CSP.NE.RINF) THEN
              DI = CSP - VJ
              IF (A(K1).GE.DI) THEN
                K0 = JPERM(J)
                IF (K0.GE.K1-6) GO TO 178
  177           CONTINUE
                  K = (K0+K1)/2
                  IF (A(K).GE.DI) THEN 
                    K1 = K
                  ELSE 
                    K0 = K
                  ENDIF
                  IF (K0.GE.K1-6) GO TO 178
                GO TO 177
  178           DO 179 K = K0+1,K1          
                  IF (A(K).LT.DI) GO TO 179
                  K1 = K - 1
                  GO TO 181
  179           CONTINUE
              ENDIF
            ENDIF
  181       IF (K1.EQ.JPERM(J)) K1 = K1 - 1
          ENDIF
 
          K0 = IP(J)
          DI = CSP - VJ
          DO 155 K = K0,K1
            I = IRN(K)
            IF (L(I).GE.LOW) GO TO 155
            DNEW = A(K) - U(I)
            IF (DNEW.GE.DI) GO TO 155
            DNEW = DNEW + VJ
            IF (DNEW.GT.D(I)) GO TO 155
            IF (IPERM(I).EQ.0) THEN
              CSP = DNEW
              ISP = K
              JSP = J
              DI = CSP - VJ
            ELSE
              IF (DNEW.GE.D(I)) GO TO 155
              D(I) = DNEW
              IF (DNEW.LE.DMIN) THEN
                IF (L(I).NE.0) THEN
                  CALL MC64FD(L(I),QLEN,M,Q,D,L,2)
                ENDIF
                LOW = LOW - 1
                Q(LOW) = I
                L(I) = LOW
              ELSE   
                IF (L(I).EQ.0) THEN
                  QLEN = QLEN + 1
                  L(I) = QLEN
                ENDIF
                CALL MC64DD(I,M,Q,D,L,2)
              ENDIF
              JJ = IPERM(I)
              OUT(JJ) = K
              PR(JJ) = J
            ENDIF
  155     CONTINUE
  150   CONTINUE
       
C If CSP = RINF, no augmenting path is found
  160   IF (CSP.EQ.RINF) GO TO 190
C Find augmenting path by tracing backward in PR; update IPERM,JPERM
        NUM = NUM + 1
        I = IRN(ISP)
        J = JSP
        IPERM(I) = J
        JPERM(J) = ISP
        DO 170 JDUM = 1,NUM
          JJ = PR(J) 
          IF (JJ.EQ.-1) GO TO 180
          K = OUT(J)
          I = IRN(K)
          IPERM(I) = JJ
          JPERM(JJ) = K
          J = JJ
  170   CONTINUE
C End of dummy loop; this point is never reached

C Reset D,L for rows in Q(LOW:UP-1)
C Reset D,L and update U for rows in Q(UP:M)
C Reset D,L for rows in binary heap Q(1:QLEN)
  180   DO 182 KK = UP,M
          I = Q(KK)
          U(I) = U(I) + D(I) - CSP
  182   CONTINUE 
  190   DO 191 KK = UP,M
          I = Q(KK)
          D(I) = RINF
          L(I) = 0
  191   CONTINUE
        DO 192 KK = LOW,UP-1
          I = Q(KK)
          D(I) = RINF
          L(I) = 0
  192   CONTINUE 
        DO 193 KK = 1,QLEN
          I = Q(KK)
          D(I) = RINF
          L(I) = 0
  193   CONTINUE

  100 CONTINUE
C End of main loop

 1000 CONTINUE
      DO 1200 J = 1,N
        K = JPERM(J)
        IF (K.NE.0) THEN
          D(J) = A(K) - U(IRN(K))
        ELSE
          D(J) = ZERO
        ENDIF
 1200 CONTINUE
      DO 1201 I = 1,M
        IF (IPERM(I).EQ.0) U(I) = ZERO
 1201 CONTINUE

C IPERM is complete if M = N and NUM = N
      IF (M.EQ.N .and. NUM.EQ.N) GO TO 2000

C Complete IPERM; L, JPERM are work arrays
      CALL MC64XD(M,N,IPERM,L,JPERM)

 2000 RETURN
      END
 
C**********************************************************************
CCCCC LAST UPDATE Tue Nov 26 03:20:26 MET 2002
      SUBROUTINE MC64ZD
     &           (M,N,IRN,LIRN,IP,LENC,IPERM,NUM,PR,ARP,CV,OUT)
      IMPLICIT NONE
C
C *** Copyright (c) 2002  I.S. Duff and J. Koster                   ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C
C PR(I) is the previous row to I in the depth first search.
C   It is used as a work array in the sorting algorithm.
C Elements (IPERM(I),I) I=1,...M  are non-zero at the end of the
C   algorithm unless N assignments have not been made.  In which case
C   (IPERM(I),I) will be zero for N-NUM entries.
C CV(I) is the most recent row extension at which column I was visited.
C ARP(I) is one less than the number of non-zeros in row I
C   which have not been scanned when looking for a cheap assignment.
C OUT(I) is one less than the number of non-zeros in row I
C   which have not been scanned during one pass through the main loop.
C
      INTEGER LIRN,M,N,NUM
      INTEGER ARP(N),CV(M),IRN(LIRN),IP(N),IPERM(M),LENC(N),OUT(N),PR(N)

      INTEGER I,II,IN1,IN2,J,J1,JORD,K,KK
C External subroutines and/or functions
      EXTERNAL MC64XD

      DO 10 I = 1,M
        CV(I) = 0
        IPERM(I) = 0
   10 CONTINUE
      DO 12 J = 1,N
        ARP(J) = LENC(J) - 1
   12 CONTINUE
      NUM = 0
C
C Main loop. Each pass round this loop either results in a new 
C assignment or gives a row with no assignment.
C
      DO 1000 JORD = 1,N
C
        J = JORD
        PR(J) = -1
        DO 70 K = 1,JORD
C Look for a cheap assignment
          IN1 = ARP(J)
          IF (IN1.LT.0) GO TO 30
          IN2 = IP(J) + LENC(J) - 1
          IN1 = IN2 - IN1
          DO 20 II = IN1,IN2
            I = IRN(II)
            IF (IPERM(I).EQ.0) GO TO 80
   20     CONTINUE
C No cheap assignment in row.
          ARP(J) = -1
C Begin looking for assignment chain starting with row J.
   30     CONTINUE
          OUT(J) = LENC(J) - 1
C Inner loop.  Extends chain by one or backtracks.
          DO 60 KK = 1,JORD
            IN1 = OUT(J)
            IF (IN1.LT.0) GO TO 50
            IN2 = IP(J) + LENC(J) - 1
            IN1 = IN2 - IN1
C Forward scan.
            DO 40 II = IN1,IN2
              I = IRN(II)
              IF (CV(I).EQ.JORD) GO TO 40
C Column I has not yet been accessed during this pass.
              J1 = J
              J = IPERM(I)
              CV(I) = JORD
              PR(J) = J1
              OUT(J1) = IN2 - II - 1
              GO TO 70
   40       CONTINUE
C Backtracking step.
   50       CONTINUE
            J = PR(J)
            IF (J.EQ.-1) GO TO 1000
   60     CONTINUE
   70   CONTINUE
C
C New assignment is made.
   80   CONTINUE
        IPERM(I) = J
        ARP(J) = IN2 - II - 1
        NUM = NUM + 1
        DO 90 K = 1,JORD
          J = PR(J)
          IF (J.EQ.-1) GO TO 1000
          II = IP(J) + LENC(J) - OUT(J) - 2
          I = IRN(II)
          IPERM(I) = J
   90   CONTINUE
C
 1000 CONTINUE


C IPERM is complete if M = N and NUM = N
      IF (M.EQ.N .and. NUM.EQ.N) GO TO 2000

C Complete IPERM; CV, ARP are work arrays
      CALL MC64XD(M,N,IPERM,CV,ARP)

 2000 RETURN
      END

C**********************************************************************
CCCCC LAST UPDATE Tue Nov 26 03:20:26 MET 2002
      SUBROUTINE MC64XD(M,N,IPERM,RW,CW)
      IMPLICIT NONE
C
C *** Copyright (c) 2002  I.S. Duff and J. Koster                   ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C
C Complete the (incomplete) row permutation in IPERM.
C
      INTEGER M,N
      INTEGER RW(M),CW(N),IPERM(M)

      INTEGER I,J,K

C If M=N, the matrix is structurally singular, complete IPERM
C If M>N, the matrix is rectangular, complete IPERM

C RW, CW are work arrays;
C Store indices of unmatched rows in RW
C Mark matched columns in CW

      DO 10 J = 1,N
        CW(J) = 0
   10 CONTINUE
      K = 0
      DO 20 I = 1,M
        IF (IPERM(I).EQ.0) THEN
          K = K + 1
          RW(K) = I
        ELSE
          J = IPERM(I)
          CW(J) = I
        ENDIF
   20 CONTINUE
      K = 0
      DO 30 J = 1,N
        IF (CW(J).NE.0) GO TO 30
        K = K + 1
        I = RW(K)
        IPERM(I) = -J
   30 CONTINUE
      DO 40 J = N+1,M
        K = K + 1
        I = RW(K)
        IPERM(I) = -J
   40 CONTINUE

      RETURN
      END

C**********************************************************************
CCCCC LAST UPDATE Tue Nov 26 03:20:26 MET 2002
      SUBROUTINE MC64YD(N,PERM,IW,INFO)
      IMPLICIT NONE
C
C *** Copyright (c) 2002  I.S. Duff and J. Koster                   ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C
C Check whether PERM(1:N) is a valid permutation vector.
C
      INTEGER N,PERM(N),IW(N),INFO 

      INTEGER I,J

      DO 20 J = 1,N
        IW(J) = 0
   20 CONTINUE
      DO 30 J = 1,N
        I = PERM(J)
        IF (I.LT.1 .OR. I.GT.N) GO TO 99
        IF (IW(I).NE.0) GO TO 99
        IW(I) = 1
   30 CONTINUE
C Successful return
      INFO = 0
      RETURN
C Error return
   99 INFO = -1
      RETURN
      END


