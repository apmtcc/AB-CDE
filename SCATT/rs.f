      DOUBLE PRECISION FUNCTION EPSLON (X)                              EPS00010
      DOUBLE PRECISION X                                                EPS00020
C                                                                       EPS00030
C     ESTIMATE UNIT ROUNDOFF IN QUANTITIES OF SIZE X.                   EPS00040
C                                                                       EPS00050
      DOUBLE PRECISION A,B,C,EPS                                        EPS00060
C                                                                       EPS00070
C     THIS PROGRAM SHOULD FUNCTION PROPERLY ON ALL SYSTEMS              EPS00080
C     SATISFYING THE FOLLOWING TWO ASSUMPTIONS,                         EPS00090
C        1.  THE BASE USED IN REPRESENTING FLOATING POINT               EPS00100
C            NUMBERS IS NOT A POWER OF THREE.                           EPS00110
C        2.  THE QUANTITY  A  IN STATEMENT 10 IS REPRESENTED TO         EPS00120
C            THE ACCURACY USED IN FLOATING POINT VARIABLES              EPS00130
C            THAT ARE STORED IN MEMORY.                                 EPS00140
C     THE STATEMENT NUMBER 10 AND THE GO TO 10 ARE INTENDED TO          EPS00150
C     FORCE OPTIMIZING COMPILERS TO GENERATE CODE SATISFYING            EPS00160
C     ASSUMPTION 2.                                                     EPS00170
C     UNDER THESE ASSUMPTIONS, IT SHOULD BE TRUE THAT,                  EPS00180
C            A  IS NOT EXACTLY EQUAL TO FOUR-THIRDS,                    EPS00190
C            B  HAS A ZERO FOR ITS LAST BIT OR DIGIT,                   EPS00200
C            C  IS NOT EXACTLY EQUAL TO ONE,                            EPS00210
C            EPS  MEASURES THE SEPARATION OF 1.0 FROM                   EPS00220
C                 THE NEXT LARGER FLOATING POINT NUMBER.                EPS00230
C     THE DEVELOPERS OF EISPACK WOULD APPRECIATE BEING INFORMED         EPS00240
C     ABOUT ANY SYSTEMS WHERE THESE ASSUMPTIONS DO NOT HOLD.            EPS00250
C                                                                       EPS00260
C     THIS VERSION DATED 4/6/83.                                        EPS00270
C                                                                       EPS00280
      A = 4.0D0/3.0D0                                                   EPS00290
   10 B = A - 1.0D0                                                     EPS00300
      C = B + B + B                                                     EPS00310
      EPS = DABS(C-1.0D0)                                               EPS00320
      IF (EPS .EQ. 0.0D0) GO TO 10                                      EPS00330
      EPSLON = EPS*DABS(X)                                              EPS00340
      RETURN                                                            EPS00350
      END                                                               EPS00360
      DOUBLE PRECISION FUNCTION PYTHAG(A,B)                             PYT00010
      DOUBLE PRECISION A,B                                              PYT00020
C                                                                       PYT00030
C     FINDS DSQRT(A**2+B**2) WITHOUT OVERFLOW OR DESTRUCTIVE UNDERFLOW  PYT00040
C                                                                       PYT00050
      DOUBLE PRECISION P,R,S,T,U                                        PYT00060
      P = DMAX1(DABS(A),DABS(B))                                        PYT00070
      IF (P .EQ. 0.0D0) GO TO 20                                        PYT00080
      R = (DMIN1(DABS(A),DABS(B))/P)**2                                 PYT00090
   10 CONTINUE                                                          PYT00100
         T = 4.0D0 + R                                                  PYT00110
         IF (T .EQ. 4.0D0) GO TO 20                                     PYT00120
         S = R/T                                                        PYT00130
         U = 1.0D0 + 2.0D0*S                                            PYT00140
         P = U*P                                                        PYT00150
         R = (S/U)**2 * R                                               PYT00160
      GO TO 10                                                          PYT00170
   20 PYTHAG = P                                                        PYT00180
      RETURN                                                            PYT00190
      END                                                               PYT00200
      SUBROUTINE RS(NM,N,A,W,MATZ,Z,FV1,FV2,IERR)                       RS 00010
C                                                                       RS 00020
      INTEGER N,NM,IERR,MATZ                                            RS 00030
      DOUBLE PRECISION A(NM,N),W(N),Z(NM,N),FV1(N),FV2(N)               RS 00040
C                                                                       RS 00050
C     THIS SUBROUTINE CALLS THE RECOMMENDED SEQUENCE OF                 RS 00060
C     SUBROUTINES FROM THE EIGENSYSTEM SUBROUTINE PACKAGE (EISPACK)     RS 00070
C     TO FIND THE EIGENVALUES AND EIGENVECTORS (IF DESIRED)             RS 00080
C     OF A REAL SYMMETRIC MATRIX.                                       RS 00090
C                                                                       RS 00100
C     ON INPUT                                                          RS 00110
C                                                                       RS 00120
C        NM  MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL    RS 00130
C        ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM            RS 00140
C        DIMENSION STATEMENT.                                           RS 00150
C                                                                       RS 00160
C        N  IS THE ORDER OF THE MATRIX  A.                              RS 00170
C                                                                       RS 00180
C        A  CONTAINS THE REAL SYMMETRIC MATRIX.                         RS 00190
C                                                                       RS 00200
C        MATZ  IS AN INTEGER VARIABLE SET EQUAL TO ZERO IF              RS 00210
C        ONLY EIGENVALUES ARE DESIRED.  OTHERWISE IT IS SET TO          RS 00220
C        ANY NON-ZERO INTEGER FOR BOTH EIGENVALUES AND EIGENVECTORS.    RS 00230
C                                                                       RS 00240
C     ON OUTPUT                                                         RS 00250
C                                                                       RS 00260
C        W  CONTAINS THE EIGENVALUES IN ASCENDING ORDER.                RS 00270
C                                                                       RS 00280
C        Z  CONTAINS THE EIGENVECTORS IF MATZ IS NOT ZERO.              RS 00290
C                                                                       RS 00300
C        IERR  IS AN INTEGER OUTPUT VARIABLE SET EQUAL TO AN ERROR      RS 00310
C           COMPLETION CODE DESCRIBED IN THE DOCUMENTATION FOR TQLRAT   RS 00320
C           AND TQL2.  THE NORMAL COMPLETION CODE IS ZERO.              RS 00330
C                                                                       RS 00340
C        FV1  AND  FV2  ARE TEMPORARY STORAGE ARRAYS.                   RS 00350
C                                                                       RS 00360
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,    RS 00370
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY RS 00380
C                                                                       RS 00390
C     THIS VERSION DATED AUGUST 1983.                                   RS 00400
C                                                                       RS 00410
C     ------------------------------------------------------------------RS 00420
C                                                                       RS 00430
      IF (N .LE. NM) GO TO 10                                           RS 00440
      IERR = 10 * N                                                     RS 00450
      GO TO 50                                                          RS 00460
C                                                                       RS 00470
   10 IF (MATZ .NE. 0) GO TO 20                                         RS 00480
C     .......... FIND EIGENVALUES ONLY ..........                       RS 00490
      CALL  TRED1(NM,N,A,W,FV1,FV2)                                     RS 00500
      CALL  TQLRAT(N,W,FV2,IERR)                                        RS 00510
      GO TO 50                                                          RS 00520
C     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........      RS 00530
   20 CALL  TRED2(NM,N,A,W,FV1,Z)                                       RS 00540
      CALL  TQL2(NM,N,W,FV1,Z,IERR)                                     RS 00550
   50 RETURN                                                            RS 00560
      END                                                               RS 00570
      SUBROUTINE TQL2(NM,N,D,E,Z,IERR)                                  TQL00010
C                                                                       TQL00020
      INTEGER I,J,K,L,M,N,II,L1,L2,NM,MML,IERR                          TQL00030
      DOUBLE PRECISION D(N),E(N),Z(NM,N)                                TQL00040
      DOUBLE PRECISION C,C2,C3,DL1,EL1,F,G,H,P,R,S,S2,TST1,TST2,PYTHAG  TQL00050
C                                                                       TQL00060
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQL2,     TQL00070
C     NUM. MATH. 11, 293-306(1968) BY BOWDLER, MARTIN, REINSCH, AND     TQL00080
C     WILKINSON.                                                        TQL00090
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).   TQL00100
C                                                                       TQL00110
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS            TQL00120
C     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE QL METHOD.               TQL00130
C     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO              TQL00140
C     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS                  TQL00150
C     FULL MATRIX TO TRIDIAGONAL FORM.                                  TQL00160
C                                                                       TQL00170
C     ON INPUT                                                          TQL00180
C                                                                       TQL00190
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL         TQL00200
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM          TQL00210
C          DIMENSION STATEMENT.                                         TQL00220
C                                                                       TQL00230
C        N IS THE ORDER OF THE MATRIX.                                  TQL00240
C                                                                       TQL00250
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.          TQL00260
C                                                                       TQL00270
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX        TQL00280
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.               TQL00290
C                                                                       TQL00300
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE           TQL00310
C          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS      TQL00320
C          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN        TQL00330
C          THE IDENTITY MATRIX.                                         TQL00340
C                                                                       TQL00350
C      ON OUTPUT                                                        TQL00360
C                                                                       TQL00370
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN          TQL00380
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT          TQL00390
C          UNORDERED FOR INDICES 1,2,...,IERR-1.                        TQL00400
C                                                                       TQL00410
C        E HAS BEEN DESTROYED.                                          TQL00420
C                                                                       TQL00430
C        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC           TQL00440
C          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,     TQL00450
C          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED       TQL00460
C          EIGENVALUES.                                                 TQL00470
C                                                                       TQL00480
C        IERR IS SET TO                                                 TQL00490
C          ZERO       FOR NORMAL RETURN,                                TQL00500
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN               TQL00510
C                     DETERMINED AFTER 30 ITERATIONS.                   TQL00520
C                                                                       TQL00530
C     CALLS PYTHAG FOR  DSQRT(A*A + B*B) .                              TQL00540
C                                                                       TQL00550
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,    TQL00560
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY TQL00570
C                                                                       TQL00580
C     THIS VERSION DATED AUGUST 1983.                                   TQL00590
C                                                                       TQL00600
C     ------------------------------------------------------------------TQL00610
C                                                                       TQL00620
      IERR = 0                                                          TQL00630
      IF (N .EQ. 1) GO TO 1001                                          TQL00640
C                                                                       TQL00650
      DO 100 I = 2, N                                                   TQL00660
  100 E(I-1) = E(I)                                                     TQL00670
C                                                                       TQL00680
      F = 0.0D0                                                         TQL00690
      TST1 = 0.0D0                                                      TQL00700
      E(N) = 0.0D0                                                      TQL00710
C                                                                       TQL00720
      DO 240 L = 1, N                                                   TQL00730
         J = 0                                                          TQL00740
         H = DABS(D(L)) + DABS(E(L))                                    TQL00750
         IF (TST1 .LT. H) TST1 = H                                      TQL00760
C     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........         TQL00770
         DO 110 M = L, N                                                TQL00780
            TST2 = TST1 + DABS(E(M))                                    TQL00790
            IF (TST2 .EQ. TST1) GO TO 120                               TQL00800
C     .......... E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT               TQL00810
C                THROUGH THE BOTTOM OF THE LOOP ..........              TQL00820
  110    CONTINUE                                                       TQL00830
C                                                                       TQL00840
  120    IF (M .EQ. L) GO TO 220                                        TQL00850
  130    IF (J .EQ. 30) GO TO 1000                                      TQL00860
         J = J + 1                                                      TQL00870
C     .......... FORM SHIFT ..........                                  TQL00880
         L1 = L + 1                                                     TQL00890
         L2 = L1 + 1                                                    TQL00900
         G = D(L)                                                       TQL00910
         P = (D(L1) - G) / (2.0D0 * E(L))                               TQL00920
         R = PYTHAG(P,1.0D0)                                            TQL00930
         D(L) = E(L) / (P + DSIGN(R,P))                                 TQL00940
         D(L1) = E(L) * (P + DSIGN(R,P))                                TQL00950
         DL1 = D(L1)                                                    TQL00960
         H = G - D(L)                                                   TQL00970
         IF (L2 .GT. N) GO TO 145                                       TQL00980
C                                                                       TQL00990
         DO 140 I = L2, N                                               TQL01000
  140    D(I) = D(I) - H                                                TQL01010
C                                                                       TQL01020
  145    F = F + H                                                      TQL01030
C     .......... QL TRANSFORMATION ..........                           TQL01040
         P = D(M)                                                       TQL01050
         C = 1.0D0                                                      TQL01060
         C2 = C                                                         TQL01070
         EL1 = E(L1)                                                    TQL01080
         S = 0.0D0                                                      TQL01090
         MML = M - L                                                    TQL01100
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........             TQL01110
         DO 200 II = 1, MML                                             TQL01120
            C3 = C2                                                     TQL01130
            C2 = C                                                      TQL01140
            S2 = S                                                      TQL01150
            I = M - II                                                  TQL01160
            G = C * E(I)                                                TQL01170
            H = C * P                                                   TQL01180
            R = PYTHAG(P,E(I))                                          TQL01190
            E(I+1) = S * R                                              TQL01200
            S = E(I) / R                                                TQL01210
            C = P / R                                                   TQL01220
            P = C * D(I) - S * G                                        TQL01230
            D(I+1) = H + S * (C * G + S * D(I))                         TQL01240
C     .......... FORM VECTOR ..........                                 TQL01250
            DO 180 K = 1, N                                             TQL01260
               H = Z(K,I+1)                                             TQL01270
               Z(K,I+1) = S * Z(K,I) + C * H                            TQL01280
               Z(K,I) = C * Z(K,I) - S * H                              TQL01290
  180       CONTINUE                                                    TQL01300
C                                                                       TQL01310
  200    CONTINUE                                                       TQL01320
C                                                                       TQL01330
         P = -S * S2 * C3 * EL1 * E(L) / DL1                            TQL01340
         E(L) = S * P                                                   TQL01350
         D(L) = C * P                                                   TQL01360
         TST2 = TST1 + DABS(E(L))                                       TQL01370
         IF (TST2 .GT. TST1) GO TO 130                                  TQL01380
  220    D(L) = D(L) + F                                                TQL01390
  240 CONTINUE                                                          TQL01400
C     .......... ORDER EIGENVALUES AND EIGENVECTORS ..........          TQL01410
      DO 300 II = 2, N                                                  TQL01420
         I = II - 1                                                     TQL01430
         K = I                                                          TQL01440
         P = D(I)                                                       TQL01450
C                                                                       TQL01460
         DO 260 J = II, N                                               TQL01470
            IF (D(J) .GE. P) GO TO 260                                  TQL01480
            K = J                                                       TQL01490
            P = D(J)                                                    TQL01500
  260    CONTINUE                                                       TQL01510
C                                                                       TQL01520
         IF (K .EQ. I) GO TO 300                                        TQL01530
         D(K) = D(I)                                                    TQL01540
         D(I) = P                                                       TQL01550
C                                                                       TQL01560
         DO 280 J = 1, N                                                TQL01570
            P = Z(J,I)                                                  TQL01580
            Z(J,I) = Z(J,K)                                             TQL01590
            Z(J,K) = P                                                  TQL01600
  280    CONTINUE                                                       TQL01610
C                                                                       TQL01620
  300 CONTINUE                                                          TQL01630
C                                                                       TQL01640
      GO TO 1001                                                        TQL01650
C     .......... SET ERROR -- NO CONVERGENCE TO AN                      TQL01660
C                EIGENVALUE AFTER 30 ITERATIONS ..........              TQL01670
 1000 IERR = L                                                          TQL01680
 1001 RETURN                                                            TQL01690
      END                                                               TQL01700
      SUBROUTINE TQLRAT(N,D,E2,IERR)                                    TQL00010
C                                                                       TQL00020
      INTEGER I,J,L,M,N,II,L1,MML,IERR                                  TQL00030
      DOUBLE PRECISION D(N),E2(N)                                       TQL00040
      DOUBLE PRECISION B,C,F,G,H,P,R,S,T,EPSLON,PYTHAG                  TQL00050
C                                                                       TQL00060
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQLRAT,   TQL00070
C     ALGORITHM 464, COMM. ACM 16, 689(1973) BY REINSCH.                TQL00080
C                                                                       TQL00090
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC              TQL00100
C     TRIDIAGONAL MATRIX BY THE RATIONAL QL METHOD.                     TQL00110
C                                                                       TQL00120
C     ON INPUT                                                          TQL00130
C                                                                       TQL00140
C        N IS THE ORDER OF THE MATRIX.                                  TQL00150
C                                                                       TQL00160
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.          TQL00170
C                                                                       TQL00180
C        E2 CONTAINS THE SQUARES OF THE SUBDIAGONAL ELEMENTS OF THE     TQL00190
C          INPUT MATRIX IN ITS LAST N-1 POSITIONS.  E2(1) IS ARBITRARY. TQL00200
C                                                                       TQL00210
C      ON OUTPUT                                                        TQL00220
C                                                                       TQL00230
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN          TQL00240
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND          TQL00250
C          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE            TQL00260
C          THE SMALLEST EIGENVALUES.                                    TQL00270
C                                                                       TQL00280
C        E2 HAS BEEN DESTROYED.                                         TQL00290
C                                                                       TQL00300
C        IERR IS SET TO                                                 TQL00310
C          ZERO       FOR NORMAL RETURN,                                TQL00320
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN               TQL00330
C                     DETERMINED AFTER 30 ITERATIONS.                   TQL00340
C                                                                       TQL00350
C     CALLS PYTHAG FOR  DSQRT(A*A + B*B) .                              TQL00360
C                                                                       TQL00370
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,    TQL00380
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY TQL00390
C                                                                       TQL00400
C     THIS VERSION DATED AUGUST 1983.                                   TQL00410
C                                                                       TQL00420
C     ------------------------------------------------------------------TQL00430
C                                                                       TQL00440
      IERR = 0                                                          TQL00450
      IF (N .EQ. 1) GO TO 1001                                          TQL00460
C                                                                       TQL00470
      DO 100 I = 2, N                                                   TQL00480
  100 E2(I-1) = E2(I)                                                   TQL00490
C                                                                       TQL00500
      F = 0.0D0                                                         TQL00510
      T = 0.0D0                                                         TQL00520
      E2(N) = 0.0D0                                                     TQL00530
C                                                                       TQL00540
      DO 290 L = 1, N                                                   TQL00550
         J = 0                                                          TQL00560
         H = DABS(D(L)) + DSQRT(E2(L))                                  TQL00570
         IF (T .GT. H) GO TO 105                                        TQL00580
         T = H                                                          TQL00590
         B = EPSLON(T)                                                  TQL00600
         C = B * B                                                      TQL00610
C     .......... LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT .......... TQL00620
  105    DO 110 M = L, N                                                TQL00630
            IF (E2(M) .LE. C) GO TO 120                                 TQL00640
C     .......... E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT              TQL00650
C                THROUGH THE BOTTOM OF THE LOOP ..........              TQL00660
  110    CONTINUE                                                       TQL00670
C                                                                       TQL00680
  120    IF (M .EQ. L) GO TO 210                                        TQL00690
  130    IF (J .EQ. 30) GO TO 1000                                      TQL00700
         J = J + 1                                                      TQL00710
C     .......... FORM SHIFT ..........                                  TQL00720
         L1 = L + 1                                                     TQL00730
         S = DSQRT(E2(L))                                               TQL00740
         G = D(L)                                                       TQL00750
         P = (D(L1) - G) / (2.0D0 * S)                                  TQL00760
         R = PYTHAG(P,1.0D0)                                            TQL00770
         D(L) = S / (P + DSIGN(R,P))                                    TQL00780
         H = G - D(L)                                                   TQL00790
C                                                                       TQL00800
         DO 140 I = L1, N                                               TQL00810
  140    D(I) = D(I) - H                                                TQL00820
C                                                                       TQL00830
         F = F + H                                                      TQL00840
C     .......... RATIONAL QL TRANSFORMATION ..........                  TQL00850
         G = D(M)                                                       TQL00860
         IF (G .EQ. 0.0D0) G = B                                        TQL00870
         H = G                                                          TQL00880
         S = 0.0D0                                                      TQL00890
         MML = M - L                                                    TQL00900
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........             TQL00910
         DO 200 II = 1, MML                                             TQL00920
            I = M - II                                                  TQL00930
            P = G * H                                                   TQL00940
            R = P + E2(I)                                               TQL00950
            E2(I+1) = S * R                                             TQL00960
            S = E2(I) / R                                               TQL00970
            D(I+1) = H + S * (H + D(I))                                 TQL00980
            G = D(I) - E2(I) / G                                        TQL00990
            IF (G .EQ. 0.0D0) G = B                                     TQL01000
            H = G * P / R                                               TQL01010
  200    CONTINUE                                                       TQL01020
C                                                                       TQL01030
         E2(L) = S * G                                                  TQL01040
         D(L) = H                                                       TQL01050
C     .......... GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST .......... TQL01060
         IF (H .EQ. 0.0D0) GO TO 210                                    TQL01070
         IF (DABS(E2(L)) .LE. DABS(C/H)) GO TO 210                      TQL01080
         E2(L) = H * E2(L)                                              TQL01090
         IF (E2(L) .NE. 0.0D0) GO TO 130                                TQL01100
  210    P = D(L) + F                                                   TQL01110
C     .......... ORDER EIGENVALUES ..........                           TQL01120
         IF (L .EQ. 1) GO TO 250                                        TQL01130
C     .......... FOR I=L STEP -1 UNTIL 2 DO -- ..........               TQL01140
         DO 230 II = 2, L                                               TQL01150
            I = L + 2 - II                                              TQL01160
            IF (P .GE. D(I-1)) GO TO 270                                TQL01170
            D(I) = D(I-1)                                               TQL01180
  230    CONTINUE                                                       TQL01190
C                                                                       TQL01200
  250    I = 1                                                          TQL01210
  270    D(I) = P                                                       TQL01220
  290 CONTINUE                                                          TQL01230
C                                                                       TQL01240
      GO TO 1001                                                        TQL01250
C     .......... SET ERROR -- NO CONVERGENCE TO AN                      TQL01260
C                EIGENVALUE AFTER 30 ITERATIONS ..........              TQL01270
 1000 IERR = L                                                          TQL01280
 1001 RETURN                                                            TQL01290
      END                                                               TQL01300
      SUBROUTINE TRED1(NM,N,A,D,E,E2)                                   TRE00010
C                                                                       TRE00020
      INTEGER I,J,K,L,N,II,NM,JP1                                       TRE00030
      DOUBLE PRECISION A(NM,N),D(N),E(N),E2(N)                          TRE00040
      DOUBLE PRECISION F,G,H,SCALE                                      TRE00050
C                                                                       TRE00060
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED1,    TRE00070
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.   TRE00080
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).   TRE00090
C                                                                       TRE00100
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX                   TRE00110
C     TO A SYMMETRIC TRIDIAGONAL MATRIX USING                           TRE00120
C     ORTHOGONAL SIMILARITY TRANSFORMATIONS.                            TRE00130
C                                                                       TRE00140
C     ON INPUT                                                          TRE00150
C                                                                       TRE00160
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL         TRE00170
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM          TRE00180
C          DIMENSION STATEMENT.                                         TRE00190
C                                                                       TRE00200
C        N IS THE ORDER OF THE MATRIX.                                  TRE00210
C                                                                       TRE00220
C        A CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE          TRE00230
C          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.               TRE00240
C                                                                       TRE00250
C     ON OUTPUT                                                         TRE00260
C                                                                       TRE00270
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANS-             TRE00280
C          FORMATIONS USED IN THE REDUCTION IN ITS STRICT LOWER         TRE00290
C          TRIANGLE.  THE FULL UPPER TRIANGLE OF A IS UNALTERED.        TRE00300
C                                                                       TRE00310
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX.    TRE00320
C                                                                       TRE00330
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL         TRE00340
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO.      TRE00350
C                                                                       TRE00360
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.    TRE00370
C          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.        TRE00380
C                                                                       TRE00390
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,    TRE00400
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY TRE00410
C                                                                       TRE00420
C     THIS VERSION DATED AUGUST 1983.                                   TRE00430
C                                                                       TRE00440
C     ------------------------------------------------------------------TRE00450
C                                                                       TRE00460
      DO 100 I = 1, N                                                   TRE00470
         D(I) = A(N,I)                                                  TRE00480
         A(N,I) = A(I,I)                                                TRE00490
  100 CONTINUE                                                          TRE00500
C     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........               TRE00510
      DO 300 II = 1, N                                                  TRE00520
         I = N + 1 - II                                                 TRE00530
         L = I - 1                                                      TRE00540
         H = 0.0D0                                                      TRE00550
         SCALE = 0.0D0                                                  TRE00560
         IF (L .LT. 1) GO TO 130                                        TRE00570
C     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........       TRE00580
         DO 120 K = 1, L                                                TRE00590
  120    SCALE = SCALE + DABS(D(K))                                     TRE00600
C                                                                       TRE00610
         IF (SCALE .NE. 0.0D0) GO TO 140                                TRE00620
C                                                                       TRE00630
         DO 125 J = 1, L                                                TRE00640
            D(J) = A(L,J)                                               TRE00650
            A(L,J) = A(I,J)                                             TRE00660
            A(I,J) = 0.0D0                                              TRE00670
  125    CONTINUE                                                       TRE00680
C                                                                       TRE00690
  130    E(I) = 0.0D0                                                   TRE00700
         E2(I) = 0.0D0                                                  TRE00710
         GO TO 300                                                      TRE00720
C                                                                       TRE00730
  140    DO 150 K = 1, L                                                TRE00740
            D(K) = D(K) / SCALE                                         TRE00750
            H = H + D(K) * D(K)                                         TRE00760
  150    CONTINUE                                                       TRE00770
C                                                                       TRE00780
         E2(I) = SCALE * SCALE * H                                      TRE00790
         F = D(L)                                                       TRE00800
         G = -DSIGN(DSQRT(H),F)                                         TRE00810
         E(I) = SCALE * G                                               TRE00820
         H = H - F * G                                                  TRE00830
         D(L) = F - G                                                   TRE00840
         IF (L .EQ. 1) GO TO 285                                        TRE00850
C     .......... FORM A*U ..........                                    TRE00860
         DO 170 J = 1, L                                                TRE00870
  170    E(J) = 0.0D0                                                   TRE00880
C                                                                       TRE00890
         DO 240 J = 1, L                                                TRE00900
            F = D(J)                                                    TRE00910
            G = E(J) + A(J,J) * F                                       TRE00920
            JP1 = J + 1                                                 TRE00930
            IF (L .LT. JP1) GO TO 220                                   TRE00940
C                                                                       TRE00950
            DO 200 K = JP1, L                                           TRE00960
               G = G + A(K,J) * D(K)                                    TRE00970
               E(K) = E(K) + A(K,J) * F                                 TRE00980
  200       CONTINUE                                                    TRE00990
C                                                                       TRE01000
  220       E(J) = G                                                    TRE01010
  240    CONTINUE                                                       TRE01020
C     .......... FORM P ..........                                      TRE01030
         F = 0.0D0                                                      TRE01040
C                                                                       TRE01050
         DO 245 J = 1, L                                                TRE01060
            E(J) = E(J) / H                                             TRE01070
            F = F + E(J) * D(J)                                         TRE01080
  245    CONTINUE                                                       TRE01090
C                                                                       TRE01100
         H = F / (H + H)                                                TRE01110
C     .......... FORM Q ..........                                      TRE01120
         DO 250 J = 1, L                                                TRE01130
  250    E(J) = E(J) - H * D(J)                                         TRE01140
C     .......... FORM REDUCED A ..........                              TRE01150
         DO 280 J = 1, L                                                TRE01160
            F = D(J)                                                    TRE01170
            G = E(J)                                                    TRE01180
C                                                                       TRE01190
            DO 260 K = J, L                                             TRE01200
  260       A(K,J) = A(K,J) - F * E(K) - G * D(K)                       TRE01210
C                                                                       TRE01220
  280    CONTINUE                                                       TRE01230
C                                                                       TRE01240
  285    DO 290 J = 1, L                                                TRE01250
            F = D(J)                                                    TRE01260
            D(J) = A(L,J)                                               TRE01270
            A(L,J) = A(I,J)                                             TRE01280
            A(I,J) = F * SCALE                                          TRE01290
  290    CONTINUE                                                       TRE01300
C                                                                       TRE01310
  300 CONTINUE                                                          TRE01320
C                                                                       TRE01330
      RETURN                                                            TRE01340
      END                                                               TRE01350
      SUBROUTINE TRED2(NM,N,A,D,E,Z)                                    TRE00010
C                                                                       TRE00020
      INTEGER I,J,K,L,N,II,NM,JP1                                       TRE00030
      DOUBLE PRECISION A(NM,N),D(N),E(N),Z(NM,N)                        TRE00040
      DOUBLE PRECISION F,G,H,HH,SCALE                                   TRE00050
C                                                                       TRE00060
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED2,    TRE00070
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.   TRE00080
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).   TRE00090
C                                                                       TRE00100
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX TO A              TRE00110
C     SYMMETRIC TRIDIAGONAL MATRIX USING AND ACCUMULATING               TRE00120
C     ORTHOGONAL SIMILARITY TRANSFORMATIONS.                            TRE00130
C                                                                       TRE00140
C     ON INPUT                                                          TRE00150
C                                                                       TRE00160
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL         TRE00170
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM          TRE00180
C          DIMENSION STATEMENT.                                         TRE00190
C                                                                       TRE00200
C        N IS THE ORDER OF THE MATRIX.                                  TRE00210
C                                                                       TRE00220
C        A CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE          TRE00230
C          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.               TRE00240
C                                                                       TRE00250
C     ON OUTPUT                                                         TRE00260
C                                                                       TRE00270
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX.    TRE00280
C                                                                       TRE00290
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL         TRE00300
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO.      TRE00310
C                                                                       TRE00320
C        Z CONTAINS THE ORTHOGONAL TRANSFORMATION MATRIX                TRE00330
C          PRODUCED IN THE REDUCTION.                                   TRE00340
C                                                                       TRE00350
C        A AND Z MAY COINCIDE.  IF DISTINCT, A IS UNALTERED.            TRE00360
C                                                                       TRE00370
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,    TRE00380
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY TRE00390
C                                                                       TRE00400
C     THIS VERSION DATED AUGUST 1983.                                   TRE00410
C                                                                       TRE00420
C     ------------------------------------------------------------------TRE00430
C                                                                       TRE00440
      DO 100 I = 1, N                                                   TRE00450
C                                                                       TRE00460
         DO 80 J = I, N                                                 TRE00470
   80    Z(J,I) = A(J,I)                                                TRE00480
C                                                                       TRE00490
         D(I) = A(N,I)                                                  TRE00500
  100 CONTINUE                                                          TRE00510
C                                                                       TRE00520
      IF (N .EQ. 1) GO TO 510                                           TRE00530
C     .......... FOR I=N STEP -1 UNTIL 2 DO -- ..........               TRE00540
      DO 300 II = 2, N                                                  TRE00550
         I = N + 2 - II                                                 TRE00560
         L = I - 1                                                      TRE00570
         H = 0.0D0                                                      TRE00580
         SCALE = 0.0D0                                                  TRE00590
         IF (L .LT. 2) GO TO 130                                        TRE00600
C     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........       TRE00610
         DO 120 K = 1, L                                                TRE00620
  120    SCALE = SCALE + DABS(D(K))                                     TRE00630
C                                                                       TRE00640
         IF (SCALE .NE. 0.0D0) GO TO 140                                TRE00650
  130    E(I) = D(L)                                                    TRE00660
C                                                                       TRE00670
         DO 135 J = 1, L                                                TRE00680
            D(J) = Z(L,J)                                               TRE00690
            Z(I,J) = 0.0D0                                              TRE00700
            Z(J,I) = 0.0D0                                              TRE00710
  135    CONTINUE                                                       TRE00720
C                                                                       TRE00730
         GO TO 290                                                      TRE00740
C                                                                       TRE00750
  140    DO 150 K = 1, L                                                TRE00760
            D(K) = D(K) / SCALE                                         TRE00770
            H = H + D(K) * D(K)                                         TRE00780
  150    CONTINUE                                                       TRE00790
C                                                                       TRE00800
         F = D(L)                                                       TRE00810
         G = -DSIGN(DSQRT(H),F)                                         TRE00820
         E(I) = SCALE * G                                               TRE00830
         H = H - F * G                                                  TRE00840
         D(L) = F - G                                                   TRE00850
C     .......... FORM A*U ..........                                    TRE00860
         DO 170 J = 1, L                                                TRE00870
  170    E(J) = 0.0D0                                                   TRE00880
C                                                                       TRE00890
         DO 240 J = 1, L                                                TRE00900
            F = D(J)                                                    TRE00910
            Z(J,I) = F                                                  TRE00920
            G = E(J) + Z(J,J) * F                                       TRE00930
            JP1 = J + 1                                                 TRE00940
            IF (L .LT. JP1) GO TO 220                                   TRE00950
C                                                                       TRE00960
            DO 200 K = JP1, L                                           TRE00970
               G = G + Z(K,J) * D(K)                                    TRE00980
               E(K) = E(K) + Z(K,J) * F                                 TRE00990
  200       CONTINUE                                                    TRE01000
C                                                                       TRE01010
  220       E(J) = G                                                    TRE01020
  240    CONTINUE                                                       TRE01030
C     .......... FORM P ..........                                      TRE01040
         F = 0.0D0                                                      TRE01050
C                                                                       TRE01060
         DO 245 J = 1, L                                                TRE01070
            E(J) = E(J) / H                                             TRE01080
            F = F + E(J) * D(J)                                         TRE01090
  245    CONTINUE                                                       TRE01100
C                                                                       TRE01110
         HH = F / (H + H)                                               TRE01120
C     .......... FORM Q ..........                                      TRE01130
         DO 250 J = 1, L                                                TRE01140
  250    E(J) = E(J) - HH * D(J)                                        TRE01150
C     .......... FORM REDUCED A ..........                              TRE01160
         DO 280 J = 1, L                                                TRE01170
            F = D(J)                                                    TRE01180
            G = E(J)                                                    TRE01190
C                                                                       TRE01200
            DO 260 K = J, L                                             TRE01210
  260       Z(K,J) = Z(K,J) - F * E(K) - G * D(K)                       TRE01220
C                                                                       TRE01230
            D(J) = Z(L,J)                                               TRE01240
            Z(I,J) = 0.0D0                                              TRE01250
  280    CONTINUE                                                       TRE01260
C                                                                       TRE01270
  290    D(I) = H                                                       TRE01280
  300 CONTINUE                                                          TRE01290
C     .......... ACCUMULATION OF TRANSFORMATION MATRICES ..........     TRE01300
      DO 500 I = 2, N                                                   TRE01310
         L = I - 1                                                      TRE01320
         Z(N,L) = Z(L,L)                                                TRE01330
         Z(L,L) = 1.0D0                                                 TRE01340
         H = D(I)                                                       TRE01350
         IF (H .EQ. 0.0D0) GO TO 380                                    TRE01360
C                                                                       TRE01370
         DO 330 K = 1, L                                                TRE01380
  330    D(K) = Z(K,I) / H                                              TRE01390
C                                                                       TRE01400
         DO 360 J = 1, L                                                TRE01410
            G = 0.0D0                                                   TRE01420
C                                                                       TRE01430
            DO 340 K = 1, L                                             TRE01440
  340       G = G + Z(K,I) * Z(K,J)                                     TRE01450
C                                                                       TRE01460
            DO 360 K = 1, L                                             TRE01470
               Z(K,J) = Z(K,J) - G * D(K)                               TRE01480
  360    CONTINUE                                                       TRE01490
C                                                                       TRE01500
  380    DO 400 K = 1, L                                                TRE01510
  400    Z(K,I) = 0.0D0                                                 TRE01520
C                                                                       TRE01530
  500 CONTINUE                                                          TRE01540
C                                                                       TRE01550
  510 DO 520 I = 1, N                                                   TRE01560
         D(I) = Z(N,I)                                                  TRE01570
         Z(N,I) = 0.0D0                                                 TRE01580
  520 CONTINUE                                                          TRE01590
C                                                                       TRE01600
      Z(N,N) = 1.0D0                                                    TRE01610
      E(1) = 0.0D0                                                      TRE01620
      RETURN                                                            TRE01630
      END                                                               TRE01640

