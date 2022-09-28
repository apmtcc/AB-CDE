!
!  ***************************************************************************
!  *                                                                         *
!  *   Program    :  rotmat2                                                 *
!  *   Function   :  reduced rotation matrix d[j,kz1,kz2](theta)             *
!  *                                                                         *
!  ***************************************************************************
       subroutine ROTMAT2(maxnod,nodea1,ja1max,kz1min,kz1max,
     &                    kz2min,kz2max,xnoda1,wnoda1,rotmata1)
       implicit none
       integer :: maxnod,nodea1,ia1,ja1,kz1,kz2,
     &            ja1max,kz1min,kz1max,kz2min,kz2max

       real*8  :: xnoda1(nodea1),wnoda1(nodea1),beta,DM,
     &  rotmata1(maxnod,kz2min:kz2max,0:ja1max,kz1min:kz1max)
c       rotmata1(maxnod,kzmin:kzmax,0:ja1max,mzmin:mzmax)
c       where kzmin,kzmax : body-fixed z-axis
c             mzmin,mzmax : space-fixed z-axis
c
       do 3000 ja1=0,ja1max
       do 3000 kz1=kz1min,kz1max
       do 3000 kz2=kz2min,kz2max
       do 3000 ia1=1,nodea1
          beta=xnoda1(ia1)
          rotmata1(ia1,kz2,ja1,kz1)=DM(ja1,kz2,kz1,beta)
     &                        *dsqrt(ja1+0.5d0)*dsqrt(wnoda1(ia1)) 
 3000  continue 
       return
       end 
c==========================================================
      function dm(ij,ik,im,beta)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     This function uses eq. (4.1.23) of Edmonds 
c     to calculate the reduced rotation matrix element 
c     d(j,k,m;beta) = <jk|exp(+i*beta*Jy/hbar)|jm>
c     -----------------------------------------------------------------
c     
      parameter (zero = 0.0d0)
      parameter (half = 0.5d0)
      parameter (one  = 1.0d0)
      parameter (two  = 2.0d0)
c
c     half integer angular momenta
c
      sj = half*nint(two*ij)
      sk = half*nint(two*ik)
      sm = half*nint(two*im)
c
c     projection ranges
c
      dm = zero
      if (sk.gt.sj .or. sk.lt.-sj)  return
      if (sm.gt.sj .or. sm.lt.-sj)  return
      if (mod(sj-sk,one) .ne. zero) return
      if (mod(sj-sm,one) .ne. zero) return      
c
c     reflection symmetries
c      
      if (sk+sm .ge. zero) then
        if (sk-sm .ge. zero) then
          tk = sk
          tm = sm
          isign = 0
        else
          tk = sm
          tm = sk
          isign = sk-sm
        endif
      else
        if (sk-sm .ge. zero) then
          tk = -sm
          tm = -sk
          isign = 0
        else
          tk = -sk
          tm = -sm
          isign = sk-sm
        endif
      endif
c
c     evaluation
c
      n = sj-tk
      ia = tk-tm
      ib = tk+tm
      a = ia
      b = ib
      beta2 = half*beta
      cosb2 = cos(beta2)
      sinb2 = sin(beta2)
      cosb = (cosb2-sinb2)*(cosb2+sinb2)
      d1 = pjacob(n,a,b,cosb)
      d2 = cosb2**ib*sinb2**ia
      d3 = d1*d2
      d4 = d3*d3
      ti = tm
      do i = 1,ia
         ti = ti+one
         d4 = d4*(sj+ti)/(sj-ti+one)
      enddo
      d4 = sqrt(d4)
      dm = sign(d4,d3)
      if (mod(isign,2) .ne. 0) dm = -dm
      return
      end

      function pjacob (n,a,b,x)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     Jacobi polynomial p(n,a,b;x)
c     Abramowitz and Stegun eq. (22.7.1)
c     ----------------------------------------------------------------- 
c
      parameter (zero = 0.0d0)
      parameter (half = 0.5d0)
      parameter (one  = 1.0d0)
      parameter (two  = 2.0d0)
c
      if (n .eq. 0) then
        fp = one
      else
        f = one
        apa = a+a
        apb = a+b
        amb = a-b
        apbamb = apb*amb
        apbp1 = apb+one
        apbp2 = apb+two
        onek = zero
        twok = zero
        fp = half*(amb+apbp2*x)
        do k = 1,n-1
          onek = onek+one
          twok = twok+two 
          a1 = (twok+two)*(onek+apbp1)*(twok+apb)
          a2 = (twok+apbp1)*apbamb
          a3 = (twok+apb)*(twok+apbp1)*(twok+apbp2)
          a4 = (twok+apa)*(onek+b)*(twok+apbp2)
          fm = f
          f = fp
          fp = ((a2+a3*x)*f-a4*fm)/a1
        enddo
      endif
      pjacob = fp
      return
      end
!
!  ***************************************************************************
!  *                                                                         *
!  *   Program    :  gauscheb                                                *
!  *   Function   :  zeroes and weights of Gauss-Chebyshev quadrature        *
!  *                                                                         *
!  ***************************************************************************
       subroutine    GAUSCHEB(kindb1,nsymb1,nodeb1,xmatb1,wmatb1)
       implicit      none
       integer       kindb1,nsymb1,nodeb1,inodb1
       real*8        xmatb1(nodeb1),wmatb1(nodeb1),b1min,b1max,delta,pi
c
       pi=dacos(-1.0d0)
       if (kindb1.eq.1) then
          b1min=  0;   b1max=2*pi;   else
          b1min=-pi;   b1max= +pi;   endif

       b1min=b1min/nsymb1
       b1max=b1max/nsymb1
       delta=(b1max-b1min)/nodeb1
       do 1200 inodb1=1,nodeb1
          xmatb1(inodb1)=b1min+(1.0d0*inodb1-0.5d0)*delta
          wmatb1(inodb1)=dsqrt(delta/(b1max-b1min))
 1200  continue
       return
       end
!
!  ***************************************************************************
!  *                                                                         *
!  *   Program    :  gausdriv                                                *
!  *   Function   :  the drive routine for gauss quadrature                  *
!  *                                                                         *
!  ***************************************************************************
       subroutine    GAUSDRIV(kind,nodes,xnode,wnode)
       implicit      none
       integer       kind,nodes,memstat,kpts
       real*8        xnode(nodes),wnode(nodes),alpha,beta,endpts(2)
       real*8,       allocatable :: dworks(:)

       allocate(dworks(nodes),stat=memstat)
       if (memstat.ne.0) stop 'failed to allocate memory in gausdriv'
c
c......for gauss-legendre and gauss-chebyshev quadrature

       if (kind.ne.1.and.kind.ne.2) then
          write(*,*)'wrong value of kind when calling GAUSDRIV'
          write(7,*)'wrong value of kind when calling GAUSDRIV'
          stop 
       endif

       kpts=0; alpha=0.0d0; beta=0.0d0

       call GAUSSQ(kind,nodes,alpha,beta,kpts,endpts,dworks,xnode,wnode)
c
c......release memory

       deallocate(dworks,stat=memstat)
       if (memstat.ne.0) stop 'failed to release memory in gausdriv'
       return
       end
c-------------------------------------------------------------------------
c
c     Gauss quadrature 
c
c-------------------------------------------------------------------------
      SUBROUTINE GAUSSQ(KIND, N, ALPHA, BETA, KPTS, ENDPTS, B, T, W)    
      IMPLICIT REAL*8 (A-H,O-Z) 
C                                                                       
C           THIS SET OF ROUTINES COMPUTES THE NODES X(I) AND WEIGHTS    
C        C(I) FOR GAUSSIAN-TYPE QUADRATURE RULES WITH PRE-ASSIGNED      
C        NODES.  THESE ARE USED WHEN ONE WISHES TO APPROXIMATE          
C                                                                       
C                 INTEGRAL (FROM A TO B)  F(X) W(X) DX                  
C                                                                       
C                              N                                         
C        BY                   SUM C  F(X )                               
C                             I=1  I    I                                
C                                                                        
C        HERE W(X) IS ONE OF SIX POSSIBLE NON-NEGATIVE WEIGHT            
C        FUNCTIONS (LISTED BELOW), AND F(X) IS THE                       
C        FUNCTION TO BE INTEGRATED.  GAUSSIAN QUADRATURE IS PARTICULARLY 
C        USEFUL ON INFINITE INTERVALS (WITH APPROPRIATE WEIGHT           
C        FUNCTIONS), SINCE THEN OTHER TECHNIQUES OFTEN FAIL.             
C                                                                        
C           ASSOCIATED WITH EACH WEIGHT FUNCTION W(X) IS A SET OF        
C        ORTHOGONAL POLYNOMIALS.  THE NODES X(I) ARE JUST THE ZEROES     
C        OF THE PROPER N-TH DEGREE POLYNOMIAL.                           
C                                                                        
C     INPUT PARAMETERS                                                   
C                                                                        
C        KIND     AN INTEGER BETWEEN 1 AND 6 GIVING THE TYPE OF          
C                 QUADRATURE RULE                                        
C                                                                        
C        KIND = 1=  LEGENDRE QUADRATURE, W(X) = 1 ON (-1, 1)             
C        KIND = 2=  CHEBYSHEV QUADRATURE OF THE FIRST KIND               
C                   W(X) = 1/SQRT(1 - X*X) ON (-1, +1)                   
C        KIND = 3=  CHEBYSHEV QUADRATURE OF THE SECOND KIND              
C                   W(X) = SQRT(1 - X*X) ON (-1, 1)                      
C        KIND = 4=  HERMITE QUADRATURE, W(X) = EXP(-X*X) ON              
C                   (-INFINITY, +INFINITY)                               
C        KIND = 5=  JACOBI QUADRATURE, W(X) = (1-X)**ALPHA * (1+X)**     
C                   BETA ON (-1, 1), ALPHA, BETA .GT. -1.                
C                   NOTE= KIND=2 AND 3 ARE A SPECIAL CASE OF THIS.       
C        KIND = 6=  GENERALIZED LAGUERRE QUADRATURE, W(X) = EXP(-X)*     
C                   X**ALPHA ON (0, +INFINITY), ALPHA .GT. -1            
C                                                                        
C        N        THE NUMBER OF POINTS USED FOR THE QUADRATURE RULE      
C        ALPHA    REAL PARAMETER USED ONLY FOR GAUSS-JACOBI AND GAUSS-   
C                 LAGUERRE QUADRATURE (OTHERWISE USE 0.).                
C        BETA     REAL PARAMETER USED ONLY FOR GAUSS-JACOBI QUADRATURE-- 
C                 (OTHERWISE USE 0.).                                    
C        KPTS     (INTEGER) NORMALLY 0, UNLESS THE LEFT OR RIGHT END-    
C                 POINT (OR BOTH) OF THE INTERVAL IS REQUIRED TO BE A    
C                 NODE (THIS IS CALLED GAUSS-RADAU OR GAUSS-LOBATTO      
C                 QUADRATURE).  THEN KPTS IS THE NUMBER OF FIXED         
C                 ENDPOINTS (1 OR 2).                                    
C        ENDPTS   REAL ARRAY OF LENGTH 2.  CONTAINS THE VALUES OF        
C                 ANY FIXED ENDPOINTS, IF KPTS = 1 OR 2.                 
C        B        REAL SCRATCH ARRAY OF LENGTH N                         
C                                                                        
C     OUTPUT PARAMETERS (BOTH ARRAYS OF LENGTH N)                        
C                                                                        
C        T        WILL CONTAIN THE DESIRED NODES X(1),,,X(N)             
C        W        WILL CONTAIN THE DESIRED WEIGHTS C(1),,,C(N)           
C                                                                        
C     SUBROUTINES REQUIRED                                               
C                                                                        
C        GBSLVE, CLASS, AND GBTQL2 ARE PROVIDED. UNDERFLOW MAY SOMETIMES 
C        OCCUR, BUT IT IS HARMLESS IF THE UNDERFLOW INTERRUPTS ARE       
C        TURNED OFF AS THEY ARE ON THIS MACHINE.                         
C                                                                        
C     ACCURACY                                                           
C                                                                        
C        THE ROUTINE WAS TESTED UP TO N = 512 FOR LEGENDRE QUADRATURE,   
C        UP TO N = 136 FOR HERMITE, UP TO N = 68 FOR LAGUERRE, AND UP    
C        TO N = 10 OR 20 IN OTHER CASES.  IN ALL BUT TWO INSTANCES,      
C        COMPARISON WITH TABLES IN REF. 3 SHOWED 12 OR MORE SIGNIFICANT  
C        DIGITS OF ACCURACY.  THE TWO EXCEPTIONS WERE THE WEIGHTS FOR    
C        HERMITE AND LAGUERRE QUADRATURE, WHERE UNDERFLOW CAUSED SOME    
C        VERY SMALL WEIGHTS TO BE SET TO ZERO.  THIS IS, OF COURSE,      
C        COMPLETELY HARMLESS.                                            
C                                                                        
C     METHOD                                                             
C                                                                        
C           THE COEFFICIENTS OF THE THREE-TERM RECURRENCE RELATION       
C        FOR THE CORRESPONDING SET OF ORTHOGONAL POLYNOMIALS ARE         
C        USED TO FORM A SYMMETRIC TRIDIAGONAL MATRIX, WHOSE              
C        EIGENVALUES (DETERMINED BY THE IMPLICIT QL-METHOD WITH          
C        SHIFTS) ARE JUST THE DESIRED NODES.  THE FIRST COMPONENTS OF    
C        THE ORTHONORMALIZED EIGENVECTORS, WHEN PROPERLY SCALED,         
C        YIELD THE WEIGHTS.  THIS TECHNIQUE IS MUCH FASTER THAN USING A  
C        ROOT-FINDER TO LOCATE THE ZEROES OF THE ORTHOGONAL POLYNOMIAL.  
C        FOR FURTHER DETAILS, SEE REF. 1.  REF. 2 CONTAINS DETAILS OF    
C        GAUSS-RADAU AND GAUSS-LOBATTO QUADRATURE ONLY.                  
C                                                                        
C     REFERENCES                                                         
C                                                                        
C        1.  GOLUB, G. H., AND WELSCH, J. H.,  CALCULATION OF GAUSSIAN   
C            QUADRATURE RULES,  MATHEMATICS OF COMPUTATION 23 (APRIL,    
C            1969), PP. 221-230.                                         
C        2.  GOLUB, G. H.,  SOME MODIFIED MATRIX EIGENVALUE PROBLEMS,    
C            SIAM REVIEW 15 (APRIL, 1973), PP. 318-334 (SECTION 7).      
C        3.  STROUD AND SECREST, GAUSSIAN QUADRATURE FORMULAS, PRENTICE- 
C            HALL, ENGLEWOOD CLIFFS, N.J., 1966.                         
C                                                                       
C     ..................................................................
C                                                                       
      REAL*8  MUZERO                                                    
      DIMENSION  B(N),T(N),W(N),ENDPTS(2)                               
C                                                                       
      CALL CLASS (KIND, N, ALPHA, BETA, B, T, MUZERO)                   
C                                                                       
C           THE MATRIX OF COEFFICIENTS IS ASSUMED TO BE SYMMETRIC.      
C           THE ARRAY T CONTAINS THE DIAGONAL ELEMENTS, THE ARRAY       
C           B THE OFF-DIAGONAL ELEMENTS.                                
C           MAKE APPROPRIATE CHANGES IN THE LOWER RIGHT 2 BY 2          
C           SUBMATRIX.                                                  
C                                                                       
      IF (KPTS.EQ.0)  GO TO 100                                         
      IF (KPTS.EQ.2)  GO TO  50                                         
C                                                                       
C           IF KPTS=1, ONLY T(N) MUST BE CHANGED                        
C                                                                       
      T(N) =GBSLVE(ENDPTS(1), N, T, B)*B(N-1)**2 + ENDPTS(1)            
      GO TO 100                                                         
C                                                                       
C           IF KPTS=2, T(N) AND B(N-1) MUST BE RECOMPUTED               
C                                                                       
   50 GAM =GBSLVE(ENDPTS(1), N, T, B)                                   
      T1 = ((ENDPTS(1) - ENDPTS(2))/(GBSLVE(ENDPTS(2), N, T, B) - GAM)) 
      B(N-1) =  SQRT(T1)                                                
      T(N) = ENDPTS(1) + GAM*T1                                         
C                                                                       
C           NOTE THAT THE INDICES OF THE ELEMENTS OF B RUN FROM 1 TO N-1
C           AND THUS THE VALUE OF B(N) IS ARBITRARY.                    
C           NOW COMPUTE THE EIGENVALUES OF THE SYMMETRIC TRIDIAGONAL    
C           MATRIX, WHICH HAS BEEN MODIFIED AS NECESSARY.               
C           THE METHOD USED IS A QL-TYPE METHOD WITH ORIGIN SHIFTING    
C                                                                       
  100 W(1) = 1.0D0                                                      
      DO 105 I = 2, N                                                   
  105    W(I) = 0.0D0                                                   
C                                                                       
      CALL GBTQL2 (N, T, B, W, IERR)                                    
      DO 110 I = 1, N                                                   
  110    W(I) = MUZERO * W(I) * W(I)                                    
C                                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C                                                                       
C                                                                       
      DOUBLE PRECISION FUNCTION GBSLVE(SHIFT, N, A, B)                                   
      IMPLICIT REAL*8 (A-H,O-Z) 
C                                                                       
C       THIS PROCEDURE PERFORMS ELIMINATION TO SOLVE FOR THE            
C       N-TH COMPONENT OF THE SOLUTION DELTA TO THE EQUATION            
C                                                                       
C             (JN - SHIFT*IDENTITY) * DELTA  = EN,                      
C                                                                       
C       WHERE EN IS THE VECTOR OF ALL ZEROES EXCEPT FOR 1 IN            
C       THE N-TH POSITION.                                              
C                                                                       
C       THE MATRIX JN IS SYMMETRIC TRIDIAGONAL, WITH DIAGONAL           
C       ELEMENTS A(I), OFF-DIAGONAL ELEMENTS B(I).  THIS EQUATION       
C       MUST BE SOLVED TO OBTAIN THE APPROPRIATE CHANGES IN THE LOWER   
C       2 BY 2 SUBMATRIX OF COEFFICIENTS FOR ORTHOGONAL POLYNOMIALS.    
C                                                                       
C                                                                       
      DIMENSION  A(N),B(N)                                              
C                                                                       
      ALPHA = A(1) - SHIFT                                              
      NM1 = N - 1                                                       
      DO 10 I = 2, NM1                                                  
   10    ALPHA = A(I) - SHIFT - B(I-1)**2/ALPHA                         
      GBSLVE = 1.0D0  /ALPHA                                            
      RETURN                                                            
      END                                                               
C                                                                       
C                                                                       
C                                                                       
      SUBROUTINE CLASS(KIND, N, ALPHA, BETA, B, A, MUZERO)              
      IMPLICIT REAL*8 (A-H,O-Z) 
C                                                                       
C           THIS PROCEDURE SUPPLIES THE COEFFICIENTS A(J), B(J) OF THE  
C        RECURRENCE RELATION                                            
C                                                                       
C             B P (X) = (X - A ) P   (X) - B   P   (X)                  
C              J J            J   J-1       J-1 J-2                     
C                                                                       
C        FOR THE VARIOUS CLASSICAL (NORMALIZED) ORTHOGONAL POLYNOMIALS, 
C        AND THE ZERO-TH MOMENT                                         
C                                                                       
C             MUZERO = INTEGRAL W(X) DX                                 
C                                                                       
C        OF THE GIVEN POLYNOMIAL   WEIGHT FUNCTION W(X).  SINCE THE     
C        POLYNOMIALS ARE ORTHONORMALIZED, THE TRIDIAGONAL MATRIX IS     
C        GUARANTEED TO BE SYMMETRIC.                                    
C                                                                       
C           THE INPUT PARAMETER ALPHA IS USED ONLY FOR LAGUERRE AND     
C        JACOBI POLYNOMIALS, AND THE PARAMETER BETA IS USED ONLY FOR    
C        JACOBI POLYNOMIALS.  THE LAGUERRE AND JACOBI POLYNOMIALS       
C        REQUIRE THE GAMMA FUNCTION.                                    
C                                                                       
C     ..................................................................
C                                                                       
      DIMENSION  A(N),B(N)                                              
      REAL*8  MUZERO                                                    
      DATA PI / 3.141592653589793D0  /                                  
C                                                                       
      NM1 = N - 1                                                       
      GO TO (10, 20, 30, 40, 50, 60), KIND                              
C                                                                       
C              KIND = 1=  LEGENDRE POLYNOMIALS P(X)                     
C              ON (-1, +1), W(X) = 1.                                   
C                                                                       
   10 MUZERO = 2.0D0                                                    
      DO 11 I = 1, NM1                                                  
         A(I) = 0.0D0                                                   
         ABI = I                                                        
   11    B(I) = ABI/ SQRT(4*ABI*ABI - 1.0D0  )                          
      A(N) = 0.0D0                                                      
      RETURN                                                            
C                                                                       
C              KIND = 2=  CHEBYSHEV POLYNOMIALS OF THE FIRST KIND T(X)  
C              ON (-1, +1), W(X) = 1 / SQRT(1 - X*X)                    
C                                                                       
   20 MUZERO = PI                                                       
      DO 21 I = 1, NM1                                                  
         A(I) = 0.0D0                                                   
   21    B(I) = 0.5D0                                                   
      B(1) =  SQRT(0.5D0  )                                             
      A(N) = 0.0D0                                                      
      RETURN                                                            
C                                                                       
C              KIND = 3=  CHEBYSHEV POLYNOMIALS OF THE SECOND KIND U(X) 
C              ON (-1, +1), W(X) = SQRT(1 - X*X)                        
C                                                                       
   30 MUZERO = PI/2.0D0                                                 
      DO 31 I = 1, NM1                                                  
         A(I) = 0.0D0                                                   
   31    B(I) = 0.5D0                                                   
      A(N) = 0.0D0                                                      
      RETURN                                                            
C                                                                       
C              KIND = 4=  HERMITE POLYNOMIALS H(X) ON (-INFINITY,       
C              +INFINITY), W(X) = EXP(-X**2)                            
C                                                                       
   40 MUZERO =  SQRT(PI)                                                
      DO 41 I = 1, NM1                                                  
         A(I) = 0.0D0                                                   
   41    B(I) =  SQRT(I/2.0D0  )                                        
      A(N) = 0.0D0                                                      
      RETURN                                                            
C                                                                       
C              KIND = 5=  JACOBI POLYNOMIALS P(ALPHA, BETA)(X) ON       
C              (-1, +1), W(X) = (1-X)**ALPHA + (1+X)**BETA, ALPHA AND   
C              BETA GREATER THAN -1                                     
C                                                                       
   50 AB = ALPHA + BETA                                                 
      ABI = 2.0D0   + AB                                                
      MUZERO = 2.0D0   ** (AB + 1.0D0  ) * DGAMMA(ALPHA + 1.0D0  ) * DGA
     VMMA(
     X BETA + 1.0D0  ) / DGAMMA(ABI)                                    
      A(1) = (BETA - ALPHA)/ABI                                         
      B(1) =  SQRT(4.0D0  *(1.0D0   + ALPHA)*(1.0D0   + BETA)/((ABI + 1.
     V0D0  )* 
     1  ABI*ABI))                                                       
      A2B2 = BETA*BETA - ALPHA*ALPHA                                    
      DO 51 I = 2, NM1                                                  
         ABI = 2.0D0  *I + AB                                           
         A(I) = A2B2/((ABI - 2.0D0  )*ABI)                              
   51    B(I) =  SQRT (4.0D0  *I*(I + ALPHA)*(I + BETA)*(I + AB)/       
     1   ((ABI*ABI - 1)*ABI*ABI))                                       
      ABI = 2.0D0  *N + AB                                              
      A(N) = A2B2/((ABI - 2.0D0  )*ABI)                                 
      RETURN                                                            
C                                                                       
C              KIND = 6=  LAGUERRE POLYNOMIALS L(ALPHA)(X) ON           
C              (0, +INFINITY), W(X) = EXP(-X) * X**ALPHA, ALPHA GREATER 
C              THAN -1.                                                 
C                                                                       
   60 MUZERO = DGAMMA(ALPHA + 1.0D0  )                                  
      DO 61 I = 1, NM1                                                  
         A(I) = 2.0D0  *I - 1.0D0   + ALPHA                             
   61    B(I) =  SQRT(I*(I + ALPHA))                                    
      A(N) = 2.0D0  *N - 1 + ALPHA                                      
      RETURN                                                            
      END                                                               
C     ------------------------------------------------------------------
C                                                                       
      SUBROUTINE GBTQL2(N, D, E, Z, IERR)                               
      IMPLICIT REAL*8 (A-H,O-Z) 
C                                                                       
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE IMTQL2,   
C     NUM. MATH. 12, 377-383(1968) BY MARTIN AND WILKINSON,             
C     AS MODIFIED IN NUM. MATH. 15, 450(1970) BY DUBRULLE.              
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 241-248(1971).   
C                                                                       
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND FIRST COMPONENTS OF THE 
C     EIGENVECTORS OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE IMPLICIT QL 
C     METHOD, AND IS ADAPTED FROM THE EISPAK ROUTINE IMTQL2             
C                                                                       
C     ON INPUT=                                                         
C                                                                       
C        N IS THE ORDER OF THE MATRIX;                                  
C                                                                       
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX;          
C                                                                       
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX        
C          IN ITS FIRST N-1 POSITIONS.  E(N) IS ARBITRARY;              
C                                                                       
C        Z CONTAINS THE FIRST ROW OF THE IDENTITY MATRIX.               
C                                                                       
C      ON OUTPUT=                                                       
C                                                                       
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN          
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT          
C          UNORDERED FOR INDICES 1, 2, ..., IERR-1;                     
C                                                                       
C        E HAS BEEN DESTROYED;                                          
C                                                                       
C        Z CONTAINS THE FIRST COMPONENTS OF THE ORTHONORMAL EIGENVECTORS
C          OF THE SYMMETRIC TRIDIAGONAL MATRIX.  IF AN ERROR EXIT IS    
C          MADE, Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED 
C          EIGENVALUES;                                                 
C                                                                       
C        IERR IS SET TO                                                 
C                                                                       
C        IERR IS SET TO                                                 
C          ZERO       FOR NORMAL RETURN,                                
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN               
C                     DETERMINED AFTER 30 ITERATIONS.                   
C                                                                       
C     ------------------------------------------------------------------
C                                                                       
      INTEGER I, J, K, L, M, N, II, MML, IERR                           
      DIMENSION  D(N),E(N),Z(N)                                         
      REAL*8  MACHEP                                                    
C                                                                       
C     ========== MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING     
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.   
C                MACHEP = 16.0D0**(-13) FOR LONG FORM ARITHMETIC        
C                ON S360 ==========                                     
       MACHEP=1.0D-14                                                   
C                                                                       
      IERR = 0                                                          
      IF (N .EQ. 1) GO TO 1001                                          
C                                                                       
      E(N) = 0.0D0                                                      
      DO 240 L = 1, N                                                   
         J = 0                                                          
C     ========== LOOK FOR SMALL SUB-DIAGONAL ELEMENT ==========         
  105    DO 110 M = L, N                                                
            IF (M .EQ. N) GO TO 120                                     
            IF ( ABS(E(M)) .LE. MACHEP * ( ABS(D(M)) +  ABS(D(M+1))))   
     X         GO TO 120                                                
  110    CONTINUE                                                       
C                                                                       
  120    P = D(L)                                                       
         IF (M .EQ. L) GO TO 240                                        
         IF (J .EQ. 30) GO TO 1000                                      
         J = J + 1                                                      
C     ========== FORM SHIFT ==========                                  
         G = (D(L+1) - P) / (2.0D0   * E(L))                            
         R =  SQRT(G*G+1.0D0  )                                         
         G = D(M) - P + E(L) / (G +  SIGN(R, G))                        
         S = 1.0D0                                                      
         C = 1.0D0                                                      
         P = 0.0D0                                                      
         MML = M - L                                                    
C     ========== FOR I=M-1 STEP -1 UNTIL L DO -- ==========             
         DO 200 II = 1, MML                                             
            I = M - II                                                  
            F = S * E(I)                                                
            B = C * E(I)                                                
            IF ( ABS(F) .LT.  ABS(G)) GO TO 150                         
            C = G / F                                                   
            R =  SQRT(C*C+1.0D0  )                                      
            E(I+1) = F * R                                              
            S = 1.0D0   / R                                             
            C = C * S                                                   
            GO TO 160                                                   
  150       S = F / G                                                   
            R =  SQRT(S*S+1.0D0  )                                      
            E(I+1) = G * R                                              
            C = 1.0D0   / R                                             
            S = S * C                                                   
  160       G = D(I+1) - P                                              
            R = (D(I) - G) * S + 2.0D0   * C * B                        
            P = S * R                                                   
            D(I+1) = G + P                                              
            G = C * R - B                                               
C     ========== FORM FIRST COMPONENT OF VECTOR ==========              
            F = Z(I+1)                                                  
            Z(I+1) = S * Z(I) + C * F                                   
            Z(I) = C * Z(I) - S * F                                     
C                                                                       
  200    CONTINUE                                                       
C                                                                       
         D(L) = D(L) - P                                                
         E(L) = G                                                       
         E(M) = 0.0D0                                                   
         GO TO 105                                                      
  240 CONTINUE                                                          
C     ========== ORDER EIGENVALUES AND EIGENVECTORS ==========          
      DO 300 II = 2, N                                                  
         I = II - 1                                                     
         K = I                                                          
         P = D(I)                                                       
C                                                                       
         DO 260 J = II, N                                               
            IF (D(J) .GE. P) GO TO 260                                  
            K = J                                                       
            P = D(J)                                                    
  260    CONTINUE                                                       
C                                                                       
         IF (K .EQ. I) GO TO 300                                        
         D(K) = D(I)                                                    
         D(I) = P                                                       
C                                                                       
         P = Z(I)                                                       
         Z(I) = Z(K)                                                    
         Z(K) = P                                                       
C                                                                       
  300 CONTINUE                                                          
C                                                                       
      GO TO 1001                                                        
C     ========== SET ERROR -- NO CONVERGENCE TO AN                      
C                EIGENVALUE AFTER 30 ITERATIONS ==========              
 1000 IERR = L                                                          
 1001 RETURN                                                            
C     ========== LAST CARD OF GBTQL2 ==========                         
      END                                                               
      DOUBLE PRECISION FUNCTION  DGAMMA(Z)                                               
      IMPLICIT REAL*8 (A-H,O-Z) 
C  THIS IS A PROCEDURE THAT EVALUATES GAMMA(Z) FOR                      
C     0 LT Z LE 3 TO 16 SIGNIFICANT FIGURES                             
C    IT IS BASED ON A CHEBYSHEV-TYPE POLYNOMIAL                         
C   APPROXIMATION GIVEN IN H. WERNER AND R. COLLINGE, MATH. COMPUT.     
C    15 (1961), PP. 195-97.                                             
C   APPROXIMATIONS TO THE GAMMA FUNCTION, ACCURATE UP TO 18 SIGNIFICANT 
C   DIGITS, MAY BE FOUND IN THE PAPER QUOTED ABOVE                      
C                                                                       
C                                                                       
C                                                                       
      DIMENSION  A(18)                                                  
C                                                                       
       A(1)=1.0D0                                                       
       A(2)=.4227843350984678D0                                         
       A(3)=.4118403304263672D0                                         
      A(4)=.0815769192502609D0                                          
      A(5)=.0742490106800904D0                                          
      A(6)=-.0002669810333484D0                                         
      A(7)=.0111540360240344D0                                          
      A(8)=-.0028525821446197D0                                         
      A(9)=.0021036287024598D0                                          
      A(10)=-.0009184843690991D0                                        
      A(11)=.0004874227944768D0                                         
      A(12)=-.0002347204018919D0                                        
      A(13)=.0001115339519666D0                                         
      A(14)=-.0000478747983834D0                                        
      A(15)=.0000175102727179D0                                         
      A(16)=-.0000049203750904D0                                        
      A(17)=.0000009199156407D0                                         
      A(18)=-.0000000839940496D0                                        
C                                                                       
C                                                                       
C                                                                       
      IF(Z.LE.1.0D0  ) GO TO 10                                         
      IF(Z.LE.2.0D0  ) GO TO 20                                         
      T=Z-2.0D0                                                         
      GO TO 30                                                          
10    T=Z                                                               
      GO TO 30                                                          
20    T=Z-1.0D0                                                         
30    P=A(18)                                                           
      DO 40 K1=1,17                                                     
      K=18-K1                                                           
      P=T*P+A(K)                                                        
40    CONTINUE                                                          
C                                                                       
      IF(Z.GT.2.0D0  ) GO TO 50                                         
      IF(Z.GT.1.0D0  ) GO TO 60                                         
      DGAMMA=P/(Z*(Z+1.0D0  ))                                          
      RETURN                                                            
60    DGAMMA=P/Z                                                        
      RETURN                                                            
50    DGAMMA=P                                                          
      RETURN                                                            
      END                                                               
c------------------------------------------------------------
        DOUBLE PRECISION FUNCTION FUNCTCG(J1,M1,J2,M2,J3)
        IMPLICIT REAL*8(A-H,O-Z)
        FUNCTCG=0
        IF (ABS(M1).GT.J1) RETURN
        IF (ABS(M2).GT.J2) RETURN
        IF (ABS(M1+M2).GT.J3) RETURN
        IF (J3.LT.ABS(J1-J2)) RETURN
        IF (J3.GT.J1+J2) RETURN

        M3=M1+M2
        M33=-M3
        FUNCTCG=(-1)**(J1+J2+M33)*F3J(J1,J2,J3,M1,M2,M33,2)
     $  *DSQRT(2.D0*J3+1.D0)
        RETURN
        END

      DOUBLE PRECISION FUNCTION F3J(JD1,JD2,JD3,MD1,MD2,MD3,ITWO)
C FROM NBS TECHNICAL NOTE 409
C
C CALCULATE 3J SYMBOL. THIS FUNCTION WORKS WITH BOTH INTEGRAL AND HALF
C INTEGRAL ARGUMENTS. IF ALL INTEGRAL ARGUMENTS, USE ITWO=2 AND JD1,ETC
C EQUAL TO THE ARGUMENTS. IF SOME HALF INTEGRAL ARGUMENTS, USE ITWO=1
C AND JD1,ETC EQUAL TO TWICE THE ARGUMENTS.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FACT3J/FL(322),NCALL
      DIMENSION MTRI(9)
      DATA ZERO,HALF,ONE/0.0D+0,0.5D+0,1.0D+0/
      DATA EPS,EPS2,EPS3,EPS4,EPS5,EPS6/1.0D-10,1.0D+30,1.0D-30,8.0D+1,
     $1.0D+10,23.02585092994046D+0/
      SAVE /FACT3J/
      J1=JD1*ITWO
      J2=JD2*ITWO
      J3=JD3*ITWO
      M1=MD1*ITWO
      M2=MD2*ITWO
      M3=MD3*ITWO
      IF(NCALL+1867)5,15,5
5     NCALL=-1867
      FL(1)=ZERO
      FL(2)=ZERO
      DO 50 N=3,322
      FN=N-1
50    FL(N)=FL(N-1)+DLOG(FN)
15    I=J1+J2-J3
      I1=I/2
      IF(I-2*I1)1000,1010,1000
1010  MTRI(1)=I1
      I=J1-J2+J3
      I1=I/2
      IF(I-2*I1)1000,1020,1000
1020  MTRI(2)=I1
      I=-J1+J2+J3
      I1=I/2
      IF(I-2*I1)1000,1030,1000
1030  MTRI(3)=I1
      IF(M1+M2+M3)1000,1040,1000
1040  I=J1+M1
      I1=I/2
      IF(I-2*I1)1000,1050,1000
1050  MTRI(4)=I1
      MTRI(5)=(J1-M1)/2
      I=J2+M2
      I1=I/2
      IF(I-2*I1)1000,1060,1000
1060  MTRI(6)=I1
      MTRI(7)=(J2-M2)/2
      I=J3+M3
      I1=I/2
      IF(I-2*I1)1000,1070,1000
1070  MTRI(8)=I1
      MTRI(9)=(J3-M3)/2
      DO 30 N=1,9
      IF(MTRI(N))1000,30,30
30    CONTINUE
      IF(J3-J2+M1)40,45,45
40    KMIN=-J3+J2-M1
      GOTO 60
45    KMIN=0
60    IF(-J3+J1+M2-KMIN)80,80,70
70    KMIN=-J3+J1+M2
80    KMIN=KMIN/2
      IF(J2-J3+M1)90,100,100
90    KMAX=J1+J2-J3
      GOTO 110
100   KMAX=J1-M1
110   IF(J2+M2-KMAX)120,130,130
120   KMAX=J2+M2
130   KMAX=KMAX/2
      MIN1=MTRI(1)-KMIN+1
      MIN2=MTRI(5)-KMIN+1
      MIN3=MTRI(6)-KMIN+1
      MIN4=(J3-J2+M1)/2+KMIN
      MIN5=(J3-J1-M2)/2+KMIN
      UK=EPS
      S=UK
      NCUT=0
      KMAX=KMAX-KMIN
      IF(KMAX)165,165,155
155   DO 160 K=1,KMAX
      UK=-UK*DFLOAT(MIN1-K)*DFLOAT(MIN2-K)*DFLOAT(MIN3-K)/(DFLOAT(KMIN+
     1K)*DFLOAT(MIN4+K)*DFLOAT(MIN5+K))
      IF(ABS(UK)-EPS2)158,157,157
157   UK=EPS*UK
      S=EPS*S
      NCUT=NCUT+1
158   IF(ABS(UK)-EPS3)165,160,160
160   S=S+UK
165   DELOG=ZERO
      DO 170 N=1,9
      NUM=MTRI(N)
170   DELOG=DELOG+FL(NUM+1)
      NUM=(J1+J2+J3)/2+2
      DELOG=HALF*(DELOG-FL(NUM))
      ULOG=-FL(KMIN+1)-FL(MIN1)-FL(MIN2)-FL(MIN3)-FL(MIN4+1)-FL(MIN5+1)
      PLOG=DELOG+ULOG
      IF(PLOG+EPS4)172,171,171
171   IF(NCUT)175,175,172
172   SIG=SIGN(ONE,S)
      S=ABS(S)
      SLOG=DLOG(S)+DFLOAT(NCUT+1)*EPS6
      F3J=SIG*EXP(SLOG+PLOG)
      GOTO 178
175   S=S*EPS5
      P=EXP(PLOG)
      F3J=P*S
178   NUM=KMIN+(J1-J2-M3)/2
      IF(MOD(NUM,2))180,190,180
180   F3J=-F3J
190   CONTINUE
      GOTO 2000
1000  F3J=ZERO
2000  RETURN
      END
      DOUBLE PRECISION FUNCTION FCOEF(J,L,JP,LP,JTOT,LAMBDA)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C     CALCULATE PERCIVAL-SEATON COEFFICIENT.
C
      DATA ZERO/0.0D+0/
      FCOEF=ZERO
      F1=F3J(J,LAMBDA,JP,0,0,0,2)
      IF(F1.EQ.ZERO)RETURN
      F2=F3J(L,LAMBDA,LP,0,0,0,2)
      IF(F2.EQ.ZERO)RETURN
      F3=F3J(L,J,JTOT,JP,LP,LAMBDA,2)
c6      F3=F6J(L,J,JTOT,JP,LP,LAMBDA,2)
      IF(F3.EQ.ZERO)RETURN
      FCOEF=F1*F2*F3*DSQRT(DFLOAT((2*J+1)*(2*JP+1)*(2*L+1)*(2*LP+1)))
      IF(MOD(J+L+LAMBDA,2).NE.0)FCOEF=-FCOEF
      RETURN
      END
c-----------------------------------------------------------------
      SUBROUTINE RBES(N,Z,ZJN,ZJNP,ZYN,ZYNP)
      IMPLICIT REAL*8 (A-H,O-Z)
C***********************************************************************
C***  THIS SUBROUTINE CALCULATES THE REGULAR RICATTI-BESSEL FUNCTION ZJN
C***  AND ITS DERIVATIVE ZJNP AND THE IRREGULAR RICATTI-BESSEL FUNCTION
C***  ZYN AND ITS DERIVATIVE ZYNP OF ORDER N FOR ARGUMENT Z.  THIS IS
C***  A COPY OF LIGHTS ROUTINE ATTRIBUTED TO R.GORDON.
C***********************************************************************
      DATA ZERO,HALF,ONE,TWO,THREE/0.D0,0.5D0,1.D0,2.D0,3.D0/
C***  EVALUATE TRIG FUNCTIONS
      SZ=DSIN(Z)
      CZ=DCOS(Z)
      ZINV=ONE/Z
      AJ=SZ
      BJ=SZ*ZINV-CZ
      AY=-CZ
      BY=-SZ-CZ*ZINV
      IF(N-1)1,2,3
C***  FUNCTIONS OF ZERO ORDER
    1 ZJN=AJ
      ZJNP=CZ
      ZYN=AY
      ZYNP=SZ
      RETURN
C***  FUNCTIONS OF ORDER ONE
    2 ZJN=BJ
      ZJNP=AJ-ZINV*BJ
      ZYN=BY
      ZYNP=AY-ZINV*BY
      RETURN
C***  FUNCTIONS OF ORDER GREATER THAN ONE
    3 XN=N
      DELF=ZINV+ZINV
      FACTOR=DELF+ZINV
C***  TEST TO SEE IF RECURRENCE CAN BE USED
      IF(Z.GT.XN)GO TO 100
C***  IN NON CLASSICAL REGION USE RECURSION RELATION FOR IRREGULAR
C***  SOLUTION AND ASCENDING POWER SERIES FOR REGULAR SOLUTION.
C***  THE COEFFICIENT,C, IN FRONT OF POWER SERIES IS CALCULATED DURING
C***  FORWARD RECURRENCE.
      ZSQ=Z*Z
      TOP=XN+XN
      XJ=THREE
      C=ZSQ/THREE
C***  IRREGULAR SOLUTION
   25 ZYN=BY*FACTOR-AY
      AY=BY
      BY=ZYN
      FACTOR=FACTOR+DELF
      XJ=XJ+TWO
      C=C*(Z/XJ)
      IF(XJ.LT.TOP)GO TO 25
      ZYNP=AY-ZINV*BY*XN
C***  REGULAR SOLUTION
      FACTOR=-HALF*ZSQ
      U=ONE
      TERM=ONE
      DU=ONE
      DTERM=ONE
      D2NP3=TOP+THREE
      DEN=D2NP3
      DFACT=ZERO
   35 DFACT=DFACT+ONE
      TERM=TERM*(FACTOR/(DFACT*DEN))
      U=U+TERM
      DEN=DEN+TWO
      DTERM=DTERM*(FACTOR/(DFACT*DEN))
      DU=DU+DTERM
C***  TEST CONVERGENCE TO SINGLE PRECISION ACCURACY
      IF(ABS(TERM).GT.1.0D-9)GO TO 35
      ZJN=U*C
      ZJNP=(XN+ONE)*ZJN*ZINV-(Z*DU/D2NP3)*C
      RETURN
C***  FOR CLASSICAL CASE USE FORWARD RECURSION FOR REGULAR
C***  AND IRREGULAR SOLUTIONS
  100 CONST=ZINV*XN
      TOP=CONST+CONST
  200 AJ=FACTOR*BJ-AJ
      AY=FACTOR*BY-AY
      FACTOR=FACTOR+DELF
      IF(FACTOR.GT.TOP)GO TO 250
      BJ=FACTOR*AJ-BJ
      BY=FACTOR*AY-BY
      FACTOR=FACTOR+DELF
      IF(FACTOR.LT.TOP)GO TO 200
      ZJN=BJ
      ZJNP=AJ-CONST*BJ
      ZYN=BY
      ZYNP=AY-CONST*BY
      RETURN
  250 ZJN=AJ
      ZJNP=BJ-CONST*AJ
      ZYN=AY
      ZYNP=BY-CONST*AY
      END
!
!  ***************************************************************************
!  *                                                                         *
!  *   Program    :  legend0                                                 *
!  *   Function   :  associate legendre function                             * 
!  *                                                                         *
!  ***************************************************************************
       subroutine LEGEND0(maxnod,nodea1,ja1max,xnoda1,wnoda1,legenda1) 
       implicit none
       integer  :: maxnod,nodea1,ja1max
       real*8   :: xnoda1(nodea1),wnoda1(nodea1),
     &             legenda1(maxnod,-ja1max:ja1max,0:ja1max) 
c
!       call DZERO1D(maxnod*(2*ja1max+1)*(ja1max+1),legenda1) 
       legenda1=0.d0
       call SPHAM(maxnod,nodea1,ja1max,legenda1,xnoda1,wnoda1) 
       return
       end 
c=========================================================================
        SUBROUTINE SPHAM(NQA12max,NQA12,J12max,YL,AN,AW)
C   This routine gives the spherical hamonic function Yjm at all quadreture
C   nodes (DSQRT(Weight) is also included).
        implicit real*8 (a-h,o-z)
        PARAMETER (MN1=300,MN2=150)
        dimension AF(MN1),YL(NQA12max,-J12max:J12max,0:J12max),
     $  AN(NQA12),AW(NQA12),YL1(MN2,MN2)

        IF(2*J12max+2.GT.MN1) STOP 'increase MN1 in subroutine SPHAM'
        call NFACTO(2*J12max+2,AF,0)

        IF(J12max+1.GT.MN2) STOP 'increase MN2 in subroutine SPHAM'

        DO 20 NQ=1,NQA12
        x=DCOS(AN(NQ))
        call S3YL(x,MN2,J12max+1,YL1,J12max+1,AF)
        DO 40 JJ=1,J12max+1
        DO 40 M1=1,JJ
        YL(NQ,M1-1,JJ-1)=YL1(JJ,M1)*DSQRT(AW(NQ))
40      YL(NQ,1-M1,JJ-1)=(-1)**(M1-1)*YL1(JJ,M1)*DSQRT(AW(NQ))
20      CONTINUE
        RETURN
        END
c------------------------------------------------------------
C In the routine S3YL, you need to provide a factorial array AF which
C can be generated by calling the routine NFACTO in your calling program
         SUBROUTINE S3YL(X,MM,N,YL,NN,AF)
         IMPLICIT REAL*8(A-H,O-Z)
         DIMENSION YL(MM,N),AF(2*MM)
C*******************************************************************
C        THIS SUBROUTINE CALCULATES SPHERICAL HARMONICS
C        YL(L+1,M+1) IS A SPHERICAL HARMONIC OF ORDER L,M AT A
C        POINT X
C        CALCULATION IS DONE BY MULTIPLYING THE RESULTS OF
C        SUBROUTINE S3LG BY A FACTOR WHICH IS
C        THE SQARE ROOT OF (L+1/2)*(L-M)!/(L+M)!
C        NOTE: THE FACTOR OF 1./(SQARE ROOT OF 2*3.14...) IS NOT
C        INCLUDED IN THIS DEFINITION
C**********************************************************************
         IF(N.GT.MM) STOP 'S3YL1'
         BETA=DACOS(X)
         NAF=2*NN
         CALL INIT2(YL,MM,NN)
         CALL S3LG(BETA,MM,N,NN,YL,AF,NAF)
         DO 3 M=1,NN
         DO 3 L=M,N
         F1=DSQRT(DFLOAT(L)-0.5D0)
    3    YL(L,M)=F1*DSQRT(AF(L-M+1)/AF(L+M-1))*YL(L,M)
         RETURN
         END
         SUBROUTINE S3LG(BETA,MMAX,NMAX,NN,P,AF,NAF)
         IMPLICIT REAL*8(A-H,O-Z)
         DIMENSION P(MMAX,NMAX),AF(NAF)
C*******************************************************************
C        P(N+1,M+1) IS THE ASSOCIATED LEGENDRE FUNCTION OF ORDER N,M
C************************************************************************
         IF((2*NN-1).GT.NAF) STOP 'S3LG'
         X=DCOS(BETA)
         Y=-DSIN(BETA)
         NDO=NMAX-2
         DO 3 M=1,NN
         MV=M-1
         A1=DFLOAT((2*MV+2)*(2*MV+1))/DFLOAT(2*(MV+1))
         P(M,M)=AF(2*MV+1)/AF(M)/2.D0**MV
         IF(MV.NE.0) P(M,M)=Y**MV*P(M,M)
         IF(NMAX.GT.M) p(M+1,M)=X*A1*P(M,M)
         DO 3 N=M,NDO
    3    P(N+2,M)=((2*N+1)*X*P(N+1,M)-(N+M-1)*P(N,M))/(N-M+2)
         RETURN
         END
      SUBROUTINE INIT2(A,MN,N2)
C*******************************************************************
C     SUBROUTINE INIT2 INITIALIZE A TWO DIMENSIONAL ARRAY.
C*********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(MN,N2)
      DO 100 J=1,N2
      DO 100 I=1,MN
  100 A(I,J)=0.D0
      RETURN
      END

      SUBROUTINE NFACTO(N,A,IOP)
      IMPLICIT REAL*8(A-H,O-Z)
C*************************************************************************
C  IF IOP=0: CALCULATES THE NFACTORIAL AN(I)=(N-1)!
C  IF IOP GREATER THAN 0: CALCULATES THE LOG OF FACTORIAL AN(I)=DLOG((N-1)!)
C*************************************************************************
      DIMENSION A(N)
      IF(N.LT.1) STOP 'N IS LESS THAN 1 IN NFACTO'
      IF(IOP.GT.0) GO TO 3
      A(1)=1.D0
      DO 1 I=2,N
      A(I)=DFLOAT(I-1)*A(I-1)
   1  CONTINUE
      RETURN
   3  CONTINUE
      A(1)=0.D0
      DO 2 I=2,N
      A(I)=A(I-1)+DLOG(DFLOAT(I-1))
   2  CONTINUE
      RETURN
      END
c==========================================================================
      subroutine lobatto ( n, x, w )

!*****************************************************************************80
!
!! LOBATTO_COMPUTE computes a Lobatto quadrature rule.
!
!  Discussion:
!
!    The integral:
!
!      Integral ( -1 <= X <= 1 ) F(X) dx
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
!
!    The quadrature rule will integrate exactly all polynomials up to
!    X**(2*N-3).
!
!    The Lobatto rule is distinguished by the fact that both endpoints
!    (-1 and 1) are always abscissas.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2007
!
!  Author:
!
!    Original MATLAB version by Greg von Winckel.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Claudio Canuto, Yousuff Hussaini, Alfio Quarteroni, Thomas Zang,
!    Spectral Methods in Fluid Dynamics,
!    Springer, 1993,
!    ISNB13: 978-3540522058,
!    LC: QA377.S676.
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996,
!    ISBN: 0-8493-2479-3,
!    LC: QA47.M315.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    N must be at least 2.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
      implicit none

      integer ( kind = 4 ) n

      integer ( kind = 4 ) i
      integer ( kind = 4 ) j
      real ( kind = 8 ) p(n,n)
      real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
      real ( kind = 8 ) tolerance
      real ( kind = 8 ) w(n)
      real ( kind = 8 ) x(n)
      real ( kind = 8 ) xold(n)

      if ( n < 2 ) then
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'LOBATTO_COMPUTE - Fatal error!'
         write ( *, '(a,i8)' ) '  Illegal value of N = ', n
         write ( *, '(a)' ) ' N must be at least 2.'
         stop
      end if

      tolerance = 100.0D+00 * epsilon ( tolerance )
!
!  Initial estimate for the abscissas is the Chebyshev-Gauss-Lobatto nodes.
!
      do i = 1, n
         x(i) = cos ( pi * real ( i - 1, kind = 8 ) 
     &        / real ( n - 1, kind = 8 ) )
      end do

      xold(1:n) = 2.0D+00

      do while ( tolerance < maxval ( abs ( x(1:n) - xold(1:n) ) ) )

        xold(1:n) = x(1:n)

        p(1:n,1) = 1.0D+00
        p(1:n,2) = x(1:n)
        do j = 2, n-1
          p(1:n,j+1) = ( real ( 2 * j - 1, kind = 8 ) * x(1:n) 
     &               * p(1:n,j) + real (   - j + 1, kind = 8 )
     &               * p(1:n,j-1) ) / real ( j, kind = 8 )
        end do

        x(1:n) = xold(1:n) - ( x(1:n) * p(1:n,n) - p(1:n,n-1) )
     &         / ( real ( n, kind = 8 ) * p(1:n,n) )
      end do

      x(1:n) = x(n:1:-1)
      w(1:n) = 2.0D+00 / ( real ( ( n - 1 ) * n, kind = 8 ) 
     &       * p(1:n,n)**2 )

      return
      end
