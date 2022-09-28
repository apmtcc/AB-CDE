        subroutine READINP
        use fileIO
        use COMPARAS
        
        implicit none
        
        integer :: i, ierr
        integer NA1PLOT, NA2PLOT, NA3PLOT
        
        EV=27.2114d0
        CM=219474.6d0
        PI=DACOS(-1.d0)
        
        read(*,*), IDPOT,VCUT, IDINIT
        read(*,*), ATOM
        read(*,*), NZ, NZVB, Zmin, Zmax   
        read(*,*), NR1, NVB1, R1min, R1max
        read(*,*), NR2, NVB2, R2min, R2max
        read(*,*), NR3, NVB3, R3min, R3max
        read(*,*), Jtot, KBZmax, parity 
        read(*,*), J1max, JS1, JD1
        read(*,*), J2max, JS2, JD2
        read(*,*), J3max, JS3, JD3
        read(*,*), J123max, J23max 
        read(*,*), NSYM
        read(*,*), Z0, R10, R20, R30, A10, A20, A30, B10, B20    
        read(*,*), istat,nstat
        read(*,*), NA1PLOT, NA2PLOT, NA3PLOT

        print*, IDPOT,VCUT,IDINIT
        print*, ATOM
        print*, NZ, NZVB, Zmin, Zmax   
        print*, NR1, NVB1, R1min, R1max
        print*, NR2, NVB2, R2min, R2max
        print*, NR3, NVB3, R3min, R3max
        print*, Jtot, KBZmax, parity 
        print*, J1max, JS1, JD1
        print*, J2max, JS2, JD2
        print*, J3max, JS3, JD3
        print*, J123max, J23max 
        print*, NSYM
        print*, Z0, R10, R20, R30, A10, A20, A30, B10, B20    
        print*, istat,nstat
        print*, NA1PLOT, NA2PLOT, NA3PLOT

        call MASS(5,ATOM,AM)

        RMS_ABCDE=(AM(1)+AM(2))*(AM(3)+AM(4)+AM(5))/
     .    (AM(1)+AM(2)+AM(3)+AM(4)+AM(5))
        RMS_AB=AM(1)*AM(2)/(AM(1)+AM(2))
        RMS_CD=AM(4)*AM(5)/(AM(4)+AM(5))
        RMS_ECD=AM(3)*(AM(4)+AM(5))/(AM(3)+AM(4)+AM(5))

        PRINT*,'RMS_ABCDE=',RMS_ABCDE,'RMS_AB=', RMS_AB
        PRINT*, 'RMS_CD=', RMS_CD, 'RMS_ECD=', RMS_ECD
        
        A10=A10*PI/180.d0
        A20=A20*PI/180.d0
        A30=A30*PI/180.d0
        B10=B10*PI/180.d0
        B20=B20*PI/180.d0
        
        VCUT=VCUT/CM
        PRINT *, 'VCUT=', VCUT

! |J1max-J2max| <= J12max <= (J1max+J2max)
! |J2max-J3max| <= L23max <= (J2max+J3max)
        if (JD1.eq.2) J1max=J1max+mod(J1max-JS1,2)
        if (JD3.eq.2) J3max=J3max+mod(J3max-JS3,2)
        
        KB1max=min(J1max, J2max)
        KM2max=min(J2max, J3max)
        print*, 'KB1max=',KB1max,'KM2max=',KM2max
        
        NA1node=J1max/JD1+1
        NA2node=(J2max+1)/JD2
        NA3node=J3max/JD3+1
        NB1node=(2*KB1max)/NSYM+1
        NB2node=2*KM2max+1
        
        if(nstat.ne.0) then
        NA1node=NA1node*2*NA1PLOT+1
        NA2node=NA2node*NA2PLOT
        NA3node=NA3node*2*NA3PLOT+1
        end if
        
        print*, 'NA1node=', NA1node, 
     .    'NA2node=', NA2node,'NA3node=', NA3node
        print*, 'NB1node=', NB1node, 
     .    'NB2node=', NB2node

        nodeax=NA1node*NA2node*NA3node
        nodebx=NB1node*NB2node
        NANODES=nodeax*nodebx
        NRGRIDS=NZVB*NVB1*NVB2*NVB3
        PRINT *, 'nodeax=', nodeax, 'nodebs=', nodebx
        PRINT *, 'NANODES=', NANODES,  'NRGRIDS', NRGRIDS

        return
        end subroutine READINP
c====================================================================
        subroutine ZBASIS
        use COMPARAS, only: NZ, NZVB, Z0, Zmin, Zmax, RMS_ABCDE
        use QUADPOINTS, only : ZQ
        use TRANSMATS, only : TMATZ
        use SMALL_ARRAYS, only : ETZ, VPOT1D, RQ0Z, HARZ, VrefZ
        
        implicit none
        integer ierr
        integer i, j, in
        real*8 DZ, PI, FAC
        
        allocate(TMATZ(NZVB,NZVB),ZQ(NZVB),ETZ(NZVB),VrefZ(NZVB),
     .  stat=ierr)
        if (NZVB.eq.1) then
          ZQ(NZVB)=Z0
          ETZ(NZVB)=0.d0
          TMATZ(NZVB,NZVB)=1.d0
          VrefZ(NZVB)=0.d0
          return
        end if
        
        DZ=(Zmax-Zmin)/(NZVB+1)
        DO I=1,NZVB
          ZQ(I)=Zmin+I*DZ
        END DO

        call SINBASIS(NZVB,Zmin,Zmax,RMS_ABCDE,TMATZ,ETZ)
        VrefZ=0.D0

        allocate(HARZ(NZ,NZ),RQ0Z(NZ),stat=ierr)
        PI=DACOS(-1.D0)
        FAC=DSQRT(2.D0/(NZ+1))
        DO I=1,NZ
        DO J=1,NZ
           HARZ(J,I)=FAC*DSIN(I*J*PI/(NZ+1))
        END DO
        END DO

        DZ=(Zmax-Zmin)/(NZ+1)
        DO I=1,NZ
          RQ0Z(I)=Zmin+I*DZ
        END DO
       
        do in=1,NZVB
          print*, in, ZQ(in)
        enddo
 
        end subroutine ZBASIS
c====================================================================
        subroutine R1BASIS
        use COMPARAS, only: NR1, NVB1, R10, R1min, R1max, RMS_AB, IDINIT
        use QUADPOINTS, only : RQ1
        use TRANSMATS, only : TMATR1
        use SMALL_ARRAYS, only : ETR1, VPOT1D, RQ0R1, HAR1, VrefR1
        
        implicit none
        integer ierr
        integer i, in
        
        allocate(TMATR1(NVB1,NVB1),RQ1(NVB1),ETR1(NVB1),VrefR1(NVB1),
     .  stat=ierr)
        if (NVB1.eq.1) then
          RQ1(NVB1)=R10
          ETR1(NVB1)=0.d0
          TMATR1(NVB1,NVB1)=1.d0
          VrefR1(NVB1)=0.d0
          return
        end if
        
        allocate(VPOT1D(NR1),RQ0R1(NR1),stat=ierr)
        call POT1D(NR1,R1min,R1max,RQ0R1,VPOT1D,'R1')

        open(28,file='r1.dat',status='unknown') 
        do in=1,NR1
          write(28,*) RQ0R1(in), VPOT1D(in)*27.2114d0
        enddo         
        close(28) 
 
       
        allocate(HAR1(NR1,NR1),stat=ierr)
        CALL VIB1D(RMS_AB,NR1,RQ0R1,VPOT1D,ETR1,HAR1,
     .  NVB1,NVB1,RQ1,RQ1,TMATR1,TMATR1,VrefR1)
        do in=1,NVB1
          print*, in, RQ1(in), VrefR1(in)*27.2114d0
        enddo
        
        deallocate(VPOT1D,stat=ierr)
        
        end subroutine R1BASIS
c====================================================================
        subroutine R2BASIS
        use COMPARAS, only: NR2, NVB2, R20, R2min, R2max, RMS_ECD,IDINIT
        use QUADPOINTS, only : RQ2
        use TRANSMATS, only : TMATR2
        use SMALL_ARRAYS, only : ETR2, VPOT1D, RQ0R2, HAR2, VrefR2
        
        implicit none
        integer i, j, ierr
        integer in
        real*8 DR, FAC, PI

        allocate(TMATR2(NVB2,NVB2),RQ2(NVB2),ETR2(NVB2),VrefR2(NVB2),
     .  stat=ierr)
        if (NVB2.eq.1) then
          RQ2(NVB2)=R20
          ETR2(NVB2)=0.d0
          TMATR2(NVB2,NVB2)=1.d0
          VrefR2(NVB2)=0.d0
          return
        end if

        allocate(VPOT1D(NR2),RQ0R2(NR2),stat=ierr)
        
        call POT1D(NR2,R2min,R2max,RQ0R2,VPOT1D,'R2')
        open(30,file='r2.dat',status='unknown') 
        do in=1,NR2
          write(30,*) RQ0R2(in), VPOT1D(in)*27.2114D0
        enddo
        close(30)

        allocate(HAR2(NR2,NR2),stat=ierr)
        CALL VIB1D(RMS_ECD,NR2,RQ0R2,VPOT1D,ETR2,HAR2,
     .  NVB2,NVB2,RQ2,RQ2,TMATR2,TMATR2,VrefR2)
         
        do in=1,NVB2
          print*, in, RQ2(in), VrefR2(in)*27.2114d0
        enddo
        
        deallocate(VPOT1D,stat=ierr)
        
        end subroutine R2BASIS
c===================================================================
        subroutine R3BASIS
        use COMPARAS, only: NR3, NVB3, R30, R3min, R3max, RMS_CD, IDINIT
        use QUADPOINTS, only : RQ3
        use TRANSMATS, only : TMATR3
        use SMALL_ARRAYS, only : ETR3, VPOT1D, RQ0R3, HAR3, VrefR3
        
        implicit none
        integer ierr
        integer i, in
        
        allocate(TMATR3(NVB3,NVB3),RQ3(NVB3),ETR3(NVB3),VrefR3(NVB3),
     .  stat=ierr)
        if (NVB3.eq.1) then
          RQ3(NVB3)=R30
          ETR3(NVB3)=0.d0
          TMATR3(NVB3,NVB3)=1.d0
          VrefR3(NVB3)=0.d0
          return
        end if
        
        allocate(VPOT1D(NR3),RQ0R3(NR3),stat=ierr)
        call POT1D(NR3,R3min,R3max,RQ0R3,VPOT1D,'R3')
        open(32,file='r3.dat',status='unknown') 
        do in=1,NR3
          write(32,*) RQ0R3(in), VPOT1D(in)*27.2114d0
        enddo
        close(32)

        allocate(HAR3(NR3,NR3),stat=ierr)
        CALL VIB1D(RMS_CD,NR3,RQ0R3,VPOT1D,ETR3,HAR3,
     .  NVB3,NVB3,RQ3,RQ3,TMATR3,TMATR3,VrefR3)
        do in=1,NVB3
          print*, in, RQ3(in), VrefR3(in)*27.2114d0
        enddo
        
        deallocate(VPOT1D,stat=ierr)
        
        end subroutine R3BASIS
c======================================================================
        SUBROUTINE VIB1D(RMS,NR,RQ0,VPOT,EIG0,HAR,NRB,NVB,RQ1,RQ2,
     $  TMAT31,TMAT32,Vref)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION RQ0(1), HAR(NR,1), Vref(1), VPOT(1), WK(NR), 
     $  DTMAT(NRB*NRB), RQ1(1), RQ2(1), EIGj(NR), COFFp(1), EIG0(1)
        REAL*8 TMAT31(NRB,1),TMAT32(NVB,1)
        
        CM=219475.D0

c solve vibrational energy
        CALL VIBSOL(NR,NR,RQ0,RMS,HAR,VPOT,EIGj)
        
        DO I=1,5
          WRITE(*,"(I8,2F14.6)") I,EIGj(I)*CM,(EIGj(I)-EIGj(1))*CM
        END DO
        
        DO I=1,NRB
          EIG0(I)=EIGj(I)
        END DO
        
        CALL PODVR(NR,NRB,HAR,DTMAT,RQ0,RQ1,RMS,EIG0,Vref)
        DO 505 I=1,NRB
        DO 505 J=1,NRB
505     TMAT31(J,I)=DTMAT((I-1)*NRB+J)

c generate PODVR for NVB
        IF(NVB >= NRB) RETURN
        
        CALL PODVR(NR,NVB,HAR,DTMAT,RQ0,RQ2,RMS,EIG0,Vref(NRB+1))
        DO 506 I=1,NVB
        DO 506 J=1,NVB
506     TMAT32(J,I)=DTMAT((I-1)*NVB+J)

        RETURN
        END
c==================================================================
        SUBROUTINE POT1D(NR,R1,R2,RQ0,VPOT,FLAG)
        implicit real*8(a-h,o-z)
        real*8 RQ0(1)
        character*2 FLAG

        DR=(R2-R1)/(NR+1)
        DO I=1,NR
          RQ0(I)=R1+DR*I
        END DO

        if (FLAG.eq.'Z ') call REFPOT1(NR,RQ0,VPOT)
        if (FLAG.eq.'R1') call REFPOT2(NR,RQ0,VPOT)
        if (FLAG.eq.'R2') call REFPOT3(NR,RQ0,VPOT)
        if (FLAG.eq.'R3') call REFPOT4(NR,RQ0,VPOT)

        END SUBROUTINE POT1D
c==================================================================
        SUBROUTINE VIBSOL(NQ1max,NQ1,RQ1,TMS,HAR,VPOT,EIG)

        !  dummy variables
        INTEGER, INTENT(IN) :: NQ1max, NQ1
        REAL(8), INTENT(IN) :: TMS, RQ1(NQ1), VPOT(NQ1)
        REAL(8), INTENT(OUT) :: EIG(NQ1), HAR(NQ1max,NQ1)

        !  subroutine variables
        INTEGER :: I, J, ID, IERR, NDVR
        REAL(8) :: AN, RL, WK1(NQ1), WK2(NQ1)

        !  function declarations
        REAL(8) :: DVRKE1, DVRKE2

c solve eigenvalue problem for one dimensional potential.
        RL = RQ1(NQ1)-RQ1(1) + 2.0_8*(RQ1(2)-RQ1(1))
        NDVR = NQ1+1

        DO I=1,NQ1
          HAR(I,I)=Vpot(I)+DVRKE1(I,NDVR,RL,TMS)
          DO J=1, I-1
            HAR(J,I)=DVRKE2(I,J,NDVR,RL,TMS)
            HAR(I,J)=HAR(J,I)
          END DO
        END DO

        ! eigensolver
        CALL RS(NQ1max,NQ1,HAR,EIG,1,HAR,WK1,WK2,IERR)  ! can be slow
        if (IERR /= 0) then
          print *, 'ERROR  DVRBASIS1 in subroutine VIBSOL'
          RETURN
        end if

        DO I = 1, NQ1
          AN = 0.0_8
          DO J = 1, NQ1
            IF(ABS(HAR(J,I)) > AN) THEN
               AN = ABS(HAR(J,I))
               ID = J
            END IF
          END DO
  
          IF(HAR(ID,I) < 0.0_8) THEN
            DO J=1,NQ1
               HAR(J,I) = -HAR(J,I)
            END DO
          END IF
        END DO  !  end of I loop
        RETURN
        END SUBROUTINE VIBSOL
c------------------------------------------------------------------------
        FUNCTION DVRKE1(I,N,RL,TMS)
        IMPLICIT REAL*8 (A-H,O-Z)
        DATA PI/3.14159265358979D0/

        DVRKE1=PI*PI/4.D0/TMS/RL/RL*((2*N*N+1.D0)/3
     $  -1.D0/DSIN(PI*I/N)**2)
        RETURN
        END
C------------------------------------------------------------------------
        FUNCTION DVRKE2(I,J,N,RL,TMS)
        IMPLICIT REAL*8 (A-H,O-Z)
        DATA PI/3.14159265358979D0/

        DVRKE2=PI*PI/4.D0/TMS/RL/RL*(-1)**(I-J)
     $*(1.D0/DSIN(PI*(I-J)/2.D0/N)**2-1.D0/DSIN(PI*(I+J)/2.D0/N)**2)
        RETURN
        END
c======================================================================
        SUBROUTINE PODVR(NR,NVB,WAVE,TMAT,RQ1,DVR,RMS,EVIB,Vref)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION WAVE(NR,1),TMAT(NVB,1),RQ1(1),DVR(1),
     $  WK1(NR),WK2(NR),Vref(1),EVIB(1)

c construct X-matrix.
        DO 20 I=1,NVB
        DO 20 J=1,I
           TMAT(J,I)=0.D0
        DO K=1,NR
        TMAT(J,I)=TMAT(J,I)+WAVE(K,J)*WAVE(K,I)*RQ1(K)
        END DO
20      TMAT(I,J)=TMAT(J,I)

        CALL RS(NVB,NVB,TMAT,DVR,1,TMAT,WK1,WK2,IERR)
        IF(IERR.NE.0) STOP 'PODVR'

        DO I=1,NVB
           IF(TMAT(1,I).LT.0.D0) THEN
           DO J=1,NVB
           TMAT(J,I)=-TMAT(J,I)
           END DO
           END IF
        END DO

c calculate the reference potentail on DVR point
        DO I=1,NVB
          WK2(I)=0.D0
        END DO

        DO N0=1,NVB

        FAC=DSQRT(2.D0/(NR+1))
        PI=DACOS(-1.D0)
        DO I=1,NR
          WK1(I)=0.D0
          DO J=1,NR
            WK1(I)=WK1(I)+WAVE(J,N0)*FAC*DSIN(I*J*PI/(NR+1))
          END DO
        END DO
        RL=RQ1(NR)-RQ1(1)+2.D0*(RQ1(2)-RQ1(1))
        FAC=DSQRT(2.D0/RL)
        R1=RQ1(1)-(RQ1(2)-RQ1(1))

        DO I=1,NVB
          WF=0.D0
          AKE=0.D0
          DO J=1,NR
            WF=WF+WK1(J)*FAC*DSIN(J*PI*(DVR(I)-R1)/RL)
            AKE=AKE+WK1(J)*FAC*DSIN(J*PI*(DVR(I)-R1)/RL)
     $          *(J*PI/RL)**2/(2.D0*RMS)
          END DO
          IF(ABS(WF).GT.WK2(I)) THEN
            WK2(I)=ABS(WF)
            Vref(I)=EVIB(N0)-AKE/WF
          END IF
        END DO

        END DO

        RETURN
        END
c================================================================
        SUBROUTINE MASS(N,ATOM,AMASS)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION AMASS(N)
        CHARACTER*1 ATOM(N)
        
        DO I=1,N
        
        IF(ATOM(I).EQ.'H' .or. ATOM(I).EQ.'h') THEN
          AMASS(I)=1837.15D0
        ELSE IF(ATOM(I).EQ.'D' .or. ATOM(I).EQ.'d') THEN
          AMASS(I)=3674.3D0
        ELSE IF(ATOM(I).EQ.'C' .or. ATOM(I).EQ.'c') THEN
          AMASS(I)=21870.D0
        ELSE IF(ATOM(I).EQ.'N' .or. ATOM(I).EQ.'n') THEN
          AMASS(I)=25515.D0
        ELSE IF(ATOM(I).EQ.'O' .or. ATOM(I).EQ.'o') THEN
          AMASS(I)=29160.D0
        ELSE IF(ATOM(I).EQ.'F' .or. ATOM(I).EQ.'f') THEN
          AMASS(I)=34627.87D0
        ELSE IF(ATOM(I).EQ.'K' .or. ATOM(I).EQ.'k') THEN
          AMASS(I)=64701.21D0
        ELSE
          PRINT*,'MASS OF', ATOM(I), 'IS NOT DEFINED'
        END IF
        
        END DO
        
        RETURN
        END
c======================================================================
        subroutine ENERGY
        use COMPARAS
        use ENGARRAYS
        use QUADPOINTS
        use SMALL_ARRAYS
        use TRANSMATS
        
        implicit none
        
        integer ierr, i, j, k, l, ms
        
        allocate(BROT_Z(NZVB),stat=ierr)
        do i=1, NZVB
          BROT_Z(i)=(2.d0*RMS_ABCDE*ZQ(i)**2.d0)
        end do
        allocate(BROT_R1(NVB1),stat=ierr)
        do i=1, NVB1
          BROT_R1(i)=(2.d0*RMS_AB*RQ1(i)**2.d0)
        end do
        allocate(BROT_R2(NVB2),stat=ierr)
        do i=1, NVB2
          BROT_R2(i)=(2.d0*RMS_ECD*RQ2(i)**2.d0)
        end do
        allocate(BROT_R3(NVB3),stat=ierr)
        do i=1, NVB3
          BROT_R3(i)=(2.d0*RMS_CD*RQ3(i)**2.d0)
        end do
        
        allocate(EIG4D(NZVB*NVB1*NVB2*NVB3),stat=ierr)
        allocate(Vref4D(NZVB*NVB1*NVB2*NVB3),stat=ierr)
        ms=0
        do l=1, NZVB
        do k=1, NVB1
        do j=1, NVB2
        do i=1, NVB3
          ms=ms+1
          EIG4D(ms)=ETZ(l)+ETR1(k)+ETR2(j)+ETR3(i)
          Vref4D(ms)=VrefZ(l)+VrefR1(k)+VrefR2(j)+VrefR3(i)
        end do
        end do
        end do
        end do
        
        allocate(CPE(NBASJJ),stat=ierr)
        do i=1, NBASJJ
          CPE(i)=Jtot*(Jtot+1.d0)+JJ123(i)*(JJ123(i)+1.d0)
        end do

301     format(I6,5(X,F19.12))
        end subroutine ENERGY
c=======================================================================
        function ITOC(I)
        character*3 ITOC
        character*1 CH

        LOW=MOD(I,10)
        IMID=MOD(I/10,10)
        IHIGH=MOD(I/100,10)

        ITOC(1:1)=CHAR(48+IHIGH)
        ITOC(2:2)=CHAR(48+IMID)
        ITOC(3:3)=CHAR(48+LOW)

        RETURN
        END
c=======================================================================
        SUBROUTINE SINBASIS(NB1,R1,R2,RMS,TMAT,ET)
        IMPLICIT REAL*8 (A-H,O-Z)
        REAL*8 TMAT(NB1,1)
        DIMENSION ET(1),RQ1(1)
        DATA PI/3.14159265358979D0/

        RL=R2-R1
        DO I=1,NB1
           ET(I)=(I*PI/RL)**2/2.D0/RMS
        END DO

        FAC=DSQRT(2.D0/(NB1+1))
        DO I=1,NB1
        DO J=1,NB1
           TMAT(J,I)=FAC*DSIN(I*J*PI/(NB1+1))
        END DO
        END DO

        RETURN
        END 
c=======================================================================
