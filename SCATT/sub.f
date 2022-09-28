        subroutine READINP
        use fileIO
        use COMPARAS
        
        implicit none
        
        integer :: i, j, NL, ierr
        
        EV=27.2114d0
        CM=219474.6d0
        PI=DACOS(-1.d0)
        
        read(*,*), IDPOT,VCUT
        read(*,*), KCUT, RCUT
        read(*,*), EIGmin, EIGmax, NE0
        read(*,*), NSTEP, NPRINT
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
        read(*,*), ESKinit, Zinit, DELTA, IDINIT, IREAD  
        read(*,*), istat1, istat2, istat3, istat4
        read(*,*), Z0, R10, R20, R30, A10, A20, A30, B10, B20    
        read(*,*), IDABS, PZABS, CZABS, ZABS, PR1ABS, CR1ABS, R1ABS

        print*, IDPOT,VCUT
        print*, KCUT, RCUT
        print*, EIGmin, EIGmax, NE0
        print*, NSTEP, NPRINT
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
        print*, ESKinit, Zinit, DELTA, IDINIT, IREAD 
        print*, istat1, istat2, istat3, istat4 
        print*, Z0, R10, R20, R30, A10, A20, A30, B10, B20    
        print*, IDABS, PZABS, CZABS, ZABS, PR1ABS, CR1ABS, R1ABS

        call MASS(5,ATOM,AM)

        RMS_ABCDE=(AM(1)+AM(2))*(AM(3)+AM(4)+AM(5))/
     .    (AM(1)+AM(2)+AM(3)+AM(4)+AM(5))
        RMS_AB=AM(1)*AM(2)/(AM(1)+AM(2))
        RMS_CD=AM(4)*AM(5)/(AM(4)+AM(5))
        RMS_ECD=AM(3)*(AM(4)+AM(5))/(AM(3)+AM(4)+AM(5))

        PRINT*,'RMS_ABCDE=',RMS_ABCDE, 'RMS_AB=', RMS_AB
        PRINT *, 'RMS_CD=', RMS_CD, 'RMS_ECD=', RMS_ECD
        
        A10=A10*PI/180.d0
        A20=A20*PI/180.d0
        A30=A30*PI/180.d0
        B10=B10*PI/180.d0
        B20=B20*PI/180.d0

        PRINT *, 'VCUT=', VCUT
        PRINT *, 'KCUT=', KCUT
        PRINT *, 'RCUT=', RCUT

        VCUT=VCUT/CM
        KCUT=KCUT/CM
        RCUT=RCUT/CM

        EIGmin=EIGmin/CM
        EIGmax=EIGmax/CM
! |J1max-J2max| <= J12max <= (J1max+J2max)
! |J2max-J3max| <= L23max <= (J2max+J3max)
        if (JD1.eq.2) J1max=J1max+mod(J1max-JS1,2)
        if (JD3.eq.2) J3max=J3max+mod(J3max-JS3,2)

        KB1max=min(J1max, J2max)
        KM2max=min(J2max, J3max)
        print*, 'KB1max=',KB1max,'KM2max=',KM2max

        NA1node=INT(J1max/JD1)+1
        NA2node=(J2max+1)/JD2
        NA3node=INT(J3max/JD3)+1
        NB1node=(2*KB1max)/NSYM+1
        NB2node=2*KM2max+1

        print*, 'NA1node=', NA1node, 'NA2node=', NA2node,
     .    'NA3node=', NA3node
        print*, 'NB1node=', NB1node, 'NB2node=', NB2node

        nodeax=NA1node*NA2node*NA3node
        nodebx=NB1node*NB2node
        NANODES=nodeax*nodebx
        NRGRIDS=NZVB*NVB1*NVB2*NVB3
        PRINT *, 'nodeax=', nodeax, 'nodebs=', nodebx
        PRINT *, 'NANODES=', NANODES, 'NRGRIDS', NRGRIDS

        return
        end subroutine READINP
c====================================================================
        subroutine ZBASIS
        use COMPARAS, only: NZ, NZVB, Z0, Zmin, Zmax, RMS_ABCDE, istat1,
     .       CM, IDINIT
        use QUADPOINTS, only : ZQ
        use TRANSMATS, only : TMATZ, TMATIZ
        use SMALL_ARRAYS, only : ETZ, VPOT1D, RQ0Z, HARZ, VrefZ,
     .       COFF_Z

        implicit none
        integer ierr
        integer i, j, in
        real*8 DZ, PI, FAC, ETIZ(NZ)

        allocate(TMATZ(NZVB,NZVB),TMATIZ(NZ,NZ),ZQ(NZVB),ETZ(NZVB),
     .   VrefZ(NZVB),COFF_Z(NZVB),stat=ierr)
        if (NZVB.eq.1) then
          ZQ(NZVB)=Z0
          ETZ(NZVB)=0.d0
          TMATZ(NZVB,NZVB)=1.d0
          VrefZ(NZVB)=0.d0
          COFF_Z(NZVB)=1.d0
          return
        end if

        DZ=(Zmax-Zmin)/(NZVB+1)
        DO I=1,NZVB
          ZQ(I)=Zmin+I*DZ
        END DO

        call SINBASIS(NZVB,Zmin,Zmax,RMS_ABCDE,TMATZ,ETZ)
        VrefZ=0.D0

        COFF_Z(1)=1.d0

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

        call SINBASIS(NZ,Zmin,ZQ(NZ+1),RMS_ABCDE,TMATIZ,ETIZ)

        do in=1,NZVB
          print*, in, ZQ(in)
        enddo


        RETURN 

        end subroutine ZBASIS
c====================================================================
        subroutine R1BASIS
        use COMPARAS, only: NR1, NVB1, R10, R1min, R1max, RMS_AB, 
     .       istat2, CM, PI, IDINIT
        use QUADPOINTS, only : RQ1
        use TRANSMATS, only : TMATR1
        use SMALL_ARRAYS, only : ETR1, VPOT1D, RQ0R1, HAR1, VrefR1,
     .       COFF_R1

        implicit none
        integer ierr
        integer i, j, in
        real*8 DR, FAC
 
        allocate(TMATR1(NVB1,NVB1),RQ1(NVB1),ETR1(NVB1),VrefR1(NVB1),
     .   COFF_R1(NVB1),stat=ierr)

        if (NVB1.eq.1) then
          RQ1(NVB1)=R10
          ETR1(NVB1)=0.d0
          TMATR1(NVB1,NVB1)=1.d0
          VrefR1(NVB1)=0.d0
          COFF_R1(NVB1)=1.d0
          return
        end if
        
        allocate(VPOT1D(NR1),RQ0R1(NR1),stat=ierr)
        call POT1D(NR1,R1min,R1max,RQ0R1,VPOT1D,'R1')
        do in=1,NR1
          write(28,*) RQ0R1(in), VPOT1D(in)-minval(VPOT1D)
        enddo         

        allocate(HAR1(NR1,NR1),stat=ierr)
        CALL VIB1D(RMS_AB,NR1,RQ0R1,VPOT1D,ETR1,HAR1,NVB1,NVB1,RQ1,RQ1,
     .     TMATR1,TMATR1,VrefR1,COFF_R1,istat2)

        do in=1,NVB1
          print*, in, RQ1(in), VrefR1(in)*CM
        enddo
        
        deallocate(VPOT1D,stat=ierr)
        
        end subroutine R1BASIS
c====================================================================
        subroutine R2BASIS
        use COMPARAS, only: NR2, NVB2, R20, R2min, R2max, RMS_ECD, 
     .       istat3, CM, IDINIT
        use QUADPOINTS, only : RQ2
        use TRANSMATS, only : TMATR2
        use SMALL_ARRAYS, only : ETR2, VPOT1D, RQ0R2, HAR2, VrefR2,
     .       COFF_R2

        implicit none
        integer i, j, ierr
        integer in
        real*8 DR, FAC, PI
        
        allocate(TMATR2(NVB2,NVB2),RQ2(NVB2),ETR2(NVB2),VrefR2(NVB2),
     .    COFF_R2(NVB2),stat=ierr)
        if (NVB2.eq.1) then
          RQ2(NVB2)=R20
          ETR2(NVB2)=0.d0
          TMATR2(NVB2,NVB2)=1.d0
          VrefR2(NVB2)=0.d0
          COFF_R2(NVB2)=1.d0
          return
        end if

        allocate(VPOT1D(NR2),RQ0R2(NR2),stat=ierr)

        if(IDINIT.eq.1) then
           call POT1D(NR2,R2min,R2max,RQ0R2,VPOT1D,'R2')
        else
           open(30,file='r2.dat',status='old')
           do i=1,NR2
             read(30,*) RQ0R2(i), VPOT1D(i)
           enddo
           close(30)
        end if

        do in=1,NR2
          write(30,*) RQ0R2(in), VPOT1D(in)-MINVAL(VPOT1D)
        enddo
        
        allocate(HAR2(NR2,NR2),stat=ierr)
        CALL VIB1D(RMS_ECD,NR2,RQ0R2,VPOT1D,ETR2,HAR2,NVB2,NVB2,RQ2,RQ2,
     .     TMATR2,TMATR2,VrefR2,COFF_R2,istat3)
         
        do in=1,NVB2
          print*, in, RQ2(in), VrefR2(in)*CM
        enddo
        
        deallocate(VPOT1D,stat=ierr)
        
        end subroutine R2BASIS
c===================================================================
        subroutine R3BASIS
        use COMPARAS, only: NR3, NVB3, R30, R3min, R3max, RMS_CD, 
     .       istat4, CM, IDINIT
        use QUADPOINTS, only : RQ3
        use TRANSMATS, only : TMATR3
        use SMALL_ARRAYS, only : ETR3, VPOT1D, RQ0R3, HAR3, VrefR3,
     .       COFF_R3

        implicit none
        integer i,ierr
        integer in
        
        allocate(TMATR3(NVB3,NVB3),RQ3(NVB3),ETR3(NVB3),VrefR3(NVB3),
     .   COFF_R3(NVB3),stat=ierr)
        if (NVB3.eq.1) then
          RQ3(NVB3)=R30
          ETR3(NVB3)=0.d0
          TMATR3(NVB3,NVB3)=1.d0
          VrefR3(NVB3)=0.d0
          COFF_R3(NVB3)=1.d0
          return
        end if
        
        allocate(VPOT1D(NR3),RQ0R3(NR3),stat=ierr)
        if(IDINIT.eq.1) then
           call POT1D(NR3,R3min,R3max,RQ0R3,VPOT1D,'R3')
        else
           open(32,file='r3.dat',status='old')
           do i=1,NR3
             read(32,*) RQ0R3(i), VPOT1D(i)
           enddo
           close(32)
        end if

        do in=1,NR3
          write(32,*) RQ0R3(in), VPOT1D(in)-MINVAL(VPOT1D)
        enddo

        allocate(HAR3(NR3,NR3),stat=ierr)
        CALL VIB1D(RMS_CD,NR3,RQ0R3,VPOT1D,ETR3,HAR3,NVB3,NVB3,RQ3,RQ3,
     .     TMATR3,TMATR3,VrefR3,COFF_R3,istat4)
        do in=1,NVB3
          print*, in, RQ3(in), VrefR3(in)*CM
        enddo
        
        deallocate(VPOT1D,stat=ierr)

        end subroutine R3BASIS
c======================================================================
        SUBROUTINE VIB1D(RMS,NR,RQ0,VPOT,EIG0,HAR,NRB,NVB,RQ1,RQ2,
     $  TMAT31,TMAT32,Vref,COFFp,ISTAT)
        IMPLICIT REAL*8 (A-H,O-Z)
        INTEGER  ISTAT
        DIMENSION RQ0(1), HAR(NR,1), Vref(1), VPOT(1), WK(NR), 
     $   DTMAT(NRB*NRB), RQ1(1), RQ2(1), EIGj(NR), COFFp(1), EIG0(1)
        REAL*8 TMAT31(NRB,1),TMAT32(NVB,1)
        REAL*8 COFF(NVB,NVB)

        CM=219475.D0

c solve vibrational energy
        CALL VIBSOL(NR,NR,RQ0,RMS,HAR,VPOT,EIGj)
      
        DO I=1,5
          WRITE(*,"(I8,2F14.6)") I,EIGj(I)*CM,(EIGj(I)-EIGj(1))*CM
        END DO

        DO I=1,NRB
          EIG0(I)=EIGj(I)
        END DO

        DO I=1,NVB
        DO Ip=1,NVB
           COFF(Ip,I)=0.D0
           DO K=1,NR
              COFF(Ip,I)=COFF(Ip,I)+HAR(K,Ip)*HAR(K,I)
           END DO
        END DO
        END DO

        DO I=1, NVB
          COFFp(I)=COFF(I,ISTAT)
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
        use COMPARAS, only : NZVB, NVB1, NVB2, NVB3, RMS_ABCDE, RMS_AB, 
     .       RMS_CD,RMS_ECD, Jtot, NBASJJ, CM
        use ENGARRAYS, only : BROT_Z, BROT_R1, BROT_R2, BROT_R3, CPE, 
     .       EIG4D, Vref4D
        use QUADPOINTS, only: ZQ, RQ1, RQ2, RQ3
        use SMALL_ARRAYS, only : ETZ, ETR1, ETR2, ETR3, VrefZ, VrefR1, 
     .       VrefR2, VrefR3
        use TRANSMATS, only  : JJ123
        
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
        function ITOC4(I)
        character*4 ITOC4
        character*1 CH

        LOW=MOD(I,10)
        IMID1=MOD(I/10,10)
        IMID2=MOD(I/100,10)
        IHIGH=MOD(I/1000,10)

        ITOC4(1:1)=CHAR(48+IHIGH)
        ITOC4(2:2)=CHAR(48+IMID2)
        ITOC4(3:3)=CHAR(48+IMID1)
        ITOC4(4:4)=CHAR(48+LOW)

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
        SUBROUTINE  FDAMP
        use QUADPOINTS, only : ZQ, RQ1
        use ENGARRAYS, only : FZABS, FR1ABS
        use COMPARAS, only : NZVB, NVB1, ZABS, CZABS, PZABS, R1ABS, 
     .    CR1ABS, PR1ABS, NZABS, NR1ABS

        implicit none
        integer I, ierr
        
        allocate(FZABS(NZVB),FR1ABS(NVB1), stat=ierr)

        call ABS_POT(PZABS,NZVB,ZQ,CZABS,ZABS,NZABS,FZABS,1)

        call ABS_POT(PR1ABS,NVB1,RQ1,CR1ABS,R1ABS,NR1ABS,FR1ABS,1)

        print*, "The parameters for damping function"
        print*, "=========================="
        PRINT*,'NZABS=',NZABS,'ZABS=',ZQ(NZABS)
        PRINT*,'NR1ABS=',NR1ABS,'R1ABS=',RQ1(NR1ABS)
        print*, "=========================="
        print*

        do I=1, NZVB
          write(212,*) I, FZABS(I)
        end do

        do I=1, NVB1
          write(213,*) I, FR1ABS(I)
        end do

        END SUBROUTINE FDAMP
c=======================================================================
        SUBROUTINE  EBASIS
        use COMPARAS, only:  PI, EIGmin, EIGmax, NE0, EV
        use QUADPOINTS, only : THETAQ, EQ
        use SCALES, only : Hplus, Hminus, Emin

        implicit none
        integer I, ierr 
        real*8 DE
        
        EIGmin=max(Emin,EIGmin)+1.0E-10

        allocate(THETAQ(NE0),EQ(NE0),stat=ierr)

        EIGmax=(EIGmax-Hplus)/Hminus
        EIGmin=(EIGmin-Hplus)/Hminus

        DE=(EIGmax-EIGmin)/(NE0-1)
        do i=1, NE0
            EQ(i)=EIGmin+(i-1)*DE
            THETAQ(i)=DACOS(EQ(i))
            write(38,*)i, EQ(i),THETAQ(i)
         end do

        RETURN 
        END SUBROUTINE  EBASIS
c==================================================================
        SUBROUTINE ABS_POT(PABS,NQ1,RQ1,CABS,RABS,NABS,FABS,ID)
c absorption potentail.

        IMPLICIT REAL*8 (A-H,O-Z)
        REAL*8 FABS(1)
        DIMENSION RQ1(NQ1)

        IF(ID.EQ.1) THEN
        RL=RQ1(NQ1)+RQ1(2)-RQ1(1)
        RABS=RL-RABS
        DO I=1,NQ1
          IF(RQ1(I).GE.RABS) THEN
            FABS(I)=DEXP(-CABS*((RQ1(I)-RABS)/(RL-RABS))**PABS)
          ELSE
            NABS=I
            FABS(I)=1.D0
          END IF
        END DO

        ELSE IF(ID.EQ.-1) THEN
          R1=RQ1(1)-(RQ1(2)-RQ1(1))
          RABS=RABS+R1
        DO I=1,NQ1
          IF(RQ1(I).GE.RABS) THEN
            FABS(I)=1.D0
          ELSE
            FABS(I)=DEXP(-CABS*((RABS-RQ1(I))/(RABS-R1))**PABS)
            NABS=I
          END IF
        END DO

        ELSE
          PRINT*,'WRONG OPTION IN WFABS'
          STOP
        END IF
        
        RETURN
        END
c======================================================
        SUBROUTINE  DNORM(NR,WF,AN)
        IMPLICIT NONE
        REAL*8  WF(NR), AN
        INTEGER NR, I

        AN=0.D0
        DO I=1, NR
           AN=AN+DABS(WF(I))**2 
        END DO

        END SUBROUTINE DNORM
c============================================================
        SUBROUTINE DSUM(NR,A1,WK1,A2,WK2)
        implicit none
        integer  NR,I
        real*8   A1, A2
        real*8   WK1(1),WK2(1)

        DO I=1, NR
          WK2(I)=A1*WK1(I)+A2*WK2(I)
        END DO

        return
        END SUBROUTINE DSUM
c------------------------------------------------------------
        SUBROUTINE DOTP(NR,WK1,WK2,DDOT)
        implicit none
        integer  NR,I
        real*8   ddot
        real*8   WK1(1),WK2(1)

        ddot=0.d0
        DO I=1, NR 
          ddot=ddot+WK1(I)*WK2(I)
        END DO

        return
        END SUBROUTINE DOTP
c------------------------------------------------------------
        SUBROUTINE DAMP(WK1)
        use COMPARAS, only : NZVB,NVB1,NVB2,NVB3,NZABS,NR1ABS,NAp,
     .                       IDABS
        use TRANSMATS, only : TMATZ, TMATR1
        use ENGARRAYS, only : FZABS, FR1ABS

        implicit none
        integer I, J, K, MS
        real*8   WK1(1), CWK(4000*NZVB)

        IF(IDABS.eq.0) RETURN

        call DTRANS(NZVB,NZVB,1,NAp*NVB1*NVB2*NVB3,WK1,CWK,TMATZ,'BG')
        call DTRANS(NVB1,NVB1,NZVB,NAp*NVB2*NVB3,WK1,CWK,TMATR1,'BG')

        DO I=NZABS, NZVB
        DO J=1, NVB1*NVB2*NVB3*NAp
          MS=(I-1)*NVB1*NVB2*NVB3*NAp+J
          WK1(MS)=WK1(MS)*FZABS(I)*FZABS(I)
        END DO
        END DO

        DO I=1, NZVB
        DO J=NR1ABS, NVB1
        DO K=1, NVB2*NVB3*NAp
          MS=(I-1)*NVB1*NVB2*NVB3*NAp+(J-1)*NVB2*NVB3*NAp+K
          WK1(MS)=WK1(MS)*FR1ABS(J)*FR1ABS(J)
        END DO
        END DO
        END DO

        call DTRANS(NVB1,NVB1,NZVB,NAp*NVB2*NVB3,WK1,CWK,TMATR1,'GB')
        call DTRANS(NZVB,NZVB,1,NAp*NVB1*NVB2*NVB3,WK1,CWK,TMATZ,'GB')

        return
        END SUBROUTINE DAMP
