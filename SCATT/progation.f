        SUBROUTINE   PREPROP
        use COMPARAS, only : KCUT, VCUT, RCUT, NZVB, NVB1, NVB2, NVB3,
     .       J1max, J23max, J3max, J123max, RMS_CD, RMS_ECD, RMS_AB, 
     .       RMS_ABCDE, Vmin, Vmax, CM 
        use SMALL_ARRAYS, only : ETZ, ETR1, ETR2, ETR3
        use QUADPOINTS, only : ZQ, RQ1, RQ2, RQ3
        use ENGARRAYS, only : Vref4D
        use SCALES, only : Emin, Emax, Hplus, Hminus
        use fileIO

        implicit none
        real*8  Rotmax, Kinetmax,Rotmin, Kinetmin 
        character*4  ITOC4
       
! for Jtot=0    
        Rotmax=j1max*(j1max+1)/(2.d0*RMS_AB*RQ1(1)**2.d0)
     .      +j23max*(j23max+1)/(2.d0*RMS_ECD*RQ2(1)**2.d0)
     .      +j3max*(j3max+1)/(2.d0*RMS_CD*RQ3(1)**2.d0)
     .      +j123max*(j123max+1)/(2.d0*RMS_ABCDE*ZQ(1)**2.d0) 
        Rotmax=min(Rotmax, Rcut)
        Rotmin=0.d0
        print *, 'Rotmax=', Rotmax*CM
        print *, 'Rotmin=', Rotmin*CM

        Kinetmax=ETZ(NZVB)+ETR1(NVB1)+ETR2(NVB2)+ETR3(NVB3)
        Kinetmin=ETZ(1)+ETR1(1)+ETR2(1)+ETR3(1)

        print *, 'Kinetmax=', Kinetmax*CM
        print *, 'Kinetmin=', Kinetmin*CM

        print *, 'Vmax=', Vmax*CM
        print *, 'Vmin=', Vmin*CM

        Emax=(Kinetmax+Rotmax+Vcut)*1.2d0     
        Emin=Kinetmin+Vmin

        if(Emin.gt.0.d0)  then
           Emin=Emin*0.8d0
        else
           Emin=Emin*1.2d0
        end if

        print *, 'Emin=', Emin*CM, 'Emax=', Emax*CM

        Hplus=(Emax+Emin)/2.d0
        Hminus=(Emax-Emin)/2.d0

        print *, 'Hplus=', Hplus*CM, 'Hminus=', Hminus*CM

        return
        END  SUBROUTINE PREPROP 
c======================================================================
        SUBROUTINE  PROPAGATION
        use COMPARAS, only : NZVB,NVB1,NVB2,NVB3,NBASJJ,NSTEP,NSTAT,
     .    NAp, NPRINT, NE0
        use WFARRAYS, only : WFINT, WK0, WK1, WKT, WCK, WKE, WKET
        use QUADPOINTS, only : ZQ, RQ1, RQ2, RQ3, A1node, A2node, 
     .       A3node, B1node, B2node
        use TRANSMATS, only : TMATZ,TMATR1, TMATR2, TMATR3
        use WFARRAYS, only :  WG
        use fileIO

        implicit none
        integer istep, i, j, k, l, m, ms, ierr, Nlarge
        real*8  AN, AN1, AN2, AN3
        real*8  tAN, tAN1, tAN2, tAN3
        real*8  WKE2D(NZVB,NVB1),tWKE2D(NZVB,NVB1) 
        real*8  pWCK(0:NSTEP-1)
        real*8, allocatable :: CWK(:)
        character*3  ITOC
        character*4  ITOC4

        allocate(WK0(NAp*NZVB*NVB1*NVB2*NVB3),
     .    WK1(NAp*NZVB*NVB1*NVB2*NVB3),
     .    WKT(NAp*NZVB*NVB1*NVB2*NVB3),stat=ierr)

        WK0=0.d0
        WK1=0.d0
        WKT=0.d0

        allocate(WCK(0:NSTEP-1),stat=ierr)
        WCK=0.d0

        allocate(WG(NE0),stat=ierr)
        WG=0.d0

!  plot the initial wave function
        NLARGE=MAX(NZVB,NVB1,NVB2,NVB3)*2000
        allocate(CWK(NLARGE*2),stat=ierr)

        call DTRANS(NZVB,NZVB,1,NAp*NVB1*NVB2*NVB3,WFINT,CWK,TMATZ,'BG')
        call DTRANS(NVB1,NVB1,NZVB,NAp*NVB2*NVB3,WFINT,CWK,TMATR1,'BG')
        call DTRANS(NVB2,NVB2,NZVB*NVB1,NAp*NVB3,WFINT,CWK,TMATR2,'BG')
        call DTRANS(NVB3,NVB3,NZVB*NVB1*NVB2,NAp,WFINT,CWK,TMATR3,'BG')
        deallocate(CWK)
      
        WKE2D=0.d0
        ms=0
        do J=1, NZVB
        do I=1, NVB1
        do K=1, NVB2
        do L=1, NVB3*NAp 
           ms=ms+1
           WKE2D(J,I)=WKE2D(J,I)+DABS(WFINT(ms))**2 
        end do
        end do
        end do
        end do

        call  catFilename('WKEN2D'//'T'//ITOC4(0))
        open(999,file=tmpFilename,status='UNKNOWN')
        do J=1, NZVB
        do I=1, NVB1
          write(999,*) ZQ(J), RQ1(I), WKE2D(J,I)
        end do
        end do
        close(999)

        allocate(CWK(NLARGE*2),stat=ierr)
        call DTRANS(NVB3,NVB3,NZVB*NVB1*NVB2,NAp,WFINT,CWK,TMATR3,'GB')
        call DTRANS(NVB2,NVB2,NZVB*NVB1,NAp*NVB3,WFINT,CWK,TMATR2,'GB')
        call DTRANS(NVB1,NVB1,NZVB,NAp*NVB2*NVB3,WFINT,CWK,TMATR1,'GB')
        call DTRANS(NZVB,NZVB,1,NAp*NVB1*NVB2*NVB3,WFINT,CWK,TMATZ,'GB')
        deallocate(CWK)

!   the first step, prepare Psi0 and Psi1

        print *, 'Prepare the first two  Chebyshev polynomial'

        call DCOPY(NAp*NZVB*NVB1*NVB2*NVB3,WFINT,1,WK0,1)          
        IF(NSTAT.gt.0) call PROJECTION(WK0,0)
        call DCOPY(NAp*NZVB*NVB1*NVB2*NVB3,WFINT,1,WK1,1)          
        IF(NSTAT.gt.0) call PROJECTION(WK1,1)
        call HMAT(WK1,WK0) 
     
        call DOTP(NAp*NZVB*NVB1*NVB2*NVB3,WK0,WFINT,WCK(0))
        call DOTP(NAp*NZVB*NVB1*NVB2*NVB3,WK1,WFINT,WCK(1))

        call DNORM(NAp*NZVB*NVB1*NVB2*NVB3,WK0,AN)
        print *, 'Norm of WK           0',  AN

        call DNORM(NAp*NZVB*NVB1*NVB2*NVB3,WK1,AN)
        print *, 'Norm of WK           1',  tAN
!  start to do recursion

        print *, 'Start to recursion'

        DO 100 ISTEP=2, NSTEP-1
          WKT=0.d0
          call DCOPY(NAp*NZVB*NVB1*NVB2*NVB3,WK1,1,WKT,1)          
          call HMAT(WK1,WKT) 
          call DAMP(WK0)
          call DSUM(NAp*NZVB*NVB1*NVB2*NVB3,-1.d0,WK0,2.d0,WK1) 
          call DOTP(NAp*NZVB*NVB1*NVB2*NVB3,WK1,WFINT,WCK(ISTEP))
          call DNORM(NAp*NZVB*NVB1*NVB2*NVB3,WK1,AN)

          print *, 'Norm of WK', istep,  AN

          call DCOPY(NAp*NZVB*NVB1*NVB2*NVB3,WKT,1,WK0,1)   

          if(MOD(ISTEP,NPRINT).eq.0) then
            call ANALYSIS(ISTEP) 
            call WRTOUT(ISTEP)
          end if

100    CONTINUE

       deallocate(WK0,WK1,WKT,WFINT)

       call  catFilename('WCK.dat')
       open(111,file=tmpFilename,status='UNKNOWN')
       do J=0, NSTEP-1 
          write(111,*) J, WCK(J)
       end do
       close(111)

        END  SUBROUTINE   PROPAGATION
C======================================================================
        SUBROUTINE WRTOUT(NSTEP)
        use WFARRAYS, only : WCK
        use fileIO

        implicit none
        integer NSTEP,J

        call  catFilename('WCK.dat')
        open(111,file=tmpFilename,status='UNKNOWN')
        do J=0, NSTEP 
           write(111,*) J, WCK(J)
        end do
        close(111)
     
        return
        END SUBROUTINE WRTOUT
