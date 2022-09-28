c-----------------------------------------------------
c this subroutine perform cosine Fourier transformation
c of the auto-correlation function in the order domain
c WCK       
c-----------------------------------------------------
        SUBROUTINE  ANALYSIS(NSTEP)
        use COMPARAS, only : PI, EV, NE0, CM
        use WFARRAYS, only : WCK, WG
        use QUADPOINTS, only : THETAQ, EQ
        use SCALES, only : Hplus, Hminus
        use fileIO

        implicit none
        integer k, j, ierr, NSTEP
        real*8  fac

        WG=0.d0
        do j=1, NE0
          WG(j)=0.d0
          do k=0, NSTEP
            fac=2.d0*dcos(k*THETAQ(j))/(PI*dsin(THETAQ(j)))
            if(k.eq.0) fac=fac/2.d0 
            WG(j)=WG(j)+fac*WCK(k) 
           end do
        end do

        call  catFilename('WG.dat')
        open(222,file=tmpFilename,status='UNKNOWN')
        do J=1, NE0
           write(222,*) (Hplus+EQ(J)*Hminus)*CM, WG(J)
        end do
        close(222)

        END SUBROUTINE ANALYSIS
c-----------------------------------------------------------------
        SUBROUTINE PROJECTION(WK,K)
        use COMPARAS, only : NZVB, NVB1, NVB2, NVB3, NBASJJ,NSTEP, 
     .       NAp, PI, EV, NE0, NSTAT, NESTAT
        use WFARRAYS, only : WKE 
        use QUADPOINTS, only : THETAQ, EQ

        implicit none
        real*8  WK(1) 
        real*8  fac
        integer K,j,i,ms 

        do i=1, NSTAT
          ms=NESTAT(i)
          fac=2.d0*dcos(k*THETAQ(ms))/(PI*dsin(THETAQ(ms)))
          do j=1, NAp*NZVB*NVB1*NVB2*NVB3
            if(k.eq.0) fac=fac/2.d0 
            WKE(j,i)=WKE(j,i)+fac*WK(j) 
          end do

        end do
        END SUBROUTINE PROJECTION
