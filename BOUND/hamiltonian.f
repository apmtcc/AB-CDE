        subroutine RVBOUND
        use COMPARAS, only : NBASJJ, NRGRIDS, istat
        use ENGARRAYS, only : VPOT9D
        
        implicit none
        
        integer nbasis, memuse, maxcyc, NROOTS, ierr
        integer i
        real*8 tolvib, eigen_eng, an, t
        
        nbasis=NBASJJ*NRGRIDS
        NROOTS=min(istat,10)
        memuse=NROOTS*4
        maxcyc=4000
        tolvib=1.0d-5

        rewind(444)

        call DSARPACK(444,nbasis,NROOTS,memuse,maxcyc,tolvib,
     .    NROOTS)

        end subroutine RVBOUND
c====================================================================
        subroutine HMAT9D(WF1, WA)
        use COMPARAS, only : NZVB,NVB1,NVB2,NVB3,NBASJJ,
     .    NANODES,NRGRIDS,NAp,J1max,J2max,J3max,J23max,
     .    KB1max,KM2max,NA1node,NA2node,NA3node,nodebx,VCUT
        use TRANSMATS, only : TMATZ, TMATR1, TMATR2, TMATR3
        use ENGARRAYS, only : EIG4D
        
        implicit none
        
        integer i, ii, j, k, l, ms, ierr, NLARGE, NLARGE2
        
        integer iz, ir1, ir2, ir3, linkPtr
        
        integer n1, n2, NR
        
        real*8 WF1(1), WA(1)
        
        real*8, allocatable :: CWK(:)
        
        real*8,allocatable :: WK2(:)
        
        real*8,allocatable :: WK3(:)
        
        real*8,allocatable :: WK4(:)
        
        real*8,allocatable :: WK5(:)
        
        real*8,allocatable :: WK6(:)
        
        real*8,allocatable :: WK7(:)

        INTEGER :: M,  MMS, NNS, ibuf
        INTEGER :: iNext
        !  parameters for column decomposition
        INTEGER :: currColLeft, nextColLeft
        
        ms=0
        do j=1, NRGRIDS
        do i=1, NAp
          ms=ms+1
          WA(ms)=WF1(ms)
        end do
        end do
        
        allocate(CWK(MAX(2000*NZVB*NVB1*NVB2*NVB3,NAp)),
     .    stat=ierr)
        
        call DTRANS(NZVB,NZVB,1,NAp*NVB1*NVB2*NVB3,WA,CWK,TMATZ,'BG')
        call DTRANS(NVB1,NVB1,NZVB,NAp*NVB2*NVB3,WA,CWK,TMATR1,'BG')
        call DTRANS(NVB2,NVB2,NZVB*NVB1,NAp*NVB3,WA,CWK,TMATR2,'BG')
        call DTRANS(NVB3,NVB3,NZVB*NVB1*NVB2,NAp,WA,CWK,TMATR3,'BG')
        
        deallocate(CWK)
        
        NLARGE=MAX(NANODES,(KB1max+1)*(J1max+1)*(J2max+1)*(J3max+1)*
     .          (J23Max+1),(2*KM2max+1)*(KB1max+1)*(J1max+1)*
     .    (J2max+1)*(J3max+1),(2*KM2max+1)*(KB1max+1)*NA1node*
     .    (J2max+1)*(J3max+1),(2*KM2max+1)*(KB1max+1)*NA1node*
     .    NA2node*(J3max+1),(2*KM2max+1)*(KB1max+1)*NA1node*
     .    NA2node*NA3node,NBASJJ)

        allocate(WK2(NLARGE),WK3(NLARGE),WK4(NLARGE),stat=ierr)
        
        NLARGE2=MAX(nodebx, (2*KM2max+1)*(KB1max+1))
        
        allocate(WK5(NLARGE2*2),WK6(NLARGE2*2),WK7(NLARGE2*2),stat=ierr)

        call RVMAT0(WA,WK2,WK3,WK4,WK5,WK6,WK7)

        deallocate(WK2,WK3,WK4)
        

        allocate(CWK(MAX(2000*NZVB*NVB1*NVB2*NVB3,NAp)),stat=ierr)

        call DTRANS(NVB3,NVB3,NZVB*NVB1*NVB2,NAp,WA,CWK,TMATR3,'GB')
        call DTRANS(NVB2,NVB2,NZVB*NVB1,NAp*NVB3,WA,CWK,TMATR2,'GB')
        call DTRANS(NVB1,NVB1,NZVB,NAp*NVB2*NVB3,WA,CWK,TMATR1,'GB')
        call DTRANS(NZVB,NZVB,1,NAp*NVB1*NVB2*NVB3,WA,CWK,TMATZ,'GB')
        
        deallocate(CWK)
       
        ms=0
        do j=1, NRGRIDS
        do i=1, NAp
          ms=ms+1
          WA(ms)=WA(ms)+WF1(ms)*EIG4D(j)
        end do
        end do
       
       end subroutine HMAT9D
c=====================================================================
!only Jtot=0, K=0
c=====================================================================
        subroutine RVMAT0(WK,WK2,WK3,WK4,WK5,WK6,WK7)
        use COMPARAS, only : NBASJJ,J1max,J2max,J3max,J23max, 
     .    KB1max,KM2max,NA1node,NA2node,NA3node,VCUT,NVB1,NVB2,
     .    NVB3,NB1node,NB2node, NRGRIDS
        use ENGARRAYS, only : BROT_Z, BROT_R1, BROT_R2, BROT_R3,
     .    CPE, LINK, Vref4D
        use TRANSMATS, only : CG1,CG2,TMATA1,TMATA2,TMATA3,JJ1,
     .    JJ23, JJ2, JJ3, JJ123

        implicit none
        
        integer iZ,iR1,iR2,iR3,j1Val,j2Val,j23Val,j3Val,j123Val,
     .    j23_min,j23_max,kbVal,kbmax,kmVal,ia,ia1,ia2,ia3,ms,ierr

        real*8 WK(NBASJJ,1), WK2(1), WK3(1), WK4(1)

        real*8 WK5(1), WK6(1), WK7(1)
        
        real*8 BROT, AN
        
        integer i, ib, iR, n1, n2, ishift
        
        ishift=(KB1max+1)*(2*KM2max+1)
        
        do iR=1, NRGRIDS

            ir3 = MOD(iR-1, NVB3) + 1
            n1 = (iR-ir3)/NVB3 + 1
            ir2 = MOD(n1-1, NVB2) + 1
            n2 = (n1-ir2)/NVB2 + 1
            ir1 = mod(n2-1, NVB1) + 1
            iz = (n2-ir1)/NVB1 + 1


! multiple centrifugal potential Jtot=0
        do ia=1, NBASJJ
          j1Val=JJ1(ia)
          j23Val=JJ23(ia)
          j3Val=JJ3(ia)
          BROT=CPE(ia)/BROT_Z(iZ)+j1Val*(j1Val+1)/BROT_R1(iR1)
     .    +j23Val*(j23Val+1)/BROT_R2(iR2)+
     .    j3Val*(j3Val+1)/BROT_R3(iR3)
          WK2(ia)=WK(ia,iR)*BROT
        end do
        
        call CG1MAT0(WK(:,iR),WK3,1)
        
        call CG2MAT0(WK3,WK4,1)
        
        call A1MAT0(WK4,WK3,1)
        
        call A2MAT0(WK3,WK4,1)
        
        call A3MAT0(WK4,WK3,1)
        
        ms=(iZ-1)*NVB1*NVB2*NVB3+(iR1-1)*NVB2*NVB3+
     .   (iR2-1)*NVB3+iR3

       do ia=1, NA1node*NA2node*NA3node
       
         if (LINK(ia,ms).eq.0) then
           do ib=1, ishift
              WK3((ia-1)*ishift+ib)=
     .          WK3((ia-1)*ishift+ib)*(VCUT-VRef4D(ms))
           end do
       
         else
            call B1B2MAT0(ms,ia,
     .       WK3((ia-1)*ishift+1:ia*ishift),
     .       WK4(1:2*NB1node*(2*KM2max+1)),
     .       WK5(1:2*NB1node*NB2node),
     .       WK6(1:NB2node*NB1node),
     .       WK7(1:2*(KB1max+1)*(2*KM2max+1)))
         end if
       
       end do
       
       call A3MAT0(WK4,WK3,-1)
       
       call A2MAT0(WK3,WK4,-1)
       
       call A1MAT0(WK4,WK3,-1)
       
       call CG2MAT0(WK3,WK4,-1)
       
       call CG1MAT0(WK(:,iR),WK3,-1)

! add centrifugal interaction
       do ia=1, NBASJJ
         WK(ia,iR)=WK2(ia)+WK(ia,iR)
       end do
       
       end do

       end subroutine RVMAT0
c========================================================================
	subroutine CG1MAT0(WK,WK3,id)
	use COMPARAS, only : NBASJJ,J1max,J2max,J3max,J23max, 
     .    KB1max
	use TRANSMATS, only : CG1, JJ1, JJ23, JJ2, JJ3, JJ123

	implicit none

	real*8 WK(NBASJJ)
	real*8 WK3(0:KB1max,0:J1max,0:J2max,0:J3max,0:J23Max)

	real*8 AN

	integer id,j1Val,j2Val,j23Val,j3Val,j123Val,
     .    j23_min,j23_max,kbVal,kbmax,kmVal,ia

! CG	
	if (id.eq.1) then

	call DZERO1D((KB1max+1)*(J1max+1)*(J2max+1)*(J3max+1)*
     .    (J23max+1), WK3)

	do ia=1, NBASJJ
	  j1Val=JJ1(ia)
          j2val=JJ2(ia)
          j3val=JJ3(ia)
	  j23Val=JJ23(ia)
	  j123Val=JJ123(ia) 
          
	  kbMax=min(j1Val, j2Val)

        do kbVal=0, kbMax

	  WK3(kbVal,j1Val,j2val,j3val,j23Val)=WK3(kbVal,j1Val,
     .      j2val,j3val,j23Val)
     .      +WK(ia)*CG1(kbVal,ia)
 
	end do
	end do

	end if

	if (id.eq.-1) then

	call DZERO1D(NBASJJ, WK)

	do ia=1, NBASJJ
	  j1Val=JJ1(ia)
          j2val=JJ2(ia)
          j3val=JJ3(ia)
	  j23Val=JJ23(ia)
	  j123Val=JJ123(ia) 
       
	  kbMax=min(j1Val, j2Val)

        do kbVal=0, kbMax

	  WK(ia)=WK(ia)+WK3(kbVal,j1Val,j2val,j3val,j23Val)
     . *CG1(kbVal,ia)
 
	end do
	end do

	end if
	end subroutine CG1MAT0
c======================================================================
	subroutine CG2MAT0(WK3,WK4,id)
	use COMPARAS, only : J1max,J2max,J3max,J23max, 
     .    KB1max,KM2max,JD1,JD2,JD3,JS1,JS2,JS3
	use TRANSMATS, only : CG2,JJ1,JJ23, JJ2, JJ3, JJ123

	implicit none

	real*8 WK3(0:KB1max,0:J1max,0:J2max,0:J3max,0:J23Max)

	real*8 WK4(-KM2max:KM2max,0:KB1max,0:J1max,0:J2max,
     .    0:J3max)

	integer id,j1Val,j2Val,j23Val,j3Val,j123Val,
     .    j23_min,j23_max,kbVal,kbmax,kmVal,ia

       real*8 an
! CG2
	if (id.eq.1) then

	call DZERO1D((2*KM2max+1)*(KB1max+1)*(J1max+1)*(J2max+1)*
     .    (J3max+1), WK4)

	do j3Val=JS3, J3max,JD3
	do j2Val=JS2, J2max,JD2
	do j1Val=JS1, J1max,JD1
	  j23_min=abs(j3val-j2val); j23_max=min(j3Val+j2Val,j23max)
	do j23Val=j23_min, j23_max
        do kbVal=0, min(j1Val, j2Val)
	do kmVal=max(-j2Val,-j3Val), min(j2Val,j3Val)
   
	   WK4(kmVal,kbVal,j1Val,j2Val,j3Val)=
     .      WK4(kmVal,kbVal,j1Val,j2Val,j3Val)
     .     +WK3(kbVal,j1Val,j2val,j3val,j23Val)
     .     *CG2(kmVal,j3Val,j23Val,j2Val)
        
        enddo
	end do
	end do
	end do
	end do
	end do

	end if

	if (id.eq.-1) then

	call DZERO1D((KB1max+1)*(J1max+1)*(J2max+1)*(J3max+1)*
     .    (J23max+1),WK3)

	do j3Val=JS3, J3max,JD3
	do j2Val=JS2, J2max,JD2
	do j1Val=JS1, J1max,JD1
	  j23_Min=abs(j3val-j2val); j23_max=min(j3Val+j2Val,j23max)
	do j23Val=j23_min, j23_max
	do kbVal=0, min(j1Val, j2Val)
	do kmVal=max(-j2Val,-j3Val), min(j2Val,j3Val)

	  WK3(kbVal,j1Val,j2val,j3val,j23Val)=
     .      WK3(kbVal,j1Val,j2val,j3val,j23Val)+
     .      WK4(kmVal,kbVal,j1Val,j2Val,j3Val)
     .      *CG2(kmVal,j3Val,j23Val,j2Val)
	
	end do
	end do
	end do
	end do
	end do
	end do

	end if

	end subroutine CG2MAT0
c====================================================================
	subroutine A1MAT0(WK4,WK5,id)
	use COMPARAS, only : J1max,J2max,J3max,J23max, 
     .    KB1max,KM2max,NA1node,JD1,JD2,JD3
	use TRANSMATS, only : TMATA1

	implicit none

	integer id, j1Val,j2Val,j3Val,kbVal,kmVal,ia1

	real*8 WK4(-KM2max:KM2max,0:KB1max,0:J1max,0:J2max,
     .    0:J3max)
	real*8 WK5(-KM2max:KM2max,0:KB1max,NA1node,
     .    0:J2max,0:J3max)
         
       real*8 an

	if (id.eq.1) then

	call DZERO1D((2*KM2max+1)*(KB1max+1)*NA1node*(J2max+1)*
     .    (J3max+1), WK5)

	do ia1=1, NA1node
	do j3Val=0, J3max
	do j2Val=0, J2max
	do j1Val=0, J1max
        do kbVal=0, KB1max
        do kmVal=-KM2max,KM2max 
            WK5(kmVal,kbVal,ia1,j2Val,j3Val)=
     .      WK5(kmVal,kbVal,ia1,j2Val,j3Val)+
     .      WK4(kmVal,kbVal,j1Val,j2Val,j3Val)
     .      *TMATA1(ia1,kbVal,j1Val)
	
	end do
	end do
	end do
	end do
	end do
        end do

	end if

	if (id.eq.-1) then

	call DZERO1D((2*KM2max+1)*(KB1max+1)*(J1max+1)*
     .    (J2max+1)*(J3max+1), WK4)

	do j3Val=0, J3max
	do j2Val=0, J2max
	do j1Val=0, J1max
        do kbVal=0, KB1max
        do kmVal=-KM2max,KM2max
        do ia1=1, NA1node
	  WK4(kmVal,kbVal,j1Val,j2Val,j3Val)=
     .      WK4(kmVal,kbVal,j1Val,j2Val,j3Val)+
     .      WK5(kmVal,kbVal,ia1,j2Val,j3Val)
     .      *TMATA1(ia1,kbVal,j1Val)
	
	end do
	end do
	end do
	end do
	end do
        end do
          
	end if

	end subroutine A1MAT0
c====================================================================
	subroutine A2MAT0(WK5, WK4, id)

	use COMPARAS, only : J1max,J2max,J3max, 
     .    KB1max,KM2max,NA1node,NA2node,JD1,JD2,JD3
	use TRANSMATS, only : TMATA2

	implicit none

	integer id, j1Val,j2Val,j3Val,kbVal,kbmax,kmVal,ia1,ia2

	real*8 AN


        real*8 WK5(-KM2max:KM2max,0:KB1max,NA1node,0:J2max,
     .    0:J3max)

	real*8 WK4(-KM2max:KM2max,0:KB1max,NA1node,NA2node,
     .    0:J3max)

	if (id.eq.1) then

	call DZERO1D((2*KM2max+1)*(KB1max+1)*NA1node*NA2node*
     .    (J3max+1), WK4)

	do j3Val=0, J3max
	do ia2=1, NA2node
	do j2Val=0, J2max
	do ia1=1, NA1node
	do kbVal=0, KB1max
	do kmVal=-KM2max, KM2max

	  WK4(kmVal,kbVal,ia1,ia2,j3Val)=
     .      WK4(kmVal,kbVal,ia1,ia2,j3Val)+
     .      WK5(kmVal,kbVal,ia1,j2Val,j3Val)
     .      *TMATA2(ia2,kmVal,j2Val,-kbVal)
	
	end do
	end do
	end do
	end do
	end do
	end do

	end if

	if (id.eq.-1) then

	call DZERO1D((2*KM2max+1)*(KB1max+1)*NA1node*(J2max+1)*
     .    (J3max+1), WK5)

	do j3Val=0, J3max
	do ia2=1, NA2node
	do j2Val=0, J2max
	do ia1=1, NA1node
	do kbVal=0, KB1max
	do kmVal=-KM2max, KM2max

	  WK5(kmVal,kbVal,ia1,j2Val,j3Val)=
     .      WK5(kmVal,kbVal,ia1,j2Val,j3Val)+
     .      WK4(kmVal,kbVal,ia1,ia2,j3Val)
     .      *TMATA2(ia2,kmVal,j2Val,-kbVal)
	
	end do
	end do
	end do
	end do
	end do
	end do

	end if

	end subroutine A2MAT0
c====================================================================
	subroutine A3MAT0(WK4, WK5, id)

	use COMPARAS, only : J1max,J2max,J3max, 
     .    KB1max,KM2max,NA1node,NA2node,NA3node,JD1,JD2,JD3
	use TRANSMATS, only : TMATA3

	implicit none

	integer id,j1Val,j2Val,j3Val,kbVal,kbmax,kmVal,ia1,ia2,ia3


        real*8 WK4(-KM2max:KM2max,0:KB1max,NA1node,NA2node,
     .    0:J3max)

	real*8 WK5(-KM2max:KM2max,0:KB1max,NA3node,NA2node,
     .    NA1node)

	real*8 AN


	if (id.eq.1) then

	call DZERO1D((2*KM2max+1)*(KB1max+1)*NA1node*NA2node*
     .    NA3node, WK5)

	do ia3=1, NA3node
	do j3Val=0, J3max
	do ia2=1, NA2node
	do ia1=1, NA1node
	do kbVal=0, KB1max
	do kmVal=-KM2max, KM2max

	  WK5(kmVal,kbVal,ia3,ia2,ia1)=
     .      WK5(kmVal,kbVal,ia3,ia2,ia1)+
     .      WK4(kmVal,kbVal,ia1,ia2,j3Val)*TMATA3(ia3,kmVal,j3Val)
	
	end do
	end do
	end do
	end do
	end do
	end do

	end if

	if (id.eq.-1) then

	call DZERO1D((2*KM2max+1)*(KB1max+1)*NA1node*NA2node*
     .    (J3max+1), WK4)

	do ia3=1, NA3node
	do j3Val=0, J3max
	do ia2=1, NA2node
	do ia1=1, NA1node
	do kbVal=0, KB1max
	do kmVal=-KM2max, KM2max

	  WK4(kmVal,kbVal,ia1,ia2,j3Val)=
     .      WK4(kmVal,kbVal,ia1,ia2,j3Val)+
     .      WK5(kmVal,kbVal,ia3,ia2,ia1)*TMATA3(ia3,kmVal,j3Val)
	
	end do
	end do
	end do
	end do
	end do
	end do

	end if

	end subroutine A3MAT0
c====================================================================
	subroutine B1B2MAT0(iR,ia,WK,WK2,WK3,WK4,WK5)
	use COMPARAS, only : KB1max,KM2max,NB1node,NB2node
	use ENGARRAYS, only : VPOT9D, POS_LIST
	use TRANSMATS, only : TMATB1,TMATB2

	implicit none

	integer iR,ia,i,j,k,ms,ierr

	real*8 WK(-kM2max:kM2max,0:kB1max)

	real*8 WK2(-kM2max:kM2max,NB1node,2)

	real*8 WK3(NB2node,NB1node,2)

	real*8 WK4(NB2node,NB1node)

	real*8 WK5(-kM2max:kM2max,0:kB1max,2)

	call DZERO1D((2*KM2max+1)*(KB1max+1)*2, WK2)

	do k=1, NB1node
	do j=0, kB1max
	do i=-kM2max, kM2max
	  WK2(i,k,1)=WK2(i,k,1)+WK(i,j)*TMATB1(k,j,1)
	  WK2(i,k,2)=WK2(i,k,2)+WK(i,j)*TMATB1(k,j,2)
	end do
	end do
	end do

	call DZERO1D(2*NB1node*NB2node,WK3)

	do k=1, NB1node
	do j=1, NB2node 
	do i=-kM2max, kM2max
	  WK3(j,k,1)=WK3(j,k,1)+WK2(i,k,1)*TMATB2(j,i,1)
	  WK3(j,k,2)=WK3(j,k,2)+WK2(i,k,2)*TMATB2(j,i,2)
	end do
	end do
	end do

	call DZERO1D(NB1node*NB2node, WK4)

	ms=0
	do j=1, NB1node
	do i=1, NB2node
	  ms=ms+1
	  WK4(i,j)=WK3(i,j,1)+WK3(i,j,2)
	  WK4(i,j)=WK4(i,j)*VPOT9D(ms,POS_LIST(ia,iR))
	end do
	end do
	
	call DZERO1D((2*KM2max+1)*(KB1max+1)*2, WK2)

	do k=1, NB1node
	do j=1, NB2node 
	do i=-kM2max, kM2max
	  WK2(i,k,1)=WK2(i,k,1)+WK4(j,k)*TMATB2(j,i,1)
	  WK2(i,k,2)=WK2(i,k,2)+WK4(j,k)*TMATB2(j,i,2)
	end do
	end do
	end do

	call DZERO1D((2*KM2max+1)*(KB1max+1)*2, WK5)

	do k=1, NB1node
	do j=0, kB1max
	do i=-kM2max, kM2max
	  WK5(i,j,1)=WK5(i,j,1)+WK2(i,k,1)*TMATB1(k,j,1)
	  WK5(i,j,2)=WK5(i,j,2)+WK2(i,k,2)*TMATB1(k,j,2)
	end do
	end do
	end do

	call DZERO1D((2*KM2max+1)*(KB1max+1), WK)

	do j=0, kB1max
	do i=-kM2max, kM2max
	  WK(i,j)=WK5(i,j,1)+WK5(i,j,2)
	end do
	end do
	
	end subroutine B1B2MAT0
c=======================================================================
	subroutine HMATDRIV(WF1, WF2)
	implicit none
	real*8 WF1(1), WF2(1)

	call HMAT9D(WF1,WF2)

	end subroutine HMATDRIV
c========================================================================
	subroutine DZERO1D(N, WK)
	implicit real*8(a-h,o-z)
	real*8 WK(1)

	do i=1, N
	  WK(i)=0.d0
	end do

	end subroutine DZERO1D
c========================================================================
	subroutine DNORM1D(N, WK, AN)
	implicit real*8(a-h,o-z)
	real*8 WK(1)

	AN=0.d0
	do i=1, N
	  AN=AN+abs(WK(i))**2.d0
	end do

	end subroutine DNORM1D
c========================================================================
