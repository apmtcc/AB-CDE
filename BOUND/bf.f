c=================================================================
c J1: the rotation number of AB
c J3: the rotation number of CD
c J23: the orbital rotation of E to CD
c J2: the coupling of J23 and J3
c J123: the coupling of J1 and J2
c=================================================================
	subroutine BASCOUP
	use COMPARAS, only : J1max, JS1, JD1, J2max, JS2, JD2,
     .  J3max, JS3, JD3, J123max, J23max, Jtot, parity, NBASJJ, 
     .  NAp, KBZmax, NKBZ  

	use TRANSMATS, only : KBZ, JJ1, JJ2, JJ3, JJ23, JJ123
	
	implicit none

	integer kzVal, j1Val, j2Val, j3Val, j23Val, j123Val, ierr

        integer in

        NBASJJ=0
	NKBZ=1
	do 10 j1Val=JS1, J1max, JD1
	do 20 j2Val=JS2, J2max, JD2
	do 30 j3Val=JS3, J3max, JD3
	do 50 j123Val=abs(j1Val-j2Val), min(J123max, j1Val+j2Val)
	do 40 j23Val=abs(j2Val-j3Val), min(J23max, j2Val+j3Val)

	  if (parity*(-1)**(Jtot+j123Val+j1Val+j23Val+J3Val).eq.(-1)) 
     .      goto 40  
          NBASJJ=NBASJJ+1
40      end do
50      end do
30      end do
20      end do
10      end do

	print*, 'NBASJJ=', NBASJJ, 'NKBZ=', NKBZ

	allocate(KBZ(NBASJJ), JJ1(NBASJJ), JJ2(NBASJJ), JJ3(NBASJJ),
     .  JJ23(NBASJJ), JJ123(NBASJJ), stat=ierr)

        NBASJJ=0
	do 100 j1Val=JS1, J1max, JD1
	do 200 j2Val=JS2, J2max, JD2
	do 300 j3Val=JS3, J3max, JD3
	do 500 j123Val=abs(j1Val-j2Val), min(J123max, j1Val+j2Val)
	do 400 j23Val=abs(j2Val-j3Val), min(J23max, j2Val+j3Val)
          if (parity*(-1)**(Jtot+j123Val+j1Val+j23Val+J3Val).eq.(-1)) 
     .      goto 400  
          NBASJJ=NBASJJ+1
	  KBZ(NBASJJ)=kzVal
	  JJ1(NBASJJ)=j1Val
	  JJ2(NBASJJ)=j2Val
	  JJ3(NBASJJ)=j3Val
	  JJ23(NBASJJ)=j23Val
	  JJ123(NBASJJ)=j123Val
400      end do
500      end do
300      end do
200      end do
100      end do
         
         NAp=NBASJJ

201      format(I8,5(X,I4))
	end subroutine BASCOUP
c==================================================================
	subroutine CGCOEF
	use TRANSMATS, only : CG1, CG2, JJ1, JJ2, JJ3, JJ23, JJ123
	use COMPARAS, only : NBASJJ, KB1max, KM2max, J1max, J2max, J3max,
     .    J23max
	implicit none
	integer j1, j2, j3, j23, j123, i, kb1, km2, ierr
	real*8 FUNCTCG
        real*8 CG
        integer ibasjj
	allocate(CG1(0:KB1max, NBASJJ),stat=ierr)
	allocate(CG2(-KM2max:KM2max,0:J3max,0:J23max,0:J2max),stat=ierr)

	CG1=0.d0

	do i=1, NBASJJ
          j1=JJ1(i)
	  j2=JJ2(i)
	  j123=JJ123(i)
	do kb1=0, KB1max
	  ! we actually deal with K=0 only
	  CG1(kb1,i)=FUNCTCG(j1,kb1,j2,-kb1,j123)
          if(kb1.gt.0) then
          CG1(kb1,i)=CG1(kb1,i)*dsqrt(2.d0)
c          write(81,801) kb1,i,CG1(kb1,i)
          endif
          
	end do
	end do

	CG2=0.d0
	do j2=0, J2max
	do j23=0, J23max
	do j3=0, J3max
	do km2=-KM2max,KM2max
           !CG2(km2,j23,j2,j3)=FUNCTCG(j3,km2,j2,0,j23)
          CG2(km2,j3,j23,j2)=FUNCTCG(j3,km2,j23,0,j2)
     .    *dsqrt((2.D0*j23+1)/(2.d0*j2+1))
c          write(82,802) j2,j23,j3,km2,CG2(km2,j3,j23,j2)
	end do
	end do
	end do
	end do
801     format(2(I4,X),F16.12)
802     format(4(I4,X),F16.12)
	end subroutine CGCOEF
c==================================================================
c the legendre polynomial for A1
c==================================================================
	subroutine A1BASIS
	use COMPARAS, only : J1max, JD1, NA1node, NSYM, A10, nstat
	use TRANSMATS, only : TMATA1
	use QUADPOINTS, only : A1node, A1weight

	implicit none
	integer i, kind, ierr
	real*8, allocatable :: WK1(:), WK2(:)

        if(nstat.ne.0) JD1=1

	allocate(A1node(NA1node),A1weight(NA1node),stat=ierr)
	if (J1max.eq.0) then
	  A1node(NA1node)=A10
	  return
	end if

	allocate(WK1(NSYM*NA1node),WK2(NSYM*NA1node),stat=ierr)
	kind=1
	call GAUSDRIV(kind,NA1node*JD1,WK1,WK2)

! AB are homonuclear molecule
! we calculate A1 within the range of [0, PI/2]
	do i=1, NA1node
	  A1node(i)=dacos(WK1(i))
	  A1weight(i)=WK2(i)*JD1
	end do

	allocate(TMATA1(NA1node,-J1max:J1max,0:J1max),stat=ierr)
        call LEGEND0(NA1node,NA1node,J1max,A1node,A1weight,TMATA1)

	end subroutine A1BASIS
c==================================================================
	subroutine A2BASIS
	use COMPARAS, only : J2max, NA2node, A20
	use TRANSMATS, only : TMATA2
	use QUADPOINTS, only : A2node, A2weight
	implicit none
	integer i, kind, ierr

	kind=1
	allocate(A2node(NA2node),A2weight(NA2node),stat=ierr)
	if (J2max.eq.0) then
	  A2node(NA2node)=A20
	  return
	end if

	call GAUSDRIV(kind,NA2node,A2node,A2weight)
   
C	call LOBATTO(NA2node,A2node,A2weight)

	do i=1, NA2node
	  A2node(i)=dacos(A2node(i))
	end do


	allocate(TMATA2(NA2node,-J2max:J2max,0:J2max,-J2max:J2max),
     &      stat=ierr)
        call ROTMAT2(NA2node,NA2node,J2max,-J2max,J2max,-J2max,
     &      J2max,A2node,A2weight,TMATA2)

	end subroutine A2BASIS
c==================================================================
	subroutine A3BASIS
	use COMPARAS, only : J3max, JD3, NA3node, NSYM, A30, nstat
	use TRANSMATS, only : TMATA3
	use QUADPOINTS, only : A3node, A3weight

	implicit none
	integer i, kind, ierr
	real*8, allocatable :: WK1(:), WK2(:)

        if(nstat.ne.0) JD3=1

	allocate(A3node(NA3node),A3weight(NA3node),stat=ierr)
	if (J3max.eq.0) then
	  A3node(NA3node)=A30
	  return
	end if

	allocate(WK1(NSYM*NA3node),WK2(NSYM*NA3node),stat=ierr)
	kind=1
	call GAUSDRIV(kind,NA3node*JD3,WK1,WK2)

! CD are homonuclear molecule
! we calculate A3 within the range of [0, PI/2]
	do i=1, NA3node
	  A3node(i)=dacos(WK1(i))
	  A3weight(i)=WK2(i)*JD3
	end do

	allocate(TMATA3(NA3node,-J3max:J3max,0:J3max),stat=ierr)
        call LEGEND0(NA3node,NA3node,J3max,A3node,A3weight,TMATA3)

	end subroutine A3BASIS
c==================================================================
	subroutine B1BASIS
	use COMPARAS, only : NB1node,KB1max,Jtot,NSYM,parity,B10
	use TRANSMATS, only : TMATB1
	use QUADPOINTS, only : B1node, B1weight
	
	implicit none
	integer ierr, kind, i, kb1, ib1 
	real*8 sinkb1, coskb1 
       
        integer ik

	allocate(B1node(NB1node), B1weight(NB1node),stat=ierr)

	if (NB1node.eq.1) then
	  B1node(NB1node)=B10
	  return
	end if
	allocate(TMATB1(NB1node,0:KB1max,2),stat=ierr)
	kind=1
        
        call GAUSCHEB(kind,nsym,NB1node,B1node,B1weight)

        do 3200 kb1=0,kb1max
          do 3000 ib1=1,NB1node
             sinkb1=dsin(kb1*B1node(ib1))*B1weight(ib1)
             coskb1=dcos(kb1*B1node(ib1))*B1weight(ib1)
             if (kb1.ne.0) then 
                coskb1=coskb1*dsqrt(2.0d0)
                sinkb1=sinkb1*dsqrt(2.0d0)
             endif
             if (parity.eq.+1) then 
                tmatb1(ib1,kb1,1)=+coskb1
                tmatb1(ib1,kb1,2)=-sinkb1
             else
                tmatb1(ib1,kb1,1)=+sinkb1
                tmatb1(ib1,kb1,2)=+coskb1
             endif
3000     continue
3200   continue

	end subroutine B1BASIS
c==================================================================
	subroutine B2BASIS
	use COMPARAS, only : NB2node, KM2max, B20
	use TRANSMATS, only : TMATB2
	use QUADPOINTS, only : B2node, B2weight
	
	implicit none
	integer ierr, kind, i, nsymb2, ma2, ib2
	real*8 sinkb2, coskb2

        integer ik

	allocate(B2node(NB2node),B2weight(NB2node),stat=ierr)
	if (NB2node.eq.1) then
	  B2node(NB2node)=B20
	  return
	end if

	kind=1
	nsymb2=1
	call GAUSCHEB(kind,nsymb2,NB2node,B2node,B2weight)

 	
	allocate(TMATB2(NB2node,-KM2max:KM2max,2),stat=ierr) !changed
        do 3200 ma2=-KM2max,KM2max
          do 3000 ib2=1,NB2node
             sinkb2=dsin(ma2*B2node(ib2))*B2weight(ib2)
             coskb2=dcos(ma2*B2node(ib2))*B2weight(ib2)
             tmatb2(ib2,ma2,1)=coskb2
             tmatb2(ib2,ma2,2)=sinkb2
 3000     continue
 3200   continue

	end subroutine B2BASIS
c=============================================================================


