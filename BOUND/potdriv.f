!
!  ***************************************************************************
!  *                                                                         *
!  *   Program    :  prepotx                                                 *
!  *   Function   :  pre-treatment of potential                              *
!  *                                                                         *
!  ***************************************************************************
       subroutine    PREPOTX
       use COMPARAS, only: Z0,R10,R20,R30,A10,A20,A30,B10,B20,
     .    EV,CM,PI,vinf

       implicit      none
       real*8  v0,va,vb,vc
       real*8  cart(3,5),cart1(3,5)
       integer i, j

       call PREPOT 

       call COORD23(Z0,R10,R20,R30,A10,A20,A30,B10,B20,cart)

!       open(100,file='tt.dat')
!       do i=1, 5
!          read(100,*)  (cart(j,i),j=1,3)
!       end do
         

       do j=1, 3
          cart1(j,5)=cart(j,1)*0.5291772d0
          cart1(j,1)=cart(j,2)*0.5291772d0
          cart1(j,2)=cart(j,3)*0.5291772d0
          cart1(j,3)=cart(j,4)*0.5291772d0
          cart1(j,4)=cart(j,5)*0.5291772d0
       end do

       v0=0.d0
       call fnh3_finn_m(cart1,v0)
       Vinf=v0

       write(20,*) '5'
       write(20,*)
       write(20,*) 'H', (cart1(j,1),j=1, 3)
       write(20,*) 'H', (cart1(j,2),j=1, 3)
       write(20,*) 'H', (cart1(j,3),j=1, 3)
       write(20,*) 'N', (cart1(j,4),j=1, 3)
       write(20,*) 'F', (cart1(j,5),j=1, 3)

       print *, 'Vinf=', Vinf/EV*CM ,'cm-1'
       
c      -----------------------------------------------------------------
       return 
       end
!
!  ***************************************************************************
!  *                                                                         *
!  *   Program    :  potdriv                                                 *
!  *   Function   :  Drive routine of potential energy                       *
!  *                                                                         *
!  ***************************************************************************
       subroutine    POTDRIV(rs,r1,r2,r3,a1,a2,a3,b1,b2,vx) 
       use COMPARAS,   only : EV, CM, PI, vinf
       use QUADPOINTS, only : ZQ
       implicit      none 
       integer       iatom,ixyz,jatom,i,j,k
       real*8        rs,r1,r2,r3,a1,a2,a3,b1,b2,vx,va,vb,vc
       real*8        cart(3,5), cart1(3,5)
       real*8        DIST,drs
       character*3 ITOC
 
       call COORD23(rs,r1,r2,r3,a1,a2,a3,b1,b2,cart)
       if(a1.gt.45.d0/180.d0*dacos(-1.d0)) then 
            vx=0.11d0
            return
       end if
            
       do j=1, 3
          cart1(j,5)=cart(j,1)*0.5291772d0
          cart1(j,1)=cart(j,2)*0.5291772d0
          cart1(j,2)=cart(j,3)*0.5291772d0
          cart1(j,3)=cart(j,4)*0.5291772d0
          cart1(j,4)=cart(j,5)*0.5291772d0
       end do

       call fnh3_finn_m(cart1,vx)

        vx=vx-Vinf

        if(vx.lt.-0.01d0)  then
            vx=3.d0
        end if

       vx=vx/EV

201    format(9(F8.4,2X))

        end 
