!
!  ***************************************************************************
!  *                                                                         *
!  *   Program    :  coord23                                                 *
!  *   Function   :  AB-CDE system, jacobi to cartesian                      *
!  *                                                                         *
!  ***************************************************************************
       subroutine    COORD23(rs,r1,r2,r3,a1,a2,a3,b1,b2,cart) 
       use COMPARAS,   only : AM
       implicit      none 
       integer       iatom
       real*8        rs,r1,r2,r3,a1,a2,a3,b1,b2,cart(3,5),cent
c      -----------------------------------------------------------------
c      AB molecule
c      -----------------------------------------------------------------
       cart(1,1)=0.0d0
       cart(2,1)=0.0d0
       cart(3,1)=+r1*am(2)/(am(1)+am(2))
       cart(1,2)=0.0d0
       cart(2,2)=0.0d0
       cart(3,2)=-r1*am(1)/(am(1)+am(2))
c      -----------------------------------------------------------------
       call ROTEULAR(2,cart(1,1),b1,a1,0.0d0)
c      -----------------------------------------------------------------
       cart(3,1)=cart(3,1)+rs
       cart(3,2)=cart(3,2)+rs
c      -----------------------------------------------------------------
c      coordinate of DE
c      -----------------------------------------------------------------
       cart(1,5)=0.0d0
       cart(2,5)=0.0d0
       cart(3,5)=0.0d0
       cart(1,4)=0.0d0
       cart(2,4)=0.0d0
       cart(3,4)=r3
c      -----------------------------------------------------------------
c      move COM to origin
c      -----------------------------------------------------------------
       cent=r3*am(4)/(am(4)+am(5))
       cart(3,4)=cart(3,4)-cent
       cart(3,5)=cart(3,5)-cent
c      -----------------------------------------------------------------
c      rotate DE
c      -----------------------------------------------------------------
       call ROTEULAR(2,cart(1,4),0.0d0,a3,0.0d0)
c      -----------------------------------------------------------------
       cart(1,3)=0.0d0
       cart(2,3)=0.0d0
       cart(3,3)=r2
c      -----------------------------------------------------------------
c      move COM to origin
c      -----------------------------------------------------------------
       cent=r2*am(3)/(am(3)+am(4)+am(5))
       cart(3,3)=cart(3,3)-cent
       cart(3,4)=cart(3,4)-cent
       cart(3,5)=cart(3,5)-cent
c      -----------------------------------------------------------------
c      rotate C-(DE)
c      -----------------------------------------------------------------
       call ROTEULAR(3,cart(1,3),0.0d0,a2,b2)
c      -----------------------------------------------------------------
       return
       end
!
!  ***************************************************************************
!  *                                                                         *
!  *   Program    :  roteular                                                *
!  *   Function   :  rotate a vector by eular angle (a,b,c)                  *
!  *                                                                         *
!  ***************************************************************************
       subroutine ROTEULAR(natom,cartes,a,b,c)
       implicit none
       integer :: natom,iatom,ixyz,jxyz
       real*8  :: a,b,c,tmatrot(3,3),cartes(3,natom),
     &            sina,cosa,sinb,cosb,sinc,cosc,temp(3)
c
       sina=dsin(a); cosa=dcos(a)
       sinb=dsin(b); cosb=dcos(b)
       sinc=dsin(c); cosc=dcos(c)
c
       tmatrot(1,1)= cosa*cosb*cosc-sina*sinc
       tmatrot(2,1)= sina*cosb*cosc+cosa*sinc
       tmatrot(3,1)=-sinb*cosc

       tmatrot(1,2)=-cosa*cosb*sinc-sina*cosc
       tmatrot(2,2)=-sina*cosb*sinc+cosa*cosc
       tmatrot(3,2)= sinb*sinc

       tmatrot(1,3)= cosa*sinb
       tmatrot(2,3)= sina*sinb
       tmatrot(3,3)= cosb
c
       do 3000 iatom=1,natom
          do 2000 ixyz=1,3
             temp(ixyz)=0.0d0
             do 1800 jxyz=1,3
                temp(ixyz)=temp(ixyz)
     &                    +tmatrot(ixyz,jxyz)*cartes(jxyz,iatom)
 1800        continue
 2000     continue
          do 2200 ixyz=1,3
             cartes(ixyz,iatom)=temp(ixyz)
 2200     continue
 3000  continue
       return
       end

