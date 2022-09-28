        subroutine POTENTIAL
        use fileIO
        use COMPARAS, only : NZVB, NVB1, NVB2, NVB3, NA1node, NA2node, 
     .    NA3node, NB1node, NB2node, NRGRIDS, VCUT, 
     .    nodeax, nodebx, IDPOT
        implicit none
        integer ierr

        if (IDPOT.eq.1) then 
          call POT9D
        end if

        if (IDPOT.eq.2) then
          call READ_POT9D
        end if

        end subroutine POTENTIAL
c==============================================================================
        subroutine POT9D 
        use fileIO
        use COMPARAS, only : NZVB, NVB1, NVB2, NVB3, NA1node, NA2node, 
     .    NA3node, NB1node, NB2node, NRGRIDS, VCUT, 
     .    nodeax, nodebx, Nblock
        use ENGARRAYS, only : LINK, Vref4D
        use QUADPOINTS, only : ZQ, RQ1, RQ2, RQ3, A1node, A2node, 
     .    A3node, B1node, B2node

        implicit none

        real*8 Z, R1, R2, R3, A1, A2, A3, B1, B2

        real*8, allocatable :: WK1(:)

        integer i, j, k, l, m, n1, n2, iQA, iQA1, iQA2, iQA3,
     $  iQB, iQB1, iQB2, NR1, NR2, ms, iBlock, jBlock, iflag, 
     $  ierr

        character*3 ITOC


        allocate(LINK(NA1node*NA2node*NA3node,NRGRIDS),stat=ierr)
        allocate(WK1(NB1node*NB2node),stat=ierr)

        print*, 'NRGRIDS=', NRGRIDS
        print*, 'NA1node*NA2node*NA3node=', NA1node*NA2node*NA3node
        print*, 'NB1node*NB2node=', NB1node*NB2node
        
! ICPU is to mimic the real taskNum
! the value of ICPU starts from 0 and ends at NCPU-1
        call catFilename('POT')
        open(100,file=tmpFilename,status='UNKNOWN',
     .    form='UNFORMATTED')

        iBlock=0
        LINK=0
        do m=1, NRGRIDS
 
            i = MOD(m-1, NVB3) + 1
            n1 = (m-i)/NVB3 + 1
            j = MOD(n1-1, NVB2) + 1
            n2 = (n1-j)/NVB2 + 1
            k = mod(n2-1, NVB1) + 1
            l = (n2-k)/NVB1 + 1

            Z=ZQ(l) 
            R1=RQ1(k) 
            R2=RQ2(j)
            R3=RQ3(i)


            do iQA=1, NA1node*NA2node*NA3node 
              iQA3=MOD(iQA-1, NA3node)+1
              n1=(iQA-iQA3)/NA3node+1
              iQA2=MOD(n1-1, NA2node)+1
              iQA1=(n1-iQA2)/NA2node+1
           
              A1=A1node(iQA1); A2=A2node(iQA2); A3=A3node(iQA3)

            do iQB=1, NB1node*NB2node
 
              iQB2=MOD(iQB-1, NB2node)+1
              iQB1=(iQB-iQB2)/NB2node+1
 
              B1=B1node(iQB1); B2=B2node(iQB2)
 
              call potdriv(Z,R1,R2,R3,A1,A2,A3,B1,B2,WK1(iQB))

            end do

            if (minval(WK1).LT.VCUT-0.0001d0) then 
              LINK(iQA,m)=LINK(iQA,m)+1
              iBlock=iBlock+1
              write(100) (Min(WK1(iQB),VCUT),
     .             iQB=1,NB1node*NB2node)
            end if 

            end do

          
        end do

        close(100)
      
        call catFilename('POT_LINK')
        open(100,file=tmpFilename,
     .    status='UNKNOWN',form='UNFORMATTED')
        write(100) iBlock  ! total potential blocks per task

        do m=1, NRGRIDS
             write(100) (LINK(i,m), i=1, nodeax)
        end do
        close(100)

        end subroutine POT9D
c=======================================================================
        subroutine READ_POT9D
        use fileIO
        use COMPARAS, only : NZVB, NVB1, NVB2, NVB3, NRGRIDS, 
     .    NA1node, NA2node, NA3node, NB1node, NB2node, VCUT,
     .    nodeax,nodebx, Vmin, Vmax, PI, CM
        use ENGARRAYS, only : LINK, VPOT9D, Vref4D, POS_LIST
        use QUADPOINTS, only : ZQ, RQ1, RQ2, RQ3, A1node, A2node, 
     .    A3node, B1node, B2node

        implicit none

        real*8 Z,R1,R2,R3,A1,A2,A3,B1,B2,pVmin,pVmax

        integer i,j,k,l,m,n1,n2,iQA,iQA1,iQA2,iQA3,
     $   iQB,iQB1,iQB2,ierr,NR1,NR2,iBlock, Nblock, item

        real*8 WK(NB2node,NB1node), tWK(NB2node,NB1node)
        
        real*8 WK1(NZVB,NVB2), tWK1(NZVB,NVB2)

        character*3 ITOC

        include 'mpif.h'

        print*, 'nodeax=', nodeax
        print*, 'NRGRIDS=', NRGRIDS
        print*, 'NA1node*NA2node*NA3node=', NA1node*NA2node*NA3node
        print*, 'NB1node*NB2node=', NB1node*NB2node
        allocate(LINK(nodeax,NRGRIDS),stat=ierr)
        LINK=0

        VMin=100.d0
        VMax=-100.d0

        call catFilename('POT_LINK')
        open(100,file=tmpFilename,status='UNKNOWN',form='UNFORMATTED')
        read(100) NBlock

        do m=1, NRGRIDS
             read(100) (LINK(i,m), i=1, nodeax)
        end do
        close(100)

        allocate(VPOT9D(nodebx,NBlock),stat=ierr)
        allocate(POS_LIST(nodeax,NRGRIDS),stat=ierr)

        call catFilename('POT')
        open(100,file=tmpFilename,status='UNKNOWN',
     .    form='UNFORMATTED')

        WK=100.d0
        WK1=100.d0

        POS_LIST=0
        iBlock=0
        do m=1, NRGRIDS
            i = MOD(m-1, NVB3) + 1
            n1 = (m-i)/NVB3 + 1
            j = MOD(n1-1, NVB2) + 1
            n2 = (n1-j)/NVB2 + 1
            k = mod(n2-1, NVB1) + 1
            l = (n2-k)/NVB1 + 1

           do iQA=1, nodeax
              if (LINK(iQA,m).eq.1) then
                iBlock=iBlock+1
                POS_LIST(iQA,m)=iBlock
                read(100) (VPOT9D(iQB,iBlock),iQB=1,nodebx)
                do iQB=1, nodebx
                  iQB2=MOD(iQB-1, NB2node)+1
                  iQB1=(iQB-iQB2)/NB2node+1
               
                  WK(iQB2, iQB1)=min(WK(iQB2, iQB1),VPOT9D(iQB,iBlock))
                  WK1(l,j)=min(WK1(l,j),VPOT9D(iQB,iBlock))

                  Vmin=min(pVmin,VPOT9D(iQB,iBLOCK))
                  Vmax=max(pVmax,VPOT9D(iQB,iBLOCK))

                  VPOT9D(iQB,iBLOCK)=VPOT9D(iQB,iBLOCK)-Vref4D(m)

                end do
              end if
            end do
         end do

        close(100)

        print *, 'Vmin=', Vmin*CM 
        print *, 'Vmax=', Vmax*CM 

        call catFilename('POT2D_B1B2')
        open(200,file=tmpFilename,status='UNKNOWN')
        do iQB1=1, NB1node
        do iQB2=1, NB2node
           write(200,*) B1node(iQB1)/PI*180.d0,B2node(iQB2)/PI*180.d0,
     .      WK(iQB2,iQB1)*CM 
        end do
        end do
        close(200)

        call catFilename('POT2D_ZR2')
        open(300,file=tmpFilename,status='UNKNOWN')
        do l=1, NZVB
        do j=1, NVB2
           write(300,*) ZQ(l),RQ2(j),WK1(l,j)*CM 
        end do
        end do
        close(300)

        end subroutine READ_POT9D
c==============================================================================
        subroutine REFPOT1(NR,RQ0,VPOT)
        use COMPARAS,only:Z0,R10,R20,R30,A10,A20,A30,B10,B20

        implicit none

        integer NR, i
        real*8 RQ0(1), VPOT(1)

        do i=1, NR
          call potdriv(RQ0(i),R10,R20,R30,A10,A20,A30,B10,B20,VPOT(i))
        end do

        end subroutine REFPOT1
c==============================================================================
        subroutine REFPOT2(NR,RQ0,VPOT)
        use COMPARAS,only:Z0,R10,R20,R30,A10,A20,A30,B10,B20

        implicit none

        integer NR, i
        real*8 RQ0(1), VPOT(1)

        do i=1, NR
          call potdriv(Z0,RQ0(i),R20,R30,A10,A20,A30,B10,B20,VPOT(i))
        end do

        end subroutine REFPOT2
c==============================================================================
        subroutine REFPOT3(NR,RQ0,VPOT)
        use COMPARAS,only:Z0,R10,R20,R30,A10,A20,A30,B10,B20

        implicit none

        integer NR, i
        real*8 RQ0(1), VPOT(1)

        do i=1, NR
          call potdriv(Z0,R10,RQ0(i),R30,A10,A20,A30,B10,B20,VPOT(i))
        end do

        end subroutine REFPOT3
c==============================================================================
        subroutine REFPOT4(NR,RQ0,VPOT)
        use COMPARAS,only:Z0,R10,R20,R30,A10,A20,A30,B10,B20

        implicit none

        integer NR, i
        real*8 RQ0(1), VPOT(1)

        do i=1, NR
          call potdriv(Z0,R10,R20,RQ0(i),A10,A20,A30,B10,B20,VPOT(i))
        end do

        end subroutine REFPOT4
c==============================================================================
