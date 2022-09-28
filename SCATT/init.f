       SUBROUTINE WFINIT
       use WFARRAYS, only : WFINT
       use COMPARAS, only : NZ, NZVB, NVB1, NVB2, NVB3, NBASJJ, NAP, 
     .    NKBZinit, KBZmin, KBZmax, J123init, IDINIT, IREAD, EV
       use TRANSMATS, only : TMATZ, TMATIZ
       use  fileIO

       implicit none

       integer maxnodes, ierr, i, j, k, l, NLARGE
       integer j12Val, j2Val, ik2, kb2Val, kbzVal, ikbz, ibasJJ
       real*8   AN, tAN
       character*3 ITOC
       real*8, allocatable :: W1(:),WTEM(:)
        
       allocate(WFINT(NAp*NZVB*NVB1*NVB2*NVB3), stat=ierr)
       WFINT=0.d0

       print *,'Read in the initial wave packet'

       call catFilename1('WFN'//itoc(IREAD))
       OPEN(300,FILE=tmpFilename1,status='unknown',
     .          form='unformatted')
       read(300)
       do i=1, NZ*NVB1*NVB2*NVB3*NAp
         read(300) WFINT(i)
       end do
       close(300)

       NLARGE=2000*MAX(NZ,NZVB)*2
       allocate(WTEM(NLARGE),stat=ierr)
       call DTRANS(NZ,NZ,1,NAp*NVB1*NVB2*NVB3,WFINT,WTEM,
     .    TMATIZ,'BG')
       call DTRANS(NZVB,NZVB,1,NAp*NVB1*NVB2*NVB3,WFINT,WTEM,
     .    TMATZ,'GB')
       deallocate(WTEM)

       call DNORM(NZVB*NVB1*NVB2*NVB3*NAp,WFINT,AN)

       print*,'norm of initial wavefunction WA is:',AN
      
       END
c=============================================================
