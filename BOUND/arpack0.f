       subroutine DSARPACK(kchanl,nbasis,nroots,memuse,maxcyc,tolvib,
     .   istat)

       implicit none 
c
c......declare variable 

       integer :: kchanl,nbasis,nroots,memuse,maxcyc,memstatus,ibasis,
     .   istat
       real*8  :: tolvib

c-----------------------------------------------------------------------
c
c     %------------------------------------------------------%
c     | Storage Declarations:                                |
c     |                                                      |
c     | The maximum dimensions for all arrays are            |
c     | set here to accommodate a problem size of            |
c     | N .le. MAXN                                          |
c     |                                                      |
c     | NEV is the number of eigenvalues requested.          |
c     |     See specifications for ARPACK usage below.       |
c     |                                                      |
c     | NCV is the largest number of basis vectors that will |
c     |     be used in the Implicitly Restarted Arnoldi      |
c     |     Process.  Work per major iteration is            |
c     |     proportional to N*NCV*NCV.                       |
c     |                                                      |
c     | You must set:                                        |
c     |                                                      |
c     | MAXN:   Maximum dimension of the A allowed.          |
c     | MAXNEV: Maximum NEV allowed.                         |
c     | MAXNCV: Maximum NCV allowed.                         |
c     %------------------------------------------------------%
c
      integer          maxn, maxnev, maxncv, ldv

c.....begining of the modified lines
c=======================================================================
c     parameter       (maxn=25600, maxnev=10, maxncv=25, ldv=maxn )
c=======================================================================

c
c     %--------------%
c     | Local Arrays |
c     %--------------%
c
c.....begining of the modified lines
c=======================================================================
c     single precision
c    &                 v(ldv,maxncv), workl(maxncv*(maxncv+8)),
c    &                 workd(3*maxn), d(maxncv,2), resid(maxn),
c    &                 ax(maxn)
c     logical          select(maxncv)

      real*8,dimension(:),allocatable   :: workl,workd,resid,ax
      real*8,dimension(:,:),allocatable :: v,d
      logical,dimension(:),allocatable  :: select
c=======================================================================

      integer          iparam(11), ipntr(11)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character        bmat*1, which*2
      integer          ido, n, nev, ncv, lworkl, info, ierr,
     &                 j, nx, ishfts, maxitr, mode1, nconv
      logical          rvec
      real*8           tol, sigma
c
c     %------------%
c     | Parameters |
c     %------------%
c
      real*8           zero
      parameter        (zero = 0.0D+0)
c  
c     %-----------------------------%
c     | BLAS & LAPACK routines used |
c     %-----------------------------%
c
      real*8           dnrm2
      external         dnrm2, daxpy
c
c     %--------------------%
c     | Intrinsic function |
c     %--------------------%
c
      intrinsic        abs
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
c     %-------------------------------------------------%
c     | The following include statement and assignments |
c     | initiate trace output from the internal         |
c     | actions of ARPACK.  See debug.doc in the        |
c     | DOCUMENTS directory for usage.  Initially, the  |
c     | most useful information will be a breakdown of  |
c     | time spent in the various stages of computation |
c     | given by setting msaupd = 1.                    |
c     %-------------------------------------------------%
c
c     include 'debug.h'
c
c\SCCS Information: @(#)
c FILE: debug.h   SID: 2.3   DATE OF SID: 11/16/95   RELEASE: 2
c
c     %---------------------------------%
c     | See debug.doc for documentation |
c     %---------------------------------%
      integer  logfil, ndigit, mgetv0,
     &         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,
     &         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,
     &         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
      common /debug/
     &         logfil, ndigit, mgetv0,
     &         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,
     &         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,
     &         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
c
      ndigit = -3
      logfil = 6
      msgets = 0
      msaitr = 0 
      msapps = 0
      msaupd = 1
      msaup2 = 0
      mseigt = 0
      mseupd = 0

c.....begining of the modified lines
c=======================================================================

C      write(*,*)nbasis,nroots,memuse,tolvib,kchanl,maxcyc
      maxn   = nbasis
      maxnev = nroots
      maxncv = memuse
      ldv    = maxn

      allocate( v(ldv,maxncv), workl(maxncv*(maxncv+8)),
     &          workd(3*maxn), d(maxncv,2), resid(maxn),
     &          ax(maxn), select(maxncv), stat=memstatus)
      if (memstatus.ne.0) stop 'failed to allocate memory in dsarpack'

c=======================================================================
c     
c     %-------------------------------------------------%
c     | The following sets dimensions for this problem. |
c     %-------------------------------------------------%
c
c.....begining of the modified lines
c=======================================================================
c     nx = 10
c     n = nx*nx

      n = nbasis
c=======================================================================
c
c     %-----------------------------------------------%
c     |                                               | 
c     | Specifications for ARPACK usage are set       | 
c     | below:                                        |
c     |                                               |
c     |    1) NEV = 4  asks for 4 eigenvalues to be   |  
c     |       computed.                               | 
c     |                                               |
c     |    2) NCV = 20 sets the length of the Arnoldi |
c     |       factorization                           |
c     |                                               |
c     |    3) This is a standard problem              |
c     |         (indicated by bmat  = 'I')            |
c     |                                               |
c     |    4) Ask for the NEV eigenvalues of          |
c     |       largest magnitude                       |
c     |         (indicated by which = 'LM')           |
c     |       See documentation in DSAUPD for the     |
c     |       other options SM, LA, SA, LI, SI.       | 
c     |                                               |
c     | Note: NEV and NCV must satisfy the following  |
c     | conditions:                                   |
c     |              NEV <= MAXNEV                    |
c     |          NEV + 1 <= NCV <= MAXNCV             |
c     %-----------------------------------------------%
c
      nev   = nroots
      ncv   = memuse

      bmat  = 'I'
      which = 'SM'
c
      if ( n .gt. maxn ) then
         print *, ' ERROR with _SSIMP: N is greater than MAXN '
         write(7,*)' ERROR with _SSIMP: N is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _SSIMP: NEV is greater than MAXNEV '
         write(7,*)' ERROR with _SSIMP: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _SSIMP: NCV is greater than MAXNCV '
         write(*,*)' ERROR with _SSIMP: NCV is greater than MAXNCV '
         go to 9000
      end if
c
c     %-----------------------------------------------------%
c     |                                                     |
c     | Specification of stopping rules and initial         |
c     | conditions before calling DSAUPD                    |
c     |                                                     |
c     | TOL  determines the stopping criterion.             |
c     |                                                     |
c     |      Expect                                         |
c     |           abs(lambdaC - lambdaT) < TOL*abs(lambdaC) |
c     |               computed   true                       |
c     |                                                     |
c     |      If TOL .le. 0,  then TOL <- macheps            |
c     |           (machine precision) is used.              |
c     |                                                     |
c     | IDO  is the REVERSE COMMUNICATION parameter         |
c     |      used to specify actions to be taken on return  |
c     |      from DSAUPD. (See usage below.)                |
c     |                                                     |
c     |      It MUST initially be set to 0 before the first |
c     |      call to DSAUPD.                                | 
c     |                                                     |
c     | INFO on entry specifies starting vector information |
c     |      and on return indicates error codes            |
c     |                                                     |
c     |      Initially, setting INFO=0 indicates that a     | 
c     |      random starting vector is requested to         |
c     |      start the ARNOLDI iteration.  Setting INFO to  |
c     |      a nonzero value on the initial call is used    |
c     |      if you want to specify your own starting       |
c     |      vector (This vector must be placed in RESID.)  | 
c     |                                                     |
c     | The work array WORKL is used in DSAUPD as           | 
c     | workspace.  Its dimension LWORKL is set as          |
c     | illustrated below.                                  |
c     |                                                     |
c     %-----------------------------------------------------%
c
      lworkl = ncv*(ncv+8)

c.....begining of the modified lines
c=======================================================================
c     tol = zero

      tol = tolvib
c=======================================================================

      info = 0
      ido = 0
c
c     %---------------------------------------------------%
c     | Specification of Algorithm Mode:                  |
c     |                                                   |
c     | This program uses the exact shift strategy        |
c     | (indicated by setting PARAM(1) = 1).              |
c     | IPARAM(3) specifies the maximum number of Arnoldi |
c     | iterations allowed.  Mode 1 of DSAUPD is used     |
c     | (IPARAM(7) = 1). All these options can be changed |
c     | by the user. For details see the documentation in |
c     | DSAUPD.                                           |
c     %---------------------------------------------------%
c
      ishfts = 1

c.....begining of the modified lines
c=======================================================================
c     maxitr = 300

      maxitr = maxcyc
c=======================================================================

      mode1 = 1
c
      iparam(1) = ishfts
c                
      iparam(3) = maxitr
c                  
      iparam(7) = mode1
c
c     %------------------------------------------------%
c     | M A I N   L O O P (Reverse communication loop) |
c     %------------------------------------------------%
c
 10   continue
c
c        %---------------------------------------------%
c        | Repeatedly call the routine DSAUPD and take | 
c        | actions indicated by parameter IDO until    |
c        | either convergence is indicated or maxitr   |
c        | has been exceeded.                          |
c        %---------------------------------------------%
c
         call dsaupd ( ido, bmat, n, which, nev, tol, resid, 
     &                 ncv, v, ldv, iparam, ipntr, workd, workl,
     &                 lworkl, info )
c
         if (ido .eq. -1 .or. ido .eq. 1) then
c
c           %--------------------------------------%
c           | Perform matrix vector multiplication |
c           |              y <--- OP*x             |
c           | The user should supply his/her own   |
c           | matrix vector multiplication routine |
c           | here that takes workd(ipntr(1)) as   |
c           | the input, and return the result to  |
c           | workd(ipntr(2)).                     |
c           %--------------------------------------%
c
c.....begining of the modified lines
c=======================================================================
c           call av (nx, workd(ipntr(1)), workd(ipntr(2)))

c            print*, '1. call HMATDRIV'

            call HMATDRIV(workd(ipntr(1)), workd(ipntr(2)))
c=======================================================================
c
c           %-----------------------------------------%
c           | L O O P   B A C K to call DSAUPD again. |
c           %-----------------------------------------%
c
            go to 10
c
         end if 
c
c     %----------------------------------------%
c     | Either we have convergence or there is |
c     | an error.                              |
c     %----------------------------------------%
c
      if ( info .lt. 0 ) then
c
c        %--------------------------%
c        | Error message. Check the |
c        | documentation in DSAUPD. |
c        %--------------------------%
c
c         print *, ' '
c         print *, ' Error with _saupd, info = ', info
c         print *, ' Check documentation in _saupd '
c         print *, ' '

c         write(7,*) ' '
c         write(7,*) ' Error with _saupd, info = ', info
c         write(7,*) ' Check documentation in _saupd '
c         write(7,*) ' '
c
      else 
c
c        %-------------------------------------------%
c        | No fatal errors occurred.                 |
c        | Post-Process using DSEUPD.                |
c        |                                           |
c        | Computed eigenvalues may be extracted.    |  
c        |                                           |
c        | Eigenvectors may be also computed now if  |
c        | desired.  (indicated by rvec = .true.)    | 
c        |                                           |
c        | The routine DSEUPD now called to do this  |
c        | post processing (Other modes may require  |
c        | more complicated post processing than     |
c        | mode1.)                                   |
c        |                                           |
c        %-------------------------------------------%
c           
          rvec = .true.
c
          call dseupd ( rvec, 'All', select, d, v, ldv, sigma, 
     &         bmat, n, which, nev, tol, resid, ncv, v, ldv, 
     &         iparam, ipntr, workd, workl, lworkl, ierr )
c
c         %----------------------------------------------%
c         | Eigenvalues are returned in the first column |
c         | of the two dimensional array D and the       |
c         | corresponding eigenvectors are returned in   |
c         | the first NCONV (=IPARAM(5)) columns of the  |
c         | two dimensional array V if requested.        |
c         | Otherwise, an orthogonal basis for the       |
c         | invariant subspace corresponding to the      |
c         | eigenvalues in D is returned in V.           |
c         %----------------------------------------------%
c
          if ( ierr .ne. 0) then
c
c            %------------------------------------%
c            | Error condition:                   |
c            | Check the documentation of DSEUPD. |
c            %------------------------------------%
c
             print *, ' '
             print *, ' Error with _seupd, info = ', ierr
             print *, ' Check the documentation of _seupd. '
             print *, ' '

             write(7,*) ' '
             write(7,*) ' Error with _seupd, info = ', ierr
             write(7,*) ' Check the documentation of _seupd. '
             write(7,*) ' '
c
          else
c
             nconv =  iparam(5)
             do 20 j=1, nconv
c
c               %---------------------------%
c               | Compute the residual norm |
c               |                           |
c               |   ||  A*x - lambda*x ||   |
c               |                           |
c               | for the NCONV accurately  |
c               | computed eigenvalues and  |
c               | eigenvectors.  (iparam(5) |
c               | indicates how many are    |
c               | accurate to the requested |
c               | tolerance)                |
c               %---------------------------%
c
c.....begining of the modified lines
c=======================================================================
c               call av(nx, v(1,j), ax)

C               print*, '2.call HMATDRIV'
                call HMATDRIV(v(1,j),ax)
c=======================================================================

                call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                d(j,2) = dnrm2(n, ax, 1)
                d(j,2) = d(j,2) / abs(d(j,1))

                if (j.le.istat) then
                  write(1000+j)d(j,1),(v(ibasis,j),ibasis=1,nbasis)
                  Write(2000,9100)j,d(j,1)*27.212,
     &                        (d(j,1)-d(1,1))*27.212 
                end if
      
C                if (j.le.istat) then
C                  do ibasis=1, nbasis
C                    write(60+j,*) v(ibasis,j)
C                  end do
C                end if

c                if (j .eq. istat) then
c                  do ibasis=1, nbasis
c                    write(227,*) abs(v(ibasis,j))**2.d0
c                  end do
c                end if

c                if (j.le.istat) then
c                call plotwv(j-1,v(1,j))
c                end if

c
 20          continue
c
c            %-----------------------------%
c            | Display computed residuals. |
c            %-----------------------------%
c
c             call dmout(6, nconv, 2, d, maxncv, -6,
c     &            'Ritz values and relative residuals')
          end if
c
c         %-------------------------------------------%
c         | Print additional convergence information. |
c         %-------------------------------------------%
c
          if ( info .eq. 1) then
c             print *, ' '
c             print *, ' Maximum number of iterations reached.'
c             print *, ' '

c             write(7,*) ' '
c             write(7,*) ' Maximum number of iterations reached.'
c             write(7,*) ' '
c          else if ( info .eq. 3) then
c             print *, ' ' 
c             print *, ' No shifts could be applied during implicit',
c     &                ' Arnoldi update, try increasing NCV.'
c             print *, ' '

c             write(7,*) ' '
c             write(7,*) ' No shifts could be applied during implicit',
c     &                  ' Arnoldi update, try increasing NCV.'
c             write(7,*) ' '
          end if      
c
c          print *, ' '
c          print *, ' _SSIMP '
c          print *, ' ====== '
c          print *, ' '
c          print *, ' Size of the matrix is ', n
c          print *, ' The number of Ritz values requested is ', nev
c          print *, ' The number of Arnoldi vectors generated',
c     &             ' (NCV) is ', ncv
c          print *, ' What portion of the spectrum: ', which
c          print *, ' The number of converged Ritz values is ', 
c     &               nconv 
c          print *, ' The number of Implicit Arnoldi update',
c     &             ' iterations taken is ', iparam(3)
c          print *, ' The number of OP*x is ', iparam(9)
c          print *, ' The convergence criterion is ', tol
c          print *, ' '

c          write(7,*) ' '
c          write(7,*) ' _SSIMP '
c          write(7,*) ' ====== '
c          write(7,*) ' '
c          write(7,*) ' Size of the matrix is ', n
c          write(7,*) ' The number of Ritz values requested is ', nev
c          write(7,*) ' The number of Arnoldi vectors generated',
c     &               ' (NCV) is ', ncv
c          write(7,*) ' What portion of the spectrum: ', which
c          write(7,*) ' The number of converged Ritz values is ',nconv
c          write(7,*) ' The number of Implicit Arnoldi update',
c     &               ' iterations taken is ', iparam(3)
c          write(7,*) ' The number of OP*x is ', iparam(9)
c          write(7,*) ' The convergence criterion is ', tol
c          write(7,*) ' '
c
      end if
c
c     %---------------------------%
c     | Done with program dssimp. |
c     %---------------------------%
c
c.....release memory 

      deallocate(v,workl,workd,d,resid,ax,select,stat=memstatus)
      if (memstatus.ne.0) stop 'failed to release memory in ssarpack'
 9000 continue
 9100 format(i10,2f16.10)
      return 
      end
