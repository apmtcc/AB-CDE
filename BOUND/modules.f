c--------------------------------------------------------------------------
        module QUADPOINTS
        implicit none
        save

        real(8), allocatable :: ZQ(:), RQ1(:), RQ2(:), RQ3(:)

        real(8), allocatable :: A1node(:), B1node(:), A2node(:),
     .    B2node(:), A3node(:)

        real(8), allocatable :: A1weight(:), B1weight(:), A2weight(:),
     .    B2weight(:), A3weight(:)

        end module QUADPOINTS
c---------------------------------------------------------------------------
        module TRANSMATS
        implicit none
        save

        real(8), allocatable::TMATZ(:,:),TMATR1(:,:),TMATR2(:,:),
     .    TMATR3(:,:)

        real(8), allocatable :: CG1(:,:),CG2(:,:,:,:), 
     .    TMATA1(:,:,:), TMATA2(:,:,:,:),
     .    TMATA3(:,:,:), TMATB1(:,:,:), TMATB2(:,:,:)        
        
        complex(16), allocatable :: CTMATB1(:,:), CTMATB2(:,:)
        
        integer, allocatable :: IDKBZ(:,:), KBZ(:),  
     .    JJ1(:), JJ2(:), JJ3(:), JJ23(:), JJ123(:)

        integer, allocatable :: idxCP1(:,:), idxCP2(:,:)

        real(8), allocatable :: CPVAL(:), CPMAT(:)
 
        end module TRANSMATS
c------------------------------------------------------------------------
        module ENGARRAYS
        implicit none
        save

        real*8, allocatable :: EIG4D(:) 

        real*8, allocatable :: CPE(:)

        real*8, allocatable :: BROT_Z(:),BROT_R1(:),BROT_R2(:),
     .                         BROT_R3(:)
        
        real*8, allocatable :: VPOT9D(:,:), Vref4D(:)
        
        integer, allocatable :: LINK(:,:)
        
        integer, allocatable :: POS_LIST(:,:)

        end module ENGARRAYS
c--------------------------------------------------------------------------
        module COMPARAS
        implicit none
        save

! atom's name and mass
! A:ATOM(1), B:ATOM(2), E:ATOM(3), C:ATOM(4), D:ATOM(5)
        character*1 :: ATOM(5) 
        real*8 AM(5)

! parameters for Z freedom, the distance between MC_AB and MC_ECD
        integer :: NZ, NZVB
        real(8) :: Zmin, Zmax

! parameters for R1 freedom, R1 is the bond distance between A and B
        integer :: NR1, NVB1
        real(8) :: R1min, R1max

! parameters for R2 freedom, R2 is the distance between E and MC_CD
        integer :: NR2, NVB2
        real(8) :: R2min, R2max

! parameters for R3 freedom, R3 is the bond distance between C and D
        integer :: NR3, NVB3
        real(8) :: R3min, R3max

! parameters for angular quadrature points
        integer :: NA1node, NA2node, NA3node, NB1node, NB2node

! parameters for rotation numbers
        integer :: Jtot, KBZmax, parity 
        integer :: J123max, J23max
        integer :: J1max, JS1, JD1
        integer :: J2max, JS2, JD2
        integer :: J3max, JS3, JD3

        integer :: NKBZ
        integer :: NBASJJ, NAp  ! total coupled rotation basis

! parameters for centrifugal potential
        integer :: NDIMCP, NDIMCP2
        integer :: NCP_COUPLED

! parameters for equlibrium conformation  
        real(8) :: Z0,R10,R20,R30,A10,A20,A30,B10,B20   
	real(8) :: vinf

        real(8) :: VCUT

! parameters for potential calculation
        integer :: NCPU, ICPU, IDPOT    

	integer :: KB1max, KM2max

	integer :: nodeax, nodebx, NANODES, NRGRIDS

	integer :: NSYM

! parameter for potential saved
	integer :: NBlock

        real(8) :: RMS_ABCDE, RMS_ECD, RMS_AB, RMS_CD

	integer :: istat, nstat, memuse, maxcyc

	real(8) :: tolVal
!parameters for analysis
        integer :: NZF,NR1F,NR2F,NR3F,NA1F,NA2F,NA3F,
     .      NB1F,NB2F, IDINIT
        
         real(8) :: EV, CM, PI

        end module COMPARAS
c-----------------------------------------------------------------------
	module SMALL_ARRAYS

	implicit none
	save

        real(8), allocatable :: ETZ(:), ETR1(:), ETR2(:), ETR3(:)
        real(8), allocatable :: VPOT1D(:),RQ0Z(:),EIG0(:),HARZ(:,:)
        real(8), allocatable :: RQ0R1(:),RQ0R2(:),RQ0R3(:)
        real(8), allocatable :: HAR1(:,:),HAR2(:,:),HAR3(:,:) 
        real(8), allocatable :: VrefZ(:)
        real(8), allocatable :: VrefR1(:)
        real(8), allocatable :: VrefR2(:)
        real(8), allocatable :: VrefR3(:)
	
        real(8), allocatable :: BROTR3(:),BROTR2(:),BROTR1(:),BROTZ(:)

        
	end module SMALL_ARRAYS
c-----------------------------------------------------------------------
        MODULE fileIO
        IMPLICIT NONE
        SAVE

        INTEGER, PARAMETER :: maxStringLength = 512
        CHARACTER(LEN=maxStringLength), PARAMETER :: 
     .  exePath = './'

      ! "scratchlocal" is a symbolically linked directory in the
      ! home directory (~/) that points to directory 
      ! in the local "/scratch" folder. 

        CHARACTER(LEN=maxStringLength) :: tmpFilename
        INTEGER, PARAMETER :: maxFileSize = 1000000000   ! less than 1 GB

        CONTAINS

        SUBROUTINE catFilename(filename)
        IMPLICIT NONE

      !  dummy variables
        CHARACTER(LEN=*), INTENT(IN) :: filename

        tmpFilename = TRIM(exePath) // filename

        END SUBROUTINE catFilename

        END MODULE fileIO
