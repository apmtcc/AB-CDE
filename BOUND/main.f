!======================================================================
!       Bound state program in the AB + CDE Jacobi coordinates 
!       Users are required to cite  the following papers：
!       1) J. Chem. Phys. 2013, 138, 124309
!       2) J. Chem. Phys. 2016, 144, 244311 
!       3) Phys. Chem. Chem. Phys. 2021, 23, 22298-22304
!       Corresponding author : Hongwei Song   Email: hwsong@wipm.ac.cn
!       Affiliation: State Key Laboratory of Magnetic Resonance and Atomic 
!        and Molecular Physics, Innovation Academy for Precision Measurement 
!        Science and Technology, Chinese Academy of Sciences, Wuhan, China
!       All Copyrights Reserved by the Original Authors.
!======================================================================
        program bound23
        use COMPARAS, only : IDPOT
        implicit none
        
        print*, 'READINP'
        call READINP

        print*, 'BASCOUP'
        call BASCOUP

        print*, 'PREPOTX' 
        call PREPOTX

        print*, 'ZBASIS'
        call ZBASIS
        
        print*, 'R1BASIS'
        call R1BASIS
        
        print*, 'R2BASIS'
        call R2BASIS
        
        print*, 'R3BASIS'
        call R3BASIS
        
        print*, 'A1BASIS'
        call A1BASIS
        
        print*, 'A2BASIS'
        call A2BASIS
        
        print*, 'A3BASIS'
        call A3BASIS
        
        print*, 'B1BASIS'
        call B1BASIS
        
        print*, 'B2BASIS'
        call B2BASIS
        
        print*, 'CGCOEF'
        call CGCOEF
        
        print*, 'ENERGY'
        call ENERGY

        print*, 'POTENTIAL'
        call POTENTIAL
        
        if (IDPOT.eq.1) goto 1000
        
        print*, 'RVBOUND'
        call RVBOUND

1000    CONTINUE 

        end program
