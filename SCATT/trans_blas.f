c=========================================================================
	SUBROUTINE DTRANS(NB,NR,N1,N2,W1,W2,TMAT,OPTION)
c ========================================
c W1 --> W1  with W2 a small working space
c ========================================
	IMPLICIT REAL*8 (A-H,O-Z)
	CHARACTER*2 OPTION
	REAL*8 TMAT(NB,1),alpha,beta,W1(N2,1),W2(2000*NR,1)
	integer omp_get_thread_num

	alpha=1.D0
	beta=0.D0

	MM=2000 
        N2D=N2

        if(NB.eq.1.and.NR.eq.1) return

	IF(OPTION .EQ. 'BG') THEN   ! GRID ---> BASIS
                                    ! sum over the first index.

	    NP=INT(N2D/MM)
	    IF(MOD(N2D,MM).NE.0) NP=NP+1

	    NBR=MAX(NB,NR)

C$OMP PARALLEL SHARED (NBR,NR,NB,N1,N2,TMAT,W1,W2,NP,alpha,beta,MM)
C$OMP&  PRIVATE (I,MS,NS,NSTART,NEND,N22,N,IP,NF,J,K)
C$OMP DO
            DO N=1,NP
	      IP=omp_get_thread_num()+1
              NSTART=(N-1)*MM
              NEND=MIN(N*MM,N2)
              N22=NEND-NSTART

	    DO I=1,N1
	     MS=(I-1)*NBR+1
	     CALL DGEMM('N','N',N22,NR,NB,alpha,W1(NSTART+1,MS),N2,TMAT,
     $       NB,beta,W2(1,IP),MM)
	
	     DO J=1,NR
	       MS=(I-1)*NBR
	       NS=(J-1)*N2
	       NF=(J-1)*MM
	     DO K=1,N22
	       W1(NSTART+K,MS+J)=W2(NF+K,IP)
	     END DO
	     END DO
	     
	    END DO
	    END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL


	ELSE IF (OPTION .EQ. 'GB') THEN    ! BASIS  --->  GRID

            NP=INT(N2D/MM)
            IF(MOD(N2D,MM).NE.0) NP=NP+1


	    NBR=MAX(NB,NR)

C$OMP PARALLEL SHARED (NBR,NR,NB,N1,N2,TMAT,W1,W2,NP,alpha,beta,MM)
C$OMP&  PRIVATE (I,MS,NS,NSTART,NEND,N22,N,IP,NF,J,K)
C$OMP DO
            DO N=1,NP
	      IP=omp_get_thread_num()+1
              NSTART=(N-1)*MM
              NEND=MIN(N*MM,N2)
              N22=NEND-NSTART

            DO I=1,N1
             MS=(I-1)*NBR+1
             CALL DGEMM('N','T',N22,NB,NR,alpha,W1(NSTART+1,MS),N2,TMAT,
     $       NB,beta,W2(1,IP),MM)

	     DO J=1,NB
               MS=(I-1)*NBR
               NS=(J-1)*N2
               NF=(J-1)*MM
             DO K=1,N22
               W1(NSTART+K,MS+J)=W2(NF+K,IP)
             END DO
             END DO

            END DO
	    END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL


	ELSE
          STOP 'OPTION IS WRONG IN DTRANS'
	END IF

	RETURN
	END 
c====================================================================
        function omp_get_max_threads()
        integer omp_get_max_threads

        omp_get_max_threads=1
        return
        end
c=================================================================
        function omp_get_thread_num()
        integer omp_get_thread_num

        omp_get_thread_num=0
        return
        end
c=======================================================================

