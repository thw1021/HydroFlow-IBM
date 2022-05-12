c========================================================================c
c          solving the trangle equations                                 c
c========================================================================c
	SUBROUTINE ATRDE(B,N,M,D,L)
	DIMENSION B(M),D(N)
	DOUBLE PRECISION B,D
	L=1
	IF (M.NE.(3*N-2)) THEN
	  L=-1
	  WRITE(*,10)
	  RETURN
	END IF
10	FORMAT(1X,'  ERR  ')
	DO 20 K=1,N-1
	  J=3*K-2
	  IF (ABS(B(J))+1.0.EQ.1.0) THEN
	    L=0
	    WRITE(*,10)
	    RETURN
	  END IF
	  B(J+1)=B(J+1)/B(J)
	  D(K)=D(K)/B(J)
	  B(J+3)=B(J+3)-B(J+2)*B(J+1)
	  D(K+1)=D(K+1)-B(J+2)*D(K)
20	CONTINUE
	IF (ABS(B(3*N-2))+1.0.EQ.1.0) THEN
	  L=0
	  WRITE(*,10)
	  RETURN
	END IF
	D(N)=D(N)/B(3*N-2)
	DO 30 K=N-1,1,-1
	  D(K)=D(K)-B(3*K-1)*D(K+1)
30	CONTINUE
	RETURN
	END
