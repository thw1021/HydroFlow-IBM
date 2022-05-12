C=====================================================================C
C                   inverse of the coefficient matrix                 c
c=====================================================================c 
	SUBROUTINE BRINV(A,N,L)
	DIMENSION A(N,N),IS(N),JS(N)
	DOUBLE PRECISION A,T,D
	L=1
	DO 100 K=1,N
	  D=0.0
	  DO 10 I=K,N
	  DO 10 J=K,N
	    IF (ABS(A(I,J)).GT.D) THEN
	      D=ABS(A(I,J))
	      IS(K)=I
	      JS(K)=J
	    END IF
10	  CONTINUE
	  IF (D+1.0.EQ.1.0) THEN
	    L=0
	    WRITE(*,20)
	    RETURN
	  END IF
20	  FORMAT(1X,'ERR**NOT INV')
	  DO 30 J=1,N
	    T=A(K,J)
	    A(K,J)=A(IS(K),J)
	    A(IS(K),J)=T
30	  CONTINUE
	  DO 40 I=1,N
	    T=A(I,K)
	    A(I,K)=A(I,JS(K))
	    A(I,JS(K))=T
40	  CONTINUE
	  A(K,K)=1/A(K,K)
	  DO 50 J=1,N
	    IF (J.NE.K) THEN
	      A(K,J)=A(K,J)*A(K,K)
	    END IF
50	  CONTINUE
	  DO 70 I=1,N
	    IF (I.NE.K) THEN
	      DO 60 J=1,N
	        IF (J.NE.K) THEN
	          A(I,J)=A(I,J)-A(I,K)*A(K,J)
	        END IF
60	      CONTINUE
	    END IF
70	  CONTINUE
	  DO 80 I=1,N
	    IF (I.NE.K) THEN
	      A(I,K)=-A(I,K)*A(K,K)
	    END IF
80	  CONTINUE
100	CONTINUE
	DO 130 K=N,1,-1
	  DO 110 J=1,N
	    T=A(K,J)
	    A(K,J)=A(JS(K),J)
	    A(JS(K),J)=T
110	  CONTINUE
	  DO 120 I=1,N
	    T=A(I,K)
	    A(I,K)=A(I,IS(K))
	    A(I,IS(K))=T
120	  CONTINUE
130	CONTINUE
	RETURN
	END
