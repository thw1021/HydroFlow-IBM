C=====================================================================C
C                   Decomposite a matrix                              c
c=====================================================================c 
	SUBROUTINE BMAQR(A,M,N,Q,L)
	DIMENSION A(M,N),Q(M,M)
	DOUBLE PRECISION A,Q,ALPHA,T,U
	IF (M.LT.N) THEN
	  L=0
	  WRITE(*,40)
	  RETURN
	END IF
40	FORMAT(1X,'  FAIL')
	DO 10 I=1,M
	DO 10 J=1,M
	  Q(I,J)=0.0
	  IF (I.EQ.J) Q(I,J)=1.0
10	CONTINUE
	NN=N
	IF (M.EQ.N) NN=M-1
	DO 200 K=1,NN
	  U=0.0
	  DO 20 I=K,M
	    IF (ABS(A(I,K)).GT.U) U=ABS(A(I,K))
20	  CONTINUE
	  ALPHA=0.0
	  DO 30 I=K,M
	    T=A(I,K)/U
	    ALPHA=ALPHA+T*T
30	  CONTINUE
	  IF (A(K,K).GT.0.0) U=-U
	  ALPHA=U*SQRT(ALPHA)
	  IF (ABS(ALPHA)+1.0.EQ.1.0) THEN
	    L=0
	    WRITE(*,40)
	    RETURN
	  END IF
	  U=SQRT(2.0*ALPHA*(ALPHA-A(K,K)))
	  IF (U+1.0.NE.1.0) THEN
	    A(K,K)=(A(K,K)-ALPHA)/U
	    DO 50 I=K+1,M
50	    A(I,K)=A(I,K)/U
	    DO 80 J=1,M
	      T=0.0
	      DO 60 L=K,M
60	      T=T+A(L,K)*Q(L,J)
	      DO 70 I=K,M
70	      Q(I,J)=Q(I,J)-2.0*T*A(I,K)
80	    CONTINUE
	    DO 110 J=K+1,N
	      T=0.0
	      DO 90 L=K,M
90	      T=T+A(L,K)*A(L,J)
	      DO 100 I=K,M
100	      A(I,J)=A(I,J)-2.0*T*A(I,K)
110	    CONTINUE
	    A(K,K)=ALPHA
	    DO 120 I=K+1,M
120	    A(I,K)=0.0
	  END IF
200	CONTINUE
	L=1
	DO 210 I=1,M-1
	DO 210 J=I+1,M
	  T=Q(I,J)
	  Q(I,J)=Q(J,I)
	  Q(J,I)=T
210	CONTINUE
	RETURN
	END

