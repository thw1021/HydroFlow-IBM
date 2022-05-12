c#=================================================================c
c         smoothing interplating techneque                         c
c==================================================================c
	SUBROUTINE GRADPC(X,Y,N,DY1,DYN,DY)
	DIMENSION X(N),Y(N),DY(N)
      DIMENSION H(N)
	DOUBLE PRECISION X,Y,DY,H,DY1,DYN,
     *                    H0,H1,BETA,ALPHA
	INTEGER N
	DY(1)=0.0
	H(1)=DY1
	H0=X(2)-X(1)
	DO 10 J=2,N-1
	  H1=X(J+1)-X(J)
	  ALPHA=H0/(H0+H1)
	  BETA=(1.0-ALPHA)*(Y(J)-Y(J-1))/H0
	  BETA=3.0*(BETA+ALPHA*(Y(J+1)-Y(J))/H1)
	  DY(J)=-ALPHA/(2.0+(1.0-ALPHA)*DY(J-1))
	  H(J)=(BETA-(1.0-ALPHA)*H(J-1))
	  H(J)=H(J)/(2.0+(1.0-ALPHA)*DY(J-1))
	  H0=H1
10	CONTINUE
	DY(N)=DYN
	DO 20 J=N-1,1,-1
20	DY(J)=DY(J)*DY(J+1)+H(J)
	RETURN
	END