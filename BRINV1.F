C=====================================================================C
C                   inverse of the coefficient matrix                 c
c=====================================================================c 
	SUBROUTINE BRINV1(A,N,L)
	DIMENSION A(L,L),A_TEMP(N,N)
	DOUBLE PRECISION A,A_TEMP
	Do I = 1, N
	Do J = 1, N
	   A_TEMP(I,J) = A(I,J)
	Enddo
	Enddo  
	Call BRINV(A_TEMP,N,M)
	Do I = 1, N
	Do J = 1, N
	   A(I,J) = A_TEMP(I,J)
	Enddo
	Enddo 
	RETURN
	END
