C############################################################################
C                     BI-CGSTAB method for equations solving                #
c                            VERSION 1.0 (28/07/2009)                       #
C                            AUTHORIZED BY ZHANG JINGXIN                    #
C                            SHANGHAI JIAO TONG UNIVERSITY                  #
C                                 SHANGHAI, CHINA                           #
c                                                                           #
c############################################################################
	Subroutine SOLVEELFSORCG
      Include './Include/OCERM_INF'
	Parameter(EPSON=1.E-20,EPSI=1.E-10,INUM=200,OMIGA = 1.25)

	Common/ELFBLK/CS(IJM,IPOLYGEN),CB(IJM),CP(IJM),X(IJM)

	Common/ELFPRE/ASL(IJM,IPOLYGEN),APL(IJM),
     &	          ASR(IJM,IPOLYGEN),APR(IJM),BBP(IJM),XP(IJM)

	Dimension R1(IJM),U1(IJM),V1(IJM),P1(IJM)
	Dimension TRACE(IJM)
c-----  initialization of the Matrix
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,K)
!$OMP DO
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
		    Do J = 1, CELL_POLYGEN(I)
	         If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
	            CS(I,J) = CS(I,J) / 
     &				   Sqrt(CP(I)) / Sqrt(CP(CELL_SIDE(I,J,2)))
	         Endif
		    Enddo
	      CB(I) = CB(I) / Sqrt(CP(I))
	   Endif
	Enddo
!$OMP END DO
!$OMP DO
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      TRACE(I) = CP(I)
	      CP(I) = 1.0
	   Endif
	Enddo
!$OMP END DO
!$OMP END PARALLEL
C===========================================================================C
C                optimazing the Matrix                                      c
c===========================================================================c
	Call PREMATRIXELFSORCG(OMIGA)                 ! Preconditioning Matrix
c-----  INITIAL VALUES	
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,K)
!$OMP DO 
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
C	      X(I) = 0.0
	      X(I)  = ELF(I)
c	      X(I) = ELF(I) - EL(I)
	   Endif
	Enddo
!$OMP END DO
!$OMP BARRIER
!$OMP DO	
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      R1(I) = 0.0
	      Do J = 1, CELL_POLYGEN(I)
	         If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
	            R1(I) = R1(I) + CS(I,J) * X(CELL_SIDE(I,J,2))
	         Endif
            Enddo
	      R1(I) = CB(I) + R1(I) - CP(I) * X(I)
	   Endif
	Enddo
!$OMP END DO
!$OMP DO
	Do I = 1,IJM
	   If(CCM(I) .EQ. 1.0) Then
	      U1(I) = 0.0
	      V1(I) = 0.0
	      BBP(I) = R1(I)
	   Endif
	Enddo
!$OMP END DO
!$OMP END PARALLEL
	Call SOLVETRIELFSORCG
	ROU = 0.0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I) REDUCTION(+:ROU)		
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      U1(I) = XP(I)
	      P1(I) = U1(I)
	      ROU = ROU + R1(I) * U1(I)
	   Endif
	Enddo
!$OMP END PARALLEL DO	
c	ROUINI = ROU
C	CALL CPU_TIME(BEGIN1)
C---------------------------------------------------------------------------C
      KNUM = 0
	ALPHA = 1.0
	AR = 1.0E6
CCC	ROU = 1.0
 10   Continue
      Do While (AR .GT. EPSI .AND. KNUM .LE. INUM) 
      KNUM = KNUM+1
C-----   AP = (A,P)
!$OMP PARALLEL DO DEFAULT(SHARED)
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      V1(I) = 0.0
	      Do J = 1, CELL_POLYGEN(I)
	         If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
	            V1(I) = V1(I) + CS(I,J) * P1(CELL_SIDE(I,J,2))
	         Endif
	      Enddo
	      V1(I) = -V1(I) + CP(I) * P1(I)
	   Endif
	Enddo
!$OMP END PARALLEL DO
C-----    BV=(P,AP)
	BV = 0.0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I) REDUCTION(+:BV)	
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      BV = BV + P1(I) * V1(I)
	   Endif
	Enddo
!$OMP END PARALLEL DO
C-----    ALPHA
	ALPHA = ROU / (BV + EPSON)
C-----    NEW X
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I)
	Do I = 1, IJM
     	   If(CCM(I) .EQ. 1.0) Then
	      X(I) = X(I) + ALPHA * P1(I)
	   Endif
	Enddo
!$OMP END PARALLEL DO
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I)
	Do I = 1, IJM
     	   If(CCM(I) .EQ. 1.0) Then
	      R1(I) = R1(I) - ALPHA * V1(I)
	      BBP(I) = R1(I)
	   Endif
	Enddo
!$OMP END PARALLEL DO
C-----    PRE-MATRIX
	Call SOLVETRIELFSORCG
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I)
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      U1(I) = XP(I)
	   Endif
	Enddo
!$OMP END PARALLEL DO
	ROUT = 0.0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I) REDUCTION(+:ROUT)
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      ROUT = ROUT + R1(I) * U1(I)
	   Endif
	Enddo
!$OMP END PARALLEL DO
	BETA = ROUT / (ROU + EPSON)
	ROU = ROUT
C-----   NEW P
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I)
	Do I = 1, IJM
     	   If(CCM(I) .EQ. 1.0) Then
	      P1(I) = U1(I) + BETA * P1(I)
	   Endif
	Enddo
!$OMP END PARALLEL DO
	AR = 0.0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I) REDUCTION(+:AR)
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      AR = AR + R1(I) * R1(I)
	   Endif
	Enddo
!$OMP END PARALLEL DO
	AR = Sqrt(AR)
	Enddo
c      If(AR .LT. EPSI .OR. KNUM .GE. INUM) Goto 1000
c	Goto 10
C---------------------------------------------------------------------------C
C1000	Continue
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I) 
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      X(I) = X(I) / Sqrt(TRACE(I))
	   Endif
	Enddo
!$OMP END PARALLEL DO
	Return
	End

c===========================================================================c
c         preconditioning the matrix for CG method                          c
c===========================================================================c
	Subroutine PREMATRIXELFSORCG(OMIGA)
      Include './Include/OCERM_INF'

	Common/ELFBLK/CS(IJM,IPOLYGEN),CB(IJM),CP(IJM),X(IJM)

	Common/ELFPRE/ASL(IJM,IPOLYGEN),APL(IJM),
     &	          ASR(IJM,IPOLYGEN),APR(IJM),BBP(IJM),XP(IJM)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J)
	Do I = 1, IJM
	   Do J = 1, CELL_POLYGEN(I)
	      ASL(I,J) = 0.0
	      ASR(I,J) = 0.0
	   Enddo
	   APL(I) = 0.0
	   APR(I) = 0.0
	Enddo
!$OMP END PARALLEL DO
C-----  L and LT
	GSCOE = Sqrt(OMIGA * (2. - OMIGA))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,K)
!$OMP DO
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      APL(I) = Sqrt(CP(I)) / GSCOE
	      APR(I) = Sqrt(CP(I)) / GSCOE
	      Do J = 1, CELL_POLYGEN(I)
	         If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
	            INDEXCELL = CELL_SIDE(I,J,2)
	            If(INDEXCELL .LT. I) Then
	               ASL(I,J) = CS(I,J) * OMIGA / GSCOE / 
     &                      Sqrt(CP(INDEXCELL))
	            Endif
C	               Do JJ = 1, CELL_POLYGEN(INDEXCELL)
C	                  If(CFM(CELL_SIDE(INDEXCELL,JJ,1)) .EQ. 1.0) Then
C	                     INDEXCELL2 = CELL_SIDE(INDEXCELL,JJ,2)
C						 If(INDEXCELL2 .EQ. I) Then
C	                        ASR(INDEXCELL,K,JJ) = 
C     &							AS(INDEXCELL,K,JJ) * OMIGA / GSCOE /
C     &                            Sqrt(AP(INDEXCELL,K))
C	                     Endif
C	                  Endif
C	               Enddo
	         Endif
	      Enddo
	   Endif
	Enddo
!$OMP END DO NOWAIT
!$OMP DO
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do J = 1, CELL_POLYGEN(I)
	         If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
	            INDEXCELL = CELL_SIDE(I,J,2)
	            If(INDEXCELL .GT. I) Then
	               ASR(I,J) = CS(I,J) * OMIGA / GSCOE / Sqrt(CP(I))
	            Endif
	         Endif
	      Enddo
	   Endif
	Enddo
!$OMP END DO NOWAIT
!$OMP END PARALLEL
	Return
	End

C############################################################################
C                     method for tri-dignal equations solving               #
c                            VERSION 1.0 (21/05/2012)                       #
C                            AUTHORIZED BY ZHANG JINGXIN                    #
C                            SHANGHAI JIAO TONG UNIVERSITY                  #
C                                 SHANGHAI, CHINA                           #
c                                                                           #
c############################################################################
	Subroutine SOLVETRIELFSORCG
      Include './Include/OCERM_INF'
	Parameter(EPSON=1.E-12,EPSI=1.E-5,INUM=200)

	Common/ELFPRE/ASL(IJM,IPOLYGEN),APL(IJM),
     &	          ASR(IJM,IPOLYGEN),APR(IJM),BBP(IJM),XP(IJM)

	Dimension R1(IJM),V1(IJM),P1(IJM)
C-----   INITIAL VALUES
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J)
!$OMP DO
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      R1(I) = 0.0
CCC	         V1(I,K) = 0.0
	      P1(I) = 0.0
	   Endif
	Enddo
!$OMP END DO
!$OMP BARRIER
C-----   FORWARD
!$OMP DO
      Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      V1(I) = 0.0
	      Do J = 1, CELL_POLYGEN(I)
	         If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
                  V1(I) = V1(I) + ASL(I,J) * P1(CELL_SIDE(I,J,2))
	         Endif
	      Enddo
	      V1(I) = BBP(I) + V1(I)
C	      P1(I) = V1(I) / APL(I)
	   Endif
	Enddo
!$OMP END DO
!$OMP DO
      Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      P1(I) = V1(I) / APL(I)
	   Endif
      ENDDO
!$OMP END DO
!$OMP BARRIER
C-----   BACKWARD
C	R1(IJM,KBM) = P1(IJM,KBM) / P1(IJM,KBM)
!$OMP DO
	Do I = IJM, 1, -1
	   If(CCM(I) .EQ. 1.0) Then
	      V1(I) = 0.0
	      Do J = 1, CELL_POLYGEN(I)
	         If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
	            V1(I) = V1(I) + ASR(I,J) * R1(CELL_SIDE(I,J,2))
	         Endif
	      Enddo
	      
	      V1(I)= P1(I) + V1(I)
C	      R1(I) = V1(I) / APR(I)
	   Endif
	Enddo
!$OMP END DO
!$OMP DO
      Do I = IJM, 1, -1
	   If(CCM(I) .EQ. 1.0) Then
	      R1(I) = V1(I) / APR(I)
	   Endif
      ENDDO
!$OMP END DO 
!$OMP DO
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      XP(I) = R1(I)
	   Endif
	Enddo
!$OMP END DO
!$OMP END PARALLEL
	Return
	End