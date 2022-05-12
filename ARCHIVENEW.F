C############################################################################
C                                                                           #
C               SAVING THE COMPUTATIONAL RESULTS                            #
C                 VERSION 1.0 (10/08/2010)                                  #
C                 AUTHORIZED BY ZHANG JINGXIN                               #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
C############################################################################
	Subroutine ARCHIVE(DAYI,HOURSTAR)
      Include './Include/OCERM_INF'

	Dimension UINS(200,100),VINS(200,100),WINS(200,100),PINS(200,100),
     &            VORXINS(200,100), VORYINS(200,100), VORZINS(200,100),
     &            TXINS(200), ELFINS(200), VISINS(200,100) 	 

      Character*10 FN
      Integer NDEM

c-----   saving the instantenous results
c	write(fn,'(i6)') nstep
c	open(111,file='./results/'//fn//'.dat',status='unknown')
c	do k = 1, kbm
c	   write(111,'(2000f8.4)') (ur(i,k),vr(i,k),wr(i,k),i=15,17985,30)
c	enddo
c	close(111)
C-----   SAVING THE RESULTS AT FIXED OBSERVATIONAL POINTS   ----------------C
	  Do N = 1, EPTS
	     I = INXE(N)
c	     ESAVE(N) = ESAVE(N) + EL(I) * SKILLI 
	     esave(n) = el(i)
	     Do K = 1, KBM
	        USAVE(N,K) = USAVE(N,K) + UR(I,K) * SKILLI
	        VSAVE(N,K) = VSAVE(N,K) + VR(I,K) * SKILLI
	        WSAVE(N,K) = WSAVE(N,K) + WR(I,K) * SKILLI
	     Enddo
	  Enddo
c	  DO K = 1, KBM
c	  DO I = 1, 200
c	     II = 35 + 70 * (I - 1)
c	     UINS(I,K) = UINS(I,K) + UR(II,K) * SKILLI
c	     VINS(I,K) = VINS(I,K) + VR(II,K) * SKILLI
c	     WINS(I,K) = WINS(I,K) + WR(II,K) * SKILLI
c	     PINS(I,K) = PINS(I,K) - DC(II) * ZZ(K) * SKILLI
c	     VORXINS(I,K) = VORXINS(I,K) + VORX(II,K) * SKILLI
c	     VORYINS(I,K) = VORYINS(I,K) + VORY(II,K) * SKILLI
c	     VORZINS(I,K) = VORZINS(I,K) + VORZ(II,K) * SKILLI
c	     VISINS(I,K) = VISINS(I,K) + (KM(II,K) - UMOL) * SKILLI
c	  ENDDO
c	  ENDDO
c	  do i = 1, 200
c	     ii = 35 + 70 * (i - 1)
c		 TXINS(I) = TXINS(I) + TBX(II) * SKILLI
c	     ELFINS(I) = ELFINS(I) + EL(II) * SKILLI
c	  enddo
cc	OPEN(1002,FILE='2-CAL-STATIC.DAT',STATUS='UNKNOWN')
cc	OPEN(1004,FILE='4-CAL-STATIC.DAT',STATUS='UNKNOWN')
cc	OPEN(1005,FILE='5-CAL-STATIC.DAT',STATUS='UNKNOWN')
cc	OPEN(1006,FILE='6-CAL-STATIC.DAT',STATUS='UNKNOWN')
cc	OPEN(1012,FILE='12-CAL-STATIC.DAT',STATUS='UNKNOWN')
cc	OPEN(1018,FILE='18-CAL-STATIC.DAT',STATUS='UNKNOWN')

cc	  DO K = 1, KBM
cc	WRITE(1002,'(10F12.4)') UR(1560,K)/0.47,
cc     &	(DC(1560)*ZZ(K)+ELF(1560)+0.12)/0.12
C
cc	WRITE(1004,'(10F12.4)') UR(2840,K)/0.47,
cc     &	(DC(2840)*ZZ(K)+ELF(2840)+0.12)/0.12

cc	WRITE(1005,'(10F12.4)') UR(3480,K)/0.47,
cc     &	(DC(3480)*ZZ(K)+ELF(3480)+0.12)/0.12

cc	WRITE(1006,'(10F12.4)') UR(4120,K)/0.47,
cc     &	(DC(4120)*ZZ(K)+ELF(4120)+0.12)/0.12

cc	WRITE(1012,'(10F12.4)') UR(7960,K)/0.47,
cc     & 	(DC(7960)*ZZ(K)+ELF(7960)+0.12)/0.12

cc	WRITE(1018,'(10F12.4)') UR(11400,K)/0.47,
cc     &	(DC(11400)*ZZ(K)+ELF(11400)+0.12)/0.12

cc	  ENDDO
cc	stop




	  If(ISKILL .NE. 0 .AND. Mod(NSTEP,ISKILL) .EQ. 0 
     &	 .AND. NSTEP .GT. ISTART) Then
	     TMIDDLE = TIME - DTI * DAYI / SKILLI / 2.
C	Write(IUTSR,*) 'TITLE="Concentration of sediment" '
C	Write(IUTSR,*) 'VARIBLES= "X","Z","U","W","BED" '

C	  WRITE(IUTSR,'(A8,A10,A4,I4,A4,I4,A10)') 
C     $	    'ZONE T=',' "000" ',', I=',400,', K=',kbm,
C     $        ', F=POINT'
     
C	     DO K = 1, KBM
C	     DO I = 1, 400
C	        II = 20 + 40 * (I - 1)
C	        WRITE(IUTSR,'(20E12.4)')  
C     &			 CXY(II,1), DC(II) * ZZ(K)+ELF(II),UR(II,K), WR(II,K),  
C     &             ELF(II), PN(II,K), TBX(II), KM(II,K) - UMOL,
C     &             VORX(II,K), VORY(II,K), VORZ(II,K),VORTICITY(ii,k)
	        
C		    WRITE(IUTSR,'(20E15.6)')  
C     &			 UR(II,K), VR(II,K), WR(II,K),  
C     &             ELF(II), PN(II,K), TBX(II), KM(II,K) - UMOL,
C     &             VORX(II,K), VORY(II,K), VORZ(II,K),VORTICITY(ii,k)

c     &			 UINS(I,K), VINS(I,K), WINS(I,K),  
c     &             ELFINS(I), PINS(I,K), TXINS(I), VISINS(I,K),
c     &             VORXINS(I,K), VORYINS(I,K), VORZINS(I,K)
C	     ENDDO
C	     ENDDO
CC	stop
	     Write(IUTSR,5600) (TMIDDLE * 24. + HOURSTAR) * 3600.,
c     &	   (ESAVE(N), (USAVE(N,K), K = 1, KBM), WSAVE(N,K), 
C     &       N = 1, EPTS)
     &		  (ESAVE(N), N = 1, EPTS) 
	       
c     &          ,ur(i,5), ur(i,10), ur(i,15), ur(i,20)
C     &      (esave(n), ur(inxe(n),5), ur(inxe(n),10), ur(inxe(n),15),
C     &		   ur(inxe(n),20),n = 1, epts)
c	     DO K = 1, KBM
c	     DO I = 1, 200
c	        UINS(I,K) = 0.0
c	        VINS(I,K) = 0.0
c	        WINS(I,K) = 0.0
c	        PINS(I,K) = 0.0
c	        VORXINS(I,K) = 0.0
c	        VORYINS(I,K) = 0.0
c	        VORZINS(I,K) = 0.0
c	        VISINS(I,K) = 0.0
c	        TXINS(I) = 0.0
c	        ELFINS(I) = 0.0
c	     ENDDO
c	     ENDDO
		 
		 DO N = 1, EPTS
	        ESAVE(N) = 0.0
	        Do K = 1, KBM
	           USAVE(N,K) = 0.0
	           VSAVE(N,K) = 0.0
	           WSAVE(N,K) = 0.0
	        Enddo
	     Enddo
	  Endif
C-----   SAVING THE COMPUTIONAL RESULTS OF THE FLUID FIELD  ----------------C
	  If(Mod(NSTEP,ISAVE) .EQ. 0.0 .OR. NSTEP .EQ. ISTART) Then
    	    If(NSTEP .EQ. ISTART) Then
	      Write(FN,'(F10.3)')  TIME * 24. + HOURSTAR
	      If(KBM .GT. 1) Then
			 IJK = IJP * KB
	         IJCELL = IJM * KBM
	         WRITE(ITF,'(A8,A4,A10,A4,A4,I8,A4,I8,A10,A20)')
     &        'ZONE T=', '"',FN,'"',' I=',IJK,' K=',IJCELL,',F=FEPOINT',
     &        ', ET=BRICK'
	         Do K = 1, KB
	            Do I = 1, IJP
	               ZP = HP(I)+ELFV(I)
	               If(ZP .LE. 0.0) Then
				      ZP = 0.01
					  ZZZ = Z(K) * ZP + (ZP - HP(I))
				   Else
				      ZZZ = Z(K)*(HP(I)+ELFV(I))+ELFV(I)
				   Endif	       
	               Write(ITF,'(20E12.4)') 
     &	 	         PXY(I,1), PXY(I,2), ZZZ, ELFV(I), PNV(I,K),
     &               URV(I,K), VRV(I,K), WRV(I,K), VORTICITYV(I,K),
     &               VORXV(I,K), VORYV(I,K), VORZV(I,K),OMIGAV(I,K), 
     &               KMV(I,K),TKEV(I,K),TDISSV(I,K)
	            Enddo
	         Enddo
	         Do K = 1, KBM 
	            Do I = 1, IJM
	               Write(ITF,'(8I8)') 
     &                 IJP * (K-1) + CELL_INDEX(I,1), 
     &                 IJP * (K-1) + CELL_INDEX(I,2), 
     &                 IJP * (K-1) + CELL_INDEX(I,3), 
     &                 IJP * (K-1) + CELL_INDEX(I,4), 
     &                 IJP * (K) + CELL_INDEX(I,1), 
     &                 IJP * (K) + CELL_INDEX(I,2), 
     &                 IJP * (K) + CELL_INDEX(I,3), 
     &                 IJP * (K) + CELL_INDEX(I,4)
	            Enddo   
	         Enddo
	      Else
C-----      2D
	         WRITE(ITF,'(A8,A4,A10,A4,A4,I8,A4,I8,A10,A20)')
     &        'ZONE T=', '"',FN,'"',' I=',IJP,' K=',IJM,',F=FEPOINT',
     &        ', ET=QUADRILATERAL'
	         Do I = 1, IJP
	            Write(ITF,'(20E12.4)') 
     &	 	       PXY(I,1), PXY(I,2), ELFV(I), PNV(I,1),
     &               URV(I,1), VRV(I,1), KMV(I,1)
	         Enddo
	         Do I = 1, IJM
	            Write(ITF,'(8I8)') CELL_INDEX(I,1), CELL_INDEX(I,2), 
     &                               CELL_INDEX(I,3), CELL_INDEX(I,4) 
	         Enddo
	      Endif
	    Else
	      Write(FN,'(F10.4)')  TIME * 24. + HOURSTAR
	      If(KBM .GT. 1) Then
		    IJK = IJP * KB
	        IJCELL = IJM * KBM
	        WRITE(ITF,'(A8,A4,A10,A4,A4,I8,A4,I8,A10,A20,A20)')
     &        'ZONE T=', '"',FN,'"',' I=',IJK,' K=',IJCELL,',F=FEPOINT',
     &        ', ET=BRICK',', D=(FECONNECT)'
	        Do K = 1, KB
	           Do I = 1, IJP
		          ZP = HP(I)+ELFV(I)
                    If(ZP .LE. 0.0) Then
				     ZP = 0.01
			  	     ZZZ = Z(K) * ZP + (ZP - HP(I))
			  	  Else
				     ZZZ = Z(K)*(HP(I)+ELFV(I))+ELFV(I)
				  Endif	       
	              Write(ITF,'(20E12.4)') 
     &	 	         PXY(I,1), PXY(I,2), ZZZ, ELFV(I), PNV(I,K),
     &               URV(I,K), VRV(I,K), WRV(I,K), VORTICITYV(I,K),
     &               VORXV(I,K), VORYV(I,K), VORZV(I,K),OMIGAV(I,K), 
     &               KMV(I,K),TKEV(I,K),TDISSV(I,K)
	           Enddo
	        Enddo
	      Else
C-----      2D
	        WRITE(ITF,'(A8,A4,A10,A4,A4,I8,A4,I8,A10,A20,A20)')
     &        'ZONE T=', '"',FN,'"',' I=',IJP,' K=',IJM,',F=FEPOINT',
     &        ', ET=QUADRILATERAL',', D=(FECONNECT)'
	        Do I = 1, IJP
	           Write(ITF,'(20E12.4)') 
     &	 	         PXY(I,1), PXY(I,2), ELFV(I), PNV(I,1),
     &                 URV(I,1), VRV(I,1), KMV(I,1)
	        Enddo
	      Endif
	    Endif
	  Endif
	Return
 5600 Format (F15.5,1000E15.5)
	End