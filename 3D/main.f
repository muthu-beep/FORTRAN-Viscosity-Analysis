      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'NTMAX.h'
      INCLUDE 'NEIBMAX.h'
      INCLUDE 'NPCH.h'
	  
      PARAMETER(NVSMAX=NTMAX/2+2, NVTMAX=120000)
      PARAMETER (NEMAX=3*(NTMAX/2) )
      
      DIMENSION FT2(NVSMAX),IN0(6,NVSMAX),IN(NEIBMAX,NVSMAX),
     >NNEI(NVSMAX),PHI(NVTMAX),XNA(3,NVTMAX),X0(2,NPCH),
     >IV(3,NTMAX),V(3,NVTMAX),V1(3,NVTMAX),VP(3,NVTMAX),U(3,NVSMAX),
     >QW(3,NVTMAX),QWP(3,NVTMAX),QM(NVSMAX),VEL(3,NVSMAX),
     >VO(NVSMAX)
	 
      DIMENSION XV(3,NTMAX),YV(3,NTMAX),ZV(3,NTMAX),XP(3),P(3)
	 
      DIMENSION L1CH(NPCH), X0CH(2,NPCH), E1CH(2,NPCH),
     >VCH(2,NVTMAXCH), QCH(2,NVTMAXCH),
     >INBCH(NPCH),INFCH(NPCH),XNWCH(2,NPCH),DSWCH(NPCH)

       dimension im(NEMAX),JM(NEMAX),NT1IJ(NEMAX),
     >K1IJ(NEMAX),NT2IJ(NEMAX),K2IJ(NEMAX),NNTI(NVSMAX)
        LOGICAL RESTR
      
      CHARACTER*11 MNAME(0:999)
      CHARACTER*11 DRNAME(0:999)
    
      REAL *8 L1CH
      INTEGER RKS
      INTEGER ALPHAR, BETA 
	  
      PARAMETER (NP=1)
	  
      COMMON / BIVEL / X0CH,L1CH,E1CH,VCH,QCH,DSWCH,XNWCH,INBCH,
     >INFCH,NVTCH

      COMMON / XP / XP,RLENG,CURMAX,FS

      COMMON/QMINIP/QMINIP(NP) /RAT/RAT(NP) 

      CALL MNAMES(MNAME)
      CALL DRNAMES(DRNAME)

      PI=DACOS(-1.D0)	  
      CA=0.8D0
      RMU=1.D0
      A=3.D0

        h=0.1d0

        X0(1,1)=-4
        X0(2,1)=3

        X0(1,2)=-4
        X0(2,2)=-2
		
		X0(1,3)=0
        X0(2,3)=-2

        X0(1,4)=0
        X0(2,4)=-h
		
        x0(1,5)=h+h*COS(5.d0/6.d0*PI)
        x0(2,5)=-h+h*SIN(5.d0/6.d0*PI)

        x0(1,6)=h+h*COS(2.d0/3.d0*PI)
        x0(2,6)=-h+h*SIN(2.d0/3.d0*PI)

        x0(1,7)=h
        x0(2,7)=0

        X0(1,8)=2-h
        X0(2,8)=0		
		
        X0(1,9)=(2-h)+h*COS(PI/3.d0)
        X0(2,9)=(-h)+h*SIN(PI/3.d0)

        x0(1,10)=(2-h)+h*COS(PI/6.d0)
        x0(2,10)=(-h)+h*SIN(PI/6.d0)

        x0(1,11)=2
        x0(2,11)=-h

        x0(1,12)=2
        x0(2,12)=-2

        x0(1,13)=6
        x0(2,13)=-2

        x0(1,14)=6
        x0(2,14)=3
		
		x0(1,15)=2
		x0(2,15)=3
		
		x0(1,16)=2
		x0(2,16)=1+h
		
		x0(1,17)=(2-h)+h*COS(11.d0*PI/6.d0)
		x0(2,17)=(1+h)+h*SIN(11.d0*PI/6.d0)
		
		x0(1,18)=(2-h)+h*COS(5.d0*PI/3.d0)
		x0(2,18)=(1+h)+h*SIN(5.d0*PI/3.d0)
		
		x0(1,19)=2-h
		x0(2,19)=1

       x0(1,20)=h
       x0(2,20)=1

       x0(1,21)=h+h*COS(4.D0*PI/3.D0)
       x0(2,21)=(1+h)+h*SIN(4.D0*PI/3.D0)

       x0(1,22)=h+h*COS(7.d0*PI/6.d0)
       x0(2,22)=(1+h)+h*SIN(7.d0*PI/6.d0)

       x0(1,23)=0
       x0(2,23)=1+h

       x0(1,24)=0
       x0(2,24)=3

        OPEN (13, FILE='2Dchan.dat')
        do ip=1,npch
        write (13,*) x0(1,ip),x0(2,ip)
        enddo
        write (13,*) x0(1,1),x0(2,1)
        CLOSE (13)

        OPEN (15, FILE='3Dchan.dat')
        WRITE (15,*)'VARIABLES = "X","Y","Z"'
        WRITE (15,322) 2*NPCH, NPCH
322     FORMAT ('ZONE N=', I5,'  E=',I4,' F=FEPOINT, ET=QUADRILATERAL')

        do ip=1,npch
        ZERO= 0.D0
        ONE=1.D0
        write (15,*) x0(1,ip),x0(2,ip),ZERO
        write (15,*) x0(1,ip),x0(2,ip),ONE
        enddo

        DO I=1,2*NPCH,2
        IF (I.EQ.2*NPCH-1)THEN
        n2=2
        n1=1
        WRITE (15,*) I, I+1, 2, 1
        ELSE
        WRITE (15,*) I, I+1,I+3, I+2
        END IF
        ENDDO

        CLOSE(15)

	  
      OPEN (10,file ='SOL2D.dat')
      READ (10,*) X0CH,L1CH,E1CH,VCH,QCH,DSWCH,XNWCH,INBCH,INFCH,NVTCH
      CLOSE(10)
	  
	  FSTEP=1.D0
c!       
          FSTEP=0.3D0
	  ALPHA=0.D0
      TOL=1.D-5
      C=0.18D0/4.D0
      COEFF1=1.D0
      COEFF2=0.1D0
      GAMMA=0.25D0

C------------------------------------------------------------------------
C      FOR RECRESTR:
        QLIM=0.09
        ALPHAR=50
        BETA=100
        COEFF= (0.144D0)**BETA
        CMESH=0.01
        ITRESTR=2000
C---------------------------------------------------------------------------------------------
      NST= 20

      OPEN(1,FILE='restart.dat', form='UNFORMATTED')
      READ(1) T,DTL,DTP,NTS,IV,NNEI,IN,NVT,NVTP,RAD,V,VP,XNA,QW,QWP,
     >VEL, IM,JM,K1IJ,K2IJ,NT1IJ,NT2IJ

      CLOSE(1)
      NVS=NTS/2+2
      NES= 3*(NTS/2)
    
      OPEN(26, FILE='TEMP1.DAT')
 
      DO 1000 KST=1,NST

      IF (MOD(KST,5).eq.1) then
      
       CALL RECRESTR(.TRUE.,1,QLIM,ALPHAR,beta,GAMMA,coeff,CMESH,
     >ITRESTR,V,XNA,QW,VEL,1,NVS,1,NTS,1,NES,IM,JM,
     >NT1IJ,K1IJ,NT2IJ,K2IJ,IN,nnei,IV,
     >QMINIP,DXMINIP,DXMAXIP,ITRESTIP, RESTR)
      print *, ' ITRESTIP=', ITRESTIP

       IF (RESTR) THEN
          WRITE (4,*) ' RECRESTR='
          DO I= 1,NVS
          DO K=1,3
          QWP(K,I)=QW(K,I)
          ENDDO
          ENDDO
          print *, ' QMINIP=', QMINIP
c          read (*, *) icont
       END IF
       END IF
 
      RKS=1
	  
      CALL VEL3D (RAD,A,ALPHA,CA,RMU,NVCONT,TOL,DTL,DTP,
     >FT2,IV,NNEI,IN,NVT,NVTP,NTS,PHI,C,V,VP,XNA,QW,QWP,U,QM)
c       print *, 'start VOLUMES'
c      CALL VOLUMES (NTS,IV, V, XNA, VOL1, VOL2, VREL, VOL)
c      stop
	 
      CALL QMINSIAD(COEFF1,COEFF2,GAMMA,IN,nnei,IV,V,
     >1,NVS,1,NTS,XNA,PHI,QM,FT2,VEL,RATMAX)
	 
      DO I=1,NVS
      VO(I)=0.D0
      ENDDO

c      VOL=0.D0

c      DO I=1,NVS
c      DO K=1,3
c      VO(I)=(V(K,I)-XP(K))*XNA(K,I)+VO(I)
c      ENDDO
c      VOL=VO(I)*PHI(I)+VOL
c      ENDDO

c      PRINT *, 'VOL=', VOL
 
c      IF (KST.EQ.1) THEN
c      VOL0=VOL
      
c      ELSE

c      FACTOR=(VOL0/VOL)**(1.D0/3.D0)
c      PRINT *, 'FACTOR=', FACTOR
c      DO I=1,NVS
c      DO K=1,3
c      V(K,I)=XP(K)+FACTOR*(V(K,I)-XP(K))
c      ENDDO
c      ENDDO

c      END IF

      WRITE(26,*) 'T=',T, ' XP(1)=',XP(1), ' XP(2)=',XP(2), ' RMS=',RMS
      WRITE(26,*) 'T=',T,' LENG=',LENG,' CURMAX=',CURMAX,' RATMAX=', RAT
c      CLOSE(26)

      USM=0.D0
      DO I=1,NVS
      USM=USM+QM(I)**2*PHI(I)
      ENDDO
      RMS=(USM/FS)**(1.d0/2.d0)

       print *, ' KST=', KST,' T=', T, ' DT=', DT
      PRINT *, 'XP=', XP(1),XP(2), 'RLENG=', RLENG
      PRINT *, 'T=', T, 'RMS=', RMS
      PRINT *, 'C/R=', CURMAX


      DX=1.D20

      DO I=1,NVS
      DO J=1,NNEI(I)

      L=IN(J,I)
      R=0.D0

      DO K=1,3
      P(K)=V(K,I)-V(K,L)  
      R=R+P(K)**2  
      ENDDO

      R=DSQRT(R) 

      IF (R.LT.DX) DX=R
 
      ENDDO
      ENDDO

      DT=FSTEP*CA*DX 
	   
      IF (RKS.EQ.1) THEN
         DO I=1,NVS
         DO K=1,3
         V(K,I)=V(K,I)+VEL(K,I)*DT
         ENDDO
         ENDDO
	  
	  DTL=DT
       
      ELSE
      
         DO I=1,NVS
         DO K=1,3
         V1(K,I)=V(K,I)
         V(K,I)=V(K,I)+VEL(K,I)*DT/2
         ENDDO
         ENDDO
	  
      CALL VEL3D (RAD,A,ALPHA,CA,RMU,NVCONT,TOL,DT/2.D0,DTP,
     >FT2,IV,NNEI,IN,NVT,NVTP,NTS,PHI,C,V,VP,XNA,QW,QWP,U,QM)
      CALL QMINSIAD(COEFF1,COEFF2,GAMMA,IN,nnei,IV,V1,
     >1,NVS,1,NTS,XNA,PHI,QM,FT2,VEL,RATMAX)

	  DO I=1,NVS
      DO K=1,3
      V(K,I)=V1(K,I)+VEL(K,I)*DT
      ENDDO
      ENDDO
	  
        DTL=DT/2.D0
      END IF

      T=T+DT
	  
      pict= 0.1d0
321   FORMAT ('ZONE N=', I5, '  E=', I5, '  F=FEPOINT, ET=TRIANGLE')
86    FORMAT (1X, 3F13.7)

      OPEN (7, FILE='MESHCURRENT.DAT')
      OPEN (8, FILE='DR2DCURRENT.DAT')
      WRITE (7,*) 'VARIABLES = "X","Y","Z"'
C      WRITE (8,*) 'VARIABLES = "X","Y"'
      WRITE (7,321)NVS, NTS
c      WRITE (8,321)NVS, NTS
      DO I=1,NVS
      WRITE (7,86) V(1,I),V(2,I),V(3,I)
      WRITE (8,86) V(1,I),V(2,I)
       ENDDO
         DO NT=1,NTS
      WRITE (7,*) IV(1,NT), IV(2,NT), IV(3,NT)
c      WRITE (8,*) IV(1,NT), IV(2,NT), IV(3,NT)
        ENDDO
      close (7)
      CLOSE (8)

      DO K=0,999
      IF(T.GE.K*PICT.AND.(T-DT).LT.(K+1.D-7)*PICT)
     >THEN
      OPEN (7, FILE=MNAME(K))
      OPEN (8, FILE=DRNAME(K))
      WRITE (7,*) 'VARIABLES = "X","Y","Z"'
C      WRITE (8,*) 'VARIABLES = "X","Y"'
      WRITE (7,321)NVS, NTS
c      WRITE (8,321)NVS, NTS
      DO I=1,NVS
      WRITE (7,86) V(1,I),V(2,I),V(3,I)
      WRITE (8,86) V(1,I),V(2,I)
       ENDDO
         DO NT=1,NTS
      WRITE (7,*) IV(1,NT), IV(2,NT), IV(3,NT)
c      WRITE (8,*) IV(1,NT), IV(2,NT), IV(3,NT)
        ENDDO
      close (7)
      CLOSE (8)
      END IF
      ENDDO

      IF (XP(1).GT.3.D0) THEN
      OPEN (17, FILE='FEND')
      WRITE(17,*)'DROP IN STRAIGHT CHANNEL'
      STOP
      END IF

      IF (XP(2).LT.-3.D0) THEN
      OPEN (17, FILE='FEND')
      WRITE(17,*)'DROP IN SIDE CHANNEL'
      STOP
      END IF 


1000  CONTINUE

      OPEN(1,FILE='restart.dat', form='UNFORMATTED')
      WRITE(1) T,DTL,DTP,NTS,IV,NNEI,IN,NVT,NVTP,RAD,V,VP,XNA,QW,QWP,
     >VEL,  IM,JM,K1IJ,K2IJ,NT1IJ,NT2IJ
      CLOSE(1)
      stop
      end

      SUBROUTINE VOLUMES (NTS,IV, V, XNA, VOL1, VOL2, VREL, VOL)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'NTMAX.h'
      PARAMETER (NVSMAX= NTMAX/2+2, NVTMAX= 120000)
      DIMENSION IV(3,*), V(3,*), XNA(3,*)   
      DIMENSION VL(3), XCOR(3), R(3,3), CC(2,3) ,P(3)

      PI=DACOS(-1.D0)  

       VL(1)=COS(PI/4.D0)
       VL(2)=SIN(PI/4.D0)
       VL(3)=0.D0

       XCOR(1)=10317.D0/10000.D0
       XCOR(2)=-317.D0/10000.D0
       XCOR(3)=0.D0

      VOL1=0.D0
      VOL2=0.D0

       DO 555 JT=1,NTS
        
       DO M=1,3
       I=IV(M,JT)

      VO=0.D0 
             
       DO K=1,3
       R(K,M)= (V(K,I)- XCOR(K))
       VO= R(K,M)*VL(K)+VO
       ENDDO

       P(M)=VO

       ENDDO

      IF (MIN(P(1),P(2),P(3)).LT.0.D0.AND.MAX(P(1),P(2),P(3)).GT.0.D0)
     >GOTO 555

       DO K=1,2
       DO J=1,3
       JJ1=IV(K+1,JT)
       JJ2=IV(1,JT)
       CC(K,J)=V(J,JJ1)-V(J,JJ2)
       ENDDO
       ENDDO

       CX1=CC(1,2)*CC(2,3)-CC(2,2)*CC(1,3)
       CY1=-CC(1,1)*CC(2,3)+CC(2,1)*CC(1,3)
       CZ1=CC(1,1)*CC(2,2)-CC(2,1)*CC(1,2)
       DS= .5D0*DSQRT(CX1**2+CY1**2+CZ1**2)


       DO M=1,3
       I=IV(M,JT)

       IF (max(P(1),P(2),P(3)).LT.0.d0 ) THEN
       P1=0.d0
       DO K=1,3
       P1=P1+(R(K,M)*XNA(K,I))
       ENDDO
       VOL1=VOL1+P1*DS/9.D0
       ELSE

CCCC       IF (min(P(1),P(2),P(3)).GT.0.d0 ) THEN
       P2=0.d0
       DO K=1,3
       P2=P2+(R(K,M)*XNA(K,I))
       ENDDO
       VOL2=VOL2+P2*DS/9.D0
       END IF

       ENDDO
       
555   CONTINUE

       VOL=VOL1+VOL2
       VREL=VOL1/VOL2
       VREL1= VOL1/VOL
       VREL2=VOL2/VOL

       PRINT *, 'VREL= ', VREL
       PRINT *, 'VREL1= ', VREL1
       PRINT *, ' VREL2= ', VREL2
       PRINT *, 'VOL= ', VOL

      RETURN
      END
