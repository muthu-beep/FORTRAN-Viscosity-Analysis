      SUBROUTINE VEL3D (RAD,A,ALPHA,CA,RMU,NVCONT,TOL,DT,DTP,
     >FT2,IV,NNEI,IN,NVT,NVTP,NTS,PHI,C,V,VP,XNA,QW,QWP,U,QM)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'NTMAX.h'
      INCLUDE 'NEIBMAX.h'
      INCLUDE 'NPCH.h'
      INCLUDE 'NRESMAX.h'
      PARAMETER (NVSMAX= NTMAX/2+2,
     >NVTMAX=120000, NPCMAX=25,NVCONTMAX=3000,NP0=NPCMAX-1)
!     CHANGE REAL TO MATCH ORIGINAL
      REAL*8 LX,LY,LZ,L1,L2,L1CH
      DIMENSION L1CH(NPCH), X0CH(2,NPCH), E1CH(2,NPCH)
      DIMENSION PHI(*),FT2(*),IN(NEIBMAX,*),FIT(3,NVSMAX),
     >CUR(NVSMAX),F(3,NVTMAX),RDX(3),RDY(3),RX1(3),RY1(3),
     >NNEI(*),PX1(3),PY1(3)
      DIMENSION V(3,*),IV(3,*),
     >XNA(3,*),VOLD(3,NVTMAX),
     >QW(3,*),QWIMP(3,NVTMAX),VINF(3,NVTMAX),IM(NVTMAX,0:NP0),
     >RMIN(NVTMAX),
     >DL(3,NVTMAX),SX(3,NVTMAX),SY(3,NVTMAX),QWOLD(3,NVTMAX),
     >XP(3),XCW(3),INB(0:NP0),INF(0:NP0),INBH(0:NP0),INFH(0:NP0),
     >INV(NVTMAX),XNW(3,NP0),CC(2,3),UINF(3,NVTMAX),DM(3,3,2)
      DIMENSION L1(NP0), L2(NP0), X0(3,NP0), N1(NP0), N2H(NP0),
     >E1(3,NP0), E2(3,NP0), RRMINI(NVSMAX)
      DIMENSION DX(3),U(3,NVSMAX),OM(3),SM(3),STR(3,3)
      DIMENSION DQW(3,NVTMAX),EA(3,NVTMAX,NRESMAX+1),
     >TM(NRESMAX,NRESMAX+1),RLAM(NRESMAX),
     >ACOEFF(NRESMAX,NRESMAX+1),RCOEFF(NRESMAX+1)
      DIMENSION X0C(2,NPCMAX),INB1(NP0),INF1(NP0),V1(2,NVCONTMAX),
     >DS1(NVCONTMAX), X(2)
      DIMENSION VP(3,*),QWP(3,*),QM(*)
      DIMENSION IFRONT(NVTMAX),IBACK(NVTMAX),IUP(NVTMAX),IDOWN(NVTMAX)
      LOGICAL SOL(NPCMAX)
      CHARACTER *7 PNAME(NP0)

      DIMENSION
     >VCH(2,NVTMAXCH), QCH(2,NVTMAXCH),
     >INBCH(NPCH),INFCH(NPCH),XNWCH(2,NPCH),DSWCH(NPCH)

      COMMON / BIVEL / X0CH,L1CH,E1CH,VCH,QCH,DSWCH,XNWCH,INBCH,
     >INFCH,NVTCH

      COMMON / XP / XP,RLENG,CURMAX,FS

      NVS=NTS/2+2
      NVSP=NVS

      LZ=2*A
      NGMRES=NRESMAX
      A0=A

      DO I=1,NVS
      PHI(I)=0.D0
      ENDDO

      FS=0.D0
      DO JT=1,NTS
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
      FS=FS+DS
      DO KJ=1,3
      I= IV(KJ,JT)
      PHI(I)=PHI(I)+DS/3.D0
      ENDDO
      ENDDO

      DO K=1,3
      XP(K)=0
      ENDDO

      DO I=1,NVS
      DO K=1,3
      XP(K)=XP(K)+V(K,I)*PHI(I)
      ENDDO
      ENDDO

      SUR=0.D0

      DO I=1,NVS
      SUR=SUR+PHI(I)
      ENDDO

      DO K=1,3
      XP(K)=XP(K)/SUR
      ENDDO

      XP(3)=0.D0

      RD=0.D0

      DO I=1,NVS
      RD1=((V(1,I)-XP(1))**2+(V(2,I)-XP(2))**2+
     >(V(3,I)-XP(3))**2)**(1.D0/2.D0)
      RD=MAX(RD1,RD)
      ENDDO
      print *, ' RD=', RD

      A1=MAX(A*RAD,1.5*RD)
      A0=A1
      nfail=0

      ANGLE=30.D0

345   CALL CONTGEN (XP,A0,ANGLE,X0CH,E1CH,L1CH,NPC,X0C,SOL)

      OPEN (13, FILE='CELL.DAT')
      DO IP=1,NPC
      WRITE (13,*) X0C(1,IP), X0C(2,IP)
      ENDDO
      CLOSE(13)

      NP=NPC-1
      WRITE (4,*) ' NP=', NP
      PRINT *, ' NP=', NP, 'NPC= ', NPC
      TOTLENG=0.D0
      DO IP=1,NP
      JP=IP+1
      IF (JP.GT.NP) JP=JP-NP
      TMP=0.D0
      DO K=1,2
      TMP=TMP+ (X0C(K,JP)-X0C(K,IP))**2
      ENDDO
      L1(IP)=DSQRT(TMP)
      DO K=1,2
      E1(K,IP)= (X0C(K,JP)-X0C(K,IP))/L1(IP)
      X0(K,IP)=X0C(K,IP)
      E2(K,IP)=0.D0
      ENDDO
      X0(3,IP)= -LZ/2.D0
      E1(3,IP)=0.D0
      E2(3,IP)=1.D0
      L2(IP)=LZ
      TOTLENG=TOTLENG+L1(IP)
      ENDDO


C      DO IP=1,NP
C      IF (L1(IP).LT. 1.001*C) THEN
C        nfail=nfail+1
C        if (nfail.eq.1) then
C        A0= 0.98*A0
C        else
C        A0= 1.025*A0
C        end if

C       PRINT *, ' RECONTOR A0=',A0, ' nfail=', nfail
C       GOTO 345
C      END IF
C      ENDDO



      WRITE (4,*) 'TOTLENG=', TOTLENG
      WRITE (4,*) ' NVCONT EST=', NVCONT
      PRINT *, 'TOTLENG=', TOTLENG

C     ------------------------------
      I=0
      DO 44 IP=1,NP
      S=0.D0
      INEW=I
      INB1(IP)=I+1
      WRITE (4,*) ' IP=', IP, ' INB1=', INB1(IP)
      PRINT *, ' IP=', IP, ' INB1=', INB1(IP)
      DO K=1,2
      X(K)=0.D0
      ENDDO

45    RR=0.D0
      DO K=1,2
      RR=RR+(X0(K,IP)+X(K)-XP(K))**2
      ENDDO
      GAP=1.D0
      H=C*GAP**ALPHA
      S=S+H
      IF (S.GT.L1(IP)) THEN
	S=S-H
	DO J=I+1,INEW
	RAT=L1(IP)/S
	DS1(J)=DS1(J)*RAT
	DO K=1,2
	V1(K,J)=X0(K,IP)+RAT*V1(K,J)+ 0.5D0*DS1(J)*E1(K,IP)
	ENDDO
	ENDDO
	GOTO 46
      ELSE
	INEW=INEW+1
	IF (INEW.GT.NVCONTMAX) THEN
	  WRITE (4,*) ' NVCONTMAX TOO SMALL'
	  PRINT *, ' NVCONTMAX TOO SMALL, INCREASE C'
	  STOP
	END IF
        DS1(INEW)=H
        DO K=1,2
        V1(K,INEW)=X(K)
        X(K)=X(K)+H*E1(K,IP)
        ENDDO
        GOTO 45
      END IF

46    N1(IP)=INEW-I
      I=INEW
      INF1(IP)=INEW
44    CONTINUE
      NVCONT=I
      WRITE (4,*) '  NVCONT=', NVCONT
      PRINT *, '  NVCONT=', NVCONT

      OPEN (13, FILE='VCELL.DAT')
      DO IP=1,NP
      PRINT *, ' IP=', IP, ' INB2=', INB1(IP), ' INF1=', INF1(IP)
      PRINT *, ' N1=', N1(IP)
      DO I=INB1(IP), INF1(IP)
      WRITE (13,*) V1(1,I),V1(2,I)
      SUM=SUM+DS1(I)
      ENDDO
      PRINT *, ' IP=', IP, ' SUM=', SUM, ' L=', L1(IP)
      ENDDO
      CLOSE(13)

      NVTOLD=NVT
      NVSOLD=NVS
      DO J=1,NVT
      DO K=1,3
      VOLD(K,J)=V(K,J)
      QWOLD(K,J)=QW(K,J)
      ENDDO
      ENDDO

321   FORMAT ('ZONE N=', I5, '  E=', I5, '  F=FEPOINT, ET=TRIANGLE')
86    FORMAT (1X, 3F13.7)
      OPEN (7, FILE='DROPMESH.DAT')
c      WRITE (7,*) 'VARIABLES = "X","Y","Z"'
      WRITE (7,*) 'VARIABLES = "X","Y"'
      WRITE (7,321)NVS, NTS
         DO I=1, NVS
c         WRITE (7,86) V(1,I),V(2,I),V(3,I)
         WRITE (7,86) V(1,I),V(2,I)
         ENDDO
         DO NT=1,NTS
         WRITE (7,*) IV(1,NT), IV(2,NT), IV(3,NT)
        ENDDO
      close (7)
      INB(0)=1
      INF(0)=NVS
      INBH(0)=1
      INFH(0)=NVS

      CALL GIM (INB,INF,V,IN,NNEI,FIT,XNA)
      CURMAX=0.d0
      DO I=1,NVS
      CUR(I)=-FIT(1,I)-FIT(3,I)
      FT2(I)=FIT(1,I)**2+FIT(3,I)**2+FIT(2,I)**2/2.D0
      curmax= max(curmax,CUR(I))
      ENDDO
      CURMAX=CURMAX*RAD

       FS=SUR
       PRINT *, ' FS=', FS

	I=NVS
	DO IP=1,NP
	  WRITE (4,*) 'IP=', IP
      PRINT *, 'IP=', IP
      CALL PNLSYMAA(I,NVTMAX,INB1(IP),INF1(IP),V1,DS1,XP,C,RAD,ALPHA,
     >L2(IP),X0(1,IP),E1(1,IP),INB(IP),INF(IP),INBH(IP),INFH(IP),
     >INV,v,XNW(1,IP),PHI,IFRONT,IBACK,IUP,IDOWN)

      WRITE (4,*) 'IP=', IP, ' INB=', INB(IP), ' INF=', INF(IP)
      PRINT *, 'IP=', IP, ' INB=', INB(IP), ' INF=', INF(IP)
      WRITE (4,*) ' '
      PRINT *, ' '
	ENDDO

	NVT=I
       WRITE (4,*)' NVT=', NVT
      PRINT *, ' NVT=', NVT

      IF (NVTOLD.EQ.0) THEN
	DTP=1.D20
        NVTP=NVT
        NVSP=NVS
	DO I=1,NVT
	DO K=1,3
	QW(K,I)=0.D0
	QWP(K,I)=0.D0
	VP(K,I)=V(K,I)
	ENDDO
	ENDDO

      ELSE
      DO I=1,NVS
      DO K=1,3
      QW(K,I)=QW(K,I)+ (QW(K,I)-QWP(K,I))*DT/DTP
      ENDDO
      ENDDO

C       -----------------------------------
	PMIN=1.D20
        DO IP=1,NP
        DO I=INBH(IP),INFH(IP)
	RRMIN=1.D20
	DO 66 J=NVSOLD+1,NVTOLD
	IF (VOLD(3,J).GT.0.D0) GOTO 66
	P=0.D0
	DO K=1,3
	P=P+(VOLD(K,J)-V(K,I))**2
	ENDDO
	IF (P.LT.RRMIN) THEN
	  RRMIN=P
	  JMIN=J
	eND IF
66      CONTINUE

        RRNEXT=1.D20
        DO 67 J=NVSOLD+1,NVTOLD
        IF (VOLD(3,J).GT.0.D0. OR. J.EQ.JMIN) GOTO 67
        P=0.D0
        DO K=1,3
        P=P+(VOLD(K,J)-V(K,I))**2
        ENDDO
        IF (P.LT.RRNEXT) THEN
          RRNEXT=P
          JNEXT=J
        END IF
67      CONTINUE
        R2=0.D0
        P=0.D0
        DO K=1,3
        R2=R2+(VOLD(K,JMIN)-VOLD(K,JNEXT))**2
        P=P+ (V(K,I)-VOLD(K,JNEXT))*(VOLD(K,JMIN)-VOLD(K,JNEXT))
        ENDDO
        P=P/R2
        PMIN=MIN(PMIN,P)

	DO K=1,3
        QW(K,I)=P*QWOLD(K,JMIN)+(1-P)*QWOLD(K,JNEXT)
        ENDDO

	RRMIN=1.D20
	DO 660 J=NVSP+1,NVTP
	IF (VP(3,J).GT.0.D0) GOTO 660
	P=0.D0
	DO K=1,3
	P=P+(VP(K,J)-V(K,I))**2
	ENDDO
	IF (P.LT.RRMIN) THEN
	  RRMIN=P
	  JMIN=J
	END IF
660     CONTINUE

        RRNEXT=1.D20

        DO 670 J=NVSP+1,NVTP
        IF (VP(3,J).GT.0.D0. OR. J.EQ.JMIN) GOTO 670
        P=0.D0
        DO K=1,3
        P=P+(VP(K,J)-V(K,I))**2
        ENDDO
        IF (P.LT.RRNEXT) THEN
          RRNEXT=P
          JNEXT=J
        END IF
670     CONTINUE
        R2=0.D0
        P=0.D0
        DO K=1,3
        R2=R2+(VP(K,JMIN)-VP(K,JNEXT))**2
        P=P+ (V(K,I)-VP(K,JNEXT))*(VP(K,JMIN)-VP(K,JNEXT))
        ENDDO
        P=P/R2
        PMIN=MIN(PMIN,P)

	DO K=1,3
        TEMP=P*QWP(K,JMIN)+(1-P)*QWP(K,JNEXT)
	QW(K,I)=QW(K,I)+ (QW(K,I)-TEMP)*DT/DTP
        ENDDO
        I1=INV(I)
        QW(1,I1)=  QW(1,I)
        QW(2,I1)=  QW(2,I)
        QW(3,I1)= -QW(3,I)
        ENDDO
        ENDDO

        IF (DABS(DT).GT. 1.D-14) THEN
	  NVTP=NVTOLD
	  NVSP=NVSOLD
	  DO I=1, NVTP
	  DO K=1,3
	  VP(K,I)=VOLD(K,I)
	  QWP(K,I)=QWOLD(K,I)
	  ENDDO
	  ENDDO
	  DTP=DT
        END IF
      END IF
      PRINT  *, ' NVT=', NVT, ' PMIN=', PMIN

        PNAME(01)='P01.DAT'
        PNAME(02)='P02.DAT'
        PNAME(03)='P03.DAT'
        PNAME(04)='P04.DAT'
        PNAME(05)='P05.DAT'
        PNAME(06)='P06.DAT'
        PNAME(07)='P07.DAT'
        PNAME(08)='P08.DAT'
        PNAME(09)='P09.DAT'
        PNAME(10)='P10.DAT'
        PNAME(11)='P11.DAT'
        PNAME(12)='P12.DAT'
        PNAME(13)='P13.DAT'
        PNAME(14)='P14.DAT'
        PNAME(15)='P15.DAT'
        PNAME(16)='P16.DAT'
        PNAME(17)='P17.DAT'
        PNAME(18)='P18.DAT'
        PNAME(19)='P19.DAT'
        PNAME(20)='P20.DAT'
        PNAME(21)='P21.DAT'
        PNAME(22)='P22.DAT'
        PNAME(23)='P23.DAT'
        PNAME(24)='P24.DAT'
        PNAME(25)='P25.DAT'


        DO IP=1,NP
        OPEN (11,FILE=PNAME(IP))
        DO I=INB(IP),INF(IP)
        X1=0.D0
        DO K=1,3
        X1=X1+(V(K,I)-X0(K,IP))*E1(K,IP)
        ENDDO
        WRITE (11,*) X1,V(3,I)-X0(3,IP)
        ENDDO
        CLOSE (11)
        ENDDO

      DO I=1,NVS
      CALL VELOCITY(V(1,I),UINF(1,I))
C      OPEN (11,FILE = 'VEL2D.DAT')
C      WRITE (11,*)' I=', I,  UINF(1,I),UINF(2,I)
C      WRITE (4,*) ' I=', I,  UINF(1,I),UINF(2,I)
!      PRINT *, ' I=', I,  UINF(1,I),UINF(2,I)
!          OPEN (11,FILE = 'VEL2D.DAT')
!          WRITE(11,126) I,T,XP(1),XP(2),XCW(1),XCW(2),XCW(3),UINF(1,I),UINF(2,I)
!126       FORMAT(I8,8F12.8)


      UINF(3,I)=0.D0
      ENDDO

      DO K=1,3
      XCW(K)=0.D0
      ENDDO
      SW=0.D0
      DO IP=1,NP
      DO I=INB(IP),INF(IP)
      DO K=1,3
      XNA(K,I)= XNW(K,IP)
      XCW(K)=XCW(K)+V(K,I)*PHI(I)
      ENDDO
      SW=SW+ PHI(I)
      ENDDO
      ENDDO

      DO K=1,3
      XCW(K)=XCW(K)/SW
      ENDDO
      WRITE (4,*) ' XCW=', XCW
      PRINT *, ' XCW=', XCW

      SBOTTOM=0.D0
      DO IP=1,NP
      TEMP=0.D0
      DO K=1,2
      TEMP=TEMP+X0C(K,IP)*XNW(K,IP)
      ENDDO
      SBOTTOM=SBOTTOM+TEMP*L1(IP)/2.D0
      ENDDO
      PRINT *, ' SBOTTOM=', SBOTTOM, ' SW=', SW
      IF (SBOTTOM.LT.0.D0) STOP


C      CALL DEFLMAT (1,NVS,XP,V,PHI,DM(1,1,1))
      CALL DEFLMAT (NVS+1,NVT,XCW,V,PHI,DM(1,1,2))

!     DYNAMICS
      PI=DACOS(-1.D0)
      PRINT *, ' START IM'
      DO IP=1,NP
      DO JP=1,NP
      DO I=INBH(IP),INFH(IP)
      RRMIN=1.D20
      DO J=INBH(JP),INFH(JP)
      P=0.D0
      DO K=1,3
      P=P+ (V(K,J)-V(K,I))**2
      ENDDO
      IF (P.LT.RRMIN) THEN
	RRMIN=P
	IM(I,JP)=J
      END IF
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      DO JP=1,NP
      DO I=1,NVS
      RRMINI(I)=1.D20
      ENDDO
      DO J=INB(JP),INF(JP)
      RRMIN=1.D20
      DO I=1,NVS
      P=0.D0
      DO K=1,3
      P=P+(V(K,J)-V(K,I))**2
      ENDDO
      IF (P.LT.RRMIN) THEN
        RRMIN=P
        IM(J,0)=I
      END IF
      IF (P.LT.RRMINI(I)) THEN
        RRMINI(I)=P
        IM(I,JP)=J
      END IF
      ENDDO
      ENDDO
      ENDDO

       FFACT=-1.D0/(4.D0*PI*CA)

!     F-CALCULATION: SELF-INTERACTION

       DO I=1,NVS
       DO K=1,3
       F(K,I)=0.D0
       ENDDO
       ENDDO

       RLENG=0.d0
       DO I=1,NVS
       DO J=I+1,NVS
       P=0.D0
       RNX=0.D0
       RNY=0.D0

       DO K=1,3
       DX(K)=V(K,J)-V(K,I)
       P=P+DX(K)**2
       RNX=RNX+DX(K)*XNA(K,J)
       RNY=RNY+DX(K)*XNA(K,I)
       ENDDO
       RLENG= max (RLENG,p)
       DO K=1,3
       RDX(K)=XNA(K,J)/DSQRT(P)
       RDY(K)=XNA(K,I)/DSQRT(P)
       ENDDO

       DO K=1,3
       RX1(K)=RNX*DX(K)
       RY1(K)=RNY*DX(K)
       ENDDO

       R= P*DSQRT(P)

       DO K=1,3
       RX1(K)=RX1(K)/R
       RY1(K)=RY1(K)/R
       ENDDO

       DO K=1,3
       PX1(K)=(RX1(K)+RDX(K))*PHI(J)
       PY1(K)=(RY1(K)+RDY(K))*PHI(I)
       ENDDO

       DO K=1,3
       F(K,I)= F(K,I)+PX1(K)*(CUR(J)-CUR(I))*FFACT
       F(K,J)= F(K,J)+PY1(K)*(CUR(I)-CUR(J))*FFACT
       ENDDO
       ENDDO
       ENDDO

       DO I=1,NVS
       DO K=1,3
       F(K,I)=F(K,I)+UINF(K,I)
       ENDDO
       ENDDO
       RLENG= dsqrt(RLENG)/RAD
       print *, ' LENG=', RLENG
c       read (*,*) icont

!      F-CALCULATION: INTERACTION BTWN THE DROP AND THE MOVING FRAME

       DO JP=1,NP

       DO I=INBH(JP),INFH(JP)

       DO K=1,3
       F(K,I)=0.D0
       ENDDO

       DO J=1,NVS
       P=0.D0
       RNX=0.D0

       DO K=1,3
       DX(K)=V(K,J)-V(K,I)
       P=P+DX(K)**2
       RNX=RNX+DX(K)*XNA(K,J)
       ENDDO

       R=DSQRT(P)

       DO K=1,3
       RDX(K)=XNA(K,J)/R
       ENDDO

       R3=R*P

       DO K=1,3
       RX1(K)=RNX*DX(K)/R3
       ENDDO

       DO K=1,3
       PX1(K)=(RX1(K)+RDX(K))*PHI(J)
       ENDDO

       IMI=IM(I,0)

       DO K=1,3
       F(K,I)=F(K,I)+PX1(K)*(CUR(J)-CUR(IMI))*FFACT
       ENDDO

       ENDDO
       ENDDO
       ENDDO

      IT=0
779   NRES=0
      IT=IT+1

      CALL BIOPER (.TRUE.,RMU,NVS,NVT,NP,INBH,INFH,INB,INF,V,X0,
     >E1,L1,E2,L2,XNW,IM,PHI,XNA,XCW,DM,SW,QW,INV,
     >IFRONT,IBACK,IUP,IDOWN,F,EA(1,1,1))

      CALL PRODMOD (EA(1,1,1), EA(1,1,1),NVT,PHI,RES0)
      RES0=DSQRT(RES0)

      DO I=1,NVT
      DO K=1,3
      EA(K,I,1)= EA(K,I,1)/RES0
      ENDDO
      ENDDO

777   DO I=1,NVT
      DO K=1,3
      QWIMP(K,I)= QW(K,I)
C     JUST FOR CONSISTENCY WITH OLD STRAIGHTFORWARD GMRES:
      IF (NRES.EQ.0) THEN
        QWIMP(K,I)= QWIMP(K,I)+RES0*EA(K,I,1)
      END IF
      DO J=1,NRES
      QWIMP(K,I)= QWIMP(K,I)+RLAM(J)*EA(K,I,J)
      ENDDO
      IF (NRES.EQ.NGMRES) THEN
	QW(K,I)= QWIMP(K,I)
      END IF
      ENDDO
      ENDDO

      RR=0.D0
      DO K=1,NRES+1
      P=0.D0
      DO J= MAX(K-1,1), NRES
      P=P+RLAM(J)*ACOEFF(J,K)
      ENDDO
      IF (K.EQ.1) P=P+RES0
      RCOEFF(K)=P
      RR=RR+P**2
      ENDDO

      CALL PRODMOD (QWIMP,QWIMP,NVT,PHI,SC)
      print *, ' RR=', rr
      ERRA=DSQRT(RR/SC)
      PRINT *, 'IT=', IT, ' ERRA=', ERRA, ' SC=', dsqrt(sc)

C     ALTERNATIVE METHOD, USING ERROR ESTIMATE ON DROP ONLY (GOOD OR NOT?)
      ERRA=0.D0
      SC=0.D0
      DO I=1,NVS
c      DO I=1,NVT
      ERR=0.D0
      P=0.D0
      DO K=1,3
      RKI=0.D0
      DO J=1,NRES+1
      RKI=RKI+RCOEFF(J)*EA(K,I,J)
      ENDDO
      ERR=ERR+RKI**2
      P=P+QWIMP(K,I)**2
      ENDDO
      ERRA=ERRA+ERR*PHI(I)
      SC=SC+P*PHI(I)
      ENDDO
      ERRA=DSQRT(ERRA/SC)

      PRINT *, 'IT=', IT, ' ERRA=', ERRA, ' SC=', dsqrt(sc)
      IF (ERRA.LT.TOL) GOTO 778

      IF (NRES.EQ.NGMRES) THEN
	GOTO 779
      ELSE
        NRES=NRES+1
      END IF
      IT=IT+1

      CALL BIOPER (.FALSE.,RMU,NVS,NVT,NP,INBH,INFH,INB,INF,
     >V,X0,E1,L1,E2,L2,
     >XNW,IM,PHI,XNA,XCW,DM,SW,EA(1,1,NRES),INV,
     >IFRONT,IBACK,IUP,IDOWN,F,EA(1,1,NRES+1))

      DO 555 J=1,NRES
      CALL PRODMOD (EA(1,1,NRES+1),EA(1,1,J),NVT,PHI,ACOEFF(NRES,J))

      DO I=1,NVT
      DO K=1,3
      EA(K,I,NRES+1)= EA(K,I,NRES+1)-ACOEFF(NRES,J)*EA(K,I,J)
      ENDDO
      ENDDO
555   CONTINUE
      CALL PRODMOD (EA(1,1,NRES+1),EA(1,1,NRES+1),NVT,PHI,P)
      P=DSQRT(P)
      ACOEFF(NRES,NRES+1)=P

      DO I=1,NVT
      DO K=1,3
      EA(K,I,NRES+1)= EA(K,I,NRES+1)/P
      ENDDO
      ENDDO

      DO I=1,NRES
      DO J=I,NRES
      TM(I,J)=0.D0
      DO K=1,I+1
      TM(I,J)= TM(I,J)+ACOEFF(I,K)*ACOEFF(J,K)
      ENDDO
      TM(J,I)=TM(I,J)
      ENDDO
      TM(I,NRES+1)= -ACOEFF(I,1)*RES0
      ENDDO

      CALL HFACT (NRES, TM, TM(1,NRES+1), RLAM)
c	CALL GMRES(NRES,NRES+1,TM)
      DO J=1,NRES
c     PRINT *, J, RLAM(J), TM(J,NRES+1)
C      RLAM(J)=TM(J,NRES+1)
      ENDDO

      GOTO 777

778   CONTINUE
      DO I=1,NVT
      DO K=1,3
      QW(K,I)=QWIMP(K,I)
      ENDDO
      ENDDO

C     CHECK:
      TST=0.D0
      DO I=1,NVS
      DO K=1,3
      TST=TST+ DABS(QW(K,I))*PHI(I)
      ENDDO
      ENDDO
c      PRINT *, ' tst=', tst
c      READ (*,*) ICONT

C      CALL DEFL(1,NVS,V,PHI,XNA,XP,DM(1,1,1),FS,QW,1,1,U,OM,DL)
C      PRINT *, ' U=', U
c------------------------------------------------------------
      DO I=1,NVS
      QM(I)=0.D0
      DO K=1,3
      U(K,I)=QW(K,I)
      QM(I)=QM(I)+QW(K,I)*XNA(K,I)
      ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE BIOPER (FIRST,RMU,NVS,NVT,NP,INBH,INFH,INB,INF,V,X0,
     >E1,L1,E2,L2,XNW,IM,PHI,XNA,XCW,DM,SW,QW,INV,
     >IFRONT,IBACK,IUP,IDOWN,F,RES)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'NTMAX.h'
      PARAMETER (NVSMAX= NTMAX/2+2,NVTMAX=120000)
      REAL*8 L1,L2
      DIMENSION PHI(*),F(3,*),RDX(3),RDY(3),RX1(3),RY1(3),
     >PX1(3),PY1(3),INB(0:*),INF(0:*)
      DIMENSION V(3,*),XNA(3,*),QW(3,*),IM(NVTMAX,0:*),
     >DL(3,NVTMAX),SX(3,NVTMAX),SY(3,NVTMAX),XCW(3),
     >INBH(0:*),INFH(0:*),INV(*),XNW(3,*),DM(3,3,2)
      DIMENSION L1(*), L2(*), X0(3,*),E1(3,*), E2(3,*)
      DIMENSION DX(3),U(3,NVSMAX),OM(3),SM(3),STR(3,3)
      DIMENSION RES(3,*), IFRONT(*),IBACK(*),IUP(*),IDOWN(*),
     >A1FIT(3,NVTMAX), A2FIT(3,NVTMAX),STR1(3,3),STR2(3,3)
      LOGICAL FIRST

      PI=DACOS(-1.D0)
      C2=3.D0/(2.D0*PI)
      DO I=1,NVT
      DO K=1,3
      DL(K,I)=0.D0
      ENDDO
      ENDDO

      IF (FIRST) THEN
        DO IP=0,NP
        DO I=INBH(IP),INFH(IP)
        DO K=1,3
        DL(K,I)=F(K,I)
        ENDDO
        ENDDO
        ENDDO
      END IF

      DO I=1,NVS
      DO K=1,3
      U(K,I)=QW(K,I)
      ENDDO
      ENDDO

C     CALCULATING LOCALLY FITTING LINEAR FORM TO QW: ---------------------
      DO IP=1,NP
      DO I=INBH(IP), INFH(IP)

      IF (IBACK(I).NE.0) THEN
	IB= IBACK(I)
      ELSE
	IB= I
      END IF

      IF (IFRONT(I).NE.0) THEN
	IFR= IFRONT(I)
      ELSE
	IFR= I
      END IF

      DX1=0.D0
      DO K=1,2
      DX1= DX1+ (V(K,IFR)-V(K,IB))*E1(K,IP)
      ENDDO
      DO K=1,3
      IF (IFR.EQ.IB) THEN
        A1FIT(K,I)=0.D0
      ELSE
        A1FIT(K,I)=(QW(K,IFR)-QW(K,IB))/DX1
      END IF
      ENDDO

      IF (IDOWN(I).NE.0) THEN
       ID= IDOWN(I)
      ELSE
       ID= I
      END IF

      IF (IUP(I).NE.0) THEN
	IU= IUP(I)
	DO K=1,3
        A2FIT(K,I)=(QW(K,IU)-QW(K,ID))/(V(3,IU)-V(3,ID))
        ENDDO
      ELSE
	DO K=1,2
C       YES!
	A2FIT(K,I)= -(QW(K,I)-QW(K,ID))/(V(3,I)+V(3,ID))
	ENDDO
	A2FIT(3,I)=(QW(3,I)+QW(3,ID))/(V(3,I)+V(3,ID))
      END IF

      do k=1,3
c      to return to older singularity subtraction, uncomment the next two lines
c      a1fit(k,i)=0.d0
c      a2fit(k,i)=0.d0
      enddo

      I1= INV(I)
      DO K=1,2
      A1FIT(K,I1)=  A1FIT(K,I)
      A2FIT(K,I1)= -A2FIT(K,I)
      ENDDO
      A1FIT(3,I1)= -A1FIT(3,I)
      A2FIT(3,I1)=  A2FIT(3,I)
      ENDDO
      ENDDO
C-------------------------------

      DO 330 JP=1,NP
      IF (INFH(JP).LT.INBH(JP)) GOTO 330
C     WALL-TO-WALL ADDBACK:
      DO 331 IP=1,NP
      IF (JP.EQ.IP) GOTO 331
      DO I=INBH(IP),INFH(IP)

      CALL RECSTRG (.FALSE.,V(1,I),X0(1,JP),E1(1,JP),L1(JP),
     >E2(1,JP),L2(JP),XNW(1,JP),STR,STR1,STR2)

      IMI= IM(I,JP)
      DO K=1,3
      DO L=1,3
      DL(K,I)=DL(K,I)+C2*STR(K,L)*QW(L,IMI)
      ENDDO
      ENDDO
      ENDDO
331   CONTINUE

C     WALL-TO-DROP ADDBACK:
      DO I=1,NVS
      CALL RECSTRG (.TRUE.,V(1,I),X0(1,JP),E1(1,JP),L1(JP),
     >E2(1,JP),L2(JP),XNW(1,JP),STR,STR1,STR2)

      DKSI1=0.D0
      DO K=1,2
      DKSI1= DKSI1+ (V(K,I)-V(K,IM(I,JP)))*E1(K,JP)
      ENDDO
      DKSI2= V(3,I)-V(3,IM(I,JP))
      IMI= IM(I,JP)
      DO K=1,3
      DO L=1,3
      DL(K,I)=DL(K,I)+C2*STR(K,L)*QW(L,IMI)+
     >C2*( (STR1(K,L)+STR(K,L)*DKSI1)*A1FIT(L,IMI) ) +
     >C2*( (STR2(K,L)+STR(K,L)*DKSI2)*A2FIT(L,IMI) )
      ENDDO
      ENDDO
      ENDDO
330   CONTINUE

!     ORDER
      CALL DEFL(NVS+1,NVT,V,PHI,XNA,XCW,DM(1,1,2),SW,QW,0,1,U,OM,DL)

C     SBOTTOM NOT IN USE ANYMORE
c!      SC= SW+ 2.D0*SBOTTOM
C      USE ONLY WITH 0 FOR NRBM
C      CALL DEFL(NVS+1,NVT,V,PHI,XNA,XCW,DM(1,1,2),SC,QW,0,1,U,OM,DL)

C      CALL DEFL(1,NVS,V,PHI,XNA,XP,DM(1,1,1),FS,QW,1,1,U,OM,DL)

C     NEW:
!     SELF-INTERACTION BETWEEN WALLS
34      DO IP=1,NP
C      WRITE (4,*) ' IP=', IP
C      PRINT *, ' IP=', IP
      DO JP=IP+1,NP
      DO I=INBH(IP), INFH(IP)
      DO K=1,3
      SY(K,I)=0.D0
      ENDDO
      ENDDO

      DO J=INBH(JP), INFH(JP)
      DO K=1,3
      SX(K,J)=0.D0
      ENDDO
      ENDDO

      DO I=INBH(IP), INFH(IP)
      IMI=IM(I,JP)
      DO J=INBH(JP), INFH(JP)
      IMJ=IM(J,IP)
      P=0.D0
      PX=0.D0
      PY=0.D0
      DO K=1,2
      DX(K)=V(K,J)-V(K,I)
      P=P+DX(K)**2
      PX=PX+ (QW(K,J)-QW(K,IMI))*DX(K)
      PY=PY+ (QW(K,I)-QW(K,IMJ))*DX(K)
      ENDDO
      R3M=V(3,J)-V(3,I)
      R3P=V(3,J)+V(3,I)
      R=P+R3M**2
      R5= 1.D0/(R**2*DSQRT(R))
      R=P+R3P**2
      R5ST= 1.D0/(R**2*DSQRT(R))

      RQX=  (PX+R3M*(QW(3,J)-QW(3,IMI)))*R5*PHI(J)
      RQXST=  (PX+R3P*(QW(3,J)+QW(3,IMI)))*R5ST*PHI(J)

      RQY= -(PY+R3M*(QW(3,I)-QW(3,IMJ)))*R5*PHI(I)
      RQYST= (-PY+R3P*(QW(3,I)+QW(3,IMJ)))*R5ST*PHI(I)

      RX=RQX+RQXST
      RY=RQY+RQYST
      DO K=1,2
      SY(K,I)=SY(K,I)+RX*DX(K)
      SX(K,J)=SX(K,J)+RY*DX(K)
      ENDDO
      SY(3,I)=SY(3,I)+RQX*R3M-RQXST*R3P
C     WAS WRONG:      SY(3,J)=SY(3,J)+RQY*R3M+RQYST*R3P
C     NOW CORRECT:
      SX(3,J)=SX(3,J)+RQY*R3M+RQYST*R3P
      ENDDO
      ENDDO

      IY=INBH(IP)
      IX=INBH(JP)

      DO I=INBH(IP),INFH(IP)
      RNX=0.D0
      DO K=1,3
      RNX=RNX+(V(K,IX)-V(K,I))*XNW(K,JP)
      ENDDO
      RNX=RNX*C2
      DO K=1,3
      DL(K,I)=DL(K,I)+ RNX*SY(K,I)
      ENDDO
      ENDDO

      DO J=INBH(JP),INFH(JP)
      RNY=0.D0
      DO K=1,3
      RNY=RNY-(V(K,IY)-V(K,J))*XNW(K,IP)
      ENDDO
      RNY=RNY*C2
      DO K=1,3
      DL(K,J)=DL(K,J)+ RNY*SX(K,J)
      ENDDO
      ENDDO
      ENDDO
      ENDDO

!     ************************************************************************
!     WALLS-TO-DROP CONTRIBUTION
	fmax=0.d0
      DO 10 JP=1,NP
      DO I=1,NVS
      DO K=1,3
      SY(K,I)=0.D0
      ENDDO
      ENDDO

      DO 20 J=INB(JP),INF(JP)
      DO 30 I=1,NVS
      IMI=IM(I,JP)
      QWRY=0.D0
      R=0.D0
      PX=0.D0
      DKSI1=0.D0
      DO K=1,2
      DKSI1=DKSI1+ (V(K,J)-V(K,IMI))*E1(K,JP)
      ENDDO
      DKSI2= V(3,J)-V(3,IMI)
      DO K=1,3
      DX(K)=V(K,J)-V(K,I)
      R=R+DX(K)**2
      PX=PX+DX(K)*(QW(K,J)-QW(K,IMI)-A1FIT(K,IMI)*DKSI1 -
     >A2FIT(K,IMI)*DKSI2)
      ENDDO
      RQX=PX*PHI(J)/(R**2*DSQRT(R))

      DO K=1,3
      SY(K,I)=SY(K,I)+RQX*DX(K)
      ENDDO
      fmax= max (fmax, dabs(rqx)*dsqrt(r) )
30     CONTINUE
20     CONTINUE

      IX=INBH(JP)
      DO I=1,NVS
      RNX=0.D0
      DO K=1,3
      RNX=RNX+(V(K,IX)-V(K,I))*XNW(K,JP)
      ENDDO
      RNX=RNX*C2
      DO K=1,3
      DL(K,I)=DL(K,I)+ RNX*SY(K,I)
      ENDDO
      ENDDO
10     CONTINUE
      print *, ' fmax=', fmax
c      read (*,*) icont

C     DROP TO MF CONTRIBUTION (U)

      UFACT=(RMU-1.D0)*3.D0/(4.D0*PI)

      DO JP=1,NP
      DO I=INBH(JP),INFH(JP)

      DO J=1,NVS
      P=0.D0
      RNX=0.D0
      UR=0.D0
      IMI=IM(I,0)

      DO K=1,3
      DX(K)=V(K,J)-V(K,I)
      P=P+DX(K)**2
      RNX=RNX+DX(K)*XNA(K,J)
      UR=UR+(U(K,J)-U(K,IMI))*DX(K)
      ENDDO

      UR=UFACT*UR/(DSQRT(P)*P**2)
      PX=UR*PHI(J)*RNX

      DO K=1,3
      DL(K,I)=DL(K,I)+PX*DX(K)
      ENDDO

      ENDDO
      ENDDO
      ENDDO

!     VELOCITY WHEN Y IS ON THE DROP

      DO I=1,NVS
      DO J=I+1,NVS
      P=0.D0
      UR=0.D0
      RNX=0.D0
      RNY=0.D0

      DO K=1,3
      DX(K)=V(K,J)-V(K,I)
      P=P+DX(K)**2
      UR=UR+(U(K,J)-U(K,I))*DX(K)
      RNX=RNX+DX(K)*XNA(K,J)
      RNY=RNY+DX(K)*XNA(K,I)
      ENDDO

      UR=UFACT*UR/(DSQRT(P)*P**2)
      PX=UR*PHI(J)*RNX
      PY=UR*PHI(I)*RNY

      DO K=1,3
      DL(K,I)=DL(K,I)+PX*DX(K)
      DL(K,J)=DL(K,J)+PY*DX(K)
      ENDDO
      ENDDO
      ENDDO

      DO I=1,NVS
      DO K=1,3
      DL(K,I)=DL(K,I)+U(K,I)*(RMU-1)/2
      ENDDO
      ENDDO

      DO IP=1,NP
      DO I=INBH(IP),INFH(IP)
      I1=INV(I)
!      PRINT *,  V(3,I), V(3,I1)
!      IF (MOD(I,10).EQ.1) READ (*,*)ICONT
      DL(1,I1)=  DL(1,I)
      DL(2,I1)=  DL(2,I)
      DL(3,I1)= -DL(3,I)
      ENDDO
      ENDDO

      DO I=1,NVT
      DO K=1,3
      QWNEW=DL(K,I)
      IF (I.GT.NVS) THEN
           QWNEW= -QWNEW
      ELSE
           QWNEW=QWNEW*2.D0/(RMU+1.D0)
      END IF

      RES(K,I)=QWNEW-QW(K,I)
      ENDDO
      ENDDO
      RETURN
      END


      SUBROUTINE DEFL(INB,INF,V,PHI,XNA,XC,DM,FS,QW,MRBM,MNORM,UAV,
     >OM,DL)
!	REAL SET TO ORIGINAL WITH SPACE
      IMPLICIT REAL *8(A-H,O-Z)
      DIMENSION QW(3,*),DL(3,*),XC(3),DM(3,3),PHI(*),XNA(3,*),
     >V(3,*),SM(3),DX(3),OM(3),UAV(3)

      DO K=1,3
      SM(K)=0.D0
      UAV(K)=0.D0
      ENDDO
      QN=0.D0
      DO 502 I=INB,INF
      P=0.D0
      DO K=1,3
      UAV(K)=UAV(K)+QW(K,I)*PHI(I)
      DX(K)=V(K,I)-XC(K)
      P=P+QW(K,I)*XNA(K,I)
      ENDDO
      QN=QN+P*PHI(I)
      SM(1)=SM(1)+(DX(2)*QW(3,I)-DX(3)*QW(2,I))*PHI(I)
      SM(2)=SM(2)+(DX(3)*QW(1,I)-DX(1)*QW(3,I))*PHI(I)
      SM(3)=SM(3)+(DX(1)*QW(2,I)-DX(2)*QW(1,I))*PHI(I)
502   CONTINUE
      QN=MNORM*QN/FS
      DO K=1,3
      UAV(K)=UAV(K)/FS
      OM(K)=0.D0
      DO L=1,3
      OM(K)=OM(K)+DM(K,L)*SM(L)
      ENDDO
      ENDDO

      DO 503 I=INB,INF
      DO K=1,3
      DX(K)=(V(K,I)-XC(K))
      ENDDO
      DL(1,I)= DL(1,I) -(UAV(1)+OM(2)*DX(3)-OM(3)*DX(2))*MRBM +
     >QN*XNA(1,I)
      DL(2,I)= DL(2,I) -(UAV(2)+OM(3)*DX(1)-OM(1)*DX(3))*MRBM +
     >QN*XNA(2,I)
      DL(3,I)= DL(3,I) -(UAV(3)+OM(1)*DX(2)-OM(2)*DX(1))*MRBM +
     >QN*XNA(3,I)
503   CONTINUE
      RETURN
      END

      SUBROUTINE DEFLMAT (INB,INF,XC,V,PHI,DM)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XC(3),V(3,*),PHI(*),T(3,6),DM(3,3)
      DO I=1,3
      DO J=1,6
       T(I,J)=0.D0
      ENDDO
       T(I,I+3)=1.D0
      ENDDO
      DO I=INB,INF
      R2=(V(1,I)-XC(1))**2 + (V(2,I)-XC(2))**2 + (V(3,I)-XC(3))**2
      DO K=1,3
      DO L=1,3
      DEL=0.D0
      IF(K.EQ.L) DEL=1.D0
      T(K,L)=T(K,L)+(DEL*R2-(V(K,I)-XC(K))*(V(L,I)-XC(L)))*PHI(I)
      ENDDO
      ENDDO
      ENDDO
      CALL GM36(3,6,T)
      DO K=1,3
      DO L=1,3
      DM(K,L)=T(K,L+3)
      ENDDO
      ENDDO
      RETURN
      END

      SUBROUTINE GM36(N,M,T)
C!	CHANGE REAL TO MATCH ORIGINAL
      REAL*8 T(3,6),R,U
      INTEGER S
      DO 1 I=1,N
      S=I
      R=T(I,I)
      DO 2 K=I+1,N
      IF(DABS(T(K,I)).GT.DABS(R)) THEN
         S=K
	 R=T(K,I)
      END IF
2     CONTINUE
      T(S,I)=T(I,I)
      DO 3 J=I+1,M
      U=T(S,J)/R
      T(S,J)=T(I,J)
      T(I,J)=U
      DO 4 K=I+1,N
4     T(K,J)=T(K,J)-T(K,I)*U
3     CONTINUE
1     CONTINUE
      DO 5 I=N-1,1,-1
      DO 5 J=N+1,M
      DO 5 K=I+1,N
      R=T(I,K)
      T(I,J)=T(I,J)-R*T(K,J)
5     CONTINUE
      RETURN
      END

      SUBROUTINE PNLSYMAA(I,NVTMAX,INB1,INF1,V1,DS1,XP,C,RAD,ALPHA,
     >L2,X0,E1,INB,INF,INBH,INFH,INV,V,XNW,PHI,IFRONT,IBACK,IUP,IDOWN)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL *8 L2,L2H
      PARAMETER (NVCONTMAX=3000)
      DIMENSION X0(3),E1(3),E2(3),V(3,*),XNW(3),INV(*),PHI(*),
     >V1(2,*), DS1(*),XP(3),X(3),DS2(NVCONTMAX),S2(NVCONTMAX),
     >IFRONT(*), IBACK(*),IUP(*),IDOWN(*)
!     ESSENTIALLY ASSUMED HERE THAT BASIS VECTOR E2 OF THE PANEL IS (0,0,1)
!     AND SYMMETRY IS ABOUT Z=0 PLANE, AND CORNER X0(3)= -L2/2
      E2(1)=0
      E2(2)=0
      E2(3)=1.D0
      N1= INF1-INB1+1
      XNW(1)=E1(2)*E2(3)-E1(3)*E2(2)
      XNW(2)=E1(3)*E2(1)-E1(1)*E2(3)
      XNW(3)=E1(1)*E2(2)-E1(2)*E2(1)
      L2H=0.5D0*L2
      I0=I
      I=0
      SN2H=0.D0
      N2HMAX=0
      DO 44 I1 = INB1,INF1
      INEW=0
      S=0.D0
      DO K=1,2
      X(K)=V1(K,I1)
      ENDDO
      X(3)=X0(3)
45    RR=0.D0
      DO K=1,3
      RR=RR+(X(K)-XP(K))**2
      ENDDO
      GAP=1.D0
      H=C*GAP**ALPHA
      S=S+H
      IF (S.GT.L2H) THEN
        S=S-H
        DO J=1,INEW
        RAT=L2H/S
        DS2(J)=DS2(J)*RAT
        DO K=1,2
        V(K,I0+I+J)=V1(K,I1)
        ENDDO
        V(3,I0+I+J)=X0(3)+RAT*S2(J)+0.5D0*DS2(J)
        PHI(I0+I+J)= DS1(I1)*DS2(J)
	M= I0+I+J

	IF (I1.NE.INB1) THEN
	  IBACK(M)=M-INEW
	ELSE
	  IBACK(M)= 0
	END IF

	IF (I1.NE.INF1) THEN
	  IFRONT(M)=M+INEW
	ELSE
	  IFRONT(M)= 0
	END IF

	IF (J.NE.1) THEN
	  IDOWN(M)= M-1
	ELSE
	  IDOWN(M)=0
	END IF

	IF (J.NE.INEW) THEN
	  IUP(M)=M+1
	ELSE
	  IUP(M)=0
	END IF

	ENDDO
        GOTO 46

      ELSE
        INEW=INEW+1
        IF (INEW.GT. NVCONTMAX) THEN
          PRINT *, ' NVCONTMAX TOO SMALL'
          STOP
        END IF

        IF (I0 + 2*(I+INEW).GT. NVTMAX) THEN
          PRINT *, ' NVTMAX TOO SMALL'
          STOP
        END IF
        DS2(INEW)=H
        S2(INEW)=S-H
        X(3)=X(3)+H
        GOTO 45
      END IF

46    N2HMAX=MAX(N2HMAX,INEW)
      SN2H=SN2H+INEW
      I=I+INEW
44    CONTINUE
      INB=I0+1
      INF=I0+2*I
      INBH=INB
      INFH=I0+I

      DO J=I0+1,I0+I
      INV(J)=J+I
      DO K=1,2
      V(K,INV(J))=V(K,J)
      ENDDO
      V(3,INV(J)) = -V(3,J)
      PHI(INV(J))=PHI(J)
      ENDDO
      I=INF
      AVN2H=SN2H/N1
      WRITE (4,*) ' N2HMAX=', N2HMAX, ' AVN2H=', AVN2H
      PRINT *, ' N2HMAX=', N2HMAX, ' AVN2H=', AVN2H
      RETURN
      END

      SUBROUTINE SPHAREA(Y,X0,DSSPH)
!     CALCULATES THE AREA OF A SPHERICAL TRIANGLE PROJECTED ON UNIT SPHERE
!     Y(J, J=1,2,3)IS THE SPHERE CENTER, AND TRIANGLE VERTICES
!     BEFORE PROJECTION ARE X0(K,J
      IMPLICIT REAL *8 (A-H,O-Z)
      DIMENSION XNA(3,3),SIDE1(3),SIDE2(3),Y(3),X0(3,*)
!     FOR THIS APPLICATION, X0(3,*) INSTEAD OF X0(3,3) IS IMPORTANT
      DO K=1,3
      P=0.D0
      DO J=1,3
      XNA(J,K)=X0(J,K)-Y(J)
      P=P+XNA(J,K)**2
      ENDDO
      P=1.D0/DSQRT(P)
      DO J=1,3
      XNA(J,K)=P*XNA(J,K)
      ENDDO
      ENDDO

      PI=DACOS(-1.D0)
      SUM=0.D0

      DO 446 K=1,3
      K1=K+1
      IF (K1.GT.3) K1=K1-3
      K2=K+2
      IF (K2.GT.3) K2=K2-3
      P1=0.D0
      P2=0.D0
      DO J=1,3
      SIDE1(J)=XNA(J,K1)-XNA(J,K)
      SIDE2(J)=XNA(J,K2)-XNA(J,K)
      P1=P1+ SIDE1(J)*XNA(J,K)
      P2=P2+ SIDE2(J)*XNA(J,K)
      ENDDO
      A1=0.D0
      A2=0.D0
      P=0.D0
      DO J=1,3
      SIDE1(J)=SIDE1(J)-P1*XNA(J,K)
      A1=A1+SIDE1(J)**2
      SIDE2(J)=SIDE2(J)-P2*XNA(J,K)
      A2=A2+SIDE2(J)**2
      P=P+SIDE1(J)*SIDE2(J)
      ENDDO
      P=P/DSQRT(A1*A2)
      P=P*(1.D0-1.D-15)
      SUM=SUM+DACOS(P)
446   CONTINUE
      DSSPH= SUM-PI
      RETURN
      END

      SUBROUTINE RECSTREX (Y,X0,E1,L1,E2,L2,XN,STR)
      IMPLICIT REAL *8 (A-H,O-Z)
      DIMENSION Y(3),X0(3),E1(3),E2(3),XN(3),STR(3,3),E(3,4),X(3,4),
     >L(4)
      REAL *8 L1,L2,L

      DO I=1,3
      DO J=1,3
      STR(I,J)=0.D0
      ENDDO
      ENDDO

      DO J=1,3
      X(J,1)=X0(J)
      E(J,1)=E1(J)

      X(J,2)=X0(J)+L1*E1(J)
      E(J,2)=E2(J)

      X(J,3)=X(J,2)+L2*E2(J)
      E(J,3)= -E1(J)

      X(J,4)=X0(J)+L2*E2(J)
      E(J,4)= -E2(J)
      ENDDO

      L(1)=L1
      L(2)=L2
      L(3)=L1
      L(4)=L2

      DO K=1,4
      CALL SEGCONT (Y,X(1,K),E(1,K),L(K),STR)
      ENDDO

      P=XN(1)*(E1(2)*E2(3)-E1(3)*E2(2))+XN(2)*
     >(E1(3)*E2(1)-E1(1)*E2(3))+ XN(3)*(E1(1)*E2(2)-E1(2)*E2(1))
      IF (P.LT.0.D0) THEN
        DO I=1,3
        DO J=1,3
	STR(I,J)=-STR(I,J)
	ENDDO
	ENDDO
      END IF

      RNX=0.D0
      DO J=1,3
      RNX=RNX+XN(J)*(X0(J)-Y(J))
      ENDDO

      CALL SPHAREA(Y,X,S1)
      DO J=1,3
      X(J,2)=X(J,3)
      X(J,3)=X(J,4)
      ENDDO
      CALL SPHAREA(Y,X,S2)
      S=(S1+S2)/3.D0
      IF(RNX.LT.0.D0) S= -S
      DO I=1,3
      STR(I,I)=STR(I,I)+S
      ENDDO
      RETURN
      END



      SUBROUTINE SEGCONT (Y,X0,DIR,L,STR)
      IMPLICIT REAL *8 (A-H,O-Z)
      REAL *8 L
      DIMENSION Y(3),X0(3),DIR(3),STR(3,3),DX(3),E(3,3)
      RR=0.D0
      P=0.D0
      DO K=1,3
      DX(K)=Y(K)-X0(K)
      RR=RR+DX(K)**2
      P=P+DX(K)*DIR(K)
      ENDDO
      H=DSQRT(RR-P**2)
      T1=-P/H
      T2= (L-P)/H
      RI1=T2/DSQRT(1.D0+T2**2) - T1/DSQRT(1.D0+T1**2)
      RI2=1.D0/DSQRT(1.D0+T1**2) -1.D0/DSQRT(1.D0+T2**2)
      RI1=RI1/3.D0
      RI2=RI2/3.D0
      E(1,3)=DIR(2)*DX(3)-DIR(3)*DX(2)
      E(2,3)=DIR(3)*DX(1)-DIR(1)*DX(3)
      E(3,3)=DIR(1)*DX(2)-DIR(2)*DX(1)
      TEMP=0.D0
      DO K=1,3
      TEMP=TEMP+E(K,3)**2
      E(K,2)=DX(K)-P*DIR(K)
      E(K,1)=DIR(K)
      ENDDO
      TEMP=1.D0/DSQRT(TEMP)
      DO K=1,3
      E(K,2)=E(K,2)/H
      E(K,3)=E(K,3)*TEMP
      ENDDO
      DO K=1,3
      DO J=1,3
      STR(K,J)=STR(K,J)+(RI2*E(J,1)-RI1*E(J,2))*E(K,3)
      ENDDO
      ENDDO
      RETURN
      END

	SUBROUTINE VELOCITY(Y,U)
	IMPLICIT REAL*8 (A-H,O-Z)
	INCLUDE 'NPCH.h'
	PARAMETER (NP=NPCH,NVTMAX=NVTMAXCH)
	COMMON / BIVEL / X0,L1,E1,V,Q,DSW,XNW,INB,INF,NVT
	DIMENSION Y(*), U(*), X0(2,NP),L1(NP),E1(2,NP),
     >	V(2,NVTMAX), Q(2,NVTMAX),STR(2,2),
     >INB(NP),INF(NP),XNW(2,NP),DSW(NP),SUM(2),DX(2),X1(2)
	REAL *8 L1

225	PI=DACOS(-1.D0)
	C2=-2.D0/PI
	DO K=1,2
	U(K)=0.D0
	ENDDO
	DO JP=1,NP
	RMIN=1.D20
	DO I=INB(JP),INF(JP)
	R=0.D0
	DO K=1,2
	R=R+(V(K,I)-Y(K))**2
	ENDDO
	IF(R.LT.RMIN) THEN
	  RMIN=R
	  JMIN=I
	ENDIF
	ENDDO
	DO K=1,2
	X1(K)=X0(K,JP)+L1(JP)*E1(K,JP)
	ENDDO
	CALL LINESTR (Y,X0(1,JP),X1,XNW(1,JP),STR)
	DO K=1,2
	DO L=1,2
	U(K)=U(K)-C2*STR(K,L)*Q(L,JMIN)
        ENDDO
        ENDDO

	 DO K=1,2
	 SUM(K)=0.D0
	 ENDDO
	 DO J=INB(JP),INF(JP)
	  P=0.D0
	  QR=0.D0
	 DO K=1,2
	  DX(K)=V(K,J)-Y(K)
	  P=P+DX(K)**2
	  QR=QR+(Q(K,J)-Q(K,JMIN))*DX(K)
	 ENDDO
	 QR=QR/(P**2)
	 DO K=1,2
	  SUM(K)=SUM(K)+QR*DX(K)
	 ENDDO
	 ENDDO
	 IX=INB(JP)
	 RNX=0.D0
	 DO K=1,2
	 RNX=RNX+(V(K,IX)-Y(K))*XNW(K,JP)
	 ENDDO
	 RNX=-RNX*DSW(JP)*C2
	 DO K=1,2
	 U(K)=U(K)+RNX*SUM(K)
	 ENDDO
	ENDDO
	RETURN
	END
!
      SUBROUTINE LINESTR (Y,X1,X2,XN,STR)
      IMPLICIT REAL *8 (A-H,O-Z)
      DIMENSION Y(2),X1(2),X2(2),XA(2),XB(2),XN(2),E(2,2),STR(2,2),
     >DPSI(2,2)
      P=0.D0
      DO K=1,2
      XA(K)=X1(K)-Y(K)
      XB(K)=X2(K)-Y(K)
      P=P+XA(K)*XN(K)
      E(K,1)=XN(K)
      ENDDO
      SIG=1.D0
      IF (P.LT.0.D0) THEN
	SIG= -SIG
	E(1,1)= -E(1,1)
        E(2,1)= -E(2,1)
      END IF
      P=XA(1)*XB(2)-XA(2)*XB(1)
      IF (P.LT.0.D0) THEN
        DO K=1,2
        TEMP=XA(K)
        XA(K)=XB(K)
        XB(K)=TEMP
        ENDDO
      END IF

      AL=0.D0
      BL=0.D0
      DO K=1,2
      AL=AL+XA(K)**2
      BL=BL+XB(K)**2
      ENDDO

      SI=(E(1,1)*XA(2)-E(2,1)*XA(1))/DSQRT(AL)
      IF (SI.GT.1.D0) SI=1.D0
      IF (SI.LT.-1.D0) SI=-1.D0
      TA=DASIN(SI)
      SI=(E(1,1)*XB(2)-E(2,1)*XB(1))/DSQRT(BL)
      IF (SI.GT.1.D0) SI=1.D0
      IF (SI.LT.-1.D0) SI=-1.D0

      TB=DASIN(SI)
      DPSI(1,1)=0.5D0*(TB-TA+0.5D0*(DSIN(2.D0*TB)-DSIN(2.D0*TA)))
      DPSI(2,2)=TB-TA-DPSI(1,1)
      DPSI(1,2)= 0.25D0*(DCOS(2.D0*TA)-DCOS(2.D0*TB))
      DPSI(2,1)=DPSI(1,2)

      E(1,2)= -E(2,1)
      E(2,2)= E(1,1)

      DO I=1,2
      DO K=1,2
      SUM=0.D0
      DO J=1,2
      DO L=1,2
      SUM=SUM+E(I,J)*E(K,L)*DPSI(J,L)
      ENDDO
      ENDDO
      STR(I,K)=SIG*SUM
      ENDDO
      ENDDO
      RETURN
      END
      SUBROUTINE CONTOR (XP,A,X0,E1,L1,NPC,X0C,SOL)
C     A FIXED VERSION OF CONTOR, TO ALLOW FOR DOUBLE-CONNECTED INTERSECTION
C     AREA OF THE SQUARE AND CHANNEL;IF THIS HAPPENS, CHOOSES THE PIECE
C     CONTAINING PARTICLE CENTER.
C       DID NOT TRY YET TO FIX SOL; NOT USED SO FAR, ANYWAY.
!     ************************************************************************
!     * XP(1),XP(2) PARTICLE CENTER LOCATION                                 *
!     * A IS HALFWIDTH OF THE SQUARE BOX (NOT TRIMMED YET) AROUND XP(1),XP(2)*
!     * X0, E1, L1 ARE ARRAYS SPECIFYING THE CHANNEL (EXACTLY LIKE IN BINEW  *
!     * AND TRAJNEW.                                                         *
!     *                                                                      *
!     * OUTPUTS:                                                             *
!     *      NPC IS THE NUMBER OF CORNER POINTS OF THE CELL (BUT THE         *
!     *      FIRST AND LAST POINT COUNTED TWICE)                             *
!     *                                                                      *
!     *     X0C IS ARRAYS OF CORNER POINTS FOR THE CELL                      *
!     *     (LISTED COUNTERCLOCKWISE)                                        *
!     *     SOL(IP) IS .TRUE. IF THE CELL EDGE ORIGINATING                   *
!     *     FROM X0(1,IP),X0(2,IP)                                           *
!     *     BELONG TO THE CHANNEL BOUNDARY, AND SOL(IP)=.FALSE. OTHERWISE.   *
!     ************************************************************************

      IMPLICIT REAL *8 (A-H,O-Z)
      INCLUDE 'NPCH.h'
      PARAMETER (NPCMAX=25, NVT=50000)
      REAL *8 L1
      LOGICAL SOL
!     XP(*) ESSENTIAL;
      DIMENSION XP(*),L1(*), X0(2,*), E1(2,*),V(2,NVT),X0SQ(2,4),
     >INB(NPCH),INF(NPCH),X0INT(2),X0C(2,*),SOL(*), X(2,4),Y(2)
      INTEGER COLOR(4)
      IF (INDF(XP,NPCH,X0).EQ.0) THEN
	PRINT *, ' PARTICLE CENTER OUTSIDE CHANNEL'
      END IF
      EPS=1.D-14*A**2
      DO IP=1,NPCMAX
      SOL(IP)=.FALSE.
      ENDDO

      S=0.D0
      DO IP=1,NPCH
      S=S+L1(IP)
      ENDDO

      I=0
      DO IP=1,NPCH
      N1=L1(IP)*NVT/S
      H1=L1(IP)/N1
      INB(IP)=I + 1
      INF(IP) = I + N1
      DO I1 = 1,N1
      I = I + 1
      DO K=1,2
      V(K,I)=X0(K,IP)+(I1-1)*H1*E1(K,IP)
      ENDDO
      ENDDO
      ENDDO
!     ************************************************************************
C      OPEN (2,FILE='V.DAT')
C      DO IP=1,NPCH
C      DO I=INB(IP),INF(IP)
C      WRITE (2,*) V(1,I),V(2,I)
C      ENDDO
C      ENDDO
C      CLOSE(2)
!      V OK
!     ************************************************************************

      X0SQ(1,1)=XP(1)-A
      X0SQ(2,1)=XP(2)-A

      X0SQ(1,2)=XP(1)+A
      X0SQ(2,2)=XP(2)-A

      X0SQ(1,3)=XP(1)+A
      X0SQ(2,3)=XP(2)+A

      X0SQ(1,4)=XP(1)-A
      X0SQ(2,4)=XP(2)+A

      OPEN (2,FILE='SQ.DAT')
      DO IP=1,4
      WRITE (2,*) X0SQ(1,IP),X0SQ(2,IP)
      ENDDO
      WRITE (2,*) X0SQ(1,1),X0SQ(2,1)
      CLOSE(2)

      IPSTART=1

51    IP=IPSTART

      NPC=0
      IPSQ0=0
      DO JP=1,4
      COLOR(JP)=0
      ENDDO
      IND=1

2     DO 1 I=INB(IP),INF(IP)
      I1=I+1
      IF (I1.GT.INF(NPCH))I1=1
      SIGOLD= -1.D20
      SIGNEW= -1.D20
      DO K=1,2
      SIGOLD=MAX(SIGOLD, DABS(V(K,I)-XP(K)) -A)
      SIGNEW=MAX(SIGNEW, DABS(V(K,I1)-XP(K)) -A)
      ENDDO
      IF (SIGNEW*SIGOLD.GT.0.D0) GOTO 1
      IF( (DABS(V(1,I)-XP(1))-A)*(DABS(V(1,I1)-XP(1))-A).LT.0.D0 ) THEN
	  K=1
      ELSE
        K=2
      END IF

      IF ( (V(K,I)-XP(K)-A)*(V(K,I1)-XP(K)-A).LT.0.D0 ) THEN
	  X0INT(K)=XP(K)+A
	  IF (K.EQ.1) THEN
	    IPSQ=2
	  ELSE
	    IPSQ=3
	  END IF
      ELSE
        X0INT(K)=XP(K)-A
	  IF (K.EQ.1) THEN
	    IPSQ=4
	  ELSE
	    IPSQ=1
	  END IF
      END IF

      M=3-K
      X0INT(M)=V(M,I)+(V(M,I1)-V(M,I))*(X0INT(K)-V(K,I))/
     >(V(K,I1)-V(K,I))
C      PRINT *, ' IP=', IP, ' K=', K, ' IPSQ=', IPSQ, ' SIGNEW=',SIGNEW
C      PRINT *, ' SIGOLD=', SIGOLD
C      PRINT *, ' IND=', IND, ' IPSQ0=', IPSQ0
C      PRINT *, ' NPC=', NPC
C      READ (*,*)  ICONT
      IF (SIGNEW.GT.0.D0 .AND. IND.EQ.0) GOTO 1
      IF (SIGNEW.GT.0.D0 .AND. IND.EQ.1) THEN
        IPSQ0=IPSQ
	IP0=IP
	NPC=NPC+1
        IF (NPC.GT.NPCMAX) THEN
	    PRINT *, ' NPCMAX TOO SMALL'
	    STOP
        END IF
	DIFF=0.D0
        DO K=1,2
        X0C(K,NPC)=X0INT(K)
        DIFF=DIFF+(X0C(K,NPC)-X0C(K,1))**2
	ENDDO
        IF (NPC.NE.1 .AND. DIFF.LT.EPS) GOTO 50
      ELSE

        IF(IPSQ0.NE.0) THEN
	    DO K=1,2
	    Y(K)=0.5D0*(X0C(K,NPC)+X0INT(K))
	    ENDDO
	    IND= INDF(Y,NPCH,X0)
	    IF(IND.EQ.0  .AND.  MOD(IP-IP0-1,NPCH).EQ.0) GOTO 1
	    IF(IPSQ.EQ.IPSQ0  .AND. IND.EQ.0) IPSQ=IPSQ+4
	    IF(IPSQ.LT.IPSQ0) IPSQ=IPSQ+4
            IND=1
	    DO 44 J=IPSQ0,IPSQ
	    J1=MOD(J-1,4)+1

          IF (J.EQ.IPSQ0) GOTO 45
	  NPC=NPC+1
          IF (NPC.GT.NPCMAX) THEN
	      PRINT *, ' NPCMAX TOO SMALL'
	      STOP
          END IF
	    DIFF=0.D0
	    DO K=1,2
	    X0C(K,NPC)=X0SQ(K,J1)
          DIFF=DIFF+(X0C(K,NPC)-X0C(K,1))**2
	    ENDDO
            IF (NPC.NE.1 .AND. DIFF.LT.EPS) GOTO 50
45          IF (COLOR(J1).EQ.1) THEN
             DO K=1,2
             Y(K)=0.5D0*(X0C(K,NPC)+X(K,J1))
             ENDDO
             IF(INDF(Y,NPCH,X0).EQ.1)  THEN
              NPC=NPC+1
              IF (NPC.GT.NPCMAX) THEN
                PRINT *, ' NPCMAX TOO SMALL'
                STOP
              END IF

              DO K=1,2
              X0C(K,NPC)=X(K,J1)
              ENDDO
              RETURN
             END IF
            END IF
44          CONTINUE

	END IF

	NPC=NPC+1
        IF (NPC.GT.NPCMAX) THEN
          PRINT *, ' NPCMAX TOO SMALL'
          STOP
        END IF
	IF (IPSQ.GT.4) IPSQ=IPSQ-4
	DIFF=0.D0
        DO K=1,2
        X0C(K,NPC)=X0INT(K)
        DIFF=DIFF+(X0C(K,NPC)-X0C(K,1))**2
	X(K,IPSQ)=X0C(K,NPC)
	ENDDO
	COLOR(IPSQ)=1
        SOL(NPC)=.TRUE.
        IF (NPC.NE.1 .AND. DIFF.LT.EPS) GOTO 50
      END IF
1     CONTINUE

      IF (SIGNEW.LT.0.D0) THEN
        NPC=NPC+1
        IF (NPC.GT.NPCMAX) THEN
          PRINT *, ' NPCMAX TOO SMALL'
          STOP
        END IF
        DIFF=0.D0
        DO K=1,2
        X0C(K,NPC)=V(K,I1)
        DIFF=DIFF+(X0C(K,NPC)-X0C(K,1))**2
        ENDDO
        SOL(NPC)=.TRUE.
        IF (NPC.NE.1 .AND. DIFF.LT.EPS) GOTO 50
      END IF

      IP=IP+1
      IF (IP.GT.NPCH) IP=1
      IF (IP.NE.IPSTART. OR. NPC.NE.0) GOTO 2

      DO K=1,2
      DO IPC=1,4
      X0C(K,IPC)=X0SQ(K,IPC)
      ENDDO
      X0C(K,5)=X0SQ(K,1)
      ENDDO
      NPC=5
      RETURN

50    IF (INDF(XP,NPC-1,X0C).EQ.0) THEN
        IPSTART=IPSTART+1
        IF (IPSTART.GT.NPCH) IPSTART=1
	GOTO 51
      END IF
      RETURN
      END

      INTEGER FUNCTION INDF(Y,NPCH,X0)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Y(2),X0(2,*),XA(2),XB(2)
      PI=DACOS(-1.D0)
      SUM=0.D0
      DO IP=1,NPCH
      IP2=IP+1
      IF (IP2.GT.NPCH) IP2=IP2-NPCH
      P=0.D0
      A=0.D0
      B=0.D0
      C=0.D0
      DO K=1,2
      XA(K)=X0(K,IP)-Y(K)
      XB(K)=X0(K,IP2)-Y(K)
      A=A+XA(K)**2
      B=B+XB(K)**2
      C=C+ (XB(K)-XA(K))**2
      P=P+XA(K)*XB(K)
      ENDDO
      PSI= (A+B-C)/(2.D0*DSQRT(A*B))
      PSI=DACOS(PSI)
      IF ( XA(1)*XB(2)-XA(2)*XB(1) .LT. 0.D0) PSI= -PSI
      SUM=SUM+PSI/(2.D0*PI)
      ENDDO

      IF(DABS(SUM).GT. 1.D-3 .AND. DABS(SUM-1.D0).GT. 1.D-3) THEN
	PRINT  *, ' ERROR IN INDF'
	STOP
      END IF
CC      PRINT *, ' SUM=',SUM
      IF (DABS(SUM).LT.1.D-3) THEN
	INDF=0
      ELSE
	INDF=1
      END IF
      RETURN
      END


      SUBROUTINE PRODMOD (X,Y,NVT,PHI,PROD)
      IMPLICIT REAL *8(A-H,O-Z)
      DIMENSION X(3,*),Y(3,*),PHI(*)
      PROD=0.D0
      DO I=1,NVT
      P=0.D0
      DO K=1,3
      P=P+X(K,I)*Y(K,I)
      ENDDO
      PROD=PROD+P*PHI(I)
      ENDDO
      RETURN
      END

      SUBROUTINE GMRES(N,M,T)
      INCLUDE 'NRESMAX.h'
      real *8 T(NRESMAX,NRESMAX+1),R,U
      INTEGER S
      DO 1 I=1,N
      S=I
      R=T(I,I)
      DO 2 K=I+1,N
      IF(DABS(T(K,I)).GT.DABS(R)) THEN
         S=K
	 R=T(K,I)
      END IF
2     CONTINUE
      T(S,I)=T(I,I)
      DO 3 J=I+1,M
      U=T(S,J)/R
      T(S,J)=T(I,J)
      T(I,J)=U
      DO 4 K=I+1,N
4     T(K,J)=T(K,J)-T(K,I)*U
3     CONTINUE
1     CONTINUE
      DO 5 I=N-1,1,-1
      DO 5 J=N+1,M
      DO 5 K=I+1,N
      R=T(I,K)
      T(I,J)=T(I,J)-R*T(K,J)
5     CONTINUE
      RETURN
      END

      SUBROUTINE HFACT (N,A,Y,X)
      IMPLICIT REAL *8( A-H,O-Z)
      INCLUDE 'NRESMAX.h'
      DIMENSION A(NRESMAX,*),B(NRESMAX,NRESMAX),Z(NRESMAX),
     >X(*),Y(*)
      DO I=1,N
      S=A(I,I)
      DO K=1,I-1
      S=S - B(K,I)**2
      ENDDO
      IF (S.LT.0.D0) THEN
	OPEN (17, FILE='FEND')
	PRINT *, ' ARG SQRT<0 IN HFACT: S=', S, ' I=', I
	WRITE (17,*) ' ARG SQRT<0 IN HFACT: S=', S, ' I=', I
	STOP
      END IF
      B(I,I)=DSQRT(S)

      DO J=I+1,N
      S=A(I,J)
      DO K=1,I-1
      S=S-B(K,I)*B(K,J)
      ENDDO
      B(I,J)= S/B(I,I)
      ENDDO
      ENDDO

C     SOLVER:
      DO I=1,N
      S=Y(I)
      DO J=1,I-1
      S=S-B(J,I)*Z(J)
      ENDDO
      Z(I)=S/B(I,I)
      ENDDO

      DO I=N,1, -1
      S=Z(I)
      DO J=I+1,N
      S= S -B(I,J)*X(J)
      ENDDO
      X(I)=S/B(I,I)
      ENDDO
      RETURN
      END

C     EXTENSION OF RECSTREX: IF ext=.TRUE., CALCULATES ADDITIONAL
C     TENSOR INTEGRALS STR1 AND STR2 (ANALYTICALLY) NEEDED FOR
C     HIGH-ORDER SINGULARITY SUBTRACTION
      SUBROUTINE RECSTRG (EXT,Y,X0,E1,L1,E2,L2,XN,STR,STR1,STR2)
      IMPLICIT REAL *8 (A-H,O-Z)
      DIMENSION Y(3),X0(3),E1(3),E2(3),XN(3),STR(3,3),E(3,4),X(3,4),
     >L(4), STR1(3,3), STR2(3,3),XM(3),R(3,5),RLT(3),D(5),H(3),
     >RR_R3L(3,3,4),R_R3S(3)

      REAL *8 L1,L2,L
      LOGICAL EXT
      DO I=1,3
      DO J=1,3
      STR(I,J)=0.D0
      ENDDO
      ENDDO

      DO J=1,3
      X(J,1)=X0(J)
      E(J,1)=E1(J)

      X(J,2)=X0(J)+L1*E1(J)
      E(J,2)=E2(J)

      X(J,3)=X(J,2)+L2*E2(J)
      E(J,3)= -E1(J)

      X(J,4)=X0(J)+L2*E2(J)
      E(J,4)= -E2(J)
      ENDDO

      L(1)=L1
      L(2)=L2
      L(3)=L1
      L(4)=L2

      DO K=1,4
      CALL SEGCONT (Y,X(1,K),E(1,K),L(K),STR)
      ENDDO

      XM(1)=E1(2)*E2(3)-E1(3)*E2(2)
      XM(2)=E1(3)*E2(1)-E1(1)*E2(3)
      XM(3)=E1(1)*E2(2)-E1(2)*E2(1)
      RNX=0.D0
      DO J=1,3
      RNX=RNX+XN(J)*(X0(J)-Y(J))
      ENDDO
      RMX=RNX

      P=XN(1)*XM(1)+ XN(2)*XM(2)+ XN(3)*XM(3)
      IF (P.LT.0.D0) THEN
	RMX= -RNX
	DO I=1,3
        DO J=1,3
	STR(I,J)=-STR(I,J)
	ENDDO
	ENDDO
      END IF

      CALL SPHAREA(Y,X,S1)
      DO J=1,3
      X(J,2)=X(J,3)
      X(J,3)=X(J,4)
      ENDDO
      CALL SPHAREA(Y,X,S2)
      OM= S1+S2
      S=(S1+S2)/3.D0
      IF(RNX.LT.0.D0) S= -S
      DO I=1,3
      STR(I,I)=STR(I,I)+S
      ENDDO

      IF (.NOT.EXT) RETURN
C     ADDL INTEGRALS OVER RECTANGLE FOR HIGH-ORDER NEAR-SINGULARITY SDUBTRACTION

      DO J=1,3
      R(J,1)=X0(J)-Y(J)
      R(J,2)= R(J,1)+L1*E1(J)
      R(J,3)= R(J,2)+L2*E2(J)
      R(J,4)= R(J,3)-L1*E1(J)
      R(J,5)=R(J,1)
      ENDDO

      DO J=1,3
      RLT(J)=0.D0
      ENDDO
      DO K=1,4
      D(K)=0.D0
      DO J=1,3
      D(K)=D(K)+ R(J,K)**2
      ENDDO
      D(K)=DSQRT(D(K))
      ENDDO
      D(5)=D(1)

      DO 1 K=1,4
      XA=0.D0
      XB=0.D0
      DO J=1,3
      XA=XA+ R(J,K)*E(J,K)
      XB=XB+ R(J,K+1)*E(J,K)
      ENDDO
      RL=DLOG ((XB+D(K+1))/(XA+D(K)))
      HH=0.D0
      DO J=1,3
      H(J)=R(J,K)-XA*E(J,K)
      HH= HH+H(J)**2
      ENDDO
      R3L=(XB/D(K+1)-XA/D(K))/HH
      XR3L= 1.D0/D(K)-1.D0/D(K+1)
      XXR3L= RL-HH*R3L

      DO I=1,3
      RLT(I)=RLT(I)+RL*E(I,K)
      DO J=I,3
      RR_R3L(I,J,K)=H(I)*H(J)*R3L+ (H(I)*E(J,K)+H(J)*E(I,K))*XR3L+
     >XXR3L*E(I,K)*E(J,K)
      RR_R3L(J,I,K)= RR_R3L(I,J,K)
      ENDDO
      ENDDO
1     CONTINUE

      R_R3S(1)= XM(2)*RLT(3)-XM(3)*RLT(2)
      R_R3S(2)= XM(3)*RLT(1)-XM(1)*RLT(3)
      R_R3S(3)= XM(1)*RLT(2)-XM(2)*RLT(1)

      DO J=1,3
      IF (RMX.GT.0.D0) THEN
	R_R3S(J)= R_R3S(J) + OM*XM(J)
      ELSE
	R_R3S(J)= R_R3S(J) - OM*XM(J)
      END IF
      ENDDO

      DO I=1,3
      DO J=I,3
      STR1(I,J)= (-RR_R3L(I,J,2) + RR_R3L(I,J,4) +
     >E1(I)*R_R3S(J)+ E1(J)*R_R3S(I) )*RNX/3.D0

      STR2(I,J)= ( RR_R3L(I,J,1) - RR_R3L(I,J,3) +
     >E2(I)*R_R3S(J)+ E2(J)*R_R3S(I) )*RNX/3.D0

      STR1(J,I)= STR1(I,J)
      STR2(J,I)= STR2(I,J)
      ENDDO
      ENDDO

      RETURN
      END


      SUBROUTINE CONTGEN(XP,A,ANGLE,x0,E1,L1,NPC,X0C,SOL)

      IMPLICIT REAL *8 (A-H,O-Z)          
      REAL *8 L1
      LOGICAL SOL         
      INCLUDE 'NPCH.h'

      DIMENSION XPROT(2), X0ROT(2,NPCH), E1ROT(2,NPCH), X0C(2,*),
     >XP(2),E1(2,*),X0(2,NPCH), L1(*),SOL(*)

      PI=DACOS(-1.D0)

      CO=COS(ANGLE*PI/180.D0)
      SI=SIN(ANGLE*PI/180.D0) 
          
      XPROT(1)= XP(1)*CO+XP(2)*SI
      XPROT(2)=-XP(1)*SI+XP(2)*CO

      DO IPC=1,NPCH
      X0ROT(1,IPC)=X0(1,IPC)*CO+X0(2,IPC)*SI
      X0ROT(2,IPC)=-X0(1,IPC)*SI+X0(2,IPC)*CO

      E1ROT(1,IPC)=E1(1,IPC)*CO+E1(2,IPC)*SI
      E1ROT(2,IPC)=-E1(1,IPC)*SI+E1(2,IPC)*CO
      ENDDO

      CALL CONTOR (XPROT,A,X0ROT,E1ROT,L1,NPC,X0C,SOL)

      DO IPC=1,NPC
      R1=X0C(1,IPC)
      R2=X0C(2,IPC)
      X0C(1,IPC)=R1*CO-R2*SI
      X0C(2,IPC)=R1*SI+R2*CO
      ENDDO
      RETURN
      END

