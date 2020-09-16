C     FOR ARBITARY NNEI and POLY
      SUBROUTINE RECRESTR(CHECK,IP,QLIM,alpha,beta,GAMMA,coeff,CMESH,
     >ITRESTR,V,XNA,QW,QWP,INB,INF,NTB,NTF,NEB,NEF,IM,JM,
     >NT1IJ,K1IJ,NT2IJ,K2IJ,IN,nnei,IV,QMINIP,DXMINIP,DXMAXIP,
     >ITRESTIP, RESTR)
      IMPLICIT REAL *8 (A-H,O-Z)
      INCLUDE 'NTMAX.h'
      include 'NEIBMAX.h'
      PARAMETER (NTD=NTMAX,nvd=NTD/2+2, NED=3*NVD, MCMAX=100)

      dimension V(3,*),XNA(3,*),INB(*),INF(*),NTB(*),
     >NTF(*),IN(neibmax,*),GRAD(3,nvd), QW(3,*),QWP(3,*),
     >IV(3,*),QMINIP(*),DXMINIP(*),DXMAXIP(*),ITRESTIP(*),
     >T(5,13),V0(3,nvd),E(3,3,nvd), nnei(*),
     >QW0(3,nvd),X(3),Y(5),FIT(5,nvd),FITQW(5,3,nvd),
     >JJMIN(nvd),RM(3),QWP0(3,nvd),FITQWP(5,3,nvd),
     >XNA0(3,NVD),CC(2,3),RIMIN(NVD),FT(3),RLA(NVD),RLA0(NVD),
     >DS(NTD), FITLA(5,nvd),PHI(NVD)

      DIMENSION  NEB(*),NEF(*),IM(*),JM(*), K1IJ(*),K2IJ(*),
     >NT1IJ(*),NT2IJ(*),QBETA(NTD),XC(3,NTD),XN(3),
     >NNE(NTD),NET(3,NTD), IM0(NED),JM0(NED),K1IJ0(NED),K2IJ0(NED),
     >NT1IJ0(NED), NT2IJ0(NED),NNEI0(NVD),IN0(NEIBMAX,NVD),IV0(3,NTD)

      DIMENSION XMIN(3),XMAX(3),MC(3),KC(3),HC(3),LL(NVD),DM(3)
      INTEGER HOC(0:MCMAX+1,0:MCMAX+1,0:MCMAX+1)
      REAL *8 LA
      INTEGER ALPHA,BETA
      LOGICAL CHECK, RESTR

c      print *, ' INB=', INB(IP), ' INF=', INF(IP)
c      print *, ' NTB=', NTB(IP), ' NTF=', NTF(IP)
c      print *, ' NEB=', NEB(IP), ' NEF=', NEF(ip) 
      RESTR=.true.
      STOT=0
      NT0=NTB(IP)-1
      DO JT=NTB(IP), NTF(IP)
      JTT=JT-NT0
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
      DS(JTT)=.5D0*DSQRT(CX1**2+CY1**2+CZ1**2)
      STOT=STOT+DS(JTT)
      ENDDO
c!    crude estimation:
       pi=dacos(-1.d0)
       rad= dsqrt(stot/(4.d0*pi))

      A2=4*STOT/(NTF(IP)-NTB(IP)+1)
      A2= A2/DSQRT(3.D0)
      AEST= DSQRT(A2)
      I0=INB(IP)-1

      DO I=INB(IP),INF(IP)
      PHI(I-I0)=0.D0
      ENDDO
      DO K=NTB(IP),NTF(IP)
      DO KJ=1,3
      I= IV(KJ,K)
      II=I-I0
      PHI(II)=PHI(II)+DS(K-NT0)/3.D0
      ENDDO
      ENDDO

C      JGG0=IGGB(IP)-1
      IS2=1
      JS2=17
C      DO JGG=IGGB(IP),IGGF(IP)
C      JJGG=JGG-JGG0
C      CALL SHELL (IS2,JS2,INBGG(JGG),INFGG(JGG),IORDGG,
C     >V,XCGG(1,JJGG),RADGG(JJGG),IMAX)
C      ENDDO


      DO I=INB(IP),INF(IP)
      II=I-I0
      DO K=1,3
      QW0(K,II)=QW(K,I)
      ENDDO
      ENDDO

c     with automatic COEFF inside
      CALL GIM (INB(IP),INF(IP),V,IN,nnei,QW,XNA)
      CURMAX=0.d0
      DO  I=INB(IP),INF(IP)
      II=I-I0
      DO K=1,3
      V0(K,II)=V(K,I)
      FT(K)=QW(K,I)
      QW(K,I)=QW0(K,II)
      QWP0(K,II)=QWP(K,I)
      XNA0(K,II)=XNA(K,I)
      E(3,K,II)=XNA(K,I)
      ENDDO
c-------------------------------
      CUR=-(FT(1)+FT(3))
      RLA(ii)=DABS(CUR)+DSQRT ((FT(1)-FT(3))**2 + FT(2)**2)
      RLA(II)=RLA(II)*RAD
      CURMAX=MAX(CURMAX,RLA(II))

      RLA(II)=  FT(1)**2+FT(3)**2+FT(2)**2/2.D0+0.001d0/RAD**2
      RLA(II)=RLA(II)**(-GAMMA/2.d0)
      RLA0(II)=RLA(II)
      ENDDO

      CONST=0.D0
      DO 1 I=INB(IP),INF(IP)
      II=I-I0
      if (PHI(II).lt. 0.d0) then
        print *, ' phi<0'
        stop
      end if
      CONST=CONST+ PHI(II)/RLA(II)**2

      II=I-I0
      A=XNA(1,I)
      B=XNA(2,I)
      C=XNA(3,I)
      IF(A*A .LE. B*B .AND. A*A .LE. C*C) THEN
      E(1,1,II)=0.D0
      E(1,2,II)=C
      E(1,3,II)=-B
      ELSE
	IF(B*B .LE. A*A .AND. B*B .LE. C*C )THEN
	 E(1,1,II)=C
	 E(1,2,II)=0.D0
	 E(1,3,II)=-A
	ELSE
	 E(1,1,II)=B
	 E(1,2,II)=-A
	 E(1,3,II)=0.D0
	END IF
      END IF
      E1L=1.D0/DSQRT(E(1,1,II)**2+E(1,2,II)**2+E(1,3,II)**2)
      DO J=1,3
      E(1,J,II)=E(1,J,II)*E1L
      ENDDO
      E(2,1,II)=E(1,2,II)*C - E(1,3,II)*B
      E(2,2,II)=E(1,3,II)*A - E(1,1,II)*C
      E(2,3,II)=E(1,1,II)*B - E(1,2,II)*A
      E2L=1.D0/DSQRT(E(2,1,II)**2+E(2,2,II)**2+E(2,3,II)**2)
      DO J=1,3
      E(2,J,II)=E(2,J,II)*E2L
      ENDDO

      DO K=1,5
      DO M=1,13
      T(K,M)=0.D0
      ENDDO
      ENDDO
      DO 2 LLL=1,nnei(i)
      J=IN(LLL,I)
      JJ=J-I0
      DO K=1,3
      X(K)=0.D0
      DO M=1,3
      X(K)=X(K)+E(K,M,II)*(V(M,J)-V(M,I))
      ENDDO
      ENDDO
      Y(1)=X(1)
      Y(2)=X(2)
      Y(3)=X(1)**2
      Y(4)=X(1)*X(2)
      Y(5)=X(2)**2

      DO K=1,5
      DO M=K,5
      T(K,M)=T(K,M)+Y(K)*Y(M)
      ENDDO
      T(K,6)=T(K,6)+Y(K)*X(3)
      DO L=1,3
      T(K,6+L)=T(K,6+L)+Y(K)*(QW(L,J)-QW(L,I))
      T(K,9+L)=T(K,9+L)+Y(K)*(QWP(L,J)-QWP(L,I))
      ENDDO
      T(K,13)=T(K,13)+Y(K)*(RLA(JJ)-RLA(II))
      ENDDO
2     CONTINUE

      DO K=1,5
      DO M=1,K-1
      T(K,M)=T(M,K)
      ENDDO
      ENDDO
      CALL GM5_13(5,13,T)
      DO K=1,5
      FIT(K,II)=T(K,6)
      DO L=1,3
      FITQW(K,L,II)=T(K,6+L)
      FITQWP(K,L,II)=T(K,9+L)
      ENDDO
      FITLA(K,II)=T(K,13)
      ENDDO
1     CONTINUE
      CONST=4.D0*CONST/(NTF(IP)-NTB(IP)+1)
      CONST=CONST/DSQRT(3.D0)

      print *, ' curmax=', curmax
ccc      if (check .and. curmax .gt.15.d0) return

c       connections backup:
	DO I=INB(IP),INF(IP)
	NNEI0(I-I0)=NNEI(I)
	DO K=1,NNEI(I)
	IN0(K,I-I0)=IN(K,I)
	ENDDO

	ENDDO

	NE0=NEB(IP)-1
	DO NE=NEB(IP),NEF(IP)
	NEE=NE-NE0
	IM0(NEE)=IM(NE)
	JM0(NEE)=JM(NE)
	K1IJ0(NEE)=K1IJ(NE)
	K2IJ0(NEE)=K2IJ(NE)
	NT1IJ0(NEE)=NT1IJ(NE)
	NT2IJ0(NEE)=NT2IJ(NE)
	ENDDO

	JT0=NTB(IP)-1
	DO JT=NTB(IP),NTF(IP)
	DO K=1,3
	IV0(K,JT-JT0)=IV(K,JT)
	ENDDO
	ENDDO

	RAT2AMIN=1.d40
        QMINMAX=0.d0
        FMIN=1.d200
	IT=0

55	IT=IT+1

	DO  I=INB(IP),INF(IP)
	II=I-I0
        if (RLA(II).lt.0.d0) then
         print *, ' RLA(II)<0'
         open (17, FILE='FEND')
         stop
        end if
	RLA(II)= 0.5D0*CONST*RLA(II)**2
	ENDDO

      if (it.gt. 150) goto 4001
c       RECONNECTION part:

      IREC=0
      I0=INB(IP)-1

      QM=1.D20
4000  DO 770 JT=NTB(IP), NTF(IP)
      DO M=1,3
      M2=M+1
      IF (M2.GT.3) M2=M2-3
      I1=IV(M,JT)
      I2=IV(M2,JT)
      RM(M)=0.D0
      DO K=1,3
      RM(M)=RM(M)+(V(K,I2)-V(K,I1))**2
      ENDDO
      ENDDO
      S=4.D0*RM(1)*RM(2)-(RM(1)+RM(2)-RM(3))**2

      SUM=RM(1)+RM(2)+RM(3)
      S= 0.25D0*DSQRT (S)
      FS=FS+S
      QM= MIN(QM,S/SUM)
      QBETA(JT-JT0) =(SUM/S) **beta

      I1=IV(1,JT)
      I2=IV(2,JT)
      I3=IV(3,JT)
      call TRNG(V(1,I1),V(1,I2),V(1,I3),XC(1,JT-JT0),RR,XN,S)
770   CONTINUE
c      print *, ' 770 complete QM=', QM, ' IT=', IT, 
c     >' QLIM=', QLIM

      IF(IT.EQ.1 .AND.QM.GT. QLIM) THEN
	RESTR=.FALSE.
	RETURN
      END IF
      DO  JT=NTB(IP), NTF(IP)
      NNE(JT-JT0)=0
      ENDDO

      DO NE=NEB(IP),NEF(IP)
      NT1=NT1IJ(NE)-JT0
      NT2=NT2IJ(NE)-JT0
      NNE(NT1)=NNE(NT1)+1
      NET(NNE(NT1),NT1)=NE
      NNE(NT2)=NNE(NT2)+1
      NET(NNE(NT2),NT2)= -NE
      ENDDO

      DO 1000 NE=NEB(IP),NEF(IP)
      I=IM(NE)
      J=JM(NE)
      K1=K1IJ(NE)
      K2=K2IJ(NE)
      NT1=NT1IJ(NE)
      NT2=NT2IJ(NE)

      i1=iv(1,nt1)+iv(2,nt1)+iv(3,nt1)
      i2=iv(1,nt2)+iv(2,nt2)+iv(3,nt2)
      if (i1.ne.i+j+k1 .or. i2.ne.i+j+k2) then
       print *, ' error in connectivities IP=', IP
       print *, ' i=', i, ' j=', j
       print *, ' k1=', k1, ' k2=', k2
       print *, iv(1,nt1),iv(2,nt1),iv(3,nt1)
       print *, iv(1,nt2),iv(2,nt2),iv(3,nt2)
       open (17, file='FEND')
       stop
      end if

      RM(1)=0.D0
      RM(2)=0.D0
      RM(3)=0.D0
      RR=0.D0
      P1=0.D0
      P2=0.D0
      DO K=1,3
      RM(1)=RM(1)+(V(K,I)-V(K,K1))**2
      RM(2)=RM(2)+(V(K,K2)-V(K,I))**2
      RM(3)=RM(3)+(V(K,K2)-V(K,K1))**2
      RR=RR+ (V(K,J)-V(K,I))**2
      P1=P1+(V(K,K1)-XC(K,NT1-JT0))* (V(K,K2)-V(K,K1))
      P2=P2+(V(K,K2)-XC(K,NT2-JT0))* (V(K,K2)-V(K,K1))
      ENDDO
      S=4.D0*RM(1)*RM(2)-(RM(1)+RM(2)-RM(3))**2
      QI =0.25D0*DSQRT (S)/(RM(1)+RM(2)+RM(3))
      RM(1)=0.D0
      RM(2)=0.D0
      DO K=1,3
      RM(1)=RM(1)+(V(K,J)-V(K,K2))**2
      RM(2)=RM(2)+(V(K,K1)-V(K,J))**2
      ENDDO
      S=4.D0*RM(1)*RM(2)-(RM(1)+RM(2)-RM(3))**2
      QJ =0.25D0*DSQRT (S)/(RM(1)+RM(2)+RM(3))
      A2=RLA(I-I0)+RLA(J-I0)
      FOLD= ((RR/A2+A2/RR)/2.D0)**ALPHA +COEFF*(QBETA(NT1-JT0) +
     >QBETA(NT2-JT0))
      A2=RLA(K1-I0)+RLA(K2-I0)
      FNEW= ((RM(3)/A2+A2/RM(3))/2.D0)**ALPHA +COEFF/QI**BETA +
     >COEFF/QJ**BETA
      IF (RM(3).GT.DABS(P1)+DABS(P2).OR. FNEW.GT.FOLD) GOTO 1000

      IF (NNEI(I).LT.6 .OR. NNEI(J).LT.6 ) THEN
C	PRINT *, ' TOO FEW NEIGHBORS'
c	OPEN (17, FILE='FEND')
c	STOP
        goto 1000
      END IF

      IF (NNEI(K1).EQ. NEIBMAX .OR. NNEI(K2).EQ.NEIBMAX) THEN
	PRINT *, ' NEIBMAX TOO SMALL'
	OPEN (17, FILE='FEND')
	STOP
      END IF

      IF (K2.LT.K1) THEN
	 IM(NE)=K2
	 JM(NE)=K1
	 K1IJ(NE)=I
	 K2IJ(NE)=J
	 IV(1,NT1)=I
	 IV(2,NT1)=K2
	 IV(3,NT1)=K1

	 IV(1,NT2)=J
	 IV(2,NT2)=K1
	 IV(3,NT2)=K2
	 NTI=NT1
	 NTJ=NT2
      ELSE
	 IM(NE)=K1
	 JM(NE)=K2
	 K1IJ(NE)=J
	 K2IJ(NE)=I
	 IV(1,NT1)=J
	 IV(2,NT1)=K1
	 IV(3,NT1)=K2

	 IV(1,NT2)=I
	 IV(2,NT2)=K2
	 IV(3,NT2)=K1
	 NTI=NT2
	 NTJ=NT1
      END IF

	 DO 33 K=1,3
	 NE1= NET(K,NT1-JT0)
	 IF (NE1.EQ.NE) GOTO 34
	 IF (NE1.GT.0 .AND. K1IJ(NE1).EQ.I) THEN
	   NT1IJ(NE1)=NTJ
	   K1IJ(NE1)=K2
	 END IF
	 IF (NE1.LT.0 .AND. K2IJ(-NE1).EQ.I) THEN
	   NT2IJ(-NE1)=NTJ
	   K2IJ(-NE1)=K2
	 END IF
	 IF (NE1.GT.0 .AND. K1IJ(NE1).EQ.J) THEN
	   NT1IJ(NE1)=NTI
	   K1IJ(NE1)=K2
	 END IF
	 IF (NE1.LT.0 .AND. K2IJ(-NE1).EQ.J) THEN
	   NT2IJ(-NE1)=NTI
	   K2IJ(-NE1)=K2
	 END IF

34	 NE1= NET(K,NT2-JT0)
	 IF(NE1.EQ. -NE) GOTO 33
	 IF (NE1.GT.0 .AND. K1IJ(NE1).EQ.I) THEN
	   NT1IJ(NE1)=NTJ
	   K1IJ(NE1)=K1
	 END IF
	 IF (NE1.LT.0 .AND. K2IJ(-NE1).EQ.I) THEN
	   NT2IJ(-NE1)=NTJ
	   K2IJ(-NE1)=K1
	 END IF
	 IF (NE1.GT.0 .AND. K1IJ(NE1).EQ.J) THEN
	   NT1IJ(NE1)=NTI
	   K1IJ(NE1)=K1
	 END IF
	 IF (NE1.LT.0 .AND. K2IJ(-NE1).EQ.J) THEN
	   NT2IJ(-NE1)=NTI
	   K2IJ(-NE1)=K1
	 END IF
33    CONTINUE

      NNEI(K1)=NNEI(K1)+1
      IN(NNEI(K1),K1)=K2
      NNEI(K2)=NNEI(K2)+1
      IN(NNEI(K2),K2)=K1

      DO K=1,NNEI(I)
      IF (IN(K,I).EQ.J) THEN
	DO L=K,NNEI(I)-1
	IN(L,I)=IN(L+1,I)
	ENDDO
	GOTO 2000
      END IF
      ENDDO
2000  NNEI(I)=NNEI(I)-1

      DO K=1,NNEI(J)
      IF (IN(K,J).EQ.I) THEN
	DO L=K,NNEI(J)-1
	IN(L,J)=IN(L+1,J)
	ENDDO
	GOTO 3000
      END IF
      ENDDO
3000  NNEI(J)=NNEI(J)-1
      IREC=IREC+1
      GOTO 4000
1000  CONTINUE
      print *, ' IREC=', IREC

4001	qmin=1.d20
        qa=0
	RMIN=1.D20
	RMAX=0.D0
        R2AMIN=1.d20
        R2AMAX=0.d0
	la=1.d40

	F=0.D0
	DO I=INB(IP),INF(IP)
	II=I-I0
	RIMIN(II)=1.D40
	DO L=1,3
	GRAD(L,II)=0.D0
	ENDDO
	ENDDO


	a2min=1.d40
        a2max=0.d0
	DO 17 I=INB(IP),INF(IP)
	II=I-I0
	DO 18 K=1,nnei(i)
	J=IN(K,I)
        IF (J.LE.I) GOTO 18
	JJ=J-I0
	A2=RLA(II)+RLA(JJ)
        A2INV=1.d0/A2
	RR=0.D0
        DO L=1,3
        RR=RR+(V(L,J)-V(L,I))**2
        ENDDO

        temp=(RR/A2+A2/RR)/2
        R1=  temp**ALPHA
        F=F+ R1
        R2=ALPHA*R1*(A2INV-A2/RR**2)/temp
        DO L=1,3
        GRAD(L,II)=GRAD(L,II) -R2*(V(L,J)-V(L,I))
        GRAD(L,J-I0)=GRAD(L,J-I0) +R2*(V(L,J)-V(L,I))
        ENDDO
        if (a2.lt.a2min) then
         a2min=a2
         rr_atmin=rr
        end if
        if (a2.gt.a2max) then
          a2max=a2
          rr_atmax=rr
        endif
        R2A=RR/A2

	RMIN=MIN(RMIN,RR)
	RMAX=MAX(RMAX,RR)
	R2AMIN=MIN(R2AMIN,R2A)
	R2AMAX=MAX(R2AMAX,R2A)

	RIMIN(II)=MIN(RIMIN(II),RR)
	RIMIN(J-I0)=MIN(RIMIN(J-I0),RR)

18	CONTINUE
17      CONTINUE
      amin=dsqrt(a2min)
      amax=dsqrt(a2max)
      rr_atmin = dsqrt (rr_atmin)
      rr_atmax = dsqrt (rr_atmax)
      R2AMIN= dsqrt(R2AMIN)
      R2AMAX= dsqrt(R2AMAX)

      rat=dsqrt(rmax/rmin)
      rat2a= r2amax/r2amin
      rat2amin=min(rat2amin,rat2a)
      if (it.eq.1) rat2a0=rat2a

      JT0=NTB(IP)-1
      DO 77 JT=NTB(IP), NTF(IP)
      DO M=1,3
      M2=M+1
      IF (M2.GT.3) M2=M2-3
      I1=IV(M,JT)
      I2=IV(M2,JT)
      RM(M)=0.D0
      DO K=1,3
      RM(M)=RM(M)+(V(K,I2)-V(K,I1))**2
      ENDDO
      ENDDO
      S=4.D0*RM(1)*RM(2)-(RM(1)+RM(2)-RM(3))**2

      SUM=RM(1)+RM(2)+RM(3)
      Q=0.25D0*DSQRT (S)/SUM
      qmin=min(qmin,q)
      qa=qa+q/(ntf(ip)-ntb(ip)+1)
      R1=COEFF/Q**BETA
      F=F+ R1
      R2= -BETA*R1/Q**2
      DO M=1,3
      D=0.125D0*(SUM-2.D0*RM(M))
      D=R2*(D/SUM**2 -2.D0*Q**2/SUM)
      M2=M+1
      IF (M2.GT.3) M2=M2-3
      I1=IV(M,JT)
      I2=IV(M2,JT)
      DO K=1,3
      GRAD(K,I1-I0)=GRAD(K,I1-I0) -D*(V(K,I2)-V(K,I1))
      GRAD(K,I2-I0)=GRAD(K,I2-I0) +D*(V(K,I2)-V(K,I1))
      ENDDO
      ENDDO
77    CONTINUE

      QMINMAX=MAX(QMINMAX,QMIN)
      IF(IT.EQ.1) QMIN0=QMIN

      if (mod(it,10).eq.1) then
         print *, ' amin=', amin, ' amax=',amax
         print *, ' rr_atmin=', rr_atmin, ' rr_atmax=', rr_atmax
         print *, ' r2amin, max=', r2amin,r2amax
	 PRINT *, ' IT=', IT, ' F=', F
         print *, ' qmin=', qmin, ' qa=', qa
         PRINT *, ' MESH RATIO=', DSQRT(RMAX/RMIN),
     >  ' rat2a=', rat2a
c       read (*,*) icont
       end if

      if (.not.check .or. it.ne.201) goto 170

      if (rat2amin.gt. rat2a0/1.15d0 .and. QMINMAX.LT.1.05d0*QMIN0)then
c      if (rat2amin.gt. rat2a0/1.05d0 .and. QMINMAX.LT.1.01d0*QMIN0)then
        do i=inb(ip),inf(ip)
        do k=1,3
        V(K,I)=V0(K,I-I0)
        enddo
        enddo

	DO I=INB(IP),INF(IP)
	NNEI(I)=NNEI0(I-I0)
	DO K=1,NNEI(I)
	IN(K,I)=IN0(K,I-I0)
	ENDDO
	ENDDO

	DO NE=NEB(IP),NEF(IP)
	NEE=NE-NE0
	IM(NE)=  IM0(NEE)
	JM(NE)=  JM0(NEE)
	K1IJ(NE)=  K1IJ0(NEE)
	K2IJ(NE)=  K2IJ0(NEE)
	NT1IJ(NE)= NT1IJ0(NEE)
	NT2IJ(NE)= NT2IJ0(NEE)
	ENDDO

	DO JT=NTB(IP),NTF(IP)
	DO K=1,3
	IV(K,JT)= IV0(K,JT-JT0)
	ENDDO
	ENDDO

	print *, ' return'
c        read (*,*) icont
        RESTR=.false.
        return
      else
        print *, ' forceful RESTRUCT IP=', IP
c        write (4,*) ' forceful RESTRUCT IP=', IP
CC        read (*,*)icont
      end if

170     DO 14 I= INB(IP),INF(IP)
	II=I-I0

	WW=0.D0
	DO K=1,3
	WW=WW+GRAD(K,II)**2
	ENDDO
	temp=cmesh*dsqrt(rimin(II)/ww)
	la=min(la,temp)
 14	CONTINUE

	DO 140 I=INB(IP),INF(IP)
	II=I-I0
	DO L=1,3
	V(L,I)= -la*GRAD(L,II)+V(L,I)
	ENDDO
 140	CONTINUE

      DO K=1,3
      XMIN(K)=  1.D20
      XMAX(K)= -1.D20
      DO I=INB(IP),INF(IP)
      II=I-I0
      XMIN(K)= MIN(XMIN(K), V(K,I), V0(K,II) )
      XMAX(K)= MAX(XMAX(K), V(K,I), V0(K,II) )
      ENDDO
      CALL GR2(IS2,JS2,DR)
      XMIN(K)=XMIN(K) -1.D-4*DR
      XMAX(K)=XMAX(K)+ 1.D-4*DR
      HC(K)=2.D0*AEST
      MC(K)= (XMAX(K)-XMIN(K))/HC(K)
      HC(K)= (XMAX(K)-XMIN(K))/MC(K)
      IF (MC(K).GT.MCMAX) THEN
	OPEN (17, FILE='FEND')
	WRITE (17,*) ' MCMAX TOO SMALL IN RECRESTR'
	STOP
      END IF
      ENDDO

      DO K1=0, MC(1)+1
      DO K2=0, MC(2)+1
      DO K3=0, MC(3)+1
      HOC(K1,K2,K3)=0
      ENDDO
      ENDDO
      ENDDO

      DO I=INB(IP),INF(IP)
      II=I-I0
      DO K=1,3
      KC(K)= (V0(K,II)-XMIN(K))/HC(K)
      KC(K)=KC(K)+1
      ENDDO
      LL(II)=HOC(KC(1),KC(2),KC(3))
      HOC(KC(1),KC(2),KC(3))=II
      ENDDO

	distmax=0
	DO 3 I=INB(IP), INF(IP)
	II=I-I0

C	CL=1.D40
C        DO  JGG=IGGB(IP),IGGF(IP)
C	JJGG=JGG-JGG0
C	P=0.D0
C	DO K=1,3
C	P=P+(XCGG(K,JJGG)-V(K,I))**2
C	ENDDO
C	CL= MIN(CL, DSQRT(P)+RADGG(JJGG))
C	ENDDO

C        DO  100 JGG=IGGB(IP),IGGF(IP)
C	JJGG=JGG-JGG0
C	P=0.D0
C        DO K=1,3
C	P=P+(XCGG(K,JJGG)-V(K,I))**2
C	ENDDO
C	IF (P.GT.(RADGG(JJGG)+CL)**2) GOTO 100
C        DO JJJ= INBGG(JGG),INFGG(JGG)
C        J=IORDGG(JJJ)
C	JJ=J-I0
C	P=0.D0
C        DO K=1,3
C        P=P+(V0(K,JJ)-V(K,I))**2
C        ENDDO
C        IF (P.LT.CL**2) THEN
C          CL=DSQRT(P)
C	  JJMIN(II)=JJ
C	END IF
C	ENDDO
C100     CONTINUE

      CL=1.D40
      DO K=1,3
      TEMP=(V(K,I)-XMIN(K))/HC(K)
      KC(K)= TEMP
      TEMP=TEMP-KC(K)
      DM(K)= HC(K)*MIN(TEMP, 1.D0-TEMP) +HC(K)
      KC(K)=KC(K)+1
      if (kc(k).le.0 .or. KC(K).gt. MC(k))  then
	print *, '  error  in recrest'
	print *,  kc(k), mc(k)
	open (17, file='FEND')
	stop
      end if
      ENDDO
      DO 100 K1=KC(1)-1, KC(1)+1
      DO 100 K2=KC(2)-1, KC(2)+1
      DO 100 K3=KC(3)-1, KC(3)+1
      JJ=HOC (K1,K2,K3)
101   IF (JJ.EQ.0) GOTO 100
      P=0.D0
      DO K=1,3
      P=P+(V0(K,JJ)-V(K,I))**2
      ENDDO
      IF (P.LT.CL) THEN
        CL=P
        JJMIN(II)=JJ
      END IF
      JJ= LL(JJ)
      GOTO 101
100   CONTINUE
      IF (DSQRT(CL).GE. MIN(DM(1),DM(2),DM(3))) THEN
	OPEN (17, FILE='FEND')
	WRITE (17,*) ' HC TOO SMALL IN RECRESTR'
	STOP
      END IF

      JJ=JJMIN(II)
      dist=0
      do m=1,3
      dist=dist+(v(m,i)-v0(m,jj))**2
      enddo
      distmax=max(distmax,dist)
      DO K=1,2
      X(K)=0.D0
      DO M=1,3
      X(K)=X(K)+E(K,M,JJ)*(V(M,I)-V0(M,JJ))
      ENDDO
      ENDDO
      Y(1)=X(1)
      Y(2)=X(2)
      Y(3)=X(1)**2
      Y(4)=X(1)*X(2)
      Y(5)=X(2)**2
      X(3)=0.D0
      RLA(II)=RLA0(JJ)
      DO K=1,5
      X(3)= X(3)+FIT(K,JJ)*Y(K)
      RLA(II)= RLA(II)+FITLA(K,JJ)*Y(K)
      ENDDO

      DO K=1,3
      V(K,I)=V0(K,JJ)
      DO L=1,3
      V(K,I)=V(K,I)+X(L)*E(L,K,JJ)
      ENDDO
      ENDDO
3     continue

cc      print *,' distmax=', dsqrt(distmax)

c        if (it.eq.1) read (*,*) icont

	if (it.lt. itrestr) goto 55

56    QMINIP(IP)=1.D40
      DO JT=NTB(IP), NTF(IP)
      DO M=1,3
      M2=M+1
      IF (M2.GT.3) M2=M2-3
      I1=IV(M,JT)
      I2=IV(M2,JT)
      RM(M)=0.D0
      DO K=1,3
      RM(M)=RM(M)+(V(K,I2)-V(K,I1))**2
      ENDDO
      ENDDO
      S=4.D0*RM(1)*RM(2)-(RM(1)+RM(2)-RM(3))**2
      SUM=RM(1)+RM(2)+RM(3)
      S=0.25D0*DSQRT (S)/SUM
      IF (S.LT.QMINIP(IP)) QMINIP(IP)=S
      ENDDO

      DO I=INB(IP),INF(IP)
      II=I-I0
      JJ=JJMIN(II)
      DO L=1,3
      XNA(L,I)=XNA0(L,JJ)
      ENDDO
      ENDDO
ccc   extra      CALL GIM (INB(IP),INF(IP),V,IN,nnei,QW,XNA)

      DO I=INB(IP),INF(IP)
      II=I-I0
      JJ=JJMIN(II)
      DO K=1,2
      X(K)=0.D0
      DO M=1,3
      X(K)=X(K)+E(K,M,JJ)*(V(M,I)-V0(M,JJ))
      ENDDO
      ENDDO
      Y(1)=X(1)
      Y(2)=X(2)
      Y(3)=X(1)**2
      Y(4)=X(1)*X(2)
      Y(5)=X(2)**2
      DO L=1,3
      QW(L,I)=QW0(L,JJ)
      QWP(L,I)=QWP0(L,JJ)
      DO K=1,5
      QW(L,I)=QW(L,I)+FITQW(K,L,JJ)*Y(K)
      QWP(L,I)=QWP(L,I)+FITQWP(K,L,JJ)*Y(K)
      ENDDO
      ENDDO
      ENDDO

      DXMINIP(IP)=dsqrt(RMIN)
      DXMAXIP(IP)=dsqrt(RMAX)
      ITRESTIP(IP)=IT
      PRINT *, ' QMIN=', QMINIP(IP)
      return
      END

      SUBROUTINE GM5_13(N,M,T)
      real *8 T(5,13),R,U
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

      SUBROUTINE TRNG(X1,X2,X3,XC,RR,XN,S)
c     FOR A TRIANGLE GIVEN BY VERTICES (X1,X2,X3), CALCULATES CENTER XC OF
c     CIRCUMSCRIBED CIRCLE, SQUARE OF THE CIRCLE RADIUS (RR), UNIT NORMAL
c     VECTOR XN TO THE TRIANGLE PLANE (ARBITRARY ORIENTATION) AND TIRANGLE
c     AREA S

      IMPLICIT REAL *8 (A-H,O-Z)
      DIMENSION X1(3),X2(3),X3(3),XC(3),E2(3),E3(3),XN(3)
      G22=0.D0
      G23=0.D0
      G33=0.D0
      DO K=1,3
      E2(K)=X2(K)-X1(K)
      E3(K)=X3(K)-X1(K)
      G22=G22+E2(K)**2
      G23=G23+E2(K)*E3(K)
      G33=G33+E3(K)**2
      ENDDO
      DET=1.D0/(G22*G33-G23**2)
      RLA=0.5D0*G33*(G22-G23)*DET
      RMU=0.5D0*G22*(G33-G23)*DET
      RR=0.D0
      DO K=1,3
      XC(K)=X1(K)+RLA*E2(K)+RMU*E3(K)
      RR=RR+(XC(K)-X1(K))**2
      ENDDO
      XN(1)=E2(2)*E3(3)-E2(3)*E3(2)
      XN(2)=E2(3)*E3(1)-E2(1)*E3(3)
      XN(3)=E2(1)*E3(2)-E2(2)*E3(1)
      TEMP=0.D0
      DO K=1,3
      TEMP=TEMP+XN(K)**2
      ENDDO
      TEMP=DSQRT(TEMP)
      DO K=1,3
      XN(K)= XN(K)/TEMP
      ENDDO
      S=0.5D0*TEMP
      RETURN
      END
