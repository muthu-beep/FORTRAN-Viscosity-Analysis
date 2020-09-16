C      FOR ARBITARRY NNEI
c      OK for POLYTRNG
cc     same as for /rh, but with smoothing in a sedimentation-like manner
cc     simplified version of qminsh (curvature adaptation term excluded;
cc    oly ratio of coeff1/coeff2 matters;
      SUBROUTINE QMINSIAD(COEFF1,COEFF2,GAMMA,IN,nnei,IV,V,
     >INB,INF,NTB,NTF,XNA,PHI,QM,FT2,VEL,RATMAX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'NEIBMAX.h'
      INCLUDE 'NTMAX.h'
      PARAMETER (NP=1, NT=NTMAX)
      PARAMETER (pi=3.14159265358979d0,NV=NT/2+2)
      DIMENSION XNA(3,*),IN(NEIBMAX,*),V(3,*),VEL(3,*),QM(*),FT2(*),
     >F(3,NV),PHI(*),DS(NT),DDS(3,3,NT),IV(3,*), RM(3),
     >G(3,NV),IN3(neibmax,NV),IV3(3,NT),DX(3),FIT(3),FICT(3),
     >XX(neibmax,3),
     >INB(*),INF(*),NTB(*),NTF(*),WEIM(3*NV),SM(3),
     >YNA(3,NV),DER(3,3,neibmax,NV),CUR(NV),E(3,3,NV),
     >FLUX(NV),nnei(*),RLA(NV)
      COMMON/SMOOTHL/SMOOTHL/ipmax/ipmax /qmin/qmin
      COMMON/QMINIP/QMINIP (NP)/RAT/RAT(NP) /RAD/RAD(NP)
      common /N3/ N3(NP)/N4/N4(NP)
      common /N5/ N5(NP)/N7/N7(NP)/N8/N8(NP)
      LOGICAL FAIL
      SMOOTHL= 1.d-5
      SMOOTHL= 0
      RATMAX=0.D0
      fluxmax=-1.d20
      fluxmin=1.d20
      fluxa=0.d0
      stot=0.d0
      itvtmax=0
      itvta=0
      qmin=1.d20
      DO 111 IP=1, NP
      QMINIP(IP)=1.D20
      RAT(IP)=0.D0
      I0=INB(IP)-1
      RMAX=0.D0
      RMIN=1.D20
      FS=0.D0
      n5(ip)=0
      n7(ip)=0
      n8(ip)=0
      DO 30 I=INB(IP),INF(IP)
C      if (nnei(i).eq.3) n3(ip)=n3(ip)+1
C      if (nnei(i).eq.4) n4(ip)=n4(ip)+1
C      if (nnei(i).eq.5) n5(ip)=n5(ip)+1
C      if (nnei(i).eq.7) n7(ip)=n7(ip)+1
C      if (nnei(i).eq.8) n8(ip)=n8(ip)+1
      II=I-I0
      DO L=1,nnei(i)
      J=IN(L,I)
      IN3(L,II)=J
	RR=0.D0
	DO K=1,3
	XX(L,K)=V(K,J)-V(K,I)
	RR=RR+(V(K,J)-V(K,I))**2
	ENDDO
	IF (RR.GT.RMAX) RMAX=RR
	IF (RR.LT.RMIN) RMIN=RR
      ENDDO
      DO K=1,3
      FICT(K)=XNA(K,I)
      ENDDO
cc      CALL MNADA(2,XX,FICT,E(1,1,II),FIT,YNA(1,II),DER(1,1,1,II),FAIL)
C        CALL MYNA(nnei(i),XX,XNA(1,I),E(1,1,II),FIT,YNA(1,II),
C     >	DER(1,1,1,II))
C        CUR(II)=-(FIT(1)+FIT(3))**3
      FS=FS+PHI(I)
30    CONTINUE

C      SUM=0.D0
C     DO I=INB(IP),INF(IP)
C      II=I-I0
C      FMAX=0.D0
C      DO L=1,nnei(i)
C      J=IN3 (L,II)
C	JJ=J-I0
C	RR=0.D0
C	FF=0.D0
C	DO K=1,3
C	RR=RR+(V(K,J)-V(K,I))**2
C        FF=FF+(YNA(K,JJ)-YNA(K,II))**2
C	ENDDO
C	FF=FF/RR
C        IF (FF.GT.FMAX) FMAX=FF
C      ENDDO
C      SUM=SUM+DSQRT(FMAX)*PHI(I)
C      ENDDO

      RAT(IP)=DSQRT (RMAX/RMIN)
      IF (RAT(IP).GT.RATMAX) RATMAX=RAT(IP)

C!      CALL MSMOOTH (INB(IP),INF(IP),IN3,nnei,V,PHI,E,CUR,FLUX)
cc!      CALL DEV (INB(IP),INF(IP),V,XNA,PHI,FS,QM,A,QA)
CCC!  do not use QA from DEV-bad idea behind
      qa=0.d0
      do i=inb(ip), inf(ip)
      qa=qa+dabs(qm(i))*phi(i)
      enddo
      qa=qa/fs
      a=rad(ip)
      FACTOR=SMOOTHL*a**5
      DO I=INB(IP), INF(IP)
      II=I-I0
        flux(ii)=flux(ii)*factor
        if (flux(ii).gt.fluxmax) fluxmax=flux(ii)
        if (flux(ii).lt.fluxmin) fluxmin=flux(ii)
C       QM(I)=QM(I)+QA*FLUX(II)
        fluxa=fluxa+dabs(flux(ii))*phi(i)
        stot=stot+phi(i)
      VN=VEL(1,I)*XNA(1,I)+VEL(2,I)*XNA(2,I)+VEL(3,I)*XNA(3,I)
      DO K=1,3
      VEL(K,I)=QM(I)*XNA(K,I)+VEL(K,I)-VN*XNA(K,I)
      G(K,II)=0.D0
      ENDDO
      ENDDO

c!!      FACTOR1=4.D0*COEFF1*(SUM/FS)**2
      FACTOR1=4.D0*COEFF1

c      print *, ' FACTOR1=', factor1
c!!      FACTOR2=COEFF2*(SUM/FS)**2
      FACTOR2=COEFF2

      RADEST=DSQRT(FS/(4.D0*PI))
      CONST=0.D0
      DO I=INB(IP),INF(IP)
      II=I-I0
      RLA(II)=  FT2(I) +0.001d0/RADEST**2
      RLA(II)=RLA(II)**(-GAMMA/2.d0)
      CONST=CONST+ PHI(I)/RLA(II)**2
      ENDDO

      CONST=4.D0*CONST/(NTF(IP)-NTB(IP)+1)
      CONST=CONST/DSQRT(3.D0)
      DO I=INB(IP),INF(IP)
      II=I-I0
      RLA(II)= 0.5D0*CONST*RLA(II)**2
      ENDDO

      NE=0
      DO 1 II=INB(IP), INF(IP)
      III=II-I0
      DO 2 J=1,nnei(ii)
      JJ=IN3(J,III)
      JJJ=JJ-I0
      IF (JJ.LT.II) GOTO 2
      NE=NE+1
      RR=0.D0
      DO K=1,3
      DX(K)=V(K,JJ)-V(K,II)
      RR=RR+DX(K)**2
      ENDDO
      A2=RLA(III)+RLA(JJJ)
      WEIM(NE)=COEFF1* (1.D0/A2-A2/RR**2)**2
2     CONTINUE
1     CONTINUE

      JT0=NTB(IP)-1
      DO JT=NTB(IP), NTF(IP)
      JTT=JT-JT0
      DO M=1,3
      IV3(M,JTT)=IV(M,JT)
      ENDDO
      DO M=1,3
      M2=M+1
      IF (M2.GT.3) M2=M2-3
      I1=IV3(M,JTT)
      I2=IV3(M2,JTT)
      RM(M)=0.D0
      DO K=1,3
      RM(M)=RM(M)+(V(K,I2)-V(K,I1))**2
      ENDDO
      ENDDO
      DS(JTT)=4.D0*RM(1)*RM(2)-(RM(1)+RM(2)-RM(3))**2
      SUM=RM(1)+RM(2)+RM(3)
      DS(JTT)=0.25D0*DSQRT (DS(JTT))/SUM
      IF (DS(JTT).LT.QMIN) QMIN=DS(JTT)
      IF (DS(JTT).LT.QMINIP(IP)) QMINIP(IP)=DS(JTT)
      DO M=1,3
      DO K=1,3
      DDS(K,M,JTT)=0.D0
      ENDDO
      ENDDO
      DO M=1,3
      D=0.125D0*(SUM-2.D0*RM(M))
      D=D/SUM**2 -2.D0*DS(JTT)**2/SUM
      M2=M+1
      IF (M2.GT.3) M2=M2-3
      I1=IV3(M,JTT)
      I2=IV3(M2,JTT)
      DO K=1,3
      DDS(K,M,JTT)=DDS(K,M,JTT)-D*(V(K,I2)-V(K,I1))
      DDS(K,M2,JTT)=DDS(K,M2,JTT)+D*(V(K,I2)-V(K,I1))
      ENDDO
      ENDDO
C     FOR SPEED:
      DS(JTT)=FACTOR2/DS(JTT)**4
      ENDDO

      FTLL=0.D0
      IT=0
200   IT=IT+1
      if (it.gt. itvtmax) itvtmax=it
      DO I=INB(IP),INF(IP)
      DO K=1,3
      F(K,I-I0)=0.D0
      ENDDO
      ENDDO

      NE=0
      FTL=0.D0
      DO 5 I=INB(IP), INF(IP)
      II=I-I0
      DO 6 L=1,nnei(i)
      J=IN3 (L,II)
      IF (J.LT.I) GOTO 6
      JJ=J-I0
      NE=NE+1
      F2=0.D0
      DO K=1,3
      DX(K)=V(K,J)-V(K,I)
      F2=F2+DX(K)*(VEL(K,J)-VEL(K,I))
      ENDDO
      TEMP=WEIM(NE)*F2
      FTL=FTL+ TEMP*F2
      F2=TEMP+TEMP
      DO K=1,3
      TEMP=F2*DX(K)
      F(K,II)=F(K,II)-TEMP
      F(K,JJ)=F(K,JJ)+TEMP
      ENDDO
6     CONTINUE
5     CONTINUE


      DO 80 JT=NTB(IP), NTF(IP)
      JTT=JT-JT0
      FF=0.D0
      DO M=1,3
      I=IV3(M,JTT)
      DO K=1,3
      FF=FF+DDS(K,M,JTT)*VEL(K,I)
      ENDDO
      ENDDO
      TEMP=DS(JTT)*FF
      FTL=FTL+TEMP*FF
      FF=TEMP+TEMP
      DO M=1,3
      I=IV3(M,JTT)
      II=I-I0
      DO K=1,3
      F(K,II)=F(K,II)+FF*DDS(K,M,JTT)
      ENDDO
      ENDDO
80    CONTINUE

      DO I=INB(IP),INF(IP)
      II=I-I0
      QQ=0.D0
      DO K=1,3
      QQ=QQ+F(K,II)*XNA(K,I)
      ENDDO
      DO K=1,3
      F(K,II)=F(K,II)-QQ*XNA(K,I)
      ENDDO
      ENDDO


      A11=0.D0
      A12=0.D0
      A22=0.D0
      B1=0.D0
      B2=0.D0
      NE=0
      DO 500 I=INB(IP), INF(IP)
      II=I-I0
      DO 600 L=1,nnei(i)
      J=IN3 (L,II)
      IF (J.LT.I) GOTO 600
      JJ=J-I0
      NE=NE+1
      F1F=0.D0
      F2F=0.D0
      F1G=0.D0
      F2G=0.D0
      DO K=1,3
      TEMP=V(K,J)-V(K,I)
      F2F=F2F+TEMP*(F(K,JJ)-F(K,II))
      F2G=F2G+TEMP*(G(K,JJ)-G(K,II))
      ENDDO
      TEMP=WEIM(NE)*F2F
      A11=A11+ TEMP*F2F
      A12=A12+ TEMP*F2G
      A22=A22+ WEIM(NE)*F2G**2
600   CONTINUE
500   CONTINUE


      DO 90 JT=NTB(IP), NTF(IP)
      JTT=JT-JT0
      FF=0.D0
      GG=0.D0
      DO M=1,3
      I=IV3(M,JTT)
      II=I-I0
      DO K=1,3
      FF=FF+DDS(K,M,JTT)*F(K,II)
      GG=GG+DDS(K,M,JTT)*G(K,II)
      ENDDO
      ENDDO
      TEMP=DS(JTT)*FF
      A11=A11+TEMP*FF
      A12=A12+TEMP*GG
      A22=A22+DS(JTT)*GG**2
90    CONTINUE

      B1=0.D0
      B2=0.D0
      DO I=INB(IP), INF(IP)
      II=I-I0
      DO K=1,3
      B1=B1+F(K,II)**2
      B2=B2+F(K,II)*G(K,II)
      ENDDO
      ENDDO
      IF (IT.EQ.1) THEN
        RLAM=-B1/(2.D0*A11)
	RMU=0.D0
      ELSE
        DET=-0.5D0/(A11*A22-A12**2)
        RLAM=(B1*A22-B2*A12)*DET
        RMU=(B2*A11-B1*A12)*DET
      END IF
      DO I=INB(IP),INF(IP)
      II=I-I0
      DO K=1,3
      G(K,II)=RLAM*F(K,II)+RMU*G(K,II)
      VEL(K,I)=VEL(K,I)+G(K,II)
      ENDDO
      ENDDO
c!
      IF (DABS(FTL-FTLL).LT. 1.D-5*FTL .and. it.gt.50) GOTO 112
c      print *, ' it=', it, ' ftl=', ftl
      FTLL=FTL
      GOTO 200
112   itvta=itvta+it
111   CONTINUE
      itvta=itvta/NP
      fluxa=fluxa/stot
      print *, ' itvtmax=', itvtmax, ' itvta=', itvta
      print *, ' fluxmin=', fluxmin, ' fluxmax=',
     >fluxmax
      print *, ' fluxa=', fluxa
      RETURN
      END
