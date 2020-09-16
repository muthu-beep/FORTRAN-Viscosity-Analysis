	SUBROUTINE VELOCITY(Y,U)
	IMPLICIT REAL*8 (A-H,O-Z)
	INCLUDE 'NPCH.h'
	PARAMETER (NP=NPCH,NVTMAX=100000)
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

