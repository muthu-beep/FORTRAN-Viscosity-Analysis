      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'NTMAX.h'
      INCLUDE 'NEIBMAX.h'
      INCLUDE 'NPCH.h'

      PARAMETER(NVSMAX=NTMAX/2+2, NVTMAX=120000)
      PARAMETER (NEMAX=3*(NTMAX/2) )

      DIMENSION FT2(NVSMAX),IN0(6,NVSMAX),IN(NEIBMAX,NVSMAX),
     >NNEI(NVSMAX),PHI(NVTMAX),XNA(3,NVTMAX),
     >IV(3,NTMAX),V(3,NVTMAX),V1(3,NVTMAX),VP(3,NVTMAX),U(3,NVSMAX),
     >QW(3,NVTMAX),QWP(3,NVTMAX),QM(NVSMAX),VEL(3,NVSMAX)

      DIMENSION XV(3,NTMAX),YV(3,NTMAX),ZV(3,NTMAX),XP(3)

      dimension im(NEMAX),JM(NEMAX),NT1IJ(NEMAX),
     >K1IJ(NEMAX),NT2IJ(NEMAX),K2IJ(NEMAX),NNTI(NVSMAX),
     >NTI(6,NVSMAX),DX(3),DY(3),XN(3)

      NTS=2160

      DTP=1.D20
      NVT=0
      RAD=0.3D0
         xp(1)=-2.d0
         xp(2)=0.5d0
         xp(3)=0.d0

c       QW=0.D0
c       QWP=0.D0

       NT1=NTS
      IF (NT1.EQ.1280) THEN
      CALL ICO (XV,YV,ZV,NT1)
      CALL REFINE (3,NT1,XV,YV,ZV)
      END IF
      IF (NT1.EQ.2160) THEN
      CALL DDC (XV,YV,ZV,NT1)
      CALL REFINE (1,NT1,XV,YV,ZV)
      CALL HLREFINE (3,NT1,XV,YV,ZV)
      END IF
      IF (NT1.EQ.3840) THEN
      CALL DDC (XV,YV,ZV,NT1)
      CALL REFINE (3,NT1,XV,YV,ZV)
      END IF
      IF (NT1.EQ.5120) THEN
      CALL ICO (XV,YV,ZV,NT1)
      CALL REFINE (4,NT1,XV,YV,ZV)
      END IF
      IF (NT1.EQ.6000) then
      CALL DDC (XV,YV,ZV,NT1)
      CALL REFINE (1,NT1,XV,YV,ZV)
      CALL HLREFINE (5,NT1,XV,YV,ZV)
      END IF
      IF (NT1.EQ.8640) THEN
      CALL DDC (XV,YV,ZV,NT1)
      CALL REFINE (2,NT1,XV,YV,ZV)
      CALL HLREFINE (3,NT1,XV,YV,ZV)
      END IF

      NVTTEMP=0
      CALL MMESHFST (NT1,XV,YV,ZV,NVTTEMP,V,IV,IN0)

C     ORIENTATION
      
      DO JT=1,NTS
      I=IV(1,JT)
      J=IV(2,JT)
      K=IV(3,JT)
      DO L=1,3
      DX(L)=V(L,J)-V(L,I)
      DY(L)=V(L,K)-V(L,I)
      ENDDO
      XN(1)=DX(2)*DY(3)-DX(3)*DY(2)
      XN(2)=DX(3)*DY(1)-DX(1)*DY(3)
      XN(3)=DX(1)*DY(2)-DX(2)*DY(1)
      P=0.D0
      DO L=1,3
      P=P+XN(L)*V(L,I)
      ENDDO
      IF (P.LT.0.D0) THEN
	IV(2,JT)=K
	IV(3,JT)=J
      END IF
      ENDDO

      NVS=NTS/2+2
C     ADDL ARRAYS NECESSARY FOR RECONNECTIONS:
      DO I=1,NVS
      NNTI(I)=0
      ENDDO
      DO JT=1,NTS
      DO K=1,3
      I=IV(K,JT)
      NNTI(I)=NNTI(I)+1
      NTI(NNTI(I),I)=JT
      ENDDO
      ENDDO

      JE=0
      DO JT=1,NTS
      DO 44 K=1,3
      K2=K+1
      IF (K2.GT.3) K2=K2-3
      I=IV(K,JT)
      J=IV(K2,JT)
      IF(I.GT.J) GOTO 44
      JE=JE+1
      IM(JE)=I
      JM(JE)=J
      NT1IJ(JE)=JT
      K3=K+2
      IF (K3.GT.3) K3=K3-3
      K1IJ(JE)=IV(K3,JT)
      DO 45 MI=1,NNTI(I)
      JT2=NTI(MI,I)
      IF (JT2.EQ.JT)GOTO 45
      DO 46 MJ=1,NNTI(J)
      IF (NTI(MJ,J).EQ. JT2) GOTO 47
46    CONTINUE
45    CONTINUE
      PRINT *, ' ERROR'
      STOP
47    NT2IJ(JE)=JT2
      DO M=1,3
      IF (IV(M,JT2).NE.I .AND. IV(M,JT2).NE.J) K2IJ(JE)= IV(M,JT2)
      ENDDO
44    CONTINUE
      ENDDO
      print *, 'ne=', je

c --------------------------------------------------------------------
      DO I=1,NVS
      DO K=1,6
      IN(K,I)=IN0(K,I)
      ENDDO
      ENDDO

      XP(3)=0.D0
      DO I=1,NVS
      DO K=1,3
      XNA(K,I)=V(K,I)
      V(K,I)=XP(K)+V(K,I)*RAD
      VP(K,I)=V(K,I)
      VEL(K,I)=0.D0
      ENDDO
      ENDDO
      DO I=1,NVS
      NNEI(I)=6
      IF (IN0(6,I).EQ.0) NNEI(I)=5
      ENDDO
      T=0.d0
      DTL=1.d20
      DTP=1.d20
      NVT=0
      NVTP=0

      OPEN(1,FILE='restart.dat', form='UNFORMATTED')
      WRITE(1) T,DTL,DTP,NTS,IV,NNEI,IN,NVT,NVTP,RAD,V,VP,XNA,QW,QWP,
     >VEL,IM,JM,K1IJ,K2IJ,NT1IJ,NT2IJ
      CLOSE(1)
      stop
      END
