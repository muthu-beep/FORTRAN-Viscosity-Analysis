c     much preferable version of JIM than mgimfast-stuff
      SUBROUTINE GIM (INB,INF,V,IN,nnei,FIT,XNA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      include 'NEIBMAX.h'
      DIMENSION FIT(3,*),XNA(3,*),IN(neibmax,*),nnei(*),
     >GRAD(3,0:3),XNOLD(3),V(3,*),
     >XN(3,0:3), XNN(3,0:3),G(3,3),FUN(0:3),T(3,4),DX(3),DIFF(3)
      logical FAILGIM
      COMMON/FAILGIM/FAILGIM
      common/ITGIMMAX/ITGIMMAX
      common/pmin/pmin
      nneimax=0
      do I=INB,INF
      NNEIMAX= MAX(NNEIMAX,nnei(i))
      ENDDO
      COEFF=1.d0
      IF(NNEIMAX.GE.10) COEFF=0.d0
      ITGIM=0
2     ITGIM =ITGIM+1
      ERRMAX=0.D0
      DO 1 I=INB, INF
      FAILGIM=.FALSE.
465   ITGRAD=0
      DO K=1,3
      XN(K,1)=XNA(K,I)
      XNOLD(K)=XNA(K,I)
      ENDDO
      CALL MGRADGIM (.FALSE.,COEFF,I,XN(1,1),V,NV,IN,nnei,XNA,
     >FIT(1,I),FUN(1),GRAD(1,1))
462   ITGRAD=ITGRAD+1
      STEP=0.05
464   P=0.D0
      DO K=1,3
       XN(K,0)=XN(K,1)-STEP*GRAD(K,1)
       P=P+XN(K,0)**2
      ENDDO
      P=1.D0/DSQRT (P)
      DIST2=0.D0
      DO K=1,3
      XNA(K,I)=XN(K,0)*P
      XN(K,0)=XNA(K,I)
      DX(K)=XN(K,1)-XN(K,0)
      DIST2=DIST2+DX(K)**2
      ENDDO

      CALL MGRADGIM (.FALSE.,COEFF,I,XN(1,0),V,NV,IN,nnei,XNA,
     >FIT(1,I),FUN(0),GRAD(1,0))
      RES=DSQRT (GRAD(1,0)**2+GRAD(2,0)**2+GRAD(3,0)**2)
      IF (RES.LT.1.D-5) GOTO 463
      
c      IF (FAIL) THEN
c        DO K=1,3
c        XN(K,1)=XN(k,0)
c        GRAD(k,1)=GRAD(k,0)
c        ENDDO
c        GOTO 462
c      END IF

      IF (FUN(0).GT.FUN(1) .OR. DIST2.GT.1.D0) THEN
        STEP=STEP*0.5D0
       GOTO 464
       END IF
       H=0.5D0*DIST2
       SQ=DSQRT (3.D0/(1.D0-0.25D0*DIST2))
       DIFF(1)=SQ*(DX(2)*XN(3,0)-DX(3)*XN(2,0))
       DIFF(2)=SQ*(DX(3)*XN(1,0)-DX(1)*XN(3,0))
       DIFF(3)=SQ*(DX(1)*XN(2,0)-DX(2)*XN(1,0))
       P2=0.D0
       P3=0.D0
       DO K=1,3
       TEMP=3.D0*H*XN(K,0)+DX(K)
       XN(K,2)=XN(K,0)+0.5D0*(DIFF(K)-TEMP)
       XN(K,3)=XN(K,0)+0.5D0*(-DIFF(K)-TEMP)
       P2=P2+XN(K,2)**2
       P3=P3+XN(K,3)**2
       ENDDO
       P2=1.D0/DSQRT (P2)
       P3=1.D0/DSQRT (P3)
       DO K=1,3
       XN(K,2)=XN(K,2)*P2
       XN(K,3)=XN(K,3)*P3
       ENDDO

       MMIN=0
       FMIN=FUN(0)
       DO M=2,3
      CALL MGRADGIM (.TRUE.,COEFF,I,XN(1,M),V,NV,IN,nnei,XNA,
     >FIT(1,I),FUN(M),GRAD(1,M))
      IF (FUN(M).LT.FMIN) THEN
        FMIN=FUN(M)
        MMIN=M
      END IF
      ENDDO
      P=0.D0
      DO K=1,3
      P=P-DX(K)*XN(K,0)
      ENDDO
      DO K=1,3
      G(K,1)=-DX(K)-P*XN(K,0)
      G(K,3)=XN(K,0)
      ENDDO
      G(1,2)=G(2,3)*G(3,1)-G(3,3)*G(2,1)
      G(2,2)=G(3,3)*G(1,1)-G(1,3)*G(3,1)
      G(3,2)=G(1,3)*G(2,1)-G(2,3)*G(1,1)
      DO M=1,2
      P=0.D0
      DO K=1,3
      P=P+G(K,M)**2
      ENDDO
      P=1.D0/DSQRT(P)
      DO K=1,3
      G(K,M)=G(K,M)*P
      ENDDO
      ENDDO

      FX=0.D0
      FY=0.D0
      DO K=1,3
      FX=FX+GRAD(K,0)*G(K,1)
      FY=FY+GRAD(K,0)*G(K,2)
      DO M=1,3
      S=0.D0
      DO L=1,3
      S=S+G(L,K)*XN(L,M)
      ENDDO
      XNN(K,M)=S
      ENDDO
      ENDDO
      DO M=1,3
      DO K=1,2
      XNN(K,M)=XNN(K,M)/XNN(3,M)
      ENDDO
      ENDDO

      DO M=1,3
      T(M,1)=XNN(1,M)**2
      T(M,2)=XNN(1,M)*XNN(2,M)
      T(M,3)=XNN(2,M)**2
      T(M,4)=FUN(M)-FUN(0)-FX*XNN(1,M)-FY*XNN(2,M)
      ENDDO
      CALL  GM34(3,4,T)
      C=T(1,4)
      D=T(2,4)
      EE=T(3,4)
      IF (C+EE.LT.0.D0 .OR. 4.D0*C*EE-D**2 .LE. 0.D0) THEN
      IF (MMIN.NE.0)
     >CALL MGRADGIM (.FALSE.,COEFF,I,XN(1,MMIN),V,NV,IN,nnei,XNA,
     >FIT(1,I),FUN(MMIN),GRAD(1,MMIN))
      FUN(1)=FMIN
	DO K=1,3
	XN(K,1)=XN(K,MMIN)
	GRAD(K,1)=GRAD(K,MMIN)
	ENDDO
        GOTO 462
      END IF
      DET=1.D0/(4.D0*C*EE-D**2)
      XNN(1,1)=(D*FY-2.D0*EE*FX)*DET
      XNN(2,1)=(D*FX-2.D0*C*FY)*DET
      P=1.D0/DSQRT (XNN(1,1)**2+XNN(2,1)**2+1.D0)
      XNN(1,1)=XNN(1,1)*P
      XNN(2,1)=XNN(2,1)*P
      XNN(3,1)=P
      DO K=1,3
      S=0.D0
      DO L=1,3
      S=S+G(K,L)*XNN(L,1)
      ENDDO
      XN(K,1)=S
      ENDDO
      P=0.D0
      DO K=1,3
      P=P+XN(K,1)**2
      ENDDO
      P=1.D0/DSQRT (P)
      tst=0.d0
      DO K=1,3
      XN(K,1)=XN(K,1)*P
      XNA(K,I)=XN(K,1)
      tst=tst+ (xn(k,1)-xn(k,mmin))**2
      ENDDO
      CALL MGRADGIM (.FALSE.,COEFF,I,XN(1,1),V,NV,IN,nnei,XNA,
     >FIT(1,I),FUN(1),GRAD(1,1))
      IF (FUN(1).GT.FMIN .or. TST.GT. DIST2) THEN
c      print *, ' unusual event in mgimsh', FUN(1),FMIN
c      print *, ' i=', i
c      open (17, file='FEND')
c       write (4,*) ' unusual event in mgimsh',FUN(1),FMIN
c      stop
c      FAILGIM= .true.
c      return

      IF (MMIN.NE.0)
     >CALL MGRADGIM (.FALSE.,COEFF,I,XN(1,MMIN),V,NV,IN,nnei,XNA,
     >FIT(1,I),FUN(MMIN),GRAD(1,MMIN))
       FUN(1)=FMIN
	DO K=1,3
	XN(K,1)=XN(K,MMIN)
	GRAD(K,1)=GRAD(K,MMIN)
	ENDDO
        GOTO 462
      END IF

      RES=DSQRT (GRAD(1,1)**2+GRAD(2,1)**2+GRAD(3,1)**2)
      IF (RES.GT.1.D-5) GOTO 462
463   continue
      ERR=0.D0
      p=0.d0
      DO K=1,3
      ERR=ERR+(XNA(K,I)-XNOLD(K))**2
      p=p+ xna(k,i)*xnold(k)
      ENDDO
      pmin= min(pmin, p)
      if (p.lt.0.4d0) then
       print *, ' p=', p, ' i=', i
       FAILGIM= .true.
       return
      end if
      IF (ERR.GT.ERRMAX) ERRMAX=ERR
1     CONTINUE
      ERRMAX=DSQRT (ERRMAX)
C      PRINT *, ' ITGIM=', ITGIM, ' ERRMAX=', ERRMAX
      IF (ERRMAX.GT. 1.D-4) GOTO 2
      itgimmax=max(itgimmax,itgim)
      RETURN
      END

      SUBROUTINE MGRADGIM (FAST,COEFF,II,XNI,V,NV,IN,nnei,XNA,
     >FIT,FUN,GRAD)
      IMPLICIT REAL *8 (A-H,O-Z)
      include 'NEIBMAX.h'
      DIMENSION V(3,*),IN(neibmax,*),XNA(3,*),T(3,4),X(neibmax,3),
     >E(3,3),Y(3,neibmax),TT(4),S(3),GRAD(*),FIT(*),Z(3),DM(neibmax),
     >XNI(*),DX(3,neibmax),nnei(*)
      LOGICAL FAST
c      common/gimcount/icount
c      icount=icount+1
      DO I=1,3
      DO J=1,4
      T(I,J)=0.D0
      ENDDO
      ENDDO
      A=XNI(1)
      B=XNI(2)
      C=XNI(3)
      E(3,1)=A
      E(3,2)=B
      E(3,3)=C
      IF(A*A .LE. B*B .AND. A*A .LE. C*C) THEN
      E(1,1)=0.D0
      E(1,2)=C
      E(1,3)=-B
      ELSE
	IF(B*B .LE. A*A .AND. B*B .LE. C*C )THEN
	 E(1,1)=C
	 E(1,2)=0.D0
	 E(1,3)=-A
	ELSE
	 E(1,1)=B
	 E(1,2)=-A
	 E(1,3)=0.D0
	END IF
      END IF

      E(2,1)=E(1,2)*C - E(1,3)*B
      E(2,2)=E(1,3)*A - E(1,1)*C
      E(2,3)=E(1,1)*B - E(1,2)*A
      DO K=1,3
      P=E(K,1)**2+E(K,2)**2+E(K,3)**2
      P=1.D0/DSQRT(P)
      DO L=1,3
      E(K,L)=E(K,L)*P
      ENDDO
      ENDDO

      DO 1 J=1, nnei(ii)
      jj=in(j,iI)
      DO K=1,3
      DX(K,J)=V(K,JJ)-V(K,II)
      ENDDO
      DO K=1,3
      Q=0.D0
      DO L=1,3
      Q=Q+E(K,L)*DX(L,J)
      ENDDO
      X(J,K)=Q
      ENDDO
      Y(1,J)=X(J,1)**2
      Y(2,J)=X(J,1)*X(J,2)
      Y(3,J)=X(J,2)**2
      DM(J)=1.D0/(Y(1,J)+Y(3,J)+X(J,3)**2)
      DO K=1,3
      TT(K)=Y(K,J)*DM(J)
      ENDDO
      TT(4)=X(J,3)*DM(J)

      DO K=1,3
      DO L=K,4
      T(K,L)=T(K,L)+TT(L)*Y(K,J)
      ENDDO
      ENDDO
1     CONTINUE

      DO K=1,3
      DO L=1, K-1
      T(K,L)=T(L,K)
      ENDDO
      ENDDO
      CALL GM34(3,4,T)
      FIT(1)=T(1,4)
      FIT(2)=T(2,4)
      FIT(3)=T(3,4)

      DO K=1,3
      Z(K)=0.D0
      S(K)=0.D0
      GRAD(K)=0.D0
      ENDDO
      FUN=0.D0
      P=DSQRT (XNA(1,II)**2+XNA(2,II)**2+XNA(3,II)**2)
      DO 2 J=1,nnei(ii)
      JJ=IN(J,II)
      F=-X(J,3)
      Q=0.D0
      DO K=1,3
      F=F+FIT(K)*Y(K,J)
      Q=Q+DX(K,J)*(XNI(K)+XNA(K,JJ))
      ENDDO
      IF (FAST) THEN
        FUN=FUN+DM(J)*(F**2+COEFF*Q**2)
	GOTO 2
      END IF
      F1=F*DM(J)
      Q1= Q*DM(J)*COEFF
      FUN=FUN+F1*F+Q1*Q
      F2=F1*X(J,3)
      DO K=1,2
      Z(K)=Z(K)-F1*X(J,K)
      S(K)=S(K)+F2*X(J,K)
      ENDDO
      DO K=1,3
      GRAD(K)= GRAD(K)+ Q1*DX(K,J)
      ENDDO
2     CONTINUE
      Z(1)=(Z(1) - 2.D0*S(1)*FIT(1) - S(2)*FIT(2))/P
      Z(2)=(Z(2)- S(1)*FIT(2) - 2.D0*S(2)*FIT(3))/P
      DO K=1,3
      DO L=1,3
      GRAD(K)=GRAD(K)+E(L,K)*Z(L)
      ENDDO
      ENDDO

      P=0.D0
      DO K=1,3
      P=P+GRAD(K)*XNA(K,II)
      ENDDO
      DO K=1,3
      GRAD(K)=2.D0*(GRAD(K)-P*XNA(K,II))
      ENDDO
      RETURN
      END

      SUBROUTINE GM34(N,M,T)
      REAL *8 T(3,4),R,U
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
