c     differes from our standard mesh.f in in being (6,*); also hlrefine
c      included
      SUBROUTINE DDC(XV,YV,ZV,NT)
      IMPLICIT REAL *8(A-H,O-Z)
      DIMENSION XV(3,*),YV(3,*),ZV(3,*),X(20), Y(20), Z(20)
      CALL ICO (XV,YV,ZV,NT)
      DO K=1,20
      XC=XV(1,K)+XV(2,K)+XV(3,K)
      YC=YV(1,K)+YV(2,K)+YV(3,K)
      ZC=ZV(1,K)+ZV(2,K)+ZV(3,K)
      D=DSQRT (XC**2+YC**2+ZC**2)
      X(K)=XC/D
      Y(K)=YC/D
      Z(K)=ZC/D
      ENDDO
      SQ=DSQRT (5.D0)
      A=4.D0/((1.D0+SQ)*DSQRT(3.D0))
      H=DSQRT (1.D0-0.25D0*A*A)
      RHO=A*DSQRT (0.25D0+SQ/10.D0)
      R=0.25D0*A*DSQRT (10.D0+22.D0/SQ)
      NT=0
      DO I=1,19
      DO 1 J=I+1,20
      DIST=(X(I)-X(J))**2+(Y(I)-Y(J))**2+(Z(I)-Z(J))**2
      IF (DABS(DIST-A*A).GT.1.D-6) GOTO 1
      AX=Y(I)*Z(J)-Z(I)*Y(J)
      AY=Z(I)*X(J)-X(I)*Z(J)
      AZ=X(I)*Y(J)-Y(I)*X(J)
      D=R*RHO/(H*DSQRT(AX**2+AY**2+AZ**2))
      AX=D*AX
      AY=D*AY
      AZ=D*AZ
      P=0.5D0*(R/H)**2
      X0=(X(I)+X(J))*P
      Y0=(Y(I)+Y(J))*P
      Z0=(Z(I)+Z(J))*P
      DO L=-1,1,2
      NT=NT+1
      XV(1,NT)=X(I)
      YV(1,NT)=Y(I)
      ZV(1,NT)=Z(I)
      XV(2,NT)=X(J)
      YV(2,NT)=Y(J)
      ZV(2,NT)=Z(J)
      XV(3,NT)=X0+L*AX
      YV(3,NT)=Y0+L*AY
      ZV(3,NT)=Z0+L*AZ
      D=DSQRT (XV(3,NT)**2+YV(3,NT)**2+ZV(3,NT)**2)
cc      PRINT *, ' D=', D, ' R=', R
      XV(3,NT)=XV(3,NT)/D
      YV(3,NT)=YV(3,NT)/D
      ZV(3,NT)=ZV(3,NT)/D
      ENDDO
1     CONTINUE
      ENDDO
      RETURN
      END
      SUBROUTINE ICO(XV,YV,ZV,NT)
      IMPLICIT REAL *8(A-H,O-Z)
      parameter (PI=3.14159265358979324D0)
      DIMENSION XV(3,*),YV(3,*),ZV(3,*)
      a=2.d0*(5.d0+dsqrt(5.d0))
      A=4.D0/DSQRT(A)
      R=5.D-1*A/DSIN(PI/5.D0)
      H=DSQRT(A*A-R*R)
      DO 1 NT=1,5
      XV(1,NT)=0.D0
      YV(1,NT)=0.D0
      ZV(1,NT)=1.D0
      FI=2.D0*PI*(NT-1)/5.D0
      XV(2,NT)=R*DCOS(FI)
      YV(2,NT)=R*DSIN(FI)
      ZV(2,NT)=1.D0-H
      FI=FI+2.D0*PI/5.D0
      XV(3,NT)=R*DCOS(FI)
      YV(3,NT)=R*DSIN(FI)
      ZV(3,NT)=1.D0-H
        DO K=2,3
	XV(K,NT+5)=XV(K,NT)
	YV(K,NT+5)=YV(K,NT)
	ZV(K,NT+5)=ZV(K,NT)
	ENDDO
      FI=FI-PI/5.D0
      XV(1,NT+5)=R*DCOS(FI)
      YV(1,NT+5)=R*DSIN(FI)
      ZV(1,NT+5)=H-1.D0
1     CONTINUE
        DO  NT=1,10
	DO K=1,3
	XV(K,NT+10)=-XV(K,NT)
	YV(K,NT+10)=-YV(K,NT)
	ZV(K,NT+10)=-ZV(K,NT)
	ENDDO
	ENDDO
      NT=20
      RETURN
      END
      SUBROUTINE REFINE (NR,NT,XV,YV,ZV)
      IMPLICIT REAL *8(A-H,O-Z)
      PARAMETER (NTMAX=160000)
      DIMENSION XV(3,*),YV(3,*),ZV(3,*),
     >XN(3,NTMAX),YN(3,NTMAX),ZN(3,NTMAX)
      common/JNK/ XN,YN,ZN
      DO 2 L=1,NR
      DO 3 N=1,NT
      IND=4*(N-1)
      DO 4 I=1,3
      M=IND+I
      XN(I,M)=XV(I,N)
      YN(I,M)=YV(I,N)
      ZN(I,M)=ZV(I,N)
        DO J=I+1,I+2
	K=J
	IF(K.GT.3) K=K-3
	XN(K,M)=(XV(I,N)+XV(K,N))
	YN(K,M)=(YV(I,N)+YV(K,N))
	ZN(K,M)=(ZV(I,N)+ZV(K,N))
        p=dsqrt(xn(k,m)**2+yn(k,m)**2+zn(k,m)**2)
        xn(k,m)=xn(k,m)/p
        yn(k,m)=yn(k,m)/p
        zn(k,m)=zn(k,m)/p
	ENDDO
4     CONTINUE
        DO K=1,2
	XN(K,M+1)=XN(K,M)
	YN(K,M+1)=YN(K,M)
	ZN(K,M+1)=ZN(K,M)
	ENDDO
      XN(3,M+1)=XN(1,M-1)
      YN(3,M+1)=YN(1,M-1)
      ZN(3,M+1)=ZN(1,M-1)
3     continue
      NT=4*NT
	DO N=1,NT
	DO K=1,3
	XV(K,N)=XN(K,N)
	YV(K,N)=YN(K,N)
	ZV(K,N)=ZN(K,N)
	ENDDO
	ENDDO
2     continue
      RETURN
      END
      SUBROUTINE MMESH(NT,XV,YV,ZV,NVT,V,IV,IN)
      IMPLICIT REAL *8(A-H,O-Z)
      DIMENSION XV(3,*),YV(3,*),ZV(3,*),
     >v(3,*),iv(3,*),in(6,*)
      inn=0
      do 300, i=1,NT
      do 200, k=1,3
      do 100, j=NVT+1,NVT+inn
      dd=(v(1,j)-xv(k,i))**2+(v(2,j)
     >-yv(k,i))**2+(v(3,j)-zv(k,i))**2
      if (dd.lt.1.d-15) goto 200
 100  continue
      inn=inn+1
      v(1,NVT+inn)=xv(k,i)
      v(2,NVT+inn)=yv(k,i)
      v(3,NVT+inn)=zv(k,i)
 200  continue
 300  continue
      do 111 n=1,NT
      do 112 k=1,3
      do 113 l=NVT+1,NVT+inn
      dd=(v(1,l)-xv(k,n))**2+(v(2,l)
     >-yv(k,n))**2+(v(3,l)-zv(k,n))**2
      if (dd.lt.1.d-15) then
         iv(k,n)=l
         goto 112
      endif
 113  continue
 112  continue
 111  continue
      DO 1 I=NVT+1,NVT+INN
      DO 1 K=1,6
      IN(K,I)=0
1     CONTINUE
      DO 4 II=1,NT
      DO 4 K=1,3
      I=IV(K,II)
      DO 2 L=K+1,K+2
      L1=L
      IF (L1.GT.3) L1=L1-3
      J=IV(L1,II)
      DO M=1,6
      IF (IN(M,I).EQ.J) GOTO 2
      IF (IN(M,I).EQ.0) GOTO 3
      ENDDO
3     IN(M,I)=J
2     CONTINUE
4     CONTINUE
      NVT=NVT+INN
      RETURN
      END


      SUBROUTINE HLREFINE(MR,NT,XV,YV,ZV)
      IMPLICIT REAL *8(A-H,O-Z)
      PARAMETER (NTMAX=160000)
      DIMENSION XV(3,*), YV(3,*), ZV(3,*),XN(3,NTMAX),YN(3,NTMAX),
     >ZN(3,NTMAX)
      COMMON /JNK/ xn,yn,zn
      IND=0
      DO 3 N=1,NT
      E2X=(XV(2,N)-XV(1,N))/MR
      E2Y=(YV(2,N)-YV(1,N))/MR
      E2Z=(ZV(2,N)-ZV(1,N))/MR
      E3X=(XV(3,N)-XV(1,N))/MR
      E3Y=(YV(3,N)-YV(1,N))/MR
      E3Z=(ZV(3,N)-ZV(1,N))/MR
      DO 4 I=0,MR-1
      DO 4 J=0,MR-1-I
      IND=IND+1
      XN(1,IND)=XV(1,N)+I*E2X+J*E3X
      YN(1,IND)=YV(1,N)+I*E2Y+J*E3Y
      ZN(1,IND)=ZV(1,N)+I*E2Z+J*E3Z
      XN(2,IND)=XN(1,IND)+E2X
      YN(2,IND)=YN(1,IND)+E2Y
      ZN(2,IND)=ZN(1,IND)+E2Z
      XN(3,IND)=XN(1,IND)+E3X
      YN(3,IND)=YN(1,IND)+E3Y
      ZN(3,IND)=ZN(1,IND)+E3Z
4     CONTINUE
      DO 5 J=0,MR-2
      DO 5 I=1,MR-1-J
      IND=IND+1
      XN(1,IND)=XV(1,N)+I*E2X+J*E3X
      YN(1,IND)=YV(1,N)+I*E2Y+J*E3Y
      ZN(1,IND)=ZV(1,N)+I*E2Z+J*E3Z
      XN(2,IND)=XN(1,IND)-E2X+E3X
      YN(2,IND)=YN(1,IND)-E2Y+E3Y
      ZN(2,IND)=ZN(1,IND)-E2Z+E3Z
      XN(3,IND)=XN(1,IND)+E3X
      YN(3,IND)=YN(1,IND)+E3Y
      ZN(3,IND)=ZN(1,IND)+E3Z
5     CONTINUE
3     CONTINUE
      NT=NT*MR**2
      DO N=1,NT
      DO K=1,3
      P=1.D0/DSQRT (XN(K,N)**2+YN(K,N)**2+ZN(K,N)**2)
      XV(K,N)=XN(K,N)*P
      YV(K,N)=YN(K,N)*P
      ZV(K,N)=ZN(K,N)*P
      ENDDO
      ENDDO
      RETURN
      END

      SUBROUTINE MMESHFST(NT,XV,YV,ZV,NVT,V,IV,IN)
C     FOR UNIT SPHERE ONLY!
      IMPLICIT REAL *8(A-H,O-Z)
      PARAMETER (NTMAX=160000,NDMAX=50)
      DIMENSION XV(3,*),YV(3,*),ZV(3,*),
     >v(3,*),iv(3,*),in(6,*),HOC(0:NDMAX,0:NDMAX,0:NDMAX),
     >LL(NTMAX/2+2)
      INTEGER HOC
      IS2=23
      JS2=13
      CALL GR2(IS2,JS2,DR)
      A=-1.D0-DR/1000.D0
      CALL GR2(IS2,JS2,DR)
      B=1.D0+DR/1000.D0
      TEMP=NT
      ND=DSQRT(TEMP)
      IF(ND.GT.NDMAX) ND=NDMAX
      DO K1=0,ND-1
      DO K2=0,ND-1
      DO K3=0,ND-1
      HOC(K1,K2,K3)=0
      ENDDO
      ENDDO
      ENDDO
      H1=ND/(B-A)

      inn=0
      do 300, i=1,NT
      do 200, k=1,3
      K1=(XV(K,I)-A)*H1
      K2=(YV(K,I)-A)*H1
      K3=(ZV(K,I)-A)*H1
      J=HOC(K1,K2,K3)
30    IF (J.EQ.0) GOTO 10
      dd=(v(1,j+NVT)-xv(k,i))**2+(v(2,j+NVT)
     >-yv(k,i))**2+(v(3,j+NVT)-zv(k,i))**2
      if (dd.lt.1.d-15) goto 20
      J=LL(J)
      GOTO 30
10    INN=INN+1
      J=INN
      LL(J)=HOC(K1,K2,K3)
      HOC(K1,K2,K3)=J
      v(1,J+NVT)=xv(k,i)
      v(2,J+NVT)=yv(k,i)
      v(3,J+NVT)=zv(k,i)
20    IV(K,I)=J+NVT
200   CONTINUE
300   CONTINUE

      DO 1 I=NVT+1,NVT+INN
      DO 1 K=1,6
      IN(K,I)=0
1     CONTINUE
      DO 4 II=1,NT
      DO 4 K=1,3
      I=IV(K,II)
      DO 2 L=K+1,K+2
      L1=L
      IF (L1.GT.3) L1=L1-3
      J=IV(L1,II)
      DO M=1,6
      IF (IN(M,I).EQ.J) GOTO 2
      IF (IN(M,I).EQ.0) GOTO 3
      ENDDO
3     IN(M,I)=J
2     CONTINUE
4     CONTINUE
      NVT=NVT+INN
      RETURN
      END

c	Generate random numbers using algorithm
c	k(n+1)=k(n)*5**17 mod 2**40
c	NDP - Fortran 1.4e
	subroutine gr2(i,j,dr)
c					k=i+j*2**31
        implicit integer*4 (a-c), integer*4 (u-z), real*8 (d-h,o-t)
        save key,mask15,mask10,two40,two9,b0,b1,b2
c				Initial settings at 1st calling
	if(key.eq.0)then
		mask15=2**15-1
		mask10=2**10-1
		two=2
		one=1
		two40=one/two**40
		two9=one/two**9
c
		ii=582758085
		jj=355
		b0=iand(ii,mask15)
		x=rshift(ii,15)
		b1=iand(x,mask15)
		x=rshift(x,15)
		b2=2*jj+x
c
		key=1
	endif
c				Decomposition 'k' onto a0,a1,a2
	a0=iand(i,mask15)
	x=rshift(i,15)
	a1=iand(x,mask15)
	x=rshift(x,15)
	a2=2*j+x
c				Multiplication 'a..' on 'b..' and getting mod 2**40
c				The result is 'c..'
	x=a0*b0
	c0=iand(x,mask15)
	x=rshift(x,15)+b1*a0
	v0=iand(x,mask15)
	v1=rshift(x,15)
	x=v0+a1*b0
	c1=iand(x,mask15)
	x=rshift(x,15)+v1+a1*b1
	y=a2*b0+b2*a0
	x=iand(x,mask10)
	y=iand(y,mask10)
	x=x+y
	c2=iand(x,mask10)
c					Making ' i' and 'j' from 'c..'
	x=c2
	j=rshift(x,1)
	x=lshift(iand(x,1),30)
	y=c1
	y=lshift(y,15)
	x=ior(x,y)
	y=c0
	i=ior(x,y)
c				Making random number 'dr'
	dr=j*two9+i*two40
	return
	end
