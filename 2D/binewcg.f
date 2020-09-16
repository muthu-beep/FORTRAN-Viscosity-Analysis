c	**************************************************************************
c	* Pinched-Flow Fractionation with 3 inputs and fixed points              *
c	* Program BIVEL                                                          *
c	* Access program TRAJ                                                    *
c	**************************************************************************
      implicit real*8 (a-h,O-Z)
      PARAMETER (NP=14,
     >NVTMAX=100000)
      DIMENSION V(2,NVTMAX),
     >PHI(NVTMAX),XNA(2,NVTMAX),
     >W(2,NVTMAX),WPR(2,NVTMAX),
     >DL(2,NVTMAX),SX(2,NVTMAX),SY(2,NVTMAX),
     >XC(2),XCW(2),INB(NP),
     >INF(NP),DSW(NP),XNW(2,NP),UINF(2,NVTMAX)
      DIMENSION L1(NP), X0(2,NP), N1(NP),E1(2,NP)
      DIMENSION DX(2),U(2),OM(2),SM(2),QF(NP),X1(2),
     >Q(2,NVTMAX),Y(2),SUM(2),jmin(nvtmax,np),STR(2,2)
      DIMENSION R(2,NVTMAX),P(2,NVTMAX),A(2,NVTMAX),
     >RT(2,NVTMAX),PT(2,NVTMAX),AT(2,NVTMAX),
     >SXT(2,NVTMAX),SYT(2,NVTMAX),WTPR(2,NVTMAX)

	REAL *8 L1
 	TOL=1.d-6
 	TOL=1.d-7
	PI=DACOS(-1.D0)
	C2=-2.D0/PI
	THETA=(PI/180)*30
	THETA2=(PI/180)*60
	ALPHA=(PI/180)*56
c	QF(1)=0.1D0
c	QF(4)=5.00D0

c       Inlet and outlet flows:

        QF(1)=1.d0
        QF(7)=-5.d0/7.d0

	QF(13)=-QF(1)-QF(7)

c     *************************************************************************
c	Routine BIVEL                                                           *
c	*************************************************************************
c	NP = Total number of Panels (for pff=13)
c	IP = Current Panel number
c	N1(IP) = The number of mesh segments on panel IP
c	L1(IP) = Length of panel IP
c	X0(1,IP), X0(2,IP) = Are coordinates of the corner point for panel
c	E1(1,IP), E1(2,IP) = Are coordinates of the unit director vector
c	                     of panel IP
c	INB(IP)= Is the first mesh node for panel number IP
c	INF(IP)= Is the last mesh node for panel IP
c	V = The array of collocation nodes (i.e. the midpoints at
c                          mesh node segments.
c	XNW(1,IP), XNW(2,IP) = Coordinates of outward unit normal to
c                          panel IP
c	DSW(IP)= The length of mesh segment for panel IP
c     UINF = The  U(infinity) that is UINF= 0 on physical walls
c	       UINF = Poiseuille flow at inlet and outlet streams
c	UINF(1,i), UINF(2,i) = Are coordinates (Uinfinity) of mesh
c	              node i.
c	C = Center of contour XCW(1),XCW(2)= Coordinates of center of
c	              contour
c	 Replace DEFLMAT with DM routine
C------References PANEL Subroutine
C------IP = 1,2,...NP
C------NP - total number of panels
C------N1(IP), N2(IP) assumes next panel number IP is divided
C------into N1 x N2 elements
C------L1(IP), L2(IP) are dimensions of panel 1 # IP
C------x0(1,IP), x0(2,IP), x0(3,IP) are coordinates at the corner point
C------of the panel # IP E1(1,IP), E1(2,IP), E1(3,IP) are coordinates
C------of unit vector for panel #IP along E2(1,IP), E2(2,IP), E2(3,IP)
C------are coordinates of unit vector for  panel # IP along 2nd axis
C------   Create commentary for complex geometries
C------ INB(IP), INF(IP) are first and last indices for the nodes on panel
C------ V(1,i), V(2,i), V(3,i) #IP for
C------ INB(IP) <= i <= INF(IP) are cartesian coordinatesfor the mesh
C------ nodes on panel # IP
C------ XNW(1,IP), XNW(2,IP), XNW(3,IP) is a unit outward normal to
C------ panel IP
C------ XNW = E1 x E2
C------ DSW(IP) is the area of mesh element on panel IP
c	**************************************************************************
c	*  New Panels for Pinched-Flow Device                                    *
c	*  Panels 1-3, Inlet, outlet panels for fluid                            *
c	*  Panel 1,2 Inlet Panels (upper (1) and lower (2))                      *
c	*  Panel 3   Outlet Panel for fluid                                      *
c	*  Panel 4-13, remainder of panels for pff device                        *
c	*  e1 moves counter-clockwise                                            *
c	*    N1...=number of elements                                            *
c	**************************************************************************
        h=0.1d0

        X0(1,1)=-6
        X0(2,1)=1

        X0(1,2)=-6
        X0(2,2)=0

        X0(1,3)=-h
        X0(2,3)=0
   
        X0(1,4)=(-h)+h*COS(PI/3.d0)
        X0(2,4)=(-h)+h*SIN(PI/3.d0)

        x0(1,5)=(-h)+h*COS(PI/6.d0)
        x0(2,5)=(-h)+h*SIN(PI/6.d0)

        x0(1,6)=0
        x0(2,6)=-h

        x0(1,7)=0
        x0(2,7)=-6

        x0(1,8)=1
        x0(2,8)=-6

        x0(1,9)=1
        x0(2,9)=-h

       x0(1,10)=(1+h)+h*COS(5.d0/6.d0*PI)
       x0(2,10)=-h+h*SIN(5.d0/6.d0*PI)

       x0(1,11)=(1+h)+h*COS(2.d0/3.d0*PI)
       x0(2,11)=-h+h*SIN(2.d0/3.d0*PI)

       x0(1,12)=1+h
       x0(2,12)=0

       x0(1,13)=7
       x0(2,13)=0

       x0(1,14)=7
       x0(2,14)=1

c       x0(1,15)=11
c       x0(2,15)=6

c       x0(1,16)=7
c       x0(2,16)=0

        OPEN (13, FILE='Tchan.dat')
        do ip=1,np
        write (13,*) x0(1,ip),x0(2,ip)
        enddo
        write (13,*) x0(1,1),x0(2,1)
        CLOSE (13)

c        stop

        OPEN (15, FILE='3Dchan.dat')
        WRITE (15,*)'VARIABLES = "X","Y","Z"'
        WRITE (15,321) 2*NP, NP
321     FORMAT ('ZONE N=', I5,'  E=',I4,' F=FEPOINT, ET=QUADRILATERAL')

        do ip=1,np
        ZERO= 0.D0
        ONE=1.D0
        write (15,*) x0(1,ip),x0(2,ip),ZERO
        write (15,*) x0(1,ip),x0(2,ip),ONE	
        enddo

        DO I=1,2*NP,2
        IF (I.EQ.2*NP-1)THEN
        n2=2
        n1=1
        WRITE (15,*) I, I+1, 2, 1
        ELSE
        WRITE (15,*) I, I+1,I+3, I+2
        END IF
        ENDDO

        CLOSE(15)

c        stop

c       For Nakashima chanel

c       OPEN (10, FILE='CHAN2D.dat')
c       READ(10,*) X0
c       CLOSE(10)

c       OPEN (13, FILE='Tchan.dat')
c       DO IP=1,NP
c       WRITE (13,*) X0(1,IP)/50.D0,X0(2,IP)/50.D0
c       ENDDO
c       CLOSE (13)

c	*************************************************************************
c	* Calculate E1 - coordinates of unit vector                             *
c	*************************************************************************
        DO IP=1,NP
        JP=IP+1
        IF (JP.GT.NP) JP=JP-NP
        TMP=0.D0
        DO K=1,2
        TMP=TMP+ (X0(K,JP)-X0(K,IP))**2
        ENDDO
c       scaled lenght of the panels
        L1(IP)=DSQRT(TMP)
        DO K=1,2
        E1(K,IP)= (X0(K,JP)-X0(K,IP))/L1(IP)
        ENDDO
        ENDDO

        do ip=1,np
        jp=ip+1
        IF (JP.GT.NP) JP=JP-NP
        temp=0.d0
        do k=1,2
        temp=temp+ e1(k,ip)*e1(k,jp)
        enddo
        print *, 'ip=', ip, ' prod=', temp, ' L1=', L1(ip)
        enddo
c	*************************************************************************
c	* comment out for now.                                                  *
c	DO IP=1,NP				                                                *
c	N1(IP)=3000                                                             *
c	ENDDO                                                                   *
c	*************************************************************************
         do ip=1,np
         n1(ip)= 1000
         enddo
         N1(1)=200
         N1(7)=200
         N1(13)=200

         N1(3)=10
         N1(4)=10
         N1(5)=10
         N1(9)=10
         N1(10)=10
         N1(11)=10

c	*************************************************************************
C	* Perform Call to PANEL subroutine                                      *
c	*************************************************************************
	I=0
	DO IP=1,NP
	Call PANEL(I,NVTMAX,N1(IP),L1(IP),X0(1,IP),
     >E1(1,IP),INB(IP),INF(IP),V,XNW(1,IP),DSW(IP))
	ENDDO
	NVT=I
	DO I=1,NVT
	DO K=1,2
	UINF(K,I)=0
	W(K,I)=0
	ENDDO
	ENDDO
c	**************************************************************************
c	* Loop for inlet and outlet panels                                       *
c	* Poiseuille flow                                                        *
c	* Fluxes QF1, QF2, QF3                                                   *
c	* QF1 + QF2 + QF3 =0                                                     *
c	* ( later - do Dimension Q(3) Q(IP)                                      *
c	* ( do dimension for Y(IP) = Y(3)                                        *
c	* QF1>0 QF2>0 QF3<0                                                      *
c	* QF(IP) = QF                                                            *
c	* L1(IP) = H                                                             *
c	* Y = Distance from channel wall                                         *
c	* H = Height of channel                                                  *
c	**************************************************************************
c	*  Give values for QF1, QF2, QF3                                         *
c	*  QF2 = 50xQF1                                                          *
c	*  QF3 = QF1 + QF2                                                       *
c	*  U = ((6*QF)/H^3)*(y(H-y))                                             *
c	*  Alternative for y:                                                    *
c	*  Y(IP)=((V(1,i)-X0(1,IP))**2)+                                         *
c	*        (V(2,i)-X0(2,IP))**2)**0.5                                      *
c	**************************************************************************
        DO IP=1,NP
        IF (IP.EQ.1 .OR. IP.EQ.7 .OR. IP.EQ.13) THEN
        DO I= INB(IP),INF(IP)
c       calculate y
        Y1=(I-INB(IP)+0.5D0)*DSW(IP)
        DO K=1,2
        UINF(K,I)=(((6.D0*(-QF(IP)))/((L1(IP))**3))*(Y1*(L1(IP)-Y1)))
     >*XNW(K,IP)
        ENDDO
        ENDDO
        END IF
        ENDDO

c-----Calculate Center of Contour (change 3->2 for 2D)
      DO K=1,2
      XCW(K)=0.D0
      ENDDO
      SW=0.D0
      DO IP=1,NP
      DO I=INB(IP),INF(IP)
      PHI(I)=DSW(IP)
      DO K=1,2
      XNA(K,I)= XNW(K,IP)
      XCW(K)=XCW(K)+V(K,I)*DSW(IP)
      ENDDO
      SW=SW+DSW(IP)
      ENDDO
      ENDDO

      DO K=1,2
      XCW(K)=XCW(K)/SW
      ENDDO
c	**************************************************************************
c	* DM routine                                                             *
c	**************************************************************************
	DM=0.D0
	DO IP=1,NP
	DO I=INB(IP),INF(IP)
	DO K=1,2
	DM=DM+((V(K,I)-XCW(K))**2)*DSW(IP)
	ENDDO
	ENDDO
	ENDDO
	DM=1.D0/DM
c	**************************************************************************
c	* Write a subroutine that calculates w' (RBP of w)                       *
c	**************************************************************************
      DO IP=1,NP
      DO JP=1,NP
      DO I=INB(IP),INF(IP)
      RRMIN=1.D20
      DO J=INB(JP),INF(JP)
      D=0.D0
      DO K=1,2
      D=D+ (V(K,J)-V(K,I))**2
      ENDDO
      IF (D.LT.RRMIN) THEN
        RRMIN=D
        JMIN(I,JP)=J
      END IF
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      DO I=1,NVT
      DO K=1,2
      W(K,I)=0.D0
      ENDDO
      ENDDO

c      read (*,*) icont
	Call RBP (XCW,V,W,PHI,NVT,DM,SW,WPR)
c	**************************************************************************
c	* Sum of :                                                               *
c	* (Uinf)k(y)                                                             *
c	* -2TKSUM  (integral of stress tensor dot n dot                          *
c	*          w(x)-w(y) dSx                                                 *
c	* -Wk(y)    (called W)                                                   *
c	* +WPRk(y)  (called WPR)                                                 *
c	* -nk(y)/SW (integral (W dot n) dS   (called WND)                        *
c     *                                                                        *
c	**************************************************************************
c	* WND Routine                                                            *
c	**************************************************************************
	WND=0.D0
	DO IP=1,NP
	DO I=INB(IP), INF(IP)
	DO K=1,2
	WND=WND+W(K,I)*XNW(K,IP)*DSW(IP)
	ENDDO
	ENDDO
	ENDDO
	WND=WND/SW
c	*************************************************************************
c	* End of calculate WND                                                  *
c	*************************************************************************
c	* Calculate W                                                           *
c	*************************************************************************
	DO IP=1,NP
	DO I=INB(IP),INF(IP)
	DO K=1,2
	R(K,I)=UINF(K,I)+WPR(K,I)-XNW(K,IP)*WND
	ENDDO
	ENDDO
	ENDDO

	DO IP=1,NP
	DO 33 JP=1,NP
	IF (JP.EQ.IP) GOTO 33
	DO K=1,2
	X1(K)=X0(K,JP)+L1(JP)*E1(K,JP)
	ENDDO
	DO I=INB(IP),INF(IP)
	CALL LINESTR (V(1,I),X0(1,JP),X1,XNW(1,JP),STR)
	DO K=1,2
	DO L=1,2
	R(K,I)=R(K,I)+C2*STR(K,L)*W(L,JMIN(I,JP))
	ENDDO
	ENDDO
        ENDDO
33      CONTINUE
        ENDDO

c	*************************************************************************
c	* SELF-INTERACTION BETWEEN WALLS                                        *
c	*************************************************************************
      DO IP=1,NP
      DO JP=IP+1,NP
      DO I=INB(IP), INF(IP)
      DO K=1,2
      SY(K,I)=0.D0
      ENDDO
      ENDDO

      DO J=INB(JP), INF(JP)
      DO K=1,2
      SX(K,J)=0.D0
      ENDDO
      ENDDO

      DO I=INB(IP), INF(IP)
      DO J=INB(JP), INF(JP)
      D=0.D0
	WRX=0.D0
	WRY=0.D0
	DO K=1,2
        DX(K)=V(K,J)-V(K,I)
	D=D+DX(K)**2
	WRX=WRX+(W(K,J)-W(K,JMIN(I,JP)))*DX(K)
	WRY=WRY-(W(K,I)-W(K,JMIN(J,IP)))*DX(K)
	ENDDO
	WRX=WRX/(D**2)
	WRY=WRY/(D**2)

      DO K=1,2
      SY(K,I)=SY(K,I)+WRX*DX(K)
      SX(K,J)=SX(K,J)+WRY*DX(K)
      ENDDO
      ENDDO
      ENDDO

      IY=INB(IP)
      IX=INB(JP)

      DO I=INB(IP),INF(IP)
      RNX=0.D0
      DO K=1,2
      RNX=RNX+(V(K,IX)-V(K,I))*XNW(K,JP)
      ENDDO
      RNX=RNX*DSW(JP)*C2
      DO K=1,2
      R(K,I)=R(K,I)+ RNX*SY(K,I)
      ENDDO
      ENDDO

      DO J=INB(JP),INF(JP)
      RNY=0.D0
      DO K=1,2
      RNY=RNY-(V(K,IY)-V(K,J))*XNW(K,IP)
      ENDDO
      RNY=RNY*DSW(IP)*C2
      DO K=1,2
      R(K,J)=R(K,J)+ RNY*SX(K,J)
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      DO I=1,NVT
      DO K=1,2
      R(K,I)=W(K,I)-R(K,I)
      P(K,I)=   R(K,I)
      RT(K,I)=  R(K,I)
      PT(K,I)=  RT(K,I)
      ENDDO
      ENDDO

      IT=0

777   IT=IT+1
      ERRMAX=0D0
      RR=0.D0
      ERRA=0.D0
      w2=0.d0
      DO I=1,NVT
      ERR=0.D0
      temp=0.d0
      DO K=1,2
      RR=RR+R(K,I)*RT(K,I)*PHI(I)
      ERR=ERR+R(K,I)**2
      temp=temp+w(k,i)**2
      ENDDO
      ERRA=ERRA+ERR*PHI(I)
      w2=w2+temp*phi(i)
      ERRMAX=MAX(ERRMAX,ERR)
      ENDDO
      W2= DSQRT(W2/SW)
      if (IT.gt.1) then
        ERRMAX=DSQRT(ERRMAX)/W2
        ERRA= DSQRT(ERRA/SW)/W2
        PRINT *,  ' IT=', IT,' ERRMAX=', ERRMAX, ' ERRA=', ERRA 
        print *, ' W2=', w2
        if (ERRA.lt. TOL)  goto 778
      end if

      DO I=1,NVT
      DO K=1,2
      A(K,I)=0.D0
      AT(K,I)=0.D0
      ENDDO
      ENDDO

      DO 44 IP=1,NP
      DO 45 JP=IP+1,NP

      IY=INB(IP)
      IX=INB(JP)

      DO I=INB(IP),INF(IP)
      RNX=0.D0
      DO K=1,2
      RNX=RNX+(V(K,IX)-V(K,I))*XNW(K,JP)
      ENDDO
      RNX=RNX*DSW(JP)*C2
      DO K=1,2
      SYT(k,i)= RNX*PT(K,I)*PHI(I)
      SY(K,I)=0.D0
      ENDDO
      ENDDO

      DO J=INB(JP),INF(JP)
      RNY=0.D0
      DO K=1,2
      RNY=RNY-(V(K,IY)-V(K,J))*XNW(K,IP)
      ENDDO
      RNY=RNY*DSW(IP)*C2
      DO K=1,2
      SXT(K,J)= RNY*PT(K,J)*PHI(J)
      SX(K,J)=0.D0
      ENDDO
      ENDDO

      DO I=INB(IP), INF(IP)
      DO J=INB(JP), INF(JP)
	D=0.D0
	WRX=0.D0
	WRY=0.D0
	SYTR=0.D0
	SXTR=0.D0
        J1= JMIN(I,JP)
        I1= JMIN(J,IP)
	DO K=1,2
        DX(K)=V(K,J)-V(K,I)
	D=D+DX(K)**2
	WRX=WRX+(P(K,J)-P(K,J1))*DX(K)
	WRY=WRY-(P(K,I)-P(K,I1))*DX(K)
	SYTR=SYTR+SYT(K,I)*DX(K)
	SXTR=SXTR+SXT(K,J)*DX(K)
	ENDDO
	D=1.D0/D**2
	WRX=WRX*D
	WRY=WRY*D
	SYTR=SYTR*D
	SXTR=SXTR*D
      DO K=1,2
      TEMP=SYTR*DX(K)
      SY(K,I)=SY(K,I)+WRX*DX(K)
      SX(K,J)=SX(K,J)+WRY*DX(K)
      AT(K,J)=AT(K,J)+TEMP
      AT(K,J1)=AT(K,J1)-TEMP
      TEMP=SXTR*DX(K)
      AT(K,I)=AT(K,I)-TEMP
      AT(K,I1)=AT(K,I1)+TEMP
      ENDDO
      ENDDO
      ENDDO

      DO I=INB(IP),INF(IP)
      RNX=0.D0
      DO K=1,2
      RNX=RNX+(V(K,IX)-V(K,I))*XNW(K,JP)
      ENDDO
      RNX=RNX*DSW(JP)*C2
      DO K=1,2
      A(K,I)=A(K,I)+ RNX*SY(K,I)
      ENDDO
      ENDDO

      DO J=INB(JP),INF(JP)
      RNY=0.D0
      DO K=1,2
      RNY=RNY-(V(K,IY)-V(K,J))*XNW(K,IP)
      ENDDO
      RNY=RNY*DSW(IP)*C2
      DO K=1,2
      A(K,J)=A(K,J)+ RNY*SX(K,J)
      ENDDO
      ENDDO

45    CONTINUE
44    CONTINUE

      DO IP=1,NP
      DO 330 JP=1,NP
      IF (JP.EQ.IP) GOTO 330
      DO K=1,2
      X1(K)=X0(K,JP)+L1(JP)*E1(K,JP)
      ENDDO
      DO I=INB(IP),INF(IP)
      CALL LINESTR (V(1,I),X0(1,JP),X1,XNW(1,JP),STR)
      J=JMIN(I,JP)
      DO K=1,2
      DO L=1,2
      A(K,I)=A(K,I)+C2*STR(K,L)*P(L,J)
      AT(L,J)=AT(L,J)+C2*STR(K,L)*PT(K,I)*PHI(I)
      ENDDO
      ENDDO
      ENDDO
330   CONTINUE
      ENDDO

      Call RBP (XCW,V,P,PHI,NVT,DM,SW,WPR)
      Call RBP (XCW,V,PT,PHI,NVT,DM,SW,WTPR)
      PRODT=0.D0
      PROD=0.D0
      DO I=1,NVT
      TEMPT=0.D0
      TEMP=0.D0
      DO K=1,2
      TEMP=TEMP+P(K,I)*XNA(K,I)
      TEMPT=TEMPT+PT(K,I)*XNA(K,I)
      ENDDO
      PROD=PROD+TEMP*PHI(I)
      PRODT=PRODT+TEMPT*PHI(I)
      ENDDO
      PROD=PROD/SW
      PRODT=PRODT/SW

      DO I=1,NVT
      DO K=1,2
      A(K,I)= P(K,I)-A(K,I) -WPR(K,I)+PROD*XNA(K,I)
      AT(K,I)= PT(K,I)-AT(K,I)/PHI(I)-WTPR(K,I)+PRODT*XNA(K,I)
      ENDDO
      ENDDO

      AL=0.D0
      DO I=1,NVT
      DO K=1,2
      AL=AL+PT(K,I)*A(K,I)*PHI(I)
      ENDDO
      ENDDO

      TE=0.D0
      DO I=1,NVT
      DO K=1,2
      TE=TE+P(K,I)*AT(K,I)*PHI(I)
      ENDDO
      ENDDO

      PRINT *, ' AL=', AL, ' TE=', TE
      AL= RR/AL
      DO I=1,NVT
      DO K=1,2
      R(K,I)= R(K,I) - AL*A(K,I)
      RT(K,I)= RT(K,I) - AL*AT(K,I)
      ENDDO
      ENDDO

      RRNEW=0.D0
      DO I=1,NVT
      DO K=1,2
      RRNEW=RRNEW+R(K,I)*RT(K,I)*PHI(I)
      ENDDO
      ENDDO

      BET=RRNEW/RR
      wtest=0.d0
      DO I=1,NVT
      DO K=1,2
      W(K,I)= W(K,I) - AL*P(K,I)
      P(K,I)= R(K,I) + BET*P(K,I)
      PT(K,I)= RT(K,I) + BET*PT(K,I)
      wtest=wtest+  dabs(w(k,i))*phi(i)
      ENDDO
      ENDDO
      print*, ' wtest=', wtest
      GOTO 777

778   continue
	CALL RBP(XCW,V,W,PHI,NVT,DM,SW,WPR)
	DO I=1,NVT
	 DO K=1,2
	  Q(K,I)=W(K,I)-WPR(K,I)/2.D0
	 ENDDO
	ENDDO
	OPEN (10, File = 'SOL2D.dat')
      WRITE (10,*)X0,L1,E1,V,Q,DSW,XNW,INB,INF,NVT
	CLOSE (10)
      STOP
      END

	SUBROUTINE RBP(XCW,V,W,PHI,NVT,DM,SW,WPR)
	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION XCW(2),V(2,*),W(2,*),PHI(*),WPR(2,*),A(2)
c	begin loop for A(K)
	DO K=1,2
	A(K)=0.D0
	DO I=1,NVT
	A(K)=A(K)+(1.0/SW)*(W(K,I)*PHI(I))
	ENDDO
	ENDDO
c	Begin loop for definition of omega
        OM=0.d0
	DO I=1,NVT
	OM=OM+((V(1,I)-XCW(1))*W(2,I)-(V(2,I)-XCW(2))*W(1,I) )*PHI(I)
	ENDDO
	OM=DM*OM
c	calculate for w, wpr
	DO I=1,NVT
	WPR(1,I)=A(1)-OM*(V(2,I)-XCW(2))
	WPR(2,I)=A(2)+OM*(V(1,I)-XCW(1))
	ENDDO
	RETURN
	END

C-----Subroutine PANEL added June 18 2007
	SUBROUTINE PANEL(I,NVTMAX,N1,L1,X0,E1,INB,INF,V,XNW,DSW)
c     change real to match original
	IMPLICIT REAL*8 (A-H,O-Z)
	REAL *8 L1
	DIMENSION X0(2),E1(2),V(2,*),XNW(2)
	H1=L1/N1
	INB=I + 1
	INF = I + N1
	DSW = H1
	XNW(1)=E1(2)
	XNW(2)=-E1(1)
	DO I1 = 1,N1
	I = I + 1
	IF (I.GT.NVTMAX) THEN
	PRINT *,'Out of Bounds - Check I GT NVTMAX'
	STOP
	ENDIF
	DO K=1,2
	V(K,I)=X0(K)+(I1-0.5D0)*H1*E1(K)
	ENDDO
	ENDDO
	RETURN
	END

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
      TA=DASIN(SI)
      SI=(E(1,1)*XB(2)-E(2,1)*XB(1))/DSQRT(BL)
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
