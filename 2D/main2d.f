      Implicit Real*8 (a-h,o-z)
      Parameter (NP=14, NVTMAX=100000)
      COMMON / BIVEL / X0,L1,E1,V,Q,DSW,XNW,INB,INF,NVT

      Dimension V(2,NVTMAX),X0(2,NP),L1(NP),Q(2,NVTMAX),
     >E1(2,NP),INB(NP),INF(NP),DSW(NP),XNW(2,NP),
     >Y1(2),U1(2)
      Real*8 L1

      Open(10,file='SOL2D.dat')
      Read(10,*) X0,L1,E1,V,Q,DSW,XNW,INB,INF,NVT
      Close(10)

      Y1(1)=-4.d0
      Y1(2)=0.5d0

      dt=1.D-2
      dt=3.d-3

      T=0.D0

      Open(11,file='streamline.dat')

      Do 1 IJK=1,150000

      Call VELOCITY(Y1,U1)


      Do I=1,2
      Y1(I)=Y1(I)+DT*U1(I)
      endDo

      if (Y1(1).GT.18.d0) STOP
    
      T=T+DT
      
      if (mod(ijk,100).eq.0) then
      print*, ijk, T, Y
      END IF
      Write(11,*) Y1(1),Y1(2)

1     Continue


      Stop
      End




