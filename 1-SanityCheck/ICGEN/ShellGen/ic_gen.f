      MODULE ic
      INTEGER (KIND=4) npbox,npshell
      REAL(KIND=4) ,dimension (:,:), ALLOCATABLE :: posshell,velcirc
      REAL(KIND=4) ,dimension (:,:), ALLOCATABLE :: poscube,velturb
      REAL(KIND=4) ,dimension (:,:), ALLOCATABLE :: posdef,veldef
      CHARACTER(20) :: box_file
      CHARACTER(20) :: ic_file
      !These next lines are the variables that are actually written
      !to the initial contitions file:
      REAL(KIND=4), dimension (:,:), ALLOCATABLE :: pos,vel
      INTEGER(KIND=4), dimension(:), ALLOCATABLE :: id
      REAL(KIND=4), dimension(:), ALLOCATABLE :: masses, u
      !These variables are for the header
      integer*4 npart_arr(0:5), nall_arr(0:5)
      real*8    massarr(0:5)
      real*8    a
      real*8    redshift
      integer*4 unused(34) 
      integer*4 flag_sfr
      integer*4 flag_feedbk
      END MODULE ic

      MODULE units
      REAL(KIND=4) pi,r_in,r_out,lbox,kpc,run,vrot,M_tot,tilt_deg
      END MODULE units

      MODULE coord_ini
      REAL(KIND=4) X0, Y0, Z0
      END MODULE coord_ini

      MODULE turb
      INTEGER (KIND=4) iseed1, iseed2
      REAL(KIND=4) vturb,nak,kmin,cnorm
      END MODULE turb

      MODULE ffta
      INTEGER (KIND=4) NM
      REAL(KIND=4) ,dimension (:,:,:), ALLOCATABLE :: avecx,avecy,avecz
      REAL(KIND=4) ,dimension (:,:,:), ALLOCATABLE :: vcurlx,vcurly,
     &vcurlz
      END MODULE ffta

      MODULE shell
      USE ic
      USE units
      USE coord_ini
      USE turb
      USE ffta
      END MODULE shell

      program shell_gen
      CALL readinput  ! read input parameters

      CALL readbox  ! read particles pos from a glass-like distr 

      CALL inparams  ! set up main parameters

      CALL findshell   ! find the pos of particles in the shell

      CALL findvcircular  ! find circ velocities

      CALL findvturb  ! find turbulent velocities

      CALL findtilt

      CALL writeic


      STOP
      END


      SUBROUTINE readbox
      USE shell
      IMPLICIT none
      real*4 smin(3),smax(3)
      INTEGER i,m,kchar,lenstr
C   reads from  given file a configuration of 10^6 particles 
C   distribution glass-like ;  box size L=1



      !open (1, file='BOX_1000k', form='unformatted')
      !open (1, file='BOX_8000k', form='unformatted')
      !open (1, file='BOX_8M', form='unformatted')
      kchar=lenstr(box_file,20)
      open (1, file=box_file(1:kchar), form='unformatted')
      !open (1, file='BOX_64M', form='unformatted')
      rewind 1
      read(1) npbox
      ALLOCATE(pos(1:3,1:npbox))
      read(1) ((pos(m,i),m=1,3),i=1,npbox)
      close(1)


      smin=1.e20
      smax=-1.e20

       DO i=1,npbox
         DO m=1,3
       smax(m)=MAX(smax(m),pos(m,i))
       smin(m)=MIN(smin(m),pos(m,i))
         END DO
       END DO

       write(6,*)' NUM P',npbox
       write(6,*)' MIN X',smin
       write(6,*)' MAX X',smax


      RETURN
      END

      SUBROUTINE readinput
      USE shell
      IMPLICIT none
C   read from   input file shell parameters : number of partciles npshell
C   inner radius r_in and outer radius r_out

      OPEN(2,FILE='INSHELL',FORM='FORMATTED')
      READ(2,*) box_file
      READ(2,*) npshell
      READ(2,*) r_in
      READ(2,*) r_out
      READ(2,*) vrot
      READ(2,*) vturb
      READ(2,*) tilt_deg
      READ(2,*) NM
      READ(2,*) M_tot
      READ(2,*) X0
      READ(2,*) Y0
      READ(2,*) Z0
      CLOSE(2)


      RETURN
      END


      SUBROUTINE inparams
      USE shell
      IMPLICIT none

      pi=4.*atan(1.)
      kpc=3.086e21
      run=kpc

      r_in=r_in*kpc
      r_out=r_out*kpc

      r_in=r_in/run
      r_out=r_out/run



      RETURN
      END


      SUBROUTINE  findshell
      USE shell
      IMPLICIT none
      real*4 volshell,sizec,avn,avcube,dxp(3),dist
      INTEGER npcube,incube,i,m,inshell
      LOGICAL testbox,testk(3)
C    setup size of box so to enclose the shell
      lbox=2.*r_out*(1.+0.02)  ! add a small space arund the shell

C    find how many particles to be extracted from the cube
      volshell=4.*pi*(r_out**3-r_in**3)/3.
      avn=float(npshell)/volshell
      npcube=avn*lbox**3   ! you need npcube parts in a cube of size lbox

      ALLOCATE(poscube(1:3,1:2*npcube))! avoid array memory violation

      IF(npcube.GT.npbox) THEN
      write(6,*)' npcube box',npcube,npbox   ! too many parts
      write(6,*)"you're requesting too much particles,"
      write(6,*)"Make a bigger box or request less particles"
      STOP
      END IF
      avcube=npbox
      sizec=(float(npcube)/avcube)**0.333333 ! extract npcube parts from the box

      incube=0
       DO i=1,npbox

       testbox=.TRUE.
       DO m=1,3
       testk(m)=.FALSE.
       IF(pos(m,i).LE.sizec) testk(m)=.TRUE.
       testbox=testbox.AND.testk(m)
       END DO

       IF(testbox) THEN
         incube=incube+1
         DO m=1,3
         poscube(m,incube)=pos(m,i)*(lbox/sizec)  ! rescale to the cube
         END DO
       END IF

       END DO

       write(6,*)' incube npcube',incube,npcube
       npcube=incube

       DEALLOCATE(pos)

       ALLOCATE(posshell(1:3,1:2*npshell))! avoid array memory violation

       inshell=0
       DO i=1,npcube

       testbox=.FALSE.
       dist=0
       DO m=1,3
       dxp(m)=poscube(m,i)-lbox/2   !   set to the cube center
       dist=dist+dxp(m)**2
       END DO

       IF(r_in**2.LE.dist.AND.dist.LE.r_out**2) testbox=.TRUE.

       IF(testbox) THEN
         inshell=inshell+1
         DO m=1,3
         posshell(m,inshell)=dxp(m)   ! origin at the cube center
         END DO
       END IF

       END DO

       write(6,*)' inchell npshell',inshell,npshell
       npshell=inshell

       ALLOCATE(velcirc(1:3,1:npshell)) ! circular velocities



      RETURN
      END


      SUBROUTINE  findvcircular  ! find circ velocities
      USE shell
      IMPLICIT none
      real*4 xi,yi,zi,vr,vt,vphi,rad,tht,phi,vx,vy,vz
      real*4 xxx
      INTEGER i


C     set-up vcirc
       velcirc=0
        vr=0
        vt=0
        vphi=vrot



        DO i=1,npshell

        xi=posshell(1,i)  ! range -lbox/2,lbox/2
        yi=posshell(2,i)
        zi=posshell(3,i)
        CALL POLAR(xi,yi,zi,rad,tht,phi,pi)
        vx=vr*sin(tht)*cos(phi)+vt*cos(tht)*cos(phi)-vphi*sin(phi)
        vy=vr*sin(tht)*sin(phi)+vt*cos(tht)*sin(phi)+vphi*cos(phi)
        vz=vr*cos(tht)-vt*sin(tht)
        velcirc(1,i)=vx
        velcirc(2,i)=vy
        velcirc(3,i)=vz
        xxx=0

        END DO


      RETURN
      END

      SUBROUTINE  findvturb  
      USE shell

      IF(vturb.EQ.0) THEN

       ALLOCATE(velturb(1:3,1:npshell)) ! turb velocities

       velturb=0

      ELSE

      CALL aknorm

      CALL setupak_2

      CALL setupvelturb_1

      END IF


      RETURN
      END

      SUBROUTINE  aknorm  ! find C
      USE shell
      IMPLICIT none
      real*4 sumint,smin,smax,sumi,volu,valint
      real*4 xxx
      INTEGER istep
      EXTERNAL fint


      nak=17./6
      cnorm=1
      kmin=2.*pi/lbox
      smin=kmin
      smax=10*smin
      volu=(lbox**3)

      sumint=0
      DO istep=1,4

      CALL qtrap(fint,smin,smax,sumi)

      smin=smax
      smax=10.*smin
      sumint=sumint+sumi
      write(6,*)' step int=',istep,sumi
      END DO

      valint=sumint*volu*4.*pi/(2.*pi)**3
C   sigma_v^2=vturb^2
      cnorm=3.*vturb**2/(2.*valint)

      RETURN
      END

      SUBROUTINE qtrap(func,a,b,s)
      INTEGER JMAX
      REAL a,b,func,s,EPS
      EXTERNAL func
      PARAMETER (EPS=1.e-4, JMAX=40)
CU    USES trapzd
      INTEGER j
      REAL olds
      IF(a.GE.b) THEN
      s=0
      RETURN
      END IF

      olds=-1.e30
      do 11 j=1,JMAX
        call trapzd(func,a,b,s,j)
        if (abs(s-olds).lt.EPS*abs(olds)) return
        olds=s
11    continue
      write(21,*)'too many steps in qtrap'
      STOP
      END

      SUBROUTINE findtilt
      USE shell
      IMPLICIT NONE
      integer*4 Nparticles, i, k
      real*4 x,y,z, vx, vy, vz
      real*4 xp,yp,zp, vxp, vyp, vzp, psi
      real*4 vxmax, vymax, vzmax, vxmin, vymin, vzmin
      real*4 xmax, ymax, zmax, xmin, ymin, zmin
      Nparticles = npshell

      ALLOCATE(posdef(1:3,1:Nparticles)) ! positions of the particles
      ALLOCATE(veldef(1:3,1:Nparticles)) ! total velocities
      posdef = 0
      veldef = 0 
      psi = tilt_deg * pi / 180.0
      write(6,*) 'tilting by', psi, 'rad'
      !Positions
      xmax = -10
      ymax = -10
      zmax = -10
      xmin = 10
      ymin = 10
      zmin = 10
      DO i=1,Nparticles
        x = posshell(1,i)
        y = posshell(2,i)
        z = posshell(3,i)
        xp =   x * cos(psi) + z * sin(psi)
        yp =   y
        zp = - x * sin(psi) + z * cos(psi)
        posdef(1,i) = xp + X0
        posdef(2,i) = yp + Y0
        posdef(3,i) = zp + Z0
        if(xp > xmax) xmax = xp
        if(yp > ymax) ymax = yp
        if(zp > zmax) zmax = zp
        if(xp < xmin) xmin = xp
        if(yp < ymin) ymin = yp
        if(zp < zmin) zmin = zp
      END DO
      write(6,*) 'x, y, z max:', xmax, ymax, zmax
      write(6,*) 'x, y, z min:', xmin, ymin, zmin

      !Give velocities
      vxmax = -10
      vymax = -10
      vzmax = -10
      vxmin = 10
      vymin = 10
      vzmin = 10
      DO i=1,Nparticles
        vx = velcirc(1,i) + velturb(1,i)
        vy = velcirc(2,i) + velturb(2,i)
        vz = velcirc(3,i) + velturb(3,i)
        vxp =   vx * cos(psi) + vz * sin(psi)
        vyp =   vy
        vzp = - vx * sin(psi) + vz * cos(psi)
        veldef(1,i) = vxp
        veldef(2,i) = vyp
        veldef(3,i) = vzp
        if(vxp > vxmax) vxmax = vxp
        if(vyp > vymax) vymax = vyp
        if(vzp > vzmax) vzmax = vzp
        if(vxp < vxmin) vxmin = vxp
        if(vyp < vymin) vymin = vyp
        if(vzp < vzmin) vzmin = vzp
      END DO
      write(6,*) 'vx, vy, vz max:', vxmax, vymax, vzmax
      write(6,*) 'vx, vy, vz min:', vxmin, vymin, vzmin
    
      RETURN
      END 

      SUBROUTINE  setupak_2  ! initialize A_k
      USE shell
      IMPLICIT none
      real*4 sumint,smin,smax,sumi,volu,valint
      real*4 xxx
      INTEGER istat

      ALLOCATE(avecx(0:NM+1,0:NM-1,0:NM-1),STAT=istat) ! allocate A comps
        IF(istat.GT.0) THEN
        write(6,*) ' alloc A too large'
        STOP
        END IF

      ALLOCATE(avecy(0:NM+1,0:NM-1,0:NM-1),STAT=istat) ! allocate A comps
        IF(istat.GT.0) THEN
        write(6,*) ' alloc A too large'
        STOP
        END IF

      ALLOCATE(avecz(0:NM+1,0:NM-1,0:NM-1),STAT=istat) ! allocate A comps
        IF(istat.GT.0) THEN
        write(6,*) ' alloc A too large'
        STOP
        END IF

        avecx=0
        avecy=0
        avecz=0

        iseed1=-2
        iseed2=-5

        CALL setupwn(avecx,NM) ! setup a Gaussian WN in x space
        CALL setupwn(avecy,NM) ! setup a Gaussian WN in x space
        CALL setupwn(avecz,NM) ! setup a Gaussian WN in x space

       CALL fftveci(avecx,NM)  ! real --> cmplx transforms
       CALL fftveci(avecy,NM)
       CALL fftveci(avecz,NM)

       CALL initak(avecx,avecy,avecz,NM)  ! setup A_K 


       CALL power(avecx,avecy,avecz,NM) ! Power spectra

       CALL fftvec(avecx,NM)  ! cmplx-->real transforms
       CALL fftvec(avecy,NM)
       CALL fftvec(avecz,NM)



      RETURN
      END

       SUBROUTINE setupwn(AR,NM)
       USE turb
       REAL AR(0:NM+1,0:NM-1,0:NM-1)
       XNM3=NM**3
       pi=4.*atan(1.)

        sumx=0

        DO k=0,nm-1
          DO j=0,nm-1
           DO i=0,nm-1
            csi1=ran1(iseed1)  ! rand number 0--1
            csi2=ran1(iseed2)  ! rand number 0--1
            YR=SQRT(-2.*ALOG(csi1)*XNM3)
            phase=2.*pi*ran1(iseed2)
            AR(i,j,k)=YR*COS(phase)
            sumx=sumx+AR(i,j,k)
           END DO
           END DO
           END DO
      
          sumx=sumx/XNM3

        DO k=0,nm-1
          DO j=0,nm-1
           DO i=0,nm-1
            AR(i,j,k)=AR(i,j,k)-sumx
           END DO
           END DO
           END DO

       AR=AR/XNM3


      RETURN
      END


       SUBROUTINE initak(avecx,avecy,avecz,NM)
       USE units
       REAL KMOD2,KMOD
      COMPLEX CF,AVECMP
      COMPLEX avecx(0:NM/2,0:NM-1,0:NM-1)
      COMPLEX avecy(0:NM/2,0:NM-1,0:NM-1)
      COMPLEX avecz(0:NM/2,0:NM-1,0:NM-1)


      NMH=NM/2
      NM1=NM-1
      NMH1=NMH+1
      pi=4.*atan(1.)
      KMIN=2.*pi/lbox

      DO 100 K=0,NM1
      KZ=K-(K/NMH1)*NM   ! kz=0...N/2,-N2/+1 ... -1

      DO 100 J=0,NM1
      JY=J-(J/NMH1)*NM

      DO 100 I=0,NM/2
      KMOD2=(I**2+JY**2+KZ**2)
      KMOD=KMIN*SQRT(KMOD2)
      PKSQROOT=VECPOT(KMOD)/SQRT(3.)  ! A_k

      avecx(i,j,k)=avecx(i,j,k)*PKSQROOT
      avecy(i,j,k)=avecy(i,j,k)*PKSQROOT
      avecz(i,j,k)=avecz(i,j,k)*PKSQROOT


100   CONTINUE

      RETURN
      END

       SUBROUTINE power(avecx,avecy,avecz,NM)
       USE units
       REAL KMOD2,KMOD,KMIN,KMAX
      COMPLEX CF,AVECMP(3)
      COMPLEX avecx(0:NM/2,0:NM-1,0:NM-1)
      COMPLEX avecy(0:NM/2,0:NM-1,0:NM-1)
      COMPLEX avecz(0:NM/2,0:NM-1,0:NM-1)
      REAL(KIND=4) ,dimension (:), ALLOCATABLE :: pbin,nbink


      NMH=NM/2
      NM1=NM-1
      NMH1=NMH+1
      pi=4.*atan(1.)
      KMIN=2.*pi/lbox
      KMAX=KMIN*NMH
      NBIN=NMH
      NBIN=NM
      ALLOCATE(PBIN(0:NBIN+1),NBINK(0:NBIN+1))



      PBIN=0
      NBINK=0

      DBIN=(KMAX-KMIN)/NBIN

      DO 100 K=0,NM1
      KZ=K-(K/NMH1)*NM   ! kz=0...N/2,-N2/+1 ... -1

      DO 100 J=0,NM1
      JY=J-(J/NMH1)*NM

      DO 100 I=0,NM/2
      KMOD2=(I**2+JY**2+KZ**2)
      KMOD=KMIN*SQRT(KMOD2)
      IF(KMOD.GE.KMAX) GO TO 100
      IF(KMOD.LT.KMIN) GO TO 100

      avecmp(1)=avecx(i,j,k)
      avecmp(2)=avecy(i,j,k)
      avecmp(3)=avecz(i,j,k)
      PKSQ=0
      DO ICOMP=1,3
      PKSQ=PKSQ+CABS(AVECMP(ICOMP))**2
      END DO
      wbin=(KMOD-KMIN)/DBIN
      jbin=wbin
      tw=wbin-jbin
      PBIN(jbin)=PBIN(jbin)+PKSQ*(1.-tw)
      PBIN(jbin+1)=PBIN(jbin+1)+PKSQ*tw
      NBINK(jbin)=NBINK(jbin)+(1.-tw)
      NBINK(jbin+1)=NBINK(jbin+1)+tw


100   CONTINUE


      DO jbin=0,NBIN
      IF(NBINK(JBIN).GT.0) THEN
      PBIN(jbin)=PBIN(jbin)/NBINK(jbin)
      END IF
      END DO

      DO jbin=0,NBIN
      xk=kmin+jbin*dbin
      ak2=vecpot(xk)**2
      write(15,900) jbin,xk,PBIN(jbin),ak2
      END DO

      xxx=0
      DEALLOCATE(PBIN,NBINK)
900   FORMAT(1x,i3,2x,f8.3,2(1x,e10.4))

      RETURN
      END


      SUBROUTINE  setupvelturb_1  ! find v turb for particles
      USE shell
      IMPLICIT none
      real*4 xi,yj,zk,hmesh,ce1,ce2,ce3,zd(3),z1,zh,zw(-1:1,3),xxx
      real*4 ax,ay,az,axx,ayy,azz,ztemp1,ztemp2,ztemp3,daxy,daxz,
     &dayz,dazx,dazy,dayx
      INTEGER i,j,k,lx,ly,lz,IX(-4:4),IY(-4:4),IZ(-4:4),MM,im,jm1,jm2,
     &kp1,km1,kp2,km2,kp3,km3,jp1,jp2,jp3,jm3,im1,im2,ip1,ip2,ip3,im3,
     &jm,km,nm1,ip

      ALLOCATE(vcurlx(0:NM-1,0:NM-1,0:NM-1)) 
      ALLOCATE(vcurly(0:NM-1,0:NM-1,0:NM-1)) 
      ALLOCATE(vcurlz(0:NM-1,0:NM-1,0:NM-1)) 

       hmesh=lbox/nm
*      6 POINTS GRADIENT
         CE1=9./(12.*HMESH)
         CE2=-3./(20.*HMESH)
         CE3=1./(60.*HMESH)
        NM1=NM-1


       DO 100 K=0,NM1
       DO 100 J=0,NM1
       DO 100 I=0,NM1

       IX(-4)=I-4
       IY(-4)=J-4
       IZ(-4)=K-4
       IF(IX(-4).LT.0) IX(-4)=IX(-4)+NM
       IF(IY(-4).LT.0) IY(-4)=IY(-4)+NM
       IF(IZ(-4).LT.0) IZ(-4)=IZ(-4)+NM

       DO MM=-3,4
       IX(MM)=MOD(IX(MM-1)+1,NM)
       IY(MM)=MOD(IY(MM-1)+1,NM)
       IZ(MM)=MOD(IZ(MM-1)+1,NM)
       END DO

        KM=0
        KP1=IZ(KM+1)                                                    
        KM1=IZ(KM-1)                                                    
        KP2=IZ(KM+2)                                                    
        KM2=IZ(KM-2)                                                  
        KP3=IZ(KM+3)
        KM3=IZ(KM-3)
        JM=0
        JP1=IY(JM+1)                                                    
        JM1=IY(JM-1)                                                    
        JP2=IY(JM+2)                                                   
        JM2=IY(JM-2)                                                   
        JP3=IY(JM+3)
        JM3=IY(JM-3)                                                    
        IM=0                                                       
        IP1=IX(IM+1)                                                    
        IM1=IX(IM-1)                                                    
        IP2=IX(IM+2)                                                    
        IM2=IX(IM-2)                                                    
        IP3=IX(IM+3)
        IM3=IX(IM-3)                                                    
C   v_x=A_z,y-A_y,z   v_y=A_x,z-A_z,x  v_z=Ay,x-A_x,y

      DAXY=CE1*(AVECX(I,JP1,K)-AVECX(I,JM1,K))+CE2*(AVECX(I,JP2,K)-
     &AVECX(I,JM2,K))+CE3*(AVECX(I,JP3,K)-AVECX(I,JM3,K))

      DAZY=CE1*(AVECZ(I,JP1,K)-AVECZ(I,JM1,K))+CE2*(AVECZ(I,JP2,K)-
     &AVECZ(I,JM2,K))+CE3*(AVECZ(I,JP3,K)-AVECZ(I,JM3,K))

      DAYZ=CE1*(AVECY(I,J,KP1)-AVECY(I,J,KM1))+CE2*(AVECY(I,J,KP2)-
     &AVECY(I,J,KM2))+CE3*(AVECY(I,J,KP3)-AVECY(I,J,KM3))

      DAXZ=CE1*(AVECX(I,J,KP1)-AVECX(I,J,KM1))+CE2*(AVECX(I,J,KP2)-
     &AVECX(I,J,KM2))+CE3*(AVECX(I,J,KP3)-AVECX(I,J,KM3))

      DAYX=CE1*(AVECY(IP1,J,K)-AVECY(IM1,J,K))+CE2*(AVECY(IP2,J,K)-
     &AVECY(IM2,J,K))+CE3*(AVECY(IP3,J,K)-AVECY(IM3,J,K))

      DAZX=CE1*(AVECZ(IP1,J,K)-AVECZ(IM1,J,K))+CE2*(AVECZ(IP2,J,K)-
     &AVECZ(IM2,J,K))+CE3*(AVECZ(IP3,J,K)-AVECZ(IM3,J,K))


        vcurlx(i,j,k)=dazy-dayz
        vcurly(i,j,k)=daxz-dazx
        vcurlz(i,j,k)=dayx-daxy


100     CONTINUE


       DEALLOCATE(avecx,avecy,avecz)

       ALLOCATE(velturb(1:3,1:npshell)) ! circular velocities

       velturb=0
       Z1=1.
       ZH=0.5


        DO ip=1,npshell

        xi=(posshell(1,ip)+lbox/2)/hmesh  ! range 0<= x <  lbox
        yj=(posshell(2,ip)+lbox/2)/hmesh
        zk=(posshell(3,ip)+lbox/2)/hmesh


        LX=XI+0.5
        LY=YJ+0.5
        LZ=ZK+0.5

        ZD(1)=XI+0.5-LX
        ZD(2)=YJ+0.5-LY
        ZD(3)=ZK+0.5-LZ

      DO 10 J=1,3                                                       
      ZW(1,J)=ZH*ZD(J)*ZD(J)                                            
      ZW(-1,J)=ZH-ZD(J)+ZW(1,J)                                         
      ZW(0,J)=Z1-ZW(-1,J)-ZW(1,J)                                       
10    CONTINUE                                                          

       IX(-1)=LX-1
       IY(-1)=LY-1
       IZ(-1)=LZ-1
       IF(IX(-1).LT.0) IX(-1)=IX(-1)+NM
       IF(IY(-1).LT.0) IY(-1)=IY(-1)+NM
       IF(IZ(-1).LT.0) IZ(-1)=IZ(-1)+NM
       DO MM=0,1
       IX(MM)=MOD(IX(MM-1)+1,NM)
       IY(MM)=MOD(IY(MM-1)+1,NM)
       IZ(MM)=MOD(IZ(MM-1)+1,NM)
       END DO


      AX=0                                                              
      AY=0                                                              
      AZ=0                                                              
                                                                        
      DO 20 KM=-1,1                                                     
        ZTEMP3=ZW(KM,3)                                                 
        K=IZ(KM)                                                        
      DO 20 JM=-1,1                                                     
        ZTEMP2=ZTEMP3*ZW(JM,2)                                          
        J=IY(JM)                                                        

      DO 20 IM=-1,1                                                     
      ZTEMP1=ZW(IM,1)*ZTEMP2                                            
                                                                        
      I=IX(IM)                                                        

      AXX=vcurlx(i,j,k)
      AYY=vcurly(i,j,k)
      AZZ=vcurlz(i,j,k)

      AX=AX+ZTEMP1*AXX
      AY=AY+ZTEMP1*AYY
      AZ=AZ+ZTEMP1*AZZ
20    CONTINUE                                                         

      velturb(1,ip)=ax
      velturb(2,ip)=ay
      velturb(3,ip)=az

      END DO

      xxx=0

      RETURN
      END



      SUBROUTINE  setupvelturb_2  ! find v turb for particles
      USE shell
      IMPLICIT none
      real*4 xi,yj,zk,hmesh,ce1,ce2,ce3,zd(3),z1,zh,zw(-1:1,3),xxx
      real*4 ax,ay,az,axx,ayy,azz,ztemp1,ztemp2,ztemp3,daxy,daxz,
     &dayz,dazx,dazy,dayx,vrotx,vroty,vrotz
      INTEGER i,j,k,lx,ly,lz,IX(-4:4),IY(-4:4),IZ(-4:4),MM,im,jm1,jm2,
     &kp1,km1,kp2,km2,kp3,km3,jp1,jp2,jp3,jm3,im1,im2,ip1,ip2,ip3,im3,
     &jm,km,nm1,ip


       hmesh=lbox/nm
*      6 POINTS GRADIENT
         CE1=9./(12.*HMESH)
         CE2=-3./(20.*HMESH)
         CE3=1./(60.*HMESH)
        NM1=NM-1



       ALLOCATE(velturb(1:3,1:npshell)) ! circular velocities

       velturb=0
       Z1=1.
       ZH=0.5


        DO ip=1,npshell

        xi=(posshell(1,ip)+lbox/2)/hmesh  ! range 0<= x <  lbox
        yj=(posshell(2,ip)+lbox/2)/hmesh
        zk=(posshell(3,ip)+lbox/2)/hmesh


        LX=XI+0.5
        LY=YJ+0.5
        LZ=ZK+0.5

        ZD(1)=XI+0.5-LX
        ZD(2)=YJ+0.5-LY
        ZD(3)=ZK+0.5-LZ

      DO 10 J=1,3                                                       
      ZW(1,J)=ZH*ZD(J)*ZD(J)                                            
      ZW(-1,J)=ZH-ZD(J)+ZW(1,J)                                         
      ZW(0,J)=Z1-ZW(-1,J)-ZW(1,J)                                       
10    CONTINUE                                                          

       IX(-4)=LX-4
       IY(-4)=LY-4
       IZ(-4)=LZ-4
       IF(IX(-4).LT.0) IX(-4)=IX(-4)+NM
       IF(IY(-4).LT.0) IY(-4)=IY(-4)+NM
       IF(IZ(-4).LT.0) IZ(-4)=IZ(-4)+NM
       DO MM=-3,4
       IX(MM)=MOD(IX(MM-1)+1,NM)
       IY(MM)=MOD(IY(MM-1)+1,NM)
       IZ(MM)=MOD(IZ(MM-1)+1,NM)
       END DO

      AX=0                                                              
      AY=0                                                              
      AZ=0                                                              

      DO 20 KM=-1,1                                                     
        ZTEMP3=ZW(KM,3)                                                 
        K=IZ(KM)
        KP1=IZ(KM+1)                                                    
        KM1=IZ(KM-1)                                                    
        KP2=IZ(KM+2)                                                    
        KM2=IZ(KM-2)                                                  
        KP3=IZ(KM+3)
        KM3=IZ(KM-3)
      DO 20 JM=-1,1                                                     
        ZTEMP2=ZTEMP3*ZW(JM,2)                                          

        J=IY(JM)
        JP1=IY(JM+1)                                                    
        JM1=IY(JM-1)                                                    
        JP2=IY(JM+2)                                                   
        JM2=IY(JM-2)                                                   
        JP3=IY(JM+3)
        JM3=IY(JM-3)                                                    

        DO 20 IM=-1,1                                                
        ZTEMP1=ZW(IM,1)*ZTEMP2                                         
        I=IX(IM)                                                       
        IP1=IX(IM+1)                                                    
        IM1=IX(IM-1)                                                    
        IP2=IX(IM+2)                                                    
        IM2=IX(IM-2)                                                    
        IP3=IX(IM+3)
        IM3=IX(IM-3)                                                    
C   v_x=A_z,y-A_y,z   v_y=A_x,z-A_z,x  v_z=Ay,x-A_x,y

      DAXY=CE1*(AVECX(I,JP1,K)-AVECX(I,JM1,K))+CE2*(AVECX(I,JP2,K)-
     &AVECX(I,JM2,K))+CE3*(AVECX(I,JP3,K)-AVECX(I,JM3,K))

      DAZY=CE1*(AVECZ(I,JP1,K)-AVECZ(I,JM1,K))+CE2*(AVECZ(I,JP2,K)-
     &AVECZ(I,JM2,K))+CE3*(AVECZ(I,JP3,K)-AVECZ(I,JM3,K))

      DAYZ=CE1*(AVECY(I,J,KP1)-AVECY(I,J,KM1))+CE2*(AVECY(I,J,KP2)-
     &AVECY(I,J,KM2))+CE3*(AVECY(I,J,KP3)-AVECY(I,J,KM3))

      DAXZ=CE1*(AVECX(I,J,KP1)-AVECX(I,J,KM1))+CE2*(AVECX(I,J,KP2)-
     &AVECX(I,J,KM2))+CE3*(AVECX(I,J,KP3)-AVECX(I,J,KM3))

      DAYX=CE1*(AVECY(IP1,J,K)-AVECY(IM1,J,K))+CE2*(AVECY(IP2,J,K)-
     &AVECY(IM2,J,K))+CE3*(AVECY(IP3,J,K)-AVECY(IM3,J,K))

      DAZX=CE1*(AVECZ(IP1,J,K)-AVECZ(IM1,J,K))+CE2*(AVECZ(IP2,J,K)-
     &AVECZ(IM2,J,K))+CE3*(AVECZ(IP3,J,K)-AVECZ(IM3,J,K))


        vrotx=dazy-dayz
        vroty=daxz-dazx
        vrotz=dayx-daxy

      AXX=vrotx
      AYY=vroty
      AZZ=vrotz

      AX=AX+ZTEMP1*AXX
      AY=AY+ZTEMP1*AYY
      AZ=AZ+ZTEMP1*AZZ
20    CONTINUE                                                         

      velturb(1,ip)=ax
      velturb(2,ip)=ay
      velturb(3,ip)=az

      END DO

      xxx=0
       DEALLOCATE(avecx,avecy,avecz)

      RETURN
      END


       SUBROUTINE fftvec(AR,NM)
       USE units
       include 'fftw3.f'
       REAL AR(1:NM+2,NM,NM)
       integer*8 plan

       NL=NM
       NN=NM

       call dfftw_plan_dft_c2r_3d(plan,NL,NM,NN,ar,ar,FFTW_ESTIMATE)
       call fftw_execute(plan)
       call fftw_destroy_plan(plan)

        RETURN
        END 


       SUBROUTINE fftveci(AR,NM)
       USE units
       include 'fftw3.f'
       REAL AR(1:NM+2,NM,NM)
       integer*8 plan

       NL=NM
       NN=NM

       call fftw_plan_dft_r2c_3d(plan,NL,NM,NN,ar,ar,FFTW_ESTIMATE)
       call fftw_execute(plan)
       call fftw_destroy_plan(plan)

        RETURN
        END 



      REAL FUNCTION vecpot(X)
      USE turb
      AK2=cnorm/(x**2+kmin**2)**nak   ! A_k**2
      vecpot=SQRT(AK2)
      RETURN
      END

      REAL FUNCTION fint(X)
      USE turb
      fint=(vecpot(x)*(x**2))**2  ! A_k^2 k^4
      RETURN
      END

      SUBROUTINE writeic
      USE shell
      IMPLICIT NONE
      integer*4 Nparticles, i, k
      real*8 m_sph
      Nparticles = npshell
      m_sph = M_tot/real(Nparticles)    

      ALLOCATE(pos(1:3,1:Nparticles)) ! positions of the particles
      ALLOCATE(vel(1:3,1:Nparticles)) ! total velocities
      ALLOCATE(id(1:Nparticles)) !Particles id's
      !ALLOCATE(masses blablabla)
      ALLOCATE(u(1:Nparticles))
      pos = 0
      vel = 0
      id = 0
      u = 0

      !Header
       npart_arr = (/ Nparticles,0,0,0,0,0/)
       nall_arr = (/Nparticles,0,0,0,0,0/)
       massarr = (/ m_sph,0d0,0d0,0d0,0d0,0d0/)
       a = 0.0
       redshift = 0.0
       flag_sfr = 0 
       flag_feedbk = 0
      
      !Positions
      DO i=1,Nparticles
       DO k=1,3
        pos(k,i)=posdef(k,i)
       END DO
      END DO

      !Give velocities
      DO i=1,Nparticles
       DO k=1,3
        vel(k,i)=veldef(k,i)
       END DO
      END DO

      !Give ids
      DO i = 1, Nparticles
       id(i) = i
      END DO

      !Give energies
      DO i = 1, Nparticles
       u(i) = 0.
      END DO
      ic_file = 'snapshot_shell'
      OPEN(12,file=ic_file,form='unformatted', convert='little_endian')
!      OPEN(13,file='ic_snap.txt',form='formatted')
      write(6,*) '*****************'
      write(6,*) 'Giving initial conditions...'
      write(6,*) 'Npart=',Nparticles
      write(6,*) 'v_rot=',vrot
      write(6,*) 'v_turb=',vturb
      write(6,*) 'tilt=',tilt_deg,'degrees'
      write(6,*) '(X0,Y0,Z0)=',X0,Y0,Z0
      write(6,*) 'IC snapshot will be written to:',ic_file
      write(6,*) '*****************'

      write(12) npart_arr, massarr , a, redshift, 
     &flag_sfr,flag_feedbk, nall_arr, unused
      write(12) ((pos(k,i),k=1,3),i=1,Nparticles)
!      do i = 1, Nparticles 
!       write(13,*) posshell(1,i), posshell(2,i), posshell(3,i) 
!      end do
      write(12) ((vel(k,i),k=1,3),i=1,Nparticles)  ! total velocities
      write(12) (id(i),i=1,Nparticles)
      write(12) (u(i),i=1,Nparticles)
      CLOSE(12) 
!      CLOSE(13) 

      RETURN
      END

      SUBROUTINE POLAR(XT,YT,ZT,R,THETA,PHI,PGR)
      R=SQRT(XT**2+YT**2+ZT**2)
      ARG=ABS(ZT)/R
      IF(ARG.GT.1) ARG=1
      IF(ARG.LT.-1) ARG=-1
      THETA=ACOS(ARG)
      IF(ZT.LT.0) THETA=PGR-THETA
      IF(XT.NE.0) GO TO 10
      IF(XT.EQ.0.AND.YT.EQ.0) PHI=0
      IF(XT.EQ.0.AND.YT.GT.0) PHI=PGR/2
      IF(XT.EQ.0.AND.YT.LT.0) PHI=1.5*PGR
      GO TO 30
10    IF(YT.NE.0) GO TO 20
      IF(XT.GT.0.AND.YT.EQ.0) PHI=0
      IF(XT.LT.0.AND.YT.EQ.0) PHI=PGR
      GO TO 30
20    PHI=ATAN(ABS(YT)/ABS(XT))
      IF(XT.LT.0.AND.YT.GT.0) PHI=PGR-PHI
      IF(XT.LT.0.AND.YT.LT.0) PHI=PGR+PHI
      IF(XT.GT.0.AND.YT.LT.0) PHI=2.*PGR-PHI
30    CONTINUE
      RETURN
      END



      SUBROUTINE trapzd(func,a,b,s,n)
      INTEGER n
      REAL a,b,s,func
      EXTERNAL func
      INTEGER it,j
      REAL del,sum,tnm,x
      if (n.eq.1) then
        s=0.5*(b-a)*(func(a)+func(b))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+del
11      continue
        s=0.5*(s+(b-a)*sum/tnm)
      endif
      return
      END


      FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END


*----------------------------------------------------------------000708-
      INTEGER FUNCTION lenstr(name,len)
*-----------------------------------------------------------------------
      CHARACTER*(*) name
      integer       len

      j = len

      do while (name(j:j) .eq. ' ')
         j = j - 1
      enddo       
      
      lenstr = j        

      if ( j .eq. 1 ) then
         write(6,*)'Warning : String begins with empty space' 
      endif

      end

