      ! BallPark code computing ballistic trajectories and formation of scoria cones
      
      MODULE mConst
      IMPLICIT NONE
      
      INTEGER, PARAMETER :: Nd=1, Nv=1000, Na=10000, Nb=Nd*Nv*Na, Ng=6000, nvar=4
      REAL, PARAMETER :: pi=3.141592653589793
      ! solvers
      INTEGER, SAVE :: ilib=1, MAXSTP=10000, OPTSTP=100
      REAL, SAVE :: eps=1.e-8, TINY=1.e-20
      ! switches
      LOGICAL, SAVE :: plotball=.false., build_cone=.true., load_cone=.false., showdist=.false., airlessan=.false.
      LOGICAL, SAVE :: avalanche=.false., backlanche=.false., checkvolume=.false., EPSL14=.false., rndsmooth=.true.
      LOGICAL, SAVE :: flatan=.false.
      ! ranges
      REAL, SAVE :: angmin=0., angmax=45., dmax=0.3, dmin=0., vmin=10., vmax=300., tmax=1000., basemax=15.e3
      REAL, SAVE :: dmean=0.04, ddev=0.3, vmean=46., vdev=0.2, angdev=31. 
      ! params
      REAL, SAVE :: grav=3.71, rhoair=0.01, rhorock=850., Cd=0.7, conetop=550., incrvol=1.e8, maxvol=1.e11
      CHARACTER(100), SAVE :: trnam='traj', cnam='cone', lconenam='cone.dat'
      
      namelist /switches/ plotball, build_cone, load_cone, lconenam, showdist, EPSL14, &
        avalanche, backlanche, checkvolume, airlessan, rndsmooth, flatan
      namelist /ranges/ angmax, angmin, dmax, dmin, vmax, vmin, tmax, basemax, &
        dmean, ddev, vmean, vdev, angdev
      namelist /params/ grav, rhoair, rhorock, Cd, conetop, incrvol, maxvol
      namelist /solvers/ ilib, MAXSTP, OPTSTP, eps
                  
      END MODULE
          
      PROGRAM Ballpark
      USE mConst                ! Shared constants
      USE nr                    ! Contains ludcmp/ludskb, ODE, and other solvers
      
      IMPLICIT NONE
      
      INTEGER :: ib,k1,k2,k3,index(3),ixhit,dist(Nb),ilay=0,landed(Nb),overshot(Nb)
      REAL :: time, dt, dx, hdid, hnext, hmax, dp, vp, ap, xscal(nvar), t3, t2, t1
      REAL :: xvec(nvar), dxvecdt(nvar), drange(Nd), vrange(Nv), angrange(Na), dhmax, dhava, voll
      REAL :: conex(Ng), coneh(Ng), layer(Ng), weight(Nb), rnd(Nb), wv, wd, wang, vollay(Nb), checkvol, danal
      REAL :: volume=0., tball=0., tcone=0.
      REAL,ALLOCATABLE :: traj(:,:,:),dcoll(:,:),vcoll(:,:),angcoll(:,:)
      LOGICAL :: hotovo=.false.
            
      open(1,file='param.in',status='old')
      read(1,switches); rewind(1);
      read(1,params); rewind(1);
      read(1,ranges); rewind(1);
      read(1,solvers); rewind(1);
      close(1)
      call LIve()

      IF(plotball) allocate(traj(5,MAXSTP,Nb))
      IF (load_cone) THEN
        call loadcone(conex,coneh)
        dx = conex(2)-conex(1)
      ELSE
        coneh(1:Ng) = 0.
        dx = basemax/Ng
        conex(1:Ng) = (/((real(k1)+0.5)*dx, k1=0,Ng-1)/)
      ENDIF
      drange(1:Nd) = (/(dmin + (real(k1)+0.5)*(dmax-dmin)/Nd, k1=0,Nd-1)/)
      vrange(1:Nv) = (/(vmin + (real(k1)+0.5)*(vmax-vmin)/Nv, k1=0,Nv-1)/)
      angrange(1:Na) = (/(angmin + (real(k1)+0.5)*(angmax-angmin)/Na, k1=0,Na-1)/)      
      IF (showdist) THEN
        allocate(dcoll(Nd,2), vcoll(Nv,2), angcoll(Na,2))
        dcoll(:,2) = 0.; vcoll(:,2) = 0.; angcoll(:,2) = 0.
        dcoll(:,1) = drange(:); vcoll(:,1) = vrange(:); angcoll(:,1) = angrange(:)
      ENDIF
      dhmax = tan(30.*pi/180.)*dx
      print *,'Parameters of the particle distribution:'
      print '(a,f5.0,a,f10.5,a,f6.1)','  vmean: ', vmean,' [m/s], vdev: ', vdev,', ang. dev: ', angdev
      print '(a,f4.0,a,f10.0,a,f8.3)',' Initial cone height ', coneh(1),' with base ', basemax, ' and dx ', dx
      print '(a,f4.0,a,f4.1,a,e20.10)',' Friction over gravity (v=',vmean,' [m/s] d=',100*dmean,' [cm]): ', &
       (3.*Cd*rhoair)/(4.*dmean*rhorock)*vmean**2 / grav
      IF (airlessan.and.flatan) THEN
        print *,'VERSION: analytical solutions (ballistics), flat base for all layers'
      ELSE IF (airlessan) THEN
        print *,'VERSION: analytical solutions (ballistics), cone gradually constructed'
      ELSE
        print *,'VERSION: numerical solutions (ballistics), cone gradually constructed'
      ENDIF
      print *
      
      DO WHILE (.not.hotovo)
          landed(:) = 0
          overshot(:) = 0
          ! Unlanded particles are reported at the outer edge of the cone
          dist(:) = Ng
          weight(:) = 0.
          call cpu_time(t1)
          IF (rndsmooth) THEN
              CALL RANDOM_NUMBER(rnd)
              rnd(:) = (rnd(:) - 0.5) * (vmax-vmin)/Nv
          ELSE
              rnd(:) = 0.
          ENDIF

            
          ! Computing ballistics
!$OMP PARALLEL DO PRIVATE(time,dt,hnext,hmax,hdid,dp,vp,ap,xvec,xscal,dxvecdt,index,ixhit,wang,wd,wv)
          DO ib=1,Nb
            index = onetothree(ib)
            dp = drange(index(1))
            vp = vrange(index(2)) + rnd(ib)
            ap = angrange(index(3))
            wang = -(ap/angdev)**2 / 2.
            IF (EPSL14) THEN
                wv = -( (log10(vp/vmean)) / vdev )**2 / 2.
                wd = -( (log10(dp/dmean)) / ddev )**2 / 2.
                weight(ib) = sin(ap*pi/180.)*exp(wang + wv + wd) * (dp**3)
            ELSE
                wv = -( (log(vp/vmean)) / vdev )**2 / 2.
                wd = -( (log(dp/dmean)) / ddev )**2 / 2.
                weight(ib) = sin(ap*pi/180.)*exp(wang + wv + wd)/(vp*dp) * (dp**3)
            ENDIF

            xvec(1:2) = (/0., coneh(1)/)
            xvec(3:4) = vp * (/ sin(ap*pi/180.), cos(ap*pi/180.) /)
            time = 0.
            ! Based on the airless analytical solution (Eq. S5 in GRL18), allowed #steps: <0.5*OPTSTP, MAXSTP>
            hnext = (2.*vp*cos(ap*pi/180.)/grav) / (10.*OPTSTP)
            hmax = 20.*hnext
            IF(plotball) print *,'Ballistic',ib,'. Diameter, velocity, angle: ', dp, vp, ap
            IF(plotball) traj(:,1,ib) = (/ time, xvec(1), xvec(2), xvec(3), xvec(4) /)

            IF (airlessan.and.flatan) THEN
                danal = (vp**2) * sin(2.*ap*pi/180.) / grav
                ixhit = int(abs(danal)/dx) + 1
                dist(ib) = ixhit
                landed(ib) = 1
            ELSE IF (airlessan) THEN
                DO k1=1,MAXSTP
                    dt = 0.5*hmax
                    time = time + dt
                    xvec(1) = xvec(3)*time
                    xvec(2) = coneh(1) + xvec(4)*time - 0.5*grav*(time**2)
                    IF(plotball) traj(:,k1+1,ib) = (/time, xvec(1), xvec(2), xvec(3), xvec(4) - grav*time/)
                    ixhit = int(abs(xvec(1))/dx) + 1
                    IF(ixhit > Ng) exit
                    IF(xvec(2) <= coneh(ixhit)) THEN
                        ! Interpolating the impact point exactly to the cone level
                        xvec(1) = xvec(1) + ((coneh(ixhit) - xvec(2)) / (xvec(4) - grav*time)) * xvec(3)
                        ixhit = int(abs(xvec(1))/dx) + 1
                        dist(ib) = ixhit
                        vollay(ib) = vollay(ib) + 2.*pi*dx*(dp**3)
                        landed(ib) = 1
                        exit
                    ENDIF
                ENDDO            
            ELSE
                DO k1=1,MAXSTP
                    dt = hnext
                    call derivs(time,xvec,dxvecdt,dp)  
                
                    ! setting a measure for computing the accuracy of the obtained ballistic curve
                    xscal(:) = abs(xvec(:)) + abs(dt*dxvecdt(:)) + TINY
                            
                    ! ADVANCING THE SOLUTION in time for each ballistic curve
                    SELECT CASE(ilib)
                    CASE(0)
                        ! explicit Euler  
                        xvec = xvec + dxvecdt*dt
                        time = time + dt
                    CASE(1)
                        ! stepper rkqs + rkck   
                        call rkqs(xvec,dxvecdt,nvar,time,dt,eps,xscal,hdid,hnext,derivs,dp)
                        IF(hnext>hmax) hnext = hmax                        
                    END SELECT
                    IF(plotball) traj(:,k1+1,ib) = (/time, xvec(1), xvec(2), xvec(3), xvec(4)/)
                    ixhit = int(abs(xvec(1))/dx) + 1
                    ! Particle has overshot basemax
                    IF (ixhit > Ng) THEN
                        overshot(ib) = 1
                        exit
                    ENDIF
                    ! Particle has landed
                    IF(xvec(2) <= coneh(ixhit)) THEN
                        dist(ib) = ixhit
                        landed(ib) = 1
                        exit
                    ENDIF
                ENDDO       ! Time-stepping loop
            ENDIF    
          ENDDO         ! Loop over different ballistic curves
!$OMP END PARALLEL DO
          call cpu_time(t2)
          tball = tball + (t2-t1)
          IF(real(Nb-sum(landed))/Nb > 0.1) stop 'More than 10% particles did not land (overshot or not fell)'
          IF(real(sum(overshot))/Nb > 0.1) stop 'More than 10% particles overshot basemax'
          
          IF (build_cone) THEN
            ilay = ilay+1
            layer(:) = 0.
            vollay(:) = 0.
            DO ib=1,Nb
                ! Constants not important, the volume of one batch is arbitrary (only the distribution with distance matters)
                layer(dist(ib)) = layer(dist(ib)) + weight(ib) / conex(dist(ib))
                ! here even the constants must be accounted, it is the actual volume added as the height was increased                                 
                vollay(ib) = vollay(ib) + 2.*pi*dx*weight(ib)
                IF (showdist) THEN
                    index = onetothree(ib)
                    dp = drange(index(1))                
                    dcoll(index(1),2) = dcoll(index(1),2) + weight(ib)/(dp**3)
                    vcoll(index(2),2) = vcoll(index(2),2) + weight(ib)
                    angcoll(index(3),2) = angcoll(index(3),2) + weight(ib)
                ENDIF
            ENDDO 

            voll = sum(vollay)
            ! The volume of ejecta in one timestep of conebuilding is such that the total volume grows by incrvol
            coneh(:) = coneh(:) + layer(:)*incrvol/voll
            volume = volume + incrvol
            IF (avalanche) THEN
                DO k1=1,Ng-1
                    IF ((coneh(k1)-coneh(k1+1))>dhmax) THEN
                        dhava = (coneh(k1)-coneh(k1+1)-dhmax)
                        coneh(k1) = coneh(k1) - dhava
                        ! The last factor accounts for the the ~r in the ring area (it preserves the total volume)
                        coneh(k1+1) = coneh(k1+1) + dhava*(conex(k1)/conex(k1+1))
                    ENDIF
                ENDDO
            ENDIF
            IF (backlanche) THEN
                DO k1=Ng,2,-1
                    IF ((coneh(k1)-coneh(k1-1))>dhmax) THEN
                        dhava = (coneh(k1)-coneh(k1-1)-dhmax)
                        coneh(k1) = coneh(k1) - dhava
                        coneh(k1-1) = coneh(k1-1) + dhava*(conex(k1)/conex(k1-1))
                    ENDIF
                ENDDO
            ENDIF
            IF (checkvolume) THEN
                checkvol = 0.
                DO k1=1,Ng
                    checkvol = checkvol + 2.*pi*conex(k1)*dx*coneh(k1)
                ENDDO
                print *,'After avalanches: ', checkvol, 'error [%]: ', (checkvol-volume)/volume*100.
            ENDIF
            IF(coneh(1) >= conetop) print *,'Desired height of the cone reached'
            IF(volume >= maxvol) print *,'Desired volume of the cone reached'
            IF((coneh(1) >= conetop).or.(volume >= maxvol)) hotovo=.true.
            print *,Nb,' ballistics computed, crater ', coneh(1),' [m] tall'
            ! Outputting cone shapes into cone.dat
            ! Shape created by shooting from flat base is baselayer, toplayer is ejected from the final cone
            IF (hotovo.or.ilay==1) THEN
                IF (ilay==1) THEN
                    open(51,file='run/baselayer.dat')
                    DO k1=1,Ng
                        write(51,*) conex(k1), layer(k1)*maxvol/voll
                    ENDDO
                    close(51)
                ENDIF
                IF (hotovo) THEN
                    open(51,file='run/cone.dat')
                    DO k1=1,Ng
                        write(51,*) conex(k1), coneh(k1)
                    ENDDO
                    close(51)
                    open(52,file='run/toplayer.dat')
                    DO k1=1,Ng
                        write(52,*) conex(k1), layer(k1)*maxvol/voll
                    ENDDO
                    close(52)                    
                ENDIF
            ENDIF
          ELSE  ! not building a cone
            hotovo = .true.
          ENDIF
          call cpu_time(t3)
          tcone = tcone + (t3-t2)

      ENDDO    ! loop until hotovo
      
      print *,'Hotovo. Total volume ',volume/1.e9, ' [km3], height ', coneh(1), ' [m]'      
      print *,'Last layer: particles landed', sum(landed),', overshot ', sum(overshot), ' / ', Nb
      IF (plotball) THEN
          print *,'Plotting the ballistic curves'
          DO ib=1,Nb
            index = onetothree(ib)
            write(trnam,'(a,i0,a,i0,a,i0)') 'traj',index(1),'_',index(2),'_',index(3)
            open(52,file='run/'//trim(trnam)//'.dat')
            DO k1=1,MAXSTP
                write(52,*) traj(:,k1,ib), drange(index(1)), vrange(index(2)), angrange(index(3))
            ENDDO
            close(52)
          ENDDO
      ENDIF
      IF (showdist) THEN
          open(66,file='run/dcoll.dat'); open(67,file='run/vcoll.dat'); open(68,file='run/angcoll.dat')
          DO k1=1,Nd
            write(66,*) dcoll(k1,1), dcoll(k1,2)
          ENDDO
          DO k2=1,Nv
            write(67,*) vcoll(k2,1), vcoll(k2,2)
          ENDDO
          DO k3=1,Na
            write(68,*) angcoll(k3,1), angcoll(k3,2)
          ENDDO
          close(66); close(67); close(68);
      ENDIF
        
      print *
      print *,'-----------------------------------------------------------'
      print *,'COMPUTATIONAL TIME ANALYSIS'
      print *,'computing ballistics  ', tball
      print *,'bulding the cone ', tcone

    CONTAINS

      SUBROUTINE derivs(time,xvec,dxvecdt,dp)
      REAL, INTENT(IN) :: time, xvec(nvar), dp
      REAL, INTENT(OUT) :: dxvecdt(nvar)
      REAL :: vamp
      ! RHS of the ballistic differential equation, xvec = (x,z,vx,vz)
      ! Eq. (1) from Broz et al., 2014, resp. Eq. 7.4b from Parfitt & Wilson, Fundamentals of Phys. Volc.
      ! dx/dt = v
      ! dv/dt = -g - (3 Cd rhoair)/(4 d rhorock) |v|v
         vamp = sqrt(xvec(3)*xvec(3) + xvec(4)*xvec(4))
         dxvecdt(1) = xvec(3)
         dxvecdt(2) = xvec(4)
         dxvecdt(3) = -(3.*Cd*rhoair)/(4.*dp*rhorock)*vamp*xvec(3)
         dxvecdt(4) = -grav - (3.*Cd*rhoair)/(4.*dp*rhorock)*vamp*xvec(4)
      END SUBROUTINE derivs
      
      FUNCTION onetothree(iin)
      INTEGER, INTENT(IN) :: iin
      INTEGER :: onetothree(3)
        onetothree(3) = mod(iin-1, Na) + 1
        onetothree(2) = mod((iin-1) / Na, Nv) + 1
        onetothree(1) = ((iin-1) / (Nv * Na)) + 1
      END FUNCTION onetothree
      
      SUBROUTINE loadcone(xcor,zcor)
      REAL,INTENT(OUT) :: xcor(Ng),zcor(Ng)
      REAL :: hlpr(2)
      INTEGER :: ierr,iread=0
        open(47,file=lconenam)
        DO WHILE(ierr>=0)
            read(47,fmt=*,iostat=ierr) hlpr(1),hlpr(2)
            IF (ierr==0) THEN
                iread = iread + 1 
                xcor(iread) = hlpr(1)
                zcor(iread) = hlpr(2)
            ENDIF
        ENDDO
        print *,'Cone read from file, iread ', iread
        IF(iread /= Ng) stop 'Cone read incorrectly'    
        close(47)
      END SUBROUTINE loadcone
      
      SUBROUTINE LIve()
        print *
        print *,'---------o-------BallPark'        
        print *,'------o------------------'
        print *,'----o------x**x----------'
        print *,'---o----x-------x--------'
        print *,'--o--x------------x------'
        print *,'-o-x---------------x-----'
        print *,'ox------------------x--VP'
        print *
        open(22,file='info.txt')
        write(22,*) "Nd =", Nd
        write(22,*) "Nv =", Nv
        write(22,*) "Na =", Na
        write(22,*) "conetop =", int(conetop)
        close(22)
      END SUBROUTINE LIve      

      END PROGRAM Ballpark
