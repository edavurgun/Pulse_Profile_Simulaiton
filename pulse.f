c     *****************************************************************
c     *           Program pulse.f
c     * ============ Version 1 =================
c     * 05/29/00: Calculates pulse profiles from a rotating compact
c     *           object, according to the method described in
c     *           Pechenick et al. 1983, ApJ, 274, 846. The rotation
c     *           is assumed to be slow, so that the exterior metric
c     *           is Schwarzschild, and photon-energy dependent 
c     *           information is taken into account. (v1.0, DP)
c     * 05/30/00: Debugged (v1.1, DP)
c     * ============ Version 2 =================
c     * 05/30/00: The function theta(x) is now calculated only once
c     *           in the beginning of the program and stored in an
c     *           array to be interpolated, improving performance.
c     *           (v2.0, DP)
c     * 06/02/00: Debugged a problem with the # of grid points 
c     *           when the size of the hot spot is very small. (DP)
c     * ============ Version 3 =================
c     * 06/02/00: The energy-dependent beaming is now read from
c     *           a user-specified file and interpolated. 
c     *           NOTE! The number of energy and mu grid points
c     *           must be known and specified a priori. (DP)
c     * ============ Version 4 =================
c     * 01/10/01: The integral over the angle x is now performed over
c     *           a linear grid but the interpolation table theta(x)
c     *           remains on a log grid. This will hopefully remove
c     *           the problem with small spot sizes and small phase angles
c     *           without requiring too many x points. (DP&FO)
c     ****************************************************************
c     * !!!NB!!!  If the spot size is small, increase the number of
c     *           grid points in the chi and phi integrals
c     ****************************************************************

 
c     ****************** MAIN PROGRAM ****************

      Implicit None

      character*40 fname

      Real*8 A0,th0,phase,a
      Integer Nph,j
      Parameter (Nph=9)                 ! # of points for pulse phase

      Real*8 Emin,Emax,E_inf
      Integer Nen,i
      Parameter (Nen=64)                 ! # of points for energy grid

      real*8 phimax
      parameter(phimax = 180.0)

      Real*8 fmax,fmin, plmax, plmin
      real*8 sum, average, power

      Real*8 Flux(Nph,Nen)              ! phase-resolved spectra
      include 'common.f'

                                        ! Initialize Parameters
      p=3.2                           ! R/2M of neutron star
      
      size=cos(70.0*pi/180.0)         ! opening of hot spot
      alpha=cos(90.0*pi/180.0)        ! angle between rot and B axes
      beta=cos(45.0*pi/180.0)         ! polar angle of observer

c      N_beam=0.0000000001d0             ! beaming exponent (not used)
      Emin=0.5                        ! min Energy (min E desired for output spec)
      Emax=6.5                        ! max Energy (max E desired for output spec)

      fname='I_mu_E_B2D14_T0.4.out'
      call Readfile(fname)              ! read input beaming
      call Calctheta                    ! calculate Theta(x) table
      write (*,*) 'Theta Calculation is done'

      do i=1,Nen-1                      ! for all photon energies
                                        ! current energy
         E_inf=10.0**(dlog10(Emin)+dfloat(i-1)/dfloat(Nen-1)
     +        *(dlog10(Emax)-dlog10(Emin)))
                                        ! energy at NS surface
         E=E_inf/dsqrt(1.0-1.0/p)
         call Int_energy                ! find grid point in input grid

         plmin=1.0+30
         plmax=0.0
         average=0.0
         do j=1,Nph                     ! for all phases
            write (*,*) 'Energy grid #: ',i,'/',Nen,'   Phase grid #: ',
     +           j,'/',Nph
                                 ! phase<pi/2 because of symmetry
            phase=phimax*dfloat(j-1)/dfloat(Nph-1)   
                                        ! calculate angle theta_0
            th0=dsqrt(1.0-alpha*alpha)*dsqrt(1.0-beta*beta)
     +           *cos(phase*pi/180.0)+alpha*beta
                                        ! calculate transfer function

            if (th0.GE.1.0) th0=0.999990
            Flux(j,i)=A0(th0) + twospot*A0(-th0)  ! one or two spots
            if (Flux(j,i).LT.plmin) plmin=Flux(j,i)
            if (Flux(j,i).GT.plmax) plmax=Flux(j,i)
            average = average + Flux(j,i)
         end do                 ! phase loop
         average = average/Nph

         power=0.0
         sum=0.0
         do j=1,Nph
            phase=phimax*dfloat(j-1)/dfloat(Nph-1)   
            if(j.eq.1.or.j.eq.Nph) then
               sum=sum+0.5*flux(j,i)*phimax/dfloat(Nph-1)
            else
               sum=sum+flux(j,i)*phimax/dfloat(Nph-1)
            endif
            power = power+(Flux(j,i)-average)**2.0

            write (18,*) phase/360.0, E_inf, Flux(j,i) !pulse prof&spec

         end do
         sum=sum/phimax
         power=dsqrt(power/Nph)/average

         write(16,*) E_inf, sum                   ! phase-avgd spec
         write(17,*) E_inf, (plmax-plmin)/(plmax+plmin) !PF vs E
         write(19,*) E_inf, power                  ! rms vs E 



      end do                    ! energy loop

      end


c     ******** Subroutine Readfile(fname) **************************
c     * 06/02/00: Reads the energy-dependent beaming from the file
c     *           given in fname. The number of grid points in both
c     *           energy and mu are assumed known. The file is assumed
c     *           to have the generic format: Energy, mu, Intensity
c     *           with the mu index varying first. (DP)
c
c     * 04/27/01: Columns 3 and 4 are read individually from the 
c                 input file and summed to get total intensity.
c                 Also, theta is read from file and the cosine computed   
c                 in this subroutine. (DP&FO)    
c     **************************************************************
      Subroutine Readfile(fname)

      Implicit None
      Integer i,j
      real*8 I1, I2, th
      character*30 fname
      include 'common.f'
      
      open (unit=1,file=fname,status='OLD')
      do i=1,Neninp                          ! Real all energies
         do j=1,Nmu                          ! and all beaming angles
            read (1,*) Energy(i), th, I1, I2
            Inp_beam(i,j) = I1 + I2
            Mu(j) = cos(th*pi/180.0)
         end do
      end do
      close(1)

*      write (*,*) 'Read Beaming File'

      end

c     ******** Subroutine Int_energy *******************************
c     * 06/02/00: Given a current energy, E, passed through the 
c     *           common block, this subroutine interpolates the
c     *           input energy grid Energy(Neninp) and returns
c     *           an index to be used in the follow-up (see Intens) 
c     *           interpolations. This is done only once per energy
c     *           point to speed up things. It uses a binary search
c     *           because the energy grid is sorted in increasing 
c     *           order and returns the index in the common block
c     *           as Int_en. (DP)
c     ***************************************************************
      Subroutine Int_energy

      Implicit None
      Real*8 intmax,intmin,intmid
      include 'common.f'

      intmin=1                           ! lower bound of current range
      intmax=Neninp                      ! upper bound of current range
      do while (intmax-intmin.GT.1)      ! while range has more than one point
         intmid=(intmax+intmin)/2        ! find mid-point
         if (Energy(nint(intmid)).GE.E) then
            intmax=intmid                ! if E>Energy(midpoint) search upper
         else
            intmin=intmid                ! if E<Energy(midpoint) search lower
         end if
      end do
      Int_en=intmax                      ! found it!
      end

c     ******** Real*8 Function Intens(thp,delta) *******************
c     * 05/29/00: Returns the specific intensity (measured locally)
c     *           that emerges from the neutron-star surface at photon
c     *           energy E (measured locally), at polar angle thp
c     *           (actually its cosine), and at an angle delta
c     *           (actually its cosine) with respect to the surface
c     *           surface normal. The polar angles are measured with
c     *           respect to the magnetic (or whatever else) axis.
c     *           The units of photon energy and of the intensity 
c     *           are arbitrary. (DP)
c     * 06/02/00: Added an interpolation routine in case the 
c     *           energy-dependent beaming is read from a file. (DP)
c     * 2002    : Added an option for dipole cooling (FO)
c     ****************************************************************
      Real*8 Function Intens(thp,delta)

      Implicit None
      Real*8 thp,delta, T, Tpole, Inten1,Inten2
      Integer intmin,intmax,intmid,int
      include 'common.f'

      if (thp.GE.size) then               ! if inside *one* hot spot

c     *** This is for one simple hot spot with analytic beaming***
c         Intens=sqrt(1.D0-delta**2.D0)**N_beam  ! use a cos^N_beam beaming
c
c     *** This is for a Heyl & Hernquist thermal cooling NS
c         Tpole=0.086d0
c         T=Tpole*(4.d0**0.8d0*thp*thp/
c     +        (3.d0*thp*thp+1.d0)**0.8d0)**0.25d0
c         Intens=E**3.d0/(dexp(E/T)-1.d0)
c         Intens=Intens*(delta**N_beam)

c     *** This is for interpolating the beaming angle
         intmin=1                          ! lower bound of current range
         intmax=Nmu                        ! upper bound of current range
         do while (intmax-intmin.GT.1)     ! while range has more than one point
            intmid=(intmax+intmin)/2       ! find mid-point
            if (Mu(intmid).GE.delta) then
               intmax=intmid               ! if d>mu(midpoint) search upper
            else
               intmin=intmid               ! if d<mu(midpoint) search lower
            end if
         end do
         Int=intmax                        ! found it!
                                           ! interpolate energies first
         Inten1=(Inp_beam(Int_en,Int)-Inp_beam(Int_en-1,Int))
     +        /(Energy(Int_en)-Energy(Int_en-1))*
     +        (E-Energy(Int_en-1))+Inp_beam(Int_en-1,Int)
         Inten2=(Inp_beam(Int_en,Int-1)-Inp_beam(Int_en-1,Int-1))
     +        /(Energy(Int_en)-Energy(Int_en-1))*
     +        (E-Energy(Int_en-1))+Inp_beam(Int_en-1,Int-1)
                                           ! now interpolate angles
         Intens=(Inten1-Inten2)/(Mu(int)-Mu(int-1))*(delta-Mu(int-1))
     +        +Inten2
c         write (*,*) delta,E,intens,inten1,inten2
c         write (*,*) int_en,int
c         write (*,*) mu(int),mu(int-1)
c         write (*,*) Inp_beam(int_en,int),Inp_beam(int_en-1,int)
c         write (*,*) Inp_beam(int_en,int-1),Inp_beam(int_en-1,int-1)
c         read (*,*)
 111     continue
      else                      ! if outside the hot spot
         Intens=0.0            ! no radiation
      end if
      
      end

c     *********** Real*8 Function Func_theta(uort,x,dum) ************
c     * 05/29/00: Integrand of the function required by Theta_int
c     *           iflag is passed through the common block and is
c     *           described in Theta_int. (DP)
c     ***************************************************************
      Real*8 Function Func_theta(uort,x,dum)

      Implicit None
      Real*8 uort,x,dum
      include 'common.f'

      if (iflag.eq.0) then                ! if x<xmax
         Func_theta=1.0/x/x-uort*uort*(1.0-2.0*uort)
         Func_theta=1.0/dsqrt(Func_theta)
      else                                ! if x=xmax
         Func_theta=-2.0*(0.50/p-uort*uort)**2.0+
     +        (1.0-1.0/p)*(1.0/p-uort*uort)
         Func_theta=2.0/dsqrt(Func_theta)
      end if

      end

c     *********** Subroutine Calctheta ******************************
c     * 05/30/00: Calculates the function theta(x), eq. 3.9 in
c     *           Pecheneck et al. 1983. Given a polar angle x at the
c     *           the observer, theta is the polar angle on the 
c     *           neutron star surface that results from tracing
c     *           backwards the light bending. The integral is
c     *           improper for x=xmax (eq.3.12 of Pechenick et al.)
c     *           and requires a special treatment. If x<xmax, then
c     *           iflag=0 and the standard integral is performed.
c     *           If x=xmax, iflag<>0 and the integral is performed
c     *           after the change of variables t=sqrt(M/R-u)
c     *           that eliminates the pole. The results are stored
c     *           in the array Thetaofx(Nchi) to be used for future
c     *           interpolation. The abcisa of the array is stored
c     *           in Xx(Nchi). (DP)
c     ***************************************************************
      Subroutine Calctheta

      Implicit None
      Real*8 x,xlog,xlogmax,xlogmin,xlogstep,xmax
      Real*8 uortmax,theta
      Integer i
      include 'common.f'
      External Func_theta

      xmax=2.0*p/dsqrt(1.0-1.0/p)    ! upper limit of x integration

      xlogmax=dlog10(xmax)              ! for log grid in (xmax-x)
      xlogmin=-3.0                     ! closest point to xmax
      xlogstep=(xlogmax-xlogmin)/dfloat(Nchi-2)  ! Nchi-1 interior points

      Xx(1)=0.0                        ! first impact parameter
      Thetaofx(1)=0.0                  ! no deflection!
      iflag=0                           ! for all interior points iflag=0
      do i=2,Nchi-1                     ! for all interior points
                                        ! integrand for i=1 is zero
                                        ! integrand for i=Npt is below

         xlog=xlogmax-dfloat(i-1)*xlogstep    ! current log value
         x=xmax-10.0**(xlog)                 ! current x value
         Xx(i)=x/xmax                   ! store abcisa
         uortmax=0.50/p                ! umax=M/R
                                        ! Integrate  
         call qtrap(Func_theta,0.0,uortmax,theta,1.0-5,x,0.0)      
         Thetaofx(i)=theta              ! store value
      end do
      iflag=1                           ! for last point
      Xx(Nchi)=1.0                     ! store absica
      uortmax=dsqrt(0.50/p)            ! tmax=sqrt(M/R)
                                        ! Integrate  
      call qtrap(Func_theta,0.0,uortmax,theta,1.0-5,x,0.0)      
      Thetaofx(Nchi)=theta
      end


c     *********** Real*8 Function Theta_int(x) **********************
c     * 05/30/00: Interpolates the arrays calculated by Calctheta
c     *           and returns the values of theta for each x.
c     *           Note that in this subroutine x is actually
c     *           x/xmax, i.e., the variable of interpolation. 
c     *           Since the array Xx(Nchi) is sorted, the subroutine
c     *           uses a binary search. (DP)
c     ***************************************************************
      Real*8 Function Theta_int(x)

      Implicit None
      Real*8 x
      Integer int,intmax,intmin,intmid
      include 'common.f'

c     ************************* old version with serial search ************
c      int=1                              ! initialize for ordinary search
c      do while (int.LT.Nchi.AND.Xx(int).LT.x)    ! are we there yet?
c         int=int+1                       ! if not increase counter
c      end do
c     *********************************************************************

      intmin=1                           ! lower bound of current range
      intmax=Nchi                        ! upper bound of current range
      do while (intmax-intmin.GT.1)      ! while range has more than one point
         intmid=(intmax+intmin)/2        ! find mid-point
         if (Xx(intmid).GE.x) then
            intmax=intmid                ! if x>XX(midpoint) search upper half
         else
            intmin=intmid                ! if x<XX(midpoint) search lower half
         end if
      end do
      int=intmax                         ! found it!
                                         ! interpolate between int-1,int
      Theta_int=(Thetaofx(int)-Thetaofx(int-1))/(Xx(int)-Xx(int-1))*
     +     (x-Xx(int-1))+Thetaofx(int-1)
      end

c     *********** Real*8 Function Func_x(x_i,xmax) ******************
c     * 05/29/00: Integrand of the function required by A0.  
c     *           It uses a simple trapezoid to evaluate the inner
c     *           integral over phi, which is often trivial. (DP)
c     ***************************************************************
      Real*8 Function Func_x(x_i,xmax,th0)
      
      Implicit None
      Real*8 x_i,xmax,th0
      Real*8 delta,thp,thNS,Theta_int
      Real*8 sum,phi_i,phimax,stepphi,Intens
      Integer j,Npt
      Parameter (Npt=200)                ! number of points for phi integral
      include 'common.f'

      delta=dsqrt(1.0-(x_i/xmax)**2.0) ! beaming angle

      if (x_i.GT.0.0) then              ! for x>0
         thNS=Theta_int(x_i/xmax)        ! theta w.r.t observer
      else
         thNS=0.0                       ! for x=0=>th_NS=0
      end if

      sum=0.0                           ! initialize sub integral

                                         ! Try to obtain limits on phi
      if (th0.EQ.1.0.or.sin(thNS).eq.0.0.OR.size.eq.0.0) then ! special
         phimax=pi                       ! only 0<phi<pi because of symmetry
      else
                                         ! find intersections with spot
         phimax=(size-cos(thNS)*th0)/sin(thNS)/dsqrt(1.0-th0*th0)

         if (phimax.GE.-1.0.AND.phimax.LE.1.0) then
            phimax=acos(phimax)
         else                            ! do all phi if limits didn't work out
            phimax=pi
         end if
      end if
      
      if (phimax.GT.0.d0) then            ! only if there are points
         stepphi=phimax/dfloat(Npt-1)     ! step in phi
         
         do j=1,Npt                       ! do only half phi>0
            phi_i=dfloat(j-1)*stepphi     ! phi w.r.t. observer
            
                                          ! cos of angle wrt. magn axis
            thp=sin(thNS)*dsqrt(1.0-th0*th0)*cos(phi_i)
     +           +cos(thNS)*th0
            
            if (j.eq.1.or.j.eq.Npt) then  ! boundary points
               sum=sum+0.50*stepphi*Intens(thp,delta)
            else                          ! interior points
               sum=sum+stepphi*Intens(thp,delta)
            end if
            
         end do                           ! phi-loop
         sum=2.0*sum                     ! for other half phi values
      end if        
      Func_x=sum*x_i                      ! integrand
      
      end

c     *********** Real*8 Function A0(th0) ****************************
c     * 05/29/00: Calculates something like the function A of eq.3.16
c     *           in Pechenick et al. 1983. th0 is the angle 
c     *           (actually its cosine) between the magnetic axis and
c     *           the line that connects the oserver with the center
c     *           of the star. The parameters of the neutron star
c     *           are passed through common blocks. E is the photon
c     *           energy in arbitrary units from the common block as well. (DP)
c     *           NB!!! the difference with Pechenick et al. (1983)
c     *                 is a factor of (1-2M/R)^{1/2} because the
c     *                 specific intensity in this code is *not*
c     *                 integrated over photon energy.
c     * 01/10/01: The integral over x is now performed on a linear grid
c     ***************************************************************
      Real*8 Function A0(th0)
            
      Implicit None
      Real*8 th0,Func_x,Theta_int
      Real*8 xmax,xstep,x_i,dx_i
      Integer i,Npt
      Parameter (Npt=3000)               ! # of grid points for x integral
      include 'common.f'

      xmax=2.0*p/dsqrt(1.0-1.0/p)    ! upper limit of x integration

      xstep=xmax/dfloat(Npt-1)          ! N-1 interior points

      A0=0.0                           ! Initialize integral
      iflag=0                           ! for all interior points iflag=0
      do i=1,Npt-1                      ! for all interior points
                                        ! integrand for i=1 is zero
                                        ! integrand for i=Npt is below

         x_i=xstep*dfloat(i-1)                  ! current x value
         if (i.GT.1) then
                                                ! x_(i+1)-x(i-1) for trapezoid
            dx_i=2.0*xstep
         else
                                                ! (x_(2)-x(1)) for trapezoid
            dx_i=xstep
         end if
         A0=A0+0.50*Func_x(x_i,xmax,th0)*dx_i  ! trapezoid rule
      end do
      iflag=1
      A0=A0+0.50*Func_x(xmax,xmax,th0)*xstep   ! point # Npt
                          !vvv NB! difference with Pechenick et al. 1983
      A0=A0*(1.0-1.0/p)**1.50*(0.50/p)**2.0   ! normalization

      end


