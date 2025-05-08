      Real*8 pi
      Parameter (pi=3.1415927d0)

      Real*8 twospot
      Parameter (twospot=0.0)      ! one spot. set twospot=1.0 for two spots

      Real*8 p                     ! p=R/2M for the neutron star
      Real*8 size                  ! cos of angular size of each hot spot
      Real*8 alpha                 ! cos of angular distance of hot spot
                                   !    from rotational equator
      Real*8 beta                  ! cos of angular distance of observer
                                   !    from rotational equator
      Real*8 nspot                 ! # of hot poles on NS surface  
      Real*8 E                     ! photon energy

      Integer Nmu,Neninp,Int_en,ntheta
      parameter(ntheta=101)
      Parameter (Nmu=31)          ! # of points for beaming angle
      Parameter (Neninp=1001)      ! # of points for input energy spectrum
      Real*8 Energy,Mu,Inp_beam,Inp_beam_int
      real*8 Inttable

      Integer iflag                ! flag for the integral in Theta_int
      Integer Nchi
      Parameter (Nchi=200)         ! number of points in x
      Real*8 Xx,Thetaofx           ! arrays to store theta(x) integration

      Common/NS/p,size,alpha,beta,nspot
      Common/Beam1/E,Int_en
      Common/Beam2/Energy(Neninp),Mu(Nmu),Inp_beam(Neninp,Nmu)
      Common/Beam3/Inp_beam_int(Nmu)
      Common/Beam4/Inttable(ntheta,Neninp,Nmu)

      Common/Flag/iflag
      Common/Lookup/Xx(Nchi),Thetaofx(Nchi)



