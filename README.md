# Pulse_Profile_Simulaiton
Fortran code for a relativistic pulse profile simulation from a rotating neutron star with one or two hot spots on its surface. The core of this code is based on the classic approach described by Pechenick, Ftaclas, &amp; Cohen (1983, ApJ, 274, 846), which models light bending and relativistic beaming in Schwarzschild geometry.

#compile 
#bash
gfortran -o pulse pulse.f common.f

#run
./pulse
