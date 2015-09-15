Program FDTD
use parameters
use allocmemo
use inputobj
use timeupdate
use output
use source
implicit none
integer::i
real*8::TF_SF_box(3,2),radius0,radius1,height,gap,centerx,centery,centerz,thickness,zpos
real*8::x0,y0,z0
real*8::molex,moley,molez,molesize1,molesize2,molesize3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Define Basic Parameters!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------!
!-------------------Read Input-------------------------------!
!------------------------------------------------------------!
!--use uniform grid only
!--define box structure size
!--define lattice size
!--define PML size
!--define incident wavelength & time 
!--define background and metal eps
call readinput(5)

!------------------------------------------------------------!
!--------------Set Default Parameters And Generate Grid------!
!------------------------------------------------------------!
!--Input  S to guarantee the Courant stability condition
!S>sqrt(1/(dx)^2+1/(dy)^2+1/(dz)^2)
call set_default_parameter(S_factor)

call calWavePara

call parameter_check

!--Material Complexity
!  0: No Metal      input eps_b
!  1: Normal Metal  input eps_b sigma
!  2: Drude Metal   input eps_b gamma_p omega_p
!  3: Single Lorentz Metal  input eps_b omega_0 gamma_p omega_p
!  4: Double Lorentz Metal  input eps_b gamma_p omega_p
!  5: Drude + Double Lorentz Metal   input eps_b gamma_p omega_p
call material_complexity(metal_t)

!--Allocate Memory
call memory()

!--Real Space Parameters
call  real_space_param(l_n,freq)

!--Set background permittivity
call background(bg_eps)

!--source info:
call source_check


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!Define Device Structure!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------------------------------------!

if( metal_type > 0) then
radius0=10.d-9		!10nm
!radius0=25.d-9
write(6,*) 'radius0',radius0*1d9,'(nm)'
size1=radius0/l_n
centerx=0.d0
centery=0.d0
centerz=0.d0

call input_metal_para('sphere',centerx,centery,centerz,0.d0,size1,0.d0,0.d0,eps_b,omg_p,gam_0,sigm_0,&
		        depsr_0_1,omega_0_1,gamma_0_1,depsr_0_2,omega_0_2, gamma_0_2,l_n)


if(Mole_type >0) then
molex=(radius0+separation)/l_n
moley=0.d0
molez=0.d0

molesize1=0.50d-9/l_n  !5A
molesize2=0.50d-9/l_n  !5A
molesize3=0.50d-9/l_n  !2.5A
molex=molex+0.5d0*molesize1

write(6,'(A,X,3(E12.5,X))') 'molexyz:    ',molex,moley,molez
write(6,'(A,X,3(E12.5,X))') 'molesize123:',molesize1,molesize2,molesize3

volfactor=dble(molesize1*lattice_x*molesize2*lattice_y*molesize3*lattice_z)
if (volfactor<1.d0) then
volfactor=1.d0
endif
write(6,*) 'volfactor',volfactor
call Molecular_struct('block',molex,moley,molez,molesize1,molesize2,molesize3,1) !xrad
endif



endif


!--assign eps ,omega and others on metal grid

call assign_metal_grid()

!--out epsilon
call out_epsilon("x",0.d0,"epsilon.x")
call out_epsilon("y",0.d0,"epsilon.y")
call out_epsilon("z",0.d0,"epsilon.z")
!!!
!-coefficients to be used in timeupdate----!
call coefficient()


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FDTD Time Propagation!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--Field initialization
call field_initialization()

!--TF_SF size
! default
TF_SF_box(1,1)= -0.260	!small
TF_SF_box(1,2)=  0.260	!large
TF_SF_box(2,1)= -0.260	!small
TF_SF_box(2,2)=  0.260	!large
TF_SF_box(3,1)= -0.260	!small
TF_SF_box(3,2)=  0.260	!large 
light_pos= floor(0.5+((0.40+zcenter)*lattice_z))
!light_pos= floor(0.5+((-0.35+zcenter)*lattice_z))
write(6,*) 'light_pos',light_pos	

!--define monitor points

!--Time propagation	

 do t =0,dd

   call   Gaussian_planewave_TF_SF(TF_SF_box,TF_SF_eps_bg,TF_SF_freq,TF_SF_to,TF_SF_tdecay)
   
   call   propagate()



!--QM Calculation   
if (lemtoqm) then
  if(mod(t,10)==0) then
    call lodestar_qmcal(molex,moley,molez)
 endif
endif

if(mod(t,interval)==0) then
call out_point('Ex',0.d0,0.d0,TF_SF_box(3,1)+0.01,0,dd,'up_Ex.dat')
call out_point('Ex',0.d0,0.d0,TF_SF_box(3,2)-0.01,0,dd,'down_Ex.dat')
!xx
call out_point('Ex',0.000d0,0.d0,0.d0,0,dd,'0_Ex.dat')
call out_point('Ex',0.025d0,0.d0,0.d0,0,dd,'1_Ex.dat')
call out_point('Ex',0.050d0,0.d0,0.d0,0,dd,'2_Ex.dat')
call out_point('Ex',0.075d0,0.d0,0.d0,0,dd,'3_Ex.dat')
call out_point('Ex',0.100d0,0.d0,0.d0,0,dd,'4_Ex.dat')
call out_point('Ex',0.125d0,0.d0,0.d0,0,dd,'5_Ex.dat')
call out_point('Ex',0.150d0,0.d0,0.d0,0,dd,'6_Ex.dat')
call out_point('Ex',0.175d0,0.d0,0.d0,0,dd,'7_Ex.dat')
call out_point('Ex',0.200d0,0.d0,0.d0,0,dd,'8_Ex.dat')
call out_point('Ex',0.225d0,0.d0,0.d0,0,dd,'9_Ex.dat')
call out_point('Ex',0.2375d0,0.d0,0.d0,0,dd,'9.5_Ex.dat')
call out_point('Ex',0.250d0,0.d0,0.d0,0,dd,'10_Ex.dat')
call out_point('Ex',0.256d0,0.d0,0.d0,0,dd,'10.24_Ex.dat')
call out_point('Ex',0.258d0,0.d0,0.d0,0,dd,'10.32_Ex.dat')
call out_point('Ex',0.275d0,0.d0,0.d0,0,dd,'11_Ex.dat')
call out_point('Ex',0.300d0,0.d0,0.d0,0,dd,'12_Ex.dat')
call out_point('Ex',0.325d0,0.d0,0.d0,0,dd,'13_Ex.dat')
call out_point('Ex',0.350d0,0.d0,0.d0,0,dd,'14_Ex.dat')
call out_point('Ex',0.375d0,0.d0,0.d0,0,dd,'15_Ex.dat')
call out_point('Ex',0.400d0,0.d0,0.d0,0,dd,'16_Ex.dat')
call out_point('Ex',0.450d0,0.d0,0.d0,0,dd,'18_Ex.dat')
!!
call out_point('Ex',0.000d0,0.d0,0.d0,0,dd,'0_Ey.dat')
call out_point('Ex',0.025d0,0.d0,0.d0,0,dd,'1_Ey.dat')
call out_point('Ex',0.050d0,0.d0,0.d0,0,dd,'2_Ey.dat')
call out_point('Ex',0.075d0,0.d0,0.d0,0,dd,'3_Ey.dat')
call out_point('Ex',0.100d0,0.d0,0.d0,0,dd,'4_Ey.dat')
call out_point('Ex',0.125d0,0.d0,0.d0,0,dd,'5_Ey.dat')
call out_point('Ex',0.150d0,0.d0,0.d0,0,dd,'6_Ey.dat')
call out_point('Ex',0.175d0,0.d0,0.d0,0,dd,'7_Ey.dat')
call out_point('Ey',0.200d0,0.d0,0.d0,0,dd,'8_Ey.dat')
call out_point('Ey',0.225d0,0.d0,0.d0,0,dd,'9_Ey.dat')
call out_point('Ey',0.2375d0,0.d0,0.d0,0,dd,'9.5_Ey.dat')
call out_point('Ey',0.250d0,0.d0,0.d0,0,dd,'10_Ey.dat')
call out_point('Ey',0.256d0,0.d0,0.d0,0,dd,'10_24_Ey.dat')
call out_point('Ey',0.256d0,0.d0,0.d0,0,dd,'10_24_Ey.dat')
call out_point('Ey',0.258d0,0.d0,0.d0,0,dd,'10.32_Ey.dat')
call out_point('Ey',0.275d0,0.d0,0.d0,0,dd,'11_Ey.dat')
call out_point('Ey',0.300d0,0.d0,0.d0,0,dd,'12_Ey.dat')
call out_point('Ey',0.325d0,0.d0,0.d0,0,dd,'13_Ey.dat')
call out_point('Ey',0.350d0,0.d0,0.d0,0,dd,'14_Ey.dat')
call out_point('Ey',0.375d0,0.d0,0.d0,0,dd,'15_Ey.dat')
call out_point('Ey',0.400d0,0.d0,0.d0,0,dd,'16_Ey.dat')
call out_point('Ey',0.450d0,0.d0,0.d0,0,dd,'18_Ey.dat')
!!yy
call out_point('Ex',0.d0,0.000d0,0.d0,0,dd,'0_yEx.dat')
call out_point('Ex',0.d0,0.025d0,0.d0,0,dd,'1_yEx.dat')
call out_point('Ex',0.d0,0.050d0,0.d0,0,dd,'2_yEx.dat')
call out_point('Ex',0.d0,0.075d0,0.d0,0,dd,'3_yEx.dat')
call out_point('Ex',0.d0,0.100d0,0.d0,0,dd,'4_yEx.dat')
call out_point('Ex',0.d0,0.125d0,0.d0,0,dd,'5_yEx.dat')
call out_point('Ex',0.d0,0.150d0,0.d0,0,dd,'6_yEx.dat')
call out_point('Ex',0.d0,0.175d0,0.d0,0,dd,'7_yEx.dat')
call out_point('Ex',0.d0,0.200d0,0.d0,0,dd,'8_yEx.dat')
call out_point('Ex',0.d0,0.225d0,0.d0,0,dd,'9_yEx.dat')
call out_point('Ex',0.d0,0.2375d0,0.d0,0,dd,'9.5_yEx.dat')
call out_point('Ex',0.d0,0.250d0,0.d0,0,dd,'10_yEx.dat')
call out_point('Ex',0.d0,0.256d0,0.d0,0,dd,'10.24_yEx.dat')
call out_point('Ex',0.d0,0.258d0,0.d0,0,dd,'10.32_yEx.dat')
call out_point('Ex',0.d0,0.275d0,0.d0,0,dd,'11_yEx.dat')
call out_point('Ex',0.d0,0.300d0,0.d0,0,dd,'12_yEx.dat')
call out_point('Ex',0.d0,0.325d0,0.d0,0,dd,'13_yEx.dat')
call out_point('Ex',0.d0,0.350d0,0.d0,0,dd,'14_yEx.dat')
call out_point('Ex',0.d0,0.375d0,0.d0,0,dd,'15_yEx.dat')
call out_point('Ex',0.d0,0.400d0,0.d0,0,dd,'16_yEx.dat')
call out_point('Ex',0.d0,0.450d0,0.d0,0,dd,'18_yEx.dat')
call out_point('Ex',0.d0,0.d0,0.d0,0,dd,'00_Ex.dat')
endif

!if(mod(t,5000)==0 .and. (dble(t)/dd)>0.07d0 .and. (dble(t)/dd)<0.15d0) then
!    call  out_plane_cut('Ex' ,'z',0.d0,'xyEx.dat')
!    call  out_plane_cut('E^2','z',0.d0,'xyE2.dat')
!    call  out_plane_cut('E^2','y',0.d0,'xzE2.dat')
!endif
enddo




write(6,*) 'Calculation Completed'


end program fdtd

