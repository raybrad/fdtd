Program FDTD
use parameters
use allocmemo
use inputobj
use timeupdate
use output
use source
implicit none
integer::i
real*8::TF_SF_box(3,2),radius0,radius1,height,gap,centerx,centery,centerz,thickness,zpos,molex,moley,molez,molesize1,molesize2,molesize3
real*8::x0,y0,z0
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
!  1: Normal Metal  input eps_b sigma
!  2: Drude Metal   input eps_b gamma_p omega_p
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

if( metal_type == 1 .or. metal_type == 2 ) then
radius0=10.d-9		!10nm
write(6,*) 'radius0',radius0*1d9,'(nm)'
size1=radius0/l_n
centerx=0.d0
centery=0.d0
centerz=0.d0

!call input_metal_para('sphere',centerx,centery,centerz,0.d0,size1,0.d0,0.d0,eps_b,omg_p,gam_0,0.d0,l_n)

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
TF_SF_box(1,1)= -0.30	!small
TF_SF_box(1,2)=  0.30	!large
TF_SF_box(2,1)= -0.30	!small
TF_SF_box(2,2)=  0.30	!large
TF_SF_box(3,1)= -0.30	!small
TF_SF_box(3,2)=  0.30	!large 
light_pos= floor(0.5+((0.30+zcenter)*lattice_z))
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

!if(mod(t,interval)==0) then
call out_point('Ex',0.d0,0.d0,TF_SF_box(3,1)+0.01,0,dd,'up_Ex.dat')
call out_point('Ex',0.d0,0.d0,TF_SF_box(3,2)-0.01,0,dd,'down_Ex.dat')

call out_point('E^2',0.275d0,0.d0,0.d0,0,dd,'11_E2.dat')
call out_point('Ex',0.275d0,0.d0,0.d0,0,dd,'11_Ex.dat')
call out_point('Ex',0.300d0,0.d0,0.d0,0,dd,'12_Ex.dat')
call out_point('Ex',0.325d0,0.d0,0.d0,0,dd,'13_Ex.dat')
call out_point('Ex',0.350d0,0.d0,0.d0,0,dd,'14_Ex.dat')
call out_point('Ex',0.375d0,0.d0,0.d0,0,dd,'15_Ex.dat')
!endif
enddo




write(6,*) 'Calculation Completed'


end program fdtd

