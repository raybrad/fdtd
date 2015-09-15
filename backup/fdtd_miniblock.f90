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
call parameter_check

!------------------------------------------------------------!
!--------------Set Default Parameters And Generate Grid------!
!------------------------------------------------------------!
!--Input  S to guarantee the Courant stability condition
!S>sqrt(1/(dx)^2+1/(dy)^2+1/(dz)^2)
call set_default_parameter(2.d0)


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
!Hex Au metal block
x0=6.d-9		!radius0=10.d-9
y0=6.d-9		!height = 6.d-9	
z0=6.d-9
write(6,*) 'radius0,height',radius0,height
!--Au disk in center
size1=x0/l_n
size2=y0/l_n
size3=z0/l_n

centerx=0.d0
centery=0.d0
centerz=0.d0

call input_metal_para('block',centerx,centery,centerz,0.d0,size1,size2,size3,eps_b,omg_p,gam_0,0.d0,l_n)

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

!-plane & point  probe


  if (mod(t,interval)==0 ) then
    call  out_plane_cut('Ex','z',0.d0,'cenEx.dat')
    call  out_plane_cut('Ey','z',0.d0,'cenEy.dat')
    call  out_plane_cut('Ez','z',0.d0,'cenEz.dat')
    call  out_plane_cut('E^2','z',0.d0,'cenE2.dat')
  end if

call out_point('Ex',0.d0,0.d0,TF_SF_box(3,1)+0.01,0,dd,'down_Ex.dat')
call out_point('Ex',0.d0,0.d0,TF_SF_box(3,2)-0.01,0,dd,'up_Ex.dat')
call out_point('Ex',0.d0,0.d0,0.d0,0,dd,'inmetal_Ex.dat')
call out_point('E^2',0.d0,0.d0,0.d0,0,dd,'inmetal_E2.dat')
!call out_line_x('E^2',-0.2d0,0.2d0,0.d0,0.d0,0,dd,'outlinexE2.dat')
call out_point('E^2',-0.1875d0,0.d0,0.d0,0,dd,'out_0.1875.dat')
call out_point('E^2',-0.18d0,0.d0,0.d0,0,dd,'out_0.18.dat')
call out_point('E^2',-0.20d0,0.d0,0.d0,0,dd,'out_0.20.dat')
enddo




write(6,*) 'Calculation Completed'


end program fdtd

