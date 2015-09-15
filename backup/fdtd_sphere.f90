Program FDTD
use parameters
use allocmemo
use inputobj
use timeupdate
use output
use source
implicit none
integer::i
real*8::TF_SF_box(3,2),radius0,radius1,height,gap,centerx,centery,centerz,thickness
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
!--define metal structure and its property------------------------------------------------------------!
!cshape,centerx,centery,centerz,centerxl,size1,size2,size3,epsilon_b,omega_p,gamma_0,lattice_n        !       
!block                                                                                                !       
!call input_metal_para("block",-0.10,0.d0,0.d0,0.d0,0.2,0.1,0.1,eps_b,omg_p,gam_0,sigm_0,l_n)         !       
!80-17-17cone                                                                                         !       
!  --------------------------\		/-------------------                                          !       	
!	                      \        /                                                              !        
!		(cx)	      /(cl)    \                                                              !        
!  --------------------------/          \-------------------                                          !       
!call input_metal_para('pyramid',-1.2d-1,0.d0,0.d0,-1.0d-2,2.0d-1,2.0d-2,2.0d-2,eps_b,omg_p,gam_0,l_n)!       
!call input_metal_para('pyramid', 1.2d-1,0.d0,0.d0, 1.0d-2,2.0d-1,2.0d-2,2.0d-2,eps_b,omg_p,gam_0,l_n)!       
!-----------------------------------------------------------------------------------------------------!

!--Metal structure
!call input_metal_para(csp,centerx1,centery1,centerz1,centerxl1,size1,size2,size3,eps_b,omg_p,gam_0,sigm_0,l_n)
!call input_metal_para(csp,centerx2,centery2,centerz2,centerxl2,size1,size2,size3,eps_b,omg_p,gam_0,sigm_0,l_n)


if( metal_type == 1 .or. metal_type == 2 ) then
!Hex Au metal cylinder
radius0=95.d-9		!center
radius1=56.d-9		!around
height =30.d-9		
gap=15.d-9	!gap=150.d-9

write(6,*) 'radius0,radius1,height,gap',radius0,radius1,height,gap
!--Au disk in center
size1=radius0/l_n
size2=height/l_n
centerz=0.d0
centerx=0.d0
centery=0.d0

call input_metal_para('sphere',0.d0,0.d0,0.d0,0.d0,size1,0.d0,0.d0,eps_b,omg_p,gam_0,0.d0,l_n)
!call input_metal_para('rod',centerx,centery,centerz,0.d0,size1,size2,0.d0,eps_b,omg_p,gam_0,0.d0,l_n)
!--Au disk around
size1=radius1/l_n
size2=height/l_n
do i=1,8
centerx=(radius0+gap+radius1)*cos(pi/4.0*(i-1))/l_n
centery=(radius0+gap+radius1)*sin(pi/4.0*(i-1))/l_n
write(6,'(A5,i,A15,F8.3,F8.3)') 'disk ',i,'centerx,centery',centerx,centery
!call input_metal_para('rod',centerx,centery,centerz,0.d0,size1,size2,0.d0,eps_b,omg_p,gam_0,0.d0,l_n)
enddo

!--Graphene layer
thickness=10.d-9
centerx=0.d0
centery=0.d0
centerz=0.d0-(height*0.5+thickness*0.5)/l_n
size1=0.4
size2=0.4
size3=thickness/l_n
!!call input_metal_para('block',centerx,centery,centerz,0.d0,size1,size2,size3,eps_b,omg_p,gam_0,0.d0,l_n)
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
size1=(TF_SF_box(1,2)-0.05)*2.d0
size2=(TF_SF_box(2,2)-0.05)*2.d0
size3=(TF_SF_box(3,2)-0.05)*2.d0
call monitor(0.d0,0.d0,0.d0,size1,size2,size3)
size1=(TF_SF_box(1,2)+0.05)*2.d0
size2=(TF_SF_box(2,2)+0.05)*2.d0
size3=(TF_SF_box(3,2)+0.05)*2.d0
call monitor(0.d0,0.d0,0.d0,size1,size2,size3)
!--Time propagation	
 do t =0,dd
!  call   Gaussian_planewave(Ef,Hf,pos,freq,tzero,del_t)
   call   Gaussian_planewave_TF_SF(TF_SF_box,TF_SF_eps_bg,TF_SF_freq,TF_SF_to,TF_SF_tdecay)
   
   call   propagate()



!--QM Calculation   
  if (lemtoqm) then
    call lodestar_qmcal(0.d0,0.d0,0.d0)
  endif

!-plane & point  probe



  if (mod(t,interval)==0 ) then
    call  out_plane('Ex', 'x',0.0d0,'xEx.dat')	
    call  out_plane('Hy', 'y',0.0d0,'yHy.dat')	
    call  out_plane('E^2','z',0.0d0,'zE2.dat')
    call  out_point_w

    size1=(TF_SF_box(1,2)-0.05)*2.d0
    size2=(TF_SF_box(2,2)-0.05)*2.d0
    size3=(TF_SF_box(3,2)-0.05)*2.d0
    call incident_intensity(TF_SF_box(1,2)-0.05,0.d0,0.d0,size1,size2)
    call Poynting_block_in('absorb.dat',0.d0,0.d0,0.d0,size1,size2,size3,0,dd)

    size1=(TF_SF_box(1,2)+0.05)*2.d0
    size2=(TF_SF_box(2,2)+0.05)*2.d0
    size3=(TF_SF_box(3,2)+0.05)*2.d0
    call Poynting_block_out('scatter.dat',0.d0,0.d0,0.d0,size1,size2,size3,0,dd)
  end if

call out_point('Ex',0.d0,0.d0,TF_SF_box(3,1)+0.05,0,dd,'up_Ex.dat')
call out_point('Ex',0.d0,0.d0,TF_SF_box(3,2)-0.05,0,dd,'down_Ex.dat')

enddo




write(6,*) 'Calculation Completed'


end program fdtd

