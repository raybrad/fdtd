module variables

!/// in fdtd.f90 ///
!module var_main
real*8:: shift
integer:: SpecN
integer:: t
integer::del_t,dd,tx,tzero,Nf  !del_t: delta t;  dd:end time
real*8 ::Real_tzero,Real_deltat,Real_totalt	!fs
real*8 ::radiusSphere,measureLambda
real*8:: SumUpper, SumLower,volfactor
integer:: TF_SF_pos(3,2),light_pos
real*8:: TF_SF_freq=0.0 !default
integer:: TF_SF_to
integer:: TF_SF_tdecay
real*8:: TF_SF_eps_bg
real*8:: separation
logical::lqmtoem,lemtoqm,triple,dofft
logical,allocatable::lmonitor(:,:,:)
!/// in struct.f90 ///
! var_struct

integer:: isize, jsize, ksize
integer:: pmlil, pmlir, pmljl, pmljr, pmlkl, pmlkr
integer:: lattice_x, lattice_y, lattice_z
real*8::l_n
real*8:: xsize, ysize, zsize
real*8:: xcenter, ycenter, zcenter
real*8:: kx, ky, kz
real*8:: orderxl, orderyl, orderzl
real*8:: orderxr, orderyr, orderzr
real*8:: sig_axl, sig_ayl, sig_azl
real*8:: sig_axr, sig_ayr, sig_azr
real*8:: ds_x, ds_y, ds_z, dt
real*8:: S_factor
real,save::planckev,pi,eo,uo,ups,light_speed,bohr,el_charge
integer:: misize, mjsize, mksize
integer:: pisize, pjsize, pksize
integer:: cisize, cjsize, cksize
integer:: metal_type=0,metal_t       !default       
integer:: Mole_type=0
real*8::bg_eps,sigm_0,eps_b,omg_p,gam_0,depsr_0_1,omega_0_1,gamma_0_1,depsr_0_2,omega_0_2, gamma_0_2
real*8::centerx1,centery1,centerz1,centerxl1,size1,size2,size3,centerx2,centery2,centerz2,centerxl2
character*10::metal_shape,Ef,Hf,sourcetyp,Gauss_typ  !mat_file
real*8::pos,freq,wav_len,lambda_min,lambda_max,df,fmin,fmax     !!wc=a/wavelength ,normalized frequency of a resonance in question
character*10::ocomp,oplane,oname
integer::interval
real*8::ovalue,amplitude
!/// in memory.f90/
!module var_memo

integer,allocatable,save::cposition(:,:,:)
real*8,allocatable,save::  epsilonx(:,:,:),epsilony(:,:,:),epsilonz(:,:,:)
!/* Input dielectric structure for Drude or single pole Lorentz. */
real*8,allocatable,save::  Lepsilon(:,:,:),Lomega(:,:,:),Lgamma(:,:,:),Lomega0(:,:,:),Lsigma(:,:,:),Molecular(:,:,:)
!/* Input dielectric structure for double pole Lorentz. */
real*8,allocatable,save::  Lomega_1(:,:,:),Lgamma_1(:,:,:),Ldepsr_1(:,:,:),Lomega_2(:,:,:),Lgamma_2(:,:,:),Ldepsr_2(:,:,:)
!/* field components */
real*8,allocatable,save::  Ex(:,:,:),Ey(:,:,:),Ez(:,:,:)
real*8,allocatable,save::  Jx(:,:,:),Jy(:,:,:),Jz(:,:,:)
real*8,allocatable,save::  Hx(:,:,:),Hy(:,:,:),Hz(:,:,:)
real*8,allocatable,save::  Dx(:,:,:),Dy(:,:,:),Dz(:,:,:)
real*8,allocatable,save::  Bx(:,:,:),By(:,:,:),Bz(:,:,:)
real*8,allocatable,save::  pEx(:,:,:),pEy(:,:,:),pEz(:,:,:) ! /* Drude or Single pole Lorentz metal extra field components. */
real*8,allocatable,save::  pJx(:,:,:),pJy(:,:,:),pJz(:,:,:)
!/* Double pole Lorentz metal extra field components. */
real*8,allocatable,save::  Jx_0(:,:,:),Jy_0(:,:,:),Jz_0(:,:,:)
real*8,allocatable,save::  Jx_1(:,:,:),Jy_1(:,:,:),Jz_1(:,:,:)
real*8,allocatable,save::  Jx_2(:,:,:),Jy_2(:,:,:),Jz_2(:,:,:)
real*8,allocatable,save::  pJx_1(:,:,:),pJy_1(:,:,:),pJz_1(:,:,:)
real*8,allocatable,save::  pJx_2(:,:,:),pJy_2(:,:,:),pJz_2(:,:,:)

!Fourier Transform	
real*8,allocatable,save::  Ex_Re(:,:,:,:),Ey_Re(:,:,:,:),Ez_Re(:,:,:,:)
real*8,allocatable,save::  Hx_Re(:,:,:,:),Hy_Re(:,:,:,:),Hz_Re(:,:,:,:)
real*8,allocatable,save::  Ex_Im(:,:,:,:),Ey_Im(:,:,:,:),Ez_Im(:,:,:,:)
real*8,allocatable,save::  Hx_Im(:,:,:,:),Hy_Im(:,:,:,:),Hz_Im(:,:,:,:)
real*8,allocatable,save::  Ex_inc_Re(:),Ex_inc_Im(:)
!/* For TF-SF formulation */
real*8,allocatable,save::  Ex_inc(:),Hy_inc(:)
real*8,allocatable,save::  Ey_inc(:),Hx_inc(:)
!/* FDTD update coefficients */
real*8,allocatable,save::  aax(:),aay(:),aaz(:)
real*8,allocatable,save::  bbx(:),bby(:),bbz(:)
real*8,allocatable,save::  ccx(:),ccy(:),ccz(:)
real*8,allocatable,save::  ddx(:,:,:),ddy(:,:,:),ddz(:,:,:)
real*8,allocatable,save::  eex(:),eey(:),eez(:)
real*8,allocatable,save::  ffx(:),ffy(:),ffz(:)
real*8,allocatable,save::  ggx(:),ggy(:),ggz(:)
real*8,allocatable,save::  hhx(:),hhy(:),hhz(:)
real*8,allocatable,save::  iix(:),iiy(:),iiz(:)
real*8,allocatable,save::  jjx(:),jjy(:),jjz(:)
real*8,allocatable,save::  kkx(:),kky(:),kkz(:)
real*8,allocatable,save::  llx(:),lly(:),llz(:)
!/* for farfield parameter calculation */
real*8,allocatable,save::  Ex_cos(:,:,:),Ex_sin(:,:,:)
real*8,allocatable,save::  Ey_cos(:,:,:),Ey_sin(:,:,:)
real*8,allocatable,save::  Hx_cos(:,:,:),Hx_sin(:,:,:)
real*8,allocatable,save::  Hy_cos(:,:,:),Hy_sin(:,:,:)

!/// in input.f90 ///
real*8:: back_epsilon

!/// in output.f90 ///
real*8:: global_W
real*8::timestep
!integer,external::nz_multiple


		!///////////// Drude model /////////////////
	        !//                                       //
		!//                    (omega_p)^2        //
                !//  eps(w) = eps_b - ------------------  //  
                !//                   w^2 + j w gamma_0   //
                !///////////////////////////////////////////
     	     !/////////////////////// Lorentz  model ////////////////////
             !//                                                       //
	     !//                         depsr_0_1*(omega_0_1)^2     //
             !//  eps(w) = eps_b + -------------------------------     //  
             !//                    omega_0_1^2 - w^2 - j w 2 gamma_0_1//
             !///////////////////////////////////////////////////////////
  !/////////////////////////////////// Double Lorentz Model /////////////////////////////////////
  !/                                                                                            /
  !/                       depsr_0_1*(omega_0_1)^2               depsr_0_2*(omega_0_2)^2    /  
  !/  eps(w) = eps_b + --------------------------------- + ---------------------------------    /
  !/                   omega_0_1^2 - w^2 - j w  2 gamma_0_1   omega_0_2^2-w^2-j w 2 gamma_0_2   /
  !//////////////////////////////////////////////////////////////////////////////////////////////
type MetalObj		   
	character*15:: cshape     	
	real*8:: centeri   
	real*8:: centerj   
	real*8:: centerk   
	real*8:: centerl   
	real*8:: size1     
	real*8:: size2     
	real*8:: size3     
	!Normal metal
	real*8:: sigma_0
	!Drude
	real*8:: epsilon_b 
	real*8:: omega_p   
	real*8:: gamma_0   
	!Single Lorentz
	real*8:: depsr_0_1
	real*8:: omega_0_1 
	real*8:: gamma_0_1 
	!Double Lorentz
	real*8:: depsr_0_2
	real*8:: omega_0_2 
	real*8:: gamma_0_2 
end type MetalObj


type MolecularObj
character*15::cshape
real*8::centeri
real*8::centerj
real*8::centerk
real*8:: size1
real*8:: size2               
real*8:: size3
end type MolecularObj


end module variables
