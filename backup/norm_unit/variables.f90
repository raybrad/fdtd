module variables

!/// in fdtd.f90 ///
!module var_main

integer:: t
integer::del_t,dd,tx,tzero,Nf   !del_t: delta t;  dd:end time
integer:: TF_SF_pos(3,2),light_pos
real*8:: TF_SF_freq=0.0 !default
integer:: TF_SF_to
integer:: TF_SF_tdecay
real*8:: TF_SF_eps_bg
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
real,save:: pi,eo,uo,ups,light_speed,bohr,el_charge
integer:: misize, mjsize, mksize
integer:: pisize, pjsize, pksize
integer:: cisize, cjsize, cksize
integer:: metal_type=0,metal_t       !default       
real*8::bg_eps,eps_b,omg_p,gam_0,sigm_0 
character*10::csp,Ef,Hf,sourcetyp  !mat_file
real*8::pos,freq,wav_len,lambda_min,lambda_max,df,fmin,fmax     !!wc=a/wavelength ,normalized frequency of a resonance in question
integer::interval
real*8::amplitude
!/// in memory.f90/
!module var_memo

integer,allocatable,save::cposition(:,:,:)
real,allocatable,save::  epsilonx(:,:,:),epsilony(:,:,:),epsilonz(:,:,:)
!/* Input dielectric structure for Drude or single pole Lorentz. */
real,allocatable,save::  Lepsilon(:,:,:),Lomega(:,:,:),Lgamma(:,:,:),Lomega0(:,:,:),Lsigma(:,:,:),Molecular(:,:,:)
!/* field components */
real,allocatable,save::  Ex(:,:,:),Ey(:,:,:),Ez(:,:,:)
real,allocatable,save::  Jx(:,:,:),Jy(:,:,:),Jz(:,:,:)
real,allocatable,save::  Hx(:,:,:),Hy(:,:,:),Hz(:,:,:)
real,allocatable,save::  Dx(:,:,:),Dy(:,:,:),Dz(:,:,:)
real,allocatable,save::  Bx(:,:,:),By(:,:,:),Bz(:,:,:)

!Fourier Transform	
real,allocatable,save::  Ex_Rw(:,:,:,:),Ey_Rw(:,:,:,:),Ez_Rw(:,:,:,:)
real,allocatable,save::  Hx_Rw(:,:,:,:),Hy_Rw(:,:,:,:),Hz_Rw(:,:,:,:)
!/* For TF-SF formulation */
real,allocatable,save::  Ex_inc(:),Hy_inc(:)
real,allocatable,save::  Ey_inc(:),Hx_inc(:)
!/* FDTD update coefficients */
real,allocatable,save::  aax(:),aay(:),aaz(:)
real,allocatable,save::  bbx(:),bby(:),bbz(:)
real,allocatable,save::  ccx(:),ccy(:),ccz(:)
real,allocatable,save::  ddx(:,:,:),ddy(:,:,:),ddz(:,:,:)
real,allocatable,save::  eex(:),eey(:),eez(:)
real,allocatable,save::  ffx(:),ffy(:),ffz(:)
real,allocatable,save::  ggx(:),ggy(:),ggz(:)
real,allocatable,save::  hhx(:),hhy(:),hhz(:)
real,allocatable,save::  iix(:),iiy(:),iiz(:)
real,allocatable,save::  jjx(:),jjy(:),jjz(:)
real,allocatable,save::  kkx(:),kky(:),kkz(:)
real,allocatable,save::  llx(:),lly(:),llz(:)

!/// in input.f90 ///
real*8:: back_epsilon

!/// in output.f90 ///
real*8:: global_W
real*8::timestep
!integer,external::nz_multiple

!!!!!!!!!!!!!type!!!!!!!
type  obj 

	character*15:: cshape
	real*8:: centeri
	real*8:: centerj
	real*8:: centerk
	real*8:: size1
	real*8:: size2
	real*8:: size3
	real*8:: epsob
!	real,allocatable:: matrix(:,:) !// contour matrix data
!	integer::   col !// matrix col num
!	integer::   row !// matrix row num
end type

!///////////// Drude model /////////////////
!//                                       //
!//                    (omega_p)^2        //
!//  eps(w) = eps_b - ------------------  //  Note: The former mobj is now generalized to Lobj.
!//                   w^2 + j w gamma_0   //
!///////////////////////////////////////////
type Lobj 

	character*15:: cshape
	real*8:: centeri
	real*8:: centerj
	real*8:: centerk
	real*8:: centerl
	real*8:: size1     !/////////////////////// Lorentz  model /////////////////
	real*8:: size2     !//                                                    //
	real*8:: size3     !//                                 (omega_p)^2        //
	real*8:: epsilon_b !//  eps(w) = eps_b + -------------------------------  //  
	real*8:: omega_p   !//                    omega_0^2 - w^2 - j w gamma_0   //
	real*8:: gamma_0   !////////////////////////////////////////////////////////
!	real*8:: omega_0   != eps_b- omega_p^2/gamma/(j w) +omega_^2/gamma/(j w +gamma)
	real*8::sigma_0
!        real,allocatable:: matrix(:,:)  !// contour matrix data
!	integer::   col !// matrix col num
!	integer::   row !// matrix row num
end type


	     !//////////////////////////////////// Double Lorentz Model //////////////////////////////////////
	     !//                                                                                            //
	     !//                             (omega_p_1)^2                 (-1)^S*(omega_p_2)^2             //      
	     !//  eps(w) = eps_b + --------------------------------- + ---------------------------------    //   
	     !//                   omega_0_1^2 - w^2 - j w gamma_0_1   omega_0_2^2 - w^2 - j w gamma_0_2    //
	     !////////////////////////////////////////////////////////////////////////////////////////////////

type Mobj
character*15::cshape
real*8::centeri
real*8::centerj
real*8::centerk
real*8:: size1
real*8:: size2
real*8:: size3
!real*8::centerl
end type



end module variables
