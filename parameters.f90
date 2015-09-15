module parameters
use variables 
!use IFPORT
implicit none


contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!Read Input!!
subroutine readinput(hFile)
integer,intent(in)::hFile
namelist /box/ xsize,ysize,zsize
namelist /lattice/ lattice_x,lattice_y,lattice_z,l_n
namelist /pml/ pmlil,pmlir,pmljl,pmljr,pmlkl,pmlkr
namelist /epsval/ bg_eps,metal_t,Mole_type,sigm_0,eps_b,omg_p,gam_0,depsr_0_1,omega_0_1,gamma_0_1,depsr_0_2,omega_0_2, gamma_0_2
namelist /input_obj/ metal_shape,centerx1,centery1,centerz1,centerxl1,size1,size2,size3,centerx2,centery2,centerz2,centerxl2,radiusSphere
namelist /sourcefield/ Ef, Hf,sourcetyp,Gauss_typ,pos,wav_len,del_t,dd,tx,tzero,lambda_min,lambda_max,triple,Nf,df,dofft,Real_tzero,Real_deltat,Real_totalt
namelist /out_res/ ocomp,oplane,ovalue,oname,interval,amplitude,measureLambda
namelist /qm_em/ lqmtoem,lemtoqm,separation
!--Default values-----------------------------------

!S_factor=2.0
S_factor=sqrt(3.d0)/0.95d0		!1.8232
!S_factor=1.8444			!to let dt=4.57*10^-4 fs
xsize=1; ysize=1; zsize=1
lattice_x=100; lattice_y=100; lattice_z=100;
l_n=100E-9 
pmlil=10; pmlir=10
pmljl=10; pmljr=10
pmlkl=10; pmlkr=10 
bg_eps=1.d0

metal_t=2	!drude metal
Mole_type=0
eps_b= 0.d0;omg_p=0.0;gam_0=0.0
sigm_0=0.0
depsr_0_1=0.d0;gamma_0_1=0.d0;omega_0_1=0.d0
depsr_0_2=0.d0;gamma_0_2=0.d0;omega_0_2=0.d0

metal_shape="sphere"
radiusSphere=10.d-9		!10nm
centery1=0.0
centerz1=0.0
centerxl1=0.0
centery2=0.0
centerz2=0.0
centerxl2=0.0

interval = 200
amplitude= 1.d0

sourcetyp=''
triple=.false.
tzero=0
tx=0

lqmtoem=.false.
lemtoqm=.false.
separation=0.d0

dofft=.false.
measureLambda=600.d-9
fmin=0.d0;fmax=0.d0
lambda_min=400d-9
lambda_max=700d-9
Nf=20
!!!define box structure size(xsize,ysize,zsize)!!!!!!!!!!!!!
rewind hFile
read(hFile,box,end=10)
10 continue

!!!!define lattice size(how many grid in each size)!!!!!!!!!
rewind hFile
read(hFile,lattice,end=20)
20 continue

!!!!define PML size: i,j,k (left,right)!!!!!!!!!!!!!!!!!!!!!

rewind hFile
read(hFile,pml,end=30)
30 continue
!!define background eps and metal eps,omega,gamma!!!!!!!!!!
rewind hFile
read(hFile,epsval,end=40)
40 continue
!!!!!!!define input (metal) shape,position!!!!!!!!!!!!!!!!!
rewind hFile
read(hFile,input_obj,end=50)
50 continue

!!!define sourcefiled parameter:norm freq,decay time,end time!!!!!!!
!!!E direction,H diretion,incident position!!!!!!!!!!!!!!!!!!!!!!!!
rewind hFile
read(hFile,sourcefield,end=60)
60 continue

rewind hFile
read(hFile,out_res,end=70)
70 continue
amplitude=amplitude/(5.291772083d-11)			!(V/bohr -> V/m)

rewind hFile
read(hFile,qm_em,end=80)
80 continue
end subroutine readinput

subroutine calWavePara
real*8::tempdt
real*8::temp1,temp2
integer::m
!if (sourcetyp(1:8)=='Gaussian') then
!temp1=0.5d0*(l_n/lambda_min+l_n/lambda_max)
!temp2=0.5d0*(l_n/lambda_min-l_n/lambda_max)
!wav_len=l_n/temp1
!del_t=S_factor*lattice_x/temp2/6.d0
!tzero=3*del_t
!dd= 6*del_t+5000\
!else
!	lambda_max=wav_len+50.d-9
!	lambda_min=wav_len-50.d-9
!endif

fmin=l_n/(lambda_max)
fmax=l_n/(lambda_min)
df=(fmax-fmin)/(Nf-1)

m=int((l_n/measureLambda-fmin)/df)+1
write(6,*) 'measureLambda position',m
if (sourcetyp(1:8)=='Gaussian' .or. sourcetyp(1:9)=='PGaussian') then
tempdt= l_n/(1d-9)/(light_speed*S_factor*ds_x*lattice_x)*1d6	!dt (fs)
tzero= int(Real_tzero/tempdt)
del_t= int(Real_deltat/tempdt)
dd= int(Real_totalt/tempdt)
endif



!normalized frequency	  100nm/400nm
freq = l_n/wav_len							 
TF_SF_freq=freq
TF_SF_to=tzero
TF_SF_tdecay=del_t
TF_SF_eps_bg=bg_eps

write(6,*) 'TF_SF_freq',TF_SF_freq
write(6,*) 'TF_SF_to',TF_SF_to
write(6,*) 'TF_SF_tdecay',TF_SF_tdecay
end subroutine calWavePara
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine parameter_check


write(6,*)'box xsize,ysize,zsize:',xsize,ysize,zsize
write(6,*)'latticex,latticey,latticez: ',lattice_x,lattice_y,lattice_z
write(6,*)'fdtd unit length l_n(m):',l_n
write(6,*)'pml(il,lr,jl,jr,kl,kr): '
write(6,*) pmlil,pmlir,pmljl,pmljr,pmlkl,pmlkr
write(6,*)'backgroud eps:',bg_eps
write(6,*)'metal parameters:',eps_b,omg_p,gam_0,sigm_0,metal_t
write(6,*) 'Mole_type',Mole_type
write(6,*)'inc wavelegth(m)',wav_len
write(6,*)'source amplitude (V/m):', amplitude							   
write(6,*) 'dofft:',dofft
write(6,*) 'metal shape',metal_shape
write(6,*)'separation:',separation
write(6,*)'light info t0,tdcay,dd:',tzero,del_t,dd
write(6,*)'light info, fmin,fmax,Nf,df',fmin,fmax,Nf,df
write(6,*) 'measureLambda',measureLambda
end subroutine parameter_check	

subroutine source_check
	
write(6,*) "dt = ",timestep,"(fs)"
write(6,*) "simulation time length:",(dble(dd)*timestep),'(fs)'
write(6,*) "max point of gaussian pulse:",(dble(3*del_t)*timestep),'(fs)'
write(6,*) "pulse frequency:",freq/abs(l_n)*light_speed*1.d-15,'(fs^-1)'						
write(6,*) "pulse period:",abs(l_n)/(freq*light_speed)*1.d15,'(fs^-1)'						
end subroutine source_check
!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!///
!
!
subroutine pml_size( il,ir,jl,jr,kl,kr)
integer::il,ir,jl,jr,kl,kr

        pmlil=il
        pmlir=ir
        pmljl=jl
        pmljr=jr
        pmlkl=kl
        pmlkr=kr
	if(il==0)   pmlil=0
	if(ir==0)   pmlir=0
	if(jl==0)   pmljl=0
	if(jr==0)   pmljr=0
	if(kl==0)   pmlkl=0
	if(kr==0)   pmlkr=0
end subroutine pml_size

subroutine set_default_parameter(S)
real*8:: S 
integer:: lcm_temp 
integer:: i, j 

       planckev=4.13566733d-15 !;//eV s
       pi=3.14159265358979323846264338327950288419716939937510 !
       eo=8.854e-12 !vacuum epsilon0 (C/V/m)  1/miu0/c^2
       uo=pi*4.0e-7 
       ups=1.0 
       light_speed=1.0/sqrt(eo*uo) 
       bohr=5.291772083d-11
       el_charge=1.602176462d-19
       write(6,*) 'set default parameter S:',S

	lcm_temp = lcm3(lattice_x, lattice_y, lattice_z) 

	ds_x = lcm_temp/lattice_x   !1
	ds_y = lcm_temp/lattice_y
	ds_z = lcm_temp/lattice_z
	!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!///
	!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!/// 
	
	S_factor = S 

	dt=1/light_speed/S_factor  !// S=Couriant parameter for numerical stability
				   !in Real space dt=dx/c/S, here we treat ds_x as 1, then dt=1/c/S, dt/ds_x=1/c/S
	orderxl = 3.5; orderyl = 3.5;  orderzl = 3.5 
	orderxr = 3.5; orderyr = 3.5;  orderzr = 3.5 

	sig_axl = 0.8; sig_ayl = 0.8;  sig_azl = 0.8 
	sig_axr = 0.8; sig_ayr = 0.8;  sig_azr = 0.8 

	kx = 1.0;  ky = 1.0;  kz =1.0 
	
	isize=floor(0.5+(xsize*lattice_x)) 
	jsize=floor(0.5+(ysize*lattice_y))
	ksize=floor(0.5+(zsize*lattice_z))
	print *,'isize',isize
	print *,'jsize',jsize
	print *,'ksize',ksize
	!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//

	xcenter=xsize/2 
	ycenter=ysize/2 
	zcenter=zsize/2   

	misize=isize 
	pisize=isize-1 
	cisize=isize/2  
	mjsize=jsize 
	pjsize=jsize-1 
	cjsize=jsize/2
	mksize=ksize 
	pksize=ksize-1 
	cksize=ksize/2  

	write(6,*)'pisize',pisize
	write(6,*)'pjsize',pjsize
	!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!///
	print *,"parameters...ok\n" 
end subroutine set_default_parameter 

subroutine set_sigma_order(oxl, oxr, oyl, oyr,ozl,ozr)
real*8::oxl,oxr,oyl,oyr,ozl,ozr
 
	orderxl = oxl; orderyl = oyl; orderzl = ozl   !//default value = 3.5
	orderxr = oxr; orderyr = oyr; orderzr = ozr   !//default value = 3.5
end subroutine set_sigma_order

subroutine set_sigma_max(axl, axr, ayl, ayr, azl, azr)
real*8::axl,axr,ayl,ayr,azl,azr
 
	sig_axl = axl; sig_ayl = ayl; sig_azl = azl;   !//default value = 1.0
	sig_axr = axr; sig_ayr = ayr; sig_azr = azr    !//default value = 1.0
end subroutine set_sigma_max

subroutine set_kappa(kappa_x,kappa_y,kappa_z)
real*8::kappa_x,kappa_y,kappa_z
 
	kx = kappa_x   !//defaut value = 1.0
	ky = kappa_y   !//defaut value = 1.0
	kz = kappa_z   !//defaut value = 1.0
end subroutine set_kappa


subroutine material_complexity(input)
 integer::input
!	/* 0: No Metal (default).
!	   1: Normal Metal.
!	   2: Drude Metal.

	metal_type = input 
	if (metal_type==0) then
		write(6,*) 'No metal'
	elseif(metal_type==1) then
		write(6,*) 'Normal metal'
	elseif(metal_type==2) then
		write(6,*) 'Drude  metal'
	elseif(metal_type==3) then
		write(6,*) 'Single Lorentz  metal'
	elseif(metal_type==4) then
		write(6,*) 'Double Lorentz  metal'
	elseif(metal_type==5) then
		write(6,*) 'Drude+Double Lorentz  metal'
	else
		write(6,*) 'metal_type > 5 Currently undealt'
		stop
	endif

end subroutine material_complexity
!
function lcm3(a,b,c)
integer::lcm3
integer::a,b,c
	lcm3= lcm2(a, lcm2(b,c))  
end function lcm3
!
function lcm2(m,n) !least commom multiple 
integer::lcm2
integer::m,n
	lcm2=( m*n/gcd2(m,n) )   !//use the relation   gcd(m,n)*lcm(m,n)=m*n
end function lcm2
!
function gcd2(m,n) !greatest common divisor
integer::gcd2
integer::m,n
integer:: temp 
integer:: r  !//remainder

	if(n>m) then !//Make m>n
	 
		temp = m 
		m = n 
		n = temp 
	endif
    do while(.true.)
      r=mod(m,n)
      if(r==0)exit
      m=n
      n=r
    enddo!//Euclid algorithm 
	gcd2=n 
end function gcd2

end module parameters
