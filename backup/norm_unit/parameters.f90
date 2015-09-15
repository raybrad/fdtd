module parameters
use variables 
implicit none


contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!Read Input!!
subroutine readinput(hFile)
integer,intent(in)::hFile
real*8::temp1,temp2,temp3
namelist /box/ xsize,ysize,zsize
namelist /lattice/ lattice_x,lattice_y,lattice_z,l_n
namelist /pml/ pmlil,pmlir,pmljl,pmljr,pmlkl,pmlkr
namelist /epsval/ bg_eps,eps_b,omg_p,gam_0,sigm_0,metal_t
namelist /sourcefield/ Ef, Hf,sourcetyp,pos,wav_len,del_t,dd,tx,tzero,lambda_min,lambda_max,triple,Nf,df,dofft
namelist /out_res/ interval,amplitude
namelist /qm_em/ lqmtoem,lemtoqm

!--Default values-----------------------------------

S_factor=2.0
xsize=1; ysize=1; zsize=1
lattice_x=100; lattice_y=100; lattice_z=100;
l_n=100E-9 
pmlil=10; pmlir=10
pmljl=10; pmljr=10
pmlkl=10; pmlkr=10 
bg_eps=1.d0

metal_t=2	!drude metal
eps_b= 0.d0;omg_p=0.0;gam_0=0.0
sigm_0=0.0

interval = 200
amplitude= 1.d0

sourcetyp=''
triple=.false.
tzero=0
tx=0

lqmtoem=.false.
lemtoqm=.false.

dofft=.true.
fmin=0.d0;fmax=0.d0
Nf=100
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

!!!define sourcefiled parameter:norm freq,decay time,end time!!!!!!!
!!!E direction,H diretion,incident position!!!!!!!!!!!!!!!!!!!!!!!!
rewind hFile
read(hFile,sourcefield,end=60)
60 continue
if (sourcetyp(1:8)=='Gaussian') then
temp1=0.5d0*(l_n/lambda_min+l_n/lambda_max)
temp2=0.5d0*(l_n/lambda_min-l_n/lambda_max)
wav_len=l_n/temp1
del_t=S_factor*lattice_x/temp2/6.d0
tzero=3*del_t
!dd= 6*del_t+5000\
else
	lambda_max=wav_len+50.d-9
	lambda_min=wav_len-50.d-9
endif

fmin=l_n/(lambda_max+50.d-9)
fmax=l_n/(lambda_min-50.d-9)
df=(fmax-fmin)/(Nf-1)


rewind hFile
read(hFile,out_res,end=70)
70 continue


rewind hFile
read(hFile,qm_em,end=80)
80 continue


!normalized frequency	  100nm/400nm
freq = l_n/wav_len							 
TF_SF_freq=freq
TF_SF_to=tzero
TF_SF_tdecay=del_t
TF_SF_eps_bg=bg_eps

end subroutine readinput
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine parameter_check


write(6,*)'box xsize,ysize,zsize:',xsize,ysize,zsize
write(6,*)'latticex,latticey,latticez: ',lattice_x,lattice_y,lattice_z
write(6,*)'fdtd unit length l_n(m):',l_n
write(6,*)'pml(il,lr,jl,jr,kl,kr): '
write(6,*) pmlil,pmlir,pmljl,pmljr,pmlkl,pmlkr
write(6,*)'backgroud eps:',bg_eps
write(6,*)'metal parameters:',eps_b,omg_p,gam_0,sigm_0,metal_t
write(6,*)'inc wavelegth(m)',wav_len
write(6,*)'source amplitude (V/m):', amplitude							   
write(6,*)'light info t0,tdcay,dd:',tzero,del_t,dd
write(6,*)'light info, fmin,fmax,Nf,df',fmin,fmax,Nf,df
end subroutine parameter_check	

subroutine source_check
	
write(6,*) "dt = ",timestep,"(fs)"
write(6,*) "simulation time length:",(dble(dd)*timestep),'(fs)'
write(6,*) "max point of gaussian pulse:",(dble(3*del_t)*timestep),'(fs)'
write(6,*) "pulse frequency:",freq/abs(l_n)*light_speed*1.d-15,'(fs^-1)'						
write(6,*) "pulse period:",abs(l_n)/(freq*light_speed)*1.d15,'(fs^-1)'						
end subroutine source_check
!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!///

subroutine set_default_parameter(S)
real*8:: S 
integer:: lcm_temp 
integer:: i, j 

       pi=3.141592 
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
	
	isize=floor(0.5+xsize*(lattice_x-1))+1 
	jsize=floor(0.5+ysize*(lattice_y-1))+1
	ksize=floor(0.5+zsize*(lattice_z-1))+1
	print *,'isize',isize
	print *,'jsize',jsize
	print *,'ksize',ksize
	!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//

	xcenter=xsize/2.d0 
	ycenter=ysize/2.d0 
	zcenter=zsize/2.d0   

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
	else
		write(6,*) 'Currently undealt'
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
