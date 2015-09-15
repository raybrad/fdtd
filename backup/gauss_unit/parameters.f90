module parameters
use variables 
implicit none


contains
!----------------------------------------------------
subroutine readinput

namelist /box/ boxsize,lattice_x,lattice_y,lattice_z
namelist /pml/ pmlil,pmlir,pmljl,pmljr,pmlkl,pmlkr
namelist /mat_para/ bg_eps,eps_b,omg_p,gam_0
namelist /source/wav_len,tdecay
namelist /qm_em/ lqmtoem,lemtoqm

!--default parameters
S=2.0						!Stability Couriant consant(S>sqrt(1/(dx)^2+1/(dy)^2+1/(dz)^2)
						!Thus dt=dx/(2*c) E~0.5H  H~0.5E
pi=3.141592 
eo=8.854e-12 					!vacuum epsilon0 (C/V/m)  1/miu0/c^2
uo=pi*4.0e-7 
ups=1.0 
light_speed=1.0/sqrt(eo*uo) 

boxsize=100.0d-9				!Size of simulation box 						
isize=100;jsize=100;ksize=100			!Grid number of x,y,z	

pmlil=10; pmlir=10				!Perfect Matched Layers parameters
pmljl=10; pmljr=10
pmlkl=10; pmlkr=10
orderxl = 3.5; orderyl = 3.5;  orderzl = 3.5 
orderxr = 3.5; orderyr = 3.5;  orderzr = 3.5 
sig_axl = 1.0; sig_ayl = 1.0;  sig_azl = 1.0 
sig_axr = 1.0; sig_ayr = 1.0;  sig_azr = 1.0 

bg_eps=1.0					!background permittivity
eps_b=1.d0;omg_p=23.89d15;gam_0=1.07162d15	!Drude model parameters (default Al)

wav_len=400.d-9					!source

lqmtoem=.false.
lemtoqm=.false.
!--read input
rewind 5
read(5,box,end=10)
10 continue

rewind 5
read(5,pml,end=20)
20 continue

rewind 5
read(5,mat_para,end=30)
30 continue

rewind 5
read(5,source,end=50)
50 continue

rewind 5
read(5,qm_em,end=60)
60 continue

!--
xcenter=0.5
ycenter=0.5
zcenter=0.5


!dx dt
dl=boxsize/(isize-1)
dt=dl/(2.d0*light_speed)

coeff=1/S						!update coeff
!normalized frequency	  100nm/400nm
freq = l_n/wav_len							 
TF_SF_freq=freq
TF_SF_to=tzero
TF_SF_tdecay=tdecay
TF_SF_eps_bg=bg_eps

write(6,*)'box size,isize,jsize,ksize:',isize,jsize,ksize
write(6,*)'pml(il,lr,jl,jr,kl,kr): '
write(6,*) pmlil,pmlir,pmljl,pmljr,pmlkl,pmlkr
write(6,*)'backgroud eps:',bg_eps
write(6,*)'metal parameters:',eps_b,omg_p,gam_0,sigm_0,metal_t

write(6,*) "parameters...ok\n" 

end subroutine readinput

end module parameters
