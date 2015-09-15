module output
use parameters
implicit none

contains
subroutine real_space_param(a_m, w_n)
real*8::a_m,a_nm,w_n

	!// a_nm is in the unit of 'nanometer' l_n
	!// w_n is the normalized frequency

	!/// Set global normalized frequency!//!//!//!//!//!//!// 
	global_W = w_n 
	!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!// (ver. 8.45)
        a_nm=a_m/(1d-9)
	write(6,*)'fdtd unit length a_m:',a_m,'(m)',a_nm,'(nm)'
	write(6,*)'normalized frequency global_W(a/lambda)',global_W
	write(6,*)'incident wavelength',a_nm/global_W,'(nm)'

	open(unit=16,file="Real_Space_Param.dat",status='replace',form='formatted') 

	write(16,*) "================================= " 	
	write(16,*) "The FDTD calculation domain  " 
	write(16,*) "================================= " 	
	write(16,*) "Lx =",xsize,"(a)"
	write(16,*) "Ly =",ysize,"(a)" 
	write(16,*) "Lz =",zsize,"(a)" 
	write(16,*) "--------------------------------- " 
	write(16,*) "Lx =",xsize*a_nm,"(nm)" 
	write(16,*) "Ly =",ysize*a_nm,"(nm)"
	write(16,*) "Lz =",zsize*a_nm,"(nm)" 
	write(16,*) " " 

	write(16,*) "================================= " 	
	write(16,*) "The PML size  " 
	write(16,*) "================================= " 	
	write(16,*) "pmlil =",pmlil,"(u)"," pmlir =",pmlir,"(u)" 
	write(16,*) "pmljl =",pmljl,"(u)"," pmljr =",pmljr,"(u)" 
	write(16,*) "pmlkl =",pmlkl,"(u)"," pmlkr =",pmlkr,"(u)" 
	write(16,*) " " 

	write(16,*) "================================= " 	
	write(16,*) "The speed of light  " 
	write(16,*) "================================= " 	
	write(16,*) "c =",light_speed,"(m/s)" 
	write(16,*) " " 

	write(16,*) "================================= " 	
	write(16,*) "The Stability parameter 'S'  " 
	write(16,*) "================================= " 	
	write(16,*) "S =",S_factor,"(1/q)" 
	write(16,*) "S =",S_factor*ds_x*lattice_x/(a_nm*1E-9),"(1/m)" 
	write(16,*) " " 

	write(16,*) "================================= " 	
	write(16,*) "The size of the FDTD grid  " 
	write(16,*) "================================= " 	
	write(16,*) "dx =",ds_x,"(q)" 
	write(16,*) "dy =",ds_y,"(q)"
	write(16,*) "dz =",ds_z,"(q)" 
	write(16,*) "--------------------------------- " 
	write(16,*) "dx =  ",a_nm/lattice_x,"(nm)" 
	write(16,*) "dy =  ",a_nm/lattice_y,"(nm)"
	write(16,*) "dz =  ",a_nm/lattice_z ,"(nm)"
	write(16,*) " " 

	write(16,*) "================================= " 	
	write(16,*) "The size of the FDTD time step  " 
	write(16,*) "================================= " 
	timestep=(a_nm)/(light_speed*S_factor*ds_x*lattice_x)*1E6
	write(16,*) "dt = ",timestep*1E-15,"(sec)" 
	write(16,*) " " 

	write(16,*) "================================= " 	
	write(16,*) "The dipole wavelength in vacuum (in FDTD) " 
	write(16,*) "================================= " 	
	write(16,*) "lambda =  ",lattice_x/w_n ,"(u)"
	write(16,*) " " 

	write(16,*) "================================= " 	
	write(16,*) "The dipole frequency (in FDTD) " 
	write(16,*) "================================= " 	
	write(16,*) "f =  ",w_n/(S_factor*ds_x*lattice_x),"(1/update)" 
	write(16,*) "T = ",(S_factor*ds_x*lattice_x)/w_n,"(update)" 
	write(16,*) " " 

	write(16,*) "================================= " 	
	write(16,*) "The Origin FFT correction factor  " 
	write(16,*) "================================= " 	
	write(16,*) "x ",S_factor*ds_x*lattice_x 
	write(16,*) " " 

	close(16) 
end subroutine real_space_param


subroutine out_epsilon(plane,value,cname)
character::plane
character*10::cname
real*8::value
integer:: i,j,k,iu 
integer:: i_range, j_range, k_range  
integer:: mz  !// index for non_uniform_grid multiplication 

	!// for Hz_parity, normal cases
	i_range = isize 
	j_range = jsize 
	k_range = ksize  
	
	iu=17

	open(unit=iu,file=cname,status='replace',form='FORMATTED') 

	if(plane=="x") then
		i=floor(0.5+(value+xcenter)*(lattice_x-1))+1 
		do k=1,k_range
				do j=1,j_range
					if(metal_type==0) then
					        write(iu,'(F12.4)',advance='no') eps(i,j,k)

	                                endif
					if(metal_type==1 .or. metal_type==2 ) then
					       if(meps(i,j,k)==0.0) then
							write(iu,'(F12.4)',advance='no') eps(i,j,k) 
						else
						     
							write(iu,'(F12.4)',advance='no') meps(i,j,k) 
						endif
					endif
				enddo
				write(iu,*)  
		enddo
	endif
	if(plane=="y") then
	 
		j=floor(0.5+(value+ycenter)*(lattice_y-1))+1 
		do k=1,k_range
				do i=1,i_range
				 
					if(metal_type==0) then
						write(iu,'(F12.4)',advance='no') eps(i,j,k) 
	                                endif
					if(metal_type==1 .or. metal_type==2 ) then
					 
						if(meps(i,j,k)==0.0) then
							write(iu,'(F12.4)',advance='no') eps(i,j,k) 
						else
							write(iu,'(F12.4)',advance='no') meps(i,j,k) 
	                                        endif
					endif
				enddo
				write(iu,*)  
		enddo
	endif

	if(plane=="z") then
	 
		k=floor(0.5+(value+zcenter)*(lattice_z-1))+1 
		do j=1,j_range
			do i=1,i_range
			 
					if(metal_type==0) then
						write(iu,'(F12.4)',advance='no') eps(i,j,k) 
					endif
					if(metal_type==1 .or. metal_type==2 ) then
					 
						if(meps(i,j,k)==0.0) then
							write(iu,'(F12.4)',advance='no') eps(i,j,k) 
						else
							write(iu,'(F12.4)',advance='no') meps(i,j,k) 
						endif
					endif
			enddo
			write(iu,*)  
		enddo
	endif
	close(iu) 
	print *,"out_epsilon...ok " 
end subroutine out_epsilon


subroutine out_plane(component,plane,value,lastname)
 
        character(len=*)::plane
	real*8::value
	character(len=*)::component,lastname
	character*20::cname
	integer:: i,j,k,iu 
	integer:: i_range, j_range, k_range  

	i_range = isize 
	j_range = jsize 
	k_range = ksize  

	iu=18
	
	write(cname,'(I7.7)') t
         cname=trim(cname)//trim(lastname) 
	open(unit=iu,file=cname,status='replace',form='formatted') 

	if(plane=="x") then
		i=floor(0.5+(value+xcenter)*(lattice_x-1))+1
		do k=1,k_range
				do j=1,j_range	
                                write(iu,'(E12.5,3x)',advance='no') grid_value(component,i,j,k) 
	                        enddo
				write(iu,*)  
		enddo
	endif

	if(plane=="y") then
	 
		j=floor(0.5+(value+ycenter)*(lattice_y-1))+1 
		do k=1,k_range
				do i=1,i_range
                                write(iu,'(E12.5,3x)',advance='no') grid_value(component,i,j,k)  
	                        !print *,'grid_value',grid_value(component,i,j,k)
	                        enddo
				write(iu,*) 
		enddo
	endif

	if(plane=="z") then
	 
		k=floor(0.5+(value+zcenter)*(lattice_z-1))+1 
		do j=1,j_range
			do i=1,i_range
                        write(iu,'(E12.5,3x)',advance='no') grid_value(component,i,j,k)  	
	                enddo
			write(iu,*)  
		enddo
	endif
	close(iu) 
	print *,"out",' ',component," ...ok " 
end subroutine out_plane

subroutine out_plane_cut(component,plane,value,lastname)
 
        character(len=*)::plane
	real*8::value
	character(len=*)::component,lastname
	character*20::cname
	integer:: i,j,k,iu 
	integer:: i_range, j_range, k_range  
	integer:: mz  !//for non_uniform_grid multiplication

	i_range = isize-pmlir 
	j_range = jsize-pmljr 
	k_range = ksize-pmlkr  

	iu=18
	
	write(cname,'(I7.7)') t
         cname=trim(cname)//trim(lastname)

	write(6,*) 'i j k range',i_range,j_range,k_range
	open(unit=iu,file=cname,status='replace',form='formatted') 

	if(plane=="x") then
		i=floor(0.5+(value+xcenter)*(lattice_x-1))+1
		do k=1+pmlkl,k_range
				do j=1+pmljl,j_range	
                                write(iu,'(E12.5,3x)',advance='no') grid_value(component,i,j,k) 
	                        enddo
				write(iu,*)  
		enddo
	endif

	if(plane=="y") then
	 
		j=floor(0.5+(value+ycenter)*(lattice_y-1))+1 
		do k=1+pmlkl,k_range
				do i=1+pmlil,i_range
                                write(iu,'(E12.5,3x)',advance='no') grid_value(component,i,j,k)  
	                        !print *,'grid_value',grid_value(component,i,j,k)
	                        enddo
				write(iu,*) 
		enddo
	endif

	if(plane=="z") then
	 
		k=floor(0.5+(value+zcenter)*(lattice_z-1))+1 
		do j=1+pmljl,j_range
			do i=1+pmlil,i_range
                        write(iu,'(E12.5,3x)',advance='no') grid_value(component,i,j,k) 	
	                enddo
			write(iu,*)  
		enddo
	endif
	close(iu) 
	print *,"out",' ',component," ...ok " 
end subroutine out_plane_cut



subroutine out_point(component,x,y,z,ti,tf,cname)
 
        character(len=*)::component
	character(len=*)::cname
        real*8::x,y,z
	integer::ti,tf
	integer:: i,j,k 

	if(ti<=t .and. t<=tf) then
       !  print *,'t ',t	
	      
		i=floor(0.5+(x+xcenter)*(lattice_x-1))+1 
		j=floor(0.5+(y+ycenter)*(lattice_y-1))+1 
		k=floor(0.5+(z+zcenter)*(lattice_z-1))+1 
		if(t==ti) then
                   open(18,file=cname,status='REPLACE')
	           close(18)
		endif
		open(unit=18,file=cname,position='APPEND') 
		write(18,'(2(E12.5,3x))') dble(t)*timestep, grid_value(component,i,j,k)    !if it's E^2 should be grid_value ^2
		close(18) 
     endif
end subroutine out_point

subroutine lodestar_qmcal(x,y,z)
! input Efield at molecular pos in  lodestar to get current at interface
real*8::x,y,z
integer::i,j,k,istat,system
i=floor(0.5+(x+xcenter)*(lattice_x-1))+1
j=floor(0.5+(y+ycenter)*(lattice_y-1))+1
k=floor(0.5+(z+zcenter)*(lattice_z-1))+1 
open(159,file='fdembound',status='REPLACE',form='formatted') !v/m *0.529*1d-10 =>eV/bohr
write(159,'(5(E20.5))')dble(t)*timestep,timestep,grid_value('Ex',i,j,k) ,grid_value('Ey',i,j,k) ,grid_value('Ez',i,j,k) 
close(159)

open(160,file='Embound.dat',position='append',form='formatted') !v/m *0.529*1d-10 =>eV/bohr
write(160,'(5(E20.5))')dble(t)*timestep,timestep,grid_value('Ex',i,j,k) ,grid_value('Ey',i,j,k) ,grid_value('Ez',i,j,k) 
close(160)
!!lodestar
istat=system('cp -f fdembound ./QMwork; cd ./QMwork;/home/hy/lodestar/bin/lodestar_half < in.td  > out.td && cp -f fdqmbound ../')
	if (istat==0) then
	  print *,'success in qm_td calculation'
	  continue
	else
	  print *,'error in in.td calculation'
	  stop
	endif


end subroutine lodestar_qmcal


function  grid_value(component,i,j,k)
character(len=*)::component
integer::i,j,k
real*8::grid_value
!           print *,component
        select case(component)
	case("Ex") 
        grid_value= Ex(i,j,k) 
	case("Ey")
	grid_value= Ey(i,j,k)
	case("Ez")
	grid_value= Ez(i,j,k)
        case("Hx")
	grid_value= Hx(i,j,k)
	case("Hy")
	grid_value= Hy(i,j,k)
	case("Hz")
	grid_value= Hz(i,j,k)
	case("Jx")
	grid_value= Jx(i,j,k)
	case("Jy")
        grid_value= Jy(i,j,k) 
	case("Jz")
	grid_value= Jz(i,j,k)
	case("J.E")
	grid_value= Ex(i,j,k)*Jx(i,j,k)+Ey(i,j,k)*Jy(i,j,k)+Ez(i,j,k)*Jz(i,j,k)
	case("E^2")
	grid_value=  (Ex(i,j,k)**2)+ (Ey(i,j,k)**2)+ (Ez(i,j,k)**2)
	case("ave_E^2")
	grid_value=  sqrt((Ex(i,j,k)**2)+ (Ey(i,j,k)**2)+ (Ez(i,j,k)**2))
!	print *,component,grid_value
	case("H^2")
	grid_value=  (Hx(i,j,k)**2)+ (Hy(i,j,k)**2)+ (Hz(i,j,k)**2) 
	case("ave_H^2")
	grid_value=  sqrt((Hx(i,j,k)**2)+ (Hy(i,j,k)**2)+ (Hz(i,j,k)**2))
	case("Sx")
	grid_value= Ey(i,j,k)*Hz(i,j,k)-Ez(i,j,k)*Hy(i,j,k) 
	case("Sy")
	grid_value= Ez(i,j,k)*Hx(i,j,k)-Ex(i,j,k)*Hz(i,j,k) 
	case("Sz")
	grid_value= Ex(i,j,k)*Hy(i,j,k)-Ey(i,j,k)*Hx(i,j,k) 
	case("LogE^2")
	grid_value= log10( (Ex(i,j,k)**2)+ (Ey(i,j,k)**2)+ (Ez(i,j,k)**2))
	case("LogH^2")
	grid_value= log10( (Hx(i,j,k)**2)+ (Hy(i,j,k)**2)+ (Hz(i,j,k)**2)) 
        case("E_Energy")
	if( meps(i,j,k)==0.0) then
	    grid_value= 0.5*eo*eps(i,j,k)*( (Ex(i,j,k)**2)+ (Ey(i,j,k)**2)+ (Ez(i,j,k)**2)) 
	else
	    grid_value= 0.5*eo*eps_m(i,j,k)*( (Ex(i,j,k)**2)+ (Ey(i,j,k)**2)+ (Ez(i,j,k)**2)) 
	endif
	case("K_Energy")
	if( meps(i,j,k)==0.0) then 
	    grid_value= 0.0 
	else
	    grid_value= 0.5*eo*eps_mK(i,j,k)*( (Ex(i,j,k)**2)+ (Ey(i,j,k)**2)+ (Ez(i,j,k)**2)) 
	endif
        case("M_Energy")	
	if(meps(i,j,k)==0.0) then 
	    grid_value= 0.5*eo*eps(i,j,k)*( (Ex(i,j,k)**2)+ (Ey(i,j,k)**2)+ (Ez(i,j,k)**2)) 
	else
	    grid_value= 0.5*eo*eps_mM(i,j,k)*( (Ex(i,j,k)**2)+ (Ey(i,j,k)**2)+ (Ez(i,j,k)**2)) 
	endif
        case("M_Energy2")
	    grid_value= 0.5*uo*( (Hx(i,j,k)**2)+ (Hy(i,j,k)**2)+ (Hz(i,j,k)**2)) 
	case("Ex2Ey2")
	    grid_value=  (Ex(i,j,k)**2)+ (Ey(i,j,k)**2) 
	case("Ey2Ez2")
	    grid_value=  (Ey(i,j,k)**2)+ (Ez(i,j,k)**2) 
	case("Ez2Ex2")
	    grid_value=  (Ez(i,j,k)**2)+ (Ex(i,j,k)**2) 
	case("LogEx2Ey2")
	    grid_value= log10( (Ex(i,j,k)**2)+ (Ey(i,j,k)**2)) 
        case("LogEy2Ez2")
	    grid_value= log10( (Ey(i,j,k)**2)+ (Ez(i,j,k)**2)) 
	case("LogEz2Ex2")
	    grid_value= log10( (Ez(i,j,k)**2)+ (Ex(i,j,k)**2)) 
	case("Hx2Hy2")
	    grid_value=  (Hx(i,j,k)**2)+ (Hy(i,j,k)**2) 
	case("Hy2Hz2")
	    grid_value=  (Hy(i,j,k)**2)+ (Hz(i,j,k)**2) 
	case("Hz2Hx2")
	    grid_value=  (Hz(i,j,k)**2)+ (Hx(i,j,k)**2) 
	case("LogHx2Hy2")
	    grid_value= log10( (Hx(i,j,k)**2)+ (Hy(i,j,k)**2)) 
	case("LogHy2Hz2")
	    grid_value= log10( (Hy(i,j,k)**2)+ (Hz(i,j,k)**2)) 
	case("LogHz2Hx2")
	    grid_value= log10( (Hz(i,j,k)**2)+ (Hx(i,j,k)**2)) 
	end select
end function grid_value


function grid_value_Rw(component,m,i,j,k)
character(len=*)::component
integer::i,j,k,m
real*8::grid_value_Rw

        select case(component)
	case("Ex") 
        grid_value_Rw= Ex_Rw(m,i,j,k) 
	case("Ey")
        grid_value_Rw= Ey_Rw(m,i,j,k) 
	case("Ez")
        grid_value_Rw= Ez_Rw(m,i,j,k) 
        case("Hx")
        grid_value_Rw= Hx_Rw(m,i,j,k) 
	case("Hy")
        grid_value_Rw= Hy_Rw(m,i,j,k) 
	case("Hz")
        grid_value_Rw= Hz_Rw(m,i,j,k) 
	case("Sx")
	grid_value_Rw= Ey_Rw(m,i,j,k)*Hz_Rw(m,i,j,k)-Ez_Rw(m,i,j,k)*Hy_Rw(m,i,j,k)
	case("Sy")
	grid_value_Rw= Ez_Rw(m,i,j,k)*Hx_Rw(m,i,j,k)-Ex_Rw(m,i,j,k)*Hz_Rw(m,i,j,k)
	case("Sz")
	grid_value_Rw= Ex_Rw(m,i,j,k)*Hy_Rw(m,i,j,k)-Ey_Rw(m,i,j,k)*Hx_Rw(m,i,j,k)
	end select
end function grid_value_Rw	


function eps(  i,  j,  k)
integer::i,j,k
real*8::eps
!?
	eps=(epsilonx(i,j,k)+epsilonx(i-1,j,k)+epsilony(i,j,k)+epsilony(i,j-1,k)+epsilonz(i,j,k)+epsilonz(i,j,k-1))/6 
!	eps=epsilonx(i,j,k) 
end function eps 


function meps(  i,  j,  k)
integer::i,j,k
real*8::meps 
	if(metal_type == 0) then
		meps=0.0
	elseif(metal_type == 1  ) then
		meps=Lsigma(i,j,k)
		
	elseif(metal_type == 2  ) then
		meps=Lepsilon(i,j,k)
	endif
end function meps

function eps_m(  i,  j,  k)
integer::i,j,k
real*8::eps_m 
real*8:: wp, w, go 


	wp = Lomega(i,j,k) 
	go = Lgamma(i,j,k) 
	w = global_W*2*pi*light_speed/(ds_x*lattice_x)  
!// w((real)
	eps_m=(Lepsilon(i,j,k) + wp*wp*(w*w-go*go)/((w*w+go*go)*(w*w+go*go)))    
end function eps_m

function eps_mK(  i,  j,  k) !//for calculting kinetic energy of electrons )
integer::i,j,k
real*8::eps_mK
real*8:: wp, w, go 


	wp = Lomega(i,j,k) 
	go = Lgamma(i,j,k) 
	w = global_W*2*3.1415926*light_speed/(ds_x*lattice_x)  

	eps_mK=(wp*wp/(w*w+go*go))    !// changed (ver. X2. 510)
end function eps_mK

function eps_mM(  i,  j,  k) !//for calculting magnetic energy 
integer::i,j,k
real*8::eps_mM 
real*8::wp, w, go 


	wp = Lomega(i,j,k) 
	go = Lgamma(i,j,k) 
	w = global_W*2*3.1415926*light_speed/(ds_x*lattice_x)  

	eps_mM=(Lepsilon(i,j,k) - wp*wp/(w*w+go*go)) 
end function eps_mM

subroutine Poynting_block_in(cname,centerx,centery,centerz,size1,size2,size3,ti,tf)
character(*)::cname	
real*8::centerx,centery,centerz,size1,size2,size3,scale_factor
integer:: i,j,k,ti,tf,m
real*8::sumxp,sumxn,sumyp,sumyn,sumzp,sumzn,sumx,sumy,sumz,lambda
integer:: x_n,x_p, y_n,y_p, z_n,z_p  !  //n:- p:+, position of the sides of the block

	x_n=(centerx + xcenter - 0.5*size1)*lattice_x
	x_p=(centerx + xcenter + 0.5*size1)*lattice_x
	y_n=(centery + ycenter - 0.5*size2)*lattice_y
	y_p=(centery + ycenter + 0.5*size2)*lattice_y
	z_n=(centerz + zcenter - 0.5*size3)*lattice_z
	z_p=(centerz + zcenter + 0.5*size3)*lattice_z
	
	scale_factor=l_n*1.d6/lattice_x

open(28,file='absorb.dat')	
do m=1,Nf
	lambda=l_n*1.d9/(fmin+df*(m-1))

	sumx=0.d0;sumy=0.d0;sumz=0.d0
	sumxp=0.d0;sumxn=0.d0;sumyp=0.d0;sumyp=0.d0;sumzp=0.d0;sumzn=0.d0
	do k=z_n,z_p,1
		do j=y_n,y_p,1 
			sumxp=sumxp+(ds_y*ds_z)*grid_value_Rw("Sx",m,x_p,j,k)
			sumxn=sumxn+(ds_y*ds_z)*grid_value_Rw("Sx",m,x_n,j,k)
		enddo
	enddo	

	do i=x_n,x_p,1
		do k=z_n,z_p,1
			sumyp=sumyp+(ds_z*ds_x)*grid_value_Rw("Sy",m,i,y_p,k)
			sumyn=sumyn+(ds_z*ds_x)*grid_value_Rw("Sy",m,i,y_n,k)
		enddo
	enddo
	do i=x_n,x_p,1
		do j=y_n,y_p,1 
			sumzp=sumzp+(ds_x*ds_y)*grid_value_Rw("Sz",m,i,j,z_p)
			sumzn=sumzn+(ds_x*ds_y)*grid_value_Rw("Sz",m,i,j,z_n)
		enddo
	enddo

	sumxp=sumxp*scale_factor*scale_factor
	sumxn=sumxn*scale_factor*scale_factor
	sumyp=sumyp*scale_factor*scale_factor
	sumyn=sumyn*scale_factor*scale_factor
	sumzp=sumzp*scale_factor*scale_factor
	sumzn=sumzn*scale_factor*scale_factor
	Sumx=(sumxp-sumxn)	
	Sumy=(sumyp-sumyn)
	Sumz=(sumzp-sumzn)

write(28,'(7(E12.5,3x))')lambda,sumxp,sumxn,sumyp,sumyn,sumzp,sumzn   

enddo
close(28)
!		if(t==ti) then
!                   open(28,file=cname,status='REPLACE')
!	           close(28)
!		endif
!		open(unit=28,file=cname,position='APPEND') 
!		write(28,'(7(E12.5,3x))') dble(t)*timestep, sumxp,sumxn,sumyp,sumyn,sumzp,sumzn   
!		write(28,'(5(E12.5,3x))') dble(t)*timestep, Sumx,Sumy,Sumz,(Sumx+Sumy+Sumz)   
!		write(28,'(2(E12.5,3x))') dble(t)*timestep,(Sumx+Sumy+Sumz)   
!		close(28) 
end subroutine Poynting_block_in

subroutine Poynting_block_out(cname,centerx,centery,centerz,size1,size2,size3,ti,tf)
character(*)::cname	
real*8::centerx,centery,centerz,size1,size2,size3
integer:: i,j,k,ti,tf
real*8::sumxp,sumxn,sumyp,sumyn,sumzp,sumzn,sumx,sumy,sumz
integer:: x_n,x_p, y_n,y_p, z_n,z_p  !  //n:- p:+, position of the sides of the block

	x_n=(centerx + xcenter - 0.5*size1)*lattice_x
	x_p=(centerx + xcenter + 0.5*size1)*lattice_x
	y_n=(centery + ycenter - 0.5*size2)*lattice_y
	y_p=(centery + ycenter + 0.5*size2)*lattice_y
	z_n=(centerz + zcenter - 0.5*size3)*lattice_z
	z_p=(centerz + zcenter + 0.5*size3)*lattice_z

	sumx=0.d0;sumy=0.d0;sumz=0.d0
	sumxp=0.d0;sumxn=0.d0;sumyp=0.d0;sumyp=0.d0;sumzp=0.d0;sumzn=0.d0


	do k=z_n,z_p,1
		do j=y_n,y_p,1 
			sumxp=sumxp+(ds_y*ds_z)*grid_value("Sx",x_p,j,k)
			sumxn=sumxn+(ds_y*ds_z)*grid_value("Sx",x_n,j,k)
		enddo
	enddo		
	do i=x_n,x_p,1
		do k=z_n,z_p,1
			sumyp=sumyp+(ds_z*ds_x)*grid_value("Sy",i,y_p,k)
			sumyn=sumyn+(ds_z*ds_x)*grid_value("Sy",i,y_n,k)
		enddo
	enddo
	do i=x_n,x_p,1
		do j=y_n,y_p,1 
			sumzp=sumzp+(ds_x*ds_y)*grid_value("Sz",i,j,z_p)
			sumzn=sumzn+(ds_x*ds_y)*grid_value("Sz",i,j,z_n)
		enddo
	enddo

	Sumx=(sumxp-sumxn)	
	Sumy=(sumyp-sumyn)
	Sumz=(sumzp-sumzn)


		if(t==ti) then
                   open(28,file=cname,status='REPLACE')
	           close(28)
		endif
		open(unit=28,file=cname,position='APPEND') 
		write(28,'(5(E12.5,3x))') dble(t)*timestep, Sumx,Sumy,Sumz,(Sumx+Sumy+Sumz)   
!		write(28,'(2(E12.5,3x))') dble(t)*timestep,(Sumx+Sumy+Sumz)   
		close(28) 
end subroutine Poynting_block_out
subroutine incident_intensity(ks,ti,tf)
real*8::sx,sy,sz
integer::ks,ti,tf

!average S_inc=|E_0*H_0|/2=sqrt(epsr)*eps0*c*E_inc^2/2=sqrt(epsr)/mu0/c*E_inc^2/2

!sx= Ey(i,j,k)*Hz(i,j,k)-Ez(i,j,k)*Hy(i,j,k)
sx=0.d0	
!sy= Ez(i,j,k)*Hx(i,j,k)-Ex(i,j,k)*Hz(i,j,k) 
sy=0.d0
!sz= Ex(i,j,k)*Hy(i,j,k)-Ey(i,j,k)*Hx(i,j,k) 
sz=Ex_inc(ks)*Hy_inc(ks)	

if(t==ti) then
        open(30,file='inc_int.dat',status='REPLACE')
	close(30)
endif
	open(unit=30,file='inc_int.dat',position='APPEND') 
	write(30,'(2(E12.5,3x))') dble(t)*timestep,(sx+sy+sz)   
	close(30) 
end subroutine incident_intensity


subroutine total_E_energy_block(centerx,centery,centerz,size1,size2,size3,ti,tf)
real*8::centerx,centery,centerz,size1,size2,size3,Eb
integer:: i,j,k,ti,tf
integer:: x_n,x_p, y_n,y_p, z_n,z_p  ! //n:- p:+, position of the sides of the block

	x_n=(centerx + xcenter - 0.5*size1)*lattice_x
	x_p=(centerx + xcenter + 0.5*size1)*lattice_x
	y_n=(centery + ycenter - 0.5*size2)*lattice_y
	y_p=(centery + ycenter + 0.5*size2)*lattice_y
	z_n=(centerz + zcenter - 0.5*size3)*lattice_z
	z_p=(centerz + zcenter + 0.5*size3)*lattice_z
	
	do i=x_n,x_p,1
		do k=z_p,z_n,-1
			do j=y_n,y_p,1 
				Eb = Eb + ds_x*ds_y*ds_z*grid_value("E_Energy",i,j,k) 
			enddo
		enddo
	enddo
		if(t==ti) then
                   open(29,file='total_eEE_b.en',status='REPLACE')
	           close(29)
		endif

	open(29,file='total_eEE_b.en',position='APPEND')
	write(29,'(2(E12.5,3x))') dble(t)*timestep,Eb   
	close(29) 
end subroutine total_E_energy_block
end module output
