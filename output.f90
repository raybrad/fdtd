module output
use parameters

!public::real_space_param,out_epsilon,out_plane,out_plane_cut,out_point,out_line_z,out_line_x,out_line_y,lodestar_qmcal,Poynting_block_in,Poynting_block_out, &
!        incident_intensity,out_point_w,total_E_energy_block

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
	!timestep=(a_m)/(light_speed*S_factor*ds_x*lattice_x)	!(s)
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

function  get_period_in_update(w_n)
integer::get_period_in_update
real*8::w_n
 
	get_period_in_update=int((S_factor*ds_x*lattice_x)/w_n)
	
end function get_period_in_update

subroutine out_epsilon(plane,value,cname)
character::plane
character*10::cname
real*8::value
integer:: i,j,k,iu 
integer:: i_range, j_range, k_range  
integer:: mz  !// index for non_uniform_grid multiplication 

	!// for Hz_parity, normal cases
	i_range = isize-1 
	j_range = jsize-1 
	k_range = ksize-1  
	
	iu=17

	open(unit=iu,file=cname,status='replace',form='FORMATTED') 

	if(plane=="x") then
		i=floor(0.5+((value+xcenter)*lattice_x)) 
		do k=0,k_range
				do j=0,j_range
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
	 
		j=floor(0.5+((value+ycenter)*lattice_y)) 
		do k=0,k_range
				do i=0,i_range
				 
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
	 
		k=floor(0.5+((value+zcenter)*lattice_z)) 
		do j=0,j_range
			do i=0,i_range
			 
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

	i_range = isize-1 
	j_range = jsize-1 
	k_range = ksize-1  

	iu=18
	
	write(cname,'(I7.7)') t
         cname=trim(cname)//trim(lastname) 
	open(unit=iu,file=cname,status='replace',form='formatted') 

	if(plane=="x") then
		i=floor(0.5+((value+xcenter)*lattice_x))
		do k=0,k_range
				do j=0,j_range	
                                write(iu,'(E20.12,3x)',advance='no') grid_value_bohr(component,i,j,k) 
	                        enddo
				write(iu,*)  
		enddo
	endif

	if(plane=="y") then
	 
		j=floor(0.5+((value+ycenter)*lattice_y)) 
		do k=0,k_range
				do i=0,i_range
                                write(iu,'(E20.12,3x)',advance='no') grid_value_bohr(component,i,j,k)  
	                        !print *,'grid_value_bohr',grid_value_bohr(component,i,j,k)
	                        enddo
				write(iu,*) 
		enddo
	endif

	if(plane=="z") then
	 
		k=floor(0.5+((value+zcenter)*lattice_z)) 
		do j=0,j_range
			do i=0,i_range
                        write(iu,'(E20.12,3x)',advance='no') grid_value_bohr(component,i,j,k)  	
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

	i_range = isize-1-pmlir 
	j_range = jsize-1-pmljr 
	k_range = ksize-1-pmlkr  

	iu=18
	
	write(cname,'(I7.7)') t
         cname=trim(cname)//trim(lastname)

	write(6,*) 'i j k range',i_range,j_range,k_range
	open(unit=iu,file=cname,status='replace',form='formatted') 

	if(plane=="x") then
		i=floor(0.5+((value+xcenter)*lattice_x))
		do k=1+pmlkl,k_range
				do j=1+pmljl,j_range	
                                write(iu,'(E20.12,3x)',advance='no') grid_value_bohr(component,i,j,k) 
	                        enddo
				write(iu,*)  
		enddo
	endif

	if(plane=="y") then
	 
		j=floor(0.5+((value+ycenter)*lattice_y)) 
		do k=1+pmlkl,k_range
				do i=1+pmlil,i_range
                                write(iu,'(E20.12,3x)',advance='no') grid_value_bohr(component,i,j,k)  
	                        !print *,'grid_value_bohr',grid_value_bohr(component,i,j,k)
	                        enddo
				write(iu,*) 
		enddo
	endif

	if(plane=="z") then
	 
		k=floor(0.5+((value+zcenter)*lattice_z)) 
		do j=1+pmljl,j_range
			do i=1+pmlil,i_range
                        write(iu,'(E20.12,3x)',advance='no') grid_value_bohr(component,i,j,k) 	
	                enddo
			write(iu,*)  
		enddo
	endif
	close(iu) 
	print *,"out",' ',component," ...ok " 
end subroutine out_plane_cut

subroutine out_plane_cut_w(component,outLambda,plane,value,lastname)
 
    character(len=*)::plane
	real*8::value
	character(len=*)::component,lastname
	character*20::cname
    real*8::outLambda
	integer:: i,j,k,iu,m
	integer:: i_range, j_range, k_range  
	integer:: mz  !//for non_uniform_grid multiplication

	i_range = isize-1-pmlir 
	j_range = jsize-1-pmljr 
	k_range = ksize-1-pmlkr  

	m=int((l_n/outLambda-fmin)/df)+1
	iu=18
	
	write(cname,'(I7.7)') t
    cname=trim(cname)//trim(lastname)

	write(6,*) 'i j k range',i_range,j_range,k_range
	open(unit=iu,file=cname,status='replace',form='formatted') 

	if(plane=="x") then
		i=floor(0.5+((value+xcenter)*lattice_x))
		do k=1+pmlkl,k_range
				do j=1+pmljl,j_range	
                                write(iu,'(E20.12,3x)',advance='no') grid_value_bohr_fft(component,m,i,j,k) 
	            enddo
				write(iu,*)  
		enddo
	endif

	if(plane=="y") then
	 
		j=floor(0.5+((value+ycenter)*lattice_y)) 
		do k=1+pmlkl,k_range
				do i=1+pmlil,i_range
                                write(iu,'(E20.12,3x)',advance='no') grid_value_bohr_fft(component,m,i,j,k)  
	                        !print *,'grid_value_bohr_fft',grid_value_bohr_fft(component,i,j,k)
	                        enddo
				write(iu,*) 
		enddo
	endif

	if(plane=="z") then
	 
		k=floor(0.5+((value+zcenter)*lattice_z)) 
		do j=1+pmljl,j_range
			do i=1+pmlil,i_range
                        write(iu,'(E20.12,3x)',advance='no') grid_value_bohr_fft(component,m,i,j,k) 	
	                enddo
			write(iu,*)  
		enddo
	endif
	close(iu) 
	print *,"out",' ',component," ...ok " 
end subroutine out_plane_cut_w



subroutine out_point(component,x,y,z,ti,tf,cname)
 
        character(len=*)::component
	character(len=*)::cname
        real*8::x,y,z
	integer::ti,tf
	integer:: i,j,k 

	if(ti<=t .and. t<=tf) then
       !  print *,'t ',t	
	      
		i=floor(0.5+((x+xcenter)*lattice_x)) 
		j=floor(0.5+((y+ycenter)*lattice_y)) 
		k=floor(0.5+((z+zcenter)*lattice_z)) 
		if(t==ti) then
                   open(18,file=cname,status='REPLACE')
	           close(18)
		endif
		open(unit=18,file=cname,position='APPEND') 
		write(18,'(2(E20.12,3x))') dble(t)*timestep, grid_value_bohr(component,i,j,k)    !if it's E^2 should be grid_value_bohr ^2
		close(18) 
     endif
end subroutine out_point

subroutine out_line_z(component,x,y,z1,z2,ti,tf,cname)
        character(len=*)::component
	character(len=*)::cname
        real*8::x,y,z1,z2
	integer::ti,tf
	integer:: i,j,k,k1,k2 

	if(ti<=t .and. t<=tf) then
       !  print *,'t ',t	
	      
		i=floor(0.5+((x+xcenter)*lattice_x)) 
		j=floor(0.5+((y+ycenter)*lattice_y)) 
		k1=floor(0.5+((z1+zcenter)*lattice_z)) 
		k2=floor(0.5+((z2+zcenter)*lattice_z)) 
		if(t==ti) then
                   open(18,file=cname,status='REPLACE')
	           close(18)
		endif
		open(unit=18,file=cname,position='APPEND')
		write(18,'(E20.12,3X)',advance='no') dble(t)*timestep
		do k=k1,k2
		write(18,'(E20.12,3x)',advance='no') grid_value_bohr(component,i,j,k)    !if it's E^2 should be grid_value_bohr ^2
		enddo
		close(18) 
     endif
end subroutine out_line_z

subroutine out_line_x(component,x1,x2,y,z,ti,tf,cname)
        character(len=*)::component
	character(len=*)::cname
        real*8::x1,x2,y,z
	integer::ti,tf
	integer::i, i1,i2,j,k 

	if(ti<=t .and. t<=tf) then
       !  print *,'t ',t	
	      
		i1=floor(0.5+((x1+xcenter)*lattice_x)) 
		i2=floor(0.5+((x2+xcenter)*lattice_x)) 
		j=floor(0.5+((y+ycenter)*lattice_y)) 
		k=floor(0.5+((z+zcenter)*lattice_z)) 
		if(t==ti) then
                   open(18,file=cname,status='REPLACE')
	           close(18)
		endif
		open(unit=18,file=cname,position='APPEND')
		write(18,'(E20.12,3X)',advance='no') dble(t)*timestep
		do i=i1,i2
		write(18,'(E20.12,3x)',advance='no') grid_value_bohr(component,i,j,k)    !if it's E^2 should be grid_value_bohr ^2
		enddo
		close(18) 
     endif
end subroutine out_line_x

subroutine out_line_y(component,x,y1,y2,z,ti,tf,cname)
        character(len=*)::component
	character(len=*)::cname
        real*8::y1,y2,z,x
	integer::ti,tf
	integer:: i,j,j1,j2,k 

	if(ti<=t .and. t<=tf) then
       !  print *,'t ',t	
	      
		i=floor(0.5+((x+xcenter)*lattice_x)) 
		j1=floor(0.5+((y1+ycenter)*lattice_y)) 
		j2=floor(0.5+((y2+ycenter)*lattice_y)) 
		k=floor(0.5+((z+zcenter)*lattice_z)) 
		if(t==ti) then
                   open(18,file=cname,status='REPLACE')
	           close(18)
		endif
		open(unit=18,file=cname,position='APPEND')
		write(18,'(E20.12,3X)',advance='no') dble(t)*timestep
		do j=j1,j2
		write(18,'(E20.12,3x)',advance='no') grid_value_bohr(component,i,j,k)    !if it's E^2 should be grid_value_bohr ^2
		enddo
		close(18) 
     endif
end subroutine out_line_y
subroutine lodestar_qmcal(x,y,z)
! input Efield at molecular pos in  lodestar to get current at interface
real*8::x,y,z
integer::i,j,k
INTEGER*4 getcwd
INTEGER*4 istat!,system
logical*4::result
CHARACTER(len=1000) :: path
CHARACTER(len=1000) :: command
i=floor(0.5+((x+xcenter)*lattice_x))
j=floor(0.5+((y+ycenter)*lattice_y))
k=floor(0.5+((z+zcenter)*lattice_z)) 
open(159,file='fdembound') !v/m *0.529*1d-10 =>eV/bohr
write(159,'(4(E20.5,X))')dble(t)*timestep,grid_value_bohr('Ex',i,j,k),grid_value_bohr('Ey',i,j,k),grid_value_bohr('Ez',i,j,k)
close(159)

open(160,file='Embound.dat',position='append',form='formatted') !v/m *0.529*1d-10 =>eV/bohr
write(160,'(4(E20.5,X))')dble(t)*timestep,grid_value_bohr('Ex',i,j,k),grid_value_bohr('Ey',i,j,k),grid_value_bohr('Ez',i,j,k)
close(160)
istat=getcwd(path)
!command="/home/hy/lodestar/bin/lodestartd <"//trim(path)//"/in.td>"//trim(path)//"/out.td" 
command='touch file.dat'
write(6,*) trim(command)
call system(trim(command),istat)
!result=systemqq(trim(command))
!write(6,*) 'result',result
!istat = GETLASTERRORQQ( )
!write(6,*) 'istat',istat
!call system(trim(command),istat) 
if (istat==0) then
!if(result) then
 print *,'success in qm calculation'
 continue
else
 print *,'  error in qm calculation',istat
 stop
endif

end subroutine lodestar_qmcal

function  grid_value_bohr(component,i,j,k)
character(len=*)::component
integer::i,j,k
real*8::grid_value_bohr
        select case(component)
	case("Ex") 
        grid_value_bohr= Ex(i,j,k)*bohr 
    case("Ex_inc")
    ks=light_pos 
    grid_value_bohr=  Ex_inc(ks)*bohr
	case("Ey")
	grid_value_bohr= Ey(i,j,k)*bohr
	case("Ez")
	grid_value_bohr= Ez(i,j,k)*bohr
        case("Hx")
	grid_value_bohr= Hx(i,j,k)*bohr
	case("Hy")
	grid_value_bohr= Hy(i,j,k)*bohr
	case("Hz")
	grid_value_bohr= Hz(i,j,k)*bohr
	case("Jx")
	grid_value_bohr= Jx(i,j,k)*bohr*bohr
	case("Jy")
        grid_value_bohr= Jy(i,j,k)*bohr*bohr 
	case("Jz")
	grid_value_bohr= Jz(i,j,k)*bohr*bohr
    case("Ex_inc^2")
    ks=light_pos 
    grid_value_bohr=  (Ex_inc(ks)*Ex_inc(ks))*bohr*bohr
	case("E^2")
	grid_value_bohr=  ((Ex(i,j,k)**2)+ (Ey(i,j,k)**2)+ (Ez(i,j,k)**2))*bohr*bohr
	case("ave_E^2")
	grid_value_bohr=  sqrt((Ex(i,j,k)**2)+ (Ey(i,j,k)**2)+ (Ez(i,j,k)**2))*bohr
	case("H^2")
	grid_value_bohr=  (Hx(i,j,k)**2)+ (Hy(i,j,k)**2)+ (Hz(i,j,k)**2)*bohr*bohr 
	case("ave_H^2")
	grid_value_bohr=  sqrt((Hx(i,j,k)**2)+ (Hy(i,j,k)**2)+ (Hz(i,j,k)**2))*bohr
	case("Sx")
	grid_value_bohr= Ey(i,j,k)*Hz(i,j,k)-Ez(i,j,k)*Hy(i,j,k)*bohr*bohr 
	case("Sy")
	grid_value_bohr= Ez(i,j,k)*Hx(i,j,k)-Ex(i,j,k)*Hz(i,j,k)*bohr*bohr 
	case("Sz")
	grid_value_bohr= Ex(i,j,k)*Hy(i,j,k)-Ey(i,j,k)*Hx(i,j,k)*bohr*bohr 
	case("LogE^2")
	grid_value_bohr= log10( (Ex(i,j,k)**2)+ (Ey(i,j,k)**2)+(Ez(i,j,k)**2))*bohr*bohr
	case("LogH^2")
	grid_value_bohr= log10( (Hx(i,j,k)**2)+ (Hy(i,j,k)**2)+(Hz(i,j,k)**2))*bohr*bohr 
        case("E_Energy")
	if( meps(i,j,k)==0.0) then
	    grid_value_bohr= 0.5*eo*eps(i,j,k)*( (Ex(i,j,k)**2)+ (Ey(i,j,k)**2)+ (Ez(i,j,k)**2))*bohr*bohr 
	else
	    grid_value_bohr= 0.5*eo*eps_m(i,j,k)*( (Ex(i,j,k)**2)+ (Ey(i,j,k)**2)+ (Ez(i,j,k)**2))*bohr*bohr 
	endif
	case("K_Energy")
	if( meps(i,j,k)==0.0) then 
	    grid_value_bohr= 0.0 
	else
	    grid_value_bohr= 0.5*eo*eps_mK(i,j,k)*( (Ex(i,j,k)**2)+ (Ey(i,j,k)**2)+ (Ez(i,j,k)**2))*bohr*bohr 
	endif
        case("M_Energy")	
	if(meps(i,j,k)==0.0) then 
	    grid_value_bohr= 0.5*eo*eps(i,j,k)*( (Ex(i,j,k)**2)+ (Ey(i,j,k)**2)+ (Ez(i,j,k)**2))*bohr*bohr 
	else
	    grid_value_bohr= 0.5*eo*eps_mM(i,j,k)*( (Ex(i,j,k)**2)+ (Ey(i,j,k)**2)+ (Ez(i,j,k)**2))*bohr*bohr 
	endif
        case("M_Energy2")
	    grid_value_bohr= 0.5*uo*( (Hx(i,j,k)**2)+ (Hy(i,j,k)**2)+ (Hz(i,j,k)**2))*bohr*bohr 
	case("Ex2Ey2")
	    grid_value_bohr=  (Ex(i,j,k)**2)+ (Ey(i,j,k)**2)*bohr*bohr 
	case("Ey2Ez2")
	    grid_value_bohr=  (Ey(i,j,k)**2)+ (Ez(i,j,k)**2) *bohr*bohr 
	case("Ez2Ex2")
	    grid_value_bohr=  (Ez(i,j,k)**2)+ (Ex(i,j,k)**2) *bohr*bohr 
	case("LogEx2Ey2")
	    grid_value_bohr= log10( (Ex(i,j,k)**2)+ (Ey(i,j,k)**2)) *bohr*bohr 
        case("LogEy2Ez2")
	    grid_value_bohr= log10( (Ey(i,j,k)**2)+ (Ez(i,j,k)**2)) *bohr*bohr 
	case("LogEz2Ex2")
	    grid_value_bohr= log10( (Ez(i,j,k)**2)+ (Ex(i,j,k)**2)) *bohr*bohr 
	case("Hx2Hy2")
	    grid_value_bohr=  (Hx(i,j,k)**2)+ (Hy(i,j,k)**2) *bohr*bohr 
	case("Hy2Hz2")
	    grid_value_bohr=  (Hy(i,j,k)**2)+ (Hz(i,j,k)**2) *bohr*bohr 
	case("Hz2Hx2")
	    grid_value_bohr=  (Hz(i,j,k)**2)+ (Hx(i,j,k)**2) *bohr*bohr 
	case("LogHx2Hy2")
	    grid_value_bohr= log10( (Hx(i,j,k)**2)+ (Hy(i,j,k)**2)) *bohr*bohr 
	case("LogHy2Hz2")
	    grid_value_bohr= log10( (Hy(i,j,k)**2)+ (Hz(i,j,k)**2)) *bohr*bohr 
	case("LogHz2Hx2")
	    grid_value_bohr= log10( (Hz(i,j,k)**2)+ (Hx(i,j,k)**2)) *bohr*bohr 
	end select
end function grid_value_bohr

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


function grid_value_bohr_fft(component,m,i,j,k)
character(len=*)::component
integer::i,j,k,m
real*8::grid_value_bohr_fft
real*8::tmp1,tmp2

select case(component)
case("Ex") 
grid_value_bohr_fft= sqrt(Ex_Re(m,i,j,k)*Ex_Re(m,i,j,k)+Ex_Im(m,i,j,k)*Ex_Im(m,i,j,k))*bohr
case("Ey")
grid_value_bohr_fft= sqrt(Ey_Re(m,i,j,k)*Ey_Re(m,i,j,k)+Ey_Im(m,i,j,k)*Ey_Im(m,i,j,k))*bohr 
case("Ez")
grid_value_bohr_fft= sqrt(Ez_Re(m,i,j,k)*Ez_Re(m,i,j,k)+Ez_Im(m,i,j,k)*Ez_Im(m,i,j,k))*bohr 
case("E^2")
grid_value_bohr_fft= ((Ex_Re(m,i,j,k)*Ex_Re(m,i,j,k)+Ex_Im(m,i,j,k)*Ex_Im(m,i,j,k))+(Ey_Re(m,i,j,k)*Ey_Re(m,i,j,k)+Ey_Im(m,i,j,k)*Ey_Im(m,i,j,k))+(Ez_Re(m,i,j,k)*Ez_Re(m,i,j,k)+Ez_Im(m,i,j,k)*Ez_Im(m,i,j,k)))*bohr*bohr
case("Ex_inc^2")
grid_value_bohr_fft=(Ex_inc_Re(m)*Ex_inc_Re(m)+Ex_inc_Im(m)*Ex_inc_Im(m))*bohr*bohr
case("GammaE2")
tmp1= ((Ex_Re(m,i,j,k)*Ex_Re(m,i,j,k)+Ex_Im(m,i,j,k)*Ex_Im(m,i,j,k))+(Ey_Re(m,i,j,k)*Ey_Re(m,i,j,k)+Ey_Im(m,i,j,k)*Ey_Im(m,i,j,k))+(Ez_Re(m,i,j,k)*Ez_Re(m,i,j,k)+Ez_Im(m,i,j,k)*Ez_Im(m,i,j,k)))*bohr*bohr
tmp2= (Ex_inc_Re(m)*Ex_inc_Re(m)+Ex_inc_Im(m)*Ex_inc_Im(m))*bohr*bohr
grid_value_bohr_fft=tmp1/tmp2
case("GammaE")
tmp1= sqrt((Ex_Re(m,i,j,k)*Ex_Re(m,i,j,k)+Ex_Im(m,i,j,k)*Ex_Im(m,i,j,k))+(Ey_Re(m,i,j,k)*Ey_Re(m,i,j,k)+Ey_Im(m,i,j,k)*Ey_Im(m,i,j,k))+(Ez_Re(m,i,j,k)*Ez_Re(m,i,j,k)+Ez_Im(m,i,j,k)*Ez_Im(m,i,j,k)))*bohr
tmp2= sqrt(Ex_inc_Re(m)*Ex_inc_Re(m)+Ex_inc_Im(m)*Ex_inc_Im(m))*bohr
grid_value_bohr_fft=tmp1/tmp2
case("Hx")
grid_value_bohr_fft= sqrt(Hx_Re(m,i,j,k)*Hx_Re(m,i,j,k)+Hx_Im(m,i,j,k)*Hx_Im(m,i,j,k))*bohr
case("Hy")
grid_value_bohr_fft= sqrt(Hy_Re(m,i,j,k)*Hy_Re(m,i,j,k)+Hy_Im(m,i,j,k)*Hy_Im(m,i,j,k))*bohr
case("Hz")
grid_value_bohr_fft= sqrt(Hz_Re(m,i,j,k)*Hz_Re(m,i,j,k)+Hz_Im(m,i,j,k)*Hz_Im(m,i,j,k))*bohr
case("ReSx")	!E*conj(H)
grid_value_bohr_fft= (Ey_Re(m,i,j,k)*Hz_Re(m,i,j,k)+Ey_Im(m,i,j,k)*Hz_Im(m,i,j,k)-Ez_Re(m,i,j,k)*Hy_Re(m,i,j,k)-Ez_Im(m,i,j,k)*Hy_Im(m,i,j,k))*bohr*bohr
case("ReSy")
grid_value_bohr_fft= (Ez_Re(m,i,j,k)*Hx_Re(m,i,j,k)+Ez_Im(m,i,j,k)*Hx_Im(m,i,j,k)-Ex_Re(m,i,j,k)*Hz_Re(m,i,j,k)-Ex_Im(m,i,j,k)*Hz_Im(m,i,j,k))*bohr*bohr
case("ReSz")
grid_value_bohr_fft= (Ex_Re(m,i,j,k)*Hy_Re(m,i,j,k)+Ex_Im(m,i,j,k)*Hy_Im(m,i,j,k)-Ey_Re(m,i,j,k)*Hx_Re(m,i,j,k)-Ey_Im(m,i,j,k)*Hx_Im(m,i,j,k))*bohr*bohr
end select
end function grid_value_bohr_fft	


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
	w = global_W*2*3.1415926*light_speed/(ds_x*lattice_x)  
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
character*20::prename
real*8::centerx,centery,centerz,size1,size2,size3
integer:: i,j,k,ti,tf,m
real*8::sumxp,sumxn,sumyp,sumyn,sumzp,sumzn,sumx,sumy,sumz,lambda,scale_factor
integer:: x_n,x_p, y_n,y_p, z_n,z_p  !  //n:- p:+, position of the sides of the block

	x_n=(centerx + xcenter - 0.5*size1)*lattice_x
	x_p=(centerx + xcenter + 0.5*size1)*lattice_x
	y_n=(centery + ycenter - 0.5*size2)*lattice_y
	y_p=(centery + ycenter + 0.5*size2)*lattice_y
	z_n=(centerz + zcenter - 0.5*size3)*lattice_z
	z_p=(centerz + zcenter + 0.5*size3)*lattice_z
	
	scale_factor=l_n*1.d6/lattice_x	!um
	
	write(prename,'(I7.7)') t
	prename=trim(prename)//cname
	open(28,file=trim(prename))	
do m=1,Nf
	lambda=l_n*1.d9/(fmin+df*(m-1))

	sumx=0.d0;sumy=0.d0;sumz=0.d0
	sumxp=0.d0;sumxn=0.d0;sumyp=0.d0;sumyp=0.d0;sumzp=0.d0;sumzn=0.d0
	do k=z_n,z_p,1
		do j=y_n,y_p,1 
			sumxp=sumxp+(ds_y*ds_z)*grid_value_bohr_fft("ReSx",m,x_p,j,k)
			sumxn=sumxn+(ds_y*ds_z)*grid_value_bohr_fft("ReSx",m,x_n,j,k)
		enddo
	enddo	

	do i=x_n,x_p,1
		do k=z_n,z_p,1
			sumyp=sumyp+(ds_z*ds_x)*grid_value_bohr_fft("ReSy",m,i,y_p,k)
			sumyn=sumyn+(ds_z*ds_x)*grid_value_bohr_fft("ReSy",m,i,y_n,k)
		enddo
	enddo
	do i=x_n,x_p,1
		do j=y_n,y_p,1 
			sumzp=sumzp+(ds_x*ds_y)*grid_value_bohr_fft("ReSz",m,i,j,z_p)
			sumzn=sumzn+(ds_x*ds_y)*grid_value_bohr_fft("ReSz",m,i,j,z_n)
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

!write(28,'(7(E20.12,3x))')lambda,sumxp,sumxn,sumyp,sumyn,sumzp,sumzn   
 write(28,'(5(E20.12,3x))') lambda,(Sumx+Sumy+Sumz),Sumx,Sumy,Sumz 
enddo
close(28)
end subroutine Poynting_block_in

subroutine Poynting_block_out(cname,centerx,centery,centerz,size1,size2,size3,ti,tf)
character(*)::cname
character*20::prename
real*8::centerx,centery,centerz,size1,size2,size3
integer:: i,j,k,ti,tf,m
real*8::sumxp,sumxn,sumyp,sumyn,sumzp,sumzn,sumx,sumy,sumz,lambda,scale_factor
integer:: x_n,x_p, y_n,y_p, z_n,z_p  !  //n:- p:+, position of the sides of the block

	x_n=(centerx + xcenter - 0.5*size1)*lattice_x
	x_p=(centerx + xcenter + 0.5*size1)*lattice_x
	y_n=(centery + ycenter - 0.5*size2)*lattice_y
	y_p=(centery + ycenter + 0.5*size2)*lattice_y
	z_n=(centerz + zcenter - 0.5*size3)*lattice_z
	z_p=(centerz + zcenter + 0.5*size3)*lattice_z
	
	scale_factor=l_n*1.d6/lattice_x	!um
	
	write(prename,'(I7.7)') t
	prename=trim(prename)//cname
	open(28,file=trim(prename))	
do m=1,Nf
	lambda=l_n*1.d9/(fmin+df*(m-1))

	sumx=0.d0;sumy=0.d0;sumz=0.d0
	sumxp=0.d0;sumxn=0.d0;sumyp=0.d0;sumyp=0.d0;sumzp=0.d0;sumzn=0.d0
	do k=z_n,z_p,1
		do j=y_n,y_p,1 
			sumxp=sumxp+(ds_y*ds_z)*grid_value_bohr_fft("ReSx",m,x_p,j,k)
			sumxn=sumxn+(ds_y*ds_z)*grid_value_bohr_fft("ReSx",m,x_n,j,k)
		enddo
	enddo	

	do i=x_n,x_p,1
		do k=z_n,z_p,1
			sumyp=sumyp+(ds_z*ds_x)*grid_value_bohr_fft("ReSy",m,i,y_p,k)
			sumyn=sumyn+(ds_z*ds_x)*grid_value_bohr_fft("ReSy",m,i,y_n,k)
		enddo
	enddo
	do i=x_n,x_p,1
		do j=y_n,y_p,1 
			sumzp=sumzp+(ds_x*ds_y)*grid_value_bohr_fft("ReSz",m,i,j,z_p)
			sumzn=sumzn+(ds_x*ds_y)*grid_value_bohr_fft("ReSz",m,i,j,z_n)
		enddo
	enddo
	
	sumxp=sumxp*scale_factor*scale_factor
	sumxn=sumxn*scale_factor*scale_factor
	sumyp=sumyp*scale_factor*scale_factor
	sumyn=sumyn*scale_factor*scale_factor
	sumzp=sumzp*scale_factor*scale_factor
	sumzn=sumzn*scale_factor*scale_factor

	Sumx= -(sumxp-sumxn)	
	Sumy= -(sumyp-sumyn)
	Sumz= -(sumzp-sumzn)

!write(28,'(7(E20.12,3x))')lambda,sumxp,sumxn,sumyp,sumyn,sumzp,sumzn   
 write(28,'(5(E20.12,3x))') lambda,(Sumx+Sumy+Sumz),Sumx,Sumy,Sumz		!um^2 
enddo
close(28)
end subroutine Poynting_block_out


subroutine incident_intensity(posz,centerx,centery,size1,size2)
real*8::sx,sy,sz,lambda,centerx,centery,size1,size2,sumzp,posz
integer::ks,m,i,j,x_n,x_p,y_n,y_p,z_p
character*20::cname

	x_n=(centerx + xcenter - 0.5*size1)*lattice_x
	x_p=(centerx + xcenter + 0.5*size1)*lattice_x
	y_n=(centery + ycenter - 0.5*size2)*lattice_y
	y_p=(centery + ycenter + 0.5*size2)*lattice_y
	z_p=(posz + zcenter)*lattice_z


	write(cname,'(I7.7)') t
	cname=trim(cname)//'inc_int.dat'
	open(unit=30,file=trim(cname)) 
do m=1,Nf
	lambda=l_n*1.d9/(fmin+df*(m-1))
	sumzp=0.d0
	do i=x_n,x_p,1
		do j=y_n,y_p,1 
			sumzp=sumzp+ds_x*ds_y*grid_value_bohr_fft("ReSz",m,i,j,z_p)
		enddo
	enddo
	write(30,'(2(E20.12,3x))') lambda,sumzp/(x_p-x_n)/(y_p-y_n)  !W/m^2 
enddo	
	close(30) 
end subroutine incident_intensity

subroutine out_point_w(component,value,lastname)
integer::m,x_p
real*8::value
character(len=*)::component,lastname
real*8::lambda
character*20::cname

x_p=floor(0.5+(value + xcenter)*lattice_x)
write(cname,'(I7.7)') t
cname=trim(cname)//trim(lastname)
open(40,file=trim(cname))
do m=1,Nf
	lambda=l_n*1.d9/(fmin+df*(m-1))
	write(40,*) lambda,grid_value_bohr_fft(component,m,x_p,isize/2,jsize/2)
enddo
close(40)
end subroutine out_point_w

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
	write(29,'(2(E20.12,3x))') dble(t)*timestep,Eb   
	close(29) 
end subroutine total_E_energy_block
end module output
