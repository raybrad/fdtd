module source
use parameters

contains
!
subroutine Gaussian_planewave(Ecomp,Hcomp,rpos,frequency,t0,tdecay)
!
character*10::Ecomp,Hcomp
real*8:: rpos,frequency
integer::t0,tdecay
integer:: i,j,k 
real*8:: phase 	
real*8::amplitude0

!example
!E~/V/m; H~A/m
!
amplitude0=amplitude/(5.291772083d-11)  !V/bohr -> V/m

!	k=floor(0.5+((rpos0+zcenter)*lattice_z)) 

	do j=0,pjsize
         do i=0,pisize
!	do j=1,pjsize-1
!	        do i=1,pisize-1
			phase = 0.0  
			call Gaussian_dipole_source_beam(Ecomp, i, j, rpos, frequency, amplitude0, phase, t0, tdecay)
			call Gaussian_dipole_source_beam(Hcomp, i, j, rpos, frequency, amplitude0, phase, t0, tdecay) 
		enddo
	enddo
end subroutine Gaussian_planewave

subroutine Gaussian_point_source(Ecomp,Hcomp,posx,posy,posz,frequency,t0,tdecay)
!
character*10::Ecomp,Hcomp
real*8:: posx,posy,posz,frequency
integer::t0,tdecay
integer:: i,j,k 
real*8:: phase 	
real*8::amplitude0

!example
!E~/V/m; H~A/m
!
amplitude0=amplitude/(5.291772083d-11)  !V/bohr -> V/m
phase = 0.0  

i=floor(0.5+((posx+xcenter)*lattice_x)) 
j=floor(0.5+((posy+ycenter)*lattice_y))

call Gaussian_dipole_source(Ecomp,i,j,posz,frequency,amplitude0,phase,t0,tdecay)
call Gaussian_dipole_source(Hcomp,i,j,posz,frequency,amplitude0,phase,t0,tdecay)

end subroutine Gaussian_point_source

subroutine Gaussian_dipole_source(component,i,j,z0,frequency,amp,phase,t0,tdecay)
character*10::component
integer::i,j
integer::t0,tdecay
real*8::z0,frequency,amp,phase
integer:: k           !// at 3*tdecay Gaussian=0.001 !// 
	k=floor(0.5+((z0+zcenter)*lattice_z)) 
	if(component=='Ex') then
		Ex(i,j,k)=Ex(i,j,k)+amp*Gauss_amp(frequency,phase,t0,tdecay)/2
		Ex(i-1,j,k)=Ex(i-1,j,k)+amp*Gauss_amp(frequency,phase,t0,tdecay)/2 
       endif
       if(component=="-Ex") then
        
       		Ex(i,j,k)=Ex(i,j,k)+amp*Gauss_amp(frequency,phase+pi,t0,tdecay)/2 
		Ex(i-1,j,k)=Ex(i-1,j,k)+amp*Gauss_amp(frequency,phase+pi,t0,tdecay)/2 
       endif
       if(component=="Ey") then
        
       		Ey(i,j,k)=Ey(i,j,k)+amp*Gauss_amp(frequency,phase,t0,tdecay)/2 
		Ey(i,j-1,k)=Ey(i,j-1,k)+amp*Gauss_amp(frequency,phase,t0,tdecay)/2 
       endif
       if(component=="-Ey") then
        
       		Ey(i,j,k)=Ey(i,j,k)+amp*Gauss_amp(frequency,phase+pi,t0,tdecay)/2 
		Ey(i,j-1,k)=Ey(i,j-1,k)+amp*Gauss_amp(frequency,phase+pi,t0,tdecay)/2 
       endif
       if(component=="Ez") then
        
       		Ez(i,j,k)=Ez(i,j,k)+amp*Gauss_amp(frequency,phase,t0,tdecay)/2 
		Ez(i,j,k-1)=Ez(i,j,k-1)+amp*Gauss_amp(frequency,phase,t0,tdecay)/2 
       endif
       if(component=="-Ez") then
        
       		Ez(i,j,k)=Ez(i,j,k)+amp*Gauss_amp(frequency,phase+pi,t0,tdecay)/2 
		Ez(i,j,k-1)=Ez(i,j,k-1)+amp*Gauss_amp(frequency,phase+pi,t0,tdecay)/2 
       endif
       if(component=="Hx") then
        
       		Hx(i,j,k)=Hx(i,j,k)+amp*Gauss_amp(frequency,phase,t0,tdecay)*sqrt(eo/uo)/4 
		Hx(i,j-1,k)=Hx(i,j-1,k)+amp*Gauss_amp(frequency,phase,t0,tdecay)*sqrt(eo/uo)/4 
		Hx(i,j,k-1)=Hx(i,j,k-1)+amp*Gauss_amp(frequency,phase,t0,tdecay)*sqrt(eo/uo)/4 
		Hx(i,j-1,k-1)=Hx(i,j-1,k-1)+amp*Gauss_amp(frequency,phase,t0,tdecay)*sqrt(eo/uo)/4 
      endif

       if(component=="-Hx") then
        
       		Hx(i,j,k)=Hx(i,j,k)+amp*Gauss_amp(frequency,phase+pi,t0,tdecay)*sqrt(eo/uo)/4 
		Hx(i,j-1,k)=Hx(i,j-1,k)+amp*Gauss_amp(frequency,phase+pi,t0,tdecay)*sqrt(eo/uo)/4 
		Hx(i,j,k-1)=Hx(i,j,k-1)+amp*Gauss_amp(frequency,phase+pi,t0,tdecay)*sqrt(eo/uo)/4 
		Hx(i,j-1,k-1)=Hx(i,j-1,k-1)+amp*Gauss_amp(frequency,phase+pi,t0,tdecay)*sqrt(eo/uo)/4 
       endif
       if(component=="Hy") then
        
		Hy(i,j,k)=Hy(i,j,k)+amp*Gauss_amp(frequency,phase,t0,tdecay)*sqrt(eo/uo)/4 
		Hy(i,j,k-1)=Hy(i,j,k-1)+amp*Gauss_amp(frequency,phase,t0,tdecay)*sqrt(eo/uo)/4 
		Hy(i-1,j,k)=Hy(i-1,j,k)+amp*Gauss_amp(frequency,phase,t0,tdecay)*sqrt(eo/uo)/4 
		Hy(i-1,j,k-1)=Hy(i-1,j,k-1)+amp*Gauss_amp(frequency,phase,t0,tdecay)*sqrt(eo/uo)/4 
       endif!
       if(component=="-Hy") then
        
       		Hy(i,j,k)=Hy(i,j,k)+amp*Gauss_amp(frequency,phase+pi,t0,tdecay)*sqrt(eo/uo)/4 
		Hy(i,j,k-1)=Hy(i,j,k-1)+amp*Gauss_amp(frequency,phase+pi,t0,tdecay)*sqrt(eo/uo)/4 
		Hy(i-1,j,k)=Hy(i-1,j,k)+amp*Gauss_amp(frequency,phase+pi,t0,tdecay)*sqrt(eo/uo)/4 
		Hy(i-1,j,k-1)=Hy(i-1,j,k-1)+amp*Gauss_amp(frequency,phase+pi,t0,tdecay)*sqrt(eo/uo)/4 
       endif
       if(component=="Hz") then
        
       		Hz(i,j,k)=Hz(i,j,k)+amp*Gauss_amp(frequency,phase,t0,tdecay)*sqrt(eo/uo)/4 
		Hz(i-1,j,k)=Hz(i-1,j,k)+amp*Gauss_amp(frequency,phase,t0,tdecay)*sqrt(eo/uo)/4 
		Hz(i,j-1,k)=Hz(i,j-1,k)+amp*Gauss_amp(frequency,phase,t0,tdecay)*sqrt(eo/uo)/4 
		Hz(i-1,j-1,k)=Hz(i-1,j-1,k)+amp*Gauss_amp(frequency,phase,t0,tdecay)*sqrt(eo/uo)/4 
       	endif
       if(component=="-Hz") then
        
       		Hz(i,j,k)=Hz(i,j,k)+amp*Gauss_amp(frequency,phase+pi,t0,tdecay)*sqrt(eo/uo)/4 
		Hz(i-1,j,k)=Hz(i-1,j,k)+amp*Gauss_amp(frequency,phase+pi,t0,tdecay)*sqrt(eo/uo)/4
		Hz(i,j-1,k)=Hz(i,j-1,k)+amp*Gauss_amp(frequency,phase+pi,t0,tdecay)*sqrt(eo/uo)/4 
		Hz(i-1,j-1,k)=Hz(i-1,j-1,k)+amp*Gauss_amp(frequency,phase+pi,t0,tdecay)*sqrt(eo/uo)/4 
      endif
end subroutine Gaussian_dipole_source
subroutine Gaussian_dipole_source_beam(component,i,j,z0,frequency,amp,phase,t0,tdecay)
character*10::component
integer::i,j
integer::t0,tdecay
real*8::z0,frequency,amp,phase
integer:: k           !// at 3*tdecay Gaussian=0.001 !// 
	k=floor(0.5+((z0+zcenter)*lattice_z)) 
	if(component=='Ex') then
		Ex(i,j,k)=Ex(i,j,k)+amp*Gauss_amp(frequency,phase,t0,tdecay)/2
	!	if (i==201) then
	!	print *,i,j,k,'Ex',Ex(i,j,k)
	!	endif
		if(i>=1) then
			Ex(i-1,j,k)=Ex(i-1,j,k)+amp*Gauss_amp(frequency,phase,t0,tdecay)/2 
		endif
       endif
       if(component=="-Ex") then
        
       	Ex(i,j,k)=Ex(i,j,k)+amp*Gauss_amp(frequency,phase+pi,t0,tdecay)/2 
		if(i>=1) then
			Ex(i-1,j,k)=Ex(i-1,j,k)+amp*Gauss_amp(frequency,phase+pi,t0,tdecay)/2 
		endif
       endif
       if(component=="Ey") then
        
       	Ey(i,j,k)=Ey(i,j,k)+amp*Gauss_amp(frequency,phase,t0,tdecay)/2 
		if(j>=1) then
			Ey(i,j-1,k)=Ey(i,j-1,k)+amp*Gauss_amp(frequency,phase,t0,tdecay)/2 
		endif
       endif
       if(component=="-Ey") then
        
       	Ey(i,j,k)=Ey(i,j,k)+amp*Gauss_amp(frequency,phase+pi,t0,tdecay)/2 
		if(j>=1) then
			Ey(i,j-1,k)=Ey(i,j-1,k)+amp*Gauss_amp(frequency,phase+pi,t0,tdecay)/2 
		endif
       endif
       if(component=="Ez") then
        
       	Ez(i,j,k)=Ez(i,j,k)+amp*Gauss_amp(frequency,phase,t0,tdecay)/2 
		if(k>=1) then
			Ez(i,j,k-1)=Ez(i,j,k-1)+amp*Gauss_amp(frequency,phase,t0,tdecay)/2 
		endif
       endif
       if(component=="-Ez") then
        
       	Ez(i,j,k)=Ez(i,j,k)+amp*Gauss_amp(frequency,phase+pi,t0,tdecay)/2 
		if(k>=1) then
			Ez(i,j,k-1)=Ez(i,j,k-1)+amp*Gauss_amp(frequency,phase+pi,t0,tdecay)/2 
		endif
       endif
       if(component=="Hx") then
        
       	Hx(i,j,k)=Hx(i,j,k)+amp*Gauss_amp(frequency,phase,t0,tdecay)*sqrt(eo/uo)/4 
		if(j>=1) then
			Hx(i,j-1,k)=Hx(i,j-1,k)+amp*Gauss_amp(frequency,phase,t0,tdecay)*sqrt(eo/uo)/4 
		endif
		if(k>=1) then
			Hx(i,j,k-1)=Hx(i,j,k-1)+amp*Gauss_amp(frequency,phase,t0,tdecay)*sqrt(eo/uo)/4 
		endif
		if(j>=1 .and. k>=1) then
			Hx(i,j-1,k-1)=Hx(i,j-1,k-1)+amp*Gauss_amp(frequency,phase,t0,tdecay)*sqrt(eo/uo)/4 
		endif
       endif
       if(component=="-Hx") then
        
       	Hx(i,j,k)=Hx(i,j,k)+amp*Gauss_amp(frequency,phase+pi,t0,tdecay)*sqrt(eo/uo)/4 
		if(j>=1) then
			Hx(i,j-1,k)=Hx(i,j-1,k)+amp*Gauss_amp(frequency,phase+pi,t0,tdecay)*sqrt(eo/uo)/4 
		endif
		if(k>=1) then
			Hx(i,j,k-1)=Hx(i,j,k-1)+amp*Gauss_amp(frequency,phase+pi,t0,tdecay)*sqrt(eo/uo)/4 
		endif
		if(j>=1 .and. k>=1) then
			Hx(i,j-1,k-1)=Hx(i,j-1,k-1)+amp*Gauss_amp(frequency,phase+pi,t0,tdecay)*sqrt(eo/uo)/4 
		endif
       endif
       if(component=="Hy") then
        
		Hy(i,j,k)=Hy(i,j,k)+amp*Gauss_amp(frequency,phase,t0,tdecay)*sqrt(eo/uo)/4 
		if(k>=1) then
			Hy(i,j,k-1)=Hy(i,j,k-1)+amp*Gauss_amp(frequency,phase,t0,tdecay)*sqrt(eo/uo)/4 
		endif
		if(i>=1) then
			Hy(i-1,j,k)=Hy(i-1,j,k)+amp*Gauss_amp(frequency,phase,t0,tdecay)*sqrt(eo/uo)/4 
		endif	
		if(i>=1 .and. k>=1) then
			Hy(i-1,j,k-1)=Hy(i-1,j,k-1)+amp*Gauss_amp(frequency,phase,t0,tdecay)*sqrt(eo/uo)/4 
		endif
       endif
       if(component=="-Hy") then
        
       	Hy(i,j,k)=Hy(i,j,k)+amp*Gauss_amp(frequency,phase+pi,t0,tdecay)*sqrt(eo/uo)/4 
		if(k>=1) then
			Hy(i,j,k-1)=Hy(i,j,k-1)+amp*Gauss_amp(frequency,phase+pi,t0,tdecay)*sqrt(eo/uo)/4 
		endif	
		if(i>=1) then
			Hy(i-1,j,k)=Hy(i-1,j,k)+amp*Gauss_amp(frequency,phase+pi,t0,tdecay)*sqrt(eo/uo)/4 
		endif	
		if(i>=1 .and. k>=1) then
			Hy(i-1,j,k-1)=Hy(i-1,j,k-1)+amp*Gauss_amp(frequency,phase+pi,t0,tdecay)*sqrt(eo/uo)/4 
		endif	
       endif
       if(component=="Hz") then
        
       	Hz(i,j,k)=Hz(i,j,k)+amp*Gauss_amp(frequency,phase,t0,tdecay)*sqrt(eo/uo)/4 
		if(i>=1) then
			Hz(i-1,j,k)=Hz(i-1,j,k)+amp*Gauss_amp(frequency,phase,t0,tdecay)*sqrt(eo/uo)/4 
		endif
		if(j>=1) then
			Hz(i,j-1,k)=Hz(i,j-1,k)+amp*Gauss_amp(frequency,phase,t0,tdecay)*sqrt(eo/uo)/4 
		endif
		if(i>=1 .and. j>=1) then
			Hz(i-1,j-1,k)=Hz(i-1,j-1,k)+amp*Gauss_amp(frequency,phase,t0,tdecay)*sqrt(eo/uo)/4 
		endif
       endif
       if(component=="-Hz") then
        
       	Hz(i,j,k)=Hz(i,j,k)+amp*Gauss_amp(frequency,phase+pi,t0,tdecay)*sqrt(eo/uo)/4 
		if(i>=1) then
			Hz(i-1,j,k)=Hz(i-1,j,k)+amp*Gauss_amp(frequency,phase+pi,t0,tdecay)*sqrt(eo/uo)/4
		endif
		if(j>=1) then
			Hz(i,j-1,k)=Hz(i,j-1,k)+amp*Gauss_amp(frequency,phase+pi,t0,tdecay)*sqrt(eo/uo)/4 
		endif
		if(i>=1 .and. j>=1) then
			Hz(i-1,j-1,k)=Hz(i-1,j-1,k)+amp*Gauss_amp(frequency,phase+pi,t0,tdecay)*sqrt(eo/uo)/4 
		endif
	endif
end subroutine Gaussian_dipole_source_beam


subroutine Gaussian_planewave_TF_SF(rpos,eps_bg,frequency,t0,tdecay)
real*8::rpos(3,2),eps_bg,frequency
integer::t0,tdecay
!//////////////////////////////////////////////////////////////////////////
!// Currently, only the case where {Ecomp = -Ex & Hcomp = Hy} is supported.
!//////////////////////////////////////////////////////////////////////////  

	TF_SF_pos(1,1)= floor(0.5+((rpos(1,1)+xcenter)*lattice_x))
	TF_SF_pos(1,2)= floor(0.5+((rpos(1,2)+xcenter)*lattice_x))
	TF_SF_pos(2,1)= floor(0.5+((rpos(2,1)+ycenter)*lattice_y))
	TF_SF_pos(2,2)= floor(0.5+((rpos(2,2)+ycenter)*lattice_y))
	TF_SF_pos(3,1)= floor(0.5+((rpos(3,1)+zcenter)*lattice_z))
	TF_SF_pos(3,2)= floor(0.5+((rpos(3,2)+zcenter)*lattice_z))
	TF_SF_freq = frequency 
	TF_SF_to = t0 
	TF_SF_tdecay = tdecay
	TF_SF_eps_bg = eps_bg
end subroutine Gaussian_planewave_TF_SF


function Gauss_amp(frequency,phase,t0,tdecay)
real*8::frequency,phase
integer::t0,tdecay
real*8::Gauss_amp
         if(t0-3*tdecay<t .and. t<t0+3*tdecay) then
		if (Gauss_typ(1:8)=='harmonic') then
			   Gauss_amp=sin(2*pi*frequency*t/S_factor/ds_x/lattice_x+phase)*exp(-1.0*((dble(t-t0)/tdecay)**2.0))
	        elseif (Gauss_typ(1:5)=='pulse') then
			   Gauss_amp=exp(-1.0*((dble(t-t0)/tdecay)**2.0))
		endif
         else
                           Gauss_amp=0.0
         endif


end function Gauss_amp

end module source
