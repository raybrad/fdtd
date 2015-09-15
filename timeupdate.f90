module timeupdate
use parameters

!The Yee algorithm(Taflove I 59-) and UPML(Taflove II 2.075-)
implicit none
real*8::tus

contains


subroutine coefficient()
integer::i,j,k


do  k=0,mksize-1
     do j=0,mjsize-1
          do i=0,misize-1
		if  (cposition(i,j,k)==1) then
		    continue      !!//Do nothing!
		
		
		else
		             ! a b c d e f-- E 
			aax(i)=(2.0*eo*kappax(dble(i))-dt*sigmax(dble(i)))/(2.0*eo*kappax(dble(i))+dt*sigmax(dble(i))) 
			aay(j)=(2.0*eo*kappay(dble(j))-dt*sigmay(dble(j)))/(2.0*eo*kappay(dble(j))+dt*sigmay(dble(j))) 
			aaz(k)=(2.0*eo*kappaz(dble(k))-dt*sigmaz(dble(k)))/(2.0*eo*kappaz(dble(k))+dt*sigmaz(dble(k))) 

			bbx(i)=(2.0*eo*dt)/(2.0*eo*kappax(dble(i))+dt*sigmax(dble(i))) 
			bby(j)=(2.0*eo*dt)/(2.0*eo*kappay(dble(j))+dt*sigmay(dble(j))) 
			bbz(k)=(2.0*eo*dt)/(2.0*eo*kappaz(dble(k))+dt*sigmaz(dble(k))) 

			ccx(i)=(2.0*eo*kappax(dble(i))-dt*sigmax(dble(i)))/(2.0*eo*kappax(dble(i))+dt*sigmax(dble(i))) 
			ccy(j)=(2.0*eo*kappay(dble(j))-dt*sigmay(dble(j)))/(2.0*eo*kappay(dble(j))+dt*sigmay(dble(j))) 
			ccz(k)=(2.0*eo*kappaz(dble(k))-dt*sigmaz(dble(k)))/(2.0*eo*kappaz(dble(k))+dt*sigmaz(dble(k))) 

			ddx(i,j,k)=2.0/(2.0*eo*kappaz(dble(k))+dt*sigmaz(dble(k)))/epsilonx(i,j,k) 
			ddy(i,j,k)=2.0/(2.0*eo*kappax(dble(i))+dt*sigmax(dble(i)))/epsilony(i,j,k) 
			ddz(i,j,k)=2.0/(2.0*eo*kappay(dble(j))+dt*sigmay(dble(j)))/epsilonz(i,j,k) 

			eex(i)=kappax(dble(i+0.5))+dt*sigmax(dble(i+0.5))/2.0/eo 
			eey(j)=kappay(dble(j+0.5))+dt*sigmay(dble(j+0.5))/2.0/eo 
			eez(k)=kappaz(dble(k+0.5))+dt*sigmaz(dble(k+0.5))/2.0/eo 
                                !  g h i j k l-- H
			ffx(i)=kappax(dble(i+0.5))-dt*sigmax(dble(i+0.5))/2.0/eo 
			ffy(j)=kappay(dble(j+0.5))-dt*sigmay(dble(j+0.5))/2.0/eo 
			ffz(k)=kappaz(dble(k+0.5))-dt*sigmaz(dble(k+0.5))/2.0/eo 

			ggx(i)=(2.0*eo*kappax(dble(i+0.5))-dt*sigmax(dble(i+0.5)))/(2.0*eo*kappax(dble(i+0.5))+dt*sigmax(dble(i+0.5))) 
			ggy(j)=(2.0*eo*kappay(dble(j+0.5))-dt*sigmay(dble(j+0.5)))/(2.0*eo*kappay(dble(j+0.5))+dt*sigmay(dble(j+0.5))) 
			ggz(k)=(2.0*eo*kappaz(dble(k+0.5))-dt*sigmaz(dble(k+0.5)))/(2.0*eo*kappaz(dble(k+0.5))+dt*sigmaz(dble(k+0.5))) 
			
			hhx(i)=(2.0*eo*dt)/(2.0*eo*kappax(dble(i+0.5))+dt*sigmax(dble(i+0.5))) 
			hhy(j)=(2.0*eo*dt)/(2.0*eo*kappay(dble(j+0.5))+dt*sigmay(dble(j+0.5))) 
			hhz(k)=(2.0*eo*dt)/(2.0*eo*kappaz(dble(k+0.5))+dt*sigmaz(dble(k+0.5))) 

			iix(i)=(2.0*eo*kappax(dble(i+0.5))-dt*sigmax(dble(i+0.5)))/(2.0*eo*kappax(dble(i+0.5))+dt*sigmax(dble(i+0.5))) 
			iiy(j)=(2.0*eo*kappay(dble(j+0.5))-dt*sigmay(dble(j+0.5)))/(2.0*eo*kappay(dble(j+0.5))+dt*sigmay(dble(j+0.5))) 
			iiz(k)=(2.0*eo*kappaz(dble(k+0.5))-dt*sigmaz(dble(k+0.5)))/(2.0*eo*kappaz(dble(k+0.5))+dt*sigmaz(dble(k+0.5))) 

			jjx(i)=(2.0*eo)/(2.0*eo*kappax(dble(i+0.5))+dt*sigmax(dble(i+0.5)))/uo/ups 
			jjy(j)=(2.0*eo)/(2.0*eo*kappay(dble(j+0.5))+dt*sigmay(dble(j+0.5)))/uo/ups 
			jjz(k)=(2.0*eo)/(2.0*eo*kappaz(dble(k+0.5))+dt*sigmaz(dble(k+0.5)))/uo/ups 

			kkx(i)=kappax(dble(i))+dt*sigmax(dble(i))/2.0/eo 
			kky(j)=kappay(dble(j))+dt*sigmay(dble(j))/2.0/eo 
			kkz(k)=kappaz(dble(k))+dt*sigmaz(dble(k))/2.0/eo 

			llx(i)=kappax(dble(i))-dt*sigmax(dble(i))/2.0/eo 
			lly(j)=kappay(dble(j))-dt*sigmay(dble(j))/2.0/eo 
			llz(k)=kappaz(dble(k))-dt*sigmaz(dble(k))/2.0/eo
		end if
	enddo
     enddo
enddo
			tus=dt/uo/ups 
	print *,"coefficient...ok" 

end subroutine coefficient

subroutine  propagate()
integer:: i,j,k
real*8::tba,t1s,t1e,t2s,t2e,t3s,t3e,t4s,t4e,t5s,t5e,t6s,t6e
!	!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//
!	!// Note 
!	!// 
!	!// 1) H-field  !! 
!	!// 2) TF-SF boundary for H
!	!// 3) E-field  !!
!	!// 4) TF-SF boundary for E
!	!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//

!	!//******* 1) real H-field ******!//
!So the default boundary is Dirichit ,since we don't update them 
!plane x=0,y=0,z=0 		  has E=0 H=0
!plane x=pisize,y=pjsize,z=pksize has E=0 H=0

do k=1,pksize-1
    do j=1,pjsize-1
          do i=1,pisize-1
		  if (cposition(i,j,k)==0) then ! PML  
		      call  pml_H_field_update(i, j, k) 
		  elseif (cposition(i,j,k)==1) then  !!//Non- PML (metal or dielectric air) 
		      call  normal_H_field_update(i, j, k)
		  else
		     continue
		  endif
         enddo
     enddo	
enddo

	
		

	!!//******* 2) TF-SF for H ******!//
	if(TF_SF_freq /= 0.0) then
	     write(6,*) 'in TF_SF calculation,step:',t
	     call H_inc_field_update_for_TF_SF() 
	     call TF_SF_correction_for_H() 
	endif


	!!//******* 3) real E-field ******!//
	
do k=1,pksize-1
    do j=1,pjsize-1
          do i=1,pisize-1
			if(cposition(i,j,k)==0) then  !// PML 
				call  pml_E_field_update(i, j, k) 
			elseif(cposition(i,j,k)==1) then!Non-PML
			   if(Lepsilon(i,j,k)/=0.0 ) then  ! Metal (Leps /=0 /=1000)
	                        if (metal_type==1) then    !normal metal
				 call  Normal_metal_E_field_update(i,j,k)
	                        elseif (metal_type==2 ) then       !drude metal
			 	 call  Drude_metal_E_field_update(i, j, k)
	                        elseif (metal_type==3 ) then       !single Lorentz metal
				 call  Lorentz_E_field_update(i,j,k)
	                        elseif (metal_type==4 ) then       !double Lorentz metal
				 call  Double_Lorentz_E_field_update(i,j,k)
	                        elseif (metal_type==5 ) then       !drude+Double Lorentz pole metal
			 	 call  Drude_double_Lorentz_E_field_update(i, j, k)
	                        endif
			   else				       !Normal dielectric material(air) (Leps =0.0)  no metal
	                       if ( Molecular(i,j,k)==1.0 ) then                                
			        call M_normal_E_field_update(i,j,k)
	                       else
				call  normal_E_field_update(i, j, k) 
	                       endif
			   endif
	                    
			else
				continue
                        endif
	 enddo
     enddo
enddo

	!!//******* 9) TF-SF for E  ******!//
	if(TF_SF_freq /= 0.0) then 
	
		call E_inc_field_update_for_TF_SF() 
		call TF_SF_correction_for_E() 
	endif
!Fourier transformed field

if(dofft) then
do k=1,pksize-1
    do j=1,pjsize-1
          do i=1,pisize-1
		    !if(lmonitor(i,j,k)) then
			call fft_field(i,j,k)	
	        !endif
	  enddo
   enddo
enddo
call fft_field_inc
endif
end subroutine propagate
			

subroutine  normal_E_field_update( i,  j,  k)
integer:: i,j,k

  Ex(i,j,k) = Ex(i,j,k) + (dt/eo/epsilonx(i,j,k))*((Hz(i,j,k)-Hz(i,j-1,k))/ds_y - (Hy(i,j,k)-Hy(i,j,k-1))/ds_z) 	
  Ey(i,j,k) = Ey(i,j,k) + (dt/eo/epsilony(i,j,k))*((Hx(i,j,k)-Hx(i,j,k-1))/ds_z - (Hz(i,j,k)-Hz(i-1,j,k))/ds_x)
  Ez(i,j,k) = Ez(i,j,k) + (dt/eo/epsilonz(i,j,k))*((Hy(i,j,k)-Hy(i-1,j,k))/ds_x - (Hx(i,j,k)-Hx(i,j-1,k))/ds_y)


end subroutine normal_E_field_update

subroutine Normal_metal_E_field_update(i,j,k)
integer::i,j,k
!dt_R/dt_F = lattice_n/(ds_x*lattice_x)
		!unit
!eo*epsilon     F/m = A/V/m*s
!sigma*dt	A/V/m*[lattice_n/(ds_x*lattice_x)*dt] =A/V/m*s
!
  Ex(i,j,k) = (2*eo*epsilonx(i,j,k)-Lsigma(i,j,k)*dt)/(2*eo*epsilonx(i,j,k)+Lsigma(i,j,k)*dt)*Ex(i,j,k) + &
    (2*dt/(2*eo*epsilonx(i,j,k)+Lsigma(i,j,k)*dt))*((Hz(i,j,k)-Hz(i,j-1,k))/ds_y - (Hy(i,j,k)-Hy(i,j,k-1))/ds_z) 	
  Ey(i,j,k) = (2*eo*epsilony(i,j,k)-Lsigma(i,j,k)*dt)/(2*eo*epsilony(i,j,k)+Lsigma(i,j,k)*dt)*Ey(i,j,k) + &
    (2*dt/(2*eo*epsilony(i,j,k)+Lsigma(i,j,k)*dt))*((Hx(i,j,k)-Hx(i,j,k-1))/ds_z - (Hz(i,j,k)-Hz(i-1,j,k))/ds_x)
  Ez(i,j,k) = (2*eo*epsilonz(i,j,k)-Lsigma(i,j,k)*dt)/(2*eo*epsilonz(i,j,k)+Lsigma(i,j,k)*dt)*Ez(i,j,k) + &
    (2*dt/(2*eo*epsilonz(i,j,k)+Lsigma(i,j,k)*dt))*((Hy(i,j,k)-Hy(i-1,j,k))/ds_x - (Hx(i,j,k)-Hx(i,j-1,k))/ds_y)
end subroutine Normal_metal_E_field_update

subroutine M_normal_E_field_update( i,  j,  k)
integer:: i,j,k
logical::iexist
real*8:: prefac,qmCurrx,qmCurry,qmCurrz
!!

  qmCurrx=0.d0
  qmCurry=0.d0
  qmCurrz=0.d0
if ( lqmtoem ) then
  inquire(file='fdqmbound',exist=iexist)
  if (iexist) then
  open(156,file='fdqmbound',status='unknown',form='formatted')
  read(156,'(3(E20.5,X))') qmCurrx,qmCurry,qmCurrz
  close(156)
!since the data from qm is e*bohr/fs, should change to el_charge*bohr*1D15 (C*m/s)
!here we choose the unit cell to be  q^3 ~(a_m/lattice_x)^3 ,a_m=l_n
!!then Jx(C/m^2/s)= (C*m/s)* lattice_x^3 /(a_m)^3    
!!     Jx(C/q^2/s)= (C*m/s)* 1q/(l_n/lattice_x) /q^3    
!!     Jx(C/q/m/s)= (C*m/s)* 1q/(l_n/lattice_x) /q^3  /((l_n/lattice_x)/1q)
!                 = (C*m/s)* (lattice_x/l_n)^2 /q
!volfactor~q^3

  prefac=el_charge*bohr*1D15*(dble(lattice_x)/l_n)**2/(volfactor)**3
  qmCurrx=qmCurrx*prefac
  qmCurry=qmCurry*prefac
  qmCurrz=qmCurrz*prefac
  write(6,*) 'prefac:',prefac
  write(6,'(A,X,3(E20.5,X))') 'qmCurr',qmCurrx,qmCurry,qmCurrz
  else
	  write(6,*) 'fdqmbound not found yet'
  endif
endif
  prefac=el_charge*bohr*1D15*(dble(lattice_x)/l_n)**2/(volfactor)**3
!  dt(1/lightspeed/S(q^-1)  (sec*q*m^-1) =>
!or dt(fs)= a_m*1E15/(C*S*ds_x*lattice_x)= dt(sec*q*m^-1)  *a_m*1E15/(ds_x*lattice_x)
!in this fashion Ex~ V/m , H~ A/m
!ds_y,ds_z,ds_x ~q ;    dt/eo/epsilonx ~ s*q/m /(C/V/m) ~ s*q*V/C
  Ex(i,j,k) = Ex(i,j,k) + (dt/eo/epsilonx(i,j,k))*((Hz(i,j,k)-Hz(i,j-1,k))/ds_y - (Hy(i,j,k)-Hy(i,j,k-1))/ds_z-qmCurrx) 	
  Ey(i,j,k) = Ey(i,j,k) + (dt/eo/epsilony(i,j,k))*((Hx(i,j,k)-Hx(i,j,k-1))/ds_z - (Hz(i,j,k)-Hz(i-1,j,k))/ds_x-qmCurry)
  Ez(i,j,k) = Ez(i,j,k) + (dt/eo/epsilonz(i,j,k))*((Hy(i,j,k)-Hy(i-1,j,k))/ds_x - (Hx(i,j,k)-Hx(i,j-1,k))/ds_y-qmCurrz) 
end subroutine M_normal_E_field_update


subroutine  Drude_metal_E_field_update(i, j, k)
integer::i,j,k
real*8:: km, bm  
real*8:: Ex_temp, Ey_temp, Ez_temp  
	km = (2-Lgamma(i,j,k)*dt)/(2+Lgamma(i,j,k)*dt) 
	bm = Lomega(i,j,k)*Lomega(i,j,k)*eo*dt/(2+Lgamma(i,j,k)*dt) 
	
	Ex_temp = Ex(i,j,k)   !// store E-field at FDTD time 'n' 
	Ey_temp = Ey(i,j,k)   !//          "               
	Ez_temp = Ez(i,j,k)   !//          "

	Ex(i,j,k) = ((2*eo*Lepsilon(i,j,k)-dt*bm)/(2*eo*Lepsilon(i,j,k)+dt*bm))*Ex(i,j,k) + (2*dt/(2*eo*Lepsilon(i,j,k)+dt*bm))*((Hz(i,j,k)-Hz(i,j-1,k))/ds_y - (Hy(i,j,k)-Hy(i,j,k-1))/ds_z - 0.5*(1+km)*Jx(i,j,k)) 
	Ey(i,j,k) = ((2*eo*Lepsilon(i,j,k)-dt*bm)/(2*eo*Lepsilon(i,j,k)+dt*bm))*Ey(i,j,k) + (2*dt/(2*eo*Lepsilon(i,j,k)+dt*bm))*((Hx(i,j,k)-Hx(i,j,k-1))/ds_z - (Hz(i,j,k)-Hz(i-1,j,k))/ds_x - 0.5*(1+km)*Jy(i,j,k)) 
	Ez(i,j,k) = ((2*eo*Lepsilon(i,j,k)-dt*bm)/(2*eo*Lepsilon(i,j,k)+dt*bm))*Ez(i,j,k) + (2*dt/(2*eo*Lepsilon(i,j,k)+dt*bm))*((Hy(i,j,k)-Hy(i-1,j,k))/ds_x - (Hx(i,j,k)-Hx(i,j-1,k))/ds_y - 0.5*(1+km)*Jz(i,j,k)) 

	Jx(i,j,k) = km*Jx(i,j,k) + bm*(Ex(i,j,k) + Ex_temp) 
	Jy(i,j,k) = km*Jy(i,j,k) + bm*(Ey(i,j,k) + Ey_temp) 
	Jz(i,j,k) = km*Jz(i,j,k) + bm*(Ez(i,j,k) + Ez_temp) 

end subroutine Drude_metal_E_field_update

subroutine Lorentz_E_field_update(i,j,k)
integer::i,j,k
real*8:: am, bm, gm, dm  
real*8:: Ex_temp, Ey_temp, Ez_temp  
real*8::Jx_temp, Jy_temp, Jz_temp 	

	am = (2 - Lomega_1(i,j,k)*Lomega_1(i,j,k)*dt*dt)/(Lgamma_1(i,j,k)*dt+1) 
	bm = (Lgamma_1(i,j,k)*dt-1)/(Lgamma_1(i,j,k)*dt+1) 
	gm = eo*Ldepsr_1(i,j,k)*Lomega_1(i,j,k)*Lomega_1(i,j,k)*dt*dt/(Lgamma_1(i,j,k)*dt+1) 
	dm = (2 *eo* Lepsilon(i,j,k) + 0.5*gm) 

	Ex_temp = Ex(i,j,k)   !// store E-field at FDTD time 'n' 
	Ey_temp = Ey(i,j,k)   !//          "               
	Ez_temp = Ez(i,j,k)   !//          "

	Jx_temp = Jx(i,j,k)   !// store J-field at FDTD time 'n' 
	Jy_temp = Jy(i,j,k)   !//          "               
	Jz_temp = Jz(i,j,k)   !//          "

	Ex(i,j,k) = (2*eo*Lepsilon(i,j,k)/dm)*Ex(i,j,k) + (0.5*gm/dm)*pEx(i,j,k) - ((am+1)*dt/dm)*Jx(i,j,k) - (bm*dt/dm)*pJx(i,j,k) + (2*dt/dm)*((Hz(i,j,k)-Hz(i,j-1,k))/ds_y - (Hy(i,j,k)-Hy(i,j,k-1))/ds_z) 
	Ey(i,j,k) = (2*eo*Lepsilon(i,j,k)/dm)*Ey(i,j,k) + (0.5*gm/dm)*pEy(i,j,k) - ((am+1)*dt/dm)*Jy(i,j,k) - (bm*dt/dm)*pJy(i,j,k) + (2*dt/dm)*((Hx(i,j,k)-Hx(i,j,k-1))/ds_z - (Hz(i,j,k)-Hz(i-1,j,k))/ds_x) 
	Ez(i,j,k) = (2*eo*Lepsilon(i,j,k)/dm)*Ez(i,j,k) + (0.5*gm/dm)*pEz(i,j,k) - ((am+1)*dt/dm)*Jz(i,j,k) - (bm*dt/dm)*pJz(i,j,k) + (2*dt/dm)*((Hy(i,j,k)-Hy(i-1,j,k))/ds_x - (Hx(i,j,k)-Hx(i,j-1,k))/ds_y) 
	
	Jx(i,j,k) = am*Jx(i,j,k) + bm*pJx(i,j,k)+ gm*(Ex(i,j,k) - pEx(i,j,k))/(2.0*dt) 
	Jy(i,j,k) = am*Jy(i,j,k) + bm*pJy(i,j,k)+ gm*(Ey(i,j,k) - pEy(i,j,k))/(2.0*dt) 
	Jz(i,j,k) = am*Jz(i,j,k) + bm*pJz(i,j,k)+ gm*(Ez(i,j,k) - pEz(i,j,k))/(2.0*dt) 

	pEx(i,j,k) = Ex_temp 
	pEy(i,j,k) = Ey_temp 
	pEz(i,j,k) = Ez_temp 

	pJx(i,j,k) = Jx_temp 
	pJy(i,j,k) = Jy_temp 
	pJz(i,j,k) = Jz_temp 
end subroutine Lorentz_E_field_update

subroutine Double_Lorentz_E_field_update(  i,   j,   k)
integer::i,j,k         
real*8:: Ex_temp, Ey_temp, Ez_temp  
real*8:: am_1, bm_1, gm_1, am_2, bm_2, gm_2,dm
real*8:: Jx_temp_1, Jy_temp_1, Jz_temp_1, Jx_temp_2, Jy_temp_2, Jz_temp_2	

	am_1 = ( 2 - Lomega_1(i,j,k)*Lomega_1(i,j,k)*dt*dt)/( 1 + dt*Lgamma_1(i,j,k))  
	bm_1 = (dt*Lgamma_1(i,j,k) - 1)/(dt*Lgamma_1(i,j,k) + 1) 
	gm_1 = eo*Ldepsr_1(i,j,k)*Lomega_1(i,j,k)*Lomega_1(i,j,k)*dt*dt/(1 + dt*Lgamma_1(i,j,k)) 

	am_2 = ( 2 - Lomega_2(i,j,k)*Lomega_2(i,j,k)*dt*dt)/( 1 + dt*Lgamma_2(i,j,k))  
	bm_2 = (dt*Lgamma_2(i,j,k) - 1)/ (dt*Lgamma_2(i,j,k) + 1) 
	gm_2 = eo*Ldepsr_2(i,j,k)*Lomega_2(i,j,k)*Lomega_2(i,j,k)*dt*dt/ (1 + dt*Lgamma_2(i,j,k)) 

	dm = (2*eo* Lepsilon(i,j,k) + gm_1/2+gm_2/2) 

	Ex_temp = Ex(i,j,k)   !// store E-field at FDTD time 'n' 
	Ey_temp = Ey(i,j,k)   !//          "               
	Ez_temp = Ez(i,j,k)   !//          "

	Jx_temp_1 = Jx_1(i,j,k)   !// store Q-field at FDTD time 'n' 
	Jy_temp_1 = Jy_1(i,j,k)   !//          "               
	Jz_temp_1 = Jz_1(i,j,k)   !//          "

	Jx_temp_2 = Jx_2(i,j,k)   !// store Q-field at FDTD time 'n' 
	Jy_temp_2 = Jy_2(i,j,k)   !//          "               
	Jz_temp_2 = Jz_2(i,j,k)   !//          "


	Ex(i,j,k) = ((2*eo*Lepsilon(i,j,k))/dm)*Ex(i,j,k)  &
		  + ((gm_1/2+gm_2/2)/dm)*pEx(i,j,k)  &
		  - ((am_1+1)*dt/dm)*Jx_1(i,j,k)-((am_2+1)*dt/dm)*Jx_2(i,j,k) &
		  - (bm_1*dt/dm)*pJx_1(i,j,k)-(bm_2*dt/dm)*pJx_2(i,j,k) &
		  + (2*dt/dm)*((Hz(i,j,k)-Hz(i,j-1,k))/ds_y - (Hy(i,j,k)-Hy(i,j,k-1))/ds_z) 

	Ey(i,j,k) = ((2*eo*Lepsilon(i,j,k))/dm)*Ey(i,j,k)  &
		  + ((gm_1/2+gm_2/2)/dm)*pEy(i,j,k)  &
		  - ((am_1+1)*dt/dm)*Jy_1(i,j,k)-((am_2+1)*dt/dm)*Jy_2(i,j,k) &
		  - (bm_1*dt/dm)*pJy_1(i,j,k)-(bm_2*dt/dm)*pJy_2(i,j,k) &
		  + (2*dt/dm)*((Hx(i,j,k)-Hx(i,j,k-1))/ds_z - (Hz(i,j,k)-Hz(i-1,j,k))/ds_x)

	Ez(i,j,k) = ((2*eo*Lepsilon(i,j,k))/dm)*Ez(i,j,k)  &
		  + ((gm_1/2+gm_2/2)/dm)*pEz(i,j,k)  &
		  - ((am_1+1)*dt/dm)*Jz_1(i,j,k)-((am_2+1)*dt/dm)*Jz_2(i,j,k) &
		  - (bm_1*dt/dm)*pJz_1(i,j,k)-(bm_2*dt/dm)*pJz_2(i,j,k) &
		  + (2*dt/dm)*((Hy(i,j,k)-Hy(i-1,j,k))/ds_x - (Hx(i,j,k)-Hx(i,j-1,k))/ds_y)


	Jx_1(i,j,k) = am_1*Jx_1(i,j,k) + bm_1* pJx_1(i,j,k)+ gm_1*(Ex(i,j,k)-pEx(i,j,k))/(2*dt)
	Jy_1(i,j,k) = am_1*Jy_1(i,j,k) + bm_1* pJy_1(i,j,k)+ gm_1*(Ey(i,j,k)-pEy(i,j,k))/(2*dt)
	Jz_1(i,j,k) = am_1*Jz_1(i,j,k) + bm_1* pJz_1(i,j,k)+ gm_1*(Ez(i,j,k)-pEz(i,j,k))/(2*dt)

	Jx_2(i,j,k) = am_2*Jx_2(i,j,k) + bm_2* pJx_2(i,j,k)+ gm_2*(Ex(i,j,k)-pEx(i,j,k))/(2*dt)
	Jy_2(i,j,k) = am_2*Jy_2(i,j,k) + bm_2* pJy_2(i,j,k)+ gm_2*(Ey(i,j,k)-pEy(i,j,k))/(2*dt)
	Jz_2(i,j,k) = am_2*Jz_2(i,j,k) + bm_2* pJz_2(i,j,k)+ gm_2*(Ez(i,j,k)-pEz(i,j,k))/(2*dt)

	pEx(i,j,k) = Ex_temp 
	pEy(i,j,k) = Ey_temp 
	pEz(i,j,k) = Ez_temp 

	pJx_1(i,j,k) = Jx_temp_1 
	pJy_1(i,j,k) = Jy_temp_1 
	pJz_1(i,j,k) = Jz_temp_1 

	pJx_2(i,j,k) = Jx_temp_2 
	pJy_2(i,j,k) = Jy_temp_2 
	pJz_2(i,j,k) = Jz_temp_2 

end subroutine Double_Lorentz_E_field_update

subroutine Drude_double_Lorentz_E_field_update(  i,   j,   k)
integer::i,j,k         
real*8:: am_0, gm_0 
real*8:: Ex_temp, Ey_temp, Ez_temp  
real*8:: am_1, bm_1, gm_1, am_2, bm_2, gm_2,dm
real*8:: Jx_temp_1, Jy_temp_1, Jz_temp_1, Jx_temp_2, Jy_temp_2, Jz_temp_2	

	am_0 = (2-Lgamma(i,j,k)*dt)/(Lgamma(i,j,k)*dt+2) !km
	gm_0 = eo*Lomega(i,j,k)*Lomega(i,j,k)*dt/(Lgamma(i,j,k)*dt+2) !bm

	am_1 = ( 2 - Lomega_1(i,j,k)*Lomega_1(i,j,k)*dt*dt)/( 1 + dt*Lgamma_1(i,j,k))  
	bm_1 = (dt*Lgamma_1(i,j,k) - 1)/(dt*Lgamma_1(i,j,k) + 1) 
	gm_1 = eo*Ldepsr_1(i,j,k)*Lomega_1(i,j,k)*Lomega_1(i,j,k)*dt*dt/(1 + dt*Lgamma_1(i,j,k)) 

	am_2 = ( 2 - Lomega_2(i,j,k)*Lomega_2(i,j,k)*dt*dt)/( 1 + dt*Lgamma_2(i,j,k))  
	bm_2 = (dt*Lgamma_2(i,j,k) - 1)/(dt*Lgamma_2(i,j,k) + 1) 
	gm_2 = eo*Ldepsr_2(i,j,k)*Lomega_2(i,j,k)*Lomega_2(i,j,k)*dt*dt/(1 + dt*Lgamma_2(i,j,k)) 

	dm = (2*eo* Lepsilon(i,j,k) + gm_0*dt+gm_1/2+gm_2/2) 
	!write(6,*) 'i,j,k',i,j,k
	!write(6,*) 'dt',dt
	!write(6,*) 'Lgamma,Lomega',Lgamma(i,j,k),Lomega(i,j,k)
	!write(6,*) 'Lgamma1,Lomega1',Lgamma_1(i,j,k),Lomega_1(i,j,k)
	!write(6,*) 'Lgamma2,Lomega2',Lgamma_2(i,j,k),Lomega_2(i,j,k)
	!write(6,*) 'am_0,gm_0',am_0,gm_0
	!write(6,*) 'am_1,bm_1,gm_1',am_1,bm_1,gm_1
	!write(6,*) 'am_2,bm_2,gm_2',am_2,bm_2,gm_2
	!write(6,*) 'dm',dm
	!write(6,*) 'c1',(2*eo*Lepsilon(i,j,k)-gm_0*dt)/dm
	!write(6,*) 'c2',(gm_1/2+gm_2/2)/dm
	!write(6,*) 'c3',((am_1+1)*dt/dm),((am_2+1)*dt/dm),((am_0+1)*dt/dm)
	!write(6,*) 'c4',(bm_1*dt/dm),(bm_2*dt/dm)
	!write(6,*) 'c5',(2*dt/dm)
!if (i==50 .and. j==50 .and. k==50) then	
!	stop
!endif
	Ex_temp = Ex(i,j,k)   !// store E-field at FDTD time 'n' 
	Ey_temp = Ey(i,j,k)   !//          "               
	Ez_temp = Ez(i,j,k)   !//          "

	Jx_temp_1 = Jx_1(i,j,k)   !// store Q-field at FDTD time 'n' 
	Jy_temp_1 = Jy_1(i,j,k)   !//          "               
	Jz_temp_1 = Jz_1(i,j,k)   !//          "

	Jx_temp_2 = Jx_2(i,j,k)   !// store Q-field at FDTD time 'n' 
	Jy_temp_2 = Jy_2(i,j,k)   !//          "               
	Jz_temp_2 = Jz_2(i,j,k)   !//          "


	Ex(i,j,k) = ((2*eo*Lepsilon(i,j,k)-gm_0*dt)/dm)*Ex(i,j,k)  &
		  + ((gm_1/2+gm_2/2)/dm)*pEx(i,j,k)  &
		  - ((am_1+1)*dt/dm)*Jx_1(i,j,k)-((am_2+1)*dt/dm)*Jx_2(i,j,k)-((am_0+1)*dt/dm)*Jx_0(i,j,k) &
		  - (bm_1*dt/dm)*pJx_1(i,j,k)-(bm_2*dt/dm)*pJx_2(i,j,k) &
		  + (2*dt/dm)*((Hz(i,j,k)-Hz(i,j-1,k))/ds_y - (Hy(i,j,k)-Hy(i,j,k-1))/ds_z) 

	Ey(i,j,k) = ((2*eo*Lepsilon(i,j,k)-gm_0*dt)/dm)*Ey(i,j,k)  &
		  + ((gm_1/2+gm_2/2)/dm)*pEy(i,j,k)  &
		  - ((am_1+1)*dt/dm)*Jy_1(i,j,k)-((am_2+1)*dt/dm)*Jy_2(i,j,k)-((am_0+1)*dt/dm)*Jy_0(i,j,k) &
		  - (bm_1*dt/dm)*pJy_1(i,j,k)-(bm_2*dt/dm)*pJy_2(i,j,k) &
		  + (2*dt/dm)*((Hx(i,j,k)-Hx(i,j,k-1))/ds_z - (Hz(i,j,k)-Hz(i-1,j,k))/ds_x)

	Ez(i,j,k) = ((2*eo*Lepsilon(i,j,k)-gm_0*dt)/dm)*Ez(i,j,k)  &
		  + ((gm_1/2+gm_2/2)/dm)*pEz(i,j,k)  &
		  - ((am_1+1)*dt/dm)*Jz_1(i,j,k)-((am_2+1)*dt/dm)*Jz_2(i,j,k)-((am_0+1)*dt/dm)*Jz_0(i,j,k) &
		  - (bm_1*dt/dm)*pJz_1(i,j,k)-(bm_2*dt/dm)*pJz_2(i,j,k) &
		  + (2*dt/dm)*((Hy(i,j,k)-Hy(i-1,j,k))/ds_x - (Hx(i,j,k)-Hx(i,j-1,k))/ds_y)
	
	Jx_0(i,j,k) = am_0*Jx_0(i,j,k) + gm_0*(Ex(i,j,k) + Ex_temp)
	Jy_0(i,j,k) = am_0*Jy_0(i,j,k) + gm_0*(Ey(i,j,k) + Ey_temp)
	Jz_0(i,j,k) = am_0*Jz_0(i,j,k) + gm_0*(Ez(i,j,k) + Ez_temp)


	Jx_1(i,j,k) = am_1*Jx_1(i,j,k) + bm_1* pJx_1(i,j,k)+ gm_1*(Ex(i,j,k)-pEx(i,j,k))/(2*dt)
	Jy_1(i,j,k) = am_1*Jy_1(i,j,k) + bm_1* pJy_1(i,j,k)+ gm_1*(Ey(i,j,k)-pEy(i,j,k))/(2*dt)
	Jz_1(i,j,k) = am_1*Jz_1(i,j,k) + bm_1* pJz_1(i,j,k)+ gm_1*(Ez(i,j,k)-pEz(i,j,k))/(2*dt)

	Jx_2(i,j,k) = am_2*Jx_2(i,j,k) + bm_2* pJx_2(i,j,k)+ gm_2*(Ex(i,j,k)-pEx(i,j,k))/(2*dt)
	Jy_2(i,j,k) = am_2*Jy_2(i,j,k) + bm_2* pJy_2(i,j,k)+ gm_2*(Ey(i,j,k)-pEy(i,j,k))/(2*dt)
	Jz_2(i,j,k) = am_2*Jz_2(i,j,k) + bm_2* pJz_2(i,j,k)+ gm_2*(Ez(i,j,k)-pEz(i,j,k))/(2*dt)

	pEx(i,j,k) = Ex_temp 
	pEy(i,j,k) = Ey_temp 
	pEz(i,j,k) = Ez_temp 

	pJx_1(i,j,k) = Jx_temp_1 
	pJy_1(i,j,k) = Jy_temp_1 
	pJz_1(i,j,k) = Jz_temp_1 

	pJx_2(i,j,k) = Jx_temp_2 
	pJy_2(i,j,k) = Jy_temp_2 
	pJz_2(i,j,k) = Jz_temp_2 

end subroutine Drude_double_Lorentz_E_field_update



subroutine pml_E_field_update(  i,   j,   k)
integer::i,j,k
real*8:: temp 

	temp=Dx(i,j,k) 
	Dx(i,j,k) = aay(j)*Dx(i,j,k) + bby(j)*((Hz(i,j,k)-Hz(i,j-1,k))/ds_y - (Hy(i,j,k)-Hy(i,j,k-1))/ds_z) 
	Ex(i,j,k) = ccz(k)*Ex(i,j,k) + ddx(i,j,k)*(eex(i)*Dx(i,j,k)-ffx(i)*temp) 
	temp=Dy(i,j,k) 
	Dy(i,j,k) = aaz(k)*Dy(i,j,k) + bbz(k)*((Hx(i,j,k)-Hx(i,j,k-1))/ds_z - (Hz(i,j,k)-Hz(i-1,j,k))/ds_x) 
	Ey(i,j,k) = ccx(i)*Ey(i,j,k) + ddy(i,j,k)*(eey(j)*Dy(i,j,k)-ffy(j)*temp) 
	temp=Dz(i,j,k) 
	Dz(i,j,k) = aax(i)*Dz(i,j,k) + bbx(i)*((Hy(i,j,k)-Hy(i-1,j,k))/ds_x - (Hx(i,j,k)-Hx(i,j-1,k))/ds_y) 
	Ez(i,j,k) = ccy(j)*Ez(i,j,k) + ddz(i,j,k)*(eez(k)*Dz(i,j,k)-ffz(k)*temp)

end subroutine pml_E_field_update

subroutine pcm_E_field_update(  i,   j,   k)
integer::i,j,k
	Ex(i,j,k) = 0.0 
	Ey(i,j,k) = 0.0 
	Ez(i,j,k) = 0.0 

end subroutine pcm_E_field_update

subroutine normal_H_field_update(  i,   j,   k)
integer::i,j,k
	Hx(i,j,k) = Hx(i,j,k) - tus*((Ez(i,j+1,k)-Ez(i,j,k))/ds_y - (Ey(i,j,k+1)-Ey(i,j,k))/ds_z) 
	Hy(i,j,k) = Hy(i,j,k) - tus*((Ex(i,j,k+1)-Ex(i,j,k))/ds_z - (Ez(i+1,j,k)-Ez(i,j,k))/ds_x) 
	Hz(i,j,k) = Hz(i,j,k) - tus*((Ey(i+1,j,k)-Ey(i,j,k))/ds_x - (Ex(i,j+1,k)-Ex(i,j,k))/ds_y) 

end subroutine normal_H_field_update

subroutine pml_H_field_update(  i,   j,   k)
integer::i,j,k
real*8:: temp 

	temp=Bx(i,j,k) 
	Bx(i,j,k) = ggy(j)*Bx(i,j,k) - hhy(j)*((Ez(i,j+1,k)-Ez(i,j,k))/ds_y - (Ey(i,j,k+1)-Ey(i,j,k))/ds_z) 
	Hx(i,j,k) = iiz(k)*Hx(i,j,k) + jjz(k)*(kkx(i)*Bx(i,j,k)-llx(i)*temp) 
	temp=By(i,j,k) 
	By(i,j,k) = ggz(k)*By(i,j,k) - hhz(k)*((Ex(i,j,k+1)-Ex(i,j,k))/ds_z - (Ez(i+1,j,k)-Ez(i,j,k))/ds_x) 
	Hy(i,j,k) = iix(i)*Hy(i,j,k) + jjx(i)*(kky(j)*By(i,j,k)-lly(j)*temp) 
	temp=Bz(i,j,k) 
	Bz(i,j,k) = ggx(i)*Bz(i,j,k) - hhx(i)*((Ey(i+1,j,k)-Ey(i,j,k))/ds_x - (Ex(i,j+1,k)-Ex(i,j,k))/ds_y) 
	Hz(i,j,k) = iiy(j)*Hz(i,j,k) + jjy(j)*(kkz(k)*Bz(i,j,k)-llz(k)*temp)


end subroutine pml_H_field_update

subroutine pcm_H_field_update(  i,   j,   k)
integer::i,j,k
	Hx(i,j,k) = 0.0 
	Hy(i,j,k) = 0.0 
	Hz(i,j,k) = 0.0 
end subroutine pcm_H_field_update

subroutine E_field_Gamma_boundary_update_x()
integer:: j,k 
	do k=1,pksize-1
	do j=1,pjsize-1
		Ex(0,j,k)=Ex(isize,j,k) 
		Ey(isize+1,j,k)=Ey(1,j,k) 
		Ez(isize+1,j,k)=Ez(1,j,k) 
	enddo
        enddo
end subroutine E_field_Gamma_boundary_update_x

subroutine E_field_Gamma_boundary_update_y()
integer:: i,k 
	do k=1,pksize-1
	do i=1,pisize-1
		Ey(i,0,k)=Ey(i,jsize,k) 
		Ex(i,jsize+1,k)=Ex(i,1,k) 
		Ez(i,jsize+1,k)=Ez(i,1,k) 
	enddo
        enddo
end subroutine E_field_Gamma_boundary_update_y

subroutine H_field_Gamma_boundary_update_x()
integer:: j,k 
	do k=1,pksize-1
	do j=1,pjsize-1
		Hx(isize+1,j,k)=Hx(1,j,k) 
		Hy(0,j,k)=Hy(isize,j,k) 
		Hz(0,j,k)=Hz(isize,j,k) 
	enddo
        enddo
end subroutine H_field_Gamma_boundary_update_x

subroutine H_field_Gamma_boundary_update_y()
integer:: i,k 
	do k=1,pksize-1
	do i=1,pisize-1
		Hy(i,jsize+1,k)=Hy(i,1,k) 
		Hx(i,0,k)=Hx(i,jsize,k) 
		Hz(i,0,k)=Hz(i,jsize,k) 
	enddo
        enddo
end subroutine H_field_Gamma_boundary_update_y


subroutine TF_SF_correction_for_H
integer:: i,j,k,il,ir,jl,jr,kl,kr

		il=TF_SF_pos(1,1)
		ir=TF_SF_pos(1,2)
		jl=TF_SF_pos(2,1)
		jr=TF_SF_pos(2,2)
		kl=TF_SF_pos(3,1)
		kr=TF_SF_pos(3,2)
		!/////// Bottom lid at k = kl		Hy(i+1/2,j,k+1/2)
		do j=jl,jr	        !(j0 ->j1)
			do i=il,ir-1    !(i0+1/2 -> i1-1/2)
				Hy(i,j,kl-1) = Hy(i,j,kl-1) + tus/ds_z*Ex_inc(kl)	!k0-1/2
			enddo
		enddo
		!/////// Top lid at k = kr
		do j=jl,jr		!(j0 ->j1)	
			do i=il,ir-1	!(i0+1/2 -> i1-1/2)
		   		Hy(i,j,kr) = Hy(i,j,kr) - tus/ds_z*Ex_inc(kr)		!k1+1/2
			enddo
		enddo

		!/////// y-walls in the left and in the right	Hz(i+1/2,j+1/2,k)
		do k=kl,kr		!(k0 ->k1)
			do i=il,ir-1    !(i0+1/2 -> i1-1/2)
				Hz(i,jl-1,k) = Hz(i,jl-1,k) - tus/ds_y*Ex_inc(k)	!j0-1/2
				Hz(i,jr,k) = Hz(i,jr,k) + tus/ds_y*Ex_inc(k)		!j1+1/2
      			!Hx~ Ez_inc
			enddo
		enddo
		!xwall Hy~Ez_inc Hz~Ey_inc
end subroutine TF_SF_correction_for_H

subroutine TF_SF_correction_for_E
integer:: i,j,k,il,ir,jl,jr,kl,kr

		il=TF_SF_pos(1,1)
		ir=TF_SF_pos(1,2)
		jl=TF_SF_pos(2,1)
		jr=TF_SF_pos(2,2)
		kl=TF_SF_pos(3,1)
		kr=TF_SF_pos(3,2)
		!/////// Bottom lid at k = kl
		do j=jl,jr		 !(j0 ->j1)
			do i=il,ir-1 	 !(i0+1/2 -> i1-1/2)
				Ex(i,j,kl) = Ex(i,j,kl) + (dt/eo/epsilonx(i,j,kl-1))/ds_z*Hy_inc(kl-1) !k0-1/2
				!Ey~Hx_inc no
			enddo
		enddo
		!/////// Top lid at k = kr  Ex(i+1/2,j,k)   
		do j=jl,jr		 !(j0 ->j1)
			do i=il,ir-1	 !(i0+1/2 -> i1-1/2)
				Ex(i,j,kr) = Ex(i,j,kr) - (dt/eo/epsilonx(i,j,kr))/ds_z*Hy_inc(kr)  !k1+1/2
				!Ey~Hx_inc no
			enddo
		enddo
		!/////// x-walls in the left and in the right Ez(i,j,k+1/2)
		do k=kl,kr-1	 	 !(k0+1/2 -> k1-1/2)
			do j=jl,jr	 !(j0 -> j1)
				Ez(il,j,k) = Ez(il,j,k) - (dt/eo/epsilonz(il-1,j,k))/ds_x*Hy_inc(k) !i0
				Ez(ir,j,k) = Ez(ir,j,k) + (dt/eo/epsilonz(ir,j,k))/ds_x*Hy_inc(k)   !i1
       			!Ey~Hz_inc no
			enddo
		enddo
		!ywall Ex~Hz_inc	 Ez~Hx_inc
end subroutine TF_SF_correction_for_E

subroutine H_inc_field_update_for_TF_SF
integer:: k, ks,nn
real*8:: temp_0, temp_1, temp_p1, temp_p0

	ks=light_pos
	temp_0 = Hy_inc(0)
	temp_1 = Hy_inc(1)
	temp_p1 = Hy_inc(pksize-1)
	temp_p0 = Hy_inc(pksize)

	do k=1,pksize-1
	       !Hy(i,j,k) = Hy(i,j,k) - tus*((Ex(i,j,k+1)-Ex(i,j,k))/ds_z - (Ez(i+1,j,k)-Ez(i,j,k))/ds_x) 
		Hy_inc(k) = Hy_inc(k) - tus*((Ex_inc(k+1)-Ex_inc(k))/ds_z)
	enddo
	
	Hy_inc(0)=temp_1+(light_speed*dt/ds_z-sqrt(TF_SF_eps_bg))/(light_speed*dt/ds_z+sqrt(TF_SF_eps_bg))*(Hy_inc(1) - temp_0)
	Hy_inc(pksize)=temp_p1 + (light_speed*dt/ds_z -sqrt(TF_SF_eps_bg))/(light_speed*dt/ds_z +sqrt(TF_SF_eps_bg))*(Hy_inc(pksize-1) - temp_p0)

!	//////////////////////////// TF_SF source amplitude //////////////////////////
	if(TF_SF_tdecay == 0 .and.  t>=TF_SF_to) then !// CW source
	  Hy_inc(ks) =amplitude*sqrt(TF_SF_eps_bg*eo/uo)*dsin( 2*pi*TF_SF_freq*(t-TF_SF_to)/(S_factor*ds_x*lattice_x) )		
	else if(TF_SF_to-3*TF_SF_tdecay<t .and. t<TF_SF_to+3*TF_SF_tdecay) then !// Gaussian pulse
     if(sourcetyp(1:8)=='Gaussian') then
	     if (.not. triple) then
	    	Hy_inc(ks) =amplitude*sqrt(TF_SF_eps_bg*eo/uo)*dsin( 2*pi*TF_SF_freq*(t-TF_SF_to)/(S_factor*ds_x*lattice_x) )*dexp(-1.d0*(dble(t-TF_SF_to)/dble(TF_SF_tdecay))**2.d0)
	  !	Hy_inc(ks) =amplitude*sqrt(TF_SF_eps_bg*eo/uo)*dsin(-pi+ 2*pi*TF_SF_freq*(t+dt/2.d0)/(S_factor*ds_x*lattice_x) )*dexp(-1.d0*(dble(t+dt/2.d0-TF_SF_to)/dble(TF_SF_tdecay))**2.d0)
	     else
	    	Hy_inc(ks)=0.d0
	        do nn=1,3
	            Hy_inc(ks)=Hy_inc(ks)+amplitude*sqrt(TF_SF_eps_bg*eo/uo)*dsin( 2*pi*TF_SF_freq*(t-TF_SF_to)/(S_factor*ds_x*lattice_x)*dble(nn*0.5))*dexp(-1.d0*(dble(t-TF_SF_to)/dble(TF_SF_tdecay))**2.d0)
	        enddo
	     endif
     elseif(sourcetyp(1:8)=='PGaussian') then
	    	Hy_inc(ks) =amplitude*sqrt(TF_SF_eps_bg*eo/uo)*dexp(-1.d0*(dble(t-TF_SF_to)/dble(TF_SF_tdecay))**2.d0)
     endif
	else 
		Hy_inc(ks) = Hy_inc(ks)
	endif
end subroutine H_inc_field_update_for_TF_SF

subroutine E_inc_field_update_for_TF_SF
integer:: k, ks,nn
real*8:: temp_0, temp_1, temp_p1, temp_p0

	ks=light_pos 

	temp_0 = Ex_inc(0)
	temp_1 = Ex_inc(1)
	temp_p1 = Ex_inc(pksize-1)
	temp_p0 = Ex_inc(pksize)

	do k=1,pksize-1
  		!Ex(i,j,k) = Ex(i,j,k) + (dt/eo/epsilonx(i,j,k))*((Hz(i,j,k)-Hz(i,j-1,k))/ds_y - (Hy(i,j,k)-Hy(i,j,k-1))/ds_z) 	
		Ex_inc(k) = Ex_inc(k) + (dt/eo/TF_SF_eps_bg)*(-((Hy_inc(k)-Hy_inc(k-1)))/ds_z)
	enddo
	!abc boundary at z (Ex Ey)
	Ex_inc(0)=temp_1+(light_speed*dt/ds_z -sqrt(TF_SF_eps_bg))/(light_speed*dt/ds_z +sqrt(TF_SF_eps_bg))*(Ex_inc(1) - temp_0)
	Ex_inc(pksize) = temp_p1 + (light_speed*dt/ds_z -sqrt(TF_SF_eps_bg))/(light_speed*dt/ds_z +sqrt(TF_SF_eps_bg))*(Ex_inc(pksize-1) - temp_p0)

	!---------------------------- TF_SF source amplitude ---------------------------could be neglected dsince already get by Hy_inc
	if(TF_SF_tdecay == 0 .and. t>=TF_SF_to) then !// CW source      !pi is the source of -Ex
		Ex_inc(ks) =amplitude*dsin(pi+2*pi*TF_SF_freq*(t-TF_SF_to)/(S_factor*ds_x*lattice_x) )		
	else if(TF_SF_to-3*TF_SF_tdecay<t .and. t<TF_SF_to+3*TF_SF_tdecay) then ! // Gaussian pulse
       if(sourcetyp(1:8)=='Gaussian') then
	       if (.not. triple) then
		Ex_inc(ks) =amplitude*dsin(pi+2*pi*TF_SF_freq*(t-TF_SF_to)/(S_factor*ds_x*lattice_x))* &
			     dexp(-1.d0*(dble(t-TF_SF_to)/dble(TF_SF_tdecay))**2.d0)
		!Ex_inc(ks) =amplitude*dsin(2*pi*TF_SF_freq*(t+dt/2.d0)/(S_factor*ds_x*lattice_x))* &
		!	     dexp(-1.d0*(dble(t+dt/2.d0-TF_SF_to)/dble(TF_SF_tdecay))**2.d0)
!		write(6,*) 'Ex_inc',dsin(pi+2*pi*TF_SF_freq*(t-TF_SF_to)/(S_factor*ds_x*lattice_x)), & 
!			    dexp(-1.d0*(dble(t-TF_SF_to)/dble(TF_SF_tdecay))**2.d0),Ex_inc(ks)
	       else
		 Ex_inc(ks)=0.d0
		 do nn=1,3
		 Ex_inc(ks)=Ex_inc(ks)+amplitude*dsin(pi+2*pi*TF_SF_freq*(t-TF_SF_to)/(S_factor*ds_x*lattice_x)&
				   *dble(nn*0.5))*dexp(-1.0*(dble(t-TF_SF_to)/dble(TF_SF_tdecay))**2.d0)
	         enddo
	       endif
        elseif(sourcetyp(1:8)=='PGaussian') then
		    Ex_inc(ks) =amplitude*dexp(-1.d0*(dble(t-TF_SF_to)/dble(TF_SF_tdecay))**2.d0)
        endif

	else 
		Ex_inc(ks) = Ex_inc(ks)
	endif

end subroutine E_inc_field_update_for_TF_SF



subroutine fft_field(i,j,k)
integer::i,j,k,m,ks
real*8::f
	
do m=1,Nf
	f=fmin+df*(m-1)
	Ex_Re(m,i,j,k)=Ex_Re(m,i,j,k)+dcos(2*pi*f*t/(S_factor*ds_x*lattice_x))*Ex(i,j,k)
	Ex_Im(m,i,j,k)=Ex_Im(m,i,j,k)+dsin(2*pi*f*t/(S_factor*ds_x*lattice_x))*Ex(i,j,k)

	Ey_Re(m,i,j,k)=Ey_Re(m,i,j,k)+dcos(2*pi*f*t/(S_factor*ds_x*lattice_x))*Ey(i,j,k)
	Ey_Im(m,i,j,k)=Ey_Im(m,i,j,k)+dsin(2*pi*f*t/(S_factor*ds_x*lattice_x))*Ey(i,j,k)

	Ez_Re(m,i,j,k)=Ez_Re(m,i,j,k)+dcos(2*pi*f*t/(S_factor*ds_x*lattice_x))*Ez(i,j,k)
	Ez_Im(m,i,j,k)=Ez_Im(m,i,j,k)+dsin(2*pi*f*t/(S_factor*ds_x*lattice_x))*Ez(i,j,k)

	Hx_Re(m,i,j,k)=Hx_Re(m,i,j,k)+dcos(2*pi*f*t/(S_factor*ds_x*lattice_x))*Hx(i,j,k)
	Hx_Im(m,i,j,k)=Hx_Im(m,i,j,k)+dsin(2*pi*f*t/(S_factor*ds_x*lattice_x))*Hx(i,j,k)

	Hy_Re(m,i,j,k)=Hy_Re(m,i,j,k)+dcos(2*pi*f*t/(S_factor*ds_x*lattice_x))*Hy(i,j,k)
	Hy_Im(m,i,j,k)=Hy_Im(m,i,j,k)+dsin(2*pi*f*t/(S_factor*ds_x*lattice_x))*Hy(i,j,k)

	Hz_Re(m,i,j,k)=Hz_Re(m,i,j,k)+dcos(2*pi*f*t/(S_factor*ds_x*lattice_x))*Hz(i,j,k)
	Hz_Im(m,i,j,k)=Hz_Im(m,i,j,k)+dsin(2*pi*f*t/(S_factor*ds_x*lattice_x))*Hz(i,j,k)
enddo
end subroutine fft_field
	
subroutine fft_field_inc
integer::m,ks
real*8::f
	
ks=light_pos 
do m=1,Nf
	f=fmin+df*(m-1)
	Ex_inc_Re(m)=Ex_inc_Re(m)+dcos(2*pi*f*t/(S_factor*ds_x*lattice_x))*Ex_inc(ks)
	Ex_inc_Im(m)=Ex_inc_Im(m)+dsin(2*pi*f*t/(S_factor*ds_x*lattice_x))*Ex_inc(ks)
enddo
end subroutine fft_field_inc
!*!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!///
real*8 function sigmax(a)
real*8::a
real*8:: sigma_maxl, sigma_maxr
!Taflov 331
!sigma_max_opt =sig_axl*(m + 1)))/(Z0*grid_dx*eps_eff)
!sig_axl=0.8~1.4(0.8 is best)  m=3~5  , eps_eff=sqrt(eps_r*mu_r)  eps_r =1.0,mu_r=1.0, Z0=sqrt(MU0/EPS0)=MU0*C=pi*4e-7*3e8=120*pi
!here grid_dx = 1, since dt/grid_dx =1/C/S    
! \sigma = (x/d)^m*sigma_max
	 !// no use of periodic boundary condition 
	        	
	        sigma_maxl=(sig_axl)*(orderxl+1)/(120*pi*sqrt(back_epsilon))
		sigma_maxr=(sig_axr)*(orderxr+1)/(120*pi*sqrt(back_epsilon))

		if(a<pmlil+0.5) then
		      sigmax=sigma_maxl*(-(a-(pmlil+0.5))/pmlil)**orderxl
		!here is tricky
		elseif(a<isize-pmlir-1) then 
		      sigmax=0.0
		else
		      sigmax=sigma_maxr*((a-(isize-pmlir-1))/pmlir)**orderxr
                endif


end function sigmax

real*8 function sigmay(a)
real*8::a
real*8:: sigma_maxl, sigma_maxr

	 ! // no use of periodic boundary condition 
		sigma_maxl=(sig_ayl)*(orderyl+1)/(120*pi*sqrt(back_epsilon))
		sigma_maxr=(sig_ayr)*(orderyr+1)/(120*pi*sqrt(back_epsilon))

		if(a<pmljl+0.5) then 
                       sigmay=sigma_maxl*(-(a-(pmljl+0.5))/pmljl)**orderyl
		elseif(a<jsize-pmljr-1) then
		       sigmay=0.0
		else
		      sigmay=sigma_maxr*((a-(jsize-pmljr-1))/pmljr)**orderyr
	        endif
end function sigmay

real*8 function sigmaz(a)
real*8::a
real*8:: sigma_maxl, sigma_maxr

	sigma_maxl=(sig_azl)*(orderzl+1)/(120*pi*sqrt(back_epsilon))
	sigma_maxr=(sig_azr)*(orderzr+1)/(120*pi*sqrt(back_epsilon))

	if(a<pmlkl+0.5) then
                    sigmaz=sigma_maxl*(-(a-(pmlkl+0.5))/pmlkl)**orderzl
	elseif(a<ksize-pmlkr-1) then
	            sigmaz=0.0
	else
                    sigmaz=sigma_maxr*((a-(ksize-pmlkr-1))/pmlkr)**orderzr

       endif
end function sigmaz

real*8 function kappax(a)
!i. optimum values of kappamax are between 7 and 20
!ii. optimum values of sigmamax are between 0.8*sigma_opt and 1.4*sigma_opt
!iii. optimum values of alphamax are 0.15 to 0.3
!\kappa = 1.0+(\kappa_max-1)*(x/d)^m
!
!3<<m<<4
real*8::a
real*8::kappa_maxl,kappa_maxr
	
	kappa_maxl=14
	kappa_maxr=14
	if(a<pmlil+0.5) then
		kappax=1.0+(kappa_maxl-1)*(-(a-(pmlil+0.5))/pmlil)**orderxl
	elseif(a<isize-pmlir-1) then 
		kappax=1.0
	else
		kappax=1.0+(kappa_maxr-1)*((a-(isize-pmlir-1))/pmlir)**orderxr
        endif
end function kappax	

real*8 function kappay(a)
real*8::a
real*8::kappa_maxl,kappa_maxr
	
	kappa_maxl=14
	kappa_maxr=14
	if(a<pmljl+0.5) then
		kappay=1.0+(kappa_maxl-1)*(-(a-(pmljl+0.5))/pmljl)**orderyl
	elseif(a<jsize-pmljr-1) then 
		kappay=1.0
	else
		kappay=1.0+(kappa_maxr-1)*((a-(jsize-pmljr-1))/pmljr)**orderyr
        endif
end function kappay	

real*8 function kappaz(a)
real*8::a
real*8::kappa_maxl,kappa_maxr
	
	kappa_maxl=14
	kappa_maxr=14
	if(a<pmlkl+0.5) then
		kappaz=1.0+(kappa_maxl-1)*(-(a-(pmlkl+0.5))/pmlkl)**orderzl
	elseif(a<ksize-pmlkr-1) then 
		kappaz=1.0
	else
		kappaz=1.0+(kappa_maxr-1)*((a-(ksize-pmlkr-1))/pmlkr)**orderzr
        endif
end function kappaz
end module timeupdate
