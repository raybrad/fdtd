module inputobj
use parameters
implicit none

integer:: Lobjn=0,Mobjn=0  !// number of Lorentz objects 
type(Lobj),save:: Lobject(0:8)  !// declaration of 1D-array, 'Lobject',for now just input maxmum two obj,if need to add more ,just change
type(Mobj),save::Mobject(0:0)
!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//
contains
subroutine background(epsilon0)
      integer:: i,j,k
      real*8::epsilon0
      back_epsilon=epsilon0
      

      do k=1,ksize
      do j=1,jsize
      do i=1,isize
             epsilonx(i,j,k)=back_epsilon
             epsilony(i,j,k)=back_epsilon
             epsilonz(i,j,k)=back_epsilon
      enddo
      enddo
      enddo
      print *,"background...ok"
end subroutine background

subroutine monitor(centerx,centery,centerz,size1,size2,size3)
real*8,intent(in)::centerx,centery,centerz,size1,size2,size3
integer::i,j,k, x_n,x_p, y_n,y_p, z_n,z_p  !  //n:- p:+, position of the sides of the block

	x_n=(centerx + xcenter - 0.5*size1)*lattice_x
	x_p=(centerx + xcenter + 0.5*size1)*lattice_x
	y_n=(centery + ycenter - 0.5*size2)*lattice_y
	y_p=(centery + ycenter + 0.5*size2)*lattice_y
	z_n=(centerz + zcenter - 0.5*size3)*lattice_z
	z_p=(centerz + zcenter + 0.5*size3)*lattice_z

	do k=z_n,z_p,1
		do j=y_n,y_p,1 
		lmonitor(x_p,j,k)=.true.
		lmonitor(x_n,j,k)=.true.
		enddo
	enddo	

	do i=x_n,x_p,1
		do k=z_n,z_p,1
		lmonitor(i,y_p,k)=.true.
		lmonitor(i,y_n,k)=.true.
		enddo
	enddo
	do i=x_n,x_p,1
		do j=y_n,y_p,1 
		lmonitor(i,j,z_n)=.true.
		lmonitor(i,j,z_p)=.true.
		enddo
	enddo

end subroutine monitor

subroutine molecular_struct(cshape,centerx,centery,centerz,size1,size2,size3,nsites)
! define region of molecular
real*8::centerx,centery,centerz,size1,size2,size3
integer::nsites
character(len=*)::cshape

Mobjn=Mobjn+1

print *,'Mobjn',Mobjn

Mobject(Mobjn-1).cshape=cshape

if (cshape=='point') then
Mobject(Mobjn-1).centeri=(centerx+xcenter)*lattice_x
Mobject(Mobjn-1).centerj=(centery+ycenter)*lattice_y
Mobject(Mobjn-1).centerk=(centerz+zcenter)*lattice_z
elseif (cshape=='block') then
Mobject(Mobjn-1).centeri=(centerx+xcenter)*lattice_x
Mobject(Mobjn-1).centerj=(centery+ycenter)*lattice_y
Mobject(Mobjn-1).centerk=(centerz+zcenter)*lattice_z
Mobject(Mobjn-1).size1=size1*lattice_x   
Mobject(Mobjn-1).size2=size2*lattice_y
Mobject(Mobjn-1).size3=size3*lattice_z



endif

print *,'make molecular structure complete',cshape
end subroutine molecular_struct
	

subroutine input_metal_para(cshape,centerx,centery,centerz,centerxl,size1,size2,size3,epsilon_b,omega_p,gamma_0,sigma_0,lattice_n)

character(len=*)::cshape              ! matrix_file(not used now)
real*8::centerx,centery,centerz,centerxl,size1,size2,size3,epsilon_b,omega_p,gamma_0,sigma_0,lattice_n	
real*8:: omega_pn, gamma_0n,sigma_0n  !//reduced frequency, NOT normalized frequency
	!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//
	!/// How to convert into the reduced unit .
	!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//
	!/*
        !     phase changeend=w_p * dt_R = w_F * dt_F 
        !     Here, w_p = angular frequency in unit of Hz
        !           dt_R = time step in unit of sec.
        !           w_F = reduced frequency in FDTD
        !           dt_F = FDTD time step 1 
        !     Remember that dt_F is NOT assumed to be 1 
        !            BUT represented in unit of (sec q m^-1) (See Eq.49 of manual)
             
        !     Now that, 
        !           dt_F = 1/(c * S(q^-1))
        !           dt_R = 1/(c * S(m^-1))
        !                = 1/(c * S(q^-1) * ds_x * lattice_x / lattice_n ) 
        !     Therefore we get,
        !           dt_R/dt_F = lattice_n/(ds_x*lattice_x)              
        !           w_F = w_p * lattice_n/(ds_x*lattice_x) 
	!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!/// */
!omega~ s^-1
	omega_pn = omega_p*lattice_n/(ds_x*lattice_x) 
	gamma_0n = gamma_0*lattice_n/(ds_x*lattice_x) 		
	sigma_0n = sigma_0*lattice_n/(ds_x*lattice_x)
	if(Lobjn==0) then !//print once!
	 
		print *,"----------------",omega_pn 	
		print *,"omega_pn = ",omega_pn 
		print *,"gamma_0n = ",gamma_0n 
		print *,"sigma_0n = ",sigma_0n 
	endif

	Lobjn=Lobjn+1 
        print *,'Lobjn',Lobjn


	Lobject(Lobjn-1).cshape=cshape

	!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//
	!//!//!//!// No Euler rotation !//!//!///
	!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//

	Lobject(Lobjn-1).epsilon_b=epsilon_b   
	Lobject(Lobjn-1).omega_p=omega_pn   !//conversion to 'FDTD' frequency
	Lobject(Lobjn-1).gamma_0=gamma_0n      !//         "      "
	Lobject(Lobjn-1).sigma_0=sigma_0n      !//         "      "

	
	if(cshape=="sphere") then!//!// not applicable for non_uniform_grid()
	 
		Lobject(Lobjn-1).centeri=(centerx+xcenter)*lattice_x 
		Lobject(Lobjn-1).centerj=(centery+ycenter)*lattice_y 
		Lobject(Lobjn-1).centerk=(centerz+zcenter)*lattice_z 
		Lobject(Lobjn-1).size1=size1*lattice_x   !// radius
	
	elseif(cshape=="block") then
	 
		Lobject(Lobjn-1).centeri=(centerx+xcenter)*lattice_x 
		Lobject(Lobjn-1).centerj=(centery+ycenter)*lattice_y 
		Lobject(Lobjn-1).centerk=(centerz+zcenter)*lattice_z   
		Lobject(Lobjn-1).size1=size1*lattice_x 
		Lobject(Lobjn-1).size2=size2*lattice_y 
		Lobject(Lobjn-1).size3=size3*lattice_z  
	elseif(cshape=='pyramid') then
		Lobject(Lobjn-1).centeri=(centerx+xcenter)*lattice_x 
		Lobject(Lobjn-1).centerj=(centery+ycenter)*lattice_y 
		Lobject(Lobjn-1).centerk=(centerz+zcenter)*lattice_z  
		Lobject(Lobjn-1).centerl=(centerxl+xcenter)*lattice_x 
		Lobject(Lobjn-1).size1=size1*lattice_x 
		Lobject(Lobjn-1).size2=size2*lattice_y 
		Lobject(Lobjn-1).size3=size3*lattice_z  
	elseif(cshape=="rod") then
	 
		Lobject(Lobjn-1).centeri=(centerx+xcenter)*lattice_x 
		Lobject(Lobjn-1).centerj=(centery+ycenter)*lattice_y 
		Lobject(Lobjn-1).centerk=(centerz+zcenter)*lattice_z 
		Lobject(Lobjn-1).size1=size1*lattice_x   
		Lobject(Lobjn-1).size2=size2*lattice_z   
	
	
	else
	      continue
	endif
end subroutine input_metal_para



subroutine assign_metal_grid()
 
integer:: i,j,k 
integer:: n 

real*8:: Le, Lw, Lr, Lo ,Ls,Mole 

	print *,"-------------------------------" 	
	print *,"making the Lorentz structure..." 

	!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!///
	!///        Medium Rules for (i,j,k)            		        !///
	!///---------------------------------           			!///
	!///  Lepsilon=Lomega= 0.0 --> dielectric      				!///
	!///  Lepsilon!=0.0 & Lepsilon/=1000 --> Drude metal			!///
	!///  Lepsilon=1000 --> PCM                      			!///
	!///  Lepsilon=0.0 & Lomega=0.0 --> medium eraser 			!///
	!///  Lsigma/=0.0		--> normal metal			!///
	!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!///
			 
	do k=1,ksize
		print *,100*k/ksize,'%'  
		do j=1,jsize
			do i=1,isize
				Le = 0.0;Lw =0.0; Lr =0.0; Lo = 0.0  !// initialization : dielectric
				Mole=0.0
				do n=0,Lobjn-1
					if(in_Lobject(n,i,j,k)==1) then
						Le=Lobject(n).epsilon_b 
						Lw=Lobject(n).omega_p 
						Lr=Lobject(n).gamma_0 
						Ls=Lobject(n).sigma_0 
						!newly added,metal could extend to some part of pml
		                		cposition(i,j,k)=1 !non-pml
					endif
				enddo

				do n=0,Mobjn-1
				     if (in_Mobject(n,i,j,k)==1) then
				         Mole=1.0
					 print *,'mole',i,j,k
			            endif
				enddo     
                               !  
                                Molecular(i,j,k)=Mole 
				Lepsilon(i,j,k) = Le    !if Lsigma /=0.0 ,Lepsilon should be set to 1.0 for metal
				Lomega(i,j,k) = Lw  
				Lgamma(i,j,k) = Lr 
				Lsigma(i,j,k) = Ls  
			enddo
		enddo	
	enddo                                                                
                                     
   !     deallocate(Lobject)                  

	print *,"assign__metal_grid...ok" 
end subroutine assign_metal_grid


function in_Mobject(n,i,j,k)
integer::n
integer::i,j,k
integer::in_Mobject
if (Mobject(n).cshape=='point') then
	if (i==floor(0.5+Mobject(n).centeri) .and. j==floor(0.5+Mobject(n).centerj) .and. k==floor(0.5+Mobject(n).centerk)) then
      		in_Mobject=1
	else
      		in_Mobject=0
	endif

elseif(Mobject(n).cshape=='block') then
 	if( Mobject(n).centeri-Mobject(n).size1/2<=i .and. i<=Mobject(n).centeri+Mobject(n).size1/2 .and. &
		Mobject(n).centerj-Mobject(n).size2/2<=j .and. j<=Mobject(n).centerj+Mobject(n).size2/2 .and.     &
	Mobject(n).centerk-0.5*Mobject(n).size3<=k .and. k<=Mobject(n).centerk+0.5*Mobject(n).size3)  then 
			in_Mobject= 1 
		else 
			in_Mobject= 0 
		endif	

endif
end function in_Mobject

function in_Lobject(n,i,j,k)
integer::n
integer::i,j,k
integer::in_Lobject
integer:: a, b  
real*8:: X, Y, Z 
integer:: temp ,temp1


	if(Lobject(n).cshape=="sphere") then
		if( (i-Lobject(n).centeri)*(i-Lobject(n).centeri)+(j-Lobject(n).centerj)*(j-Lobject(n).centerj)+(k-Lobject(n).centerk)*(k-Lobject(n).centerk)<=Lobject(n).size1*Lobject(n).size1 ) then
		     in_Lobject= 1 
		else 
		     in_Lobject= 0 
		endif

	elseif(Lobject(n).cshape=="block") then
	 
		if( Lobject(n).centeri-Lobject(n).size1/2<=i .and. i<=Lobject(n).centeri+Lobject(n).size1/2 .and. &
		Lobject(n).centerj-Lobject(n).size2/2<=j .and. j<=Lobject(n).centerj+Lobject(n).size2/2 .and.     &
		Lobject(n).centerk-Lobject(n).size3/2<=k .and. k<=Lobject(n).centerk+Lobject(n).size3/2 ) then 
			in_Lobject= 1 
		else 
			in_Lobject= 0 
		endif	

	elseif(Lobject(n).cshape=="pyramid") then
	  	if (Lobject(n).centeri-Lobject(n).size1/2<=i .and. i<=Lobject(n).centeri+Lobject(n).size1/2) then  
			if( Lobject(n).centerj-Lobject(n).size2/2<=j .and. j<=Lobject(n).centerj+Lobject(n).size2/2 .and.     &
			Lobject(n).centerk-Lobject(n).size3/2<=k .and. k<=Lobject(n).centerk+Lobject(n).size3/2 ) then 
			in_Lobject= 1 
			else 
			in_Lobject= 0 
			endif
	  	elseif (i<= Lobject(n).centerl .and. i>=Lobject(n).centeri+Lobject(n).size1/2) then
	  	temp=(Lobject(n).centerl-i)*Lobject(n).size2/(Lobject(n).centerl-Lobject(n).centeri-Lobject(n).size1/2)    !tempsize2(y)
	  	temp1=(Lobject(n).centerl-i)*Lobject(n).size3/(Lobject(n).centerl-Lobject(n).centeri-Lobject(n).size1/2) 
	    !    print *,'x',temp,temp1
	        	if(Lobject(n).centerj-temp/2<=j .and.  Lobject(n).centerj+temp/2>=j .and. Lobject(n).centerj-temp1/2<=k .and. Lobject(n).centerj+temp1/2>=k ) then
		    	in_Lobject= 1
			else
			in_Lobject= 0
			endif
	  	elseif (i>= Lobject(n).centerl .and. i<Lobject(n).centeri-Lobject(n).size1/2) then
	  	temp=(Lobject(n).centerl-i)*Lobject(n).size2/(Lobject(n).centerl-Lobject(n).centeri+Lobject(n).size1/2)    !tempsize2(y)
	  	temp1=(Lobject(n).centerl-i)*Lobject(n).size3/(Lobject(n).centerl-Lobject(n).centeri+Lobject(n).size1/2)  
	        	if(Lobject(n).centerj-temp/2<=j .and.  Lobject(n).centerj+temp/2>=j .and. Lobject(n).centerj-temp1/2<=k .and. Lobject(n).centerj+temp1/2>=k ) then
		    	in_Lobject= 1
			else
			in_Lobject= 0
			endif
	          else
	 		in_Lobject=0

		  endif
	elseif(Lobject(n).cshape=="rod") then
	 
		if( (i-Lobject(n).centeri)*(i-Lobject(n).centeri)+(j-Lobject(n).centerj)*(j-Lobject(n).centerj)<=Lobject(n).size1*Lobject(n).size1 .and. (Lobject(n).centerk-0.5*Lobject(n).size2)<=k .and. k<=(Lobject(n).centerk+0.5*Lobject(n).size2) ) then 
		     in_Lobject= 1 
		else
		     in_Lobject= 0 
		endif


	else
	      continue
	endif
end function in_Lobject


end module inputobj 
