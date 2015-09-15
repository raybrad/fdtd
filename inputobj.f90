module inputobj
use parameters
implicit none

integer:: MetalObjn=0,MolecluarObjn=0  !// number of Lorentz objects 
type(MetalObj),save:: MetalObject(0:8)  !// declaration of 1D-array, 'MetalObject',for now just input maxmum two obj,if need to add more ,just change
type(MolecularObj),save::MolecularObject(0:0)
!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//
contains
subroutine background(epsilon0)
      integer:: i,j,k
      real*8::epsilon0
      back_epsilon=epsilon0
      

      do k=0,mksize-1
      do j=0,mjsize-1
      do i=0,misize-1
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

MolecluarObjn=MolecluarObjn+1

print *,'MolecluarObjn',MolecluarObjn

MolecularObject(MolecluarObjn-1).cshape=cshape

if (cshape=='point') then
MolecularObject(MolecluarObjn-1).centeri=(centerx+xcenter)*lattice_x
MolecularObject(MolecluarObjn-1).centerj=(centery+ycenter)*lattice_y
MolecularObject(MolecluarObjn-1).centerk=(centerz+zcenter)*lattice_z
elseif (cshape=='block') then
MolecularObject(MolecluarObjn-1).centeri=(centerx+xcenter)*lattice_x
MolecularObject(MolecluarObjn-1).centerj=(centery+ycenter)*lattice_y
MolecularObject(MolecluarObjn-1).centerk=(centerz+zcenter)*lattice_z
MolecularObject(MolecluarObjn-1).size1=size1*lattice_x   
MolecularObject(MolecluarObjn-1).size2=size2*lattice_y
MolecularObject(MolecluarObjn-1).size3=size3*lattice_z



endif

print *,'make molecular structure complete',cshape
end subroutine molecular_struct
	

subroutine input_metal_para(cshape,centerx,centery,centerz,centerxl,size1,size2,size3,epsilon_b,omega_p,gamma_0,sigma_0,depsr_0_1,omega_0_1,gamma_0_1,depsr_0_2,omega_0_2,gamma_0_2,lattice_n)
character(len=*)::cshape              ! matrix_file(not used now)
real*8::centerx,centery,centerz,centerxl,size1,size2,size3,epsilon_b,omega_p,gamma_0,sigma_0,lattice_n	
real*8::depsr_0_1,omega_0_1,gamma_0_1,depsr_0_2,omega_0_2,gamma_0_2
real*8::omega_pn, gamma_0n,sigma_0n  !//reduced frequency, NOT normalized frequency
real*8::omega_0n_1,gamma_0n_1,omega_0n_2,gamma_0n_2  
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
!omega~ eV -> rad/s^-1  ->  ..
	omega_pn = omega_p*2.d0*PI/PLANCKEV*lattice_n/(ds_x*lattice_x) 
	gamma_0n = gamma_0*2.d0*PI/PLANCKEV*lattice_n/(ds_x*lattice_x) 		
	sigma_0n = sigma_0*lattice_n/(ds_x*lattice_x)		!??
	
	omega_0n_1 = omega_0_1*2.d0*PI/PLANCKEV*lattice_n/(ds_x*lattice_x) 
	gamma_0n_1 = gamma_0_1*2.d0*PI/PLANCKEV*lattice_n/(ds_x*lattice_x) 	
	omega_0n_2 = omega_0_2*2.d0*PI/PLANCKEV*lattice_n/(ds_x*lattice_x) 
	gamma_0n_2 = gamma_0_2*2.d0*PI/PLANCKEV*lattice_n/(ds_x*lattice_x) 
	if(MetalObjn==0) then !//print once!
	 
		print *,"----------------" 	
		print *,"epsilon_b",epsilon_b
		print *,"omega_pn = ",omega_pn 
		print *,"gamma_0n = ",gamma_0n 
		print *,"sigma_0n = ",sigma_0n 
		print *,"depsr_0_1 = ",depsr_0_1 
		print *,"omega_0n_1 = ",omega_0n_1 
		print *,"gamma_0n_1 = ",gamma_0n_1 
		print *,"depsr_0_2 = ",depsr_0_2 
		print *,"omega_0n_2 = ",omega_0n_2 
		print *,"gamma_0n_2 = ",gamma_0n_2 
	endif

	MetalObjn=MetalObjn+1 
        print *,'MetalObjn',MetalObjn


	MetalObject(MetalObjn-1).cshape=cshape

	!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//
	!//!//!//!// No Euler rotation !//!//!///
	!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//

	MetalObject(MetalObjn-1).epsilon_b=epsilon_b   
	MetalObject(MetalObjn-1).omega_p=omega_pn   !//conversion to 'FDTD' frequency
	MetalObject(MetalObjn-1).gamma_0=gamma_0n      !//         "      "
	MetalObject(MetalObjn-1).sigma_0=sigma_0n      !//         "      "
	MetalObject(MetalObjn-1).depsr_0_1=depsr_0_1   !//conversion to 'FDTD' frequency
	MetalObject(MetalObjn-1).omega_0_1=omega_0n_1 
	MetalObject(MetalObjn-1).gamma_0_1=gamma_0n_1   
	MetalObject(MetalObjn-1).depsr_0_2=depsr_0_2   !//conversion to 'FDTD' frequency
	MetalObject(MetalObjn-1).omega_0_2=omega_0n_2 
	MetalObject(MetalObjn-1).gamma_0_2=gamma_0n_2   

	
	if(cshape=="sphere") then!//!// not applicable for non_uniform_grid()
	 
		MetalObject(MetalObjn-1).centeri=(centerx+xcenter)*lattice_x 
		MetalObject(MetalObjn-1).centerj=(centery+ycenter)*lattice_y 
		MetalObject(MetalObjn-1).centerk=(centerz+zcenter)*lattice_z 
		MetalObject(MetalObjn-1).size1=size1*lattice_x   !// radius
	
	elseif(cshape=="block") then
	 
		MetalObject(MetalObjn-1).centeri=(centerx+xcenter)*lattice_x 
		MetalObject(MetalObjn-1).centerj=(centery+ycenter)*lattice_y 
		MetalObject(MetalObjn-1).centerk=(centerz+zcenter)*lattice_z   
		MetalObject(MetalObjn-1).size1=size1*lattice_x 
		MetalObject(MetalObjn-1).size2=size2*lattice_y 
		MetalObject(MetalObjn-1).size3=size3*lattice_z  
	elseif(cshape=='pyramid') then
		MetalObject(MetalObjn-1).centeri=(centerx+xcenter)*lattice_x 
		MetalObject(MetalObjn-1).centerj=(centery+ycenter)*lattice_y 
		MetalObject(MetalObjn-1).centerk=(centerz+zcenter)*lattice_z  
		MetalObject(MetalObjn-1).centerl=(centerxl+xcenter)*lattice_x 
		MetalObject(MetalObjn-1).size1=size1*lattice_x 
		MetalObject(MetalObjn-1).size2=size2*lattice_y 
		MetalObject(MetalObjn-1).size3=size3*lattice_z  
	elseif(cshape=="rod") then
	 
		MetalObject(MetalObjn-1).centeri=(centerx+xcenter)*lattice_x 
		MetalObject(MetalObjn-1).centerj=(centery+ycenter)*lattice_y 
		MetalObject(MetalObjn-1).centerk=(centerz+zcenter)*lattice_z 
		MetalObject(MetalObjn-1).size1=size1*lattice_x   
		MetalObject(MetalObjn-1).size2=size2*lattice_z   
	
	
	else
	      continue
	endif
end subroutine input_metal_para



subroutine assign_metal_grid()
 
integer:: i,j,k 
integer:: n 

real*8:: Le, Lw, Lr, Lo ,Ls,Mole 
real*8:: Lepsr_1, Lr_1, Lo_1, Lepsr_2, Lr_2, Lo_2 

	print *,"-------------------------------" 	
	print *,"making structure..." 

	!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!///
	!///        Medium Rules for (i,j,k)            		        !///
	!///---------------------------------           			!///
	!///  Lepsilon=Lomega= 0.0 --> dielectric      				!///
	!///  Lepsilon!=0.0 & Lepsilon/=1000 --> Drude metal			!///
	!///  Lepsilon=1000 --> PCM                      			!///
	!///  Lepsilon=0.0 & Lomega=0.0 --> medium eraser 			!///
	!///  Lsigma/=0.0		--> normal metal			!///
	!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!//!///
			 
	do k=0,mksize-1
		print *,100*k/mksize,'%'  
		do j=0,mjsize-1
			do i=0,misize-1
				Le = 0.0;Lw =0.0; Lr =0.0; Lo = 0.0  !// initialization : dielectric
	                        Lepsr_1 =0.0; Lr_1 =0.0; Lo_1 =0.0
	                        Lepsr_2 =0.0; Lr_2 =0.0; Lo_2 =0.0  !// initialization : dielectric
				Mole=0.0
				do n=0,MetalObjn-1
					if(in_MetalObject(n,i,j,k)==1) then
						Le=MetalObject(n).epsilon_b 
						Lw=MetalObject(n).omega_p 
						Lr=MetalObject(n).gamma_0 
						Ls=MetalObject(n).sigma_0 
					if(metal_type==3 .or. metal_type==4 .or. metal_type==5) then	
						Lepsr_1=MetalObject(n).depsr_0_1 
						Lr_1=MetalObject(n).gamma_0_1 
						Lo_1=MetalObject(n).omega_0_1 
						Lepsr_2=MetalObject(n).depsr_0_2 
						Lr_2=MetalObject(n).gamma_0_2 
						Lo_2=MetalObject(n).omega_0_2 
					endif
						!newly added,metal could extend to some part of pml
		                		cposition(i,j,k)=1 !non-pml
					endif
			
				enddo

				do n=0,MolecluarObjn-1
				     if (in_MolecularObject(n,i,j,k)==1) then
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
				
				if(metal_type==3 .or. metal_type==4 .or. metal_type==5) then	
				Ldepsr_1(i,j,k) = Lepsr_1  
				Lgamma_1(i,j,k) = Lr_1 
				Lomega_1(i,j,k) = Lo_1  
				Ldepsr_2(i,j,k) = Lepsr_2  
				Lgamma_2(i,j,k) = Lr_2 
				Lomega_2(i,j,k) = Lo_2  
				endif
			enddo
		enddo	
	enddo                                                                
                                     
   !     deallocate(MetalObject)                  

	print *,"assign__metal_grid...ok" 
end subroutine assign_metal_grid

function in_MolecularObject(n,i,j,k)
integer::n
integer::i,j,k
integer::in_MolecularObject
if (MolecularObject(n).cshape=='point') then
	if (i==floor(0.5+MolecularObject(n).centeri) .and. j==floor(0.5+MolecularObject(n).centerj) .and. k==floor(0.5+MolecularObject(n).centerk)) then
      		in_MolecularObject=1
	else
      		in_MolecularObject=0
	endif

elseif(MolecularObject(n).cshape=='block') then
 	if( MolecularObject(n).centeri-MolecularObject(n).size1/2<=i .and. i<=MolecularObject(n).centeri+MolecularObject(n).size1/2 .and. &
		MolecularObject(n).centerj-MolecularObject(n).size2/2<=j .and. j<=MolecularObject(n).centerj+MolecularObject(n).size2/2 .and.     &
	MolecularObject(n).centerk-0.5*MolecularObject(n).size3<=k .and. k<=MolecularObject(n).centerk+0.5*MolecularObject(n).size3)  then 
			in_MolecularObject= 1 
		else 
			in_MolecularObject= 0 
		endif	

endif
end function in_MolecularObject

function in_MetalObject(n,i,j,k)
integer::n
integer::i,j,k
integer::in_MetalObject
integer:: a, b  
real*8:: X, Y, Z 
integer:: temp ,temp1


	if(MetalObject(n).cshape=="sphere") then
		if( (i-MetalObject(n).centeri)*(i-MetalObject(n).centeri)+(j-MetalObject(n).centerj)*(j-MetalObject(n).centerj)+(k-MetalObject(n).centerk)*(k-MetalObject(n).centerk)<=MetalObject(n).size1*MetalObject(n).size1 ) then
		     in_MetalObject= 1 
		else 
		     in_MetalObject= 0 
		endif

	elseif(MetalObject(n).cshape=="block") then
	 
		if( MetalObject(n).centeri-MetalObject(n).size1/2<=i .and. i<=MetalObject(n).centeri+MetalObject(n).size1/2 .and. &
		MetalObject(n).centerj-MetalObject(n).size2/2<=j .and. j<=MetalObject(n).centerj+MetalObject(n).size2/2 .and.     &
		MetalObject(n).centerk-MetalObject(n).size3/2<=k .and. k<=MetalObject(n).centerk+MetalObject(n).size3/2 ) then 
			in_MetalObject= 1 
		else 
			in_MetalObject= 0 
		endif	

	elseif(MetalObject(n).cshape=="pyramid") then
	  	if (MetalObject(n).centeri-MetalObject(n).size1/2<=i .and. i<=MetalObject(n).centeri+MetalObject(n).size1/2) then  
			if( MetalObject(n).centerj-MetalObject(n).size2/2<=j .and. j<=MetalObject(n).centerj+MetalObject(n).size2/2 .and.     &
			MetalObject(n).centerk-MetalObject(n).size3/2<=k .and. k<=MetalObject(n).centerk+MetalObject(n).size3/2 ) then 
			in_MetalObject= 1 
			else 
			in_MetalObject= 0 
			endif
	  	elseif (i<= MetalObject(n).centerl .and. i>=MetalObject(n).centeri+MetalObject(n).size1/2) then
	  	temp=(MetalObject(n).centerl-i)*MetalObject(n).size2/(MetalObject(n).centerl-MetalObject(n).centeri-MetalObject(n).size1/2)    !tempsize2(y)
	  	temp1=(MetalObject(n).centerl-i)*MetalObject(n).size3/(MetalObject(n).centerl-MetalObject(n).centeri-MetalObject(n).size1/2) 
	    !    print *,'x',temp,temp1
	        	if(MetalObject(n).centerj-temp/2<=j .and.  MetalObject(n).centerj+temp/2>=j .and. MetalObject(n).centerj-temp1/2<=k .and. MetalObject(n).centerj+temp1/2>=k ) then
		    	in_MetalObject= 1
			else
			in_MetalObject= 0
			endif
	  	elseif (i>= MetalObject(n).centerl .and. i<MetalObject(n).centeri-MetalObject(n).size1/2) then
	  	temp=(MetalObject(n).centerl-i)*MetalObject(n).size2/(MetalObject(n).centerl-MetalObject(n).centeri+MetalObject(n).size1/2)    !tempsize2(y)
	  	temp1=(MetalObject(n).centerl-i)*MetalObject(n).size3/(MetalObject(n).centerl-MetalObject(n).centeri+MetalObject(n).size1/2)  
	        	if(MetalObject(n).centerj-temp/2<=j .and.  MetalObject(n).centerj+temp/2>=j .and. MetalObject(n).centerj-temp1/2<=k .and. MetalObject(n).centerj+temp1/2>=k ) then
		    	in_MetalObject= 1
			else
			in_MetalObject= 0
			endif
	          else
	 		in_MetalObject=0

		  endif
	elseif(MetalObject(n).cshape=="rod") then
	 
		if( (i-MetalObject(n).centeri)*(i-MetalObject(n).centeri)+(j-MetalObject(n).centerj)*(j-MetalObject(n).centerj)<=MetalObject(n).size1*MetalObject(n).size1 .and. (MetalObject(n).centerk-0.5*MetalObject(n).size2)<=k .and. k<=(MetalObject(n).centerk+0.5*MetalObject(n).size2) ) then 
		     in_MetalObject= 1 
		else
		     in_MetalObject= 0 
		endif


	else
	      continue
	endif
end function in_MetalObject


end module inputobj 
