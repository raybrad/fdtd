module allocmemo

use parameters

implicit none
contains
subroutine memory()
 
	integer:: i,j,k
	real*8::value
allocate(epsilonx(0:misize-1,0:mjsize-1,0:mksize-1),&
	epsilony(0:misize-1,0:mjsize-1,0:mksize-1),& 
	epsilonz(0:misize-1,0:mjsize-1,0:mksize-1),&

	ddx (0:misize-1,0:mjsize-1,0:mksize-1),& 
	ddy (0:misize-1,0:mjsize-1,0:mksize-1),& 
	ddz (0:misize-1,0:mjsize-1,0:mksize-1)) 

	epsilonx=0.d0;epsilony=0.d0;epsilonz=0.d0
	ddx=0.d0;ddy=0.d0;ddz=0.d0
! For TF-SF formulation 
	if(TF_SF_freq /= 0.0) then
	 
         allocate(Ex_inc(0:mksize-1),&
	           Hy_inc(0:mksize-1),Ey_inc(0:mksize-1),Hx_inc(0:mksize-1)) 
		   Ex_inc=0.d0
		   Hy_inc=0.d0
	endif

	if(dofft) then
	allocate(Ex_Re(Nf,0:misize-1,0:mjsize-1,0:mksize-1),Ey_Re(Nf,0:misize-1,0:mjsize-1,0:mksize-1), & 
		 Ez_Re(Nf,0:misize-1,0:mjsize-1,0:mksize-1),Hx_Re(Nf,0:misize-1,0:mjsize-1,0:mksize-1), &
		 Hy_Re(Nf,0:misize-1,0:mjsize-1,0:mksize-1),Hz_Re(Nf,0:misize-1,0:mjsize-1,0:mksize-1), &
		 Ex_Im(Nf,0:misize-1,0:mjsize-1,0:mksize-1),Ey_Im(Nf,0:misize-1,0:mjsize-1,0:mksize-1), & 
		 Ez_Im(Nf,0:misize-1,0:mjsize-1,0:mksize-1),Hx_Im(Nf,0:misize-1,0:mjsize-1,0:mksize-1), &
		 Hy_Im(Nf,0:misize-1,0:mjsize-1,0:mksize-1),Hz_Im(Nf,0:misize-1,0:mjsize-1,0:mksize-1), &
         Ex_inc_Re(Nf),Ex_inc_Im(Nf))

		 Ex_Re=0.d0;Ey_Re=0.d0;Ez_Re=0.d0
		 Ex_Im=0.d0;Ey_Im=0.d0;Ez_Im=0.d0
		 Hx_Re=0.d0;Hy_Re=0.d0;Hz_Re=0.d0
		 Hx_Im=0.d0;Hy_Im=0.d0;Hz_Im=0.d0
         Ex_inc_Re=0.d0;Ex_inc_Im=0.d0
	endif

	allocate(lmonitor(0:misize-1,0:mjsize-1,0:mksize-1))
	lmonitor=.false.
	! Normal or Drude metal or Single Lorentz or Drude+ Double Lorentz
allocate(	Lepsilon(0:misize-1,0:mjsize-1,0:mksize-1),& 
		Lomega(0:misize-1,0:mjsize-1,0:mksize-1),& 
		Lgamma(0:misize-1,0:mjsize-1,0:mksize-1),&
		Lsigma(0:misize-1,0:mjsize-1,0:mksize-1),&
		Molecular(0:misize-1,0:mjsize-1,0:mksize-1))
		
		Lepsilon=0.0
		Lomega=0.0
		Lgamma=0.0
		Lsigma=0.0
		Molecular=0.0
	!/*Double Lorentz or Drude +  Double Lorentz pole. */
	if(metal_type==3 .or.  metal_type ==4 .or. metal_type==5 ) then
	 
allocate(	Ldepsr_1(0:misize-1,0:mjsize-1,0:mksize-1),& 
		Lomega_1(0:misize-1,0:mjsize-1,0:mksize-1),& 
		Lgamma_1(0:misize-1,0:mjsize-1,0:mksize-1),& 
		Ldepsr_2(0:misize-1,0:mjsize-1,0:mksize-1),& 
		Lomega_2(0:misize-1,0:mjsize-1,0:mksize-1),& 
		Lgamma_2(0:misize-1,0:mjsize-1,0:mksize-1))

		Ldepsr_1=0.0
		Lomega_1=0.0
		Lgamma_1=0.0
		Ldepsr_2=0.0
		Lomega_2=0.0
		Lgamma_2=0.0
	endif

allocate(Ex(0:misize-1,0:mjsize-1,0:mksize-1),& 
	 Ey(0:misize-1,0:mjsize-1,0:mksize-1),& 
	 Ez(0:misize-1,0:mjsize-1,0:mksize-1),& 
	 Hx(0:misize-1,0:mjsize-1,0:mksize-1),& 
         Hy(0:misize-1,0:mjsize-1,0:mksize-1),& 
	 Hz(0:misize-1,0:mjsize-1,0:mksize-1),& 
	
         Dx(0:misize-1,0:mjsize-1,0:mksize-1),& 
	 Dy(0:misize-1,0:mjsize-1,0:mksize-1),& 
	 Dz(0:misize-1,0:mjsize-1,0:mksize-1),& 
	
	 Bx(0:misize-1,0:mjsize-1,0:mksize-1),& 
	 By(0:misize-1,0:mjsize-1,0:mksize-1),& 
	 Bz(0:misize-1,0:mjsize-1,0:mksize-1)) 

	!/* Drude or single pole Lorentzian (non re-interpreted). */
	if(metal_type==2 .or. metal_type == 3) then
	 
allocate(	Jx(0:misize-1,0:mjsize-1,0:mksize-1),& 
		Jy(0:misize-1,0:mjsize-1,0:mksize-1),& 
		Jz(0:misize-1,0:mjsize-1,0:mksize-1))
		Jx=0.0;Jy=0.0;Jz=0.0
	endif

	if(metal_type == 3) then
allocate(	pJx (0:misize-1,0:mjsize-1,0:mksize-1),& 
		pJy (0:misize-1,0:mjsize-1,0:mksize-1),& 
		pJz (0:misize-1,0:mjsize-1,0:mksize-1))

		pJx=0.0;pJy=0.0;pJz=0.0
	endif

	if(metal_type == 3  .or. metal_type==5) then

allocate(	pEx (0:misize-1,0:mjsize-1,0:mksize-1),& 
		pEy (0:misize-1,0:mjsize-1,0:mksize-1),& 
		pEz (0:misize-1,0:mjsize-1,0:mksize-1))
		
		pEx=0.0;pEy=0.0;pEz=0.0
	endif

	if( metal_type == 4 .or. metal_type == 5) then
allocate(     Jx_1(0:misize-1,0:mjsize-1,0:mksize-1),& 
	      Jy_1(0:misize-1,0:mjsize-1,0:mksize-1),& 
              Jz_1(0:misize-1,0:mjsize-1,0:mksize-1),&
	      Jx_2(0:misize-1,0:mjsize-1,0:mksize-1),& 
	      Jy_2(0:misize-1,0:mjsize-1,0:mksize-1),& 
	      Jz_2(0:misize-1,0:mjsize-1,0:mksize-1),&

	      pJx_1(0:misize-1,0:mjsize-1,0:mksize-1),& 
              pJy_1(0:misize-1,0:mjsize-1,0:mksize-1),& 
	      pJz_1(0:misize-1,0:mjsize-1,0:mksize-1),&

	      pJx_2(0:misize-1,0:mjsize-1,0:mksize-1),& 
              pJy_2(0:misize-1,0:mjsize-1,0:mksize-1),& 
	      pJz_2(0:misize-1,0:mjsize-1,0:mksize-1))

	      Jx_1=0.0;Jy_1=0.0;Jz_1=0.0
	      Jx_2=0.0;Jy_2=0.0;Jz_2=0.0
	     pJx_1=0.0;pJy_1=0.0;pJz_1=0.0
	     pJx_2=0.0;pJy_2=0.0;pJz_2=0.0
	endif	
	if(metal_type==5) then
allocate(     Jx_0(0:misize-1,0:mjsize-1,0:mksize-1),& 
	      Jy_0(0:misize-1,0:mjsize-1,0:mksize-1),& 
              Jz_0(0:misize-1,0:mjsize-1,0:mksize-1))
	      Jx_0=0.0;Jy_0=0.0;Jz_0=0.0
	endif
allocate(cposition(0:misize-1,0:mjsize-1,0:mksize-1))

	do k=0,mksize-1
		do j=0,mjsize-1
			do i=0,misize-1
		if((pmlil+1)<=i .and. i<=(isize-pmlir-2) .and. (pmljl+1)<=j .and. j<=(jsize-pmljr-2) .and. (pmlkl+1)<=k .and. k<=(ksize-pmlkr-2)) then
		                cposition(i,j,k)=1 !non-pml
		else
		                cposition(i,j,k)=0 !pml
		endif
		        enddo
		enddo
	enddo
	write(6,*) 'pml x',0,pmlil,isize-pmlir-1,misize-1
	write(6,*) 'pml y',0,pmljl,jsize-pmljr-1,mjsize-1
	write(6,*) 'pml z',0,pmlkl,ksize-pmlkr-1,mksize-1

allocate(aax(0:misize-1),& 
 	 aay(0:mjsize-1),& 
	 aaz(0:mksize-1),& 
	 bbx(0:misize-1),& 
	 bby(0:mjsize-1),& 
	 bbz(0:mksize-1),& 
	 ccx(0:misize-1),& 
	 ccy(0:mjsize-1),& 
	 ccz(0:mksize-1),& 
	 eex(0:misize-1),& 
	 eey(0:mjsize-1),& 
	 eez(0:mksize-1),& 
	 ffx(0:misize-1),& 
	 ffy(0:mjsize-1),& 
	 ffz(0:mksize-1),& 
	 ggx(0:misize-1),& 
	 ggy(0:mjsize-1),& 
	 ggz(0:mksize-1),& 
	 hhx(0:misize-1),& 
	 hhy(0:mjsize-1),& 
	 hhz(0:mksize-1),& 
	 iix(0:misize-1),& 
	 iiy(0:mjsize-1),& 
	 iiz(0:mksize-1),& 
	 jjx(0:misize-1),& 
	 jjy(0:mjsize-1),& 
	 jjz(0:mksize-1),& 
	 kkx(0:misize-1),& 
	 kky(0:mjsize-1),& 
	 kkz(0:mksize-1),& 
	 llx(0:misize-1),& 
	 lly(0:mjsize-1),& 
	 llz(0:mksize-1)) 
	
	 aax=0.d0;aay=0.d0;aaz=0.d0
	 bbx=0.d0;bby=0.d0;bbz=0.d0
	 ccx=0.d0;ccy=0.d0;ccz=0.d0
	 eex=0.d0;eey=0.d0;eez=0.d0
	 ffx=0.d0;ffy=0.d0;ffz=0.d0
	 ggx=0.d0;ggy=0.d0;ggz=0.d0
	 hhx=0.d0;hhy=0.d0;hhz=0.d0
	 iix=0.d0;iiy=0.d0;iiz=0.d0
	 jjx=0.d0;jjy=0.d0;jjz=0.d0
	 kkx=0.d0;kky=0.d0;kkz=0.d0
	 llx=0.d0;lly=0.d0;llz=0.d0


	print *,"memory...ok" 
end subroutine memory

subroutine field_initialization()
 
	Ex=0.0;Ey=0.0;Ez=0.0 
	Hx=0.0;Hy=0.0;Hz=0.0 
	Dx=0.0;Dy=0.0;Dz=0.0 
	Bx=0.0;By=0.0;Bz=0.0 
	
	Ex_Re=0.d0;Ey_Re=0.d0;Ez_Re=0.d0
	Hx_Re=0.d0;Hy_Re=0.d0;Hz_Re=0.d0
end subroutine field_initialization

end module allocmemo
