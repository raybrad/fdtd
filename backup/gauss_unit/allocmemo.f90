module allocmemo

use parameters

implicit none
contains
subroutine memory()
 
integer:: i,j,k
real*8::value

allocate(epsilonx(isize,jsize,ksize),epsilony(isize,jsize,ksize),epsilonz(isize,jsize,ksize),&
	 ddx(isize,jsize,ksize),ddy(isize,jsize,ksize),ddz(isize,jsize,ksize)) 

!For TF-SF formulation 
allocate(Ex_inc(ksize),Hy_inc(ksize),Ey_inc(ksize),Hx_inc(ksize))

!Drude metal
allocate(Lepsilon(isize,jsize,ksize),Lomega(isize,jsize,ksize),Lgamma(isize,jsize,ksize),Molecular(isize,jsize,ksize))

!Field	
allocate(Ex(isize,jsize,ksize),Ey(isize,jsize,ksize),Ez(isize,jsize,ksize),& 
	 Hx(isize,jsize,ksize),Hy(isize,jsize,ksize),Hz(isize,jsize,ksize),& 
         Dx(isize,jsize,ksize),Dy(isize,jsize,ksize),Dz(isize,jsize,ksize),& 
	 Bx(isize,jsize,ksize),By(isize,jsize,ksize),Bz(isize,jsize,ksize),& 
	 Jx(isize,jsize,ksize),Jy(isize,jsize,ksize),Jz(isize,jsize,ksize))

!pml,non-pml distribution	
allocate(cposition(isize,jsize,ksize))

do k=1,ksize
  do j=1,jsize
    do i=1,isize
     if((pmlil+1)<=i .and. i<=(isize-pmlir) .and. (pmljl+1)<=j .and. j<=(jsize-pmljr) .and. (pmlkl+1)<=k .and. k<=(ksize-pmlkr)) then
                 cposition(i,j,k)=1 !non-pml
     else
                 cposition(i,j,k)=0 !pml
     endif
    enddo
  enddo
enddo
write(6,*) 'pml x',1,pmlil,isize-pmlir+1,isize
write(6,*) 'pml y',1,pmljl,jsize-pmljr+1,jsize
write(6,*) 'pml z',1,pmlkl,ksize-pmlkr+1,ksize

allocate(aax(isize),aay(jsize),aaz(ksize),bbx(isize),bby(jsize),bbz(ksize),ccx(isize),ccy(jsize),ccz(ksize),& 
	 eex(isize),eey(jsize),eez(ksize),ffx(isize),ffy(jsize),ffz(ksize),ggx(isize),ggy(jsize),ggz(ksize),& 
	 hhx(isize),hhy(jsize),hhz(ksize),iix(isize),iiy(jsize),iiz(ksize),jjx(isize),jjy(jsize),jjz(ksize),& 
	 kkx(isize),kky(jsize),kkz(ksize),llx(isize),lly(jsize),llz(ksize)) 



write(6,*)"memory...ok" 
end subroutine memory

subroutine field_initialization()
 
	Ex=0.0;Ey=0.0;Ez=0.0 
	Hx=0.0;Hy=0.0;Hz=0.0 
	Dx=0.0;Dy=0.0;Dz=0.0 
	Bx=0.0;By=0.0;Bz=0.0 
	Jx=0.0;Jy=0.0;Jz=0.0 

end subroutine field_initialization

end module allocmemo
