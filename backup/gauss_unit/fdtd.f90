!--3D FDTD method based on Taflov's book
!--Boundary:UPML and Dirichit boundary 
!--Source:Ex/Hy planewave along z direction,implemented by TF/SF
!--Units: choose E = sqrt(eps0/mu0)*E0(standard units) 
!--		 H = H0(standard units)

program fdtd
implicit none

!-------------------------------------------Parameters Setup--------------------------------------------------
!--default parameters

TF_SF_xl= floor(0.5+((rpos(1,1)+xcenter)*lattice_x))
TF_SF_xr= floor(0.5+((rpos(1,2)+xcenter)*lattice_x))
TF_SF_yl= floor(0.5+((rpos(2,1)+ycenter)*lattice_y))
TF_SF_yr= floor(0.5+((rpos(2,2)+ycenter)*lattice_y))
TF_SF_kl= floor(0.5+((rpos(3,1)+zcenter)*lattice_z))
TF_SF_kr= floor(0.5+((rpos(3,2)+zcenter)*lattice_z))

interval = 200
			

