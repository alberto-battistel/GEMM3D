PROGRAM grid3D_irregular

IMPLICIT NONE 
INTEGER, PARAMETER	:: dp = kind(1.d0)
REAL(dp), PARAMETER	:: pi = 4*datan(1.d0), mi0 = 4.d0*pi*1.d-7
COMPLEX(dp), PARAMETER	:: ci = dcmplx(0.d0, 1.d0)

REAL(dp), ALLOCATABLE, DIMENSION(:)	:: Esph, Xhet, Yhet, Zhet, Xinhet, Yinhet, Zinhet, rho_het, rho_1D
REAL(dp), ALLOCATABLE, DIMENSION(:)	:: rho_e,  ip1_het, ip2_het, ip3_het, ip4_het, pip5_het, pip6_het
REAL(dp), ALLOCATABLE, DIMENSION(:)	:: x_node, y_node, z_node, x_obs, y_obs, z_obs
REAL(dp)				:: f, dx_obs, dy_obs, dz_obs
REAL(dp)	         	        :: xi_grid, xf_grid, yi_grid, yf_grid, zi_grid, zf_grid 


INTEGER , ALLOCATABLE, DIMENSION(:)	:: Ninh, NXinhet, NYinhet, NZinhet, indobs, eleobs, no_obs
INTEGER , ALLOCATABLE, DIMENSION(:,:)	:: melem, Marestas, mnosaresta, Mblcsvznhs
INTEGER , ALLOCATABLE, DIMENSION(:)	:: vbordas, vlabel, nelebloc, elebloc, vnobs, label, IPELE
INTEGER 				:: Nelem, Nnode, Narestas, Nbordas, Nh, Nhet, Pinhet
INTEGER 				:: NXhet, NYhet, NZhet, Nobs, Nlinhas, INPUTIP

 CHARACTER(LEN=100)			:: nome_arqs
 CHARACTER(LEN=300)			:: nome_path


 call init_dados_grid3D( xi_grid, xf_grid, yi_grid, yf_grid, zi_grid, zf_grid, Nh, Esph, Ninh, INPUTIP, &
			 ip1_het, ip2_het, ip3_het, ip4_het, pip5_het, pip6_het, rho_1D, Nhet, rho_het, x_obs, y_obs, z_obs, Xhet, &
			 Zhet, Yhet, Nxhet, Nyhet, Nzhet, Pinhet, nome_arqs, nome_path )
 if( Pinhet == 1 )then
	 call init_nodein_grid3D( Xhet, Zhet, Yhet, Xinhet, Yinhet, Zinhet, NXhet, NYhet, &
							 NZhet, NXinhet, NYinhet, NZinhet )
 end if

 call grid3D( xi_grid, xf_grid, yi_grid, yf_grid, zi_grid, zf_grid, Nh, Esph, Ninh, Nhet, Xhet, Yhet, &
		Zhet, NXhet, NYhet, NZhet, Xinhet, Yinhet, Zinhet, NXinhet, NYinhet, NZinhet, x_obs, y_obs, &
		z_obs, INPUTIP, IPELE,  ip1_het, ip2_het, ip3_het, ip4_het, pip5_het, pip6_het, Nobs, Pinhet, rho_het, rho_1D, nelebloc, &
		elebloc, indobs, eleobs, no_obs, melem, rho_e, vlabel, label, Marestas, Mblcsvznhs, mnosaresta, vbordas, &
		Nelem, Nnode, Narestas, Nbordas, x_node, y_node, z_node, nome_arqs, nome_path )
 

 deallocate( rho_het, rho_1D, rho_e, Marestas, mnosaresta, Mblcsvznhs )
 deallocate( vlabel, vbordas, Melem, x_node, y_node, z_node )!
 deallocate( nelebloc, elebloc, indobs, eleobs, no_obs, IPELE )! 
 if( INPUTIP == 1 .or. INPUTIP == 2  )then
	 deallocate( ip1_het, ip2_het, ip3_het )
	 if( INPUTIP == 1  )deallocate( ip4_het )
	 if( INPUTIP == 3  )deallocate( ip4_het, pip5_het, pip6_het )
 end if

CONTAINS

!====================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!====================================================================================================

SUBROUTINE init_dados_grid3D( xi_grid, xf_grid, yi_grid, yf_grid, zi_grid, zf_grid, Nh, Esph, Ninh, INPUTIP, &
			      pip1, pip2, pip3, pip4, pip5, pip6, rho_1D, Nhet, rho_het, x_obs, y_obs, z_obs, Xhet, &
			      Zhet, Yhet, Nxhet, Nyhet, Nzhet, Pinhet, nome_arqs, nome_path )
IMPLICIT NONE
REAL(dp)	         	   , INTENT(out) :: xi_grid, xf_grid, yi_grid, yf_grid, zi_grid, zf_grid 
REAL(dp), ALLOCATABLE, DIMENSION(:), INTENT(out) :: Esph, Xhet, Yhet, Zhet, rho_het, rho_1D
REAL(dp), ALLOCATABLE, DIMENSION(:), INTENT(out) :: x_obs, y_obs, z_obs, pip1, pip2, pip3, pip4, pip5, pip6
INTEGER ,			     INTENT(out) :: NXhet, NYhet, NZhet, Nh, Nhet, Pinhet, INPUTIP
INTEGER , ALLOCATABLE, DIMENSION(:), INTENT(out) :: Ninh
 CHARACTER(LEN=100)		   , INTENT(out) :: nome_arqs
 CHARACTER(LEN=300)		   , INTENT(out) :: nome_path

REAL(dp) :: xi, yi, zi
INTEGER	 :: i, j, k(2)

open( unit=20 , file='model_input.dat', status='old' , action='read' )

read(20,*) nome_path
read(20,*) nome_arqs

read(20,*) xi_grid, xf_grid
read(20,*) yi_grid, yf_grid
read(20,*) zi_grid, zf_grid

read(20,*) Nh
if( Nh > 1 )then

	Allocate( Esph(Nh-1), rho_1D(Nh), Ninh(Nh-1) )
	do i=1, Nh
		read(20,*) rho_1D(i)
	end do
	do i=1, Nh-1
		read(20,*) Esph(i)
	end do
	do i=1, Nh-1
		read(20,*) Ninh(i)
	end do

else if( Nh == 1 )then

	Allocate( Esph(Nh), rho_1D(Nh), Ninh(Nh) )
	read(20,*) rho_1D(1)
	read(20,*) Esph(1)
	read(20,*) Ninh(1)

end if


read(20,*) Nobs
allocate( x_obs(Nobs), y_obs(Nobs), z_obs(Nobs) )

x_obs = 0.d0
y_obs = 0.d0
z_obs = 0.d0

do i=1, Nobs
	read(20,*) x_obs(i), y_obs(i), z_obs(i)
end do

read(20,*) NXhet, NYhet, NZhet, Pinhet, INPUTIP




if( Pinhet == 0 )then

	allocate( Xhet(NXhet), Yhet(NYhet), Zhet(NZhet) )


	Xhet = 0.d0
	Yhet = 0.d0
	Zhet = 0.d0

	do i=1, NXhet
		read(20,*) Xhet(i) 
	end do
	
	do i=1, NYhet
		read(20,*) Yhet(i)
	end do
	
	do i=1, NZhet
		read(20,*) Zhet(i)
	end do

elseif( Pinhet == 1 )then

	allocate( Xhet(NXhet), Yhet(NYhet), Zhet(NZhet) )
	allocate( NXinhet(NXhet-1), NYinhet(NYhet-1), NZinhet(NZhet-1) )

	Xhet = 0.d0
	Yhet = 0.d0
	Zhet = 0.d0

	NXinhet = 0
	NXinhet = 0
	NXinhet = 0

	do i=1, NXhet-1
		read(20,*) Xhet(i), NXinhet(i) 
	end do
	read(20,*) Xhet(NXhet)

	do i=1, NYhet-1
		read(20,*) Yhet(i), NYinhet(i) 
	end do
	read(20,*) Yhet(NYhet)

	do i=1, NZhet-1
		read(20,*) Zhet(i), NZinhet(i) 
	end do
	read(20,*) Zhet(NZhet)

end if

Nhet = (NXhet-1)*(NYhet-1)*(NZhet-1)

allocate( rho_het(Nhet) )
rho_het = 0.d0


if( INPUTIP == 0 )then

	do i=1, Nhet
		read(20,*) rho_het(i)
	end do	
elseif( INPUTIP == 1 )then
	allocate(  pip1(Nhet), pip2(Nhet), pip3(Nhet), pip4(Nhet) )
	pip1 = 0.d0
	pip2 = 0.d0
	pip3 = 0.d0
	pip4 = 0.d0
	do i=1, Nhet
		read(20,*) rho_het(i), pip1(i), pip2(i), pip3(i), pip4(i)
	end do	

elseif( INPUTIP == 2 )then
	allocate(  pip1(Nhet), pip2(Nhet), pip3(Nhet) )
	pip1 = 0.d0
	pip2 = 0.d0
	pip3 = 0.d0
	do i=1, Nhet
		read(20,*) rho_het(i), pip1(i), pip2(i), pip3(i)
	end do	
elseif( INPUTIP == 3 )then
	allocate(  pip1(Nhet), pip2(Nhet), pip3(Nhet), pip4(Nhet), pip5(Nhet), pip6(Nhet) )
	pip1 = 0.d0
	pip2 = 0.d0
	pip3 = 0.d0
	pip4 = 0.d0
	pip5 = 0.d0
	pip6 = 0.d0
	do i=1, Nhet
		read(20,*) rho_het(i), pip1(i), pip2(i), pip3(i), pip4(i), pip5(i), pip6(i)
	end do	
end if

 close(20)

END SUBROUTINE init_dados_grid3D

!====================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!====================================================================================================

SUBROUTINE init_nodein_grid3D( Xhet, Zhet, Yhet, Xinhet, Yinhet, Zinhet, NXhet, NYhet, &
							NZhet, NXinhet, NYinhet, NZinhet )

REAL(dp), ALLOCATABLE, DIMENSION(:), INTENT(out) :: Xinhet, Yinhet, Zinhet
REAL(dp), 	       DIMENSION(:), INTENT(in)  :: Xhet, Yhet, Zhet
INTEGER ,	       DIMENSION(:), INTENT(in)  :: NXinhet, NYinhet, NZinhet
INTEGER ,			     INTENT(in)  :: NXhet, NYhet, NZhet

REAL(dp) :: dxin, dyin, dzin
INTEGER	 :: i, j, k, cont(2)

	cont = 0

	allocate( Xinhet(sum(NXinhet)), Yinhet(sum(NYinhet)), Zinhet(sum(NZinhet)) )

	print*, 'Xinhet'
	cont(1) = 1	
	do i=1, NXhet-1

		cont(2) = sum(NXinhet(1:i))

		k = 0
		dxin = (Xhet(i+1)-Xhet(i))/(NXinhet(i)+1)
		do j=cont(2), cont(1), -1
			k = k + 1 
			Xinhet(j) = Xhet(i+1) - k*dxin
		end do
		print*, i, Xinhet(cont(1):cont(2))	
		cont(1) = cont(2) + 1
	end do

	print*, 'Yinhet'
	cont(1) = 1
	do i=1, NYhet-1

		cont(2) = sum(NYinhet(1:i))

		k = 0
		dyin = (Yhet(i+1)-Yhet(i))/(NYinhet(i)+1)
		do j=cont(2), cont(1),-1
			k = k + 1 
			Yinhet(j) = Yhet(i+1) - k*dyin
		end do
		print*, i, Yinhet(cont(1):cont(2))
		cont(1) = cont(2) + 1
	end do

	print*, 'Zinhet'
	cont(1) = 1
	do i=1, NZhet-1

		cont(2) = sum(NZinhet(1:i))

		k = 0
		dzin = (Zhet(i+1)-Zhet(i))/(NZinhet(i)+1)
		do j=cont(2), cont(1), -1
			k = k + 1 
			Zinhet(j) = Zhet(i+1) - k*dzin
		end do
		print*, i, Zinhet(cont(1):cont(2))
		cont(1) = cont(2) + 1
	end do

END SUBROUTINE init_nodein_grid3D

!====================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!====================================================================================================

SUBROUTINE grid3D( xi_grid, xf_grid, yi_grid, yf_grid, zi_grid, zf_grid, Nh, Esph, Ninh, Nhet, Xhet, &
		   Yhet, Zhet, NXhet, NYhet, NZhet, Xinhet, Yinhet, Zinhet, NXinhet, NYinhet, NZinhet, &
		   x_obs, y_obs, z_obs, INPUTIP, IPELE, ip1_het, ip2_het, ip3_het, ip4_het, pip5_het, pip6_het, Nobs, in_out, &
		   rho_het, rho_1D, nelebloc, elebloc, indobs, eleobs, no_obs, melem, rho_e, vlabel, &
		   label, Marestas, Mblcsvznhs, mnosaresta, vbordas, Nelem, Nnode, Narestas, Nbordas, x_node, y_node, &
		   z_node, nome_arqs, nome_path )

IMPLICIT NONE 
 REAL(dp), ALLOCATABLE, DIMENSION(:)  , INTENT(out) :: rho_e
 REAL(dp), ALLOCATABLE, DIMENSION(:)  , INTENT(out) :: x_node, y_node, z_node
 INTEGER , ALLOCATABLE, DIMENSION(:,:), INTENT(out) :: melem, Marestas, mnosaresta, Mblcsvznhs
 INTEGER , ALLOCATABLE, DIMENSION(:)  , INTENT(out) :: vbordas, vlabel, nelebloc, IPELE
 INTEGER , ALLOCATABLE, DIMENSION(:)  , INTENT(out) :: elebloc, indobs, eleobs, no_obs
 INTEGER ,	      		        INTENT(out) :: Nelem, Nnode, Narestas, Nbordas 

 REAL(dp),	         INTENT(in)	:: xi_grid, xf_grid, yi_grid, yf_grid, zi_grid, zf_grid
 REAL(dp), DIMENSION(:), INTENT(in)	:: x_obs, y_obs, z_obs, ip1_het, ip2_het, ip3_het, ip4_het, pip5_het, pip6_het
 REAL(dp), DIMENSION(:), INTENT(in)	:: Esph, Xhet, Yhet, Zhet, Xinhet, Yinhet, Zinhet, rho_het, rho_1D
 INTEGER,	         INTENT(in)	:: NXhet, NYhet, NZhet, Nh, Nhet, Nobs, in_out, INPUTIP
 INTEGER,  DIMENSION(:), INTENT(in)	:: Ninh, NXinhet, NYinhet, NZinhet, label
 CHARACTER(LEN=100)    , INTENT(in)	:: nome_arqs
 CHARACTER(LEN=300)    , INTENT(in)	:: nome_path

 REAL(dp), ALLOCATABLE, DIMENSION(:) :: x, y, z


 REAL(dp)		:: xi_body, xf_body, yi_body, yf_body, zf_body, zi_body
 REAL(dp)		:: dx(2), dy(2), dz(5), w, rho_ar
 REAL(dp)		:: xc_e, yc_e, zc_e
 COMPLEX(dp)		:: skin, skin_ar
 INTEGER		:: Nx, Ny, Nz, Npx, Npy, Npz, i, j, k, l, m, n, cont
 INTEGER		:: indx, indy, indz, borda(3), num(16), nskin, nele_b, max_ind

 rho_ar = 1.d12


 xi_body = Xhet(1)
 yi_body = Yhet(1)
 zi_body = Zhet(1)

 xf_body = Xhet(Nxhet)
 yf_body = Yhet(Nyhet) 
 zf_body = Zhet(Nzhet)

!========================== vector of coordinate X ===========================================

 call VCOORDXY( Xhet, Xinhet, xi_grid, xf_grid, xi_body, xf_body, x_obs, & 
					NXinhet, Nx, Nobs, NXhet, in_out, x )

!========================== vector of coordinate Y ===========================================

 call VCOORDXY( Yhet, Yinhet, yi_grid, yf_grid, yi_body, yf_body, y_obs, & 
					NYinhet, Ny, Nobs, NYhet, in_out, y )

!========================== vector of coordinate z ===========================================

 call VCOORDZ_new( Zhet, Zinhet, Esph, zi_grid, zf_grid, z_obs, &
		     NZinhet, Ninh, Nz, Nobs, NZhet, in_out, Nh, z)

!=======================================================================================
!  Starts the element vector per inversion block.
!=======================================================================================

 call INITVEC_ELEM_BLOC( Nhet, NXhet, NYhet, NZhet, Xhet, Yhet, Zhet, x, y, z, &
					 Nx, Ny, Nz, nelebloc, elebloc )

!=======================================================================================
!  Initiates the element vector per node of the observations.
!=======================================================================================

 call INITVEC_ELEM_OBS( Nobs, indobs, eleobs )

!=======================================================================================
!  array of nodes.
!=======================================================================================

 call Mnode( Nx, Ny, Nz, Nobs, x, y, z, x_obs, y_obs, z_obs, xi_grid, xf_grid, yi_grid, &
		yf_grid, zi_grid, zf_grid, vlabel, label, no_obs, x_node, y_node, z_node, Nnode ) 

!=======================================================================================
!  array of nodes per element and elements per inverse block.
!=======================================================================================

 call Mele_Mbloc_Mobs( Nx, Ny, Nz, Nelem, Nhet, NXhet, NYhet, NZhet, Nh, no_obs,  INPUTIP, & 
		       IPELE,  ip1_het, x_node, y_node, z_node, Xhet, Yhet, Zhet, rho_1D, &
		       rho_ar, rho_e, rho_het, Esph, Melem, nelebloc, elebloc, eleobs, indobs )

!=======================================================================================
!  Matrix of edges per element.
!=======================================================================================

 call Maresta( Nx, Ny, Nz, Nelem, Marestas )

!=======================================================================================
!  Matrix of neighboring cells by cell.
!=======================================================================================

 call Mblocs_vznhs_cell( Nhet, NXhet-1, NYhet-1, NZhet-1, Mblcsvznhs )

!=======================================================================================
! nodes from the edge.
!=======================================================================================

 call arestas( Nx, Ny, Nz, Narestas, mnosaresta )

!=======================================================================================
! Edge on boundary Vector.
!=======================================================================================

 call vfronteira( NX, Ny, Nz, Nbordas, vbordas )

!=======================================================================================
! Vector of edges for observations and derivatives.
!=======================================================================================

 call escrev_saida( x_node, y_node, z_node, vlabel, Melem, Marestas, Mblcsvznhs, mnosaresta, vbordas, &
		    elebloc, eleobs, nelebloc, indobs, no_obs, rho_het, rho_e, INPUTIP, &
		    IPELE, ip1_het, ip2_het, ip3_het, ip4_het, pip5_het, pip6_het, Nnode, Nelem, Narestas, &
		    Nbordas, Nobs, Nhet, nome_arqs, nome_path )

!=======================================================================================
! release local memory. 
!=======================================================================================

print*,'============================================================================'
print*,'============================================================================'
print*,' Number of points along the x-axis ', size(x)
print*, x
print*,'============================================================================'
print*,' Number of points along the y-axis ', size(y)
print*, y
print*,'============================================================================'
print*,' Number of points along the z-axis ', size(z)
print*, z
print*,'============================================================================'
print*,'============================================================================'
print*, 'Nº measures    -->', Nobs
print*, 'nodes per axis -->', Nx, Ny, Nz
print*, 'total of node  -->', Nnode
print*, 'Nº elements    -->', Nelem
print*, 'Nº Edges       -->', Narestas
print*,'============================================================================'
print*, 'xi_grid, xf_grid ', xi_grid, xf_grid
print*, 'yi_grid, yf_grid ', yi_grid, yf_grid
print*, 'zi_grid, zf_grid ', zi_grid, zf_grid
print*,'============================================================================'
print*, 'xi_body, xf_body ', xi_body, xf_body 
print*, 'yi_body, yf_body ', yi_body, yf_body
print*, 'zi_body, zf_body ', zi_body, zf_body
print*,'============================================================================'

deallocate( x, y, z )

END SUBROUTINE grid3D

!=======================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=======================================================================================

SUBROUTINE VCOORDXY( Xhet, Xinhet, xi_grid, xf_grid, xi_body, xf_body, x_obs, & 
					NXinhet, Nx, Nobs, NXhet, in_out, x )
IMPLICIT NONE
REAL(dp), ALLOCATABLE, DIMENSION(:), INTENT(out)   :: x
REAL(dp),	       DIMENSION(:), INTENT(in)    :: Xhet, Xinhet, x_obs
REAL(dp),			     INTENT(in)    :: xi_grid, xf_grid, xi_body, xf_body
INTEGER ,	       DIMENSION(:), INTENT(in)    :: NXinhet
INTEGER ,			     INTENT(in)    :: Nobs, NXhet, in_out
INTEGER ,			     INTENT(inout) :: Nx

INTEGER , ALLOCATABLE, DIMENSION(:) :: ind_aux
REAL(dp), ALLOCATABLE, DIMENSION(:) :: x_aux, x_aux1, x_aux2, xtemp
REAL(dp) :: val_aux, dx_obs_aux, an(4), q(2), soma(2), interv(2), dxbody, x1, x2, tol=1.d-10
INTEGER	 :: i, j, k, indx, Nx_aux(2), ndx(2), cont(2), Ndxbody

if( xi_grid > xi_body )then
  print*,'One of the initial coordinates that defines the boundary'
  print*,'of Dirichlet is greater than that of heterogeneity'		
  stop
elseif( xf_grid < xf_body)then
  print*,'One of the final coordinates that defines the boundary'
  print*,'of Dirichlet is lower than that of heterogeneity'		
  stop
elseif( x_obs(1) < xi_grid .or. x_obs(Nobs) > xf_grid )then
 print*,'There are coordinates of observations outside'
 print*,'of the GRID region'
 stop 
end if

!=================================================================
! 
! Defining the vector with the coordinates of the INVESRSION GRID
! and their respective user-defined discretizations.
! 
!=================================================================

if( in_out == 0 )then

	Nx = NXhet
	allocate( xtemp(Nx) )
	xtemp = 0.d0

	indx = indx + 1
	xtemp(indx:indx+NXhet-1) = Xhet(:)
	an(2) = Xhet(2)-Xhet(1)

else if( in_out == 1 )then

	Nx = NXhet + sum( NXinhet )
	allocate( xtemp(Nx) )
	xtemp = 0.d0

	indx = indx + 1 
	xtemp(indx) = Xhet(1)
	k = 1
	j = 0 
	do i=1, NXhet-1
		
		j = j + NXinhet(i) 
		xtemp(indx+1:indx+NXinhet(i)) = Xinhet(k:j)
		indx = indx + NXinhet(i) + 1
		xtemp(indx) = Xhet(i+1)
		k = j + 1

	end do

end if

!=================================================================
! 
! Defining the vector with the coordinates of the external 
! discretization to the INVERSION GRID, as well as the final 
! vector with all the coordinates along the horizontal direction
! in which it is being discretized.
! 
!=================================================================

	an(1) = xtemp(2)-xtemp(1) 
	an(2) = xtemp(Nx)-xtemp(Nx-1) 

if( dabs(xi_grid-xi_body) < 5.d-2*dabs(an(1)) .and. & 
    dabs(xf_grid-xf_body) < 5.d-2*dabs(an(2)) )then

	Nx = Nx + Nobs
	allocate(x_aux(Nx))
	x_aux(1:Nobs)    = x_obs
	x_aux(Nobs+1:Nx) = xtemp
	deallocate(xtemp)

	call ordenar_x( x_aux, Nx )
	call eliminar_valrept( x_aux, x_obs, x, Nx )
	deallocate( x_aux )
	return

elseif( dabs(xi_grid-xi_body) < 5.d-2*an(1) )then		

	q(1)  = 1.2d0
	q(2)  = 1.5d0
	call discretizar_interv( xf_grid, xf_body, an(2), q, x_aux2 ) 

	do j=1, size(x_aux2)
		x_aux2(j) = x_aux2(j) + xf_body
	end do
	Nx_aux(1) = Nx + Nobs
	Nx        = Nx_aux(1) + size(x_aux2)

	allocate( x_aux(Nx) )
	x_aux(1:Nobs)	        = x_obs
	x_aux(Nobs+1:Nx_aux(1)) = xtemp
	x_aux(Nx_aux(1)+1:Nx)   = x_aux2

	call ordenar_x( x_aux, Nx )
	call eliminar_valrept( x_aux, x_obs, x, Nx )
	deallocate( x_aux )
	return

elseif( dabs(xf_grid-xf_body) < 5.d-2*an(2) )then		

	q(1)  = 1.2d0
	q(2)  = 1.5d0
	call discretizar_interv( xi_grid, xi_body, an(1), q, x_aux1 ) 

	do j=1, size(x_aux1)
		x_aux1(j) = -x_aux1(j) + xi_body
	end do

	Nx_aux(1) = Nx+Nobs
	Nx        = Nx_aux(1) + size(x_aux1)

	allocate( x_aux(Nx) )
	x_aux(1:Nobs)	        = x_obs
	x_aux(Nobs+1:Nx_aux(1)) = xtemp
	x_aux(Nx_aux(1)+1:Nx)   = x_aux1

	call ordenar_x( x_aux, Nx )
	call eliminar_valrept( x_aux, x_obs, x, Nx )
	deallocate( x_aux )
	return
else

!	print*,' OK!'
	q(1)  = 1.2d0
	q(2)  = 1.5d0
	call discretizar_interv( xi_grid, xi_body, an(1), q, x_aux1 ) 
	call discretizar_interv( xf_grid, xf_body, an(2), q, x_aux2 ) 

	do j=1, size(x_aux1)
		x_aux1(j) = -x_aux1(j) + xi_body
	end do

	do j=1, size(x_aux2)
		x_aux2(j) = x_aux2(j) + xf_body
	end do

	Nx_aux(1) = size(x_aux1) + Nobs
	Nx_aux(2) = Nx + Nx_aux(1)
	Nx        = Nx_aux(2) + size(x_aux2)

	allocate( x_aux(Nx) )

	x_aux(1:Nobs)                = x_obs
	x_aux(Nobs+1:Nx_aux(1))      = x_aux1
	x_aux(Nx_aux(1)+1:Nx_aux(2)) = xtemp
	x_aux(Nx_aux(2)+1:Nx)	     = x_aux2

	call ordenar_x( x_aux, Nx )
	call eliminar_valrept( x_aux, x_obs, x, Nx )
	deallocate( x_aux )
	return

end if

END SUBROUTINE VCOORDXY

!=======================================================================================
!=======================================================================================

	SUBROUTINE discretizar_interv( x_grid, x_body, an1, q, xaux )
	IMPLICIT NONE
	REAL(dp), ALLOCATABLE, DIMENSION(:), INTENT(out) :: xaux
	REAL(dp) 			   , INTENT(in)  :: x_grid, x_body, an1, q(:)
	
	REAL(dp), ALLOCATABLE, DIMENSION(:) :: xaux1, xaux2
	REAL(dp) :: soma(2), interv(2), an(2)
	INTEGER  :: i, ndx(2)
	
		an(1) = an1
		soma  = 0.d0
		ndx(1) = 0
		interv(1) =  2.d-1*dabs( x_grid-x_body )
		interv(2) = (1.d0 - 2.d-1)*dabs( x_grid-x_body )
	
		do 
			soma(1) = soma(1) + an(1)*q(1)**ndx(1)
			if( soma(1) < interv(1) )then
	
				if( dabs(soma(1)-soma(2)) < dabs(interv(1)-soma(1)) )then				 
					ndx(1)  = ndx(1) + 1	
		 			soma(2) = soma(1)
				elseif( ndx(1) == 0 )then
	
					ndx(1) = 1
					allocate( xaux1(1) )
					xaux1(1) = interv(1)	
					an(2)    = interv(1)
					exit						
				else
					ndx(1)  = ndx(1) + 1	
					allocate( xaux1(ndx(1)) )
					if( ndx(1) == 2 )then
						xaux1(1) = an(1)
						xaux1(2) = interv(1)
					else
						xaux1(1) = an(1)
						do i=2, ndx(1)-1
							xaux1(i) = xaux1(i-1) + an(1)*q(1)**(i-1)
						end do
						xaux1(ndx(1)) = interv(1)
						an(2)    = dabs(xaux1(ndx(1)) - xaux1(ndx(1)-1))
					end if
					exit
				end if
			else
				if( .not. allocated(xaux1) )then

					ndx(1)  = ndx(1) + 1	
					allocate( xaux1(ndx(1)) )
					if( ndx(1) == 2 )then
						xaux1(1) = an(1)
						xaux1(2) = interv(1)
					else
						xaux1(1) = an(1)
						do i=2, ndx(1)-1
							xaux1(i) = xaux1(i-1) + an(1)*q(1)**(i-1)
						end do
						xaux1(ndx(1)) = interv(1)
						an(2)    = dabs(xaux1(ndx(1)) - xaux1(ndx(1)-1))
					end if
					exit
				end if
			end if	
		end do
		
		soma = 0
		ndx(2) = 0
		do 
			soma(1) = soma(1) + an(2)*q(2)**ndx(2)
			if( soma(1) < interv(2) )then
	
				if( dabs(soma(1)-soma(2)) < dabs(interv(2)-soma(1)) )then				 
					ndx(2)  = ndx(2) + 1	
		 			soma(2) = soma(1)
				elseif( ndx(2) == 0 )then
	
					ndx(2) = 1
					allocate( xaux2(1) )
					xaux2(1) = interv(2)
					exit				
				else
					ndx(2)  = ndx(2) + 1	
					allocate( xaux2(ndx(2)) )
					if( ndx(2) == 2)then
						xaux2(1) = xaux1(ndx(1))+an(2)
						xaux2(2) = sum(interv)
					else
						xaux2(1) = xaux1(ndx(1)) + an(2)
						do i=2, ndx(2)-1
							xaux2(i) = xaux2(i-1) + an(2)*q(2)**(i-1)
						end do
						xaux2(ndx(2)) = sum(interv)

					end if	
					exit
				end if
			else
				if( .not. allocated(xaux2) )then
					ndx(2)  = ndx(2) + 1	
					allocate( xaux2(ndx(2)) )
					if( ndx(2) == 2)then
						xaux2(1) = xaux1(ndx(1))+an(2)
						xaux2(2) = sum(interv)
					else
						xaux2(1) = xaux1(ndx(1)) + an(2)
						do i=2, ndx(2)-1
							xaux2(i) = xaux2(i-1) + an(2)*q(2)**(i-1)
						end do
						xaux2(ndx(2)) = sum(interv)

					end if	
					exit
				end if
			end if
	
		end do
	
		allocate( xaux(ndx(1)+ndx(2)) )
	
		xaux(1:ndx(1)) = xaux1
		xaux(ndx(1)+1:ndx(1)+ndx(2)) = xaux2

		deallocate( xaux1, xaux2 ) 
	
	END SUBROUTINE

!=======================================================================================
!=======================================================================================

	SUBROUTINE ordenar_x( x1, Nx1 )
	IMPLICIT NONE
	REAL(dp), DIMENSION(:), INTENT(inout) :: x1
	INTEGER		      , INTENT(inout) :: Nx1

	REAL(dp) :: X_aux(Nx1), val_aux
	INTEGER  :: i, j

	X_aux = x1
	do i=1, Nx1-1
		val_aux = 0
		do j=1+i, Nx1
	
			if( X_aux(j) < X_aux(i) )then
				val_aux  = X_aux(i)
				X_aux(i) = X_aux(j)
				X_aux(j) = val_aux
			end if
		end do
	end do
	X1 = x_aux
	X_aux = x1
	
	END SUBROUTINE

!=======================================================================================
!=======================================================================================

	SUBROUTINE eliminar_valrept( x1, obs, x2, Nx1 )
	IMPLICIT NONE
	REAL(dp), ALLOCATABLE, DIMENSION(:), INTENT(inout) :: x2
	REAL(dp), DIMENSION(:)		   , INTENT(inout) :: x1
	REAL(dp), DIMENSION(:)		   , INTENT(in)    :: obs	
	INTEGER		      		   , INTENT(inout) :: Nx1

	INTEGER, ALLOCATABLE, DIMENSION(:) :: ind_aux
	INTEGER  :: i, k, cont(2), sizeobs

	allocate(ind_aux(Nx1))	
	sizeobs = size(obs)
	cont = 0
	ind_aux = 0

	do i=1, Nx1-1
	
		if( x1(i) == x1(i+1) )then
			cont(1) = cont(1) + 1
			ind_aux(i+1) = 1
		else if( x1(i+1) < 0.d0 )then
			if( x1(i)/x1(i+1) < (1.d0+1.d-4) )then
				cont(1) = cont(1) + 1
				do k=1, sizeobs
					if( x1(i) == obs(k) )then
						ind_aux(i+1) = 1
						exit
					else if( x1(i+1) == obs(k) )then
						ind_aux(i)   = 1
						exit
					else
						ind_aux(i+1) = 1
					end if	
				end do			
			end if
		else if( x1(i+1) == 0.d0 .and. x1(i) < 0.d0 )then
			if( dabs(x1(i)) < (1.d-4) )then
				cont(1) = cont(1) + 1
				do k=1, sizeobs
					if( x1(i) == obs(k) )then
						ind_aux(i+1) = 1
						exit
					else if( x1(i+1) == obs(k) )then
						ind_aux(i)   = 1
						exit
					else
						ind_aux(i+1) = 1
					end if	
				end do
			end if	
		else if( x1(i+1) > 0.d0 .and. x1(i) < 0.d0 )then
			if( x1(i+1) < (1.d-4) .and. dabs(x1(i)) < (1.d-4) )then
				cont(1) = cont(1) + 1
				do k=1, sizeobs
					if( x1(i) == obs(k) )then
						ind_aux(i+1) = 1
						exit
					else if( x1(i+1) == obs(k) )then
						ind_aux(i)   = 1
						exit
					else
						ind_aux(i+1) = 1
					end if	
				end do
			end if
		else if( x1(i+1) > 0.d0 .and. x1(i) > 0.d0 )then
			if( x1(i+1)/x1(i) < (1.d0+1.d-4) )then
				cont(1) = cont(1) + 1
				do k=1, sizeobs
					if( x1(i) == obs(k) )then
						ind_aux(i+1) = 1
						exit
					else if( x1(i+1) == obs(k) )then
						ind_aux(i)   = 1
						exit
					else
						ind_aux(i+1) = 1
					end if	
				end do			
			end if
		else if( x1(i+1) > 0.d0 .and. x1(i) == 0.d0 )then
			if( x1(i+1) < 1.d-4 )then
				cont(1) = cont(1) + 1
				do k=1, sizeobs
					if( x1(i)   == obs(k) )then
						ind_aux(i+1) = 1
						exit
					else if( x1(i+1) == obs(k) )then
						ind_aux(i)   = 1
						exit
					else
						ind_aux(i+1) = 1
					end if	
				end do			
			end if
	
		end if

	end do


	allocate( x2( Nx1-cont(1) ) )
	do i=1, Nx1
		if( ind_aux(i) == 0 )then
		cont(2) = cont(2) + 1
		x2(cont(2)) = x1(i) 
		endif
	end do	

	Nx1 = Nx1-cont(1)
	deallocate(ind_aux)
	END SUBROUTINE

!=======================================================================================
!=======================================================================================

!=======================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=======================================================================================

SUBROUTINE VCOORDZ_new( Zhet, Zinhet, Esph, zi_grid, zf_grid, z_obs, &
			NZinhet, Ninh, Nz, Nobs, NZhet, in_out, Nh, z)

IMPLICIT NONE
REAL(dp), ALLOCATABLE, DIMENSION(:), INTENT(out)   :: z
REAL(dp), 	       DIMENSION(:), INTENT(in)    :: Zhet, Zinhet, Esph, z_obs
REAL(dp),			     INTENT(in)    :: zi_grid, zf_grid
INTEGER , 			     INTENT(inout) :: Nz
INTEGER , 	       DIMENSION(:), INTENT(in)    :: NZinhet, Ninh
INTEGER , 			     INTENT(in)    :: NZhet, Nobs, in_out, Nh!, nskin

REAL(dp), ALLOCATABLE, DIMENSION(:) :: z_aux1, z_aux2, z_aux3, z_aux4
INTEGER , ALLOCATABLE, DIMENSION(:) :: ind_aux

REAL(dp) :: val_aux, dz_aux, dz_obs, q(3), an(2), soma(2), dzbody, zi_body, zf_body, tol0=1.d-5, tol1=1.d-3!, tol2=1.d-5
INTEGER	 :: i, j, k(2), indz, Nz_aux(2), cont(2), ndz(3), cnt, Ndzbody

Ndzbody = 2
dzbody  = Zhet(NZhet) - Zhet(1)

!=================================================================
! Defining the vector with the coordinates of the layers (strata) 
! and their respective discretizations defined by the user.
!=================================================================

if( Nh > 1 )then

	Nz_aux(1) = Nh + sum( Ninh )
	allocate( z_aux1(Nz_aux(1)) )
	z_aux1 = 0.d0

	indz = 1
	do i=1, Nh-1

		dz_aux = Esph(i)/(Ninh(i)+1)
		do j=1, Ninh(i) + 1
			indz = indz + 1
			z_aux1(indz) = z_aux1(indz-1) + dz_aux
		end do
	end do

else

	Nz_aux(1) = Nh
	allocate( z_aux1(Nz_aux(1)) )
	z_aux1 = 0.d0

end if


!=================================================================
! Defining the vector with the coordinates of the INVERSION GRID 
! and their respective user-defined discretizations.
!=================================================================

if( in_out == 0 )then

	Nz_aux(2) = NZhet 
	allocate( z_aux2(Nz_aux(2)) )
	z_aux2 = 0.d0
	z_aux2 = Zhet

else if( in_out == 1 )then 

	Nz_aux(2) = NZhet + sum(NZinhet)
	allocate( z_aux2(Nz_aux(2)) )
	z_aux2 = 0.d0

	indz = 0
	k(1) = 1
	do i = 1, NZhet-1
	
		indz = indz + 1
		z_aux2(indz) = Zhet(i)
		if( NZinhet(i) == 0 )cycle
		k(2) = sum(NZinhet(1:i))
		do j=k(1), k(2)
			indz = indz + 1
			z_aux2(indz) = Zinhet(j)
		end do
		k(1) = k(2)+1
	
	end do
	z_aux2(Nz_aux(2)) = Zhet(NZhet)

end if

 Nz = sum(Nz_aux)
 allocate( z_aux3(Nz+Nobs), ind_aux(Nz+Nobs) )
 z_aux3 = 0.d0
 z_aux3(1:Nz_aux(1))    = z_aux1
 z_aux3(Nz_aux(1)+1:Nz) = z_aux2
 z_aux3(Nz+1:Nz+Nobs)   = z_obs

 Nz= Nz+Nobs

 call ordenar_x( z_aux3, Nz )
 call eliminar_valrept( z_aux3, z_obs, z, Nz )

 deallocate( z_aux1, z_aux2, z_aux3 )
 allocate( z_aux1(Nz) )

 z_aux1 = z 
 deallocate( z )

!===============================================================================
! adding points of the air, after the z of heterog and in the final 
! region of the mesh.
!===============================================================================

 q(1)  = 1.8d0
 q(2)  = 2.0d0
 an(1) = z_aux1(2)
 call discretizar_interv( zi_grid, z_aux1(1), an(1), q, z_aux2 ) 
 
 do j=1, size(z_aux2)
	z_aux2(j) = -z_aux2(j) + z_aux1(1)
 end do

 q(1)  = 1.5d0
 q(2)  = 2.0d0
 an(2) = z_aux1(Nz)-z_aux1(Nz-1)
 call discretizar_interv( z_aux1(Nz), zf_grid, an(2), q, z_aux3 ) 

 do j=1, size(z_aux3)
	z_aux3(j) = z_aux3(j) + z_aux1(Nz)
 end do

 allocate( z_aux4(Nz + size(z_aux2) + size(z_aux3)) )
 z_aux4(1:Nz)                 = z_aux1
 z_aux4(Nz+1:Nz+size(z_aux2)) = z_aux2
 z_aux4(Nz+size(z_aux2)+1:Nz+size(z_aux2)+size(z_aux3)) = z_aux3
 Nz = Nz + size(z_aux2) + size(z_aux3)
 deallocate( z_aux1, z_aux2, z_aux3 ) 

 call ordenar_x( z_aux4, Nz )
 call eliminar_valrept( z_aux4, z_obs, z, Nz )
 deallocate( z_aux4 )

END SUBROUTINE VCOORDZ_new

!=======================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=======================================================================================
SUBROUTINE INITVEC_ELEM_BLOC( Nhet, NXhet, NYhet, NZhet, Xhet, Yhet, Zhet, x, y, z, &
							Nx, Ny, Nz, nelebloc, elebloc )
 	
IMPLICIT NONE
INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(out) :: nelebloc, elebloc
REAL(dp), 	      DIMENSION(:), INTENT(in)  :: Xhet, Yhet, Zhet, x, y, z
INTEGER,			    INTENT(in)  :: NXhet, NYhet, NZhet
INTEGER,			    INTENT(in)  :: Nx, Ny, Nz, Nhet

INTEGER, ALLOCATABLE, DIMENSION(:) :: Nxin, Nyin, Nzin
INTEGER :: i, j, k, cont

allocate( Nxin(NXhet-1), Nyin(NYhet-1), Nzin(NZhet-1 ), nelebloc(Nhet+1) )


k = 1
do i=1, NXhet-1
	cont = 0
	do j=k, Nx
		k=j
		if( Xhet(i) < x(j) .and. x(j) < Xhet(i+1) )then
			cont = cont + 1		
		end if	
		if( x(j) > Xhet(i+1) )exit
	end do
	Nxin(i) = cont+1
end do

k = 1
do i=1, NYhet-1
	cont = 0
	do j=k, Ny
		k=j
		if( Yhet(i) < y(j) .and. y(j) < Yhet(i+1) )then
			cont = cont + 1		
		end if	
		if( y(j) > Yhet(i+1) )exit
	end do
	Nyin(i) = cont+1
end do

k = 1
do i=1, NZhet-1
	cont = 0
	do j=k, Nz
		k=j
		if( Zhet(i) < z(j) .and. z(j) < Zhet(i+1) )then
			cont = cont + 1		
		end if	
		if( z(j) > Zhet(i+1) )exit
	end do
	Nzin(i) = cont+1
end do

nelebloc = 0
cont	 = 1
do i=1, NZhet-1
	do j=1, NYhet-1
		do k=1, NXhet-1
			cont = cont + 1
			nelebloc(cont) = Nxin(k)*Nyin(j)*Nzin(i)
		end do
	end do
end do

nelebloc(1)=1
do i=2, Nhet+1
	nelebloc(i) = sum(nelebloc(i-1:i))
end do

cont = nelebloc(Nhet+1)-1
allocate( elebloc(cont) )
elebloc = 0

END SUBROUTINE INITVEC_ELEM_BLOC

!=======================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=======================================================================================

SUBROUTINE INITVEC_ELEM_OBS( Nobs, indobs, eleobs )
IMPLICIT NONE
INTEGER , ALLOCATABLE, DIMENSION(:), INTENT(out) :: indobs, eleobs
INTEGER ,	      		     INTENT(in)	 :: Nobs

INTEGER	:: i, j, neleobs, num

neleobs = 8
num = Nobs*neleobs
allocate( indobs(Nobs+1), eleobs(num) )

indobs = 8
eleobs = 0

indobs(1)=1
do i=2, Nobs+1
	indobs(i) = sum(indobs(i-1:i))
end do

END SUBROUTINE INITVEC_ELEM_OBS

!=======================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=======================================================================================


SUBROUTINE Mnode( Nx, Ny, Nz, Nobs, x, y, z, x_obs, y_obs, z_obs, xi_grid, xf_grid, yi_grid, &
		  yf_grid, zi_grid, zf_grid, vlabel, label, no_obs, x_node, y_node, z_node, Nnode ) 

IMPLICIT NONE
REAL(dp), ALLOCATABLE, DIMENSION(:), INTENT(out) :: x_node, y_node, z_node
INTEGER , ALLOCATABLE, DIMENSION(:), INTENT(out) :: no_obs, vlabel
INTEGER ,	      		     INTENT(out) :: Nnode

REAL(dp),	       DIMENSION(:), INTENT(in)  :: x, y, z, x_obs, y_obs, z_obs
REAL(dp),	       		     INTENT(in)  :: xi_grid, xf_grid, yi_grid
REAL(dp),	       		     INTENT(in)  :: yf_grid, zi_grid, zf_grid
INTEGER ,	      		     INTENT(in)  :: Nx, Ny, Nz, Nobs, label(:)

INTEGER	:: l, i, j, k, p, cont(2)

Nnode = Nx*Ny*Nz
allocate( x_node(Nnode), y_node(Nnode), z_node(Nnode), vlabel(Nnode), no_obs(Nobs) )

x_node = 0.d0
y_node = 0.d0
z_node = 0.d0
no_obs = 0
vlabel = 0
cont   = 0
p      = 0

do l=1, nz
	do j=1, ny
		do i=1, nx
			cont(1) = cont(1) + 1
			 if( x(i) == xi_grid .or. x(i) == xf_grid .or. &
 		 	     y(j) == yi_grid .or. y(j) == yf_grid .or. & 
		    	     z(l) == zi_grid .or. z(l) == zf_grid )then

				x_node(cont(1)) = x(i) ; y_node(cont(1)) = y(j) ; z_node(cont(1)) = z(l) ; vlabel(cont(1)) = 1
			else

				x_node(cont(1)) = x(i) ; y_node(cont(1)) = y(j) ; z_node(cont(1)) = z(l) ; vlabel(cont(1)) = 0

			end if
			
			if( cont(2) == Nobs )cycle
			do k = p + 1, Nobs

				if( x(i)==x_obs(k) .and. y(j)==y_obs(k) .and. z(l)==z_obs(k) )then
					no_obs(k) = cont(1)
					cont(2) = cont(2) + 1						
				end if
			end do

		end do
	end do
	
end do

END SUBROUTINE Mnode

!=======================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=======================================================================================

SUBROUTINE Mele_Mbloc_Mobs( Nx, Ny, Nz, Nelem, Nhet, NXhet, NYhet, NZhet, Nh, no_obs, INPUTIP, & 
		            IPELE,  ip_het, x_node, y_node, z_node, Xhet, Yhet, Zhet, rho_1D, &
			    rho_ar, rho_e, rho_het, Esph, Melem, nelebloc, elebloc, eleobs, indobs )

IMPLICIT NONE
REAL(dp), ALLOCATABLE, DIMENSION(:)  , INTENT(out) :: rho_e
INTEGER , ALLOCATABLE, DIMENSION(:,:), INTENT(out) :: Melem
INTEGER ,	      		       INTENT(out) :: Nelem
INTEGER , ALLOCATABLE, DIMENSION(:),   INTENT(out) :: IPELE

REAL(dp),	       DIMENSION(:),   INTENT(in)    :: ip_het
REAL(dp),	       DIMENSION(:),   INTENT(in)    :: x_node, y_node, z_node, rho_1D
REAL(dp),	       DIMENSION(:),   INTENT(in)    :: Xhet, Yhet, Zhet, rho_het, Esph
REAL(dp),	       		       INTENT(inout) :: rho_ar

INTEGER ,	       DIMENSION(:),   INTENT(inout) :: elebloc, eleobs, nelebloc, indobs, no_obs
INTEGER ,	      		       INTENT(in)    :: Nx, Ny, Nz, NXhet, INPUTIP
INTEGER ,	      		       INTENT(in)    :: NYhet, NZhet, Nhet, Nh


INTEGER , ALLOCATABLE, DIMENSION(:) :: no_obs_aux
REAL(dp)			    :: xc_e, yc_e, zc_e, zinf, zsup, NXNY, NXhNYh
REAL(dp)			    :: xi_body, xf_body, yi_body
REAL(dp)	 		    :: yf_body, zi_body, zf_body
INTEGER				    :: l, i, j, n, m, k, p, cont, cont1, blocj

Nelem = (Nx-1)*(Ny-1)*(Nz-1)
allocate( rho_e(Nelem), Melem(Nelem, 8), no_obs_aux(Nobs) )
allocate( IPELE(Nelem) )

no_obs_aux = 0
rho_e   = 0.d0
Melem   = 0
elebloc = 0
cont    = 0
cont1   = 1
NXNY    = Nx*Ny
NXhNYh  = (NXhet-1)*(NYhet-1)
IPELE   = 0

do l=1, Nz-1
	do j=1, Ny-1
		do i=1, Nx-1

			!=========================================================================
			! Mapping nos by elements.
			!=========================================================================

			cont = cont + 1
			Melem(cont,1) = (l-1)*NXNY + (j-1)*Nx + i
			Melem(cont,2) = (l-1)*NXNY + (j-1)*Nx + i + 1
			Melem(cont,3) = (l-1)*NXNY + j*Nx + i + 1
			Melem(cont,4) = (l-1)*NXNY + j*Nx + i 
			Melem(cont,5) = l*NXNY + (j-1)*Nx + i 
			Melem(cont,6) = l*NXNY + (j-1)*Nx + i + 1
			Melem(cont,7) = l*NXNY + j*Nx + i + 1
			Melem(cont,8) = l*NXNY + j*Nx + i
		
			!=========================================================================
			! Mapping the elements that share points of observation.
			!=========================================================================

			if( cont1 <= Nobs )then
				do m=1, Nobs
					do n=1, 8
						
						if( Melem(cont,n) == no_obs(m) )then

							if( no_obs_aux(m) == 1 )exit!Go to 20
							no_obs_aux(m) = 1

							eleobs(indobs(m))   = cont
							eleobs(indobs(m)+1) = cont + 1
							eleobs(indobs(m)+2) = cont + Nx-1
							eleobs(indobs(m)+3) = cont + Nx
							eleobs(indobs(m)+4) = cont + (Nx-1)*(Ny-1)
							eleobs(indobs(m)+5) = cont + (Nx-1)*(Ny-1) + 1
							eleobs(indobs(m)+6) = cont + (Nx-1)*(Ny-1) + Nx - 1
							eleobs(indobs(m)+7) = cont + (Nx-1)*(Ny-1) + Nx
							cont1 = cont1 + 1
							do k=1, Nobs
								if( no_obs(m) == no_obs(k) .and. k /= m )then
									eleobs(indobs(k):indobs(k)+7) = eleobs(indobs(m):indobs(m)+7) 
									no_obs_aux(k) = 1
									cont1 = cont1 + 1
								end if	
							end do
							exit	
						end if
					end do
				20 end do
			end if

			!=============================================================================
			! Mapping the elements belonging to their respective Inversion Blocks as well 
                        ! as mapping the resistivity of each Grid element.
			!=============================================================================

			xc_e = sum( x_node(Melem(cont,:)) )/8
			yc_e = sum( y_node(Melem(cont,:)) )/8
			zc_e = sum( z_node(Melem(cont,:)) )/8
		
			xi_body = Xhet(1)	; xf_body = Xhet(NXhet)
			yi_body = Yhet(1)	; yf_body = Yhet(NYhet) 
			zi_body = Zhet(1)	; zf_body = Zhet(NZhet)  
		
			if( xc_e > xi_body .and. xc_e < xf_body .and. &
			    yc_e > yi_body .and. yc_e < yf_body .and. &
			    zc_e > zi_body .and. zc_e < zf_body )then
		
				do m=1, NZhet-1				
		
					zi_body = Zhet(m)
					zf_body = Zhet(m+1)
					if( zc_e < zi_body .or. zc_e > zf_body )cycle
		
					do n=1, NYhet-1						
		
						yi_body = Yhet(n)
						yf_body = Yhet(n+1)
						if( yc_e < yi_body .or. yc_e > yf_body )cycle
		
						do k=1, NXhet-1						
		
							xi_body = Xhet(k)
							xf_body = Xhet(k+1)
		
							if( xc_e > xi_body .and. xc_e < xf_body )then

								blocj	    = (m-1)*NXhNYh + &
									      ((NXhet-1)*(n-1)) + k
								rho_e(cont) = rho_het(blocj)

								if( INPUTIP > 0 )then
									if( ip_het(blocj) /= 0.d0 )then
										IPELE(cont) = 1
									end if
								end if
								do p=nelebloc(blocj), nelebloc(blocj+1)-1
									if( elebloc(p) == 0 )then
										elebloc(p) = cont
										go to 10
									end if
								end do
					
							end if
						end do
					end do
				end do
			
			else if( zc_e < 0.d0 )then
				rho_e(cont) = rho_ar
			else if( zc_e > 0.d0 )then

				if( Nh > 1 )then
					zinf = 0.d0					
					do k=1, Nh-1
						zsup = sum(Esph(1:k))
						if( zc_e > zinf .and. zc_e < zsup )then		
							rho_e(cont) = rho_1D(k)
							go to 10
						else
							zinf = zsup
						end if
					end do
					if( rho_e(cont) == 0.d0)rho_e(cont) = rho_1D(Nh)

				else
					rho_e(cont) = rho_1D(Nh)		
				end if

			10 end if
		end do
	end do
end do

END SUBROUTINE Mele_Mbloc_Mobs

!=======================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=======================================================================================

SUBROUTINE Maresta( Nx, Ny, Nz, Nelem, Marestas )

IMPLICIT NONE 
INTEGER, ALLOCATABLE, DIMENSION(:,:), INTENT(out) :: Marestas
INTEGER, 			      INTENT(IN)  :: Nx, Ny, Nz, Nelem  

INTEGER :: i, j, l, cont, num(2)

allocate( Marestas(Nelem,12) )
cont = 0

do l=1, Nz-1
	do j=1, Ny-1
		do i=1, Nx-1

			cont = cont + 1		
			Marestas(cont,1) = (l-1)*(Nx-1)*Ny + (j-1)*(Nx-1) + i
			Marestas(cont,2) = (l-1)*(Nx-1)*Ny + j*(Nx-1) + i
			Marestas(cont,3) = l*(Nx-1)*Ny + (j-1)*(Nx-1) + i 
			Marestas(cont,4) = l*(Nx-1)*Ny + j*(Nx-1) + i 

			Marestas(cont,5) = (l-1)*Nx*(Ny-1) + (i-1)*(Ny-1) + j 
			Marestas(cont,6) = l*Nx*(Ny-1) + (i-1)*(Ny-1) + j
			Marestas(cont,7) = (l-1)*Nx*(Ny-1) + i*(Ny-1) + j 
			Marestas(cont,8) = l*Nx*(Ny-1) + i*(Ny-1) + j

			Marestas(cont,9)  = (l-1)*Nx*Ny + (j-1)*Nx + i 
			Marestas(cont,10) = (l-1)*Nx*Ny + (j-1)*Nx + i + 1
			Marestas(cont,11) = (l-1)*Nx*Ny + j*Nx + i 
			Marestas(cont,12) = (l-1)*Nx*Ny + j*Nx + i + 1
		end do
	end do
end do

num(1) = (Nx-1)*Ny*Nz
num(2) = num(1) + Nx*(Ny-1)*Nz

Marestas(:,5:8) = Marestas(:,5:8) + (Nx-1)*Ny*Nz
Marestas(:,9:12) = Marestas(:,9:12) + (Nx-1)*Ny*Nz + Nx*(Ny-1)*Nz

END SUBROUTINE Maresta

!=======================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=======================================================================================

SUBROUTINE Mblocs_vznhs_cell( Nhet, Nhetx, Nhety, Nhetz, Mblcsvznhs )
IMPLICIT NONE
 INTEGER, ALLOCATABLE, DIMENSION(:,:), INTENT(out)  :: Mblcsvznhs
 INTEGER                             , INTENT(in)   :: Nhetx, Nhety, Nhetz, Nhet

 INTEGER :: k, i, j, indc, indc_aux


 ALLOCATE( Mblcsvznhs(Nhet,6) )

 Mblcsvznhs = 0
 indc_aux   = Nhetx*Nhety 
 print*, indc_aux


if( Nhetx > 1 .and. Nhety > 1 .and. Nhetz > 1 )then

	 indc = 0
	 do k=1, Nhetz-1
		do j=1, Nhety
			do i=1, Nhetx
				indc = indc + 1
	 			Mblcsvznhs( indc_aux+indc, 5) = indc
	 			Mblcsvznhs( indc, 6)          = indc_aux + indc
			end do
		end do
	 end do

	 indc = 0
	 do k=1, Nhetz
		do j=1, Nhety-1
			do i=1, Nhetx
				indc = indc + 1
	 			Mblcsvznhs( indc+Nhetx, 3) = indc
	 			Mblcsvznhs( indc, 4)       = Nhetx + indc	
			end do
		end do
		indc = indc + Nhetx
	 end do

	 indc_aux = 1
	 indc = 0
	 do k=1, Nhetz
		do j=1, Nhety

	       	        indc = indc + 1
			Mblcsvznhs( indc, 2) = indc + 1
			if( (Nhetx-1) > 1)then
				do i=2, Nhetx-1
					indc = indc + 1
		 			Mblcsvznhs( indc, 1) = indc - 1
		 			Mblcsvznhs( indc, 2) = indc + 1
			 	end do
			end if
	       	        indc = indc + 1
			Mblcsvznhs( indc, 1) = indc - 1
		end do
	 end do

elseif( Nhetx == 1 .and. Nhety > 1 .and. Nhetz > 1 )then

	 indc = 0
	 do k=1, Nhetz-1
		do j=1, Nhety
			indc = indc + 1
	 		Mblcsvznhs( indc_aux+indc, 5) = indc
	 		Mblcsvznhs( indc, 6)          = indc_aux + indc
		end do
	 end do

	 indc = 0
	 do k=1, Nhetz
		do j=1, Nhety-1
			indc = indc + 1
	 		Mblcsvznhs( indc+Nhetx, 3) = indc
	 		Mblcsvznhs( indc, 4)       = Nhetx + indc	
		end do
		indc = indc + Nhetx
	 end do


elseif(  Nhetx > 1 .and. Nhety == 1 .and. Nhetz > 1  )then

	 indc = 0
	 do k=1, Nhetz-1
		do i=1, Nhetx
			indc = indc + 1
			Mblcsvznhs( indc_aux+indc, 5) = indc
			Mblcsvznhs( indc, 6)          = indc_aux + indc
		end do
	 end do

	 indc_aux = 1
	 indc = 0
	 do k=1, Nhetz
       	        indc = indc + 1
		Mblcsvznhs( indc, 2) = indc + 1
		if( (Nhetx-1) > 1)then
			do i=2, Nhetx-1
				indc = indc + 1
	 			Mblcsvznhs( indc, 1) = indc - 1
	 			Mblcsvznhs( indc, 2) = indc + 1
		 	end do
		end if
       	        indc = indc + 1
		Mblcsvznhs( indc, 1) = indc - 1
	 end do


elseif( Nhetx > 1 .and. Nhety > 1 .and. Nhetz == 1 )then


	 indc = 0
	do j=1, Nhety-1
		do i=1, Nhetx
			indc = indc + 1
 			Mblcsvznhs( indc+Nhetx, 3) = indc
 			Mblcsvznhs( indc, 4)       = Nhetx + indc	
		end do
	end do
	indc = indc + Nhetx


	indc_aux = 1
	indc = 0
	do j=1, Nhety

       	        indc = indc + 1
		Mblcsvznhs( indc, 2) = indc + 1
		if( (Nhetx-1) > 1)then
			do i=2, Nhetx-1
				indc = indc + 1
	 			Mblcsvznhs( indc, 1) = indc - 1
	 			Mblcsvznhs( indc, 2) = indc + 1
		 	end do
		end if
       	        indc = indc + 1
		Mblcsvznhs( indc, 1) = indc - 1
	end do

elseif( Nhetx > 1 .and. Nhety == 1 .and. Nhetz == 1 )then

	 indc_aux = 1
         indc = 1
	 Mblcsvznhs( indc, 2) = indc + 1
	if( (Nhetx-1) > 1)then
		do i=2, Nhetx-1
			indc = indc + 1
			Mblcsvznhs( indc, 1) = indc - 1
			Mblcsvznhs( indc, 2) = indc + 1
	 	end do
	end if
        indc = indc + 1
	Mblcsvznhs( indc, 1) = indc - 1

elseif( Nhetx == 1 .and. Nhety > 1 .and. Nhetz == 1 )then

	 indc = 0
	 do k=1, Nhetz
		do j=1, Nhety-1
			indc = indc + 1
 			Mblcsvznhs( indc+Nhetx, 3) = indc
 			Mblcsvznhs( indc, 4)       = Nhetx + indc	
		end do
		indc = indc + Nhetx
	 end do

elseif( Nhetx == 1 .and. Nhety == 1 .and. Nhetz > 1 )then

	 indc = 0
	 do k=1, Nhetz-1
		indc = indc + 1
		Mblcsvznhs( indc_aux+indc, 5) = indc
		Mblcsvznhs( indc, 6)          = indc_aux + indc
	 end do

end if

 END SUBROUTINE Mblocs_vznhs_cell

!=======================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=======================================================================================

SUBROUTINE arestas( Nx, Ny, Nz, Narestas, mnosaresta )

IMPLICIT NONE
INTEGER, ALLOCATABLE, DIMENSION(:,:), INTENT(out) :: mnosaresta
INTEGER, 			      INTENT(out) :: Narestas
INTEGER, 			      INTENT(in)  :: Nx, Ny, Nz

INTEGER :: l, i, j, borda(3), num(2)

Narestas = (Nx-1)*Ny*Nz + Nx*(Ny-1)*Nz + Ny*Nx*(Nz-1)
allocate( mnosaresta(Narestas,2) )

num(1) = (Nx-1)*Ny*Nz
num(2) = num(1) + Nx*(Ny-1)*Nz

do l=1, Nz
	do  j=1, Ny
		do i=1, Nx-1

			borda(1) = (l-1)*(Nx-1)*Ny + (j-1)*(Nx-1) + i
			mnosaresta( borda(1),1) = borda(1) + (j-1) + (l-1)*Ny
			mnosaresta( borda(1),2) = borda(1) + (j-1) + 1 + (l-1)*Ny

			if( j < Ny )then

				borda(2) = (l-1)*(Ny-1)*Nx + (i-1)*(Ny-1) + j + num(1)			
				mnosaresta( borda(2),1) = borda(1) + (j-1) + (l-1)*Ny
				mnosaresta( borda(2),2) = (l-1)*Nx*Ny + j*Nx + i

				if( i == Nx-1 )then
					borda(2) = (l-1)*(Ny-1)*Nx + i*(Ny-1) + j + num(1)			
					mnosaresta( borda(2),1) = (l-1)*Nx*Ny + j*Nx
					mnosaresta( borda(2),2) = (l-1)*Nx*Ny + (j+1)*Nx
				end if
			end if			
		end do

		if( l < Nz )then
			do i=1, Nx

				borda(1) = (l-1)*Nx*Ny + (j-1)*Nx + i
				borda(3) = (l-1)*Nx*Ny + (j-1)*Nx + i + num(2)							
				mnosaresta( borda(3),1) = borda(1) 
				mnosaresta( borda(3),2) = borda(1) + Nx*Ny 

			end do
		end if

	end do

end do

END SUBROUTINE arestas

!=======================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=======================================================================================

SUBROUTINE vfronteira( NX, Ny, Nz, Nbordas, vbordas )

IMPLICIT NONE
INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(out) :: vbordas
INTEGER,			    INTENT(out) :: Nbordas
INTEGER,			    INTENT(in)  :: Nx, Ny, Nz

INTEGER	:: l, i, j, cont, num(2)

Nbordas  = 2*( (Nx-1)*Nz + (Nx-1)*(Ny-2) + &
	       (Ny-1)*Nz + (Ny-1)*(Nx-2) + &
	       (Nz-1)*Nx + (Nz-1)*(Ny-2) )

allocate( vbordas(Nbordas) )

num(1) = (Nx-1)*Ny*Nz
num(2) = num(1) + Nx*(Ny-1)*Nz

cont = 0
do l=1, Nz !Nz-1
	do j=1, Ny !Ny-1
		do i=1, Nx-1
			
			if( l==1 .or. l==Nz )then
				cont = cont + 1
				vbordas(cont) = (l-1)*(Nx-1)*Ny + (j-1)*(Nx-1) + i
			else if( j==1 .or. j==Ny )then
				cont = cont + 1
				vbordas(cont) = (l-1)*(Nx-1)*Ny + (j-1)*(Nx-1) + i
			end if

		end do
	end do
end do

do l=1, Nz
	do i=1, Nx !Nx-1
		do j=1, Ny-1
			
			if( l==1 .or. l==Nz )then
				cont = cont + 1
				vbordas(cont) = num(1) + (l-1)*Nx*(Ny-1) + (i-1)*(Ny-1) + j
!				print*, num(1)
			else if( i==1 .or. i==Nx )then
				cont = cont + 1
				vbordas(cont) = num(1) + (l-1)*Nx*(Ny-1) + (i-1)*(Ny-1) + j
			end if

		end do
	end do
end do

do l=1, Nz-1
	do j=1, Ny
		do i=1, Nx !Nx-1
			
			if( j==1 .or. j==Ny )then
				cont = cont + 1
				vbordas(cont) = num(2) + (l-1)*Nx*Ny + (j-1)*Nx + i
			else if( i==1 .or. i==Nx )then
				cont = cont + 1
				vbordas(cont) = num(2) + (l-1)*Nx*Ny + (j-1)*Nx + i
			end if

		end do
	end do
end do

END SUBROUTINE vfronteira

!=======================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=======================================================================================

SUBROUTINE escrev_saida( x_node, y_node, z_node, vlabel, Melem, Marestas, Mblcsvznhs, mnosaresta, vbordas, &
			 elebloc, eleobs, nelebloc, indobs, no_obs, rhohet, rho_e, INPUTIP, &
		    	 IPELE, ip1_het, ip2_het, ip3_het, ip4_het, pip5_het, pip6_het, Nnode, Nelem, Narestas, &
		    	 Nbordas, Nobs, Nhet, nome_arqs, nome_path )


IMPLICIT NONE
REAL(dp),	DIMENSION(:)  , INTENT(in) :: rho_e, ip1_het, ip2_het, ip3_het, ip4_het, pip5_het, pip6_het
REAL(dp),	DIMENSION(:)  , INTENT(in) :: x_node, y_node, z_node, rhohet
INTEGER,	DIMENSION(:,:), INTENT(in) :: Melem, Marestas, Mblcsvznhs, mnosaresta
INTEGER,	DIMENSION(:)  , INTENT(in) :: vlabel, vbordas, elebloc, eleobs, IPELE
INTEGER,	DIMENSION(:)  , INTENT(in) :: nelebloc, indobs, no_obs
INTEGER			      , INTENT(in) :: Nnode, Nelem, Narestas, Nbordas, Nobs, Nhet, INPUTIP
 CHARACTER(LEN=100)	      , INTENT(in) :: nome_arqs, nome_path

INTEGER	:: i


open( unit=30  , file=trim(nome_path)//'/grid3D/'//trim(nome_arqs)//'.node'  , status='replace' , action='write' )
open( unit=40  , file=trim(nome_path)//'/grid3D/'//trim(nome_arqs)//'.ele'   , status='replace' , action='write' )
open( unit=50  , file=trim(nome_path)//'/grid3D/'//trim(nome_arqs)//'.rste'  , status='replace' , action='write' )
open( unit=60  , file=trim(nome_path)//'/grid3D/'//trim(nome_arqs)//'.rsts'  , status='replace' , action='write' )
open( unit=70  , file=trim(nome_path)//'/grid3D/'//trim(nome_arqs)//'.rstb'  , status='replace' , action='write' )
open( unit=80  , file=trim(nome_path)//'/grid3D/'//trim(nome_arqs)//'.eblc'  , status='replace' , action='write' )
open( unit=90  , file=trim(nome_path)//'/grid3D/'//trim(nome_arqs)//'.eobs'  , status='replace' , action='write' )
open( unit=100 , file=trim(nome_path)//'/grid3D/'//trim(nome_arqs)//'.obs'   , status='replace' , action='write' )
open( unit=115 , file=trim(nome_path)//'/grid3D/'//trim(nome_arqs)//'.blcvz' , status='replace' , action='write' )

!======================================================================================
!======================================================================================

write(30,*) Nnode
30 format( 2x, 3d25.15, 1I4 )
do i=1, Nnode
	write(30,30) x_node(i), y_node(i), z_node(i), vlabel(i)
end do

!======================================================================================
!======================================================================================

write(40,*) Nelem
40 format( 2x, 8I10, 1d25.15, 1I4 )

if( INPUTIP == 0 )then
	do i=1, Nelem
		write(40,40) Melem(i,1:8), rho_e(i)
	end do
else
	do i=1, Nelem
		write(40,40) Melem(i,1:8), rho_e(i), IPELE(i)
	end do
end if

!======================================================================================
!======================================================================================

write(50,*) Nelem
50 format( 2x, 12I10 )
do i=1, Nelem
	write(50,50) Marestas(i,1:12)
end do


write(60,*) Narestas
do i=1, Narestas
	write(60,*) mnosaresta(i,1:2)
end do

write(70,*) Nbordas
do i=1, Nbordas
	write(70,*) vbordas(i)
end do


!======================================================================================
!======================================================================================

60 format( 2x, 1I10, 7d30.15  )
write(80,*) Nhet+1
if( INPUTIP == 0 )then
	do i=1, Nhet
		write(80,60) nelebloc(i), rhohet(i)
	end do	
elseif( INPUTIP == 1 )then
	do i=1, Nhet
		write(80,60) nelebloc(i), rhohet(i), ip1_het(i), ip2_het(i), ip3_het(i), ip4_het(i)
	end do
elseif( INPUTIP == 2 )then
	do i=1, Nhet
		write(80,60) nelebloc(i), rhohet(i), ip1_het(i), ip2_het(i), ip3_het(i)
	end do
elseif( INPUTIP == 3 )then
	do i=1, Nhet
		write(80,60) nelebloc(i), rhohet(i), ip1_het(i), ip2_het(i), ip3_het(i), ip4_het(i), pip5_het(i), pip6_het(i)
	end do

end if

do i=1, Nhet
	write(115,40) Mblcsvznhs(i,1:6)
end do


!======================================================================================
!======================================================================================

write(80,60) nelebloc(Nhet+1)
write(80,60) nelebloc(Nhet+1)-1
do i=1, nelebloc(Nhet+1)-1
	write(80,60) elebloc(i)
end do

write(90,*) Nobs+1
do i=1, Nobs+1
	write(90,*) indobs(i)
end do

write(90,*) indobs(Nobs+1)-1
do i=1, indobs(Nobs+1)-1
	write(90,*) eleobs(i)
end do

!======================================================================================
!======================================================================================

write(100,*) Nobs
100 format( 2x, 1I10, 3d25.15 )
do i=1, Nobs
	write(100,100) no_obs(i), x_node(no_obs(i)), y_node(no_obs(i)), z_node(no_obs(i))
end do

!======================================================================================
!======================================================================================

 close(30)
 close(40)
 close(50)
 close(60)
 close(70)
 close(80)
 close(90)
 close(100)

END SUBROUTINE escrev_saida

!=======================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=======================================================================================

END PROGRAM grid3D_irregular
