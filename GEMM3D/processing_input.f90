MODULE ler_escrever_arq
!=================================================================================================
				USE G_variables 
!=================================================================================================
CONTAINS
!=================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=================================================================================================

SUBROUTINE inici_entradas
IMPLICIT NONE


 CHARACTER(len=100)	:: arq_node, arq_ele, arq_borda, arq_arestas, arq_nos_are
 CHARACTER(len=100)	:: arq_ele_bloc, arq_ele_obs, arq_obs
 INTEGER		:: i, j, k


 call leitura_dados_entrada( ngl, arq_node, arq_ele, arq_borda, arq_arestas, &
			     arq_nos_are, arq_ele_bloc, arq_ele_obs, arq_obs )

 call ler_grid_3D( arq_node, arq_ele, arq_borda, arq_arestas, arq_nos_are )

 call  print_dados_1D2D3D_font( nnosmalha, Nelem, narestas )

 call dados_obs_df( arq_ele_bloc, arq_ele_obs, arq_obs )

 call esclha_filtro
 
END SUBROUTINE inici_entradas

!============================================================================================
! 8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!============================================================================================

SUBROUTINE  esclha_filtro

 ! (3)is the parameter that sets the Kong filter as default with 241 abscissas and weights

   autr_fltJ0 = 3
   numpJ0     = 241
   autr_fltJ1 = 3 
   numpJ1     = 241

END SUBROUTINE  esclha_filtro

!============================================================================================
! 8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!============================================================================================

SUBROUTINE leitura_dados_entrada( ngl, arq_node, arq_ele, arq_borda, arq_arestas, &
				   arq_nos_are, arq_ele_bloc, arq_ele_obs, arq_obs )
IMPLICIT NONE
 INTEGER, INTENT(out)             :: ngl
 CHARACTER(len=100), INTENT(out)  :: arq_node, arq_ele, arq_borda, arq_arestas
 CHARACTER(len=100), INTENT(out)  :: arq_nos_are, arq_ele_bloc, arq_ele_obs, arq_obs

 CHARACTER(len=50) :: nomevar
 INTEGER:: i

	OPEN( unit=100, file='input_data.dat', status='old', action='read' )

	read(100,*) arq_node
	read(100,*) arq_ele
	read(100,*) arq_arestas
	read(100,*) arq_borda
	read(100,*) arq_nos_are
	read(100,*) arq_ele_bloc
	read(100,*) arq_ele_obs
	read(100,*) arq_obs

!	read(100,*) Nfreq, rbch
	read(100,*) ISOURCE, SDIMENS

	read(100,*) Nfreq
	Allocate( vetf(Nfreq) )
 	nomevar = 'Vector of frequencies'
	call alloc_rdvet( vetf, Nfreq, nomevar )

	do i=1, Nfreq
		read(100,*)vetf(i)
	end do


	read(100,*) Ntransm	

 	nomevar = 'Transmitter coordinate vector on x-axis'
	call alloc_rdvet( cxT, Ntransm, nomevar )
 	nomevar = 'Transmitter coordinate vector on y-axis'
	call alloc_rdvet( cyT, Ntransm, nomevar )
 	nomevar = 'Transmitter coordinate vector on z-axis'
	call alloc_rdvet( czT, Ntransm, nomevar )

	do i=1, Ntransm
		read(100,*) cxT(i), cyT(i), czT(i)
	end do


	read(100,*) NC
	if( NC > 1 )then

	 	nomevar = 'resistivity vector of the 1D model'
		call alloc_rdvet( roj, NC, nomevar )
	 	nomevar = 'relative permittivity vector of the 1D model'
		call alloc_rdvet( mirj, NC, nomevar )
	 	nomevar = 'Thick layer vector'
		call alloc_rdvet( hj, NC-1, nomevar )
 
		do i=1, NC
			read(100,*) roj(i)
		end do
		do i=1, NC
			read(100,*) mirj(i)
		end do
	
		do i=1, NC-1
			read(100,*) hj(i)
		end do

	else if( NC == 1 )then
	
	 	nomevar = 'resistivity vector of the 1D model'
		call alloc_rdvet( roj, NC, nomevar )
	 	nomevar = 'relative permittivity vector of the 1D model'
		call alloc_rdvet( mirj, NC, nomevar )
	 	nomevar = 'Thick layer vector'
		call alloc_rdvet( hj, NC, nomevar )

		read(100,*) roj(1)
		read(100,*) mirj(1)
		read(100,*) hj(1)
	
	end if

	read(100,*) Ncorpos, INPUTIP
 	nomevar = 'resistivity vector of the 3D model'
	call alloc_rdvet( rojc, Ncorpos, nomevar )
	if( INPUTIP == 1 )then

	 	nomevar = 'chargeability vector of the 3D model'
		call alloc_rdvet( CARGBLD, Ncorpos, nomevar )
	 	nomevar = 'Relaxation time vector of the 3D model'
		call alloc_rdvet( TMRLX, Ncorpos, nomevar )
	 	nomevar = 'VPortion fraction vector of the 3D model'
		call alloc_rdvet( FCPIP, Ncorpos, nomevar )
	 	nomevar = 'Electrochemical parameter vector of the 3D model'
		call alloc_rdvet( PEQMC, Ncorpos, nomevar )

	elseif( INPUTIP ==2 )then

	 	nomevar = 'chargeability vector of the 3D model'
		call alloc_rdvet( CARGBLD, Ncorpos, nomevar )
	 	nomevar = 'Relaxation time vector of the 3D model'
		call alloc_rdvet( TMRLX, Ncorpos, nomevar )
	 	nomevar = 'frequency dependency parameter vector of the 3D model'
		call alloc_rdvet( CCFREQ, Ncorpos, nomevar )

	elseif( INPUTIP ==3 )then

	 	nomevar = 'chargeability vector of the 3D model'
		call alloc_rdvet( CARGBLD, Ncorpos, nomevar )
	 	nomevar = 'relaxation time vector of the 3D model'
		call alloc_rdvet( TMRLX, Ncorpos, nomevar )
	 	nomevar = 'frequency dependency parameter vector of the 3D model'
		call alloc_rdvet( CCFREQ, Ncorpos, nomevar )

	 	nomevar = '2nd chargeability vector of the 3D model'
		call alloc_rdvet( CARGBLD2, Ncorpos, nomevar )
	 	nomevar = '2nd relaxation time vector of the 3D model'
		call alloc_rdvet( TMRLX2, Ncorpos, nomevar )
	 	nomevar = '2nd frequency dependency parameter vector of the 3D model'
		call alloc_rdvet( CCFREQ2, Ncorpos, nomevar )

	endif

	rojc(:)	= 0.d0
	ngl 	= 1

	call zizf_kmadaj

	close(unit=100)

END SUBROUTINE leitura_dados_entrada

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE ler_grid_3D( arq_node, arq_ele, arq_borda, arq_arestas, arq_nos_are )
 IMPLICIT NONE
 CHARACTER(len=100), INTENT(in)	:: arq_node, arq_ele, arq_borda, arq_arestas, arq_nos_are

 CHARACTER(len=50) :: nomevar
 INTEGER:: l,k

	OPEN( unit=200, file=trim(pathr)//'grid3D/'//arq_node, status='old', action='read' )
	OPEN( unit=250, file=trim(pathr)//'grid3D/'//arq_ele, status='old', action='read' )
	OPEN( unit=300, file=trim(pathr)//'grid3D/'//arq_arestas, status='old', action='read' )
	OPEN( unit=350, file=trim(pathr)//'grid3D/'//arq_borda, status='old', action='read' )
	OPEN( unit=400, file=trim(pathr)//'grid3D/'//arq_nos_are, status='old', action='read' )

	read(200,*) nnosmalha 
	read(250,*) Nelem
	read(300,*) Nelem
	read(350,*) nborda
	read(400,*) Narestas

	nomevar = 'GRID array of nodes'
	call  alloc_rdmat( Mnos, nnosmalha, 3, nomevar )
	nomevar = 'Matrix of elements'
	call  alloc_intmat( Melem, Nelem, 8, nomevar )
	nomevar = 'Element Edges Matrix'
	call  alloc_intmat( Marestas, Nelem, 12, nomevar )
	nomevar = 'Grid Edges Matrix'
	call  alloc_intmat( Mnosaresta, Narestas, 2, nomevar )	
	nomevar = 'edge boundary vector'
	call  alloc_intvet( vbordas, nborda, nomevar )
	nomevar = 'resistivity vector of the model'
	call  alloc_rdvet( vrho_e, nelem, nomevar )

	do l=1, nnosmalha
		read(200,*) (Mnos(l,k), k=1,3)
	end do

	if( INPUTIP == 0 )then
		do l=1, nelem
			read(250,*) (Melem(l,k), k=1,8), vrho_e(l)
		end do
	else	
		nomevar = 'index IP vector of the elements'
		call alloc_intvet( IPELE, Nelem, nomevar )
		do l=1, nelem
			read(250,*) (Melem(l,k), k=1,8), vrho_e(l), IPELE(l)
		end do
	end if

	do l=1, nelem
		read(300,*) (Marestas(l,k), k=1,12)
	end do
	do l=1,nborda
		read(350,*) vbordas(l)
	end do
	do l=1, Narestas
		read(400,*) (mnosaresta(l,k), k=1,2)
	end do

	CLOSE(200) ; CLOSE (250)
	CLOSE(300) ; CLOSE (350)
	CLOSE(400)

END SUBROUTINE ler_grid_3D

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE dados_obs_df( arq_ele_bloc, arq_ele_obs, arq_obs )
IMPLICIT NONE

 CHARACTER(len=100), INTENT(in)	    :: arq_ele_bloc, arq_ele_obs, arq_obs

 INTEGER, ALLOCATABLE, DIMENSION(:) :: inde, indnos
 INTEGER 			    :: i, j, k, cnt(2), nloc
 CHARACTER(len=50)		    :: nomevar


open( unit=1300 , file=trim(pathr)//'grid3D/'//arq_ele_bloc, status='old', action='read' )
open( unit=1400 , file=trim(pathr)//'grid3D/'//arq_ele_obs , status='old', action='read' )
open( unit=1500 , file=trim(pathr)//'grid3D/'//arq_obs     , status='old', action='read' )

read(1300,*) Nbloc
read(1400,*) Nobs

!==========================================================================
!==========================================================================

nomevar = 'number of elements per observation vector'
call  alloc_intvet( indobs, Nobs, nomevar )
nomevar = 'number of elements per inversion GRID block'
call  alloc_intvet( nelebloc, Nbloc, nomevar )


if( INPUTIP == 0 )then
	do i=1, Nbloc-1
		read(1300,*) nelebloc(i), rojc(i)
	end do
elseif( INPUTIP == 1 )then
	do i=1, Nbloc-1
		read(1300,*) nelebloc(i), rojc(i), CARGBLD(i), TMRLX(i), FCPIP(i), PEQMC(i)
	end do
elseif( INPUTIP == 2 )then
	do i=1, Nbloc-1
		read(1300,*) nelebloc(i), rojc(i), CARGBLD(i), TMRLX(i), CCFREQ(i)
	end do
elseif( INPUTIP == 3 )then
	do i=1, Nbloc-1
		read(1300,*) nelebloc(i), rojc(i), CARGBLD(i), TMRLX(i), CCFREQ(i), CARGBLD2(i), TMRLX2(i), CCFREQ2(i)
	end do
end if

read(1300,*) nelebloc(Nbloc)
read(1300,*)  size_eblc
nomevar = 'index vector of the elements per inversion GRID block'
call  alloc_intvet( elebloc, size_eblc, nomevar )

do i=1, size_eblc
	read(1300,*) elebloc(i)
end do

!==========================================================================
!==========================================================================

do i=1, Nobs
	read(1400,*) indobs(i)
end do

read(1400,*)  size_eobs
nomevar = 'index vector of elements per observation point'
call  alloc_intvet( eleobs, size_eobs, nomevar )

do i=1, size_eobs
	read(1400,*) eleobs(i)
end do

read(1500,*) Nobs
allocate( no_obs(Nobs) ) 
nomevar = 'index vector of elements of the observation point'
call  alloc_intvet( no_obs, Nobs, nomevar )

do i=1, Nobs
	read(1500,*) no_obs(i)
end do

Nbloc = Nbloc - 1

 close( 1300 )
 close( 1400 )
 close( 1500 )

END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE READ_primarios_Transm( FDP )
IMPLICIT NONE 
 CHARACTER(LEN=20), INTENT(in) :: FDP

 REAL(dp), allocatable, dimension(:)  :: ExpR_aux, EypR_aux, EzpR_aux
 REAL(dp), allocatable, dimension(:)  :: ExpI_aux, EypI_aux, EzpI_aux    
 INTEGER  :: i, j, k, n, unid, notnun(3)

 CHARACTER(len=50)		    :: nomevar

n = 8*(nelebloc( Nbloc+1)-1)

nomevar = 'vector of the real part of the field Ex 1D of the Transmitter'
call  alloc_rdvet( ExpR_aux, n, nomevar )
nomevar = 'vector of the real part of the field Ey 1D of the do Transmitter'
call  alloc_rdvet( EypR_aux, n, nomevar )
nomevar = 'vector of the real part of the field Ez 1D of the Transmitter'
call  alloc_rdvet( EzpR_aux, n, nomevar )
nomevar = 'vector of the imaginary part of the field Ex 1D of the Transmitter'
call  alloc_rdvet( ExpI_aux, n, nomevar )
nomevar = 'vector of the imaginary part of the field Ey 1D of the Transmitter'
call  alloc_rdvet( EypI_aux, n, nomevar )
nomevar = 'vector of the imaginary part of the field Ez 1D of the Transmitter'
call  alloc_rdvet( EzpI_aux, n, nomevar )



IF( .NOT. ALLOCATED( Exp_trans ) ) THEN

	nomevar = 'Field Vector Ex 1D Transmitter'
	call  alloc_cdvet( Exp_trans, n, nomevar )
	nomevar = 'Field Vector Ey 1D Transmitter'
	call  alloc_cdvet( Eyp_trans, n, nomevar )
	nomevar = 'Field Vector Ez 1D Transmitter'
	call  alloc_cdvet( Ezp_trans, n, nomevar )

END IF

    if( FDP == 'ExDP' )then
	unid = 4000
elseif( FDP == 'EyDP' )then
	unid = 4100
elseif( FDP == 'EzDP' )then
	unid = 4200
elseif( FDP == 'HxDP' )then
	unid = 4300
elseif( FDP == 'HyDP' )then
	unid = 4400
elseif( FDP == 'HzDP' )then
	unid = 4500
end if


Exp_trans = (0.d0,0.d0)
Eyp_trans = (0.d0,0.d0)
Ezp_trans = (0.d0,0.d0)

ExpR_aux = 0.d0
EypR_aux = 0.d0
EzpR_aux = 0.d0

ExpI_aux = 0.d0
EypI_aux = 0.d0
EzpI_aux = 0.d0

notnun  = 0

read(unid) ExpR_aux
read(unid) ExpI_aux
read(unid) EypR_aux
read(unid) EypI_aux
read(unid) EzpR_aux
read(unid) EzpI_aux 

do j=nelebloc(1), n

	if( isnan( ExpR_aux(j) ) .or. isnan( ExpI_aux(j) ) )then
		notnun(1) = notnun(1) + 1
	end if
	if( isnan( EypR_aux(j) ) .or. isnan( EypI_aux(j) ) )then
		notnun(2) = notnun(2) + 1
	end if
	if( isnan( EzpR_aux(j) ) .or. isnan( EzpI_aux(j) ) )then
		notnun(3) = notnun(3) + 1
	end if
		
end do


	if( notnun(1) /= 0 )then
		print*, "Number of NaN in the evaluation of the transmitter's primary Ex:"
		print*, notnun(1)
		print*, ' '
	end if
	if( notnun(2) /= 0 )then
		print*, "Number of NaN in the evaluation of the transmitter's primary Ey:"
		print*, notnun(2)
		print*, ' '
	end if
	if( notnun(3) /= 0 )then
		print*, "Number of NaN in the evaluation of the transmitter's primary Ez:"
		print*, notnun(3)
		print*, ' '
	end if

	if( notnun(1)/= 0 .or. notnun(2)/= 0 .or. notnun(3)/= 0)stop

	Exp_trans = ExpR_aux + ci*ExpI_aux
	Eyp_trans = EypR_aux + ci*EypI_aux
	Ezp_trans = EzpR_aux + ci*EzpI_aux

deallocate( ExpR_aux, EypR_aux, EzpR_aux )
deallocate( ExpI_aux, EypI_aux, EzpI_aux )

END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE READ_primarios_Transm_BCH
IMPLICIT NONE 
 REAL(dp), allocatable, dimension(:)  :: ExpR_aux, EypR_aux
 REAL(dp), allocatable, dimension(:)  :: ExpI_aux, EypI_aux
 INTEGER                              :: i, j, k, n, unid, notnun(2)

 CHARACTER(len=50)                    :: nomevar

n = 8*(nelebloc( Nbloc+1)-1)

nomevar = 'vector of the real part of the field Ex 1D of the Transmitter'
call  alloc_rdvet( ExpR_aux, n, nomevar )
nomevar = 'vector of the real part of the field Ey 1D of the Transmitter'
call  alloc_rdvet( EypR_aux, n, nomevar )
nomevar = 'vector of the imaginary part of the field Ex 1D of the Transmitter'
call  alloc_rdvet( ExpI_aux, n, nomevar )
nomevar = 'vector of the imaginary part of the field Ey 1D of the Transmitter'
call  alloc_rdvet( EypI_aux, n, nomevar )

IF( .NOT. ALLOCATED( Exp_trans ) ) THEN
	nomevar = 'Field Vector Ex 1D Transmitter'
	call  alloc_cdvet( Exp_trans, n, nomevar )
	nomevar = 'Field Vector Ey 1D Transmitter'
	call  alloc_cdvet( Eyp_trans, n, nomevar )
END IF

Exp_trans = (0.d0,0.d0)
Eyp_trans = (0.d0,0.d0)

ExpR_aux = 0.d0
EypR_aux = 0.d0

ExpI_aux = 0.d0
EypI_aux = 0.d0

notnun  = 0

unid=4600

read(unid) ExpR_aux
read(unid) ExpI_aux
read(unid) EypR_aux
read(unid) EypI_aux

do j=nelebloc(1), n

	if( isnan( ExpR_aux(j) ) .or. isnan( ExpI_aux(j) ) )then
		notnun(1) = notnun(1) + 1
	end if
	if( isnan( EypR_aux(j) ) .or. isnan( EypI_aux(j) ) )then
		notnun(2) = notnun(2) + 1
	end if
	
end do

	if( notnun(1) /= 0 )then
		print*, "Number of NaN in the evaluation of the transmitter's primary Ex:"
		print*, notnun(1)
		print*, ' '
	end if
	if( notnun(2) /= 0 )then
		print*, "Number of NaN in the evaluation of the transmitter's primary Ep:"
		print*, notnun(2)
		print*, ' '
	end if

	if( notnun(1)/= 0 .or. notnun(2)/= 0 )stop

	Exp_trans = ExpR_aux + ci*ExpI_aux
	Eyp_trans = EypR_aux + ci*EypI_aux

deallocate( ExpR_aux, EypR_aux )
deallocate( ExpI_aux, EypI_aux )

END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE READ_primarios_ADJ_MT( FDP )
IMPLICIT NONE 
 CHARACTER(LEN=20), INTENT(in) :: FDP

 REAL(dp), allocatable, dimension(:)  :: ExpR_aux, EypR_aux, EzpR_aux
 REAL(dp), allocatable, dimension(:)  :: ExpI_aux, EypI_aux, EzpI_aux    
 INTEGER  :: i, j, k, n, unid, notnun(3)

 CHARACTER(len=50)		    :: nomevar

n = 8*(nelebloc( Nbloc+1)-1)

nomevar = 'vector of the real part of the field Ex 1D adjoint'
call  alloc_rdvet( ExpR_aux, n, nomevar )
nomevar = 'vector of the real part of the field Ey 1D adjoint'
call  alloc_rdvet( EypR_aux, n, nomevar )
nomevar = 'vector of the real part of the field Ez 1D adjoint'
call  alloc_rdvet( EzpR_aux, n, nomevar )
nomevar = 'vector of the imaginary part of the field Ex 1D adjoint'
call  alloc_rdvet( ExpI_aux, n, nomevar )
nomevar = 'vector of the imaginary part of the field Ey 1D adjoint'
call  alloc_rdvet( EypI_aux, n, nomevar )
nomevar = 'vector of the imaginary part of the field Ez 1D adjoint'
call  alloc_rdvet( EzpI_aux, n, nomevar )

IF( .NOT. ALLOCATED( Exp_adj ) ) THEN

	nomevar = 'Field Vector Ex 1D Adjoint'
	call  alloc_cdvet( Exp_adj, n, nomevar )
	nomevar = 'Field Vector Ey 1D Adjoint'
	call  alloc_cdvet( Eyp_adj, n, nomevar )
	nomevar = 'Field Vector Ez 1D Adjoint'
	call  alloc_cdvet( Ezp_adj, n, nomevar )
END IF

    if( FDP == 'ExDP' )then
	unid = 3000
elseif( FDP == 'EyDP' )then
	unid = 3100
elseif( FDP == 'EzDP' )then
	unid = 3200
elseif( FDP == 'HxDP' )then
	unid = 3300
elseif( FDP == 'HyDP' )then
	unid = 3400
elseif( FDP == 'HzDP' )then
	unid = 3500
end if


Exp_adj = (0.d0,0.d0)
Eyp_adj = (0.d0,0.d0)
Ezp_adj = (0.d0,0.d0)

ExpR_aux = 0.d0
EypR_aux = 0.d0
EzpR_aux = 0.d0

ExpI_aux = 0.d0
EypI_aux = 0.d0
EzpI_aux = 0.d0

notnun  = 0

read(unid) ExpR_aux
read(unid) ExpI_aux
read(unid) EypR_aux
read(unid) EypI_aux
read(unid) EzpR_aux
read(unid) EzpI_aux 

do j=nelebloc(1), n

			
	if( isnan( ExpR_aux(j) ) .or. isnan( ExpI_aux(j) ) )then
		notnun(1) = notnun(1) + 1
	end if
	if( isnan( EypR_aux(j) ) .or. isnan( EypI_aux(j) ) )then
		notnun(2) = notnun(2) + 1
	end if
	if( isnan( EzpR_aux(j) ) .or. isnan( EzpI_aux(j) ) )then
		notnun(3) = notnun(3) + 1
	end if
		
end do

	if( notnun(1) /= 0 )then
		print*, "Number of NaN in the evaluation of the receiver's primary Ex:"
		print*, notnun(1)
		print*, ' '
	end if
	if( notnun(2) /= 0 )then
		print*,  "Number of NaN in the evaluation of the receiver's primary Ey:"
		print*, notnun(2)
		print*, ' '
	end if
	if( notnun(3) /= 0 )then
		print*,  "Number of NaN in the evaluation of the receiver's primary Ez:"
		print*, notnun(3)
		print*, ' '
	end if

	if( notnun(1)/= 0 .or. notnun(2)/= 0 .or. notnun(3)/= 0)stop

	Exp_adj = ExpR_aux + ci*ExpI_aux
	Eyp_adj = EypR_aux + ci*EypI_aux
	Ezp_adj = EzpR_aux + ci*EzpI_aux

deallocate( ExpR_aux, EypR_aux, EzpR_aux )
deallocate( ExpI_aux, EypI_aux, EzpI_aux )

END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE READ_primarios_ADJ( FDP )
IMPLICIT NONE 
 CHARACTER(LEN=20), INTENT(in) :: FDP

 REAL(dp), allocatable, dimension(:)  :: ExpR_aux, EypR_aux, EzpR_aux
 REAL(dp), allocatable, dimension(:)  :: ExpI_aux, EypI_aux, EzpI_aux    
 INTEGER  :: i, j, k, n, unid, notnun(3)

 CHARACTER(len=50)		    :: nomevar

n = 8*(nelebloc( Nbloc+1)-1)

nomevar = 'vector of the real part of the field Ex 1D adjoint'
call  alloc_rdvet( ExpR_aux, n, nomevar )
nomevar = 'vector of the real part of the field Ey 1D adjoint'
call  alloc_rdvet( EypR_aux, n, nomevar )
nomevar = 'vector of the real part of the field Ez 1D adjoint'
call  alloc_rdvet( EzpR_aux, n, nomevar )
nomevar = 'vector of the imaginary part of the field Ex 1D adjoint'
call  alloc_rdvet( ExpI_aux, n, nomevar )
nomevar = 'vector of the imaginary part of the field Ey 1D adjoint'
call  alloc_rdvet( EypI_aux, n, nomevar )
nomevar = 'vector of the imaginary part of the field Ez 1D adjoint'
call  alloc_rdvet( EzpI_aux, n, nomevar )

IF( .NOT. ALLOCATED( Exp_adj ) ) THEN

	nomevar = 'Field Vector Ex 1D Adjoint'
	call  alloc_cdvet( Exp_adj, n, nomevar )
	nomevar = 'Field Vector Ey 1D Adjoint'
	call  alloc_cdvet( Eyp_adj, n, nomevar )
	nomevar = 'Field Vector Ez 1D Adjoint'
	call  alloc_cdvet( Ezp_adj, n, nomevar )
END IF


    if( FDP == 'ExDP' )then
	unid = 3010
elseif( FDP == 'EyDP' )then
	unid = 3110
elseif( FDP == 'EzDP' )then
	unid = 3210
elseif( FDP == 'HxDP' )then
	unid = 3310
elseif( FDP == 'HyDP' )then
	unid = 3410
elseif( FDP == 'HzDP' )then
	unid = 3510
end if

Exp_adj = (0.d0,0.d0)
Eyp_adj = (0.d0,0.d0)
Ezp_adj = (0.d0,0.d0)

ExpR_aux = 0.d0
EypR_aux = 0.d0
EzpR_aux = 0.d0

ExpI_aux = 0.d0
EypI_aux = 0.d0
EzpI_aux = 0.d0

notnun  = 0

read(unid) ExpR_aux
read(unid) ExpI_aux
read(unid) EypR_aux
read(unid) EypI_aux
read(unid) EzpR_aux
read(unid) EzpI_aux 

do j=nelebloc(1), n
			
	if( isnan( ExpR_aux(j) ) .or. isnan( ExpI_aux(j) ) )then
		notnun(1) = notnun(1) + 1
	end if
	if( isnan( EypR_aux(j) ) .or. isnan( EypI_aux(j) ) )then
		notnun(2) = notnun(2) + 1
	end if
	if( isnan( EzpR_aux(j) ) .or. isnan( EzpI_aux(j) ) )then
		notnun(3) = notnun(3) + 1
	end if
		
end do

	if( notnun(1) /= 0 )then
		print*, "Number of NaN in the evaluation of the receiver's primary Ex:"
		print*, notnun(1)
		print*, ' '
	end if
	if( notnun(2) /= 0 )then
		print*, "Number of NaN in the evaluation of the receiver's primary Ey:"
		print*, notnun(2)
		print*, ' '
	end if
	if( notnun(3) /= 0 )then
		print*, "Number of NaN in the evaluation of the receiver's primary Ez:"
		print*, notnun(3)
		print*, ' '
	end if

	if( notnun(1)/= 0 .or. notnun(2)/= 0 .or. notnun(3)/= 0)stop

	Exp_adj = ExpR_aux + ci*ExpI_aux
	Eyp_adj = EypR_aux + ci*EypI_aux
	Ezp_adj = EzpR_aux + ci*EzpI_aux

deallocate( ExpR_aux, EypR_aux, EzpR_aux )
deallocate( ExpI_aux, EypI_aux, EzpI_aux )

END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE w_x_gaussleg_QWE
IMPLICIT NONE
 INTEGER:: i, n
 CHARACTER(len=20)  :: doc1
 CHARACTER(len=50)  :: nomevar

	doc1  = 'xw_gaussleg.dat' 
	open( unit=800, file=doc1, DEFAULTFILE=trim(pathr)//'1D_font_orig/', status='old', action='read' )

	read(800,*) n

	nomevar = 'vector of abscissa for quadrature of Gauss-Legendree'
	call  alloc_rdvet( xGssLg, n, nomevar )


	nomevar = 'vector of Gauss-Legendre quadrature weights'
	call  alloc_rdvet( wGssLg, n, nomevar )

	do i=1, n
		read(800,*) xGssLg(i), wGssLg(i)
	end do

	close(800)

END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE Zeros_J0J1
IMPLICIT NONE
 INTEGER		   :: i, nzrs_J0, nzrs_J1
 CHARACTER(len=20) :: doc1
 CHARACTER(len=50) :: nomevar

	doc1  = 'zeros_J0.dat' 
	open( unit=900, file=doc1, DEFAULTFILE=trim(pathr)//'1D_font_orig/', status='old', action='read' )


	doc1  = 'zeros_J1.dat' 
	open( unit=950, file=doc1, DEFAULTFILE=trim(pathr)//'1D_font_orig/', status='old', action='read' )

	read(900,*) nzrs_J0 ;  read(950,*) nzrs_J1

	nomevar = 'zeros vector of Bessel Function J0'
	call  alloc_rdvet( zrs_J0, nzrs_J0, nomevar )

	nomevar = 'zeros vector of Bessel Function J1'
	call  alloc_rdvet( zrs_J1, nzrs_J1, nomevar )

	do i=1, nzrs_J0
		read(900,*) zrs_J0(i)
	end do
	do i=1, nzrs_J1
		read(950,*) zrs_J1(i)
	end do

	close(900) ; CLOSE(950)

END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE zizf_kmadaj

IMPLICIT NONE
INTEGER:: l
	
	IF( NC > 1 )THEN
		allocate(zkmadaj(0:NC-1))
		zkmadaj(0) = 0.d0
		DO  l=1, NC-1
			zkmadaj(l) = sum(hj(1:l))
		END DO
	END IF

END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE alloc_intmat( mtrz, diml, dimc, NOMMAT )
 INTEGER, ALLOCATABLE, DIMENSION(:,:), INTENT(out) :: mtrz
 INTEGER, 			       INTENT(in)  :: diml, dimc
 CHARACTER(LEN=50),		       INTENT(in)  :: NOMMAT 

 INTEGER :: IERR

	ALLOCATE( mtrz(diml, dimc), STAT = IERR )
	IF ( IERR == 0 )THEN
		mtrz =0.d0
	ELSE IF ( IERR /= 0 ) THEN
		WRITE(*,*)' ERROR AT ALLOCATION ', TRIM(NOMMAT)
		STOP
	END IF

END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE alloc_rdmat( mtrz, diml, dimc, NOMMAT )
 REAL(dp), ALLOCATABLE, DIMENSION(:,:), INTENT(out) :: mtrz
 INTEGER, 			        INTENT(in)  :: diml, dimc
 CHARACTER(LEN=50),		        INTENT(in)  :: NOMMAT 

 INTEGER :: IERR

	ALLOCATE( mtrz(diml, dimc), STAT = IERR )
	IF ( IERR == 0 )THEN
		mtrz = 0.d0
	ELSE IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR AT ALLOCATION ', TRIM(NOMMAT)
		STOP
	END IF

END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE alloc_cdmat( mtrz, diml, dimc, NOMMAT )
 COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:), INTENT(out) :: mtrz
 INTEGER, 			      	   INTENT(in)  :: diml, dimc
 CHARACTER(LEN=50),		       	   INTENT(in)  :: NOMMAT 

 INTEGER :: IERR

	ALLOCATE( mtrz(diml, dimc), STAT = IERR )
	IF ( IERR == 0 )THEN
		mtrz = ( 0.d0, 0.d0 )
	ELSE IF ( IERR /= 0 ) THEN
		WRITE(*,*)' ERROR AT ALLOCATION ', TRIM(NOMMAT)
		STOP
	END IF

END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE alloc_intVet( vet, dimv, NOMVET  )

 INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(out) :: vet
 INTEGER,			     INTENT(in)  :: dimv
 CHARACTER(LEN=50),		     INTENT(in)  :: NOMVET 

 INTEGER ::  IERR

	ALLOCATE( vet(dimv), STAT = IERR )
	IF ( IERR == 0 )THEN
		vet = 0
	ELSE IF ( IERR /= 0 ) THEN
		WRITE(*,*)' ERROR AT ALLOCATION ', TRIM(NOMVET)
		STOP
	END IF

END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================


SUBROUTINE alloc_rdvet( vet, dimv, NOMVET )
 REAL(dp), ALLOCATABLE, DIMENSION(:), INTENT(out) :: vet
 INTEGER ,			      INTENT(in)  :: dimv
 CHARACTER(LEN=50),		      INTENT(in)  :: NOMVET 

 INTEGER ::  IERR

	ALLOCATE( vet(dimv), STAT = IERR )
	IF ( IERR == 0 )THEN
		vet = 0.d0
	ELSE IF ( IERR /= 0 ) THEN
		WRITE(*,*)' ERROR AT ALLOCATION ', TRIM(NOMVET)
		STOP
	END IF

END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE alloc_cdvet( vet, dimv, NOMVET )
COMPLEX(dp), ALLOCATABLE, DIMENSION(:), INTENT(out) :: vet
 INTEGER ,			        INTENT(in)  :: dimv
 CHARACTER(LEN=50),		        INTENT(in)  :: NOMVET 

 INTEGER ::  IERR

	ALLOCATE( vet(dimv), STAT = IERR )
	IF ( IERR == 0 )THEN
		vet = ( 0.d0, 0d0 )
	ELSE IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR AT ALLOCATION ', TRIM(NOMVET)
		STOP
	END IF


END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE print_dados_1D2D3D_font( nnosmalha, Nelem, narestas )

implicit none
integer, intent(in)	:: nnosmalha, Nelem, narestas 

integer:: i
print*,'============================================================================'
Print*,'				Source Data:'
print*,'============================================================================'
Print*,"' Source's frequencies: "
print*,'	'
print*, vetf
print*,'	'
Print*,' Electric current from the loop or dipole(s): '
print*,'	'
print*, ISOURCE
print*,'	'
Print*,' Radius the loop or length of the dipole(s): '
print*,'	'
print*, FDIMENS
print*,'	'
print*,'============================================================================'
Print*,'				1D model data:'
print*,'============================================================================'
print*,' Number of layers: ', NC
print*,' Layer resistivities 1,2, ... NC: '
	do i=1, NC
		print*,'					', roj(i)
	end do
		print*,'	'
print*,' Thickness of layers 1,2, ... NC:'
	do i=1, NC
		print*,'					', hj(i)
	end do
print*,'============================================================================'
print*,'============================================================================'
Print*,'				2D model data:'
print*,'============================================================================'
print*,' Number of global nodes	    :', nnosmalha
print*,' Number of Elements	    :', Nelem
print*,' Number of edges   	    :', narestas
print*,' Number of heterogeneities  :', NCorpos
print*,'============================================================================'

END SUBROUTINE print_dados_1D2D3D_font

! ================================================================================================================
! 8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
! ================================================================================================================

END MODULE ler_escrever_arq
