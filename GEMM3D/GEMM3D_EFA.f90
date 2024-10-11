MODULE GEMM3D

USE G_variables
USE ler_escrever_arq
USE mdlgm1_1D
USE mdlgm2_1D
USE DADOS_SAIDA

CONTAINS

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================

SUBROUTINE inic_entrada( path_raiz )
 CHARACTER(LEN=300), INTENT(in) :: path_raiz

       pathr = path_raiz

	call inici_entradas
	call coluna_ptr_3D( Narestas )

END SUBROUTINE

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================

SUBROUTINE ATUALIZANDO_vrho_e( rhohet )
IMPLICIT NONE
 REAL(dp), DIMENSION(:), INTENT(in) :: rhohet

 INTEGER :: i, j

do i=1, Nbloc
	do j=nelebloc(i), nelebloc(i+1)-1 
		vrho_e(elebloc(j)) = rhohet(i)
	end do
end do

rojc = rhohet

END SUBROUTINE

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================

SUBROUTINE MODELO_SINTET_MT
 IMPLICIT NONE

 INTEGER  :: i
 
 print*,'======================================================================='
 print*,'STARTING 3D DATA SIMULATION'
 print*,'	'

 call openRW_arq_MT
 do i=1, Nfreq
	
	call MODELO_DIRETO_MT( 1, vetf(i) )

	call escrever_resultados_MT3D( Es_obs_1, Hs_obs_1, Es_obs_2, Hs_obs_2, Et_obs_1, Ht_obs_1, &
					Et_obs_2, Ht_obs_2, rho_ap_xx, rho_ap_yy, rho_ap_xy, &
					rho_ap_yx, fase_xx, fase_yy, fase_xy, fase_yx, Nobs )

	call MODELO_DIRETO_MT( 2, vetf(i) )

 end do
 call close_arq_MT

 print*,'3D DATA SIMULATION FINISHED'
 print*,'======================================================================='

END SUBROUTINE

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================

SUBROUTINE MODELO_DIRETO_MT( etapa, valf )
 IMPLICIT NONE
 REAL(dp), INTENT(in) :: valf
 INTEGER , INTENT(in) :: etapa

 INTEGER  :: modo(2)
 

 if( etapa == 1 )then

	 modo(1) = 1
	 modo(2) = 2

	call INICIANDO_PARDISO

	call montagem_matriz_EFA( valf )

	call MT_3D_EFA( valf )

	call liberando_memory_PARDISO( ES_1 )

	call campo_H_E_obs_MT( valf, modo(1), Es_1, Es_obs_1, Hs_obs_1, Et_obs_1, Ht_obs_1 )

	call campo_H_E_obs_MT( valf, modo(2), Es_2, Es_obs_2, Hs_obs_2, Et_obs_2, Ht_obs_2 )

	call rho_ap_fase(valf)

	nnz = nnz_aux

 elseif( etapa == 2 )then

	deallocate( Es_1, Es_2, val_nz, row_ind, Es_obs_1, Hs_obs_1, Es_obs_2, Hs_obs_2, &
			   Et_obs_1, Ht_obs_1, Et_obs_2, Ht_obs_2, rho_ap_xx, rho_ap_yy, &
	  		   rho_ap_xy, rho_ap_yx, fase_xx, fase_yy, fase_xy, fase_yx )

 end if

END SUBROUTINE

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================


SUBROUTINE MODELO_SINTET_MT_CSIG
 IMPLICIT NONE

 INTEGER  :: i
 
print*,'============================================================================'
 print*,'STARTING DATA SIMULATION 3D MT WITH SIGMA (w)'
 print*,'	'

 call openRW_arq_MT
 do i=1, Nfreq
	
	call MODELO_DIRETO_MT_CSIG( 1, vetf(i) )

	call escrever_resultados_MT3D( Es_obs_1, Hs_obs_1, Es_obs_2, Hs_obs_2, Et_obs_1, Ht_obs_1, &
					Et_obs_2, Ht_obs_2, rho_ap_xx, rho_ap_yy, rho_ap_xy, &
					rho_ap_yx, fase_xx, fase_yy, fase_xy, fase_yx, Nobs )

	call MODELO_DIRETO_MT_CSIG( 2, vetf(i) )

 end do
 call close_arq_MT

 print*,'3D MT DATA SIMULATION WITH SIGMA (w) FINISHED'
print*,'============================================================================'

END SUBROUTINE

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================

SUBROUTINE MODELO_DIRETO_MT_CSIG( etapa, valf )
 IMPLICIT NONE
 REAL(dp), INTENT(in) :: valf
 INTEGER , INTENT(in) :: etapa

 INTEGER  :: modo(2)
 

 if( etapa == 1 )then

	 modo(1) = 1
	 modo(2) = 2

	call INICIANDO_PARDISO

	call montagem_matriz_EFA_CSIG( valf )

	call MT_3D_EFA_CSIG( valf )

	call liberando_memory_PARDISO( ES_1 )

	call campo_H_E_obs_MT( valf, modo(1), Es_1, Es_obs_1, Hs_obs_1, Et_obs_1, Ht_obs_1 )

	call campo_H_E_obs_MT( valf, modo(2), Es_2, Es_obs_2, Hs_obs_2, Et_obs_2, Ht_obs_2 )

	call rho_ap_fase(valf)

	nnz = nnz_aux

 elseif( etapa == 2 )then

	deallocate( Es_1, Es_2, val_nz, row_ind, Es_obs_1, Hs_obs_1, Es_obs_2, Hs_obs_2, &
			   Et_obs_1, Ht_obs_1, Et_obs_2, Ht_obs_2, rho_ap_xx, rho_ap_yy, &
	  		   rho_ap_xy, rho_ap_yx, fase_xx, fase_yy, fase_xy, fase_yx )

 end if

END SUBROUTINE

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================

SUBROUTINE MODELO_SINTET_DPEH( FDP, Trans, nproc )
 IMPLICIT NONE
 CHARACTER(LEN=20), INTENT(in) :: FDP
 INTEGER	  , INTENT(in) :: Trans, nproc

 INTEGER  :: i, l


 if( Trans == 1 )then
	print*,'============================================================================'
	print*,'	'
	print*,'CALCULATING 1D FIELD FOR THE DIPOLE TRANSMITTER: ',trim(FDP)
	call DPEH1D( FDP, nproc )
	print*,'CALCULATION OF DIPOLE 1D FIELD TRANSMITTER FINISHED'		
	print*,'	'
	print*,'============================================================================'
 else							
	print*,'============================================================================'
	call openRW_arq_DPEM1D( FDP )
 endif

 print*,'============================================================================'
 print*,'INITIATING 3D SIMULATION OF THE DIPOLE FIELDS: ',trim(FDP)
 print*,'	'

 call openRW_arq_DPEH( FDP )
 do i=1, Nfreq

	call INICIANDO_PARDISO

	call montagem_matriz_EFA( vetf(i) )

	do l=1, Ntransm	
	
		call READ_primarios_Transm( FDP )
	
		call DPEH_3D_EFA_Transm( FDP, vetf(i), i, l, Es_1 )
		
		call campo_H_E_obs_DPEH( vetf(i), cxT(l), cyT(l), czT(l), FDP, Es_1, Es_obs_DP, Hs_obs_DP, Et_obs_DP, Ht_obs_DP )

		call escrever_resultados_DPEH( Es_obs_DP, Hs_obs_DP, Et_obs_DP, Ht_obs_DP, Nobs )
	
		deallocate( Es_obs_DP, Hs_obs_DP, Et_obs_DP, Ht_obs_DP  )	

	end do
	
	call liberando_memory_PARDISO( Es_1 )

	deallocate( val_nz, row_ind, Es_1 )

	nnz = nnz_aux

 end do
 call close_arq_DPEH

 call close_arq_DPEM1D( FDP )

 print*,'3D SIMULATION OF THE DIPOLE FIELDS:',trim(FDP),'FINISHED'
print*,'============================================================================'


END SUBROUTINE

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================

SUBROUTINE MODELO_SINTET_DPEH_CSIG( FDP, Trans, nproc )
 IMPLICIT NONE
 CHARACTER(LEN=20), INTENT(in) :: FDP
 INTEGER	  , INTENT(in) :: Trans, nproc

 INTEGER  :: i, l


 if( Trans == 1 )then
	print*,'============================================================================'
	print*,'	'
	print*,'CALCULATING 1D FIELD FOR THE DIPOLE TRANSMITTER: ',trim(FDP),'WITH SIGMA(w)'
	call DPEH1D( FDP, nproc )
	print*,'CALCULATION OF DIPOLE 1D FIELD TRANSMITTER FINISHED'		
	print*,'	'
	print*,'============================================================================'
 else							
	print*,'============================================================================'
	call openRW_arq_DPEM1D( FDP )
 endif

 print*,'============================================================================'
 print*,'INITIATING 3D SIMULATION OF THE DIPOLE FIELDS: ',trim(FDP),'WITH SIGMA(w)'
 print*,'	'

 call openRW_arq_DPEH( FDP )
 do i=1, Nfreq

	call INICIANDO_PARDISO

	call montagem_matriz_EFA_CSIG( vetf(i) )

	do l=1, Ntransm	
	
		call READ_primarios_Transm( FDP )
	
		call DPEH_3D_EFA_Transm_CSIG( FDP, vetf(i), i, l, Es_1 )

		
		call campo_H_E_obs_DPEH( vetf(i), cxT(l), cyT(l), czT(l), FDP, Es_1, Es_obs_DP, Hs_obs_DP, Et_obs_DP, Ht_obs_DP )

		call escrever_resultados_DPEH( Es_obs_DP, Hs_obs_DP, Et_obs_DP, Ht_obs_DP, Nobs )
	
		deallocate( Es_obs_DP, Hs_obs_DP, Et_obs_DP, Ht_obs_DP  )	

	end do
	
	call liberando_memory_PARDISO( Es_1 )

	deallocate( val_nz, row_ind, Es_1 )

	nnz = nnz_aux

 end do
 call close_arq_DPEH

 call close_arq_DPEM1D( FDP )

 print*,'3D SIMULATION OF THE DIPOLE FIELDS: ',trim(FDP),'FINISHED'
print*,'============================================================================'


END SUBROUTINE

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================


SUBROUTINE MODELO_SINTET_CSAMT( armz1D, nproc )
 IMPLICIT NONE
 INTEGER, INTENT(in) :: armz1D, nproc


 CHARACTER(LEN=20) :: FDP1, FDP2
 INTEGER	   :: i, l

 
 FDP1 = 'ExDP'
 FDP2 = 'EyDP'
 if( armz1D == 1 )then
	print*,'============================================================================'
	print*,'	'
	print*,'CALCULATING 1D FIELD FOR THE DIPOLES TRANSMITTER: ', trim(FDP1),' ',trim(FDP2)
	call DPEH1D( FDP1, nproc )
	call DPEH1D( FDP2, nproc )
	print*,'CALCULATION OF DIPOLES 1D FIELD TRANSMITTER FINISHED'		
	print*,'	'
	print*,'============================================================================'
 else							
	print*,'============================================================================'
	call openRW_arq_DPEM1D( FDP1 )
	call openRW_arq_DPEM1D( FDP2 )
 endif
 
 print*,'============================================================================'
 print*,'INITIATING 3D SIMULATION OF THE CSAMT FIELDS'
 print*,'	'


 call openRW_arq_CSAMT
 do i=1, Nfreq

	call INICIANDO_PARDISO

	call montagem_matriz_EFA( vetf(i) )

	do l=1, Ntransm	
	
		call READ_primarios_Transm( FDP1 )
	
		call DPEH_3D_EFA_Transm( FDP1, vetf(i), i, l, Es_1 )
		
		call campo_H_E_obs_DPEH( vetf(i), cxT(l), cyT(l), czT(l), FDP1, Es_1, Es_obs_1, Hs_obs_1, Et_obs_1, Ht_obs_1 )

		call READ_primarios_Transm( FDP2 )
	
		call DPEH_3D_EFA_Transm( FDP2, vetf(i), i, 2, Es_2 )
		
		call campo_H_E_obs_DPEH( vetf(i), cxT(l), cyT(l), czT(l), FDP2, Es_2, Es_obs_2, Hs_obs_2, Et_obs_2, Ht_obs_2 )

                call rho_ap_fase(vetf(i))

		call escrever_resultados_CSAMT3D( Es_obs_1, Hs_obs_1, Es_obs_2, Hs_obs_2, Et_obs_1, Ht_obs_1, &
						  Et_obs_2, Ht_obs_2, rho_ap_xx, rho_ap_yy, rho_ap_xy, &
						  rho_ap_yx, fase_xx, fase_yy, fase_xy, fase_yx, Nobs )

		deallocate( Es_obs_1, Hs_obs_1, Es_obs_2, Hs_obs_2, Et_obs_1, Ht_obs_1, Et_obs_2, Ht_obs_2 )
                deallocate( rho_ap_xx, rho_ap_yy, rho_ap_xy, rho_ap_yx, fase_xx, fase_yy, fase_xy, fase_yx )

	end do
	
	call liberando_memory_PARDISO( Es_1 )

	deallocate( val_nz, row_ind, Es_1, Es_2 )

	nnz = nnz_aux

 end do
 call close_arq_CSAMT
 call close_arq_DPEM1D( FDP1 )
 call close_arq_DPEM1D( FDP2 )

 print*,'CSAMT 3D DATA SIMULATION FINISHED'
 print*,'============================================================================'

END SUBROUTINE

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================

SUBROUTINE MODELO_SINTET_CSAMT_CSIG( armz1D, nproc )
 IMPLICIT NONE
 INTEGER, INTENT(in) :: armz1D, nproc


 CHARACTER(LEN=20) :: FDP1, FDP2
 INTEGER	   :: i, l

 
 FDP1 = 'ExDP'
 FDP2 = 'EyDP'
 if( armz1D == 1 )then
	print*,'============================================================================'
	print*,'	'
	print*,'CALCULATING 1D FIELD FOR THE DIPOLES TRANSMITTER: ', trim(FDP1),trim(FDP2)
	call DPEH1D( FDP1, nproc )
	call DPEH1D( FDP2, nproc )
	print*,'CALCULATION OF DIPOLES 1D FIELD TRANSMITTER FINISHED'		
	print*,'	'
	print*,'============================================================================'
 else							
	print*,'============================================================================'
	call openRW_arq_DPEM1D( FDP1 )
	call openRW_arq_DPEM1D( FDP2 )
 endif
 
 print*,'======================================================================='
 print*,'STARTING 3D CSAMT DATA SIMULATION WITH COMPLEX CONDUCTIVITY'
 print*,'	'


 call openRW_arq_CSAMT
 do i=1, Nfreq

	call INICIANDO_PARDISO

	call montagem_matriz_EFA_CSIG( vetf(i) )

	do l=1, Ntransm	
	
		call READ_primarios_Transm( FDP1 )
	
		call DPEH_3D_EFA_Transm_CSIG( FDP1, vetf(i), i, l, Es_1 )
		
		call campo_H_E_obs_DPEH( vetf(i), cxT(l), cyT(l), czT(l), FDP1, Es_1, Es_obs_1, Hs_obs_1, Et_obs_1, Ht_obs_1 )

		call READ_primarios_Transm( FDP2 )
	
		call DPEH_3D_EFA_Transm_CSIG( FDP2, vetf(i), i, 2, Es_2 )
		
		call campo_H_E_obs_DPEH( vetf(i), cxT(l), cyT(l), czT(l), FDP2, Es_2, Es_obs_2, Hs_obs_2, Et_obs_2, Ht_obs_2 )

                call rho_ap_fase(vetf(i))

		call escrever_resultados_MT3D( Es_obs_1, Hs_obs_1, Es_obs_2, Hs_obs_2, Et_obs_1, Ht_obs_1, &
						Et_obs_2, Ht_obs_2, rho_ap_xx, rho_ap_yy, rho_ap_xy, &
						rho_ap_yx, fase_xx, fase_yy, fase_xy, fase_yx, Nobs )

		deallocate( Es_obs_1, Hs_obs_1, Es_obs_2, Hs_obs_2, Et_obs_1, Ht_obs_1, Et_obs_2, Ht_obs_2 )
                deallocate( rho_ap_xx, rho_ap_yy, rho_ap_xy, rho_ap_yx, fase_xx, fase_yy, fase_xy, fase_yx )
	end do
	
	call liberando_memory_PARDISO( Es_1 )

	deallocate( val_nz, row_ind, Es_1, Es_2 )

	nnz = nnz_aux

 end do
 call close_arq_CSAMT
 call close_arq_DPEM1D( FDP1 )
 call close_arq_DPEM1D( FDP2 )

 print*,'STARTING 3D CSAMT DATA SIMULATION WITH COMPLEX CONDUCTIVITY FINISHED'
 print*,'======================================================================='

END SUBROUTINE

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================

SUBROUTINE DIRETO_CSAMT( i, valf, FDP1, FDP2 )
 IMPLICIT NONE
 CHARACTER(len=20), INTENT(in) :: FDP1, FDP2
 REAL(dp)         , INTENT(in) :: valf
 INTEGER          , INTENT(in) :: i 

 INTEGER :: l

  do l=1, Ntransm	

	call READ_primarios_Transm( FDP1 )

	call DPEH_3D_EFA_Transm( FDP1, valf, i, l, Es_1 )

	call campo_H_E_obs_DPEH( valf, cxT(l), cyT(l), czT(l), FDP1, Es_1, Es_obs_1, Hs_obs_1, Et_obs_1, Ht_obs_1 )

	call READ_primarios_Transm( FDP2 )

	call DPEH_3D_EFA_Transm( FDP2, valf, i, 2, Es_2 )
		
	call campo_H_E_obs_DPEH( valf, cxT(l), cyT(l), czT(l), FDP2, Es_2, Es_obs_2, Hs_obs_2, Et_obs_2, Ht_obs_2 )

        call rho_ap_fase(vetf(i))

	deallocate( Es_obs_1, Hs_obs_1, Es_obs_2, Hs_obs_2, Et_obs_1, Ht_obs_1, Et_obs_2, Ht_obs_2 )

 end do

 END SUBROUTINE 


!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================

SUBROUTINE MODELO_SINTET_HCL( Trans, nproc )
 IMPLICIT NONE
 INTEGER, INTENT(in) :: Trans, nproc

 INTEGER 	   :: i, l
 CHARACTER(len=20) :: FDP
 REAL(dp)	   :: t1, t2
 
 call w_x_gaussleg_QWE		! abscissas e pesos Para a Integração com QWE
 call Zeros_J0J1		! Zeros reais J0 e J1 p/ definir os subintervalos de aplicação da quadratura
 FDP = 'HCL'

 if( Trans == 1 )then
	print*,'============================================================================'
	print*,'	'
	print*,'CALCULATING 1D FIELD FOR THE TRANSMITTER: ',trim(FDP)

	call BCH1D( nproc )

	print*,'CALCULATING 1D FIELD FOR THE TRANSMITTER FINISHED'		
	print*,'	'
	print*,'============================================================================'
 else							
	print*,'============================================================================'
	call openRW_arq_DPEM1D( FDP )
 endif
 print*,'======================================================================='
 print*,'STARTING 3D SIMULATION OF HCL FIELDS'
 print*,'	'

 call openRW_arq_BCH

 do i=1, Nfreq

	call INICIANDO_PARDISO

	call montagem_matriz_EFA( vetf(i) )

	do l=1, Ntransm

		call READ_primarios_Transm_BCH

		call BCH_3D_EFA( vetf(i), i, l, Es_1, cxT(l), cyT(l), czT(l)) 
	
		call campo_H_E_obs_DPEH( vetf(i), cxT(l), cyT(l), czT(l), FDP, Es_1, Es_obs_BCH, Hs_obs_BCH, Et_obs_BCH, Ht_obs_BCH )
	
		call escrever_resultados_BCH3D( Es_obs_BCH, Hs_obs_BCH, Et_obs_BCH, Ht_obs_BCH, Nobs )

		deallocate( Es_obs_BCH, Hs_obs_BCH, Et_obs_BCH, Ht_obs_BCH  )

	end do

	nnz = nnz_aux

	call liberando_memory_PARDISO( Es_1 )
	deallocate( val_nz, row_ind, Es_1 )

 end do

 call close_arq_DPEM1D( FDP )
 call close_arq_BCH
 print*,'3D SIMULATION OF THE HCL FIELDS FINISHED'
 print*,'======================================================================='


END SUBROUTINE

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================

SUBROUTINE MODELO_SINTET_HCL_CSIG( Trans, nproc )
 IMPLICIT NONE
 INTEGER, INTENT(in) :: Trans, nproc

 INTEGER  :: i, l, K, readoff(2)
 CHARACTER(len=20) :: FDP
	
 REAL(dp) :: t1, t2
 
 call w_x_gaussleg_QWE		! abscissas e pesos Para a Integração com QWE
 call Zeros_J0J1		! Zeros reais J0 e J1 p/ definir os subintervalos de aplicação da quadratura
 FDP = 'HCL'

 if( Trans == 1 )then
	print*,'============================================================================'
	print*,'	'
	print*,'CALCULATING 1D FIELD FOR THE TRANSMITTER: ',trim(FDP)
	call BCH1D( nproc )
	print*,'CALCULATING 1D FIELD FOR THE TRANSMITTER FINISHED'
	print*,'	'
	print*,'============================================================================'
 else							
	print*,'============================================================================'
	call openRW_arq_DPEM1D( FDP )
 endif
 print*,'======================================================================='
 print*,'STARTING 3D SIMULATION OF HCL FIELDS WITH SIGMA(w)'
 print*,'	'

 call openRW_arq_BCH

 do i=1, Nfreq

	call INICIANDO_PARDISO

	call montagem_matriz_EFA_CSIG( vetf(i) )

	do l=1, Ntransm

		call READ_primarios_Transm_BCH

		call BCH_3D_EFA_CSIG( vetf(i), i, l, Es_1, cxT(l), cyT(l), czT(l)) 
	
		call campo_H_E_obs_DPEH( vetf(i), cxT(l), cyT(l), czT(l), FDP, Es_1, Es_obs_BCH, Hs_obs_BCH, Et_obs_BCH, Ht_obs_BCH )
	
		call escrever_resultados_BCH3D( Es_obs_BCH, Hs_obs_BCH, Et_obs_BCH, Ht_obs_BCH, Nobs )

		deallocate( Es_obs_BCH, Hs_obs_BCH, Et_obs_BCH, Ht_obs_BCH  )

	end do

	nnz = nnz_aux
	call liberando_memory_PARDISO( Es_1 )
	deallocate( val_nz, row_ind, Es_1 )

 end do

 call close_arq_DPEM1D( FDP )
 call close_arq_BCH
 print*,'3D SIMULATION OF THE HCL FIELDS FINISHED'
 print*,'======================================================================='


END SUBROUTINE

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================

SUBROUTINE JACOB_DPEM3D( DPT, DPR, fadj, Trans, nproc, MSB, Yestm )
 IMPLICIT NONE
 REAL(dp), DIMENSION(:,:), INTENT(inout) :: MSB
 REAL(dp), DIMENSION(:)  , INTENT(inout) :: Yestm
 INTEGER		 , INTENT(in)    :: fadj, Trans, nproc
 CHARACTER(LEN=20)	 , INTENT(in)    :: DPT, DPR

 INTEGER				:: Licamp, Lfcamp, Liphi, Lfphi, nL
 INTEGER				:: i, k, j, l, etapa, nL2


 if( fadj == 1 )then
	print*,'==========================================================='
	print*,'	'
	print*,'CALCULATING THE 1D FIELD FOR THE ADJOINT DIPOLAR SOURCE:', trim(DPR)
	print*,'	'
	call DPEH1D_adj( DPR, nproc )
	print*,'CALCULATING OF THE 1D FIELD FOR THE ADJOINT SOURCE FINISHED'
	print*,'	'
	print*,'==========================================================='
 end if

 if( Trans == 1 )then

	IF( DPT /= 'HCL' )THEN
		print*,'==========================================================='
		print*,'	'
		print*,'CALCULATING 1D FIELD FOR THE TRANSMITTER:',trim(DPT)
		print*,'	'
		call DPEH1D( DPT, nproc )
		print*,'CALCULATING 1D FIELD FOR THE TRANSMITTER FINISHED'		
		print*,'	'
		print*,'==========================================================='
	ELSE

		call w_x_gaussleg_QWE		! abscissas e pesos Para a Integração com QWE
		call Zeros_J0J1			! Zeros reais J0 e J1 p/ definir os subintervalos de aplicação da quadratura

		print*,'==========================================================='
		print*,'	'
		print*,'CALCULATING 1D FIELD FOR THE TRANSMITTER:',trim(DPT)
		print*,'	'
		call BCH1D( nproc )
		print*,'CALCULATING 1D FIELD FOR THE TRANSMITTER FINISHED'		
		print*,'	'
		print*,'==========================================================='

	END IF

 else							
	print*,'==========================================================='
	call openRW_arq_DPEM1D( DPT )
	if( DPT == 'HCL' )then
		call w_x_gaussleg_QWE		! abscissas e pesos Para a Integração com QWE
		call Zeros_J0J1			! Zeros reais J0 e J1 p/ definir os subintervalos de aplicação da quadratura
	 end if

 endif


 print*,'	'
 Print*,'STARTING THE CALCULATING OF THE JACOBIAN - DIPOLE'
 print*,'	'
 nL = Nobs*Nfreq*Ntransm
 MSB = 0.d0
 Yestm = 0.d0
 do i=1, Nfreq

	etapa = 1
	call SENSIB_DPEM3D( DPT, DPR, etapa, cxT(1), cyT(1), czT(1), 1, vetf(i), &
				  i, MSB(Licamp:Lfcamp,:), MSB(Liphi:Lfphi,:) )

	nL2 = (i-1)*Ntransm*Nobs
	do l=1, Ntransm

		Licamp = nL2 + (l-1)*Nobs+1
		Lfcamp = nL2 + l*Nobs
		Liphi  = nL+ nL2 + (l-1)*Nobs+1
		Lfphi  = nL + nL2 + l*Nobs

		etapa = 2
		call SENSIB_DPEM3D( DPT, DPR, etapa, cxT(l), cyT(l), czT(l), l+1, vetf(i), &
					  i, MSB(Licamp:Lfcamp,:), MSB(Liphi:Lfphi,:) )

		etapa = 3
		call SENSIB_DPEM3D( DPT, DPR, etapa, cxT(l), cyT(l), czT(l), l+1, vetf(i), &
					  i, MSB(Licamp:Lfcamp,:), MSB(Liphi:Lfphi,:) )

		  if( DPR == 'ExDP' )then

			Yestm(Licamp:Lfcamp) = log(cdabs(Et_obs_1(:,1)))
			Yestm( Liphi:Liphi ) = (180.d0/pi)*datan2( aimag(Et_obs_1(:,1)), &
								   dreal(Et_obs_1(:,1)) )
		elseif( DPR == 'EyDP' )then

			Yestm(Licamp:Lfcamp) = log(cdabs(Et_obs_1(:,2))) 
			Yestm( Liphi:Liphi ) = (180.d0/pi)*datan2( aimag(Et_obs_1(:,2)), &
								   dreal(Et_obs_1(:,2)) )
		elseif( DPR == 'EzDP' )then

			Yestm(Licamp:Lfcamp) = log(cdabs(Et_obs_1(:,3))) 
			Yestm( Liphi:Liphi ) = (180.d0/pi)*datan2( aimag(Et_obs_1(:,3)), &
								   dreal(Et_obs_1(:,3)) )

		elseif( DPR == 'HxDP' )then

			Yestm(Licamp:Lfcamp) = log(cdabs(Ht_obs_1(:,1)))
			Yestm( Liphi:Liphi ) = (180.d0/pi)*datan2( aimag(Ht_obs_1(:,1)), &
								   dreal(Ht_obs_1(:,1)) )
		elseif( DPR == 'HyDP' )then

			Yestm(Licamp:Lfcamp) = log(cdabs(Ht_obs_1(:,2))) 
			Yestm( Liphi:Liphi ) = (180.d0/pi)*datan2( aimag(Ht_obs_1(:,2)), &
								   dreal(Ht_obs_1(:,2)) )
		elseif( DPR == 'HzDP' )then

			Yestm(Licamp:Lfcamp) = log(cdabs(Ht_obs_1(:,3))) 
			Yestm( Liphi:Liphi ) = (180.d0/pi)*datan2( aimag(Ht_obs_1(:,3)), &
								   dreal(Ht_obs_1(:,3)) )
		end if
	end do	
	etapa = 4
	call SENSIB_DPEM3D( DPT, DPR, etapa, cxT(1), cyT(1), czT(1), 0, vetf(i), &
				  i, MSB(Licamp:Lfcamp,:), MSB(Liphi:Lfphi,:) )

 end do
 call close_arq_DPEM1D( DPT )

 Print*,'CALCULATING OF THE JACOBIAN FINISHED - DIPOLE'
 print*,'	'
 print*,'==========================================================='

 END SUBROUTINE

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================

SUBROUTINE SENSIB_DPEM3D( DPT, DPR, etapa, Tx, Ty, Tz, indt, f, indf, MSBcamp, MSBfase )
	
 IMPLICIT NONE
 REAL(dp), DIMENSION(:,:), INTENT(out) :: MSBcamp, MSBfase
 REAL(dp),		   INTENT(in)  :: f, Tx, Ty, Tz
 INTEGER ,		   INTENT(in)  :: etapa, indf, indt
 CHARACTER(LEN=20),	   INTENT(in)  :: DPT, DPR

 INTEGER :: IERR

 if( etapa == 1 )then

	call INICIANDO_PARDISO

	call montagem_matriz_EFA( f )

 else if( etapa == 2 )then

	IF( DPT/='HCL')THEN

		call READ_primarios_Transm( DPT )
		call DPEH_3D_EFA_Transm( DPT, f, indf, indt, Es_1 )
       ELSE

	       call READ_primarios_Transm_BCH
               call BCH_3D_EFA( f, indf, indt, Es_1 )
       END IF


	call campo_H_E_obs_DPEH( f, Tx, Ty, Tz, DPT, Es_1, Es_obs_1, Hs_obs_1, Et_obs_1, Ht_obs_1 )

 else if( etapa == 3 )then

	call openRW_arq_DPEM1D_adj( DPR )

	call SENSIB_DPEM( DPT, DPR, f, indf, MSBcamp, MSBfase ) 

	call close_arq_DPEM1D_adj( DPR )

 else if( etapa == 4 )then

	call liberando_memory_PARDISO( Es_1 )

	nnz = nnz_aux

	deallocate(  val_nz, row_ind, Es_1, Es_obs_1, Hs_obs_1, Et_obs_1, Ht_obs_1 )	

 end if


END SUBROUTINE

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================

SUBROUTINE SENSIB_DPEM( DPT, DPR, f, indf, MSBcamp, MSBfase )
 IMPLICIT NONE
 REAL(dp), DIMENSION(:,:), INTENT(out) :: MSBcamp, MSBfase
 REAL(dp),	               INTENT(in)  :: f
 CHARACTER(LEN=20),	       INTENT(in)  :: DPT, DPR
 INTEGER, 	      	       INTENT(in)  :: indf

 COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:) :: MSF1 
 COMPLEX(dp) :: zeta
 INTEGER     :: i, j, IERR

	ALLOCATE( MSF1(Nobs,Nbloc), STAT = IERR )
		IF ( IERR == 0 )THEN
		MSF1 = (0.d0,0.d0)
	ELSE IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN ALLOCATING MATRIX (MSF1)'
		WRITE(*,*)'ON THE SUBROUTINE SENSIB_DPEM'
		STOP
	END IF

	ALLOCATE( Es_adj_1(Narestas*ngl), STAT = IERR )
		IF ( IERR == 0 )THEN
		Es_adj_1 = (0.d0,0.d0)
	ELSE IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN ALLOCATING VECTOR (Es_adj_1)'
		WRITE(*,*)'ON THE SUBROUTINE SENSIB_DPEM'
		STOP
	END IF

!	print*,'Iniciando Sensibilidade da Componente :',trim(DPR)

	if( indf > 1 )then
		do i=1, indf-1
			do j=1, Nobs
				call READ_primarios_ADJ( DPR )! Leitura  do primario : DPR
			end do
		end do
	end if

	do j=1, Nobs

		call READ_primarios_ADJ( DPR )! Leitura  do primario : DPR

		call DPEH_3D_EFA_Fadj( f, Es_adj_1 )! Calculando Sensibiidade da componente : DPR

		call SENSIB_DPR_3( DPT, DPR, f, Es_adj_1, j, MSF1 )!; print*,'Final Sensibilidade ',trim(DPR),' obs',j
	end do

	if( DPR == 'HxDP' .or. DPR == 'HyDP' .or. DPR == 'HzDP' )then
		zeta = (0.d0,1.d0)*2.d0*pi*f*mi0
		MSF1 = MSF1/(-zeta)
	end if

	call MSB_camp_fase( DPR, MSF1, MSBcamp, MSBfase )

	DEALLOCATE( MSF1 , STAT = IERR )
	IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN DEALLOCATING MATRIX (MSF1)'
		WRITE(*,*)'ON THE SUBROUTINE SENSIB_DPEM'
		STOP
	END IF

	DEALLOCATE( Es_adj_1 , STAT = IERR )
	IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN DEALLOCATING VECTOR (Es_adj_1)'
		WRITE(*,*)'ON THE SUBROUTINE SENSIB_DPEM'
		STOP
	END IF

END SUBROUTINE

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================

SUBROUTINE JACOB_DPEM3D_CSIG( DPT, DPR, fadj, Trans, nproc, MSB, Yestm )
 IMPLICIT NONE
 REAL(dp), DIMENSION(:,:), INTENT(inout) :: MSB
 REAL(dp), DIMENSION(:)  , INTENT(inout) :: Yestm
 INTEGER		 , INTENT(in)    :: fadj, Trans, nproc
 CHARACTER(LEN=20)	 , INTENT(in)    :: DPT, DPR

 INTEGER				:: Licamp, Lfcamp, Liphi, Lfphi, nL
 INTEGER				:: i, k, j, l, etapa, nL2


 if( fadj == 1 )then
	print*,'==========================================================='
	print*,'	'
	print*,'CALCULATING THE 1D FIELD FOR THE ADJOINT DIPOLAR SOURCE:', trim(DPR)
	print*,'	'
	call DPEH1D_adj( DPR, nproc )
	print*,'CALCULATING OF THE 1D FIELD FOR THE ADJOINT SOURCE FINISHED'
	print*,'	'
	print*,'==========================================================='
 end if

 if( Trans == 1 )then

	IF( DPT /= 'HCL' )THEN
		print*,'==========================================================='
		print*,'	'
		print*,'CALCULATING 1D FIELD FOR THE TRANSMITTER:',trim(DPT)
		print*,'	'
		call DPEH1D( DPT, nproc )
		print*,'CALCULATING 1D FIELD FOR THE TRANSMITTER FINISHED'		
		print*,'	'
		print*,'==========================================================='
	ELSE

		call w_x_gaussleg_QWE		! abscissas e pesos Para a Integração com QWE
		call Zeros_J0J1			! Zeros reais J0 e J1 p/ definir os subintervalos de aplicação da quadratura

		print*,'==========================================================='
		print*,'	'
		print*,'CALCULATING 1D FIELD FOR THE TRANSMITTER:',trim(DPT)
		print*,'	'
		call BCH1D( nproc )
		print*,'CALCULATING 1D FIELD FOR THE TRANSMITTER FINISHED'		
		print*,'	'
		print*,'==========================================================='

	END IF

 else							
	print*,'==========================================================='
	call openRW_arq_DPEM1D( DPT )
	if( DPT == 'HCL' )then
		call w_x_gaussleg_QWE		! abscissas e pesos Para a Integração com QWE
		call Zeros_J0J1			! Zeros reais J0 e J1 p/ definir os subintervalos de aplicação da quadratura
	 end if

 endif


 print*,'	'
 Print*,'STARTING THE CALCULATING OF THE JACOBIAN - DIPOLE WITH SIGMA(w)'
 print*,'	'
 nL = Nobs*Nfreq*Ntransm
 MSB = 0.d0
 Yestm = 0.d0
 do i=1, Nfreq

!	print*,'FREQUÊNCIA', i
	etapa = 1
	call SENSIB_DPEM3D_CSIG( DPT, DPR, etapa, cxT(1), cyT(1), czT(1), 1, vetf(i), &
				  i, MSB(Licamp:Lfcamp,:), MSB(Liphi:Lfphi,:) )

	nL2 = (i-1)*Ntransm*Nobs
	do l=1, Ntransm

!		print*,'TRANSMISSOR', l

		Licamp = nL2 + (l-1)*Nobs+1
		Lfcamp = nL2 + l*Nobs
		Liphi  = nL+ nL2 + (l-1)*Nobs+1
		Lfphi  = nL + nL2 + l*Nobs

		etapa = 2
		call SENSIB_DPEM3D_CSIG( DPT, DPR, etapa, cxT(l), cyT(l), czT(l), l+1, vetf(i), &
					  i, MSB(Licamp:Lfcamp,:), MSB(Liphi:Lfphi,:) )

		etapa = 3
		call SENSIB_DPEM3D_CSIG( DPT, DPR, etapa, cxT(l), cyT(l), czT(l), l+1, vetf(i), &
					  i, MSB(Licamp:Lfcamp,:), MSB(Liphi:Lfphi,:) )

		  if( DPR == 'ExDP' )then

			Yestm(Licamp:Lfcamp) = log(cdabs(Et_obs_1(:,1)))
			Yestm( Liphi:Liphi ) = (180.d0/pi)*datan2( aimag(Et_obs_1(:,1)), &
								   dreal(Et_obs_1(:,1)) )
		elseif( DPR == 'EyDP' )then

			Yestm(Licamp:Lfcamp) = log(cdabs(Et_obs_1(:,2))) 
			Yestm( Liphi:Liphi ) = (180.d0/pi)*datan2( aimag(Et_obs_1(:,2)), &
								   dreal(Et_obs_1(:,2)) )
		elseif( DPR == 'EzDP' )then

			Yestm(Licamp:Lfcamp) = log(cdabs(Et_obs_1(:,3))) 
			Yestm( Liphi:Liphi ) = (180.d0/pi)*datan2( aimag(Et_obs_1(:,3)), &
								   dreal(Et_obs_1(:,3)) )

		elseif( DPR == 'HxDP' )then

			Yestm(Licamp:Lfcamp) = log(cdabs(Ht_obs_1(:,1)))
			Yestm( Liphi:Liphi ) = (180.d0/pi)*datan2( aimag(Ht_obs_1(:,1)), &
								   dreal(Ht_obs_1(:,1)) )
		elseif( DPR == 'HyDP' )then

			Yestm(Licamp:Lfcamp) = log(cdabs(Ht_obs_1(:,2))) 
			Yestm( Liphi:Liphi ) = (180.d0/pi)*datan2( aimag(Ht_obs_1(:,2)), &
								   dreal(Ht_obs_1(:,2)) )
		elseif( DPR == 'HzDP' )then

			Yestm(Licamp:Lfcamp) = log(cdabs(Ht_obs_1(:,3))) 
			Yestm( Liphi:Liphi ) = (180.d0/pi)*datan2( aimag(Ht_obs_1(:,3)), &
								   dreal(Ht_obs_1(:,3)) )
		end if
	end do	
	etapa = 4
	call SENSIB_DPEM3D_CSIG( DPT, DPR, etapa, cxT(1), cyT(1), czT(1), 0, vetf(i), &
				  i, MSB(Licamp:Lfcamp,:), MSB(Liphi:Lfphi,:) )

 end do
 call close_arq_DPEM1D( DPT )

 Print*,'CALCULATING OF THE JACOBIAN FINISHED - DIPOLE '
 print*,'	'
 print*,'==========================================================='

 END SUBROUTINE

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================

SUBROUTINE SENSIB_DPEM3D_CSIG( DPT, DPR, etapa, Tx, Ty, Tz, indt, f, indf, MSBcamp, MSBfase )
	
 IMPLICIT NONE
 REAL(dp), DIMENSION(:,:), INTENT(out) :: MSBcamp, MSBfase
 REAL(dp),		   INTENT(in)  :: f, Tx, Ty, Tz
 INTEGER ,		   INTENT(in)  :: etapa, indf, indt
 CHARACTER(LEN=20),	   INTENT(in)  :: DPT, DPR

 INTEGER :: IERR

 if( etapa == 1 )then

	call INICIANDO_PARDISO

	call montagem_matriz_EFA_CSIG( f )


 else if( etapa == 2 )then

	IF( DPT/='HCL')THEN

		call READ_primarios_Transm( DPT )
		call DPEH_3D_EFA_Transm_CSIG( DPT, f, indf, indt, Es_1 )
       ELSE

	       call READ_primarios_Transm_BCH
               call BCH_3D_EFA_CSIG( f, indf, indt, Es_1 )
       END IF

	call campo_H_E_obs_DPEH( f, Tx, Ty, Tz, DPT, Es_1, Es_obs_1, Hs_obs_1, Et_obs_1, Ht_obs_1 )

 else if( etapa == 3 )then

	call openRW_arq_DPEM1D_adj( DPR )

	call SENSIB_DPEM_CSIG( DPT, DPR, f, indf, MSBcamp, MSBfase ) 

	call close_arq_DPEM1D_adj( DPR )

 else if( etapa == 4 )then

	call liberando_memory_PARDISO( Es_1 )

	nnz = nnz_aux

	deallocate(  val_nz, row_ind, Es_1, Es_obs_1, Hs_obs_1, Et_obs_1, Ht_obs_1 )	

 end if


END SUBROUTINE

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================

SUBROUTINE SENSIB_DPEM_CSIG( DPT, DPR, f, indf, MSBcamp, MSBfase )
 IMPLICIT NONE
 REAL(dp), DIMENSION(:,:), INTENT(out) :: MSBcamp, MSBfase
 REAL(dp),	               INTENT(in)  :: f
 CHARACTER(LEN=20),	       INTENT(in)  :: DPT, DPR
 INTEGER, 	      	       INTENT(in)  :: indf

 COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:) :: MSF1 
 COMPLEX(dp) :: zeta
 INTEGER     :: i, j, IERR

	ALLOCATE( MSF1(Nobs,Nbloc), STAT = IERR )
		IF ( IERR == 0 )THEN
		MSF1 = (0.d0,0.d0)
	ELSE IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN ALLOCATING MATRIX (MSF1)'
		WRITE(*,*)'ON THE SUBROUTINE SENSIB_DPEM_CSIG'
		STOP
	END IF

	ALLOCATE( Es_adj_1(Narestas*ngl), STAT = IERR )
		IF ( IERR == 0 )THEN
		Es_adj_1 = (0.d0,0.d0)
	ELSE IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN ALLOCATING VECTOR (Es_adj_1)'
		WRITE(*,*)'ON THE SUBROUTINE SENSIB_DPEM_CSIG'
		STOP
	END IF

!	print*,'Iniciando Sensibilidade da Componente :',trim(DPR)

	if( indf > 1 )then
		do i=1, indf-1
			do j=1, Nobs
				call READ_primarios_ADJ( DPR )! Leitura  do primario : DPR
			end do
		end do
	end if

	do j=1, Nobs

		call READ_primarios_ADJ( DPR )! Leitura  do primario : DPR

		call DPEH_3D_EFA_Fadj_CSIG( f, Es_adj_1 )! Calculando Sensibiidade da componente : DPR

		call SENSIB_DPR_3( DPT, DPR, f, Es_adj_1, j, MSF1 )!; print*,'Final Sensibilidade ',trim(DPR),' obs',j
	end do

	if( DPR == 'HxDP' .or. DPR == 'HyDP' .or. DPR == 'HzDP' )then
		zeta = (0.d0,1.d0)*2.d0*pi*f*mi0
		MSF1 = MSF1/(-zeta)
	end if

	call MSB_camp_fase( DPR, MSF1, MSBcamp, MSBfase )

	DEALLOCATE( MSF1 , STAT = IERR )
	IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN DEALLOCATING MATRIX (MSF1)'
		WRITE(*,*)'ON THE SUBROUTINE SENSIB_DPEM_CSIG'
		STOP
	END IF

	DEALLOCATE( Es_adj_1 , STAT = IERR )
	IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN DEALLOCATING VECTOR (Es_adj_1)'
		WRITE(*,*)'ON THE SUBROUTINE SENSIB_DPEM_CSIG'
		STOP
	END IF

END SUBROUTINE

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================

SUBROUTINE SENSIB_DPR_3( DPT, DPR, f, Es_adj, lin, MSF1 )

IMPLICIT NONE
 COMPLEX(dp), DIMENSION(:,:), INTENT(inout) :: MSF1
 COMPLEX(dp), DIMENSION(:)  , INTENT(in)    :: Es_adj
 REAL(dp)   ,		      INTENT(in)    :: f
 INTEGER    ,		      INTENT(in)    :: lin
 CHARACTER(LEN=20), 	      INTENT(in)    :: DPT, DPR

COMPLEX(dp), DIMENSION(12)	:: Ep_adj, E_adj, E_f1
COMPLEX(dp), DIMENSION(8)	:: campx_D, campy_D, campz_D, campx_A, campy_A, campz_A
REAL(dp)   , DIMENSION(12,12)	:: MF
REAL(dp)   , DIMENSION(4,4)	:: ml

COMPLEX(dp), DIMENSION(8)	:: Exp_D, Eyp_D, Ezp_D, Ex_D, Ey_D, Ez_D
COMPLEX(dp), DIMENSION(8)	:: Exp_A, Eyp_A, Ezp_A, Ex_A, Ey_A, Ez_A

COMPLEX(dp)	:: inte1, Ex, Ey, Ez, Hx, Hy, Hz, soma, neta0, zeta, inte(3)
REAL(dp)	:: dsx, dsy, dsz, Iw, w, sigmas(NC), l(3)
REAL(dp)	:: x(2), y(2), z(2), xrec, yrec, zrec, Rinte(3), Iinte(3)!, Rinte2(3), Iinte2(3)
INTEGER		:: cont, i, j, k, p, nosG_e(8), aresta(12), inde, ncoef, nloc


campx_D = (0.d0,0.d0)
campy_D = (0.d0,0.d0)
campz_D = (0.d0,0.d0)

campx_A = (0.d0,0.d0)
campy_A = (0.d0,0.d0)
campz_A = (0.d0,0.d0)

Ex_D = (0.d0,0.d0)
Ey_D = (0.d0,0.d0)
Ez_D = (0.d0,0.d0)

Ex_A = (0.d0,0.d0)
Ey_A = (0.d0,0.d0)
Ez_A = (0.d0,0.d0)

Iw  = 1.d0
dsx = 1.d0
dsy = 1.d0
dsz = 1.d0

xrec = Mnos(no_obs(lin),1)
yrec = Mnos(no_obs(lin),2)
zrec = Mnos(no_obs(lin),3)

w = (pi+pi)*f
sigmas = 1.d0/roj
neta0 = 1.d-12 
zeta = (0.d0,1.d0)*w*mi0

do k=1, Nbloc

	inte1 = (0.d0,0.d0)

	do j=nelebloc(k), nelebloc(k+1)-1

		inde	  = elebloc(j)
		nosG_e(:) = Melem(inde,:)
		aresta(:) = Marestas(inde,:)
      
		x(1:2) = Mnos(nosG_e(1:2),1)
		y(1) = Mnos(nosG_e(1),2)
		y(2) = Mnos(nosG_e(4),2)
		z(1) = Mnos(nosG_e(1),3)
		z(2) = Mnos(nosG_e(5),3)

		p = j+j+j+j + j+j+j+j

		Exp_D(:) = Exp_trans(p-7:p)
		Eyp_D(:) = Eyp_trans(p-7:p)
		Ezp_D(:) = Ezp_trans(p-7:p)

		Exp_A(:) = Exp_adj(p-7:p)
		Eyp_A(:) = Eyp_adj(p-7:p)
		Ezp_A(:) = Ezp_adj(p-7:p)

		do i=1, 4

			E_f1(i)   = Es_1(aresta(i))
			E_f1(i+4) = Es_1(aresta(i+4)) 
			E_f1(i+8) = Es_1(aresta(i+8))
			
			E_adj(i)   = Es_adj(aresta(i))
			E_adj(i+4) = Es_adj(aresta(i+4))
			E_adj(i+8) = Es_adj(aresta(i+8))

		end do


		do i=1, 2

			Ex_D(i+i-1:i+i)   = E_f1(i+i-1)
			Ex_D(3+i+i:4+i+i) = E_f1(i+i)

			Ex_A(i+i-1:i+i)   = E_adj(i+i-1)
			Ex_A(3+i+i:4+i+i) = E_adj(i+i)

		end do
	
		Ey_D(1) = E_f1(5)				; Ey_A(1) = E_adj(5)
		Ey_D(5) = Ey_D(1)    				; Ey_A(5) = Ey_A(1)
			
		Ey_D(2) = E_f1(7)				; Ey_A(2) = E_adj(7)
		Ey_D(6) = Ey_D(2)  				; Ey_A(6) = Ey_A(2)  
			
		Ey_D(4) = E_f1(6)				; Ey_A(4) = E_adj(6)
		Ey_D(8) = Ey_D(4)   				; Ey_A(8) = Ey_A(4) 
			
		Ey_D(3) = E_f1(8)				; Ey_A(3) = E_adj(8)
		Ey_D(7) = Ey_D(3)   				; Ey_A(7) = Ey_A(3) 
		
		
		Ez_D(1) = E_f1(9)				; Ez_A(1) = E_adj(9)
		Ez_D(4) = Ez_D(1)  				; Ez_A(4) = Ez_A(1)
			
		Ez_D(2) = E_f1(10)				; Ez_A(2) = E_adj(10)
		Ez_D(3) = Ez_D(2)  				; Ez_A(3) = Ez_A(2)
				
		Ez_D(5) = E_f1(11)				; Ez_A(5) = E_adj(11)
		Ez_D(8) = Ez_D(5)  				; Ez_A(8) = Ez_A(5)
 			
		Ez_D(6) = E_f1(12)				; Ez_A(6) = E_adj(12)
		Ez_D(7) = Ez_D(6)  				; Ez_A(7) = Ez_A(6)

		
		Ex_D(1:2) = Ex_D(1:2) + Exp_D(1:2)		; Ex_A(1:2) = Ex_A(1:2) + Exp_A(1:2)
		Ex_D(7:8) = Ex_D(7:8) + Exp_D(7:8)		; Ex_A(7:8) = Ex_A(7:8) + Exp_A(7:8)

		Ex_D(3) = Ex_D(3) + Exp_D(6)			; Ex_A(3) = Ex_A(3) + Exp_A(6)
		Ex_D(4) = Ex_D(4) + Exp_D(5)			; Ex_A(4) = Ex_A(4) + Exp_A(5)
		Ex_D(5) = Ex_D(5) + Exp_D(4)			; Ex_A(5) = Ex_A(5) + Exp_A(4)
		Ex_D(6) = Ex_D(6) + Exp_D(3)			; Ex_A(6) = Ex_A(6) + Exp_A(3)

		Ey_D(1:2) = Ey_D(1:2) + Eyp_D(1:2)		; Ey_A(1:2) = Ey_A(1:2) + Eyp_A(1:2)
		Ey_D(7:8) = Ey_D(7:8) + Eyp_D(7:8)		; Ey_A(7:8) = Ey_A(7:8) + Eyp_A(7:8)

		Ey_D(3) = Ey_D(3) + Eyp_D(6)			; Ey_A(3) = Ey_A(3) + Eyp_A(6)
		Ey_D(4) = Ey_D(4) + Eyp_D(5)			; Ey_A(4) = Ey_A(4) + Eyp_A(5)
		Ey_D(5) = Ey_D(5) + Eyp_D(4)			; Ey_A(5) = Ey_A(5) + Eyp_A(4)
		Ey_D(6) = Ey_D(6) + Eyp_D(3)			; Ey_A(6) = Ey_A(6) + Eyp_A(3)

		Ez_D(1:2) = Ez_D(1:2) + Ezp_D(1:2)		; Ez_A(1:2) = Ez_A(1:2) + Ezp_A(1:2)
		Ez_D(7:8) = Ez_D(7:8) + Ezp_D(7:8)		; Ez_A(7:8) = Ez_A(7:8) + Ezp_A(7:8)

		Ez_D(3) = Ez_D(3) + Ezp_D(6)			; Ez_A(3) = Ez_A(3) + Ezp_A(6)
		Ez_D(4) = Ez_D(4) + Ezp_D(5)			; Ez_A(4) = Ez_A(4) + Ezp_A(5)
		Ez_D(5) = Ez_D(5) + Ezp_D(4)			; Ez_A(5) = Ez_A(5) + Ezp_A(4)
		Ez_D(6) = Ez_D(6) + Ezp_D(3)			; Ez_A(6) = Ez_A(6) + Ezp_A(3)

		inte(1) = IntVolBlc( Ex_D, Ex_A, x, y, z )
		inte(2) = IntVolBlc( Ey_D, Ey_A, x, y, z )
		inte(3) = IntVolBlc( Ez_D, Ez_A, x, y, z )

		inte1 = inte1 + sum(inte(:)) 
	
	end do
	
	MSF1(lin,k) = inte1
end do

deallocate( Exp_adj, Eyp_adj, Ezp_adj )			

END SUBROUTINE

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================

SUBROUTINE MSB_camp_fase( DPR, MSE_1, MSBcamp, MSBfase )
IMPLICIT NONE
 REAL(dp)   , DIMENSION(:,:), INTENT(out) :: MSBcamp, MSBfase
 COMPLEX(dp), DIMENSION(:,:), INTENT(in)  :: MSE_1
 CHARACTER(LEN=20), 	      INTENT(in)  :: DPR

 COMPLEX(dp),ALLOCATABLE, DIMENSION(:) :: EHt
 COMPLEX(dp)	:: Ex_dsigma 
 INTEGER	:: i, j, IERR

 MSBcamp = 0.d0
 MSBfase = 0.d0

 ALLOCATE( EHt(Nobs), STAT = IERR )
	IF ( IERR == 0 )THEN
	EHt = (0.d0,0.d0)
 ELSE IF ( IERR /= 0 ) THEN
	WRITE(*,*)'ERROR IN DEALLOCATING VECTOR (EHt)'
	WRITE(*,*)'ON THE SUBROUTINE MSB_camp_fase'
	STOP
 END IF

     if( DPR == 'ExDP' )then
	EHt  = Et_obs_1(:,1)
 else if( DPR == 'HxDP' )then       
	EHt  = Ht_obs_1(:,1)
 else if( DPR == 'EyDP' )then
	EHt  = Et_obs_1(:,2)
 else if( DPR == 'HyDP')then
	EHt  = Ht_obs_1(:,2)
 else if( DPR == 'EzDP' )then
	EHt  = Et_obs_1(:,3)
 else if( DPR == 'HzDP')then
	EHt  = Ht_obs_1(:,3)
 end if

 do j=1, Nbloc
	do i=1, Nobs 

		Ex_dsigma    = MSE_1(i,j)/rojc(j)/EHt(i)
		MSBcamp(i,j) = REAL( Ex_dsigma )
		MSBfase(i,j) = AIMAG( Ex_dsigma )

	end do
 end do

 DEALLOCATE( EHt , STAT = IERR )
 IF ( IERR /= 0 ) THEN
	WRITE(*,*)'ERROR IN DEALLOCATING VECTOR (EHt)'
	WRITE(*,*)'ON THE SUBROUTINE MSB_camp_fase'
	STOP
 END IF


END SUBROUTINE

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================

SUBROUTINE JACOB_CSAMT3D( locadj, armz1D , fadj, nproc, MSB, Yestm )
 IMPLICIT NONE
 REAL(dp), DIMENSION(:,:), INTENT(inout) :: MSB
 REAL(dp), DIMENSION(:)	 , INTENT(inout) :: Yestm
 INTEGER		 , INTENT(in)    :: locadj, armz1D, fadj, nproc
 
 CHARACTER(LEN=20)			:: elem_Zij, FDP1, FDP2
 INTEGER				:: i, k, j, etapa
 INTEGER				:: Lirhoxy, Lfrhoxy, Lirhoyx, Lfrhoyx, nL, nL2 
 INTEGER				:: Liphixy, Lfphixy, Liphiyx, Lfphiyx

 
 FDP1 = 'ExDP'
 FDP2 = 'EyDP'

 if( armz1D == 1 )then
	print*,'============================================================================'
	print*,'	'
	print*,'CALCULATING 1D FIELD FOR TRANSMITTER DIPOLES: ', trim(FDP1), trim(FDP2)
	print*,'	'
	call DPEH1D( FDP1, nproc )
	call DPEH1D( FDP2, nproc )
	print*,'CALCULATING 1D FIELD FOR TRANSMITTER DIPOLES FINISHED'		
	print*,'	'
	print*,'============================================================================'
 else							
	print*,'============================================================================'
	call openRW_arq_DPEM1D( FDP1 )
	call openRW_arq_DPEM1D( FDP2 )
 endif

 if( fadj == 1 )then
	print*,'==========================================================='
	print*,'	'
	print*,'CALCULATING 1D FIELD FOR ADJOINT SOURCES - CSAMT:'
	print*,'	'
	call DPEH1D_ADJ_CSAMT( nproc )
	print*,'CALCULATING 1D FIELD FOR ADJOINT SOURCES FINISHED.'
	print*,	'	'
	print*,'==========================================================='
 else		
	print*,'==========================================================='
	if( locadj == 0 )then
		call openRW_arq_1D_adj_CSAMT	
	elseif( locadj == 1 )then
		call openRW_arq_1D_adj_MT	
	endif
	
 end if

 print*,'	'
 Print*,'STARTING THE JACOBIAN CALCULATION - CSAMT'
 print*,'	'
 nL  = Nobs*Nfreq
 nL2 = 2*Nobs*Nfreq
 MSB = 0.d0
 Yestm = 0.d0
 do i=1, Nfreq

	Lirhoxy =  (i-1)*Nobs+1
	Lfrhoxy =   i*Nobs
	Liphixy = nL2 + (i-1)*Nobs+1
	Lfphixy = nL2 + i*Nobs

	Lirhoyx = nL + (i-1)*Nobs+1
	Lfrhoyx = nL +  i*Nobs
	Liphiyx = nL + nL2 + (i-1)*Nobs+1
	Lfphiyx = nL + nL2 +  i*Nobs

	elem_Zij = 'xy'
	etapa = 1
	call SENSIB_CSAMT3D( elem_Zij , etapa, i, vetf(i), MSB(Lirhoxy:Lfrhoxy,:), MSB(Liphixy:Lfphixy,:) )
!	print*,'Etapa 1 OK!'
	
	etapa = 2
	call SENSIB_CSAMT3D( elem_Zij , etapa, i, vetf(i), MSB(Lirhoxy:Lfrhoxy,:), MSB(Liphixy:Lfphixy,:) )
!	print*,'Etapa 2 OK!'
	
	etapa = 3
	call SENSIB_CSAMT3D( elem_Zij , etapa, i, vetf(i), MSB(Lirhoxy:Lfrhoxy,:), MSB(Liphixy:Lfphixy,:) )
!	print*,'Etapa 3 xy OK!'

	etapa = 3
	elem_Zij ='yx'
	call SENSIB_CSAMT3D( elem_Zij , etapa, i, vetf(i), MSB(Lirhoyx:Lfrhoyx,:), MSB(Liphiyx:Lfphiyx,:) )
!	print*,'Etapa 3 yx OK!'

	Yestm( Lirhoxy:Lfrhoxy ) = rho_ap_xy
	Yestm( Lirhoyx:Lfrhoyx ) = rho_ap_yx
	Yestm( Liphixy:Lfphixy ) = fase_xy
	Yestm( Liphiyx:Lfphiyx ) = fase_yx

	etapa = 4
	call SENSIB_CSAMT3D( elem_Zij, etapa, i, vetf(i), MSB(Lirhoyx:Lfrhoyx,:), MSB(Liphiyx:Lfphiyx,:) )
!	print*,'Etapa 4 OK!'
	
	nnz = nnz_aux
 end do

 call close_arq_DPEM1D( FDP1 )
 call close_arq_DPEM1D( FDP2 )

if( locadj == 0 )then
  call close_arq_1D_adj_CSAMT
else
  call close_arq_1D_adj_MT
endif

 Print*,'JACOBIAN CALCULATION FINISHED - CSAMT'
 print*,'	'
 print*,'==========================================================='



 END SUBROUTINE

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================

SUBROUTINE SENSIB_CSAMT3D( Z_ij, etapa, indf, f, MSBrho, MSBfase )
	
 IMPLICIT NONE
 REAL(dp), DIMENSION(:,:), INTENT(out) :: MSBrho, MSBfase
 REAL(dp),		   INTENT(in)  :: f
 INTEGER ,		   INTENT(in)  :: etapa, indf
 CHARACTER(LEN=20),	   INTENT(in)  :: Z_ij

 CHARACTER(LEN=50) :: nomevar
 CHARACTER(LEN=20) :: FDP1, FDP2
 INTEGER :: n, IERR


 FDP1 = 'ExDP'
 FDP2 = 'EyDP'
 n = 8*(nelebloc( Nbloc+1)-1)

if( etapa == 1 )then

	call INICIANDO_PARDISO

	call montagem_matriz_EFA( f )

	call READ_primarios_Transm( FDP1 )

	call DPEH_3D_EFA_Transm( FDP1, f, indf, 1, Es_1 )
		
	call campo_H_E_obs_DPEH( f, cxT(1), cyT(1), czT(1), FDP1, Es_1, Es_obs_1, Hs_obs_1, Et_obs_1, Ht_obs_1 )

	nomevar = 'vector for the field Ex 1D of the first Transmitter for the CSAMT'
	call  alloc_cdvet( Exp_trans1, n, nomevar )
	nomevar = 'vector for the field Ey 1D of the first Transmitter for the CSAMT'
	call  alloc_cdvet( Eyp_trans1, n, nomevar )
	nomevar = 'vector for the field Ez 1D of the first Transmitter for the CSAMT'
	call  alloc_cdvet( Ezp_trans1, n, nomevar )

        Exp_trans1 = Exp_trans
        Eyp_trans1 = Eyp_trans
        Ezp_trans1 = Ezp_trans

	call READ_primarios_Transm( FDP2 )
	
	call DPEH_3D_EFA_Transm( FDP2, f, indf, 2, Es_2 )
		
	call campo_H_E_obs_DPEH( f, cxT(1), cyT(1), czT(1), FDP2, Es_2, Es_obs_2, Hs_obs_2, Et_obs_2, Ht_obs_2 )

	nomevar = 'vector for the field Ex 1D of the second Transmitter for the CSAMT'
	call  alloc_cdvet( Exp_trans2, n, nomevar )
	nomevar = 'vector for the field Ey 1D of the second Transmitter for the CSAMT'
	call  alloc_cdvet( Eyp_trans2, n, nomevar )
	nomevar = 'vector for the field Ez 1D of the second Transmitter for the CSAMT'
	call  alloc_cdvet( Ezp_trans2, n, nomevar )

        Exp_trans2 = Exp_trans
        Eyp_trans2 = Eyp_trans
        Ezp_trans2 = Ezp_trans

	DEALLOCATE( Exp_trans , STAT = IERR )
	IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN DEALLOCATING VECTOR (Exp_trans)'
		STOP
	END IF
	DEALLOCATE( Eyp_trans , STAT = IERR )
	IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN DEALLOCATING VECTOR (Eyp_trans)'
		STOP
	END IF
	DEALLOCATE( Ezp_trans , STAT = IERR )
	IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN DEALLOCATING VECTOR (Ezp_trans)'
		STOP
	END IF

	call rho_ap_fase( f )

else if( etapa == 2 )then

	MT_CSAMT = 'CSAMT'
	call SENSIB_H( f )

else if( etapa == 3 )then

	MT_CSAMT = 'CSAMT'
	call SENSIB_E( Z_ij, f, MSBrho, MSBfase )

else if( etapa == 4 )then

	call liberando_memory_PARDISO( Es_1 )

	call liberando_variavs
end if

END SUBROUTINE SENSIB_CSAMT3D

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================

SUBROUTINE JACOB_CSAMT3D_CSIG( locadj, armz1D , fadj, nproc, MSB, Yestm )
 IMPLICIT NONE
 REAL(dp), DIMENSION(:,:), INTENT(inout) :: MSB
 REAL(dp), DIMENSION(:)	 , INTENT(inout) :: Yestm
 INTEGER		 , INTENT(in)    :: locadj, armz1D, fadj, nproc
 
 CHARACTER(LEN=20)			:: elem_Zij, FDP1, FDP2
 INTEGER				:: i, k, j, etapa
 INTEGER				:: Lirhoxy, Lfrhoxy, Lirhoyx, Lfrhoyx, nL, nL2 
 INTEGER				:: Liphixy, Lfphixy, Liphiyx, Lfphiyx

 
 FDP1 = 'ExDP'
 FDP2 = 'EyDP'

 if( armz1D == 1 )then
	print*,'============================================================================'
	print*,'	'
	print*,'CALCULATING 1D FIELD FOR TRANSMITTER DIPOLES: ', trim(FDP1), trim(FDP2)
	call DPEH1D( FDP1, nproc )
	call DPEH1D( FDP2, nproc )
	print*,'CALCULATING 1D FIELD FOR TRANSMITTER DIPOLES FINISHED'		
	print*,'	'
	print*,'============================================================================'
 else							
	print*,'============================================================================'
	call openRW_arq_DPEM1D( FDP1 )
	call openRW_arq_DPEM1D( FDP2 )
 endif

 if( fadj == 1 )then
	print*,'==========================================================='
	print*,'	'
	print*,'CALCULATING 1D FIELD FOR ADJOINT SOURCES - CSAMT:'
	print*,'	'
	call DPEH1D_ADJ_CSAMT( nproc )
	print*,'CALCULATING 1D FIELD FOR ADJOINT SOURCES FINISHED.'
	print*,	'	'
	print*,'==========================================================='
 else		
	print*,'==========================================================='
	if( locadj == 0 )then
		call openRW_arq_1D_adj_CSAMT	
	elseif( locadj == 1 )then
		call openRW_arq_1D_adj_MT	
	endif
	
 end if

 print*,'	'
 Print*,'STARTING THE JACOBIAN CALCULATION WITH SIGMA(w) - CSAMT '
 print*,'	'
 nL  = Nobs*Nfreq
 nL2 = 2*Nobs*Nfreq
 MSB = 0.d0
 Yestm = 0.d0
 do i=1, Nfreq

	Lirhoxy =  (i-1)*Nobs+1
	Lfrhoxy =   i*Nobs
	Liphixy = nL2 + (i-1)*Nobs+1
	Lfphixy = nL2 + i*Nobs

	Lirhoyx = nL + (i-1)*Nobs+1
	Lfrhoyx = nL +  i*Nobs
	Liphiyx = nL + nL2 + (i-1)*Nobs+1
	Lfphiyx = nL + nL2 +  i*Nobs

	elem_Zij = 'xy'
	etapa = 1
	call SENSIB_CSAMT3D_CSIG( elem_Zij , etapa, i, vetf(i), MSB(Lirhoxy:Lfrhoxy,:), MSB(Liphixy:Lfphixy,:) )
!	print*,'Etapa 1 OK!'
	
	etapa = 2
	call SENSIB_CSAMT3D_CSIG( elem_Zij , etapa, i, vetf(i), MSB(Lirhoxy:Lfrhoxy,:), MSB(Liphixy:Lfphixy,:) )
!	print*,'Etapa 2 OK!'
	
	etapa = 3
	call SENSIB_CSAMT3D_CSIG( elem_Zij , etapa, i, vetf(i), MSB(Lirhoxy:Lfrhoxy,:), MSB(Liphixy:Lfphixy,:) )
!	print*,'Etapa 3 xy OK!'

	etapa = 3
	elem_Zij ='yx'
	call SENSIB_CSAMT3D_CSIG( elem_Zij , etapa, i, vetf(i), MSB(Lirhoyx:Lfrhoyx,:), MSB(Liphiyx:Lfphiyx,:) )
!	print*,'Etapa 3 yx OK!'

	Yestm( Lirhoxy:Lfrhoxy ) = rho_ap_xy
	Yestm( Lirhoyx:Lfrhoyx ) = rho_ap_yx
	Yestm( Liphixy:Lfphixy ) = fase_xy
	Yestm( Liphiyx:Lfphiyx ) = fase_yx

	etapa = 4
	call SENSIB_CSAMT3D( elem_Zij, etapa, i, vetf(i), MSB(Lirhoyx:Lfrhoyx,:), MSB(Liphiyx:Lfphiyx,:) )
!	print*,'Etapa 4 OK!'
	
	nnz = nnz_aux
 end do

 call close_arq_DPEM1D( FDP1 )
 call close_arq_DPEM1D( FDP2 )

if( locadj == 0 )then
  call close_arq_1D_adj_CSAMT
else
  call close_arq_1D_adj_MT
endif

 Print*,'JACOBIAN CALCULATION FINISHED - CSAMT'
 print*,'	'
 print*,'==========================================================='



 END SUBROUTINE

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================

SUBROUTINE SENSIB_CSAMT3D_CSIG( Z_ij, etapa, indf, f, MSBrho, MSBfase )
	
 IMPLICIT NONE
 REAL(dp), DIMENSION(:,:), INTENT(out) :: MSBrho, MSBfase
 REAL(dp),		   INTENT(in)  :: f
 INTEGER ,		   INTENT(in)  :: etapa, indf
 CHARACTER(LEN=20),	   INTENT(in)  :: Z_ij

 CHARACTER(LEN=50) :: nomevar
 CHARACTER(LEN=20) :: FDP1, FDP2
 INTEGER :: n, IERR


 FDP1 = 'ExDP'
 FDP2 = 'EyDP'
 n = 8*(nelebloc( Nbloc+1)-1)

if( etapa == 1 )then

	call INICIANDO_PARDISO

	call montagem_matriz_EFA_CSIG( f )

	call READ_primarios_Transm( FDP1 )

	call DPEH_3D_EFA_Transm_CSIG( FDP1, f, indf, 1, Es_1 )

	call campo_H_E_obs_DPEH( f, cxT(1), cyT(1), czT(1), FDP1, Es_1, Es_obs_1, Hs_obs_1, Et_obs_1, Ht_obs_1 )

	nomevar = 'vector for the field Ex 1D of the first Transmitter for the CSAMT'
	call  alloc_cdvet( Exp_trans1, n, nomevar )
	nomevar = 'vector for the field Ey 1D of the first Transmitter for the CSAMT'
	call  alloc_cdvet( Eyp_trans1, n, nomevar )
	nomevar = 'vector for the field Ez 1D of the first Transmitter for the CSAMT'
	call  alloc_cdvet( Ezp_trans1, n, nomevar )

        Exp_trans1 = Exp_trans
        Eyp_trans1 = Eyp_trans
        Ezp_trans1 = Ezp_trans

	call READ_primarios_Transm( FDP2 )
	
	call DPEH_3D_EFA_Transm_CSIG( FDP2, f, indf, 2, Es_2 )
		
	call campo_H_E_obs_DPEH( f, cxT(1), cyT(1), czT(1), FDP2, Es_2, Es_obs_2, Hs_obs_2, Et_obs_2, Ht_obs_2 )

	nomevar = 'vector for the field Ex 1D of the second Transmitter for the CSAMT'
	call  alloc_cdvet( Exp_trans2, n, nomevar )
	nomevar = 'vector for the field Ey 1D of the second Transmitter for the CSAMT'
	call  alloc_cdvet( Eyp_trans2, n, nomevar )
	nomevar = 'vector for the field Ez 1D of the second Transmitter for the CSAMT'
	call  alloc_cdvet( Ezp_trans2, n, nomevar )

        Exp_trans2 = Exp_trans
        Eyp_trans2 = Eyp_trans
        Ezp_trans2 = Ezp_trans

	DEALLOCATE( Exp_trans , STAT = IERR )
	IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN DEALLOCATING VECTOR (Exp_trans)'
		STOP
	END IF
	DEALLOCATE( Eyp_trans , STAT = IERR )
	IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN DEALLOCATING VECTOR (Eyp_trans)'
		STOP
	END IF
	DEALLOCATE( Ezp_trans , STAT = IERR )
	IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN DEALLOCATING VECTOR (Ezp_trans)'
		STOP
	END IF

	call rho_ap_fase( f )

else if( etapa == 2 )then

	MT_CSAMT = 'CSAMT'
	call SENSIB_H( f )

else if( etapa == 3 )then

	MT_CSAMT = 'CSAMT'
	call SENSIB_E( Z_ij, f, MSBrho, MSBfase )

else if( etapa == 4 )then

	call liberando_memory_PARDISO( Es_1 )

	call liberando_variavs
end if

END SUBROUTINE

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================


SUBROUTINE JACOB_MT3D( locadj, fadj, nproc, MSB, Yestm )
 IMPLICIT NONE
 REAL(dp), DIMENSION(:,:), INTENT(inout) :: MSB
 REAL(dp), DIMENSION(:)	 , INTENT(inout) :: Yestm
 INTEGER		 , INTENT(in)    :: locadj, fadj, nproc
 
 CHARACTER(LEN=20)			:: elem_Zij
 INTEGER				:: i, k, j, etapa
 INTEGER				:: Lirhoxy, Lfrhoxy, Lirhoyx, Lfrhoyx, nL, nL2 
 INTEGER				:: Liphixy, Lfphixy, Liphiyx, Lfphiyx

 if( fadj == 1 )then
	print*,'==========================================================='
	print*,'	'
	print*,'CALCULATING 1D FIELD FOR ADJOINT SOURCES - MT:'
	print*,'	'
	call DPEH1D_ADJ_MT( nproc )
	print*,'CALCULATING 1D FIELD FOR ADJOINT SOURCES FINISHED - MT'
	print*,	'	'
	print*,'==========================================================='
 else		
	print*,'==========================================================='
	if( locadj == 0 )then
		call openRW_arq_1D_adj_MT	
	else
		call openRW_arq_1D_adj_CSAMT	
	endif
 end if

 print*,'	'
 Print*,'STARTING THE CALCULATION JACOBIAN - MT'
 print*,'	'
 nL  = Nobs*Nfreq
 nL2 = 2*Nobs*Nfreq
 MSB = 0.d0
 Yestm = 0.d0
 do i=1, Nfreq

	Lirhoxy =  (i-1)*Nobs+1
	Lfrhoxy =   i*Nobs
	Liphixy = nL2 + (i-1)*Nobs+1
	Lfphixy = nL2 + i*Nobs

	Lirhoyx = nL + (i-1)*Nobs+1
	Lfrhoyx = nL +  i*Nobs
	Liphiyx = nL + nL2 + (i-1)*Nobs+1
	Lfphiyx = nL + nL2 +  i*Nobs

	elem_Zij = 'xy'
	etapa = 1
	call SENSIB_MT3D( elem_Zij , etapa, vetf(i), MSB(Lirhoxy:Lfrhoxy,:), MSB(Liphixy:Lfphixy,:) )
	!print*,'Etapa 1 OK!'
	
	etapa = 2
	call SENSIB_MT3D( elem_Zij , etapa, vetf(i), MSB(Lirhoxy:Lfrhoxy,:), MSB(Liphixy:Lfphixy,:) )
	!print*,'Etapa 2 OK!'
	
	etapa = 3
	call SENSIB_MT3D( elem_Zij , etapa, vetf(i), MSB(Lirhoxy:Lfrhoxy,:), MSB(Liphixy:Lfphixy,:) )
	!print*,'Etapa 3 OK!'

	etapa = 3
	elem_Zij ='yx'
	call SENSIB_MT3D( elem_Zij , etapa, vetf(i), MSB(Lirhoyx:Lfrhoyx,:), MSB(Liphiyx:Lfphiyx,:) )
	!print*,'Etapa 3 OK!'

	Yestm( Lirhoxy:Lfrhoxy ) = rho_ap_xy
	Yestm( Lirhoyx:Lfrhoyx ) = rho_ap_yx
	Yestm( Liphixy:Lfphixy ) = fase_xy
	Yestm( Liphiyx:Lfphiyx ) = fase_yx

	etapa = 4
	call SENSIB_MT3D( elem_Zij, etapa, vetf(i), MSB(Lirhoyx:Lfrhoyx,:), MSB(Liphiyx:Lfphiyx,:) )
	!print*,'Etapa 4 OK!'
	
	nnz = nnz_aux
 end do

if( locadj == 0 )then
  call close_arq_1D_adj_MT
else
  call close_arq_1D_adj_CSAMT
endif


 Print*,'CALCULATION JACOBIAN FINISHED - MT'
 print*,'	'
 print*,'==========================================================='



 END SUBROUTINE

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================

SUBROUTINE SENSIB_MT3D( Z_ij, etapa, f, MSBrho, MSBfase )

 IMPLICIT NONE
 REAL(dp), DIMENSION(:,:), INTENT(out) :: MSBrho, MSBfase
 REAL(dp),		   INTENT(in)  :: f
 INTEGER ,		   INTENT(in)  :: etapa
 CHARACTER(LEN=20),	   INTENT(in)  :: Z_ij

 INTEGER :: modo, IERR

if( etapa == 1 )then

	call INICIANDO_PARDISO

	call montagem_matriz_EFA( f )

	call MT_3D_EFA( f )

	modo = 1
	call campo_H_E_obs_MT( f, modo, Es_1, Es_obs_1, Hs_obs_1, Et_obs_1, Ht_obs_1 )

	modo = 2
	call campo_H_E_obs_MT( f, modo, Es_2, Es_obs_2, Hs_obs_2, Et_obs_2, Ht_obs_2 )

	call rho_ap_fase( f )

else if( etapa == 2 )then

	MT_CSAMT = 'MT'
	call SENSIB_H( f )

else if( etapa == 3 )then

	MT_CSAMT = 'MT'
	call SENSIB_E( Z_ij, f, MSBrho, MSBfase )

else if( etapa == 4 )then

	call liberando_memory_PARDISO( Es_1 )

	call liberando_variavs
end if

END SUBROUTINE SENSIB_MT3D

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================

SUBROUTINE JACOB_MT3D_CSIG( locadj, fadj, nproc, MSB, Yestm )
 IMPLICIT NONE
 REAL(dp), DIMENSION(:,:), INTENT(inout) :: MSB
 REAL(dp), DIMENSION(:)	 , INTENT(inout) :: Yestm
 INTEGER		 , INTENT(in)    :: locadj, fadj, nproc
 
 CHARACTER(LEN=20)			:: elem_Zij
 INTEGER				:: i, k, j, etapa
 INTEGER				:: Lirhoxy, Lfrhoxy, Lirhoyx, Lfrhoyx, nL, nL2 
 INTEGER				:: Liphixy, Lfphixy, Liphiyx, Lfphiyx


 if( fadj == 1 )then
	print*,'==========================================================='
	print*,'	'
	print*,'CALCULATING 1D FIELD FOR ADJOINT SOURCES - MT:'
	print*,'	'
	call DPEH1D_ADJ_MT( nproc )
	print*,'CALCULATING 1D FIELD FOR ADJOINT SOURCES FINISHED - MT'
	print*,	'	'
	print*,'==========================================================='
 else		
	print*,'==========================================================='
	if( locadj == 0 )then
		call openRW_arq_1D_adj_MT	
	else
		call openRW_arq_1D_adj_CSAMT	
	endif
 end if

 print*,'	'
 Print*,'STARTING THE CALCULTION JACOBIAN WITH SIGMA(w) - MT'
 print*,'	'
 nL  = Nobs*Nfreq
 nL2 = 2*Nobs*Nfreq
 MSB = 0.d0
 Yestm = 0.d0
 do i=1, Nfreq

	Lirhoxy =  (i-1)*Nobs+1
	Lfrhoxy =   i*Nobs
	Liphixy = nL2 + (i-1)*Nobs+1
	Lfphixy = nL2 + i*Nobs

	Lirhoyx = nL + (i-1)*Nobs+1
	Lfrhoyx = nL +  i*Nobs
	Liphiyx = nL + nL2 + (i-1)*Nobs+1
	Lfphiyx = nL + nL2 +  i*Nobs

	elem_Zij = 'xy'
	etapa = 1
	call SENSIB_MT3D_CSIG( elem_Zij , etapa, vetf(i), MSB(Lirhoxy:Lfrhoxy,:), MSB(Liphixy:Lfphixy,:) )
	!print*,'Etapa 1 OK!'
	
	etapa = 2
	call SENSIB_MT3D_CSIG( elem_Zij , etapa, vetf(i), MSB(Lirhoxy:Lfrhoxy,:), MSB(Liphixy:Lfphixy,:) )
	!print*,'Etapa 2 OK!'
	
	etapa = 3
	call SENSIB_MT3D_CSIG( elem_Zij , etapa, vetf(i), MSB(Lirhoxy:Lfrhoxy,:), MSB(Liphixy:Lfphixy,:) )
	!print*,'Etapa 3 OK!'

	etapa = 3
	elem_Zij ='yx'
	call SENSIB_MT3D_CSIG( elem_Zij , etapa, vetf(i), MSB(Lirhoyx:Lfrhoyx,:), MSB(Liphiyx:Lfphiyx,:) )
	!print*,'Etapa 3 OK!'

	Yestm( Lirhoxy:Lfrhoxy ) = rho_ap_xy
	Yestm( Lirhoyx:Lfrhoyx ) = rho_ap_yx
	Yestm( Liphixy:Lfphixy ) = fase_xy
	Yestm( Liphiyx:Lfphiyx ) = fase_yx

	etapa = 4
	call SENSIB_MT3D_CSIG( elem_Zij, etapa, vetf(i), MSB(Lirhoyx:Lfrhoyx,:), MSB(Liphiyx:Lfphiyx,:) )
	!print*,'Etapa 4 OK!'
	
	nnz = nnz_aux
 end do

if( locadj == 0 )then
  call close_arq_1D_adj_MT
else
  call close_arq_1D_adj_CSAMT
endif


 Print*,'CALCULTION JACOBIAN FINISHED - MT'
 print*,'	'
 print*,'==========================================================='



 END SUBROUTINE

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================

SUBROUTINE SENSIB_MT3D_CSIG( Z_ij, etapa, f, MSBrho, MSBfase )

 IMPLICIT NONE
 REAL(dp), DIMENSION(:,:), INTENT(out) :: MSBrho, MSBfase
 REAL(dp),		   INTENT(in)  :: f
 INTEGER ,		   INTENT(in)  :: etapa
 CHARACTER(LEN=20),	   INTENT(in)  :: Z_ij

 INTEGER :: modo, IERR

if( etapa == 1 )then

	call INICIANDO_PARDISO

	call montagem_matriz_EFA_CSIG( f )

	call MT_3D_EFA_CSIG( f )

	modo = 1
	call campo_H_E_obs_MT( f, modo, Es_1, Es_obs_1, Hs_obs_1, Et_obs_1, Ht_obs_1 )

	modo = 2
	call campo_H_E_obs_MT( f, modo, Es_2, Es_obs_2, Hs_obs_2, Et_obs_2, Ht_obs_2 )

	call rho_ap_fase( f )

else if( etapa == 2 )then

	MT_CSAMT = 'MT'
	call SENSIB_H( f )

else if( etapa == 3 )then

	MT_CSAMT = 'MT'
	call SENSIB_E( Z_ij, f, MSBrho, MSBfase )

else if( etapa == 4 )then

	call liberando_memory_PARDISO( Es_1 )

	call liberando_variavs
end if

END SUBROUTINE

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================

SUBROUTINE liberando_variavs
INTEGER :: IERR

	DEALLOCATE( val_nz , STAT = IERR )
	IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN DEALLOCATING VECTOR (val_nz)'
		STOP
	END IF

	DEALLOCATE( row_ind , STAT = IERR )
	IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN DEALLOCATING VECTOR (row_ind)'
		STOP
	END IF

	DEALLOCATE( rho_ap_xx , STAT = IERR )
	IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN DEALLOCATING VECTOR (rho_ap_xx)'
		STOP
	END IF

	DEALLOCATE( rho_ap_yy , STAT = IERR )
	IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN DEALLOCATING VECTOR (rho_ap_yy)'
		STOP
	END IF
	DEALLOCATE( rho_ap_xy , STAT = IERR )
	IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN DEALLOCATING VECTOR (rho_ap_xy)'
		STOP
	END IF
	DEALLOCATE( rho_ap_yx , STAT = IERR )
	IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN DEALLOCATING VECTOR (rho_ap_yx)'
		STOP
	END IF

	DEALLOCATE( fase_xx , STAT = IERR )
	IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN DEALLOCATING VECTOR (fase_xx)'
		STOP
	END IF

	DEALLOCATE( fase_yy , STAT = IERR )
	IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN DEALLOCATING VECTOR (fase_yy)'
		STOP
	END IF
	DEALLOCATE( fase_xy , STAT = IERR )
	IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN DEALLOCATING VECTOR (fase_xy)'
		STOP
	END IF
	DEALLOCATE( fase_yx , STAT = IERR )
	IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN DEALLOCATING VECTOR (fase_yx)'
		STOP
	END IF

	DEALLOCATE( Es_1 , STAT = IERR )
	IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN DEALLOCATING VECTOR (Es_1)'
		STOP
	END IF

	DEALLOCATE( Es_2 , STAT = IERR )
	IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN DEALLOCATING VECTOR (Es_2)'
		STOP
	END IF

	DEALLOCATE( Es_obs_1 , STAT = IERR )
	IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN DEALLOCATING VECTOR (Es_obs_1)'
		STOP
	END IF

	DEALLOCATE( Es_obs_2 , STAT = IERR )
	IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN DEALLOCATING VECTOR (Es_obs_2)'
		STOP
	END IF

	DEALLOCATE( Et_obs_1 , STAT = IERR )
	IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN DEALLOCATING VECTOR (Et_obs_1)'
		STOP
	END IF

	DEALLOCATE( Et_obs_2 , STAT = IERR )
	IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN DEALLOCATING VECTOR (Et_obs_2)'
		STOP
	END IF

	DEALLOCATE( Hs_obs_1 , STAT = IERR )
	IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN DEALLOCATING VECTOR (Hs_obs_1)'
		STOP
	END IF

	DEALLOCATE( Hs_obs_2 , STAT = IERR )
	IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN DEALLOCATING VECTOR (Hs_obs_2)'
		STOP
	END IF

	DEALLOCATE( Ht_obs_1 , STAT = IERR )
	IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN DEALLOCATING VECTOR (Ht_obs_1)'
		STOP
	END IF

	DEALLOCATE( Ht_obs_2 , STAT = IERR )
	IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN DEALLOCATING VECTOR (Ht_obs_2)'
		STOP
	END IF

	DEALLOCATE( MSHx1 , STAT = IERR )
	IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN DEALLOCATING MATRIX (MSHx1)'
		STOP
	END IF

	DEALLOCATE( MSHx2 , STAT = IERR )
	IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN DEALLOCATING MATRIX (MSHx2)'
		STOP
	END IF

	DEALLOCATE( MSHy1 , STAT = IERR )
	IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN DEALLOCATING MATRIX (MSHx3)'
		STOP
	END IF

	DEALLOCATE( MSHy2 , STAT = IERR )
	IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN DEALLOCATING MATRIX (MSHx4)'
		STOP
	END IF

END SUBROUTINE

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================

SUBROUTINE SENSIB_H( f )
 IMPLICIT NONE
 REAL(dp) , INTENT(in)	:: f
 
 COMPLEX(dp) 	   :: zeta
 CHARACTER(LEN=20) :: FDP
 INTEGER	   :: i, j, IERR

	IF( .NOT. ALLOCATED( MSHx1 ) ) THEN
		allocate( MSHx1(Nobs,Nbloc), MSHx2(Nobs,Nbloc), MSHy1(Nobs,Nbloc), MSHy2(Nobs,Nbloc) )
	
		MSHx1 = (0.d0,0.d0)
		MSHx2 = (0.d0,0.d0)
		MSHy1 = (0.d0,0.d0)
		MSHy2 = (0.d0,0.d0)	
	ELSE	
		
		MSHx1 = (0.d0,0.d0)
		MSHx2 = (0.d0,0.d0)
		MSHy1 = (0.d0,0.d0)
		MSHy2 = (0.d0,0.d0)
	
	END IF	
	
		ALLOCATE( Es_adj_3(Narestas*ngl), STAT = IERR )
			IF ( IERR == 0 )THEN
			Es_adj_3 = (0.d0,0.d0)
		ELSE IF	( IERR /= 0 ) THEN
			WRITE(*,*)'ERROR IN ALLOCAR THE VECTOR (Es_adj_3)'
			WRITE(*,*)'ON THE SUBROUTINE SENSIB_H'
			STOP
		END IF
		ALLOCATE( Es_adj_4(Narestas*ngl), STAT = IERR )
			IF ( IERR == 0 )THEN
			Es_adj_4 = (0.d0,0.d0)
		ELSE IF ( IERR /= 0 ) THEN
			WRITE(*,*)'ERROR IN DEALLOCATING VECTOR (Es_adj_4)'
			WRITE(*,*)'ON THE SUBROUTINE SENSIB_H'
			STOP
		END IF

!	print*,'Iniciando sensibilidade das componente Hx e Hy'

	IF( MT_CSAMT == 'MT')THEN	

		do j=1, Nobs

!			print*,'Iniciando sensibilidade da componente Hx(j,:)' , j
			FDP = 'HxDP'
			call READ_primarios_ADJ_MT( fDP ); ! ' Leitura do primario DPHx '
		
			call DPEH_3D_EFA_Fadj( f, Es_adj_3 )! ' SOLUÇÂO do DPHx '
		
			call SENSIB_E_H_1( FDP, f, Es_adj_3, j, MSHx1, MSHx2 )!; print*,' Sensibilidade Hx(j,:)-(fim),j=',j 
		
!			print*,'Iniciando sensibilidade da componente Hy(j,:)' 
			FDP = 'HyDP'
			call READ_primarios_ADJ_MT( FDP )! ' Leitura do primario DPHy '
		
			call DPEH_3D_EFA_Fadj( f , Es_adj_4 )! ' SOLUÇÂO do DPHy '
	
			call SENSIB_E_H_1( FDP, f, Es_adj_4, j, MSHy1, MSHy2 ); !print*,' Sensibilidade Hy(j,:)-(fim),j=',j
	
		end do

	ELSEIF(	MT_CSAMT == 'CSAMT' )THEN

		do j=1, Nobs

!			print*,'Iniciando sensibilidade da componente Hx(j,:)' , j
			FDP = 'HxDP'
			call READ_primarios_ADJ_MT( FDP ); !print*, ' Leitura do primario DPHx '
		
			call DPEH_3D_EFA_Fadj( f, Es_adj_3 ); !print*, ' SOLUÇÂO do DPHx '
		
			call SENSIB_E_H_CSAMT( FDP, f, Es_adj_3, j, MSHx1, MSHx2 ); !print*,' Sensibilidade Hx(j,:)-(fim),j=',j 
		
!			print*,'Iniciando sensibilidade da componente Hy(j,:)' 
			FDP = 'HyDP'
			call READ_primarios_ADJ_MT( FDP )! ' Leitura do primario DPHy '
		
			call DPEH_3D_EFA_Fadj( f , Es_adj_4 )! ' SOLUÇÂO do DPHy '
	
			call SENSIB_E_H_CSAMT( FDP, f, Es_adj_4, j, MSHy1, MSHy2 ); !print*,' Sensibilidade Hy(j,:)-(fim),j=',j
	
		end do
	END IF

		DEALLOCATE( Es_adj_3 , STAT = IERR )
		IF ( IERR /= 0 ) THEN
			WRITE(*,*)'ERROR IN DEALLOCATING VECTOR (Es_adj_3)'
			WRITE(*,*)'ON THE SUBROUTINE SENSIB_H'
			STOP
		END IF
		DEALLOCATE( Es_adj_4 , STAT = IERR )
		IF ( IERR /= 0 ) THEN
			WRITE(*,*)'ERROR IN DEALLOCATING VECTOR (Es_adj_4)'
			WRITE(*,*)'ON THE SUBROUTINE SENSIB_H'
			STOP
		END IF

	zeta = (0.d0,1.d0)*2.d0*pi*f*mi0
	MSHx1 = MSHx1/(-zeta)
	MSHx2 = MSHx2/(-zeta)
	MSHy1 = MSHy1/(-zeta)
	MSHy2 = MSHy2/(-zeta)

	END SUBROUTINE 

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================


SUBROUTINE SENSIB_H_CSIG( f )
 IMPLICIT NONE
 REAL(dp) , INTENT(in)	:: f
 
 COMPLEX(dp) 	   :: zeta
 CHARACTER(LEN=20) :: FDP
 INTEGER	   :: i, j, IERR

	IF( .NOT. ALLOCATED( MSHx1 ) ) THEN
		allocate( MSHx1(Nobs,Nbloc), MSHx2(Nobs,Nbloc), MSHy1(Nobs,Nbloc), MSHy2(Nobs,Nbloc) )
	
		MSHx1 = (0.d0,0.d0)
		MSHx2 = (0.d0,0.d0)
		MSHy1 = (0.d0,0.d0)
		MSHy2 = (0.d0,0.d0)	
	ELSE	
		
		MSHx1 = (0.d0,0.d0)
		MSHx2 = (0.d0,0.d0)
		MSHy1 = (0.d0,0.d0)
		MSHy2 = (0.d0,0.d0)
	
	END IF	
	
		ALLOCATE( Es_adj_3(Narestas*ngl), STAT = IERR )
			IF ( IERR == 0 )THEN
			Es_adj_3 = (0.d0,0.d0)
		ELSE IF	( IERR /= 0 ) THEN
			WRITE(*,*)'ERROR IN ALLOCATING VECTOR (Es_adj_3)'
			WRITE(*,*)'ON THE SUBROUTINE SENSIB_H'
			STOP
		END IF
		ALLOCATE( Es_adj_4(Narestas*ngl), STAT = IERR )
			IF ( IERR == 0 )THEN
			Es_adj_4 = (0.d0,0.d0)
		ELSE IF ( IERR /= 0 ) THEN
			WRITE(*,*)'ERROR IN ALLOCATING VECTOR (Es_adj_4)'
			WRITE(*,*)'ON THE SUBROUTINE SENSIB_H'
			STOP
		END IF

!	print*,'Iniciando sensibilidade das componente Hx e Hy'

	IF( MT_CSAMT == 'MT')THEN	

		do j=1, Nobs

!			print*,'Iniciando sensibilidade da componente Hx(j,:)' , j
			FDP = 'HxDP'
			call READ_primarios_ADJ_MT( fDP ); ! ' Leitura do primario DPHx '
		
			call DPEH_3D_EFA_Fadj_CSIG( f, Es_adj_3 )! ' SOLUÇÂO do DPHx '
		
			call SENSIB_E_H_1( FDP, f, Es_adj_3, j, MSHx1, MSHx2 )!; print*,' Sensibilidade Hx(j,:)-(fim),j=',j 
		
!			print*,'Iniciando sensibilidade da componente Hy(j,:)' 
			FDP = 'HyDP'
			call READ_primarios_ADJ_MT( FDP )! ' Leitura do primario DPHy '
		
			call DPEH_3D_EFA_Fadj_CSIG( f , Es_adj_4 )! ' SOLUÇÂO do DPHy '
	
			call SENSIB_E_H_1( FDP, f, Es_adj_4, j, MSHy1, MSHy2 ); !print*,' Sensibilidade Hy(j,:)-(fim),j=',j
	
		end do

	ELSEIF(	MT_CSAMT == 'CSAMT' )THEN

		do j=1, Nobs

!			print*,'Iniciando sensibilidade da componente Hx(j,:)' , j
			FDP = 'HxDP'
			call READ_primarios_ADJ_MT( FDP )!; print*, ' Leitura do primario DPHx '
		
			call DPEH_3D_EFA_Fadj_CSIG( f, Es_adj_3 )!; print*, ' SOLUÇÂO do DPHx '
		
			call SENSIB_E_H_CSAMT( FDP, f, Es_adj_3, j, MSHx1, MSHx2 )!; print*,' Sensibilidade Hx(j,:)-(fim),j=',j 
		
!			print*,'Iniciando sensibilidade da componente Hy(j,:)' 
			FDP = 'HyDP'
			call READ_primarios_ADJ_MT( FDP )! ' Leitura do primario DPHy '
		
			call DPEH_3D_EFA_Fadj_CSIG( f , Es_adj_4 )! ' SOLUÇÂO do DPHy '
	
			call SENSIB_E_H_CSAMT( FDP, f, Es_adj_4, j, MSHy1, MSHy2 ); !print*,' Sensibilidade Hy(j,:)-(fim),j=',j
	
		end do
	END IF

		DEALLOCATE( Es_adj_3 , STAT = IERR )
		IF ( IERR /= 0 ) THEN
			WRITE(*,*)'ERROR IN DEALLOCATING VECTOR (Es_adj_3)'
			WRITE(*,*)'ON THE SUBROUTINE SENSIB_H'
			STOP
		END IF
		DEALLOCATE( Es_adj_4 , STAT = IERR )
		IF ( IERR /= 0 ) THEN
			WRITE(*,*)'ERROR IN DEALLOCATING VECTOR (Es_adj_4)'
			WRITE(*,*)'ON THE SUBROUTINE SENSIB_H'
			STOP
		END IF

	zeta = (0.d0,1.d0)*2.d0*pi*f*mi0
	MSHx1 = MSHx1/(-zeta)
	MSHx2 = MSHx2/(-zeta)
	MSHy1 = MSHy1/(-zeta)
	MSHy2 = MSHy2/(-zeta)

	END SUBROUTINE 

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================


SUBROUTINE SENSIB_E( Z_ij, f, MSBcamp, MSBfase )
 IMPLICIT NONE
 REAL(dp),  DIMENSION(:,:), INTENT(out) :: MSBcamp, MSBfase
 REAL(dp),		    INTENT(in)  :: f
 CHARACTER(LEN=20), 	    INTENT(in)  :: Z_ij

 CHARACTER(LEN=20) :: FDP
 INTEGER :: i, j, IERR

 20 format( 2x, 2d25.16)
 if( Z_ij == 'xy' )then
 
	IF( .NOT. ALLOCATED( MSEx1 ) ) THEN	
		allocate( MSEx1(Nobs,Nbloc), MSEx2(Nobs,Nbloc) )
		MSEx1 = (0.d0,0.d0)
		MSEx2 = (0.d0,0.d0)
	ELSE
		MSEx1 = (0.d0,0.d0)
		MSEx2 = (0.d0,0.d0)
	END IF

	ALLOCATE( Es_adj_1(Narestas*ngl), STAT = IERR )
		IF ( IERR == 0 )THEN
		Es_adj_1 = (0.d0,0.d0)
	ELSE IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN ALLOCATING VECTOR (Es_adj_1)'
		WRITE(*,*)'ON THE SUBROUTINE SENSIB_E'
		STOP
	END IF

!	print*,'Iniciando sensibilidade da componente Ex'
	FDP   = 'ExDP'
	IF( MT_CSAMT == 'MT')THEN

		do j=1, Nobs
	
			call READ_primarios_ADJ_MT( FDP )! ' Leitura do primario DPEx '
		
			call DPEH_3D_EFA_Fadj( f, Es_adj_1 )! print*,' SOLUÇÂO do DPEx ' ,j
	
			call SENSIB_E_H_1( FDP, f, Es_adj_1, j, MSEx1, MSEx2 )!' Sensibilidade Ex(j,:)-(fim),j=',j
	
		end do

	ELSEIF( MT_CSAMT == 'CSAMT')THEN

		do j=1, Nobs
	
			call READ_primarios_ADJ_MT( FDP )!;print*,  ' Leitura do primario DPEx '
		
			call DPEH_3D_EFA_Fadj( f, Es_adj_1 )!;print*,' SOLUÇÂO do DPEx ' ,j
	
			call SENSIB_E_H_CSAMT( FDP, f, Es_adj_1, j, MSEx1, MSEx2 )!;print*, ' Sensibilidade Ex(j,:)-(fim),j=',j
	
		end do
	ENDIF


	DEALLOCATE( Es_adj_1 , STAT = IERR )
	IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN DEALLOCATING VECTOR (Es_adj_1)'
		WRITE(*,*)'ON THE SUBROUTINE SENSIB_E'
		STOP
	END IF
      

	call SENSIB_RHOap_FASE( Z_ij, MSEx1, MSEx2, MSBcamp, MSBfase )
     
	deallocate( MSEx1, MSEx2 )


 else if( Z_ij == 'yx' )then

	IF( .NOT. ALLOCATED( MSEy1 ) ) THEN	
		allocate( MSEy1(Nobs,Nbloc), MSEy2(Nobs,Nbloc) )
		MSEy1 = (0.d0,0.d0)
		MSEy2 = (0.d0,0.d0)
	ELSE
		MSEy1 = (0.d0,0.d0)
		MSEy2 = (0.d0,0.d0)
	END IF

	ALLOCATE( Es_adj_2(Narestas*ngl), STAT = IERR )
		IF ( IERR == 0 )THEN
		Es_adj_2 = (0.d0,0.d0)
	ELSE IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN ALLOCATING VECTOR (Es_adj_2)'
		WRITE(*,*)'ON THE SUBROUTINE SENSIB_E'
		STOP
	END IF

!	print*,'Iniciando sensibilidade da componente Ey'
	FDP   = 'EyDP'
	IF( MT_CSAMT == 'MT' )THEN

		do j=1, Nobs
	
			call READ_primarios_ADJ_MT( FDP )! ' Leitura do primario DPEy '
	
			call DPEH_3D_EFA_Fadj( f, Es_adj_2 )! print*, ' SOLUÇÂO do DPEy ', j
	
			call SENSIB_E_H_1( FDP, f, Es_adj_2, j, MSEy1, MSEy2 )!' Sensibilidade Ey(j,:)-(fim),j=',j
	
		end do

	ELSEIF( MT_CSAMT == 'CSAMT' )THEN

		do j=1, Nobs
	
			call READ_primarios_ADJ_MT( FDP )!;print*, ' Leitura do primario DPEy '
	
			call DPEH_3D_EFA_Fadj( f, Es_adj_2 )!;print*, ' SOLUÇÂO do DPEy ', j
	
			call SENSIB_E_H_CSAMT( FDP, f, Es_adj_2, j, MSEy1, MSEy2 )!;print*,' Sensibilidade Ey(j,:)-(fim),j=',j
	
		end do
	ENDIF

	DEALLOCATE( Es_adj_2 , STAT = IERR )
	IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN DEALLOCATING VECTOR (Es_adj_2)'
		WRITE(*,*)'ON THE SUBROUTINE SENSIB_E'
		STOP
	END IF

	call SENSIB_RHOap_FASE( Z_ij, MSEy1, MSEy2, MSBcamp, MSBfase )

	deallocate( MSEy1, MSEy2 )

 end if

END SUBROUTINE

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================


SUBROUTINE SENSIB_E_CSIG( Z_ij, f, MSBcamp, MSBfase )
 IMPLICIT NONE
 REAL(dp),  DIMENSION(:,:), INTENT(out) :: MSBcamp, MSBfase
 REAL(dp),		    INTENT(in)  :: f
 CHARACTER(LEN=20), 	    INTENT(in)  :: Z_ij

 CHARACTER(LEN=20) :: FDP
 INTEGER :: i, j, IERR

 20 format( 2x, 2d25.16)
 if( Z_ij == 'xy' )then
 
	IF( .NOT. ALLOCATED( MSEx1 ) ) THEN	
		allocate( MSEx1(Nobs,Nbloc), MSEx2(Nobs,Nbloc) )
		MSEx1 = (0.d0,0.d0)
		MSEx2 = (0.d0,0.d0)
	ELSE
		MSEx1 = (0.d0,0.d0)
		MSEx2 = (0.d0,0.d0)
	END IF

	ALLOCATE( Es_adj_1(Narestas*ngl), STAT = IERR )
		IF ( IERR == 0 )THEN
		Es_adj_1 = (0.d0,0.d0)
	ELSE IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN ALLOCATING VECTOR (Es_adj_1)'
		WRITE(*,*)'ON THE SUBROUTINE SENSIB_E'
		STOP
	END IF

!	print*,'Iniciando sensibilidade da componente Ex'
	FDP   = 'ExDP'
	IF( MT_CSAMT == 'MT')THEN

		do j=1, Nobs
	
			call READ_primarios_ADJ_MT( FDP )! ' Leitura do primario DPEx '
		
			call DPEH_3D_EFA_Fadj_CSIG( f, Es_adj_1 )! print*,' SOLUÇÂO do DPEx ' ,j
	
			call SENSIB_E_H_1( FDP, f, Es_adj_1, j, MSEx1, MSEx2 )!' Sensibilidade Ex(j,:)-(fim),j=',j
	
		end do

	ELSEIF( MT_CSAMT == 'CSAMT')THEN

		do j=1, Nobs
	
			call READ_primarios_ADJ_MT( FDP )!;print*,  ' Leitura do primario DPEx '
		
			call DPEH_3D_EFA_Fadj_CSIG( f, Es_adj_1 )!;print*,' SOLUÇÂO do DPEx ' ,j
	
			call SENSIB_E_H_CSAMT( FDP, f, Es_adj_1, j, MSEx1, MSEx2 )!;print*, ' Sensibilidade Ex(j,:)-(fim),j=',j
	
		end do
	ENDIF


	DEALLOCATE( Es_adj_1 , STAT = IERR )
	IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN DEALLOCATING VECTOR (Es_adj_1)'
		WRITE(*,*)'ON THE SUBROUTINE SENSIB_E'
		STOP
	END IF
      

	call SENSIB_RHOap_FASE( Z_ij, MSEx1, MSEx2, MSBcamp, MSBfase )
     
	deallocate( MSEx1, MSEx2 )


 else if( Z_ij == 'yx' )then

	IF( .NOT. ALLOCATED( MSEy1 ) ) THEN	
		allocate( MSEy1(Nobs,Nbloc), MSEy2(Nobs,Nbloc) )
		MSEy1 = (0.d0,0.d0)
		MSEy2 = (0.d0,0.d0)
	ELSE
		MSEy1 = (0.d0,0.d0)
		MSEy2 = (0.d0,0.d0)
	END IF

	ALLOCATE( Es_adj_2(Narestas*ngl), STAT = IERR )
		IF ( IERR == 0 )THEN
		Es_adj_2 = (0.d0,0.d0)
	ELSE IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN ALLOCATING VECTOR (Es_adj_2)'
		WRITE(*,*)'ON THE SUBROUTINE SENSIB_E'
		STOP
	END IF

	!print*,'Iniciando sensibilidade da componente Ey'
	FDP   = 'EyDP'
	IF( MT_CSAMT == 'MT' )THEN

		do j=1, Nobs
	
			call READ_primarios_ADJ_MT( FDP )! ' Leitura do primario DPEy '
	
			call DPEH_3D_EFA_Fadj_CSIG( f, Es_adj_2 )! print*, ' SOLUÇÂO do DPEy ', j
	
			call SENSIB_E_H_1( FDP, f, Es_adj_2, j, MSEy1, MSEy2 )!' Sensibilidade Ey(j,:)-(fim),j=',j
	
		end do

	ELSEIF( MT_CSAMT == 'CSAMT' )THEN

		do j=1, Nobs
	
			call READ_primarios_ADJ_MT( FDP )!;print*, ' Leitura do primario DPEy '
	
			call DPEH_3D_EFA_Fadj_CSIG( f, Es_adj_2 )!;print*, ' SOLUÇÂO do DPEy ', j
	
			call SENSIB_E_H_CSAMT( FDP, f, Es_adj_2, j, MSEy1, MSEy2 )!;print*,' Sensibilidade Ey(j,:)-(fim),j=',j
	
		end do
	ENDIF

	DEALLOCATE( Es_adj_2 , STAT = IERR )
	IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN DEALLOCATING VECTOR (Es_adj_2)'
		WRITE(*,*)'ON THE SUBROUTINE SENSIB_E'
		STOP
	END IF

	call SENSIB_RHOap_FASE( Z_ij, MSEy1, MSEy2, MSBcamp, MSBfase )

	deallocate( MSEy1, MSEy2 )

 end if

END SUBROUTINE

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================

SUBROUTINE SENSIB_E_H_1( DPR, f, Es_adj, lin, MSF1, MSF2 )

IMPLICIT NONE
COMPLEX(dp), DIMENSION(:,:), INTENT(inout) :: MSF1, MSF2
COMPLEX(dp), DIMENSION(:)  , INTENT(in)	   :: Es_adj
REAL(dp)   ,		     INTENT(in)    :: f
INTEGER    ,		     INTENT(in)    :: lin

 CHARACTER(LEN=20), 	      		INTENT(in)  :: DPR

COMPLEX(dp), DIMENSION(12)	:: Ep_adj, E_adj, E_f1, E_f2
REAL(dp)   , DIMENSION(12,12)	:: MF
REAL(dp)   , DIMENSION(4,4)	:: ml

COMPLEX(dp), DIMENSION(8)	:: Ex1_D, Ey1_D, Ez1_D
COMPLEX(dp), DIMENSION(8)	:: Ex2_D, Ey2_D, Ez2_D
COMPLEX(dp), DIMENSION(8)	:: Exp_D, Eyp_D
COMPLEX(dp), DIMENSION(8)	:: Exp_A, Eyp_A, Ezp_A
COMPLEX(dp), DIMENSION(8)	:: Ex_A , Ey_A , Ez_A

COMPLEX(dp)	:: inte1, inte2, inte_1(3), inte_2(3), EyHy, HxEx
COMPLEX(dp)	:: Ex, Ey, Ez, Hx, Hy, Hz, soma, neta0, zeta
REAL(dp)	:: dsx, dsy, dsz, Iw, w, sigmas(NC)
REAL(dp)	:: x(2), y(2), z(2), xrec, yrec, zrec
REAL(dp)	:: xno(8), yno(8), zno(8)
INTEGER		:: i, j, k, p, nosG_e(8), aresta(12), inde



E_f1   = (0.d0,0.d0)
E_f2   = (0.d0,0.d0)
Ep_adj = (0.d0,0.d0)
E_adj  = (0.d0,0.d0)

Ex1_D  = (0.d0,0.d0)
Ey1_D  = (0.d0,0.d0)
Ez1_D  = (0.d0,0.d0)

Ex2_D  = (0.d0,0.d0)
Ey2_D  = (0.d0,0.d0)
Ez2_D  = (0.d0,0.d0)

Exp_D  = (0.d0,0.d0)
Eyp_D  = (0.d0,0.d0)

Exp_A  = (0.d0,0.d0)
Eyp_A  = (0.d0,0.d0)
Ezp_A  = (0.d0,0.d0)

Ex_A  = (0.d0,0.d0)
Ey_A  = (0.d0,0.d0)
Ez_A  = (0.d0,0.d0)

Iw  = 1.d0
dsx = 1.d0
dsy = 1.d0
dsz = 1.d0

w = (pi+pi)*f
sigmas = 1.d0/roj
neta0 = 1.d-12 
zeta = (0.d0,1.d0)*w*mi0

xrec = Mnos(no_obs(lin),1)
yrec = Mnos(no_obs(lin),2)
zrec = Mnos(no_obs(lin),3)

do k=1, Nbloc

	inte1 = (0.d0,0.d0)
	inte2 = (0.d0,0.d0)

	do j=nelebloc(k), nelebloc(k+1)-1

		inde = elebloc(j)
		nosG_e(:) = Melem(inde,:)
		aresta(:) = Marestas(inde,:)

		x(1:2) = Mnos(nosG_e(1:2),1)
		y(1) = Mnos(nosG_e(1),2)
		y(2) = Mnos(nosG_e(4),2)
		z(1) = Mnos(nosG_e(1),3)
		z(2) = Mnos(nosG_e(5),3)

!		xno(:) = Mnos(nosG_e(:),1)
!		yno(:) = Mnos(nosG_e(:),2)
		zno(:) = Mnos(nosG_e(:),3)
		
		do i=1, 4

			E_f1(i)   = Es_1(aresta(i))
			E_f1(i+4) = Es_1(aresta(i+4))
			E_f1(i+8) = Es_1(aresta(i+8))

			E_f2(i)   = Es_2(aresta(i))			
			E_f2(i+4) = Es_2(aresta(i+4))			
			E_f2(i+8) = Es_2(aresta(i+8))			
			
			E_adj(i)   = Es_adj(aresta(i))
			E_adj(i+4) = Es_adj(aresta(i+4))
			E_adj(i+8) = Es_adj(aresta(i+8))

		end do

		do i=1, 2

			Ex1_D(i+i-1:i+i)   = E_f1(i+i-1)
			Ex1_D(3+i+i:4+i+i) = E_f1(i+i)

			Ex2_D(i+i-1:i+i)   = E_f2(i+i-1)
			Ex2_D(3+i+i:4+i+i) = E_f2(i+i)

			Ex_A(i+i-1:i+i)   = E_adj(i+i-1)
			Ex_A(3+i+i:4+i+i) = E_adj(i+i)

		end do
	
		Ey1_D(1) = E_f1(5)	; Ey2_D(1) = E_f2(5)	; Ey_A(1) = E_adj(5)
		Ey1_D(5) = E_f1(5)	; Ey2_D(5) = E_f2(5)	; Ey_A(5) = E_adj(5)
			
		Ey1_D(2) = E_f1(7)	; Ey2_D(2) = E_f2(7)	; Ey_A(2) = E_adj(7)
		Ey1_D(6) = E_f1(7)	; Ey2_D(6) = E_f2(7)	; Ey_A(6) = E_adj(7) 
			
		Ey1_D(4) = E_f1(6)	; Ey2_D(4) = E_f2(6)	; Ey_A(4) = E_adj(6)
		Ey1_D(8) = E_f1(6)	; Ey2_D(8) = E_f2(6)	; Ey_A(8) = E_adj(6)
			
		Ey1_D(3) = E_f1(8)	; Ey2_D(3) = E_f2(8)	; Ey_A(3) = E_adj(8)
		Ey1_D(7) = E_f1(8)	; Ey2_D(7) = E_f2(8)	; Ey_A(7) = E_adj(8)
		
		
		Ez1_D(1) = E_f1(9)	; Ez2_D(1) = E_f2(9)	; Ez_A(1) = E_adj(9)
		Ez1_D(4) = E_f1(9)	; Ez2_D(4) = E_f2(9) 	; Ez_A(4) = E_adj(9)
			
		Ez1_D(2) = E_f1(10)	; Ez2_D(2) = E_f2(10)	; Ez_A(2) = E_adj(10)
		Ez1_D(3) = E_f1(10)	; Ez2_D(3) = E_f2(10) 	; Ez_A(3) = E_adj(10)
				
		Ez1_D(5) = E_f1(11)	; Ez2_D(5) = E_f2(11)	; Ez_A(5) = E_adj(11)
		Ez1_D(8) = E_f1(11)	; Ez2_D(8) = E_f2(11) 	; Ez_A(8) = E_adj(11)
 			
		Ez1_D(6) = E_f1(12)	; Ez2_D(6) = E_f2(12)	; Ez_A(6) = E_adj(12)
		Ez1_D(7) = E_f1(12)	; Ez2_D(7) = E_f2(12) 	; Ez_A(7) = E_adj(12)

               ! PRIMÁRIO FONTE ORIGINAL.
		do i=1, 8
			call EyHy_TETM( zno(i), f, 2, EyHy, Exp_D(i) )
			call EyHy_TETM( zno(i), f, 1, Eyp_D(i), HxEx )
		end do

		! MT-Polarização 1.

		Ex1_D(1:2) = Ex1_D(1:2) + Exp_D(1:2)
		Ex1_D(7:8) = Ex1_D(7:8) + Exp_D(7:8)

		Ex1_D(3) = Ex1_D(3) + Exp_D(6) 
		Ex1_D(4) = Ex1_D(4) + Exp_D(5) 
		Ex1_D(5) = Ex1_D(5) + Exp_D(4) 
		Ex1_D(6) = Ex1_D(6) + Exp_D(3) 

		Ey1_D(1:2) = Ey1_D(1:2) + Eyp_D(1:2)
		Ey1_D(7:8) = Ey1_D(7:8) + Eyp_D(7:8)

		Ey1_D(3) = Ey1_D(3) + Eyp_D(6)
		Ey1_D(4) = Ey1_D(4) + Eyp_D(5)
		Ey1_D(5) = Ey1_D(5) + Eyp_D(4)
		Ey1_D(6) = Ey1_D(6) + Eyp_D(3)

		! MT-Polarização 2.

		Ex2_D(1:2) = Ex2_D(1:2) + (-1.d0*Eyp_D(1:2))
		Ex2_D(7:8) = Ex2_D(7:8) + (-1.d0*Eyp_D(7:8))

		Ex2_D(3) = Ex2_D(3) + (-1.d0*Eyp_D(6))
		Ex2_D(4) = Ex2_D(4) + (-1.d0*Eyp_D(5))
		Ex2_D(5) = Ex2_D(5) + (-1.d0*Eyp_D(4))
		Ex2_D(6) = Ex2_D(6) + (-1.d0*Eyp_D(3))

		Ey2_D(1:2) = Ey2_D(1:2) + Exp_D(1:2)
		Ey2_D(7:8) = Ey2_D(7:8) + Exp_D(7:8)

		Ey2_D(3) = Ey2_D(3) + Exp_D(6)
		Ey2_D(4) = Ey2_D(4) + Exp_D(5)
		Ey2_D(5) = Ey2_D(5) + Exp_D(4)
		Ey2_D(6) = Ey2_D(6) + Exp_D(3)

		! Campo adjunto.
                ! PRIMÁRIO FONTE ADJUNTA.

		p = j+j+j+j + j+j+j+j

		Exp_A(:) = Exp_adj(p-7:p)
		Eyp_A(:) = Eyp_adj(p-7:p)
		Ezp_A(:) = Ezp_adj(p-7:p)


		Ex_A(1:2) = Ex_A(1:2) + Exp_A(1:2)
		Ex_A(7:8) = Ex_A(7:8) + Exp_A(7:8)

		Ex_A(3) = Ex_A(3) + Exp_A(6)
		Ex_A(4) = Ex_A(4) + Exp_A(5)
		Ex_A(5) = Ex_A(5) + Exp_A(4)
		Ex_A(6) = Ex_A(6) + Exp_A(3)

		Ey_A(1:2) = Ey_A(1:2) + Eyp_A(1:2)
		Ey_A(7:8) = Ey_A(7:8) + Eyp_A(7:8)

		Ey_A(3) = Ey_A(3) + Eyp_A(6)
		Ey_A(4) = Ey_A(4) + Eyp_A(5)
		Ey_A(5) = Ey_A(5) + Eyp_A(4)
		Ey_A(6) = Ey_A(6) + Eyp_A(3)

		Ez_A(1:2) = Ez_A(1:2) + Ezp_A(1:2)
		Ez_A(7:8) = Ez_A(7:8) + Ezp_A(7:8)

		Ez_A(3) = Ez_A(3) + Ezp_A(6)
		Ez_A(4) = Ez_A(4) + Ezp_A(5)
		Ez_A(5) = Ez_A(5) + Ezp_A(4)
		Ez_A(6) = Ez_A(6) + Ezp_A(3)

		! CALCULANDO A SENSIBILIDADE DO j-ésimo ELEMENTO.
		
		inte_1(1) = IntVolBlc( Ex1_D, Ex_A, x, y, z )
		inte_1(2) = IntVolBlc( Ey1_D, Ey_A, x, y, z )
		inte_1(3) = IntVolBlc( Ez1_D, Ez_A, x, y, z )

		inte_2(1) = IntVolBlc( Ex2_D, Ex_A, x, y, z )
		inte_2(2) = IntVolBlc( Ey2_D, Ey_A, x, y, z )
		inte_2(3) = IntVolBlc( Ez2_D, Ez_A, x, y, z )

		inte1 = inte1 + sum(inte_1(:)) 
		inte2 = inte2 + sum(inte_2(:)) 

	end do
	
	MSF1(lin,k) = inte1
	MSF2(lin,k) = inte2

end do

deallocate( Exp_adj, Eyp_adj, Ezp_adj )			

END SUBROUTINE

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================

SUBROUTINE SENSIB_E_H_CSAMT( DPR, f, Es_adj, lin, MSF1, MSF2 )

IMPLICIT NONE
COMPLEX(dp), DIMENSION(:,:), INTENT(inout) :: MSF1, MSF2
COMPLEX(dp), DIMENSION(:)  , INTENT(in)	   :: Es_adj
REAL(dp)   ,		     INTENT(in)    :: f
INTEGER    ,		     INTENT(in)    :: lin

 CHARACTER(LEN=20), 	      		INTENT(in)  :: DPR

COMPLEX(dp), DIMENSION(12)	:: Ep_adj, E_adj, E_f1, E_f2
REAL(dp)   , DIMENSION(12,12)	:: MF
REAL(dp)   , DIMENSION(4,4)	:: ml

COMPLEX(dp), DIMENSION(8)	:: Ex1_D, Ey1_D, Ez1_D
COMPLEX(dp), DIMENSION(8)	:: Ex2_D, Ey2_D, Ez2_D
COMPLEX(dp), DIMENSION(8)	:: Exp_D, Eyp_D, Ezp_D
COMPLEX(dp), DIMENSION(8)	:: Exp_A, Eyp_A, Ezp_A
COMPLEX(dp), DIMENSION(8)	:: Ex_A , Ey_A , Ez_A

COMPLEX(dp)	:: inte1, inte2, inte_1(3), inte_2(3), EyHy, HxEx
COMPLEX(dp)	:: Ex, Ey, Ez, Hx, Hy, Hz, soma, neta0, zeta
REAL(dp)	:: dsx, dsy, dsz, Iw, w, sigmas(NC)
REAL(dp)	:: x(2), y(2), z(2), xrec, yrec, zrec
REAL(dp)	:: xno(8), yno(8), zno(8)
INTEGER		:: i, j, k, p, nosG_e(8), aresta(12), inde



E_f1   = (0.d0,0.d0)
E_f2   = (0.d0,0.d0)
Ep_adj = (0.d0,0.d0)
E_adj  = (0.d0,0.d0)

Ex1_D  = (0.d0,0.d0)
Ey1_D  = (0.d0,0.d0)
Ez1_D  = (0.d0,0.d0)

Ex2_D  = (0.d0,0.d0)
Ey2_D  = (0.d0,0.d0)
Ez2_D  = (0.d0,0.d0)

Exp_D  = (0.d0,0.d0)
Eyp_D  = (0.d0,0.d0)

Exp_A  = (0.d0,0.d0)
Eyp_A  = (0.d0,0.d0)
Ezp_A  = (0.d0,0.d0)

Ex_A  = (0.d0,0.d0)
Ey_A  = (0.d0,0.d0)
Ez_A  = (0.d0,0.d0)

Iw  = 1.d0
dsx = 1.d0
dsy = 1.d0
dsz = 1.d0

w = (pi+pi)*f
sigmas = 1.d0/roj
neta0 = 1.d-12 
zeta = (0.d0,1.d0)*w*mi0

xrec = Mnos(no_obs(lin),1)
yrec = Mnos(no_obs(lin),2)
zrec = Mnos(no_obs(lin),3)

do k=1, Nbloc

	inte1 = (0.d0,0.d0)
	inte2 = (0.d0,0.d0)

	do j=nelebloc(k), nelebloc(k+1)-1

		inde = elebloc(j)
		nosG_e(:) = Melem(inde,:)
		aresta(:) = Marestas(inde,:)

		x(1:2) = Mnos(nosG_e(1:2),1)
		y(1) = Mnos(nosG_e(1),2)
		y(2) = Mnos(nosG_e(4),2)
		z(1) = Mnos(nosG_e(1),3)
		z(2) = Mnos(nosG_e(5),3)

		do i=1, 4

			E_f1(i)   = Es_1(aresta(i))
			E_f1(i+4) = Es_1(aresta(i+4))
			E_f1(i+8) = Es_1(aresta(i+8))

			E_f2(i)   = Es_2(aresta(i))			
			E_f2(i+4) = Es_2(aresta(i+4))			
			E_f2(i+8) = Es_2(aresta(i+8))			
			
			E_adj(i)   = Es_adj(aresta(i))
			E_adj(i+4) = Es_adj(aresta(i+4))
			E_adj(i+8) = Es_adj(aresta(i+8))

		end do

		do i=1, 2

			Ex1_D(i+i-1:i+i)   = E_f1(i+i-1)
			Ex1_D(3+i+i:4+i+i) = E_f1(i+i)

			Ex2_D(i+i-1:i+i)   = E_f2(i+i-1)
			Ex2_D(3+i+i:4+i+i) = E_f2(i+i)

			Ex_A(i+i-1:i+i)   = E_adj(i+i-1)
			Ex_A(3+i+i:4+i+i) = E_adj(i+i)

		end do
	
		Ey1_D(1) = E_f1(5)	; Ey2_D(1) = E_f2(5)	; Ey_A(1) = E_adj(5)
		Ey1_D(5) = E_f1(5)	; Ey2_D(5) = E_f2(5)	; Ey_A(5) = E_adj(5)
			
		Ey1_D(2) = E_f1(7)	; Ey2_D(2) = E_f2(7)	; Ey_A(2) = E_adj(7)
		Ey1_D(6) = E_f1(7)	; Ey2_D(6) = E_f2(7)	; Ey_A(6) = E_adj(7) 
			
		Ey1_D(4) = E_f1(6)	; Ey2_D(4) = E_f2(6)	; Ey_A(4) = E_adj(6)
		Ey1_D(8) = E_f1(6)	; Ey2_D(8) = E_f2(6)	; Ey_A(8) = E_adj(6)
			
		Ey1_D(3) = E_f1(8)	; Ey2_D(3) = E_f2(8)	; Ey_A(3) = E_adj(8)
		Ey1_D(7) = E_f1(8)	; Ey2_D(7) = E_f2(8)	; Ey_A(7) = E_adj(8)
		
		
		Ez1_D(1) = E_f1(9)	; Ez2_D(1) = E_f2(9)	; Ez_A(1) = E_adj(9)
		Ez1_D(4) = E_f1(9)	; Ez2_D(4) = E_f2(9) 	; Ez_A(4) = E_adj(9)
			
		Ez1_D(2) = E_f1(10)	; Ez2_D(2) = E_f2(10)	; Ez_A(2) = E_adj(10)
		Ez1_D(3) = E_f1(10)	; Ez2_D(3) = E_f2(10) 	; Ez_A(3) = E_adj(10)
				
		Ez1_D(5) = E_f1(11)	; Ez2_D(5) = E_f2(11)	; Ez_A(5) = E_adj(11)
		Ez1_D(8) = E_f1(11)	; Ez2_D(8) = E_f2(11) 	; Ez_A(8) = E_adj(11)
 			
		Ez1_D(6) = E_f1(12)	; Ez2_D(6) = E_f2(12)	; Ez_A(6) = E_adj(12)
		Ez1_D(7) = E_f1(12)	; Ez2_D(7) = E_f2(12) 	; Ez_A(7) = E_adj(12)


		p = j+j+j+j + j+j+j+j


               ! CSAMT - Polarização 1.

		Exp_D(:) = Exp_trans1(p-7:p)
		Eyp_D(:) = Eyp_trans1(p-7:p)
		Ezp_D(:) = Ezp_trans1(p-7:p)

		Ex1_D(1:2) = Ex1_D(1:2) + Exp_D(1:2)
		Ex1_D(7:8) = Ex1_D(7:8) + Exp_D(7:8)

		Ex1_D(3) = Ex1_D(3) + Exp_D(6) 
		Ex1_D(4) = Ex1_D(4) + Exp_D(5) 
		Ex1_D(5) = Ex1_D(5) + Exp_D(4) 
		Ex1_D(6) = Ex1_D(6) + Exp_D(3) 

		Ey1_D(1:2) = Ey1_D(1:2) + Eyp_D(1:2)
		Ey1_D(7:8) = Ey1_D(7:8) + Eyp_D(7:8)

		Ey1_D(3) = Ey1_D(3) + Eyp_D(6)
		Ey1_D(4) = Ey1_D(4) + Eyp_D(5)
		Ey1_D(5) = Ey1_D(5) + Eyp_D(4)
		Ey1_D(6) = Ey1_D(6) + Eyp_D(3)

		Ez1_D(1:2) = Ez1_D(1:2) + Ezp_D(1:2)
		Ez1_D(7:8) = Ez1_D(7:8) + Ezp_D(7:8)

		Ez1_D(3) = Ez1_D(3) + Ezp_D(6)
		Ez1_D(4) = Ez1_D(4) + Ezp_D(5)
		Ez1_D(5) = Ez1_D(5) + Ezp_D(4)
		Ez1_D(6) = Ez1_D(6) + Ezp_D(3)


               ! CSAMT - Polarização 2.

		Exp_D(:) = Exp_trans2(p-7:p)
		Eyp_D(:) = Eyp_trans2(p-7:p)
		Ezp_D(:) = Ezp_trans2(p-7:p)

		Ex2_D(1:2) = Ex2_D(1:2) + Eyp_D(1:2)
		Ex2_D(7:8) = Ex2_D(7:8) + Eyp_D(7:8)

		Ex2_D(3) = Ex2_D(3) + Exp_D(6)
		Ex2_D(4) = Ex2_D(4) + Exp_D(5)
		Ex2_D(5) = Ex2_D(5) + Exp_D(4)
		Ex2_D(6) = Ex2_D(6) + Exp_D(3)

		Ey2_D(1:2) = Ey2_D(1:2) + Exp_D(1:2)
		Ey2_D(7:8) = Ey2_D(7:8) + Exp_D(7:8)

		Ey2_D(3) = Ey2_D(3) + Eyp_D(6)
		Ey2_D(4) = Ey2_D(4) + Eyp_D(5)
		Ey2_D(5) = Ey2_D(5) + Eyp_D(4)
		Ey2_D(6) = Ey2_D(6) + Eyp_D(3)

		Ez2_D(1:2) = Ez2_D(1:2) + Ezp_D(1:2)
		Ez2_D(7:8) = Ez2_D(7:8) + Ezp_D(7:8)

		Ez2_D(3) = Ez2_D(3) + Ezp_D(6)
		Ez2_D(4) = Ez2_D(4) + Ezp_D(5)
		Ez2_D(5) = Ez2_D(5) + Ezp_D(4)
		Ez2_D(6) = Ez2_D(6) + Ezp_D(3)


		! CAMPO ADJUNTO TOTAL.

                ! PRIMÁRIO FONTE ADJUNTA.
		Exp_A(:) = Exp_adj(p-7:p)
		Eyp_A(:) = Eyp_adj(p-7:p)
		Ezp_A(:) = Ezp_adj(p-7:p)


		Ex_A(1:2) = Ex_A(1:2) + Exp_A(1:2)
		Ex_A(7:8) = Ex_A(7:8) + Exp_A(7:8)

		Ex_A(3) = Ex_A(3) + Exp_A(6)
		Ex_A(4) = Ex_A(4) + Exp_A(5)
		Ex_A(5) = Ex_A(5) + Exp_A(4)
		Ex_A(6) = Ex_A(6) + Exp_A(3)

		Ey_A(1:2) = Ey_A(1:2) + Eyp_A(1:2)
		Ey_A(7:8) = Ey_A(7:8) + Eyp_A(7:8)

		Ey_A(3) = Ey_A(3) + Eyp_A(6)
		Ey_A(4) = Ey_A(4) + Eyp_A(5)
		Ey_A(5) = Ey_A(5) + Eyp_A(4)
		Ey_A(6) = Ey_A(6) + Eyp_A(3)


		Ez_A(1:2) = Ez_A(1:2) + Ezp_A(1:2)
		Ez_A(7:8) = Ez_A(7:8) + Ezp_A(7:8)

		Ez_A(3) = Ez_A(3) + Ezp_A(6)
		Ez_A(4) = Ez_A(4) + Ezp_A(5)
		Ez_A(5) = Ez_A(5) + Ezp_A(4)
		Ez_A(6) = Ez_A(6) + Ezp_A(3)

		
		inte_1(1) = IntVolBlc( Ex1_D, Ex_A, x, y, z )
		inte_1(2) = IntVolBlc( Ey1_D, Ey_A, x, y, z )
		inte_1(3) = IntVolBlc( Ez1_D, Ez_A, x, y, z )

		inte_2(1) = IntVolBlc( Ex2_D, Ex_A, x, y, z )
		inte_2(2) = IntVolBlc( Ey2_D, Ey_A, x, y, z )
		inte_2(3) = IntVolBlc( Ez2_D, Ez_A, x, y, z )

		inte1 = inte1 + sum(inte_1(:)) 
		inte2 = inte2 + sum(inte_2(:)) 

	end do
	
	MSF1(lin,k) = inte1
	MSF2(lin,k) = inte2

end do

deallocate( Exp_adj, Eyp_adj, Ezp_adj )			

END SUBROUTINE

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================

SUBROUTINE SENSIB_RHOap_FASE( Z_ij, MSE_1, MSE_2, MSBrho, MSBfase )
IMPLICIT NONE
 REAL(dp)   , DIMENSION(:,:), INTENT(out) :: MSBrho, MSBfase
 COMPLEX(dp), DIMENSION(:,:), INTENT(in)  :: MSE_1, MSE_2
 CHARACTER(LEN=20), 	      INTENT(in)  :: Z_ij

 COMPLEX(dp)	:: det_H, det_HE, dHdsk, dHEdsk, H, HE
 INTEGER	:: i, j

 MSBrho  = 0.d0
 MSBfase = 0.d0


 if( Z_ij == 'xy' )then

	do i=1, Nobs 

		H  = Ht_obs_1(i,1)*Ht_obs_2(i,2)-Ht_obs_1(i,2)*Ht_obs_2(i,1)
		HE = Ht_obs_1(i,1)*Et_obs_2(i,1)-Et_obs_1(i,1)*Ht_obs_2(i,1)

		do j=1, Nbloc
			
			dHdsk  = ( MSHx1(i,j)*Ht_obs_2(i,2) + MSHy2(i,j)*Ht_obs_1(i,1) ) - &
				 ( MSHy1(i,j)*Ht_obs_2(i,1) + MSHx2(i,j)*Ht_obs_1(i,2) )
	
			dHEdsk = ( MSHx1(i,j)*Et_obs_2(i,1) + MSE_2(i,j)*Ht_obs_1(i,1) ) - &
				 ( MSE_1(i,j)*Ht_obs_2(i,1) + MSHx2(i,j)*Et_obs_1(i,1) )
	
			MSBrho(i,j)  = 2.d0*REAL( dHEdsk/HE - dHdsk/H )/rojc(j)  
			MSBfase(i,j) = AIMAG( dHEdsk/HE - dHdsk/H )/rojc(j)

		end do
	end do

 else if( Z_ij == 'yx' )then	

	do i=1, Nobs 

		H  = Ht_obs_1(i,1)*Ht_obs_2(i,2)-Ht_obs_1(i,2)*Ht_obs_2(i,1)
		HE = Et_obs_1(i,2)*Ht_obs_2(i,2)-Et_obs_2(i,2)*Ht_obs_1(i,2)

		do j=1, Nbloc
				
			dHdsk  = ( MSHx1(i,j)*Ht_obs_2(i,2) + MSHy2(i,j)*Ht_obs_1(i,1) ) - &
				 ( MSHy1(i,j)*Ht_obs_2(i,1) + MSHx2(i,j)*Ht_obs_1(i,2) )
	
			dHEdsk = ( MSE_1(i,j)*Ht_obs_2(i,2) + Et_obs_1(i,2)*MSHy2(i,j) ) - &
				 ( MSE_2(i,j)*Ht_obs_1(i,2) + Et_obs_1(i,1)*MSHy1(i,j) )
		
			MSBrho(i,j)  = 2.d0*REAL( dHEdsk/HE - dHdsk/H )/rojc(j)
	
			MSBfase(i,j) = AIMAG( dHEdsk/HE - dHdsk/H )/rojc(j)

		end do
	end do
 end if



END SUBROUTINE

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================

complex(dp) function IntVolBlc( D, A, x, y, z )
!==============================================
!	D : Campo Direto
! 	A : Campo Adjunto
!==============================================
implicit none
complex(dp),intent(in) :: D(8) , A(8)
real(dp),intent(in)    :: x(2),y(2),z(2)


complex(dp) :: D1,D2,D3,D4,D5,D6,D7,D8
complex(dp) :: A1,A2,A3,A4,A5,A6,A7,A8
real(dp)    :: x1,x2,y1,y2,z1,z2

!%Função de integração de volume, sobre a expressão de interpolação
!%trilinear sobre os 8 vértices de um bloco.
!%
!%                   1------------2                   -----> x
!%                  /\           /\                 /|
!%                 / \          / \                / |
!%                /  \         /  \               /  |
!%               /   \        /   \              /   V
!%              /    4-------/----3              y   z
!%             5------------6    /
!%             \   /        \   /
!%             \  /         \  /
!%             \ /          \ /
!%             \/           \/
!%             8------------7 
!%

x1 = x(1); x2 = x(2); y1 = y(1); y2 = y(2); z1 = z(1); z2 = z(2);
!%
D1 = D(1); D2 = D(2); D3 = D(3); D4 = D(4);
D5 = D(5); D6 = D(6); D7 = D(7); D8 = D(8);
!%
A1 = A(1); A2 = A(2); A3 = A(3); A4 = A(4);
A5 = A(5); A6 = A(6); A7 = A(7); A8 = A(8);
!%

IntVolBlc = -( ( 2 * A3 * D1 + 4 * A4 * D1 + 4 * A5 * D1 + 2 * A6 * D1 + &
              A7 * D1 + 2 * A8 * D1 + 4 * A3 * D2 + 2 * A4 * D2 + 2 * &
              A5 * D2 + 4 * A6 * D2 + 2 * A7 * D2 + A8 * D2 + 8 * A3 * &
              D3 + 4 * A4 * D3 + A5 * D3 + 2 * A6 * D3 + 4 * A7 * D3 + &
              2 * A8 * D3 + 4 * A3 * D4 + 8 * A4 * D4 + 2 * A5 * D4 + &
              A6 * D4 + 2 * A7 * D4 + 4 * A8 * D4 + A3 * D5 + 2 * A4 * &
              D5 + 8 * A5 * D5 + 4 * A6 * D5 + 2 * A7 * D5 + 4 * A8 * &
              D5 + 2 * A3 * D6 + A4 * D6 + 4 * A5 * D6 + 8 * A6 * D6 + &
              4 * A7 * D6 + 2 * A8 * D6 + 4 * A3 * D7 + 2 * A4 * D7 + &
              2 * A5 * D7 + 4 * A6 * D7 + 8 * A7 * D7 + 4 * A8 * D7 + &
              2 * ( A3 + 2 * A4 + 2 * A5 + A6 + 2 * A7 + 4 * A8 ) * D8 + &
              A2 * ( 2 * ( 2 * D1 + 4 * D2 + 2 * D3 + D4 + D5 + 2 * D6 + &
              D7 ) + D8 ) + A1 * ( 8 * D1 + 4 * D2 + 2 * D3 + 4 * D4 + &
              4 * D5 + 2 * D6 + D7 + 2 * D8 ) ) * ( x1 - x2 ) * ( y1 - &
              y2 ) * ( z1 - z2 ) ) / 216;


end function


!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================


SUBROUTINE INICIANDO_PARDISO

 IMPLICIT NONE

	INTEGER	:: n, phase, error

	
	!================================== STARTING PARDISO ==================================

	!
	! mtype  : Define the matrix type
	!
	! maxfct : Maximum number of factors with identical sparsity 
	!
	! structure that must be kept in memory at the same time.
	!
	! mnum : Indicates the actual matrix for the solution phase
	!	
	! n     : Number of equations in the sparse linear systems of equation
	!	  logo, n=Narestas*ngl ou n=nnode*ngl.
	!
	! a     : Array that Contains the non-zero elements of the coefficient matrix 
	!         'A' corresponding to the indices in ja. Portanto, a = val_nnz.

	! ja    : Array ja(*) contains column indices of the sparse matrix A stored 
	!         in compressed sparse row format. Como A = A', ja = row_ind.
	!
	! ia    : Array, dimension (n+1). For i≤n, ia(i) points to the first column 
	!	  index of row i in the array ja in compressed sparse row format. 
	!	  Assim sendo, ia = co_ptr.
	!
	! perm  : Holds the permutation vector of size n or specifies elements used 
	!	  for computing a partial solution
	!
	! nrhs  : Number of right-hand sides that need to be solved for.
	!
	! b     : Array, dimension (n, nrhs). On entry, contains the right-hand side 
	!	  vector/matrix B, which is placed in memory contiguously.
	!
	! msgvl : Message level information. If msglvl = 0 then pardiso generates no 
	!	  output, if msglvl = 1 the solver prints statistical information to 
	!	  the screen.
	!
	! x     :
	!
	! error : The tabled error indicator.
	!
	! phase : Controls the execution of the solver.

	! mtype = 3 define a Complex and structurally symmetric matrix
	! mtype = 6 define a Complex and symmetric matrix
	!

	allocate( pt(64) )
!	solverp = 0
!	mtype = 3
	mtype = 6 

	! Sstarting variable pt (Solver internal data address pointer)
	! e a variável iparm (Based on this value of mtype pardisoinit
	! sets default values for the iparm array.)

!	stop
	call pardisoinit( pt, mtype, iparm )

	maxfct = 1
	mnum = 1
	nrhs = 1
	msglvl = 1
	n = Narestas*ngl

	allocate( perm(n) )


	!'Fill-in reducing ordering for the input matrix. '
	!iparm(2)=', iparm(2)

	!'Pivoting perturbation: This parameter instructs Intel '
	!'MKL PARDISO how to handle small pivots or zero pivots'	 
	! iparm(10)

	!'iparm(11): Intel MKL PARDISO uses a maximum weight matching '
	!'algorithm to permute large elements on the diagonal and to  '
	!'scale so that the diagonal elements are equal to 1 and the '
	!'absolute values of the off-diagonal entries are < .OR. = to 1'
	!'iparm(11)=', iparm(11)

	!'iparm(13): Improved accuracy using (non-) symmetric weighted matching.'
	!'value 0 its disable matching. Default for symmetric indefinite matrices'
	!'iparm(13)=', iparm(13)

	!'The parameter iparm(21): Pivoting for symmetric indefinite matrices. '
	!', iparm(21)

	!'iparm(24) Parallel factorization control. This algorithm generally improves'
	!'scalability in case of parallel factorization on many threads (more than eight).'
	!'iparm(24)==',iparm(24)

	!'iparm(25): Parallel forward/backward solve control. 0 for (par) and 1 for (seq).'
	!'iparm(25)=', iparm(25)

	!'iparm(34): Optimal number of threads for conditional numerical'
	!'reproducibility (CNR) mode. But setting up iparm(60) to 0 in' 
	!'order to use the in-core version and produce numerically repeatable results.'

	!'iparm(63): Size of the minimum OOC memory for numerical factorization and solution.'

END SUBROUTINE

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================

SUBROUTINE FATORACAO_MG_PARDISO

IMPLICIT NONE 
	 COMPLEX(dp), ALLOCATABLE, DIMENSION(:) :: solv
	 REAL(dp) :: tt1, tt2
	 INTEGER :: n, phase, error
	
	 n = Narestas*ngl
	 allocate( solv(n) )
	 solv = (0.d0,0.d0) 

	!===================================== phase 1 ========================================

	! Phase      : Fill-reduction analysis and symbolic factorization ( vary options )
	! phase = 11 : Analysis
	phase = 11 

	call pardiso( pt, maxfct, mnum, mtype, phase, n, val_nz, col_ptr, row_ind, &
				  perm, nrhs, iparm, msglvl, vfonte, solv, error )

!	print*,'======================================================================='
!	if (error .eq. 0) then
!		write (*,*) 'Analysis succeeded'
!	else
!		write(*,*) 'error from Analysis = ', error
!	endif
!	print*,'======================================================================='

	!===================================== phase 2 ========================================

	! phase = 22 : Numerical factorization	
	phase = 22 

	call pardiso( pt, maxfct, mnum, mtype, phase, n, val_nz, col_ptr, row_ind, &
				  perm, nrhs, iparm, msglvl, vfonte, solv, error )

!	print*,'======================================================================='
!	if (error .eq. 0) then
!		write (*,*) 'Factorization succeeded'
!	else
!		write(*,*) 'error from factorization = ', error
!	endif
!	print*,'======================================================================='


	deallocate( solv )

END SUBROUTINE

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================

SUBROUTINE liberando_memory_PARDISO( solv )

IMPLICIT NONE
 	COMPLEX(dp), DIMENSION(:), INTENT(inout):: solv
 
	INTEGER	:: n, phase, error


	n=Narestas*ngl


	!===================================== phase 4 ========================================

	! phase = -1 : Release all internal memory for all matrices
 	phase = -1 

	call pardiso( pt, maxfct, mnum, mtype, phase, n, val_nz, col_ptr, row_ind, &
				  perm, nrhs, iparm, msglvl, vfonte, solv, error )

!	print*,'======================================================================='
!	if (error .eq. 0) then
!		write (*,*) 'Release all internal memory succeeded'
!	else
!		write(*,*) 'error from Release all internal memory = ', error
!	endif
!	print*,'======================================================================='

	!=================================== END PARDISO ======================================
	deallocate( perm, pt, vfonte )

END SUBROUTINE

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================

SUBROUTINE MT_3D_EFA( valf )

IMPLICIT NONE
	REAL(dp), INTENT(in)  :: valf


	INTEGER		  :: modo(2), n, phase, error, i, notnun(1), IERR
	REAL(dp)	  :: t1, t2, tt1, tt2, valp(2)
	CHARACTER(LEN=30) :: modeling, polarizacao! mudar p/ global


	modo(1)	= 1
	modo(2)	= 2
	n = Narestas*ngl

	ALLOCATE( Es_1(n), STAT = IERR )
		IF ( IERR == 0 )THEN
		Es_1 = (0.d0,0.d0)
	ELSE IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN ALLOCATING VECTOR (Es_1)'
		WRITE(*,*)'ON THE SUBROUTINE MT_3D_EFA'
		STOP
	END IF

	ALLOCATE( Es_2(n), STAT = IERR )
		IF ( IERR == 0 )THEN
		Es_2 = (0.d0,0.d0)
	ELSE IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN ALLOCATING VECTOR (Es_2)'
		WRITE(*,*)'ON THE SUBROUTINE MT_3D_EFA'
		STOP
	END IF

	modeling = 'MT_TE'
	call fonte_MT3D_EFA( modeling, valf )

	call Cond_Fronteira_DH_3D( ngl, vbordas, nborda, nobs )

	call FATORACAO_MG_PARDISO	

	!===================================== phase 3 ========================================

	! phase = 33 :Solve, iterative refinement
 	phase = 33 

	call pardiso( pt, maxfct, mnum, mtype, phase, n, val_nz, col_ptr, row_ind, &
				  perm, nrhs, iparm, msglvl, vfonte, Es_1, error )

!	print*,'======================================================================='

	do i=1, Narestas
		if( isnan( real(Es_1(i)) ) .or. isnan( aimag(Es_1(i)) ) )then
			notnun(1) = notnun(1) + 1
		end if
	end do
	if( notnun(1) >= 1 )then
		print*, notnun(1),'Not a number present in the solution for 1st polarization';stop
	end if

!	if (error .eq. 0) then
!		write (*,*) 'Solve, iterative refinement succeeded'
!	else
!		write(*,*) 'error from Solve, iterative refinement = ', error;stop
!	endif
!	print*,'======================================================================='

	!================================== FIM - phase 3 ====================================

	deallocate( vfonte )

	modeling = 'MT_TM'
	call fonte_MT3D_EFA( modeling, valf )

	call Cond_Fronteira_DH_3D( ngl, vbordas, nborda, nobs )

	!===================================== phase 3 ========================================
	
	! phase = 33 :Solve, iterative refinement
 	phase = 33 
	
	call pardiso( pt, maxfct, mnum, mtype, phase, n, val_nz, col_ptr, row_ind, &
				  perm, nrhs, iparm, msglvl, vfonte, Es_2, error )
	
	do i=1, Narestas
		if( isnan( real(Es_2(i)) ) .or. isnan( aimag(Es_2(i)) ) )then
			notnun(1) = notnun(1) + 1
		end if
	end do
	if( notnun(1) >= 1 )then
		print*, notnun(1),'Not a number present in the solution for 2nd polarization' ;stop
	end if
!	print*,'======================================================================='
!	if (error .eq. 0) then
!		write (*,*) 'Solve, iterative refinement succeeded'
!	else
!		write(*,*) 'error from Solve, iterative refinement = ', error;stop
!	endif
!	print*,'======================================================================='
	!================================== FIM - phase 3 ====================================

END SUBROUTINE MT_3D_EFA

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================

SUBROUTINE MT_3D_EFA_CSIG( valf )

IMPLICIT NONE
	REAL(dp), INTENT(in)  :: valf


	INTEGER		  :: modo(2), n, phase, error, i, notnun(1), IERR
	REAL(dp)	  :: t1, t2, tt1, tt2, valp(2)
	CHARACTER(LEN=30) :: modeling, polarizacao! mudar p/ global


	modo(1)	= 1
	modo(2)	= 2
	n = Narestas*ngl

	ALLOCATE( Es_1(n), STAT = IERR )
		IF ( IERR == 0 )THEN
		Es_1 = (0.d0,0.d0)
	ELSE IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN ALLOCATING VECTOR (Es_1)'
		WRITE(*,*)'ON THE SUBROUTINE MT_3D_EFA'
		STOP
	END IF

	ALLOCATE( Es_2(n), STAT = IERR )
		IF ( IERR == 0 )THEN
		Es_2 = (0.d0,0.d0)
	ELSE IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN ALLOCATING VECTOR (Es_2)'
		WRITE(*,*)'ON THE SUBROUTINE MT_3D_EFA'
		STOP
	END IF

	modeling = 'MT_TE'
	call fonte_MT3D_EFA_CSIG( modeling, valf )

	call Cond_Fronteira_DH_3D( ngl, vbordas, nborda, nobs )

	call FATORACAO_MG_PARDISO	

	!===================================== phase 3 ========================================

	! phase = 33 :Solve, iterative refinement
 	phase = 33 

	call pardiso( pt, maxfct, mnum, mtype, phase, n, val_nz, col_ptr, row_ind, &
				  perm, nrhs, iparm, msglvl, vfonte, Es_1, error )

!	print*,'======================================================================='

	do i=1, Narestas
		if( isnan( real(Es_1(i)) ) .or. isnan( aimag(Es_1(i)) ) )then
			notnun(1) = notnun(1) + 1
		end if
	end do
	if( notnun(1) >= 1 )then
		print*, notnun(1),'Not a number present in the solution for 1st polarization';stop
	end if

!	if (error .eq. 0) then
!		write (*,*) 'Solve, iterative refinement succeeded'
!	else
!		write(*,*) 'error from Solve, iterative refinement = ', error;stop
!	endif
!	print*,'======================================================================='

	!================================== FIM - phase 3 ====================================

	deallocate( vfonte )

	modeling = 'MT_TM'

	call fonte_MT3D_EFA_CSIG( modeling, valf )

	call Cond_Fronteira_DH_3D( ngl, vbordas, nborda, nobs )

	!===================================== phase 3 ========================================
	
	! phase = 33 :Solve, iterative refinement
 	phase = 33 
	
	call pardiso( pt, maxfct, mnum, mtype, phase, n, val_nz, col_ptr, row_ind, &
				  perm, nrhs, iparm, msglvl, vfonte, Es_2, error )
	
	do i=1, Narestas
		if( isnan( real(Es_2(i)) ) .or. isnan( aimag(Es_2(i)) ) )then
			notnun(1) = notnun(1) + 1
		end if
	end do
	if( notnun(1) >= 1 )then
		print*, notnun(1),'Not a number present in the solution for 2nd polarization' ;stop
	end if
!	print*,'======================================================================='
	if (error /= 0) then
		write(*,*) 'error from Solve, iterative refinement = ', error;stop
	endif
!	print*,'======================================================================='
	!================================== FIM - phase 3 ====================================

END SUBROUTINE MT_3D_EFA_CSIG

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================

SUBROUTINE DPEH_3D_EFA( FDP, valf, indf, indt, Tx, Ty, Tz, solv )
IMPLICIT NONE

	COMPLEX(dp), ALLOCATABLE, DIMENSION(:) , INTENT(out) :: solv
	REAL(dp)   ,				 INTENT(in)  :: valf, Tx, Ty, Tz
	INTEGER    ,				 INTENT(in)  :: indt, indf
	CHARACTER(LEN=20),			 INTENT(in)  :: FDP	

	REAL(dp)				:: t1, t2, tt1, tt2
	INTEGER					:: n, phase, error, i, notnun(1), IERR
    INTEGER, save 				:: indfold
	
	call fonte_DPEH_EFA_1( FDP, valf, Tx, Ty, Tz )

	call Cond_Fronteira_DH_3D( ngl, vbordas, nborda, nobs )

	if( indf == 1 .and. indt == 1 )then
   		indfold = indf
		call FATORACAO_MG_PARDISO
	elseif( indf > indfold )then
   		indfold = indf
		call FATORACAO_MG_PARDISO
	end if		


	n = Narestas*ngl
	ALLOCATE( solv(n), STAT = IERR )
		IF ( IERR == 0 )THEN
		solv = (0.d0,0.d0)
	ELSE IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN ALLOCATING VECTOR (solv)'
		WRITE(*,*)'ON THE SUBROUTINE DPEH_3D_EFA'
		STOP
	END IF

	!===================================== phase 3 ========================================
	! phase = 33 :Solve, iterative refinement
	phase = 33 
	call pardiso( pt, maxfct, mnum, mtype, phase, n, val_nz, col_ptr, row_ind, &
			  perm, nrhs, iparm, msglvl, vfonte, solv, error )

!	print*,'======================================================================='
	do i=1, Narestas
		if( isnan( real(solv(i)) ) .or. isnan( aimag(solv(i)) ) )then
			notnun(1) = notnun(1) + 1
		end if
	end do
	if( notnun(1) >= 1 )then
		print*, notnun(1),'Not a number present in the solution';stop
	end if
!	print*,'======================================================================='

	if (error /= 0) then
		write(*,*) 'error from Solve, iterative refinement = ', error;stop
	endif
!	print*,'======================================================================='

	!================================== FIM - phase 3 ====================================

END SUBROUTINE

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================


SUBROUTINE BCH_3D_EFA( valf, indf, indt, solv )
IMPLICIT NONE

	COMPLEX(dp), ALLOCATABLE, DIMENSION(:) , INTENT(out) :: solv
	REAL(dp)   ,			         INTENT(in) :: valf
	INTEGER    ,			         INTENT(in)  :: indt, indf

	INTEGER         :: i, n, phase, error, IERR, notnun(1)
        INTEGER, save   :: indfold

!	print*,' Inici... Montagem do Vetor Font...!'	
	call fonte_BCH_EFA_EParmz( valf )

	call Cond_Fronteira_DH_3D( ngl, vbordas, nborda, nobs )
!	print*,' Vfonte e Cond... de Vetor Font...OK!'


	if( indf == 1 .and. indt == 1 )then
   		indfold = indf
		call FATORACAO_MG_PARDISO
	elseif( indf > indfold )then
   		indfold = indf
		call FATORACAO_MG_PARDISO
	end if		


	n = Narestas*ngl
	ALLOCATE( solv(n), STAT = IERR )
		IF ( IERR == 0 )THEN
		solv = (0.d0,0.d0)
	ELSE IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN ALLOCATING VECTOR (solv)'
		WRITE(*,*)'ON THE SUBROUTINE BCH_3D_EFA'
		STOP
	END IF


	!===================================== phase 3 ========================================
	! phase = 33 :Solve, iterative refinement
	phase = 33 
	call pardiso( pt, maxfct, mnum, mtype, phase, n, val_nz, col_ptr, row_ind, &
			  perm, nrhs, iparm, msglvl, vfonte, solv, error )

!	print*,'======================================================================='
	do i=1, Narestas
		if( isnan( real(solv(i)) ) .or. isnan( aimag(solv(i)) ) )then
			notnun(1) = notnun(1) + 1
		end if
	end do
	if( notnun(1) >= 1 )then
		print*, notnun(1),'Not a number present in the solution';stop
	end if
!	print*,'======================================================================='

	if (error /= 0) then
		write(*,*) 'error from Solve, iterative refinement = ', error
	endif
!	print*,'======================================================================='



END SUBROUTINE

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================

SUBROUTINE BCH_3D_EFA_CSIG( valf, indf, indt, solv )
IMPLICIT NONE

	COMPLEX(dp), ALLOCATABLE, DIMENSION(:) , INTENT(out) :: solv
	REAL(dp)   ,			         INTENT(in)  :: valf
	INTEGER    ,			         INTENT(in)  :: indt, indf

	INTEGER         :: i, n, phase, error, IERR, notnun(1)
        INTEGER, save   :: indfold

!	print*,' Inici... Montagem do Vetor Font...!'	
	call fonte_BCH_EFA_EParmz_CSIG( valf )

	call Cond_Fronteira_DH_3D( ngl, vbordas, nborda, nobs )
!	print*,' Vfonte e Cond... de Vetor Font...OK!'

	if( indf == 1 .and. indt == 1 )then
   		indfold = indf
		call FATORACAO_MG_PARDISO
	elseif( indf > indfold )then
   		indfold = indf
		call FATORACAO_MG_PARDISO
	end if		


	n = Narestas*ngl
	ALLOCATE( solv(n), STAT = IERR )
		IF ( IERR == 0 )THEN
		solv = (0.d0,0.d0)
	ELSE IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN ALLOCATING VECTOR (solv)'
		WRITE(*,*)'ON THE SUBROUTINE BCH_3D_EFA'
		STOP
	END IF


	!===================================== phase 3 ========================================
	! phase = 33 :Solve, iterative refinement
	phase = 33 
	call pardiso( pt, maxfct, mnum, mtype, phase, n, val_nz, col_ptr, row_ind, &
			  perm, nrhs, iparm, msglvl, vfonte, solv, error )

	print*,'======================================================================='
	do i=1, Narestas
		if( isnan( real(solv(i)) ) .or. isnan( aimag(solv(i)) ) )then
			notnun(1) = notnun(1) + 1
		end if
	end do
	if( notnun(1) >= 1 )then
		print*, notnun(1),'Not a number present in the solution';stop
	end if
	print*,'======================================================================='

	if (error .eq. 0) then
		write (*,*) 'Solve, iterative refinement succeeded'
	else
		write(*,*) 'error from Solve, iterative refinement = ', error
	endif
	print*,'======================================================================='

	!======================================================================================

	!print*,'======================================================================='
	!print*,' '
	!print*,'DPolo 3D realizado com sucesso para Frequência f=',valf,'Hz'
	!print*,' '
	!print*,'======================================================================='


END SUBROUTINE

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================

SUBROUTINE DPEH_3D_EFA_Transm( FDP, valf, indf, indt, solv )
IMPLICIT NONE

	COMPLEX(dp), ALLOCATABLE, DIMENSION(:) , INTENT(out) :: solv
	REAL(dp)  ,				 INTENT(in)  :: valf
	INTEGER   ,				 INTENT(in)  :: indf, indt 
	CHARACTER(LEN=20),			 INTENT(in)  :: FDP	

	REAL(dp)				:: t1, t2, tt1, tt2
	INTEGER					:: n, phase, error, i, notnun(1),IERR
   	INTEGER, save :: indfold


	call fonte_DPEH_EFA_Transm( valf )

	call Cond_Fronteira_DH_3D( ngl, vbordas, nborda, nobs )!'Vetor fonte e Condição de Fronteira'

	if( indf == 1 .and. indt == 1 )then
   		indfold = indf
		call FATORACAO_MG_PARDISO
	elseif( indf > indfold )then
   		indfold = indf
		call FATORACAO_MG_PARDISO
	end if		

	n = Narestas*ngl
	ALLOCATE( solv(n), STAT = IERR )
		IF ( IERR == 0 )THEN
		solv = (0.d0,0.d0)
	ELSE IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN ALLOCATING VECTOR (solv)'
		WRITE(*,*)'ON THE SUBROUTINE DPEH_3D_EFA_Transm'
		STOP
	END IF

	!===================================== phase 3 ========================================
	! phase = 33 :Solve, iterative refinement
	phase = 33 
	call pardiso( pt, maxfct, mnum, mtype, phase, n, val_nz, col_ptr, row_ind, &
			  perm, nrhs, iparm, msglvl, vfonte, solv, error )

!	print*,'======================================================================='

	do i=1, Narestas
		if( isnan( real(solv(i)) ) .or. isnan( aimag(solv(i)) ) )then
			notnun(1) = notnun(1) + 1
		end if
	end do
	if( notnun(1) >= 1 )then
		print*, notnun(1),'Not a number present in the solution';stop
	end if
	
!	print*,'======================================================================='
!	if (error .eq. 0) then
!		write (*,*) 'Solve, iterative refinement succeeded'
!	else
!		write(*,*) 'error from Solve, iterative refinement = ', error
!	endif
!	print*,'======================================================================='

	!================================== FIM - phase 3 ====================================

END SUBROUTINE

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================


SUBROUTINE DPEH_3D_EFA_Transm_CSIG( FDP, valf, indf, indt, solv )
IMPLICIT NONE

	COMPLEX(dp), ALLOCATABLE, DIMENSION(:) , INTENT(out) :: solv
	REAL(dp)  ,				 INTENT(in)  :: valf
	INTEGER   ,				 INTENT(in)  :: indf, indt
	CHARACTER(LEN=20),			 INTENT(in)  :: FDP	

	REAL(dp)				:: t1, t2, tt1, tt2
	INTEGER					:: n, phase, error, i, notnun(1),IERR
        INTEGER, save :: indfold


	call fonte_DPEH_EFA_Transm_CSIG( valf )

	call Cond_Fronteira_DH_3D( ngl, vbordas, nborda, nobs )!'Vetor fonte e Condição de Fronteira'

	if( indf == 1 .and. indt == 1 )then
   		indfold = indf
		call FATORACAO_MG_PARDISO
	elseif( indf > indfold )then
   		indfold = indf
		call FATORACAO_MG_PARDISO
	end if		


	n = Narestas*ngl
	ALLOCATE( solv(n), STAT = IERR )
		IF ( IERR == 0 )THEN
		solv = (0.d0,0.d0)
	ELSE IF ( IERR /= 0 ) THEN
		WRITE(*,*)'ERROR IN ALLOCATING VECTOR (solv)'
		WRITE(*,*)'ON THE SUBROUTINE DPEH_3D_EFA_Transm'
		STOP
	END IF

	!===================================== phase 3 ========================================
	! phase = 33 :Solve, iterative refinement
	phase = 33 
	call pardiso( pt, maxfct, mnum, mtype, phase, n, val_nz, col_ptr, row_ind, &
			  perm, nrhs, iparm, msglvl, vfonte, solv, error )

!	print*,'======================================================================='
	do i=1, Narestas
		if( isnan( real(solv(i)) ) .or. isnan( aimag(solv(i)) ) )then
			notnun(1) = notnun(1) + 1
		end if
	end do
	if( notnun(1) >= 1 )then
		print*, notnun(1),'Not a number present in the solution';stop
	end if
!	print*,'======================================================================='
	if (error /= 0) then
		write(*,*) 'error from Solve, iterative refinement = ', error
	endif
!	print*,'======================================================================='

	!================================== FIM - phase 3 ====================================

END SUBROUTINE

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================

SUBROUTINE DPEH_3D_EFA_Fadj( valf, solv)

IMPLICIT NONE

	COMPLEX(dp), DIMENSION(:) , INTENT(out)	:: solv
	REAL(dp)   ,		    INTENT(in)	:: valf

	REAL(dp)				:: t1, t2, tt1, tt2
	INTEGER					:: n, phase, error
	

	call fonte_DPEH_EFA_Recep( valf )

	call Cond_Fronteira_DH_3D( ngl, vbordas, nborda, nobs )

	n = Narestas*ngl
	solv   = (0.d0,0.d0)

	!===================================== phase 3 ========================================

	! phase = 33 :Solve, iterative refinement
 	phase = 33 

	call pardiso( pt, maxfct, mnum, mtype, phase, n, val_nz, col_ptr, row_ind, &
				  perm, nrhs, iparm, msglvl, vfonte, solv, error )

!	print*,'======================================================================='
	do i=1, Narestas
		if( isnan( real(solv(i)) ) .or. isnan( aimag(solv(i)) ) )then
			notnun(1) = notnun(1) + 1
		end if
	end do
	if( notnun(1) >= 1 )then
		print*, notnun(1),'Not a number present in the solution';stop
	end if
!	print*,'======================================================================='
	if (error /= 0) then
		write(*,*) 'error from Solve, iterative refinement = ', error
	endif
!	print*,'======================================================================='

	!================================== FIM - phase 3 ====================================

END SUBROUTINE

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================


SUBROUTINE DPEH_3D_EFA_Fadj_CSIG( valf, solv)

IMPLICIT NONE

	COMPLEX(dp), DIMENSION(:) , INTENT(out)	:: solv
	REAL(dp)   ,		    INTENT(in)	:: valf

	REAL(dp)				:: t1, t2, tt1, tt2
	INTEGER					:: n, phase, error
	

	call fonte_DPEH_EFA_Recep_CSIG( valf )

	call Cond_Fronteira_DH_3D( ngl, vbordas, nborda, nobs )

	n = Narestas*ngl
	solv   = (0.d0,0.d0)

	!===================================== phase 3 ========================================

	! phase = 33 :Solve, iterative refinement
 	phase = 33 

	call pardiso( pt, maxfct, mnum, mtype, phase, n, val_nz, col_ptr, row_ind, &
				  perm, nrhs, iparm, msglvl, vfonte, solv, error )

!	print*,'======================================================================='
	do i=1, Narestas
		if( isnan( real(solv(i)) ) .or. isnan( aimag(solv(i)) ) )then
			notnun(1) = notnun(1) + 1
		end if
	end do
	if( notnun(1) >= 1 )then
		print*, notnun(1),'Not a number present in the solution';stop
	end if
!	print*,'======================================================================='

	if (error /= 0) then
		write(*,*) 'error from Solve, iterative refinement = ', error
	endif
!	print*,'======================================================================='

	!================================== FIM - phase 3 ====================================

END SUBROUTINE

!=========================================================================================================
!888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=========================================================================================================

SUBROUTINE  montagem_matriz_EFA( valf )
IMPLICIT NONE
REAL(dp), INTENT(in) :: valf

COMPLEX(dp), DIMENSION(12,12)	:: mloc
REAL(dp), DIMENSION(12,12)	:: ME, MF
REAL(dp), DIMENSION(4,4)	:: mk1, mk2, mk3, ml
REAL(dp) 			:: l(3), x(8), y(8), z(8), sigma_e, w
INTEGER				:: nosG_e(8), areG_e(12), nares_e, i, k

w = 2.d0*pi*valf
ME   = 0.D0
MF   = 0.D0
mloc = (0.d0, 0.d0 )
nares_e = 12

call matriz_k1( mk1 )
call matriz_k2( mk2 )
call matriz_k3( mk3 )
call matriz_ml( ml )

allocate( val_nz(nnz), row_ind(nnz) ) 

val_nz = (0.d0, 0.d0 )
row_ind = 0

do i=1, Nelem

	sigma_e = 1.d0/vrho_e(i)

	nosG_e(:) = Melem(i,1:8)
	areG_e(:) = Marestas(i,:)
	x(:) = Mnos(nosG_e(:),1)
	y(:) = Mnos(nosG_e(:),2)
	z(:) = Mnos(nosG_e(:),3)

	l(1) = dabs(x(2)-x(1))
	l(2) = dabs(y(4)-y(1))
	l(3) = dabs(z(5)-z(1))  

	ME(1:4,1:4)   = (l(1)/6.d0)*( mk1*l(3)/l(2) + mk2*l(2)/l(3) )
	ME(5:8,5:8)   = (l(2)/6.d0)*( mk1*l(1)/l(3) + mk2*l(3)/l(1) )
	ME(9:12,9:12) = (l(3)/6.d0)*( mk1*l(2)/l(1) + mk2*l(1)/l(2) )

	ME(1:4,5:8)  = (-l(3)/6.d0)* mk3
	ME(1:4,9:12) = (-l(2)/6.d0)*transpose(mk3)
	ME(5:8,9:12) = (-l(1)/6.d0)* mk3
	ME(5:8,1:4)  = (-l(3)/6.d0)*transpose(mk3)
	ME(9:12,1:4) = (-l(2)/6.d0)* mk3
	ME(9:12,5:8) = (-l(1)/6.d0)*transpose(mk3)

	MF(1:4,1:4)   = l(1)*l(2)*l(3)*Ml/36.d0
	MF(5:8,5:8)   = l(1)*l(2)*l(3)*Ml/36.d0
	MF(9:12,9:12) = l(1)*l(2)*l(3)*Ml/36.d0

	mloc = ME + ci*w*mi0*sigma_e*MF
	call val_nnz_msim( mloc, areG_e, ngl, nares_e )

end do

	call eliminando_zeros_cm( Narestas, nnz, ngl )
	call ordenacao( Narestas, ngl )

END SUBROUTINE montagem_matriz_EFA

!=================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=================================================================================================

SUBROUTINE  montagem_matriz_EFA_CSIG( valf )
IMPLICIT NONE
REAL(dp), INTENT(in) :: valf

COMPLEX(dp), DIMENSION(12,12)	:: mloc
REAL(dp), DIMENSION(12,12)	:: ME, MF
REAL(dp), DIMENSION(4,4)	:: mk1, mk2, mk3, ml
REAL(dp) 			:: l(3), x(8), y(8), z(8), sigma_e,w
COMPLEX(dp)			:: sigmac_e
INTEGER				:: nosG_e(8), areG_e(12), nares_e, i, k, p

w = 2.d0*pi*valf
ME   = 0.D0
MF   = 0.D0
mloc = (0.d0, 0.d0 )
nares_e = 12
sigma_e = (0.d0, 0.d0 )
!
call matriz_k1( mk1 )
call matriz_k2( mk2 )
call matriz_k3( mk3 )
call matriz_ml( ml )

allocate( val_nz(nnz), row_ind(nnz) ) 

val_nz = (0.d0, 0.d0 )
row_ind = 0

do i=1, Nelem

	sigma_e = 1.d0/vrho_e(i)

	nosG_e(:) = Melem(i,1:8)
	areG_e(:) = Marestas(i,:)
	x(:) = Mnos(nosG_e(:),1)
	y(:) = Mnos(nosG_e(:),2)
	z(:) = Mnos(nosG_e(:),3)

	l(1) = dabs(x(2)-x(1))
	l(2) = dabs(y(4)-y(1))
	l(3) = dabs(z(5)-z(1))  

	ME(1:4,1:4)   = (l(1)/6.d0)*( mk1*l(3)/l(2) + mk2*l(2)/l(3) )
	ME(5:8,5:8)   = (l(2)/6.d0)*( mk1*l(1)/l(3) + mk2*l(3)/l(1) )
	ME(9:12,9:12) = (l(3)/6.d0)*( mk1*l(2)/l(1) + mk2*l(1)/l(2) )

	ME(1:4,5:8)  = (-l(3)/6.d0)*mk3
	ME(1:4,9:12) = (-l(2)/6.d0)*transpose(mk3)
	ME(5:8,9:12) = (-l(1)/6.d0)*mk3
	ME(5:8,1:4)  = (-l(3)/6.d0)*transpose(mk3)
	ME(9:12,1:4) = (-l(2)/6.d0)*mk3
	ME(9:12,5:8) = (-l(1)/6.d0)*transpose(mk3)

	MF(1:4,1:4)   = l(1)*l(2)*l(3)*Ml/36.d0
	MF(5:8,5:8)   = l(1)*l(2)*l(3)*Ml/36.d0
	MF(9:12,9:12) = l(1)*l(2)*l(3)*Ml/36.d0



	if( IPELE(i) == 0 )then

		mloc = ME + ci*w*mi0*sigma_e*MF
		call val_nnz_msim( mloc, areG_e, ngl, nares_e )

	elseif( IPELE(i) == 1 )then

		if( INPUTIP == 1 )then	

			do k=1, Nbloc
				do p=nelebloc(k), nelebloc(k+1)-1
					if( elebloc(p) == i )then
						sigmac_e = sig_iw_Dias( sigma_e, w, CARGBLD(k), TMRLX(k), FCPIP(k), PEQMC(k) )
						go to 10
					end if
				end do
			end do 


		elseif( INPUTIP == 2 )then


			do k=1, Nbloc
				do p=nelebloc(k), nelebloc(k+1)-1
					if( elebloc(p) == i )then
						sigmac_e = sig_iw_cole( sigma_e, w, CARGBLD(k), TMRLX(k), CCFREQ(k) )
						go to 10

					end if
				end do
			end do 

		elseif( INPUTIP == 3 )then


			do k=1, Nbloc
				do p=nelebloc(k), nelebloc(k+1)-1
					if( elebloc(p) == i )then
						sigmac_e = sig_iw_Multcole( sigma_e, w, CARGBLD(k), TMRLX(k), CCFREQ(k), &
											CARGBLD2(k), TMRLX2(k), CCFREQ2(k) )

						go to 10
					end if
				end do
			end do 

		10 end if 

		mloc    = ME + ci*w*mi0*sigmac_e*MF
		call val_nnz_msim( mloc, areG_e, ngl, nares_e )
	end if

end do

	call eliminando_zeros_cm( Narestas, nnz, ngl )
	call ordenacao( Narestas, ngl )

END SUBROUTINE montagem_matriz_EFA_CSIG

!=================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=================================================================================================

SUBROUTINE fonte_MT3D_EFA( modlgm, valf )
IMPLICIT NONE
 REAL(dp),		 INTENT(in)	:: valf
 CHARACTER(LEN=20), INTENT(in)		:: modlgm

 COMPLEX(dp), DIMENSION(12)	:: campo_F, vloc
 COMPLEX(dp)			:: Ephi, EyHy, HxEx, Ex, Ey, Hx, Hy
 REAL(dp), DIMENSION(12,12)	:: MF
 REAL(dp), DIMENSION(4,4)	:: ml
 REAL(dp) 			:: l(3), x(8), y(8), z(8), sigma_e, rho_e, Delt_sigma_e
 REAL(dp) 			:: xm(12), ym(12), zm(12), w
 INTEGER			:: nosG_e(8), areG_e(12), nnos_e, meio, i, j, k, md(2), notnun(1)

IF( .NOT. ALLOCATED( vfonte) ) THEN
	ALLOCATE( vfonte(ngl*Narestas) )
	vfonte  = ( 0.d0, 0.d0 )
ELSE
	vfonte  = ( 0.d0, 0.d0 )
END IF

call matriz_ML( ml )

w = (pi+pi)*valf
vfonte = (0.d0, 0.d0)
campo_F = (0.d0, 0.d0)
vloc	= (0.d0, 0.d0)
MF	= 0.d0
nnos_e  = 8

do i=1, nelem

	rho_e = vrho_e(i) 
	nosG_e(:) = Melem(i,:)
	areG_e(:) = Marestas(i,:)
	x(:) = Mnos(nosG_e(:),1)
	y(:) = Mnos(nosG_e(:),2)
	z(:) = Mnos(nosG_e(:),3)
	l(1) = dabs(x(2)-x(1))
	l(2) = dabs(y(4)-y(1))
	l(3) = dabs(z(5)-z(1))  

	MF(1:4,1:4)   = l(1)*l(2)*l(3)*Ml/36.d0
	MF(5:8,5:8)   = l(1)*l(2)*l(3)*Ml/36.d0
	MF(9:12,9:12) = l(1)*l(2)*l(3)*Ml/36.d0

	do j=1, Nbloc

		if( rho_e == rojc(j) )then

			sigma_e = 1.d0/rho_e
			Delt_sigma_e = Deltaneta( sigma_e, nnos_e, z ) 

			if( Delt_sigma_e == 0.d0 ) exit

			xm(1)  = x(1) + x(2) ; ym(1)  = y(1) + y(2) ; zm(1)  = z(1) + z(2)
			xm(2)  = x(3) + x(4) ; ym(2)  = y(3) + y(4) ; zm(2)  = z(3) + z(4)
			xm(3)  = x(5) + x(6) ; ym(3)  = y(5) + y(6) ; zm(3)  = z(5) + z(6)
			xm(4)  = x(7) + x(8) ; ym(4)  = y(7) + y(8) ; zm(4)  = z(7) + z(8)

			xm(5)  = x(1) + x(4) ; ym(5)  = y(1) + y(4) ; zm(5)  = z(1) + z(4)
			xm(6)  = x(5) + x(8) ; ym(6)  = y(5) + y(8) ; zm(6)  = z(5) + z(8)
			xm(7)  = x(2) + x(3) ; ym(7)  = y(2) + y(3) ; zm(7)  = z(2) + z(3)
			xm(8)  = x(6) + x(7) ; ym(8)  = y(6) + y(7) ; zm(8)  = z(6) + z(7)

			xm(9)  = x(1) + x(5) ; ym(9)  = y(1) + y(5) ; zm(9)  = z(1) + z(5)
			xm(10) = x(2) + x(6) ; ym(10) = y(2) + y(6) ; zm(10) = z(2) + z(6)
			xm(11) = x(4) + x(8) ; ym(11) = y(4) + y(8) ; zm(11) = z(4) + z(8)
			xm(12) = x(3) + x(7) ; ym(12) = y(3) + y(7) ; zm(12) = z(3) + z(7)

			xm = xm/2.d0 ; ym = ym/2.d0 ; zm = zm/2.d0


			if( modlgm == 'MT_TE' )then

				md(1) = 1
				md(2) = 2

				do k=1, 4!noele

					call EyHy_TETM( zm(k), valf, md(2), EyHy, HxEx )
					campo_F(k)   = HxEx
					call EyHy_TETM( zm(k+4), valf, md(1), EyHy, HxEx )
					campo_F(k+4) = EyHy
					campo_F(k+8) = (0.d0,0.d0)
				end do

			else if( modlgm == 'MT_TM' )then

				md(1) = 1
				md(2) = 2

				do k=1, 4!noele
					call EyHy_TETM( zm(k), valf, md(1), EyHy, HxEx )
					campo_F(k)   = -1.d0*EyHy
					call EyHy_TETM( zm(k+4), valf, md(2), EyHy, HxEx )
					campo_F(k+4) = HxEx
					campo_F(k+8) = (0.d0,0d0) 
				end do

			end if

			vloc = ci*w*mi0*Delt_sigma_e*matmul(MF,campo_F)*(-1.d0)

			do k=1, 12*ngl
				vfonte(areG_e(k)) = vfonte(areG_e(k)) + vloc(k)

				if( isnan( real(vloc(k)) ) .or. isnan( aimag(vloc(k)) ) )then
					notnun(1) = notnun(1) + 1
					print*, 'Not a number present in the source vector in the element:', i;stop
				end if
			enddo

			exit
			
		end if
	end do
end do

END SUBROUTINE fonte_MT3D_EFA

!=================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=================================================================================================

SUBROUTINE fonte_MT3D_EFA_CSIG( modlgm, valf )
IMPLICIT NONE
 REAL(dp),		 INTENT(in)	:: valf
 CHARACTER(LEN=20), INTENT(in)		:: modlgm

 COMPLEX(dp), DIMENSION(12)	:: campo_F, vloc
 COMPLEX(dp)			:: Ephi, EyHy, HxEx, Ex, Ey, Hx, Hy, Delt_Csigma_e
 REAL(dp), DIMENSION(12,12)	:: MF
 REAL(dp), DIMENSION(4,4)	:: ml
 REAL(dp) 			:: l(3), x(8), y(8), z(8), sigma_e, rho_e, Delt_sigma_e
 REAL(dp) 			:: xm(12), ym(12), zm(12), w
 INTEGER			:: nosG_e(8), areG_e(12), nnos_e, meio, i, j, k, md(2), notnun(1)

IF( .NOT. ALLOCATED( vfonte) ) THEN
	ALLOCATE( vfonte(ngl*Narestas) )
	vfonte  = ( 0.d0, 0.d0 )
ELSE
	vfonte  = ( 0.d0, 0.d0 )
END IF

call matriz_ML( ml )

w = (pi+pi)*valf
vfonte = (0.d0, 0.d0)
campo_F = (0.d0, 0.d0)
vloc	= (0.d0, 0.d0)
MF	= 0.d0
nnos_e  = 8

do i=1, nelem

	rho_e = vrho_e(i) 
	nosG_e(:) = Melem(i,:)
	areG_e(:) = Marestas(i,:)
	x(:) = Mnos(nosG_e(:),1)
	y(:) = Mnos(nosG_e(:),2)
	z(:) = Mnos(nosG_e(:),3)
	l(1) = dabs(x(2)-x(1))
	l(2) = dabs(y(4)-y(1))
	l(3) = dabs(z(5)-z(1))  

	MF(1:4,1:4)   = l(1)*l(2)*l(3)*Ml/36.d0
	MF(5:8,5:8)   = l(1)*l(2)*l(3)*Ml/36.d0
	MF(9:12,9:12) = l(1)*l(2)*l(3)*Ml/36.d0

	do j=1, Nbloc

		if( rho_e == rojc(j) )then

			sigma_e = 1.d0/rho_e

			if( IPELE(i) == 1 )then
				Delt_csigma_e = CDeltaneta( w, sigma_e, nnos_e, j, z ) 
				if( CDABS(Delt_Csigma_e) == 0.d0  ) exit
			elseif( IPELE(i) == 0 )then
				Delt_sigma_e = Deltaneta( sigma_e, nnos_e, z ) 
				if( Delt_sigma_e == 0.d0  ) exit
			end if

			
			xm(1)  = x(1) + x(2) ; ym(1)  = y(1) + y(2) ; zm(1)  = z(1) + z(2)
			xm(2)  = x(3) + x(4) ; ym(2)  = y(3) + y(4) ; zm(2)  = z(3) + z(4)
			xm(3)  = x(5) + x(6) ; ym(3)  = y(5) + y(6) ; zm(3)  = z(5) + z(6)
			xm(4)  = x(7) + x(8) ; ym(4)  = y(7) + y(8) ; zm(4)  = z(7) + z(8)

			xm(5)  = x(1) + x(4) ; ym(5)  = y(1) + y(4) ; zm(5)  = z(1) + z(4)
			xm(6)  = x(5) + x(8) ; ym(6)  = y(5) + y(8) ; zm(6)  = z(5) + z(8)
			xm(7)  = x(2) + x(3) ; ym(7)  = y(2) + y(3) ; zm(7)  = z(2) + z(3)
			xm(8)  = x(6) + x(7) ; ym(8)  = y(6) + y(7) ; zm(8)  = z(6) + z(7)

			xm(9)  = x(1) + x(5) ; ym(9)  = y(1) + y(5) ; zm(9)  = z(1) + z(5)
			xm(10) = x(2) + x(6) ; ym(10) = y(2) + y(6) ; zm(10) = z(2) + z(6)
			xm(11) = x(4) + x(8) ; ym(11) = y(4) + y(8) ; zm(11) = z(4) + z(8)
			xm(12) = x(3) + x(7) ; ym(12) = y(3) + y(7) ; zm(12) = z(3) + z(7)

			xm = xm/2.d0 ; ym = ym/2.d0 ; zm = zm/2.d0


			if( modlgm == 'MT_TE' )then

				md(1) = 1
				md(2) = 2

				do k=1, 4!noele

					call EyHy_TETM( zm(k), valf, md(2), EyHy, HxEx )
					campo_F(k)   = HxEx
					call EyHy_TETM( zm(k+4), valf, md(1), EyHy, HxEx )
					campo_F(k+4) = EyHy
					campo_F(k+8) = (0.d0,0.d0)
				end do

			else if( modlgm == 'MT_TM' )then

				md(1) = 1
				md(2) = 2

				do k=1, 4!noele
					call EyHy_TETM( zm(k), valf, md(1), EyHy, HxEx )
					campo_F(k)   = -1.d0*EyHy
					call EyHy_TETM( zm(k+4), valf, md(2), EyHy, HxEx )
					campo_F(k+4) = HxEx
					campo_F(k+8) = (0.d0,0d0) 
				end do

			end if

			if( IPELE(i) == 1 )then
				vloc = ci*w*mi0*Delt_csigma_e*matmul(MF,campo_F)*(-1.d0)
			elseif( IPELE(i) == 0 )then
				vloc = ci*w*mi0*Delt_sigma_e*matmul(MF,campo_F)*(-1.d0)
			end if


			do k=1, 12*ngl
				vfonte(areG_e(k)) = vfonte(areG_e(k)) + vloc(k)

				if( isnan( real(vloc(k)) ) .or. isnan( aimag(vloc(k)) ) )then
					notnun(1) = notnun(1) + 1
					print*, 'Not a number present in the source vector in the element:', i;stop
				end if
			enddo

			exit
			
		end if
	end do
end do

END SUBROUTINE fonte_MT3D_EFA_CSIG

!=================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=================================================================================================

SUBROUTINE fonte_BCH_EFA( valf, Tx, Ty, Tz )
IMPLICIT NONE
 REAL(dp),		    INTENT(in)	:: valf, Tx, Ty, Tz

 COMPLEX(dp), DIMENSION(12)	:: campo_F, vloc, FF
 COMPLEX(dp)			:: Ex, Ey, Ez, Ephi, zeta

 REAL(dp), DIMENSION(12,12)	:: MF
 REAL(dp), DIMENSION(4,4)	:: ml
 REAL(dp) 				    :: l(3), x(8), y(8), z(8), sigma_e, rho_e, Delt_sigma_e
 REAL(dp) 				    :: xm(12), ym(12), zm(12), w
 INTEGER				    :: nosG_e(8), areG_e(12), nnos_e, meio, i, j, k, p, bb

IF( .NOT. ALLOCATED( vfonte) ) THEN
	ALLOCATE( vfonte(ngl*Narestas) )
	vfonte  = ( 0.d0, 0.d0 )
ELSE
	vfonte  = ( 0.d0, 0.d0 )
END IF

call matriz_ML( ml )

w = (pi+pi)*valf

vfonte  = (0.d0, 0.d0)
campo_F = (0.d0, 0.d0)
vloc	= (0.d0, 0.d0)
Ex	    = (0.d0, 0.d0)
Ey	    = (0.d0, 0.d0)
Ez   	= (0.d0, 0.d0)

MF	    = 0.d0
nnos_e  = 8
zeta    = ci*w*mi0

do i=1, nelem

	rho_e = vrho_e(i) 
	nosG_e(:) = Melem(i,:)
	areG_e(:) = Marestas(i,:)
	x(:) = Mnos(nosG_e(:),1)
	y(:) = Mnos(nosG_e(:),2)
	z(:) = Mnos(nosG_e(:),3)
	l(1) = dabs(x(2)-x(1))
	l(2) = dabs(y(4)-y(1))
	l(3) = dabs(z(5)-z(1))  

	MF(1:4,1:4)   = l(1)*l(2)*l(3)*Ml/36.d0
	MF(5:8,5:8)   = l(1)*l(2)*l(3)*Ml/36.d0
	MF(9:12,9:12) = l(1)*l(2)*l(3)*Ml/36.d0

	do j=1, Nbloc

		if( rho_e == rojc(j) )then
			sigma_e = 1.d0/rho_e
			Delt_sigma_e = Deltaneta( sigma_e, nnos_e, z ) 

			if( Delt_sigma_e == 0.d0 ) exit

			xm(1)  = x(1) + x(2) ; ym(1)  = y(1) + y(2) ; zm(1)  = z(1) + z(2)
			xm(2)  = x(3) + x(4) ; ym(2)  = y(3) + y(4) ; zm(2)  = z(3) + z(4)
			xm(3)  = x(5) + x(6) ; ym(3)  = y(5) + y(6) ; zm(3)  = z(5) + z(6)
			xm(4)  = x(8) + x(7) ; ym(4)  = y(8) + y(7) ; zm(4)  = z(8) + z(7)

			xm(5)  = x(1) + x(4) ; ym(5)  = y(1) + y(4) ; zm(5)  = z(1) + z(4)
			xm(6)  = x(5) + x(8) ; ym(6)  = y(5) + y(8) ; zm(6)  = z(5) + z(8)
			xm(7)  = x(2) + x(3) ; ym(7)  = y(2) + y(3) ; zm(7)  = z(2) + z(3)
			xm(8)  = x(6) + x(7) ; ym(8)  = y(6) + y(7) ; zm(8)  = z(6) + z(7)

			xm(9)  = x(1) + x(5) ; ym(9)  = y(1) + y(5) ; zm(9)  = z(1) + z(5)
			xm(10) = x(2) + x(6) ; ym(10) = y(2) + y(6) ; zm(10) = z(2) + z(6)
			xm(11) = x(4) + x(8) ; ym(11) = y(4) + y(8) ; zm(11) = z(4) + z(8)
			xm(12) = x(3) + x(7) ; ym(12) = y(3) + y(7) ; zm(12) = z(3) + z(7)

			xm = xm/2.d0 ; ym = ym/2.d0 ; zm = zm/2.d0


			do k=1, 4!noele

				call BCH1DQWE_campo_E( valf, Tx, Ty, Tz, xm(k), ym(k), zm(k), Ephi, Ex, Ey )
				campo_F(k)   = Ex
				
				call BCH1DQWE_campo_E( valf, Tx, Ty, Tz, xm(k+4), ym(k+4), zm(k+4), Ephi, Ex, Ey )
				campo_F(k+4) = Ey

				campo_F(k+8) = (0.d0,0.d0)

			end do

			vloc = zeta*Delt_sigma_e*matmul(MF,campo_F)*(-1.d0)

			do k=1, 12*ngl
				vfonte(areG_e(k)) = vfonte(areG_e(k)) + vloc(k)

				if( isnan( real(campo_F(k)) ) .or. isnan( aimag(campo_F(k)) ) )then
					print*, 'Not a number present in the source vector in the element:', i;stop
				end if
			enddo
			exit
			
		end if
	end do
end do

END SUBROUTINE fonte_BCH_EFA

!=================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=================================================================================================

SUBROUTINE fonte_BCH_EFA_CSIG( valf, Tx, Ty, Tz )
IMPLICIT NONE
 REAL(dp),		    INTENT(in)	:: valf, Tx, Ty, Tz

 COMPLEX(dp), DIMENSION(12)	:: campo_F, vloc, FF
 COMPLEX(dp)		        :: Ex, Ey, Ez, Ephi, zeta
 COMPLEX(dp)			:: Delt_csigma_e 

 REAL(dp), DIMENSION(12,12)	:: MF
 REAL(dp), DIMENSION(4,4)	:: ml
 REAL(dp)                       :: l(3), x(8), y(8), z(8), sigma_e, rho_e, Delt_sigma_e
 REAL(dp)                       :: xm(12), ym(12), zm(12), w
 INTEGER                        :: nosG_e(8), areG_e(12), nnos_e, meio, i, j, k, p, bb

IF( .NOT. ALLOCATED( vfonte ) ) THEN
	ALLOCATE( vfonte(ngl*Narestas) )
	vfonte  = ( 0.d0, 0.d0 )
ELSE
	vfonte  = ( 0.d0, 0.d0 )
END IF

call matriz_ML( ml )

w = (pi+pi)*valf

vfonte  = (0.d0, 0.d0)
campo_F = (0.d0, 0.d0)
vloc	= (0.d0, 0.d0)
Ex	= (0.d0, 0.d0)
Ey	= (0.d0, 0.d0)
Ez   	= (0.d0, 0.d0)

MF	= 0.d0
nnos_e  = 8
zeta    = ci*w*mi0

do i=1, nelem

	rho_e = vrho_e(i) 
	nosG_e(:) = Melem(i,:)
	areG_e(:) = Marestas(i,:)
	x(:) = Mnos(nosG_e(:),1)
	y(:) = Mnos(nosG_e(:),2)
	z(:) = Mnos(nosG_e(:),3)
	l(1) = dabs(x(2)-x(1))
	l(2) = dabs(y(4)-y(1))
	l(3) = dabs(z(5)-z(1))  

	MF(1:4,1:4)   = l(1)*l(2)*l(3)*Ml/36.d0
	MF(5:8,5:8)   = l(1)*l(2)*l(3)*Ml/36.d0
	MF(9:12,9:12) = l(1)*l(2)*l(3)*Ml/36.d0

	do j=1, Nbloc

		if( rho_e == rojc(j) )then

			sigma_e = 1.d0/rho_e

			if( IPELE(i) == 1 )then
				Delt_csigma_e = CDeltaneta( w, sigma_e, nnos_e, j, z ) 
				if( CDABS(Delt_Csigma_e) == 0.d0  ) exit
			elseif( IPELE(i) == 0 )then
				Delt_sigma_e = Deltaneta( sigma_e, nnos_e, z ) 
				if( Delt_sigma_e == 0.d0  ) exit
			end if

			xm(1)  = x(1) + x(2) ; ym(1)  = y(1) + y(2) ; zm(1)  = z(1) + z(2)
			xm(2)  = x(3) + x(4) ; ym(2)  = y(3) + y(4) ; zm(2)  = z(3) + z(4)
			xm(3)  = x(5) + x(6) ; ym(3)  = y(5) + y(6) ; zm(3)  = z(5) + z(6)
			xm(4)  = x(8) + x(7) ; ym(4)  = y(8) + y(7) ; zm(4)  = z(8) + z(7)

			xm(5)  = x(1) + x(4) ; ym(5)  = y(1) + y(4) ; zm(5)  = z(1) + z(4)
			xm(6)  = x(5) + x(8) ; ym(6)  = y(5) + y(8) ; zm(6)  = z(5) + z(8)
			xm(7)  = x(2) + x(3) ; ym(7)  = y(2) + y(3) ; zm(7)  = z(2) + z(3)
			xm(8)  = x(6) + x(7) ; ym(8)  = y(6) + y(7) ; zm(8)  = z(6) + z(7)

			xm(9)  = x(1) + x(5) ; ym(9)  = y(1) + y(5) ; zm(9)  = z(1) + z(5)
			xm(10) = x(2) + x(6) ; ym(10) = y(2) + y(6) ; zm(10) = z(2) + z(6)
			xm(11) = x(4) + x(8) ; ym(11) = y(4) + y(8) ; zm(11) = z(4) + z(8)
			xm(12) = x(3) + x(7) ; ym(12) = y(3) + y(7) ; zm(12) = z(3) + z(7)

			xm = xm/2.d0 ; ym = ym/2.d0 ; zm = zm/2.d0

			do k=1, 4!noele

				call BCH1DQWE_campo_E( valf, Tx, Ty, Tz, xm(k), ym(k), zm(k), Ephi, Ex, Ey )
				campo_F(k)   = Ex
				
				call BCH1DQWE_campo_E( valf, Tx, Ty, Tz, xm(k+4), ym(k+4), zm(k+4), Ephi, Ex, Ey )
				campo_F(k+4) = Ey

				campo_F(k+8) = (0.d0,0.d0)

			end do

			if( IPELE(i) == 1 )then
				vloc = ci*w*mi0*Delt_csigma_e*matmul(MF,campo_F)*(-1.d0)
			elseif( IPELE(i) == 0 )then
				vloc = ci*w*mi0*Delt_sigma_e*matmul(MF,campo_F)*(-1.d0)
			end if

			do k=1, 12*ngl
				vfonte(areG_e(k)) = vfonte(areG_e(k)) + vloc(k)

				if( isnan( real(campo_F(k)) ) .or. isnan( aimag(campo_F(k)) ) )then
					print*, 'Not a number present in the source vector in the element:', i;stop
				end if
			enddo
			exit
			
		end if
	end do
end do

END SUBROUTINE fonte_BCH_EFA_CSIG

!=================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=================================================================================================


SUBROUTINE fonte_BCH_EFA_EParmz( valf )
IMPLICIT NONE
 REAL(dp),	    INTENT(in)	:: valf

 COMPLEX(dp), DIMENSION(12)	:: campo_F, vloc, FF
 COMPLEX(dp), DIMENSION(8)	:: Exp_aux, Eyp_aux
 COMPLEX(dp)			:: Ephi, zeta

 REAL(dp), DIMENSION(12,12)	:: MF
 REAL(dp), DIMENSION(4,4)	:: ml
 REAL(dp) 				    :: l(3), x(8), y(8), z(8), sigma_e, rho_e, Delt_sigma_e
 REAL(dp) 				    :: w
 INTEGER				    :: nosG_e(8), areG_e(12), nnos_e, meio, i, j, k, p, bb

IF( .NOT. ALLOCATED( vfonte) ) THEN
	ALLOCATE( vfonte(ngl*Narestas) )
	vfonte  = ( 0.d0, 0.d0 )
ELSE
	vfonte  = ( 0.d0, 0.d0 )
END IF

call matriz_ML( ml )

w = (pi+pi)*valf

vfonte  = (0.d0, 0.d0)
campo_F = (0.d0, 0.d0)
vloc	= (0.d0, 0.d0)

MF	= 0.d0
nnos_e  = 8
zeta    = ci*w*mi0

do i=1, nelem

	rho_e = vrho_e(i) 
	nosG_e(:) = Melem(i,:)
	areG_e(:) = Marestas(i,:)
	x(:) = Mnos(nosG_e(:),1)
	y(:) = Mnos(nosG_e(:),2)
	z(:) = Mnos(nosG_e(:),3)
	l(1) = dabs(x(2)-x(1))
	l(2) = dabs(y(4)-y(1))
	l(3) = dabs(z(5)-z(1))  

	MF(1:4,1:4)   = l(1)*l(2)*l(3)*Ml/36.d0
	MF(5:8,5:8)   = l(1)*l(2)*l(3)*Ml/36.d0
	MF(9:12,9:12) = l(1)*l(2)*l(3)*Ml/36.d0

	do j=1, Nbloc

		if( rho_e == rojc(j) )then

			sigma_e = 1.d0/rho_e
			Delt_sigma_e = Deltaneta( sigma_e, nnos_e, z ) 

			if( Delt_sigma_e == 0.d0 ) exit


			do k=nelebloc(j), nelebloc(Nbloc+1)-1 		
				
				if( elebloc(k) == i )then
					p = k+k+k+k + k+k+k+k
					exit
				else if( k == nelebloc(j+1)-1 )then		
					Go to 10
				end if
			end do
		
			Exp_aux(:) = Exp_trans(p-7:p)
			Eyp_aux(:) = Eyp_trans(p-7:p)

			campo_F(1) = (Exp_aux(1) + Exp_aux(2))/2.d0
			campo_F(2) = (Exp_aux(3) + Exp_aux(4))/2.d0
			campo_F(3) = (Exp_aux(5) + Exp_aux(6))/2.d0
			campo_F(4) = (Exp_aux(7) + Exp_aux(8))/2.d0

			campo_F(5) = (Eyp_aux(1) + Eyp_aux(4))/2.d0
			campo_F(6) = (Eyp_aux(5) + Eyp_aux(8))/2.d0
			campo_F(7) = (Eyp_aux(2) + Eyp_aux(3))/2.d0
			campo_F(8) = (Eyp_aux(6) + Eyp_aux(7))/2.d0

			vloc = zeta*Delt_sigma_e*matmul(MF,campo_F)*(-1.d0)
			
			do k=1, 12*ngl
				vfonte(areG_e(k)) = vfonte(areG_e(k)) + vloc(k)
			enddo
			exit
			
		end if
	10 end do
end do

END SUBROUTINE fonte_BCH_EFA_EParmz

!=================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=================================================================================================

SUBROUTINE fonte_BCH_EFA_EParmz_CSIG( valf )
IMPLICIT NONE
 REAL(dp),		    INTENT(in)	:: valf

 COMPLEX(dp), DIMENSION(12)	:: campo_F, vloc, FF
 COMPLEX(dp), DIMENSION(8)	:: Exp_aux, Eyp_aux
 COMPLEX(dp)			:: Delt_csigma_e, zeta 

 REAL(dp), DIMENSION(12,12)	:: MF
 REAL(dp), DIMENSION(4,4)	:: ml
 REAL(dp)                       :: l(3), x(8), y(8), z(8), sigma_e, rho_e, Delt_sigma_e
 REAL(dp)                       :: xm(12), ym(12), zm(12), w
 INTEGER                        :: nosG_e(8), areG_e(12), nnos_e, meio, i, j, k, p, bb

IF( .NOT. ALLOCATED( vfonte) ) THEN
	ALLOCATE( vfonte(ngl*Narestas) )
	vfonte  = ( 0.d0, 0.d0 )
ELSE
	vfonte  = ( 0.d0, 0.d0 )
END IF

call matriz_ML( ml )

w = (pi+pi)*valf

vfonte  = (0.d0, 0.d0)
campo_F = (0.d0, 0.d0)
vloc	= (0.d0, 0.d0)

MF	    = 0.d0
nnos_e  = 8
zeta    = ci*w*mi0
Delt_csigma_e = (0.D0,0.D0)
Delt_sigma_e = 0.D0

do i=1, nelem

	rho_e = vrho_e(i) 
	nosG_e(:) = Melem(i,:)
	areG_e(:) = Marestas(i,:)
	x(:) = Mnos(nosG_e(:),1)
	y(:) = Mnos(nosG_e(:),2)
	z(:) = Mnos(nosG_e(:),3)
	l(1) = dabs(x(2)-x(1))
	l(2) = dabs(y(4)-y(1))
	l(3) = dabs(z(5)-z(1))  

	MF(1:4,1:4)   = l(1)*l(2)*l(3)*Ml/36.d0
	MF(5:8,5:8)   = l(1)*l(2)*l(3)*Ml/36.d0
	MF(9:12,9:12) = l(1)*l(2)*l(3)*Ml/36.d0

	do j=1, Nbloc

		if( rho_e == rojc(j) )then

			sigma_e = 1.d0/rho_e

			if( IPELE(i) == 1 )then
				Delt_csigma_e = CDeltaneta( w, sigma_e, nnos_e, j, z )
				if( CDABS(Delt_Csigma_e) == 0.d0  ) exit
			elseif( IPELE(i) == 0 )then
				Delt_sigma_e = Deltaneta( sigma_e, nnos_e, z ) 
				if( Delt_sigma_e == 0.d0  ) exit
			end if

			do k=nelebloc(j), nelebloc(Nbloc+1)-1 		
				
				if( elebloc(k) == i )then
					p = k+k+k+k + k+k+k+k
					exit
				else if( k == nelebloc(j+1)-1 )then		
					Go to 10
				end if
			end do
		
			Exp_aux(:) = Exp_trans(p-7:p)
			Eyp_aux(:) = Eyp_trans(p-7:p)

			campo_F(1) = (Exp_aux(1) + Exp_aux(2))/2.d0
			campo_F(2) = (Exp_aux(3) + Exp_aux(4))/2.d0
			campo_F(3) = (Exp_aux(5) + Exp_aux(6))/2.d0
			campo_F(4) = (Exp_aux(7) + Exp_aux(8))/2.d0

			campo_F(5) = (Eyp_aux(1) + Eyp_aux(4))/2.d0
			campo_F(6) = (Eyp_aux(5) + Eyp_aux(8))/2.d0
			campo_F(7) = (Eyp_aux(2) + Eyp_aux(3))/2.d0
			campo_F(8) = (Eyp_aux(6) + Eyp_aux(7))/2.d0

			if( IPELE(i) == 1 )then
				vloc = ci*w*mi0*Delt_csigma_e*matmul(MF,campo_F)*(-1.d0)
			elseif( IPELE(i) == 0 )then
				vloc = ci*w*mi0*Delt_sigma_e*matmul(MF,campo_F)*(-1.d0)
			end if

			do k=1, 12*ngl
				vfonte(areG_e(k)) = vfonte(areG_e(k)) + vloc(k)
			enddo
			exit
			
		end if
	10 end do
end do

END SUBROUTINE fonte_BCH_EFA_EParmz_CSIG

!=================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=================================================================================================

SUBROUTINE fonte_DPEH_EFA_Transm( valf )
IMPLICIT NONE
 REAL(dp),	    INTENT(in)	:: valf

 COMPLEX(dp), DIMENSION(12)	:: campo_F, vloc
 COMPLEX(dp), DIMENSION(8)	:: Exp_aux, Eyp_aux, Ezp_aux
 REAL(dp), DIMENSION(12,12)	:: MF
 REAL(dp), DIMENSION(4,4)	:: ml
 REAL(dp) 			:: w, l(3), x(8), y(8), z(8), sigma_e, rho_e, Delt_sigma_e
 INTEGER			:: nosG_e(8), areG_e(12), nnos_e, meio, i, j, k, p

IF( .NOT. ALLOCATED( vfonte) ) THEN
	ALLOCATE( vfonte(ngl*Narestas) )
	vfonte  = ( 0.d0, 0.d0 )
ELSE
	vfonte  = ( 0.d0, 0.d0 )
END IF

call matriz_ML( ml )

w = (pi+pi)*valf
campo_F = (0.d0, 0.d0)
vloc	= (0.d0, 0.d0)
MF	= 0.d0
nnos_e  = 8

do i=1, nelem

	rho_e = vrho_e(i) 
	nosG_e(:) = Melem(i,:)
	areG_e(:) = Marestas(i,:)
	x(:) = Mnos(nosG_e(:),1)
	y(:) = Mnos(nosG_e(:),2)
	z(:) = Mnos(nosG_e(:),3)
	l(1) = dabs(x(2)-x(1))
	l(2) = dabs(y(4)-y(1))
	l(3) = dabs(z(5)-z(1))  

	MF(1:4,1:4)   = l(1)*l(2)*l(3)*Ml/36.d0
	MF(5:8,5:8)   = l(1)*l(2)*l(3)*Ml/36.d0
	MF(9:12,9:12) = l(1)*l(2)*l(3)*Ml/36.d0

	do j=1, Nbloc

		if( rho_e == rojc(j) )then

			sigma_e = 1.d0/rho_e
			Delt_sigma_e = Deltaneta( sigma_e, nnos_e, z ) 

			if( Delt_sigma_e == 0.d0 ) exit


			do k=nelebloc(j), nelebloc(Nbloc+1)-1 		
				
				if( elebloc(k) == i )then
					p = k+k+k+k + k+k+k+k
					exit
				else if( k == nelebloc(j+1)-1 )then		
					Go to 10
				end if
			end do
		
			Exp_aux(:) = Exp_trans(p-7:p)
			Eyp_aux(:) = Eyp_trans(p-7:p)
			Ezp_aux(:) = Ezp_trans(p-7:p)

			campo_F(1) = (Exp_aux(1) + Exp_aux(2))/2.d0
			campo_F(2) = (Exp_aux(3) + Exp_aux(4))/2.d0
			campo_F(3) = (Exp_aux(5) + Exp_aux(6))/2.d0
			campo_F(4) = (Exp_aux(7) + Exp_aux(8))/2.d0

			campo_F(5) = (Eyp_aux(1) + Eyp_aux(4))/2.d0
			campo_F(6) = (Eyp_aux(5) + Eyp_aux(8))/2.d0
			campo_F(7) = (Eyp_aux(2) + Eyp_aux(3))/2.d0
			campo_F(8) = (Eyp_aux(6) + Eyp_aux(7))/2.d0

			campo_F(9)  = (Ezp_aux(1) + Ezp_aux(5))/2.d0
			campo_F(10) = (Ezp_aux(2) + Ezp_aux(6))/2.d0
			campo_F(11) = (Ezp_aux(4) + Ezp_aux(8))/2.d0
			campo_F(12) = (Ezp_aux(3) + Ezp_aux(7))/2.d0

			vloc = ci*w*mi0*Delt_sigma_e*matmul(MF,campo_F)*(-1.d0)
			
			do k=1, 12*ngl
				vfonte(areG_e(k)) = vfonte(areG_e(k)) + vloc(k)
			enddo
			exit
			
		end if
	10 end do
end do

END SUBROUTINE fonte_DPEH_EFA_Transm

!=================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=================================================================================================

SUBROUTINE fonte_DPEH_EFA_Transm_CSIG( valf )
IMPLICIT NONE
 REAL(dp),	    INTENT(in)	:: valf

 COMPLEX(dp), DIMENSION(12)	:: campo_F, vloc
 COMPLEX(dp), DIMENSION(8)	:: Exp_aux, Eyp_aux, Ezp_aux
 COMPLEX(dp)			:: Delt_csigma_e 
 REAL(dp), DIMENSION(12,12)	:: MF
 REAL(dp), DIMENSION(4,4)	:: ml
 REAL(dp) 			:: w, l(3), x(8), y(8), z(8), sigma_e, rho_e, Delt_sigma_e
 INTEGER			:: nosG_e(8), areG_e(12), nnos_e, meio, i, j, k, p

IF( .NOT. ALLOCATED( vfonte) ) THEN
	ALLOCATE( vfonte(ngl*Narestas) )
	vfonte  = ( 0.d0, 0.d0 )
ELSE
	vfonte  = ( 0.d0, 0.d0 )
END IF

call matriz_ML( ml )

w = (pi+pi)*valf
campo_F = (0.d0, 0.d0)
vloc	= (0.d0, 0.d0)
MF	= 0.d0
nnos_e  = 8
Delt_Csigma_e = (0.D0,0.D0)
Delt_sigma_e = 0.D0

do i=1, nelem

	rho_e = vrho_e(i) 
	nosG_e(:) = Melem(i,:)
	areG_e(:) = Marestas(i,:)
	x(:) = Mnos(nosG_e(:),1)
	y(:) = Mnos(nosG_e(:),2)
	z(:) = Mnos(nosG_e(:),3)
	l(1) = dabs(x(2)-x(1))
	l(2) = dabs(y(4)-y(1))
	l(3) = dabs(z(5)-z(1))  

	MF(1:4,1:4)   = l(1)*l(2)*l(3)*Ml/36.d0
	MF(5:8,5:8)   = l(1)*l(2)*l(3)*Ml/36.d0
	MF(9:12,9:12) = l(1)*l(2)*l(3)*Ml/36.d0

	do j=1, Nbloc

		if( rho_e == rojc(j) )then

			sigma_e = 1.d0/rho_e

			if( IPELE(i) == 1 )then
				Delt_csigma_e = CDeltaneta( w, sigma_e, nnos_e, j, z )
				if( CDABS(Delt_Csigma_e) == 0.d0  ) exit
			elseif( IPELE(i) == 0 )then
				Delt_sigma_e = Deltaneta( sigma_e, nnos_e, z ) 
				if( Delt_sigma_e == 0.d0  ) exit
			end if

			do k=nelebloc(j), nelebloc(Nbloc+1)-1 		
				
				if( elebloc(k) == i )then
					p = k+k+k+k + k+k+k+k
					exit
				else if( k == nelebloc(j+1)-1 )then		
					Go to 10
				end if
			end do
		
			Exp_aux(:) = Exp_trans(p-7:p)
			Eyp_aux(:) = Eyp_trans(p-7:p)
			Ezp_aux(:) = Ezp_trans(p-7:p)

			campo_F(1) = (Exp_aux(1) + Exp_aux(2))/2.d0
			campo_F(2) = (Exp_aux(3) + Exp_aux(4))/2.d0
			campo_F(3) = (Exp_aux(5) + Exp_aux(6))/2.d0
			campo_F(4) = (Exp_aux(7) + Exp_aux(8))/2.d0

			campo_F(5) = (Eyp_aux(1) + Eyp_aux(4))/2.d0
			campo_F(6) = (Eyp_aux(5) + Eyp_aux(8))/2.d0
			campo_F(7) = (Eyp_aux(2) + Eyp_aux(3))/2.d0
			campo_F(8) = (Eyp_aux(6) + Eyp_aux(7))/2.d0

			campo_F(9)  = (Ezp_aux(1) + Ezp_aux(5))/2.d0
			campo_F(10) = (Ezp_aux(2) + Ezp_aux(6))/2.d0
			campo_F(11) = (Ezp_aux(4) + Ezp_aux(8))/2.d0
			campo_F(12) = (Ezp_aux(3) + Ezp_aux(7))/2.d0

			if( IPELE(i) == 1 )then
				vloc = ci*w*mi0*Delt_csigma_e*matmul(MF,campo_F)*(-1.d0)
			elseif( IPELE(i) == 0 )then
				vloc = ci*w*mi0*Delt_sigma_e*matmul(MF,campo_F)*(-1.d0)
			end if

			do k=1, 12*ngl
				vfonte(areG_e(k)) = vfonte(areG_e(k)) + vloc(k)
			enddo
			exit
			
		end if
	10 end do
end do

END SUBROUTINE fonte_DPEH_EFA_Transm_CSIG

!=================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=================================================================================================

SUBROUTINE fonte_DPEH_EFA_recep( valf )
IMPLICIT NONE
 REAL(dp),	    INTENT(in)	:: valf

 COMPLEX(dp), DIMENSION(12)	:: campo_F, vloc
 COMPLEX(dp), DIMENSION(8)	:: Exp_aux, Eyp_aux, Ezp_aux
 REAL(dp), DIMENSION(12,12)	:: MF
 REAL(dp), DIMENSION(4,4)	:: ml
 REAL(dp) 			:: w, l(3), x(8), y(8), z(8), sigma_e, rho_e, Delt_sigma_e
 INTEGER			:: nosG_e(8), areG_e(12), nnos_e, meio, i, j, k, p

IF( .NOT. ALLOCATED( vfonte) ) THEN
	ALLOCATE( vfonte(ngl*Narestas) )
	vfonte  = ( 0.d0, 0.d0 )
ELSE
	vfonte  = ( 0.d0, 0.d0 )
END IF

call matriz_ML( ml )

w = (pi+pi)*valf
campo_F = (0.d0, 0.d0)
vloc	= (0.d0, 0.d0)
MF	= 0.d0
nnos_e  = 8

do i=1, nelem

	rho_e = vrho_e(i) 
	nosG_e(:) = Melem(i,:)
	areG_e(:) = Marestas(i,:)
	x(:) = Mnos(nosG_e(:),1)
	y(:) = Mnos(nosG_e(:),2)
	z(:) = Mnos(nosG_e(:),3)
	l(1) = dabs(x(2)-x(1))
	l(2) = dabs(y(4)-y(1))
	l(3) = dabs(z(5)-z(1))  

	MF(1:4,1:4)   = l(1)*l(2)*l(3)*Ml/36.d0
	MF(5:8,5:8)   = l(1)*l(2)*l(3)*Ml/36.d0
	MF(9:12,9:12) = l(1)*l(2)*l(3)*Ml/36.d0

	do j=1, Nbloc

		if( rho_e == rojc(j) )then

			sigma_e = 1.d0/rho_e
			Delt_sigma_e = Deltaneta( sigma_e, nnos_e, z ) 

			if( Delt_sigma_e == 0.d0 ) exit


			do k=nelebloc(j), nelebloc(Nbloc+1)-1 		
				
				if( elebloc(k) == i )then
					p = k+k+k+k + k+k+k+k
					exit
				else if( k == nelebloc(j+1)-1 )then		
					Go to 10
				end if
			end do

			Exp_aux(:) = Exp_adj(p-7:p)
			Eyp_aux(:) = Eyp_adj(p-7:p)
			Ezp_aux(:) = Ezp_adj(p-7:p)

			campo_F(1) = (Exp_aux(1) + Exp_aux(2))/2.d0
			campo_F(2) = (Exp_aux(3) + Exp_aux(4))/2.d0
			campo_F(3) = (Exp_aux(5) + Exp_aux(6))/2.d0
			campo_F(4) = (Exp_aux(7) + Exp_aux(8))/2.d0

			campo_F(5) = (Eyp_aux(1) + Eyp_aux(4))/2.d0
			campo_F(6) = (Eyp_aux(5) + Eyp_aux(8))/2.d0
			campo_F(7) = (Eyp_aux(2) + Eyp_aux(3))/2.d0
			campo_F(8) = (Eyp_aux(6) + Eyp_aux(7))/2.d0

			campo_F(9)  = (Ezp_aux(1) + Ezp_aux(5))/2.d0
			campo_F(10) = (Ezp_aux(2) + Ezp_aux(6))/2.d0
			campo_F(11) = (Ezp_aux(4) + Ezp_aux(8))/2.d0
			campo_F(12) = (Ezp_aux(3) + Ezp_aux(7))/2.d0


			vloc = ci*w*mi0*Delt_sigma_e*matmul(MF,campo_F)*(-1.d0)

			do k=1, 12*ngl
				vfonte(areG_e(k)) = vfonte(areG_e(k)) + vloc(k)
			enddo
			exit
			
		end if
	10 end do
end do

END SUBROUTINE fonte_DPEH_EFA_recep

!=================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=================================================================================================

SUBROUTINE fonte_DPEH_EFA_recep_CSIG( valf )
IMPLICIT NONE
 REAL(dp),	    INTENT(in)	:: valf

 COMPLEX(dp), DIMENSION(12)	:: campo_F, vloc
 COMPLEX(dp), DIMENSION(8)	:: Exp_aux, Eyp_aux, Ezp_aux
 COMPLEX(dp)			:: Delt_csigma_e 
 REAL(dp), DIMENSION(12,12)	:: MF
 REAL(dp), DIMENSION(4,4)	:: ml
 REAL(dp) 			:: w, l(3), x(8), y(8), z(8), sigma_e, rho_e, Delt_sigma_e
 INTEGER			:: nosG_e(8), areG_e(12), nnos_e, meio, i, j, k, p

IF( .NOT. ALLOCATED( vfonte) ) THEN
	ALLOCATE( vfonte(ngl*Narestas) )
	vfonte  = ( 0.d0, 0.d0 )
ELSE
	vfonte  = ( 0.d0, 0.d0 )
END IF

call matriz_ML( ml )

w = (pi+pi)*valf
campo_F = (0.d0, 0.d0)
vloc	= (0.d0, 0.d0)
MF	= 0.d0
nnos_e  = 8

do i=1, nelem

	rho_e = vrho_e(i) 
	nosG_e(:) = Melem(i,:)
	areG_e(:) = Marestas(i,:)
	x(:) = Mnos(nosG_e(:),1)
	y(:) = Mnos(nosG_e(:),2)
	z(:) = Mnos(nosG_e(:),3)
	l(1) = dabs(x(2)-x(1))
	l(2) = dabs(y(4)-y(1))
	l(3) = dabs(z(5)-z(1))  

	MF(1:4,1:4)   = l(1)*l(2)*l(3)*Ml/36.d0
	MF(5:8,5:8)   = l(1)*l(2)*l(3)*Ml/36.d0
	MF(9:12,9:12) = l(1)*l(2)*l(3)*Ml/36.d0

	do j=1, Nbloc

		if( rho_e == rojc(j) )then

			sigma_e = 1.d0/rho_e
			if( IPELE(i) == 1 )then
				Delt_csigma_e = CDeltaneta( w, sigma_e, nnos_e, j, z ) 
				if( CDABS(Delt_Csigma_e) == 0.d0  ) exit
			elseif( IPELE(i) == 0 )then
				Delt_sigma_e = Deltaneta( sigma_e, nnos_e, z ) 
				if( Delt_sigma_e == 0.d0  ) exit
			end if

			if( Delt_sigma_e == 0.d0 ) exit


			do k=nelebloc(j), nelebloc(Nbloc+1)-1 		
				
				if( elebloc(k) == i )then
					p = k+k+k+k + k+k+k+k
					exit
				else if( k == nelebloc(j+1)-1 )then		
					Go to 10
				end if
			end do

			Exp_aux(:) = Exp_adj(p-7:p)
			Eyp_aux(:) = Eyp_adj(p-7:p)
			Ezp_aux(:) = Ezp_adj(p-7:p)

			campo_F(1) = (Exp_aux(1) + Exp_aux(2))/2.d0
			campo_F(2) = (Exp_aux(3) + Exp_aux(4))/2.d0
			campo_F(3) = (Exp_aux(5) + Exp_aux(6))/2.d0
			campo_F(4) = (Exp_aux(7) + Exp_aux(8))/2.d0

			campo_F(5) = (Eyp_aux(1) + Eyp_aux(4))/2.d0
			campo_F(6) = (Eyp_aux(5) + Eyp_aux(8))/2.d0
			campo_F(7) = (Eyp_aux(2) + Eyp_aux(3))/2.d0
			campo_F(8) = (Eyp_aux(6) + Eyp_aux(7))/2.d0

			campo_F(9)  = (Ezp_aux(1) + Ezp_aux(5))/2.d0
			campo_F(10) = (Ezp_aux(2) + Ezp_aux(6))/2.d0
			campo_F(11) = (Ezp_aux(4) + Ezp_aux(8))/2.d0
			campo_F(12) = (Ezp_aux(3) + Ezp_aux(7))/2.d0


			if( IPELE(i) == 1 )then
				vloc = ci*w*mi0*Delt_csigma_e*matmul(MF,campo_F)*(-1.d0)
			elseif( IPELE(i) == 0 )then
				vloc = ci*w*mi0*Delt_sigma_e*matmul(MF,campo_F)*(-1.d0)
			end if


			do k=1, 12*ngl
				vfonte(areG_e(k)) = vfonte(areG_e(k)) + vloc(k)
			enddo
			exit
			
		end if
	10 end do
end do

END SUBROUTINE fonte_DPEH_EFA_recep_CSIG

!=================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=================================================================================================

FUNCTION Deltaneta( neta_e,  nnose, z ) 
IMPLICIT NONE
REAL(dp), INTENT(in)	:: neta_e, z(:)
INTEGER, INTENT(in)	:: nnose

REAL(dp):: zbc, Deltaneta
INTEGER	:: k 

zbc = sum(z(1: nnose))/nnose

if( NC == 1 )then
	if( zbc <= 0.d0 )then
		Deltaneta = (neta_e - 1.d-12)
	else
		Deltaneta = (neta_e - 1.d0/roj(1))
	end if
else 				
	if( zbc <= 0.d0 )then
		Deltaneta = (neta_e - 1.d-12)			
	else			
		do k=1, NC-1	
			if( zbc > zkmadaj(k-1) .and. zbc < zkmadaj(k) )then
				Deltaneta = (neta_e - 1.d0/roj(k))
			else if( zbc > zkmadaj(NC-1) )then
				Deltaneta = (neta_e - 1.d0/roj(NC))	
			end if
		end do				
	end if
end if

END FUNCTION Deltaneta

!=================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=================================================================================================

COMPLEX(dp) FUNCTION CDeltaneta( w, neta_e,  nnose, indblc, z ) 
IMPLICIT NONE
REAL(dp), INTENT(in)	:: w, neta_e, z(:)
INTEGER, INTENT(in)	:: nnose, indblc

REAL(dp)    :: zbc
COMPLEX(dp) :: sigmac_e
INTEGER	    :: k, p

zbc = sum(z(1: nnose))/nnose

if( NC == 1 )then
	if( zbc <= 0.d0 )then
		CDeltaneta = (neta_e - 1.d-12)
	else
		
		if( INPUTIP == 1 )then

			sigmac_e = sig_iw_Dias( neta_e, w, CARGBLD(indblc), TMRLX(indblc), FCPIP(indblc), PEQMC(indblc) )

		elseif( INPUTIP == 2 )then

			sigmac_e = sig_iw_cole( neta_e, w, CARGBLD(indblc), TMRLX(indblc), CCFREQ(indblc) )

		elseif( INPUTIP == 3 )then

			sigmac_e = sig_iw_Multcole( neta_e, w, CARGBLD(indblc), TMRLX(indblc), CCFREQ(indblc), &
								CARGBLD2(indblc), TMRLX2(indblc), CCFREQ2(indblc) )
		end if

		CDeltaneta = sigmac_e - 1.d0/roj(1)
	end if
else 				
	if( zbc <= 0.d0 )then
		CDeltaneta = (neta_e - 1.d-12)			
	else			
		do k=1, NC-1	
			if( zbc > zkmadaj(k-1) .and. zbc < zkmadaj(k) )then

			
				if( INPUTIP == 1 )then	
					sigmac_e = sig_iw_Dias( neta_e, w, CARGBLD(indblc), TMRLX(indblc), FCPIP(indblc), PEQMC(indblc) )
				elseif( INPUTIP == 2 )then
					sigmac_e = sig_iw_cole( neta_e, w, CARGBLD(indblc), TMRLX(indblc), CCFREQ(indblc) )
				elseif( INPUTIP == 3 )then
					sigmac_e = sig_iw_Multcole( neta_e, w, CARGBLD(indblc), TMRLX(indblc), CCFREQ(indblc), &
									       CARGBLD2(indblc), TMRLX2(indblc), CCFREQ2(indblc) )
				end if

				CDeltaneta = ( sigmac_e - 1.d0/roj(k) )

			else if( zbc > zkmadaj(NC-1) )then

				if( INPUTIP == 1 )then	
					sigmac_e = sig_iw_Dias( neta_e, w, CARGBLD(indblc), TMRLX(indblc), FCPIP(indblc), PEQMC(indblc) )
				elseif( INPUTIP == 2 )then
					sigmac_e = sig_iw_cole( neta_e, w, CARGBLD(indblc), TMRLX(indblc), CCFREQ(indblc) )
				elseif( INPUTIP == 3 )then
					sigmac_e = sig_iw_Multcole( neta_e, w, CARGBLD(indblc), TMRLX(indblc), CCFREQ(indblc), &
									       CARGBLD2(indblc), TMRLX2(indblc), CCFREQ2(indblc) )
				end if

				CDeltaneta = sigmac_e - 1.d0/roj(NC)
			end if
		end do				
	end if
end if

END FUNCTION CDeltaneta

!=================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=================================================================================================

SUBROUTINE coluna_ptr_3D( neqtn )
IMPLICIT NONE
INTEGER, INTENT(in)         		:: neqtn

INTEGER::  aux(2), i

allocate( col_ptr( ngl*neqtn+1 )) ; col_ptr = 0 ;

col_ptr = 40*ngl

nnz = sum(col_ptr(:))
aux(1) = col_ptr(1)
col_ptr(1) = 1

do i=2, ngl*neqtn
	aux(2) = col_ptr(i)
	col_ptr(i) = aux(1) + 1
	aux(1) = aux(1) + aux(2)
end do

col_ptr(ngl*neqtn+1) = aux(1) + 1

nnz_aux = nnz

END SUBROUTINE

!=================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=================================================================================================

SUBROUTINE val_nnz_mestsim( mloc, nos_e, ngl, nnos_e)
COMPLEX(dp), DIMENSION(:,:), INTENT(in)	:: mloc
INTEGER, DIMENSION(:), INTENT(in)	:: nos_e
INTEGER, INTENT(in)		 	:: ngl, nnos_e

INTEGER :: i, j, k, dim_ij

dim_ij = nnos_e*ngl

do j=1, dim_ij
	do i=1, dim_ij
		do k=col_ptr(nos_e(j)), col_ptr(nos_e(j)+1)-1
			if( row_ind(k) == nos_e(i) )then
				val_nz(k) = val_nz(k) + mloc(i,j)
				exit
			else if( row_ind(k) == 0 )then
				val_nz(k) = mloc(i,j)
				row_ind(k)= nos_e(i)
				exit
			else
				cycle				
			end if		
		end do
	end do
end do

END SUBROUTINE

!=================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=================================================================================================

SUBROUTINE val_nnz_msim( mloc, nos_e, ngl, nnos_e)
COMPLEX(dp), DIMENSION(:,:), INTENT(in)	:: mloc
INTEGER, DIMENSION(:), INTENT(in)	:: nos_e
INTEGER, INTENT(in)		 	:: ngl, nnos_e

INTEGER :: i, j, k, dim_ij

dim_ij = nnos_e*ngl

do j=1, dim_ij
	do i=1, dim_ij
		do k=col_ptr(nos_e(j)), col_ptr(nos_e(j)+1)-1
			if( row_ind(k) == nos_e(i) )then
				val_nz(k) = val_nz(k) + mloc(i,j)
				exit
			else if( row_ind(k) == 0 .and. nos_e(i) >= nos_e(j) )then
				val_nz(k) = mloc(i,j)
				row_ind(k)= nos_e(i)
				exit
			else
				cycle				
			end if		
		end do
	end do
end do

END SUBROUTINE

!=================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=================================================================================================

SUBROUTINE eliminando_zeros_cm( neqtn, nnz, ngl )
IMPLICIT NONE
INTEGER, INTENT(in)	:: neqtn
INTEGER, INTENT(inout)  :: ngl, nnz

INTEGER    , ALLOCATABLE, DIMENSION(:) :: row_ind_aux, col_ptr_aux
COMPLEX(dp), ALLOCATABLE, DIMENSION(:)  :: val_nz_aux

INTEGER:: l(2), j, k, nzeros, ng, ncol

nzeros=0 ; l(2)=0 ; l(1)=1 ; ncol = neqtn*ngl
ng = 0

do j=1, ncol

	l(1) = col_ptr(j)
	l(2) = col_ptr(j+1)-1

	do k=l(1), l(2)

		if( k > nnz )then
			print*, 'col', j
			print*, 'k > ', nnz
		end if

		if( cdabs(val_nz(k)) /= 0.d0 )then
			ng = ng + 1
		end if
	end do
end do

allocate( val_nz_aux( ng ), row_ind_aux( ng ), col_ptr_aux(ncol+1) )
val_nz_aux = 0
row_ind_aux = 0
col_ptr_aux = 1

nnz=0 ; nzeros = 0

do j=1, ncol
	l(1) = col_ptr(j)
	l(2) = col_ptr(j+1)-1
	do k=l(1), l(2)
		if( val_nz(k) /= (0.d0, 0.d0) )then
			nnz = nnz + 1
			val_nz_aux(nnz) = val_nz(k)
			row_ind_aux(nnz) = row_ind(k)
			col_ptr_aux(j+1) = col_ptr(j+1) - nzeros
		else
			nzeros = nzeros + 1
			col_ptr_aux(j+1) = col_ptr(j+1) - nzeros
		end if
	end do
end do

deallocate( row_ind, val_nz )

allocate( val_nz( nnz ), row_ind( nnz ) )

val_nz = val_nz_aux
row_ind = row_ind_aux
col_ptr = col_ptr_aux

deallocate( val_nz_aux, row_ind_aux, col_ptr_aux )

END SUBROUTINE

!=================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=================================================================================================

subroutine ordenacao( neqtn, ngl )

!=================================================================================================
! 
! Program that orders the column vector compressed in increasing direction
! 
! 
! 
! neqtn :-------------> Number of system variables
! ngl :---------------> No. of degrees of freedom of the problem
! 
!=================================================================================================

IMPLICIT NONE
integer, intent(in)		     :: neqtn, ngl

integer :: l, k, ntemp, min_ind_tmp
integer :: max_ind_tmp, min_ind, max_ind
complex(dp) :: val_tmp(2)

integer, allocatable, dimension(:):: ind_temp, ind_aux
complex(dp), allocatable, dimension(:):: val_temp
integer	:: ndim, nk

ndim = 1
nk = neqtn*ngl

do l=1, nk

	ntemp = col_ptr(l+1) - col_ptr(l)
	allocate( val_temp(ntemp), ind_temp(ntemp) )
	val_temp(1:ntemp) = (0.d0, 0.d0) 
	ind_temp(1:ntemp) = 0
	val_temp(1:ntemp) = val_nz( col_ptr(l):col_ptr(l+1)-1 )
	ind_temp(1:ntemp) = row_ind( col_ptr(l):col_ptr(l+1)-1 )

	do k=1, int(ntemp/2)

		min_ind = minloc( ind_temp(k:ntemp-(k-1)), ndim )+(K-1)

		min_ind_tmp = ind_temp(k) 
		ind_temp(k) = ind_temp(min_ind)
		ind_temp(min_ind) = min_ind_tmp

		val_tmp(1) = val_temp(k)
		val_temp(k) = val_temp(min_ind)
		val_temp(min_ind) = val_tmp(1)

		max_ind = maxloc( ind_temp(k:ntemp-(k-1)), ndim )+(K-1)

		max_ind_tmp = ind_temp(ntemp-(k-1)) 
		ind_temp(ntemp-(k-1)) = ind_temp(max_ind)
		ind_temp(max_ind) = max_ind_tmp

		val_tmp(2) = val_temp(ntemp-(k-1))
		val_temp(ntemp-(k-1)) = val_temp(max_ind)
		val_temp(max_ind) = val_tmp(2)
	end do
	val_nz( col_ptr(l):col_ptr(l+1)-1 ) = val_temp 
	row_ind( col_ptr(l):col_ptr(l+1)-1 )= ind_temp
	deallocate( val_temp, ind_temp )
end do

end subroutine

!======================================================================================
!88888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!======================================================================================	

SUBROUTINE Cond_Fronteira_DH_3D( ngl, vnos_borda1, nborda, nobs )
IMPLICIT NONE
 INTEGER, INTENT(in)				:: ngl, nborda, nobs
 INTEGER, DIMENSION(:), INTENT(in)		:: vnos_borda1

 INTEGER	:: i, j, k

do i=1, nborda
	do k=1, ngl 
		do j=col_ptr(ngl*vnos_borda1(i)-(k-1)), col_ptr(ngl*vnos_borda1(i)-(k-1)+1)-1
			if( row_ind(j) == ngl*vnos_borda1(i)-(k-1) )then
				val_nz(j) = 1.d30
			end if
			vfonte(ngl*vnos_borda1(i)-(k-1))=0.d0
		end do
	end do
end do

END SUBROUTINE Cond_Fronteira_DH_3D

!=================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=================================================================================================

SUBROUTINE matriz_k1( mk1 )
IMPLICIT NONE 
REAL(dp), DIMENSION(:,:), INTENT(out):: mk1
 

	mk1(1,:) = [ 2.d0, -2.d0,  1.d0, -1.d0 ]
	mk1(2,:) = [-2.d0,  2.d0, -1.d0,  1.d0 ] 
	mk1(3,:) = [ 1.d0, -1.d0,  2.d0, -2.d0 ] 
	mk1(4,:) = [-1.d0,  1.d0, -2.d0,  2.d0 ]  


END SUBROUTINE matriz_k1

!=================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=================================================================================================

SUBROUTINE matriz_k2( mk2 )
IMPLICIT NONE 
REAL(dp), DIMENSION(:,:), INTENT(out):: mk2

	mk2(1,:) = [ 2.d0,  1.d0, -2.d0, -1.d0 ]
	mk2(2,:) = [ 1.d0,  2.d0, -1.d0, -2.d0 ] 
	mk2(3,:) = [-2.d0, -1.d0,  2.d0,  1.d0 ] 
	mk2(4,:) = [-1.d0, -2.d0,  1.d0,  2.d0 ]  

END SUBROUTINE matriz_k2

!=================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=================================================================================================

SUBROUTINE matriz_k3( mk3 )
IMPLICIT NONE 
REAL(dp), DIMENSION(:,:), INTENT(out):: mk3

	mk3(1,:) = [ 2.d0,  1.d0, -2.d0, -1.d0 ]
	mk3(2,:) = [-2.d0, -1.d0,  2.d0,  1.d0 ] 
	mk3(3,:) = [ 1.d0,  2.d0, -1.d0, -2.d0 ] 
	mk3(4,:) = [-1.d0, -2.d0,  1.d0,  2.d0 ]  

END SUBROUTINE matriz_k3

!=================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=================================================================================================

SUBROUTINE matriz_Ml( ml )
IMPLICIT NONE 
REAL(dp), DIMENSION(:,:), INTENT(out):: ml

	ml(1,:) = [ 4.d0, 2.d0, 2.d0, 1.d0 ]
	ml(2,:) = [ 2.d0, 4.d0, 1.d0, 2.d0 ] 
	ml(3,:) = [ 2.d0, 1.d0, 4.d0, 2.d0 ] 
	ml(4,:) = [ 1.d0, 2.d0, 2.d0, 4.d0 ]  

END SUBROUTINE matriz_Ml

!=================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=================================================================================================

COMPLEX(DP) FUNCTION sig_iw_cole( sig, omg, m, tal, ccf )
!=============================================================
! 
! COLE-COLE'S PARÂMETERS:
! 
! sig ---- True conductivity at the limit when x = infinity.
! 
! m ------- Is the polarization coefficient or chargeability.
! 
! tal ----- Time constant or relaxation time of the pair
!           electric layer of Helmholtz
! 
! ccf ----- It is the frequency correlation coefficient
!           commonly expressed by the letter 'C'.
! 
! f ------- Frequency.
! 
!=============================================================
IMPLICIT NONE 
REAL(dp), INTENT(in):: sig, m, tal, ccf, omg

REAL(dp) :: sig0, sig00

!========================================================================
!                 Se Sig representar a condutividade
!========================================================================
! 
	 sig0 = sig*(1-m)
 	 sig_iw_cole = sig0/( 1-m*( 1-(1/(1+(ci*omg*tal)**ccf)) ) )
! 
!========================================================================

END FUNCTION sig_iw_cole


!=================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=================================================================================================

COMPLEX(DP) FUNCTION sig_iw_Multcole( sig0, omg, m1, tal1, ccf1, m2, tal2, ccf2 )
!=============================================================
! 
! COLE-COLE'S PARÂMETERS:
! 
! sig0 --------- Condutividade na frequência zero corrente d.c
!	         quando w = 0.
! 
! m1,m2 -------- Are the polarization coefficient or chargeability.
! 
! tal1,tal2 ---- Time constant or relaxation time of the pair
! 
! 
! ccf1, ccf2 --- It is the frequency correlation coefficient
!                commonly expressed by the letter 'C'.
! 
! f ------------ Frequency.
! 
!=============================================================
IMPLICIT NONE 
REAL(dp), INTENT(in):: sig0, omg, m1, tal1, ccf1, m2, tal2, ccf2

COMPLEX(dp) :: rhow, rho0


	rhow =  ( 1.d0 - m1/(1.d0+(1.d0-m1)*(ci*omg*Tal1)**CCf1) ) * &
		( 1.d0 - m2/(1.d0+(1.d0-m2)*(ci*omg*Tal2)**CCf2) ) *sig0

	sig_iw_Multcole = rhow


END FUNCTION sig_iw_Multcole


!=================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=================================================================================================

COMPLEX(DP) FUNCTION sig_iw_Dias( sig0, omg, m, tal, det, eta )
!=============================================================
! 
! DIAS'S PARÂMETERS:
! 
! sig0 ---- Condutividade na frequência zero corrente d.c
!	    quando w = 0.
! 
! m ------- Is the polarization coefficient or chargeability.
! 
! tal ----- Time constant or relaxation time of the pair
!           electric layer of Helmholtz
! 
! det ----- Fraction pore length affected by polarization.
! 
! eta ----- Electrochemical Parameter.
! 
!=============================================================
IMPLICIT NONE 
REAL(dp), INTENT(in):: sig0, m, tal, det, eta, omg

REAL(dp)    :: alf, lmb, lmb2, bta
COMPLEX(dp) :: iomg, valaux, mu, sig0aux


        iomg   = ci*omg
	valaux = omg*tal*(1+(eta/cdsqrt(iomg)))
	mu   = ci*valaux
	lmb  = 1 + mu 
	lmb2 = 1 + (1-det)*mu
	alf  = m*(1-det)/(1-m)
	bta  = 1/(eta*det)
	

	sig0aux = sig0*(1-m)
	sig_iw_Dias = sig0aux*( 1 + alf*(lmb*bta*cdsqrt(iomg))/(1+lmb2*bta*cdsqrt(iomg)) )


END FUNCTION sig_iw_Dias

!=================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=================================================================================================

END MODULE GEMM3D
