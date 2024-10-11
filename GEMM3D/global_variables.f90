include 'mkl_pardiso.f90'
MODULE G_variables

	USE MKL_PARDISO
!================================================================================================================
!===========================       Declaration of global variables  	    =====================================
!================================================================================================================

	INTEGER    , PARAMETER	:: dp = kind(1.d0)
	COMPLEX(dp), PARAMETER	:: ci=(0.d0,1.d0)
	REAL(dp)   , PARAMETER	:: pi = 4.d0*datan(1.d0)
	REAL(dp)   , PARAMETER	:: mi0 = 4.d-7*pi, epson = 8.85d-12
	REAL(dp)   , PARAMETER	:: tolrel=10.d-6, tolabs=10.d-30, small = 1.d-15

	COMPLEX(dp), ALLOCATABLE, DIMENSION(:)	:: vfonte, val_nz, Es_1, Es_2

	COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:):: Es_obs_1, Hs_obs_1, Et_obs_1, Ht_obs_1, Et_obs_1aux, Et_obs_2aux
	COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:):: Es_obs_2, Hs_obs_2, Et_obs_2, Ht_obs_2, Ht_obs_1aux, Ht_obs_2aux

	COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:):: Es_obs_DP , Hs_obs_DP , Et_obs_DP , Ht_obs_DP
	COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:):: Es_obs_BCH, Hs_obs_BCH, Et_obs_BCH, Ht_obs_BCH

	COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:):: MSEx1, MSEx2, MSEy1, MSEy2
	COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:):: MSHx1, MSHx2, MSHy1, MSHy2

	COMPLEX(dp), ALLOCATABLE, DIMENSION(:)	:: Es_adj_1, Hs_adj_1, Et_adj_1, Ht_adj_1
	COMPLEX(dp), ALLOCATABLE, DIMENSION(:)	:: Es_adj_2, Hs_adj_2, Et_adj_2, Ht_adj_2
	COMPLEX(dp), ALLOCATABLE, DIMENSION(:)	:: Es_adj_3, Hs_adj_3, Et_adj_3, Ht_adj_3
	COMPLEX(dp), ALLOCATABLE, DIMENSION(:)	:: Es_adj_4, Hs_adj_4, Et_adj_4, Ht_adj_4

	COMPLEX(dp), ALLOCATABLE, DIMENSION(:)	:: Exp_obs_adj, Eyp_obs_adj, Ezp_obs_adj
	COMPLEX(dp), ALLOCATABLE, DIMENSION(:)	:: Exp_trans , Eyp_trans , Ezp_trans
	COMPLEX(dp), ALLOCATABLE, DIMENSION(:)	:: Exp_trans1, Eyp_trans1, Ezp_trans1
	COMPLEX(dp), ALLOCATABLE, DIMENSION(:)	:: Exp_trans2, Eyp_trans2, Ezp_trans2
	COMPLEX(dp), ALLOCATABLE, DIMENSION(:)	:: Exp_adj, Eyp_adj, Ezp_adj



	REAL(dp), ALLOCATABLE, DIMENSION(:,:)	:: Mnos
	REAL(dp), ALLOCATABLE, DIMENSION(:)	:: hj, roj, mirj, rojc, zkmadaj
	REAL(dp), ALLOCATABLE, DIMENSION(:)	:: vetf, cxT, cyT, czT, zrs_J0, zrs_J1, xGssLg, wGssLg

	REAL(dp), ALLOCATABLE, DIMENSION(:) 	:: rho_ap_xx, rho_ap_yy, rho_ap_xy, rho_ap_yx
	REAL(dp), ALLOCATABLE, DIMENSION(:) 	:: fase_xx, fase_yy, fase_xy, fase_yx, vrho_e

	REAL(dp), ALLOCATABLE, DIMENSION(:) 	:: rho_ap_xy_aux, rho_ap_yx_aux
	REAL(dp), ALLOCATABLE, DIMENSION(:) 	:: fase_xy_aux, fase_yx_aux

	REAL(dp), ALLOCATABLE, DIMENSION(:) 	:: CARGBLD, TMRLX, FCPIP, PEQMC, CCFREQ, CARGBLD2, TMRLX2, CCFREQ2 

	INTEGER, ALLOCATABLE, DIMENSION(:,:)	:: Marestas, Mnosaresta, Melem
	INTEGER, ALLOCATABLE, DIMENSION(:)	:: row_ind, col_ptr, vbordas, IPELE
	INTEGER, ALLOCATABLE, DIMENSION(:)	:: nelebloc, elebloc, eleobs, indobs, no_obs

	REAL(dp)				:: ISOURCE, SDIMENS

	INTEGER					:: Nnode, NC, nfreq, nnosmalha, Nelem, Narestas, nborda
	INTEGER					:: Ncorpos, ngl, ny, nf, NP, Nhet, INPUTIP
	INTEGER					:: Nobs, Nbloc, size_eblc, size_eobs, Ntransm
	INTEGER					:: nnz, nnz_aux, autr_fltJ0, numpJ0, autr_fltJ1, numpJ1

	CHARACTER(LEN=300)                      :: pathr, MT_CSAMT		

!======================================  Input variables for pardiso calls  =====================================

	TYPE(MKL_PARDISO_HANDLE), ALLOCATABLE, DIMENSION(:)  :: pt
	INTEGER, ALLOCATABLE, DIMENSION(:)                   :: perm

	INTEGER		:: nrhs, mtype, maxfct
	INTEGER		:: mnum, msglvl, iparm(64)

!================================================================================================================

END MODULE G_variables
