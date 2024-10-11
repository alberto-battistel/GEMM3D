PROGRAM MAIN_ME
!=======================================================================================

 USE GEMM3D

!=======================================================================================
IMPLICIT NONE
 REAL(dp), ALLOCATABLE, DIMENSION(:,:)	:: MSB
 REAL(dp), ALLOCATABLE, DIMENSION(:)	:: Ax
 INTEGER				:: IERR, etapa,  nite, nL, k, i, j
 INTEGER				:: adj1D, Trans1D, nproc, nlin, locadj
 
 CHARACTER(LEN=20)			:: DPT, DPR, modlagm 
 CHARACTER(LEN=300)			:: path_atual

	!Global variables that you will use to allocate the sensitivity matrix 
        !and the estimated data vector:

	!Nobs     --> Integer that has the number of observations
	!Nfreq    --> Integer that has the number of frequencies
	!Ntransm  --> Integer that has the number of transmitters
	!Nbloc	  --> Integer that has the number of Inversion cells (or Heterogeneity)

	! SECTION 1: Starts input and the filter for the 1D modeling of the dipoles.	

		! Entering the path/GEMM3D/
		 path_atual = '~/Desktop/GEMM3D/'	
		 call  inic_entrada( path_atual )

	! SECTION 2: If you want to change the filter.
		
		!Exemple: Changing to the 101-point Karry Key filter
		!autr_fltJ0 = 4
		!numpJ0     = 101
		!autr_fltJ1 = 4 
		!numpJ1     = 101
	 

	! SECTION 3: Examples to simulate the direct models of MT, CSAMT, Circular Coil (Coplanar Array) and dipoles.

                !===============================================================================
		! Ex( MT )
		! Routine that simulates MT 3D results and stores in the folder /MT_fields

                      ! call MODELO_SINTET_MT
                      ! call MODELO_SINTET_MT_CSIG

                !===============================================================================
		! Ex( CSAMT )
		! Routine that simulates CSAMT 3D results and stores in the folder /CSAMT_fields

                      ! adj1D = 1
                      ! nproc = 6
                      ! call MODELO_SINTET_CSAMT( adj1D, nproc )
                      ! call MODELO_SINTET_CSAMT_CSIG( adj1D, nproc )

                !===============================================================================
		! Ex( DPs )
		!Routine that simulates the fields of the dipoles, in this case, the Electric 
                !in the x-direction and stores in the folder / DPEM / DPEx fields

                      ! DPT = 'ExDP'
                      ! adj1D = 1
                      ! nproc = 6
                      ! call MODELO_SINTET_DPEH( DPT, adj1D, nproc )
                      ! call MODELO_SINTET_DPEH_CSIG( DPT, adj1D, n;proc )

                !===============================================================================
		! Ex( HCL )
		! Routine that simulates 3D HCL results and stores in folder / HCL_fields.

                !      Trans1D = 1
                !      nproc   = 6
                !      call MODELO_SINTET_HCL( Trans1D, nproc )
                !      call MODELO_SINTET_HCL_CSIG( Trans1D, nproc )

                !===============================================================================


	! SECTION 4: Examples for the Jacobian calculation( or sensitivity matrix ).

                !===============================================================================
		modlagm = 'MT'     ! If MT uncomment and comment the other two
!		modlagm = 'CSAMT'  ! If CSAMT uncomment and comment the other two
!		modlagm = 'DPEM'   ! If DPEM uncomment and comment the other two



		!Defining the number of lines of the sensitivity matrix for MT and dipoles
		if( modlagm == 'MT' .or. modlagm == 'CSAMT' )then	
	
			nlin = 4*Nobs*Nfreq
	
		elseif( modlagm == 'DPEM' )then
	
			nlin = 2*Nobs*Nfreq*Ntransm
		end if
	

		IF( modlagm == 'MT' .or.  modlagm == 'CSAMT' .or. modlagm == 'DPEM' )THEN
	
			!Allocating the sensitivity matrix for MT or dipoles
			ALLOCATE( MSB(nlin,Nbloc), STAT = IERR )
				IF ( IERR == 0 )THEN
				MSB = 0.d0
			ELSE IF ( IERR /= 0 ) THEN
				WRITE(*,*)'ERROR IN THE ALLOCATION OF THE SENSITIVITY MATRIX (MSB)'
				STOP
			END IF
			ALLOCATE( Ax(nlin), STAT = IERR )
				IF ( IERR == 0 )THEN
				Ax = 0.d0
			ELSE IF ( IERR /= 0 ) THEN
				WRITE(*,*)'ERROR IN THE ALLOCATION OF THE VECTOR (Ax) ESTIMATED '
				STOP
			END IF

		END IF
	
		20 format( 2x, 2d25.16 )


                !===============================================================================
                !                          calculating the sensitivity matrix
                !===============================================================================

                !===============================================================================
		!Example - MT: 
                !===============================================================================

		if( modlagm == 'MT' )then	
	
	
			call open_arq_sensib_MT

			locadj  = 0  ! informs whether or not the program will read the adjunct sources already calculated and  				     ! stored from a previous CSAMT simulation for the same model. If / = 0 YES, if = 0 NO. 	                                 ! Therefore, note that if locadj == 1, adj1D must be started by 0.
                                     ! But if adj1D = 1, the value of locadj will be ignored.
			adj1D   = 1
			nproc   = 6
	
			nite = 1
			do i=1,nite
	
				if( i > 1 )then
					adj1D   = 0
				end if			
	
				call JACOB_MT3D( locadj, adj1D, nproc, MSB, Ax )
				!call JACOB_MT3D_CSIG( locadj, adj1D, nproc, MSB, Ax )
			
				nL = 2*Nobs*Nfreq
				do k=1, nL		
					do j=1, Nbloc
						write(3600,20) MSB(k,j)
						write(3700,20) MSB(nL+k,j) 
					end do
				end do
			end do
			
		end if

                !===============================================================================
		!Example - CSAMT:
                !===============================================================================


		if( modlagm == 'CSAMT' )then	

			call open_arq_sensib_CSAMT

			Trans1D = 1
			adj1D   = 1
			locadj  = 0  ! informs whether or not the program will read the adjunct sources already calculated and   					     ! stored from a previous MT simulation for the same model. If / = 0 YES, if = 0 NO.
                                     !  Therefore, note that if locadj == 1, adj1D must be started by 0.
                                     ! But if adj1D = 1, the value of locadj will be ignored.
			nproc   = 6
			nite = 1

			do i=1,nite
	
				if( i > 1 )then
					adj1D   = 0
				end if			

				call JACOB_CSAMT3D( locadj, Trans1D, adj1D, nproc, MSB, Ax )
				!call JACOB_CSAMT3D_CSIG( locadj, Trans1D, adj1D, nproc, MSB, Ax )

				nL = 2*Nobs*Nfreq
				do k=1, nL		
					do j=1, Nbloc
						write(3600,20) MSB(k,j)
						write(3700,20) MSB(nL+k,j) 
					end do
				end do
			end do

		end if

                !===============================================================================
		!Example - DPs 
                !===============================================================================

		if( modlagm == 'DPEM' )then
	
	
			DPT = 'ExDP'
			DPR = 'ExDP'
			Trans1D = 1
			adj1D   = 1
			nproc   = 6
			
			call open_arq_sensib_DPEM( DPT )
	
			nite = 1
			do i=1,nite
	
				if( i > 1 )then
					Trans1D = 0
					adj1D   = 0
				end if			
				
				call JACOB_DPEM3D( DPT, DPR, adj1D, Trans1D, nproc, MSB, Ax )
				!call JACOB_DPEM3D_CSIG( DPT, DPR, adj1D, Trans1D, nproc, MSB, Ax )
		
				nL = Nobs*Nfreq*Ntransm 

				 do k=1, nL		
					do j=1, Nbloc
						write(3610,20) MSB(k,j)
						write(3710,20) MSB(nL+k,j)
					end do
				
				 end do	
			end do	
	
		end if

                !===============================================================================
                !                       Deallocate: MSB e Ax 
                !===============================================================================

		IF( modlagm == 'MT' .or. modlagm == 'DPEM' .or. modlagm == 'CSAMT' )THEN
	
			! MT, CSAMT or dipoles
			DEALLOCATE( MSB, STAT = IERR )
			IF ( IERR /= 0 ) THEN
				WRITE(*,*)'MEMORY RELEASE ERROR OF THE SENSITIVITY MATRIX (MSB).'
				WRITE(*,*)'ERROR =',IERR
				STOP
			END IF
			DEALLOCATE( Ax, STAT = IERR )
			IF ( IERR /= 0 ) THEN
				WRITE(*,*)'MEMORY RELEASE ERROR OF THE VECTOR (Ax).'
				WRITE(*,*)'ERROR =',IERR
				STOP
			END IF
		END IF
                !===============================================================================
                !===============================================================================


END PROGRAM MAIN_ME
