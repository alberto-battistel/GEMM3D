MODULE DADOS_SAIDA

USE G_variables
USE ler_escrever_arq
USE mdlgm1_1D
USE mdlgm2_1D

CONTAINS

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE DPEH1D_ADJ_MT( nprcss )
IMPLICIT NONE	
 INTEGER, INTENT(in)			:: nprcss
 REAL(dp)   , ALLOCATABLE, DIMENSION(:) :: sigmas
 COMPLEX(dp), ALLOCATABLE, DIMENSION(:) :: Ex_dpex, Ey_dpex, Ez_dpex
 COMPLEX(dp), ALLOCATABLE, DIMENSION(:) :: Ex_dpey, Ey_dpey, Ez_dpey
 COMPLEX(dp), ALLOCATABLE, DIMENSION(:) :: Ex_dpez, Ey_dpez, Ez_dpez

 COMPLEX(dp), ALLOCATABLE, DIMENSION(:) :: Ex_dphx, Ey_dphx, Ez_dphx
 COMPLEX(dp), ALLOCATABLE, DIMENSION(:) :: Ex_dphy, Ey_dphy, Ez_dphy
 COMPLEX(dp), ALLOCATABLE, DIMENSION(:) :: Ex_dphz, Ey_dphz, Ez_dphz

 COMPLEX(dp)	:: Hx, Hy, Hz
 COMPLEX(dp)	:: Ex_aux, Ey_aux, Ez_aux, z_aux, neta0, zeta
 REAL(dp)	:: w, Iw, dsx, dsy, dsz, Tx, Ty, Tz, t1, t2
 REAL(dp)	:: x(8), y(8), z(8)
 INTEGER	:: i, j, k, p, l, inde, nosG_e(8), areG_e(12), estts(6)
 INTEGER	:: num, num2, id, OMP_GET_THREAD_NUM, IERR 

 CHARACTER(len=50)                    :: nomevar

call openRW_arq_1D_adj_MT

allocate( sigmas(NC))

sigmas = 1.d0/roj

Iw  = ISOURCE
dsx = SDIMENS
dsy = SDIMENS
dsz = SDIMENS

neta0 = 1.d-12  !(0.d0,1.d0)*omega*eps    !

!======================================================================================================
!					 Primário do Receptor
!======================================================================================================

num  = (nelebloc(Nbloc+1)-1)*8
num2 = num
num = 0

nomevar = 'Ex_dpex field vector'
call alloc_cdvet( Ex_dpex, num2, nomevar )
nomevar = 'Ey_dpex field vector '
call alloc_cdvet( Ey_dpex, num2, nomevar )
nomevar = 'Ez_dpex field vector '
call alloc_cdvet( Ez_dpex, num2, nomevar )

nomevar = 'Ex_dpey field vector'
call alloc_cdvet( Ex_dpey, num2, nomevar )
nomevar = 'Ey_dpey field vector'
call alloc_cdvet( Ey_dpey, num2, nomevar )
nomevar = 'Ez_dpey field vector'
call alloc_cdvet( Ez_dpey, num2, nomevar )

nomevar = 'Ex_dphx field vector'
call alloc_cdvet( Ex_dphx, num2, nomevar )
nomevar = 'Ey_dphx field vector'
call alloc_cdvet( Ey_dphx, num2, nomevar )
nomevar = 'Ez_dphx field vector'
call alloc_cdvet( Ez_dphx, num2, nomevar )

nomevar = 'Ex_dphy field vector'
call alloc_cdvet( Ex_dphy, num2, nomevar )
nomevar = 'Ey_dphy field vector'
call alloc_cdvet( Ey_dphy, num2, nomevar )
nomevar = 'Ez_dphy field vector'
call alloc_cdvet( Ez_dphy, num2, nomevar )


 do i=1, Nfreq

	w = (pi+pi)*vetf(i)
	zeta = (0.d0,1.d0)*w*mi0

	do j=1, Nobs

		Tx  = Mnos(no_obs(j),1)
		Ty  = Mnos(no_obs(j),2)
		Tz  = Mnos(no_obs(j),3)
		
		do l=nelebloc(1), nelebloc(Nbloc+1)-1

			inde = elebloc(l)
			nosG_e(:) = Melem(inde,:)
			areG_e(:) = Marestas(inde,:)
			x(:) = Mnos(nosG_e(:),1)
			y(:) = Mnos(nosG_e(:),2)
			z(:) = Mnos(nosG_e(:),3)

			CALL OMP_SET_NUM_THREADS( nprcss )

			!==================================================================================
			!==================================================================================

		        !$OMP PARALLEL PRIVATE(p, num) 
			!$OMP DO
			do p=1, 8

				num  = (l-1)*8 + p 				
				call dehx_xyz_loops( autr_fltJ0, numpJ0, autr_fltJ1, numpJ1, Tx, Ty, Iw, dsx, Tz, &
						     NC, hj, sigmas, neta0, zeta, x(p), y(p), z(p), Ex_dpex(num), &
									   Ey_dpex(num), Ez_dpex(num), Hx, Hy, Hz )
			end do
		        !$OMP END  DO
			!$OMP END PARALLEL

			!==================================================================================
			!==================================================================================

		        !$OMP PARALLEL PRIVATE(p, num)
			!$OMP DO
			do p=1, 8

				num  = (l-1)*8 + p 
				call dehy_xyz_loops(  autr_fltJ0, numpJ0, autr_fltJ1, numpJ1, Tx, Ty, Iw, dsy, Tz, &
						      NC, hj, sigmas, neta0, zeta, x(p), y(p), z(p), Ex_dpey(num), &
									    Ey_dpey(num), Ez_dpey(num), Hx, Hy, Hz )
			end do
		        !$OMP END DO
			!$OMP END PARALLEL

			!==================================================================================
			!==================================================================================

		        !$OMP PARALLEL PRIVATE(p, num)
			!$OMP DO
			do p=1, 8

				num  = (l-1)*8 + p 
				call dmhx_xyz_loops( autr_fltJ0, numpJ0, autr_fltJ1, numpJ1, Tx, Ty, Iw, dsx, Tz, &
						     NC, hj, sigmas, neta0, zeta, x(p), y(p), z(p), Ex_dphx(num), &
									   Ey_dphx(num), Ez_dphx(num), Hx, Hy, Hz )
			end do
		        !$OMP END  DO
			!$OMP END PARALLEL

			!==================================================================================
			!==================================================================================

		        !$OMP PARALLEL PRIVATE(p, num)
			!$OMP DO
			do p=1, 8

				num  = (l-1)*8 + p 
				call dmhx_xyz_loops( autr_fltJ0, numpJ0, autr_fltJ1, numpJ1, Ty, -Tx, Iw, dsx, Tz, &
						     NC, hj, sigmas, neta0, zeta, y(p), -x(p), z(p), Ex_dphy(num), &
									    Ey_dphy(num), Ez_dphy(num), Hx, Hy, Hz )
			end do
		        !$OMP END DO
			!$OMP END PARALLEL

			!==================================================================================
			!==================================================================================

		end do

 
		call armaznr_1Dadj_MT( num2, zeta, Ex_dpex, Ey_dpex, Ez_dpex, Ex_dpey, Ey_dpey, Ez_dpey, &
					           Ex_dphx, Ey_dphx, Ez_dphx, Ex_dphy, Ey_dphy, Ez_dphy)
 
		Ex_dpex=(0.d0,0.d0); Ey_dpex=(0.d0,0.d0); Ez_dpex=(0.d0,0.d0)
		Ex_dpey=(0.d0,0.d0); Ey_dpey=(0.d0,0.d0); Ez_dpey=(0.d0,0.d0)
		
		Ex_dphx=(0.d0,0.d0); Ey_dphx=(0.d0,0.d0); Ez_dphx=(0.d0,0.d0)
		Ex_dphy=(0.d0,0.d0); Ey_dphy=(0.d0,0.d0); Ez_dphy=(0.d0,0.d0)

	end do	
 end do

 call close_arq_1D_adj_MT
 call openRW_arq_1D_adj_MT

deallocate( Ex_dpex, Ey_dpex, Ez_dpex )
deallocate( Ex_dpey, Ey_dpey, Ez_dpey )
deallocate( Ex_dphx, Ey_dphx, Ez_dphx )
deallocate( Ex_dphy, Ey_dphy, Ez_dphy )

END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

subroutine armaznr_1Dadj_MT( num, zeta, Ex_dpex, Ey_dpex, Ez_dpex, Ex_dpey, Ey_dpey, Ez_dpey, &
				        Ex_dphx, Ey_dphx, Ez_dphx, Ex_dphy, Ey_dphy, Ez_dphy)

 COMPLEX(dp), DIMENSION(:), INTENT(in) :: Ex_dpex, Ey_dpex, Ez_dpex, Ex_dpey, Ey_dpey, Ez_dpey
 COMPLEX(dp), DIMENSION(:), INTENT(in) :: Ex_dphx, Ey_dphx, Ez_dphx, Ex_dphy, Ey_dphy, Ez_dphy
 COMPLEX(dp)		  , INTENT(in) :: zeta
 INTEGER		  , INTENT(in) :: num

 INTEGER :: estts(6)

		write(unit=3000, IOSTAT=estts(1)) (dreal(Ex_dpex(l)),l=1,num)
		write(unit=3000, IOSTAT=estts(2)) (aimag(Ex_dpex(l)),l=1,num)
		write(unit=3000, IOSTAT=estts(3)) (dreal(Ey_dpex(l)),l=1,num)
		write(unit=3000, IOSTAT=estts(4)) (aimag(Ey_dpex(l)),l=1,num)
		write(unit=3000, IOSTAT=estts(5)) (dreal(Ez_dpex(l)),l=1,num)
		write(unit=3000, IOSTAT=estts(6)) (aimag(Ez_dpex(l)),l=1,num)
		IF( estts(1) /= 0 )THEN
			PRINT*,'An error occurred while writing to the 3000 drive file'
			PRINT*,'when writing the real part of the vector Ex_dpex'
		END IF
		IF( estts(2) /= 0 )THEN
			PRINT*,'An error occurred while writing to the 3000 drive file'
			PRINT*,'when writing the imaginary part of the vector Ex_dpex' 	 
		END IF
		IF( estts(3) /= 0 )THEN
			PRINT*,'An error occurred while writing to the 3000 drive file'
			PRINT*,'when writing the real part of the vector Ey_dpex'
		END IF
		IF( estts(4) /= 0 )THEN
			PRINT*,'An error occurred while writing to the 3000 drive file'
			PRINT*,'when writing the imaginary part of the vector Ey_dpex' 	 
		END IF
		IF( estts(5) /= 0 )THEN
			PRINT*,'An error occurred while writing to the 3000 drive file'
			PRINT*,'when writing the real part of the vector Ez_dpex'
		END IF
		IF( estts(6) /= 0 )THEN
			PRINT*,'An error occurred while writing to the 3000 drive file'
			PRINT*,'when writing the imaginary part of the vector Ez_dpex' 	 
		END IF


		write(unit=3100, IOSTAT=estts(1)) (dreal(Ex_dpey(l)),l=1,num)
		write(unit=3100, IOSTAT=estts(2)) (aimag(Ex_dpey(l)),l=1,num)
		write(unit=3100, IOSTAT=estts(3)) (dreal(Ey_dpey(l)),l=1,num)
		write(unit=3100, IOSTAT=estts(4)) (aimag(Ey_dpey(l)),l=1,num)
		write(unit=3100, IOSTAT=estts(5)) (dreal(Ez_dpey(l)),l=1,num)
		write(unit=3100, IOSTAT=estts(6)) (aimag(Ez_dpey(l)),l=1,num)
		IF( estts(1) /= 0 )THEN
			PRINT*,'An error occurred while writing to the 3100 drive file'
			PRINT*,'when writing the real part of the vector Ex_dpey'
		END IF
		IF( estts(2) /= 0 )THEN
			PRINT*,'An error occurred while writing to the 3100 drive file'
			PRINT*,'when writing the imaginary part of the vector Ex_dpey' 	 
		END IF
		IF( estts(3) /= 0 )THEN
			PRINT*,'An error occurred while writing to the 3100 drive file'
			PRINT*,'when writing the real part of the vector Ey_dpey'
		END IF
		IF( estts(4) /= 0 )THEN
			PRINT*,'An error occurred while writing to the 3100 drive file'
			PRINT*,'when writing the imaginary part of the vector Ey_dpey' 	 
		END IF
		IF( estts(5) /= 0 )THEN
			PRINT*,'An error occurred while writing to the 3100 drive file'
			PRINT*,'when writing the real part of the vector Ez_dpey'
		END IF
		IF( estts(6) /= 0 )THEN
			PRINT*,'An error occurred while writing to the 3100 drive file'
			PRINT*,'when writing the imaginary part of the vector Ez_dpey' 	 
		END IF


		write(unit=3300, IOSTAT=estts(1)) (dreal(Ex_dphx(l)),l=1,num)
		write(unit=3300, IOSTAT=estts(2)) (aimag(Ex_dphx(l)),l=1,num)
		write(unit=3300, IOSTAT=estts(3)) (dreal(Ey_dphx(l)),l=1,num)
		write(unit=3300, IOSTAT=estts(4)) (aimag(Ey_dphx(l)),l=1,num)
		write(unit=3300, IOSTAT=estts(5)) (dreal(Ez_dphx(l)),l=1,num)
		write(unit=3300, IOSTAT=estts(6)) (aimag(Ez_dphx(l)),l=1,num)
		IF( estts(1) /= 0 )THEN
			PRINT*,'An error occurred while writing to the 3300 drive file'
			PRINT*,'when writing the real part of the vector Ex_dphx'
		END IF
		IF( estts(2) /= 0 )THEN
			PRINT*,'An error occurred while writing to the 3300 drive file'
			PRINT*,'when writing the imaginary part of the vector Ex_dphx' 	 
		END IF
		IF( estts(3) /= 0 )THEN
			PRINT*,'An error occurred while writing to the 3300 drive file'
			PRINT*,'when writing the real part of the vector Ey_dphx'
		END IF
		IF( estts(4) /= 0 )THEN
			PRINT*,'An error occurred while writing to the 3300 drive file'
			PRINT*,'when writing the imaginary part of the vector Ey_dphx' 	 
		END IF
		IF( estts(5) /= 0 )THEN
			PRINT*,'An error occurred while writing to the 3300 drive file'
			PRINT*,'when writing the real part of the vector Ez_dphx'
		END IF
		IF( estts(6) /= 0 )THEN
			PRINT*,'An error occurred while writing to the 3300 drive file'
			PRINT*,'when writing the imaginary part of the vector Ez_dphx' 	 
		END IF


		write(unit=3400, IOSTAT=estts(3)) (dreal(-Ey_dphy(l)),l=1,num)
		write(unit=3400, IOSTAT=estts(4)) (aimag(-Ey_dphy(l)),l=1,num)
		write(unit=3400, IOSTAT=estts(1)) (dreal( Ex_dphy(l)),l=1,num)
		write(unit=3400, IOSTAT=estts(2)) (aimag( Ex_dphy(l)),l=1,num)
		write(unit=3400, IOSTAT=estts(5)) (dreal( Ez_dphy(l)),l=1,num)
		write(unit=3400, IOSTAT=estts(6)) (aimag( Ez_dphy(l)),l=1,num)
		IF( estts(1) /= 0 )THEN
			PRINT*,'An error occurred while writing to the 3400 drive file'
			PRINT*,'when writing the real part of the vector Ex_dphy'
		END IF
		IF( estts(2) /= 0 )THEN
			PRINT*,'An error occurred while writing to the 3400 drive file'
			PRINT*,'when writing the imaginary part of the vector Ex_dphy' 	 
		END IF
		IF( estts(3) /= 0 )THEN
			PRINT*,'An error occurred while writing to the 3400 drive file'
			PRINT*,'when writing the real part of the vector Ey_dphy'
		END IF
		IF( estts(4) /= 0 )THEN
			PRINT*,'An error occurred while writing to the 3400 drive file'
			PRINT*,'when writing the imaginary part of the vector Ey_dphy' 	 
		END IF
		IF( estts(5) /= 0 )THEN
			PRINT*,'An error occurred while writing to the 3400 drive file'
			PRINT*,'when writing the real part of the vector Ez_dphy'
		END IF
		IF( estts(6) /= 0 )THEN
			PRINT*,'An error occurred while writing to the 3400 drive file'
			PRINT*,'when writing the imaginary part of the vector Ez_dphy' 	 
		END IF

end subroutine

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE DPEH1D_ADJ_CSAMT( nprcss )
IMPLICIT NONE	
 INTEGER, INTENT(in)			:: nprcss
 REAL(dp)   , ALLOCATABLE, DIMENSION(:) :: sigmas
 COMPLEX(dp), ALLOCATABLE, DIMENSION(:) :: Ex_dpex, Ey_dpex, Ez_dpex
 COMPLEX(dp), ALLOCATABLE, DIMENSION(:) :: Ex_dpey, Ey_dpey, Ez_dpey
 COMPLEX(dp), ALLOCATABLE, DIMENSION(:) :: Ex_dpez, Ey_dpez, Ez_dpez

 COMPLEX(dp), ALLOCATABLE, DIMENSION(:) :: Ex_dphx, Ey_dphx, Ez_dphx
 COMPLEX(dp), ALLOCATABLE, DIMENSION(:) :: Ex_dphy, Ey_dphy, Ez_dphy
 COMPLEX(dp), ALLOCATABLE, DIMENSION(:) :: Ex_dphz, Ey_dphz, Ez_dphz

 COMPLEX(dp)	   :: Hx, Hy, Hz
 COMPLEX(dp)	   :: Ex_aux, Ey_aux, Ez_aux, z_aux, neta0, zeta
 REAL(dp)	   :: w, Iw, dsx, dsy, dsz, Tx, Ty, Tz, t1, t2
 REAL(dp)	   :: x(8), y(8), z(8)
 INTEGER	   :: i, j, k, p, l, inde, nosG_e(8), areG_e(12)
 INTEGER	   :: num, num2, id, OMP_GET_THREAD_NUM 

 CHARACTER(len=50) :: nomevar

call openRW_arq_1D_adj_CSAMT

allocate( sigmas(NC))

sigmas = 1.d0/roj

Iw  = ISOURCE
dsx = SDIMENS
dsy = SDIMENS
dsz = SDIMENS

neta0 = 1.d-12  !(0.d0,1.d0)*omega*eps    !

!======================================================================================================
!					 Primário do Receptor
!======================================================================================================

num = (nelebloc(Nbloc+1)-1)*8
num2 = num
num = 0

nomevar = 'Ex_dpex field vector'
call alloc_cdvet( Ex_dpex, num2, nomevar )
nomevar = 'Ey_dpex field vector '
call alloc_cdvet( Ey_dpex, num2, nomevar )
nomevar = 'Ez_dpex field vector '
call alloc_cdvet( Ez_dpex, num2, nomevar )

nomevar = 'Ex_dpey field vector'
call alloc_cdvet( Ex_dpey, num2, nomevar )
nomevar = 'Ey_dpey field vector'
call alloc_cdvet( Ey_dpey, num2, nomevar )
nomevar = 'Ez_dpey field vector'
call alloc_cdvet( Ez_dpey, num2, nomevar )

nomevar = 'Ex_dphx field vector'
call alloc_cdvet( Ex_dphx, num2, nomevar )
nomevar = 'Ey_dphx field vector'
call alloc_cdvet( Ey_dphx, num2, nomevar )
nomevar = 'Ez_dphx field vector'
call alloc_cdvet( Ez_dphx, num2, nomevar )

nomevar = 'Ex_dphy field vector'
call alloc_cdvet( Ex_dphy, num2, nomevar )
nomevar = 'Ey_dphy field vector'
call alloc_cdvet( Ey_dphy, num2, nomevar )
nomevar = 'Ez_dphy field vector'
call alloc_cdvet( Ez_dphy, num2, nomevar )


 do i=1, Nfreq

	w = (pi+pi)*vetf(i)
	zeta = (0.d0,1.d0)*w*mi0

	do j=1, Nobs

		Tx  = Mnos(no_obs(j),1)
		Ty  = Mnos(no_obs(j),2)
		Tz  = Mnos(no_obs(j),3)
		
		do l=nelebloc(1), nelebloc(Nbloc+1)-1

			inde = elebloc(l)
			nosG_e(:) = Melem(inde,:)
			areG_e(:) = Marestas(inde,:)
			x(:) = Mnos(nosG_e(:),1)
			y(:) = Mnos(nosG_e(:),2)
			z(:) = Mnos(nosG_e(:),3)

			CALL OMP_SET_NUM_THREADS( nprcss )

			!==================================================================================
			!==================================================================================

			!$OMP PARALLEL PRIVATE(p, num) 
			!$OMP DO
			do p=1, 8

				num  = (l-1)*8 + p 				
				call dehx_xyz_loops( autr_fltJ0, numpJ0, autr_fltJ1, numpJ1, Tx, Ty, Iw, dsx, Tz, &
						     NC, hj, sigmas, neta0, zeta, x(p), y(p), z(p), Ex_dpex(num), &
									   Ey_dpex(num), Ez_dpex(num), Hx, Hy, Hz )
			end do
			!$OMP END PARALLEL

			!==================================================================================
			!==================================================================================

			!$OMP PARALLEL PRIVATE(p, num) 
			!$OMP DO
			do p=1, 8

				num  = (l-1)*8 + p 
				call dehy_xyz_loops(  autr_fltJ0, numpJ0, autr_fltJ1, numpJ1, Tx, Ty, Iw, dsy, Tz, &
						      NC, hj, sigmas, neta0, zeta, x(p), y(p), z(p), Ex_dpey(num), &
									    Ey_dpey(num), Ez_dpey(num), Hx, Hy, Hz )
			end do
                        !$OMP END  DO
			!$OMP END PARALLEL

			!==================================================================================
			!==================================================================================

			!$OMP PARALLEL PRIVATE(p, num) 
			!$OMP DO
			do p=1, 8

				num  = (l-1)*8 + p 
				call dmhx_xyz_loops( autr_fltJ0, numpJ0, autr_fltJ1, numpJ1, Tx, Ty, Iw, dsx, Tz, &
						     NC, hj, sigmas, neta0, zeta, x(p), y(p), z(p), Ex_dphx(num), &
									   Ey_dphx(num), Ez_dphx(num), Hx, Hy, Hz )
			end do
                        !$OMP END  DO
			!$OMP END PARALLEL

			!==================================================================================
			!==================================================================================

			!$OMP PARALLEL PRIVATE(p, num) 
			!$OMP DO
			do p=1, 8

				num  = (l-1)*8 + p 
				call dmhx_xyz_loops( autr_fltJ0, numpJ0, autr_fltJ1, numpJ1, Ty, -Tx, Iw, dsx, Tz, &
						     NC, hj, sigmas, neta0, zeta, y(p), -x(p), z(p), Ex_dphy(num), &
									    Ey_dphy(num), Ez_dphy(num), Hx, Hy, Hz )
			end do
                        !$OMP END  DO
			!$OMP END PARALLEL

			!==================================================================================
			!==================================================================================

		end do

		call armaznr_1Dadj_MT( num2, zeta, Ex_dpex, Ey_dpex, Ez_dpex, Ex_dpey, Ey_dpey, Ez_dpey, &
					           Ex_dphx, Ey_dphx, Ez_dphx, Ex_dphy, Ey_dphy, Ez_dphy)

		Ex_dpex=(0.d0,0.d0); Ey_dpex=(0.d0,0.d0); Ez_dpex=(0.d0,0.d0)
		Ex_dpey=(0.d0,0.d0); Ey_dpey=(0.d0,0.d0); Ez_dpey=(0.d0,0.d0)
		
		Ex_dphx=(0.d0,0.d0); Ey_dphx=(0.d0,0.d0); Ez_dphx=(0.d0,0.d0)
		Ex_dphy=(0.d0,0.d0); Ey_dphy=(0.d0,0.d0); Ez_dphy=(0.d0,0.d0)

	end do	
 end do

 call close_arq_1D_adj_CSAMT
 call openRW_arq_1D_adj_CSAMT

deallocate( Ex_dpex, Ey_dpex, Ez_dpex )
deallocate( Ex_dpey, Ey_dpey, Ez_dpey )
deallocate( Ex_dphx, Ey_dphx, Ez_dphx )
deallocate( Ex_dphy, Ey_dphy, Ez_dphy )

END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE DPEH1D_adj( DPEH, nprcss )
IMPLICIT NONE	
 INTEGER, INTENT(in)			:: nprcss	
 CHARACTER(LEN=20), INTENT(in)	:: DPEH

 REAL(dp)   , ALLOCATABLE, DIMENSION(:) :: sigmas
 COMPLEX(dp), ALLOCATABLE, DIMENSION(:) :: Exp_dp, Eyp_dp, Ezp_dp

 COMPLEX(dp)	:: Hx, Hy, Hz
 COMPLEX(dp)	:: Ex_aux, Ey_aux, Ez_aux, z_aux, neta0, zeta
 REAL(dp)	:: w, Iw, dsx, dsy, dsz, Tx, Ty, Tz, t1, t2
 REAL(dp)	:: x(8), y(8), z(8)
 INTEGER	:: i, j, k, p, l, inde, nosG_e(8), areG_e(12)
 INTEGER	:: num, num2, id, OMP_GET_THREAD_NUM , estts(6)

 CHARACTER(len=50) :: nomevar

call openRW_arq_DPEM1D_adj( DPEH )

allocate( sigmas(NC))

sigmas = 1.d0/roj

Iw  = ISOURCE
dsx = SDIMENS
dsy = SDIMENS
dsz = SDIMENS

neta0 = 1.d-12!  (0.d0,1.d0)*omega*eps    !

num  = (nelebloc(Nbloc+1)-1)*8
num2 = num
num  = 0

nomevar = 'field vector Exp_dp'
call alloc_cdvet( Exp_dp, num2, nomevar )
nomevar = 'field vector Eyp_dp'
call alloc_cdvet( Eyp_dp, num2, nomevar )
nomevar = 'field vector Ezp_dp'
call alloc_cdvet( Ezp_dp, num2, nomevar )


!======================================================================================================
!					 Primário do Receptor
!======================================================================================================


if( DPEH == 'ExDP' )then

	 do i=1, Nfreq

		w = (pi+pi)*vetf(i)
		zeta = (0.d0,1.d0)*w*mi0
	
		do j=1, Nobs

			Tx  = Mnos(no_obs(j),1)
			Ty  = Mnos(no_obs(j),2)
			Tz  = Mnos(no_obs(j),3)
	
			do l=nelebloc(1), nelebloc(Nbloc+1)-1
	
				inde = elebloc(l)
				nosG_e(:) = Melem(inde,:)
				areG_e(:) = Marestas(inde,:)
				x(:) = Mnos(nosG_e(:),1)
				y(:) = Mnos(nosG_e(:),2)
				z(:) = Mnos(nosG_e(:),3)
	
				CALL OMP_SET_NUM_THREADS( nprcss )
	
				!==================================================================================
				!==================================================================================
	
				!$OMP PARALLEL PRIVATE(p, num)
				!$OMP DO
				do p=1, 8
	
					num  = (l-1)*8 + p 				
					call dehx_xyz_loops( autr_fltJ0, numpJ0, autr_fltJ1, numpJ1, Tx, Ty, Iw, dsx, Tz, &
							      NC, hj, sigmas, neta0, zeta, x(p), y(p), z(p), Exp_dp(num), &
								                     Eyp_dp(num), Ezp_dp(num), Hx, Hy, Hz )
				end do
                                !$OMP END  DO
				!$OMP END PARALLEL
	
				!==================================================================================
				!==================================================================================
	
			end do

			write(unit=3010, IOSTAT=estts(1)) (dreal(Exp_dp(l)),l=1,num2)
			write(unit=3010, IOSTAT=estts(2)) (aimag(Exp_dp(l)),l=1,num2)	

			write(unit=3010, IOSTAT=estts(3)) (dreal(Eyp_dp(l)),l=1,num2)	
			write(unit=3010, IOSTAT=estts(4)) (aimag(Eyp_dp(l)),l=1,num2)	

			write(unit=3010, IOSTAT=estts(5)) (dreal(Ezp_dp(l)),l=1,num2)	
			write(unit=3010, IOSTAT=estts(6)) (aimag(Ezp_dp(l)),l=1,num2)	
	
			IF( estts(1) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 3010 drive file'
				PRINT*,'when writing the real part of the vector Exp_dp'
			END IF
			IF( estts(2) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 3010 drive file'
				PRINT*,'when writing the imaginary part of the vector Exp_dp' 	 
			END IF
			IF( estts(3) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 3010 drive file'
				PRINT*,'when writing the real part of the vector Eyp_dp'
			END IF
			IF( estts(4) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 3010 drive file'
				PRINT*,'when writing the imaginary part of the vector Eyp_dp' 	 
			END IF
			IF( estts(5) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 3010 drive file'
				PRINT*,'when writing the real part of the vector Ezp_dp'
			END IF
			IF( estts(6) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 3010 drive file'
				PRINT*,'when writing the imaginary part of the vector Ezp_dp' 	 
			END IF


			Exp_dp = (0.d0,0.d0)
			Eyp_dp = (0.d0,0.d0)
			Ezp_dp = (0.d0,0.d0)
	
		end do
	 end do

elseif( DPEH == 'EyDP' )then

	 do i=1, Nfreq
	
		w = (pi+pi)*vetf(i)
		zeta = (0.d0,1.d0)*w*mi0
	
		do j=1, Nobs
	
			Tx  = Mnos(no_obs(j),1)
			Ty  = Mnos(no_obs(j),2)
			Tz  = Mnos(no_obs(j),3)
	
			do l=nelebloc(1), nelebloc(Nbloc+1)-1
	
				inde = elebloc(l)
				nosG_e(:) = Melem(inde,:)
				areG_e(:) = Marestas(inde,:)
				x(:) = Mnos(nosG_e(:),1)
				y(:) = Mnos(nosG_e(:),2)
				z(:) = Mnos(nosG_e(:),3)
	
				CALL OMP_SET_NUM_THREADS( nprcss )
	
				!==================================================================================
				!==================================================================================
	
				!$OMP PARALLEL PRIVATE(p, num) 
				!$OMP DO
				do p=1, 8
	
					num  = (l-1)*8 + p 
					call dehy_xyz_loops(  autr_fltJ0, numpJ0, autr_fltJ1, numpJ1, Tx, Ty, Iw, dsy, Tz, &
							       NC, hj, sigmas, neta0, zeta, x(p), y(p), z(p), Exp_dp(num), &
								                      Eyp_dp(num), Ezp_dp(num), Hx, Hy, Hz )
				end do
                                !$OMP END  DO
				!$OMP END PARALLEL
	
				!==================================================================================
				!==================================================================================
	
			end do


			write(unit=3110, IOSTAT=estts(1)) (dreal(Exp_dp(l)),l=1,num2)
			write(unit=3110, IOSTAT=estts(2)) (aimag(Exp_dp(l)),l=1,num2)	

			write(unit=3110, IOSTAT=estts(3)) (dreal(Eyp_dp(l)),l=1,num2)	
			write(unit=3110, IOSTAT=estts(4)) (aimag(Eyp_dp(l)),l=1,num2)	

			write(unit=3110, IOSTAT=estts(5)) (dreal(Ezp_dp(l)),l=1,num2)	
			write(unit=3110, IOSTAT=estts(6)) (aimag(Ezp_dp(l)),l=1,num2)	

			IF( estts(1) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 3110 drive file'
				PRINT*,'when writing the real part of the vector Exp_dp'
			END IF
			IF( estts(2) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 3110 drive file'
				PRINT*,'when writing the imaginary part of the vector Exp_dp' 	 
			END IF
			IF( estts(3) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 3110 drive file'
				PRINT*,'when writing the real part of the vector Eyp_dp'
			END IF
			IF( estts(4) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 3110 drive file'
				PRINT*,'when writing the imaginary part of the vector Eyp_dp' 	 
			END IF
			IF( estts(5) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 3110 drive file'
				PRINT*,'when writing the real part of the vector Ezp_dp'
			END IF
			IF( estts(6) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 3110 drive file'
				PRINT*,'when writing the imaginary part of the vector Ezp_dp' 	 
			END IF
	
			Exp_dp = (0.d0,0.d0)
			Eyp_dp = (0.d0,0.d0)
			Ezp_dp = (0.d0,0.d0)
			
		end do
	 end do


elseif( DPEH == 'EzDP' )then

	 do i=1, Nfreq
	
		w = (pi+pi)*vetf(i)
		zeta = (0.d0,1.d0)*w*mi0
	
		do j=1, Nobs
	
			Tx  = Mnos(no_obs(j),1)
			Ty  = Mnos(no_obs(j),2)
			Tz  = Mnos(no_obs(j),3)
	
			do l=nelebloc(1), nelebloc(Nbloc+1)-1
	
				inde = elebloc(l)
				nosG_e(:) = Melem(inde,:)
				areG_e(:) = Marestas(inde,:)
				x(:) = Mnos(nosG_e(:),1)
				y(:) = Mnos(nosG_e(:),2)
				z(:) = Mnos(nosG_e(:),3)
	
				CALL OMP_SET_NUM_THREADS( nprcss )
	
				!==================================================================================
				!==================================================================================
	
				!$OMP PARALLEL PRIVATE(p, num) 
				!$OMP DO
				do p=1, 8
	
					num  = (l-1)*8 + p 
					call dev_xyz_loops(  autr_fltJ0, numpJ0, autr_fltJ1, numpJ1, Tx, Ty, Iw, dsz, Tz, &
							      NC, hj, sigmas, neta0, zeta, x(p), y(p), z(p), Exp_dp(num), &
								                     Eyp_dp(num), Ezp_dp(num), Hx, Hy, Hz )
				end do
                                !$OMP END  DO
				!$OMP END PARALLEL
	
				!==================================================================================
				!==================================================================================
	
			end do

			write(unit=3210, IOSTAT=estts(1)) (dreal(Exp_dp(l)),l=1,num2)
			write(unit=3210, IOSTAT=estts(2)) (aimag(Exp_dp(l)),l=1,num2)	

			write(unit=3210, IOSTAT=estts(3)) (dreal(Eyp_dp(l)),l=1,num2)	
			write(unit=3210, IOSTAT=estts(4)) (aimag(Eyp_dp(l)),l=1,num2)	

			write(unit=3210, IOSTAT=estts(5)) (dreal(Ezp_dp(l)),l=1,num2)	
			write(unit=3210, IOSTAT=estts(6)) (aimag(Ezp_dp(l)),l=1,num2)	
	
			IF( estts(1) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 3210 drive file'
				PRINT*,'when writing the real part of the vector Exp_dp'
			END IF
			IF( estts(2) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 3210 drive file'
				PRINT*,'when writing the imaginary part of the vector Exp_dp' 	 
			END IF
			IF( estts(3) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 3210 drive file'
				PRINT*,'when writing the real part of the vector Eyp_dp'
			END IF
			IF( estts(4) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 3210 drive file'
				PRINT*,'when writing the imaginary part of the vector Eyp_dp' 	 
			END IF
			IF( estts(5) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 3210 drive file'
				PRINT*,'when writing the real part of the vector Ezp_dp'
			END IF
			IF( estts(6) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 3210 drive file'
				PRINT*,'when writing the imaginary part of the vector Ezp_dp' 	 
			END IF

			Exp_dp = (0.d0,0.d0)
			Eyp_dp = (0.d0,0.d0)
			Ezp_dp = (0.d0,0.d0)
			
		end do
	 end do

elseif( DPEH == 'HxDP' )then

	 do i=1, Nfreq

		w = (pi+pi)*vetf(i)
		zeta = (0.d0,1.d0)*w*mi0
	
		do j=1, Nobs
	
			Tx  = Mnos(no_obs(j),1)
			Ty  = Mnos(no_obs(j),2)
			Tz  = Mnos(no_obs(j),3)
	
			do l=nelebloc(1), nelebloc(Nbloc+1)-1
	
				inde = elebloc(l)
				nosG_e(:) = Melem(inde,:)
				areG_e(:) = Marestas(inde,:)
				x(:) = Mnos(nosG_e(:),1)
				y(:) = Mnos(nosG_e(:),2)
				z(:) = Mnos(nosG_e(:),3)
	
				CALL OMP_SET_NUM_THREADS( nprcss )
	
				!==================================================================================
				!==================================================================================

				!$OMP PARALLEL PRIVATE(p, num) 
				!$OMP DO
				do p=1, 8
	
					num  = (l-1)*8 + p 
					call dmhx_xyz_loops( autr_fltJ0, numpJ0, autr_fltJ1, numpJ1, Tx, Ty, Iw, dsx, Tz, &
							      NC, hj, sigmas, neta0, zeta, x(p), y(p), z(p), Exp_dp(num), &
								                     Eyp_dp(num), Ezp_dp(num), Hx, Hy, Hz )
				end do
                                !$OMP END  DO
				!$OMP END PARALLEL
	
				!==================================================================================
				!==================================================================================
	
			end do

			write(unit=3310, IOSTAT=estts(1)) (dreal(Exp_dp(l)),l=1,num2)
			write(unit=3310, IOSTAT=estts(2)) (aimag(Exp_dp(l)),l=1,num2)	

			write(unit=3310, IOSTAT=estts(3)) (dreal(Eyp_dp(l)),l=1,num2)	
			write(unit=3310, IOSTAT=estts(4)) (aimag(Eyp_dp(l)),l=1,num2)	

			write(unit=3310, IOSTAT=estts(5)) (dreal(Ezp_dp(l)),l=1,num2)	
			write(unit=3310, IOSTAT=estts(6)) (aimag(Ezp_dp(l)),l=1,num2)	
			IF( estts(1) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 3310 drive file'
				PRINT*,'when writing the real part of the vector Exp_dp'
			END IF
			IF( estts(2) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 3310 drive file'
				PRINT*,'when writing the imaginary part of the vector Exp_dp' 	 
			END IF
			IF( estts(3) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 3310 drive file'
				PRINT*,'when writing the real part of the vector Eyp_dp'
			END IF
			IF( estts(4) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 3310 drive file'
				PRINT*,'when writing the imaginary part of the vector Eyp_dp' 	 
			END IF
			IF( estts(5) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 3310 drive file'
				PRINT*,'when writing the real part of the vector Ezp_dp'
			END IF
			IF( estts(6) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 3310 drive file'
				PRINT*,'when writing the imaginary part of the vector Ezp_dp' 	 
			END IF
	
			Exp_dp = (0.d0,0.d0)
			Eyp_dp = (0.d0,0.d0)
			Ezp_dp = (0.d0,0.d0)

		end do
	 end do	

elseif( DPEH == 'HyDP' )then
	
	 do i=1, Nfreq
	
		w = (pi+pi)*vetf(i)
		zeta = (0.d0,1.d0)*w*mi0
	
		do j=1, Nobs
	
			Tx  = Mnos(no_obs(j),1)
			Ty  = Mnos(no_obs(j),2)
			Tz  = Mnos(no_obs(j),3)
	
			do l=nelebloc(1), nelebloc(Nbloc+1)-1
	
				inde = elebloc(l)
				nosG_e(:) = Melem(inde,:)
				areG_e(:) = Marestas(inde,:)
				x(:) = Mnos(nosG_e(:),1)
				y(:) = Mnos(nosG_e(:),2)
				z(:) = Mnos(nosG_e(:),3)
	
				CALL OMP_SET_NUM_THREADS( nprcss )
	
				!==================================================================================
				!==================================================================================
	
				!$OMP PARALLEL PRIVATE(p, num) 
				!$OMP DO
				do p=1, 8
	
					num  = (l-1)*8 + p 
					call dmhx_xyz_loops( autr_fltJ0, numpJ0, autr_fltJ1, numpJ1, Ty, -Tx, Iw, dsx, Tz, &
							      NC, hj, sigmas, neta0, zeta, y(p), -x(p), z(p), Exp_dp(num), &
								                      Eyp_dp(num), Ezp_dp(num), Hx, Hy, Hz )
				end do
                                !$OMP END  DO
				!$OMP END PARALLEL
	
				!==================================================================================
				!==================================================================================
	
			end do

			write(unit=3410, IOSTAT=estts(3)) (dreal(-Eyp_dp(l)),l=1,num2)	
			write(unit=3410, IOSTAT=estts(4)) (aimag(-Eyp_dp(l)),l=1,num2)	

			write(unit=3410, IOSTAT=estts(1)) (dreal( Exp_dp(l)),l=1,num2)
			write(unit=3410, IOSTAT=estts(2)) (aimag( Exp_dp(l)),l=1,num2)	

			write(unit=3410, IOSTAT=estts(5)) (dreal( Ezp_dp(l)),l=1,num2)	
			write(unit=3410, IOSTAT=estts(6)) (aimag( Ezp_dp(l)),l=1,num2)
			IF( estts(1) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 3410 drive file'
				PRINT*,'when writing the real part of the vector Exp_dp'
			END IF
			IF( estts(2) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 3410 drive file'
				PRINT*,'when writing the imaginary part of the vector Exp_dp' 	 
			END IF
			IF( estts(3) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 3410 drive file'
				PRINT*,'when writing the real part of the vector Eyp_dp'
			END IF
			IF( estts(4) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 3410 drive file'
				PRINT*,'when writing the imaginary part of the vector Eyp_dp' 	 
			END IF
			IF( estts(5) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 3410 drive file'
				PRINT*,'when writing the real part of the vector Ezp_dp'
			END IF
			IF( estts(6) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 3410 drive file'
				PRINT*,'when writing the imaginary part of the vector Ezp_dp' 	 
			END IF
	
			Exp_dp = (0.d0,0.d0)
			Eyp_dp = (0.d0,0.d0)
			Ezp_dp = (0.d0,0.d0)
	
		end do
	 end do

elseif( DPEH == 'HzDP' )then

	 do i=1, Nfreq

		w = (pi+pi)*vetf(i)
		zeta = (0.d0,1.d0)*w*mi0
	
		do j=1, Nobs
	
			Tx  = Mnos(no_obs(j),1)
			Ty  = Mnos(no_obs(j),2)
			Tz  = Mnos(no_obs(j),3)
	
			do l=nelebloc(1), nelebloc(Nbloc+1)-1
	
				inde = elebloc(l)
				nosG_e(:) = Melem(inde,:)
				areG_e(:) = Marestas(inde,:)
				x(:) = Mnos(nosG_e(:),1)
				y(:) = Mnos(nosG_e(:),2)
				z(:) = Mnos(nosG_e(:),3)
	
				CALL OMP_SET_NUM_THREADS( nprcss )
	
				!==================================================================================
				!==================================================================================

				!$OMP PARALLEL PRIVATE(p, num) 
				!$OMP DO
				do p=1, 8
	
					num  = (l-1)*8 + p 
					call dmv_xyz_loops( autr_fltJ0, numpJ0, autr_fltJ1, numpJ1, Tx, Ty, dsz, Tz, &
							    NC, hj, sigmas, neta0, zeta, x(p), y(p), z(p), Exp_dp(num), &
								                  Eyp_dp(num), Ezp_dp(num), Hx, Hy, Hz )
				end do	
                                !$OMP END  DO	
				!$OMP END PARALLEL
	
				!==================================================================================
				!==================================================================================
	
			end do

			write(unit=3510, IOSTAT=estts(1)) (dreal(Exp_dp(l)),l=1,num2)
			write(unit=3510, IOSTAT=estts(2)) (aimag(Exp_dp(l)),l=1,num2)	
	
			write(unit=3510, IOSTAT=estts(3)) (dreal(Eyp_dp(l)),l=1,num2)	
			write(unit=3510, IOSTAT=estts(4)) (aimag(Eyp_dp(l)),l=1,num2)	
	
			write(unit=3510, IOSTAT=estts(5)) (dreal(Ezp_dp(l)),l=1,num2)	
			write(unit=3510, IOSTAT=estts(6)) (aimag(Ezp_dp(l)),l=1,num2)	
			IF( estts(1) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 3510 drive file'
				PRINT*,'when writing the real part of the vector Exp_dp'
			END IF
			IF( estts(2) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 3510 drive file'
				PRINT*,'when writing the imaginary part of the vector Exp_dp' 	 
			END IF
			IF( estts(3) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 3510 drive file'
				PRINT*,'when writing the real part of the vector Eyp_dp'
			END IF
			IF( estts(4) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 3510 drive file'
				PRINT*,'when writing the imaginary part of the vector Eyp_dp' 	 
			END IF
			IF( estts(5) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 3510 drive file'
				PRINT*,'when writing the real part of the vector Ezp_dp'
			END IF
			IF( estts(6) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 3510 drive file'
				PRINT*,'when writing the imaginary part of the vector Ezp_dp' 	 
			END IF

			Exp_dp = (0.d0,0.d0)
			Eyp_dp = (0.d0,0.d0)
			Ezp_dp = (0.d0,0.d0)
	
		end do
	 end do	

end if

 deallocate( Exp_dp, Eyp_dp, Ezp_dp )
 call close_arq_DPEM1D_adj( DPEH )

END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE DPEH1D( DPEH, nprcss )
IMPLICIT NONE	
 INTEGER, INTENT(in)			:: nprcss	
 CHARACTER(LEN=20), INTENT(in)		:: DPEH
 REAL(dp)   , ALLOCATABLE, DIMENSION(:) :: sigmas
 COMPLEX(dp), ALLOCATABLE, DIMENSION(:) :: Exp_dp, Eyp_dp, Ezp_dp

 COMPLEX(dp)	:: Hx, Hy, Hz
 COMPLEX(dp)	:: Ex_aux, Ey_aux, Ez_aux, z_aux, neta0, zeta
 REAL(dp)	:: w, Iw, dsx, dsy, dsz, Tx, Ty, Tz, t1, t2
 REAL(dp)	:: x(8), y(8), z(8)
 INTEGER	:: i, j, k, p, l, inde, nosG_e(8), areG_e(12)
 INTEGER	:: num, num2, id, OMP_GET_THREAD_NUM,  IERR, estts(6) 

 CHARACTER(len=50) :: nomevar

 call openRW_arq_DPEM1D( DPEH )

allocate( sigmas(NC))

sigmas = 1.d0/roj

Iw  = ISOURCE
dsx = SDIMENS
dsy = SDIMENS
dsz = SDIMENS

neta0 = 1.d-12  !(0.d0,1.d0)*omega*eps    !

num = (nelebloc(Nbloc+1)-1)*8
num2 = num
num  = 0


if( DPEH /= 'HCL')then

	nomevar = 'field vector Exp_dp'
	call alloc_cdvet( Exp_dp, num2, nomevar )
	nomevar = 'field vector Eyp_dp'
	call alloc_cdvet( Eyp_dp, num2, nomevar )
	nomevar = 'field vector Ezp_dp'
	call alloc_cdvet( Ezp_dp, num2, nomevar )
else
	nomevar = 'field vector Exp_dp'
	call alloc_cdvet( Exp_dp, num2, nomevar )
	nomevar = 'field vector Eyp_dp'
	call alloc_cdvet( Eyp_dp, num2, nomevar )
endif

!======================================================================================================
!					 Primário do Transmissor
!======================================================================================================


if( DPEH == 'ExDP' )then

	 do i=1, Nfreq
	
		w = (pi+pi)*vetf(i)
		zeta = (0.d0,1.d0)*w*mi0
	
		do j=1, Ntransm

			Tx  = cxT(j)
			Ty  = cyT(j)
			Tz  = czT(j)
	
			do l=nelebloc(1), nelebloc(Nbloc+1)-1; 
		
				inde = elebloc(l)
				nosG_e(:) = Melem(inde,:)
				areG_e(:) = Marestas(inde,:)
				x(:) = Mnos(nosG_e(:),1)
				y(:) = Mnos(nosG_e(:),2)
				z(:) = Mnos(nosG_e(:),3)
	
				CALL OMP_SET_NUM_THREADS( nprcss )
	
				!==================================================================================
				!==================================================================================
				!$OMP PARALLEL PRIVATE(p, num) 
				!$OMP DO
				do p=1, 8
	
					num  = (l-1)*8 + p 				
					call dehx_xyz_loops( autr_fltJ0, numpJ0, autr_fltJ1, numpJ1, Tx, Ty, Iw, dsx, Tz, &
							      NC, hj, sigmas, neta0, zeta, x(p), y(p), z(p), Exp_dp(num), &
								                     Eyp_dp(num), Ezp_dp(num), Hx, Hy, Hz )
				end do
                                !$OMP END  DO
				!$OMP END PARALLEL
	
				!==================================================================================
				!==================================================================================
	
			end do

			write(unit=4000, IOSTAT=estts(1)) (dreal(Exp_dp(l)),l=1,num2)
			write(unit=4000, IOSTAT=estts(2)) (aimag(Exp_dp(l)),l=1,num2)	

			write(unit=4000, IOSTAT=estts(3)) (dreal(Eyp_dp(l)),l=1,num2)	
			write(unit=4000, IOSTAT=estts(4)) (aimag(Eyp_dp(l)),l=1,num2)	

			write(unit=4000, IOSTAT=estts(5)) (dreal(Ezp_dp(l)),l=1,num2)	
			write(unit=4000, IOSTAT=estts(6)) (aimag(Ezp_dp(l)),l=1,num2)	

			IF( estts(1) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 4000 drive file'
				PRINT*,'when writing the real part of the vector Exp_dp'
			END IF
			IF( estts(2) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 4000 drive file'
				PRINT*,'when writing the imaginary part of the vector Exp_dp' 	 
			END IF
			IF( estts(3) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 4000 drive file'
				PRINT*,'when writing the real part of the vector Eyp_dp'
			END IF
			IF( estts(4) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 4000 drive file'
				PRINT*,'when writing the imaginary part of the vector Eyp_dp' 	 
			END IF
			IF( estts(5) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 4000 drive file'
				PRINT*,'when writing the real part of the vector Ezp_dp'
			END IF
			IF( estts(6) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 4000 drive file'
				PRINT*,'when writing the imaginary part of the vector Ezp_dp' 	 
			END IF

			Exp_dp = (0.d0,0.d0)
			Eyp_dp = (0.d0,0.d0)
			Ezp_dp = (0.d0,0.d0)
	
		end do
	 end do

elseif( DPEH == 'EyDP' )then

	 do i=1, Nfreq

		w = (pi+pi)*vetf(i)
		zeta = (0.d0,1.d0)*w*mi0
	
		do j=1, Ntransm

			Tx  = cxT(j)
			Ty  = cyT(j)
			Tz  = czT(j)
	
			do l=nelebloc(1), nelebloc(Nbloc+1)-1;
	
				inde = elebloc(l)
				nosG_e(:) = Melem(inde,:)
				areG_e(:) = Marestas(inde,:)
				x(:) = Mnos(nosG_e(:),1)
				y(:) = Mnos(nosG_e(:),2)
				z(:) = Mnos(nosG_e(:),3)
	
				CALL OMP_SET_NUM_THREADS( nprcss )
	
				!==================================================================================
				!==================================================================================
	
				!$OMP PARALLEL PRIVATE(p, num) 
				!$OMP DO
				do p=1, 8
	
					num  = (l-1)*8 + p 
					call dehy_xyz_loops(  autr_fltJ0, numpJ0, autr_fltJ1, numpJ1, Tx, Ty, Iw, dsy, Tz, &
							       NC, hj, sigmas, neta0, zeta, x(p), y(p), z(p), Exp_dp(num), &
								                      Eyp_dp(num), Ezp_dp(num), Hx, Hy, Hz )
				end do
                                !$OMP END  DO
				!$OMP END PARALLEL
	
				!==================================================================================
				!==================================================================================
	
			end do

			write(unit=4100, IOSTAT=estts(1)) (dreal(Exp_dp(l)),l=1,num2)
			write(unit=4100, IOSTAT=estts(2)) (aimag(Exp_dp(l)),l=1,num2)	

			write(unit=4100, IOSTAT=estts(3)) (dreal(Eyp_dp(l)),l=1,num2)	
			write(unit=4100, IOSTAT=estts(4)) (aimag(Eyp_dp(l)),l=1,num2)	

			write(unit=4100, IOSTAT=estts(5)) (dreal(Ezp_dp(l)),l=1,num2)	
			write(unit=4100, IOSTAT=estts(6)) (aimag(Ezp_dp(l)),l=1,num2)	

			IF( estts(1) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 4100 drive file'
				PRINT*,'when writing the real part of the vector Exp_dp'
			END IF
			IF( estts(2) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 4100 drive file'
				PRINT*,'when writing the imaginary part of the vector Exp_dp' 	 
			END IF
			IF( estts(3) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 4100 drive file'
				PRINT*,'when writing the real part of the vector Eyp_dp'
			END IF
			IF( estts(4) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 4100 drive file'
				PRINT*,'when writing the imaginary part of the vector Eyp_dp' 	 
			END IF
			IF( estts(5) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 4100 drive file'
				PRINT*,'when writing the real part of the vector Ezp_dp'
			END IF
			IF( estts(6) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 4100 drive file'
				PRINT*,'when writing the imaginary part of the vector Ezp_dp' 	 
			END IF

			Exp_dp = (0.d0,0.d0)
			Eyp_dp = (0.d0,0.d0)
			Ezp_dp = (0.d0,0.d0)
			
		end do
	 end do


elseif( DPEH == 'EzDP' )then

	 do i=1, Nfreq

		w = (pi+pi)*vetf(i)
		zeta = (0.d0,1.d0)*w*mi0
	
		do j=1, Ntransm
	
			Tx  = cxT(j)
			Ty  = cyT(j)
			Tz  = czT(j)
	
			do l=nelebloc(1), nelebloc(Nbloc+1)-1; 
	
				inde = elebloc(l)
				nosG_e(:) = Melem(inde,:)
				areG_e(:) = Marestas(inde,:)
				x(:) = Mnos(nosG_e(:),1)
				y(:) = Mnos(nosG_e(:),2)
				z(:) = Mnos(nosG_e(:),3)
	
				CALL OMP_SET_NUM_THREADS( nprcss )
	
				!==================================================================================
				!==================================================================================
	
				!$OMP PARALLEL PRIVATE(p, num) 
				!$OMP DO
				do p=1, 8
	
					num  = (l-1)*8 + p 
					call dev_xyz_loops( autr_fltJ0, numpJ0, autr_fltJ1, numpJ1, Tx, Ty, Iw, dsy, Tz, &
							     NC, hj, sigmas, neta0, zeta, x(p), y(p), z(p), Exp_dp(num), &
										  Eyp_dp(num), Ezp_dp(num), Hx, Hy, Hz )
				end do
                                !$OMP END  DO
				!$OMP END PARALLEL
	
				!==================================================================================
				!==================================================================================
	
			end do
	
			write(unit=4200, IOSTAT=estts(1)) (dreal(Exp_dp(l)),l=1,num2)
			write(unit=4200, IOSTAT=estts(2)) (aimag(Exp_dp(l)),l=1,num2)	

			write(unit=4200, IOSTAT=estts(3)) (dreal(Eyp_dp(l)),l=1,num2)	
			write(unit=4200, IOSTAT=estts(4)) (aimag(Eyp_dp(l)),l=1,num2)	

			write(unit=4200, IOSTAT=estts(5)) (dreal(Ezp_dp(l)),l=1,num2)	
			write(unit=4200, IOSTAT=estts(6)) (aimag(Ezp_dp(l)),l=1,num2)	

			IF( estts(1) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 4200 drive file'
				PRINT*,'when writing the real part of the vector Exp_dp'
			END IF
			IF( estts(2) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 4200 drive file'
				PRINT*,'when writing the imaginary part of the vector Exp_dp' 	 
			END IF
			IF( estts(3) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 4200 drive file'
				PRINT*,'when writing the real part of the vector Eyp_dp'
			END IF
			IF( estts(4) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 4200 drive file'
				PRINT*,'when writing the imaginary part of the vector Eyp_dp' 	 
			END IF
			IF( estts(5) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 4200 drive file'
				PRINT*,'when writing the real part of the vector Ezp_dp'
			END IF
			IF( estts(6) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 4200 drive file'
				PRINT*,'when writing the imaginary part of the vector Ezp_dp' 	 
			END IF

			Exp_dp = (0.d0,0.d0)
			Eyp_dp = (0.d0,0.d0)
			Ezp_dp = (0.d0,0.d0)
			
		end do
	 end do


elseif( DPEH == 'HxDP' )then

	 do i=1, Nfreq

		w = (pi+pi)*vetf(i)
		zeta = (0.d0,1.d0)*w*mi0
	
		do j=1, Ntransm

			Tx  = cxT(j)
			Ty  = cyT(j)
			Tz  = czT(j)
	
			do l=nelebloc(1), nelebloc(Nbloc+1)-1;
	
				inde = elebloc(l)
				nosG_e(:) = Melem(inde,:)
				areG_e(:) = Marestas(inde,:)
				x(:) = Mnos(nosG_e(:),1)
				y(:) = Mnos(nosG_e(:),2)
				z(:) = Mnos(nosG_e(:),3)
	
				CALL OMP_SET_NUM_THREADS( nprcss )
	
				!==================================================================================
				!==================================================================================

				!$OMP PARALLEL PRIVATE(p, num) 
				!$OMP DO
				do p=1, 8
	
					num  = (l-1)*8 + p 
					call dmhx_xyz_loops( autr_fltJ0, numpJ0, autr_fltJ1, numpJ1, Tx, Ty, Iw, dsx, Tz, &
							      NC, hj, sigmas, neta0, zeta, x(p), y(p), z(p), Exp_dp(num), &
								                     Eyp_dp(num), Ezp_dp(num), Hx, Hy, Hz )
				end do
                                !$OMP END  DO
				!$OMP END PARALLEL
	
				!==================================================================================
				!==================================================================================
	
			end do

			write(unit=4300, IOSTAT=estts(1)) (dreal(Exp_dp(l)),l=1,num2)
			write(unit=4300, IOSTAT=estts(2)) (aimag(Exp_dp(l)),l=1,num2)	

			write(unit=4300, IOSTAT=estts(3)) (dreal(Eyp_dp(l)),l=1,num2)	
			write(unit=4300, IOSTAT=estts(4)) (aimag(Eyp_dp(l)),l=1,num2)	

			write(unit=4300, IOSTAT=estts(5)) (dreal(Ezp_dp(l)),l=1,num2)	
			write(unit=4300, IOSTAT=estts(6)) (aimag(Ezp_dp(l)),l=1,num2)	

			IF( estts(1) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 4300 drive file'
				PRINT*,'when writing the real part of the vector Exp_dp'
			END IF
			IF( estts(2) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 4300 drive file'
				PRINT*,'when writing the imaginary part of the vector Exp_dp' 	 
			END IF
			IF( estts(3) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 4300 drive file'
				PRINT*,'when writing the real part of the vector Eyp_dp'
			END IF
			IF( estts(4) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 4300 drive file'
				PRINT*,'when writing the imaginary part of the vector Eyp_dp' 	 
			END IF
			IF( estts(5) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 4300 drive file'
				PRINT*,'when writing the real part of the vector Ezp_dp'
			END IF
			IF( estts(6) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 4300 drive file'
				PRINT*,'when writing the imaginary part of the vector Ezp_dp' 	 
			END IF

			Exp_dp = (0.d0,0.d0)
			Eyp_dp = (0.d0,0.d0)
			Ezp_dp = (0.d0,0.d0)

		end do
	 end do	

elseif( DPEH == 'HyDP' )then
	
	 do i=1, Nfreq
	
		w = (pi+pi)*vetf(i)
		zeta = (0.d0,1.d0)*w*mi0
	
		do j=1, Ntransm
	
			Tx  = cxT(j)
			Ty  = cyT(j)
			Tz  = czT(j)
	
			do l=nelebloc(1), nelebloc(Nbloc+1)-1; 
	
				inde = elebloc(l)
				nosG_e(:) = Melem(inde,:)
				areG_e(:) = Marestas(inde,:)
				x(:) = Mnos(nosG_e(:),1)
				y(:) = Mnos(nosG_e(:),2)
				z(:) = Mnos(nosG_e(:),3)
	
				CALL OMP_SET_NUM_THREADS( nprcss )
	
				!==================================================================================
				!==================================================================================
	
				!$OMP PARALLEL PRIVATE(p, num) 
				!$OMP DO
				do p=1, 8
	
					num  = (l-1)*8 + p 
					call dmhx_xyz_loops( autr_fltJ0, numpJ0, autr_fltJ1, numpJ1, Ty, -Tx, Iw, dsx, Tz, &
							      NC, hj, sigmas, neta0, zeta, y(p), -x(p), z(p), Exp_dp(num), &
								                      Eyp_dp(num), Ezp_dp(num), Hx, Hy, Hz )
				end do
                                !$OMP END  DO
				!$OMP END PARALLEL
	
				!==================================================================================
				!==================================================================================
	
			end do
	
			write(unit=4400, IOSTAT=estts(3)) (dreal(-Eyp_dp(l)),l=1,num2)	
			write(unit=4400, IOSTAT=estts(4)) (aimag(-Eyp_dp(l)),l=1,num2)	

			write(unit=4400, IOSTAT=estts(1)) (dreal( Exp_dp(l)),l=1,num2)
			write(unit=4400, IOSTAT=estts(2)) (aimag( Exp_dp(l)),l=1,num2)	

			write(unit=4400, IOSTAT=estts(5)) (dreal( Ezp_dp(l)),l=1,num2)	
			write(unit=4400, IOSTAT=estts(6)) (aimag( Ezp_dp(l)),l=1,num2)	

			IF( estts(1) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 4400 drive file'
				PRINT*,'when writing the real part of the vector Exp_dp'
			END IF
			IF( estts(2) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 4400 drive file'
				PRINT*,'when writing the imaginary part of the vector Exp_dp' 	 
			END IF
			IF( estts(3) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 4400 drive file'
				PRINT*,'when writing the real part of the vector Eyp_dp'
			END IF
			IF( estts(4) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 4400 drive file'
				PRINT*,'when writing the imaginary part of the vector Eyp_dp' 	 
			END IF
			IF( estts(5) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 4400 drive file'
				PRINT*,'when writing the real part of the vector Ezp_dp'
			END IF
			IF( estts(6) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 4400 drive file'
				PRINT*,'when writing the imaginary part of the vector Ezp_dp' 	 
			END IF
	
			Exp_dp = (0.d0,0.d0)
			Eyp_dp = (0.d0,0.d0)
			Ezp_dp = (0.d0,0.d0)
	
		end do
	 end do

elseif( DPEH == 'HzDP' )then

	 do i=1, Nfreq

		w = (pi+pi)*vetf(i)
		zeta = (0.d0,1.d0)*w*mi0
	
		do j=1, Ntransm
	
			Tx  = cxT(j)
			Ty  = cyT(j)
			Tz  = czT(j)
	
			do l=nelebloc(1), nelebloc(Nbloc+1)-1; 
	
				inde = elebloc(l)
				nosG_e(:) = Melem(inde,:)
				areG_e(:) = Marestas(inde,:)
				x(:) = Mnos(nosG_e(:),1)
				y(:) = Mnos(nosG_e(:),2)
				z(:) = Mnos(nosG_e(:),3)
	
				CALL OMP_SET_NUM_THREADS( nprcss )
	
				!==================================================================================
				!==================================================================================

				!$OMP PARALLEL PRIVATE(p, num)
				!$OMP DO
				do p=1, 8
	
					num  = (l-1)*8 + p 
					call dmv_xyz_loops( autr_fltJ0, numpJ0, autr_fltJ1, numpJ1, Tx, Ty, dsz, Tz, &
							    NC, hj, sigmas, neta0, zeta, x(p), y(p), z(p), Exp_dp(num), &
										Eyp_dp(num), Ezp_dp(num), Hx, Hy, Hz )
				end do
                                !$OMP END  DO
				!$OMP END PARALLEL
	
				!==================================================================================
				!==================================================================================
	
			end do

			write(unit=4500, IOSTAT=estts(1)) (dreal(Exp_dp(l)),l=1,num2)
			write(unit=4500, IOSTAT=estts(2)) (aimag(Exp_dp(l)),l=1,num2)	

			write(unit=4500, IOSTAT=estts(3)) (dreal(Eyp_dp(l)),l=1,num2)	
			write(unit=4500, IOSTAT=estts(4)) (aimag(Eyp_dp(l)),l=1,num2)	

			write(unit=4500, IOSTAT=estts(5)) (dreal(Ezp_dp(l)),l=1,num2)	
			write(unit=4500, IOSTAT=estts(6)) (aimag(Ezp_dp(l)),l=1,num2)	

			IF( estts(1) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 4500 drive file'
				PRINT*,'when writing the real part of the vector Exp_dp'
			END IF
			IF( estts(2) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 4500 drive file'
				PRINT*,'when writing the imaginary part of the vector Exp_dp' 	 
			END IF
			IF( estts(3) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 4500 drive file'
				PRINT*,'when writing the real part of the vector Eyp_dp'
			END IF
			IF( estts(4) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 4500 drive file'
				PRINT*,'when writing the imaginary part of the vector Eyp_dp' 	 
			END IF
			IF( estts(5) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 4500 drive file'
				PRINT*,'when writing the real part of the vector Ezp_dp'
			END IF
			IF( estts(6) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 4500 drive file'
				PRINT*,'when writing the imaginary part of the vector Ezp_dp' 	 
			END IF

			Exp_dp = (0.d0,0.d0)
			Eyp_dp = (0.d0,0.d0)
			Ezp_dp = (0.d0,0.d0)

		end do
	 end do	

elseif( DPEH == 'HCL' )then


         Ex_aux  = (0.d0,0.d0)
	 do i=1, Nfreq

		w = vetf(i)
	
		do j=1, Ntransm
	
			Tx  = cxT(j)
			Ty  = cyT(j)
			Tz  = czT(j)
	
			do l=nelebloc(1), nelebloc(Nbloc+1)-1
	
				inde = elebloc(l)
				nosG_e(:) = Melem(inde,:)
				areG_e(:) = Marestas(inde,:)
				x(:) = Mnos(nosG_e(:),1)
				y(:) = Mnos(nosG_e(:),2)
				z(:) = Mnos(nosG_e(:),3)

				!==================================================================================
				!==================================================================================

				do p=1, 8
					num  = (l-1)*8 + p 
					call BCH1DQWE_pll_campo_E( nprcss, w, Tx, Ty, Tz, x(p), y(p), z(p), & 
								            Ex_aux, Exp_dp(num), Eyp_dp(num) )
				end do

				!==================================================================================
				!==================================================================================

			end do

			write(unit=4600, IOSTAT=estts(1)) (dreal(Exp_dp(l)),l=1,num2)
			write(unit=4600, IOSTAT=estts(2)) (aimag(Exp_dp(l)),l=1,num2)	

			write(unit=4600, IOSTAT=estts(3)) (dreal(Eyp_dp(l)),l=1,num2)	
			write(unit=4600, IOSTAT=estts(4)) (aimag(Eyp_dp(l)),l=1,num2)	

			IF( estts(1) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 4600 drive file'
				PRINT*,'when writing the real part of the vector Exp_dp'
			END IF
			IF( estts(2) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 4600 drive file'
				PRINT*,'when writing the imaginary part of the vector Exp_dp' 	 
			END IF
			IF( estts(3) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 4600 drive file'
				PRINT*,'when writing the real part of the vector Eyp_dp'
			END IF
			IF( estts(4) /= 0 )THEN
				PRINT*,'An error occurred while writing to the 4600 drive file'
				PRINT*,'when writing the imaginary part of the vector Eyp_dp' 	 
			END IF

			Exp_dp = (0.d0,0.d0)
			Eyp_dp = (0.d0,0.d0)


		end do
	 end do	
end if

 deallocate( Exp_dp, Eyp_dp, Ezp_dp )

 call close_arq_DPEM1D( DPEH )
 call openRW_arq_DPEM1D( DPEH )

END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE BCH1D( nprcss )
IMPLICIT NONE	
 INTEGER, INTENT(in)			:: nprcss	
 REAL(dp)   , ALLOCATABLE, DIMENSION(:) :: sigmas
 COMPLEX(dp), ALLOCATABLE, DIMENSION(:) :: Exp_dp, Eyp_dp

 COMPLEX(dp)	:: Hx, Hy, Hz
 COMPLEX(dp)	:: Ex_aux, neta0, zeta
 REAL(dp)	:: w, Tx, Ty, Tz
 REAL(dp)	:: x(8), y(8), z(8)
 INTEGER	:: i, j, p, l, inde, nosG_e(8), areG_e(12)
 INTEGER	:: IERR , num, num2, id, proc_num, thread_num
 INTEGER	:: OMP_GET_NUM_PROCS, notnun(1), Nprcss0, estts(4)
 LOGICAL        :: arqexist, arqopened
 character(len=20) :: DPEH, arqpad, arqwrite, arqacess

 CHARACTER(len=50) :: nomevar

 DPEH = 'HCL'
 call openRW_arq_DPEM1D( DPEH )
 
 allocate( sigmas(NC) )

 sigmas = 1.d0/roj
 neta0 = 1.d-12
 Ex_aux  = (0.d0,0.d0)

 num  = (nelebloc(Nbloc+1)-1)*8
 num2 = num
 num  = 0
 
 nomevar = 'field vector Exp_dp'
 call alloc_cdvet( Exp_dp, num2, nomevar )
 nomevar = 'field vector Eyp_dp'
 call alloc_cdvet( Eyp_dp, num2, nomevar )

!======================================================================================================
!					 Primário do Transmissor
!======================================================================================================

 do i=1, Nfreq
	w = vetf(i)

	do j=1, Ntransm

		Tx  = cxT(j)
		Ty  = cyT(j)
		Tz  = czT(j)

		do l=nelebloc(1), nelebloc(Nbloc+1)-1!

			inde = elebloc(l)
			nosG_e(:) = Melem(inde,:)
			areG_e(:) = Marestas(inde,:)
			x(:) = Mnos(nosG_e(:),1)
			y(:) = Mnos(nosG_e(:),2)
			z(:) = Mnos(nosG_e(:),3)

			!==================================================================================
			!==================================================================================

			do p=1, 8
				num  = (l-1)*8 + p 
				call BCH1DQWE_pll_campo_E( nprcss, w, Tx, Ty, Tz, x(p), y(p), z(p), Ex_aux, Exp_dp(num), Eyp_dp(num) )

			end do

			!==================================================================================
			!==================================================================================

		end do

		write(unit=4600, IOSTAT=estts(1)) (dreal(Exp_dp(l)),l=1,num2)
		write(unit=4600, IOSTAT=estts(2)) (aimag(Exp_dp(l)),l=1,num2)
		write(unit=4600, IOSTAT=estts(3)) (dreal(Eyp_dp(l)),l=1,num2)
		write(unit=4600, IOSTAT=estts(4)) (aimag(Eyp_dp(l)),l=1,num2)

		IF( estts(1) /= 0 )THEN
			PRINT*,'An error occurred while writing to the 4600 drive file'
			PRINT*,'when writing the real part of the vector Exp_dp'
		END IF
		IF( estts(2) /= 0 )THEN
			PRINT*,'An error occurred while writing to the 4600 drive file'
			PRINT*,'when writing the imaginary part of the vector Exp_dp' 	 
		END IF
		IF( estts(3) /= 0 )THEN
			PRINT*,'An error occurred while writing to the 4600 drive file'
			PRINT*,'when writing the real part of the vector Eyp_dp'
		END IF
		IF( estts(4) /= 0 )THEN
			PRINT*,'An error occurred while writing to the 4600 drive file'
			PRINT*,'when writing the imaginary part of the vector Eyp_dp' 	 
		END IF

		Exp_dp = (0.d0,0.d0)
		Eyp_dp = (0.d0,0.d0)

	end do
 end do	


 deallocate( Exp_dp, Eyp_dp)


 call close_arq_DPEM1D( DPEH )
 call openRW_arq_DPEM1D( DPEH )

END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE campo_H_E_obs_MT( valf, md, vsol, Es_obs, Hs_obs, Et_obs, Ht_obs )
IMPLICIT NONE
 COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:), INTENT(out) :: Es_obs, Hs_obs, Et_obs, Ht_obs
 COMPLEX(dp), 		   DIMENSION(:)  , INTENT(in)  :: vsol
 REAL(dp)  			    	 , INTENT(in)  :: valf
 INTEGER                                 , INTENT(in)  :: md

 COMPLEX(dp), DIMENSION(Nobs,3)	:: Ep_obs, Hp_obs
 COMPLEX(dp)			:: EyHy, HxEX, Ephi, Ex, Ey, Hr, Hx, Hy, Hz
 REAL(dp)			:: x, y, z, w
 INTEGER			:: i, j, N

N = Nobs

allocate( Es_obs(N,3), Hs_obs(N,3), Et_obs(N,3), Ht_obs(N,3)  )

Ep_obs = (0.d0, 0.d0)
Hp_obs = (0.d0, 0.d0)

Es_obs = (0.d0, 0.d0)
Hs_obs = (0.d0, 0.d0)

Et_obs = (0.d0, 0.d0)
Ht_obs = (0.d0, 0.d0)
w = (pi+pi)*valf

call campo_Hs_obs( vsol, Hs_obs, w )
				

if( md == 1 )then

	do i=1, N

		z = Mnos(no_obs(i),3)
		call EyHy_TETM( z, valf, 1, EyHy, HxEx ) ! parâmetro interiro indica o modo TE-(Ey, Hx) Usando TE+TM 
		Ep_obs(i,2) = EyHy
		Hp_obs(i,1) = HxEX

		call EyHy_TETM( z, valf, 2, EyHy, HxEx ) ! parâmetro interiro indica o modo TM-(Ey, Hx) Usando TE+TM 
		Ep_obs(i,1) = HxEx
		Hp_obs(i,2) = EyHy

	end do

else if( md == 2 )then

	do i=1, N

		z = Mnos(no_obs(i),3)
		call EyHy_TETM( z, valf, 1, EyHy, HxEx ) ! parâmetro inteiro que indica o modo TE-(Ey, Hx) Usando TE+TM 90º Ant-Horário
		Ep_obs(i,1) = -1.d0*EyHy
		Hp_obs(i,2) = HxEX
	
		call EyHy_TETM( z, valf, 2, EyHy, HxEx ) ! parâmetro inteiro que indica o modo TM-(Ey, Hx) Usando TE+TM 90º Ant-Horário
		Ep_obs(i,2) = HxEx
		Hp_obs(i,1) = -1.d0*EyHy
	end do

end if


do i=1, N


	Es_obs(i,1) = ( vsol(Marestas(eleobs(indobs(i)),4)) + &
			vsol(Marestas(eleobs(indobs(i)+1),4)) )/2.d0

	Es_obs(i,2) = ( vsol(Marestas(eleobs(indobs(i)),8)) + &
			vsol(Marestas(eleobs(indobs(i)+2),8)) )/2.d0

	Es_obs(i,3) = ( vsol(Marestas(eleobs(indobs(i)),12)) + &
			vsol(Marestas(eleobs(indobs(i)+4),12)) )/2.d0

	Et_obs(i,:) = Ep_obs(i,:) + Es_obs(i,:) 
	Ht_obs(i,:) = Hp_obs(i,:) + Hs_obs(i,:) 

end do


END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE campo_H_E_obs_DPEH( valf, Tx, Ty, Tz, FDP, vsol, Es_obs, Hs_obs, Et_obs, Ht_obs )
IMPLICIT NONE
 COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:), INTENT(out)	:: Es_obs, Hs_obs, Et_obs, Ht_obs
 COMPLEX(dp), 		   DIMENSION(:)  , INTENT(in)	:: vsol
 REAL(dp), 		 		   INTENT(in)	:: valf, Tx, Ty, Tz
 CHARACTER(LEN=20), 			   INTENT(in)	:: FDP		

 COMPLEX(dp), DIMENSION(Nobs,3)	:: Ep_obs, Hp_obs
 COMPLEX(dp)					:: Ephi, Ex, Ey, Ez, Hr, Hx, Hy, Hz, neta0, zeta
 REAL(dp)						:: x, y, z, dsx, dsy, dsz, Iw, w, sigmas(NC)
 INTEGER						:: i, j, N

N = Nobs
Iw  = ISOURCE
dsx = SDIMENS
dsy = SDIMENS
dsz = SDIMENS

w = (pi+pi)*valf
sigmas = 1.d0/roj
neta0  = 1.d-12  !(0.d0,1.d0)*omega*eps    !
zeta   = (0.d0,1.d0)*w*mi0

allocate( Es_obs(N,3), Hs_obs(N,3), Et_obs(N,3), Ht_obs(N,3)  )

Ep_obs = (0.d0, 0.d0)
Hp_obs = (0.d0, 0.d0)

Es_obs = (0.d0, 0.d0)
Hs_obs = (0.d0, 0.d0)

Et_obs = (0.d0, 0.d0)
Ht_obs = (0.d0, 0.d0)

call campo_Hs_obs( vsol, Hs_obs, w )
				

if( FDP == 'ExDP' )then

	do i=1, N

		z = Mnos(no_obs(i),3)
		y = Mnos(no_obs(i),2)
		x = Mnos(no_obs(i),1)

		call dehx_xyz_loops( autr_fltJ0, numpJ0, autr_fltJ1, numpJ1, Tx, Ty, Iw, dsx, Tz, NC, &
					    hj, sigmas, neta0, zeta, x, y, z, Ex, Ey, Ez , Hx, Hy, Hz )

		Ep_obs(i,1) = Ex
		Ep_obs(i,2) = Ey
		Ep_obs(i,3) = Ez
		Hp_obs(i,1) = Hx
		Hp_obs(i,2) = Hy
		Hp_obs(i,3) = Hz

	end do

elseif( FDP == 'EyDP' )then

	do i=1, N

		z = Mnos(no_obs(i),3)
		y = Mnos(no_obs(i),2)
		x = Mnos(no_obs(i),1)

		call dehy_xyz_loops( autr_fltJ0, numpJ0, autr_fltJ1, numpJ1, Tx, Ty, Iw, dsy, Tz, NC, &
					   hj, sigmas, neta0, zeta, x, y, z, Ex, Ey, Ez , Hx, Hy, Hz )

		Ep_obs(i,1) = Ex
		Ep_obs(i,2) = Ey
		Ep_obs(i,3) = Ez
		Hp_obs(i,1) = Hx
		Hp_obs(i,2) = Hy
		Hp_obs(i,3) = Hz

	end do

elseif( FDP == 'EzDP' )then

	do i=1, N

		z = Mnos(no_obs(i),3)
		y = Mnos(no_obs(i),2)
		x = Mnos(no_obs(i),1)

		call dev_xyz_loops( autr_fltJ0, numpJ0, autr_fltJ1, numpJ1, Tx, Ty, Iw, dsz, Tz, NC, &
							     hj, sigmas, neta0, zeta, x, y, z, Ex, Ey, Ez , Hx, Hy, Hz )

		Ep_obs(i,1) = Ex
		Ep_obs(i,2) = Ey
		Ep_obs(i,3) = Ez
		Hp_obs(i,1) = Hx
		Hp_obs(i,2) = Hy
		Hp_obs(i,3) = Hz

	end do

elseif( FDP == 'HxDP' )then

	do i=1, N

		z = Mnos(no_obs(i),3)
		y = Mnos(no_obs(i),2)
		x = Mnos(no_obs(i),1)

		call dmhx_xyz_loops( autr_fltJ0, numpJ0, autr_fltJ1, numpJ1, Tx, Ty, Iw, dsx, Tz, NC, &
							     hj, sigmas, neta0, zeta, x, y, z, Ex, Ey, Ez , Hx, Hy, Hz )

		Ep_obs(i,1) = Ex
		Ep_obs(i,2) = Ey
		Ep_obs(i,3) = Ez
		Hp_obs(i,1) = Hx
		Hp_obs(i,2) = Hy
		Hp_obs(i,3) = Hz

	end do

elseif( FDP == 'HyDP' )then

	do i=1, N

		z = Mnos(no_obs(i),3)
		y = Mnos(no_obs(i),2)
		x = Mnos(no_obs(i),1)

		call dmhx_xyz_loops( autr_fltJ0, numpJ0, autr_fltJ1, numpJ1, Ty, -Tx, Iw, dsy, Tz, NC, &
							     hj, sigmas, neta0, zeta, y, -x, z, Ey, Ex, Ez , Hy, Hx, Hz )
		Ex = -Ex
		Hx = -Hx

		Ep_obs(i,1) = Ex
		Ep_obs(i,2) = Ey
		Ep_obs(i,3) = Ez
		Hp_obs(i,1) = Hx
		Hp_obs(i,2) = Hy
		Hp_obs(i,3) = Hz

	end do

elseif( FDP == 'HzDP' )then

	do i=1, N

		z = Mnos(no_obs(i),3)
		y = Mnos(no_obs(i),2)
		x = Mnos(no_obs(i),1)

		call dmv_xyz_loops( autr_fltJ0, numpJ0, autr_fltJ1, numpJ1, Tx, Ty, dsz, Tz, NC, &
							     hj, sigmas, neta0, zeta, x, y, z, Ex, Ey, Ez , Hx, Hy, Hz )

		Ep_obs(i,1) = Ex
		Ep_obs(i,2) = Ey
		Ep_obs(i,3) = Ez
		Hp_obs(i,1) = Hx
		Hp_obs(i,2) = Hy
		Hp_obs(i,3) = Hz

	end do

elseif( FDP == 'HCL' )then

	do i=1, N

		z = Mnos(no_obs(i),3)
		y = Mnos(no_obs(i),2)
		x = Mnos(no_obs(i),1)

		call BCH1DQWE_campo_E( valf, Tx, Ty, Tz, x, y, z, Ephi, Ex, Ey )
		call BCH1DQWE_campo_H( valf, Tx, Ty, Tz, x, y, z, Hr, Hx, Hy, Hz )

		Ep_obs(i,1) = Ex
		Ep_obs(i,2) = Ey
		Ep_obs(i,3) = Ez
		Hp_obs(i,1) = Hx
		Hp_obs(i,2) = Hy
		Hp_obs(i,3) = Hz

	end do

end if



do i=1, N

	Es_obs(i,1) = ( vsol(Marestas(eleobs(indobs(i)),4)) + &
			vsol(Marestas(eleobs(indobs(i)+1),4)) )/2.d0

	Es_obs(i,2) = ( vsol(Marestas(eleobs(indobs(i)),8)) + &
			vsol(Marestas(eleobs(indobs(i)+2),8)) )/2.d0

	Es_obs(i,3) = ( vsol(Marestas(eleobs(indobs(i)),12)) + &
			vsol(Marestas(eleobs(indobs(i)+4),12)) )/2.d0

	Et_obs(i,:) =   Ep_obs(i,:) + Es_obs(i,:)
	Ht_obs(i,:) =   Hp_obs(i,:) + Hs_obs(i,:)		

end do

END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE campo_Hs_obs( vsol, Hs_obs, w )
IMPLICIT NONE
COMPLEX(dp), DIMENSION(:,:), INTENT(out) :: Hs_obs
COMPLEX(dp), DIMENSION(:),   INTENT(in)	 :: vsol
REAL(dp),		     INTENT(in)  :: w	

COMPLEX(dp), ALLOCATABLE, DIMENSION(:)	:: dfxdy, dfxdz, dfydx, dfydz, dfzdx, dfzdy
REAL(dp)				:: lx(2), ly(2), lz(2)
INTEGER					:: i, N, aresta_x(10), aresta_y(10), aresta_z(10)

 N = Nobs

allocate( dfxdy(N), dfxdz(N), dfydx(N), dfydz(N), dfzdx(N), dfzdy(N) )

!=====================================================================================
! calculando as derivadas por diferença finita ( Razão incremental )
!=====================================================================================

do i=1, N


	lx(1) = ( Mnos(Mnosaresta(Marestas(eleobs(indobs(i)),1),2),1) - &
		  Mnos(Mnosaresta(Marestas(eleobs(indobs(i)),1),1),1) )! /2.d0

	lx(2) = ( Mnos(Mnosaresta(Marestas(eleobs(indobs(i)+1),1),2),1) - &
		  Mnos(Mnosaresta(Marestas(eleobs(indobs(i)+1),1),1),1) )! /2.d0

	ly(1) = ( Mnos(Mnosaresta(Marestas(eleobs(indobs(i)),5),2),2) - &
		  Mnos(Mnosaresta(Marestas(eleobs(indobs(i)),5),1),2) )! /2.d0

	ly(2) = ( Mnos(Mnosaresta(Marestas(eleobs(indobs(i)+2),5),2),2) - &
		  Mnos(Mnosaresta(Marestas(eleobs(indobs(i)+2),5),1),2) )! /2.d0

	lz(1) = ( Mnos(Mnosaresta(Marestas(eleobs(indobs(i)),9),2),3) - &
		  Mnos(Mnosaresta(Marestas(eleobs(indobs(i)),9),1),3) )! /2.d0

	lz(2) = ( Mnos(Mnosaresta(Marestas(eleobs(indobs(i)+4),9),2),3) - &
		  Mnos(Mnosaresta(Marestas(eleobs(indobs(i)+4),9),1),3) )! /2.d0



	aresta_x(1)  = Marestas(eleobs(indobs(i)),4)
	aresta_x(2)  = Marestas(eleobs(indobs(i)+1),4)
	aresta_x(3)  = Marestas(eleobs(indobs(i)),3)
	aresta_x(4)  = Marestas(eleobs(indobs(i)+2),4)
	aresta_x(5)  = Marestas(eleobs(indobs(i)+1),3)
	aresta_x(6)  = Marestas(eleobs(indobs(i)+3),4)
	aresta_x(7)  = Marestas(eleobs(indobs(i)),2)
	aresta_x(8)  = Marestas(eleobs(indobs(i)+4),4)
	aresta_x(9)  = Marestas(eleobs(indobs(i)+1),2)
	aresta_x(10) = Marestas(eleobs(indobs(i)+5),4)
	
	aresta_y(1)  = Marestas(eleobs(indobs(i)),8)
	aresta_y(2)  = Marestas(eleobs(indobs(i)+2),8)
	aresta_y(3)  = Marestas(eleobs(indobs(i)),6)
	aresta_y(4)  = Marestas(eleobs(indobs(i)+1),8)
	aresta_y(5)  = Marestas(eleobs(indobs(i)+2),6)
	aresta_y(6)  = Marestas(eleobs(indobs(i)+3),8)
	aresta_y(7)  = Marestas(eleobs(indobs(i)),7)
	aresta_y(8)  = Marestas(eleobs(indobs(i)+4),8)
	aresta_y(9)  = Marestas(eleobs(indobs(i)+2),7)
	aresta_y(10) = Marestas(eleobs(indobs(i)+6),8)

	aresta_z(1)  = Marestas(eleobs(indobs(i)),12)
	aresta_z(2)  = Marestas(eleobs(indobs(i)+4),12)
	aresta_z(3)  = Marestas(eleobs(indobs(i)),11)
	aresta_z(4)  = Marestas(eleobs(indobs(i)+1),12)
	aresta_z(5)  = Marestas(eleobs(indobs(i)+4),11)
	aresta_z(6)  = Marestas(eleobs(indobs(i)+5),12)
	aresta_z(7)  = Marestas(eleobs(indobs(i)),10)
	aresta_z(8)  = Marestas(eleobs(indobs(i)+2),12)
	aresta_z(9)  = Marestas(eleobs(indobs(i)+4),10)
	aresta_z(10) = Marestas(eleobs(indobs(i)+6),12)



	dfxdy(i) = ( (vsol(aresta_x(1)) - vsol(aresta_x(3)) )/ly(1) + &
			     (vsol(aresta_x(4)) - vsol(aresta_x(1)) )/ly(2) )/4.d0 + & 
			   ( (vsol(aresta_x(2)) - vsol(aresta_x(5)) )/ly(1) + &
			     (vsol(aresta_x(6)) - vsol(aresta_x(2)) )/ly(2) )/4.d0


	dfxdz(i) = ( (vsol(aresta_x(1)) - vsol(aresta_x(7)) )/lz(1) + &
			     (vsol(aresta_x(8)) - vsol(aresta_x(1)) )/lz(2) )/4.d0 + & 
			   ( (vsol(aresta_x(2)) - vsol(aresta_x(9)) )/lz(1) + &
			     (vsol(aresta_x(10)) - vsol(aresta_x(2)) )/lz(2) )/4.d0


	dfydx(i) = ( (vsol(aresta_y(1)) - vsol(aresta_y(3)) )/lx(1) + &
			     (vsol(aresta_y(4)) - vsol(aresta_y(1)) )/lx(2) )/4.d0 + & 
			   ( (vsol(aresta_y(2)) - vsol(aresta_y(5)) )/lx(1) + &
			     (vsol(aresta_y(6)) - vsol(aresta_y(2)) )/lx(2) )/4.d0 


	dfydz(i) = ( (vsol(aresta_y(1)) - vsol(aresta_y(7)) )/lz(1) + &
			     (vsol(aresta_y(8)) - vsol(aresta_y(1)) )/lz(2) )/4.d0 + & 
			   ( (vsol(aresta_y(2)) - vsol(aresta_y(9)) )/lz(1) + &
			     (vsol(aresta_y(10)) - vsol(aresta_y(2)) )/lz(2) )/4.d0 


	dfzdx(i) = ( (vsol(aresta_z(1)) - vsol(aresta_z(3)) )/lx(1) + &
			     (vsol(aresta_z(4)) - vsol(aresta_z(1)) )/lx(2) )/4.d0 + & 
			   ( (vsol(aresta_z(2)) - vsol(aresta_z(5)) )/lx(1) + &
			     (vsol(aresta_z(6)) - vsol(aresta_z(2)) )/lx(2) )/4.d0 
	

	dfzdy(i) = ( (vsol(aresta_z(1)) - vsol(aresta_z(7)) )/ly(1) + &
			     (vsol(aresta_z(8)) - vsol(aresta_z(1)) )/ly(2) )/4.d0 + & 
			   ( (vsol(aresta_z(2)) - vsol(aresta_z(9)) )/ly(1) + &
			     (vsol(aresta_z(10)) - vsol(aresta_z(2)) )/ly(2) )/4.d0 

end do

Hs_obs(:,1) = dfzdy - dfydz 
Hs_obs(:,2) = dfxdz - dfzdx
Hs_obs(:,3) = dfydx - dfxdy

Hs_obs = -1.d0*Hs_obs/ci/w/mi0

END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE rho_ap_fase( f)
IMPLICIT NONE
 REAL(dp), INTENT(in) :: f

 COMPLEX(dp)	:: DET(5)
 INTEGER	:: i, N
 REAL(dp)	:: w

 w = (pi+pi)*f

N = Nobs

allocate( rho_ap_xx(N), fase_xx(N), rho_ap_yy(N), fase_yy(N) )
allocate( rho_ap_xy(N), fase_xy(N), rho_ap_yx(N), fase_yx(N) )

rho_ap_xx = 0.d0	; fase_xx = 0.d0
rho_ap_yy = 0.d0	; fase_yy = 0.d0
rho_ap_xy = 0.d0	; fase_xy = 0.d0
rho_ap_yx = 0.d0	; fase_yx = 0.d0

do i=1, N

	DET(1)	= Ht_obs_1(i,1)*Ht_obs_2(i,2) - Ht_obs_1(i,2)*Ht_obs_2(i,1)

	DET(2)	= Et_obs_1(i,1)*Ht_obs_2(i,2) - Ht_obs_1(i,2)*Et_obs_2(i,1)
	DET(3)	= Ht_obs_1(i,1)*Et_obs_2(i,1) - Et_obs_1(i,1)*Ht_obs_2(i,1)
	DET(4)	= Et_obs_1(i,2)*Ht_obs_2(i,2) - Ht_obs_1(i,2)*Et_obs_2(i,2)
	DET(5)	= Ht_obs_1(i,1)*Et_obs_2(i,2) - Et_obs_1(i,2)*Ht_obs_2(i,1)

	rho_ap_xx(i) = cdabs( DET(2)/DET(1) )**2/w/mi0
	rho_ap_yy(i) = cdabs( DET(5)/DET(1) )**2/w/mi0

	rho_ap_xy(i) = cdabs( DET(3)/DET(1) )**2/w/mi0
	rho_ap_yx(i) = cdabs( DET(4)/DET(1) )**2/w/mi0

	fase_xx(i) = datan2(aimag(DET(2)/DET(1)), dreal(DET(2)/DET(1)))*180.d0/pi
	fase_yy(i) = datan2(aimag(DET(5)/DET(1)), dreal(DET(5)/DET(1)))*180.d0/pi

	fase_xy(i) = datan2(aimag(DET(3)/DET(1)), dreal(DET(3)/DET(1)))*180.d0/pi
	fase_yx(i) = datan2(aimag(DET(4)/DET(1)), dreal(DET(4)/DET(1)))*180.d0/pi

end do

END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE openRW_arq_MT
 character(len=80)  :: msgm
 integer            :: estts

	open( unit=1100, file=trim(pathr)//'/fields_MT/field_Es_1.dat', &
					status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
		IF( estts /= 0 )THEN
			PRINT*,TRIM(MSGM)
			STOP	
		END IF 
	open( unit=1200, file=trim(pathr)//'/fields_MT/field_Hs_1.dat', &
					status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
		IF( estts /= 0 )THEN
			PRINT*,TRIM(MSGM)
			STOP	
		END IF 
	open( unit=1300, file=trim(pathr)//'/fields_MT/field_Es_2.dat', &
					status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
		IF( estts /= 0 )THEN
			PRINT*,TRIM(MSGM)
			STOP	
		END IF 
	open( unit=1400, file=trim(pathr)//'/fields_MT/field_Hs_2.dat', &
					status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
		IF( estts /= 0 )THEN
			PRINT*,TRIM(MSGM)
			STOP	
		END IF 
	open( unit=1500, file=trim(pathr)//'/fields_MT/field_Et_1.dat', &
					status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
		IF( estts /= 0 )THEN
			PRINT*,TRIM(MSGM)
			STOP	
		END IF 
	open( unit=1600, file=trim(pathr)//'/fields_MT/field_Ht_1.dat', &
					status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
		IF( estts /= 0 )THEN
			PRINT*,TRIM(MSGM)
			STOP	
		END IF 
	open( unit=1700, file=trim(pathr)//'/fields_MT/field_Et_2.dat', &
					status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
		IF( estts /= 0 )THEN
			PRINT*,TRIM(MSGM)
			STOP	
		END IF 
	open( unit=1800, file=trim(pathr)//'/fields_MT/field_Ht_2.dat', &
					status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
		IF( estts /= 0 )THEN
			PRINT*,TRIM(MSGM)
			STOP	
		END IF 
	open( unit=1900, file=trim(pathr)//'/fields_MT/rho_ap_phse.dat', &
					status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
		IF( estts /= 0 )THEN
			PRINT*,TRIM(MSGM)
			STOP	
		END IF 

END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE close_arq_MT

	close( 1100 )
	close( 1200 )
	close( 1300 )
	close( 1400 )
	close( 1500 )
	close( 1600 )
	close( 1700 )
	close( 1800 )
	close( 1900 )

END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================



SUBROUTINE openRW_arq_CSAMT
 character(len=80)  :: msgm
 integer            :: estts

	open( unit=1110, file=trim(pathr)//'/fields_CSAMT/field_Es_1.dat', &
					status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
		IF( estts /= 0 )THEN
			PRINT*,TRIM(MSGM)
			STOP	
		END IF 
	open( unit=1210, file=trim(pathr)//'/fields_CSAMT/field_Hs_1.dat', &
					status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
		IF( estts /= 0 )THEN
			PRINT*,TRIM(MSGM)
			STOP	
		END IF
	open( unit=1310, file=trim(pathr)//'/fields_CSAMT/field_Es_2.dat', &
					status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
		IF( estts /= 0 )THEN
			PRINT*,TRIM(MSGM)
			STOP	
		END IF
	open( unit=1410, file=trim(pathr)//'/fields_CSAMT/field_Hs_2.dat', &
					status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
		IF( estts /= 0 )THEN
			PRINT*,TRIM(MSGM)
			STOP	
		END IF
	open( unit=1510, file=trim(pathr)//'/fields_CSAMT/field_Et_1.dat', &
					status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
		IF( estts /= 0 )THEN
			PRINT*,TRIM(MSGM)
			STOP	
		END IF
	open( unit=1610, file=trim(pathr)//'/fields_CSAMT/field_Ht_1.dat', &
					status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
		IF( estts /= 0 )THEN
			PRINT*,TRIM(MSGM)
			STOP	
		END IF
	open( unit=1710, file=trim(pathr)//'/fields_CSAMT/field_Et_2.dat', &
					status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
		IF( estts /= 0 )THEN
			PRINT*,TRIM(MSGM)
			STOP	
		END IF
	open( unit=1810, file=trim(pathr)//'/fields_CSAMT/field_Ht_2.dat', &
					status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
		IF( estts /= 0 )THEN
			PRINT*,TRIM(MSGM)
			STOP	
		END IF
	open( unit=1910, file=trim(pathr)//'/fields_CSAMT/rho_ap_phase.dat', &
					status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
		IF( estts /= 0 )THEN
			PRINT*,TRIM(MSGM)
			STOP	
		END IF

END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE close_arq_CSAMT

	close( 1110 )
	close( 1210 )
	close( 1310 )
	close( 1410 )
	close( 1510 )
	close( 1610 )
	close( 1710 )
	close( 1810 )
	close( 1910 )

END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE openRW_arq_DPEH( DPEH )
IMPLICIT NONE
 CHARACTER(LEN=20), INTENT(in)		:: DPEH

 character(len=80)  :: msgm
 integer            :: estts

	if( DPEH == 'ExDP' )then


		open( unit=2100, file=trim(pathr)//'fields_DPEM/ExDP/field_Es_DIP.dat', &
							status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
		IF( estts /= 0 )THEN
			PRINT*,TRIM(MSGM)
			STOP	
		END IF

		open( unit=2200, file=trim(pathr)//'fields_DPEM/ExDP/field_Hs_DIP.dat', &
							status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
		IF( estts /= 0 )THEN
			PRINT*,TRIM(MSGM)
			STOP	
		END IF

		open( unit=2300, file=trim(pathr)//'fields_DPEM/ExDP/field_Et_DIP.dat', &
							status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
		IF( estts /= 0 )THEN
			PRINT*,TRIM(MSGM)
			STOP	
		END IF

		open( unit=2400, file=trim(pathr)//'fields_DPEM/ExDP/field_Ht_DIP.dat', &
							status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
		IF( estts /= 0 )THEN
			PRINT*,TRIM(MSGM)
			STOP	
		END IF


	elseif( DPEH == 'EyDP' )then


		open( unit=2100, file=trim(pathr)//'/fields_DPEM/EyDP/field_Es_DIP.dat', &
							status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
		IF( estts /= 0 )THEN
			PRINT*,TRIM(MSGM)
			STOP	
		END IF

		open( unit=2200, file=trim(pathr)//'/fields_DPEM/EyDP/field_Hs_DIP.dat', &
							status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
		IF( estts /= 0 )THEN
			PRINT*,TRIM(MSGM)
			STOP	
		END IF

		open( unit=2300, file=trim(pathr)//'/fields_DPEM/EyDP/field_Et_DIP.dat', &
							status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
		IF( estts /= 0 )THEN
			PRINT*,TRIM(MSGM)
			STOP	
		END IF

		open( unit=2400, file=trim(pathr)//'/fields_DPEM/EyDP/field_Ht_DIP.dat', &
							status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
		IF( estts /= 0 )THEN
			PRINT*,TRIM(MSGM)
			STOP	
		END IF

	elseif( DPEH == 'EzDP' )then


		open( unit=2100, file=trim(pathr)//'/fields_DPEM/EzDP/field_Es_DIP.dat', &
							status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
		IF( estts /= 0 )THEN
			PRINT*,TRIM(MSGM)
			STOP	
		END IF

		open( unit=2200, file=trim(pathr)//'/fields_DPEM/EzDP/field_Hs_DIP.dat', &
							status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
		IF( estts /= 0 )THEN
			PRINT*,TRIM(MSGM)
			STOP	
		END IF

		open( unit=2300, file=trim(pathr)//'/fields_DPEM/EzDP/field_Et_DIP.dat', &
							status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
		IF( estts /= 0 )THEN
			PRINT*,TRIM(MSGM)
			STOP	
		END IF

		open( unit=2400, file=trim(pathr)//'/fields_DPEM/EzDP/field_Ht_DIP.dat', &
							status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
		IF( estts /= 0 )THEN
			PRINT*,TRIM(MSGM)
			READ(*,*)	
		END IF


	elseif( DPEH == 'HxDP' )then


		open( unit=2100, file=trim(pathr)//'/fields_DPEM/DPHx/campo_Es_DIP.dat', &
							status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
		IF( estts /= 0 )THEN
			PRINT*,TRIM(MSGM)
			STOP	
		END IF

		open( unit=2200, file=trim(pathr)//'/fields_DPEM/DPHx/campo_Hs_DIP.dat', &
							status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
		IF( estts /= 0 )THEN
			PRINT*,TRIM(MSGM)
			STOP	
		END IF

		open( unit=2300, file=trim(pathr)//'/fields_DPEM/DPHx/campo_Et_DIP.dat', &
							status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
		IF( estts /= 0 )THEN
			PRINT*,TRIM(MSGM)
			STOP	
		END IF

		open( unit=2400, file=trim(pathr)//'/fields_DPEM/DPHx/campo_Ht_DIP.dat', &
							status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
		IF( estts /= 0 )THEN
			PRINT*,TRIM(MSGM)
			STOP	
		END IF


	elseif( DPEH == 'HyDP' )then


		open( unit=2100, file=trim(pathr)//'/fields_DPEM/DPHy/campo_Es_DIP.dat', &
							status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
		IF( estts /= 0 )THEN
			PRINT*,TRIM(MSGM)
			STOP	
		END IF

		open( unit=2200, file=trim(pathr)//'/fields_DPEM/DPHy/campo_Hs_DIP.dat', &
							status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
		IF( estts /= 0 )THEN
			PRINT*,TRIM(MSGM)
			STOP	
		END IF

		open( unit=2300, file=trim(pathr)//'/fields_DPEM/DPHy/campo_Et_DIP.dat', &
							status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
		IF( estts /= 0 )THEN
			PRINT*,TRIM(MSGM)
			STOP	
		END IF

		open( unit=2400, file=trim(pathr)//'/fields_DPEM/DPHy/campo_Ht_DIP.dat', &
							status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
		IF( estts /= 0 )THEN
			PRINT*,TRIM(MSGM)
			STOP	
		END IF


	elseif( DPEH == 'HzDP' )then


		open( unit=2100, file=trim(pathr)//'/fields_DPEM/DPHz/campo_Es_DIP.dat', &
							status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm ) 
		IF( estts /= 0 )THEN
			PRINT*,TRIM(MSGM)
			STOP	
		END IF

		open( unit=2200, file=trim(pathr)//'/fields_DPEM/DPHz/campo_Hs_DIP.dat', &
							status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm ) 
		IF( estts /= 0 )THEN
			PRINT*,TRIM(MSGM)
			STOP	
		END IF

		open( unit=2300, file=trim(pathr)//'/fields_DPEM/DPHz/campo_Et_DIP.dat', &
							status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
		IF( estts /= 0 )THEN
			PRINT*,TRIM(MSGM)
			STOP	
		END IF

		open( unit=2400, file=trim(pathr)//'/fields_DPEM/DPHz/campo_Ht_DIP.dat', &
							status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
		IF( estts /= 0 )THEN
			PRINT*,TRIM(MSGM)
			STOP	
		END IF

	end if


END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE close_arq_DPEH

	close( 2100 )
	close( 2200 )
	close( 2300 )
	close( 2400 )
 
END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE openRW_arq_BCH
 character(len=80)  :: msgm
 integer            :: estts

	open( unit=2110, file=trim(pathr)//'fields_DPEM/HCL/fields_Es_HCL.dat', &
			                     status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
	IF( estts /= 0 )THEN
		PRINT*,TRIM(MSGM)
		STOP	
	END IF


	open( unit=2210, file=trim(pathr)//'fields_DPEM/HCL/field_Hs_HCL.dat', &
			                     status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
	IF( estts /= 0 )THEN
		PRINT*,TRIM(MSGM)
		STOP	
	END IF


	open( unit=2310, file=trim(pathr)//'fields_DPEM/HCL/field_Et_HCL.dat', &
			                     status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
	IF( estts /= 0 )THEN
		PRINT*,TRIM(MSGM)
		STOP	
	END IF


	open( unit=2410, file=trim(pathr)//'fields_DPEM/HCL/field_Ht_HCL.dat', &
			                     status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
	IF( estts /= 0 )THEN
		PRINT*,TRIM(MSGM)
		STOP	
	END IF


END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================


SUBROUTINE close_arq_BCH

	close( 2110 )
	close( 2210 )
	close( 2310 )
	close( 2410 )
 
END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE openRW_arq_1D_adj_MT
IMPLICIT NONE
 character(len=20)  :: doc1, doc2, doc3, doc4
 character(len=300) :: path
 character(len=80)  :: msgm
 integer            :: estts


	path = '1D_adj_sourc/MT'
	path = trim(pathr)//trim(path)

	doc1  = 'Ep_adj_dpREx.dat' 
	open( unit=3000, file=doc1, DEFAULTFILE=path, &
			 form='unformatted', status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
	IF( estts /= 0 )THEN
		PRINT*,TRIM(MSGM)
		STOP	
	END IF


	doc2 = 'Ep_adj_dpREy.dat' 
	open( unit=3100, file=doc2, DEFAULTFILE=path, &
			 form='unformatted', status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
	IF( estts /= 0 )THEN
		PRINT*,TRIM(MSGM)
		STOP	
	END IF


	doc3 = 'Ep_adj_dpRHx.dat' 
	open( unit=3300, file=doc3, DEFAULTFILE=path, &
			 form='unformatted', status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
	IF( estts /= 0 )THEN
		PRINT*,TRIM(MSGM)
		STOP	
	END IF


	doc4 = 'Ep_adj_dpRHy.dat' 
	open( unit=3400, file=doc4, DEFAULTFILE=path, &
			 form='unformatted', status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
	IF( estts /= 0 )THEN
		PRINT*,TRIM(MSGM)
		STOP	
	END IF



END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE close_arq_1D_adj_MT

	close( 3000 )
	close( 3100 )
	close( 3300 )
	close( 3400 )

END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE openRW_arq_1D_adj_CSAMT
IMPLICIT NONE
 character(len=20)  :: doc1, doc2, doc3, doc4
 character(len=300) :: path
 character(len=80)  :: msgm
 integer            :: estts


	path = '1D_adj_sourc/CSAMT'
	path = trim(pathr)//trim(path)


	doc1  = 'Ep_adj_dpREx.dat' 
	open( unit=3000, file=doc1, DEFAULTFILE=path, &
			 form='unformatted', status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
	IF( estts /= 0 )THEN
		PRINT*,TRIM(MSGM)
		STOP	
	END IF


	doc2 = 'Ep_adj_dpREy.dat' 
	open( unit=3100, file=doc2, DEFAULTFILE=path, &
			 form='unformatted', status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
	IF( estts /= 0 )THEN
		PRINT*,TRIM(MSGM)
		STOP	
	END IF


	doc3 = 'Ep_adj_dpRHx.dat' 
	open( unit=3300, file=doc3, DEFAULTFILE=path, &
			 form='unformatted', status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
	IF( estts /= 0 )THEN
		PRINT*,TRIM(MSGM)
		STOP	
	END IF


	doc4 = 'Ep_adj_dpRHy.dat' 
	open( unit=3400, file=doc4, DEFAULTFILE=path, &
			 form='unformatted', status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )

	IF( estts /= 0 )THEN
		PRINT*,TRIM(MSGM)
		STOP
	END IF


END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE close_arq_1D_adj_CSAMT

	close( 3000 )
	close( 3100 )
	close( 3300 )
	close( 3400 )

END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE openRW_arq_DPEM1D_adj( DPEH )
IMPLICIT NONE
 CHARACTER(LEN=20), INTENT(in)		:: DPEH

 character(len=20)  :: doc
 character(len=300) :: path
 character(len=80)  :: msgm
 integer            :: unid, estts

 

	path = '1D_adj_sourc/DPEM'
	path = trim(pathr)//trim(path)

	if( DPEH == 'ExDP' )then

		doc  = 'Ep_adj_dpREx.dat'
		unid = 3010 

	elseif( DPEH == 'EyDP' )then

		doc  = 'Ep_adj_dpREy.dat' 
		unid = 3110

	elseif( DPEH == 'EzDP' )then

		doc  = 'Ep_adj_dpREz.dat' 
		unid = 3210

	elseif( DPEH == 'HxDP' )then

		doc  = 'Ep_adj_dpRHx.dat' 
		unid = 3310

	elseif( DPEH == 'HyDP' )then

		doc  = 'Ep_adj_dpRHy.dat' 
		unid = 3410
	elseif( DPEH == 'HzDP' )then

		doc  = 'Ep_adj_dpRHz.dat' 
		unid = 3510
	end if

	open( unit=unid, file=doc, DEFAULTFILE=path, &
			 form='unformatted', status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )

	IF( estts /= 0 )THEN
		PRINT*,TRIM(MSGM)
		STOP	
	END IF


END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================


SUBROUTINE close_arq_DPEM1D_adj( DPEH )
IMPLICIT NONE
 CHARACTER(LEN=20), INTENT(in)		:: DPEH

	if( DPEH == 'ExDP' )then
		close( 3010 )
	elseif( DPEH == 'EyDP' )then
		close( 3110 )
	elseif( DPEH == 'EzDP' )then
		close( 3210 )
	elseif( DPEH == 'HxDP' )then
		close( 3310 )	
	elseif( DPEH == 'HyDP' )then
		close( 3410 )
	elseif( DPEH == 'HzDP' )then
		close( 3510 )
	end if	

END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE openRW_arq_DPEM1D( DPEH )
IMPLICIT NONE
 CHARACTER(LEN=20), INTENT(in)		:: DPEH
 
 character(len=20)  :: doc
 character(len=80)  :: msgm
 character(len=300) :: path
 integer            :: unid, estts

	path = '1D_orig_sourc'
	path = trim(pathr)//trim(path)

	if( DPEH == 'ExDP' )then

		doc  = 'Ep_adj_dpTEx.dat' 
		unid  = 4000

	elseif( DPEH == 'EyDP' )then

		doc = 'Ep_adj_dpTEy.dat' 
		unid  = 4100

	elseif( DPEH == 'EzDP' )then

		doc = 'Ep_adj_dpTEz.dat' 
		unid  = 4200

	elseif( DPEH == 'HxDP' )then

		doc = 'Ep_adj_dpTHx.dat' 
		unid  = 4300

	elseif( DPEH == 'HyDP' )then

		doc = 'Ep_adj_dpTHy.dat' 
		unid  = 4400

	elseif( DPEH == 'HzDP' )then

		doc = 'Ep_adj_dpTHz.dat' 
		unid  = 4500

	elseif( DPEH == 'HCL' )then

		doc = 'Ep_adj_BCH.dat' 
		unid  = 4600

	end if

	open( unit=unid, file=doc, DEFAULTFILE=path, &
			 form='unformatted', status='old', action='readwrite', IOSTAT=estts, IOMSG=msgm )
	IF( estts /= 0 )THEN
		PRINT*,TRIM(MSGM)
		STOP	
	END IF

END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================


SUBROUTINE close_arq_DPEM1D( DPEH )
IMPLICIT NONE
 CHARACTER(LEN=20), INTENT(in)		:: DPEH

	if( DPEH == 'ExDP' )then
		close( 4000 )
	elseif( DPEH == 'EyDP' )then
		close( 4100 )
	elseif( DPEH == 'EzDP' )then
		close( 4200 )
	elseif( DPEH == 'HxDP' )then
		close( 4300 )	
	elseif( DPEH == 'HyDP' )then
		close( 4400 )
	elseif( DPEH == 'HzDP' )then
		close( 4500 )
	elseif( DPEH == 'HCL' )then
		close( 4600 )
	end if	

END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================


SUBROUTINE open_arq_sensib_MT		
IMPLICIT NONE
 character(len=20) :: doc1, doc2
 character(len=300) :: path
 character(len=80)  :: msgm
 integer            :: estts

	path = trim(pathr)//'sensib_MT'

	doc1 = 'msb_rho_adj.dat'
	open( unit=3600, file=doc1, DEFAULTFILE=path, status='old', action='write', IOSTAT=estts, IOMSG=msgm )
	IF( estts /= 0 )THEN
		PRINT*,TRIM(MSGM)
		STOP	
	END IF 


	doc2 = 'msb_phase_adj.dat'
	open( unit=3700, file=doc2, DEFAULTFILE=path, status='old', action='write', IOSTAT=estts, IOMSG=msgm )
	IF( estts /= 0 )THEN
		PRINT*,TRIM(MSGM)
		STOP	
	END IF

END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE open_arq_sensib_CSAMT	
IMPLICIT NONE
 character(len=20) :: doc1, doc2
 character(len=300) :: path
 character(len=80)  :: msgm
 integer            :: estts

	path = trim(pathr)//'sensib_CSAMT'

	doc1 = 'msb_rho_adj.dat'
	open( unit=3600, file=doc1, DEFAULTFILE=path, status='old', action='write', IOSTAT=estts, IOMSG=msgm )
	IF( estts /= 0 )THEN
		PRINT*,TRIM(MSGM)
		STOP	
	END IF

	doc2 = 'msb_phase_adj.dat'
	open( unit=3700, file=doc2, DEFAULTFILE=path, status='old', action='write', IOSTAT=estts, IOMSG=msgm )
	IF( estts /= 0 )THEN
		PRINT*,TRIM(MSGM)
		STOP	
	END IF

END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE open_arq_sensib_DPEM( DPEH )
IMPLICIT NONE
 CHARACTER(LEN=20), INTENT(in)		:: DPEH		

 character(len=20)  :: doc1, doc2
 character(len=300) :: path
 character(len=80)  :: msgm
 integer            :: estts

	if( DPEH == 'ExDP' )then
		path = trim(pathr)//'sensib_DPEM/ExDP'
	elseif( DPEH == 'EyDP' )then
		path = trim(pathr)//'sensib_DPEM/ExDP'
	elseif( DPEH == 'EzDP' )then
		path = trim(pathr)//'sensib_DPEM/EzDP'
	elseif( DPEH == 'HxDP' )then
		path = trim(pathr)//'sensib_DPEM/HxDP'
	elseif( DPEH == 'HyDP' )then
		path = trim(pathr)//'sensib_DPEM/HyDP'
	elseif( DPEH == 'HzDP' )then
		path = trim(pathr)//'sensib_DPEM/HzDP'
	elseif( DPEH == 'HCL' )then
		path = trim(pathr)//'sensib_DPEM/HCL'
	end if

	doc1 = 'msb_EH_adj.dat'
	doc2 = 'msb_fase_adj.dat'

	open( unit=3610, file=doc1, DEFAULTFILE=path, status='old', action='write', IOSTAT=estts, IOMSG=msgm )
	IF( estts /= 0 )THEN
		PRINT*,TRIM(MSGM)
		STOP	
	END IF
	open( unit=3710, file=doc2, DEFAULTFILE=path, status='old', action='write', IOSTAT=estts, IOMSG=msgm )
	IF( estts /= 0 )THEN
		PRINT*,TRIM(MSGM)
		STOP	
	END IF	


END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE escrever_resultados_MT3D( Es_obs_1, Hs_obs_1, Es_obs_2, Hs_obs_2,  Et_obs_1, Ht_obs_1, Et_obs_2, Ht_obs_2, &
				    rho_ap_xx, rho_ap_yy, rho_ap_xy, rho_ap_yx, fase_xx, fase_yy, fase_xy, fase_yx, N )
IMPLICIT NONE
 COMPLEX(dp), DIMENSION(:,:), INTENT(in) :: Es_obs_1, Hs_obs_1, Es_obs_2, Hs_obs_2
 COMPLEX(dp), DIMENSION(:,:), INTENT(in) :: Et_obs_1, Ht_obs_1, Et_obs_2, Ht_obs_2
 REAL(dp)   , DIMENSION(:)  , INTENT(in) :: rho_ap_xx, rho_ap_yy, rho_ap_xy, rho_ap_yx
 REAL(dp)   , DIMENSION(:)  , INTENT(in) :: fase_xx, fase_yy, fase_xy, fase_yx
 INTEGER    ,                 INTENT(in) :: N

 INTEGER	:: i

40 format( 2x, 6d25.16 )
60 format( 2x, 8d25.16 )

do i=1, N
	write(1100,40) dreal(Es_obs_1(i,1)), aimag(Es_obs_1(i,1)), dreal(Es_obs_1(i,2)), &
			aimag(Es_obs_1(i,2)), dreal(Es_obs_1(i,3)), aimag(Es_obs_1(i,3))
	write(1200,40) dreal(Hs_obs_1(i,1)), aimag(Hs_obs_1(i,1)), dreal(Hs_obs_1(i,2)), &
			aimag(Hs_obs_1(i,2)), dreal(Hs_obs_1(i,3)), aimag(Hs_obs_1(i,3))
	write(1300,40) dreal(Es_obs_2(i,1)), aimag(Es_obs_2(i,1)), dreal(Es_obs_2(i,2)), &
			aimag(Es_obs_2(i,2)), dreal(Es_obs_2(i,3)), aimag(Es_obs_2(i,3))
	write(1400,40) dreal(Hs_obs_2(i,1)), aimag(Hs_obs_2(i,1)), dreal(Hs_obs_2(i,2)), &
			aimag(Hs_obs_2(i,2)), dreal(Hs_obs_2(i,3)), aimag(Hs_obs_2(i,3))

	write(1500,40) dreal(Et_obs_1(i,1)), aimag(Et_obs_1(i,1)), dreal(Et_obs_1(i,2)), &
			aimag(Et_obs_1(i,2)), dreal(Et_obs_1(i,3)), aimag(Et_obs_1(i,3))
	write(1600,40) dreal(Ht_obs_1(i,1)), aimag(Ht_obs_1(i,1)), dreal(Ht_obs_1(i,2)), &
			aimag(Ht_obs_1(i,2)), dreal(Ht_obs_1(i,3)), aimag(Ht_obs_1(i,3))
	write(1700,40) dreal(Et_obs_2(i,1)), aimag(Et_obs_2(i,1)), dreal(Et_obs_2(i,2)), &
			aimag(Et_obs_2(i,2)), dreal(Et_obs_2(i,3)), aimag(Et_obs_2(i,3))
	write(1800,40) dreal(Ht_obs_2(i,1)), aimag(Ht_obs_2(i,1)), dreal(Ht_obs_2(i,2)), &
			aimag(Ht_obs_2(i,2)), dreal(Ht_obs_2(i,3)), aimag(Ht_obs_2(i,3))

end do

do i=1, N
	write(1900,60)  rho_ap_xx(i), fase_xx(i), rho_ap_yy(i), fase_yy(i), &
			rho_ap_xy(i), fase_xy(i), rho_ap_yx(i), fase_yx(i)
end do

END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE escrever_resultados_CSAMT3D( Es_obs_1, Hs_obs_1, Es_obs_2, Hs_obs_2,  Et_obs_1, Ht_obs_1, Et_obs_2, Ht_obs_2, &
				        rho_ap_xx, rho_ap_yy, rho_ap_xy, rho_ap_yx, fase_xx, fase_yy, fase_xy, fase_yx, N )
IMPLICIT NONE
 COMPLEX(dp), DIMENSION(:,:), INTENT(in) :: Es_obs_1, Hs_obs_1, Es_obs_2, Hs_obs_2
 COMPLEX(dp), DIMENSION(:,:), INTENT(in) :: Et_obs_1, Ht_obs_1, Et_obs_2, Ht_obs_2
 REAL(dp)   , DIMENSION(:)  , INTENT(in) :: rho_ap_xx, rho_ap_yy, rho_ap_xy, rho_ap_yx
 REAL(dp)   , DIMENSION(:)  , INTENT(in) :: fase_xx, fase_yy, fase_xy, fase_yx
 INTEGER    ,                 INTENT(in) :: N

 INTEGER	:: i

40 format( 2x, 6d25.16 )
60 format( 2x, 8d25.16 )

do i=1, N
	write(1110,40) dreal(Es_obs_1(i,1)), aimag(Es_obs_1(i,1)), dreal(Es_obs_1(i,2)), &
			aimag(Es_obs_1(i,2)), dreal(Es_obs_1(i,3)), aimag(Es_obs_1(i,3))
	write(1210,40) dreal(Hs_obs_1(i,1)), aimag(Hs_obs_1(i,1)), dreal(Hs_obs_1(i,2)), &
			aimag(Hs_obs_1(i,2)), dreal(Hs_obs_1(i,3)), aimag(Hs_obs_1(i,3))
	write(1310,40) dreal(Es_obs_2(i,1)), aimag(Es_obs_2(i,1)), dreal(Es_obs_2(i,2)), &
			aimag(Es_obs_2(i,2)), dreal(Es_obs_2(i,3)), aimag(Es_obs_2(i,3))
	write(1410,40) dreal(Hs_obs_2(i,1)), aimag(Hs_obs_2(i,1)), dreal(Hs_obs_2(i,2)), &
			aimag(Hs_obs_2(i,2)), dreal(Hs_obs_2(i,3)), aimag(Hs_obs_2(i,3))

	write(1510,40) dreal(Et_obs_1(i,1)), aimag(Et_obs_1(i,1)), dreal(Et_obs_1(i,2)), &
			aimag(Et_obs_1(i,2)), dreal(Et_obs_1(i,3)), aimag(Et_obs_1(i,3))
	write(1610,40) dreal(Ht_obs_1(i,1)), aimag(Ht_obs_1(i,1)), dreal(Ht_obs_1(i,2)), &
			aimag(Ht_obs_1(i,2)), dreal(Ht_obs_1(i,3)), aimag(Ht_obs_1(i,3))
	write(1710,40) dreal(Et_obs_2(i,1)), aimag(Et_obs_2(i,1)), dreal(Et_obs_2(i,2)), &
			aimag(Et_obs_2(i,2)), dreal(Et_obs_2(i,3)), aimag(Et_obs_2(i,3))
	write(1810,40) dreal(Ht_obs_2(i,1)), aimag(Ht_obs_2(i,1)), dreal(Ht_obs_2(i,2)), &
			aimag(Ht_obs_2(i,2)), dreal(Ht_obs_2(i,3)), aimag(Ht_obs_2(i,3))

end do

do i=1, N
	write(1910,60)  rho_ap_xx(i), fase_xx(i), rho_ap_yy(i), fase_yy(i), &
			rho_ap_xy(i), fase_xy(i), rho_ap_yx(i), fase_yx(i)
end do

END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE escrever_resultados_DPEH( Es_obs_1, Hs_obs_1, Et_obs_1, Ht_obs_1, N )
IMPLICIT NONE
 COMPLEX(dp), DIMENSION(:,:), INTENT(in) :: Es_obs_1, Hs_obs_1, Et_obs_1, Ht_obs_1
 INTEGER, 		     		  INTENT(in) :: N

 INTEGER	:: i

40 format( 2x, 6d25.16 )
60 format( 2x, 8d25.16 )

do i=1, N
	write(2100,40) dreal(Es_obs_1(i,1)), aimag(Es_obs_1(i,1)), dreal(Es_obs_1(i,2)), &
			aimag(Es_obs_1(i,2)), dreal(Es_obs_1(i,3)), aimag(Es_obs_1(i,3))
	write(2200,40) dreal(Hs_obs_1(i,1)), aimag(Hs_obs_1(i,1)), dreal(Hs_obs_1(i,2)), &
			aimag(Hs_obs_1(i,2)), dreal(Hs_obs_1(i,3)), aimag(Hs_obs_1(i,3))

	write(2300,40) dreal(Et_obs_1(i,1)), aimag(Et_obs_1(i,1)), dreal(Et_obs_1(i,2)), &
			aimag(Et_obs_1(i,2)), dreal(Et_obs_1(i,3)), aimag(Et_obs_1(i,3))
	write(2400,40) dreal(Ht_obs_1(i,1)), aimag(Ht_obs_1(i,1)), dreal(Ht_obs_1(i,2)), &
			aimag(Ht_obs_1(i,2)), dreal(Ht_obs_1(i,3)), aimag(Ht_obs_1(i,3))

end do

END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE escrever_resultados_BCH3D( Es_obs_1, Hs_obs_1, Et_obs_1, Ht_obs_1, N )
IMPLICIT NONE
COMPLEX(dp), DIMENSION(:,:), INTENT(in)	:: Es_obs_1, Hs_obs_1, Et_obs_1, Ht_obs_1
INTEGER,					 INTENT(in)	:: N

INTEGER	:: i

40 format( 2x, 6d25.16 )
60 format( 2x, 8d25.16 )

do i=1, N
	write(2110,40) dreal(Es_obs_1(i,1)), aimag(Es_obs_1(i,1)), dreal(Es_obs_1(i,2)), &
			aimag(Es_obs_1(i,2)), dreal(Es_obs_1(i,3)), aimag(Es_obs_1(i,3))
	write(2210,40) dreal(Hs_obs_1(i,1)), aimag(Hs_obs_1(i,1)), dreal(Hs_obs_1(i,2)), &
			aimag(Hs_obs_1(i,2)), dreal(Hs_obs_1(i,3)), aimag(Hs_obs_1(i,3))

	write(2310,40) dreal(Et_obs_1(i,1)), aimag(Et_obs_1(i,1)), dreal(Et_obs_1(i,2)), &
			aimag(Et_obs_1(i,2)), dreal(Et_obs_1(i,3)), aimag(Et_obs_1(i,3))
	write(2410,40) dreal(Ht_obs_1(i,1)), aimag(Ht_obs_1(i,1)), dreal(Ht_obs_1(i,2)), &
			aimag(Ht_obs_1(i,2)), dreal(Ht_obs_1(i,3)), aimag(Ht_obs_1(i,3))

end do

END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

END MODULE dados_saida
