MODULE mdlgm2_1D
!================================================================================================================

	USE ifport 
	USE variaveis_G

!================================================================================================================

        REAL(dp), PRIVATE :: zsup, zinf, h0, a, cz, w, r, IS
        INTEGER , PRIVATE :: kmada

!================================================================================================================

CONTAINS

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE BCH1DQWE_campo_E( valf, Tx, Ty, Tz, x, y, z, Ephip, Ex_p, Eyp )
IMPLICIT NONE
COMPLEX(dp), INTENT(out)	:: Ephip, Ex_p, Eyp
REAL(dp), INTENT(in)		:: valf, Tx, Ty, Tz, x, y, z

INTEGER	:: meio, ERRO

 r  = dsqrt((Tx-x)**2+(Ty-y)**2)
 if( r <= 1.d-10)then
  r= 1.d-1
 end if

 IS = ISOURCE
 a  = SDIMENS
 h0 = Tz
 cz = z
 w  = (pi+pi)*valf
 Ephip = (0.d0,0.d0)

call regiao_coord_z( meio, z, kmada )

if( meio == 1 )then

	!===========================================================================================
          call algoritmoepsilon_zerosreais( kernelephi_1, zrs_J1, tolrel, tolabs, xGssLg, wGssLg, Ephip , ERRO )
	!===========================================================================================

else if( meio == 2 )then  

	!===========================================================================================    
          call algoritmoepsilon_zerosreais( kernelephi_2, zrs_J1, tolrel, tolabs, xGssLg, wGssLg, Ephip , ERRO )
	!===========================================================================================	 				

else if( meio == 3 )then  

	!===========================================================================================
          call algoritmoepsilon_zerosreais( kernelephi_3, zrs_J1, tolrel, tolabs, xGssLg, wGssLg, Ephip , ERRO )
	!===========================================================================================

else if( meio == 4 )then  

	!===========================================================================================
         call algoritmoepsilon_zerosreais( kernelephi_4, zrs_J1, tolrel, tolabs, xGssLg, wGssLg, Ephip , ERRO )
	!===========================================================================================			

end if


Ex_p = -(y-Ty)*Ephip/r
Eyp   = (x-Tx)*Ephip/r

IF( ERRO == 1 )THEN
      print*,'=================================================================='
      print*,'The epsilon algorithm did not achieve the desired convergence with'
      print*,'number of intervals and pre-set quadrature points'
      print*,'=================================================================='
      STOP
END IF

END SUBROUTINE BCH1DQWE_campo_E

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================
 
SUBROUTINE BCH1DQWE_campo_H( valf, Tx, Ty, Tz, x, y, z, Hrp, Hxp, Hyp, Hzp )
IMPLICIT NONE
COMPLEX(dp), INTENT(out)	::  Hrp, Hxp, Hyp, Hzp
REAL(dp), INTENT(in)		::  valf, Tx, Ty, Tz, x, y, z

INTEGER	 :: meio, ERRO

 r  = dsqrt((Tx-x)**2+(Ty-y)**2)
 if( r <= 1.d-10)then
  r= 1.d-1
 end if

 IS = ISOURCE
 a  = SDIMENS
 h0 = Tz
 cz = z
 w  = (pi+pi)*valf

call regiao_coord_z( meio, z, kmada )

if( meio == 1 )then

	!===========================================================================================
	  call algoritmoepsilon_zerosreais( kernelhr_1, zrs_J1, tolrel, tolabs, xGssLg, wGssLg, Hrp , ERRO )
	  call algoritmoepsilon_zerosreais( kernelhz_1, zrs_J0, tolrel, tolabs, xGssLg, wGssLg, Hzp , ERRO )
	!===========================================================================================

else if( meio == 2 )then  

	!===========================================================================================    
	  call algoritmoepsilon_zerosreais( kernelhr_2, zrs_J1, tolrel, tolabs, xGssLg, wGssLg, Hrp , ERRO )
	  call algoritmoepsilon_zerosreais( kernelhz_2, zrs_J0, tolrel, tolabs, xGssLg, wGssLg, Hzp , ERRO )
	!==========================================================================================	 				

else if( meio == 3 )then  


	!===========================================================================================
	  call algoritmoepsilon_zerosreais( kernelhr_3, zrs_J1, tolrel, tolabs, xGssLg, wGssLg, Hrp , ERRO )
	  call algoritmoepsilon_zerosreais( kernelhz_3, zrs_J0, tolrel, tolabs, xGssLg, wGssLg, Hzp , ERRO )
	!===========================================================================================

else if( meio == 4 )then  

	!===========================================================================================
	  call algoritmoepsilon_zerosreais( kernelhr_4, zrs_J1, tolrel, tolabs, xGssLg, wGssLg, Hrp , ERRO )
	  call algoritmoepsilon_zerosreais( kernelhz_4, zrs_J0, tolrel, tolabs, xGssLg, wGssLg, Hzp , ERRO )
	!===========================================================================================			!

end if


Hxp = (x-Tx)*Hrp/r
Hyp = (y-Ty)*Hrp/r
Hzp = -Hzp/ci/w/mi0

IF( ERRO == 1 )THEN
      print*,'=================================================================='
      print*,'The epsilon algorithm did not achieve the desired convergence with'
      print*,'number of intervals and pre-set quadrature points'
      print*,'=================================================================='
      STOP
END IF
			
END SUBROUTINE BCH1DQWE_campo_H

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE BCH1DQWE_pll_campo_E( nprocss, valf, Tx, Ty, Tz, x, y, z, Ephip, Ex_p, Eyp )
IMPLICIT NONE
COMPLEX(dp), INTENT(out)    :: Ephip, Ex_p, Eyp
REAL(dp), INTENT(in)	    :: valf, Tx, Ty, Tz, x, y, z
INTEGER , INTENT(in)        :: nprocss 

INTEGER	:: meio, ERRO

 r  = dsqrt((Tx-x)**2+(Ty-y)**2)
 if( r <= 1.d-10)then
  r= 1.d-1
 end if

 IS = ISOURCE
 a  = SDIMENS
 h0 = Tz
 cz = z
 w  = (pi+pi)*valf
 Ephip = (0.d0,0.d0)

call regiao_coord_z( meio, z, kmada )

if( meio == 1 )then

	!===========================================================================================
          call algoritmoepsilon_zerosreais_pll( nprocss, kernelephi_1, zrs_J1, tolrel, tolabs, xGssLg, wGssLg, Ephip, ERRO )
	!===========================================================================================

else if( meio == 2 )then  

	!===========================================================================================    
          call algoritmoepsilon_zerosreais_pll( nprocss, kernelephi_2, zrs_J1, tolrel, tolabs, xGssLg, wGssLg, Ephip, ERRO )
	!===========================================================================================	 				

else if( meio == 3 )then  

	!===========================================================================================
          call algoritmoepsilon_zerosreais_pll( nprocss, kernelephi_3, zrs_J1, tolrel, tolabs, xGssLg, wGssLg, Ephip, ERRO )
	!===========================================================================================

else if( meio == 4 )then  

	!===========================================================================================
          call algoritmoepsilon_zerosreais_pll( nprocss, kernelephi_4, zrs_J1, tolrel, tolabs, xGssLg, wGssLg, Ephip, ERRO )
	!===========================================================================================			

end if


Ex_p = -(y-Ty)*Ephip/r
Eyp   = (x-Tx)*Ephip/r

IF( ERRO == 1 )THEN
      print*,'=================================================================='
      print*,'The epsilon algorithm did not achieve the desired convergence with'
      print*,'number of intervals and pre-set quadrature points'
      print*,'=================================================================='
      STOP	
END IF

END SUBROUTINE BCH1DQWE_pll_campo_E

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================
 
SUBROUTINE BCH1DQWE_pll_campo_H( nprocss, valf, Tx, Ty, Tz, x, y, z, Hrp, Hxp, Hyp, Hzp )
IMPLICIT NONE
COMPLEX(dp), INTENT(out)   :: Hrp, Hxp, Hyp, Hzp
REAL(dp), INTENT(in)	   :: valf, Tx, Ty, Tz, x, y, z
INTEGER , INTENT(in)	   :: nprocss	

INTEGER	 :: meio, ERRO

 r  = dsqrt((Tx-x)**2+(Ty-y)**2)
 if( r <= 1.d-10)then
  r= 1.d-1
 end if

 IS = ISOURCE
 a  = SDIMENS
 h0 = Tz
 cz = z
 w  = (pi+pi)*valf

call regiao_coord_z( meio, z, kmada )

if( meio == 1 )then

	!===========================================================================================
	  call algoritmoepsilon_zerosreais_pll( nprocss, kernelhr_1, zrs_J1, tolrel, tolabs, xGssLg, wGssLg, Hrp, ERRO )
	  call algoritmoepsilon_zerosreais_pll( nprocss, kernelhz_1, zrs_J0, tolrel, tolabs, xGssLg, wGssLg, Hzp, ERRO )
	!===========================================================================================

else if( meio == 2 )then  

	!===========================================================================================    
	  call algoritmoepsilon_zerosreais_pll( nprocss, kernelhr_2, zrs_J1, tolrel, tolabs, xGssLg, wGssLg, Hrp, ERRO )
	  call algoritmoepsilon_zerosreais_pll( nprocss, kernelhz_2, zrs_J0, tolrel, tolabs, xGssLg, wGssLg, Hzp, ERRO )
	!==========================================================================================	 				

else if( meio == 3 )then  


	!===========================================================================================
	  call algoritmoepsilon_zerosreais_pll( nprocss, kernelhr_3, zrs_J1, tolrel, tolabs, xGssLg, wGssLg, Hrp, ERRO )
	  call algoritmoepsilon_zerosreais_pll( nprocss, kernelhz_3, zrs_J0, tolrel, tolabs, xGssLg, wGssLg, Hzp, ERRO )
	!===========================================================================================

else if( meio == 4 )then  

	!===========================================================================================
	  call algoritmoepsilon_zerosreais_pll( nprocss, kernelhr_4, zrs_J1, tolrel, tolabs, xGssLg, wGssLg, Hrp, ERRO )
	  call algoritmoepsilon_zerosreais_pll( nprocss, kernelhz_4, zrs_J0, tolrel, tolabs, xGssLg, wGssLg, Hzp, ERRO )
	!===========================================================================================			!

end if


Hxp = (x-Tx)*Hrp/r
Hyp = (y-Ty)*Hrp/r
Hzp = -Hzp/ci/w/mi0

IF( ERRO == 1 )THEN
      print*,'=================================================================='
      print*,'The epsilon algorithm did not achieve the desired convergence with'
      print*,'number of intervals and pre-set quadrature points'
      print*,'=================================================================='
      STOP	
END IF
			
END SUBROUTINE BCH1DQWE_pll_campo_H

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE regiao_coord_z( meio, z, kmada )
IMPLICIT NONE
INTEGER, INTENT(out)	:: meio, kmada
REAL(dp), INTENT(in)	:: z

INTEGER	 :: cont, i

!======================================================================================
if( NC == 1 )then
	if( z < h0 )then
		meio = 1  ;  kmada = 0 
	else if( z > h0 .and. z < 0.d0 )then 
		meio = 2  ;  kmada = 0 
	else if( z >= 0.d0 )then
		meio = 4  ;  kmada = 1 ;  zinf = 0.d0
	end if
else
	if( z < h0 .and. h0 < 0.d0 )then
		 meio = 1  ;  kmada = 0
	else if( z >= h0 .and. z < 0.d0 )then 
		meio = 2  ;  kmada = 0
	else if( z >= 0.d0 )then
	
		cont = 0 ;	zinf = 0.d0 ;	   zsup = hj(1)
		do i=1, NC      

			meio = 3 ; kmada = i
			if( zinf == sum(hj(1:NC-1)) )then
				meio = 4
			end if
			if( zinf <= z .and. z < zsup )exit
		        if( i < NC )then
				zinf = zinf + hj(i) 
				zsup = zsup + hj(i+1)            	
			end if
		end do
	end if
end if

END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

COMPLEX(dp) FUNCTION kernelephi_1( g )
IMPLICIT NONE
REAL(dp), INTENT(in):: g 
 
   kernelephi_1 = Tej(g , kmada)*( cdexp(uj(g, kmada)*(cz-h0)) + Rtej(g, kmada)* &
							cdexp(uj(g,kmada)*(cz+h0)) )*g*dbesj1(g)

END FUNCTION kernelephi_1

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

COMPLEX(dp) FUNCTION kernelephi_2( g )
IMPLICIT NONE
REAL(dp), INTENT(in):: g
 
   kernelephi_2 = Tej(g,kmada)*( cdexp(-uj(g, kmada)*(cz-h0)) + Rtej(g,kmada)* &
						cdexp(uj(g,kmada)*(cz+h0)) )*g*dbesj1(g)
   
END FUNCTION kernelephi_2

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

COMPLEX(dp) FUNCTION kernelephi_3( g )
IMPLICIT NONE
REAL(dp), INTENT(in):: g 

   kernelephi_3 = Tej(g,kmada)*( cdexp(-uj(g,kmada)*(cz-zinf)) + Rtej(g,kmada)* &
					cdexp(uj(g,kmada)*(cz-zsup-hj(kmada))) )*g*dbesj1(g)
 
END FUNCTION kernelephi_3

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

COMPLEX(dp) FUNCTION kernelephi_4( g )
IMPLICIT NONE
REAL(dp), INTENT(in):: g
 
   kernelephi_4 = Tej(g,kmada)*cdexp(-uj(g,kmada)*(cz-zinf))*g*dbesj1(g)

END FUNCTION kernelephi_4

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

COMPLEX(dp) FUNCTION kernelhr_1( g )
IMPLICIT NONE
REAL(dp), INTENT(in):: g 
 
   kernelhr_1 = Yintj(g,kmada)*Tej(g,kmada)*( cdexp(uj(g,kmada)*(cz-h0)) + Rtej(g,kmada)* &
					cdexp(uj(g,kmada)*(cz+h0)) )*g*dbesj1(g)

END FUNCTION kernelhr_1

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

COMPLEX(dp) FUNCTION kernelhr_2( g )
IMPLICIT NONE
REAL(dp), INTENT(in):: g 

   kernelhr_2 = -Yintj(g,kmada)*Tej(g,kmada)*( cdexp(-uj(g,kmada)*(cz-h0)) - Rtej(g,kmada)* &
					cdexp(uj(g,kmada)*(cz+h0)) )*g*dbesj1(g)

END FUNCTION kernelhr_2

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

COMPLEX(dp) FUNCTION kernelhr_3( g )
IMPLICIT NONE
REAL(dp), INTENT(in):: g 

   kernelhr_3 = -Yintj(g,kmada)*Tej(g,kmada)*( cdexp(-uj(g,kmada)*(cz-zinf)) - Rtej(g,kmada)* & 
	 				cdexp(uj(g,kmada)*(cz-zsup-hj(kmada))) )*g*dbesj1(g)



END FUNCTION kernelhr_3

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

COMPLEX(dp) FUNCTION kernelhr_4( g )
IMPLICIT NONE
REAL(dp), INTENT(in):: g 

   kernelhr_4 = -Yintj(g,kmada)*Tej(g,kmada)*cdexp(-uj(g,kmada)*(cz-zinf))*g*dbesj1(g)


END FUNCTION kernelhr_4

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

COMPLEX(dp) FUNCTION kernelhz_1( g )
IMPLICIT NONE
REAL(dp), INTENT(in):: g 
 
   kernelhz_1 = Tej_Hz(g , kmada)*( cdexp(uj(g, kmada)*(cz-h0)) + Rtej(g, kmada)* &
							cdexp(uj(g,kmada)*(cz+h0)) )*g**2

END FUNCTION kernelhz_1

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

COMPLEX(dp) FUNCTION kernelhz_2( g )
IMPLICIT NONE
REAL(dp), INTENT(in):: g
 
   kernelhz_2 = Tej_Hz(g,kmada)*( cdexp(-uj(g, kmada)*(cz-h0)) + Rtej(g,kmada)* &
						cdexp(uj(g,kmada)*(cz+h0)) )*g**2

END FUNCTION kernelhz_2

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

COMPLEX(dp) FUNCTION kernelhz_3( g )
IMPLICIT NONE
REAL(dp), INTENT(in):: g 
 
   kernelhz_3 = Tej_Hz(g,kmada)*( cdexp(-uj(g,kmada)*(cz-zinf)) + Rtej(g,kmada)* &
					cdexp(uj(g,kmada)*(cz-zsup-hj(kmada))) )*g**2
END FUNCTION kernelhz_3

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

COMPLEX(dp) FUNCTION kernelhz_4( g )
IMPLICIT NONE
REAL(dp), INTENT(in):: g 
 
   kernelhz_4 = Tej_Hz(g,kmada)*cdexp(-uj(g,kmada)*(cz-zinf))*g**2

END FUNCTION kernelhz_4

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

COMPLEX(dp) FUNCTION Uj( g, kmdj )

IMPLICIT NONE
INTEGER, INTENT(in):: kmdj
REAL(dp), INTENT(in):: g

COMPLEX(dp):: zetaj, netaj
REAL(8):: gra

if( r >= a )then
        gra = g/r
else
	gra = g/a
end if

if( kmdj == 0 )then
	zetaj = ci*w*mi0
	netaj = ci*w*epson
	uj = cdsqrt( (gra)**2.d0 + netaj*zetaj )         
else
	zetaj = dcmplx( 0.d0, w*mi0 )
	netaj = dcmplx( 1.d0/roj(kmdj), 0.d0 )
	uj = cdsqrt( (gra)**2 + netaj*zetaj )
end if

END FUNCTION

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

COMPLEX(dp) FUNCTION Yintj( g, kmdj )
IMPLICIT NONE
REAL(dp), INTENT(in):: g
INTEGER, INTENT(in):: kmdj
COMPLEX(dp):: zetaj

	zetaj = dcmplx( 0.d0, w*mi0 ) 
	Yintj = uj(g, kmdj)/zetaj

END FUNCTION
!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

COMPLEX(dp) FUNCTION Yapj( g, kmdj )
IMPLICIT NONE
INTEGER, INTENT(in):: kmdj
REAL(dp), INTENT(in):: g
COMPLEX(dp), ALLOCATABLE, DIMENSION(:):: adm_apj
INTEGER:: j

allocate( adm_apj(NC) )

	adm_apj(NC) = Yintj( g, NC )

if( kmdj < NC ) then

	do j=NC-1, kmdj, -1
		adm_apj(j) = Yintj(g,j)*( adm_apj(j+1) + Yintj(g,j)*cdtanh(uj(g,j)*hj(j)) )/ & 
			                  ( Yintj(g,j) + adm_apj(j+1)*cdtanh(uj(g,j)*hj(j)) )
	end do
	Yapj = adm_apj(kmdj) 
else
	Yapj = adm_apj(NC) 
end if

deallocate(adm_apj)

END FUNCTION

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

COMPLEX(dp) FUNCTION Rtej( g, kmdj )
IMPLICIT NONE
INTEGER, INTENT(in):: kmdj
REAL(dp), INTENT(in):: g

if( NC == 1 )then       
	Rtej = (Yintj(g,0)-Yintj(g,1))/(Yintj(g,0)+Yintj(g,1))
else
	Rtej = (Yintj(g,kmdj)-Yapj(g,kmdj+1)) / (Yintj(g,kmdj)+Yapj(g,kmdj+1))
end if

END FUNCTION

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

COMPLEX(dp) FUNCTION Tej( g, kmdj )
IMPLICIT NONE
INTEGER, INTENT(in):: kmdj
REAL(dp), INTENT(in):: g
COMPLEX(dp), ALLOCATABLE, DIMENSION(:):: Tej_temp
COMPLEX(dp):: zeta0

INTEGER:: j
REAL(dp)::ra

allocate( Tej_temp(0:NC) )

zeta0=(0.d0, 1.d0)*w*mi0
if( r >= a )then
        ra = a/r
	Tej_temp(0) = -a*IS*zeta0*dbesj1(g*ra)/2.d0/uj(g,0)/r**2
else
	ra = r/a
	Tej_temp(0) = -a*IS*zeta0*dbesj1(g*ra)/2.d0/uj(g,0)/a**2
end if
if( kmdj == 1 .and. NC == 1 ) then

	Tej_temp(kmdj) = 2.d0*Yintj(g, kmdj-1)*Tej_temp(kmdj-1) / (Yintj(g,kmdj-1)+Yintj(g,kmdj))

else if( NC > 1 ) then

    do j=1, kmdj
	if( j == 1 ) then
		Tej_temp(1) = Tej_temp(0)*( 1.d0+Rtej(g,0) ) / ( 1.d0+Rtej(g,1)*cdexp(-2.d0*uj(g,1)*hj(1)) )
        else if( j > 1 .and. j < NC )then
 	        Tej_temp(j) = Tej_temp(j-1)*( 1d0+Rtej(g,j-1) )*cdexp(-uj(g,j-1)*hj(j-1)) / &
					    ( 1.d0+Rtej(g,j)*cdexp(-2.d0*uj(g,j)*hj(j)) )
	else if( j == NC )then
		Tej_temp(j) = Tej_temp(j-1)*( 1.d0+Rtej(g,j-1) )*cdexp(-uj(g,j-1)*hj(j-1))
	end if
   end do

end if

	Tej = Tej_temp(kmdj)

deallocate(Tej_temp)

END FUNCTION

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

COMPLEX(dp) FUNCTION Tej_Hz( g, kmdj )
IMPLICIT NONE
INTEGER, INTENT(in):: kmdj
REAL(dp), INTENT(in):: g
COMPLEX(dp), ALLOCATABLE, DIMENSION(:):: Tej_temp
COMPLEX(dp):: zeta0

INTEGER:: j
REAL(dp)::ra

allocate( Tej_temp(0:NC) )

zeta0=(0.d0, 1.d0)*w*mi0
if( r >= a )then
        ra = a/r
	Tej_temp(0) = -a*IS*zeta0*dbesj1(g*ra)*dbesj0(g)/2.d0/uj(g,0)/r**3
else
	ra = r/a
	Tej_temp(0) = -a*IS*zeta0*dbesj1(g)*dbesj0(g*ra)/2.d0/uj(g,0)/a**3
end if
if( kmdj == 1 .and. NC == 1 ) then

	Tej_temp(kmdj) = 2.d0*Yintj(g, kmdj-1)*Tej_temp(kmdj-1) / (Yintj(g,kmdj-1)+Yintj(g,kmdj))

else if( NC > 1 ) then

    do j=1, kmdj
	if( j == 1 ) then
		Tej_temp(1) = Tej_temp(0)*( 1.d0+Rtej(g,0) ) / ( 1.d0+Rtej(g,1)*cdexp(-2.d0*uj(g,1)*hj(1)) )
        else if( j > 1 .and. j < NC )then
 	        Tej_temp(j) = Tej_temp(j-1)*( 1d0+Rtej(g,j-1) )*cdexp(-uj(g,j-1)*hj(j-1)) / &
					    ( 1.d0+Rtej(g,j)*cdexp(-2.d0*uj(g,j)*hj(j)) )
	else if( j == NC )then
		Tej_temp(j) = Tej_temp(j-1)*( 1.d0+Rtej(g,j-1) )*cdexp(-uj(g,j-1)*hj(j-1))
	end if
   end do

end if

	Tej_Hz = Tej_temp(kmdj)

deallocate(Tej_temp)

END FUNCTION

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

COMPLEX(dp) FUNCTION cdtanh(z)
IMPLICIT NONE
COMPLEX(dp):: z

    cdtanh = (1.d0-cdexp(-(z+z))) / (1.d0+cdexp(-(z+z)))

END FUNCTION

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

subroutine algoritmoepsilon_zerosreais( f, zeros, tol_rel, tol_abs, x, w, s, ERRO )
!================================================================================================================
! 
!   INPUT:
! 
!   f       : integrating function that must be written with integration variable the argument of the Bessel function.
! 
!   zeros   : zeros to be considered. The number of zeros also determines the maximum number of intervals to be 
!             considered in extrapolation. If there is no convergence, increasing the number of zeros may be the 
!             solution. In general I did not need more than 1000 intervals. 
! 
!   tol_rel : relative tolerance between the partial sums. In general do not adopt greater than 1.d-6; In mCSEM 
!             good results are obtained from 1.d-9, or lower.
! 
!   tol_abs : absolute tolerance between the partial sums. I've always used 1.d-30 without problems;
! 
!   x       : abscissas of Gauss-Legendre;
!   w       :Gauss-Legendre weights.
! 
!   output:
! 
!   s       : integration result.
! 
!================================================================================================================
    implicit none
	complex(dp), external :: f
    real(dp), intent(in) :: zeros(:), x(:), w(:)
    real(dp), intent(in) :: tol_rel, tol_abs
    complex(dp), intent(out) :: s
    INTEGER, INTENT(out) :: ERRO 

	integer(4), parameter :: nini = 1
	real(dp), parameter :: small = 1.d-300
	real(dp), parameter :: big = 1.d+300
	real(dp), parameter :: x0 = 1.d-20
    integer(4) :: ntermos, i, j, n
    integer(4) :: nmaxintrvl, sit
    complex(dp), dimension(:), allocatable :: E, extrap
	real(dp), dimension(:), allocatable::  er_abs, er_rel, zero
	
	complex(dp) :: soma, res, sn_1, aux2, aux1, dif
	real(dp) :: x_a, x_b

	ERRO = 0
	nmaxintrvl = size(zeros)
    ntermos = nmaxintrvl - nini - 1
    allocate( E( ntermos ), extrap( ntermos ), er_abs( ntermos ), er_rel( ntermos ) )
	E = 0.d0
	extrap = 0.d0
	er_abs = 0.d0

    soma = ( 0.d0, 0.d0 )
	res = (0.d0,0.d0)
	x_a = x0
	x_b = zeros(1)
    do i = 1, nini
        call quadGauLeg( f, x_a, x_b, x, w, res )
        soma = soma + res
		x_a = zeros( i )
		x_b = zeros( i + 1 )
    end do

    s = 0.d0
    sn_1 = 0.d0
    do i = nini + 1, ntermos

        sn_1 = sn_1 + s
        call quadGauLeg( f, x_a, x_b, x, w, res )        
		x_a = zeros( i )
		x_b = zeros( i + 1 )
        
        n = i - nini
        E ( n + 1 ) = E ( n ) + res
        
        aux2 = 0.d0
        do j = n + 1, 2, -1
            aux1 = aux2
            aux2 = E ( j - 1 )
            dif = E ( j ) - aux2
            if ( cdabs( dif ) < small ) then
                E ( j - 1 ) = big
            else
                E ( j - 1 ) = aux1 + 1.d0 / dif
            end if
        end do
        
        if ( mod( n, 2 )  == 0 ) then
            extrap(n) = E ( 1 ) + soma
        else
            extrap(n) = E ( 2 ) + soma
        end if

        if ( n > 1) then
            er_abs(n) = cdabs( extrap(n) - extrap(n-1) )
            er_rel(n) = er_abs(n) / cdabs( extrap(n) )
            	s = extrap(n)
            if ( er_rel(n) < ( tol_rel + tol_abs / cdabs( extrap(n) ) ) ) then
            	s = extrap(n)
                exit
            else if ( i == ntermos ) then
		ERRO = 1
		print*,'Erro relativo percentual',er_rel(n)
            end if
        end if

    end do

end subroutine algoritmoepsilon_zerosreais

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

subroutine algoritmoepsilon_zerosreais_pll( nprocss, f, zeros, tol_rel, tol_abs, x, w, s, ERRO )
!================================================================================================================
! 
!   INPUT:
! 
!   f       : integrating function that must be written with integration variable the argument of the Bessel function.
! 
!   zeros   : zeros to be considered. The number of zeros also determines the maximum number of intervals to be 
!             considered in extrapolation. If there is no convergence, increasing the number of zeros may be the 
!             solution. In general I did not need more than 1000 intervals. 
! 
!   tol_rel : relative tolerance between the partial sums. In general do not adopt greater than 1.d-6; In mCSEM 
!             good results are obtained from 1.d-9, or lower.
! 
!   tol_abs : absolute tolerance between the partial sums. I've always used 1.d-30 without problems;
! 
!   x       : abscissas of Gauss-Legendre;
!   w       :Gauss-Legendre weights.
! 
!   output:
! 
!   s       : integration result.
! 
!================================================================================================================
    implicit none
    complex(dp), external    :: f
    real(dp), intent(in)     :: zeros(:), x(:), w(:)
    real(dp), intent(in)     :: tol_rel, tol_abs
    complex(dp), intent(out) :: s
    INTEGER, INTENT(out)     :: ERRO
    INTEGER, INTENT(in)      :: nprocss  

	integer(4), parameter :: nini = 1
	real(dp), parameter :: small = 1.d-300
	real(dp), parameter :: big = 1.d+300
	real(dp), parameter :: x0 = 1.d-20
    integer(4) :: ntermos, i, j, n
    integer(4) :: nmaxintrvl, sit
    complex(dp), dimension(:), allocatable :: E, extrap
	real(dp), dimension(:), allocatable::  er_abs, er_rel, zero
	
	complex(dp) :: soma, res, sn_1, aux2, aux1, dif
	real(dp) :: x_a, x_b

	ERRO = 0
	nmaxintrvl = size(zeros)
    ntermos = nmaxintrvl - nini - 1
    allocate( E( ntermos ), extrap( ntermos ), er_abs( ntermos ), er_rel( ntermos ) )
	E = 0.d0
	extrap = 0.d0
	er_abs = 0.d0

    soma = ( 0.d0, 0.d0 )
	res = (0.d0,0.d0)
	x_a = x0
	x_b = zeros(1)
    do i = 1, nini
        call quadGauLeg_pll( nprocss,f , x_a, x_b, x, w, res )
        soma = soma + res
		x_a = zeros( i )
		x_b = zeros( i + 1 )
    end do

    s = 0.d0
    sn_1 = 0.d0
    do i = nini + 1, ntermos

        sn_1 = sn_1 + s
        call quadGauLeg_pll( nprocss,f , x_a, x_b, x, w, res )        
		x_a = zeros( i )
		x_b = zeros( i + 1 )
        
        n = i - nini
        E ( n + 1 ) = E ( n ) + res
        
        aux2 = 0.d0
        do j = n + 1, 2, -1
            aux1 = aux2
            aux2 = E ( j - 1 )
            dif = E ( j ) - aux2
            if ( cdabs( dif ) < small ) then
                E ( j - 1 ) = big
            else
                E ( j - 1 ) = aux1 + 1.d0 / dif
            end if
        end do
        
        if ( mod( n, 2 )  == 0 ) then
            extrap(n) = E ( 1 ) + soma
        else
            extrap(n) = E ( 2 ) + soma
        end if

        if ( n > 1) then
            er_abs(n) = cdabs( extrap(n) - extrap(n-1) )
            er_rel(n) = er_abs(n) / cdabs( extrap(n) )
            	s = extrap(n)
            if ( er_rel(n) < ( tol_rel + tol_abs / cdabs( extrap(n) ) ) ) then
            	s = extrap(n)
                exit
            else if ( i == ntermos ) then
		ERRO = 1
		print*,'Erro relativo percentual',er_rel(n)
            end if
        end if

    end do

end subroutine algoritmoepsilon_zerosreais_pll

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

subroutine quadGauLeg( func, xa, xb, x, wt, s )
    implicit none
    real(dp), intent(in) :: xa, xb, x(:), wt(:)
    complex(dp), intent(out) :: s
    complex(dp), external :: func
	
    integer(4) :: i, n
    complex(dp) :: sx
    real(dp) :: hx, mx
    
    n = size(x)
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=132
!	Variable change for the interval [a, b] to become [-1,1].
!	We use the functions x = mx * t + hx ([-1,1] ---> [a, b]) that are given by:
	mx = 0.5d0 * ( xb - xa )	
	hx = 0.5d0 * ( xa + xb )
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=132
!	With this change of variable the integral, in 1D, becomes:
!	     b		  1
!           /		/	   			  
!	   |	       | b-a   [ b-a      b+a  ]
!	   | f(z) dz = | --- f| ---  t +  --- | dt	,
!	   |	       |  2   [  2	 2   ]
!	  /	       /		   
!	 a	       -1								
! 
!       then the integrand becomes:
! 
!	mx*f(mx*t+hx)
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=132
	sx = (0.d0,0.d0)
	do i = 1, n
        	sx = sx + func( mx * x( i ) + hx ) * wt( i )
        end do
	s = mx * sx
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=132
	end subroutine quadGauLeg

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

subroutine quadGauLeg_pll( nprocss, func, xa, xb, x, wt, s )
    implicit none
    real(dp)   , intent(in)  :: xa, xb, x(:), wt(:)
    complex(dp), intent(out) :: s
    complex(dp), external    :: func
    integer    ,intent(in)   :: nprocss   	

    integer(4) :: i, n
    complex(dp) :: sx
    real(dp) :: hx, mx
    
    n = size(x)
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=132
!	Variable change for the interval [a, b] to become [-1,1].
!	We use the functions x = mx * t + hx ([-1,1] ---> [a, b]) that are given by:
	mx = 0.5d0 * ( xb - xa )	
	hx = 0.5d0 * ( xa + xb )
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=132
!	With this change of variable the integral, in 1D, becomes:
!	     b		  1
!           /		/	   			  
!	   |	       | b-a   [ b-a      b+a  ]
!	   | f(z) dz = | --- f| ---  t +  --- | dt	,
!	   |	       |  2   [  2	 2   ]
!	  /	       /		   
!	 a	       -1								, então nosso integrando passa a ser:
!       then the integrand becomes:
! 
!	mx*f(mx*t+hx)
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=132
	sx = (0.d0,0.d0)

 	call OMP_SET_NUM_THREADS( nprocss )
        !$omp parallel  PRIVATE(i) REDUCTION( +:sx )
        !$omp DO  
		do i = 1, n
			!$OMP ATOMIC UPDATE
	        	sx = sx + func( mx * x( i ) + hx ) * wt( i )
	        end do
        !$omp end DO
        !$omp end parallel
	s = mx * sx
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=132
	end subroutine quadGauLeg_pll

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE EyHy_TETM( z, freq, md, EyHy, HxEx )
IMPLICIT NONE
INTEGER, INTENT(in):: md
REAL(dp), INTENT(in):: freq, z
COMPLEX(dp), INTENT(out):: EyHy, HxEx

COMPLEX(dp), ALLOCATABLE, DIMENSION(:):: uj, TEHj, YZapj, YZintj, Rtemj
REAL(dp):: interv_z, w
COMPLEX(dp):: k
INTEGER:: cn, meio, N

w = (pi+pi)*freq 

	if( NC > 1 )then
		CALL ident_kmada( z, hj, meio, NC, cn )
		CAll val_Uj_ZYint( freq, uj, YZintj, roj, mirj, NC, md )
		CAll ZY_int_apar( uj, YZintj, YZapj, hj, NC )
		CAll coef_reflexTETM( NC, YZintj, YZapj, Rtemj )
		CAll coef_transEjHj( uj, Rtemj, TEHj, hj, NC )
	end if

	if( NC == 1 )then

		k = cdsqrt(-ci*w*mi0*(1.d0/roj(1) + ci*w*epson) ) 
		EyHy = cdexp(-ci*k*z)
		if( md==1 )then
			HxEx =  (-k/w/mi0)*EyHy
		else if( md==2 )then
			HxEx =  ci*k*roj(1)*EyHy

		end if

	else if( NC > 1 .and. cn < NC )then

		interv_z = sum(hj(1:cn))

		EyHy = TEHj(cn)*(cdexp(-uj(cn)*(z-interv_z)) + Rtemj(cn)*cdexp(uj(cn)*(z-interv_z)))
		if( md == 1 )then
			YZintj(cn) = -YZintj(cn)
		end if
		HxEx = YZintj(cn)*TEHj(cn)*(cdexp(-uj(cn)*(z-interv_z)) - Rtemj(cn)*cdexp(uj(cn)*(z-interv_z)) )

	else if( cn == NC )then

		interv_z = sum(hj(1:cn-1))
		EyHy = TEHj(cn)*cdexp(-uj(cn)*(z-interv_z))
		if( md == 1 )then
			YZintj(cn) = -YZintj(cn)
		end if
		HxEx = YZintj(cn)*TEHj(cn)*EyHy

	end if

END SUBROUTINE

!======================================================================================
!88888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!======================================================================================

SUBROUTINE ident_kmada( z, hj, meio, NC, kmada )
IMPLICIT NONE
REAL(dp), INTENT(in):: z
REAL(dp), DIMENSION(:), INTENT(in):: hj
INTEGER, INTENT(inout):: NC, meio, kmada

REAL(dp):: zinf, zsup, h0=0.d0 
INTEGER:: cont, i

!======================================================================================
if( z <= h0 .and. h0 <= 0 )then
	 meio = 1  ;  kmada = 0
else if( z > h0 .and. z <= 0 )then 
	meio = 2  ;  kmada = 0
else if( z >= 0 )then

	cont = 0 ;	zinf = 0.d0 ;	   zsup = hj(1)
	do i=1, NC            

		meio = 3 ; kmada = i
		if( zinf < z .and. z <= zsup )exit
	        if( i < NC-1 )then
			zinf = zinf + hj(i) 
			zsup = zsup + hj(i+1)            	
		end if
	end do
end if
!======================================================================================

END SUBROUTINE

!================================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!================================================================================================================

SUBROUTINE val_Uj_ZYint( freq, uj, YZintj, roj, mirj, N , md)

IMPLICIT NONE
INTEGER:: N , md
REAL(dp), INTENT(in):: freq
REAL(dp), DIMENSION(:), INTENT(in):: roj, mirj
COMPLEX(dp), ALLOCATABLE, DIMENSION(:), INTENT(out):: uj, YZintj 

COMPLEX(dp) zeta, neta0
real(dp):: w

allocate( uj(0:N), YZintj(0:N) )

w = 8.d0*datan(1.d0)*freq

neta0 = (0.d0,1.d0)*w*epson    !
zeta  = (0.d0,1.d0)*w*mi0

uj(0)   = cdsqrt(zeta*neta0)
uj(1:N) = (ci+1)*dsqrt(w*mi0*mirj(1:N)/roj(1:N)/2.d0)

if( md == 1 )then
	YZintj(0)   = uj(0)/zeta ! TE
	YZintj(1:N) = uj(1:N)/w/mi0/mirj(1:N)/ci ! TE
else if( md == 2 )then
	YZintj(0)   = uj(0)/neta0
	YZintj(1:N) = uj(1:N)*roj(1:N)	! TM
end if

END SUBROUTINE

!=================================================================================================================
!88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=================================================================================================================

SUBROUTINE ZY_int_apar( uj, YZintj, YZapj, hj, N )

IMPLICIT NONE
INTEGER, INTENT(in):: N
REAL(dp), DIMENSION(:), INTENT(in):: hj
COMPLEX(dp), 		INTENT(in) :: uj(0:N), YZintj(0:N)
COMPLEX(dp), ALLOCATABLE, DIMENSION(:), INTENT(out):: YZapj

INTEGER:: j

allocate( YZapj(N) )
YZapj(N) = YZintj(N)

if (N > 1) then
	do j=N-1, 1, -1
	        YZapj(j) = YZintj(j)*( YZapj(j+1)+ YZintj(j)*arctanh(uj(j)*hj(j)) ) / &
			             ( YZintj(j) + YZapj(j+1)*arctanh(uj(j)*hj(j)) )  
	end do
end if
END SUBROUTINE

!=================================================================================================================
!88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=================================================================================================================

SUBROUTINE coef_reflexTETM( N, YZintj, YZapj, Rj )

IMPLICIT NONE
INTEGER    ,		 		INTENT(in)  :: N
COMPLEX(dp), 				INTENT(in)  :: YZapj(N), YZintj(0:N)
COMPLEX(dp), ALLOCATABLE, DIMENSION(:), INTENT(out) :: Rj        
INTEGER:: j

if (N > 1) then
allocate( Rj(0:N-1) )

	do j=0,N-1 
		Rj(j) = ( YZintj(j)-YZapj(j+1) ) / (YZintj(j) + YZapj(j+1))
	end do
end if
END SUBROUTINE

!=================================================================================================================
!88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=================================================================================================================

SUBROUTINE coef_transEjHj( uj, Rj, Tj, hj, N )
IMPLICIT NONE
INTEGER,		INTENT(in) :: N
REAL(dp), DIMENSION(:), INTENT(in) :: hj
COMPLEX(dp), 		INTENT(in) :: uj(0:N), Rj(0:N)

COMPLEX(dp), ALLOCATABLE, DIMENSION(:), INTENT(out):: Tj        
INTEGER:: j

allocate( Tj(0:N) )

if(N == 1)then
	Tj(N) = 1.0d0
else
	Tj(0) = 1.d0/(1.0d0 + Rj(0))
	do j=1, N-1
		Tj(j) = Tj(j-1)*(1+Rj(j-1))*cdexp(-uj(j)*hj(j)) /(1.d0+Rj(j)*cdexp(-2.d0*uj(j)*hj(j)))
	end do
	Tj(N) = Tj(N-1)*(1+Rj(N-1))
end if

END SUBROUTINE

!=================================================================================================================
!88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=================================================================================================================

COMPLEX FUNCTION arctanh(z)
IMPLICIT NONE
COMPLEX(dp):: z

	arctanh = (1.d0-cdexp(-2.d0*z)) / (1.d0+cdexp(-2.d0*z))

END FUNCTION

!=================================================================================================================
!88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=================================================================================================================

END MODULE mdlgm2_1D
