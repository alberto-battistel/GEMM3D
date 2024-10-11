module  mdlgm1_1D
!======================================================================================================
!
! Autor: The subroutines included in this module were implemented and generously hosted by Professor 
!
!        Dr. Valdelírio da Silva e Silva universidade federal do Pará(UFPA) in Brazil.
!
!======================================================================================================
contains
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
	subroutine dehx_xyz_loops( idtfcd_cJ0, nJ0, idtfcd_cJ1, nJ1, Tx, Ty, Iw, dsx, h0, n, esp, condut, &
						 neta, zeta, cx, cy, z, Ex_p, Ey_p, Ez_p, Hx_p, Hy_p, Hz_p)
	implicit none
	integer(4), intent(in)  :: n, idtfcd_cJ0, nJ0, idtfcd_cJ1, nJ1
	real(8)   , intent(in)  :: Tx, Ty, Iw, dsx, h0, esp(:), condut(1:n), cx, cy, z !,neta
	complex*16, intent(in)  :: zeta,neta
	complex*16, intent(out) :: Ex_p, Ey_p, Ez_p, Hx_p, Hy_p, Hz_p

	real(8),parameter::pi=3.141592653589793238462643383279502884197d0
	integer(4)::i,j,k,camad,camadT,filtro,ident_fJ0,ident_fJ1
	real(8)::x,y,r,kr,k_r
	real(8),dimension(:),allocatable::h,krJ0,krJ1,w_J0,w_J1,prof

!	Para uso de loops:
	complex*16,dimension(:),allocatable::wvnb2,AMdwJ0,AMdwJ1,AMupJ0,AMupJ1,FEdwJ0,FEdwJ1,FEupJ0,FEupJ1
	complex*16,dimension(:,:),allocatable::uJ0,uJ1,AdmIntJ0,AdmIntJ1,ImpIntJ0,ImpIntJ1,uhJ0,uhJ1,tghJ0,tghJ1
	complex*16,dimension(:,:),allocatable::AdmApdwJ0,AdmApdwJ1,ImpApdwJ0,ImpApdwJ1,AdmApupJ0,AdmApupJ1,ImpApupJ0,ImpApupJ1
	complex*16,dimension(:,:),allocatable::RTEdwJ0,RTEdwJ1,RTMdwJ0,RTMdwJ1,RTEupJ0,RTEupJ1,RTMupJ0,RTMupJ1
	complex*16,dimension(:,:),allocatable::TMdwJ0,TMdwJ1,TEdwJ0,TEdwJ1,TMupJ0,TMupJ1,TEupJ0,TEupJ1

	complex*16,dimension(:),allocatable::Ktmdz_J0,Ktmdz_J1,Ktm_J0,Ktm_J1,Kte_J0,Kte_J1,Ktedz_J0,Ktedz_J1
	complex*16,dimension(:),allocatable::kernelExJ0,kernelExJ1,kernelEyJ0,kernelEyJ1,kernelEzJ1
	complex*16,dimension(:),allocatable::kernelHxJ0,kernelHxJ1,kernelHyJ0,kernelHyJ1,kernelHzJ1

	complex*16::kerEx_J1,kerEx_J0,kerEy_J1,kerEy_J0,kerEz_J1
	complex*16::kerHx_J1,kerHx_J0,kerHy_J1,kerHy_J0,kerHz_J1

	if (cx==Tx) then
	x=1.d0!*Tx/dabs(Tx)
	else
	x=cx-Tx
	end if
	y=cy-Ty
	r = dsqrt(x**2 + y**2)

	allocate(h(0:n),prof(-1:n))
	if (size(esp)==n) then
		h(0)=0.d0
		h(1:n)=esp
	else
		h(0)=0.d0
		h(1:n-1)=esp
		h(n)=1.d300
	end if
!	criando um novo vetor de profundidades que se adeque à qualquer situação patológica
	prof(-1)=-1.d300
	prof(0)=0.d0
	if (n > 1) then
	 prof(1) = h(1)
	 if (n > 2) then
	  do k=2,n-1
	   prof(k) = prof(k-1) + h(k)
	  end do
	 end if
	end if
	prof(n)=1.d300

!para descobrir em que camada está a observação
	if (z < 0.d0) then
		camad=0
	else if (z >= prof(n-1)) then
		camad=n
	else
	do i=n-1,1,-1
	if (z >= prof(i-1)) then
		camad=i
		exit
	end if
	end do
	end if

!para descobrir em que camada está o transmissor
	if (h0 < 0.d0) then
		camadT = 0
	else if (h0 >= prof(n-1)) then
		camadT = n
	else
	do j=n-1,1,-1
	if (h0 >= prof(j-1)) then
		camadT = j
		exit
	end if
	end do
	end if

!!	write(*,*)'Entre com o criador dos filtros J0: Rijo(0), Frayzer(1), Guptasarma(2), Kong(3) ou Key(4)'
!!	read(*,*)idtfcd_cJ0
	filtro=0	!esta variável direciona o uso de filtros J0 e J1 em vez de seno e cosseno
!!	
!!	As variáveis idtfcd_cJ0 e idtfcd_cJ1 são do tipo inteiro 
!!      e estão associadas aos seguinte desenvovedor do filtro 
!!
!!      PARA idtfcd_cJ0 ou _cJ1 (filtros J0 e J1)::  

!!	idtfcd_cJ0 ou _cJ1 = 0 ---> Rijo 
!!	idtfcd_cJ0 ou _cJ1 = 1 ---> Frayzer
!!	idtfcd_cJ0 ou _cJ1 = 2 ---> Guptasarma
!!	idtfcd_cJ0 ou _cJ1 = 3 ---> Kong
!!	idtfcd_cJ0 ou _cJ1 = 4 ---> Key


!!	PARA o de Rijo       ---> nJ0 = 61 		;  nJ1 = 47 
!!	PARA o de Frayzer    ---> nJ0 = 37  		;  nJ1 = 27
!!	PARA o de Guptasarma ---> nJ0 = 61 ou 120       ;  nJ1 = 47 ou 140.
!!	PARA o de Kong       ---> nJ0 = 61 ou 241       ;  nJ1 = 61 ou 241
!!	PARA o de Key	     ---> nJ0 = 101, 201 ou 401 ;  nJ1 = 101, 201 ou 401

!	idtfcd_cJ0 = 4
	ident_fJ0 = 0
!	nJ0 = 401
!	idtfcd_cJ1 = 4
	ident_fJ1 = 1
!	nJ1 = 401

	allocate(KrJ0(nJ0),KrJ1(nJ1),w_J0(nJ0),w_J1(nJ1))

	call constfiltro(filtro,idtfcd_cJ0,ident_fJ0,nJ0,r,KrJ0,w_J0)
	call constfiltro(filtro,idtfcd_cJ1,ident_fJ1,nJ1,r,KrJ1,w_J1)

 	allocate(wvnb2(0:n),uJ0(nJ0,0:n),uJ1(nJ1,0:n),AdmIntJ0(nJ0,0:n),AdmIntJ1(nJ1,0:n),ImpIntJ0(nJ0,0:n),ImpIntJ1(nJ1,0:n))
 	allocate(uhJ0(nJ0,0:n),uhJ1(nJ1,0:n),tghJ0(nJ0,0:n),tghJ1(nJ1,0:n))
! 
 	do i=0,n
 		if (i==0) then
 		wvnb2(i)=-zeta*neta
 		uJ0(:,i)=cdsqrt(krJ0*krJ0-wvnb2(i))
 		uJ1(:,i)=cdsqrt(krJ1*krJ1-wvnb2(i))
 		AdmIntJ0(:,i)=uJ0(:,i)/zeta
 		AdmIntJ1(:,i)=uJ1(:,i)/zeta
 		ImpIntJ0(:,i)=uJ0(:,i)/neta
 		ImpIntJ1(:,i)=uJ1(:,i)/neta
 		uhJ0(:,i)=uJ0(:,i)*h(i)
 		uhJ1(:,i)=uJ1(:,i)*h(i)
 		tghJ0(:,i)=(1.d0-cdexp(-2.d0*uhJ0(:,i)))/(1.d0+cdexp(-2.d0*uhJ0(:,i)))
 		tghJ1(:,i)=(1.d0-cdexp(-2.d0*uhJ1(:,i)))/(1.d0+cdexp(-2.d0*uhJ1(:,i)))
 		else
 		wvnb2(i) = -zeta*condut(i)
 		uJ0(:,i)=cdsqrt(krJ0*krJ0-wvnb2(i))
 		uJ1(:,i)=cdsqrt(krJ1*krJ1-wvnb2(i))
 		AdmIntJ0(:,i)=uJ0(:,i)/zeta
 		AdmIntJ1(:,i)=uJ1(:,i)/zeta
 		ImpIntJ0(:,i)=uJ0(:,i)/condut(i)
 		ImpIntJ1(:,i)=uJ1(:,i)/condut(i)
 		uhJ0(:,i)=uJ0(:,i)*h(i)
 		uhJ1(:,i)=uJ1(:,i)*h(i)
 		tghJ0(:,i)=(1.d0-cdexp(-2.d0*uhJ0(:,i)))/(1.d0+cdexp(-2.d0*uhJ0(:,i)))
 		tghJ1(:,i)=(1.d0-cdexp(-2.d0*uhJ1(:,i)))/(1.d0+cdexp(-2.d0*uhJ1(:,i)))
 		end if
 	end do
 
 	allocate(AdmApdwJ0(nJ0,1:n),AdmApdwJ1(nJ1,1:n),ImpApdwJ0(nJ0,1:n),ImpApdwJ1(nJ1,1:n))
 	allocate(RTEdwJ0(nJ0,0:n),RTEdwJ1(nJ1,0:n),RTMdwJ0(nJ0,0:n),RTMdwJ1(nJ1,0:n))
 
 	do i=n,1,-1
 		if (i==n) then
 		AdmApdwJ0(:,i)=AdmIntJ0(:,i)
 		AdmApdwJ1(:,i)=AdmIntJ1(:,i)
 		ImpApdwJ0(:,i)=ImpIntJ0(:,i)
 		ImpApdwJ1(:,i)=ImpIntJ1(:,i)
 		RTEdwJ0(:,i)=(0.d0,0.d0)
 		RTEdwJ1(:,i)=(0.d0,0.d0)
 		RTMdwJ0(:,i)=(0.d0,0.d0)
 		RTMdwJ1(:,i)=(0.d0,0.d0)
 		else
 		AdmApdwJ0(:,i)=AdmIntJ0(:,i)*(AdmApdwJ0(:,i+1)+AdmIntJ0(:,i)* &
 			tghJ0(:,i))/(AdmIntJ0(:,i)+AdmApdwJ0(:,i+1)*tghJ0(:,i))
 		AdmApdwJ1(:,i)=AdmIntJ1(:,i)*(AdmApdwJ1(:,i+1)+AdmIntJ1(:,i)* &
 			tghJ1(:,i))/(AdmIntJ1(:,i)+AdmApdwJ1(:,i+1)*tghJ1(:,i))
 		ImpApdwJ0(:,i)=ImpIntJ0(:,i)*(ImpApdwJ0(:,i+1)+ImpIntJ0(:,i)* &
 			tghJ0(:,i))/(ImpIntJ0(:,i)+ImpApdwJ0(:,i+1)*tghJ0(:,i))
 		ImpApdwJ1(:,i)=ImpIntJ1(:,i)*(ImpApdwJ1(:,i+1)+ImpIntJ1(:,i)* &
 			tghJ1(:,i))/(ImpIntJ1(:,i)+ImpApdwJ1(:,i+1)*tghJ1(:,i))
 		RTEdwJ0(:,i)=(AdmIntJ0(:,i)-AdmApdwJ0(:,i+1))/(AdmIntJ0(:,i)+AdmApdwJ0(:,i+1))
 		RTEdwJ1(:,i)=(AdmIntJ1(:,i)-AdmApdwJ1(:,i+1))/(AdmIntJ1(:,i)+AdmApdwJ1(:,i+1))
 		RTMdwJ0(:,i)=(ImpIntJ0(:,i)-ImpApdwJ0(:,i+1))/(ImpIntJ0(:,i)+ImpApdwJ0(:,i+1))
 		RTMdwJ1(:,i)=(ImpIntJ1(:,i)-ImpApdwJ1(:,i+1))/(ImpIntJ1(:,i)+ImpApdwJ1(:,i+1))
 		end if	
 	end do
 		RTEdwJ0(:,0)=(AdmIntJ0(:,0)-AdmApdwJ0(:,1))/(AdmIntJ0(:,0)+AdmApdwJ0(:,1))
 		RTEdwJ1(:,0)=(AdmIntJ1(:,0)-AdmApdwJ1(:,1))/(AdmIntJ1(:,0)+AdmApdwJ1(:,1))
 		RTMdwJ0(:,0)=(ImpIntJ0(:,0)-ImpApdwJ0(:,1))/(ImpIntJ0(:,0)+ImpApdwJ0(:,1))
 		RTMdwJ1(:,0)=(ImpIntJ1(:,0)-ImpApdwJ1(:,1))/(ImpIntJ1(:,0)+ImpApdwJ1(:,1))
 
 	allocate(AdmApupJ0(nJ0,0:n-1),AdmApupJ1(nJ1,0:n-1),ImpApupJ0(nJ0,0:n-1),ImpApupJ1(nJ1,0:n-1))
 	allocate(RTEupJ0(nJ0,0:n),RTEupJ1(nJ1,0:n),RTMupJ0(nJ0,0:n),RTMupJ1(nJ1,0:n))
 
 	do i=0,n-1
 		if (i==0) then
 		AdmApupJ0(:,i)=AdmIntJ0(:,i)
 		AdmApupJ1(:,i)=AdmIntJ1(:,i)
 		ImpApupJ0(:,i)=ImpIntJ0(:,i)
 		ImpApupJ1(:,i)=ImpIntJ1(:,i)
 		RTEupJ0(:,i)=(0.d0,0.d0)
 		RTEupJ1(:,i)=(0.d0,0.d0)
 		RTMupJ0(:,i)=(0.d0,0.d0)
 		RTMupJ1(:,i)=(0.d0,0.d0)
 		else
 		AdmApupJ0(:,i)=AdmIntJ0(:,i)*(AdmApupJ0(:,i-1)+AdmIntJ0(:,i)* &
 			tghJ0(:,i))/(AdmIntJ0(:,i)+AdmApupJ0(:,i-1)*tghJ0(:,i))
 		AdmApupJ1(:,i)=AdmIntJ1(:,i)*(AdmApupJ1(:,i-1)+AdmIntJ1(:,i)* &
 			tghJ1(:,i))/(AdmIntJ1(:,i)+AdmApupJ1(:,i-1)*tghJ1(:,i))
 		ImpApupJ0(:,i)=ImpIntJ0(:,i)*(ImpApupJ0(:,i-1)+ImpIntJ0(:,i)* &
 			tghJ0(:,i))/(ImpIntJ0(:,i)+ImpApupJ0(:,i-1)*tghJ0(:,i))
 		ImpApupJ1(:,i)=ImpIntJ1(:,i)*(ImpApupJ1(:,i-1)+ImpIntJ1(:,i)* &
 			tghJ1(:,i))/(ImpIntJ1(:,i)+ImpApupJ1(:,i-1)*tghJ1(:,i))
 		RTEupJ0(:,i)=(AdmIntJ0(:,i)-AdmApupJ0(:,i-1))/(AdmIntJ0(:,i)+AdmApupJ0(:,i-1))
 		RTEupJ1(:,i)=(AdmIntJ1(:,i)-AdmApupJ1(:,i-1))/(AdmIntJ1(:,i)+AdmApupJ1(:,i-1))
 		RTMupJ0(:,i)=(ImpIntJ0(:,i)-ImpApupJ0(:,i-1))/(ImpIntJ0(:,i)+ImpApupJ0(:,i-1))
 		RTMupJ1(:,i)=(ImpIntJ1(:,i)-ImpApupJ1(:,i-1))/(ImpIntJ1(:,i)+ImpApupJ1(:,i-1))
 		end if
 	end do
 		RTEupJ0(:,n)=(AdmIntJ0(:,n)-AdmApupJ0(:,n-1))/(AdmIntJ0(:,n)+AdmApupJ0(:,n-1))
 		RTEupJ1(:,n)=(AdmIntJ1(:,n)-AdmApupJ1(:,n-1))/(AdmIntJ1(:,n)+AdmApupJ1(:,n-1))
 		RTMupJ0(:,n)=(ImpIntJ0(:,n)-ImpApupJ0(:,n-1))/(ImpIntJ0(:,n)+ImpApupJ0(:,n-1))
 		RTMupJ1(:,n)=(ImpIntJ1(:,n)-ImpApupJ1(:,n-1))/(ImpIntJ1(:,n)+ImpApupJ1(:,n-1))
 
 	allocate(AMdwJ0(nJ0),AMdwJ1(nJ1),AMupJ0(nJ0),AMupJ1(nJ1),FEdwJ0(nJ0),FEdwJ1(nJ1),FEupJ0(nJ0),FEupJ1(nJ1))

 	AMdwJ0=(cdexp(-uJ0(:,camadT)*(prof(camadT)-h0)) - RTMupJ0(:,camadT)* &
				cdexp(uJ0(:,camadT)*(prof(camadT-1)-h(camadT)-h0))) / &
				(1.d0-RTMupJ0(:,camadT)*RTMdwJ0(:,camadT)*cdexp(-2.d0*uhJ0(:,camadT)))
 	AMdwJ1=(cdexp(-uJ1(:,camadT)*(prof(camadT)-h0)) - RTMupJ1(:,camadT)* &
 				cdexp(uJ1(:,camadT)*(prof(camadT-1)-h(camadT)-h0))) / &
 				(1.d0-RTMupJ1(:,camadT)*RTMdwJ1(:,camadT)*cdexp(-2.d0*uhJ1(:,camadT)))
 	AMupJ0=(cdexp(uJ0(:,camadT)*(prof(camadT-1)-h0))-RTMdwJ0(:,camadT)* &
 				cdexp(-uJ0(:,camadT)*(prof(camadT)+h(camadT)-h0))) / &
 				(1.d0-RTMupJ0(:,camadT)*RTMdwJ0(:,camadT)*cdexp(-2.d0*uhJ0(:,camadT)))
 	AMupJ1=(cdexp(uJ1(:,camadT)*(prof(camadT-1)-h0))-RTMdwJ1(:,camadT)* &
 				cdexp(-uJ1(:,camadT)*(prof(camadT)+h(camadT)-h0))) / &
 				(1.d0-RTMupJ1(:,camadT)*RTMdwJ1(:,camadT)*cdexp(-2.d0*uhJ1(:,camadT)))
 	FEdwJ0=(cdexp(-uJ0(:,camadT)*(prof(camadT)-h0))+RTEupJ0(:,camadT)* &
 				cdexp(uJ0(:,camadT)*(prof(camadT-1)-h(camadT)-h0))) / &
 				(1.d0-RTEupJ0(:,camadT)*RTEdwJ0(:,camadT)*cdexp(-2.d0*uhJ0(:,camadT)))
 	FEdwJ1=(cdexp(-uJ1(:,camadT)*(prof(camadT)-h0))+RTEupJ1(:,camadT)* &
 				cdexp(uJ1(:,camadT)*(prof(camadT-1)-h(camadT)-h0))) / &
 				(1.d0-RTEupJ1(:,camadT)*RTEdwJ1(:,camadT)*cdexp(-2.d0*uhJ1(:,camadT)))
 	FEupJ0=(cdexp(uJ0(:,camadT)*(prof(camadT-1)-h0))+RTEdwJ0(:,camadT)* &
 				cdexp(-uJ0(:,camadT)*(prof(camadT)+h(camadT)-h0))) / &
 				(1.d0-RTEupJ0(:,camadT)*RTEdwJ0(:,camadT)*cdexp(-2.d0*uhJ0(:,camadT)))
 	FEupJ1=(cdexp(uJ1(:,camadT)*(prof(camadT-1)-h0))+RTEdwJ1(:,camadT)* &
 				cdexp(-uJ1(:,camadT)*(prof(camadT)+h(camadT)-h0))) / &
 				(1.d0-RTEupJ1(:,camadT)*RTEdwJ1(:,camadT)*cdexp(-2.d0*uhJ1(:,camadT)))

 	if (camad > camadT) then
 	allocate(TMdwJ0(nJ0,camadT:camad),TMdwJ1(nJ1,camadT:camad),TEdwJ0(nJ0,camadT:camad),TEdwJ1(nJ1,camadT:camad))
 		do j=camadT,camad
 			if (j == camadT) then
 			TMdwJ0(:,j)= - Iw * dsx / 2.d0
 			TMdwJ1(:,j)= - Iw * dsx / 2.d0
 			TEdwJ0(:,j)= - Iw * dsx / (2.d0 * AdmIntJ0(:,camadT))
 			TEdwJ1(:,j)= - Iw * dsx / (2.d0 * AdmIntJ1(:,camadT))
 			elseif (j == (camadT + 1) .and. j == n) then
 			TMdwJ0(:,j)=TMdwJ0(:,j-1)*(cdexp(-uJ0(:,camadT)*(prof(camadT)-h0)) - &
				RTMupJ0(:,camadT)*AMupJ0(:)*cdexp(-uhJ0(:,camadT)) + &
				RTMdwJ0(:,camadT)*AMdwJ0(:))
 			TMdwJ1(:,j)=TMdwJ1(:,j-1)*(cdexp(-uJ1(:,camadT)*(prof(camadT)-h0)) - &
				RTMupJ1(:,camadT)*AMupJ1(:)*cdexp(-uhJ1(:,camadT)) + &
				RTMdwJ1(:,camadT)*AMdwJ1(:))
 			TEdwJ0(:,j)=TEdwJ0(:,j-1)*(cdexp(-uJ0(:,camadT)*(prof(camadT)-h0)) + &
				RTEupJ0(:,camadT)*FEupJ0(:)*cdexp(-uhJ0(:,camadT)) + &
				RTEdwJ0(:,camadT)*FEdwJ0(:))
			TEdwJ1(:,j)=TEdwJ1(:,j-1)*(cdexp(-uJ1(:,camadT)*(prof(camadT)-h0)) + &
				RTEupJ1(:,camadT)*FEupJ1(:)*cdexp(-uhJ1(:,camadT)) + &
				RTEdwJ1(:,camadT)*FEdwJ1(:))
 			elseif (j == (camadT + 1) .and. j /= n) then
 			TMdwJ0(:,j)=TMdwJ0(:,j-1)*(cdexp(-uJ0(:,camadT)*(prof(camadT)-h0)) - &
				RTMupJ0(:,camadT)*AMupJ0(:)*cdexp(-uhJ0(:,camadT)) + &
				RTMdwJ0(:,camadT)*AMdwJ0(:)) / (1.d0 + RTMdwJ0(:,j)*cdexp(-2.d0*uhJ0(:,j)))
 			TMdwJ1(:,j)=TMdwJ1(:,j-1)*(cdexp(-uJ1(:,camadT)*(prof(camadT)-h0)) - &
				RTMupJ1(:,camadT)*AMupJ1(:)*cdexp(-uhJ1(:,camadT)) + &
				RTMdwJ1(:,camadT)*AMdwJ1(:)) / (1.d0 + RTMdwJ1(:,j)*cdexp(-2.d0*uhJ1(:,j)))
 			TEdwJ0(:,j)=TEdwJ0(:,j-1)*(cdexp(-uJ0(:,camadT)*(prof(camadT)-h0)) + &
				RTEupJ0(:,camadT)*FEupJ0(:)*cdexp(-uhJ0(:,camadT)) + &
				RTEdwJ0(:,camadT)*FEdwJ0(:)) / (1.d0 + RTEdwJ0(:,j)*cdexp(-2.d0*uhJ0(:,j)))
 			TEdwJ1(:,j)=TEdwJ1(:,j-1)*(cdexp(-uJ1(:,camadT)*(prof(camadT)-h0)) + &
				RTEupJ1(:,camadT)*FEupJ1(:)*cdexp(-uhJ1(:,camadT)) + &
				RTEdwJ1(:,camadT)*FEdwJ1(:)) / (1.d0 + RTEdwJ1(:,j)*cdexp(-2.d0*uhJ1(:,j)))
 			elseif (j /= n) then
 			TMdwJ0(:,j)=TMdwJ0(:,j-1)*(1.d0+RTMdwJ0(:,j-1))*cdexp(-uhJ0(:,j-1))/ &
 						(1.d0 + RTMdwJ0(:,j)*cdexp(-2.d0*uhJ0(:,j)))
 			TMdwJ1(:,j)=TMdwJ1(:,j-1)*(1.d0+RTMdwJ1(:,j-1))*cdexp(-uhJ1(:,j-1))/ &
 						(1.d0 + RTMdwJ1(:,j)*cdexp(-2.d0*uhJ1(:,j)))
 			TEdwJ0(:,j)=TEdwJ0(:,j-1)*(1.d0+RTEdwJ0(:,j-1))*cdexp(-uhJ0(:,j-1))/ &
 						(1.d0 + RTEdwJ0(:,j)*cdexp(-2.d0*uhJ0(:,j)))
 			TEdwJ1(:,j)=TEdwJ1(:,j-1)*(1.d0+RTEdwJ1(:,j-1))*cdexp(-uhJ1(:,j-1))/ &
 						(1.d0 + RTEdwJ1(:,j)*cdexp(-2.d0*uhJ1(:,j)))
 			elseif (j==n) then
 			TMdwJ0(:,j)=TMdwJ0(:,j-1)*(1.d0+RTMdwJ0(:,j-1))*cdexp(-uhJ0(:,j-1))
 			TMdwJ1(:,j)=TMdwJ1(:,j-1)*(1.d0+RTMdwJ1(:,j-1))*cdexp(-uhJ1(:,j-1))
 			TEdwJ0(:,j)=TEdwJ0(:,j-1)*(1.d0+RTEdwJ0(:,j-1))*cdexp(-uhJ0(:,j-1))
 			TEdwJ1(:,j)=TEdwJ1(:,j-1)*(1.d0+RTEdwJ1(:,j-1))*cdexp(-uhJ1(:,j-1))
 			end if
 		end do
 	elseif (camad < camadT) then
 	allocate(TMupJ0(nJ0,camad:camadT),TMupJ1(nJ1,camad:camadT),TEupJ0(nJ0,camad:camadT),TEupJ1(nJ1,camad:camadT))
 		do j=camadT,camad,-1
 			if (j == camadT) then
 			TMupJ0(:,j)= Iw * dsx / 2.d0
 			TMupJ1(:,j)= Iw * dsx / 2.d0
 			TEupJ0(:,j)= - Iw * dsx / (2.d0 * AdmIntJ0(:,camadT))
 			TEupJ1(:,j)= - Iw * dsx / (2.d0 * AdmIntJ1(:,camadT))
 			elseif (j == (camadT - 1) .and. j == 0) then
 			TMupJ0(:,j)=TMupJ0(:,j+1)*(cdexp(-uJ0(:,camadT)*h0) + &
				RTMupJ0(:,camadT)*AMupJ0(:) - RTMdwJ0(:,camadT)*AMdwJ0(:)* &
				cdexp(-uhJ0(:,camadT)))
 			TMupJ1(:,j)=TMupJ1(:,j+1)*(cdexp(-uJ1(:,camadT)*h0) + &
				RTMupJ1(:,camadT)*AMupJ1(:) - RTMdwJ1(:,camadT)*AMdwJ1(:)* &
				cdexp(-uhJ1(:,camadT)))
 			TEupJ0(:,j)=TEupJ0(:,j+1)*(cdexp(-uJ0(:,camadT)*h0) + &
				RTEupJ0(:,camadT)*FEupJ0(:) + RTEdwJ0(:,camadT)*FEdwJ0(:)* &
				cdexp(-uhJ0(:,camadT)))
 			TEupJ1(:,j)=TEupJ1(:,j+1)*(cdexp(-uJ1(:,camadT)*h0) + &
				RTEupJ1(:,camadT)*FEupJ1(:) + RTEdwJ1(:,camadT)*FEdwJ1(:)* &
				cdexp(-uhJ1(:,camadT)))
 			elseif (j == (camadT - 1) .and. j /= 0) then
 			TMupJ0(:,j)=TMupJ0(:,j+1)*(cdexp(uJ0(:,camadT)*(prof(camadT-1)-h0)) + &
				RTMupJ0(:,camadT)*AMupJ0(:) - RTMdwJ0(:,camadT)*AMdwJ0(:)* &
				cdexp(-uhJ0(:,camadT))) / (1.d0+RTMupJ0(:,j)*cdexp(-2.d0*uhJ0(:,j)))
 			TMupJ1(:,j)=TMupJ1(:,j+1)*(cdexp(uJ1(:,camadT)*(prof(camadT-1)-h0)) + &
				RTMupJ1(:,camadT)*AMupJ1(:) - RTMdwJ1(:,camadT)*AMdwJ1(:)* &
				cdexp(-uhJ1(:,camadT))) / (1.d0+RTMupJ1(:,j)*cdexp(-2.d0*uhJ1(:,j)))
 			TEupJ0(:,j)=TEupJ0(:,j+1)*(cdexp(uJ0(:,camadT)*(prof(camadT-1)-h0)) + &
				RTEupJ0(:,camadT)*FEupJ0(:) + RTEdwJ0(:,camadT)*FEdwJ0(:)* &
				cdexp(-uhJ0(:,camadT))) / (1.d0+RTEupJ0(:,j)*cdexp(-2.d0*uhJ0(:,j)))
 			TEupJ1(:,j)=TEupJ1(:,j+1)*(cdexp(uJ1(:,camadT)*(prof(camadT-1)-h0)) + &
				RTEupJ1(:,camadT)*FEupJ1(:) + RTEdwJ1(:,camadT)*FEdwJ1(:)* &
				cdexp(-uhJ1(:,camadT))) / (1.d0+RTEupJ1(:,j)*cdexp(-2.d0*uhJ1(:,j)))
 			elseif (j /= 0) then
 			TMupJ0(:,j)=TMupJ0(:,j+1)*(1.d0+RTMupJ0(:,j+1))*cdexp(-uhJ0(:,j+1)) / &
				(1.d0 + RTMupJ0(:,j)*cdexp(-2.d0*uhJ0(:,j)))
 			TMupJ1(:,j)=TMupJ1(:,j+1)*(1.d0+RTMupJ1(:,j+1))*cdexp(-uhJ1(:,j+1)) / &
				(1.d0 + RTMupJ1(:,j)*cdexp(-2.d0*uhJ1(:,j)))
 			TEupJ0(:,j)=TEupJ0(:,j+1)*(1.d0+RTEupJ0(:,j+1))*cdexp(-uhJ0(:,j+1)) / &
				(1.d0 + RTEupJ0(:,j)*cdexp(-2.d0*uhJ0(:,j)))
 			TEupJ1(:,j)=TEupJ1(:,j+1)*(1.d0+RTEupJ1(:,j+1))*cdexp(-uhJ1(:,j+1)) / &
				(1.d0 + RTEupJ1(:,j)*cdexp(-2.d0*uhJ1(:,j)))
 			elseif (j == 0) then
 			TMupJ0(:,j)=TMupJ0(:,1)*(1.d0+RTMupJ0(:,1))*cdexp(-uhJ0(:,1))
 			TMupJ1(:,j)=TMupJ1(:,1)*(1.d0+RTMupJ1(:,1))*cdexp(-uhJ1(:,1))
 			TEupJ0(:,j)=TEupJ0(:,1)*(1.d0+RTEupJ0(:,1))*cdexp(-uhJ0(:,1))
 			TEupJ1(:,j)=TEupJ1(:,1)*(1.d0+RTEupJ1(:,1))*cdexp(-uhJ1(:,1))
 			end if
 		end do
 	else
 	allocate(TMdwJ0(nJ0,camadT:camad),TMdwJ1(nJ1,camadT:camad),TEdwJ0(nJ0,camadT:camad),TEdwJ1(nJ1,camadT:camad))
 	allocate(TMupJ0(nJ0,camad:camadT),TMupJ1(nJ1,camad:camadT),TEupJ0(nJ0,camad:camadT),TEupJ1(nJ1,camad:camadT))
 			TMdwJ0(:,camad)= - Iw * dsx / 2.d0
 			TMdwJ1(:,camad)= - Iw * dsx / 2.d0
 			TEdwJ0(:,camad)= - Iw * dsx / (2.d0 * AdmIntJ0(:,camadT))
 			TEdwJ1(:,camad)= - Iw * dsx / (2.d0 * AdmIntJ1(:,camadT))
 			TMupJ0(:,camad)= Iw * dsx / 2.d0
 			TMupJ1(:,camad)= Iw * dsx / 2.d0
 			TEupJ0(:,camad)= - Iw * dsx / (2.d0 * AdmIntJ0(:,camadT))
 			TEupJ1(:,camad)= - Iw * dsx / (2.d0 * AdmIntJ1(:,camadT))
 	end if

 	allocate(Ktmdz_J0(nJ0),Ktmdz_J1(nJ1),Ktm_J0(nJ0),Ktm_J1(nJ1),Kte_J0(nJ0),Kte_J1(nJ1),Ktedz_J0(nJ0),Ktedz_J1(nJ1))
 	allocate(kernelExJ0(nJ0),kernelExJ1(nJ1),kernelEyJ0(nJ0),kernelEyJ1(nJ1),kernelEzJ1(nJ1))
 	allocate(kernelHxJ0(nJ0),kernelHxJ1(nJ1),kernelHyJ0(nJ0),kernelHyJ1(nJ1),kernelHzJ1(nJ1))
 	if (camad == 0 .and. camadT /= 0) then
 		Ktmdz_J0=(ImpIntJ0(:,0)*TMupJ0(:,0)*cdexp(uJ0(:,0)*z))*w_J0(:)
 		Ktmdz_J1=(ImpIntJ1(:,0)*TMupJ1(:,0)*cdexp(uJ1(:,0)*z))*w_J1(:)
 		Kte_J0=(TEupJ0(:,0)*cdexp(uJ0(:,0)*z))*w_J0(:)
 		Kte_J1=(TEupJ1(:,0)*cdexp(uJ1(:,0)*z))*w_J1(:)
 		Ktm_J0=(TMupJ0(:,0)*cdexp(uJ0(:,0)*z))*w_J0(:)
 		Ktm_J1=(TMupJ1(:,0)*cdexp(uJ1(:,0)*z))*w_J1(:)
 		Ktedz_J0=(AdmIntJ0(:,0)*TEupJ0(:,0)*cdexp(uJ0(:,0)*z))*w_J0(:)
 		Ktedz_J1=(AdmIntJ1(:,0)*TEupJ1(:,0)*cdexp(uJ1(:,0)*z))*w_J1(:)

 		kernelExJ1 = (2.d0*x**2/r**3 - 1.d0/r) * Ktmdz_J1 - (2.d0*y**2/r**3 - 1.d0/r) * Kte_J1
 		kernelExJ0 = (x**2/r**2 * Ktmdz_J0 - y**2/r**2 * Kte_J0) * krJ0
 		Ex_p = (sum(kernelExJ1) - sum(kernelExJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelEyJ1 = 2.d0*x*y/r**3 * (Ktmdz_J1 + Kte_J1)
 		kernelEyJ0 = x*y/r**2 * (Ktmdz_J0 + Kte_J0) * krJ0
 		Ey_p = (sum(kernelEyJ1) - sum(kernelEyJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelEzJ1 = x/(r*neta) * Ktm_J1 * krJ1 * krJ1
 		Ez_p = - sum(kernelEzJ1) / (2.d0*pi*r)				!este último r é decorrente do uso dos filtros

 		kernelHxJ1 = 2.d0*x*y/r**3 * (Ktm_J1 + Ktedz_J1)
 		kernelHxJ0 = x*y/r**2 * (Ktm_J0 + Ktedz_J0) * krJ0
 		Hx_p = (sum(kernelHxJ1) - sum(kernelHxJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros

 		kernelHyJ1 = (1.d0/r - 2.d0*x**2/r**3) * Ktm_J1 - (1.d0/r - 2.d0*y**2/r**3) * Ktedz_J1
 		kernelHyJ0 = (x**2/r**2 * Ktm_J0 - y**2/r**2 * Ktedz_J0) * krJ0
 		Hy_p = (sum(kernelHyJ1) + sum(kernelHyJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelHzJ1 = y/(r*zeta) * Kte_J1 * krJ1 * krJ1
 		Hz_p = - sum(kernelHzJ1) / (2.d0*pi*r)				!este último r é decorrente do uso dos filtros
 
 	elseif (camad < camadT) then !camada k
 		ktmdz_J0=(ImpIntJ0(:,camad)*TMupJ0(:,camad)*(cdexp(uJ0(:,camad)*(z-prof(camad))) - &
 			RTMupJ0(:,camad)*cdexp(-uJ0(:,camad)*(z-prof(camad-1)+h(camad)))))*w_J0(:)
 		ktmdz_J1=(ImpIntJ1(:,camad)*TMupJ1(:,camad)*(cdexp(uJ1(:,camad)*(z-prof(camad))) - &
 			RTMupJ1(:,camad)*cdexp(-uJ1(:,camad)*(z-prof(camad-1)+h(camad)))))*w_J1(:)
 		Kte_J0=(TEupJ0(:,camad)*(cdexp(uJ0(:,camad)*(z-prof(camad))) + &
 			 RTEupJ0(:,camad)*cdexp(-uJ0(:,camad)*(z-prof(camad-1)+h(camad)))))*w_J0(:)
 		Kte_J1=(TEupJ1(:,camad)*(cdexp(uJ1(:,camad)*(z-prof(camad))) + &
 			 RTEupJ1(:,camad)*cdexp(-uJ1(:,camad)*(z-prof(camad-1)+h(camad)))))*w_J1(:)
 		ktm_J0=(TMupJ0(:,camad)*(cdexp(uJ0(:,camad)*(z-prof(camad))) + &
 			RTMupJ0(:,camad)*cdexp(-uJ0(:,camad)*(z-prof(camad-1)+h(camad)))))*w_J0(:)
 		ktm_J1=(TMupJ1(:,camad)*(cdexp(uJ1(:,camad)*(z-prof(camad))) + &
 			RTMupJ1(:,camad)*cdexp(-uJ1(:,camad)*(z-prof(camad-1)+h(camad)))))*w_J1(:)
 		Ktedz_J0=(AdmIntJ0(:,camad)*TEupJ0(:,camad)*(cdexp(uJ0(:,camad)*(z-prof(camad))) - &
 			RTEupJ0(:,camad)*cdexp(-uJ0(:,camad)*(z-prof(camad-1)+h(camad)))))*w_J0(:)
 		Ktedz_J1=(AdmIntJ1(:,camad)*TEupJ1(:,camad)*(cdexp(uJ1(:,camad)*(z-prof(camad))) - &
 			RTEupJ1(:,camad)*cdexp(-uJ1(:,camad)*(z-prof(camad-1)+h(camad)))))*w_J1(:)
 
 		kernelExJ1 = (2.d0*x**2/r**3 - 1.d0/r) * Ktmdz_J1 - (2.d0*y**2/r**3 - 1.d0/r) * Kte_J1
 		kernelExJ0 = (x**2/r**2 * Ktmdz_J0 - y**2/r**2 * Kte_J0) * krJ0
 		Ex_p = (sum(kernelExJ1) - sum(kernelExJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelEyJ1 = 2.d0*x*y/r**3 * (Ktmdz_J1 + Kte_J1)
 		kernelEyJ0 = x*y/r**2 * (Ktmdz_J0 + Kte_J0) * krJ0
 		Ey_p = (sum(kernelEyJ1) - sum(kernelEyJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelEzJ1 = x/(r*condut(camad)) * Ktm_J1 * krJ1 * krJ1
 		Ez_p = - sum(kernelEzJ1) / (2.d0*pi*r)				!este último r é decorrente do uso dos filtros
 
 		kernelHxJ1 = 2.d0*x*y/r**3 * (Ktm_J1 + Ktedz_J1)
 		kernelHxJ0 = (x*y/r**2 * (Ktm_J0 + Ktedz_J0)) * krJ0
 		Hx_p = (sum(kernelHxJ1) - sum(kernelHxJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelHyJ1 = (1.d0/r - 2.d0*x**2/r**3) * Ktm_J1 - (1.d0/r - 2.d0*y**2/r**3) * Ktedz_J1
 		kernelHyJ0 = (x**2/r**2 * Ktm_J0 - y**2/r**2 * Ktedz_J0) * krJ0
 		Hy_p = (sum(kernelHyJ1) + sum(kernelHyJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelHzJ1 = y/(r*zeta) * Kte_J1 * krJ1 * krJ1
 		Hz_p = - sum(kernelHzJ1) / (2.d0*pi*r)				!este último r é decorrente do uso dos filtros
 
 	elseif (camad == camadT .and. z <= h0) then	!na mesma camada do transmissor mas acima dele
		Ktmdz_J0=(ImpIntJ0(:,camad)*TMupJ0(:,camad)*(cdexp(uJ0(:,camad)*(z-h0)) - &
 			RTMupJ0(:,camad)*AMupJ0(:)*cdexp(-uJ0(:,camad)*(z-prof(camad-1))) - &
 			RTMdwJ0(:,camad)*AMdwJ0(:)*cdexp(uJ0(:,camad)*(z-prof(camad)))))*w_J0(:)
 		Ktmdz_J1=(ImpIntJ1(:,camad)*TMupJ1(:,camad)*(cdexp(uJ1(:,camad)*(z-h0)) - &
 			RTMupJ1(:,camad)*AMupJ1(:)*cdexp(-uJ1(:,camad)*(z-prof(camad-1))) - &
 			RTMdwJ1(:,camad)*AMdwJ1(:)*cdexp(uJ1(:,camad)*(z-prof(camad)))))*w_J1(:)
 		Kte_J0=(TEupJ0(:,camad)*(cdexp(uJ0(:,camad)*(z-h0)) + &
 			RTEupJ0(:,camad)*FEupJ0(:)*cdexp(-uJ0(:,camad)*(z-prof(camad-1))) + &
 			RTEdwJ0(:,camad)*FEdwJ0(:)*cdexp(uJ0(:,camad)*(z-prof(camad)))))*w_J0(:)
 		Kte_J1=(TEupJ1(:,camad)*(cdexp(uJ1(:,camad)*(z-h0)) + &
 			RTEupJ1(:,camad)*FEupJ1(:)*cdexp(-uJ1(:,camad)*(z-prof(camad-1))) + &
 			RTEdwJ1(:,camad)*FEdwJ1(:)*cdexp(uJ1(:,camad)*(z-prof(camad)))))*w_J1(:)
 		Ktm_J0=(TMupJ0(:,camad)*(cdexp(uJ0(:,camad)*(z-h0)) + &
 			RTMupJ0(:,camad)*AMupJ0(:)*cdexp(-uJ0(:,camad)*(z-prof(camad-1))) - &
 			RTMdwJ0(:,camad)*AMdwJ0(:)*cdexp(uJ0(:,camad)*(z-prof(camad)))))*w_J0(:)
 		Ktm_J1=(TMupJ1(:,camad)*(cdexp(uJ1(:,camad)*(z-h0)) + &
 			RTMupJ1(:,camad)*AMupJ1(:)*cdexp(-uJ1(:,camad)*(z-prof(camad-1))) - &
 			RTMdwJ1(:,camad)*AMdwJ1(:)*cdexp(uJ1(:,camad)*(z-prof(camad)))))*w_J1(:)
 		Ktedz_J0=(AdmIntJ0(:,camad)*TEupJ0(:,camad)*(cdexp(uJ0(:,camad)*(z-h0)) - &
 			RTEupJ0(:,camad)*FEupJ0(:)*cdexp(-uJ0(:,camad)*(z-prof(camad-1))) + &
 			RTEdwJ0(:,camad)*FEdwJ0(:)*cdexp(uJ0(:,camad)*(z-prof(camad)))))*w_J0(:)
 		Ktedz_J1=(AdmIntJ1(:,camad)*TEupJ1(:,camad)*(cdexp(uJ1(:,camad)*(z-h0)) - &
 			RTEupJ1(:,camad)*FEupJ1(:)*cdexp(-uJ1(:,camad)*(z-prof(camad-1))) + &
 			RTEdwJ1(:,camad)*FEdwJ1(:)*cdexp(uJ1(:,camad)*(z-prof(camad)))))*w_J1(:)

 		kernelExJ1 = (2.d0*x**2/r**3 - 1.d0/r) * Ktmdz_J1 - (2.d0*y**2/r**3 - 1.d0/r) * Kte_J1
 		kernelExJ0 = (x**2/r**2 * Ktmdz_J0 - y**2/r**2 * Kte_J0) * krJ0
 		Ex_p = (sum(kernelExJ1) - sum(kernelExJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelEyJ1 = 2.d0*x*y/r**3 * (Ktmdz_J1 + Kte_J1)
 		kernelEyJ0 = x*y/r**2 * (Ktmdz_J0 + Kte_J0) * krJ0
 		Ey_p = (sum(kernelEyJ1) - sum(kernelEyJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
		if (camad /= 0) then
 		kernelEzJ1 = x/(r*condut(camad)) * Ktm_J1 * krJ1 * krJ1
		else
 		kernelEzJ1 = x/(r*neta) * Ktm_J1 * krJ1 * krJ1
		end if
 		Ez_p = - sum(kernelEzJ1) / (2.d0*pi*r)						!este último r é decorrente do uso dos filtros
 
 		kernelHxJ1 = 2.d0*x*y/r**3 * (Ktm_J1 + Ktedz_J1)
 		kernelHxJ0 = x*y/r**2 * (Ktm_J0 + Ktedz_J0) * krJ0
 		Hx_p = (sum(kernelHxJ1) - sum(kernelHxJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelHyJ1 = (1.d0/r - 2.d0*x**2/r**3) * Ktm_J1 - (1.d0/r - 2.d0*y**2/r**3) * Ktedz_J1
 		kernelHyJ0 = (x**2/r**2 * Ktm_J0 - y**2/r**2 * Ktedz_J0) * krJ0
 		Hy_p = (sum(kernelHyJ1) + sum(kernelHyJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelHzJ1 = y/(r*zeta) * Kte_J1 * krJ1 * krJ1
 		Hz_p = - sum(kernelHzJ1) / (2.d0*pi*r)						!este último r é decorrente do uso dos filtros
 
 	elseif (camad == camadT .and. z > h0) then	!na mesma camada do transmissor mas abaixo dele
 		Ktmdz_J0=(ImpIntJ0(:,camad)*TMdwJ0(:,camad)*(cdexp(-uJ0(:,camad)*(z-h0)) - &
 			RTMupJ0(:,camad)*AMupJ0(:)*cdexp(-uJ0(:,camad)*(z-prof(camad-1))) - &
 			RTMdwJ0(:,camad)*AMdwJ0(:)*cdexp(uJ0(:,camad)*(z-prof(camad)))))*w_J0(:)
 		Ktmdz_J1=(ImpIntJ1(:,camad)*TMdwJ1(:,camad)*(cdexp(-uJ1(:,camad)*(z-h0)) - &
 			RTMupJ1(:,camad)*AMupJ1(:)*cdexp(-uJ1(:,camad)*(z-prof(camad-1))) - &
 			RTMdwJ1(:,camad)*AMdwJ1(:)*cdexp(uJ1(:,camad)*(z-prof(camad)))))*w_J1(:)
 		Kte_J0=(TEdwJ0(:,camad)*(cdexp(-uJ0(:,camad)*(z-h0)) + &
 			RTEupJ0(:,camad)*FEupJ0(:)*cdexp(-uJ0(:,camad)*(z-prof(camad-1))) + &
 			RTEdwJ0(:,camad)*FEdwJ0(:)*cdexp(uJ0(:,camad)*(z-prof(camad)))))*w_J0(:)
 		Kte_J1=(TEdwJ1(:,camad)*(cdexp(-uJ1(:,camad)*(z-h0)) + &
 			RTEupJ1(:,camad)*FEupJ1(:)*cdexp(-uJ1(:,camad)*(z-prof(camad-1))) + &
 			RTEdwJ1(:,camad)*FEdwJ1(:)*cdexp(uJ1(:,camad)*(z-prof(camad)))))*w_J1(:)
		Ktm_J0=(TMdwJ0(:,camad)*(cdexp(-uJ0(:,camad)*(z-h0)) - &
 			RTMupJ0(:,camad)*AMupJ0(:)*cdexp(-uJ0(:,camad)*(z-prof(camad-1))) + &
 			RTMdwJ0(:,camad)*AMdwJ0(:)*cdexp(uJ0(:,camad)*(z-prof(camad)))))*w_J0(:)
 		Ktm_J1=(TMdwJ1(:,camad)*(cdexp(-uJ1(:,camad)*(z-h0)) - &
 			RTMupJ1(:,camad)*AMupJ1(:)*cdexp(-uJ1(:,camad)*(z-prof(camad-1))) + &
 			RTMdwJ1(:,camad)*AMdwJ1(:)*cdexp(uJ1(:,camad)*(z-prof(camad)))))*w_J1(:)
 		Ktedz_J0=(AdmIntJ0(:,camad)*TEdwJ0(:,camad)*(cdexp(-uJ0(:,camad)*(z-h0)) + &
 			RTEupJ0(:,camad)*FEupJ0(:)*cdexp(-uJ0(:,camad)*(z-prof(camad-1))) - &
 			RTEdwJ0(:,camad)*FEdwJ0(:)*cdexp(uJ0(:,camad)*(z-prof(camad)))))*w_J0(:)
 		Ktedz_J1=(AdmIntJ1(:,camad)*TEdwJ1(:,camad)*(cdexp(-uJ1(:,camad)*(z-h0)) + &
 			RTEupJ1(:,camad)*FEupJ1(:)*cdexp(-uJ1(:,camad)*(z-prof(camad-1))) - &
 			RTEdwJ1(:,camad)*FEdwJ1(:)*cdexp(uJ1(:,camad)*(z-prof(camad)))))*w_J1(:)
 
 		kernelExJ1 = (1.d0/r - 2.d0*x**2/r**3) * Ktmdz_J1 + (1.d0/r - 2.d0*y**2/r**3) * Kte_J1
 		kernelExJ0 = (x**2/r**2 * Ktmdz_J0 + y**2/r**2 * Kte_J0) * krJ0
 		Ex_p = (sum(kernelExJ1) + sum(kernelExJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelEyJ1 = 2.d0*x*y/r**3 * (-Ktmdz_J1 + Kte_J1)
 		kernelEyJ0 = x*y/r**2 * (-Ktmdz_J0 + Kte_J0) * krJ0
 		Ey_p = (sum(kernelEyJ1) - sum(kernelEyJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
		if (camad /= 0) then
 		kernelEzJ1 = x/(r*condut(camad)) * Ktm_J1 * krJ1 * krJ1
		else
 		kernelEzJ1 = x/(r*neta) * Ktm_J1 * krJ1 * krJ1
		end if
 		Ez_p = - sum(kernelEzJ1) / (2.d0*pi*r)						!este último r é decorrente do uso dos filtros
 
 		kernelHxJ1 = 2.d0*x*y/r**3 * (Ktm_J1 - Ktedz_J1)
 		kernelHxJ0 = (x*y/r**2 * (Ktm_J0 - Ktedz_J0)) * krJ0
 		Hx_p = (sum(kernelHxJ1) - sum(kernelHxJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelHyJ1 = (1.d0/r - 2.d0*x**2/r**3) * Ktm_J1 + (1.d0/r - 2.d0*y**2/r**3) * Ktedz_J1
 		kernelHyJ0 = (x**2/r**2 * Ktm_J0 + y**2/r**2 * Ktedz_J0) * krJ0
 		Hy_p = (sum(kernelHyJ1) + sum(kernelHyJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelHzJ1 = y/(r*zeta) * Kte_J1 * krJ1 * krJ1
 		Hz_p = - sum(kernelHzJ1) / (2.d0*pi*r)						!este último r é decorrente do uso dos filtros
 
 	elseif (camad > camadT .and. camad /= n) then !camada j
 		Ktmdz_J0=(ImpIntJ0(:,camad)*TMdwJ0(:,camad)*(cdexp(-uJ0(:,camad)*(z-prof(camad-1))) - &
 			RTMdwJ0(:,camad)*cdexp(uJ0(:,camad)*(z-prof(camad)-h(camad)))))*w_J0(:)
 		Ktmdz_J1=(ImpIntJ1(:,camad)*TMdwJ1(:,camad)*(cdexp(-uJ1(:,camad)*(z-prof(camad-1))) - &
 			RTMdwJ1(:,camad)*cdexp(uJ1(:,camad)*(z-prof(camad)-h(camad)))))*w_J1(:)
 		Kte_J0=(TEdwJ0(:,camad)*(cdexp(-uJ0(:,camad)*(z-prof(camad-1))) + &
 			RTEdwJ0(:,camad)*cdexp(uJ0(:,camad)*(z-prof(camad)-h(camad)))))*w_J0(:)
 		Kte_J1=(TEdwJ1(:,camad)*(cdexp(-uJ1(:,camad)*(z-prof(camad-1))) + &
 			RTEdwJ1(:,camad)*cdexp(uJ1(:,camad)*(z-prof(camad)-h(camad)))))*w_J1(:)
 		Ktm_J0=(TMdwJ0(:,camad)*(cdexp(-uJ0(:,camad)*(z-prof(camad-1))) + &
 			RTMdwJ0(:,camad)*cdexp(uJ0(:,camad)*(z-prof(camad)-h(camad)))))*w_J0(:)
 		Ktm_J1=(TMdwJ1(:,camad)*(cdexp(-uJ1(:,camad)*(z-prof(camad-1))) + &
 			RTMdwJ1(:,camad)*cdexp(uJ1(:,camad)*(z-prof(camad)-h(camad)))))*w_J1(:)
 		Ktedz_J0=(AdmIntJ0(:,camad)*TEdwJ0(:,camad)*(cdexp(-uJ0(:,camad)*(z-prof(camad-1))) - &
 			RTEdwJ0(:,camad)*cdexp(uJ0(:,camad)*(z-prof(camad)-h(camad)))))*w_J0(:)
 		Ktedz_J1=(AdmIntJ1(:,camad)*TEdwJ1(:,camad)*(cdexp(-uJ1(:,camad)*(z-prof(camad-1))) - &
 			RTEdwJ1(:,camad)*cdexp(uJ1(:,camad)*(z-prof(camad)-h(camad)))))*w_J1(:)
 
 		kernelExJ1 = (1.d0/r - 2.d0*x**2/r**3) * Ktmdz_J1 + (1.d0/r - 2.d0*y**2/r**3) * Kte_J1
 		kernelExJ0 = (x**2/r**2 * Ktmdz_J0 + y**2/r**2 * Kte_J0) * krJ0
 		Ex_p = (sum(kernelExJ1) + sum(kernelExJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelEyJ1 = 2.d0*x*y/r**3 * (-Ktmdz_J1 + Kte_J1)
 		kernelEyJ0 = x*y/r**2 * (-Ktmdz_J0 + Kte_J0) * krJ0
 		Ey_p = (sum(kernelEyJ1) - sum(kernelEyJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelEzJ1 = x/(r*condut(camad)) * Ktm_J1 * krJ1 * krJ1
 		Ez_p = - sum(kernelEzJ1) / (2.d0*pi*r)						!este último r é decorrente do uso dos filtros
 
 		kernelHxJ1 = 2.d0*x*y/r**3 * (Ktm_J1 - Ktedz_J1)
 		kernelHxJ0 = x*y/r**2 * (Ktm_J0 - Ktedz_J0) * krJ0
 		Hx_p = (sum(kernelHxJ1) - sum(kernelHxJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelHyJ1 = (1.d0/r - 2.d0*x**2/r**3) * Ktm_J1 + (1.d0/r - 2.d0*y**2/r**3) * Ktedz_J1
 		kernelHyJ0 = (x**2/r**2 * Ktm_J0 + y**2/r**2 * Ktedz_J0) * krJ0
 		Hy_p = (sum(kernelHyJ1) + sum(kernelHyJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelHzJ1 = y/(r*zeta) * Kte_J1 * krJ1 * krJ1
 		Hz_p = - sum(kernelHzJ1) / (2.d0*pi*r)						!este último r é decorrente do uso dos filtros
 	else	!camada n
 		Ktmdz_J0=(ImpIntJ0(:,n)*TMdwJ0(:,n)*cdexp(-uJ0(:,n)*(z-prof(n-1))))*w_J0(:)
 		Ktmdz_J1=(ImpIntJ1(:,n)*TMdwJ1(:,n)*cdexp(-uJ1(:,n)*(z-prof(n-1))))*w_J1(:)
 		Kte_J0=(TEdwJ0(:,n)*cdexp(-uJ0(:,n)*(z-prof(n-1))))*w_J0(:)
 		Kte_J1=(TEdwJ1(:,n)*cdexp(-uJ1(:,n)*(z-prof(n-1))))*w_J1(:)
 		Ktm_J0=(TMdwJ0(:,n)*cdexp(-uJ0(:,n)*(z-prof(n-1))))*w_J0(:)
 		Ktm_J1=(TMdwJ1(:,n)*cdexp(-uJ1(:,n)*(z-prof(n-1))))*w_J1(:)
 		Ktedz_J0=(AdmIntJ0(:,n)*TEdwJ0(:,n)*cdexp(-uJ0(:,n)*(z-prof(n-1))))*w_J0(:)
 		Ktedz_J1=(AdmIntJ1(:,n)*TEdwJ1(:,n)*cdexp(-uJ1(:,n)*(z-prof(n-1))))*w_J1(:)
 
 		kernelExJ1 = (1.d0/r - 2.d0*x**2/r**3) * Ktmdz_J1 + (1.d0/r - 2.d0*y**2/r**3) * Kte_J1
 		kernelExJ0 = (x**2/r**2 * Ktmdz_J0 + y**2/r**2 * Kte_J0) * krJ0
 		Ex_p = (sum(kernelExJ1) + sum(kernelExJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelEyJ1 = 2.d0*x*y/r**3 * (-Ktmdz_J1 + Kte_J1)
 		kernelEyJ0 = x*y/r**2 * (-Ktmdz_J0 + Kte_J0) * krJ0
 		Ey_p = (sum(kernelEyJ1) - sum(kernelEyJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelEzJ1 = x/(r*condut(n)) * Ktm_J1 * krJ1 * krJ1
 		Ez_p = - sum(kernelEzJ1) / (2.d0*pi*r)						!este último r é decorrente do uso dos filtros
 
 		kernelHxJ1 = 2.d0*x*y/r**3 * (Ktm_J1 - Ktedz_J1)
 		kernelHxJ0 = x*y/r**2 * (Ktm_J0 - Ktedz_J0) * krJ0
 		Hx_p = (sum(kernelHxJ1) - sum(kernelHxJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelHyJ1 = (1.d0/r - 2.d0*x**2/r**3) * Ktm_J1 + (1.d0/r - 2.d0*y**2/r**3) * Ktedz_J1
 		kernelHyJ0 = (x**2/r**2 * Ktm_J0 + y**2/r**2 * Ktedz_J0) * krJ0
 		Hy_p = (sum(kernelHyJ1) + sum(kernelHyJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelHzJ1 = y/(r*zeta) * Kte_J1 * krJ1 * krJ1
 		Hz_p = - sum(kernelHzJ1) / (2.d0*pi*r)						!este último r é decorrente do uso dos filtros
 
 	end if
 
 	deallocate(h,KrJ0,KrJ1,w_J0,w_J1)
 	deallocate(wvnb2,uJ0,uJ1,AdmIntJ0,AdmIntJ1,ImpIntJ0,ImpIntJ1,uhJ0,uhJ1,tghJ0,tghJ1)
 	deallocate(AdmApdwJ0,AdmApdwJ1,ImpApdwJ0,ImpApdwJ1,RTEdwJ0,RTEdwJ1,RTMdwJ0,RTMdwJ1)
 	deallocate(AdmApupJ0,AdmApupJ1,ImpApupJ0,ImpApupJ1,RTEupJ0,RTEupJ1,RTMupJ0,RTMupJ1)
 	deallocate(AMdwJ0,AMdwJ1,AMupJ0,AMupJ1,FEdwJ0,FEdwJ1,FEupJ0,FEupJ1)
 	deallocate(Ktmdz_J0,Ktmdz_J1,Ktm_J0,Ktm_J1,Kte_J0,Kte_J1,Ktedz_J0,Ktedz_J1)
 	deallocate(kernelExJ0,kernelExJ1,kernelEyJ0,kernelEyJ1,kernelEzJ1)
 	deallocate(kernelHxJ0,kernelHxJ1,kernelHyJ0,kernelHyJ1,kernelHzJ1)
	end subroutine dehx_xyz_loops
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
	subroutine dehy_xyz_loops( idtfcd_cJ0, nJ0, idtfcd_cJ1, nJ1, Tx, Ty, Iw, dsy, h0, n, esp, condut, & 
						neta, zeta, cx, cy, z, Ex_p, Ey_p, Ez_p, Hx_p, Hy_p, Hz_p)
	implicit none
	integer(4),intent(in)::n, idtfcd_cJ0, nJ0, idtfcd_cJ1, nJ1
	real(8),intent(in)::Tx,Ty,Iw,dsy,h0,esp(:),condut(1:n),cx,cy,z !,neta
	complex*16,intent(in)::zeta,neta
	complex*16,intent(out)::Ex_p,Ey_p,Ez_p,Hx_p,Hy_p,Hz_p

	real(8),parameter::pi=3.141592653589793238462643383279502884197d0
	integer(4)::i,j,k,camad,camadT,filtro,ident_fJ0,ident_fJ1
	real(8)::x,y,r,kr,k_r
	real(8),dimension(:),allocatable::h,krJ0,krJ1,w_J0,w_J1,prof

!	Para uso de loops:
	complex*16,dimension(:),allocatable::wvnb2,AMdwJ0,AMdwJ1,AMupJ0,AMupJ1,FEdwJ0,FEdwJ1,FEupJ0,FEupJ1
	complex*16,dimension(:,:),allocatable::uJ0,uJ1,AdmIntJ0,AdmIntJ1,ImpIntJ0,ImpIntJ1,uhJ0,uhJ1,tghJ0,tghJ1
	complex*16,dimension(:,:),allocatable::AdmApdwJ0,AdmApdwJ1,ImpApdwJ0,ImpApdwJ1,AdmApupJ0,AdmApupJ1,ImpApupJ0,ImpApupJ1
	complex*16,dimension(:,:),allocatable::RTEdwJ0,RTEdwJ1,RTMdwJ0,RTMdwJ1,RTEupJ0,RTEupJ1,RTMupJ0,RTMupJ1
	complex*16,dimension(:,:),allocatable::TMdwJ0,TMdwJ1,TEdwJ0,TEdwJ1,TMupJ0,TMupJ1,TEupJ0,TEupJ1

	complex*16,dimension(:),allocatable::Ktmdz_J0,Ktmdz_J1,Ktm_J0,Ktm_J1,Kte_J0,Kte_J1,Ktedz_J0,Ktedz_J1
	complex*16,dimension(:),allocatable::kernelExJ0,kernelExJ1,kernelEyJ0,kernelEyJ1,kernelEzJ1
	complex*16,dimension(:),allocatable::kernelHxJ0,kernelHxJ1,kernelHyJ0,kernelHyJ1,kernelHzJ1

	complex*16::kerEx_J1,kerEx_J0,kerEy_J1,kerEy_J0,kerEz_J1
	complex*16::kerHx_J1,kerHx_J0,kerHy_J1,kerHy_J0,kerHz_J1

	if (cx==Tx) then
	x=1.d0
	else
	x=cx-Tx
	end if
	y=cy-Ty
	r = dsqrt(x**2 + y**2)

	allocate(h(0:n),prof(-1:n))
	if (size(esp)==n) then
		h(0)=0.d0
		h(1:n)=esp
	else
		h(0)=0.d0
		h(1:n-1)=esp
		h(n)=1.d300
	end if
!	criando um novo vetor de profundidades que se adeque à qualquer situação patológica
	prof(-1)=-1.d300
	prof(0)=0.d0
	if (n > 1) then
	 prof(1) = h(1)
	 if (n > 2) then
	  do k=2,n-1
	   prof(k) = prof(k-1) + h(k)
	  end do
	 end if
	end if
	prof(n)=1.d300

!para descobrir em que camada está a observação
	if (z <= 0.d0) then
		camad=0
	else if (z > prof(n-1)) then
		camad=n
	else
	do i=n-1,1,-1
	if (z > prof(i-1)) then
		camad=i
		exit
	end if
	end do
	end if

!para descobrir em que camada está o transmissor
	if (h0 <= 0.d0) then
		camadT = 0
	else if (h0 > prof(n-1)) then
		camadT = n
	else
	do j=n-1,1,-1
	if (h0 > prof(j-1)) then
		camadT = j
		exit
	end if
	end do
	end if

!!	write(*,*)'Entre com o criador dos filtros J0: Rijo(0), Frayzer(1), Guptasarma(2), Kong(3) ou Key(4)'
!!	read(*,*)idtfcd_cJ0
	filtro=0	!esta variável direciona o uso de filtros J0 e J1 em vez de seno e cosseno
!!	call identfiltro(filtro,idtfcd_cJ0,ident_fJ0,nJ0)
!!	write(*,*)'Entre com o criador dos filtros J1: Rijo(0), Frayzer(1), Guptasarma(2), Kong(3) ou Key(4)'
!!	read(*,*)idtfcd_cJ1
!!	call identfiltro(filtro,idtfcd_cJ1,ident_fJ1,nJ1)

!	idtfcd_cJ0 = 4
	ident_fJ0 = 0
!	nJ0 = 101
!	idtfcd_cJ1 = 4
	ident_fJ1 = 1
!	nJ1 = 101

	allocate(KrJ0(nJ0),KrJ1(nJ1),w_J0(nJ0),w_J1(nJ1))

	call constfiltro(filtro,idtfcd_cJ0,ident_fJ0,nJ0,r,KrJ0,w_J0)
	call constfiltro(filtro,idtfcd_cJ1,ident_fJ1,nJ1,r,KrJ1,w_J1)

 	allocate(wvnb2(0:n),uJ0(nJ0,0:n),uJ1(nJ1,0:n),AdmIntJ0(nJ0,0:n),AdmIntJ1(nJ1,0:n),ImpIntJ0(nJ0,0:n),ImpIntJ1(nJ1,0:n))
 	allocate(uhJ0(nJ0,0:n),uhJ1(nJ1,0:n),tghJ0(nJ0,0:n),tghJ1(nJ1,0:n))
! 
 	do i=0,n
 		if (i==0) then
 		wvnb2(i)=-zeta*neta
 		uJ0(:,i)=cdsqrt(krJ0*krJ0-wvnb2(i))
 		uJ1(:,i)=cdsqrt(krJ1*krJ1-wvnb2(i))
 		AdmIntJ0(:,i)=uJ0(:,i)/zeta
 		AdmIntJ1(:,i)=uJ1(:,i)/zeta
 		ImpIntJ0(:,i)=uJ0(:,i)/neta
 		ImpIntJ1(:,i)=uJ1(:,i)/neta
 		uhJ0(:,i)=uJ0(:,i)*h(i)
 		uhJ1(:,i)=uJ1(:,i)*h(i)
 		tghJ0(:,i)=(1.d0-cdexp(-2.d0*uhJ0(:,i)))/(1.d0+cdexp(-2.d0*uhJ0(:,i)))
 		tghJ1(:,i)=(1.d0-cdexp(-2.d0*uhJ1(:,i)))/(1.d0+cdexp(-2.d0*uhJ1(:,i)))
 		else
 		wvnb2(i) = -zeta*condut(i)
 		uJ0(:,i)=cdsqrt(krJ0*krJ0-wvnb2(i))
 		uJ1(:,i)=cdsqrt(krJ1*krJ1-wvnb2(i))
 		AdmIntJ0(:,i)=uJ0(:,i)/zeta
 		AdmIntJ1(:,i)=uJ1(:,i)/zeta
 		ImpIntJ0(:,i)=uJ0(:,i)/condut(i)
 		ImpIntJ1(:,i)=uJ1(:,i)/condut(i)
 		uhJ0(:,i)=uJ0(:,i)*h(i)
 		uhJ1(:,i)=uJ1(:,i)*h(i)
 		tghJ0(:,i)=(1.d0-cdexp(-2.d0*uhJ0(:,i)))/(1.d0+cdexp(-2.d0*uhJ0(:,i)))
 		tghJ1(:,i)=(1.d0-cdexp(-2.d0*uhJ1(:,i)))/(1.d0+cdexp(-2.d0*uhJ1(:,i)))
 		end if
 	end do
 
 	allocate(AdmApdwJ0(nJ0,1:n),AdmApdwJ1(nJ1,1:n),ImpApdwJ0(nJ0,1:n),ImpApdwJ1(nJ1,1:n))
 	allocate(RTEdwJ0(nJ0,0:n),RTEdwJ1(nJ1,0:n),RTMdwJ0(nJ0,0:n),RTMdwJ1(nJ1,0:n))
 
 	do i=n,1,-1
 		if (i==n) then
 		AdmApdwJ0(:,i)=AdmIntJ0(:,i)
 		AdmApdwJ1(:,i)=AdmIntJ1(:,i)
 		ImpApdwJ0(:,i)=ImpIntJ0(:,i)
 		ImpApdwJ1(:,i)=ImpIntJ1(:,i)
 		RTEdwJ0(:,i)=(0.d0,0.d0)
 		RTEdwJ1(:,i)=(0.d0,0.d0)
 		RTMdwJ0(:,i)=(0.d0,0.d0)
 		RTMdwJ1(:,i)=(0.d0,0.d0)
 		else
 		AdmApdwJ0(:,i)=AdmIntJ0(:,i)*(AdmApdwJ0(:,i+1)+AdmIntJ0(:,i)* &
 			tghJ0(:,i))/(AdmIntJ0(:,i)+AdmApdwJ0(:,i+1)*tghJ0(:,i))
 		AdmApdwJ1(:,i)=AdmIntJ1(:,i)*(AdmApdwJ1(:,i+1)+AdmIntJ1(:,i)* &
 			tghJ1(:,i))/(AdmIntJ1(:,i)+AdmApdwJ1(:,i+1)*tghJ1(:,i))
 		ImpApdwJ0(:,i)=ImpIntJ0(:,i)*(ImpApdwJ0(:,i+1)+ImpIntJ0(:,i)* &
 			tghJ0(:,i))/(ImpIntJ0(:,i)+ImpApdwJ0(:,i+1)*tghJ0(:,i))
 		ImpApdwJ1(:,i)=ImpIntJ1(:,i)*(ImpApdwJ1(:,i+1)+ImpIntJ1(:,i)* &
 			tghJ1(:,i))/(ImpIntJ1(:,i)+ImpApdwJ1(:,i+1)*tghJ1(:,i))
 		RTEdwJ0(:,i)=(AdmIntJ0(:,i)-AdmApdwJ0(:,i+1))/(AdmIntJ0(:,i)+AdmApdwJ0(:,i+1))
 		RTEdwJ1(:,i)=(AdmIntJ1(:,i)-AdmApdwJ1(:,i+1))/(AdmIntJ1(:,i)+AdmApdwJ1(:,i+1))
 		RTMdwJ0(:,i)=(ImpIntJ0(:,i)-ImpApdwJ0(:,i+1))/(ImpIntJ0(:,i)+ImpApdwJ0(:,i+1))
 		RTMdwJ1(:,i)=(ImpIntJ1(:,i)-ImpApdwJ1(:,i+1))/(ImpIntJ1(:,i)+ImpApdwJ1(:,i+1))
 		end if	
 	end do
 		RTEdwJ0(:,0)=(AdmIntJ0(:,0)-AdmApdwJ0(:,1))/(AdmIntJ0(:,0)+AdmApdwJ0(:,1))
 		RTEdwJ1(:,0)=(AdmIntJ1(:,0)-AdmApdwJ1(:,1))/(AdmIntJ1(:,0)+AdmApdwJ1(:,1))
 		RTMdwJ0(:,0)=(ImpIntJ0(:,0)-ImpApdwJ0(:,1))/(ImpIntJ0(:,0)+ImpApdwJ0(:,1))
 		RTMdwJ1(:,0)=(ImpIntJ1(:,0)-ImpApdwJ1(:,1))/(ImpIntJ1(:,0)+ImpApdwJ1(:,1))
 
 	allocate(AdmApupJ0(nJ0,0:n-1),AdmApupJ1(nJ1,0:n-1),ImpApupJ0(nJ0,0:n-1),ImpApupJ1(nJ1,0:n-1))
 	allocate(RTEupJ0(nJ0,0:n),RTEupJ1(nJ1,0:n),RTMupJ0(nJ0,0:n),RTMupJ1(nJ1,0:n))
 
 	do i=0,n-1
 		if (i==0) then
 		AdmApupJ0(:,i)=AdmIntJ0(:,i)
 		AdmApupJ1(:,i)=AdmIntJ1(:,i)
 		ImpApupJ0(:,i)=ImpIntJ0(:,i)
 		ImpApupJ1(:,i)=ImpIntJ1(:,i)
 		RTEupJ0(:,i)=(0.d0,0.d0)
 		RTEupJ1(:,i)=(0.d0,0.d0)
 		RTMupJ0(:,i)=(0.d0,0.d0)
 		RTMupJ1(:,i)=(0.d0,0.d0)
 		else
 		AdmApupJ0(:,i)=AdmIntJ0(:,i)*(AdmApupJ0(:,i-1)+AdmIntJ0(:,i)* &
 			tghJ0(:,i))/(AdmIntJ0(:,i)+AdmApupJ0(:,i-1)*tghJ0(:,i))
 		AdmApupJ1(:,i)=AdmIntJ1(:,i)*(AdmApupJ1(:,i-1)+AdmIntJ1(:,i)* &
 			tghJ1(:,i))/(AdmIntJ1(:,i)+AdmApupJ1(:,i-1)*tghJ1(:,i))
 		ImpApupJ0(:,i)=ImpIntJ0(:,i)*(ImpApupJ0(:,i-1)+ImpIntJ0(:,i)* &
 			tghJ0(:,i))/(ImpIntJ0(:,i)+ImpApupJ0(:,i-1)*tghJ0(:,i))
 		ImpApupJ1(:,i)=ImpIntJ1(:,i)*(ImpApupJ1(:,i-1)+ImpIntJ1(:,i)* &
 			tghJ1(:,i))/(ImpIntJ1(:,i)+ImpApupJ1(:,i-1)*tghJ1(:,i))
 		RTEupJ0(:,i)=(AdmIntJ0(:,i)-AdmApupJ0(:,i-1))/(AdmIntJ0(:,i)+AdmApupJ0(:,i-1))
 		RTEupJ1(:,i)=(AdmIntJ1(:,i)-AdmApupJ1(:,i-1))/(AdmIntJ1(:,i)+AdmApupJ1(:,i-1))
 		RTMupJ0(:,i)=(ImpIntJ0(:,i)-ImpApupJ0(:,i-1))/(ImpIntJ0(:,i)+ImpApupJ0(:,i-1))
 		RTMupJ1(:,i)=(ImpIntJ1(:,i)-ImpApupJ1(:,i-1))/(ImpIntJ1(:,i)+ImpApupJ1(:,i-1))
 		end if
 	end do
 		RTEupJ0(:,n)=(AdmIntJ0(:,n)-AdmApupJ0(:,n-1))/(AdmIntJ0(:,n)+AdmApupJ0(:,n-1))
 		RTEupJ1(:,n)=(AdmIntJ1(:,n)-AdmApupJ1(:,n-1))/(AdmIntJ1(:,n)+AdmApupJ1(:,n-1))
 		RTMupJ0(:,n)=(ImpIntJ0(:,n)-ImpApupJ0(:,n-1))/(ImpIntJ0(:,n)+ImpApupJ0(:,n-1))
 		RTMupJ1(:,n)=(ImpIntJ1(:,n)-ImpApupJ1(:,n-1))/(ImpIntJ1(:,n)+ImpApupJ1(:,n-1))
 
 	allocate(AMdwJ0(nJ0),AMdwJ1(nJ1),AMupJ0(nJ0),AMupJ1(nJ1),FEdwJ0(nJ0),FEdwJ1(nJ1),FEupJ0(nJ0),FEupJ1(nJ1))

 	AMdwJ0=(cdexp(-uJ0(:,camadT)*(prof(camadT)-h0)) - RTMupJ0(:,camadT)* &
				cdexp(uJ0(:,camadT)*(prof(camadT-1)-h(camadT)-h0))) / &
				(1.d0-RTMupJ0(:,camadT)*RTMdwJ0(:,camadT)*cdexp(-2.d0*uhJ0(:,camadT)))
 	AMdwJ1=(cdexp(-uJ1(:,camadT)*(prof(camadT)-h0)) - RTMupJ1(:,camadT)* &
 				cdexp(uJ1(:,camadT)*(prof(camadT-1)-h(camadT)-h0))) / &
 				(1.d0-RTMupJ1(:,camadT)*RTMdwJ1(:,camadT)*cdexp(-2.d0*uhJ1(:,camadT)))
 	AMupJ0=(cdexp(uJ0(:,camadT)*(prof(camadT-1)-h0))-RTMdwJ0(:,camadT)* &
 				cdexp(-uJ0(:,camadT)*(prof(camadT)+h(camadT)-h0))) / &
 				(1.d0-RTMupJ0(:,camadT)*RTMdwJ0(:,camadT)*cdexp(-2.d0*uhJ0(:,camadT)))
 	AMupJ1=(cdexp(uJ1(:,camadT)*(prof(camadT-1)-h0))-RTMdwJ1(:,camadT)* &
 				cdexp(-uJ1(:,camadT)*(prof(camadT)+h(camadT)-h0))) / &
 				(1.d0-RTMupJ1(:,camadT)*RTMdwJ1(:,camadT)*cdexp(-2.d0*uhJ1(:,camadT)))
 	FEdwJ0=(cdexp(-uJ0(:,camadT)*(prof(camadT)-h0))+RTEupJ0(:,camadT)* &
 				cdexp(uJ0(:,camadT)*(prof(camadT-1)-h(camadT)-h0))) / &
 				(1.d0-RTEupJ0(:,camadT)*RTEdwJ0(:,camadT)*cdexp(-2.d0*uhJ0(:,camadT)))
 	FEdwJ1=(cdexp(-uJ1(:,camadT)*(prof(camadT)-h0))+RTEupJ1(:,camadT)* &
 				cdexp(uJ1(:,camadT)*(prof(camadT-1)-h(camadT)-h0))) / &
 				(1.d0-RTEupJ1(:,camadT)*RTEdwJ1(:,camadT)*cdexp(-2.d0*uhJ1(:,camadT)))
 	FEupJ0=(cdexp(uJ0(:,camadT)*(prof(camadT-1)-h0))+RTEdwJ0(:,camadT)* &
 				cdexp(-uJ0(:,camadT)*(prof(camadT)+h(camadT)-h0))) / &
 				(1.d0-RTEupJ0(:,camadT)*RTEdwJ0(:,camadT)*cdexp(-2.d0*uhJ0(:,camadT)))
 	FEupJ1=(cdexp(uJ1(:,camadT)*(prof(camadT-1)-h0))+RTEdwJ1(:,camadT)* &
 				cdexp(-uJ1(:,camadT)*(prof(camadT)+h(camadT)-h0))) / &
 				(1.d0-RTEupJ1(:,camadT)*RTEdwJ1(:,camadT)*cdexp(-2.d0*uhJ1(:,camadT)))

 	if (camad > camadT) then
 	allocate(TMdwJ0(nJ0,camadT:camad),TMdwJ1(nJ1,camadT:camad),TEdwJ0(nJ0,camadT:camad),TEdwJ1(nJ1,camadT:camad))
 		do j=camadT,camad
 			if (j == camadT) then
 			TMdwJ0(:,j)= - Iw * dsy / 2.d0
 			TMdwJ1(:,j)= - Iw * dsy / 2.d0
 			TEdwJ0(:,j)= - Iw * dsy / (2.d0 * AdmIntJ0(:,camadT))
 			TEdwJ1(:,j)= - Iw * dsy / (2.d0 * AdmIntJ1(:,camadT))
 			elseif (j == (camadT + 1) .and. j == n) then
 			TMdwJ0(:,j)=TMdwJ0(:,j-1)*(cdexp(-uJ0(:,camadT)*(prof(camadT)-h0)) - &
				RTMupJ0(:,camadT)*AMupJ0(:)*cdexp(-uhJ0(:,camadT)) + &
				RTMdwJ0(:,camadT)*AMdwJ0(:))
 			TMdwJ1(:,j)=TMdwJ1(:,j-1)*(cdexp(-uJ1(:,camadT)*(prof(camadT)-h0)) - &
				RTMupJ1(:,camadT)*AMupJ1(:)*cdexp(-uhJ1(:,camadT)) + &
				RTMdwJ1(:,camadT)*AMdwJ1(:))
 			TEdwJ0(:,j)=TEdwJ0(:,j-1)*(cdexp(-uJ0(:,camadT)*(prof(camadT)-h0)) + &
				RTEupJ0(:,camadT)*FEupJ0(:)*cdexp(-uhJ0(:,camadT)) + &
				RTEdwJ0(:,camadT)*FEdwJ0(:))
			TEdwJ1(:,j)=TEdwJ1(:,j-1)*(cdexp(-uJ1(:,camadT)*(prof(camadT)-h0)) + &
				RTEupJ1(:,camadT)*FEupJ1(:)*cdexp(-uhJ1(:,camadT)) + &
				RTEdwJ1(:,camadT)*FEdwJ1(:))
 			elseif (j == (camadT + 1) .and. j /= n) then
 			TMdwJ0(:,j)=TMdwJ0(:,j-1)*(cdexp(-uJ0(:,camadT)*(prof(camadT)-h0)) - &
				RTMupJ0(:,camadT)*AMupJ0(:)*cdexp(-uhJ0(:,camadT)) + &
				RTMdwJ0(:,camadT)*AMdwJ0(:)) / (1.d0 + RTMdwJ0(:,j)*cdexp(-2.d0*uhJ0(:,j)))
 			TMdwJ1(:,j)=TMdwJ1(:,j-1)*(cdexp(-uJ1(:,camadT)*(prof(camadT)-h0)) - &
				RTMupJ1(:,camadT)*AMupJ1(:)*cdexp(-uhJ1(:,camadT)) + &
				RTMdwJ1(:,camadT)*AMdwJ1(:)) / (1.d0 + RTMdwJ1(:,j)*cdexp(-2.d0*uhJ1(:,j)))
 			TEdwJ0(:,j)=TEdwJ0(:,j-1)*(cdexp(-uJ0(:,camadT)*(prof(camadT)-h0)) + &
				RTEupJ0(:,camadT)*FEupJ0(:)*cdexp(-uhJ0(:,camadT)) + &
				RTEdwJ0(:,camadT)*FEdwJ0(:)) / (1.d0 + RTEdwJ0(:,j)*cdexp(-2.d0*uhJ0(:,j)))
 			TEdwJ1(:,j)=TEdwJ1(:,j-1)*(cdexp(-uJ1(:,camadT)*(prof(camadT)-h0)) + &
				RTEupJ1(:,camadT)*FEupJ1(:)*cdexp(-uhJ1(:,camadT)) + &
				RTEdwJ1(:,camadT)*FEdwJ1(:)) / (1.d0 + RTEdwJ1(:,j)*cdexp(-2.d0*uhJ1(:,j)))
 			elseif (j /= n) then
 			TMdwJ0(:,j)=TMdwJ0(:,j-1)*(1.d0+RTMdwJ0(:,j-1))*cdexp(-uhJ0(:,j-1))/ &
 						(1.d0 + RTMdwJ0(:,j)*cdexp(-2.d0*uhJ0(:,j)))
 			TMdwJ1(:,j)=TMdwJ1(:,j-1)*(1.d0+RTMdwJ1(:,j-1))*cdexp(-uhJ1(:,j-1))/ &
 						(1.d0 + RTMdwJ1(:,j)*cdexp(-2.d0*uhJ1(:,j)))
 			TEdwJ0(:,j)=TEdwJ0(:,j-1)*(1.d0+RTEdwJ0(:,j-1))*cdexp(-uhJ0(:,j-1))/ &
 						(1.d0 + RTEdwJ0(:,j)*cdexp(-2.d0*uhJ0(:,j)))
 			TEdwJ1(:,j)=TEdwJ1(:,j-1)*(1.d0+RTEdwJ1(:,j-1))*cdexp(-uhJ1(:,j-1))/ &
 						(1.d0 + RTEdwJ1(:,j)*cdexp(-2.d0*uhJ1(:,j)))
 			elseif (j==n) then
 			TMdwJ0(:,j)=TMdwJ0(:,j-1)*(1.d0+RTMdwJ0(:,j-1))*cdexp(-uhJ0(:,j-1))
 			TMdwJ1(:,j)=TMdwJ1(:,j-1)*(1.d0+RTMdwJ1(:,j-1))*cdexp(-uhJ1(:,j-1))
 			TEdwJ0(:,j)=TEdwJ0(:,j-1)*(1.d0+RTEdwJ0(:,j-1))*cdexp(-uhJ0(:,j-1))
 			TEdwJ1(:,j)=TEdwJ1(:,j-1)*(1.d0+RTEdwJ1(:,j-1))*cdexp(-uhJ1(:,j-1))
 			end if
 		end do
 	elseif (camad < camadT) then
 	allocate(TMupJ0(nJ0,camad:camadT),TMupJ1(nJ1,camad:camadT),TEupJ0(nJ0,camad:camadT),TEupJ1(nJ1,camad:camadT))
 		do j=camadT,camad,-1
 			if (j == camadT) then
 			TMupJ0(:,j)= Iw * dsy / 2.d0
 			TMupJ1(:,j)= Iw * dsy / 2.d0
 			TEupJ0(:,j)= - Iw * dsy / (2.d0 * AdmIntJ0(:,camadT))
 			TEupJ1(:,j)= - Iw * dsy / (2.d0 * AdmIntJ1(:,camadT))
 			elseif (j == (camadT - 1) .and. j == 0) then
 			TMupJ0(:,j)=TMupJ0(:,j+1)*(cdexp(-uJ0(:,camadT)*h0) + &
				RTMupJ0(:,camadT)*AMupJ0(:) - RTMdwJ0(:,camadT)*AMdwJ0(:)* &
				cdexp(-uhJ0(:,camadT)))
 			TMupJ1(:,j)=TMupJ1(:,j+1)*(cdexp(-uJ1(:,camadT)*h0) + &
				RTMupJ1(:,camadT)*AMupJ1(:) - RTMdwJ1(:,camadT)*AMdwJ1(:)* &
				cdexp(-uhJ1(:,camadT)))
 			TEupJ0(:,j)=TEupJ0(:,j+1)*(cdexp(-uJ0(:,camadT)*h0) + &
				RTEupJ0(:,camadT)*FEupJ0(:) + RTEdwJ0(:,camadT)*FEdwJ0(:)* &
				cdexp(-uhJ0(:,camadT)))
 			TEupJ1(:,j)=TEupJ1(:,j+1)*(cdexp(-uJ1(:,camadT)*h0) + &
				RTEupJ1(:,camadT)*FEupJ1(:) + RTEdwJ1(:,camadT)*FEdwJ1(:)* &
				cdexp(-uhJ1(:,camadT)))
 			elseif (j == (camadT - 1) .and. j /= 0) then
 			TMupJ0(:,j)=TMupJ0(:,j+1)*(cdexp(uJ0(:,camadT)*(prof(camadT-1)-h0)) + &
				RTMupJ0(:,camadT)*AMupJ0(:) - RTMdwJ0(:,camadT)*AMdwJ0(:)* &
				cdexp(-uhJ0(:,camadT))) / (1.d0+RTMupJ0(:,j)*cdexp(-2.d0*uhJ0(:,j)))
 			TMupJ1(:,j)=TMupJ1(:,j+1)*(cdexp(uJ1(:,camadT)*(prof(camadT-1)-h0)) + &
				RTMupJ1(:,camadT)*AMupJ1(:) - RTMdwJ1(:,camadT)*AMdwJ1(:)* &
				cdexp(-uhJ1(:,camadT))) / (1.d0+RTMupJ1(:,j)*cdexp(-2.d0*uhJ1(:,j)))
 			TEupJ0(:,j)=TEupJ0(:,j+1)*(cdexp(uJ0(:,camadT)*(prof(camadT-1)-h0)) + &
				RTEupJ0(:,camadT)*FEupJ0(:) + RTEdwJ0(:,camadT)*FEdwJ0(:)* &
				cdexp(-uhJ0(:,camadT))) / (1.d0+RTEupJ0(:,j)*cdexp(-2.d0*uhJ0(:,j)))
 			TEupJ1(:,j)=TEupJ1(:,j+1)*(cdexp(uJ1(:,camadT)*(prof(camadT-1)-h0)) + &
				RTEupJ1(:,camadT)*FEupJ1(:) + RTEdwJ1(:,camadT)*FEdwJ1(:)* &
				cdexp(-uhJ1(:,camadT))) / (1.d0+RTEupJ1(:,j)*cdexp(-2.d0*uhJ1(:,j)))
 			elseif (j /= 0) then
 			TMupJ0(:,j)=TMupJ0(:,j+1)*(1.d0+RTMupJ0(:,j+1))*cdexp(-uhJ0(:,j+1)) / &
				(1.d0 + RTMupJ0(:,j)*cdexp(-2.d0*uhJ0(:,j)))
 			TMupJ1(:,j)=TMupJ1(:,j+1)*(1.d0+RTMupJ1(:,j+1))*cdexp(-uhJ1(:,j+1)) / &
				(1.d0 + RTMupJ1(:,j)*cdexp(-2.d0*uhJ1(:,j)))
 			TEupJ0(:,j)=TEupJ0(:,j+1)*(1.d0+RTEupJ0(:,j+1))*cdexp(-uhJ0(:,j+1)) / &
				(1.d0 + RTEupJ0(:,j)*cdexp(-2.d0*uhJ0(:,j)))
 			TEupJ1(:,j)=TEupJ1(:,j+1)*(1.d0+RTEupJ1(:,j+1))*cdexp(-uhJ1(:,j+1)) / &
				(1.d0 + RTEupJ1(:,j)*cdexp(-2.d0*uhJ1(:,j)))
 			elseif (j == 0) then
 			TMupJ0(:,j)=TMupJ0(:,1)*(1.d0+RTMupJ0(:,1))*cdexp(-uhJ0(:,1))
 			TMupJ1(:,j)=TMupJ1(:,1)*(1.d0+RTMupJ1(:,1))*cdexp(-uhJ1(:,1))
 			TEupJ0(:,j)=TEupJ0(:,1)*(1.d0+RTEupJ0(:,1))*cdexp(-uhJ0(:,1))
 			TEupJ1(:,j)=TEupJ1(:,1)*(1.d0+RTEupJ1(:,1))*cdexp(-uhJ1(:,1))
 			end if
 		end do
 	else
 	allocate(TMdwJ0(nJ0,camadT:camad),TMdwJ1(nJ1,camadT:camad),TEdwJ0(nJ0,camadT:camad),TEdwJ1(nJ1,camadT:camad))
 	allocate(TMupJ0(nJ0,camad:camadT),TMupJ1(nJ1,camad:camadT),TEupJ0(nJ0,camad:camadT),TEupJ1(nJ1,camad:camadT))
 			TMdwJ0(:,camad)= - Iw * dsy / 2.d0
 			TMdwJ1(:,camad)= - Iw * dsy / 2.d0
 			TEdwJ0(:,camad)= - Iw * dsy / (2.d0 * AdmIntJ0(:,camadT))
 			TEdwJ1(:,camad)= - Iw * dsy / (2.d0 * AdmIntJ1(:,camadT))
 			TMupJ0(:,camad)= Iw * dsy / 2.d0
 			TMupJ1(:,camad)= Iw * dsy / 2.d0
 			TEupJ0(:,camad)= - Iw * dsy / (2.d0 * AdmIntJ0(:,camadT))
 			TEupJ1(:,camad)= - Iw * dsy / (2.d0 * AdmIntJ1(:,camadT))
 	end if

 	allocate(Ktmdz_J0(nJ0),Ktmdz_J1(nJ1),Ktm_J0(nJ0),Ktm_J1(nJ1),Kte_J0(nJ0),Kte_J1(nJ1),Ktedz_J0(nJ0),Ktedz_J1(nJ1))
 	allocate(kernelExJ0(nJ0),kernelExJ1(nJ1),kernelEyJ0(nJ0),kernelEyJ1(nJ1),kernelEzJ1(nJ1))
 	allocate(kernelHxJ0(nJ0),kernelHxJ1(nJ1),kernelHyJ0(nJ0),kernelHyJ1(nJ1),kernelHzJ1(nJ1))
 	if (camad == 0 .and. camadT /= 0) then
 		Ktmdz_J0=(ImpIntJ0(:,0)*TMupJ0(:,0)*cdexp(uJ0(:,0)*z))*w_J0(:)
 		Ktmdz_J1=(ImpIntJ1(:,0)*TMupJ1(:,0)*cdexp(uJ1(:,0)*z))*w_J1(:)
 		Kte_J0=(TEupJ0(:,0)*cdexp(uJ0(:,0)*z))*w_J0(:)
 		Kte_J1=(TEupJ1(:,0)*cdexp(uJ1(:,0)*z))*w_J1(:)
 		Ktm_J0=(TMupJ0(:,0)*cdexp(uJ0(:,0)*z))*w_J0(:)
 		Ktm_J1=(TMupJ1(:,0)*cdexp(uJ1(:,0)*z))*w_J1(:)
 		Ktedz_J0=(AdmIntJ0(:,0)*TEupJ0(:,0)*cdexp(uJ0(:,0)*z))*w_J0(:)
 		Ktedz_J1=(AdmIntJ1(:,0)*TEupJ1(:,0)*cdexp(uJ1(:,0)*z))*w_J1(:)

 		kernelExJ1 = 2.d0*x*y/r**3 * (Ktmdz_J1 + Kte_J1)
 		kernelExJ0 = x*y/r**2 * (Ktmdz_J0 + Kte_J0) * krJ0
 		Ex_p = (sum(kernelExJ1) - sum(kernelExJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelEyJ1 = (2*y*y/r**3 - 1.d0/r) * Ktmdz_J1 - (2*x*x/r**3 - 1.d0/r) * Kte_J1
 		kernelEyJ0 = (y*y/r**2 * Ktmdz_J0 - x*x/r**2 * Kte_J0) * krJ0
 		Ey_p = (sum(kernelEyJ1) - sum(kernelEyJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelEzJ1 = y/(r*neta) * Ktm_J1 * krJ1 * krJ1
 		Ez_p = - sum(kernelEzJ1) / (2.d0*pi*r)				!este último r é decorrente do uso dos filtros

 		kernelHxJ1 = (2*y*y/r**3 - 1/r) * Ktm_J1 - (2*x*x/r**3 - 1/r) * Ktedz_J1
 		kernelHxJ0 = (y*y/r**2 * Ktm_J0 - x*x/r**2 * Ktedz_J0) * krJ0
 		Hx_p = (sum(kernelHxJ1) - sum(kernelHxJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros

 		kernelHyJ1 = 2*x*y/r**3 * (Ktm_J1 + Ktedz_J1)
 		kernelHyJ0 = x*y/r**2 * (Ktm_J0 + Ktedz_J0) * krJ0
 		Hy_p = - (sum(kernelHyJ1) - sum(kernelHyJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelHzJ1 = x/(r*zeta) * Kte_J1 * krJ1 * krJ1
 		Hz_p = sum(kernelHzJ1) / (2.d0*pi*r)				!este último r é decorrente do uso dos filtros
 
 	elseif (camad < camadT) then !camada k
 		ktmdz_J0=(ImpIntJ0(:,camad)*TMupJ0(:,camad)*(cdexp(uJ0(:,camad)*(z-prof(camad))) - &
 			RTMupJ0(:,camad)*cdexp(-uJ0(:,camad)*(z-prof(camad-1)+h(camad)))))*w_J0(:)
 		ktmdz_J1=(ImpIntJ1(:,camad)*TMupJ1(:,camad)*(cdexp(uJ1(:,camad)*(z-prof(camad))) - &
 			RTMupJ1(:,camad)*cdexp(-uJ1(:,camad)*(z-prof(camad-1)+h(camad)))))*w_J1(:)
 		Kte_J0=(TEupJ0(:,camad)*(cdexp(uJ0(:,camad)*(z-prof(camad))) + &
 			 RTEupJ0(:,camad)*cdexp(-uJ0(:,camad)*(z-prof(camad-1)+h(camad)))))*w_J0(:)
 		Kte_J1=(TEupJ1(:,camad)*(cdexp(uJ1(:,camad)*(z-prof(camad))) + &
 			 RTEupJ1(:,camad)*cdexp(-uJ1(:,camad)*(z-prof(camad-1)+h(camad)))))*w_J1(:)
 		ktm_J0=(TMupJ0(:,camad)*(cdexp(uJ0(:,camad)*(z-prof(camad))) + &
 			RTMupJ0(:,camad)*cdexp(-uJ0(:,camad)*(z-prof(camad-1)+h(camad)))))*w_J0(:)
 		ktm_J1=(TMupJ1(:,camad)*(cdexp(uJ1(:,camad)*(z-prof(camad))) + &
 			RTMupJ1(:,camad)*cdexp(-uJ1(:,camad)*(z-prof(camad-1)+h(camad)))))*w_J1(:)
 		Ktedz_J0=(AdmIntJ0(:,camad)*TEupJ0(:,camad)*(cdexp(uJ0(:,camad)*(z-prof(camad))) - &
 			RTEupJ0(:,camad)*cdexp(-uJ0(:,camad)*(z-prof(camad-1)+h(camad)))))*w_J0(:)
 		Ktedz_J1=(AdmIntJ1(:,camad)*TEupJ1(:,camad)*(cdexp(uJ1(:,camad)*(z-prof(camad))) - &
 			RTEupJ1(:,camad)*cdexp(-uJ1(:,camad)*(z-prof(camad-1)+h(camad)))))*w_J1(:)
 
 		kernelExJ1 = 2.d0*x*y/r**3 * (Ktmdz_J1 + Kte_J1)
 		kernelExJ0 = x*y/r**2 * (Ktmdz_J0 + Kte_J0) * krJ0
 		Ex_p = (sum(kernelExJ1) - sum(kernelExJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelEyJ1 = (2*y*y/r**3 - 1.d0/r) * Ktmdz_J1 - (2*x*x/r**3 - 1.d0/r) * Kte_J1
 		kernelEyJ0 = (y*y/r**2 * Ktmdz_J0 - x*x/r**2 * Kte_J0) * krJ0
 		Ey_p = (sum(kernelEyJ1) - sum(kernelEyJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelEzJ1 = y/(r*condut(camad)) * Ktm_J1 * krJ1 * krJ1
 		Ez_p = - sum(kernelEzJ1) / (2.d0*pi*r)				!este último r é decorrente do uso dos filtros
 
 		kernelHxJ1 = (2*y*y/r**3 - 1/r) * Ktm_J1 - (2*x*x/r**3 - 1/r) * Ktedz_J1
 		kernelHxJ0 = (y*y/r**2 * Ktm_J0 - x*x/r**2 * Ktedz_J0) * krJ0
 		Hx_p = (sum(kernelHxJ1) - sum(kernelHxJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelHyJ1 = 2*x*y/r**3 * (Ktm_J1 + Ktedz_J1)
 		kernelHyJ0 = x*y/r**2 * (Ktm_J0 + Ktedz_J0) * krJ0
 		Hy_p = - (sum(kernelHyJ1) - sum(kernelHyJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelHzJ1 = x/(r*zeta) * Kte_J1 * krJ1 * krJ1
 		Hz_p = sum(kernelHzJ1) / (2.d0*pi*r)				!este último r é decorrente do uso dos filtros
 
 	elseif (camad == camadT .and. z <= h0) then	!na mesma camada do transmissor mas acima dele
		Ktmdz_J0=(ImpIntJ0(:,camad)*TMupJ0(:,camad)*(cdexp(uJ0(:,camad)*(z-h0)) - &
 			RTMupJ0(:,camad)*AMupJ0(:)*cdexp(-uJ0(:,camad)*(z-prof(camad-1))) - &
 			RTMdwJ0(:,camad)*AMdwJ0(:)*cdexp(uJ0(:,camad)*(z-prof(camad)))))*w_J0(:)
 		Ktmdz_J1=(ImpIntJ1(:,camad)*TMupJ1(:,camad)*(cdexp(uJ1(:,camad)*(z-h0)) - &
 			RTMupJ1(:,camad)*AMupJ1(:)*cdexp(-uJ1(:,camad)*(z-prof(camad-1))) - &
 			RTMdwJ1(:,camad)*AMdwJ1(:)*cdexp(uJ1(:,camad)*(z-prof(camad)))))*w_J1(:)
 		Kte_J0=(TEupJ0(:,camad)*(cdexp(uJ0(:,camad)*(z-h0)) + &
 			RTEupJ0(:,camad)*FEupJ0(:)*cdexp(-uJ0(:,camad)*(z-prof(camad-1))) + &
 			RTEdwJ0(:,camad)*FEdwJ0(:)*cdexp(uJ0(:,camad)*(z-prof(camad)))))*w_J0(:)
 		Kte_J1=(TEupJ1(:,camad)*(cdexp(uJ1(:,camad)*(z-h0)) + &
 			RTEupJ1(:,camad)*FEupJ1(:)*cdexp(-uJ1(:,camad)*(z-prof(camad-1))) + &
 			RTEdwJ1(:,camad)*FEdwJ1(:)*cdexp(uJ1(:,camad)*(z-prof(camad)))))*w_J1(:)
 		Ktm_J0=(TMupJ0(:,camad)*(cdexp(uJ0(:,camad)*(z-h0)) + &
 			RTMupJ0(:,camad)*AMupJ0(:)*cdexp(-uJ0(:,camad)*(z-prof(camad-1))) - &
 			RTMdwJ0(:,camad)*AMdwJ0(:)*cdexp(uJ0(:,camad)*(z-prof(camad)))))*w_J0(:)
 		Ktm_J1=(TMupJ1(:,camad)*(cdexp(uJ1(:,camad)*(z-h0)) + &
 			RTMupJ1(:,camad)*AMupJ1(:)*cdexp(-uJ1(:,camad)*(z-prof(camad-1))) - &
 			RTMdwJ1(:,camad)*AMdwJ1(:)*cdexp(uJ1(:,camad)*(z-prof(camad)))))*w_J1(:)
 		Ktedz_J0=(AdmIntJ0(:,camad)*TEupJ0(:,camad)*(cdexp(uJ0(:,camad)*(z-h0)) - &
 			RTEupJ0(:,camad)*FEupJ0(:)*cdexp(-uJ0(:,camad)*(z-prof(camad-1))) + &
 			RTEdwJ0(:,camad)*FEdwJ0(:)*cdexp(uJ0(:,camad)*(z-prof(camad)))))*w_J0(:)
 		Ktedz_J1=(AdmIntJ1(:,camad)*TEupJ1(:,camad)*(cdexp(uJ1(:,camad)*(z-h0)) - &
 			RTEupJ1(:,camad)*FEupJ1(:)*cdexp(-uJ1(:,camad)*(z-prof(camad-1))) + &
 			RTEdwJ1(:,camad)*FEdwJ1(:)*cdexp(uJ1(:,camad)*(z-prof(camad)))))*w_J1(:)

 		kernelExJ1 = 2.d0*x*y/r**3 * (Ktmdz_J1 + Kte_J1)
 		kernelExJ0 = x*y/r**2 * (Ktmdz_J0 + Kte_J0) * krJ0
 		Ex_p = (sum(kernelExJ1) - sum(kernelExJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelEyJ1 = (2*y*y/r**3 - 1.d0/r) * Ktmdz_J1 - (2*x*x/r**3 - 1.d0/r) * Kte_J1
 		kernelEyJ0 = (y*y/r**2 * Ktmdz_J0 - x*x/r**2 * Kte_J0) * krJ0
 		Ey_p = (sum(kernelEyJ1) - sum(kernelEyJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
		if (camad /= 0) then
 		kernelEzJ1 = y/(r*condut(camad)) * Ktm_J1 * krJ1 * krJ1
		else
 		kernelEzJ1 = y/(r*neta) * Ktm_J1 * krJ1 * krJ1
		end if
 		Ez_p = - sum(kernelEzJ1) / (2.d0*pi*r)						!este último r é decorrente do uso dos filtros
 
 		kernelHxJ1 = (2*y*y/r**3 - 1/r) * Ktm_J1 - (2*x*x/r**3 - 1/r) * Ktedz_J1
 		kernelHxJ0 = (y*y/r**2 * Ktm_J0 - x*x/r**2 * Ktedz_J0) * krJ0
 		Hx_p = (sum(kernelHxJ1) - sum(kernelHxJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelHyJ1 = 2*x*y/r**3 * (Ktm_J1 + Ktedz_J1)
 		kernelHyJ0 = x*y/r**2 * (Ktm_J0 + Ktedz_J0) * krJ0
 		Hy_p = - (sum(kernelHyJ1) - sum(kernelHyJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelHzJ1 = x/(r*zeta) * Kte_J1 * krJ1 * krJ1
 		Hz_p = sum(kernelHzJ1) / (2.d0*pi*r)				!este último r é decorrente do uso dos filtros
 
 	elseif (camad == camadT .and. z > h0) then	!na mesma camada do transmissor mas abaixo dele
 		Ktmdz_J0=(ImpIntJ0(:,camad)*TMdwJ0(:,camad)*(cdexp(-uJ0(:,camad)*(z-h0)) - &
 			RTMupJ0(:,camad)*AMupJ0(:)*cdexp(-uJ0(:,camad)*(z-prof(camad-1))) - &
 			RTMdwJ0(:,camad)*AMdwJ0(:)*cdexp(uJ0(:,camad)*(z-prof(camad)))))*w_J0(:)
 		Ktmdz_J1=(ImpIntJ1(:,camad)*TMdwJ1(:,camad)*(cdexp(-uJ1(:,camad)*(z-h0)) - &
 			RTMupJ1(:,camad)*AMupJ1(:)*cdexp(-uJ1(:,camad)*(z-prof(camad-1))) - &
 			RTMdwJ1(:,camad)*AMdwJ1(:)*cdexp(uJ1(:,camad)*(z-prof(camad)))))*w_J1(:)
 		Kte_J0=(TEdwJ0(:,camad)*(cdexp(-uJ0(:,camad)*(z-h0)) + &
 			RTEupJ0(:,camad)*FEupJ0(:)*cdexp(-uJ0(:,camad)*(z-prof(camad-1))) + &
 			RTEdwJ0(:,camad)*FEdwJ0(:)*cdexp(uJ0(:,camad)*(z-prof(camad)))))*w_J0(:)
 		Kte_J1=(TEdwJ1(:,camad)*(cdexp(-uJ1(:,camad)*(z-h0)) + &
 			RTEupJ1(:,camad)*FEupJ1(:)*cdexp(-uJ1(:,camad)*(z-prof(camad-1))) + &
 			RTEdwJ1(:,camad)*FEdwJ1(:)*cdexp(uJ1(:,camad)*(z-prof(camad)))))*w_J1(:)
		Ktm_J0=(TMdwJ0(:,camad)*(cdexp(-uJ0(:,camad)*(z-h0)) - &
 			RTMupJ0(:,camad)*AMupJ0(:)*cdexp(-uJ0(:,camad)*(z-prof(camad-1))) + &
 			RTMdwJ0(:,camad)*AMdwJ0(:)*cdexp(uJ0(:,camad)*(z-prof(camad)))))*w_J0(:)
 		Ktm_J1=(TMdwJ1(:,camad)*(cdexp(-uJ1(:,camad)*(z-h0)) - &
 			RTMupJ1(:,camad)*AMupJ1(:)*cdexp(-uJ1(:,camad)*(z-prof(camad-1))) + &
 			RTMdwJ1(:,camad)*AMdwJ1(:)*cdexp(uJ1(:,camad)*(z-prof(camad)))))*w_J1(:)
 		Ktedz_J0=(AdmIntJ0(:,camad)*TEdwJ0(:,camad)*(cdexp(-uJ0(:,camad)*(z-h0)) + &
 			RTEupJ0(:,camad)*FEupJ0(:)*cdexp(-uJ0(:,camad)*(z-prof(camad-1))) - &
 			RTEdwJ0(:,camad)*FEdwJ0(:)*cdexp(uJ0(:,camad)*(z-prof(camad)))))*w_J0(:)
 		Ktedz_J1=(AdmIntJ1(:,camad)*TEdwJ1(:,camad)*(cdexp(-uJ1(:,camad)*(z-h0)) + &
 			RTEupJ1(:,camad)*FEupJ1(:)*cdexp(-uJ1(:,camad)*(z-prof(camad-1))) - &
 			RTEdwJ1(:,camad)*FEdwJ1(:)*cdexp(uJ1(:,camad)*(z-prof(camad)))))*w_J1(:)
 
 		kernelExJ1 = 2.d0*x*y/r**3 * (-Ktmdz_J1 + Kte_J1)
 		kernelExJ0 = x*y/r**2 * (-Ktmdz_J0 + Kte_J0) * krJ0
 		Ex_p = (sum(kernelExJ1) - sum(kernelExJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelEyJ1 = (1.d0/r - 2*y*y/r**3) * Ktmdz_J1 + (1.d0/r - 2*x*x/r**3) * Kte_J1
 		kernelEyJ0 = (y*y/r**2 * Ktmdz_J0 + x*x/r**2 * Kte_J0) * krJ0
 		Ey_p = (sum(kernelEyJ1) + sum(kernelEyJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
		if (camad /= 0) then
 		kernelEzJ1 = y/(r*condut(camad)) * Ktm_J1 * krJ1 * krJ1
		else
 		kernelEzJ1 = y/(r*neta) * Ktm_J1 * krJ1 * krJ1
		end if
 		Ez_p = - sum(kernelEzJ1) / (2.d0*pi*r)						!este último r é decorrente do uso dos filtros
 
 		kernelHxJ1 = (2*y*y/r**3 - 1/r) * Ktm_J1 + (2*x*x/r**3 - 1/r) * Ktedz_J1
 		kernelHxJ0 = (y*y/r**2 * Ktm_J0 + x*x/r**2 * Ktedz_J0) * krJ0
 		Hx_p = (sum(kernelHxJ1) - sum(kernelHxJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelHyJ1 = 2*x*y/r**3 * (Ktm_J1 - Ktedz_J1)
 		kernelHyJ0 = x*y/r**2 * (Ktm_J0 - Ktedz_J0) * krJ0
 		Hy_p = - (sum(kernelHyJ1) - sum(kernelHyJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelHzJ1 = x/(r*zeta) * Kte_J1 * krJ1 * krJ1
 		Hz_p = sum(kernelHzJ1) / (2.d0*pi*r)				!este último r é decorrente do uso dos filtros
 
 	elseif (camad > camadT .and. camad /= n) then !camada j
 		Ktmdz_J0=(ImpIntJ0(:,camad)*TMdwJ0(:,camad)*(cdexp(-uJ0(:,camad)*(z-prof(camad-1))) - &
 			RTMdwJ0(:,camad)*cdexp(uJ0(:,camad)*(z-prof(camad)-h(camad)))))*w_J0(:)
 		Ktmdz_J1=(ImpIntJ1(:,camad)*TMdwJ1(:,camad)*(cdexp(-uJ1(:,camad)*(z-prof(camad-1))) - &
 			RTMdwJ1(:,camad)*cdexp(uJ1(:,camad)*(z-prof(camad)-h(camad)))))*w_J1(:)
 		Kte_J0=(TEdwJ0(:,camad)*(cdexp(-uJ0(:,camad)*(z-prof(camad-1))) + &
 			RTEdwJ0(:,camad)*cdexp(uJ0(:,camad)*(z-prof(camad)-h(camad)))))*w_J0(:)
 		Kte_J1=(TEdwJ1(:,camad)*(cdexp(-uJ1(:,camad)*(z-prof(camad-1))) + &
 			RTEdwJ1(:,camad)*cdexp(uJ1(:,camad)*(z-prof(camad)-h(camad)))))*w_J1(:)
 		Ktm_J0=(TMdwJ0(:,camad)*(cdexp(-uJ0(:,camad)*(z-prof(camad-1))) + &
 			RTMdwJ0(:,camad)*cdexp(uJ0(:,camad)*(z-prof(camad)-h(camad)))))*w_J0(:)
 		Ktm_J1=(TMdwJ1(:,camad)*(cdexp(-uJ1(:,camad)*(z-prof(camad-1))) + &
 			RTMdwJ1(:,camad)*cdexp(uJ1(:,camad)*(z-prof(camad)-h(camad)))))*w_J1(:)
 		Ktedz_J0=(AdmIntJ0(:,camad)*TEdwJ0(:,camad)*(cdexp(-uJ0(:,camad)*(z-prof(camad-1))) - &
 			RTEdwJ0(:,camad)*cdexp(uJ0(:,camad)*(z-prof(camad)-h(camad)))))*w_J0(:)
 		Ktedz_J1=(AdmIntJ1(:,camad)*TEdwJ1(:,camad)*(cdexp(-uJ1(:,camad)*(z-prof(camad-1))) - &
 			RTEdwJ1(:,camad)*cdexp(uJ1(:,camad)*(z-prof(camad)-h(camad)))))*w_J1(:)
 
 		kernelExJ1 = 2.d0*x*y/r**3 * (-Ktmdz_J1 + Kte_J1)
 		kernelExJ0 = x*y/r**2 * (-Ktmdz_J0 + Kte_J0) * krJ0
 		Ex_p = (sum(kernelExJ1) - sum(kernelExJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelEyJ1 = (1.d0/r - 2*y*y/r**3) * Ktmdz_J1 + (1.d0/r - 2*x*x/r**3) * Kte_J1
 		kernelEyJ0 = (y*y/r**2 * Ktmdz_J0 + x*x/r**2 * Kte_J0) * krJ0
 		Ey_p = (sum(kernelEyJ1) + sum(kernelEyJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelEzJ1 = y/(r*condut(camad)) * Ktm_J1 * krJ1 * krJ1
 		Ez_p = - sum(kernelEzJ1) / (2.d0*pi*r)				!este último r é decorrente do uso dos filtros
 
 		kernelHxJ1 = (2*y*y/r**3 - 1/r) * Ktm_J1 + (2*x*x/r**3 - 1/r) * Ktedz_J1
 		kernelHxJ0 = (y*y/r**2 * Ktm_J0 + x*x/r**2 * Ktedz_J0) * krJ0
 		Hx_p = (sum(kernelHxJ1) - sum(kernelHxJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelHyJ1 = 2*x*y/r**3 * (Ktm_J1 - Ktedz_J1)
 		kernelHyJ0 = x*y/r**2 * (Ktm_J0 - Ktedz_J0) * krJ0
 		Hy_p = - (sum(kernelHyJ1) - sum(kernelHyJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelHzJ1 = x/(r*zeta) * Kte_J1 * krJ1 * krJ1
 		Hz_p = sum(kernelHzJ1) / (2.d0*pi*r)				!este último r é decorrente do uso dos filtros

 	else	!camada n
 		Ktmdz_J0=(ImpIntJ0(:,n)*TMdwJ0(:,n)*cdexp(-uJ0(:,n)*(z-prof(n-1))))*w_J0(:)
 		Ktmdz_J1=(ImpIntJ1(:,n)*TMdwJ1(:,n)*cdexp(-uJ1(:,n)*(z-prof(n-1))))*w_J1(:)
 		Kte_J0=(TEdwJ0(:,n)*cdexp(-uJ0(:,n)*(z-prof(n-1))))*w_J0(:)
 		Kte_J1=(TEdwJ1(:,n)*cdexp(-uJ1(:,n)*(z-prof(n-1))))*w_J1(:)
 		Ktm_J0=(TMdwJ0(:,n)*cdexp(-uJ0(:,n)*(z-prof(n-1))))*w_J0(:)
 		Ktm_J1=(TMdwJ1(:,n)*cdexp(-uJ1(:,n)*(z-prof(n-1))))*w_J1(:)
 		Ktedz_J0=(AdmIntJ0(:,n)*TEdwJ0(:,n)*cdexp(-uJ0(:,n)*(z-prof(n-1))))*w_J0(:)
 		Ktedz_J1=(AdmIntJ1(:,n)*TEdwJ1(:,n)*cdexp(-uJ1(:,n)*(z-prof(n-1))))*w_J1(:)
 
 		kernelExJ1 = 2.d0*x*y/r**3 * (-Ktmdz_J1 + Kte_J1)
 		kernelExJ0 = x*y/r**2 * (-Ktmdz_J0 + Kte_J0) * krJ0
 		Ex_p = (sum(kernelExJ1) - sum(kernelExJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelEyJ1 = (1.d0/r - 2*y*y/r**3) * Ktmdz_J1 + (1.d0/r - 2*x*x/r**3) * Kte_J1
 		kernelEyJ0 = (y*y/r**2 * Ktmdz_J0 + x*x/r**2 * Kte_J0) * krJ0
 		Ey_p = (sum(kernelEyJ1) + sum(kernelEyJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelEzJ1 = y/(r*condut(camad)) * Ktm_J1 * krJ1 * krJ1
 		Ez_p = - sum(kernelEzJ1) / (2.d0*pi*r)				!este último r é decorrente do uso dos filtros
 
 		kernelHxJ1 = (2*y*y/r**3 - 1/r) * Ktm_J1 + (2*x*x/r**3 - 1/r) * Ktedz_J1
 		kernelHxJ0 = (y*y/r**2 * Ktm_J0 + x*x/r**2 * Ktedz_J0) * krJ0
 		Hx_p = (sum(kernelHxJ1) - sum(kernelHxJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelHyJ1 = 2*x*y/r**3 * (Ktm_J1 - Ktedz_J1)
 		kernelHyJ0 = x*y/r**2 * (Ktm_J0 - Ktedz_J0) * krJ0
 		Hy_p = - (sum(kernelHyJ1) - sum(kernelHyJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelHzJ1 = x/(r*zeta) * Kte_J1 * krJ1 * krJ1
 		Hz_p = sum(kernelHzJ1) / (2.d0*pi*r)				!este último r é decorrente do uso dos filtros
 
 	end if
 
 	deallocate(h,KrJ0,KrJ1,w_J0,w_J1)
 	deallocate(wvnb2,uJ0,uJ1,AdmIntJ0,AdmIntJ1,ImpIntJ0,ImpIntJ1,uhJ0,uhJ1,tghJ0,tghJ1)
 	deallocate(AdmApdwJ0,AdmApdwJ1,ImpApdwJ0,ImpApdwJ1,RTEdwJ0,RTEdwJ1,RTMdwJ0,RTMdwJ1)
 	deallocate(AdmApupJ0,AdmApupJ1,ImpApupJ0,ImpApupJ1,RTEupJ0,RTEupJ1,RTMupJ0,RTMupJ1)
 	deallocate(AMdwJ0,AMdwJ1,AMupJ0,AMupJ1,FEdwJ0,FEdwJ1,FEupJ0,FEupJ1)
 	deallocate(Ktmdz_J0,Ktmdz_J1,Ktm_J0,Ktm_J1,Kte_J0,Kte_J1,Ktedz_J0,Ktedz_J1)
 	deallocate(kernelExJ0,kernelExJ1,kernelEyJ0,kernelEyJ1,kernelEzJ1)
 	deallocate(kernelHxJ0,kernelHxJ1,kernelHyJ0,kernelHyJ1,kernelHzJ1)
	end subroutine dehy_xyz_loops
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

	subroutine dev_xyz_loops( idtfcd_cJ0, nJ0, idtfcd_cJ1, nJ1, Tx, Ty, Iw, dsz, h0, n, esp, condut, &
						 neta, zeta, cx, cy, z, Ex_p, Ey_p, Ez_p, Hx_p, Hy_p, Hz_p)
	implicit none
	integer(4),intent(in)::n, idtfcd_cJ0, nJ0, idtfcd_cJ1, nJ1
	real(8),intent(in)::Tx,Ty,Iw,dsz,h0,esp(:),condut(1:n),cx,cy,z !,neta
	complex*16,intent(in)::zeta,neta
	complex*16,intent(out)::Ex_p,Ey_p,Ez_p,Hx_p,Hy_p,Hz_p

	real(8),parameter::pi=3.141592653589793238462643383279502884197d0
	integer(4)::i,j,k,camad,camadT,filtro,ident_fJ0,ident_fJ1
	real(8)::x,y,r,kr,k_r
	real(8),dimension(:),allocatable::h,krJ0,krJ1,w_J0,w_J1,prof

!	Para uso de loops:
	complex*16,dimension(:),allocatable::wvnb2,AMdwJ0,AMdwJ1,AMupJ0,AMupJ1
	complex*16,dimension(:,:),allocatable::uJ0,uJ1,ImpIntJ0,ImpIntJ1,uhJ0,uhJ1,tghJ0,tghJ1
	complex*16,dimension(:,:),allocatable::ImpApdwJ0,ImpApdwJ1,ImpApupJ0,ImpApupJ1
	complex*16,dimension(:,:),allocatable::RTMdwJ0,RTMdwJ1,RTMupJ0,RTMupJ1
	complex*16,dimension(:,:),allocatable::TMdwJ0,TMdwJ1,TMupJ0,TMupJ1

	complex*16,dimension(:),allocatable::Ktmdz_J1,Ktm_J0,Ktm_J1
	complex*16,dimension(:),allocatable::kernelExJ1,kernelEyJ1,kernelEzJ0
	complex*16,dimension(:),allocatable::kernelHxJ1,kernelHyJ1

	Hz_p = (0.d0,0.d0)

	if (cx==Tx) then
	x=1.d0
	else
	x=cx-Tx
	end if
	y=cy-Ty
	r = dsqrt(x**2 + y**2)

	allocate(h(0:n),prof(-1:n))
	if (size(esp)==n) then
		h(0)=0.d0
		h(1:n)=esp
	else
		h(0)=0.d0
		h(1:n-1)=esp
		h(n)=1.d300
	end if
!	criando um novo vetor de profundidades que se adeque à qualquer situação patológica
	prof(-1)=-1.d300
	prof(0)=0.d0
	if (n > 1) then
	 prof(1) = h(1)
	 if (n > 2) then
	  do k=2,n-1
	   prof(k) = prof(k-1) + h(k)
	  end do
	 end if
	end if
	prof(n)=1.d300

!para descobrir em que camada está a observação
	if (z <= 0.d0) then
		camad=0
	else if (z > prof(n-1)) then
		camad=n
	else
	do i=n-1,1,-1
	if (z > prof(i-1)) then
		camad=i
		exit
	end if
	end do
	end if

!para descobrir em que camada está o transmissor
	if (h0 <= 0.d0) then
		camadT = 0
	else if (h0 > prof(n-1)) then
		camadT = n
	else
	do j=n-1,1,-1
	if (h0 > prof(j-1)) then
		camadT = j
		exit
	end if
	end do
	end if

!!	write(*,*)'Entre com o criador dos filtros J0: Rijo(0), Frayzer(1), Guptasarma(2), Kong(3) ou Key(4)'
!!	read(*,*)idtfcd_cJ0
	filtro=0	!esta variável direciona o uso de filtros J0 e J1 em vez de seno e cosseno
!!	call identfiltro(filtro,idtfcd_cJ0,ident_fJ0,nJ0)
!!	write(*,*)'Entre com o criador dos filtros J1: Rijo(0), Frayzer(1), Guptasarma(2), Kong(3) ou Key(4)'
!!	read(*,*)idtfcd_cJ1
!!	call identfiltro(filtro,idtfcd_cJ1,ident_fJ1,nJ1)

!	idtfcd_cJ0 = 3
	ident_fJ0 = 0
!	nJ0 = 241
!	idtfcd_cJ1 = 3
	ident_fJ1 = 1
!	nJ1 = 241

	allocate(KrJ0(nJ0),KrJ1(nJ1),w_J0(nJ0),w_J1(nJ1))

	call constfiltro(filtro,idtfcd_cJ0,ident_fJ0,nJ0,r,KrJ0,w_J0)
	call constfiltro(filtro,idtfcd_cJ1,ident_fJ1,nJ1,r,KrJ1,w_J1)

 	allocate(wvnb2(0:n),uJ0(nJ0,0:n),uJ1(nJ1,0:n),ImpIntJ0(nJ0,0:n),ImpIntJ1(nJ1,0:n))
 	allocate(uhJ0(nJ0,0:n),uhJ1(nJ1,0:n),tghJ0(nJ0,0:n),tghJ1(nJ1,0:n))
! 
 	do i=0,n
 		if (i==0) then
 		wvnb2(i)=-zeta*neta
 		uJ0(:,i)=cdsqrt(krJ0*krJ0-wvnb2(i))
 		uJ1(:,i)=cdsqrt(krJ1*krJ1-wvnb2(i))
 		ImpIntJ0(:,i)=uJ0(:,i)/neta
 		ImpIntJ1(:,i)=uJ1(:,i)/neta
 		uhJ0(:,i)=uJ0(:,i)*h(i)
 		uhJ1(:,i)=uJ1(:,i)*h(i)
 		tghJ0(:,i)=(1.d0-cdexp(-2.d0*uhJ0(:,i)))/(1.d0+cdexp(-2.d0*uhJ0(:,i)))
 		tghJ1(:,i)=(1.d0-cdexp(-2.d0*uhJ1(:,i)))/(1.d0+cdexp(-2.d0*uhJ1(:,i)))
 		else
 		wvnb2(i) = -zeta*condut(i)
 		uJ0(:,i)=cdsqrt(krJ0*krJ0-wvnb2(i))
 		uJ1(:,i)=cdsqrt(krJ1*krJ1-wvnb2(i))
 		ImpIntJ0(:,i)=uJ0(:,i)/condut(i)
 		ImpIntJ1(:,i)=uJ1(:,i)/condut(i)
 		uhJ0(:,i)=uJ0(:,i)*h(i)
 		uhJ1(:,i)=uJ1(:,i)*h(i)
 		tghJ0(:,i)=(1.d0-cdexp(-2.d0*uhJ0(:,i)))/(1.d0+cdexp(-2.d0*uhJ0(:,i)))
 		tghJ1(:,i)=(1.d0-cdexp(-2.d0*uhJ1(:,i)))/(1.d0+cdexp(-2.d0*uhJ1(:,i)))
 		end if
 	end do
 
 	allocate(ImpApdwJ0(nJ0,1:n),ImpApdwJ1(nJ1,1:n),RTMdwJ0(nJ0,0:n),RTMdwJ1(nJ1,0:n))
 
 	do i=n,1,-1
 		if (i==n) then
 		ImpApdwJ0(:,i)=ImpIntJ0(:,i)
 		ImpApdwJ1(:,i)=ImpIntJ1(:,i)
 		RTMdwJ0(:,i)=(0.d0,0.d0)
 		RTMdwJ1(:,i)=(0.d0,0.d0)
 		else
 		ImpApdwJ0(:,i)=ImpIntJ0(:,i)*(ImpApdwJ0(:,i+1)+ImpIntJ0(:,i)* &
 			tghJ0(:,i))/(ImpIntJ0(:,i)+ImpApdwJ0(:,i+1)*tghJ0(:,i))
 		ImpApdwJ1(:,i)=ImpIntJ1(:,i)*(ImpApdwJ1(:,i+1)+ImpIntJ1(:,i)* &
 			tghJ1(:,i))/(ImpIntJ1(:,i)+ImpApdwJ1(:,i+1)*tghJ1(:,i))
 		RTMdwJ0(:,i)=(ImpIntJ0(:,i)-ImpApdwJ0(:,i+1))/(ImpIntJ0(:,i)+ImpApdwJ0(:,i+1))
 		RTMdwJ1(:,i)=(ImpIntJ1(:,i)-ImpApdwJ1(:,i+1))/(ImpIntJ1(:,i)+ImpApdwJ1(:,i+1))
 		end if	
 	end do
 		RTMdwJ0(:,0)=(ImpIntJ0(:,0)-ImpApdwJ0(:,1))/(ImpIntJ0(:,0)+ImpApdwJ0(:,1))
 		RTMdwJ1(:,0)=(ImpIntJ1(:,0)-ImpApdwJ1(:,1))/(ImpIntJ1(:,0)+ImpApdwJ1(:,1))
 
 	allocate(ImpApupJ0(nJ0,0:n-1),ImpApupJ1(nJ1,0:n-1),RTMupJ0(nJ0,0:n),RTMupJ1(nJ1,0:n))
 
 	do i=0,n-1
 		if (i==0) then
 		ImpApupJ0(:,i)=ImpIntJ0(:,i)
 		ImpApupJ1(:,i)=ImpIntJ1(:,i)
 		RTMupJ0(:,i)=(0.d0,0.d0)
 		RTMupJ1(:,i)=(0.d0,0.d0)
 		else
 		ImpApupJ0(:,i)=ImpIntJ0(:,i)*(ImpApupJ0(:,i-1)+ImpIntJ0(:,i)* &
 			tghJ0(:,i))/(ImpIntJ0(:,i)+ImpApupJ0(:,i-1)*tghJ0(:,i))
 		ImpApupJ1(:,i)=ImpIntJ1(:,i)*(ImpApupJ1(:,i-1)+ImpIntJ1(:,i)* &
 			tghJ1(:,i))/(ImpIntJ1(:,i)+ImpApupJ1(:,i-1)*tghJ1(:,i))
 		RTMupJ0(:,i)=(ImpIntJ0(:,i)-ImpApupJ0(:,i-1))/(ImpIntJ0(:,i)+ImpApupJ0(:,i-1))
 		RTMupJ1(:,i)=(ImpIntJ1(:,i)-ImpApupJ1(:,i-1))/(ImpIntJ1(:,i)+ImpApupJ1(:,i-1))
 		end if
 	end do
 		RTMupJ0(:,n)=(ImpIntJ0(:,n)-ImpApupJ0(:,n-1))/(ImpIntJ0(:,n)+ImpApupJ0(:,n-1))
 		RTMupJ1(:,n)=(ImpIntJ1(:,n)-ImpApupJ1(:,n-1))/(ImpIntJ1(:,n)+ImpApupJ1(:,n-1))
 
 	allocate(AMdwJ0(nJ0),AMdwJ1(nJ1),AMupJ0(nJ0),AMupJ1(nJ1))

 	AMdwJ0=(cdexp(-uJ0(:,camadT)*(prof(camadT)-h0)) + RTMupJ0(:,camadT)* &
				cdexp(uJ0(:,camadT)*(prof(camadT-1)-h(camadT)-h0))) / &
				(1.d0-RTMupJ0(:,camadT)*RTMdwJ0(:,camadT)*cdexp(-2.d0*uhJ0(:,camadT)))
 	AMdwJ1=(cdexp(-uJ1(:,camadT)*(prof(camadT)-h0)) + RTMupJ1(:,camadT)* &
 				cdexp(uJ1(:,camadT)*(prof(camadT-1)-h(camadT)-h0))) / &
 				(1.d0-RTMupJ1(:,camadT)*RTMdwJ1(:,camadT)*cdexp(-2.d0*uhJ1(:,camadT)))
 	AMupJ0=(cdexp(uJ0(:,camadT)*(prof(camadT-1)-h0))+RTMdwJ0(:,camadT)* &
 				cdexp(-uJ0(:,camadT)*(prof(camadT)+h(camadT)-h0))) / &
 				(1.d0-RTMupJ0(:,camadT)*RTMdwJ0(:,camadT)*cdexp(-2.d0*uhJ0(:,camadT)))
 	AMupJ1=(cdexp(uJ1(:,camadT)*(prof(camadT-1)-h0))+RTMdwJ1(:,camadT)* &
 				cdexp(-uJ1(:,camadT)*(prof(camadT)+h(camadT)-h0))) / &
 				(1.d0-RTMupJ1(:,camadT)*RTMdwJ1(:,camadT)*cdexp(-2.d0*uhJ1(:,camadT)))

 	if (camad > camadT) then
 	allocate(TMdwJ0(nJ0,camadT:camad),TMdwJ1(nJ1,camadT:camad))
 		do j=camadT,camad
 			if (j == camadT) then
 			TMdwJ0(:,j)= Iw * dsz / (2.d0*uJ0(:,camadT))
 			TMdwJ1(:,j)= Iw * dsz / (2.d0*uJ1(:,camadT))
 			elseif (j == (camadT + 1) .and. j == n) then
 			TMdwJ0(:,j)=TMdwJ0(:,j-1)*(cdexp(-uJ0(:,camadT)*(prof(camadT)-h0)) + &
				RTMupJ0(:,camadT)*AMupJ0(:)*cdexp(-uhJ0(:,camadT)) + &
				RTMdwJ0(:,camadT)*AMdwJ0(:))
 			TMdwJ1(:,j)=TMdwJ1(:,j-1)*(cdexp(-uJ1(:,camadT)*(prof(camadT)-h0)) + &
				RTMupJ1(:,camadT)*AMupJ1(:)*cdexp(-uhJ1(:,camadT)) + &
				RTMdwJ1(:,camadT)*AMdwJ1(:))
 			elseif (j == (camadT + 1) .and. j /= n) then
 			TMdwJ0(:,j)=TMdwJ0(:,j-1)*(cdexp(-uJ0(:,camadT)*(prof(camadT)-h0)) + &
				RTMupJ0(:,camadT)*AMupJ0(:)*cdexp(-uhJ0(:,camadT)) + &
				RTMdwJ0(:,camadT)*AMdwJ0(:)) / (1.d0 + RTMdwJ0(:,j)*cdexp(-2.d0*uhJ0(:,j)))
 			TMdwJ1(:,j)=TMdwJ1(:,j-1)*(cdexp(-uJ1(:,camadT)*(prof(camadT)-h0)) + &
				RTMupJ1(:,camadT)*AMupJ1(:)*cdexp(-uhJ1(:,camadT)) + &
				RTMdwJ1(:,camadT)*AMdwJ1(:)) / (1.d0 + RTMdwJ1(:,j)*cdexp(-2.d0*uhJ1(:,j)))
 			elseif (j /= n) then
 			TMdwJ0(:,j)=TMdwJ0(:,j-1)*(1.d0+RTMdwJ0(:,j-1))*cdexp(-uhJ0(:,j-1))/ &
 						(1.d0 + RTMdwJ0(:,j)*cdexp(-2.d0*uhJ0(:,j)))
 			TMdwJ1(:,j)=TMdwJ1(:,j-1)*(1.d0+RTMdwJ1(:,j-1))*cdexp(-uhJ1(:,j-1))/ &
 						(1.d0 + RTMdwJ1(:,j)*cdexp(-2.d0*uhJ1(:,j)))
 			elseif (j==n) then
 			TMdwJ0(:,j)=TMdwJ0(:,j-1)*(1.d0+RTMdwJ0(:,j-1))*cdexp(-uhJ0(:,j-1))
 			TMdwJ1(:,j)=TMdwJ1(:,j-1)*(1.d0+RTMdwJ1(:,j-1))*cdexp(-uhJ1(:,j-1))
 			end if
 		end do
 	elseif (camad < camadT) then
 	allocate(TMupJ0(nJ0,camad:camadT),TMupJ1(nJ1,camad:camadT))
 		do j=camadT,camad,-1
 			if (j == camadT) then
 			TMupJ0(:,j)= Iw * dsz / (2.d0*uJ0(:,camadT))
 			TMupJ1(:,j)= Iw * dsz / (2.d0*uJ1(:,camadT))
 			elseif (j == (camadT - 1) .and. j == 0) then
 			TMupJ0(:,j)=TMupJ0(:,j+1)*(cdexp(-uJ0(:,camadT)*h0) + &
				RTMupJ0(:,camadT)*AMupJ0(:) + RTMdwJ0(:,camadT)*AMdwJ0(:)* &
				cdexp(-uhJ0(:,camadT)))
 			TMupJ1(:,j)=TMupJ1(:,j+1)*(cdexp(-uJ1(:,camadT)*h0) + &
				RTMupJ1(:,camadT)*AMupJ1(:) + RTMdwJ1(:,camadT)*AMdwJ1(:)* &
				cdexp(-uhJ1(:,camadT)))
 			elseif (j == (camadT - 1) .and. j /= 0) then
 			TMupJ0(:,j)=TMupJ0(:,j+1)*(cdexp(uJ0(:,camadT)*(prof(camadT-1)-h0)) + &
				RTMupJ0(:,camadT)*AMupJ0(:) + RTMdwJ0(:,camadT)*AMdwJ0(:)* &
				cdexp(-uhJ0(:,camadT))) / (1.d0+RTMupJ0(:,j)*cdexp(-2.d0*uhJ0(:,j)))
 			TMupJ1(:,j)=TMupJ1(:,j+1)*(cdexp(uJ1(:,camadT)*(prof(camadT-1)-h0)) + &
				RTMupJ1(:,camadT)*AMupJ1(:) + RTMdwJ1(:,camadT)*AMdwJ1(:)* &
				cdexp(-uhJ1(:,camadT))) / (1.d0+RTMupJ1(:,j)*cdexp(-2.d0*uhJ1(:,j)))
 			elseif (j /= 0) then
 			TMupJ0(:,j)=TMupJ0(:,j+1)*(1.d0+RTMupJ0(:,j+1))*cdexp(-uhJ0(:,j+1)) / &
				(1.d0 + RTMupJ0(:,j)*cdexp(-2.d0*uhJ0(:,j)))
 			TMupJ1(:,j)=TMupJ1(:,j+1)*(1.d0+RTMupJ1(:,j+1))*cdexp(-uhJ1(:,j+1)) / &
				(1.d0 + RTMupJ1(:,j)*cdexp(-2.d0*uhJ1(:,j)))
 			elseif (j == 0) then
 			TMupJ0(:,j)=TMupJ0(:,1)*(1.d0+RTMupJ0(:,1))*cdexp(-uhJ0(:,1))
 			TMupJ1(:,j)=TMupJ1(:,1)*(1.d0+RTMupJ1(:,1))*cdexp(-uhJ1(:,1))
 			end if
 		end do
 	else
 	allocate(TMdwJ0(nJ0,camadT:camad),TMdwJ1(nJ1,camadT:camad))
 	allocate(TMupJ0(nJ0,camad:camadT),TMupJ1(nJ1,camad:camadT))
 			TMdwJ0(:,camad)= Iw * dsz / (2.d0*uJ0(:,camadT))
 			TMdwJ1(:,camad)= Iw * dsz / (2.d0*uJ1(:,camadT))
 			TMupJ0(:,camad)= Iw * dsz / (2.d0*uJ0(:,camadT))
 			TMupJ1(:,camad)= Iw * dsz / (2.d0*uJ1(:,camadT))
 	end if

 	allocate(Ktmdz_J1(nJ1),Ktm_J0(nJ0),Ktm_J1(nJ1))
 	allocate(kernelExJ1(nJ1),kernelEyJ1(nJ1),kernelEzJ0(nJ0))
 	allocate(kernelHxJ1(nJ1),kernelHyJ1(nJ1))
 	if (camad == 0 .and. camadT /= 0) then
 		Ktmdz_J1=(ImpIntJ1(:,0)*TMupJ1(:,0)*cdexp(uJ1(:,0)*z))*w_J1(:)
 		Ktm_J0=(TMupJ0(:,0)*cdexp(uJ0(:,0)*z))*w_J0(:)
 		Ktm_J1=(TMupJ1(:,0)*cdexp(uJ1(:,0)*z))*w_J1(:)

 		kernelExJ1 = x/r * Ktmdz_J1 * krJ1 * krJ1
 		Ex_p = - sum(kernelExJ1) / (2.d0*pi*r)		!este último r é decorrente do uso dos filtros
 
 		kernelEyJ1 = y/r * Ktmdz_J1 * krJ1 * krJ1
 		Ey_p = - sum(kernelEyJ1) / (2.d0*pi*r)		!este último r é decorrente do uso dos filtros
 
 		kernelEzJ0 = Ktm_J0 * krJ0 * krJ0 * krJ0
 		Ez_p = sum(kernelEzJ0) / (2.d0*pi*r*neta)	!este último r é decorrente do uso dos filtros

 		kernelHxJ1 = y/r * Ktm_J1 * krJ1 * krJ1
 		Hx_p = - sum(kernelHxJ1) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros

 		kernelHyJ1 = x/r * Ktm_J1 * krJ1 * krJ1
 		Hy_p = sum(kernelHyJ1) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 	elseif (camad < camadT) then !camada k
 		ktmdz_J1=(ImpIntJ1(:,camad)*TMupJ1(:,camad)*(cdexp(uJ1(:,camad)*(z-prof(camad))) - &
 			RTMupJ1(:,camad)*cdexp(-uJ1(:,camad)*(z-prof(camad-1)+h(camad)))))*w_J1(:)
 		ktm_J0=(TMupJ0(:,camad)*(cdexp(uJ0(:,camad)*(z-prof(camad))) + &
 			RTMupJ0(:,camad)*cdexp(-uJ0(:,camad)*(z-prof(camad-1)+h(camad)))))*w_J0(:)
 		ktm_J1=(TMupJ1(:,camad)*(cdexp(uJ1(:,camad)*(z-prof(camad))) + &
 			RTMupJ1(:,camad)*cdexp(-uJ1(:,camad)*(z-prof(camad-1)+h(camad)))))*w_J1(:)
 
 		kernelExJ1 = x/r * Ktmdz_J1 * krJ1 * krJ1
 		Ex_p = - sum(kernelExJ1) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelEyJ1 = y/r * Ktmdz_J1 * krJ1 * krJ1
 		Ey_p = - sum(kernelEyJ1) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelEzJ0 = Ktm_J0 * krJ0 * krJ0 * krJ0
 		Ez_p = sum(kernelEzJ0) / (2.d0*pi*r*condut(camad))	!este último r é decorrente do uso dos filtros
 
 		kernelHxJ1 = y/r * Ktm_J1 * krJ1 * krJ1
 		Hx_p = - sum(kernelHxJ1) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelHyJ1 = x/r * Ktm_J1 * krJ1 * krJ1
 		Hy_p = sum(kernelHyJ1) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 	elseif (camad == camadT .and. z <= h0) then	!na mesma camada do transmissor mas acima dele
 		Ktmdz_J1=(ImpIntJ1(:,camad)*TMupJ1(:,camad)*(cdexp(uJ1(:,camad)*(z-h0)) - &
 			RTMupJ1(:,camad)*AMupJ1(:)*cdexp(-uJ1(:,camad)*(z-prof(camad-1))) + &
 			RTMdwJ1(:,camad)*AMdwJ1(:)*cdexp(uJ1(:,camad)*(z-prof(camad)))))*w_J1(:)
 		Ktm_J0=(TMupJ0(:,camad)*(cdexp(uJ0(:,camad)*(z-h0)) + &
 			RTMupJ0(:,camad)*AMupJ0(:)*cdexp(-uJ0(:,camad)*(z-prof(camad-1))) + &
 			RTMdwJ0(:,camad)*AMdwJ0(:)*cdexp(uJ0(:,camad)*(z-prof(camad)))))*w_J0(:)
 		Ktm_J1=(TMupJ1(:,camad)*(cdexp(uJ1(:,camad)*(z-h0)) + &
 			RTMupJ1(:,camad)*AMupJ1(:)*cdexp(-uJ1(:,camad)*(z-prof(camad-1))) + &
 			RTMdwJ1(:,camad)*AMdwJ1(:)*cdexp(uJ1(:,camad)*(z-prof(camad)))))*w_J1(:)

 		kernelExJ1 = x/r * Ktmdz_J1 * krJ1 * krJ1
 		Ex_p = - sum(kernelExJ1) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelEyJ1 = y/r * Ktmdz_J1 * krJ1 * krJ1
 		Ey_p = - sum(kernelEyJ1) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
		if (camad /= 0) then
 		kernelEzJ0 = Ktm_J0 * krJ0 * krJ0 * krJ0
 		Ez_p = sum(kernelEzJ0) / (2.d0*pi*r*condut(camad))	!este último r é decorrente do uso dos filtros
		else
 		kernelEzJ0 = Ktm_J0 * krJ0 * krJ0 * krJ0
 		Ez_p = sum(kernelEzJ0) / (2.d0*pi*r*neta)	!este último r é decorrente do uso dos filtros
		end if
 
 		kernelHxJ1 = y/r * Ktm_J1 * krJ1 * krJ1
 		Hx_p = - sum(kernelHxJ1) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelHyJ1 = x/r * Ktm_J1 * krJ1 * krJ1
 		Hy_p = sum(kernelHyJ1) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 	elseif (camad == camadT .and. z > h0) then	!na mesma camada do transmissor mas abaixo dele
 		Ktmdz_J1=(ImpIntJ1(:,camad)*TMdwJ1(:,camad)*(cdexp(-uJ1(:,camad)*(z-h0)) + &
 			RTMupJ1(:,camad)*AMupJ1(:)*cdexp(-uJ1(:,camad)*(z-prof(camad-1))) - &
 			RTMdwJ1(:,camad)*AMdwJ1(:)*cdexp(uJ1(:,camad)*(z-prof(camad)))))*w_J1(:)
		Ktm_J0=(TMdwJ0(:,camad)*(cdexp(-uJ0(:,camad)*(z-h0)) + &
 			RTMupJ0(:,camad)*AMupJ0(:)*cdexp(-uJ0(:,camad)*(z-prof(camad-1))) + &
 			RTMdwJ0(:,camad)*AMdwJ0(:)*cdexp(uJ0(:,camad)*(z-prof(camad)))))*w_J0(:)
 		Ktm_J1=(TMdwJ1(:,camad)*(cdexp(-uJ1(:,camad)*(z-h0)) + &
 			RTMupJ1(:,camad)*AMupJ1(:)*cdexp(-uJ1(:,camad)*(z-prof(camad-1))) + &
 			RTMdwJ1(:,camad)*AMdwJ1(:)*cdexp(uJ1(:,camad)*(z-prof(camad)))))*w_J1(:)
 
 		kernelExJ1 = x/r * Ktmdz_J1 * krJ1 * krJ1
 		Ex_p = sum(kernelExJ1) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelEyJ1 = y/r * Ktmdz_J1 * krJ1 * krJ1
 		Ey_p = sum(kernelEyJ1) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
		if (camad /= 0) then
 		kernelEzJ0 = Ktm_J0 * krJ0 * krJ0 * krJ0
 		Ez_p = sum(kernelEzJ0) / (2.d0*pi*r*condut(camad))	!este último r é decorrente do uso dos filtros
		else
 		kernelEzJ0 = Ktm_J0 * krJ0 * krJ0 * krJ0
 		Ez_p = sum(kernelEzJ0) / (2.d0*pi*r*neta)	!este último r é decorrente do uso dos filtros
		end if
 
 		kernelHxJ1 = y/r * Ktm_J1 * krJ1 * krJ1
 		Hx_p = - sum(kernelHxJ1) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelHyJ1 = x/r * Ktm_J1 * krJ1 * krJ1
 		Hy_p = sum(kernelHyJ1) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 	elseif (camad > camadT .and. camad /= n) then !camada j
 		Ktmdz_J1=(ImpIntJ1(:,camad)*TMdwJ1(:,camad)*(cdexp(-uJ1(:,camad)*(z-prof(camad-1))) - &
 			RTMdwJ1(:,camad)*cdexp(uJ1(:,camad)*(z-prof(camad)-h(camad)))))*w_J1(:)
 		Ktm_J0=(TMdwJ0(:,camad)*(cdexp(-uJ0(:,camad)*(z-prof(camad-1))) + &
 			RTMdwJ0(:,camad)*cdexp(uJ0(:,camad)*(z-prof(camad)-h(camad)))))*w_J0(:)
 		Ktm_J1=(TMdwJ1(:,camad)*(cdexp(-uJ1(:,camad)*(z-prof(camad-1))) + &
 			RTMdwJ1(:,camad)*cdexp(uJ1(:,camad)*(z-prof(camad)-h(camad)))))*w_J1(:)

 		kernelExJ1 = x/r * Ktmdz_J1 * krJ1 * krJ1
 		Ex_p = sum(kernelExJ1) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelEyJ1 = y/r * Ktmdz_J1 * krJ1 * krJ1
 		Ey_p = sum(kernelEyJ1) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelEzJ0 = Ktm_J0 * krJ0 * krJ0 * krJ0
 		Ez_p = sum(kernelEzJ0) / (2.d0*pi*r*condut(camad))	!este último r é decorrente do uso dos filtros
 
 		kernelHxJ1 = y/r * Ktm_J1 * krJ1 * krJ1
 		Hx_p = - sum(kernelHxJ1) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelHyJ1 = x/r * Ktm_J1 * krJ1 * krJ1
 		Hy_p = sum(kernelHyJ1) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 	else	!camada n
		Ktmdz_J1=(ImpIntJ1(:,n)*TMdwJ1(:,n)*cdexp(-uJ1(:,n)*(z-prof(n-1))))*w_J1(:)
		Ktm_J0=(TMdwJ0(:,n)*cdexp(-uJ0(:,n)*(z-prof(n-1))))*w_J0(:)
		Ktm_J1=(TMdwJ1(:,n)*cdexp(-uJ1(:,n)*(z-prof(n-1))))*w_J1(:)

 		kernelExJ1 = x/r * Ktmdz_J1 * krJ1 * krJ1
 		Ex_p = sum(kernelExJ1) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros

 		kernelEyJ1 = y/r * Ktmdz_J1 * krJ1 * krJ1
 		Ey_p = sum(kernelEyJ1) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros

 		kernelEzJ0 = Ktm_J0 * krJ0 * krJ0 * krJ0
 		Ez_p = sum(kernelEzJ0) / (2.d0*pi*r*condut(camad))	!este último r é decorrente do uso dos filtros

 		kernelHxJ1 = y/r * Ktm_J1 * krJ1 * krJ1
 		Hx_p = - sum(kernelHxJ1) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros

 		kernelHyJ1 = x/r * Ktm_J1 * krJ1 * krJ1
 		Hy_p = sum(kernelHyJ1) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros

 	end if
 
	deallocate(h,KrJ0,KrJ1,w_J0,w_J1)
	deallocate(wvnb2,uJ0,uJ1,ImpIntJ0,ImpIntJ1,uhJ0,uhJ1,tghJ0,tghJ1)
	deallocate(ImpApdwJ0,ImpApdwJ1,RTMdwJ0,RTMdwJ1)
	deallocate(ImpApupJ0,ImpApupJ1,RTMupJ0,RTMupJ1)
	deallocate(AMdwJ0,AMdwJ1,AMupJ0,AMupJ1)
	deallocate(Ktmdz_J1,Ktm_J0,Ktm_J1)
	deallocate(kernelExJ1,kernelEyJ1,kernelEzJ0)
	deallocate(kernelHxJ1,kernelHyJ1)
	end subroutine dev_xyz_loops
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

	subroutine dmhx_xyz_loops( idtfcd_cJ0, nJ0, idtfcd_cJ1, nJ1, Tx, Ty, Iw, mx, h0, n, esp, condut, &
						neta, zeta, cx, cy, z, Ex_p, Ey_p, Ez_p, Hx_p, Hy_p, Hz_p)
	implicit none
	integer(4),intent(in)::n, idtfcd_cJ0, nJ0, idtfcd_cJ1, nJ1
	real(8),intent(in)::Tx,Ty,Iw,mx,h0,esp(:),condut(1:n),cx,cy,z  !,neta
	complex*16,intent(in)::zeta,neta
	complex*16,intent(out)::Ex_p,Ey_p,Ez_p,Hx_p,Hy_p,Hz_p

	real(8),parameter::pi=3.141592653589793238462643383279502884197d0
	integer(4)::i,j,k,camad,camadT,filtro,ident_fJ0,ident_fJ1
	real(8)::x,y,r,kr,k_r
	real(8),dimension(:),allocatable::h,krJ0,krJ1,w_J0,w_J1,prof

!	Para uso de loops:
	complex*16,dimension(:),allocatable::wvnb2,AMdwJ0,AMdwJ1,AMupJ0,AMupJ1,FEdwJ0,FEdwJ1,FEupJ0,FEupJ1
	complex*16,dimension(:,:),allocatable::uJ0,uJ1,AdmIntJ0,AdmIntJ1,ImpIntJ0,ImpIntJ1,uhJ0,uhJ1,tghJ0,tghJ1
	complex*16,dimension(:,:),allocatable::AdmApdwJ0,AdmApdwJ1,ImpApdwJ0,ImpApdwJ1,AdmApupJ0,AdmApupJ1,ImpApupJ0,ImpApupJ1
	complex*16,dimension(:,:),allocatable::RTEdwJ0,RTEdwJ1,RTMdwJ0,RTMdwJ1,RTEupJ0,RTEupJ1,RTMupJ0,RTMupJ1
	complex*16,dimension(:,:),allocatable::TMdwJ0,TMdwJ1,TEdwJ0,TEdwJ1,TMupJ0,TMupJ1,TEupJ0,TEupJ1

	complex*16,dimension(:),allocatable::Ktmdz_J0,Ktmdz_J1,Ktm_J0,Ktm_J1,Kte_J0,Kte_J1,Ktedz_J0,Ktedz_J1
	complex*16,dimension(:),allocatable::kernelExJ0,kernelExJ1,kernelEyJ0,kernelEyJ1,kernelEzJ1
	complex*16,dimension(:),allocatable::kernelHxJ0,kernelHxJ1,kernelHyJ0,kernelHyJ1,kernelHzJ1

	complex*16::kerEx_J1,kerEx_J0,kerEy_J1,kerEy_J0,kerEz_J1
	complex*16::kerHx_J1,kerHx_J0,kerHy_J1,kerHy_J0,kerHz_J1

	if (cx==Tx) then
	x=1.d0!Tx/dabs(Tx)  !0.1*Tx/dabs(Tx)
	else
	x=cx-Tx
	end if
	y=cy-Ty
	r = dsqrt(x**2 + y**2)

	allocate(h(0:n),prof(-1:n))
	if (size(esp)==n) then
		h(0)=0.d0
		h(1:n)=esp
	else
		h(0)=0.d0
		h(1:n-1)=esp
		h(n)=1.d300
	end if
!	criando um novo vetor de profundidades que se adeque à qualquer situação patológica
	prof(-1)=-1.d300
	prof(0)=0.d0
	if (n > 1) then
	 prof(1) = h(1)
	 if (n > 2) then
	  do k=2,n-1
	   prof(k) = prof(k-1) + h(k)
	  end do
	 end if
	end if
	prof(n)=1.d300

!para descobrir em que camada está a observação
	if (z <= 0.d0) then
		camad=0
	else if (z > prof(n-1)) then
		camad=n
	else
	do i=n-1,1,-1
	if (z > prof(i-1)) then
		camad=i
		exit
	end if
	end do
	end if

!para descobrir em que camada está o transmissor
	if (h0 <= 0.d0) then
		camadT = 0
	else if (h0 > prof(n-1)) then
		camadT = n
	else
	do j=n-1,1,-1
	if (h0 > prof(j-1)) then
		camadT = j
		exit
	end if
	end do
	end if

!!	write(*,*)'Entre com o criador dos filtros J0: Rijo(0), Frayzer(1), Guptasarma(2), Kong(3) ou Key(4)'
!!	read(*,*)idtfcd_cJ0
	filtro=0	!esta variável direciona o uso de filtros J0 e J1 em vez de seno e cosseno
!!	call identfiltro(filtro,idtfcd_cJ0,ident_fJ0,nJ0)
!!	write(*,*)'Entre com o criador dos filtros J1: Rijo(0), Frayzer(1), Guptasarma(2), Kong(3) ou Key(4)'
!!	read(*,*)idtfcd_cJ1
!!	call identfiltro(filtro,idtfcd_cJ1,ident_fJ1,nJ1)

!	idtfcd_cJ0 = 4
	ident_fJ0 = 0
!	nJ0 = 101
!	idtfcd_cJ1 = 4
	ident_fJ1 = 1
!	nJ1 = 101

	allocate(KrJ0(nJ0),KrJ1(nJ1),w_J0(nJ0),w_J1(nJ1))

	call constfiltro(filtro,idtfcd_cJ0,ident_fJ0,nJ0,r,KrJ0,w_J0)
	call constfiltro(filtro,idtfcd_cJ1,ident_fJ1,nJ1,r,KrJ1,w_J1)

 	allocate(wvnb2(0:n),uJ0(nJ0,0:n),uJ1(nJ1,0:n),AdmIntJ0(nJ0,0:n),AdmIntJ1(nJ1,0:n),ImpIntJ0(nJ0,0:n),ImpIntJ1(nJ1,0:n))
 	allocate(uhJ0(nJ0,0:n),uhJ1(nJ1,0:n),tghJ0(nJ0,0:n),tghJ1(nJ1,0:n))
! 
 	do i=0,n
 		if (i==0) then
 		wvnb2(i)=-zeta*neta
 		uJ0(:,i)=cdsqrt(krJ0*krJ0-wvnb2(i))
 		uJ1(:,i)=cdsqrt(krJ1*krJ1-wvnb2(i))
 		AdmIntJ0(:,i)=uJ0(:,i)/zeta
 		AdmIntJ1(:,i)=uJ1(:,i)/zeta
 		ImpIntJ0(:,i)=uJ0(:,i)/neta
 		ImpIntJ1(:,i)=uJ1(:,i)/neta
 		uhJ0(:,i)=uJ0(:,i)*h(i)
 		uhJ1(:,i)=uJ1(:,i)*h(i)
 		tghJ0(:,i)=(1.d0-cdexp(-2.d0*uhJ0(:,i)))/(1.d0+cdexp(-2.d0*uhJ0(:,i)))
 		tghJ1(:,i)=(1.d0-cdexp(-2.d0*uhJ1(:,i)))/(1.d0+cdexp(-2.d0*uhJ1(:,i)))
 		else
 		wvnb2(i) = -zeta*condut(i)
 		uJ0(:,i)=cdsqrt(krJ0*krJ0-wvnb2(i))
 		uJ1(:,i)=cdsqrt(krJ1*krJ1-wvnb2(i))
 		AdmIntJ0(:,i)=uJ0(:,i)/zeta
 		AdmIntJ1(:,i)=uJ1(:,i)/zeta
 		ImpIntJ0(:,i)=uJ0(:,i)/condut(i)
 		ImpIntJ1(:,i)=uJ1(:,i)/condut(i)
 		uhJ0(:,i)=uJ0(:,i)*h(i)
 		uhJ1(:,i)=uJ1(:,i)*h(i)
 		tghJ0(:,i)=(1.d0-cdexp(-2.d0*uhJ0(:,i)))/(1.d0+cdexp(-2.d0*uhJ0(:,i)))
 		tghJ1(:,i)=(1.d0-cdexp(-2.d0*uhJ1(:,i)))/(1.d0+cdexp(-2.d0*uhJ1(:,i)))
 		end if
 	end do
 
 	allocate(AdmApdwJ0(nJ0,1:n),AdmApdwJ1(nJ1,1:n),ImpApdwJ0(nJ0,1:n),ImpApdwJ1(nJ1,1:n))
 	allocate(RTEdwJ0(nJ0,0:n),RTEdwJ1(nJ1,0:n),RTMdwJ0(nJ0,0:n),RTMdwJ1(nJ1,0:n))
 
 	do i=n,1,-1
 		if (i==n) then
 		AdmApdwJ0(:,i)=AdmIntJ0(:,i)
 		AdmApdwJ1(:,i)=AdmIntJ1(:,i)
 		ImpApdwJ0(:,i)=ImpIntJ0(:,i)
 		ImpApdwJ1(:,i)=ImpIntJ1(:,i)
 		RTEdwJ0(:,i)=(0.d0,0.d0)
 		RTEdwJ1(:,i)=(0.d0,0.d0)
 		RTMdwJ0(:,i)=(0.d0,0.d0)
 		RTMdwJ1(:,i)=(0.d0,0.d0)
 		else
 		AdmApdwJ0(:,i)=AdmIntJ0(:,i)*(AdmApdwJ0(:,i+1)+AdmIntJ0(:,i)* &
 			tghJ0(:,i))/(AdmIntJ0(:,i)+AdmApdwJ0(:,i+1)*tghJ0(:,i))
 		AdmApdwJ1(:,i)=AdmIntJ1(:,i)*(AdmApdwJ1(:,i+1)+AdmIntJ1(:,i)* &
 			tghJ1(:,i))/(AdmIntJ1(:,i)+AdmApdwJ1(:,i+1)*tghJ1(:,i))
 		ImpApdwJ0(:,i)=ImpIntJ0(:,i)*(ImpApdwJ0(:,i+1)+ImpIntJ0(:,i)* &
 			tghJ0(:,i))/(ImpIntJ0(:,i)+ImpApdwJ0(:,i+1)*tghJ0(:,i))
 		ImpApdwJ1(:,i)=ImpIntJ1(:,i)*(ImpApdwJ1(:,i+1)+ImpIntJ1(:,i)* &
 			tghJ1(:,i))/(ImpIntJ1(:,i)+ImpApdwJ1(:,i+1)*tghJ1(:,i))
 		RTEdwJ0(:,i)=(AdmIntJ0(:,i)-AdmApdwJ0(:,i+1))/(AdmIntJ0(:,i)+AdmApdwJ0(:,i+1))
 		RTEdwJ1(:,i)=(AdmIntJ1(:,i)-AdmApdwJ1(:,i+1))/(AdmIntJ1(:,i)+AdmApdwJ1(:,i+1))
 		RTMdwJ0(:,i)=(ImpIntJ0(:,i)-ImpApdwJ0(:,i+1))/(ImpIntJ0(:,i)+ImpApdwJ0(:,i+1))
 		RTMdwJ1(:,i)=(ImpIntJ1(:,i)-ImpApdwJ1(:,i+1))/(ImpIntJ1(:,i)+ImpApdwJ1(:,i+1))
 		end if	
 	end do
 		RTEdwJ0(:,0)=(AdmIntJ0(:,0)-AdmApdwJ0(:,1))/(AdmIntJ0(:,0)+AdmApdwJ0(:,1))
 		RTEdwJ1(:,0)=(AdmIntJ1(:,0)-AdmApdwJ1(:,1))/(AdmIntJ1(:,0)+AdmApdwJ1(:,1))
 		RTMdwJ0(:,0)=(ImpIntJ0(:,0)-ImpApdwJ0(:,1))/(ImpIntJ0(:,0)+ImpApdwJ0(:,1))
 		RTMdwJ1(:,0)=(ImpIntJ1(:,0)-ImpApdwJ1(:,1))/(ImpIntJ1(:,0)+ImpApdwJ1(:,1))
 
 	allocate(AdmApupJ0(nJ0,0:n-1),AdmApupJ1(nJ1,0:n-1),ImpApupJ0(nJ0,0:n-1),ImpApupJ1(nJ1,0:n-1))
 	allocate(RTEupJ0(nJ0,0:n),RTEupJ1(nJ1,0:n),RTMupJ0(nJ0,0:n),RTMupJ1(nJ1,0:n))
 
 	do i=0,n-1
 		if (i==0) then
 		AdmApupJ0(:,i)=AdmIntJ0(:,i)
 		AdmApupJ1(:,i)=AdmIntJ1(:,i)
 		ImpApupJ0(:,i)=ImpIntJ0(:,i)
 		ImpApupJ1(:,i)=ImpIntJ1(:,i)
 		RTEupJ0(:,i)=(0.d0,0.d0)
 		RTEupJ1(:,i)=(0.d0,0.d0)
 		RTMupJ0(:,i)=(0.d0,0.d0)
 		RTMupJ1(:,i)=(0.d0,0.d0)
 		else
 		AdmApupJ0(:,i)=AdmIntJ0(:,i)*(AdmApupJ0(:,i-1)+AdmIntJ0(:,i)* &
 			tghJ0(:,i))/(AdmIntJ0(:,i)+AdmApupJ0(:,i-1)*tghJ0(:,i))
 		AdmApupJ1(:,i)=AdmIntJ1(:,i)*(AdmApupJ1(:,i-1)+AdmIntJ1(:,i)* &
 			tghJ1(:,i))/(AdmIntJ1(:,i)+AdmApupJ1(:,i-1)*tghJ1(:,i))
 		ImpApupJ0(:,i)=ImpIntJ0(:,i)*(ImpApupJ0(:,i-1)+ImpIntJ0(:,i)* &
 			tghJ0(:,i))/(ImpIntJ0(:,i)+ImpApupJ0(:,i-1)*tghJ0(:,i))
 		ImpApupJ1(:,i)=ImpIntJ1(:,i)*(ImpApupJ1(:,i-1)+ImpIntJ1(:,i)* &
 			tghJ1(:,i))/(ImpIntJ1(:,i)+ImpApupJ1(:,i-1)*tghJ1(:,i))
 		RTEupJ0(:,i)=(AdmIntJ0(:,i)-AdmApupJ0(:,i-1))/(AdmIntJ0(:,i)+AdmApupJ0(:,i-1))
 		RTEupJ1(:,i)=(AdmIntJ1(:,i)-AdmApupJ1(:,i-1))/(AdmIntJ1(:,i)+AdmApupJ1(:,i-1))
 		RTMupJ0(:,i)=(ImpIntJ0(:,i)-ImpApupJ0(:,i-1))/(ImpIntJ0(:,i)+ImpApupJ0(:,i-1))
 		RTMupJ1(:,i)=(ImpIntJ1(:,i)-ImpApupJ1(:,i-1))/(ImpIntJ1(:,i)+ImpApupJ1(:,i-1))
 		end if
 	end do
 		RTEupJ0(:,n)=(AdmIntJ0(:,n)-AdmApupJ0(:,n-1))/(AdmIntJ0(:,n)+AdmApupJ0(:,n-1))
 		RTEupJ1(:,n)=(AdmIntJ1(:,n)-AdmApupJ1(:,n-1))/(AdmIntJ1(:,n)+AdmApupJ1(:,n-1))
 		RTMupJ0(:,n)=(ImpIntJ0(:,n)-ImpApupJ0(:,n-1))/(ImpIntJ0(:,n)+ImpApupJ0(:,n-1))
 		RTMupJ1(:,n)=(ImpIntJ1(:,n)-ImpApupJ1(:,n-1))/(ImpIntJ1(:,n)+ImpApupJ1(:,n-1))
 
 	allocate(AMdwJ0(nJ0),AMdwJ1(nJ1),AMupJ0(nJ0),AMupJ1(nJ1),FEdwJ0(nJ0),FEdwJ1(nJ1),FEupJ0(nJ0),FEupJ1(nJ1))

 	AMdwJ0=(cdexp(-uJ0(:,camadT)*(prof(camadT)-h0)) + RTMupJ0(:,camadT)* &
				cdexp(uJ0(:,camadT)*(prof(camadT-1)-h(camadT)-h0))) / &
				(1.d0 - RTMupJ0(:,camadT)*RTMdwJ0(:,camadT)*cdexp(-2.d0*uhJ0(:,camadT)))
 	AMdwJ1=(cdexp(-uJ1(:,camadT)*(prof(camadT)-h0)) + RTMupJ1(:,camadT)* &
 				cdexp(uJ1(:,camadT)*(prof(camadT-1)-h(camadT)-h0))) / &
 				(1.d0 - RTMupJ1(:,camadT)*RTMdwJ1(:,camadT)*cdexp(-2.d0*uhJ1(:,camadT)))

 	AMupJ0=(cdexp(uJ0(:,camadT)*(prof(camadT-1)-h0)) + RTMdwJ0(:,camadT)* &
 				cdexp(-uJ0(:,camadT)*(prof(camadT)+h(camadT)-h0))) / &
 				(1.d0 - RTMupJ0(:,camadT)*RTMdwJ0(:,camadT)*cdexp(-2.d0*uhJ0(:,camadT)))
 	AMupJ1=(cdexp(uJ1(:,camadT)*(prof(camadT-1)-h0)) + RTMdwJ1(:,camadT)* &
 				cdexp(-uJ1(:,camadT)*(prof(camadT)+h(camadT)-h0))) / &
 				(1.d0 - RTMupJ1(:,camadT)*RTMdwJ1(:,camadT)*cdexp(-2.d0*uhJ1(:,camadT)))

 	FEdwJ0=(cdexp(-uJ0(:,camadT)*(prof(camadT)-h0)) - RTEupJ0(:,camadT)* &
 				cdexp(uJ0(:,camadT)*(prof(camadT-1)-h(camadT)-h0))) / &
 				(1.d0 - RTEupJ0(:,camadT)*RTEdwJ0(:,camadT)*cdexp(-2.d0*uhJ0(:,camadT)))
 	FEdwJ1=(cdexp(-uJ1(:,camadT)*(prof(camadT)-h0)) - RTEupJ1(:,camadT)* &
 				cdexp(uJ1(:,camadT)*(prof(camadT-1)-h(camadT)-h0))) / &
 				(1.d0 - RTEupJ1(:,camadT)*RTEdwJ1(:,camadT)*cdexp(-2.d0*uhJ1(:,camadT)))

 	FEupJ0=(cdexp(uJ0(:,camadT)*(prof(camadT-1)-h0)) - RTEdwJ0(:,camadT)* &
 				cdexp(-uJ0(:,camadT)*(prof(camadT)+h(camadT)-h0))) / &
 				(1.d0 - RTEupJ0(:,camadT)*RTEdwJ0(:,camadT)*cdexp(-2.d0*uhJ0(:,camadT)))
 	FEupJ1=(cdexp(uJ1(:,camadT)*(prof(camadT-1)-h0)) - RTEdwJ1(:,camadT)* &
 				cdexp(-uJ1(:,camadT)*(prof(camadT)+h(camadT)-h0))) / &
 				(1.d0 - RTEupJ1(:,camadT)*RTEdwJ1(:,camadT)*cdexp(-2.d0*uhJ1(:,camadT)))

 	if (camad > camadT) then
 	allocate(TMdwJ0(nJ0,camadT:camad),TMdwJ1(nJ1,camadT:camad),TEdwJ0(nJ0,camadT:camad),TEdwJ1(nJ1,camadT:camad))
 		do j=camadT,camad
 			if (j == camadT) then

 			TMdwJ0(:,j)= zeta * mx / (2.d0 * ImpIntJ0(:,camadT))
 			TMdwJ1(:,j)= zeta * mx / (2.d0 * ImpIntJ1(:,camadT))
 			TEdwJ0(:,j)= - zeta * mx / 2.d0
 			TEdwJ1(:,j)= - zeta * mx / 2.d0

 			elseif (j == (camadT + 1) .and. j == n) then

 			TMdwJ0(:,j)=TMdwJ0(:,j-1)*(cdexp(-uJ0(:,camadT)*(prof(camadT)-h0)) + &
				RTMupJ0(:,camadT)*AMupJ0(:)*cdexp(-uhJ0(:,camadT)) + &
				RTMdwJ0(:,camadT)*AMdwJ0(:))
 			TMdwJ1(:,j)=TMdwJ1(:,j-1)*(cdexp(-uJ1(:,camadT)*(prof(camadT)-h0)) + &
				RTMupJ1(:,camadT)*AMupJ1(:)*cdexp(-uhJ1(:,camadT)) + &
				RTMdwJ1(:,camadT)*AMdwJ1(:))
 			TEdwJ0(:,j)=TEdwJ0(:,j-1)*(cdexp(-uJ0(:,camadT)*(prof(camadT)-h0)) - &
				RTEupJ0(:,camadT)*FEupJ0(:)*cdexp(-uhJ0(:,camadT)) + &
				RTEdwJ0(:,camadT)*FEdwJ0(:))
			TEdwJ1(:,j)=TEdwJ1(:,j-1)*(cdexp(-uJ1(:,camadT)*(prof(camadT)-h0)) - &
				RTEupJ1(:,camadT)*FEupJ1(:)*cdexp(-uhJ1(:,camadT)) + &
				RTEdwJ1(:,camadT)*FEdwJ1(:))

 			elseif (j == (camadT + 1) .and. j /= n) then

 			TMdwJ0(:,j)=TMdwJ0(:,j-1)*(cdexp(-uJ0(:,camadT)*(prof(camadT)-h0)) + &
				RTMupJ0(:,camadT)*AMupJ0(:)*cdexp(-uhJ0(:,camadT)) + &
				RTMdwJ0(:,camadT)*AMdwJ0(:)) / (1.d0 + RTMdwJ0(:,j)*cdexp(-2.d0*uhJ0(:,j)))
 			TMdwJ1(:,j)=TMdwJ1(:,j-1)*(cdexp(-uJ1(:,camadT)*(prof(camadT)-h0)) + &
				RTMupJ1(:,camadT)*AMupJ1(:)*cdexp(-uhJ1(:,camadT)) + &
				RTMdwJ1(:,camadT)*AMdwJ1(:)) / (1.d0 + RTMdwJ1(:,j)*cdexp(-2.d0*uhJ1(:,j)))
 			TEdwJ0(:,j)=TEdwJ0(:,j-1)*(cdexp(-uJ0(:,camadT)*(prof(camadT)-h0)) - &
				RTEupJ0(:,camadT)*FEupJ0(:)*cdexp(-uhJ0(:,camadT)) + &
				RTEdwJ0(:,camadT)*FEdwJ0(:)) / (1.d0 + RTEdwJ0(:,j)*cdexp(-2.d0*uhJ0(:,j)))
 			TEdwJ1(:,j)=TEdwJ1(:,j-1)*(cdexp(-uJ1(:,camadT)*(prof(camadT)-h0)) - &
				RTEupJ1(:,camadT)*FEupJ1(:)*cdexp(-uhJ1(:,camadT)) + &
				RTEdwJ1(:,camadT)*FEdwJ1(:)) / (1.d0 + RTEdwJ1(:,j)*cdexp(-2.d0*uhJ1(:,j)))

 			elseif (j /= n) then

 			TMdwJ0(:,j)=TMdwJ0(:,j-1)*(1.d0+RTMdwJ0(:,j-1))*cdexp(-uhJ0(:,j-1))/ &
 						(1.d0 + RTMdwJ0(:,j)*cdexp(-2.d0*uhJ0(:,j)))
 			TMdwJ1(:,j)=TMdwJ1(:,j-1)*(1.d0+RTMdwJ1(:,j-1))*cdexp(-uhJ1(:,j-1))/ &
 						(1.d0 + RTMdwJ1(:,j)*cdexp(-2.d0*uhJ1(:,j)))
 			TEdwJ0(:,j)=TEdwJ0(:,j-1)*(1.d0+RTEdwJ0(:,j-1))*cdexp(-uhJ0(:,j-1))/ &
 						(1.d0 + RTEdwJ0(:,j)*cdexp(-2.d0*uhJ0(:,j)))
 			TEdwJ1(:,j)=TEdwJ1(:,j-1)*(1.d0+RTEdwJ1(:,j-1))*cdexp(-uhJ1(:,j-1))/ &
 						(1.d0 + RTEdwJ1(:,j)*cdexp(-2.d0*uhJ1(:,j)))

 			elseif (j==n) then

 			TMdwJ0(:,j)=TMdwJ0(:,j-1)*(1.d0+RTMdwJ0(:,j-1))*cdexp(-uhJ0(:,j-1))
 			TMdwJ1(:,j)=TMdwJ1(:,j-1)*(1.d0+RTMdwJ1(:,j-1))*cdexp(-uhJ1(:,j-1))
 			TEdwJ0(:,j)=TEdwJ0(:,j-1)*(1.d0+RTEdwJ0(:,j-1))*cdexp(-uhJ0(:,j-1))
 			TEdwJ1(:,j)=TEdwJ1(:,j-1)*(1.d0+RTEdwJ1(:,j-1))*cdexp(-uhJ1(:,j-1))

 			end if
 		end do
 	elseif (camad < camadT) then
 	allocate(TMupJ0(nJ0,camad:camadT),TMupJ1(nJ1,camad:camadT),TEupJ0(nJ0,camad:camadT),TEupJ1(nJ1,camad:camadT))
 		do j=camadT,camad,-1
 			if (j == camadT) then

 			TMupJ0(:,j) = zeta * mx / (2.d0 * ImpIntJ0(:,camadT))
 			TMupJ1(:,j) = zeta * mx / (2.d0 * ImpIntJ1(:,camadT))
 			TEupJ0(:,j) = zeta * mx / 2.d0
 			TEupJ1(:,j) = zeta * mx / 2.d0

 			elseif (j == (camadT - 1) .and. j == 0) then

 			TMupJ0(:,j)=TMupJ0(:,j+1)*(cdexp(-uJ0(:,camadT)*h0) + &
				RTMupJ0(:,camadT)*AMupJ0(:) + RTMdwJ0(:,camadT)*AMdwJ0(:)* &
				cdexp(-uhJ0(:,camadT)))
 			TMupJ1(:,j)=TMupJ1(:,j+1)*(cdexp(-uJ1(:,camadT)*h0) + &
				RTMupJ1(:,camadT)*AMupJ1(:) + RTMdwJ1(:,camadT)*AMdwJ1(:)* &
				cdexp(-uhJ1(:,camadT)))
 			TEupJ0(:,j)=TEupJ0(:,j+1)*(cdexp(-uJ0(:,camadT)*h0) + &
				RTEupJ0(:,camadT)*FEupJ0(:) - RTEdwJ0(:,camadT)*FEdwJ0(:)* &
				cdexp(-uhJ0(:,camadT)))
 			TEupJ1(:,j)=TEupJ1(:,j+1)*(cdexp(-uJ1(:,camadT)*h0) + &
				RTEupJ1(:,camadT)*FEupJ1(:) - RTEdwJ1(:,camadT)*FEdwJ1(:)* &
				cdexp(-uhJ1(:,camadT)))

 			elseif (j == (camadT - 1) .and. j /= 0) then

 			TMupJ0(:,j)=TMupJ0(:,j+1)*(cdexp(uJ0(:,camadT)*(prof(camadT-1)-h0)) + &
				RTMupJ0(:,camadT)*AMupJ0(:) + RTMdwJ0(:,camadT)*AMdwJ0(:)* &
				cdexp(-uhJ0(:,camadT))) / (1.d0+RTMupJ0(:,j)*cdexp(-2.d0*uhJ0(:,j)))
 			TMupJ1(:,j)=TMupJ1(:,j+1)*(cdexp(uJ1(:,camadT)*(prof(camadT-1)-h0)) + &
				RTMupJ1(:,camadT)*AMupJ1(:) + RTMdwJ1(:,camadT)*AMdwJ1(:)* &
				cdexp(-uhJ1(:,camadT))) / (1.d0+RTMupJ1(:,j)*cdexp(-2.d0*uhJ1(:,j)))
 			TEupJ0(:,j)=TEupJ0(:,j+1)*(cdexp(uJ0(:,camadT)*(prof(camadT-1)-h0)) + &
				RTEupJ0(:,camadT)*FEupJ0(:) - RTEdwJ0(:,camadT)*FEdwJ0(:)* &
				cdexp(-uhJ0(:,camadT))) / (1.d0+RTEupJ0(:,j)*cdexp(-2.d0*uhJ0(:,j)))
 			TEupJ1(:,j)=TEupJ1(:,j+1)*(cdexp(uJ1(:,camadT)*(prof(camadT-1)-h0)) + &
				RTEupJ1(:,camadT)*FEupJ1(:) - RTEdwJ1(:,camadT)*FEdwJ1(:)* &
				cdexp(-uhJ1(:,camadT))) / (1.d0+RTEupJ1(:,j)*cdexp(-2.d0*uhJ1(:,j)))

 			elseif (j /= 0) then

 			TMupJ0(:,j)=TMupJ0(:,j+1)*(1.d0+RTMupJ0(:,j+1))*cdexp(-uhJ0(:,j+1)) / &
				(1.d0 + RTMupJ0(:,j)*cdexp(-2.d0*uhJ0(:,j)))
 			TMupJ1(:,j)=TMupJ1(:,j+1)*(1.d0+RTMupJ1(:,j+1))*cdexp(-uhJ1(:,j+1)) / &
				(1.d0 + RTMupJ1(:,j)*cdexp(-2.d0*uhJ1(:,j)))
 			TEupJ0(:,j)=TEupJ0(:,j+1)*(1.d0+RTEupJ0(:,j+1))*cdexp(-uhJ0(:,j+1)) / &
				(1.d0 + RTEupJ0(:,j)*cdexp(-2.d0*uhJ0(:,j)))
 			TEupJ1(:,j)=TEupJ1(:,j+1)*(1.d0+RTEupJ1(:,j+1))*cdexp(-uhJ1(:,j+1)) / &
				(1.d0 + RTEupJ1(:,j)*cdexp(-2.d0*uhJ1(:,j)))

 			elseif (j == 0) then

 			TMupJ0(:,j)=TMupJ0(:,1)*(1.d0+RTMupJ0(:,1))*cdexp(-uhJ0(:,1))
 			TMupJ1(:,j)=TMupJ1(:,1)*(1.d0+RTMupJ1(:,1))*cdexp(-uhJ1(:,1))
 			TEupJ0(:,j)=TEupJ0(:,1)*(1.d0+RTEupJ0(:,1))*cdexp(-uhJ0(:,1))
 			TEupJ1(:,j)=TEupJ1(:,1)*(1.d0+RTEupJ1(:,1))*cdexp(-uhJ1(:,1))

 			end if
 		end do
 	else
 	allocate(TMdwJ0(nJ0,camadT:camad),TMdwJ1(nJ1,camadT:camad),TEdwJ0(nJ0,camadT:camad),TEdwJ1(nJ1,camadT:camad))
 	allocate(TMupJ0(nJ0,camad:camadT),TMupJ1(nJ1,camad:camadT),TEupJ0(nJ0,camad:camadT),TEupJ1(nJ1,camad:camadT))

 			TMdwJ0(:,camad) = zeta * mx / (2.d0 * ImpIntJ0(:,camadT))
 			TMdwJ1(:,camad) = zeta * mx / (2.d0 * ImpIntJ1(:,camadT))
 			
 			TEdwJ0(:,camad) = - zeta * mx / 2.d0
 			TEdwJ1(:,camad) = - zeta * mx / 2.d0
 			
 			TMupJ0(:,camad) = TMdwJ0(:,camad)
 			TMupJ1(:,camad) = TMdwJ1(:,camad)
 			
 			TEupJ0(:,camad) = - TEdwJ0(:,camad)
 			TEupJ1(:,camad) = - TEdwJ1(:,camad)

 	end if

 	allocate(Ktmdz_J0(nJ0),Ktmdz_J1(nJ1),Ktm_J0(nJ0),Ktm_J1(nJ1),Kte_J0(nJ0),Kte_J1(nJ1),Ktedz_J0(nJ0),Ktedz_J1(nJ1))
 	allocate(kernelExJ0(nJ0),kernelExJ1(nJ1),kernelEyJ0(nJ0),kernelEyJ1(nJ1),kernelEzJ1(nJ1))
 	allocate(kernelHxJ0(nJ0),kernelHxJ1(nJ1),kernelHyJ0(nJ0),kernelHyJ1(nJ1),kernelHzJ1(nJ1))
 	if (camad == 0 .and. camadT /= 0) then

 		Ktmdz_J0=(ImpIntJ0(:,0)*TMupJ0(:,0)*cdexp(uJ0(:,0)*z))*w_J0(:)
 		Ktmdz_J1=(ImpIntJ1(:,0)*TMupJ1(:,0)*cdexp(uJ1(:,0)*z))*w_J1(:)
 		Kte_J0=(TEupJ0(:,0)*cdexp(uJ0(:,0)*z))*w_J0(:)
 		Kte_J1=(TEupJ1(:,0)*cdexp(uJ1(:,0)*z))*w_J1(:)

 		Ktm_J0=(TMupJ0(:,0)*cdexp(uJ0(:,0)*z))*w_J0(:)
 		Ktm_J1=(TMupJ1(:,0)*cdexp(uJ1(:,0)*z))*w_J1(:)
 		Ktedz_J0=(AdmIntJ0(:,0)*TEupJ0(:,0)*cdexp(uJ0(:,0)*z))*w_J0(:)
 		Ktedz_J1=(AdmIntJ1(:,0)*TEupJ1(:,0)*cdexp(uJ1(:,0)*z))*w_J1(:)

 		kernelExJ1 = x * y * (Ktmdz_J1 - Kte_J1) / (r**3)
 		kernelExJ0 = x * y * (Ktmdz_J0 - Kte_J0) * krJ0 / (2.d0 * r**2)
 		Ex_p = (sum(kernelExJ1) - sum(kernelExJ0)) / (pi * r)	!este último r é decorrente do uso dos filtros
 
 		kernelEyJ1 = (1.d0/r - 2.d0*y**2/r**3) * Ktmdz_J1 + (1.d0/r - 2.d0*x**2/r**3) * Kte_J1
 		kernelEyJ0 = (y * y * Ktmdz_J0 + x * x * Kte_J0) * krJ0 / r**2
 		Ey_p = -(sum(kernelEyJ1) + sum(kernelEyJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelEzJ1 = (y / (r*neta)) * Ktm_J1 * krJ1 * krJ1
 		Ez_p = - sum(kernelEzJ1) / (2.d0*pi*r)				!este último r é decorrente do uso dos filtros

 		kernelHxJ1 = (1.d0/r - 2.d0*y**2/r**3) * Ktm_J1 + (1.d0/r - 2.d0*x**2/r**3) * Ktedz_J1
 		kernelHxJ0 = (y * y * Ktm_J0 + x * x * Ktedz_J0) * krJ0 / r**2
 		Hx_p = -(sum(kernelHxJ1) + sum(kernelHxJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros

 		kernelHyJ1 = x * y * (Ktm_J1 - Ktedz_J1) / r**3
 		kernelHyJ0 = x * y * (Ktm_J0 - Ktedz_J0) * krJ0 / (2.d0 * r * r)
 		Hy_p = (-sum(kernelHyJ1) + sum(kernelHyJ0)) / (pi * r)	!este último r é decorrente do uso dos filtros
 
 		kernelHzJ1 = x/(r*zeta) * Kte_J1 * krJ1 * krJ1
 		Hz_p = - sum(kernelHzJ1) / (2.d0*pi*r)				!este último r é decorrente do uso dos filtros
 
 	elseif (camad < camadT) then !camada k

 		ktmdz_J0=(ImpIntJ0(:,camad)*TMupJ0(:,camad)*(cdexp(uJ0(:,camad)*(z-prof(camad))) - &
 			RTMupJ0(:,camad)*cdexp(-uJ0(:,camad)*(z-prof(camad-1)+h(camad)))))*w_J0(:)
 		ktmdz_J1=(ImpIntJ1(:,camad)*TMupJ1(:,camad)*(cdexp(uJ1(:,camad)*(z-prof(camad))) - &
 			RTMupJ1(:,camad)*cdexp(-uJ1(:,camad)*(z-prof(camad-1)+h(camad)))))*w_J1(:)
 		Kte_J0=(TEupJ0(:,camad)*(cdexp(uJ0(:,camad)*(z-prof(camad))) + &
 			 RTEupJ0(:,camad)*cdexp(-uJ0(:,camad)*(z-prof(camad-1)+h(camad)))))*w_J0(:)
 		Kte_J1=(TEupJ1(:,camad)*(cdexp(uJ1(:,camad)*(z-prof(camad))) + &
 			 RTEupJ1(:,camad)*cdexp(-uJ1(:,camad)*(z-prof(camad-1)+h(camad)))))*w_J1(:)

 		ktm_J0=(TMupJ0(:,camad)*(cdexp(uJ0(:,camad)*(z-prof(camad))) + &
 			RTMupJ0(:,camad)*cdexp(-uJ0(:,camad)*(z-prof(camad-1)+h(camad)))))*w_J0(:)
 		ktm_J1=(TMupJ1(:,camad)*(cdexp(uJ1(:,camad)*(z-prof(camad))) + &
 			RTMupJ1(:,camad)*cdexp(-uJ1(:,camad)*(z-prof(camad-1)+h(camad)))))*w_J1(:)
 		Ktedz_J0=(AdmIntJ0(:,camad)*TEupJ0(:,camad)*(cdexp(uJ0(:,camad)*(z-prof(camad))) - &
 			RTEupJ0(:,camad)*cdexp(-uJ0(:,camad)*(z-prof(camad-1)+h(camad)))))*w_J0(:)
 		Ktedz_J1=(AdmIntJ1(:,camad)*TEupJ1(:,camad)*(cdexp(uJ1(:,camad)*(z-prof(camad))) - &
 			RTEupJ1(:,camad)*cdexp(-uJ1(:,camad)*(z-prof(camad-1)+h(camad)))))*w_J1(:)
 
 		kernelExJ1 = x * y * (Ktmdz_J1 - Kte_J1) / (r * r * r)
 		kernelExJ0 = x * y * (Ktmdz_J0 - Kte_J0) * krJ0 / (2.d0 * r * r)
 		Ex_p = (sum(kernelExJ1) - sum(kernelExJ0)) / (pi * r)	!este último r é decorrente do uso dos filtros
 
 		kernelEyJ1 = (1.d0/r - 2.d0*y**2/r**3) * Ktmdz_J1 + (1.d0/r - 2.d0*x**2/r**3) * Kte_J1
 		kernelEyJ0 = (y * y * Ktmdz_J0 + x * x * Kte_J0) * krJ0 / r**2
 		Ey_p = -(sum(kernelEyJ1) + sum(kernelEyJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelEzJ1 = (y/(r*condut(camad))) * Ktm_J1 * krJ1 * krJ1
 		Ez_p = -sum(kernelEzJ1) / (2.d0*pi*r)				!este último r é decorrente do uso dos filtros
 
 		kernelHxJ1 = (1.d0/r - 2.d0*y**2/r**3) * Ktm_J1 + (1.d0/r - 2.d0*x**2/r**3) * Ktedz_J1
 		kernelHxJ0 = (y * y * Ktm_J0 + x * x * Ktedz_J0) * krJ0 / r**2
 		Hx_p = -(sum(kernelHxJ1) + sum(kernelHxJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelHyJ1 = x * y * (Ktm_J1 - Ktedz_J1) / r**3
 		kernelHyJ0 = x * y * (Ktm_J0 - Ktedz_J0) * krJ0 / (2.d0 * r * r)
 		Hy_p = (-sum(kernelHyJ1) + sum(kernelHyJ0)) / (pi * r)	!este último r é decorrente do uso dos filtros
 
 		kernelHzJ1 = x/(r*zeta) * Kte_J1 * krJ1 * krJ1
 		Hz_p = - sum(kernelHzJ1) / (2.d0*pi*r)				!este último r é decorrente do uso dos filtros
 
 	elseif (camad == camadT .and. z <= h0) then	!na mesma camada do transmissor mas acima dele

		Ktmdz_J0=(ImpIntJ0(:,camad)*TMupJ0(:,camad)*(cdexp(uJ0(:,camad)*(z-h0)) - &
 			RTMupJ0(:,camad)*AMupJ0(:)*cdexp(-uJ0(:,camad)*(z-prof(camad-1))) + &
 			RTMdwJ0(:,camad)*AMdwJ0(:)*cdexp(uJ0(:,camad)*(z-prof(camad)))))*w_J0(:)
 		Ktmdz_J1=(ImpIntJ1(:,camad)*TMupJ1(:,camad)*(cdexp(uJ1(:,camad)*(z-h0)) - &
 			RTMupJ1(:,camad)*AMupJ1(:)*cdexp(-uJ1(:,camad)*(z-prof(camad-1))) + &
 			RTMdwJ1(:,camad)*AMdwJ1(:)*cdexp(uJ1(:,camad)*(z-prof(camad)))))*w_J1(:)
 		Kte_J0=(TEupJ0(:,camad)*(cdexp(uJ0(:,camad)*(z-h0)) + &
 			RTEupJ0(:,camad)*FEupJ0(:)*cdexp(-uJ0(:,camad)*(z-prof(camad-1))) - &
 			RTEdwJ0(:,camad)*FEdwJ0(:)*cdexp(uJ0(:,camad)*(z-prof(camad)))))*w_J0(:)
 		Kte_J1=(TEupJ1(:,camad)*(cdexp(uJ1(:,camad)*(z-h0)) + &
 			RTEupJ1(:,camad)*FEupJ1(:)*cdexp(-uJ1(:,camad)*(z-prof(camad-1))) - &
 			RTEdwJ1(:,camad)*FEdwJ1(:)*cdexp(uJ1(:,camad)*(z-prof(camad)))))*w_J1(:)

 		Ktm_J0=(TMupJ0(:,camad)*(cdexp(uJ0(:,camad)*(z-h0)) + &
 			RTMupJ0(:,camad)*AMupJ0(:)*cdexp(-uJ0(:,camad)*(z-prof(camad-1))) + &
 			RTMdwJ0(:,camad)*AMdwJ0(:)*cdexp(uJ0(:,camad)*(z-prof(camad)))))*w_J0(:)
 		Ktm_J1=(TMupJ1(:,camad)*(cdexp(uJ1(:,camad)*(z-h0)) + &
 			RTMupJ1(:,camad)*AMupJ1(:)*cdexp(-uJ1(:,camad)*(z-prof(camad-1))) + &
 			RTMdwJ1(:,camad)*AMdwJ1(:)*cdexp(uJ1(:,camad)*(z-prof(camad)))))*w_J1(:)
 		Ktedz_J0=(AdmIntJ0(:,camad)*TEupJ0(:,camad)*(cdexp(uJ0(:,camad)*(z-h0)) - &
 			RTEupJ0(:,camad)*FEupJ0(:)*cdexp(-uJ0(:,camad)*(z-prof(camad-1))) - &
 			RTEdwJ0(:,camad)*FEdwJ0(:)*cdexp(uJ0(:,camad)*(z-prof(camad)))))*w_J0(:)
 		Ktedz_J1=(AdmIntJ1(:,camad)*TEupJ1(:,camad)*(cdexp(uJ1(:,camad)*(z-h0)) - &
 			RTEupJ1(:,camad)*FEupJ1(:)*cdexp(-uJ1(:,camad)*(z-prof(camad-1))) - &
 			RTEdwJ1(:,camad)*FEdwJ1(:)*cdexp(uJ1(:,camad)*(z-prof(camad)))))*w_J1(:)

 		kernelExJ1 = x * y * (Ktmdz_J1 - Kte_J1) / (r * r * r)
 		kernelExJ0 = x * y * (Ktmdz_J0 - Kte_J0) * krJ0 / (2.d0 * r * r)
 		Ex_p = (sum(kernelExJ1) - sum(kernelExJ0)) / (pi * r)	!este último r é decorrente do uso dos filtros
 
 		kernelEyJ1 = (1.d0/r - 2.d0*y**2/r**3) * Ktmdz_J1 + (1.d0/r - 2.d0*x**2/r**3) * Kte_J1
 		kernelEyJ0 = (y * y * Ktmdz_J0 + x * x * Kte_J0) * krJ0 / r**2
 		Ey_p = -(sum(kernelEyJ1) + sum(kernelEyJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
		if (camad /= 0) then
 		kernelEzJ1 = (y/(r*condut(camad))) * Ktm_J1 * krJ1 * krJ1
		else
 		kernelEzJ1 = (y/(r*neta)) * Ktm_J1 * krJ1 * krJ1
		end if
 		Ez_p = -sum(kernelEzJ1) / (2.d0*pi*r)				!este último r é decorrente do uso dos filtros
 
 		kernelHxJ1 = (1.d0/r - 2.d0*y**2/r**3) * Ktm_J1 + (1.d0/r - 2.d0*x**2/r**3) * Ktedz_J1
 		kernelHxJ0 = (y * y * Ktm_J0 + x * x * Ktedz_J0) * krJ0 / r**2
 		Hx_p = -(sum(kernelHxJ1) + sum(kernelHxJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelHyJ1 = x * y * (Ktm_J1 - Ktedz_J1) / r**3
 		kernelHyJ0 = x * y * (Ktm_J0 - Ktedz_J0) * krJ0 / (2.d0 * r * r)
 		Hy_p = (-sum(kernelHyJ1) + sum(kernelHyJ0)) / (pi * r)	!este último r é decorrente do uso dos filtros
 
 		kernelHzJ1 = x/(r*zeta) * Kte_J1 * krJ1 * krJ1
 		Hz_p = - sum(kernelHzJ1) / (2.d0*pi*r)				!este último r é decorrente do uso dos filtros
 
 	elseif (camad == camadT .and. z > h0) then	!na mesma camada do transmissor mas abaixo dele

 		Ktmdz_J0=(ImpIntJ0(:,camad)*TMdwJ0(:,camad)*(cdexp(-uJ0(:,camad)*(z-h0)) + &
 			RTMupJ0(:,camad)*AMupJ0(:)*cdexp(-uJ0(:,camad)*(z-prof(camad-1))) - &
 			RTMdwJ0(:,camad)*AMdwJ0(:)*cdexp(uJ0(:,camad)*(z-prof(camad)))))*w_J0(:)
 		Ktmdz_J1=(ImpIntJ1(:,camad)*TMdwJ1(:,camad)*(cdexp(-uJ1(:,camad)*(z-h0)) + &
 			RTMupJ1(:,camad)*AMupJ1(:)*cdexp(-uJ1(:,camad)*(z-prof(camad-1))) - &
 			RTMdwJ1(:,camad)*AMdwJ1(:)*cdexp(uJ1(:,camad)*(z-prof(camad)))))*w_J1(:)

 		Kte_J0=(TEdwJ0(:,camad)*(cdexp(-uJ0(:,camad)*(z-h0)) - &
 			RTEupJ0(:,camad)*FEupJ0(:)*cdexp(-uJ0(:,camad)*(z-prof(camad-1))) + &
 			RTEdwJ0(:,camad)*FEdwJ0(:)*cdexp(uJ0(:,camad)*(z-prof(camad)))))*w_J0(:)
 		Kte_J1=(TEdwJ1(:,camad)*(cdexp(-uJ1(:,camad)*(z-h0)) - &
 			RTEupJ1(:,camad)*FEupJ1(:)*cdexp(-uJ1(:,camad)*(z-prof(camad-1))) + &
 			RTEdwJ1(:,camad)*FEdwJ1(:)*cdexp(uJ1(:,camad)*(z-prof(camad)))))*w_J1(:)

		Ktm_J0=(TMdwJ0(:,camad)*(cdexp(-uJ0(:,camad)*(z-h0)) + &
 			RTMupJ0(:,camad)*AMupJ0(:)*cdexp(-uJ0(:,camad)*(z-prof(camad-1))) + &
 			RTMdwJ0(:,camad)*AMdwJ0(:)*cdexp(uJ0(:,camad)*(z-prof(camad)))))*w_J0(:)
 		Ktm_J1=(TMdwJ1(:,camad)*(cdexp(-uJ1(:,camad)*(z-h0)) + &
 			RTMupJ1(:,camad)*AMupJ1(:)*cdexp(-uJ1(:,camad)*(z-prof(camad-1))) + &
 			RTMdwJ1(:,camad)*AMdwJ1(:)*cdexp(uJ1(:,camad)*(z-prof(camad)))))*w_J1(:)

 		Ktedz_J0=(AdmIntJ0(:,camad)*TEdwJ0(:,camad)*(cdexp(-uJ0(:,camad)*(z-h0)) - &
 			RTEupJ0(:,camad)*FEupJ0(:)*cdexp(-uJ0(:,camad)*(z-prof(camad-1))) - &
 			RTEdwJ0(:,camad)*FEdwJ0(:)*cdexp(uJ0(:,camad)*(z-prof(camad)))))*w_J0(:)
 		Ktedz_J1=(AdmIntJ1(:,camad)*TEdwJ1(:,camad)*(cdexp(-uJ1(:,camad)*(z-h0)) - &
 			RTEupJ1(:,camad)*FEupJ1(:)*cdexp(-uJ1(:,camad)*(z-prof(camad-1))) - &
 			RTEdwJ1(:,camad)*FEdwJ1(:)*cdexp(uJ1(:,camad)*(z-prof(camad)))))*w_J1(:)
 
 		kernelExJ1 = x * y * (Ktmdz_J1 + Kte_J1) / (r * r * r)
 		kernelExJ0 = x * y * (Ktmdz_J0 + Kte_J0) * krJ0 / (2.d0 * r * r)
 		Ex_p = (-sum(kernelExJ1) + sum(kernelExJ0)) / (pi * r)	!este último r é decorrente do uso dos filtros
 
 		kernelEyJ1 = (1.d0/r - 2.d0*y**2/r**3) * Ktmdz_J1 - (1.d0/r - 2.d0*x**2/r**3) * Kte_J1
 		kernelEyJ0 = (y * y * Ktmdz_J0 - x * x * Kte_J0) * krJ0 / r**2
 		Ey_p = (sum(kernelEyJ1) + sum(kernelEyJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros

		if (camad /= 0) then
 		kernelEzJ1 = (y/(r*condut(camad))) * Ktm_J1 * krJ1 * krJ1
		else
 		kernelEzJ1 = (y/(r*neta)) * Ktm_J1 * krJ1 * krJ1
		end if
 		Ez_p = -sum(kernelEzJ1) / (2.d0*pi*r)				!este último r é decorrente do uso dos filtros
 
 		kernelHxJ1 = (1.d0/r - 2.d0*y**2/r**3) * Ktm_J1 - (1.d0/r - 2.d0*x**2/r**3) * Ktedz_J1
 		kernelHxJ0 = (y * y * Ktm_J0 - x * x * Ktedz_J0) * krJ0 / r**2
 		Hx_p = -(sum(kernelHxJ1) + sum(kernelHxJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelHyJ1 = x * y * (Ktm_J1 + Ktedz_J1) / r**3
 		kernelHyJ0 = x * y * (Ktm_J0 + Ktedz_J0) * krJ0 / (2.d0 * r * r)
 		Hy_p = (-sum(kernelHyJ1) + sum(kernelHyJ0)) / (pi * r)	!este último r é decorrente do uso dos filtros
 
 		kernelHzJ1 = x/(r*zeta) * Kte_J1 * krJ1 * krJ1
 		Hz_p = - sum(kernelHzJ1) / (2.d0*pi*r)				!este último r é decorrente do uso dos filtros
 
 	elseif (camad > camadT .and. camad /= n) then !camada j

 		Ktmdz_J0=(ImpIntJ0(:,camad)*TMdwJ0(:,camad)*(cdexp(-uJ0(:,camad)*(z-prof(camad-1))) - &
 			RTMdwJ0(:,camad)*cdexp(uJ0(:,camad)*(z-prof(camad)-h(camad)))))*w_J0(:)
 		Ktmdz_J1=(ImpIntJ1(:,camad)*TMdwJ1(:,camad)*(cdexp(-uJ1(:,camad)*(z-prof(camad-1))) - &
 			RTMdwJ1(:,camad)*cdexp(uJ1(:,camad)*(z-prof(camad)-h(camad)))))*w_J1(:)

 		Kte_J0=(TEdwJ0(:,camad)*(cdexp(-uJ0(:,camad)*(z-prof(camad-1))) + &
 			RTEdwJ0(:,camad)*cdexp(uJ0(:,camad)*(z-prof(camad)-h(camad)))))*w_J0(:)
 		Kte_J1=(TEdwJ1(:,camad)*(cdexp(-uJ1(:,camad)*(z-prof(camad-1))) + &
 			RTEdwJ1(:,camad)*cdexp(uJ1(:,camad)*(z-prof(camad)-h(camad)))))*w_J1(:)

 		Ktm_J0=(TMdwJ0(:,camad)*(cdexp(-uJ0(:,camad)*(z-prof(camad-1))) + &
 			RTMdwJ0(:,camad)*cdexp(uJ0(:,camad)*(z-prof(camad)-h(camad)))))*w_J0(:)
 		Ktm_J1=(TMdwJ1(:,camad)*(cdexp(-uJ1(:,camad)*(z-prof(camad-1))) + &
 			RTMdwJ1(:,camad)*cdexp(uJ1(:,camad)*(z-prof(camad)-h(camad)))))*w_J1(:)

 		Ktedz_J0=(AdmIntJ0(:,camad)*TEdwJ0(:,camad)*(cdexp(-uJ0(:,camad)*(z-prof(camad-1))) - &
 			RTEdwJ0(:,camad)*cdexp(uJ0(:,camad)*(z-prof(camad)-h(camad)))))*w_J0(:)
 		Ktedz_J1=(AdmIntJ1(:,camad)*TEdwJ1(:,camad)*(cdexp(-uJ1(:,camad)*(z-prof(camad-1))) - &
 			RTEdwJ1(:,camad)*cdexp(uJ1(:,camad)*(z-prof(camad)-h(camad)))))*w_J1(:)
 
 		kernelExJ1 = x * y * (Ktmdz_J1 + Kte_J1) / (r * r * r)
 		kernelExJ0 = x * y * (Ktmdz_J0 + Kte_J0) * krJ0 / (2.d0 * r * r)
 		Ex_p = (-sum(kernelExJ1) + sum(kernelExJ0)) / (pi * r)	!este último r é decorrente do uso dos filtros
 
 		kernelEyJ1 = (1.d0/r - 2.d0*y**2/r**3) * Ktmdz_J1 - (1.d0/r - 2.d0*x**2/r**3) * Kte_J1
 		kernelEyJ0 = (y * y * Ktmdz_J0 - x * x * Kte_J0) * krJ0 / r**2
 		Ey_p = (sum(kernelEyJ1) + sum(kernelEyJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelEzJ1 = (y/(r*condut(camad))) * Ktm_J1 * krJ1 * krJ1
 		Ez_p = -sum(kernelEzJ1) / (2.d0*pi*r)				!este último r é decorrente do uso dos filtros
 
 		kernelHxJ1 = (1.d0/r - 2.d0*y**2/r**3) * Ktm_J1 - (1.d0/r - 2.d0*x**2/r**3) * Ktedz_J1
 		kernelHxJ0 = (y * y * Ktm_J0 - x * x * Ktedz_J0) * krJ0 / r**2
 		Hx_p = -(sum(kernelHxJ1) + sum(kernelHxJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelHyJ1 = x * y * (Ktm_J1 + Ktedz_J1) / r**3
 		kernelHyJ0 = x * y * (Ktm_J0 + Ktedz_J0) * krJ0 / (2.d0 * r * r)
 		Hy_p = (-sum(kernelHyJ1) + sum(kernelHyJ0)) / (pi * r)	!este último r é decorrente do uso dos filtros
 
 		kernelHzJ1 = x/(r*zeta) * Kte_J1 * krJ1 * krJ1
 		Hz_p = - sum(kernelHzJ1) / (2.d0*pi*r)				!este último r é decorrente do uso dos filtros

 	else	!camada n

 		Ktmdz_J0=(ImpIntJ0(:,n)*TMdwJ0(:,n)*cdexp(-uJ0(:,n)*(z-prof(n-1))))*w_J0(:)
 		Ktmdz_J1=(ImpIntJ1(:,n)*TMdwJ1(:,n)*cdexp(-uJ1(:,n)*(z-prof(n-1))))*w_J1(:)

 		Kte_J0=(TEdwJ0(:,n)*cdexp(-uJ0(:,n)*(z-prof(n-1))))*w_J0(:)
 		Kte_J1=(TEdwJ1(:,n)*cdexp(-uJ1(:,n)*(z-prof(n-1))))*w_J1(:)

 		Ktm_J0=(TMdwJ0(:,n)*cdexp(-uJ0(:,n)*(z-prof(n-1))))*w_J0(:)
 		Ktm_J1=(TMdwJ1(:,n)*cdexp(-uJ1(:,n)*(z-prof(n-1))))*w_J1(:)

 		Ktedz_J0=(AdmIntJ0(:,n)*TEdwJ0(:,n)*cdexp(-uJ0(:,n)*(z-prof(n-1))))*w_J0(:)
 		Ktedz_J1=(AdmIntJ1(:,n)*TEdwJ1(:,n)*cdexp(-uJ1(:,n)*(z-prof(n-1))))*w_J1(:)
 
 		kernelExJ1 = x * y * (Ktmdz_J1 + Kte_J1) / (r * r * r)
 		kernelExJ0 = x * y * (Ktmdz_J0 + Kte_J0) * krJ0 / (2.d0 * r * r)
 		Ex_p = (-sum(kernelExJ1) + sum(kernelExJ0)) / (pi * r)	!este último r é decorrente do uso dos filtros
 
 		kernelEyJ1 = (1.d0/r - 2.d0*y**2/r**3) * Ktmdz_J1 - (1.d0/r - 2.d0*x**2/r**3) * Kte_J1
 		kernelEyJ0 = (y * y * Ktmdz_J0 - x * x * Kte_J0) * krJ0 / r**2
 		Ey_p = (sum(kernelEyJ1) + sum(kernelEyJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelEzJ1 = (y/(r*condut(n))) * Ktm_J1 * krJ1 * krJ1
 		Ez_p = -sum(kernelEzJ1) / (2.d0*pi*r)				!este último r é decorrente do uso dos filtros
 
 		kernelHxJ1 = (1.d0/r - 2.d0*y**2/r**3) * Ktm_J1 - (1.d0/r - 2.d0*x**2/r**3) * Ktedz_J1
 		kernelHxJ0 = (y * y * Ktm_J0 - x * x * Ktedz_J0) * krJ0 / r**2
 		Hx_p = -(sum(kernelHxJ1) + sum(kernelHxJ0)) / (2.d0*pi*r)	!este último r é decorrente do uso dos filtros
 
 		kernelHyJ1 = x * y * (Ktm_J1 + Ktedz_J1) / r**3
 		kernelHyJ0 = x * y * (Ktm_J0 + Ktedz_J0) * krJ0 / (2.d0 * r * r)
 		Hy_p = (-sum(kernelHyJ1) + sum(kernelHyJ0)) / (pi * r)	!este último r é decorrente do uso dos filtros
 
 		kernelHzJ1 = x/(r*zeta) * Kte_J1 * krJ1 * krJ1
 		Hz_p = - sum(kernelHzJ1) / (2.d0*pi*r)				!este último r é decorrente do uso dos filtros
 
 	end if
 
 	deallocate(h,KrJ0,KrJ1,w_J0,w_J1)
 	deallocate(wvnb2,uJ0,uJ1,AdmIntJ0,AdmIntJ1,ImpIntJ0,ImpIntJ1,uhJ0,uhJ1,tghJ0,tghJ1)
 	deallocate(AdmApdwJ0,AdmApdwJ1,ImpApdwJ0,ImpApdwJ1,RTEdwJ0,RTEdwJ1,RTMdwJ0,RTMdwJ1)
 	deallocate(AdmApupJ0,AdmApupJ1,ImpApupJ0,ImpApupJ1,RTEupJ0,RTEupJ1,RTMupJ0,RTMupJ1)
 	deallocate(AMdwJ0,AMdwJ1,AMupJ0,AMupJ1,FEdwJ0,FEdwJ1,FEupJ0,FEupJ1)
 	deallocate(Ktmdz_J0,Ktmdz_J1,Ktm_J0,Ktm_J1,Kte_J0,Kte_J1,Ktedz_J0,Ktedz_J1)
 	deallocate(kernelExJ0,kernelExJ1,kernelEyJ0,kernelEyJ1,kernelEzJ1)
 	deallocate(kernelHxJ0,kernelHxJ1,kernelHyJ0,kernelHyJ1,kernelHzJ1)
	end subroutine dmhx_xyz_loops
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
subroutine dmv_xyz_loops( idtfcd_cJ0, nJ0, idtfcd_cJ1, nJ1, Tx, Ty, mz, h0, n, esp, condut, &
				neta, zeta, cx, cy, z, Ex_p, Ey_p, Ez_p, Hx_p, Hy_p, Hz_p )
	implicit none
	integer(4), intent(in) :: n, idtfcd_cJ0, nJ0, idtfcd_cJ1, nJ1
	real(8), intent(in) :: Tx, Ty, mz, h0, esp(:), condut(1:n), cx, cy, z
	complex*16, intent(in) :: zeta, neta
	complex*16, intent(out) :: Ex_p, Ey_p, Ez_p, Hx_p, Hy_p, Hz_p

	real(8), parameter :: pi = 3.141592653589793238462643383279502884197d0
	integer(4) :: i, j, k, camad, camadT, filtro, ident_fJ0, ident_fJ1
	real(8) :: x, y, r, kr, k_r
	real(8), dimension(:), allocatable :: h, krJ0, krJ1, w_J0, w_J1, prof

!	Para uso de loops:
	complex*16, dimension(:), allocatable :: wvnb2, FEdwJ0, FEdwJ1, FEupJ0, FEupJ1
	complex*16, dimension(:,:), allocatable :: uJ0, uJ1, AdmIntJ0, AdmIntJ1, uhJ0, uhJ1, tghJ0, tghJ1
	complex*16, dimension(:,:), allocatable :: AdmApdwJ0, AdmApdwJ1, AdmApupJ0, AdmApupJ1
	complex*16, dimension(:,:), allocatable :: RTEdwJ0, RTEdwJ1, RTEupJ0, RTEupJ1
	complex*16, dimension(:,:), allocatable :: TEdwJ0, TEdwJ1, TEupJ0, TEupJ1

	complex*16, dimension(:), allocatable :: Kte_J0, Kte_J1, Ktedz_J1
	complex*16, dimension(:), allocatable :: kernelExJ0, kernelExJ1, kernelEyJ0, kernelEyJ1
	complex*16, dimension(:), allocatable :: kernelHxJ0, kernelHxJ1, kernelHyJ0, kernelHyJ1, kernelHzJ1

	complex*16 :: kerEx_J1, kerEx_J0, kerEy_J1, kerEy_J0
	complex*16 :: kerHx_J1, kerHx_J0, kerHy_J1, kerHy_J0, kerHz_J1

	Ez_p = (0.d0,0.d0)

	if ( cx == Tx ) then
        x = dsign( 1.d0, Tx )
	else
        x = cx - Tx
	end if

    y = cy - Ty
	r = dsqrt( x ** 2 + y ** 2 )

	allocate( h(0 : n), prof(-1 : n) )
	if ( size(esp) == n ) then
		h(0) = 0.d0
		h(1 : n) = esp
	else
		h(0) = 0.d0
		h(1 : n - 1) = esp
		h(n) = 1.d300
	end if
!	criando um novo vetor de profundidades que se adeque à qualquer situação patológica
	prof(-1) = -1.d300
	prof(0) = 0.d0
	if ( n > 1 ) then
        prof(1) = h(1)
        if ( n > 2 ) then
            do k = 2, n - 1
                prof(k) = prof(k - 1) + h(k)
            end do
        end if
	end if
	prof(n) = 1.d300

!para descobrir em que camada está a observação
	if ( z <= 0.d0 ) then
		camad = 0
	else if ( z > prof(n - 1) ) then
		camad = n
	else
        do i = n - 1, 1, -1
            if ( z > prof(i - 1) ) then
                camad = i
                exit
            end if
        end do
	end if

!para descobrir em que camada está o transmissor
	if ( h0 <= 0.d0 ) then
        camadT = 0
	else if ( h0 > prof(n - 1) ) then
		camadT = n
	else
        do j = n - 1, 1, -1
            if ( h0 > prof(j - 1) ) then
                camadT = j
                exit
            end if
        end do
	end if

!!	write(*,*)'Entre com o criador dos filtros J0: Rijo(0), Frayzer(1), Guptasarma(2), Kong(3) ou Key(4)'
!!	read(*,*)idtfcd_cJ0
	filtro = 0	!esta variável direciona o uso de filtros J0 e J1 em vez de seno e cosseno
!!	call identfiltro(filtro,idtfcd_cJ0,ident_fJ0,nJ0)
!!	write(*,*)'Entre com o criador dos filtros J1: Rijo(0), Frayzer(1), Guptasarma(2), Kong(3) ou Key(4)'
!!	read(*,*)idtfcd_cJ1
!!	call identfiltro(filtro,idtfcd_cJ1,ident_fJ1,nJ1)

!	idtfcd_cJ0 = 3
	ident_fJ0 = 0
!	nJ0 = 241
!	idtfcd_cJ1 = 3
	ident_fJ1 = 1
!	nJ1 = 241

	allocate( KrJ0(nJ0), KrJ1(nJ1), w_J0(nJ0), w_J1(nJ1) )

	call constfiltro( filtro, idtfcd_cJ0, ident_fJ0, nJ0, r, KrJ0, w_J0 )
	call constfiltro( filtro, idtfcd_cJ1, ident_fJ1, nJ1, r, KrJ1, w_J1 )

 	allocate( wvnb2(0 : n), uJ0(nJ0, 0 : n), uJ1(nJ1, 0 : n), AdmIntJ0(nJ0, 0 : n), AdmIntJ1(nJ1, 0 : n) )
 	allocate( uhJ0(nJ0, 0 : n), uhJ1(nJ1, 0 : n), tghJ0(nJ0, 0 : n), tghJ1(nJ1, 0 : n) )
! 
 	do i = 0, n
        if ( i == 0 ) then
            wvnb2(i) = -zeta * neta
            uJ0(:,i) = cdsqrt( krJ0 * krJ0 - wvnb2(i) )
            uJ1(:,i) = cdsqrt( krJ1 * krJ1 - wvnb2(i) )
            AdmIntJ0(:,i) = uJ0(:,i) / zeta
            AdmIntJ1(:,i) = uJ1(:,i) / zeta
            uhJ0(:,i) = uJ0(:,i) * h(i)
            uhJ1(:,i) = uJ1(:,i) * h(i)
            tghJ0(:,i) = ( 1.d0 - cdexp( -2.d0 * uhJ0(:,i) ) ) / ( 1.d0 + cdexp( -2.d0 * uhJ0(:,i) ) )
            tghJ1(:,i) = ( 1.d0 - cdexp( -2.d0 * uhJ1(:,i) ) ) / ( 1.d0 + cdexp( -2.d0 * uhJ1(:,i) ) )
 		else
            wvnb2(i) = -zeta * condut(i)
            uJ0(:,i) = cdsqrt( krJ0 * krJ0 - wvnb2(i) )
            uJ1(:,i) = cdsqrt( krJ1 * krJ1 - wvnb2(i) )
            AdmIntJ0(:,i) = uJ0(:,i) / zeta
            AdmIntJ1(:,i) = uJ1(:,i) / zeta
            uhJ0(:,i) = uJ0(:,i) * h(i)
            uhJ1(:,i) = uJ1(:,i) * h(i)
            tghJ0(:,i) = ( 1.d0 - cdexp( -2.d0 * uhJ0(:,i) ) ) / ( 1.d0 + cdexp( -2.d0 * uhJ0(:,i) ) )
            tghJ1(:,i) = ( 1.d0 - cdexp( -2.d0 * uhJ1(:,i) ) ) / ( 1.d0 + cdexp( -2.d0 * uhJ1(:,i) ) )
 		end if
 	end do
 
 	allocate( AdmApdwJ0(nJ0, 1 : n), AdmApdwJ1(nJ1, 1 : n) )
 	allocate( RTEdwJ0(nJ0, 0 : n), RTEdwJ1(nJ1, 0 : n) )
 
 	do i = n, 1, -1
        if ( i == n ) then
            AdmApdwJ0(:,i) = AdmIntJ0(:,i)
            AdmApdwJ1(:,i) = AdmIntJ1(:,i)
            RTEdwJ0(:,i) = (0.d0, 0.d0)
            RTEdwJ1(:,i) = (0.d0, 0.d0)
 		else
            AdmApdwJ0(:,i) = AdmIntJ0(:,i) * ( AdmApdwJ0(:,i + 1) + AdmIntJ0(:,i) * &
            	tghJ0(:,i) ) / ( AdmIntJ0(:,i) + AdmApdwJ0(:,i + 1) * tghJ0(:,i) )
            AdmApdwJ1(:,i) = AdmIntJ1(:,i) * ( AdmApdwJ1(:,i + 1) + AdmIntJ1(:,i) * &
            	tghJ1(:,i) ) / ( AdmIntJ1(:,i) + AdmApdwJ1(:,i + 1) * tghJ1(:,i) )
            RTEdwJ0(:,i) = ( AdmIntJ0(:,i) - AdmApdwJ0(:,i + 1) ) / ( AdmIntJ0(:,i) + AdmApdwJ0(:,i + 1) )
            RTEdwJ1(:,i) = ( AdmIntJ1(:,i) - AdmApdwJ1(:,i + 1) ) / ( AdmIntJ1(:,i) + AdmApdwJ1(:,i + 1) )
 		end if	
 	end do
    RTEdwJ0(:,0) = ( AdmIntJ0(:,0) - AdmApdwJ0(:,1) ) / ( AdmIntJ0(:,0) + AdmApdwJ0(:,1) )
    RTEdwJ1(:,0) = ( AdmIntJ1(:,0) - AdmApdwJ1(:,1) ) / ( AdmIntJ1(:,0) + AdmApdwJ1(:,1) )
 
 	allocate( AdmApupJ0(nJ0,0 : n - 1), AdmApupJ1(nJ1, 0 : n - 1) )
 	allocate( RTEupJ0(nJ0, 0 : n), RTEupJ1(nJ1,0 : n) )
 
 	do i = 0, n - 1
        if ( i == 0 ) then
            AdmApupJ0(:,i) = AdmIntJ0(:,i)
            AdmApupJ1(:,i) = AdmIntJ1(:,i)
            RTEupJ0(:,i) = (0.d0,0.d0)
            RTEupJ1(:,i) = (0.d0,0.d0)
 		else
            AdmApupJ0(:,i) = AdmIntJ0(:,i) * ( AdmApupJ0(:,i - 1) + AdmIntJ0(:,i) * &
            	tghJ0(:,i) ) / ( AdmIntJ0(:,i) + AdmApupJ0(:,i - 1) * tghJ0(:,i) )
            AdmApupJ1(:,i) = AdmIntJ1(:,i) * ( AdmApupJ1(:,i - 1) + AdmIntJ1(:,i) * &
            	tghJ1(:,i) ) / ( AdmIntJ1(:,i) + AdmApupJ1(:,i - 1) * tghJ1(:,i) )
            RTEupJ0(:,i) = ( AdmIntJ0(:,i) - AdmApupJ0(:,i - 1) ) / ( AdmIntJ0(:,i) + AdmApupJ0(:,i - 1) )
            RTEupJ1(:,i) = ( AdmIntJ1(:,i) - AdmApupJ1(:,i - 1) ) / ( AdmIntJ1(:,i) + AdmApupJ1(:,i - 1) )
 		end if
 	end do
    RTEupJ0(:,n) = ( AdmIntJ0(:,n) - AdmApupJ0(:,n - 1) ) / ( AdmIntJ0(:,n) + AdmApupJ0(:,n - 1) )
    RTEupJ1(:,n) = ( AdmIntJ1(:,n) - AdmApupJ1(:,n - 1) ) / ( AdmIntJ1(:,n) + AdmApupJ1(:,n - 1) )

 	allocate( FEdwJ0(nJ0), FEdwJ1(nJ1), FEupJ0(nJ0), FEupJ1(nJ1) )

 	FEdwJ0 = ( cdexp( -uJ0(:,camadT) * ( prof(camadT) - h0) ) + RTEupJ0(:,camadT) * &
            cdexp( uJ0(:,camadT) * ( prof(camadT - 1) - h(camadT) - h0 ) ) ) / &
            ( 1.d0 - RTEupJ0(:,camadT) * RTEdwJ0(:,camadT) * cdexp( -2.d0 * uhJ0(:,camadT) ) )
 	FEdwJ1 = ( cdexp( -uJ1(:,camadT) * ( prof(camadT) - h0) ) + RTEupJ1(:,camadT) * &
            cdexp( uJ1(:,camadT) * ( prof(camadT - 1) - h(camadT) - h0 ) ) ) / &
            ( 1.d0 - RTEupJ1(:,camadT) * RTEdwJ1(:,camadT) * cdexp( -2.d0 * uhJ1(:,camadT) ) )

 	FEupJ0 = ( cdexp( uJ0(:,camadT) * ( prof(camadT - 1) - h0 ) ) + RTEdwJ0(:,camadT) * &
            cdexp( -uJ0(:,camadT) * ( prof(camadT) + h(camadT) - h0 ) ) ) / &
            ( 1.d0 - RTEupJ0(:,camadT) * RTEdwJ0(:,camadT) * cdexp( -2.d0 * uhJ0(:,camadT) ) )
 	FEupJ1 = ( cdexp( uJ1(:,camadT) * ( prof(camadT - 1) - h0 ) ) + RTEdwJ1(:,camadT) * &
            cdexp( -uJ1(:,camadT) * ( prof(camadT) + h(camadT) - h0 ) ) ) / &
            ( 1.d0 - RTEupJ1(:,camadT) * RTEdwJ1(:,camadT) * cdexp( -2.d0 * uhJ1(:,camadT) ) )

 	if ( camad > camadT ) then
        allocate( TEdwJ0(nJ0,camadT : camad), TEdwJ1(nJ1,camadT : camad) )
        do j = camadT, camad
            if ( j == camadT ) then

                TEdwJ0(:,j)= zeta * mz / ( 2.d0 * uJ0(:,j) )
                TEdwJ1(:,j)= zeta * mz / ( 2.d0 * uJ1(:,j) )

 			else if ( j == (camadT + 1) .and. j == n ) then

                TEdwJ0(:,j) = TEdwJ0(:,j - 1) * ( cdexp( -uJ0(:,camadT) * ( prof(camadT) - h0 ) ) + &
                        RTEupJ0(:,camadT) * FEupJ0(:) * cdexp( -uhJ0(:,camadT) ) + RTEdwJ0(:,camadT) * FEdwJ0(:) )
                TEdwJ1(:,j) = TEdwJ1(:,j - 1) * ( cdexp( -uJ1(:,camadT) * ( prof(camadT) - h0) ) + &
                        RTEupJ1(:,camadT) * FEupJ1(:) * cdexp( -uhJ1(:,camadT) ) + RTEdwJ1(:,camadT) * FEdwJ1(:) )

 			else if ( j == (camadT + 1) .and. j /= n ) then

                TEdwJ0(:,j) = TEdwJ0(:,j - 1) * ( cdexp( -uJ0(:,camadT) * ( prof(camadT) - h0 ) ) + &
                        RTEupJ0(:,camadT) * FEupJ0(:) * cdexp( -uhJ0(:,camadT) ) + &
                        RTEdwJ0(:,camadT) * FEdwJ0(:) ) / ( 1.d0 + RTEdwJ0(:,j) * cdexp( -2.d0 * uhJ0(:,j) ) )
                TEdwJ1(:,j) = TEdwJ1(:,j - 1) * ( cdexp( -uJ1(:,camadT) * ( prof(camadT) - h0 ) ) + &
                        RTEupJ1(:,camadT) * FEupJ1(:) * cdexp( -uhJ1(:,camadT) ) + &
                        RTEdwJ1(:,camadT) * FEdwJ1(:) ) / ( 1.d0 + RTEdwJ1(:,j) * cdexp( -2.d0 * uhJ1(:,j) ) )

 			else if ( j /= n ) then

                TEdwJ0(:,j) = TEdwJ0(:,j - 1) * ( 1.d0 + RTEdwJ0(:,j - 1) ) * cdexp( -uhJ0(:,j - 1) ) / &
                        ( 1.d0 + RTEdwJ0(:,j) * cdexp( -2.d0 * uhJ0(:,j) ) )
                TEdwJ1(:,j) = TEdwJ1(:,j - 1) * ( 1.d0 + RTEdwJ1(:,j - 1) ) * cdexp( -uhJ1(:,j - 1) ) / &
                        ( 1.d0 + RTEdwJ1(:,j) * cdexp( -2.d0 * uhJ1(:,j) ) )

 			else if ( j == n ) then

                TEdwJ0(:,j) = TEdwJ0(:,j - 1) * ( 1.d0 + RTEdwJ0(:,j - 1) ) * cdexp( -uhJ0(:,j - 1) )
                TEdwJ1(:,j) = TEdwJ1(:,j - 1) * ( 1.d0 + RTEdwJ1(:,j - 1) ) * cdexp( -uhJ1(:,j - 1) )

 			end if
 		end do
 	else if ( camad < camadT ) then
        allocate( TEupJ0(nJ0,camad : camadT), TEupJ1(nJ1,camad : camadT) )
        do j = camadT, camad, -1
            if ( j == camadT ) then

                TEupJ0(:,j) = zeta * mz / ( 2.d0 * uJ0(:,j) )
                TEupJ1(:,j) = zeta * mz / ( 2.d0 * uJ1(:,j) )

 			else if ( j == (camadT - 1) .and. j == 0 ) then

                TEupJ0(:,j) = TEupJ0(:,j + 1) * ( cdexp( -uJ0(:,camadT) * h0 ) + &
                    RTEupJ0(:,camadT) * FEupJ0(:) + RTEdwJ0(:,camadT) * FEdwJ0(:) * cdexp( -uhJ0(:,camadT) ) )
                TEupJ1(:,j) = TEupJ1(:,j + 1) * ( cdexp( -uJ1(:,camadT) * h0 ) + &
                    RTEupJ1(:,camadT) * FEupJ1(:) + RTEdwJ1(:,camadT) * FEdwJ1(:) * cdexp( -uhJ1(:,camadT) ) )

 			else if ( j == (camadT - 1) .and. j /= 0 ) then

                TEupJ0(:,j) = TEupJ0(:,j + 1) * ( cdexp(uJ0(:,camadT) * ( prof(camadT - 1) - h0 ) ) + &
                    RTEupJ0(:,camadT) * FEupJ0(:) + RTEdwJ0(:,camadT) * FEdwJ0(:) * &
                    cdexp( -uhJ0(:,camadT) ) ) / ( 1.d0 + RTEupJ0(:,j) * cdexp( -2.d0 * uhJ0(:,j) ) )
                TEupJ1(:,j) = TEupJ1(:,j + 1) * ( cdexp( uJ1(:,camadT) * ( prof(camadT - 1) - h0 ) ) + &
                    RTEupJ1(:,camadT) * FEupJ1(:) + RTEdwJ1(:,camadT) * FEdwJ1(:) * &
                    cdexp( -uhJ1(:,camadT) ) ) / ( 1.d0 + RTEupJ1(:,j) * cdexp( -2.d0 * uhJ1(:,j) ) )

 			else if ( j /= 0 ) then

                TEupJ0(:,j) = TEupJ0(:,j + 1) * ( 1.d0 + RTEupJ0(:,j + 1) ) * cdexp( -uhJ0(:,j + 1) ) / &
                    ( 1.d0 + RTEupJ0(:,j) * cdexp( -2.d0 * uhJ0(:,j) ) )
                TEupJ1(:,j) = TEupJ1(:,j + 1) * ( 1.d0 + RTEupJ1(:,j + 1) ) * cdexp( -uhJ1(:,j + 1) ) / &
                    ( 1.d0 + RTEupJ1(:,j) * cdexp( -2.d0 * uhJ1(:,j) ) )

 			else if ( j == 0 ) then

                TEupJ0(:,j) = TEupJ0(:,1) * ( 1.d0 + RTEupJ0(:,1) ) * cdexp( -uhJ0(:,1) )
                TEupJ1(:,j) = TEupJ1(:,1) * ( 1.d0 + RTEupJ1(:,1) ) * cdexp( -uhJ1(:,1) )

 			end if
 		end do
 	else
 	
        allocate( TEdwJ0(nJ0,camadT : camad), TEdwJ1(nJ1,camadT : camad) )
        allocate( TEupJ0(nJ0,camad : camadT), TEupJ1(nJ1,camad : camadT) )

		TEdwJ0(:,camad) = zeta * mz / ( 2.d0 * uJ0(:,camadT) )
		TEdwJ1(:,camad) = zeta * mz / ( 2.d0 * uJ1(:,camadT) )
 			
		TEupJ0(:,camad) = TEdwJ0(:,camad)
		TEupJ1(:,camad) = TEdwJ1(:,camad)

 	end if

 	allocate( Kte_J0(nJ0), Kte_J1(nJ1), Ktedz_J1(nJ1) )
 	allocate( kernelExJ0(nJ0), kernelExJ1(nJ1), kernelEyJ0(nJ0), kernelEyJ1(nJ1) )
 	allocate( kernelHxJ0(nJ0), kernelHxJ1(nJ1), kernelHyJ0(nJ0), kernelHyJ1(nJ1), kernelHzJ1(nJ1) )
 	if ( camad == 0 .and. camadT /= 0 ) then

 		Kte_J0 = TEupJ0(:,0) * cdexp( uJ0(:,0) * z ) * w_J0(:)
 		Kte_J1 = TEupJ1(:,0) * cdexp( uJ1(:,0) * z ) * w_J1(:)
 		Ktedz_J1 = AdmIntJ1(:,0) * Kte_J1

 		kernelExJ1 = y * Kte_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
 		Ex_p = sum( kernelExJ1 ) / r     !este último r é decorrente do uso dos filtros
 
 		kernelEyJ1 = - x * Kte_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
 		Ey_p = sum( kernelEyJ1 ) / r     !este último r é decorrente do uso dos filtros
 
 		kernelHxJ1 = - x * Ktedz_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
 		Hx_p = sum( kernelHxJ1 ) / r     !este último r é decorrente do uso dos filtros

 		kernelHyJ1 = - y * Ktedz_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
 		Hy_p = sum( kernelHyJ1 ) / r     !este último r é decorrente do uso dos filtros
 
 		kernelHzJ1 = Kte_J0 * krJ0 * krJ0 * krJ0 / ( 2.d0 * pi * zeta )
 		Hz_p = sum( kernelHzJ1 ) / r     !este último r é decorrente do uso dos filtros
 
    else if ( camad < camadT ) then !camada k

 		Kte_J0 = ( TEupJ0(:,camad) * ( cdexp( uJ0(:,camad) * ( z - prof(camad) ) ) + &
                RTEupJ0(:,camad) * cdexp( -uJ0(:,camad) * ( z - prof(camad - 1) + h(camad) ) ) ) ) * w_J0(:)
 		Kte_J1 = ( TEupJ1(:,camad) * ( cdexp( uJ1(:,camad) * ( z - prof(camad) ) ) + &
                RTEupJ1(:,camad) * cdexp( -uJ1(:,camad) * ( z - prof(camad - 1) + h(camad) ) ) ) ) * w_J1(:)
 		Ktedz_J1 = ( AdmIntJ1(:,camad) * TEupJ1(:,camad) * ( cdexp( uJ1(:,camad) * ( z - prof(camad) ) ) - &
                RTEupJ1(:,camad) * cdexp( -uJ1(:,camad) * ( z - prof(camad - 1) + h(camad) ) ) ) ) * w_J1(:)
 
 		kernelExJ1 = y * Kte_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
 		Ex_p = sum( kernelExJ1 ) / r     !este último r é decorrente do uso dos filtros
 
 		kernelEyJ1 = - x * Kte_J1 * krJ1 ** 2 / ( 2.d0 * pi * r )
 		Ey_p = sum( kernelEyJ1 ) / r     !este último r é decorrente do uso dos filtros
 
 		kernelHxJ1 = - x * Ktedz_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
 		Hx_p = sum( kernelHxJ1 ) / r     !este último r é decorrente do uso dos filtros

 		kernelHyJ1 = - y * Ktedz_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
 		Hy_p = sum( kernelHyJ1 ) / r     !este último r é decorrente do uso dos filtros
 
 		kernelHzJ1 = Kte_J0 * krJ0 * krJ0 * krJ0 / ( 2.d0 * pi * zeta )
 		Hz_p = sum( kernelHzJ1 ) / r     !este último r é decorrente do uso dos filtros
 
 	else if ( camad == camadT .and. z <= h0 ) then	!na mesma camada do transmissor mas acima dele

 		Kte_J0 = ( TEupJ0(:,camad) * ( cdexp( uJ0(:,camad) * ( z - h0 ) ) + &
                RTEupJ0(:,camad) * FEupJ0(:) * cdexp( -uJ0(:,camad) * ( z - prof(camad - 1) ) ) + &
                RTEdwJ0(:,camad) * FEdwJ0(:) * cdexp( uJ0(:,camad) * ( z - prof(camad) ) ) ) ) * w_J0(:)
 		Kte_J1 = ( TEupJ1(:,camad) * ( cdexp( uJ1 (:,camad) * ( z - h0 ) ) + &
                RTEupJ1(:,camad) * FEupJ1(:) * cdexp( -uJ1(:,camad) * ( z - prof(camad - 1) ) ) + &
                RTEdwJ1(:,camad) * FEdwJ1(:) * cdexp( uJ1(:,camad) * ( z - prof(camad) ) ) ) ) * w_J1(:)
		Ktedz_J1 = ( AdmIntJ1(:,camad) * TEupJ1(:,camad) * ( cdexp( uJ1(:,camad) * ( z - h0 ) ) - &
                RTEupJ1(:,camad) * FEupJ1(:) * cdexp( -uJ1(:,camad) * ( z - prof(camad - 1) ) ) + &
                RTEdwJ1(:,camad) * FEdwJ1(:) * cdexp( uJ1(:,camad) * ( z - prof(camad) ) ) ) ) * w_J1(:)

 		kernelExJ1 = y * Kte_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
 		Ex_p = sum( kernelExJ1 ) / r     !este último r é decorrente do uso dos filtros
 
 		kernelEyJ1 = - x * Kte_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
 		Ey_p = sum( kernelEyJ1 ) / r     !este último r é decorrente do uso dos filtros
 
 		kernelHxJ1 = - x * Ktedz_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
 		Hx_p = sum( kernelHxJ1 ) / r     !este último r é decorrente do uso dos filtros

 		kernelHyJ1 = - y * Ktedz_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
 		Hy_p = sum( kernelHyJ1 ) / r     !este último r é decorrente do uso dos filtros
 
 		kernelHzJ1 = Kte_J0 * krJ0 * krJ0 * krJ0 / ( 2.d0 * pi * zeta )
 		Hz_p = sum( kernelHzJ1 ) / r     !este último r é decorrente do uso dos filtros
 
 	else if ( camad == camadT .and. z > h0 ) then	!na mesma camada do transmissor mas abaixo dele

 		Kte_J0 = ( TEdwJ0(:,camad) * ( cdexp( -uJ0(:,camad) * ( z - h0 ) ) + &
                RTEupJ0(:,camad) * FEupJ0(:) * cdexp( -uJ0(:,camad) * ( z - prof(camad - 1) ) ) + &
                RTEdwJ0(:,camad) * FEdwJ0(:) * cdexp( uJ0(:,camad) * ( z - prof(camad) ) ) ) ) * w_J0(:)
 		Kte_J1 = ( TEdwJ1(:,camad) * ( cdexp( -uJ1(:,camad) * ( z - h0 ) ) + &
                RTEupJ1(:,camad) * FEupJ1(:) * cdexp( -uJ1(:,camad) * ( z - prof(camad - 1) ) ) + &
                RTEdwJ1(:,camad) * FEdwJ1(:) * cdexp( uJ1(:,camad) * ( z - prof(camad) ) ) ) ) * w_J1(:)
 		Ktedz_J1 = ( - AdmIntJ1(:,camad) * TEdwJ1(:,camad) * ( cdexp( -uJ1(:,camad) * ( z - h0 ) ) + &
                RTEupJ1(:,camad) * FEupJ1(:) * cdexp( -uJ1(:,camad) * ( z - prof(camad - 1) ) ) - &
                RTEdwJ1(:,camad) * FEdwJ1(:) * cdexp( uJ1(:,camad) * ( z - prof(camad) ) ) ) ) * w_J1(:)
 
 		kernelExJ1 = y * Kte_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
 		Ex_p = sum( kernelExJ1 ) / r     !este último r é decorrente do uso dos filtros
 
 		kernelEyJ1 = - x * Kte_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
 		Ey_p = sum( kernelEyJ1 ) / r     !este último r é decorrente do uso dos filtros
 
 		kernelHxJ1 = - x * Ktedz_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
 		Hx_p = sum( kernelHxJ1 ) / r     !este último r é decorrente do uso dos filtros

 		kernelHyJ1 = - y * Ktedz_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
 		Hy_p = sum( kernelHyJ1 ) / r     !este último r é decorrente do uso dos filtros
 
 		kernelHzJ1 = Kte_J0 * krJ0 * krJ0 * krJ0 / ( 2.d0 * pi * zeta )
 		Hz_p = sum( kernelHzJ1 ) / r     !este último r é decorrente do uso dos filtros
 
 	else if ( camad > camadT .and. camad /= n ) then !camada j

 		Kte_J0 = ( TEdwJ0(:,camad) * ( cdexp( -uJ0(:,camad) * ( z - prof(camad - 1) ) ) + &
                RTEdwJ0(:,camad) * cdexp( uJ0(:,camad) * ( z - prof(camad) - h(camad) ) ) ) ) * w_J0(:)
 		Kte_J1 = ( TEdwJ1(:,camad) * ( cdexp( -uJ1(:,camad) * ( z - prof(camad - 1) ) ) + &
                RTEdwJ1(:,camad) * cdexp( uJ1(:,camad) * ( z - prof(camad) - h(camad) ) ) ) ) * w_J1(:)
 		Ktedz_J1 = ( - AdmIntJ1(:,camad) * TEdwJ1(:,camad) * ( cdexp( -uJ1(:,camad) * ( z - prof(camad - 1) ) ) - &
                RTEdwJ1(:,camad) * cdexp( uJ1(:,camad) * ( z - prof(camad) - h(camad) ) ) ) ) * w_J1(:)
 
 		kernelExJ1 = y * Kte_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
 		Ex_p = sum( kernelExJ1 ) / r     !este último r é decorrente do uso dos filtros
 
 		kernelEyJ1 = - x * Kte_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
 		Ey_p = sum( kernelEyJ1 ) / r     !este último r é decorrente do uso dos filtros
 
 		kernelHxJ1 = - x * Ktedz_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
 		Hx_p = sum( kernelHxJ1 ) / r     !este último r é decorrente do uso dos filtros

 		kernelHyJ1 = - y * Ktedz_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
 		Hy_p = sum( kernelHyJ1 ) / r     !este último r é decorrente do uso dos filtros
 
 		kernelHzJ1 = Kte_J0 * krJ0 * krJ0 * krJ0 / ( 2.d0 * pi * zeta )
 		Hz_p = sum( kernelHzJ1 ) / r     !este último r é decorrente do uso dos filtros

 	else	!camada n

 		Kte_J0 = ( TEdwJ0(:,n) * cdexp( -uJ0(:,n) * ( z - prof(n - 1) ) ) ) * w_J0(:)
 		Kte_J1 = ( TEdwJ1(:,n) * cdexp( -uJ1(:,n) * ( z - prof(n - 1) ) ) ) * w_J1(:)
 		Ktedz_J1=( - AdmIntJ1(:,n) * TEdwJ1(:,n) * cdexp( -uJ1(:,n) * ( z - prof(n - 1) ) ) ) * w_J1(:)
 
 		kernelExJ1 = y * Kte_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
 		Ex_p = sum( kernelExJ1 ) / r     !este último r é decorrente do uso dos filtros
 
 		kernelEyJ1 = - x * Kte_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
 		Ey_p = sum( kernelEyJ1 ) / r     !este último r é decorrente do uso dos filtros
 
 		kernelHxJ1 = - x * Ktedz_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
 		Hx_p = sum( kernelHxJ1 ) / r     !este último r é decorrente do uso dos filtros

 		kernelHyJ1 = - y * Ktedz_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
 		Hy_p = sum( kernelHyJ1 ) / r     !este último r é decorrente do uso dos filtros
 
 		kernelHzJ1 = Kte_J0 * krJ0 * krJ0 * krJ0 / ( 2.d0 * pi * zeta )
 		Hz_p = sum( kernelHzJ1 ) / r     !este último r é decorrente do uso dos filtros
 
 	end if
 
 	deallocate( h, KrJ0, KrJ1, w_J0, w_J1 )
 	deallocate( wvnb2, uJ0, uJ1, AdmIntJ0, AdmIntJ1, uhJ0, uhJ1, tghJ0, tghJ1 )
 	deallocate( AdmApdwJ0, AdmApdwJ1, RTEdwJ0, RTEdwJ1 )
 	deallocate( AdmApupJ0, AdmApupJ1, RTEupJ0, RTEupJ1 )
 	deallocate( FEdwJ0, FEdwJ1, FEupJ0, FEupJ1 )
 	deallocate( Kte_J0, Kte_J1, Ktedz_J1 )
 	deallocate( kernelExJ0, kernelExJ1, kernelEyJ0, kernelEyJ1 )
 	deallocate( kernelHxJ0, kernelHxJ1, kernelHyJ0, kernelHyJ1, kernelHzJ1 )
end subroutine dmv_xyz_loops
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

	subroutine constfiltro(tipofiltro,criador,idtfcd_f,np_f,r,Kr,w)

	implicit none
	real(8),intent(in)::r
	integer(4),intent(in)::tipofiltro,criador,idtfcd_f,np_f
	real(8),intent(out)::Kr(np_f),w(np_f)
	real(8)::abs_f(np_f)

	real(8),dimension(:),allocatable::lr_abs_J061,lr_pes_J061,lr_abs_J147,lr_pes_J147
	real(8),dimension(:),allocatable::fl_abs_J037,fl_pes_J037,fl_abs_J127,fl_pes_J127,fl_abs_J219,fl_pes_J219
	real(8),dimension(:),allocatable::gs_abs_J061,gs_pes_J061,gs_abs_J0120,gs_pes_J0120
	real(8),dimension(:),allocatable::gs_abs_J147,gs_pes_J147,gs_abs_J1140,gs_pes_J1140
	real(8),dimension(:),allocatable::k_abs_J0161,k_pes_J061,k_pes_J161,k_abs_J01241,k_pes_J0241,k_pes_J1241
	real(8),dimension(:),allocatable::kk_abs_J01101,kk_pes_J0101,kk_pes_J1101,kk_abs_J01201,kk_pes_J0201,kk_pes_J1201
	real(8),dimension(:),allocatable::kk_abs_J01401,kk_pes_J0401,kk_pes_J1401

	real(8),dimension(:),allocatable::fl_abs_s19,fl_pes_sen19,fl_abs_s30,fl_pes_sen30,fl_abs_s40,fl_pes_sen40
	real(8),dimension(:),allocatable::fl_abs_c19,fl_pes_cos19,fl_abs_c30,fl_pes_cos30,fl_abs_c40,fl_pes_cos40
	real(8),dimension(:),allocatable::kk_abs_81,kk_pes_sen81,kk_pes_cos81,kk_abs_241,kk_pes_sen241,kk_pes_cos241
	real(8),dimension(:),allocatable::kk_abs_601,kk_pes_sen601,kk_pes_cos601
!	tipofiltro é uma variável que identifica se os filtros a serem usados serão para transformadas de Hankel (0) ou de Fourier (1)
	if (tipofiltro == 0) then !Para uso de filtros J0, J1 ou J2
		if (criador==0) then !Para o uso dos filtros de Rijo
			allocate(lr_abs_J061(61),lr_pes_J061(61),lr_abs_J147(47),lr_pes_J147(47))
			call J0J1Rijo(lr_abs_J061,lr_pes_J061,lr_abs_J147,lr_pes_J147)
			if (idtfcd_f == 0) then
				abs_f=lr_abs_J061
				w=lr_pes_J061
			elseif (idtfcd_f == 1) then
				abs_f=lr_abs_J147
				w=lr_pes_J147
			end if
			deallocate(lr_abs_J061,lr_pes_J061,lr_abs_J147,lr_pes_J147)
		kr=dexp(abs_f)/dabs(r)
		elseif (criador==1) then !Para o uso dos filtros de Frayzer
			allocate(fl_abs_J037(37),fl_pes_J037(37),fl_abs_J127(27),fl_pes_J127(27),fl_abs_J219(19),fl_pes_J219(19))
			call J0J1J2Frayzer(fl_abs_J037,fl_pes_J037,fl_abs_J127,fl_pes_J127,fl_abs_J219,fl_pes_J219)
			if (idtfcd_f == 0) then
				abs_f=fl_abs_J037
				w=fl_pes_J037
			elseif (idtfcd_f == 1) then
				abs_f=fl_abs_J127
				w=fl_pes_J127
			end if
			deallocate(fl_abs_J037,fl_pes_J037,fl_abs_J127,fl_pes_J127,fl_abs_J219,fl_pes_J219)
		kr=dexp(abs_f)/dabs(r)
		elseif (criador==2) then !Para o uso dos filtros de Guptasarma
			allocate(gs_abs_J061(61),gs_pes_J061(61),gs_abs_J0120(120),gs_pes_J0120(120),gs_abs_J147(47),gs_pes_J147(47))
			allocate(gs_abs_J1140(140),gs_pes_J1140(140))
			call J0J1GS(gs_abs_J061,gs_pes_J061,gs_abs_J0120,gs_pes_J0120,gs_abs_J147,gs_pes_J147,gs_abs_J1140,gs_pes_J1140)
			if (idtfcd_f == 0) then
				if (np_f == 61) then
					abs_f=gs_abs_J061
					w=gs_pes_J061
				elseif (np_f == 120) then
					abs_f=gs_abs_J0120
					w=gs_pes_J0120
				end if
			elseif (idtfcd_f == 1) then
				if (np_f == 47) then
					abs_f=gs_abs_J147
					w=gs_pes_J147
				elseif (np_f == 140) then
					abs_f=gs_abs_J1140
					w=gs_pes_J1140
				end if
			end if
			deallocate(gs_abs_J061,gs_pes_J061,gs_abs_J0120,gs_pes_J0120,gs_abs_J147,gs_pes_J147,gs_abs_J1140,gs_pes_J1140)
		kr=(10.d0**(abs_f))/dabs(r)
		elseif (criador==3) then !Para o uso dos filtros de Kong
			allocate(k_abs_J0161(61),k_pes_J061(61),k_pes_J161(61),k_abs_J01241(241),k_pes_J0241(241),k_pes_J1241(241))
			call J0J1Kong(k_abs_J0161,k_pes_J061,k_pes_J161,k_abs_J01241,k_pes_J0241,k_pes_J1241)
			if (idtfcd_f == 0) then
				if (np_f == 61) then
					abs_f=k_abs_J0161
					w=k_pes_J061
				elseif (np_f == 241) then
					abs_f=k_abs_J01241
					w=k_pes_J0241
				end if
			elseif (idtfcd_f == 1) then
				if (np_f == 61) then
					abs_f=k_abs_J0161
					w=k_pes_J161
				elseif (np_f == 241) then
					abs_f=k_abs_J01241
					w=k_pes_J1241
				end if
			end if
			deallocate(k_abs_J0161,k_pes_J061,k_pes_J161,k_abs_J01241,k_pes_J0241,k_pes_J1241)
		kr=abs_f/dabs(r)
		elseif (criador==4) then !Para o uso dos filtros de Key
			allocate(kk_abs_J01101(101),kk_pes_J0101(101),kk_pes_J1101(101),kk_abs_J01201(201),kk_pes_J0201(201))
			allocate(kk_pes_J1201(201),kk_abs_J01401(401),kk_pes_J0401(401),kk_pes_J1401(401))
			call J0J1Key(kk_abs_J01101,kk_pes_J0101,kk_pes_J1101,kk_abs_J01201,kk_pes_J0201,kk_pes_J1201, &
					kk_abs_J01401,kk_pes_J0401,kk_pes_J1401)
			if (idtfcd_f == 0) then
				if (np_f == 101) then
					abs_f=kk_abs_J01101
					w=kk_pes_J0101
				elseif (np_f == 201) then
					abs_f=kk_abs_J01201
					w=kk_pes_J0201
				elseif (np_f == 401) then
					abs_f=kk_abs_J01401
					w=kk_pes_J0401
				end if
			elseif (idtfcd_f == 1) then
				if (np_f == 101) then
					abs_f=kk_abs_J01101
					w=kk_pes_J1101
				elseif (np_f == 201) then
					abs_f=kk_abs_J01201
					w=kk_pes_J1201
				elseif (np_f == 401) then
					abs_f=kk_abs_J01401
					w=kk_pes_J1401
				end if
			end if
			deallocate(kk_abs_J01101,kk_pes_J0101,kk_pes_J1101,kk_abs_J01201,kk_pes_J0201)
			deallocate(kk_pes_J1201,kk_abs_J01401,kk_pes_J0401,kk_pes_J1401)
		kr=abs_f/dabs(r) 
		end if
!	elseif (tipofiltro == 1) then !Para uso de filtros seno ou cosseno
	else
		print*,'especificacao incorreta do tipo de filtro'
		stop
	end if

	end subroutine

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!Atenção: Os filtros de Kong e de Key que tem os valores identificados como de abscissas,
!	na verdade já são os valores de exp(s1+(i-1)*del) ou 10**(s1+(i-1)*del) (dependendo da técnica usada pelo autor)

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
	subroutine J0J1Rijo(lr_abs_J061,lr_pes_J061,lr_abs_J147,lr_pes_J147)
	implicit none
	real(8),intent(out)::lr_abs_J061(61),lr_pes_J061(61),lr_abs_J147(47),lr_pes_J147(47)
	integer(4)::ii
	real(8)::s1RJ061,delRJ061,s1RJ147,delRJ147
!	Filtro J0 do Rijo (61 pontos)
	s1RJ061=-11.702888735142237d0
	delRJ061=0.2685696197447527d0
	lr_abs_J061=(/s1RJ061+(/(ii,ii=0,60)/)*delRJ061/)
	lr_pes_J061=(/3.30220475766d-4, -1.18223623458d-3, 2.01879495264d-3,&
			-2.13218719891d-3, 1.60839063172d-3, -9.09156346708d-4,&
			4.37889252738d-4, -1.55298878782d-4, 7.98411962729d-5,&
			4.37268394072d-6, 3.94253441247d-5, 4.02675924344d-5,&
			5.66053344653d-5, 7.25774926389d-5, 9.55412535465d-5,&
			1.24699163157d-4, 1.63262166579d-4, 2.13477133718d-4,&
			2.79304232173d-4, 3.65312787897d-4, 4.77899413107d-4,&
			6.25100170825d-4, 8.17726956451d-4, 1.06961339341d-3,&
			1.39920928148d-3, 1.83020380399d-3, 2.39417015791d-3,&
			3.13158560774d-3, 4.09654426763d-3, 5.35807925630d-3,&
			7.00889482693d-3, 9.16637526490d-3, 1.19891721272d-2,&
			1.56755740646d-2, 2.04953856060d-2, 2.67778388247d-2,&
			3.49719672729d-2, 4.55975312615d-2, 5.93498881451d-2,&
			7.69179091244d-2, 9.91094769804d-2, 1.26166963993d-1,&
			1.57616825575d-1, 1.89707800260d-1, 2.13804195282d-1,&
			2.08669340316d-1, 1.40250562745d-1, -3.65385242807d-2,&
			-2.98004010732d-1, -4.2189814924d-1, 5.94373771266d-2,&
			5.29621428353d-1, -4.41362405166d-1, 1.90355040550d-1,&
			-6.19966386785d-2, 1.87255115744d-2, -5.68736766738d-3,&
			1.68263510609d-3, -4.38587145792d-4, 8.59117336292d-5,&
			-9.15853765160d-6 /)

!	Filtro J1 do Rijo (47 pontos)
	s1RJ147=-7.024684869538879d0
	delRJ147=0.2546636319446458d0
	lr_abs_J147=(/s1RJ147+(/(ii,ii=0,46)/)*delRJ147/)
	lr_pes_J147=(/3.17926147465d-6, -9.73811660718d-6, 1.64866227408d-5, &
			-1.81501261160d-5, 1.87556556369d-5, -1.46550406038d-5, &
			1.53799733803d-5, -6.95628273934d-6, 1.41881555665d-5,  &
			3.41445665537d-6, 2.13941715512d-5, 2.34962369042d-5, &
			4.84340283290d-5, 7.33732978590d-5, 1.27703784430d-4, &
			2.08120025730d-4, 3.49803898913d-4, 5.79107814687d-4, &
			9.65887918451d-4, 1.60401273703d-3, 2.66903777685d-3, &
			4.43111590040d-3, 7.35631696247d-3, 1.21782796293d-2, &
			2.01097829218d-2, 3.30096953061d-2, 5.37143591532d-2, &
			8.60516613299d-2, 1.34267607144d-1, 2.00125033067d-1, &
			2.74027505792d-1, 3.18168749246d-1, 2.41655667461d-1, &
			-5.40549161658d-2, -4.46912952135d-1, -1.92231885629d-1, &
			5.52376753950d-1, -3.57429049025d-1, 1.41510519002d-1, &
			-4.61421935309d-2, 1.48273761923d-2, -5.07479209193d-3, &
			1.83829713749d-3, -6.67742804324d-4, 2.21277518118d-4, &
			-5.66248732755d-5, 7.88229202853d-6/)
	end subroutine J0J1Rijo

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

	subroutine J0J1J2Frayzer(fl_abs_J037,fl_pes_J037,fl_abs_J127,fl_pes_J127,fl_abs_J219,fl_pes_J219)
	implicit none
	real(8),intent(out)::fl_abs_J037(37),fl_pes_J037(37),fl_abs_J127(27),fl_pes_J127(27),fl_abs_J219(19),fl_pes_J219(19)

	integer(4)::ii
	real(8)::s1_J037,del_J037,s1_J127,del_J127,s1_J219,del_J219
!	Filtro J0. 37 pontos (Frayzer)
	s1_J037 = -8.7915450115593092d0
	del_J037 = 0.34398729728929939d0
	fl_abs_J037=(/s1_J037+(/(ii,ii=0,36)/)*del_J037/)
	fl_pes_J037=(/0.00238335617914d0, -0.00923451277729d0, 0.02079895052544d0, -0.03338702409972d0, &
			0.04477739890139d0, -0.05167765846960d0, 0.05587472556026d0, -0.05512276066568d0, &
			0.05444557860139d0, -0.04892210690497d0, 0.04732656195005d0, -0.03860491533605d0, &
			0.03920960657489d0, -0.02644804816704d0, 0.03260982461445d0, -0.01228397713755d0, &
			0.02956411797486d0, 0.00595004660943d0, 0.03317515498364d0, 0.03274470853033d0, &
			0.04933875502497d0, 0.07661759525931d0, 0.08883650000582d0, 0.15133386760606d0, &
			0.16507829454944d0, 0.26169301627952d0, 0.24688255898384d0, 0.27750176445175d0, &
			-0.00155554395068d0, -0.33630611758063d0, -0.50651480451307d0, 0.59211451923489d0, &
			-0.23450184447820d0, 0.05446237320377d0, -0.00928958955398d0, 0.00119894072586d0, &
			-0.00008971645558d0 /)
!	Filtro J1. 27pontos (Frayzer)
	s1_J127=-4.61195591517857d0
	del_J127=0.33288944128788d0
	fl_abs_J127=(/s1_J127+(/(ii,ii=0,26)/)*del_J127/)
	fl_pes_J127=(/0.00022403915103d0, -0.00063831184603d0, 0.00113167818280d0, -0.00081481684918d0, &
			0.00065948253261d0, 0.00015110345657d0, 0.00165282210988d0, 0.00080418132518d0, &
			0.00380611316285d0, 0.00619808024688d0, 0.01404440767965d0, 0.02189018716297d0, &
			0.04964931228922d0, 0.08370610919229d0, 0.16660930791854d0, 0.25725196330244d0, &
			0.40730714915038d0, 0.34895491488489d0, -0.04241399715350d0, -0.65781916295135d0, &
			0.46094644953685d0, -0.15431387506167d0, 0.03858981403104d0, -0.00949545480093d0, &
			0.00231033358900d0, -0.00043488073049d0, 0.00004191435466d0 /)
!	Filtro J2. 19 pontos (Frayzer)
	s1_J219=-1.55050598185024d0
	del_J219=0.18634701555858d0
	fl_abs_J219=(/s1_J219+(/(ii,ii=0,18)/)*del_J219/)
	fl_pes_J219=(/0.00799552054539d0, -0.02611146942531d0, 0.03908786880947d0, -0.00917038216156d0, &
			-0.03460228618293d0, 0.01193259214300d0, 0.06058316111850d0, 0.02489345333357d0, &
			-0.05843364695515d0, -0.04064913677478d0, 0.12862574793110d0, 0.25252742978437d0, &
			0.11110399521022d0, -0.11496118208512d0, 0.10097209623246d0, 0.80339619274433d0, &
			0.82739761272033d0, -1.49990618049328d0, 0.41518056949677d0 /)

	end subroutine J0J1J2Frayzer

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

	subroutine J0J1GS(gs_abs_J061,gs_pes_J061,gs_abs_J0120,gs_pes_J0120,gs_abs_J147,gs_pes_J147,gs_abs_J1140,gs_pes_J1140)
	implicit none
	real(8),intent(out)::gs_abs_J061(61),gs_pes_J061(61),gs_abs_J0120(120),gs_pes_J0120(120)
	real(8),intent(out)::gs_abs_J147(47),gs_pes_J147(47),gs_abs_J1140(140),gs_pes_J1140(140)

	integer(4)::ii
	real(8)::s1_J061,del_J061,s1_J0120,del_J0120,s1_J147,del_J147,s1_J1140,del_J1140
!	Filtro J0. 61 pontos (Guptasarma e Singh, 1997)
	s1_J061=-5.0825d0
	del_J061=1.16638303862d-1
	gs_abs_J061=(/s1_J061+(/(ii,ii=0,60)/)*del_J061/)
	gs_pes_J061=(/3.30220475766d-4,-1.18223623458d-3,2.01879495264d-3,-2.13218719891d-3, &
			1.60839063172d-3,-9.09156346708d-4,4.37889252738d-4,-1.55298878782d-4, &
			7.98411962729d-5,4.37268394072d-6,3.94253441247d-5,4.02675924344d-5, &
			5.66053344653d-5,7.25774926389d-5,9.55412535465d-5,1.24699163157d-4, &
			1.63262166579d-4,2.13477133718d-4,2.79304232173d-4,3.65312787897d-4, &
			4.77899413107d-4,6.25100170825d-4,8.17726956451d-4,1.06961339341d-3, &
			1.39920928148d-3,1.83020380399d-3,2.39417015791d-3,3.13158560774d-3, &
			4.09654426763d-3,5.35807925630d-3,7.00889482693d-3,9.16637526490d-3, &
			1.19891721272d-2,1.56755740646d-2,2.04953856060d-2,2.67778388247d-2, &
			3.49719672729d-2,4.55975312615d-2,5.93498881451d-2,7.69179091244d-2, &
			9.91094769804d-2,1.26166963993d-1,1.57616825575d-1,1.89707800260d-1, &
			2.13804195282d-1,2.08669340316d-1,1.40250562745d-1,-3.65385242807d-2, &
			-2.98004010732d-1,-4.21898149249d-1,5.94373771266d-2,5.29621428353d-1, &
			-4.41362405166d-1,1.90355040550d-1,-6.19966386785d-2,1.87255115744d-2, &
			-5.68736766738d-3,1.68263510609d-3,-4.38587145792d-4,8.59117336292d-5, &
			-9.15853765160d-6/)
!	Filtro J0. 120 pontos (Guptasarma e Singh, 1997)
	s1_J0120=-8.3885d0
	del_J0120=9.04226468670d-2
	gs_abs_J0120=(/s1_J0120+(/(ii,ii=0,119)/)*del_J0120/)
	gs_pes_J0120=(/9.62801364263d-7,-5.02069203805d-6,1.25268783953d-5,-1.99324417376d-5, &
			2.29149033546d-5,-2.04737583809d-5,1.49952002937d-5,-9.37502840980d-6, &
			5.20156955323d-6,-2.62939890538d-6,1.26550848081d-6,-5.73156151923d-7, &
			2.76281274155d-7,-1.09963734387d-7,7.38038330280d-8,-9.31614600001d-9, &
			3.87247135578d-8,2.10303178461d-8,4.10556513877d-8,4.13077946246d-8, &
			5.68828741789d-8,6.59543638130d-8,8.40811858728d-8,1.01532550003d-7, &
			1.26437360082d-7,1.54733678097d-7,1.91218582499d-7,2.35008851918d-7, &
			2.89750329490d-7,3.56550504341d-7,4.39299297826d-7,5.40794544880d-7, &
			6.66136379541d-7,8.20175040653d-7,1.01015545059d-6,1.24384500153d-6, &
			1.53187399787d-6,1.88633707689d-6,2.32307100992d-6,2.86067883258d-6, &
			3.52293208580d-6,4.33827546442d-6,5.34253613351d-6,6.57906223200d-6, &
			8.10198829111d-6,9.97723263578d-6,1.22867312381d-5,1.51305855976d-5, &
			1.86329431672d-5,2.29456891669d-5,2.82570465155d-5,3.47973610445d-5, &
			4.28521099371d-5,5.27705217882d-5,6.49856943660d-5,8.00269662180d-5, &
			9.85515408752d-5,1.21361571831d-4,1.49454562334d-4,1.84045784500d-4, &
			2.26649641428d-4,2.79106748890d-4,3.43716968725d-4,4.23267056591d-4, &
			5.21251001943d-4,6.41886194381d-4,7.90483105615d-4,9.73420647376d-4, &
			1.19877439042d-3,1.47618560844d-3,1.81794224454d-3,2.23860214971d-3, &
			2.75687537633d-3,3.39471308297d-3,4.18062141752d-3,5.14762977308d-3, &
			6.33918155348d-3,7.80480111772d-3,9.61064602702d-3,1.18304971234d-2, &
			1.45647517743d-2,1.79219149417d-2,2.20527911163d-2,2.71124775541d-2, &
			3.33214363101d-2,4.08864842127d-2,5.01074356716d-2,6.12084049407d-2, &
			7.45146949048d-2,9.00780900611d-2,1.07940155413d-1,1.27267746478d-1, &
			1.46676027814d-1,1.62254276550d-1,1.68045766353d-1,1.52383204788d-1, &
			1.01214136498d-1,-2.44389126667d-3,-1.54078468398d-1,-3.03214415655d-1, &
			-2.97674373379d-1,7.93541259524d-3,4.26273267393d-1,1.00032384844d-1, &
			-4.94117404043d-1,3.92604878741d-1,-1.90111691178d-1,7.43654896362d-2, &
			-2.78508428343d-2,1.09992061155d-2,-4.69798719697d-3,2.12587632706d-3, &
			-9.81986734159d-4,4.44992546836d-4,-1.89983519162d-4,7.31024164292d-5, &
			-2.40057837293d-5,6.23096824846d-6,-1.12363896552d-6,1.04470606055d-7/)
!	Filtro J1. 47 pontos (Guptasarma e Singh, 1997)
	s1_J147=-3.05078187595d0
	del_J147=1.10599010095d-1
	gs_abs_J147=(/s1_J147+(/(ii,ii=0,46)/)*del_J147/)
	gs_pes_J147=(/3.17926147465d-6,-9.73811660718d-6,1.64866227408d-5,-1.81501261160d-5, &
			1.87556556369d-5,-1.46550406038d-5,1.53799733803d-5,-6.95628273934d-6, &
			1.41881555665d-5,3.41445665537d-6,2.13941715512d-5,2.34962369042d-5, &
			4.84340283290d-5,7.33732978590d-5,1.27703784430d-4,2.08120025730d-4, &
			3.49803898913d-4,5.79107814687d-4,9.65887918451d-4,1.60401273703d-3, &
			2.66903777685d-3,4.43111590040d-3,7.35631696247d-3,1.21782796293d-2, &
			2.01097829218d-2,3.30096953061d-2,5.37143591532d-2,8.60516613299d-2, &
			1.34267607144d-1,2.00125033067d-1,2.74027505792d-1,3.18168749246d-1, &
			2.41655667461d-1,-5.40549161658d-2,-4.46912952135d-1,-1.92231885629d-1, &
			5.52376753950d-1,-3.57429049025d-1,1.41510519002d-1,-4.61421935309d-2, &
			1.48273761923d-2,-5.07479209193d-3,1.83829713749d-3,-6.67742804324d-4, &
			2.21277518118d-4,-5.66248732755d-5,7.88229202853d-6/)
!	Filtro J1. 140 pontos (Guptasarma e Singh, 1997)
	s1_J1140=-7.91001919d0
	del_J1140=8.79671439570d-2
	gs_abs_J1140=(/s1_J1140+(/(ii,ii=0,139)/)*del_J1140/)
	gs_pes_J1140=(/-6.76671159511d-14,3.39808396836d-13,-7.43411889153d-13,8.93613024469d-13, &
			-5.47341591896d-13,-5.84920181906d-14,5.20780672883d-13,-6.92656254606d-13, &
			6.88908045074d-13,-6.39910528298d-13,5.82098912530d-13,-4.84912700478d-13, &
			3.54684337858d-13,-2.10855291368d-13,1.00452749275d-13,5.58449957721d-15, &
			-5.67206735175d-14,1.09107856853d-13,-6.04067500756d-14,8.84512134731d-14, &
			2.22321981827d-14,8.38072239207d-14,1.23647835900d-13,1.44351787234d-13, &
			2.94276480713d-13,3.39965995918d-13,6.17024672340d-13,8.25310217692d-13, &
			1.32560792613d-12,1.90949961267d-12,2.93458179767d-12,4.33454210095d-12, &
			6.55863288798d-12,9.78324910827d-12,1.47126365223d-11,2.20240108708d-11, &
			3.30577485691d-11,4.95377381480d-11,7.43047574433d-11,1.11400535181d-10, &
			1.67052734516d-10,2.50470107577d-10,3.75597211630d-10,5.63165204681d-10, &
			8.44458166896d-10,1.26621795331d-9,1.89866561359d-9,2.84693620927d-9, &
			4.26886170263d-9,6.40104325574d-9,9.59798498616d-9,1.43918931885d-8, &
			2.15798696769d-8,3.23584600810d-8,4.85195105813d-8,7.27538583183d-8, &
			1.09090191748d-7,1.63577866557d-7,2.45275193920d-7,3.67784458730d-7, &
			5.51470341585d-7,8.26916206192d-7,1.23991037294d-6,1.85921554669d-6, &
			2.78777669034d-6,4.18019870272d-6,6.26794044911d-6,9.39858833064d-6, &
			1.40925408889d-5,2.11312291505d-5,3.16846342900d-5,4.75093313246d-5, &
			7.12354794719d-5,1.06810848460d-4,1.60146590551d-4,2.40110903628d-4, &
			3.59981158972d-4,5.39658308918d-4,8.08925141201d-4,1.21234066243d-3, &
			1.81650387595d-3,2.72068483151d-3,4.07274689463d-3,6.09135552241d-3, &
			9.09940027636d-3,1.35660714813d-2,2.01692550906d-2,2.98534800308d-2, &
			4.39060697220d-2,6.39211368217d-2,9.16763946228d-2,1.28368795114d-1, &
			1.73241920046d-1,2.19830379079d-1,2.51193131178d-1,2.32380049895d-1, &
			1.17121080205d-1,-1.17252913088d-1,-3.52148528535d-1,-2.71162871370d-1, &
			2.91134747110d-1,3.17192840623d-1,-4.93075681595d-1,3.11223091821d-1, &
			-1.36044122543d-1,5.12141261934d-2,-1.90806300761d-2,7.57044398633d-3, &
			-3.25432753751d-3,1.49774676371d-3,-7.24569558272d-4,3.62792644965d-4, &
			-1.85907973641d-4,9.67201396593d-5,-5.07744171678d-5,2.67510121456d-5, &
			-1.40667136728d-5,7.33363699547d-6,-3.75638767050d-6,1.86344211280d-6, &
			-8.71623576811d-7,3.61028200288d-7,-1.05847108097d-7,-1.51569361490d-8, &
			6.67633241420d-8,-8.33741579804d-8,8.31065906136d-8,-7.53457009758d-8, &
			6.48057680299d-8,-5.37558016587d-8,4.32436265303d-8,-3.37262648712d-8, &
			2.53558687098d-8,-1.81287021528d-8,1.20228328586d-8,-7.10898040664d-9, &
			3.53667004588d-9,-1.36030600198d-9,3.52544249042d-10,-4.53719284366d-11/)
	end subroutine J0J1GS

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

	subroutine J0J1Kong(k_abs_J0161,k_pes_J061,k_pes_J161,k_abs_J01241,k_pes_J0241,k_pes_J1241)
!	Cuidado! Os valores identificados como abscissas, na verdade já são os valores de exp(s1+(i-1)*del) ou 10**(s1+(i-1)*del)
	implicit none
	real(8),intent(out)::k_abs_J0161(61),k_pes_J061(61),k_pes_J161(61),k_abs_J01241(241),k_pes_J0241(241),k_pes_J1241(241)

    k_abs_J0161(1:3)   = (/ 2.3517745856009100D-02,  2.6649097336355482D-02,  3.0197383422318501D-02 /)
    k_abs_J0161(4:6)   = (/ 3.4218118311666032D-02,  3.8774207831722009D-02,  4.3936933623407420D-02 /)
    k_abs_J0161(7:9)   = (/ 4.9787068367863938D-02,  5.6416139503777350D-02,  6.3927861206707570D-02 /)
    k_abs_J0161(10:12) = (/ 7.2439757034251456D-02,  8.2084998623898800D-02,  9.3014489210663506D-02 /)
    k_abs_J0161(13:15) = (/ 1.0539922456186430D-01,  1.1943296826671961D-01,  1.3533528323661270D-01 /)
    k_abs_J0161(16:18) = (/ 1.5335496684492850D-01,  1.7377394345044520D-01,  1.9691167520419400D-01 /)
    k_abs_J0161(19:21) = (/ 2.2313016014842979D-01,  2.5283959580474641D-01,  2.8650479686019009D-01 /)
    k_abs_J0161(22:24) = (/ 3.2465246735834979D-01,  3.6787944117144239D-01,  4.1686201967850839D-01 /)
    k_abs_J0161(25:27) = (/ 4.7236655274101469D-01,  5.3526142851899028D-01,  6.0653065971263342D-01 /)
    k_abs_J0161(28:30) = (/ 6.8728927879097224D-01,  7.7880078307140488D-01,  8.8249690258459546D-01 /)
    k_abs_J0161(31:33) = (/ 1.0000000000000000D+00,  1.1331484530668261D+00,  1.2840254166877421D+00 /)
    k_abs_J0161(34:36) = (/ 1.4549914146182010D+00,  1.6487212707001280D+00,  1.8682459574322221D+00 /)
    k_abs_J0161(37:39) = (/ 2.1170000166126748D+00,  2.3988752939670981D+00,  2.7182818284590451D+00 /)
    k_abs_J0161(40:42) = (/ 3.0802168489180310D+00,  3.4903429574618419D+00,  3.9550767229205772D+00 /)
    k_abs_J0161(43:45) = (/ 4.4816890703380636D+00,  5.0784190371800806D+00,  5.7546026760057307D+00 /)
    k_abs_J0161(46:48) = (/ 6.5208191203301116D+00,  7.3890560989306504D+00,  8.3728974881272649D+00 /)
    k_abs_J0161(49:51) = (/ 9.4877358363585262D+00,  1.0751013186076360D+01,  1.2182493960703470D+01 /)
    k_abs_J0161(52:54) = (/ 1.3804574186067100D+01,  1.5642631884188170D+01,  1.7725424121461639D+01 /)
    k_abs_J0161(55:57) = (/ 2.0085536923187671D+01,  2.2759895093526730D+01,  2.5790339917193059D+01 /)
    k_abs_J0161(58:60) = (/ 2.9224283781234941D+01,  3.3115451958692312D+01,  3.7524723159601002D+01 /)
    k_abs_J0161(61:61) = (/ 4.2521082000062783D+01/)
    
    k_pes_J061(1:3)   = (/  1.4463210615326699D+02, -1.1066222143752420D+03,  3.7030010025325978D+03 /)
    k_pes_J061(4:6)   = (/ -6.8968188464424520D+03,  7.1663544112656937D+03, -2.4507884783377681D+03 /)
    k_pes_J061(7:9)   = (/ -4.0166567754046082D+03,  6.8623845298546094D+03, -5.0013321011775661D+03 /)
    k_pes_J061(10:12) = (/  2.1291291365196648D+03, -1.3845222435542289D+03,  2.1661554291595580D+03 /)
    k_pes_J061(13:15) = (/ -2.2260393789657141D+03,  8.0317156013986391D+02,  1.0142221718890841D+03 /)
    k_pes_J061(16:18) = (/ -1.9350455051432630D+03,  1.6601169447226580D+03, -7.5159684285420133D+02 /)
    k_pes_J061(19:21) = (/ -9.0315984178183285D+01,  5.0705574889546148D+02, -5.1207646422722519D+02 /)
    k_pes_J061(22:24) = (/  2.9722959494490038D+02, -5.0248319908072993D+01, -1.2290725861955920D+02 /)
    k_pes_J061(25:27) = (/  1.9695244755899429D+02, -1.9175679966946601D+02,  1.4211755630338590D+02 /)
    k_pes_J061(28:30) = (/ -7.7463216543224149D+01,  1.7638009334931201D+01,  2.8855056499202671D+01 /)
    k_pes_J061(31:33) = (/ -5.9225643887809561D+01,  7.5987941373668960D+01, -8.1687962781233580D+01 /)
    k_pes_J061(34:36) = (/  8.0599209238447102D+01, -7.4895905328771619D+01,  6.7516291538794434D+01 /)
    k_pes_J061(37:39) = (/ -5.9325033647358048D+01,  5.1617042242841528D+01, -4.4664967446820263D+01 /)
    k_pes_J061(40:42) = (/  3.8366152052928278D+01, -3.3308787868993100D+01,  2.8278671651033459D+01 /)
    k_pes_J061(43:45) = (/ -2.4505863388620480D+01,  2.0469632532079750D+01, -1.7074034940700429D+01 /)
    k_pes_J061(46:48) = (/  1.4206119215530070D+01, -1.0904435643084650D+01,  8.7518389425802283D+00 /)
    k_pes_J061(49:51) = (/ -6.7721665239085622D+00,  4.5096884588095891D+00, -3.2704247166629590D+00 /)
    k_pes_J061(52:54) = (/  2.6827195063720430D+00, -1.8406031821386459D+00,  9.1586697140412443D-01 /)
    k_pes_J061(55:57) = (/ -3.2436011485890798D-01,  8.0675176189581893D-02, -1.2881307195759690D-02 /)
    k_pes_J061(58:60) = (/  7.0489137468452920D-04,  2.3846917590855061D-04, -6.9102205995825531D-05 /)
    k_pes_J061(61:61) = (/  6.7792635718095777D-06 /)
    
    k_pes_J161(1:3) = (/  4.6440396425864918D+01, -4.5034239857914162D+02,  1.7723440076223640D+03 /)
    k_pes_J161(4:6) = (/ -3.7559735516994660D+03,  4.4736494009764137D+03, -2.2476603569606068D+03 /)
    k_pes_J161(7:9) = (/ -1.5219842155931799D+03,  3.4904608559273802D+03, -2.4814243247472318D+03 /)
    k_pes_J161(10:12) = (/  5.7328164634108396D+02,  5.3132044837659631D-01,  6.8895205008006235D+02 /)
    k_pes_J161(13:15) = (/ -1.2012013872160269D+03,  7.9679138423597340D+02,  4.9874460187939818D+01 /)
    k_pes_J161(16:18) = (/ -5.6367338332457007D+02,  4.7971936503711203D+02, -5.8979702298044558D+01 /)
    k_pes_J161(19:21) = (/ -3.1935800954986922D+02,  4.5762551999442371D+02, -3.7239927283248380D+02 /)
    k_pes_J161(22:24) = (/  1.8255852885279569D+02, -2.3504740340815669D-01, -1.1588151583545380D+02 /)
    k_pes_J161(25:27) = (/  1.5740956677133170D+02, -1.4334746114883359D+02,  9.9857411013284818D+01 /)
    k_pes_J161(28:30) = (/ -4.8246322019171487D+01,  2.0371404343057380D+00,  3.3003938094974323D+01 /)
    k_pes_J161(31:33) = (/ -5.5476151884197712D+01,  6.7354852323852583D+01, -7.0735403363284121D+01 /)
    k_pes_J161(34:36) = (/  6.8872932663164747D+01, -6.3272750944993042D+01,  5.6501568721817442D+01 /)
    k_pes_J161(37:39) = (/ -4.8706577819918110D+01,  4.1737211284663481D+01, -3.4776621242200903D+01 /)
    k_pes_J161(40:42) = (/  2.9161717578906430D+01, -2.3886749056000909D+01,  1.9554007583544220D+01 /)
    k_pes_J161(43:45) = (/ -1.5966397353366460D+01,  1.2429310210239199D+01, -1.0139180791868180D+01 /)
    k_pes_J161(46:48) = (/  7.4716493393871861D+00, -5.5509479014742613D+00,  4.3380799768234208D+00 /)
    k_pes_J161(49:51) = (/ -2.5911516181746550D+00,  1.6300524630626780D+00, -1.4041567266387460D+00 /)
    k_pes_J161(52:54) = (/  7.5225141726873213D-01,  4.6808777208492733D-02, -3.6630197849601159D-01 /)
    k_pes_J161(55:57) = (/  2.8948389902792782D-01, -1.3705521898064801D-01,  4.6292091649913013D-02 /)
    k_pes_J161(58:60) = (/ -1.1721281347435180D-02,  2.2002397354029149D-03, -2.8146036357227600D-04 /)
    k_pes_J161(61:61) = (/  1.8788896009128770D-05 /)
    
    k_abs_J01241(1:3)   = (/ 4.0973497897978643D-04,   4.3725238042414673D-04,    4.6661782370309847D-04 /)
    k_abs_J01241(4:6)   = (/ 4.9795542150327349D-04,   5.3139762179825294D-04,    5.6708576763830356D-04 /)
    k_abs_J01241(7:9)   = (/ 6.0517069453505309D-04,   6.4581336796593124D-04,    6.8918556369279311D-04 /)
    k_abs_J01241(10:12) = (/ 7.3547059377007120D-04,   7.8486408131093165D-04,    8.3757478728596971D-04 /)
    k_abs_J01241(13:15) = (/ 8.9382549284889267D-04,   9.5385394091834482D-04,    1.0179138409954376D-03 /)
    k_abs_J01241(16:18) = (/ 1.0862759414638579D-03,   1.1592291739045914D-03,    1.2370818742617167D-03 /)
    k_abs_J01241(19:21) = (/ 1.3201630860205026D-03,   1.4088239509056637D-03,    1.5034391929775724D-03 /)
    k_abs_J01241(22:24) = (/ 1.6044087023989032D-03,   1.7121592255655228D-03,    1.8271461687448972D-03 /)
    k_abs_J01241(25:27) = (/ 1.9498555228451206D-03,   2.0808059174495293D-03,    2.2205508127982939D-03 /)
    k_abs_J01241(28:30) = (/ 2.3696808389813539D-03,   2.5288262922292556D-03,    2.6986597988524685D-03 /)
    k_abs_J01241(31:33) = (/ 2.8798991580882404D-03,   3.0733103758703055D-03,    3.2797109023435700D-03 /)
    k_abs_J01241(34:36) = (/ 3.4999730868071656D-03,   3.7350278646880674D-03,    3.9858686921282905D-03 /)
    k_abs_J01241(37:39) = (/ 4.2535557448151254D-03,   4.5392203988006852D-03,    4.8440700122489664D-03 /)
    k_abs_J01241(40:42) = (/ 5.1693930283203291D-03,   5.5165644207607716D-03,    5.8870515052116155D-03 /)
    k_abs_J01241(43:45) = (/ 6.2824201408011177D-03,   6.7043413482289257D-03,    7.1545983723145792D-03 /)
    k_abs_J01241(46:48) = (/ 7.6350942188599617D-03,   8.1478596976799818D-03,    8.6950620057954561D-03 /)
    k_abs_J01241(49:51) = (/ 9.2790138870647437D-03,   9.9021834069673818D-03,    1.0567204383852655D-02 /)
    k_abs_J01241(52:54) = (/ 1.1276887520740558D-02,   1.2034232284723775D-02,    1.2842439584178571D-02 /)
    k_abs_J01241(55:57) = (/ 1.3704925297364945D-02,   1.4625334709594208D-02,    1.5607557919982831D-02 /)
    k_abs_J01241(58:60) = (/ 1.6655746282908664D-02,   1.7774329953659442D-02,    1.8968036612429941D-02 /)
    k_abs_J01241(61:63) = (/ 2.0241911445804381D-02,   2.1601338470175833D-02,    2.3052063287225571D-02 /)
    k_abs_J01241(64:66) = (/ 2.4600217367638302D-02,   2.6252343965687961D-02,    2.8015425774221808D-02 /)
    k_abs_J01241(67:69) = (/ 2.9896914436926308D-02,   3.1904762042607962D-02,    3.4047454734599344D-02 /)
    k_abs_J01241(70:72) = (/ 3.6334048577339996D-02,   3.8774207831722009D-02,    4.1378245800970381D-02 /)
    k_abs_J01241(73:75) = (/ 4.4157168419692860D-02,   4.7122720770327912D-02,    5.0287436723591865D-02 /)
    k_abs_J01241(76:78) = (/ 5.3664691912730114D-02,   5.7268760265467331D-02,    6.1114874332588359D-02 /)
    k_abs_J01241(79:81) = (/ 6.5219289668127525D-02,   6.9599353533269015D-02,    7.4273578214333877D-02 /)
    k_abs_J01241(82:84) = (/ 7.9261719264731550D-02,   8.4584859001564691D-02,    9.0265495609784266D-02 /)
    k_abs_J01241(85:87) = (/ 9.6327638230493035D-02,   1.0279690843528640D-01,    1.0970064851551141D-01 /)
    k_abs_J01241(88:90) = (/ 1.1706803704412637D-01,   1.2493021219858241D-01,    1.3332040336594936D-01 /)
    k_abs_J01241(91:93) = (/ 1.4227407158651356D-01,   1.5182905942943059D-01,    1.6202575093388075D-01 /)
    k_abs_J01241(94:96) = (/ 1.7290724229171636D-01,   1.8451952399298926D-01,    1.9691167520419406D-01 /)
    k_abs_J01241(97:99) = (/ 2.1013607120076472D-01,   2.2424860473053532D-01,    2.3930892224375450D-01 /)
    k_abs_J01241(100:102) = (/ 2.5538067598807768D-01,   2.7253179303401259D-01,    2.9083476236785155D-01 /)
    k_abs_J01241(103:105) = (/ 3.1036694126548503D-01,   3.3121088224198103D-01,    3.5345468195878016D-01 /)
    k_abs_J01241(106:108) = (/ 3.7719235356315689D-01,   4.0252422403363597D-01,    4.2955735821073915D-01 /)
    k_abs_J01241(109:111) = (/ 4.5840601130522352D-01,   4.8919211179633149D-01,    5.2204577676101604D-01 /)
    k_abs_J01241(112:114) = (/ 5.5710586181217392D-01,   5.9452054797019427D-01,    6.3444796794822822D-01 /)
    k_abs_J01241(115:117) = (/ 6.7705687449816465D-01,   7.2252735364207221D-01,    7.7105158580356625D-01 /)
    k_abs_J01241(118:120) = (/ 8.2283465805601841D-01,   8.7809543092056130D-01,    9.3706746337740343D-01 /)
    k_abs_J01241(121:123) = (/ 1.0000000000000000D+00,   1.0671590243841926D+00,    1.1388283833246218D+00 /)
    k_abs_J01241(124:126) = (/ 1.2153109864897309D+00,   1.2969300866657718D+00,    1.3840306459807514D+00 /)
    k_abs_J01241(127:129) = (/ 1.4769807938826425D+00,   1.5761733830339912D+00,    1.6820276496988864D+00 /)
    k_abs_J01241(130:132) = (/ 1.7949909856399000D+00,   1.9155408290138962D+00,    2.0441866822585570D+00 /)
    k_abs_J01241(133:135) = (/ 2.1814722654982011D+00,   2.3279778145702346D+00,    2.4843225333848169D+00 /)
    k_abs_J01241(136:138) = (/ 2.6511672109826070D+00,   2.8292170143515598D+00,    3.0192244688065686D+00 /)
    k_abs_J01241(139:141) = (/ 3.2219926385284996D+00,   3.4383785207051250D+00,    3.6692966676192444D+00 /)
    k_abs_J01241(142:144) = (/ 3.9157230519927220D+00,   4.1786991919232470D+00,    4.4593365528478257D+00 /)
    k_abs_J01241(145:147) = (/ 4.7588212451378542D+00,   5.0784190371800815D+00,    5.4194807051312059D+00 /)
    k_abs_J01241(148:150) = (/ 5.7834477419567749D+00,   6.1718584498835538D+00,    6.5863544420150681D+00 /)
    k_abs_J01241(151:153) = (/ 7.0286875805892945D+00,   7.5007273812029629D+00,    8.0044689142963534D+00 /)
    k_abs_J01241(154:156) = (/ 8.5420412372940930D+00,   9.1157163930403051D+00,    9.7279190125598838D+00 /)
    k_abs_J01241(157:159) = (/ 1.0381236562731843D+01,   1.1078430282186428D+01,    1.1822446851646363D+01 /)
    k_abs_J01241(160:162) = (/ 1.2616430848036902D+01,   1.3463738035001692D+01,    1.4367949545996751D+01 /)
    k_abs_J01241(163:165) = (/ 1.5332887019907195D+01,   1.6362628753157214D+01,    1.7461526936579997D+01 /)
    k_abs_J01241(166:168) = (/ 1.8634226049899006D+01,   1.9885682491564729D+01,    2.1221185526912038D+01 /)
    k_abs_J01241(169:171) = (/ 2.2646379643175401D+01,   2.4167288405845099D+01,    2.5790339917193062D+01 /)
    k_abs_J01241(172:174) = (/ 2.7522393984568446D+01,   2.9370771113289432D+01,    3.1343283446669389D+01 /)
    k_abs_J01241(175:177) = (/ 3.3448267783944921D+01,   3.5694620815655881D+01,    3.8091836725399020D+01 /)
    k_abs_J01241(178:180) = (/ 4.0650047316878776D+01,   4.3380064835851620D+01,    4.6293427667950432D+01 /)
    k_abs_J01241(181:183) = (/ 4.9402449105530195D+01,   5.2720269389647328D+01,    5.6260911247127851D+01 /)
    k_abs_J01241(184:186) = (/ 6.0039339157450577D+01,   6.4071522599936642D+01,    6.8374503548558152D+01 /)
    k_abs_J01241(187:189) = (/ 7.2966468499632811D+01,   7.7866825336828100D+01,    8.3096285358343764D+01 /)
    k_abs_J01241(190:192) = (/ 8.8676950812960627D+01,   9.4632408314924064D+01,    1.0098782853248096D+02 /)
    k_abs_J01241(193:195) = (/ 1.0777007257140046D+02,   1.1500780550310940D+02,    1.2273161751726525D+02 /)
    k_abs_J01241(196:198) = (/ 1.3097415321081860D+02,   1.3977024956000301D+02,    1.4915708315838788D+02 /)
    k_abs_J01241(199:201) = (/ 1.5917432734329714D+02,   1.6986431987468302D+02,    1.8127224187515122D+02 /)
    k_abs_J01241(202:204) = (/ 1.9344630878742183D+02,   2.0643797415630826D+02,    2.2030214709649516D+02 /)
    k_abs_J01241(205:207) = (/ 2.3509742436523857D+02,   2.5088633802084459D+02,    2.6773561971364717D+02 /)
    k_abs_J01241(208:210) = (/ 2.8571648272651305D+02,   3.0490492295690876D+02,    3.2538204011263201D+02 /)
    k_abs_J01241(211:213) = (/ 3.4723438047873475D+02,   3.7055430270433595D+02,    3.9544036815532411D+02 /)
    k_abs_J01241(214:216) = (/ 4.2199775748276141D+02,   4.5033871516762099D+02,    4.8058302392070897D+02 /)
    k_abs_J01241(217:219) = (/ 5.1285851094282907D+02,   5.4730158818487951D+02,    5.8405782889129489D+02 /)
    k_abs_J01241(220:222) = (/ 6.2328258286358414D+02,   6.6514163304436181D+02,    7.0981189619693009D+02 /)
    k_abs_J01241(223:225) = (/ 7.5748217064180938D+02,   8.0835393421053413D+02,    8.6264219578923701D+02 /)
    k_abs_J01241(226:228) = (/ 9.2057640405108020D+02,   9.8240141721825944D+02,    1.0483785379522853D+03 /)
    k_abs_J01241(229:231) = (/ 1.1187866177464875D+03,   1.1939232354884318D+03,    1.2741059551734540D+03 /)
    k_abs_J01241(232:234) = (/ 1.3596736680849924D+03,   1.4509880251144575D+03,    1.5484349652742915D+03 /)
    k_abs_J01241(235:237) = (/ 1.6524263468644833D+03,   1.7634016881866382D+03,    1.8818300251626902D+03 /)
    k_abs_J01241(238:240) = (/ 2.0082118937094979D+03,   2.1430814452477584D+03,    2.2870087042864643D+03 /)
    k_abs_J01241(241:241) = (/ 2.4406019776245007D+03 /)
      
    k_pes_J0241(1:3) = (/  2.0521734894828349D+01, -1.8902686822627982D+02,  7.3615711333621391D+02 /)
    k_pes_J0241(4:6) = (/ -1.5743636401162100D+03,  1.9044172285355407D+03, -9.3177949828687565D+02 /)
    k_pes_J0241(7:9) = (/ -8.4756674787751410D+02,  1.8896580243120454D+03, -1.4364973619770130D+03 /)
    k_pes_J0241(10:12) = (/  3.9999475225222255D+02, -1.3542454844021287D+02,  9.2532241451272546D+02 /)
    k_pes_J0241(13:15) = (/ -1.9710056537324476D+03,  2.4446279171484089D+03, -2.1906123203102729D+03 /)
    k_pes_J0241(16:18) = (/  1.5850478978775775D+03, -1.0505328230896389D+03,  7.6984613632493506D+02 /)
    k_pes_J0241(19:21) = (/ -7.0040976495766233D+02,  7.1606276337453505D+02, -7.1676220487815579D+02 /)
    k_pes_J0241(22:24) = (/  6.6234563164374151D+02, -5.5792724680287802D+02,  4.2758743543890813D+02 /)
    k_pes_J0241(25:27) = (/ -2.9556591300691770D+02,  1.7825665706981184D+02, -8.3414675513698882D+01 /)
    k_pes_J0241(28:30) = (/  1.2417017098470755D+01,  3.7033552727335611D+01, -6.8808606022416711D+01 /)
    k_pes_J0241(31:33) = (/  8.7059849275066583D+01, -9.5546538383371242D+01,  9.7381161313521631D+01 /)
    k_pes_J0241(34:36) = (/ -9.4989208296887298D+01,  9.0188673544147917D+01, -8.4292718636942382D+01 /)
    k_pes_J0241(37:39) = (/  7.8226832173348868D+01, -7.2614408878060885D+01,  6.7859682477117360D+01 /)
    k_pes_J0241(40:42) = (/ -6.4196278592205076D+01,  6.1741547025682522D+01, -6.0519141387284385D+01 /)
    k_pes_J0241(43:45) = (/  6.0491769334480402D+01, -6.1560425911383518D+01,  6.3581206383691445D+01 /)
    k_pes_J0241(46:48) = (/ -6.6354594873531482D+01,  6.9646829420226268D+01, -7.3188103112902326D+01 /)
    k_pes_J0241(49:51) = (/  7.6708283429767803D+01, -7.9940731778724427D+01,  8.2660409124105158D+01 /)
    k_pes_J0241(52:54) = (/ -8.4678217919647295D+01,  8.5870333056460041D+01, -8.6157409171557106D+01 /)
    k_pes_J0241(55:57) = (/  8.5526114160418373D+01, -8.3997300873449348D+01,  8.1646015415233720D+01 /)
    k_pes_J0241(58:60) = (/ -7.8563495286564020D+01,  7.4879819626555161D+01, -7.0721125046688840D+01 /)
    k_pes_J0241(61:63) = (/  6.6236029868394425D+01, -6.1547931592552032D+01,  5.6787391190038285D+01 /)
    k_pes_J0241(64:66) = (/ -5.2042250889790630D+01,  4.7401124144337359D+01, -4.2904131591225969D+01 /)
    k_pes_J0241(67:69) = (/  3.8598898252452933D+01, -3.4488213955995406D+01,  3.0594494176351635D+01 /)
    k_pes_J0241(70:72) = (/ -2.6898802810010366D+01,  2.3412439481509015D+01, -2.0106495314306827D+01 /)
    k_pes_J0241(73:75) = (/  1.6993408313817032D+01, -1.4047029053019282D+01,  1.1294767677459204D+01 /)
    k_pes_J0241(76:78) = (/ -8.7239710357940154D+00,  6.3824509509500222D+00, -4.2667780201601095D+00 /)
    k_pes_J0241(79:81) = (/  2.4325051404370148D+00, -8.6492951505166771D-01, -3.9467042614078562D-01 /)
    k_pes_J0241(82:84) = (/  1.3928518582065126D+00, -2.1151897072040806D+00,  2.6442660480353983D+00 /)
    k_pes_J0241(85:87) = (/ -2.9851029846728716D+00,  3.2414793692894390D+00, -3.4159478759422446D+00 /)
    k_pes_J0241(88:90) = (/  3.6114813897567299D+00, -3.8065481471148583D+00,  4.0881401855578137D+00 /)
    k_pes_J0241(91:93) = (/ -4.3995650940379534D+00,  4.8090112941999275D+00, -5.2256078193681690D+00 /)
    k_pes_J0241(94:96) = (/  5.7072204037420695D+00, -6.1390449324831593D+00,  6.5839615678399346D+00 /)
    k_pes_J0241(97:99) = (/ -6.9164266653291033D+00,  7.2188877597238257D+00, -7.3636288680068356D+00 /)
    k_pes_J0241(100:102) = (/  7.4607479507817942D+00, -7.3815887701068963D+00,  7.2656732754121194D+00 /)
    k_pes_J0241(103:105) = (/ -6.9791586417386533D+00,  6.6886849969172548D+00, -6.2474308336775373D+00 /)
    k_pes_J0241(106:108) = (/  5.8445619106015023D+00, -5.3113072582819267D+00,  4.8557012512460096D+00 /)
    k_pes_J0241(109:111) = (/ -4.2807487698944051D+00,  3.8139742889863348D+00, -3.2285888930306363D+00 /)
    k_pes_J0241(112:114) = (/  2.7762320886023009D+00, -2.2014921362632469D+00,  1.7852627321359549D+00 /)
    k_pes_J0241(115:117) = (/ -1.2441166658179819D+00,  8.9066883548526044D-01, -4.1204139236605647D-01 /)
    k_pes_J0241(118:120) = (/  1.5134230126345752D-01,  2.3623794557833316D-01, -3.8120367816755335D-01 /)
    k_pes_J0241(121:123) = (/  6.6106381680388782D-01, -6.8548466976773659D-01,  8.6073950596627435D-01 /)
    k_pes_J0241(124:126) = (/ -7.8244282877711779D-01,  8.7517842089752496D-01, -7.3021880955390983D-01 /)
    k_pes_J0241(127:129) = (/  7.7262232832472633D-01, -6.0480748494904057D-01,  6.2658579242671697D-01 /)
    k_pes_J0241(130:132) = (/ -4.7604045939428874D-01,  4.9376650403683020D-01, -3.8947933764667242D-01 /)
    k_pes_J0241(133:135) = (/  4.0212565692397595D-01, -3.6042458247357878D-01,  3.5204231781815948D-01 /)
    k_pes_J0241(136:138) = (/ -3.7930843668720610D-01,  3.2733689195982146D-01, -4.2315508593901219D-01 /)
    k_pes_J0241(139:141) = (/  3.0989441332277490D-01, -4.6647831068082313D-01,  2.9127018318261688D-01 /)
    k_pes_J0241(142:144) = (/ -4.8705387071535239D-01,  2.7698618486052751D-01, -4.6732646495222491D-01 /)
    k_pes_J0241(145:147) = (/  2.8247063840665948D-01, -3.9734195755492102D-01,  3.2057830222839734D-01 /)
    k_pes_J0241(148:150) = (/ -2.8520009926766782D-01,  3.8137726386396886D-01, -1.7293887127458049D-01 /)
    k_pes_J0241(151:153) = (/  4.1451464616516831D-01, -1.3721835010891198D-01,  3.4710931937928907D-01 /)
    k_pes_J0241(154:156) = (/ -2.3427502149669574D-01,  1.7531537381430198D-01, -3.8202053727848895D-01 /)
    k_pes_J0241(157:159) = (/  6.9544164681199228D-02, -3.4148035441124835D-01,  2.2691549707658903D-01 /)
    k_pes_J0241(160:162) = (/ -7.7194405672412439D-02,  3.9072641204282460D-01, -8.3147250278941501D-02 /)
    k_pes_J0241(163:165) = (/  9.2847198484461116D-02, -3.6217207710792204D-01,  1.0822725873639430D-02 /)
    k_pes_J0241(166:168) = (/  5.1804674714713058D-03,  3.1786543576503390D-01, -7.5374320952788773D-02 /)
    k_pes_J0241(169:171) = (/ -1.9068492965240180D-01, -8.3243368421653877D-02,  2.8203125407893237D-01 /)
    k_pes_J0241(172:174) = (/ -1.6118506235275606D-02, -2.6963245268577724D-01,  1.9777266987246947D-01 /)
    k_pes_J0241(175:177) = (/  8.3472930934463230D-02, -2.7664291725690993D-01,  2.8131050087465398D-01 /)
    k_pes_J0241(178:180) = (/ -1.7190500226396740D-01,  4.7929281626787730D-02,  4.0609103740114338D-02 /)
    k_pes_J0241(181:183) = (/ -8.8052363812875192D-02,  1.0655635530077956D-01, -1.0911039889769197D-01 /)
    k_pes_J0241(184:186) = (/  1.0429389492393636D-01, -9.6633502650203901D-02,  8.8193660973234708D-02 /)
    k_pes_J0241(187:189) = (/ -7.9818272272792340D-02,  7.1826372817129266D-02, -6.4338746209327960D-02 /)
    k_pes_J0241(190:192) = (/  5.7410843811670373D-02, -5.1077264811879450D-02,  4.5360450730759908D-02 /)
    k_pes_J0241(193:195) = (/ -4.0267658416644959D-02,  3.5786523769444523D-02, -3.1883184387833892D-02 /)
    k_pes_J0241(196:198) = (/  2.8503987928848736D-02, -2.5580391582215669D-02,  2.3035937972928848D-02 /)
    k_pes_J0241(199:201) = (/ -2.0793869607217731D-02,  1.8783956187539304D-02, -1.6947463041594817D-02 /)
    k_pes_J0241(202:204) = (/  1.5239810059902827D-02, -1.3631136080120022D-02,  1.2105408954537602D-02 /)
    k_pes_J0241(205:207) = (/ -1.0658744979296539D-02,  9.2973012802827770D-03, -8.0347421903655229D-03 /)
    k_pes_J0241(208:210) = (/  6.8891191289946248D-03, -5.8791440085348362D-03,  5.0201744968206546D-03 /)
    k_pes_J0241(211:213) = (/ -4.3205477515280673D-03,  3.7790104974130003D-03, -3.3838316346466528D-03 /)
    k_pes_J0241(214:216) = (/  3.1138040173687830D-03, -2.9408759849158255D-03,  2.8337622799129800D-03 /)
    k_pes_J0241(217:219) = (/ -2.7617063348841968D-03,  2.6976528553792770D-03, -2.6203758209448719D-03 /)
    k_pes_J0241(220:222) = (/  2.5154519816519300D-03, -2.3752432671978732D-03,  2.1981931889228849D-03 /)
    k_pes_J0241(223:225) = (/ -1.9877598903184075D-03,  1.7512378319139862D-03, -1.4986007967920582D-03 /)
    k_pes_J0241(226:228) = (/  1.2413783399553321D-03, -9.9151066433110242D-04,  7.6014790547529541D-04 /)
    k_pes_J0241(229:231) = (/ -5.5645455379022023D-04,  3.8658922048399793D-04, -2.5308541294203246D-04 /)
    k_pes_J0241(232:234) = (/  1.5481736123335846D-04, -8.7598356988947244D-05,  4.5276882194925013D-05 /)
    k_pes_J0241(235:237) = (/ -2.1044884609057091D-05,  8.6194148443616213D-06, -3.0265483432492266D-06 /)
    k_pes_J0241(238:240) = (/  8.7586378588036060D-07, -1.9631579670478718D-07,  3.0405955678339503D-08 /)
    k_pes_J0241(241:241) = (/ -2.4535953018971818D-09 /)  
     
    k_pes_J1241(1:3) = (/ -6.8036776043707992D+00,  1.2311367914708828D+02, -6.0393880274694629D+02 /)
    k_pes_J1241(4:6) = (/  1.4366739793969700D+03, -1.8485584122626103D+03,  1.0200389392810666D+03 /)
    k_pes_J1241(7:9) = (/  5.5616172944428206D+02, -1.3603216652028611D+03,  7.0167261135195611D+02 /)
    k_pes_J1241(10:12) = (/  5.3402931472384910D+02, -1.1429297666385887D+03,  9.0089651035093164D+02 /)
    k_pes_J1241(13:15) = (/ -4.1234322841018411D+02,  2.0384051761852257D+02, -2.8956552394430940D+02 /)
    k_pes_J1241(16:18) = (/  4.0377046372586426D+02, -3.6323786432883537D+02,  1.8298774315085831D+02 /)
    k_pes_J1241(19:21) = (/  1.9289907486739775D+01, -1.4580899521803201D+02,  1.6824962609150825D+02 /)
    k_pes_J1241(22:24) = (/ -1.1236235241859872D+02,  2.2347944322515357D+01,  6.4941329832499832D+01 /)
    k_pes_J1241(25:27) = (/ -1.2923900641568204D+02,  1.6503157425386775D+02, -1.7577683834104877D+02 /)
    k_pes_J1241(28:30) = (/  1.6854431620724952D+02, -1.5058653072004714D+02,  1.2774003008692418D+02 /)
    k_pes_J1241(31:33) = (/ -1.0402698302316379D+02,  8.1858707646356891D+01, -6.2444337865410809D+01 /)
    k_pes_J1241(34:36) = (/  4.6197726342228336D+01, -3.3061362052896769D+01,  2.2733224772330665D+01 /)
    k_pes_J1241(37:39) = (/ -1.4810952273048434D+01,  8.8755599393231179D+00, -4.5352205696513836D+00 /)
    k_pes_J1241(40:42) = (/  1.4449219838043461D+00,  6.8743711356557613D-01, -2.1019399621163908D+00 /)
    k_pes_J1241(43:45) = (/  2.9910240715100782D+00, -3.5052807631070255D+00,  3.7597471619400826D+00 /)
    k_pes_J1241(46:48) = (/ -3.8399078959004260D+00,  3.8072695188293357D+00, -3.7042465619688505D+00 /)
    k_pes_J1241(49:51) = (/  3.5586759844836791D+00, -3.3877703451391699D+00,  3.2017608368233237D+00 /)
    k_pes_J1241(52:54) = (/ -3.0067309393330288D+00,  2.8068645418014970D+00, -2.6056123347708602D+00 /)
    k_pes_J1241(55:57) = (/  2.4063835683466084D+00, -2.2124155853618990D+00,  2.0267333591783983D+00 /)
    k_pes_J1241(58:60) = (/ -1.8516694378892813D+00,  1.6889204130272601D+00, -1.5392291077389748D+00 /)
    k_pes_J1241(61:63) = (/  1.4027911155907105D+00, -1.2791049665324250D+00,  1.1676349053464947D+00 /)
    k_pes_J1241(64:66) = (/ -1.0675953584302726D+00,  9.7863461730665857D-01, -9.0027354444851670D-01 /)
    k_pes_J1241(67:69) = (/  8.3249709182475107D-01, -7.7482174407455295D-01,  7.2701622202684046D-01 /)
    k_pes_J1241(70:72) = (/ -6.8800213341352723D-01,  6.5698937864345264D-01, -6.3225165202126643D-01 /)
    k_pes_J1241(73:75) = (/  6.1273246350052379D-01, -5.9649977526660647D-01,  5.8279442137494575D-01 /)
    k_pes_J1241(76:78) = (/ -5.6994253522747218D-01,  5.5790401373575793D-01, -5.4545052833277607D-01 /)
    k_pes_J1241(79:81) = (/  5.3336442347695479D-01, -5.2069797861987954D-01,  5.0887851227988856D-01 /)
    k_pes_J1241(82:84) = (/ -4.9685712521558462D-01,  4.8646743314063406D-01, -4.7619096686921952D-01 /)
    k_pes_J1241(85:87) = (/  4.6816529269778301D-01, -4.6015156373071614D-01,  4.5466141540832178D-01 /)
    k_pes_J1241(88:90) = (/ -4.4857363825706975D-01,  4.4499201335812283D-01, -4.3980181069680696D-01 /)
    k_pes_J1241(91:93) = (/  4.3704907382079500D-01, -4.3150722634129624D-01,  4.2859194645669874D-01 /)
    k_pes_J1241(94:96) = (/ -4.2173773973311074D-01,  4.1817784119266610D-01, -4.0960162465920091D-01 /)
    k_pes_J1241(97:99) = (/  4.0553804034464430D-01, -3.9532815510851277D-01,  3.9136208758411489D-01 /)
    k_pes_J1241(100:102) = (/ -3.7983009113722693D-01,  3.7671762005172676D-01, -3.6406330003120191D-01 /)
    k_pes_J1241(103:105) = (/  3.6244935256730459D-01, -3.4854450330689657D-01,  3.4888232036420985D-01 /)
    k_pes_J1241(106:108) = (/ -3.3325449310423266D-01,  3.3591600342224537D-01, -3.1786716813473143D-01 /)
    k_pes_J1241(109:111) = (/  3.2333450491673743D-01, -3.0206709086775535D-01,  3.1108664637173294D-01 /)
    k_pes_J1241(112:114) = (/ -2.8574426655786522D-01,  2.9938562076378405D-01, -2.6899179718067823D-01 /)
    k_pes_J1241(115:117) = (/  2.8862849386853856D-01, -2.5197407940080729D-01,  2.7924443578937680D-01 /)
    k_pes_J1241(118:120) = (/ -2.3479292380166417D-01,  2.7158086404171844D-01, -2.1743177724854762D-01 /)
    k_pes_J1241(121:123) = (/  2.6585734716833231D-01, -1.9977886402875333D-01,  2.6215186968647208D-01 /)
    k_pes_J1241(124:126) = (/ -1.8169913895160320D-01,  2.6037686710257568D-01, -1.6314391449016330D-01 /)
    k_pes_J1241(127:129) = (/  2.6022362236599678D-01, -1.4430969552882361D-01,  2.6106112640293244D-01 /)
    k_pes_J1241(130:132) = (/ -1.2585051878275111D-01,  2.6176557155826263D-01, -1.0912157374630198D-01 /)
    k_pes_J1241(133:135) = (/  2.6047541087804182D-01, -9.6410349758274366D-02,  2.5434822704405619D-01 /)
    k_pes_J1241(136:138) = (/ -9.1058907636613021D-02,  2.3952590152331843D-01, -9.7233907394681532D-02 /)
    k_pes_J1241(139:141) = (/  2.1169127037476246D-01, -1.1886530901829440D-01,  1.6780342130647322D-01 /)
    k_pes_J1241(142:144) = (/ -1.5708700322500274D-01,  1.0962587924131692D-01, -2.0578176345887339D-01 /)
    k_pes_J1241(145:147) = (/  4.8885661268589009D-02, -2.4642892106977957D-01,  1.1194065592765716D-02 /)
    k_pes_J1241(148:150) = (/ -2.4747581882772479D-01,  3.0611266302904498D-02, -1.7906330306689583D-01 /)
    k_pes_J1241(151:153) = (/  1.2284053143305627D-01, -5.2283741701101033D-02,  2.3805370442772930D-01 /)
    k_pes_J1241(154:156) = (/  3.8699307341758679D-02,  2.4690725918834208D-01, -3.8424035688006668D-02 /)
    k_pes_J1241(157:159) = (/  7.1201119703791560D-02, -2.3889263981505893D-01, -7.7114690286621504D-02 /)
    k_pes_J1241(160:162) = (/ -2.1836754494419724D-01,  1.2315414896678636D-01,  7.5338270471735747D-02 /)
    k_pes_J1241(163:165) = (/  2.7696240054976468D-01, -6.8197423407317256D-02, -9.2154357023181005D-02 /)
    k_pes_J1241(166:168) = (/ -2.7364306064385413D-01,  1.4776915895856893D-01,  1.5391397605873250D-01 /)
    k_pes_J1241(169:171) = (/  1.1202666775507018D-01, -3.3501468763636766D-01,  5.3448876160720617D-02 /)
    k_pes_J1241(172:174) = (/  1.8824613952112992D-01, -4.6814043305072224D-03, -3.0449707752313021D-01 /)
    k_pes_J1241(175:177) = (/  3.7651614016005935D-01, -2.0743158205566595D-01, -1.1420479303107758D-02 /)
    k_pes_J1241(178:180) = (/  1.4728687125987736D-01, -1.8431611574819040D-01,  1.6270659034524407D-01 /)
    k_pes_J1241(181:183) = (/ -1.2222143749935742D-01,  8.4181157909789264D-02, -5.5281360387327068D-02 /)
    k_pes_J1241(184:186) = (/  3.5292561498256686D-02, -2.1995586173467183D-02,  1.3230592337602392D-02 /)
    k_pes_J1241(187:189) = (/ -7.4157047320195287D-03,  3.5088522519723777D-03, -8.5018538826118518D-04 /)
    k_pes_J1241(190:192) = (/ -9.7550322159493529D-04,  2.2326305963338822D-03, -3.0933635050624424D-03 /)
    k_pes_J1241(193:195) = (/  3.6729552029512069D-03, -4.0509173928643906D-03,  4.2837535729585508D-03 /)
    k_pes_J1241(196:198) = (/ -4.4126872566431202D-03,  4.4684173716055568D-03, -4.4741053242745125D-03 /)
    k_pes_J1241(199:201) = (/  4.4473069057176040D-03, -4.4012645981419867D-03,  4.3457940693551285D-03 /)
    k_pes_J1241(202:204) = (/ -4.2878952176510369D-03,  4.2321726027587606D-03, -4.1811395672186466D-03 /)
    k_pes_J1241(205:207) = (/  4.1354762911556932D-03, -4.0942937154137878D-03,  4.0554223224550864D-03 /)
    k_pes_J1241(208:210) = (/ -4.0157136359124266D-03,  3.9713297394453244D-03, -3.9180044788809047D-03 /)
    k_pes_J1241(211:213) = (/  3.8512768724064841D-03, -3.7667091650295829D-03,  3.6601079883520934D-03 /)
    k_pes_J1241(214:216) = (/ -3.5277760443730628D-03,  3.3668346001750932D-03, -3.1756548702066185D-03 /)
    k_pes_J1241(217:219) = (/  2.9543911277564146D-03, -2.7055140265840505D-03,  2.4341413354852411D-03 /)
    k_pes_J1241(220:222) = (/ -2.1479350621587340D-03,  1.8564376110443456D-03, -1.5699285271508928D-03 /)
    k_pes_J1241(223:225) = (/  1.2980859369581127D-03, -1.0488156933979571D-03,  8.2752974192331252D-04 /)
    k_pes_J1241(226:228) = (/ -6.3697315746016306D-04,  4.7751718968428873D-04, -3.4772734444372060D-04 /)
    k_pes_J1241(229:231) = (/  2.4499928240303252D-04, -1.6610663484386119D-04,  1.0758513633705277D-04 /)
    k_pes_J1241(232:234) = (/ -6.5954676827986505D-05,  3.7834534141493396D-05, -2.0025414475625305D-05 /)
    k_pes_J1241(235:237) = (/  9.6121169858568779D-06, -4.0939692455506115D-06,  1.5034842918522733D-06 /)
    k_pes_J1241(238:240) = (/ -4.5718370169166360D-07,  1.0805490301555136D-07, -1.7684735514194769D-08 /)
    k_pes_J1241(241:241) = (/ 1.5081394131412287D-09 /)	
	end subroutine J0J1Kong

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

	subroutine J0J1Key(kk_abs_J01101,kk_pes_J0101,kk_pes_J1101,kk_abs_J01201,kk_pes_J0201,kk_pes_J1201, &
				kk_abs_J01401,kk_pes_J0401,kk_pes_J1401)
!	Cuidado! Os valores identificados como abscissas, na verdade já são os valores de exp(s1+(i-1)*del) ou 10**(s1+(i-1)*del)
!	(dependendo da técnica usada pelo criador do filtro)
	real(8),intent(out)::kk_abs_J01101(101),kk_pes_J0101(101),kk_pes_J1101(101)
	real(8),intent(out)::kk_abs_J01201(201),kk_pes_J0201(201),kk_pes_J1201(201)
	real(8),intent(out)::kk_abs_J01401(401),kk_pes_J0401(401),kk_pes_J1401(401)

    kk_abs_J01101(   1:   3) = (/  3.182780796509667D-03,  3.570677233218250D-03,  4.005847942090417D-03 /)
    kk_abs_J01101(   4:   6) = (/  4.494054401183452D-03,  5.041760259690979D-03,  5.656216913953108D-03 /)
    kk_abs_J01101(   7:   9) = (/  6.345559512909110D-03,  7.118914664064660D-03,  7.986521265955502D-03 /)
    kk_abs_J01101(  10:  12) = (/  8.959866066878890D-03,  1.005183574463358D-02,  1.127688752074056D-02 /)
    kk_abs_J01101(  13:  15) = (/  1.265124056800531D-02,  1.419309074557768D-02,  1.592285150451168D-02 /)
    kk_abs_J01101(  16:  18) = (/  1.786342415331403D-02,  2.004050106168401D-02,  2.248290581673552D-02 /)
    kk_abs_J01101(  19:  21) = (/  2.522297483522721D-02,  2.829698548418685D-02,  3.174563637806794D-02 /)
    kk_abs_J01101(  22:  24) = (/  3.561458621137186D-02,  3.995505826065390D-02,  4.482451855926687D-02 /)
    kk_abs_J01101(  25:  27) = (/  5.028743672359186D-02,  5.641613950377735D-02,  6.329176835964070D-02 /)
    kk_abs_J01101(  28:  30) = (/  7.100535373963698D-02,  7.965902028589801D-02,  8.936733892175319D-02 /)
    kk_abs_J01101(  31:  33) = (/  1.002588437228037D-01,  1.124777336542896D-01,  1.261857817050387D-01 /)
    kk_abs_J01101(  34:  36) = (/  1.415644766941340D-01,  1.588174261069207D-01,  1.781730517728984D-01 /)
    kk_abs_J01101(  37:  39) = (/  1.998876140751445D-01,  2.242486047305353D-01,  2.515785530597565D-01 /)
    kk_abs_J01101(  40:  42) = (/  2.822392961405233D-01,  3.166367693790532D-01,  3.552263809249515D-01 /)
    kk_abs_J01101(  43:  45) = (/  3.985190410845142D-01,  4.470879265593564D-01,  5.015760690660555D-01 /)
    kk_abs_J01101(  46:  48) = (/  5.627048688069557D-01,  6.312836455069260D-01,  7.082203534678000D-01 /)
    kk_abs_J01101(  49:  51) = (/  7.945336025033340D-01,  8.913661439068313D-01,  1.000000000000000D+00 /)
    kk_abs_J01101(  52:  54) = (/  1.121873437571938D+00,  1.258600009929478D+00,  1.411989919667659D+00 /)
    kk_abs_J01101(  55:  57) = (/  1.584073984994482D+00,  1.777130526914038D+00,  1.993715533243082D+00 /)
    kk_abs_J01101(  58:  60) = (/  2.236696498819987D+00,  2.509290389936297D+00,  2.815106235624064D+00 /)
    kk_abs_J01101(  61:  63) = (/  3.158192909689768D+00,  3.543092736108982D+00,  3.974901627494749D+00 /)
    kk_abs_J01101(  64:  66) = (/  4.459336552847826D+00,  5.002811227833588D+00,  5.612521029693157D+00 /)
    kk_abs_J01101(  67:  69) = (/  6.296538261026657D+00,  7.063919023701211D+00,  7.924823117849490D+00 /)
    kk_abs_J01101(  70:  72) = (/  8.890648553371371D+00,  9.974182454814724D+00,  1.118977035755271D+01 /)
    kk_abs_J01101(  73:  75) = (/  1.255350613666823D+01,  1.408344508312441D+01,  1.579984294826040D+01 /)
    kk_abs_J01101(  76:  78) = (/  1.772542412146164D+01,  1.988568249156473D+01,  2.230921897527583D+01 /)
    kk_abs_J01101(  79:  81) = (/  2.502812018133782D+01,  2.807838322380105D+01,  3.150039230874794D+01 /)
    kk_abs_J01101(  82:  84) = (/  3.533945340427969D+01,  3.964639407257260D+01,  4.447823640552875D+01 /)
    kk_abs_J01101(  85:  87) = (/  4.989895197340787D+01,  5.598030878164415D+01,  6.280282144920171D+01 /)
    kk_abs_J01101(  88:  90) = (/  7.045681718843254D+01,  7.904363169956446D+01,  8.867695081296063D+01 /)
    kk_abs_J01101(  91:  93) = (/  9.948431564193386D+01,  1.116088111737080D+02,  1.252109606547652D+02 /)
    kk_abs_J01101(  94:  96) = (/  1.404708508514462D+02,  1.575905163233670D+02,  1.767966142764323D+02 /)
    kk_abs_J01101(  97:  99) = (/  1.983434254093812D+02,  2.225162204838159D+02,  2.496350371896939D+02 /)
    kk_abs_J01101( 100: 101) = (/   2.800589173104007D+02,  3.141906602856942D+02 /)

    kk_pes_J0101(   1:   3) = (/  1.761499629714979D+00, -1.261104080802587D+01,  4.589440735380221D+01 /)
    kk_pes_J0101(   4:   6) = (/ -1.137201615271146D+02,  2.171373845748458D+02, -3.427418757004401D+02 /)
    kk_pes_J0101(   7:   9) = (/  4.682970092296048D+02, -5.721810218704758D+02,  6.404705872914975D+02 /)
    kk_pes_J0101(  10:  12) = (/ -6.689903305418517D+02,  6.614840943418809D+02, -6.261458042620394D+02 /)
    kk_pes_J0101(  13:  15) = (/  5.724742839763362D+02, -5.091188391804840D+02,  4.429121852279642D+02 /)
    kk_pes_J0101(  16:  18) = (/ -3.786064826795332D+02,  3.191588258655895D+02, -2.660734942256455D+02 /)
    kk_pes_J0101(  19:  21) = (/  2.198909349916423D+02, -1.804655347160418D+02,  1.473251045730838D+02 /)
    kk_pes_J0101(  22:  24) = (/ -1.197617622162473D+02,  9.706708172663220D+01, -7.848237715099181D+01 /)
    kk_pes_J0101(  25:  27) = (/  6.337970122621439D+01, -5.112377096515730D+01,  4.125319165745542D+01 /)
    kk_pes_J0101(  28:  30) = (/ -3.328022871692455D+01,  2.690546850964645D+01, -2.176077775919447D+01 /)
    kk_pes_J0101(  31:  33) = (/  1.767900438275361D+01, -1.437298500741315D+01,  1.177821110514306D+01 /)
    kk_pes_J0101(  34:  36) = (/ -9.653214292533669D+00,  8.014188621856462D+00, -6.636716121370467D+00 /)
    kk_pes_J0101(  37:  39) = (/  5.605112419659892D+00, -4.690136885102554D+00,  4.040369173672801D+00 /)
    kk_pes_J0101(  40:  42) = (/ -3.404103887583860D+00,  2.995747642598741D+00, -2.524176278751386D+00 /)
    kk_pes_J0101(  43:  45) = (/  2.274303867444810D+00, -1.898578152750098D+00,  1.759903177266298D+00 /)
    kk_pes_J0101(  46:  48) = (/ -1.438255316175359D+00,  1.383718859899911D+00, -1.090411693791455D+00 /)
    kk_pes_J0101(  49:  51) = (/  1.103187499775871D+00, -8.232645366851894D-01,  8.897253309567960D-01 /)
    kk_pes_J0101(  52:  54) = (/ -6.181710614766990D-01,  7.214321035775552D-01, -4.660523646863802D-01 /)
    kk_pes_J0101(  55:  57) = (/  5.784691305884657D-01, -3.658095952962029D-01,  4.407748104824833D-01 /)
    kk_pes_J0101(  58:  60) = (/ -3.217298446491200D-01,  2.909042759264825D-01, -3.338253349563184D-01 /)
    kk_pes_J0101(  61:  63) = (/  1.293953529828388D-01, -3.718639185374572D-01,  8.698186744359777D-03 /)
    kk_pes_J0101(  64:  66) = (/ -3.389828291018994D-01,  4.963953772877858D-02, -1.147290346722791D-01 /)
    kk_pes_J0101(  67:  69) = (/  2.818501136331604D-01,  1.360129862128020D-01,  2.683528605579231D-01 /)
    kk_pes_J0101(  70:  72) = (/ -1.505769588191184D-01, -2.031064772010313D-01, -2.437687529694006D-01 /)
    kk_pes_J0101(  73:  75) = (/  2.876770799753536D-01,  2.011329588874316D-01, -2.275463527816749D-01 /)
    kk_pes_J0101(  76:  78) = (/ -1.749917271598988D-01,  4.306671898041470D-01, -3.992878212851492D-01 /)
    kk_pes_J0101(  79:  81) = (/  2.550550974806911D-01, -1.336748872779371D-01,  6.423008702535922D-02 /)
    kk_pes_J0101(  82:  84) = (/ -3.108808914792973D-02,  1.639821648563433D-02, -9.842274279418138D-03 /)
    kk_pes_J0101(  85:  87) = (/  6.685502920904526D-03, -4.955006305606947D-03,  3.851240446211632D-03 /)
    kk_pes_J0101(  88:  90) = (/ -3.047184325211406D-03,  2.405935290382968D-03, -1.869259018281964D-03 /)
    kk_pes_J0101(  91:  93) = (/  1.412821020240033D-03, -1.027340635931714D-03,  7.099124273512882D-04 /)
    kk_pes_J0101(  94:  96) = (/ -4.593396038963651D-04,  2.731571829548821D-04, -1.457155899557221D-04 /)
    kk_pes_J0101(  97:  99) = (/  6.748422483369748D-05, -2.590402870915282D-05,  7.676168023731497D-06 /)
    kk_pes_J0101( 100: 101) = (/  -1.549433571510781D-06,  1.584912057504873D-07 /)

    kk_pes_J1101(   1:   3) = (/  3.382421751388876D-04, -2.080089427136974D-03,  6.616516890510206D-03 /)
    kk_pes_J1101(   4:   6) = (/ -1.444442106541980D-02,  2.445603365138576D-02, -3.439020151920317D-02 /)
    kk_pes_J1101(   7:   9) = (/  4.205277853503258D-02, -4.616883971957859D-02,  4.665044563723574D-02 /)
    kk_pes_J1101(  10:  12) = (/ -4.417942425707860D-02,  3.983124989286448D-02, -3.456241855974430D-02 /)
    kk_pes_J1101(  13:  15) = (/  2.919015247012813D-02, -2.413480064473489D-02,  1.973163931127073D-02 /)
    kk_pes_J1101(  16:  18) = (/ -1.595784764119277D-02,  1.292753867241687D-02, -1.040069021237603D-02 /)
    kk_pes_J1101(  19:  21) = (/  8.496859564365926D-03, -6.857608112687733D-03,  5.734515396281575D-03 /)
    kk_pes_J1101(  22:  24) = (/ -4.635710275375116D-03,  4.034369992101164D-03, -3.217728935403965D-03 /)
    kk_pes_J1101(  25:  27) = (/  2.990802607299357D-03, -2.269398404454681D-03,  2.363714769044985D-03 /)
    kk_pes_J1101(  28:  30) = (/ -1.583785021876379D-03,  2.025371816097564D-03, -1.027952174606392D-03 /)
    kk_pes_J1101(  31:  33) = (/  1.922153595502174D-03, -5.043945522132147D-04,  2.054567905217619D-03 /)
    kk_pes_J1101(  34:  36) = (/  7.633835394740574D-05,  2.471902278652444D-03,  8.190345324336234D-04 /)
    kk_pes_J1101(  37:  39) = (/  3.277371965548541D-03,  1.870779100759767D-03,  4.642867014031185D-03 /)
    kk_pes_J1101(  40:  42) = (/  3.454727180711420D-03,  6.839165720545689D-03,  5.911273173907802D-03 /)
    kk_pes_J1101(  43:  45) = (/  1.028566504388743D-02,  9.752018620581887D-03,  1.561636209584471D-02 /)
    kk_pes_J1101(  46:  48) = (/  1.572789508971819D-02,  2.375132082749520D-02,  2.489168113861482D-02 /)
    kk_pes_J1101(  49:  51) = (/  3.593610268802203D-02,  3.858602166700044D-02,  5.363612015669736D-02 /)
    kk_pes_J1101(  52:  54) = (/  5.817318184012853D-02,  7.800216974165830D-02,  8.408281670359681D-02 /)
    kk_pes_J1101(  55:  57) = (/  1.083042686943279D-01,  1.133712057428050D-01,  1.383564865699443D-01 /)
    kk_pes_J1101(  58:  60) = (/  1.348125772346644D-01,  1.504073961323417D-01,  1.224262032649067D-01 /)
    kk_pes_J1101(  61:  63) = (/  1.106446266767018D-01,  3.772825296930285D-02, -1.353204196275998D-02 /)
    kk_pes_J1101(  64:  66) = (/ -1.266046406313744D-01, -1.730074119319237D-01, -2.327661941998823D-01 /)
    kk_pes_J1101(  67:  69) = (/ -1.312677137162940D-01, -1.104509755532312D-02,  2.294237929047981D-01 /)
    kk_pes_J1101(  70:  72) = (/  2.327988918717147D-01,  9.149367519642598D-02, -3.004964869169849D-01 /)
    kk_pes_J1101(  73:  75) = (/ -1.797040516746428D-01,  1.817308988818578D-01,  2.935346562382805D-01 /)
    kk_pes_J1101(  76:  78) = (/ -4.502465801331593D-01,  2.151483212376926D-01,  4.924585434461976D-02 /)
    kk_pes_J1101(  79:  81) = (/ -1.679842097056792D-01,  1.754997616383300D-01, -1.428351980172456D-01 /)
    kk_pes_J1101(  82:  84) = (/  1.076274808049898D-01, -8.030885482867917D-02,  6.079630210587774D-02 /)
    kk_pes_J1101(  85:  87) = (/ -4.687267027003955D-02,  3.664115351116358D-02, -2.884253032421485D-02 /)
    kk_pes_J1101(  88:  90) = (/  2.269535305564210D-02, -1.772060472957579D-02,  1.362455706705945D-02 /)
    kk_pes_J1101(  91:  93) = (/ -1.022825775681532D-02,  7.424720522409015D-03, -5.150709666353805D-03 /)
    kk_pes_J1101(  94:  96) = (/  3.365543012321504D-03, -2.033399376196629D-03,  1.108895838202336D-03 /)
    kk_pes_J1101(  97:  99) = (/ -5.283266194475448D-04,  2.099779511098949D-04, -6.484566005000160D-05 /)
    kk_pes_J1101( 100: 101) = (/   1.372984806560325D-05, -1.482645851470641D-06 /)

    kk_abs_J01201(   1:   3) = (/  6.112527611295728D-04,  6.582011330626792D-04,  7.087554594672159D-04 /)
    kk_abs_J01201(   4:   6) = (/  7.631927021868981D-04,  8.218110956199024D-04,  8.849317805958155D-04 /)
    kk_abs_J01201(   7:   9) = (/  9.529005637454616D-04,  1.026089812002297D-03,  1.104900492614410D-03 /)
    kk_abs_J01201(  10:  12) = (/  1.189764369843323D-03,  1.281146370384213D-03,  1.379547130466501D-03 /)
    kk_abs_J01201(  13:  15) = (/  1.485505738589109D-03,  1.599602688916441D-03,  1.722463061515276D-03 /)
    kk_abs_J01201(  16:  18) = (/  1.854759946855503D-03,  1.997218133335814D-03,  2.150618078036461D-03 /)
    kk_abs_J01201(  19:  21) = (/  2.315800182452862D-03,  2.493669396634629D-03,  2.685200176953821D-03 /)
    kk_abs_J01201(  22:  24) = (/  2.891441824663527D-03,  3.113524234494096D-03,  3.352664084780651D-03 /)
    kk_abs_J01201(  25:  27) = (/  3.610171503034560D-03,  3.887457243476130D-03,  4.186040415850672D-03 /)
    kk_abs_J01201(  28:  30) = (/  4.507556807870229D-03,  4.853767846875483D-03,  5.226570249814254D-03 /)
    kk_abs_J01201(  31:  33) = (/  5.628006414404065D-03,  6.060275608406676D-03,  6.525746018315052D-03 /)
    kk_abs_J01201(  34:  36) = (/  7.026967723461506D-03,  7.566686666625634D-03,  8.147859697679989D-03 /)
    kk_abs_J01201(  37:  39) = (/  8.773670772690183D-03,  9.447548397216066D-03,  1.017318440937716D-02 /)
    kk_abs_J01201(  40:  42) = (/  1.095455420558538D-02,  1.179593851975157D-02,  1.270194687528350D-02 /)
    kk_abs_J01201(  43:  45) = (/  1.367754283835671D-02,  1.472807121080859D-02,  1.585928731163164D-02 /)
    kk_abs_J01201(  46:  48) = (/  1.707738850748481D-02,  1.838904816496257D-02,  1.980145221062947D-02 /)
    kk_abs_J01201(  49:  51) = (/  2.132233849911398D-02,  2.296003920493998D-02,  2.472352647033940D-02 /)
    kk_abs_J01201(  52:  54) = (/  2.662246155912740D-02,  2.866724778592989D-02,  3.086908751073575D-02 /)
    kk_abs_J01201(  55:  57) = (/  3.324004351101861D-02,  3.579310506765532D-02,  3.854225912669247D-02 /)
    kk_abs_J01201(  58:  60) = (/  4.150256692682124D-02,  4.469024651236346D-02,  4.812276158381668D-02 /)
    kk_abs_J01201(  61:  63) = (/  5.181891717272583D-02,  5.579896266503613D-02,  6.008470273734039D-02 /)
    kk_abs_J01201(  64:  66) = (/  6.469961681378547D-02,  6.966898769808187D-02,  7.502004008532698D-02 /)
    kk_abs_J01201(  67:  69) = (/  8.078208971247929D-02,  8.698670396460390D-02,  9.366787481677047D-02 /)
    kk_abs_J01201(  70:  72) = (/  1.008622050590664D-01,  1.086091088249580D-01,  1.169510275215943D-01 /)
    kk_abs_J01201(  73:  75) = (/  1.259336623450285D-01,  1.356062246541897D-01,  1.460217055752807D-01 /)
    kk_abs_J01201(  76:  78) = (/  1.572371663136276D-01,  1.693140507634552D-01,  1.823185221282210D-01 /)
    kk_abs_J01201(  79:  81) = (/  1.963218253956815D-01,  2.114006776535105D-01,  2.276376883838127D-01 /)
    kk_abs_J01201(  82:  84) = (/  2.451218120391174D-01,  2.639488353792869D-01,  2.842219022392174D-01 /)
    kk_abs_J01201(  85:  87) = (/  3.060520786022707D-01,  3.295589610751891D-01,  3.548713320980274D-01 /)
    kk_abs_J01201(  88:  90) = (/  3.821278654786651D-01,  4.114778861171706D-01,  4.430821880821678D-01 /)
    kk_abs_J01201(  91:  93) = (/  4.771139155210344D-01,  5.137595112299984D-01,  5.532197380808739D-01 /)
    kk_abs_J01201(  94:  96) = (/  5.957107789003212D-01,  6.414654208273198D-01,  6.907343306373547D-01 /)
    kk_abs_J01201(  97:  99) = (/  7.437874280201796D-01,  8.009153643346592D-01,  8.624311149420454D-01 /)
    kk_abs_J01201( 100: 102) = (/  9.286716938412872D-01,  1.000000000000000D+00,  1.076806805496220D+00 /)
    kk_abs_J01201( 103: 105) = (/  1.159512896362974D+00,  1.248571377864284D+00,  1.344470156832053D+00 /)
    kk_abs_J01201( 106: 108) = (/  1.447734614663324D+00,  1.558930485621915D+00,  1.678666956213205D+00 /)
    kk_abs_J01201( 109: 111) = (/  1.807600002612004D+00,  1.946435984427591D+00,  2.095935514494364D+00 /)
    kk_abs_J01201( 112: 114) = (/  2.256917625888753D+00,  2.430264259001381D+00,  2.616925093246914D+00 /)
    kk_abs_J01201( 115: 117) = (/  2.817922749882108D+00,  3.034358394435675D+00,  3.267417769442918D+00 /)
    kk_abs_J01201( 118: 120) = (/  3.518377690535413D+00,  3.788613041474606D+00,  4.079604306451588D+00 /)
    kk_abs_J01201( 121: 123) = (/  4.392945680918757D+00,  4.730353805388542D+00,  5.093677170047323D+00 /)
    kk_abs_J01201( 124: 126) = (/  5.484906241707685D+00,  5.906184368579528D+00,  6.359819522601831D+00 /)
    kk_abs_J01201( 127: 129) = (/  6.848296943665372D+00,  7.374292754997836D+00,  7.940688624303139D+00 /)
    kk_abs_J01201( 130: 132) = (/  8.550587550976035D+00,  9.207330865882248D+00,  9.914516536837411D+00 /)
    kk_abs_J01201( 133: 135) = (/  1.067601888007134D+01,  1.149600978566695D+01,  1.237898157325731D+01 /)
    kk_abs_J01201( 136: 138) = (/  1.332977160319577D+01,  1.435358877803146D+01,  1.545604207947845D+01 /)
    kk_abs_J01201( 139: 141) = (/  1.664317129721834D+01,  1.792148011788406D+01,  1.929797175550276D+01 /)
    kk_abs_J01201( 142: 144) = (/  2.078018731859920D+01,  2.237624712415386D+01,  2.409489518475410D+01 /)
    kk_abs_J01201( 145: 147) = (/  2.594554711266131D+01,  2.793834170323650D+01,  3.008419648032392D+01 /)
    kk_abs_J01201( 148: 150) = (/  3.239486750789821D+01,  3.488301379565316D+01,  3.756226665137786D+01 /)
    kk_abs_J01201( 151: 153) = (/  4.044730436006738D+01,  4.355393259889750D+01,  4.689917102861648D+01 /)
    kk_abs_J01201( 154: 156) = (/  5.050134653574537D+01,  5.438019363641357D+01,  5.855696259189233D+01 /)
    kk_abs_J01201( 157: 159) = (/  6.305453582813728D+01,  6.789755329714343D+01,  7.311254746690634D+01 /)
    kk_abs_J01201( 160: 162) = (/  7.872808867953015D+01,  8.477494167382795D+01,  9.128623412992303D+01 /)
    kk_abs_J01201( 163: 165) = (/  9.829763815922249D+01,  1.058475657340557D+02,  1.139773791276396D+02 /)
    kk_abs_J01201( 166: 168) = (/  1.227316175172651D+02,  1.321582409921502D+02,  1.423088933027569D+02 /)
    kk_abs_J01201( 169: 171) = (/  1.532391847910440D+02,  1.650089970516890D+02,  1.776828109933644D+02 /)
    kk_abs_J01201( 172: 174) = (/  1.913300600973533D+02,  2.060255108088307D+02,  2.218496721447841D+02 /)
    kk_abs_J01201( 175: 177) = (/  2.388892367626087D+02,  2.572375559057747D+02,  2.769951508285525D+02 /)
    kk_abs_J01201( 178: 180) = (/  2.982702635016372D+02,  3.211794496157137D+02,  3.458482171317310D+02 /)
    kk_abs_J01201( 181: 183) = (/  3.724117138761822D+02,  4.010154679483840D+02,  4.318161849960710D+02 /)
    kk_abs_J01201( 184: 186) = (/  4.649826067271839D+02,  5.006964353612040D+02,  5.391533290846429D+02 /)
    kk_abs_J01201( 187: 189) = (/  5.805639739642864D+02,  6.251552381906738D+02,  6.731714149753277D+02 /)
    kk_abs_J01201( 190: 192) = (/  7.248755609109528D+02,  7.805509371268035D+02,  8.405025611345947D+02 /)
    kk_abs_J01201( 193: 195) = (/  9.050588778667341D+02,  9.745735590616712D+02,  1.049427440854279D+03 /)
    kk_abs_J01201( 196: 198) = (/  1.130030610186370D+03,  1.216824651467729D+03,  1.310285065796017D+03 /)
    kk_abs_J01201( 199: 201) = (/  1.410923875989213D+03,  1.519292431702289D+03,  1.635984429995926D+03 /)

    kk_pes_J0201(   1:   3) = (/  1.104702818216325D-01, -3.002860174042879D-01,  0.000000000000000D+00 /)
    kk_pes_J0201(   4:   6) = (/  9.304611998317738D-01, -1.237989456789895D+00,  0.000000000000000D+00 /)
    kk_pes_J0201(   7:   9) = (/  1.522782496923861D+00, -1.481262207187227D+00,  0.000000000000000D+00 /)
    kk_pes_J0201(  10:  12) = (/  1.200938682543965D+00, -1.042356281617017D+00,  0.000000000000000D+00 /)
    kk_pes_J0201(  13:  15) = (/  8.186482154233182D-01, -7.837979439107681D-01,  0.000000000000000D+00 /)
    kk_pes_J0201(  16:  18) = (/  1.072087935017888D+00, -2.017973482744253D+00,  2.638551914228436D+00 /)
    kk_pes_J0201(  19:  21) = (/ -2.916973810417140D+00,  2.932547717791041D+00, -2.782030176901836D+00 /)
    kk_pes_J0201(  22:  24) = (/  2.548634287940637D+00, -2.286132540834447D+00,  2.028668029433413D+00 /)
    kk_pes_J0201(  25:  27) = (/ -1.790466549754032D+00,  1.579778206432123D+00, -1.395100794535786D+00 /)
    kk_pes_J0201(  28:  30) = (/  1.236899777150533D+00, -1.099633354001948D+00,  9.829180606149369D-01 /)
    kk_pes_J0201(  31:  33) = (/ -8.808127495830721D-01,  7.939417322911940D-01, -7.166399975760631D-01 /)
    kk_pes_J0201(  34:  36) = (/  6.508052453577512D-01, -5.909087763650980D-01,  5.400346833052113D-01 /)
    kk_pes_J0201(  37:  39) = (/ -4.925062292020488D-01,  4.525012829058494D-01, -4.139282151279450D-01 /)
    kk_pes_J0201(  40:  42) = (/  3.820303825432799D-01, -3.500694920581675D-01,  3.243907465091445D-01 /)
    kk_pes_J0201(  43:  45) = (/ -2.973991010648227D-01,  2.766317890637754D-01, -2.534315625567141D-01 /)
    kk_pes_J0201(  46:  48) = (/  2.366617712595763D-01, -2.163875214294563D-01,  2.029716420595088D-01 /)
    kk_pes_J0201(  49:  51) = (/ -1.849663778885546D-01,  1.744470283910267D-01, -1.581901786171230D-01 /)
    kk_pes_J0201(  52:  54) = (/  1.502386381629129D-01, -1.352956168139171D-01,  1.296732438173242D-01 /)
    kk_pes_J0201(  55:  57) = (/ -1.156619860385481D-01,  1.121981827706036D-01, -9.877037098650596D-02 /)
    kk_pes_J0201(  58:  60) = (/  9.735307747326948D-02, -8.418324552747586D-02,  8.475347956447533D-02 /)
    kk_pes_J0201(  61:  63) = (/ -7.152829800825916D-02,  7.407480907190452D-02, -6.048336872780390D-02 /)
    kk_pes_J0201(  64:  66) = (/  6.504162604072715D-02, -5.077074526100943D-02,  5.742706238192050D-02 /)
    kk_pes_J0201(  67:  69) = (/ -4.215804851336973D-02,  5.105276831205472D-02, -3.445402816016690D-02 /)
    kk_pes_J0201(  70:  72) = (/  4.578139059809960D-02, -2.749842403920171D-02,  4.150734318146184D-02 /)
    kk_pes_J0201(  73:  75) = (/ -2.115449634951816D-02,  3.815224456264026D-02, -1.530500719575409D-02 /)
    kk_pes_J0201(  76:  78) = (/  3.566006962242232D-02, -9.844572094042197D-03,  3.398735588004622D-02 /)
    kk_pes_J0201(  79:  81) = (/ -4.669634342048067D-03,  3.309505748308261D-02,  3.248775057265341D-04 /)
    kk_pes_J0201(  82:  84) = (/  3.294860285661111D-02,  5.240451270870069D-03,  3.352191931272241D-02 /)
    kk_pes_J0201(  85:  87) = (/  1.016839266998432D-02,  3.479748010810809D-02,  1.518640694860538D-02 /)
    kk_pes_J0201(  88:  90) = (/  3.675989077470884D-02,  2.035623099573585D-02,  3.938537294158406D-02 /)
    kk_pes_J0201(  91:  93) = (/  2.571720552065411D-02,  4.262879025301092D-02,  3.127265047233792D-02 /)
    kk_pes_J0201(  94:  96) = (/  4.640586621922704D-02,  3.696789312241264D-02,  5.056664152720505D-02 /)
    kk_pes_J0201(  97:  99) = (/  4.265736230160079D-02,  5.485681544583107D-02,  4.805697373905592D-02 /)
    kk_pes_J0201( 100: 102) = (/  5.886326726197704D-02,  5.267934237797620D-02,  6.194081754102766D-02 /)
    kk_pes_J0201( 103: 105) = (/  5.575263672089661D-02,  6.312348308951768D-02,  5.613171205140025D-02 /)
    kk_pes_J0201( 106: 108) = (/  6.103877063939782D-02,  5.223073389072554D-02,  5.387340754058757D-02 /)
    kk_pes_J0201( 109: 111) = (/  4.205003959253396D-02,  3.949537291826706D-02,  2.344464896766128D-02 /)
    kk_pes_J0201( 112: 114) = (/  1.592975900866907D-02, -5.112454956546657D-03, -1.750424757311479D-02 /)
    kk_pes_J0201( 115: 117) = (/ -4.293028583599693D-02, -5.806536835195215D-02, -8.443837687273714D-02 /)
    kk_pes_J0201( 118: 120) = (/ -9.651641994648259D-02, -1.160353562902818D-01, -1.145716151033802D-01 /)
    kk_pes_J0201( 121: 123) = (/ -1.152746912280352D-01, -8.753115163764941D-02, -5.917367601312370D-02 /)
    kk_pes_J0201( 124: 126) = (/ -6.660378399098311D-04,  4.966494483361900D-02,  1.157206357429209D-01 /)
    kk_pes_J0201( 127: 129) = (/  1.448966611107768D-01,  1.588998330585965D-01,  1.031590212340854D-01 /)
    kk_pes_J0201( 130: 132) = (/  2.389465737393602D-02, -1.007600802615937D-01, -1.701081679147513D-01 /)
    kk_pes_J0201( 133: 135) = (/ -1.806596339730873D-01, -5.044220811489127D-02,  1.067558829070489D-01 /)
    kk_pes_J0201( 136: 138) = (/  2.243843205181916D-01,  1.124861192136151D-01, -1.060750025290042D-01 /)
    kk_pes_J0201( 139: 141) = (/ -2.494689576049384D-01, -2.532563907981785D-02,  2.313093162590108D-01 /)
    kk_pes_J0201( 142: 144) = (/  1.265041517428334D-01, -2.716895564501260D-01, -5.741673301698751D-02 /)
    kk_pes_J0201( 145: 147) = (/  2.828355798553848D-01, -1.151529184187685D-01, -1.895238214387882D-01 /)
    kk_pes_J0201( 148: 150) = (/  3.498981818372285D-01, -3.283594487966644D-01,  2.228689935519562D-01 /)
    kk_pes_J0201( 151: 153) = (/ -1.191882703384384D-01,  4.914136861200422D-02, -1.085046330669215D-02 /)
    kk_pes_J0201( 154: 156) = (/ -6.989673166059397D-03,  1.396225098602933D-02, -1.583247654053119D-02 /)
    kk_pes_J0201( 157: 159) = (/  1.553798887482501D-02, -1.443534141636538D-02,  1.310934487425667D-02 /)
    kk_pes_J0201( 160: 162) = (/ -1.179562561923210D-02,  1.057909167982074D-02, -9.482129155265809D-03 /)
    kk_pes_J0201( 163: 165) = (/  8.502636346417821D-03, -7.630076138978968D-03,  6.852043148244981D-03 /)
    kk_pes_J0201( 166: 168) = (/ -6.156761038065365D-03,  5.533849408416496D-03, -4.974383010461074D-03 /)
    kk_pes_J0201( 169: 171) = (/  4.470699747473056D-03, -4.016161870832697D-03,  3.604955260613827D-03 /)
    kk_pes_J0201( 172: 174) = (/ -3.231951863754572D-03,  2.892630190613161D-03, -2.583034645430386D-03 /)
    kk_pes_J0201( 175: 177) = (/  2.299750969077188D-03, -2.039879838855092D-03,  1.801000017171617D-03 /)
    kk_pes_J0201( 178: 180) = (/ -1.581121784044337D-03,  1.378637896797070D-03, -1.192281463551114D-03 /)
    kk_pes_J0201( 181: 183) = (/  1.021096380613436D-03, -8.644169540225837D-04,  7.218438967094129D-04 /)
    kk_pes_J0201( 184: 186) = (/ -5.932006323351060D-04,  4.784596555937276D-04, -3.776406369869407D-04 /)
    kk_pes_J0201( 187: 189) = (/  2.906937486164167D-04, -2.173880970439557D-04,  1.572247021446439D-04 /)
    kk_pes_J0201( 190: 192) = (/ -1.093880480994978D-04,  7.274278828090175D-05, -4.587508579334320D-05 /)
    kk_pes_J0201( 193: 195) = (/  2.717284254234203D-05, -1.493569664311702D-05,  7.502072444202401D-06 /)
    kk_pes_J0201( 196: 198) = (/ -3.374980009202976D-06,  1.323028566570774D-06, -4.342813347631881D-07 /)
    kk_pes_J0201( 199: 201) = (/  1.120529389225740D-07, -2.023679217396607D-08,  1.923133952677995D-09 /)

    kk_pes_J1201(   1:   3) = (/  1.289633927144360D-05, -4.692852957012775D-05,  5.712407500240781D-05 /)
    kk_pes_J1201(   4:   6) = (/  0.000000000000000D+00, -5.401898356504400D-05,  0.000000000000000D+00 /)
    kk_pes_J1201(   7:   9) = (/  1.163813605855986D-04, -1.341585158986474D-04,  0.000000000000000D+00 /)
    kk_pes_J1201(  10:  12) = (/  1.563529882751311D-04, -1.701932285008713D-04,  0.000000000000000D+00 /)
    kk_pes_J1201(  13:  15) = (/  2.685212722282496D-04, -5.148623386148907D-04,  6.653519984942213D-04 /)
    kk_pes_J1201(  16:  18) = (/ -7.072232255490351D-04,  6.684047762347991D-04, -5.847964888682289D-04 /)
    kk_pes_J1201(  19:  21) = (/  4.877023048974497D-04, -3.939313465485314D-04,  3.136345180982827D-04 /)
    kk_pes_J1201(  22:  24) = (/ -2.470064062519923D-04,  1.953997902093073D-04, -1.539602296158165D-04 /)
    kk_pes_J1201(  25:  27) = (/  1.234529832565310D-04, -9.817008106125012D-05,  8.058804362211881D-05 /)
    kk_pes_J1201(  28:  30) = (/ -6.451459975546781D-05,  5.456430430893415D-05, -4.344419043404446D-05 /)
    kk_pes_J1201(  31:  33) = (/  3.815523049837875D-05, -2.949393655677912D-05,  2.732187176767718D-05 /)
    kk_pes_J1201(  34:  36) = (/ -1.959842972127203D-05,  1.982911726363689D-05, -1.200904791302960D-05 /)
    kk_pes_J1201(  37:  39) = (/  1.443251606646564D-05, -5.683761857979987D-06,  1.043594555177340D-05 /)
    kk_pes_J1201(  40:  42) = (/  4.513240229146117D-08,  7.457609318727016D-06,  5.653216751529018D-06 /)
    kk_pes_J1201(  43:  45) = (/  5.308688404222969D-06,  1.152796570665543D-05,  3.932961600888582D-06 /)
    kk_pes_J1201(  46:  48) = (/  1.803352844385906D-05,  3.380683574927984D-06,  2.555716029426899D-05 /)
    kk_pes_J1201(  49:  51) = (/  3.803729181402945D-06,  3.454895491783953D-05,  5.466139820523614D-06 /)
    kk_pes_J1201(  52:  54) = (/  4.556234770039426D-05,  8.768101084525483D-06,  5.930097379335330D-05 /)
    kk_pes_J1201(  55:  57) = (/  1.428416302046287D-05,  7.667672615961499D-05,  2.281866040881824D-05 /)
    kk_pes_J1201(  58:  60) = (/  9.888436530447273D-05,  3.548288867698162D-05,  1.274995405373852D-04 /)
    kk_pes_J1201(  61:  63) = (/  5.380041779170104D-05,  1.646090070216299D-04,  7.984956415148433D-05 /)
    kk_pes_J1201(  64:  66) = (/  2.129842798716441D-04,  1.164554629053559D-04,  2.763133928114338D-04 /)
    kk_pes_J1201(  67:  69) = (/  1.674485387417859D-04,  3.595100112861982D-04,  2.380118276339850D-04 /)
    kk_pes_J1201(  70:  72) = (/  4.691253417077570D-04,  3.351463657538389D-04,  6.138970952511212D-04 /)
    kk_pes_J1201(  73:  75) = (/  4.682919161509869D-04,  8.054805703670692D-04,  6.501519986304682D-04 /)
    kk_pes_J1201(  76:  78) = (/  1.059417913885959D-03,  8.977881651498095D-04,  1.396414494631849D-03 /)
    kk_pes_J1201(  79:  81) = (/  1.234064772364374D-03,  1.844009097024346D-03,  1.689538327387063D-03 /)
    kk_pes_J1201(  82:  84) = (/  2.438744274863229D-03,  2.304895516259244D-03,  3.228952488337877D-03 /)
    kk_pes_J1201(  85:  87) = (/  3.134049181497694D-03,  4.278260457055239D-03,  4.247981356211152D-03 /)
    kk_pes_J1201(  88:  90) = (/  5.669864127830826D-03,  5.739341687299852D-03,  7.511500876566557D-03 /)
    kk_pes_J1201(  91:  93) = (/  7.727614768189288D-03,  9.940768446104526D-03,  1.036426807004136D-02 /)
    kk_pes_J1201(  94:  96) = (/  1.312987922331795D-02,  1.383652144093020D-02,  1.728788292640136D-02 /)
    kk_pes_J1201(  97:  99) = (/  1.836695924406547D-02,  2.265651442167874D-02,  2.420372071766719D-02 /)
    kk_pes_J1201( 100: 102) = (/  2.949262366900291D-02,  3.159191631258840D-02,  3.802499854417571D-02 /)
    kk_pes_J1201( 103: 105) = (/  4.071065391868849D-02,  4.836588871048378D-02,  5.155150426795294D-02 /)
    kk_pes_J1201( 106: 108) = (/  6.034835565661112D-02,  6.370530818484274D-02,  7.325374787602650D-02 /)
    kk_pes_J1201( 109: 111) = (/  7.602246037354336D-02,  8.540182884134039D-02,  8.613756797984258D-02 /)
    kk_pes_J1201( 112: 114) = (/  9.363045212001810D-02,  8.994780552031496D-02,  9.285438346694046D-02 /)
    kk_pes_J1201( 115: 117) = (/  8.138538169597652D-02,  7.625701418518774D-02,  5.331840281941071D-02 /)
    kk_pes_J1201( 118: 120) = (/  3.728410637826844D-02,  1.088918303862871D-03, -2.482185562667677D-02 /)
    kk_pes_J1201( 121: 123) = (/ -6.966985038573692D-02, -9.508370436510773D-02, -1.320457796877950D-01 /)
    kk_pes_J1201( 124: 126) = (/ -1.333188913518803D-01, -1.345679831656020D-01, -8.403605511993056D-02 /)
    kk_pes_J1201( 127: 129) = (/ -3.270153309655380D-02,  6.374700977990792D-02,  1.243720031600276D-01 /)
    kk_pes_J1201( 130: 132) = (/  1.818967705583269D-01,  1.375740652175917D-01,  5.811873304105276D-02 /)
    kk_pes_J1201( 133: 135) = (/ -1.043087758627017D-01, -1.815595823147628D-01, -1.819917129979687D-01 /)
    kk_pes_J1201( 136: 138) = (/  1.631402693428312D-02,  1.752375632474410D-01,  2.152308945759099D-01 /)
    kk_pes_J1201( 139: 141) = (/ -6.135834728025826D-02, -2.260012732339317D-01, -1.058818826360049D-01 /)
    kk_pes_J1201( 142: 144) = (/  2.751720328057493D-01,  7.282078436121274D-02, -2.453904197631155D-01 /)
    kk_pes_J1201( 145: 147) = (/ -2.981552944151773D-02,  3.295321158355015D-01, -3.395688755289581D-01 /)
    kk_pes_J1201( 148: 150) = (/  1.437347757121746D-01,  5.569137331864311D-02, -1.618975349034848D-01 /)
    kk_pes_J1201( 151: 153) = (/  1.840881618330268D-01, -1.638475650685530D-01,  1.323098979187242D-01 /)
    kk_pes_J1201( 154: 156) = (/ -1.035844311836129D-01,  8.132823321716304D-02, -6.502108525281801D-02 /)
    kk_pes_J1201( 157: 159) = (/  5.316632459343718D-02, -4.441483689628730D-02,  3.779017256454915D-02 /)
    kk_pes_J1201( 160: 162) = (/ -3.263631595413427D-02,  2.852104339616387D-02, -2.515785440611284D-02 /)
    kk_pes_J1201( 163: 165) = (/  2.235339482483744D-02, -1.997416046894760D-02,  1.792574037444691D-02 /)
    kk_pes_J1201( 166: 168) = (/ -1.613980618851269D-02,  1.456582802493779D-02, -1.316569266298338D-02 /)
    kk_pes_J1201( 169: 171) = (/  1.191013283575226D-02, -1.077631048875747D-02,  9.746151537603939D-03 /)
    kk_pes_J1201( 172: 174) = (/ -8.805179856684300D-03,  7.941688885035929D-03, -7.146145811536250D-03 /)
    kk_pes_J1201( 175: 177) = (/  6.410759357461669D-03, -5.729165010828318D-03,  5.096195694120956D-03 /)
    kk_pes_J1201( 178: 180) = (/ -4.507714364593528D-03,  3.960490050623528D-03, -3.452101711144719D-03 /)
    kk_pes_J1201( 181: 183) = (/  2.980856037463172D-03, -2.545706551739226D-03,  2.146162475389078D-03 /)
    kk_pes_J1201( 184: 186) = (/ -1.782177147988815D-03,  1.454007747414243D-03, -1.162041366288416D-03 /)
    kk_pes_J1201( 187: 189) = (/  9.065878213110874D-04, -6.876473801626548D-04,  5.046718335603190D-04 /)
    kk_pes_J1201( 190: 192) = (/ -3.563489117729183D-04,  2.404501739881472D-04, -1.537865721775565D-04 /)
    kk_pes_J1201( 193: 195) = (/  9.230834105005988D-05, -5.136252914208434D-05,  2.608318212185880D-05 /)
    kk_pes_J1201( 196: 198) = (/ -1.184499677743133D-05,  4.678587718668116D-06, -1.544064223784288D-06 /)
    kk_pes_J1201( 199: 201) = (/  3.995921803021977D-07, -7.219416019819105D-08,  6.844864105885603D-09 /)

    kk_abs_J01401(   1:   3) = (/  6.825603376334870D-08,  7.375625734276621D-08,  7.969970121723440D-08 /)
    kk_abs_J01401(   4:   6) = (/  8.612208106759941D-08,  9.306199062400369D-08,  1.005611335855232D-07 /)
    kk_abs_J01401(   7:   9) = (/  1.086645742284082D-07,  1.174210082088922D-07,  1.268830551878997D-07 /)
    kk_abs_J01401(  10:  12) = (/  1.371075750361032D-07,  1.481560095194920D-07,  1.600947515187244D-07 /)
    kk_abs_J01401(  13:  15) = (/  1.729955440010014D-07,  1.869359111419852D-07,  2.019996241884814D-07 /)
    kk_abs_J01401(  16:  18) = (/  2.182772048613791D-07,  2.358664693239236D-07,  2.548731159841640D-07 /)
    kk_abs_J01401(  19:  21) = (/  2.754113606638370D-07,  2.976046229505788D-07,  3.215862678579245D-07 /)
    kk_abs_J01401(  22:  24) = (/  3.475004072499321D-07,  3.755027658463905D-07,  4.057616170126568D-07 /)
    kk_abs_J01401(  25:  27) = (/  4.384587939575313D-07,  4.737907824157171D-07,  5.119699013810637D-07 /)
    kk_abs_J01401(  28:  30) = (/  5.532255789859379D-07,  5.978057311938045D-07,  6.459782515899233D-07 /)
    kk_abs_J01401(  31:  33) = (/  6.980326212227157D-07,  7.542816481697429D-07,  8.150633472817832D-07 /)
    kk_abs_J01401(  34:  36) = (/  8.807429714008947D-07,  9.517152062585649D-07,  1.028406542243639D-06 /)
    kk_abs_J01401(  37:  39) = (/  1.111277837292621D-06,  1.200827086303379D-06,  1.297592413714413D-06 /)
    kk_abs_J01401(  40:  42) = (/  1.402155307232813D-06,  1.515144112143249D-06,  1.637237807196195D-06 /)
    kk_abs_J01401(  43:  45) = (/  1.769170084765623D-06,  1.911733759794935D-06,  2.065785534025608D-06 /)
    kk_abs_J01401(  46:  48) = (/  2.232251144137991D-06,  2.412130924740802D-06,  2.606505819638747D-06 /)
    kk_abs_J01401(  49:  51) = (/  2.816543877501467D-06,  3.043507270968020D-06,  3.288759881366484D-06 /)
    kk_abs_J01401(  52:  54) = (/  3.553775494627146D-06,  3.840146657640719D-06,  4.149594248281688D-06 /)
    kk_abs_J01401(  55:  57) = (/  4.483977816605428D-06,  4.845306759362147D-06,  5.235752394978102D-06 /)
    kk_abs_J01401(  58:  60) = (/  5.657661011565692D-06,  6.113567966371402D-06,  6.606212921388751D-06 /)
    kk_abs_J01401(  61:  63) = (/  7.138556306690845D-06,  7.713797110415130D-06,  8.335392102304759D-06 /)
    kk_abs_J01401(  64:  66) = (/  9.007076606325931D-06,  9.732886947188948D-06,  1.051718470566020D-05 /)
    kk_abs_J01401(  67:  69) = (/  1.136468292842132D-05,  1.228047444997721D-05,  1.327006249680666D-05 /)
    kk_abs_J01401(  70:  72) = (/  1.433939375766396D-05,  1.549489411875887D-05,  1.674350727855744D-05 /)
    kk_abs_J01401(  73:  75) = (/  1.809273647424976D-05,  1.955068957062923D-05,  2.112612778233479D-05 /)
    kk_abs_J01401(  76:  78) = (/  2.282851832224015D-05,  2.466809129236741D-05,  2.665590115919790D-05 /)
    kk_abs_J01401(  79:  81) = (/  2.880389318280080D-05,  3.112497519896879D-05,  3.363309518571897D-05 /)
    kk_abs_J01401(  82:  84) = (/  3.634332508027547D-05,  3.927195135021125D-05,  4.243657286301518D-05 /)
    kk_abs_J01401(  85:  87) = (/  4.585620664220731D-05,  4.955140214551744D-05,  5.354436475185514D-05 /)
    kk_abs_J01401(  88:  90) = (/  5.785908919913498D-05,  6.252150377482026D-05,  6.755962612566159D-05 /)
    kk_abs_J01401(  91:  93) = (/  7.300373162293296D-05,  7.888653529491450D-05,  8.524338841989974D-05 /)
    kk_abs_J01401(  94:  96) = (/  9.211249096110730D-05,  9.953512112007229D-05,  1.075558833879609D-04 /)
    kk_abs_J01401(  97:  99) = (/  1.162229765854152D-04,  1.255884835016470D-04,  1.357086838732945D-04 /)
    kk_abs_J01401( 100: 102) = (/  1.466443925838176D-04,  1.584613251157513D-04,  1.712304924519203D-04 /)
    kk_abs_J01401( 103: 105) = (/  1.850286277986747D-04,  1.999386476954356D-04,  2.160501502814794D-04 /)
    kk_abs_J01401( 106: 108) = (/  2.334599537141681D-04,  2.522726779741279D-04,  2.726013735535827D-04 /)
    kk_abs_J01401( 109: 111) = (/  2.945682008058004D-04,  3.183051640380327D-04,  3.439549047593057D-04 /)
    kk_abs_J01401( 112: 114) = (/  3.716715588498817D-04,  4.016216828033581D-04,  4.339852546074165D-04 /)
    kk_abs_J01401( 115: 117) = (/  4.689567552777785D-04,  5.067463375445843D-04,  5.475810887141261D-04 /)
    kk_abs_J01401( 118: 120) = (/  5.917063952947988D-04,  6.393874175876612D-04,  6.909106831027882D-04 /)
    kk_abs_J01401( 121: 123) = (/  7.465858083766792D-04,  8.067473595375511D-04,  8.717568627991351D-04 /)
    kk_abs_J01401( 124: 126) = (/  9.420049769645596D-04,  1.017913840995439D-03,  1.099939610753318D-03 /)
    kk_abs_J01401( 127: 129) = (/  1.188575200157420D-03,  1.284353243230985D-03,  1.387849294835929D-03 /)
    kk_abs_J01401( 130: 132) = (/  1.499685289329845D-03,  1.620533277929307D-03,  1.751119467238237D-03 /)
    kk_abs_J01401( 133: 135) = (/  1.892228583209940D-03,  2.044708586766894D-03,  2.209475769415756D-03 /)
    kk_abs_J01401( 136: 138) = (/  2.387520259478373D-03,  2.579911972027180D-03,  2.787807038279695D-03 /)
    kk_abs_J01401( 139: 141) = (/  3.012454753087960D-03,  3.255205082272187D-03,  3.517516774912128D-03 /)
    kk_abs_J01401( 142: 144) = (/  3.800966129344984D-03,  4.107256465546973D-03,  4.438228360820604D-03 /)
    kk_abs_J01401( 145: 147) = (/  4.795870710296421D-03,  5.182332678714725D-03,  5.599936615308509D-03 /)
    kk_abs_J01401( 148: 150) = (/  6.051192009396498D-03,  6.538810570549064D-03,  7.065722523947535D-03 /)
    kk_abs_J01401( 151: 153) = (/  7.635094218859962D-03,  8.250347156047146D-03,  8.915178548439553D-03 /)
    kk_abs_J01401( 154: 156) = (/  9.633583538639471D-03,  1.040987920675908D-02,  1.124873051286366D-02 /)
    kk_abs_J01401( 157: 159) = (/  1.215517832991493D-02,  1.313466973567140D-02,  1.419309074557768D-02 /)
    kk_abs_J01401( 160: 162) = (/  1.533680168334325D-02,  1.657267540176125D-02,  1.790813858344623D-02 /)
    kk_abs_J01401( 163: 165) = (/  1.935121636967760D-02,  2.091058058553484D-02,  2.259560185112186D-02 /)
    kk_abs_J01401( 166: 168) = (/  2.441640589203004D-02,  2.638393438742415D-02,  2.851001072140295D-02 /)
    kk_abs_J01401( 169: 171) = (/  3.080741103275108D-02,  3.328994099003863D-02,  3.597251875342965D-02 /)
    kk_abs_J01401( 172: 174) = (/  3.887126462173843D-02,  4.200359790344555D-02,  4.538834159379689D-02 /)
    kk_abs_J01401( 175: 177) = (/  4.904583548701673D-02,  5.299805840335580D-02,  5.726876026546736D-02 /)
    kk_abs_J01401( 178: 180) = (/  6.188360481779279D-02,  6.687032384659380D-02,  7.225888382737922D-02 /)
    kk_abs_J01401( 181: 183) = (/  7.808166600115317D-02,  8.437366096160970D-02,  9.117267892259785D-02 /)
    kk_abs_J01401( 184: 186) = (/  9.851957692940828D-02,  1.064585043792528D-01,  1.150371683263324D-01 /)
    kk_abs_J01401( 187: 189) = (/  1.243071201657794D-01,  1.343240654192323D-01,  1.451481984836237D-01 /)
    kk_abs_J01401( 190: 192) = (/  1.568445643547723D-01,  1.694834494994701D-01,  1.831408042249150D-01 /)
    kk_abs_J01401( 193: 195) = (/  1.978986990836147D-01,  2.138458180564171D-01,  2.310779914773302D-01 /)
    kk_abs_J01401( 196: 198) = (/  2.496987719026136D-01,  2.698200563846869D-01,  2.915627588902595D-01 /)
    kk_abs_J01401( 199: 201) = (/  3.150575369034133D-01,  3.404455765799855D-01,  3.678794411714423D-01 /)
    kk_abs_J01401( 202: 204) = (/  3.975239878166446D-01,  4.295573582107391D-01,  4.641720491043618D-01 /)
    kk_abs_J01401( 205: 207) = (/  5.015760690660556D-01,  5.419941884591870D-01,  5.856692901447937D-01 /)
    kk_abs_J01401( 208: 210) = (/  6.328638290270813D-01,  6.838614092123558D-01,  7.389684882589442D-01 /)
    kk_abs_J01401( 211: 213) = (/  7.985162187593771D-01,  8.628624383213754D-01,  9.323938199059482D-01 /)
    kk_abs_J01401( 214: 216) = (/  1.007528195444534D+00,  1.088717066698399D+00,  1.176448318448691D+00 /)
    kk_abs_J01401( 217: 219) = (/  1.271249150321405D+00,  1.373689244865351D+00,  1.484384190920914D+00 /)
    kk_abs_J01401( 220: 222) = (/  1.603999182851514D+00,  1.733253017867395D+00,  1.872922415462686D+00 /)
    kk_abs_J01401( 223: 225) = (/  2.023846684922348D+00,  2.186932768947246D+00,  2.363160693705795D+00 /)
    kk_abs_J01401( 226: 228) = (/  2.553589458062927D+00,  2.759363397376282D+00,  2.981719060101298D+00 /)
    kk_abs_J01401( 229: 231) = (/  3.221992638528500D+00,  3.481627998306172D+00,  3.762185354999911D+00 /)
    kk_abs_J01401( 232: 234) = (/  4.065350649828702D+00,  4.392945680918757D+00,  4.746939050956379D+00 /)
    kk_abs_J01401( 235: 237) = (/  5.129457997027160D+00,  5.542801173730021D+00,  5.989452466383113D+00 /)
    kk_abs_J01401( 238: 240) = (/  6.472095917328692D+00,  6.993631855032968D+00,  7.557194322904825D+00 /)
    kk_abs_J01401( 241: 243) = (/  8.166169912567650D+00,  8.824218114758272D+00,  9.535293310146759D+00 /)
    kk_abs_J01401( 244: 246) = (/  1.030366853222556D+01,  1.113396114506531D+01,  1.203116059024153D+01 /)
    kk_abs_J01401( 247: 249) = (/  1.300065836967064D+01,  1.404828044452976D+01,  1.518032224495389D+01 /)
    kk_abs_J01401( 250: 252) = (/  1.640358650089262D+01,  1.772542412146164D+01,  1.915377836844383D+01 /)
    kk_abs_J01401( 253: 255) = (/  2.069723258938951D+01,  2.236506179715657D+01,  2.416728840584508D+01 /)
    kk_abs_J01401( 256: 258) = (/  2.611474245805797D+01,  2.821912670540861D+01,  3.049308693336139D+01 /)
    kk_abs_J01401( 259: 261) = (/  3.295028795300461D+01,  3.560549571641010D+01,  3.847466604903214D+01 /)
    kk_abs_J01401( 262: 264) = (/  4.157504053236070D+01,  4.492525011301342D+01,  4.854542706087923D+01 /)
    kk_abs_J01401( 265: 267) = (/  5.245732594909905D+01,  5.668445438288377D+01,  6.125221426275198D+01 /)
    kk_abs_J01401( 268: 270) = (/  6.618805443107456D+01,  7.152163561921920D+01,  7.728500868650357D+01 /)
    kk_abs_J01401( 271: 273) = (/  8.351280722204120D+01,  9.024245566687466D+01,  9.751439420705401D+01 /)
    kk_abs_J01401( 274: 276) = (/  1.053723217891024D+02,  1.138634587182101D+02,  1.230388304171765D+02 /)
    kk_abs_J01401( 277: 279) = (/  1.329535740512827D+02,  1.436672698616796D+02,  1.552442991983601D+02 /)
    kk_abs_J01401( 280: 282) = (/  1.677542314042285D+02,  1.812722418751512D+02,  1.958795638082189D+02 /)
    kk_abs_J01401( 283: 285) = (/  2.116639763528940D+02,  2.287203320984661D+02,  2.471511270676237D+02 /)
    kk_abs_J01401( 286: 288) = (/  2.670671166413821D+02,  2.885879811166153D+02,  3.118430448957007D+02 /)
    kk_abs_J01401( 289: 291) = (/  3.369720536300714D+02,  3.641260139877280D+02,  3.934681010910957D+02 /)
    kk_abs_J01401( 292: 294) = (/  4.251746390782468D+02,  4.594361606799342D+02,  4.964585521797127D+02 /)
    kk_abs_J01401( 295: 297) = (/  5.364642906374987D+02,  5.796937808113659D+02,  6.264067998114887D+02 /)
    kk_abs_J01401( 298: 300) = (/  6.768840581675209D+02,  7.314288866902692D+02,  7.903690592644506D+02 /)
    kk_abs_J01401( 301: 303) = (/  8.540587625261516D+02,  9.228807242613036D+02,  9.972485133152526D+02 /)
    kk_abs_J01401( 304: 306) = (/  1.077609024834175D+03,  1.164445165772805D+03,  1.258278756806354D+03 /)
    kk_abs_J01401( 307: 309) = (/  1.359673668084992D+03,  1.469239207674401D+03,  1.587633783044448D+03 /)
    kk_abs_J01401( 310: 312) = (/  1.715568857608799D+03,  1.853813226091298D+03,  2.003197634410938D+03 /)
    kk_abs_J01401( 313: 315) = (/  2.164619771847478D+03,  2.339049665486888D+03,  2.527535509363266D+03 /)
    kk_abs_J01401( 316: 318) = (/  2.731209963326043D+03,  2.951296959483918D+03,  3.189119057127294D+03 /)
    kk_abs_J01401( 319: 321) = (/  3.446105390326745D+03,  3.723800255966753D+03,  4.023872393822313D+03 /)
    kk_abs_J01401( 322: 324) = (/  4.348125014444879D+03,  4.698506635117757D+03,  5.077122788996905D+03 /)
    kk_abs_J01401( 325: 327) = (/  5.486248677800499D+03,  5.928342844080489D+03,  6.406061945236226D+03 /)
    kk_abs_J01401( 328: 330) = (/  6.922276718051178D+03,  7.480089229687659D+03,  8.082851518805124D+03 /)
    kk_abs_J01401( 331: 333) = (/  8.734185738821492D+03,  9.438005924363450D+03,  1.019854151170577D+04 /)
    kk_abs_J01401( 334: 336) = (/  1.102036275454031D+04,  1.190840818780434D+04,  1.286801430460541D+04 /)
    kk_abs_J01401( 337: 339) = (/  1.390494762457918D+04,  1.502543934638713D+04,  1.623622279258970D+04 /)
    kk_abs_J01401( 340: 342) = (/  1.754457387191111D+04,  1.895835480204384D+04,  2.048606135573400D+04 /)
    kk_abs_J01401( 343: 345) = (/  2.213687391406209D+04,  2.392071263371084D+04,  2.584829705973495D+04 /)
    kk_abs_J01401( 346: 348) = (/  2.793121054206043D+04,  3.018196984280969D+04,  3.261410035274024D+04 /)
    kk_abs_J01401( 349: 351) = (/  3.524221736879158D+04,  3.808211392115859D+04,  4.115085567766677D+04 /)
    kk_abs_J01401( 352: 354) = (/  4.446688349575319D+04,  4.805012423831561D+04,  5.192211051935027D+04 /)
    kk_abs_J01401( 355: 357) = (/  5.610611009895962D+04,  6.062726570529947D+04,  6.551274612369082D+04 /)
    kk_abs_J01401( 358: 360) = (/  7.079190946082873D+04,  7.649647956518646D+04,  8.266073666376984D+04 /)
    kk_abs_J01401( 361: 363) = (/  8.932172336080558D+04,  9.651946723626507D+04,  1.042972213818742D+05 /)
    kk_abs_J01401( 364: 366) = (/  1.127017243200503D+05,  1.217834808676890D+05,  1.315970656325818D+05 /)
    kk_abs_J01401( 367: 369) = (/  1.422014509662509D+05,  1.536603613439579D+05,  1.660426563014429D+05 /)
    kk_abs_J01401( 370: 372) = (/  1.794227442295626D+05,  1.938810295134221D+05,  2.095043957029773D+05 /)
    kk_abs_J01401( 373: 375) = (/  2.263867276186048D+05,  2.446294755291000D+05,  2.643422647924012D+05 /)
    kk_abs_J01401( 376: 378) = (/  2.856435546225249D+05,  3.086613499414067D+05,  3.335339705933576D+05 /)
    kk_abs_J01401( 379: 381) = (/  3.604108825445382D+05,  3.894535960623372D+05,  4.208366362720539D+05 /)
    kk_abs_J01401( 382: 384) = (/  4.547485919232073D+05,  4.913932486677662D+05,  5.309908136604783D+05 /)
    kk_abs_J01401( 385: 387) = (/  5.737792388402272D+05,  6.200156508443463D+05,  6.699778961486300D+05 /)
    kk_abs_J01401( 388: 390) = (/  7.239662107181761D+05,  7.823050242024113D+05,  8.453449095161883D+05 /)
    kk_abs_J01401( 391: 393) = (/  9.134646895224825D+05,  9.870737134762700D+05,  1.066614316909347D+06 /)
    kk_abs_J01401( 394: 396) = (/  1.152564479738165D+06,  1.245440698567906D+06,  1.345801090453257D+06 /)
    kk_abs_J01401( 397: 399) = (/  1.454248746767143D+06,  1.571435357331702D+06,  1.698065126589807D+06 /)
    kk_abs_J01401( 400: 401) = (/   1.834899005350438D+06,  1.982759263537569D+06 /)

    kk_pes_J0401(   1:   3) = (/  0.000000000000000D+00,  0.000000000000000D+00,  9.824913062932289D-08 /)
    kk_pes_J0401(   4:   6) = (/  0.000000000000000D+00,  0.000000000000000D+00,  0.000000000000000D+00 /)
    kk_pes_J0401(   7:   9) = (/  0.000000000000000D+00,  0.000000000000000D+00,  0.000000000000000D+00 /)
    kk_pes_J0401(  10:  12) = (/  1.587400592121677D-08,  0.000000000000000D+00,  0.000000000000000D+00 /)
    kk_pes_J0401(  13:  15) = (/  0.000000000000000D+00,  8.566005894708157D-08,  0.000000000000000D+00 /)
    kk_pes_J0401(  16:  18) = (/  0.000000000000000D+00,  0.000000000000000D+00,  4.185852284487868D-08 /)
    kk_pes_J0401(  19:  21) = (/  0.000000000000000D+00,  6.042919228413999D-08,  0.000000000000000D+00 /)
    kk_pes_J0401(  22:  24) = (/  0.000000000000000D+00,  9.223778633040080D-08,  0.000000000000000D+00 /)
    kk_pes_J0401(  25:  27) = (/  0.000000000000000D+00,  9.605835817155985D-08,  0.000000000000000D+00 /)
    kk_pes_J0401(  28:  30) = (/  7.119233736108770D-08,  0.000000000000000D+00,  1.084534194934211D-07 /)
    kk_pes_J0401(  31:  33) = (/  2.978122322309773D-08,  0.000000000000000D+00,  2.229467165084615D-07 /)
    kk_pes_J0401(  34:  36) = (/ -1.762389145190478D-07,  3.656694527170595D-07, -2.190824497653861D-07 /)
    kk_pes_J0401(  37:  39) = (/  3.591123314594687D-07, -1.330459122900614D-07,  2.693078682355359D-07 /)
    kk_pes_J0401(  40:  42) = (/  0.000000000000000D+00,  1.681753537860265D-07,  1.292951481169703D-07 /)
    kk_pes_J0401(  43:  45) = (/  8.732830659371336D-08,  2.394635732137388D-07,  3.263161048145093D-08 /)
    kk_pes_J0401(  46:  48) = (/  3.320472817306875D-07,  0.000000000000000D+00,  4.140474794253468D-07 /)
    kk_pes_J0401(  49:  51) = (/ -1.691089584764466D-08,  4.929867154878337D-07, -2.351486069903310D-08 /)
    kk_pes_J0401(  52:  54) = (/  5.748749293789135D-07, -2.296118919509815D-08,  6.634488608548862D-07 /)
    kk_pes_J0401(  55:  57) = (/ -1.578522272463715D-08,  7.602932506949016D-07, -3.950291386691670D-10 /)
    kk_pes_J0401(  58:  60) = (/  8.657561597587320D-07,  2.577126803166016D-08,  9.803339051961460D-07 /)
    kk_pes_J0401(  61:  63) = (/  6.493462588145047D-08,  1.105842139938396D-06,  1.182703157662404D-07 /)
    kk_pes_J0401(  64:  66) = (/  1.245710190589679D-06,  1.861786291659876D-07,  1.404457760420360D-06 /)
    kk_pes_J0401(  67:  69) = (/  2.690715579072592D-07,  1.587041404534934D-06,  3.680016417418645D-07 /)
    kk_pes_J0401(  70:  72) = (/  1.798553595216537D-06,  4.849140624691158D-07,  2.044272068810943D-06 /)
    kk_pes_J0401(  73:  75) = (/  6.226640605691353D-07,  2.329872825482119D-06,  7.849774436589454D-07 /)
    kk_pes_J0401(  76:  78) = (/  2.661676063084173D-06,  9.764401550293115D-07,  3.046886487275798D-06 /)
    kk_pes_J0401(  79:  81) = (/  1.202531755803571D-06,  3.493832828910101D-06,  1.469704682136747D-06 /)
    kk_pes_J0401(  82:  84) = (/  4.012202093180585D-06,  1.785529348521925D-06,  4.613260893326240D-06 /)
    kk_pes_J0401(  85:  87) = (/  2.158900926480479D-06,  5.310102517657744D-06,  2.600268757751806D-06 /)
    kk_pes_J0401(  88:  90) = (/  6.117954165362566D-06,  3.121898205359148D-06,  7.054521021164678D-06 /)
    kk_pes_J0401(  91:  93) = (/  3.738203205552756D-06,  8.140366436119160D-06,  4.466138097026229D-06 /)
    kk_pes_J0401(  94:  96) = (/  9.399370737137856D-06,  5.325632087195220D-06,  1.085928513468306D-05 /)
    kk_pes_J0401(  97:  99) = (/  6.340092174159829D-06,  1.255236713273140D-05,  7.537013084429620D-06 /)
    kk_pes_J0401( 100: 102) = (/  1.451609696381449D-05,  8.948708367122905D-06,  1.679400873721405D-05 /)
    kk_pes_J0401( 103: 105) = (/  1.061315399682074D-05,  1.943668323451985D-05,  1.257494571773497D-05 /)
    kk_pes_J0401( 106: 108) = (/  2.250293226239968D-05,  1.488640124393654D-05,  2.606118359025671D-05 /)
    kk_pes_J0401( 109: 111) = (/  1.760885939507792D-05,  3.019107498684298D-05,  2.081422548042667D-05 /)
    kk_pes_J0401( 112: 114) = (/  3.498528720118145D-05,  2.458680024846277D-05,  4.055166655188445D-05 /)
    kk_pes_J0401( 115: 117) = (/  2.902542214383600D-05,  4.701571256263205D-05,  3.424594329627977D-05 /)
    kk_pes_J0401( 118: 120) = (/  5.452351542272840D-05,  4.038409635240394D-05,  6.324519057873240D-05 /)
    kk_pes_J0401( 121: 123) = (/  4.759885125947362D-05,  7.337888648387832D-05,  5.607631517546404D-05 /)
    kk_pes_J0401( 124: 126) = (/  8.515549372910889D-05,  6.603428945751648D-05,  9.884406328190645D-05 /)
    kk_pes_J0401( 127: 129) = (/  7.772775341712903D-05,  1.147579427511257D-04,  9.145538883362232D-05 /)
    kk_pes_J0401( 130: 132) = (/  1.332619721204058D-04,  1.075670047677990D-04,  1.547811559467743D-04 /)
    kk_pes_J0401( 133: 135) = (/  1.264719996713977D-04,  1.798108274980929D-04,  1.486493783884393D-04 /)
    kk_pes_J0401( 136: 138) = (/  2.089282897741272D-04,  1.746596606473501D-04,  2.428062685411137D-04 /)
    kk_pes_J0401( 139: 141) = (/  2.051588438598003D-04,  2.822286254590445D-04,  2.409146476076457D-04 /)
    kk_pes_J0401( 142: 144) = (/  3.281087927133202D-04,  2.828254175118406D-04,  3.815112277320206D-04 /)
    kk_pes_J0401( 145: 147) = (/  3.319424735329088D-04,  4.436760084416005D-04,  3.894966102002753D-04 /)
    kk_pes_J0401( 148: 150) = (/  5.160474740163103D-04,  4.569285549031484D-04,  6.003086167616667D-04 /)
    kk_pes_J0401( 151: 153) = (/  5.359233992792353D-04,  6.984220304753900D-04,  6.284506339729007D-04 /)
    kk_pes_J0401( 154: 156) = (/  8.126770850250102D-04,  7.368119832657103D-04,  9.457437019690870D-04 /)
    kk_pes_J0401( 157: 159) = (/  8.636982308314989D-04,  1.100734372093935D-03,  1.012255913625649D-03 /)
    kk_pes_J0401( 160: 162) = (/  1.281276054621223D-03,  1.186165424832822D-03,  1.491593633007820D-03 /)
    kk_pes_J0401( 163: 165) = (/  1.389731874050673D-03,  1.736607751791357D-03,  1.627989098002781D-03 /)
    kk_pes_J0401( 166: 168) = (/  2.022051249754969D-03,  1.906816944033515D-03,  2.354607536102306D-03 /)
    kk_pes_J0401( 169: 171) = (/  2.233075556223528D-03,  2.742069502777252D-03,  2.614764820469374D-03 /)
    kk_pes_J0401( 172: 174) = (/  3.193516622456554D-03,  3.061214433236700D-03,  3.719512717457636D-03 /)
    kk_pes_J0401( 175: 177) = (/  3.583305236990679D-03,  4.332330454742194D-03,  4.193717895916075D-03 /)
    kk_pes_J0401( 178: 180) = (/  5.046212230367296D-03,  4.907200820233952D-03,  5.877674969863473D-03 /)
    kk_pes_J0401( 181: 183) = (/  5.740851469703553D-03,  6.845854941792223D-03,  6.714414341214463D-03 /)
    kk_pes_J0401( 184: 186) = (/  7.972865827753512D-03,  7.850607650897248D-03,  9.284131207419897D-03 /)
    kk_pes_J0401( 187: 189) = (/  9.175463215382743D-03,  1.080866736231174D-02,  1.071861578305127D-02 /)
    kk_pes_J0401( 190: 192) = (/  1.257928364139867D-02,  1.251344881977102D-02,  1.463261489338917D-02 /)
    kk_pes_J0401( 193: 195) = (/  1.459697067058206D-02,  1.700882559413141D-02,  1.700922269488248D-02 /)
    kk_pes_J0401( 196: 198) = (/  1.975072507757058D-02,  1.979189868224237D-02,  2.290187743495209D-02 /)
    kk_pes_J0401( 199: 201) = (/  2.298566772167304D-02,  2.650305108874451D-02,  2.662540459598260D-02 /)
    kk_pes_J0401( 202: 204) = (/  3.058599251785477D-02,  3.073210898400032D-02,  3.516298150334437D-02 /)
    kk_pes_J0401( 205: 207) = (/  3.529966137617352D-02,  4.020993030999442D-02,  4.027366645135865D-02 /)
    kk_pes_J0401( 208: 210) = (/  4.563986763942291D-02,  4.551858910431959D-02,  5.126253569381094D-02 /)
    kk_pes_J0401( 211: 213) = (/  5.076838045690198D-02,  5.672499952578218D-02,  5.555543753580153D-02 /)
    kk_pes_J0401( 214: 216) = (/  6.142864225226326D-02,  5.911480957167643D-02,  6.442218515556448D-02 /)
    kk_pes_J0401( 217: 219) = (/  6.026861996000471D-02,  6.428390136013452D-02,  5.731690882023245D-02 /)
    kk_pes_J0401( 220: 222) = (/  5.903878467303038D-02,  4.800846542566535D-02,  4.622345209535847D-02 /)
    kk_pes_J0401( 223: 225) = (/  2.975631767989811D-02,  2.332973372901844D-02,  4.075631045367046D-04 /)
    kk_pes_J0401( 226: 228) = (/ -1.097998927420109D-02, -3.996770310336364D-02, -5.439690105073641D-02 /)
    kk_pes_J0401( 229: 231) = (/ -8.586172077870871D-02, -9.710630904000793D-02, -1.221886100159675D-01 /)
    kk_pes_J0401( 232: 234) = (/ -1.181757180654116D-01, -1.226608581391981D-01, -8.819223320765970D-02 /)
    kk_pes_J0401( 235: 237) = (/ -5.976246996708556D-02,  1.003867346106077D-02,  6.143837055743897D-02 /)
    kk_pes_J0401( 238: 240) = (/  1.358126160861001D-01,  1.546019688831242D-01,  1.597454471983037D-01 /)
    kk_pes_J0401( 241: 243) = (/  7.378168765574926D-02, -2.155543820940989D-02, -1.576840287359735D-01 /)
    kk_pes_J0401( 244: 246) = (/ -1.867618104598368D-01, -1.378646829885588D-01,  6.349123295842654D-02 /)
    kk_pes_J0401( 247: 249) = (/  1.950187860384398D-01,  1.936424769502024D-01, -7.533434525407418D-02 /)
    kk_pes_J0401( 250: 252) = (/ -2.318813953601980D-01, -1.072893080046963D-01,  2.602429230035626D-01 /)
    kk_pes_J0401( 253: 255) = (/  1.145512544892857D-01, -2.400299276793784D-01, -8.781414190173367D-02 /)
    kk_pes_J0401( 256: 258) = (/  3.482215332524013D-01, -2.648837943057846D-01,  1.177899549516786D-02 /)
    kk_pes_J0401( 259: 261) = (/  1.809751035116563D-01, -2.458202901044509D-01,  2.258610969476037D-01 /)
    kk_pes_J0401( 262: 264) = (/ -1.773155121294782D-01,  1.311747401438518D-01, -9.665134793648206D-02 /)
    kk_pes_J0401( 265: 267) = (/  7.308596000248034D-02, -5.738192263561636D-02,  4.678168146095997D-02 /)
    kk_pes_J0401( 268: 270) = (/ -3.939779158727641D-02,  3.405264070131274D-02, -3.003070542957877D-02 /)
    kk_pes_J0401( 271: 273) = (/  2.689468371509007D-02, -2.437148059023670D-02,  2.228571230862690D-02 /)
    kk_pes_J0401( 274: 276) = (/ -2.052141594646854D-02,  1.899976167605663D-02, -1.766579040557837D-02 /)
    kk_pes_J0401( 277: 279) = (/  1.648030895118207D-02, -1.541479785692305D-02,  1.444812372540324D-02 /)
    kk_pes_J0401( 280: 282) = (/ -1.356436078983561D-02,  1.275131275416127D-02, -1.199948916264833D-02 /)
    kk_pes_J0401( 283: 285) = (/  1.130138518277636D-02, -1.065096922990164D-02,  1.004331590462835D-02 /)
    kk_pes_J0401( 286: 288) = (/ -9.474341780324707D-03,  8.940614242462378D-03, -8.439212039530420D-03 /)
    kk_pes_J0401( 289: 291) = (/  7.967622109909870D-03, -7.523661440170574D-03,  7.105415845769664D-03 /)
    kk_pes_J0401( 292: 294) = (/ -6.711190184736050D-03,  6.339466716099865D-03, -5.988869875713372D-03 /)
    kk_pes_J0401( 295: 297) = (/  5.658136700741345D-03, -5.346092874868108D-03,  5.051634994011519D-03 /)
    kk_pes_J0401( 298: 300) = (/ -4.773719694137324D-03,  4.511359360882715D-03, -4.263622706076311D-03 /)
    kk_pes_J0401( 301: 303) = (/  4.029637517014100D-03, -3.808592909460376D-03,  3.599739234584605D-03 /)
    kk_pes_J0401( 304: 306) = (/ -3.402384890213290D-03,  3.215890296536024D-03, -3.039660023908146D-03 /)
    kk_pes_J0401( 307: 309) = (/  2.873134442084575D-03, -2.715782351627641D-03,  2.567095892455205D-03 /)
    kk_pes_J0401( 310: 312) = (/ -2.426588488831973D-03,  2.293795626711503D-03, -2.168277194783892D-03 /)
    kk_pes_J0401( 313: 315) = (/  2.049619594802782D-03, -1.937436245806754D-03,  1.831366196380340D-03 /)
    kk_pes_J0401( 316: 318) = (/ -1.731071555824182D-03,  1.636234772852003D-03, -1.546556467549061D-03 /)
    kk_pes_J0401( 319: 321) = (/  1.461754019926524D-03, -1.381560792781408D-03,  1.305725738551140D-03 /)
    kk_pes_J0401( 322: 324) = (/ -1.234013089394381D-03,  1.166201823788569D-03, -1.102084687325930D-03 /)
    kk_pes_J0401( 325: 327) = (/  1.041466730622445D-03, -9.841635538957642D-04,  9.299996161562829D-04 /)
    kk_pes_J0401( 328: 330) = (/ -8.788069763395656D-04,  8.304246517808708D-04, -7.846985037183290D-04 /)
    kk_pes_J0401( 331: 333) = (/  7.414813684955915D-04, -7.006331617591514D-04,  6.620208464608574D-04 /)
    kk_pes_J0401( 334: 336) = (/ -6.255183200732592D-04,  5.910063209799918D-04, -5.583723812133077D-04 /)
    kk_pes_J0401( 337: 339) = (/  5.275107580030028D-04, -4.983222495996592D-04,  4.707138540343669D-04 /)
    kk_pes_J0401( 340: 342) = (/ -4.445983047168279D-04,  4.198935530394733D-04, -3.965222590204658D-04 /)
    kk_pes_J0401( 343: 345) = (/  3.744113329795447D-04, -3.534915670228068D-04,  3.336973877578390D-04 /)
    kk_pes_J0401( 346: 348) = (/ -3.149667279395557D-04,  2.972409648555553D-04, -2.804648423393001D-04 /)
    kk_pes_J0401( 349: 351) = (/  2.645863035783382D-04, -2.495562056240159D-04,  2.353279391742063D-04 /)
    kk_pes_J0401( 352: 354) = (/ -2.218570165506096D-04,  2.091007054227963D-04, -1.970177752442844D-04 /)
    kk_pes_J0401( 355: 357) = (/  1.855683950070826D-04, -1.747141837505830D-04,  1.644183743725741D-04 /)
    kk_pes_J0401( 358: 360) = (/ -1.546460137084759D-04,  1.453641062486393D-04, -1.365416355278548D-04 /)
    kk_pes_J0401( 361: 363) = (/  1.281494594613588D-04, -1.201601331637577D-04,  1.125477274134751D-04 /)
    kk_pes_J0401( 364: 366) = (/ -1.052876878881750D-04,  9.835675748306749D-05, -9.173298412621861D-05 /)
    kk_pes_J0401( 367: 369) = (/  8.539584183369774D-05, -7.932647157887555D-05,  7.350799507002284D-05 /)
    kk_pes_J0401( 370: 372) = (/ -6.792579804522193D-05,  6.256766460994959D-05, -5.742369885284320D-05 /)
    kk_pes_J0401( 373: 375) = (/  5.248607676312167D-05, -4.774876610278635D-05,  4.320737034558452D-05 /)
    kk_pes_J0401( 376: 378) = (/ -3.885918466176108D-05,  3.470344555454732D-05, -3.074167696093530D-05 /)
    kk_pes_J0401( 379: 381) = (/  2.697802016509551D-05, -2.341946443717172D-05,  2.007591620962982D-05 /)
    kk_pes_J0401( 382: 384) = (/ -1.696002345446848D-05,  1.408663272782194D-05, -1.147175200077509D-05 /)
    kk_pes_J0401( 385: 387) = (/  9.130950853154427D-06, -7.077245367140833D-06,  5.318671603956403D-06 /)
    kk_pes_J0401( 388: 390) = (/ -3.855913768030704D-06,  2.680452779662483D-06, -1.773679649172923D-06 /)
    kk_pes_J0401( 391: 393) = (/  1.107284379526082D-06, -6.450375038180866D-07,  3.458519700889004D-07 /)
    kk_pes_J0401( 394: 396) = (/ -1.677370250259620D-07,  7.195813987460663D-08, -2.650936383269536D-08 /)
    kk_pes_J0401( 397: 399) = (/  8.053026529021296D-09, -1.901548822685746D-09,  3.176144531356312D-10 /)
    kk_pes_J0401( 400: 401) = (/  -3.147399480514734D-11,  1.181734681811153D-12 /)
    kk_pes_J1401(   1:   3) = (/  0.000000000000000D+00,  0.000000000000000D+00,  0.000000000000000D+00 /)
    kk_pes_J1401(   4:   6) = (/ -1.761153377678605D-10,  0.000000000000000D+00,  0.000000000000000D+00 /)
    kk_pes_J1401(   7:   9) = (/  0.000000000000000D+00,  0.000000000000000D+00,  0.000000000000000D+00 /)
    kk_pes_J1401(  10:  12) = (/  0.000000000000000D+00,  3.440216255617151D-10,  0.000000000000000D+00 /)
    kk_pes_J1401(  13:  15) = (/  0.000000000000000D+00,  0.000000000000000D+00, -5.012130889133619D-10 /)
    kk_pes_J1401(  16:  18) = (/  0.000000000000000D+00,  0.000000000000000D+00,  5.096249946359556D-10 /)
    kk_pes_J1401(  19:  21) = (/  0.000000000000000D+00,  0.000000000000000D+00, -5.137927435894773D-10 /)
    kk_pes_J1401(  22:  24) = (/  0.000000000000000D+00,  4.901252024472790D-10,  0.000000000000000D+00 /)
    kk_pes_J1401(  25:  27) = (/ -3.412499448385044D-10,  0.000000000000000D+00,  2.331695173457025D-10 /)
    kk_pes_J1401(  28:  30) = (/  0.000000000000000D+00, -2.490134046395816D-10,  1.954047749927753D-10 /)
    kk_pes_J1401(  31:  33) = (/  0.000000000000000D+00, -8.103005437723991D-11,  0.000000000000000D+00 /)
    kk_pes_J1401(  34:  36) = (/  1.279088845538237D-10, -1.851480627564174D-10,  1.332990147779435D-10 /)
    kk_pes_J1401(  37:  39) = (/  0.000000000000000D+00, -1.624963592448070D-10,  3.098241711700410D-10 /)
    kk_pes_J1401(  40:  42) = (/ -4.149316491270542D-10,  4.685626300171698D-10, -4.712327624043219D-10 /)
    kk_pes_J1401(  43:  45) = (/  4.305062986549322D-10, -3.542442130313303D-10,  2.522911525991727D-10 /)
    kk_pes_J1401(  46:  48) = (/ -1.313966839434664D-10,  0.000000000000000D+00,  1.374783013515889D-10 /)
    kk_pes_J1401(  49:  51) = (/ -2.736921677833286D-10,  4.061582817302266D-10, -5.279722941941932D-10 /)
    kk_pes_J1401(  52:  54) = (/  6.386355968571155D-10, -7.315427954629247D-10,  8.087932062353638D-10 /)
    kk_pes_J1401(  55:  57) = (/ -8.640922282673907D-10,  9.028429079409167D-10, -9.186974584267662D-10 /)
    kk_pes_J1401(  58:  60) = (/  9.208119531836431D-10, -9.017596435374981D-10,  8.746151932449153D-10 /)
    kk_pes_J1401(  61:  63) = (/ -8.291880695503223D-10,  7.826366656707373D-10, -7.199240415978938D-10 /)
    kk_pes_J1401(  64:  66) = (/  6.630127031117425D-10, -5.898849811863955D-10,  5.290634200040055D-10 /)
    kk_pes_J1401(  67:  69) = (/ -4.493240826470479D-10,  3.886261259968761D-10, -3.037437366057724D-10 /)
    kk_pes_J1401(  70:  72) = (/  2.460094930171363D-10, -1.562353890563393D-10,  1.042595244297130D-10 /)
    kk_pes_J1401(  73:  75) = (/ -9.277742086714826D-12, -3.359801793948191D-11,  1.346759112142224D-10 /)
    kk_pes_J1401(  76:  78) = (/ -1.644581113778544D-10,  2.735615970185790D-10, -2.856750577185324D-10 /)
    kk_pes_J1401(  79:  81) = (/  4.064174138388151D-10, -3.956626260990240D-10,  5.339967379460912D-10 /)
    kk_pes_J1401(  82:  84) = (/ -4.942289770979441D-10,  6.590721699414525D-10, -5.825238873978475D-10 /)
    kk_pes_J1401(  85:  87) = (/  7.864186241426929D-10, -6.626363292105279D-10,  9.225771435176782D-10 /)
    kk_pes_J1401(  88:  90) = (/ -7.369977848939469D-10,  1.075704953982971D-09, -8.079491312165849D-10 /)
    kk_pes_J1401(  91:  93) = (/  1.255853770105183D-09, -8.775982413533787D-10,  1.475677833588153D-09 /)
    kk_pes_J1401(  94:  96) = (/ -9.478165158795805D-10,  1.751454564451092D-09, -1.020150878388706D-09 /)
    kk_pes_J1401(  97:  99) = (/  2.104270247990516D-09, -1.095430340362651D-09,  2.561356147221457D-09 /)
    kk_pes_J1401( 100: 102) = (/ -1.172968828002544D-09,  3.157686036019299D-09, -1.249340900707534D-09 /)
    kk_pes_J1401( 103: 105) = (/  3.938151654389299D-09, -1.316804538118081D-09,  4.960584143866221D-09 /)
    kk_pes_J1401( 106: 108) = (/ -1.361050365193285D-09,  6.299594474808514D-09, -1.357732636276717D-09 /)
    kk_pes_J1401( 109: 111) = (/  8.051372567815041D-09, -1.267357938643338D-09,  1.034000380570114D-08 /)
    kk_pes_J1401( 112: 114) = (/ -1.028333203364199D-09,  1.332639788605693D-08, -5.481128505786461D-10 /)
    kk_pes_J1401( 115: 117) = (/  1.722132657425943D-08,  3.079572462647345D-10,  2.230389235352557D-08 /)
    kk_pes_J1401( 118: 120) = (/  1.731612035551629D-09,  2.894640331467214D-08,  3.992876253106220D-09 /)
    kk_pes_J1401( 121: 123) = (/  3.764713519770196D-08,  7.470342606008511D-09,  4.907355253852564D-08 /)
    kk_pes_J1401( 124: 126) = (/  1.269390990642245D-08,  6.411984689966026D-08,  2.040424455749986D-08 /)
    kk_pes_J1401( 127: 129) = (/  8.398456134765634D-08,  3.163423793932673D-08,  1.102763405195791D-07 /)
    kk_pes_J1401( 130: 132) = (/  4.782016944498328D-08,  1.451576997806275D-07,  7.095396726449770D-08 /)
    kk_pes_J1401( 133: 135) = (/  1.915405053682995D-07,  1.037907010433357D-07,  2.533530424909100D-07 /)
    kk_pes_J1401( 136: 138) = (/  1.501316966404849D-07,  3.359015461416407D-07,  2.152163676698347D-07 /)
    kk_pes_J1401( 139: 141) = (/  4.463548364795692D-07,  3.062649070955455D-07,  5.943990986479385D-07 /)
    kk_pes_J1401( 142: 144) = (/  4.332201278554725D-07,  7.931323845043337D-07,  6.097551748617194D-07 /)
    kk_pes_J1401( 145: 147) = (/  1.060288415924516D-06,  8.546425300312143D-07,  1.419911983644722D-06 /)
    kk_pes_J1401( 148: 150) = (/  1.193611725125921D-06,  1.904650353313021D-06,  1.661883236417129D-06 /)
    kk_pes_J1401( 151: 153) = (/  2.558861699766762D-06,  2.307654632812588D-06,  3.442811157821135D-06 /)
    kk_pes_J1401( 154: 156) = (/  3.196895656510743D-06,  4.638358992076794D-06,  4.419902312928682D-06 /)
    kk_pes_J1401( 157: 159) = (/  6.256715142637839D-06,  6.100218940817017D-06,  8.449032402277294D-06 /)
    kk_pes_J1401( 160: 162) = (/  8.406776377975507D-06,  1.142087099410692D-05,  1.157042201253195D-05 /)
    kk_pes_J1401( 163: 165) = (/  1.545191683298927D-05,  1.590647868284089D-05,  2.092279678518197D-05 /)
    kk_pes_J1401( 166: 168) = (/  2.184558773386808D-05,  2.835151124366015D-05,  2.997585693259872D-05 /)
    kk_pes_J1401( 169: 171) = (/  3.844299184771377D-05,  4.110034139799833D-05,  5.215663114685968D-05 /)
    kk_pes_J1401( 172: 174) = (/  5.631528820209043D-05,  7.079837414439344D-05,  7.711659421807641D-05 /)
    kk_pes_J1401( 175: 177) = (/  9.614618225427727D-05,  1.055447661767286D-04,  1.306206751127902D-04 /)
    kk_pes_J1401( 178: 180) = (/  1.443823581310115D-04,  1.775170452103382D-04,  1.974225796898412D-04 /)
    kk_pes_J1401( 181: 183) = (/  2.413201291349723D-04,  2.698342899158978D-04,  3.281318482161766D-04 /)
    kk_pes_J1401( 184: 186) = (/  3.686575352518240D-04,  4.462498902719190D-04,  5.034750551777089D-04 /)
    kk_pes_J1401( 187: 189) = (/  6.069499176331384D-04,  6.873187128242988D-04,  8.255410973454176D-04 /)
    kk_pes_J1401( 190: 192) = (/  9.378875323431232D-04,  1.122784104067260D-03,  1.279178262961256D-03 /)
    kk_pes_J1401( 193: 195) = (/  1.526781204453462D-03,  1.743655859738250D-03,  2.075473273444583D-03 /)
    kk_pes_J1401( 196: 198) = (/  2.375110960283700D-03,  2.819903510936917D-03,  3.232361399449634D-03 /)
    kk_pes_J1401( 199: 201) = (/  3.828410021593598D-03,  4.393943043291972D-03,  5.191863182780206D-03 /)
    kk_pes_J1401( 202: 204) = (/  5.963857250198409D-03,  7.029925656135096D-03,  8.078221794986575D-03 /)
    kk_pes_J1401( 205: 207) = (/  9.497983356558860D-03,  1.091218414985025D-02,  1.279369144280816D-02 /)
    kk_pes_J1401( 208: 210) = (/  1.468546702846075D-02,  1.716070017133079D-02,  1.966301954392782D-02 /)
    kk_pes_J1401( 211: 213) = (/  2.288457156694697D-02,  2.614381914454830D-02,  3.027139461924632D-02 /)
    kk_pes_J1401( 214: 216) = (/  3.442503605732608D-02,  3.959217858802792D-02,  4.471956265829849D-02 /)
    kk_pes_J1401( 217: 219) = (/  5.096501947760914D-02,  5.699213341078931D-02,  6.413328819589428D-02 /)
    kk_pes_J1401( 220: 222) = (/  7.066611680136001D-02,  7.808858155732867D-02,  8.415253392588408D-02 /)
    kk_pes_J1401( 223: 225) = (/  9.050430770970036D-02,  9.420119432022898D-02,  9.704525747507635D-02 /)
    kk_pes_J1401( 226: 228) = (/  9.525207334318243D-02,  9.091107970583592D-02,  7.941214687854993D-02 /)
    kk_pes_J1401( 229: 231) = (/  6.360896948482617D-02,  3.851999787085855D-02,  8.930119389229475D-03 /)
    kk_pes_J1401( 232: 234) = (/ -2.931232114631221D-02, -6.745472857764939D-02, -1.071387129122266D-01 /)
    kk_pes_J1401( 235: 237) = (/ -1.337787672525563D-01, -1.467487065654186D-01, -1.280178925639679D-01 /)
    kk_pes_J1401( 238: 240) = (/ -8.266114319901335D-02, -1.829018947010989D-03,  8.538653035086596D-02 /)
    kk_pes_J1401( 241: 243) = (/  1.634948424549636D-01,  1.768283682900049D-01,  1.165478617534897D-01 /)
    kk_pes_J1401( 244: 246) = (/ -2.909864536607605D-02, -1.650763395456089D-01, -2.108189429067250D-01 /)
    kk_pes_J1401( 247: 249) = (/ -6.548403763878953D-02,  1.474231756411919D-01,  2.387075880611132D-01 /)
    kk_pes_J1401( 250: 252) = (/  8.364392500807047D-03, -2.375301184090862D-01, -1.324045834042435D-01 /)
    kk_pes_J1401( 253: 255) = (/  2.623684401711879D-01,  9.558257428232823D-02, -2.795911533969242D-01 /)
    kk_pes_J1401( 256: 258) = (/  4.155824322494574D-02,  2.715606356255319D-01, -3.742520323983440D-01 /)
    kk_pes_J1401( 259: 261) = (/  2.853293391721109D-01, -1.417457299548064D-01,  3.117732200130105D-02 /)
    kk_pes_J1401( 262: 264) = (/  2.889458199693983D-02, -5.278084724642508D-02,  5.771601780075224D-02 /)
    kk_pes_J1401( 265: 267) = (/ -5.477531833561626D-02,  4.939546060797920D-02, -4.378947822887341D-02 /)
    kk_pes_J1401( 268: 270) = (/  3.869238569366965D-02, -3.425716943762586D-02,  3.044093551672354D-02 /)
    kk_pes_J1401( 271: 273) = (/ -2.715226194815648D-02,  2.430080341386173D-02, -2.181052920166688D-02 /)
    kk_pes_J1401( 274: 276) = (/  1.962073026812987D-02, -1.768361271316385D-02,  1.596144663779621D-02 /)
    kk_pes_J1401( 277: 279) = (/ -1.442412760579961D-02,  1.304729403368628D-02, -1.181093160074335D-02 /)
    kk_pes_J1401( 280: 282) = (/  1.069835428430446D-02, -9.695462302234538D-03,  8.790199228345667D-03 /)
    kk_pes_J1401( 283: 285) = (/ -7.972150985921739D-03,  7.232245468946166D-03, -6.562523603657329D-03 /)
    kk_pes_J1401( 286: 288) = (/  5.955961653845067D-03, -5.406331000039078D-03,  4.908085823543291D-03 /)
    kk_pes_J1401( 289: 291) = (/ -4.456271717308975D-03,  4.046450048794231D-03, -3.674634444757500D-03 /)
    kk_pes_J1401( 292: 294) = (/  3.337237040467819D-03, -3.031022902436713D-03,  2.753071298464685D-03 /)
    kk_pes_J1401( 295: 297) = (/ -2.500742536477930D-03,  2.271649199484839D-03, -2.063630842184971D-03 /)
    kk_pes_J1401( 298: 300) = (/  1.874731507609227D-03, -1.703179664969160D-03,  1.547370305159833D-03 /)
    kk_pes_J1401( 301: 303) = (/ -1.405848965175527D-03,  1.277297436350464D-03, -1.160520897100170D-03 /)
    kk_pes_J1401( 304: 306) = (/  1.054436232968270D-03, -9.580613717996333D-04,  8.705055444165467D-04 /)
    kk_pes_J1401( 307: 309) = (/ -7.909604401097324D-04,  7.186922398324291D-04, -6.530344915884751D-04 /)
    kk_pes_J1401( 310: 312) = (/  5.933817682686722D-04, -5.391840290971849D-04,  4.899415886335732D-04 /)
    kk_pes_J1401( 313: 315) = (/ -4.452005847911536D-04,  4.045488440443307D-04, -3.676120749678280D-04 /)
    kk_pes_J1401( 316: 318) = (/  3.340503647933679D-04, -3.035549812866991D-04,  2.758454813737592D-04 /)
    kk_pes_J1401( 319: 321) = (/ -2.506671088369273D-04,  2.277884460323346D-04, -2.069992810457240D-04 /)
    kk_pes_J1401( 322: 324) = (/  1.881086606885530D-04, -1.709431124465889D-04,  1.553450267830473D-04 /)
    kk_pes_J1401( 325: 327) = (/ -1.411711924592810D-04,  1.282914747613330D-04, -1.165876238497402D-04 /)
    kk_pes_J1401( 328: 330) = (/  1.059521991589964D-04, -9.628759515338238D-05,  8.750515384295626D-05 /)
    kk_pes_J1401( 331: 333) = (/ -7.952435126636077D-05,  7.227204913633232D-05, -6.568180805443277D-05 /)
    kk_pes_J1401( 334: 336) = (/  5.969326263717584D-05, -5.425155895475095D-05,  4.930685077009869D-05 /)
    kk_pes_J1401( 337: 339) = (/ -4.481384639404849D-05,  4.073139627769274D-05, -3.702211350786090D-05 /)
    kk_pes_J1401( 340: 342) = (/  3.365202295484529D-05, -3.059023779502484D-05,  2.780866410092890D-05 /)
    kk_pes_J1401( 343: 345) = (/ -2.528173545505838D-05,  2.298617982129317D-05, -2.090081950209839D-05 /)
    kk_pes_J1401( 346: 348) = (/  1.900640153317426D-05, -1.728545097076178D-05,  1.572213560023288D-05 /)
    kk_pes_J1401( 349: 351) = (/ -1.430213090870262D-05,  1.301247986370528D-05, -1.184145017024185D-05 /)
    kk_pes_J1401( 352: 354) = (/  1.077839734116676D-05, -9.813642618247219D-06,  8.938371729153144D-06 /)
    kk_pes_J1401( 355: 357) = (/ -8.144555804633981D-06,  7.424890780058305D-06, -6.772748023670692D-06 /)
    kk_pes_J1401( 358: 360) = (/  6.182129388552279D-06, -5.647624656180153D-06,  5.164374046764935D-06 /)
    kk_pes_J1401( 361: 363) = (/ -4.728037195336446D-06,  4.334761689817321D-06, -3.981136442356169D-06 /)
    kk_pes_J1401( 364: 366) = (/  3.664118637205433D-06, -3.380939384592361D-06,  3.129009529723152D-06 /)
    kk_pes_J1401( 367: 369) = (/ -2.905847492954934D-06,  2.709034327844308D-06, -2.536183084429876D-06 /)
    kk_pes_J1401( 370: 372) = (/  2.384905895111026D-06, -2.252772087223033D-06,  2.137259365710216D-06 /)
    kk_pes_J1401( 373: 375) = (/ -2.035698487395405D-06,  1.945208220753021D-06, -1.862627051737732D-06 /)
    kk_pes_J1401( 376: 378) = (/  1.784472814922380D-06, -1.706983978589542D-06,  1.626297101583022D-06 /)
    kk_pes_J1401( 379: 381) = (/ -1.538789746921996D-06,  1.441575114294373D-06, -1.333077139855850D-06 /)
    kk_pes_J1401( 382: 384) = (/  1.213536734497111D-06, -1.085223736705093D-06,  9.521451230126823D-07 /)
    kk_pes_J1401( 385: 387) = (/ -8.192296419664399D-07,  6.912547372495382D-07, -5.719422069867315D-07 /)
    kk_pes_J1401( 388: 390) = (/  4.635461014353267D-07, -3.669835113430816D-07,  2.823229600705464D-07 /)
    kk_pes_J1401( 391: 393) = (/ -2.093541854410099D-07,  1.479865821668501D-07, -9.831343971854163D-08 /)
    kk_pes_J1401( 394: 396) = (/  6.033173947686275D-08, -3.349527235750482D-08,  1.640107598073682D-08 /)
    kk_pes_J1401( 397: 399) = (/ -6.856420858188184D-09,  2.340536901536035D-09, -6.097093750212437D-10 /)
    kk_pes_J1401( 400: 401) = (/   1.074334228952679D-10, -9.576603341064843D-12 /)	

	end subroutine J0J1Key

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

end module mdlgm1_1D
