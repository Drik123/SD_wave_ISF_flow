program spread_depolarize

  integer, parameter :: DP = SELECTED_REAL_KIND(14)
  integer, parameter :: Nx = 26
  integer, parameter :: Ny = 26
  integer, parameter :: Nx_fin = 101
  integer, parameter :: Ny_fin = 101
  integer, parameter :: Nt = 600001
  integer, parameter :: Nt_fin = 1000

  real(KIND=DP) :: dx,dy,dt,dt2,RR,Farad,Vi,f,S,temp,phi,Imax,D_K,D_Na,eps,Cm,Vtot
  real(KIND=DP) :: gK_ghk_dr,gK_ghk_a,gk_hh,gNa_ghk_T,gNa_ghk_P,gNa_hh,ghh
  real(KIND=DP) :: pi,alpha0,freq,lambda,kwave,omega,beta,Am,Ap,Ae,Aw,An,As,weight
  real(KIND=DP) :: res1,Ng,toler,r,lp,L,mdss,kh,delP,T,Vscale,tort,Diffp,Pe,Pe2,Le3,D,utotss
  integer :: it,i,j,k,p,q,file_step,last_step,itinit,itmax,isConverged1
  real(KIND=DP), dimension(Nx) :: x
  real(KIND=DP), dimension(Ny) :: y
  real(KIND=DP), dimension(Ny,Nx) :: K_e,K_i,Na_e,Na_i,Em,K_tot,E_K,E_Na
  real(KIND=DP), dimension(Ny,Nx) :: I_K,I_leak,I_pump,I_K_dr,I_K_a,I_K_leak,I_K_pump
  real(KIND=DP), dimension(Ny,Nx) :: I_Na_T,I_Na_P,I_Na_leak,I_Na_pump,I_Na,I_tot
  real(KIND=DP), dimension(Ny,Nx) :: mkdr,hkdr,alphamkdr,betamkdr,alphahkdr,betahkdr
  real(KIND=DP), dimension(Ny,Nx) :: mka,hka,alphamka,betamka,alphahka,betahka
  real(KIND=DP), dimension(Ny,Nx) :: mnat,hnat,alphamnat,betamnat,alphahnat,betahnat
  real(KIND=DP), dimension(Ny,Nx) :: mnap,hnap,alphamnap,betamnap,alphahnap,betahnap
  real(KIND=DP), dimension(Ny,Nx) :: alpha, dalphadt
  real(KIND=DP), dimension(Ny_fin,Nx_fin) :: alpha_fin_old,alpha_fin,dalphadt_fin,perm
  real(KIND=DP), dimension(Nt) :: tt
  real(KIND=DP), dimension(Nt_fin)::tw
  real(KIND=DP), dimension(14,Ny,Nx):: z, dzdt
  
  real(KIND=DP), dimension(Nx_fin) :: x_fin
  real(KIND=DP), dimension(Ny_fin) :: y_fin,mds,utots
  real(KIND=DP), dimension(Ny_fin,Nx_fin)::PP,Pnew,u,v,xm,ym,diff1,md,utot
  real(KIND=DP), dimension(Ny_fin,Nx_fin)::S1
  real(KIND=DP), dimension(Nx_fin*Ny_fin)::Pread,Pread2
  
  EXTERNAL rk4
  EXTERNAl SystemDerives

  real :: start, finish
  CHARACTER(40) numchr

  file_step = 20000
  last_step = 1000*file_step ! previously was 1800, changed to 1k

  dx=0.001
  dt=5e-3!dt=2.5e-4 !in ms (Changed from 5e-3 to 2.5e-4)

  x(1)=0.0
  y(1)=0.0
  tt(1)=0.0

  do i=2,Nx
  x(i)=x(i-1)+dx
  end do

  do i=2,Ny
  y(i)=y(i-1)+dx
  end do

  do i=2,Nt
  tt(i)=tt(i-1)+dt
  end do

!potassium goldman-hodgkins-katz conductance
gK_ghk_dr = 100e-5 !in S/cm^2
gK_ghk_a = 10e-5
!sodium goldman-hodgkins-katz conductance
gNa_ghk_T = 100e-5
gNa_ghk_P = 2e-5
!sodium hodgkin-huxley condectance
gNa_hh = 2e-5
!sodium hodgkin-huxley condectance
gk_hh = 7e-5
!general leak conductance
ghh = 20e-5

Farad = 96.485 !Coloumb/mM
RR = 8.51 !universal gas constant mv Coloumb/mM K

Vi = 2160e-12 !cm^3
S = 1586e-8 !cm^2
f = 0.15
temp = 298 ! K
phi = (RR*temp)/Farad !mV
!diffusion coeffs
D_K = 1.96e-5 !cm^2/s
D_Na = 1.33e-5 !cm^2/s

Imax = 0.013 !in mA/(cm^2)

Cm = 7.5e-7 !membrane capacitance in s/\Omega cm^2

eps = 0.5e-3 !in cm

!initiate e - extracellular, i - intracellular

!define the system matrix
do j=1,Ny
do i=1,Nx
z(1,j,i) = 141.46
z(2,j,i) = 3.86 + 150.0*exp(-((x(i)+y(j)-0.04)**2)/(2*(eps**2)))
z(3,j,i) = 9.82
z(4,j,i) = 133.45
z(5,j,i) = -70 !in mV
!gates
z(6,j,i) = 0.12 !mka
z(7,j,i) = 0.12 !hka
z(8,j,i) = 1.3e-3 !mkdr
z(9,j,i) = 0.1 !doesn't matter, hkdr
z(10,j,i) = 5e-3 !mnat
z(11,j,i) = 0.995 !hnat
z(12,j,i) = 1.3e-2 !mnap
z(13,j,i) = 0.1 !hnap
z(14,j,i) = 0.2 !alpha
end do
end do

do j=1,Ny
do i=1,Nx
Na_e(j,i) = z(1,j,i)
K_e(j,i) = z(2,j,i) !in mM
Na_i(j,i) = z(3,j,i)
K_i(j,i) = z(4,j,i)
Em(j,i) = z(5,j,i) !in mV
!gates
mka(j,i) = z(6,j,i)
hka(j,i) = z(7,j,i)
mkdr(j,i) = z(8,j,i)
hkdr(j,i) = z(9,j,i)
mnat(j,i) = z(10,j,i)
hnat(j,i) = z(11,j,i)
mnap(j,i) = z(12,j,i)
hnap(j,i) = z(13,j,i)
!gates
alpha(j,i) = z(14,j,i)
end do
end do

! Begin - semi preamble of 2nd part
pi=4.0*atan(1.0)
kh = 200.0 !!hydraulic conductivity in (\mu m)^2 mmHg^-1 s^-1
alpha0=0.2
delP = 0.4 !!pressure difference in mmHg
L=250.0 !!length between arteriole and venule in \mu m
r=15.0/L !!Vasculature radius
lp=5.0/L !! PVS length
T=(L**2)/(delP*kh) !!time in seconds
Vscale=((delP*kh)/L)*60 !! velocity in \mu m/min
D=2.3e-10 !!diffusivity of A \beta in m^2 s^-1
tort=2.04 !!tortuosity
Diffp=delP*kh*1.0e-12
Pe=(Diffp*(tort**2))/D
Pe2 = Pe*1.26 !for oligomers Pe = Pex2^(1/3)
Pe3 = Pe*1.5876

Ng=Nx_fin*Ny_fin
x_fin(1)=0.0
y_fin(1)=0.0
!grid spacing
dx=0.01
dy=0.01
dt2=0.1!dt2=0.005 ! Controls the total time (changed from 0.1 to 0.005)

do i=2,Nx_fin
x_fin(i)=x_fin(i-1)+dx
end do

do i=2,Ny_fin
y_fin(i)=y_fin(i-1)+dy
end do

do i=2,Nt_fin
tw(i)=tw(i-1)+dt2
end do

!!!!!!!!!!!!!!!!!!!!mesh_grid
  do j=1,Ny_fin
    do i=1,Nx_fin
      xm(j,i)=x_fin(i)
    end do
  end do
!!!!!!!!!!!!!!!!!
  do j=1,Ny_fin
    do i=1,Nx_fin
      ym(j,i)=y_fin(j)
    end do
  end do
!!!!!alpha for delta waves

freq = 1.0 !1 Hz
lambda = 1000.0 !1000 \mu m or 1 mm
Am = 0.1

omega = 2*pi*freq*T
kwave = (2*pi/lambda)*L
! End - semi preamble of 2nd part

call cpu_time(start)

! Time loop starts
do n=1,last_step
call SystemDerives(z,dzdt)

call rk4(z,dzdt,dt,z,SystemDerives)

do j=1,Ny
do i=1,Nx
Na_e(j,i) = z(1,j,i)
K_e(j,i) = z(2,j,i) !in mM
Na_i(j,i) = z(3,j,i)
K_i(j,i) = z(4,j,i)
Em(j,i) = z(5,j,i) !in mV
!gates
alpha(j,i) = z(14,j,i)
end do
end do
!write out
if (mod(n,(file_step)) == 0) then

print *,'time index for evolution is', n
write(numchr,*) (n/(file_step))
open (16,file="K_e"//trim(adjustl(numchr))//".txt")
write (16,*) K_e
close (16)

!write(numchr,*) (n/file_step)
!open (16,file="K_i"//trim(adjustl(numchr))//".txt")
!write (16,*) K_i
!close (16)

!write(numchr,*) (n/file_step)
!open (16,file="Na_e"//trim(adjustl(numchr))//".txt")
!write (16,*) Na_e
!close (16)

!write(numchr,*) (n/file_step)
!open (16,file="Na_i"//trim(adjustl(numchr))//".txt")
!write (16,*) Na_i
!close (16)


!write(numchr,*) (n/file_step)
open (16,file="Em"//trim(adjustl(numchr))//".txt")
write (16,*) Em
close (16)

!write(numchr,*) (n/file_step)
!open (16,file="dalphadt"//trim(adjustl(numchr))//".txt")
!write (16,*) dalphadt
!close (16)

!write(numchr,*) (n/file_step)
!open (16,file="alpha"//trim(adjustl(numchr))//".txt")
!write (16,*) alpha
!close (16)

alpha_fin_old = alpha_fin

call interp2_manual(alpha, dalphadt, alpha_fin, dalphadt_fin)

!write(numchr,*) (n/file_step)
open (16,file="alpha_fin"//trim(adjustl(numchr))//".txt")
write (16,*) alpha_fin
close (16)

call make_dalphadtfin(alpha_fin_old, alpha_fin, dalphadt_fin)

!write(numchr,*) (n/file_step)
open (16,file="dalphadt_fin"//trim(adjustl(numchr))//".txt")
write (16,*) dalphadt_fin
close (16)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! insert the 2nd loop after here
do j=1,Ny_fin
	do i=1,Nx_fin
		dalphadt_fin(j,i)=dalphadt_fin(j,i)*T 
	end do 
end do ! this was written to adjust dalphadt after reading it (idk why)
!!!!!!!!!!!!!!!!!!!!!!
PP=0.0
Pnew=0.0
u=0.0
v=0.0
md=0.0
weight=1.0
diff1=0.0
beta=0.1
!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!sources and sinks - sets boundary condition on pressure at PVS

do j=1,Ny_fin
     do i=1,Nx_fin
     	PP(j,i)=0.0_DP
     	perm(j,i) = (alpha_fin(j,i)**3/(1.0_DP - alpha_fin(j,i))**2) / &
     	        (alpha0**3/(1.0_DP - alpha0)**2)
     	S1(j,i)=beta*dalphadt_fin(j,i)
	     if ((xm(j,i)**2 + (ym(j,i)-0.5)**2) < (r+lp)**2) then
	       PP(j,i)=1.0_DP
	     end if
     end do 
end do

 	 !!!!!!!!!!!!!!!!!!!!!!!! gausss_seidel parameters
 	  !Ap = 2.0*((1.0/dx**2) + (1.0/dy**2))
 	  !Ae = 1.0/(dx**2)
 	  !Aw = 1.0/dx**2
 	  !An = 1.0/dy**2
 	  !As = 1.0/dy**2
 	  itinit=0
 	  itmax=300000
    	
 	  do it=itinit,itmax
 	 !!call gausss_seidel (only changes from orgnl code are x_fin and y_fin)
	 call darcy_flow(r, lp, i, j, Ng, weight, x_fin, y_fin, Ap, Ae, Aw, An, As, &
	                 S1, Pnew, PP, diff1, res1, dx, dy, u, v, md, mds, mdss,       &
	                 xm, ym, utot, utots, utotss, perm)
    	
 	   if ( mod(it,100000) == 0) then
 	   print *,'residual for gausss_seidel is', res1,'at iteration', it, 'and mean mass divergence', mdss
 	   end if
 	 !!end iteration loop
 	  end do
  	 !!write flow
      !do j=1,Ny_fin
      !do i=1,Nx_fin
      !   u(j,i) = 80.0*u(j,i) ! Remove the 80s here !!!
      !   v(j,i) = 80.0*v(j,i)
      !end do
      !end do
 print *,'n=', n, ' and file_step=', file_step 
  close(15) 
 ! write(numchr,*) (n/(file_step/100))
  open(15,file="P"//trim(adjustl(numchr))//".txt")
  write(15,*) PP
  close(15)

  open(15,file="u"//trim(adjustl(numchr))//".txt")
  write(15,*) u
  close(15)

  open(15,file="v"//trim(adjustl(numchr))//".txt")
  write(15,*) v
  close(15)

  !open(15,file="source"//trim(adjustl(numchr))//".txt")
  !write(15,*) S1
  !close(15)

end if


! end time loop
end do

! Start the second code
!file_step = 1
!last_step = 900000
! Semi-preamble of 2nd part

! start of the second loop (in the second code)
!do n=1,Nt_fin

	! Read and nondimensionalize the dalphadt
	!write(numchr,*) (n)
	!open(14,file="dalphadt_fin"//trim(adjustl(numchr))//".txt")
	!read(14,*) Pread
	!close(14)

	!dalphadt_fin = reshape(Pread,(/Ny_fin, Nx_fin/))

	!do j=1,Ny_fin
		!do i=1,Nx_fin
			!dalphadt_fin(j,i)=dalphadt_fin(j,i)*T 
		!end do 
	!end do
	
	! Read the alpha
     !write(numchr,*) (n)
	   !open(14,file="alpha_fin"//trim(adjustl(numchr))//".txt")
     !read(14,*) Pread2
     !close(14)

     !alpha_fin = reshape(Pread2,(/Ny_fin, Nx_fin/))    ! we dont need to read alpha and dalphadt if we join the loops
	
       !!!!!!!!!!!!!!!!!initialize
 	  ! remaining of the 2nd loop after reading alpha and dalphadt

!end do

!!! End sequence
call cpu_time(finish)

print *, "The real CPU time in seconds is ", finish-start

end program spread_depolarize ! End of main

subroutine SystemDerives(z,dzdt)

  implicit none

  integer, parameter :: DP = SELECTED_REAL_KIND(14)
  integer, parameter :: Nx = 26
  integer, parameter :: Ny = 26
  integer, parameter :: Nt = 600001

  real(KIND=DP) :: dx,dt,R,Farad,Vi,f,S,temp,phi,Imax,D_K,D_Na,eps,Cm,Vtot
  real(KIND=DP) :: gK_ghk_dr,gK_ghk_a,gk_hh,gNa_ghk_T,gNa_ghk_P,gNa_hh,ghh
  integer :: it,i,j,k,p,q,file_step,last_step
  real(KIND=DP), dimension(Nx) :: x
  real(KIND=DP), dimension(Ny) :: y
  real(KIND=DP), dimension(Ny,Nx) :: K_e,K_i,Na_e,Na_i,Em,K_tot,E_K,E_Na
  real(KIND=DP), dimension(Ny,Nx) :: I_K,I_leak,I_pump,I_K_dr,I_K_a,I_K_leak,I_K_pump
  real(KIND=DP), dimension(Ny,Nx) :: I_Na_T,I_Na_P,I_Na_leak,I_Na_pump,I_Na,I_tot
  real(KIND=DP), dimension(Ny,Nx) :: mkdr,hkdr,alphamkdr,betamkdr,alphahkdr,betahkdr
  real(KIND=DP), dimension(Ny,Nx) :: mka,hka,alphamka,betamka,alphahka,betahka
  real(KIND=DP), dimension(Ny,Nx) :: mnat,hnat,alphamnat,betamnat,alphahnat,betahnat
  real(KIND=DP), dimension(Ny,Nx) :: mnap,hnap,alphamnap,betamnap,alphahnap,betahnap
  real(KIND=DP), dimension(Ny,Nx) :: dmkadt,dhkadt,dmkdrdt,dhkdrdt,dmnatdt,dhnatdt,dmnapdt,dhnapdt
   real(KIND=DP), dimension(Ny,Nx) :: alpha, dalphadt,Na_ir,Na_er,K_ir,K_er,pi_e,pi_i,alphastar
  real(KIND=DP), dimension(Nt) :: t
  real(KIND=DP), dimension(14,Ny,Nx):: z, dzdt


  dx=0.001
  dt=5e-3
  x(1)=0.0
  y(1)=0.0
  t(1)=0.0

  do i=2,Nx
  x(i)=x(i-1)+dx
  end do

  do i=2,Ny
  y(i)=y(i-1)+dx
  end do

  do i=2,Nt
  t(i)=t(i-1)+dt
  end do

!potassium goldman-hodgkins-katz conductance
gK_ghk_dr = 100e-5 !in S/cm^2
gK_ghk_a = 10e-5
!sodium goldman-hodgkins-katz conductance
gNa_ghk_T = 100e-5
gNa_ghk_P = 2e-5
!sodium hodgkin-huxley condectance
gNa_hh = 2e-5
!sodium hodgkin-huxley condectance
gk_hh = 7e-5
!general leak conductance
ghh = 20e-5

  Farad = 96.485 !Coloumb/mM
  R = 8.51 !universal gas constant mv Coloumb/mM K

  Vi = 2160e-12 !cm^3
  Vtot = 1.15*Vi !cm^3
  S = 1586e-8 !cm^2
  f = 0.15
  temp = 298 ! K
  phi = (R*temp)/Farad !mV
  !diffusion coeffs
  D_K = 1.96e-5 !cm^2/s
  D_Na = 1.33e-5 !cm^2/s

  Imax = 0.013 !in mA/(cm^2)

  Cm = 7.5e-4 !membrane capacitance in ms/\Omega cm^2

  eps = 0.5e-3 !in cm

  do j=1,Ny
  do i=1,Nx
  Na_e(j,i) = z(1,j,i)
  K_e(j,i) = z(2,j,i) !in mM
  Na_i(j,i) = z(3,j,i)
  K_i(j,i) = z(4,j,i)
  Em(j,i) = z(5,j,i) !in mV
  !gates
  mka(j,i) = z(6,j,i)
  hka(j,i) = z(7,j,i)
  mkdr(j,i) = z(8,j,i)
  hkdr(j,i) = z(9,j,i)
  mnat(j,i) = z(10,j,i)
  hnat(j,i) = z(11,j,i)
  mnap(j,i) = z(12,j,i)
  hnap(j,i) = z(13,j,i)
  !gates
  alpha(j,i) = z(14,j,i)
  end do
  end do


  !!!!!!!!!!!!!!!!!!!
  do j=1,Ny
  do i=1,Nx
  pi_e(j,i) = Na_e(j,i) + K_e(j,i) + 3.86+141.6
  pi_i(j,i) = Na_i(j,i) + K_i(j,i) + 9.82+133.45
  end do
  end do
!!!!!!!!!dalphadt
  do j=1,Ny
  do i=1,Nx
  alphastar(j,i) = 1.0 - 0.8*(1.21 - 0.21*exp((pi_e(j,i)-pi_i(j,i))/20)) !compare with Eq. 10 in Wei et al., but with alpha0=0.2
  end do
  end do

  !dalphadt(i) = -(Pf*Vw*S/Vtot)*(A(i)-(0.8/(1-alpha(i)))*Ar(i))
     !dalphadt(i)=0
  do j=1,Ny   
  do i=1,Nx   
     dalphadt(j,i)=(alphastar(j,i) - alpha(j,i))/250   !Wei,Ullah_et_al. J. NeuroSci
  end do
  end do


  !pump currents
  do j=1,Ny
  do i=1,Nx
  I_pump(j,i) = Imax/(((1 + 2.0/K_e(j,i))**2)*((1 + 7.7/Na_i(j,i))**3))
  end do
  end do

  do j=1,Ny
  do i=1,Nx
  I_K_pump(j,i) = -2.0*I_pump(j,i)
  I_Na_pump(j,i) = 3.0*I_pump(j,i)
  end do
  end do
 !potassium transient current
  call ka(p,q,alphamka,betamka,alphahka,betahka,x,Em)
  call gates(dmkadt,mka,alphamka,betamka)
  call gates(dhkadt,hka,alphahka,betahka)
  call GHK(p,q,Farad,phi,gK_ghk_a,I_K_a,mka,hka,K_e,K_i,Em,x)

  !delayed rectified potassium
  call kdr(p,q,alphamkdr,betamkdr,alphahkdr,betahkdr,x,Em)
  call gates(dmkdrdt,mkdr,alphamkdr,betamkdr)
  call gates(dhkdrdt,hkdr,alphahkdr,betahkdr)
  call GHK(p,q,Farad,phi,gK_ghk_dr,I_K_dr,mkdr,hkdr,K_e,K_i,Em,x)
  !potassium leak current
  do j=1,Ny
  do i=1,Nx
  E_K(j,i) = phi*log(K_e(j,i)/K_i(j,i))
  end do
  end do
  call HH(p,q,gk_hh,phi,K_e,K_i,Em,I_K_leak,E_K,x)

  !total potassium current
  do j=1,Ny
  do i=1,Nx
  I_K(j,i)=I_K_dr(j,i)+I_K_a(j,i)+I_K_leak(j,i)+I_K_pump(j,i)
  end do
  end do


  !fast transient sodium current
  call nat(p,q,alphamnat,betamnat,alphahnat,betahnat,x,Em)
  call gates(dmnatdt,mnat,alphamnat,betamnat)
  call gates(dhnatdt,hnat,alphahnat,betahnat)
  call GHK(p,q,Farad,phi,gNa_ghk_T,I_Na_T,mnat,hnat,Na_e,Na_i,Em,x)

  !persistant sodium current
  call nap(p,q,alphamnap,betamnap,alphahnap,betahnap,x,Em)
  call gates(dmnapdt,mnap,alphamnap,betamnap)
  call gates(dhnapdt,hnap,alphahnap,betahnap)
  call GHK(p,q,Farad,phi,gNa_ghk_P,I_Na_P,mnap,hnap,Na_e,Na_i,Em,x)

  !sodium leak current
  do j=1,Ny
  do i=1,Nx
  E_Na(j,i) = phi*log(Na_e(j,i)/Na_i(j,i))
  end do
  end do

  call HH(p,q,gNa_hh,phi,Na_e,Na_i,Em,I_Na_leak,E_Na,x)

  !total sodium current
  do j=1,Ny
  do i=1,Nx
  I_Na(j,i) = I_Na_T(j,i)+I_Na_P(j,i)+I_Na_leak(j,i)+I_Na_pump(j,i)
  end do
  end do

  !general leak current
  do j=1,Ny
  do i=1,Nx
  I_leak(j,i) = ghh*(Em(j,i)+70.0)
  end do
  end do

  !total current across membrane
  do j=1,Ny
  do i=1,Nx
  I_tot(j,i) = I_Na(j,i) + I_K(j,i) + I_leak(j,i)
  end do
  end do

  do j=2,Ny-1
  do i=2,Nx-1
  dzdt(1,j,i) = ((S/(Farad*Vtot*alpha(j,i)))*I_Na(j,i))*1e-3 + &
 (D_Na*(Na_e(j,i+1)-2*Na_e(j,i)+Na_e(j,i-1))/(dx**2))*1e-3 +&
 (D_Na*(Na_e(j+1,i)-2*Na_e(j,i)+Na_e(j-1,i))/(dx**2))*1e-3+&

 ((D_Na/alpha(j,i))*(Na_e(j,i+1)-Na_e(j,i-1))*(alpha(j,i+1)-alpha(j,i-1))/(4*dx**2))*1e-3+&

 ((D_Na/alpha(j,i))*(Na_e(j+1,i)-Na_e(j-1,i))*(alpha(j+1,i)-alpha(j-1,i))/(4*dx**2))*1e-3-&

  (Na_e(j,i)/alpha(j,i))*dalphadt(j,i)

  dzdt(2,j,i) = ((S/(Farad*Vtot*alpha(j,i)))*I_K(j,i))*1e-3 + &
 (D_K*(K_e(j,i+1)-2*K_e(j,i)+K_e(j,i-1))/(dx**2))*1e-3 +&
 (D_K*(K_e(j+1,i)-2*K_e(j,i)+K_e(j-1,i))/(dx**2))*1e-3 +&

 ((D_K/alpha(j,i))*(K_e(j,i+1)-K_e(j,i-1))*(alpha(j,i+1)-alpha(j,i-1))/(4*dx**2))*1e-3+&

 ((D_K/alpha(j,i))*(K_e(j+1,i)-K_e(j-1,i))*(alpha(j+1,i)-alpha(j-1,i))/(4*dx**2))*1e-3-&

(K_e(j,i)/alpha(j,i))*dalphadt(j,i)
  end do
  end do

!boundary conditions of Na_e, K_e
!xwalls
do j=2,Ny-1

  dzdt(1,j,1) = ((S/(Farad*Vtot*alpha(j,1)))*I_Na(j,1))*1e-3 + &
 (D_Na*(2*Na_e(j,2)-2*Na_e(j,1))/(dx**2))*1e-3 +&
 (D_Na*(Na_e(j+1,1)-2*Na_e(j,1)+Na_e(j-1,1))/(dx**2))*1e-3+&
 ((D_Na/alpha(j,1))*(Na_e(j+1,1)-Na_e(j-1,1))*(alpha(j+1,1)-alpha(j-1,1))/(4*dx**2))*1e-3-&
 (Na_e(j,1)/alpha(j,1))*dalphadt(j,1) !Na at x=0

 dzdt(1,j,Nx) = ((S/(Farad*Vtot*alpha(j,Nx)))*I_Na(j,Nx))*1e-3 + &
 (D_Na*(2*Na_e(j,Nx-1)-2*Na_e(j,Nx))/(dx**2))*1e-3 +&
 (D_Na*(Na_e(j+1,Nx)-2*Na_e(j,Nx)+Na_e(j-1,Nx))/(dx**2))*1e-3+&
 ((D_Na/alpha(j,Nx))*(Na_e(j+1,Nx)-Na_e(j-1,Nx))*(alpha(j+1,Nx)-alpha(j-1,Nx))/(4*dx**2))*1e-3-&
 (Na_e(j,Nx)/alpha(j,Nx))*dalphadt(j,Nx) !Na at x=Nx

 dzdt(2,j,1) = ((S/(Farad*Vtot*alpha(j,1)))*I_K(j,1))*1e-3 + &
 (D_K*(2*K_e(j,2)-2*K_e(j,1))/(dx**2))*1e-3 +&
 (D_K*(K_e(j+1,1)-2*K_e(j,1)+K_e(j-1,1))/(dx**2))*1e-3 +&
 ((D_K/alpha(j,1))*(K_e(j+1,1)-K_e(j-1,1))*(alpha(j+1,1)-alpha(j-1,1))/(4*dx**2))*1e-3-&
 (K_e(j,1)/alpha(j,1))*dalphadt(j,1) !K at x=0

 dzdt(2,j,Nx) = ((S/(Farad*Vtot*alpha(j,Nx)))*I_K(j,Nx))*1e-3 + &
 (D_K*(2*K_e(j,Nx-1)-2*K_e(j,Nx))/(dx**2))*1e-3 +&
 (D_K*(K_e(j+1,Nx)-2*K_e(j,Nx)+K_e(j-1,Nx))/(dx**2))*1e-3 +&
 ((D_K/alpha(j,Nx))*(K_e(j+1,Nx)-K_e(j-1,Nx))*(alpha(j+1,Nx)-alpha(j-1,Nx))/(4*dx**2))*1e-3-&
 (K_e(j,Nx)/alpha(j,Nx))*dalphadt(j,Nx) !K at x=Nx

end do

!ywalls
do i=2,Nx-1
  
  dzdt(1,1,i) = ((S/(Farad*Vtot*alpha(1,i)))*I_Na(1,i))*1e-3 + &
 (D_Na*(Na_e(1,i+1)-2*Na_e(1,i)+Na_e(1,i-1))/(dx**2))*1e-3 +&
 (D_Na*(2*Na_e(2,i)-2*Na_e(1,i))/(dx**2))*1e-3+&
  ((D_Na/alpha(1,i))*(Na_e(1,i+1)-Na_e(1,i-1))*(alpha(1,i+1)-alpha(1,i-1))/(4*dx**2))*1e-3-&
  (Na_e(1,i)/alpha(1,i))*dalphadt(1,i) !Na at y=0

  dzdt(1,Nx,i) = ((S/(Farad*Vtot*alpha(Nx,i)))*I_Na(Nx,i))*1e-3 + &
 (D_Na*(Na_e(Nx,i+1)-2*Na_e(Nx,i)+Na_e(Nx,i-1))/(dx**2))*1e-3 +&
 (D_Na*(2*Na_e(Nx-1,i)-2*Na_e(Nx,i))/(dx**2))*1e-3+&
  ((D_Na/alpha(Nx,i))*(Na_e(Nx,i+1)-Na_e(Nx,i-1))*(alpha(Nx,i+1)-alpha(Nx,i-1))/(4*dx**2))*1e-3-&
  (Na_e(Nx,i)/alpha(Nx,i))*dalphadt(Nx,i) !Na at y=Nx


  dzdt(2,1,i) = ((S/(Farad*Vtot*alpha(1,i)))*I_K(1,i))*1e-3 + &
 (D_K*(K_e(1,i+1)-2*K_e(1,i)+K_e(1,i-1))/(dx**2))*1e-3 +&
 (D_K*(2*K_e(2,i)-2*K_e(1,i))/(dx**2))*1e-3 +&
 ((D_K/alpha(1,i))*(K_e(1,i+1)-K_e(1,i-1))*(alpha(1,i+1)-alpha(1,i-1))/(4*dx**2))*1e-3-&
 (K_e(1,i)/alpha(1,i))*dalphadt(1,i) !K at y=0

  dzdt(2,Ny,i) = ((S/(Farad*Vtot*alpha(Ny,i)))*I_K(Ny,i))*1e-3 + &
 (D_K*(K_e(Ny,i+1)-2*K_e(Ny,i)+K_e(Ny,i-1))/(dx**2))*1e-3 +&
 (D_K*(2*K_e(Ny-1,i)-2*K_e(Ny,i))/(dx**2))*1e-3 +&
 ((D_K/alpha(Ny,i))*(K_e(Ny,i+1)-K_e(Ny,i-1))*(alpha(Ny,i+1)-alpha(Ny,i-1))/(4*dx**2))*1e-3-&
 (K_e(Ny,i)/alpha(Ny,i))*dalphadt(Ny,i) !K at y=Ny

end do

!!Na corners

dzdt(1,1,1) = ((S/(Farad*Vtot*alpha(1,1)))*I_Na(1,1))*1e-3 + &
 (D_Na*(2*Na_e(1,2)-2*Na_e(1,1))/(dx**2))*1e-3 +&
 (D_Na*(2*Na_e(2,1)-2*Na_e(1,1))/(dx**2))*1e-3-&
  (Na_e(1,1)/alpha(1,1))*dalphadt(1,1)

  dzdt(1,Ny,Nx) = ((S/(Farad*Vtot*alpha(Ny,Nx)))*I_Na(Ny,Nx))*1e-3 + &
 (D_Na*(2*Na_e(Ny,Nx-1)-2*Na_e(Ny,Nx))/(dx**2))*1e-3 +&
 (D_Na*(2*Na_e(Ny-1,Nx)-2*Na_e(Ny,Nx))/(dx**2))*1e-3-&
  (Na_e(Ny,Nx)/alpha(Ny,Nx))*dalphadt(Ny,Nx)

  dzdt(1,Ny,1) = ((S/(Farad*Vtot*alpha(Ny,1)))*I_Na(Ny,1))*1e-3 + &
 (D_Na*(2*Na_e(Ny,2)-2*Na_e(Ny,1))/(dx**2))*1e-3 +&
 (D_Na*(2*Na_e(Ny-1,1)-2*Na_e(Ny,1))/(dx**2))*1e-3-&
  (Na_e(Ny,1)/alpha(Ny,1))*dalphadt(Ny,1)

  dzdt(1,1,Nx) = ((S/(Farad*Vtot*alpha(1,Nx)))*I_Na(1,Nx))*1e-3 + &
 (D_Na*(2*Na_e(1,Nx-1)-2*Na_e(1,Nx))/(dx**2))*1e-3 +&
 (D_Na*(2*Na_e(2,Nx)-2*Na_e(1,Nx))/(dx**2))*1e-3-&
  (Na_e(1,Nx)/alpha(1,Nx))*dalphadt(1,Nx)


!! K_e corners
  
  dzdt(2,1,1) = ((S/(Farad*Vtot*alpha(1,1)))*I_K(1,1))*1e-3 + &
 (D_K*(2*K_e(1,2)-2*K_e(1,1))/(dx**2))*1e-3 +&
 (D_K*(2*K_e(2,1)-2*K_e(1,1))/(dx**2))*1e-3 -&
 (K_e(1,1)/alpha(1,1))*dalphadt(1,1)

 dzdt(2,Ny,Nx) = ((S/(Farad*Vtot*alpha(Ny,Nx)))*I_K(Ny,Nx))*1e-3 + &
 (D_K*(2*K_e(Ny,Nx-1)-2*K_e(Ny,Nx))/(dx**2))*1e-3 +&
 (D_K*(2*K_e(Ny-1,Nx)-2*K_e(Ny,Nx))/(dx**2))*1e-3 -&
 (K_e(Ny,Nx)/alpha(Ny,Nx))*dalphadt(Ny,Nx)

 dzdt(2,Ny,1) = ((S/(Farad*Vtot*alpha(Ny,1)))*I_K(Ny,1))*1e-3 + &
 (D_K*(2*K_e(Ny,2)-2*K_e(Ny,1))/(dx**2))*1e-3 +&
 (D_K*(2*K_e(Ny-1,1)-2*K_e(Ny,1))/(dx**2))*1e-3 -&
 (K_e(Ny,1)/alpha(Ny,1))*dalphadt(Ny,1)

 dzdt(2,1,Nx) = ((S/(Farad*Vtot*alpha(1,Nx)))*I_K(1,Nx))*1e-3 + &
 (D_K*(2*K_e(1,Nx-1)-2*K_e(1,Nx))/(dx**2))*1e-3 +&
 (D_K*(2*K_e(2,Nx)-2*K_e(1,Nx))/(dx**2))*1e-3 -&
 (K_e(1,Nx)/alpha(1,Nx))*dalphadt(1,Nx)


!!

!!update with alpha
  do j=1,Ny
  do i=1,Nx
  dzdt(3,j,i) = - ((S/(Farad*Vtot*(1-alpha(j,i))))*I_Na(j,i))*1e-3+&
  (Na_i(j,i)/(1-alpha(j,i)))*dalphadt(j,i)
  dzdt(4,j,i) = - ((S/(Farad*Vtot*(1-alpha(j,i))))*I_K(j,i))*1e-3+&
  (K_i(j,i)/(1-alpha(j,i)))*dalphadt(j,i)
  dzdt(5,j,i) = - (I_tot(j,i))/Cm

  dzdt(6,j,i) = dmkadt(j,i)
  dzdt(7,j,i) = dhkadt(j,i)
  dzdt(8,j,i) = dmkdrdt(j,i)
  dzdt(9,j,i) = dhkdrdt(j,i)
  dzdt(10,j,i) = dmnatdt(j,i)
  dzdt(11,j,i) = dhnatdt(j,i)
  dzdt(12,j,i) = dmnapdt(j,i)
  dzdt(13,j,i) = dhnapdt(j,i)
  !alpha
  dzdt(14,j,i) = dalphadt(j,i)
  end do
  end do

end subroutine SystemDerives

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!subroutines
!!!subroutine kdr
subroutine kdr(p,q,alpham,betam,alphah,betah,x,Em)

implicit none

integer, parameter :: DP = SELECTED_REAL_KIND(14)
integer, parameter :: Nx = 26
integer, parameter :: Ny = 26

integer :: i,j,k,p,q
real(KIND=DP), dimension(Ny,Nx) :: alpham,alphah,betam,betah,Em
real(KIND=DP), dimension(Nx) :: x

p = 2
q = 0

do j=1,Ny
do i=1,Nx
alpham(j,i) = 0.016*(Em(j,i)+34.9)/(1-exp(-(0.2*Em(j,i) + 6.98)))
betam(j,i) = 0.25*exp(-(0.025*Em(j,i) + 1.25))

alphah(j,i) = 0.0
betah(j,i) = 0.0
end do
end do

end subroutine kdr
!!!subroutine ka
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ka(p,q,alpham,betam,alphah,betah,x,Em)

implicit none

integer, parameter :: DP = SELECTED_REAL_KIND(14)
integer, parameter :: Nx = 26
integer, parameter :: Ny = 26

integer :: i,j,k,p,q
real(KIND=DP), dimension(Ny,Nx) :: alpham,alphah,betam,betah,Em
real(KIND=DP), dimension(Nx) :: x

p = 2
q = 1

do j=1,Ny
do i=1,Nx
  alpham(j,i) = 0.02*(Em(j,i)+56.9)/(1-exp(-(0.1*Em(j,i)+5.69)))
  betam(j,i) = 0.0175*(Em(j,i)+29.9)/(exp(0.1*Em(j,i)+2.99)-1)

  alphah(j,i) = 0.016*exp(-(0.056*Em(j,i)+4.61))
  betah(j,i) = 0.5/(1+exp(-(0.2*Em(j,i)+11.98)))
end do
end do

end subroutine ka
!!!subroutine nap
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine nap(p,q,alpham,betam,alphah,betah,x,Em)

implicit none

integer, parameter :: DP = SELECTED_REAL_KIND(14)
integer, parameter :: Nx = 26
integer, parameter :: Ny = 26

integer :: i,j,k,p,q
real(KIND=DP), dimension(Ny,Nx) :: alpham,alphah,betam,betah,Em
real(KIND=DP), dimension(Nx) :: x

p = 2
q = 1

do j=1,Ny
do i=1,Nx
alpham(j,i) = 1/(6*(1+exp(-(0.143*Em(j,i)+5.67))))
betam(j,i) = exp(-(0.143*Em(j,i)+5.67))*alpham(j,i)

alphah(j,i) = (5.12e-8)*exp(-(0.056*Em(j,i)+2.94))
betah(j,i) = (1.6e-6)/(1+exp(-(0.2*Em(j,i)+8)))
end do
end do

end subroutine nap
!!!subroutine nat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine nat(p,q,alpham,betam,alphah,betah,x,Em)

implicit none

integer, parameter :: DP = SELECTED_REAL_KIND(14)
integer, parameter :: Nx = 26
integer, parameter :: Ny = 26

integer :: i,j,k,p,q
real(KIND=DP), dimension(Ny,Nx) :: alpham,alphah,betam,betah,Em
real(KIND=DP), dimension(Nx) :: x

p = 3
q = 1

do j=1,Ny
do i=1,Nx
  alpham(j,i) = 0.32*(Em(j,i)+51.9)/(1-exp(-(0.25*Em(j,i)+12.975)))
  betam(j,i) = 0.28*(Em(j,i)+24.89)/(exp(0.2*Em(j,i)+4.978)-1)

  alphah(j,i) = 0.128*exp(-(0.056*Em(j,i)+2.94))
  betah(j,i) = 4/(1+exp(-(0.2*Em(j,i)+6)))
end do
end do

end subroutine nat

subroutine GHK(p,q,Farad,phi,gion_ghk,I_GHK,m,h,ion_e,ion_i,Em,x)

implicit none

integer, parameter :: DP = SELECTED_REAL_KIND(14)
integer, parameter :: Nx = 26
integer, parameter :: Ny = 26

integer :: i,k,p,q,j
real(KIND=DP) :: Farad,phi,gion_ghk
real(KIND=DP), dimension(Nx) :: x
real(KIND=DP), dimension(Ny,Nx) :: I_GHK,m,h,ion_e,ion_i,Em

do j=1,Ny
do i=1,Nx
I_GHK(j,i) = gion_ghk*(m(j,i)**p)*(h(j,i)**q)*(Farad*Em(j,i)*(ion_i(j,i) - &
(exp(-Em(j,i)/phi))*ion_e(j,i)))/(phi*(1 - exp(-Em(j,i)/phi)))
end do
end do

end subroutine GHK

subroutine HH(p,q,gion_hh,phi,ion_e,ion_i,Em,I_HH,E_ion,x)

  implicit none

  integer, parameter :: DP = SELECTED_REAL_KIND(14)
  integer, parameter :: Nx = 26
  integer, parameter :: Ny = 26

  integer :: i,j,k,p,q
  real(KIND=DP) :: gion_hh,phi
  real(KIND=DP), dimension(Ny,Nx) :: ion_e,ion_i,Em,I_HH, E_ion
  real(KIND=DP), dimension(Nx) :: x

  do j=1,Ny
  do i=1,Nx
  I_HH(j,i) = gion_hh*(Em(j,i) - E_ion(j,i))
  end do
 end do

end subroutine HH

subroutine gates(dmdt,m,alpham,betam)

  implicit none

  integer, parameter :: DP = SELECTED_REAL_KIND(14)
  integer, parameter :: Nx = 26
  integer, parameter :: Ny = 26

  integer :: i,j,k,p,q
  real(KIND=DP), dimension(Ny,Nx) :: m,alpham,betam,dmdt
  real(KIND=DP), dimension(Nx) :: x

  do j=1,Ny
  do i=1,Nx
    dmdt(j,i) = (alpham(j,i)*(1.0-m(j,i))-betam(j,i)*m(j,i))
  end do
  end do

end subroutine gates

!!!!!!!!!!!!!!!!!!!!!rk4
subroutine rk4(z,dzdt,h,zout,derivs)
  integer, parameter :: DP = SELECTED_REAL_KIND(14)
  integer, parameter :: Nx = 26
  integer, parameter :: Ny = 26
  EXTERNAL derivs
  INTEGER i,j,k
  real(KIND=DP), dimension(14,Ny,Nx):: z, dzdt, dzt, dzm, zt, zout
  real(KIND=DP) :: h,hh,h6


  hh=h*0.5d0
  h6=h/6.0d0

  do i = 1,14     ! First step.
  do k=1,Ny
  do j=1,Nx
     zt(i,k,j)=z(i,k,j)+hh*dzdt(i,k,j)
  end do
  end do
  end do
  call derivs(zt,dzt)    !  Secondstep.
  do i = 1,14
  do k=1,Ny
  do j=1,Nx
     zt(i,k,j)=z(i,k,j)+hh*dzt(i,k,j)
  end do
  end do
  end do

  call derivs(zt,dzm)    !  Thirdstep.
  do i = 1,14
  do k=1,Ny
  do j=1,Nx
     zt(i,k,j)=z(i,k,j)+h*dzm(i,k,j)
     dzm(i,k,j)=dzt(i,k,j)+dzm(i,k,j)
  end do
  end do
  end do

  call derivs(zt,dzt)   !  Fourthstep.
  do i = 1,14
  do k=1,Ny
  do j=1,Nx        !  Accumulate increments with proper weights.
     zout(i,k,j)=z(i,k,j)+h6*(dzdt(i,k,j)+dzt(i,k,j)+2.0d0*dzm(i,k,j))
  end do
  end do
  end do
!return

end subroutine rk4

subroutine interp2_manual(alpha, dalphadt, alpha_fin, dalphadt_fin) ! Hard coded for 26 interpolated to 101 points (OR A+1 to 4A+1)
     integer, parameter :: DP = SELECTED_REAL_KIND(14)
     integer, parameter :: Nx = 26
     integer, parameter :: Ny = 26
     integer, parameter :: Nx_fin = 101
     integer, parameter :: Ny_fin = 101
	integer :: i,j,k,k1,k2
	real(KIND=DP), dimension(Ny,Nx) :: alpha, dalphadt
	real(KIND=DP), dimension(Ny_fin,Nx_fin) :: alpha_fin, dalphadt_fin
	real(KIND=DP) :: s,t
	
	do j=1,Ny-1
		do i=1,Nx-1
			do k1=3,0,-1
				s=real(3-k1)/4.0
				do k2=3,0,-1
					t=real(3-k2)/4.0
					
alpha_fin(4*j-k2,4*i-k1) = (1.0-s)*(1.0-t)*alpha(j,i) &
											+ s*(1.0-t)*alpha(j,i+1) &
											+ (1.0-s)*t*alpha(j+1,i) &
											+ s*t*alpha(j+1,i+1)
!dalphadt_fin(4*j-k2,4*i-k1) = (1-s)*(1-t)*dalphadt(j,i) &
!		 									+ s*(1-t)*dalphadt(j,i+1) &
!											+ (1-s)*t*dalphadt(j+1,i) &
!											+ s*t*dalphadt(j+1,i+1)

				end do
			end do
		end do
	end do
	! the last row and column
	alpha_fin(Ny_fin, Nx_fin)    = alpha(Ny, Nx)
	dalphadt_fin(Ny_fin, Nx_fin) = dalphadt(Ny, Nx)
	do i = 1, Nx-1
	   do k1 = 3, 0, -1
	      s = (3.0d0 - dble(k1)) / 4.0d0     ! s = 0,0.25,0.5,0.75
	      alpha_fin(101, 4*i - k1) = (1.0d0 - s)*alpha(Ny, i) + s*alpha(Ny, i+1)
	      dalphadt_fin(101, 4*i - k1) = (1.0d0 - s)*dalphadt(Ny, i) + s*dalphadt(Ny, i+1)
	   end do
	end do
	alpha_fin(Ny_fin, Nx_fin)    = alpha(Ny, Nx)
	dalphadt_fin(Ny_fin, Nx_fin) = dalphadt(Ny, Nx)
	
end subroutine interp2_manual

subroutine make_dalphadtfin(alpha_fin_old, alpha_fin, dalphadt_fin)
    integer, parameter :: DP = SELECTED_REAL_KIND(14)
    integer, parameter :: Nx_fin = 101
    integer, parameter :: Ny_fin = 101
    integer, parameter :: Nt_fin = 1000
    integer :: i,j
    real(KIND=DP), dimension(Ny_fin,Nx_fin) :: alpha_fin_old, alpha_fin, dalphadt_fin

    ! (alpha(k)-alpha(k-1))/0.1 at every x and y
    ! the 0.1 is due to the dt in 2nd part of the code
  
      do i=1, Nx_fin
        do j=1, Ny_fin
          dalphadt_fin(j,i) = (alpha_fin(j,i) - alpha_fin_old(j,i))/0.1
        end do
      end do
  

end subroutine make_dalphadtfin

!Subroutine darcy flow from the second code)
subroutine darcy_flow(r,lp,i,j,Ng,weight,x,y,Ap,Ae,Aw,An,As,S1,Pnew,P,diff1,res1,&
                dx,dy,u,v,md,mds,mdss,xm,ym,utot,utots,utotss,perm)
implicit none

integer, parameter :: DP = SELECTED_REAL_KIND(14)
integer, parameter :: Nx = 101
integer, parameter :: Ny = 101
integer, parameter :: Nt = 1000

real(KIND=DP) :: dx,dy,Ap,Ae,Aw,An,As,weight,res1,Ng,toler,r,lp,L,mdss,utotss,kE,kW,kN,kS,res_sum
integer :: i,j,k
real(KIND=DP), dimension(Nx)::x
real(KIND=DP), dimension(Ny)::y,mds,utots
real(KIND=DP), dimension(Ny,Nx)::P,Pnew,u,v,S1,xm,ym,diff1,md,utot,perm

do j=2,Ny-1
  do i=2,Nx-1
	  kE = 0.5_DP*(perm(j,i+1) + perm(j,i))
	  kW = 0.5_DP*(perm(j,i-1) + perm(j,i))
	  kN = 0.5_DP*(perm(j+1,i) + perm(j,i))
	  kS = 0.5_DP*(perm(j-1,i) + perm(j,i))

	  Ae = kE/dx**2
	  Aw = kW/dx**2
	  An = kN/dy**2
	  As = kS/dy**2
	  Ap = Ae + Aw + An + As
	  
	  Pnew(j,i) = (1.0/Ap)*((1.0-weight)*Ap*P(j,i) &
	  + weight*(As*P(j-1,i)+An*P(j+1,i)+Ae*P(j,i+1)+Aw*P(j,i-1))&
	  -S1(j,i))
  end do
end do

do j=2,Ny-1
  do i=2,Nx-1
	  kE = 0.5_DP*(perm(j,i+1) + perm(j,i))
	  kW = 0.5_DP*(perm(j,i-1) + perm(j,i))
	  kN = 0.5_DP*(perm(j+1,i) + perm(j,i))
	  kS = 0.5_DP*(perm(j-1,i) + perm(j,i))
	  diff1(j,i) = ( &
	     ( kE*(P(j,i+1)-P(j,i)) - kW*(P(j,i)-P(j,i-1)) )/dx**2 &
	   + ( kN*(P(j+1,i)-P(j,i)) - kS*(P(j,i)-P(j-1,i)) )/dy**2 &
	   - S1(j,i) )**2
	   res_sum = res_sum + diff1(j,i)
 end do
end do

!!Neumann boundary conditions - noflux
!!east-west
do j=2,Ny-1
	kE = 0.5_DP*(perm(j,2) + perm(j,1))
	kN = 0.5_DP*(perm(j+1,1) + perm(j,1))
	kS = 0.5_DP*(perm(j-1,1) + perm(j,1))

	Ae = kE/dx**2
	Aw = Ae      ! ghost/mirror no-flux treatment
	An = kN/dy**2
	As = kS/dy**2
	Ap = Ae + Aw + An + As
	
 	Pnew(j,1) = (1.0/Ap)*((Ae+Aw)*Pnew(j,2) + An*Pnew(j+1,1)&
  	+ As*Pnew(j-1,1)-S1(j,1))
  
	kW = 0.5_DP*(perm(j,Nx-1) + perm(j,Nx))
	kN = 0.5_DP*(perm(j+1,Nx) + perm(j,Nx))
	kS = 0.5_DP*(perm(j-1,Nx) + perm(j,Nx))
	
	Aw = kW/dx**2
	Ae = Aw      ! ghost/mirror no-flux treatment
	An = kN/dy**2
	As = kS/dy**2
	Ap = Ae + Aw + An + As
	
	Pnew(j,Nx) = (1/Ap)*((Ae+Aw)*Pnew(j,Nx-1)&
	+ An*Pnew(j+1,Nx) + As*Pnew(j-1,Nx)-S1(j,Nx))
end do

!!north-south
do i=2,Nx-1
	kE = 0.5_DP*(perm(1,i+1) + perm(1,i))
	kW = 0.5_DP*(perm(1,i-1) + perm(1,i))
	kN = 0.5_DP*(perm(2,i) + perm(1,i))

	Ae = kE/dx**2
	Aw = kW/dx**2
	An = kN/dy**2
	As = An      ! ghost/mirror no-flux treatment
	Ap = Ae + Aw + An + As
 	
	Pnew(1,i)=(1.0/Ap)*(Ae*Pnew(1,i+1)+Aw*Pnew(1,i-1)&
  	+ (An+As)*Pnew(2,i)-S1(1,i))
	
	kE = 0.5_DP*(perm(Ny,i+1) + perm(Ny,i))
	kW = 0.5_DP*(perm(Ny,i-1) + perm(Ny,i))
	kS = 0.5_DP*(perm(Ny-1,i) + perm(Ny,i))

	Ae = kE/dx**2
	Aw = kW/dx**2
	As = kS/dy**2
	An = As      ! ghost/mirror no-flux treatment
	Ap = Ae + Aw + An + As
	
     Pnew(Ny,i)=(1.0/Ap)*(Ae*Pnew(Ny,i+1)&
     +Aw*Pnew(Ny,i-1) + (An+As)*Pnew(Ny-1,i)-S1(Ny,i))
end do

!!corners
kE = 0.5_DP*(perm(1,2) + perm(1,1))
kN = 0.5_DP*(perm(2,1) + perm(1,1))
Ae = kE/dx**2
Aw = Ae      ! mirror/no-flux in west direction
An = kN/dy**2
As = An      ! mirror/no-flux in south direction
Ap = Ae + Aw + An + As
Pnew(1,1)=(1.0/Ap)*((Ae+Aw)*Pnew(1,2)&
 + (An+As)*Pnew(2,1)-S1(1,1))
 
 kW = 0.5_DP*(perm(Ny,Nx-1) + perm(Ny,Nx))
 kS = 0.5_DP*(perm(Ny-1,Nx) + perm(Ny,Nx))
 Aw = kW/dx**2
 Ae = Aw      ! mirror/no-flux in east direction
 As = kS/dy**2
 An = As      ! mirror/no-flux in north direction
 Ap = Ae + Aw + An + As
 Pnew(Ny,Nx)=(1.0/Ap)*((Ae+Aw)*Pnew(Ny,Nx-1)&
 + (An+As)*Pnew(Ny-1,Nx)-S1(Ny,Nx))
 
 kW = 0.5_DP*(perm(1,Nx-1) + perm(1,Nx))
 kN = 0.5_DP*(perm(2,Nx) + perm(1,Nx))
 Aw = kW/dx**2
 Ae = Aw      ! mirror/no-flux in east direction
 An = kN/dy**2
 As = An      ! mirror/no-flux in south direction
 Ap = Ae + Aw + An + As
 Pnew(1,Nx)=(1.0/Ap)*((Ae+Aw)*Pnew(1,Nx-1)&
 + (An+As)*Pnew(2,Nx)-S1(1,Nx))
 
 kE = 0.5_DP*(perm(Ny,2) + perm(Ny,1))
 kS = 0.5_DP*(perm(Ny-1,1) + perm(Ny,1))
 Ae = kE/dx**2
 Aw = Ae      ! mirror/no-flux in west direction
 As = kS/dy**2
 An = As      ! mirror/no-flux in north direction
 Ap = Ae + Aw + An + As
 Pnew(Ny,1)=(1.0/Ap)*((Ae+Aw)*Pnew(Ny,2)&
 + (An+As)*Pnew(Ny-1,1)-S1(Ny,1))

!!fixed sources and sinks

!!!!!!!!!!!!!!!!!!!!sources and sinks - sets boundary condition on pressure at PVS
!!!!!!!!!!
do j=1,Ny
	do i=1,Nx
      	!PP(j,i)=0.0
      	if ((xm(j,i)**2 + (ym(j,i)-0.5)**2) < (r+lp)**2) then
        	  Pnew(j,i)=1.0
      	else if (((xm(j,i)-1)**2 + (ym(j,i)-0.5)**2) < (r+lp)**2) then
        	  Pnew(j,i)=0.0
      	else
        	  Pnew(j,i)=Pnew(j,i)
          end if
	end do
end do
	
res1 = ((1/Ng)*(sum(diff1)**0.5))
!Pnew = Pnew - sum(Pnew)/size(Pnew)
P=Pnew

!!velocities
do j=1,Ny
    do i=2,Nx-1
        !u(j,i) = -(P(j,i+1)-P(j,i-1))/(2*dx)
	   u(j,i) = -perm(j,i)*(P(j,i+1)-P(j,i-1))/(2.0_DP*dx)
    end do
end do

do j=2,Ny-1
    do i=1,Nx
        !v(j,i) = -(P(j+1,i)-P(j-1,i))/(2*dy)
	   v(j,i) = -perm(j,i)*(P(j+1,i)-P(j-1,i))/(2.0_DP*dy)
    end do
end do

do j=2,Ny-1
    do i=2,Nx-1
     md(j,i) = (u(j,i+1)-u(j,i-1))/(2*dx) + (v(j+1,i)-v(j-1,i))/(2*dy)
   end do
end do

do j=2,Ny-1
    do i=2,Nx-1
     utot(j,i) = (u(j,i)**2 + v(j,i)**2)**(0.5)
   end do
end do

do j=1,Ny
  mds(j) = sum(md(j,:))/(size(md,2))
end do

mdss = sum(mds)/(size(mds,1))

do j=1,Ny
  utots(j) = sum(utot(j,:))/(size(utot,2))
end do

utotss = sum(utots)/(size(utots,1))
!!
end subroutine darcy_flow
!!!!!!!!!!!!!!!!!!!!!!!!!!subroutine advective transport	
	
