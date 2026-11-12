program darcy_interstitium

  integer, parameter :: DP = SELECTED_REAL_KIND(14)
  integer, parameter :: Nx = 101
  integer, parameter :: Ny = 101
  integer, parameter :: Nt = 1000

  real(KIND=DP) :: dx,dy,dt,Ap,Ae,Aw,An,As,weight,res1,res2,Ng,toler,r,lp,L,mdss
  real(KIND=DP) :: kh,delP,T,Vscale,tort,Diffp,Pe,Pe2,Le3,D,utotss
  real(KIND=DP) :: pi,alpha0,freq,lambda,kwave,omega,beta,Am,b
  integer :: it,itinit,itmax,isConverged1,i,j,k,file_step,last_step
  real(KIND=DP), dimension(Nx)::x
  real(KIND=DP), dimension(Ny)::y,mds,utots
  real(KIND=DP), dimension(Ny,Nx)::P,Pnew,u,v,S1,xm,ym,diff1,md,utot,perm
  real(KIND=DP), dimension(Nt)::tw
  real(KIND=DP), dimension(Nt,Ny,Nx)::alpha,dalphadt
  real(KIND=DP), dimension(Nx*Ny)::Pread
  real :: start, finish
  CHARACTER(40) numchr

  call cpu_time(start)


  file_step = 1
  last_step = 900000
  !!parameters
  !!change these as needed
  pi=4.0*atan(1.0)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
  !!!!!!!!!!!!!!!!!!!!!!!!!!grid and time step
  Ng=Nx*Ny
  x(1)=0.0
  y(1)=0.0
  !grid spacing
  dx=0.01
  dy=0.01
  dt=0.00002
  

  do i=2,Nx
  x(i)=x(i-1)+dx
  end do

  do i=2,Ny
  y(i)=y(i-1)+dy
  end do

  do i=2,Nt
  tw(i)=tw(i-1)+dt
  end do
!!!!!!!!!!!!!!!!!!!!mesh_grid
  do j=1,Ny
    do i=1,Nx
      xm(j,i)=x(i)
    end do
  end do
!!!!!!!!!!!!!!!!!
  do j=1,Ny
    do i=1,Nx
      ym(j,i)=y(j)
    end do
  end do
!!!!!alpha for delta waves

freq = 1.0 !1 Hz
lambda = 1000.0 !1000 \mu m or 1 mm
Am = 0.1

omega = 2*pi*freq*T
kwave = (2*pi/lambda)*L

    do k=1,Nt
    do j=1,Ny
    do i=1,Nx

    alpha(k,j,i)=alpha0*(1+Am*SIN(omega*tw(k)-kwave*y(j)-kwave*x(i)))

    dalphadt(k,j,i)=alpha0*omega*Am*COS(omega*tw(k)-kwave*y(j)-kwave*x(i))

    end do
    end do
    end do

!    print *,'this is dalphadt', dalphadt
   

! start time loop
    do k=1,Nt 


   !!!!!!!!!!!!!!!!!initialize
 P=0.0
 Pnew=0.0
 u=0.0
 v=0.0
 md=0.0
 weight=1.0
 diff1=0.0
 b=0.9
 !beta = (1-b)*(alpha0**3)/(kh*(1-alpha0)**2) !Dimensional beta
beta=0.1
 !!!!!!!!!!!!!!!!!!!
!!
!!!!!!!!!!!!!!!!!!!!sources and sinks - sets boundary condition on pressure at PVS
do j=1,Ny
do i=1,Nx

	perm(j,i) = (alpha(k,j,i)**3/(1.0_DP - alpha(k,j,i))**2) / &
	        (alpha0**3/(1.0_DP - alpha0)**2)
	if ((xm(j,i)**2 + (ym(j,i)-0.5)**2) < (r+lp)**2) then
	        P(j,i)=0.0
	else if (((xm(j,i)-1)**2 + (ym(j,i)-0.5)**2) < (r+lp)**2) then
	        P(j,i)=0.0
	else
	        P(j,i)=0.0
	end if

end do
end do

do j=1,Ny
do i=1,Nx

!S1(j,i)=((beta*alpha0**2)/(1-alph0**2))*dalphadt(k,j,i)*((1-alpha(k,j,i))**2)/(alpha(k,j,i)**3)
!S1(j,i) = ((1-b)*(alpha0**3)*(1-alpha(k,j,i))**2)*dalphadt(k,j,i)/(((1-alpha0)**2)*alpha(k,j,i)**3)
S1(j,i) = (beta/80)*dalphadt(k,j,i)*((1-alpha(k,j,i))**2)/(alpha(k,j,i)**3)

end do 
end do 
!!!!!!!!!!!!!!!!!!!!!!!! gausss_seidel parameters
 !Ap = 2.0*((1.0/dx**2) + (1.0/dy**2))
 !Ae = 1.0/(dx**2)
 !Aw = 1.0/dx**2
 !An = 1.0/dy**2
 !As = 1.0/dy**2

 itinit=0
 itmax=500000

 do it=itinit,itmax
!!call gausss_seidel
  call darcy_flow(r,lp,i,j,Ng,weight,x,y,Ap,Ae,Aw,An,As,S1,Pnew,P,diff1,&
  res1,res2,dx,dy,u,v,md,mds,mdss,xm,ym,utot,utots,utotss,perm)

  if ( mod(it,100000) == 0) then
  print *,'residual_1 is', res1,' residual_2 is', res2,'at iteration', it,&
   'and mean mass divergence', mdss
  end if
!!end iteration loop
 end do
!!write flow
    !do j=1,Ny
    !do i=1,Nx
    !   u(j,i) = u(j,i)
    !   v(j,i) = v(j,i)
    !end do
    !end do
 
 
if ( mod(k,file_step) == 0) then
print *,'time of wave travel is', k
   write(numchr,*) (k/file_step)
open(15,file="P"//trim(adjustl(numchr))//".txt")
write(15,*) P
close(15)

open(15,file="u"//trim(adjustl(numchr))//".txt")
write(15,*) u
close(15)

open(15,file="v"//trim(adjustl(numchr))//".txt")
write(15,*) v
close(15)

open(15,file="dalpha"//trim(adjustl(numchr))//".txt")
write(15,*) S1
close(15)
end if


    end do

!

call cpu_time(finish)

print *, "The real CPU time in seconds is ", finish-start

end program darcy_interstitium



subroutine darcy_flow(r,lp,i,j,Ng,weight,x,y,Ap,Ae,Aw,An,As,S1,Pnew,P,&
	diff1,res1,res2,dx,dy,u,v,md,mds,mdss,xm,ym,utot,utots,utotss,perm)
implicit none

integer, parameter :: DP = SELECTED_REAL_KIND(14)
integer, parameter :: Nx = 101
integer, parameter :: Ny = 101
integer, parameter :: Nt = 1000

real(KIND=DP) :: dx,dy,Ap,Ae,Aw,An,As,weight,res_sum,res1,res2,Ng,toler,r,lp,L,mdss,utotss
real(KIND=DP) :: kE,kW,kN,kS
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

res_sum=0.0_DP
do j=2,Ny-1
  do i=2,Nx-1
	  kE = 0.5_DP*(perm(j,i+1) + perm(j,i))
	  kW = 0.5_DP*(perm(j,i-1) + perm(j,i))
	  kN = 0.5_DP*(perm(j+1,i) + perm(j,i))
	  kS = 0.5_DP*(perm(j-1,i) + perm(j,i))
	  diff1(j,i) =((1.0/Ap)*(As*P(j-1,i)+An*P(j+1,i)+Ae*P(j,i+1)+Aw*P(j,i-1)&
	  -Ap*P(j,i)-S1(j,i)))**2
	  res_sum = res_sum + diff1(j,i)
 end do
end do
res2 = ((1/Ng)*res_sum)**0.5
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

if((xm(j,i)**2 + (ym(j,i)-0.5)**2) < (r+lp)**2) then
        Pnew(j,i)=0.0
else if (((xm(j,i)-1)**2 + (ym(j,i)-0.5)**2) < (r+lp)**2) then
        Pnew(j,i)=0.0
else
        Pnew(j,i)=Pnew(j,i)
end if

end do
end do

res1 = ((1/Ng)*(sum(diff1)))**0.5
!!update
!Pnew(1,1)=0.0_DP
P=Pnew

!!velocities
do j=1,Ny
    do i=2,Nx-1
        u(j,i) = -perm(j,i)*(P(j,i+1)-P(j,i-1))/(2*dx)
    end do
end do

do j=2,Ny-1
    do i=1,Nx
        v(j,i) = -perm(j,i)*(P(j+1,i)-P(j-1,i))/(2*dy)
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
