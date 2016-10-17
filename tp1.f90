program TP
  implicit none
  integer,parameter   :: Np=500000
  real                :: Tint,R,T,dt,tau,time,tmax,T0
  integer(kind=4)     :: seed
  real, dimension(Np) :: Ep,Ep0
  integer             :: i,n

  call srand(seed)
  
  Tint = 0.
  R    = 280.
  tau  = 5.
  dt   = 1.e-2
  T    = 300.
  
  time = 0.
  tmax = 3*tau
  
  do i=1,Np
     call generateEnergy(100000.,200000.,Ep(i))
     Tint = Tint + Ep(i)
  end do
  Tint = Tint/(Np*R)
  T0 = Tint
  Ep0 = Ep

  open(unit=10,file="energy.dat")

  do while(time<tmax)
     time = time + dt
     
     call collision()
     
     write(10,*) time, moyenne(Ep)/R, moyenne(Ep)/R-T-exp(-2*time/tau)*(T0-T)
     
     Ep0 = Ep

  end do
  
  close(10)

  call histogram(Ep,100)

contains
  subroutine generateEnergy(min, max, E)
    real,intent(in)  :: min,max
    real,intent(out) :: E
    
    E = min + rand()*(max-min)
    
  end subroutine 

  subroutine collision()!rho,T1,part_m,Npart_m) 
    !real                       :: rho
    !real,dimension(Np)         :: T1
    !integer                    :: Npart_m
    !integer,dimension(Npart_m) :: part_m 
    real :: sig

    do i=1,Np
       sig = gaussian(0.,1.)
       Ep(i)=(Ep0(i)+R*T*dt*(1+sig**2)/tau+2*sqrt(dt*R*T*Ep0(i)/tau)*sig)&
         /(1+2*dt/tau)
    end do


  end subroutine collision

  function moyenne(Tab)
    real, dimension(:) :: Tab
    real               :: moyenne
    integer :: i
    
    moyenne = 0.
    do i=1,size(Tab)
       moyenne = moyenne + Tab(i)
    enddo
    moyenne = moyenne / size(Tab) 
  end function moyenne

  function gaussian(esp, var)
    real :: esp, var
    real, parameter :: pi = 3.14159265359
    real :: a, b
    real :: gaussian

    a = 1-rand()
    b = rand()
    
    gaussian = sqrt(-2*log(a)) * cos(2*pi*b)*var + esp

  end function gaussian

  subroutine histogram(Tab,nb_int)
    real,dimension(:) :: Tab
    integer :: nb_int, i, j
    real    :: max, min
    real,dimension(:,:),allocatable :: Val

    allocate(Val(nb_int,2))

    Val=0.

    min=Tab(1)
    max=Tab(1)
    do i=2,size(Tab)
       if (Tab(i)<min)then
          min=Tab(i)
       elseif(Tab(i)>max)then
          max=Tab(i)
       end if
    end do


    do i=1,nb_int
       Val(i,1)=min + i*(max-min)/nb_int
    end do


    do i=1,size(Tab)
       do j=1,nb_int
          if(Tab(i)<Val(j,1))then
             Val(j,2)=Val(j,2)+1
             exit
          endif
       end do
    end do

    open(unit=11,file="histogram.dat")
    do i=1,nb_int
       write(11,*) Val(i,1), Val(i,2)
    end do
    close(11)
    deallocate(Val)
  end subroutine histogram


end program TP
