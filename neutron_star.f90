! Fortran 2008 implementation of a neutron-star
! equation-of-state and TOV solver.
!
! Input:
!   CSV file with header:
!   rho,X1,X08,X06,X04,X02,X00001
!
!   rho  : baryon number density [fm^-3]
!   X... : energy per baryon [MeV]
!
! Outputs:
!   P_Yp.dat       : rho, P [MeV/fm^3], Yp
!   mr_table.csv   : central density [g/cm^3], central pressure [erg/cm^3],
!                    mass [Msun], radius [km]
!   MR_curve.dat   : R [km], M [Msun], central density [g/cm^3],
!                    central pressure [erg/cm^3]
!   P_and_Yp.gp    : gnuplot script for pressure and proton fraction
!   MR_curve.gp    : gnuplot script for mass-radius results
!
! Compile:
!   gfortran -O2 -std=f2008 neutron_star.f90 -o neutron_star
!
! Run:
!   ./neutron_star data/example_eos.csv
!
! Plot:
!   gnuplot P_and_Yp.gp
!   gnuplot MR_curve.gp

module constants
  implicit none
  integer, parameter :: dp = kind(1.0d0)

  real(dp), parameter :: hbarc = 197.3269804_dp
  real(dp), parameter :: mn = 939.5654133_dp
  real(dp), parameter :: mp = 938.2720881_dp
  real(dp), parameter :: me = 0.51099895_dp
  real(dp), parameter :: mmu = 105.6583755_dp

  real(dp), parameter :: mevfm3_to_ergcm3 = 1.602176634d33
  real(dp), parameter :: mb = 1.66053907d-24
  real(dp), parameter :: G = 6.67430d-8
  real(dp), parameter :: c = 2.99792458d10
  real(dp), parameter :: Msun = 1.98892d33
  real(dp), parameter :: pi = 3.1415926535897932384626433832795_dp
end module constants


module csv_module
  use constants
  implicit none
contains

  subroutine replace_commas(line)
    character(len=*), intent(inout) :: line
    integer :: i
    do i = 1, len_trim(line)
       if (line(i:i) == ',') line(i:i) = ' '
    end do
  end subroutine replace_commas


  subroutine count_data_rows(filename, n)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: n
    integer :: u, ios
    character(len=4096) :: line
    real(dp) :: x1, x2, x3, x4, x5, x6, x7

    n = 0
    open(newunit=u, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'ERROR: cannot open input CSV.'

    read(u,'(A)',iostat=ios) line
    if (ios /= 0) error stop 'ERROR: empty input CSV.'

    do
       read(u,'(A)',iostat=ios) line
       if (ios /= 0) exit
       if (len_trim(line) == 0) cycle
       call replace_commas(line)
       read(line,*,iostat=ios) x1,x2,x3,x4,x5,x6,x7
       if (ios == 0) n = n + 1
    end do
    close(u)

    if (n < 4) error stop 'ERROR: at least four valid EOS rows are required.'
  end subroutine count_data_rows


  subroutine read_eos_csv(filename, rho, e)
    character(len=*), intent(in) :: filename
    real(dp), allocatable, intent(out) :: rho(:), e(:,:)
    integer :: n, u, i, ios
    character(len=4096) :: line

    call count_data_rows(filename, n)
    allocate(rho(n), e(6,n))

    open(newunit=u, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'ERROR: cannot open input CSV.'

    read(u,'(A)',iostat=ios) line
    i = 0

    do
       read(u,'(A)',iostat=ios) line
       if (ios /= 0) exit
       if (len_trim(line) == 0) cycle

       call replace_commas(line)
       read(line,*,iostat=ios) rho(i+1), e(1,i+1), e(2,i+1), e(3,i+1), &
                                e(4,i+1), e(5,i+1), e(6,i+1)
       if (ios == 0) i = i + 1
    end do
    close(u)

    if (i /= n) error stop 'ERROR: inconsistent CSV parsing.'
    if (any(rho <= 0.0_dp)) error stop 'ERROR: rho must be positive.'
  end subroutine read_eos_csv
end module csv_module


module spline_module
  use constants
  implicit none

  type spline_t
     integer :: n = 0
     real(dp), allocatable :: x(:), y(:), y2(:)
  end type spline_t

contains

  subroutine spline_build_natural(x,y,s)
    real(dp), intent(in) :: x(:),y(:)
    type(spline_t), intent(out) :: s
    integer :: n,i,k
    real(dp), allocatable :: u(:)
    real(dp) :: sig,p

    n=size(x)
    if (size(y)/=n .or. n<2) error stop 'ERROR: invalid spline input.'

    allocate(s%x(n),s%y(n),s%y2(n),u(n))
    s%n=n
    s%x=x
    s%y=y
    s%y2=0.0_dp
    u=0.0_dp

    do i=2,n-1
       if (x(i+1)<=x(i) .or. x(i)<=x(i-1)) &
            error stop 'ERROR: spline x-grid must be strictly increasing.'
       sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
       p=sig*s%y2(i-1)+2.0_dp
       s%y2(i)=(sig-1.0_dp)/p
       u(i)=(6.0_dp*((y(i+1)-y(i))/(x(i+1)-x(i)) - &
            (y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1)) - &
            sig*u(i-1))/p
    end do

    do k=n-1,1,-1
       s%y2(k)=s%y2(k)*s%y2(k+1)+u(k)
    end do

    deallocate(u)
  end subroutine spline_build_natural


  subroutine spline_build_notaknot(x,y,s)
    real(dp), intent(in) :: x(:),y(:)
    type(spline_t), intent(out) :: s
    integer :: n,i,j,k,piv
    real(dp), allocatable :: A(:,:),b(:),h(:)
    real(dp) :: maxv,tmp,fac

    n=size(x)
    if (size(y)/=n .or. n<2) error stop 'ERROR: invalid spline input.'

    allocate(s%x(n),s%y(n),s%y2(n))
    s%n=n
    s%x=x
    s%y=y
    s%y2=0.0_dp

    if (n==2) return

    allocate(A(n,n),b(n),h(n-1))

    do i=1,n-1
       h(i)=x(i+1)-x(i)
       if (h(i)<=0.0_dp) error stop 'ERROR: spline x-grid must be strictly increasing.'
    end do

    A=0.0_dp
    b=0.0_dp

    A(1,1)=-1.0_dp/h(1)
    A(1,2)=1.0_dp/h(1)+1.0_dp/h(2)
    A(1,3)=-1.0_dp/h(2)

    do i=2,n-1
       A(i,i-1)=h(i-1)
       A(i,i)=2.0_dp*(h(i-1)+h(i))
       A(i,i+1)=h(i)
       b(i)=6.0_dp*((y(i+1)-y(i))/h(i)-(y(i)-y(i-1))/h(i-1))
    end do

    A(n,n-2)=1.0_dp/h(n-2)
    A(n,n-1)=-(1.0_dp/h(n-2)+1.0_dp/h(n-1))
    A(n,n)=1.0_dp/h(n-1)

    do k=1,n-1
       piv=k
       maxv=abs(A(k,k))

       do i=k+1,n
          if (abs(A(i,k))>maxv) then
             maxv=abs(A(i,k))
             piv=i
          end if
       end do

       if (maxv<1.0d-30) error stop 'ERROR: singular not-a-knot spline system.'

       if (piv/=k) then
          do j=k,n
             tmp=A(k,j)
             A(k,j)=A(piv,j)
             A(piv,j)=tmp
          end do
          tmp=b(k)
          b(k)=b(piv)
          b(piv)=tmp
       end if

       do i=k+1,n
          fac=A(i,k)/A(k,k)
          A(i,k)=0.0_dp
          A(i,k+1:n)=A(i,k+1:n)-fac*A(k,k+1:n)
          b(i)=b(i)-fac*b(k)
       end do
    end do

    do i=n,1,-1
       if (i<n) then
          s%y2(i)=(b(i)-sum(A(i,i+1:n)*s%y2(i+1:n)))/A(i,i)
       else
          s%y2(i)=b(i)/A(i,i)
       end if
    end do

    deallocate(A,b,h)
  end subroutine spline_build_notaknot


  function spline_eval(s,xx) result(v)
    type(spline_t), intent(in) :: s
    real(dp), intent(in) :: xx
    real(dp) :: v,h,a,b
    integer :: klo,khi,k

    if (xx<s%x(1) .or. xx>s%x(s%n)) &
         error stop 'ERROR: spline evaluation outside its interpolation range.'

    klo=1
    khi=s%n

    do while(khi-klo>1)
       k=(khi+klo)/2
       if(s%x(k)>xx) then
          khi=k
       else
          klo=k
       end if
    end do

    h=s%x(khi)-s%x(klo)
    a=(s%x(khi)-xx)/h
    b=(xx-s%x(klo))/h

    v=a*s%y(klo)+b*s%y(khi)+ &
      ((a**3-a)*s%y2(klo)+(b**3-b)*s%y2(khi))*h*h/6.0_dp
  end function spline_eval


  function spline_deriv(s,xx) result(v)
    type(spline_t), intent(in) :: s
    real(dp), intent(in) :: xx
    real(dp) :: v,h,a,b
    integer :: klo,khi,k

    if (xx<s%x(1) .or. xx>s%x(s%n)) &
         error stop 'ERROR: spline derivative outside its interpolation range.'

    klo=1
    khi=s%n

    do while(khi-klo>1)
       k=(khi+klo)/2
       if(s%x(k)>xx) then
          khi=k
       else
          klo=k
       end if
    end do

    h=s%x(khi)-s%x(klo)
    a=(s%x(khi)-xx)/h
    b=(xx-s%x(klo))/h

    v=(s%y(khi)-s%y(klo))/h + &
      ((-3.0_dp*a*a+1.0_dp)*s%y2(klo) + &
       (3.0_dp*b*b-1.0_dp)*s%y2(khi))*h/6.0_dp
  end function spline_deriv
end module spline_module


module eos_module
  use constants
  use spline_module
  implicit none

  type eos_data
     integer :: n = 0
     real(dp), allocatable :: rho(:), delta(:)
     type(spline_t), allocatable :: espline(:)
  end type eos_data

contains

  subroutine build_eos(rho, e, eos)
    real(dp), intent(in) :: rho(:), e(:,:)
    type(eos_data), intent(out) :: eos
    real(dp), parameter :: xratio(6) = &
         [1.0_dp,0.8_dp,0.6_dp,0.4_dp,0.2_dp,0.0001_dp]
    real(dp) :: y(size(rho))
    integer :: j

    eos%n = size(rho)
    allocate(eos%rho(eos%n), eos%delta(6), eos%espline(6))
    eos%rho = rho

    do j=1,6
       eos%delta(j)=(1.0_dp-xratio(j))/(1.0_dp+xratio(j))
       y=e(j,:)
       call spline_build_natural(rho,y,eos%espline(j))
    end do
  end subroutine build_eos


  subroutine eA_fit_coeffs(eos,r,a0,a2,a4)
    type(eos_data), intent(in) :: eos
    real(dp), intent(in) :: r
    real(dp), intent(out) :: a0,a2,a4
    real(dp) :: M(3,3),b(3),d,ev,A(3,4)
    integer :: j

    M=0.0_dp
    b=0.0_dp

    do j=1,6
       d=eos%delta(j)
       ev=spline_eval(eos%espline(j),r)

       M(1,1)=M(1,1)+1.0_dp
       M(1,2)=M(1,2)+d**2
       M(1,3)=M(1,3)+d**4
       M(2,1)=M(2,1)+d**2
       M(2,2)=M(2,2)+d**4
       M(2,3)=M(2,3)+d**6
       M(3,1)=M(3,1)+d**4
       M(3,2)=M(3,2)+d**6
       M(3,3)=M(3,3)+d**8

       b(1)=b(1)+ev
       b(2)=b(2)+d**2*ev
       b(3)=b(3)+d**4*ev
    end do

    A(:,1:3)=M
    A(:,4)=b
    call solve3(A)

    a0=A(1,4)
    a2=A(2,4)
    a4=A(3,4)
  end subroutine eA_fit_coeffs


  subroutine solve3(A)
    real(dp), intent(inout) :: A(3,4)
    integer :: i,j,k,piv
    real(dp) :: maxv,tmp,fac

    do k=1,2
       piv=k
       maxv=abs(A(k,k))

       do i=k+1,3
          if (abs(A(i,k))>maxv) then
             maxv=abs(A(i,k))
             piv=i
          end if
       end do

       if (maxv<1.0d-14) error stop 'ERROR: singular 3x3 least-squares matrix.'

       if (piv/=k) then
          do j=k,4
             tmp=A(k,j)
             A(k,j)=A(piv,j)
             A(piv,j)=tmp
          end do
       end if

       do i=k+1,3
          fac=A(i,k)/A(k,k)
          A(i,k:4)=A(i,k:4)-fac*A(k,k:4)
       end do
    end do

    if (abs(A(3,3))<1.0d-14) &
         error stop 'ERROR: singular least-squares matrix.'

    A(3,4)=A(3,4)/A(3,3)
    A(2,4)=(A(2,4)-A(2,3)*A(3,4))/A(2,2)
    A(1,4)=(A(1,4)-A(1,2)*A(2,4)-A(1,3)*A(3,4))/A(1,1)
  end subroutine solve3


  function eA_of_delta(eos,r,d) result(v)
    type(eos_data), intent(in) :: eos
    real(dp), intent(in) :: r,d
    real(dp) :: v,a0,a2,a4

    call eA_fit_coeffs(eos,r,a0,a2,a4)
    v=a0+a2*d**2+a4*d**4
  end function eA_of_delta


  function deA_ddelta(eos,r,d) result(v)
    type(eos_data), intent(in) :: eos
    real(dp), intent(in) :: r,d
    real(dp) :: v,a0,a2,a4

    call eA_fit_coeffs(eos,r,a0,a2,a4)
    v=2.0_dp*a2*d+4.0_dp*a4*d**3
  end function deA_ddelta


  function fermi_kf(n) result(kf)
    real(dp), intent(in) :: n
    real(dp) :: kf

    if(n>0.0_dp) then
       kf=(3.0_dp*pi*pi*n)**(1.0_dp/3.0_dp)
    else
       kf=0.0_dp
    end if
  end function fermi_kf


  function n_from_mu(mu,mass) result(n)
    real(dp), intent(in) :: mu,mass
    real(dp) :: n,kf

    if(mu<=mass) then
       n=0.0_dp
    else
       kf=sqrt(mu*mu-mass*mass)/hbarc
       n=kf**3/(3.0_dp*pi*pi)
    end if
  end function n_from_mu


  subroutine lepton_eps_press(kf,mass,eps,pres)
    real(dp), intent(in) :: kf,mass
    real(dp), intent(out) :: eps,pres
    real(dp) :: mf,mu,x,logterm

    if(kf<=0.0_dp) then
       eps=0.0_dp
       pres=0.0_dp
       return
    end if

    mf=mass/hbarc
    mu=sqrt(kf*kf+mf*mf)
    x=kf
    logterm=log((x+mu)/mf)

    eps=(1.0_dp/(8.0_dp*pi*pi))* &
        (x*mu*(2.0_dp*x*x+mf*mf)-mf**4*logterm)*hbarc

    pres=(1.0_dp/(24.0_dp*pi*pi))* &
         (x*mu*(2.0_dp*x*x-3.0_dp*mf*mf)+ &
          3.0_dp*mf**4*logterm)*hbarc
  end subroutine lepton_eps_press


  function charge_balance(eos,r,d) result(f)
    type(eos_data), intent(in) :: eos
    real(dp), intent(in) :: r,d
    real(dp) :: f,Yp,np,mue,ne,nmu

    Yp=(1.0_dp-d)/2.0_dp
    np=r*Yp
    mue=max(2.0_dp*deA_ddelta(eos,r,d),0.0_dp)
    ne=n_from_mu(mue,me)
    nmu=n_from_mu(mue,mmu)

    f=np-ne-nmu
  end function charge_balance


  function beta_equilibrium_delta(eos,r) result(droot)
    type(eos_data), intent(in) :: eos
    real(dp), intent(in) :: r
    real(dp) :: droot,lo,hi,fl,fh,mid,fm
    integer :: iter

    lo=1.0d-6
    hi=0.999_dp
    fl=charge_balance(eos,r,lo)
    fh=charge_balance(eos,r,hi)

    if(fl*fh>0.0_dp) then
       if(abs(fl)<abs(fh)) then
          droot=lo
       else
          droot=hi
       end if
       return
    end if

    do iter=1,200
       mid=0.5_dp*(lo+hi)
       fm=charge_balance(eos,r,mid)

       if(abs(fm)<1.0d-12 .or. &
          abs(hi-lo)<1.0d-10*max(1.0_dp,abs(mid))) exit

       if(fl*fm<=0.0_dp) then
          hi=mid
          fh=fm
       else
          lo=mid
          fl=fm
       end if
    end do

    droot=0.5_dp*(lo+hi)
  end function beta_equilibrium_delta


  subroutine build_core_eos(eos,npoints,rhog,epscore,Yp,Pcore,thr_rho,thr_Yp)
    type(eos_data), intent(in) :: eos
    integer, intent(in) :: npoints
    real(dp), allocatable, intent(out) :: rhog(:),epscore(:),Yp(:),Pcore(:)
    real(dp), intent(out) :: thr_rho,thr_Yp
    type(spline_t) :: seps
    real(dp), allocatable :: deps(:)
    real(dp) :: d,yn,ypp,eint,mue,ne,nmu,epse,epsmu,pdummy
    real(dp) :: kfn,kfp,kfe
    integer :: i
    logical :: durca_found

    allocate(rhog(npoints),epscore(npoints),Yp(npoints),Pcore(npoints),deps(npoints))

    durca_found=.false.
    thr_rho=-1.0_dp
    thr_Yp=-1.0_dp

    do i=1,npoints
       rhog(i)=eos%rho(1)+real(i-1,dp)* &
            (eos%rho(eos%n)-eos%rho(1))/real(npoints-1,dp)

       d=beta_equilibrium_delta(eos,rhog(i))
       ypp=(1.0_dp-d)/2.0_dp
       yn=(1.0_dp+d)/2.0_dp
       eint=eA_of_delta(eos,rhog(i),d)

       mue=max(2.0_dp*deA_ddelta(eos,rhog(i),d),0.0_dp)
       ne=n_from_mu(mue,me)
       nmu=n_from_mu(mue,mmu)

       call lepton_eps_press(fermi_kf(ne),me,epse,pdummy)
       call lepton_eps_press(fermi_kf(nmu),mmu,epsmu,pdummy)

       epscore(i)=rhog(i)*(yn*mn+ypp*mp+eint)+epse+epsmu
       Yp(i)=ypp

       kfn=fermi_kf(rhog(i)*yn)
       kfp=fermi_kf(rhog(i)*ypp)
       kfe=fermi_kf(ne)

       if(.not.durca_found .and. (kfn<=kfp+kfe)) then
          thr_rho=rhog(i)
          thr_Yp=ypp
          write(*,'(A,F10.5,A,F10.5,A)') &
               '>> Direct Urca opens at rho = ',rhog(i), &
               ' fm^-3, Yp = ',ypp,' (with muons included)'
          durca_found=.true.
       end if
    end do

    call spline_build_notaknot(rhog,epscore,seps)

    do i=1,npoints
       deps(i)=spline_deriv(seps,rhog(i))
       Pcore(i)=rhog(i)*deps(i)-epscore(i)
    end do

    deallocate(deps)
  end subroutine build_core_eos
end module eos_module


module tov_module
  use constants
  use spline_module
  implicit none

  type full_eos_t
     integer :: n=0
     real(dp), allocatable :: rho(:),eps(:),P(:)
     type(spline_t) :: eps_log_spline
     type(spline_t) :: P_logrho_spline
  end type full_eos_t

contains

  subroutine build_full_eos(rhog,epscore,Pcore,eosfull)
    real(dp), intent(in) :: rhog(:),epscore(:),Pcore(:)
    type(full_eos_t), intent(out) :: eosfull
    integer :: ncore,nlow,i,nfull,nkeep,j,k
    real(dp) :: r0,P0,eps0,gammaeff,Klow
    real(dp), allocatable :: rlow(:),Plow(:),elow(:)
    real(dp), allocatable :: rm(:),em(:),Pm(:)
    real(dp) :: r1n,r2n,P1n,P2n
    real(dp) :: r1_,r2_,P1_,P2_,e2,epred,int1,int2
    logical, allocatable :: good(:)

    ncore=size(rhog)

    r0=rhog(1)*1.0d39*mb
    P0=Pcore(1)*mevfm3_to_ergcm3
    eps0=epscore(1)*mevfm3_to_ergcm3

    if(P0<=0.0_dp) error stop 'ERROR: first EOS pressure must be positive.'

    r1n=rhog(1)
    r2n=rhog(4)
    P1n=Pcore(1)
    P2n=Pcore(4)

    if(P2n<=P1n .or. P1n<=0.0_dp) &
         error stop 'ERROR: low-density EOS requires positive increasing pressure.'

    gammaeff=log(P2n/P1n)/log(r2n/r1n)
    Klow=P0/r0**gammaeff

    nlow=500
    allocate(rlow(nlow),Plow(nlow),elow(nlow))

    do i=1,nlow
       rlow(i)=1.0d4*(r0/1.0d4)**(real(i-1,dp)/real(nlow-1,dp))
       Plow(i)=Klow*rlow(i)**gammaeff
    end do

    elow(nlow)=eps0

    do i=nlow-1,1,-1
       r1_=rlow(i)
       r2_=rlow(i+1)
       P1_=Plow(i)
       P2_=Plow(i+1)
       e2=elow(i+1)

       int2=(e2+P2_)/r2_
       epred=e2-int2*(r2_-r1_)
       int1=(epred+P1_)/r1_
       elow(i)=e2-0.5_dp*(int1+int2)*(r2_-r1_)
    end do

    nfull=(nlow-1)+ncore

    allocate(rm(nfull),em(nfull),Pm(nfull))
    rm(1:nlow-1)=rlow(1:nlow-1)
    em(1:nlow-1)=elow(1:nlow-1)
    Pm(1:nlow-1)=Plow(1:nlow-1)

    rm(nlow:nfull)=rhog*1.0d39*mb
    em(nlow:nfull)=epscore*mevfm3_to_ergcm3
    Pm(nlow:nfull)=Pcore*mevfm3_to_ergcm3

    allocate(good(nfull))
    good=.false.
    good(1)=.true.

    do i=2,nfull
       good(i)=(Pm(i)>Pm(i-1))
    end do

    nkeep=count(good)

    if(nkeep<4) error stop 'ERROR: insufficient monotonic EOS points.'

    allocate(eosfull%rho(nkeep),eosfull%eps(nkeep),eosfull%P(nkeep))
    eosfull%n=nkeep

    k=0
    do j=1,nfull
       if(good(j)) then
          k=k+1
          eosfull%rho(k)=rm(j)
          eosfull%eps(k)=em(j)
          eosfull%P(k)=Pm(j)
       end if
    end do

    call spline_build_notaknot(log(eosfull%P),log(eosfull%eps), &
                               eosfull%eps_log_spline)
    call spline_build_notaknot(log(eosfull%rho),log(eosfull%P), &
                               eosfull%P_logrho_spline)

    deallocate(rlow,Plow,elow,rm,em,Pm,good)
  end subroutine build_full_eos


  function eps_from_P(eos,Pval) result(eps)
    type(full_eos_t), intent(in) :: eos
    real(dp), intent(in) :: Pval
    real(dp) :: eps

    if(Pval<=eos%P(1)) then
       eps=eos%eps(1)*(Pval/eos%P(1))
    else if(Pval>=eos%P(eos%n)) then
       eps=exp(spline_eval(eos%eps_log_spline,log(eos%P(eos%n))))
    else
       eps=exp(spline_eval(eos%eps_log_spline,log(Pval)))
    end if
  end function eps_from_P


  function central_pressure(eos,rhoc) result(Pc)
    type(full_eos_t), intent(in) :: eos
    real(dp), intent(in) :: rhoc
    real(dp) :: Pc,lr

    lr=log(rhoc)

    if(lr<log(eos%rho(1)) .or. lr>log(eos%rho(eos%n))) &
         error stop 'ERROR: central density outside EOS range.'

    Pc=exp(spline_eval(eos%P_logrho_spline,lr))
  end function central_pressure


  subroutine tov_rhs(eos,r,m,p,dm,dpdr)
    type(full_eos_t), intent(in) :: eos
    real(dp), intent(in) :: r,m,p
    real(dp), intent(out) :: dm,dpdr
    real(dp) :: e,den,num

    if(p<=eos%P(1) .or. r<1.0d-6) then
       dm=0.0_dp
       dpdr=0.0_dp
       return
    end if

    e=eps_from_P(eos,p)/c**2
    dm=4.0_dp*pi*r*r*e

    num=G*(e+p/c**2)*(m+4.0_dp*pi*r**3*p/c**2)
    den=r*r*(1.0_dp-2.0_dp*G*m/(r*c**2))

    dpdr=-num/den
  end subroutine tov_rhs


  subroutine rk45_dp_step(eos,r,h,m,p,m5,p5,errm,errp)
    type(full_eos_t), intent(in) :: eos
    real(dp), intent(in) :: r,h,m,p
    real(dp), intent(out) :: m5,p5,errm,errp
    real(dp) :: k1m,k1p,k2m,k2p,k3m,k3p,k4m,k4p
    real(dp) :: k5m,k5p,k6m,k6p,k7m,k7p
    real(dp) :: mt,pt

    call tov_rhs(eos,r,m,p,k1m,k1p)

    mt=m+h*(1.0_dp/5.0_dp*k1m)
    pt=p+h*(1.0_dp/5.0_dp*k1p)
    call tov_rhs(eos,r+h/5.0_dp,mt,pt,k2m,k2p)

    mt=m+h*(3.0_dp/40.0_dp*k1m+9.0_dp/40.0_dp*k2m)
    pt=p+h*(3.0_dp/40.0_dp*k1p+9.0_dp/40.0_dp*k2p)
    call tov_rhs(eos,r+3.0_dp*h/10.0_dp,mt,pt,k3m,k3p)

    mt=m+h*(44.0_dp/45.0_dp*k1m-56.0_dp/15.0_dp*k2m+ &
            32.0_dp/9.0_dp*k3m)
    pt=p+h*(44.0_dp/45.0_dp*k1p-56.0_dp/15.0_dp*k2p+ &
            32.0_dp/9.0_dp*k3p)
    call tov_rhs(eos,r+4.0_dp*h/5.0_dp,mt,pt,k4m,k4p)

    mt=m+h*(19372.0_dp/6561.0_dp*k1m-25360.0_dp/2187.0_dp*k2m+ &
            64448.0_dp/6561.0_dp*k3m-212.0_dp/729.0_dp*k4m)
    pt=p+h*(19372.0_dp/6561.0_dp*k1p-25360.0_dp/2187.0_dp*k2p+ &
            64448.0_dp/6561.0_dp*k3p-212.0_dp/729.0_dp*k4p)
    call tov_rhs(eos,r+8.0_dp*h/9.0_dp,mt,pt,k5m,k5p)

    mt=m+h*(9017.0_dp/3168.0_dp*k1m-355.0_dp/33.0_dp*k2m+ &
            46732.0_dp/5247.0_dp*k3m+49.0_dp/176.0_dp*k4m- &
            5103.0_dp/18656.0_dp*k5m)
    pt=p+h*(9017.0_dp/3168.0_dp*k1p-355.0_dp/33.0_dp*k2p+ &
            46732.0_dp/5247.0_dp*k3p+49.0_dp/176.0_dp*k4p- &
            5103.0_dp/18656.0_dp*k5p)
    call tov_rhs(eos,r+h,mt,pt,k6m,k6p)

    mt=m+h*(35.0_dp/384.0_dp*k1m+500.0_dp/1113.0_dp*k3m+ &
            125.0_dp/192.0_dp*k4m-2187.0_dp/6784.0_dp*k5m+ &
            11.0_dp/84.0_dp*k6m)
    pt=p+h*(35.0_dp/384.0_dp*k1p+500.0_dp/1113.0_dp*k3p+ &
            125.0_dp/192.0_dp*k4p-2187.0_dp/6784.0_dp*k5p+ &
            11.0_dp/84.0_dp*k6p)
    call tov_rhs(eos,r+h,mt,pt,k7m,k7p)

    m5=m+h*(35.0_dp/384.0_dp*k1m+500.0_dp/1113.0_dp*k3m+ &
            125.0_dp/192.0_dp*k4m-2187.0_dp/6784.0_dp*k5m+ &
            11.0_dp/84.0_dp*k6m)

    p5=p+h*(35.0_dp/384.0_dp*k1p+500.0_dp/1113.0_dp*k3p+ &
            125.0_dp/192.0_dp*k4p-2187.0_dp/6784.0_dp*k5p+ &
            11.0_dp/84.0_dp*k6p)

    errm=h*((35.0_dp/384.0_dp-5179.0_dp/57600.0_dp)*k1m+ &
           (500.0_dp/1113.0_dp-7571.0_dp/16695.0_dp)*k3m+ &
           (125.0_dp/192.0_dp-393.0_dp/640.0_dp)*k4m+ &
           (-2187.0_dp/6784.0_dp+92097.0_dp/339200.0_dp)*k5m+ &
           (11.0_dp/84.0_dp-187.0_dp/2100.0_dp)*k6m- &
           1.0_dp/40.0_dp*k7m)

    errp=h*((35.0_dp/384.0_dp-5179.0_dp/57600.0_dp)*k1p+ &
           (500.0_dp/1113.0_dp-7571.0_dp/16695.0_dp)*k3p+ &
           (125.0_dp/192.0_dp-393.0_dp/640.0_dp)*k4p+ &
           (-2187.0_dp/6784.0_dp+92097.0_dp/339200.0_dp)*k5p+ &
           (11.0_dp/84.0_dp-187.0_dp/2100.0_dp)*k6p- &
           1.0_dp/40.0_dp*k7p)
  end subroutine rk45_dp_step


  subroutine integrate_star(eos,Pc,Msolar,Rkm)
    type(full_eos_t), intent(in) :: eos
    real(dp), intent(in) :: Pc
    real(dp), intent(out) :: Msolar,Rkm
    real(dp), parameter :: rstart=1.0_dp
    real(dp), parameter :: rmax=3.0d6
    real(dp), parameter :: hmax=1.0d3
    real(dp), parameter :: rtol=1.0d-8
    real(dp) :: r,h,m,p,mnew,pnew,em,ep,scale_m,scale_p,errnorm,fac
    real(dp) :: rs,ms,ps,rt,mt,pt,fs,ft,theta,r_event,m_event
    real(dp) :: pmin,surf
    integer :: iter,maxiter

    pmin=eos%P(1)
    surf=pmin*1.0001_dp

    r=rstart
    m=4.0_dp*pi*r**3*(eps_from_P(eos,Pc)/c**2)/3.0_dp
    p=Pc
    h=hmax
    maxiter=100000

    do iter=1,maxiter
       if(p<=surf .or. r>=rmax) exit

       h=min(h,hmax,rmax-r)

       call rk45_dp_step(eos,r,h,m,p,mnew,pnew,em,ep)

       scale_m=1.0d10+rtol*max(abs(m),abs(mnew))
       scale_p=pmin*1.0d-6+rtol*max(abs(p),abs(pnew))

       errnorm=sqrt(0.5_dp*((em/scale_m)**2+(ep/scale_p)**2))

       if(errnorm<=1.0_dp .or. h<=1.0d-6) then
          rs=r
          ms=m
          ps=p
          rt=r+h
          mt=mnew
          pt=pnew

          if(ps>surf .and. pt<=surf) then
             fs=ps-surf
             ft=pt-surf
             theta=fs/(fs-ft)
             r_event=rs+theta*(rt-rs)
             m_event=ms+theta*(mt-ms)

             Msolar=m_event/Msun
             Rkm=r_event/1.0d5
             return
          end if

          r=rt
          m=mt
          p=pt

          if(1.0_dp-2.0_dp*G*m/(r*c**2)<=0.0_dp) exit
       end if

       if(errnorm==0.0_dp) then
          fac=5.0_dp
       else
          fac=0.9_dp*errnorm**(-0.2_dp)
          fac=max(0.2_dp,min(5.0_dp,fac))
       end if

       h=min(hmax,h*fac)
    end do

    Msolar=m/Msun
    Rkm=r/1.0d5
  end subroutine integrate_star


  subroutine solve_tov_sequence(eos,rho_core_min,nstars,mr)
    type(full_eos_t), intent(in) :: eos
    real(dp), intent(in) :: rho_core_min
    integer, intent(in) :: nstars
    real(dp), allocatable, intent(out) :: mr(:,:)
    real(dp) :: rhoc,frac,Pc,M,R
    integer :: i

    allocate(mr(nstars,4))

    do i=1,nstars
       frac=real(i-1,dp)/real(nstars-1,dp)
       rhoc=(rho_core_min*1.001_dp)* &
            ((eos%rho(eos%n)*0.999_dp)/(rho_core_min*1.001_dp))**frac

       Pc=central_pressure(eos,rhoc)
       call integrate_star(eos,Pc,M,R)

       mr(i,:)=[rhoc,Pc,M,R]
    end do
  end subroutine solve_tov_sequence
end module tov_module


program neutron_star_pipeline
  use constants
  use csv_module
  use eos_module
  use tov_module
  implicit none

  character(len=1024) :: filename
  integer :: nargs,n,npoints,i,ios
  real(dp), allocatable :: rho(:),etable(:,:)
  real(dp), allocatable :: rhog(:),epscore(:),Yp(:),Pcore(:)
  real(dp), allocatable :: mr(:,:)
  real(dp) :: rho_core_min,thr_rho,thr_Yp
  integer :: imax,iu
  real(dp) :: R14,Mmax

  nargs=command_argument_count()

  if(nargs>=1) then
     call get_command_argument(1,filename)
  else
     filename='data/example_eos.csv'
  end if

  write(*,'(A)') 'Reading EOS table: '//trim(filename)
  call read_eos_csv(trim(filename),rho,etable)

  n=size(rho)
  write(*,'(A,I0)') 'EOS rows: ',n

  block
     type(eos_data) :: eos
     type(full_eos_t) :: eosfull

     call build_eos(rho,etable,eos)

     npoints=121
     call build_core_eos(eos,npoints,rhog,epscore,Yp,Pcore,thr_rho,thr_Yp)

     open(newunit=iu,file='P_Yp.dat',status='replace',action='write',iostat=ios)
     if(ios/=0) error stop 'ERROR: cannot create P_Yp.dat.'

     write(iu,'(A)') '# rho(fm^-3)  P(MeV/fm^3)  Yp'

     do i=1,npoints
        write(iu,'(3(ES24.15,1X))') rhog(i),Pcore(i),Yp(i)
     end do
     close(iu)

     call build_full_eos(rhog,epscore,Pcore,eosfull)
     rho_core_min=rhog(1)*1.0d39*mb

     call solve_tov_sequence(eosfull,rho_core_min,80,mr)

     imax=maxloc(mr(:,3),dim=1)
     Mmax=mr(imax,3)

     write(*,'(A,F10.5,A,ES12.5,A,F10.4,A)') &
          'Mmax = ',Mmax,' Msun at rho_c = ',mr(imax,1), &
          ' g/cm^3, R = ',mr(imax,4),' km'

     call interpolate_radius_at_mass(mr,imax,1.4_dp,R14)

     if(R14>0.0_dp) then
        write(*,'(A,F10.4,A)') 'R(1.4 Msun) = ',R14,' km'
     end if

     open(newunit=iu,file='mr_table.csv',status='replace',action='write',iostat=ios)
     if(ios/=0) error stop 'ERROR: cannot create mr_table.csv.'

     write(iu,'(A)') 'rho_c[g/cm3],P_c[erg/cm3],M[Msun],R[km]'

     do i=1,size(mr,1)
        write(iu,'(4(ES24.15,1X))') mr(i,1),mr(i,2),mr(i,3),mr(i,4)
     end do
     close(iu)

     open(newunit=iu,file='MR_curve.dat',status='replace',action='write',iostat=ios)
     if(ios/=0) error stop 'ERROR: cannot create MR_curve.dat.'

     write(iu,'(A)') '# R(km) M(Msun) rho_c(g/cm3) Pc(erg/cm3)'

     do i=1,size(mr,1)
        write(iu,'(4(ES24.15,1X))') mr(i,4),mr(i,3),mr(i,1),mr(i,2)
     end do
     close(iu)

     call write_gnuplot_scripts(thr_Yp)

     deallocate(eos%rho,eos%delta,eos%espline)
     deallocate(eosfull%rho,eosfull%eps,eosfull%P)
  end block

  deallocate(rho,etable,rhog,epscore,Yp,Pcore,mr)

  write(*,'(A)') 'Done.'
  write(*,'(A)') 'Run: gnuplot P_and_Yp.gp'
  write(*,'(A)') 'Run: gnuplot MR_curve.gp'

contains

  subroutine interpolate_radius_at_mass(mr,imax,target,R)
    real(dp), intent(in) :: mr(:,:),target
    integer, intent(in) :: imax
    real(dp), intent(out) :: R
    integer :: i
    real(dp) :: m1,m2,r1,r2

    R=0.0_dp

    do i=1,imax-1
       m1=mr(i,3)
       m2=mr(i+1,3)

       if((target-m1)*(target-m2)<=0.0_dp .and. m2/=m1) then
          r1=mr(i,4)
          r2=mr(i+1,4)
          R=r1+(target-m1)*(r2-r1)/(m2-m1)
          return
       end if
    end do

    write(*,'(A)') 'WARNING: M=1.4 Msun is outside the stable sequence.'
  end subroutine interpolate_radius_at_mass


  subroutine write_gnuplot_scripts(threshold_yp)
    real(dp), intent(in) :: threshold_yp
    integer :: u,ios

    open(newunit=u,file='P_and_Yp.gp',status='replace',action='write',iostat=ios)
    if(ios/=0) error stop 'ERROR: cannot create P_and_Yp.gp.'

    write(u,'(A)') 'set terminal pngcairo size 1050,420 enhanced'
    write(u,'(A)') 'set output "P_and_Yp.png"'
    write(u,'(A)') 'set multiplot layout 1,2'
    write(u,'(A)') 'set xlabel "rho (fm^{-3})"'
    write(u,'(A)') 'set ylabel "P (MeV fm^{-3})"'
    write(u,'(A)') 'set title "Pressure in beta equilibrium"'
    write(u,'(A)') 'plot "P_Yp.dat" using 1:2 with lines lw 2 notitle'

    write(u,'(A)') 'set xlabel "rho (fm^{-3})"'
    write(u,'(A)') 'set ylabel "Yp = Z/A"'
    write(u,'(A)') 'set title "Beta-equilibrium proton fraction"'

    if(threshold_yp>0.0_dp) then
       write(u,'(A,F8.5,A,F8.5,A)') &
            'set arrow from graph 0,first ',threshold_yp, &
            ' to graph 1,first ',threshold_yp,' nohead dt 2 lc rgb "black"'
    end if

    write(u,'(A)') 'plot "P_Yp.dat" using 1:3 with lines lw 2 notitle'
    write(u,'(A)') 'unset multiplot'
    close(u)

    open(newunit=u,file='MR_curve.gp',status='replace',action='write',iostat=ios)
    if(ios/=0) error stop 'ERROR: cannot create MR_curve.gp.'

    write(u,'(A)') 'set terminal pngcairo size 1100,450 enhanced'
    write(u,'(A)') 'set output "MR_curve.png"'
    write(u,'(A)') 'set multiplot layout 1,2'

    write(u,'(A)') 'set xlabel "R (km)"'
    write(u,'(A)') 'set ylabel "M (M_{sun})"'
    write(u,'(A)') 'set title "Neutron-star mass-radius relation"'
    write(u,'(A)') 'set yrange [0:2.4]'
    write(u,'(A)') 'plot "MR_curve.dat" using 1:2 with linespoints pt 7 ps 0.45 notitle'

    write(u,'(A)') 'set xlabel "Central density rho_c (g/cm^3)"'
    write(u,'(A)') 'set ylabel "M (M_{sun})"'
    write(u,'(A)') 'set title "Mass vs central density"'
    write(u,'(A)') 'set xrange [0:2.0e15]'
    write(u,'(A)') 'set xtics ("1x10^{15}" 1.0e15, "2x10^{15}" 2.0e15)'
    write(u,'(A)') 'plot "MR_curve.dat" using 3:2 with linespoints pt 7 ps 0.45 notitle'

    write(u,'(A)') 'unset multiplot'
    close(u)
  end subroutine write_gnuplot_scripts

end program neutron_star_pipeline
