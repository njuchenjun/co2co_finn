#define pesdir "./"
! --------------------------------------------------------------------
! fit to UCCSD(T)-F12a/avtz energy
! input C-O1 C-O2, O1-O2 in Bohr, or Cartesian of COO in Bohr
! return energy in Hartree, with v(OH+HCl)=0.d0
! --------------------------------------------------------------------
! program test
! implicit none
! real*8 r(3),v0,v
! do while(.true.)
!   read(*,*,end=99) r,v0
!   call pot_co2(1,r,v); v=v*27.21138505d0
!   write(*,'(10f10.5)') r,v0,v,v0-v
! enddo
! 99 continue
! end program
! --------------------------------------------------------------------
subroutine pot_co2_bba(r1,r2,th,v)
  implicit none
  real*8,intent(in) :: r1,r2,th ! bohr , degree
  real*8,intent(out) :: v ! hartree
  real*8 :: dist(3,1),pot(1)
  real*8,parameter :: pi=dacos(-1.d0)
  dist(1,1)=r1
  dist(2,1)=r2
  dist(3,1)=dsqrt(r1**2+r2**2-2.d0*r1*r2*dcos(th/180.d0*pi))
  call pot_co2(1,dist,pot)
  v=pot(1)
  return
end subroutine
! --------------------------------------------------------------------
subroutine pot_co2_cart(x,v)
  implicit none
  real*8,intent(in)  :: x(3,3)
  real*8,intent(out) :: v
  real*8 :: r(3),d(3)
  d=x(1:3,1)-x(1:3,2); r(1)=dsqrt(dot_product(d,d))
  d=x(1:3,1)-x(1:3,3); r(2)=dsqrt(dot_product(d,d))
  d=x(1:3,2)-x(1:3,3); r(3)=dsqrt(dot_product(d,d))
  call pot_co2(1,r,v)
  write(444,*) x,v
  return
end subroutine pot_co2_cart
! --------------------------------------------------------------------
subroutine pot_co2_cart1(x,vco2,dvco2)
  implicit none
  real*8,intent(in)  :: x(3,3)  !OOC
  real*8,intent(out) :: vco2,dvco2(3)
  real*8 :: v(7)
  real*8 :: x_ooc(3,3)
  real*8 :: r(3),d(3),rd(3),rc(3,7)
  real*8,parameter :: dr=1.d-2
  real*8,parameter :: cm=219474.5d0
  integer :: l

  x_ooc=x(:,[3,1,2])
  d=x_ooc(1:3,1)-x_ooc(1:3,2); r(1)=dsqrt(dot_product(d,d))
  d=x_ooc(1:3,1)-x_ooc(1:3,3); r(2)=dsqrt(dot_product(d,d))
  d=x_ooc(1:3,2)-x_ooc(1:3,3); r(3)=dsqrt(dot_product(d,d))

  !write(441,*) x,r

  call pot_co2(1,r,v(1))
  do l=1,3
     rc(:,l*2)=r
     rc(l,l*2)=rc(l,l*2)-dr
     rd=rc(:,l*2)
     call pot_co2(1,rd,v(l*2))
     rc(:,l*2+1)=r
     rc(l,l*2+1)=rc(l,l*2+1)+dr
     rd=rc(:,l*2+1)
     call pot_co2(1,rd,v(l*2+1))
   end do

  !call pot_co2(1,rc,v) 
  
  vco2=v(1) ! vco2 : a.u

  dvco2(1)=(v(3)-v(2))/2.0d0/dr
  dvco2(2)=(v(5)-v(4))/2.0d0/dr
  dvco2(3)=(v(7)-v(6))/2.0d0/dr

  !write(442,*) v,dvco2

  return
end subroutine pot_co2_cart1
! --------------------------------------------------------------------
subroutine pot_co2(n,r,v)
  implicit none
  integer,intent(in) :: n
  real*8,intent(in)  :: r(3,n)
  real*8,intent(out) :: v(n)
  real*8 :: rx(3,n),v1(n),v2(n),w1(n),w2(n)
  integer i
  real*8 :: logsig,x
  logsig(x)=1.d0/(1.d0+dexp(-x))

  rx=r; call co2_nnfit_36(3,n,rx,v1) ! version 2021/10/11, 3947 points
  rx(1,1:n)=r(2,1:n); rx(2,1:n)=r(1,1:n); call co2_nnfit_36(3,n,rx,v2)

  do i=1,n
    w1(i)=logsig(20.d0*(r(2,i)-r(1,i)))
    w2(i)=1.d0-w1(i)
  enddo
  v=v1*w1+v2*w2
  v=v/27.21138505d0
  return
contains
! --------------------------------------------------------------------
  subroutine nsim(r0,v,idv,dv,n0,n1,n2,rg,w1,b1,w2,b2,w3,b3,vg)
  implicit none
  integer,intent(in) :: n0,n1,n2,idv
  real*8,intent(in)  :: r0(n0),rg(2,n0),vg(2)
  real*8,intent(in)  :: w1(n0,n1),b1(n1)
  real*8,intent(in)  :: w2(n1,n2),b2(n2)
  real*8,intent(in)  :: w3(n2),b3
  real*8,intent(out) :: v,dv(n0)
  integer :: i,j,k
  real*8  :: r(n0),rgg(n0),vgg,ax(n1),bx(n2)
  real*8  :: dvtm,rtmp,rt1(n1),rt2(n2)
  real*8,external :: ddot
  v=0.d0
  r=r0
  vgg=vg(2)-vg(1)
  ! mapminmax [-1,1]
  do i=1,n0
    rgg(i)=rg(2,i)-rg(1,i)
    r(i)=2.d0*(r(i)-rg(1,i))/rgg(i)-1.d0
  end do
  ! 1st layer
  rt1=b1
  call dgemv('t',n0,n1,1.d0,w1,n0,r,1,1.d0,rt1,1)
  ax=dtanh(rt1)
  ! 2nd layer
  rt2=b2
  call dgemv('t',n1,n2,1.d0,w2,n1,ax,1,1.d0,rt2,1)
  bx=dtanh(rt2)
  ! output layer
  v=b3+ddot(n2,w3,1,bx,1)
  !reverse map
  v=vgg*(v+1.d0)/2+vg(1)
  if(idv.ne.1) return
  ! calculate first derivatives, dv(i)=dv/dr(i)
  dv=0.d0
  do i=1,n0
    do k=1,n2
      dvtm=0.d0
      do j=1,n1
        dvtm=dvtm+w2(j,k)*w1(i,j)*(1-ax(j)**2)
      enddo
      dv(i)=dv(i)+w3(k)*dvtm*(1-bx(k)**2)
    enddo
    dv(i)=dv(i)*vgg/rgg(i)
  enddo
  return
  end subroutine nsim

! --------------------------------------------------------------------
!  simulate a neural network with two hidden layers
!     n0-n1(tansig)-n2(tansig)-1(pureline)
!  blas routines used in this subroutine: dgemm, dgemv
!  the syntax keeps the same with nsim, except input nt at first
!  Chen Jun, 2013/11/26
!
  subroutine nsimx(nt,r,v,idv,dv,n0,n1,n2,rg,w1,b1,w2,b2,w3,b3,vg)
!
!  description of parameters
!  nt  : total number of points in r,v and dv
!  r   : input of nn, r(n0,nt), r(1:n0,i) is the i-th input
!  v   : output of nn, v(nt), v(i) is the i-th input
!  idv : if idv==0, calc. v only, else calc. v and derivatives dv
!  dv  : derivatives of dv/dr, dv(n0,nt), dv(j,i)=d v(i) / d r(j,i)
!  n0  : neurons of the input layer, dimension of r
!  n1  : neurons of the first hidden layer
!  n2  : neurons of the second hidden layer
!  rg  : input ranges of the training set, rg(2,n0)
!  w1  : weights connect input layer with 1-st hidden layer
!  b1  : biases of 1-st hidden layer
!  w2  : weights connect 1-st hidden layer with 2-nd hidden layer
!  b2  : biases of 2-nd hidden layer
!  w3  : weights connect 2-nd hidden layer with output layer
!  b3  : biase if output layer
!  vg  : output range of the training set, vg(2)
!
  implicit none
  integer,intent(in) :: nt,idv,n0,n1,n2
  real*8,intent(in) :: r(n0,nt)
  real*8,intent(out) :: v(nt),dv(n0,nt)
!
  real*8,intent(in) :: w1(n0,n1),b1(n1),w2(n1,n2),b2(n2),w3(n2),b3
  real*8,intent(in) :: rg(2,n0),vg(2)
!
  real*8 :: x0(n0,nt),x1(nt,n1),x2(nt,n2),tmp,rgg(n0),vgg
  integer :: i,j,k,n
!
  v=0.d0
  dv=0.d0
  do i=1,n0
    rgg(i)=rg(2,i)-rg(1,i)
  enddo
  vgg=vg(2)-vg(1)
!
  x0=r
!
  do i=1,n0
    x0(i,1:nt)=(x0(i,1:nt)-rg(1,i))/rgg(i)*2.d0-1.d0
  enddo
!
  do i=1,n1
   x1(1:nt,i)=b1(i)
  enddo
!
!  x1=trans(x0(n0,nt))*w1(n0,n1)+x1(nt,n1)
  call dgemm('t','n',nt,n1,n0,1.d0,x0,n0,w1,n0,1.d0,x1,nt)
!
  x1=dtanh(x1)
!
  do i=1,n2
   x2(1:nt,i)=b2(i)
  enddo
!
!  x2=x1(nt,n1)*w2(n1,n2)+x2(nt,n2)
  call dgemm('n','n',nt,n2,n1,1.d0,x1,nt,w2,n1,1.d0,x2,nt)
!
  x2=dtanh(x2)
!
  v(1:nt)=b3
!
!  v=x2(nt,n2)*w3(n2)+v(nt)
  call dgemv('n',nt,n2,1.d0,x2,nt,w3,1,1.d0,v,1)
!
  v=(v+1.d0)/2.d0*vgg+vg(1)
!
  if (idv.eq.0) return
!  calculate first derivatives
  do n=1,nt
    do i=1,n0
      do k=1,n2
        tmp=0.d0
        do j=1,n1
          tmp=tmp+w2(j,k)*w1(i,j)*(1-x1(n,j)**2)
        enddo
        dv(i,n)=dv(i,n)+w3(k)*tmp*(1-x2(n,k)**2)
      enddo
      dv(i,n)=dv(i,n)*vgg/(rgg(i))
    enddo
  enddo
  return
  end subroutine nsimx
! --------------------------------------------------------------------
  subroutine cart2dist3(c,r)
  implicit none
  integer,parameter :: natom=3
  real*8,intent(in) :: c(3,natom)
  real*8,intent(out) :: r(natom*(natom-1)/2)
  integer i,j,k,id
  id=0
  do i=1,natom-1
  do j=i+1,natom
    id=id+1
    r(id)=0.d0
    do k=1,3
    r(id)=r(id)+(c(k,i)-c(k,j))**2
    enddo
    r(id)=dsqrt(r(id))
  enddo
  enddo
  return
  end subroutine cart2dist3
! --------------------------------------------------------------------
  subroutine co2_nnfit_36(ndim,ntot,r0,vx)                       !@@@@
  implicit none
  integer,intent(in) :: ndim,ntot
  real*8,intent(in) :: r0(ndim,ntot)
  real*8,intent(out) :: vx(ntot)

  integer,parameter :: s0=3,s1=15,s2=15                          !@@@@
  character*99,save :: wfile=pesdir//'pes_co2_w36.dat' !@@@@

  real*8 :: dv(ndim,ntot),dvx(ndim),tmp
  integer :: i
  integer,save :: init=0
  integer,save :: fid
  real*8,save :: w1(s0,s1),b1(s1),w2(s1,s2),b2(s2),w3(s2),b3
  real*8,save :: rg(2,s0),vg(2)
  if(ndim.ne.s0) stop "ndim .ne. s0"
  if (init.eq.0) then
    fid=123456
    open(fid,file=trim(wfile),status='old',action='read')
    read(fid,*) w1,b1,w2,b2,w3,b3,rg,vg
    close(fid)
    init=1
  endif
  if(ntot.ge.24) then
    call nsimx(ntot,r0,vx,0,dv,s0,s1,s2,rg,w1,b1,w2,b2,w3,b3,vg)
  else
    do i=1,ntot
    call nsim(r0(1,i),vx(i),0,dvx,s0,s1,s2,rg,w1,b1,w2,b2,w3,b3,vg)
    enddo
  endif
  return
  end subroutine
! --------------------------------------------------------------------
end subroutine pot_co2
! --------------------------------------------------------------------
