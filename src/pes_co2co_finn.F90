#define pesdir "./"
!program test
!  implicit none
!  real*8 :: r(6),v,dvdr(6),vco,dvco
!  !! r: bohr, O1-C1, O1-O2, O1-C2, C1-O2, C1-C2, O2-C2
!  !! v: cm-1;  dvdr: cm-1/bohr
!  r=[2.709248d0,10.463631d0,13.177895d0,8.340613d0,11.211281d0,3.161466d0]
!  call pes_ococ(r,v,1,dvdr)
!  print*,v
!  print*,dvdr
!  call potco(r(1),vco,dvco)
!  print*,vco,dvco
!end program

subroutine pes_co2co(x,v,vco2,vco,vint,dvdr,dvco2,dvco)
  implicit none
  real*8,intent(in) :: x(3,5) !! Cartesian of O O C O C, Angstrom
  real*8,intent(out) :: v,vco2,vco,vint,dvdr(10),dvco,dvco2(3)
  real*8 :: x_o2c_oc(3,5),r(10),rs(10,21),ps(13,21),rco,p(13),g(10)
  real*8 :: r1(10),r2(10),fim(13)
  real*8 :: vfita(21)
  real*8,parameter :: cm=219474.5d0
  real*8,parameter :: ang=0.5291772083d0
  real*8,parameter :: ev=27.21138505d0
  real*8,parameter :: dr=0.01d0
  integer :: ipes
  integer :: i,j,k,l
  real*8 :: va
  
  !!! input coordinate order : OOC-OC 
  !x_o2c_oc=x(:,[2,3,1,5,4]) ! OOC-OC Angstrom

  call cart2dist5(x,r)

  k=1
  r2=r
  do l=1,10,1
     g(l)=1.0d0/r2(l)
  end do
  call fi_a2bcd(g,p) ! r: O2C + OC
  call pipnn_vdw_o2c_oc_fit01_a(13,1,p,vfita(k))
  do i=1,10,1
     do j=-1,1,2
        k=k+1
        r1(1:10)=r(1:10)
        r1(i)=r1(i)+j*dr
        r2(1:10)=r1(1:10)
        do l=1,10,1
           g(l)=1.0d0/r2(l)
        end do
        call fi_a2bcd(g,p) ! r: O2C + OC
        call pipnn_vdw_o2c_oc_fit01_a(13,1,p,vfita(k))
     end do
  end do

  vint=vfita(1)  ! Vint cm-1

  dvdr(1)=(vfita(3)-vfita(2))/2.d0/dr ! cm-1/Angstrom
  dvdr(2)=(vfita(5)-vfita(4))/2.d0/dr ! cm-1/Angstrom
  dvdr(3)=(vfita(7)-vfita(6))/2.d0/dr ! cm-1/Angstrom
  dvdr(4)=(vfita(9)-vfita(8))/2.d0/dr ! cm-1/Angstrom
  dvdr(5)=(vfita(11)-vfita(10))/2.d0/dr ! cm-1/Angstrom
  dvdr(6)=(vfita(13)-vfita(12))/2.d0/dr ! cm-1/Angstrom
  dvdr(7)=(vfita(15)-vfita(14))/2.d0/dr ! cm-1/Angstrom
  dvdr(8)=(vfita(17)-vfita(16))/2.d0/dr ! cm-1/Angstrom
  dvdr(9)=(vfita(19)-vfita(18))/2.d0/dr ! cm-1/Angstrom
  dvdr(10)=(vfita(21)-vfita(20))/2.d0/dr ! cm-1/Angstrom
 
  call pot_co2_cart1(x(:,1:3)/ang,vco2,dvco2)  ! bohr; a.u

  rco=dsqrt(dot_product(x(:,4)-x(:,5),x(:,4)-x(:,5)))/ang ! bohr
  call potco(rco,vco,dvco) ! cm-1, cm-1/bohr

  v=vint+vco2*cm+vco ! cm-1

  return
end subroutine pes_co2co

subroutine fi_a2bcd(r,p)
  ! compute the fundamental invariants
  ! atom order : A A B C D
  ! r (input) real(kind=8) the internuclear distances
  ! p (output) real(kind=8) the fundamental invariants
  implicit none
  real(kind=8),intent(in) :: r(10)
  real(kind=8),intent(out) :: p(13)
  p(1)=r(2)+r(5)
  p(2)=r(3)+r(6)
  p(3)=r(4)+r(7)
  p(4)=dsqrt(r(2)**2+r(5)**2)
  p(5)=dsqrt(r(3)**2+r(6)**2)
  p(6)=dsqrt(r(4)**2+r(7)**2)
  p(7)=dsqrt(r(2)*r(3)+r(5)*r(6))
  p(8)=dsqrt(r(2)*r(4)+r(5)*r(7))
  p(9)=dsqrt(r(3)*r(4)+r(6)*r(7))
  p(10)=r(1)
  p(11)=r(8)
  p(12)=r(9)
  p(13)=r(10)
end subroutine fi_a2bcd

subroutine pipnn_vdw_o2c_oc_fit01_a(ndim,ntot,r0,vx)
  implicit none
  integer,intent(in) :: ndim,ntot
  real*8,intent(in) :: r0(ndim,ntot)
  real*8,intent(out) :: vx(ntot)

  integer,parameter :: s0=13,s1=45,s2=45
  character*99,save :: wfile=pesdir//"pes_co2co_w05.dat"
  real*8,save :: w1(s0,s1),b1(s1),w2(s1,s2),b2(s2),w3(s2),b3,rg(2,s0),vg(2)

  real*8 :: dv(ndim,ntot)
  integer :: i
  integer,save :: init=0
  integer,save :: fid
  if(ndim.ne.s0) stop "ndim .ne. s0"
  if (init.eq.0) then
    init=1
    fid=123456
    open(fid,file=trim(wfile),status='old',action='read')
    read(fid,*) w1,b1,w2,b2,w3,b3,rg,vg
    close(fid)
  endif
  if(ntot.ge.24) then
    call nsimx(ntot,r0,vx,0,dv,s0,s1,s2,rg,w1,b1,w2,b2,w3,b3,vg)
  else
    do i=1,ntot
    call nsim(r0(1,i),vx(i),0,dv(:,i),s0,s1,s2,rg,w1,b1,w2,b2,w3,b3,vg)
    enddo
  endif
  return
end subroutine
! --------------------------------------------------------------------
subroutine potco(r,v,dv)
  ! 2020/4/21  C-O AE-CCSD(T)(short range) MRCI(long range) / ACV5Z Energy in cm-1
  ! updated: 2020/6/25, dv/dr
  implicit none
  real*8,intent(in) :: r !! bohr
  real*8,intent(out) :: v,dv !! cm-1, cm-1/bohr
  real*8 :: x,y1j
  integer j
  integer,parameter :: ns=8
  real*8,parameter :: w1(ns)=[  2.5731334302226245d+00,&
                               -4.1635475296365723d+00,&
                                6.2652688635816842d+00,&
                               -8.1228824615517947d+00,&
                                2.5634824563114837d+01,&
                               -1.6643878666848999d+00,&
                                1.2863481593265147d+01,&
                               -5.1395051186435685d+00]
  real*8,parameter :: b1(ns)=[ -1.6421783363943476d+00,&
                                1.5364774081764951d+00,&
                               -1.1338455332512607d+00,&
                               -1.4448051810696709d+00,&
                                7.5217573991947644d+00,&
                               -1.4005229119446290d+00,&
                                1.1053854378210930d+01,&
                               -5.9299626180269485d+00]
  real*8,parameter :: w2(ns)=[ -8.8664820626030757d-03,&
                                7.8571245473773067d-03,&
                               -3.6411047563733342d-03,&
                               -4.0358215533209145d-03,&
                                9.6640587889626451d-04,&
                               -1.4325782866595651d+00,&
                                1.2002568907875554d-02,&
                                8.3983298757280007d+00]
  real*8,parameter :: b2=6.8970075393140338d+00
  real*8,parameter :: ra=1.4000000000000000d+00,rb=7.0000000000000000d+00
  real*8,parameter :: va=9.8087308941410573d-02,vb=1.9558422718340193d+05

  v=0.d0
  dv=0.d0
  x=2.d0*(r-ra)/(rb-ra)-1.d0
  do j=1,ns
    y1j=tanh(b1(j)+w1(j)*x)
    v=v+w2(j)*y1j
    dv=dv+w2(j)*(1.d0-y1j**2)*w1(j)
  enddo
  v=v+b2
  v=(v+1.d0)*(vb-va)/2.d0+va
  v=v+0.560096850315234d0
  dv=dv*(vb-va)/(rb-ra)
  return
end subroutine potco
! --------------------------------------------------------------------
subroutine nsim(r0,v,idv,dv,n0,n1,n2,rg,w1,b1,w2,b2,w3,b3,vg)
  ! ---- simulate a neural network with two hidden layers
  ! ----      n0-n1-n2-1
  ! ---- blas routines used in this subroutine: dgemv, ddot
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
subroutine nsimx(nt,r,v,idv,dv,n0,n1,n2,rg,w1,b1,w2,b2,w3,b3,vg)
  ! ---- simulate a neural network with two hidden layers
  ! ----    n0-n1(tansig)-n2(tansig)-1(pureline)
  ! ---- blas routines used in this subroutine: dgemm, dgemv
  ! ---- the syntax keeps the same with nsim, except input nt at first
  ! ---- Chen Jun, 2013/11/26
  !
  ! ---- description of parameters
  ! ---- nt  : total number of points in r,v and dv
  ! ---- r   : input of nn, r(n0,nt), r(1:n0,i) is the i-th input
  ! ---- v   : output of nn, v(nt), v(i) is the i-th input
  ! ---- idv : if idv==0, calc. v only, else calc. v and derivatives dv
  ! ---- dv  : derivatives of dv/dr, dv(n0,nt), dv(j,i)=d v(i) / d r(j,i)
  ! ---- n0  : neurons of the input layer, dimension of r
  ! ---- n1  : neurons of the first hidden layer
  ! ---- n2  : neurons of the second hidden layer
  ! ---- rg  : input ranges of the training set, rg(2,n0)
  ! ---- w1  : weights connect input layer with 1-st hidden layer
  ! ---- b1  : biases of 1-st hidden layer
  ! ---- w2  : weights connect 1-st hidden layer with 2-nd hidden layer
  ! ---- b2  : biases of 2-nd hidden layer
  ! ---- w3  : weights connect 2-nd hidden layer with output layer
  ! ---- b3  : biase of output layer
  ! ---- vg  : output range of the training set, vg(2)
  implicit none
  integer,intent(in) :: nt,idv,n0,n1,n2
  real*8,intent(in) :: r(n0,nt)
  real*8,intent(out) :: v(nt),dv(n0,nt)
  real*8,intent(in) :: w1(n0,n1),b1(n1),w2(n1,n2),b2(n2),w3(n2),b3
  real*8,intent(in) :: rg(2,n0),vg(2)
  real*8 :: x0(n0,nt),x1(nt,n1),x2(nt,n2),tmp,rgg(n0),vgg
  integer :: i,j,k,n

  v=0.d0
  do i=1,n0
    rgg(i)=rg(2,i)-rg(1,i)
  enddo
  vgg=vg(2)-vg(1)

  x0=r

  do i=1,n0
    x0(i,1:nt)=(x0(i,1:nt)-rg(1,i))/rgg(i)*2.d0-1.d0
  enddo

  do i=1,n1
   x1(1:nt,i)=b1(i)
  enddo

  ! ---- x1=trans(x0(n0,nt))*w1(n0,n1)+x1(nt,n1)
  !call gpu_dgemm('t','n',nt,n1,n0,1.d0,x0,n0,w1,n0,1.d0,x1,nt)
  call dgemm('t','n',nt,n1,n0,1.d0,x0,n0,w1,n0,1.d0,x1,nt)

  x1=dtanh(x1)

  do i=1,n2
   x2(1:nt,i)=b2(i)
  enddo

  ! ---- x2=x1(nt,n1)*w2(n1,n2)+x2(nt,n2)
  !call gpu_dgemm('n','n',nt,n2,n1,1.d0,x1,nt,w2,n1,1.d0,x2,nt)
  call dgemm('n','n',nt,n2,n1,1.d0,x1,nt,w2,n1,1.d0,x2,nt)

  x2=dtanh(x2)

  v(1:nt)=b3

  ! ---- v=x2(nt,n2)*w3(n2)+v(nt)
  call dgemv('n',nt,n2,1.d0,x2,nt,w3,1,1.d0,v,1)

  v=(v+1.d0)/2.d0*vgg+vg(1)

  if (idv.ne.1) return
  ! ---- calculate first derivatives
  dv=0.d0
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
subroutine cart2dist5(c,r)
  implicit none
  integer,parameter :: natom=5
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
end subroutine

