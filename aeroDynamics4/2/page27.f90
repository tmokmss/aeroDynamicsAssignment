! Solve Riemann Problem for Linear Hyperbolic System

program page27
  implicit none
  integer, parameter :: size = 300
  integer, parameter :: dim = 2
  integer :: mx
  real :: dx, xmin, xmax
  real, dimension(-5:size+5) :: x
  real, dimension(dim) :: ram
  real, dimension(dim, dim) :: a, r, ri
  real, dimension(-5:size+5) :: us, vs, ui, vi
  real, dimension(dim, -5:size+5) :: alfa, q, qold, flux
  real, dimension(-5:size+5) :: u, v

  integer :: i, n, nlast
  real :: bbb, ccc, ddd, det
  real :: ul, vl, ur, vr, cfl, dt, a1l, a2l, a1r, a2r
  real :: dq1, dq2, alfa1, alfa2, f1lt, f2lt, f1rt, f2rt

  ! set grid
  mx = 201
  xmin = -1.0
  xmax = 1.0
  dx = (xmax-xmin)/float(mx-1)
  do i=-2,mx+3
    x(i) = xmin+dx*float(i-1)
  enddo

  ! set parameters
  a(1, 1) = 5
  a(1, 2) = 2
  a(2, 1) = 10
  a(2, 2) = 4

  !eigenvalues & vectors
  bbb = -(a(1,1)+a(2,2))
  ccc = a(1,1)*a(2,2)-a(1,2)*a(2,1)
  ddd = bbb**2 - 4.0*ccc
  ram(1) = (-bbb-sqrt(ddd))/2.0
  ram(2) = (-bbb+sqrt(ddd))/2.0
  r(1,1) = 1.0
  r(1,2) = 1.0
  r(2,1) = -(a(1,1)-ram(1))/a(1,2)
  r(2,2) = -(a(1,1)-ram(2))/a(1,2)
  det = r(1,1)*r(2,2)-r(1,2)*r(2,1)
  ri(1,1) = r(2,2)/det
  ri(2,1) = -r(2,1)/det
  ri(1,2) = -r(1,2)/det
  ri(2,2) = r(1,1)/det
  write(*,*) r(1,1), r(1,2)
  write(*,*) r(2,1), r(2,2)
  write(*,*) ram
  write(*,*) ri(1,1), ri(1,2)
  write(*,*) ri(2,1), ri(2,2)
  
  ul = 1.0
  vl = 5.0
  ur = 3.0
  vr = 2.0
  cfl = 0.5
  nlast = 40
  dt = cfl*dx/amax1(ram(1),ram(2))

  ! solve exact solution
  a1l = ri(1,1)*ul+ri(1,2)*vl
  a2l = ri(2,1)*ul+ri(2,2)*vl
  a1r = ri(1,1)*ur+ri(1,2)*vr
  a2r = ri(2,1)*ur+ri(2,2)*vr
  do i=1, mx
    alfa(1,i) = a1l
    if (x(i)>ram(1)*dt*float(nlast)) alfa(1,i) = a1r
    alfa(2,i) = a2l
    if (x(i)>ram(2)*dt*float(nlast)) alfa(2,i) = a2r
    us(i) = r(1,1)*alfa(1,i)+r(1,2)*alfa(2,i)
    vs(i) = r(2,1)*alfa(1,i)+r(2,2)*alfa(2,i)
  enddo

  ! numerical solutions
  ! initial condition
  do i=1, mx
    if (x(i)<=0.0) then
      u(i) = ul
      v(i) = vl
    else
      u(i) = ur
      v(i) = vr
    endif
    q(1,i) = u(i)
    q(2,i) = v(i)
    ui(i) = u(i)
    vi(i) = v(i)
  enddo

  ! Roe method
  do n=1, nlast
    do i=1, mx-1
      dq1 = q(1,i+1)-q(1,i)
      dq2 = q(2,i+1)-q(2,i)
      alfa1 = ri(1,1)*dq1+ri(1,2)*dq2
      alfa2 = ri(2,1)*dq1+ri(2,2)*dq2
      f1lt = a(1,1)*u(i)+a(1,2)*v(i)
      f2lt = a(2,1)*u(i)+a(2,2)*v(i)
      f1rt = a(1,1)*u(i+1)+a(1,2)*v(i+1)
      f2rt = a(2,1)*u(i+1)+a(2,2)*v(i+1)
      flux(1,i) = 0.5*(f1lt+f1rt) &
                - 0.5*alfa1*abs(ram(1))*r(1,1) &
                - 0.5*alfa2*abs(ram(2))*r(1,2)
      flux(2,i) = 0.5*(f2lt+f2rt) &
                - 0.5*alfa1*abs(ram(1))*r(2,1) &
                - 0.5*alfa2*abs(ram(2))*r(2,2)
    enddo
    do i=2, mx-1
      q(1,i) = q(1,i)-(dt/dx)*(flux(1,i)-flux(1,i-1))
      q(2,i) = q(2,i)-(dt/dx)*(flux(2,i)-flux(2,i-1))
      u(i) = q(1,i)
      v(i) = q(2,i)
    enddo
  enddo
 
  ! write results
  open(unit=10,file='config27.txt', action='readwrite')
  write(10,*) 'mx: ', mx
  write(10,*) 'cfl: ', cfl
  write(10,*) 'dx: ', dx
  write(10,*) 'dt: ', dt
  write(10,*) 'nlast: ', nlast
  write(10,*) 'ul: ', ul
  write(10,*) 'vl: ', vl
  write(10,*) 'ur: ', ur
  write(10,*) 'vr: ', vr
  close(unit=10)
  
  open(unit=11,file='result27.csv', action='readwrite')
  write(11,*) 'x, u(initial), v, u(Exact), v, u(Comp), v'
  do i=1, mx
    write(11, *) x(i),',',ui(i),',',vi(i),',',us(i),',',vs(i),',',u(i),',',v(i)
  enddo
  close(unit=11)

end program
