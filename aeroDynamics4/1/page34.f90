! 1D linear-hyperbolic problem
! Check Symmetric-TVD scheme

program page34
  implicit none
  integer, parameter :: size = 100
  real, dimension(-1:size) :: x, yi, ye, ys1, ys2, ys3, flux, work
  real :: c, cfl, dx, dt, travel
  integer :: mx, nlast

  c = 1.0
  mx = 51
  nlast = 30
  cfl = 0.5
  dx = 0.02
  dt = cfl*dx
  travel = dt*c*float(nlast)

  ! set grid
  do i=-1,mx+1
    x(i) = dx*float(i-1)
  enddo

  ! set initial condition
  do i = -1,mx+1
    if(x(i) < 0.5) then
      yi(i) = 1.0
    else
      yi(i) = 0.0
    endif
  enddo

  ! set exact condition
  do i=-1,mx+1
    if (x(i)<0.5+travel) then
      ye(i) = 1.0
    else
      ye(i) = 0.0
    endif
  enddo

  ! solve scheme by minmod1
  do i=-1, mx+1
    ys1(i) = yi(i)
  enddo

  do n=1, nlast
    do i=0, mx-1
      dm = ys1(i)-ys1(i-1)
      d0 = ys1(i+1)-ys1(i)
      dp = ys1(i+2)-ys1(i+1)
      s = sign(1.0, dm)
      q = s*amax1(0.0, amin1(s*dm,s*d0,s*dp))
      ac = abs(c)
      if (ac < ecp) ac = (c**2+ecp**2)*0.5/ecp
      f = -(dt*c**2/dx*q+ac*(d0-q))
      flux(i) = 0.5*(c*ys1(i)+c*ys1(i+1)+f)
    enddo
    do i=1, mx-1
      work(i) = ys1(i)-(dt/dx)*(flux(i)-flux(i-1))
    enddo
    ys1(1:mx-1) = work(1:mx-1)
  enddo

end program
