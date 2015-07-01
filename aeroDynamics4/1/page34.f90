! 1D linear-hyperbolic problem
! Check Symmetric-TVD scheme

program page34
  implicit none
  integer, parameter :: size = 520
  real, dimension(-1:size) :: x, yi, ye, ys1, ys2, ys3, flux, work
  real :: c, cfl, dx, dt, travel
  integer :: mx, nlast

  real :: dm,d0,dp,s,q,ac,f,ecp,sb1,sb2,scp
  integer :: i,n

  c = 1.0
  mx = size-1
  nlast = 30
  cfl = 0.5
  dx = 0.002
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
  ys1(-1:mx+1) = yi(-1:mx+1)

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

  ! solve scheme by minmod2
  ys2(-1:mx+1) = yi(-1:mx+1)

  do n=1, nlast
    do i=0, mx-1
      dm = ys2(i)-ys2(i-1)
      d0 = ys2(i+1)-ys2(i)
      dp = ys2(i+2)-ys2(i+1)
      s = sign(1.0, dm)
      q = s*amax1(0.0, amin1(2.*s*dm, 2.*s*d0, 2.*s*dp, 0.5*s*(dm+dp)))
      ac =abs(c)
      if (ac < ecp) ac = (c**2+ecp**2)*0.5/ecp
      f = -(dt*c**2/dx*q+ac*(d0-q))
      flux(i) = 0.5*(c*ys2(i)+c*ys2(i+1)+f)
    enddo
    do i=1,mx-1
      work(i) = ys2(i)-cfl*(flux(i)-flux(i-1))
    enddo
    ys2(1:mx-1) = work(1:mx-1)
  enddo

  ! solve scheme by superbee
  ys3(-1:mx+1) = yi(-1:mx+1)
  do n=1, nlast
    do i=0,mx-1
      dm = ys3(i)-ys3(i-1)
      d0 = ys3(i+1)-ys3(i)
      dp = ys3(i+2)-ys3(i+1)
      s = sign(1.0, dm)
      sb1 = s*amax1(0.0, amin1(2.*abs(dm),s*d0), amin1(abs(dm), 2.*s*d0))
      s = sign(1.0, d0)
      sb2 = s*amax1(0.0, amin1(2.*abs(d0), s*dp), amin1(abs(d0), 2.*s*dp))
      q = sb1+sb2-d0
      ac = abs(c)
      if (ac<scp) ac = (c**2+ecp**2)*0.5/ecp
      f = -(dt*c**2/dx*q+ac*(d0-q))
      flux(i) = 0.5*(c*ys3(i)+c*ys3(i+1)+f)
    enddo
    do i=1, mx-1
      work(i) = ys3(i) - cfl*(flux(i)-flux(i-1))
    enddo
    ys3(1:mx-1) = work(1:mx-1)
  enddo
  
  ! write results
  open(unit=10,file='config34.txt', action='readwrite')
  write(10,*) 'mx: ', mx
  write(10,*) 'cfl: ', cfl
  write(10,*) 'nlast: ', nlast
  close(unit=10)
  
  open(unit=11,file='result34.csv', action='readwrite')
  write(11,*) 'j, initial, exact, minmod1, minmod2, superbee'
  do i=0, mx
    write(11, *) x(i),',',yi(i),',',ye(i),',',ys1(i),',',ys2(i),',',ys3(i)
  enddo
  close(unit=11)

end program
