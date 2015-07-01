! 1D hyperbolic problem
! 2015/06/21

program page8
  implicit none
  integer, parameter :: size = 10
  real, dimension(0:size) ::  x, yi, ye, ys, work
  
  integer, parameter :: mx = size
  integer, parameter :: nlast = 5
  real, parameter :: cfl = 0.5

  integer i, n, j
  real :: travel = real(nlast)*cfl
  integer :: centerx = mx / 2
  ! set grid
  do i=0,mx
    x(i) = real(i)
  end do
  
  ! set initial condition
  do i=0,mx
    if (x(i) < centerx) then
      yi(i) = 1.0
    else
      yi(i)=0.0
    endif
  enddo

  ! set exact solution
  do i=0, mx
    if (x(i) < centerx+travel) then
      ye(i) = 1.0
    else
      ye(i) = 0.0
    endif
  enddo

  ! numerical solution
  ys = yi

  do n=1,nlast
    do i=1,mx-1
      ! scheme #1
      !work(i) = ys(i)-0.5*cfl*(ys(i+1)-ys(i-1))
      ! scheme #2
      ! work(i) = ys(i)-cfl*(ys(i)-ys(i-1))
      ! scheme #3 Lax-Wendroff
      work(i) = 0.5*cfl*(1.+cfl)*ys(i-1)+(1.-cfl**2)*ys(i) - 0.5*cfl*(1.-cfl)*ys(i+1)
    enddo
    ys(1:mx-1) = work(1:mx-1)
  enddo

  ! write results
  open(unit=10,file='config.txt', action='readwrite')
  write(10,*) 'mx: ', mx
  write(10,*) 'cfl: ', cfl
  write(10,*) 'nlast: ', nlast
  close(unit=10)
  
  open(unit=11,file='result.csv', action='readwrite')
  write(11,*) 'j, initial, exact, numerical'
  do i=0, mx
    write(11, *) x(i),',',yi(i),',',ye(i),',',ys(i)
  enddo
  close(unit=11)
  print *, ys

end program
