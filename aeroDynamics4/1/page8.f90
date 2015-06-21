! 1D hyperbolic problem
! 2015/06/21

program page8
  implicit none
  integer, parameter :: size = 100
  real, dimension(0:size) ::  x, yi, ye, ys, work
  
  integer, parameter :: mx = size
  integer, parameter :: nlast = 4000000
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
      print *, i
    else
      ye(i) = 0.0
    endif
  enddo

  ! numerical solution
  ys = yi

  do n=0,nlast
    do j=i,mx-1
      ! scheme #1
      work(i) = ys(i)-0.5*cfl*(ys(i+1)-ys(i-1))
      ! scheme #2
      ! work(i) = ys(i)-cfl*(ys(i)-ys(i-1))
      ! scheme #3 Lax-Wendroff
      !work(i) = 0.5*cfl*(1.+cfl)*ys(i-1)+(1.-cfl**2)*ys(i) - 0.5*cfl*(1.-cfl)*ys(i+1)
    enddo
    !work = (/ () /)
    ys(1:mx-1) = work(1:mx-1)
  enddo

  print *, ys

end program
