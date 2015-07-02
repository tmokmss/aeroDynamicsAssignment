! Shock Tube Problem by solving 1D Euler eqs.
! for parfect gas

program page38
  implicit none
  integer, parameter :: num = 300, dime = 3
  real :: dx, cfl, dt, time, eps, g, rhol, ul, pl, rhor, &
          ur, pr, rgas, cvgas, ecp
  real, dimension(-5:num+5) :: x, rho, u, p, e
  real, dimension(dime, -5:num+5) :: q, qold, flux
  real, dimension(:,:), allocatable :: xt, rhot, ut, pt, et
  integer :: i, n, nlast, mx
  real :: atmpa, platm, pratm, tl, tr, xmax, xmin, tmsec, ustar

  ! set grid
  mx = num
  xmin = -2.0
  xmax = 2.0
  dx = (xmax-xmin)/float(mx-1)
  do i=-2,mx+3
    x(i) = xmin+dx*float(i-1)
  enddo

  ! set parameters
  atmpa = 1.013e+05
  eps = 1.0e-06
  g = 1.4
  rgas = 8.314e+03/28.96
  cvgas = rgas/(g-1.0)
  ul = 0
  tl = 400
  platm = 5
  pl = platm*atmpa
  rhol = pl/(rgas*tl)
  ur = 0
  tr = 300
  pratm = 1
  pr = pratm*atmpa
  rhor = pr/(rgas*tr)
  tmsec = 4
  time = tmsec*1.0e-03
  cfl = 0.5
  ecp = 0.125
  ustar = sqrt(g*rgas*tl)

  ! numerical solutions
  dt = cfl*dx/sqrt(g*rgas*amax1(tl,tr))
  nlast = int(time/dt)
  time = dt*float(nlast)

  !initial condition 
  do i=1,mx 
    if (x(i) < 0.0) then
      rho(i) = rhol
      u(i) = ul
      p(i) = pl
      e(i) = p(i)/((g-1.0)*rho(i))
    else
      rho(i) = rhor
      u(i) = ur
      p(i) = pr
      e(i) = p(i)/((g-1.0)*rho(i))
    endif
    q(1,i) = rho(i)
    q(2,i) = rho(i)*u(i)
    q(3,i) = rho(i)*(e(i)+0.5*u(i)**2)
  enddo
 
  allocate(rhot(1:nlast, -5:num+5), &
           ut(1:nlast, -5:num+5), &
           pt(1:nlast, -5:num+5), &
           et(1:nlast, -5:num+5))

  do n=1,nlast
    qold(1:3, 1:mx) = q(1:3, 1:mx)
    call calflx
    do i=4,mx-3
      q(1,i) = qold(1,i)-(dt/dx)*(flux(1,i)-flux(1,i-1))
      q(2,i) = qold(2,i)-(dt/dx)*(flux(2,i)-flux(2,i-1))
      q(3,i) = qold(3,i)-(dt/dx)*(flux(3,i)-flux(3,i-1))
      rho(i) = q(1,i)
      u(i) = q(2,i)/rho(i)
      p(i) = (g-1.0)*(q(3,i)-0.5*rho(i)*u(i)**2)
      e(i) = p(i)/((g-1.0)*rho(i))
    enddo
    rhot(n, :) = rho(:)
    ut(n, :) = u(:)
    pt(n, :) = p(:)
    et(n, :) = e(:)
  enddo

  ! write results
  open(unit=10,file='config38.txt', action='readwrite')
  write(10,*) 'mx: ', mx
  write(10,*) 'cfl: ', cfl
  write(10,*) 'dx: ', dx
  write(10,*) 'dt: ', dt
  write(10,*) 'nlast: ', nlast
  write(10,*) 'time: ', time
  write(10,*) 'g: ', g
  write(10,*) 'rhol: ', rhol
  write(10,*) 'ul: ', ul
  write(10,*) 'pl: ', pl
  write(10,*) 'rhor: ', rhor
  write(10,*) 'ur: ', ur
  write(10,*) 'pr: ', pr
  close(unit=10)
  
  open(unit=11,file='result38.csv', action='readwrite')
  write(11,*) 'x(m), density, u, p, T'
  do i=0, mx
    write(11, *) x(i),',',rho(i),',',u(i),',',p(i)/atmpa,',',e(i)/cvgas
  enddo
  close(unit=11)

  open(unit=12,file='ut38.csv', action='readwrite')
  do n=1, nlast
    write(12, fmt='(F15.7)',advance='no') n*dt
    write(12, fmt='(a)',advance='no') ','
    do i=0, mx-1
      write(12, fmt='(F15.7)',advance='no') ut(n,i)
      write(12, fmt='(a)', advance='no') ','
    enddo
    write(12, *) ut(n,mx)
  enddo
  close(unit=12)
  
  open(unit=13,file='rhot38.csv', action='readwrite')
  do n=1, nlast
    write(13, fmt='(F15.7)',advance='no') n*dt
    write(13, fmt='(a)',advance='no') ','
    do i=0, mx-1
      write(13, fmt='(F15.7)',advance='no') rhot(n,i)
      write(13, fmt='(a)', advance='no') ','
    enddo
    write(13, *) rhot(n,mx)
  enddo
  close(unit=13)

  open(unit=14,file='pt38.csv', action='readwrite')
  do n=1, nlast
    write(14, fmt='(F15.7)',advance='no') n*dt
    write(14, fmt='(a)',advance='no') ','
    do i=0, mx-1
      write(14, fmt='(F15.7)',advance='no') pt(n,i)
      write(14, fmt='(a)', advance='no') ','
    enddo
    write(14, *) pt(n,mx)
  enddo
  close(unit=14)

  open(unit=15,file='et38.csv', action='readwrite')
  do n=1, nlast
    write(15, fmt='(F15.7)',advance='no') n*dt
    write(15, fmt='(a)',advance='no') ','
    do i=0, mx-1
      write(15, fmt='(F15.7)',advance='no') et(n,i)
      write(15, fmt='(a)', advance='no') ','
    enddo
    write(15, *) et(n,mx)
  enddo
  close(unit=15)
  !call saveCsvt(et, 'et382.csv', 16)

contains
  subroutine saveCsvt(array, fname, unitnum)
    real, dimension(:, :), intent(in) :: array
    character(*), intent(in) :: fname
    integer, intent(in) :: unitnum

    open(unit=unitnum,file=fname, action='readwrite')
    do n=1, nlast
      write(unitnum, fmt='(F15.7)',advance='no') n*dt
      write(unitnum, fmt='(a)',advance='no') ','
      do i=0, mx-1
        write(unitnum, fmt='(F15.7)',advance='no') array(n,i)
        write(unitnum, fmt='(a)', advance='no') ','
      enddo
      write(unitnum, *) array(n,mx)
    enddo
    close(unit=unitnum)
  end subroutine

  subroutine calflx()
    real, dimension(dime, -5:num+5) :: ram, alfa
    real, dimension(dime, dime, -5:num+5) :: veck
    real :: q1lt, q1rt, q2lt, q2rt, q3lt, q3rt, dq1, dq2, dq3, &
            rlt, ult, plt, hlt, rrt, urt, prt, hrt, ubar, hbar, abar, &
            b2, b1, ecpx, f1lt, f1rt, f2lt, f2rt, f3lt, f3rt, ff1, ff2, &
            ph1, ph2, ph3, qqq, sss, tvd1, tvd2, tvd3
    do i=1,mx-1
      q1lt = q(1,i)
      q1rt = q(1,i+1)
      q2lt = q(2,i)
      q2rt = q(2,i+1)
      q3lt = q(3,i)
      q3rt = q(3,i+1)
      dq1 = q(1,i+1)-q(1,i)
      dq2 = q(2,i+1)-q(2,i)
      dq3 = q(3,i+1)-q(3,i)
      rlt = q1lt
      ult = q2lt/rlt
      plt = (g-1.0)*(q3lt-0.5*rlt*ult**2)
      hlt = (q3lt+plt)/rlt
      rrt = q1rt
      urt = q2rt/q1rt
      prt = (g-1.0)*(q3rt-0.5*rrt*urt**2)
      hrt = (q3rt+prt)/rrt
     
      ubar = (sqrt(rlt)*ult+sqrt(rrt)*urt)/(sqrt(rlt)+sqrt(rrt))
      hbar = (sqrt(rlt)*hlt+sqrt(rrt)*hrt)/(sqrt(rlt)+sqrt(rrt))
      abar = sqrt(amax1((g-1.0)*(hbar-0.5*ubar**2), &
                        amin1(g*plt/rlt, g*prt/rrt)))
      ram(1,i) = ubar-abar
      veck(1,1,i) = 1.0
      veck(2,1,i) = ubar-abar
      veck(3,1,i) = hbar-ubar*abar
      ram(2,i) = ubar
      veck(1,2,i) = 1.0
      veck(2,2,i) = ubar
      veck(3,2,i) = 0.5*ubar**2
      ram(3,i) = ubar+abar
      veck(1,3,i) = 1.0
      veck(2,3,i) = ubar+abar
      veck(3,3,i) = hbar+ubar*abar
      b2 = (g-1.0)/abar**2
      b1 = b2*ubar**2/2.0
      alfa(1,i) = 0.5*(b1+ubar/abar)*dq1+0.5*(-b2*ubar-1.0/abar)*dq2+0.5*b2*dq3
      alfa(2,i) = (1.0-b1)*dq1+b2*ubar*dq2-b2*dq3
      alfa(3,i) = 0.5*(b1-ubar/abar)*dq1+0.5*(-b2*ubar+1.0/abar)*dq2+0.5*b2*dq3
    enddo
    
    ! numerical flux
    do i=3,mx-3
      ecpx = ecp*amax1(abs(ram(1,i)), abs(ram(2,i)), abs(ram(3,i)))
      ff1 = (dt/dx)*ram(1,i)**2
      ff2 = abs(ram(1,i))
      if (ff2<ecpx) ff2 = (ff2**2+ecpx**2)*0.5/ecpx
      sss = sign(1.0, alfa(1,i))
      qqq = sss*amax1(0.0, amin1(sss*2.0*alfa(1,i-1) &
                                ,sss*2.0*alfa(1,i) &
                                ,sss*2.0*alfa(1,i+1) &
                                ,sss*0.5*(alfa(1,i-1)+alfa(1,i+1))))
      ph1 = -ff1*qqq-ff2*(alfa(1,i)-qqq)
      
      ff1 = (dt/dx)*ram(2,i)**2
      ff2 = abs(ram(2,i))
      if (ff2<ecpx) ff2 = (ff2**2+ecpx**2)*0.5/ecpx
      sss = sign(1.0, alfa(2,i))
      qqq = sss*amax1(0.0, amin1(sss*2.0*alfa(2,i-1) &
                                ,sss*2.0*alfa(2,i) &
                                ,sss*2.0*alfa(2,i+1) &
                                ,sss*0.5*(alfa(2,i-1)+alfa(2,i+1))))
      ph2 = -ff1*qqq-ff2*(alfa(2,i)-qqq)

      ff1 = (dt/dx)*ram(3,i)**2 
      ff2 = abs(ram(3,i))
      if (ff2<ecpx) ff2 = (ff2**2+ecpx**2)*0.5/ecpx
      sss = sign(1.0, alfa(3,i))
      qqq = sss*amax1(0.0, amin1(sss*2.0*alfa(3,i-1) &
                                ,sss*2.0*alfa(3,i) &
                                ,sss*2.0*alfa(3,i+1) &
                                ,sss*0.5*(alfa(3,i-1)+alfa(3,i+1))))
      ph3 = -ff1*qqq-ff2*(alfa(3,i)-qqq)

      tvd1 = veck(1,1,i)*ph1+veck(1,2,i)*ph2+veck(1,3,i)*ph3
      tvd2 = veck(2,1,i)*ph1+veck(2,2,i)*ph2+veck(2,3,i)*ph3
      tvd3 = veck(3,1,i)*ph1+veck(3,2,i)*ph2+veck(3,3,i)*ph3
      q1lt = q(1,i)
      q1rt = q(1,i+1)
      q2lt = q(2,i)
      q2rt = q(2,i+1)
      q3lt = q(3,i)
      q3rt = q(3,i+1)

      rlt = q1lt
      ult = q2lt/rlt
      plt = (g-1.0)*(q3lt-0.5*rlt*ult**2)
      hlt = (q3lt+plt)/rlt
      rrt = q1rt
      urt = q2rt/q1rt
      prt = (g-1.0)*(q3rt-0.5*rrt*urt**2)
      hrt = (q3rt+prt)/rrt
      f1lt = rlt*ult
      f2lt = rlt*ult**2+plt
      f3lt = ult*(q3lt+plt)
      f1rt = rrt*urt
      f2rt = rrt*urt**2+prt
      f3rt = urt*(q3rt+prt)
      flux(1,i) = 0.5*(f1lt+f1rt+tvd1)
      flux(2,i) = 0.5*(f2lt+f2rt+tvd2)
      flux(3,i) = 0.5*(f3lt+f3rt+tvd3)
    enddo
  end subroutine
end program
