program xv2el
   use, intrinsic                :: iso_fortran_env
   implicit none
   integer, parameter            :: P = real64
   real(P)                       :: R_norm, V_norm, h_norm, Rdot, PI, G, msun, mu
   real(P)                       :: a, e, inc, omega, varpi, f
   real(P)                       :: sin_omega, cos_omega, sin_varpi, cos_varpi, sin_f, cos_f
   real(P), dimension(3)         :: R, V, h
   real(P), dimension(1001)      :: t, px, py, pz, vx, vy, vz
   real(P), dimension(1001,7)    :: input_data, output_data
   integer                       :: count_start, count_rate, count_end, i, j
   character(len = 200)          :: file_in, file_out, path
   character(len=*), parameter   :: fmt = '(ES23.16,10(",",ES23.16,:))'
   call system_clock(count_start, count_rate)
   
   path = '/Users/carlislewishard/Documents/Classes/Year_4/Numerical_Dynamics/Module3/&
   module3-carlislewishard/data/'

   file_in = trim(path) // trim('id000004-XV.csv')
   file_out = trim(path) // trim('id000004-EL.csv')
   open(unit = 11, file = file_in, status = 'old')
   open(unit = 12, file = file_out, status = 'replace')
   read(11,*) !Skip the first line
   do i = 1,1001
      read(11,*) input_data(i,:)
   end do

   close(11)

   !Write the output header
   write(12,*) "t,a,e,inc,omega,varpi,f"

   !Calculate the constants and define the output units
   PI = 4.0_P * atan(1.0_P)
   G = 4.0_P * PI**2 
   msun = 1.0_P
   mu = G * msun

   !Name the input data
   t(:) = input_data(:,1)
   px(:) = input_data(:,2)
   py(:) = input_data(:,3)
   pz(:) = input_data(:,4)
   vx(:) = input_data(:,5)
   vy(:) = input_data(:,6)
   vz(:) = input_data(:,7)

   do j = 1,1001  

      !Calculate R 
      R(1) = px(j) !AU
      R(2) = py(j) !AU
      R(3) = pz(j) !AU
      R_norm = norm2(R(:)) 

      !Calculate V 
      V(1) = vx(j) * 365.25_P !AU/day -> AU/year
      V(2) = vy(j) * 365.25_P !AU/day -> AU/year
      V(3) = vz(j) * 365.25_P !AU/day -> AU/year
      V_norm = norm2(V(:))

      !Calculate h
      h(1) = (R(2) * V(3)) - (R(3) * V(2))
      h(2) = (R(3) * V(1)) - (R(1) * V(3))
      h(3) = (R(1) * V(2)) - (R(2) * V(1))
      h_norm = norm2(h(:))

      !Calculate Rdot
      Rdot = sign(V_norm**2 - (h_norm / R_norm)**2, dot_product(R(:), V(:)))

      !Calculate a (semi-major axis)
      a = ((2.0_P / R_norm) - (V_norm**2 / mu))**(-1)

      !Calculate e (eccentricity)
      e = sqrt(1.0_P - (h_norm**2 / (a * mu)))

      !Calculate inc (inclination)
      inc = acos(h(3) / h_norm)

      !Calculate omega (longitude of the ascending note)
      sin_omega = sign(h(1), h(3)) / (h_norm * sin(inc))
      cos_omega = sign(-h(2), h(3)) / (h_norm * sin(inc))
      omega = atan2(sin_omega, cos_omega)

      !Calculate f (true anomaly)
      sin_f = a * (1.0_P - e**2) / (h_norm * e) * Rdot
      cos_f = ((a * (1.0_P - e**2) / R_norm) - 1.0_P) / e
      f = atan2(sin_f, cos_f)

      !Calculate varpi (argument of periapsis)
      sin_varpi = R(3) / (R_norm * sin(inc))
      cos_varpi = ((R(1) / R_norm) + (sin(inc) * sin_varpi * cos(inc))) / cos(omega)
      varpi = atan2(sin_varpi, cos_varpi) - omega

      !Name the output data
      output_data(j,1) = t(j)
      output_data(j,2) = a
      output_data(j,3) = e
      output_data(j,4) = inc
      output_data(j,5) = omega
      output_data(j,6) = varpi
      output_data(j,7) = f

      !Write to the output file
      write(12,fmt) output_data(j,:)

   end do

   close(12)

   call system_clock(count_end, count_rate)
   write(*,*) "INTEGRATION TIME (S) : ", ((count_end  - count_start) / (count_rate * 1.0_P))

end program xv2el