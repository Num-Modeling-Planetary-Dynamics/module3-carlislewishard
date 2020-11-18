program xv2el
   use, intrinsic                :: iso_fortran_env
   implicit none
   integer, parameter            :: P = real64
   real(P)                       :: R_norm, V_norm, h_norm, Rdot
   real(P)                       :: sin_omega, cos_omega, sin_varpi, cos_varpi, sin_f, cos_f
   real(P), dimension(3)         :: R, V, h
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

   do j = 1,1001  
      !Copy t (time)
      output_data(j,1) = input_data(j,1)

      !Calculate R 
      R(1) = input_data(j,2)
      R(2) = input_data(j,3)
      R(3) = input_data(j,4)
      R_norm = norm2(R(:)) 

      !Calculate V 
      V(1) = input_data(j,5)
      V(2) = input_data(j,6)
      V(3) = input_data(j,7)
      V_norm = norm2(V(:)) !/ 365.25_P !AU/days

      !Calculate h
      h(1) = (R(2) * V(3)) - (R(3) * V(2))
      h(2) = (R(3) * V(1)) - (R(1) * V(3))
      h(3) = (R(1) * V(2)) - (R(2) * V(1))
      h_norm = norm2(h(:))

      !Calculate Rdot
      Rdot = sign(V_norm**2 - (h_norm / R_norm)**2, dot_product(R(:), V(:)))

      !Calculate a (semi-major axis)
      output_data(j,2) = ((2.0_P / R_norm) - V_norm**2)**(-1)

      !Calculate e (eccentricity)
      output_data(j,3) = sqrt(1.0_P - (h_norm**2 / output_data(j,2)))

      !Calculate inc (inclination)
      output_data(j,4) = acos(h(3) / h_norm)

      !Calculate omega (longitude of the ascending note)
      sin_omega = sign(h(1), h(3)) / (h_norm * sin(output_data(j,4)))
      cos_omega = sign(-h(2), h(3)) / (h_norm * sin(output_data(j,4)))
      output_data(j,5) = atan2(sin_omega, cos_omega)

      !Calculate f (true anomaly)
      sin_f = output_data(j,2) * (1.0_P - output_data(j,3)**2) / (h_norm * output_data(j,3)) * Rdot
      cos_f = ((output_data(j,2) * (1.0_P - output_data(j,3)**2) / R_norm) - 1.0_P) / output_data(j,3)
      output_data(j,7) = atan2(sin_f, cos_f)

      !Calculate varpi (argument of periapsis)
      sin_varpi = R(3) / (R_norm * sin(output_data(j,4)))
      cos_varpi = ((R(1) / R_norm) + (sin(output_data(j,5)) * sin_varpi * cos(output_data(j,4)))) / cos(output_data(j,5))
      output_data(j,6) = atan2(sin_varpi, cos_varpi) - output_data(j,5)

      !Write to the output file
      write(12,fmt) output_data(j,:)

   end do

   close(12)

   call system_clock(count_end, count_rate)
   write(*,*) "INTEGRATION TIME (S) : ", ((count_end  - count_start) / (count_rate * 1.0_P))

end program xv2el