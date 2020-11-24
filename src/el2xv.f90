program el2xv
   use, intrinsic                :: iso_fortran_env
   implicit none
   integer, parameter            :: P = real64
   real(P)                       :: R_norm, V_norm, h_norm, PI, G, msun, mu
   real(P)                       :: px, py, pz, vx, vy, vz, rdot
   real(P), dimension(3)         :: R, V, h
   real(P), dimension(1001)      :: t, a, e, inc, omega, varpi, f, littleomega
   real(P), dimension(1001,7)    :: input_data, output_data
   integer                       :: count_start, count_rate, count_end, i, j
   character(len = 200)          :: file_in, file_out, path
   character(len=*), parameter   :: fmt = '(ES23.16,10(",",ES23.16,:))'
   call system_clock(count_start, count_rate)
   
   path = '/Users/carlislewishard/Documents/Classes/Year_4/Numerical_Dynamics/Module3/&
   module3-carlislewishard/data/'

   file_in = trim(path) // trim('id000010-EL.csv')
   file_out = trim(path) // trim('id000010-XV-new.csv')
   open(unit = 11, file = file_in, status = 'old')
   open(unit = 12, file = file_out, status = 'replace')
   read(11,*) !Skip the first line
   do i = 1,1001
      read(11,*) input_data(i,:)
   end do

   close(11)

   !Write the output header
   write(12,*) "t,px,py,pz,vx,vy,vz"

   !Calculate the constants and define the input units
   PI = 4.0_P * atan(1.0_P)
   G = 4.0_P * PI**2 
   msun = 1.0_P
   mu = G * msun

   !Name the input data
   t(:) = input_data(:,1)
   a(:) = input_data(:,2)
   e(:) = input_data(:,3)
   inc(:) = input_data(:,4)
   omega(:) = input_data(:,5)
   varpi(:) = input_data(:,6)
   f(:) = input_data(:,7)

   do j = 1,1001 

      !Calculate littleomega
      littleomega(j) = varpi(j) - omega(j)

      !Calculate Rx, Ry, Rz
      R_norm = (a(j) * (1.0_P - e(j)**2)) / ((e(j) * cos(f(j))) + 1.0_P)
      rdot = (a(j) * e(j) * sin(f(j))) / sqrt(1.0_P - e(j)**2)

      px = R_norm * ((cos(omega(j)) * cos(littleomega(j) + f(j))) - (sin(omega(j)) * sin(littleomega(j) + f(j)) * cos(inc(j)))) 
      py = R_norm * ((sin(omega(j)) * cos(littleomega(j) + f(j))) + (cos(omega(j)) * sin(littleomega(j) + f(j)) * cos(inc(j))))
      pz = R_norm * sin(littleomega(j) + f(j)) * sin(inc(j))

      !Calculate hx, hy, hz
      h_norm = sqrt(mu * a(j) * (1.0_P - e(j)**2))
      h(3) = h_norm * cos(inc(j))
      h(1) = sign(h_norm * sin(inc(j)) * sin(omega(j)), h(3))
      h(2) = sign(h_norm * sin(inc(j)) * cos(omega(j)), h(3))

      !Calculate vx, vy, vz
      V_norm = sqrt(mu * ((2.0_P / R_norm) - (1.0_P / a(j))))

      vx = ((rdot * cos(omega(j)) * cos(littleomega(j) + f(j))) - &
         (rdot * sin(omega(j)) * sin(littleomega(j) + f(j)) * cos(inc(j)))) 
      vy = ((rdot * sin(omega(j)) * cos(littleomega(j) + f(j))) + &
         (rdot * cos(omega(j)) * sin(littleomega(j) + f(j)) * cos(inc(j)))) 
      vz = (rdot * sin(littleomega(j) + f(j)) * sin(inc(j))) 

      if (j == 1) then
         write(*,*) 'el2xv', vx, vy, vz, V_norm
      end if 

      !Name the output data
      output_data(j,1) = t(j)
      output_data(j,2) = px
      output_data(j,3) = py
      output_data(j,4) = pz
      output_data(j,5) = vx / 365.25_P
      output_data(j,6) = vy / 365.25_P
      output_data(j,7) = vz / 365.25_P

      !Write to the output file
      write(12,fmt) output_data(j,:)

   end do

   close(12)

   call system_clock(count_end, count_rate)
   write(*,*) "INTEGRATION TIME (S) : ", ((count_end  - count_start) / (count_rate * 1.0_P))

end program el2xv