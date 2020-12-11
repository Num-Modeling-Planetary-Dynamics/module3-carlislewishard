program leapfrog
   use, intrinsic                :: iso_fortran_env
   implicit none
   integer, parameter            :: P = real64
   integer                       :: dt_max, count_start, count_rate, count_end, i
   real(P)                       :: PI, G, m_sun, h, tmax, t
   real(P), dimension(9,7)       :: input_data
   real(P), dimension(9)         :: r_x_old, r_y_old, r_z_old, v_x_old, v_y_old, v_z_old
   real(P), dimension(9)         :: r_x_new, r_y_new, r_z_new, v_x_new, v_y_new, v_z_new
   real(P), dimension(9)         :: r_x_new_half, r_y_new_half, r_z_new_half
   real(P), dimension(9)         :: func_x, func_y, func_z, r_norm3
   character(len = 200)          :: mer_file_in, ven_file_in, ear_file_in, mar_file_in
   character(len = 200)          :: jup_file_in, sat_file_in, ura_file_in, nep_file_in, plu_file_in
   character(len = 200)          :: mer_file_out, ven_file_out, ear_file_out, mar_file_out
   character(len = 200)          :: jup_file_out, sat_file_out, ura_file_out, nep_file_out, plu_file_out
   character(len = 200)          :: path
   character(len=*), parameter   :: fmt = '(ES23.16,10(",",ES23.16,:))'

   call system_clock(count_start, count_rate)

   path = '/Users/carlislewishard/Documents/Classes/Year_4/Numerical_Dynamics/Module3/&
   module3-carlislewishard/data/'

   mer_file_in = trim(path) // trim('id000002-XV.csv')
   ven_file_in = trim(path) // trim('id000003-XV.csv')
   ear_file_in = trim(path) // trim('id000004-XV.csv')
   mar_file_in = trim(path) // trim('id000005-XV.csv')
   jup_file_in = trim(path) // trim('id000006-XV.csv')
   sat_file_in = trim(path) // trim('id000007-XV.csv')
   ura_file_in = trim(path) // trim('id000008-XV.csv')
   nep_file_in = trim(path) // trim('id000009-XV.csv')
   plu_file_in = trim(path) // trim('id000010-XV.csv')

   mer_file_out = trim(path) // trim('id000002-XV-LF.csv')
   ven_file_out = trim(path) // trim('id000003-XV-LF.csv')
   ear_file_out = trim(path) // trim('id000004-XV-LF.csv')
   mar_file_out = trim(path) // trim('id000005-XV-LF.csv')
   jup_file_out = trim(path) // trim('id000006-XV-LF.csv')
   sat_file_out = trim(path) // trim('id000007-XV-LF.csv')
   ura_file_out = trim(path) // trim('id000008-XV-LF.csv')
   nep_file_out = trim(path) // trim('id000009-XV-LF.csv')
   plu_file_out = trim(path) // trim('id000010-XV-LF.csv')

   open(unit = 11, file = mer_file_in, status = 'old')
   open(unit = 12, file = ven_file_in, status = 'old')
   open(unit = 13, file = ear_file_in, status = 'old')
   open(unit = 14, file = mar_file_in, status = 'old')
   open(unit = 15, file = jup_file_in, status = 'old')
   open(unit = 16, file = sat_file_in, status = 'old')
   open(unit = 17, file = ura_file_in, status = 'old')
   open(unit = 18, file = nep_file_in, status = 'old')
   open(unit = 19, file = plu_file_in, status = 'old')

   open(unit = 21, file = mer_file_out, status = 'replace')
   open(unit = 22, file = ven_file_out, status = 'replace')
   open(unit = 23, file = ear_file_out, status = 'replace')
   open(unit = 24, file = mar_file_out, status = 'replace')
   open(unit = 25, file = jup_file_out, status = 'replace')
   open(unit = 26, file = sat_file_out, status = 'replace')
   open(unit = 27, file = ura_file_out, status = 'replace')
   open(unit = 28, file = nep_file_out, status = 'replace')
   open(unit = 29, file = plu_file_out, status = 'replace')

   !Skip the header line
   read(11,*)
   read(12,*)
   read(13,*)
   read(14,*)
   read(15,*)
   read(16,*)
   read(17,*)
   read(18,*)
   read(19,*)

   !Read the first line of data
   read(11,*) input_data(1,:)
   read(12,*) input_data(2,:)
   read(13,*) input_data(3,:)
   read(14,*) input_data(4,:)
   read(15,*) input_data(5,:)
   read(16,*) input_data(6,:)
   read(17,*) input_data(7,:)
   read(18,*) input_data(8,:)
   read(19,*) input_data(9,:)

   close(11)
   close(12)
   close(13)
   close(14)
   close(15)
   close(16)
   close(17)
   close(18)
   close(19)

   write(21,*) "t,px,py,pz,vx,vy,vz"
   write(22,*) "t,px,py,pz,vx,vy,vz"
   write(23,*) "t,px,py,pz,vx,vy,vz"
   write(24,*) "t,px,py,pz,vx,vy,vz"
   write(25,*) "t,px,py,pz,vx,vy,vz"
   write(26,*) "t,px,py,pz,vx,vy,vz"
   write(27,*) "t,px,py,pz,vx,vy,vz"
   write(28,*) "t,px,py,pz,vx,vy,vz"
   write(29,*) "t,px,py,pz,vx,vy,vz"

   PI = 4.0_P * atan(1.0_P)                                                               !pi
   G = 4.0_P * PI**2                                                                      !gravitational constant in solar masses
   m_sun = 1.0_P                                                                          !mass of the sun in solar masses     

   h = 0.008_P                                                                            !timestep in years
   tmax = 1000000.0_P                                                                     !total simulation time in years
   dt_max = int(tmax / h)                                                                 !total number of timesteps

   input_data(:,5) = input_data(:,5) * 365.25_P                                           !AU/day -> AU/year
   input_data(:,6) = input_data(:,6) * 365.25_P                                           !AU/day -> AU/year
   input_data(:,7) = input_data(:,7) * 365.25_P                                           !AU/day -> AU/year

   write(21,fmt) 0.0, input_data(1,2), input_data(1,3), input_data(1,4), input_data(1,5), input_data(1,6), input_data(1,7)
   write(22,fmt) 0.0, input_data(2,2), input_data(2,3), input_data(2,4), input_data(2,5), input_data(2,6), input_data(2,7)
   write(23,fmt) 0.0, input_data(3,2), input_data(3,3), input_data(3,4), input_data(3,5), input_data(3,6), input_data(3,7)
   write(24,fmt) 0.0, input_data(4,2), input_data(4,3), input_data(4,4), input_data(4,5), input_data(4,6), input_data(4,7)
   write(25,fmt) 0.0, input_data(5,2), input_data(5,3), input_data(5,4), input_data(5,5), input_data(5,6), input_data(5,7)
   write(26,fmt) 0.0, input_data(6,2), input_data(6,3), input_data(6,4), input_data(6,5), input_data(6,6), input_data(6,7)
   write(27,fmt) 0.0, input_data(7,2), input_data(7,3), input_data(7,4), input_data(7,5), input_data(7,6), input_data(7,7)
   write(28,fmt) 0.0, input_data(8,2), input_data(8,3), input_data(8,4), input_data(8,5), input_data(8,6), input_data(8,7)
   write(29,fmt) 0.0, input_data(9,2), input_data(9,3), input_data(9,4), input_data(9,5), input_data(9,6), input_data(9,7)

   r_x_old(:) = input_data(:,2)
   r_y_old(:) = input_data(:,3)
   r_z_old(:) = input_data(:,4)
   v_x_old(:) = input_data(:,5)
   v_y_old(:) = input_data(:,6)
   v_z_old(:) = input_data(:,7)

   do i = 1, dt_max

      !something is wrong here 
      r_x_new_half(:) = r_x_old(:) + (0.5_P * h * v_x_old(:))
      r_y_new_half(:) = r_y_old(:) + (0.5_P * h * v_y_old(:))
      r_z_new_half(:) = r_z_old(:) + (0.5_P * h * v_z_old(:))

      !r_norm3(:) = (r_x_new_half(:)**2 + r_y_new_half(:)**2 + r_z_new_half(:)**2)**(3.0_P/2.0_P)
      r_norm3(:) = (r_x_old(:)**2 + r_y_old(:)**2 + r_z_old(:)**2)**(3.0_P/2.0_P)

      !func_x(:) = - (G * m_sun) * (r_x_new_half(:) / r_norm3(:))
      !func_y(:) = - (G * m_sun) * (r_y_new_half(:) / r_norm3(:))
      !func_z(:) = - (G * m_sun) * (r_z_new_half(:) / r_norm3(:))

      func_x(:) = - (G * m_sun) * (r_x_old(:) / r_norm3(:))
      func_y(:) = - (G * m_sun) * (r_y_old(:) / r_norm3(:))
      func_z(:) = - (G * m_sun) * (r_z_old(:) / r_norm3(:))

      v_x_new(:) = v_x_old(:) + (h * func_x(:))
      v_y_new(:) = v_y_old(:) + (h * func_y(:))
      v_z_new(:) = v_z_old(:) + (h * func_z(:))

      r_x_new(:) = r_x_old(:) + (0.5_P * h * v_x_new(:))
      r_y_new(:) = r_y_old(:) + (0.5_P * h * v_y_new(:))
      r_z_new(:) = r_z_old(:) + (0.5_P * h * v_z_new(:))

      t = i * h

      if (mod(t, 10000.0_P) == 0) then
         write(*,*) "Time = ", t
         write(21,fmt) t, r_x_new(1), r_y_new(1), r_z_new(1), v_x_new(1), v_y_new(1), v_z_new(1)
         write(22,fmt) t, r_x_new(2), r_y_new(2), r_z_new(2), v_x_new(2), v_y_new(2), v_z_new(2)
         write(23,fmt) t, r_x_new(3), r_y_new(3), r_z_new(3), v_x_new(3), v_y_new(3), v_z_new(3)
         write(24,fmt) t, r_x_new(4), r_y_new(4), r_z_new(4), v_x_new(4), v_y_new(4), v_z_new(4)
         write(25,fmt) t, r_x_new(5), r_y_new(5), r_z_new(5), v_x_new(5), v_y_new(5), v_z_new(5)
         write(26,fmt) t, r_x_new(6), r_y_new(6), r_z_new(6), v_x_new(6), v_y_new(6), v_z_new(6)
         write(27,fmt) t, r_x_new(7), r_y_new(7), r_z_new(7), v_x_new(7), v_y_new(7), v_z_new(7)
         write(28,fmt) t, r_x_new(8), r_y_new(8), r_z_new(8), v_x_new(8), v_y_new(8), v_z_new(8)
         write(29,fmt) t, r_x_new(9), r_y_new(9), r_z_new(9), v_x_new(9), v_y_new(9), v_z_new(9)
      end if

      r_x_old(:) = r_x_new(:)
      r_y_old(:) = r_y_new(:)
      r_z_old(:) = r_z_new(:)
      v_x_old(:) = v_x_new(:) 
      v_y_old(:) = v_y_new(:)
      v_z_old(:) = v_z_new(:)

   end do

   close(21)
   close(22)
   close(23)
   close(24)
   close(25)
   close(26)
   close(27)
   close(28)
   close(29)

   call system_clock(count_end, count_rate)
   write(*,*) "INTEGRATION TIME (S) : ", ((count_end - count_start) / (count_rate * 1.0_P))

end program leapfrog