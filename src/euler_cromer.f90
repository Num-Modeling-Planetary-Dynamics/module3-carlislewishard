program euler_cromer
   use, intrinsic                :: iso_fortran_env
   implicit none
   integer, parameter            :: P = real64
   integer                       :: dt_max, count_start, count_rate, count_end, i
   real(P)                       :: PI, G, m_sun, m_n, m_p, h, tmax, t
   real(P)                       :: mag_r_n_0, mag_v_n_0, mag_r_p_0, mag_v_p_0
   real(P)                       :: energy_n_0, angmom_n_0, energy_p_0, angmom_p_0
   real(P)                       :: energy_n, energy_p, angmom_n, angmom_p
   real(P), dimension(7)         :: n_input_data, p_input_data
   real(P)                       :: r_n_x_old, r_n_y_old, r_n_z_old, v_n_x_old, v_n_y_old, v_n_z_old
   real(P)                       :: r_p_x_old, r_p_y_old, r_p_z_old, v_p_x_old, v_p_y_old, v_p_z_old
   real(P)                       :: r_n_x_new, r_n_y_new, r_n_z_new, v_n_x_new, v_n_y_new, v_n_z_new
   real(P)                       :: r_p_x_new, r_p_y_new, r_p_z_new, v_p_x_new, v_p_y_new, v_p_z_new
   real(P)                       :: mag_r_n, mag_v_n, mag_r_p, mag_v_p, dist_np
   real(P)                       :: v_n_x_kep, v_n_y_kep, v_n_z_kep, v_n_x_int, v_n_y_int, v_n_z_int
   real(P)                       :: v_p_x_kep, v_p_y_kep, v_p_z_kep, v_p_x_int, v_p_y_int, v_p_z_int
   character(len = 200)          :: n_file_in, p_file_in, n_file_out, p_file_out, path
   character(len=*), parameter   :: fmt = '(ES23.16,10(",",ES23.16,:))'

   call system_clock(count_start, count_rate)
   
   path = '/Users/carlislewishard/Documents/Classes/Year_4/Numerical_Dynamics/Module3/&
   module3-carlislewishard/data/'

   n_file_in = trim(path) // trim('id000009-XV.csv')
   p_file_in = trim(path) // trim('id000010-XV.csv')
   n_file_out = trim(path) // trim('id000009-XV-EC.csv')
   p_file_out = trim(path) // trim('id000010-XV-EC.csv')

   open(unit = 11, file = n_file_in, status = 'old')
   open(unit = 12, file = p_file_in, status = 'old')
   open(unit = 13, file = n_file_out, status = 'replace')
   open(unit = 14, file = p_file_out, status = 'replace')

   read(11,*) !Skip the first line
   read(12,*) !Skip the first line

   read(11,*) n_input_data(:)
   read(12,*) p_input_data(:)

   close(11)
   close(12)

   write(13,*) "t,px,py,pz,vx,vy,vz,energy,angmom"
   write(14,*) "t,px,py,pz,vx,vy,vz,energy,angmom"

   PI = 4.0_P * atan(1.0_P)                                                               !pi
   G = 4.0_P * PI**2                                                                      !gravitational constant in solar masses
   m_sun = 1.0_P                                                                          !mass of the sun in solar masses
   m_n = 0.00005149_P                                                                     !mass of neptune in solar masses
   m_p = 0.00000000658086572_P                                                            !mass of pluto in solar masses

   h = 1.0_P                                                                              !timestep in years
   tmax = 10.0_P                                                                      !total simulation time in years
   dt_max = int(tmax / h)                                                                 !total number of timesteps

   mag_r_n_0 = sqrt(n_input_data(2)**2 + n_input_data(3)**2 + n_input_data(4)**2)         !|r| for neptune at t=0
   mag_v_n_0 = sqrt(n_input_data(5)**2 + n_input_data(6)**2 + n_input_data(7)**2)         !|v| for neptune at t=0

   mag_r_p_0 = sqrt(p_input_data(2)**2 + p_input_data(3)**2 + p_input_data(4)**2)         !|r| for pluto at t=0
   mag_v_p_0 = sqrt(p_input_data(5)**2 + p_input_data(6)**2 + p_input_data(7)**2)         !|v| for pluto at t=0

   energy_n_0 = 0.5_P * mag_v_n_0**2 - ((G * m_sun) / mag_r_n_0)                          !energy from vis viva for neptune at t=0
   angmom_n_0 = (n_input_data(3) * n_input_data(7) - n_input_data(4) * n_input_data(6)) + &
      (n_input_data(4) * n_input_data(5) - n_input_data(2) * n_input_data(7)) + &
      (n_input_data(2) * n_input_data(6) - n_input_data(3) * n_input_data(5))             !angmom from r cros v for neptune at t=0

   energy_p_0 = 0.5_P * mag_v_p_0**2 - ((G * m_sun) / mag_r_p_0)                          !energy from vis viva for pluto at t=0
   angmom_p_0 = (p_input_data(3) * p_input_data(7) - p_input_data(4) * p_input_data(6)) + &
      (p_input_data(4) * p_input_data(5) - p_input_data(2) * p_input_data(7)) + &
      (p_input_data(2) * p_input_data(6) - p_input_data(3) * p_input_data(5))             !angmom from r cros v for pluto at t=0

   write(13,fmt) 0.0, n_input_data(2), n_input_data(3), n_input_data(4), &
      n_input_data(5), n_input_data(6), n_input_data(7), energy_n_0, angmom_n_0
   write(14,fmt) 0.0, p_input_data(2), p_input_data(3), p_input_data(4), &
      p_input_data(5), p_input_data(6), p_input_data(7), energy_p_0, angmom_p_0

   r_n_x_old = n_input_data(2)
   r_n_y_old = n_input_data(3)
   r_n_z_old = n_input_data(4)
   v_n_x_old = n_input_data(5)
   v_n_y_old = n_input_data(6)
   v_n_z_old = n_input_data(7)

   r_p_x_old = p_input_data(2)
   r_p_y_old = p_input_data(3)
   r_p_z_old = p_input_data(4)
   v_p_x_old = p_input_data(5)
   v_p_y_old = p_input_data(6)
   v_p_z_old = p_input_data(7)

   do i = 1, dt_max

      mag_r_n = sqrt(r_n_x_old**2 + r_n_y_old**2 + r_n_z_old**2)
      mag_r_p = sqrt(r_p_x_old**2 + r_p_y_old**2 + r_p_z_old**2)
      dist_np = sqrt((r_n_x_old - r_p_x_old)**2 + (r_n_y_old - r_p_y_old)**2 + (r_n_z_old - r_p_z_old)**2)

      v_n_x_kep = - G * m_sun * (r_n_x_old / (mag_r_n)**3)
      v_n_x_int = - (G * (m_p + m_n) * (r_n_x_old - r_p_x_old)) / dist_np**3
      v_n_x_new = v_n_x_old + (h * (v_n_x_kep + v_n_x_int))

      v_n_y_kep = - G * m_sun * (r_n_y_old / (mag_r_n)**3)
      v_n_y_int = - (G * (m_p + m_n) * (r_n_y_old - r_p_y_old)) / dist_np**3
      v_n_y_new = v_n_y_old + (h * (v_n_y_kep + v_n_y_int))

      v_n_z_kep = - G * m_sun * (r_n_z_old / (mag_r_n)**3)
      v_n_z_int = - (G * (m_p + m_n) * (r_n_z_old - r_p_z_old)) / dist_np**3
      v_n_z_new = v_n_z_old + (h * (v_n_z_kep + v_n_z_int))

      v_p_x_kep = - G * m_sun * (r_p_x_old / (mag_r_p)**3)
      v_p_x_int = - (G * (m_p + m_n) * (r_p_x_old - r_n_x_old)) / dist_np**3
      v_p_x_new = v_p_x_old + (h * (v_p_x_kep + v_p_x_int))

      v_p_y_kep = - G * m_sun * (r_p_y_old / (mag_r_p)**3)
      v_p_y_int = - (G * (m_p + m_n) * (r_p_y_old - r_n_y_old)) / dist_np**3
      v_p_y_new = v_p_y_old + (h * (v_p_y_kep + v_p_y_int))

      v_p_z_kep = - G * m_sun * (r_p_z_old / (mag_r_p)**3)
      v_p_z_int = - (G * (m_p + m_n) * (r_p_z_old - r_n_z_old)) / dist_np**3
      v_p_z_new = v_p_z_old + (h * (v_p_z_kep + v_p_z_int))

      r_n_x_new = r_n_x_old + (v_n_x_new * h)
      r_n_y_new = r_n_y_old + (v_n_y_new * h)
      r_n_z_new = r_n_z_old + (v_n_z_new * h)

      r_p_x_new = r_p_x_old + (v_p_x_new * h)
      r_p_y_new = r_p_y_old + (v_p_y_new * h)
      r_p_z_new = r_p_z_old + (v_p_z_new * h)

      mag_r_n = sqrt(r_n_x_new**2 + r_n_y_new**2 + r_n_z_new**2)
      mag_r_p = sqrt(r_p_x_new**2 + r_p_y_new**2 + r_p_z_new**2)

      mag_v_n = sqrt(v_n_x_new**2 + v_n_y_new**2 + v_n_z_new**2)
      mag_v_p = sqrt(v_p_x_new**2 + v_p_y_new**2 + v_p_z_new**2)

      energy_n = 0.5_P * mag_v_n**2 - ((G * m_sun) / mag_r_n)
      energy_p = 0.5_P * mag_v_p**2 - ((G * m_sun) / mag_r_p)

      angmom_n = (r_n_y_new * v_n_z_new - r_n_z_new * v_n_y_new) + &
         (r_n_z_new * v_n_x_new - r_n_x_new * v_n_z_new) + &
         (r_n_x_new * v_n_y_new - r_n_y_new * v_n_x_new) 

      angmom_p = (r_p_y_new * v_p_z_new - r_p_z_new * v_p_y_new) + &
         (r_p_z_new * v_p_x_new - r_p_x_new * v_p_z_new) + &
         (r_p_x_new * v_p_y_new - r_p_y_new * v_p_x_new) 

      t = i * h
      write(*,*) "Time = ", t

      write(13,fmt) t, r_n_x_old, r_n_y_old, r_n_z_old, v_n_x_old, v_n_y_old, v_n_z_old, energy_n, angmom_n
      write(14,fmt) t, r_p_x_old, r_p_y_old, r_p_z_old, v_p_x_old, v_p_y_old, v_p_z_old, energy_p, angmom_p

      r_n_x_old = r_n_x_new
      r_n_y_old = r_n_y_new
      r_n_z_old = r_n_z_new

      v_n_x_old = v_n_x_new
      v_n_y_old = v_n_y_new
      v_n_z_old = v_n_z_new

      r_p_x_old = r_p_x_new
      r_p_y_old = r_p_y_new
      r_p_z_old = r_p_z_new

      v_p_x_old = v_p_x_new
      v_p_y_old = v_p_y_new
      v_p_z_old = v_p_z_new

   end do 

   close(13)
   close(14)

   call system_clock(count_end, count_rate)
   write(*,*) "INTEGRATION TIME (S) : ", ((count_end - count_start) / (count_rate * 1.0_P))

end program euler_cromer