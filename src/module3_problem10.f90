program module3_problem10
   use, intrinsic                :: iso_fortran_env
   implicit none
   integer, parameter            :: P = real64
   integer                       :: dt_max, count_start, count_rate, count_end, i, j
   real(P)                       :: PI, G, m_sun, m_nep, m_plu, m_sys, h, tmax, t, t_half
   real(P)                       :: nep_r_norm, nep_v_norm, nep_h_norm, plu_r_norm, plu_v_norm, plu_h_norm 
   real(P)                       :: nep_a, nep_e, nep_n, nep_M_new, plu_a, plu_e, plu_n, plu_M_new
   real(P)                       :: nep_func, nep_func_1d, nep_func_2d, nep_func_3d
   real(P)                       :: plu_func, plu_func_1d, plu_func_2d, plu_func_3d
   real(P)                       :: nep_del1, nep_del2, nep_del3, plu_del1, plu_del2, plu_del3
   real(P)                       :: nep_E_new, plu_E_new, nep_r_norm3, plu_r_norm3
   real(P)                       :: nep_func_x, nep_func_y, nep_func_z, plu_func_x, plu_func_y, plu_func_z
   real(P), dimension(2,7)       :: input_data
   real(P), dimension(3)         :: vb_sun, vb_nep, vb_plu, vh_sun, vh_nep, vh_plu
   real(P), dimension(3)         :: nep_r_drift, plu_r_drift, nep_v_kick, plu_v_kick
   real(P), dimension(6)         :: nep_old, plu_old
   real(P), dimension(11)        :: nep_E_val, plu_E_val
   character(len = 200)          :: path, nep_file_in, plu_file_in, nep_file_out, plu_file_out
   character(len=*), parameter   :: fmt = '(ES23.16,10(",",ES23.16,:))'

   call system_clock(count_start, count_rate)

   path = '/Users/carlislewishard/Documents/Classes/Year_4/Numerical_Dynamics/Module3/&
   module3-carlislewishard/data/'

   nep_file_in = trim(path) // trim('id000009-XV.csv')
   plu_file_in = trim(path) // trim('id000010-XV.csv')

   nep_file_out = trim(path) // trim('id000009-XV-INT.csv')
   plu_file_out = trim(path) // trim('id000010-XV-INT.csv')

   open(unit = 11, file = nep_file_in, status = 'old')
   open(unit = 12, file = plu_file_in, status = 'old')

   open(unit = 21, file = nep_file_out, status = 'replace')
   open(unit = 22, file = plu_file_out, status = 'replace')

   !Skip the header line
   read(11,*)
   read(12,*)

   !Read the first line of data
   read(11,*) input_data(1,:)
   read(12,*) input_data(2,:)

   close(11)
   close(12)

   PI = 4.0_P * atan(1.0_P)                                                               !pi
   G = 4.0_P * PI**2                                                                      !gravitational constant in solar masses
   m_sun = 1.0_P
   m_nep = 0.00005149_P
   m_plu = 6.58086572E-9_P                                                                !mass of the sun in solar masses     
   m_sys = m_sun + m_nep + m_plu

   h = 0.008_P                                                                            !timestep in years
   tmax = 100000.0_P                                                                      !total simulation time in years
   dt_max = int(tmax / h)                                                                 !total number of timesteps

   input_data(:,5) = input_data(:,5) * 365.25_P                                           !AU/day -> AU/year
   input_data(:,6) = input_data(:,6) * 365.25_P                                           !AU/day -> AU/year
   input_data(:,7) = input_data(:,7) * 365.25_P                                           !AU/day -> AU/year

   !Write t=0 positions and velocities in original coordinate frame
   write(21,fmt) 0.0, input_data(1,2), input_data(1,3), input_data(1,4), input_data(1,5), input_data(1,6), input_data(1,7)
   write(22,fmt) 0.0, input_data(2,2), input_data(2,3), input_data(2,4), input_data(2,5), input_data(2,6), input_data(2,7)

   !Switch to the velocities to barycentric coordinates
   vb_sun(:) = - ((m_nep * input_data(1,5:7)) + (m_plu * input_data(2,5:7))) / m_sys
   vb_nep(:) = input_data(1,5:7) + vb_sun(:)
   vb_plu(:) = input_data(2,5:7) + vb_sun(:)

   nep_old(1) = input_data(1,2)
   nep_old(2) = input_data(1,3)
   nep_old(3) = input_data(1,4)
   nep_old(4) = vb_nep(1)
   nep_old(5) = vb_nep(2)
   nep_old(6) = vb_nep(3)

   plu_old(1) = input_data(2,2)
   plu_old(2) = input_data(2,3)
   plu_old(3) = input_data(2,4)
   plu_old(4) = vb_plu(1)
   plu_old(5) = vb_plu(2)
   plu_old(6) = vb_plu(3)

   do i = 1, dt_max

      t = i * h
      t_half = 0.5_P * t

      !Switch from XV->EL
      nep_r_norm = sqrt(nep_old(1)**2 + nep_old(2)**2 + nep_old(3)**2)
      nep_v_norm = sqrt(nep_old(4)**2 + nep_old(5)**2 + nep_old(6)**2)
      nep_h_norm = sqrt(((nep_old(2) * nep_old(6)) - (nep_old(3) * nep_old(5)))**2 &
         + ((nep_old(3) * nep_old(4)) - (nep_old(1) * nep_old(6)))**2 &
         + ((nep_old(1) * nep_old(5)) - (nep_old(2) * nep_old(4)))**2)
      nep_a = ((2.0_P / nep_r_norm) - (nep_v_norm**2 / (G * m_sun)))**(-1)
      nep_e = sqrt(1.0_P - (nep_h_norm**2 / (nep_a * G * m_sun)))
      nep_n = (nep_a**3 / (G * m_sun))**(-0.5_P)
      nep_M_new = mod(nep_n * t_half, 2.0_P * PI) 
      nep_E_val(1) = acos((1.0_P - (nep_r_norm / nep_a)) / nep_e)
      if ((nep_E_val(1) + PI) > (2.0_P * PI)) then
         nep_E_val(1) = nep_E_val(1) - PI
      end if 
      plu_r_norm = sqrt(plu_old(1)**2 + plu_old(2)**2 + plu_old(3)**2)
      plu_v_norm = sqrt(plu_old(4)**2 + plu_old(5)**2 + plu_old(6)**2)
      plu_h_norm = sqrt(((plu_old(2) * plu_old(6)) - (plu_old(3) * plu_old(5)))**2 &
         + ((plu_old(3) * plu_old(4)) - (plu_old(1) * plu_old(6)))**2 &
         + ((plu_old(1) * plu_old(5)) - (plu_old(2) * plu_old(4)))**2)
      plu_a = ((2.0_P / plu_r_norm) - (plu_v_norm**2 / (G * m_sun)))**(-1)
      plu_e = sqrt(1.0_P - (plu_h_norm**2 / (plu_a * G * m_sun)))
      plu_n = (plu_a**3 / (G * m_sun))**(-0.5_P)
      plu_M_new = mod(plu_n * t_half, 2.0_P * PI) 
      plu_E_val(1) = acos((1.0_P - (plu_r_norm / plu_a)) / plu_e)
      if ((plu_E_val(1) + PI) > (2.0_P * PI)) then
         plu_E_val(1) = plu_E_val(1) - PI
      end if 

      !Drift positions using Danby 
      do j = 1, 10
         !Calculate the first, second, and third derivatives of Kepler's Equation
         nep_func = nep_E_val(j) - (nep_e * sin(nep_E_val(j))) - nep_M_new
         nep_func_1d = 1.0_P - (nep_e * cos(nep_E_val(j)))
         nep_func_2d = nep_e * sin(nep_E_val(j))
         nep_func_3d = nep_e * cos(nep_E_val(j))

         plu_func = plu_E_val(j) - (plu_e * sin(plu_E_val(j))) - plu_M_new
         plu_func_1d = 1.0_P - (plu_e * cos(plu_E_val(j)))
         plu_func_2d = plu_e * sin(plu_E_val(j))
         plu_func_3d = plu_e * cos(plu_E_val(j))

         !Calculate the corresponding delta functions
         nep_del1 = - nep_func / nep_func_1d
         nep_del2 = - nep_func / (nep_func_1d + (0.5_P * nep_del1 * nep_func_2d))
         nep_del3 = - nep_func / (nep_func_1d + (0.5_P * nep_del2 * nep_func_2d) &
            + ((1.0_P/6.0_P) * nep_del2**2.0_P * nep_func_3d))

         plu_del1 = - plu_func / plu_func_1d
         plu_del2 = - plu_func / (plu_func_1d + (0.5_P * plu_del1 * plu_func_2d))
         plu_del3 = - plu_func / (plu_func_1d + (0.5_P * plu_del2 * plu_func_2d) &
            + ((1.0_P/6.0_P) * plu_del2**2.0_P * plu_func_3d))

         !Calculate the eccentric anomaly in radians
         nep_E_val(j+1) = (nep_E_val(j) + nep_del3) 

         plu_E_val(j+1) = (plu_E_val(j) + plu_del3)        
      end do

      nep_E_new = nep_E_val(11)

      plu_E_new = plu_E_val(11)

      !Calculate the new positions
      nep_r_drift(1) = nep_a * (cos(nep_E_new) - nep_e)
      nep_r_drift(2) = nep_a * sqrt(1.0_P - nep_e**2.0_P) * sin(nep_E_new)
      nep_r_drift(3) = nep_old(3) + (0.5_P * h * nep_old(6))

      plu_r_drift(1) = plu_a * (cos(plu_E_new) - plu_e)
      plu_r_drift(2) = plu_a * sqrt(1.0_P - plu_e**2.0_P) * sin(plu_E_new)
      plu_r_drift(3) = plu_old(3) + (0.5_P * h * plu_old(6))

      !v_sys(:) = ((m_nep * nep_old(4:6)) + (m_plu * plu_old(4:6))) / m_sun
      !nep_r_drift(1:3) = nep_old(1:3) + (v_sys(:) * dt)
      !plu_r_drift(1:3) = plu_old(1:3) + (v_sys(:) * dt)

      !Kick velocities using leapfrog
      nep_r_norm3 = (nep_r_drift(1)**2 + nep_r_drift(2)**2 + nep_r_drift(3)**2)**(3.0_P/2.0_P)

      plu_r_norm3 = (plu_r_drift(1)**2 + plu_r_drift(2)**2 + plu_r_drift(3)**2)**(3.0_P/2.0_P)

      nep_func_x = - (G * m_sun) * (nep_r_drift(1) / nep_r_norm3)
      nep_func_y = - (G * m_sun) * (nep_r_drift(2) / nep_r_norm3)
      nep_func_z = - (G * m_sun) * (nep_r_drift(3) / nep_r_norm3)

      plu_func_x = - (G * m_sun) * (plu_r_drift(1) / plu_r_norm3)
      plu_func_y = - (G * m_sun) * (plu_r_drift(2) / plu_r_norm3)
      plu_func_z = - (G * m_sun) * (plu_r_drift(3) / plu_r_norm3)

      nep_v_kick(1) = nep_old(4) + (h * nep_func_x)
      nep_v_kick(2) = nep_old(5) + (h * nep_func_y)
      nep_v_kick(3) = nep_old(6) + (h * nep_func_z)

      plu_v_kick(1) = plu_old(4) + (h * plu_func_x)
      plu_v_kick(2) = plu_old(5) + (h * plu_func_y)
      plu_v_kick(3) = plu_old(6) + (h * plu_func_z)

      !nep_v_kick(1:3) = nep_old(4:6) + ( * dt)

      !Drift positions using Danby
      nep_r_norm = sqrt(nep_old(1)**2 + nep_old(2)**2 + nep_old(3)**2)
      nep_v_norm = sqrt(nep_v_kick(1)**2 + nep_v_kick(2)**2 + nep_v_kick(3)**2)
      nep_h_norm = sqrt(((nep_old(2) * nep_v_kick(3)) - (nep_old(3) * nep_v_kick(2)))**2 &
         + ((nep_old(3) * nep_v_kick(1)) - (nep_old(1) * nep_v_kick(3)))**2 &
         + ((nep_old(1) * nep_v_kick(3)) - (nep_old(2) * nep_v_kick(1)))**2)
      nep_a = ((2.0_P / nep_r_norm) - (nep_v_norm**2 / (G * m_sun)))**(-1)
      nep_e = sqrt(1.0_P - (nep_h_norm**2 / (nep_a * G * m_sun)))
      nep_n = (nep_a**3 / (G * m_sun))**(-0.5_P)
      nep_M_new = mod(nep_n * t, 2.0_P * PI) 
      nep_E_val(1) = acos((1.0_P - (nep_r_norm / nep_a)) / nep_e)
      if ((nep_E_val(1) + PI) > (2.0_P * PI)) then
         nep_E_val(1) = nep_E_val(1) - PI
      end if 

      plu_r_norm = sqrt(plu_old(1)**2 + plu_old(2)**2 + plu_old(3)**2)
      plu_v_norm = sqrt(plu_v_kick(1)**2 + plu_v_kick(2)**2 + plu_v_kick(3)**2)
      plu_h_norm = sqrt(((plu_old(2) * plu_v_kick(3)) - (plu_old(3) * plu_v_kick(2)))**2 &
         + ((plu_old(3) * plu_v_kick(1)) - (plu_old(1) * plu_v_kick(3)))**2 &
         + ((plu_old(1) * plu_v_kick(3)) - (plu_old(2) * plu_v_kick(1)))**2)
      plu_a = ((2.0_P / plu_r_norm) - (plu_v_norm**2 / (G * m_sun)))**(-1)
      plu_e = sqrt(1.0_P - (plu_h_norm**2 / (plu_a * G * m_sun)))
      plu_n = (plu_a**3 / (G * m_sun))**(-0.5_P)
      plu_M_new = mod(plu_n * t, 2.0_P * PI) 
      plu_E_val(1) = acos((1.0_P - (plu_r_norm / plu_a)) / plu_e)
      if ((plu_E_val(1) + PI) > (2.0_P * PI)) then
         plu_E_val(1) = plu_E_val(1) - PI
      end if 
 
      do j = 1, 10
         !Calculate the first, second, and third derivatives of Kepler's Equation
         nep_func = nep_E_val(j) - (nep_e * sin(nep_E_val(j))) - nep_M_new
         nep_func_1d = 1.0_P - (nep_e * cos(nep_E_val(j)))
         nep_func_2d = nep_e * sin(nep_E_val(j))
         nep_func_3d = nep_e * cos(nep_E_val(j))

         plu_func = plu_E_val(j) - (plu_e * sin(plu_E_val(j))) - plu_M_new
         plu_func_1d = 1.0_P - (plu_e * cos(plu_E_val(j)))
         plu_func_2d = plu_e * sin(plu_E_val(j))
         plu_func_3d = plu_e * cos(plu_E_val(j))

         !Calculate the corresponding delta functions
         nep_del1 = - nep_func / nep_func_1d
         nep_del2 = - nep_func / (nep_func_1d + (0.5_P * nep_del1 * nep_func_2d))
         nep_del3 = - nep_func / (nep_func_1d + (0.5_P * nep_del2 * nep_func_2d) &
            + ((1.0_P/6.0_P) * nep_del2**2.0_P * nep_func_3d))

         plu_del1 = - plu_func / plu_func_1d
         plu_del2 = - plu_func / (plu_func_1d + (0.5_P * plu_del1 * plu_func_2d))
         plu_del3 = - plu_func / (plu_func_1d + (0.5_P * plu_del2 * plu_func_2d) &
            + ((1.0_P/6.0_P) * plu_del2**2.0_P * plu_func_3d))

         !Calculate the eccentric anomaly
         nep_E_val(j+1) = (nep_E_val(j) + nep_del3) !eccentric anomaly in radians

         plu_E_val(j+1) = (plu_E_val(j) + plu_del3) !eccentric anomaly in radians         
      end do

      nep_E_new = nep_E_val(11)

      plu_E_new = plu_E_val(11)

      !Calculate the new positions
      nep_r_drift(1) = nep_a * (cos(nep_E_new) - nep_e)
      nep_r_drift(2) = nep_a * sqrt(1.0_P - nep_e**2.0_P) * sin(nep_E_new)
      nep_r_drift(3) = nep_old(3) + (0.5_P * h * nep_v_kick(3))

      plu_r_drift(1) = plu_a * (cos(plu_E_new) - plu_e)
      plu_r_drift(2) = plu_a * sqrt(1.0_P - plu_e**2.0_P) * sin(plu_E_new)
      plu_r_drift(3) = plu_old(3) + (0.5_P * h * plu_v_kick(3))

      !Switch back to original coordinates
      vh_sun(:) = - ((m_nep * nep_v_kick(:)) + (m_plu * plu_v_kick(:))) / m_sun
      vh_nep(:) = nep_v_kick(:) - vb_sun(:)
      vh_plu(:) = plu_v_kick(:) - vb_sun(:)

      !Write new positions and velocities in original coordinate frame
      if (mod(t, 10000.0_P) == 0) then
         write(*,*) "Time = ", t
         write(21,fmt) t, nep_r_drift(1), nep_r_drift(2), nep_r_drift(3), vh_nep(1), &
            vh_nep(2), vh_nep(3)
         write(22,fmt) t, plu_r_drift(1), plu_r_drift(2), plu_r_drift(3), vh_plu(1), &
            vh_plu(2), vh_plu(3)
      end if

      nep_old(1:3) = nep_r_drift(:)
      nep_old(4:6) = nep_v_kick(:)
      plu_old(1:3) = plu_r_drift(:)
      plu_old(4:6) = plu_v_kick(:)

   end do

   close(21)
   close(22)

   call system_clock(count_end, count_rate)
   write(*,*) "INTEGRATION TIME (S) : ", ((count_end - count_start) / (count_rate * 1.0_P))

end program module3_problem10