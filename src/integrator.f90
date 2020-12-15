module integrator
use, intrinsic                :: iso_fortran_env
implicit none
integer, parameter                         :: P = real64
real(P), parameter                         :: PI = 3.14159265_P
real(P), parameter                         :: G = 39.47841760_P
real(P), parameter                         :: m_sun = 1.0_P
real(P), parameter                         :: m_nep = 0.00005149_P
real(P), parameter                         :: m_plu = 6.58086572E-9_P  
real(P), parameter                         :: dt = 0.008_P      
real(P), parameter                         :: tmax = 35.17_P !100000.0_P
character(len=*), parameter                :: fmt = '(ES23.16,10(",",ES23.16,:))'

contains 

   subroutine read_input(path, in1, in2, out1, out2, input_data)
      implicit none
      character(len = 200), intent(in)     :: path, in1, in2, out1, out2
      real(P), dimension(2,7), intent(out) :: input_data
      character(len = 200)                 :: file_in1, file_in2, file_out1, file_out2

      file_in1 = trim(path) // trim(in1)
      file_in2 = trim(path) // trim(in2)
      file_out1 = trim(path) // trim(out1)
      file_out2 = trim(path) // trim(out2)
      open(unit = 11, file = file_in1, status = 'old')
      open(unit = 12, file = file_in2, status = 'old')
      open(unit = 21, file = file_out1, status = 'replace')
      open(unit = 22, file = file_out2, status = 'replace')
      !Skip the header line
      read(11,*)
      read(12,*)
      !Read the first line of data
      read(11,*) input_data(1,:)
      read(12,*) input_data(2,:)
      close(11)
      close(12)
      return
   end subroutine read_input

   subroutine helio2bary(input_data, vb_out)
      implicit none
      real(P), dimension(2,7), intent(in)  :: input_data
      real(P), dimension(2,3), intent(out) :: vb_out
      real(P)                              :: m_sys 
      real(P), dimension(3)                :: vb_sun
      
      m_sys = m_sun + m_nep + m_plu
      vb_sun(:) = - ((m_nep * input_data(1,5:7)) + (m_plu * input_data(2,5:7))) / m_sys
      vb_out(1,:) = input_data(1,5:7) + vb_sun(:)
      vb_out(2,:) = input_data(2,5:7) + vb_sun(:)
      return
   end subroutine helio2bary

   subroutine xv2E0(pl, t, a, e, n, M, E0) 
      implicit none
      real(P), dimension(6), intent(in)    :: pl
      real(P), intent(in)                  :: t
      real(P), intent(out)                 :: a, e, n, M, E0
      real(P)                              :: r_norm, v_norm, h_norm
      
      r_norm = sqrt(pl(1)**2 + pl(2)**2 + pl(3)**2)
      v_norm = sqrt(pl(4)**2 + pl(5)**2 + pl(6)**2)
      h_norm = sqrt(((pl(2) * pl(6)) - (pl(3) * pl(5)))**2 &
         + ((pl(3) * pl(4)) - (pl(1) * pl(6)))**2 &
         + ((pl(1) * pl(5)) - (pl(2) * pl(4)))**2)
      a = ((2.0_P / r_norm) - (v_norm**2 / (G * m_sun)))**(-1)
      e = sqrt(1.0_P - (h_norm**2 / (a * G * m_sun)))
      n = (a**3 / (G * m_sun))**(-0.5_P)
      M = mod(n * t, 2.0_P * PI) 
      E0 = acos((1.0_P - (r_norm / a)) / e)
      if ((E0 + PI) > (2.0_P * PI)) then
         E0 = E0 - PI
      end if 

      return
   end subroutine xv2E0

   subroutine sun_drift(pl1, pl2)
      implicit none
      real(P), dimension(6), intent(inout) :: pl1, pl2
      real(P), dimension(3)                :: v_sys

      v_sys(:) = ((m_nep * pl1(4:6)) + (m_plu * pl2(4:6))) / m_sun
      pl1(1:3) = pl1(1:3) + (v_sys(:) * dt)
      pl2(1:3) = pl2(1:3) + (v_sys(:) * dt)
      
      return
   end subroutine sun_drift

   subroutine pl_kick(pl1, pl2)
      implicit none
      real(P), dimension(6), intent(inout) :: pl1, pl2
      real(P), dimension(3)                :: acc_sys
      real(P)                              :: r, rx, ry, rz
      
      rx = pl1(1) - pl2(1)
      ry = pl1(2) - pl2(2)
      rz = pl1(3) - pl2(3)
      r = sqrt(rx**2 + ry**2 + rz**2)
      acc_sys(:) = (G * m_nep * m_plu) / ((m_nep + m_plu) * r**2)
      pl1(4:6) = pl1(4:6) + (acc_sys(:) * dt)
      pl2(4:6) = pl2(4:6) + (acc_sys(:) * dt)
      return
   end subroutine pl_kick

   subroutine pl_drift(pl, a, e, M, E0)
      implicit none
      real(P), dimension(6), intent(inout) :: pl
      real(P), intent(in)                  :: a, e, M, E0
      integer                              :: j 
      real(P)                              :: func, func_1d, func_2d, func_3d, del1, del2, del3
      real(P), dimension(11)               :: E_val
      
      E_val(1) = E0
      !Drift positions using Danby 
      do j = 1, 10
         !Calculate the first, second, and third derivatives of Kepler's Equation
         func = E_val(j) - (e * sin(E_val(j))) - M
         func_1d = 1.0_P - (e * cos(E_val(j)))
         func_2d = e * sin(E_val(j))
         func_3d = e * cos(E_val(j))
         !Calculate the corresponding delta functions
         del1 = - func / func_1d
         del2 = - func / (func_1d + (0.5_P * del1 * func_2d))
         del3 = - func / (func_1d + (0.5_P * del2 * func_2d) + ((1.0_P/6.0_P) * del2**2.0_P * func_3d))
         !Calculate the eccentric anomaly in radians
         E_val(j+1) = (E_val(j) + del3) 
      end do

      write(*,*) "ORBEL", a, e, M, E0

      !Calculate the new positions
      pl(1) = a * (cos(E_val(11)) - e)
      pl(2) = a * sqrt(1.0_P - e**2.0_P) * sin(E_val(11))
      pl(3) = pl(3) + (0.5_P * dt * pl(6))
      return
   end subroutine pl_drift

   subroutine bary2helio(pl1, pl2, vh_out)
      implicit none
      real(P), dimension(6), intent(in)    :: pl1
      real(P), dimension(6), intent(in)    :: pl2
      real(P), dimension(2,3), intent(out) :: vh_out
      real(P), dimension(3)                :: vh_sun

      vh_sun(:) = - ((m_nep * pl1(4:6)) + (m_plu * pl2(4:6))) / m_sun
      vh_out(1,:) = pl1(4:6) - vh_sun(:)
      vh_out(2,:) = pl2(4:6) - vh_sun(:)
      
      return
   end subroutine bary2helio

end module integrator 

program neptune_pluto
   use integrator 
   implicit none
   integer                              :: dt_max, count_start, count_rate, count_end, i
   character(len = 200)                 :: path, nep_file_in, plu_file_in, nep_file_out, plu_file_out
   real(P), dimension(2,7)              :: input_data
   real(P), dimension(2,3)              :: vb_out, vh_out
   real(P), dimension(6)                :: nep, plu
   real(P)                              :: nep_a, nep_e, nep_n, nep_M, nep_E0, plu_a, plu_e, plu_n, plu_M, plu_E0, t, t_half

   call system_clock(count_start, count_rate)

   path = '/Users/carlislewishard/Documents/Classes/Year_4/Numerical_Dynamics/Module3/&
   module3-carlislewishard/data/'

   nep_file_in = 'id000009-XV.csv'
   plu_file_in = 'id000010-XV.csv'

   nep_file_out = 'id000009-XV-INT.csv'
   plu_file_out = 'id000010-XV-INT.csv'

   call read_input(path, nep_file_in, plu_file_in, nep_file_out, plu_file_out, input_data)

   !Convert the velocities from AU/day -> AU/year 
   input_data(:,5:7) = input_data(:,5:7) * 365.25_P    

   !Write t=0 positions and velocities in original coordinate frame
   write(21,fmt) 0.0, input_data(1,2), input_data(1,3), input_data(1,4), input_data(1,5), input_data(1,6), input_data(1,7)
   write(22,fmt) 0.0, input_data(2,2), input_data(2,3), input_data(2,4), input_data(2,5), input_data(2,6), input_data(2,7) 

   call helio2bary(input_data, vb_out)

   nep(1) = input_data(1,2)
   nep(2) = input_data(1,3)
   nep(3) = input_data(1,4)
   nep(4) = vb_out(1,1)
   nep(5) = vb_out(1,2)
   nep(6) = vb_out(1,3)

   plu(1) = input_data(2,2)
   plu(2) = input_data(2,3)
   plu(3) = input_data(2,4)
   plu(4) = vb_out(2,1)
   plu(5) = vb_out(2,2)
   plu(6) = vb_out(2,3)    

   t = 0.0_P 

   dt_max = int(tmax / dt) 

   do i = 1, dt_max

      t = i * dt
      t_half = 0.5_P * t

      !if (mod(t, 10000.0_P) == 0) then
         write(*,*) '#1 NEP', nep(:)
         !write(*,*) '#1 PLU', plu(:)
      !end if

      call sun_drift(nep, plu)

      !if (mod(t, 10000.0_P) == 0) then
         write(*,*) '#2 NEP', nep(:)
         !write(*,*) '#2 PLU', plu(:)
      !end if

      call pl_kick(nep, plu)

      !if (mod(t, 10000.0_P) == 0) then
         write(*,*) '#3 NEP', nep(:)
         !write(*,*) '#3 PLU', plu(:)
      !end if

      call xv2E0(nep, t, nep_a, nep_e, nep_n, nep_M, nep_E0)
      call xv2E0(plu, t, plu_a, plu_e, plu_n, plu_M, plu_E0)

      !if (mod(t, 10000.0_P) == 0) then
         write(*,*) '#4 NEP', nep(:)
         !write(*,*) '#4 PLU', plu(:)
      !end if

      call pl_drift(nep, nep_a, nep_e, nep_M, nep_E0)
      call pl_drift(plu, plu_a, plu_e, plu_M, plu_E0)

      !if (mod(t, 10000.0_P) == 0) then
         write(*,*) '#5 NEP', nep(:)
         !write(*,*) '#5 PLU', plu(:)
      !end if

      call pl_kick(nep, plu)

      !if (mod(t, 10000.0_P) == 0) then
         write(*,*) '#6 NEP', nep(:)
         !write(*,*) '#6 PLU', plu(:)
      !end if

      call sun_drift(nep, plu)

      !if (mod(t, 10000.0_P) == 0) then
         write(*,*) '#7 NEP', nep(:)
         !write(*,*) '#7 PLU', plu(:)
      !end if

      call bary2helio(nep, plu, vh_out)

      !if (mod(t, 10000.0_P) == 0) then
         write(*,*) '#8 NEP', nep(:)
         !write(*,*) '#8 PLU', plu(:)
      !end if

      !Write new positions and velocities in original coordinate frame
      if (mod(t, 10000.0_P) == 0) then
         write(*,*) "Time = ", t
         write(21,fmt) t, nep(1), nep(2), nep(3), vh_out(1,1), vh_out(1,2), vh_out(1,3)
         write(22,fmt) t, plu(1), plu(2), plu(3), vh_out(2,1), vh_out(2,2), vh_out(2,3)
      end if

      !if (mod(t, 10000.0_P) == 0) then
         write(*,*) '#9 NEP', nep(:)
         !write(*,*) '#9 PLU', plu(:)
      !end if

   end do

   close(21)
   close(22)

   call system_clock(count_end, count_rate)
   write(*,*) "INTEGRATION TIME (S) : ", ((count_end - count_start) / (count_rate * 1.0_P))

end program neptune_pluto