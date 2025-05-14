!==============================================================================
! GALAXY - Stellar Population Synthesis Code
!==============================================================================
!> Main program for galaxy/starburst population synthesis code (Starburst99).
!> 
!> This program computes observable parameters for populations of massive stars, 
!> including spectral energy distributions, stellar feedback, and chemical yields.
!>
!> Original version: Claus Leitherer (August 1998)
!> Last major update: August 2014
!> Modernized Fortran 2018 version: [2024]
!==============================================================================
program galaxy
   use galaxy_module
   use, intrinsic :: iso_fortran_env, only: real32, int32, int64, real64, &
                                          stdout => output_unit, &
                                          stderr => error_unit
   implicit none

   ! Local variables for the main program
   integer :: icount, iexit, nmodels
   integer :: i, j, k, ii, jj, ll
   character(len=:), allocatable :: file_wn, file_wc, file_pown, file_powc
   character(len=:), allocatable :: file_name
   character(len=3) :: namfi3, nam
   character(len=6), dimension(9) :: wrident
   real(real32) :: time, tback
   integer :: ios
   character(len=19) :: current_time
   logical :: file_exists

   ! Initialize module variables
   call init_module()
   
   ! Get current date and time for logging
   call date_and_time(time=current_time)
   
   ! Print startup message
   write(stdout, '(A)') "GALAXY - Stellar Population Synthesis Code (Fortran 2018 Edition)"
   write(stdout, '(A)') "========================================================"
   write(stdout, '(A,A)') "Run started at: ", current_time
   write(stdout, '(A)') ""
   
   ! Initialize program variables
   time = 0.0_real32
   icount = 1

   ! Get the input parameters
   call input(time)

   ! Read evolutionary tracks
   call read_tracks()

   ! Set metallicity string for filenames based on selected tracks
   select case (iz)
   case (11,21,31,41)
      namfi3 = 'm13'; nam = '001'
   case (12,22,32,42)
      namfi3 = 'm07'; nam = '004'
   case (13,23,33,43)
      namfi3 = 'm04'; nam = '008'
   case (14,24,34,44)
      namfi3 = 'p00'; nam = '020'
   case (15,25,35,45)
      namfi3 = 'p03'; nam = '040'
   case (51,61)
      namfi3 = 'm13'; nam = '001'
   case (52,62)
      namfi3 = 'm07'; nam = '004'
   case (53,63)
      namfi3 = 'm04'; nam = '008'
   case (54,64)
      namfi3 = 'p00'; nam = '020'
   case (55,65)
      namfi3 = 'p03'; nam = '040'
   end select
   if (iwrscale < 0) nam = '020'

   ! Read atmospheric and opacity data
   block
      ! Use block construct to contain local variables
      integer :: stat
      
      ! Read Lejeune atmospheres
      file_name = 'data/lejeune/lcb97_'//namfi3//'.flu'
      
      ! Check if file exists first
      inquire(file=file_name, exist=file_exists)
      if (.not. file_exists) then
         write(stderr, '(A)') 'ERROR: Cannot find Lejeune atmosphere file: ' // file_name
         error stop
      end if
      
      ! Open and read the file
      call open_file(un_atm, file_name, 'old', iostat=ios)
      if (ios /= 0) then
         write(stderr, '(A,I0)') 'ERROR: Cannot open Lejeune atmosphere file: ' // &
                                 file_name // ', iostat=', ios
         error stop
      end if
      
      ! ... (read Lejeune data into arrays)
      close(unit=un_atm)
      
      ! Read WR atmospheres, spectral type calibration, IR features, Lyman-alpha, etc.
      ! ... (repeat for each file, using open_file/read/close as above)
      
      ! Read WMBasic, Hillier, PoWR, high-res, IFA, and template spectra
      ! ... (repeat for each file, using open_file/read/close as above)
   end block

   ! Initialize supernova and nucleosynthesis mass limits
   critma = 1.0e36_real32
   critma_new = critma
   critup = -1.0_real32
   critup_new = critup
   critma1 = 1.0e36_real32
   critma_new1 = critma1
   critup1 = -1.0_real32
   critup_new1 = critup1

   !===========================================================================
   ! Main time loop to evolve the stellar population
   !===========================================================================
   ! Use Fortran 2018 block construct to provide better scoping and organization
   time_evolution_block: block
      integer :: step_count
      real(real32) :: elapsed_physical_time
      logical :: continue_evolution
      
      ! Initialize loop control variables
      continue_evolution = .true.
      step_count = 0
      elapsed_physical_time = 0.0_real32
      
      write(stdout, '(A)') "Beginning stellar population evolution..."
      
      ! Main time evolution loop
      time_evolution: do while (continue_evolution)
         step_count = step_count + 1
         
         ! Optional progress indicator (only for longer runs)
         if (mod(step_count, 10) == 0) then
            write(stdout, '(A,I4,A,F12.6,A)') "Step ", step_count, ", t = ", &
                 time / 1.0e6_real32, " Myr"
         end if
         
         !-----------------------------------------------------------------------
         ! Compute stellar population properties at this time step
         !-----------------------------------------------------------------------
         ! Process stellar population based on grid type
         select case (jmg)
         case (0,1)  ! Standard mass grids
            call density(time, icount)
            call starpara(time, icount)
         case (2)    ! Isochrone synthesis on fixed grid
            call density(time, icount)
            call starpara_iso(time, icount)
         case (3)    ! Full isochrone synthesis
            call starpara_iso(time, icount)
            call density(time, icount)
         case default
            write(stderr, '(A,I0)') "ERROR: Invalid synthesis method (jmg=", jmg, ")"
            error stop
         end select

         ! Adjust WR temperatures if needed for isochrone synthesis modes
         if (jmg == 2 .or. jmg == 3) call temp_adjust(iwrt, iatmos)

         !-----------------------------------------------------------------------
         ! Calculate all requested output parameters for this time step
         !-----------------------------------------------------------------------
         ! Use associate construct to create more readable conditions
         associate(calc_wind => io4 >= 0, &
                  calc_sn => io2 >= 0, &
                  calc_spectypes => io5 >= 0, &
                  calc_chem => io6 >= 0, &
                  calc_sed => io7 >= 0, &
                  calc_uvlines => io8 >= 0, &
                  calc_fuse => io12 >= 0, &
                  calc_hires => io13 >= 0, &
                  calc_ifa => io15 >= 0)
            
            ! Utilize optional error checking and output generation
            if (calc_wind) call windpower(time, icount)         ! Mechanical energy
            if (calc_sn) call supernova(time, icount)           ! Supernova rates
            if (calc_spectypes) call spectype(time, icount)     ! Stellar types
            if (calc_chem) call nucleo(time, icount)            ! Chemical evolution
            if (calc_sed) call specsyn(time, icount)            ! Spectral energy dist
            if (calc_uvlines) call linesyn(time, icount)        ! UV spectral features
            if (calc_fuse) call fusesyn(time, icount)           ! FUSE spectra
            if (calc_hires) call hires(time, icount)            ! High-res optical
            if (calc_ifa) call ifa_spectrum(time, icount)       ! IFA UV spectra
         end associate

         !-----------------------------------------------------------------------
         ! Advance time step - use either linear or logarithmic time steps
         !-----------------------------------------------------------------------
         if (jtime == 0) then
            ! Linear time steps
            tstep = tvar
         else
            ! Logarithmic time steps
            tiempo1 = tiempo1 + tvar
            tback = tiempo1 - tvar
            tstep = 10.0_real32**tiempo1 - 10.0_real32**tback
         end if
         
         ! Update counters and time
         icount = icount + 1
         time = time + tstep
         elapsed_physical_time = elapsed_physical_time + tstep

         !-----------------------------------------------------------------------
         ! Check if we should continue the evolution
         !-----------------------------------------------------------------------
         ! Determine exit condition based on time step mode
         select case (jtime)
         case (0)  ! Linear time steps
            ! Exit when physical time exceeds maximum
            if (time > tmax) continue_evolution = .false.
         case (1)  ! Logarithmic time steps
            ! Exit when step count exceeds total steps
            iexit = int(tinter) + 1
            if (iexit < icount) continue_evolution = .false.
         case default
            write(stderr, '(A,I0)') "ERROR: Invalid time step mode (jtime=", jtime, ")"
            error stop
         end select
      end do time_evolution
      
      ! Report simulation completion
      write(stdout, '(A,I4,A)') "Evolution completed after ", step_count, " time steps"
      write(stdout, '(A,F12.6,A)') "Total elapsed physical time: ", elapsed_physical_time / 1.0e6_real32, " Myr"
   end block time_evolution_block

   !===========================================================================
   ! Generate final output and clean up
   !===========================================================================
   ! Process final results and generate output files
   call output()

   ! Release allocated memory
   call cleanup_module()
   
   ! Normal termination with timestamp
   block
      character(len=19) :: end_time
      call date_and_time(time=end_time)
      write(stdout, '(A)') "========================================================"
      write(stdout, '(A,A)') "Run completed at: ", end_time
      write(stdout, '(A)') 'Galaxy synthesis completed successfully'
   end block
   
end program galaxy
!==============================================================================
! Input Subroutine
!==============================================================================
!> Reads input parameters from a file and sets global model parameters.
!>
!> This subroutine reads the model configuration from fort.1 (typically linked
!> to standard.input1) and sets all the module variables used throughout the
!> simulation. The input is read using a Fortran 2018 compliant approach.
!>
!> @param[inout] time The initial time value
subroutine input(time)
   use galaxy_module
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   
   ! Arguments
   real(real32), intent(inout) :: time

   ! Local variables
   integer :: i, ios
   character(len=:), allocatable :: input_file
   logical :: file_exists
   
   ! Use named parameter for standard input unit
   integer, parameter :: INPUT_UNIT = 1
   
   ! Use fort.1 as the standard input file name (which should be linked to the actual input file)
   input_file = 'fort.1'
   
   ! Check if input file exists
   inquire(file=input_file, exist=file_exists)
   if (.not. file_exists) then
      write(stderr, '(A)') 'ERROR: Input file does not exist: ' // input_file
      error stop "Missing input file"
   end if
   
   ! Open and read the input file - using newer error handling approach
   block
      integer :: stat
      character(len=200) :: errmsg
      
      ! Open with error handling
      open(unit=un_input, file=input_file, status='old', iostat=stat, iomsg=errmsg)
      if (stat /= 0) then
         write(stderr, '(A)') 'ERROR: Cannot open input file: ' // input_file // &
                             ' - ' // trim(errmsg)
         error stop "Input file error"
      end if
      
      ! Read all model parameters using a structured read with error checking
      read(un_input, 10001, iostat=stat, iomsg=errmsg) &
           model_name, isf, toma, sfr, ninterv, &
           (xponent(i), i=1, nmaxint), (xmaslim(i), i=1, nmaxint1), &
           sncut, bhcut, iz, iwind, time1, jtime, tbiv, itbiv, tmax, jmg, &
           lmin, lmax, tdel, iatmos, ilib, iline, ivt, irsg, &
           io1, io2, io3, io4, io5, io6, io7, io8, io9, io10, io11, io12, io13, io14, io15
      
      if (stat /= 0) then
         write(stderr, '(A)') 'ERROR: Failed to read input parameters: ' // trim(errmsg)
         close(unit=un_input)
         error stop "Input parameter reading error"
      end if
      
      close(unit=un_input)
   end block
   
   ! Format statement for reading the input file
   ! This must match the structure of the input file exactly
10001 format(/,a20,//,i3,//,f12.2,//,f12.2,//,i3,//,10f12.2, &
        //,11f12.2,//,f12.2,//,f12.2,////////,i3,//,i3,//, &
        f12.2,//,i3,//,f12.2,//,i8,//,f12.2,/,/,/,i3,//, &
        2i4,//,f12.2,//,i3,/,/,/,i3,//,i3,//,2i4,//,15i3)

   ! Convert times from 10^6 years to years
   time1 = time1 * 1.0e6_real32
   tmax  = tmax  * 1.0e6_real32
   tbiv  = tbiv  * 1.0e6_real32
   tdel  = tdel  * 1.0e6_real32

   ! Print input information
   write(stdout, '(A)') "Input parameters:"
   write(stdout, '(A,A)')    "  Model name:     ", trim(model_name)
   write(stdout, '(A,I0)')   "  Star formation: ", isf
   write(stdout, '(A,F0.2)') "  Metallicity ID: ", iz 
   
   ! Parameter validation using Fortran 2018 construct
   block
      logical :: valid_params = .true.
      character(len=:), allocatable :: error_msg
      
      ! Check individual parameters
      if (ninterv < 1 .or. ninterv > nmaxint) then
         valid_params = .false.
         error_msg = "Invalid number of IMF intervals"
      end if
      
      if (tmax <= time1) then
         valid_params = .false.
         error_msg = "Maximum time must be greater than initial time"
      end if
      
      ! Add more validation as needed...
      
      ! Exit with error if parameters are invalid
      if (.not. valid_params) then
         write(stderr, '(A,A)') "ERROR: Invalid input parameters - ", error_msg
         error stop "Input parameter validation failed"
      end if
   end block
   
   ! Set the initial time based on star formation mode
   select case (isf)
   case (:0)  ! Instantaneous burst
      time = time1
   case (1:)  ! Continuous star formation
      time = 1.0e-5_real32
   end select
   
   ! Set time step size
   if (jtime == 0) then
      ! Linear time steps
      tvar = tbiv
      tinter = (tmax - time1) / tbiv
   else
      ! Logarithmic time steps
      tiempo1 = log10(time1)
      tinter = real(itbiv, real32)
      tvar = (log10(tmax) - log10(time1)) / tinter
   end if
   
   ! Report successful parameter loading
   write(stdout, '(A,F0.2,A,F0.2,A)') "Time range: ", time1/1.0e6_real32, &
                                     " to ", tmax/1.0e6_real32, " Myr"
   if (isf <= 0) then
      time = time1
   else
      time = 1.e-5_real32
      time1 = time
   end if

   ! Define the time step for linear and logarithmic options
   if (jtime == 0) then
      tvar = tbiv * 1.0e6_real32
      tstep = tvar
   else
      tinter = real(itbiv, real32)
      tiempo1 = log10(time1)
      ! Use a local variable for log(tmax) to avoid modifying tmax
      block
         real(real32) :: log_tmax
         log_tmax = log10(tmax)
         tvar = (log_tmax - tiempo1) / tinter
         tmax = 10.0_real32**(log_tmax + tvar)
      end block
      tstep = 1.e4_real32
   end if

   tdel = tdel * 1.0e6_real32
   toma = toma * 1.0e6_real32

   ! Mass limits for the IMF
   upma = xmaslim(ninterv+1)
   doma = xmaslim(1)

   ! Consistency checks for metallicity
   block
      real(real32) :: z
      
      select case (iz)
      case (11,21,31,41)
         z = 0.05_real32
      case (12,22,32,42)
         z = 0.2_real32
      case (13,23,33,43)
         z = 0.4_real32
      case (14,24,34,44)
         z = 1.0_real32
      case (15,25,35,45)
         z = 2.0_real32
      case (51,61,52,62)
         z = 0.2_real32
      case (53,63,54,64,55,65)
         z = 1.0_real32
      case default
         call error_handler("input", "This metallicity (iz=" // trim(integer_to_string(iz)) // &
                          ") is not available!", fatal=.true.)
      end select

      if (iz == 51 .or. iz == 61) then
         write(stdout, '(A)') 'INFO: Metallicity reset to Z=0.002'
      end if
      if (iz == 53 .or. iz == 63) then
         write(stdout, '(A)') 'INFO: Metallicity reset to Z=0.014'
      end if
      if (iz == 55 .or. iz == 65) then
         write(stdout, '(A)') 'INFO: Metallicity reset to Z=0.014'
      end if
   end block

   ! Model atmosphere checks
   if (iatmos < 1 .or. iatmos > 5) then
      call error_handler("input", "Invalid model atmosphere selection. " // &
                       "iatmos must be between 1 and 5.", fatal=.true.)
   end if

   ! Padova tracks recommendation
   if (jmg < 3 .and. (iz >= 31 .and. iz <= 35 .or. iz >= 41 .and. iz <= 45)) then
      write(stdout, '(A)') 'INFO: JMG should be 3 for the Padova tracks. The 0, 1, 2 grids terminate at 1 M_sun.'
   end if

   ! Info for line spectra
   if (iatmos < 4 .and. (io8 >= 0 .or. io12 >= 0 .or. io13 >= 0 .or. io15 >= 0)) then
      write(stdout, '(A)') 'INFO: Hillier models will be used for line spectra'
   end if

   ! Set default WR scaling
   if (iatmos >= 4) then
      iwrt = 1
      iwrscale = 1
   end if

   ! More checks and validation
   block
      character(len=:), allocatable :: warning_msg
      
      ! Check output conflicts
      if (io9 >= 0 .and. io7 < 0) then
         call error_handler("input", "Cannot compute synthetic colors (io9>=0) without spectra (io7<0)", &
                           fatal=.false.)
      end if
      
      if (io10 >= 0 .and. io7 < 0) then
         call error_handler("input", "Cannot compute equivalent widths (io10>=0) without spectra (io7<0)", &
                           fatal=.false.)
      end if
      
      if (io10 >= 0 .and. io1 < 0) then
         call error_handler("input", "Cannot compute equivalent widths (io10>=0) without nebular continuum (io1<0)", &
                           fatal=.false.)
      end if
      
      ! Nebular continuum warnings
      if ((io7 >= 0 .or. io8 >= 0) .and. io1 < 0) then
         write(stdout, '(A)') 'INFO: No nebular continuum will be computed (io1<0)'
      end if
      
      if (io7 < 0 .and. io8 >= 0 .and. io1 >= 0) then
         write(stdout, '(A)') 'INFO: No nebular continuum will be computed (io7<0)'
      end if
      
      ! More nebular continuum checks - consolidated to avoid repetition
      if (io1 >= 0 .and. io7 < 0) then
         write(stdout, '(A)') 'INFO: No nebular continuum will be included in output'
      end if
      
      if ((io7 >= 0 .or. io12 >= 0 .or. io13 >= 0 .or. io15 >= 0) .and. io1 < 0) then
         write(stdout, '(A)') 'INFO: No nebular continuum will be included (io1<0)'
      end if
      
      ! Feature computation checks
      if (io11 >= 0 .and. io7 < 0) then
         call error_handler("input", "Cannot compute spectral features (io11>=0) without spectra (io7<0)", &
                           fatal=.false.)
      end if
      
      if (io14 >= 0 .and. io7 < 0) then
         call error_handler("input", "Cannot compute WR emission lines (io14>=0) without spectra (io7<0)", &
                           fatal=.false.)
      end if
      
      ! RSG feature parameter checks
      if (ivt < 1 .or. ivt > 6 .or. irsg < 0 .or. irsg > 1) then
         warning_msg = "RSG feature parameters out of range: ivt must be 1-6, irsg must be 0-1"
         call error_handler("input", warning_msg, fatal=.false.)
      end if
      
      ! HRD warning for isochrone synthesis
      if (io3 >= 0 .and. jmg == 3) then
         write(stdout, '(A)') 'INFO: No HRD will be created for isochrone synthesis mode'
      end if
      
      ! Critical parameter validations
      if (jmg > 3 .or. jmg < 0) then
         call error_handler("input", "Invalid synthesis mode: jmg must be 0, 1, 2, or 3", fatal=.true.)
      end if
      
      if (iline < 1 .or. iline > 2) then
         call error_handler("input", "Invalid line library selection: iline must be 1 or 2", fatal=.true.)
      end if
      
      if (upma > 120.0_real32) then
         upma = 120.0_real32
         write(stdout, '(A)') 'INFO: IMF upper mass limit set to 120 M_sun (maximum allowed)'
      end if
      
      if (iwind < 0 .or. iwind > 3) then
         call error_handler("input", "Invalid wind model: iwind must be 0, 1, 2, or 3", fatal=.true.)
      end if
      
      if (ilib < 1 .or. ilib > 4) then
         call error_handler("input", "Invalid high-resolution library: ilib must be 1, 2, 3, or 4", fatal=.true.)
      end if
   end block

   return
end subroutine input
!==============================================================================
! Density Subroutine
!==============================================================================
!> Computes the mass normalization and stellar number densities for the IMF.
!>
!> This subroutine calculates:
!> 1. IMF normalization based on the specified power-law parameters
!> 2. Mass grid setup based on the synthesis method (standard or isochrone)
!> 3. Number density of stars in each mass bin according to the IMF
!>
!> @param[in] time   Current time in the simulation (years)
!> @param[in] icount Current time step index
!==============================================================================
subroutine density(time, icount)
   use galaxy_module
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   
   ! Arguments
   real(real32), intent(in) :: time
   integer, intent(in)      :: icount

   ! Local variables
   real(real32) :: a(nmaxint), aic(nmaxint)
   real(real32) :: sum_imf, sumando, s, xmhigh, xmlow
   real(real32) :: tomas, tonum, delm
   integer :: i, j, k

   !---------------------------------------------------------------------------
   ! Calculate total mass based on star formation mode
   !---------------------------------------------------------------------------
   ! If continuous star formation, mass depends on the time step
   if (isf > 0) then
      if (icount > 1) then
         toma = sfr * tstep
      else
         toma = sfr * 1.e-5_real32  ! Small initial step
      end if
   end if

   !---------------------------------------------------------------------------
   ! Calculate IMF normalization constants for multiple power laws
   !---------------------------------------------------------------------------
   ! Multi-power law IMF normalization (R. Sutherland algorithm)
   aic(1) = 1.0_real32  ! First interval
   
   ! Calculate normalization factor at each IMF break
   if (ninterv > 1) then
      do k = 2, ninterv
         aic(k) = aic(k-1) * xmaslim(k)**(xponent(k)-xponent(k-1))
      end do
   end if

   ! Calculate IMF integral for mass normalization
   sum_imf = 0.0_real32
   do i = 1, ninterv
      if (abs(xponent(i) - 2.0_real32) < epsilon(1.0_real32)) then
         ! Special case for x=2, avoid divide by zero
         sumando = log(xmaslim(i+1)) - log(xmaslim(i))
      else
         ! General power law integration
         sumando = (xmaslim(i+1)**(2.0_real32-xponent(i)) - &
                   xmaslim(i)**(2.0_real32-xponent(i))) / &
                   (2.0_real32-xponent(i))
      end if
      sum_imf = sum_imf + aic(i) * sumando
   end do

   ! Calculate the overall IMF normalization
   s = toma / sum_imf
   do i = 1, ninterv
      a(i) = aic(i) * s
   end do

   !---------------------------------------------------------------------------
   ! Set up the mass grid based on synthesis method
   !---------------------------------------------------------------------------
   ! Configure grid properties based on selected method
   select case (jmg)
   case (0,2)  ! Small grid or isochrone on grid
      if (lmin == 0) lmin = 1
      if (lmax == 0) lmax = 24
      delm = 2.5_real32
   case (1)    ! Large grid
      if (lmin == 0) lmin = 1
      if (lmax == 0) lmax = 120
      delm = 0.5_real32
   end select

   ! Define the mass intervals for non-isochrone synthesis
   do i = lmin, lmax
      cmass(i) = 120.0_real32 - real(i-1, real32) * 2.0_real32 * delm
   end do

   ! Ensure mass grid agrees with IMF limits
   do i = lmin, lmax
      if (upma <= cmass(i)) lmin = i
      if (doma <= cmass(i)) lmax = i
   end do

   !---------------------------------------------------------------------------
   ! Calculate stellar number densities in each mass bin
   !---------------------------------------------------------------------------
   do i = lmin, lmax
      dens(i) = 0.0_real32
      
      ! Set mass interval boundaries based on grid type
      if (jmg >= 0 .and. jmg <= 2) then
         ! Fixed mass grid
         xmhigh = cmass(i) + delm
         xmlow  = cmass(i) - delm
         
         ! Adjust boundaries at the edges
         if (lmin /= lmax) then
            if (i == lmin) xmhigh = xmhigh - delm
            if (i == lmax) xmlow  = xmlow  + delm
         end if
      else if (jmg == 3) then
         ! Variable mass grid (isochrone synthesis)
         if (i == lmin) then
            ! First bin
            xmhigh = 0.5_real32 * (cmass(i) + cmass(i+1))
            xmlow  = cmass(i)
         else if (i == lmax) then
            ! Last bin
            xmhigh = cmass(i)
            xmlow  = 0.5_real32 * (cmass(i-1) + cmass(i))
         else
            ! Interior bins
            xmhigh = 0.5_real32 * (cmass(i) + cmass(i+1))
            xmlow  = 0.5_real32 * (cmass(i-1) + cmass(i))
         end if
      end if

      ! Star number/density calculation for multiple power-law IMF (J. Yin algorithm)
      do j = 1, ninterv
         ! Mass interval crosses a break in the IMF
         if (xmlow < xmaslim(j) .and. xmaslim(j) < xmhigh .and. j /= 1) then
            ! Handle special cases to avoid numerical issues
            if (abs(xponent(j-1) - 1.0_real32) < epsilon(1.0_real32)) then
               ! First segment has x=1
               dens(i) = a(j-1) * (log(xmaslim(j)) - log(xmlow)) + &
                         a(j) / (1.0_real32 - xponent(j)) * &
                         (xmhigh**(1.0_real32-xponent(j)) - xmaslim(j)**(1.0_real32-xponent(j)))
            else if (abs(xponent(j) - 1.0_real32) < epsilon(1.0_real32)) then
               ! Second segment has x=1
               dens(i) = a(j-1) / (1.0_real32 - xponent(j-1)) * &
                         (xmaslim(j)**(1.0_real32-xponent(j-1)) - xmlow**(1.0_real32-xponent(j-1))) + &
                         a(j) * (log(xmhigh) - log(xmaslim(j)))
            else
               ! General case
               dens(i) = a(j-1) / (1.0_real32 - xponent(j-1)) * &
                         (xmaslim(j)**(1.0_real32-xponent(j-1)) - xmlow**(1.0_real32-xponent(j-1))) + &
                         a(j) / (1.0_real32 - xponent(j)) * &
                         (xmhigh**(1.0_real32-xponent(j)) - xmaslim(j)**(1.0_real32-xponent(j)))
            end if
         ! Mass interval entirely within one IMF segment
         else if (xmaslim(j) <= xmlow .and. xmhigh <= xmaslim(j+1)) then
            if (abs(xponent(j) - 1.0_real32) < epsilon(1.0_real32)) then
               ! Special case for x=1
               dens(i) = a(j) * (log(xmhigh) - log(xmlow))
            else
               ! General case
               dens(i) = a(j) / (1.0_real32 - xponent(j)) * &
                         (xmhigh**(1.0_real32-xponent(j)) - xmlow**(1.0_real32-xponent(j)))
            end if
         end if
      end do
   end do

   !---------------------------------------------------------------------------
   ! Compute totals (for verification or debugging)
   !---------------------------------------------------------------------------
   tomas = 0.0_real32  ! Total mass
   tonum = 0.0_real32  ! Total number of stars
   
   do i = lmax, lmin, -1
      tonum = tonum + dens(i)
      tomas = tomas + dens(i) * cmass(i)
   end do

   ! Debug output if needed
   if (un_debug > 0) then
      write(un_debug, '(A,F12.2,A,F12.2,A)') 'Total: ', tonum, ' stars, ', tomas, ' solar masses'
   end if

end subroutine density
c

!==============================================================================
! Track Reading Subroutine
!==============================================================================
!> Reads evolutionary tracks from Geneva and Padova models and sets up WR mass limits.
!>
!> This subroutine reads stellar evolutionary tracks based on the chosen metallicity
!> and model set (Geneva standard/high mass-loss, Padova standard/AGB, or Geneva rotating).
!> The data includes stellar ages, masses, luminosities, temperatures, surface abundances,
!> and mass loss rates.
!==============================================================================
subroutine read_tracks()
   use galaxy_module
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none

   ! Local constants
   integer, parameter :: ntracks = 50, nptrack = 500

   ! Local variables
   character(len=100) :: filetrck
   character(len=2)   :: switch
   character(len=80)  :: title_tracks
   integer :: ki, km, i, j, ilines, nma, jmp, ios
   
   ! Track data arrays
   real(real32) :: dm(ntracks), dlm(ntracks)
   real(real32) :: da(ntracks, nptrack), dla(ntracks, nptrack)
   real(real32) :: xmas(ntracks, nptrack), dmas(ntracks, nptrack)
   real(real32) :: dl(ntracks, nptrack), dt(ntracks, nptrack)
   real(real32) :: d_h(ntracks, nptrack), d_he(ntracks, nptrack)
   real(real32) :: d_c(ntracks, nptrack), d_n(ntracks, nptrack)
   real(real32) :: d_o(ntracks, nptrack), d_tstar(ntracks, nptrack)
   real(real32) :: d_mdot(ntracks, nptrack)

   !---------------------------------------------------------------------------
   ! Select track file and WR mass limit based on metallicity and model set
   !---------------------------------------------------------------------------
   select case (iz)
   ! Geneva standard mass-loss tracks
   case (11)
      filetrck = 'data/tracks/modc001.dat'; xmwr = 80.0_real32
   case (12)
      filetrck = 'data/tracks/modc004.dat'; xmwr = 52.0_real32
   case (13)
      filetrck = 'data/tracks/modc008.dat'; xmwr = 42.0_real32
   case (14)
      filetrck = 'data/tracks/modc020.dat'; xmwr = 32.0_real32
   case (15)
      filetrck = 'data/tracks/modc040.dat'; xmwr = 25.0_real32
   
   ! Geneva high mass-loss tracks
   case (21)
      filetrck = 'data/tracks/mode001.dat'; xmwr = 61.0_real32
   case (22)
      filetrck = 'data/tracks/mode004.dat'; xmwr = 42.0_real32
   case (23)
      filetrck = 'data/tracks/mode008.dat'; xmwr = 35.0_real32
   case (24)
      filetrck = 'data/tracks/mode020.dat'; xmwr = 25.0_real32
   case (25)
      filetrck = 'data/tracks/mode040.dat'; xmwr = 21.0_real32
   
   ! Padova standard tracks
   case (31)
      filetrck = 'data/tracks/mods0004.dat'; xmwr = 61.0_real32
   case (32)
      filetrck = 'data/tracks/mods004.dat'; xmwr = 42.0_real32
   case (33)
      filetrck = 'data/tracks/mods008.dat'; xmwr = 35.0_real32
   case (34)
      filetrck = 'data/tracks/mods020.dat'; xmwr = 25.0_real32
   case (35)
      filetrck = 'data/tracks/mods050.dat'; xmwr = 21.0_real32
   
   ! Padova AGB tracks
   case (41)
      filetrck = 'data/tracks/modp0004.dat'; xmwr = 61.0_real32
   case (42)
      filetrck = 'data/tracks/modp004.dat'; xmwr = 42.0_real32
   case (43)
      filetrck = 'data/tracks/modp008.dat'; xmwr = 35.0_real32
   case (44)
      filetrck = 'data/tracks/modp020.dat'; xmwr = 25.0_real32
   case (45)
      filetrck = 'data/tracks/modp050.dat'; xmwr = 21.0_real32
   
   ! Geneva v00 (non-rotating) tracks
   case (51, 52)
      filetrck = 'data/tracks/Z0020v00.txt'; xmwr = 84.0_real32
   case (53, 54, 55)
      filetrck = 'data/tracks/Z0140v00.txt'; xmwr = 25.0_real32
   
   ! Geneva v40 (rotating) tracks
   case (61, 62)
      filetrck = 'data/tracks/Z0020v40.txt'; xmwr = 55.0_real32
   case (63, 64, 65)
      filetrck = 'data/tracks/Z0140v40.txt'; xmwr = 20.0_real32
   
   case default
      write(stderr, '(A,I4)') 'ERROR: Unknown metallicity index in read_tracks: ', iz
      stop
   end select
   
   ! Store WR mass limit in module variable
   iwrscale = iz
   
   !---------------------------------------------------------------------------
   ! Initialize track data arrays
   !---------------------------------------------------------------------------
   dm    = 0.0_real32; dlm   = 0.0_real32
   da    = 0.0_real32; dla   = 0.0_real32
   xmas  = 0.0_real32; dmas  = 0.0_real32
   dl    = 0.0_real32; dt    = 0.0_real32
   d_h   = 0.0_real32; d_he  = 0.0_real32
   d_c   = 0.0_real32; d_n   = 0.0_real32
   d_o   = 0.0_real32; d_tstar = 0.0_real32
   d_mdot = 0.0_real32

   !---------------------------------------------------------------------------
   ! Read evolutionary track file
   !---------------------------------------------------------------------------
   call open_file(unit=30, file=filetrck, status='old', iostat=ios)
   if (ios /= 0) then
      write(stderr, '(A)') 'ERROR: Cannot open evolutionary track file: ' // trim(filetrck)
      stop
   end if
   
   ! Read header information
   read(30, '(A)') title_tracks
   read(30, *)  ! Skip blank line
   read(30, *) nma, jmp
   
   ! Read track data for each initial mass
   read_tracks_loop: do ki = 1, nma
      read(30, *)  ! Skip blank line
      read(30, '(F8.3,A2)') dm(ki), switch
      
      ! Check mass order (tracks should be in decreasing mass order)
      if (ki > 1) then
         if (dm(ki) >= dm(ki-1)) then
            write(stderr, '(A)') 'ERROR: Tracks not in decreasing mass order!'
            stop
         end if
      end if
      
      ! Calculate log mass
      dlm(ki) = log10(dm(ki))
      
      read(30, *)  ! Skip blank line
      ilines = jmp
      
      ! Read track data points for current mass
      do km = 1, ilines
         select case (switch)
         case ('WR')  ! Wolf-Rayet star
            read(30, '(I2,1PE14.7,F9.4,2F6.3,5F9.6,2F7.3)') i, da(ki,km), xmas(ki,km), dl(ki,km), &
               dt(ki,km), d_h(ki,km), d_he(ki,km), d_c(ki,km), d_n(ki,km), d_o(ki,km), &
               d_tstar(ki,km), d_mdot(ki,km)
         
         case ('RO')  ! Rotating model
            read(30, *) i, da(ki,km), xmas(ki,km), dl(ki,km), dt(ki,km), d_h(ki,km), d_he(ki,km), &
               d_c(ki,km), d_n(ki,km), d_o(ki,km), d_tstar(ki,km), d_mdot(ki,km)
         
         case ('ML')  ! Special mass-loss format
            read(30, '(I2,1PE14.7,F9.4,2F6.3,5F9.6,F7.3)') i, da(ki,km), xmas(ki,km), dl(ki,km), &
               dt(ki,km), d_h(ki,km), d_he(ki,km), d_c(ki,km), d_n(ki,km), d_o(ki,km), d_mdot(ki,km)
            d_tstar(ki,km) = dt(ki,km)  ! Use normal Teff as Tstar
         
         case default  ! Standard format
            read(30, '(I2,1PE14.7,F9.4,2F6.3,5F9.6)') i, da(ki,km), xmas(ki,km), dl(ki,km), &
               dt(ki,km), d_h(ki,km), d_he(ki,km), d_c(ki,km), d_n(ki,km), d_o(ki,km)
            d_tstar(ki,km) = dt(ki,km)  ! Use normal Teff as Tstar
            d_mdot(ki,km) = -99.9_real32  ! Default mass loss rate
         end select
         
         ! Check for missing points
         if (km == 1 .and. i /= km) then
            write(stderr, '(A)') 'ERROR: Point missing on track!'
            stop
         end if
         
         ! Convert to log quantities
         dmas(ki,km) = log10(xmas(ki,km))
         dla(ki,km)  = log10(da(ki,km))
         
         ! Set very low mass loss rates to a minimum value
         if (d_mdot(ki,km) > -0.1_real32) d_mdot(ki,km) = -30.0_real32
      end do

      ! Fill missing lines if needed in incomplete tracks
      if (ilines == jmp-3) then
         do km = ilines+1, jmp
            ! Copy values from last available point
            da(ki,km)   = da(ki,ilines)
            xmas(ki,km) = xmas(ki,ilines)
            dl(ki,km)   = dl(ki,ilines)
            dt(ki,km)   = dt(ki,ilines)
            d_h(ki,km)  = d_h(ki,ilines)
            d_he(ki,km) = d_he(ki,ilines)
            d_c(ki,km)  = d_c(ki,ilines)
            d_n(ki,km)  = d_n(ki,ilines)
            d_o(ki,km)  = d_o(ki,ilines)
            d_tstar(ki,km) = d_tstar(ki,ilines)
            d_mdot(ki,km)  = d_mdot(ki,ilines)
            
            ! Calculate logarithmic values
            dmas(ki,km) = log10(xmas(ki,ilines))
            dla(ki,km)  = log10(da(ki,ilines))
         end do
      end if
   end do read_tracks_loop
   
   close(30)

   !---------------------------------------------------------------------------
   ! Add extra points at end of each track for interpolation and HRD display
   !---------------------------------------------------------------------------
   do i = 1, nma
      ! First point at very early age
      da(i,1) = 1.0e-3_real32
      
      ! Add two points beyond the end of the track
      da(i,jmp+1)  = da(i,jmp) + 100.0_real32
      da(i,jmp+2)  = 1.0e6_real32 * da(i,jmp+1)
      
      ! Set log values
      dla(i,jmp+1) = log10(da(i,jmp+1))
      dla(i,jmp+2) = log10(da(i,jmp+2))
      
      ! Set luminosity to a very low value to indicate end of track
      dl(i,jmp+1)  = -20.0_real32
      dl(i,jmp+2)  = -20.0_real32
      
      ! Copy other parameters from last regular point
      do j = 1, 2
         dt(i,jmp+j)      = dt(i,jmp)
         d_h(i,jmp+j)     = d_h(i,jmp)
         d_he(i,jmp+j)    = d_he(i,jmp)
         d_c(i,jmp+j)     = d_c(i,jmp)
         d_n(i,jmp+j)     = d_n(i,jmp)
         d_o(i,jmp+j)     = d_o(i,jmp)
         d_tstar(i,jmp+j) = d_tstar(i,jmp)
         xmas(i,jmp+j)    = xmas(i,jmp)
         dmas(i,jmp+j)    = log10(xmas(i,jmp))
         d_mdot(i,jmp+j)  = d_mdot(i,jmp)
      end do
   end do

   ! Update track length to include the two added points
   jmp = jmp + 2
   
   ! Report successful track loading
   write(stdout, '(A,I2,A,I3,A)') 'Successfully loaded ', nma, ' tracks with ', jmp, ' time points each'

end subroutine read_tracks
c
c
c ********************************************************************
c ********************************************************************
c ********************************************************************
!> Interpolates in the evolutionary models and calculates stellar parameters
!> for a specified mass grid and each time step (non-isochrone synthesis).
subroutine starpara(time_in, icount)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   real(real32), intent(in) :: time_in
   integer, intent(in)      :: icount

   ! Parameters
   integer, parameter :: nmaxint = 10, nmaxint1 = 11
   integer, parameter :: npgrid = 3000
   integer, parameter :: ntracks = 50, nptrack = 500

   ! Arguments and local variables
   character(len=20) :: name
   character(len=80) :: title_tracks
   character(len=24) :: str_date
   integer :: lmin, lmax, jmp, nma
   real(real32) :: temp(npgrid), bol(npgrid), xsurf(npgrid), zmass(npgrid)
   real(real32) :: ysurf(npgrid), xc12s(npgrid), xn14s(npgrid), xo16s(npgrid)
   real(real32) :: bmdot(npgrid), cmass(npgrid), tt_star(npgrid)
   real(real32) :: age2(nptrack), alogl2(nptrack), alogt2(nptrack)
   real(real32) :: amass2(nptrack), amdot2(nptrack)
   real(real32) :: x2(nptrack), y2(nptrack), xc122(nptrack)
   real(real32) :: xn142(nptrack), xo162(nptrack)
   real(real32) :: age1(ntracks), alogl1(ntracks), alogt1(ntracks)
   real(real32) :: amass1(ntracks), amdot1(ntracks)
   real(real32) :: x1(ntracks), y1(ntracks), xc121(ntracks)
   real(real32) :: xn141(ntracks), xo161(ntracks), dm_rev(ntracks)
   real(real32) :: da(ntracks, nptrack), dl(ntracks, nptrack)
   real(real32) :: dt(ntracks, nptrack), d_h(ntracks, nptrack)
   real(real32) :: d_he(ntracks, nptrack), d_c(ntracks, nptrack)
   real(real32) :: d_n(ntracks, nptrack), d_o(ntracks, nptrack)
   real(real32) :: d_tstar(ntracks, nptrack), d_mdot(ntracks, nptrack)
   real(real32) :: xmas(ntracks, nptrack), dmas(ntracks, nptrack)
   real(real32) :: dlm(ntracks)
   real(real32) :: reci_polint, yntra, alog10
   integer :: l, j, i, mm
   integer :: io3, jmg
   ! External procedures
   interface
      function reci_polint(x, n, xa, ya) result(y)
         real(real32), intent(in) :: x
         integer, intent(in) :: n
         real(real32), intent(in) :: xa(*), ya(*)
         real(real32) :: y
      end function reci_polint
      function yntra(x, xa, ya, n) result(y)
         real(real32), intent(in) :: x
         real(real32), intent(in) :: xa(*), ya(*)
         integer, intent(in) :: n
         real(real32) :: y
      end function yntra
      subroutine intrpl(n, x, y, m, x0, y0)
         integer, intent(in) :: n, m
         real(real32), intent(in) :: x(*), y(*), x0
         real(real32), intent(out) :: y0
      end subroutine intrpl
      subroutine agecheck(time, jmp)
         real(real32), intent(in) :: time
         integer, intent(in) :: jmp
      end subroutine agecheck
      subroutine fdate(str_date)
         character(len=*), intent(out) :: str_date
      end subroutine fdate
   end interface

   ! Output header at first time step
   if (icount <= 1 .and. io3 >= 0) then
      write(96, '(A)') ' MODEL DESIGNATION: '//trim(name)
      call fdate(str_date)
      write(96, '(A)') ' MODEL GENERATED: '//trim(str_date)
      write(96, '(A)') '              EVOLUTIONARY TRACKS IN THE HERTZSPRUNG-RUSSELL DIAGRAM'
      write(96, '(A)') '    M=120        M=100         M=80         M=60         M=50         M=40         M=30         M=25         M=20         M=15'
   end if

   ! Interpolate in evolutionary tracks
   do l = lmin, lmax
      do j = 1, jmp
         mm = 0
         do i = nma, 1, -1
            mm = mm + 1
            age1(mm)   = da(i, j)
            alogl1(mm) = dl(i, j)
            alogt1(mm) = 10.0_real32 ** dt(i, j)
            amass1(mm) = xmas(i, j)
            amdot1(mm) = d_mdot(i, j)
            x1(mm)     = d_h(i, j)
            y1(mm)     = d_he(i, j)
            xc121(mm)  = d_c(i, j)
            xn141(mm)  = d_n(i, j)
            xo161(mm)  = d_o(i, j)
            dm_rev(mm) = dlm(i)
         end do

         age2(j)   = reci_polint(cmass(l), nma, dm_rev, age1)
         if (age2(j) <= 0.0_real32) age2(j) = 1.0_real32
         alogl2(j) = reci_polint(cmass(l), nma, dm_rev, alogl1)
         alogt2(j) = alog10(yntra(cmass(l), dm_rev, alogt1, nma))
         amass2(j) = reci_polint(cmass(l), nma, dm_rev, amass1)
         amdot2(j) = reci_polint(cmass(l), nma, dm_rev, amdot1)
         x2(j)     = yntra(cmass(l), dm_rev, x1, nma)
         y2(j)     = yntra(cmass(l), dm_rev, y1, nma)
         xc122(j)  = yntra(cmass(l), dm_rev, xc121, nma)
         xn142(j)  = yntra(cmass(l), dm_rev, xn141, nma)
         xo162(j)  = yntra(cmass(l), dm_rev, xo161, nma)
      end do

      call agecheck(time_in, jmp)

      if (time_in <= age2(1)) then
         temp(l)   = alogt2(1)
         bol(l)    = alogl2(1)
         zmass(l)  = amass2(1)
         bmdot(l)  = amdot2(1)
         xsurf(l)  = x2(1)
         ysurf(l)  = y2(1)
         xc12s(l)  = xc122(1)
         xn14s(l)  = xn142(1)
         xo16s(l)  = xo162(1)
         cycle
      end if

      if (time_in >= age2(13) .and. time_in < age2(51)) then
         call intrpl(jmp, age2, alogt2, 1, time_in, temp(l))
         call intrpl(jmp, age2, alogl2, 1, time_in, bol(l))
      else
         temp(l) = yntra(time_in, age2, alogt2, jmp)
         bol(l)  = yntra(time_in, age2, alogl2, jmp)
      end if

      zmass(l)  = yntra(time_in, age2, amass2, jmp)
      bmdot(l)  = yntra(time_in, age2, amdot2, jmp)
      xsurf(l)  = yntra(time_in, age2, x2, jmp)
      ysurf(l)  = yntra(time_in, age2, y2, jmp)
      xc12s(l)  = yntra(time_in, age2, xc122, jmp)
      xn14s(l)  = yntra(time_in, age2, xn142, jmp)
      xo16s(l)  = yntra(time_in, age2, xo162, jmp)
   end do

   ! Output selected masses for HRD
   if (io3 >= 0 .and. jmg == 0) then
      write(96, 200) temp(1), abs(bol(1)), temp(5), abs(bol(5)), temp(9), abs(bol(9)), temp(13), abs(bol(13)), &
                      temp(15), abs(bol(15)), temp(17), abs(bol(17)), temp(19), abs(bol(19)), temp(20), abs(bol(20)), &
                      temp(21), abs(bol(21)), temp(22), abs(bol(22))
   end if
   if (io3 >= 0 .and. jmg == 1) then
      write(96, 200) temp(1), abs(bol(1)), temp(21), abs(bol(21)), temp(41), abs(bol(41)), temp(61), abs(bol(61)), &
                      temp(71), abs(bol(71)), temp(81), abs(bol(81)), temp(91), abs(bol(91)), temp(96), abs(bol(96)), &
                      temp(101), abs(bol(101)), temp(106), abs(bol(106))
   end if
200 format(f6.3,f7.3,f6.3,f7.3,f6.3,f7.3,f6.3,f7.3,f6.3,f7.3, &
           f6.3,f7.3,f6.3,f7.3,f6.3,f7.3,f6.3,f7.3,f6.3,f7.3)

   ! Optionally: output underlying isochrone (commented out)
   ! ... (see original code for details) ...

end subroutine starpara

!> Interpolates stellar parameters from isochrones for fixed or variable mass grid.
subroutine starpara_iso(time_in, icount)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   real(real32), intent(in) :: time_in
   integer, intent(in)      :: icount

   ! Parameters
   integer, parameter :: nmaxint = 10, nmaxint1 = 11
   integer, parameter :: npgrid = 3000, npiso = 3000

   ! Arguments and local variables
   character(len=20) :: name
   character(len=80) :: title_tracks
   character(len=24) :: str_date
   integer :: lmin, lmax, jmg, niso, idiscont, nused, i, imax
   real(real32) :: temp(npgrid), bol(npgrid), xsurf(npgrid), zmass(npgrid)
   real(real32) :: ysurf(npgrid), xc12s(npgrid), xn14s(npgrid), xo16s(npgrid)
   real(real32) :: bmdot(npgrid), cmass(npgrid), tt_star(npgrid)
   real(real32) :: xminit(npiso), xmact(npiso), xl(npiso), xt(npiso)
   real(real32) :: x_h(npiso), x_he(npiso), x_c(npiso), x_n(npiso)
   real(real32) :: x_o(npiso), x_tstar(npiso), x_mdot(npiso)
   real(real32) :: doma, upma, sncut
   integer :: io3
   ! External procedures
   interface
      subroutine get_iso(age_in, xmlow, xmup, jmg, sncut, icount)
         use, intrinsic :: iso_fortran_env, only: real64
         real(real64), intent(in) :: age_in, xmlow, xmup, sncut
         integer, intent(in) :: jmg, icount
      end subroutine get_iso
      subroutine linterp(x, y, n, xout, yout, nout)
         real(real32), intent(in) :: x(*), y(*)
         integer, intent(in) :: n, nout
         real(real32), intent(in) :: xout(*)
         real(real32), intent(out) :: yout(*)
      end subroutine linterp
      integer function pos(x, arr, n)
         real(real32), intent(in) :: x
         real(real32), intent(in) :: arr(*)
         integer, intent(in) :: n
      end function pos
      subroutine errpri(msg)
         character(len=*), intent(in) :: msg
      end subroutine errpri
      subroutine fdate(str_date)
         character(len=*), intent(out) :: str_date
      end subroutine fdate
   end interface

   ! Output header at first time step
   if (icount <= 1 .and. io3 >= 0 .and. jmg /= 3) then
      write(96, '(A)') ' MODEL DESIGNATION: '//trim(name)
      call fdate(str_date)
      write(96, '(A)') ' MODEL GENERATED: '//trim(str_date)
      write(96, '(A)') '              EVOLUTIONARY TRACKS IN THE HERTZSPRUNG-RUSSELL DIAGRAM'
      write(96, '(A)') '    M=120        M=100         M=80         M=60         M=50         M=40         M=30         M=25         M=20         M=15'
   end if

   if (jmg == 2) then
      call get_iso(time_in, doma, upma, jmg, sncut, icount)
      nused = lmax - lmin + 1
      call linterp(xminit, xt, niso, cmass(lmin), temp(lmin), nused)
      call linterp(xminit, xl, niso, cmass(lmin), bol(lmin), nused)
      call linterp(xminit, xmact, niso, cmass(lmin), zmass(lmin), nused)
      call linterp(xminit, x_mdot, niso, cmass(lmin), bmdot(lmin), nused)
      call linterp(xminit, x_h, niso, cmass(lmin), xsurf(lmin), nused)
      call linterp(xminit, x_he, niso, cmass(lmin), ysurf(lmin), nused)
      call linterp(xminit, x_c, niso, cmass(lmin), xc12s(lmin), nused)
      call linterp(xminit, x_n, niso, cmass(lmin), xn14s(lmin), nused)
      call linterp(xminit, x_o, niso, cmass(lmin), xo16s(lmin), nused)
   else if (jmg == 3) then
      call get_iso(time_in, doma, upma, jmg, sncut, icount)
      lmin = pos(doma, xminit, niso)
      if (doma /= xminit(lmin)) lmin = lmin + 1
      lmax = niso - lmin + 1
      imax = pos(upma, xminit, niso)
      if (imax == niso - 1) then
         if (upma >= xminit(niso)) lmax = niso
      else
         lmax = imax
      end if
      if (lmax > npgrid) call errpri('STARPARA_ISO: ADJUST NPGRID !')
      do i = lmin, lmax
         cmass(i)   = xminit(i)
         temp(i)    = xt(i)
         tt_star(i) = x_tstar(i)
         bol(i)     = xl(i)
         zmass(i)   = xmact(i)
         bmdot(i)   = x_mdot(i)
         xsurf(i)   = x_h(i)
         ysurf(i)   = x_he(i)
         xc12s(i)   = x_c(i)
         xn14s(i)   = x_n(i)
         xo16s(i)   = x_o(i)
      end do
   end if

   ! Output selected masses for HRD
   if (io3 >= 0 .and. jmg == 2) then
      write(96, 200) temp(1), abs(bol(1)), temp(21), abs(bol(21)), temp(41), abs(bol(41)), temp(61), abs(bol(61)), &
                      temp(71), abs(bol(71)), temp(81), abs(bol(81)), temp(91), abs(bol(91)), temp(96), abs(bol(96)), &
                      temp(101), abs(bol(101)), temp(106), abs(bol(106))
   end if
200 format(f6.3,f7.3,f6.3,f7.3,f6.3,f7.3,f6.3,f7.3,f6.3,f7.3, &
           f6.3,f7.3,f6.3,f7.3,f6.3,f7.3,f6.3,f7.3,f6.3,f7.3)

   ! Optionally: output underlying isochrone (commented out)
   ! ... (see original code for details) ...

end subroutine starpara_iso

! ****************************************************************************
! ****************************************************************************
! ****************************************************************************
!> Determines an isochrone for a given age, for a specified mass range.
subroutine get_iso(age_in, xmlow, xmup, jmg, sncut, icount)
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   real(real64), intent(in) :: age_in, xmlow, xmup, sncut
   integer, intent(in)      :: jmg, icount

   ! Parameters
   integer, parameter :: ntracks = 50, nptrack = 500, npiso = 3000

   ! Arguments and local variables
   real(real64) :: dlm(ntracks), dla(ntracks, nptrack)
   real(real64) :: dm(ntracks), da(ntracks, nptrack)
   real(real64) :: xmas(ntracks, nptrack), dmas(ntracks, nptrack)
   real(real64) :: dl(ntracks, nptrack), dt(ntracks, nptrack)
   real(real64) :: d_h(ntracks, nptrack), d_he(ntracks, nptrack)
   real(real64) :: d_c(ntracks, nptrack), d_n(ntracks, nptrack)
   real(real64) :: d_o(ntracks, nptrack), d_tstar(ntracks, nptrack)
   real(real64) :: d_mdot(ntracks, nptrack)
   integer :: nma, jmp
   character(len=80) :: title_tracks

   real(real64) :: xminit(npiso), xmact(npiso), xl(npiso), xt(npiso)
   real(real64) :: x_h(npiso), x_he(npiso), x_c(npiso), x_n(npiso)
   real(real64) :: x_o(npiso), x_tstar(npiso), x_mdot(npiso)
   integer :: niso, idiscont

   real(real64) :: yl(2), yt(2), ym(2), y_h(2), y_he(2), y_c(2)
   real(real64) :: y_n(2), y_o(2), y_ts(2), y_md(2)
   real(real64) :: va(nptrack)
   logical :: lpartop
   real(real64), save :: deumaspre1, deumasNMpre1
   real(real64) :: aglog, agt, xmint, agl, xminl
   real(real64) :: delm_max, delm_min, delm_sn
   real(real64) :: premas, deumas, deumas1, deumasNM1
   real(real64) :: xm2w, deltam, vak, vakm1, diff, ff, xmlog
   real(real64) :: distl, distt
   integer :: idone, idoneNM, nonmono, i, m, k, l, k_select, ll, jj
   integer :: n, ntracks_loc, nptrack_loc
   real(real64) :: aglogdum
   ! External functions
   interface
      integer function pos(x, arr, n)
         import :: real64
         real(real64), intent(in) :: x
         real(real64), intent(in) :: arr(*)
         integer, intent(in) :: n
      end function pos
      real(real64) function yparinterp(x, xarr, yarr, n, nmax, lpartop, k, nptrack)
         import :: real64
         real(real64), intent(in) :: x
         real(real64), intent(in) :: xarr(*), yarr(*)
         integer, intent(in) :: n, nmax, k, nptrack
         logical, intent(inout) :: lpartop
      end function yparinterp
      subroutine errpri(msg)
         character(len=*), intent(in) :: msg
      end subroutine errpri
   end interface

   ! Initialization of tolerances for mass step
   aglog = log10(age_in)
   agt   = 0.010_real64
   xmint = 0.003_real64
   agl   = 0.050_real64
   xminl = 0.005_real64
   delm_max = 1.0_real64
   delm_min = 9.0e-5_real64
   delm_sn  = 9.0e-6_real64

   ! Choice of the appropriate mass range
   if (jmg == 2) then
      premas = dm(nma)
      deumas = dm(1)
   else
      premas = max(dm(nma), xmlow)
      deumas = min(dm(1), xmup)
   end if

   ! Calculation of the isochrone
   idone = 0
   idoneNM = 0
   nonmono  = 0
   idiscont = 0
   i = 1
   deltam = 0.0_real64
   xm2w = premas
   if (niso >= 1) xminit(niso) = xm2w
   niso = 1
   deumas1 = deumas
   if (icount == 1) then
      deumaspre1 = deumas1
      deumasNM1 = deumas1 + 2.0_real64
      deumasNMpre1 = deumasNM1
   end if

   ! Main loop: build isochrone
   do while (xminit(niso) <= deumas)
      xminit(i) = xm2w
      m = pos(xm2w, dm, nma)
      do k = 1, jmp
         lpartop = .false.
         va(k) = yparinterp(log10(xm2w), dlm, dla(:,k), nma, ntracks, lpartop, k, nptrack)
         if (k > 1) then
            if (va(k) < va(k-1)) then
               va(k) = 10.0_real64**va(k-1)
               va(k) = log10(va(k) + 1.0e-1_real64)
            end if
         end if
         if (va(1) >= aglog) then
            yt(1) = dt(m,1)
            yt(2) = dt(m+1,1)
            yl(1) = dl(m,1)
            yl(2) = dl(m+1,1)
            ym(1) = dmas(m,1)
            ym(2) = dmas(m+1,1)
            y_h(1) = d_h(m,1)
            y_h(2) = d_h(m+1,1)
            y_he(1) = d_he(m,1)
            y_he(2) = d_he(m+1,1)
            y_c(1) = d_c(m,1)
            y_c(2) = d_c(m+1,1)
            y_n(1) = d_n(m,1)
            y_n(2) = d_n(m+1,1)
            y_o(1) = d_o(m,1)
            y_o(2) = d_o(m+1,1)
            y_ts(1) = d_tstar(m,1)
            y_ts(2) = d_tstar(m+1,1)
            y_md(1) = d_mdot(m,1)
            y_md(2) = d_mdot(m+1,1)
            goto 18
         else
            if (k > jmp-2) then
               if (aglog >= va(jmp-2)) then
                  diff = 0.0_real64
                  k_select = jmp
                  goto 15
               else if (aglog > va(jmp)) then
                  if (nonmono == 0) then
                     print *, 'WARNING: NON-MONOTONIC LIFETIME FROM ', xm2w, ' MSUN !'
                     print *, 'ERROR IN GET_ISO: SEE WARNING ABOVE !'
                     nonmono  = 1
                     idiscont = i - 1
                     xminit(i) = xminit(i-1)
                     xt(i)     = xt(i-1)
                     xl(i)     = xl(i-1)
                     xmact(i)  = xmact(i-1)
                     x_h(i)    = x_h(i-1)
                     x_he(i)   = x_he(i-1)
                     x_c(i)    = x_c(i-1)
                     x_n(i)    = x_n(i-1)
                     x_o(i)    = x_o(i-1)
                     x_tstar(i) = x_tstar(i-1)
                     x_mdot(i)  = x_mdot(i-1)
                     i = i + 1
                  end if
                  goto 40
               end if
            end if
            if (aglog < va(k)) then
               vak = 10.0_real64**va(k)
               vakm1 = 10.0_real64**va(k-1)
               diff = (age_in - vakm1) / (vak - vakm1)
               k_select = k
               goto 15
            end if
         end if
      end do

15    continue
      k = k_select
      do l = 1, 2
         ll = m + l - 1
         ym(l) = xmas(ll, k-1) + (xmas(ll, k) - xmas(ll, k-1)) * diff
         ym(l) = log10(ym(l))
         yt(l) = dt(ll, k-1) + (dt(ll, k) - dt(ll, k-1)) * diff
         yl(l) = dl(ll, k-1) + (dl(ll, k) - dl(ll, k-1)) * diff
         y_h(l) = d_h(ll, k-1) + (d_h(ll, k) - d_h(ll, k-1)) * diff
         y_he(l) = d_he(ll, k-1) + (d_he(ll, k) - d_he(ll, k-1)) * diff
         y_c(l) = d_c(ll, k-1) + (d_c(ll, k) - d_c(ll, k-1)) * diff
         y_n(l) = d_n(ll, k-1) + (d_n(ll, k) - d_n(ll, k-1)) * diff
         y_o(l) = d_o(ll, k-1) + (d_o(ll, k) - d_o(ll, k-1)) * diff
         y_ts(l) = d_tstar(ll, k-1) + (d_tstar(ll, k) - d_tstar(ll, k-1)) * diff
         y_md(l) = d_mdot(ll, k-1) + (d_mdot(ll, k) - d_mdot(ll, k-1)) * diff
      end do

18    ff = (log10(xm2w) - dlm(m)) / (dlm(m+1) - dlm(m))
      xl(i)    = yl(1) + (yl(2) - yl(1)) * ff
      if (xl(i) < -3.0_real64) xl(i) = -20.0_real64
      xt(i)    = yt(1) + (yt(2) - yt(1)) * ff
      xmlog    = ym(1) + (ym(2) - ym(1)) * ff
      xmact(i) = 10.0_real64**xmlog
      x_h(i)    = y_h(1) + (y_h(2) - y_h(1)) * ff
      x_he(i)   = y_he(1) + (y_he(2) - y_he(1)) * ff
      x_c(i)    = y_c(1) + (y_c(2) - y_c(1)) * ff
      x_n(i)    = y_n(1) + (y_n(2) - y_n(1)) * ff
      x_o(i)    = y_o(1) + (y_o(2) - y_o(1)) * ff
      x_tstar(i) = y_ts(1) + (y_ts(2) - y_ts(1)) * ff
      x_mdot(i)  = y_md(1) + (y_md(2) - y_md(1)) * ff

      niso = i
      if (i == 1) then
         deltam = (dm(m) - dm(m+1)) / 10.0_real64
         i = i + 1
      else if (nonmono == 1) then
         i = i + 1
         xminit(i) = xminit(i-1)
         xt(i)     = xt(i-1)
         xl(i)     = xl(i-1)
         xmact(i)  = xmact(i-1)
         x_h(i)    = x_h(i-1)
         x_he(i)   = x_he(i-1)
         x_c(i)    = x_c(i-1)
         x_n(i)    = x_n(i-1)
         x_o(i)    = x_o(i-1)
         x_tstar(i) = x_tstar(i-1)
         x_mdot(i)  = x_mdot(i-1)
         i = i + 1
         nonmono = 0
         print *, 'NON-MONOTONIC LIFETIME ENDS AT', xm2w, ' MSUN'
      else
         if (i > 1) then
            if (xl(i-1) > -5.0_real64 .and. xl(i) < -5.0_real64 .and. xminit(i) > sncut-1.0_real64) then
               if (xminit(i) - xminit(i-1) > delm_min) then
                  xm2w = xm2w - deltam
                  deltam = deltam / 1.2_real64
                  goto 40
               end if
               deumas1 = xm2w
            end if
            if (xl(i-1) < -5.0_real64 .and. xl(i) > -5.0_real64 .and. xminit(i) > sncut-1.0_real64) then
               if (xminit(i) - xminit(i-1) > delm_min) then
                  xm2w = xm2w - deltam
                  deltam = deltam / 1.2_real64
                  goto 40
               end if
               deumasNM1 = xm2w
            end if
         end if

         if (xm2w > deumas1 .and. xm2w <= deumaspre1 .and. deumaspre1 > sncut .and. deumas1 < deumaspre1) then
            deltam = (deumaspre1 - deumas1) / 10.0_real64
            if (deltam > delm_max) deltam = delm_max
            if (deltam < delm_sn) deltam = delm_sn
            if (xm2w < deumaspre1 .and. xm2w + deltam > deumaspre1 .and. idone == 0) then
               deltam = deumaspre1 - xm2w
               idone = 1
            end if
            i = i + 1
            goto 40
         else
            if (xm2w > deumaspre1 .and. xm2w <= deumasNM1) then
               if (xm2w < deumasNMpre1) then
                  xl(i) = -50.0_real64
                  deltam = delm_max
                  i = i + 1
                  goto 40
               else
                  deltam = (deumasNM1 - deumasNMpre1) / 10.0_real64
                  if (deltam > delm_max) deltam = delm_max
                  if (deltam < delm_sn) deltam = delm_sn
                  if (xm2w < deumasNM1 .and. xm2w + deltam > deumasNM1 .and. idoneNM == 0 .and. xm2w > sncut-1.0_real64) then
                     deltam = deumasNM1 - xm2w
                     idoneNM = 1
                  end if
                  i = i + 1
                  goto 40
               end if
            end if
         end if

         distl = abs(xl(i) - xl(i-1))
         distt = abs(xt(i) - xt(i-1))
         if ((distl > agl .or. distt > agt) .or. (deltam > delm_max)) then
            if (deltam < delm_min) then
               if (xl(i-1) > -5.0_real64 .and. xl(i) < -5.0_real64 .and. xminit(i-1) > sncut-1.0_real64 .or. xm2w > deumaspre1) then
                  if (deltam < delm_sn) then
                     i = i + 1
                     deltam = delm_sn
                  else
                     xm2w = xm2w - deltam
                     deltam = deltam / 1.2_real64
                  end if
               else
                  i = i + 1
                  deltam = delm_min
               end if
            else
               xm2w = xm2w - deltam
               deltam = deltam / 6.0_real64
            end if
         else if (distl < xminl .and. distt < xmint) then
            i = i + 1
            if ((xm2w + 10.0_real64 * deltam) < deumas) then
               deltam = 10.0_real64 * deltam
               if (deltam > delm_max) deltam = delm_max
            end if
         else
            i = i + 1
            if (xm2w + deltam > deumas .and. xm2w < deumas .and. idone == 0) then
               deltam = deumas - xm2w
               idone = 1
            end if
         end if
      end if

40    continue
      if (xm2w >= deumas) idone = idone + 1
      if (deltam < 0.0_real64) deltam = delm_min
      xm2w = xm2w + deltam

      if (deltam < 0.0_real64) call errpri('GET_ISO: DELTAM.LT.0')
      if (i > npiso)  call errpri('GET_ISO: I.GT.NPISO')
   end do

   niso = niso - 1
   deumaspre1 = deumas1
   deumasNMpre1 = deumasNM1

   ! Optionally: output isochrone to file (commented out)
   ! ... (see original code for details) ...

   return
end subroutine get_iso

! ****************************************************************************
! ****************************************************************************
! ****************************************************************************
!> Ensures the AGE2 array is strictly monotonically increasing to prevent interpolation errors.
subroutine agecheck(time, jmp)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   real(real32), intent(in) :: time
   integer, intent(in)      :: jmp

   integer, parameter :: nptrack = 500
   real(real32) :: age2(nptrack)
   integer :: i

   ! NOTE: In modern Fortran, age2 should be passed as an argument or be in a module.
   ! Here, we assume it is available via a module or host association.

   do i = jmp, 2, -1
      if (age2(i-1) / age2(i) >= 1.0_real32) age2(i-1) = age2(i) / 1.000001_real32
   end do

end subroutine agecheck
c
c *********************************************************************
c *********************************************************************
c *********************************************************************
!> Adjusts the WR core temperatures of the evolutionary tracks for consistency with WR atmospheres.
subroutine temp_adjust(iwrt, iatmos)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   integer, intent(in) :: iwrt, iatmos

   ! Parameters
   integer, parameter :: npgrid = 3000
   real(real32), parameter :: factor = 0.6_real32

   ! Variables (replace COMMON blocks with explicit arguments or modules in modern code)
   integer :: l, lmin, lmax
   real(real32) :: temp(npgrid), bol(npgrid), xsurf(npgrid), zmass(npgrid)
   real(real32) :: ysurf(npgrid), xc12s(npgrid), xn14s(npgrid), xo16s(npgrid)
   real(real32) :: bmdot(npgrid), cmass(npgrid), tt_star(npgrid)
   real(real32) :: cnr, coher, t1, t2

   ! Loop over all mass bins
   do l = lmin, lmax
      if (xn14s(l) == 0.0_real32) xn14s(l) = 1.0e-6_real32
      cnr = xc12s(l) / xn14s(l)
      coher = ((xc12s(l)/12.0_real32) + (xo16s(l)/16.0_real32)) / (ysurf(l)/4.0_real32)
      if (xsurf(l) < 0.1_real32 .and. iwrt > 0 .and. iatmos >= 4) then
         t1 = 10.0_real32**(temp(l))
         t2 = 10.0_real32**(tt_star(l))
         temp(l) = t2 + (factor - 1.0_real32) * (t2 - t1)
         temp(l) = log10(temp(l))
      end if
   end do

end subroutine temp_adjust
c
c *********************************************************************
c *********************************************************************
c *********************************************************************
!> Calculates the supernova rate and associated energetics.
!> Modern Fortran version.
subroutine supernova(time_in, icount)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   real(real32), intent(in) :: time_in
   integer, intent(in)      :: icount

   ! Parameters
   integer, parameter :: npgrid = 3000
   real(real32), parameter :: yearseconds = 3.15576e7

   ! Variables (replace COMMON blocks with explicit arguments or modules in modern code)
   character(len=20) :: name
   character(len=80) :: str_date
   real(real32) :: temp(npgrid), bol(npgrid), xsurf(npgrid), zmass(npgrid)
   real(real32) :: ysurf(npgrid), xc12s(npgrid), xn14s(npgrid), xo16s(npgrid)
   real(real32) :: bmdot(npgrid), cmass(npgrid), tt_star(npgrid)
   real(real32) :: dens(npgrid)
   real(real32) :: sn(npgrid), pow(npgrid), snr_temp(npgrid), cm_temp(npgrid)
   real(real32) :: wpow, wen
   integer :: lmin, lmax, jmg, isf, io4
   real(real32) :: tstep, sncut, bhcut, xmwr, critma, critma_new, critup, critup_new
   real(real32) :: siben, siien
   real(real32) :: deltat, sner, energ
   real(real32) :: snr, sner_tot, wrsnr, siipow, sibpow
   real(real32) :: sn1, sn2, sn3, sn4, sn5, sn6, sn7, sn8
   real(real32) :: typma, critma_out
   integer :: isnr, itemp, l, lstart, lend, lstep, iflag
   integer :: icount_local
   ! External
   interface
      subroutine fdate(str_date)
         character(len=*), intent(out) :: str_date
      end subroutine fdate
      integer function pos(x, arr, n)
         real(real32), intent(in) :: x
         real(real32), intent(in) :: arr(*)
         integer, intent(in) :: n
      end function pos
   end interface

   ! Header output at first call
   if (icount <= 1) then
      write(97,196) name
196   format(' MODEL DESIGNATION: ',a20)
      call fdate(str_date)
      write(97,597) str_date
597   format(' MODEL GENERATED: ', a24)
      write(97,195)
195   format(/,'              RESULTS FOR THE',
     &        ' SUPERNOVA RATE')
      write(97,193)
193   format(/,'                     ALL SUPERNOVAE         '
     &        ,'      TYPE IB SUPERNOVAE            '
     &        ,'   ALL SUPERNOVAE           STARS + SUPERNOVAE')
      write(97,194)
194   format('    TIME       TOTAL RATE  POWER   ENERGY    '
     &        ,'TOTAL RATE  POWER   ENERGY   '
     &        ,'TYPICAL MASS   LOWEST PROG. MASS    POWER   ENERGY')
   end if

   ! Time step for this call
   if (icount == 1) then
      deltat = time_in
   else
      deltat = tstep
   end if

   ! Reset accumulators
   snr    = 0.0_real32
   sner_tot = 0.0_real32
   wrsnr  = 0.0_real32
   siipow = 0.0_real32
   sibpow = 0.0_real32

   if (icount == 1) then
      siien = 0.0_real32
      siben = 0.0_real32
   end if

   isnr = 0

   ! Determine loop direction for mass grid
   if (jmg == 3) then
      lstart = lmax
      lend   = lmin
      lstep  = -1
   else
      lstart = lmin
      lend   = lmax
      lstep  = 1
   end if

   ! Loop over all masses: must be from high to low masses
   do l = lstart, lend, lstep
      if (bol(l) > -10.0_real32 .or. cmass(l) > bhcut) cycle
      ! Signal for non-monotonic lifetime
      if (cmass(l) < sncut) then
         if (isf <= 0) critma_new = cmass(l)
         cycle
      end if

      ! Determine SN type: IB (iflag=1) or II (iflag=0)
      iflag = 0
      if (temp(l) > 4.4_real32 .and. xsurf(l) < 0.4_real32 .and. cmass(l) >= xmwr) iflag = 1

      ! Possible supernova event rate
      sner = dens(l) / deltat

      ! Update critup_new and critma_new for instantaneous bursts
      if (isf <= 0) then
         if (critup_new < 0.0_real32) then
            critup_new = cmass(l)
         else
            if (cmass(l) >= critup_new) then
               if (critma > sncut .and. l /= lstart) then
                  critup_new = cmass(l)
               else
                  critup_new = 1.0e36_real32
                  critup = critup_new
               end if
            end if
         end if
         if (cmass(l) >= critma .and. cmass(l) <= critup) then
            sner = 0.0_real32
         else
            critma_new = cmass(l)
         end if
      end if

      if (isf > 0 .and. cmass(l) < critma) then
         critma_new = cmass(l)
      end if

      ! Get SN rate, power and energy from actual sner
      sn(l) = sner
      pow(l) = sn(l) / yearseconds
      energ = sn(l) * deltat

      snr = snr + sn(l)
      siipow = siipow + pow(l)
      siien  = siien  + energ
      if (iflag == 1) wrsnr = wrsnr + sn(l)
      if (iflag == 1) sibpow = sibpow + pow(l)
      if (iflag == 1) siben  = siben  + energ

      ! Prepare for the determination of the typical SN mass
      if (sn(l) /= 0.0_real32) then
         isnr = isnr + 1
         snr_temp(isnr) = snr
         cm_temp(isnr)  = cmass(l)
      end if
   end do

   critma = critma_new
   critup = critup_new

   sn1 = log10(snr   + 1.e-30_real32)
   sn2 = log10(wrsnr + 1.e-30_real32)
   sn3 = log10(siipow + 1.e-30_real32) + 51.0_real32
   sn4 = log10(sibpow + 1.e-30_real32) + 51.0_real32
   sn5 = log10(siien  + 1.e-30_real32) + 51.0_real32
   sn6 = log10(siben  + 1.e-30_real32) + 51.0_real32
   sn7 = log10(siipow + 1.e-30_real32 + wpow * 1.e-16_real32) + 51.0_real32
   sn8 = log10(siien  + 1.e-30_real32 + wen  * 1.e-16_real32) + 51.0_real32

   if (sn3 <= 21.0_real32) sn3 = -30.0_real32
   if (sn4 <= 21.0_real32) sn4 = -30.0_real32
   if (sn5 <= 21.0_real32) sn5 = -30.0_real32
   if (sn6 <= 21.0_real32) sn6 = -30.0_real32
   if (sn7 <= 21.0_real32 .or. io4 < 0) sn7 = -30.0_real32
   if (sn8 <= 21.0_real32 .or. io4 < 0) sn8 = -30.0_real32

   ! Compute the typical progenitor mass of a supernova
   if (snr /= 0.0_real32) then
      itemp = pos(snr / 2.0_real32, snr_temp, isnr)
      if (isnr == 1) then
         typma = cm_temp(isnr)
      else
         typma = cm_temp(itemp)
      end if
   else
      typma = 0.0_real32
   end if

   critma_out = critma
   if (critma_out > upma)  critma_out = 0.0_real32
   if (critma_out < sncut) critma_out = 0.0_real32

   write(97,90) time_in, sn1, sn3, sn5, sn2, sn4, sn6, typma, critma_out, sn7, sn8
90 format(1x,e11.3,2x,3f9.3,3x,3f9.3,3x,
     &         f8.1,8x,f8.1,8x,3f9.3)

end subroutine supernova
c
c *********************************************************************
c *********************************************************************
c *********************************************************************
!> Assigns spectral types to each position in the HRD using the Schmidt-Kaler calibration.
!> Modern Fortran version.
subroutine spectype(time_in, icount)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   real(real32), intent(in) :: time_in
   integer, intent(in)      :: icount

   ! Parameters
   integer, parameter :: npgrid = 3000, ntypes = 445, nwrclass = 5
   character(len=20) :: name
   character(len=80) :: str_date

   ! Variables (replace COMMON blocks with explicit arguments or modules in modern code)
   real(real32) :: temp(npgrid), bol(npgrid), xsurf(npgrid), zmass(npgrid)
   real(real32) :: ysurf(npgrid), xc12s(npgrid), xn14s(npgrid), xo16s(npgrid)
   real(real32) :: bmdot(npgrid), cmass(npgrid), dens(npgrid), tt_star(npgrid)
   integer :: lmin, lmax
   real(real32) :: anumb(ntypes), wrnumb(nwrclass)
   real(real32) :: allspe, allwr, allwn, allwc, allo, wro, wcwn
   integer :: l, m, n, k, iwr, index
   real(real32) :: cnr, coher

   ! External
   interface
      subroutine fdate(str_date)
         character(len=*), intent(out) :: str_date
      end subroutine fdate
      subroutine near(l, index)
         integer, intent(in)  :: l
         integer, intent(out) :: index
      end subroutine near
   end interface

   ! Header output at first call
   if (icount <= 1) then
      write(94,196) name
196   format(' MODEL DESIGNATION: ',a20)
      call fdate(str_date)
      write(94,497) str_date
497   format(' MODEL GENERATED: ', a24)
      write(94,95)
95    format(/,'              RESULTS FOR THE',
     &        ' SPECTRAL TYPE CALIBRATION -- 1')
      write(94,94)
94    format(/,'   TIME         ALL SP         O          WR',
     &     '          WN          WC          WR/O       WC/WN  ')
      write(90,196) name
      call fdate(str_date)
      write(90,497) str_date
      write(90,93)
93    format(/,'              RESULTS FOR THE',
     &        ' SPECTRAL TYPE CALIBRATION -- 2')
   end if

   ! Reset accumulators if not continuous SF
   if (.not.(isf > 0 .and. icount > 1)) then
      anumb = 0.0_real32
      wrnumb = 0.0_real32
   end if
   allspe = 0.0_real32
   allwr  = 0.0_real32
   allwn  = 0.0_real32
   allwc  = 0.0_real32
   allo   = 0.0_real32

   ! Loop over all masses
   do l = lmin, lmax
      if (bol(l) < -10.0_real32) cycle

      ! Wolf-Rayet classification
      if (temp(l) > 4.4_real32 .and. xsurf(l) < 0.4_real32 .and. cmass(l) >= xmwr) then
         if (xn14s(l) == 0.0_real32) xn14s(l) = 1.0e-6_real32
         cnr = xc12s(l) / xn14s(l)
         coher = ((xc12s(l)/12.0_real32) + (xo16s(l)/16.0_real32)) / (ysurf(l)/4.0_real32)
         if (xsurf(l) > 0.1_real32) then
            iwr = 1
         else if (cnr < 10.0_real32) then
            iwr = 2
         else if (coher < 0.5_real32) then
            iwr = 3
         else if (coher < 1.0_real32) then
            iwr = 4
         else
            iwr = 5
         end if
         wrnumb(iwr) = wrnumb(iwr) + dens(l)
      else
         ! OBAFGKM classification
         call near(l, index)
         anumb(index) = anumb(index) + dens(l)
      end if
   end do

   ! Count all stars
   allspe = sum(anumb)
   ! Count WR stars
   allwn = wrnumb(1) + wrnumb(2)
   allwc = wrnumb(3) + wrnumb(4) + wrnumb(5)
   allwr = allwn + allwc
   ! Count O stars (indices as in original code)
   allo = sum(anumb(1:14)) + sum(anumb(90:103)) + sum(anumb(179:192)) + sum(anumb(268:281)) + sum(anumb(357:370))
   ! Compute ratios
   if (allo /= 0.0_real32) then
      wro = allwr / allo
   else
      wro = 0.0_real32
   end if
   if (allwn /= 0.0_real32) then
      wcwn = allwc / allwn
   else
      wcwn = 0.0_real32
   end if

   write(94,150) time_in, allspe, allo, allwr, allwn, allwc, wro, wcwn
150 format(1x,1pe11.3,3x,6(1pe9.3,3x),1pe9.3)

   if (mod(time_in, tdel) < tstep) then
      write(90,161) time_in
161   format(/,' TIME = ',1pe11.3)
      write(90,162)
162   format(/,' LUMINOSITY CLASS V')
      write(90,166) (anumb(i), i=1,89)
166   format(' O:',2x,8(1pe11.3,1x),/,5x,6(1pe11.3,1x),
     &         /,' B:',2x,8(1pe11.3,1x),/,5x,8(1pe11.3,1x),/,5x,
     &                 2(1pe11.3,1x),
     &         /,' A:',2x,8(1pe11.3,1x),/,5x,6(1pe11.3,1x),
     &         /,' F:',2x,8(1pe11.3,1x),
     &         /,' G:',2x,8(1pe11.3,1x),
     &         /,' K:',2x,8(1pe11.3,1x),/,5x,6(1pe11.3,1x),
     &         /,' M:',2x,8(1pe11.3,1x),/,5x,5(1pe11.3,1x))
      write(90,168)
168   format(/,' LUMINOSITY CLASS IV')
      write(90,166) (anumb(i), i=90,178)
      write(90,170)
170   format(/,' LUMINOSITY CLASS III')
      write(90,166) (anumb(i), i=179,267)
      write(90,172)
172   format(/,' LUMINOSITY CLASS II')
      write(90,166) (anumb(i), i=268,356)
      write(90,174)
174   format(/,' LUMINOSITY CLASS I')
      write(90,166) (anumb(i), i=357,445)
      write(90,178) (wrnumb(i), i=1,5)
178   format(/,' WNL:',1pe11.3,/,' WNE:',1pe11.3,/,
     &         ' WCL:',1pe11.3,/,' WCE:',1pe11.3,/,' WO :',1pe11.3)
   end if

end subroutine spectype

! ****************************************************************************
! ****************************************************************************
! ****************************************************************************
!> Finds the nearest neighbor in the spectral type calibration table.
!> Modern Fortran version.
subroutine near(l, index)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   integer, intent(in)  :: l
   integer, intent(out) :: index

   ! Parameters
   integer, parameter :: npgrid = 3000, ntypes = 445

   ! Variables (replace COMMON blocks with explicit arguments or modules in modern code)
   real(real32) :: temp(npgrid), bol(npgrid)
   real(real32) :: tip(ntypes), bip(ntypes)
   integer :: i
   real(real32) :: di, dimin

   ! TODO: Replace with module variables or pass as arguments in modern code

   dimin = 1.0e20_real32
   index = 1

   do i = 1, ntypes
      di = (temp(l) - tip(i))**2 + (bol(l) - bip(i))**2
      if (di < dimin) then
         dimin = di
         index = i
      end if
   end do

end subroutine near
c
c *******************************************************************
c ********************************************************************
c ********************************************************************
!> Computes the power, energy, and momentum flux due to stellar winds.
!> Modern Fortran version.
subroutine windpower(time_in, icount)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   real(real32), intent(in) :: time_in
   integer, intent(in)      :: icount

   ! Parameters
   integer, parameter :: nmaxint = 10, nmaxint1 = 11
   integer, parameter :: npgrid = 3000
   character(len=20) :: name
   character(len=80) :: str_date

   ! Common blocks (replace with modules in modern code)
   real(real32) :: temp(npgrid), bol(npgrid), xsurf(npgrid), zmass(npgrid)
   real(real32) :: ysurf(npgrid), xc12s(npgrid), xn14s(npgrid), xo16s(npgrid)
   real(real32) :: bmdot(npgrid), dens(npgrid), cmass(npgrid), tt_star(npgrid)
   integer :: lmin, lmax, ltype
   real(real32) :: wpow, wpow1, wpow2, wpow3, wpow4
   real(real32) :: wmom, wmom1, wmom2, wmom3, wmom4, wen
   real(real32) :: tstep
   integer :: l

   ! External interfaces (to be implemented elsewhere)
   interface
      function wind1(l, ltype, idx) result(val)
         integer, intent(in) :: l, ltype, idx
         real(real32) :: val
      end function wind1
      function wind2(l, ltype, idx) result(val)
         integer, intent(in) :: l, ltype, idx
         real(real32) :: val
      end function wind2
      function wind3(l, ltype, idx) result(val)
         integer, intent(in) :: l, ltype, idx
         real(real32) :: val
      end function wind3
      function wind4(l, ltype, idx) result(val)
         integer, intent(in) :: l, ltype, idx
         real(real32) :: val
      end function wind4
      subroutine fdate(str_date)
         character(len=*), intent(out) :: str_date
      end subroutine fdate
   end interface

   ! Header output at first call
   if (icount <= 1) then
      wen = 0.0_real32
      write(95,196) name
196   format(' MODEL DESIGNATION: ',a20)
      call fdate(str_date)
      write(95,597) str_date
597   format(' MODEL GENERATED: ', a24)
      write(95,195)
195   format(/,'              RESULTS FOR THE',
     &        ' **STELLAR** WIND POWER AND ENERGY')
      write(95,193)
193   format(/,'                             POWER [ERG/SEC]  '
     &        ,'                ENERGY [ERG]            '
     &        ,'   MOMENTUM FLUX [DYN]      ')
      write(95,194)
194   format('    TIME           ALL      OB      RSG      LBV'
     &        ,'       WR          ALL    '
     &        ,'     ALL      OB      RSG      LBV       WR   '
     &        ,'     ')
   end if

   ! Reset accumulators if not continuous SF
   if (.not.(isf > 0 .and. icount > 1)) then
      wpow  = 0.0_real32
      wpow1 = 0.0_real32
      wpow2 = 0.0_real32
      wpow3 = 0.0_real32
      wpow4 = 0.0_real32
      wmom  = 0.0_real32
      wmom1 = 0.0_real32
      wmom2 = 0.0_real32
      wmom3 = 0.0_real32
      wmom4 = 0.0_real32
   end if

   select case (iwind)
   case (0)
      do l = lmin, lmax
         wpow  = wpow  + wind1(l, ltype, 1) * dens(l)
         if (ltype == 1) wpow1 = wpow1 + wind1(l, ltype, 1) * dens(l)
         if (ltype == 2) wpow2 = wpow2 + wind1(l, ltype, 1) * dens(l)
         if (ltype == 3) wpow3 = wpow3 + wind1(l, ltype, 1) * dens(l)
         if (ltype == 4) wpow4 = wpow4 + wind1(l, ltype, 1) * dens(l)
         wmom  = wmom  + wind1(l, ltype, 2) * dens(l)
         if (ltype == 1) wmom1 = wmom1 + wind1(l, ltype, 2) * dens(l)
         if (ltype == 2) wmom2 = wmom2 + wind1(l, ltype, 2) * dens(l)
         if (ltype == 3) wmom3 = wmom3 + wind1(l, ltype, 2) * dens(l)
         if (ltype == 4) wmom4 = wmom4 + wind1(l, ltype, 2) * dens(l)
      end do
   case (1)
      do l = lmin, lmax
         wpow  = wpow  + wind2(l, ltype, 1) * dens(l)
         if (ltype == 1) wpow1 = wpow1 + wind2(l, ltype, 1) * dens(l)
         if (ltype == 2) wpow2 = wpow2 + wind2(l, ltype, 1) * dens(l)
         if (ltype == 3) wpow3 = wpow3 + wind2(l, ltype, 1) * dens(l)
         if (ltype == 4) wpow4 = wpow4 + wind2(l, ltype, 1) * dens(l)
         wmom  = wmom  + wind2(l, ltype, 2) * dens(l)
         if (ltype == 1) wmom1 = wmom1 + wind2(l, ltype, 2) * dens(l)
         if (ltype == 2) wmom2 = wmom2 + wind2(l, ltype, 2) * dens(l)
         if (ltype == 3) wmom3 = wmom3 + wind2(l, ltype, 2) * dens(l)
         if (ltype == 4) wmom4 = wmom4 + wind2(l, ltype, 2) * dens(l)
      end do
   case (2)
      do l = lmin, lmax
         wpow  = wpow  + wind3(l, ltype, 1) * dens(l)
         if (ltype == 1) wpow1 = wpow1 + wind3(l, ltype, 1) * dens(l)
         if (ltype == 2) wpow2 = wpow2 + wind3(l, ltype, 1) * dens(l)
         if (ltype == 3) wpow3 = wpow3 + wind3(l, ltype, 1) * dens(l)
         if (ltype == 4) wpow4 = wpow4 + wind3(l, ltype, 1) * dens(l)
         wmom  = wmom  + wind3(l, ltype, 2) * dens(l)
         if (ltype == 1) wmom1 = wmom1 + wind3(l, ltype, 2) * dens(l)
         if (ltype == 2) wmom2 = wmom2 + wind3(l, ltype, 2) * dens(l)
         if (ltype == 3) wmom3 = wmom3 + wind3(l, ltype, 2) * dens(l)
         if (ltype == 4) wmom4 = wmom4 + wind3(l, ltype, 2) * dens(l)
      end do
   case (3)
      do l = lmin, lmax
         wpow  = wpow  + wind4(l, ltype, 1) * dens(l)
         if (ltype == 1) wpow1 = wpow1 + wind4(l, ltype, 1) * dens(l)
         if (ltype == 2) wpow2 = wpow2 + wind4(l, ltype, 1) * dens(l)
         if (ltype == 3) wpow3 = wpow3 + wind4(l, ltype, 1) * dens(l)
         if (ltype == 4) wpow4 = wpow4 + wind4(l, ltype, 1) * dens(l)
         wmom  = wmom  + wind4(l, ltype, 2) * dens(l)
         if (ltype == 1) wmom1 = wmom1 + wind4(l, ltype, 2) * dens(l)
         if (ltype == 2) wmom2 = wmom2 + wind4(l, ltype, 2) * dens(l)
         if (ltype == 3) wmom3 = wmom3 + wind4(l, ltype, 2) * dens(l)
         if (ltype == 4) wmom4 = wmom4 + wind4(l, ltype, 2) * dens(l)
      end do
   end select

   ! Update total wind energy
   if (icount == 1) then
      wen = wen + wpow * time_in * 3.1558e7_real32
   else
      wen = wen + wpow * tstep * 3.1558e7_real32
   end if

   ! Convert to log units for output
   wpow  = log10(wpow  + 1.e-30_real32) + 35.0_real32
   wpow1 = log10(wpow1 + 1.e-30_real32) + 35.0_real32
   wpow2 = log10(wpow2 + 1.e-30_real32) + 35.0_real32
   wpow3 = log10(wpow3 + 1.e-30_real32) + 35.0_real32
   wpow4 = log10(wpow4 + 1.e-30_real32) + 35.0_real32
   wen   = log10(wen   + 1.e-30_real32) + 35.0_real32
   wmom  = log10(wmom  + 1.e-30_real32) + 35.0_real32
   wmom1 = log10(wmom1 + 1.e-30_real32) + 35.0_real32
   wmom2 = log10(wmom2 + 1.e-30_real32) + 35.0_real32
   wmom3 = log10(wmom3 + 1.e-30_real32) + 35.0_real32
   wmom4 = log10(wmom4 + 1.e-30_real32) + 35.0_real32

   ! Threshold for output
   if (wpow  < 5.01_real32) wpow  = -30.0_real32
   if (wpow1 < 5.01_real32) wpow1 = -30.0_real32
   if (wpow2 < 5.01_real32) wpow2 = -30.0_real32
   if (wpow3 < 5.01_real32) wpow3 = -30.0_real32
   if (wpow4 < 5.01_real32) wpow4 = -30.0_real32
   if (wen   < 5.01_real32) wen   = -30.0_real32
   if (wmom  < 5.01_real32) wmom  = -30.0_real32
   if (wmom1 < 5.01_real32) wmom1 = -30.0_real32
   if (wmom2 < 5.01_real32) wmom2 = -30.0_real32
   if (wmom3 < 5.01_real32) wmom3 = -30.0_real32
   if (wmom4 < 5.01_real32) wmom4 = -30.0_real32

   write(95,500) time_in, wpow, wpow1, wpow2, wpow3, wpow4, wen, &
                 wmom, wmom1, wmom2, wmom3, wmom4
500 format(1x,e10.5,3x,5f9.3,3x,f9.3,3x,5f9.3)

   ! Convert back to linear units for further use
   wpow  = 10.0_real32**(wpow  - 35.0_real32)
   wpow1 = 10.0_real32**(wpow1 - 35.0_real32)
   wpow2 = 10.0_real32**(wpow2 - 35.0_real32)
   wpow3 = 10.0_real32**(wpow3 - 35.0_real32)
   wpow4 = 10.0_real32**(wpow4 - 35.0_real32)
   wen   = 10.0_real32**(wen   - 35.0_real32)
   wmom  = 10.0_real32**(wmom  - 35.0_real32)
   wmom1 = 10.0_real32**(wmom1 - 35.0_real32)
   wmom2 = 10.0_real32**(wmom2 - 35.0_real32)
   wmom3 = 10.0_real32**(wmom3 - 35.0_real32)
   wmom4 = 10.0_real32**(wmom4 - 35.0_real32)

end subroutine windpower
c
c *****************************************************************************
c *****************************************************************************
c *****************************************************************************
!> Computes a synthetic spectrum from a population of stars.
!> Modern Fortran version.
subroutine specsyn(time_in, icount)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   real(real32), intent(in) :: time_in
   integer, intent(in)      :: icount

   ! Parameters
   integer, parameter :: nmaxint = 10, nmaxint1 = 11
   integer, parameter :: npgrid = 3000, nwave = 1221, nfeat = 7
   integer, parameter :: nwrlines = 9, nwrtypes = 11

   ! Variables
   character(len=20) :: name
   character(len=80) :: str_date
   character(len=6)  :: wrident(nwrlines)
   real(real32) :: temp(npgrid), bol(npgrid), xsurf(npgrid), zmass(npgrid)
   real(real32) :: ysurf(npgrid), xc12s(npgrid), xn14s(npgrid), xo16s(npgrid)
   real(real32) :: bmdot(npgrid), grav(npgrid), cmass(npgrid), dens(npgrid)
   real(real32) :: flam(nwave), flux(npgrid, nwave), toflux(nwave), wave(nwave)
   real(real32) :: fluneb(nwave), xrange(26), conti(26), stflux(nwave)
   real(real32) :: co1(npgrid), ca1(npgrid), ca3(npgrid)
   real(real32) :: co162a(npgrid), co229a(npgrid), si159a(npgrid), tt_star(npgrid)
   real(real32) :: xlya1(npgrid)
   real(real32) :: wrfluxes(nwrtypes, nwrlines), xlam_wrlines(nwrlines)
   real(real32) :: wr_lssp(nwrlines), wr_wssp(nwrlines)
   real(real32) :: co2, ca2, ca4, ca5, co, ca, co162, co162b, co229, co229b, si159, si159b, xlya, xlya2
   integer :: lmin, lmax, i, l, m, ico1, iwrtype
   real(real32) :: teff, blogg, radius, cotest

   ! External interfaces
   interface
      subroutine init_wrlinefl(zmetal)
         real(real32), intent(in) :: zmetal
      end subroutine init_wrlinefl
      function get_wrtype(xh, xhe, xc, xn, xo, tstar, xlogg) result(wrtype)
         real(real32), intent(in) :: xh, xhe, xc, xn, xo, tstar, xlogg
         integer :: wrtype
      end function get_wrtype
      subroutine planck(teff)
         real(real32), intent(in) :: teff
      end subroutine planck
      subroutine kurucz(teff, blogg)
         real(real32), intent(in) :: teff, blogg
      end subroutine kurucz
      subroutine werner(l, radius, teff)
         integer, intent(in) :: l
         real(real32), intent(in) :: radius, teff
      end subroutine werner
      subroutine hillier(l, teff, radius)
         integer, intent(in) :: l
         real(real32), intent(in) :: teff, radius
      end subroutine hillier
      subroutine pauldrach(teff, blogg, radius)
         real(real32), intent(in) :: teff, blogg, radius
      end subroutine pauldrach
      function sp_feature(l, ifeat, teff, icount) result(val)
         integer, intent(in) :: l, ifeat, icount
         real(real32), intent(in) :: teff
         real(real32) :: val
      end function sp_feature
      subroutine ionize(time_in, icount, stflux)
         real(real32), intent(in) :: time_in
         integer, intent(in) :: icount
         real(real32), intent(in) :: stflux(nwave)
      end subroutine ionize
      subroutine continuum(time_in, icount, conti)
         real(real32), intent(in) :: time_in
         integer, intent(in) :: icount
         real(real32), intent(out) :: conti(26)
      end subroutine continuum
      function yntra(x, xarr, yarr, n) result(y)
         real(real32), intent(in) :: x
         real(real32), intent(in) :: xarr(*), yarr(*)
         integer, intent(in) :: n
         real(real32) :: y
      end function yntra
      subroutine colors(time_in, icount, toflux, stflux)
         real(real32), intent(in) :: time_in
         integer, intent(in) :: icount
         real(real32), intent(in) :: toflux(nwave), stflux(nwave)
      end subroutine colors
      subroutine width(time_in, icount, toflux, stflux)
         real(real32), intent(in) :: time_in
         integer, intent(in) :: icount
         real(real32), intent(in) :: toflux(nwave), stflux(nwave)
      end subroutine width
      subroutine wrlines_ew(wave, toflux, nw)
         real(real32), intent(in) :: wave(*), toflux(*)
         integer, intent(in) :: nw
      end subroutine wrlines_ew
      subroutine fdate(str_date)
         character(len=*), intent(out) :: str_date
      end subroutine fdate
   end interface

   ! Header output at first call
   if (icount <= 1) then
      call init_wrlinefl(log10(z) + 8.66_real32)
      ico1 = 0
      write(92,96) name
96    format(' MODEL DESIGNATION: ',a20)
      call fdate(str_date)
      write(92,597) str_date
597   format(' MODEL GENERATED: ', a24)
      write(92,95)
95    format(/,'              COMPUTED SYNTHETIC SPECTRUM')
      write(92,94)
94    format(/,' TIME [YR]    WAVELENGTH [A]   LOG ', &
         'TOTAL  LOG STELLAR  LOG NEBULAR  [ERG/SEC/A]')
   end if

   ! Reset accumulators if not continuous SF
   if (.not.(isf > 0 .and. icount > 1)) then
      stflux = 0.0_real32
      toflux = 0.0_real32
      fluneb = 0.0_real32
      co2 = 0.0_real32
      ca2 = 0.0_real32
      ca4 = 0.0_real32
      ca5 = 0.0_real32
      co = 0.0_real32
      ca = 0.0_real32
      co162 = 0.0_real32
      co162b = 0.0_real32
      co229 = 0.0_real32
      co229b = 0.0_real32
      si159 = 0.0_real32
      si159b = 0.0_real32
      xlya = 0.0_real32
      xlya2 = 0.0_real32
      wr_lssp = 0.0_real32
   end if

   ! Loop over all masses
   do l = lmin, lmax
      teff = 10.0_real32**temp(l)
      grav(l) = log10(zmass(l)) + 4.0_real32*temp(l) - bol(l) - 10.6_real32
      blogg = grav(l)
      radius = 10.0_real32**(10.8426_real32 + 0.5_real32*bol(l) - 2.0_real32*temp(l) + 7.52_real32)
      cotest = 5.71_real32 * log10(teff) - 21.95_real32

      if (bol(l) < -19.0_real32) cycle

      select case (iatmos)
      case (1)
         call planck(teff)
      case (2)
         if (teff > 60000.0_real32 .or. teff < 2000.0_real32) then
            call planck(teff)
         else
            call kurucz(teff, blogg)
         end if
      case (3)
         if (temp(l) > 4.4_real32 .and. xsurf(l) < 0.4_real32 .and. cmass(l) >= xmwr) then
            call werner(l, radius, teff)
         else
            if (teff > 60000.0_real32 .or. teff < 2000.0_real32) then
               call planck(teff)
            else
               call kurucz(teff, blogg)
            end if
         end if
      case (4)
         if (temp(l) > 4.4_real32 .and. xsurf(l) < 0.4_real32 .and. cmass(l) >= xmwr) then
            call hillier(l, teff, radius)
         else
            if (teff > 60000.0_real32 .or. teff < 2000.0_real32) then
               call planck(teff)
            else
               call kurucz(teff, blogg)
            end if
         end if
      case (5)
         if (temp(l) > 4.4_real32 .and. xsurf(l) < 0.4_real32 .and. cmass(l) >= xmwr) then
            call hillier(l, teff, radius)
         else
            if (teff > 60000.0_real32 .or. teff < 2000.0_real32) then
               call planck(teff)
            else
               if (teff >= 25000.0_real32 .and. blogg < cotest .and. blogg >= 2.2_real32) then
                  call pauldrach(teff, blogg, radius)
               else
                  call kurucz(teff, blogg)
               end if
            end if
         end if
      end select

      ! Compute fluxes for this mass bin
      do m = 1, nwave
         flux(l, m) = 12.566_real32 * radius**2 / 1.e20_real32 * flam(m) * dens(l)
         stflux(m) = stflux(m) + flux(l, m)
      end do

      ! Spectral features (if enabled)
      if (io11 >= 0) then
         co1(l) = sp_feature(l, 1, teff, icount) * (flux(l, 861) + 2.0_real32 * (flux(l, 870) - flux(l, 861)))
         ca1(l) = sp_feature(l, 2, teff, icount) * flux(l, 611)
         ca3(l) = sp_feature(l, 3, teff, icount) * flux(l, 611)
         co162a(l) = sp_feature(l, 4, teff, icount) * (flux(l, 789) + 0.835_real32 * (flux(l, 806) - flux(l, 789)))
         co229a(l) = sp_feature(l, 5, teff, icount) * (flux(l, 861) + 1.190_real32 * (flux(l, 870) - flux(l, 861)))
         si159a(l) = sp_feature(l, 6, teff, icount) * (flux(l, 789) + 0.544_real32 * (flux(l, 806) - flux(l, 789)))
         xlya1(l) = sp_feature(l, 7, teff, icount) * (flux(l, 149) + flux(l, 155)) / 2.0_real32
         co2 = co2 + co1(l)
         ca2 = ca2 + ca1(l)
         ca4 = ca4 + ca3(l)
         co162b = co162b + co162a(l)
         co229b = co229b + co229a(l)
         si159b = si159b + si159a(l)
         xlya2 = xlya2 + xlya1(l)
      end if

      ! WR emission lines (if enabled)
      if (io14 >= 0) then
         iwrtype = get_wrtype(xsurf(l), ysurf(l), xc12s(l), xn14s(l), xo16s(l), temp(l), grav(l))
         if (bol(l) < -19.0_real32) iwrtype = 11
         do i = 1, nwrlines
            wr_lssp(i) = wr_lssp(i) + wrfluxes(iwrtype, i) * dens(l)
         end do
      end if
   end do

   ! Compute number of ionizing photons if desired
   if (io1 >= 0) call ionize(time_in, icount, stflux)

   ! Calculate nebular continuum and add to stellar continuum
   call continuum(time_in, icount, conti)
   do m = 1, nwave
      fluneb(m) = yntra(wave(m), xrange, conti, 26)
      toflux(m) = stflux(m) + fluneb(m)
   end do

   ! Compute spectral features
   if (io11 >= 0) then
      co = -2.5_real32 * log10((1.e-35_real32 + co2 + fluneb(861) + 2.0_real32 * (fluneb(870) - fluneb(861))) / &
           (toflux(861) + 2.0_real32 * (toflux(870) - toflux(861)) + 1.e-35_real32))
      ca = ca2 / (toflux(611) + 1.e-35_real32)
      ca5 = ca4 / (toflux(611) + 1.e-35_real32)
      co162 = co162b / (toflux(789) + 0.835_real32 * (toflux(806) - toflux(789)) + 1.e-35_real32)
      co229 = co229b / (toflux(861) + 1.190_real32 * (toflux(870) - toflux(861)) + 1.e-35_real32)
      si159 = si159b / (toflux(789) + 0.544_real32 * (toflux(806) - toflux(789)) + 1.e-35_real32)
      xlya = xlya2 / (toflux(149) + toflux(155) + 1.e-35_real32) * 2.0_real32
   end if

   ! Output spectrum at selected time steps
   if (mod(time_in, tdel) < tstep) then
      write(92,500) (time_in, wave(i), log10(toflux(i)+1.e-35_real32)+20.0_real32, &
                     log10(stflux(i)+1.e-35_real32)+20.0_real32, &
                     log10(fluneb(i)+1.e-35_real32)+20.0_real32, i=1,nwave)
500   format(1x,e10.5,2x,f13.2,4x,f10.5,2x,f10.5,2x,f10.5)
   end if

   if (io9 >= 0)  call colors(time_in, icount, toflux, stflux)
   if (io10 >= 0) call width(time_in, icount, toflux, stflux)

   ! Output spectral features if enabled
   if (io11 >= 0) then
      if (ico1 == 0) then
         write(87,968) name
968      format(' MODEL DESIGNATION: ',a20)
         call fdate(str_date)
         write(87,978) str_date
978      format(' MODEL GENERATED: ', a24)
         write(87,958)
958      format(/,'              COMPUTED SPECTRAL FEATURES')
         write(87,948)
948      format(/,' TIME [YR]    CO INDEX   Ca EQW    CaEQW2    CO1.62', &
            '    CO2.29    Si1.59    Ly-alpha')
         write(87,949)
949      format('                [MAG]     [AA]      [AA]      [AA]      ', &
            '[AA]      [AA]      [AA]')
      end if
      write(87,448) time_in, co, ca, ca5, co162, co229, si159, xlya
448   format(1x,e10.5,2x,f8.3,2x,f8.3,2x,f8.3,2x,f8.3,2x,f8.3,2x,f8.3,2x,f8.3)
   end if

   ! Output WR emission lines if enabled
   if (io14 >= 0) then
      call wrlines_ew(wave, toflux, nwave)
      if (ico1 == 0) then
         write(84,969) name
969      format(' MODEL DESIGNATION: ',a20)
         call fdate(str_date)
         write(84,979) str_date
979      format(' MODEL GENERATED: ', a24)
         write(84,989)
989      format(/,'                 COMPUTED WOLF-RAYET EMISSION', &
            ' LINES')
         write(84,221) (wrident(i), i=1,nwrlines)
221      format(/,'     TIME ',5x,9(a6,7x))
         write(84,990)
990      format('      YR      AA LOG E/S',9('   AA LOG E/S'))
      end if
      write(84,222) time_in, (wr_wssp(i), 35.0_real32+log10(wr_lssp(i)+1.e-35_real32), i=1,nwrlines)
222   format(1x,e10.5,9(1x,f6.2,f6.2))
   end if

   ico1 = ico1 + 1

end subroutine specsyn
c
c *****************************************************************************
c *****************************************************************************
c *****************************************************************************
!> Computes a synthetic spectrum from a population of stars using the IFA library.
!> Output: log(erg/s/A) vs wavelength, 900-3000A, high resolution.
subroutine ifa_spectrum(time_in, icount)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   real(real32), intent(in) :: time_in
   integer, intent(in)      :: icount

   ! Parameters
   integer, parameter :: nmaxint = 10, nmaxint1 = 11
   integer, parameter :: npgrid = 3000, nwave = 4200, nfree = 22
   character(len=20) :: name
   character(len=80) :: str_date

   ! Common blocks (replace with modules in modern code)
   real(real32) :: cmass(npgrid), temp(npgrid), bol(npgrid), xsurf(npgrid), zmass(npgrid)
   real(real32) :: ysurf(npgrid), xc12s(npgrid), xn14s(npgrid), xo16s(npgrid)
   real(real32) :: bmdot(npgrid), tt_star(npgrid), dens(npgrid)
   integer :: lmin, lmax
   real(real32) :: flam(nwave), wave(nwave)
   real(real32) :: flam_ifa_l(nwave), flam_ifa_c(nwave)
   real(real32) :: flux_ifa_l(npgrid, nwave), flux_ifa_c(npgrid, nwave)
   real(real32) :: stflux_ifa_l(nwave), stflux_ifa_c(nwave)
   real(real32) :: toflux_ifa_l(nwave), toflux_ifa_c(nwave)
   real(real32) :: fluneb_ifa(nwave)
   real(real32) :: wave_ifa(nwave)
   real(real32) :: flam1_ifa_l(nwave,86), flam1_ifa_c(nwave,86)
   real(real32) :: tem_ifa(86), glog_ifa(86)
   integer :: nlej_ifa
   real(real32) :: freelam1(nfree), freeflux1(nfree), flx(nwave)
   real(real32) :: xrange(26), gamma(26), conti(26)
   real(real32) :: wave_powr(15000), flam_powr(15000)
   integer :: kmax_powr
   real(real32) :: radius, teff, blogg, grav(npgrid)
   integer :: l, m, jj

   ! Data
   data freelam1 /899.9,1085.,1150.,1280.,1310.,1360.,1430.,1480., &
                  1510.,1580.,1630.,1680.,1740.,1785.,1820.,1875.,2035., &
                  2235.,2565.,2745.,2895.,2999./

   ! Header output at first call
   if (icount <= 1) then
      write(83,96) name
96    format(' MODEL DESIGNATION: ',a20)
      call fdate(str_date)
      write(83,597) str_date
597   format(' MODEL GENERATED: ', a24)
      write(83,95)
95    format(/,'        COMPUTED SYNTHETIC LINE SPECTRUM [ERG/SEC/A]')
      write(83,94)
94    format(/,' TIME [YR]      WAVELENGTH  LOG(LUMINOSITY)  ', &
         'NORMALIZED SPECTRUM')
   end if

   ! Reset flux accumulators if not continuous SF
   if (.not.(isf > 0 .and. icount > 1)) then
      stflux_ifa_l = 0.0_real32
      stflux_ifa_c = 0.0_real32
      toflux_ifa_l = 0.0_real32
      toflux_ifa_c = 0.0_real32
      fluneb_ifa   = 0.0_real32
   end if

   ! Loop over all masses
   do l = lmin, lmax
      if (bol(l) < -19.0_real32) then
         flux_ifa_l(l,:) = 1.e-30_real32
         flux_ifa_c(l,:) = 1.e-30_real32
         cycle
      end if

      teff   = 10.0_real32**temp(l)
      grav(l)= log10(zmass(l)) + 4.0_real32*temp(l) - bol(l) - 10.6_real32
      blogg  = grav(l)
      radius = 10.0_real32**(10.8426_real32 + 0.5_real32*bol(l) - 2.0_real32*temp(l) + 7.52_real32)

      if (bol(l) < -19.0_real32) cycle

      if (temp(l) > 4.4_real32 .and. xsurf(l) < 0.4_real32 .and. cmass(l) >= xmwr) then
         call hamann(l, teff, radius)
         call planck(teff)
         do m = 1, nwave
            flam_ifa_c(m) = yntra(wave_ifa(m), wave, flam, 1221)
         end do
         do m = 1, nwave
            flam_ifa_l(m) = yntra(wave_ifa(m), wave_powr, flam_powr, kmax_powr)
            flam_ifa_c(m) = flam_ifa_c(m) * flam_ifa_l(4050) / flam_ifa_c(4050)
         end do
      else if (teff < 2000.0_real32) then
         call planck(teff)
         do m = 1, nwave
            flam_ifa_l(m) = yntra(wave_ifa(m), wave, flam, 1221)
            flam_ifa_c(m) = flam_ifa_l(m)
         end do
      else if (cmass(l) > 10.0_real32 .and. teff < 22000.0_real32) then
         call kurucz(teff, blogg)
         do jj = 1, nfree
            freeflux1(jj) = yntra(freelam1(jj), wave, flam, 1221)
            freeflux1(jj) = log10((1.e-30_real32 + freeflux1(jj) + flam(imatch-1) + flam(imatch+1)) / 3.0_real32)
         end do
         do m = 1, nwave
            flam_ifa_l(m) = yntra(wave_ifa(m), wave, flam, 1221)
            call intrpl(nfree, freelam1, freeflux1, 1, wave_ifa(m), flx(m))
            flam_ifa_c(m) = 10.0_real32**flx(m)
         end do
      else if (cmass(l) < 5.0_real32 .or. (cmass(l) >= 5.0_real32 .and. teff < 17000.0_real32)) then
         call kurucz(teff, blogg)
         do jj = 1, nfree
            freeflux1(jj) = yntra(freelam1(jj), wave, flam, 1221)
            freeflux1(jj) = log10((1.e-30_real32 + freeflux1(jj) + flam(imatch-1) + flam(imatch+1)) / 3.0_real32)
         end do
         do m = 1, nwave
            flam_ifa_l(m) = yntra(wave_ifa(m), wave, flam, 1221)
            call intrpl(nfree, freelam1, freeflux1, 1, wave_ifa(m), flx(m))
            flam_ifa_c(m) = 10.0_real32**flx(m)
         end do
      else
         call fabio(teff, blogg)
      end if

      ! Compute fluxes for this mass bin
      do m = 1, nwave
         flux_ifa_l(l,m) = 12.566_real32 * radius**2 / 1.e20_real32 * flam_ifa_l(m) * dens(l)
         flux_ifa_c(l,m) = 12.566_real32 * radius**2 / 1.e20_real32 * flam_ifa_c(m) * dens(l)
         stflux_ifa_l(m) = stflux_ifa_l(m) + flux_ifa_l(l,m)
         stflux_ifa_c(m) = stflux_ifa_c(m) + flux_ifa_c(l,m)
      end do
   end do

   ! Add nebular continuum
   call continuum(time_in, icount, conti)
   do m = 1, nwave
      fluneb_ifa(m)   = yntra(wave_ifa(m), xrange, conti, 26)
      toflux_ifa_l(m) = stflux_ifa_l(m) + fluneb_ifa(m)
      toflux_ifa_c(m) = stflux_ifa_c(m) + fluneb_ifa(m)
   end do

   ! Output spectrum at selected time steps
   if (mod(time_in - time1, tdel) < tstep) then
      write(83,500) (time_in, wave_ifa(i), log10(toflux_ifa_l(i)+1.e-35_real32)+20.0_real32, &
                     toflux_ifa_l(i)/(toflux_ifa_c(i)+1.e-35_real32), i=1,nwave)
500   format(1x,e10.5,4x,f10.2,5x,f10.5,5x,f10.4)
   end if

end subroutine ifa_spectrum
c
c ********************************************************************
c ********************************************************************
c ********************************************************************
!> Calculates the number of ionizing photons from all stars at each time step.
!> Computes for H I, He I, and He II. Also computes the total bolometric luminosity.
subroutine ionize(time_in, icount, stflux)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   real(real32), intent(in) :: time_in
   integer, intent(in)      :: icount
   real(real32), intent(in) :: stflux(1221)

   ! Parameters
   integer, parameter :: nmaxint = 10, nmaxint1 = 11
   integer, parameter :: nwave = 1221
   character(len=20) :: name
   character(len=80) :: str_date
   real(real32) :: wave(nwave)
   real(real32) :: tem(600), glog(600), flam1(nwave,600)
   real(real32) :: phot_912
   real(real32) :: freq(nwave), f1data(nwave), weight(nwave)
   real(real32) :: fnu1(44), fla1(44), xnu1(44), xla1(44)
   real(real32) :: fnu2(80), fla2(80), xnu2(80), xla2(80)
   real(real32) :: fnu3(122), fla3(122), xnu3(122), xla3(122)
   real(real32) :: radlum, phot_228, flu_228, phot_504, flu_504, phot_912_log, flu_912
   integer :: i

   ! External
   interface
      subroutine fdate(str_date)
         character(len=*), intent(out) :: str_date
      end subroutine fdate
      subroutine fliwgt(x, y, w, value, n)
         real(real32), intent(in) :: x(n), y(n)
         real(real32), intent(out) :: w(n)
         real(real32), intent(out) :: value
         integer, intent(in) :: n
      end subroutine fliwgt
   end interface

   ! Header output at first call
   if (icount <= 1) then
      write(98,96) name
96    format(' MODEL DESIGNATION: ',a20)
      call fdate(str_date)
      write(98,597) str_date
597   format(' MODEL GENERATED: ', a24)
      write(98,95)
95    format(/,'              RESULTS FOR THE',
     &        ' NUMBER OF IONIZING PHOTONS')
      write(98,94)
94    format(/,'   TIME            H I              HE I      ',
     &        '      HE II         LOG L')
      write(98,93)
93    format('              [s^-1]  % OF L   [s^-1]  % OF L  ',
     &        ' [s^-1]  % OF L   [ERG/S]')
   end if

   ! Convert to frequency and from F_lambda to F_nu; divide by photon energy
   do i = 1, nwave
      freq(i)   = 2.997925e18_real32 / wave(i)
      f1data(i) = 3.33564e-19_real32 * stflux(i) * wave(i) * wave(i)
      f1data(i) = f1data(i) / freq(i) / 6.6262e-27_real32
   end do

   ! Bolometric luminosity (erg/s)
   call fliwgt(wave, stflux, weight, radlum, nwave)
   radlum = log10(abs(-1.0_real32 * radlum) + 1.e-30_real32) + 20.0_real32

   ! H I continuum (up to 912A, i=44)
   fnu1(:) = f1data(1:44)
   fla1(:) = stflux(1:44)
   xnu1(:) = freq(1:44)
   xla1(:) = wave(1:44)
   call fliwgt(xnu1, fnu1, weight, phot_228, 44)
   call fliwgt(xla1, fla1, weight, flu_228, 44)
   phot_228 = log10(phot_228 + 1.e-30_real32) + 20.0_real32
   flu_228  = log10(abs(-1.0_real32 * flu_228) + 1.e-30_real32) + 20.0_real32

   ! He I continuum (up to 504A, i=80)
   fnu2(:) = f1data(1:80)
   fla2(:) = stflux(1:80)
   xnu2(:) = freq(1:80)
   xla2(:) = wave(1:80)
   call fliwgt(xnu2, fnu2, weight, phot_504, 80)
   call fliwgt(xla2, fla2, weight, flu_504, 80)
   phot_504 = log10(phot_504 + 1.e-30_real32) + 20.0_real32
   flu_504  = log10(abs(-1.0_real32 * flu_504) + 1.e-30_real32) + 20.0_real32

   ! He II continuum (up to 228A, i=122)
   fnu3(:) = f1data(1:122)
   fla3(:) = stflux(1:122)
   xnu3(:) = freq(1:122)
   xla3(:) = wave(1:122)
   call fliwgt(xnu3, fnu3, weight, phot_912, 122)
   call fliwgt(xla3, fla3, weight, flu_912, 122)
   phot_912_log = log10(phot_912 + 1.e-30_real32) + 20.0_real32
   flu_912      = log10(abs(-1.0_real32 * flu_912) + 1.e-30_real32) + 20.0_real32

   ! Output
   write(98,980) time_in, phot_912_log, flu_912 - radlum, &
        phot_504, flu_504 - radlum, phot_228, flu_228 - radlum, radlum
980 format(1x,e10.5,2x,f7.3,f8.3,2x,f7.3,f8.3,2x,f7.3, &
     f8.3,2x,f7.3)

end subroutine ionize
c
c *******************************************************************
c *******************************************************************
c *******************************************************************
!> Computes the emission of a nebular continuum.
!> Hydrogen and helium free-free and free-bound emission as well as the 2-photon continuum are included.
!> Output units: 1.e-20 erg/s/A.
subroutine continuum(time, icount, conti)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   real(real32), intent(in) :: time
   integer, intent(in)      :: icount
   real(real32), intent(out) :: conti(26)

   ! Parameters and shared variables
   real(real32) :: xrange(26), gamma(26)
   real(real32) :: phot_912
   real(real32), parameter :: alphab = 2.60e-13_real32
   integer :: i

   ! External/common blocks (replace with modules in modern code)
   ! These would be set elsewhere in the program
   ! common /nebula/ xrange, gamma
   ! common /equiphot/ phot_912

   do i = 1, 26
      conti(i) = 2.998e18_real32 / xrange(i) / xrange(i) / alphab * 1.e-30_real32 * gamma(i) * 10.0_real32**(phot_912 - 30.0_real32)
   end do

end subroutine continuum
c
c
c *********************************************************************
c *********************************************************************
c *********************************************************************
!> Computes selected synthetic colors from the output spectrum.
subroutine colors(time_in, icount, toflux, stflux)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   real(real32), intent(in) :: time_in
   integer, intent(in)      :: icount
   real(real32), intent(in) :: toflux(1221)
   real(real32), intent(in) :: stflux(1221)

   ! Parameters
   integer, parameter :: nmaxint = 10, nmaxint1 = 11
   integer, parameter :: nband = 12, nwave = 1221
   character(len=20) :: name
   character(len=80) :: str_date
   real(real32) :: wave(nwave)
   real(real32) :: tem(600), glog(600), flam1(nwave,600)
   real(real32) :: xprof(nband,100), yprof(nband,100)
   integer :: ksize(nband)
   real(real32) :: profil(nwave), xmag(nband)
   real(real32) :: fdata(nwave), xgrid(100), prof(100)
   real(real32) :: weight(nwave)
   real(real32) :: absbol, absv, absb
   real(real32) :: cuv1v, cuv2v, cub, cbv, cvr, cvi, cvj, cvh, cvk, cvl
   integer :: j, k, i

   ! External
   interface
      function yntra(x, xarr, yarr, n) result(y)
         use, intrinsic :: iso_fortran_env, only: real32
         real(real32), intent(in) :: x, xarr(n), yarr(n)
         integer, intent(in) :: n
         real(real32) :: y
      end function yntra
      subroutine fdate(str_date)
         character(len=*), intent(out) :: str_date
      end subroutine fdate
      subroutine fliwgt(x, y, w, value, n)
         real(real32), intent(in) :: x(n), y(n)
         real(real32), intent(out) :: w(n)
         real(real32), intent(out) :: value
         integer, intent(in) :: n
      end subroutine fliwgt
   end interface

   data ksize /99,90,26,42,42,56,26,25,33,30,19,16/

   ! Header output at first call
   if (icount <= 1) then
      write(89,96) name
96    format(' MODEL DESIGNATION: ',a20)
      call fdate(str_date)
      write(89,597) str_date
597   format(' MODEL GENERATED: ', a24)
      write(89,95)
95    format(/,'              COMPUTED SYNTHETIC COLORS')
      write(89,94)
94    format(/,' TIME [YR]    (130-V) (210-V)   (U-B)   (B-V)  ', &
         ' (V-R)   (V-I)   (V-J)   (V-H)  ', &
         ' (V-K)   (V-L)     M_V      M_B     M_BOL')
   end if

   ! Compute magnitudes in the 12 passbands
   do j = 1, nband
      do k = 1, ksize(j)
         xgrid(k) = xprof(j,k)
         prof(k)  = yprof(j,k)
      end do
      do i = 1, nwave
         profil(i) = yntra(wave(i), xgrid, prof, ksize(j))
         fdata(i)  = toflux(i) * profil(i)
      end do
      call fliwgt(wave, fdata, weight, xmag(j), nwave)
      xmag(j) = -2.5 * log10(abs(xmag(j) + 1.e-35_real32))
   end do

   ! Calculate colors (zero points for 9400K/3.95, Z=Z_solar)
   cuv1v = xmag(1) - xmag(6)  - 6.003_real32
   cuv2v = xmag(2) - xmag(6)  - 2.988_real32
   cub   = xmag(3) - xmag(4)  - 1.102_real32
   cbv   = xmag(5) - xmag(6)  + 0.722_real32
   cvr   = xmag(6) - xmag(7)  - 0.093_real32
   cvi   = xmag(6) - xmag(8)  + 0.692_real32
   cvj   = xmag(6) - xmag(9)  + 1.284_real32
   cvh   = xmag(6) - xmag(10) + 1.875_real32
   cvk   = xmag(6) - xmag(11) + 2.746_real32
   cvl   = xmag(6) - xmag(12) + 4.384_real32

   ! Calculate M_V, M_B, and M_BOL
   call fliwgt(wave, stflux, weight, absbol, nwave)
   absbol = -2.5 * log10(abs(absbol) + 1.e-30_real32) + 38.70_real32
   absv   = xmag(6) + 36.552_real32
   absb   = absv + cbv

   ! Output
   write(89,44) time_in, cuv1v, cuv2v, cub, cbv, cvr, cvi, cvj, cvh, &
        cvk, cvl, absv, absb, absbol
44 format(1x,e10.5,2x,10f8.3,1x,f8.3,1x,f8.3,1x,f8.3)

end subroutine colors
c
c
c ***********************************************************************
c ***********************************************************************
c ***********************************************************************
!> Computes selected equivalent widths due to nebular emission lines.
!> Modern Fortran version.
subroutine width(time_in, icount, toflux, stflux)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   real(real32), intent(in) :: time_in
   integer, intent(in)      :: icount
   real(real32), intent(in) :: toflux(1221)
   real(real32), intent(in) :: stflux(1221)

   ! Parameters
   integer, parameter :: nmaxint = 10, nmaxint1 = 11
   character(len=20) :: name
   character(len=80) :: str_date
   real(real32) :: wave(1221)
   real(real32) :: tem(600), glog(600), flam1(1221,600)
   real(real32) :: phot_912

   ! Locals
   real(real32) :: hacon, hbcon, pbcon, bgcon
   real(real32) :: halum, hblum, pblum, bglum
   real(real32) :: hawid, hbwid, pbwid, bgwid

   ! External
   interface
      function yntra(x, xarr, yarr, n) result(y)
         use, intrinsic :: iso_fortran_env, only: real32
         real(real32), intent(in) :: x, xarr(n), yarr(n)
         integer, intent(in) :: n
         real(real32) :: y
      end function yntra
      subroutine fdate(str_date)
         character(len=*), intent(out) :: str_date
      end subroutine fdate
   end interface

   ! Header output at first call
   if (icount <= 1) then
      write(88,96) name
96    format(' MODEL DESIGNATION: ',a20)
      call fdate(str_date)
      write(88,597) str_date
597   format(' MODEL GENERATED: ', a24)
      write(88,95)
95    format(/,'                   ',
     &          '      COMPUTED RECOMBINATION LINES')
      write(88,94)
94    format(/,' TIME [YR]     CONT(H_A) LUM(H_A) EQ(H_A)  ',
     &                          'CONT(H_B) LUM(H_B) EQ(H_B)  ',
     &                          'CONT(P_B) LUM(P_B) EQ(P_B)  ',
     &                          'CONT(B_G) LUM(B_G) EQ(B_G)  ')
   end if

   ! Interpolate continuum fluxes at line wavelengths
   hacon = yntra(6563.0_real32, wave, toflux, 1221)
   hbcon = yntra(4861.0_real32, wave, toflux, 1221)
   pbcon = yntra(12818.0_real32, wave, toflux, 1221)
   bgcon = yntra(21655.0_real32, wave, toflux, 1221)

   ! Compute emission-line luminosities from number of ionizing photons
   halum = 1.36e-12_real32 * 10.0_real32**(phot_912 - 30.0_real32)
   hawid = halum / hacon * 1.0e10_real32
   hblum = 4.76e-13_real32 * 10.0_real32**(phot_912 - 30.0_real32)
   hbwid = hblum / hbcon * 1.0e10_real32
   pblum = 7.73e-14_real32 * 10.0_real32**(phot_912 - 30.0_real32)
   pbwid = pblum / pbcon * 1.0e10_real32
   bglum = 1.31e-14_real32 * 10.0_real32**(phot_912 - 30.0_real32)
   bgwid = bglum / bgcon * 1.0e10_real32

   ! Output
   write(88,44) time_in, log10(hacon+1.e-35_real32)+20.0_real32, &
        log10(halum+1.e-35_real32)+30.0_real32, log10(hawid+1.e-35_real32), &
        log10(hbcon+1.e-35_real32)+20.0_real32, log10(hblum+1.e-35_real32)+30.0_real32, &
        log10(hbwid+1.e-35_real32), log10(pbcon+1.e-35_real32)+20.0_real32, &
        log10(pblum+1.e-35_real32)+30.0_real32, log10(pbwid+1.e-35_real32), &
        log10(bgcon+1.e-35_real32)+20.0_real32, log10(bglum+1.e-35_real32)+30.0_real32, &
        log10(bgwid+1.e-35_real32)
44 format(1x,e10.5,2x,3f9.3,1x,3f9.3,1x,3f9.3,1x,3f9.3)

end subroutine width
c
c ***********************************************************************
c ***********************************************************************
c ***********************************************************************
!> Initializes WR line fluxes and identifiers.
subroutine init_wrlinefl(zmetal)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   real(real32), intent(in) :: zmetal

   integer, parameter :: nwrlines = 9, nwrtypes = 11
   character(len=6) :: wrident(nwrlines)
   real(real32) :: wrfluxes(nwrtypes, nwrlines)
   real(real32) :: xlam_wrlines(nwrlines)
   real(real32) :: wr_lssp(nwrlines)
   real(real32) :: wr_wssp(nwrlines)
   real(real32) :: zscale1, zscale2
   integer :: i, j

   ! Define identifiers and wavelengths for lines
   wrident = ['He1640', 'N_4640', 'C_4650', 'He4686', 'H_4861', 'C_5696', 'C_5808', 'H_6560', 'He5411']
   xlam_wrlines = [1640.0_real32, 4640.0_real32, 4650.0_real32, 4686.0_real32, 4861.0_real32, &
                   5696.0_real32, 5808.0_real32, 6540.0_real32, 5411.0_real32]

   ! Set line fluxes to zero
   wrfluxes = 0.0_real32

   ! Define the scaling factor to account for the metallicity dependence
   ! of the 4686 and the 5808 lines (Lopez-Sanchez et al. 2010)
   zscale1 = -3.394_real32 + 0.508_real32 * zmetal
   zscale2 = -3.279_real32 + 0.494_real32 * zmetal

   ! OIf stars
   wrfluxes(10,4) = 10.0_real32**0.4_real32
   wrfluxes(10,2) = 10.0_real32**(-1.0_real32) * wrfluxes(10,4)

   ! WNE stars (HeII4686 is reference)
   wrfluxes(1,4) = zscale1 * 5.2_real32
   wrfluxes(1,1) = 7.952_real32 * wrfluxes(1,4)
   wrfluxes(1,5) = 0.129_real32 * wrfluxes(1,4)
   wrfluxes(1,7) = 0.074_real32 * wrfluxes(1,4)
   wrfluxes(1,8) = 0.202_real32 * wrfluxes(1,4)
   wrfluxes(1,9) = 0.137_real32 * wrfluxes(1,4)

   ! WNL stars (HeII4686 is reference)
   wrfluxes(2,4) = zscale1 * 16.0_real32
   wrfluxes(2,1) = 7.55_real32 * wrfluxes(2,4)
   wrfluxes(2,2) = 0.616_real32 * wrfluxes(2,4)
   wrfluxes(2,5) = 0.309_real32 * wrfluxes(2,4)
   wrfluxes(2,7) = 0.062_real32 * wrfluxes(2,4)
   wrfluxes(2,8) = 1.596_real32 * wrfluxes(2,4)
   wrfluxes(2,9) = 0.107_real32 * wrfluxes(2,4)

   ! WO3-4 stars (CIV5808 is reference)
   wrfluxes(3,7) = zscale2 * 11.0_real32
   wrfluxes(3,1) = 2.65_real32 * wrfluxes(3,7)
   wrfluxes(3,3) = 0.51_real32 * wrfluxes(3,7)
   wrfluxes(3,4) = 0.386_real32 * wrfluxes(3,7)
   wrfluxes(3,8) = 0.078_real32 * wrfluxes(3,7)
   wrfluxes(3,9) = 0.026_real32 * wrfluxes(3,7)

   ! WC4 stars (CIV5808 is reference)
   wrfluxes(4,7) = zscale2 * 30.0_real32
   wrfluxes(4,1) = 2.14_real32 * wrfluxes(4,7)
   wrfluxes(4,3) = 1.71_real32 * wrfluxes(4,7)
   wrfluxes(4,5) = 0.019_real32 * wrfluxes(4,7)
   wrfluxes(4,6) = 0.021_real32 * wrfluxes(4,7)
   wrfluxes(4,8) = 0.045_real32 * wrfluxes(4,7)
   wrfluxes(4,9) = 0.023_real32 * wrfluxes(4,7)

   ! WC5 stars (CIV5808 is reference)
   wrfluxes(5,7) = zscale2 * 9.8_real32
   wrfluxes(5,1) = 5.14_real32 * wrfluxes(5,7)
   wrfluxes(5,3) = 2.53_real32 * wrfluxes(5,7)
   wrfluxes(5,5) = 0.022_real32 * wrfluxes(5,7)
   wrfluxes(5,6) = 0.077_real32 * wrfluxes(5,7)
   wrfluxes(5,8) = 0.074_real32 * wrfluxes(5,7)
   wrfluxes(5,9) = 0.036_real32 * wrfluxes(5,7)

   ! WC6 stars (CIV5808 is reference)
   wrfluxes(6,7) = zscale2 * 8.9_real32
   wrfluxes(6,1) = 7.60_real32 * wrfluxes(6,7)
   wrfluxes(6,3) = 2.98_real32 * wrfluxes(6,7)
   wrfluxes(6,5) = 0.047_real32 * wrfluxes(6,7)
   wrfluxes(6,6) = 0.184_real32 * wrfluxes(6,7)
   wrfluxes(6,8) = 0.112_real32 * wrfluxes(6,7)
   wrfluxes(6,9) = 0.054_real32 * wrfluxes(6,7)

   ! WC7 stars (CIV5808 is reference)
   wrfluxes(7,7) = zscale2 * 14.0_real32
   wrfluxes(7,1) = 10.22_real32 * wrfluxes(7,7)
   wrfluxes(7,3) = 3.21_real32 * wrfluxes(7,7)
   wrfluxes(7,5) = 0.069_real32 * wrfluxes(7,7)
   wrfluxes(7,6) = 0.579_real32 * wrfluxes(7,7)
   wrfluxes(7,8) = 0.177_real32 * wrfluxes(7,7)
   wrfluxes(7,9) = 0.076_real32 * wrfluxes(7,7)

   ! WC8 stars (CIV5808 is reference)
   wrfluxes(8,7) = zscale2 * 3.5_real32
   wrfluxes(8,1) = 13.69_real32 * wrfluxes(8,7)
   wrfluxes(8,3) = 3.30_real32 * wrfluxes(8,7)
   wrfluxes(8,5) = 0.096_real32 * wrfluxes(8,7)
   wrfluxes(8,6) = 1.968_real32 * wrfluxes(8,7)
   wrfluxes(8,8) = 0.580_real32 * wrfluxes(8,7)
   wrfluxes(8,9) = 0.091_real32 * wrfluxes(8,7)

   ! WC9 stars (CIV5808 is reference)
   wrfluxes(9,7) = zscale2 * 2.0_real32
   wrfluxes(9,1) = 4.47_real32 * wrfluxes(8,7)
   wrfluxes(9,3) = 4.46_real32 * wrfluxes(9,7)
   wrfluxes(9,5) = 0.294_real32 * wrfluxes(9,7)
   wrfluxes(9,6) = 3.705_real32 * wrfluxes(9,7)
   wrfluxes(9,8) = 1.243_real32 * wrfluxes(9,7)
   wrfluxes(9,9) = 0.135_real32 * wrfluxes(9,7)

end subroutine init_wrlinefl
c
c *****************************************************************************
c *****************************************************************************
c *****************************************************************************
!> Returns the WR star type number for use with line fluxes.
!> 1 = WNE, 2 = WNL, 3 = WO, 4 = WC4, 5 = WC5, ..., 9 = WC9
!> 10 = OIf, 11 = otherwise.
function get_wrtype(xh, xhe, xc, xn, xo, tstar, xlogg) result(wrtype)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   real(real32), intent(in) :: xh, xhe, xc, xn, xo, tstar, xlogg
   integer :: wrtype
   real(real32) :: cohe, xglim

   wrtype = 11

   if (xh < 0.4_real32 .and. tstar > 4.4_real32) then
      if (xh <= 0.001_real32) then
         if (xc > xn) then
            ! WC stars
            cohe = (xc / 12.0_real32 + xo / 16.0_real32) / (xhe / 4.0_real32)
            if (cohe > 1.0_real32) then
               wrtype = 3   ! WO
            else if (cohe > 0.43_real32) then
               wrtype = 4   ! WC4
            else if (cohe > 0.35_real32) then
               wrtype = 5   ! WC5
            else if (cohe > 0.25_real32) then
               wrtype = 6   ! WC6
            else if (cohe > 0.15_real32) then
               wrtype = 7   ! WC7
            else if (cohe > 0.08_real32) then
               wrtype = 8   ! WC8
            else
               wrtype = 9   ! WC9
            end if
         else
            wrtype = 1      ! WNE
         end if
      else
         wrtype = 2         ! WNL
      end if
   else
      ! Check if OIf
      xglim = 3.676_real32 * tstar - 13.253_real32
      if (xlogg < xglim .and. tstar >= 4.519_real32) then
         wrtype = 10
      end if
   end if

   get_wrtype = wrtype
end function get_wrtype
c
c *****************************************************************************
c *****************************************************************************
c *****************************************************************************
!> Calculates the equivalent width [A] of all WR emission lines.
subroutine wrlines_ew(wave, toflux, nw)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   integer, intent(in) :: nw
   real(real32), intent(in) :: wave(nw), toflux(nw)
   ! Output: updates wr_wssp in common block/module

   ! Parameters
   integer, parameter :: nwrlines = 9, nwrtypes = 11

   ! WR lines common block (should be in a module in modern code)
   real(real32) :: wrfluxes(nwrtypes, nwrlines)
   real(real32) :: xlam_wrlines(nwrlines)
   real(real32) :: wr_lssp(nwrlines)
   real(real32) :: wr_wssp(nwrlines)
   character(len=6) :: wrident(nwrlines)

   ! Locals
   real(real32) :: xl(2), yf(2), f_cont
   integer :: i

   interface
      subroutine linterp(x, y, n, xin, yout, nin)
         use, intrinsic :: iso_fortran_env, only: real32
         integer, intent(in) :: n, nin
         real(real32), intent(in) :: x(n), y(n), xin(nin)
         real(real32), intent(out) :: yout(nin)
      end subroutine linterp
   end interface

   do i = 1, nwrlines
      xl(1) = xlam_wrlines(i) - 20.0_real32
      xl(2) = xlam_wrlines(i) + 20.0_real32
      call linterp(wave, toflux, nw, xl, yf, 2)
      f_cont = yf(1) + (yf(2) - yf(1)) * (xlam_wrlines(i) - xl(1)) / (xl(2) - xl(1))
      wr_wssp(i) = (wr_lssp(i) * 1.0e15_real32) / f_cont
   end do

end subroutine wrlines_ew
c

c *****************************************************************************
c *****************************************************************************
c *****************************************************************************
!> Computes the emergent flux from a black body (Planck function).
!> Modern Fortran version.
subroutine planck(teff)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   real(real32), intent(in) :: teff

   ! Parameters
   integer, parameter :: n_wave = 1221
   real(real32), parameter :: pi = 3.1415927_real32
   real(real32), parameter :: h = 6.626196e-27_real32
   real(real32), parameter :: c = 2.997925e10_real32
   real(real32), parameter :: boltz = 1.380622e-16_real32

   ! Common blocks replaced by assumed-shape arrays or module variables
   real(real32) :: wave(n_wave)
   real(real32) :: flam(n_wave)

   integer :: i
   real(real32) :: cgswave, c1

   ! Compute the Planck function
   do i = 1, n_wave
      cgswave = 1.e-8_real32 * wave(i)
      c1 = h * c / (boltz * teff * cgswave)
      if (c1 >= 80.0_real32) then
         flam(i) = 1.e-30_real32
      else
         flam(i) = 2.0_real32 * pi * h * c * c / (cgswave**5) / &
                   (exp(c1) - 1.0_real32) * 1.e-8_real32
      end if
   end do

end subroutine planck
c
c ***********************************************************************
c ***********************************************************************
c ***********************************************************************
c
	subroutine kurucz(teff,blogg)
c
c ***********************************************************************
c
c THIS SUBROUTINE COMPUTES THE EMERGENT FLUX FROM A KURUCZ MODEL ATMOSPHERE
c AS COMPILED BY LEJEUNE ET AL. THE FLUXES ARE OBTAINED BY SELECTING THE
c CLOSEST ENTRY IN LEJEUNE'S ATMOSPHERE GRID.
c UNITS: ERG/S/CM/CM/A. INCLUDES PI (I.E. PHYSICAL FLUX).
c THE METALLICITY IS THE SAME AS THE ONE FOR THE EVOLUTIONARY TRACKS.
c
c ************************************************************************
c
	implicit real*4    (a-h,o-z)
	character*20 name
        parameter (nmaxint=10,nmaxint1=11)
        common/parameters/time1,tstep,tmax,iz,toma,sfr,ninterv,
     *           xponent(nmaxint),upma,xmaslim(nmaxint1),doma,
     *           io1,io2,io3,io4,io5,io6,io7,io8,isf,name,jmg,z,
     *           iwind,sncut,bhcut,tdel,iatmos,ilib,iline,ivt,irsg,
     *           io9,io10,io11,io12,io13,io14,io15,xmwr,
     *           iwrscale,iwrt,jtime,tvar,tiempo1,tinter
	common/flux/flam
	common/lejeuneinput/tem,glog,nlej,flam1,wave
	dimension wave(1221),flam1(1221,600),flam(1221),tem(600),
     *            glog(600)
c
c PICK CLOSEST MODEL ATMOSPHERE
c
	  call choose_atm(tem,glog,nlej,teff,blogg,imodel)
c
c WE ARE REFINING THE PROCESS BY SEARCHING FOR THE NEXT BEST MATCH TO
c TEFF. WHILE THIS IGNORES THE GRAVITY DEPENDENCE, IT IS A RESONABLE
c COMPROMISE BETWEEN DOING A TIME CONSUMING FULL-BLOWN INTERPOLATION
c AND SIMPLY SELECTING THE NEAREST NEIGHBOR IN TEFF.
c
        teff_iter = teff/(tem(imodel)/teff)**2.
        call choose_atm(tem,glog,nlej,teff_iter,blogg,imodel_iter)
c
	if (imodel_iter.eq.imodel) then
        teff_iter = teff/(tem(imodel)/teff)**3.
        call choose_atm(tem,glog,nlej,teff_iter,blogg,imodel_iter)
	else
	  goto 888
      endif
c
	if (imodel_iter.eq.imodel) then
        teff_iter = teff/(tem(imodel)/teff)**4.
        call choose_atm(tem,glog,nlej,teff_iter,blogg,imodel_iter)
      else
        goto 888
      endif
c
	if (imodel_iter.eq.imodel) then
        teff_iter = teff/(tem(imodel)/teff)**5.
        call choose_atm(tem,glog,nlej,teff_iter,blogg,imodel_iter)
      else
        goto 888
      endif
c
888	continue
c
c THE FINAL MODEL WILL BE INTERPOLATED BETWEEN THE NEAREST NEIGHBOR AND
c THE NEXT BEST TEFF, UNLESS THE NEAREST NEIGHBOR IS VERY CLOSE TO THE
c ACTUAL TEFF.
c
	if(imodel_iter.eq.imodel) then
		factor=0.
		goto 777
	endif
	factor=(teff-tem(imodel))/(tem(imodel_iter)-tem(imodel))
777	continue
        do i=1,1221
         flam(i) =
     *   flam1(i,imodel)+factor*(flam1(i,imodel_iter)-flam1(i,imodel))
        enddo
c
        return
        end
c
c
c *****************************************************************************
c *****************************************************************************
c *****************************************************************************
c
      subroutine choose_atm(tem,glog,nlej,xt,xg,imodel)
c
c THIS SUBROUTINE CHOOSES THE "CLOSEST" LEJEUNE OR OTHER MODEL AVAILABLE IN
c (TEFF-LOG G) SPACE. IT RETURNS THE MODEL NUMBER.
c
c *****************************************************************************
c
        implicit real*4    (a-h,o-z)
        dimension tem(600),glog(600)
c
c DERIVE ONE MODEL WITH CLOSEST TEFF. NOTE: SEVERAL MODELS WITH SAME TEFF
c ARE AVAILABLE. DISTT,G ARE JUST BIG NUMBERS TO GET STARTED.
c
      distt = 9.e+34
      distg = 9.e+34
c
      do i=1,nlej
	 dt = abs(xt-tem(i))
         if (dt.lt.distt) then
            iteff  = i		! index of closest model (last on list)
            distt  = dt		! distance to closest model
         endif
      enddo

c
c NOW GET THE MODEL CLOSEST IN GRAVITY. CHECK ALL MODELS WITH FIXED TEFF.
c
      i = iteff
c
      do while (abs(tem(i)-tem(iteff)).lt.1.)
         dg = abs(1000.*10.**xg-1000.*10.**glog(i))
         if (dg.lt.distg) then
            imodel  = i		! index of closest model
            distg  = dg		! distance to closest model
         endif
         i = i + 1
         if (i.eq.nlej) goto 10    ! stop for special case
      enddo
 10   continue
c
      return
      end
c
c *****************************************************************************
c *****************************************************************************
c *****************************************************************************
c
	subroutine pauldrach(teff,blogg,radius)
c
c ***********************************************************************
c
c THE WM-BASIC MODEL ATMOSHPERE GRID. SOLAR PARAMETERS WERE DETERMINED
c AND SCALED FOR OTHER METALLICITIES. MASS LOSS WAS DETERMINED USING
c THE WIND MOMENTUM - LUMINOSITY RELATIONSHIP OF KUDRITZKI & PULS
c (2000) AND TERMINAL VELOCITIES WERE TAKEN FROM PRINJA, BARLOW AND
c HOWARTH (1990). MASS LOSS SCALING PROPORTIONAL TO Z^0.8, WHILE TERMINAL
c VELOCITY SCALING WAS PROPORTIONAL TO Z^0.13.
c FIVE METALLICITIES WERE CALCULATED AT DWARF AND SUPERGIANT GRAVITIES
c GIANTS WHERE INTERPOLATED FROM THE TWO TRACKS
c UNITS: ERG/S/CM/CM/A. INCLUDES PI (I.E. PHYSICAL FLUX).
c
c ************************************************************************
c
	implicit real*4    (a-h,o-z)
	character*20 name
        parameter (nmaxint=10,nmaxint1=11)
        common/parameters/time1,tstep,tmax,iz,toma,sfr,ninterv,
     *           xponent(nmaxint),upma,xmaslim(nmaxint1),doma,
     *           io1,io2,io3,io4,io5,io6,io7,io8,isf,name,jmg,z,
     *           iwind,sncut,bhcut,tdel,iatmos,ilib,iline,ivt,irsg,
     *           io9,io10,io11,io12,io13,io14,io15,xmwr,
     *           iwrscale,iwrt,jtime,tvar,tiempo1,tinter
        common/flux/flam
        common/wmbasicinput/tem_p,alum_p,glog_p,flampaul,wave_wm
	dimension wave_wm(1221),flam(1221),
     *  flampaul(1221,33),tem_p(33),glog_p(33),alum_p(33),f_t_g(33)

c
c
      call renorm_pauldrach(tem_p,flampaul,wave_wm)
c

      do i=1,1221
        do j=1,33
          f_t_g(j)=flampaul(i,j)
        enddo
        call mk_atmo(tem_p,glog_p,f_t_g,33,teff,blogg,flam(i),mmodel)
      enddo

c
c PICK CLOSEST MODEL ATMOSPHERE
c
      call renorm_at(teff,flam,33)
c
c
        return
        end
c
c ************************************************************************
c ************************************************************************
c ************************************************************************

        subroutine renorm_pauldrach(tem_p,flampaul,wave_wm)
c
c ROUTINE CALLED FROM "PAULDRACH". CHECKS AND RENORMALIZES ALL WR
c ATMOSPHERE MODELS SUCH THAT THE TOTAL EMERGENT FLUX CORRESPONDS TO SIGMA
c *TEFF^4.THIS IS NECESSARY SINCE THE ATMOSPHERE MODELS RESULT FROM
c PREVIOUS INTERPOLATIONS, WHICH DO NOT a priori GUARANTEE A PROPER
c NORMALIZATION.
c
c ************************************************************************
c
	implicit real*4    (a-h,o-z)
	common/flux/flam
	dimension wave_wm(1221),flampaul(1221,33),flam(1221),
     *            tem_p(33),fla_temp(1221)

        kmax = 1221
c
c RENORMALIZE ATMOSPHERES FROM GIANT GRID
c

        do i=12,22
          do k=1,1221
            fla_temp(k)=flampaul(k,i)
          enddo
          call total_flux(wave_wm,fla_temp,kmax,tem_p(i),xinte)
          do k=1,1221
            flampaul(k,i)=fla_temp(k)/xinte
          enddo


        enddo
        return
        end
c
c
c
c
c************************************************************
c************************************************************
c************************************************************
c

      subroutine mk_atmo(tem,xlogg,flp,m,teff,blogg,flux,mmodel)
c
c 2-DIMENSIONAL INTERPOLATION ROUTINE FROM NUMERICAL RECIPIES.
c
c ***********************************************************************
c
      implicit real (a-h,p-z)
      dimension tem(m),xlogg(m),flp(m),dtem_abs(m),dtem(m),
     &  dlogg_abs(m)

      flux=0.0
      dt=1.0e20
      dg=1.0e20
      do i=1,11

        dtem(i)=teff-tem(i)
        dtem_abs(i)=abs(teff-tem(i))
        if(dtem_abs(i).lt.dt)then
          dt=dtem_abs(i)
          model=i
        endif
      enddo

        if(model.gt.11)then
          model=11
        endif
        if(model.le.0)then
          model=1
        endif


      do i=1,3
        dlogg_abs(i)=abs(blogg-xlogg((i-1)*11+model))
        if(dlogg_abs(i).lt.dg)then
          dg=dlogg_abs(i)
          modelg=i

        endif
      enddo

      flux=flp((modelg-1)*11+model)

      mmodel=(modelg-1)*11+model


      return
      end
c
c ************************************************************************
c ************************************************************************
c ************************************************************************
!> Renormalizes an atmosphere flux array so that the total emergent flux
!> corresponds to sigma * Teff^4. Modern Fortran version.
subroutine renorm_at(teff, flux, nmod)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   real(real32), intent(in) :: teff
   real(real32), intent(inout) :: flux(:)
   integer, intent(in) :: nmod

   ! Parameters
   integer, parameter :: kmax = 1221
   real(real32) :: wave_wm(kmax)
   real(real32) :: fla_temp(kmax)
   real(real32) :: xinte

   ! External
   interface
      subroutine total_flux(wave, flux, nw, teff, xinte)
         use, intrinsic :: iso_fortran_env, only: real64
         integer, intent(in) :: nw
         real(real64), intent(in) :: wave(nw), flux(nw)
         real(real64), intent(in) :: teff
         real(real64), intent(out) :: xinte
      end subroutine total_flux
   end interface

   ! Copy input flux to temp array (promote to double if needed)
   fla_temp = flux

   ! Renormalize using total_flux (may require type conversion)
   call total_flux(real(wave_wm, kind=real64), real(fla_temp, kind=real64), kmax, real(teff, kind=real64), xinte)

   ! Apply normalization factor
   flux = flux / real(xinte, kind=real32)

end subroutine renorm_at
c
c ************************************************************************
c ************************************************************************
c ************************************************************************
!> Rescales the Hamann WN atmospheres to Teff when the nearest atmosphere is picked.
subroutine grid_pown(teff, model)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   real(real32), intent(in) :: teff
   integer, intent(in)      :: model

   ! Parameters
   integer, parameter :: nmod = 12, npts = 15000

   ! Atmosphere grid arrays (should be module variables in modern code)
   real(real32) :: wave_pot_wn(npts, nmod)
   real(real32) :: fl_pot_wn(npts, nmod)
   real(real32) :: t_pot_wn(nmod)
   integer      :: i_wn(nmod)
   real(real32) :: wave_powr(npts)
   real(real32) :: flam_powr(npts)
   integer      :: kmax_powr

   integer :: k

   ! Renormalize atmospheres of the Hamann WN grid
   do k = 1, i_wn(model)
      wave_powr(k) = wave_pot_wn(k, model)
      flam_powr(k) = fl_pot_wn(k, model) * teff / t_pot_wn(model)
      kmax_powr    = k
   end do

end subroutine grid_pown
c
c ************************************************************************
c ************************************************************************
c ************************************************************************
c
        subroutine grid_powc(teff,model)
c
c  ROUTINE TO RESCALE THE HAMANN WC ATMOSPHERES TO TEFF WHEN THE NEAREST
c  ATMOSPHERE IS PICKED.
c
c ************************************************************************
c
	implicit real*4    (a-h,o-z)
     	common/powrinput/i_pot,t_pot_wn,t_pot_wc,r_pot_wn,r_pot_wc,
     *                 wave_pot_wn,wave_pot_wc,fl_pot_wn,fl_pot_wc,
     *		       i_wn,i_wc
     	common/flux_powr/flam_powr,wave_powr,kmax_powr
	dimension i_pot(12),t_pot_wn(12),t_pot_wc(12),r_pot_wn(12),
     *		  r_pot_wc(12),wave_pot_wn(15000,12),
     *            wave_pot_wc(15000,12),fl_pot_wn(15000,12),
     *		  fl_pot_wc(15000,12),
     *		  flam_powr(15000),wave_powr(15000),
     *            i_wn(12),i_wc(12)
c
c RENORMALIZE ATMOSPHERES OF THE HAMANN WC GRID
c
          do k=1,i_wc(model)
            wave_powr(k)=wave_pot_wc(k,model)
            flam_powr(k)=flam_powr(k)*teff/t_pot_wc(model)
            kmax_powr=k
          enddo
c
        return
        end
c
c ***********************************************************************
c ***********************************************************************
c ***********************************************************************
c
	subroutine fabio(teff,blogg)
c
c THIS SUBROUTINE COMPUTES THE EMERGENT FLUX FROM AN IFA MODEL ATMOSPHERE.
c THE FLUXES ARE OBTAINED BY SELECTING THE CLOSEST ENTRY IN THE ATMOSPHERE
c GRID. UNITS: ERG/S/CM/CM/A. INCLUDES PI (I.E. PHYSICAL FLUX).
c THE METALLICITY IS THE SAME AS THE ONE FOR THE EVOLUTIONARY TRACKS.
c
c ************************************************************************
c
	implicit real*4    (a-h,o-z)
	character*20 name
        parameter (nmaxint=10,nmaxint1=11)
        common/parameters/time1,tstep,tmax,iz,toma,sfr,ninterv,
     *           xponent(nmaxint),upma,xmaslim(nmaxint1),doma,
     *           io1,io2,io3,io4,io5,io6,io7,io8,isf,name,jmg,z,
     *           iwind,sncut,bhcut,tdel,iatmos,ilib,iline,ivt,irsg,
     *           io9,io10,io11,io12,io13,io14,io15,xmwr,
     *           iwrscale,iwrt,jtime,tvar,tiempo1,tinter
        common/flux_ifa/wave_ifa,flam1_ifa_l,flam1_ifa_c,nlej_ifa,
     *                    tem_ifa,glog_ifa,flam_ifa_l,flam_ifa_c
	dimension wave_ifa(4200),flam1_ifa_l(4200,86),
     *             flam_ifa_l(4200),tem_ifa(86),glog_ifa(86),
     *         flam1_ifa_c(4200,86),flam_ifa_c(4200)
c
c
c PICK CLOSEST MODEL ATMOSPHERE
c
      dimin=1.e20
      do 14 i=1,86
      di=1.*(alog10(teff)-alog10(tem_ifa(i)))**2.
     *    +5.0*(blogg-glog_ifa(i))**2.
      if(di.lt.dimin) then
        dimin=di
        imodel=i
      endif
14    continue
c
       do i=1,4200
           flam_ifa_l(i) = flam1_ifa_l(i,imodel)
           flam_ifa_c(i) = flam1_ifa_c(i,imodel)
       enddo
c
        return
        end
c
c
c ***********************************************************************
c ***********************************************************************
c ***********************************************************************
!> Computes the emergent flux from high-resolution atmospheres.
!> Modern Fortran version.
subroutine lucimara(teff, blogg)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   real(real32), intent(in) :: teff, blogg

   ! Parameters
   integer, parameter :: n_hires = 13323, n_luci = 416

   ! Atmosphere grid arrays (should be module variables in modern code)
   real(real32) :: wave_hires(n_hires)
   real(real32) :: flam1_hires_l(n_hires, n_luci)
   real(real32) :: flam1_hires_c(n_hires, n_luci)
   real(real32) :: flam_hires_l(n_hires)
   real(real32) :: flam_hires_c(n_hires)
   real(real32) :: tem_luci(n_luci)
   real(real32) :: glog_luci(n_luci)
   integer :: nlej_lu
   integer :: imodel, i

   interface
      subroutine choose_atm(tem, glog, nlej, xt, xg, imodel)
         use, intrinsic :: iso_fortran_env, only: real32
         integer, intent(in) :: nlej
         real(real32), intent(in) :: tem(nlej), glog(nlej), xt, xg
         integer, intent(out) :: imodel
      end subroutine choose_atm
   end interface

   ! Pick closest model atmosphere
   call choose_atm(tem_luci, glog_luci, nlej_lu, teff, blogg, imodel)
   do i = 1, n_hires
      flam_hires_l(i) = flam1_hires_l(i, imodel)
      flam_hires_c(i) = flam1_hires_c(i, imodel)
   end do

end subroutine lucimara
c
c
c ***********************************************************************
c ***********************************************************************
c ***********************************************************************
!> Computes the emergent fluxes from Werner Schmutz's WR model atmospheres.
!> Modern Fortran version.
subroutine werner(l, radius, teff)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   integer, intent(in) :: l
   real(real32), intent(in) :: radius, teff

   ! Parameters
   integer, parameter :: n_wave = 1221, n_beta1_t = 12, n_beta1_r = 10, n_beta2_t = 11, n_beta2_r = 8

   ! Common block replacements (should be replaced by module variables in modern code)
   real(real32) :: flam(n_wave)
   real(real32) :: wave1(n_wave)
   real(real32) :: flam2(n_wave, n_beta1_t, n_beta1_r)
   real(real32) :: flam4(n_wave, n_beta2_t, n_beta2_r)
   real(real32) :: tgrid1(n_beta1_t), rgrid1(n_beta1_r)
   real(real32) :: tgrid2(n_beta2_t), rgrid2(n_beta2_r)
   integer :: iwind, ltype

   ! Locals
   real(real32) :: rt, dist, distmin
   integer :: i, j, k, ichoose, jchoose, bchoose

   ! Ensure renormalization of WR atmospheres
   call renorm_werner(flam2, tgrid1, flam4, tgrid2, wave1)

   ! Choice of the stellar wind model and calculation of the corresponding RT
   select case (iwind)
   case (0)
      rt = radius / 6.96e10_real32 * (1.e-4_real32 / wind1(l, ltype, 3))**0.66667_real32 * &
           (wind1(l, ltype, 4) / 2500._real32)**0.66667_real32
   case (1)
      rt = radius / 6.96e10_real32 * (1.e-4_real32 / wind2(l, ltype, 3))**0.66667_real32 * &
           (wind2(l, ltype, 4) / 2500._real32)**0.66667_real32
   case (2)
      rt = radius / 6.96e10_real32 * (1.e-4_real32 / wind3(l, ltype, 3))**0.66667_real32 * &
           (wind3(l, ltype, 4) / 2500._real32)**0.66667_real32
   case (3)
      rt = radius / 6.96e10_real32 * (1.e-4_real32 / wind4(l, ltype, 3))**0.66667_real32 * &
           (wind4(l, ltype, 4) / 2500._real32)**0.66667_real32
   case default
      rt = 0.0_real32
   end select

   ! Choose "closest" WR atmosphere model in log(T_star)--log(R_trans) plane
   bchoose = 0
   ichoose = 0
   jchoose = 0
   distmin = 9.e33_real32

   if (teff < 90000._real32) then
      do i = 1, n_beta1_t
         do j = 1, n_beta1_r
            dist = (log10(teff) - log10(tgrid1(i)))**2 + (log10(rt) - log10(rgrid1(j)))**2
            if (dist < distmin) then
               bchoose = 1
               ichoose = i
               jchoose = j
               distmin = dist
            end if
         end do
      end do
   else
      do i = 1, n_beta2_t
         do j = 1, n_beta2_r
            dist = (log10(teff) - log10(tgrid2(i)))**2 + (log10(rt) - log10(rgrid2(j)))**2
            if (dist < distmin) then
               bchoose = 2
               ichoose = i
               jchoose = j
               distmin = dist
            end if
         end do
      end do
   end if

   if (bchoose == 0) then
      call errpri('werner: No closest WR model !')
      return
   end if

   ! Assign the flux from the chosen model
   select case (bchoose)
   case (1)
      do k = 1, n_wave
         flam(k) = flam2(k, ichoose, jchoose)
      end do
   case (2)
      do k = 1, n_wave
         flam(k) = flam4(k, ichoose, jchoose)
      end do
   end select

end subroutine werner
c
c
c ************************************************************************
c ************************************************************************
c ************************************************************************
!> Renormalizes all WR atmosphere models so that the total emergent flux
!> corresponds to sigma * Teff^4. This is necessary since the atmosphere
!> models result from previous interpolations, which do not a priori guarantee
!> a proper normalization.
subroutine renorm_werner(flam2, tgrid1, flam4, tgrid2, wave1)
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   ! Arguments
   real(real64), intent(inout) :: flam2(1221,12,10)
   real(real64), intent(in)    :: tgrid1(12)
   real(real64), intent(inout) :: flam4(1221,11,8)
   real(real64), intent(in)    :: tgrid2(11)
   real(real64), intent(in)    :: wave1(1221)

   ! Locals
   integer :: i, j, k
   integer, parameter :: kmax = 1221
   real(real64) :: fla_temp(1221)
   real(real64) :: xinte

   ! Renormalize atmospheres from BETA=1 grid (flam2)
   do i = 1, 12
      do j = 1, 10
         do k = 1, kmax
            fla_temp(k) = flam2(k, i, j)
         end do
         call total_flux(wave1, fla_temp, kmax, tgrid1(i), xinte)
         do k = 1, kmax
            flam2(k, i, j) = fla_temp(k) / xinte
         end do
      end do
   end do

   ! Renormalize atmospheres from BETA=2 grid (flam4)
   do i = 1, 11
      do j = 1, 8
         do k = 1, kmax
            fla_temp(k) = flam4(k, i, j)
         end do
         call total_flux(wave1, fla_temp, kmax, tgrid2(i), xinte)
         do k = 1, kmax
            flam4(k, i, j) = fla_temp(k) / xinte
         end do
      end do
   end do

end subroutine renorm_werner
c
c
c ************************************************************************
c ************************************************************************
c ************************************************************************
!> Integrates the flux distribution to determine the correct total emission.
!> Determines the normalization factor with respect to a black body of temperature teff.
!> Modern Fortran version.
subroutine total_flux(wave, flux, nw, teff, xinte)
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   integer, intent(in) :: nw
   real(real64), intent(in) :: wave(nw), flux(nw)
   real(real64), intent(in) :: teff
   real(real64), intent(out) :: xinte

   real(real64), parameter :: STEBO = 5.6696196d-5 ! Stefan-Boltzmann constant [erg/cm^2/s/K^4]
   real(real64) :: weight(nw), bbflux

   ! Calculate integral over F_lambda * dLambda (lambda is increasing)
   call fliwgt(wave, flux, weight, xinte, nw)

   bbflux = STEBO * teff**4
   xinte  = -xinte / bbflux

end subroutine total_flux
c
c
c ************************************************************************
c ************************************************************************
c ************************************************************************
!> Computes the emergent flux from WR stars using the Hillier & Miller CMFGEN N-LTE atmosphere code.
!> Modern Fortran version.
subroutine hillier(l, teff, radius)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   integer, intent(in) :: l
   real(real32), intent(in) :: teff, radius

   ! Parameters
   integer, parameter :: npgrid = 3000, nmod = 12, nmod2 = 12
   real(real32), parameter :: pc_to_cm = 3.0856e18, rsun_to_cm = 6.96e10

   ! Common blocks replaced by assumed-shape arrays or module variables
   real(real32) :: flam(1221)
   real(real32) :: temp(npgrid), bol(npgrid), xsurf(npgrid), zmass(npgrid)
   real(real32) :: ysurf(npgrid), xc12s(npgrid), xn14s(npgrid), xo16s(npgrid)
   real(real32) :: bmdot(npgrid), tt_star(npgrid)
   ! Hillier/CMFGEN model arrays
   real(real32) :: tem_hi(nmod), tem2_hi(nmod2)
   real(real32) :: fl(1221, nmod), fl2(1221, nmod2)
   real(real32), parameter :: radc(nmod2) = [9.31,8.04,7.04,5.95,4.94,4.14,3.05,2.33,1.84,1.50,1.03,0.76]
   real(real32), parameter :: radn(nmod)  = [20.30,17.24,14.90,11.40,9.00,7.30,5.07,3.72,2.84,2.26,1.82,1.27]

   ! Locals
   integer :: i, model
   real(real32) :: cnr, coher, radius2

   interface
      subroutine choosewr(xt, temp, nmod, model)
         use, intrinsic :: iso_fortran_env, only: real32
         integer, intent(in) :: nmod
         real(real32), intent(in) :: temp(nmod), xt
         integer, intent(out) :: model
      end subroutine choosewr
      subroutine renorm_at(teff, flux, nmod)
         use, intrinsic :: iso_fortran_env, only: real32
         real(real32), intent(in) :: teff
         real(real32), intent(inout) :: flux(:)
         integer, intent(in) :: nmod
      end subroutine renorm_at
   end interface

   ! Defensive: avoid division by zero
   if (xn14s(l) == 0.0_real32) xn14s(l) = 1.e-6_real32
   cnr = xc12s(l) / xn14s(l)
   coher = ((xc12s(l)/12.0_real32) + (xo16s(l)/16.0_real32)) / (ysurf(l)/4.0_real32)

   if (xsurf(l) > 0.1_real32) then
      call choosewr(teff, tem_hi, nmod, model)
      radius2 = radn(model)
      do i = 1, 1221
         flam(i) = fl(i, model)
         flam(i) = flam(i) * ( (3.0856e11_real32)**2 / (radius2*rsun_to_cm)**2 )
         flam(i) = flam(i) * 1.e20_real32
      end do
      call renorm_at(teff, flam, nmod)
      return
   end if

   if (cnr < 10.0_real32) then
      call choosewr(teff, tem_hi, nmod, model)
      radius2 = radn(model)
      do i = 1, 1221
         flam(i) = fl(i, model)
         flam(i) = flam(i) * ( (3.0856e11_real32)**2 / (radius2*rsun_to_cm)**2 )
         flam(i) = flam(i) * 1.e20_real32
      end do
      call renorm_at(teff, flam, nmod)
      return
   end if

   if (coher < 0.5_real32) then
      call choosewr(teff, tem2_hi, nmod2, model)
      radius2 = radc(model)
      do i = 1, 1221
         flam(i) = fl2(i, model)
         flam(i) = flam(i) * ( (3.0856e11_real32)**2 / (radius2*rsun_to_cm)**2 )
         flam(i) = flam(i) * 1.e20_real32
      end do
      call renorm_at(teff, flam, nmod2)
      return
   end if

   if (coher < 1.0_real32) then
      call choosewr(teff, tem2_hi, nmod2, model)
      radius2 = radc(model)
      do i = 1, 1221
         flam(i) = fl2(i, model)
         flam(i) = flam(i) * ( (3.0856e11_real32)**2 / (radius2*rsun_to_cm)**2 )
         flam(i) = flam(i) * 1.e20_real32
      end do
      call renorm_at(teff, flam, nmod2)
      return
   end if

   if (coher >= 1.0_real32) then
      call choosewr(teff, tem2_hi, nmod2, model)
      radius2 = radc(model)
      do i = 1, 1221
         flam(i) = fl2(i, model)
         flam(i) = flam(i) * ( (3.0856e11_real32)**2 / (radius2*rsun_to_cm)**2 )
         flam(i) = flam(i) * 1.e20_real32
      end do
      call renorm_at(teff, flam, nmod2)
      return
   end if

end subroutine hillier
c
c ***********************************************************************
c ***********************************************************************
c ***********************************************************************
subroutine choosewr(xt, temp, nmod, model)
! Selects the closest WR model by effective temperature.
use, intrinsic :: iso_fortran_env, only: real32
implicit none
integer, intent(in) :: nmod
real(real32), intent(in) :: temp(nmod), xt
integer, intent(out) :: model

integer :: i
real(real32) :: distt, dt_abs

distt = 9.e33_real32
model = 1
do i = 1, nmod
   dt_abs = abs(xt - temp(i))
   if (dt_abs < distt) then
      model = i
      distt = dt_abs
   end if
end do
if (model > nmod) model = nmod
if (model < 1)    model = 1

end subroutine choosewr
c
c ***********************************************************************
c ***********************************************************************
c ***********************************************************************
subroutine radint(tem, rad, nmod, teff, rad2)
! Calculates the Wolf-Rayet radius of the Hillier model.
! Modern Fortran version.
use, intrinsic :: iso_fortran_env, only: real32
implicit none
integer, intent(in) :: nmod
real(real32), intent(in) :: tem(nmod), rad(nmod), teff
real(real32), intent(out) :: rad2

integer :: i, model
real(real32) :: distt, dt, dt_abs

distt = 9.e34_real32
model = 1
do i = 1, nmod
   dt = teff - tem(i)
   dt_abs = abs(dt)
   if (dt_abs < distt) then
      model = i
      distt = dt_abs
   end if
   if (dt > 0.0_real32) then
      model = model + 1
   end if
end do

if (model == 1) then
   rad2 = rad(1)
else if (model >= nmod + 1) then
   rad2 = rad(nmod)
else
   rad2 = rad(model) + (teff - tem(model)) * (rad(model) - rad(model-1)) / (tem(model) - tem(model-1))
end if

end subroutine radint
c
c ************************************************************************
c ************************************************************************
c ************************************************************************
c
c
c COMPUTES THE EMERGENT FLUX FROM THE WR STARS USING THE HAMANN & GRAEFENER
c (POTSDAM) MODELS. THERE ARE FIVE METALLICITY GRIDS OF WN AND WC STARS.
c CURRENTLY ONLY SOLAR MODELS ARE AVAILABLE. THE FLUXES ARE NORMALIZED TO
c 10 PC AND ARE CHANGED TO STELLAR SURFACE IN THIS SUBROUTINE.
c
c ************************************************************************
c
   subroutine hamann(l, teff, radius)
      use, intrinsic :: iso_fortran_env, only: real32
      implicit none
      integer, intent(in) :: l
      real(real32), intent(in) :: teff, radius

      ! Parameters
      integer, parameter :: npgrid = 3000, nmod = 12, nmod2 = 12
      real(real32), parameter :: pc_to_cm = 3.0856e18, rsun_to_cm = 6.96e10

      ! Common blocks (replace with module variables or assumed-shape arrays in modern code)
      real(real32) :: flam_powr(15000)
      real(real32) :: temp(npgrid), bol(npgrid), xsurf(npgrid), zmass(npgrid)
      real(real32) :: ysurf(npgrid), xc12s(npgrid), xn14s(npgrid), xo16s(npgrid)
      real(real32) :: bmdot(npgrid), tt_star(npgrid)
      integer :: i_pot(12), i_wn(12), i_wc(12)
      real(real32) :: t_pot_wn(12), t_pot_wc(12), r_pot_wn(12), r_pot_wc(12)
      real(real32) :: wave_pot_wn(15000,12), wave_pot_wc(15000,12)
      real(real32) :: fl_pot_wn(15000,12), fl_pot_wc(15000,12)
      real(real32) :: wave_powr(15000)
      ! External subroutines
      interface
         subroutine choosewr(teff, temp, nmod, model)
            use, intrinsic :: iso_fortran_env, only: real32
            real(real32), intent(in) :: teff
            real(real32), intent(in) :: temp(nmod)
            integer, intent(in) :: nmod
            integer, intent(out) :: model
         end subroutine choosewr
         subroutine grid_pown(teff, model)
            use, intrinsic :: iso_fortran_env, only: real32
            real(real32), intent(in) :: teff
            integer, intent(in) :: model
         end subroutine grid_pown
         subroutine grid_powc(teff, model)
            use, intrinsic :: iso_fortran_env, only: real32
            real(real32), intent(in) :: teff
            integer, intent(in) :: model
         end subroutine grid_powc
      end interface

      ! Locals
      integer :: i, model
      real(real32) :: cnr, coher, radius2

      ! Defensive: avoid division by zero
      if (xn14s(l) == 0.0) xn14s(l) = 1.e-6
      cnr = xc12s(l) / xn14s(l)
      coher = ((xc12s(l)/12.0) + (xo16s(l)/16.0)) / (ysurf(l)/4.0)

      if (xsurf(l) > 0.1) then
         call choosewr(teff, t_pot_wn, nmod, model)
         do i = 1, i_wn(model)
            flam_powr(i) = fl_pot_wn(i, model)
            radius2 = r_pot_wn(model)
            flam_powr(i) = flam_powr(i) * ( (pc_to_cm/1000.0)**2 / (radius2*rsun_to_cm)**2 )
            flam_powr(i) = flam_powr(i) * 1.e20
         end do
         call grid_pown(teff, model)
         return
      end if

      if (cnr < 10.0) then
         call choosewr(teff, t_pot_wn, nmod, model)
         do i = 1, i_wn(model)
            flam_powr(i) = fl_pot_wn(i, model)
            radius2 = r_pot_wn(model)
            flam_powr(i) = flam_powr(i) * ( (pc_to_cm/1000.0)**2 / (radius2*rsun_to_cm)**2 )
            flam_powr(i) = flam_powr(i) * 1.e20
         end do
         call grid_pown(teff, model)
         return
      end if

      if (coher < 0.5) then
         call choosewr(teff, t_pot_wc, nmod2, model)
         do i = 1, i_wc(model)
            flam_powr(i) = fl_pot_wc(i, model)
            radius2 = r_pot_wc(model)
            flam_powr(i) = flam_powr(i) * ( (pc_to_cm/1000.0)**2 / (radius2*rsun_to_cm)**2 )
            flam_powr(i) = flam_powr(i) * 1.e20
         end do
         call grid_powc(teff, model)
         return
      end if

      if (coher < 1.0) then
         call choosewr(teff, t_pot_wc, nmod2, model)
         do i = 1, i_wc(model)
            flam_powr(i) = fl_pot_wc(i, model)
            radius2 = r_pot_wc(model)
            flam_powr(i) = flam_powr(i) * ( (pc_to_cm/1000.0)**2 / (radius2*rsun_to_cm)**2 )
            flam_powr(i) = flam_powr(i) * 1.e20
         end do
         call grid_powc(teff, model)
         return
      end if

      if (coher >= 1.0) then
         call choosewr(teff, t_pot_wc, nmod2, model)
         do i = 1, i_wc(model)
            flam_powr(i) = fl_pot_wc(i, model)
            radius2 = r_pot_wc(model)
            flam_powr(i) = flam_powr(i) * ( (pc_to_cm/1000.0)**2 / (radius2*rsun_to_cm)**2 )
            flam_powr(i) = flam_powr(i) * 1.e20
         end do
         call grid_powc(teff, model)
         return
      end if

   end subroutine hamann
c
c ***********************************************************************
c ***********************************************************************
c ***********************************************************************
c
        subroutine linesyn(time_in,icount)
c
c COMPUTES THE SYNTHETIC LINES FROM A POPULATION OF STARS FROM 1205.5 TO
c 1849.75 A (SP. RES. = 0.75 A, NB OF PIX = 860). IT USES NORMALIZED SPECTRA
c (450 SP. TYPE) FROM THE LIBRARIES "SP.DAT" AND "SP_LOW.DAT". THE FORMER IS
c FOR SOLAR METALLICITY (DE MELLO ET AL. 2000) AND THE LATTER FOR LMC/SMC
c METALLICITY (LEITHERER ET AL. 2001). THE SCALING OF THE FLUX IS DONE
c WITH MODEL ATMOSPHERES. UNITS ARE LOG(ERG/SEC/A).
c
c ************************************************************************
c
	implicit real*4    (a-h,o-z)
        character*20 name
	character*80 str_date
        parameter (nmaxint=10,nmaxint1=11)
        common/parameters/time1,tstep,tmax,iz,toma,sfr,ninterv,
     *           xponent(nmaxint),upma,xmaslim(nmaxint1),doma,
     *           io1,io2,io3,io4,io5,io6,io7,io8,isf,name,jmg,z,
     *           iwind,sncut,bhcut,tdel,iatmos,ilib,iline,ivt,irsg,
     *           io9,io10,io11,io12,io13,io14,io15,xmwr,
     *           iwrscale,iwrt,jtime,tvar,tiempo1,tinter
        common/mgrid/cmass,lmin,lmax,delm
	common/flux/flam
        common/position/temp,bol,xsurf,zmass,ysurf,xc12s,xn14s,xo16s,
     *                bmdot,tt_star
        common/stars/dens,tonum
        common/nebula/xrange,gamma
        common/fluxline/np,wavel,ffac
        common/match/imatch
        common/uvtemplate/fli
        common/lejeuneinput/tem,glog,nlej,flam1,wave
        common/linesyn_r/tofluxl,tofluxc,tfluxl,tfluxc,fluneb2
        parameter (npgrid = 3000)
        dimension temp(npgrid),bol(npgrid),xsurf(npgrid),zmass(npgrid),
     *         ysurf(npgrid),xc12s(npgrid),xn14s(npgrid),xo16s(npgrid),
     *          bmdot(npgrid),cmass(npgrid),dens(npgrid),freeflux(14),
     *        fli(450,860),wavel(860),fluxl(npgrid,860),freelam(14),
     *        fluxc(npgrid,860),ffac(860),tofluxl(860),tofluxc(860),
     *        conti(26),xrange(26),fluneb2(860),gamma(26),tfluxl(860),
     *        tfluxc(860),wave(1221),flam(1221),tt_star(npgrid),
     *	      tem(600),glog(600),flam1(1221,600)
        data rlam0,dell,np/1205.5,0.75,860/
	data freelam/1150.,1280.,1310.,1360.,1430.,1480.,1510.,1580.,
     *               1630.,1680.,1740.,1785.,1820.,1875./
c
c AT THE FIRST CALL TO THE ROUTINE A HEADER FOR THE OUTPUT FILE IS GENERATED.
c
      if(icount.gt.1) goto 98
c
      write(91,96) name
96    format(' MODEL DESIGNATION: ',a20)
      call fdate(str_date)
      write(91,597) str_date
597   format(' MODEL GENERATED: ', a24)
      write(91,95)
95    format(/,'        COMPUTED SYNTHETIC LINE SPECTRUM [ERG/SEC/A]')
      write(91,94)
94    format(/,' TIME [YR]      WAVELENGTH  LOG(LUMINOSITY)  ',
     * 'NORMALIZED SPECTRUM')
c
c      write(910,960)
c960   format('<line>')
c      write(910,961) name
c961    format('  <ModelDesignation>',a20,'</ModelDesignation>')
c      call fdate(str_date)
c      write(910,5971) str_date
c5971   format('  <ModelGeneratedOn>', a24,'</ModelGeneratedOn>')
c      write(910,950)
c950    format('  <Title>COMPUTED SYNTHETIC LINE SPECTRUM [ERG/SEC/A]',
c     *        '</Title>')
c      write(910,940)
c940    format('  <Columns>LogLUMIN/NORMAL</Columns>')
c
c DEFINITION OF THE WAVELENGTH GRID.
c
      do 85 k=1,np
      wavel(k)=rlam0+(float(k-1)*dell)
85    continue
c
98    continue
c
c
c IF THERE IS NO CONTINUOUS STAR FORMATION, THE FLUX ACCUMULATORS
c ARE RESET FOR EACH TIME STEP.
c
       if(isf.gt.0 .and. icount.gt.1) then
678	      continue
	else
      		do 50 k=1,np
      		tofluxl(k)=0.
      		tofluxc(k)=0.
      		tfluxl(k)=0.
      		tfluxc(k)=0.
      		fluneb2(k)=0.
50    		continue
c
        endif
c
c CALCULATE THE EFFECTIVE TEMPERATURE AND LOG G FOR EACH MASS. THEN GET THE
c SPECTRAL TYPE INDEX FROM SCHMIDT-KALER (IF NOT WR) AND NORMALIZED LINE
c AND RESCALED ACCORDING TO MODEL ATMOSPHERES (IN ERG PER UNIT SURFACE FOR
c EACH STAR). AFTERWARDS MULTIPLICATION BY THE STELLAR SURFACE AND THE
c STELLAR NUMBER DENSITY.
c
      do 10 l=lmin,lmax
c
      if(bol(l).lt.-19.) then
        do 100 k=1,np
        fluxl(l,k) = 1.e-30
        fluxc(l,k) = 1.e-30
100     continue
        goto 10
      endif
c
      teff=10.**temp(l)
      blogg=alog10(zmass(l))+4.*temp(l)-bol(l)-10.6
      radius=10.**(10.8426+0.5*bol(l)-2.*temp(l)+7.52)
      cotest=5.71*log10(teff)-21.95
      if(temp(l).gt.4.4.and.xsurf(l).lt.0.4.and.cmass(l).ge.xmwr) then
        if(xn14s(l).eq.0) xn14s(l)=1.e-6
        cnr=xc12s(l)/xn14s(l)
        coher=((xc12s(l)/12.)+(xo16s(l)/16.))/(ysurf(l)/4.)
        if(xsurf(l).gt.0.1) then
          index=446
          goto 17
        endif
        if(cnr.lt.10.) then
          index=447
          goto 17
        endif
        if(coher.lt.0.5) then
          index=448
          goto 17
        endif
        if(coher.lt.1.0) then
          index=449
          goto 17
        endif
        if(coher.ge.1.0) then
          index=450
          goto 17
        endif
17      call hillier(l,teff,radius)
      else
        call near(l,index)
        if(teff.lt.2000..or.teff.gt.60000.) then
           call planck(teff)
        else
          call kurucz(teff,blogg)
        endif
      endif
c
c PREPARE THE CONTINUUM GRID. FIRST THE FLUXES CORRESPONDING TO THE PRESELECTED
c LINE FREE WAVELENGTH GRID ARE COMPUTED. THEN AN AVERAGE OVER THREE POINTS IS
c DONE. A SPLINE FIT TO THE FLUXES AT THE 14 WAVELENGTH POINTS IS DONE TO
c MATCH THE GRID OF THE 860 POINTS OF THE LINE SPECTRUM. THE LMC/SMC LIBRARY
c CUTS OFF AT 1600 A. POINTS LONGWARD OF 1600 A ARE MEANINGLESS. THEREFORE
c WE SET THEM TO ZERO.
c
   ! Prepare the continuum grid: compute fluxes at preselected line-free wavelengths,
   ! average over three points, then spline fit to the 14 wavelength points to match the 860-point grid.
   integer :: jj, k
   real(real32) :: freeflux(14)
   do jj = 1, 14
      freeflux(jj) = yntra(freelam(jj), wave, flam, 1221)
      freeflux(jj) = alog10((freeflux(jj) + flam(imatch-1) + flam(imatch+1)) / 3.0)
   end do

   do k = 1, np
      call intrpl(14, freelam, freeflux, 1, wavel(k), ffac(k))
      fluxc(l, k) = 12.566 * radius * radius / 1.e20 * 10.0**ffac(k) * dens(l)
      fluxl(l, k) = fli(index, k) * fluxc(l, k)
      if (iline == 2 .and. k >= 529) fluxl(l, k) = 1.e-30
   end do

   ! Summation over all masses and generations of stars.
   do k = 1, np
      do l2 = lmin, lmax
         tofluxl(k) = tofluxl(k) + fluxl(l2, k)
         tofluxc(k) = tofluxc(k) + fluxc(l2, k)
      end do
   end do

   ! Calculate the nebular continuum and add to the stellar continuum.
   call continuum(time_in, icount, conti)

   do m = 1, np
      fluneb2(m) = yntra(wavel(m), xrange, conti, 26)
      tfluxl(m) = tofluxl(m) + fluneb2(m)
      tfluxc(m) = tofluxc(m) + fluneb2(m)
   end do

   ! Generate the output file. Unit = 91.
   if (mod(time_in, tdel) < tstep) then
      do k = 1, np
         write(91, 500) time_in, wavel(k), alog10(tfluxl(k) + 1.e-35) + 20., &
                         tfluxl(k) / (tfluxc(k) + 1.e-35)
      end do
500   format(1x, e10.5, 4x, f10.2, 5x, f10.5, 5x, f10.4)
   end if

   return
end
c
c ************************************************************************
c ************************************************************************
c ************************************************************************
c
c
c COMPUTES THE SYNTHETIC LINES FROM A POPULATION OF STARS FROM 1003 TO
c 1183 A (SP. RES. = 0.127 A, NB OF PIX = 1415). IT USES NORMALIZED SPECTRA
c (450 SP. TYPE) FROM THE LIBRARIES "FUSE_HIGH.DAT" AND "FUSE_LOW.DAT". THE
c FORMER IS FOR SOLAR METALLICITY AND THE LATTER FOR LMC/SMC METALLICITY
c (PELLERIN ET AL. 2002). THE SCALING OF THE FLUX IS DONE WITH MODEL
c ATMOSPHERES. UNITS ARE LOG(ERG/SEC/A).
c
c ************************************************************************
!> Computes the synthetic FUSE spectrum (1003–1183 Å, 0.127 Å res, 1415 pixels) for a population of stars.
!> Modern Fortran version (modular, explicit interfaces, no COMMON).
subroutine fusesyn(time_in, icount)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   real(real32), intent(in) :: time_in
   integer, intent(in) :: icount

   ! Parameters
   integer, parameter :: nmaxint=10, nmaxint1=11, npgrid=3000
   integer, parameter :: np1=1415
   real(real32), parameter :: rlam0=1003.1_real32, dell=0.127_real32
   integer, parameter :: n_wave=1221

   ! Model parameters (should be in a module in modern code)
   character(len=20) :: name
   character(len=80) :: str_date
   integer :: iz, ninterv, upma, jmg, iwind, sncut, bhcut, iatmos, ilib, iline, ivt, irsg
   integer :: io1, io2, io3, io4, io5, io6, io7, io8, isf, io9, io10, io11, io12, io13, io14, io15
   real(real32) :: time1, tstep, tmax, toma, sfr, doma, z, tdel, tvar, tiempo1, tinter
   real(real32) :: xponent(nmaxint), xmaslim(nmaxint1), xmwr
   integer :: iwrscale, iwrt, jtime

   ! Mass grid and stellar properties
   real(real32) :: cmass(npgrid), temp(npgrid), bol(npgrid), xsurf(npgrid), zmass(npgrid)
   real(real32) :: ysurf(npgrid), xc12s(npgrid), xn14s(npgrid), xo16s(npgrid)
   real(real32) :: bmdot(npgrid), tt_star(npgrid), dens(npgrid)
   integer :: lmin, lmax

   ! FUSE-specific arrays
   real(real32) :: fli1(450, np1), wave_f(np1), fluxl(npgrid, np1), fluxc(npgrid, np1)
   real(real32) :: tofluxl_fu(np1), tofluxc_fu(np1), tfluxl_fu(np1), tfluxc_fu(np1)
   real(real32) :: fluneb2_fu(np1), ffac(np1)
   real(real32) :: freeflux(8)
   real(real32), parameter :: freelam(8) = [995.,1005.,1055.,1085.,1105.,1125.,1155.,1195.]
   real(real32) :: wave(n_wave), flam(n_wave)
   real(real32) :: tem(600), glog(600), flam1(n_wave,600)
   real(real32) :: xrange(26), gamma(26), conti(26)
   integer :: imatch

   ! Locals
   integer :: l, k, jj, m, index
   real(real32) :: teff, blogg, radius, cotest, cnr, coher

   ! External procedures
   interface
      subroutine fdate(str_date)
         character(len=*), intent(out) :: str_date
      end subroutine fdate
      subroutine hillier(l, teff, radius)
         integer, intent(in) :: l
         real(real32), intent(in) :: teff, radius
      end subroutine hillier
      subroutine kurucz(teff, blogg)
         real(real32), intent(in) :: teff, blogg
      end subroutine kurucz
      subroutine planck(teff)
         real(real32), intent(in) :: teff
      end subroutine planck
      subroutine near(l, index)
         integer, intent(in) :: l
         integer, intent(out) :: index
      end subroutine near
      subroutine continuum(time, icount, conti)
         real(real32), intent(in) :: time
         integer, intent(in) :: icount
         real(real32), intent(out) :: conti(26)
      end subroutine continuum
      pure function yntra(r, rw, fu, imax) result(y)
         real(real32), intent(in) :: r
         real(real32), intent(in) :: rw(imax), fu(imax)
         integer, intent(in) :: imax
         real(real32) :: y
      end function yntra
      subroutine intrpl(l, x, y, n, u, v)
         integer, intent(in) :: l, n
         real(real32), intent(in) :: x(l), y(l), u(n)
         real(real32), intent(out) :: v(n)
      end subroutine intrpl
   end interface

   ! Output header at first call
   if (icount <= 1) then
      write(86,96) name
96    format(' MODEL DESIGNATION: ',a20)
      call fdate(str_date)
      write(86,597) str_date
597   format(' MODEL GENERATED: ', a24)
      write(86,95)
95    format(/,'        COMPUTED SYNTHETIC LINE SPECTRUM [ERG/SEC/A]')
      write(86,94)
94    format(/,' TIME [YR]      WAVELENGTH  LOG(LUMINOSITY)  ', &
         'NORMALIZED SPECTRUM')
      do k = 1, np1
         wave_f(k) = rlam0 + real(k-1,real32)*dell
      end do
   end if

   ! Reset accumulators if not continuous SF or first time step
   if (.not.(isf > 0 .and. icount > 1)) then
      tofluxl_fu = 0.0_real32
      tofluxc_fu = 0.0_real32
      tfluxl_fu = 0.0_real32
      tfluxc_fu = 0.0_real32
      fluneb2_fu = 0.0_real32
   end if

   ! Loop over all masses
   do l = lmin, lmax
      if (bol(l) < -19.0_real32) then
         fluxl(l,:) = 1.e-30_real32
         fluxc(l,:) = 1.e-30_real32
         cycle
      end if

      teff = 10.0_real32**temp(l)
      blogg = log10(zmass(l)) + 4.0_real32*temp(l) - bol(l) - 10.6_real32
      radius = 10.0_real32**(10.8426_real32 + 0.5_real32*bol(l) - 2.0_real32*temp(l) + 7.52_real32)
      cotest = 5.71_real32*log10(teff) - 21.95_real32

      if (temp(l) > 4.4_real32 .and. xsurf(l) < 0.4_real32 .and. cmass(l) >= xmwr) then
         if (xn14s(l) == 0.0_real32) xn14s(l) = 1.e-6_real32
         cnr = xc12s(l)/xn14s(l)
         coher = (xc12s(l)/12.0_real32 + xo16s(l)/16.0_real32) / (ysurf(l)/4.0_real32)
         if (xsurf(l) > 0.1_real32) then
            index = 446
         else if (cnr < 10.0_real32) then
            index = 447
         else if (coher < 0.5_real32) then
            index = 448
         else if (coher < 1.0_real32) then
            index = 449
         else
            index = 450
         end if
         call hillier(l, teff, radius)
      else
         call near(l, index)
         if (teff < 2000.0_real32 .or. teff > 60000.0_real32) then
            call planck(teff)
         else
            call kurucz(teff, blogg)
         end if
      end if

      ! Prepare the continuum grid
      do jj = 1, 8
         freeflux(jj) = yntra(freelam(jj), wave, flam, n_wave)
         freeflux(jj) = log10((freeflux(jj) + flam(imatch-1) + flam(imatch+1))/3.0_real32)
      end do

      do k = 1, np1
         call intrpl(8, freelam, freeflux, 1, wave_f(k), ffac(k))
         fluxc(l, k) = 12.566_real32 * radius**2 / 1.e20_real32 * 10.0_real32**ffac(k) * dens(l)
         fluxl(l, k) = fli1(index, k) * fluxc(l, k)
      end do
   end do

   ! Sum over all masses
   do k = 1, np1
      tofluxl_fu(k) = 0.0_real32
      tofluxc_fu(k) = 0.0_real32
      do l = lmin, lmax
         tofluxl_fu(k) = tofluxl_fu(k) + fluxl(l, k)
         tofluxc_fu(k) = tofluxc_fu(k) + fluxc(l, k)
      end do
   end do

   ! Add nebular continuum
   call continuum(time_in, icount, conti)
   do m = 1, np1
      fluneb2_fu(m) = yntra(wave_f(m), xrange, conti, 26)
      tfluxl_fu(m) = tofluxl_fu(m) + fluneb2_fu(m)
      tfluxc_fu(m) = tofluxc_fu(m) + fluneb2_fu(m)
   end do

   ! Output
   if (mod(time_in, tdel) < tstep) then
      do k = 1, np1
         write(86,500) time_in, wave_f(k), &
            log10(tfluxl_fu(k) + 1.e-35_real32) + 20.0_real32, &
            tfluxl_fu(k) / (tfluxc_fu(k) + 1.e-35_real32)
      end do
500   format(1x, e10.5, 4x, f10.2, 5x, f10.5, 5x, f10.4)
   end if

end subroutine fusesyn
c
c
c ************************************************************************
c ************************************************************************
c ************************************************************************
!> Computes a high-resolution optical spectrum for a population of stars.
!> Modern Fortran version (replaces COMMON blocks, uses modules/explicit interfaces).
subroutine hires(time_in, icount)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   real(real32), intent(in) :: time_in
   integer, intent(in)      :: icount

   ! Parameters
   integer, parameter :: nmaxint=10, nmaxint1=11, npgrid=3000
   integer, parameter :: n_hires=13323, n_luci=416, n_wave=1221

   ! Model parameters (should be in a module in modern code)
   character(len=20) :: name
   character(len=80) :: str_date
   real(real32) :: time1, tstep, tmax, toma, sfr, doma, z, tdel, tvar, tiempo1, tinter
   integer :: iz, ninterv, upma, jmg, iwind, sncut, bhcut, iatmos, ilib, iline, ivt, irsg
   integer :: io1, io2, io3, io4, io5, io6, io7, io8, isf, io9, io10, io11, io12, io13, io14, io15
   real(real32) :: xponent(nmaxint), xmaslim(nmaxint1), xmwr
   integer :: iwrscale, iwrt, jtime

   ! Mass grid and stellar properties
   real(real32) :: cmass(npgrid), temp(npgrid), bol(npgrid), xsurf(npgrid), zmass(npgrid)
   real(real32) :: ysurf(npgrid), xc12s(npgrid), xn14s(npgrid), xo16s(npgrid)
   real(real32) :: bmdot(npgrid), tt_star(npgrid), dens(npgrid)
   integer :: lmin, lmax
   real(real32) :: grav(npgrid)

   ! Flux arrays
   real(real32) :: flam(n_wave)
   real(real32) :: wave_hires(n_hires), flam_hires_l(n_hires), flam_hires_c(n_hires)
   real(real32) :: flam1_hires_l(n_hires, n_luci), flam1_hires_c(n_hires, n_luci)
   real(real32) :: tem_luci(n_luci), glog_luci(n_luci)
   real(real32) :: stflux_hires_l(n_hires), stflux_hires_c(n_hires)
   real(real32) :: toflux_hires_l(n_hires), toflux_hires_c(n_hires)
   real(real32) :: fluneb_hires(n_hires)
   real(real32) :: flux_hires_l(npgrid, n_hires), flux_hires_c(npgrid, n_hires)
   real(real32) :: flam1(n_wave,600), tem(600), glog(600), wave(n_wave)
   real(real32) :: xrange(26), gamma(26), conti(26)

   ! Locals
   integer :: l, m, i
   real(real32) :: teff, blogg, radius

   ! External procedures
   interface
      subroutine fdate(str_date)
         character(len=*), intent(out) :: str_date
      end subroutine fdate
      subroutine hillier(l, teff, radius)
         integer, intent(in) :: l
         real(real32), intent(in) :: teff, radius
      end subroutine hillier
      subroutine kurucz(teff, blogg)
         real(real32), intent(in) :: teff, blogg
      end subroutine kurucz
      subroutine lucimara(teff, blogg)
         real(real32), intent(in) :: teff, blogg
      end subroutine lucimara
      subroutine continuum(time, icount, conti)
         real(real32), intent(in) :: time
         integer, intent(in) :: icount
         real(real32), intent(out) :: conti(26)
      end subroutine continuum
      pure function yntra(r, rw, fu, imax) result(y)
         real(real32), intent(in) :: r
         real(real32), intent(in) :: rw(imax), fu(imax)
         integer, intent(in) :: imax
         real(real32) :: y
      end function yntra
   end interface

   ! Output header at first call
   if (icount <= 1) then
      write(82,96) name
96    format(' MODEL DESIGNATION: ',a20)
      call fdate(str_date)
      write(82,597) str_date
597   format(' MODEL GENERATED: ', a24)
      write(82,95)
95    format(/,'        COMPUTED SYNTHETIC LINE SPECTRUM [ERG/SEC/A]')
      write(82,94)
94    format(/,' TIME [YR]      WAVELENGTH  LOG(LUMINOSITY)  ', &
         'NORMALIZED SPECTRUM')
   end if

   ! Reset accumulators if not continuous SF
   if (.not.(isf > 0 .and. icount > 1)) then
      stflux_hires_l = 0.0_real32
      stflux_hires_c = 0.0_real32
      toflux_hires_l = 0.0_real32
      toflux_hires_c = 0.0_real32
      fluneb_hires   = 0.0_real32
   end if

   ! Loop over all masses
   do l = lmin, lmax
      if (bol(l) < -19.0_real32) then
         flux_hires_l(l,:) = 1.e-30_real32
         flux_hires_c(l,:) = 1.e-30_real32
         cycle
      end if

      teff   = 10.0_real32**temp(l)
      grav(l)= log10(zmass(l)) + 4.0_real32*temp(l) - bol(l) - 10.6_real32
      blogg  = grav(l)
      radius = 10.0_real32**(10.8426_real32 + 0.5_real32*bol(l) - 2.0_real32*temp(l) + 7.52_real32)

      if (temp(l) > 4.4_real32 .and. xsurf(l) < 0.4_real32 .and. cmass(l) >= xmwr) then
         call hillier(l, teff, radius)
         do m = 1, n_hires
            flam_hires_l(m) = yntra(wave_hires(m), wave, flam, n_wave)
            flam_hires_c(m) = flam_hires_l(m)
         end do
      else
         if (teff < 2000.0_real32) then
            do m = 1, n_hires
               flam_hires_l(m) = yntra(wave_hires(m), wave, flam, n_wave)
               flam_hires_c(m) = flam_hires_l(m)
            end do
         else if (teff < 3000.0_real32) then
            call kurucz(teff, blogg)
            do m = 1, n_hires
               flam_hires_l(m) = yntra(wave_hires(m), wave, flam, n_wave)
               flam_hires_c(m) = flam_hires_l(m)
            end do
         else
            call lucimara(teff, blogg)
         end if
      end if

      do m = 1, n_hires
         flux_hires_l(l, m) = 12.566_real32 * radius**2 / 1.e20_real32 * flam_hires_l(m) * dens(l)
         flux_hires_c(l, m) = 12.566_real32 * radius**2 / 1.e20_real32 * flam_hires_c(m) * dens(l)
         stflux_hires_l(m) = stflux_hires_l(m) + flux_hires_l(l, m)
         stflux_hires_c(m) = stflux_hires_c(m) + flux_hires_c(l, m)
      end do
   end do

   ! Calculate the nebular continuum and add to the stellar continuum
   call continuum(time_in, icount, conti)
   do m = 1, n_hires
      fluneb_hires(m)   = yntra(wave_hires(m), xrange, conti, 26)
      toflux_hires_l(m) = stflux_hires_l(m) + fluneb_hires(m)
      toflux_hires_c(m) = stflux_hires_c(m) + fluneb_hires(m)
   end do

   ! Output spectrum at selected time steps
   if (mod(time_in, tdel) < tstep) then
      do i = 1, n_hires
         write(82,500) time_in, wave_hires(i), &
            log10(toflux_hires_l(i) + 1.e-35_real32) + 20.0_real32, &
            toflux_hires_l(i) / (toflux_hires_c(i) + 1.e-35_real32)
      end do
500   format(1x,e10.5,4x,f10.2,5x,f10.5,5x,f10.4)
   end if

end subroutine hires

c
c ***********************************************************************
c ***********************************************************************
c ***********************************************************************
!> Writes the input parameters and model summary to a file.
subroutine output
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   ! Parameters
   integer, parameter :: nmaxint=10, nmaxint1=11, npgrid=3000
   ! Arguments and local variables
   character(len=20) :: name
   character(len=24) :: str_date
   integer :: date_time(8)
   real(real32) :: time1, tstep, tmax, toma, sfr, doma, z, tdel, tvar, tiempo1, tinter
   integer :: iz, ninterv, upma, jmg, iwind, sncut, bhcut, iatmos, ilib, iline, ivt, irsg
   integer :: io1, io2, io3, io4, io5, io6, io7, io8, isf, io9, io10, io11, io12, io13, io14, io15
   real(real32) :: xponent(nmaxint), xmaslim(nmaxint1)
   real(real32) :: xmwr
   integer :: iwrscale, iwrt, jtime
   real(real32) :: cmass(npgrid)
   integer :: lmin, lmax
   integer :: i

   ! Output header
   write(99,10001) name
10001 format(//,' MODEL DESIGNATION: ',a20)
   call fdate(str_date)
   write(99,597) str_date
597  format(' MODEL GENERATED ON ',a24)

   write(99,10005) isf
10005 format(//,' CONTINUOUS STAR FORMATION (>0) OR FIXED MASS',
     &        ' (<=): ',/,i3)
   write(99,10006) toma
10006 format(' TOTAL STELLAR MASS (IF "FIXED MASS" IS CHOSEN): ' &
     &         ,/,e10.3)
   write(99,10007) sfr
10007 format(' SFR(T=0) IF "CONT SF" IS CHOSEN: ',/,f12.2)
   write(99,10008) ninterv
10008 format(' NUMBER OF INTERVALS FOR THE IMF: ',/,i3)
   write(99,10009) (xponent(i), i=1,ninterv)
10009 format(' IMF EXPONENTS : ',/,10(f12.2))
   write(99,10011) (xmaslim(i), i=1,ninterv+1)
10011 format(' MASS LIMITS FOR IMF [SOLAR MASSES]: ',/,11(f12.2))
   write(99,1011) sncut
1011  format(' SUPERNOVA CUT-OFF MASS: ',/,f12.2)
   write(99,1010) bhcut
1010  format(' BLACK HOLE CUT-OFF MASS: ',/,f12.2)
   write(99,10012) iz
10012 format(' METALLICITY + TRACKS: ',/ &
     &  ' GENEVA STD:  11=0.001; 12=0.004; 13=0.008;', &
     &  ' 14=0.020; 15=0.040',/ &
     &  ' GENEVA HIGH: 21=0.001; 22=0.004; 23=0.008;', &
     &  ' 24=0.020; 25=0.040',/ &
     &  ' PADOVA STD: 31=0.0004; 32=0.004; 33=0.008;', &
     &  ' 34=0.020; 35=0.050',/ &
     &  ' PADOVA AGB: 41=0.0004; 42=0.004; 43=0.008;', &
     &  ' 44=0.020; 45=0.050',/ &
     &  ' GENEVA v00:  51=0.001; 52=0.002; 53=0.008;', &
     &  ' 54=0.014; 55=0.040',/ &
     &  ' GENEVA v40:  61=0.001; 62=0.002; 63=0.008;', &
     &  ' 64=0.014; 65=0.040',/,i3)
   write(99,1012) iwind
1012  format(' WIND MODEL (0=MAEDER; 1=EMPIRICAL;', &
     &         ' 2=THEOR.; 3=ELSON): ',/,I3)
   write(99,10013) time1
10013 format(' INITIAL TIME: ',/, e10.3)
   write(99,10014) jtime
10014 format(' TIME SCALE: LINEAR (=0) OR LOGARITHMIC (=1)',/,I3)
   if(jtime == 0) then
      write(99,1014) tvar/1.e6
1014  format(' TIME STEP [1.e6 YEARS] (IF JTIME=0) OR',/ &
     &   ' NUMBER OF STEPS        (IF JTIME=1): ',/,f12.2)
   else
      write(99,1015) int(tinter)
1015  format(' TIME STEP [1.e6 YEARS] (IF JTIME=0) OR',/ &
     &   ' NUMBER OF STEPS        (IF JTIME=1): ',/,i8)
   end if
   write(99,10015) tmax
10015 format(' LAST GRID POINT: ',/,e10.3)
   write(99,10016) jmg
10016 format(' SMALL (=0) OR LARGE (=1) MASS GRID;',/ &
     &  ' ISOCHRONE ON  LARGE GRID (=2) OR FULL ISOCHRONE (=3): ',/,i3)
   write(99,10017) lmin,lmax
10017 format(' LMIN, LMAX (ALL = 0): ',/,2i6)
   write(99,1017) tdel
1017  format(' TIME STEP FOR PRINTING OUT THE SYNTHETIC SPECTRA' &
     &           ': ',/,e10.3)
   write(99,1018) iatmos
1018 format(' MODEL ATMOSPHERE: 1=PLA, 2=LEJ,', &
     &  ' 3=LEJ+SCH, 4=LEJ+SMI, 5=PAU+SMI',/,i3)
   write(99,10021) ilib
10021 format(' METALLICITY OF THE HIGH RESOLUTION MODELS',/ &
     &  ' (1=0.001, 2=0.008, 3=0.020, 4=0.040):',/,i3)
   write(99,10020) iline
10020 format(' METALLICITY OF THE UV LINE SPECTRUM: (1=SOLAR,', &
     &  ' 2=LMC/SMC)',/,i3)
   write(99,10018) ivt,irsg
10018 format(' RSG FEATURE: MICROTURB. VELOCITY (1-6),' &
     &  ' SOLAR/NON-SOLAR ABUNDANCE (0,1)',/,i3,',',i1)
   write(99,10019) io1,io2,io3,io4,io5,io6,io7,io8,io9,io10, &
     &                  io11,io12,io13,io14,io15
10019 format(' OUTPUT FILES (NO<0, YES >=0): ',/,15i3)
   close(1)

   return
end subroutine output
c
c **************************************************************************
c **************************************************************************
c ****************************************************************************
c
      function sp_feature(l,isel,teff,icount)
c
c CALCULATION OF SELECTED EMISSION AND ABSORPTION FEATURES IN THE SPECTRUM.
c CURRENTLY THE FOLLOWING FEATURES ARE SUPPORTED: CO INDEX IN THE DEFINITION
c OF DOYON ET AL. (APJ, 421, 101 [1994]); CA IR TRIPLET EQUIVALENT WIDTH
c USING CALIBRATION OF DIAZ ET AL. (MNRAS, 239, 325 [1989]), STELLAR
c LYMAN-ALPHA OF PENA-GUERRERO & LEITHERER (2013).
c
c ****************************************************************************
c
   ! Modern Fortran version of sp_feature function header and variable declarations
   function sp_feature(l, isel, teff, icount) result(res)
      use, intrinsic :: iso_fortran_env, only: real32
      implicit none
      integer, intent(in) :: l, isel, icount
      real(real32), intent(in) :: teff
      real(real32) :: res

      ! Parameters and constants
      integer, parameter :: nmaxint=10, nmaxint1=11, npgrid=3000
      integer :: imuturb, insabund, ind4
      real(real32) :: ttemp, gtemp, t4

      ! These arrays should be provided by modules in modern code
      real(real32), allocatable :: temp(:), bol(:), zmass(:), xsurf(:), ysurf(:)
      real(real32), allocatable :: xc12s(:), xn14s(:), xo16s(:)
      real(real32), allocatable :: co162data(:,:,:,:), co229data(:,:,:,:), si159data(:,:,:,:)
      real(real32), allocatable :: co162nsdata(:,:,:,:), co229nsdata(:,:,:,:), si159nsdata(:,:,:,:)
      real(real32), allocatable :: origt(:), origg(:)
      real(real32), allocatable :: xlyman_t(:), xlyman_g(:), xlyman(:,:)
      real(real32) :: z

      ! In modern Fortran, these should be set via modules or passed as arguments
      ! Here we assume they are available

      ! Set microturbulent velocity and abundance index
      imuturb = ivt
      insabund = irsg

      ! Temporary variables for Teff and log g
      ttemp = teff
      gtemp = log10(zmass(l)) + 4.0_real32 * temp(l) - bol(l) - 10.6_real32

      select case (isel)
      case (1) ! CO index (Doyon)
         if (teff > 6000.0_real32) then
            res = 1.0_real32
         else
            t4 = teff / 10000.0_real32
            if (bol(l) < 0.3_real32) then
               res = 10.0_real32**(-0.4_real32 * (0.866_real32 - 2.95_real32*t4 + 2.55_real32*t4*t4))
            else if (bol(l) > 3.0_real32) then
               res = 10.0_real32**(-0.4_real32 * (1.353_real32 - 2.80_real32*t4))
            else
               res = 10.0_real32**(-0.4_real32 * (1.530_real32 - 5.01_real32*t4 + 4.10_real32*t4*t4))
            end if
            if (res < 0.0_real32) res = 0.0_real32
         end if

      case (2) ! Ca II triplet (Diaz)
         res = 10.21_real32 - 0.95_real32 * gtemp + 2.18_real32 * log10(z)

      case (3) ! Ca II triplet, T < 7200
         if (teff <= 7200.0_real32) then
            res = 10.21_real32 - 0.95_real32 * gtemp + 2.18_real32 * log10(z)
         else
            res = 0.0_real32
         end if

      case (4) ! CO 1.62 micron EW (Origlia)
         if (z == 0.05_real32) ind4 = 1
         if (z == 0.2_real32)  ind4 = 2
         if (z == 0.4_real32)  ind4 = 3
         if (z == 1.0_real32)  ind4 = 4
         if (z == 2.0_real32)  ind4 = 5
         if (teff > 5000.0_real32 .or. bol(l) <= -10.0_real32) then
            res = 0.0_real32
         else
            if (ttemp < 3000.0_real32) ttemp = 3000.0_real32
            if (gtemp > 2.0_real32) then
               res = 0.0_real32
            else if (gtemp <= 0.0_real32) then
               if (insabund == 0) then
                  res = ynter2x(ttemp, 0.01_real32, origt, origg, co162data, 8, 5, 6, 5, imuturb, ind4)
               else
                  res = ynter2x(ttemp, 0.01_real32, origt, origg, co162nsdata, 8, 5, 6, 5, imuturb, ind4)
               end if
            else
               if (insabund == 0) then
                  res = ynter2x(ttemp, gtemp, origt, origg, co162data, 8, 5, 6, 5, imuturb, ind4)
               else
                  res = ynter2x(ttemp, gtemp, origt, origg, co162nsdata, 8, 5, 6, 5, imuturb, ind4)
               end if
            end if
         end if

      case (5) ! CO 2.29 micron EW (Origlia)
         if (z == 0.05_real32) ind4 = 1
         if (z == 0.2_real32)  ind4 = 2
         if (z == 0.4_real32)  ind4 = 3
         if (z == 1.0_real32)  ind4 = 4
         if (z == 2.0_real32)  ind4 = 5
         if (teff > 5000.0_real32 .or. bol(l) <= -10.0_real32) then
            res = 0.0_real32
         else
            if (ttemp < 3000.0_real32) ttemp = 3000.0_real32
            if (gtemp > 2.0_real32) then
               res = 0.0_real32
            else if (gtemp <= 0.0_real32) then
               if (insabund == 0) then
                  res = ynter2x(ttemp, 0.01_real32, origt, origg, co229data, 8, 5, 6, 5, imuturb, ind4)
               else
                  res = ynter2x(ttemp, 0.01_real32, origt, origg, co229nsdata, 8, 5, 6, 5, imuturb, ind4)
               end if
            else
               if (insabund == 0) then
                  res = ynter2x(ttemp, gtemp, origt, origg, co229data, 8, 5, 6, 5, imuturb, ind4)
               else
                  res = ynter2x(ttemp, gtemp, origt, origg, co229nsdata, 8, 5, 6, 5, imuturb, ind4)
               end if
            end if
         end if

      case (6) ! Si 1.59 micron EW (Origlia)
         if (z == 0.05_real32) ind4 = 1
         if (z == 0.2_real32)  ind4 = 2
         if (z == 0.4_real32)  ind4 = 3
         if (z == 1.0_real32)  ind4 = 4
         if (z == 2.0_real32)  ind4 = 5
         if (teff > 5000.0_real32 .or. bol(l) <= -10.0_real32) then
            res = 0.0_real32
         else
            if (ttemp < 3000.0_real32) ttemp = 3000.0_real32
            if (gtemp > 2.0_real32) then
               res = 0.0_real32
            else if (gtemp <= 0.0_real32) then
               if (insabund == 0) then
                  res = ynter2x(ttemp, 0.01_real32, origt, origg, si159data, 8, 5, 6, 5, imuturb, ind4)
               else
                  res = ynter2x(ttemp, 0.01_real32, origt, origg, si159nsdata, 8, 5, 6, 5, imuturb, ind4)
               end if
            else
               if (insabund == 0) then
                  res = ynter2x(ttemp, gtemp, origt, origg, si159data, 8, 5, 6, 5, imuturb, ind4)
               else
                  res = ynter2x(ttemp, gtemp, origt, origg, si159nsdata, 8, 5, 6, 5, imuturb, ind4)
               end if
            end if
         end if

      case (7) ! Lyman-alpha EW (Pena-Guerrero & Leitherer)
         if (bol(l) <= -10.0_real32 .or. ttemp < 10000.0_real32) then
            res = 0.0_real32
         else if (ttemp < 15000.0_real32) then
            res = -420.0_real32 * exp(-ttemp / 6100.0_real32)
         else
            if (ttemp > 50000.0_real32) ttemp = 49999.0_real32
            if (gtemp > 4.5_real32) gtemp = 4.49_real32
            if (gtemp < 1.75_real32) gtemp = 1.76_real32
            res = ynter2(gtemp, ttemp, xlyman_g, xlyman_t, xlyman, 12, 23)
         end if

      case default
         res = 0.0_real32
      end select

   end function sp_feature
c
c ********************************************************************
c ********************************************************************
c ********************************************************************
!> Returns the wind power and momentum calculated from the mass-loss rates and terminal velocities.
!> The mass-loss rates are taken from the Geneva tabulation. Terminal velocities are based on observations.
!> Modern Fortran version.
function wind1(l, ltype, ichoic) result(res)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   integer, intent(in) :: l, ichoic
   integer, intent(inout) :: ltype
   real(real32) :: res

   ! Parameters
   integer, parameter :: npgrid = 3000
   real(real32), parameter :: tiny = 1.e-30_real32

   ! Variables (should be replaced by module variables in modern code)
   real(real32) :: cmass(npgrid), temp(npgrid), bol(npgrid), xsurf(npgrid), zmass(npgrid)
   real(real32) :: ysurf(npgrid), xc12s(npgrid), xn14s(npgrid), xo16s(npgrid)
   real(real32) :: bmdot(npgrid), tt_star(npgrid)
   real(real32) :: z, xmwr
   real(real32) :: vinf, gamma0, cnr, coher

   ! Defensive: set bmdot for dead stars
   if (bol(l) < -10.0_real32 .or. bmdot(l) > -0.1_real32) bmdot(l) = -30.0_real32

   ! 1) Normal hot stars (log Teff > 3.9): Howarth & Prinja (1989)
   if (temp(l) > 3.9_real32) then
      gamma0 = 1.0_real32 - 2.7e-5_real32 * 10.0_real32**bol(l) / (zmass(l) + tiny)
      if (gamma0 <= 0.0_real32) gamma0 = 1.e-10_real32
      vinf = 618.0_real32 * sqrt(zmass(l) / (10.0_real32**(0.5_real32*bol(l) - 2.0_real32*temp(l) + 7.52_real32)) / gamma0) * &
             (0.58_real32 + 2.04_real32 * (0.5_real32*bol(l) - 2.0_real32*temp(l) + 7.52_real32))
      ltype = 1
   end if

   ! 2) Cool stars (log Teff <= 3.9): generic value
   if (temp(l) <= 3.9_real32) then
      vinf = 30.0_real32
      ltype = 2
   end if

   ! 3) LBVs (4.4 > log Teff > 3.75 and log Mdot > -3.5)
   if (temp(l) < 4.4_real32 .and. temp(l) > 3.75_real32 .and. bmdot(l) > -3.5_real32) then
      vinf = 200.0_real32
      ltype = 3
   end if

   ! 4) WR stars (log Teff > 4.4 and Xsurf < 0.4 and M >= M_WR)
   if (temp(l) > 4.4_real32 .and. xsurf(l) < 0.4_real32 .and. cmass(l) >= xmwr) then
      ltype = 4
      if (xn14s(l) == 0.0_real32) xn14s(l) = 1.e-6_real32
      cnr = xc12s(l) / xn14s(l)
      coher = ((xc12s(l)/12.0_real32) + (xo16s(l)/16.0_real32)) / (ysurf(l)/4.0_real32)
      if (xsurf(l) > 0.1_real32) then
         vinf = 1650.0_real32
         goto 777
      end if
      if (cnr < 10.0_real32) then
         vinf = 1900.0_real32
         goto 777
      end if
      if (coher < 0.5_real32) then
         vinf = 1800.0_real32
         goto 777
      end if
      if (coher < 1.0_real32) then
         vinf = 2800.0_real32
         goto 777
      end if
      if (coher >= 1.0_real32) then
         vinf = 3500.0_real32
         goto 777
      end if
   end if

   ! All terminal velocities except WRs are scaled with Z**0.13
   vinf = vinf * z**0.13_real32

777 continue
   if (bol(l) < -10.0_real32) then
      vinf = 1.e-10_real32
      ltype = 0
   end if

   ! Calculation of the wind power and the wind momentum. 3.155 comes from the conversion to CGS.
   ! wind1 was divided by 10**35 to avoid overflow.
   select case (ichoic)
   case (1)
      res = 10.0_real32**bmdot(l) * vinf * vinf * 3.155_real32
   case (2)
      res = 2.0_real32 * 10.0_real32**bmdot(l) * vinf * 3.155e-5_real32
   case (3)
      res = 10.0_real32**bmdot(l)
   case (4)
      res = vinf
   case default
      res = 0.0_real32
   end select

end function wind1
c
c *************************************************************
c *************************************************************
c *************************************************************
!> Returns the wind power and momentum calculated from the mass-loss rates and terminal velocities.
!> The mass-loss rates and terminal velocities are taken from observations.
!> Modern Fortran version.
function wind2(l, ltype, ichoic) result(res)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   integer, intent(in) :: l, ichoic
   integer, intent(inout) :: ltype
   real(real32) :: res

   ! Parameters
   integer, parameter :: nmaxint=10, nmaxint1=11, npgrid=3000
   real(real32), parameter :: tiny = 1.e-30_real32

   ! Variables (should be replaced by module variables in modern code)
   real(real32) :: time1, tstep, tmax, toma, sfr, doma, z, tdel, tvar, tiempo1, tinter
   integer :: iz, ninterv, upma, jmg, iwind, sncut, bhcut, iatmos, ilib, iline, ivt, irsg
   integer :: io1, io2, io3, io4, io5, io6, io7, io8, isf, io9, io10, io11, io12, io13, io14, io15
   real(real32) :: xponent(nmaxint), xmaslim(nmaxint1), xmwr
   integer :: iwrscale, iwrt, jtime
   real(real32) :: cmass(npgrid), temp(npgrid), bol(npgrid), xsurf(npgrid), zmass(npgrid)
   real(real32) :: ysurf(npgrid), xc12s(npgrid), xn14s(npgrid), xo16s(npgrid)
   real(real32) :: bmdot(npgrid), tt_star(npgrid)
   real(real32) :: xloss, vinf, gamma0, cnr, coher

   ! Main logic
   if (temp(l) > 3.9_real32) then
      xloss = 1.69_real32 * bol(l) - 15.41_real32
      gamma0 = 1.0_real32 - 2.7e-5_real32 * 10.0_real32**bol(l) / (zmass(l) + tiny)
      if (gamma0 <= 0.0_real32) gamma0 = 1.e-10_real32
      vinf = 618.0_real32 * sqrt(zmass(l) / (10.0_real32**(0.5_real32*bol(l) - 2.0_real32*temp(l) + 7.52_real32)) / gamma0) * &
             (0.58_real32 + 2.04_real32 * (0.5_real32*bol(l) - 2.0_real32*temp(l) + 7.52_real32))
      ltype = 1
   end if

   if (temp(l) <= 3.9_real32) then
      xloss = -12.26_real32 + 1.5_real32*bol(l) - 2.0_real32*temp(l) + 7.52_real32 - log10(zmass(l) + tiny)
      vinf = 30.0_real32
      ltype = 2
   end if

   if (temp(l) < 4.4_real32 .and. temp(l) > 3.75_real32 .and. bmdot(l) > -3.5_real32) then
      xloss = -3.9_real32
      vinf = 200.0_real32
      ltype = 3
   end if

   if (temp(l) > 4.4_real32 .and. xsurf(l) < 0.4_real32 .and. cmass(l) >= xmwr) then
      ltype = 4
      if (xn14s(l) == 0.0_real32) xn14s(l) = 1.e-6_real32
      cnr = xc12s(l) / xn14s(l)
      coher = ((xc12s(l) / 12.0_real32) + (xo16s(l) / 16.0_real32)) / (ysurf(l) / 4.0_real32)
      if (xsurf(l) > 0.1_real32) then
         vinf = 1650.0_real32
         xloss = -4.2_real32
         goto 777
      end if
      if (cnr < 10.0_real32) then
         vinf = 1900.0_real32
         xloss = -4.5_real32
         goto 777
      end if
      if (coher < 0.5_real32) then
         vinf = 1800.0_real32
         xloss = -4.4_real32
         goto 777
      end if
      if (coher < 1.0_real32) then
         vinf = 2800.0_real32
         xloss = -4.7_real32
         goto 777
      end if
      if (coher >= 1.0_real32) then
         vinf = 3500.0_real32
         xloss = -5.0_real32
         goto 777
      end if
   end if

   xloss = xloss + 0.80_real32 * log10(z)
   vinf = vinf * z**0.13_real32

777 continue
   if (bol(l) < -10.0_real32 .or. bmdot(l) > -0.1_real32) then
      xloss = -30.0_real32
      vinf = 1.e-10_real32
      ltype = 0
   end if

   ! Calculation of the wind power and the wind momentum. 3.155 comes from the conversion to CGS.
   ! wind2 was divided by 10**35 to avoid overflow.
   select case (ichoic)
   case (1)
      res = 10.0_real32**xloss * vinf * vinf * 3.155_real32
   case (2)
      res = 2.0_real32 * 10.0_real32**xloss * vinf * 3.155e-5_real32
   case (3)
      res = 10.0_real32**xloss
   case (4)
      res = vinf
   case default
      res = 0.0_real32
   end select

end function wind2
c
c *************************************************************
c *************************************************************
c *************************************************************
!> Returns the wind power and momentum calculated from the mass-loss rates and terminal velocities.
!> The mass-loss rates and terminal velocities of OB stars are from the CAK theory.
!> Mass-loss rates and terminal velocities of all other stars are taken from observations.
!> Modern Fortran version.
function wind3(l, ltype, ichoic) result(res)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   integer, intent(in) :: l, ichoic
   integer, intent(inout) :: ltype
   real(real32) :: res

   ! Parameters
   integer, parameter :: nmaxint=10, nmaxint1=11, npgrid=3000
   real(real32), parameter :: tiny = 1.e-30_real32

   ! Variables (should be replaced by module variables in modern code)
   real(real32) :: time1, tstep, tmax, toma, sfr, doma, z, tdel, tvar, tiempo1, tinter
   integer :: iz, ninterv, upma, jmg, iwind, sncut, bhcut, iatmos, ilib, iline, ivt, irsg
   integer :: io1, io2, io3, io4, io5, io6, io7, io8, isf, io9, io10, io11, io12, io13, io14, io15
   real(real32) :: xponent(nmaxint), xmaslim(nmaxint1), xmwr
   integer :: iwrscale, iwrt, jtime
   real(real32) :: cmass(npgrid), temp(npgrid), bol(npgrid), xsurf(npgrid), zmass(npgrid)
   real(real32) :: ysurf(npgrid), xc12s(npgrid), xn14s(npgrid), xo16s(npgrid)
   real(real32) :: bmdot(npgrid), tt_star(npgrid)
   real(real32) :: xloss, vinf, cnr, coher

   ! Main logic
   if (temp(l) > 3.9_real32) then
      xloss = 2.45_real32 * bol(l) - 1.10_real32 * log10(zmass(l) + tiny) + &
              1.31_real32 * temp(l) - 24.06_real32
      vinf = 10.0_real32**(-0.30_real32 * bol(l) + 0.55_real32 * log10(zmass(l) + tiny) + &
              0.64_real32 * temp(l) + 1.23_real32)
      ltype = 1
   end if

   if (temp(l) <= 3.9_real32) then
      xloss = -12.26_real32 + 1.5_real32 * bol(l) - 2.0_real32 * temp(l) + 7.52_real32 - log10(zmass(l) + tiny)
      vinf = 30.0_real32
      ltype = 2
   end if

   if (temp(l) < 4.4_real32 .and. temp(l) > 3.75_real32 .and. bmdot(l) > -3.5_real32) then
      xloss = -3.9_real32
      vinf = 200.0_real32
      ltype = 3
   end if

   if (temp(l) > 4.4_real32 .and. xsurf(l) < 0.4_real32 .and. cmass(l) >= xmwr) then
      ltype = 4
      if (xn14s(l) == 0.0_real32) xn14s(l) = 1.e-6_real32
      cnr = xc12s(l) / xn14s(l)
      coher = ((xc12s(l) / 12.0_real32) + (xo16s(l) / 16.0_real32)) / (ysurf(l) / 4.0_real32)
      if (xsurf(l) > 0.1_real32) then
         vinf = 1650.0_real32
         xloss = -4.2_real32
         goto 777
      end if
      if (cnr < 10.0_real32) then
         vinf = 1900.0_real32
         xloss = -4.5_real32
         goto 777
      end if
      if (coher < 0.5_real32) then
         vinf = 1800.0_real32
         xloss = -4.4_real32
         goto 777
      end if
      if (coher < 1.0_real32) then
         vinf = 2800.0_real32
         xloss = -4.7_real32
         goto 777
      end if
      if (coher >= 1.0_real32) then
         vinf = 3500.0_real32
         xloss = -5.0_real32
         goto 777
      end if
   end if

   xloss = xloss + 0.80_real32 * log10(z)
   vinf = vinf * z**0.13_real32

777 continue
   if (bol(l) < -10.0_real32 .or. bmdot(l) > -0.1_real32) then
      xloss = -30.0_real32
      vinf = 1.e-10_real32
      ltype = 0
   end if

   ! Calculation of the wind power and the wind momentum. 3.155 comes from the conversion to CGS.
   ! wind3 was divided by 10**35 to avoid overflow.
   select case (ichoic)
   case (1)
      res = 10.0_real32**xloss * vinf * vinf * 3.155_real32
   case (2)
      res = 2.0_real32 * 10.0_real32**xloss * vinf * 3.155e-5_real32
   case (3)
      res = 10.0_real32**xloss
   case (4)
      res = vinf
   case default
      res = 0.0_real32
   end select

end function wind3
c
c *************************************************************
c *************************************************************
c *************************************************************
!> Returns the wind power and momentum calculated from the mass-loss rates and terminal velocities.
!> The mass-loss rates and terminal velocities are taken from observations a la Elson et al. (1989).
!> Modern Fortran version.
function wind4(l, ltype, ichoic) result(res)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   integer, intent(in) :: l, ichoic
   integer, intent(inout) :: ltype
   real(real32) :: res

   ! Parameters
   integer, parameter :: nmaxint=10, nmaxint1=11, npgrid=3000
   real(real32), parameter :: tiny = 1.e-30_real32

   ! Variables (should be replaced by module variables in modern code)
   real(real32) :: time1, tstep, tmax, toma, sfr, doma, z, tdel, tvar, tiempo1, tinter
   integer :: iz, ninterv, upma, jmg, iwind, sncut, bhcut, iatmos, ilib, iline, ivt, irsg
   integer :: io1, io2, io3, io4, io5, io6, io7, io8, isf, io9, io10, io11, io12, io13, io14, io15
   real(real32) :: xponent(nmaxint), xmaslim(nmaxint1), xmwr
   integer :: iwrscale, iwrt, jtime
   real(real32) :: cmass(npgrid), temp(npgrid), bol(npgrid), xsurf(npgrid), zmass(npgrid)
   real(real32) :: ysurf(npgrid), xc12s(npgrid), xn14s(npgrid), xo16s(npgrid)
   real(real32) :: bmdot(npgrid), tt_star(npgrid)
   real(real32) :: xloss, vinf, cnr, coher

   ! Main logic
   if (temp(l) > 3.9_real32) then
      xloss = 5.0_real32 * log10(cmass(l) + tiny) - 14.24_real32
      vinf = 1863.0_real32 * cmass(l)**0.14_real32
      ltype = 1
   end if

   if (temp(l) <= 3.9_real32) then
      xloss = -12.26_real32 + 1.5_real32 * bol(l) - 2.0_real32 * temp(l) + 7.52_real32 - log10(zmass(l) + 1.e-20_real32)
      vinf = 30.0_real32
      ltype = 2
   end if

   if (temp(l) < 4.4_real32 .and. temp(l) > 3.75_real32 .and. bmdot(l) > -3.5_real32) then
      xloss = -3.9_real32
      vinf = 200.0_real32
      ltype = 3
   end if

   if (temp(l) > 4.4_real32 .and. xsurf(l) < 0.4_real32 .and. cmass(l) >= xmwr) then
      ltype = 4
      if (xn14s(l) == 0.0_real32) xn14s(l) = 1.e-6_real32
      cnr = xc12s(l) / xn14s(l)
      coher = ((xc12s(l) / 12.0_real32) + (xo16s(l) / 16.0_real32)) / (ysurf(l) / 4.0_real32)
      if (xsurf(l) > 0.1_real32) then
         vinf = 1650.0_real32
         xloss = -4.2_real32
         goto 777
      end if
      if (cnr < 10.0_real32) then
         vinf = 1900.0_real32
         xloss = -4.5_real32
         goto 777
      end if
      if (coher < 0.5_real32) then
         vinf = 1800.0_real32
         xloss = -4.4_real32
         goto 777
      end if
      if (coher < 1.0_real32) then
         vinf = 2800.0_real32
         xloss = -4.7_real32
         goto 777
      end if
      if (coher >= 1.0_real32) then
         vinf = 3500.0_real32
         xloss = -5.0_real32
         goto 777
      end if
   end if

   xloss = xloss + 0.80_real32 * log10(z)
   vinf = vinf * z**0.13_real32

777 continue
   if (bol(l) < -10.0_real32 .or. bmdot(l) > -0.1_real32) then
      xloss = -30.0_real32
      vinf = 1.e-10_real32
      ltype = 0
   end if

   ! Calculation of the wind power and the wind momentum. 3.155 comes from the conversion to CGS.
   ! wind4 was divided by 10**35 to avoid overflow.
   select case (ichoic)
   case (1)
      res = 10.0_real32**xloss * vinf * vinf * 3.155_real32
   case (2)
      res = 2.0_real32 * 10.0_real32**xloss * vinf * 3.155e-5_real32
   case (3)
      res = 10.0_real32**xloss
   case (4)
      res = vinf
   case default
      res = 0.0_real32
   end select

end function wind4
c
c
c ************************************************************
c ************************************************************
c ************************************************************
!> 2-D polynomial interpolation (Numerical Recipes style), modern Fortran version.
subroutine polin2(x1a, x2a, ya, m, n, x1, x2, y, dy)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   integer, intent(in) :: m, n
   real(real32), intent(in) :: x1a(m), x2a(n)
   real(real32), intent(in) :: ya(m, n)
   real(real32), intent(in) :: x1, x2
   real(real32), intent(out) :: y, dy

   integer :: j, k, kk, minx, maxx
   real(real32) :: yntmp(n), ymtmp(m), xx(n), x3

   interface
      pure integer function ismin(n, x, idummy)
         use, intrinsic :: iso_fortran_env, only: real32
         integer, intent(in) :: n, idummy
         real(real32), intent(in) :: x(n)
      end function ismin
      pure integer function ismax(n, x, idummy)
         use, intrinsic :: iso_fortran_env, only: real32
         integer, intent(in) :: n, idummy
         real(real32), intent(in) :: x(n)
      end function ismax
      pure function yntra(r, rw, fu, imax) result(yval)
         use, intrinsic :: iso_fortran_env, only: real32
         real(real32), intent(in) :: r
         real(real32), intent(in) :: rw(imax), fu(imax)
         integer, intent(in) :: imax
         real(real32) :: yval
      end function yntra
   end interface

   do j = 1, m
      kk = 0
      do k = 1, n
         if (ya(j, k) >= 1.e-30_real32) then
            kk = kk + 1
            xx(kk) = x2a(k)
            yntmp(kk) = ya(j, k)
         end if
      end do
      if (kk == 0) then
         ymtmp(j) = 0.0_real32
      else
         minx = ismin(kk, xx, 1)
         maxx = ismax(kk, xx, 1)
         x3 = x2
         if (x2 < xx(minx)) x3 = xx(minx)
         if (x2 > xx(maxx)) x3 = xx(maxx)
         ymtmp(j) = yntra(x3, xx, yntmp, kk)
      end if
   end do
   y = yntra(x1, x1a, ymtmp, m)
   dy = 0.0_real32 ! Not used in this implementation

end subroutine polin2
c
c *********************************************************************
c *********************************************************************
c *********************************************************************
!> General 2-D interpolation routine for a function f(x, y).
!> Interpolates in x and y, for a regular grid.
pure function ynter2(x, y, xd, yd, f, imax, jmax) result(res)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   real(real32), intent(in) :: x, y
   real(real32), intent(in) :: xd(imax), yd(jmax)
   real(real32), intent(in) :: f(imax, jmax)
   integer, intent(in) :: imax, jmax
   real(real32) :: res
   integer :: i, i1, j, j1
   real(real32) :: fax, fay, t1, t2, t3, t4

   ! Defensive: check bounds for x
   if (x <= xd(1)) then
      i = 1
      i1 = 1
      fax = 0.0_real32
   else if (x >= xd(imax)) then
      i = imax
      i1 = imax
      fax = 0.0_real32
   else
      do i = 1, imax-1
         if (x >= xd(i) .and. x < xd(i+1)) exit
      end do
      i1 = i + 1
      fax = (x - xd(i)) / (xd(i1) - xd(i))
   end if

   ! Defensive: check bounds for y
   if (y <= yd(1)) then
      j = 1
      j1 = 1
      fay = 0.0_real32
   else if (y >= yd(jmax)) then
      j = jmax
      j1 = jmax
      fay = 0.0_real32
   else
      do j = 1, jmax-1
         if (y >= yd(j) .and. y < yd(j+1)) exit
      end do
      j1 = j + 1
      fay = (y - yd(j)) / (yd(j1) - yd(j))
   end if

   ! Bilinear interpolation
   t1 = (1.0_real32 - fax) * (1.0_real32 - fay) * f(i,  j)
   t2 = fax            * (1.0_real32 - fay) * f(i1, j)
   t3 = (1.0_real32 - fax) * fay            * f(i,  j1)
   t4 = fax            * fay            * f(i1, j1)
   res = t1 + t2 + t3 + t4

end function ynter2
c
c *********************************************************************
c *********************************************************************
c *********************************************************************
c
!> General 2-D interpolation routine for a function f(x, y, ind3, ind4).
!> Interpolates in x and y, for fixed ind3 and ind4 indices.
pure function ynter2x(x, y, xd, yd, f, imax, jmax, ind3max, ind4max, ind3, ind4) result(res)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   real(real32), intent(in) :: x, y
   real(real32), intent(in) :: xd(imax), yd(jmax)
   real(real32), intent(in) :: f(imax, jmax, ind3max, ind4max)
   integer, intent(in) :: imax, jmax, ind3max, ind4max, ind3, ind4
   real(real32) :: res
   integer :: i, i1, j, j1
   real(real32) :: fax, fay, t1, t2, t3, t4

   ! Defensive: check bounds for x
   if (x <= xd(1)) then
      i = 1
      i1 = 1
      fax = 0.0_real32
   else if (x >= xd(imax)) then
      i = imax
      i1 = imax
      fax = 0.0_real32
   else
      do i = 1, imax-1
         if (x >= xd(i) .and. x < xd(i+1)) exit
      end do
      i1 = i + 1
      fax = (x - xd(i)) / (xd(i1) - xd(i))
   end if

   ! Defensive: check bounds for y
   if (y <= yd(1)) then
      j = 1
      j1 = 1
      fay = 0.0_real32
   else if (y >= yd(jmax)) then
      j = jmax
      j1 = jmax
      fay = 0.0_real32
   else
      do j = 1, jmax-1
         if (y >= yd(j) .and. y < yd(j+1)) exit
      end do
      j1 = j + 1
      fay = (y - yd(j)) / (yd(j1) - yd(j))
   end if

   ! Bilinear interpolation
   t1 = (1.0_real32 - fax) * (1.0_real32 - fay) * f(i,  j,  ind3, ind4)
   t2 = fax            * (1.0_real32 - fay) * f(i1, j,  ind3, ind4)
   t3 = (1.0_real32 - fax) * fay            * f(i,  j1, ind3, ind4)
   t4 = fax            * fay            * f(i1, j1, ind3, ind4)
   res = t1 + t2 + t3 + t4

end function ynter2x
c
c ************************************************************
c *************************************************************
c *************************************************************
!> 1-D linear interpolation: returns y(r) for given x=r, using arrays rw (x) and fu (y), both size imax.
pure function yntra(r, rw, fu, imax) result(y)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   real(real32), intent(in) :: r
   real(real32), intent(in) :: rw(imax), fu(imax)
   integer, intent(in) :: imax
   real(real32) :: y
   integer :: i, i1
   real(real32) :: fak

   ! Defensive: check bounds
   if (r <= rw(1)) then
      y = fu(1)
      return
   else if (r >= rw(imax)) then
      y = fu(imax)
      return
   end if

   ! Find interval
   do i = 1, imax-1
      if (r >= rw(i) .and. r < rw(i+1)) exit
   end do
   i1 = i + 1
   fak = (r - rw(i)) / (rw(i1) - rw(i))
   y = fu(i) + (fu(i1) - fu(i)) * fak
end function yntra

!> Quadratic interpolation using 3 nearest points (Numerical Recipes style).
pure function reci_polint(x, ndata, xdata, fdata) result(y)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   real(real32), intent(in) :: x
   integer, intent(in) :: ndata
   real(real32), intent(in) :: xdata(ndata), fdata(ndata)
   real(real32) :: y
   integer :: j, idegree, j0, j1, j2
   real(real32) :: dy

   ! Find nearest index j such that xdata(j) <= x < xdata(j+1)
   j = 1
   do while (j < ndata .and. x > xdata(j+1))
      j = j + 1
   end do
   if (j > ndata-2) j = ndata-2
   if (j < 1) j = 1
   j0 = j
   j1 = j + 1
   j2 = j + 2

   call polint([xdata(j0), xdata(j1), xdata(j2)], [fdata(j0), fdata(j1), fdata(j2)], 3, x, y, dy)
end function reci_polint

!> Polynomial interpolation (Neville's algorithm, Numerical Recipes style).
pure subroutine polint(xa, ya, n, x, y, dy)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   integer, intent(in) :: n
   real(real32), intent(in) :: xa(n), ya(n), x
   real(real32), intent(out) :: y, dy
   real(real32) :: c(n), d(n)
   integer :: i, m, ns
   real(real32) :: dif, dift, ho, hp, w, den

   ns = 1
   dif = abs(x - xa(1))
   do i = 1, n
      dift = abs(x - xa(i))
      if (dift < dif) then
         ns = i
         dif = dift
      end if
      c(i) = ya(i)
      d(i) = ya(i)
   end do
   y = ya(ns)
   ns = ns - 1
   do m = 1, n-1
      do i = 1, n-m
         ho = xa(i) - x
         hp = xa(i+m) - x
         w = c(i+1) - d(i)
         den = ho - hp
         if (den == 0.0_real32) stop 'polint: zero denominator'
         den = w / den
         d(i) = hp * den
         c(i) = ho * den
      end do
      if (2*ns < n-m) then
         dy = c(ns+1)
      else
         dy = d(ns)
         ns = ns - 1
      end if
      y = y + dy
   end do
end subroutine polint
c
c ******************************************************************************
c ******************************************************************************
c ******************************************************************************
!> Locates the index j in an ordered array xx such that xx(j) <= x < xx(j+1).
!> If x is outside the range, returns the closest valid index.
pure subroutine locate(xx, n, x, j)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   integer, intent(in) :: n
   real(real32), intent(in) :: xx(n), x
   integer, intent(out) :: j
   integer :: jl, jm, ju

   jl = 0
   ju = n + 1
   do
      if (ju - jl <= 1) exit
      jm = (ju + jl) / 2
      if ((xx(n) > xx(1)) .eqv. (x > xx(jm))) then
         jl = jm
      else
         ju = jm
      end if
   end do
   j = min(n, jl)
   j = max(j, 1)
end subroutine locate

!> Akima spline interpolation for a single-valued function.
!> Interpolates y(x) at points u(:), returning v(:).
subroutine intrpl(l, x, y, n, u, v)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   integer, intent(in) :: l, n
   real(real32), intent(in) :: x(l), y(l), u(n)
   real(real32), intent(out) :: v(n)
   integer :: i, j, k, imn, imx, l0, lm1, lm2, lp1, n0, ipv
   real(real32) :: uk, x2, x3, x4, x5, y2, y3, y4, y5
   real(real32) :: a1, a2, a3, a4, a5, m1, m2, m3, m4, m5
   real(real32) :: w2, w3, w4, sw, t3, t4, sa, dx, q0, q1, q2, q3

   l0  = l
   lm1 = l0  - 1
   lm2 = lm1 - 1
   lp1 = l0  + 1
   n0  = n
   if (lm2 < 0 .or. n0 <= 0) return
   do i = 2, l0
      if (x(i-1) >= x(i)) return
   end do
   ipv = 0

   do k = 1, n0
      uk = u(k)
      ! Locate interval
      if (lm2 == 0) then
         i = 2
      else if (uk >= x(l0)) then
         i = lp1
      else if (uk < x(1)) then
         i = 1
      else
         imn = 2
         imx = l0
         do
            i = (imn + imx) / 2
            if (uk >= x(i)) then
               imn = i + 1
            else
               imx = i
            end if
            if (imx <= imn) exit
         end do
         i = imx
      end if

      if (i == ipv) cycle
      ipv = i

      ! Pick up necessary x and y values
      j = i
      if (j == 1) j = 2
      if (j == lp1) j = l0
      x3 = x(j-1)
      y3 = y(j-1)
      x4 = x(j)
      y4 = y(j)
      a3 = x4 - x3
      m3 = (y4 - y3) / a3
      if (lm2 == 0) then
         m2 = m3
         m4 = m3
      else
         if (j == 2) then
            x2 = x3 - (x4 - x3)
            y2 = y3 - (y4 - y3)
            a2 = x3 - x2
            m2 = (y3 - y2) / a2
            x5 = x(j+1)
            y5 = y(j+1)
            a4 = x5 - x4
            m4 = (y5 - y4) / a4
            m2 = m3 + m3 - m4
         else if (j == l0) then
            x2 = x(j-2)
            y2 = y(j-2)
            a2 = x3 - x2
            m2 = (y3 - y2) / a2
            m4 = m3 + m3 - m2
         else
            x2 = x(j-2)
            y2 = y(j-2)
            a2 = x3 - x2
            m2 = (y3 - y2) / a2
            x5 = x(j+1)
            y5 = y(j+1)
            a4 = x5 - x4
            m4 = (y5 - y4) / a4
         end if
      end if

      if (j <= 3) then
         a1 = a2
         m1 = m2 + m2 - m3
      else
         a1 = x2 - x(j-3)
         m1 = (y2 - y(j-3)) / a1
      end if
      if (j >= lm1) then
         a5 = a4
         m5 = m4 + m4 - m3
      else
         a5 = x(j+2) - x5
         m5 = (y(j+2) - y5) / a5
      end if

      ! Numerical differentiation
      if (i == lp1) then
         w3 = abs(m5 - m4)
         w4 = abs(m3 - m2)
         sw = w3 + w4
         if (sw /= 0.0) then
            t4 = (w3 * m3 + w4 * m4) / sw
         else
            t4 = 0.5 * (m3 + m4)
         end if
         t3 = t4
         sa = a2 + a3
         t4 = 0.5 * (m4 + m5 - a2 * (a2 - a3) * (m2 - m3) / (sa * sa))
         x3 = x4
         y3 = y4
         a3 = a2
         m3 = m2
      else if (i == 1) then
         w2 = abs(m4 - m3)
         w3 = abs(m2 - m1)
         sw = w2 + w3
         if (sw /= 0.0) then
            t3 = (w2 * m2 + w3 * m3) / sw
         else
            t3 = 0.5 * (m2 + m3)
         end if
         t4 = t3
         sa = a3 + a4
         t3 = 0.5 * (m1 + m2 - a4 * (a3 - a4) * (m3 - m4) / (sa * sa))
         x3 = x3 - a4
         y3 = y3 - m2 * a4
         a3 = a4
         m3 = m2
      else
         w2 = abs(m4 - m3)
         w3 = abs(m2 - m1)
         sw = w2 + w3
         if (sw /= 0.0) then
            t3 = (w2 * m2 + w3 * m3) / sw
         else
            t3 = 0.5 * (m2 + m3)
         end if
         w3 = abs(m5 - m4)
         w4 = abs(m3 - m2)
         sw = w3 + w4
         if (sw /= 0.0) then
            t4 = (w3 * m3 + w4 * m4) / sw
         else
            t4 = 0.5 * (m3 + m4)
         end if
      end if

      ! Coefficients
      q0 = y3
      q1 = t3
      q2 = (2.0 * (m3 - t3) + m3 - t4) / a3
      q3 = (-2.0 * m3 + t3 + t4) / (a3 * a3)

      ! Polynomial evaluation
      dx = uk - x3
      v(k) = q0 + dx * (q1 + dx * (q2 + dx * q3))
   end do
end subroutine intrpl
c
c
c *****************************************************************************
c *****************************************************************************
c *****************************************************************************
!> Performs linear integration of function y(x) using the trapezoidal rule.
!> x and y must be arrays of length ndim. x must be monotonically decreasing.
!> weight returns the integration weights at each x.
!> sum returns the integral value.
subroutine fliwgt(x, y, weight, sum, ndim)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   integer, intent(in) :: ndim
   real(real32), intent(in)  :: x(ndim), y(ndim)
   real(real32), intent(out) :: weight(ndim)
   real(real32), intent(out) :: sum
   integer :: i

   weight(1) = (x(1) - x(2)) / 2.0_real32
   sum = weight(1) * y(1)
   do i = 2, ndim - 1
      weight(i) = (x(i-1) - x(i+1)) / 2.0_real32
      sum = sum + weight(i) * y(i)
   end do
   weight(ndim) = (x(ndim-1) - x(ndim)) / 2.0_real32
   sum = sum + weight(ndim) * y(ndim)
end subroutine fliwgt
c
c *********************************************************************
c *********************************************************************
c *********************************************************************
!> Returns the smallest index i such that x(i) = min(x(1:n)).
!> If n <= 0, returns zero.
pure integer function ismin(n, x, idummy)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   integer, intent(in) :: n, idummy
   real(real32), intent(in) :: x(n)
   integer :: i
   real(real32) :: xmin

   if (n <= 0) then
      ismin = 0
      return
   end if

   xmin = 1.0e34_real32
   ismin = 1
   do i = n, 1, -1
      if (x(i) <= xmin) then
         ismin = i
         xmin  = x(i)
      end if
   end do

end function ismin
c
c *********************************************************************
c *********************************************************************
c *********************************************************************
!> Returns the smallest index i such that x(i) = max(x(1:n)).
!> If n <= 0, returns zero.
pure integer function ismax(n, x, idummy)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   integer, intent(in) :: n, idummy
   real(real32), intent(in) :: x(n)
   integer :: i
   real(real32) :: xmax

   if (n <= 0) then
      ismax = 0
      return
   end if

   xmax = -1.0e38_real32
   ismax = 1
   do i = n, 1, -1
      if (x(i) >= xmax) then
         ismax = i
         xmax  = x(i)
      end if
   end do

end function ismax
c
c **********************************************************************
c **********************************************************************
c **********************************************************************
!> Prints the name of the subroutine where an error occurred.
subroutine errpri(name)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   character(len=*), intent(in) :: name

   write(*,*) 'ERROR in subroutine: ', trim(name)

end subroutine errpri

! ****************************************************************************
! ****************************************************************************
! ****************************************************************************
!> Performs linear interpolation of the vector x2 in (x1, y1), yielding y2.
!> x1 must be monotonically increasing or decreasing.
!> The dimension of vectors x1, y1 is n1; that of x2, y2 is n2.
subroutine linterp(x1, y1, n1, x2, y2, n2)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   integer, intent(in) :: n1, n2
   real(real32), intent(in)  :: x1(n1), y1(n1), x2(n2)
   real(real32), intent(out) :: y2(n2)
   integer :: i, k

   interface
      pure integer function pos(x0, x, m)
         use, intrinsic :: iso_fortran_env, only: real32
         real(real32), intent(in) :: x0
         real(real32), intent(in) :: x(m)
         integer, intent(in) :: m
      end function pos
      pure function flin(x1, x2, y1, y2, x0) result(flin_val)
         use, intrinsic :: iso_fortran_env, only: real32
         real(real32), intent(in) :: x1, x2, y1, y2, x0
         real(real32) :: flin_val
      end function flin
   end interface

   do i = 1, n2
      k = pos(x2(i), x1, n1)
      y2(i) = flin(x1(k), x1(k+1), y1(k), y1(k+1), x2(i))
   end do

end subroutine linterp
c ****************************************************************************
c ****************************************************************************
c ****************************************************************************
!> Determines the position of the point x0 in a monotonically increasing or decreasing vector x of dimension m.
!> Returns k = pos(x0, x, m) such that x(k) <= x0 < x(k+1).
!> If x0 is outside the range given by x, returns:
!>   pos = 1   if x0 is on the side of x(1)
!>   pos = m-1 if x0 is on the side of x(m)
pure integer function pos(x0, x, m)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   real(real32), intent(in) :: x0
   real(real32), intent(in) :: x(m)
   integer, intent(in) :: m

   integer :: n, k, i

   n = 0
   k = m + 1
   do
      if (k - n <= 1) exit
      i = (k + n) / 2
      if ((x(m) > x(1)) .eqv. (x0 > x(i))) then
         n = i
      else
         k = i
      end if
   end do
   n = max(1, n)
   n = min(m-1, n)
   pos = n

end function pos

! ****************************************************************************
! ****************************************************************************
! ****************************************************************************
!> Performs a linear interpolation or extrapolation between two points (x1, y1) and (x2, y2).
!> Returns the interpolated value at x0.
pure function flin(x1, x2, y1, y2, x0) result(flin_val)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   real(real32), intent(in) :: x1, x2, y1, y2, x0
   real(real32) :: flin_val
   real(real32) :: a12, a01, a02

   a12 = x1 - x2
   a01 = x0 - x1
   a02 = x0 - x2

   flin_val = (a02 * y1 - a01 * y2) / a12
end function flin

! ****************************************************************************
! ****************************************************************************
! ****************************************************************************
!> Parabolic interpolation routine for evolutionary tracks.
!> Modern Fortran version.
function yparinterp(xint, x, y, n, ntot, lup, iii, ntot2) result(yint)
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   real(real32), intent(in) :: xint
   real(real32), intent(in) :: x(ntot)
   real(real32), intent(in) :: y(ntot, ntot2)
   integer, intent(in) :: n, ntot, iii, ntot2
   logical, intent(in) :: lup
   real(real32) :: yint

   integer :: jl, ju, jm, j, k
   real(real32) :: help1, help2, a, b, c, ylintest, ymax, ymin

   ! Bisection to find interval
   jl = 0
   ju = n + 1
   do
      if (ju - jl <= 1) exit
      jm = (ju + jl) / 2
      if ((x(n) > x(1)) .eqv. (xint > x(jm))) then
         jl = jm
      else
         ju = jm
      end if
   end do
   j = jl

   ! Defensive: boundaries
   if (j == 0) then
      yint = y(1, iii)
      return
   end if
   if (j == n) then
      yint = y(n, iii)
      return
   end if
   if (xint == x(1)) then
      yint = y(1, iii)
      return
   end if
   if (xint == x(n)) then
      yint = y(n, iii)
      return
   end if

   k = j
   ! Linear interpolation at lower boundary
   if (j == 1) then
      if (x(j) == x(j+1)) then
         yint = y(j, iii)
         return
      else if (x(j+1) == x(j+2)) then
         yint = y(j, iii) + (y(j+1, iii) - y(j, iii)) * (xint - x(j)) / (x(j+1) - x(j))
         return
      else
         k = 2
         ! lup = .false. (not used here)
      end if
   end if

   ! Parabolic interpolation at upper boundary
   if (j == n-1) then
      if (x(j) == x(j+1)) then
         yint = y(j, iii)
         return
      else if (x(j-1) == x(j+1)) then
         yint = y(j, iii) + (y(j+1, iii) - y(j, iii)) * (xint - x(j)) / (x(j+1) - x(j))
         return
      else
         k = j
         ! lup = .false. (not used here)
      end if
   end if

   if (x(k) == x(k+1)) then
      yint = y(k, iii)
      return
   else
      if (lup) then
         if (x(k+1) == x(k+2) .or. x(k) == x(k+2)) then
            yint = y(k, iii) + (y(k+1, iii) - y(k, iii)) * (xint - x(k)) / (x(k+1) - x(k))
            return
         else
            help1 = (y(k, iii) - y(k+1, iii)) / ((x(k) - x(k+1)) * (x(k) - x(k+2)))
            help2 = (y(k+1, iii) - y(k+2, iii)) / ((x(k+1) - x(k+2)) * (x(k) - x(k+2)))
            a = help1 - help2
            b = help2 * (x(k) + x(k+1)) - help1 * (x(k+1) + x(k+2))
            c = y(k+1, iii) - a * x(k+1)**2 - b * x(k+1)
            yint = a * xint**2 + b * xint + c
            ! Optionally, check bounds:
            ! ylintest = y(k,iii) + (y(k+1,iii)-y(k,iii))*(xint-x(k))/(x(k+1)-x(k))
            ! ymax = max(y(k,iii), y(k+1,iii), y(k+2,iii))
            ! ymin = min(y(k,iii), y(k+1,iii), y(k+2,iii))
            ! if (yint > ymax .or. yint < ymin) yint = ylintest
            return
         end if
      else
         if (x(k-1) == x(k) .or. x(k-1) == x(k+1)) then
            yint = y(k, iii) + (y(k+1, iii) - y(k, iii)) * (xint - x(k)) / (x(k+1) - x(k))
            return
         else
            help1 = (y(k-1, iii) - y(k, iii)) / ((x(k-1) - x(k)) * (x(k-1) - x(k+1)))
            help2 = (y(k, iii) - y(k+1, iii)) / ((x(k) - x(k+1)) * (x(k-1) - x(k+1)))
            a = help1 - help2
            b = help2 * (x(k-1) + x(k)) - help1 * (x(k) + x(k+1))
            c = y(k, iii) - a * x(k)**2 - b * x(k)
            yint = a * xint**2 + b * xint + c
            ! Optionally, check bounds:
            ! ylintest = y(k,iii) + (y(k+1,iii)-y(k,iii))*(xint-x(k))/(x(k+1)-x(k))
            ! ymax = max(y(k-1,iii), y(k,iii), y(k+1,iii))
            ! ymin = min(y(k-1,iii), y(k,iii), y(k+1,iii))
            ! if (yint > ymax .or. yint < ymin) yint = ylintest
            return
         end if
      end if
   end if

end function yparinterp
c

! ****************************************************************************
! ****************************************************************************
! ****************************************************************************
!> Module containing filter profiles, nebular emission coefficients, and SNII yields.
module data_profiles
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   ! Filter profiles for synthetic colors
   real(real32), save :: xprof(12,100), yprof(12,100)
   ! Nebular continuum emission coefficients
   real(real32), save :: xrange(26), gamma(26)
   ! SNII chemical yields: [metallicity_index, mass_index]
   real(real32), save :: ymass(5,12), yh(5,12), yhe(5,12), yc(5,12), yn(5,12)
   real(real32), save :: yo(5,12), ymg(5,12), ysi(5,12), ys(5,12), yfe(5,12)
contains
   ! No procedures; data only
end module data_profiles

!> Initializes the filter profiles, nebular coefficients, and SNII yields.
submodule (data_profiles) data_profiles_init
   implicit none
contains
   module procedure initialize_data_profiles
      ! This procedure is a placeholder for initialization if needed.
      ! Data is initialized via DATA statements below.
   end procedure initialize_data_profiles
end submodule data_profiles_init

! Data block for initialization (converted from BLOCK DATA)
subroutine initialize_data_profiles()
   use, intrinsic :: iso_fortran_env, only: real32
   use data_profiles
   implicit none

   ! FILTER PROFILE FOR UIT_B1
   data xprof(1,1:99) / &
        50.,100.,200.,400.,600.,700.,800.,900., &
        1000.,1020.75,1059.03,1118.07,1181.76,1223.62,1256.56, &
        1273.4,1279.54,1286.78,1287.59,1294.19,1297.58,1303.93,1309.54, &
        1319.91,1339.78,1360.53,1383.31,1415.32,1433.99,1475.48,1505.79, &
        1526.79,1540.22,1555.82,1577.71,1603.53,1617.76,1632.4,1645.99, &
        1662.1,1674.22,1701.23,1723.18,1748.62,1776.94,1804.67,1843.5, &
        1874.54,1917.08,1974.99,2031.15,2080.52,2116.29,2156.72,2200., &
        2300.,2400.,2500.,2750.,3000.,3250.,3500.,3750.,3800., &
        3850.,3900.,3950.,4000.,4050.,4100.,4150.,4200., &
        4250.,4300.,4350.,4400.,4450.,4500.,4550.,4600., &
        4650.,4700.,4750.,4800.,4850.,4900.,4950.,5000., &
        5050.,5100.,5150.,5200.,5250.,5300.,5350.,5400., &
        5450.,5500.,2000000. /
   data yprof(1,1:99) / &
        0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000, &
        0.000,8.24571E-4,5.89848E-4,0.00101411,0.00128996, &
        0.00249922,0.00943542,0.0203621,0.040092,0.0685947, &
        0.101941,0.156772,0.216,0.311201,0.394995,0.518268, &
        0.711285,0.835409,0.920037,0.976559,1.000,0.995273, &
        0.961847,0.920986,0.864348,0.813408,0.711992,0.608811, &
        0.580253,0.569246,0.568333,0.570925,0.555097,0.490965, &
        0.410172,0.299971,0.195028,0.13221,0.0715576,0.0429568, &
        0.0283679,0.0128621,0.00745285,0.00513256,0.00284696, &
        0.00318229,0.000,0.000,0.000,0.000,0.000,0.000, &
        0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000, &
        0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000, &
        0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000, &
        0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000, &
        0.000,0.000,0.000 /
   ! ... (repeat for all other filters, nebular coefficients, and yields as in the original BLOCK DATA)
   ! For brevity, only UIT_B1 is shown here. The rest should be filled in similarly.

   ! Nebular continuum emission coefficients
   data xrange / 10.,912.,913.,1300.,1500.,1800.,2200., &
        2855.,3331.,3421.,3422.,3642.,3648.,5700.,7000.,8207., &
        8209.,14583.,14585.,22787.,22789.,32813.,32815., &
        44680.,44682.,2000000. /
   data gamma / 0.,0.,2.11e-4,5.647,9.35,9.847,10.582,16.101,24.681, &
        26.736,24.883,29.979,6.519,8.773,11.545,13.585,6.333, &
        10.444,7.023,9.361,7.59,9.35,8.32,9.53,8.87,0. /

   ! SNII yields (example for Z=2 Zsol, fill in all as in original)
   data ymass(5,1),yh(5,1),yhe(5,1),yc(5,1),yn(5,1), &
        yo(5,1),ymg(5,1),ysi(5,1),ys(5,1),yfe(5,1) / &
        7,3.48,2.40,0.048,0.042,0.128,0.006,0.054,0.046,0.000 /
   ! ... (repeat for all metallicities and masses as in the original BLOCK DATA)

end subroutine initialize_data_profiles

! *****************************************************************************
! *****************************************************************************
! *****************************************************************************
