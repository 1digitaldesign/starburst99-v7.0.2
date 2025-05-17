!> Galaxy code module for shared data and constants
!>
!> This module provides shared variables, constants, and utility functions
!> for the galaxy stellar population synthesis code, using modern Fortran 2018
!> features to replace outdated COMMON blocks and implicit global variables.
module galaxy_module
   use, intrinsic :: iso_fortran_env, only: real32, real64, int32, int64, &
                                          stdout => output_unit, &
                                          stderr => error_unit
   implicit none
   private  ! Default everything to private
   
   ! Explicitly export public entities
   public :: init_module, cleanup_module, open_file, linear_interp
   public :: exp10, integer_to_string
   
   ! Export all necessary variables and types
   public :: pi, solar_mass, solar_lum, year_in_sec, k_boltz, h_planck, c_light
   public :: nmaxint, nmaxint1, np, np1, npgrid
   public :: un_input, un_output, un_spectrum, un_quanta, un_snr, un_power, un_sptyp
   public :: un_yield, un_uvline, un_color, un_atm, un_debug, un_wrline, un_hires, un_fuse
   public :: model_name, isf, toma, sfr, ninterv, xponent, xmaslim, sncut, bhcut
   public :: iz, iwind, time1, jtime, tbiv, itbiv, tmax, jmg, lmin, lmax, tdel
   public :: iatmos, ilib, iline, ivt, irsg
   public :: io1, io2, io3, io4, io5, io6, io7, io8, io9, io10, io11, io12, io13, io14, io15
   public :: tvar, tinter, tiempo1, tstep, upma, doma
   public :: iwrt, iwrscale, xmwr
   public :: critma, critma_new, critup, critup_new, critma1, critma_new1, critup1, critup_new1
   public :: cmass, dens, wavel, spectra
   public :: track_data, tracks, model_parameters
   public :: wind_power, sn_rates, sp_type_counts, element_yields, uv_lines, fuv_lines, hires_lines
   
   !---------------------------------------------------------------------------
   ! Constants 
   !---------------------------------------------------------------------------
   ! Mathematical constants
   real(real32), parameter :: pi = 3.14159265358979323846_real32
   
   ! Astronomical constants
   real(real32), parameter :: solar_mass = 1.989e33_real32     ! g
   real(real32), parameter :: solar_lum = 3.826e33_real32      ! erg/s
   real(real32), parameter :: year_in_sec = 3.1557e7_real32    ! seconds
   real(real32), parameter :: parsec = 3.0857e18_real32        ! cm
   real(real32), parameter :: lsun_mw = 4.0e10_real32          ! L_sun

   ! Atomic constants
   real(real32), parameter :: k_boltz = 1.380649e-16_real32    ! erg/K
   real(real32), parameter :: h_planck = 6.62607015e-27_real32 ! erg-s
   real(real32), parameter :: c_light = 2.99792458e10_real32   ! cm/s
   real(real32), parameter :: sigma_sb = 5.670374419e-5_real32 ! erg cm^-2 s^-1 K^-4

   !---------------------------------------------------------------------------
   ! Module-level parameters 
   !---------------------------------------------------------------------------
   integer, parameter :: nmaxint = 10     ! Maximum number of IMF intervals
   integer, parameter :: nmaxint1 = 11    ! Maximum number of IMF mass boundaries
   integer, parameter :: np = 860         ! Size of standard wavelength grid
   integer, parameter :: np1 = 1415       ! Size of extended wavelength grid
   integer, parameter :: npgrid = 3000    ! Size of mass grid
   
   !---------------------------------------------------------------------------
   ! File unit numbers
   !---------------------------------------------------------------------------
   ! Instead of using hard-coded unit numbers, we define them here
   integer, parameter :: un_input = 10     ! Input parameter file
   integer, parameter :: un_output = 11    ! Main output file
   integer, parameter :: un_spectrum = 12  ! Spectrum output
   integer, parameter :: un_quanta = 13    ! Ionizing photon output
   integer, parameter :: un_snr = 14       ! Supernova rate output
   integer, parameter :: un_power = 15     ! Wind power output
   integer, parameter :: un_sptyp = 16     ! Spectral types output
   integer, parameter :: un_yield = 17     ! Chemical yields output
   integer, parameter :: un_uvline = 18    ! UV line spectrum output
   integer, parameter :: un_color = 19     ! Color index output
   integer, parameter :: un_atm = 20       ! Atmosphere model data
   integer, parameter :: un_debug = 21     ! Debug output
   integer, parameter :: un_wrline = 22    ! WR line output
   integer, parameter :: un_hires = 23     ! High-resolution output
   integer, parameter :: un_fuse = 24      ! FUSE output

   !---------------------------------------------------------------------------
   ! Derived types for main data structures
   !---------------------------------------------------------------------------
   !> Type containing all model configuration parameters
   type :: model_parameters
      ! General model parameters
      character(len=:), allocatable :: name    ! Model designation
      integer :: sf_mode                        ! Star formation mode: >0 continuous, <=0 fixed mass
      real(real32) :: total_mass               ! Total stellar mass (10^6 solar masses)
      real(real32) :: sf_rate                  ! Star formation rate (solar masses per year)
      
      ! IMF parameters
      integer :: num_intervals                 ! Number of intervals for the IMF
      real(real32), allocatable :: exponents(:)  ! IMF exponents
      real(real32), allocatable :: mass_limits(:) ! Mass boundaries for IMF (solar masses)
      real(real32) :: sn_cutoff                ! Supernova cut-off mass (solar masses)
      real(real32) :: bh_cutoff                ! Black hole cut-off mass (solar masses)
      
      ! Model selection parameters
      integer :: metallicity_id                ! Metallicity + tracks identifier
      integer :: wind_model                    ! Wind model (0: evolution, 1: emp., 2: theor., 3: Elson)
      
      ! Time step parameters
      real(real32) :: initial_time             ! Initial time (years)
      integer :: time_scale                    ! 0=linear, 1=logarithmic
      real(real32) :: time_step                ! Time step (years) if time_scale=0
      integer :: num_steps                     ! Number of steps if time_scale=1
      real(real32) :: max_time                 ! Last grid point (years)
      
      ! Atmosphere and library options
      integer :: atmosphere_model              ! Atmosphere model options
      integer :: hires_metallicity             ! Metallicity of high-resolution models
      integer :: uv_library                    ! Library for UV line spectrum (1=solar, 2=LMC/SMC)
      
      ! Output control
      logical, allocatable :: outputs(:)       ! Output flags for various outputs
   end type model_parameters
   
   !> Type for stellar evolutionary track data
   type :: track_data
      ! Basic track information
      integer :: num_masses                    ! Number of initial masses
      integer :: num_points                    ! Number of time points per track
      real(real32) :: metallicity              ! Z value for this track set
      character(len=:), allocatable :: source  ! Source of track data (Geneva, Padova, etc.)
      
      ! Main track arrays
      real(real32), allocatable :: init_mass(:)       ! Initial stellar mass (solar masses)
      real(real32), allocatable :: log_init_mass(:)   ! Log of initial mass
      
      ! Per-track, per-time point data
      real(real32), allocatable :: age(:,:)           ! Stellar age (years)
      real(real32), allocatable :: log_age(:,:)       ! Log of age
      real(real32), allocatable :: mass(:,:)          ! Current mass (solar masses)
      real(real32), allocatable :: log_mass(:,:)      ! Log of current mass
      real(real32), allocatable :: log_lum(:,:)       ! Log luminosity (solar units)
      real(real32), allocatable :: log_teff(:,:)      ! Log effective temperature (K)
      real(real32), allocatable :: mdot(:,:)          ! Mass loss rate (solar masses/yr)
      
      ! Surface abundances
      real(real32), allocatable :: h_frac(:,:)        ! Surface hydrogen mass fraction
      real(real32), allocatable :: he_frac(:,:)       ! Surface helium mass fraction
      real(real32), allocatable :: c_frac(:,:)        ! Surface carbon mass fraction
      real(real32), allocatable :: n_frac(:,:)        ! Surface nitrogen mass fraction
      real(real32), allocatable :: o_frac(:,:)        ! Surface oxygen mass fraction
      
   contains
      procedure :: init => initialize_track_data
      procedure :: cleanup => cleanup_track_data
      procedure :: get_mass_index => find_mass_index
      procedure :: interpolate_in_time => interp_track_time
   end type track_data
   
   !---------------------------------------------------------------------------
   ! Model configuration variables  
   !---------------------------------------------------------------------------
   ! General model parameters
   ! These will be migrated to model_parameters type in future refactoring
   character(len=20) :: model_name          ! Model designation
   integer :: isf                           ! Star formation mode: >0 continuous, <=0 fixed mass
   real(real32) :: toma                     ! Total stellar mass (10^6 solar masses)
   real(real32) :: sfr                      ! Star formation rate (solar masses per year)
   
   ! IMF parameters
   integer :: ninterv                       ! Number of intervals for the IMF
   real(real32) :: xponent(nmaxint)         ! IMF exponents
   real(real32) :: xmaslim(nmaxint1)        ! Mass boundaries for IMF (solar masses)
   real(real32) :: sncut                    ! Supernova cut-off mass (solar masses)
   real(real32) :: bhcut                    ! Black hole cut-off mass (solar masses)
   
   ! Metallicity, tracks, and wind model
   integer :: iz                            ! Metallicity + tracks identifier
   integer :: iwind                         ! Wind model (0: evolution, 1: emp., 2: theor., 3: Elson)
   
   ! Time step parameters
   real(real32) :: time1                    ! Initial time (10^6 years)
   integer :: jtime                         ! Time scale: 0=linear, 1=logarithmic
   real(real32) :: tbiv                     ! Time step (10^6 years) if jtime=0
   integer :: itbiv                         ! Number of steps if jtime=1
   real(real32) :: tmax                     ! Last grid point (10^6 years)
   
   ! Grid parameters
   ! Synthesis method selection:
   integer :: jmg                           ! 0: small mass grid
                                            ! 1: large mass grid 
                                            ! 2: isochrone on large grid
                                            ! 3: full isochrone
   integer :: lmin, lmax                    ! Min and max indices of evolutionary tracks
   real(real32) :: tdel                     ! Time step for printing spectra (10^6 years)
   
   ! Atmosphere and library options
   integer :: iatmos                        ! Atmosphere model options
   integer :: ilib                          ! Metallicity of high-resolution models
   integer :: iline                         ! Library for UV line spectrum (1=solar, 2=LMC/SMC)
   integer :: ivt, irsg                     ! RSG feature parameters
   
   ! Output options
   integer :: io1, io2, io3, io4, io5       ! Output flags for various outputs
   integer :: io6, io7, io8, io9, io10      ! More output flags
   integer :: io11, io12, io13, io14, io15  ! More output flags
   
   ! Derived time variables
   real(real32) :: tvar                     ! Current time step size
   real(real32) :: tinter                   ! Number of time steps
   real(real32) :: tiempo1                  ! Current time (log scale)
   real(real32) :: tstep                    ! Current time step size adjusted for log/linear scale
   real(real32) :: upma, doma               ! Upper and lower mass limits

   ! WR related parameters
   integer :: iwrt                          ! WR temperature adjustment method
   integer :: iwrscale                      ! WR scaling parameter
   real(real32) :: xmwr                     ! Minimum mass for WR stars
   
   ! SN/nucleosynthesis mass limits
   real(real32) :: critma                   ! Critical mass for supernovae
   real(real32) :: critma_new               ! Updated critical mass
   real(real32) :: critup                   ! Upper critical mass
   real(real32) :: critup_new               ! Updated upper critical mass
   real(real32) :: critma1                  ! Secondary critical mass
   real(real32) :: critma_new1              ! Updated secondary critical mass
   real(real32) :: critup1                  ! Secondary upper critical mass
   real(real32) :: critup_new1              ! Updated secondary upper critical mass

   ! Additional variables used in arrays 
   ! (could be expanded based on the full code structure)
   real(real32), allocatable :: cmass(:)    ! Grid of stellar masses
   real(real32), allocatable :: dens(:)     ! Number density of stars per mass bin
   real(real32), allocatable :: wavel(:)    ! Wavelength grid
   real(real32), allocatable :: spectra(:,:)! Spectral data arrays
   
   ! Arrays for physical outputs
   real(real32), allocatable :: wind_power(:)       ! Wind power per mass bin
   real(real32), allocatable :: sn_rates(:)         ! Supernova rates per mass bin
   real(real32), allocatable :: sp_type_counts(:)   ! Counts of different spectral types
   real(real32), allocatable :: element_yields(:)   ! Chemical yields for different elements
   real(real32), allocatable :: uv_lines(:)         ! UV line strengths
   real(real32), allocatable :: fuv_lines(:)        ! FUV line strengths
   real(real32), allocatable :: hires_lines(:)      ! High-resolution optical line strengths

   ! Track data for multiple metallicities
   type(track_data), allocatable :: tracks(:)
   
   ! Enhanced error message variable
   character(len=:), allocatable :: error_message
   
   !---------------------------------------------------------------------------
   ! Abstract interfaces for procedures
   !---------------------------------------------------------------------------
   abstract interface
      function interp_func(x, y, xval) result(yval)
         import :: real32
         real(real32), intent(in) :: x(:), y(:)
         real(real32), intent(in) :: xval
         real(real32) :: yval
      end function interp_func
   end interface
   
   !---------------------------------------------------------------------------
   ! Interface blocks for utility functions
   !---------------------------------------------------------------------------
   interface
      module subroutine error_handler(routine, message, fatal)
         character(len=*), intent(in) :: routine
         character(len=*), intent(in) :: message
         logical, intent(in), optional :: fatal
      end subroutine error_handler
      
      pure module function flin(x1, y1, x2, y2, x) result(y)
         real(real32), intent(in) :: x1, y1, x2, y2, x
         real(real32) :: y
      end function flin
   end interface

contains

   !---------------------------------------------------------------------------
   ! Type-bound procedures for track_data
   !---------------------------------------------------------------------------
   
   !> Initialize track_data structure
   subroutine initialize_track_data(this, num_masses, num_points, metallicity, source)
      class(track_data), intent(inout) :: this
      integer, intent(in) :: num_masses, num_points
      real(real32), intent(in) :: metallicity
      character(len=*), intent(in), optional :: source
      
      ! Set basic information
      this%num_masses = num_masses
      this%num_points = num_points
      this%metallicity = metallicity
      
      ! Set source if provided
      if (present(source)) then
         this%source = source
      else
         this%source = "Unknown"
      end if
      
      ! Allocate arrays
      if (allocated(this%init_mass)) call this%cleanup()
      
      allocate(this%init_mass(num_masses))
      allocate(this%log_init_mass(num_masses))
      
      allocate(this%age(num_masses, num_points))
      allocate(this%log_age(num_masses, num_points))
      allocate(this%mass(num_masses, num_points))
      allocate(this%log_mass(num_masses, num_points))
      allocate(this%log_lum(num_masses, num_points))
      allocate(this%log_teff(num_masses, num_points))
      allocate(this%mdot(num_masses, num_points))
      
      allocate(this%h_frac(num_masses, num_points))
      allocate(this%he_frac(num_masses, num_points))
      allocate(this%c_frac(num_masses, num_points))
      allocate(this%n_frac(num_masses, num_points))
      allocate(this%o_frac(num_masses, num_points))
      
      ! Initialize arrays
      this%init_mass = 0.0_real32
      this%log_init_mass = 0.0_real32
      this%age = 0.0_real32
      this%log_age = 0.0_real32
      this%mass = 0.0_real32
      this%log_mass = 0.0_real32
      this%log_lum = 0.0_real32
      this%log_teff = 0.0_real32
      this%mdot = 0.0_real32
      this%h_frac = 0.0_real32
      this%he_frac = 0.0_real32
      this%c_frac = 0.0_real32
      this%n_frac = 0.0_real32
      this%o_frac = 0.0_real32
   end subroutine initialize_track_data
   
   !> Clean up track_data allocatable arrays
   subroutine cleanup_track_data(this)
      class(track_data), intent(inout) :: this
      
      if (allocated(this%init_mass)) deallocate(this%init_mass)
      if (allocated(this%log_init_mass)) deallocate(this%log_init_mass)
      if (allocated(this%age)) deallocate(this%age)
      if (allocated(this%log_age)) deallocate(this%log_age)
      if (allocated(this%mass)) deallocate(this%mass)
      if (allocated(this%log_mass)) deallocate(this%log_mass)
      if (allocated(this%log_lum)) deallocate(this%log_lum)
      if (allocated(this%log_teff)) deallocate(this%log_teff)
      if (allocated(this%mdot)) deallocate(this%mdot)
      if (allocated(this%h_frac)) deallocate(this%h_frac)
      if (allocated(this%he_frac)) deallocate(this%he_frac)
      if (allocated(this%c_frac)) deallocate(this%c_frac)
      if (allocated(this%n_frac)) deallocate(this%n_frac)
      if (allocated(this%o_frac)) deallocate(this%o_frac)
      if (allocated(this%source)) deallocate(this%source)
   end subroutine cleanup_track_data
   
   !> Find the index of the closest mass track
   pure function find_mass_index(this, mass) result(idx)
      class(track_data), intent(in) :: this
      real(real32), intent(in) :: mass
      integer :: idx
      integer :: i
      real(real32) :: log_mass, min_diff, diff
      
      log_mass = log10(mass)
      idx = 1
      min_diff = abs(this%log_init_mass(1) - log_mass)
      
      do i = 2, this%num_masses
         diff = abs(this%log_init_mass(i) - log_mass)
         if (diff < min_diff) then
            min_diff = diff
            idx = i
         end if
      end do
   end function find_mass_index
   
   !> Interpolate track properties for a given time
   !> Interpolate track properties for a given time
   !> 
   !> @param[in] this The track_data instance
   !> @param[in] mass_idx The index of the mass track to use
   !> @param[in] time The time point to interpolate at (in years)
   !> @return Interpolated stellar properties at the given time
   function interp_track_time(this, mass_idx, time) result(props)
      class(track_data), intent(in) :: this
      integer, intent(in) :: mass_idx
      real(real32), intent(in) :: time
      
      ! Return structure with interpolated properties
      type :: track_properties
         real(real32) :: mass
         real(real32) :: log_lum
         real(real32) :: log_teff
         real(real32) :: mdot
         real(real32) :: h_frac
         real(real32) :: he_frac
         real(real32) :: c_frac
         real(real32) :: n_frac
         real(real32) :: o_frac
      end type track_properties
      
      type(track_properties) :: props
      real(real32) :: log_time, t1, t2
      integer :: i, idx
      
      ! Find time indices
      log_time = log10(time)
      idx = 1
      
      ! Find index where time is between time points
      do i = 1, this%num_points-1
         if (this%log_age(mass_idx, i) <= log_time .and. &
             this%log_age(mass_idx, i+1) > log_time) then
            idx = i
            exit
         end if
      end do
      
      ! Special case for times outside track range
      if (log_time < this%log_age(mass_idx, 1)) then
         idx = 1
      else if (log_time > this%log_age(mass_idx, this%num_points)) then
         idx = this%num_points - 1
      end if
      
      ! Linear interpolation parameter
      if (abs(this%log_age(mass_idx, idx+1) - this%log_age(mass_idx, idx)) < epsilon(1.0_real32)) then
         t1 = 0.0_real32
      else
         t1 = (log_time - this%log_age(mass_idx, idx)) / &
              (this%log_age(mass_idx, idx+1) - this%log_age(mass_idx, idx))
      end if
      t2 = 1.0_real32 - t1
      
      ! Interpolate all properties
      props%mass = 10.0_real32 ** (t2 * this%log_mass(mass_idx, idx) + &
                                  t1 * this%log_mass(mass_idx, idx+1))
      props%log_lum = t2 * this%log_lum(mass_idx, idx) + &
                      t1 * this%log_lum(mass_idx, idx+1)
      props%log_teff = t2 * this%log_teff(mass_idx, idx) + &
                       t1 * this%log_teff(mass_idx, idx+1)
      props%mdot = t2 * this%mdot(mass_idx, idx) + &
                   t1 * this%mdot(mass_idx, idx+1)
      props%h_frac = t2 * this%h_frac(mass_idx, idx) + &
                     t1 * this%h_frac(mass_idx, idx+1)
      props%he_frac = t2 * this%he_frac(mass_idx, idx) + &
                      t1 * this%he_frac(mass_idx, idx+1)
      props%c_frac = t2 * this%c_frac(mass_idx, idx) + &
                     t1 * this%c_frac(mass_idx, idx+1)
      props%n_frac = t2 * this%n_frac(mass_idx, idx) + &
                     t1 * this%n_frac(mass_idx, idx+1)
      props%o_frac = t2 * this%o_frac(mass_idx, idx) + &
                     t1 * this%o_frac(mass_idx, idx+1)
   end function interp_track_time

   !---------------------------------------------------------------------------
   ! Utility functions
   !---------------------------------------------------------------------------
   
   !> Base-10 exponential function (10^x)
   !> This is a convenience function for 10^x calculation
   pure elemental function exp10(x) result(y)
      real(real32), intent(in) :: x
      real(real32) :: y
      
      y = 10.0_real32 ** x
   end function exp10
   
   !---------------------------------------------------------------------------
   ! Utility subroutines
   !---------------------------------------------------------------------------
   
   !> Initialize the module variables
   subroutine init_module()
      ! Allocate arrays based on parameters
      if (allocated(cmass)) deallocate(cmass)
      if (allocated(dens)) deallocate(dens)
      if (allocated(wavel)) deallocate(wavel)
      if (allocated(spectra)) deallocate(spectra)
      
      allocate(cmass(npgrid), source=0.0_real32)
      allocate(dens(npgrid), source=0.0_real32)
      allocate(wavel(np), source=0.0_real32)
      allocate(spectra(np, 3), source=0.0_real32)  ! stellar, nebular, total
   end subroutine init_module
   
   !> Clean up allocatable arrays when program ends
   subroutine cleanup_module()
      if (allocated(cmass)) deallocate(cmass)
      if (allocated(dens)) deallocate(dens)
      if (allocated(wavel)) deallocate(wavel)
      if (allocated(spectra)) deallocate(spectra)
      
      ! Clean up physical output arrays
      if (allocated(wind_power)) deallocate(wind_power)
      if (allocated(sn_rates)) deallocate(sn_rates)
      if (allocated(sp_type_counts)) deallocate(sp_type_counts)
      if (allocated(element_yields)) deallocate(element_yields)
      if (allocated(uv_lines)) deallocate(uv_lines)
      if (allocated(fuv_lines)) deallocate(fuv_lines)
      if (allocated(hires_lines)) deallocate(hires_lines)
      
      if (allocated(tracks)) then
         call cleanup_tracks()
         deallocate(tracks)
      end if
      if (allocated(error_message)) deallocate(error_message)
   end subroutine cleanup_module
   
   !> Clean up track data
   subroutine cleanup_tracks()
      integer :: i
      
      if (allocated(tracks)) then
         do i = 1, size(tracks)
            call tracks(i)%cleanup()
         end do
      end if
   end subroutine cleanup_tracks
   
   !> Safely open a file with error handling
   !>
   !> @param[in] unit The Fortran logical unit number to use
   !> @param[in] filename The name of the file to open
   !> @param[in] status File status ('old', 'new', 'replace', etc.)
   !> @param[in] form Optional file format ('formatted' or 'unformatted')
   !> @param[out] iostat Optional I/O status code (if not present, errors are fatal)
   subroutine open_file(unit, filename, status, form, iostat)
      integer, intent(in) :: unit
      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: status
      character(len=*), intent(in), optional :: form
      integer, intent(out), optional :: iostat
      
      integer :: ios
      character(len=:), allocatable :: form_loc
      character(len=:), allocatable :: message
      character(len=200) :: errmsg
      
      ! Default to formatted
      form_loc = 'formatted'
      if (present(form)) form_loc = form
      
      ! Open with condition handling
      if (form_loc == 'formatted') then
         open(unit=unit, file=filename, status=status, iostat=ios, &
              iomsg=errmsg)
      else
         open(unit=unit, file=filename, status=status, form=form_loc, &
              iostat=ios, iomsg=errmsg)
      end if
      
      ! Handle result
      if (ios == 0) then
         ! Normal return - no error
         if (present(iostat)) iostat = 0
      else
         ! Error occurred
         if (present(iostat)) then
            iostat = ios
         else
            ! Format error message
            message = "Error opening file: " // trim(filename) // &
                     " (status=" // trim(status) // ", iostat=" // &
                     trim(integer_to_string(ios)) // ") - " // trim(errmsg)
            
            ! Report error through error handler
            call error_handler("open_file", message, fatal=.true.)
         end if
      end if
   end subroutine open_file
   
   !> Convert integer to string using internal write
   pure function integer_to_string(i) result(str)
      integer, intent(in) :: i
      character(len=:), allocatable :: str
      character(len=20) :: temp
      
      write(temp, '(I0)') i
      str = trim(temp)
   end function integer_to_string
   
   !> Linear interpolation function (wrapper for flin to maintain API compatibility)
   pure elemental function linear_interp(x1, y1, x2, y2, x) result(y)
      real(real32), intent(in) :: x1, y1, x2, y2, x
      real(real32) :: y
      
      ! Call the submodule implementation
      y = flin(x1, y1, x2, y2, x)
   end function linear_interp

end module galaxy_module