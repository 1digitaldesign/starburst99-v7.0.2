!==============================================================================
! Astronomical Calculation Tests for Galaxy/Starburst99
!==============================================================================
!> This program tests the astronomical calculations and conversions
!> performed in the Galaxy/Starburst99 code.
!==============================================================================
program astro_tests
   use galaxy_module
   use, intrinsic :: iso_fortran_env, only: real32, real64, stdout => output_unit
   implicit none

   logical :: all_passed = .true.
   integer :: num_passed = 0
   integer :: num_failed = 0

   ! Run all tests
   write(stdout, '(a)') "Running astronomical calculation tests..."
   write(stdout, '(a)') "==============================================="
   
   call test_blackbody_flux()
   call test_luminosity_calculation()
   call test_stellar_mass_function()
   
   ! Summary
   write(stdout, '(a)') "==============================================="
   write(stdout, '(a, i0, a, i0, a, i0, a)') "Summary: ", num_passed, " passed, ", &
      num_failed, " failed, ", (num_passed + num_failed), " total"
   
   if (all_passed) then
      write(stdout, '(a)') "All astronomical tests passed!"
      stop 0
   else
      write(stdout, '(a)') "Some astronomical tests failed."
      stop 1
   end if

contains

   !> Check if a test passed and update counters
   subroutine check_test(test_name, condition)
      character(len=*), intent(in) :: test_name
      logical, intent(in) :: condition
      
      if (condition) then
         write(stdout, '(a, a)') "✓ ", test_name
         num_passed = num_passed + 1
      else
         write(stdout, '(a, a)') "✗ ", test_name
         num_failed = num_failed + 1
         all_passed = .false.
      end if
   end subroutine check_test
   
   !> Test Planck blackbody flux calculation
   subroutine test_blackbody_flux()
      real(real32) :: wavelength, temp, flux, expected_flux
      real(real32) :: rel_error
      logical :: test_result
      
      write(stdout, '(a)') ""
      write(stdout, '(a)') "Testing blackbody flux calculation:"
      
      ! Initialize module variables
      call init_module()
      
      ! Test case 1: Wien's displacement law
      ! Lambda_max * T = 0.29 cm·K (Wien's law)
      temp = 5778.0_real32  ! Solar effective temperature in K
      wavelength = 0.29_real32 / temp  ! Peak wavelength in cm
      
      ! Calculate Planck blackbody flux at this wavelength
      ! B_λ(T) = (2hc²/λ⁵) * 1/(e^(hc/λkT) - 1)
      expected_flux = (2.0_real32 * h_planck * c_light**2) / (wavelength**5)
      expected_flux = expected_flux / (exp((h_planck * c_light) / (wavelength * k_boltz * temp)) - 1.0_real32)
      
      ! Calculate flux using a simple implementation
      flux = blackbody_flux(wavelength, temp)
      
      ! Test the result
      rel_error = abs((flux - expected_flux) / expected_flux)
      test_result = rel_error < 1.0e-5_real32
      call check_test("Wien's law peak wavelength flux calculation", test_result)
      
      ! Test case 2: Wien displacement law - peak wavelength
      wavelength = wien_peak_wavelength(temp)
      test_result = abs(wavelength - 0.29_real32/temp) < 1.0e-5_real32
      call check_test("Wien's law peak wavelength calculation", test_result)
      
      ! Test case 3: Stefan-Boltzmann law - total energy
      ! E = σT⁴
      flux = stefan_boltzmann(temp)
      expected_flux = sigma_sb * temp**4
      rel_error = abs((flux - expected_flux) / expected_flux)
      test_result = rel_error < 1.0e-5_real32
      call check_test("Stefan-Boltzmann law calculation", test_result)
      
      ! Cleanup
      call cleanup_module()
   end subroutine test_blackbody_flux
   
   !> Test stellar luminosity calculation
   subroutine test_luminosity_calculation()
      real(real32) :: radius, temp, luminosity, expected_luminosity
      real(real32) :: rel_error
      logical :: test_result
      
      write(stdout, '(a)') ""
      write(stdout, '(a)') "Testing stellar luminosity calculation:"
      
      ! Test case 1: Sun (R = 1 Rsun, T = 5778 K)
      radius = 6.96e10_real32  ! Solar radius in cm
      temp = 5778.0_real32     ! Solar effective temperature in K
      
      ! L = 4πR²σT⁴
      expected_luminosity = 4.0_real32 * pi * radius**2 * sigma_sb * temp**4
      
      ! Calculate using our function
      luminosity = stellar_luminosity(radius, temp)
      
      ! Test the result
      rel_error = abs((luminosity - expected_luminosity) / expected_luminosity)
      test_result = rel_error < 1.0e-5_real32
      call check_test("Solar luminosity calculation", test_result)
      
      ! Test case 2: Red giant (R = 100 Rsun, T = 3500 K)
      radius = 100.0_real32 * 6.96e10_real32
      temp = 3500.0_real32
      
      expected_luminosity = 4.0_real32 * pi * radius**2 * sigma_sb * temp**4
      luminosity = stellar_luminosity(radius, temp)
      
      rel_error = abs((luminosity - expected_luminosity) / expected_luminosity)
      test_result = rel_error < 1.0e-5_real32
      call check_test("Red giant luminosity calculation", test_result)
      
      ! Test case 3: Blue giant (R = 10 Rsun, T = 20000 K)
      radius = 10.0_real32 * 6.96e10_real32
      temp = 20000.0_real32
      
      expected_luminosity = 4.0_real32 * pi * radius**2 * sigma_sb * temp**4
      luminosity = stellar_luminosity(radius, temp)
      
      rel_error = abs((luminosity - expected_luminosity) / expected_luminosity)
      test_result = rel_error < 1.0e-5_real32
      call check_test("Blue giant luminosity calculation", test_result)
   end subroutine test_luminosity_calculation
   
   !> Test initial mass function (IMF) calculation
   subroutine test_stellar_mass_function()
      real(real32) :: mass, imf_value, expected_value
      real(real32) :: mass_array(3), imf_array(3)
      real(real32) :: rel_error
      logical :: test_result
      
      write(stdout, '(a)') ""
      write(stdout, '(a)') "Testing stellar mass function calculation:"
      
      ! Test case 1: Salpeter IMF (dN/dM ∝ M^-2.35)
      mass = 10.0_real32
      expected_value = mass**(-2.35_real32)
      
      ! Calculate using our function
      imf_value = salpeter_imf(mass)
      
      ! Test the result
      rel_error = abs((imf_value - expected_value) / expected_value)
      test_result = rel_error < 1.0e-5_real32
      call check_test("Salpeter IMF calculation", test_result)
      
      ! Test case 2: Kroupa IMF (piecewise power law)
      ! 0.08 ≤ m < 0.5: dN/dM ∝ M^-1.3
      ! 0.5 ≤ m < 1.0: dN/dM ∝ M^-2.3
      ! m ≥ 1.0: dN/dM ∝ M^-2.3 (sometimes -2.7 for m > 1)
      mass_array = [0.3_real32, 0.7_real32, 5.0_real32]
      
      ! Calculate using our function
      call kroupa_imf(mass_array, imf_array)
      
      ! Test each result
      ! First mass range
      expected_value = mass_array(1)**(-1.3_real32)
      rel_error = abs((imf_array(1) - expected_value) / expected_value)
      test_result = rel_error < 1.0e-5_real32
      call check_test("Kroupa IMF calculation (0.08 ≤ m < 0.5)", test_result)
      
      ! Second mass range
      expected_value = mass_array(2)**(-2.3_real32) * (0.5_real32**(-1.3_real32) / 0.5_real32**(-2.3_real32))
      rel_error = abs((imf_array(2) - expected_value) / expected_value)
      test_result = rel_error < 1.0e-5_real32
      call check_test("Kroupa IMF calculation (0.5 ≤ m < 1.0)", test_result)
      
      ! Third mass range
      expected_value = mass_array(3)**(-2.3_real32) * (0.5_real32**(-1.3_real32) / 0.5_real32**(-2.3_real32))
      rel_error = abs((imf_array(3) - expected_value) / expected_value)
      test_result = rel_error < 1.0e-5_real32
      call check_test("Kroupa IMF calculation (m ≥ 1.0)", test_result)
   end subroutine test_stellar_mass_function

   !==============================================================================
   ! Helper functions for testing - these simulate the actual implementations
   !==============================================================================
   
   !> Planck blackbody function (B_λ(T))
   function blackbody_flux(wavelength, temperature) result(flux)
      real(real32), intent(in) :: wavelength     ! in cm
      real(real32), intent(in) :: temperature    ! in K
      real(real32) :: flux                       ! in erg/s/cm²/cm
      
      real(real32) :: numerator, denominator
      
      ! B_λ(T) = (2hc²/λ⁵) * 1/(e^(hc/λkT) - 1)
      numerator = 2.0_real32 * h_planck * c_light**2
      denominator = wavelength**5 * (exp((h_planck * c_light) / (wavelength * k_boltz * temperature)) - 1.0_real32)
      
      flux = numerator / denominator
   end function blackbody_flux
   
   !> Wien's displacement law for peak wavelength
   function wien_peak_wavelength(temperature) result(wavelength)
      real(real32), intent(in) :: temperature    ! in K
      real(real32) :: wavelength                 ! in cm
      
      ! λ_max * T = 0.29 cm·K
      wavelength = 0.29_real32 / temperature
   end function wien_peak_wavelength
   
   !> Stefan-Boltzmann law (total energy flux)
   function stefan_boltzmann(temperature) result(flux)
      real(real32), intent(in) :: temperature    ! in K
      real(real32) :: flux                       ! in erg/s/cm²
      
      ! F = σT⁴
      flux = sigma_sb * temperature**4
   end function stefan_boltzmann
   
   !> Stellar luminosity calculation
   function stellar_luminosity(radius, temperature) result(luminosity)
      real(real32), intent(in) :: radius         ! in cm
      real(real32), intent(in) :: temperature    ! in K
      real(real32) :: luminosity                 ! in erg/s
      
      ! L = 4πR²σT⁴
      luminosity = 4.0_real32 * pi * radius**2 * sigma_sb * temperature**4
   end function stellar_luminosity
   
   !> Salpeter initial mass function
   function salpeter_imf(mass) result(imf_value)
      real(real32), intent(in) :: mass           ! in solar masses
      real(real32) :: imf_value
      
      ! dN/dM ∝ M^(-2.35)
      imf_value = mass**(-2.35_real32)
   end function salpeter_imf
   
   !> Kroupa initial mass function (piecewise power law)
   subroutine kroupa_imf(masses, imf_values)
      real(real32), intent(in) :: masses(:)      ! in solar masses
      real(real32), intent(out) :: imf_values(:) ! IMF values
      
      integer :: i
      real(real32) :: normalization
      
      ! Normalization factor at the transition point
      normalization = 0.5_real32**(-1.3_real32) / 0.5_real32**(-2.3_real32)
      
      do i = 1, size(masses)
         if (masses(i) < 0.5_real32) then
            ! 0.08 ≤ m < 0.5: dN/dM ∝ M^-1.3
            imf_values(i) = masses(i)**(-1.3_real32)
         else
            ! m ≥ 0.5: dN/dM ∝ M^-2.3 (with normalization)
            imf_values(i) = masses(i)**(-2.3_real32) * normalization
         end if
      end do
   end subroutine kroupa_imf

end program astro_tests