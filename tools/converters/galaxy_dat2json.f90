!==============================================================================
! DAT to JSON Converter Module
!==============================================================================
!> Module for converting GALAXY/Starburst99 data files to JSON format
module galaxy_dat2json
    use, intrinsic :: iso_fortran_env, only: real32, real64, int32, int64, &
                                           stdout => output_unit, &
                                           stderr => error_unit
    implicit none
    private  ! Default everything to private
    
    ! Public procedures
    public :: convert_irfeatures, convert_tracks, convert_modfiles
    
    ! Constants for file handling
    integer, parameter :: MAX_LINE_LEN = 1024
    integer, parameter :: MAX_COLUMNS = 100
    
contains

    !> Convert irfeatures.dat file to JSON
    subroutine convert_irfeatures(input_file, output_file)
        character(len=*), intent(in) :: input_file
        character(len=*), intent(in) :: output_file
        
        ! Local variables
        integer :: io_stat, unit_in, unit_out
        character(len=MAX_LINE_LEN) :: line
        character(len=:), allocatable :: errmsg
        real(real32), allocatable :: wavelengths(:), extinctions(:)
        real(real32), allocatable :: data_matrix(:,:)
        integer :: i, j, num_wavelengths, num_extinctions
        
        ! Open input file
        open(newunit=unit_in, file=input_file, status='old', action='read', &
             iostat=io_stat, iomsg=errmsg)
        if (io_stat /= 0) then
            write(stderr, '(a)') 'Error opening input file: ' // trim(errmsg)
            return
        end if
        
        ! Read first line containing wavelengths
        read(unit_in, *, iostat=io_stat) ! Skip first line
        if (io_stat /= 0) then
            write(stderr, '(a)') 'Error reading wavelengths: ' // trim(errmsg)
            close(unit_in)
            return
        end if
        
        ! Count number of wavelengths (columns in first row)
        rewind(unit_in)
        read(unit_in, '(a)') line
        ! Count the number of values by parsing the line
        num_wavelengths = 0
        do i = 1, len_trim(line)
            if (i == 1 .and. line(i:i) /= ' ') then
                num_wavelengths = num_wavelengths + 1
            else if (i > 1 .and. line(i:i) /= ' ' .and. line(i-1:i-1) == ' ') then
                num_wavelengths = num_wavelengths + 1
            end if
        end do
        
        allocate(wavelengths(num_wavelengths))
        read(line, *, iostat=io_stat) wavelengths
        
        ! Read second line containing extinction values
        read(unit_in, '(a)') line
        ! Count the number of values by parsing the line
        num_extinctions = 0
        do i = 1, len_trim(line)
            if (i == 1 .and. line(i:i) /= ' ') then
                num_extinctions = num_extinctions + 1
            else if (i > 1 .and. line(i:i) /= ' ' .and. line(i-1:i-1) == ' ') then
                num_extinctions = num_extinctions + 1
            end if
        end do
        
        allocate(extinctions(num_extinctions))
        read(line, *, iostat=io_stat) extinctions
        
        ! Allocate data matrix
        allocate(data_matrix(num_wavelengths, num_extinctions))
        
        ! Read data values - make sure we read the right number of lines
        do i = 1, num_wavelengths
            read(unit_in, '(a)', iostat=io_stat) line
            if (io_stat /= 0) then
                write(stderr, '(a,i0)') 'Error reading data line ', i
                exit
            end if
            
            ! Parse each number in the line
            read(line, *, iostat=io_stat) (data_matrix(i,j), j=1, num_extinctions)
            if (io_stat /= 0) then
                write(stderr, '(a,i0)') 'Error parsing values in line ', i
                exit
            end if
        end do
        
        close(unit_in)
        
        ! Write JSON output
        open(newunit=unit_out, file=output_file, status='replace', action='write', &
             iostat=io_stat, iomsg=errmsg)
        if (io_stat /= 0) then
            write(stderr, '(a)') 'Error opening output file: ' // trim(errmsg)
            return
        end if
        
        ! Write JSON header
        write(unit_out, '(a)') '{'
        
        ! Write wavelengths array
        write(unit_out, '(a)') '  "wavelengths": ['
        do i = 1, num_wavelengths
            if (i < num_wavelengths) then
                write(unit_out, '(f10.2, a)') wavelengths(i), ','
            else
                write(unit_out, '(f10.2)') wavelengths(i)
            end if
        end do
        write(unit_out, '(a)') '  ],'
        
        ! Write extinctions array
        write(unit_out, '(a)') '  "extinctions": ['
        do i = 1, num_extinctions
            if (i < num_extinctions) then
                write(unit_out, '(f10.2, a)') extinctions(i), ','
            else
                write(unit_out, '(f10.2)') extinctions(i)
            end if
        end do
        write(unit_out, '(a)') '  ],'
        
        ! Write data matrix
        write(unit_out, '(a)') '  "data": ['
        do i = 1, num_wavelengths
            write(unit_out, '(a)') '    ['
            do j = 1, num_extinctions
                if (j < num_extinctions) then
                    write(unit_out, '(f10.3, a)') data_matrix(i,j), ','
                else
                    write(unit_out, '(f10.3)') data_matrix(i,j)
                end if
            end do
            
            if (i < num_wavelengths) then
                write(unit_out, '(a)') '    ],'
            else
                write(unit_out, '(a)') '    ]'
            end if
        end do
        write(unit_out, '(a)') '  ]'
        
        ! Close JSON document
        write(unit_out, '(a)') '}'
        
        close(unit_out)
    end subroutine convert_irfeatures

    !> Convert evolutionary track files to JSON
    subroutine convert_tracks(input_file, output_file)
        character(len=*), intent(in) :: input_file
        character(len=*), intent(in) :: output_file
        
        ! Local variables
        integer :: io_stat, unit_in, unit_out
        character(len=MAX_LINE_LEN) :: line, track_name
        character(len=:), allocatable :: errmsg
        integer :: i, j, num_rows, num_cols
        real(real64), allocatable :: track_data(:,:)
        
        ! Open input file
        open(newunit=unit_in, file=input_file, status='old', action='read', &
             iostat=io_stat, iomsg=errmsg)
        if (io_stat /= 0) then
            write(stderr, '(a)') 'Error opening input file: ' // trim(errmsg)
            return
        end if
        
        ! Read track name (first line)
        read(unit_in, '(a)', iostat=io_stat) track_name
        if (io_stat /= 0) then
            write(stderr, '(a)') 'Error reading track name: ' // trim(errmsg)
            close(unit_in)
            return
        end if
        
        ! Read dimensions (second line)
        read(unit_in, *, iostat=io_stat) num_cols, num_rows
        if (io_stat /= 0) then
            write(stderr, '(a)') 'Error reading dimensions: ' // trim(errmsg)
            close(unit_in)
            return
        end if
        
        ! Skip metadata lines (typically 3 more lines)
        do i = 1, 3
            read(unit_in, '(a)', iostat=io_stat) line
            if (io_stat /= 0) exit
        end do
        
        ! Allocate track data matrix
        allocate(track_data(num_rows, num_cols))
        
        ! Read track data
        do i = 1, num_rows
            read(unit_in, *, iostat=io_stat) j, (track_data(i,j), j=1, num_cols)
            if (io_stat /= 0) exit
        end do
        
        close(unit_in)
        
        ! Write JSON output
        open(newunit=unit_out, file=output_file, status='replace', action='write', &
             iostat=io_stat, iomsg=errmsg)
        if (io_stat /= 0) then
            write(stderr, '(a)') 'Error opening output file: ' // trim(errmsg)
            return
        end if
        
        ! Write JSON header
        write(unit_out, '(a)') '{'
        write(unit_out, '(a, a, a)') '  "name": "', trim(track_name), '",'
        write(unit_out, '(a, i0, a)') '  "columns": ', num_cols, ','
        write(unit_out, '(a, i0, a)') '  "rows": ', num_rows, ','
        
        ! Write track data as array of objects
        write(unit_out, '(a)') '  "tracks": ['
        do i = 1, num_rows
            write(unit_out, '(a)') '    {'
            write(unit_out, '(a, i0, a)') '      "index": ', i, ','
            
            ! Write parameters as an array
            write(unit_out, '(a)') '      "parameters": ['
            do j = 1, num_cols
                if (j < num_cols) then
                    write(unit_out, '(es16.8, a)') track_data(i,j), ','
                else
                    write(unit_out, '(es16.8)') track_data(i,j)
                end if
            end do
            write(unit_out, '(a)') '      ]'
            
            if (i < num_rows) then
                write(unit_out, '(a)') '    },'
            else
                write(unit_out, '(a)') '    }'
            end if
        end do
        write(unit_out, '(a)') '  ]'
        
        ! Close JSON document
        write(unit_out, '(a)') '}'
        
        close(unit_out)
    end subroutine convert_tracks

    !> Convert model files (mod*.dat) to JSON
    subroutine convert_modfiles(input_file, output_file)
        character(len=*), intent(in) :: input_file
        character(len=*), intent(in) :: output_file
        
        ! Local variables
        integer :: io_stat, unit_in, unit_out
        character(len=MAX_LINE_LEN) :: line
        character(len=:), allocatable :: errmsg
        integer :: i, j, num_rows, num_cols
        real(real64), allocatable :: mod_data(:,:)
        
        ! Open input file
        open(newunit=unit_in, file=input_file, status='old', action='read', &
             iostat=io_stat, iomsg=errmsg)
        if (io_stat /= 0) then
            write(stderr, '(a)') 'Error opening input file: ' // trim(errmsg)
            return
        end if
        
        ! Count number of rows and columns
        num_rows = 0
        num_cols = 0
        
        ! Read first line to determine number of columns
        read(unit_in, '(a)', iostat=io_stat) line
        if (io_stat == 0) then
            ! Count columns by reading values from the first line
            read(line, *, iostat=io_stat) (i, j=1, MAX_COLUMNS)
            num_cols = j - 1  ! Subtract 1 because the loop increments one more time
        end if
        
        ! Count rows
        rewind(unit_in)
        do
            read(unit_in, '(a)', iostat=io_stat) line
            if (io_stat /= 0) exit
            num_rows = num_rows + 1
        end do
        
        ! Allocate data array
        allocate(mod_data(num_rows, num_cols))
        
        ! Read data
        rewind(unit_in)
        do i = 1, num_rows
            read(unit_in, *, iostat=io_stat) (mod_data(i,j), j=1, num_cols)
            if (io_stat /= 0) exit
        end do
        
        close(unit_in)
        
        ! Write JSON output
        open(newunit=unit_out, file=output_file, status='replace', action='write', &
             iostat=io_stat, iomsg=errmsg)
        if (io_stat /= 0) then
            write(stderr, '(a)') 'Error opening output file: ' // trim(errmsg)
            return
        end if
        
        ! Write JSON
        write(unit_out, '(a)') '{'
        write(unit_out, '(a)') '  "data": ['
        
        do i = 1, num_rows
            write(unit_out, '(a)') '    ['
            do j = 1, num_cols
                if (j < num_cols) then
                    write(unit_out, '(es16.8, a)') mod_data(i,j), ','
                else
                    write(unit_out, '(es16.8)') mod_data(i,j)
                end if
            end do
            
            if (i < num_rows) then
                write(unit_out, '(a)') '    ],'
            else
                write(unit_out, '(a)') '    ]'
            end if
        end do
        
        write(unit_out, '(a)') '  ]'
        write(unit_out, '(a)') '}'
        
        close(unit_out)
    end subroutine convert_modfiles

end module galaxy_dat2json