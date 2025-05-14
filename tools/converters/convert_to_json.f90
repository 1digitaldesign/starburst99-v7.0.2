!==============================================================================
! DAT to JSON Converter Program
!==============================================================================
!> Program to convert GALAXY/Starburst99 data files to JSON format
program convert_to_json
    use galaxy_dat2json
    use, intrinsic :: iso_fortran_env, only: stdout => output_unit, &
                                           stderr => error_unit
    implicit none
    
    ! Local variables
    character(len=1024) :: input_file, output_file
    character(len=10) :: file_type
    integer :: num_args, i, status
    character(len=256) :: arg
    logical :: help_requested = .false.
    
    ! Process command line arguments
    num_args = command_argument_count()
    
    if (num_args < 2) then
        help_requested = .true.
    else
        ! Get file type
        call get_command_argument(1, file_type, status=status)
        if (status /= 0) then
            write(stderr, '(a)') 'Error: Invalid file type argument'
            help_requested = .true.
        end if
        
        ! Get input file
        call get_command_argument(2, input_file, status=status)
        if (status /= 0) then
            write(stderr, '(a)') 'Error: Invalid input file argument'
            help_requested = .true.
        end if
        
        ! Get output file (optional)
        if (num_args >= 3) then
            call get_command_argument(3, output_file, status=status)
            if (status /= 0) then
                write(stderr, '(a)') 'Error: Invalid output file argument'
                help_requested = .true.
            end if
        else
            ! Default output file name (replace .dat with .json)
            output_file = trim(input_file)
            i = index(output_file, '.dat')
            if (i > 0) then
                output_file = output_file(1:i-1) // '.json'
            else
                output_file = trim(output_file) // '.json'
            end if
        end if
    end if
    
    ! Display help if requested or no arguments
    if (help_requested) then
        write(stdout, '(a)') 'Usage: convert_to_json <type> <input_file> [output_file]'
        write(stdout, '(a)') ''
        write(stdout, '(a)') 'Convert GALAXY/Starburst99 data files to JSON format'
        write(stdout, '(a)') ''
        write(stdout, '(a)') 'Arguments:'
        write(stdout, '(a)') '  type        File type: irfeatures, tracks, modfiles'
        write(stdout, '(a)') '  input_file  Input data file path'
        write(stdout, '(a)') '  output_file Output JSON file path (optional)'
        write(stdout, '(a)') ''
        write(stdout, '(a)') 'Examples:'
        write(stdout, '(a)') '  convert_to_json irfeatures auxil/irfeatures.dat'
        write(stdout, '(a)') '  convert_to_json tracks tracks/Z0020v00.txt output.json'
        write(stdout, '(a)') '  convert_to_json modfiles tracks/modc001.dat'
        stop
    end if
    
    ! Convert files based on type
    select case (trim(file_type))
        case ('irfeatures')
            call convert_irfeatures(trim(input_file), trim(output_file))
            
        case ('tracks')
            call convert_tracks(trim(input_file), trim(output_file))
            
        case ('modfiles')
            call convert_modfiles(trim(input_file), trim(output_file))
            
        case default
            write(stderr, '(a)') 'Error: Unknown file type. Must be one of: irfeatures, tracks, modfiles'
            stop 1
    end select
    
    write(stdout, '(a, a, a)') 'Successfully converted ', trim(input_file), ' to JSON format'
    write(stdout, '(a, a)') 'Output file: ', trim(output_file)
    
end program convert_to_json