program solution_converter
    use mesh, only: mesh_t
    use field, only: field_t
    use solution_io, only: solution_io_t
    use vtk_converter, only: vtk_converter_t
    implicit none
    
    ! Converter objects
    type(solution_io_t) :: sol_io
    type(vtk_converter_t) :: vtk_conv
    
    ! Command line arguments
    character(len=256) :: input_file, output_file, output_dir
    character(len=10) :: format
    integer :: num_args
    logical :: output_file_specified = .false.
    logical :: output_dir_specified = .false.
    logical :: format_specified = .false.
    
    ! Get command line arguments
    num_args = command_argument_count()
    
    ! If no arguments, print usage and exit
    if (num_args < 1) then
        call print_usage()
        stop
    end if
    
    ! Get input file (required)
    call get_command_argument(1, input_file)
    
    ! Process optional arguments
    if (num_args > 1) then
        call process_optional_args(2, num_args, output_file, output_dir, format, &
                                  output_file_specified, output_dir_specified, format_specified)
    end if
    
    ! Initialize converter objects
    call sol_io%init()
    
    ! Set output directory for VTK converter if specified
    if (output_dir_specified) then
        call vtk_conv%init(output_dir=output_dir)
    else
        call vtk_conv%init()
    end if
    
    ! Convert solution file
    if (output_file_specified) then
        call vtk_conv%convert_solution(input_file, output_file, sol_io)
    else
        call vtk_conv%convert_solution(input_file, solution_io=sol_io)
    end if
    
    print *, "Conversion completed successfully."
    
contains

    !> Print usage information
    subroutine print_usage()
        print *, "cuCFD Solution Converter"
        print *, "Usage: solution_converter input_file [options]"
        print *, ""
        print *, "Convert cuCFD solution files to VTK format for visualization."
        print *, ""
        print *, "Arguments:"
        print *, "  input_file             Path to input solution file (.bin)"
        print *, ""
        print *, "Options:"
        print *, "  -o, --output FILE      Specify output file name"
        print *, "  -d, --dir DIRECTORY    Specify output directory"
        print *, "  -f, --format FORMAT    Specify output format (vtk or vtu, default: vtk)"
        print *, ""
        print *, "Examples:"
        print *, "  solution_converter results/solution_000100.bin"
        print *, "  solution_converter results/solution_000100.bin -o temperature.vtk"
        print *, "  solution_converter results/solution_000100.bin -d paraview_data -f vtu"
    end subroutine print_usage
    
    !> Process optional command line arguments
    subroutine process_optional_args(start_idx, end_idx, output_file, output_dir, format, &
                                    output_file_specified, output_dir_specified, format_specified)
        integer, intent(in) :: start_idx, end_idx
        character(len=*), intent(out) :: output_file, output_dir, format
        logical, intent(out) :: output_file_specified, output_dir_specified, format_specified
        
        integer :: i
        character(len=256) :: arg
        
        i = start_idx
        do while (i <= end_idx)
            call get_command_argument(i, arg)
            
            select case (trim(arg))
                case ('-o', '--output')
                    ! Output file
                    if (i + 1 <= end_idx) then
                        i = i + 1
                        call get_command_argument(i, output_file)
                        output_file_specified = .true.
                    else
                        print *, "Error: Missing value for option ", trim(arg)
                        call print_usage()
                        stop
                    end if
                    
                case ('-d', '--dir')
                    ! Output directory
                    if (i + 1 <= end_idx) then
                        i = i + 1
                        call get_command_argument(i, output_dir)
                        output_dir_specified = .true.
                    else
                        print *, "Error: Missing value for option ", trim(arg)
                        call print_usage()
                        stop
                    end if
                    
                case ('-f', '--format')
                    ! Output format
                    if (i + 1 <= end_idx) then
                        i = i + 1
                        call get_command_argument(i, format)
                        format_specified = .true.
                        
                        ! Validate format
                        if (trim(format) /= 'vtk' .and. trim(format) /= 'vtu') then
                            print *, "Error: Unsupported format '", trim(format), "'. Must be 'vtk' or 'vtu'."
                            call print_usage()
                            stop
                        end if
                    else
                        print *, "Error: Missing value for option ", trim(arg)
                        call print_usage()
                        stop
                    end if
                    
                case default
                    print *, "Warning: Ignoring unknown option: ", trim(arg)
            end select
            
            i = i + 1
        end do
        
        ! Set default format if not specified
        if (.not. format_specified) then
            format = 'vtk'
        end if
    end subroutine process_optional_args

end program solution_converter 