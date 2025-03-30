module solution_io
    use mesh, only: mesh_t, WP
    use field, only: field_t, FIELD_TYPE_CELL, FIELD_TYPE_FACE, FIELD_TYPE_VERTEX
    implicit none
    
    !> Solution I/O data structure
    type :: solution_io_t
        character(len=256) :: output_dir = "results"   ! Output directory
        character(len=64)  :: prefix = "solution"      ! Filename prefix
        integer :: output_frequency = 100              ! Output frequency
        logical :: write_fields = .true.               ! Whether to write fields
        
    contains
        procedure :: init => solution_io_init
        procedure :: write_solution => solution_io_write_solution
        procedure :: read_solution => solution_io_read_solution
        procedure :: write_field => solution_io_write_field
        procedure :: read_field => solution_io_read_field
        procedure :: create_directory => solution_io_create_directory
    end type solution_io_t
    
contains

    !> Initialize solution I/O
    subroutine solution_io_init(this, output_dir, prefix, output_frequency)
        class(solution_io_t), intent(inout) :: this
        character(len=*), intent(in), optional :: output_dir
        character(len=*), intent(in), optional :: prefix
        integer, intent(in), optional :: output_frequency
        
        ! Set output directory
        if (present(output_dir)) then
            this%output_dir = output_dir
        end if
        
        ! Set filename prefix
        if (present(prefix)) then
            this%prefix = prefix
        end if
        
        ! Set output frequency
        if (present(output_frequency)) then
            this%output_frequency = output_frequency
        end if
        
        ! Create output directory if it doesn't exist
        call this%create_directory(this%output_dir)
    end subroutine solution_io_init
    
    !> Create directory if it doesn't exist 
    subroutine solution_io_create_directory(this, dir_path)
        class(solution_io_t), intent(inout) :: this
        character(len=*), intent(in) :: dir_path
        logical :: dir_exists
        integer :: stat
        
        print *, "DEBUG: Creating directory if it doesn't exist: '", trim(dir_path), "'"
        
        ! Check if directory exists
        inquire(file=trim(dir_path)//'/.', exist=dir_exists)
        
        ! Create directory if it doesn't exist
        if (.not. dir_exists) then
            print *, "DEBUG: Directory does not exist, creating: ", trim(dir_path)
            ! Use system call for directory creation - platform dependent
#ifdef _WIN32
            call execute_command_line('mkdir "' // trim(dir_path) // '"', wait=.true., exitstat=stat)
#else
            call execute_command_line('mkdir -p "' // trim(dir_path) // '"', wait=.true., exitstat=stat)
#endif
            if (stat /= 0) then
                print *, "Warning: Failed to create directory: ", trim(dir_path)
            else
                print *, "DEBUG: Directory created successfully: ", trim(dir_path)
            end if
        else
            print *, "DEBUG: Directory already exists: ", trim(dir_path)
        end if
    end subroutine solution_io_create_directory
    
    !> Write solution to file
    subroutine solution_io_write_solution(this, mesh, fields, iteration, time)
        class(solution_io_t), intent(inout) :: this
        type(mesh_t), intent(in) :: mesh
        type(field_t), intent(in) :: fields(:)
        integer, intent(in) :: iteration
        real(WP), intent(in) :: time
        integer :: i, unit_num, ios, comment_pos
        character(len=256) :: filename, filepath, header_file
        character(len=256) :: clean_dir, clean_prefix
        logical :: dir_exists
        
        print *, "DEBUG: Attempting to write solution for iteration ", iteration
        
        ! Clean up the output directory string to remove quotes, comments, and line endings
        clean_dir = trim(adjustl(this%output_dir))
        
        ! Remove surrounding quotes if present
        if (len_trim(clean_dir) >= 2) then
            if (clean_dir(1:1) == '"' .and. clean_dir(len_trim(clean_dir):len_trim(clean_dir)) == '"') then
                clean_dir = clean_dir(2:len_trim(clean_dir)-1)
            end if
        end if
        
        ! Remove any comment portion (text after #)
        comment_pos = index(clean_dir, '#')
        if (comment_pos > 0) then
            clean_dir = clean_dir(1:comment_pos-1)
        end if
        
        ! Remove any trailing newlines or carriage returns
        do while (len_trim(clean_dir) > 0 .and. (clean_dir(len_trim(clean_dir):len_trim(clean_dir)) == char(10) .or. &
                                                 clean_dir(len_trim(clean_dir):len_trim(clean_dir)) == char(13)))
            clean_dir = clean_dir(1:len_trim(clean_dir)-1)
        end do
        
        clean_dir = trim(adjustl(clean_dir))
        clean_prefix = trim(this%prefix)
        
        print *, "DEBUG: Using cleaned output directory: '", trim(clean_dir), "'"
        
        ! Make sure the directory exists
        call this%create_directory(clean_dir)
        
        ! Create filename with iteration number
        filename = trim(clean_prefix) // "_" // trim(int2str(iteration)) // ".bin"
        filepath = trim(clean_dir) // "/" // trim(filename)
        
        ! Create header file path
        header_file = trim(clean_dir) // "/" // trim(clean_prefix) // "_" // trim(int2str(iteration)) // ".hdr"
        
        print *, "DEBUG: Header file path: '", trim(header_file), "'"
        
        ! Check if we can write to the directory
        inquire(file=trim(clean_dir)//'/.', exist=dir_exists)
        if (.not. dir_exists) then
            print *, "Error: Directory does not exist: ", trim(clean_dir)
            print *, "Creating directory..."
            call this%create_directory(clean_dir)
        end if
        
        ! Try to create a test file to verify write permissions
        open(newunit=unit_num, file=trim(clean_dir)//"/test_write.tmp", status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            print *, "Error: Could not create test file in directory: ", trim(clean_dir), " (IOSTAT=", ios, ")"
            print *, "Please check directory permissions."
            return
        else
            close(unit_num, status='delete')
            print *, "DEBUG: Write permission test passed for directory: ", trim(clean_dir)
        end if
        
        open(newunit=unit_num, file=trim(header_file), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            print *, "Error: Could not open header file for writing: ", trim(header_file), " (IOSTAT=", ios, ")"
            return
        else
            print *, "DEBUG: Opened header file: ", trim(header_file)
        end if
        
        ! Write header info
        write(unit_num, *) "cuCFD Solution File"
        write(unit_num, *) "Version: 1.0"
        write(unit_num, *) "Iteration: ", iteration
        write(unit_num, *) "Time: ", time
        write(unit_num, *) "Mesh cells: ", mesh%topo%ncells
        write(unit_num, *) "Mesh faces: ", mesh%topo%nfaces
        write(unit_num, *) "Mesh vertices: ", mesh%topo%nvertices
        write(unit_num, *) "Number of fields: ", size(fields)
        
        ! Write field information
        do i = 1, size(fields)
            write(unit_num, *) "Field ", i, ": ", trim(fields(i)%name), ", Type: ", fields(i)%field_type, &
                              ", Components: ", fields(i)%ncomp
        end do
        
        ! Close header file
        close(unit_num)
        
        ! Open binary file for writing
        open(newunit=unit_num, file=trim(filepath), status='replace', form='unformatted', &
             access='stream', action='write', iostat=ios)
        if (ios /= 0) then
            print *, "Error: Could not open solution file for writing: ", trim(filepath), " (IOSTAT=", ios, ")"
            return
        else
            print *, "DEBUG: Opened solution file: ", trim(filepath)
        end if
        
        ! Write file signature and version
        write(unit_num) 'CUCFD001'  ! 8-byte signature
        
        ! Write solution metadata
        write(unit_num) iteration
        write(unit_num) time
        write(unit_num) size(fields)
        
        ! Write mesh information
        write(unit_num) mesh%topo%ncells
        write(unit_num) mesh%topo%nfaces
        write(unit_num) mesh%topo%nvertices
        
        ! Write mesh dimensions
        write(unit_num) mesh%nx, mesh%ny, mesh%nz
        write(unit_num) mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, mesh%zmin, mesh%zmax
        
        ! Write fields
        do i = 1, size(fields)
            ! First check if field is on device, copy from device if needed
            if (fields(i)%on_device) then
                print *, "DEBUG: Field ", trim(fields(i)%name), " is on device, should be copied first"
                ! We need a non-const reference to the field to call copy_from_device
                ! For now, let's assume fields are already copied from device if needed
                ! In a real implementation, this would need a better solution
            end if
            
            ! Write field metadata
            write(unit_num) len_trim(fields(i)%name)
            write(unit_num) fields(i)%name(1:len_trim(fields(i)%name))
            write(unit_num) fields(i)%field_type
            write(unit_num) fields(i)%ncomp
            write(unit_num) fields(i)%nentries
            
            ! Write field data
            write(unit_num) fields(i)%data
        end do
        
        ! Close file
        close(unit_num)
        
        print *, "Solution written to: ", trim(filepath)
        print *, "DEBUG: Finished writing solution for iteration ", iteration
    end subroutine solution_io_write_solution
    
    !> Read solution from file
    subroutine solution_io_read_solution(this, mesh, fields, iteration, time, filepath)
        class(solution_io_t), intent(inout) :: this
        type(mesh_t), intent(in) :: mesh
        type(field_t), allocatable, intent(inout) :: fields(:)
        integer, intent(out) :: iteration
        real(WP), intent(out) :: time
        character(len=*), intent(in) :: filepath
        integer :: i, unit_num, ios, nfields
        integer :: ncells, nfaces, nvertices
        integer :: field_type, ncomp, nentries, name_len
        character(len=64) :: field_name
        character(len=8) :: signature
        
        ! Open binary file for reading
        open(newunit=unit_num, file=trim(filepath), status='old', form='unformatted', &
             access='stream', action='read', iostat=ios)
        if (ios /= 0) then
            print *, "Error: Could not open solution file for reading: ", trim(filepath)
            return
        end if
        
        ! Read and verify file signature
        read(unit_num) signature
        if (signature /= 'CUCFD001') then
            print *, "Error: Invalid file signature in: ", trim(filepath)
            close(unit_num)
            return
        end if
        
        ! Read solution metadata
        read(unit_num) iteration
        read(unit_num) time
        read(unit_num) nfields
        
        ! Read and verify mesh information
        read(unit_num) ncells, nfaces, nvertices
        
        ! Verify the mesh dimensions match
        if (ncells /= mesh%topo%ncells .or. nfaces /= mesh%topo%nfaces .or. nvertices /= mesh%topo%nvertices) then
            print *, "Error: Mesh dimensions in file do not match current mesh."
            print *, "File: cells=", ncells, ", faces=", nfaces, ", vertices=", nvertices
            print *, "Mesh: cells=", mesh%topo%ncells, ", faces=", mesh%topo%nfaces, ", vertices=", mesh%topo%nvertices
            close(unit_num)
            return
        end if
        
        ! Skip mesh dimensions in the file
        read(unit_num) ! Skip nx, ny, nz
        read(unit_num) ! Skip xmin, xmax, ymin, ymax, zmin, zmax
        
        ! Allocate fields if needed
        if (allocated(fields)) then
            if (size(fields) /= nfields) then
                deallocate(fields)
                allocate(fields(nfields))
            end if
        else
            allocate(fields(nfields))
        end if
        
        ! Read fields
        do i = 1, nfields
            ! Read field metadata
            read(unit_num) name_len
            field_name = ''  ! Initialize to empty string
            read(unit_num) field_name(1:name_len)
            read(unit_num) field_type
            read(unit_num) ncomp
            read(unit_num) nentries
            
            ! Initialize the field
            call fields(i)%init(mesh, field_name(1:name_len), field_type, ncomp)
            
            ! Verify that field dimensions match
            if (fields(i)%nentries /= nentries) then
                print *, "Warning: Field '", trim(field_name), "' has different number of entries. Expected: ", &
                          nentries, ", got: ", fields(i)%nentries
            end if
            
            ! Read field data
            read(unit_num) fields(i)%data
        end do
        
        ! Close file
        close(unit_num)
        
        print *, "Solution read from: ", trim(filepath)
    end subroutine solution_io_read_solution
    
    !> Write a single field to file
    subroutine solution_io_write_field(this, field, iteration, filepath)
        class(solution_io_t), intent(inout) :: this
        type(field_t), intent(in) :: field
        integer, intent(in) :: iteration
        character(len=*), intent(in), optional :: filepath
        character(len=256) :: actual_filepath
        integer :: unit_num, ios
        
        ! Determine the file path
        if (present(filepath)) then
            actual_filepath = filepath
        else
            write(actual_filepath, '(A,"/",A,"_",I0.6,"_",A,".bin")') &
                trim(this%output_dir), trim(this%prefix), iteration, trim(field%name)
        end if
        
        ! Open binary file for writing
        open(newunit=unit_num, file=trim(actual_filepath), status='replace', form='unformatted', &
             access='stream', action='write', iostat=ios)
        if (ios /= 0) then
            print *, "Error: Could not open field file for writing: ", trim(actual_filepath)
            return
        end if
        
        ! Write file signature and version
        write(unit_num) 'CUFIELD1'  ! 8-byte signature
        
        ! Write field metadata
        write(unit_num) len_trim(field%name)
        write(unit_num) field%name(1:len_trim(field%name))
        write(unit_num) field%field_type
        write(unit_num) field%ncomp
        write(unit_num) field%nentries
        
        ! Write field data
        write(unit_num) field%data
        
        ! Close file
        close(unit_num)
        
        print *, "Field written to: ", trim(actual_filepath)
    end subroutine solution_io_write_field
    
    !> Read a single field from file
    subroutine solution_io_read_field(this, field, mesh, filepath)
        class(solution_io_t), intent(inout) :: this
        type(field_t), intent(inout) :: field
        type(mesh_t), intent(in) :: mesh
        character(len=*), intent(in) :: filepath
        integer :: unit_num, ios
        integer :: field_type, ncomp, nentries, name_len
        character(len=64) :: field_name
        character(len=8) :: signature
        
        ! Open binary file for reading
        open(newunit=unit_num, file=trim(filepath), status='old', form='unformatted', &
             access='stream', action='read', iostat=ios)
        if (ios /= 0) then
            print *, "Error: Could not open field file for reading: ", trim(filepath)
            return
        end if
        
        ! Read and verify file signature
        read(unit_num) signature
        if (signature /= 'CUFIELD1') then
            print *, "Error: Invalid field file signature in: ", trim(filepath)
            close(unit_num)
            return
        end if
        
        ! Read field metadata
        read(unit_num) name_len
        field_name = ''  ! Initialize to empty string
        read(unit_num) field_name(1:name_len)
        read(unit_num) field_type
        read(unit_num) ncomp
        read(unit_num) nentries
        
        ! Initialize the field
        call field%init(mesh, field_name(1:name_len), field_type, ncomp)
        
        ! Verify that field dimensions match
        if (field%nentries /= nentries) then
            print *, "Warning: Field '", trim(field_name), "' has different number of entries. Expected: ", &
                      nentries, ", got: ", field%nentries
        end if
        
        ! Read field data
        read(unit_num) field%data
        
        ! Close file
        close(unit_num)
        
        print *, "Field read from: ", trim(filepath)
    end subroutine solution_io_read_field

    !> Convert integer to string
    function int2str(n) result(str)
        integer, intent(in) :: n
        character(len=20) :: str
        
        write(str, '(I0)') n
    end function int2str

end module solution_io 