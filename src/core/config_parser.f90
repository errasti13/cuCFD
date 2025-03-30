module config_parser
    use mesh, only: WP
    implicit none
    
    ! Configuration data structure
    type :: config_section_t
        character(len=64) :: name = ""                ! Section name
        integer :: num_items = 0                       ! Number of key-value pairs
        character(len=64), allocatable :: keys(:)      ! Keys
        character(len=256), allocatable :: values(:)   ! Values
        type(config_section_t), pointer :: subsections(:) => null()  ! Nested sections - must be pointer for recursive type
        integer :: num_subsections = 0                ! Number of subsections
        integer :: parent_idx = 0                     ! Index in parent's subsections array
    end type config_section_t
    
    ! Main configuration data structure
    type :: config_t
        type(config_section_t) :: root                ! Root section
        logical :: is_loaded = .false.                 ! Whether config has been loaded
        character(len=512) :: filename = ""            ! Config filename
        
    contains
        procedure :: load => config_load
        procedure :: get_string => config_get_string
        procedure :: get_real => config_get_real
        procedure :: get_integer => config_get_integer
        procedure :: get_logical => config_get_logical
        procedure :: get_real_array => config_get_real_array
        procedure :: get_integer_array => config_get_integer_array
        procedure :: get_section => config_get_section
        procedure :: has_key => config_has_key
        procedure :: has_section => config_has_section
        procedure :: print => config_print
        procedure :: cleanup => config_cleanup
    end type config_t
    
contains

    !> Load configuration from file
    subroutine config_load(this, filename)
        class(config_t), intent(inout) :: this
        character(len=*), intent(in) :: filename
        integer :: unit, ios, i, line_num
        character(len=512) :: line, clean_line
        
        ! Initialize root section
        this%root%name = "root"
        
        ! Store filename
        this%filename = filename
        
        ! Open configuration file
        open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, "Error: Could not open configuration file: ", trim(filename)
            return
        end if
        
        print *, "Reading configuration from: ", trim(filename)
        
        ! Parse file line by line and build configuration tree
        line_num = 0
        call parse_section(unit, this%root, 0, line_num)
        
        ! Close file
        close(unit)
        
        ! Mark as loaded
        this%is_loaded = .true.
        
        ! Print configuration tree
        ! call this%print()
    end subroutine config_load
    
    !> Parse a section of the configuration file
    recursive subroutine parse_section(unit, section, indent_level, line_num)
        integer, intent(in) :: unit
        type(config_section_t), intent(inout) :: section
        integer, intent(in) :: indent_level
        integer, intent(inout) :: line_num
        
        character(len=512) :: line, clean_line
        integer :: ios, colon_pos, i, j, actual_indent, next_indent
        character(len=64) :: key
        character(len=256) :: value
        logical :: is_section, is_array_item
        
        ! Temporary storage for key-value pairs
        character(len=64), allocatable :: temp_keys(:)
        character(len=256), allocatable :: temp_values(:)
        
        ! Temporary storage for subsections
        type(config_section_t), pointer :: temp_subsections(:) => null()
        
        ! Read file until end or until we encounter a line with lower indentation
        do
            ! Read line
            read(unit, '(A)', iostat=ios) line
            line_num = line_num + 1
            
            ! Exit if end of file
            if (ios /= 0) exit
            
            ! Clean up line
            clean_line = trim(adjustl(line))
            
            ! Skip empty lines and comments
            if (len_trim(clean_line) == 0 .or. clean_line(1:1) == '#') then
                cycle
            end if
            
            ! Determine indentation level
            actual_indent = 0
            do i = 1, len(line)
                if (line(i:i) == ' ') then
                    actual_indent = actual_indent + 1
                else
                    exit
                end if
            end do
            
            ! Convert spaces to indent level (assuming 2 spaces per level)
            actual_indent = actual_indent / 2
            
            ! If indentation is less than current level, we're done with this section
            if (actual_indent < indent_level) then
                ! Backtrack to previous line
                backspace(unit)
                line_num = line_num - 1
                exit
            end if
            
            ! Skip if indentation doesn't match current level or next level
            if (actual_indent > indent_level + 1) then
                cycle
            end if
            
            ! Find colon position (key-value separator)
            colon_pos = index(clean_line, ':')
            
            ! Skip line if no colon found
            if (colon_pos == 0) then
                cycle
            end if
            
            ! Extract key and value
            key = trim(adjustl(clean_line(1:colon_pos-1)))
            value = ""
            
            ! Check if line defines a subsection
            is_section = .false.
            if (colon_pos == len_trim(clean_line)) then
                ! This is a section header
                is_section = .true.
            else
                ! Extract value
                value = trim(adjustl(clean_line(colon_pos+1:)))
                
                ! Remove comments from value
                i = index(value, '#')
                if (i > 0) then
                    value = trim(adjustl(value(1:i-1)))
                end if
                
                ! Remove quotes if present
                if (len_trim(value) >= 2) then
                    if ((value(1:1) == '"' .and. value(len_trim(value):len_trim(value)) == '"') .or. &
                        (value(1:1) == "'" .and. value(len_trim(value):len_trim(value)) == "'")) then
                        value = value(2:len_trim(value)-1)
                    end if
                end if
                
            end if
            
            ! Process based on whether it's a section or key-value pair
            if (is_section) then
                ! Add subsection
                if (.not. associated(section%subsections)) then
                    allocate(section%subsections(1))
                    section%num_subsections = 1
                else
                    ! Resize subsections array
                    allocate(temp_subsections(section%num_subsections + 1))
                    temp_subsections(1:section%num_subsections) = section%subsections
                    
                    ! Deallocate old array and assign new one
                    deallocate(section%subsections)
                    section%subsections => temp_subsections
                    nullify(temp_subsections)
                    
                    section%num_subsections = section%num_subsections + 1
                end if
                
                ! Initialize new subsection
                section%subsections(section%num_subsections)%name = key
                section%subsections(section%num_subsections)%parent_idx = section%num_subsections
                
                ! Parse subsection recursively
                call parse_section(unit, section%subsections(section%num_subsections), indent_level + 1, line_num)
            else
                ! Add key-value pair to current section
                if (.not. allocated(section%keys)) then
                    allocate(section%keys(1))
                    allocate(section%values(1))
                    section%num_items = 1
                else
                    ! Resize key-value arrays
                    allocate(temp_keys(section%num_items + 1))
                    allocate(temp_values(section%num_items + 1))
                    temp_keys(1:section%num_items) = section%keys
                    temp_values(1:section%num_items) = section%values
                    call move_alloc(temp_keys, section%keys)
                    call move_alloc(temp_values, section%values)
                    section%num_items = section%num_items + 1
                end if
                
                ! Store key-value pair
                section%keys(section%num_items) = key
                section%values(section%num_items) = value
            end if
        end do
    end subroutine parse_section
    
    !> Get string value from configuration
    function config_get_string(this, key_path, default_value) result(value)
        class(config_t), intent(in) :: this
        character(len=*), intent(in) :: key_path
        character(len=*), intent(in), optional :: default_value
        character(len=256) :: value
        
        type(config_section_t), pointer :: section
        character(len=64) :: section_path(10), key
        integer :: num_sections, i
        logical :: need_cleanup = .false.
        
        ! Initialize with default value
        if (present(default_value)) then
            value = default_value
        else
            value = ""
        end if
        
        ! Check if configuration is loaded
        if (.not. this%is_loaded) then
            return
        end if
        
        ! Split key path into sections
        call split_key_path(key_path, section_path, num_sections, key)
        
        ! Navigate to section
        nullify(section)
        
        ! Special case for root section
        if (num_sections == 0) then
            allocate(section)
            section = this%root
            need_cleanup = .true.
        else
            ! For non-root sections, navigate through the structure
            do i = 1, this%root%num_subsections
                if (trim(this%root%subsections(i)%name) == trim(section_path(1))) then
                    if (num_sections == 1) then
                        section => this%root%subsections(i)
                    else
                        call navigate_to_section(this%root%subsections(i), section_path(2:), num_sections-1, section)
                    end if
                    exit
                end if
            end do
        end if
        
        ! If section not found, return default
        if (.not. associated(section)) return
        
        ! Find key in section
        do i = 1, section%num_items
            if (trim(section%keys(i)) == trim(key)) then
                value = section%values(i)
                exit
            end if
        end do
        
        ! Cleanup if we allocated a copy
        if (need_cleanup) then
            deallocate(section)
        end if
        nullify(section)
    end function config_get_string
    
    !> Get real value from configuration
    function config_get_real(this, key_path, default_value) result(value)
        class(config_t), intent(in) :: this
        character(len=*), intent(in) :: key_path
        real(WP), intent(in), optional :: default_value
        real(WP) :: value
        
        character(len=256) :: str_value
        integer :: ios
        
        ! Initialize with default value
        if (present(default_value)) then
            value = default_value
        else
            value = 0.0_WP
        end if
        
        ! Get string value
        str_value = this%get_string(key_path)
        
        ! Convert to real
        if (len_trim(str_value) > 0) then
            read(str_value, *, iostat=ios) value
            if (ios /= 0) then
                print *, "Warning: Could not convert '", trim(str_value), "' to real for key ", trim(key_path)
                if (present(default_value)) value = default_value
            end if
        end if
    end function config_get_real
    
    !> Get integer value from configuration
    function config_get_integer(this, key_path, default_value) result(value)
        class(config_t), intent(in) :: this
        character(len=*), intent(in) :: key_path
        integer, intent(in), optional :: default_value
        integer :: value
        
        character(len=256) :: str_value
        integer :: ios
        
        ! Initialize with default value
        if (present(default_value)) then
            value = default_value
        else
            value = 0
        end if
        
        ! Get string value
        str_value = this%get_string(key_path)
        
        ! Convert to integer
        if (len_trim(str_value) > 0) then
            read(str_value, *, iostat=ios) value
            if (ios /= 0) then
                print *, "Warning: Could not convert '", trim(str_value), "' to integer for key ", trim(key_path)
                if (present(default_value)) value = default_value
            end if
        end if
    end function config_get_integer
    
    !> Get logical value from configuration
    function config_get_logical(this, key_path, default_value) result(value)
        class(config_t), intent(in) :: this
        character(len=*), intent(in) :: key_path
        logical, intent(in), optional :: default_value
        logical :: value
        
        character(len=256) :: str_value
        
        ! Initialize with default value
        if (present(default_value)) then
            value = default_value
        else
            value = .false.
        end if
        
        ! Get string value
        str_value = this%get_string(key_path)
        
        ! Convert to logical
        if (len_trim(str_value) > 0) then
            ! Convert various string representations to logical
            if (str_value == "true" .or. str_value == "yes" .or. str_value == "on" .or. &
                str_value == "True" .or. str_value == "Yes" .or. str_value == "On" .or. &
                str_value == "TRUE" .or. str_value == "YES" .or. str_value == "ON" .or. &
                str_value == "1" .or. str_value == ".true.") then
                value = .true.
            else if (str_value == "false" .or. str_value == "no" .or. str_value == "off" .or. &
                    str_value == "False" .or. str_value == "No" .or. str_value == "Off" .or. &
                    str_value == "FALSE" .or. str_value == "NO" .or. str_value == "OFF" .or. &
                    str_value == "0" .or. str_value == ".false.") then
                value = .false.
            else
                print *, "Warning: Could not convert '", trim(str_value), "' to logical for key ", trim(key_path)
                if (present(default_value)) value = default_value
            end if
        end if
    end function config_get_logical
    
    !> Get real array from configuration
    function config_get_real_array(this, key_path, default_value) result(array)
        class(config_t), intent(in) :: this
        character(len=*), intent(in) :: key_path
        real(WP), intent(in), optional :: default_value(:)
        real(WP), allocatable :: array(:)
        
        character(len=256) :: str_value
        character(len=256) :: array_str
        integer :: i, j, k, num_values, ios
        
        ! Get string value
        str_value = this%get_string(key_path)
        
        ! Check if array notation [a, b, c]
        if (len_trim(str_value) > 2 .and. str_value(1:1) == '[' .and. &
            str_value(len_trim(str_value):len_trim(str_value)) == ']') then
            
            ! Remove brackets
            array_str = str_value(2:len_trim(str_value)-1)
            
            ! Count number of values (count commas + 1)
            num_values = 1
            do i = 1, len_trim(array_str)
                if (array_str(i:i) == ',') num_values = num_values + 1
            end do
            
            ! Allocate array
            allocate(array(num_values))
            
            ! Parse values
            j = 1
            k = 1
            do i = 1, len_trim(array_str)
                if (array_str(i:i) == ',' .or. i == len_trim(array_str)) then
                    ! Extract value
                    if (i == len_trim(array_str)) then
                        read(array_str(j:i), *, iostat=ios) array(k)
                    else
                        read(array_str(j:i-1), *, iostat=ios) array(k)
                    end if
                    
                    ! Check for errors
                    if (ios /= 0) then
                        print *, "Warning: Could not parse array element for key ", trim(key_path)
                        array(k) = 0.0_WP
                    end if
                    
                    ! Move to next value
                    j = i + 1
                    k = k + 1
                end if
            end do
        else
            ! Single value
            if (len_trim(str_value) > 0) then
                allocate(array(1))
                read(str_value, *, iostat=ios) array(1)
                if (ios /= 0) then
                    print *, "Warning: Could not convert '", trim(str_value), "' to real for key ", trim(key_path)
                    array(1) = 0.0_WP
                end if
            else if (present(default_value)) then
                ! Use default value
                allocate(array(size(default_value)))
                array = default_value
            else
                ! Empty array
                allocate(array(0))
            end if
        end if
    end function config_get_real_array
    
    !> Get integer array from configuration
    function config_get_integer_array(this, key_path, default_value) result(array)
        class(config_t), intent(in) :: this
        character(len=*), intent(in) :: key_path
        integer, intent(in), optional :: default_value(:)
        integer, allocatable :: array(:)
        
        character(len=256) :: str_value
        character(len=256) :: array_str
        integer :: i, j, k, num_values, ios
        
        ! Get string value
        str_value = this%get_string(key_path)
        
        ! Check if array notation [a, b, c]
        if (len_trim(str_value) > 2 .and. str_value(1:1) == '[' .and. &
            str_value(len_trim(str_value):len_trim(str_value)) == ']') then
            
            ! Remove brackets
            array_str = str_value(2:len_trim(str_value)-1)
            
            ! Count number of values (count commas + 1)
            num_values = 1
            do i = 1, len_trim(array_str)
                if (array_str(i:i) == ',') num_values = num_values + 1
            end do
            
            ! Allocate array
            allocate(array(num_values))
            
            ! Parse values
            j = 1
            k = 1
            do i = 1, len_trim(array_str)
                if (array_str(i:i) == ',' .or. i == len_trim(array_str)) then
                    ! Extract value
                    if (i == len_trim(array_str)) then
                        read(array_str(j:i), *, iostat=ios) array(k)
                    else
                        read(array_str(j:i-1), *, iostat=ios) array(k)
                    end if
                    
                    ! Check for errors
                    if (ios /= 0) then
                        print *, "Warning: Could not parse array element for key ", trim(key_path)
                        array(k) = 0
                    end if
                    
                    ! Move to next value
                    j = i + 1
                    k = k + 1
                end if
            end do
        else
            ! Single value
            if (len_trim(str_value) > 0) then
                allocate(array(1))
                read(str_value, *, iostat=ios) array(1)
                if (ios /= 0) then
                    print *, "Warning: Could not convert '", trim(str_value), "' to integer for key ", trim(key_path)
                    array(1) = 0
                end if
            else if (present(default_value)) then
                ! Use default value
                allocate(array(size(default_value)))
                array = default_value
            else
                ! Empty array
                allocate(array(0))
            end if
        end if
    end function config_get_integer_array
    
    !> Get configuration section
    subroutine config_get_section(this, section_path, section)
        class(config_t), intent(in) :: this
        character(len=*), intent(in) :: section_path
        type(config_section_t), pointer, intent(out) :: section
        
        character(len=64) :: path_parts(10)
        integer :: num_parts, i
        
        ! Initialize
        nullify(section)
        
        ! Check if configuration is loaded
        if (.not. this%is_loaded) then
            print *, "Warning: Configuration not loaded. Cannot get section ", trim(section_path)
            return
        end if
        
        ! Special case for root section
        if (len_trim(section_path) == 0 .or. section_path == "root") then
            ! We can't directly assign a non-pointer to a pointer, so need to create a temporary solution
            ! This is a workaround since we can't use 'target' attribute
            allocate(section)
            section = this%root
            return
        end if
        
        ! Split section path
        call split_section_path(section_path, path_parts, num_parts)
        
        ! For non-root sections, navigate through the structure
        ! We need to modify the navigate_to_section logic to handle this special case
        if (num_parts > 0) then
            ! Look for the first section in root's subsections
            do i = 1, this%root%num_subsections
                if (trim(this%root%subsections(i)%name) == trim(path_parts(1))) then
                    if (num_parts == 1) then
                        section => this%root%subsections(i)
                    else
                        call navigate_to_section(this%root%subsections(i), path_parts(2:), num_parts-1, section)
                    end if
                    return
                end if
            end do
        end if
    end subroutine config_get_section
    
    !> Check if key exists in configuration
    function config_has_key(this, key_path) result(exists)
        class(config_t), intent(in) :: this
        character(len=*), intent(in) :: key_path
        logical :: exists
        
        type(config_section_t), pointer :: section
        character(len=64) :: section_path(10), key
        integer :: num_sections, i
        
        ! Initialize
        exists = .false.
        
        ! Check if configuration is loaded
        if (.not. this%is_loaded) return
        
        ! Split key path into sections
        call split_key_path(key_path, section_path, num_sections, key)
        
        ! Navigate to section
        nullify(section)
        
        ! Special case for root section
        if (num_sections == 0) then
            ! We need to allocate and copy for root
            allocate(section)
            section = this%root
        else
            ! For non-root sections, navigate through the structure
            do i = 1, this%root%num_subsections
                if (trim(this%root%subsections(i)%name) == trim(section_path(1))) then
                    if (num_sections == 1) then
                        section => this%root%subsections(i)
                    else
                        call navigate_to_section(this%root%subsections(i), section_path(2:), num_sections-1, section)
                    end if
                    exit
                end if
            end do
        end if
        
        ! If section not found, return false
        if (.not. associated(section)) return
        
        ! Find key in section
        do i = 1, section%num_items
            if (trim(section%keys(i)) == trim(key)) then
                exists = .true.
                exit
            end if
        end do
        
        ! Cleanup if we allocated a copy
        if (num_sections == 0 .and. associated(section)) then
            deallocate(section)
            nullify(section)
        end if
    end function config_has_key
    
    !> Check if section exists in configuration
    function config_has_section(this, section_path) result(exists)
        class(config_t), intent(in) :: this
        character(len=*), intent(in) :: section_path
        logical :: exists
        
        type(config_section_t), pointer :: section
        
        ! Initialize
        exists = .false.
        
        ! Check if configuration is loaded
        if (.not. this%is_loaded) return
        
        ! Special case for root section
        if (len_trim(section_path) == 0 .or. section_path == "root") then
            exists = .true.
            return
        end if
        
        ! Try to get section
        nullify(section)
        call this%get_section(section_path, section)
        
        ! Check if section exists
        exists = associated(section)
        
        ! Cleanup section pointer if needed
        if (associated(section)) then
            ! We only need to deallocate if it's a copy of the root
            ! Pointers to subsections don't need deallocation
            if (section_path == "root" .or. len_trim(section_path) == 0) then
                deallocate(section)
            end if
            nullify(section)
        end if
    end function config_has_section
    
    !> Print configuration tree
    subroutine config_print(this)
        class(config_t), intent(in) :: this
        
        ! Check if configuration is loaded
        if (.not. this%is_loaded) then
            print *, "Configuration not loaded."
            return
        end if
        
        ! Print root section recursively
        print *, ""
        print *, "Configuration from file: ", trim(this%filename)
        print *, "------------------------------------------------------------"
        call print_section(this%root, 0)
        print *, "------------------------------------------------------------"
        print *, ""
    end subroutine config_print
    
    !> Print configuration section recursively
    recursive subroutine print_section(section, indent)
        type(config_section_t), intent(in) :: section
        integer, intent(in) :: indent
        
        character(len=256) :: indent_str
        integer :: i
        
        ! Create indentation string
        indent_str = ""
        do i = 1, indent
            indent_str = trim(indent_str) // "  "
        end do
        
        ! Print section name
        if (trim(section%name) /= "root") then
            print *, trim(indent_str), "[", trim(section%name), "]"
        end if
        
        ! Print key-value pairs
        do i = 1, section%num_items
            print *, trim(indent_str), "  ", trim(section%keys(i)), " = ", trim(section%values(i))
        end do
        
        ! Print subsections
        do i = 1, section%num_subsections
            if (i > 1 .or. section%num_items > 0) print *, "" ! Add space for better readability
            call print_section(section%subsections(i), indent + 1)
        end do
    end subroutine print_section
    
    !> Cleanup configuration
    subroutine config_cleanup(this)
        class(config_t), intent(inout) :: this
        
        ! Cleanup root section recursively
        call cleanup_section(this%root)
        
        ! Reset state
        this%is_loaded = .false.
        this%filename = ""
    end subroutine config_cleanup
    
    !> Cleanup configuration section recursively
    recursive subroutine cleanup_section(section)
        type(config_section_t), intent(inout) :: section
        integer :: i
        
        ! Cleanup subsections
        if (associated(section%subsections)) then
            do i = 1, section%num_subsections
                call cleanup_section(section%subsections(i))
            end do
            deallocate(section%subsections)
            nullify(section%subsections)
        end if
        
        ! Cleanup key-value arrays
        if (allocated(section%keys)) deallocate(section%keys)
        if (allocated(section%values)) deallocate(section%values)
        
        ! Reset state
        section%num_items = 0
        section%num_subsections = 0
    end subroutine cleanup_section
    
    !> Split key path into section path and key
    subroutine split_key_path(key_path, section_path, num_sections, key)
        character(len=*), intent(in) :: key_path
        character(len=*), intent(out) :: section_path(:)
        integer, intent(out) :: num_sections
        character(len=*), intent(out) :: key
        
        integer :: i, j, last_dot
        
        ! Find last dot in key path
        last_dot = 0
        do i = 1, len_trim(key_path)
            if (key_path(i:i) == '.') last_dot = i
        end do
        
        ! If no dot, key is in root section
        if (last_dot == 0) then
            num_sections = 0
            key = key_path
            return
        end if
        
        ! Extract key
        key = key_path(last_dot+1:)
        
        ! Split section path on dots
        num_sections = 1
        j = 1
        do i = 1, last_dot - 1
            if (key_path(i:i) == '.') then
                section_path(num_sections) = key_path(j:i-1)
                num_sections = num_sections + 1
                j = i + 1
            end if
        end do
        
        ! Add last section
        section_path(num_sections) = key_path(j:last_dot-1)
    end subroutine split_key_path
    
    !> Split section path into parts
    subroutine split_section_path(section_path, path_parts, num_parts)
        character(len=*), intent(in) :: section_path
        character(len=*), intent(out) :: path_parts(:)
        integer, intent(out) :: num_parts
        
        integer :: i, j
        
        ! Split section path on dots
        num_parts = 1
        j = 1
        do i = 1, len_trim(section_path)
            if (section_path(i:i) == '.') then
                path_parts(num_parts) = section_path(j:i-1)
                num_parts = num_parts + 1
                j = i + 1
            end if
        end do
        
        ! Add last part
        path_parts(num_parts) = section_path(j:)
    end subroutine split_section_path
    
    !> Navigate to section in configuration tree
    recursive subroutine navigate_to_section(current, section_path, num_sections, result_section)
        type(config_section_t), intent(in) :: current
        character(len=*), intent(in) :: section_path(:)
        integer, intent(in) :: num_sections
        type(config_section_t), pointer, intent(out) :: result_section
        
        integer :: i
        
        ! If no more sections to navigate, we've found the right section
        if (num_sections == 0) then
            ! We can't directly assign a non-pointer to a pointer, so allocate and copy
            allocate(result_section)
            result_section = current
            return
        end if
        
        ! Find subsection with matching name
        do i = 1, current%num_subsections
            if (trim(current%subsections(i)%name) == trim(section_path(1))) then
                ! Navigate to next section
                if (num_sections == 1) then
                    result_section => current%subsections(i)
                else
                    call navigate_to_section(current%subsections(i), section_path(2:), num_sections - 1, result_section)
                end if
                return
            end if
        end do
        
        ! Section not found, return null
        nullify(result_section)
    end subroutine navigate_to_section

end module config_parser 