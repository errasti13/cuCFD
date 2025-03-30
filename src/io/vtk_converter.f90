module vtk_converter
    use mesh, only: mesh_t, WP
    use field, only: field_t, FIELD_TYPE_CELL, FIELD_TYPE_FACE, FIELD_TYPE_VERTEX
    use solution_io, only: solution_io_t
    implicit none
    
    !> VTK converter data structure
    type :: vtk_converter_t
        character(len=256) :: output_dir = "vtk_output"  ! Output directory
        character(len=64)  :: prefix = "solution"        ! Filename prefix
        
    contains
        procedure :: init => vtk_converter_init
        procedure :: convert_solution => vtk_converter_convert_solution
        procedure :: write_vtk_file => vtk_converter_write_vtk_file
        procedure :: write_vtu_file => vtk_converter_write_vtu_file
        procedure :: create_directory => vtk_converter_create_directory
    end type vtk_converter_t
    
contains

    !> Initialize VTK converter
    subroutine vtk_converter_init(this, output_dir, prefix)
        class(vtk_converter_t), intent(inout) :: this
        character(len=*), intent(in), optional :: output_dir
        character(len=*), intent(in), optional :: prefix
        
        ! Set output directory
        if (present(output_dir)) then
            this%output_dir = output_dir
        end if
        
        ! Set filename prefix
        if (present(prefix)) then
            this%prefix = prefix
        end if
        
        ! Create output directory if it doesn't exist
        call this%create_directory(this%output_dir)
    end subroutine vtk_converter_init
    
    !> Create directory if it doesn't exist
    subroutine vtk_converter_create_directory(this, dir_path)
        class(vtk_converter_t), intent(inout) :: this
        character(len=*), intent(in) :: dir_path
        logical :: dir_exists
        integer :: stat
        
        ! Check if directory exists
        inquire(file=trim(dir_path)//'/.', exist=dir_exists)
        
        ! Create directory if it doesn't exist
        if (.not. dir_exists) then
            ! Use system call for directory creation - platform dependent
#ifdef _WIN32
            call execute_command_line('mkdir "' // trim(dir_path) // '"', wait=.true., exitstat=stat)
#else
            call execute_command_line('mkdir -p "' // trim(dir_path) // '"', wait=.true., exitstat=stat)
#endif
            if (stat /= 0) then
                print *, "Warning: Failed to create directory: ", trim(dir_path)
            end if
        end if
    end subroutine vtk_converter_create_directory
    
    !> Convert solution file to VTK format
    subroutine vtk_converter_convert_solution(this, solution_file, vtk_file, solution_io)
        class(vtk_converter_t), intent(inout) :: this
        character(len=*), intent(in) :: solution_file
        character(len=*), intent(in), optional :: vtk_file
        type(solution_io_t), intent(inout), optional :: solution_io
        
        type(mesh_t) :: mesh
        type(field_t), allocatable :: fields(:)
        type(solution_io_t) :: local_io
        character(len=256) :: output_file
        integer :: iteration
        real(WP) :: time
        
        ! Initialize mesh
        call mesh%init()
        
        ! Create solution I/O object if not provided
        if (.not. present(solution_io)) then
            call local_io%init()
        else
            local_io = solution_io
        end if
        
        ! Determine output file path
        if (present(vtk_file)) then
            output_file = vtk_file
        else
            ! Extract iteration number from filename (assuming format like solution_000100.bin)
            ! This is a simplification - in a real implementation you'd parse the filename more robustly
            output_file = trim(this%output_dir) // '/' // trim(this%prefix) // '_' // &
                          trim(solution_file(index(solution_file, '_')+1:index(solution_file, '.bin')-1)) // '.vtk'
        end if
        
        ! Read solution file
        call local_io%read_solution(mesh, fields, iteration, time, solution_file)
        
        ! Write VTK file
        call this%write_vtk_file(mesh, fields, output_file, iteration, time)
        
        ! Clean up
        if (allocated(fields)) then
            do iteration = 1, size(fields)
                call fields(iteration)%cleanup()
            end do
            deallocate(fields)
        end if
        call mesh%cleanup()
        
        print *, "Solution converted to VTK: ", trim(output_file)
    end subroutine vtk_converter_convert_solution
    
    !> Write VTK legacy format file
    subroutine vtk_converter_write_vtk_file(this, mesh, fields, vtk_file, iteration, time)
        class(vtk_converter_t), intent(inout) :: this
        type(mesh_t), intent(in) :: mesh
        type(field_t), intent(in) :: fields(:)
        character(len=*), intent(in) :: vtk_file
        integer, intent(in) :: iteration
        real(WP), intent(in) :: time
        
        integer :: unit_num, ios, i, j, k, idx
        
        ! Open VTK file for writing
        open(newunit=unit_num, file=trim(vtk_file), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            print *, "Error: Could not open VTK file for writing: ", trim(vtk_file)
            return
        end if
        
        ! Write VTK header
        write(unit_num, '(A)') '# vtk DataFile Version 3.0'
        write(unit_num, '(A,I0,A,F10.6)') 'cuCFD Solution - Iteration: ', iteration, ', Time: ', time
        write(unit_num, '(A)') 'ASCII'
        write(unit_num, '(A)') 'DATASET STRUCTURED_GRID'
        write(unit_num, '(A,3(1X,I0))') 'DIMENSIONS', mesh%nx+1, mesh%ny+1, mesh%nz+1
        write(unit_num, '(A,1X,I0,1X,A)') 'POINTS', (mesh%nx+1)*(mesh%ny+1)*(mesh%nz+1), 'float'
        
        ! Write vertex coordinates
        do k = 0, mesh%nz
            do j = 0, mesh%ny
                do i = 0, mesh%nx
                    ! Calculate vertex coordinates based on domain dimensions
                    write(unit_num, '(3(F12.6,1X))') &
                        mesh%xmin + i * (mesh%xmax - mesh%xmin) / mesh%nx, &
                        mesh%ymin + j * (mesh%ymax - mesh%ymin) / mesh%ny, &
                        mesh%zmin + k * (mesh%zmax - mesh%zmin) / mesh%nz
                end do
            end do
        end do
        
        ! Write cell data
        write(unit_num, '(A,1X,I0)') 'CELL_DATA', mesh%nx * mesh%ny * mesh%nz
        
        ! Write fields
        do i = 1, size(fields)
            ! Only write cell-centered fields
            if (fields(i)%field_type == FIELD_TYPE_CELL) then
                ! Write field data based on number of components
                select case (fields(i)%ncomp)
                    case (1)
                        ! Scalar field
                        write(unit_num, '(A,1X,A,1X,A)') 'SCALARS', trim(fields(i)%name), 'float'
                        write(unit_num, '(A)') 'LOOKUP_TABLE default'
                        
                        ! Write scalar data
                        do k = 1, mesh%nz
                            do j = 1, mesh%ny
                                do i = 1, mesh%nx
                                    idx = i + (j-1)*mesh%nx + (k-1)*mesh%nx*mesh%ny
                                    write(unit_num, '(F12.6)') fields(i)%data(idx)
                                end do
                            end do
                        end do
                        
                    case (3)
                        ! Vector field
                        write(unit_num, '(A,1X,A,1X,A)') 'VECTORS', trim(fields(i)%name), 'float'
                        
                        ! Write vector data
                        do k = 1, mesh%nz
                            do j = 1, mesh%ny
                                do i = 1, mesh%nx
                                    idx = i + (j-1)*mesh%nx + (k-1)*mesh%nx*mesh%ny
                                    write(unit_num, '(3(F12.6,1X))') &
                                        fields(i)%data((idx-1)*3+1), &
                                        fields(i)%data((idx-1)*3+2), &
                                        fields(i)%data((idx-1)*3+3)
                                end do
                            end do
                        end do
                        
                    case default
                        print *, "Warning: Field '", trim(fields(i)%name), "' has unsupported number of components: ", &
                                  fields(i)%ncomp, ". Skipping."
                end select
            end if
        end do
        
        ! Close file
        close(unit_num)
    end subroutine vtk_converter_write_vtk_file
    
    !> Write VTK XML UnstructuredGrid format file (.vtu)
    subroutine vtk_converter_write_vtu_file(this, mesh, fields, vtu_file, iteration, time)
        class(vtk_converter_t), intent(inout) :: this
        type(mesh_t), intent(in) :: mesh
        type(field_t), intent(in) :: fields(:)
        character(len=*), intent(in) :: vtu_file
        integer, intent(in) :: iteration
        real(WP), intent(in) :: time
        
        integer :: unit_num, ios, i, j, k, idx, cell_idx
        
        ! Open VTU file for writing
        open(newunit=unit_num, file=trim(vtu_file), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            print *, "Error: Could not open VTU file for writing: ", trim(vtu_file)
            return
        end if
        
        ! Write XML header
        write(unit_num, '(A)') '<?xml version="1.0"?>'
        write(unit_num, '(A)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
        write(unit_num, '(A)') '  <UnstructuredGrid>'
        
        ! Write piece information
        write(unit_num, '(A,I0,A,I0,A)') '    <Piece NumberOfPoints="', &
                                         (mesh%nx+1) * (mesh%ny+1) * (mesh%nz+1), &
                                         '" NumberOfCells="', &
                                         mesh%nx * mesh%ny * mesh%nz, '">'
        
        ! Write point data
        write(unit_num, '(A)') '      <Points>'
        write(unit_num, '(A)') '        <DataArray type="Float32" NumberOfComponents="3" format="ascii">'
        
        ! Write vertex coordinates
        do k = 0, mesh%nz
            do j = 0, mesh%ny
                do i = 0, mesh%nx
                    ! Calculate vertex coordinates based on domain dimensions
                    write(unit_num, '(3(F12.6,1X))') &
                        mesh%xmin + i * (mesh%xmax - mesh%xmin) / mesh%nx, &
                        mesh%ymin + j * (mesh%ymax - mesh%ymin) / mesh%ny, &
                        mesh%zmin + k * (mesh%zmax - mesh%zmin) / mesh%nz
                end do
            end do
        end do
        
        write(unit_num, '(A)') '        </DataArray>'
        write(unit_num, '(A)') '      </Points>'
        
        ! Write cells (connectivity, offsets, types)
        write(unit_num, '(A)') '      <Cells>'
        
        ! Connectivity
        write(unit_num, '(A)') '        <DataArray type="Int32" Name="connectivity" format="ascii">'
        
        ! Write hexahedron connectivity for each cell
        do k = 1, mesh%nz
            do j = 1, mesh%ny
                do i = 1, mesh%nx
                    ! Calculate vertex indices for this cell
                    ! Vertices are ordered in VTK as:
                    ! 0: (i,j,k), 1: (i+1,j,k), 2: (i+1,j+1,k), 3: (i,j+1,k)
                    ! 4: (i,j,k+1), 5: (i+1,j,k+1), 6: (i+1,j+1,k+1), 7: (i,j+1,k+1)
                    idx = (i-1) + (j-1)*(mesh%nx+1) + (k-1)*(mesh%nx+1)*(mesh%ny+1)
                    
                    write(unit_num, '(8(I0,1X))') &
                        idx, idx+1, idx+mesh%nx+2, idx+mesh%nx+1, &
                        idx+(mesh%nx+1)*(mesh%ny+1), idx+1+(mesh%nx+1)*(mesh%ny+1), &
                        idx+mesh%nx+2+(mesh%nx+1)*(mesh%ny+1), idx+mesh%nx+1+(mesh%nx+1)*(mesh%ny+1)
                end do
            end do
        end do
        
        write(unit_num, '(A)') '        </DataArray>'
        
        ! Offsets
        write(unit_num, '(A)') '        <DataArray type="Int32" Name="offsets" format="ascii">'
        
        ! Write offset for each cell (8 vertices per cell)
        do i = 1, mesh%nx * mesh%ny * mesh%nz
            write(unit_num, '(I0)') i * 8
        end do
        
        write(unit_num, '(A)') '        </DataArray>'
        
        ! Cell types
        write(unit_num, '(A)') '        <DataArray type="UInt8" Name="types" format="ascii">'
        
        ! Write cell type for each cell (12 = VTK_HEXAHEDRON)
        do i = 1, mesh%nx * mesh%ny * mesh%nz
            write(unit_num, '(I0)') 12
        end do
        
        write(unit_num, '(A)') '        </DataArray>'
        write(unit_num, '(A)') '      </Cells>'
        
        ! Write cell data
        write(unit_num, '(A)') '      <CellData>'
        
        ! Add iteration and time as field data
        write(unit_num, '(A)') '        <DataArray type="Int32" Name="Iteration" NumberOfComponents="1" format="ascii">'
        do i = 1, mesh%nx * mesh%ny * mesh%nz
            write(unit_num, '(I0)') iteration
        end do
        write(unit_num, '(A)') '        </DataArray>'
        
        write(unit_num, '(A)') '        <DataArray type="Float32" Name="Time" NumberOfComponents="1" format="ascii">'
        do i = 1, mesh%nx * mesh%ny * mesh%nz
            write(unit_num, '(F12.6)') time
        end do
        write(unit_num, '(A)') '        </DataArray>'
        
        ! Write fields
        do i = 1, size(fields)
            ! Only write cell-centered fields
            if (fields(i)%field_type == FIELD_TYPE_CELL) then
                ! Write field data based on number of components
                select case (fields(i)%ncomp)
                    case (1)
                        ! Scalar field
                        write(unit_num, '(A,A,A)') '        <DataArray type="Float32" Name="', &
                                                  trim(fields(i)%name), &
                                                  '" NumberOfComponents="1" format="ascii">'
                        
                        ! Write scalar data
                        do k = 1, mesh%nz
                            do j = 1, mesh%ny
                                do i = 1, mesh%nx
                                    cell_idx = i + (j-1)*mesh%nx + (k-1)*mesh%nx*mesh%ny
                                    write(unit_num, '(F12.6)') fields(i)%data(cell_idx)
                                end do
                            end do
                        end do
                        
                    case (3)
                        ! Vector field
                        write(unit_num, '(A,A,A)') '        <DataArray type="Float32" Name="', &
                                                  trim(fields(i)%name), &
                                                  '" NumberOfComponents="3" format="ascii">'
                        
                        ! Write vector data
                        do k = 1, mesh%nz
                            do j = 1, mesh%ny
                                do i = 1, mesh%nx
                                    cell_idx = i + (j-1)*mesh%nx + (k-1)*mesh%nx*mesh%ny
                                    write(unit_num, '(3(F12.6,1X))') &
                                        fields(i)%data((cell_idx-1)*3+1), &
                                        fields(i)%data((cell_idx-1)*3+2), &
                                        fields(i)%data((cell_idx-1)*3+3)
                                end do
                            end do
                        end do
                        
                    case default
                        print *, "Warning: Field '", trim(fields(i)%name), "' has unsupported number of components: ", &
                                  fields(i)%ncomp, ". Skipping."
                end select
                
                write(unit_num, '(A)') '        </DataArray>'
            end if
        end do
        
        write(unit_num, '(A)') '      </CellData>'
        write(unit_num, '(A)') '    </Piece>'
        write(unit_num, '(A)') '  </UnstructuredGrid>'
        write(unit_num, '(A)') '</VTKFile>'
        
        ! Close file
        close(unit_num)
    end subroutine vtk_converter_write_vtu_file

end module vtk_converter 