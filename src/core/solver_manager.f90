module solver_manager
    use omp_lib
    use mesh, only: mesh_t, WP
    use field, only: field_t, FIELD_TYPE_CELL
    use boundary, only: boundary_manager_t
    use solution_io, only: solution_io_t
    use config_parser, only: config_t, config_section_t
    
    ! Import solver modules - use preprocessor to make these dependencies optional
    ! This allows for more flexible compilation order
#ifdef HAVE_DIFFUSION
    use diffusion, only: diffusion_solver_t
#endif
#ifdef HAVE_INCOMPRESSIBLE_FLOW
    use incompressible_flow, only: incompressible_flow_solver_t
#endif
    implicit none
    
    ! Constants for solver types
    integer, parameter :: SOLVER_TYPE_DIFFUSION = 1
    integer, parameter :: SOLVER_TYPE_INCOMPRESSIBLE_FLOW = 2
    
    !> Solver manager data structure
    type :: solver_manager_t
        ! Mesh
        type(mesh_t), pointer :: mesh => null()
        
        ! Boundary condition manager
        type(boundary_manager_t), pointer :: bc_manager => null()
        
        ! Fields storage
        integer :: num_fields = 0
        type(field_t), allocatable :: fields(:)
        
        ! Solver type and instances
        integer :: solver_type = 0
        
#ifdef HAVE_DIFFUSION
        type(diffusion_solver_t) :: diffusion_solver
#endif
#ifdef HAVE_INCOMPRESSIBLE_FLOW
        type(incompressible_flow_solver_t) :: flow_solver
#endif
        
        ! Settings
        real(WP) :: dt = 0.001_WP
        real(WP) :: end_time = 1.0_WP
        integer :: max_iterations = 1000
        
        ! I/O control
        logical :: do_write_solution = .true.
        integer :: output_frequency = 10
        character(len=256) :: output_dir = "results"
        character(len=64) :: case_name = "simulation"
        
        ! Simulation state
        real(WP) :: current_time = 0.0_WP
        integer :: current_iteration = 0
        logical :: simulation_finished = .false.
        
        ! GPU acceleration
        logical :: use_gpu = .false.
        
    contains
        procedure :: init => solver_manager_init
        procedure :: setup_from_config => solver_manager_setup_from_config
        procedure :: create_mesh => solver_manager_create_mesh
        procedure :: create_fields => solver_manager_create_fields
        procedure :: setup_boundary_conditions => solver_manager_setup_boundary_conditions
        procedure :: setup_solver => solver_manager_setup_solver
        procedure :: advance => solver_manager_advance
        procedure :: write_solution => solver_manager_write_solution
        procedure :: cleanup => solver_manager_cleanup
    end type solver_manager_t
    
contains

    !> Initialize the solver manager
    subroutine solver_manager_init(this, config_file)
        class(solver_manager_t), intent(inout) :: this
        character(len=*), intent(in) :: config_file
        
        ! Check if GPU is available
        this%use_gpu = (omp_get_num_devices() > 0)
        
        ! Setup the solver manager from configuration file
        call this%setup_from_config(config_file)
    end subroutine solver_manager_init
    
    !> Setup solver manager from configuration file
    subroutine solver_manager_setup_from_config(this, config_file)
        class(solver_manager_t), intent(inout) :: this
        character(len=*), intent(in) :: config_file
        character(len=256) :: solver_type_str
        type(config_t) :: config
        
        print *, "Setting up solver from configuration file: ", trim(config_file)
        
        ! Load configuration
        call config%load(config_file)
        
        ! Create mesh based on configuration
        call this%create_mesh(config_file)
        
        ! Determine solver type from configuration
        ! Check if solver type is specified in config
        if (config%has_key("solver.type")) then
            solver_type_str = config%get_string("solver.type", "")
            
            if (trim(solver_type_str) == "diffusion") then
                this%solver_type = SOLVER_TYPE_DIFFUSION
            else if (trim(solver_type_str) == "incompressible_flow") then
                this%solver_type = SOLVER_TYPE_INCOMPRESSIBLE_FLOW
            else
                ! Default to diffusion if unknown type
                print *, "Warning: Unknown solver type '", trim(solver_type_str), "'. Defaulting to diffusion."
                solver_type_str = "diffusion"
                this%solver_type = SOLVER_TYPE_DIFFUSION
            end if
        end if
        
        print *, "Solver type: ", trim(solver_type_str)
        
        ! Create fields based on configuration
        call this%create_fields(config_file)
        
        ! Setup boundary conditions based on configuration
        call this%setup_boundary_conditions(config_file)
        
        ! Setup solver based on type
        call this%setup_solver(config_file)
        
        ! Clean up config
        call config%cleanup()
    end subroutine solver_manager_setup_from_config
    
    !> Create mesh based on configuration
    subroutine solver_manager_create_mesh(this, config_file)
        class(solver_manager_t), intent(inout) :: this
        character(len=*), intent(in) :: config_file
        integer :: nx, ny, nz
        real(WP) :: xmin, xmax, ymin, ymax, zmin, zmax
        type(config_t) :: config
        integer, allocatable :: mesh_dimensions(:)
        real(WP), allocatable :: domain_bounds(:)
        
        ! Load configuration
        call config%load(config_file)
        
        ! Read mesh dimensions from config
        if (config%has_key("mesh.dimensions")) then
            mesh_dimensions = config%get_integer_array("mesh.dimensions")
            if (size(mesh_dimensions) >= 3) then
                nx = mesh_dimensions(1)
                ny = mesh_dimensions(2)
                nz = mesh_dimensions(3)
            else
                print *, "Warning: mesh.dimensions array should have 3 values. Using defaults."
            end if
        else
            print *, "Using default mesh dimensions: ", nx, ny, nz
        end if
        
        ! Read domain bounds from config
        if (config%has_key("mesh.domain")) then
            domain_bounds = config%get_real_array("mesh.domain")
            if (size(domain_bounds) >= 6) then
                xmin = domain_bounds(1)
                xmax = domain_bounds(2)
                ymin = domain_bounds(3)
                ymax = domain_bounds(4)
                zmin = domain_bounds(5)
                zmax = domain_bounds(6)
            else
                print *, "Warning: mesh.domain array should have 6 values. Using defaults."
            end if
        else
            print *, "Using default domain bounds: ", xmin, xmax, ymin, ymax, zmin, zmax
        end if
        
        ! Allocate mesh
        allocate(this%mesh)
        
        ! Initialize mesh
        call this%mesh%init()
        
        ! Create cartesian mesh
        call this%mesh%create_cartesian(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax)
        
        ! Optimize for GPU if available
        if (this%use_gpu) then
            call this%mesh%optimize_for_gpu(8, 8, 8)
        end if
        
        print *, "Created mesh with ", this%mesh%topo%ncells, " cells"
        
        ! Clean up
        call config%cleanup()
    end subroutine solver_manager_create_mesh
    
    !> Create fields based on configuration
    subroutine solver_manager_create_fields(this, config_file)
        class(solver_manager_t), intent(inout) :: this
        character(len=*), intent(in) :: config_file
        integer :: i
        
        ! Different fields are created depending on the solver type
        select case (this%solver_type)
            case (SOLVER_TYPE_DIFFUSION)
                ! For diffusion solver, we need temperature, diffusivity, and source
                this%num_fields = 3
                allocate(this%fields(this%num_fields))
                
                ! Initialize fields
                call this%fields(1)%init(this%mesh, "temperature", FIELD_TYPE_CELL, 1)
                call this%fields(2)%init(this%mesh, "diffusivity", FIELD_TYPE_CELL, 1)
                call this%fields(3)%init(this%mesh, "source", FIELD_TYPE_CELL, 1)
                
                ! Set initial values
                call this%fields(1)%set_uniform([0.0_WP])  ! Initial temperature
                call this%fields(2)%set_uniform([1.0_WP])  ! Uniform diffusivity
                call this%fields(3)%set_uniform([0.0_WP])  ! Zero source term
                
            case (SOLVER_TYPE_INCOMPRESSIBLE_FLOW)
                ! For incompressible flow solver, we need velocity components, pressure, and viscosity
                this%num_fields = 5
                allocate(this%fields(this%num_fields))
                
                ! Initialize fields
                call this%fields(1)%init(this%mesh, "u", FIELD_TYPE_CELL, 1)      ! X velocity
                call this%fields(2)%init(this%mesh, "v", FIELD_TYPE_CELL, 1)      ! Y velocity
                call this%fields(3)%init(this%mesh, "w", FIELD_TYPE_CELL, 1)      ! Z velocity
                call this%fields(4)%init(this%mesh, "p", FIELD_TYPE_CELL, 1)      ! Pressure
                call this%fields(5)%init(this%mesh, "viscosity", FIELD_TYPE_CELL, 1)  ! Viscosity
                
                ! Set initial values
                call this%fields(1)%set_uniform([0.0_WP])  ! Initial u velocity
                call this%fields(2)%set_uniform([0.0_WP])  ! Initial v velocity
                call this%fields(3)%set_uniform([0.0_WP])  ! Initial w velocity
                call this%fields(4)%set_uniform([0.0_WP])  ! Initial pressure
                call this%fields(5)%set_uniform([0.01_WP]) ! Default viscosity (Re=100 for 1m domain with U=1)
            
            case default
                print *, "Error: Unknown solver type. Cannot create fields."
        end select
        
        ! Copy fields to device if using GPU
        if (this%use_gpu) then
            do i = 1, this%num_fields
                call this%fields(i)%copy_to_device()
            end do
        end if
    end subroutine solver_manager_create_fields
    
    !> Setup boundary conditions based on configuration
    subroutine solver_manager_setup_boundary_conditions(this, config_file)
        class(solver_manager_t), intent(inout) :: this
        character(len=*), intent(in) :: config_file
        integer :: BC_DIRICHLET, BC_NEUMANN, BC_WALL, BC_SYMMETRY
        type(config_t) :: config
        type(config_section_t), pointer :: bc_section => null()
        character(len=64) :: wall_name, bc_type
        character(len=256) :: velocity_str
        real(WP) :: velocity(3)
        integer :: i, j
        
        ! Get boundary condition types from boundary module
        BC_DIRICHLET = 1
        BC_NEUMANN = 2
        BC_SYMMETRY = 3
        BC_WALL = 5
        
        ! Allocate boundary condition manager
        allocate(this%bc_manager)
        
        ! Initialize boundary condition manager
        call this%bc_manager%init(this%mesh)
        
        ! Load configuration
        call config%load(config_file)
        
        ! Try to read boundary conditions from config file
        if (config%has_section("boundary_conditions")) then
            print *, "Setting up boundary conditions from config file..."
            
            ! Get boundary conditions section
            call config%get_section("boundary_conditions", bc_section)
            
            if (associated(bc_section)) then
                print *, "Found boundary conditions in config file:"
                
                ! Loop through the walls (north, south, east, west, etc.)
                do i = 1, bc_section%num_subsections
                    print *, "Processing boundary condition section: ", i
                    wall_name = bc_section%subsections(i)%name
                    bc_type = ""
                    velocity = [0.0_WP, 0.0_WP, 0.0_WP]
                    
                    ! Get properties for this wall
                    do j = 1, bc_section%subsections(i)%num_items
                        if (trim(bc_section%subsections(i)%keys(j)) == "type") then
                            bc_type = trim(bc_section%subsections(i)%values(j))
                        end if
                        
                        if (trim(bc_section%subsections(i)%keys(j)) == "velocity") then
                            velocity_str = trim(bc_section%subsections(i)%values(j))
                            ! Parse velocity string like [1.0, 0.0, 0.0]
                            read(velocity_str(2:len_trim(velocity_str)-1), *) velocity
                        end if
                    end do
                    
                    print *, "Wall: ", trim(wall_name), " Type: ", trim(bc_type), " Velocity: ", velocity
                    
                    ! Map wall name to index
                    if (trim(wall_name) == "north") then
                        if (trim(bc_type) == "wall") then
                            call this%bc_manager%add_bc("u", 4, BC_WALL, [velocity(1)])
                            call this%bc_manager%add_bc("v", 4, BC_WALL, [velocity(2)])
                            call this%bc_manager%add_bc("w", 4, BC_WALL, [velocity(3)])
                            call this%bc_manager%add_bc("p", 4, BC_NEUMANN, [0.0_WP])
                        else if (trim(bc_type) == "symmetry") then
                            call this%bc_manager%add_bc("u", 4, BC_SYMMETRY, [0.0_WP])
                            call this%bc_manager%add_bc("v", 4, BC_SYMMETRY, [0.0_WP])
                            call this%bc_manager%add_bc("w", 4, BC_SYMMETRY, [0.0_WP])
                            call this%bc_manager%add_bc("p", 4, BC_SYMMETRY, [0.0_WP])
                        end if
                    else if (trim(wall_name) == "south") then
                        if (trim(bc_type) == "wall") then
                            call this%bc_manager%add_bc("u", 3, BC_WALL, [velocity(1)])
                            call this%bc_manager%add_bc("v", 3, BC_WALL, [velocity(2)])
                            call this%bc_manager%add_bc("w", 3, BC_WALL, [velocity(3)])
                            call this%bc_manager%add_bc("p", 3, BC_NEUMANN, [0.0_WP])
                        else if (trim(bc_type) == "symmetry") then
                            call this%bc_manager%add_bc("u", 3, BC_SYMMETRY, [0.0_WP])
                            call this%bc_manager%add_bc("v", 3, BC_SYMMETRY, [0.0_WP])
                            call this%bc_manager%add_bc("w", 3, BC_SYMMETRY, [0.0_WP])
                            call this%bc_manager%add_bc("p", 3, BC_SYMMETRY, [0.0_WP])
                        end if
                    else if (trim(wall_name) == "east") then
                        if (trim(bc_type) == "wall") then
                            call this%bc_manager%add_bc("u", 2, BC_WALL, [velocity(1)])
                            call this%bc_manager%add_bc("v", 2, BC_WALL, [velocity(2)])
                            call this%bc_manager%add_bc("w", 2, BC_WALL, [velocity(3)])
                            call this%bc_manager%add_bc("p", 2, BC_NEUMANN, [0.0_WP])
                        else if (trim(bc_type) == "symmetry") then
                            call this%bc_manager%add_bc("u", 2, BC_SYMMETRY, [0.0_WP])
                            call this%bc_manager%add_bc("v", 2, BC_SYMMETRY, [0.0_WP])
                            call this%bc_manager%add_bc("w", 2, BC_SYMMETRY, [0.0_WP])
                            call this%bc_manager%add_bc("p", 2, BC_SYMMETRY, [0.0_WP])
                        end if
                    else if (trim(wall_name) == "west") then
                        if (trim(bc_type) == "wall") then
                            call this%bc_manager%add_bc("u", 1, BC_WALL, [velocity(1)])
                            call this%bc_manager%add_bc("v", 1, BC_WALL, [velocity(2)])
                            call this%bc_manager%add_bc("w", 1, BC_WALL, [velocity(3)])
                            call this%bc_manager%add_bc("p", 1, BC_NEUMANN, [0.0_WP])
                        else if (trim(bc_type) == "symmetry") then
                            call this%bc_manager%add_bc("u", 1, BC_SYMMETRY, [0.0_WP])
                            call this%bc_manager%add_bc("v", 1, BC_SYMMETRY, [0.0_WP])
                            call this%bc_manager%add_bc("w", 1, BC_SYMMETRY, [0.0_WP])
                            call this%bc_manager%add_bc("p", 1, BC_SYMMETRY, [0.0_WP])
                        end if
                    else if (trim(wall_name) == "top") then
                        if (trim(bc_type) == "wall") then
                            call this%bc_manager%add_bc("u", 6, BC_WALL, [velocity(1)])
                            call this%bc_manager%add_bc("v", 6, BC_WALL, [velocity(2)])
                            call this%bc_manager%add_bc("w", 6, BC_WALL, [velocity(3)])
                            call this%bc_manager%add_bc("p", 6, BC_NEUMANN, [0.0_WP])
                        else if (trim(bc_type) == "symmetry") then
                            call this%bc_manager%add_bc("u", 6, BC_SYMMETRY, [0.0_WP])
                            call this%bc_manager%add_bc("v", 6, BC_SYMMETRY, [0.0_WP])
                            call this%bc_manager%add_bc("w", 6, BC_SYMMETRY, [0.0_WP])
                            call this%bc_manager%add_bc("p", 6, BC_SYMMETRY, [0.0_WP])
                        end if
                    else if (trim(wall_name) == "bottom") then
                        if (trim(bc_type) == "wall") then
                            call this%bc_manager%add_bc("u", 5, BC_WALL, [velocity(1)])
                            call this%bc_manager%add_bc("v", 5, BC_WALL, [velocity(2)])
                            call this%bc_manager%add_bc("w", 5, BC_WALL, [velocity(3)])
                            call this%bc_manager%add_bc("p", 5, BC_NEUMANN, [0.0_WP])
                        else if (trim(bc_type) == "symmetry") then
                            call this%bc_manager%add_bc("u", 5, BC_SYMMETRY, [0.0_WP])
                            call this%bc_manager%add_bc("v", 5, BC_SYMMETRY, [0.0_WP])
                            call this%bc_manager%add_bc("w", 5, BC_SYMMETRY, [0.0_WP])
                            call this%bc_manager%add_bc("p", 5, BC_SYMMETRY, [0.0_WP])
                        end if
                    end if
                end do
                
                print *, "Boundary conditions from config file applied successfully"
                ! Clean up
                call config%cleanup()
                return
            end if
        end if
        
        ! If we get here, use default boundary conditions
        print *, "Warning: Using default boundary conditions - config file boundary conditions not applied"
        
        ! Setup boundary conditions based on solver type
        select case (this%solver_type)
            case (SOLVER_TYPE_DIFFUSION)
                ! For diffusion solver, we set up temperature boundary conditions
                call this%bc_manager%add_bc("temperature", 1, BC_DIRICHLET, [1.0_WP])  ! West = hot
                call this%bc_manager%add_bc("temperature", 2, BC_DIRICHLET, [0.0_WP])  ! East = cold
                call this%bc_manager%add_bc("temperature", 3, BC_NEUMANN, [0.0_WP])    ! South = insulated
                call this%bc_manager%add_bc("temperature", 4, BC_NEUMANN, [0.0_WP])    ! North = insulated
                call this%bc_manager%add_bc("temperature", 5, BC_NEUMANN, [0.0_WP])    ! Bottom = insulated
                call this%bc_manager%add_bc("temperature", 6, BC_NEUMANN, [0.0_WP])    ! Top = insulated
                
            case (SOLVER_TYPE_INCOMPRESSIBLE_FLOW)
                ! For incompressible flow solver, we set up velocity and pressure boundary conditions
                ! This is for a lid-driven cavity case
                
                ! Velocities
                call this%bc_manager%add_bc("u", 1, BC_WALL, [0.0_WP])  ! West wall (no-slip)
                call this%bc_manager%add_bc("u", 2, BC_WALL, [0.0_WP])  ! East wall (no-slip)
                call this%bc_manager%add_bc("u", 3, BC_WALL, [0.0_WP])  ! South wall (no-slip)
                call this%bc_manager%add_bc("u", 4, BC_WALL, [1.0_WP])  ! North wall (moving lid)
                call this%bc_manager%add_bc("u", 5, BC_SYMMETRY, [0.0_WP])  ! Bottom (symmetry)
                call this%bc_manager%add_bc("u", 6, BC_SYMMETRY, [0.0_WP])  ! Top (symmetry)
                
                call this%bc_manager%add_bc("v", 1, BC_WALL, [0.0_WP])  ! West wall (no-slip)
                call this%bc_manager%add_bc("v", 2, BC_WALL, [0.0_WP])  ! East wall (no-slip)
                call this%bc_manager%add_bc("v", 3, BC_WALL, [0.0_WP])  ! South wall (no-slip)
                call this%bc_manager%add_bc("v", 4, BC_WALL, [0.0_WP])  ! North wall (no-slip)
                call this%bc_manager%add_bc("v", 5, BC_SYMMETRY, [0.0_WP])  ! Bottom (symmetry)
                call this%bc_manager%add_bc("v", 6, BC_SYMMETRY, [0.0_WP])  ! Top (symmetry)
                
                call this%bc_manager%add_bc("w", 1, BC_WALL, [0.0_WP])  ! West wall (no-slip)
                call this%bc_manager%add_bc("w", 2, BC_WALL, [0.0_WP])  ! East wall (no-slip)
                call this%bc_manager%add_bc("w", 3, BC_WALL, [0.0_WP])  ! South wall (no-slip)
                call this%bc_manager%add_bc("w", 4, BC_WALL, [0.0_WP])  ! North wall (no-slip)
                call this%bc_manager%add_bc("w", 5, BC_SYMMETRY, [0.0_WP])  ! Bottom (symmetry)
                call this%bc_manager%add_bc("w", 6, BC_SYMMETRY, [0.0_WP])  ! Top (symmetry)
                
                ! Pressure
                call this%bc_manager%add_bc("p", 1, BC_NEUMANN, [0.0_WP])  ! West wall (zero gradient)
                call this%bc_manager%add_bc("p", 2, BC_NEUMANN, [0.0_WP])  ! East wall (zero gradient)
                call this%bc_manager%add_bc("p", 3, BC_NEUMANN, [0.0_WP])  ! South wall (zero gradient)
                call this%bc_manager%add_bc("p", 4, BC_NEUMANN, [0.0_WP])  ! North wall (zero gradient)
                call this%bc_manager%add_bc("p", 5, BC_SYMMETRY, [0.0_WP])  ! Bottom (symmetry)
                call this%bc_manager%add_bc("p", 6, BC_SYMMETRY, [0.0_WP])  ! Top (symmetry)
                
            case default
                print *, "Error: Unknown solver type. Cannot set up boundary conditions."
        end select
        
        print *, "Boundary conditions set up for solver"
    end subroutine solver_manager_setup_boundary_conditions
    
    !> Setup solver based on configuration
    subroutine solver_manager_setup_solver(this, config_file)
        class(solver_manager_t), intent(inout) :: this
        character(len=*), intent(in) :: config_file
        type(config_t) :: config
        integer :: i, j  ! Loop counters for boundary conditions
        type(config_section_t), pointer :: bc_section => null()
        
        ! Load configuration
        call config%load(config_file)
        
        ! Set solver parameters
        this%dt = 0.001_WP  ! Default time step
        this%end_time = 1.0_WP  ! Default end time
        this%max_iterations = 1000  ! Default max iterations
        
        ! Override with values from config
        if (config%has_section("solver")) then
            
            ! Time step
            if (config%has_key("solver.time_step")) then
                this%dt = config%get_real("solver.time_step", this%dt)
            end if
            
            ! End time
            if (config%has_key("solver.end_time")) then
                this%end_time = config%get_real("solver.end_time", this%end_time)
            end if
            
            ! Solver type
            if (config%has_key("solver.type")) then
                this%solver_type = config%get_integer("solver.type", this%solver_type)
            end if
            
            ! Max iterations (if present)
            if (config%has_key("solver.max_iterations")) then
                this%max_iterations = config%get_integer("solver.max_iterations", this%max_iterations)
            end if
        else
            print *, "No solver section found in config file, using defaults or filename heuristics"
            ! Default values
            this%dt = 0.001_WP
            this%end_time = 1.0_WP
            this%max_iterations = 1000
            this%solver_type = SOLVER_TYPE_INCOMPRESSIBLE_FLOW

            print *, "Using default solver settings"
            print *, "  * Time step: ", this%dt
            print *, "  * End time: ", this%end_time
            print *, "  * Max iterations: ", this%max_iterations
            print *, "  * Solver type: ", this%solver_type
        end if

        
        ! Read output settings from config file
        if (config%has_section("output")) then
            print *, "DEBUG: Found output section in config file"
            
            ! Read output directory
            if (config%has_key("output.directory")) then
                this%output_dir = config%get_string("output.directory", this%output_dir)
            end if
            
            ! Read case name
            if (config%has_key("output.case_name")) then
                this%case_name = config%get_string("output.case_name", this%case_name)
            end if
            
            ! Read output frequency
            if (config%has_key("output.frequency")) then
                this%output_frequency = config%get_integer("output.frequency", this%output_frequency)
            end if
        else
            print *, "No output section found in config file, using defaults"
        end if
        
        ! Recalculate max_iterations based on end_time and dt
        this%max_iterations = ceiling(this%end_time / this%dt)
        print *, "Final max_iterations (based on end_time/dt): ", this%max_iterations
        
        ! Initialize the appropriate solver based on type
        select case (this%solver_type)
            case (SOLVER_TYPE_DIFFUSION)
#ifdef HAVE_DIFFUSION
                ! Initialize diffusion solver
                call this%diffusion_solver%init(this%mesh, this%dt, this%max_iterations, 1.0e-6_WP)
                
                ! Setup fields for diffusion solver
                call this%diffusion_solver%setup_fields(this%fields(1), this%fields(2), this%fields(3))
#else
                print *, "Error: Diffusion solver not available. Recompile with diffusion support."
                stop
#endif
                
            case (SOLVER_TYPE_INCOMPRESSIBLE_FLOW)
#ifdef HAVE_INCOMPRESSIBLE_FLOW
                ! Initialize incompressible flow solver
                call this%flow_solver%init(this%mesh, this%dt, this%end_time, this%max_iterations, &
                                         1.0e-6_WP, 1.0e-5_WP, "SIMPLE")
                
                ! Setup fields for incompressible flow solver
                call this%flow_solver%setup_fields(this%fields(1), this%fields(2), this%fields(3), &
                                                this%fields(4), this%fields(5))
                
                ! Set boundary condition manager
                call this%flow_solver%setup_boundary_manager(this%bc_manager)
#else
                print *, "Error: Incompressible flow solver not available. Recompile with flow solver support."
                stop
#endif
                
            case default
                print *, "Error: Unknown solver type. Cannot set up solver."
        end select
        
        ! Clean up
        call config%cleanup()
    end subroutine solver_manager_setup_solver
    
    !> Advance solution by one time step
    subroutine solver_manager_advance(this)
        class(solver_manager_t), intent(inout) :: this
        
        ! Advance solution based on solver type
        select case (this%solver_type)
            case (SOLVER_TYPE_DIFFUSION)
#ifdef HAVE_DIFFUSION
                call this%diffusion_solver%advance()
#else
                print *, "Error: Diffusion solver not available."
                stop
#endif
                
            case (SOLVER_TYPE_INCOMPRESSIBLE_FLOW)
#ifdef HAVE_INCOMPRESSIBLE_FLOW
                call this%flow_solver%advance()
#else
                print *, "Error: Incompressible flow solver not available."
                stop
#endif
                
            case default
                print *, "Error: Unknown solver type. Cannot advance solution."
        end select
        
        ! Update simulation state
        this%current_iteration = this%current_iteration + 1
        this%current_time = this%current_time + this%dt
        
        ! Check if simulation is finished
        if (this%current_time >= this%end_time .or. this%current_iteration >= this%max_iterations) then
            this%simulation_finished = .true.
        end if
    end subroutine solver_manager_advance
    
    !> Write solution to file
    subroutine solver_manager_write_solution(this, solution_io)
        class(solver_manager_t), intent(inout) :: this
        type(solution_io_t), intent(inout) :: solution_io
        integer :: i
        
        ! Keep track of the original solution_io settings
        character(len=256) :: orig_output_dir
        character(len=64) :: orig_prefix
        integer :: orig_frequency
        
        ! Save original values
        orig_output_dir = solution_io%output_dir
        orig_prefix = solution_io%prefix
        orig_frequency = solution_io%output_frequency
        
        ! Copy fields from device if using GPU
        if (this%use_gpu) then
            do i = 1, this%num_fields
                call this%fields(i)%copy_from_device()
            end do
        end if
        
        call solution_io%write_solution(this%mesh, this%fields, &
                                      this%current_iteration, this%current_time)
    end subroutine solver_manager_write_solution
    
    !> Cleanup
    subroutine solver_manager_cleanup(this)
        class(solver_manager_t), intent(inout) :: this
        integer :: i
        
        ! Cleanup fields
        do i = 1, this%num_fields
            call this%fields(i)%cleanup()
        end do
        
        ! Deallocate fields array
        if (allocated(this%fields)) then
            deallocate(this%fields)
        end if
        
        ! Cleanup solver based on type
        select case (this%solver_type)
            case (SOLVER_TYPE_DIFFUSION)
#ifdef HAVE_DIFFUSION
                call this%diffusion_solver%cleanup()
#endif
                
            case (SOLVER_TYPE_INCOMPRESSIBLE_FLOW)
#ifdef HAVE_INCOMPRESSIBLE_FLOW
                call this%flow_solver%cleanup()
#endif
                
            case default
                ! Nothing to do
        end select
        
        ! Cleanup boundary condition manager
        if (associated(this%bc_manager)) then
            call this%bc_manager%cleanup()
            deallocate(this%bc_manager)
            this%bc_manager => null()
        end if
        
        ! Cleanup mesh
        if (associated(this%mesh)) then
            call this%mesh%cleanup()
            deallocate(this%mesh)
            this%mesh => null()
        end if
    end subroutine solver_manager_cleanup

end module solver_manager 