module incompressible_flow
    use omp_lib
    use mesh, only: mesh_t, WP
    use field, only: field_t, FIELD_TYPE_CELL, FIELD_TYPE_FACE
    use boundary, only: boundary_manager_t
    implicit none
    
    ! Constants for the pressure-velocity coupling algorithm
    integer, parameter :: SIMPLE = 1
    integer, parameter :: PISO = 2
    integer, parameter :: PIMPLE = 3
    
    !> Incompressible flow solver type
    type :: incompressible_flow_solver_t
        ! Fields
        type(field_t), pointer :: u => null()           ! X-velocity
        type(field_t), pointer :: v => null()           ! Y-velocity
        type(field_t), pointer :: w => null()           ! Z-velocity
        type(field_t), pointer :: p => null()           ! Pressure
        type(field_t), pointer :: viscosity => null()   ! Kinematic viscosity
        
        ! Internal fields for computation
        type(field_t) :: u_star        ! Intermediate x-velocity
        type(field_t) :: v_star        ! Intermediate y-velocity
        type(field_t) :: w_star        ! Intermediate z-velocity
        type(field_t) :: p_prime       ! Pressure correction
        type(field_t) :: mass_flux     ! Mass flux at faces
        type(field_t) :: face_u        ! Face-interpolated x-velocity
        type(field_t) :: face_v        ! Face-interpolated y-velocity
        type(field_t) :: face_w        ! Face-interpolated z-velocity
        
        ! Boundary conditions
        type(boundary_manager_t), pointer :: bc_manager => null()
        
        ! Mesh
        type(mesh_t), pointer :: mesh => null()
        
        ! Solver settings
        real(WP) :: dt                         ! Time step
        real(WP) :: end_time                   ! End time
        real(WP) :: time                       ! Current time
        integer :: max_iterations              ! Maximum number of outer iterations
        integer :: max_inner_iterations        ! Maximum number of inner iterations
        real(WP) :: tolerance_velocity         ! Convergence tolerance for velocity
        real(WP) :: tolerance_pressure         ! Convergence tolerance for pressure
        integer :: pressure_correction_type    ! Pressure-velocity coupling algorithm
        
        ! Residual tracking
        real(WP) :: u_residual = 0.0_WP        ! Residual for x-momentum equation
        real(WP) :: v_residual = 0.0_WP        ! Residual for y-momentum equation
        real(WP) :: w_residual = 0.0_WP        ! Residual for z-momentum equation
        real(WP) :: p_residual = 0.0_WP        ! Residual for pressure equation
        
        ! Algorithm control
        real(WP) :: alpha_u = 0.7_WP           ! Velocity under-relaxation factor
        real(WP) :: alpha_p = 0.3_WP           ! Pressure under-relaxation factor
        
        ! GPU acceleration
        logical :: use_gpu = .false.           ! Will be set based on device availability
        
    contains
        procedure :: init => incompressible_flow_solver_init
        procedure :: setup_fields => incompressible_flow_solver_setup_fields
        procedure :: setup_boundary_manager => incompressible_flow_solver_setup_boundary_manager
        procedure :: advance => incompressible_flow_solver_advance
        procedure :: momentum_predictor => incompressible_flow_solver_momentum_predictor
        procedure :: pressure_correction => incompressible_flow_solver_pressure_correction
        procedure :: momentum_corrector => incompressible_flow_solver_momentum_corrector
        procedure :: update_mass_flux => incompressible_flow_solver_update_mass_flux
        procedure :: check_convergence => incompressible_flow_solver_check_convergence
        procedure :: apply_boundary_conditions => incompressible_flow_solver_apply_boundary_conditions
        procedure :: get_residuals => incompressible_flow_solver_get_residuals
        procedure :: cleanup => incompressible_flow_solver_cleanup
        procedure :: init_fields => incompressible_flow_solver_init_fields
    end type incompressible_flow_solver_t
    
contains

    !> Initialize the incompressible flow solver
    subroutine incompressible_flow_solver_init(this, mesh, dt, end_time, max_iter, &
                                               tol_vel, tol_p, p_corr_type)
        class(incompressible_flow_solver_t), intent(inout) :: this
        type(mesh_t), target, intent(in) :: mesh
        real(WP), intent(in) :: dt
        real(WP), intent(in) :: end_time
        integer, intent(in) :: max_iter
        real(WP), intent(in) :: tol_vel
        real(WP), intent(in) :: tol_p
        character(len=*), intent(in) :: p_corr_type
        
        ! Store mesh and parameters
        this%mesh => mesh
        this%dt = dt
        this%end_time = end_time
        this%max_iterations = max_iter
        this%tolerance_velocity = tol_vel
        this%tolerance_pressure = tol_p
        
        ! Set pressure-velocity coupling method
        select case (trim(p_corr_type))
            case ("SIMPLE")
                this%pressure_correction_type = SIMPLE
            case ("PISO")
                this%pressure_correction_type = PISO
            case ("PIMPLE")
                this%pressure_correction_type = PIMPLE
            case default
                print *, "Warning: Unknown pressure correction algorithm: ", trim(p_corr_type)
                print *, "Defaulting to SIMPLE"
                this%pressure_correction_type = SIMPLE
        end select
        
        ! Initialize inner iterations based on chosen algorithm
        select case (this%pressure_correction_type)
            case (SIMPLE)
                this%max_inner_iterations = 1
            case (PISO)
                this%max_inner_iterations = 2
            case (PIMPLE)
                this%max_inner_iterations = 3
            case default
                this%max_inner_iterations = 1
        end select
        
        ! Set current time to zero
        this%time = 0.0_WP
        
        ! Check if GPU is available
        this%use_gpu = (omp_get_num_devices() > 0)
    end subroutine incompressible_flow_solver_init
    
    !> Setup fields for the incompressible flow solver
    subroutine incompressible_flow_solver_setup_fields(this, u_field, v_field, w_field, p_field, viscosity_field)
        class(incompressible_flow_solver_t), intent(inout) :: this
        type(field_t), target, intent(inout) :: u_field, v_field, w_field, p_field, viscosity_field
        
        ! Set field pointers
        this%u => u_field
        this%v => v_field
        this%w => w_field
        this%p => p_field
        this%viscosity => viscosity_field
        
        ! Initialize internal fields
        call this%init_fields()
    end subroutine incompressible_flow_solver_setup_fields
    
    !> Set the boundary condition manager
    subroutine incompressible_flow_solver_setup_boundary_manager(this, bc_manager)
        class(incompressible_flow_solver_t), intent(inout) :: this
        type(boundary_manager_t), target, intent(in) :: bc_manager
        
        ! Store boundary condition manager
        this%bc_manager => bc_manager
        
        ! Now that we have boundary conditions, apply them to the fields
        call this%apply_boundary_conditions()
        
        ! Update mass flux now that boundaries are established
        call this%update_mass_flux()
        
        print *, "Boundary conditions applied to the flow fields"
    end subroutine incompressible_flow_solver_setup_boundary_manager
    
    !> Apply boundary conditions to all fields
    subroutine incompressible_flow_solver_apply_boundary_conditions(this)
        class(incompressible_flow_solver_t), intent(inout) :: this
        integer :: i, j, k, id, nx, ny, nz
        
        ! Apply boundary conditions via boundary manager if available
        if (associated(this%bc_manager)) then
            call this%bc_manager%apply(this%u)
            call this%bc_manager%apply(this%v)
            call this%bc_manager%apply(this%w)
            call this%bc_manager%apply(this%p)
        else
            ! Skip boundary condition application if bc_manager isn't set yet
            ! This might happen during initialization
            print *, "Warning: Boundary condition manager not set yet. Skipping boundary conditions."
            return
        end if
        
        ! Additional explicit enforcement for lid-driven cavity
        ! This ensures the top lid boundary condition is properly applied
        nx = this%mesh%nx
        ny = this%mesh%ny
        nz = this%mesh%nz
        
        ! Apply north wall (lid) velocity: u=1.0, v=0.0, w=0.0
        j = ny  ! Top row
        do k = 1, nz
            do i = 1, nx
                id = i + (j-1)*nx + (k-1)*nx*ny
                this%u%data(id) = 1.0_WP  ! Moving lid with u=1.0
                this%v%data(id) = 0.0_WP  ! No vertical velocity
                this%w%data(id) = 0.0_WP  ! No z-direction velocity
            end do
        end do
        
        ! Ensure all other walls have no-slip conditions (u=v=w=0)
        ! South wall (j=1)
        j = 1
        do k = 1, nz
            do i = 1, nx
                id = i + (j-1)*nx + (k-1)*nx*ny
                this%u%data(id) = 0.0_WP
                this%v%data(id) = 0.0_WP
                this%w%data(id) = 0.0_WP
            end do
        end do
        
        ! West wall (i=1)
        i = 1
        do k = 1, nz
            do j = 1, ny
                id = i + (j-1)*nx + (k-1)*nx*ny
                this%u%data(id) = 0.0_WP
                this%v%data(id) = 0.0_WP
                this%w%data(id) = 0.0_WP
            end do
        end do
        
        ! East wall (i=nx)
        i = nx
        do k = 1, nz
            do j = 1, ny
                id = i + (j-1)*nx + (k-1)*nx*ny
                this%u%data(id) = 0.0_WP
                this%v%data(id) = 0.0_WP
                this%w%data(id) = 0.0_WP
            end do
        end do
    end subroutine incompressible_flow_solver_apply_boundary_conditions
    
    !> Advance the solution by one time step
    subroutine incompressible_flow_solver_advance(this)
        class(incompressible_flow_solver_t), intent(inout) :: this
        integer :: outer_iter, inner_iter
        logical :: converged
        real(WP) :: max_u, max_v, max_p
        integer :: i, j, k, id, nx, ny, nz, ncells
        
        ! Initialize convergence flag
        converged = .false.
        
        ! Get dimensions
        ncells = this%mesh%topo%ncells
        nx = this%mesh%nx
        ny = this%mesh%ny
        nz = this%mesh%nz

        ! Debug: Add small perturbation to kickstart flow if first iteration 
        ! This is for the first iteration of the whole simulation
        if (this%time < this%dt*0.5) then
            ! Add slight perturbation to internal cells' velocity to kickstart flow
            do k = 1, nz
                do j = 1, ny
                    do i = 1, nx
                        id = i + (j-1)*nx + (k-1)*nx*ny
                        ! Skip boundaries
                        if (i > 1 .and. i < nx .and. j > 1 .and. j < ny) then
                            ! Add small velocity perturbation based on position
                            this%u%data(id) = this%u%data(id) + 0.01_WP * sin(real(j,WP)/real(ny,WP) * 3.14159_WP)
                            this%v%data(id) = this%v%data(id) + 0.01_WP * cos(real(i,WP)/real(nx,WP) * 3.14159_WP)
                        end if
                    end do
                end do
            end do
            
            ! Apply boundary conditions to ensure correct values at walls
            call this%apply_boundary_conditions()
        end if
        
        ! Make sure boundary conditions are applied correctly from the start
        call this%apply_boundary_conditions()
        
        ! SIMPLE algorithm outer iteration loop
        do outer_iter = 1, this%max_iterations
            ! 1. Momentum predictor step
            call this%momentum_predictor()
            
            ! Apply boundary conditions after momentum predictor to enforce lid velocity
            call this%apply_boundary_conditions()
            
            ! Inner correction loop (PISO/PIMPLE would do multiple corrections)
            do inner_iter = 1, this%max_inner_iterations
                ! 2. Pressure correction step
                call this%pressure_correction()
                
                ! 3. Momentum correction step
                call this%momentum_corrector()
                
                ! 4. Update mass flux at faces
                call this%update_mass_flux()
            end do
            
            ! 5. Apply boundary conditions
            call this%apply_boundary_conditions()
            
            ! 6. Check for convergence
            call this%check_convergence(converged)
            if (converged) then
                exit
            end if
        end do
        
        ! Update time
        this%time = this%time + this%dt
    end subroutine incompressible_flow_solver_advance
    
    !> Momentum predictor step (solving momentum equations for intermediate velocity)
    subroutine incompressible_flow_solver_momentum_predictor(this)
        class(incompressible_flow_solver_t), intent(inout) :: this
        integer :: i, j, k, id, ncells, nx, ny, nz
        real(WP) :: viscous_term_u, viscous_term_v, viscous_term_w
        real(WP) :: convective_term_u, convective_term_v, convective_term_w
        real(WP) :: pressure_grad_u, pressure_grad_v, pressure_grad_w
        real(WP) :: dx, dy, dz, visc
        real(WP) :: u_e, u_w, u_n, u_s, u_t, u_b   ! East, west, north, south, top, bottom values
        real(WP) :: v_e, v_w, v_n, v_s, v_t, v_b   ! for neighboring cells
        real(WP) :: w_e, w_w, w_n, w_s, w_t, w_b
        real(WP) :: p_e, p_w, p_n, p_s, p_t, p_b
        
        ! Number of cells and grid dimensions
        ncells = this%mesh%topo%ncells
        nx = this%mesh%nx
        ny = this%mesh%ny
        nz = this%mesh%nz
        
        ! Grid spacing (assuming uniform grid)
        dx = (this%mesh%xmax - this%mesh%xmin) / real(nx, WP)
        dy = (this%mesh%ymax - this%mesh%ymin) / real(ny, WP)
        dz = (this%mesh%zmax - this%mesh%zmin) / real(nz, WP)
        
        ! Copy current velocity fields to intermediate fields
        this%u_star%data = this%u%data
        this%v_star%data = this%v%data
        this%w_star%data = this%w%data
        
        ! If using GPU, update device data
        if (this%use_gpu) then
            !$omp target update to(this%u_star%data(1:ncells), this%v_star%data(1:ncells), this%w_star%data(1:ncells))
        end if
        
        ! Compute terms for momentum equations
        ! Note: This is a simplified approach using a structured grid assumption
        ! A real CFD code would handle unstructured grids more explicitly
        !$omp parallel do private(i, j, k, id, visc, &
        !$omp                     u_e, u_w, u_n, u_s, u_t, u_b, &
        !$omp                     v_e, v_w, v_n, v_s, v_t, v_b, &
        !$omp                     w_e, w_w, w_n, w_s, w_t, w_b, &
        !$omp                     p_e, p_w, p_n, p_s, p_t, p_b, &
        !$omp                     viscous_term_u, viscous_term_v, viscous_term_w, &
        !$omp                     convective_term_u, convective_term_v, convective_term_w, &
        !$omp                     pressure_grad_u, pressure_grad_v, pressure_grad_w)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    ! Convert 3D index to 1D cell ID
                    id = i + (j-1)*nx + (k-1)*nx*ny
                    
                    ! Get viscosity for this cell
                    visc = this%viscosity%data(id)
                    
                    ! Define indices for neighboring cells with boundary wrapping
                    ! East neighbor (i+1)
                    if (i < nx) then
                        u_e = this%u%data(id + 1)
                        v_e = this%v%data(id + 1)
                        w_e = this%w%data(id + 1)
                        p_e = this%p%data(id + 1)
                    else
                        ! At eastern boundary - use boundary condition instead
                        u_e = this%u%data(id)  ! Zero gradient for now
                        v_e = this%v%data(id)
                        w_e = this%w%data(id)
                        p_e = this%p%data(id)
                    end if
                    
                    ! West neighbor (i-1)
                    if (i > 1) then
                        u_w = this%u%data(id - 1)
                        v_w = this%v%data(id - 1)
                        w_w = this%w%data(id - 1)
                        p_w = this%p%data(id - 1)
                    else
                        ! At western boundary
                        u_w = this%u%data(id)
                        v_w = this%v%data(id)
                        w_w = this%w%data(id)
                        p_w = this%p%data(id)
                    end if
                    
                    ! North neighbor (j+1)
                    if (j < ny) then
                        u_n = this%u%data(id + nx)
                        v_n = this%v%data(id + nx)
                        w_n = this%w%data(id + nx)
                        p_n = this%p%data(id + nx)
                    else
                        ! At northern boundary - for lid-driven cavity, north wall is moving
                        u_n = 1.0_WP  ! Lid velocity in x-direction
                        v_n = 0.0_WP
                        w_n = 0.0_WP
                        p_n = this%p%data(id)
                    end if
                    
                    ! South neighbor (j-1)
                    if (j > 1) then
                        u_s = this%u%data(id - nx)
                        v_s = this%v%data(id - nx)
                        w_s = this%w%data(id - nx)
                        p_s = this%p%data(id - nx)
                    else
                        ! At southern boundary
                        u_s = 0.0_WP  ! No-slip wall
                        v_s = 0.0_WP
                        w_s = 0.0_WP
                        p_s = this%p%data(id)
                    end if
                    
                    ! Top neighbor (k+1) - for 2D problems, we'll ignore this
                    if (k < nz) then
                        u_t = this%u%data(id + nx*ny)
                        v_t = this%v%data(id + nx*ny)
                        w_t = this%w%data(id + nx*ny)
                        p_t = this%p%data(id + nx*ny)
                    else
                        u_t = this%u%data(id)
                        v_t = this%v%data(id)
                        w_t = this%w%data(id)
                        p_t = this%p%data(id)
                    end if
                    
                    ! Bottom neighbor (k-1)
                    if (k > 1) then
                        u_b = this%u%data(id - nx*ny)
                        v_b = this%v%data(id - nx*ny)
                        w_b = this%w%data(id - nx*ny)
                        p_b = this%p%data(id - nx*ny)
                    else
                        u_b = this%u%data(id)
                        v_b = this%v%data(id)
                        w_b = this%w%data(id)
                        p_b = this%p%data(id)
                    end if
                    
                    ! Compute viscous terms (Laplacian of velocity)
                    ! Using a simple second-order central difference scheme
                    viscous_term_u = visc * ( &
                        (u_e - 2.0_WP*this%u%data(id) + u_w) / (dx*dx) + &
                        (u_n - 2.0_WP*this%u%data(id) + u_s) / (dy*dy) + &
                        (u_t - 2.0_WP*this%u%data(id) + u_b) / (dz*dz) )
                    
                    viscous_term_v = visc * ( &
                        (v_e - 2.0_WP*this%v%data(id) + v_w) / (dx*dx) + &
                        (v_n - 2.0_WP*this%v%data(id) + v_s) / (dy*dy) + &
                        (v_t - 2.0_WP*this%v%data(id) + v_b) / (dz*dz) )
                    
                    viscous_term_w = visc * ( &
                        (w_e - 2.0_WP*this%w%data(id) + w_w) / (dx*dx) + &
                        (w_n - 2.0_WP*this%w%data(id) + w_s) / (dy*dy) + &
                        (w_t - 2.0_WP*this%w%data(id) + w_b) / (dz*dz) )
                    
                    ! Simplified convective terms (first-order upwind)
                    ! For lid-driven cavity, this term drives the flow
                    convective_term_u = this%u%data(id) * (this%u%data(id) - u_w) / dx + &
                                       this%v%data(id) * (this%u%data(id) - u_s) / dy
                    
                    convective_term_v = this%u%data(id) * (this%v%data(id) - v_w) / dx + &
                                       this%v%data(id) * (this%v%data(id) - v_s) / dy
                    
                    convective_term_w = 0.0_WP  ! 2D flow for lid-driven cavity
                    
                    ! Pressure gradient terms
                    pressure_grad_u = (p_e - p_w) / (2.0_WP * dx)
                    pressure_grad_v = (p_n - p_s) / (2.0_WP * dy)
                    pressure_grad_w = 0.0_WP  ! 2D flow
                    
                    ! Update u_star with under-relaxation
                    this%u_star%data(id) = this%u%data(id) + &
                                        this%dt * (viscous_term_u - convective_term_u - pressure_grad_u)
                    
                    ! Apply under-relaxation
                    this%u_star%data(id) = (1.0_WP - this%alpha_u) * this%u%data(id) + &
                                         this%alpha_u * this%u_star%data(id)
                    
                    ! Update v_star with under-relaxation
                    this%v_star%data(id) = this%v%data(id) + &
                                        this%dt * (viscous_term_v - convective_term_v - pressure_grad_v)
                    
                    ! Apply under-relaxation
                    this%v_star%data(id) = (1.0_WP - this%alpha_u) * this%v%data(id) + &
                                         this%alpha_u * this%v_star%data(id)
                    
                    ! Update w_star with under-relaxation
                    this%w_star%data(id) = this%w%data(id) + &
                                        this%dt * (viscous_term_w - convective_term_w - pressure_grad_w)
                    
                    ! Apply under-relaxation
                    this%w_star%data(id) = (1.0_WP - this%alpha_u) * this%w%data(id) + &
                                         this%alpha_u * this%w_star%data(id)
                end do
            end do
        end do
        !$omp end parallel do
        
        ! If using GPU, copy results back to host
        if (this%use_gpu) then
            !$omp target update from(this%u_star%data(1:ncells), this%v_star%data(1:ncells), this%w_star%data(1:ncells))
        end if
    end subroutine incompressible_flow_solver_momentum_predictor
    
    !> Pressure correction step (solving pressure correction equation)
    subroutine incompressible_flow_solver_pressure_correction(this)
        class(incompressible_flow_solver_t), intent(inout) :: this
        integer :: i, j, k, id, iter, ncells, nx, ny, nz
        real(WP) :: dx, dy, dz, source, residual, max_residual
        real(WP) :: p_e, p_w, p_n, p_s, p_t, p_b
        real(WP) :: u_e, u_w, v_n, v_s  ! Face velocities for divergence calculation
        real(WP) :: div_u, coef
        real(WP), allocatable :: p_prime_old(:)
        integer, parameter :: MAX_PRESSURE_ITER = 50
        real(WP), parameter :: PRESSURE_TOLERANCE = 1.0e-4_WP
        real(WP) :: max_p  ! For pressure perturbation check
        
        ! Number of cells and grid dimensions
        ncells = this%mesh%topo%ncells
        nx = this%mesh%nx
        ny = this%mesh%ny
        nz = this%mesh%nz
        
        ! Grid spacing (assuming uniform grid)
        dx = (this%mesh%xmax - this%mesh%xmin) / real(nx, WP)
        dy = (this%mesh%ymax - this%mesh%ymin) / real(ny, WP)
        dz = (this%mesh%zmax - this%mesh%zmin) / real(nz, WP)
        
        ! Check if pressure is all zeros, and perturb if needed to kickstart
        max_p = 0.0_WP
        do i = 1, ncells
            max_p = max(max_p, abs(this%p%data(i)))
        end do
        
        if (max_p < 1.0e-10_WP .and. this%time < this%dt*2.0_WP) then
            do k = 1, nz
                do j = 2, ny-1
                    do i = 2, nx-1
                        id = i + (j-1)*nx + (k-1)*nx*ny
                        ! Small pressure perturbation in the center region
                        if (i > nx/4 .and. i < 3*nx/4 .and. j > ny/4 .and. j < 3*ny/4) then
                            this%p%data(id) = 0.001_WP * sin(real(i,WP)/real(nx,WP) * 3.14159_WP) * &
                                             cos(real(j,WP)/real(ny,WP) * 3.14159_WP)
                        end if
                    end do
                end do
            end do
        end if
        
        ! Initialize pressure correction to zero
        this%p_prime%data = 0.0_WP
        
        ! Allocate temporary array for Jacobi iteration
        allocate(p_prime_old(ncells))
        
        ! Jacobi iteration for pressure correction equation
        do iter = 1, MAX_PRESSURE_ITER
            ! Save current solution
            p_prime_old = this%p_prime%data
            
            ! Reset maximum residual
            max_residual = 0.0_WP
            
            ! Solve pressure correction equation
            !$omp parallel do private(i, j, k, id, p_e, p_w, p_n, p_s, p_t, p_b, &
            !$omp                     u_e, u_w, v_n, v_s, div_u, source, coef, residual) &
            !$omp reduction(max:max_residual)
            do k = 1, nz
                do j = 1, ny
                    do i = 1, nx
                        ! Convert 3D index to 1D cell ID
                        id = i + (j-1)*nx + (k-1)*nx*ny
                        
                        ! Skip boundary cells for simplicity
                        if (i > 1 .and. i < nx .and. j > 1 .and. j < ny) then
                            ! Get neighbor cell pressure corrections
                            p_e = p_prime_old(id + 1)    ! East
                            p_w = p_prime_old(id - 1)    ! West
                            p_n = p_prime_old(id + nx)   ! North
                            p_s = p_prime_old(id - nx)   ! South
                            
                            if (k < nz) then
                                p_t = p_prime_old(id + nx*ny) ! Top
                            else
                                p_t = p_prime_old(id)        ! Use cell value at boundary
                            end if
                            
                            if (k > 1) then
                                p_b = p_prime_old(id - nx*ny) ! Bottom
                            else
                                p_b = p_prime_old(id)        ! Use cell value at boundary
                            end if
                            
                            ! Better computation of face velocities for divergence calculation
                            u_e = 0.5_WP * (this%u_star%data(id) + this%u_star%data(id+1))    ! u at east face
                            u_w = 0.5_WP * (this%u_star%data(id) + this%u_star%data(id-1))    ! u at west face
                            v_n = 0.5_WP * (this%v_star%data(id) + this%v_star%data(id+nx))   ! v at north face
                            v_s = 0.5_WP * (this%v_star%data(id) + this%v_star%data(id-nx))   ! v at south face
                            
                            ! Compute velocity divergence using face velocities
                            div_u = (u_e - u_w) / dx + (v_n - v_s) / dy
                            
                            ! Source term scaled by density/dt
                            source = div_u / this%dt
                            
                            ! Coefficient for central cell
                            coef = 2.0_WP / (dx*dx) + 2.0_WP / (dy*dy) + 2.0_WP / (dz*dz)
                            
                            ! Update pressure correction using Jacobi iteration
                            this%p_prime%data(id) = (source + &
                                p_e / (dx*dx) + p_w / (dx*dx) + &
                                p_n / (dy*dy) + p_s / (dy*dy) + &
                                p_t / (dz*dz) + p_b / (dz*dz)) / coef
                            
                            ! Compute local residual
                            residual = abs(this%p_prime%data(id) - p_prime_old(id))
                            max_residual = max(max_residual, residual)
                        end if
                    end do
                end do
            end do
            !$omp end parallel do
            
            ! Add direct enforcement of boundary conditions for p_prime
            ! Set zero pressure gradient at the boundaries
            do k = 1, nz
                do i = 1, nx
                    ! South wall (j=1)
                    this%p_prime%data(i + 0*nx + (k-1)*nx*ny) = this%p_prime%data(i + 1*nx + (k-1)*nx*ny)
                    
                    ! North wall (j=ny)
                    this%p_prime%data(i + (ny-1)*nx + (k-1)*nx*ny) = this%p_prime%data(i + (ny-2)*nx + (k-1)*nx*ny)
                end do
                
                do j = 1, ny
                    ! West wall (i=1)
                    this%p_prime%data(1 + (j-1)*nx + (k-1)*nx*ny) = this%p_prime%data(2 + (j-1)*nx + (k-1)*nx*ny)
                    
                    ! East wall (i=nx)
                    this%p_prime%data(nx + (j-1)*nx + (k-1)*nx*ny) = this%p_prime%data(nx-1 + (j-1)*nx + (k-1)*nx*ny)
                end do
            end do
            
            ! Check for convergence
            if (max_residual < PRESSURE_TOLERANCE) then
                exit
            end if
        end do
        
        ! Apply pressure correction to pressure field with under-relaxation
        do i = 1, ncells
            this%p%data(i) = this%p%data(i) + this%alpha_p * this%p_prime%data(i)
        end do
        
        ! Free temporary array
        deallocate(p_prime_old)
        
        ! If using GPU, copy results back to host
        if (this%use_gpu) then
            !$omp target update from(this%p_prime%data(1:ncells), this%p%data(1:ncells))
        end if
    end subroutine incompressible_flow_solver_pressure_correction
    
    !> Momentum corrector step (correcting velocity based on pressure correction)
    subroutine incompressible_flow_solver_momentum_corrector(this)
        class(incompressible_flow_solver_t), intent(inout) :: this
        integer :: i, j, k, id, ncells, nx, ny, nz
        real(WP) :: pressure_corr_grad_u, pressure_corr_grad_v, pressure_corr_grad_w
        real(WP) :: dx, dy, dz
        real(WP) :: p_e, p_w, p_n, p_s, p_t, p_b
        
        ! Number of cells and grid dimensions
        ncells = this%mesh%topo%ncells
        nx = this%mesh%nx
        ny = this%mesh%ny
        nz = this%mesh%nz
        
        ! Grid spacing (assuming uniform grid)
        dx = (this%mesh%xmax - this%mesh%xmin) / real(nx, WP)
        dy = (this%mesh%ymax - this%mesh%ymin) / real(ny, WP)
        dz = (this%mesh%zmax - this%mesh%zmin) / real(nz, WP)
        
        ! OpenMP parallel loop for correcting velocities
        !$omp parallel do private(i, j, k, id, p_e, p_w, p_n, p_s, p_t, p_b, &
        !$omp                     pressure_corr_grad_u, pressure_corr_grad_v, pressure_corr_grad_w)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    ! Convert 3D index to 1D cell ID
                    id = i + (j-1)*nx + (k-1)*nx*ny
                    
                    ! Skip boundary cells
                    if (i > 1 .and. i < nx .and. j > 1 .and. j < ny) then
                        ! Get pressure correction at neighboring cells for gradient calculation
                        p_e = this%p_prime%data(id + 1)    ! East
                        p_w = this%p_prime%data(id - 1)    ! West
                        p_n = this%p_prime%data(id + nx)   ! North
                        p_s = this%p_prime%data(id - nx)   ! South
                        
                        if (k < nz) then
                            p_t = this%p_prime%data(id + nx*ny) ! Top
                        else
                            p_t = this%p_prime%data(id)        ! Use cell value at boundary
                        end if
                        
                        if (k > 1) then
                            p_b = this%p_prime%data(id - nx*ny) ! Bottom
                        else
                            p_b = this%p_prime%data(id)        ! Use cell value at boundary
                        end if
                        
                        ! Compute pressure correction gradients using central differences
                        pressure_corr_grad_u = (p_e - p_w) / (2.0_WP * dx)
                        pressure_corr_grad_v = (p_n - p_s) / (2.0_WP * dy)
                        pressure_corr_grad_w = (p_t - p_b) / (2.0_WP * dz)
                        
                        ! Correct velocities using pressure correction gradient
                        ! Note: Pressure gradient is divided by density (assumed 1.0 here)
                        this%u%data(id) = this%u_star%data(id) - this%dt * pressure_corr_grad_u
                        this%v%data(id) = this%v_star%data(id) - this%dt * pressure_corr_grad_v
                        this%w%data(id) = this%w_star%data(id) - this%dt * pressure_corr_grad_w
                    else
                        ! For boundary cells, just use the predicted velocities
                        this%u%data(id) = this%u_star%data(id)
                        this%v%data(id) = this%v_star%data(id)
                        this%w%data(id) = this%w_star%data(id)
                    end if
                end do
            end do
        end do
        !$omp end parallel do
        
        ! If using GPU, update device data
        if (this%use_gpu) then
            !$omp target update to(this%u%data(1:ncells), this%v%data(1:ncells), this%w%data(1:ncells))
        end if
    end subroutine incompressible_flow_solver_momentum_corrector
    
    !> Update mass flux at faces (for continuity equation)
    subroutine incompressible_flow_solver_update_mass_flux(this)
        class(incompressible_flow_solver_t), intent(inout) :: this
        integer :: i, j, k, face_id, id, id_m1
        integer :: ncells, nfaces, nx, ny, nz
        real(WP) :: dx, dy, dz, face_area
        
        ! Ensure fields are properly set up
        if (.not. associated(this%u) .or. .not. associated(this%v) .or. &
            .not. associated(this%w) .or. .not. associated(this%mesh)) then
            print *, "Warning: Fields not properly set up in update_mass_flux. Skipping calculation."
            return
        end if
        
        ! Get mesh dimensions
        ncells = this%mesh%topo%ncells
        nfaces = this%mesh%topo%nfaces
        nx = this%mesh%nx
        ny = this%mesh%ny
        nz = this%mesh%nz
        
        ! Grid spacing
        dx = (this%mesh%xmax - this%mesh%xmin) / real(nx, WP)
        dy = (this%mesh%ymax - this%mesh%ymin) / real(ny, WP)
        dz = (this%mesh%zmax - this%mesh%zmin) / real(nz, WP)
        
        ! Reset face velocities and mass flux
        this%face_u%data = 0.0_WP
        this%face_v%data = 0.0_WP
        this%face_w%data = 0.0_WP
        this%mass_flux%data = 0.0_WP
        
        ! Calculate face velocities and mass flux
        ! For structured grid, we'll process internal x-faces, y-faces, and z-faces separately
        
        ! Process x-faces (constant i planes)
        face_id = 1
        do k = 1, nz
            do j = 1, ny
                do i = 0, nx  ! One more face than cells in each direction
                    if (i > 0 .and. i < nx) then
                        ! Internal face - interpolate from adjacent cells
                        id = i + (j-1)*nx + (k-1)*nx*ny  ! Right cell ID
                        
                        ! Linear interpolation of u velocity to face
                        this%face_u%data(face_id) = 0.5_WP * (this%u%data(id) + this%u%data(id-1))
                        
                        ! Face area for flux calculation (dy*dz)
                        face_area = dy * dz
                        
                        ! Mass flux = density * velocity * area (density=1 for incompressible)
                        this%mass_flux%data(face_id) = this%face_u%data(face_id) * face_area
                    else if (i == 0) then
                        ! West boundary face
                        ! No need to compute id, just directly set boundary values
                        this%face_u%data(face_id) = 0.0_WP  ! No-slip for lid-driven cavity
                        face_area = dy * dz
                        this%mass_flux%data(face_id) = this%face_u%data(face_id) * face_area
                    else if (i == nx) then
                        ! East boundary face
                        id = nx + (j-1)*nx + (k-1)*nx*ny  ! Last cell ID
                        this%face_u%data(face_id) = 0.0_WP  ! No-slip for lid-driven cavity
                        face_area = dy * dz
                        this%mass_flux%data(face_id) = this%face_u%data(face_id) * face_area
                    end if
                    
                    face_id = face_id + 1
                end do
            end do
        end do
        
        ! Process y-faces (constant j planes)
        do k = 1, nz
            do j = 0, ny  ! One more face than cells in each direction
                do i = 1, nx
                    if (j > 0 .and. j < ny) then
                        ! Internal face - interpolate from adjacent cells
                        id = i + (j-1)*nx + (k-1)*nx*ny  ! Upper cell ID
                        
                        ! Linear interpolation of v velocity to face
                        this%face_v%data(face_id) = 0.5_WP * (this%v%data(id) + this%v%data(id-nx))
                        
                        ! Face area for flux calculation (dx*dz)
                        face_area = dx * dz
                        
                        ! Mass flux = density * velocity * area (density=1 for incompressible)
                        this%mass_flux%data(face_id) = this%face_v%data(face_id) * face_area
                    else if (j == 0) then
                        ! South boundary face
                        ! No need to compute id, just directly set boundary values
                        this%face_v%data(face_id) = 0.0_WP  ! No-slip for lid-driven cavity
                        face_area = dx * dz
                        this%mass_flux%data(face_id) = this%face_v%data(face_id) * face_area
                    else if (j == ny) then
                        ! North boundary face - moving lid with u=1.0
                        id = i + (ny-1)*nx + (k-1)*nx*ny  ! Last cell ID in this column
                        ! For the lid-driven cavity, the north boundary has u=1.0
                        if (i >= 1 .and. i <= nx) then
                            this%face_v%data(face_id) = 0.0_WP  ! No vertical velocity at the lid
                        end if
                        face_area = dx * dz
                        this%mass_flux%data(face_id) = this%face_v%data(face_id) * face_area
                    end if
                    
                    face_id = face_id + 1
                end do
            end do
        end do
        
        ! Process z-faces (constant k planes) - simplified for 2D
        do k = 0, nz  ! One more face than cells in each direction
            do j = 1, ny
                do i = 1, nx
                    if (k > 0 .and. k < nz) then
                        ! Internal face - interpolate from adjacent cells
                        id = i + (j-1)*nx + (k-1)*nx*ny  ! Top cell ID
                        
                        ! Linear interpolation of w velocity to face
                        this%face_w%data(face_id) = 0.5_WP * (this%w%data(id) + this%w%data(id-nx*ny))
                        
                        ! Face area for flux calculation (dx*dy)
                        face_area = dx * dy
                        
                        ! Mass flux = density * velocity * area (density=1 for incompressible)
                        this%mass_flux%data(face_id) = this%face_w%data(face_id) * face_area
                    else if (k == 0) then
                        ! Bottom boundary face - symmetry plane for 2D simulation
                        ! No need to compute id, just directly set boundary values
                        this%face_w%data(face_id) = 0.0_WP
                        face_area = dx * dy
                        this%mass_flux%data(face_id) = 0.0_WP
                    else if (k == nz) then
                        ! Top boundary face - symmetry plane for 2D simulation
                        this%face_w%data(face_id) = 0.0_WP
                        face_area = dx * dy
                        this%mass_flux%data(face_id) = 0.0_WP
                    end if
                    
                    face_id = face_id + 1
                end do
            end do
        end do
        
        ! If using GPU, update device data
        if (this%use_gpu) then
            !$omp target update to(this%face_u%data(1:nfaces), this%face_v%data(1:nfaces), &
            !$omp                   this%face_w%data(1:nfaces), this%mass_flux%data(1:nfaces))
        end if
    end subroutine incompressible_flow_solver_update_mass_flux
    
    !> Check convergence
    subroutine incompressible_flow_solver_check_convergence(this, converged)
        class(incompressible_flow_solver_t), intent(inout) :: this
        logical, intent(out) :: converged
        integer :: i, ncells
        real(WP) :: residual_sum_u, residual_sum_v, residual_sum_w, residual_sum_p
        real(WP) :: diff_u, diff_v, diff_w, diff_p
        
        ! Number of cells
        ncells = this%mesh%topo%ncells
        
        ! Initialize residual sums
        residual_sum_u = 0.0_WP
        residual_sum_v = 0.0_WP
        residual_sum_w = 0.0_WP
        residual_sum_p = 0.0_WP
        
        ! Compute residuals based on difference between current and intermediate values
        do i = 1, ncells
            ! Calculate differences
            diff_u = abs(this%u%data(i) - this%u_star%data(i))
            diff_v = abs(this%v%data(i) - this%v_star%data(i))
            diff_w = abs(this%w%data(i) - this%w_star%data(i))
            diff_p = abs(this%p_prime%data(i))
            
            ! Accumulate squared differences
            residual_sum_u = residual_sum_u + diff_u * diff_u
            residual_sum_v = residual_sum_v + diff_v * diff_v
            residual_sum_w = residual_sum_w + diff_w * diff_w
            residual_sum_p = residual_sum_p + diff_p * diff_p
        end do
        
        ! Compute root-mean-square residuals
        this%u_residual = sqrt(residual_sum_u / real(ncells, WP))
        this%v_residual = sqrt(residual_sum_v / real(ncells, WP))
        this%w_residual = sqrt(residual_sum_w / real(ncells, WP))
        this%p_residual = sqrt(residual_sum_p / real(ncells, WP))
        
        ! Check if all residuals are below tolerance
        converged = (this%u_residual < this%tolerance_velocity) .and. &
                   (this%v_residual < this%tolerance_velocity) .and. &
                   (this%w_residual < this%tolerance_velocity) .and. &
                   (this%p_residual < this%tolerance_pressure)
    end subroutine incompressible_flow_solver_check_convergence
    
    !> Cleanup
    subroutine incompressible_flow_solver_cleanup(this)
        class(incompressible_flow_solver_t), intent(inout) :: this
        
        ! Cleanup internal fields
        call this%u_star%cleanup()
        call this%v_star%cleanup()
        call this%w_star%cleanup()
        call this%p_prime%cleanup()
        call this%mass_flux%cleanup()
        call this%face_u%cleanup()
        call this%face_v%cleanup()
        call this%face_w%cleanup()
        
        ! Reset pointers
        this%u => null()
        this%v => null()
        this%w => null()
        this%p => null()
        this%viscosity => null()
        this%bc_manager => null()
        this%mesh => null()
    end subroutine incompressible_flow_solver_cleanup

    !> Get residuals for velocity components (u, v, w) and pressure (p)
    subroutine incompressible_flow_solver_get_residuals(this, residuals)
        class(incompressible_flow_solver_t), intent(in) :: this
        real(WP), intent(out) :: residuals(:)
        
        ! Check if array is large enough
        if (size(residuals) >= 4) then
            ! Return residuals in correct order
            residuals(1) = this%u_residual
            residuals(2) = this%v_residual
            residuals(3) = this%w_residual
            residuals(4) = this%p_residual
        end if
    end subroutine incompressible_flow_solver_get_residuals

    !> Initialize the fields with small perturbations to help start the flow
    subroutine incompressible_flow_solver_init_fields(this)
        class(incompressible_flow_solver_t), intent(inout) :: this
        integer :: i, j, k, id, nx, ny, nz
        real(WP) :: x, y, z, dx, dy, dz
        real(WP) :: perturbation_magnitude = 0.05_WP  ! 5% of lid velocity
        
        print *, "Initializing fields for incompressible flow solver..."
        
        ! Initialize internal fields
        call this%u_star%init(this%mesh, "u_star", FIELD_TYPE_CELL, 1)
        call this%v_star%init(this%mesh, "v_star", FIELD_TYPE_CELL, 1)
        call this%w_star%init(this%mesh, "w_star", FIELD_TYPE_CELL, 1)
        call this%p_prime%init(this%mesh, "p_prime", FIELD_TYPE_CELL, 1)
        call this%mass_flux%init(this%mesh, "mass_flux", FIELD_TYPE_FACE, 1)
        call this%face_u%init(this%mesh, "face_u", FIELD_TYPE_FACE, 1)
        call this%face_v%init(this%mesh, "face_v", FIELD_TYPE_FACE, 1)
        call this%face_w%init(this%mesh, "face_w", FIELD_TYPE_FACE, 1)
        
        ! Copy to device if using GPU
        if (this%use_gpu) then
            call this%u_star%copy_to_device()
            call this%v_star%copy_to_device()
            call this%w_star%copy_to_device()
            call this%p_prime%copy_to_device()
            call this%mass_flux%copy_to_device()
            call this%face_u%copy_to_device()
            call this%face_v%copy_to_device()
            call this%face_w%copy_to_device()
        end if
        
        print *, "Adding small perturbations to velocity field..."
        
        ! Get grid dimensions
        nx = this%mesh%nx
        ny = this%mesh%ny
        nz = this%mesh%nz
        
        ! Get cell sizes
        dx = (this%mesh%xmax - this%mesh%xmin) / real(nx, WP)
        dy = (this%mesh%ymax - this%mesh%ymin) / real(ny, WP)
        dz = (this%mesh%zmax - this%mesh%zmin) / real(nz, WP)
        
        ! Add small perturbations to internal cells
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    id = i + (j-1)*nx + (k-1)*nx*ny
                    
                    ! Calculate cell center coordinates
                    x = this%mesh%xmin + (i - 0.5_WP) * dx
                    y = this%mesh%ymin + (j - 0.5_WP) * dy
                    
                    ! Add sinusoidal perturbation to u velocity (stronger near the lid)
                    this%u%data(id) = perturbation_magnitude * sin(2.0_WP * 3.14159_WP * x) * (y / this%mesh%ymax)
                    
                    ! Add sinusoidal perturbation to v velocity
                    this%v%data(id) = 0.5_WP * perturbation_magnitude * sin(2.0_WP * 3.14159_WP * y) * (1.0_WP - y / this%mesh%ymax)
                end do
            end do
        end do
        
        ! Note: boundary conditions will be applied after bc_manager is set up
        print *, "Initial field setup complete."
    end subroutine incompressible_flow_solver_init_fields

end module incompressible_flow 