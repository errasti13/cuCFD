module diffusion
    use omp_lib
    use mesh, only: mesh_t, WP
    use field, only: field_t, FIELD_TYPE_CELL, FIELD_TYPE_FACE
    implicit none
    
    !> Diffusion solver data structure
    type :: diffusion_solver_t
        ! Fields
        type(field_t), pointer :: temperature => null()  ! Temperature field
        type(field_t), pointer :: diffusivity => null()  ! Diffusivity field
        type(field_t), pointer :: source => null()       ! Source term
        
        ! Internal fields for computation
        type(field_t) :: flux                ! Face fluxes
        type(field_t) :: flux_div            ! Flux divergence
        
        ! Solver configuration
        real(WP) :: dt                       ! Time step
        integer :: max_iterations            ! Maximum number of iterations
        real(WP) :: convergence_tolerance    ! Convergence tolerance
        
        ! Residual tracking
        real(WP) :: current_residual = 0.0_WP  ! Current residual for temperature equation
        
        ! Mesh pointer
        type(mesh_t), pointer :: mesh => null()
        
        ! GPU acceleration flag
        logical :: use_gpu = .false.  ! Will be set based on device availability
        
    contains
        procedure :: init => diffusion_solver_init
        procedure :: setup_fields => diffusion_solver_setup_fields
        procedure :: compute_fluxes => diffusion_solver_compute_fluxes
        procedure :: compute_divergence => diffusion_solver_compute_divergence
        procedure :: advance => diffusion_solver_advance
        procedure :: get_residual => diffusion_solver_get_residual
        procedure :: cleanup => diffusion_solver_cleanup
    end type diffusion_solver_t
    
contains

    !> Initialize the diffusion solver
    subroutine diffusion_solver_init(this, mesh, dt, max_iter, tolerance)
        class(diffusion_solver_t), intent(inout) :: this
        type(mesh_t), target, intent(in) :: mesh
        real(WP), intent(in) :: dt
        integer, intent(in) :: max_iter
        real(WP), intent(in) :: tolerance
        
        ! Store mesh and solver parameters
        this%mesh => mesh
        this%dt = dt
        this%max_iterations = max_iter
        this%convergence_tolerance = tolerance
        
        ! Check if GPU is available
        this%use_gpu = (omp_get_num_devices() > 0)
        
        ! Print GPU usage status
        if (this%use_gpu) then
            print *, "Diffusion solver: Using GPU acceleration"
        else
            print *, "Diffusion solver: Using CPU only"
        end if
    end subroutine diffusion_solver_init
    
    !> Setup fields for diffusion solver
    subroutine diffusion_solver_setup_fields(this, temperature, diffusivity, source)
        class(diffusion_solver_t), intent(inout) :: this
        type(field_t), target, intent(in) :: temperature
        type(field_t), target, intent(in) :: diffusivity
        type(field_t), target, intent(in) :: source
        
        ! Store field pointers
        this%temperature => temperature
        this%diffusivity => diffusivity
        this%source => source
        
        ! Initialize internal fields
        call this%flux%init(this%mesh, "flux", FIELD_TYPE_CELL, 1)
        call this%flux_div%init(this%mesh, "flux_divergence", FIELD_TYPE_CELL, 1)
        
        ! Copy to device if using GPU
        if (this%use_gpu) then
            call this%flux%copy_to_device()
            call this%flux_div%copy_to_device()
        end if
    end subroutine diffusion_solver_setup_fields
    
    !> Compute diffusive fluxes at cell faces
    subroutine diffusion_solver_compute_fluxes(this)
        class(diffusion_solver_t), intent(inout) :: this
        integer :: i, ncells
        
        ! Skip if not all fields are available
        if (.not. associated(this%temperature) .or. &
            .not. associated(this%diffusivity)) then
            print *, "Error: Fields not set for diffusion solver"
            return
        end if
        
        ! Number of cells
        ncells = this%mesh%topo%ncells

        ! Use CPU version for now, but with optimized code paths
        ! This avoids complex OpenMP offload issues that are causing the data transfer errors
        
        ! Simple diffusion operation
        do i = 1, ncells
            ! Simple placeholder - just compute some heat diffusion
            ! In a real solver, we would use proper face-cell connectivity
            ! But for this demo, we'll use simplifications
            this%flux%data(i) = -1.0_WP * this%diffusivity%data(i) * this%temperature%data(i)
        end do
        
        ! If using GPU, copy the updated data back to the device
        if (this%use_gpu .and. this%flux%on_device) then
            ! Explicit data update
            !$omp target update to(this%flux%data(1:ncells))
        end if
    end subroutine diffusion_solver_compute_fluxes
    
    !> Compute divergence of fluxes for each cell
    subroutine diffusion_solver_compute_divergence(this)
        class(diffusion_solver_t), intent(inout) :: this
        integer :: i, ncells
        
        ! Number of cells
        ncells = this%mesh%topo%ncells
        
        ! Simple divergence calculation
        do i = 1, ncells
            ! Simple placeholder - add some diffusion effect
            this%flux_div%data(i) = 0.1_WP * this%flux%data(i)
        end do
        
        ! If using GPU, copy the updated data back to the device
        if (this%use_gpu .and. this%flux_div%on_device) then
            ! Explicit data update
            !$omp target update to(this%flux_div%data(1:ncells))
        end if
    end subroutine diffusion_solver_compute_divergence
    
    !> Advance solution by one time step
    subroutine diffusion_solver_advance(this)
        class(diffusion_solver_t), intent(inout) :: this
        integer :: i, ncells
        real(WP) :: local_dt ! Local variable for dt
        real(WP) :: residual_sum
        
        ! Compute fluxes at cell faces
        call this%compute_fluxes()
        
        ! Compute divergence of fluxes
        call this%compute_divergence()
        
        ! Number of cells and local dt
        ncells = this%mesh%topo%ncells
        local_dt = this%dt ! Copy dt to local variable
        
        ! Calculate residual before updating the solution
        residual_sum = 0.0_WP
        do i = 1, ncells
            residual_sum = residual_sum + (this%flux_div%data(i) + this%source%data(i))**2
        end do
        this%current_residual = sqrt(residual_sum / real(ncells, WP))
        
        ! Update temperature field
        if (this%use_gpu .and. this%temperature%on_device .and. &
            this%flux_div%on_device .and. this%source%on_device) then
            
            ! First check data sizes to ensure they match
            if (allocated(this%temperature%data) .and. allocated(this%flux_div%data) .and. allocated(this%source%data)) then
                print *, "Temperature data size:", size(this%temperature%data)
                print *, "Flux divergence data size:", size(this%flux_div%data)
                print *, "Source data size:", size(this%source%data)
            end if

            ! Combine data map and kernel launch into a single directive
            ! Explicitly map all required data, use firstprivate for local_dt
            ! Simplify mapping, relying more on implicit rules based on usage
            !$omp target teams distribute parallel do map(tofrom: this%temperature%data(1:ncells)) &            
            !$omp                                    firstprivate(local_dt) private(i)
            do i = 1, ncells
                this%temperature%data(i) = this%temperature%data(i) + &
                                         local_dt * (this%flux_div%data(i) + this%source%data(i))
            end do
            !$omp end target teams distribute parallel do
            
        else
            ! CPU computation
            do i = 1, ncells
                this%temperature%data(i) = this%temperature%data(i) + &
                                         this%dt * (this%flux_div%data(i) + this%source%data(i))
            end do
        end if
        
        ! Update boundary conditions
        call this%temperature%update_boundaries()
    end subroutine diffusion_solver_advance
    
    !> Get residual of the temperature equation
    function diffusion_solver_get_residual(this) result(residual)
        class(diffusion_solver_t), intent(in) :: this
        real(WP) :: residual
        
        ! Return stored residual
        residual = this%current_residual
    end function diffusion_solver_get_residual
    
    !> Clean up diffusion solver
    subroutine diffusion_solver_cleanup(this)
        class(diffusion_solver_t), intent(inout) :: this
        
        ! Clean up internal fields
        call this%flux%cleanup()
        call this%flux_div%cleanup()
        
        ! Remove field pointers
        this%temperature => null()
        this%diffusivity => null()
        this%source => null()
        
        ! Remove mesh pointer
        this%mesh => null()
    end subroutine diffusion_solver_cleanup

end module diffusion 