program main
    use omp_lib
    use mesh, only: mesh_t, WP
    use field, only: field_t
    use solution_io, only: solution_io_t
    use config_parser, only: config_t
    use solver_manager, only: solver_manager_t
    implicit none
    
    ! Configuration
    type(config_t) :: config
    
    ! Solver manager
    type(solver_manager_t) :: solver
    
    ! I/O
    type(solution_io_t) :: solution_io
    
    ! Timing
    real(WP) :: start_time, end_time, sim_time
    
    ! Command line arguments
    integer :: num_args
    character(len=256) :: config_file
    
    print *, "==================================================="
    print *, "cuCFD - CUDA-accelerated Computational Fluid Dynamics Solver"
    print *, "==================================================="
    print *, ""
    
    ! Process command line arguments
    num_args = command_argument_count()
    if (num_args < 1) then
        print *, "Error: No configuration file specified."
        print *, "Usage: cucfd <config_file>"
        stop
    end if
    
    ! Get configuration file name
    call get_command_argument(1, config_file)
    
    ! Load configuration
    call config%load(config_file)
    
    ! Initialize solver manager with configuration
    call solver%init(config_file)
    
    ! Initialize solution I/O
    call solution_io%init(solver%output_dir, solver%case_name, solver%output_frequency)
    
    ! Print simulation information
    print *, "Starting simulation with:"
    print *, "  * Time step: ", solver%dt
    print *, "  * End time: ", solver%end_time
    print *, "  * Max iterations: ", solver%max_iterations
    print *, "  * Output directory: ", trim(solver%output_dir)
    print *, "  * GPU acceleration: ", solver%use_gpu
    print *, ""
    
    ! Write initial solution
    if (solver%do_write_solution) then
        call solver%write_solution(solution_io)
    end if
    
    ! Start timing
    start_time = omp_get_wtime()
    
    ! Main simulation loop
    do while (.not. solver%simulation_finished)
        ! Advance solution by one time step
        call solver%advance()
        
        ! Print progress every 10% of iterations
        if (mod(solver%current_iteration, max(1, solver%max_iterations/10)) == 0) then
            print '(A,F6.1,A,I0,A,I0,A)', "  Progress: ", &
                100.0 * solver%current_time / real(solver%end_time), "% (", &
                solver%current_iteration, "/", solver%max_iterations, ")"
        end if
        
        ! Write solution if output frequency reached
        if (solver%do_write_solution .and. mod(solver%current_iteration, solver%output_frequency) == 0) then
            call solver%write_solution(solution_io)
        end if
    end do
    
    ! Stop timing
    end_time = omp_get_wtime()
    sim_time = end_time - start_time
    
    ! Write final solution if not already written
    if (solver%do_write_solution .and. mod(solver%current_iteration, solver%output_frequency) /= 0) then
        call solver%write_solution(solution_io)
    end if
    
    ! Print results
    print *, ""
    print *, "Simulation completed"
    print *, "Time steps: ", solver%current_iteration
    print *, "Execution time: ", sim_time, " seconds"
    
    ! Cleanup
    call solver%cleanup()
    call config%cleanup()
    
    print *, "Done."
    
end program main 