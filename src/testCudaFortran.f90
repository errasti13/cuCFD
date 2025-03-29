program openmp_gpu_vs_cpu
    use omp_lib
    use gpu_utils
    implicit none

    integer, parameter :: NUM_SIZES = 4
    integer, parameter :: ITERATIONS = 10
    integer, parameter :: WARMUP_ITER = 3
    integer, parameter :: BLOCK_SIZE = 256
    integer, dimension(NUM_SIZES) :: sizes = (/10000000, 50000000, 100000000, 200000000/)
    real, allocatable :: a(:), b(:), c(:), c_cpu(:)
    integer :: i, j, n, iter
    double precision :: t1, t2, t3, t4, cpu_total, gpu_total, speedup
    logical :: is_correct

    ! Set number of CPU threads
    call omp_set_num_threads(16)

    ! Initialize GPU
    call init_gpu()

    ! Print system information
    print '(A,I0)', "Number of CPU threads: ", omp_get_num_threads()
    print '(A)', "GPU execution enabled: ", is_gpu_available()
    print '(A)', ""

    ! Print header
    print '(A)', "Size        CPU Time(s)    GPU Time(s)    Speedup    Match?"
    print '(A)', "--------------------------------------------------------"

    do j = 1, NUM_SIZES
        n = sizes(j)
        
        allocate(a(n), b(n), c(n), c_cpu(n))
        
        !$omp parallel do simd
        do i = 1, n
            a(i) = sin(real(i))
            b(i) = cos(real(i))
        end do
        !$omp end parallel do simd

        if (is_gpu_available()) then
            !$omp target enter data map(to:a,b) map(alloc:c)
            do iter = 1, WARMUP_ITER
                !$omp target teams distribute parallel do num_teams(n/BLOCK_SIZE) thread_limit(BLOCK_SIZE)
                do i = 1, n
                    c(i) = sqrt(a(i)**2 + b(i)**2) * sin(a(i) + b(i))
                end do
                !$omp end target teams distribute parallel do
            end do
        endif

        cpu_total = 0.0d0
        gpu_total = 0.0d0

        do iter = 1, ITERATIONS
            ! CPU computation
            t1 = omp_get_wtime()
            !$omp parallel do simd
            do i = 1, n
                c_cpu(i) = sqrt(a(i)**2 + b(i)**2) * sin(a(i) + b(i))
            end do
            !$omp end parallel do simd
            t2 = omp_get_wtime()
            cpu_total = cpu_total + (t2 - t1)

            if (is_gpu_available()) then
                t3 = omp_get_wtime()
                !$omp target teams distribute parallel do num_teams(n/BLOCK_SIZE) thread_limit(BLOCK_SIZE)
                do i = 1, n
                    c(i) = sqrt(a(i)**2 + b(i)**2) * sin(a(i) + b(i))
                end do
                !$omp end target teams distribute parallel do
                t4 = omp_get_wtime()
                gpu_total = gpu_total + (t4 - t3)
            else
                gpu_total = 0.0d0
            endif
        end do

        if (is_gpu_available()) then
            !$omp target exit data map(from:c) map(delete:a,b)
        endif

        ! Calculate average times and speedup
        cpu_total = cpu_total / ITERATIONS
        gpu_total = gpu_total / ITERATIONS
        speedup = cpu_total / gpu_total
        is_correct = all(abs(c - c_cpu) < 1.0e-5)

        ! Print results in tabular format
        print '(I9,2F14.6,F11.2,L8)', n, cpu_total, gpu_total, speedup, is_correct

        deallocate(a, b, c, c_cpu)
    end do

    ! Cleanup
    call cleanup_gpu()

end program openmp_gpu_vs_cpu