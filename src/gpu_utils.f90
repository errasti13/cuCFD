module gpu_utils
    use omp_lib
    implicit none
    private
    public :: init_gpu, cleanup_gpu, is_gpu_available

    logical, public :: gpu_available = .false.
    
contains
    subroutine init_gpu()
        integer :: num_devices
        num_devices = omp_get_num_devices()
        gpu_available = (num_devices > 0)
        
        if (gpu_available) then
            !$omp target
            !$omp end target
            if (omp_get_num_devices() == 0) then
                print '(A)', "WARNING: GPU detected but not compatible"
                gpu_available = .false.
            else
                print '(A,I0)', "GPU initialized. Available devices: ", num_devices
            end if
        end if
    end subroutine

    logical function is_gpu_available()
        is_gpu_available = gpu_available
    end function

    subroutine cleanup_gpu()
        if (gpu_available) then
            !$omp target exit data
        end if
    end subroutine
end module gpu_utils
