module field
    use omp_lib
    use mesh, only: mesh_t, WP
    implicit none
    
    !> Field types
    integer, parameter :: FIELD_TYPE_CELL = 1    ! Cell-centered field
    integer, parameter :: FIELD_TYPE_FACE = 2    ! Face-centered field
    integer, parameter :: FIELD_TYPE_VERTEX = 3  ! Vertex-centered field
    
    !> Field data structure
    type :: field_t
        character(len=64) :: name            ! Field name
        integer :: field_type                ! Field type (cell/face/vertex)
        integer :: ncomp                     ! Number of components
        integer :: nentries                  ! Number of entries (cells/faces/vertices)
        logical :: on_device = .false.       ! Flag indicating if field is on device
        
        ! Field data (flattened for better memory access patterns)
        real(WP), allocatable :: data(:)     ! (nentries * ncomp)
        
        ! Pointer to mesh
        type(mesh_t), pointer :: mesh => null()
        
    contains
        procedure :: init => field_init
        procedure :: set_uniform => field_set_uniform
        procedure :: copy_to_device => field_copy_to_device
        procedure :: copy_from_device => field_copy_from_device
        procedure :: update_boundaries => field_update_boundaries
        procedure :: cleanup => field_cleanup
    end type field_t
    
contains

    !> Initialize a field
    subroutine field_init(this, mesh, name, field_type, ncomp)
        class(field_t), intent(inout) :: this
        type(mesh_t), target, intent(in) :: mesh
        character(len=*), intent(in) :: name
        integer, intent(in) :: field_type
        integer, intent(in) :: ncomp
        
        this%name = name
        this%field_type = field_type
        this%ncomp = ncomp
        this%mesh => mesh
        
        ! Set number of entries based on field type
        select case (field_type)
        case (FIELD_TYPE_CELL)
            this%nentries = mesh%topo%ncells
        case (FIELD_TYPE_FACE)
            this%nentries = mesh%topo%nfaces
        case (FIELD_TYPE_VERTEX)
            this%nentries = mesh%topo%nvertices
        case default
            print *, "Error: Unknown field type ", field_type
            stop
        end select
        
        ! Allocate field data
        allocate(this%data(this%nentries * this%ncomp))
        
        ! Initialize to zero
        this%data = 0.0_WP
    end subroutine field_init
    
    !> Set uniform value for field
    subroutine field_set_uniform(this, values)
        class(field_t), intent(inout) :: this
        real(WP), intent(in) :: values(this%ncomp)
        integer :: i, j, idx
        
        ! Set uniform value for each component
        do i = 1, this%nentries
            do j = 1, this%ncomp
                idx = (i-1)*this%ncomp + j
                this%data(idx) = values(j)
            end do
        end do
        
        ! Update field on device if needed
        if (this%on_device) then
            call this%copy_to_device()
        end if
    end subroutine field_set_uniform
    
    !> Copy field to device
    subroutine field_copy_to_device(this)
        class(field_t), intent(inout) :: this
        
        ! Only copy if we have OpenMP target device
        if (omp_get_num_devices() > 0) then
            !$omp target enter data map(to:this%data)
            this%on_device = .true.
        else
            this%on_device = .false.
        end if
    end subroutine field_copy_to_device
    
    !> Copy field from device
    subroutine field_copy_from_device(this)
        class(field_t), intent(inout) :: this
        
        ! Only copy if field is on device
        if (this%on_device) then
            !$omp target exit data map(from:this%data)
        end if
    end subroutine field_copy_from_device
    
    !> Update boundary values for a field
    subroutine field_update_boundaries(this)
        class(field_t), intent(inout) :: this
        
        ! Only needed for cell-centered fields
        if (this%field_type /= FIELD_TYPE_CELL) return
        
        ! TODO: Implement boundary condition application
        ! Will need to loop over boundary patches and apply
        ! appropriate conditions for each patch type
        
    end subroutine field_update_boundaries
    
    !> Clean up field
    subroutine field_cleanup(this)
        class(field_t), intent(inout) :: this
        
        ! Copy from device if needed
        if (this%on_device) then
            call this%copy_from_device()
        end if
        
        ! Free memory
        if (allocated(this%data)) deallocate(this%data)
        
        ! Remove mesh reference
        this%mesh => null()
    end subroutine field_cleanup

end module field 