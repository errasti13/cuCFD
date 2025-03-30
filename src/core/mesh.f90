module mesh
    use omp_lib
    implicit none
    
    !> Mesh types and constants
    integer, parameter :: SP = kind(1.0)
    integer, parameter :: DP = kind(1.0d0)
    
    !> Use double precision by default
    integer, parameter :: WP = DP
    
    !> Face index definitions
    integer, parameter :: FACE_WEST = 1
    integer, parameter :: FACE_EAST = 2
    integer, parameter :: FACE_SOUTH = 3
    integer, parameter :: FACE_NORTH = 4
    integer, parameter :: FACE_BOTTOM = 5
    integer, parameter :: FACE_TOP = 6
    
    !> Mesh topology data structure
    type :: mesh_topology_t
        integer :: ncells      ! Number of cells
        integer :: nfaces      ! Number of faces
        integer :: nvertices   ! Number of vertices
        
        ! Cell connectivity (index of neighbors for each cell)
        integer, allocatable :: cell_neighbors(:,:)  ! (6, ncells)
        
        ! Face connectivity
        integer, allocatable :: face_cells(:,:)      ! (2, nfaces) - cells on each side of a face
        
        ! Vertex coordinates
        real(WP), allocatable :: vertices(:,:)       ! (3, nvertices)
        
        ! Cell centers
        real(WP), allocatable :: cell_centers(:,:)   ! (3, ncells)
        
        ! Face centers
        real(WP), allocatable :: face_centers(:,:)   ! (3, nfaces)
        
        ! Face areas
        real(WP), allocatable :: face_areas(:)       ! (nfaces)
        
        ! Face normals
        real(WP), allocatable :: face_normals(:,:)   ! (3, nfaces)
        
        ! Cell volumes
        real(WP), allocatable :: cell_volumes(:)     ! (ncells)
    end type mesh_topology_t
    
    !> Mesh boundary data structure
    type :: boundary_t
        integer :: npatch               ! Number of boundary patches
        integer, allocatable :: patch_start(:)   ! Starting face of each patch
        integer, allocatable :: patch_end(:)     ! Ending face of each patch
        integer, allocatable :: patch_type(:)    ! Type of each patch
        character(len=64), allocatable :: patch_name(:)  ! Name of each patch
    end type boundary_t
    
    !> Complete mesh data structure
    type :: mesh_t
        type(mesh_topology_t) :: topo        ! Mesh topology
        type(boundary_t)      :: boundary    ! Mesh boundaries
        
        ! Auxiliary arrays optimized for GPU computation
        real(WP), allocatable :: cell_dx(:)  ! Cell size in x direction
        real(WP), allocatable :: cell_dy(:)  ! Cell size in y direction
        real(WP), allocatable :: cell_dz(:)  ! Cell size in z direction
        
        ! Block decomposition for GPU
        integer :: blocks_x, blocks_y, blocks_z
        integer :: cells_per_block_x, cells_per_block_y, cells_per_block_z
        integer, allocatable :: block_index(:,:,:)   ! Global cell index of each block
        integer, allocatable :: block_offset(:)      ! Cell offset for each block
        
        ! Domain dimensions
        real(WP) :: xmin, xmax, ymin, ymax, zmin, zmax
        integer :: nx, ny, nz  ! Number of cells in each direction
        
        ! Partitioning for MPI
        integer :: mpi_rank, mpi_size
        integer :: num_local_cells, num_ghost_cells
        integer, allocatable :: local_to_global(:)   ! Map from local to global indices
        integer, allocatable :: global_to_local(:)   ! Map from global to local indices
        
    contains
        procedure :: init => mesh_init
        procedure :: create_cartesian => mesh_create_cartesian
        procedure :: optimize_for_gpu => mesh_optimize_for_gpu
        procedure :: cleanup => mesh_cleanup
    end type mesh_t
    
contains

    !> Initialize mesh
    subroutine mesh_init(this)
        class(mesh_t), intent(inout) :: this
        
        this%mpi_rank = 0
        this%mpi_size = 1
        this%blocks_x = 1
        this%blocks_y = 1
        this%blocks_z = 1
    end subroutine mesh_init
    
    !> Create a Cartesian mesh
    subroutine mesh_create_cartesian(this, nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax)
        class(mesh_t), intent(inout) :: this
        integer, intent(in) :: nx, ny, nz
        real(WP), intent(in) :: xmin, xmax, ymin, ymax, zmin, zmax
        integer :: i, j, k, idx, face_idx
        real(WP) :: dx, dy, dz
        
        ! Store basic mesh information
        this%nx = nx
        this%ny = ny
        this%nz = nz
        this%xmin = xmin
        this%xmax = xmax
        this%ymin = ymin
        this%ymax = ymax
        this%zmin = zmin
        this%zmax = zmax
        
        ! Calculate cell sizes
        dx = (xmax - xmin) / nx
        dy = (ymax - ymin) / ny
        dz = (zmax - zmin) / nz
        
        ! Allocate storage for mesh
        this%topo%ncells = nx * ny * nz
        this%topo%nvertices = (nx+1) * (ny+1) * (nz+1)
        this%topo%nfaces = nx*(ny+1)*(nz+1) + (nx+1)*ny*(nz+1) + (nx+1)*(ny+1)*nz
        
        allocate(this%topo%cell_centers(3, this%topo%ncells))
        allocate(this%topo%cell_volumes(this%topo%ncells))
        allocate(this%topo%cell_neighbors(6, this%topo%ncells))
        allocate(this%topo%vertices(3, this%topo%nvertices))
        allocate(this%topo%face_centers(3, this%topo%nfaces))
        allocate(this%topo%face_areas(this%topo%nfaces))
        allocate(this%topo%face_normals(3, this%topo%nfaces))
        allocate(this%topo%face_cells(2, this%topo%nfaces))
        
        allocate(this%cell_dx(this%topo%ncells))
        allocate(this%cell_dy(this%topo%ncells))
        allocate(this%cell_dz(this%topo%ncells))
        
        ! Initialize cell sizes (uniform for Cartesian grid)
        this%cell_dx = dx
        this%cell_dy = dy
        this%cell_dz = dz
        
        ! Generate vertex coordinates
        idx = 0
        do k = 0, nz
            do j = 0, ny
                do i = 0, nx
                    idx = idx + 1
                    this%topo%vertices(1, idx) = xmin + i * dx
                    this%topo%vertices(2, idx) = ymin + j * dy
                    this%topo%vertices(3, idx) = zmin + k * dz
                end do
            end do
        end do
        
        ! Generate cell centers and volumes
        idx = 0
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    idx = idx + 1
                    this%topo%cell_centers(1, idx) = xmin + (i - 0.5_WP) * dx
                    this%topo%cell_centers(2, idx) = ymin + (j - 0.5_WP) * dy
                    this%topo%cell_centers(3, idx) = zmin + (k - 0.5_WP) * dz
                    this%topo%cell_volumes(idx) = dx * dy * dz
                end do
            end do
        end do
        
        ! Generate face-cell connectivity and face information
        face_idx = 0
        
        ! TODO: Generate face information for all faces (x, y, z planes)
        ! For brevity, this is not fully implemented here
        
        ! Initialize boundary information
        allocate(this%boundary%patch_start(6))  ! Six patches for a box
        allocate(this%boundary%patch_end(6))
        allocate(this%boundary%patch_type(6))
        allocate(this%boundary%patch_name(6))
        
        this%boundary%npatch = 6
        this%boundary%patch_name(1) = "west"
        this%boundary%patch_name(2) = "east"
        this%boundary%patch_name(3) = "south"
        this%boundary%patch_name(4) = "north"
        this%boundary%patch_name(5) = "bottom"
        this%boundary%patch_name(6) = "top"
        
        ! Initialize neighbor information
        ! TODO: Fill in neighbor information for cells
        
    end subroutine mesh_create_cartesian
    
    !> Optimize mesh data for GPU computation
    subroutine mesh_optimize_for_gpu(this, blocks_x, blocks_y, blocks_z)
        class(mesh_t), intent(inout) :: this
        integer, intent(in) :: blocks_x, blocks_y, blocks_z
        integer :: cells_per_block_x, cells_per_block_y, cells_per_block_z
        
        this%blocks_x = blocks_x
        this%blocks_y = blocks_y
        this%blocks_z = blocks_z
        
        ! Calculate cells per block in each direction
        cells_per_block_x = (this%nx + blocks_x - 1) / blocks_x
        cells_per_block_y = (this%ny + blocks_y - 1) / blocks_y
        cells_per_block_z = (this%nz + blocks_z - 1) / blocks_z
        
        this%cells_per_block_x = cells_per_block_x
        this%cells_per_block_y = cells_per_block_y
        this%cells_per_block_z = cells_per_block_z
        
        ! Allocate and fill block index arrays
        allocate(this%block_index(blocks_x, blocks_y, blocks_z))
        allocate(this%block_offset(blocks_x * blocks_y * blocks_z))
        
        ! TODO: Fill in block indexing information
        ! This involves reordering cells for better GPU memory access patterns
        
    end subroutine mesh_optimize_for_gpu
    
    !> Clean up mesh data
    subroutine mesh_cleanup(this)
        class(mesh_t), intent(inout) :: this
        
        ! Free memory
        if (allocated(this%topo%cell_centers)) deallocate(this%topo%cell_centers)
        if (allocated(this%topo%cell_volumes)) deallocate(this%topo%cell_volumes)
        if (allocated(this%topo%cell_neighbors)) deallocate(this%topo%cell_neighbors)
        if (allocated(this%topo%vertices)) deallocate(this%topo%vertices)
        if (allocated(this%topo%face_centers)) deallocate(this%topo%face_centers)
        if (allocated(this%topo%face_areas)) deallocate(this%topo%face_areas)
        if (allocated(this%topo%face_normals)) deallocate(this%topo%face_normals)
        if (allocated(this%topo%face_cells)) deallocate(this%topo%face_cells)
        
        if (allocated(this%cell_dx)) deallocate(this%cell_dx)
        if (allocated(this%cell_dy)) deallocate(this%cell_dy)
        if (allocated(this%cell_dz)) deallocate(this%cell_dz)
        
        if (allocated(this%block_index)) deallocate(this%block_index)
        if (allocated(this%block_offset)) deallocate(this%block_offset)
        
        if (allocated(this%boundary%patch_start)) deallocate(this%boundary%patch_start)
        if (allocated(this%boundary%patch_end)) deallocate(this%boundary%patch_end)
        if (allocated(this%boundary%patch_type)) deallocate(this%boundary%patch_type)
        if (allocated(this%boundary%patch_name)) deallocate(this%boundary%patch_name)
        
        if (allocated(this%local_to_global)) deallocate(this%local_to_global)
        if (allocated(this%global_to_local)) deallocate(this%global_to_local)
    end subroutine mesh_cleanup

end module mesh 