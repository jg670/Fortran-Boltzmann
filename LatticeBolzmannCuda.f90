program LatticeBoltzmannCuda
    use lattice_boltzmann_cuda
    use mpi
    implicit none

    integer :: ierr, rank, comm_size
    integer :: dims(2), periods(2), coords(2)
    integer :: comm_cart
    integer :: neighbor_n, neighbor_s, neighbor_e, neighbor_w
    integer :: args_count, io_status, interval_len, interval_num
    integer :: norm_width, norm_height, mod_width, mod_height
    character(len=20) :: width_arg, height_arg, coarray_dim_arg, interval_length_arg, num_intervals_arg

    ! Lattice arrays
    real(8), allocatable :: lattice_even(:,:,:), lattice_odd(:,:,:)

    ! Initialize MPI
    call MPI_Init(ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, comm_size, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

    args_count = command_argument_count()
    if (args_count < 5) then
        if (rank == 0) print *, "Missing arguments: [width height grid_x interval_length num_intervals]"
        call MPI_Finalize(ierr)
        stop
    end if

    ! Parse arguments (assuming parsing succeeds for brevity)
    call get_command_argument(1, width_arg); read(width_arg, *) global_width
    call get_command_argument(2, height_arg); read(height_arg, *) global_height
    call get_command_argument(3, coarray_dim_arg); read(coarray_dim_arg, *) coarray_dimensions
    call get_command_argument(4, interval_length_arg); read(interval_length_arg, *) interval_len
    call get_command_argument(5, num_intervals_arg); read(num_intervals_arg, *) interval_num

    ! 2D Cartesian Topology setup
    dims(1) = coarray_dimensions
    dims(2) = comm_size / coarray_dimensions
    periods = [0, 0] ! Non-periodic globally by default

    call MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, .false., comm_cart, ierr)
    call MPI_Comm_rank(comm_cart, rank, ierr)
    call MPI_Cart_coords(comm_cart, rank, 2, coords, ierr)

    ! Discover neighbors (Shift for X and Y dimensions)
    call MPI_Cart_shift(comm_cart, 0, 1, neighbor_w, neighbor_e, ierr)
    call MPI_Cart_shift(comm_cart, 1, 1, neighbor_n, neighbor_s, ierr)

    ! Domain sizing (coords are 0-indexed)
    norm_width = global_width / dims(1)
    mod_width = mod(global_width, dims(1))
    if (coords(1) == dims(1) - 1) then
        instance_width = norm_width + mod_width + 2
    else
        instance_width = norm_width + 2
    end if

    norm_height = global_height / dims(2)
    mod_height = mod(global_height, dims(2))
    if (coords(2) == dims(2) - 1) then
        instance_height = norm_height + mod_height + 2
    else
        instance_height = norm_height + 2
    end if

    allocate(lattice_even(directions, instance_width, instance_height))
    allocate(lattice_odd(directions, instance_width, instance_height))

    ! Initialize lattice on CPU
    call populate_lattice_sliding_lid_parallel(lattice_even)
    lattice_odd = lattice_even
    
    call MPI_Barrier(comm_cart, ierr)

    ! Map arrays to GPU memory
    !$acc enter data copyin(lattice_even, lattice_odd)

    call output_results(lattice_even, lattice_odd, interval_len, interval_num, comm_cart, neighbor_n, neighbor_s, neighbor_e, neighbor_w)

    ! Clean up GPU memory
    !$acc exit data copyout(lattice_even)
    !$acc exit data delete(lattice_odd)

    call MPI_Finalize(ierr)
    
    ! ... (Keep gather_and_write adapted for MPI_Send/Recv to Rank 0) ...

    contains

        subroutine do_parallel_performance_test(lattice_even, lattice_odd, comm_cart, n_n, n_s, n_e, n_w, coords, dims)
            use mpi
            real(8), intent(inout) :: lattice_even(directions, instance_width, instance_height)
            real(8), intent(inout) :: lattice_odd(directions, instance_width, instance_height)
            integer, intent(in) :: comm_cart, n_n, n_s, n_e, n_w
            integer, intent(in) :: coords(2), dims(2)

            integer :: rank, ierr, step_idx
            real(8) :: start_time, end_time, t_elapsed, mlups
            integer(8) :: total_nodes, total_steps

            call MPI_Comm_rank(comm_cart, rank, ierr)

            total_steps = 1000
            ! Total fluid nodes is exactly the global physical domain
            total_nodes = int(global_width, 8) * int(global_height, 8)

            ! CRITICAL: Ensure all MPI processes are synced and GPU queues are empty before timing
            !$acc wait
            call MPI_Barrier(comm_cart, ierr)
            
            start_time = MPI_Wtime()

            ! Main simulation loop
            do step_idx = 1, total_steps
                if (mod(step_idx, 2) /= 0) then
                    call perform_one_time_step_fast(lattice_even, lattice_odd, comm_cart, n_n, n_s, n_e, n_w, coords, dims)
                else
                    call perform_one_time_step_fast(lattice_odd, lattice_even, comm_cart, n_n, n_s, n_e, n_w, coords, dims)
                end if
            end do

            ! CRITICAL: Wait for the GPU to finish the final step before stopping the clock
            !$acc wait
            call MPI_Barrier(comm_cart, ierr)
            
            end_time = MPI_Wtime()

            ! Calculate elapsed time and MLUPS (only on Rank 0 to avoid duplicate printing)
            if (rank == 0) then
                t_elapsed = end_time - start_time
                mlups = (real(total_nodes, 8) * real(total_steps, 8)) / (t_elapsed * 1.0e6_8)
                
                print *, "Total time (s): ", t_elapsed
                print *, "MLUPS: ", mlups
            end if

        end subroutine do_parallel_performance_test

end program LatticeBoltzmannCuda

module lattice_boltzmann_cuda
    implicit none

    type :: velocity
        real(8) :: x
        real(8) :: y
    end type velocity

    integer :: global_width = 15, global_height = 15, instance_height = 17, instance_width = 17
    integer, parameter :: directions = 9
    real(8) :: omega = 1
    real(8), parameter :: shift_directions_x(directions) = [0.0,  1.0,  0.0, -1.0,  0.0,  1.0, -1.0, -1.0,  1.0]
    real(8), parameter :: shift_directions_y(directions) = [0.0,  0.0,  1.0,  0.0, -1.0,  1.0,  1.0, -1.0, -1.0]
    real(8), parameter :: weights(directions) = [4.0_8/9.0_8, 1.0_8/9.0_8, 1.0_8/9.0_8, 1.0_8/9.0_8, 1.0_8/9.0_8, 1.0_8/36.0_8, 1.0_8/36.0_8, 1.0_8/36.0_8, 1.0_8/36.0_8]

    real(8), parameter :: pi = 4.0_8 * atan(1.0_8)
    real(8), parameter :: epsilon = 0.01
    
    type(velocity) :: wall_speed = velocity(0.1_8, 0.0_8)

    real(8), parameter :: poiseuille_in_pressure = 0.3, poiseuille_out_pressure = 0.29

    integer :: coarray_dimensions = 1

    ! Encoding for boundary configuration. 0 is normal periodic boundaries, 1 is Couette flow, 2 is Poiseuille flow, 3 is sliding lid !
    integer :: boundary_configuration = 0

    contains

        subroutine perform_one_time_step_fast(lattice_in, lattice_out, comm_cart, n_n, n_s, n_e, n_w, coords, dims)
            real(8), intent(inout) :: lattice_in(directions, instance_width, instance_height)
            real(8), intent(inout) :: lattice_out(directions, instance_width, instance_height)
            integer, intent(in) :: comm_cart, n_n, n_s, n_e, n_w
            integer, intent(in) :: coords(2), dims(2)

            ! 1. Update all ghost nodes via CUDA-aware MPI
            call sync_ghost_nodes(lattice_in, comm_cart, n_n, n_s, n_e, n_w)

            ! 2. Apply boundary conditions (Needs to know its place in the global grid)
            call sliding_lid_boundary_parallel(lattice_in, coords, dims)

            ! 3. Fused Kernel (Calculates on the GPU)
            call stream_and_collide(lattice_in, lattice_out)

        end subroutine perform_one_time_step_fast

        subroutine sync_ghost_nodes(lattice, comm, n_n, n_s, n_e, n_w)
            use mpi
            real(8), intent(inout) :: lattice(directions, instance_width, instance_height)
            integer, intent(in) :: comm, n_n, n_s, n_e, n_w
            
            integer :: req(8), req_count
            integer :: stat(MPI_STATUS_SIZE, 8), ierr
            
            req_count = 0

            !$acc host_data use_device(lattice)

            ! ==========================================
            ! PHASE 1: X-Direction (East/West) Sync
            ! ==========================================
            
            ! Send Left Inner (2) to West, Receive from West into Left Ghost (1)
            if (n_w /= MPI_PROC_NULL) then
                req_count = req_count + 1
                call MPI_Isend(lattice(1, 2, 1), directions*instance_height, MPI_DOUBLE_PRECISION, n_w, 0, comm, req(req_count), ierr)
                req_count = req_count + 1
                call MPI_Irecv(lattice(1, 1, 1), directions*instance_height, MPI_DOUBLE_PRECISION, n_w, 1, comm, req(req_count), ierr)
            end if

            ! Send Right Inner (W-1) to East, Receive from East into Right Ghost (W)
            if (n_e /= MPI_PROC_NULL) then
                req_count = req_count + 1
                call MPI_Isend(lattice(1, instance_width - 1, 1), directions*instance_height, MPI_DOUBLE_PRECISION, n_e, 1, comm, req(req_count), ierr)
                req_count = req_count + 1
                call MPI_Irecv(lattice(1, instance_width, 1), directions*instance_height, MPI_DOUBLE_PRECISION, n_e, 0, comm, req(req_count), ierr)
            end if

            if (req_count > 0) call MPI_Waitall(req_count, req(1:req_count), stat, ierr)

            ! ==========================================
            ! PHASE 2: Y-Direction (North/South) Sync
            ! ==========================================
            req_count = 0

            ! Send Bottom Inner (2) to North, Receive from North into Bottom Ghost (1)
            ! (Assuming index 1 is North / Y-axis orientation)
            if (n_n /= MPI_PROC_NULL) then
                req_count = req_count + 1
                call MPI_Isend(lattice(1, 1, 2), directions*instance_width, MPI_DOUBLE_PRECISION, n_n, 2, comm, req(req_count), ierr)
                req_count = req_count + 1
                call MPI_Irecv(lattice(1, 1, 1), directions*instance_width, MPI_DOUBLE_PRECISION, n_n, 3, comm, req(req_count), ierr)
            end if

            ! Send Top Inner (H-1) to South, Receive from South into Top Ghost (H)
            if (n_s /= MPI_PROC_NULL) then
                req_count = req_count + 1
                call MPI_Isend(lattice(1, 1, instance_height - 1), directions*instance_width, MPI_DOUBLE_PRECISION, n_s, 3, comm, req(req_count), ierr)
                req_count = req_count + 1
                call MPI_Irecv(lattice(1, 1, instance_height), directions*instance_width, MPI_DOUBLE_PRECISION, n_s, 2, comm, req(req_count), ierr)
            end if

            if (req_count > 0) call MPI_Waitall(req_count, req(1:req_count), stat, ierr)

            !$acc end host_data
        end subroutine sync_ghost_nodes

        subroutine stream_and_collide(lattice_in, lattice_out)
            real(8), intent(in) :: lattice_in(directions, instance_width, instance_height)
            real(8), intent(inout) :: lattice_out(directions, instance_width, instance_height)
            
            integer :: j, k
            real(8) :: f0, f1, f2, f3, f4, f5, f6, f7, f8
            real(8) :: density, vel_x, vel_y, u_sq, c_dot_u

            ! Execute on GPU, collapsing loops for maximum parallelism
            !$acc parallel loop collapse(2) present(lattice_in, lattice_out) private(f0, f1, f2, f3, f4, f5, f6, f7, f8, density, vel_x, vel_y, u_sq, c_dot_u)
            do k = 2, instance_height - 1
                do j = 2, instance_width - 1
                    
                    ! 1. PULL STREAMING
                    f0 = lattice_in(1, j, k)
                    f1 = lattice_in(2, j - 1, k)
                    ! ... [Keep original f2 through f8 reads] ...
                    
                    ! 2. MACROSCOPIC VARIABLES
                    density = f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8
                    vel_x = (f1 - f3 + f5 - f6 - f7 + f8) / density
                    vel_y = (f2 - f4 + f5 + f6 - f7 - f8) / density
                    u_sq = vel_x**2 + vel_y**2

                    ! 3. COLLISION
                    lattice_out(1, j, k) = f0 - omega * (f0 - weights(1) * density * (1.0_8 - 1.5_8 * u_sq))
                    ! ... [Keep original lattice_out(2) through lattice_out(9) math] ...
                    
                end do
            end do
        end subroutine stream_and_collide

        subroutine sliding_lid_boundary_parallel(lattice, coords, dims)
            real(8), intent(inout) :: lattice(directions, instance_width, instance_height)
            integer, intent(in) :: coords(2), dims(2)

            ! Handle left border (MPI coords are 0-indexed)
            if (coords(1) == 0) then
                !$acc kernels present(lattice)
                lattice(2,1,:) = lattice(4,2,:)
                lattice(6,1,2:instance_height) = lattice(8,2,1:instance_height-1)
                lattice(9,1,1:instance_height-1) = lattice(7,2,2:instance_height)
                !$acc end kernels
            end if

            ! Handle right border
            if (coords(1) == dims(1) - 1) then
                !$acc kernels present(lattice)
                lattice(4,instance_width,:) = lattice(2, instance_width-1,:)
                lattice(7,instance_width,2:instance_height) = lattice(9,instance_width-1,1:instance_height-1)
                lattice(8,instance_width,1:instance_height-1) = lattice(6,instance_width-1,2:instance_height)
                !$acc end kernels
            end if

            ! Handle bottom border (assuming coords(2) == dims(2)-1 is the bottom)
            if (coords(2) == dims(2) - 1) then
                !$acc kernels present(lattice)
                lattice(3,:,instance_height) = lattice(5,:,instance_height-1)
                lattice(6,1:instance_width-1,instance_height) = lattice(8,2:instance_width,instance_height-1)
                lattice(7,2:instance_width,instance_height) = lattice(9,1:instance_width-1,instance_height-1)
                !$acc end kernels
            end if

            ! Handle top border (moving lid)
            if (coords(2) == 0) then
                !$acc kernels present(lattice)
                lattice(5,:,1) = lattice(3,:,2)
                
                lattice(8,2:instance_width,1) = lattice(6,1:instance_width-1,2) - &
                    2.0_8 * weights(6) * sum(lattice(:, 1:instance_width-1, 2), dim=1) * &
                    ((shift_directions_x(6) * wall_speed%x) / (1.0_8 / 3.0_8))
                    
                lattice(9,1:instance_width-1,1) = lattice(7,2:instance_width,2) - &
                    2.0_8 * weights(7) * sum(lattice(:, 2:instance_width, 2), dim=1) * &
                    ((shift_directions_x(7) * wall_speed%x) / (1.0_8 / 3.0_8))
                !$acc end kernels
            end if

        end subroutine sliding_lid_boundary_parallel

end module lattice_boltzmann_cuda