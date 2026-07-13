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
            
            ! Contiguous buffers for X-direction
            real(8) :: send_w(directions, instance_height), recv_w(directions, instance_height)
            real(8) :: send_e(directions, instance_height), recv_e(directions, instance_height)
            integer :: d, k
            
            req_count = 0

            ! Create buffers on the GPU
            !$acc data create(send_w, recv_w, send_e, recv_e)

            ! === PACK X-BUFFERS ON GPU ===
            !$acc parallel loop collapse(2) present(lattice, send_w, send_e)
            do k = 1, instance_height
                do d = 1, directions
                    send_w(d, k) = lattice(d, 2, k)
                    send_e(d, k) = lattice(d, instance_width - 1, k)
                end do
            end do

            ! Uncomment when running on the cluster and comment out flags below !$acc host_data use_device(lattice, send_w, recv_w, send_e, recv_e)
            !$acc update host(send_w, send_e)
            !$acc update host(lattice(1:directions, 1:instance_width, 2:2))
            !$acc update host(lattice(1:directions, 1:instance_width, instance_height-1:instance_height-1))

            ! === PHASE 1: X-Direction Sync ===
            if (n_w /= MPI_PROC_NULL) then
                req_count = req_count + 1
                call MPI_Isend(send_w, directions*instance_height, MPI_DOUBLE_PRECISION, n_w, 0, comm, req(req_count), ierr)
                req_count = req_count + 1
                call MPI_Irecv(recv_w, directions*instance_height, MPI_DOUBLE_PRECISION, n_w, 1, comm, req(req_count), ierr)
            end if

            if (n_e /= MPI_PROC_NULL) then
                req_count = req_count + 1
                call MPI_Isend(send_e, directions*instance_height, MPI_DOUBLE_PRECISION, n_e, 1, comm, req(req_count), ierr)
                req_count = req_count + 1
                call MPI_Irecv(recv_e, directions*instance_height, MPI_DOUBLE_PRECISION, n_e, 0, comm, req(req_count), ierr)
            end if

            !if (req_count > 0) call MPI_Waitall(req_count, req(1:req_count), stat, ierr)

            ! === PHASE 2: Y-Direction Sync (Contiguous, no buffers needed) ===
            !req_count = 0

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
            
            !!$acc end host_data
            !$acc update device(recv_w, recv_e)
            !$acc update device(lattice(1:directions, 1:instance_width, 1:1))
            !$acc update device(lattice(1:directions, 1:instance_width, instance_height:instance_height))

            ! === UNPACK X-BUFFERS ON GPU ===
            !$acc parallel loop collapse(2) present(lattice, recv_w, recv_e)
            do k = 1, instance_height
                do d = 1, directions
                    if (n_w /= MPI_PROC_NULL) lattice(d, 1, k) = recv_w(d, k)
                    if (n_e /= MPI_PROC_NULL) lattice(d, instance_width, k) = recv_e(d, k)
                end do
            end do

            !$acc end data
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
                    f2 = lattice_in(3, j, k + 1)     ! x: 0, y: 1   -> Pull South
                    f3 = lattice_in(4, j + 1, k)     ! x:-1, y: 0   -> Pull East
                    f4 = lattice_in(5, j, k - 1)     ! x: 0, y:-1   -> Pull North
                    f5 = lattice_in(6, j - 1, k + 1) ! x: 1, y: 1   -> Pull South-West
                    f6 = lattice_in(7, j + 1, k + 1) ! x:-1, y: 1   -> Pull South-East
                    f7 = lattice_in(8, j + 1, k - 1) ! x:-1, y:-1   -> Pull North-East
                    f8 = lattice_in(9, j - 1, k - 1) ! x: 1, y:-1   -> Pull North-West
                    
                    ! 2. MACROSCOPIC VARIABLES
                    density = f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8
                    vel_x = (f1 - f3 + f5 - f6 - f7 + f8) / density
                    vel_y = (f2 - f4 + f5 + f6 - f7 - f8) / density
                    u_sq = vel_x**2 + vel_y**2

                    ! 3. COLLISION
                    lattice_out(1, j, k) = f0 - omega * (f0 - weights(1) * density * (1.0_8 - 1.5_8 * u_sq))

                    c_dot_u = vel_x
                    lattice_out(2, j, k) = f1 - omega * (f1 - weights(2) * density * (1.0_8 + 3.0_8 * c_dot_u + 4.5_8 * c_dot_u**2 - 1.5_8 * u_sq))
                    
                    c_dot_u = vel_y
                    lattice_out(3, j, k) = f2 - omega * (f2 - weights(3) * density * (1.0_8 + 3.0_8 * c_dot_u + 4.5_8 * c_dot_u**2 - 1.5_8 * u_sq))
                    
                    c_dot_u = -vel_x
                    lattice_out(4, j, k) = f3 - omega * (f3 - weights(4) * density * (1.0_8 + 3.0_8 * c_dot_u + 4.5_8 * c_dot_u**2 - 1.5_8 * u_sq))
                    
                    c_dot_u = -vel_y
                    lattice_out(5, j, k) = f4 - omega * (f4 - weights(5) * density * (1.0_8 + 3.0_8 * c_dot_u + 4.5_8 * c_dot_u**2 - 1.5_8 * u_sq))
                    
                    c_dot_u = vel_x + vel_y
                    lattice_out(6, j, k) = f5 - omega * (f5 - weights(6) * density * (1.0_8 + 3.0_8 * c_dot_u + 4.5_8 * c_dot_u**2 - 1.5_8 * u_sq))
                    
                    c_dot_u = -vel_x + vel_y
                    lattice_out(7, j, k) = f6 - omega * (f6 - weights(7) * density * (1.0_8 + 3.0_8 * c_dot_u + 4.5_8 * c_dot_u**2 - 1.5_8 * u_sq))
                    
                    c_dot_u = -vel_x - vel_y
                    lattice_out(8, j, k) = f7 - omega * (f7 - weights(8) * density * (1.0_8 + 3.0_8 * c_dot_u + 4.5_8 * c_dot_u**2 - 1.5_8 * u_sq))
                    
                    c_dot_u = vel_x - vel_y
                    lattice_out(9, j, k) = f8 - omega * (f8 - weights(9) * density * (1.0_8 + 3.0_8 * c_dot_u + 4.5_8 * c_dot_u**2 - 1.5_8 * u_sq))
                    
                end do
            end do
        end subroutine stream_and_collide

        subroutine sliding_lid_boundary_parallel(lattice, coords, dims)
            real(8), intent(inout) :: lattice(directions, instance_width, instance_height)
            integer, intent(in) :: coords(2), dims(2)

            integer :: j

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
                ! Direction 5
                !$acc parallel loop present(lattice)
                do j = 1, instance_width
                    lattice(5, j, 1) = lattice(3, j, 2)
                end do
                
               ! Direction 8 (Pulls from j, writes to j+1)
                !$acc parallel loop present(lattice)
                do j = 1, instance_width - 1
                    lattice(8, j+1, 1) = lattice(6, j, 2) - 2.0_8 * weights(6) * &
                        (lattice(1, j, 2) + lattice(2, j, 2) + lattice(3, j, 2) + &
                         lattice(4, j, 2) + lattice(5, j, 2) + lattice(6, j, 2) + &
                         lattice(7, j, 2) + lattice(8, j, 2) + lattice(9, j, 2)) * &
                        ((shift_directions_x(6) * wall_speed%x) / (1.0_8 / 3.0_8))
                end do
                    
                ! Direction 9 (Pulls from j+1, writes to j)
                !$acc parallel loop present(lattice)
                do j = 1, instance_width - 1
                    lattice(9, j, 1) = lattice(7, j+1, 2) - 2.0_8 * weights(7) * &
                        (lattice(1, j+1, 2) + lattice(2, j+1, 2) + lattice(3, j+1, 2) + &
                         lattice(4, j+1, 2) + lattice(5, j+1, 2) + lattice(6, j+1, 2) + &
                         lattice(7, j+1, 2) + lattice(8, j+1, 2) + lattice(9, j+1, 2)) * &
                        ((shift_directions_x(7) * wall_speed%x) / (1.0_8 / 3.0_8))
                end do
            end if

        end subroutine sliding_lid_boundary_parallel

        subroutine populate_lattice_sliding_lid_parallel(lattice)

            real(8), intent(inout) :: lattice(directions, instance_width, instance_height)

            real(8) :: density_arr(instance_width, instance_height)
            type(velocity) :: velocity_arr(instance_width, instance_height)

            density_arr = 1.0_8
            velocity_arr = velocity(0.0_8, 0.0_8)

            call set_lid_velocity_given_reynolds(1000.0_8)

            lattice = calculate_equilibrium(density_arr, velocity_arr)

        end subroutine populate_lattice_sliding_lid_parallel

        subroutine set_lid_velocity_given_reynolds(reynolds_number)

            real(8), intent(in) ::  reynolds_number
            
            wall_speed = velocity((reynolds_number * calculate_analytical_viscosity()) / (coarray_dimensions * (instance_width - 2.0_8)), 0.0_8)

        end subroutine set_lid_velocity_given_reynolds

        function calculate_equilibrium(density_arr, velocity_arr)

            real(8), intent(in) :: density_arr(instance_width, instance_height)
            type(velocity), intent(in) :: velocity_arr(instance_width, instance_height)

            real(8) :: calculate_equilibrium(directions, instance_width, instance_height)

            integer :: i
            real(8) :: u_squared(instance_width, instance_height), velocity_sum(instance_width, instance_height)

            calculate_equilibrium = 0.0_8

            u_squared(2:instance_width-1, 2:instance_height-1) = velocity_arr(2:instance_width-1, 2:instance_height-1)%x**2 + velocity_arr(2:instance_width-1, 2:instance_height-1)%y**2

            do i=1, directions
                velocity_sum(2:instance_width-1, 2:instance_height-1) = shift_directions_x(i) * velocity_arr(2:instance_width-1, 2:instance_height-1)%x + shift_directions_y(i) * velocity_arr(2:instance_width-1, 2:instance_height-1)%y

                calculate_equilibrium(i,2:instance_width-1, 2:instance_height-1) = weights(i) * density_arr(2:instance_width-1, 2:instance_height-1) * (1.0_8 + 3.0_8 * velocity_sum(2:instance_width-1, 2:instance_height-1) + (4.5_8 * velocity_sum(2:instance_width-1, 2:instance_height-1)**2) - 1.5_8 * u_squared(2:instance_width-1, 2:instance_height-1))
            end do

        end function calculate_equilibrium

        function calculate_analytical_viscosity()

            real(8) :: calculate_analytical_viscosity
            
            calculate_analytical_viscosity = (1.0_8 / 3.0_8) * ((1.0_8 / omega) - 0.5_8)

        end function calculate_analytical_viscosity

end module lattice_boltzmann_cuda

program LatticeBoltzmannCudaMain
    use lattice_boltzmann_cuda
    use mpi
    implicit none

    integer :: ierr, rank, comm_size
    integer :: dims(2), coords(2)
    integer :: comm_cart
    integer :: neighbor_n, neighbor_s, neighbor_e, neighbor_w
    integer :: args_count, io_status, interval_len, interval_num
    integer :: norm_width, norm_height, mod_width, mod_height
    character(len=20) :: width_arg, height_arg, coarray_dim_arg, interval_length_arg, num_intervals_arg
    logical :: periods(2)

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
    periods = [.false., .false.] ! Non-periodic globally by default

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

    call output_results(lattice_even, lattice_odd, interval_len, interval_num, comm_cart, neighbor_n, neighbor_s, neighbor_e, neighbor_w, coords, dims)

    ! Clean up GPU memory
    !$acc exit data copyout(lattice_even)
    !$acc exit data delete(lattice_odd)

    call MPI_Finalize(ierr)

    contains

        subroutine output_results(lattice_even, lattice_odd, interval_length, num_intervals, comm_cart, n_n, n_s, n_e, n_w, coords, dims)
            use mpi
            real(8), intent(inout) :: lattice_even(directions, instance_width, instance_height)
            real(8), intent(inout) :: lattice_odd(directions, instance_width, instance_height)
            integer, intent(in) :: interval_length, num_intervals
            integer, intent(in) :: comm_cart, n_n, n_s, n_e, n_w, coords(2), dims(2)

            integer :: interval, sub_interval, ierr
            integer(8) :: total_steps

            ! Initial lattice output (step 0)
            call gather_and_write(lattice_even, 0, comm_cart, dims)

            total_steps = 0

            do interval = 1, num_intervals
                do sub_interval = 1, interval_length
                    total_steps = total_steps + 1
                    
                    if (mod(total_steps, 2) /= 0) then
                        call perform_one_time_step_fast(lattice_even, lattice_odd, comm_cart, n_n, n_s, n_e, n_w, coords, dims)
                    else 
                        call perform_one_time_step_fast(lattice_odd, lattice_even, comm_cart, n_n, n_s, n_e, n_w, coords, dims)
                    end if
                end do

                ! CRITICAL: Wait for GPU queues to empty and sync MPI before gathering
                !$acc wait
                call MPI_Barrier(comm_cart, ierr)

                if (mod(total_steps, 2) /= 0) then
                    call gather_and_write(lattice_odd, interval * interval_length, comm_cart, dims)
                else 
                    call gather_and_write(lattice_even, interval * interval_length, comm_cart, dims)
                end if
            end do
        end subroutine output_results


        subroutine gather_and_write(lattice, step_num, comm, dims)
            use mpi
            real(8), intent(inout) :: lattice(directions, instance_width, instance_height)
            integer, intent(in) :: step_num, comm, dims(2)
            
            real(8), allocatable :: global_density(:,:)
            type(velocity), allocatable :: global_velocity(:,:)
            
            integer :: rank, comm_size, remote_rank, ierr, max_w, max_h
            integer :: remote_coords(2), rx, ry, rw, rh, rx_start, ry_start, j, k
            integer :: stat(MPI_STATUS_SIZE)
            character(len=50) :: file_name
            character(len=20) :: interval_str
            integer :: d
            
            call MPI_Comm_rank(comm, rank, ierr)
            call MPI_Comm_size(comm, comm_size, ierr)

            ! 1. Bring local GPU data down to the CPU Host
            !$acc update host(lattice)
            
            if (rank == 0) then
                allocate(global_density(global_width, global_height))
                allocate(global_velocity(global_width, global_height))
                
                ! Assemble the decomposed domain
                do remote_rank = 0, comm_size - 1
                    
                    ! Get the cartesian coordinates of the remote rank
                    call MPI_Cart_coords(comm, remote_rank, 2, remote_coords, ierr)
                    rx = remote_coords(1) ! 0-indexed
                    ry = remote_coords(2) ! 0-indexed
                    
                    ! Calculate exact local bounds for this specific chunk
                    if (rx == dims(1) - 1) then
                        rw = (global_width / dims(1)) + mod(global_width, dims(1))
                    else
                        rw = global_width / dims(1)
                    end if
                    
                    if (ry == dims(2) - 1) then
                        rh = (global_height / dims(2)) + mod(global_height, dims(2))
                    else
                        rh = global_height / dims(2)
                    end if
                    
                    ! Map local chunk to global starting coordinates
                    rx_start = (rx) * (global_width / dims(1)) + 1
                    ry_start = (ry) * (global_height / dims(2)) + 1
                    
                    ! CRITICAL FIX: Create a dynamically sized block so the receive buffer 
                    ! perfectly matches the exact memory stride of the incoming chunk.
                    block
                        real(8) :: exact_buffer(directions, rw + 2, rh + 2)
                        real(8) :: local_density(rw + 2, rh + 2)
                        type(velocity) :: local_velocity(rw + 2, rh + 2)
                    
                        ! Pull the lattice data via MPI
                        if (remote_rank == 0) then
                            exact_buffer = lattice(:, 1:rw+2, 1:rh+2)
                        else
                            call MPI_Recv(exact_buffer, directions * (rw+2) * (rh+2), MPI_DOUBLE_PRECISION, &
                                          remote_rank, 0, comm, stat, ierr)
                        end if
                        
                        ! Inline Density
                        local_density(1:rw+2, 1:rh+2) = 1.0_8
                        local_density(2:rw+1, 2:rh+1) = sum(exact_buffer(:, 2:rw+1, 2:rh+1), dim=1)
                        
                        ! Inline Velocity
                        local_velocity(1:rw+2, 1:rh+2)%x = 0.0_8
                        local_velocity(1:rw+2, 1:rh+2)%y = 0.0_8
                        
                        do d = 1, directions
                            local_velocity(2:rw+1, 2:rh+1)%x = local_velocity(2:rw+1, 2:rh+1)%x + &
                                exact_buffer(d, 2:rw+1, 2:rh+1) * shift_directions_x(d)
                            local_velocity(2:rw+1, 2:rh+1)%y = local_velocity(2:rw+1, 2:rh+1)%y + &
                                exact_buffer(d, 2:rw+1, 2:rh+1) * shift_directions_y(d)
                        end do
                        
                        local_velocity(2:rw+1, 2:rh+1)%x = local_velocity(2:rw+1, 2:rh+1)%x / local_density(2:rw+1, 2:rh+1)
                        local_velocity(2:rw+1, 2:rh+1)%y = local_velocity(2:rw+1, 2:rh+1)%y / local_density(2:rw+1, 2:rh+1)
                        
                        ! Stitch into the global field
                        do j = 1, rw
                            do k = 1, rh
                                global_density(rx_start + j - 1, ry_start + k - 1) = local_density(j + 1, k + 1)
                                global_velocity(rx_start + j - 1, ry_start + k - 1) = local_velocity(j + 1, k + 1)
                            end do
                        end do
                    end block
                end do
                
                ! Write the fully assembled domain to file
                write(interval_str, '(I0)') step_num
                file_name = "./Visualization/output-" // trim(adjustl(interval_str)) // "-sliding-lid.txt"
                open(1, file=file_name, status="replace", action="write")
                do j = 1, global_width
                    do k = 1, global_height
                        ! CRITICAL FIX: Explicit format string to guarantee clean, identical spacing
                        write(1, '(I0, A, I0, A, G15.7, A, G15.7, A, G15.7)') &
                            j, ", ", k, ", ", global_density(j,k), ", ", global_velocity(j,k)%x, ", ", global_velocity(j,k)%y
                    end do
                end do
                close(1)
                
                deallocate(global_density, global_velocity)
                
            else
                ! Worker Ranks simply send their chunk to Rank 0
                call MPI_Send(lattice, directions * instance_width * instance_height, MPI_DOUBLE_PRECISION, 0, 0, comm, ierr)
            end if
            
        end subroutine gather_and_write

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

end program LatticeBoltzmannCudaMain