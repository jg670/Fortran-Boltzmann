program LatticeBoltzmannMain
    use lattice_boltzmann
    implicit none

    integer :: total_images, args_count, io_status, interval_len, interval_num
    integer :: current_img, current_img_x, current_img_y, img_grid_x, img_grid_y
    integer :: norm_width, norm_height, mod_width, mod_height
    character(len=20) :: width_arg, height_arg, coarray_dim_arg, interval_length_arg, num_intervals_arg

    ! Lattice arrays !
    real(8), allocatable :: lattice_even(:,:,:)[:,:], lattice_odd(:,:,:)[:,:]

    ! Get info about image structure !
    total_images = num_images()
    current_img = this_image()


    !======================================================!
    !                   Argument Handling                  !
    !======================================================!

    args_count = command_argument_count()

    if (args_count < 5) then
        print *, "Missing arguments, should be [width height img_dim interval_length num_intervals]"
        stop
    end if

    call get_command_argument(1, width_arg)
    call get_command_argument(2, height_arg)
    call get_command_argument(3, coarray_dim_arg)
    call get_command_argument(4, interval_length_arg)
    call get_command_argument(5, num_intervals_arg)

    read(width_arg, *, iostat=io_status) global_width
    if (io_status /= 0) then
        print *, "Argument '", trim(width_arg), "' is not a valid integer."
        stop
    end if

    read(height_arg, *, iostat=io_status) global_height
    if (io_status /= 0) then
        print *, "Argument '", trim(height_arg), "' is not a valid integer."
        stop
    end if

    read(coarray_dim_arg, *, iostat=io_status) coarray_dimensions
    if (io_status /= 0) then
        print *, "Argument '", trim(coarray_dim_arg), "' is not a valid integer."
        stop
    end if

    read(interval_length_arg, *, iostat=io_status) interval_len
    if (io_status /= 0) then
        print *, "Argument '", trim(interval_length_arg), "' is not a valid integer."
        stop
    end if

    read(num_intervals_arg, *, iostat=io_status) interval_num
    if (io_status /= 0) then
        print *, "Argument '", trim(num_intervals_arg), "' is not a valid integer."
        stop
    end if


    !======================================================!
    !           Splitting Nodes Amongst Images             !
    !======================================================!

    ! Image grid dimensions !
    img_grid_x = coarray_dimensions
    img_grid_y = (total_images + img_grid_x - 1) / img_grid_x

    ! Current image coordinates !
    current_img_x = mod(current_img - 1, img_grid_x) + 1
    current_img_y = (current_img - 1) / img_grid_x + 1

    ! Determine width of image (and get remainder if not evenly divisible) !
    norm_width = global_width / img_grid_x
    mod_width = mod(global_width, img_grid_x)

    ! If current image is on the border, then add remainder to its width !
    if (current_img_x == img_grid_x) then
        instance_width = norm_width + mod_width + 2
    else
        instance_width = norm_width + 2
    end if

    ! Do the same operation to height !
    norm_height = global_height / img_grid_y
    mod_height = mod(global_height, img_grid_y)

    if (current_img_y == img_grid_y) then
        instance_height = norm_height + mod_height + 2
    else
        instance_height = norm_height + 2
    end if

    ! Allocate local image arrays !
    allocate(lattice_even(directions, instance_width, instance_height)[coarray_dimensions,*])
    allocate(lattice_odd(directions, instance_width, instance_height)[coarray_dimensions,*])

    
    !======================================================!
    !                  Run the simulation                  !
    !======================================================!

    !!! Uncomment which initial conditions to use for simulation, USE ONLY ONE !!!

    ! IMPORTANT NOTE: The following simulations do not support multiple images !
    !call populate_lattice_random(lattice_initial)
    !call populate_lattice_dense_center(lattice_initial)
    !call populate_lattice_shear_wave(lattice_initial)
    !call populate_lattice_couette(lattice_initial)
    !call populate_lattice_poiseuille(lattice_initial)
    !call populate_lattice_sliding_lid(lattice_initial)

    ! This simulation is the only one that supports using multiple images !
    call populate_lattice_sliding_lid_parallel(lattice_even)

    lattice_odd = lattice_even
    
    sync all

    !!! Uncomment which test to perform (if any) !!!

    ! Runs simulation with different omega values [0.2 -> 1.8] for 100 time steps, outputting the x velocity every 10 time steps along with the step number !
    ! To be used in tandem with the code in visualization jupyter notebook to compare with the analytical solution !
    ! IMPORTANT NOTE: This test uses hardcoded values and only works with a 30 x 30 grid and one image !
    !call do_shear_wave_decay_test(lattice_initial)

    ! Runs the simulation for 1000 time steps without any file I/O and calculates the MLUPS, printing the results !
    !call do_parallel_performance_test(lattice_even, lattice_odd)

    ! Runs the simulation for interval_len * interval_num time steps, outputting the lattice state every lattice_len steps interval_num times !
    call output_results(lattice_even, lattice_odd, interval_len, interval_num)

    contains

        ! Runs the simulation for interval_len * interval_num time steps, outputting the lattice state every lattice_len steps interval_num times !
        subroutine output_results(lattice_even, lattice_odd, interval_length, num_intervals)

            ! Input variable declarations !
            real(8), intent(inout) :: lattice_even(directions, instance_width, instance_height)[coarray_dimensions,*]
            real(8), intent(inout) :: lattice_odd(directions, instance_width, instance_height)[coarray_dimensions,*]
            integer, intent(in) :: interval_length, num_intervals

            ! Iterators !
            integer :: interval, sub_interval
            integer(8) :: total_steps

            ! Initial lattice output (step 0) !
            call gather_and_write(lattice_even, 0)

            total_steps = 0

            ! Interval handling loop !
            do interval=1, num_intervals
                do sub_interval=1, interval_length
                    total_steps = total_steps + 1
                    
                    ! If total_steps is odd, use lattice_even as input data and lattice_odd as output lattice, if total_steps is even then vice versa !
                    if (mod(total_steps, 2) /= 0) then
                        call perform_one_time_step_fast(lattice_even, lattice_odd)
                    else 
                        call perform_one_time_step_fast(lattice_odd, lattice_even)
                    end if
                end do

                ! Ensure all images finish calculating before image 1 starts reading !
                sync all

                ! If total_steps is odd, output lattice_odd, otherwise output lattice_even !
                if (mod(total_steps, 2) /= 0) then
                    call gather_and_write(lattice_odd, interval * interval_length)
                else 
                    call gather_and_write(lattice_even, interval * interval_length)
                end if

            end do
        end subroutine output_results

        ! Helper subroutine to collect and write lattice data to output file !
        subroutine gather_and_write(lattice, step_num)

            ! Input variable declarations !
            real(8), intent(in) :: lattice(directions, instance_width, instance_height)[coarray_dimensions,*]
            integer, intent(in) :: step_num
            
            ! Global lattice aggregates !
            real(8), allocatable :: global_density(:,:)
            type(velocity), allocatable :: global_velocity(:,:)
            
            ! Image variables !
            integer :: img, remote_img, rx, ry, rx_start, ry_start, rw, rh
            
            ! Iterators !
            integer :: j, k

            ! File I/O variables !
            character(:), allocatable :: file_name
            character(len=20) :: interval_str
            
            ! Local image lattices !
            real(8) :: remote_lattice(directions, instance_width, instance_height)
            real(8) :: local_density(instance_width, instance_height)
            type(velocity) :: local_velocity(instance_width, instance_height)
            
            img = this_image()
            
            ! Only run on first image (prevents race conditions / duplicate outputs) !
            if (img == 1) then
                allocate(global_density(global_width, global_height))
                allocate(global_velocity(global_width, global_height))

                ! Loop over all images !
                do remote_img = 1, total_images

                    ! Determine image position in image grid !
                    rx = mod(remote_img - 1, img_grid_x) + 1
                    ry = (remote_img - 1) / img_grid_x + 1
                    
                    ! Figure out size of the lattice on this image !
                    if (rx == img_grid_x) then
                        rw = norm_width + mod_width
                    else
                        rw = norm_width
                    end if
                    
                    if (ry == img_grid_y) then
                        rh = norm_height + mod_height
                    else
                        rh = norm_height
                    end if
                    
                    ! Map (image) local lattice coords onto global lattice !
                    rx_start = (rx - 1) * norm_width + 1
                    ry_start = (ry - 1) * norm_height + 1
                    
                    ! Get the lattice from the image !
                    remote_lattice = lattice(:,:,:)[rx, ry]
                    
                    ! Calculate density and velocity of image's lattice !
                    local_density = calculate_density_array(remote_lattice)
                    local_velocity = calculate_average_velocity_array(remote_lattice, local_density)
                    
                    ! Add image data to the global array !
                    do j = 1, rw
                        do k = 1, rh
                            global_density(rx_start + j - 1, ry_start + k - 1) = local_density(j + 1, k + 1)
                            global_velocity(rx_start + j - 1, ry_start + k - 1) = local_velocity(j + 1, k + 1)
                        end do
                    end do
                end do
                
                ! Write the global arrays to file (NOTE: files are currently appended with "-sliding-lid", change this for different simulations so the python code knows which files to use) !
                write(interval_str, '(I0)') step_num
                file_name = "./Visualization/output-" // trim(adjustl(interval_str)) // "-sliding-lid.txt"
                open(1, file=file_name, status="replace", action="write")
                do j = 1, global_width
                    do k = 1, global_height
                        write(1, *) j, ", ", k, ", ", global_density(j,k), ", ", global_velocity(j,k)%x, ", ", global_velocity(j,k)%y
                    end do
                end do
                close(1)
                
                deallocate(global_density, global_velocity)
            end if
        end subroutine gather_and_write

        ! Runs simulation with different omega values [0.2 -> 1.8] for 100 time steps, outputting the x velocity every 10 time steps along with the step number !
        ! To be used in tandem with the code in visualization jupyter notebook to compare with the analytical solution !
        ! IMPORTANT NOTE: This test uses hardcoded values and only works with a 30 x 30 grid and one image !
        subroutine do_shear_wave_decay_test(lattice)

            ! Input variable declarations !
            real(8), intent(inout) :: lattice(directions, instance_width, instance_height)[*]
            
            real(8), dimension(9) :: omega_vals = [0.2_8, 0.4_8, 0.6_8, 0.8_8, 1.0_8, 1.2_8, 1.4_8, 1.6_8, 1.8_8]
            integer :: o_idx, i, d
            real(8) :: point_density, point_momentum_x
            real(8) :: target_u
            
            integer :: decay_unit = 25
            character(len=50) :: filename
            integer(8) :: total_steps = 1000
            
            ! Target analytical global coordinate
            integer :: target_global_y = 8 
            
            ! Array indices accounting for a 1-cell thick ghost boundary
            integer :: target_array_y, target_array_x
            
            target_array_y = target_global_y + 1
            target_array_x = 2 ! Global x = 1

            do o_idx = 1, size(omega_vals)
                omega = omega_vals(o_idx)
                call populate_lattice_shear_wave(lattice)

                ! No need to check for current_img == 1
                write(filename, '("./visualization/shear_decay_omega_", F3.1, ".txt")') omega
                open(unit=decay_unit, file=trim(filename), status="replace")
                write(decay_unit, *) "TimeStep Velocity_X"

                do i = 1, total_steps
                    call perform_one_time_step(lattice)
                    
                    target_u = 0.0_8
                    point_density = 0.0_8
                    point_momentum_x = 0.0_8
                    
                    ! Calculate moments directly at the specific node
                    do d = 1, directions
                        point_density = point_density + lattice(d, target_array_x, target_array_y)
                        point_momentum_x = point_momentum_x + lattice(d, target_array_x, target_array_y) * shift_directions_x(d)
                    end do
                    
                    if (point_density > 0.0_8) then
                        target_u = abs(point_momentum_x / point_density)
                    end if

                    ! No reduction (co_max) needed, just directly write the output
                    if (mod(i, 10) == 0) then 
                        write(decay_unit, *) i, target_u
                    end if
                end do

                close(decay_unit)
                print *, "Completed analytical decay test for omega = ", omega
            end do
            
        end subroutine do_shear_wave_decay_test

        subroutine do_parallel_performance_test(lattice_even, lattice_odd)

            real(8), intent(inout) :: lattice_even(directions, instance_width, instance_height)[coarray_dimensions,*]
            real(8), intent(inout) :: lattice_odd(directions, instance_width, instance_height)[coarray_dimensions,*]


            integer :: start_time, end_time, rate
            real(8) :: t_elapsed, mlups
            integer(8) :: total_nodes, total_steps
            integer :: img, step_idx ! Declared step_idx

            total_steps = 1000

            ! Total fluid nodes is exactly the global physical domain
            total_nodes = int(global_width, 8) * int(global_height, 8) 

            ! Start timer !
            call system_clock(count=start_time, count_rate=rate)

            ! Main simulation loop !
            do step_idx = 1, total_steps
                ! Pass the dummy argument 'lattice', not 'lattice_initial'
                if (mod(step_idx, 2) /= 0) then
                    call perform_one_time_step_fast(lattice_even, lattice_odd)
                else
                    call perform_one_time_step_fast(lattice_odd, lattice_even)
                end if
            end do

            ! Stop timer !
            call system_clock(count=end_time, count_rate=rate)

            ! Calculate elapsed time and MLUPS (only on image 1 to avoid duplicate printing) !
            img = this_image()
            if (img == 1) then
                t_elapsed = real(end_time - start_time, 8) / real(rate, 8)
                mlups = (real(total_nodes, 8) * real(total_steps, 8)) / (t_elapsed * 1.0e6_8)
                
                print *, "Total time (s): ", t_elapsed
                print *, "MLUPS: ", mlups
            end if

            sync all

        end subroutine do_parallel_performance_test
        
end program LatticeBoltzmannMain