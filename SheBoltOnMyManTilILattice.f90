program SheBoltOnMyManTilILattice
    use lattice_boltzmann
    implicit none

    integer :: total_images, args_count, io_status, interval_len, interval_num
    integer :: current_img, current_img_x, current_img_y, img_grid_x, img_grid_y
    integer :: norm_width, norm_height, mod_width, mod_height
    character(len=20) :: width_arg, height_arg, coarray_dim_arg, interval_length_arg, num_intervals_arg

    ! Lattice arrays !
    real(8), allocatable :: lattice_initial(:,:,:)[:,:]

    total_images = num_images()
    current_img = this_image()

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

    img_grid_x = coarray_dimensions
    img_grid_y = (total_images + img_grid_x - 1) / img_grid_x

    current_img_x = mod(current_img - 1, img_grid_x) + 1
    current_img_y = (current_img - 1) / img_grid_x + 1

    norm_width = global_width / img_grid_x
    mod_width = mod(global_width, img_grid_x)

    if (current_img_x == img_grid_x) then
        instance_width = norm_width + mod_width + 2
    else
        instance_width = norm_width + 2
    end if

    norm_height = global_height / img_grid_y
    mod_height = mod(global_height, img_grid_y)

    if (current_img_y == img_grid_y) then
        instance_height = norm_height + mod_height + 2
    else
        instance_height = norm_height + 2
    end if

    allocate(lattice_initial(directions, instance_width, instance_height)[coarray_dimensions,*])

    ! Intial lattice !
    !call populate_lattice_random(lattice_initial)
    !call populate_lattice_dense_center(lattice_initial)
    !call populate_lattice_shear_wave(lattice_initial)
    !call populate_lattice_couette(lattice_initial)
    !call populate_lattice_poiseuille(lattice_initial)
    !call populate_lattice_sliding_lid(lattice_initial)
    call populate_lattice_sliding_lid_parallel(lattice_initial)
    
    sync all

    !call do_shear_wave_decay_test( lattice_initial)
    call do_parallel_performance_test(lattice_initial)
    !call output_results(lattice_initial, interval_len, interval_num)

    contains

        subroutine output_results(lattice, interval_length, num_intervals)
            real(8), intent(inout) :: lattice(directions, instance_width, instance_height)[coarray_dimensions,*]
            integer, intent(in) :: interval_length, num_intervals
            integer :: interval, sub_interval

            ! Initial lattice output (step 0)
            call gather_and_write(lattice, 0)

            do interval=1, num_intervals
                do sub_interval=1, interval_length
                    call perform_one_time_step(lattice)
                end do

                ! CRITICAL: Ensure all images finish calculating before image 1 starts reading
                sync all

                call gather_and_write(lattice, interval * interval_length)
            end do
        end subroutine output_results

        ! Sibling subroutine to handle the pulling and writing
        subroutine gather_and_write(lattice, step_num)
            real(8), intent(in) :: lattice(directions, instance_width, instance_height)[coarray_dimensions,*]
            integer, intent(in) :: step_num
            
            ! Global arrays that only image 1 needs to construct the full field
            real(8), allocatable :: global_density(:,:)
            type(velocity), allocatable :: global_velocity(:,:)
            
            integer :: img, remote_img, rx, ry, rx_start, ry_start, rw, rh, j, k
            character(:), allocatable :: file_name
            character(len=20) :: interval_str
            
            ! Buffers scaled to your instance bounds to match your module functions
            real(8) :: remote_lattice(directions, instance_width, instance_height)
            real(8) :: local_density(instance_width, instance_height)
            type(velocity) :: local_velocity(instance_width, instance_height)
            
            img = this_image()
            
            if (img == 1) then
                allocate(global_density(global_width, global_height))
                allocate(global_velocity(global_width, global_height))

                ! Assemble the decomposed domain by looping over all images
                do remote_img = 1, total_images
                    ! Determine the remote image's location in the coarray grid
                    rx = mod(remote_img - 1, img_grid_x) + 1
                    ry = (remote_img - 1) / img_grid_x + 1
                    
                    ! Calculate local bounds for this specific chunk (ignoring ghost nodes)
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
                    
                    ! Map local chunk to global starting coordinates
                    rx_start = (rx - 1) * norm_width + 1
                    ry_start = (ry - 1) * norm_height + 1
                    
                    ! Pull the lattice data from the remote image
                    remote_lattice = lattice(:,:,:)[rx, ry]
                    
                    ! Calculate macroscopic variables using module functions
                    local_density = calculate_density_array(remote_lattice)
                    local_velocity = calculate_average_velocity_array(remote_lattice, local_density)
                    
                    ! Stitch the valid data (ignoring ghost boundaries 1 and N) into the global field
                    do j = 1, rw
                        do k = 1, rh
                            global_density(rx_start + j - 1, ry_start + k - 1) = local_density(j + 1, k + 1)
                            global_velocity(rx_start + j - 1, ry_start + k - 1) = local_velocity(j + 1, k + 1)
                        end do
                    end do
                end do
                
                ! Write the fully assembled domain to file
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

        subroutine do_shear_wave_decay_test(lattice)
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

        subroutine do_parallel_performance_test(lattice)
            ! Changed to intent(inout) because the simulation modifies it
            real(8), intent(inout) :: lattice(directions, instance_width, instance_height)[coarray_dimensions,*]

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
                call perform_one_time_step(lattice)
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
        
end program SheBoltOnMyManTilILattice