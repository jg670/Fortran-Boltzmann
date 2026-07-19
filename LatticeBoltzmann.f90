module lattice_boltzmann
    implicit none

    ! Type used for velocity vectors !
    type :: velocity
        real(8) :: x
        real(8) :: y
    end type velocity

    ! Lattice dimensions, will be overwritten by command line arguments !
    integer :: global_width = 15, global_height = 15, instance_height = 17, instance_width = 17

    ! Simulation parameters !
    integer, parameter :: directions = 9
    real(8) :: omega = 1
    real(8), parameter :: shift_directions_x(directions) = [0.0,  1.0,  0.0, -1.0,  0.0,  1.0, -1.0, -1.0,  1.0]
    real(8), parameter :: shift_directions_y(directions) = [0.0,  0.0,  1.0,  0.0, -1.0,  1.0,  1.0, -1.0, -1.0]
    real(8), parameter :: weights(directions) = [4.0_8/9.0_8, 1.0_8/9.0_8, 1.0_8/9.0_8, 1.0_8/9.0_8, 1.0_8/9.0_8, 1.0_8/36.0_8, 1.0_8/36.0_8, 1.0_8/36.0_8, 1.0_8/36.0_8]
    real(8), parameter :: pi = 4.0_8 * atan(1.0_8)
    real(8), parameter :: epsilon = 0.01
    
    ! Wall velocity, used when simulation has a moving wall !
    type(velocity) :: wall_speed = velocity(0.1_8, 0.0_8)

    ! Input and output pressure used for Poiseuille flow !
    real(8), parameter :: poiseuille_in_pressure = 0.3, poiseuille_out_pressure = 0.29

    ! Width of coarray grid !
    integer :: coarray_dimensions = 1

    ! Encoding for boundary configuration. 0 is normal periodic boundaries, 1 is Couette flow, 2 is Poiseuille flow, 3 is sliding lid !
    integer :: boundary_configuration = 0

    contains

        subroutine perform_one_time_step(lattice_in, lattice_out)

            real(8), intent(inout) :: lattice_in(directions, instance_width, instance_height)[coarray_dimensions,*]
            real(8), intent(inout) :: lattice_out(directions, instance_width, instance_height)[coarray_dimensions,*]

            ! 1. Update all ghost nodes in the input array so the Pull scheme has valid data [cite: 12]
            call sync_ghost_nodes(lattice_in)

            ! 2. Apply appropriate boundary conditions to ghost nodes [cite: 11, 56, 57]
            if (boundary_configuration == 0) then
                call periodic_boundary(lattice_in)
            else if (boundary_configuration == 1) then
                call couette_boundary(lattice_in)
            else if (boundary_configuration == 2) then
                call poiseuille_boundary(lattice_in)
            else if (boundary_configuration == 3) then
                call sliding_lid_boundary(lattice_in)
            end if

            ! 3. Fused Kernel (Pulls from lattice_in, computes, writes to lattice_out) [cite: 13]
            call stream_and_collide(lattice_in, lattice_out)

        end subroutine perform_one_time_step

        ! Synchronize ghost boundaries across all images
        subroutine sync_ghost_nodes(lattice)
            
            real(8), intent(inout) :: lattice(directions, instance_width, instance_height)[coarray_dimensions,*]
            
            integer :: img(2), total_images, img_grid_y

            img = this_image(lattice)
            total_images = num_images()
            img_grid_y = (total_images + coarray_dimensions - 1) / coarray_dimensions

            ! ==========================================
            ! PHASE 1: X-Direction (East/West) Sync
            ! ==========================================
            
            ! Send left inner edge to Western neighbor's right ghost edge
            if (img(1) /= 1) then
                lattice(:, instance_width, :)[img(1) - 1, img(2)] = lattice(:, 2, :)
            end if

            ! Send right inner edge to Eastern neighbor's left ghost edge
            if (img(1) /= coarray_dimensions) then
                lattice(:, 1, :)[img(1) + 1, img(2)] = lattice(:, instance_width - 1, :)
            end if

            ! CRITICAL: Sync here so X-data is fully received before Y-data is sent.
            ! This ensures diagonal particles correctly propagate into the corner ghost nodes.
            sync all

            ! ==========================================
            ! PHASE 2: Y-Direction (North/South) Sync
            ! ==========================================
            
            ! Send bottom inner edge to Southern neighbor's top ghost edge
            if (img(2) /= img_grid_y) then
                lattice(:, :, 1)[img(1), img(2) + 1] = lattice(:, :, instance_height - 1)
            end if

            ! Send top inner edge to Northern neighbor's bottom ghost edge
            if (img(2) /= 1) then
                lattice(:, :, instance_height)[img(1), img(2) - 1] = lattice(:, :, 2)
            end if

            ! Ensure all ghost nodes (including corners) are ready before streaming
            sync all

        end subroutine sync_ghost_nodes

        subroutine stream_and_collide(lattice_in, lattice_out)

            real(8), intent(in) :: lattice_in(directions, instance_width, instance_height)
            real(8), intent(inout) :: lattice_out(directions, instance_width, instance_height)

            integer :: j, k
            real(8) :: f0, f1, f2, f3, f4, f5, f6, f7, f8
            real(8) :: density, vel_x, vel_y, u_sq, c_dot_u

            do k = 2, instance_height - 1
                do j = 2, instance_width - 1
                    
                    ! 1. PULL STREAMING
                    ! Read directly from neighbors based on your shift_directions arrays
                    f0 = lattice_in(1, j, k)
                    f1 = lattice_in(2, j - 1, k)     ! x: 1, y: 0   -> Pull West
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

                    ! 3. COLLISION (Inline Equilibrium and Relaxation)
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

        subroutine sliding_lid_boundary(lattice)

            real(8), intent(inout) :: lattice(directions, instance_width, instance_height)[coarray_dimensions,*]

            integer :: img(2), total_images, img_grid_y

            img = this_image(lattice)
            total_images = num_images()
            img_grid_y = (total_images + coarray_dimensions - 1) / coarray_dimensions

            ! Handle left border !
            if (img(1) == 1) then
                ! Left side bounce-back boundary (indexes are +1 since fortran is 1-indexed) !
                lattice(2,1,:) = lattice(4,2,:)
                lattice(6,1,2:instance_height) = lattice(8,2,1:instance_height-1)
                lattice(9,1,1:instance_height-1) = lattice(7,2,2:instance_height)
            end if

            ! Handle right border !
            if (img(1) == coarray_dimensions) then
                ! Right side bounce-back boundary (indexes are +1 since fortran is 1-indexed) !
                lattice(4,instance_width,:) = lattice(2, instance_width-1,:)
                lattice(7,instance_width,2:instance_height) = lattice(9,instance_width-1,1:instance_height-1)
                lattice(8,instance_width,1:instance_height-1) = lattice(6,instance_width-1,2:instance_height)
            
            end if

            ! Handle bottom border !
            if (img(2) == img_grid_y) then
                ! Bottom bounce-back boundary (indexes are +1 since fortran is 1-indexed) !
                lattice(3,:,instance_height) = lattice(5,:,instance_height-1)
                lattice(6,1:instance_width-1,instance_height) = lattice(8,2:instance_width,instance_height-1)
                lattice(7,2:instance_width,instance_height) = lattice(9,1:instance_width-1,instance_height-1)

            end if

            ! Handle top border !
            if (img(2) == 1) then
                ! Top moving boundary (y velocity is 0, so nothing is added) !
                lattice(5,:,1) = lattice(3,:,2)
                lattice(8,2:instance_width,1) = lattice(6,1:instance_width-1,2) - 2.0_8 * weights(6) * sum(lattice(:, 1:instance_width-1, 2), dim=1) * ((shift_directions_x(6) * wall_speed%x) / (1.0_8 / 3.0_8))
                lattice(9,1:instance_width-1,1) = lattice(7,2:instance_width,2) - 2.0_8 * weights(7) * sum(lattice(:, 2:instance_width, 2), dim=1) * ((shift_directions_x(7) * wall_speed%x) / (1.0_8 / 3.0_8))

            end if

        end subroutine sliding_lid_boundary

        subroutine periodic_boundary(lattice)

            real(8), intent(inout) :: lattice(directions, instance_width, instance_height)[coarray_dimensions,*]
            
            integer :: img(2), total_images, img_grid_y

            img = this_image(lattice)
            total_images = num_images()
            img_grid_y = (total_images + coarray_dimensions - 1) / coarray_dimensions

            ! ==========================================
            ! PHASE 1: Global X-Direction Wrap-Around
            ! ==========================================
            
            ! If on the West edge of the global grid, pull from the East edge
            if (img(1) == 1) then
                lattice(:, 1, :) = lattice(:, instance_width - 1, :)[coarray_dimensions, img(2)]
            end if

            ! If on the East edge of the global grid, pull from the West edge
            if (img(1) == coarray_dimensions) then
                lattice(:, instance_width, :) = lattice(:, 2, :)[1, img(2)]
            end if

            ! Sync X-data so diagonal particles wrapping around the corners are ready for the Y-pass
            sync all

            ! ==========================================
            ! PHASE 2: Global Y-Direction Wrap-Around
            ! ==========================================
            
            ! If on the North edge of the global grid, pull from the South edge
            if (img(2) == 1) then
                lattice(:, :, 1) = lattice(:, :, instance_height - 1)[img(1), img_grid_y]
            end if

            ! If on the South edge of the global grid, pull from the North edge
            if (img(2) == img_grid_y) then
                lattice(:, :, instance_height) = lattice(:, :, 2)[img(1), 1]
            end if

            ! Final sync to ensure all ghost nodes are populated before the fused kernel runs
            sync all

        end subroutine periodic_boundary

        subroutine couette_boundary(lattice)
            real(8), intent(inout) :: lattice(directions, instance_width, instance_height)[coarray_dimensions,*]
            
            integer :: img(2), total_images, img_grid_y
            real(8) :: average_density

            img = this_image(lattice)
            total_images = num_images()
            img_grid_y = (total_images + coarray_dimensions - 1) / coarray_dimensions
            
            average_density = sum(lattice) / (instance_width * instance_height)

            ! Handle global periodic side walls across coarray boundaries
            if (img(1) == 1) then
                lattice(:,1,:) = lattice(:, instance_width-1, :)[coarray_dimensions, img(2)]
            end if
            if (img(1) == coarray_dimensions) then
                lattice(:,instance_width,:) = lattice(:, 2, :)[1, img(2)]
            end if
            
            sync all ! Ensure cross-image periodic assignments resolve before pulling

            ! Handle bottom bounce-back border [cite: 125, 126]
            if (img(2) == img_grid_y) then
                lattice(3,:,instance_height) = lattice(5,:,instance_height-1)
                lattice(6,1:instance_width-1,instance_height) = lattice(8,2:instance_width,instance_height-1)
                lattice(7,2:instance_width,instance_height) = lattice(9,1:instance_width-1,instance_height-1)
            end if

            ! Handle top moving boundary [cite: 127, 129]
            if (img(2) == 1) then
                lattice(5,:,1) = lattice(3,:,2)
                lattice(8,2:instance_width,1) = lattice(6,1:instance_width-1,2) - 2.0_8 * weights(6) * average_density * ((shift_directions_x(6) * wall_speed%x) / (1.0_8 / 3.0_8))
                lattice(9,1:instance_width-1,1) = lattice(7,2:instance_width,2) - 2.0_8 * weights(7) * average_density * ((shift_directions_x(7) * wall_speed%x) / (1.0_8 / 3.0_8))
            end if
        end subroutine couette_boundary

        subroutine poiseuille_boundary(lattice)
            real(8), intent(inout) :: lattice(directions, instance_width, instance_height)[coarray_dimensions,*]
            
            integer :: img(2), total_images, img_grid_y
            real(8) :: density_in, density_out
            real(8) :: density_arr(instance_width, instance_height)
            type(velocity) :: velocity_arr(instance_width, instance_height)
            real(8) :: f_star_east_minus_feq(directions, 1, instance_height), f_star_west_minus_feq(directions, 1, instance_height)
            real(8) :: equilibriums(directions, instance_width, instance_height)

            img = this_image(lattice)
            total_images = num_images()
            img_grid_y = (total_images + coarray_dimensions - 1) / coarray_dimensions

            ! Pressure and macroscopic calculations [cite: 98]
            density_in = poiseuille_in_pressure / (1.0_8 / 3.0_8)
            density_out = poiseuille_out_pressure / (1.0_8 / 3.0_8)
            
            density_arr = calculate_density_array(lattice)
            velocity_arr = calculate_average_velocity_array(lattice, density_arr)
            equilibriums = calculate_equilibrium(density_arr, velocity_arr)

            if (img(1) == coarray_dimensions) then
                f_star_east_minus_feq(:,1,:) = lattice(:, instance_width-1, :) - equilibriums(:, instance_width-1, :)
            end if
            
            if (img(1) == 1) then
                f_star_west_minus_feq(:,1,:) = lattice(:, 2, :) - equilibriums(:, 2, :)
            end if

            ! Calculate target velocities and densities on the boundaries [cite: 100, 102]
            if (img(1) == coarray_dimensions) then
                velocity_arr(instance_width-1,:)%x = velocity_arr(instance_width-1,:)%x * (density_arr(instance_width-1,:) / density_in)
                velocity_arr(instance_width-1,:)%y = velocity_arr(instance_width-1,:)%y * (density_arr(instance_width-1,:) / density_in)
                density_arr(instance_width-1,:) = density_in
            end if
            
            if (img(1) == 1) then
                velocity_arr(2,:)%x = velocity_arr(2,:)%x * (density_arr(2,:) / density_out)
                velocity_arr(2,:)%y = velocity_arr(2,:)%y * (density_arr(2,:) / density_out)
                density_arr(2,:) = density_out
            end if

            ! Re-evaluate equilibrium with updated target densities [cite: 102]
            equilibriums = calculate_equilibrium(density_arr, velocity_arr)

            ! Apply periodic boundary with pressure differential across images 
            if (img(1) == coarray_dimensions) then
                ! Push East boundary to the West-most image's left ghost node
                lattice(:, 1, :)[1, img(2)] = equilibriums(:, instance_width-1, :) + f_star_east_minus_feq(:,1,:)
            end if

            if (img(1) == 1) then
                ! Push West boundary to the East-most image's right ghost node
                lattice(:, instance_width, :)[coarray_dimensions, img(2)] = equilibriums(:, 2, :) + f_star_west_minus_feq(:,1,:)
            end if
            
            sync all

            ! Top bounce-back [cite: 107, 108]
            if (img(2) == 1) then
                lattice(5,:,1) = lattice(3,:,2)
                lattice(8,2:instance_width,1) = lattice(6,1:instance_width-1,2)
                lattice(9,1:instance_width-1,1) = lattice(7,2:instance_width,2)
            end if

            ! Bottom bounce-back [cite: 105, 106]
            if (img(2) == img_grid_y) then
                lattice(3,:,instance_height) = lattice(5,:,instance_height-1)
                lattice(6,1:instance_width-1,instance_height) = lattice(8,2:instance_width,instance_height-1)
                lattice(7,2:instance_width,instance_height) = lattice(9,1:instance_width-1,instance_height-1)
            end if
        end subroutine poiseuille_boundary

        function calculate_density_array(lattice)

            real(8), intent(in) :: lattice(directions, instance_width, instance_height)

            real(8) :: calculate_density_array(instance_width, instance_height)

            calculate_density_array = 1.0_8

            calculate_density_array(2:instance_width-1, 2:instance_height-1) = sum(lattice(:, 2:instance_width-1, 2:instance_height-1), dim=1)

        end function calculate_density_array

        function calculate_average_velocity_array(lattice, density_arr)

            real(8), intent(in) :: lattice(directions, instance_width, instance_height)
            real(8), intent(in) :: density_arr(instance_width, instance_height)

            type(velocity) :: calculate_average_velocity_array(instance_width, instance_height)
            integer :: i

            calculate_average_velocity_array%x = 0.0_8
            calculate_average_velocity_array%y = 0.0_8

            do i=1, directions
                calculate_average_velocity_array(2:instance_width-1, 2:instance_height-1)%x = calculate_average_velocity_array(2:instance_width-1, 2:instance_height-1)%x + lattice(i,2:instance_width-1, 2:instance_height-1) * shift_directions_x(i)
                calculate_average_velocity_array(2:instance_width-1, 2:instance_height-1)%y = calculate_average_velocity_array(2:instance_width-1, 2:instance_height-1)%y + lattice(i,2:instance_width-1, 2:instance_height-1) * shift_directions_y(i)
            end do
            calculate_average_velocity_array(2:instance_width-1, 2:instance_height-1)%x = calculate_average_velocity_array(2:instance_width-1, 2:instance_height-1)%x / density_arr(2:instance_width-1, 2:instance_height-1)
            calculate_average_velocity_array(2:instance_width-1, 2:instance_height-1)%y = calculate_average_velocity_array(2:instance_width-1, 2:instance_height-1)%y / density_arr(2:instance_width-1, 2:instance_height-1)

        end function calculate_average_velocity_array

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

        ! Initialize lattice with (semi)-random values !
        subroutine populate_lattice_random(lattice)

            real(8), intent(inout) :: lattice(directions, instance_width, instance_height)
            integer :: i, j, k
            real(8) :: u

            call random_seed()

            ! populate with random numbers !
            do j=1, instance_width
                do k=1, instance_height
                    do i=1, directions
                        call random_number(u)
                        lattice(i,j,k) = 0.07 + (u * 0.02)
                    end do
                end do
            end do

        end subroutine populate_lattice_random

        ! Initialize lattice with dense center !
        subroutine populate_lattice_dense_center(lattice)

            real(8), intent(inout) :: lattice(directions, instance_width, instance_height)[coarray_dimensions,*]
            integer :: j, k, global_j, global_k
            integer :: img(2)
            real(8) :: epsilon1 = 0.01
            real(8) :: global_nodes

            img = this_image(lattice)

            ! on one point set f(0) to 4/9 + 4/9 * epsilon, f(1,2,3,4) to 1/9 + 1/9 * epsilon, f(5,6,7,8) to 1/36 + 1/36 + epsilon. now roh = 1+epsilon
            ! all other points subtract 4/9 * 1/(m*n-1) * epsilon

            ! populate with dense center !
            global_nodes = real(global_width * global_height, 8)

            do j=1, instance_width
                do k=1, instance_height
                    
                    ! Map local ghost/inner nodes to global coordinates
                    global_j = (img(1) - 1) * (instance_width - 2) + (j - 1)
                    global_k = (img(2) - 1) * (instance_height - 2) + (k - 1)

                    ! Check against the global center
                    if (global_j == global_width/2 .and. global_k == global_height/2) then
                        lattice(1,j,k) = weights(1) + weights(1) * epsilon1
                        lattice(2:5,j,k) = weights(2) + weights(2) * epsilon1
                        lattice(6:,j,k) = weights(6) + weights(6) * epsilon1
                    else
                        lattice(1,j,k) = weights(1) - weights(1) * epsilon1 * (1.0_8 / (global_nodes - 1.0_8))
                        lattice(2:5,j,k) = weights(2) - weights(2) * epsilon1 * (1.0_8 / (global_nodes - 1.0_8))
                        lattice(6:,j,k) = weights(6) - weights(6) * epsilon1  * (1.0_8 / (global_nodes - 1.0_8))
                    end if
                    
                end do
            end do

        end subroutine populate_lattice_dense_center

        subroutine populate_lattice_shear_wave(lattice)

            real(8), intent(inout) :: lattice(directions, instance_width, instance_height)[coarray_dimensions,*]

            real(8) :: density_arr(instance_width, instance_height)
            type(velocity) :: velocity_arr(instance_width, instance_height)

            integer :: y_pos, global_k
            integer :: img(2)

            img = this_image(lattice)

            density_arr = 1.0_8
            velocity_arr%x = 0.0_8
            velocity_arr%y = 0.0_8

            do y_pos=2, instance_height-1
                global_k = (img(2) - 1) * (instance_height - 2) + y_pos
                velocity_arr(:,y_pos)%x = epsilon * sin((2.0_8 * pi * (global_k-1)) / (global_height))
            end do

            lattice = calculate_equilibrium(density_arr, velocity_arr)

        end subroutine populate_lattice_shear_wave

        subroutine populate_lattice_couette(lattice)

            real(8), intent(inout) :: lattice(directions, instance_width, instance_height)

            real(8) :: density_arr(instance_width, instance_height)
            type(velocity) :: velocity_arr(instance_width, instance_height)

            boundary_configuration = 1

            density_arr = 1.0_8
            velocity_arr = velocity(0.0_8, 0.0_8)

            lattice = calculate_equilibrium(density_arr, velocity_arr)

        end subroutine populate_lattice_couette

        subroutine populate_lattice_poiseuille(lattice)

            real(8), intent(inout) :: lattice(directions, instance_width, instance_height)

            real(8) :: density_arr(instance_width, instance_height)
            type(velocity) :: velocity_arr(instance_width, instance_height)

            boundary_configuration = 2

            density_arr = 1.0_8
            velocity_arr = velocity(0.0_8, 0.0_8)

            lattice = calculate_equilibrium(density_arr, velocity_arr)

        end subroutine populate_lattice_poiseuille

        subroutine populate_lattice_sliding_lid(lattice)

            real(8), intent(inout) :: lattice(directions, instance_width, instance_height)

            real(8) :: density_arr(instance_width, instance_height)
            type(velocity) :: velocity_arr(instance_width, instance_height)

            boundary_configuration = 3

            density_arr = 1.0_8
            velocity_arr = velocity(0.0_8, 0.0_8)

            call set_lid_velocity_given_reynolds(1000.0_8)

            lattice = calculate_equilibrium(density_arr, velocity_arr)

        end subroutine populate_lattice_sliding_lid

        function calculate_analytical_viscosity()

            real(8) :: calculate_analytical_viscosity
            
            calculate_analytical_viscosity = (1.0_8 / 3.0_8) * ((1.0_8 / omega) - 0.5_8)

        end function calculate_analytical_viscosity

        subroutine check_viscosity(intial_a, new_a, time_step)

            real(8), intent(in) :: intial_a, new_a
            integer, intent(in) :: time_step

            real(8) :: log_term, calculated_viscosity, k, analytical_viscosity

            k = (2.0_8 * pi) / instance_height
            log_term = log(new_a / intial_a)
            calculated_viscosity = log_term * (-1.0_8 / (k**2 * time_step))

            analytical_viscosity = calculate_analytical_viscosity()

            ! write(interval_str, '(I0)') step_num
            ! file_name = "./Visualization/output-" // trim(adjustl(interval_str)) // "-.txt"
            ! open(1, file=file_name, status="replace", action="write")
            ! do j = 1, global_width
            !     do k = 1, global_height
            !         write(1, *) j, ", ", k, ", ", global_density(j,k), ", ", global_velocity(j,k)%x, ", ", global_velocity(j,k)%y
            !     end do
            ! end do
            ! close(1)

            print *, "Calculated viscosity: ", calculated_viscosity
            print *, "Analytical viscosity: ", analytical_viscosity
            print *, "Difference: ", calculated_viscosity - analytical_viscosity


        end subroutine check_viscosity

        subroutine check_poiseuille_velocity(point_density, point_speed, y_pos)

            real(8), intent(in) :: point_density, point_speed
            integer :: y_pos

            real(8) :: analytical_viscosity, mu, analytical_speed

            analytical_viscosity = calculate_analytical_viscosity()

            mu = point_density * analytical_viscosity

            analytical_speed = (-1.0_8 / (2.0_8 * mu)) * ((poiseuille_out_pressure - poiseuille_in_pressure) / instance_width) * (y_pos - 1.5_8) * ((instance_height - 0.5_8) - y_pos)

            print *, "Calculated speed: ", point_speed
            print *, "Analytical: ", analytical_speed
            print *, "Difference: ", point_speed - analytical_speed

        end subroutine check_poiseuille_velocity

        subroutine set_lid_velocity_given_reynolds(reynolds_number)

            real(8), intent(in) ::  reynolds_number
            
            wall_speed = velocity((reynolds_number * calculate_analytical_viscosity()) / (coarray_dimensions * (instance_width - 2.0_8)), 0.0_8)

        end subroutine set_lid_velocity_given_reynolds
        
end module lattice_boltzmann