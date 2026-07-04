module lattice_boltzmann
    implicit none

    type :: velocity
        real(8) :: x
        real(8) :: y
    end type velocity

    integer, parameter :: height = 32, width = 32, directions = 9
    real(8), parameter :: omega = 1
    real(8), parameter :: shift_directions_x(directions) = [0.0,  1.0,  0.0, -1.0,  0.0,  1.0, -1.0, -1.0,  1.0]
    real(8), parameter :: shift_directions_y(directions) = [0.0,  0.0,  1.0,  0.0, -1.0,  1.0,  1.0, -1.0, -1.0]
    real(8), parameter :: weights(directions) = [4.0_8/9.0_8, 1.0_8/9.0_8, 1.0_8/9.0_8, 1.0_8/9.0_8, 1.0_8/9.0_8, 1.0_8/36.0_8, 1.0_8/36.0_8, 1.0_8/36.0_8, 1.0_8/36.0_8]

    real(8), parameter :: pi = 4.0_8 * atan(1.0_8)
    real(8), parameter :: epsilon = 0.01
    
    type(velocity), parameter :: wall_speed = velocity(1.0_8, 0.0_8)

    real(8), parameter :: poiseuille_in_pressure = 0.3, poiseuille_out_pressure = 0.29

    ! encoding for boundary configuration. 0 is normal periodic boundaries, 1 is Couette flow, 2 is Poiseuille flow, and 3 is sliding lid !
    integer :: boundary_configuration = 0

    contains

        function perform_one_time_step(lattice)

            real(8), intent(in) :: lattice(directions, width, height)

            real(8) :: perform_one_time_step(directions, width, height)

            ! HPC with GPUs says collision first and most places seemed to agree so that is how I implemented it. !
            perform_one_time_step = collision_step(lattice)

            ! Perform boundary conditions now to avoid potential errors with bad data leaking into the simulation from ghost nodes when doing the first step !
            if (boundary_configuration == 0) then
                perform_one_time_step = periodic_boundary(perform_one_time_step)

            else if (boundary_configuration == 1) then
                perform_one_time_step = couette_boundary(perform_one_time_step)

            else if (boundary_configuration == 2) then
                perform_one_time_step = poiseuille_boundary(perform_one_time_step)

            else if (boundary_configuration == 3) then
                perform_one_time_step = sliding_lid_boundary(perform_one_time_step)

            end if

            perform_one_time_step = streaming_step(perform_one_time_step)

        end function perform_one_time_step

        ! Perform a streaming step, return the new lattice !
        function streaming_step(lattice)

            real(8), intent(in) :: lattice(directions, width, height)

            real(8) :: streaming_step(directions, width, height)
            integer :: i 

            ! Perform the array shifts !
            do i=1, directions
                streaming_step(i,:,:) = cshift(lattice(i,:,:), int(-shift_directions_x(i)), 1)
                streaming_step(i,:,:) = cshift(streaming_step(i,:,:), int(shift_directions_y(i)), 2)
            end do

        end function streaming_step

        ! Perform a collision step, return the new lattice !
        function collision_step(lattice)

            real(8), intent(in) :: lattice(directions, width, height)

            real(8) :: collision_step(directions, width, height)
            integer :: i
            real(8) :: density_arr(width, height)
            type(velocity) :: velocity_arr(width, height)
            real(8) :: equilibriums(directions, width, height)
            
            density_arr = calculate_density_array(lattice)

            velocity_arr = calculate_average_velocity_array(lattice, density_arr)

            equilibriums = calculate_equilibrium(density_arr, velocity_arr)

            do i=1, directions
                collision_step(i,:,:) = lattice(i,:,:) + omega * (equilibriums(i,:,:) - lattice(i,:,:))
            end do

        end function collision_step

        function periodic_boundary(lattice)

            real(8), intent(in) :: lattice(directions, width, height)

            real(8) :: periodic_boundary(directions, width, height)

            periodic_boundary = lattice

            periodic_boundary(:,1,:) = periodic_boundary(:,width-1,:)
            periodic_boundary(:,width,:) = periodic_boundary(:,2,:)
            periodic_boundary(:,:,1) = periodic_boundary(:,:,height-1)
            periodic_boundary(:,:,height) = periodic_boundary(:,:,2)

        end function periodic_boundary

        function couette_boundary(lattice)

            real(8), intent(in) :: lattice(directions, width, height)

            real(8) :: couette_boundary(directions, width, height)
            real(8) :: average_density

            couette_boundary = lattice

            ! Side walls with periodic boundary conditions !
            couette_boundary(:,1,:) = couette_boundary(:, width-1,:)
            couette_boundary(:,width,:) = couette_boundary(:,2,:)

            ! Bottom bounce-back boundary (indexes are +1 since fortran is 1-indexed) !
            couette_boundary(3,:,height) = couette_boundary(5,:,height-1)
            couette_boundary(6,1:width-1,height) = couette_boundary(8,2:width,height-1)
            couette_boundary(7,2:width,height) = couette_boundary(9,1:width-1,height-1)

            average_density = sum(lattice) / (width * height)

            ! Top moving boundary (y velocity is 0, so nothing is added) !
            couette_boundary(5,:,1) = couette_boundary(3,:,2)
            couette_boundary(8,2:width,1) = couette_boundary(6,1:width-1,2) - 2.0_8 * weights(6) * average_density * ((shift_directions_x(6) * wall_speed%x) / (1.0_8 / 3.0_8))
            couette_boundary(9,1:width-1,1) = couette_boundary(7,2:width,2) - 2.0_8 * weights(7) * average_density * ((shift_directions_x(7) * wall_speed%x) / (1.0_8 / 3.0_8))

        end function couette_boundary

        function poiseuille_boundary(lattice)

            real(8), intent(in) :: lattice(directions, width, height)

            real(8) :: poiseuille_boundary(directions, width, height)

            real(8) :: density_in, density_out
            real(8) :: density_arr(width, height)
            type(velocity) :: velocity_arr(width, height)
            real(8) :: f_star_west_minus_feq(directions, 1, height), f_star_east_minus_feq(directions, 1, height)
            real(8) :: equilibriums(directions, width, height)

            poiseuille_boundary = lattice

            ! Pressure calculations !
            density_in = poiseuille_in_pressure / (1.0_8 / 3.0_8)
            density_out = poiseuille_out_pressure / (1.0_8 / 3.0_8)

            density_arr = calculate_density_array(lattice)
            velocity_arr = calculate_average_velocity_array(lattice, density_arr)

            equilibriums = calculate_equilibrium(density_arr, velocity_arr)

            f_star_east_minus_feq(:,1,:) = lattice(:,width-1,:) - equilibriums(:,width-1,:)
            f_star_west_minus_feq(:,1,:) = lattice(:,2,:) - equilibriums(:,2,:)

            ! Update velocities to match new pressure !
            velocity_arr(width-1,:)%x = velocity_arr(width-1,:)%x * (density_arr(width-1,:) / density_in)
            velocity_arr(width-1,:)%y = velocity_arr(width-1,:)%y * (density_arr(width-1,:) / density_in)
        
            velocity_arr(2,:)%x = velocity_arr(2,:)%x * (density_arr(2,:) / density_out)
            velocity_arr(2,:)%y = velocity_arr(2,:)%y * (density_arr(2,:) / density_out)

            ! Set density array east border to density_in and west border to density_out !
            density_arr(width-1,:) = density_in
            density_arr(2,:) = density_out

            equilibriums = calculate_equilibrium(density_arr, velocity_arr)

            ! Side walls with periodic boundary conditions and pressure differential !
            poiseuille_boundary(:,1,:) = equilibriums(:,width-1,:) + f_star_east_minus_feq(:,1,:)
            poiseuille_boundary(:,width,:) = equilibriums(:,2,:) + f_star_west_minus_feq(:,1,:)


            ! Bottom bounce-back boundary (indexes are +1 since fortran is 1-indexed) !
            poiseuille_boundary(3,:,height) = poiseuille_boundary(5,:,height-1)
            poiseuille_boundary(6,1:width-1,height) = poiseuille_boundary(8,2:width,height-1)
            poiseuille_boundary(7,2:width,height) = poiseuille_boundary(9,1:width-1,height-1)

            ! Top bounce-back boundary (indexes are +1 since fortran is 1-indexed) !
            poiseuille_boundary(5,:,1) = poiseuille_boundary(3,:,2)
            poiseuille_boundary(8,2:width,1) = poiseuille_boundary(6,1:width-1,2)
            poiseuille_boundary(9,1:width-1,1) = poiseuille_boundary(7,2:width,2)

        end function poiseuille_boundary

        function sliding_lid_boundary(lattice)

            real(8), intent(in) ::  lattice(directions, width, height)

            real(8) :: sliding_lid_boundary(directions, width, height)

            real(8) :: average_density

            sliding_lid_boundary = lattice

            average_density = sum(lattice) / (width * height)

            ! Top moving boundary (y velocity is 0, so nothing is added) !
            sliding_lid_boundary(5,:,1) = sliding_lid_boundary(3,:,2)
            sliding_lid_boundary(8,2:width,1) = sliding_lid_boundary(6,1:width-1,2) - 2.0_8 * weights(6) * average_density * ((shift_directions_x(6) * wall_speed%x) / (1.0_8 / 3.0_8))
            sliding_lid_boundary(9,1:width-1,1) = sliding_lid_boundary(7,2:width,2) - 2.0_8 * weights(7) * average_density * ((shift_directions_x(7) * wall_speed%x) / (1.0_8 / 3.0_8))

            ! Bottom bounce-back boundary (indexes are +1 since fortran is 1-indexed) !
            sliding_lid_boundary(3,:,height) = sliding_lid_boundary(5,:,height-1)
            sliding_lid_boundary(6,1:width-1,height) = sliding_lid_boundary(8,2:width,height-1)
            sliding_lid_boundary(7,2:width,height) = sliding_lid_boundary(9,1:width-1,height-1)

            ! Left side bounce-back boundary (indexes are +1 since fortran is 1-indexed) !
            sliding_lid_boundary(2,1,:) = sliding_lid_boundary(4,2,:)
            sliding_lid_boundary(6,1,2:height) = sliding_lid_boundary(8,2,1:height-1)
            sliding_lid_boundary(9,1,1:height-1) = sliding_lid_boundary(7,2,2:height)

            ! Right side bounce-back boundary (indexes are +1 since fortran is 1-indexed) !
            sliding_lid_boundary(4,width,:) = sliding_lid_boundary(2, width-1,:)
            sliding_lid_boundary(7,width,2:height) = sliding_lid_boundary(9,width-1,1:height-1)
            sliding_lid_boundary(8,width,1:height-1) = sliding_lid_boundary(6,width-1,2:height)

        end function sliding_lid_boundary

        function calculate_density_array(lattice)

            real(8), intent(in) :: lattice(directions, width, height)

            real(8) :: calculate_density_array(width, height)

            calculate_density_array = sum(lattice, dim=1)

        end function calculate_density_array

        function calculate_average_velocity_array(lattice, density_arr)

            real(8), intent(in) :: lattice(directions, width, height)
            real(8), intent(in) :: density_arr(width, height)

            type(velocity) :: calculate_average_velocity_array(width, height)
            integer :: i

            calculate_average_velocity_array%x = 0.0_8
            calculate_average_velocity_array%y = 0.0_8

            do i=1, directions
                calculate_average_velocity_array%x = calculate_average_velocity_array%x + lattice(i,:,:) * shift_directions_x(i)
                calculate_average_velocity_array%y = calculate_average_velocity_array%y + lattice(i,:,:) * shift_directions_y(i)
            end do
            calculate_average_velocity_array%x = calculate_average_velocity_array%x / density_arr
            calculate_average_velocity_array%y = calculate_average_velocity_array%y / density_arr

        end function calculate_average_velocity_array

        function calculate_equilibrium(density_arr, velocity_arr)

            real(8), intent(in) :: density_arr(width, height)
            type(velocity), intent(in) :: velocity_arr(width, height)

            real(8) :: calculate_equilibrium(directions, width, height)

            integer :: i
            real(8) :: u_squared(width, height), velocity_sum(width, height)

            u_squared = velocity_arr%x**2 + velocity_arr%y**2

            do i=1, directions
                velocity_sum = shift_directions_x(i) * velocity_arr%x + shift_directions_y(i) * velocity_arr%y

                calculate_equilibrium(i,:,:) = weights(i) * density_arr * (1.0_8 + 3.0_8 * velocity_sum + (4.5_8 * velocity_sum**2) - 1.5_8 * u_squared)
            end do

        end function calculate_equilibrium

        ! Initialize lattice with (semi)-random values !
        subroutine populate_lattice_random(lattice)

            real(8), intent(inout) :: lattice(directions, width, height)
            integer :: i, j, k
            real(8) :: u

            call random_seed()

            ! populate with random numbers !
            do j=1, width
                do k=1, height
                    do i=1, directions
                        call random_number(u)
                        lattice(i,j,k) = 0.07 + (u * 0.02)
                    end do
                end do
            end do

        end subroutine populate_lattice_random

        ! Initialize lattice with dense center !
        subroutine populate_lattice_dense_center(lattice)

            real(8), intent(inout) :: lattice(directions, width, height)
            integer :: j, k
            real(8) :: epsilon1 = 0.01

            ! on one point set f(0) to 4/9 + 4/9 * epsilon, f(1,2,3,4) to 1/9 + 1/9 * epsilon, f(5,6,7,8) to 1/36 + 1/36 + epsilon. now roh = 1+epsilon
            ! all other points subtract 4/9 * 1/(m*n-1) * epsilon

            ! populate with dense center !
            do j=1, width
                do k=1, height
                    if (j == width/2 .and. k == height/2) then
                        lattice(1,j,k) = weights(1) + weights(1) * epsilon1
                        lattice(2:5,j,k) = weights(2) + weights(2) * epsilon1
                        lattice(6:,j,k) = weights(6) + weights(6) * epsilon1
                    else
                        lattice(1,j,k) = weights(1) - weights(1) * epsilon1 * (1.0 / ((width * height) - 1.0))
                        lattice(2:5,j,k) = weights(2) - weights(2) * epsilon1 * (1.0 / ((width * height) - 1.0))
                        lattice(6:,j,k) = weights(6) - weights(6) * epsilon1  * (1.0 / ((width * height) - 1.0))
                    end if
                end do
            end do

        end subroutine populate_lattice_dense_center

        subroutine populate_lattice_shear_wave(lattice)

            real(8), intent(inout) :: lattice(directions, width, height)

            real(8) :: density_arr(width, height)
            type(velocity) :: velocity_arr(width, height)

            integer :: y_pos

            density_arr = 1.0_8
            velocity_arr%y = 0.0_8

            do y_pos=2, height-1
                velocity_arr(:,y_pos)%x = epsilon * sin((2.0_8 * pi * (y_pos-1)) / (height-2))
            end do

            lattice = calculate_equilibrium(density_arr, velocity_arr)

        end subroutine populate_lattice_shear_wave

        subroutine populate_lattice_couette(lattice)

            real(8), intent(inout) :: lattice(directions, width, height)

            real(8) :: density_arr(width, height)
            type(velocity) :: velocity_arr(width, height)

            boundary_configuration = 1

            density_arr = 1.0_8
            velocity_arr = velocity(0.0_8, 0.0_8)

            lattice = calculate_equilibrium(density_arr, velocity_arr)

        end subroutine populate_lattice_couette

        subroutine populate_lattice_poiseuille(lattice)

            real(8), intent(inout) :: lattice(directions, width, height)

            real(8) :: density_arr(width, height)
            type(velocity) :: velocity_arr(width, height)

            boundary_configuration = 2

            density_arr = 1.0_8
            velocity_arr = velocity(0.0_8, 0.0_8)

            lattice = calculate_equilibrium(density_arr, velocity_arr)

        end subroutine populate_lattice_poiseuille

        subroutine populate_lattice_sliding_lid(lattice)

            real(8), intent(inout) :: lattice(directions, width, height)

            real(8) :: density_arr(width, height)
            type(velocity) :: velocity_arr(width, height)

            boundary_configuration = 3

            density_arr = 1.0_8
            velocity_arr = velocity(0.0_8, 0.0_8)

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

            k = (2.0_8 * pi) / height
            log_term = log(new_a / intial_a)
            calculated_viscosity = log_term * (-1.0_8 / (k**2 * time_step))

            analytical_viscosity = calculate_analytical_viscosity()

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

            analytical_speed = (-1.0_8 / (2.0_8 * mu)) * ((poiseuille_out_pressure - poiseuille_in_pressure) / width) * (y_pos - 1.5_8) * ((height - 0.5_8) - y_pos)

            print *, "Calculated speed: ", point_speed
            print *, "Analytical: ", analytical_speed
            print *, "Difference: ", point_speed - analytical_speed

        end subroutine check_poiseuille_velocity
        
end module lattice_boltzmann