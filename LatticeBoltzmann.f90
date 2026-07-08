module lattice_boltzmann
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

    ! Encoding for boundary configuration. 0 is normal periodic boundaries, 1 is Couette flow, 2 is Poiseuille flow, 3 is sliding lid , 4 is for parallel sliding lid !
    integer :: boundary_configuration = 0

    contains

        subroutine perform_one_time_step(lattice)

            real(8), intent(inout) :: lattice(directions, instance_width, instance_height)[coarray_dimensions,*]

            ! HPC with GPUs says collision first and most places seemed to agree so that is how I implemented it. !
            call collision_step(lattice)

            ! Perform boundary conditions now to avoid potential errors with bad data leaking into the simulation from ghost nodes when doing the first step !
            if (boundary_configuration == 0) then
                call periodic_boundary(lattice)

            else if (boundary_configuration == 1) then
                call couette_boundary(lattice)

            else if (boundary_configuration == 2) then
                call poiseuille_boundary(lattice)

            else if (boundary_configuration == 3) then
                call sliding_lid_boundary(lattice)

            else if (boundary_configuration == 4) then
                call sliding_lid_boundary_parallel(lattice)

            end if

            call streaming_step(lattice)

        end subroutine perform_one_time_step

        ! Perform a streaming step, return the new lattice !
        subroutine streaming_step(lattice)

            real(8), intent(inout) :: lattice(directions, instance_width, instance_height)

            integer :: i 

            ! Perform the array shifts !
            do i=1, directions
                lattice(i,:,:) = cshift(lattice(i,:,:), int(-shift_directions_x(i)), 1)
                lattice(i,:,:) = cshift(lattice(i,:,:), int(shift_directions_y(i)), 2)
            end do

        end subroutine streaming_step

        ! Perform a collision step, return the new lattice !
        subroutine collision_step(lattice)

            real(8), intent(inout) :: lattice(directions, instance_width, instance_height)

            integer :: i
            real(8) :: density_arr(instance_width, instance_height)
            type(velocity) :: velocity_arr(instance_width, instance_height)
            real(8) :: equilibriums(directions, instance_width, instance_height)
            
            density_arr = calculate_density_array(lattice)

            velocity_arr = calculate_average_velocity_array(lattice, density_arr)

            equilibriums = calculate_equilibrium(density_arr, velocity_arr)

            do i=1, directions
                lattice(i,2:instance_width-1, 2:instance_height-1) = lattice(i,2:instance_width-1, 2:instance_height-1) + omega * (equilibriums(i,2:instance_width-1, 2:instance_height-1) - lattice(i,2:instance_width-1, 2:instance_height-1))
            end do

        end subroutine collision_step

        subroutine periodic_boundary(lattice)

            real(8), intent(inout) :: lattice(directions, instance_width, instance_height)

            lattice(:,1,:) = lattice(:,instance_width-1,:)
            lattice(:,instance_width,:) = lattice(:,2,:)
            lattice(:,:,1) = lattice(:,:,instance_height-1)
            lattice(:,:,instance_height) = lattice(:,:,2)

        end subroutine periodic_boundary

        subroutine couette_boundary(lattice)

            real(8), intent(inout) :: lattice(directions, instance_width, instance_height)

            real(8) :: average_density

            average_density = sum(lattice) / (instance_width * instance_height)

            ! Side walls with periodic boundary conditions !
            lattice(:,1,:) = lattice(:, instance_width-1,:)
            lattice(:,instance_width,:) = lattice(:,2,:)

            ! Bottom bounce-back boundary (indexes are +1 since fortran is 1-indexed) !
            lattice(3,:,instance_height) = lattice(5,:,instance_height-1)
            lattice(6,1:instance_width-1,instance_height) = lattice(8,2:instance_width,instance_height-1)
            lattice(7,2:instance_width,instance_height) = lattice(9,1:instance_width-1,instance_height-1)

            ! Top moving boundary (y velocity is 0, so nothing is added) !
            lattice(5,:,1) = lattice(3,:,2)
            lattice(8,2:instance_width,1) = lattice(6,1:instance_width-1,2) - 2.0_8 * weights(6) * average_density * ((shift_directions_x(6) * wall_speed%x) / (1.0_8 / 3.0_8))
            lattice(9,1:instance_width-1,1) = lattice(7,2:instance_width,2) - 2.0_8 * weights(7) * average_density * ((shift_directions_x(7) * wall_speed%x) / (1.0_8 / 3.0_8))

        end subroutine couette_boundary

        subroutine poiseuille_boundary(lattice)

            real(8), intent(inout) :: lattice(directions, instance_width, instance_height)

            real(8) :: density_in, density_out
            real(8) :: density_arr(instance_width, instance_height)
            type(velocity) :: velocity_arr(instance_width, instance_height)
            real(8) :: f_star_west_minus_feq(directions, 1, instance_height), f_star_east_minus_feq(directions, 1, instance_height)
            real(8) :: equilibriums(directions, instance_width, instance_height)

            ! Pressure calculations !
            density_in = poiseuille_in_pressure / (1.0_8 / 3.0_8)
            density_out = poiseuille_out_pressure / (1.0_8 / 3.0_8)

            density_arr = calculate_density_array(lattice)
            velocity_arr = calculate_average_velocity_array(lattice, density_arr)

            equilibriums = calculate_equilibrium(density_arr, velocity_arr)

            f_star_east_minus_feq(:,1,:) = lattice(:,instance_width-1,:) - equilibriums(:,instance_width-1,:)
            f_star_west_minus_feq(:,1,:) = lattice(:,2,:) - equilibriums(:,2,:)

            ! Update velocities to match new pressure !
            velocity_arr(instance_width-1,:)%x = velocity_arr(instance_width-1,:)%x * (density_arr(instance_width-1,:) / density_in)
            velocity_arr(instance_width-1,:)%y = velocity_arr(instance_width-1,:)%y * (density_arr(instance_width-1,:) / density_in)
        
            velocity_arr(2,:)%x = velocity_arr(2,:)%x * (density_arr(2,:) / density_out)
            velocity_arr(2,:)%y = velocity_arr(2,:)%y * (density_arr(2,:) / density_out)

            ! Set density array east border to density_in and west border to density_out !
            density_arr(instance_width-1,:) = density_in
            density_arr(2,:) = density_out

            equilibriums = calculate_equilibrium(density_arr, velocity_arr)

            ! Side walls with periodic boundary conditions and pressure differential !
            lattice(:,1,:) = equilibriums(:,instance_width-1,:) + f_star_east_minus_feq(:,1,:)
            lattice(:,instance_width,:) = equilibriums(:,2,:) + f_star_west_minus_feq(:,1,:)


            ! Bottom bounce-back boundary (indexes are +1 since fortran is 1-indexed) !
            lattice(3,:,instance_height) = lattice(5,:,instance_height-1)
            lattice(6,1:instance_width-1,instance_height) = lattice(8,2:instance_width,instance_height-1)
            lattice(7,2:instance_width,instance_height) = lattice(9,1:instance_width-1,instance_height-1)

            ! Top bounce-back boundary (indexes are +1 since fortran is 1-indexed) !
            lattice(5,:,1) = lattice(3,:,2)
            lattice(8,2:instance_width,1) = lattice(6,1:instance_width-1,2)
            lattice(9,1:instance_width-1,1) = lattice(7,2:instance_width,2)

        end subroutine poiseuille_boundary

        subroutine sliding_lid_boundary(lattice)

            real(8), intent(inout) ::  lattice(directions, instance_width, instance_height)

            ! Bottom bounce-back boundary (indexes are +1 since fortran is 1-indexed) !
            lattice(3,:,instance_height) = lattice(5,:,instance_height-1)
            lattice(6,1:instance_width-1,instance_height) = lattice(8,2:instance_width,instance_height-1)
            lattice(7,2:instance_width,instance_height) = lattice(9,1:instance_width-1,instance_height-1)

            ! Left side bounce-back boundary (indexes are +1 since fortran is 1-indexed) !
            lattice(2,1,:) = lattice(4,2,:)
            lattice(6,1,2:instance_height) = lattice(8,2,1:instance_height-1)
            lattice(9,1,1:instance_height-1) = lattice(7,2,2:instance_height)

            ! Right side bounce-back boundary (indexes are +1 since fortran is 1-indexed) !
            lattice(4,instance_width,:) = lattice(2, instance_width-1,:)
            lattice(7,instance_width,2:instance_height) = lattice(9,instance_width-1,1:instance_height-1)
            lattice(8,instance_width,1:instance_height-1) = lattice(6,instance_width-1,2:instance_height)

            ! Top moving boundary (y velocity is 0, so nothing is added) !
            lattice(5,:,1) = lattice(3,:,2)
            lattice(8,2:instance_width,1) = lattice(6,1:instance_width-1,2) - 2.0_8 * weights(6) * sum(lattice(:, 1:instance_width-1, 2), dim=1) * ((shift_directions_x(6) * wall_speed%x) / (1.0_8 / 3.0_8))
            lattice(9,1:instance_width-1,1) = lattice(7,2:instance_width,2) - 2.0_8 * weights(7) * sum(lattice(:, 2:instance_width, 2), dim=1) * ((shift_directions_x(7) * wall_speed%x) / (1.0_8 / 3.0_8))

        end subroutine sliding_lid_boundary

        subroutine sliding_lid_boundary_parallel(lattice)

            real(8), intent(inout) :: lattice(directions, instance_width, instance_height)[coarray_dimensions,*]

            integer :: img(2), total_images

            img = this_image(lattice)
            total_images = num_images()

            ! Handle left border !
            if (img(1) /= 1) then
                lattice(:,instance_width,:)[img(1) - 1, img(2)] = lattice(:,2,:)

            else 
                ! Left side bounce-back boundary (indexes are +1 since fortran is 1-indexed) !
                lattice(2,1,:) = lattice(4,2,:)
                lattice(6,1,2:instance_height) = lattice(8,2,1:instance_height-1)
                lattice(9,1,1:instance_height-1) = lattice(7,2,2:instance_height)

            end if

            ! Handle right border !
            if (img(1) /= coarray_dimensions) then
                lattice(:,1,:)[img(1) + 1, img(2)] = lattice(:,instance_width-1,:)

            else 
                ! Right side bounce-back boundary (indexes are +1 since fortran is 1-indexed) !
                lattice(4,instance_width,:) = lattice(2, instance_width-1,:)
                lattice(7,instance_width,2:instance_height) = lattice(9,instance_width-1,1:instance_height-1)
                lattice(8,instance_width,1:instance_height-1) = lattice(6,instance_width-1,2:instance_height)
            
            end if

            sync all

            ! Handle bottom border !
            if (img(2) /= (total_images / coarray_dimensions)) then
                lattice(:,:,1)[img(1), img(2) + 1] = lattice(:,:,instance_height-1)

            else 
                ! Bottom bounce-back boundary (indexes are +1 since fortran is 1-indexed) !
                lattice(3,:,instance_height) = lattice(5,:,instance_height-1)
                lattice(6,1:instance_width-1,instance_height) = lattice(8,2:instance_width,instance_height-1)
                lattice(7,2:instance_width,instance_height) = lattice(9,1:instance_width-1,instance_height-1)

            end if

            ! Handle top border !
            if (img(2) /= 1) then
                lattice(:,:,instance_height)[img(1), img(2) - 1] = lattice(:,:,2)

            else 
                ! Top moving boundary (y velocity is 0, so nothing is added) !
                lattice(5,:,1) = lattice(3,:,2)
                lattice(8,2:instance_width,1) = lattice(6,1:instance_width-1,2) - 2.0_8 * weights(6) * sum(lattice(:, 1:instance_width-1, 2), dim=1) * ((shift_directions_x(6) * wall_speed%x) / (1.0_8 / 3.0_8))
                lattice(9,1:instance_width-1,1) = lattice(7,2:instance_width,2) - 2.0_8 * weights(7) * sum(lattice(:, 2:instance_width, 2), dim=1) * ((shift_directions_x(7) * wall_speed%x) / (1.0_8 / 3.0_8))

            end if

            sync all

        end subroutine sliding_lid_boundary_parallel

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

            real(8), intent(inout) :: lattice(directions, instance_width, instance_height)
            integer :: j, k
            real(8) :: epsilon1 = 0.01

            ! on one point set f(0) to 4/9 + 4/9 * epsilon, f(1,2,3,4) to 1/9 + 1/9 * epsilon, f(5,6,7,8) to 1/36 + 1/36 + epsilon. now roh = 1+epsilon
            ! all other points subtract 4/9 * 1/(m*n-1) * epsilon

            ! populate with dense center !
            do j=1, instance_width
                do k=1, instance_height
                    if (j == instance_width/2 .and. k == instance_height/2) then
                        lattice(1,j,k) = weights(1) + weights(1) * epsilon1
                        lattice(2:5,j,k) = weights(2) + weights(2) * epsilon1
                        lattice(6:,j,k) = weights(6) + weights(6) * epsilon1
                    else
                        lattice(1,j,k) = weights(1) - weights(1) * epsilon1 * (1.0_8 / ((instance_width * instance_height) - 1.0_8))
                        lattice(2:5,j,k) = weights(2) - weights(2) * epsilon1 * (1.0_8 / ((instance_width * instance_height) - 1.0_8))
                        lattice(6:,j,k) = weights(6) - weights(6) * epsilon1  * (1.0_8 / ((instance_width * instance_height) - 1.0_8))
                    end if
                end do
            end do

        end subroutine populate_lattice_dense_center

        subroutine populate_lattice_shear_wave(lattice)

            real(8), intent(inout) :: lattice(directions, instance_width, instance_height)

            real(8) :: density_arr(instance_width, instance_height)
            type(velocity) :: velocity_arr(instance_width, instance_height)

            integer :: y_pos

            density_arr = 1.0_8
            velocity_arr%y = 0.0_8

            do y_pos=2, instance_height-1
                velocity_arr(:,y_pos)%x = epsilon * sin((2.0_8 * pi * (y_pos-1)) / (instance_height-2))
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

        subroutine populate_lattice_sliding_lid_parallel(lattice)

            real(8), intent(inout) :: lattice(directions, instance_width, instance_height)

            real(8) :: density_arr(instance_width, instance_height)
            type(velocity) :: velocity_arr(instance_width, instance_height)
            integer :: img

            boundary_configuration = 4

            img = this_image()

            density_arr = 1.0_8
            velocity_arr = velocity(0.0_8, 0.0_8)

            call set_lid_velocity_given_reynolds(1000.0_8)

            lattice = calculate_equilibrium(density_arr, velocity_arr)

        end subroutine populate_lattice_sliding_lid_parallel

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