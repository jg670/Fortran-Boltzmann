module lattice_boltzmann
    implicit none

    type :: velocity
        real(8) :: x
        real(8) :: y
    end type velocity

    integer, parameter :: height = 30, width = 30, directions = 9
    real(8), parameter :: omega = 1
    real(8), parameter :: shift_directions_x(directions) = [0.0,  1.0,  0.0, -1.0,  0.0,  1.0, -1.0, -1.0,  1.0]
    real(8), parameter :: shift_directions_y(directions) = [0.0,  0.0,  1.0,  0.0, -1.0,  1.0,  1.0, -1.0, -1.0]
    real(8), parameter :: weights(directions) = [4.0_8/9.0_8, 1.0_8/9.0_8, 1.0_8/9.0_8, 1.0_8/9.0_8, 1.0_8/9.0_8, 1.0_8/36.0_8, 1.0_8/36.0_8, 1.0_8/36.0_8, 1.0_8/36.0_8]

    real(8), parameter :: pi = 4.0_8 * atan(1.0_8)
    real(8), parameter :: epsilon = 0.01

    contains

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

            do y_pos=1, height
                velocity_arr(:,y_pos)%x = epsilon * sin((2.0_8 * pi * (y_pos-1)) / height)
            end do

            lattice = calculate_equilibrium(density_arr, velocity_arr)

        end subroutine populate_lattice_shear_wave
        
end module lattice_boltzmann