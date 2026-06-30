module lattice_boltzmann
  implicit none

  type :: velocity
    real(8) :: x
    real(8) :: y
  end type velocity

  integer, parameter :: height = 10, width = 15, directions = 9
  real(8), parameter :: omega = 1.0
  real(8), parameter :: shift_directions_x(directions) = [0.0,  1.0,  0.0, -1.0,  0.0,  1.0, -1.0, -1.0,  1.0]
  real(8), parameter :: shift_directions_y(directions) = [0.0,  0.0,  1.0,  0.0, -1.0,  1.0,  1.0, -1.0, -1.0]
  real(8), parameter :: weights(directions) = [4.0_8/9.0_8, 1.0_8/9.0_8, 1.0_8/9.0_8, 1.0_8/9.0_8, 1.0_8/9.0_8, 1.0_8/36.0_8, 1.0_8/36.0_8, 1.0_8/36.0_8, 1.0_8/36.0_8]

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
      integer :: i, j, k

      real(8) :: average_density, equilibrium, velocity_sum
      type(velocity) :: average_velocity

      do j=1, width
        do k=1, height
          average_density = sum(lattice(:,j,k))
          average_velocity%x = sum(lattice(:,j,k) * shift_directions_x(:)) / average_density
          average_velocity%y = sum(lattice(:,j,k) * shift_directions_y(:)) / average_density
          do i=1, directions
            velocity_sum = shift_directions_x(i) * average_velocity%x + shift_directions_y(i) * average_velocity%y
            equilibrium = weights(i) * average_density * (1.0 + 3.0 * velocity_sum + (9.0 * velocity_sum**2) / 2.0 - 1.5 * (average_velocity%x**2 + average_velocity%y**2))
            collision_step(i,j,k) = lattice(i,j,k) + omega * (equilibrium - lattice(i,j,k))
          end do
        end do
      end do            

    end function collision_step

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
      integer :: i, j, k
      real(8) :: epsilon = 0.01

      ! on one point set f(0) to 4/9 + 4/9 * epsilon, f(1,2,3,4) to 1/9 + 1/9 * epsilon, f(5,6,7,8) to 1/36 + 1/36 + epsilon. now roh = 1+epsilon
      ! all other points subtract 4/9 * 1/(m*n-1) * epsilon

      ! populate with dense center !
      do j=1, width
        do k=1, height
          if (j == width/2 .and. k == height/2) then
            lattice(1,j,k) = weights(1) + weights(1) * epsilon
            lattice(2:5,j,k) = weights(2) + weights(2) * epsilon
            lattice(6:,j,k) = weights(6) + weights(6) * epsilon
          else
            lattice(1,j,k) = weights(1) - weights(1) * epsilon * (1.0 / ((width * height) - 1.0))
            lattice(2:5,j,k) = weights(2) - weights(2) * epsilon * (1.0 / ((width * height) - 1.0))
            lattice(6:,j,k) = weights(6) - weights(6) * epsilon  * (1.0 / ((width * height) - 1.0))
          end if
        end do
      end do

    end subroutine populate_lattice_dense_center

end module lattice_boltzmann