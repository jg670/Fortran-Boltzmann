module lattice_output_operations
    use lattice_boltzmann

    implicit none

    contains

        ! More or less just for visualization !
        function calculate_density_array(lattice)

            real(8), intent(in) :: lattice(directions, width, height)

            real(8) :: calculate_density_array(width, height)
            integer :: j, k

            ! Calculate density by summing over all directions !
            do j=1, width
                do k=1, height
                    calculate_density_array(j,k) = sum(lattice(:,j,k))
                end do
            end do

        end function calculate_density_array

        ! More or less just for visualization !
        function calculate_average_velocity_array(lattice, density_arr)

            real(8), intent(in) :: lattice(directions, width, height)
            real(8), intent(in) :: density_arr(width, height)

            type(velocity) :: calculate_average_velocity_array(width, height)
            integer :: i, j, k

            do j=1, width
                do k=1, height
                    ! Average x velocity !
                    calculate_average_velocity_array(j,k)%x = sum(lattice(:,j,k) * shift_directions_x(:)) / density_arr(j,k)
                    ! Average y velocity !
                    calculate_average_velocity_array(j,k)%y = sum(lattice(:,j,k) * shift_directions_y(:)) / density_arr(j,k)
                end do
            end do

        end function calculate_average_velocity_array

end module lattice_output_operations

program SheBoltOnMyManTilILattice
    use lattice_output_operations
    implicit none

    ! Image data !
    integer :: img, img_dim, img_coords(2)

    ! Loop indexes !
    integer :: i, j, k

    ! Lattice arrays !
    real(8), allocatable :: lattice(:,:,:)[:,:], lattice_new(:,:,:)[:,:], density(:,:)[:,:], density_new(:,:)[:,:]
    type(velocity), allocatable :: average_velocity(:,:)[:,:], average_velocity_new(:,:)[:,:]

    ! Format for printing the lattice !
    !character(len=50) :: format = '(15(f0.2," "),"  ",15(f0.2," "))'
    character(len=50) :: format = '(15(f0.4," "))'

    img_dim = 1

    allocate(lattice(directions, width, height)[img_dim,*])
    allocate(lattice_new(directions, width, height)[img_dim,*])
    allocate(density(width, height)[img_dim,*])
    allocate(density_new(width, height)[img_dim,*])
    allocate(average_velocity(width, height)[img_dim,*])
    allocate(average_velocity_new(width, height)[img_dim,*])

    img = this_image()
    ! img_coords = this_image(lattice)
    ! print *, "Image ", img, " coordinates: (", img_coords, ")"

    ! Populate lattice !
    !call populate_lattice_random(lattice, width, height, directions)
    call populate_lattice_dense_center(lattice)

    ! Perform one streaming step !
    lattice_new = streaming_step(lattice)

    ! Perform one collision step !
    lattice_new = collision_step(lattice_new)

    ! Calculate density arrays !
    density = calculate_density_array(lattice)

    ! Calculate average velocity arrays !
    ! average_velocity = calculate_average_velocity_array(lattice, density, width, height, directions, shift_directions_x, shift_directions_y)
    ! average_velocity_new = calculate_average_velocity_array(lattice_new, density_new, width, height, directions, shift_directions_x, shift_directions_y)

    ! Print initial lattice !
    if (img .eq. 1) then
        print *, " "
        print *, "Initial lattice:"
        do j=1, height
            write(unit=6, fmt=format) density(:,j)
        end do
        print *, " "
        print *, "Initial density: ", sum(lattice)
    end if
    sync all

    do i=1, 1
        do j=1, 100
            ! Perform one streaming step !
            lattice_new = streaming_step(lattice_new)

            ! Perform one collision step !
            lattice_new = collision_step(lattice_new)
        end do

        ! Calculate density array !
        density_new = calculate_density_array(lattice_new)

        ! Print density lattices after each step !
        if (img .eq. 1) then
            print *, " "
            print *, "After ", i*j, " steps:"
            do j=1, height
                write(unit=6, fmt=format) density_new(:,j)
            end do
            print *, " "
            print *, "Density diff: ", sum(lattice_new) - sum(lattice)
        end if
        sync all

        do j=1, 100000
            ! Perform one streaming step !
            lattice_new = streaming_step(lattice_new)

            ! Perform one collision step !
            lattice_new = collision_step(lattice_new)
        end do

        ! Calculate density array !
        density_new = calculate_density_array(lattice_new)

        ! Print density lattices after each step !
        if (img .eq. 1) then
            print *, " "
            print *, "After ", i*j, " steps:"
            do j=1, height
                write(unit=6, fmt=format) density_new(:,j)
            end do
            print *, " "
            print *, "Density diff: ", sum(lattice_new) - sum(lattice)
        end if
        sync all
    end do

end program SheBoltOnMyManTilILattice