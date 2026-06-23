module lattice_operations
    contains
        subroutine populate_lattice(lattice, width, height, directions)
            implicit none

            integer, intent(in) :: width, height, directions
            real, intent(inout) :: lattice(directions, width, height)

            integer :: i, j, k
            real :: u

            ! populate with random numbers !
            do i=1, directions
                do j=1, width
                    do k=1, height
                        call random_number(u)
                        lattice(i,j,k) = u
                    end do
                end do
            end do

        end subroutine populate_lattice

        function streaming_step(lattice, width, height, directions)
            implicit none

            integer, intent(in) :: width, height, directions
            real, intent(in) :: lattice(directions, width, height)

            real :: streaming_step(directions, width, height)
            integer :: i, j, k
            integer :: shift_directions_x(directions), shift_directions_y(directions)

            ! Initialize shift directions for D2Q9 model !
            shift_directions_x = [0,  1,  0, -1,  0,  1, -1, -1,  1]
            shift_directions_y = [0,  0,  1,  0, -1,  1,  1, -1, -1] 

            ! Perform the array shifts !
            do i=1, directions
                streaming_step(i,:,:) = cshift(lattice(i,:,:), -shift_directions_x(i), 1)
                streaming_step(i,:,:) = cshift(streaming_step(i,:,:), shift_directions_y(i), 2)
            end do
        end function streaming_step

        function calculate_density_array(lattice, width, height, directions)
            implicit none

            integer, intent(in) :: width, height, directions
            real, intent(in) :: lattice(directions, width, height)

            real :: calculate_density_array(width, height)
            integer :: i, j, k

            ! Calculate density by summing over all directions !
            do j=1, width
                do k=1, height
                    calculate_density_array(j,k) = sum(lattice(:,j,k))
                end do
            end do

        end function calculate_density_array

        ! function handle_ghost_nodes(lattice, width, height, directions)
        !     implicit none

        !     integer, intent(in) :: width, height, directions
        !     real, intent(in) :: lattice(directions, width, height)
        !     real :: handle_ghost_nodes(directions, width, height)

        !     integer :: i

            
            
        ! end function handle_ghost_nodes

end module lattice_operations

program SheBoltOnMyManTilILattice
    use lattice_operations
    implicit none

    ! Dimensions of the lattice !
    integer :: height, width, directions

    ! Image data !
    integer :: img, img_dim, img_coords(2)

    ! Loop indexes !
    integer :: i, j, k

    ! Lattice arrays !
    real, target, allocatable :: lattice(:,:,:)[:,:], lattice_new(:,:,:)[:,:], density(:,:)[:,:], density_new(:,:)[:,:]

    ! Format for printing the lattice !
    character(len=50) :: format = '(15(f0.2," "),"  ",15(f0.2," "))'

    ! Pointers for ghost nodes !
    real, pointer, dimension(:,:) :: north_ptr, south_ptr
    real, pointer, dimension(:,:) :: east_ptr, west_ptr

    ! Initialize dimensions !
    height = 15
    width = 15
    directions = 9
    img_dim = 2

    allocate(lattice(directions, width, height)[img_dim,*])
    allocate(lattice_new(directions, width, height)[img_dim,*])
    allocate(density(width, height)[img_dim,*])
    allocate(density_new(width, height)[img_dim,*])
    allocate(north_ptr(directions, width))
    allocate(south_ptr(directions, width))
    allocate(east_ptr(directions, height))
    allocate(west_ptr(directions, height))

    img = this_image()
    ! img_coords = this_image(lattice)
    ! print *, "Image ", img, " coordinates: (", img_coords, ")"

    ! Populate lattice !
    call populate_lattice(lattice, width, height, directions)

    ! Populate pointers for ghost nodes !
    ! north_ptr => lattice(:, :, 1)[1, 1]
    ! south_ptr => lattice(:, :, height)[1, 1]
    ! west_ptr  => lattice(:, 1, :)[1, 1]
    ! east_ptr  => lattice(:, width, :)[1, 1]

    ! Perform one streaming step !
    lattice_new = streaming_step(lattice, width, height, directions)

    ! Calculate density arrays !
    density = calculate_density_array(lattice, width, height, directions)
    density_new = calculate_density_array(lattice_new, width, height, directions)

    ! ! Print initial and subsequent lattice !
    ! print *, " "
    ! if (img .eq. 1) then
    !     do i=1, directions
    !         print *, "Direction: ", i-1
    !         print *, "Before streaming step:"
    !         do j=1, height
    !             write(unit=6, fmt=format) lattice(i,:,j)[1,1], lattice(i,:,j)[2,1]
    !         end do
    !         print *, " "
    !         do j=1, height
    !             write(unit=6, fmt=format) lattice(i,:,j)[1,2], lattice(i,:,j)[2,2]
    !         end do
    !         print *, " "
    !         print *, "After streaming step:"
    !         do j=1, height
    !             write(unit=6, fmt=format) lattice_new(i,:,j)[1,1], lattice_new(i,:,j)[2,1]
    !         end do
    !         print *, " "
    !         do j=1, height
    !             write(unit=6, fmt=format) lattice_new(i,:,j)[1,2], lattice_new(i,:,j)[2,2]
    !         end do
    !         print *, " "
    !     end do
    ! end if
    ! sync all

    ! Print initial and subsequent density lattices !
    print *, " "
    if (img .eq. 1) then
        print *, "Before streaming step:"
        do j=1, height
            write(unit=6, fmt=format) density(:,j)[1,1], density(:,j)[2,1]
        end do
        print *, " "
        do j=1, height
            write(unit=6, fmt=format) density(:,j)[1,2], density(:,j)[2,2]
        end do
        print *, " "
        print *, "After streaming step:"
        do j=1, height
            write(unit=6, fmt=format) density_new(:,j)[1,1], density_new(:,j)[2,1]
        end do
        print *, " "
        do j=1, height
            write(unit=6, fmt=format) density_new(:,j)[1,2], density_new(:,j)[2,2]
        end do
        print *, " "
    end if
    sync all

end program SheBoltOnMyManTilILattice