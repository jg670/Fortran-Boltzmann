program SheBoltOnMyManTilILattice
    use lattice_boltzmann
    implicit none

    ! Image data !
    integer :: img, img_dim, img_coords(2)

    ! Loop indexes !
    integer :: i, j, k

    ! Lattice arrays !
    real(8), allocatable :: lattice_initial(:,:,:)[:,:], lattice_1_step(:,:,:)[:,:], lattice_100_steps(:,:,:)[:,:], lattice_100000_steps(:,:,:)[:,:]
    real (8), allocatable :: density_initial(:,:)[:,:], density_1_step(:,:)[:,:], density_100_steps(:,:)[:,:], density_100000_steps(:,:)[:,:]
    type(velocity), allocatable :: average_velocity_initial(:,:)[:,:], average_velocity_1_step(:,:)[:,:], average_velocity_100_steps(:,:)[:,:], average_velocity_100000_steps(:,:)[:,:]

    ! Format for printing the lattice !
    !character(len=50) :: format = '(15(f0.2," "),"  ",15(f0.2," "))'
    character(len=50) :: format = '(15(f0.4," "))'

    img_dim = 1

    allocate(lattice_initial(directions, width, height)[img_dim,*])
    allocate(lattice_1_step(directions, width, height)[img_dim,*])
    allocate(lattice_100_steps(directions, width, height)[img_dim,*])
    allocate(lattice_100000_steps(directions, width, height)[img_dim,*])
    allocate(density_initial(width, height)[img_dim,*])
    allocate(density_1_step(width, height)[img_dim,*])
    allocate(density_100_steps(width, height)[img_dim,*])
    allocate(density_100000_steps(width, height)[img_dim,*])
    allocate(average_velocity_initial(width, height)[img_dim,*])
    allocate(average_velocity_1_step(width, height)[img_dim,*])
    allocate(average_velocity_100_steps(width, height)[img_dim,*])
    allocate(average_velocity_100000_steps(width, height)[img_dim,*])

    img = this_image()
    ! img_coords = this_image(lattice)
    ! print *, "Image ", img, " coordinates: (", img_coords, ")"


    ! Intial lattice !
    !call populate_lattice_random(lattice_initial)
    call populate_lattice_dense_center(lattice_initial)
    density_initial = calculate_density_array(lattice_initial)
    average_velocity_initial = calculate_average_velocity_array(lattice_initial, density_initial)    

    ! Print initial lattice !
    if (img .eq. 1) then
        print *, " "
        print *, "Initial lattice:"
        do k=1, height
            write(unit=6, fmt=format) density_initial(:,k)
        end do
        print *, " "
        print *, "Initial density: ", sum(lattice_initial)

        ! Output file !
        open(1, file="C:\Users\jackg\OneDrive\Desktop\Fortran-Project\Visualization\output-0.txt", status="replace", action="write")
        do j=1, width
            do k=1, height
                write(1, *) j, ", ", k, ", ", density_initial(j,k), ", ", average_velocity_initial(j,k)%x, ", ", average_velocity_initial(j,k)%y
            end do
        end do
        close(1)
    end if
    sync all


    ! Lattice after 1 step !
    lattice_1_step = streaming_step(lattice_initial)
    lattice_1_step = collision_step(lattice_1_step)
    density_1_step = calculate_density_array(lattice_1_step)
    average_velocity_1_step = calculate_average_velocity_array(lattice_1_step, density_1_step)

    ! Print lattice after 1 step !
    if (img .eq. 1) then
        print *, " "
        print *, "Lattice after 1 step:"
        do k=1, height
            write(unit=6, fmt=format) density_1_step(:,k)
        end do
        print *, " "
        print *, "Lattice density after 1 step: ", sum(lattice_1_step)

        ! Output file !
        open(1, file="C:\Users\jackg\OneDrive\Desktop\Fortran-Project\Visualization\output-1.txt", status="replace", action="write")
        do j=1, width
            do k=1, height
                write(1, *) j, ", ", k, ", ", density_1_step(j,k), ", ", average_velocity_1_step(j,k)%x, ", ", average_velocity_1_step(j,k)%y
            end do
        end do
        close(1)
    end if
    sync all


    ! Lattice after 100 steps !
    lattice_100_steps = lattice_1_step
    do i=1, 99
        lattice_100_steps = streaming_step(lattice_100_steps)
        lattice_100_steps = collision_step(lattice_100_steps)
    end do
    density_100_steps = calculate_density_array(lattice_100_steps)
    average_velocity_100_steps = calculate_average_velocity_array(lattice_100_steps, density_100_steps)

    ! Print lattice after 100 steps !
    if (img .eq. 1) then
        print *, " "
        print *, "Lattice after 100 steps"
        do k=1, height
            write(unit=6, fmt=format) density_100_steps(:,k)
        end do
        print *, " "
        print *, "Lattice density after 100 steps: ", sum(lattice_100_steps)

        ! Output file !
        open(1, file="C:\Users\jackg\OneDrive\Desktop\Fortran-Project\Visualization\output-100.txt", status="replace", action="write")
        do j=1, width
            do k=1, height
                write(1, *) j, ", ", k, ", ", density_100_steps(j,k), ", ", average_velocity_100_steps(j,k)%x, ", ", average_velocity_100_steps(j,k)%y
            end do
        end do
        close(1)
    end if
    sync all


    ! Lattice after 100000 steps !
    lattice_100000_steps = lattice_100_steps
    do i=1, 99900
        lattice_100000_steps = streaming_step(lattice_100000_steps)
        lattice_100000_steps = collision_step(lattice_100000_steps)
    end do
    density_100000_steps = calculate_density_array(lattice_100000_steps)
    average_velocity_100000_steps = calculate_average_velocity_array(lattice_100000_steps, density_100000_steps)

    ! Print lattice after 100000 steps !
    if (img .eq. 1) then
        print *, " "
        print *, "Lattice after 100000 steps: "
        do k=1, height
            write(unit=6, fmt=format) density_100000_steps(:,k)
        end do
        print *, " "
        print *, "Lattice density after 100000 steps: ", sum(lattice_100000_steps)

        ! Output file !
        open(1, file="C:\Users\jackg\OneDrive\Desktop\Fortran-Project\Visualization\output-100000.txt", status="replace", action="write")
        do j=1, width
            do k=1, height
                write(1, *) j, ", ", k, ", ", density_100000_steps(j,k), ", ", average_velocity_100000_steps(j,k)%x, ", ", average_velocity_100000_steps(j,k)%y
            end do
        end do
        close(1)
    end if
    sync all

end program SheBoltOnMyManTilILattice