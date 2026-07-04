program SheBoltOnMyManTilILattice
    use lattice_boltzmann
    implicit none

    ! Image data !
    integer :: img, img_dim !, img_coords(2)


    ! Lattice arrays !
    real(8), allocatable :: lattice_initial(:,:,:)[:,:]
    
    ! Format for printing the lattice !
    ! character(len=50) :: format = '(15(f0.2," "),"  ",15(f0.2," "))'
    ! character(len=50) :: format = '(15(f0.4," "))'

    img_dim = 1

    allocate(lattice_initial(directions, width, height)[img_dim,*])

    ! img = this_image()
    ! img_coords = this_image(lattice)
    ! print *, "Image ", img, " coordinates: (", img_coords, ")"


    ! Intial lattice !
    !call populate_lattice_random(lattice_initial)
    !call populate_lattice_dense_center(lattice_initial)
    !call populate_lattice_shear_wave(lattice_initial)
    !call populate_lattice_couette(lattice_initial)
    call populate_lattice_poiseuille(lattice_initial)
    
    call output_results(lattice_initial, 25, 1000)

end program SheBoltOnMyManTilILattice

subroutine output_results(lattice, interval_length, num_intervals)
    use lattice_boltzmann

    real(8), intent(inout) :: lattice(directions, width, height)
    integer, intent(in) :: interval_length, num_intervals

    integer :: interval, sub_interval, img
    real(8) :: density_arr(width, height)
    type(velocity) :: velocity_arr(width, height)
    character(:), allocatable :: file_name
    character(len=20) :: interval_str

    !real(8) :: intial_a

    img = this_image()

    if (img .eq. 1) then

        ! Intial lattice !
        density_arr = calculate_density_array(lattice)
        velocity_arr = calculate_average_velocity_array(lattice, density_arr)

        !intial_a = velocity_arr(15,12)%x

        open(1, file="C:\Users\jackg\OneDrive\Desktop\Fortran-Project\Visualization\output-0-poiseuille.txt", status="replace", action="write")
            do j=2, width-1
                do k=2, height-1
                    write(1, *) j-1, ", ", k-1, ", ", density_arr(j,k), ", ", velocity_arr(j,k)%x, ", ", velocity_arr(j,k)%y
                end do
            end do
        close(1)

        do interval=1, num_intervals
            do sub_interval=1, interval_length
                lattice = perform_one_time_step(lattice)
            end do

            density_arr = calculate_density_array(lattice)
            velocity_arr = calculate_average_velocity_array(lattice, density_arr)

            !call check_viscosity(intial_a, velocity_arr(15,12)%x, interval * interval_length)
            !call check_poiseuille_velocity(density_arr(15,12), velocity_arr(15,12)%x, 12)

            ! Output file !
            write(interval_str, '(I0)') interval * interval_length
            file_name = "C:\Users\jackg\OneDrive\Desktop\Fortran-Project\Visualization\output-" // trim(adjustl(interval_str)) // "-poiseuille.txt"
            open(1, file=file_name, status="replace", action="write")
            do j=2, width-1
                do k=2, height-1
                    write(1, *) j-1, ", ", k-1, ", ", density_arr(j,k), ", ", velocity_arr(j,k)%x, ", ", velocity_arr(j,k)%y
                end do
            end do
            close(1)            
        end do
    end if
    sync all

end subroutine output_results