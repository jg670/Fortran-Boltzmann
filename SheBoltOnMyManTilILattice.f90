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

            !call check_viscosity(lattice, intial_a, velocity_arr(15,12)%x, interval * interval_length)

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

subroutine check_viscosity(lattice, intial_a, new_a, time_step)
    use lattice_boltzmann

    real(8), intent(in) ::  lattice(directions, width, height)
    real(8), intent(in) :: intial_a, new_a
    integer, intent(in) :: time_step

    real(8) :: log_term, calculated_viscosity, k, analytical_viscosity

    k = (2.0_8 * pi) / height
    log_term = log(new_a / intial_a)
    calculated_viscosity = log_term * (-1.0_8 / (k**2 * time_step))

    analytical_viscosity = (1.0_8 / 3.0_8) * ((1.0_8 / omega) - 0.5_8)

    print *, "Calculated viscosity: ", calculated_viscosity
    print *, "Analytical viscosity: ", analytical_viscosity
    print *, "Difference: ", calculated_viscosity - analytical_viscosity


end subroutine check_viscosity