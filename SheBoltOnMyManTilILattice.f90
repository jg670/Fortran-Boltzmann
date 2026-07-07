program SheBoltOnMyManTilILattice
    use lattice_boltzmann
    implicit none

    ! Image data !
    integer :: total_images, args_count, io_status
    character(len=20) :: width_arg, height_arg, coarray_dim_arg

    ! Lattice arrays !
    real(8), allocatable :: lattice_initial(:,:,:)[:,:]

    total_images = num_images()

    args_count = command_argument_count()

    if (args_count < 3) then
        print *, "Missing arguments"
        stop
    end if

    call get_command_argument(1, width_arg)
    call get_command_argument(2, height_arg)
    call get_command_argument(3, coarray_dim_arg)

    read(width_arg, *, iostat=io_status) global_width
    if(io_status /= 0) then
        print *, "Argument '", trim(width_arg), "' is not a valid integer."
    end if

    read(height_arg, *, iostat=io_status) global_height
    if(io_status /= 0) then
        print *, "Argument '", trim(height_arg), "' is not a valid integer."
    end if

    read(coarray_dim_arg, *, iostat=io_status) coarray_dimensions
    if(io_status /= 0) then
        print *, "Argument '", trim(coarray_dim_arg), "' is not a valid integer."
    end if

    instance_width = global_width / coarray_dimensions + 2
    instance_height = global_height / (total_images / coarray_dimensions) + 2

    allocate(lattice_initial(directions, instance_width, instance_height)[coarray_dimensions,*])

    ! Intial lattice !
    !call populate_lattice_random(lattice_initial)
    !call populate_lattice_dense_center(lattice_initial)
    !call populate_lattice_shear_wave(lattice_initial)
    !call populate_lattice_couette(lattice_initial)
    !call populate_lattice_poiseuille(lattice_initial)
    !call populate_lattice_sliding_lid(lattice_initial)
    call populate_lattice_sliding_lid_parallel(lattice_initial)
    
    sync all

    call output_results_parallel(lattice_initial)

    contains

        subroutine output_results(lattice, interval_length, num_intervals)

            real(8), intent(inout) :: lattice(directions, instance_width, instance_height)[coarray_dimensions,*]
            integer, intent(in) :: interval_length, num_intervals

            integer :: interval, sub_interval, img, j, k
            real(8) :: density_arr(instance_width, instance_height)
            type(velocity) :: velocity_arr(instance_width, instance_height)
            character(:), allocatable :: file_name
            character(len=20) :: interval_str

            !real(8) :: intial_a

            img = this_image()

            if (img .eq. 1) then

                ! Intial lattice !
                density_arr = calculate_density_array(lattice)
                velocity_arr = calculate_average_velocity_array(lattice, density_arr)

                !intial_a = velocity_arr(15,12)%x

                open(1, file="C:\Users\jackg\OneDrive\Desktop\Fortran-Project\Visualization\output-0-sliding-lid.txt", status="replace", action="write")
                    do j=2, instance_width-1
                        do k=2, instance_height-1
                            write(1, *) j-1, ", ", k-1, ", ", density_arr(j,k), ", ", velocity_arr(j,k)%x, ", ", velocity_arr(j,k)%y
                        end do
                    end do
                close(1)

                do interval=1, num_intervals
                    do sub_interval=1, interval_length
                        call perform_one_time_step(lattice)
                    end do

                    density_arr = calculate_density_array(lattice)
                    velocity_arr = calculate_average_velocity_array(lattice, density_arr)

                    !call check_viscosity(intial_a, velocity_arr(15,12)%x, interval * interval_length)
                    !call check_poiseuille_velocity(density_arr(15,12), velocity_arr(15,12)%x, 12)

                    ! Output file !
                    write(interval_str, '(I0)') interval * interval_length
                    file_name = "C:\Users\jackg\OneDrive\Desktop\Fortran-Project\Visualization\output-" // trim(adjustl(interval_str)) // "-sliding-lid.txt"
                    open(1, file=file_name, status="replace", action="write")
                    do j=2, instance_width-1
                        do k=2, instance_height-1
                            write(1, *) j-1, ", ", k-1, ", ", density_arr(j,k), ", ", velocity_arr(j,k)%x, ", ", velocity_arr(j,k)%y
                        end do
                    end do
                    close(1)            
                end do
            end if
            sync all

        end subroutine output_results

        subroutine output_results_parallel(lattice)
            real(8), intent(inout) :: lattice(directions, instance_width, instance_height)[coarray_dimensions,*]

            integer :: img, j
            character(len=50) :: format = '(15(f0.4," "),"  ",15(f0.4," "))'

            img = this_image()

            if (img .eq. 1) then
                write(6,*) "Intial Lattice"
                do j=2,instance_height-1
                    write(unit=6, fmt=format) lattice(2,2:instance_width-1,j)[1,1], lattice(2,2:instance_width-1,j)[2,1]
                end do
                print *, ""
                do j=2,instance_height-1
                    write(unit=6, fmt=format) lattice(2,2:instance_width-1,j)[1,2], lattice(2,2:instance_width-1,j)[2,2]
                end do
                print *, ""
                print *, ""
            end if
            sync all

            call perform_one_time_step(lattice)

            if (img .eq. 1) then
                write(6,*) "Lattice after 1 stream right"
                do j=2,instance_height-1
                    write(unit=6, fmt=format) lattice(2,2:instance_width-1,j)[1,1], lattice(2,2:instance_width-1,j)[2,1]
                end do
                print *, ""
                do j=2,instance_height-1
                    write(unit=6, fmt=format) lattice(2,2:instance_width-1,j)[1,2], lattice(2,2:instance_width-1,j)[2,2]
                end do
                print *, ""
                print *, ""
            end if
            sync all

            call perform_one_time_step(lattice)

            if (img .eq. 1) then
                write(6,*) "Lattice after 2 stream right"
                do j=2,instance_height-1
                    write(unit=6, fmt=format) lattice(2,2:instance_width-1,j)[1,1], lattice(2,2:instance_width-1,j)[2,1]
                end do
                print *, ""
                do j=2,instance_height-1
                    write(unit=6, fmt=format) lattice(2,2:instance_width-1,j)[1,2], lattice(2,2:instance_width-1,j)[2,2]
                end do
                print *, ""
                print *, ""
            end if
            sync all

        end subroutine output_results_parallel

end program SheBoltOnMyManTilILattice