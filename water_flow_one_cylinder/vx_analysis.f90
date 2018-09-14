program vx_analysis
    use iso_c_binding, only: c_int, c_long_long, c_double
    use iso_fortran_env, only: iostat_end
    implicit none
    character(len=1024) :: infile, outfile
    integer(c_int) :: num_chunks, num_columns, values_in_chunk, atoms_in_chunk
    integer(c_long_long) :: num_atoms, timestep, start_timestep

    integer(c_int) :: triclinic, boundary_conditions(2,3)
    real(c_double) :: boundary(2,3), angles(3)

    real(c_double), allocatable :: atom_data(:,:), vx_hist(:)
    integer(c_int), allocatable :: count_hist(:)
    integer(c_int) :: num_bins
    real(c_double) :: bin_size = 0, y_fraction, Lx, Ly, ymin, ymax

    integer :: u, i, j, fstat = 0
    integer :: x_idx, vx_idx, bin_idx

    namelist /input/ infile, start_timestep, outfile, y_fraction, num_bins, x_idx, vx_idx
    read(*,nml=input)

    allocate(vx_hist(num_bins), count_hist(num_bins), atom_data(1,1))
    vx_hist(:) = 0
    count_hist(:) = 0

    open(newunit=u, file=trim(infile), access="stream", action="read", iostat=fstat)
    if (fstat /= 0) then
        write(*,*) "Cannot open " // infile
        error stop
    end if

    steploop: do
        read(u, iostat=fstat) timestep, num_atoms, triclinic, boundary_conditions, boundary
        if (fstat == iostat_end) then
            exit steploop
        end if

        write(*,*) "Reading timestep", timestep

        if (triclinic > 0) then
            read(u, iostat=fstat) angles
            if (fstat == iostat_end) then
                exit steploop
            end if
        end if

        if (bin_size == 0) then
            Lx = boundary(2,1) - boundary(1,1)
            Ly = boundary(2,2) - boundary(1,2)
            bin_size = Lx/num_bins
            ymin = (0.5-y_fraction/2)*Ly
            ymax = (0.5+y_fraction/2)*Ly
        end if

        read(u, iostat=fstat) num_columns, num_chunks
        if (fstat == iostat_end) then
            exit steploop
        end if

        chunkloop: do i = 1, num_chunks
            read(u, iostat=fstat) values_in_chunk
            if (fstat == iostat_end) then
                exit steploop
            end if
            atoms_in_chunk = values_in_chunk/num_columns

            if(size(atom_data,2) < atoms_in_chunk) then
                deallocate(atom_data)
                allocate(atom_data(num_columns, atoms_in_chunk))
            end if

            read(u, iostat=fstat) atom_data(:,1:atoms_in_chunk)
            if (fstat == iostat_end) then
                exit steploop
            end if

            if (timestep < start_timestep) cycle chunkloop

            atomloop: do j = 1, atoms_in_chunk
                associate(x  => atom_data(x_idx, j), &
                          y  => atom_data(x_idx+1, j), &
                          vx => atom_data(vx_idx, j))
                    if (y < ymin .or. y > ymax) cycle atomloop
                    bin_idx = int(x/bin_size) + 1
                    bin_idx = min(bin_idx, num_bins)
                    bin_idx = max(bin_idx, 1)
                    count_hist(bin_idx) = count_hist(bin_idx) + 1
                    vx_hist(bin_idx) = vx_hist(bin_idx) + vx
                end associate
            end do atomloop
        end do chunkloop
    end do steploop
    where(count_hist > 0)
        vx_hist = vx_hist/count_hist
    end where
    close(u)

    open(newunit=u, file=trim(outfile), status="replace")
    do i = 1, num_bins
        write(u,*) (i-0.5d0)*bin_size, vx_hist(i), abs((i-0.5d0)*bin_size-Lx/2)
    end do
    close(u)
end program
