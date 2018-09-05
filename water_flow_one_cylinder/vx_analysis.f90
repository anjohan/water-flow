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
    integer(c_int) :: num_bins = 1000
    real(c_double) :: bin_size = 0

    integer :: u, i, j, fstat = 0
    integer :: x_idx = 3, vx_idx = 6, bin_idx

    namelist /input/ infile, start_timestep, outfile
    read(*,nml=input)

    allocate(vx_hist(num_bins), count_hist(num_bins), atom_data(1,1))
    vx_hist(:) = 0
    count_hist(:) = 0

    open(newunit=u, file=trim(infile), access="stream", action="read", iostat=fstat)
    if (fstat /= 0) then
        error stop "Cannot open " // infile
    end if

    steploop: do
        read(u, iostat=fstat) timestep
        if (fstat == iostat_end) then
            exit steploop
        end if

        read(u) num_atoms, triclinic, boundary_conditions, boundary
        if (triclinic > 0) then
            read(u) angles
        end if

        if (bin_size == 0) then
            bin_size = (boundary(2,1) - boundary(1,1))/num_bins
        end if

        read(u) num_columns, num_chunks

        chunkloop: do i = 1, num_chunks
            read(u) values_in_chunk
            atoms_in_chunk = values_in_chunk/num_columns

            if(size(atom_data,2) < atoms_in_chunk) then
                deallocate(atom_data)
                allocate(atom_data(num_columns, atoms_in_chunk))
            end if

            read(u) atom_data(:,1:atoms_in_chunk)

            if (timestep < start_timestep) cycle chunkloop

            atomloop: do j = 1, atoms_in_chunk
                associate(x  => atom_data(x_idx, j), &
                          vx => atom_data(vx_idx, j))
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
        write(u,*) (i-0.5d0)*bin_size, vx_hist(i)
    end do
    close(u)
end program
