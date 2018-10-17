program vx_profile_analysis
    use iso_fortran_env, only: iostat_end, real64
    implicit none

    character(len=1024) :: inbase, Fstring, Ffmt, outbase
    character(len=:), allocatable :: filename, fmtstring

    integer :: i, j, u

    integer :: timestep, num_chunks, num_atoms, fstat, bin_id, num_columns = 6
    integer :: x_idx = 2, y_idx = 3, count_idx = 5, vx_idx = 6
    real(real64) :: F, ymin, ymax, bin_size, radius
    real(real64), allocatable :: vx_hists(:,:), values(:,:), count_hists(:,:), x(:), Fs(:), &
                                 vx_front(:,:), vx_rear(:,:)

    real(real64) :: y_fraction, Lx, Ly, Fmin, Fstep
    integer :: num_bins, numFs, start_timestep, start_bin_asymm

    namelist /input/ inbase, start_timestep, outbase, y_fraction, radius, num_bins, &
                     Lx, Ly, Fmin, Fstep, numFs, Ffmt

    read(*, nml=input)

    bin_size = Lx/num_bins
    ymin = (0.5-y_fraction/2)*Ly
    ymax = (0.5+y_fraction/2)*Ly
    x = [(bin_size*(i-0.5), i = 1, num_bins)]

    Fs = [(Fmin + i*Fstep, i = 0, numFs-1)]

    allocate(count_hists(num_bins, numFs), vx_hists(num_bins, numFs), values(1,num_columns))
    vx_hists(:,:) = 0
    count_hists(:,:) = 0

    !$omp parallel do default(firstprivate) shared(count_hists, vx_hists)
    forces: do i = 1, numFs
        F = Fs(i)
        Fstring(:) = " "
        write(Fstring, fmt=Ffmt) F

        filename = trim(inbase) // trim(Fstring) // ".profile"

        open(newunit=u, file=filename, status="old")

        ! skip header
        do j = 1, 3
            read(u,*)
        end do

        timesteps: do
            read(u, *, iostat=fstat) timestep, num_chunks, num_atoms

            if (fstat == iostat_end) exit timesteps
            write(*,*) "Reading timestep", timestep, "from ", filename

            if (size(values, 1) < num_chunks) then
                deallocate(values)
                allocate(values(num_columns, num_chunks))
            end if

            read(u,*) values(:,1:num_chunks)


            if (timestep < start_timestep) cycle timesteps

            chunks: do j = 1, num_chunks
                associate(x => values(x_idx,j), y=> values(y_idx, j), &
                          vx => values(vx_idx, j), num_atoms => values(count_idx,j))
                    if (y < ymin .or. y > ymax) cycle chunks

                    bin_id = int(x/bin_size) + 1
                    bin_id = min(bin_id, num_bins)
                    bin_id = max(bin_id, 1)

                    vx_hists(bin_id,i) = vx_hists(bin_id,i) + vx
                    count_hists(bin_id,i) = count_hists(bin_id,i) + num_atoms
                end associate
            end do chunks

        end do timesteps
        close(u)
    end do forces
    !$omp end parallel do

    where (count_hists > 0) vx_hists(:,:) = vx_hists(:,:)/count_hists(:,:)

    vx_front = vx_hists(num_bins/2:1:-1, :)
    vx_rear  = vx_hists(num_bins/2+1:, :)

    open(newunit=u, file=trim(outbase)//"_vx_vs_x.dat", status="replace")
    fmtstring = "('x ',*('F='," // Ffmt // ",'eV/Å',:,x))"
    write(u, fmt=fmtstring) (Fs(i) , i = 1, numFs)
    do i = 1, num_bins
        write(u, fmt="(*(f0.6,:,x))") x(i), (vx_hists(i,j), j = 1, numFs)
    end do
    close(u)

    open(newunit=u, file=trim(outbase)//"_vx_vs_x_wrapped.dat", status="replace")
    fmtstring = "('x ',*('F='," // Ffmt // ",'eV/Å,infront',x,'F='," &
                // Ffmt // ",'eV/Å,behind',:,x))"
    write(u, fmt=fmtstring) ([Fs(i),Fs(i)] , i = 1, numFs)
    do i = 0,num_bins/2-1
        write(u, fmt="(*(f0.6,:,x))") abs(x(num_bins/2+i)-Lx/2), &
                                      (vx_rear(i+1,j), &
                                       vx_front(i+1,j), j = 1, numFs)
    end do
    close(u)

    start_bin_asymm = ceiling(radius/bin_size)

    open(newunit=u, file=trim(outbase) // "_asymmetry.dat", status="replace")
    write(u, *) "F[eV/Å] <vx>[Å/ps] relative_2norm_diff relative_integral_diff"
    do i = 1, numFs
        associate(front => vx_front(start_bin_asymm:,i), &
                  rear => vx_rear(start_bin_asymm:,i))
            write(u, *) Fs(i), &
                        sum(front + rear) / (size(front) + size(rear)), &
                        norm2(front - rear) / norm2(rear), &
                        sum(front - rear) / sum(rear)
        end associate
    end do
    close(u)

end program
