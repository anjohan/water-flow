program vx_profile_analysis
    use iso_fortran_env, only: iostat_end, real64, int64
    implicit none

    character(len=1024) :: inbase, Fstring, Ffmt, outbase
    character(len=:), allocatable :: filename, fmtstring

    integer :: i, j, k, l, u

    integer(int64) :: step
    integer(int64), allocatable :: steps(:)
    integer :: Nx, Ny, Nz, fstat
    real(real64) :: F, ymin, ymax, zmin, zmax, xmin, xmax,  radius, &
                    final_v_sum2, final_v_sum
    real(real64), allocatable :: values(:,:,:,:), Fs(:), &
                                 final_v(:,:,:,:), correlations(:,:), rel_diffs(:,:)

    real(real64) :: y_fraction, Lx, Ly, Fmin, Fstep
    integer :: numFs, start_timestep, ave_x, num_steps

    namelist /input/ inbase, start_timestep, outbase, y_fraction, radius, &
                     Fmin, Fstep, numFs, Ffmt, ave_x

    if (this_image() == 1) then
        read(*, nml=input)
    end if

    call co_broadcast(inbase, 1)
    call co_broadcast(start_timestep, 1)
    call co_broadcast(outbase, 1)
    call co_broadcast(y_fraction, 1)
    call co_broadcast(radius, 1)
    call co_broadcast(Fmin, 1)
    call co_broadcast(Fstep, 1)
    call co_broadcast(numFs, 1)
    call co_broadcast(Ffmt, 1)
    call co_broadcast(ave_x, 1)

    Fs = [(Fmin + i*Fstep, i = 0, numFs-1)]

    forces: do l = 1, numFs
        if (mod(l-1, num_images()) /= this_image()-1) cycle forces
        F = Fs(l)
        Fstring(:) = " "
        write(Fstring, fmt=Ffmt) F

        filename = trim(inbase) // trim(Fstring) // ".bin"

        open(newunit=u, file=filename, access="stream", status="old")

        read(u) Nx, Ny, Nz, xmin, xmax, ymin, ymax, zmin, zmax
        if (.not. allocated(values)) allocate(values(4, Nz, Ny, Nx))

        num_steps = 0
        timesteps: do
            read(u, iostat=fstat) step

            write(*,"(a,x,i0,x,a,x,i0,x,a)") "Image", this_image(), "reading timestep", &
                                             step, "from " // filename

            if (fstat == iostat_end) exit timesteps
            if (fstat /= 0) error stop

            read(u, iostat=fstat) values
            if (fstat /= 0) error stop
            num_steps = num_steps+1
        end do timesteps

        final_v = values(2:4,:,:,:)
        final_v_sum = sum(final_v)
        final_v_sum2 = sum(final_v**2)

        if (.not. allocated(correlations)) then
            allocate(correlations(num_steps, numFs))
            correlations(:,:) = 0
            allocate(rel_diffs, source=correlations)
        end if

        if (.not. allocated(steps)) allocate(steps(num_steps))

        close(u)

        open(newunit=u, file=filename, access="stream", status="old")

        read(u) Nx, Ny, Nz, xmin, xmax, ymin, ymax, zmin, zmax

        do i = 1, num_steps
            read(u) steps(i), values

            write(*,"(a,x,i0,x,a,x,i0,x,a)") "Image", this_image(), "reading timestep", &
                                             steps(i), "from " // filename

            correlations(i,l) = sum(final_v(:,:,:,:)*values(2:4,:,:,:))/final_v_sum2
            rel_diffs(i,l) = sum(final_v - values(2:4,:,:,:))/final_v_sum
            if (this_image() == 1) write(*,*) correlations(i,l)
        end do

        close(u)

    end do forces
    call co_sum(correlations)
    call co_sum(rel_diffs)

    deallocate(values)

    if (this_image() == 1) then
        open(newunit=u, file=trim(outbase)//"_corr_final.dat", status="replace")
        fmtstring = "('step ',*('F='," // Ffmt // ",'eV/Å',:,x))"
        write(u, fmt=fmtstring) Fs

        do i = 1, num_steps
            write(u, *) steps(i), correlations(i,:)
        end do

        close(u)

        open(newunit=u, file=trim(outbase)//"_rel_diff.dat", status="replace")
        write(u, fmt=fmtstring) Fs

        do i = 1, num_steps
            write(u, *) steps(i), rel_diffs(i,:)
        end do

        close(u)
    end if
end program
