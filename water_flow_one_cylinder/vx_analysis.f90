program vx_profile_analysis
    use iso_fortran_env, only: iostat_end, real64, int64
    implicit none

    character(len=1024) :: inbase, Fstring, Ffmt, outbase
    character(len=:), allocatable :: filename, fmtstring

    integer :: i, j, k, l, u

    integer(int64) :: step
    integer :: Nx, Ny, Nz, fstat
    real(real64) :: F, ymin, ymax, zmin, zmax, xmin, xmax,  radius
    real(real64), allocatable :: v(:,:,:,:,:), values(:,:,:,:), Fs(:), &
                                 counts(:,:,:,:)

    real(real64) :: y_fraction, Lx, Ly, Fmin, Fstep
    integer :: numFs, start_timestep, ave_x

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
        if (.not. allocated(counts)) then
            allocate(counts(Nz, Ny, Nx, numFs))
            counts(:,:,:,:) = 0
        end if

        if (.not. allocated(v)) then
            allocate(v(3, Nz, Ny, Nx, numFs))
            v(:,:,:,:,:) = 0
        end if

        steps: do
            read(u, iostat=fstat) step

            if (fstat == iostat_end) exit steps
            if (fstat /= 0) error stop

            read(u, iostat=fstat) values
            if (fstat /= 0) error stop

            if (step < start_timestep) cycle steps

            write(*,"(a,x,i0,x,a,x,i0,x,a)") "Image", this_image(), "reading timestep", &
                                             step, "from " // filename

            counts(:,:,:,l) = counts(:,:,:,l) + values(1,:,:,:)

            do i = 1, Nx
                do j = 1, Ny
                    do k = 1, Nz
                        v(:, k, j, i, l) = v(:, k, j, i, l) + values(1, k, j, i) * values(2:4, k, j, i)
                    end do
                end do
            end do
        end do steps

        do i = 1, Nx
            do j = 1, Ny
                do k = 1, Nz
                if (counts(k,j,i,l) /= 0) v(:, k,j,i,l) = v(:, k,j,i,l) &
                                                           / counts(k,j,i,l)
                end do
            end do
        end do
    end do forces
    call co_sum(v)
    call co_sum(counts)

    deallocate(values)
    if (this_image() /= 1) deallocate(v, counts)

    if (this_image() == 1) then
    block
        real(real64) :: dx, dy ,dz, tmp
        integer :: start_bin_asymm, ymin_idx, ymax_idx
        real(real64), allocatable :: vx_ave_z(:,:,:), vx_ave_yz(:,:), &
                                     vx_front(:,:), vx_rear(:,:), counts_ave_z(:,:,:)

        dx = (xmax - xmin)/Nx
        dy = (ymax - ymin)/Ny
        dz = (zmax - zmin)/Nz

        start_bin_asymm = ceiling(radius/dx)
        ymin_idx = floor((0.5d0 - y_fraction/2)*(ymax-ymin)/dy)
        ymax_idx = floor((0.5d0 + y_fraction/2)*(ymax-ymin)/dy)

        allocate(vx_ave_z(Ny,Nx,numFs), &
                 vx_ave_yz(Nx, numFs))

        counts_ave_z = sum(counts(:,:,:,:), dim=1)

        where (counts_ave_z /= 0)
            vx_ave_z = sum(v(1,:,:,:,:)*counts(:,:,:,:), dim=1)/counts_ave_z
        else where
            vx_ave_z = 0
        end where

        do l = 1, numFs
            do i = 1, Nx
                tmp = sum(counts_ave_z(ymin_idx:ymax_idx,i,l))
                if (tmp == 0) then
                    vx_ave_yz(i, l) = 0
                    cycle
                end if

                vx_ave_yz(i, l) = sum(vx_ave_z(ymin_idx:ymax_idx,i,l) &
                                      * counts_ave_z(ymin_idx:ymax_idx,i,l)) &
                                / tmp
            end do
        end do

        vx_front = vx_ave_yz(Nx/2-start_bin_asymm:2:-1, :)
        vx_rear = vx_ave_yz(Nx/2+1+start_bin_asymm:Nx-1, :)
        open(newunit=u, file=trim(outbase)//"_vx_vs_x.dat", status="replace")
        fmtstring = "('x ',*('F='," // Ffmt // ",'eV/Å',:,x))"
        write(u, fmt=fmtstring) Fs
        do i = 2, Nx-2
            if (i < Nx/2-start_bin_asymm .or. i > Nx/2+start_bin_asymm) then
                write(u, fmt="(*(f0.6,:,x))") (i-1)*dx, (vx_ave_yz(i,j), j = 1, numFs)
            else
                write(u, fmt="(*(f0.6,:,x))") (i-1)*dx, (0.0d0, j = 1, numFs)
            end if
        end do
        close(u)

        open(newunit=u, file=trim(outbase)//"_vx_vs_x_wrapped.dat", status="replace")
        fmtstring = "('x ',*('F='," // Ffmt // ",'eV/Å,infront',x,'F='," &
            // Ffmt // ",'eV/Å,behind',:,x))"
        write(u, fmt=fmtstring) ([Fs(i),Fs(i)] , i = 1, numFs)
        do i = 1, size(vx_front, 1)
            write(u, fmt="(*(f0.6,:,x))") radius+(i-1)*dx, &
                                          (vx_rear(i,j), &
                                           vx_front(i,j), j = 1, numFs)
        end do
        close(u)

        open(newunit=u, file=trim(outbase) // "_asymmetry.dat", status="replace")
        write(u, *) "F[eV/Å] <vx>[Å/ps] relative_2norm_diff relative_integral_diff"
        do i = 1, numFs
            write(u, *) Fs(i), &
                sum(vx_front(:,i) + vx_rear(:,i)) / (size(vx_front(:,i)) + size(vx_rear(:,i))), &
                norm2(vx_front(:,i) - vx_rear(:,i)) / norm2(vx_rear(:,i)), &
                sum(vx_front(:,i) - vx_rear(:,i)) / sum(vx_rear(:,i))
            write(*,*) norm2(vx_rear(:,i)), sum(vx_rear(:,i))
        end do
        close(u)
    end block
    end if
end program
