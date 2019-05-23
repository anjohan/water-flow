program analysis
    use iso_fortran_env, only: real64

    implicit none

    integer, parameter :: numFs = 10, num_bins_x = 400, num_bins_tot = 400*400*50

    real(real64), parameter :: a = 7.16d0, r = 4*a, Ly=30*a, ymin=0.45d0*Ly, ymax=0.55d0*Ly

    real(real64) :: Fs(numFs), vx(num_bins_x, numFs)[*], v(6, num_bins_tot)

    integer :: vx_counts(num_bins_x)

    character(len=1024) :: fname

    integer :: iu, ou, i, Fi, Li, xi
    real(real64) :: Lx, dx, cyl_xmin, cyl_xmax


    Fs = [(0.00001d0*i, i = 1, numFs)]
    vx(:,:) = 0

    vx_counts(:) = 0

    lengths: do Li = 1, 1

        Lx = 30*a*Li
        dx = Lx/num_bins_x

        cyl_xmin = Lx/2 - r
        cyl_xmax = Lx/2 + r

        forces: do Fi = 1, numFs

            if (mod(Fi-1, num_images()) /= this_image()-1) cycle forces

            fname(:) = " "
            write(fname, '("./continuum_results/length_index_",i0,"_force_index_",i0)') Li, Fi

            open(newunit=iu, file=trim(fname), action="read")

            write(*,*) "Image ", this_image(), " reading ", trim(fname)

            ! skip header
            do i = 1, 9
                read(iu, *)
            end do

            read(iu, *) v
            close(iu)

            ! Ångstrøm
            v(1:3, :) = v(1:3, :) * 1.0d10

            vx_counts(:) = 0

            do i = 1, num_bins_tot
                associate(x => v(1,i), y => v(2,i))
                    if ((x < cyl_xmin .or. x > cyl_xmax) .and. y > ymin .and. y < ymax) then
                        xi = floor(x / dx) + 1
                        vx_counts(xi) = vx_counts(xi) + 1
                        vx(xi, Fi) = vx(xi, Fi) + v(4, i)
                    end if
                end associate
            end do

            where (vx_counts > 0)
                vx(:, Fi) = vx(:, Fi) / vx_counts(:)
            end where

        end do forces

        sync all
        call co_sum(vx)
        if (this_image() /= 1) cycle lengths

        fname(:) = " "
        write(fname, '("/home/anders/master/data/continuum/lx",i0,".dat")') 30*Li

        open(newunit=ou, file=trim(fname), action="write")

        write(ou, '("x",*(:,x,"F=0",f0.5,"eV/Å"))') Fs

        do i = 1, num_bins_x
            write(ou, "(*(:,x,f0.6))") (i-1)*dx, vx(i,:)
        end do

        close(ou)
    end do lengths

end program
