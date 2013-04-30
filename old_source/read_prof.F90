subroutine read_prof(filename,radius,rho,eint,np,ipos)
    implicit none
    integer i, u, cnt
    integer, intent(in) :: np
    integer, intent(out) :: ipos
    real, intent(out) :: radius(np), rho(np), eint(np)
    real :: nradius(np), nrho(np)
    real :: frac(np)
    character(50), intent(in) :: filename
    parameter (u=20)

    open(u, file=filename, status='old')

    read(u,*) ipos
    if (ipos .gt. np) then
        write(*,*) 'Error: ipos = ', ipos ,' is larger than np = ', np
        stop
    endif

    do i=1, ipos
        read(u, *) radius(i), rho(i), frac(i), eint(i)
    enddo
    rho = 10**rho

    nrho(1) = rho(1)
    nradius(1) = radius(1)
    cnt = 1
    do i=2, ipos
        if (rho(i) .ne. rho(i-1) .and. radius(i) .ne. radius(i-1)) then 
            cnt = cnt + 1
            nradius(cnt) = radius(i)
            nrho(cnt) = rho(i)
        endif
    enddo
    ipos = cnt
    radius = nradius
    rho = nrho

    close(u)
    return
end
