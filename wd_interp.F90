subroutine wd_interp(don_mass, don_pos, acc_pos, noz_pos, don_dist, acc_rate)

    use Simulation_data 
    use gr_mpoleData, ONLY: Xcm, Ycm, Zcm, Mtot
    use Grid_interface, ONLY: Grid_getMinCellSize
    use PhysicalConstants_interface, ONLY : PhysicalConstants_get
    use Driver_data, ONLY: dr_simTime
    implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

    double precision, intent(out) :: don_mass
    double precision, intent(out) :: don_dist, acc_rate
    double precision, intent(out) :: don_pos(3), acc_pos(3), noz_pos(3)
    double precision              :: d12, ind_frac, L1, min_cell_size, G
    double precision              :: noz_thick = NGUARD
    integer                       :: i

    don_mass = sim_donorMass
    call Grid_getMinCellSize(min_cell_size)
    call PhysicalConstants_get("Newton", G)

    do i = 1, sim_TableRows-1
        don_mass = don_mass - (sim_table(i,12)+sim_table(i+1,12))/2.*(sim_table(i+1,1)-sim_table(i,1))
        if (dr_simTime .le. sim_table(i+1,1)) then
            ind_frac = (dr_simTime - sim_table(i,1))/(sim_table(i+1,1) - sim_table(i,1))
            acc_pos(1) = (sim_table(i,6) + (sim_table(i+1,6) - sim_table(i,6))*ind_frac)
            acc_pos(2) = (sim_table(i,7) + (sim_table(i+1,7) - sim_table(i,7))*ind_frac)
            acc_pos(3) = (sim_table(i,8) + (sim_table(i+1,8) - sim_table(i,8))*ind_frac)
            don_pos(1) = Xcm - acc_pos(1) + (sim_table(i,3) + (sim_table(i+1,3) - sim_table(i,3))*ind_frac)
            don_pos(2) = Ycm - acc_pos(2) + (sim_table(i,4) + (sim_table(i+1,4) - sim_table(i,4))*ind_frac)
            don_pos(3) = Zcm - acc_pos(3) + (sim_table(i,5) + (sim_table(i+1,5) - sim_table(i,5))*ind_frac)
            d12 = sqrt((Xcm - don_pos(1))**2. + (Ycm - don_pos(2))**2. + (Zcm - don_pos(3))**2.)
            noz_pos(1) = L1(don_mass, Mtot, d12, 1./sqrt(d12**3. / G / (don_mass + Mtot))) + min_cell_size + don_pos(1)
            noz_pos(2) = Ycm - acc_pos(2)
            noz_pos(3) = Zcm - acc_pos(3)
            don_dist = sqrt(((sim_table(i,3) + (sim_table(i+1,3) - sim_table(i,3))*ind_frac) - acc_pos(1))**2. + &
                ((sim_table(i,4) + (sim_table(i+1,4) - sim_table(i,4))*ind_frac) - acc_pos(2))**2. + &
                ((sim_table(i,5) + (sim_table(i+1,5) - sim_table(i,5))*ind_frac) - acc_pos(3))**2.)
            acc_rate = 2.*(sim_table(i,12) + (sim_table(i+1,12) - sim_table(i,12)) * ind_frac)/noz_thick !Used for thick nozzle
            !write(*,*) 'interp output', i, dr_simTime, sim_table(i+1,1), acc_rate, noz_thick
            exit
        endif
    enddo

end subroutine wd_interp

double precision function L1(m_donor, m_acc, r, omega)
    use PhysicalConstants_interface, ONLY : PhysicalConstants_get
    implicit none

    double precision, intent(in) :: m_donor, m_acc, r, omega
    double precision             :: mu, G, f, df
    double precision, parameter  :: tol = 1.0e-6
    integer                      :: i = 0
    call PhysicalConstants_get('Newton',G)
    !write(*,*) 'm_donor, m_acc, r, omega', m_donor, m_acc, r, omega
    mu = m_acc/(m_acc + m_donor)
    L1 = sqrt(m_donor)/(sqrt(m_donor)+sqrt(m_acc))*r !Ignore centrifugal force to obtain initial guess.
    f = omega**2.*(L1 - r*mu) !We know first term is zero
    !write(*,*) 'fterms', G*(m_acc/(L1-r)**2. - m_donor/L1**2.), omega**2.*(L1 - r*mu)
    !write(*,*) 'G', G
    !write(*,*) 'L1, f', L1, f
    do while (abs(f) .ge. tol)
        df = (-2*G*(m_acc*L1**3. - m_donor*(L1-r)**3.)) / (L1*(L1-r))**3.+omega**2.
        !L1 = L1 - sign(min(abs(f/df), abs(L1/10.)), f/df)
        L1 = L1 - f/df
        f = G*(m_acc/(L1-r)**2. - m_donor/L1**2.) + omega**2.*(L1 - r*mu)
        !write(*,*) 'L1, f, df', L1, f, df
        i = i + 1
    enddo
    return
end function L1
