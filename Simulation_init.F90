!!****if* source/Simulation/SimulationMain/WDAccrection/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init(integer(IN) :: myPE)
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get rotuine for initialization.
!!
!! ARGUMENTS
!!
!!   myPE -   current processor number
!!
!! PARAMETERS
!!
!!  sim_pAmbient       Initial ambient pressure
!!  sim_rhoAmbient     Initial ambient density
!!  sim_xctr           Explosion center coordinates
!!  sim_yctr           Explosion center coordinates
!!  sim_zctr           Explosion center coordinates
!!  sim_nsubzones      Number of `sub-zones' in cells for applying 1d profile
!!
!!***

subroutine Simulation_init(myPE)

    use Simulation_data 
    use RuntimeParameters_interface, ONLY : RuntimeParameters_get
    use Logfile_interface
    use Grid_interface, ONLY : Grid_getMinCellSize

    implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

    integer, intent(in) :: myPE
    integer             :: i, j, k, ierr, istat, ii, jj, kk, gdim, jLo, jHi, n
    double precision    :: start_t, rho0_cgs, massacc, wd_mass_tot, min_grid, sumRho, don_mass, don_dist
    double precision    :: xDist, yDist, zDist, xCoord, yCoord, zCoord, frac, grid_center, dist
    double precision    :: don_pos(3), acc_pos(3), noz_pos(3), acc_rate
    character(len=32), dimension(1,2) :: block_buff
    character(len=32) :: int_to_str
    double precision,allocatable,dimension(:,:,:) :: grid_3d
    logical             :: calc_3d = .false.

    sim_pi = PI
    call RuntimeParameters_get('sim_tAmbient', sim_tAmbient)
    call RuntimeParameters_get('sim_rhoAmbient', sim_rhoAmbient)
    call RuntimeParameters_get('sim_nsubzones',sim_nSubZones)
    call RuntimeParameters_get('sim_xctr',sim_xCenter)
    call RuntimeParameters_get('sim_yctr',sim_yCenter)
    call RuntimeParameters_get('sim_zctr',sim_zCenter)
    call RuntimeParameters_get('smallx', sim_smallX)
    call RuntimeParameters_get('smlrho', sim_smallRho)
    call RuntimeParameters_get('smallp', sim_smallP)
    call RuntimeParameters_get('xmin',sim_xMin)
    call RuntimeParameters_get('ymin',sim_yMin)
    call RuntimeParameters_get('zmin',sim_zMin)
    call RuntimeParameters_get('xmax',sim_xMax)
    call RuntimeParameters_get('ymax',sim_yMax)
    call RuntimeParameters_get('zmax',sim_zMax)
    call RuntimeParameters_get('tinitial',sim_tInitial)
    call RuntimeParameters_get('sim_tNozz',sim_tNozz)
    call RuntimeParameters_get('sim_tRelax',sim_tRelax)
    call RuntimeParameters_get('sim_relaxRate',sim_relaxRate)
    call RuntimeParameters_get('sim_accMass',sim_accMass)
    call RuntimeParameters_get('sim_accTemp',sim_accTemp)
    call RuntimeParameters_get('sim_donorMass',sim_donorMass)
    call RuntimeParameters_get('sim_donorTemp',sim_donorTemp)
    call RuntimeParameters_get('sim_dataFile',sim_dataFile)

    if (sim_nSubZones .le. 1) sim_nSubZones = 2

    sim_inSubZones = 1./real(sim_nSubZones)
    sim_inSubzm1   = 1./real(sim_nSubZones-1)
    sim_inszd      = sim_inSubZones**NDIM

    !if (sim_tInitial + sim_tRelax .lt. sim_tNozz) sim_tRelax = sim_tNozz - sim_tInitial
    if (myPE .eq. MASTER_PE) then
        !Read in evolution data.
        call read_table_dims(sim_dataFile, sim_tableRows, sim_tableCols)
        allocate(sim_table(sim_tableRows, sim_tableCols))
        call read_table(sim_dataFile, sim_table, sim_tableRows, sim_tableCols)
        start_t = sim_table(1,1)
        sim_table(:,1) = sim_table(:,1) - start_t + sim_tNozz
        sim_table(:,2:11) = 1.0d9*sim_table(:,2:11)
        !sim_table(:,12) = 2.0d33*sim_table(:,12)

        !Initialize accretor WD
        wd_xn = sim_smallX
        wd_xn(C12_SPEC) = 0.5
        wd_xn(O16_SPEC) = 0.5

        !how precise do we want the mass?
        massacc = 1.d-5

        !guess central density (cgs) 
        call Logfile_stampMessage(myPE,'[Simulation Init] Guessing Central Density')
        call wd_guess_rho_0(sim_accMass,rho0_cgs)
        !print *, 'rho0_cgs', rho0_cgs

        call wd_interp(don_mass, don_pos, acc_pos, noz_pos, don_dist, acc_rate)

        write(*,*) 'don_dist, don_mass', don_dist, don_mass
        !loop until mass is correct
        wd_mass_tot = 0.0
        n = 0
        call Grid_getMinCellSize(min_grid)
        do while (abs(wd_mass_tot-sim_accMass)/sim_accMass.gt.massacc)
            if (n .gt. 0) rho0_cgs = rho0_cgs * 10**((log10(sim_accMass)-log10(wd_mass_tot)))
            !print *, 'rho0_cgs', rho0_cgs
            write (int_to_str, '(i11)') n
            call Logfile_stampMessage(myPE,'[Simulation Init] Generating Accretor Profile, Iteration #' // &
                trim(adjustl(int_to_str)))

            call wd_intout(wd_xn,rho0_cgs,sim_accTemp,wd_mass_tot,ipos,radius,rhop,don_dist,don_mass)

            write (int_to_str, '(e10.4)') radius(ipos)
            call Logfile_stampMessage(myPE,'Radius: ' // trim(adjustl(int_to_str)))
            write (int_to_str, '(e10.4)') wd_mass_tot
            call Logfile_stampMessage(myPE,'Mass: ' // trim(adjustl(int_to_str)))
            write (int_to_str, '(i10)') ipos
            call Logfile_stampMessage(myPE,'Profile length: ' // trim(adjustl(int_to_str)))
            write (int_to_str, '(e10.4)') abs(wd_mass_tot-sim_accMass)/sim_accMass
            call Logfile_stampMessage(myPE,'Mass Error: ' // trim(adjustl(int_to_str)))
            write (int_to_str, '(e10.4)') rhop(1)
            call Logfile_stampMessage(myPE,'Central Density: ' // trim(adjustl(int_to_str)))

            if (n .gt. 10000) then
                print *, 'wd mass did not converge!'
                stop
            endif
            n = n + 1
        enddo

        massacc = 1.d-3
        wd_mass_tot = 0.0
        n = 0
        do while (abs(wd_mass_tot-sim_accMass)/sim_accMass.gt.massacc)
            if (n .gt. 0) rho0_cgs = rho0_cgs * 10**((log10(sim_accMass)-log10(wd_mass_tot)))
            !print *, 'rho0_cgs', rho0_cgs
            write (int_to_str, '(i11)') n
            call Logfile_stampMessage(myPE,'[Simulation Init] Generating Accretor in 3D, Iteration #' // &
                trim(adjustl(int_to_str)))

            call wd_intout(wd_xn,rho0_cgs,sim_accTemp,wd_mass_tot,ipos,radius,rhop,don_dist,don_mass)

            gdim = 2*ceiling(radius(ipos)/min_grid)
            grid_center = gdim/2*min_grid
            allocate(grid_3d(0:gdim-1,0:gdim-1,0:gdim-1),stat=istat)

            do k = 0, gdim-1
                zCoord = min_grid*k
                do j = 0, gdim-1
                    yCoord = min_grid*j
                    do i = 0, gdim-1
                        xCoord = min_grid*i
                        sumRho = 0.0
                        do kk = 0, sim_nSubZones-1
                           zDist = zCoord + (kk*sim_inSubzm1)*min_grid - grid_center
                           do jj = 0, sim_nSubZones-1
                              yDist = yCoord + (jj*sim_inSubzm1)*min_grid - grid_center
                              do ii = 0, sim_nSubZones-1
                                 xDist = xCoord + (ii*sim_inSubzm1)*min_grid - grid_center
                                 dist = sqrt(xDist**2 + yDist**2 + zDist**2)
                                 if (dist .gt. radius(ipos)) cycle
                                 call sim_find(radius, ipos, dist, jLo)
                                 if (jLo .eq. 0) then
                                    jLo = 1
                                    jHi = 1
                                    frac = 0.
                                 else if (jLo .eq. ipos) then
                                    jHi = ipos
                                    frac = 0.
                                 else
                                    jHi = jLo + 1
                                    frac = (dist - radius(jLo)) / & 
                                         (radius(jHi)-radius(jLo))
                                 endif
                                 sumRho = sumRho + rhop(jLo) + frac*(rhop(jHi) - rhop(jLo))
                              enddo
                           enddo
                        enddo
                        grid_3d(i,j,k) = sumRho
                    enddo
                enddo
            enddo
            wd_mass_tot = sum(grid_3d)*min_grid**3.0*sim_inszd
            deallocate(grid_3d)

            write (int_to_str, '(e10.4)') radius(ipos)
            call Logfile_stampMessage(myPE,'Radius: ' // trim(adjustl(int_to_str)))
            write (int_to_str, '(e10.4)') wd_mass_tot
            call Logfile_stampMessage(myPE,'Mass: ' // trim(adjustl(int_to_str)))
            write (int_to_str, '(i10)') ipos
            call Logfile_stampMessage(myPE,'Profile length: ' // trim(adjustl(int_to_str)))
            write (int_to_str, '(e10.4)') abs(wd_mass_tot-sim_accMass)/sim_accMass
            call Logfile_stampMessage(myPE,'Mass Error: ' // trim(adjustl(int_to_str)))
            write (int_to_str, '(e10.4)') rhop(1)
            call Logfile_stampMessage(myPE,'Central Density: ' // trim(adjustl(int_to_str)))

            if (n .gt. 10000) then
                print *, 'wd mass did not converge!'
                stop
            endif
            n = n + 1
        enddo

    endif

    call MPI_BCAST(radius, np, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(rhop, np, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(eint, np, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(wd_xn, NSPECIES, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(ipos, 1, FLASH_INTEGER, MASTER_PE, MPI_COMM_WORLD, ierr)

    wd_radius = radius(ipos)

    call MPI_BCAST(sim_tableRows, 1, FLASH_INTEGER, MASTER_PE, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(sim_tableCols, 1, FLASH_INTEGER, MASTER_PE, MPI_COMM_WORLD, ierr)
    if (myPE .ne. MASTER_PE) allocate(sim_table(sim_tableRows, sim_tableCols))
    call MPI_BCAST(sim_table, sim_tableRows*sim_tableCols, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)

end subroutine Simulation_init
