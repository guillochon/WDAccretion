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

    implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

    integer, intent(in) :: myPE
    integer             :: i, ierr
    double precision    :: start_t

    sim_pi = PI
    call RuntimeParameters_get('sim_pAmbient', sim_pAmbient)
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

    if (sim_nSubZones .le. 1) sim_nSubZones = 2

    sim_inSubZones = 1./real(sim_nSubZones)
    sim_inSubzm1   = 1./real(sim_nSubZones-1)
    sim_inszd      = sim_inSubZones**NDIM

    if (sim_tInitial + sim_tRelax .lt. sim_tNozz) sim_tRelax = sim_tNozz - sim_tInitial
    if (myPE .eq. MASTER_PE) then
        call read_prof('WD_hires_helm_prof_no_coul_0.90_CaO.dat', radius, rhop, eint, np, ipos)
        call read_table_dims('quantities_interp.dat', sim_tableRows, sim_tableCols)
        allocate(sim_table(sim_tableRows, sim_tableCols))
        call read_table('quantities_interp.dat', sim_table, sim_tableRows, sim_tableCols)
        start_t = sim_table(1,1)
        sim_table(:,1) = sim_table(:,1) - start_t + sim_tNozz
        sim_table(:,2:11) = 1.0d9*sim_table(:,2:11)
        sim_table(:,12) = 2.0d33*sim_table(:,12)
    endif

    radius(ipos + 1) = 2*radius(ipos) - radius(ipos - 1)
    rhop(ipos + 1) = sim_rhoAmbient
    ipos = ipos + 1

    call MPI_BCAST(radius, np, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(rhop, np, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(eint, np, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(ipos, 1, FLASH_INTEGER, MASTER_PE, MPI_COMM_WORLD, ierr)

    call MPI_BCAST(sim_tableRows, 1, FLASH_INTEGER, MASTER_PE, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(sim_tableCols, 1, FLASH_INTEGER, MASTER_PE, MPI_COMM_WORLD, ierr)
    if (myPE .ne. MASTER_PE) allocate(sim_table(sim_tableRows, sim_tableCols))
    call MPI_BCAST(sim_table, sim_tableRows*sim_tableCols, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)

end subroutine Simulation_init
