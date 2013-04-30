!!****if* source/Driver/DriverMain/Driver_sourceTerms
!!
!! NAME
!!
!!  Driver_sourceTerms
!!
!! SYNOPSIS
!!
!!  Driver_sourceTerms(integer(IN)::blockCount,
!!                     integer(IN)::blockList(blockCount),
!!                     real(IN) :: dt)
!!
!! DESCRIPTION
!!
!!  Driver for source terms. Instead of calling all these routines 
!!  from Driver_evolveFlash we call Driver_sourceTerms which then
!!  makes the calls to Cool, Burn, Heat and Stir.  If a unit is not
!!  included in the simulation, the routine will be a stub and return
!!  without doing anything.
!! 
!!
!! ARGUMENTS
!!  blockCount   : The number of blocks in the list
!!  blockList    : The list of blocks on which to apply the stirring operator
!!  dt           : the current timestep
!!
!!***



subroutine Driver_sourceTerms(blockCount, blockList, dt) 

    use Driver_data, ONLY: dr_simTime, dr_initialSimTime
    !use Flame_interface, ONLY : Flame_step
    use Stir_interface, ONLY : Stir
    use Heat_interface, ONLY : Heat
    use Burn_interface, ONLY : Burn
    use Cool_interface, ONLY : Cool
    use Simulation_data, ONLY: sim_smallX, sim_rhoAmbient, sim_tRelax, wd_radius, sim_accTemp, sim_donorTemp, rhop, &
        sim_tInitial, ipos, sim_accMass, sim_relaxRate
    use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr,&
        Grid_getCellCoords, Grid_putPointData, Grid_getMinCellSize
    use Eos_interface, ONLY : Eos_wrapped, Eos
    use gr_mpoleData, ONLY: Xcm, Ycm, Zcm, Mtot
    use PhysicalConstants_interface, ONLY : PhysicalConstants_get
    use RuntimeParameters_interface, ONLY : RuntimeParameters_mapStrToInt, RuntimeParameters_get
    implicit none
#include "Eos.h"
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  
    integer, parameter :: ngp = 10000
    real, intent(IN)    :: dt
    integer, intent(IN) :: blockCount
    integer, dimension(blockCount), intent(IN):: blockList
    
    integer  ::  i, j, k, l, put, lb, ierr
    real     ::  xx, yy, zz, local_cell_size, min_cell_size, c2, rnozz, cell_dens
    real     ::  nnozz, G, falling_time, local_xnozz, relax_rate
    real     ::  new_mass, mass_sum, xcenter, ycenter, zcenter
    integer  ::  istat
  
    real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
    integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
    integer :: sizeX,sizeY,sizeZ
    real, dimension(SPECIES_BEGIN:SPECIES_END) :: xn
    real, dimension(:,:,:,:),pointer :: solnData
    real, dimension(:,:,:,:),allocatable :: velCopy
    real, dimension(EOS_NUM) :: eosData
    real, dimension(MDIM) :: avg_vel, new_vel, tot_avg_vel
    integer,dimension(MDIM) :: axis
  
    logical :: gcell = .true.
    logical :: changedGrid
  
    !nozzle variables
    integer cnt,ri,tot_cnt
    real :: rho(0:ngp),dr,delt,yzdist,theta
    real :: m_nozz_tot,r(0:ngp),rho0_cgs,massacc,mrate,mrate_adj,accel_frac
    real :: rhoi, don_pos(3), acc_pos(3), noz_pos(3)
    
    real :: ke, don_mass, don_dist, adj_cfactor, d12, L1
    real :: donor_size = 1.0d9
    real :: donor_dens = 2.0d2
    real :: noz_thick = NGUARD !Nozzle thickness (in grid cells)
    
    xn = sim_smallX
    xn(HE4_SPEC) = 1.0
  
    call gr_mpoleCenterOfMass(DENS_VAR)
    call wd_interp(don_mass, don_pos, acc_pos, noz_pos, don_dist, mrate)

    if (dr_simTime .lt. sim_tInitial + sim_tRelax) then
        call PhysicalConstants_get("Newton", G)
        !adj_cfactor = sqrt(1.0-(dt*(1.0-sim_relaxRate)))
        relax_rate = (dr_simTime - sim_tInitial)/sim_tRelax*(1.0 - sim_relaxRate) + sim_relaxRate 
        !adj_cfactor = exp(-dt/sqrt(wd_radius**3./(2.*G*sim_accMass)))
        do lb = 1, blockCount
            call Grid_getBlkIndexLimits(blockList(lb),blkLimits,blkLimitsGC)
            sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
            allocate(xCoord(sizeX),stat=istat)
            sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
            allocate(yCoord(sizeY),stat=istat)
            sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
            allocate(zCoord(sizeZ),stat=istat)
  
            if (NDIM == 3) call Grid_getCellCoords&
                                (KAXIS, blockList(lb), CENTER, gcell, zCoord, sizeZ)
            if (NDIM >= 2) call Grid_getCellCoords&
                                (JAXIS, blockList(lb), CENTER,gcell, yCoord, sizeY)
            call Grid_getCellCoords(IAXIS, blockList(lb), CENTER, gcell, xCoord, sizeX)
            call Grid_getBlkPtr(blockList(lb),solnData)
            allocate(velCopy(VELX_VAR:VELZ_VAR,&
                             blkLimits(LOW, IAXIS):blkLimits(HIGH, IAXIS), &
                             blkLimits(LOW, JAXIS):blkLimits(HIGH, JAXIS), &
                             blkLimits(LOW, KAXIS):blkLimits(HIGH, KAXIS)),stat=istat)
            velCopy(VELX_VAR:VELZ_VAR,:,:,:) = solnData(VELX_VAR:VELZ_VAR,blkLimits(LOW, IAXIS):blkLimits(HIGH, IAXIS),&
                                                                          blkLimits(LOW, JAXIS):blkLimits(HIGH, JAXIS),&
                                                                          blkLimits(LOW, KAXIS):blkLimits(HIGH, KAXIS))
            do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
                do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
                    do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
                        !Arnett 1994
                        !velCopy(VELX_VAR,i,j,k) = 0.5*solnData(VELX_VAR,i,j,k) + 0.25*solnData(VELX_VAR,i-1,j,k) &
                        !    + 0.25*solnData(VELX_VAR,i+1,j,k)
                        !velCopy(VELY_VAR,i,j,k) = 0.5*solnData(VELY_VAR,i,j,k) + 0.25*solnData(VELY_VAR,i,j-1,k) &
                        !    + 0.25*solnData(VELY_VAR,i,j+1,k)
                        !velCopy(VELZ_VAR,i,j,k) = 0.5*solnData(VELZ_VAR,i,j,k) + 0.25*solnData(VELZ_VAR,i,j,k-1) &
                        !    + 0.25*solnData(VELZ_VAR,i,j,k+1)

                        !reduce by constant factor
                        !velCopy(VELX_VAR,i,j,k) = velCopy(VELX_VAR,i,j,k)*sim_relaxRate
                        !velCopy(VELY_VAR,i,j,k) = velCopy(VELY_VAR,i,j,k)*sim_relaxRate
                        !velCopy(VELZ_VAR,i,j,k) = velCopy(VELZ_VAR,i,j,k)*sim_relaxRate
                        velCopy(VELX_VAR,i,j,k) = velCopy(VELX_VAR,i,j,k)*relax_rate
                        velCopy(VELY_VAR,i,j,k) = velCopy(VELY_VAR,i,j,k)*relax_rate
                        velCopy(VELZ_VAR,i,j,k) = velCopy(VELZ_VAR,i,j,k)*relax_rate
                        
                        !!Arnett 1994
                        !if (solnData(DENS_VAR,i,j,k) .gt. donor_dens) then
                        !    velCopy(VELX_VAR,i,j,k) = 0.5*solnData(VELX_VAR,i,j,k)
                        !    if (solnData(DENS_VAR,i-1,j,k) .gt. donor_dens) &
                        !        velCopy(VELX_VAR,i,j,k) = velCopy(VELX_VAR,i,j,k) + 0.25*solnData(VELX_VAR,i-1,j,k)
                        !    if (solnData(DENS_VAR,i+1,j,k) .gt. donor_dens) &
                        !        velCopy(VELX_VAR,i,j,k) = velCopy(VELX_VAR,i,j,k) + 0.25*solnData(VELX_VAR,i+1,j,k)
                        !    velCopy(VELY_VAR,i,j,k) = 0.5*solnData(VELY_VAR,i,j,k)
                        !    if (solnData(DENS_VAR,i,j-1,k) .gt. donor_dens) &
                        !        velCopy(VELY_VAR,i,j,k) = velCopy(VELY_VAR,i,j,k) + 0.25*solnData(VELY_VAR,i,j-1,k)
                        !    if (solnData(DENS_VAR,i,j+1,k) .gt. donor_dens) &
                        !        velCopy(VELY_VAR,i,j,k) = velCopy(VELY_VAR,i,j,k) + 0.25*solnData(VELY_VAR,i,j+1,k)
                        !    velCopy(VELZ_VAR,i,j,k) = 0.5*solnData(VELZ_VAR,i,j,k)
                        !    if (solnData(DENS_VAR,i,j,k-1) .gt. donor_dens) &
                        !        velCopy(VELZ_VAR,i,j,k) = velCopy(VELZ_VAR,i,j,k) + 0.25*solnData(VELZ_VAR,i,j,k-1)
                        !    if (solnData(DENS_VAR,i,j,k+1) .gt. donor_dens) &
                        !        velCopy(VELZ_VAR,i,j,k) = velCopy(VELZ_VAR,i,j,k) + 0.25*solnData(VELZ_VAR,i,j,k+1)

                        if (solnData(DENS_VAR,i,j,k) .gt. rhop(ipos)) then
                            solnData(TEMP_VAR,i,j,k) = sim_accTemp
                        endif
                        eosData(EOS_DENS) = solnData(DENS_VAR,i,j,k)
                        eosData(EOS_TEMP) = solnData(TEMP_VAR,i,j,k)
                        call Eos(MODE_DENS_TEMP,1,eosData,solnData(SPECIES_BEGIN:SPECIES_END,i,j,k))
                        solnData(EINT_VAR,i,j,k) = eosData(EOS_EINT)
                        solnData(PRES_VAR,i,j,k) = eosData(EOS_PRES)
                        solnData(GAMC_VAR,i,j,k) = eosData(EOS_GAMC)
                        solnData(ENER_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) + &
                            0.5*(velCopy(VELX_VAR,i,j,k)**2. + velCopy(VELY_VAR,i,j,k)**2. + velCopy(VELZ_VAR,i,j,k)**2.)
                        !endif
                    enddo
                enddo
            enddo
  
            solnData(VELX_VAR:VELZ_VAR,blkLimits(LOW, IAXIS):blkLimits(HIGH, IAXIS),&
                                       blkLimits(LOW, JAXIS):blkLimits(HIGH, JAXIS),&
                                       blkLimits(LOW, KAXIS):blkLimits(HIGH, KAXIS)) = velCopy(VELX_VAR:VELZ_VAR,:,:,:)
            deallocate(velCopy)
  
            do k = blkLimitsGC(LOW, KAXIS), blkLimitsGC(HIGH, KAXIS)
                zz = zCoord(k)
                do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)
                    yy = yCoord(j)
                    do i = blkLimitsGC(LOW, IAXIS), blkLimitsGC(HIGH, IAXIS)
                        xx = xCoord(i)
                        !Reset region near donor sink particle
                        if (sqrt((xx-don_pos(1))**2.+(yy-don_pos(2))**2.+(zz-don_pos(3))**2.) .lt. donor_size) then
                            solnData(DENS_VAR,i,j,k) = donor_dens
                            solnData(TEMP_VAR,i,j,k) = sim_donorTemp
                            solnData(SPECIES_BEGIN:SPECIES_END,i,j,k) = xn
                            solnData(VELX_VAR,i,j,k) = 0.0d0
                            solnData(VELY_VAR,i,j,k) = 0.0d0
                            solnData(VELZ_VAR,i,j,k) = 0.0d0
                            eosData(EOS_DENS) = solnData(DENS_VAR,i,j,k)
                            eosData(EOS_TEMP) = solnData(TEMP_VAR,i,j,k)
                            call Eos(MODE_DENS_TEMP,1,eosData,solnData(SPECIES_BEGIN:SPECIES_END,i,j,k))
                            solnData(EINT_VAR,i,j,k) = eosData(EOS_EINT)
                            solnData(PRES_VAR,i,j,k) = eosData(EOS_PRES)
                            solnData(GAMC_VAR,i,j,k) = eosData(EOS_GAMC)
                            solnData(ENER_VAR,i,j,k) = solnData(EINT_VAR,i,j,k)
                        endif
                    enddo
                enddo
            enddo
            call Grid_releaseBlkPtr(blockList(lb), solnData)
            deallocate(xCoord)
            deallocate(yCoord)
            deallocate(zCoord)
            call Eos_wrapped(MODE_DENS_TEMP,blkLimitsGC,blockList(lb))
        enddo
    else
        call PhysicalConstants_get('Newton',G)
        call RuntimeParameters_get('sim_xctr',xcenter)
        call RuntimeParameters_get('sim_yctr',ycenter)
        call RuntimeParameters_get('sim_zctr',zcenter)
        call Grid_getMinCellSize(min_cell_size)
  
        ! calculate the velocity of the center of the star so it can be substracted from all grid cells.
        avg_vel = 0
        cnt = 0
        do lb = 1, blockCount
            call Grid_getBlkIndexLimits(blockList(lb),blkLimits,blkLimitsGC)
            sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
            allocate(xCoord(sizeX),stat=istat)
            sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
            allocate(yCoord(sizeY),stat=istat)
            sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
            allocate(zCoord(sizeZ),stat=istat)
  
            if (NDIM == 3) call Grid_getCellCoords&
                                (KAXIS, blockList(lb), CENTER, gcell, zCoord, sizeZ)
            if (NDIM >= 2) call Grid_getCellCoords&
                                (JAXIS, blockList(lb), CENTER, gcell, yCoord, sizeY)
            call Grid_getCellCoords(IAXIS, blockList(lb), CENTER, gcell, xCoord, sizeX)
  
            call Grid_getBlkPtr(blockList(lb),solnData)
            do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
                zz = zCoord(k)
                do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
                    yy = yCoord(j)
                    do i = blkLimits(LOW,IAXIS), blkLimits(HIGH, IAXIS)
                        xx = xCoord(i)
                        if (sqrt((xx-Xcm)**2.+(yy-Ycm)**2.+(zz-Zcm)**2.) .le. 5*min_cell_size) then
                            avg_vel(1) = avg_vel(1) + solnData(VELX_VAR,i,j,k)
                            avg_vel(2) = avg_vel(2) + solnData(VELY_VAR,i,j,k)
                            avg_vel(3) = avg_vel(3) + solnData(VELZ_VAR,i,j,k)
                            cnt = cnt + 1
                        endif
                    enddo
                enddo
            enddo
            call Grid_releaseBlkPtr(blockList(lb), solnData)
            deallocate(xCoord)
            deallocate(yCoord)
            deallocate(zCoord)
        enddo
  
        call MPI_ALLREDUCE(cnt, tot_cnt, 1, FLASH_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_ALLREDUCE(avg_vel, tot_avg_vel, MDIM, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
        if (tot_cnt .eq. 0) call Driver_abortFlash("Zero cells in averaging volume!")
        tot_avg_vel = tot_avg_vel / tot_cnt
  
        !how precise do we want the mass?
        massacc = 1.d-5
  
        !the following makes the profile shape independent of time-step size
        !not the best way, this is resolution dependent!
        falling_time = sqrt(2*((Xcm - noz_pos(1))**2. + (Ycm - noz_pos(2))**2. + (Zcm - noz_pos(3))**2.)*noz_thick/2.*min_cell_size/G/Mtot)
        !print *, 'falling_time, Xcm, Ycm, Zcm, noz_pos(1,2,3), Mtot, mrate', falling_time, Xcm, Ycm, Zcm, &
        !    noz_pos(1), noz_pos(2), noz_pos(3), Mtot, mrate
        !include donor when calculating falling time
        !falling_time = sqrt(2*min_cell_size/G/&
        !    abs(Mtot/(sqrt((Xcm-noz_pos(1))**2.+(Ycm-noz_pos(2))**2.+(Zcm-noz_pos(3))**2.)-min_cell_size)**2.-&
        !    don_mass/(sqrt((don_pos(1)-noz_pos(1))**2.+(don_pos(2)-noz_pos(2))**2.+(don_pos(3)-noz_pos(3))**2.)+min_cell_size)**2.))
        mrate_adj = mrate * falling_time
  
        !guess central density (cgs) 
        call nozz_guess_rho_0(mrate_adj,rho0_cgs)
        !print *, 'rho0_cgs, mrate_adj', rho0_cgs, mrate_adj
  
        !loop until mass is correct
        call nozz_intout(xn,rho0_cgs,sim_donorTemp,m_nozz_tot,don_mass,cnt,r,rho,noz_pos(1),noz_pos(2),noz_pos(3),don_pos(1),don_pos(2),don_pos(3))
        i = 0
        do while (abs(m_nozz_tot-mrate_adj)/mrate_adj.gt.massacc)
            rho0_cgs = rho0_cgs * 10**((log10(mrate_adj)-log10(m_nozz_tot))/2.)
            !print *, 'rho0_cgs', rho0_cgs
            call nozz_intout(xn,rho0_cgs,sim_donorTemp,m_nozz_tot,don_mass,cnt,r,rho,noz_pos(1),noz_pos(2),noz_pos(3),don_pos(1),don_pos(2),don_pos(3))
            if (i .gt. 100) then
                print *, 'nozzle mass did not converge!'
                stop
            endif
            i = i + 1
        enddo
  
        !adjust rho so that mass is added proportional to time-step size
        rho = rho*dt/falling_time
  
        !radius of nozzle
        rnozz = r(cnt)
        delt = rnozz/dble(cnt)
  
        do lb = 1, blockCount
            call Grid_getBlkIndexLimits(blockList(lb),blkLimits,blkLimitsGC)
            sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
            allocate(xCoord(sizeX),stat=istat)
            sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
            allocate(yCoord(sizeY),stat=istat)
            sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
            allocate(zCoord(sizeZ),stat=istat)
  
            if (NDIM == 3) call Grid_getCellCoords&
                                (KAXIS, blockList(lb), CENTER, gcell, zCoord, sizeZ)
            if (NDIM >= 2) call Grid_getCellCoords&
                                (JAXIS, blockList(lb), CENTER,gcell, yCoord, sizeY)
            call Grid_getCellCoords(IAXIS, blockList(lb), CENTER, gcell, xCoord, sizeX)
            local_cell_size = abs(xCoord(2) - xCoord(1))
            c2 = local_cell_size/2.
  
            !if ((noz_pos(1) .lt. xCoord(1) - c2) .or. (noz_pos(1) .gt. xCoord(sizeX) + c2)) cycle
            !if ((yCoord(sizeY) - c2 .lt. noz_pos(2) - rnozz) .or. (yCoord(1) + c2 .gt. noz_pos(2) + rnozz)) cycle
            !if ((zCoord(sizeZ) - c2 .lt. noz_pos(3) - rnozz) .or. (zCoord(1) + c2 .gt. noz_pos(3) + rnozz)) cycle
  
            changedGrid = .false.
            call Grid_getBlkPtr(blockList(lb),solnData)
            do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
               zz = zCoord(k)
               do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)
                  yy = yCoord(j)
                  theta = atan2(yy-noz_pos(2),zz-noz_pos(3))
                  yzdist = sqrt((yy-noz_pos(2))**2 + (zz-noz_pos(3))**2)
                  !Have to figure out WHY this is true!
                  local_xnozz = noz_pos(1) - yzdist/tan(1.0)
  
                  do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH, IAXIS)
                     xx = xCoord(i)
  
                     !Subtract the COM's velocity from all grid cells (keeps WD in center of grid)
                     solnData(VELX_VAR,i,j,k) = solnData(VELX_VAR,i,j,k) - tot_avg_vel(1)
                     solnData(VELY_VAR,i,j,k) = solnData(VELY_VAR,i,j,k) - tot_avg_vel(2)
                     solnData(VELZ_VAR,i,j,k) = solnData(VELZ_VAR,i,j,k) - tot_avg_vel(3)

                     !Add a fictitious force to move the center of mass to the center of the simulation
                     solnData(VELX_VAR,i,j,k) = solnData(VELX_VAR,i,j,k) + (xcenter - Xcm)/dt/10.
                     solnData(VELY_VAR,i,j,k) = solnData(VELY_VAR,i,j,k) + (ycenter - Ycm)/dt/10
                     solnData(VELZ_VAR,i,j,k) = solnData(VELZ_VAR,i,j,k) + (zcenter - Zcm)/dt/10.
  
                     if (yzdist .le. rnozz .and. .not. ((zz .lt. noz_pos(3) - rnozz) .or. (zz .gt. noz_pos(3) + rnozz) .or. &
                         (yy .lt. noz_pos(2) - rnozz) .or. (yy .gt. noz_pos(2) + rnozz) .or. &
                         (xx .lt. local_xnozz - (noz_thick/2.+0.5)*min_cell_size) .or. &
                         (xx .gt. local_xnozz + (noz_thick/2.+0.5)*min_cell_size))) then
                     !if (yzdist .le. rnozz .and. .not. ((zz .lt. noz_pos(3) - rnozz) .or. (zz .gt. noz_pos(3) + rnozz) .or. &
                     !    (yy .lt. noz_pos(2) - rnozz) .or. (yy .gt. noz_pos(2) + rnozz) .or. (xx .lt. local_xnozz - c2) .or. &
                     !    (xx .gt. local_xnozz + c2))) then
                         changedGrid = .true.
                         ri = int(floor(yzdist/delt))
                         dr = yzdist/delt - ri
                         if (ri .ge. cnt - 1) then
                             new_mass = rho(ri) + (sim_rhoAmbient-rho(ri))*dr
                         else
                             new_mass = rho(ri) + (rho(ri+1)-rho(ri))*dr
                         endif
  
                         !We want to add mass to the nozzle without destroying the composition of the stuff that's
                         !already there...
                         !if (xx .lt. local_xnozz + c2 + (noz_thick - 1.0)*min_cell_size) solnData(DENS_VAR,i,j,k) = 0.0
                         mass_sum = solnData(DENS_VAR,i,j,k)+new_mass
                         solnData(SPECIES_BEGIN:SPECIES_END,i,j,k) = &
                             (solnData(DENS_VAR,i,j,k)*solnData(SPECIES_BEGIN:SPECIES_END,i,j,k)+new_mass*xn)/mass_sum
                         solnData(SPECIES_BEGIN:SPECIES_END,i,j,k) = solnData(SPECIES_BEGIN:SPECIES_END,i,j,k)/&
                             sum(solnData(SPECIES_BEGIN:SPECIES_END,i,j,k))
                         !do l = 1, NSPECIES
                         !    solnData(SPECIES_BEGIN+l-1,i,j,k) = min(max(solnData(SPECIES_BEGIN+l-1,i,j,k),sim_smallX),1.0)
                         !enddo
                         if (xx .lt. local_xnozz) then
                             solnData(VELX_VAR,i,j,k) = solnData(DENS_VAR,i,j,k)*solnData(VELX_VAR,i,j,k)/mass_sum
                             solnData(VELY_VAR,i,j,k) = solnData(DENS_VAR,i,j,k)*solnData(VELY_VAR,i,j,k)/mass_sum
                             solnData(VELZ_VAR,i,j,k) = solnData(DENS_VAR,i,j,k)*solnData(VELZ_VAR,i,j,k)/mass_sum
                         else
                             accel_frac = (xx-local_xnozz+c2)/min_cell_size
                             new_vel(1) = G*Mtot/(local_xnozz+min_cell_size-Xcm)**2.*falling_time*accel_frac
                             new_vel(2) = G*yzdist*(Mtot/((local_xnozz+min_cell_size-Xcm)**2.+yzdist**2.)**1.5 + &
                                 don_mass/((local_xnozz-don_pos(1)-min_cell_size)**2.+yzdist**2.)**1.5)*falling_time*accel_frac
                             new_vel(3) = -new_vel(2)*cos(theta)
                             new_vel(2) = -new_vel(2)*sin(theta)
                             !solnData(VELX_VAR:VELZ_VAR,i,j,k) = new_vel(1:3)
                             solnData(VELX_VAR,i,j,k) = (solnData(DENS_VAR,i,j,k)*solnData(VELX_VAR,i,j,k) + new_mass*new_vel(1))/mass_sum
                             solnData(VELY_VAR,i,j,k) = (solnData(DENS_VAR,i,j,k)*solnData(VELY_VAR,i,j,k) + new_mass*new_vel(2))/mass_sum
                             solnData(VELZ_VAR,i,j,k) = (solnData(DENS_VAR,i,j,k)*solnData(VELZ_VAR,i,j,k) + new_mass*new_vel(3))/mass_sum
                         endif
  
                         if ((abs(min_cell_size/solnData(VELX_VAR,i,j,k)) .lt. 1e-4) .or. &
                             (abs(min_cell_size/solnData(VELY_VAR,i,j,k)) .lt. 1e-4) .or. &
                             (abs(min_cell_size/solnData(VELZ_VAR,i,j,k)) .lt. 1e-4)) then
                             write(*,*), 'min_cell_size', min_cell_size
                             write(*,*), 'new_vel', new_vel(1), new_vel(2), new_vel(3)
                             write(*,*), 'xx, c2', xx, c2
                             write(*,*), 'dens, vel', solnData(DENS_VAR,i,j,k), solnData(VELX_VAR,i,j,k), solnData(VELY_VAR,i,j,k), solnData(VELZ_VAR,i,j,k)
                             write(*,*), 'G, Mtot, don_mass, local_xnozz, accel_frac', G, Mtot, don_mass, local_xnozz, accel_frac
                             write(*,*), 'falling_time, don_pos(1), Xcm, yzdist', falling_time, don_pos(1), Xcm, yzdist
                             write(*,*), 'new_mass, mass_sum, theta', new_mass, mass_sum, theta
                             write(*,*), 'species_data', solnData(SPECIES_BEGIN:SPECIES_END,i,j,k)
                             call Driver_abortFlash("Nozzle velocities too large!")
                         endif
                             
                         solnData(DENS_VAR,i,j,k) = mass_sum
                         solnData(TEMP_VAR,i,j,k) = sim_donorTemp
                         eosData(EOS_DENS) = solnData(DENS_VAR,i,j,k)
                         eosData(EOS_TEMP) = solnData(TEMP_VAR,i,j,k)
                         call Eos(MODE_DENS_TEMP,1,eosData,solnData(SPECIES_BEGIN:SPECIES_END,i,j,k))
                         solnData(EINT_VAR,i,j,k) = eosData(EOS_EINT)
                         solnData(TEMP_VAR,i,j,k) = eosData(EOS_TEMP)
                         !solnData(PRES_VAR,i,j,k) = eosData(EOS_PRES)
                     endif
  
                     !Reset region near donor sink particle
                     if (sqrt((xx-don_pos(1))**2.+(yy-don_pos(2))**2.+(zz-don_pos(3))**2.) .lt. donor_size) then
                         changedGrid = .true.
                         solnData(DENS_VAR,i,j,k) = donor_dens
                         solnData(TEMP_VAR,i,j,k) = sim_donorTemp
                         solnData(SPECIES_BEGIN:SPECIES_END,i,j,k) = xn
                         eosData(EOS_DENS) = solnData(DENS_VAR,i,j,k)
                         eosData(EOS_TEMP) = solnData(TEMP_VAR,i,j,k)
                         call Eos(MODE_DENS_TEMP,1,eosData,xn)
                         solnData(TEMP_VAR,i,j,k) = eosData(EOS_TEMP)
                         solnData(EINT_VAR,i,j,k) = eosData(EOS_EINT)
                         !solnData(PRES_VAR,i,j,k) = eosData(EOS_PRES)
                         solnData(VELX_VAR,i,j,k) = 0.0d0
                         solnData(VELY_VAR,i,j,k) = 0.0d0
                         solnData(VELZ_VAR,i,j,k) = 0.0d0
                     endif
  
                     solnData(ENER_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) + &
                         0.5*(solnData(VELX_VAR,i,j,k)**2. + solnData(VELY_VAR,i,j,k)**2. + solnData(VELZ_VAR,i,j,k)**2.)
                  enddo
               enddo
            enddo
            if (changedGrid) call Eos_wrapped(MODE_DENS_TEMP,blkLimitsGC,blockList(lb))
            call Grid_releaseBlkPtr(blockList(lb), solnData)
            deallocate(xCoord)
            deallocate(yCoord)
            deallocate(zCoord)
        enddo
  
        call Stir(blockCount, blockList, dt) 
        !call Flame_step(blockCount, blockList, dt)
        call Burn(blockCount, blockList, dt) 
        call Heat(blockCount, blockList, dt, dr_simTime) 
        call Cool(blockCount, blockList, dt, dr_simTime)
    endif
  
  
    return
end subroutine Driver_sourceTerms
