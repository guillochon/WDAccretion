!!****if* source/physics/Gravity/GravityMain/Poisson/Gravity_accelOneRow
!!
!! NAME
!!
!!  Gravity_accelOneRow
!!
!!
!! SYNOPSIS
!!
!!  Gravity_accelOneRow(integer(2),intent(in):: pos, 
!!                      integer, intent(in) :: sweepDir, 
!!                      integer, intent(in) :: blockID, 
!!                      integer, intent(in) :: numCells, 
!!                      real(numCells),intent(out) :: grav, 
!!                      integer, intent(in),optional :: potentialIndex)
!!                      
!!                      
!!
!! DESCRIPTION
!!
!!  Compute components of the zone-averaged gravitational
!!  acceleration on all mesh blocks.  Either a single component
!!  of the acceleration or all three can be computed.
!!
!!  This routine computes the gravitational acceleration for a row
!!  of zones in a specified direction in a given block. First-order
!!  finite-volume differencing is used everywhere.  Formulae based
!!  on long stencils (usually high-order) may produce differences
!!  at the block boundaries for siblings as hydro solver may require
!!  several valid guard cells (e.g., PPM with parabolic
!!  interpolation for force terms needs 3 valid guard cells). Not
!!  providing such valid data may result in violation of conservation. 
!!
!! ARGUMENTS
!!
!!  pos     -       Row indices transverse to the sweep direction
!!  sweepDir   -       The sweep direction:  test against sweep_x,
!!                                 sweep_y, and sweep_z
!!  blockID   -     The local identifier of the block to work on
!!  grav()   -       Array to receive result
!!  numCells -       Number of cells to update in grav array
!!  potentialIndex      -  if specified,  Variable # to take as potential.
!!                         Default is GPOT_VAR for the potential stored in the
!!                         gpot slot of unk, which should correspond to the
!!                         potential at the current timestep.
!!
!!
!!***


subroutine Gravity_accelOneRow (pos, sweepDir, blockID, numCells, grav, potentialIndex)

  use Driver_data, ONLY: dr_simTime
  use Grid_interface, ONLY : Grid_getBlkPhysicalSize, Grid_getBlkPtr, &
    Grid_releaseBlkPtr, Grid_getCellCoords, Grid_getBlkIndexLimits
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use gr_mpoleData, ONLY: Mtot, Xcm, Ycm, Zcm
  use Simulation_data, ONLY: sim_table, sim_tableRows, sim_tableCols, &
      sim_tNozz, sim_tRelax, sim_donorMass, sim_tInitial
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Gravity_data, ONLY: grv_thresh
  implicit none

#include "constants.h"
#include "Flash.h"

  integer, dimension(2), intent(in) :: pos
  integer, intent(in)               :: sweepDir, blockID,  numCells
  real, intent(inout)               :: grav(numCells)
  integer, intent(in),optional        :: potentialIndex
  real            :: blockSize(MDIM)
  real, pointer, dimension(:,:,:,:) :: solnVec
  
  integer         :: ii, iimin, iimax, lb
  real            :: gpot(numCells), velx(numCells), vely(numCells), dens(numCells), delxinv
  real, parameter :: onesixth = 1.e0/6.e0
  real            :: don_pos(3), acc_pos(3), noz_pos(3)
  integer         :: potVar
  
#ifdef FIXEDBLOCKSIZE
  real,dimension(GRID_KHI_GC) :: zCenter
  real,dimension(GRID_JHI_GC) :: yCenter
  real,dimension(GRID_IHI_GC) :: xCenter
#else
  real,allocatable,dimension(:) ::xCenter,yCenter,zCenter
#endif
  real :: dr32, tmpdr32

  integer :: sizeX,sizeY,sizez

  integer :: j,k
  logical :: gcell = .true., useCoriolis = .false.
  real :: t, newton
  real :: don_mass, don_dist, acc_rate
  real :: msum, m1, m2, d1, d12, theta, r, omega, ind_frac
  !==================================================
  
  

  call Grid_getBlkPhysicalSize(blockID, blockSize)
  
  
  call Grid_getBlkPtr(blockID, solnVec)

!! IF a variable index is explicitly specified, assume that as the potential
!! otherwise use the default current potential GPOT_VAR  
  if(present(potentialIndex)) then
     potVar=potentialIndex
  else
     potVar=GPOT_VAR
  end if

  iimin   = 1
  iimax   = numCells


  !Get row of potential values and compute inverse of zone spacing  
  if (sweepDir == SWEEP_X) then                     ! x-direction
     delxinv = real(NXB) / blockSize(IAXIS)
     
     gpot(:) = solnVec(potVar,:,pos(1),pos(2))
     dens(:) = solnVec(DENS_VAR,:,pos(1),pos(2))
     velx(:) = solnVec(VELX_VAR,:,pos(1),pos(2))
     vely(:) = solnVec(VELY_VAR,:,pos(1),pos(2))
  elseif (sweepDir == SWEEP_Y) then                 ! y-direction
     delxinv = real(NYB) / blockSize(JAXIS)
     
     gpot(:) = solnVec(potVar,pos(1),:,pos(2))
     dens(:) = solnVec(DENS_VAR,pos(1),:,pos(2))
     velx(:) = solnVec(VELX_VAR,pos(1),:,pos(2))
     vely(:) = solnVec(VELY_VAR,pos(1),:,pos(2))
  else                                            ! z-direction
     delxinv = real(NZB) / blockSize(KAXIS)
     
     gpot(:) = solnVec(potVar,pos(1),pos(2),:)
     dens(:) = solnVec(DENS_VAR,pos(1),pos(2),:)
     velx(:) = solnVec(VELX_VAR,pos(1),pos(2),:)
     vely(:) = solnVec(VELY_VAR,pos(1),pos(2),:)
  endif
  
  !-------------------------------------------------------------------------------
  
  !               Compute gravitational acceleration
  
  
  !**************** first-order differences
  !                 preserves conservation
  
  delxinv = 0.5e0 * delxinv
  
  do ii = iimin+1, iimax-1
     grav(ii) = delxinv * (gpot(ii-1) - gpot(ii+1))
  enddo
  
  grav(iimin) = grav(iimin+1)     ! this is invalid data - must not be used
  grav(iimax) = grav(iimax-1)
  
  
  call Grid_releaseBlkPtr(blockID, solnVec)
  
  !Now include the point mass
  j=pos(1)
  k=pos(2)
#ifndef FIXEDBLOCKSIZE
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  sizeX=blkLimitsGC(HIGH,IAXIS)
  sizeY=blkLimitsGC(HIGH,JAXIS)
  sizeZ=blkLimitsGC(HIGH,KAXIS)
  allocate(xCenter(sizeX))
  allocate(yCenter(sizeY))
  allocate(zCenter(sizeZ))
#else
  sizeX=GRID_IHI_GC
  sizeY=GRID_JHI_GC
  sizeZ=GRID_KHI_GC
#endif
  call PhysicalConstants_get("Newton", newton)
  call wd_interp(don_mass, don_pos, acc_pos, noz_pos, don_dist, acc_rate)

  if (dr_simTime .ge. sim_tInitial + sim_tRelax) then
      useCoriolis = .true.
  else
      useCoriolis = .false.
  endif

  m1 = don_mass
  m2 = Mtot
  d12 = sqrt((Xcm - don_pos(1))**2. + (Ycm - don_pos(2))**2. + (Zcm - don_pos(3))**2.)
  msum = m1 + m2
  d1 = m2 / msum * d12 !center of mass coordinate of binary along line between centers of masses of the stars
  d1 = d1*cos(atan(abs(Zcm - don_pos(3))/sqrt((Xcm - don_pos(1))**2.+(Ycm - don_pos(2))**2.))) !distance of barycenter from donor in XY plane
  omega = 1. / sqrt(d12**3. / newton / msum)
  
  zCenter = 0.
  yCenter = 0.
  if (NDIM == 3) then 
     call Grid_getCellCoords(KAXIS, blockID, CENTER, gcell, zCenter, sizeZ)
     zCenter = zCenter - don_pos(3)
  endif
  if (NDIM >= 2) then
     call Grid_getCellCoords(JAXIS, blockID, CENTER, gcell, yCenter, sizeY)
     yCenter = yCenter - don_pos(2)
  endif
  call Grid_getCellCoords(IAXIS, blockID, CENTER, gcell, xCenter, sizeX)
  xCenter = xCenter - don_pos(1)
  

  if (sweepDir .eq. SWEEP_X) then                       ! x-component

     tmpdr32 = yCenter(j)*yCenter(j) + zCenter(k)*zCenter(k) 

     do ii = 1, numCells

        dr32 = sqrt(xCenter(ii)*xCenter(ii) + tmpdr32)
        dr32 = dr32*dr32*dr32

        theta = atan(yCenter(ii)/(xCenter(ii) - d1))
        r = sqrt((xCenter(ii) - d1)**2. + yCenter(ii)**2.)
        grav(ii) = grav(ii) - newton*m1*xCenter(ii)/dr32 + omega**2. * r * cos(theta)
        if (useCoriolis) grav(ii) = grav(ii) + 2.*omega*vely(ii)
     enddo


  else if (sweepDir .eq. SWEEP_Y) then          ! y-component

     tmpdr32 = xCenter(j)*xCenter(j) + zCenter(k)*zCenter(k) 

     do ii = 1, numCells
        
        dr32 = sqrt(yCenter(ii)*yCenter(ii) + tmpdr32)
        dr32 = dr32*dr32*dr32

        theta = atan(yCenter(ii)/(xCenter(ii) - d1))
        r = sqrt((xCenter(ii) - d1)**2. + yCenter(ii)**2.)
        grav(ii) = grav(ii) - newton*m1*yCenter(ii)/dr32 + omega**2. * r * sin(theta)
        if (useCoriolis) grav(ii) = grav(ii) - 2.*omega*velx(ii)
     enddo

  else if (sweepDir .eq. SWEEP_Z) then          ! z-component

     tmpdr32 = xCenter(j)*xCenter(j) + yCenter(k)*yCenter(k) 

     do ii = 1, numCells
        
        dr32 = sqrt(zCenter(ii)*zCenter(ii) + tmpdr32)           
        dr32 = dr32*dr32*dr32
        
        grav(ii) = grav(ii) - newton*m1*zCenter(ii)/dr32
     enddo

  endif

!==============================================================================
#ifndef FIXEDBLOCKSIZE
  deallocate(xCenter)
  deallocate(yCenter)
  deallocate(zCenter)
#endif
  do ii = 1, numCells                             
      if (dens(ii) .le. grv_thresh) grav(ii) = 0.0
  enddo                                           

  return
   
end subroutine Gravity_accelOneRow


