!!****if* source/Simulation/SimulationMain/Sedov/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!! 
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer :: blockId, 
!!                       integer :: myPE)
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.  This version sets up the Sedov spherical
!!  explosion problem.
!!
!!  References:  Sedov, L. I., 1959, Similarity and Dimensional Methods
!!                 in Mechanics (New York:  Academic)
!!
!!               Landau, L. D., & Lifshitz, E. M., 1987, Fluid Mechanics,
!!                 2d ed. (Oxford:  Pergamon)
!!
!! ARGUMENTS
!!
!!  blockId -        The number of the block to initialize
!!  myPE -          current processor number
!!
!! PARAMETERS
!!
!!  sim_pAmbient       Initial ambient pressure
!!  sim_rhoAmbient     Initial ambient density
!!  sim_expEnergy      Explosion energy (distributed over 2^dimen central zones)
!!  sim_rInit          Radial position of inner edge of grid (for 1D )
!!  sim_xctr           Explosion center coordinates
!!  sim_yctr           Explosion center coordinates
!!  sim_zctr           Explosion center coordinates
!!  sim_nsubzones      Number of `sub-zones` in cells for applying 1d profile
!!
!!
!!***

subroutine Simulation_initBlock (blockId, myPE)

  use Simulation_data, ONLY: sim_xMax, sim_xMin, sim_yMax, sim_yMin, sim_zMax, sim_zMin, &
     &  sim_tInitial, sim_gamma, sim_pAmbient, sim_rhoAmbient, &
     &  sim_smallX, sim_smallRho, sim_smallP, &
     &  sim_nSubZones, sim_xCenter, sim_yCenter, sim_zCenter, sim_inSubzm1, sim_inszd, &
     &  radius, rhop, eint, ipos
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr,&
    Grid_getCellCoords, Grid_putPointData
  use PhysicalConstants_interface, ONLY: PhysicalConstants_get
  use Multispecies_interface, ONLY:  Multispecies_getSumFrac, Multispecies_getSumInv, Multispecies_getAvg
  
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"
  
  integer,intent(IN) ::  blockId
  integer,intent(IN) ::  myPE
  
  integer  ::  i, j, k, n, jLo, jHi
  integer  ::  ii, jj, kk, put
  real     ::  distInv, xDist, yDist, zDist
  real     ::  sumRho, sumE
  real     ::  vel, diagonal
  real     ::  xx, dxx, yy, dyy, zz, dzz, frac, theta, phi
  real     ::  vx, vy, vz, p, rho, e, ek, t, mp, kb, G, vtot
  real     ::  dist, gam
  logical  ::  validGeom
  integer  ::  istat

  real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  integer :: sizeX,sizeY,sizeZ
  integer,dimension(MDIM) :: axis
  real, dimension(SPECIES_BEGIN:SPECIES_END) :: xn = 0.0
  real, dimension(:,:,:,:),pointer :: solnData

  logical :: gcell = .true.

  double precision  mtot,mu
  integer mode
  !
  !  Construct the radial samples needed for the initialization.
  !
  
  call PhysicalConstants_get("proton mass", mp)
  call PhysicalConstants_get("Boltzmann", kb)
  call PhysicalConstants_get("Newton", G)
  
  !xn(HE4_SPEC) = 1.0
  xn(C12_SPEC) = 0.5
  xn(O16_SPEC) = 0.5
  
  ! get the coordinate information for the current block from the database

  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
  allocate(xCoord(sizeX),stat=istat)
  sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
  allocate(yCoord(sizeY),stat=istat)
  sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
  allocate(zCoord(sizeZ),stat=istat)

  if (NDIM == 3) call Grid_getCellCoords&
                      (KAXIS, blockId, CENTER, gcell, zCoord, sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords&
                      (JAXIS, blockId, CENTER,gcell, yCoord, sizeY)
  call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, xCoord, sizeX)
  !
  !     For each cell
  !  
#ifdef FL_NON_PERMANENT_GUARDCELLS
  call Grid_getBlkPtr(blockId,solnData)
#endif
  do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
     ! Find a real difference between z's if problem is >= 3D
     if (NDIM > 2) then
        if (k .eq. 1) then
           dzz = zCoord(2) - zCoord(1) 
        else
           dzz = zCoord(k) - zCoord(k-1) 
        endif
     ! Otherwise this problem is <= 2D, so dzz is meaningless
     else
       dzz = 0.0
     endif
     zz = zCoord(k)
     
     do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)
        ! Find a real difference between y's if problem is >= 2D
        if (NDIM > 1) then
           if (j .eq. 1) then
              dyy = yCoord(2) - yCoord(1) 
           else
              dyy = yCoord(j) - yCoord(j-1) 
           endif
        ! Otherwise this problem is <= 1D, so dyy is meaningless
        else
          dyy = 0.0
        endif
        yy = yCoord(j)
        

        do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH, IAXIS)
           xx = xCoord(i)
           if (i .eq. 1) then
              dxx = xCoord(2) - xCoord(1) 
           else
              dxx = xCoord(i) - xCoord(i-1) 
           endif
           
           sumRho = 0.
           sumE   = 0.
           
           !
           !       Break the cell into sim_nSubZones^NDIM sub-zones, and look up the
           !       appropriate quantities along the 1d profile for that subzone.  
           !
           !       Have the final values for the zone be equal to the average of
           !       the subzone values.
           ! 

           do kk = 0, (sim_nSubZones-1)*K3D
              zz    = zCoord(k) + (kk*sim_inSubzm1-.5)*dzz 
              zDist = (zz - sim_zCenter) * K3D
              
              do jj = 0, (sim_nSubZones-1)*K2D
                 yy    = yCoord(j) + (jj*sim_inSubzm1-.5)*dyy
                 yDist = (yy - sim_yCenter) * K2D
                 
                 do ii = 0, (sim_nSubZones-1)
                    xx    = xCoord(i) + (ii*sim_inSubzm1-.5)*dxx
                    xDist = xx - sim_xCenter
                    
                    dist    = sqrt( xDist**2 + yDist**2 + zDist**2 )
                    distInv = 1. / max( dist, 1.E-10 )
                    call sim_find (radius, ipos, dist, jLo)
                    !
                    !  a point at `dist' is frac-way between jLo and jHi.   We do a
                    !  linear interpolation of the quantities at jLo and jHi and sum those.
                    ! 
                    if (jLo .eq. 0) then
                       jLo = 1
                       jHi = 1
                       frac = 0.
                    else if (jLo .eq. ipos) then
                       jLo = ipos
                       jHi = ipos
                       frac = 0.
                    else
                       jHi = jLo + 1
                       frac = (dist - radius(jLo)) / & 
                            (radius(jHi)-radius(jLo))
                    endif
                    ! 
                    !   Now total these quantities.   Note that  v is a radial velocity; 
                    !   we multiply by the tangents of the appropriate angles to get
                    !   the projections in the x, y and z directions.
                    !
                    sumE = sumE +  & 
                         eint(jLo) + frac*(eint(jHi) - eint(jLo))
                    
                    sumRho = sumRho + & 
                         rhop(jLo) + frac*(rhop(jHi) - rhop(jLo))
                 enddo
              enddo
           enddo
           
           xx = xCoord(i) - sim_xCenter
           yy = yCoord(j) - sim_yCenter
           zz = zCoord(k) - sim_zCenter
           dist = sqrt(xx**2 + yy**2 + zz**2)
           if (dist .le. radius(ipos)) then
              rho = max (sumRho * sim_inszd, sim_rhoAmbient)
           else
              rho = sim_rhoAmbient
           endif

           !p   = max (sumE   * sim_inszd * rho, sim_pAmbient)
           !if (dist .le. radius(ipos)) then
              vx = 0.0d0
              vy = 0.0d0
              vz = 0.0d0
           !else
           !   vtot = sqrt(g * 0.9 * 2.e33 / dist)
           !   theta = acos(zz/dist)
           !   phi = atan(yy/xx)
           !   vx = vtot*cos(phi)*sin(theta)
           !   vy = vtot*sin(phi)*sin(theta)
           !   vz = vtot*cos(theta)
           !endif

           ek  = 0.5*(vx*vx + vy*vy + vz*vz)
           t = 5.0e5
           e   = e + ek
           !e   = max (e, sim_smallP)
           
           axis(IAXIS)=i
           axis(JAXIS)=j
           axis(KAXIS)=k


#ifdef FL_NON_PERMANENT_GUARDCELLS
           solnData(DENS_VAR,i,j,k)=rho
           solnData(TEMP_VAR,i,j,k)=t
           solnData(VELX_VAR,i,j,k)=vx
           solnData(VELY_VAR,i,j,k)=vy
           solnData(VELZ_VAR,i,j,k)=vz
           do put=1,NSPECIES
              if (xn(SPECIES_BEGIN+put-1) == 0.0) xn(SPECIES_BEGIN+put-1) = sim_smallX
              solnData(SPECIES_BEGIN+put-1,i,j,k)=xn(SPECIES_BEGIN+put-1)
           enddo
#else
           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rho)
           call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, axis, t)    
           call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, vx)
           call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, vy)
           call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, axis, vz)
           do put=1,NSPECIES
              if (xn(SPECIES_BEGIN+put-1) == 0.0) xn(SPECIES_BEGIN+put-1) = sim_smallX
              call Grid_putPointData(blockID,CENTER,SPECIES_BEGIN+put-1,&
                   EXTERIOR,axis,xn(SPECIES_BEGIN+put-1))
           enddo
#endif
        enddo
     enddo
  enddo
#ifdef FL_NON_PERMANENT_GUARDCELLS
  call Grid_releaseBlkPtr(blockID, solnData)
#endif
  deallocate(xCoord)
  deallocate(yCoord)
  deallocate(zCoord)
  return
end subroutine Simulation_initBlock


!******************************************************************************

!  Routine:     sim_find()

!  Description: Given a monotonically increasing table x(nn) and a test value
!               x0, return the index i of the largest table value less than
!               or equal to x0 (or 0 if x0 < x(1)).  Use binary search.

subroutine sim_find (x, nn, x0, i)

  implicit none

! Arguments, LBR guessed intent on these
  integer, intent(IN) :: nn
  integer, intent(OUT):: i
  real, intent(IN)    :: x(nn), x0

! local variables
  integer  il, ir, im

  if (x0 .lt. x(1)) then

     i = 0

  elseif (x0 .gt. x(nn)) then

     i = nn

  else

     il = 1
     ir = nn
10   if (ir .eq. il+1) goto 20
     im = (il + ir) / 2
     if (x(im) .gt. x0) then
        ir = im
     else
        il = im
     endif
     goto 10
20   i = il

  endif

  return
end subroutine sim_find
