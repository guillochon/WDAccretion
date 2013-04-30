!!****if* source/Grid/GridMain/paramesh/Grid_markRefineDerefine
!!
!! NAME
!!  Grid_markRefineDerefine
!!
!! SYNOPSIS
!!
!!  Grid_markRefineDerefine(integer(IN) :: myPE)
!!  
!! DESCRIPTION 
!!  Mark blocks for refinement or derefinement
!!  This routine is used with AMR only where individual 
!!  blocks are marked for refinement or derefinement based upon
!!  some refinement criterion. The Uniform Grid does not need
!!  this routine, and uses the stub.
!!
!! ARGUMENTS
!!  myPE : my processor id.
!! 
!! NOTES
!!
!! Every unit uses a few unit scope variables that are
!! accessible to all routines within the unit, but not to the
!! routines outside the unit. For Grid unit these variables begin with "gr_"
!! like, gr_myPE or gr_eosMode, and are stored in fortran
!! module Grid_data (in file Grid_data.F90). The other variables
!! are local to the specific routines and do not have the prefix "gr_"
!!
!!
!!***

subroutine Grid_markRefineDerefine(MyPE)

  use Grid_data, ONLY : gr_refine_cutoff, gr_derefine_cutoff,&
                        gr_refine_filter,gr_refine_val_cutoff,&
                        gr_numRefineVars,gr_refine_var,gr_refineOnParticleCount,&
                        gr_blkList
  use tree, ONLY : newchild, refine, derefine, stay, nodetype,&
       lrefine,lrefine_max, parent, nchild,child
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_getListOfBlocks,&
                             Grid_releaseBlkPtr,Grid_fillGuardCells
  use Driver_data, ONLY: dr_simTime, dr_restart, dr_initialSimTime
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Simulation_data, ONLY: sim_tRelax, wd_radius, sim_donorMass, &
      sim_tableRows, sim_tableCols, sim_table, sim_tInitial, sim_tNozz, &
      sim_xCenter, sim_yCenter, sim_zCenter
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get

  implicit none

#include "constants.h"
#include "Flash_mpi.h"
#include "Flash.h"

  integer,intent(IN) :: MyPE
  
  real, dimension(:,:,:,:), pointer :: solnData
  real :: ref_cut,deref_cut,ref_filter,ref_val_cut
  real :: xcenter,ycenter,zcenter,don_dist,acc_rate
  real :: maxptime,shrinkdur,box,zbox,don_mass,nozz_dist
  real :: don_pos(3), acc_pos(3), noz_pos(3)
  integer       :: l,iref,blkCount,lb,i,j
  logical :: doEos=.true.
  integer,parameter :: maskSize = NUNK_VARS+NDIM*NFACE_VARS
  logical,dimension(maskSize) :: gcMask
  real :: maxvals(MAXBLOCKS),maxvals_parent(MAXBLOCKS)
  integer :: nsend,nrecv,ierr
  integer :: reqr(MAXBLOCKS),reqs(MAXBLOCKS)
  integer :: statr(MPI_STATUS_SIZE,MAXBLOCKS)
  integer :: stats(MPI_STATUS_SIZE,MAXBLOCKS)
  integer :: cur_max_refine

  ! that are implemented in this file need values in guardcells

  !if ((.not. dr_restart) .or. dr_simTime > dr_initialSimTime + 0.1) then
      gcMask=.false.
      do i = 1,gr_numRefineVars
         iref = gr_refine_var(i)
         if (iref > 0) gcMask(iref) = .TRUE.
      end do

      call Grid_fillGuardCells(MyPE,CENTER_FACES,ALLDIR,doEos=.true.,&
           maskSize=maskSize, mask=gcMask, makeMaskConsistent=.true.)

      newchild(:) = .FALSE.
      refine(:)   = .FALSE.
      derefine(:) = .FALSE.
      stay(:)     = .FALSE.

      !do l = 1,gr_numRefineVars
      !   iref = gr_refine_var(l)
      !   ref_cut = gr_refine_cutoff(l)
      !   deref_cut = gr_derefine_cutoff(l)
      !   ref_filter = gr_refine_filter(l)
      !   ref_val_cut = gr_refine_val_cutoff(l)
      !   call gr_markRefineDerefine(MyPE,iref,ref_cut,deref_cut,ref_filter,&
      !       ref_val_cut)
      !end do

      call Grid_getListOfBlocks(ACTIVE_BLKS, gr_blkList,blkCount)

      cur_max_refine = 0
      do i = 1, blkCount
         lb = gr_blkList(i)
         if (nodetype(lb) == LEAF .and. lrefine(lb) .gt. cur_max_refine) then
            cur_max_refine = lrefine(lb)
         endif
      enddo

      call gr_markInRectangle(sim_xCenter-4d9,sim_xCenter,sim_yCenter-1.d9,&                                               
           sim_yCenter+1.d9,sim_zCenter-1.d9,sim_zCenter+1.d9,lrefine_max-2,0)                                              
      if (dr_simTime .ge. sim_tInitial + sim_tNozz) then
          call wd_interp(don_mass, don_pos, acc_pos, noz_pos, don_dist, acc_rate)
          call gr_markInRadius(sim_xCenter,sim_yCenter,sim_zCenter,abs(sim_xCenter-noz_pos(1)),lrefine_max-1,0)                                  
      endif
      if (dr_simTime .ge. sim_tInitial + sim_tRelax - 1.0) &                                           
          call gr_markInRectangle(sim_xCenter-2.5d9,sim_xCenter,sim_yCenter-6.d8,&                                         
               sim_yCenter+6.d8,sim_zCenter-6.d8,sim_zCenter+6.d8,lrefine_max,0) !Refine 1 second before the nozzle turns on
      call gr_markInRadius(sim_xCenter,sim_yCenter,sim_zCenter,1.3*wd_radius,lrefine_max,0)                                
  !endif
    
  !call gr_markInRadius(xcenter,ycenter,zcenter,5.75e8,lrefine_max-1,1)
  !!call gr_markInRadius(xcenter,ycenter,zcenter,2.0e7,lrefine_max,0)

  !if (t .eq. 0.0) then
  !    !write(*,*) 'entered sphere refine'
  !    call gr_markInRadius(xcenter,ycenter,zcenter,1.5e9,lrefine_max)
  !else
  !    do l = 1,gr_numRefineVars
  !       ref_val_cut = gr_refine_val_cutoff(l)
  !       call gr_markVarThreshold(iref,ref_val_cut,0,0)
  !    enddo

  !    do i = 1, blkCount
  !       lb = gr_blkList(i)
  !       call Grid_getBlkPtr(lb,solnData,CENTER)
  !       maxvals(lb) = maxval(solnData(iref,:,:,:))
  !       call Grid_releaseBlkPtr(lb,solnData)
  !    end do

! !    Communicate maxvals of parents to their leaf children.
! !    Maximally refined children collect messages from parents.

  !    maxvals_parent(:) = 0.0
  !    nrecv = 0
  !    do i = 1, blkCount
  !       lb = gr_blkList(i)
  !       if (nodetype(lb) == LEAF .AND. lrefine(lb) == lrefine_max) then
  !          if(parent(1,lb).gt.-1) then
  !             if (parent(2,lb).ne.MyPE) then
  !                nrecv = nrecv + 1
  !                call MPI_IRecv(maxvals_parent(lb),1, &
  !                   FLASH_REAL, &
  !                   parent(2,lb), &
  !                   lb, &
  !                   MPI_COMM_WORLD, &
  !                   reqr(nrecv), &
  !                   ierr)
  !             else
  !                maxvals_parent(lb) = maxvals(parent(1,lb))
  !             end if
  !          end if
  !       end if
  !    end do

  !    ! parents send maxvals to children

  !    nsend = 0
  !    do i = 1, blkCount
  !       lb = gr_blkList(i)
  !       if (nodetype(lb) == PARENT .AND. lrefine(lb) == lrefine_max-1) then
  !          do j = 1,nchild
  !             if(child(1,j,lb).gt.-1) then
  !                if (child(2,j,lb).ne.MyPE) then
  !                   nsend = nsend + 1
  !                   call MPI_ISend(maxvals(lb), &
  !                      1, &
  !                      FLASH_REAL, &
  !                      child(2,j,lb), &  ! PE TO SEND TO
  !                      child(1,j,lb), &  ! THIS IS THE TAG
  !                      MPI_COMM_WORLD, &
  !                      reqs(nsend), &
  !                      ierr)
  !                end if
  !             end if
  !          end do
  !       end if
  !    end do

  !    if (nsend.gt.0) then
  !       call MPI_Waitall (nsend, reqs, stats, ierr)
  !    end if
  !    if (nrecv.gt.0) then
  !       call MPI_Waitall (nrecv, reqr, statr, ierr)
  !    end if

!!!      maxvals_parent(:) = 0.0  ! <-- uncomment line for previous behavior
  !    do i = 1, blkCount
  !       lb = gr_blkList(i)
  !       if (nodetype(lb) == LEAF) then
  !          if (maxvals(lb) < gr_refine_val_cutoff(1)) then
  !             refine(lb)   = .false.
! !              if (maxvals_parent(lb) < gr_refine_val_cutoff .AND. .NOT. stay(lb)) derefine(lb)   = .true.
  !             if (maxvals_parent(lb) < gr_refine_val_cutoff(1)) derefine(lb)   = .true.
  !          endif
  !       end if
  !    enddo

  !endif

  return
end subroutine Grid_markRefineDerefine

