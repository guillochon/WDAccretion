!!****if* source/Simulation/SimulationMain/Sedov/Simulation_data
!!
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  use Simulation_data 
!!
!!  DESCRIPTION
!!
!!  Stores the local data for Simulation setup: Sedov
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
!!
!!***

module Simulation_data

  implicit none
#include "constants.h"
#include "Flash.h"

  !! *** Runtime Parameters *** !!

  double precision, save    :: sim_tAmbient, sim_rhoAmbient
  double precision, save    :: sim_gamma, sim_xCenter, sim_yCenter, sim_zCenter
  double precision, save    :: sim_smallX, sim_smallRho, sim_smallP, sim_pi
  double precision, save    :: sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin, sim_zMax
  integer, save             :: sim_nSubZones
  integer, save             :: sim_tableRows, sim_tableCols
  double precision, save    :: sim_tInitial
  double precision, dimension(:,:), allocatable, save :: sim_table
  character(len=MAX_STRING_LENGTH), save :: sim_dataFile

  !! *** Variables pertaining to this Simulation *** !!

  integer, parameter        :: np = 100000
  double precision, save    :: sim_inSubZones, sim_inSubzm1
  double precision, save    :: sim_inszd, sim_tNozz, sim_tRelax, sim_relaxRate
  double precision, save    :: sim_accMass, sim_accTemp, sim_donorMass, sim_donorTemp, wd_radius

  double precision, save    :: radius(np), rhop(np), eint(np)
  double precision, save    :: wd_xn(SPECIES_BEGIN:SPECIES_END)
  integer, save             :: ipos

end module Simulation_data
