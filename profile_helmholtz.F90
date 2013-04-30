subroutine wd_intout(xn,rho0_cgs,T0,wd_mass_tot,cnt,r,rho,donor_dist,donor_mass)
!***************************************************************
!                                                              *
!     integrate outward until pressure gets negative,          *
!     give back total mass to correct guess for central        *
!     density; 15.6.98                                         *
!     Version for Helmholtz-EOS; SKR 22.1.2003                 *
!                                                              *
!***************************************************************
      
      use Eos_interface, ONLY : Eos
      use Grid_interface, ONLY : Grid_getMinCellSize
      use PhysicalConstants_interface, ONLY : PhysicalConstants_get
      use RuntimeParameters_interface, ONLY : RuntimeParameters_mapStrToInt, RuntimeParameters_get
      implicit none
#include "Eos.h"
#include "constants.h"
#include "Flash.h"

      integer ngp
      parameter (ngp=100000)
      double precision, intent(in) :: T0, donor_dist, donor_mass
      double precision, intent(in), dimension(NSPECIES) :: xn
      double precision, intent(out) :: wd_mass_tot, r(ngp), rho(ngp)
      integer, intent(out) :: cnt

      integer i
      double precision rho_cgs_cut
      double precision m(ngp),pr(ngp),dr,ri,pri,dthick
      double precision rhoi,rho0_cgs
      double precision p0,deltm,gg
      double precision, dimension(EOS_NUM) :: eosData

      call RuntimeParameters_get('sim_rhoAmbient',rho_cgs_cut)
      call PhysicalConstants_get("Newton", gg)
      call Grid_getMinCellSize(dthick)

      eosData(EOS_DENS) = rho0_cgs
      eosData(EOS_TEMP) = T0
      call Eos(MODE_DENS_TEMP,1,eosData,xn)

      p0 = eosData(EOS_PRES)
      pr(1) = p0
      rho(1) = rho0_cgs
      m(1) = 0.
      cnt = 0
      !write(*,*) 'central pr, rho, m', pr(0), rho(0), m(0)
!
!--give grid size
!     
      dr = dthick/100.
      
!     calculate masses and densities at grid points
      do i=2,ngp
         cnt= cnt+1
         ri= (i-1)*dr
         r(i)= ri
        
!     mass inside grid point i
         deltm= 4.*PI*ri**2*rho(i-1)*dr
         m(i)= deltm + m(i-1)

!     pressure at grid point i
         !Include gravity from companion
         !pr(i)=  pr(i-1) - gg*rho(i-1)*dr*(m(i)/ri**2. - donor_mass/(donor_dist-ri)**2.)
         pr(i)=  pr(i-1) - gg*rho(i-1)*dr*m(i)/ri**2.
         if(pr(i).lt.0.) exit
        
!     solve EOS for corresponding density 
         pri= pr(i)
         rhoi= rho(i-1)
         call calcrho_helm(pri,T0,xn,rhoi)
         if(rhoi.lt.rho_cgs_cut) exit

         rho(i)= rhoi
      enddo
   
      wd_mass_tot= m(cnt)

      return
end


subroutine nozz_intout(xn,rho0_cgs,T0,m_nozz_tot,donor_mass,cnt,r,rho,xnozz,ynozz,znozz,grv_ptxpos,grv_ptypos,grv_ptzpos)
!***************************************************************
!                                                              *
!     integrate outward until pressure gets negative,          *
!     give back total mass to correct guess for central        *
!     density; 15.6.98                                         *
!     Version for Helmholtz-EOS; SKR 22.1.2003                 *
!                                                              *
!***************************************************************
      
      use Eos_interface, ONLY : Eos
      use Grid_interface, ONLY : Grid_getMinCellSize
      use gr_mpoleData, ONLY: Mtot, Xcm, Ycm, Zcm
      use PhysicalConstants_interface, ONLY : PhysicalConstants_get
      use RuntimeParameters_interface, ONLY : RuntimeParameters_mapStrToInt, RuntimeParameters_get
      implicit none
#include "Eos.h"
#include "constants.h"
#include "Flash.h"

      integer ngp
      parameter (ngp=10000)
      double precision, intent(in) :: xnozz, ynozz, znozz, grv_ptxpos, grv_ptypos, grv_ptzpos
      double precision, intent(in) :: T0, donor_mass
      double precision, intent(in), dimension(NSPECIES) :: xn
      double precision, intent(out) :: m_nozz_tot, r(ngp), rho(ngp)
      integer, intent(out) :: cnt

      integer i
      double precision rho_cgs_cut
      double precision m(ngp),pr(ngp),dr,ri,pri,dthick
      double precision rhoi,rho0_cgs,d12,d_db,omega2,theta,r_bary
      double precision p0,deltm,d1,d2,local_d1,local_d2,m1,m2,gg,accretor_mass
      double precision, dimension(EOS_NUM) :: eosData

!     fix temperature to T0
      
      call RuntimeParameters_get('sim_rhoAmbient',rho_cgs_cut)
      call PhysicalConstants_get("Newton", gg)
      call Grid_getMinCellSize(dthick)

      eosData(EOS_DENS) = rho0_cgs
      eosData(EOS_TEMP) = T0
      call Eos(MODE_DENS_TEMP,1,eosData,xn)

      p0 = eosData(EOS_PRES)
      pr(1) = p0
      rho(1) = rho0_cgs
      m(1) = 0.
      cnt = 0
      !write(*,*) 'central pr, rho, m', pr(0), rho(0), m(0)

      m1 = donor_mass
      m2 = Mtot
      d1 = sqrt((xnozz - grv_ptxpos)**2. + (ynozz - grv_ptypos)**2. + (znozz - grv_ptzpos)**2.)
      d2 = sqrt((xnozz - Xcm)**2. + (ynozz - Ycm)**2. + (znozz - Zcm)**2.)

      d12 = sqrt((Xcm - grv_ptxpos)**2. + (Ycm - grv_ptypos)**2. + (Zcm - grv_ptzpos)**2.)
      d_db = m2 / (m1 + m2) * d12 !center of mass coordinate of binary along line between centers of masses of the stars
      d_db = d_db*cos(atan(abs(Zcm - grv_ptzpos)/sqrt((Xcm - grv_ptxpos)**2.+(Ycm - grv_ptypos)**2.))) !distance of barycenter from donor in XY plane
      omega2 = 1. / (d12**3. / gg / (m1 + m2))
!
!--give grid size (cgs); for normal mass WDs dr=1.e4 is a good choice
!     
      dr = dthick/100.
      
!     calculate masses and densities at grid points
      !L1_pot = -gg*(m1/d1+m2/d2)-0.5*omega**2.*(d1 - d_db)**2.
      do i=2,ngp
         cnt= cnt+1
         ri= (i-1)*dr
         !x_it = xnozz
         !do while (abs(pot-L1_pot)/L1_pot .gt. 1e-6)
         !    d1_it = sqrt((x - grv_ptxpos)**2. + (ynozz + ri - grv_ptypos)**2. + (znozz - grv_ptzpos)**2.)
         !    d2_it = sqrt((x - Xcm)**2. + (ynozz + ri - Ycm)**2. + (znozz - Zcm)**2.)
         !    d_db_it = sqrt(
         !    pot = -gg*(m1/d1_it + m2/d2_it)-0.5*omega**2.*((d1_it-d_db)**2. + ri**2.)
         !    gpot = -
         !enddo
         local_d1 = d1 - ri/tan(1.0)
         local_d2 = d2 + ri/tan(1.0)
         theta = atan(ri/(local_d1 - d_db))
         r_bary = sqrt((local_d1 - d_db)**2. + ri**2.)
         r(i)= ri
        
!     mass inside grid point i
         deltm = 2.*PI*ri*rho(i-1)*dr*dthick
         m(i)= deltm + m(i-1)

!     pressure at grid point i
         pr(i) = pr(i-1) - gg*rho(i-1)*dr*ri* &
             (m1/(local_d1**2.+ri**2.)**1.5 + m2/(local_d2**2.+ri**2.)**1.5) + &
             rho(i-1)*dr*omega2*r_bary*sin(theta)
!         pr(i)=  pr(i-1) - gg*rho(i-1)*m(i)*dr/ri**2 
         if(pr(i).lt.0.) exit
        
!     solve EOS for corresponding density 
         pri= pr(i)
         rhoi= rho(i-1)
         call calcrho_helm(pri,T0,xn,rhoi)
         if(rhoi.lt.rho_cgs_cut) exit

         rho(i)= rhoi
      enddo
   
      m_nozz_tot= m(cnt)

      return
end


subroutine calcrho_helm(presi,Ti,xn,rhoi)
!************************************************************
!                                                           *
!     given temperature and desired pressure calculate      *
!     necessary rho; SKR 22.2.2003                          *
!                                                           *
!************************************************************

      use Eos_interface, ONLY : Eos
      implicit none
#include "Eos.h"
#include "constants.h"
#include "Flash.h"

      integer maxit
      double precision rho_tol,rhomin
      parameter(rho_tol=1.d-8,rhomin=1.d-3,maxit=100)
      double precision, intent(in) :: presi,Ti
      double precision, intent(in), dimension(NSPECIES) :: xn
      double precision, intent(out) :: rhoi

      integer it
      double precision abar,zbar,rho,pdes,rho_new,error
      double precision pres,dpd

      double precision, dimension(EOS_NUM) :: eosData
      logical, dimension(EOS_VARS+1:EOS_NUM) :: mask = .false.

      mask(EOS_DPD) = .true.

!     desired pressure 
      pdes= presi

      eosData(EOS_TEMP) = Ti
      eosData(EOS_DENS) = rhoi
      call Eos(MODE_DENS_TEMP,1,eosData,xn,mask)
      pres = eosData(EOS_PRES) 
      dpd = eosData(EOS_DPD)
      rho = rhoi
      
!..here the pipeline is only 1 element long

!     now do iteration
      do it=1,maxit
         !print *, 'pres, dpd, rho, pdes', pres, dpd, rho, pdes
         rho_new = rho - (pres-pdes)/dpd

!     Dont allow the pressure to change by more than an 
!     order of magnitude in a single iteration
         if(rho_new.gt.10.e0*rho) rho_new = 10.e0*rho
         if(rho_new.lt.0.1e0*rho) rho_new = 0.1e0*rho

!     Compute the error
         error = abs((rho_new-rho)/rho)
         
!     Store the new temperature
         rho = rho_new
         eosData(EOS_DENS) = rho

         if (rho.le.rhomin) then
            rho = rhomin
            error = 0.1*rho_tol
         endif
         if (error.lt.rho_tol) exit

!     call HELMHOLTZ-EOS (all output quantities are stored in'vector_eos.dek')
         call Eos(MODE_DENS_TEMP,1,eosData,xn,mask)
         pres = eosData(EOS_PRES) 
         dpd = eosData(EOS_DPD)
      enddo

      if(it.ge.maxit)then
         write(*,*) 'calc_rho_helm: no convergence after', maxit, 'iterations'
         write(*,*) 'current values are (rho, pres, dpd, error): '
         write(*,*) rho, pres, dpd, error
         write(*,*) 'xn'
         write(*,*) xn
         stop
      endif

!     If converged:
      rhoi = rho_new

      return 
end

subroutine wd_guess_rho_0(mdes,rho0_cgs)

!******************************************************************
!                                                                 *
!     provide guess value for central density of WD of mass mdes  *
!     SKR 22.7.2005                                               *
!                                                                 *
!******************************************************************

      implicit none

      integer ilo,ihi,i
      double precision, intent(in) :: mdes
      double precision, intent(out) :: rho0_cgs
      double precision mgrid(5),lrhogrid(5)     
      double precision dlrhodm,lrho0_cgs

!     mass grid
      mgrid(1)=0.02*2d33
      mgrid(2)=0.2*2d33
      mgrid(3)=0.6*2d33
      mgrid(4)=1.2*2d33
      mgrid(5)=1.35*2d33

!     log-density grid
      lrhogrid(1)=3.26
      lrhogrid(2)=5.28
      lrhogrid(3)=6.51
      lrhogrid(4)=8.07
      lrhogrid(5)=8.88

!     get lower index
      
      if (mdes .lt. mgrid(1)) then
          lrho0_cgs = lrhogrid(1)
      elseif (mdes .gt. mgrid(size(mgrid))) then
          lrho0_cgs = lrhogrid(5)
      else
          do i = 1, size(mgrid)-1
              if (mdes .gt. mgrid(i)) then
                  dlrhodm= (lrhogrid(i+1)-lrhogrid(i))/(mgrid(i+1)-mgrid(i))
                  lrho0_cgs= lrhogrid(i)+ dlrhodm*(mdes-mgrid(i))
              endif
          enddo
      endif
      rho0_cgs= 10**lrho0_cgs

      return
end


subroutine nozz_guess_rho_0(mdes,rho0_cgs)

!******************************************************************
!                                                                 *
!     provide guess value for central density of WD of mass mdes  *
!     SKR 22.7.2005                                               *
!                                                                 *
!******************************************************************

      implicit none

      integer ilo,ihi,i
      double precision, intent(in) :: mdes
      double precision, intent(out) :: rho0_cgs
      double precision mgrid(6),lrhogrid(6)     
      double precision dlrhodm,lrho0_cgs

!     mass grid
      mgrid(1)=1.d-12*2d33
      mgrid(2)=1.d-10*2d33
      mgrid(3)=1.d-8*2d33
      mgrid(4)=1.d-6*2d33
      mgrid(5)=1.d-4*2d33
      mgrid(6)=1.d-2*2d33

!     log-density grid
      lrhogrid(1)=0.15
      lrhogrid(2)=1.40
      lrhogrid(3)=2.61
      lrhogrid(4)=3.81
      lrhogrid(5)=5.00
      lrhogrid(6)=6.1

!     get lower index
      
      if (mdes .lt. mgrid(1)) then
          lrho0_cgs = log10(mdes/2d33)/2. + 7.
      elseif (mdes .gt. mgrid(size(mgrid))) then
          lrho0_cgs = log10(mdes/2d33)/2. + 7.
      else
          do i = 1, size(mgrid)-1
              if (mdes .gt. mgrid(i)) then
                  dlrhodm= (lrhogrid(i+1)-lrhogrid(i))/(mgrid(i+1)-mgrid(i))
                  lrho0_cgs= lrhogrid(i)+ dlrhodm*(mdes-mgrid(i))
              endif
          enddo
      endif
      rho0_cgs= 10**lrho0_cgs

      return
end

