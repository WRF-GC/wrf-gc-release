!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   WRF-GCHP
!   GEOS-Chem High Performance-powered Chemistry Add-On for WRF Model
!
!   WRF & GCHP are (c) their original authors.
!   WRF-GCHP coupling layer (WGCL) is (c) Atmospheric Chemistry and Climate Group, Peking University
!
!   Developed by Haipeng Lin <linhaipeng@pku.edu.cn>, Xu Feng, 2018-01
!   Peking University, School of Physics
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   Codename Pumpkin: Abstracted Bindings for Chemistry-to-WRF
!
!   This Chemical Interface (chem/) is written after comprehensive study of
!   the original chem_driver.f from WRF-Chem v3.6.1
!   which is (c) their respective authors.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  MODULE: module_gocart_aerosols
!  DESCRIPTION: A dummy gocart aerosols subroutine to satisfy WRF_CHEM = 1
!               triggering all sorts of weird compilation links in WRF.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE module_gocart_aerosols
   USE module_model_constants, only:mwdry
   INTEGER, PARAMETER ::NBC1 = 1, NOC1 = 2, NBC2 = 3, NOC2 = 4

CONTAINS
   subroutine gocart_aerosols_driver(ktau, dt, config_flags, t_phy, moist, &
                                     chem, rho_phy, dz8w, p8w, dx, g, &
                                     ids, ide, jds, jde, kds, kde, &
                                     ims, ime, jms, jme, kms, kme, &
                                     its, ite, jts, jte, kts, kte)
      use module_configure
      use module_state_description
      implicit none
      type(grid_config_rec_type), intent(in)    :: config_flags

      integer, intent(in) :: ktau, &
                             ids, ide, jds, jde, kds, kde, &
                             ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte
      real, dimension(ims:ime, kms:kme, jms:jme, num_moist), &
         intent(in) :: moist
      real, dimension(ims:ime, kms:kme, jms:jme, num_chem), &
         intent(inout) :: chem
      real, dimension(ims:ime, kms:kme, jms:jme), &
         intent(in) :: t_phy, dz8w, p8w, rho_phy
      real, intent(in) :: dt, dx, g

      !! Stub
   end subroutine gocart_aerosols_driver

   subroutine sum_pm_gocart( &
      alt, chem, pm2_5_dry, pm2_5_dry_ec, pm10, &
      ids, ide, jds, jde, kds, kde, &
      ims, ime, jms, jme, kms, kme, &
      its, ite, jts, jte, kts, kte)

      USE module_configure
      USE module_state_description
      implicit none

      integer, intent(in) :: ids, ide, jds, jde, kds, kde, &
                             ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte
      real, dimension(ims:ime, kms:kme, jms:jme), &
         intent(inout) :: pm2_5_dry, pm2_5_dry_ec, pm10
      real, dimension(ims:ime, kms:kme, jms:jme), &
         intent(in) :: alt
      real, dimension(ims:ime, kms:kme, jms:jme, num_chem), &
         intent(in) :: chem

      !! Stub
   end subroutine sum_pm_gocart

END MODULE module_gocart_aerosols
