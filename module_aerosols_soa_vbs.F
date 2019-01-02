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
!  MODULE: module_aerosols_soa_vbs
!  DESCRIPTION: A dummy SOA/VBS aerosols subroutine to satisfy WRF_CHEM = 1
!               triggering all sorts of weird compilation links in WRF.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define cw_species_are_in_registry

module module_aerosols_soa_vbs
  implicit none

  contains

  ! sum_pm_soa_vbs
  ! req. by dyn_em/start_em, stubbed below
  subroutine sum_pm_soa_vbs(alt, chem, h2oaj, h2oai, &
      pm2_5_dry, pm2_5_water, pm2_5_dry_ec, pm10, dust_opt, &
      ids, ide, jds, jde, kds, kde, &
      ims, ime, jms, jme, kms, kme, &
      its, ite, jts, jte, kts, kte)
      use module_state_description, only: num_chem

      integer, intent(in) :: dust_opt, &
                             ids, ide, jds, jde, kds, kde, &
                             ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte

      real, dimension(ims:ime, kms:kme, jms:jme, num_chem), &
         intent(in) :: chem

      real, dimension(ims:ime, kms:kme, jms:jme), &
         intent(in) :: alt, h2oaj, h2oai

      real, dimension(ims:ime, kms:kme, jms:jme), &
         intent(out) :: pm2_5_dry, pm2_5_water, pm2_5_dry_ec, pm10

   end subroutine sum_pm_soa_vbs

end module module_aerosols_soa_vbs