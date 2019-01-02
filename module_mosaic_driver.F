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
!  MODULE: module_mosaic_driver
!  DESCRIPTION: A dummy MOSAIC driver to satisfy WRF_CHEM = 1 requirement when
!               compiling WRF (namely, routines in em_real)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module module_mosaic_driver

#define cw_species_are_in_registry

contains
   subroutine mosaic_aerchem_driver(id, curr_secs, ktau, dtstep, ktauc, dtstepc, config_flags, &
      t_phy, rho_phy, p_phy, &
      moist, chem, vbs_nbin, &
      ids, ide, jds, jde, kds, kde, &
      ims, ime, jms, jme, kms, kme, &
      its, ite, jts, jte, kts, kte)

      use module_configure, only:grid_config_rec_type, &
         p_qv, &
         p_so2, p_ho2, p_so4aj, p_corn, p_hcl, p_mtf, &
         p_so4_a01, p_water_a01, p_num_a01, &
         p_so4_a04, p_water_a04, p_num_a04

      use module_state_description, only: num_moist, num_chem

      implicit none

      integer, intent(in) :: id, ktau, ktauc, &
         ids, ide, jds, jde, kds, kde, &
         ims, ime, jms, jme, kms, kme, &
         its, ite, jts, jte, kts, kte, &
         vbs_nbin(1)
      !   id - domain index
      !   ktau - time step number
      !   ktauc - gas and aerosol chemistry time step number

      !   [ids:ide, kds:kde, jds:jde] - spatial (x,z,y) indices for "domain"
      !   [ims:ime, kms:kme, jms:jme] - spatial (x,z,y) indices for "memory"
      !   Most arrays that are arguments to chem_driver
      !   are dimensioned with these spatial indices.
      !   [its:ite, kts:kte, jts:jte] - spatial (x,z,y) indices for "tile"
      !   chem_driver and routines under it do calculations
      !   over these spatial indices.

      real(kind=8), intent(in) :: curr_secs
      real, intent(in) :: dtstep, dtstepc
      !   dtstep - main model time step (s)
      !   dtstepc - time step for gas and aerosol chemistry(s)

      real, intent(in), &
         dimension(ims:ime, kms:kme, jms:jme) :: &
         t_phy, rho_phy, p_phy
      !   t_phy - temperature (K)
      !   rho_phy - air density (kg/m^3)
      !   p_phy - air pressure (Pa)

      real, intent(in), &
         dimension(ims:ime, kms:kme, jms:jme, 1:num_moist) :: moist
      !   moist - mixing ratios of moisture species (water vapor,
      !   cloud water, ...) (kg/kg for mass species, #/kg for number species)

      real, intent(inout), &
         dimension(ims:ime, kms:kme, jms:jme, 1:num_chem) :: chem
      !   chem - mixing ratios of trace gase (ppm) and aerosol species
      !   (ug/kg for mass species, #/kg for number species)

      type(grid_config_rec_type), intent(in) :: config_flags
      !   config_flags - configuration and control parameters

      return
   end subroutine mosaic_aerchem_driver


   ! sum_pm_mosaic stub
   ! Called by start_em.F in CBMZ_MOSAIC_ ... options
   subroutine sum_pm_mosaic(alt, chem, &
      pm2_5_dry, pm2_5_water, pm2_5_dry_ec, pm10, &
      ids, ide, jds, jde, kds, kde, &
      ims, ime, jms, jme, kms, kme, &
      its, ite, jts, jte, kts, kte)

      use module_state_description, only: num_chem
      implicit none

      integer, intent(in) :: &
         ids, ide, jds, jde, kds, kde, &
         ims, ime, jms, jme, kms, kme, &
         its, ite, jts, jte, kts, kte

      real, dimension(ims:ime, kms:kme, jms:jme), &
         intent(in) :: alt

      real, dimension(ims:ime, kms:kme, jms:jme, num_chem), &
         intent(in) :: chem

      real, dimension(ims:ime, kms:kme, jms:jme), &
         intent(out) :: pm2_5_dry, pm2_5_water, pm2_5_dry_ec, pm10
   end subroutine sum_pm_mosaic

   subroutine aerchem_debug_dump(iflag, iclm, jclm, dtchem)
      implicit none
      integer iflag, iclm, jclm
      real dtchem

      return
   end subroutine aerchem_debug_dump

end module module_mosaic_driver
