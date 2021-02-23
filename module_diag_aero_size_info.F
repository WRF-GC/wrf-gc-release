!---------------------------------------------------------------------------
! This module is used to diagnose the number concentrations and mass 
! concnetrations in each size bin based on the MOSAIC 4 (default)or 8 bins. 
!---------------------------------------------------------------------------
  module module_diag_aero_size_info
    
    USE module_state_description
    USE module_model_constants
    USE module_configure

    contains

    subroutine diag_aero_wrfgc_aercld_ptrs(num_chem, is_aerosol, config_flags)

    ! This subroutine is used to provide a means of organizing and accessing
    ! aerosol species in the "chem" array by their chemical component, 
    ! size bin (or mode), "type", and "phase" based on the module_data_geoschem_asect.

    !            USE module_state_description
    !            USE module_model_constants
    !            USE module_configure
                USE module_data_gigc_asect

                implicit none
                integer, intent(in)  :: num_chem
		logical, intent(out) :: is_aerosol(num_chem)
                type (grid_config_rec_type), intent(in) :: config_flags
        
        ! local variables	
  		real    :: dlo, dhi, xlo, xhi, dx
                integer :: iphase, isize, itype, n
                
                nphase_aer   = 2
                ai_phase     = 1
                cw_phase     = 2
                
                ntype_aer    = 1
                nsize_aer(1) = 4
                ncomp_aer(1) = 10
                
                ! densities of chemical components
                itype = 1
                dens_aer(1, itype) = dens_so4_aer_gc
                dens_aer(2, itype) = dens_nit_aer_gc
                dens_aer(3, itype) = dens_nh4_aer_gc
                dens_aer(4, itype) = dens_oc_aer_gc
                dens_aer(5, itype) = dens_oc_aer_gc
                dens_aer(6, itype) = dens_bc_aer_gc
                dens_aer(7, itype) = dens_bc_aer_gc
                dens_aer(8, itype) = dens_seas_aer_gc
                dens_aer(9, itype) = dens_dst2_aer_gc
                dens_aer(10,itype) = dens_soas_aer_gc 

                ! molecular weight of chemical components
                mw_aer(1, itype) = mw_so4_aer
                mw_aer(2, itype) = mw_nit_aer
                mw_aer(3, itype) = mw_nh4_aer
                mw_aer(4, itype) = mw_oc_aer
                mw_aer(5, itype) = mw_oc_aer
                mw_aer(6, itype) = mw_bc_aer
                mw_aer(7, itype) = mw_bc_aer
                mw_aer(8, itype) = mw_seas_aer
                mw_aer(9, itype) = mw_dst_aer
                mw_aer(10,itype) = mw_soas_aer

                ! hygroscopicity of chemical components
                hygro_aer(1, itype) = hygro_so4_aer
                hygro_aer(2, itype) = hygro_nit_aer
                hygro_aer(3, itype) = hygro_nh4_aer
                hygro_aer(4, itype) = hygro_oc_aer
                hygro_aer(5, itype) = hygro_oc_aer
                hygro_aer(6, itype) = hygro_bc_aer
                hygro_aer(7, itype) = hygro_bc_aer
                hygro_aer(8, itype) = hygro_seas_aer
                hygro_aer(9, itype) = hygro_dst_aer
                hygro_aer(10,itype) = hygro_soas_aer

                dlo   = 0.0390625 * 1.0E-4 ! unit:cm
                dhi   = 10.0 * 1.0E-4      ! unit:cm
                xlo   = log( dlo )
                xhi   = log( dhi )
                dx    = (xhi - xlo)/nsize_aer(itype)
                do isize = 1, nsize_aer(itype)
                   dlo_sect(isize, itype)   = exp( xlo + dx*(isize-1) )
                   dhi_sect(isize, itype)   = exp( xlo + dx*isize )
                   dcen_sect(isize, itype)  = sqrt( dlo_sect(isize, itype) * dhi_sect(isize, itype))
                   sigmag_aer(isize, itype) = (dhi_sect(isize, itype)/dlo_sect(isize, itype))**0.289
                end do


                itype = 1
                isize = 1
                numptr_aer( isize, itype, ai_phase)    = p_diag_num_a1
                massptr_aer(1, isize, itype, ai_phase) = p_diag_so4_a1
                massptr_aer(2, isize, itype, ai_phase) = p_diag_nit_a1
                massptr_aer(3, isize, itype, ai_phase) = p_diag_nh4_a1
                massptr_aer(4, isize, itype, ai_phase) = p_diag_ocpi_a1
                massptr_aer(5, isize, itype, ai_phase) = p_diag_ocpo_a1
                massptr_aer(6, isize, itype, ai_phase) = p_diag_bcpi_a1
                massptr_aer(7, isize, itype, ai_phase) = p_diag_bcpo_a1
                massptr_aer(8, isize, itype, ai_phase) = p_diag_seas_a1
                massptr_aer(9, isize, itype, ai_phase) = p_diag_dst_a1
                massptr_aer(10,isize, itype, ai_phase) = p_diag_soas_a1
                
                isize = 2
                numptr_aer( isize, itype, ai_phase)    = p_diag_num_a2
                massptr_aer(1, isize, itype, ai_phase) = p_diag_so4_a2
                massptr_aer(2, isize, itype, ai_phase) = p_diag_nit_a2
                massptr_aer(3, isize, itype, ai_phase) = p_diag_nh4_a2
                massptr_aer(4, isize, itype, ai_phase) = p_diag_ocpi_a2
                massptr_aer(5, isize, itype, ai_phase) = p_diag_ocpo_a2
                massptr_aer(6, isize, itype, ai_phase) = p_diag_bcpi_a2
                massptr_aer(7, isize, itype, ai_phase) = p_diag_bcpo_a2
                massptr_aer(8, isize, itype, ai_phase) = p_diag_seas_a2
                massptr_aer(9, isize, itype, ai_phase) = p_diag_dst_a2
                massptr_aer(10,isize, itype, ai_phase) = p_diag_soas_a2

                isize = 3
                numptr_aer( isize, itype, ai_phase)    = p_diag_num_a3
                massptr_aer(1, isize, itype, ai_phase) = p_diag_so4_a3
                massptr_aer(2, isize, itype, ai_phase) = p_diag_nit_a3
                massptr_aer(3, isize, itype, ai_phase) = p_diag_nh4_a3
                massptr_aer(4, isize, itype, ai_phase) = p_diag_ocpi_a3
                massptr_aer(5, isize, itype, ai_phase) = p_diag_ocpo_a3
                massptr_aer(6, isize, itype, ai_phase) = p_diag_bcpi_a3
                massptr_aer(7, isize, itype, ai_phase) = p_diag_bcpo_a3
                massptr_aer(8, isize, itype, ai_phase) = p_diag_seas_a3
                massptr_aer(9, isize, itype, ai_phase) = p_diag_dst_a3
                massptr_aer(10,isize, itype, ai_phase) = p_diag_soas_a3
                
                isize = 4
                numptr_aer( isize, itype, ai_phase)    = p_diag_num_a4
                massptr_aer(1, isize, itype, ai_phase) = p_diag_so4_a4
                massptr_aer(2, isize, itype, ai_phase) = p_diag_nit_a4
                massptr_aer(3, isize, itype, ai_phase) = p_diag_nh4_a4
                massptr_aer(4, isize, itype, ai_phase) = p_diag_ocpi_a4
                massptr_aer(5, isize, itype, ai_phase) = p_diag_ocpo_a4
                massptr_aer(6, isize, itype, ai_phase) = p_diag_bcpi_a4
                massptr_aer(7, isize, itype, ai_phase) = p_diag_bcpo_a4
                massptr_aer(8, isize, itype, ai_phase) = p_diag_seas_a4
                massptr_aer(9, isize, itype, ai_phase) = p_diag_dst_a4
                massptr_aer(10,isize, itype, ai_phase) = p_diag_soas_a4

                isize = 1
                numptr_aer( isize, itype, cw_phase)    = p_diag_num_cw1
                massptr_aer(1, isize, itype, cw_phase) = p_diag_so4_cw1
                massptr_aer(2, isize, itype, cw_phase) = p_diag_nit_cw1
                massptr_aer(3, isize, itype, cw_phase) = p_diag_nh4_cw1
                massptr_aer(4, isize, itype, cw_phase) = p_diag_ocpi_cw1
                massptr_aer(5, isize, itype, cw_phase) = p_diag_ocpo_cw1
                massptr_aer(6, isize, itype, cw_phase) = p_diag_bcpi_cw1
                massptr_aer(7, isize, itype, cw_phase) = p_diag_bcpo_cw1
                massptr_aer(8, isize, itype, cw_phase) = p_diag_seas_cw1
                massptr_aer(9, isize, itype, cw_phase) = p_diag_dst_cw1
                massptr_aer(10,isize, itype, cw_phase) = p_diag_soas_cw1

                isize = 2
                numptr_aer( isize, itype, cw_phase)    = p_diag_num_cw2
                massptr_aer(1, isize, itype, cw_phase) = p_diag_so4_cw2
                massptr_aer(2, isize, itype, cw_phase) = p_diag_nit_cw2
                massptr_aer(3, isize, itype, cw_phase) = p_diag_nh4_cw2
                massptr_aer(4, isize, itype, cw_phase) = p_diag_ocpi_cw2
                massptr_aer(5, isize, itype, cw_phase) = p_diag_ocpo_cw2
                massptr_aer(6, isize, itype, cw_phase) = p_diag_bcpi_cw2
                massptr_aer(7, isize, itype, cw_phase) = p_diag_bcpo_cw2
                massptr_aer(8, isize, itype, cw_phase) = p_diag_seas_cw2
                massptr_aer(9, isize, itype, cw_phase) = p_diag_dst_cw2
                massptr_aer(10,isize, itype, cw_phase) = p_diag_soas_cw2

                isize = 3
                numptr_aer( isize, itype, cw_phase)    = p_diag_num_cw3
                massptr_aer(1, isize, itype, cw_phase) = p_diag_so4_cw3
                massptr_aer(2, isize, itype, cw_phase) = p_diag_nit_cw3
                massptr_aer(3, isize, itype, cw_phase) = p_diag_nh4_cw3
                massptr_aer(4, isize, itype, cw_phase) = p_diag_ocpi_cw3
                massptr_aer(5, isize, itype, cw_phase) = p_diag_ocpo_cw3
                massptr_aer(6, isize, itype, cw_phase) = p_diag_bcpi_cw3
                massptr_aer(7, isize, itype, cw_phase) = p_diag_bcpo_cw3
                massptr_aer(8, isize, itype, cw_phase) = p_diag_seas_cw3
                massptr_aer(9, isize, itype, cw_phase) = p_diag_dst_cw3
                massptr_aer(10,isize, itype, cw_phase) = p_diag_soas_cw3

                isize = 4
                numptr_aer( isize, itype, cw_phase)    = p_diag_num_cw4
                massptr_aer(1, isize, itype, cw_phase) = p_diag_so4_cw4
                massptr_aer(2, isize, itype, cw_phase) = p_diag_nit_cw4
                massptr_aer(3, isize, itype, cw_phase) = p_diag_nh4_cw4
                massptr_aer(4, isize, itype, cw_phase) = p_diag_ocpi_cw4
                massptr_aer(5, isize, itype, cw_phase) = p_diag_ocpo_cw4
                massptr_aer(6, isize, itype, cw_phase) = p_diag_bcpi_cw4
                massptr_aer(7, isize, itype, cw_phase) = p_diag_bcpo_cw4
                massptr_aer(8, isize, itype, cw_phase) = p_diag_seas_cw4
                massptr_aer(9, isize, itype, cw_phase) = p_diag_dst_cw4
                massptr_aer(10,isize, itype, cw_phase) = p_diag_soas_cw4
                
                ! No Aerosol water
                ! do itype = 1, ntype_aer
                !  do isize = 1, nsize_aer(n)
                !     waterptr_aer(isize, itype) = -1
                !  end do
                ! end do

                ! Aerosol water
                itype = 1
                waterptr_aer(1, itype) = p_diag_water_a1
                waterptr_aer(2, itype) = p_diag_water_a2
                waterptr_aer(3, itype) = p_diag_water_a3
                waterptr_aer(4, itype) = p_diag_water_a4

        end subroutine diag_aero_wrfgc_aercld_ptrs


    subroutine diag_aero_size_info(nbin_o, chem, num_chem, relhum, is_gc,    &
                                  ids, ide, jds, jde, kds, kde,             &
                                  ims, ime, jms, jme, kms, kme,             &
                                  its, ite, jts, jte, kts, kte              )
    ! This subroutine is used to diagnose the mass/num concentrations in each size bin.

        USE module_data_gigc_asect
        USE module_data_sorgam, only: dginin, dginia, dginic, sginin, sginia, sginic
        USE module_data_gocart_dust, only: ra_dust, rb_dust

        implicit none

        integer, intent(in)                :: nbin_o
        real, dimension(ims:ime, kms:kme, jms:jme), intent(in)                         :: relhum
	logical, intent(in)                :: is_gc
        integer, intent(in)                :: num_chem
	real, dimension(ims:ime, kms:kme, jms:jme, num_chem),intent(inout)             :: chem
        integer, intent(in)                :: ids, ide, jds, jde, kds, kde,  &
                                              ims, ime, jms, jme, kms, kme,  &
                                              its, ite, jts, jte, kts, kte
    ! local variables
        real(kind = 8)            :: dens_h2o, mole_dryair
        real(kind = 8)            :: mass_so4i, mass_so4j, mass_niti, mass_nitj,         &
                                     mass_nh4i, mass_nh4j, mass_ocpii, mass_ocpij,       &
                                     mass_ocpoi, mass_ocpoj, mass_bcpii, mass_bcpij,     &
                                     mass_bcpoi, mass_bcpoj, mass_salj, mass_salc,       &
                                     mass_soasi, mass_soasj
        real(kind = 8)            :: mass_so4_gc, mass_nit_gc, mass_nh4_gc, mass_ocpi_gc, mass_ocpo_gc, &
                                     mass_bcpi_gc, mass_bcpo_gc, mass_dusttmp, mass_soas_gc
        real(kind = 8)            :: vol_h2otmp, vol_wet_aer, volaer_tmp

        real(kind = 8)            :: conv2so4, conv2nit, conv2nh4, conv2ocbc,   &
                                     conv2dst, conv2seas, conv2soas
        real(kind = 8)            :: dlogoc, dhigoc
 
	real, dimension(1:nbin_o) :: xnum_secti, xnum_sectj, xnum_sectc, xnum_sect_sna_oc, xnum_sect_bc, xnum_sect_sala, xnum_sect_salc
	real, dimension(1:nbin_o) :: xmas_secti, xmas_sectj, xmas_sectc, xmas_sect_sna_oc, xmas_sect_bc, xmas_sect_sala, xmas_sect_salc
	real, dimension(1:nbin_o) :: xdia_cm
	real, parameter           :: FRAC2Aitken = 0.10 ! Fraction of modal mass in Aitken mode
	real, parameter           :: ndust1 = 4         ! Number of dust bins of GEOS-Chem
	real                      :: dustfrc_gigc4bin(ndust1, nbin_o) ! GEOS-Chem dust size distribution - mass fracs in MOSAIC 4-bins
        real                      :: dgnum_um, dlo_um, dhi_um, duma, sixpi, dlo, dhi, dxbin, relh_frc, dgnum_sna_oc, dgnum_bc, sigma_sna_oc, sigma_bc, &
                                     dgnum_sala, dgnum_salc, sigma_sala, sigma_salc
        real                      :: dlo_sectm(nbin_o), dhi_sectm(nbin_o)
	real                      :: pi
        parameter (pi = 3.141592653589)
        integer                   :: i, k, j, isize, itype, icomp, m, n
                
        sixpi    = 6.0/pi
        dlo_um   = 0.0390625  ! low diameter limit (um)
        dhi_um   = 10.0       ! high diameter limit (um)
        duma     = 1.0
                
        dlo      = dlo_um*1.0E-6  ! low diameter limit (m)
        dhi      = dhi_um*1.0E-6  ! high diameter limit (m)
        dxbin    = (log(dhi) - log(dlo))/nbin_o
                do n = 1, nbin_o
                   dlo_sectm(n) = exp( log(dlo) + dxbin * (n-1) )
                   dhi_sectm(n) = exp( log(dlo) + dxbin * n )
                   xdia_cm(n)   = 0.5 * (dlo_sectm(n) + dhi_sectm(n)) * 1.0E2 ! unit:cm
                end do
                
                ! Dust bin mass fraction
                    dustfrc_gigc4bin = 0.
                        do m = 1, ndust1
                            dlogoc = ra_dust(m) * 2.E-6 ! low diameter limit (m)
                            dhigoc = rb_dust(m) * 2.E-6 ! high diameter limit (m)
                            do n = 1, nbin_o
                                dustfrc_gigc4bin(m, n) = max(DBLE(0.),min(DBLE(log(dhi_sectm(n))),log(dhigoc))- &
                                             max(log(dlogoc),DBLE(log(dlo_sectm(n)))) )/(log(dhigoc)-log(dlogoc))
                            end do
                        end do

        ! water densities in g/cm^3
            dens_h2o  = 1.0
        ! dry air mole weight g/mol
            mole_dryair = 28.97

        ! units:
        ! * mass     - ug/kg (dry air)
        ! * number   - #/kg (dry air)
        ! * volume   - cm^3/cm^3
        ! * diameter - cm
        
            do j = jts, jte
            do k = kts, kte
            do i = its, ite
            mass_so4i   = 0.0
            mass_so4j   = 0.0
            mass_niti   = 0.0
            mass_nitj   = 0.0
            mass_nh4i   = 0.0
            mass_nh4j   = 0.0
            mass_ocpii  = 0.0
            mass_ocpij  = 0.0
            mass_ocpoi  = 0.0
            mass_ocpoj  = 0.0
            mass_bcpii  = 0.0
            mass_bcpij  = 0.0
            mass_bcpoi  = 0.0
            mass_bcpoj  = 0.0
            mass_salj   = 0.0
            mass_salc   = 0.0
            mass_soasi  = 0.0
            mass_soasj  = 0.0
            mass_so4_gc = 0.0
            mass_nit_gc = 0.0
            mass_nh4_gc = 0.0
            mass_ocpi_gc= 0.0
            mass_ocpo_gc= 0.0
            mass_bcpi_gc= 0.0
            mass_bcpo_gc= 0.0
            mass_soas_gc= 0.0

        ! convert ppmv sulfate to ug/kg
            conv2so4 = mw_so4_aer/mole_dryair * 1.0E3
        ! convert ppmv nitrate to ug/kg		
            conv2nit  = mw_nit_aer/mole_dryair * 1.0E3
        ! convert ppmv ammonium to ug/kg		
            conv2nh4  = mw_nh4_aer/mole_dryair * 1.0E3
        ! convert ppmv oc/bc to ug/kg		
            conv2ocbc = mw_oc_aer/mole_dryair * 1.0E3
        ! convert ppmv dust to ug/kg		
            conv2dst  = mw_dst_aer/mole_dryair * 1.0E3
        ! convert ppmv sea salt ug/kg		
            conv2seas = mw_seas_aer/mole_dryair * 1.0E3
        ! convert ppmv simple soa ug/kg
            conv2soas = mw_soas_aer/mole_dryair * 1.0E3

            if (is_gc) then
                mass_so4_gc  = chem(i, k, j, p_so4) * conv2so4
                mass_nit_gc  = chem(i, k, j, p_nit) * conv2nit
                mass_nh4_gc  = chem(i, k, j, p_nh4) * conv2nh4
                mass_ocpi_gc = chem(i, k, j, p_ocpi) * conv2ocbc
                mass_ocpo_gc = chem(i, k, j, p_ocpo) * conv2ocbc
                mass_bcpi_gc = chem(i, k, j, p_bcpi) * conv2ocbc
                mass_bcpo_gc = chem(i, k, j, p_bcpo) * conv2ocbc
                mass_soas_gc = chem(i, k, j, p_soas) * conv2soas

            else
           ! Put the part of aerosol mass into Aitken mode 
                mass_so4i  = FRAC2Aitken * chem(i, k, j, p_so4) * conv2so4
                mass_niti  = FRAC2Aitken * chem(i, k, j, p_nit) * conv2nit
                mass_nh4i  = FRAC2Aitken * chem(i, k, j, p_nh4) * conv2nh4
                mass_ocpii = FRAC2Aitken * chem(i, k, j, p_ocpi) * conv2ocbc
                mass_ocpoi = FRAC2Aitken * chem(i, k, j, p_ocpo) * conv2ocbc
                mass_bcpii = FRAC2Aitken * chem(i, k, j, p_bcpi) * conv2ocbc
                mass_bcpoi = FRAC2Aitken * chem(i, k, j, p_bcpo) * conv2ocbc
                mass_soasi = FRAC2Aitken * chem(i, k, j, p_soas) * conv2soas

    ! Put the rest of aerosol mass into Accumulation mode
                mass_so4j  = (1.0-FRAC2Aitken) * chem(i, k, j, p_so4) * conv2so4
                mass_nitj  = (1.0-FRAC2Aitken) * chem(i, k, j, p_nit) * conv2nit
                mass_nh4j  = (1.0-FRAC2Aitken) * chem(i, k, j, p_nh4) * conv2nh4
                mass_ocpij = (1.0-FRAC2Aitken) * chem(i, k, j, p_ocpi) * conv2ocbc
                mass_ocpoj = (1.0-FRAC2Aitken) * chem(i, k, j, p_ocpo) * conv2ocbc
                mass_bcpij = (1.0-FRAC2Aitken) * chem(i, k, j, p_bcpi) * conv2ocbc
                mass_bcpoj = (1.0-FRAC2Aitken) * chem(i, k, j, p_bcpo) * conv2ocbc
                mass_soasj = (1.0-FRAC2Aitken) * chem(i, k, j, p_soas) * conv2soas

            end if

    ! sea salt in accumulation and coarse mode
                mass_salj = chem(i, k, j, p_sala) * conv2seas ! SALA represents the accumulation mode sea salt aerosol
                mass_salc = chem(i, k, j, p_salc) * conv2seas ! SALC represents the coarse mode sea salt aerosol

        ! call sect02_new to calculate the mass and number for each section
        ! sect02_new expects input in um
        ! size distribution is different for gc species
        ! secti is for aitken mode
        ! sectj is for accum mode
        ! sectc is for coarse mode
           if (is_gc) then
              dgnum_sna_oc = 0.07*2 ! unit:um
              sigma_sna_oc = 1.6  
              call sect02_new(dgnum_sna_oc, sigma_sna_oc, duma, nbin_o, dlo_um, dhi_um, &
                              xnum_sect_sna_oc, xmas_sect_sna_oc)
              dgnum_bc     = 0.02*2 ! unit:um
              sigma_bc     = 1.6
              call sect02_new(dgnum_bc, sigma_bc, duma, nbin_o, dlo_um, dhi_um,   &
                              xnum_sect_bc, xmas_sect_bc)
              itype = 1
              icomp = 1 ! SO4 (ug/kg-dryair)
                        do isize = 1, nsize_aer(itype)
                           chem(i, k, j, massptr_aer(icomp, isize, itype, ai_phase)) = mass_so4_gc * xmas_sect_sna_oc(isize)
                        end do
                        icomp = 2 ! NIT (ug/kg-dryair)
                        do isize = 1, nsize_aer(itype)
                           chem(i, k, j, massptr_aer(icomp, isize, itype, ai_phase)) = mass_nit_gc * xmas_sect_sna_oc(isize)
                        end do
                        icomp = 3 ! NH4 (ug/kg-dryair)
                        do isize = 1, nsize_aer(itype)
                           chem(i, k, j, massptr_aer(icomp, isize, itype, ai_phase)) = mass_nh4_gc * xmas_sect_sna_oc(isize)
                        end do
                        icomp = 4 ! OCPI (ug/kg-dryair)
                        do isize = 1, nsize_aer(itype)
                           chem(i, k, j, massptr_aer(icomp, isize, itype, ai_phase)) = mass_ocpi_gc * xmas_sect_sna_oc(isize)
                        end do
                        icomp = 5 ! OCPO (ug/kg-dryair)
                        do isize = 1, nsize_aer(itype)
                           chem(i, k, j, massptr_aer(icomp, isize, itype, ai_phase)) = mass_ocpo_gc * xmas_sect_sna_oc(isize)
                        end do
                        icomp = 6 ! BCPI (ug/kg-dryair)
                        do isize = 1, nsize_aer(itype)
                           chem(i, k, j, massptr_aer(icomp, isize, itype, ai_phase)) = mass_bcpi_gc * xmas_sect_bc(isize)
                        end do
                        icomp = 7 ! BCPO (ug/kg-dryair)
                        do isize = 1, nsize_aer(itype)
                           chem(i, k, j, massptr_aer(icomp, isize, itype, ai_phase)) = mass_bcpo_gc * xmas_sect_bc(isize)
                        end do
                        icomp = 10 ! SOAS (ug/kg-dryair)
                        do isize = 1, nsize_aer(itype)
                           chem(i, k, j, massptr_aer(icomp, isize, itype, ai_phase)) = mass_soas_gc * xmas_sect_sna_oc(isize)
                        end do

                else
                dgnum_um = dginin * 1.E6 
                call sect02_new(dgnum_um, sginin, duma, nbin_o, dlo_um, dhi_um,     &
                                xnum_secti, xmas_secti)

                dgnum_um = dginia * 1.E6
                call sect02_new(dgnum_um, sginia, duma, nbin_o, dlo_um, dhi_um,     &
                            xnum_sectj, xmas_sectj)

                dgnum_um = dginic * 1.E6
                call sect02_new(dgnum_um, sginic, duma, nbin_o, dlo_um, dhi_um,     &
                                xnum_sectc, xmas_sectc)

                        itype = 1
                        icomp = 1 ! S04 (ug/kg-dryair)
                        do isize = 1, nsize_aer(itype)
                           chem(i, k, j, massptr_aer(icomp, isize, itype, ai_phase)) = mass_so4i * xmas_secti(isize) + mass_so4j * xmas_sectj(isize)
                        end do
                        icomp = 2 ! NIT (ug/kg-dryair)
                        do isize = 1, nsize_aer(itype)
                           chem(i, k, j, massptr_aer(icomp, isize, itype, ai_phase)) = mass_niti * xmas_secti(isize) + mass_nitj * xmas_sectj(isize)
                        end do
                        icomp = 3 ! NH4 (ug/kg-dryair)
                        do isize = 1, nsize_aer(itype)
                           chem(i, k, j, massptr_aer(icomp, isize, itype,ai_phase)) = mass_nh4i * xmas_secti(isize) + mass_nh4j * xmas_sectj(isize)
                        end do
                        icomp = 4 ! OCPI (ug/kg-dryair)
                        do isize = 1, nsize_aer(itype)
                           chem(i, k, j, massptr_aer(icomp, isize, itype, ai_phase)) = mass_ocpii * xmas_secti(isize) + mass_ocpij * xmas_sectj(isize)
                        end do
                        icomp = 5 ! OCPO (ug/kg-dryair)
                        do isize = 1, nsize_aer(itype)
                           chem(i, k, j, massptr_aer(icomp, isize, itype, ai_phase)) = mass_ocpoi * xmas_secti(isize) + mass_ocpoj * xmas_sectj(isize)
                        end do
                        icomp = 6 ! BCPI (ug/kg-dryair)
                        do isize = 1, nsize_aer(itype)
                           chem(i, k, j, massptr_aer(icomp, isize, itype, ai_phase)) = mass_bcpii * xmas_secti(isize) + mass_bcpij * xmas_sectj(isize)
                        end do
                        icomp = 7 ! BCPO (ug/kg-dryair)
                        do isize = 1, nsize_aer(itype)
                           chem(i, k, j, massptr_aer(icomp, isize, itype, ai_phase)) = mass_bcpoi * xmas_secti(isize) + mass_bcpoj * xmas_sectj(isize)
                        end do
                        icomp = 10 ! SOAS (ug/kg-dryair)
                        do isize = 1, nsize_aer(itype)
                           chem(i, k, j, massptr_aer(icomp, isize, itype, ai_phase)) = mass_soasi * xmas_secti(isize) + mass_soasj * xmas_sectj(isize)
                        end do

            end if

                ! Sea salt aerosol
                dgnum_sala = 0.09*2 ! unit:um
                sigma_sala = 1.5
                call sect02_new(dgnum_sala, sigma_sala, duma, nbin_o, dlo_um, dhi_um, &
                              xnum_sect_sala, xmas_sect_sala)
                dgnum_salc = 0.4*2 ! unit:um
                sigma_salc = 1.8
                call sect02_new(dgnum_salc, sigma_salc, duma, nbin_o, dlo_um, dhi_um, &
                              xnum_sect_salc, xmas_sect_salc)


                icomp = 8 ! SEA-SALT (ug/kg-dryair)
                do isize = 1, nsize_aer(itype)
                   chem(i, k, j, massptr_aer(icomp, isize, itype, ai_phase)) = mass_salj * xmas_sect_sala(isize) + mass_salc * xmas_sect_salc(isize)
                end do
                 
                ! Dust aerosol
                   icomp = 9 ! Dust (ug/kg-dryair) 
                   do isize = 1, nsize_aer(itype) 
                      n = 0
                      mass_dusttmp = 0.0
                      do m = p_dst1, p_dst4
                         n = n + 1
                         mass_dusttmp = mass_dusttmp + dustfrc_gigc4bin(n, isize)*chem(i, k, j, m)
                      end do
                      chem(i, k, j, massptr_aer(icomp, isize, itype, ai_phase)) = mass_dusttmp * conv2dst
                   end do

                ! Aerosol water (ug/kg-dryair) 
                    relh_frc = amin1(0.9, relhum(i, k, j)) ! Put in fractional relative humidity, max of 0.9 here
                    itype = 1
                    do isize = 1, nsize_aer(itype)
                       vol_h2otmp = 0.0
                       do icomp = 1, ncomp_aer(itype)
                          vol_h2otmp = vol_h2otmp + chem(i, k, j, massptr_aer(icomp, isize, itype, ai_phase)) * &                                                                                         hygro_aer(icomp, itype) / dens_aer(icomp, itype)
                       end do
                          vol_h2otmp = relh_frc * vol_h2otmp / (1. - relh_frc)
                          chem(i, k, j, waterptr_aer(isize, itype)) = vol_h2otmp * dens_h2o
                    end do
                
                ! Aerosol number concentrations in each size bin
                    do isize = 1, nsize_aer(itype)
                       volaer_tmp = 0.0
                       do icomp = 1, ncomp_aer(itype)
                          volaer_tmp = volaer_tmp + chem(i, k, j, massptr_aer(icomp, isize, itype, ai_phase)) / &
                                       dens_aer(icomp, itype)
                       end do
                          vol_wet_aer = volaer_tmp + chem(i, k, j, waterptr_aer(isize, itype))/dens_h2o
                          chem(i, k, j, numptr_aer(isize, itype, ai_phase)) = vol_wet_aer * 1.0E-6 * sixpi / (xdia_cm(isize)*xdia_cm(isize)*xdia_cm(isize))
                    end do

    end do ! end i loop
    end do ! end k loop
    end do ! end j loop


        end subroutine diag_aero_size_info

  
    subroutine sect02_new(dgnum_um, sigmag, duma, nbin, dlo_um, dhi_um, &
                          xnum_sect, xmas_sect)
  
    ! This subroutine is based on the subroutine sect02 in the WRF-Chem 
    ! chem/module_optical_averaging.F
    ! user specifies a single log-normal mode and a set of section boundaries
    ! prog calculates mass and number for each section
    
        implicit none
        real, dimension(nbin), intent(out) :: xnum_sect, xmas_sect
        integer                            :: n, nbin
        real                               :: dgnum, dgnum_um, dhi, dhi_um,  &
                                              dlo, dlo_um, duma, dumfrac,    &
                                              dx, sigmag, sumnum, summas,    &
                                              sx, sxroot2, thi, tlo, x0, x3, &
                                              xhi, xlo, xmtot, xntot
        real                               :: dlo_sect(nbin), dhi_sect(nbin)  
        real                               :: pi
        parameter (pi = 3.141592653589)

            xmtot = duma
            xntot = duma

            dlo = dlo_um*1.0E-4
            dhi = dhi_um*1.0E-4
            xlo = log( dlo )
            xhi = log( dhi )
            dx  = (xhi - xlo)/nbin
        do n = 1, nbin
            dlo_sect(n) = exp( xlo + dx*(n-1) )
            dhi_sect(n) = exp( xlo + dx*n )
        end do

        dgnum = dgnum_um*1.0E-4
        sx = log( sigmag )
        x0 = log( dgnum )
        x3 = x0 + 3.*sx*sx
        sxroot2 = sx * sqrt( 2.0 )
        sumnum = 0.
        summas = 0.
            do n = 1, nbin
            xlo = log( dlo_sect(n) )
            xhi = log( dhi_sect(n) )
            tlo = (xlo - x0)/sxroot2
            thi = (xhi - x0)/sxroot2
            if (tlo .le. 0.) then
                dumfrac = 0.5*( erfc_num_recipes(-thi) - erfc_num_recipes(-tlo) )
            else
                dumfrac = 0.5*( erfc_num_recipes(tlo) - erfc_num_recipes(thi) )
            end if
            xnum_sect(n) = xntot*dumfrac
            tlo = (xlo - x3)/sxroot2
            thi = (xhi - x3)/sxroot2
            if (tlo .le. 0.) then
                dumfrac = 0.5*( erfc_num_recipes(-thi) - erfc_num_recipes(-tlo) )
            else
                dumfrac = 0.5*( erfc_num_recipes(tlo) - erfc_num_recipes(thi) )
            end if
            xmas_sect(n) = xmtot*dumfrac
            sumnum = sumnum + xnum_sect(n)
            summas = summas + xmas_sect(n)
        end do
    end subroutine sect02_new

!-----------------------------------------------------------------------
    real function erfc_num_recipes( x )
!
!   from press et al, numerical recipes, 1990, page 164
!
        implicit none
        real x
        double precision erfc_dbl, dum, t, z
        z = abs(x)
        t = 1.0/(1.0 + 0.5*z)
        dum =  ( -z*z - 1.26551223 + t*(1.00002368 + t*(0.37409196 +   &
          t*(0.09678418 + t*(-0.18628806 + t*(0.27886807 +   &
                                           t*(-1.13520398 +   &
          t*(1.48851587 + t*(-0.82215223 + t*0.17087277 )))))))))
        erfc_dbl = t * exp(dum)
        if (x .lt. 0.0) erfc_dbl = 2.0d0 - erfc_dbl
        erfc_num_recipes = erfc_dbl
        return

    end function erfc_num_recipes

  end module module_diag_aero_size_info
