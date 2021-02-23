module module_data_gigc_asect

implicit none

!-----------------------------------------------------------------------
!
!   The variables in this module provide a means of organizing and accessing
!   aerosol species in the "chem" array by their chemical component, 
!   size bin (or mode), "type", and "phase"
!
!   Their purpose is to allow flexible coding of process modules, 
!   compared to "hard-coding" using the chem array p_xxx indices
!   (e.g., p_diag_so4_a1, p_diag_so4_a2, ...; p_diag_num_a1, ...)
!
!-----------------------------------------------------------------------
!
!   maxd_atype = maximum allowable number of aerosol types
!   maxd_asize = maximum allowable number of aerosol size bins
!   maxd_acomp = maximum allowable number of chemical components
!       in each aerosol size bin
!   maxd_aphase = maximum allowable number of aerosol phases 
!       (gas, cloud, ice, rain, ...)
!
!   ntype_aer = number of aerosol types
!       The current GEOS-Chem aerosol scheme has only 1 type, with the 
!       implicit assumption of internal mixing.
!
!   nphase_aer = number of aerosol phases
!
!   ai_phase = phase (p) index for interstitial (unactivated) aerosol particles
!   cw_phase = phase (p) index for aerosol particles in cloud water
!   ci_phase = phase (p) index for aerosol particles in cloud ice
!   rn_phase = phase (p) index for aerosol particles in rain
!   sn_phase = phase (p) index for aerosol particles in snow
!   gr_phase = phase (p) index for aerosol particles in graupel
!   [Note:  the value of "xx_phase" will be between 1 and nphase_aer 
!       for phases that are active in a simulation.  The others
!       will have non-positive values.]
!   nsize_aer(t) = number of aerosol size bins for aerosol type t
!
!   ncomp_aer(t) = number of "regular" chemical components for aerosol type t
!   massptr_aer(c,s,t,p) = the position/index in the chem array for mixing- 
!       ratio for chemical component c, size bin s, type t, and phase p.
!
!   lptr_so4_aer(s,t,p) = the position/index in the chem array for mixing-
!       ratio for sulfate for aerosol size bin s, type t, and phase p
!   (similar lptr's are defined for no3, cl, msa, co3, 
!       nh4, na, ca, oin, oc, bc, ...)
!   [Note:  the massptr_aer allow you to loop over all species of 
!       an aerosol type.  The lptr_so4_aer, etc., allow you to access
!       a specific chemical component.]
!   numptr_aer(s,t,p) = the position/index in the chem array for mixing- 
!       ratio of particle number for size bin s, type t, and phase p.
!   dens_aer(c,t) = dry density (g/cm^3) of aerosol chemical component 
!       c of type t
!   mw_aer(c,t) = molecular weight (g/mole) of aerosol chemical component 
!       c of type t
!   hygro_aer(c,t) = bulk hygroscopicity (--) at dilute conditions 
!       (RH near 100%) of aerosol chemical component c of type t
!   [For definition of bulk hygroscopicity, 
!       see Abdul-Razzak and Ghan, 2004, J Geophys Res, V105, p. 6837-6844.]
!   dlo_sect(s,t) = 1-particle diameter (cm) at lower boundary of section m
!   dhi_sect(s,t) = 1-particle diameter (cm) at upper boundary of section m
!-----------------------------------------------------------------------

        integer, parameter :: maxd_atype  = 1
        integer, parameter :: maxd_asize  = 4
        integer, parameter :: maxd_acomp  = 10
        integer, parameter :: maxd_aphase = 2

        integer, save :: ai_phase = -999888777
        integer, save :: cw_phase = -999888777
        integer, save :: msectional = 1

        integer, save :: ntype_aer = 1                ! number of types		
        integer, save :: nphase_aer = 2               ! number of phases
        integer, save :: nsize_aer( maxd_atype )      ! number of size bins
        integer, save :: ncomp_aer( maxd_atype )      ! number of chemical components
        integer, save :: massptr_aer( maxd_acomp, maxd_asize, maxd_atype, maxd_aphase ), & 
                                                      ! index for mixing ratio
                         numptr_aer( maxd_asize, maxd_atype, maxd_aphase ), & ! index for the number mixing ratio                                      
                         waterptr_aer( maxd_asize, maxd_atype )  ! index for aerosol water
 
		real, save    :: dens_aer( maxd_acomp, maxd_atype),   &  ! unit: g/cm^3
                                 mw_aer( maxd_acomp, maxd_atype),     &  ! unit: g/mol
                                 hygro_aer( maxd_acomp, maxd_atype),  &  ! unitless
                                 dcen_sect( maxd_asize, maxd_atype),  &  ! unit: cm
                                 dlo_sect( maxd_asize, maxd_atype),   &  ! unit: cm
                                 dhi_sect( maxd_asize, maxd_atype),   &  ! unit: cm
                                 sigmag_aer(maxd_asize, maxd_atype)
 
        ! dry densities (g/cm^3)
        real, parameter :: dens_so4_aer_gc  = 1.7
        real, parameter :: dens_nit_aer_gc  = 1.8
        real, parameter :: dens_nh4_aer_gc  = 1.8   
        real, parameter :: dens_oc_aer_gc   = 1.3
        real, parameter :: dens_bc_aer_gc   = 1.8
	real, parameter :: dens_seas_aer_gc = 2.2
	real, parameter :: dens_dst1_aer_gc = 2.5
	real, parameter :: dens_dst2_aer_gc = 2.65
	real, parameter :: dens_dst3_aer_gc = 2.65
	real, parameter :: dens_dst4_aer_gc = 2.65
        real, parameter :: dens_soas_aer_gc = 1.5 ! dens based on Murphy and Pandis, ES&T, 2009
	real, parameter :: dens_water_aer_gc= 1.0

        ! molecular weight (g/mol)
	real, parameter :: mw_so4_aer       = 96.0
        real, parameter :: mw_nit_aer       = 62.0
	real, parameter :: mw_nh4_aer       = 18.0
	real, parameter :: mw_oc_aer        = 12.0
	real, parameter :: mw_bc_aer        = 12.0
	real, parameter :: mw_seas_aer      = 31.4
	real, parameter :: mw_dst_aer       = 29.0
        real, parameter :: mw_soas_aer      = 150.0
	real, parameter :: mw_water_aer     = 18.0

        ! hygroscopicity (dimensionless)
	real, parameter :: hygro_so4_aer    = 0.5
	real, parameter :: hygro_nit_aer    = 0.5
	real, parameter :: hygro_nh4_aer    = 0.5
	real, parameter :: hygro_oc_aer     = 0.2
	real, parameter :: hygro_bc_aer     = 1.E-6
	real, parameter :: hygro_seas_aer   = 1.16
	real, parameter :: hygro_dst_aer    = 0.14
        real, parameter :: hygro_soas_aer   = 0.14

end module module_data_gigc_asect
