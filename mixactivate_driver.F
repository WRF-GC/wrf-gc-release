!**********************************************************************************  
!
      SUBROUTINE mixactivate_driver(id, ktau, dtstep, config_flags,         &
                                    rho, t_phy, w, cldfra, cldfra_old,  &
                                    z, dz8w, p8w, t8w, exch_h, moist,       &
                                    scalar, chem, ccn1, ccn2,               &
                                    ccn3, ccn4, ccn5, ccn6, nsource,        & 
                                    ids,ide, jds,jde, kds,kde,              &
                                    ims,ime, jms,jme, kms,kme,              &
                                    its,ite, jts,jte, kts,kte               )

!------------------------------------------------------------------------
   USE module_configure
   USE module_state_description
   USE module_model_constants
   USE module_domain_type, only : domain
   USE module_mixactivate_wrappers, only: wrfgc_mixactivate

   IMPLICIT NONE
   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags
   INTEGER,      INTENT(IN   ) :: id,                                  &
                                  ids,ide, jds,jde, kds,kde,           &
                                  ims,ime, jms,jme, kms,kme,           &
                                  its,ite, jts,jte, kts,kte
   INTEGER,      INTENT(IN   ) :: ktau
   REAL,         INTENT(IN   ) :: dtstep
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_moist), INTENT(IN ) :: moist
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_scalar), INTENT(INOUT) :: scalar
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ), INTENT(INOUT ) ::  chem
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN ) :: rho, t_phy, dz8w,  &
                                                                p8w, t8w, exch_h, w,   &
                                                                z
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(INOUT) :: cldfra, cldfra_old
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(OUT) :: ccn1, ccn2, ccn3,      &
                                                                ccn4, ccn5, ccn6,      &
                                                                nsource
   
                                                                
! LOCAL VARIABLES
   INTEGER                                      :: idrydep_onoff
   REAL, DIMENSION( its:ite, jts:jte, num_chem) :: ddvel
   REAL, DIMENSION( ims:ime, jms:jme, num_chem) :: qsrflx
  ! REAL, DIMENSION( ims:ime, kms:kme, jms:jme) :: rho

! arrays that hold the aerosol optical properties
   
  !  where( alt(its:ite, kts:kte, jts:jte) /= 0.)
  !       dryrho_phy(its:ite, kts:kte, jts:jte) = 1./alt(its:ite, kts:kte, jts:jte)
  !  elsewhere
  !       dryrho_phy(its:ite, kts:kte, jts:jte) = 0.
  !  end where
    
    ddvel(:, :, :) = 0.0
    idrydep_onoff = 0

    call wrfgc_mixactivate(id, ktau, dtstep, config_flags, idrydep_onoff,    &
                           rho, t_phy, w, cldfra, cldfra_old,         &
                           ddvel, z, dz8w, p8w, t8w, exch_h,                 &
                           moist(ims, kms, jms, P_QV),                       & 
                           moist(ims, kms, jms, P_QC),                       &
                           moist(ims, kms, jms, P_QI),                       &
                           scalar(ims, kms, jms,P_QNDROP),                   & 
                           f_qc, f_qi, chem,                                 &
                           ccn1, ccn2, ccn3, ccn4, ccn5, ccn6, nsource,      &
                           qsrflx,                                           &
                           ids,ide, jds,jde, kds,kde,                        &
                           ims,ime, jms,jme, kms,kme,                        &
                           its,ite, jts,jte, kts,kte                         )

           
END SUBROUTINE mixactivate_driver
