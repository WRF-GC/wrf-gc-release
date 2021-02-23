!**********************************************************************************  
! This computer software was prepared by Battelle Memorial Institute, hereinafter
! the Contractor, under Contract No. DE-AC05-76RL0 1830 with the Department of 
! Energy (DOE). NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY,
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! MOSAIC module: see module_mosaic_driver.F for references and terms of use
!**********************************************************************************  
!
! WRF-chem V3.0 : Original version of optical_driver written by Jerome Fast (PNNL)
!                 and James Barnard (PNNL)
!
!WRF:MODEL_LAYER:CHEMISTRY
!
      SUBROUTINE optical_driver(id,curr_secs,dtstep,config_flags,haveaer,&
               is_gc, chem,dz8w,alt,relhum,                              &
               tauaer1,tauaer2,tauaer3,tauaer4,                          &
               !czhao
               extaer1,extaer2,extaer3,extaer4,                          &
               gaer1,gaer2,gaer3,gaer4,                                  &
               waer1,waer2,waer3,waer4,                                  &
               bscoef1,bscoef2,bscoef3,bscoef4,                          &
               l2aer,l3aer,l4aer,l5aer,l6aer,l7aer,                      & 
               extaerlw1,extaerlw2,extaerlw3,extaerlw4,extaerlw5,extaerlw6, &
               extaerlw7,extaerlw8,extaerlw9,extaerlw10,extaerlw11,extaerlw12, &
               extaerlw13,extaerlw14,extaerlw15,extaerlw16,  &
               tauaerlw1,tauaerlw2,tauaerlw3,tauaerlw4,tauaerlw5,tauaerlw6, & 
               tauaerlw7,tauaerlw8,tauaerlw9,tauaerlw10,tauaerlw11,tauaerlw12, & 
               tauaerlw13,tauaerlw14,tauaerlw15,tauaerlw16,  & 
               ids,ide, jds,jde, kds,kde,                                &
               ims,ime, jms,jme, kms,kme,                                &
               its,ite, jts,jte, kts,kte                                 )

!------------------------------------------------------------------------
   USE module_configure
   USE module_state_description
   USE module_model_constants
   USE module_optical_averaging
   USE module_data_rrtmgaeropt, only: nswbands,nlwbands
   use infnan,                 only: inf

   IMPLICIT NONE
   INTEGER,      INTENT(IN   ) :: id,                                  &
                                  ids,ide, jds,jde, kds,kde,           &
                                  ims,ime, jms,jme, kms,kme,           &
                                  its,ite, jts,jte, kts,kte
   REAL(KIND=8), INTENT(IN   ) :: curr_secs
   REAL,         INTENT(IN   ) :: dtstep
!
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),             &
         INTENT(IN ) ::  chem
!
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                       &  
         INTENT(IN ) ::  relhum, dz8w, alt
   LOGICAL, INTENT(IN ) :: is_gc
!
! arrays that hold the aerosol optical properties
!
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                       &  
         INTENT(INOUT ) ::                                             &
           tauaer1, tauaer2, tauaer3, tauaer4,                         &
           !czhao
           extaer1, extaer2, extaer3, extaer4,                         &
           gaer1, gaer2, gaer3, gaer4,                                 &
           waer1, waer2, waer3, waer4,                                 &
           bscoef1, bscoef2, bscoef3, bscoef4                              
   !for rrtmg shortwave and longwave --czhao
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme),                &
         INTENT(INOUT ) :: extaerlw1,extaerlw2,extaerlw3,extaerlw4,extaerlw5, & 
                           extaerlw6,extaerlw7,extaerlw8,extaerlw9,extaerlw10, & 
                           extaerlw11,extaerlw12,extaerlw13,extaerlw14,extaerlw15, &
                           extaerlw16,   & 
                           tauaerlw1,tauaerlw2,tauaerlw3,tauaerlw4,tauaerlw5, &
                           tauaerlw6,tauaerlw7,tauaerlw8,tauaerlw9,tauaerlw10, & 
                           tauaerlw11,tauaerlw12,tauaerlw13,tauaerlw14,tauaerlw15, &
                           tauaerlw16
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme,1:4) ::  & 
         tauaersw,extaersw,gaersw,waersw,bscoefsw 
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, 1:4 ),                  &  
         INTENT(INOUT ) ::                                             &
           l2aer, l3aer, l4aer, l5aer, l6aer, l7aer

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme,1:nlwbands) ::  &
         extaerlw,tauaerlw 

   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   character*100 msg
   integer lunerr 

   LOGICAL, INTENT(IN) :: haveaer
!         
! local variables
!
      logical processingAerosols
      integer nbin_o
      integer option_method, option_mie

!-----------------------------------------------------------------
! compute only if simulating aerosols and aer_ra_feedback=1

!  IF (config_flags%aer_ra_feedback .eq. 0) THEN
!        call wrf_debug(15,'no feedback, return from optical driver')
!    return
!  ENDIF
   select case (config_flags%chem_opt)
   case (233)
      processingAerosols = .true.
      call wrf_debug(15,'optical driver: process aerosols true')
   case default
      processingAerosols = .false.
      call wrf_debug(15,'optical driver: process aerosols false')
   end select

  if( processingAerosols ) then
!
! select aerosol optical property option
! VOLUME: volume averaging of refractive indicies
! * for MADE/SORGAM, assume same 8 size bins as MOSAIC by default
! SHELL: shell-core approach, placeholder
!
   select case (config_flags%chem_opt)
   case (233)
     nbin_o = 8
   end select
!
     call wrf_debug(15,'optical averaging')
     aer_op_opt_select: SELECT CASE(config_flags%aer_op_opt)
     CASE (VOLUME_APPROX)
       option_method=1
       option_mie=1
     CASE (MAXWELL_APPROX)
       option_method=2
       option_mie=1
     CASE (VOLUME_EXACT)
       option_method=1
       option_mie=2
     CASE (MAXWELL_EXACT)
       option_method=2
       option_mie=2
     CASE (SHELL_EXACT)
       option_method=3
       option_mie=2
     CASE DEFAULT
        if( config_flags%aer_op_opt > 0 ) then
           call wrf_message('WARNING: Invalid aer_op_opt. Defaulting to VOLUME_APPROX.')
           option_method=1
           option_mie=1
        end if
     END SELECT aer_op_opt_select

     if( config_flags%aer_op_opt > 0 ) then
        call wrf_debug(15,'optical driver: call optical averaging')
!        lunerr=-1

        !BSINGH(PNNL)- The followingvariables had an undefined behavior at the boundaries.
        tauaersw(:,:,:,:) = inf 
        extaersw(:,:,:,:) = inf
        gaersw(:,:,:,:)   = inf
        waersw(:,:,:,:)   = inf
        bscoefsw(:,:,:,:) = inf
        !long wave
        extaerlw(:,:,:,:) = inf
        tauaerlw(:,:,:,:) = inf

        call optical_averaging(id,curr_secs,dtstep,config_flags,     &
             nbin_o,haveaer,option_method,option_mie,chem,dz8w,alt,  &
             relhum,is_gc,                                           &
             tauaersw,extaersw,gaersw,waersw,bscoefsw,               &
             l2aer,l3aer,l4aer,l5aer,l6aer,l7aer,                    &
             tauaerlw,extaerlw,                                      &
             ids,ide, jds,jde, kds,kde,                              &
             ims,ime, jms,jme, kms,kme,                              &
             its,ite, jts,jte, kts,kte                               )

             !short wave
             tauaer1=tauaersw(:,:,:,1)
             tauaer2=tauaersw(:,:,:,2)
             tauaer3=tauaersw(:,:,:,3)
             tauaer4=tauaersw(:,:,:,4)
             extaer1=extaersw(:,:,:,1)
             extaer2=extaersw(:,:,:,2)
             extaer3=extaersw(:,:,:,3)
             extaer4=extaersw(:,:,:,4)
             gaer1=gaersw(:,:,:,1)
             gaer2=gaersw(:,:,:,2)
             gaer3=gaersw(:,:,:,3)
             gaer4=gaersw(:,:,:,4)
             waer1=waersw(:,:,:,1)
             waer2=waersw(:,:,:,2)
             waer3=waersw(:,:,:,3)
             waer4=waersw(:,:,:,4)
             bscoef1=bscoefsw(:,:,:,1)
             bscoef2=bscoefsw(:,:,:,2)
             bscoef3=bscoefsw(:,:,:,3)
             bscoef4=bscoefsw(:,:,:,4)
             !long wave
             extaerlw1=extaerlw(:,:,:,1)
             extaerlw2=extaerlw(:,:,:,2)
             extaerlw3=extaerlw(:,:,:,3)
             extaerlw4=extaerlw(:,:,:,4)
             extaerlw5=extaerlw(:,:,:,5)
             extaerlw6=extaerlw(:,:,:,6)
             extaerlw7=extaerlw(:,:,:,7)
             extaerlw8=extaerlw(:,:,:,8)
             extaerlw9=extaerlw(:,:,:,9)
             extaerlw10=extaerlw(:,:,:,10)
             extaerlw11=extaerlw(:,:,:,11)
             extaerlw12=extaerlw(:,:,:,12)
             extaerlw13=extaerlw(:,:,:,13)
             extaerlw14=extaerlw(:,:,:,14)
             extaerlw15=extaerlw(:,:,:,15)
             extaerlw16=extaerlw(:,:,:,16)
             tauaerlw1=tauaerlw(:,:,:,1)
             tauaerlw2=tauaerlw(:,:,:,2)
             tauaerlw3=tauaerlw(:,:,:,3)
             tauaerlw4=tauaerlw(:,:,:,4)
             tauaerlw5=tauaerlw(:,:,:,5)
             tauaerlw6=tauaerlw(:,:,:,6)
             tauaerlw7=tauaerlw(:,:,:,7)
             tauaerlw8=tauaerlw(:,:,:,8)
             tauaerlw9=tauaerlw(:,:,:,9)
             tauaerlw10=tauaerlw(:,:,:,10)
             tauaerlw11=tauaerlw(:,:,:,11)
             tauaerlw12=tauaerlw(:,:,:,12)
             tauaerlw13=tauaerlw(:,:,:,13)
             tauaerlw14=tauaerlw(:,:,:,14)
             tauaerlw15=tauaerlw(:,:,:,15)
             tauaerlw16=tauaerlw(:,:,:,16)
        call wrf_debug(15,'optical driver: after call optical averaging')
     else
        !If aer_op_opt==0 then the optical arrays are already set to
        !zero in chemics_init so there will not be a problem if the
        !user has selected aer_ra_feedback=1.
     end if
!
   endif
   return

END SUBROUTINE optical_driver
