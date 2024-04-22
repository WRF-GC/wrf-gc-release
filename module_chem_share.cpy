module module_chem_share
   use module_state_description

#if ( WRF_CHEM == 1 )
   contains

   integer function get_last_gas(chem_opt)
      implicit none
      integer, intent(in) :: chem_opt

      ! determine the index of the last gas species, which depends
      ! upon the gas mechanism.
      !
      ! only supports WRF-GC.
      get_last_gas = 0
      if(chem_opt .ge. 1) then
         get_last_gas = p_xyle
      endif
   end function get_last_gas
#endif
end module module_chem_share
