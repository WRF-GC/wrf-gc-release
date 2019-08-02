!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               __          _______  ______       _____  _____                 !
!               \ \        / /  __ \|  ____|     / ____|/ ____|                !
!                \ \  /\  / /| |__) | |__ ______| |  __| |                     !
!                 \ \/  \/ / |  _  /|  __|______| | |_ | |                     !
!                  \  /\  /  | | \ \| |         | |__| | |____                 !
!                   \/  \/   |_|  \_\_|          \_____|\_____|                !
!                                                                              !
!----------------------- ALPHA VERSION, v0.9 (20190802) -----------------------!
!
! WRF-GC: GEOS-Chem High Performance-powered Chemistry Add-On for WRF Model
! Developed by Haipeng Lin <hplin@g.harvard.edu>, Xu Feng <fengx7@pku.edu.cn>
!    January 2018, Peking University, Dept of Atmospheric and Oceanic Sciences
!    Correspondence to: Tzung-May Fu <fuzm@sustech.edu.cn>
!
! ALPHA INFORMATION:
!    WRF-GC Alpha (version 0.9) is experimental. Please notify the authors of
!    any bugs, suggestions and feature requests through email or the GitHub
!    repository.
!
! COPYRIGHT STATEMENT:
!    Permission is hereby granted, free of charge, to any person obtaining a copy
!   of this software and associated documentation files (the "Software"), to 
!   use, copy, modify the Software, and to permit persons to whom the Software is
!   furnished to do so, subject to the following conditions:
!
!   - The above copyright notice and this permission notice shall be included in all
!   copies or substantial portions of the Software.
! 
!   - The Software, modified in part or in full may not be redistributed without
!   express permission from the copyright holder(s).
! 
!   Except as contained in this notice or in attribution, the name of the WRF-GC model
!   shall not be used as an endorsement for distributing modified copies of the
!   Software without prior written permission from the copyright holder(s).
! 
!   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
!   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
!   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
!   SOFTWARE.
! 
!  WRF and the GEOS-Chem model, GCHP are (c) their original authors.
!
!------------------------------------------------------------------------------- 
module GIGC_Convert_State_Mod

    use module_domain
    use module_dm
    use module_model_constants
    use module_state_description
    use module_chem_utilities, only: calc_zenithb, getpar, get_cloud_optical_depth
    use module_convection_prep

    use GIGC_Chunk_Mod
    use PRECISION_MOD
    use Species_Mod,    only: MISSING_INT
    use ERROR_MOD,      only: IT_IS_NAN
    use Input_Opt_Mod,  only: OptInput
    use State_Chm_Mod,  only: ChmState, Ind_
    use State_Met_Mod,  only: MetState
    use State_Diag_Mod, only: DgnState

    implicit none
    private

    public :: GIGC_IdxSetup
    public :: GIGC_Get_WRF
    public :: GIGC_Set_WRF
    public :: GIGC_Chunk_Cleanup
    integer :: gi_no, gi_o3, gi_pan, gi_co, gi_alk4, gi_isop, gi_hno3, gi_h2o2, gi_acet, gi_mek, gi_ald2, gi_rcho, gi_mvk, gi_macr, gi_npmn, &
               gi_ppn, gi_r4n2, gi_prpe, gi_c3h8, gi_ch2o, gi_c2h6, gi_n2o5, gi_hno4, gi_mp, gi_dms, gi_so2, gi_so4, gi_so4s, gi_msa, &
               gi_nh3, gi_nh4, gi_nit, gi_nits, gi_bcpi, gi_ocpi, gi_bcpo, gi_ocpo, gi_dst1, gi_dst2, gi_dst3, gi_dst4, gi_sala, gi_salc
    integer :: gi_br2, gi_br, gi_bro, gi_hobr, gi_hbr, gi_brno2, gi_brno3, gi_chbr3, gi_ch2br2, gi_ch3br, gi_mpn, gi_isopnd, gi_isopnb, gi_moba, &
               gi_propnn, gi_hac, gi_glyc, gi_mvkn, gi_macrn, gi_map, gi_no2, gi_no3, gi_hno2, gi_benz, gi_tolu, gi_xyle, gi_mtpa, gi_limo
    integer :: gi_mtpo, gi_soap, gi_soas, gi_eoh, gi_mgly, gi_glyx, gi_acta, gi_hpald, gi_dhdn, gi_ethln, gi_hcooh, gi_iepoxa, gi_iepoxb, &
               gi_iepoxd, gi_isn1, gi_ripa, gi_ripb, gi_ripd, gi_imae, gi_soaie, gi_soame, gi_soagx, gi_soamg, gi_lvoc, gi_lvocoa, gi_isn1og
    integer :: gi_isn1oa, gi_monits, gi_monitu, gi_honit, gi_ionita, gi_monita, gi_indiol, gi_ipmn, gi_hc187, gi_n2o, gi_ocs, gi_ch4, gi_brcl, &
               gi_hcl, gi_ccl4, gi_ch3cl, gi_ch3ccl3, gi_cfc113, gi_cfc114, gi_cfc115, gi_hcfc123, gi_hcfc141b, gi_hcfc142b, gi_cfc11, gi_cfc12, &
               gi_hcfc22, gi_h1211, gi_h1301, gi_h2402, gi_cl, gi_clo, gi_hocl, gi_clno3, gi_clno2, gi_cloo, gi_oclo, gi_cl2, gi_cl2o2, gi_h2o
    integer :: gi_brsala, gi_brsalc, gi_chcl3, gi_ch2cl2, gi_ch3i, gi_ch2i2, gi_ch2icl, gi_ch2ibr, gi_hoi, gi_i2, gi_ibr, gi_icl, gi_i, gi_io, &
               gi_hi, gi_oio, gi_ino, gi_iono, gi_iono2, gi_i2o2, gi_i2o3, gi_i2o4, gi_isala, gi_isalc, gi_aeri, gi_pfe, gi_co2, gi_lbro2h, &
               gi_lbro2n, gi_lisopoh, gi_lisopno3, gi_ltro2h, gi_ltro2n, gi_lxro2h, gi_lxro2n, gi_pyac, gi_so4h1, gi_so4h2, gi_pox, gi_lox
    integer :: gi_pco, gi_lco, gi_pso4, gi_lch4, gi_ph2o2, gi_dhdc, gi_pmnn, gi_tro2, gi_bro2, gi_n, gi_xro2, gi_imao3, gi_etp, gi_ra3p, gi_rb3p, &
               gi_rp, gi_pp, gi_prpn, gi_pip, gi_r4p, gi_dhpcarp, gi_iap, gi_hpc52o2, gi_vrp, gi_atooh, gi_dhmob, gi_isnp, gi_maop, gi_mrp, gi_inpn
    integer :: gi_mobaoo, gi_diboo, gi_h, gi_macrno2, gi_roh, gi_isnohoo, gi_mvkoo, gi_isnoob, gi_gaoo, gi_ch3choo, gi_iepoxoo, gi_mglyoo, gi_prn1, &
               gi_mgloo, gi_a3o2, gi_man2, gi_isnooa, gi_po2, gi_isopndo2, gi_b3o2, gi_macroo, gi_r4n1, gi_maopo2, gi_ato2, gi_rco3, gi_limo2, &
               gi_olnd, gi_olnn, gi_ko2, gi_isopnbo2, gi_pio2, gi_hc5oo, gi_nmao3, gi_vro2, gi_ch2oo, gi_hc5, gi_rio2, gi_eto2, gi_ino2, gi_r4o2
    integer :: gi_mro2, gi_mco3, gi_mo2, gi_o1d, gi_o, gi_ho2, gi_oh, gi_h2, gi_moh, gi_n2, gi_o2, gi_rcooh
contains
    subroutine GIGC_IdxSetup(am_I_Root)
        logical, intent(in) :: am_I_Root

        gi_no = IND_('NO')
        gi_o3 = IND_('O3')
        gi_pan = IND_('PAN')
        gi_co = IND_('CO')
        gi_alk4 = IND_('ALK4')
        gi_isop = IND_('ISOP')
        gi_hno3 = IND_('HNO3')
        gi_h2o2 = IND_('H2O2')
        gi_acet = IND_('ACET')
        gi_mek = IND_('MEK')
        gi_ald2 = IND_('ALD2')
        gi_rcho = IND_('RCHO')
        gi_mvk = IND_('MVK')
        gi_macr = IND_('MACR')
        gi_npmn = IND_('NPMN')
        gi_ppn = IND_('PPN')
        gi_r4n2 = IND_('R4N2')
        gi_prpe = IND_('PRPE')
        gi_c3h8 = IND_('C3H8')
        gi_ch2o = IND_('CH2O')
        gi_c2h6 = IND_('C2H6')
        gi_n2o5 = IND_('N2O5')
        gi_hno4 = IND_('HNO4')
        gi_mp = IND_('MP')
        gi_dms = IND_('DMS')
        gi_so2 = IND_('SO2')
        gi_so4 = IND_('SO4')
        gi_so4s = IND_('SO4S')
        gi_msa = IND_('MSA')
        gi_nh3 = IND_('NH3')
        gi_nh4 = IND_('NH4')
        gi_nit = IND_('NIT')
        gi_nits = IND_('NITS')
        gi_bcpi = IND_('BCPI')
        gi_ocpi = IND_('OCPI')
        gi_bcpo = IND_('BCPO')
        gi_ocpo = IND_('OCPO')
        gi_dst1 = IND_('DST1')
        gi_dst2 = IND_('DST2')
        gi_dst3 = IND_('DST3')
        gi_dst4 = IND_('DST4')
        gi_sala = IND_('SALA')
        gi_salc = IND_('SALC')
        gi_br2 = IND_('BR2')
        gi_br = IND_('BR')
        gi_bro = IND_('BRO')
        gi_hobr = IND_('HOBR')
        gi_hbr = IND_('HBR')
        gi_brno2 = IND_('BRNO2')
        gi_brno3 = IND_('BRNO3')
        gi_chbr3 = IND_('CHBR3')
        gi_ch2br2 = IND_('CH2BR2')
        gi_ch3br = IND_('CH3BR')
        gi_mpn = IND_('MPN')
        gi_isopnd = IND_('ISOPND')
        gi_isopnb = IND_('ISOPNB')
        gi_moba = IND_('MOBA')
        gi_propnn = IND_('PROPNN')
        gi_hac = IND_('HAC')
        gi_glyc = IND_('GLYC')
        gi_mvkn = IND_('MVKN')
        gi_macrn = IND_('MACRN')
        gi_map = IND_('MAP')
        gi_no2 = IND_('NO2')
        gi_no3 = IND_('NO3')
        gi_hno2 = IND_('HNO2')
        gi_benz = IND_('BENZ')
        gi_tolu = IND_('TOLU')
        gi_xyle = IND_('XYLE')
        gi_mtpa = IND_('MTPA')
        gi_limo = IND_('LIMO')
        gi_mtpo = IND_('MTPO')
        gi_soap = IND_('SOAP')
        gi_soas = IND_('SOAS')
        gi_eoh = IND_('EOH')
        gi_mgly = IND_('MGLY')
        gi_glyx = IND_('GLYX')
        gi_acta = IND_('ACTA')
        gi_hpald = IND_('HPALD')
        gi_dhdn = IND_('DHDN')
        gi_ethln = IND_('ETHLN')
        gi_hcooh = IND_('HCOOH')
        gi_iepoxa = IND_('IEPOXA')
        gi_iepoxb = IND_('IEPOXB')
        gi_iepoxd = IND_('IEPOXD')
        gi_isn1 = IND_('ISN1')
        gi_ripa = IND_('RIPA')
        gi_ripb = IND_('RIPB')
        gi_ripd = IND_('RIPD')
        gi_imae = IND_('IMAE')
        gi_soaie = IND_('SOAIE')
        gi_soame = IND_('SOAME')
        gi_soagx = IND_('SOAGX')
        gi_soamg = IND_('SOAMG')
        gi_lvoc = IND_('LVOC')
        gi_lvocoa = IND_('LVOCOA')
        gi_isn1og = IND_('ISN1OG')
        gi_isn1oa = IND_('ISN1OA')
        gi_monits = IND_('MONITS')
        gi_monitu = IND_('MONITU')
        gi_honit = IND_('HONIT')
        gi_ionita = IND_('IONITA')
        gi_monita = IND_('MONITA')
        gi_indiol = IND_('INDIOL')
        gi_ipmn = IND_('IPMN')
        gi_hc187 = IND_('HC187')
        gi_n2o = IND_('N2O')
        gi_ocs = IND_('OCS')
        gi_ch4 = IND_('CH4')
        gi_brcl = IND_('BRCL')
        gi_hcl = IND_('HCL')
        gi_ccl4 = IND_('CCL4')
        gi_ch3cl = IND_('CH3CL')
        gi_ch3ccl3 = IND_('CH3CCL3')
        gi_cfc113 = IND_('CFC113')
        gi_cfc114 = IND_('CFC114')
        gi_cfc115 = IND_('CFC115')
        gi_hcfc123 = IND_('HCFC123')
        gi_hcfc141b = IND_('HCFC141B')
        gi_hcfc142b = IND_('HCFC142B')
        gi_cfc11 = IND_('CFC11')
        gi_cfc12 = IND_('CFC12')
        gi_hcfc22 = IND_('HCFC22')
        gi_h1211 = IND_('H1211')
        gi_h1301 = IND_('H1301')
        gi_h2402 = IND_('H2402')
        gi_cl = IND_('CL')
        gi_clo = IND_('CLO')
        gi_hocl = IND_('HOCL')
        gi_clno3 = IND_('CLNO3')
        gi_clno2 = IND_('CLNO2')
        gi_cloo = IND_('CLOO')
        gi_oclo = IND_('OCLO')
        gi_cl2 = IND_('CL2')
        gi_cl2o2 = IND_('CL2O2')
        gi_h2o = IND_('H2O')
        gi_brsala = IND_('BRSALA')
        gi_brsalc = IND_('BRSALC')
        gi_chcl3 = IND_('CHCL3')
        gi_ch2cl2 = IND_('CH2CL2')
        gi_ch3i = IND_('CH3I')
        gi_ch2i2 = IND_('CH2I2')
        gi_ch2icl = IND_('CH2ICL')
        gi_ch2ibr = IND_('CH2IBR')
        gi_hoi = IND_('HOI')
        gi_i2 = IND_('I2')
        gi_ibr = IND_('IBR')
        gi_icl = IND_('ICL')
        gi_i = IND_('I')
        gi_io = IND_('IO')
        gi_hi = IND_('HI')
        gi_oio = IND_('OIO')
        gi_ino = IND_('INO')
        gi_iono = IND_('IONO')
        gi_iono2 = IND_('IONO2')
        gi_i2o2 = IND_('I2O2')
        gi_i2o3 = IND_('I2O3')
        gi_i2o4 = IND_('I2O4')
        gi_isala = IND_('ISALA')
        gi_isalc = IND_('ISALC')
        gi_aeri = IND_('AERI')
        gi_pfe = IND_('PFE')
        gi_co2 = IND_('CO2')
        gi_lbro2h = IND_('LBRO2H')
        gi_lbro2n = IND_('LBRO2N')
        gi_lisopoh = IND_('LISOPOH')
        gi_lisopno3 = IND_('LISOPNO3')
        gi_ltro2h = IND_('LTRO2H')
        gi_ltro2n = IND_('LTRO2N')
        gi_lxro2h = IND_('LXRO2H')
        gi_lxro2n = IND_('LXRO2N')
        gi_pyac = IND_('PYAC')
        gi_so4h1 = IND_('SO4H1')
        gi_so4h2 = IND_('SO4H2')
        gi_pox = IND_('POX')
        gi_lox = IND_('LOX')
        gi_pco = IND_('PCO')
        gi_lco = IND_('LCO')
        gi_pso4 = IND_('PSO4')
        gi_lch4 = IND_('LCH4')
        gi_ph2o2 = IND_('PH2O2')
        gi_dhdc = IND_('DHDC')
        gi_pmnn = IND_('PMNN')
        gi_tro2 = IND_('TRO2')
        gi_bro2 = IND_('BRO2')
        gi_n = IND_('N')
        gi_xro2 = IND_('XRO2')
        gi_imao3 = IND_('IMAO3')
        gi_etp = IND_('ETP')
        gi_ra3p = IND_('RA3P')
        gi_rb3p = IND_('RB3P')
        gi_rp = IND_('RP')
        gi_pp = IND_('PP')
        gi_prpn = IND_('PRPN')
        gi_pip = IND_('PIP')
        gi_r4p = IND_('R4P')
        gi_dhpcarp = IND_('DHPCARP')
        gi_iap = IND_('IAP')
        gi_hpc52o2 = IND_('HPC52O2')
        gi_vrp = IND_('VRP')
        gi_atooh = IND_('ATOOH')
        gi_dhmob = IND_('DHMOB')
        gi_isnp = IND_('ISNP')
        gi_maop = IND_('MAOP')
        gi_mrp = IND_('MRP')
        gi_inpn = IND_('INPN')
        gi_mobaoo = IND_('MOBAOO')
        gi_diboo = IND_('DIBOO')
        gi_h = IND_('H')
        gi_macrno2 = IND_('MACRNO2')
        gi_roh = IND_('ROH')
        gi_isnohoo = IND_('ISNOHOO')
        gi_mvkoo = IND_('MVKOO')
        gi_isnoob = IND_('ISNOOB')
        gi_gaoo = IND_('GAOO')
        gi_ch3choo = IND_('CH3CHOO')
        gi_iepoxoo = IND_('IEPOXOO')
        gi_mglyoo = IND_('MGLYOO')
        gi_prn1 = IND_('PRN1')
        gi_mgloo = IND_('MGLOO')
        gi_a3o2 = IND_('A3O2')
        gi_man2 = IND_('MAN2')
        gi_isnooa = IND_('ISNOOA')
        gi_po2 = IND_('PO2')
        gi_isopndo2 = IND_('ISOPNDO2')
        gi_b3o2 = IND_('B3O2')
        gi_macroo = IND_('MACROO')
        gi_r4n1 = IND_('R4N1')
        gi_maopo2 = IND_('MAOPO2')
        gi_ato2 = IND_('ATO2')
        gi_rco3 = IND_('RCO3')
        gi_limo2 = IND_('LIMO2')
        gi_olnd = IND_('OLND')
        gi_olnn = IND_('OLNN')
        gi_ko2 = IND_('KO2')
        gi_isopnbo2 = IND_('ISOPNBO2')
        gi_pio2 = IND_('PIO2')
        gi_hc5oo = IND_('HC5OO')
        gi_nmao3 = IND_('NMAO3')
        gi_vro2 = IND_('VRO2')
        gi_ch2oo = IND_('CH2OO')
        gi_hc5 = IND_('HC5')
        gi_rio2 = IND_('RIO2')
        gi_eto2 = IND_('ETO2')
        gi_ino2 = IND_('INO2')
        gi_r4o2 = IND_('R4O2')
        gi_mro2 = IND_('MRO2')
        gi_mco3 = IND_('MCO3')
        gi_mo2 = IND_('MO2')
        gi_o1d = IND_('O1D')
        gi_o = IND_('O')
        gi_ho2 = IND_('HO2')
        gi_oh = IND_('OH')
        gi_h2 = IND_('H2')
        gi_moh = IND_('MOH')
        gi_n2 = IND_('N2')
        gi_o2 = IND_('O2')
        gi_rcooh = IND_('RCOOH')
    end subroutine GIGC_IdxSetup
    subroutine GIGC_Get_WRF(am_I_Root, &
        config_flags, grid, &
        num_chem, chem, num_scalar, scalar, num_moist, moist, &
        dz8w, p8w, t8w, pi_phy, &
        f_qc, f_qi, f_qs, &
        its, ite, jts, jte, &
        ide, jde, &
        kts, kte, &
        ids, jds, kds, &
        curr_secs, &
        Input_Opt, State_Met, State_Chm, &
        WRFGC_Phys_Time)
        use PRESSURE_MOD,  only: ACCEPT_EXTERNAL_APBP

        implicit none

        logical, intent(in) :: am_I_Root
        type(grid_config_rec_type), intent(in) :: config_flags

        type(domain), target :: grid
        integer, intent(in)  :: num_chem, num_scalar, num_moist ! # of WRF chem, moist variables. Check "dummy_new_args.inc" & "actual_new_args.inc" in WRF compile.
        real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33, num_chem), intent(in)   :: chem
        real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33, num_scalar), intent(in) :: scalar
        real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33, num_moist), intent(in)  :: moist
        real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33), intent(in)  :: dz8w, p8w, pi_phy, t8w
        integer, intent(in) :: its, ite, jts, jte, ide, jde, kts, kte, ids, jds, kds

        real(KIND=8) :: curr_secs
        type(OptInput), intent(inout) :: Input_Opt
        type(MetState), intent(inout) :: State_Met
        type(ChmState), intent(inout) :: State_Chm

        real(KIND=8), intent(inout)   :: WRFGC_Phys_Time  ! Cumulative time spent in physics/parameterization calculation time

        integer :: debug_level      ! WRF Debug Level
        integer :: GEOS_CHEM_RC     ! GEOS-Chem Return Code (0 = Success)

        integer :: IM, II           ! # of lons on this PET
        integer :: JM, JJ           ! # of lats on this PET
        integer :: LM               ! # of levs on this PET
        integer :: i, j, k          ! Loop variable over xts, xte (x = i, j, k)
        integer :: N

        real :: zen, pardb, pardif  ! For computing %SUNCOSmid, %PARDB, %PARDF
        real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33) :: cmfmc, pmflxrain, pmflxsnow    ! For computing %CMFMC, %PFLCU, %PFICU
        real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33) :: dqrcu  ! For computing %DQRCU
        integer, dimension(grid%sm31:grid%em31, grid%sm33:grid%em33) :: cloud_bot ! For computing %DQRCU
        real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33) :: dtrainu, dtraind, dtrain ! For computing %DTRAIN
        real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33) :: reevapcn ! For computing %REEVAPCN
        logical, dimension(grid%sm31:grid%em31, grid%sm33:grid%em33) :: cu_act_flag ! For cu_ntiedtke
        real, dimension(grid%sm31:grid%em31, grid%sm33:grid%em33) :: raincv_diag, pratec_diag ! For cu_ntiedtke
        real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33) :: rqccuten_diag, rqicuten_diag, rqvcuten_diag, &
                                                                                          rthcuten_diag, rucuten_diag, rvcuten_diag ! For preventing re-writing the values of grid%
        real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33) :: taucldc, taucldi, optd !For computing %TAUCLW, %TAUCLI, %OPTD
        logical, optional, intent(in) :: f_qc, f_qi, f_qs ! For computing %TAUCLW, %TAUCLI, %OPTD        

        real(fp), allocatable :: Ap(:)
        real(fp), allocatable :: Bp(:)

        character*256 :: mminlu_loc

        real(KIND=8)          :: WRFGC_Time_Temp_Start, WRFGC_Time_Temp_End

        logical, save :: FIRST = .true.      ! Is first run?

        call nl_get_debug_level(1, debug_level)
        call set_wrf_debug_level(debug_level)

        call nl_get_mminlu(1, mminlu_loc) 

        IM = ite - its + 1
        JM = jte - jts + 1
        LM = kte - kts + 1


        WRFGC_Time_Temp_Start = MPI_Wtime()

        if(config_flags%cu_physics .eq. 16) then
        
            raincv_diag   = grid%raincv
            pratec_diag   = grid%pratec
            rqccuten_diag = grid%rqccuten
            rqicuten_diag = grid%rqicuten
            rqvcuten_diag = grid%rqvcuten
            rthcuten_diag = grid%rthcuten
            rucuten_diag  = grid%rucuten
            rvcuten_diag  = grid%rvcuten
        
            call cu_ntiedtke(                                         &
                 grid%dt, grid%itimestep, grid%stepcu                 &
                ,raincv_diag, pratec_diag, grid%qfx, grid%hfx         &
                ,grid%u_phy, grid%v_phy, grid%w_2, grid%t_phy         &
                ,moist(grid%sm31,grid%sm32,grid%sm33,p_qv)            &
                ,moist(grid%sm31,grid%sm32,grid%sm33,p_qc)            &
                ,moist(grid%sm31,grid%sm32,grid%sm33,p_qi)            &
                ,pi_phy, grid%rho                                     &
                ,grid%rqvften, grid%rthften                           &
                ,dz8w, grid%p_hyd, grid%p_hyd_w, grid%xland           & 
                ,cu_act_flag, grid%dx                                 &
                ,ids, ide, jds, jde, kds, kte+1                       &
                ,grid%sm31, grid%em31, grid%sm33, grid%em33           &
                ,grid%sm32, grid%em32                                 &
                ,its, ite, jts, jte, kts, kte                         &
                ,rthcuten_diag, rqvcuten_diag, rqccuten_diag          &
                ,rqicuten_diag, rucuten_diag, rvcuten_diag            &
                ,pmflxrain, pmflxsnow, cmfmc, dqrcu, cloud_bot        &
                ,dtrainu, dtraind, dtrain, reevapcn                   )
        endif
        
        if(config_flags%ra_lw_physics .eq. 4) then
           call get_cloud_optical_depth(t8w, dz8w                     & 
               ,grid%t_phy, grid%p_hyd_w                              &
               ,moist(grid%sm31,grid%sm32,grid%sm33,p_qc)             &
               ,moist(grid%sm31,grid%sm32,grid%sm33,p_qi)             &
               ,moist(grid%sm31,grid%sm32,grid%sm33,p_qs)             &
               ,taucldc, taucldi, optd                                &
               ,f_qc, f_qi, f_qs                                      &
               ,grid%sm31, grid%em31, grid%sm33, grid%em33            &
               ,grid%sm32, grid%em32                                  &
               ,its, ite, jts, jte, kts, kte                          )
        endif
        
        if(config_flags%hybrid_opt .eq. 2) then
            allocate(Ap(LM+1))
            allocate(Bp(LM+1))

            do k = kts, kte
                Ap(k) = (grid%c4f(k) + (1 - grid%c3f(k)) * grid%p_top) * .01_fp
                Bp(k) = grid%c3f(k)
            enddo

            Ap(LM+1) = grid%p_top * .01_fp
            Bp(LM+1) = 0


            call ACCEPT_EXTERNAL_APBP(am_I_Root, Ap, Bp, GEOS_CHEM_RC)

            deallocate(Ap)
            deallocate(Bp)
        else
            write(6, *) "====================================================="
            write(6, *) "|               W A R N I N G (WRF-GC)              |"
            write(6, *) "====================================================="
            write(6, *) "  GEOS-CHEM USES A VERTICAL HYBRID-SIGMA GRID.       "
            write(6, *) "  WRF MUST BE CONFIGURED TO USE THIS GRID using the  "
            write(6, *) "  namelist option &dynamics: hybrid_opt = 2, and     "
            write(6, *) "  RECOMPILING with ./configure -hyb.                 "
            write(6, *) "                                                     "
            write(6, *) "  WE WERE *NOT* ABLE TO DETECT THIS IN YOUR WRF-GC   "
            write(6, *) "  CONFIGURATION, WHICH MEANS THE VERTICAL LEVELS MAY "
            write(6, *) "  BE INACCURATE AND OUTRIGHT WRONG. PLEASE CHECK.    "
            write(6, *) "====================================================="
        endif

        WRFGC_Time_Temp_End = MPI_Wtime()
        WRFGC_Phys_Time = WRFGC_Phys_Time + (WRFGC_Time_Temp_End - WRFGC_Time_Temp_Start)

        call wrf_debug(100, "GIGC_Convert_State_Mod before State_Met% Conversion")

        do j = jts, jte
        do i = its, ite
            II = i - its + 1
            JJ = j - jts + 1

            State_Met%ALBD(II, JJ) = grid%albedo(i, j)

            State_Met%CLDFRC(II, JJ) = grid%CLDT(i, j)

            
            State_Met%EFLUX(II, JJ) = grid%LH(i, j)

            State_Met%FRCLND(II, JJ) = State_Met%FRLAND(II, JJ)

            if(grid%lu_mask(i, j) .eq. 0 .and. grid%lakemask(i, j) .eq. 0 .and. grid%snowh(i, j) .eq. 0) then
                 State_Met%FRLAND(II, JJ) = 1.0_fp
                 State_Met%FROCEAN(II, JJ) = 0.0_fp
                 State_Met%FRSNO(II, JJ) = 0.0_fp
                 State_Met%FRLAKE(II, JJ) = 0.0_fp
                 State_Met%FRLANDIC(II, JJ) = 0.0_fp
            else if(grid%lu_mask(i, j) .eq. 1) then
                 State_Met%FRLAND(II, JJ) = 0.0_fp
                 State_Met%FROCEAN(II, JJ) = 1.0_fp
                 State_Met%FRSNO(II, JJ) = 0.0_fp
                 State_Met%FRLAKE(II, JJ) = 0.0_fp
                 State_Met%FRLANDIC(II, JJ) = 0.0_fp
            else if(grid%lakemask(i, j) .eq. 1) then
                 State_Met%FRLAND(II, JJ) = 0.0_fp
                 State_Met%FROCEAN(II, JJ) = 0.0_fp
                 State_Met%FRSNO(II, JJ) = 0.0_fp
                 State_Met%FRLAKE(II, JJ) = 1.0_fp
                 State_Met%FRLANDIC(II, JJ) = 0.0_fp
            else if(grid%snowh(i, j) .ge. 0) then
                 State_Met%FRLAND(II, JJ) = 0.0_fp
                 State_Met%FROCEAN(II, JJ) = 0.0_fp
                 State_Met%FRSNO(II, JJ) = 1.0_fp
                 State_Met%FRLAKE(II, JJ) = 0.0_fp
                 State_Met%FRLANDIC(II, JJ) = 1.0_fp
            endif

            State_Met%FRSEAICE(II, JJ) = grid%icefrac_gc(i, j)

            if(mminlu_loc .eq. 'USGS') then
                select case(int(grid%lu_index(i, j)))
                case (1, 31, 32, 33) ! Urban, residential
                    State_Met%LandTypeFrac(II, JJ, 2 ) = 1.0_fp
                    State_Met%XLAI        (II, JJ, 2 ) = grid%lai(i, j)
                case (2, 3, 4, 5, 6) ! Dryland Cropland, Cr/Grassland Mosaic, Cr/Woodland Mosaic
                    State_Met%LandTypeFrac(II, JJ, 39) = 1.0_fp ! Cool Irrigated Cropland
                    State_Met%XLAI        (II, JJ, 39) = grid%lai(i, j)
                case (7, 8, 9) ! Grassland, Shrubland, Mixed Sh/Gr
                    State_Met%LandTypeFrac(II, JJ, 43) = 1.0_fp ! Cold Grassland
                    State_Met%XLAI        (II, JJ, 43) = grid%lai(i, j)
                case (10) ! Savanna
                    State_Met%LandTypeFrac(II, JJ, 44) = 1.0_fp ! Woody Savanna
                    State_Met%XLAI        (II, JJ, 44) = grid%lai(i, j)
                case (11) ! Deciduous Broadleaf Forest
                    State_Met%LandTypeFrac(II, JJ, 6 ) = 1.0_fp
                    State_Met%XLAI        (II, JJ, 6 ) = grid%lai(i, j)
                case (12) ! Deciduous Needleleaf Forest
                    State_Met%LandTypeFrac(II, JJ, 5 ) = 1.0_fp
                    State_Met%XLAI        (II, JJ, 5 ) = grid%lai(i, j)
                case (13) ! Evergreen Broadleaf Forest
                    State_Met%LandTypeFrac(II, JJ, 7 ) = 1.0_fp
                    State_Met%XLAI        (II, JJ, 7 ) = grid%lai(i, j)
                case (14) ! Evergreen Needleleaf Forest
                    State_Met%LandTypeFrac(II, JJ, 4 ) = 1.0_fp ! Coniferous Forest
                    State_Met%XLAI        (II, JJ, 4 ) = grid%lai(i, j)
                case (15) ! Mixed forest
                    State_Met%LandTypeFrac(II, JJ, 25) = 1.0_fp
                    State_Met%XLAI        (II, JJ, 25) = grid%lai(i, j)
                case (16) ! Water bodies
                    State_Met%LandTypeFrac(II, JJ, 1 ) = 1.0_fp ! Inland/Sea Water
                    State_Met%XLAI        (II, JJ, 1 ) = grid%lai(i, j)
                case (18) ! Wooded Wetland
                    State_Met%LandTypeFrac(II, JJ, 14) = 1.0_fp ! Wooded Wet Swamp
                    State_Met%XLAI        (II, JJ, 14) = grid%lai(i, j)
                case (19) ! Barren or Sparsely Vegetated
                    State_Met%LandTypeFrac(II, JJ, 3 ) = 1.0_fp ! Low Sparse Grassland
                    State_Met%XLAI        (II, JJ, 3 ) = grid%lai(i, j)
                case (21) ! Wooded Tundra
                    State_Met%LandTypeFrac(II, JJ, 64) = 1.0_fp ! Wooden Tundra
                    State_Met%XLAI        (II, JJ, 64) = grid%lai(i, j)
                case (23) ! Bare Ground Tundra
                    State_Met%LandTypeFrac(II, JJ, 54) = 1.0_fp ! Barren Tundra
                    State_Met%XLAI        (II, JJ, 54) = grid%lai(i, j)
                case (24) ! Snow or Ice
                    State_Met%LandTypeFrac(II, JJ, 13) = 1.0_fp ! Glacier Ice
                    State_Met%XLAI        (II, JJ, 13) = grid%lai(i, j)
                case (25) ! Playa
                    State_Met%LandTypeFrac(II, JJ, 8 ) = 1.0_fp ! Bare Desert
                    State_Met%XLAI        (II, JJ, 8 ) = grid%lai(i, j)
                case default
                    State_Met%LandTypeFrac(II, JJ, 2 ) = 1.0_fp
                    State_Met%XLAI        (II, JJ, 2 ) = grid%lai(i, j)
                end select
            else if(mminlu_loc .eq. 'MODIFIED_IGBP_MODIS_NOAH') then

            endif














            
            State_Met%GWETROOT(II, JJ) = grid%sm100200(i, j)
            
            State_Met%GWETTOP(II, JJ) = grid%sm000010(i, j)
            
            State_Met%HFLUX(II, JJ) = grid%hfx(i, j)

            State_Met%LAI(II, JJ) = grid%lai(i, j)

            State_Met%MODISLAI(II, JJ) = grid%lai(i, j)

            
            if(grid%lu_mask(i, j) .eq. 0) then
                State_Met%LWI(II, JJ) = 1
            else if(grid%lu_mask(i, j) .eq. 1) then
                State_Met%LWI(II, JJ) = 0
            endif
            
            State_Met%PBLH(II, JJ) = grid%pblh(i, j)

            State_Met%PBL_TOP_L(II, JJ) = grid%kpbl(i, j)

            State_Met%PHIS(II, JJ) = grid%ph_2(i, kts, j) + grid%phb(i, kts, j)
            
            State_Met%PRECANV(II, JJ) = (grid%snowncv(i, j) + grid%graupelncv(i, j) + grid%hailncv(i, j))/grid%dt
            
            State_Met%PRECCON(II, JJ) = grid%pratec(i, j)
            
            State_Met%PRECTOT(II, JJ) = (grid%rainncv(i, j) + grid%snowncv(i, j) + grid%graupelncv(i, j) + grid%hailncv(i, j))/grid%dt + grid%pratec(i, j)
            
            State_Met%PRECLSC(II, JJ) = grid%rainncv(i, j)/grid%dt

            State_Met%PS1_WET(II, JJ) = .01*grid%psfc(i, j)
            State_Met%PS2_WET(II, JJ) = .01*(grid%psfc(i, j) + grid%dpsdt(i, j) * grid%dt)
            State_Met%PSC2_WET(II, JJ) = .01*grid%psfc(i, j)

            
            State_Met%PS1_DRY(II, JJ) = .01*grid%psfc(i, j) - 28.97e-3_fp*max(moist(i, kts, j, p_qv) * (1._fp / (1._fp + moist(i, kts, j, p_qv))), 1.0e-30_fp) /&
                                    (18.01e+0_fp + 10.96e-3_fp * max(moist(i, kts, j, p_qv) * (1._fp / (1._fp + moist(i, kts, j, p_qv))), 1.0e-30_fp)) * grid%psfc(i, j)
            State_Met%PS2_DRY(II, JJ) = State_Met%PS1_DRY(II, JJ) + .01 * grid%dpsdt(i, j) * grid%dt
            State_Met%PSC2_DRY(II, JJ) = .01*grid%psfc(i, j) - 28.97e-3_fp*max(moist(i, kts, j, p_qv) * (1._fp / (1._fp + moist(i, kts, j, p_qv))), 1.0e-30_fp) /&
                                    (18.01e+0_fp + 10.96e-3_fp * max(moist(i, kts, j, p_qv) * (1._fp / (1._fp + moist(i, kts, j, p_qv))), 1.0e-30_fp)) * grid%psfc(i, j)
            



            State_Met%SEAICE00(II, JJ) = 1.0_fp

            State_Met%SLP(II, JJ) = grid%pslv_gc(i, j) / 100.


            State_Met%SNODP(II, JJ) = grid%snowh(i, j)

            State_Met%SNOMAS(II, JJ) = grid%acsnow(i, j)
            
            State_Met%SUNCOS(II, JJ) = grid%coszen(i, j)
            
            call calc_zenithb(grid%xlat(i, j), grid%xlong(i, j)*(-1.), grid%julday, grid%dt, grid%gmt, curr_secs, zen)
            State_Met%SUNCOSmid(II, JJ) = cos(zen)

            call getpar(grid%swdown(i, j), 0.01*(grid%p(i, kts, j)+grid%pb(i, kts, j)), zen, pardb, pardif)
            State_Met%PARDR(II, JJ) = pardb
            State_Met%PARDF(II, JJ) = pardif

            State_Met%SWGDN(II, JJ) = grid%swdown(i, j)
            
            if (.not. FIRST) then
                State_Met%TO3(II, JJ) = 0.0_fp
                do k = kts, kte
                   State_Met%TO3(II, JJ) =  State_Met%TO3(II, JJ) + chem(i, k, j, p_o3) * State_Met%AD(II, JJ, k) / (grid%dx * grid%dy / grid%msftx(i, j) / grid%msfty(i, j)) * 22.4e+0_fp / 28.97e+0_fp * 0.1_fp
                enddo
            endif

            State_Met%TROPP(II, JJ) = grid%tropo_p(i, j) / 100.

            State_Met%TS(II, JJ) = grid%t2(i, j)

            State_Met%TSKIN(II, JJ) = grid%tsk(i, j)
            
            State_Met%U10M(II, JJ) = grid%u10(i, j)

            State_Met%USTAR(II, JJ) = grid%ust(i, j)
            
            State_Met%UVALBEDO(II, JJ) = grid%albedo(i, j)
            
            State_Met%V10M(II, JJ) = grid%v10(i, j)

            State_Met%Z0(II, JJ) = grid%znt(i, j)
            

            State_Met%PFLCU(II, JJ, LM+1) = 0._fp
            State_Met%PFICU(II, JJ, LM+1) = 0._fp
            State_Met%PFLLSAN(II, JJ, LM+1) = 0._fp
            State_Met%PFILSAN(II, JJ, LM+1) = 0._fp
            
            State_Met%PEDGE(II, JJ, kts) = 0.01_fp*grid%psfc(i, j)

            State_Met%PEDGE(II, JJ, LM + 1) = 0.01_fp*grid%p_top
        enddo
        enddo

        call wrf_debug(100, "GIGC_Convert_State_Mod after 2-D State_Met% Conversion")

        do k = kts, kte
        do j = jts, jte
        do i = its, ite
            II = i - its + 1
            JJ = j - jts + 1


            State_Met%AREA_M2(II, JJ, 1) = grid%dx * grid%dy / grid%msftx(i, j) / grid%msfty(i, j)

            if(k .ne. kts) then ! Start from kts + 1
                State_Met%PEDGE(II, JJ, k) = .01_fp*(grid%c3f(k)*(grid%psfc(i, j) - grid%p_top) + grid%c4f(k) + grid%p_top)
            endif

            State_Met%CLDF(II, JJ, k) = grid%cldfra(i, k, j)

            if(config_flags%cu_physics .eq. 10) then
                State_Met%CMFMC(II, JJ, k) = grid%mfup_cup(i, k, j)
            else if(config_flags%cu_physics .eq. 7) then
                State_Met%CMFMC(II, JJ, k) = grid%cmfmcdzm(i, k, j)
            else if(config_flags%cu_physics .eq. 16) then
                State_Met%CMFMC(II, JJ, k) = cmfmc(i, k, j)
            endif

            State_Met%T(II, JJ, k) = (grid%t_2(i, k, j) + t0)*((grid%p(i, k, j) + grid%pb(i, k, j))/p1000mb)**rcp

            if(config_flags%cu_physics .eq. 7) then
                State_Met%DQRCU(II, JJ, k) = grid%zmntprpd(i, k, j)
            else if(config_flags%cu_physics .eq. 16) then
                State_Met%DQRCU(II, JJ, k) = dqrcu(i, k, j)
            endif
            
            if(config_flags%mp_physics .eq. 10 .or. config_flags%mp_physics .eq. 8) then
                State_Met%DQRLSAN(II, JJ, k) = grid%rainprod(i, k, j)
            else if(config_flags%mp_physics .eq. 11) then
                State_Met%DQRLSAN(II, JJ, k) = grid%prain3d(i, k, j)
            endif
            
            if(config_flags%cu_physics .eq. 7) then
                State_Met%DTRAIN(II, JJ, k) = grid%du3d(i, k, j)*grid%dp3d(i, k, j)*100._fp/9.8_fp
                if(State_Met%DTRAIN(II, JJ, k) > 0.5) then
                    State_Met%DTRAIN(II, JJ, k) = 0.5
                endif
            else if(config_flags%cu_physics .eq. 16) then
                State_Met%DTRAIN(II, JJ, k) = dtrain(i, k, j)
            endif

            if(State_Met%AIRVOL(II, JJ, k) .ne. 0) then
                State_Met%OMEGA(II, JJ, k) = -9.8_fp*State_Met%AD(II, JJ, k)/State_Met%AIRVOL(II, JJ, k)*grid%w_2(i, k, j)
            endif

            State_Met%QI(II, JJ, k) = moist(i, k, j, p_qi)
            State_Met%QL(II, JJ, k) = moist(i, k, j, p_qc)

            if(config_flags%cu_physics .eq. 7) then
                State_Met%REEVAPCN(II, JJ, k) = abs(grid%evapcdp3d(i, k, j))
            else if(config_flags%cu_physics .eq. 16) then
                State_Met%REEVAPCN(II, JJ, k) = abs(reevapcn(i, k, j)*State_Met%AD(II, JJ, k)/State_Met%AREA_M2(II, JJ, k))
            endif
            
            if(config_flags%mp_physics .eq. 10) then
                State_Met%REEVAPLS(II, JJ, k) = abs(grid%evapprod(i, k, j))
            else if(config_flags%mp_physics .eq. 11) then
                State_Met%REEVAPLS(II, JJ, k) = abs(grid%nevapr3d(i, k, j))
            endif

            State_Met%RH(II, JJ, k) = max(.1, min(.95, moist(i, k, j, p_qv)/ &
                                   (3.80*exp(17.27*(State_Met%T(II, JJ, k) - 273.)/ &
                                             (State_Met%T(II, JJ, k) - 36.))/(.01*(grid%p(i, k, j) + grid%pb(i, k, j)))))) * 100

            State_Met%SPHU(II, JJ, k) = max(moist(i, k, j, p_qv) * (1._fp / (1._fp + moist(i, k, j, p_qv))), 1.0e-30_fp)
            State_Met%SPHU1(II, JJ, k) = State_Met%SPHU(II, JJ, k)
            if(config_flags%cu_physics .eq. 7) then
                State_Met%SPHU2(II, JJ, k) = State_Met%SPHU(II, JJ, k) + grid%zmdq(i, k, j)*1.0e3_fp*grid%dt
            else ! Unsupported physics option
                State_Met%SPHU2(II, JJ, k) = State_Met%SPHU(II, JJ, k)
            endif

            if(i .ne. ide) then
                State_Met%U(II, JJ, k) = 0.5*(grid%u_2(i, k, j) + grid%u_2(i + 1, k, j))
            else
                State_Met%U(II, JJ, k) = grid%u_2(i, k, j)
            endif

            if(j .ne. jde) then
                State_Met%V(II, JJ, k) = 0.5*(grid%v_2(i, k, j) + grid%v_2(i, k, j + 1))
            else
                State_Met%V(II, JJ, k) = grid%v_2(i, k, j)
            endif

            if(config_flags%ra_lw_physics .eq. 3 .or. config_flags%ra_lw_physics .eq. 5) then
                State_Met%TAUCLI(II, JJ, k) = grid%taucldi(i, k, j)
                State_Met%TAUCLW(II, JJ, k) = grid%taucldc(i, k, j)
                State_Met%OPTD(II, JJ, k) = State_Met%TAUCLI(II, JJ, k) + State_Met%TAUCLW(II, JJ, k)
            else if(config_flags%ra_lw_physics .eq. 4) then
                State_Met%TAUCLI(II, JJ, k) = taucldi(i, k, j)
                State_Met%TAUCLW(II, JJ, k) = taucldc(i, k, j)
                State_Met%OPTD(II, JJ, k) = optd(i, k, j)
            endif

            if(config_flags%cu_physics .eq. 7) then
                State_Met%PFLCU(II, JJ, k) = grid%zmflxprc(i, k, j)
                State_Met%PFICU(II, JJ, k) = grid%zmflxsnw(i, k, j)
            else if(config_flags%cu_physics .eq. 16) then
                State_Met%PFLCU(II, JJ, k) = pmflxrain(i, k, j)
                State_Met%PFICU(II, JJ, k) = pmflxsnow(i, k, j)
            endif

            State_Met%PFLLSAN(II, JJ, k) = grid%precr(i, k, j)
            State_Met%PFILSAN(II, JJ, k) = grid%preci(i, k, j) + grid%precs(i, k, j) + grid%precg(i, k, j)

            State_Met%TMPU1(II, JJ, k) = State_Met%T(II, JJ, k)
            State_Met%TMPU2(II, JJ, k) = State_Met%T(II, JJ, k)

            if(State_Met%AIRVOL(II, JJ, k) .ne. 0) then
                State_Met%UPDVVEL(II, JJ, k) = -9.8_fp*State_Met%AD(II, JJ, k)/State_Met%AIRVOL(II, JJ, k)*grid%w_2(i, k, j)
            endif

            grid%diagikj0(i, k, j) = State_Met%CMFMC(II, JJ, k)
            grid%diagikj1(i, k, j) = State_Met%DTRAIN(II, JJ, k)
            grid%diagikj2(i, k, j) = State_Met%PFLCU(II, JJ, k)
            grid%diagikj3(i, k, j) = State_Met%PFLLSAN(II, JJ, k)
            grid%diagikj4(i, k, j) = State_Met%DQRCU(II, JJ, k)
            grid%diagikj5(i, k, j) = State_Met%AD(II, JJ, k)
            grid%diagikj6(i, k, j) = State_Met%AIRVOL(II, JJ, k)
            grid%diagikj7(i, k, j) = State_Met%UPDVVEL(II, JJ, k)
            grid%diagij0(i, j) = State_Met%PRECTOT(II, JJ)
            grid%diagij1(i, j) = State_Met%PRECCON(II, JJ)
            grid%diagij2(i, j) = State_Met%PRECLSC(II, JJ)
            grid%diagij3(i, j) = State_Met%TO3(II, JJ)

            if(k .eq. kte) then
                if(minval(State_Met%PEDGE(II, JJ, 1:LM)) .ge. State_Met%TROPP(II, JJ)) then
                    State_Met%TROPP(II, JJ) = minval(State_Met%PEDGE(II, JJ, 1:LM))
                    write(6, *) "GIGC_Convert_State_Mod fixed T-PAUSE for GEOS-CHEM-CGM, set T-PAUSE to LLPAR"
                endif
            endif

            if(.not. FIRST) then
                State_Chm%Species(II, JJ, k, gi_no) = chem(i, k, j, p_no) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_o3) = chem(i, k, j, p_o3) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_pan) = chem(i, k, j, p_pan) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_co) = chem(i, k, j, p_co) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_alk4) = chem(i, k, j, p_alk4) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_isop) = chem(i, k, j, p_isoprene) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hno3) = chem(i, k, j, p_hno3) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_h2o2) = chem(i, k, j, p_h2o2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_acet) = chem(i, k, j, p_acet) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_mek) = chem(i, k, j, p_mek) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ald2) = chem(i, k, j, p_ald2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_rcho) = chem(i, k, j, p_rcho) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_mvk) = chem(i, k, j, p_mvk) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_macr) = chem(i, k, j, p_macr) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_npmn) = chem(i, k, j, p_npmn) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ppn) = chem(i, k, j, p_ppn) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_r4n2) = chem(i, k, j, p_r4n2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_prpe) = chem(i, k, j, p_prpe) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_c3h8) = chem(i, k, j, p_c3h8) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ch2o) = chem(i, k, j, p_hcho) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_c2h6) = chem(i, k, j, p_c2h6) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_n2o5) = chem(i, k, j, p_n2o5) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hno4) = chem(i, k, j, p_hno4) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_mp) = chem(i, k, j, p_ch3ooh) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_dms) = chem(i, k, j, p_dms) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_so2) = chem(i, k, j, p_so2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_so4) = chem(i, k, j, p_so4) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_so4s) = chem(i, k, j, p_so4s) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_msa) = chem(i, k, j, p_msa) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_nh3) = chem(i, k, j, p_nh3) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_nh4) = chem(i, k, j, p_nh4) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_nit) = chem(i, k, j, p_nit) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_nits) = chem(i, k, j, p_nits) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_bcpi) = chem(i, k, j, p_bcpi) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ocpi) = chem(i, k, j, p_ocpi) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_bcpo) = chem(i, k, j, p_bcpo) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ocpo) = chem(i, k, j, p_ocpo) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_dst1) = chem(i, k, j, p_dst1) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_dst2) = chem(i, k, j, p_dst2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_dst3) = chem(i, k, j, p_dst3) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_dst4) = chem(i, k, j, p_dst4) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_sala) = chem(i, k, j, p_sala) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_salc) = chem(i, k, j, p_salc) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_br2) = chem(i, k, j, p_br2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_br) = chem(i, k, j, p_br) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_bro) = chem(i, k, j, p_bro) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hobr) = chem(i, k, j, p_hobr) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hbr) = chem(i, k, j, p_hbr) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_brno2) = chem(i, k, j, p_brno2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_brno3) = chem(i, k, j, p_brno3) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_chbr3) = chem(i, k, j, p_chbr3) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ch2br2) = chem(i, k, j, p_ch2br2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ch3br) = chem(i, k, j, p_ch3br) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_mpn) = chem(i, k, j, p_mpn) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_isopnd) = chem(i, k, j, p_isopnd) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_isopnb) = chem(i, k, j, p_isopnb) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_moba) = chem(i, k, j, p_moba) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_propnn) = chem(i, k, j, p_propnn) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hac) = chem(i, k, j, p_hac) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_glyc) = chem(i, k, j, p_glyc) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_mvkn) = chem(i, k, j, p_mvkn) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_macrn) = chem(i, k, j, p_macrn) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_map) = chem(i, k, j, p_map) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_no2) = chem(i, k, j, p_no2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_no3) = chem(i, k, j, p_no3) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hno2) = chem(i, k, j, p_hno2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_benz) = chem(i, k, j, p_benzene) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_tolu) = chem(i, k, j, p_toluene) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_xyle) = chem(i, k, j, p_xylenes) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_mtpa) = chem(i, k, j, p_mtpa) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_limo) = chem(i, k, j, p_limon) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_mtpo) = chem(i, k, j, p_mtpo) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_soap) = chem(i, k, j, p_soap) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_soas) = chem(i, k, j, p_soas) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_eoh) = chem(i, k, j, p_eoh) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_mgly) = chem(i, k, j, p_mgly) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_glyx) = chem(i, k, j, p_glyx) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_acta) = chem(i, k, j, p_acta) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hpald) = chem(i, k, j, p_hpald) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_dhdn) = chem(i, k, j, p_dhdn) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ethln) = chem(i, k, j, p_ethln) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hcooh) = chem(i, k, j, p_hcooh) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_iepoxa) = chem(i, k, j, p_iepoxa) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_iepoxb) = chem(i, k, j, p_iepoxb) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_iepoxd) = chem(i, k, j, p_iepoxd) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_isn1) = chem(i, k, j, p_isn1) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ripa) = chem(i, k, j, p_ripa) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ripb) = chem(i, k, j, p_ripb) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ripd) = chem(i, k, j, p_ripd) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_imae) = chem(i, k, j, p_imae) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_soaie) = chem(i, k, j, p_soaie) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_soame) = chem(i, k, j, p_soame) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_soagx) = chem(i, k, j, p_soagx) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_soamg) = chem(i, k, j, p_soamg) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_lvoc) = chem(i, k, j, p_lvoc) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_lvocoa) = chem(i, k, j, p_lvocoa) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_isn1og) = chem(i, k, j, p_isn1og) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_isn1oa) = chem(i, k, j, p_isn1oa) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_monits) = chem(i, k, j, p_monits) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_monitu) = chem(i, k, j, p_monitu) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_honit) = chem(i, k, j, p_honit) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ionita) = chem(i, k, j, p_ionita) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_monita) = chem(i, k, j, p_monita) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_indiol) = chem(i, k, j, p_indiol) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ipmn) = chem(i, k, j, p_ipmn) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hc187) = chem(i, k, j, p_hc187) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_n2o) = chem(i, k, j, p_n2o) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ocs) = chem(i, k, j, p_ocs) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ch4) = chem(i, k, j, p_ch4) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_brcl) = chem(i, k, j, p_brcl) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hcl) = chem(i, k, j, p_hcl) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ccl4) = chem(i, k, j, p_ccl4) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ch3cl) = chem(i, k, j, p_ch3cl) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ch3ccl3) = chem(i, k, j, p_ch3ccl3) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_cfc113) = chem(i, k, j, p_cfc113) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_cfc114) = chem(i, k, j, p_cfc114) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_cfc115) = chem(i, k, j, p_cfc115) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hcfc123) = chem(i, k, j, p_hcfc123) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hcfc141b) = chem(i, k, j, p_hcfc141b) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hcfc142b) = chem(i, k, j, p_hcfc142b) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_cfc11) = chem(i, k, j, p_cfc11) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_cfc12) = chem(i, k, j, p_cfc12) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hcfc22) = chem(i, k, j, p_hcfc22) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_h1211) = chem(i, k, j, p_h1211) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_h1301) = chem(i, k, j, p_h1301) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_h2402) = chem(i, k, j, p_h2402) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_cl) = chem(i, k, j, p_cl) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_clo) = chem(i, k, j, p_clo) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hocl) = chem(i, k, j, p_hocl) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_clno3) = chem(i, k, j, p_clno3) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_clno2) = chem(i, k, j, p_clno2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_cloo) = chem(i, k, j, p_cloo) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_oclo) = chem(i, k, j, p_oclo) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_cl2) = chem(i, k, j, p_cl2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_cl2o2) = chem(i, k, j, p_cl2o2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_h2o) = chem(i, k, j, p_h2o) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_brsala) = chem(i, k, j, p_brsala) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_brsalc) = chem(i, k, j, p_brsalc) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_chcl3) = chem(i, k, j, p_chcl3) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ch2cl2) = chem(i, k, j, p_ch2cl2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ch3i) = chem(i, k, j, p_ch3i) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ch2i2) = chem(i, k, j, p_ch2i2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ch2icl) = chem(i, k, j, p_ch2icl) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ch2ibr) = chem(i, k, j, p_ch2ibr) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hoi) = chem(i, k, j, p_hoi) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_i2) = chem(i, k, j, p_i2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ibr) = chem(i, k, j, p_ibr) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_icl) = chem(i, k, j, p_icl) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_i) = chem(i, k, j, p_i) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_io) = chem(i, k, j, p_io) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hi) = chem(i, k, j, p_hi) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_oio) = chem(i, k, j, p_oio) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ino) = chem(i, k, j, p_ino) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_iono) = chem(i, k, j, p_iono) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_iono2) = chem(i, k, j, p_iono2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_i2o2) = chem(i, k, j, p_i2o2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_i2o3) = chem(i, k, j, p_i2o3) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_i2o4) = chem(i, k, j, p_i2o4) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_isala) = chem(i, k, j, p_isala) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_isalc) = chem(i, k, j, p_isalc) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_aeri) = chem(i, k, j, p_aeri) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_pfe) = chem(i, k, j, p_pfe) * 1.0e-6_fp

            endif
        enddo ! Loop over i
        enddo ! Loop over j
        enddo ! Loop over k

        call wrf_debug(100, "GIGC_Convert_State_Mod after 3-D State_Met% Conversion")

        if(config_flags%have_bcs_chem .and. FIRST) then
            call wrf_debug(1, "GIGC_Convert_State_Mod updating chemistry concentrations in GEOS-Chem from WRF IC values")
            do k = kts, kte
            do j = jts, jte
            do i = its, ite
                II = i - its + 1
                JJ = j - jts + 1

                State_Chm%Species(II, JJ, k, gi_no) = chem(i, k, j, p_no) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_o3) = chem(i, k, j, p_o3) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_pan) = chem(i, k, j, p_pan) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_co) = chem(i, k, j, p_co) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_alk4) = chem(i, k, j, p_alk4) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_isop) = chem(i, k, j, p_isoprene) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hno3) = chem(i, k, j, p_hno3) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_h2o2) = chem(i, k, j, p_h2o2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_acet) = chem(i, k, j, p_acet) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_mek) = chem(i, k, j, p_mek) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ald2) = chem(i, k, j, p_ald2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_mvk) = chem(i, k, j, p_mvk) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_macr) = chem(i, k, j, p_macr) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_c3h8) = chem(i, k, j, p_c3h8) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ch2o) = chem(i, k, j, p_hcho) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_c2h6) = chem(i, k, j, p_c2h6) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_n2o5) = chem(i, k, j, p_n2o5) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hno4) = chem(i, k, j, p_hno4) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_mp) = chem(i, k, j, p_ch3ooh) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_dms) = chem(i, k, j, p_dms) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_so2) = chem(i, k, j, p_so2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_so4) = chem(i, k, j, p_so4) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_nh3) = chem(i, k, j, p_nh3) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_nh4) = chem(i, k, j, p_nh4) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_no2) = chem(i, k, j, p_no2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_no3) = chem(i, k, j, p_no3) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_tolu) = chem(i, k, j, p_toluene) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_limo) = chem(i, k, j, p_limon) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_mgly) = chem(i, k, j, p_mgly) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ch4) = chem(i, k, j, p_ch4) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hc5) = chem(i, k, j, p_hc5) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ho2) = chem(i, k, j, p_ho2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_oh) = chem(i, k, j, p_oh) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_prpe) = chem(i, k, j, p_prpe) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_nit) = chem(i, k, j, p_nit) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_bcpi) = chem(i, k, j, p_bcpi) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ocpi) = chem(i, k, j, p_ocpi) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_bcpo) = chem(i, k, j, p_bcpo) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ocpo) = chem(i, k, j, p_ocpo) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_dst1) = chem(i, k, j, p_dst1) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_dst2) = chem(i, k, j, p_dst2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_dst3) = chem(i, k, j, p_dst3) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_dst4) = chem(i, k, j, p_dst4) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_sala) = chem(i, k, j, p_sala) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_salc) = chem(i, k, j, p_salc) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hac) = chem(i, k, j, p_hac) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_map) = chem(i, k, j, p_map) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_eoh) = chem(i, k, j, p_eoh) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_acta) = chem(i, k, j, p_acta) * 1.0e-6_fp
            enddo
            enddo
            enddo
        endif

        call wrf_debug(100, "GIGC_Convert_State_Mod after State Conversion")

        if(debug_level .ge. 999 .and. FIRST) then
            write(6, *) "%%%%%% GIGC_Convert_State_Mod Chemistry State Export %%%%%%"
            write(6, *) "Units: ", State_Chm%Spc_Units
            do N = 1, State_Chm%nSpecies
                WRITE(6, *) "N:", N, State_Chm%SpcData(N)%Info%Name, &
                            "Is_Adv:", State_Chm%SpcData(N)%Info%Is_Advected
            enddo
        endif

        FIRST = .false.

    end subroutine GIGC_Get_WRF
    subroutine GIGC_Set_WRF(am_I_Root, &
        config_flags, grid, &
        num_chem, chem, num_scalar, scalar, num_moist, moist, &
        its, ite, jts, jte, &
        ide, jde, &
        kts, kte, &
        Input_Opt, State_Met, State_Chm, State_Diag)
        use aerosol_mod, only:   PM25
        implicit none

        logical, intent(in) :: am_I_Root
        type(grid_config_rec_type), intent(in) :: config_flags
        type(domain), target :: grid

        integer, intent(in) :: num_chem, num_scalar, num_moist ! # of WRF chem, moist variables. Check "dummy_new_args.inc" & "actual_new_args.inc" in WRF compile.
        real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33, num_chem), intent(inout)   :: chem
        real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33, num_scalar), intent(inout) :: scalar
        real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33, num_moist), intent(inout)  :: moist

        integer, intent(in) :: its, ite, jts, jte, ide, jde, kts, kte

        type(OptInput), intent(inout) :: Input_Opt
        type(MetState), intent(inout) :: State_Met
        type(ChmState), intent(inout) :: State_Chm
        type(DgnState), intent(inout) :: State_Diag
        integer :: debug_level      ! WRF Debug Level
        integer :: N

        integer :: IM, II           ! # of lons on this PET
        integer :: JM, JJ           ! # of lats on this PET
        integer :: LM               ! # of levs on this PET
        integer :: i, j, k          ! Loop variable over xts, xte (x = i, j, k)

        call nl_get_debug_level(1, debug_level)
        call set_wrf_debug_level(debug_level)

        IM = ite - its + 1
        JM = jte - jts + 1
        LM = kte - kts + 1





        where(State_Chm%Species < 0.0e0)
            State_Chm%Species = 1.0e-36
        endwhere

        call wrf_debug(1, "GIGC_Convert_State_Mod updating chemistry concs/diags in WRF from GEOS-Chem values")

        do j = jts, jte
        do k = kts, kte
        do i = its, ite

            II = i - its + 1
            JJ = j - jts + 1

            if(IT_IS_NAN(State_Chm%Species(II, JJ, k, 7))) then
                State_Chm%Species(II, JJ, k, 8) = 1.0e-36_fp
            endif

            if(IT_IS_NAN(State_Chm%Species(II, JJ, k, 8))) then
                State_Chm%Species(II, JJ, k, 8) = 1.0e-36_fp
            endif

            if(IT_IS_NAN(State_Chm%Species(II, JJ, k, 26))) then
                State_Chm%Species(II, JJ, k, 26) = 1.0e-36_fp
            endif

            chem(i, k, j, p_no) = State_Chm%Species(II, JJ, k, gi_no) * 1.0e+6_fp
            chem(i, k, j, p_o3) = State_Chm%Species(II, JJ, k, gi_o3) * 1.0e+6_fp
            chem(i, k, j, p_pan) = State_Chm%Species(II, JJ, k, gi_pan) * 1.0e+6_fp
            chem(i, k, j, p_co) = State_Chm%Species(II, JJ, k, gi_co) * 1.0e+6_fp
            chem(i, k, j, p_alk4) = State_Chm%Species(II, JJ, k, gi_alk4) * 1.0e+6_fp
            chem(i, k, j, p_isoprene) = State_Chm%Species(II, JJ, k, gi_isop) * 1.0e+6_fp
            chem(i, k, j, p_hno3) = State_Chm%Species(II, JJ, k, gi_hno3) * 1.0e+6_fp
            chem(i, k, j, p_h2o2) = State_Chm%Species(II, JJ, k, gi_h2o2) * 1.0e+6_fp
            chem(i, k, j, p_acet) = State_Chm%Species(II, JJ, k, gi_acet) * 1.0e+6_fp
            chem(i, k, j, p_mek) = State_Chm%Species(II, JJ, k, gi_mek) * 1.0e+6_fp
            chem(i, k, j, p_ald2) = State_Chm%Species(II, JJ, k, gi_ald2) * 1.0e+6_fp
            chem(i, k, j, p_rcho) = State_Chm%Species(II, JJ, k, gi_rcho) * 1.0e+6_fp
            chem(i, k, j, p_mvk) = State_Chm%Species(II, JJ, k, gi_mvk) * 1.0e+6_fp
            chem(i, k, j, p_macr) = State_Chm%Species(II, JJ, k, gi_macr) * 1.0e+6_fp
            chem(i, k, j, p_npmn) = State_Chm%Species(II, JJ, k, gi_npmn) * 1.0e+6_fp
            chem(i, k, j, p_ppn) = State_Chm%Species(II, JJ, k, gi_ppn) * 1.0e+6_fp
            chem(i, k, j, p_r4n2) = State_Chm%Species(II, JJ, k, gi_r4n2) * 1.0e+6_fp
            chem(i, k, j, p_prpe) = State_Chm%Species(II, JJ, k, gi_prpe) * 1.0e+6_fp
            chem(i, k, j, p_c3h8) = State_Chm%Species(II, JJ, k, gi_c3h8) * 1.0e+6_fp
            chem(i, k, j, p_hcho) = State_Chm%Species(II, JJ, k, gi_ch2o) * 1.0e+6_fp
            chem(i, k, j, p_c2h6) = State_Chm%Species(II, JJ, k, gi_c2h6) * 1.0e+6_fp
            chem(i, k, j, p_n2o5) = State_Chm%Species(II, JJ, k, gi_n2o5) * 1.0e+6_fp
            chem(i, k, j, p_hno4) = State_Chm%Species(II, JJ, k, gi_hno4) * 1.0e+6_fp
            chem(i, k, j, p_ch3ooh) = State_Chm%Species(II, JJ, k, gi_mp) * 1.0e+6_fp
            chem(i, k, j, p_dms) = State_Chm%Species(II, JJ, k, gi_dms) * 1.0e+6_fp
            chem(i, k, j, p_so2) = State_Chm%Species(II, JJ, k, gi_so2) * 1.0e+6_fp
            chem(i, k, j, p_so4) = State_Chm%Species(II, JJ, k, gi_so4) * 1.0e+6_fp
            chem(i, k, j, p_so4s) = State_Chm%Species(II, JJ, k, gi_so4s) * 1.0e+6_fp
            chem(i, k, j, p_msa) = State_Chm%Species(II, JJ, k, gi_msa) * 1.0e+6_fp
            chem(i, k, j, p_nh3) = State_Chm%Species(II, JJ, k, gi_nh3) * 1.0e+6_fp
            chem(i, k, j, p_nh4) = State_Chm%Species(II, JJ, k, gi_nh4) * 1.0e+6_fp
            chem(i, k, j, p_nit) = State_Chm%Species(II, JJ, k, gi_nit) * 1.0e+6_fp
            chem(i, k, j, p_nits) = State_Chm%Species(II, JJ, k, gi_nits) * 1.0e+6_fp
            chem(i, k, j, p_bcpi) = State_Chm%Species(II, JJ, k, gi_bcpi) * 1.0e+6_fp
            chem(i, k, j, p_ocpi) = State_Chm%Species(II, JJ, k, gi_ocpi) * 1.0e+6_fp
            chem(i, k, j, p_bcpo) = State_Chm%Species(II, JJ, k, gi_bcpo) * 1.0e+6_fp
            chem(i, k, j, p_ocpo) = State_Chm%Species(II, JJ, k, gi_ocpo) * 1.0e+6_fp
            chem(i, k, j, p_dst1) = State_Chm%Species(II, JJ, k, gi_dst1) * 1.0e+6_fp
            chem(i, k, j, p_dst2) = State_Chm%Species(II, JJ, k, gi_dst2) * 1.0e+6_fp
            chem(i, k, j, p_dst3) = State_Chm%Species(II, JJ, k, gi_dst3) * 1.0e+6_fp
            chem(i, k, j, p_dst4) = State_Chm%Species(II, JJ, k, gi_dst4) * 1.0e+6_fp
            chem(i, k, j, p_sala) = State_Chm%Species(II, JJ, k, gi_sala) * 1.0e+6_fp
            chem(i, k, j, p_salc) = State_Chm%Species(II, JJ, k, gi_salc) * 1.0e+6_fp
            chem(i, k, j, p_br2) = State_Chm%Species(II, JJ, k, gi_br2) * 1.0e+6_fp
            chem(i, k, j, p_br) = State_Chm%Species(II, JJ, k, gi_br) * 1.0e+6_fp
            chem(i, k, j, p_bro) = State_Chm%Species(II, JJ, k, gi_bro) * 1.0e+6_fp
            chem(i, k, j, p_hobr) = State_Chm%Species(II, JJ, k, gi_hobr) * 1.0e+6_fp
            chem(i, k, j, p_hbr) = State_Chm%Species(II, JJ, k, gi_hbr) * 1.0e+6_fp
            chem(i, k, j, p_brno2) = State_Chm%Species(II, JJ, k, gi_brno2) * 1.0e+6_fp
            chem(i, k, j, p_brno3) = State_Chm%Species(II, JJ, k, gi_brno3) * 1.0e+6_fp
            chem(i, k, j, p_chbr3) = State_Chm%Species(II, JJ, k, gi_chbr3) * 1.0e+6_fp
            chem(i, k, j, p_ch2br2) = State_Chm%Species(II, JJ, k, gi_ch2br2) * 1.0e+6_fp
            chem(i, k, j, p_ch3br) = State_Chm%Species(II, JJ, k, gi_ch3br) * 1.0e+6_fp
            chem(i, k, j, p_mpn) = State_Chm%Species(II, JJ, k, gi_mpn) * 1.0e+6_fp
            chem(i, k, j, p_isopnd) = State_Chm%Species(II, JJ, k, gi_isopnd) * 1.0e+6_fp
            chem(i, k, j, p_isopnb) = State_Chm%Species(II, JJ, k, gi_isopnb) * 1.0e+6_fp
            chem(i, k, j, p_moba) = State_Chm%Species(II, JJ, k, gi_moba) * 1.0e+6_fp
            chem(i, k, j, p_propnn) = State_Chm%Species(II, JJ, k, gi_propnn) * 1.0e+6_fp
            chem(i, k, j, p_hac) = State_Chm%Species(II, JJ, k, gi_hac) * 1.0e+6_fp
            chem(i, k, j, p_glyc) = State_Chm%Species(II, JJ, k, gi_glyc) * 1.0e+6_fp
            chem(i, k, j, p_mvkn) = State_Chm%Species(II, JJ, k, gi_mvkn) * 1.0e+6_fp
            chem(i, k, j, p_macrn) = State_Chm%Species(II, JJ, k, gi_macrn) * 1.0e+6_fp
            chem(i, k, j, p_map) = State_Chm%Species(II, JJ, k, gi_map) * 1.0e+6_fp
            chem(i, k, j, p_no2) = State_Chm%Species(II, JJ, k, gi_no2) * 1.0e+6_fp
            chem(i, k, j, p_no3) = State_Chm%Species(II, JJ, k, gi_no3) * 1.0e+6_fp
            chem(i, k, j, p_hno2) = State_Chm%Species(II, JJ, k, gi_hno2) * 1.0e+6_fp
            chem(i, k, j, p_benzene) = State_Chm%Species(II, JJ, k, gi_benz) * 1.0e+6_fp
            chem(i, k, j, p_toluene) = State_Chm%Species(II, JJ, k, gi_tolu) * 1.0e+6_fp
            chem(i, k, j, p_xylenes) = State_Chm%Species(II, JJ, k, gi_xyle) * 1.0e+6_fp
            chem(i, k, j, p_mtpa) = State_Chm%Species(II, JJ, k, gi_mtpa) * 1.0e+6_fp
            chem(i, k, j, p_limon) = State_Chm%Species(II, JJ, k, gi_limo) * 1.0e+6_fp
            chem(i, k, j, p_mtpo) = State_Chm%Species(II, JJ, k, gi_mtpo) * 1.0e+6_fp
            chem(i, k, j, p_soap) = State_Chm%Species(II, JJ, k, gi_soap) * 1.0e+6_fp
            chem(i, k, j, p_soas) = State_Chm%Species(II, JJ, k, gi_soas) * 1.0e+6_fp
            chem(i, k, j, p_eoh) = State_Chm%Species(II, JJ, k, gi_eoh) * 1.0e+6_fp
            chem(i, k, j, p_mgly) = State_Chm%Species(II, JJ, k, gi_mgly) * 1.0e+6_fp
            chem(i, k, j, p_glyx) = State_Chm%Species(II, JJ, k, gi_glyx) * 1.0e+6_fp
            chem(i, k, j, p_acta) = State_Chm%Species(II, JJ, k, gi_acta) * 1.0e+6_fp
            chem(i, k, j, p_hpald) = State_Chm%Species(II, JJ, k, gi_hpald) * 1.0e+6_fp
            chem(i, k, j, p_dhdn) = State_Chm%Species(II, JJ, k, gi_dhdn) * 1.0e+6_fp
            chem(i, k, j, p_ethln) = State_Chm%Species(II, JJ, k, gi_ethln) * 1.0e+6_fp
            chem(i, k, j, p_hcooh) = State_Chm%Species(II, JJ, k, gi_hcooh) * 1.0e+6_fp
            chem(i, k, j, p_iepoxa) = State_Chm%Species(II, JJ, k, gi_iepoxa) * 1.0e+6_fp
            chem(i, k, j, p_iepoxb) = State_Chm%Species(II, JJ, k, gi_iepoxb) * 1.0e+6_fp
            chem(i, k, j, p_iepoxd) = State_Chm%Species(II, JJ, k, gi_iepoxd) * 1.0e+6_fp
            chem(i, k, j, p_isn1) = State_Chm%Species(II, JJ, k, gi_isn1) * 1.0e+6_fp
            chem(i, k, j, p_ripa) = State_Chm%Species(II, JJ, k, gi_ripa) * 1.0e+6_fp
            chem(i, k, j, p_ripb) = State_Chm%Species(II, JJ, k, gi_ripb) * 1.0e+6_fp
            chem(i, k, j, p_ripd) = State_Chm%Species(II, JJ, k, gi_ripd) * 1.0e+6_fp
            chem(i, k, j, p_imae) = State_Chm%Species(II, JJ, k, gi_imae) * 1.0e+6_fp
            chem(i, k, j, p_soaie) = State_Chm%Species(II, JJ, k, gi_soaie) * 1.0e+6_fp
            chem(i, k, j, p_soame) = State_Chm%Species(II, JJ, k, gi_soame) * 1.0e+6_fp
            chem(i, k, j, p_soagx) = State_Chm%Species(II, JJ, k, gi_soagx) * 1.0e+6_fp
            chem(i, k, j, p_soamg) = State_Chm%Species(II, JJ, k, gi_soamg) * 1.0e+6_fp
            chem(i, k, j, p_lvoc) = State_Chm%Species(II, JJ, k, gi_lvoc) * 1.0e+6_fp
            chem(i, k, j, p_lvocoa) = State_Chm%Species(II, JJ, k, gi_lvocoa) * 1.0e+6_fp
            chem(i, k, j, p_isn1og) = State_Chm%Species(II, JJ, k, gi_isn1og) * 1.0e+6_fp
            chem(i, k, j, p_isn1oa) = State_Chm%Species(II, JJ, k, gi_isn1oa) * 1.0e+6_fp
            chem(i, k, j, p_monits) = State_Chm%Species(II, JJ, k, gi_monits) * 1.0e+6_fp
            chem(i, k, j, p_monitu) = State_Chm%Species(II, JJ, k, gi_monitu) * 1.0e+6_fp
            chem(i, k, j, p_honit) = State_Chm%Species(II, JJ, k, gi_honit) * 1.0e+6_fp
            chem(i, k, j, p_ionita) = State_Chm%Species(II, JJ, k, gi_ionita) * 1.0e+6_fp
            chem(i, k, j, p_monita) = State_Chm%Species(II, JJ, k, gi_monita) * 1.0e+6_fp
            chem(i, k, j, p_indiol) = State_Chm%Species(II, JJ, k, gi_indiol) * 1.0e+6_fp
            chem(i, k, j, p_ipmn) = State_Chm%Species(II, JJ, k, gi_ipmn) * 1.0e+6_fp
            chem(i, k, j, p_hc187) = State_Chm%Species(II, JJ, k, gi_hc187) * 1.0e+6_fp
            chem(i, k, j, p_n2o) = State_Chm%Species(II, JJ, k, gi_n2o) * 1.0e+6_fp
            chem(i, k, j, p_ocs) = State_Chm%Species(II, JJ, k, gi_ocs) * 1.0e+6_fp
            chem(i, k, j, p_ch4) = State_Chm%Species(II, JJ, k, gi_ch4) * 1.0e+6_fp
            chem(i, k, j, p_brcl) = State_Chm%Species(II, JJ, k, gi_brcl) * 1.0e+6_fp
            chem(i, k, j, p_hcl) = State_Chm%Species(II, JJ, k, gi_hcl) * 1.0e+6_fp
            chem(i, k, j, p_ccl4) = State_Chm%Species(II, JJ, k, gi_ccl4) * 1.0e+6_fp
            chem(i, k, j, p_ch3cl) = State_Chm%Species(II, JJ, k, gi_ch3cl) * 1.0e+6_fp
            chem(i, k, j, p_ch3ccl3) = State_Chm%Species(II, JJ, k, gi_ch3ccl3) * 1.0e+6_fp
            chem(i, k, j, p_cfc113) = State_Chm%Species(II, JJ, k, gi_cfc113) * 1.0e+6_fp
            chem(i, k, j, p_cfc114) = State_Chm%Species(II, JJ, k, gi_cfc114) * 1.0e+6_fp
            chem(i, k, j, p_cfc115) = State_Chm%Species(II, JJ, k, gi_cfc115) * 1.0e+6_fp
            chem(i, k, j, p_hcfc123) = State_Chm%Species(II, JJ, k, gi_hcfc123) * 1.0e+6_fp
            chem(i, k, j, p_hcfc141b) = State_Chm%Species(II, JJ, k, gi_hcfc141b) * 1.0e+6_fp
            chem(i, k, j, p_hcfc142b) = State_Chm%Species(II, JJ, k, gi_hcfc142b) * 1.0e+6_fp
            chem(i, k, j, p_cfc11) = State_Chm%Species(II, JJ, k, gi_cfc11) * 1.0e+6_fp
            chem(i, k, j, p_cfc12) = State_Chm%Species(II, JJ, k, gi_cfc12) * 1.0e+6_fp
            chem(i, k, j, p_hcfc22) = State_Chm%Species(II, JJ, k, gi_hcfc22) * 1.0e+6_fp
            chem(i, k, j, p_h1211) = State_Chm%Species(II, JJ, k, gi_h1211) * 1.0e+6_fp
            chem(i, k, j, p_h1301) = State_Chm%Species(II, JJ, k, gi_h1301) * 1.0e+6_fp
            chem(i, k, j, p_h2402) = State_Chm%Species(II, JJ, k, gi_h2402) * 1.0e+6_fp
            chem(i, k, j, p_cl) = State_Chm%Species(II, JJ, k, gi_cl) * 1.0e+6_fp
            chem(i, k, j, p_clo) = State_Chm%Species(II, JJ, k, gi_clo) * 1.0e+6_fp
            chem(i, k, j, p_hocl) = State_Chm%Species(II, JJ, k, gi_hocl) * 1.0e+6_fp
            chem(i, k, j, p_clno3) = State_Chm%Species(II, JJ, k, gi_clno3) * 1.0e+6_fp
            chem(i, k, j, p_clno2) = State_Chm%Species(II, JJ, k, gi_clno2) * 1.0e+6_fp
            chem(i, k, j, p_cloo) = State_Chm%Species(II, JJ, k, gi_cloo) * 1.0e+6_fp
            chem(i, k, j, p_oclo) = State_Chm%Species(II, JJ, k, gi_oclo) * 1.0e+6_fp
            chem(i, k, j, p_cl2) = State_Chm%Species(II, JJ, k, gi_cl2) * 1.0e+6_fp
            chem(i, k, j, p_cl2o2) = State_Chm%Species(II, JJ, k, gi_cl2o2) * 1.0e+6_fp
            chem(i, k, j, p_h2o) = State_Chm%Species(II, JJ, k, gi_h2o) * 1.0e+6_fp
            chem(i, k, j, p_brsala) = State_Chm%Species(II, JJ, k, gi_brsala) * 1.0e+6_fp
            chem(i, k, j, p_brsalc) = State_Chm%Species(II, JJ, k, gi_brsalc) * 1.0e+6_fp
            chem(i, k, j, p_chcl3) = State_Chm%Species(II, JJ, k, gi_chcl3) * 1.0e+6_fp
            chem(i, k, j, p_ch2cl2) = State_Chm%Species(II, JJ, k, gi_ch2cl2) * 1.0e+6_fp
            chem(i, k, j, p_ch3i) = State_Chm%Species(II, JJ, k, gi_ch3i) * 1.0e+6_fp
            chem(i, k, j, p_ch2i2) = State_Chm%Species(II, JJ, k, gi_ch2i2) * 1.0e+6_fp
            chem(i, k, j, p_ch2icl) = State_Chm%Species(II, JJ, k, gi_ch2icl) * 1.0e+6_fp
            chem(i, k, j, p_ch2ibr) = State_Chm%Species(II, JJ, k, gi_ch2ibr) * 1.0e+6_fp
            chem(i, k, j, p_hoi) = State_Chm%Species(II, JJ, k, gi_hoi) * 1.0e+6_fp
            chem(i, k, j, p_i2) = State_Chm%Species(II, JJ, k, gi_i2) * 1.0e+6_fp
            chem(i, k, j, p_ibr) = State_Chm%Species(II, JJ, k, gi_ibr) * 1.0e+6_fp
            chem(i, k, j, p_icl) = State_Chm%Species(II, JJ, k, gi_icl) * 1.0e+6_fp
            chem(i, k, j, p_i) = State_Chm%Species(II, JJ, k, gi_i) * 1.0e+6_fp
            chem(i, k, j, p_io) = State_Chm%Species(II, JJ, k, gi_io) * 1.0e+6_fp
            chem(i, k, j, p_hi) = State_Chm%Species(II, JJ, k, gi_hi) * 1.0e+6_fp
            chem(i, k, j, p_oio) = State_Chm%Species(II, JJ, k, gi_oio) * 1.0e+6_fp
            chem(i, k, j, p_ino) = State_Chm%Species(II, JJ, k, gi_ino) * 1.0e+6_fp
            chem(i, k, j, p_iono) = State_Chm%Species(II, JJ, k, gi_iono) * 1.0e+6_fp
            chem(i, k, j, p_iono2) = State_Chm%Species(II, JJ, k, gi_iono2) * 1.0e+6_fp
            chem(i, k, j, p_i2o2) = State_Chm%Species(II, JJ, k, gi_i2o2) * 1.0e+6_fp
            chem(i, k, j, p_i2o3) = State_Chm%Species(II, JJ, k, gi_i2o3) * 1.0e+6_fp
            chem(i, k, j, p_i2o4) = State_Chm%Species(II, JJ, k, gi_i2o4) * 1.0e+6_fp
            chem(i, k, j, p_isala) = State_Chm%Species(II, JJ, k, gi_isala) * 1.0e+6_fp
            chem(i, k, j, p_isalc) = State_Chm%Species(II, JJ, k, gi_isalc) * 1.0e+6_fp
            chem(i, k, j, p_aeri) = State_Chm%Species(II, JJ, k, gi_aeri) * 1.0e+6_fp
            chem(i, k, j, p_pfe) = State_Chm%Species(II, JJ, k, gi_pfe) * 1.0e+6_fp
            chem(i, k, j, p_co2) = State_Chm%Species(II, JJ, k, gi_co2) * 1.0e+6_fp
            chem(i, k, j, p_lbro2h) = State_Chm%Species(II, JJ, k, gi_lbro2h) * 1.0e+6_fp
            chem(i, k, j, p_lbro2n) = State_Chm%Species(II, JJ, k, gi_lbro2n) * 1.0e+6_fp
            chem(i, k, j, p_lisopoh) = State_Chm%Species(II, JJ, k, gi_lisopoh) * 1.0e+6_fp
            chem(i, k, j, p_lisopno3) = State_Chm%Species(II, JJ, k, gi_lisopno3) * 1.0e+6_fp
            chem(i, k, j, p_ltro2h) = State_Chm%Species(II, JJ, k, gi_ltro2h) * 1.0e+6_fp
            chem(i, k, j, p_ltro2n) = State_Chm%Species(II, JJ, k, gi_ltro2n) * 1.0e+6_fp
            chem(i, k, j, p_lxro2h) = State_Chm%Species(II, JJ, k, gi_lxro2h) * 1.0e+6_fp
            chem(i, k, j, p_lxro2n) = State_Chm%Species(II, JJ, k, gi_lxro2n) * 1.0e+6_fp
            chem(i, k, j, p_pyac) = State_Chm%Species(II, JJ, k, gi_pyac) * 1.0e+6_fp
            chem(i, k, j, p_so4h1) = State_Chm%Species(II, JJ, k, gi_so4h1) * 1.0e+6_fp
            chem(i, k, j, p_so4h2) = State_Chm%Species(II, JJ, k, gi_so4h2) * 1.0e+6_fp
            chem(i, k, j, p_pox) = State_Chm%Species(II, JJ, k, gi_pox) * 1.0e+6_fp
            chem(i, k, j, p_lox) = State_Chm%Species(II, JJ, k, gi_lox) * 1.0e+6_fp
            chem(i, k, j, p_pco) = State_Chm%Species(II, JJ, k, gi_pco) * 1.0e+6_fp
            chem(i, k, j, p_lco) = State_Chm%Species(II, JJ, k, gi_lco) * 1.0e+6_fp
            chem(i, k, j, p_pso4) = State_Chm%Species(II, JJ, k, gi_pso4) * 1.0e+6_fp
            chem(i, k, j, p_lch4) = State_Chm%Species(II, JJ, k, gi_lch4) * 1.0e+6_fp
            chem(i, k, j, p_ph2o2) = State_Chm%Species(II, JJ, k, gi_ph2o2) * 1.0e+6_fp
            chem(i, k, j, p_dhdc) = State_Chm%Species(II, JJ, k, gi_dhdc) * 1.0e+6_fp
            chem(i, k, j, p_pmnn) = State_Chm%Species(II, JJ, k, gi_pmnn) * 1.0e+6_fp
            chem(i, k, j, p_tro2) = State_Chm%Species(II, JJ, k, gi_tro2) * 1.0e+6_fp
            chem(i, k, j, p_bro2) = State_Chm%Species(II, JJ, k, gi_bro2) * 1.0e+6_fp
            chem(i, k, j, p_n) = State_Chm%Species(II, JJ, k, gi_n) * 1.0e+6_fp
            chem(i, k, j, p_xro2) = State_Chm%Species(II, JJ, k, gi_xro2) * 1.0e+6_fp
            chem(i, k, j, p_imao3) = State_Chm%Species(II, JJ, k, gi_imao3) * 1.0e+6_fp
            chem(i, k, j, p_etp) = State_Chm%Species(II, JJ, k, gi_etp) * 1.0e+6_fp
            chem(i, k, j, p_ra3p) = State_Chm%Species(II, JJ, k, gi_ra3p) * 1.0e+6_fp
            chem(i, k, j, p_rb3p) = State_Chm%Species(II, JJ, k, gi_rb3p) * 1.0e+6_fp
            chem(i, k, j, p_rp) = State_Chm%Species(II, JJ, k, gi_rp) * 1.0e+6_fp
            chem(i, k, j, p_pp) = State_Chm%Species(II, JJ, k, gi_pp) * 1.0e+6_fp
            chem(i, k, j, p_prpn) = State_Chm%Species(II, JJ, k, gi_prpn) * 1.0e+6_fp
            chem(i, k, j, p_pip) = State_Chm%Species(II, JJ, k, gi_pip) * 1.0e+6_fp
            chem(i, k, j, p_r4p) = State_Chm%Species(II, JJ, k, gi_r4p) * 1.0e+6_fp
            chem(i, k, j, p_dhpcarp) = State_Chm%Species(II, JJ, k, gi_dhpcarp) * 1.0e+6_fp
            chem(i, k, j, p_iap) = State_Chm%Species(II, JJ, k, gi_iap) * 1.0e+6_fp
            chem(i, k, j, p_hpc52o2) = State_Chm%Species(II, JJ, k, gi_hpc52o2) * 1.0e+6_fp
            chem(i, k, j, p_vrp) = State_Chm%Species(II, JJ, k, gi_vrp) * 1.0e+6_fp
            chem(i, k, j, p_atooh) = State_Chm%Species(II, JJ, k, gi_atooh) * 1.0e+6_fp
            chem(i, k, j, p_dhmob) = State_Chm%Species(II, JJ, k, gi_dhmob) * 1.0e+6_fp
            chem(i, k, j, p_isnp) = State_Chm%Species(II, JJ, k, gi_isnp) * 1.0e+6_fp
            chem(i, k, j, p_maop) = State_Chm%Species(II, JJ, k, gi_maop) * 1.0e+6_fp
            chem(i, k, j, p_mrp) = State_Chm%Species(II, JJ, k, gi_mrp) * 1.0e+6_fp
            chem(i, k, j, p_inpn) = State_Chm%Species(II, JJ, k, gi_inpn) * 1.0e+6_fp
            chem(i, k, j, p_mobaoo) = State_Chm%Species(II, JJ, k, gi_mobaoo) * 1.0e+6_fp
            chem(i, k, j, p_diboo) = State_Chm%Species(II, JJ, k, gi_diboo) * 1.0e+6_fp
            chem(i, k, j, p_h) = State_Chm%Species(II, JJ, k, gi_h) * 1.0e+6_fp
            chem(i, k, j, p_macrno2) = State_Chm%Species(II, JJ, k, gi_macrno2) * 1.0e+6_fp
            chem(i, k, j, p_roh) = State_Chm%Species(II, JJ, k, gi_roh) * 1.0e+6_fp
            chem(i, k, j, p_isnohoo) = State_Chm%Species(II, JJ, k, gi_isnohoo) * 1.0e+6_fp
            chem(i, k, j, p_mvkoo) = State_Chm%Species(II, JJ, k, gi_mvkoo) * 1.0e+6_fp
            chem(i, k, j, p_isnoob) = State_Chm%Species(II, JJ, k, gi_isnoob) * 1.0e+6_fp
            chem(i, k, j, p_gaoo) = State_Chm%Species(II, JJ, k, gi_gaoo) * 1.0e+6_fp
            chem(i, k, j, p_ch3choo) = State_Chm%Species(II, JJ, k, gi_ch3choo) * 1.0e+6_fp
            chem(i, k, j, p_iepoxoo) = State_Chm%Species(II, JJ, k, gi_iepoxoo) * 1.0e+6_fp
            chem(i, k, j, p_mglyoo) = State_Chm%Species(II, JJ, k, gi_mglyoo) * 1.0e+6_fp
            chem(i, k, j, p_prn1) = State_Chm%Species(II, JJ, k, gi_prn1) * 1.0e+6_fp
            chem(i, k, j, p_mgloo) = State_Chm%Species(II, JJ, k, gi_mgloo) * 1.0e+6_fp
            chem(i, k, j, p_a3o2) = State_Chm%Species(II, JJ, k, gi_a3o2) * 1.0e+6_fp
            chem(i, k, j, p_man2) = State_Chm%Species(II, JJ, k, gi_man2) * 1.0e+6_fp
            chem(i, k, j, p_isnooa) = State_Chm%Species(II, JJ, k, gi_isnooa) * 1.0e+6_fp
            chem(i, k, j, p_po2) = State_Chm%Species(II, JJ, k, gi_po2) * 1.0e+6_fp
            chem(i, k, j, p_isopndo2) = State_Chm%Species(II, JJ, k, gi_isopndo2) * 1.0e+6_fp
            chem(i, k, j, p_b3o2) = State_Chm%Species(II, JJ, k, gi_b3o2) * 1.0e+6_fp
            chem(i, k, j, p_macroo) = State_Chm%Species(II, JJ, k, gi_macroo) * 1.0e+6_fp
            chem(i, k, j, p_r4n1) = State_Chm%Species(II, JJ, k, gi_r4n1) * 1.0e+6_fp
            chem(i, k, j, p_maopo2) = State_Chm%Species(II, JJ, k, gi_maopo2) * 1.0e+6_fp
            chem(i, k, j, p_ato2) = State_Chm%Species(II, JJ, k, gi_ato2) * 1.0e+6_fp
            chem(i, k, j, p_rco3) = State_Chm%Species(II, JJ, k, gi_rco3) * 1.0e+6_fp
            chem(i, k, j, p_limo2) = State_Chm%Species(II, JJ, k, gi_limo2) * 1.0e+6_fp
            chem(i, k, j, p_olnd) = State_Chm%Species(II, JJ, k, gi_olnd) * 1.0e+6_fp
            chem(i, k, j, p_olnn) = State_Chm%Species(II, JJ, k, gi_olnn) * 1.0e+6_fp
            chem(i, k, j, p_ko2) = State_Chm%Species(II, JJ, k, gi_ko2) * 1.0e+6_fp
            chem(i, k, j, p_isopnbo2) = State_Chm%Species(II, JJ, k, gi_isopnbo2) * 1.0e+6_fp
            chem(i, k, j, p_pio2) = State_Chm%Species(II, JJ, k, gi_pio2) * 1.0e+6_fp
            chem(i, k, j, p_hc5oo) = State_Chm%Species(II, JJ, k, gi_hc5oo) * 1.0e+6_fp
            chem(i, k, j, p_nmao3) = State_Chm%Species(II, JJ, k, gi_nmao3) * 1.0e+6_fp
            chem(i, k, j, p_vro2) = State_Chm%Species(II, JJ, k, gi_vro2) * 1.0e+6_fp
            chem(i, k, j, p_ch2oo) = State_Chm%Species(II, JJ, k, gi_ch2oo) * 1.0e+6_fp
            chem(i, k, j, p_hc5) = State_Chm%Species(II, JJ, k, gi_hc5) * 1.0e+6_fp
            chem(i, k, j, p_rio2) = State_Chm%Species(II, JJ, k, gi_rio2) * 1.0e+6_fp
            chem(i, k, j, p_eto2) = State_Chm%Species(II, JJ, k, gi_eto2) * 1.0e+6_fp
            chem(i, k, j, p_ino2) = State_Chm%Species(II, JJ, k, gi_ino2) * 1.0e+6_fp
            chem(i, k, j, p_r4o2) = State_Chm%Species(II, JJ, k, gi_r4o2) * 1.0e+6_fp
            chem(i, k, j, p_mro2) = State_Chm%Species(II, JJ, k, gi_mro2) * 1.0e+6_fp
            chem(i, k, j, p_mco3) = State_Chm%Species(II, JJ, k, gi_mco3) * 1.0e+6_fp
            chem(i, k, j, p_mo2) = State_Chm%Species(II, JJ, k, gi_mo2) * 1.0e+6_fp
            chem(i, k, j, p_o1d) = State_Chm%Species(II, JJ, k, gi_o1d) * 1.0e+6_fp
            chem(i, k, j, p_o) = State_Chm%Species(II, JJ, k, gi_o) * 1.0e+6_fp
            chem(i, k, j, p_ho2) = State_Chm%Species(II, JJ, k, gi_ho2) * 1.0e+6_fp
            chem(i, k, j, p_oh) = State_Chm%Species(II, JJ, k, gi_oh) * 1.0e+6_fp
            chem(i, k, j, p_h2) = State_Chm%Species(II, JJ, k, gi_h2) * 1.0e+6_fp
            chem(i, k, j, p_moh) = State_Chm%Species(II, JJ, k, gi_moh) * 1.0e+6_fp
            chem(i, k, j, p_n2) = State_Chm%Species(II, JJ, k, gi_n2) * 1.0e+6_fp
            chem(i, k, j, p_o2) = State_Chm%Species(II, JJ, k, gi_o2) * 1.0e+6_fp
            chem(i, k, j, p_rcooh) = State_Chm%Species(II, JJ, k, gi_rcooh) * 1.0e+6_fp

            if(allocated(PM25)) then
                grid%pm2_5_dry(i, k, j) = PM25(II, JJ, k) * 1.0e+9_fp ! kg/m^3 -> \mug/m^3
            endif








            if(config_flags%gc_diagn_spc_n0 .ne. 0 .and. config_flags%gc_diagn_spc_n0 .le. State_Chm%nSpecies) then
                grid%cldcnvflx_n0(i, k, j) = State_Diag%CloudConvFlux(II, JJ, k, config_flags%gc_diagn_spc_n0)


                grid%gcemisdrydep_full_n0(i, j) = State_Diag%BudgetEmisDryDepFull(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n0)%Info%AdvectID)
                grid%gcemisdrydep_trop_n0(i, j) = State_Diag%BudgetEmisDryDepTrop(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n0)%Info%AdvectID)
                grid%gcemisdrydep_pbl_n0(i, j)  = State_Diag%BudgetEmisDryDepPBL(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n0)%Info%AdvectID)

                grid%gcmixing_full_n0(i, j) = State_Diag%BudgetMixingFull(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n0)%Info%AdvectID)
                grid%gcmixing_trop_n0(i, j) = State_Diag%BudgetMixingTrop(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n0)%Info%AdvectID)
                grid%gcmixing_pbl_n0(i, j)  = State_Diag%BudgetMixingPBL(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n0)%Info%AdvectID)

                grid%gcconv_full_n0(i, j) = State_Diag%BudgetConvectionFull(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n0)%Info%AdvectID)
                grid%gcconv_trop_n0(i, j) = State_Diag%BudgetConvectionTrop(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n0)%Info%AdvectID)
                grid%gcconv_pbl_n0(i, j)  = State_Diag%BudgetConvectionPBL(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n0)%Info%AdvectID)

                grid%gcchem_full_n0(i, j) = State_Diag%BudgetChemistryFull(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n0)%Info%AdvectID)
                grid%gcchem_trop_n0(i, j) = State_Diag%BudgetChemistryTrop(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n0)%Info%AdvectID)
                grid%gcchem_pbl_n0(i, j)  = State_Diag%BudgetChemistryPBL(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n0)%Info%AdvectID)

                if(State_Chm%SpcData(config_flags%gc_diagn_spc_n0)%Info%WetDepId .ne. MISSING_INT) then
                    grid%gcwetdep_full_n0(i, j) = State_Diag%BudgetWetDepFull(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n0)%Info%WetDepId)
                    grid%gcwetdep_trop_n0(i, j) = State_Diag%BudgetWetDepTrop(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n0)%Info%WetDepId)
                    grid%gcwetdep_pbl_n0(i, j)  = State_Diag%BudgetWetDepPBL(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n0)%Info%WetDepId)
                    grid%wetlscnvfrc_n0(i, k, j) = State_Diag%WetLossConvFrac(II, JJ, k, &
                                                   State_Chm%SpcData(config_flags%gc_diagn_spc_n0)%Info%WetDepId)
                    grid%wetlscnv_n0(i, k, j) = State_Diag%WetLossConv(II, JJ, k, &
                                                   State_Chm%SpcData(config_flags%gc_diagn_spc_n0)%Info%WetDepId)

                    grid%wetlossls_n0(i, k, j) = State_Diag%WetLossLS(II, JJ, k, &
                                                   State_Chm%SpcData(config_flags%gc_diagn_spc_n0)%Info%WetDepId)
                endif
            endif

            if(config_flags%gc_diagn_spc_n1 .ne. 0 .and. config_flags%gc_diagn_spc_n1 .le. State_Chm%nSpecies) then
                grid%cldcnvflx_n1(i, k, j) = State_Diag%CloudConvFlux(II, JJ, k, config_flags%gc_diagn_spc_n1)

                grid%gcemisdrydep_full_n1(i, j) = State_Diag%BudgetEmisDryDepFull(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n1)%Info%AdvectID)
                grid%gcemisdrydep_trop_n1(i, j) = State_Diag%BudgetEmisDryDepTrop(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n1)%Info%AdvectID)
                grid%gcemisdrydep_pbl_n1(i, j)  = State_Diag%BudgetEmisDryDepPBL(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n1)%Info%AdvectID)

                grid%gcmixing_full_n1(i, j) = State_Diag%BudgetMixingFull(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n1)%Info%AdvectID)
                grid%gcmixing_trop_n1(i, j) = State_Diag%BudgetMixingTrop(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n1)%Info%AdvectID)
                grid%gcmixing_pbl_n1(i, j)  = State_Diag%BudgetMixingPBL(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n1)%Info%AdvectID)

                grid%gcconv_full_n1(i, j) = State_Diag%BudgetConvectionFull(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n1)%Info%AdvectID)
                grid%gcconv_trop_n1(i, j) = State_Diag%BudgetConvectionTrop(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n1)%Info%AdvectID)
                grid%gcconv_pbl_n1(i, j)  = State_Diag%BudgetConvectionPBL(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n1)%Info%AdvectID)

                grid%gcchem_full_n1(i, j) = State_Diag%BudgetChemistryFull(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n1)%Info%AdvectID)
                grid%gcchem_trop_n1(i, j) = State_Diag%BudgetChemistryTrop(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n1)%Info%AdvectID)
                grid%gcchem_pbl_n1(i, j)  = State_Diag%BudgetChemistryPBL(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n1)%Info%AdvectID)


                if(State_Chm%SpcData(config_flags%gc_diagn_spc_n1)%Info%WetDepId .ne. MISSING_INT) then
                    grid%gcwetdep_full_n1(i, j) = State_Diag%BudgetWetDepFull(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n1)%Info%WetDepId)
                    grid%gcwetdep_trop_n1(i, j) = State_Diag%BudgetWetDepTrop(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n1)%Info%WetDepId)
                    grid%gcwetdep_pbl_n1(i, j)  = State_Diag%BudgetWetDepPBL(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n1)%Info%WetDepId)

                    grid%wetlscnvfrc_n1(i, k, j) = State_Diag%WetLossConvFrac(II, JJ, k, &
                                                   State_Chm%SpcData(config_flags%gc_diagn_spc_n1)%Info%WetDepId)
                    grid%wetlscnv_n1(i, k, j) = State_Diag%WetLossConv(II, JJ, k, &
                                                   State_Chm%SpcData(config_flags%gc_diagn_spc_n1)%Info%WetDepId)

                    grid%wetlossls_n1(i, k, j) = State_Diag%WetLossLS(II, JJ, k, &
                                                   State_Chm%SpcData(config_flags%gc_diagn_spc_n1)%Info%WetDepId)
                endif
            endif

            if(config_flags%gc_diagn_spc_n2 .ne. 0 .and. config_flags%gc_diagn_spc_n2 .le. State_Chm%nSpecies) then
                grid%cldcnvflx_n2(i, k, j) = State_Diag%CloudConvFlux(II, JJ, k, config_flags%gc_diagn_spc_n2)

                grid%gcemisdrydep_full_n2(i, j) = State_Diag%BudgetEmisDryDepFull(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n2)%Info%AdvectID)
                grid%gcemisdrydep_trop_n2(i, j) = State_Diag%BudgetEmisDryDepTrop(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n2)%Info%AdvectID)
                grid%gcemisdrydep_pbl_n2(i, j)  = State_Diag%BudgetEmisDryDepPBL(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n2)%Info%AdvectID)

                grid%gcmixing_full_n2(i, j) = State_Diag%BudgetMixingFull(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n2)%Info%AdvectID)
                grid%gcmixing_trop_n2(i, j) = State_Diag%BudgetMixingTrop(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n2)%Info%AdvectID)
                grid%gcmixing_pbl_n2(i, j)  = State_Diag%BudgetMixingPBL(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n2)%Info%AdvectID)

                grid%gcconv_full_n2(i, j) = State_Diag%BudgetConvectionFull(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n2)%Info%AdvectID)
                grid%gcconv_trop_n2(i, j) = State_Diag%BudgetConvectionTrop(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n2)%Info%AdvectID)
                grid%gcconv_pbl_n2(i, j)  = State_Diag%BudgetConvectionPBL(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n2)%Info%AdvectID)

                grid%gcchem_full_n2(i, j) = State_Diag%BudgetChemistryFull(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n2)%Info%AdvectID)
                grid%gcchem_trop_n2(i, j) = State_Diag%BudgetChemistryTrop(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n2)%Info%AdvectID)
                grid%gcchem_pbl_n2(i, j)  = State_Diag%BudgetChemistryPBL(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n2)%Info%AdvectID)

                if(State_Chm%SpcData(config_flags%gc_diagn_spc_n2)%Info%WetDepId .ne. MISSING_INT) then
                    grid%gcwetdep_full_n2(i, j) = State_Diag%BudgetWetDepFull(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n2)%Info%WetDepId)
                    grid%gcwetdep_trop_n2(i, j) = State_Diag%BudgetWetDepTrop(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n2)%Info%WetDepId)
                    grid%gcwetdep_pbl_n2(i, j)  = State_Diag%BudgetWetDepPBL(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n2)%Info%WetDepId)

                    grid%wetlscnvfrc_n2(i, k, j) = State_Diag%WetLossConvFrac(II, JJ, k, &
                                                   State_Chm%SpcData(config_flags%gc_diagn_spc_n2)%Info%WetDepId)
                    grid%wetlscnv_n2(i, k, j) = State_Diag%WetLossConv(II, JJ, k, &
                                                   State_Chm%SpcData(config_flags%gc_diagn_spc_n2)%Info%WetDepId)

                    grid%wetlossls_n2(i, k, j) = State_Diag%WetLossLS(II, JJ, k, &
                                                   State_Chm%SpcData(config_flags%gc_diagn_spc_n2)%Info%WetDepId)
                endif
            endif

            if(config_flags%gc_diagn_spc_n3 .ne. 0 .and. config_flags%gc_diagn_spc_n3 .le. State_Chm%nSpecies) then
                grid%cldcnvflx_n3(i, k, j) = State_Diag%CloudConvFlux(II, JJ, k, config_flags%gc_diagn_spc_n3)

                grid%gcemisdrydep_full_n3(i, j) = State_Diag%BudgetEmisDryDepFull(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n3)%Info%AdvectID)
                grid%gcemisdrydep_trop_n3(i, j) = State_Diag%BudgetEmisDryDepTrop(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n3)%Info%AdvectID)
                grid%gcemisdrydep_pbl_n3(i, j)  = State_Diag%BudgetEmisDryDepPBL(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n3)%Info%AdvectID)

                grid%gcmixing_full_n3(i, j) = State_Diag%BudgetMixingFull(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n3)%Info%AdvectID)
                grid%gcmixing_trop_n3(i, j) = State_Diag%BudgetMixingTrop(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n3)%Info%AdvectID)
                grid%gcmixing_pbl_n3(i, j)  = State_Diag%BudgetMixingPBL(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n3)%Info%AdvectID)

                grid%gcconv_full_n3(i, j) = State_Diag%BudgetConvectionFull(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n3)%Info%AdvectID)
                grid%gcconv_trop_n3(i, j) = State_Diag%BudgetConvectionTrop(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n3)%Info%AdvectID)
                grid%gcconv_pbl_n3(i, j)  = State_Diag%BudgetConvectionPBL(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n3)%Info%AdvectID)

                grid%gcchem_full_n3(i, j) = State_Diag%BudgetChemistryFull(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n3)%Info%AdvectID)
                grid%gcchem_trop_n3(i, j) = State_Diag%BudgetChemistryTrop(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n3)%Info%AdvectID)
                grid%gcchem_pbl_n3(i, j)  = State_Diag%BudgetChemistryPBL(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n3)%Info%AdvectID)

                if(State_Chm%SpcData(config_flags%gc_diagn_spc_n3)%Info%WetDepId .ne. MISSING_INT) then
                    grid%gcwetdep_full_n3(i, j) = State_Diag%BudgetWetDepFull(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n3)%Info%WetDepId)
                    grid%gcwetdep_trop_n3(i, j) = State_Diag%BudgetWetDepTrop(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n3)%Info%WetDepId)
                    grid%gcwetdep_pbl_n3(i, j)  = State_Diag%BudgetWetDepPBL(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n3)%Info%WetDepId)

                    grid%wetlscnvfrc_n3(i, k, j) = State_Diag%WetLossConvFrac(II, JJ, k, &
                                                   State_Chm%SpcData(config_flags%gc_diagn_spc_n3)%Info%WetDepId)
                    grid%wetlscnv_n3(i, k, j) = State_Diag%WetLossConv(II, JJ, k, &
                                                   State_Chm%SpcData(config_flags%gc_diagn_spc_n3)%Info%WetDepId)

                    grid%wetlossls_n3(i, k, j) = State_Diag%WetLossLS(II, JJ, k, &
                                                   State_Chm%SpcData(config_flags%gc_diagn_spc_n3)%Info%WetDepId)
                endif
            endif
        enddo
        enddo
        enddo


        if(debug_level .ge. 100) then
            WRITE(6, *) "%%%% REGISTERED STATE_CHM SPECIES DUMP %%%%"
            do N = 1, MIN(State_Chm%nSpecies, 15)
                WRITE(6, *) "N:", N, "Name:", State_Chm%SpcData(N)%Info%Name
                WRITE(6, *) "Value at PET (1,1,1):", State_Chm%Species(1, 1, 1, N)
            enddo
        endif

        if(debug_level .ge. 1 .and. debug_level .lt. 100) then
            WRITE(6, *) "%%%% State_Chm IMPORTANT SPECIES DUMP - DEBUG WRF|GCHP %%%%"
            WRITE(6, *) "At I = 1, J = 1, L = 1, units: ", State_Chm%Spc_Units
            WRITE(6, *) "O3:", State_Chm%Species(1, 1, 1, IND_('O3'))
            WRITE(6, *) "CO:", State_Chm%Species(1, 1, 1, IND_('CO'))
            WRITE(6, *) "NH3:", State_Chm%Species(1, 1, 1, IND_('NH3'))
            WRITE(6, *) "NO2:", State_Chm%Species(1, 1, 1, IND_('NO2'))

            WRITE(6, *) "Sums at L = 1"
            WRITE(6, *) "O3:", SUM(State_Chm%Species(:, :, 1, IND_('O3')))
            WRITE(6, *) "CO:", SUM(State_Chm%Species(:, :, 1, IND_('CO')))
            WRITE(6, *) "NH3:", SUM(State_Chm%Species(:, :, 1, IND_('NH3')))

            WRITE(6, *) "%%%% State_Diag DUMP - DEBUG WRF|GC %%%%"
            WRITE(6, *) "At I = 1, J = 1 (L = 1) - diagnostic states"
            WRITE(6, *) "Species#", config_flags%gc_diagn_spc_n0, " ",        &
                                    config_flags%gc_diagn_spc_n1, " ",        &
                                    config_flags%gc_diagn_spc_n2, " ",        &
                                    config_flags%gc_diagn_spc_n3
            WRITE(6, *) "Species Names -"
            IF(config_flags%gc_diagn_spc_n0 .ne. 0 .and. config_flags%gc_diagn_spc_n0 .le. State_Chm%nSpecies) THEN
                WRITE(6, *) "0: ", State_Chm%SpcData(config_flags%gc_diagn_spc_n0)%Info%Name
            ENDIF
            IF(config_flags%gc_diagn_spc_n1 .ne. 0 .and. config_flags%gc_diagn_spc_n1 .le. State_Chm%nSpecies) THEN
                WRITE(6, *) "1: ", State_Chm%SpcData(config_flags%gc_diagn_spc_n1)%Info%Name
            ENDIF
            IF(config_flags%gc_diagn_spc_n2 .ne. 0 .and. config_flags%gc_diagn_spc_n2 .le. State_Chm%nSpecies) THEN
                WRITE(6, *) "2: ", State_Chm%SpcData(config_flags%gc_diagn_spc_n2)%Info%Name
            ENDIF
            IF(config_flags%gc_diagn_spc_n3 .ne. 0 .and. config_flags%gc_diagn_spc_n3 .le. State_Chm%nSpecies) THEN
                WRITE(6, *) "3: ", State_Chm%SpcData(config_flags%gc_diagn_spc_n3)%Info%Name
            ENDIF
            WRITE(6, *) "Budget Diags:", State_Diag%Archive_Budget,           &
                                         State_Diag%Archive_BudgetEmisDryDep, &
                                         State_Diag%Archive_BudgetMixing,     &
                                         State_Diag%Archive_BudgetConvection, &
                                         State_Diag%Archive_BudgetChemistry,  &
                                         State_Diag%Archive_BudgetWetDep
            WRITE(6, *) "(Glo|EmisDryDep|Mix|Conv|Chem|WetDep)"
            WRITE(6, *) "Values - FULL for N0", config_flags%gc_diagn_spc_n0
            WRITE(6, *) State_Diag%BudgetEmisDryDepFull(1, 1, config_flags%gc_diagn_spc_n0)
            WRITE(6, *) State_Diag%BudgetMixingFull(1, 1, config_flags%gc_diagn_spc_n0)
            WRITE(6, *) State_Diag%BudgetConvectionFull(1, 1, config_flags%gc_diagn_spc_n0)
            WRITE(6, *) State_Diag%BudgetChemistryFull(1, 1, config_flags%gc_diagn_spc_n0)
            WRITE(6, *) State_Diag%BudgetWetDepFull(1, 1, config_flags%gc_diagn_spc_n0)
            WRITE(6, *) "-- end of wrf|gc temporary diagnostics debug --"
        endif

    end subroutine GIGC_Set_WRF
    subroutine GIGC_Chunk_Cleanup(config_flags, kts, kte, Input_Opt, State_Met, State_Chm)
        implicit none

        type(grid_config_rec_type), intent(in) :: config_flags
        integer, intent(in) :: kts, kte
        type(OptInput), intent(inout) :: Input_Opt
        type(MetState), intent(inout) :: State_Met
        type(ChmState), intent(inout) :: State_Chm



    end subroutine GIGC_Chunk_Cleanup

end module GIGC_Convert_State_Mod
