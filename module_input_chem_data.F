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
!  MODULE: module_input_chem_data
!  DESCRIPTION: Input Chemistry Data Module for "Pumpkin" Abstraction Layer
!               Satisfying all WRF Calls through stubbing or redirection to external
!               parameters, as of WRF v3.9.1.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module module_input_chem_data
   use module_io_domain
   use module_domain
   use module_data_sorgam, ONLY : conmin, rgasuniv, epsilc, grav
   use module_get_file_names, only: eligible_file_name, number_of_eligible_files, unix_ls

   implicit none

   ! last_chem_time
   ! req. by chem/chem_driver for timestep calculations.
   type(WRFU_Time), dimension(max_domains) :: last_chem_time

   ! loop index
   integer :: k_loop

   ! lo
   ! Number of chemicals in initial profile
   integer :: lo

   ! logg
   ! Number of final chemical species
   integer :: logg

   ! kx
   ! Number of vertical levels in temp profile
   integer :: kx

   ! kxm1
   integer :: kxm1

   ! Values for kx, lo, logg from WRF-Chem 3.9.1
   parameter(kx=16, kxm1=kx-1, lo=34, logg=350)

   ! zfa
   ! profile heights in meters, req. by input_chem_profile
   real, dimension(kx + 1) :: zfa
   data ZFA/0., 85., 212., 385., 603., 960., 1430., 2010., &
      2850., 4010., 5340., 6900., 8510., 10200., 12100., 16000., &
      21000./

   ! dens
   real, dimension(kx) :: dens
   data dens/2.738E+18, 5.220E+18, 7.427E+18, 9.202E+18, &
      1.109E+19, 1.313E+19, 1.525E+19, 1.736E+19, &
      1.926E+19, 2.074E+19, 2.188E+19, 2.279E+19, &
      2.342E+19, 2.384E+19, 2.414E+19, 2.434E+19/

   ! iref
   ! The reference index.
   integer, dimension(logg) :: iref

   ! fracref
   ! Reference fraction corresponding to iref.
   ! e.g. for WRF-Chem RADM, species # 1 is SO2. iref(1) = 12, xl(12, k) = so2 profile.
   ! this is converted from the NALROM chemistry model, according to wrf-chem 3.9.1 source.
   real, dimension(logg) :: fracref

   ! xl
   ! req. by input_chem_profile
   real, dimension(lo, kx) :: xl

   !wig: To match the xl profile to the correct species, match WRF's p_<species>
   !     flag with iref(p_<species>-1) to get the value of the first index in xl,
   !     e.g. p_o3=6, iref(6-1)=1, so xl(1,:) is the ozone profile.
   !     See gasprofile_init_pnnl for an explination of what height
   !     each index represents.
   DATA(xl(1, k_loop), k_loop=1, kx) &
      /1.68E-07, 1.68E-07, 5.79E-08, 5.24E-08, 5.26E-08, &
      5.16E-08, 4.83E-08, 4.50E-08, 4.16E-08, 3.80E-08, 3.56E-08, &
      3.35E-08, 3.15E-08, 3.08E-08, 3.06E-08, 3.00E-08/

   DATA(xl(2, k_loop), k_loop=1, kx) &
      /4.06E-10, 4.06E-10, 2.16E-10, 1.37E-10, 9.47E-11, &
      6.95E-11, 5.31E-11, 4.19E-11, 3.46E-11, 3.01E-11, 2.71E-11, &
      2.50E-11, 2.35E-11, 2.26E-11, 2.20E-11, 2.16E-11/

   DATA(xl(3, k_loop), k_loop=1, kx) &
      /9.84E-10, 9.84E-10, 5.66E-10, 4.24E-10, 3.26E-10, &
      2.06E-10, 1.12E-10, 7.33E-11, 7.03E-11, 7.52E-11, 7.96E-11, &
      7.56E-11, 7.27E-11, 7.07E-11, 7.00E-11, 7.00E-11/

   DATA(xl(4, k_loop), k_loop=1, kx) &
      /8.15E-10, 8.15E-10, 8.15E-10, 8.15E-10, 8.15E-10, &
      8.65E-10, 1.07E-09, 1.35E-09, 1.47E-09, 1.47E-09, 1.47E-09, &
      1.47E-09, 1.45E-09, 1.43E-09, 1.40E-09, 1.38E-09/

   DATA(xl(5, k_loop), k_loop=1, kx) &
      /4.16E-10, 4.16E-10, 4.16E-10, 4.16E-10, 4.16E-10, &
      4.46E-10, 5.57E-10, 1.11E-09, 1.63E-09, 1.63E-09, 1.63E-09, &
      1.63E-09, 1.61E-09, 1.59E-09, 1.57E-09, 1.54E-09/

   !  CO is 70 ppbv at top, 80 throughout troposphere
   DATA(xl(6, k_loop), k_loop=1, kx)/7.00E-08, kxm1*8.00E-08/

   DATA(xl(7, k_loop), k_loop=1, kx) &
      /8.33E-29, 8.33E-29, 8.33E-29, 8.33E-29, 8.33E-29, &
      1.33E-28, 3.54E-28, 1.85E-28, 1.29E-29, 1.03E-30, 1.72E-31, &
      7.56E-32, 1.22E-31, 2.14E-31, 2.76E-31, 2.88E-31/

   DATA(xl(8, k_loop), k_loop=1, kx) &
      /9.17E-11, 9.17E-11, 9.17E-11, 9.17E-11, 9.17E-11, &
      1.03E-10, 1.55E-10, 2.68E-10, 4.47E-10, 4.59E-10, 4.72E-10, &
      4.91E-10, 5.05E-10, 5.13E-10, 5.14E-10, 5.11E-10/
   DATA(xl(9, k_loop), k_loop=1, kx) &
      /7.10E-12, 7.10E-12, 7.10E-12, 7.10E-12, 7.10E-12, &
      7.36E-12, 1.02E-11, 2.03E-11, 2.98E-11, 3.01E-11, 3.05E-11, &
      3.08E-11, 3.08E-11, 3.06E-11, 3.03E-11, 2.99E-11/
   DATA(xl(10, k_loop), k_loop=1, kx) &
      /4.00E-11, 4.00E-11, 4.00E-11, 3.27E-11, 2.51E-11, &
      2.61E-11, 2.20E-11, 1.69E-11, 1.60E-11, 1.47E-11, 1.37E-11, &
      1.30E-11, 1.24E-11, 1.20E-11, 1.18E-11, 1.17E-11/
   DATA(xl(11, k_loop), k_loop=1, kx) &
      /1.15E-16, 1.15E-16, 2.46E-15, 2.30E-14, 1.38E-13, &
      6.25E-13, 2.31E-12, 7.32E-12, 1.87E-11, 3.68E-11, 6.10E-11, &
      9.05E-11, 1.22E-10, 1.50E-10, 1.70E-10, 1.85E-10/
   DATA(xl(12, k_loop), k_loop=1, kx) &
      /1.00E-10, 1.00E-10, 1.00E-10, 1.00E-10, 1.00E-10, &
      1.00E-10, 1.00E-10, 1.00E-10, 1.00E-10, 1.00E-10, 1.00E-10, &
      1.00E-10, 1.00E-10, 1.00E-10, 1.00E-10, 1.00E-10/
   DATA(xl(13, k_loop), k_loop=1, kx) &
      /1.26E-11, 1.26E-11, 2.02E-11, 2.50E-11, 3.02E-11, &
      4.28E-11, 6.62E-11, 1.08E-10, 1.54E-10, 2.15E-10, 2.67E-10, &
      3.24E-10, 3.67E-10, 3.97E-10, 4.16E-10, 4.31E-10/
   DATA(xl(14, k_loop), k_loop=1, kx) &
      /1.15E-16, 1.15E-16, 2.46E-15, 2.30E-14, 1.38E-13, &
      6.25E-13, 2.31E-12, 7.32E-12, 1.87E-11, 3.68E-11, 6.10E-11, &
      9.05E-11, 1.22E-10, 1.50E-10, 1.70E-10, 1.85E-10/
   DATA(xl(15, k_loop), k_loop=1, kx) &
      /1.00E-20, 1.00E-20, 6.18E-20, 4.18E-18, 1.23E-16, &
      2.13E-15, 2.50E-14, 2.21E-13, 1.30E-12, 4.66E-12, 1.21E-11, &
      2.54E-11, 4.47E-11, 6.63E-11, 8.37E-11, 9.76E-11/
   DATA(xl(16, k_loop), k_loop=1, kx) &
      /1.23E-11, 1.23E-11, 1.23E-11, 1.23E-11, 1.23E-11, &
      1.20E-11, 9.43E-12, 3.97E-12, 1.19E-12, 1.11E-12, 9.93E-13, &
      8.66E-13, 7.78E-13, 7.26E-13, 7.04E-13, 6.88E-13/
   DATA(xl(17, k_loop), k_loop=1, kx) &
      /1.43E-12, 1.43E-12, 1.43E-12, 1.43E-12, 1.43E-12, &
      1.50E-12, 2.64E-12, 8.90E-12, 1.29E-11, 1.30E-11, 1.32E-11, &
      1.32E-11, 1.31E-11, 1.30E-11, 1.29E-11, 1.27E-11/
   DATA(xl(18, k_loop), k_loop=1, kx) &
      /3.61E-13, 3.61E-13, 3.61E-13, 3.61E-13, 3.61E-13, &
      3.58E-13, 5.22E-13, 1.75E-12, 2.59E-12, 2.62E-12, 2.64E-12, &
      2.66E-12, 2.65E-12, 2.62E-12, 2.60E-12, 2.57E-12/
   DATA(xl(19, k_loop), k_loop=1, kx) &
      /5.00E-11, 5.00E-11, 5.00E-11, 5.00E-11, 5.00E-11, &
      5.00E-11, 5.00E-11, 5.00E-11, 5.00E-11, 5.00E-11, 5.00E-11, &
      5.00E-11, 5.00E-11, 5.00E-11, 5.00E-11, 5.00E-11/

   DATA(xl(20, k_loop), k_loop=1, kx)/kx*1.E-20/
   DATA(xl(21, k_loop), k_loop=1, kx)/kx*1.E-20/
   DATA(xl(22, k_loop), k_loop=1, kx)/kx*1.E-20/
   DATA(xl(23, k_loop), k_loop=1, kx)/kx*1.E-20/
   DATA(xl(24, k_loop), k_loop=1, kx)/kx*1.E-20/
   DATA(xl(25, k_loop), k_loop=1, kx)/kx*1.E-20/

   ! Propane - Gregory PEM-West A 25 ppt median marine boundary layer
   DATA(xl(26, k_loop), k_loop=1, kx) &
      /5.00E-13, 1.24E-12, 2.21E-12, 3.27E-12, 4.71E-12, &
      6.64E-12, 9.06E-12, 1.19E-11, 1.47E-11, 1.72E-11, &
      1.93E-11, 2.11E-11, 2.24E-11, 2.34E-11, 2.42E-11, 2.48E-11/
   ! Acetylene - Gregory PEM-West A 53 ppt median marine boundary layer
   DATA(xl(27, k_loop), k_loop=1, kx) &
      /1.00E-12, 2.48E-12, 4.42E-12, 6.53E-12, 9.42E-12, &
      1.33E-11, 1.81E-11, 2.37E-11, 2.95E-11, 3.44E-11, &
      3.85E-11, 4.22E-11, 4.49E-11, 4.69E-11, 4.84E-11, 4.95E-11/
   ! OH
   DATA(xl(28, k_loop), k_loop=1, kx) &
      /9.80E+06, 9.80E+06, 4.89E+06, 2.42E+06, 1.37E+06, &
      9.18E+05, 7.29E+05, 6.26E+05, 5.01E+05, 4.33E+05, 4.05E+05, &
      3.27E+05, 2.54E+05, 2.03E+05, 1.74E+05, 1.52E+05/
   ! HO2
   DATA(xl(29, k_loop), k_loop=1, kx) &
      /5.74E+07, 5.74E+07, 7.42E+07, 8.38E+07, 8.87E+07, &
      9.76E+07, 1.15E+08, 1.34E+08, 1.46E+08, 1.44E+08, 1.40E+08, &
      1.36E+08, 1.31E+08, 1.28E+08, 1.26E+08, 1.26E+08/
   ! NO3+N2O5
   DATA(xl(30, k_loop), k_loop=1, kx) &
      /5.52E+05, 5.52E+05, 3.04E+05, 2.68E+05, 2.32E+05, &
      1.66E+05, 1.57E+05, 1.72E+05, 1.98E+05, 2.22E+05, 2.43E+05, &
      2.75E+05, 3.00E+05, 3.18E+05, 3.32E+05, 3.39E+05/
   ! HO2NO2
   DATA(xl(31, k_loop), k_loop=1, kx) &
      /7.25E+07, 7.25E+07, 6.36E+07, 5.55E+07, 4.94E+07, &
      3.66E+07, 2.01E+07, 9.57E+06, 4.75E+06, 2.37E+06, 1.62E+06, &
      9.86E+05, 7.05E+05, 5.63E+05, 4.86E+05, 4.41E+05/
   ! Sum of RO2 &
   DATA(xl(32, k_loop), k_loop=1, kx) &
      /9.14E+06, 9.14E+06, 1.46E+07, 2.14E+07, 2.76E+07, &
      3.62E+07, 5.47E+07, 1.19E+08, 2.05E+08, 2.25E+08, 2.39E+08, &
      2.58E+08, 2.82E+08, 2.99E+08, 3.08E+08, 3.15E+08/
   ! O3 <--This is not the O3 used for RADM2 or CBMZ (wig)
   DATA(xl(33, k_loop), k_loop=1, kx) &
      /8.36E+11, 8.36E+11, 4.26E+11, 4.96E+11, 6.05E+11, &
      6.93E+11, 7.40E+11, 7.74E+11, 7.82E+11, 7.75E+11, 7.69E+11, &
      7.59E+11, 7.54E+11, 7.50E+11, 7.47E+11, 7.45E+11/
   ! NOx (NO+NO2)
   DATA(xl(34, k_loop), k_loop=1, kx) &
      /1.94E+09, 1.94E+09, 1.53E+09, 1.24E+09, 1.04E+09, &
      8.96E+08, 7.94E+08, 7.11E+08, 6.44E+08, 6.00E+08, 5.70E+08, &
      5.49E+08, 5.35E+08, 5.28E+08, 5.24E+08, 5.23E+08/

contains
   ! get_last_gas
   ! Get the index of the last gas species depending on mechanism. req. by chem/chem_driver

   ! Finished, note that gas registry is:
   ! a3o2,acet,acta,ald2,alk4,ato2,atooh,b3o2,benzene,br,br2,brcl,brno2,brno3,bro,bro2,brsala,brsalc,c2h6,c3h8,c4hvp1,c4hvp2,ccl4,cfc11,cfc113,cfc114,cfc115,cfc12,ch2br2,ch2cl2,ch2i2,ch2ibr,ch2icl,ch2o,ch2oo,ch3br,ch3ccl3,ch3choo,ch3cl,ch3i,ch4,chbr3,chcl3,cl,cl2,cl2o2,clno2,clno3,clo,cloo,co,co2,eoh,ethln,etno3,eto2,etp,glyc,glyx,h,h1211,h1301,h2,h2402,h2o,h2o2,hac,hbr,hc5a,hcfc123,hcfc141b,hcfc142b,hcfc22,hcl,hcooh,hi,hmhp,hmml,hno2,hno3,hno4,ho2,hobr,hocl,hoi,honit,hpald1,hpald1oo,hpald2,hpald2oo,hpald3,hpald4,hpethnl,i,i2,i2o2,i2o3,i2o4,ibr,iche,ichoo,icl,icn,icnoo,icpdh,idc,idchp,idhdp,idhnboo,idhndoo1,idhndoo2,idhpe,idn,idnoo,iepoxa,iepoxaoo,iepoxb,iepoxboo,iepoxd,ihn1,ihn2,ihn3,ihn4,ihoo1,ihoo4,ihpnboo,ihpndoo,ihpoo1,ihpoo2,ihpoo3,ina,ino,ino2b,ino2d,inpb,inpd,io,iono,iono2,iprno3,isopnoo1,isopnoo2,isoprene,itcn,ithn,ko2,lbro2h,lbro2n,lch4,lco,limo2,limon,lisopno3,lisopoh,lox,ltro2h,ltro2n,lvoc,lvocoa,lxro2h,lxro2n,macr,macr1oo,macr1ooh,macrno2,map,mco3,mcrdh,mcrenol,mcrhn,mcrhnb,mcrhp,mcrohoo,mek,meno3,mgly,mo2,moh,monits,monitu,mp,mpan,mpn,mtpa,mtpo,mvk,mvkdh,mvkhc,mvkhcb,mvkhp,mvkn,mvkohoo,mvkpc,n,n2,n2o,n2o5,nh3,no,no2,no3,nprno3,o,o1d,o2,o3,oclo,ocs,oh,oio,olnd,olnn,othro2,pan,pco,pfe,ph2o2,pio2,pip,po2,pox,pp,ppn,prn1,propnn,prpe,prpn,pso4,pyac,r4n1,r4n2,r4o2,r4p,ra3p,rb3p,rcho,rco3,rcooh,ripa,ripb,ripc,ripd,roh,rp,so2,so4h1,so4h2,toluene,tro2,xro2,xylenes
   ! 
   ! The non-gas registry is:
   ! aeri,bcpi,bcpo,dms,dst1,dst2,dst3,dst4,indiol,ionita,isala,isalc,monita,msa,nh4,nit,nits,ocpi,ocpo,sala,salc,so4,so4s,soagx,soaie,soap,soas

   integer function get_last_gas(chem_opt)
      implicit none
      integer, intent(in) :: chem_opt

      select case (chem_opt)
      case (0)
         get_last_gas = 0
      case (1)
         get_last_gas = p_ho2
      case (233) ! GEOS-Chem
         get_last_gas = p_xylenes

      case default
         call wrf_error_fatal("Pumpkin module_input_chem_data::get_last_gas: could not decipher chem_opt value")

      end select

   end function get_last_gas

   ! setup_gasprofile_maps
   ! Sets up the cross reference mapping indices and fractional
   ! apportionment of the default species profiles for use with ICs and BCs.
   ! req. by chem/chemics_init
   !!!! CHEMISTRY DEVELOPERS: YOU MUST UPDATE THIS FOR THE RIGHT INDEXES DEP. ON YOUR REGISTRY.
   subroutine setup_gasprofile_maps(chem_opt, numgas)
      integer, intent(in) :: chem_opt, numgas
      select case (chem_opt)
      case (1)
         !! CHEMISTRY DEVELOPERS: CHEMISTRY NEEDS TO BE ADDED HERE.
         !! Either by including a new module file (good practice) or hardcode it...
         ! call setup_gasprofile_map_geoschem
      end select
   end subroutine setup_gasprofile_maps

   ! an example for setup_gasprofile_maps is set below.
   ! you should preferably add this to a module file!
   ! this is setup_gasprofile_map_radm_racm from WRF-Chem 3.6.1
   ! (c) original authors
   ! subroutine setup_gasprofile_map_radm_racm

   !    iref(:) = 7 !default value
   !    iref(1:41) = (/12, 19, 2, 2, 1, 3, 4, 9, 8, 5, 5, 32, 6, 6, 6, 30, 30, 10, 26, 13, 11, 6, 6, &
   !                   14, 15, 15, 23, 23, 32, 16, 23, 31, 17, 23, 23, 23, 23, 23, 7, 28, 29/)

   !    fracref(:) = 1. !default value
   !    fracref(1:41) = (/1., 1., .75, .25, 1., 1., 1., 1., 1., 1., &
   !                      .5, .5, 6.25E-4, 7.5E-4, 6.25E-5, .1, &
   !                      .9, 1., 1., 1., 1., 8.E-3, 1., 1., 1., .5, &
   !                      1., 1., .5, 1., 1., 1., 1., 1., 1., 1., 1., &
   !                      1., 1., 1., 1./)

   !    ggnam(:) = 'JUNK' !default value
   !    ggnam(1:41) = (/'SO2 ', 'SULF', 'NO2 ', 'NO  ', 'O3  ', 'HNO3', &
   !                    'H2O2', 'ALD ', 'HCHO', 'OP1 ', 'OP2 ', 'PAA ', &
   !                    'ORA1', 'ORA2', 'NH3 ', 'N2O5', 'NO3 ', 'PAN ', &
   !                    'HC3 ', 'HC5 ', 'HC8 ', 'ETH ', 'CO  ', 'OL2 ', &
   !                    'OLT ', 'OLI ', 'TOL ', 'XYL ', 'ACO3', 'TPAN', &
   !                    'HONO', 'HNO4', 'KET ', 'GLY ', 'MGLY', 'DCB ', &
   !                    'ONIT', 'CSL ', 'ISO ', 'HO  ', 'HO2 '/)

   ! end subroutine setup_gasprofile_map_radm_racm

   ! input_chem_profile
   ! req. by real_em
   ! Based on original WRF-Chem module_input_chem_data::input_chem_profile, (c) original authors
   subroutine input_chem_profile(grid)
      implicit none
      type(domain) :: grid

      integer :: i, j, k, &
                 ids, ide, jds, jde, kds, kde, &
                 ims, ime, jms, jme, kms, kme, &
                 ips, ipe, jps, jpe, kps, kpe
      integer :: fid, ierr, numgas
      integer :: debug_level

      REAL, ALLOCATABLE, DIMENSION(:, :, :) :: si_zsigf, si_zsig

      ! Mean gravitational acceleration (m/s^2)
      real grav
      parameter (grav=9.80622)

      ! Get grid dimensions
      call get_ijk_from_grid(grid, &
                             ids, ide, jds, jde, kds, kde, &
                             ims, ime, jms, jme, kms, kme, &
                             ips, ipe, jps, jpe, kps, kpe)

      ! Get scalar grid point heights
      allocate (si_zsigf(ims:ime, kms:kme, jms:jme))
      allocate (si_zsig(ims:ime, kms:kme, jms:jme))

      si_zsigf = (grid%ph_1 + grid%phb)/grav

      do k = 1, kde - 1
         si_zsig(:, k, :) = 0.5*(si_zsigf(:, k, :) + si_zsigf(:, k + 1, :))
      enddo
      si_zsig(:, kde, :) = 0.5*(3.*si_zsigf(:, kde, :) - si_zsigf(:, kde - 1, :))

      ! Determine the index of the last gas species
      numgas = get_last_gas(grid%chem_opt)

      ! Setup the cross reference mappings between the default profiles and
      ! the gas mechanism species (wig, 2-May-2007)
      call setup_gasprofile_maps(grid%chem_opt, numgas)

      ! Interpolate the chemistry data to the grid. These values should typically
      ! be set to match the values in bdy_chem_value_tracer so that the boundaries
      ! and interior match each other.
      if (grid%chem_opt == 0) then
         ! this is set as an example from the original module_input_chem_data
         ! grid%chem(ims:ime, kms:kme, jms:jme, 1:numgas) = 0.
      else
         call make_chem_profile(ims, ime, jms, jme, kms, kme, num_chem, numgas, &
                                grid%chem_opt, si_zsig, grid%chem)
      end if

      call wrf_debug(100, ' input_chem_profile: exit subroutine ')

      deallocate (si_zsigf); deallocate (si_zsig)
      return
   end subroutine input_chem_profile

   ! make_chem_profile
   ! req. by input_chem_profile, in turn required by real_em
   ! Based on original WRF-Chem module_input_chem_data::input_chem_profile, (c) original authors
   subroutine make_chem_profile(nx1, nx2, ny1, ny2, nz1, nz2, nch, numgas, &
                                chem_opt, zgrid, chem)
      implicit none

      integer, intent(in) :: nx1, ny1, nz1, nx2, ny2, nz2
      integer, intent(in) :: nch, numgas, chem_opt
      real, dimension(nx1:nx2, nz1:nz2, ny1:ny2) :: zgrid

      integer :: i, j, k, l, is
      real, dimension(nx1:nx2, nz1:kx, ny1:ny2, lo + 1) :: chprof
      real, dimension(nx1:nx2, nz1:kx, ny1:ny2) :: zprof

      real, dimension(nx1:nx2, nz1:nz2, ny1:ny2, nch) :: chem
      real, dimension(nx1:nx2, nz1:nz2, ny1:ny2, lo) :: stor

      if (nch .NE. num_chem) then
         call wrf_error_fatal("Pumpkin module_input_chem_data: wrong number of chemical species.")
      endif

      ! Vertically flip the chemistry data as it is given top down and
      ! heights are bottom up. Fill temp 3D chemical and profile array,
      ! keep chem slot 1 open as vinterp_chem assumes there is no data.
      do j = ny1, ny2
         do k = 1, kx
            do i = nx1, nx2
               chprof(i, k, j, 2:lo + 1) = xl(1:lo, kx - k + 1)
               zprof(i, k, j) = 0.5*(zfa(k) + zfa(k + 1))
            enddo
         enddo
      enddo

      ! Change number concentrations to mixing ratios for short-lived NALROM species
      do k = 1, kx
         chprof(:, k, :, lo - 5:lo + 1) = chprof(:, k, :, lo - 5:lo + 1)/dens(k)
      end do

      ! Interpolate temp 3D chemical and profile array to WRF grid
      call vinterp_chem(nx1, nx2, ny1, ny2, nz1, kx, nz2, lo, zprof, zgrid, &
                        chprof, chem, .false.)

      ! place interpolated data into temp storage array
      stor(nx1:nx2, nz1:nz2, ny1:ny2, 1:lo) = chem(nx1:nx2, nz1:nz2, ny1:ny2, 2:lo + 1)

      ! Here is where the chemistry profile is constructed
      chem(nx1:nx2, nz1:nz2, ny1:ny2, 1) = -999.

      do l = 2, numgas
         is = iref(l - 1)
         do j = ny1, ny2
            do k = nz1, nz2
               do i = nx1, nx2
                  chem(i, k, j, l) = fracref(l - 1) * stor(i, k, j, is) * 1.E6
               enddo
            enddo
         enddo
      enddo

      return
   end subroutine make_chem_profile

   ! vinterp_chem
   ! Interpolates columns of chemistry data from one set of height surfaces to another.
   ! req. by make_chem_profile
   ! WRF-Chem v3.6.1, (c) original authors
   subroutine vinterp_chem(nx1, nx2, ny1, ny2, nz1, nz_in, nz_out, nch, z_in, z_out, &
                           data_in, data_out, extrapolate)
      integer, intent(in)                :: nx1, nx2
      integer, intent(in)                :: ny1, ny2
      integer, intent(in)                :: nz1
      integer, intent(in)                :: nz_in
      integer, intent(in)                :: nz_out
      integer, intent(in)                :: nch
      real, intent(in)                   :: z_in(nx1:nx2, nz1:nz_in, ny1:ny2)
      real, intent(in)                   :: z_out(nx1:nx2, nz1:nz_out, ny1:ny2)
      real, intent(in)                   :: data_in(nx1:nx2, nz1:nz_in, ny1:ny2, nch)
      real, intent(out)                  :: data_out(nx1:nx2, nz1:nz_out, ny1:ny2, nch)
      logical, intent(in)                :: extrapolate

      integer                            :: i, j, l
      integer                            :: k, kk
      real                               :: desired_z
      real                               :: dvaldz
      real                               :: wgt0

      ! Loop over the number of chemical species
      chem_loop: DO l = 2, nch
         data_out(:, :, :, l) = -99999.9
         DO j = ny1, ny2
            DO i = nx1, nx2
               output_loop: DO k = nz1, nz_out
                  desired_z = z_out(i, k, j)
                  IF (desired_z .LT. z_in(i, 1, j)) THEN
                     IF ((desired_z - z_in(i, 1, j)) .LT. 0.0001) THEN
                        data_out(i, k, j, l) = data_in(i, 1, j, l)
                     ELSE
                        IF (extrapolate) THEN
                           ! Extrapolate downward because desired height level is below
                           ! the lowest level in our input data.  Extrapolate using simple
                           ! 1st derivative of value with respect to height for the bottom 2
                           ! input layers.

                           ! Add a check to make sure we are not using the gradient of
                           ! a very thin layer
                           IF ((z_in(i, 1, j) - z_in(i, 2, j)) .GT. 0.001) THEN
                              dvaldz = (data_in(i, 1, j, l) - data_in(i, 2, j, l))/ &
                                       (z_in(i, 1, j) - z_in(i, 2, j))
                           ELSE
                              dvaldz = (data_in(i, 1, j, l) - data_in(i, 3, j, l))/ &
                                       (z_in(i, 1, j) - z_in(i, 3, j))
                           ENDIF
                           data_out(i, k, j, l) = MAX(data_in(i, 1, j, l) + &
                                                      dvaldz*(desired_z - z_in(i, 1, j)), 0.)
                        ELSE
                           data_out(i, k, j, l) = data_in(i, 1, j, l)
                        ENDIF
                     ENDIF
                  ELSE IF (desired_z .GT. z_in(i, nz_in, j)) THEN
                     IF ((z_in(i, nz_in, j) - desired_z) .LT. 0.0001) THEN
                        data_out(i, k, j, l) = data_in(i, nz_in, j, l)
                     ELSE
                        IF (extrapolate) THEN
                           ! Extrapolate upward
                           IF ((z_in(i, nz_in - 1, j) - z_in(i, nz_in, j)) .GT. 0.0005) THEN
                              dvaldz = (data_in(i, nz_in, j, l) - data_in(i, nz_in - 1, j, l))/ &
                                       (z_in(i, nz_in, j) - z_in(i, nz_in - 1, j))
                           ELSE
                              dvaldz = (data_in(i, nz_in, j, l) - data_in(i, nz_in - 2, j, l))/ &
                                       (z_in(i, nz_in, j) - z_in(i, nz_in - 2, j))
                           ENDIF
                           data_out(i, k, j, l) = MAX(data_in(i, nz_in, j, l) + &
                                                      dvaldz*(desired_z - z_in(i, nz_in, j)), 0.)
                        ELSE
                           data_out(i, k, j, l) = data_in(i, nz_in, j, l)
                        ENDIF
                     ENDIF
                  ELSE
                     ! We can trap between two levels and linearly interpolate
                     input_loop: DO kk = 1, nz_in - 1
                        IF (desired_z .EQ. z_in(i, kk, j)) THEN
                           data_out(i, k, j, l) = data_in(i, kk, j, l)
                           EXIT input_loop
                        ELSE IF (desired_z .EQ. z_in(i, kk + 1, j)) THEN
                           data_out(i, k, j, l) = data_in(i, kk + 1, j, l)
                           EXIT input_loop
                        ELSE IF ((desired_z .GT. z_in(i, kk, j)) .AND. &
                                 (desired_z .LT. z_in(i, kk + 1, j))) THEN
                           wgt0 = (desired_z - z_in(i, kk + 1, j))/ &
                                  (z_in(i, kk, j) - z_in(i, kk + 1, j))
                           data_out(i, k, j, l) = MAX(wgt0*data_in(i, kk, j, l) + &
                                                      (1.-wgt0)*data_in(i, kk + 1, j, l), 0.)
                           EXIT input_loop
                        ENDIF
                     ENDDO input_loop
                  ENDIF
               ENDDO output_loop
            ENDDO
         ENDDO
      ENDDO chem_loop

      return
   end subroutine vinterp_chem

   ! flow_dep_bdy_chem
   ! sets zero gradient conditions for outflow and a set profile value
   ! for inflow in the boundary specified region. Note that field must be unstaggered.
   ! The velocities, u and v, will only be used to check their sign (coupled vels OK)
   ! spec_zone is the width of the outer specified b.c.s that are set here.
   !
   ! req. by solve_em
   subroutine flow_dep_bdy_chem(chem, &
                                chem_bxs, chem_btxs, &
                                chem_bxe, chem_btxe, &
                                chem_bys, chem_btys, &
                                chem_bye, chem_btye, &
                                dt, spec_bdy_width, z, have_bcs_chem, &
                                u, v, config_flags, alt, &
                                t, pb, p, t0, p1000mb, rcp, ph, phb, g, &
                                spec_zone, ic, &
                                ids, ide, jds, jde, kds, kde, &
                                ims, ime, jms, jme, kms, kme, &
                                ips, ipe, jps, jpe, kps, kpe, &
                                its, ite, jts, jte, kts, kte)
      implicit none
      INTEGER, INTENT(IN)    :: ids, ide, jds, jde, kds, kde
      INTEGER, INTENT(IN)    :: ims, ime, jms, jme, kms, kme
      INTEGER, INTENT(IN)    :: ips, ipe, jps, jpe, kps, kpe
      INTEGER, INTENT(IN)    :: its, ite, jts, jte, kts, kte
      INTEGER, INTENT(IN)    :: spec_zone, spec_bdy_width, ic
      REAL, INTENT(IN) :: dt

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT) :: chem
      REAL, DIMENSION(jms:jme, kds:kde, spec_bdy_width), INTENT(IN) :: chem_bxs, chem_bxe, chem_btxs, chem_btxe
      REAL, DIMENSION(ims:ime, kds:kde, spec_bdy_width), INTENT(IN) :: chem_bys, chem_bye, chem_btys, chem_btye
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: z
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: alt
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: u
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: v
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(IN) :: &
         ph, phb, t, pb, p
      real, INTENT(IN) :: g, rcp, t0, p1000mb
      TYPE(grid_config_rec_type) config_flags

      INTEGER    :: i, j, k, numgas
      INTEGER    :: ibs, ibe, jbs, jbe, itf, jtf, ktf
      INTEGER    :: i_inner, j_inner
      INTEGER    :: b_dist
      integer    :: itestbc, i_bdy_method
      real tempfac, convfac
      real       :: chem_bv_def
      logical    :: have_bcs_chem

      ! from module_data_sorgam
      real conmin
      parameter(conmin=1.E-16)

      chem_bv_def = conmin
      numgas = get_last_gas(config_flags%chem_opt)
      itestbc = 0
      if (p_nu0 .gt. 1) itestbc = 1
      ibs = ids
      ibe = ide - 1
      itf = min(ite, ide - 1)
      jbs = jds
      jbe = jde - 1
      jtf = min(jte, jde - 1)
      ktf = kde - 1

      ! Set boundary transfer method.
      ! This varies by mechanism and aerosol model, but we use a simplified version
      ! of the original WRF-Chem code. (hplin, 6/19/19)
      i_bdy_method = 0
      if (have_bcs_chem) then
         i_bdy_method = 6
      endif

      ! Below actual code from module_input_chem_data of original WRF-Chem
      ! (c) original authors

      ! Below statement is from bdy_chem_value_gcm (chem, chem_b, chem_bt, dt, ic)
      ! from WRF-Chem -- this is used for i_bdy_method = 6 when have_bcs_chem
      ! is enabled. (hplin, 6/19/19)
      !
      ! chem = max(epsilc, chem_b + chem_bt * dt)

      if (jts - jbs .lt. spec_zone) then
      ! y-start boundary
         do j = jts, min(jtf, jbs + spec_zone - 1)
            b_dist = j - jbs
            do k = kts, ktf
               do i = max(its, b_dist + ibs), min(itf, ibe - b_dist)
                  i_inner = max(i, ibs + spec_zone)
                  i_inner = min(i_inner, ibe - spec_zone)
                  if (v(i, k, j) .lt. 0.) then
                     chem(i, k, j) = chem(i_inner, k, jbs + spec_zone)
                  else
                     if(i_bdy_method .eq. 6) then
                        chem(i, k, j) = max(epsilc, chem_bys(i, k, 1) + chem_btys(i, k, 1) * dt)
                     else
                        ! See default value above
                        chem(i, k, j) = chem_bv_def
                     endif
                  endif
               enddo
            enddo
         enddo
      endif
      if (jbe - jtf .lt. spec_zone) then
      ! y-end boundary
         do j = max(jts, jbe - spec_zone + 1), jtf
            b_dist = jbe - j
            do k = kts, ktf
               do i = max(its, b_dist + ibs), min(itf, ibe - b_dist)
                  i_inner = max(i, ibs + spec_zone)
                  i_inner = min(i_inner, ibe - spec_zone)
                  if (v(i, k, j + 1) .gt. 0.) then
                     chem(i, k, j) = chem(i_inner, k, jbe - spec_zone)
                  else
                     if(i_bdy_method .eq. 6) then
                        chem(i, k, j) = max(epsilc, chem_bye(i, k, 1) + chem_btye(i, k, 1) * dt)
                     else
                        ! See default value above
                        chem(i, k, j) = chem_bv_def
                     endif
                  endif
               enddo
            enddo
         enddo
      endif

      if (its - ibs .lt. spec_zone) then
      ! x-start boundary
         do i = its, min(itf, ibs + spec_zone - 1)
            b_dist = i - ibs
            do k = kts, ktf
               do j = max(jts, b_dist + jbs + 1), min(jtf, jbe - b_dist - 1)
                  j_inner = max(j, jbs + spec_zone)
                  j_inner = min(j_inner, jbe - spec_zone)
                  if (u(i, k, j) .lt. 0.) then
                     chem(i, k, j) = chem(ibs + spec_zone, k, j_inner)
                  else
                     if(i_bdy_method .eq. 6) then
                        chem(i, k, j) = max(epsilc, chem_bxs(j, k, 1) + chem_btxs(j, k, 1) * dt)
                     else
                        ! See default value above
                        chem(i, k, j) = chem_bv_def
                     endif
                  endif
               enddo
            enddo
         enddo
      endif

      if (ibe - itf .lt. spec_zone) then
      ! x-end boundary
         do i = max(its, ibe - spec_zone + 1), itf
            b_dist = ibe - i
            do k = kts, ktf
               do j = max(jts, b_dist + jbs + 1), min(jtf, jbe - b_dist - 1)
                  j_inner = max(j, jbs + spec_zone)
                  j_inner = min(j_inner, jbe - spec_zone)
                  if (u(i + 1, k, j) .gt. 0.) then
                     chem(i, k, j) = chem(ibe - spec_zone, k, j_inner)
                  else
                     if(i_bdy_method .eq. 6) then
                        chem(i, k, j) = max(epsilc, chem_bxe(j, k, 1) + chem_btxe(j, k, 1) * dt)
                     else
                        ! See default value above
                        chem(i, k, j) = chem_bv_def
                     endif
                  endif
               enddo
            enddo
         enddo
      endif
   end subroutine flow_dep_bdy_chem

   ! chem_dbg
   ! Chemistry Debug function req. by chem/chem_driver
   !!!! CHEMISTRY DEVELOPERS: YOU WILL LIKELY NEED TO UPDATE THIS SPECIES LIST FOR DEBUGGING.
   !!!! THIS IS SEPARATE FROM REGISTRY.CHEM
#ifdef CHEM_DBG_I
   subroutine chem_dbg(i, j, k, dtstep, itimestep, &
                       dz8w, t_phy, p_phy, rho_phy, chem, emis_ant, &
                       ids, ide, jds, jde, kds, kde, &
                       ims, ime, jms, jme, kms, kme, &
                       its, ite, jts, jte, kts, kte, &
                       kemit, &
                       ph_macr, ph_o31d, ph_o33p, ph_no2, ph_no3o2, ph_no3o, ph_hno2, &
                       ph_hno3, ph_hno4, ph_h2o2, ph_ch2or, ph_ch2om, ph_ch3cho, &
                       ph_ch3coch3, ph_ch3coc2h5, ph_hcocho, ph_ch3cocho, &
                       ph_hcochest, ph_ch3o2h, ph_ch3coo2h, ph_ch3ono2, ph_hcochob, ph_n2o5, &
                       ph_o2)

      implicit none
      integer, intent(in) :: i, j, k, &
                             ids, ide, jds, jde, kds, kde, &
                             ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte, &
                             kemit
      real, intent(in) :: dtstep
      integer, intent(in) :: itimestep
      real, dimension(ims:ime, kms:kme, jms:jme, num_chem), intent(inout) :: chem
      real, dimension(ims:ime, kms:kme, jms:jme), intent(in) :: dz8w, t_phy, p_phy, rho_phy
      real, dimension(ims:ime, kms:kemit, jms:jme, num_emis_ant), &
         intent(in) :: emis_ant
      real, dimension(ims:ime, kms:kme, jms:jme), &
         intent(in), optional :: &
         ph_macr, ph_o31d, ph_o33p, ph_no2, ph_no3o2, ph_no3o, ph_hno2, &
         ph_hno3, ph_hno4, ph_h2o2, ph_ch2or, ph_ch2om, ph_ch3cho, &
         ph_ch3coch3, ph_ch3coc2h5, ph_hcocho, ph_ch3cocho, &
         ph_hcochest, ph_ch3o2h, ph_ch3coo2h, ph_ch3ono2, ph_hcochob, ph_n2o5, &
         ph_o2

      integer :: n

      print *, "itimestep =", itimestep
      print *, "MET DATA AT (i,k,j):", i, k, j
      print *, "t_phy,p_phy,rho_phy=", t_phy(i, k, j), p_phy(i, k, j), rho_phy(i, k, j)

      print *, "CHEM_DBG PRINT (PPM or ug/m^3) AT (i,k,j):", i, k, j
      do n = 1, num_chem
         print *, n, chem(i, k, j, n)
      end do

      print*
   end subroutine chem_dbg
#endif

end module module_input_chem_data
