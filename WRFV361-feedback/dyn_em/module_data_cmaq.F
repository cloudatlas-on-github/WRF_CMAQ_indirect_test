
MODULE module_data_cmaq


!*******************************************************************                                                                *
!*  parameters and variables for aerosol-cloud interactions        *                                                                 *
!    (7/20/2010, Shaocai Yu)
!*******************************************************************
!
!   maxd_atype = maximum allowable number of aerosol types
!   maxd_asize = maximum allowable number of aerosol size bins
!   maxd_acomp = maximum allowable number of chemical components
!	in each aerosol size bin
!   maxd_aphase = maximum allowable number of aerosol phases (gas, cloud, ice, rain, ...)
!
!   ntype_aer = number of aerosol types = t
!   nsize_aer(t) = number of aerosol size bins for aerosol type t. each bin w/ same set of components
!**********************************************************************
!  For the CMAQ case:
! aitken mode has one set of species, so are treated as isize=1 of itype=1
! accum mode has one set of species, so are treated as isize=1 of itype=2
! coarse mode has one set of species, so is  treated as isize=1 of itype=3
!**********************************************************************
!            ntype_aer = 3
!	     nsize_aer(1)=1  ! nsize_aer(ntype_aer)
!	     nsize_aer(2)=1
!            nsize_aer(3)=1
!   nphase_aer = number of aerosol phases
!
!   msectional - if positive, moving-center sectional code is utilized,
!	and each mode is actually a section.
!   maerosolincw - if positive, both unactivated/interstitial and activated
!	aerosol species are simulated.  if zero/negative, only the
!	unactivated are simulated.
!
!   ncomp_aer(t) = number of chemical components for aerosol type t
!   ncomp_aer_nontracer(t) = number of "non-tracer" chemical components while in gchm code
!   mastercompptr_aer(c,t) = mastercomp type/i.d. for chemical component c
!	(1=sulfate, others to be defined) and aerosol type t.
!   massptr_aer(c,s,t,p) = gchm r-array index for the mixing ratio
!	(moles-x/mole-air) for chemical component c in size bin s for type t and phase p
!
!   waterptr_aer(s,t) = mixing ratio (moles-water/mole-air) for water
!       associated with aerosol size bin s and type t
!   hygroptr_aer(s,t) = gchm r-array index for the bulk hygroscopicity of the size bin and type
!   numptr_aer(s,t,p) = gchm r-array index for the number mixing ratio
!	(particles/mole-air) for aerosol size bin s, type t, and phase p
!       If zero or negative, then number is not being simulated.
!
!   mprognum_aer(s,t,p) - if positive, number mixing-ratio for size s, type t,
!       and phase p will be prognosed.  Otherwise, no.
!
!   ntot_mastercomp_aer = number of aerosol chemical components defined
!   dens_mastercomp_aer(mc) = dry density (g/cm^3) of aerosol master chemical component type c
!   mw_mastercomp_aer(mc) = molecular weight of aerosol master chemical component type mc
!   name_mastercomp_aer(mc) = name of aerosol master chemical component type mc
!   mc=mastercompptr_aer(c,t)
!   dens_aer(c,t) = dry density (g/cm^3) of aerosol chemical component type c and type t
!   mw_aer(c,t) = molecular weight of aerosol chemical component type c and type t
!   name_aer(c,t) = name of aerosol chemical component type c and type t
!
!   lptr_so4_aer(s,t,p) = gchm r-array index for the
!	mixing ratio for sulfate associated with aerosol size bin s, type t, and phase p
!   (similar for msa, oc, bc, nacl, dust)
!
!-----------------------------------------------------------------------
!
!   volumcen_sect(s,t)= volume (cm^3) at center of section m
!   volumlo_sect(s,t) = volume (cm^3) at lower boundary of section m
!   volumhi_sect(s,t) = volume (cm^3) at upper boundary of section m
!
!   dlo_sect(s,t) = diameter (cm) at lower boundary of section m
!   dhi_sect(s,t) = diameter (cm) at upper boundary of section m
!   dcen_sect(s,t) = volume arithmetic-mean diameter (cm) of section m
!	(corresponds to volumcen_sect == 0.5*(volumlo_sect + volumhi_sect)
!
!-----------------------------------------------------------------------

	integer, parameter :: maxd_atype = 3
	integer, parameter :: maxd_asize = 1
	integer, parameter :: maxd_acomp = 43
	integer, parameter :: maxd_aphase = 2
	integer, save :: ai_phase ! interstitial phase of aerosol
	integer, save :: cw_phase ! cloud water phase of aerosol
	integer, save :: ci_phase ! cloud ice  phase of aerosol
	integer, save :: cr_phase ! rain  phase of aerosol
	integer, save :: cs_phase ! snow  phase of aerosol
	integer, save :: cg_phase ! graupel phase of aerosol

	integer, save :: ntype_aer = 3 ! number of types
!	integer, save :: ntot_mastercomp_aer = 29 ! number of master components
	integer, save :: nphase_aer = 2 ! number of phases (ai, cw only here)

        REAL, PARAMETER :: pirs=3.14159265358979324
        REAL, PARAMETER :: f6dpi=6.0/pirs

	integer, save ::   &
      	  msectional, maerosolincw,   &
      	  nsize_aer( maxd_atype ),   & ! number of size bins
      	  ncomp_aer( maxd_atype ),   & ! number of chemical components
      	  ncomp_aer_nontracer( maxd_atype ),   &
          mastercompptr_aer(maxd_acomp, maxd_atype), &   !  mastercomp index
      	  massptr_aer( maxd_acomp, maxd_asize, maxd_atype, maxd_aphase ), & ! index for mixing ratio
      	  waterptr_aer( maxd_asize, maxd_atype ), & ! index for aerosol water
      	  hygroptr_aer( maxd_asize, maxd_atype ), & ! index for aerosol hygroscopicity
      	  numptr_aer( maxd_asize, maxd_atype, maxd_aphase ), & ! index for the number mixing ratio
          mprognum_aer(maxd_asize,maxd_atype,maxd_aphase)

	real, save ::   &
          dens_aer( maxd_acomp, maxd_atype ),   &
          dens_mastercomp_aer( maxd_acomp ),   &
      	  mw_mastercomp_aer( maxd_acomp ), &
      	  mw_aer( maxd_acomp, maxd_atype ),  &
      	  hygro_mastercomp_aer( maxd_acomp ), &
      	  hygro_aer( maxd_acomp, maxd_atype )
	character*10, save ::   &
      	  name_mastercomp_aer( maxd_acomp ), &
      	  name_aer( maxd_acomp, maxd_atype )

	real, save ::   &
          volumcen_sect( maxd_asize, maxd_atype ),   &
          volumlo_sect( maxd_asize, maxd_atype ),   &
          volumhi_sect( maxd_asize, maxd_atype ),   &
          dcen_sect( maxd_asize, maxd_atype ),   &
          dlo_sect( maxd_asize, maxd_atype ),   &
          dhi_sect( maxd_asize, maxd_atype ),   &
	  sigmag_aer(maxd_asize, maxd_atype)

	integer, save ::                     &
	  lptr_aso4_aer( maxd_asize,maxd_atype,maxd_aphase),     &
	  lptr_ano3_aer( maxd_asize,maxd_atype,maxd_aphase ),    &
	  lptr_anh4_aer( maxd_asize,maxd_atype,maxd_aphase),     &	
	  lptr_orgalk_aer(maxd_asize,maxd_atype,maxd_aphase),    &
	  lptr_orgxyl1_aer(maxd_asize,maxd_atype,maxd_aphase),   &
	  lptr_orgxyl2_aer(maxd_asize,maxd_atype,maxd_aphase),   &
	  lptr_orgxyl3_aer(maxd_asize,maxd_atype,maxd_aphase),   &
	  lptr_orgtol1_aer(maxd_asize,maxd_atype,maxd_aphase),   &
	  lptr_orgtol2_aer(maxd_asize,maxd_atype,maxd_aphase),   &	
	  lptr_orgtol3_aer(maxd_asize,maxd_atype,maxd_aphase),   &	
	  lptr_orgbnz1_aer(maxd_asize,maxd_atype,maxd_aphase),   &	
	  lptr_orgbnz2_aer(maxd_asize,maxd_atype,maxd_aphase),   &
	  lptr_orgbnz3_aer(maxd_asize,maxd_atype,maxd_aphase),   &	
	  lptr_orgtrp1_aer(maxd_asize,maxd_atype,maxd_aphase),   &	
	  lptr_orgtrp2_aer(maxd_asize,maxd_atype,maxd_aphase),   &
	  lptr_orgiso1_aer(maxd_asize,maxd_atype,maxd_aphase),   &	
	  lptr_orgiso2_aer(maxd_asize,maxd_atype,maxd_aphase),   &	
	  lptr_orgsqt_aer(maxd_asize,maxd_atype,maxd_aphase),    &
          lptr_orgISO3_aer(maxd_asize,maxd_atype,maxd_aphase),   & 
	  lptr_orgOLGA_aer(maxd_asize,maxd_atype,maxd_aphase),   & 
	  lptr_orgOLGB_aer(maxd_asize,maxd_atype,maxd_aphase),   & 		  
	  lptr_aorgc_aer(maxd_asize,maxd_atype,maxd_aphase),     &
	  lptr_aorgpa_aer(maxd_asize,maxd_atype,maxd_aphase),    &	  	
	  lptr_aec_aer(maxd_asize,maxd_atype,maxd_aphase),       &	
	  lptr_aothr_aer(maxd_asize,maxd_atype,maxd_aphase),       &		
	  lptr_ana_aer(maxd_asize,maxd_atype,maxd_aphase),       &
	  lptr_acl_aer(maxd_asize,maxd_atype,maxd_aphase),       &
	  lptr_acors_aer(maxd_asize,maxd_atype,maxd_aphase),     &
	  lptr_asoil_aer(maxd_asize,maxd_atype,maxd_aphase) 	

	logical, save ::                     &
      	  do_cloudchem_aer(maxd_asize,maxd_atype)


!   molecular weights (g/mol)
        real, parameter :: mw_aso4_aer    = 96.0
        real, parameter :: mw_ano3_aer    = 62.0
        real, parameter :: mw_anh4_aer    = 18.0
        real, parameter :: mw_AALK_aer    = 150.0	
        real, parameter :: mw_AXYL_aer    = 192.0 	
        real, parameter :: mw_ATOL_aer    = 168.0	
        real, parameter :: mw_ABNZ_aer    = 144.0	
        real, parameter :: mw_ATRP_aer    = 168.0	
        real, parameter :: mw_AISO_aer    = 96.0	
        real, parameter :: mw_ASQT_aer    = 378.0	
        real, parameter :: mw_AISO3_aer   = 162.0	
        real, parameter :: mw_AOLGA_aer   = 176.4	
        real, parameter :: mw_AOLGB_aer   = 252.0	
        real, parameter :: mw_aorgc_aer   = 177.0
        real, parameter :: mw_aorgpa_aer  = 220.0	
        real, parameter :: mw_aec_aer     = 12.0
        real, parameter :: mw_aothr_aer   = 200.0
        real, parameter :: mw_ana_aer     = 23.0
        real, parameter :: mw_acl_aer     = 35.0
        real, parameter :: mw_acors_aer   = 100.0
        real, parameter :: mw_asoil_aer   = 100.0
	real, parameter :: mw_water_aer   = 18.016
	  
!   dry densities (g/cm3)

        real, parameter :: dens_aso4_aer    = 1.8       ! = rhoso4
        real, parameter :: dens_ano3_aer    = 1.8       ! = rhono3
        real, parameter :: dens_anh4_aer    = 1.8       ! = rhonh4
        real, parameter :: dens_AALK_aer    = 2.0	! cmaq	
        real, parameter :: dens_AXYL_aer    = 2.0	! cmaq
        real, parameter :: dens_ATOL_aer    = 2.0	! cmaq
        real, parameter :: dens_ABNZ_aer    = 2.0	! cmaq
        real, parameter :: dens_ATRP_aer    = 2.0	! cmaq
        real, parameter :: dens_AISO_aer    = 2.0	! cmaq
        real, parameter :: dens_ASQT_aer    = 2.0	! cmaq
        real, parameter :: dens_AISO3_aer   = 2.0	! cmaq
        real, parameter :: dens_AOLGA_aer   = 2.0	! cmaq
        real, parameter :: dens_AOLGB_aer   = 2.0	! cmaq
        real, parameter :: dens_aorgc_aer   = 2.0       ! cmaq
        real, parameter :: dens_aorgpa_aer  = 2.0       ! cmaq	
        real, parameter :: dens_aec_aer     = 2.2       ! cmaq
        real, parameter :: dens_aothr_aer   = 2.2       ! cmaq
        real, parameter :: dens_ana_aer     = 2.2       ! cmaq
        real, parameter :: dens_acl_aer     = 2.2       ! cmaq
        real, parameter :: dens_acors_aer   = 2.2       ! cmaq
        real, parameter :: dens_asoil_aer   = 2.6       ! cmaq
!   water density (g/cm3)
	real, parameter :: dens_water_aer   = 1.0

!   hygroscopicity (dimensionless)

        real, parameter :: hygro_aso4_aer    = 0.50
        real, parameter :: hygro_ano3_aer    = 0.50
        real, parameter :: hygro_anh4_aer    = 0.50
        real, parameter :: hygro_AALK_aer    = 0.14		
        real, parameter :: hygro_AXYL_aer    = 0.14	
        real, parameter :: hygro_ATOL_aer    = 0.14	
        real, parameter :: hygro_ABNZ_aer    = 0.14	
        real, parameter :: hygro_ATRP_aer    = 0.14	
        real, parameter :: hygro_AISO_aer    = 0.14	
        real, parameter :: hygro_ASQT_aer    = 0.14	
        real, parameter :: hygro_AISO3_aer   = 0.14	
        real, parameter :: hygro_AOLGA_aer   = 0.20	
        real, parameter :: hygro_AOLGB_aer   = 0.20	 
        real, parameter :: hygro_aorgc_aer   = 0.50      
        real, parameter :: hygro_aorgpa_aer  = 0.14		
        real, parameter :: hygro_aec_aer     = 1.e-6    
        real, parameter :: hygro_aothr_aer   = 0.10      
        real, parameter :: hygro_ana_aer     = 1.16     
        real, parameter :: hygro_acl_aer     = 1.16     
        real, parameter :: hygro_acors_aer   = 0.03      
        real, parameter :: hygro_asoil_aer   = 0.03      


END Module module_data_cmaq
