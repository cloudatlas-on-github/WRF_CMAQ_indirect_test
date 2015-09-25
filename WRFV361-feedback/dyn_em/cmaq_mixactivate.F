!----------------------------------------------------------------------
! This module contains interface wrapper routines to couple the cmaq aerosol
!  with mixactivate in the physics directory (Shaocai Yu, 7/10/2011). 
!----------------------------------------------------------------------

      subroutine cmaq_mixactivate (                                    &
           id, ktau, dtstep, config_flags, idrydep_onoff,              &
           rho_phy, p_phy,t_phy, w, cldfra, cldfra_old,                &
           ddvel, z, dz8w, p_at_w, t_at_w, exch_h,                     &
           qv, qc, qi,qs,f_ice_phy,f_rain_phy,qndrop3d,                &
           qnice3d,f_qndrop,f_qnice,                                   &
           rndst1,rndst2,rndst3,rndst4,nacon1,nacon2,nacon3,nacon4,    &
	   f_qv,f_qc, f_qi,f_qs,                                       &
	   ccn1, ccn2, ccn3, ccn4, ccn5, ccn6, nsource,                &
           dgn_1,dgn_2,dgn_3,sig_1,sig_2,sig_3,                        &
           ae_mass_01, ae_mass_02, ae_mass_03, ae_mass_04, ae_mass_05, &
           ae_mass_06, ae_mass_07, ae_mass_08, ae_mass_09, ae_mass_10, &
           ae_mass_11, ae_mass_12, ae_mass_13, ae_mass_14, ae_mass_15, &
           ae_mass_16, ae_mass_17, ae_mass_18, ae_mass_19, ae_mass_20, &
           ae_mass_21, ae_mass_22, ae_mass_23, ae_mass_24, ae_mass_25, &
           ae_mass_26, ae_mass_27, ae_mass_28, ae_mass_29, ae_mass_30, &
           ae_mass_31, ae_mass_32, ae_mass_33, ae_mass_34, ae_mass_35, &
           ae_mass_36, ae_mass_37, ae_mass_38, ae_mass_39, ae_mass_40, &
           ae_mass_41, ae_mass_42, ae_mass_43,                         &
           ae_num_1, ae_num_2,ae_num_3,                                &
           ids,ide, jds,jde, kds,kde,                                  &
           ims,ime, jms,jme, kms,kme,                                  &
           its,ite, jts,jte, kts,kte,                                  &
!----------------------------------------------------------------------
!Added by R Atlas for thompson indirect effect
           xland, dx, dy, qnum_wf, f_qnum_wf, qnum_if, f_qnum_if,       &
           hygro_1, hygro_2, hygro_3,                                  &
           morrison_indirect, thompson_indirect )
!----------------------------------------------------------------------
	   	
        use module_configure, only: grid_config_rec_type
	use module_state_description, only:  num_chem
	use module_data_cmaq
        use module_aerosols_cmaq_init_aercld
!	use module_Cmixactivate, only: Cmixactivate
	use module_Cmixactivate, only: Cmixactivate,cmaq_icenucleat

!-----------------------------------------------------
!    ..../frame/module_configure.F
!    ..../frame/module_state_description.F 
!    ..../dyn_em/module_data_cmaq.F
!    ..../dyn_em/module_aerosols_cmaq_init_aercld.F
!    ..../phys/module_Cmixactivate.F
!----------------------------------------------------	

!--------------------------------------	
! wrapper to call mixactivate for cmaq description of aerosol
!  REVISION HISTORY:
!
!   10 March 2010 Shaocai Yu:  create the original code
!   7 October 2011 Shaocai Yu: add ice parameterization code
!-----------------------------------------------------------

implicit none

!-------------------------------------------------------------------	
! Definitions for each variable
!-------------------------------------------------------------------------
!--id                                     !WRF grid%id
!--ktau                                   !WRF time step count
!--dtstep                                 !time step for microphysics (s)
!--config_flags                           !configuration and control parameters
!--idrydep_onoff    
!--rho_phy(ims:ime, kms:kme, jms:jme)     ! density of air at mid-level  (kg/m3)
!--t_phy	                          ! potential temperature    (K)
!--w                                      ! vertical velocity (m/s)
!--cldfra                                 !current cloud fraction
!--cldfra_old                             !cloud fraction on previous time step        
!--ddvel( its:ite, jts:jte, num_chem )    ! dry deposition velocity  (m/s)
!--z                                      !3D height with lowest level(Height above sea level, m)  
!--dz8w(ims:ime, kms:kme, jms:jme)        ! layer thickness between full levels (m) 
!--p_at_w(ims:ime, kms:kme, jms:jme)      ! 3D pressure at full levels (Pa)
!--t_at_w(ims:ime, kms:kme, jms:jme)      ! temperature at layer interfaces 
!--exch_h                                 !vertical diffusivity (m2/s)           
!--qv                                     !water vapor mixing ratio (kg/kg)
!--qc                                     !cloud water mixing ratio (kg/kg)
!--qi                                     !cloud ice mixing ratio (kg/kg)
!--qndrop3d                               !droplet number mixing ratio (#/kg)
!--f_qc                                   !LOGICAL, OPTIONAL :: f_qc
!--f_qi      L                            !OGICAL, OPTIONAL :: f_qi          
!--ccn1,ccn2,ccn3,ccn4,ccn5,ccn6          ! number conc of aerosols activated at supersat
!                                         !(/0.02,0.05,0.1,0.2,0.5,1.0/)
!--nsource(ims:ime,kms:kme,jms:jme)       ! droplet number source (#/kg/s)   
!--dgn(ims:ime, kms:kme, jms:jme, 3)      !CMAQ diameter for each mode (cm)
!--sig(ims:ime, kms:kme, jms:jme, 3)      !CMAQ standard deviation for each mode
!--ae_mass(ims:ime, kms:kme, jms:jme, 28) !CMAQ mass concentration for each species (mass mixing ratio, kg/kg-air)
!--ae_num (ims:ime, kms:kme, jms:jme, 3)  !CMAQ number concentration for each mode(#/kg-air)                                                                    
!-- ids                                   !start index for i in domain
!-- ide                                   !end index for i in domain
!-- jds                                   !start index for j in domain
!-- jde                                   !end index for j in domain
!-- kds                                   !start index for k in domain
!-- kde                                   !end index for k in domain 
!-- ims                                   !start index for i in memory
!-- ime                                   !end index for i in memory
!-- jms                                   !start index for j in memory
!-- jme                                   !end index for j in memory
!-- kms                                   !start index for k in memory
!-- kme                                   !end index for k in memory
!-- its                                   !start index for i in tile
!-- ite                                   !end index for i in tile
!-- jts                                   !start index for j in tile
!-- jte                                   !end index for j in tile
!-- kts                                   !start index for k in tile
!-- kte                                   !end index for k in tile                      
!----------------------------------------------------------------------------------------------

! David: use this temporarily
	integer, parameter :: loc_num_chem = 93 !93  !P_H2Ocorn    !500

!   subr arguments

!	integer, intent(in) ::                  &
!		ids, ide, jds, jde, kds, kde,   &
!		ims, ime, jms, jme, kms, kme,   &
!		its, ite, jts, jte, kts, kte

	integer, intent(in) ::                  &
		id, ktau,                       &
		ids, ide, jds, jde, kds, kde,   &
		ims, ime, jms, jme, kms, kme,   &
		its, ite, jts, jte, kts, kte,   &
                idrydep_onoff
	real, intent(in) :: dtstep, dx, dy
	real, intent(in),   &
		dimension( ims:ime, kms:kme, jms:jme ) :: &
		rho_phy, p_phy,t_phy, w,   &
		z, dz8w, p_at_w, t_at_w, exch_h
	real, intent(inout),   &
		dimension( ims:ime, kms:kme, jms:jme ) :: cldfra, cldfra_old
	real,    &
		dimension( its:ite, jts:jte, loc_num_chem ) :: ddvel
	real, intent(inout),   &
		dimension( ims:ime, kms:kme, jms:jme ) :: &
		qv, qc, qi, qs
        LOGICAL, intent(in) :: f_qc, f_qi
	real, intent(inout),   &
		dimension( ims:ime, kms:kme, jms:jme  ) :: &
		qndrop3d
	real, save, allocatable :: chem (:,:,:,:)
        real, intent(out), dimension(ims:ime,kms:kme,jms:jme) :: nsource, &
	     ccn1,ccn2,ccn3,ccn4,ccn5,ccn6  ! number conc of aerosols activated at supersat
	type(grid_config_rec_type), intent(in) :: config_flags
!--------------------------------------------------------For test purpose
	real,    &
		dimension( ims:ime, kms:kme, jms:jme  ) :: &
		Tqnice3d
	real,    &
		dimension( ims:ime, kms:kme, jms:jme ) :: Tcldfra, Tcldfra_old
	real, &
		dimension( ims:ime, kms:kme, jms:jme ) :: Tf_ice_phy,Tf_rain_phy

!--------------------------------------------------------------------
! variables for CMAQ aerosol information

        REAL, DIMENSION (ims:ime, kms:kme, jms:jme ),         &
              INTENT(INout)  ::  dgn_1,dgn_2,dgn_3,sig_1,sig_2,sig_3,ae_num_1,ae_num_2,ae_num_3
        REAL, DIMENSION (ims:ime, kms:kme, jms:jme ),         &
              INTENT(INout)  ::  &
           ae_mass_01, ae_mass_02, ae_mass_03, ae_mass_04, ae_mass_05, &
           ae_mass_06, ae_mass_07, ae_mass_08, ae_mass_09, ae_mass_10, &
           ae_mass_11, ae_mass_12, ae_mass_13, ae_mass_14, ae_mass_15, &
           ae_mass_16, ae_mass_17, ae_mass_18, ae_mass_19, ae_mass_20, &
           ae_mass_21, ae_mass_22, ae_mass_23, ae_mass_24, ae_mass_25, &
           ae_mass_26, ae_mass_27, ae_mass_28, ae_mass_29, ae_mass_30, &
           ae_mass_31, ae_mass_32, ae_mass_33, ae_mass_34, ae_mass_35, &
           ae_mass_36, ae_mass_37, ae_mass_38, ae_mass_39, ae_mass_40, &
           ae_mass_41, ae_mass_42, ae_mass_43

        REAL, DIMENSION (ims:ime, kms:kme, jms:jme, 43 ) :: ae_mass

!	real, PARAMETER :: pirs=3.14159265358979324
        INTEGER, PARAMETER :: max_chem =100	     
        REAL, DIMENSION(ims:ime,jms:jme), INTENT(IN):: xland

!-----------------------------------------------------------------	
! dgn(its:ite,kts:kte,jts:jte,3): mean diameter for nuclei (1),accumulation (1),coarse(1) mode (cm)
! sig(its:ite,kts:kte,jts:jte,3): standard deivation for nuclei (1),accumulation (1),coarse(3) mode
! ae_num(its:ite,kts:kte,jts:jte,3) : #/kg-air;
! ae_mass(its:ite,kts:kte,jts:jte,28):ug/kg-air (mixing ratio)
!------------------------------------------------------------------     
       REAL, DIMENSION (ims:ime, kms:kme, jms:jme, maxd_asize, maxd_atype):: &
          cvolumcen_sect,cvolumlo_sect,cvolumhi_sect,cdcen_sect,   &
          cdlo_sect, cdhi_sect,csigmag_aer
       
       real hygro( its:ite, kts:kte, jts:jte,maxd_asize, maxd_atype )	  	  
       real sumhygro,sumvol	     
       REAL, DIMENSION (ims:ime, kms:kme, jms:jme):: cdp_meanvol_tmp
!log(sginin),log(sginia),log(sginic)   
!----------------------------------------------------------------                        
       REAL, DIMENSION (ims:ime, kms:kme, jms:jme, 3 )::       &
             xxlsg    
!-----------------------------------------------------------------                    
! log(sginin ) ** 2,log(sginina) ** 2, log(sginina) ** 2     
!--------------------------------------------------------------------                     
       REAL, DIMENSION (ims:ime, kms:kme, jms:jme, 3 )::      &
             l2sig     
	     
!-------------------------------------------------------------------
!   initialize pointers used by aerosol-cloud-interaction routines
!     Currently, for all local variables (12/03/2009)
!       like SUBROUTINE aerosols_sorgam_init_aercld_ptrs 
!       in module_aerosol_sorgam.F
!-------------------------------------------------------------------    
    
!    INTEGER loc_num_chem
    LOGICAL is_aerosol(loc_num_chem)
    
! local vars    
	real qsrflx(ims:ime, jms:jme, loc_num_chem) ! dry deposition flux of aerosol

!     integer iphase, isize, itype, l, ll, n, p1st
     integer iphase, isize, itype, ll, n, p1st,p,m,l,j,k,i
     REAL dp_meanvol_tmp
     
     INTEGER , PARAMETER :: PARAM_FIRST_SCALAR = 1
     integer ptr, error_stat

     real, allocatable, save :: cfrac_num (:, :, :, :)
     real, allocatable, save :: cfrac_mass (:, :, :, :)
     integer :: stat
     logical, save :: firstime = .true.

!--------------for ice parameterization 
      real, intent(inout),   &
		dimension( ims:ime, kms:kme, jms:jme ) :: f_ice_phy,f_rain_phy
      LOGICAL, intent(in) :: f_qndrop,f_qnice, f_qv,f_qs
      real, intent(inout), dimension( ims:ime, kms:kme, jms:jme  ) :: &
		qnice3d     
      real, dimension(ims:ime,kms:kme,jms:jme) ::  &
        nihf3d,niimm3d,nidep3d,nimey3d
      real, intent(inout),dimension(ims:ime,kms:kme,jms:jme) ::  &
        rndst1,rndst2,rndst3,rndst4,nacon1,nacon2,nacon3,nacon4  ! dust for contact freezing 
!----------------------------------------------------------------------------
!Added by R Atlas for thompson indirect effect
       real, intent(out), dimension(its:ite, kts:kte, jts:jte) :: &
                hygro_1, hygro_2, hygro_3
       real, intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: &
                qnum_wf, qnum_if
       logical, intent(in) :: f_qnum_wf, f_qnum_if     
       logical, intent(in) :: morrison_indirect, thompson_indirect     
       real :: ii
!-----------------------------------------------------------
 
!   param_first_scalar is from module_state_description (???)   
!        nphase_aer = 1
!     	if(p_so4cwj.ge. param_first_scalar) then
!                nphase_aer = 2
!     	endif
!-------Shaocai------------------------------------
!  initial mean diameter for nuclei mode [ m ]    
!      REAL dginin
!      PARAMETER (dginin=0.01E-6)
!
!  initial mean diameter for accumulation mode [ m ]
!      REAL dginia
!      PARAMETER (dginia=0.07E-6)
!
! initial mean diameter for coarse mode [ m ]  
!      REAL dginic
!      PARAMETER (dginic=1.0E-6)
!  initial sigma-G for nucleimode                 
!      REAL sginin
!      PARAMETER (sginin=1.70)
!
!  initial sigma-G for accumulation mode          
!      REAL sginia
!      PARAMETER (sginia=2.00)
!
! initial sigma-G for coarse mode               
!     REAL sginic
!      PARAMETER (sginic=2.5)
!    ae_mass (ug/kg-air) ae_num (#/kg-air)
!----------------------------------------------------          
!              dgn(:,:,:,1)= 0.01E-6
!              dgn(:,:,:,2)= 0.07E-6
!              dgn(:,:,:,3)= 1.0E-6	      
!	      sig(:,:,:,1)= 1.02
!	      sig(:,:,:,2)= 1.02	      
!	      sig(:,:,:,3)= 2.20	      
!	      ae_num(:,:,:,:)=100.e6  !1000.e6
!              ae_mass(:,:,:,:)=10.0

!--------------------need initialize qndrop3d(:,:,:) to be zero because we create cloud drop
!from total CMAQ-aerosols instead of just interstital aerosol at this stage (Shaocai 3/15/2010)
!---------------------------------------------------------

!              dgn_2(:,:,:)= dgn_2(:,:,:)*2.0
!              dgn_3(:,:,:)= dgn_3(:,:,:)*2.0
	      
	      
           ae_mass (:,:,:, 1) = ae_mass_01
           ae_mass (:,:,:, 2) = ae_mass_02
           ae_mass (:,:,:, 3) = ae_mass_03
           ae_mass (:,:,:, 4) = ae_mass_04
           ae_mass (:,:,:, 5) = ae_mass_05
           ae_mass (:,:,:, 6) = ae_mass_06
           ae_mass (:,:,:, 7) = ae_mass_07
           ae_mass (:,:,:, 8) = ae_mass_08
           ae_mass (:,:,:, 9) = ae_mass_09
           ae_mass (:,:,:,10) = ae_mass_10
           ae_mass (:,:,:,11) = ae_mass_11
           ae_mass (:,:,:,12) = ae_mass_12
           ae_mass (:,:,:,13) = ae_mass_13
           ae_mass (:,:,:,14) = ae_mass_14
           ae_mass (:,:,:,15) = ae_mass_15
           ae_mass (:,:,:,16) = ae_mass_16
           ae_mass (:,:,:,17) = ae_mass_17
           ae_mass (:,:,:,18) = ae_mass_18
           ae_mass (:,:,:,19) = ae_mass_19
           ae_mass (:,:,:,20) = ae_mass_20
           ae_mass (:,:,:,21) = ae_mass_21
           ae_mass (:,:,:,22) = ae_mass_22
           ae_mass (:,:,:,23) = ae_mass_23
           ae_mass (:,:,:,24) = ae_mass_24
           ae_mass (:,:,:,25) = ae_mass_25
           ae_mass (:,:,:,26) = ae_mass_26
           ae_mass (:,:,:,27) = ae_mass_27
           ae_mass (:,:,:,28) = ae_mass_28
           ae_mass (:,:,:,29) = ae_mass_29
           ae_mass (:,:,:,30) = ae_mass_30
           ae_mass (:,:,:,31) = ae_mass_31
           ae_mass (:,:,:,32) = ae_mass_32
           ae_mass (:,:,:,33) = ae_mass_33
           ae_mass (:,:,:,34) = ae_mass_34
           ae_mass (:,:,:,35) = ae_mass_35
           ae_mass (:,:,:,36) = ae_mass_36
           ae_mass (:,:,:,37) = ae_mass_37
           ae_mass (:,:,:,38) = ae_mass_38
           ae_mass (:,:,:,39) = ae_mass_39
           ae_mass (:,:,:,40) = ae_mass_40
           ae_mass (:,:,:,41) = ae_mass_41
           ae_mass (:,:,:,42) = ae_mass_42
           ae_mass (:,:,:,43) = ae_mass_43

           if (firstime) then
              firstime = .false.
              allocate (cfrac_num (ims:ime, kms:kme, jms:jme, 200), stat=stat)
              allocate (cfrac_mass (ims:ime, kms:kme, jms:jme, 200), stat=stat)
              cfrac_num = 1.0
              cfrac_mass = 1.0
           end if

!               qndrop3d(:,:,:)=0.0

	   	  
!             goto 999

        if (.not. allocated(chem)) then
           allocate (chem( ims:ime, kms:kme, jms:jme, &
                     loc_num_chem), stat=error_stat)
!                    (nphase_aer + maxd_aphase * sum(ncomp_aer)) * ntype_aer * sum(nsize_aer) ) , stat=error_stat)
        end if

	ai_phase=-999888777
	cw_phase=-999888777
	ci_phase=-999888777
	cr_phase=-999888777
	cs_phase=-999888777
	cg_phase=-999888777
	if(nphase_aer>=1)ai_phase=1
	if(nphase_aer>=2)cw_phase=2
	if(nphase_aer>=3)cr_phase=3
	if(nphase_aer>=4)ci_phase=4
	if(nphase_aer>=5)cw_phase=5
	if(nphase_aer>=6)cg_phase=6

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
!--------------------------------------------------------
	     nsize_aer(1)=1  ! nsize_aer(ntype_aer)
	     nsize_aer(2)=1
             nsize_aer(3)=1

	msectional = 0
	maerosolincw = 0
!--------------------------------------------------------
! maerosolincw - if positive, both unactivated/interstitial and activated
!	aerosol species are simulated.  if zero/negative, only the
!	unactivated are simulated. (???)	
!--------------------------------------------------	
#if defined ( cw_species_are_in_registry )
	maerosolincw = 1
#endif

        write(6,*)'maerosolincw,msectional=',maerosolincw,msectional
	
	name_mastercomp_aer( 1) = 'sulfate'
	dens_mastercomp_aer( 1) = dens_aso4_aer
	mw_mastercomp_aer(   1) =   mw_aso4_aer
	hygro_mastercomp_aer(1) = hygro_aso4_aer

	name_mastercomp_aer( 2) = 'nitrate'
	dens_mastercomp_aer( 2) = dens_ano3_aer
	mw_mastercomp_aer(   2) =   mw_ano3_aer
	hygro_mastercomp_aer(2) = hygro_ano3_aer

	name_mastercomp_aer( 3) = 'ammonium'
	dens_mastercomp_aer( 3) = dens_anh4_aer
	mw_mastercomp_aer(   3) =   mw_anh4_aer
	hygro_mastercomp_aer(3) = hygro_anh4_aer

	name_mastercomp_aer( 4) = 'orgalk'
	dens_mastercomp_aer( 4) = dens_AALK_aer
	mw_mastercomp_aer(   4) =   mw_AALK_aer
	hygro_mastercomp_aer(4) = hygro_AALK_aer
	
	name_mastercomp_aer( 5) = 'orgxyl1'
	dens_mastercomp_aer( 5) = dens_AXYL_aer
	mw_mastercomp_aer(   5) =   mw_AXYL_aer
	hygro_mastercomp_aer(5) = hygro_AXYL_aer
	
	name_mastercomp_aer( 6) = 'orgxyl2'
	dens_mastercomp_aer( 6) = dens_AXYL_aer
	mw_mastercomp_aer(   6) =   mw_AXYL_aer
	hygro_mastercomp_aer(6) = hygro_AXYL_aer
	
	name_mastercomp_aer( 7) = 'orgxyl3'
	dens_mastercomp_aer( 7) = dens_AXYL_aer
	mw_mastercomp_aer(   7) =   mw_AXYL_aer
	hygro_mastercomp_aer(7) = hygro_AXYL_aer

	name_mastercomp_aer( 8) = 'orgtol1'
	dens_mastercomp_aer( 8) = dens_ATOL_aer
	mw_mastercomp_aer(   8) =   mw_ATOL_aer
	hygro_mastercomp_aer(8) = hygro_ATOL_aer

	name_mastercomp_aer( 9) = 'orgtol2'
	dens_mastercomp_aer( 9) = dens_ATOL_aer
	mw_mastercomp_aer(   9) =   mw_ATOL_aer
	hygro_mastercomp_aer(9) = hygro_ATOL_aer
	
	name_mastercomp_aer( 10) = 'orgtol3'
	dens_mastercomp_aer( 10) = dens_ATOL_aer
	mw_mastercomp_aer(   10) =   mw_ATOL_aer
	hygro_mastercomp_aer(10) = hygro_ATOL_aer

	name_mastercomp_aer( 11) = 'orgbnz1'
	dens_mastercomp_aer( 11) = dens_ABNZ_aer
	mw_mastercomp_aer(   11) =   mw_ABNZ_aer
	hygro_mastercomp_aer(11) = hygro_ABNZ_aer
	
	name_mastercomp_aer( 12) = 'orgbnz2'
	dens_mastercomp_aer( 12) = dens_ABNZ_aer
	mw_mastercomp_aer(   12) =   mw_ABNZ_aer
	hygro_mastercomp_aer(12) = hygro_ABNZ_aer
			
	name_mastercomp_aer( 13) = 'orgbnz3'
	dens_mastercomp_aer( 13) = dens_ABNZ_aer
	mw_mastercomp_aer(   13) =   mw_ABNZ_aer
	hygro_mastercomp_aer(13) = hygro_ABNZ_aer
	
	name_mastercomp_aer( 14) = 'orgtrp1'
	dens_mastercomp_aer( 14) = dens_ATRP_aer
	mw_mastercomp_aer(   14) =   mw_ATRP_aer
	hygro_mastercomp_aer(14) = hygro_ATRP_aer
	
	name_mastercomp_aer( 15) = 'orgtrp2'
	dens_mastercomp_aer( 15) = dens_ATRP_aer
	mw_mastercomp_aer(   15) =   mw_ATRP_aer
	hygro_mastercomp_aer(15) = hygro_ATRP_aer		

	name_mastercomp_aer( 16) = 'orgiso1'
	dens_mastercomp_aer( 16) = dens_AISO_aer
	mw_mastercomp_aer(   16) =   mw_AISO_aer
	hygro_mastercomp_aer(16) = hygro_AISO_aer

	name_mastercomp_aer( 17) = 'orgiso2'
	dens_mastercomp_aer( 17) = dens_AISO_aer
	mw_mastercomp_aer(   17) =   mw_AISO_aer
	hygro_mastercomp_aer(17) = hygro_AISO_aer
	
	name_mastercomp_aer( 18) = 'orgsqt'
	dens_mastercomp_aer( 18) = dens_ASQT_aer
	mw_mastercomp_aer(   18) =   mw_ASQT_aer
	hygro_mastercomp_aer(18) = hygro_ASQT_aer	

	name_mastercomp_aer( 19) = 'orgISO3'
	dens_mastercomp_aer( 19) = dens_AISO3_aer
	mw_mastercomp_aer(   19) =   mw_AISO3_aer
	hygro_mastercomp_aer(19) = hygro_AISO3_aer

	name_mastercomp_aer( 20) = 'orgOLGA'
	dens_mastercomp_aer( 20) = dens_AOLGA_aer
	mw_mastercomp_aer(   20) =   mw_AOLGA_aer
	hygro_mastercomp_aer(20) = hygro_AOLGA_aer

	name_mastercomp_aer( 21) = 'orgOLGB'
	dens_mastercomp_aer( 21) = dens_AOLGB_aer
	mw_mastercomp_aer(   21) =   mw_AOLGB_aer
	hygro_mastercomp_aer(21) = hygro_AOLGB_aer
		
	name_mastercomp_aer( 22) = 'cloud_SOA'
	dens_mastercomp_aer( 22) = dens_aorgc_aer
	mw_mastercomp_aer(   22) =   mw_aorgc_aer
	hygro_mastercomp_aer(22) = hygro_aorgc_aer

	name_mastercomp_aer( 23) = 'primary_org'
	dens_mastercomp_aer( 23) = dens_aorgpa_aer
	mw_mastercomp_aer(   23) =   mw_aorgpa_aer
	hygro_mastercomp_aer(23) = hygro_aorgpa_aer
	
	name_mastercomp_aer( 24) = 'aec'
	dens_mastercomp_aer( 24) = dens_aec_aer
	mw_mastercomp_aer(   24) =   mw_aec_aer
	hygro_mastercomp_aer(24) = hygro_aec_aer
	
	name_mastercomp_aer( 25) = 'unspeci_aothr'
	dens_mastercomp_aer( 25) = dens_aothr_aer
	mw_mastercomp_aer(   25) =   mw_aothr_aer
	hygro_mastercomp_aer(25) = hygro_aothr_aer	

  	name_mastercomp_aer( 26) = 'sodium'
  	dens_mastercomp_aer( 26) = dens_ana_aer
  	mw_mastercomp_aer(   26) =   mw_ana_aer
  	hygro_mastercomp_aer(26) = hygro_ana_aer
  
  	name_mastercomp_aer( 27) = 'chloride'
  	dens_mastercomp_aer( 27) = dens_acl_aer
  	mw_mastercomp_aer(   27) =   mw_acl_aer
  	hygro_mastercomp_aer(27) = hygro_acl_aer
	
	name_mastercomp_aer( 28) = 'coarse'
	dens_mastercomp_aer( 28) = dens_acors_aer
	mw_mastercomp_aer(   28) =   mw_acors_aer
	hygro_mastercomp_aer(28) = hygro_acors_aer

	name_mastercomp_aer( 29) = 'asoil'
	dens_mastercomp_aer( 29) = dens_asoil_aer
	mw_mastercomp_aer(   29) =  mw_asoil_aer
	hygro_mastercomp_aer(29) = hygro_asoil_aer

	lptr_aso4_aer(    :,:,:) = 0
	lptr_ano3_aer(    :,:,:) = 0
	lptr_anh4_aer(    :,:,:) = 0	
	lptr_orgalk_aer(  :,:,:) = 0
	lptr_orgxyl1_aer( :,:,:) = 0
	lptr_orgxyl2_aer( :,:,:) = 0
	lptr_orgxyl3_aer( :,:,:) = 0
	lptr_orgtol1_aer( :,:,:) = 0
	lptr_orgtol2_aer( :,:,:) = 0	
	lptr_orgtol3_aer( :,:,:) = 0	
	lptr_orgbnz1_aer( :,:,:) = 0	
	lptr_orgbnz2_aer( :,:,:) = 0
	lptr_orgbnz3_aer( :,:,:) = 0	
	lptr_orgtrp1_aer( :,:,:) = 0	
	lptr_orgtrp2_aer( :,:,:) = 0
	lptr_orgiso1_aer( :,:,:) = 0	
	lptr_orgiso2_aer( :,:,:) = 0	
	lptr_orgsqt_aer(  :,:,:) = 0
	lptr_orgISO3_aer( :,:,:) = 0	
	lptr_orgOLGA_aer( :,:,:) = 0	
	lptr_orgOLGB_aer( :,:,:) = 0				
	lptr_aorgc_aer(   :,:,:) = 0
	lptr_aorgpa_aer(  :,:,:) = 0		
	lptr_aec_aer(     :,:,:) = 0	
	lptr_aothr_aer(   :,:,:) = 0		
	lptr_ana_aer(     :,:,:) = 0
	lptr_acl_aer(     :,:,:) = 0
	lptr_acors_aer(   :,:,:) = 0
	lptr_asoil_aer(   :,:,:) = 0	
	numptr_aer(       :,:,:) = 0

	do_cloudchem_aer(:,:) = .false.

! Aitken mode
	itype = 1
	isize = 1
	ncomp_aer(itype) = 8
!-----------------------------------------------------------
!        p_aso4ai etc. should be defined 
!            in module_state_description.F  (???)
!-----------------------------------------------------------	
	numptr_aer(     isize,itype,ai_phase) = p_nu0
	lptr_aso4_aer(  isize,itype,ai_phase) = p_aso4ai
	lptr_ano3_aer(  isize,itype,ai_phase) = p_ano3ai
	lptr_anh4_aer(  isize,itype,ai_phase) = p_anh4ai
	lptr_aorgpa_aer(isize,itype,ai_phase) = p_aorgpaai
	lptr_aec_aer(   isize,itype,ai_phase) = p_aecai
	lptr_aothr_aer( isize,itype,ai_phase) = p_aothrai		
	lptr_ana_aer(   isize,itype,ai_phase) = p_anaai
	lptr_acl_aer(   isize,itype,ai_phase) = p_aclai
	
! aerosol in cloud water
        if(cw_phase.gt.0)then	
	  numptr_aer(     isize,itype,cw_phase) = p_nu0cw
	  lptr_aso4_aer(  isize,itype,cw_phase) = p_aso4cwi
	  lptr_ano3_aer(  isize,itype,cw_phase) = p_ano3cwi
	  lptr_anh4_aer(  isize,itype,cw_phase) = p_anh4cwi
	  lptr_aorgpa_aer(isize,itype,cw_phase) = p_aorgpacwi	
	  lptr_aec_aer(   isize,itype,cw_phase) = p_aeccwi
	  lptr_aothr_aer( isize,itype,cw_phase) = p_aothrcwi		
	  lptr_ana_aer(   isize,itype,cw_phase) = p_anacwi
	  lptr_acl_aer(   isize,itype,cw_phase) = p_aclcwi
	  do_cloudchem_aer(isize,itype) = .true.
	endif

! Accumulation mode
	itype = 2
	isize = 1
	ncomp_aer(itype) = 28

	numptr_aer(      isize,itype,ai_phase) = p_ac0
	lptr_aso4_aer(   isize,itype,ai_phase) = p_aso4aj
	lptr_ano3_aer(   isize,itype,ai_phase) = p_ano3aj
	lptr_anh4_aer(   isize,itype,ai_phase) = p_anh4aj
	lptr_orgalk_aer( isize,itype,ai_phase) = p_aalkaj
	lptr_orgxyl1_aer(isize,itype,ai_phase) = p_axyl1aj
	lptr_orgxyl2_aer(isize,itype,ai_phase) = p_axyl2aj
	lptr_orgxyl3_aer(isize,itype,ai_phase) = p_axyl3aj
	lptr_orgtol1_aer(isize,itype,ai_phase) = p_atol1aj
	lptr_orgtol2_aer(isize,itype,ai_phase) = p_atol2aj	
	lptr_orgtol3_aer(isize,itype,ai_phase) = p_atol3aj	
	lptr_orgbnz1_aer(isize,itype,ai_phase) = p_abnz1aj	
	lptr_orgbnz2_aer(isize,itype,ai_phase) = p_abnz2aj
	lptr_orgbnz3_aer(isize,itype,ai_phase) = p_abnz3aj	
	lptr_orgtrp1_aer(isize,itype,ai_phase) = p_atrp1aj	
	lptr_orgtrp2_aer(isize,itype,ai_phase) = p_atrp2aj
	lptr_orgiso1_aer(isize,itype,ai_phase) = p_aiso1aj	
	lptr_orgiso2_aer(isize,itype,ai_phase) = p_aiso2aj	
	lptr_orgsqt_aer( isize,itype,ai_phase) = p_asqtaj	
	lptr_orgISO3_aer(isize,itype,ai_phase) = p_aISO3aj	
	lptr_orgOLGA_aer(isize,itype,ai_phase) = p_aOLGAaj
	lptr_orgOLGB_aer(isize,itype,ai_phase) = p_aOLGBaj			
	lptr_aorgc_aer(  isize,itype,ai_phase) = p_aorgcaj
	lptr_aorgpa_aer( isize,itype,ai_phase) = p_aorgpaaj		
	lptr_aec_aer(    isize,itype,ai_phase) = p_aecaj
	lptr_aothr_aer(  isize,itype,ai_phase) = p_aothraj		
	lptr_ana_aer(    isize,itype,ai_phase) = p_anaaj
	lptr_acl_aer(    isize,itype,ai_phase) = p_aclaj
	lptr_asoil_aer(  isize,itype,ai_phase) = p_asoilaj
		
! aerosol in cloud water
        if(cw_phase.gt.0)then	
	  numptr_aer(      isize,itype,cw_phase) = p_ac0cw
	  lptr_aso4_aer(   isize,itype,cw_phase) = p_aso4cwj
	  lptr_ano3_aer(   isize,itype,cw_phase) = p_ano3cwj
	  lptr_anh4_aer(   isize,itype,cw_phase) = p_anh4cwj
  	  lptr_orgalk_aer( isize,itype,cw_phase) = p_aalkcwj
	  lptr_orgxyl1_aer(isize,itype,cw_phase) = p_axyl1cwj
	  lptr_orgxyl2_aer(isize,itype,cw_phase) = p_axyl2cwj
	  lptr_orgxyl3_aer(isize,itype,cw_phase) = p_axyl3cwj
	  lptr_orgtol1_aer(isize,itype,cw_phase) = p_atol1cwj
	  lptr_orgtol2_aer(isize,itype,cw_phase) = p_atol2cwj	
	  lptr_orgtol3_aer(isize,itype,cw_phase) = p_atol3cwj	
	  lptr_orgbnz1_aer(isize,itype,cw_phase) = p_abnz1cwj	
	  lptr_orgbnz2_aer(isize,itype,cw_phase) = p_abnz2cwj
	  lptr_orgbnz3_aer(isize,itype,cw_phase) = p_abnz3cwj	
	  lptr_orgtrp1_aer(isize,itype,cw_phase) = p_atrp1cwj	
	  lptr_orgtrp2_aer(isize,itype,cw_phase) = p_atrp2cwj
	  lptr_orgiso1_aer(isize,itype,cw_phase) = p_aiso1cwj	
	  lptr_orgiso2_aer(isize,itype,cw_phase) = p_aiso2cwj	
	  lptr_orgsqt_aer( isize,itype,cw_phase) = p_asqtcwj
	  lptr_orgISO3_aer(isize,itype,cw_phase) = p_aISO3cwj
	  lptr_orgOLGA_aer(isize,itype,cw_phase) = p_aOLGAcwj
	  lptr_orgOLGB_aer(isize,itype,cw_phase) = p_aOLGBcwj	  	  	  
	  lptr_aorgc_aer(  isize,itype,cw_phase) = p_aorgccwj
	  lptr_aorgpa_aer( isize,itype,cw_phase) = p_aorgpacwj	  	
	  lptr_aec_aer(    isize,itype,cw_phase) = p_aeccwj
	  lptr_aothr_aer(  isize,itype,cw_phase) = p_aothrcwj		
	  lptr_ana_aer(    isize,itype,cw_phase) = p_anacwj
	  lptr_acl_aer(    isize,itype,cw_phase) = p_aclcwj
	  lptr_asoil_aer(  isize,itype,cw_phase) = p_asoilcwj          	  
	  do_cloudchem_aer(isize,itype) = .true.	
	endif

! coarse mode
	itype = 3
	isize = 1
	ncomp_aer(itype) = 7
	
	numptr_aer(     isize,itype,ai_phase) = p_corn
	lptr_aso4_aer(  isize,itype,ai_phase) = p_aso4ak
	lptr_ano3_aer(  isize,itype,ai_phase) = p_ano3ak
	lptr_anh4_aer(  isize,itype,ai_phase) = p_anh4ak
	lptr_ana_aer(   isize,itype,ai_phase) = p_anaak
	lptr_acl_aer(   isize,itype,ai_phase) = p_aclak
	lptr_acors_aer( isize,itype,ai_phase) = p_acorsak
	lptr_asoil_aer( isize,itype,ai_phase) = p_asoilak
				
! aerosol in cloud water
        if(cw_phase.gt.0)then	
	  numptr_aer(     isize,itype,cw_phase) = p_corncw
	  lptr_aso4_aer(  isize,itype,cw_phase) = p_aso4cwk
	  lptr_ano3_aer(  isize,itype,cw_phase) = p_ano3cwk
	  lptr_anh4_aer(  isize,itype,cw_phase) = p_anh4cwk
	  lptr_ana_aer(   isize,itype,cw_phase) = p_anacwk
	  lptr_acl_aer(   isize,itype,cw_phase) = p_aclcwk
	  lptr_acors_aer( isize,itype,cw_phase) = p_acorscwk
	  lptr_asoil_aer( isize,itype,cw_phase) = p_asoilcwk		
! no cloudchem for coarse mode because it has no so4/nh4/no3 species
!yes in cmaq as we have so4 etc. (Shaocai Yu)
	  do_cloudchem_aer(isize,itype) = .true.
	endif

	massptr_aer(:,:,:,:) = -999888777
	mastercompptr_aer(:,:) = -999888777

!------------------------------------------
!     param_first_scalar (???)=1
!------------------------------------------

	p1st = param_first_scalar

	do iphase=1,nphase_aer  !1,2
	do itype=1,ntype_aer    !1,3
	do n = 1, nsize_aer(itype)  !1,1
	    ll = 0
	    if (lptr_aso4_aer(n,itype,iphase) .ge. p1st) then
		ll = ll + 1
		massptr_aer(ll,n,itype,iphase) = lptr_aso4_aer(n,itype,iphase)
		mastercompptr_aer(ll,itype) = 1
	    end if
	    if (lptr_ano3_aer(n,itype,iphase) .ge. p1st) then
		ll = ll + 1
		massptr_aer(ll,n,itype,iphase) = lptr_ano3_aer(n,itype,iphase)
		mastercompptr_aer(ll,itype) = 2
	    end if
	    if (lptr_anh4_aer(n,itype,iphase) .ge. p1st) then
		ll = ll + 1
		massptr_aer(ll,n,itype,iphase) = lptr_anh4_aer(n,itype,iphase)
		mastercompptr_aer(ll,itype) = 3
	    end if
	    if (lptr_orgalk_aer(n,itype,iphase) .ge. p1st) then
		ll = ll + 1
		massptr_aer(ll,n,itype,iphase) = lptr_orgalk_aer(n,itype,iphase)
		mastercompptr_aer(ll,itype) = 4
	    end if	    
	    if (lptr_orgxyl1_aer(n,itype,iphase) .ge. p1st) then
		ll = ll + 1
		massptr_aer(ll,n,itype,iphase) = lptr_orgxyl1_aer(n,itype,iphase)
		mastercompptr_aer(ll,itype) = 5
	    end if
	    if (lptr_orgxyl2_aer(n,itype,iphase) .ge. p1st) then
		ll = ll + 1
		massptr_aer(ll,n,itype,iphase) = lptr_orgxyl2_aer(n,itype,iphase)
		mastercompptr_aer(ll,itype) = 6
	    end if
	    if (lptr_orgxyl3_aer(n,itype,iphase) .ge. p1st) then
		ll = ll + 1
		massptr_aer(ll,n,itype,iphase) = lptr_orgxyl3_aer(n,itype,iphase)
		mastercompptr_aer(ll,itype) = 7
	    end if
	    if (lptr_orgtol1_aer(n,itype,iphase) .ge. p1st) then
		ll = ll + 1
		massptr_aer(ll,n,itype,iphase) = lptr_orgtol1_aer(n,itype,iphase)
		mastercompptr_aer(ll,itype) = 8
	    end if
	    if (lptr_orgtol2_aer(n,itype,iphase) .ge. p1st) then
		ll = ll + 1
		massptr_aer(ll,n,itype,iphase) = lptr_orgtol2_aer(n,itype,iphase)
		mastercompptr_aer(ll,itype) = 9
	    end if
	    if (lptr_orgtol3_aer(n,itype,iphase) .ge. p1st) then
		ll = ll + 1
		massptr_aer(ll,n,itype,iphase) = lptr_orgtol3_aer(n,itype,iphase)
		mastercompptr_aer(ll,itype) = 10
	    end if
	    if (lptr_orgbnz1_aer(n,itype,iphase) .ge. p1st) then
		ll = ll + 1
		massptr_aer(ll,n,itype,iphase) = lptr_orgbnz1_aer(n,itype,iphase)
		mastercompptr_aer(ll,itype) = 11
	    end if	    
	    if (lptr_orgbnz2_aer(n,itype,iphase) .ge. p1st) then
		ll = ll + 1
		massptr_aer(ll,n,itype,iphase) = lptr_orgbnz2_aer(n,itype,iphase)
		mastercompptr_aer(ll,itype) = 12
	    end if	    
	    if (lptr_orgbnz3_aer(n,itype,iphase) .ge. p1st) then
		ll = ll + 1
		massptr_aer(ll,n,itype,iphase) = lptr_orgbnz3_aer(n,itype,iphase)
		mastercompptr_aer(ll,itype) = 13
	    end if	    
	    if (lptr_orgtrp1_aer(n,itype,iphase) .ge. p1st) then
		ll = ll + 1
		massptr_aer(ll,n,itype,iphase) = lptr_orgtrp1_aer(n,itype,iphase)
		mastercompptr_aer(ll,itype) = 14
	    end if
	    if (lptr_orgtrp2_aer(n,itype,iphase) .ge. p1st) then
		ll = ll + 1
		massptr_aer(ll,n,itype,iphase) = lptr_orgtrp2_aer(n,itype,iphase)
		mastercompptr_aer(ll,itype) = 15
	    end if
	    if (lptr_orgiso1_aer(n,itype,iphase) .ge. p1st) then
		ll = ll + 1
		massptr_aer(ll,n,itype,iphase) = lptr_orgiso1_aer(n,itype,iphase)
		mastercompptr_aer(ll,itype) = 16
	    end if	
	    if (lptr_orgiso2_aer(n,itype,iphase) .ge. p1st) then
		ll = ll + 1
		massptr_aer(ll,n,itype,iphase) = lptr_orgiso2_aer(n,itype,iphase)
		mastercompptr_aer(ll,itype) = 17
	    end if	    
	    if (lptr_orgsqt_aer(n,itype,iphase) .ge. p1st) then
		ll = ll + 1
		massptr_aer(ll,n,itype,iphase) = lptr_orgsqt_aer(n,itype,iphase)
		mastercompptr_aer(ll,itype) = 18
	    end if	    	    
	    if (lptr_orgISO3_aer(n,itype,iphase) .ge. p1st) then
		ll = ll + 1
		massptr_aer(ll,n,itype,iphase) = lptr_orgISO3_aer(n,itype,iphase)
		mastercompptr_aer(ll,itype) = 19
	    end if
	    if (lptr_orgOLGA_aer(n,itype,iphase) .ge. p1st) then
		ll = ll + 1
		massptr_aer(ll,n,itype,iphase) = lptr_orgOLGA_aer(n,itype,iphase)
		mastercompptr_aer(ll,itype) = 20
	    end if
	    if (lptr_orgOLGB_aer(n,itype,iphase) .ge. p1st) then
		ll = ll + 1
		massptr_aer(ll,n,itype,iphase) = lptr_orgOLGB_aer(n,itype,iphase)
		mastercompptr_aer(ll,itype) = 21
	    end if
	    if (lptr_aorgc_aer(n,itype,iphase) .ge. p1st) then
		ll = ll + 1
		massptr_aer(ll,n,itype,iphase) = lptr_aorgc_aer(n,itype,iphase)
		mastercompptr_aer(ll,itype) = 22
	    end if	    
	    if (lptr_aorgpa_aer(n,itype,iphase) .ge. p1st) then
		ll = ll + 1
		massptr_aer(ll,n,itype,iphase) = lptr_aorgpa_aer(n,itype,iphase)
		mastercompptr_aer(ll,itype) = 23
	    end if	    	    
	    if (lptr_aec_aer(n,itype,iphase) .ge. p1st) then
		ll = ll + 1
		massptr_aer(ll,n,itype,iphase) = lptr_aec_aer(n,itype,iphase)
		mastercompptr_aer(ll,itype) = 24
	    end if
	    if (lptr_aothr_aer(n,itype,iphase) .ge. p1st) then
		ll = ll + 1
		massptr_aer(ll,n,itype,iphase) = lptr_aothr_aer(n,itype,iphase)
		mastercompptr_aer(ll,itype) = 25
	    end if
	    if (lptr_ana_aer(n,itype,iphase) .ge. p1st) then
		ll = ll + 1
		massptr_aer(ll,n,itype,iphase) = lptr_ana_aer(n,itype,iphase)
		mastercompptr_aer(ll,itype) = 26
	    end if
	    if (lptr_acl_aer(n,itype,iphase) .ge. p1st) then
		ll = ll + 1
		massptr_aer(ll,n,itype,iphase) = lptr_acl_aer(n,itype,iphase)
		mastercompptr_aer(ll,itype) = 27				    	    
	    if (lptr_acors_aer(n,itype,iphase) .ge. p1st) then
		ll = ll + 1
		massptr_aer(ll,n,itype,iphase) = lptr_acors_aer(n,itype,iphase)
		mastercompptr_aer(ll,itype) = 28
	    end if
	    if (lptr_asoil_aer(n,itype,iphase) .ge. p1st) then
		ll = ll + 1
		massptr_aer(ll,n,itype,iphase) = lptr_asoil_aer(n,itype,iphase)
		mastercompptr_aer(ll,itype) = 29
	    end if
		
            endif
	    ncomp_aer_nontracer(itype) = ll

	    ncomp_aer(itype) = ll

	    mprognum_aer(n,itype,iphase) = 0
	    if (numptr_aer(n,itype,iphase) .ge. p1st) then
		mprognum_aer(n,itype,iphase) = 1
	    end if
!            write(6,*)'massptr_aer(ll,n,itype,iphase)',massptr_aer(ll,n,itype,iphase),&
!	               ll,n,itype,iphase
!            write(6,*)'mastercompptr_aer(ll,itype)',mastercompptr_aer(ll,itype),ll,itype
	    
	    
! if (maxval(massptr_aer(:,n,itype,iphase)) > 42) then
! print *, ' ==d== massptr_aer ', n, itype, iphase, minval(massptr_aer(:,n,itype,iphase)), maxval(massptr_aer(:,n,itype,iphase)), maxloc(massptr_aer(:,n,itype,iphase))
! end if

	end do ! size
	end do ! type
	end do ! phase

        waterptr_aer(:,:) = 0

	do itype=1,ntype_aer      !ntype_aer = 3
	do ll=1,ncomp_aer(itype)
	   dens_aer(ll,itype) = dens_mastercomp_aer(mastercompptr_aer(ll,itype))
	   mw_aer(ll,itype) = mw_mastercomp_aer(mastercompptr_aer(ll,itype))
	   hygro_aer(ll,itype) = hygro_mastercomp_aer(mastercompptr_aer(ll,itype))
	   name_aer(ll,itype) = name_mastercomp_aer(mastercompptr_aer(ll,itype))
	  
!           write(6,*)'itype,ll,mastercompptr_aer(ll,itype),dens_aer(ll,itype),   &
!	         mw_aer(ll,itype),hygro_aer(ll,itype),name_aer(ll,itype)'
!           write(6,*)itype,ll,mastercompptr_aer(ll,itype),dens_aer(ll,itype),   &
!	         mw_aer(ll,itype),hygro_aer(ll,itype),name_aer(ll,itype)	   
	   
	   
	end do
	end do

	is_aerosol(:) = .false.
	do iphase=1,nphase_aer  ! 2
	do itype=1,ntype_aer    !   3
	do n = 1, nsize_aer(itype)   ! 1
	    do ll = 1, ncomp_aer(itype)  ! 8, 28, 7
	      is_aerosol(massptr_aer(ll,n,itype,iphase))=.true.
	    end do
	      is_aerosol(numptr_aer(n,itype,iphase))=.true.
	end do ! size
	end do ! type
	end do ! phase
	
    
!------------------------------------------------------------------        
! for sectional
!    the dhi/dlo_sect are the upper/lower bounds for 
!       mean-volume diameter for a section/bin
! for modal
!    they should be set to reasonable upper/lower
!       bounds for mean-volume diameters of each modes
!    they are primarily used to put reasonable bounds
!       on number (in relation to mass/volume)
!    the dcen_sect are used by initwet for the impaction scavenging
!       lookup tables, and should represent a "base" mean-volume diameter
! dp_meanvol_tmp (below) is the made-sorgam default initial value 
!    for mean-volume diameter (in cm)
! terminology:  (pi/6) * (mean-volume diameter)**3 ==
!	(volume mixing ratio of section/mode)/(number mixing ratio)
!---------------------------------------------------------------------

!*** Compute these once and they will all  be saved in COMMON
        xxlsg(:,:,:,1) = log(sig_1(:,:,:))
        xxlsg(:,:,:,2) = log(sig_2(:,:,:))
        xxlsg(:,:,:,3) = log(sig_3(:,:,:))
        l2sig(:,:,:,:) = xxlsg(:,:,:,:)**2
 	  
        cdhi_sect(:,:,:,:,:) = 0.0
        cdlo_sect(:,:,:,:,:) = 0.0
!--------------------------------------------------
!       aitken mode
!-------------------------------------------------
        itype = 1
        isize = 1
        csigmag_aer(:,:,:,isize,itype)=sig_1(:,:,:)
        cdp_meanvol_tmp(:,:,:) = 1.0e2*dgn_1(:,:,:)*exp(1.5*l2sig(:,:,:,1)) ! aitken
        cdcen_sect(:,:,:,isize,itype) = cdp_meanvol_tmp(:,:,:)
        cdhi_sect(:,:,:,isize,itype)  = cdp_meanvol_tmp(:,:,:)*4.0
        cdlo_sect(:,:,:,isize,itype)  = cdp_meanvol_tmp(:,:,:)/4.0	
	
if ((its <= 17) .and. (17 <= ite) .and. (jts <= 20) .and. (20 <= jte)) then
	   write(6, *) 'aitken===sig_1(17,2,20),dgn_1(17,2,20), &
	    isize,itype,csigmag_aer,cdcen_sect,cdhi_sect,cdlo_sect', & 
	   sig_1(17,2,20),dgn_1(17,2,20),isize,itype,csigmag_aer(17,2,20,isize,itype), & 
	       cdcen_sect(17,2,20,isize,itype),&
	       cdhi_sect(17,2,20,isize,itype),&
	       cdlo_sect(17,2,20,isize,itype)
end if
	       
!--------------------------------------------------
!       accum mode
!------------------------------------------------- 		
        itype = 2
        isize = 1
        csigmag_aer(:,:,:,isize,itype)=sig_2(:,:,:)
        cdp_meanvol_tmp(:,:,:) = 1.0e2*dgn_2(:,:,:)*exp(1.5*l2sig(:,:,:,2)) ! accum
        cdcen_sect(:,:,:,isize,itype) = cdp_meanvol_tmp(:,:,:)
        cdhi_sect(:,:,:,isize,itype)  = cdp_meanvol_tmp(:,:,:)*4.0
        cdlo_sect(:,:,:,isize,itype)  = cdp_meanvol_tmp(:,:,:)/4.0
	
if ((its <= 17) .and. (17 <= ite) .and. (jts <= 20) .and. (20 <= jte)) then
	   write(6, *) 'accum===sig_2(17,2,20),dgn_2(17,2,20), &
	      isize,itype,csigmag_aer,cdcen_sect,cdhi_sect,cdlo_sect', & 
	      sig_2(17,2,20),dgn_2(17,2,20),isize,itype,csigmag_aer(17,2,20,isize,itype), & 
	       cdcen_sect(17,2,20,isize,itype),&
	       cdhi_sect(17,2,20,isize,itype),&
	       cdlo_sect(17,2,20,isize,itype)	
end if
            	
!--------------------------------------------------
!       coarse mode
!-------------------------------------------------	
        itype = 3
        isize = 1
        csigmag_aer(:,:,:,isize,itype)=sig_3(:,:,:)
        cdp_meanvol_tmp(:,:,:) = 1.0e2*dgn_3(:,:,:)*exp(1.5*l2sig(:,:,:,3)) ! coarse
        cdcen_sect(:,:,:,isize,itype) = cdp_meanvol_tmp(:,:,:)
        cdhi_sect(:,:,:,isize,itype)  = cdp_meanvol_tmp(:,:,:)*4.0
        cdlo_sect(:,:,:,isize,itype)  = cdp_meanvol_tmp(:,:,:)/4.0

if ((its <= 17) .and. (17 <= ite) .and. (jts <= 20) .and. (20 <= jte)) then
	   write(6, *) 'coarse===,sig_3(17,2,20),dgn_3(17,2,20),&
	      isize,itype,csigmag_aer,cdcen_sect,cdhi_sect,cdlo_sect', & 
	      sig_3(17,2,20),dgn_3(17,2,20),isize,itype,csigmag_aer(17,2,20,isize,itype), & 
	       cdcen_sect(17,2,20,isize,itype),&
	       cdhi_sect(17,2,20,isize,itype),&
	       cdlo_sect(17,2,20,isize,itype)
end if

        do itype = 1, ntype_aer
        do isize = 1, nsize_aer(itype)
           cvolumcen_sect(:,:,:,isize,itype) = (pirs/6.0)*(cdcen_sect(:,:,:,isize,itype)**3)
           cvolumlo_sect(:,:,:,isize,itype)  = (pirs/6.0)*(cdlo_sect(:,:,:,isize,itype)**3)
           cvolumhi_sect(:,:,:,isize,itype)  = (pirs/6.0)*(cdhi_sect(:,:,:,isize,itype)**3)
        end do
        end do    
    
!---End for initialization for aerosol-cloud interaction

!             goto 999   ! (OK before here)

        do p=1,maxd_aphase  !maxd_aphase=2, interstitial and cloud-activited
        do n=1,ntype_aer !ntype_aer= 3 
        do m=1,nsize_aer(n)
!	   write(6,*)'numptr_aer( m, n, p )',numptr_aer( m, n, p ),m, n, p	
	   do l=1,ncomp_aer(n)  !ncomp_aer(1)=8,ncomp_aer(2)=28,ncomp_aer(3)=7	
! 		write(6,*)'massptr_aer(l, m, n, p)',massptr_aer(l, m, n, p),l,m,n,p
!		write(6,*)'its:ite,kts:kte,jts:jte',its,ite,kts,kte,jts,jte
!		write(6,*)'ims:ime,kms:kme,jms:jme',ims,ime,kms,kme,jms,jme
	 end do ! ncomp_aer(n)
	 end do ! size
	 end do ! type
	 end do ! phase

        do p=1,maxd_aphase  !maxd_aphase=2, interstitial and cloud-activited
        do n=1,ntype_aer !ntype_aer= 3 
        do m=1,nsize_aer(n)
	   do l=1,ncomp_aer(n)  !ncomp_aer(1)=8,ncomp_aer(2)=28,ncomp_aer(3)=7
     
!-----------------------------------------------------------------------------------
!	      massptr_aer(l, m, n, p)=ptr   !  need do something here???? for cmaq
!-------------------------------------------------------------------------------------
! maer is ug/kg-air (mixing ratio);  naer is #/kg-air;  dgnum is cm;  dens_aer is g/cm3
! 1.e6 factor converts g to ug
! 	      maer= 1.0e6 * naer * dens_aer(l,n) * ( (3.1416/6.) *   &
!                   (dgnum_aer(m,n)**3) * exp( 4.5*((log(sigmag_aer(m,n)))**2) ) )
!------------------------------------------------------------------------------------
!	      if(p.eq.ai_phase)then
!	         chem(its:ite,kts:kte,jts:jte,ptr)=ae_mass(its:ite,kts:kte,jts:jte,ptr)
!	      else
!	         chem(its:ite,kts:kte,jts:jte,ptr)=0.   ???????
!	      endif
!---------------------------------------------------------------------------------	      
	      if(p.eq.ai_phase)then
	      

               do  j=jts,jte 
               do  k=kts,kte
               do  i=its,ite	

!         if((massptr_aer(l, m, n, p).eq.1).and.(n.eq.2).and.(m.eq.2))then
!           write(*,2999) ' ==d== cfma ', i, j, k,n,l,cfrac_mass(i,k,j,massptr_aer(l, m, n, p))
!         endif
	 
 2999   format(A13,8x,I3,4x,I3,4x,I3,4x,I4,4x,I4,4x,E12.6)


	       		      
! if (massptr_aer(l, m, n, p) > 42) then
! print *, ' ==d== Cmixactivate chem 00 ', massptr_aer(l, m, n, p), l, m, n, p
! end if
	         chem(i,k,j,massptr_aer(l, m, n, p))=          &    !???
		       ae_mass(i,k,j,massptr_aer(l, m, n, p))  &
		       *cfrac_mass(i,k,j,massptr_aer(l, m, n, p))
		       
!                if(i.eq.152.and.k.eq.23.and.j.eq.189.and.massptr_aer(l, m, n, p).eq.4)then
!		  write(6,*)'Here in cmaq_mixa',i,k,j 		 
! 		  write(6,*)'massptr_aer(l, m, n, p)',massptr_aer(l, m, n, p),l,m,n,p
!		  write(6,*)'chem',chem(17,2,20,massptr_aer(l, m, n, p)),chem(i,k,j,massptr_aer(l, m, n, p))
!		  write(6,*)'ae_mass',ae_mass(17,2,20,massptr_aer(l, m, n, p))
!		  write(6,*)'cfrac_mass(i,k,j,massptr_aer(l, m, n, p)',   &
!		       cfrac_mass(i,k,j,massptr_aer(l, m, n, p)),l, m, n, p
		       
!		  write(6,*)'its:ite,kts:kte,jts:jte',its,ite,kts,kte,jts,jte
!		  write(6,*)'ims:ime,kms:kme,jms:jme',ims,ime,kms,kme,jms,jme
!                endif
		enddo
		enddo
		enddo
 
		
!		stop

!        goto 999   !  (OK)
	      else
	      
	      
!	         chem(its:ite,kts:kte,jts:jte,massptr_aer(l, m, n, p))=0.    ! cloud activated
		 
               do  j=jts,jte 
               do  k=kts,kte
               do  i=its,ite	
	       		      
	         chem(i,k,j,massptr_aer(l, m, n, p))=          &           !???
		       ae_mass(i,k,j,massptr_aer(l, m, n, p)-p_corn)  &
		       *(1.0-cfrac_mass(i,k,j,massptr_aer(l, m, n, p)-p_corn))

				       
!                if(i.eq.17.and.k.eq.2.and.j.eq.20)then		 
! 		  write(6,*)'massptr_aer(l, m, n, p)',massptr_aer(l, m, n, p),l,m,n,p
!		  write(6,*)'chem',chem(17,2,20,massptr_aer(l, m, n, p))
!		  write(6,*)'ae_mass',ae_mass(17,2,20,massptr_aer(l, m, n, p))
!		  write(6,*)'cfrac_mass(i,k,j,massptr_aer(l, m, n, p)',  &
!		      cfrac_mass(i,k,j,massptr_aer(l, m, n, p)),l, m, n, p
!		  write(6,*)'its:ite,kts:kte,jts:jte',its,ite,kts,kte,jts,jte
!		  write(6,*)'ims:ime,kms:kme,jms:jme',ims,ime,kms,kme,jms,jme
!                endif
		enddo
		enddo
		enddo		 
		 
	      endif	      
!        goto 999   !  	(OK)      
	      
	      
	 end do ! ncomp_aer(n)
	 end do ! size
	 end do ! type
	 end do ! phase
	 
!       ptr=0

!         goto 999 ! (ok)

       do p=1,nphase_aer
       do n=1,ntype_aer   !ntype_aer= 3
       do m=1,nsize_aer(n)
!-------------------------------------------------------------------------       
!          numptr_aer( m, n, p )=ptr   !need do something here ??????
!	  if(p.eq.ai_phase)then
!	     chem(its:ite,kts:kte,jts:jte,ptr)=ae_num(its:ite,kts:kte,jts:jte,ptr)  !interstitial 
!	  else
!	     chem(its:ite,kts:kte,jts:jte,ptr)=0.    ! cloud activated
!	  endif
!--------------------------------------------------------------------------------
          if(p.eq.ai_phase)then
	  	 
            do  j=jts,jte 
            do  k=kts,kte
            do  i=its,ite	
!         if(n.eq.2)then
!           write(*,1999) ' ==d== cfnu ', i, j, k,n,cfrac_num(i,k,j,numptr_aer( m, n, p ))
!         endif

 1999   format(A13,8x,I3,4x,I3,4x,I3,4x,I4,4x,E12.6)
	  	  
             if (n == 1) then
                chem(i,k,j,numptr_aer( m, n, p ))=ae_num_1(i,k,j) &
		    *cfrac_num(i,k,j,numptr_aer( m, n, p ))  !interstitial
             else if (n == 2) then
                chem(i,k,j,numptr_aer( m, n, p ))=ae_num_2(i,k,j) &
		    *cfrac_num(i,k,j,numptr_aer( m, n, p ))  !interstitial
             else if (n == 3) then
                chem(i,k,j,numptr_aer( m, n, p ))=ae_num_3(i,k,j) &
		    *cfrac_num(i,k,j,numptr_aer( m, n, p ))  !interstitial
             end if
	     
!             if(i.eq.17.and.k.eq.2.and.j.eq.20)then
!	        write(6,*)'p',p
! 		write(6,*)'chem,ae_num_1 ',chem(17,2,20,numptr_aer( m, n, p )),ae_num_1(17,2,20)
! 		write(6,*)'chem,ae_num_2 ',chem(17,2,20,numptr_aer( m, n, p )),ae_num_2(17,2,20)
! 		write(6,*)'chem,ae_num_3 ',chem(17,2,20,numptr_aer( m, n, p )),ae_num_3(17,2,20)	       
!                write(6,*)'cfrac_num',cfrac_num(i,k,j,numptr_aer( m, n, p ))
!	     endif
	       
            enddo
            enddo		
            enddo
		
          else
!             chem(its:ite,kts:kte,jts:jte,numptr_aer( m, n, p ))=0.    ! cloud activated
	     
            do  j=jts,jte 
            do  k=kts,kte
            do  i=its,ite	
	  	  
             if (n == 1) then
                chem(i,k,j,numptr_aer( m, n, p ))=ae_num_1(i,k,j) &
		    *(1.0-cfrac_num(i,k,j,numptr_aer( m, n, p )-p_corn))  !cloud activated
             else if (n == 2) then
                chem(i,k,j,numptr_aer( m, n, p ))=ae_num_2(i,k,j) &
		    *(1.0-cfrac_num(i,k,j,numptr_aer( m, n, p )-p_corn))  !cloud activated
             else if (n == 3) then
                chem(i,k,j,numptr_aer( m, n, p ))=ae_num_3(i,k,j) &
		    *(1.0-cfrac_num(i,k,j,numptr_aer( m, n, p )-p_corn))  !cloud activated
             end if
!             if(i.eq.17.and.k.eq.2.and.j.eq.20)then
!	        write(6,*)'p',p
! 		write(6,*)'chem,ae_num_1 ',chem(17,2,20,numptr_aer( m, n, p )),ae_num_1(17,2,20)
! 		write(6,*)'chem,ae_num_2 ',chem(17,2,20,numptr_aer( m, n, p )),ae_num_2(17,2,20)
! 		write(6,*)'chem,ae_num_3 ',chem(17,2,20,numptr_aer( m, n, p )),ae_num_3(17,2,20)	       
!                write(6,*)'cfrac_num',cfrac_num(i,k,j,numptr_aer( m, n, p ))
!	     endif
	     
            enddo
            enddo		
            enddo	     
	     
	     
          endif
 		
	  	  
	end do ! size
	end do ! type
	end do ! phase	 
	
!           goto 999 ! (ok) 
	   
	   
	   	
       do n=1,ntype_aer
       do m=1,nsize_aer(n)      
          if(m.eq.1)then
	    waterptr_aer(m, n)= P_H2Onu
	  else if(m.eq.2)then 
	    waterptr_aer(m, n)= P_H2Oac	  
	  else if(m.eq.3)then
	    waterptr_aer(m, n)= P_H2Ocorn	  
	  endif
	  
!	  waterptr_aer(m, n)=-1
	end do ! size
	end do ! type
	ddvel(its:ite,jts:jte,:)=0.	 	 
	 
!-------------------------------------------------------
  
   
! calculate volume-weighted bulk hygroscopicity for each type and size
!-------------------------------------------------------------------
!hygro_aer(l,n): dimensionless, hydroscopicity for each component
!chem(i,k,j,massptr_aer(l,m,n,ai_phase)): mass mixing ratio, ug/kg-air (mixing ratio);  
!dens_aer(l,n):  dry aerosol densities (g/cm3)
!---------------------------------------------------------------------

      do 100 j=jts,jte 
      do 100 k=kts,kte
      do 100 i=its,ite

       do n=1,ntype_aer   !ntype_aer= 3
       do m=1,nsize_aer(n)
	       sumhygro=0
	       sumvol=0
	       do l=1,ncomp_aer(n)   !ncomp_aer(1)=8,ncomp_aer(2)=28,ncomp_aer(3)=7
	          sumhygro = sumhygro+hygro_aer(l,n)*   &
                   chem(i,k,j,massptr_aer(l,m,n,ai_phase))/dens_aer(l,n)
	          sumvol = sumvol+chem(i,k,j,massptr_aer(l,m,n,ai_phase))/dens_aer(l,n)
		  
!		  if((i.eq.17).and.(j.eq.20).and.(k.eq.2))then
! 		    write(6,*)'sumhygro',hygro_aer(l,n),i,k,j,massptr_aer(l,m,n,ai_phase)
!		    write(6,*)'chem(i,k,j,massptr_aer(l,m,n,ai_phase))',  &
!		          chem(i,k,j,massptr_aer(l,m,n,ai_phase))
!		    write(6,*)'dens_aer(l,n)',dens_aer(l,n)
!		  endif
		  
	       end do ! comp
               hygro(i,k,j,m,n)=sumhygro/sumvol
	       
	       if((i.eq.17).and.(j.eq.20).and.(k.eq.2))then
!	       if(hygro(i,k,j,m,n).eq.0.0)then
		  write(6,*)'in ---sumhygro,sumvol',sumhygro,sumvol
		  write(6,*)'in ---hygro(i,k,j,m,n)',hygro(i,k,j,m,n),i,k,j,m,n,sumhygro,sumvol	       
	       endif	       
	       
	end do ! size
	end do ! type
  100 continue

      hygro_1=hygro(:,:,:,1,1)
      hygro_2=hygro(:,:,:,1,2)
      hygro_3=hygro(:,:,:,1,3)

!            write(6,*)'in CMAQ_mixactivate before cmaq_icenucleat'  
!             goto 999 ! (ok)
	           
!-------------------------------------------------------------------------
!   CMAQ aerosol ice nucleation
!-------Shaocai
!       do i=its,ite 
!do k=kds+1,kde 
!       do j=jts,jte
	
!	     f_rain_phy(i,k,j)=0.2
!	     f_ice_phy(i,k,j)=0.2

!          if(f_ice_phy(i,k,j).lt.1.0E-10)then
!	     f_ice_phy(i,k,j)=0.0
!	  endif
	  
!          if(qnice3d(i,k,j).eq.0.0)then
!	     f_ice_phy(i,k,j)=0.0
!	  endif
!	  
!          if(f_rain_phy(i,k,j).lt.1.0E-10)then
!	     f_rain_phy(i,k,j)=0.0
!	  endif
!	  	  	      
!          if(f_ice_phy(i,k,j).ne.0.0)then
!	   if(i.eq.its.and.k.eq.(kds+1).and.j.eq.jts)then
!          write(6, *) 'in cmaq_mixactivate before calling cmaq_icenucleat  &
!	           qndrop3d,i,k,j',    &
!                       qndrop3d(i,k,j),i,k,j  
		       
!	   write(6, *) 'msectional',msectional 
!           write(6, *) 'chem',chem 
!	   write(6, *) 'loc_num_chem',loc_num_chem 
!	   write(6, *) 'qv, qc, qi, qndrop3d',qv, qc, qi,qs, qnice3d  
!           write(6, *) 't_phy, w, ddvel, idrydep_onoff',t_phy, w, ddvel, idrydep_onoff  
!           write(6, *) 'maxd_acomp, maxd_asize',maxd_acomp, maxd_asize 
!	   write(6, *) 'maxd_atype, maxd_aphase',maxd_atype, maxd_aphase
!           write(6, *) 'ncomp_aer, nsize_aer',ncomp_aer, nsize_aer 
!	   write(6, *) 'ntype_aer, nphase_aer',ntype_aer, nphase_aer 
!           write(6, *) 'numptr_aer, massptr_aer',numptr_aer, massptr_aer 
!	   write(6, *) 'cdlo_sect, cdhi_sect',cdlo_sect, cdhi_sect 
!	   write(6, *) 'csigmag_aer, cdcen_sect',csigmag_aer, cdcen_sect
!           write(6, *) 'dens_aer, mw_aer',dens_aer, mw_aer
!            write(6, *) 'waterptr_aer, hygro',waterptr_aer, hygro 
!	   write(6, *) 'ai_phase, cw_phase',ai_phase, cw_phase
!           write(6, *) 'nihf3d,niimm3d,nidep3d,nimey3d',nihf3d,niimm3d,nidep3d,nimey3d 
!           write(6, *) 'rndst1,rndst2,rndst3,rndst4',rndst1,rndst2,rndst3,rndst4
!	   write(6, *) 'nacon1,nacon2,nacon3,nacon4',nacon1,nacon2,nacon3,nacon4	   	   
!           write(6, *) 'ids,ide, jds,jde, kds,kde',ids,ide, jds,jde, kds,kde                          
!           write(6, *) 'ims,ime, jms,jme, kms,kme',ims,ime, jms,jme, kms,kme                         
!           write(6, *) 'its,ite, jts,jte, kts,kte',its,ite, jts,jte, kts,kte                         
!           write(6, *) 'rho_phy, z, dz8w',rho_phy, z, dz8w 
!	   write(6, *) 'p_at_w, t_at_w, exch_h',p_at_w, t_at_w, exch_h     
!           write(6, *) 'cldfra, cldfra_old, qsrflx',cldfra, cldfra_old, qsrflx                      
!	   write(6, *) 'ccn1, ccn2, ccn3, ccn4',ccn1, ccn2, ccn3, ccn4 
!	   write(6, *) 'ccn5, ccn6, nsource',ccn5, ccn6, nsource      
!	   write(6, *) 'id, ktau, dtstep',id, ktau, dtstep
!           write(6, *) 'f_qc, f_qi',f_qc, f_qi,f_qs,config_flags%ra_sw_physics,f_ice_phy(i,k,j)
!	   write(6, *) ' f_rain_phy',f_rain_phy(i,k,j),qnice3d(i,k,j),config_flags%ra_sw_physics                
!          endif
!	  endif
!enddo
!enddo
!enddo
 
!       write(6,*)'chem(its,kts,jts,2) before cmaq_icenucleat',chem(its,kts,jts,2)
	
!	write(6,*)'qnice3d',qnice3d
!	write(6,*)'cldfra',cldfra 
!	write(6,*)'cldfra_old',cldfra_old
!	 Tqnice3d = qnice3d
!	 Tcldfra=cldfra
!	 Tcldfra_old=cldfra_old
!         Tf_ice_phy=f_ice_phy
!	 Tf_rain_phy=f_rain_phy

      call cmaq_icenucleat (msectional,chem,loc_num_chem,qv,qc,qi,qs,qnice3d,  &
           t_phy,p_phy,w,maxd_acomp,maxd_asize, maxd_atype, maxd_aphase,   &
           ncomp_aer, nsize_aer, ntype_aer, nphase_aer,  &
           numptr_aer, massptr_aer, cdlo_sect, cdhi_sect, csigmag_aer, cdcen_sect,  &
           dens_aer, mw_aer,waterptr_aer, hygro,  ai_phase, cw_phase,  &
            !czhao 
           nihf3d,niimm3d,nidep3d,nimey3d,  &
           rndst1,rndst2,rndst3,rndst4,nacon1,nacon2,nacon3,nacon4,  &
           dgn_1,dgn_2,dgn_3,   &
           ids,ide, jds,jde, kds,kde,                            &
           ims,ime, jms,jme, kms,kme,                            &
           its,ite, jts,jte, kts,kte,                             &
           rho_phy, z, dz8w, p_at_w, t_at_w, exch_h,      &
           cldfra, cldfra_old,id, ktau, dtstep, &
           f_qv,f_qc,f_qi,f_qs,f_ice_phy,f_rain_phy, &
!----------------------------------------------------------------------------
!Added by R Atlas for thompson indirect effect
           config_flags%ra_sw_physics, config_flags%icloud, xland, &
           dx, dy, qnum_wf, qnum_if, &
           morrison_indirect, thompson_indirect)

!-------------------------------------------------------------------------

!          write(6,*)'in CMAQ_mixactivate after cmaq_icenucleat'
! 	   write(6,*)'qnice3d',qnice3d
!	   write(6,*)'cldfra',cldfra 
!	   write(6,*)'cldfra_old',cldfra_old
!   	 qnice3d = Tqnice3d
!	 cldfra=Tcldfra
!	 cldfra_old=Tcldfra_old
!         f_ice_phy=Tf_ice_phy
!	 f_rain_phy=Tf_rain_phy
!---------------------------------------------------------------	
!--------------------------------------------------------------------      

if (morrison_indirect) then ! R Atlas

      call Cmixactivate(  msectional, &
           chem, loc_num_chem, qv, qc, qi, qndrop3d,   &
           t_phy, w, ddvel, idrydep_onoff,  &
           maxd_acomp, maxd_asize, maxd_atype, maxd_aphase,   &
           ncomp_aer, nsize_aer, ntype_aer, nphase_aer,  &
           numptr_aer, massptr_aer, cdlo_sect, cdhi_sect, csigmag_aer, cdcen_sect,  &
           dens_aer, mw_aer,           &
           waterptr_aer, hygro,  ai_phase, cw_phase,                 &
           ids,ide, jds,jde, kds,kde,                            &
           ims,ime, jms,jme, kms,kme,                            &
           its,ite, jts,jte, kts,kte,                            &
           rho_phy, z, dz8w, p_at_w, t_at_w, exch_h,      &
           cldfra, cldfra_old, qsrflx,                      &
	       ccn1, ccn2, ccn3, ccn4, ccn5, ccn6, nsource,       &
	       id, ktau, dtstep, &
           f_qc, f_qi, cfrac_num, cfrac_mass )

end if

 999 continue
 
!	     ccn1(:,:,:)=100.0
!	     ccn2(:,:,:)=100.0
!	     ccn3(:,:,:)=100.0
!	     ccn4(:,:,:)=100.0
!	     ccn5(:,:,:)=100.0
!	     ccn6(:,:,:)=100.0 

  

 10999 continue

!R Atlas ------------------------ 
      hygro_1=hygro(:,:,:,1,1)
      hygro_2=hygro(:,:,:,1,2)
      hygro_3=hygro(:,:,:,1,3)
!--------------------------------
      end subroutine cmaq_mixactivate

!END MODULE module_mixactivate_wrappers
