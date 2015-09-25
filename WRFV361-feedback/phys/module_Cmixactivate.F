!************************************************************************
! This module is for CMAQ aerosol output because aerosol size distributions 
!        in CMAQ are changing with time
! (Shaocai Yu) (12/18/2009)
! 
!************************************************************************

MODULE module_Cmixactivate
! PRIVATE
PUBLIC cmaq_icenucleat
CONTAINS

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!   nov-04 sg ! replaced amode with aer and expanded aerosol dimension to include type and phase

! 06-nov-2005 rce - grid_id & ktau added to arg list
! 25-apr-2006 rce - dens_aer is (g/cm3), NOT (kg/m3)
subroutine Cmixactivate(  msectional,            &
           chem, num_chem, qv, qc, qi, qndrop3d,         &
           temp, w, ddvel, idrydep_onoff,  &
           maxd_acomp, maxd_asize, maxd_atype, maxd_aphase,   &
           ncomp_aer, nsize_aer, ntype_aer, nphase_aer,  &
           numptr_aer, massptr_aer, cdlo_sect, cdhi_sect, csigmag_aer, cdpvolmean_aer,  &
           dens_aer, mw_aer,               &
           waterptr_aer, hygro, ai_phase, cw_phase,              &
           ids,ide, jds,jde, kds,kde,                            &
           ims,ime, jms,jme, kms,kme,                            &
           its,ite, jts,jte, kts,kte,                            &
           rho, zm, dz8w, p_at_w, t_at_w, kvh,      &
           cldfra, cldfra_old, qsrflx,          &
           ccn1, ccn2, ccn3, ccn4, ccn5, ccn6, nsource,       &
           grid_id, ktau, dtstep, &
           f_qc, f_qi, cfrac_num, cfrac_mass)

!     vertical diffusion and nucleation of cloud droplets
!     assume cloud presence controlled by cloud fraction
!     doesn't distinguish between warm, cold clouds

  USE module_model_constants, only: g, rhowater, xlv, cp, rvovrd, r_d, r_v, mwdry, ep_2
  USE module_radiation_driver, only: cal_cldfra2
  USE module_mixactivate

  implicit none

!     input

  INTEGER, intent(in) ::         grid_id, ktau
  INTEGER, intent(in) ::         num_chem
  integer, intent(in) ::         ids,ide, jds,jde, kds,kde,    &
                                 ims,ime, jms,jme, kms,kme,    &
                                 its,ite, jts,jte, kts,kte

  integer, intent(in) :: maxd_aphase, nphase_aer, maxd_atype, ntype_aer
  integer, intent(in) :: maxd_asize, maxd_acomp, nsize_aer(maxd_atype)
  integer, intent(in) ::   &
       ncomp_aer( maxd_atype  ),   &
       massptr_aer( maxd_acomp, maxd_asize, maxd_atype, maxd_aphase ),   &
       waterptr_aer( maxd_asize, maxd_atype ),   &
       numptr_aer( maxd_asize, maxd_atype, maxd_aphase), &
       ai_phase, cw_phase
  integer, intent(in) :: msectional ! 1 for sectional, 0 for modal
  integer, intent(in) :: idrydep_onoff
  
  REAL, intent(in), DIMENSION( ims:ime, kms:kme, jms:jme,maxd_asize, maxd_atype ) :: &
       cdlo_sect,cdhi_sect,csigmag_aer,cdpvolmean_aer
!----------------------------------------------------------------------------
!       cdlo_sect       ! minimum size of section (cm)
!       cdhi_sect       ! maximum size of section (cm)
!       csigmag_aer     ! geometric standard deviation of aerosol size dist
!       cdpvolmean_aer  ! mean-volume diameter (cm) of mode    
!----------------------------------------------------------------------------  
!   local varibles here for CMAQ case (Shaocai Yu)
!------------------------------------------------------------------------- 
  real ::                       &
       dlo_sect( maxd_asize, maxd_atype ),   & ! minimum size of section (cm)
       dhi_sect( maxd_asize, maxd_atype ),   & ! maximum size of section (cm)
       sigmag_aer(maxd_asize, maxd_atype),   & ! geometric standard deviation of aerosol size dist
       dpvolmean_aer(maxd_asize, maxd_atype)   ! mean-volume diameter (cm) of mode
!------------------------------------------------------------------------------
       
   REAL, intent(in) :: &      
       dens_aer( maxd_acomp, maxd_atype),    & ! density (g/cm3) of material
       mw_aer( maxd_acomp, maxd_atype)       ! molecular weight (g/mole)

! terminology:  (pi/6) * (mean-volume diameter)**3 ==
!       (volume mixing ratio of section/mode)/(number mixing ratio)

! David: change the last dimenstion to :
  REAL, intent(inout), DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ) :: &
       chem ! aerosol molar mixing ratio (ug/kg or #/kg)

  REAL, intent(in), DIMENSION( ims:ime, kms:kme, jms:jme ) :: &
       qv, qc, qi ! water species (vapor, cloud drops, cloud ice) mixing ratio (g/g)

  LOGICAL, OPTIONAL :: f_qc, f_qi

  REAL, intent(inout), DIMENSION( ims:ime, kms:kme, jms:jme ) :: &
       qndrop3d    ! water species mixing ratio (g/g)

  real, intent(in) :: dtstep             ! time step for microphysics (s)
  real, intent(in) :: temp(ims:ime, kms:kme, jms:jme)    ! temperature (K)
  real, intent(in) :: w(ims:ime, kms:kme, jms:jme)   ! vertical velocity (m/s)
  real, intent(in) :: rho(ims:ime, kms:kme, jms:jme)    ! density at mid-level  (kg/m3)
  REAL, intent(in) :: ddvel( its:ite, jts:jte, num_chem ) ! deposition velocity  (m/s)
  real, intent(in) :: zm(ims:ime, kms:kme, jms:jme)     ! geopotential height of level (m)
  real, intent(in) :: dz8w(ims:ime, kms:kme, jms:jme) ! layer thickness (m)
  real, intent(in) :: p_at_w(ims:ime, kms:kme, jms:jme) ! pressure at layer interface (Pa)
  real, intent(in) :: t_at_w(ims:ime, kms:kme, jms:jme) ! temperature at layer interface (K)
  real, intent(in) :: kvh(ims:ime, kms:kme, jms:jme)    ! vertical diffusivity (m2/s)
  real, intent(inout) :: cldfra_old(ims:ime, kms:kme, jms:jme)! cloud fraction on previous time step
  real, intent(inout) :: cldfra(ims:ime, kms:kme, jms:jme)    ! cloud fraction
  real, intent(in) :: hygro( its:ite, kts:kte, jts:jte, maxd_asize, maxd_atype ) ! bulk hygroscopicity   &

  REAL, intent(out), DIMENSION( ims:ime, jms:jme, num_chem ) ::   qsrflx ! dry deposition rate for aerosol
  real, intent(out), dimension(ims:ime,kms:kme,jms:jme) :: nsource, &  ! droplet number source (#/kg/s)
       ccn1,ccn2,ccn3,ccn4,ccn5,ccn6  ! number conc of aerosols activated at supersat

  real, intent(inout) :: cfrac_num (ims:ime, kms:kme, jms:jme, 200)
  real, intent(inout) :: cfrac_mass (ims:ime, kms:kme, jms:jme, 200)

!--------------------Local storage-------------------------------------
!
  real :: dgnum_aer(maxd_asize, maxd_atype) ! median diameter (cm) of number distrib of mode
  real :: qndrop(kms:kme)      ! cloud droplet number mixing ratio (#/kg)
  real :: lcldfra(kms:kme)     ! liquid cloud fraction
  real :: lcldfra_old(kms:kme) ! liquid cloud fraction for previous timestep
  real :: wtke(kms:kme)        ! turbulent vertical velocity at base of layer k (m2/s)
  real zn(kms:kme)             ! g/pdel (m2/g) for layer
  real zs(kms:kme)             ! inverse of distance between levels (m)
  real, parameter :: zkmin = 0.01
  real, parameter :: zkmax = 100.
  real cs(kms:kme)             ! air density (kg/m3) at layer center
  real csbot(kms:kme)          ! air density (kg/m3) at layer bottom
  real csbot_cscen(kms:kme)    ! csbot(k)/cs(k)
  real dz(kms:kme)             ! geometric thickness of layers (m)

  real wdiab                   ! diabatic vertical velocity
!      real, parameter :: wmixmin = 0.1 ! minimum turbulence vertical velocity (m/s)
  real, parameter :: wmixmin = 0.2 ! minimum turbulence vertical velocity (m/s)
!      real, parameter :: wmixmin = 1.0 ! minimum turbulence vertical velocity (m/s)
  real :: qndrop_new(kms:kme)  ! droplet number nucleated on cloud boundaries
  real :: ekd(kms:kme)         ! diffusivity for droplets (m2/s)
  real :: ekk(kms:kme)         ! density*diffusivity for droplets (kg/m3 m2/s)
  real :: srcn(kms:kme)        ! droplet source rate (/s)
  real, parameter :: sq2pi = 2.5066282746
  real dtinv

  integer km1,kp1
  real wbar,wmix,wmin,wmax
  real dum
  real tmpa, tmpb, tmpc, tmpc1, tmpc2, tmpd, tmpe, tmpf
  real tmpcourno
  real dact
  real fluxntot         ! (#/cm2/s)
  real fac_srflx
  real depvel_drop, depvel_tmp
  real, parameter :: depvel_uplimit = 1.0 ! upper limit for dep vels (m/s)
  real :: surfrate(num_chem) ! surface exchange rate (/s)
  real surfratemax      ! max surfrate for all species treated here
  real surfrate_drop    ! surfade exchange rate for droplelts
  real dtmin,tinv,dtt
  integer nsubmix,nsubmix_bnd
  integer i,j,k,m,n,nsub
  real dtmix
  real alogarg
  real qcld
  real pi
  integer nnew,nsav,ntemp
  real :: overlapp(kms:kme),overlapm(kms:kme) ! cloud overlap
  real ::  ekkp(kms:kme),ekkm(kms:kme) ! zn*zs*density*diffusivity
!  integer, save :: count_submix(100)=0 ! wig: Note that this is a no-no for tile threads with OMP

  integer lnum,lnumcw,l,lmass,lmasscw,lsfc,lsfccw,ltype,lsig,lwater
  integer :: ntype(maxd_asize)

  real ::  naerosol(maxd_asize, maxd_atype)    ! interstitial aerosol number conc (/m3)
  real ::  naerosolcw(maxd_asize, maxd_atype)  ! activated number conc (/m3)
  real ::   maerosol(maxd_acomp,maxd_asize, maxd_atype)   ! interstit mass conc (kg/m3)
  real ::   maerosolcw(maxd_acomp,maxd_asize, maxd_atype) ! activated mass conc (kg/m3)
  real ::   maerosol_tot(maxd_asize, maxd_atype)     ! species-total interstit mass conc (kg/m3)
  real ::   maerosol_totcw(maxd_asize, maxd_atype)   ! species-total activated mass conc (kg/m3)
  real ::   vaerosol(maxd_asize, maxd_atype) ! interstit+activated aerosol volume conc (m3/m3)
  real ::   vaerosolcw(maxd_asize, maxd_atype) ! activated aerosol volume conc (m3/m3)
  real ::   raercol(kms:kme,num_chem,2) ! aerosol mass, number mixing ratios
  real ::   source(kms:kme) !

  real ::   fn(maxd_asize, maxd_atype)         ! activation fraction for aerosol number
  real ::   fs(maxd_asize, maxd_atype)         ! activation fraction for aerosol sfcarea
  real ::   fm(maxd_asize, maxd_atype)         ! activation fraction for aerosol mass
  integer ::   ncomp(maxd_atype)

  real ::   fluxn(maxd_asize, maxd_atype)      ! number  activation fraction flux (m/s)
  real ::   fluxs(maxd_asize, maxd_atype)      ! sfcarea activation fraction flux (m/s)
  real ::   fluxm(maxd_asize, maxd_atype)      ! mass    activation fraction flux (m/s)
  real ::   flux_fullact(kms:kme)              ! 100%    activation fraction flux (m/s)
!     note:  activation fraction fluxes are defined as
!     fluxn = [flux of activated aero. number into cloud (#/m2/s)]
!           / [aero. number conc. in updraft, just below cloudbase (#/m3)]

  real :: nact(kms:kme,maxd_asize, maxd_atype)  ! fractional aero. number  activation rate (/s)
  real :: mact(kms:kme,maxd_asize, maxd_atype)  ! fractional aero. mass    activation rate (/s)
  real :: npv(maxd_asize, maxd_atype) ! number per volume concentration (/m3)
  real scale

  real :: hygro_aer(maxd_asize, maxd_atype)  ! hygroscopicity of aerosol mode
  real :: exp45logsig     ! exp(4.5*alogsig**2)
  real :: alogsig(maxd_asize, maxd_atype) ! natl log of geometric standard dev of aerosol
  integer, parameter :: psat=6  ! number of supersaturations to calc ccn concentration
  real ccn(kts:kte,psat)        ! number conc of aerosols activated at supersat
  real, parameter :: supersat(psat)= &! supersaturation (%) to determine ccn concentration
       (/0.02,0.05,0.1,0.2,0.5,1.0/)
  real super(psat) ! supersaturation
  real, parameter :: surften = 0.076 ! surface tension of water w/respect to air (N/m)
  real :: ccnfact(psat,maxd_asize, maxd_atype)
  real :: amcube(maxd_asize, maxd_atype) ! cube of dry mode radius (m)
  real :: argfactor(maxd_asize, maxd_atype)
  real aten ! surface tension parameter
  real t0 ! reference temperature
  real sm ! critical supersaturation
  real arg

  integer,parameter :: icheck_colmass = 0
           ! icheck_colmass > 0 turns on mass/number conservation checking
           ! values of 1, 10, 100 produce less to more diagnostics
  integer :: colmass_worst_ij( 2, 0:maxd_acomp, maxd_asize, maxd_atype )
  integer :: colmass_maxworst_i(3)
  real :: colmass_bgn( 0:maxd_acomp, maxd_asize, maxd_atype, maxd_aphase )
  real :: colmass_end( 0:maxd_acomp, maxd_asize, maxd_atype, maxd_aphase )
  real :: colmass_sfc( 0:maxd_acomp, maxd_asize, maxd_atype, maxd_aphase )
  real :: colmass_worst( 0:maxd_acomp, maxd_asize, maxd_atype )
  real :: colmass_maxworst_r
  real :: rhodz( kts:kte ), rhodzsum

!!$#if (defined AIX)
!!$#define ERF erf
!!$#define ERFC erfc
!!$#else
!!$#define ERF erf
!!$    real erf
!!$#define ERFC erfc
!!$    real erfc
!!$#endif

  character*8, parameter :: ccn_name(psat)=(/'CCN1','CCN2','CCN3','CCN4','CCN5','CCN6'/)

  integer :: stat

! real, external :: erf_alt, erfc_num_recipes


! if ((its <= 146) .and. (146 <= ite) .and. (jts <= 210) .and. (210 <= jte)) then
!          write(6, *) 'in Cmixactivate at beiginning Here22============ '		   
!          write(6, *)'qndrop3d(146,22,210)====  ', qndrop3d(146,22,210)
! end if

  colmass_worst(:,:,:) = 0.0
  colmass_worst_ij(:,:,:,:) = -1

!  cfrac_num (:,:,:,:) =0.5
!  cfrac_mass (:,:,:,:) =0.5
  
  arg = 1.0
  if (abs(0.8427-ERF_ALT(arg))/0.8427>0.001) then
     write (6,*) 'erf_alt(1.0) = ',ERF_ALT(arg)
     call wrf_error_fatal('dropmixnuc: Error function error')
  endif
  arg = 0.0
  if (ERF_ALT(arg) /= 0.0) then
     write (6,*) 'erf_alt(0.0) = ',ERF_ALT(arg)
     call wrf_error_fatal('dropmixnuc: Error function error')
  endif

  pi = 4.*atan(1.0)
  dtinv=1./dtstep

  depvel_drop =  0.1 ! prescribed here rather than getting it from dry_dep_driver
  if (idrydep_onoff .le. 0) depvel_drop =  0.0
  depvel_drop =  min(depvel_drop,depvel_uplimit)

!---------this has to be done within K-loop because sigmag_aer(m,n) changes for CMAQ case 
!   except ncomp(n) here (Shaocai Yu)
!  do n=1,ntype_aer
!     do m=1,nsize_aer(n)
!        ncomp(n)=ncomp_aer(n)
!        alogsig(m,n)=alog(sigmag_aer(m,n))
!        dgnum_aer(m,n) = dpvolmean_aer(m,n) * exp( -1.5*alogsig(m,n)*alogsig(m,n) )
!    print *,'sigmag_aer,dgnum_aer=',sigmag_aer(m,n),dgnum_aer(m,n)
!        ! npv is used only if number is diagnosed from volume
!        npv(m,n)=6./(pi*(0.01*dgnum_aer(m,n))**3*exp(4.5*alogsig(m,n)*alogsig(m,n)))
!     end do
!  end do
!---------------------------------------------------------------------------------

  do n=1,ntype_aer
        ncomp(n)=ncomp_aer(n)
  end do

  t0=273.15   !wig, 1-Mar-2009: Added .15
  aten=2.*surften/(r_v*t0*rhowater)
  super(:)=0.01*supersat(:)

!---------this has to be done within K-loop because sigmag_aer(m,n) changes for CMAQ case
!   (Shaocai Yu)  
!  do n=1,ntype_aer
!     do m=1,nsize_aer(n)
!        exp45logsig=exp(4.5*alogsig(m,n)*alogsig(m,n))
!        argfactor(m,n)=2./(3.*sqrt(2.)*alogsig(m,n))
!        amcube(m,n)=3./(4.*pi*exp45logsig*npv(m,n))
!     enddo
!  enddo
!----------------------------------------------------------------------

!------Shaocai	
! if ((its <= 146) .and. (146 <= ite) .and. (jts <= 210) .and. (210 <= jte)) then
!        write(6,*)'in mixactivate before call cal_cldfra2  '
! write(6,*)'qndrop3d(146,22,210)',qndrop3d(146,22,210)  
! write(6,*)'cldfra(146,22,210)',cldfra(146,22,210)
! end if

!-----------------------------------------------------------


  IF( PRESENT(F_QC) .AND. PRESENT ( F_QI ) ) THEN
     CALL cal_cldfra2(CLDFRA,qc,qi,f_qc,f_qi,     &   
          ids,ide, jds,jde, kds,kde,              &
          ims,ime, jms,jme, kms,kme,              &
          its,ite, jts,jte, kts,kte               )
!------Shaocai	
! if ((its <= 146) .and. (146 <= ite) .and. (jts <= 210) .and. (210 <= jte)) then
!          write(6,*)'in mixactivate just after call cal_cldfra  '
! write(6,*)'qndrop3d(146,22,210)',qndrop3d(146,22,210)  
! write(6,*)'cldfra(146,22,210)',cldfra(146,22,210)
! end if

!-----------------------------------------------------------	

  END IF

!------Shaocai	
! if ((its <= 146) .and. (146 <= ite) .and. (jts <= 210) .and. (210 <= jte)) then
!        write(6,*)'in mixactivate after call cal_cldfra  '
!        write(6,*)'qndrop3d(146,22,210)',qndrop3d(146,22,210)  
!        write(6,*)'cldfra(146,22,210)',cldfra(146,22,210)
!        write(6,*)'before loop,jts,jte,its,ite,kts,kte',jts,jte,its,ite,kts,kte 
! end if

!-----------------------------------------------------------
!  for test only (Shaocai)
     do j=jts,jte
     do i=its,ite
     do k=kts,kte
     do n=1,ntype_aer
        do m=1,nsize_aer(n)
           lnum=numptr_aer(m,n,ai_phase)
           lnumcw=numptr_aer(m,n,cw_phase)
!             if((chem(i,k,j,lnum)+chem(i,k,j,lnumcw)) == 0.0) then
!	        write(6,*)' heree1==d== 00 chem(i,k,j,lnum)',chem(i,k,j,lnumcw),chem(i,k,j,lnum),i,k,j,lnum
!             endif
           do l=1,ncomp(n)	      
              lmass=massptr_aer(l,m,n,ai_phase)
              lmasscw=massptr_aer(l,m,n,cw_phase)
              
!             if((chem(i,k,j,lmass)+chem(i,k,j,lmasscw)) == 0.0) then
!                write(6,*)' Heree2==d== 00 chem(i,k,j,lmass)',chem(i,k,j,lmass),chem(i,k,j,lmasscw),i,k,j,lmass, lmasscw
!             endif
      enddo
      enddo
      enddo
      enddo
      enddo      
      enddo      
!--------------------------------------------------------------------


  qsrflx(its:ite,jts:jte,:) = 0.

!     start loop over columns

OVERALL_MAIN_J_LOOP: do j=jts,jte
OVERALL_MAIN_I_LOOP: do i=its,ite

!        load number nucleated into qndrop on cloud boundaries

! initialization for current i .........................................

     do k=kts+1,kte
	    zs(k)=1./(zm(i,k,j)-zm(i,k-1,j))
	 enddo
	 zs(kts)=zs(kts+1)
! David: removed the next line
    zs(kte+1)=0.

     do k=kts,kte
!!$	    if(qndrop3d(i,k,j).lt.-10.e6.or.qndrop3d(i,k,j).gt.1.E20)then
!!$!	       call wrf_error_fatal("1")
!!$	    endif
        if(f_qi)then
           qcld=qc(i,k,j)+qi(i,k,j)
        else
           qcld=qc(i,k,j)
        endif
        if(qcld.lt.-1..or.qcld.gt.1.)then
           write(6,'(a,g12.2,a,3i5)')'qcld=',qcld,' for i,k,j=',i,k,j
           call wrf_error_fatal("1")
        endif
        if(qcld.gt.1.e-20)then
           lcldfra(k)=cldfra(i,k,j)*qc(i,k,j)/qcld
           lcldfra_old(k)=cldfra_old(i,k,j)*qc(i,k,j)/qcld
	   
!------Shaocai
!        if(i.eq.146.and.k.eq.22.and.j.eq.210)then
!        write(6,*)'in mixactivate after if(qcld.gt.1.e-20)then'
! write(6,*)'qndrop3d(i,k,j), i,k,j',qndrop3d(i,k,j),i,k,j 
! write(6,*)'cldfra(i,k,j),qc(i,k,j),qcld',cldfra(i,k,j),qc(i,k,j),qcld
! write(6,*)'cldfra_old(i,k,j),qc(i,k,j),qcld',cldfra_old(i,k,j),qc(i,k,j),qcld
! write(6,*)'lcldfra(k),lcldfra_old(k)',lcldfra(k),lcldfra_old(k)
! endif
!---------------------------------------		   
	   
	   
        else
           lcldfra(k)=0.
           lcldfra_old(k)=0.
        endif
        qndrop(k)=qndrop3d(i,k,j)
	
!------Shaocai
!        if(i.eq.146.and.k.eq.22.and.j.eq.210)then	 
!        write(6,*)'in mixactivate after qndrop(k)=qndrop3d(i,k,j)then'
! write(6,*)'qndrop3d(i,k,j),qndrop(k), i,k,j',qndrop3d(i,k,j),qndrop(k),i,k,j 
! write(6,*)'cldfra(i,k,j),qc(i,k,j),qcld',cldfra(i,k,j),qc(i,k,j),qcld
! write(6,*)'cldfra_old(i,k,j),qc(i,k,j),qcld',cldfra_old(i,k,j),qc(i,k,j),qcld
! write(6,*)'lcldfra(k),lcldfra_old(k)',lcldfra(k),lcldfra_old(k)
! endif
!---------------------------------------		
	
!	    qndrop(k)=1.e5
        cs(k)=rho(i,k,j) ! air density (kg/m3)
        dz(k)=dz8w(i,k,j)
!---------------------------------Shaocai
!if(dz(k).eq.0.0)then
!  write(6,*)'dz8w(i,k,j)',dz8w(i,k,j),i,k,j
!endif
!--------------------------------------------------
        do n=1,ntype_aer
           do m=1,nsize_aer(n)
              nact(k,m,n)=0.
              mact(k,m,n)=0.
           enddo
        enddo
        zn(k)=1./(cs(k)*dz(k))
        if(k>kts)then
           ekd(k)=kvh(i,k,j)
           ekd(k)=max(ekd(k),zkmin)
           ekd(k)=min(ekd(k),zkmax)
        else
           ekd(k)=0
        endif
!           diagnose subgrid vertical velocity from diffusivity
        if(k.eq.kts)then
           wtke(k)=sq2pi*depvel_drop
!               wtke(k)=sq2pi*kvh(i,k,j)
!               wtke(k)=max(wtke(k),wmixmin)
        else
           wtke(k)=sq2pi*ekd(k)/dz(k)
        endif
        wtke(k)=max(wtke(k),wmixmin)
        nsource(i,k,j)=0.
     enddo
! David: removed the next three lines
    nsource(i,kte+1,j) = 0.
    qndrop(kte+1)      = 0.
    zn(kte+1)          = 0.

     do k = kts+1, kte
        tmpa = dz(k-1) ; tmpb = dz(k)
        tmpc = tmpa/(tmpa + tmpb)
        csbot(k) = cs(k-1)*(1.0-tmpc) + cs(k)*tmpc
        csbot_cscen(k) = csbot(k)/cs(k)
     end do
     csbot(kts) = cs(kts)
     csbot_cscen(kts) = 1.0
! David: removed the next three lines
    csbot(kte+1) = cs(kte)
    csbot_cscen(kte+1) = 1.0

    !  calculate surface rate and mass mixing ratio for aerosol

     surfratemax = 0.0
     nsav=1
     nnew=2
     surfrate_drop=depvel_drop/dz(kts)
     surfratemax = max( surfratemax, surfrate_drop )
     do n=1,ntype_aer
        do m=1,nsize_aer(n)
           lnum=numptr_aer(m,n,ai_phase)
           lnumcw=numptr_aer(m,n,cw_phase)
           if(lnum>0)then
              depvel_tmp = max( 0.0, min( ddvel(i,j,lnum), depvel_uplimit ) )
              surfrate(lnum)=depvel_tmp/dz(kts)
              surfrate(lnumcw)=surfrate_drop
              surfratemax = max( surfratemax, surfrate(lnum) )
!             scale = 1000./mwdry ! moles/kg
              scale = 1.
              raercol(kts:kte,lnumcw,nsav)=chem(i,kts:kte,j,lnumcw)*scale ! #/kg
              raercol(kts:kte,lnum,nsav)=chem(i,kts:kte,j,lnum)*scale

!------Shaocai
!         if(i.eq.177.and.j.eq.206)then	 
!	 write(6,*)'if(lnum>0)then'
!	 write(6,*)'raercol(13,lnumcw,nsav)',raercol(13,lnumcw,nsav)
!         write(6,*)'raercol(13,lnum,nsav)',raercol(13,lnum,nsav)
!	 write(6,*)'kts:kte,lnum,lnumcw,nsav',kts,kte,lnum,lnumcw,nsav
!	 endif
!-----------------------------------------------	      
	      
	      
           endif
           do l=1,ncomp(n)
              lmass=massptr_aer(l,m,n,ai_phase)
              lmasscw=massptr_aer(l,m,n,cw_phase)
!             scale = mw_aer(l,n)/mwdry
              scale = 1.e-9 ! kg/ug
              depvel_tmp = max( 0.0, min( ddvel(i,j,lmass), depvel_uplimit ) )
              surfrate(lmass)=depvel_tmp/dz(kts)
              surfrate(lmasscw)=surfrate_drop
              surfratemax = max( surfratemax, surfrate(lmass) )
              raercol(kts:kte,lmasscw,nsav)=chem(i,kts:kte,j,lmasscw)*scale ! kg/kg
              raercol(kts:kte,lmass,nsav)=chem(i,kts:kte,j,lmass)*scale ! kg/kg

!--------------Shaocai
	      do k=kts,kte
		if((raercol(k,lmasscw,nsav)+raercol(k,lmass,nsav)).eq.0)then
		    raercol(k,lmasscw,nsav) = 1.e-30
		    raercol(k,lmass,nsav)   = 1.e-30
!-------------------------------------------------------------------------------
!                  1.0e-9*1.0e-30=1.0e-39===0 in here because of very small number (Shaocai)
!		  write(6,*)'000Here1,raercol(k,lmasscw,nsav)',raercol(k,lmasscw,nsav)+raercol(k,lmass,nsav),&
!		     k,lmasscw,lmass,nsav, kts,kte		    
!                  write(6,*)'chem(i,k,j,lmasscw)',chem(i,k,j,lmasscw),chem(i,k,j,lmass),i,k,j
!----------------------------------------------------------
		  
		endif		     
              enddo
!------------------------------------------------------------

!------Shaocai
!        if(i.eq.146.and.j.eq.210.and.l.eq.1)then
! write(6,*)'after do l=1,ncomp(n)',1,ncomp(n)	
! write(6,*)'raercol(13,lmasscw,nsav)',raercol(13,lmasscw,nsav)
!        write(6,*)'raercol(13,lmass,nsav)',raercol(13,lmass,nsav)
! write(6,*)'kts:kte,lmass,lmasscw,nsav',kts,kte,lmass,lmasscw,nsav
! endif
!-----------------------------------------------

           enddo

           lwater=waterptr_aer(m,n)
           if(lwater>0)then
              depvel_tmp = max( 0.0, min( ddvel(i,j,lwater), depvel_uplimit ) )
              surfrate(lwater)=depvel_tmp/dz(kts)
              surfratemax = max( surfratemax, surfrate(lwater) )
              raercol(kts:kte,lwater,nsav)=chem(i,kts:kte,j,lwater) ! don't bother to convert units,
             ! because it doesn't contribute to aerosol mass
           endif
        enddo ! size
     enddo ! type


! mass conservation checking
     if (icheck_colmass > 0) then

! calc initial column burdens
        colmass_bgn(:,:,:,:) = 0.0
        colmass_end(:,:,:,:) = 0.0
        colmass_sfc(:,:,:,:) = 0.0
        rhodz(kts:kte) = 1.0/zn(kts:kte)
        rhodzsum = sum( rhodz(kts:kte) )
        do n=1,ntype_aer
           do m=1,nsize_aer(n)
              lnum=numptr_aer(m,n,ai_phase)
              lnumcw=numptr_aer(m,n,cw_phase)
              if(lnum>0)then
                 colmass_bgn(0,m,n,1) = sum( chem(i,kts:kte,j,lnum  )*rhodz(kts:kte) )
                 colmass_bgn(0,m,n,2) = sum( chem(i,kts:kte,j,lnumcw)*rhodz(kts:kte) )
              endif
              do l=1,ncomp(n)
                 lmass=massptr_aer(l,m,n,ai_phase)
                 lmasscw=massptr_aer(l,m,n,cw_phase)
                 colmass_bgn(l,m,n,1) = sum( chem(i,kts:kte,j,lmass  )*rhodz(kts:kte) )
                 colmass_bgn(l,m,n,2) = sum( chem(i,kts:kte,j,lmasscw)*rhodz(kts:kte) )
              enddo
           enddo ! size
        enddo ! type
     endif ! (icheck_colmass > 0)


!        droplet nucleation/aerosol activation

! k-loop for growing/shrinking cloud calcs .............................
GROW_SHRINK_MAIN_K_LOOP: do k=kts,kte
        km1=max0(k-1,1)
        kp1=min0(k+1,kde-1)


!       if(lcldfra(k)-lcldfra_old(k).gt.0.01)then   ! this line is the "old" criterion
!       go to 10

!       growing cloud PLUS
!       upwards vertical advection when lcldfra(k-1) < lcldfra(k)
!
! tmpc1 = cloud fraction increase from previous time step
        tmpc1 = max( (lcldfra(k)-lcldfra_old(k)), 0.0 )
        if (k > kts) then
! tmpc2 = fraction of layer for which vertical advection from below
!             (over dtstep) displaces cloudy air with clear air
!       = (courant number using upwards w at layer bottom)*(difference in cloud fraction)
           tmpcourno = dtstep*max(w(i,k,j),0.0)/dz(k)
           tmpc2 = max( (lcldfra(k)-lcldfra(km1)), 0.0 ) * tmpcourno
           tmpc2 = min( tmpc2, 1.0 )
!          tmpc2 = 0.0   ! this turns off the vertical advect part
        else
           tmpc2 = 0.0
        endif

        if ((tmpc1 > 0.001) .or. (tmpc2 > 0.001)) then
	
!------Shaocai
!        if(i.eq.146.and.k.eq.22.and.j.eq.210)then	 
! write(6,*)'after if ((tmpc1 > 0.001) .or. (tmpc2 > 0.001))'
! write(6,*)'cldfra(i,k,j),cldfra_old(i,k,j)',cldfra(i,k,j),cldfra_old(i,k,j)
!        write(6,*)'tmpc1,tmpc2,0.001',tmpc1,tmpc2,0.001
! write(6,*)'lcldfra(k),lcldfra_old(k)',lcldfra(k),lcldfra_old(k)
! endif
!---------------------------------------

!                wmix=wtke(k)
           wbar=w(i,k,j)+wtke(k)     !wtke(k) =infinite
!	   write(6,*)'wbar,w(i,k,j),wtke(k),i,k,j',wbar,w(i,k,j),wtke(k),i,k,j
           wmix=0.
           wmin=0.
! 06-nov-2005 rce - increase wmax from 10 to 50 (deep convective clouds)
           wmax=50.
           wdiab=0

!---------------------------------------------------------------------------------
!       Put CMAQ aerosol for each cell to dlo_sect etc. (Shaocai Yu)
!----------------------------------------------------------------------------------       
         do n=1,ntype_aer
         do m=1,nsize_aer(n)
            dlo_sect(m,n)     =cdlo_sect(i,k,j,m,n)
	    dhi_sect(m,n)     =cdhi_sect(i,k,j,m,n)   
	    sigmag_aer(m,n)   =csigmag_aer(i,k,j,m,n)
	    dpvolmean_aer(m,n)=cdpvolmean_aer(i,k,j,m,n) 
	    
            alogsig(m,n)=alog(sigmag_aer(m,n))
            dgnum_aer(m,n) = dpvolmean_aer(m,n) * exp( -1.5*alogsig(m,n)*alogsig(m,n) )
	    
        if((i.eq.146.and.k.eq.22.and.j.eq.210))then	    	    
           print *,'sigmag_aer,dgnum_aer=',sigmag_aer(m,n),dgnum_aer(m,n),m,n,i,k,j
           print *,'ntype_aer,nsize_aer(n)=',ntype_aer,nsize_aer(n),n
	   print *,'ncomp(n)',ncomp(n)
	   print *,'alogsig(m,n)',alogsig(m,n),m,n
	   print *,'dgnum_aer(m,n)',dgnum_aer(m,n)
	endif
	
! npv is used only if number is diagnosed from volume
            npv(m,n)=6./(pi*(0.01*dgnum_aer(m,n))**3*exp(4.5*alogsig(m,n)*alogsig(m,n)))
         end do
         end do	
	    
  do n=1,ntype_aer
     do m=1,nsize_aer(n)
        exp45logsig=exp(4.5*alogsig(m,n)*alogsig(m,n))
        argfactor(m,n)=2./(3.*sqrt(2.)*alogsig(m,n))
        amcube(m,n)=3./(4.*pi*exp45logsig*npv(m,n))
        if((i.eq.146.and.k.eq.22.and.j.eq.210))then	    	    
          print *,'exp45logsig',exp45logsig
          print *,'argfactor(m,n)',argfactor(m,n),m,n
          print *,'amcube(m,n)',amcube(m,n),m,n	  	   
        endif	
	
     enddo
  enddo	  
!----------------------------------------------------------------------------------------

!                load aerosol properties, assuming external mixtures
           do n=1,ntype_aer
              do m=1,nsize_aer(n)
!------Shaocai
!       if((i.eq.146.and.k.eq.22.and.j.eq.210))then
!        write(6,*)'======in mixactivate before calling loadaer'
!        write(6,*)'n,m=',n,m
!        write(6,*)'raercol(1,1,nsav),k,kms,kme',raercol(1,1,nsav),k,kms,kme
!        write(6,*)'num_chem,cs(k), npv(m,n)',num_chem,cs(k), npv(m,n)    
!        write(6,*)'dlo_sect(m,n),dhi_sect(m,n)',dlo_sect(m,n),dhi_sect(m,n)             
!        write(6,*)'maxd_acomp, ncomp(n)',maxd_acomp, ncomp(n) 
!        write(6,*)'grid_id, ktau, i, j, m, n',grid_id, ktau, i, j, m, n   
!        write(6,*)'numptr_aer(m,n,ai_phase),numptr_aer(m,n,cw_phase)',&
!        numptr_aer(m,n,ai_phase),numptr_aer(m,n,cw_phase)  
!        write(6,*)'dens_aer(1,n)', dens_aer(1,n)   
!        write(6,*)'massptr_aer(1,m,n,ai_phase), massptr_aer(1,m,n,cw_phase)', &
!        massptr_aer(1,m,n,ai_phase), massptr_aer(1,m,n,cw_phase)  
!        write(6,*)'maerosol(1,m,n), maerosolcw(1,m,n)', &
!        maerosol(1,m,n), maerosolcw(1,m,n)          
!        write(6,*)'maerosol_tot(m,n), maerosol_totcw(m,n)', &
!        maerosol_tot(m,n), maerosol_totcw(m,n)      
!        write(6,*)'naerosol(m,n), naerosolcw(m,n)',naerosol(m,n), naerosolcw(m,n)                  
!        write(6,*)'vaerosol(m,n), vaerosolcw(m,n)',vaerosol(m,n), vaerosolcw(m,n)
!        write(6,*)'hygro_aer(m,n)',hygro_aer(m,n)	 
!        endif
!---------------------------------------	     	      
	      
                 call loadaer(raercol(1,1,nsav),k,kms,kme,num_chem,    &
                      cs(k), npv(m,n), dlo_sect(m,n),dhi_sect(m,n),             &
                      maxd_acomp, ncomp(n), &
                      grid_id, ktau, i, j, m, n,   &
                      numptr_aer(m,n,ai_phase),numptr_aer(m,n,cw_phase),  &
                      dens_aer(1,n),    &
                      massptr_aer(1,m,n,ai_phase), massptr_aer(1,m,n,cw_phase),  &
                      maerosol(1,m,n), maerosolcw(1,m,n),          &
                      maerosol_tot(m,n), maerosol_totcw(m,n),      &
                      naerosol(m,n), naerosolcw(m,n),                  &
                      vaerosol(m,n), vaerosolcw(m,n) )

                 hygro_aer(m,n)=hygro(i,k,j,m,n)

!------Shaocai
!               if(hygro(i,k,j,m,n).eq.0.0)then
!	  write(6,*)'here1 hygro_aer(m,n)=',hygro_aer(m,n),hygro(i,k,j,m,n),i,k,j,m,n
!	endif
!       if((i.eq.146.and.k.eq.22.and.j.eq.210))then		
!        write(6,*)'======in mixactivate after calling loadaer'
!        write(6,*)'n,m=',n,m
!        write(6,*)'raercol(1,1,nsav),k,kms,kme',raercol(1,1,nsav),k,kms,kme
!        write(6,*)'num_chem,cs(k), npv(m,n)',num_chem,cs(k), npv(m,n)    
!        write(6,*)'dlo_sect(m,n),dhi_sect(m,n)',dlo_sect(m,n),dhi_sect(m,n)             
!        write(6,*)'maxd_acomp, ncomp(n)',maxd_acomp, ncomp(n) 
!        write(6,*)'grid_id, ktau, i, j, m, n',grid_id, ktau, i, j, m, n   
!        write(6,*)'numptr_aer(m,n,ai_phase),numptr_aer(m,n,cw_phase)',&
!        numptr_aer(m,n,ai_phase),numptr_aer(m,n,cw_phase)  
!        write(6,*)'dens_aer(1,n)', dens_aer(1,n)   
!        write(6,*)'massptr_aer(1,m,n,ai_phase), massptr_aer(1,m,n,cw_phase)', &
!        massptr_aer(1,m,n,ai_phase), massptr_aer(1,m,n,cw_phase)  
!        write(6,*)'maerosol(1,m,n), maerosolcw(1,m,n)', &
!        maerosol(1,m,n), maerosolcw(1,m,n)          
!        write(6,*)'maerosol_tot(m,n), maerosol_totcw(m,n)', &
!        maerosol_tot(m,n), maerosol_totcw(m,n)      
!        write(6,*)'naerosol(m,n), naerosolcw(m,n)',naerosol(m,n), naerosolcw(m,n)                  
!        write(6,*)'vaerosol(m,n), vaerosolcw(m,n)',vaerosol(m,n), vaerosolcw(m,n)
!        write(6,*)'hygro_aer(m,n)',hygro_aer(m,n)
!        endif
!---------------------------------------			 
		 
		 
		 
		 
              enddo
           enddo

! 06-nov-2005 rce - grid_id & ktau added to arg list

!------Shaocai
!       if((i.eq.146.and.k.eq.22.and.j.eq.210))then
!         if(qndrop(k).gt. 0.0) then
!        write(6,*)'in mixactivate before calling activate-firsttime'
! write(6,*)'qndrop(k), i,k,j',qndrop(k),i,k,j 
! write(6,*)'lcldfra(k),lcldfra_old(k)',lcldfra(k),lcldfra_old(k)
!	 endif
! endif
!---------------------------------------

!------Shaocai
!       if((i.eq.146.and.k.eq.22.and.j.eq.210))then
!        write(6,*)'======in mixactivate before calling activate-firsttime'
! write(6,*)'wbar,wmix,wdiab,wmin,wmax',wbar,wmix,wdiab,wmin,wmax
! write(6,*)'temp(i,k,j),cs(k)',temp(i,k,j),cs(k)
! write(6,*)'msectional, maxd_atype=',msectional, maxd_atype 
! write(6,*)'ntype_aer, maxd_asize, nsize_aer=',ntype_aer, maxd_asize, nsize_aer
! write(6,*)'naerosol, vaerosol',naerosol, vaerosol,i,k,j
! write(6,*)'raercol(1,1,nsav)',raercol(1,1,nsav)
! write(6,*)'dlo_sect,dhi_sect,sigmag_aer=',dlo_sect,dhi_sect,sigmag_aer	 
! write(6,*)'hygro_aer(m,n)',hygro_aer	  	 
! write(6,*)'fn(1,1),fs(1,1),fm(1,1)=',fn,fs,fm          
!        write(6,*)'fluxn,fluxs,fluxm',fluxn,fluxs,fluxm
! endif
!---------------------------------------

           call activate(wbar,wmix,wdiab,wmin,wmax,temp(i,k,j),cs(k), &
                msectional, maxd_atype, ntype_aer, maxd_asize, nsize_aer,    &
                naerosol, vaerosol,  &
                dlo_sect,dhi_sect,sigmag_aer,hygro_aer,              &
                fn,fs,fm,fluxn,fluxs,fluxm,flux_fullact(k), grid_id, ktau, i, j, k )

!------Shaocai
!        do n = 1,ntype_aer
!        do m = 1,nsize_aer(n)
! if(naerosol(m,n).gt.0.0)then
!       if((i.eq.146.and.k.eq.22.and.j.eq.210))then
!        write(6,*)'======in mixactivate just after calling activate-firsttime'
! write(6,*)'wbar,wmix,wdiab,wmin,wmax',wbar,wmix,wdiab,wmin,wmax
! write(6,*)'temp(i,k,j),cs(k),i,k,j',temp(i,k,j),cs(k),i,k,j
! write(6,*)'msectional, maxd_atype=',msectional, maxd_atype 
! write(6,*)'ntype_aer, maxd_asize, nsize_aer=',ntype_aer, maxd_asize, nsize_aer
! write(6,*)'naerosol, vaerosol',naerosol, vaerosol,i,k,j
! write(6,*)'raercol(1,1,nsav)',raercol(1,1,nsav)
! write(6,*)'dlo_sect,dhi_sect,sigmag_aer=',dlo_sect,dhi_sect,sigmag_aer	
! write(6,*)'hygro_aer(m,n)',hygro_aer	  
! write(6,*)'fn(1,1),fs(1,1),fm(1,1)=',n,m,fn,fs,fm          
!        write(6,*)'fluxn,fluxs,fluxm',fluxn,fluxs,fluxm
! endif
! endif
! enddo
! enddo
!---------------------------------------


           do n = 1,ntype_aer
              do m = 1,nsize_aer(n)
                 lnum   = numptr_aer(m,n,ai_phase)
                 lnumcw = numptr_aer(m,n,cw_phase)
                 if (tmpc1 > 0.0) then
                    dact = tmpc1*fn(m,n)*raercol(k,lnum,nsav) ! interstitial only
                 else
                    dact = 0.0
                 endif
                 if (tmpc2 > 0.0) then
                    dact = dact + tmpc2*fn(m,n)*raercol(km1,lnum,nsav) ! interstitial only
                 endif
                 dact = min( dact, 0.99*raercol(k,lnum,nsav) )
                 raercol(k,lnumcw,nsav) = raercol(k,lnumcw,nsav)+dact
                 raercol(k,lnum,  nsav) = raercol(k,lnum,  nsav)-dact
                 qndrop(k) = qndrop(k)+dact
                 nsource(i,k,j) = nsource(i,k,j)+dact*dtinv
                 do l = 1,ncomp(n)
                    lmass   = massptr_aer(l,m,n,ai_phase)
                    lmasscw = massptr_aer(l,m,n,cw_phase)
                    if (tmpc1 > 0.0) then
                       dact = tmpc1*fm(m,n)*raercol(k,lmass,nsav) ! interstitial only
		       
!		       write(6,*)'dact',tmpc1,fm(m,n),raercol(k,lmass,nsav)
!		       stop
                    else
                       dact = 0.0
                    endif
                    if (tmpc2 > 0.0) then
                       dact = dact + tmpc2*fm(m,n)*raercol(km1,lmass,nsav) ! interstitial only
                    endif
                    dact = min( dact, 0.99*raercol(k,lmass,nsav) )
                    raercol(k,lmasscw,nsav)  =  raercol(k,lmasscw,nsav)+dact
                    raercol(k,lmass,  nsav)  =  raercol(k,lmass,  nsav)-dact
!--------------Shaocai
		if((raercol(k,lmasscw,nsav)+raercol(k,lmass,nsav)).eq.0)then
		    raercol(k,lmasscw,nsav) = 1.e-30
		    raercol(k,lmass,nsav)   = 1.e-30		
		
!-----------------------------------------------------		
!		  write(6,*)'000Here2,raercol(k,lmasscw,nsav)',raercol(k,lmasscw,nsav)+raercol(k,lmass,nsav),&
!		     k,lmasscw,lmass,nsav, kts,kte
!-----------------------------------------------------
		endif		     
!------------------------------------------------------------		    
		    
		    
                 enddo
              enddo
           enddo
!   10 continue

!------Shaocai	
!        if(i.eq.146.and.k.eq.22.and.j.eq.210)then	   
!         if(qndrop(k).gt. 0.0) then	
!        write(6,*)'in mixactivate after calling activate-firsttime'
! write(6,*)'qndrop(k), i,k,j',qndrop(k),i,k,j 
! write(6,*)'fn(1,1),fs(1,1),fm(1,1)=',fn,fs,fm 
! write(6,*)'naerosol, vaerosol',naerosol, vaerosol,i,k,j
! write(6,*)'raercol(1,1,nsav)',raercol(1,1,nsav)	          
! write(6,*)'lcldfra(k),lcldfra_old(k)',lcldfra(k),lcldfra_old(k)
!	 endif	
! endif  
!-----------------------------------------------



        endif   ! ((tmpc1 > 0.001) .or. (tmpc2 > 0.001))


        if(lcldfra(k) < lcldfra_old(k) .and. lcldfra_old(k) > 1.e-20)then   ! this line is the "old" criterion
!         go to 20

!       shrinking cloud ......................................................

!                droplet loss in decaying cloud
           nsource(i,k,j)=nsource(i,k,j)+qndrop(k)*(lcldfra(k)-lcldfra_old(k))*dtinv
           qndrop(k)=qndrop(k)*(1.+lcldfra(k)-lcldfra_old(k))
!                 convert activated aerosol to interstitial in decaying cloud

           tmpc = (lcldfra_old(k)-lcldfra(k))/lcldfra_old(k)
           do n=1,ntype_aer
              do m=1,nsize_aer(n)
                 lnum=numptr_aer(m,n,ai_phase)
                 lnumcw=numptr_aer(m,n,cw_phase)
                 if(lnum.gt.0)then
                    dact=raercol(k,lnumcw,nsav)*tmpc
                    raercol(k,lnumcw,nsav)=raercol(k,lnumcw,nsav)+dact
                    raercol(k,lnum,nsav)=raercol(k,lnum,nsav)-dact	    
                 endif
                 do l=1,ncomp(n)
                    lmass=massptr_aer(l,m,n,ai_phase)
                    lmasscw=massptr_aer(l,m,n,cw_phase)
                    dact=raercol(k,lmasscw,nsav)*tmpc
                    raercol(k,lmasscw,nsav)=raercol(k,lmasscw,nsav)+dact
                    raercol(k,lmass,nsav)=raercol(k,lmass,nsav)-dact
!--------------Shaocai
		if((raercol(k,lmasscw,nsav)+raercol(k,lmass,nsav)).eq.0)then
		    raercol(k,lmasscw,nsav) = 1.e-30
		    raercol(k,lmass,nsav)   = 1.e-30		
		
!		  write(6,*)'000Here3,raercol(k,lmasscw,nsav)',raercol(k,lmasscw,nsav)+raercol(k,lmass,nsav),&
!		     k,lmasscw,lmass,nsav, kts,kte
!---------------------------------------------
		endif		     
!------------------------------------------------------------		    
                 enddo
              enddo
           enddo
!             20 continue
        endif

     enddo GROW_SHRINK_MAIN_K_LOOP
! end of k-loop for growing/shrinking cloud calcs ......................



! ......................................................................
! start of main k-loop for calc of old cloud activation tendencies ..........
! this loop does "set up" for the nsubmix loop
!
! rce-comment
!    changed this part of code to use current cloud fraction (lcldfra) exclusively

OLD_CLOUD_MAIN_K_LOOP: do k=kts,kte
        km1=max0(k-1,kts)
        kp1=min0(k+1,kde-1)
        flux_fullact(k) = 0.0
        if(lcldfra(k).gt.0.01)then
!          go to 30

!               old cloud
              if(lcldfra(k)-lcldfra(km1).gt.0.01.or.k.eq.kts)then

!                   interior cloud
!                   cloud base

                 wdiab=0
                 wmix=wtke(k) ! spectrum of updrafts
                 wbar=w(i,k,j) ! spectrum of updrafts
!                    wmix=0. ! single updraft
!               wbar=wtke(k) ! single updraft
! 06-nov-2005 rce - increase wmax from 10 to 50 (deep convective clouds)
                 wmax=50.
                 ekd(k)=wtke(k)*dz(k)/sq2pi
                 alogarg=max(1.e-20,1/lcldfra(k)-1.)
                 wmin=wbar+wmix*0.25*sq2pi*alog(alogarg)
		 
!---------------------------------------------------------------------------------
!       Put CMAQ aerosol for each cell to dlo_sect etc. (Shaocai Yu)
!----------------------------------------------------------------------------------       
         do n=1,ntype_aer
         do m=1,nsize_aer(n)
            dlo_sect(m,n)     =cdlo_sect(i,k,j,m,n)
	    dhi_sect(m,n)     =cdhi_sect(i,k,j,m,n)   
	    sigmag_aer(m,n)   =csigmag_aer(i,k,j,m,n)
	    dpvolmean_aer(m,n)=cdpvolmean_aer(i,k,j,m,n) 
	    
            alogsig(m,n)=alog(sigmag_aer(m,n))
            dgnum_aer(m,n) = dpvolmean_aer(m,n) * exp( -1.5*alogsig(m,n)*alogsig(m,n) )
	    
        if((i.eq.146.and.k.eq.22.and.j.eq.210).and.(n.eq.1.and.m.eq.1))then	    
	    
           print *,'sigmag_aer,dgnum_aer=',sigmag_aer(m,n),dgnum_aer(m,n),m,n,i,k,j
	
	endif 
	 	 
! npv is used only if number is diagnosed from volume
            npv(m,n)=6./(pi*(0.01*dgnum_aer(m,n))**3*exp(4.5*alogsig(m,n)*alogsig(m,n)))
         end do
         end do	
	   
  do n=1,ntype_aer
     do m=1,nsize_aer(n)
        exp45logsig=exp(4.5*alogsig(m,n)*alogsig(m,n))
        argfactor(m,n)=2./(3.*sqrt(2.)*alogsig(m,n))
        amcube(m,n)=3./(4.*pi*exp45logsig*npv(m,n))
     enddo
  enddo	 	  
!----------------------------------------------------------------------------------------

                 do n=1,ntype_aer
                    do m=1,nsize_aer(n)
                       call loadaer(raercol(1,1,nsav),km1,kms,kme,num_chem,    &
                            cs(k), npv(m,n),dlo_sect(m,n),dhi_sect(m,n),               &
                            maxd_acomp, ncomp(n), &
                            grid_id, ktau, i, j, m, n,   &
                            numptr_aer(m,n,ai_phase),numptr_aer(m,n,cw_phase),  &
                            dens_aer(1,n),   &
                            massptr_aer(1,m,n,ai_phase), massptr_aer(1,m,n,cw_phase),  &
                            maerosol(1,m,n), maerosolcw(1,m,n),          &
                            maerosol_tot(m,n), maerosol_totcw(m,n),      &
                            naerosol(m,n), naerosolcw(m,n),                  &
                            vaerosol(m,n), vaerosolcw(m,n) )

                       hygro_aer(m,n)=hygro(i,k,j,m,n)

                    enddo
                 enddo
!          print *,'old cloud wbar,wmix=',wbar,wmix

                 call activate(wbar,wmix,wdiab,wmin,wmax,temp(i,k,j),cs(k), &
                      msectional, maxd_atype, ntype_aer, maxd_asize, nsize_aer,    &
                      naerosol, vaerosol,  &
                      dlo_sect,dhi_sect, sigmag_aer,hygro_aer,                    &
                      fn,fs,fm,fluxn,fluxs,fluxm,flux_fullact(k), grid_id, ktau, i, j, k )
                 
! rce-comment
!    the activation-fraction fluxes (fluxn, fluxm) from subr activate assume that
!       wbar << wmix, which is valid for global-model scale but not mesoscale
!    for wrf-chem application, divide these by flux_fullact to get a 
!       "flux-weighted-average" activation fraction, then multiply by (ekd(k)*zs(k)) 
!       which is the local "turbulent vertical-mixing velocity"
                 if (k > kts) then
                    if (flux_fullact(k) > 1.0e-20) then
                       tmpa = ekd(k)*zs(k)
                       tmpf = flux_fullact(k)
                       do n=1,ntype_aer
                       do m=1,nsize_aer(n)
                          tmpb = max( fluxn(m,n), 0.0 ) / max( fluxn(m,n), tmpf )
                          fluxn(m,n) = tmpa*tmpb
                          tmpb = max( fluxm(m,n), 0.0 ) / max( fluxm(m,n), tmpf )
                          fluxm(m,n) = tmpa*tmpb
                       enddo
                       enddo
                    else
                       fluxn(:,:) = 0.0
                       fluxm(:,:) = 0.0
                    endif
                 endif

                 if(k.gt.kts)then
                    tmpc = lcldfra(k)-lcldfra(km1)
                 else
                    tmpc=lcldfra(k)
                 endif
! rce-comment
!    flux of activated mass into layer k (in kg/m2/s)
!       = "actmassflux" = dumc*fluxm*raercol(kp1,lmass)*csbot(k)
!    source of activated mass (in kg/kg/s) = flux divergence
!       = actmassflux/(cs(i,k)*dz(i,k))
!    so need factor of csbot_cscen = csbot(k)/cs(i,k)
!                tmpe=1./(dz(k))
                 tmpe = csbot_cscen(k)/(dz(k))
                 fluxntot=0.
                 do n=1,ntype_aer
                 do m=1,nsize_aer(n)
                    fluxn(m,n)=fluxn(m,n)*tmpc
!                   fluxs(m,n)=fluxs(m,n)*tmpc
                    fluxm(m,n)=fluxm(m,n)*tmpc
                    lnum=numptr_aer(m,n,ai_phase)
                    fluxntot=fluxntot+fluxn(m,n)*raercol(km1,lnum,nsav)
!             print *,'fn=',fn(m,n),' for m,n=',m,n
!             print *,'old cloud tmpc=',tmpc,' fn=',fn(m,n),' for m,n=',m,n
                    nact(k,m,n)=nact(k,m,n)+fluxn(m,n)*tmpe
                    mact(k,m,n)=mact(k,m,n)+fluxm(m,n)*tmpe
                 enddo
                 enddo
                 flux_fullact(k) = flux_fullact(k)*tmpe
                 nsource(i,k,j)=nsource(i,k,j)+fluxntot*zs(k)
                 fluxntot=fluxntot*cs(k)
              endif
!       30 continue

        else
!       go to 40
!              no cloud
           if(qndrop(k).gt.10000.e6)then
              print *,'i,k,j,lcldfra,qndrop=',i,k,j,lcldfra(k),qndrop(k)
              print *,'cldfra,ql,qi',cldfra(i,k,j),qc(i,k,j),qi(i,k,j)
           endif
           nsource(i,k,j)=nsource(i,k,j)-qndrop(k)*dtinv
           qndrop(k)=0.
!              convert activated aerosol to interstitial in decaying cloud
           do n=1,ntype_aer
              do m=1,nsize_aer(n)
                 lnum=numptr_aer(m,n,ai_phase)
                 lnumcw=numptr_aer(m,n,cw_phase)
                 if(lnum.gt.0)then
                    raercol(k,lnum,nsav)=raercol(k,lnum,nsav)+raercol(k,lnumcw,nsav)
                    raercol(k,lnumcw,nsav)=0.
                 endif
                 do l=1,ncomp(n)
                    lmass=massptr_aer(l,m,n,ai_phase)
                    lmasscw=massptr_aer(l,m,n,cw_phase)
                    raercol(k,lmass,nsav)=raercol(k,lmass,nsav)+raercol(k,lmasscw,nsav)
                    raercol(k,lmasscw,nsav)=0.
!--------------Shaocai
		if((raercol(k,lmasscw,nsav)+raercol(k,lmass,nsav)).eq.0)then
		    raercol(k,lmasscw,nsav) = 1.e-30
		    raercol(k,lmass,nsav)   = 1.e-30		
		
!		  write(6,*)'000Here4,raercol(k,lmasscw,nsav)',raercol(k,lmasscw,nsav)+raercol(k,lmass,nsav),&
!		     k,lmasscw,lmass,nsav, kts,kte
		endif		     
!------------------------------------------------------------		    
		    
                 enddo
              enddo
           enddo
!      40 continue
        endif

     enddo OLD_CLOUD_MAIN_K_LOOP

!    cycle OVERALL_MAIN_I_LOOP
!------Shaocai	
!   if(i.eq.146.and.j.eq.210)then	   
!         if(qndrop(13).gt. 0.0) then	
!        write(6,*)'in mixactivate here2 after enddo OLD_CLOUD_MAIN_K_LOOP '
! write(6,*)'qndrop(13), i,k,j',qndrop(13),i,13,j  
!	 write(6,*)'fn(1,1),fs(1,1),fm(1,1)=',fn,fs,fm         
!	 write(6,*)'dlo_sect,dhi_sect,sigmag_aer=',dlo_sect,dhi_sect,sigmag_aer
!	 write(6,*)'msectional, maxd_atype=',msectional, maxd_atype 
!	 write(6,*)'ntype_aer, maxd_asize, nsize_aer=',ntype_aer, maxd_asize, nsize_aer 
! write(6,*)'lcldfra(13),lcldfra_old(13)',lcldfra(13),lcldfra_old(13)
!	 endif
! endif
!-----------------------------------------------------------

!    switch nsav, nnew so that nnew is the updated aerosol

     ntemp=nsav
     nsav=nnew
     nnew=ntemp

!    load new droplets in layers above, below clouds

     dtmin=dtstep
     ekk(kts)=0.0
! rce-comment -- ekd(k) is eddy-diffusivity at k/k-1 interface
!   want ekk(k) = ekd(k) * (density at k/k-1 interface)
     do k=kts+1,kte
        ekk(k)=ekd(k)*csbot(k)
     enddo
! David: remove the following line
    ekk(kte+1)=0.0
! David: reduce the ending index by 1
!     do k=kts,kte-1

    do k=kts,kte
        ekkp(k)=zn(k)*ekk(k+1)*zs(k+1)
        ekkm(k)=zn(k)*ekk(k)*zs(k)
        tinv=ekkp(k)+ekkm(k)
        if(k.eq.kts)tinv=tinv+surfratemax
        if(tinv.gt.1.e-6)then
           dtt=1./tinv
           dtmin=min(dtmin,dtt)
        endif
     enddo
     dtmix=0.9*dtmin
     nsubmix=dtstep/dtmix+1
     if(nsubmix>100)then
        nsubmix_bnd=100
     else
        nsubmix_bnd=nsubmix
     endif
!     count_submix(nsubmix_bnd)=count_submix(nsubmix_bnd)+1
     dtmix=dtstep/nsubmix
     fac_srflx = -1.0/(zn(1)*nsubmix)
     
     do k=kts,kte
        kp1=min(k+1,kde-1)
        km1=max(k-1,1)
        if(lcldfra(kp1).gt.0)then
           overlapp(k)=min(lcldfra(k)/lcldfra(kp1),1.)
        else
           overlapp(k)=1.
        endif
        if(lcldfra(km1).gt.0)then
           overlapm(k)=min(lcldfra(k)/lcldfra(km1),1.)
        else
           overlapm(k)=1.
        endif
     enddo

!------Shaocai	
!   if(i.eq.146.and.j.eq.210)then
!         if(qndrop(13).gt. 0.0) then	
!        write(6,*)'in mixactivate here1 '
! write(6,*)'qndrop(13), i,k,j',qndrop(13),i,13,j  
!	 write(6,*)'fn(1,1),fs(1,1),fm(1,1)=',fn,fs,fm         
! write(6,*)'lcldfra(13),lcldfra_old(13)',lcldfra(13),lcldfra_old(13)
!	 endif
! endif
!-----------------------------------------------------------

! ......................................................................
! start of nsubmix-loop for calc of old cloud activation tendencies ....
OLD_CLOUD_NSUBMIX_LOOP: do nsub=1,nsubmix
        qndrop_new(kts:kte)=qndrop(kts:kte)
!           switch nsav, nnew so that nsav is the updated aerosol
        ntemp=nsav
        nsav=nnew
        nnew=ntemp
        srcn(:)=0.0
        do n=1,ntype_aer
           do m=1,nsize_aer(n)
              lnum=numptr_aer(m,n,ai_phase)
!              update droplet source
! rce-comment - activation source in layer k involves particles from k-1
!             srcn(kts  :kte)=srcn(kts  :kte)+nact(kts  :kte,m,n)*(raercol(kts:kte  ,lnum,nsav))
              srcn(kts+1:kte)=srcn(kts+1:kte)+nact(kts+1:kte,m,n)*(raercol(kts:kte-1,lnum,nsav))
! rce-comment - new formulation for k=kts should be implemented
              srcn(kts      )=srcn(kts      )+nact(kts      ,m,n)*(raercol(kts      ,lnum,nsav))
           enddo
        enddo
        call explmix(qndrop,srcn,ekkp,ekkm,overlapp,overlapm,   &
             qndrop_new,surfrate_drop,kms,kme,kts,kte,dtmix,.false.)
	     
!------Shaocai	
!   if(i.eq.146.and.j.eq.210)then
!         if(qndrop(13).gt. 0.0) then	
!        write(6,*)'in mixactivate after call explmix () '
! write(6,*)'qndrop(13), i,k,j',qndrop(13),i,13,j  
!	 write(6,*)'fn(1,1),fs(1,1),fm(1,1)=',fn,fs,fm         
!	 write(6,*)'dlo_sect,dhi_sect,sigmag_aer=',dlo_sect,dhi_sect,sigmag_aer
!	 write(6,*)'msectional, maxd_atype=',msectional, maxd_atype 
!	 write(6,*)'ntype_aer, maxd_asize, nsize_aer=',ntype_aer, maxd_asize, nsize_aer 
! write(6,*)'lcldfra(13),lcldfra_old(13)',lcldfra(13),lcldfra_old(13)
!	 endif
! endif
!-----------------------------------------------------------	     
	     
	     
	     
        do n=1,ntype_aer
           do m=1,nsize_aer(n)
              lnum=numptr_aer(m,n,ai_phase)
              lnumcw=numptr_aer(m,n,cw_phase)
              if(lnum>0)then
! rce-comment - activation source in layer k involves particles from k-1
!                source(kts  :kte)= nact(kts  :kte,m,n)*(raercol(kts:kte  ,lnum,nsav))
                 source(kts+1:kte)= nact(kts+1:kte,m,n)*(raercol(kts:kte-1,lnum,nsav))
! rce-comment - new formulation for k=kts should be implemented
                 source(kts      )= nact(kts      ,m,n)*(raercol(kts      ,lnum,nsav))
                 call explmix(raercol(1,lnumcw,nnew),source,ekkp,ekkm,overlapp,overlapm, &
                      raercol(1,lnumcw,nsav),surfrate(lnumcw),kms,kme,kts,kte,dtmix,&
                      .false.)
                 call explmix(raercol(1,lnum,nnew),source,ekkp,ekkm,overlapp,overlapm,  &
                      raercol(1,lnum,nsav),surfrate(lnum),kms,kme,kts,kte,dtmix, &
                      .true.,raercol(1,lnumcw,nsav))
                 qsrflx(i,j,lnum) = qsrflx(i,j,lnum) + fac_srflx*            &
                      raercol(kts,lnum,nsav)*surfrate(lnum)
                 qsrflx(i,j,lnumcw) = qsrflx(i,j,lnumcw) + fac_srflx*        &
                      raercol(kts,lnumcw,nsav)*surfrate(lnumcw)
                 if (icheck_colmass > 0) then
                    tmpf = dtmix*rhodz(kts)
                    colmass_sfc(0,m,n,1) = colmass_sfc(0,m,n,1) &
                          + raercol(kts,lnum  ,nsav)*surfrate(lnum  )*tmpf
                    colmass_sfc(0,m,n,2) = colmass_sfc(0,m,n,2) &
                          + raercol(kts,lnumcw,nsav)*surfrate(lnumcw)*tmpf
                 endif
              endif
              do l=1,ncomp(n)
                 lmass=massptr_aer(l,m,n,ai_phase)
                 lmasscw=massptr_aer(l,m,n,cw_phase)
! rce-comment - activation source in layer k involves particles from k-1
!                source(kts  :kte)= mact(kts  :kte,m,n)*(raercol(kts:kte  ,lmass,nsav))
                 source(kts+1:kte)= mact(kts+1:kte,m,n)*(raercol(kts:kte-1,lmass,nsav))
! rce-comment - new formulation for k=kts should be implemented
                 source(kts      )= mact(kts      ,m,n)*(raercol(kts      ,lmass,nsav))
                 call explmix(raercol(1,lmasscw,nnew),source,ekkp,ekkm,overlapp,overlapm, &
                      raercol(1,lmasscw,nsav),surfrate(lmasscw),kms,kme,kts,kte,dtmix,  &
                      .false.)
                 call explmix(raercol(1,lmass,nnew),source,ekkp,ekkm,overlapp,overlapm,  &
                      raercol(1,lmass,nsav),surfrate(lmass),kms,kme,kts,kte,dtmix,  &
                      .true.,raercol(1,lmasscw,nsav))
                 qsrflx(i,j,lmass) = qsrflx(i,j,lmass) + fac_srflx*          &
                      raercol(kts,lmass,nsav)*surfrate(lmass)
                 qsrflx(i,j,lmasscw) = qsrflx(i,j,lmasscw) + fac_srflx*      &
                      raercol(kts,lmasscw,nsav)*surfrate(lmasscw)
                 if (icheck_colmass > 0) then
                    ! colmass_sfc calculation
                    !    colmass_bgn/end = bgn/end column burden = sum.over.k.of{ rho(k)*dz(k)*chem(k,l) }
                    !    colmass_sfc = surface loss over dtstep
                    !       = sum.over.nsubmix.substeps{ depvel(l)*rho(kts)*chem(kts,l)*dtmix }
                    !    surfrate(l) = depvel(l)/dz(kts) so need to multiply by dz(kts)
                    !    for mass, raercol(k,l) = chem(k,l)*1.0e-9, so need to multiply by 1.0e9
                    tmpf = dtmix*rhodz(kts)*1.0e9
                    colmass_sfc(l,m,n,1) = colmass_sfc(l,m,n,1) &
                          + raercol(kts,lmass  ,nsav)*surfrate(lmass  )*tmpf
                    colmass_sfc(l,m,n,2) = colmass_sfc(l,m,n,2) &
				  + raercol(kts,lmasscw,nsav)*surfrate(lmasscw)*tmpf
			 endif
		      enddo
		      lwater=waterptr_aer(m,n)  ! aerosol water
		      if(lwater>0)then
			 source(:)=0.
			 call explmix(   raercol(1,lwater,nnew),source,ekkp,ekkm,overlapp,overlapm,   &
			      raercol(1,lwater,nsav),surfrate(lwater),kms,kme,kts,kte,dtmix,  &
			      .true.,source)
		      endif
		   enddo ! size
		enddo ! type

	     enddo OLD_CLOUD_NSUBMIX_LOOP

	!------Shaocai	
!	   if(i.eq.146.and.j.eq.210)then
	!         if(qndrop(13).gt. 0.0) then	
!	 write(6,*)'in mixactivate after enddo OLD_CLOUD_NSUBMIX_LOOP'
!	 write(6,*)'qndrop(13), i,k,j',qndrop(13),i,13,j  
	!	 write(6,*)'fn(1,1),fs(1,1),fm(1,1)=',fn,fs,fm         
	!	 write(6,*)'dlo_sect,dhi_sect,sigmag_aer=',dlo_sect,dhi_sect,sigmag_aer
	!	 write(6,*)'msectional, maxd_atype=',msectional, maxd_atype 
	!	 write(6,*)'ntype_aer, maxd_asize, nsize_aer=',ntype_aer, maxd_asize, nsize_aer 
!	 write(6,*)'lcldfra(13),lcldfra_old(13)',lcldfra(13),lcldfra_old(13)
	!	 endif
!	 endif
	!-----------------------------------------------------------

	!    cycle OVERALL_MAIN_I_LOOP

	!        evaporate particles again if no cloud

	     do k=kts,kte
		if(lcldfra(k).eq.0.)then

	!              no cloud

		   qndrop(k)=0.
	!              convert activated aerosol to interstitial in decaying cloud
		   do n=1,ntype_aer
		      do m=1,nsize_aer(n)
			 lnum=numptr_aer(m,n,ai_phase)
			 lnumcw=numptr_aer(m,n,cw_phase)
			 if(lnum.gt.0)then
			    raercol(k,lnum,nnew)=raercol(k,lnum,nnew)+raercol(k,lnumcw,nnew)
			    raercol(k,lnumcw,nnew)=0.
			 endif
			 do l=1,ncomp(n)
			    lmass=massptr_aer(l,m,n,ai_phase)
                    lmasscw=massptr_aer(l,m,n,cw_phase)
                    raercol(k,lmass,nnew)=raercol(k,lmass,nnew)+raercol(k,lmasscw,nnew)
                    raercol(k,lmasscw,nnew)=0.
!--------------Shaocai
		if((raercol(k,lmasscw,nnew)+raercol(k,lmass,nnew)).eq.0)then

		    raercol(k,lmasscw,nnew) = 1.e-30
		    raercol(k,lmass,nnew)   = 1.e-30		
!		  write(6,*)'000Here5,raercol(k,lmasscw,nnew)',raercol(k,lmasscw,nnew)+raercol(k,lmass,nnew),&
!		     k,lmasscw,lmass,nnew, kts,kte
!------------------------------------------------
		endif		     
!------------------------------------------------------------		    
		    
		    
                 enddo
              enddo
           enddo
        endif
     enddo

!------Shaocai	
!        if(i.eq.146.and.j.eq.210)then
!         if(qndrop(13).gt. 0.0) then	
!        write(6,*)'in mixactivate after evaporate '
!        write(6,*)'qndrop(13), i,k,j',qndrop(13),i,13,j  
!	 write(6,*)'fn(1,1),fs(1,1),fm(1,1)=',fn,fs,fm         
!        write(6,*)'lcldfra(13),lcldfra_old(13)',lcldfra(13),lcldfra_old(13)
!	 endif
!       endif
!-----------------------------------------------------------

!         cycle OVERALL_MAIN_I_LOOP

!        droplet number

     do k=kts,kte
!       if(lcldfra(k).gt.0.1)then
!           write(6,'(a,3i5,f12.1)')'i,j,k,qndrop=',i,j,k,qndrop(k)
!       endif
!       if(qndrop(k).lt.-10.e6.or.qndrop(k).gt.1.e12)then
!          write(6,'(a,g12.2,a,3i5)')'after qndrop=',qndrop(k),' for i,k,j=',i,k,j
!       endif

        qndrop3d(i,k,j) = max(qndrop(k),1.e-6)

!-------------------------Shaocai Yu
!         if((qndrop3d(i,k,j).ge.1.0).and.(qc(i,k,j).gt.0.0).and.(k.lt.4)) then
!	    write(6,*)'problem here in Cmixactivate====d'
!	    write(6,*)'qndrop3d(i,k,j),qc(i,k,j)',qndrop3d(i,k,j),qc(i,k,j),i,k,j
!	    
!	 
!	 endif
!----------------------------------------------	


        if(qndrop3d(i,k,j).lt.-10.e6.or.qndrop3d(i,k,j).gt.1.E20)then
           write(6,'(a,g12.2,a,3i5)')'after qndrop=',qndrop3d(i,k,j),' for i,k,j=',i,k,j
        endif
        if(qc(i,k,j).lt.-1..or.qc(i,k,j).gt.1.)then
           write(6,'(a,g12.2,a,3i5)')'qc=',qc(i,k,j),' for i,k,j=',i,k,j
           call wrf_error_fatal("1")
        endif
        if(qi(i,k,j).lt.-1..or.qi(i,k,j).gt.1.)then
           write(6,'(a,g12.2,a,3i5)')'qi=',qi(i,k,j),' for i,k,j=',i,k,j
           call wrf_error_fatal("1")
        endif
        if(qv(i,k,j).lt.-1..or.qv(i,k,j).gt.1.)then
           write(6,'(a,g12.2,a,3i5)')'qv=',qv(i,k,j),' for i,k,j=',i,k,j
           call wrf_error_fatal("1")
        endif
        cldfra_old(i,k,j) = cldfra(i,k,j)
!       if(k.gt.6.and.k.lt.11)cldfra_old(i,k,j)=1.

!------Shaocai	
!        if(i.eq.146.and.k.eq.22.and.j.eq.210)then	   
!         if(qndrop(k).gt. 0.0) then	
!        write(6,*)'in mixactivate hereee '
! write(6,*)'qndrop(k), i,k,j',qndrop(k),i,k,j  
!	 write(6,*)'fn(1,1),fs(1,1),fm(1,1)=',fn,fs,fm         
! write(6,*)'lcldfra(k),lcldfra_old(k)',lcldfra(k),lcldfra_old(k)
!	 endif
! endif
!-----------------------------------------------------------

     enddo

!    cycle OVERALL_MAIN_I_LOOP

!        update chem and convert back to mole/mole

!   print *, ' ==d== 00000 e ', size(numptr_aer,1), size(numptr_aer,2), size(numptr_aer,3), ntype_aer, nsize_aer, ai_phase, cw_phase

     ccn(:,:) = 0.
     do n=1,ntype_aer
        do m=1,nsize_aer(n)
           lnum=numptr_aer(m,n,ai_phase)
           lnumcw=numptr_aer(m,n,cw_phase)
           if(lnum.gt.0)then
              !          scale=mwdry*0.001
              scale = 1.
              chem(i,kts:kte,j,lnumcw)= raercol(kts:kte,lnumcw,nnew)*scale
              chem(i,kts:kte,j,lnum)= raercol(kts:kte,lnum,nnew)*scale
	      
!--------------------------Shaocai
	      do k=kts,kte	
		if((raercol(k,lnumcw,nnew)+raercol(k,lnum,nnew)).eq.0)then
		  cfrac_num(i,k,j,lnum)=0.5  
		else
	          cfrac_num(i,k,j,lnum)       &
		    =raercol(k,lnum,nnew)/(raercol(k,lnumcw,nnew)+raercol(k,lnum,nnew))		
		endif
	      enddo
!----------------------------------------------

           endif
           do l=1,ncomp(n)
              lmass=massptr_aer(l,m,n,ai_phase)
              lmasscw=massptr_aer(l,m,n,cw_phase)
!          scale = mwdry/mw_aer(l,n)
              scale = 1.e9
              chem(i,kts:kte,j,lmasscw)=raercol(kts:kte,lmasscw,nnew)*scale ! ug/kg
              chem(i,kts:kte,j,lmass)=raercol(kts:kte,lmass,nnew)*scale ! ug/kg
	      
!--------------------------Shaocai	      
	      do k=kts,kte
		if((raercol(k,lmasscw,nnew)+raercol(k,lmass,nnew)).eq.0)then
!----------------------------------------------------------------
!		  write(*,*)'cfrac_mass(i,k,j,lmass)',raercol(k,lmasscw,nnew)+raercol(k,lmass,nnew),&
!		     cfrac_mass(i,k,j,lmass),i,k,j,lmass,k, lnumcw,nnew,kts, kte
!----------------------------------------------------------------		     
		  cfrac_mass(i,k,j,lmass)= 0.5  
		else		
	          cfrac_mass(i,k,j,lmass)       &
		    =raercol(k,lmass,nnew)/(raercol(k,lmasscw,nnew)+raercol(k,lmass,nnew))		
		   
!                  stop
		endif
		
	      enddo	      
!--------------------------	
      
           enddo
           lwater=waterptr_aer(m,n)
           if(lwater>0)chem(i,kts:kte,j,lwater)=raercol(kts:kte,lwater,nnew) ! don't convert units
           do k=kts,kte
!---------------------------------------------------------------------------------
!       Put CMAQ aerosol for each cell to dlo_sect etc. (Shaocai Yu)
!----------------------------------------------------------------------------------       
            dlo_sect(m,n)     =cdlo_sect(i,k,j,m,n)
	    dhi_sect(m,n)     =cdhi_sect(i,k,j,m,n)   
	    sigmag_aer(m,n)   =csigmag_aer(i,k,j,m,n)
	    dpvolmean_aer(m,n)=cdpvolmean_aer(i,k,j,m,n) 
	    
            alogsig(m,n)=alog(sigmag_aer(m,n))
            dgnum_aer(m,n) = dpvolmean_aer(m,n) * exp( -1.5*alogsig(m,n)*alogsig(m,n) )
	    
        if((i.eq.146.and.k.eq.22.and.j.eq.210).and.(n.eq.1.and.m.eq.1))then	    
	    
           print *,'sigmag_aer,dgnum_aer=',sigmag_aer(m,n),dgnum_aer(m,n),m,n,i,k,j
	
	endif 
	 	 
! npv is used only if number is diagnosed from volume
            npv(m,n)=6./(pi*(0.01*dgnum_aer(m,n))**3*exp(4.5*alogsig(m,n)*alogsig(m,n)))
	   
        exp45logsig=exp(4.5*alogsig(m,n)*alogsig(m,n))
        argfactor(m,n)=2./(3.*sqrt(2.)*alogsig(m,n))
        amcube(m,n)=3./(4.*pi*exp45logsig*npv(m,n))
!----------------------------------------------------------------------------------------	   
	   
              sm=2.*aten*sqrt(aten/(27.*hygro(i,k,j,m,n)*amcube(m,n)))
              do l=1,psat
                 arg=argfactor(m,n)*log(sm/super(l))
                 if(arg<2)then
                    if(arg<-2)then
                       ccnfact(l,m,n)=1.e-6 ! convert from #/m3 to #/cm3
                    else
                       ccnfact(l,m,n)=1.e-6*0.5*ERFC_NUM_RECIPES(arg)
                    endif
                 else
                    ccnfact(l,m,n) = 0.
                 endif
!                 ccn concentration as diagnostic
!                 assume same hygroscopicity and ccnfact for cloud-phase and aerosol phase particles

!------Shaocai	
!         if(i.eq.146.and.k.eq.22.and.j.eq.210)then	   
!         write(6,*)'in mixactivate before CCN (k,l)============ i,k,j,m,n',i,k,j,m,n,ntype_aer,nsize_aer(n)
!	 write(6,*)'aten,hygro(i,k,j,m,n),amcube(m,n)',aten,hygro(i,k,j,m,n),amcube(m,n)
!	 write(6,*)'arg,argfactor(m,n),sm,super(l),ERFC_NUM_RECIPES(arg)',arg,argfactor(m,n),sm,super(l),ERFC_NUM_RECIPES(arg)
!	 write(6,*)'qndrop(k), i,k,j,lnum,nnew',qndrop(k),i,k,j,lnum,nnew  
!	 write(6,*)'k,l,ccn(k,l),raercol(k,lnum,nnew)',k,l,ccn(k,l),raercol(k,lnum,nnew)
!         write(6,*)'raercol(k,lnumcw,nnew),cs(k),ccnfact(l,m,n)', &
!	        raercol(k,lnumcw,nnew),cs(k),ccnfact(l,m,n)        
!	 write(6,*)'super(l),lcldfra(k),lcldfra_old(k)',super(l),lcldfra(k),lcldfra_old(k)
!	 endif
!-----------------------------------------------------------


                 ccn(k,l)=ccn(k,l)+(raercol(k,lnum,nnew)+raercol(k,lnumcw,nnew))*cs(k)*ccnfact(l,m,n)
!------Shaocai	
!         if(i.eq.146.and.k.eq.22.and.j.eq.210)then	   
!         write(6,*)'in mixactivate after CCN (k,l)============= i,k,j,m,n',i,k,j,m,n
!	 write(6,*)'qndrop(k), i,k,j,lnum,nnew',qndrop(k),i,k,j,lnum,nnew  
!	 write(6,*)'k,l,ccn(k,l),raercol(k,lnum,nnew)',k,l,ccn(k,l),raercol(k,lnum,nnew)
!         write(6,*)'raercol(k,lnumcw,nnew),cs(k),ccnfact(l,m,n)', &
!	        raercol(k,lnumcw,nnew),cs(k),ccnfact(l,m,n)        
!	 write(6,*)'super(l),lcldfra(k),lcldfra_old(k)',super(l),lcldfra(k),lcldfra_old(k)
!	 endif
!-----------------------------------------------------------


              enddo
           enddo
        enddo
     enddo
     do l=1,psat
        !wig, 22-Nov-2006: added vertical bounds to prevent out-of-bounds at top
        if(l.eq.1)ccn1(i,kts:kte,j)=ccn(:,l)
        if(l.eq.2)ccn2(i,kts:kte,j)=ccn(:,l)
        if(l.eq.3)ccn3(i,kts:kte,j)=ccn(:,l)
        if(l.eq.4)ccn4(i,kts:kte,j)=ccn(:,l)
        if(l.eq.5)ccn5(i,kts:kte,j)=ccn(:,l)
        if(l.eq.6)ccn6(i,kts:kte,j)=ccn(:,l)
!------Shaocai	
!        if(i.eq.146.and.j.eq.210)then	   
!          write(6,*)'in mixactivate after do l=1,psat '
!	  write(6,*)'l,ccn(2,l),i,j',l,ccn(2,l),i,j
!	endif
!-----------------------------------------------------------	
	
     end do

! mass conservation checking
     if (icheck_colmass > 0) then
! calc final column burdens
        do n=1,ntype_aer
        do m=1,nsize_aer(n)
           lnum=numptr_aer(m,n,ai_phase)
           lnumcw=numptr_aer(m,n,cw_phase)
           if(lnum>0)then
              colmass_end(0,m,n,1) = sum( chem(i,kts:kte,j,lnum  )*rhodz(kts:kte) )
              colmass_end(0,m,n,2) = sum( chem(i,kts:kte,j,lnumcw)*rhodz(kts:kte) )
           endif
           do l=1,ncomp(n)
              lmass=massptr_aer(l,m,n,ai_phase)
              lmasscw=massptr_aer(l,m,n,cw_phase)
              colmass_end(l,m,n,1) = sum( chem(i,kts:kte,j,lmass  )*rhodz(kts:kte) )
              colmass_end(l,m,n,2) = sum( chem(i,kts:kte,j,lmasscw)*rhodz(kts:kte) )
           enddo
        enddo ! size
        enddo ! type
! calc burden change errors for each interstitial/activated pair
        do n=1,ntype_aer
        do m=1,nsize_aer(n)
           do l=0,ncomp(n)
              ! tmpa & tmpb = beginning & ending column burden divided by rhodzsum,
              !             = beginning & ending column-mean mixing ratios
              ! tmpc = loss to surface divided by rhodzsum,
              tmpa = ( colmass_bgn(l,m,n,1) + colmass_bgn(l,m,n,2) )/rhodzsum
              tmpb = ( colmass_end(l,m,n,1) + colmass_end(l,m,n,2) )/rhodzsum
              tmpc = ( colmass_sfc(l,m,n,1) + colmass_sfc(l,m,n,2) )/rhodzsum

              ! tmpd = ((final burden) + (sfc loss)) - (initial burden)
              !      = burden change error
              tmpd = (tmpb + tmpc) - tmpa
              tmpe = max( tmpa, 1.0e-20 )

              ! tmpf = (burden change error) / (initial burden)
              if (abs(tmpd) < 1.0e5*tmpe) then
                 tmpf = tmpd/tmpe
              else if (tmpf < 0.0) then
                 tmpf = -1.0e5
              else
                 tmpf = 1.0e5
              end if
              if (abs(tmpf) > abs(colmass_worst(l,m,n))) then
                 colmass_worst(l,m,n) = tmpf
                 colmass_worst_ij(1,l,m,n) = i
                 colmass_worst_ij(2,l,m,n) = j
              endif
           enddo
        enddo ! size
        enddo ! type
     endif ! (icheck_colmass > 0)


     enddo OVERALL_MAIN_I_LOOP ! end of main loop over i
     enddo OVERALL_MAIN_J_LOOP ! end of main loop over j


! mass conservation checking
     if (icheck_colmass > 0) then
        if (icheck_colmass >= 100) write(*,'(a)') &
             'mixactivate colmass worst errors bgn - type, size, comp, err, i, j'
        colmass_maxworst_r = 0.0
        colmass_maxworst_i(:) = -1
        do n=1,ntype_aer
        do m=1,nsize_aer(n)
           do l=0,ncomp(n)
              if (icheck_colmass >= 100) &
                 write(*,'(3i3,1p,e10.2,2i4)') n, m, l, &
                 colmass_worst(l,m,n), colmass_worst_ij(1:2,l,m,n) 
              if (abs(colmass_worst(l,m,n)) > abs(colmass_maxworst_r)) then
                 colmass_maxworst_r = colmass_worst(l,m,n) 
                 colmass_maxworst_i(1) = n
                 colmass_maxworst_i(2) = m
                 colmass_maxworst_i(3) = l
              end if
           enddo
        enddo ! size
        enddo ! type
        if ((icheck_colmass >= 10) .or. (abs(colmass_maxworst_r) >= 1.0e-6)) &
             write(*,'(a,3i3,1p,e10.2)') 'mixactivate colmass maxworst', &
             colmass_maxworst_i(1:3), colmass_maxworst_r
     endif ! (icheck_colmass > 0)

     return
   end subroutine Cmixactivate

!----------------------------------------------------------------------
subroutine cmaq_icenucleat (msectional,chem,num_chem,qv,qc,qi,qs,qnice3d,  &
           temp,p_phy,w,maxd_acomp,maxd_asize, maxd_atype, maxd_aphase,   &
           ncomp_aer, nsize_aer, ntype_aer, nphase_aer,  &
           numptr_aer, massptr_aer, dlo_sect, dhi_sect, sigmag_aer, dpvolmean,  &
           dens_aer, mw_aer,waterptr_aer, hygro,  ai_phase, cw_phase,  &
           !czhao 
           nihf3d,niimm3d,nidep3d,nimey3d,  &
           rndst1,rndst2,rndst3,rndst4,nacon1,nacon2,nacon3,nacon4,  &
           dgnuc,dgacc,dgcor,    &
           ids,ide, jds,jde, kds,kde,                            &
           ims,ime, jms,jme, kms,kme,                            &
           its,ite, jts,jte, kts,kte,                             &
	   rho, zm, dz8w, p_at_w, t_at_w, kvh,      &
           cldfra, cldfra_old,grid_id, ktau, dtstep, &
           f_qv,f_qc,f_qi,f_qs,f_ice_phy,f_rain_phy,sw_physics, &
!----------------------------------------------------------------------------
!Added by R Atlas for thompson indirect effect and icloud=3
           icloud, xland, dx, dy, num_wf, num_if, &
           morrison_indirect, thompson_indirect)
!******************************************************************
!ice nucleation is called in module_mixactivate_wrappers.F in chem/ 
! Subroutine to Compute Aerosol-Ice nucleation rate 
! * Author: Chun Zhao 
! * Originators of parts of code: Chun Zhao, Xiaohong Liu
! Last update: Jan 2011
! madificated by Shaocai Yu for CMAQ aerosol simulation (5/20/2011)
!
! Contact:
! Chun Zhao, PhD
! Scientist
! Pacific Northwest National Laboratory
! P.O. Box 999, MSIN K9-24
! Richland, WA, 99352
! Phone: (509) 371-6372
! Email: Chun.Zhao@pnl.gov
!
! Please report any bugs or problems to Chun Zhao
!
! Terms of Use:
!  1) Users are requested to consult the primary author prior to
!     modifying this module or incorporating it or its submodules in
!     another code. This is meant to ensure that the any linkages and/or
!     assumptions will not adversely affect the operation of this module.
!  2) Users preparing publications resulting from the usage of this code
!     are requested to cite one or more of the references below
!     (depending on the application) for proper acknowledgement.
!
! References: 
! * Liu and Penner (2005), .. , Meteorol. Z., 14, 499-514.
! * Liu et al. (2007), Inclusion of ice microphsics in the NCAR Community Atmospheric Model Version 3 (CAM3), J. Climate, 20, 4526-4547. 
!******************************************************************

  USE module_model_constants, only: g, rhowater, xlv, cp, rvovrd, r_d, r_v, mwdry, ep_2
  USE module_radiation_driver, only: cal_cldfra1,cal_cldfra2,cal_cldfra3
  USE module_state_description, ONLY : RRTMSCHEME,RRTMG_SWSCHEME,SWRADSCHEME   &
                                      ,GSFCSWSCHEME,GFDLSWSCHEME,CAMSWSCHEME

  implicit none

!     input
  INTEGER, intent(in) ::         grid_id,ktau
  INTEGER, intent(in) ::         num_chem
  integer, intent(in) ::         ids,ide, jds,jde, kds,kde,    &
                                 ims,ime, jms,jme, kms,kme,    &
                                 its,ite, jts,jte, kts,kte
				 
  integer, intent(in) :: maxd_aphase, nphase_aer, maxd_atype, ntype_aer
  integer, intent(in) :: maxd_asize, maxd_acomp, nsize_aer(maxd_atype)
  integer, intent(in) ::   &
       ncomp_aer( maxd_atype  ),   &
       massptr_aer( maxd_acomp, maxd_asize, maxd_atype, maxd_aphase ),   &
       waterptr_aer( maxd_asize, maxd_atype ),   &
       numptr_aer( maxd_asize, maxd_atype, maxd_aphase), &
       ai_phase, cw_phase
  integer, intent(in) :: msectional ! 1 for sectional, 0 for modal
  real, intent(in)  :: dens_aer( maxd_acomp, maxd_atype), & ! density (g/cm3)of material
       mw_aer( maxd_acomp, maxd_atype)       ! molecular weight (g/mole)
  real,intent(in),dimension(ims:ime, kms:kme,jms:jme,maxd_asize,maxd_atype) :: &
       dlo_sect,   & ! minimum size of section (cm)
       dhi_sect,   & ! maximum size of section (cm)
       sigmag_aer,   & ! geometric standard deviation of aerosol size dist
       dpvolmean   ! mean-volume diameter (cm) of mode
  real dpvolmean_aer(maxd_asize, maxd_atype)    ! mean-surface diameter (cm) of mode
  real dpsrfmean_aer(maxd_asize, maxd_atype)    ! mean-surface diameter (cm) of mode
! terminology:  (pi/6) * (mean-volume diameter)**3 ==(volume mixing ratio of section/mode)/(number mixing ratio)

  REAL, intent(in), DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ) :: &
       chem ! aerosol molar mixing ratio (ug/kg or #/kg)

  real raercol(kms:kme,num_chem) ! aerosol mass, number mixing ratios (ug/m3 or #/m3)
  real scale0

  REAL, intent(inout), DIMENSION( ims:ime, kms:kme, jms:jme ) :: &
       qv, qc, qi,qs ! water species (vapor, cloud drops, cloud ice) mixing ratio (g/g)

  REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),INTENT(IN   ):: &
      F_ICE_PHY, F_RAIN_PHY
      
  LOGICAL, OPTIONAL :: f_qv,f_qc, f_qi,f_qs

  INTEGER, INTENT(IN)      ::  sw_physics

  REAL, intent(inout), DIMENSION( ims:ime, kms:kme, jms:jme ) :: &
       qnice3d    ! ice nucleation number concentration (#/kg)

  REAL, intent(inout), DIMENSION( ims:ime, kms:kme, jms:jme ) :: &
       nihf3d,niimm3d,nidep3d,nimey3d

  !czhao
  REAL, intent(out), DIMENSION( ims:ime, kms:kme, jms:jme ) :: &
   nacon1,nacon2,nacon3,nacon4   ! number in 4 dust bins for contact freezing 
  REAL, intent(out), DIMENSION( ims:ime, kms:kme, jms:jme ) :: &
   rndst1,rndst2,rndst3,rndst4   ! radius of 4 dust bins for contact freezing 
  REAL, intent(in), DIMENSION( ims:ime, kms:kme, jms:jme ) :: &
   dgnuc,dgacc,dgcor 
   
  REAL, DIMENSION( 1:4 ) ::  nacon1d, rndst1d

  real, intent(in) :: dtstep             ! time step for microphysics (s)
  real, intent(in) :: temp(ims:ime, kms:kme, jms:jme)    ! temperature (K)
  real, intent(in) :: p_phy(ims:ime, kms:kme, jms:jme)    ! pressure 
  real, intent(in) :: w(ims:ime, kms:kme, jms:jme)   ! vertical velocity (m/s)
  real, intent(in) :: rho(ims:ime, kms:kme, jms:jme)    ! density at mid-level  (kg/m3)
  real, intent(in) :: zm(ims:ime, kms:kme, jms:jme)     ! geopotential height of level (m)
  real, intent(in) :: dz8w(ims:ime, kms:kme, jms:jme) ! layer thickness (m)
  real, intent(in) :: p_at_w(ims:ime, kms:kme, jms:jme) ! pressure at layer interface (Pa)
  real, intent(in) :: t_at_w(ims:ime, kms:kme, jms:jme) ! temperature at layer interface (K)
  real, intent(in) :: kvh(ims:ime, kms:kme, jms:jme)    ! vertical diffusivity (m2/s)
  real, intent(inout) :: cldfra_old(ims:ime, kms:kme, jms:jme)! cloud fraction on previous time step
  real, intent(inout) :: cldfra(ims:ime, kms:kme, jms:jme)    ! cloud fraction
  real, intent(in) :: hygro( its:ite, kts:kte, jts:jte, maxd_asize, maxd_atype ) ! bulk hygroscopicity   &
!--------------------Local storage-------------------------------------
!
  real :: dgnum_m(maxd_asize, maxd_atype) ! median diameter (m) of number distrib of mode
  real :: qnice(kms:kme)      ! ice number mixing ratio (#/kg)
  real :: nihf2,niimm2,nidep2,nimey2,dum2 !dummy variables for activated ice nuclei by diffferent processes
  real :: icldfra(kms:kme)     ! ice cloud fraction
  real :: wtke(kms:kme)        ! turbulent vertical velocity at base of layer k (m/s)
  real :: wsubi(kms:kme)       ! sub grid vertical velocity for ice nucleation (m/s)
  real zn(kms:kme)             ! g/pdel (m2/g) for layer
  real zs(kms:kme)             ! inverse of distance between levels (m)
  real, parameter :: zkmin = 0.01
  real, parameter :: zkmax = 100.
  real cs(kms:kme)             ! air density (kg/m3) at layer center
  real csbot(kms:kme)          ! air density (kg/m3) at layer bottom
  real csbot_cscen(kms:kme)    ! csbot(k)/cs(k)
  real dz(kms:kme)             ! geometric thickness of layers (m)
  
  real wdiab                   ! diabatic vertical velocity
  real, parameter :: wmixmin = 0.001 ! minimum turbulence vertical velocity (m/s)
  real, parameter :: wmixmax = 0.2 ! maximum turbulence vertical velocity (m/s)
  real :: ekd(kms:kme)         ! diffusivity for ice (m2/s)
  real :: ekk(kms:kme)         ! density*diffusivity for ice (kg/m3 m2/s)
  real, parameter :: sq2pi = 2.5066282746

  real wbar,wmix,wmin,wmax
  real dum
  integer i,j,k,m,n,nsub
  real qcld
  real pi
  real ::  ekkp(kms:kme),ekkm(kms:kme) ! zn*zs*density*diffusivity

  integer l,ltype
  integer :: ntype(maxd_asize)

  integer ncomp(maxd_atype)
  integer lnum,lmass

  real :: hygro_aer(maxd_asize, maxd_atype)  ! hygroscopicity of aerosol mode
  real :: exp45logsig     ! exp(4.5*alogsig**2)
  real :: alogsig(maxd_asize, maxd_atype) ! natl log of geometric standard dev of aerosol

  real :: qsat            ! liquid-ice weighted sat mixing rat (kg/kg)
  real :: esat            ! liquid-ice weighted sat vapor press (pa)
! real :: eslq(kts:kte)    ! liquid sat vapor pressure (pa)
! real :: esic(kts:kte)    ! ice sat vapor pressure (pa)
  real :: gammas          ! parameter for cond/evap of cloud water
  real, parameter :: tmelt = 273.15     ! freezing T of fresh water (K)
  real :: relhum1d(kts:kte)
  real :: nfice(kts:kte)  !fice variable
  real :: dumfice
  real, parameter :: qsmall = 1.e-18   !min mixing ratio  
  real, parameter :: mincld = 1.e-7   !min cloud 

  CHARACTER (LEN=200) :: message

!----------------------------------------------------------------------------
!Added by R Atlas for thompson indirect effect

  integer, intent(in) :: icloud
  real, dimension(ims:ime,jms:jme), intent(in):: xland
  real, intent(in) :: dx, dy
  real, intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: &
                num_wf, num_if
  logical, intent(in) :: morrison_indirect, thompson_indirect
  real :: num_wf1d, num_if1d
!----------------------------------------------------------------------------

 !    write(6,*)'inside cmaq_icenucleat here1'
!    write(6,*)'chem(its,kts,jts,2) in cmaq_icenucleat',chem(its,kts,jts,2)

!    write(6,*)'inside cmaq_icenucleat here1a',254,85,26,num_chem,chem(254,26,85,1:num_chem)

      raercol = 0.0   

OVERALL_MAIN_J_LOOP: do j=jts,jte
OVERALL_MAIN_I_LOOP: do i=its,ite

     do k=kts+1,kte
      zs(k)=1./(zm(i,k,j)-zm(i,k-1,j))
     enddo
      zs(kts)=zs(kts+1)
      zs(kte+1)=0.

      qnice(kts:kte)    = 0.
!---------------------------------------------------------------
!diagnose ice cloud fraction 
!---------------------------------------------------------------

     cldfra_select: SELECT CASE(sw_physics)

     CASE (RRTMG_SWSCHEME)
 
       IF ( ICLOUD .LT. 3 ) THEN  ! R Atlas
 
          IF (PRESENT(F_QC) .AND. PRESENT ( F_QI ) ) THEN
             CALL cal_cldfra1(CLDFRA,qv,qc,qi,qs,          &
                   F_QV,F_QC,F_QI,F_QS,temp,p_phy,         &
                   F_ICE_PHY,F_RAIN_PHY,                   &
                   ids,ide, jds,jde, kds,kde,              &
                   ims,ime, jms,jme, kms,kme,              &
                   its,ite, jts,jte, kts,kte               )
          ENDIF

!----------------------------------------------------------------------------
!Added by R Atlas for use with thompson cloud fraction scheme

         ELSE IF ( ICLOUD == 3 ) THEN

             CALL cal_cldfra3(CLDFRA, qv, qc, qi, qs,                  &
       &              p_phy,temp,rho, XLAND, dx, dy,                   &
       &              ids,ide, jds,jde, kds,kde,                       &
       &              ims,ime, jms,jme, kms,kme,                       &
       &              its,ite, jts,jte, kts,kte)

       END IF

!-------------------------------------------------------------------------

     CASE DEFAULT

        IF (PRESENT(F_QC) .AND. PRESENT ( F_QI ) ) THEN
          CALL cal_cldfra2(CLDFRA,qc,qi,f_qc,f_qi,  &
            ids,ide, jds,jde, kds,kde,              &
            ims,ime, jms,jme, kms,kme,              &
            its,ite, jts,jte, kts,kte               )
        ENDIF

     END SELECT cldfra_select

     do k=kts,kte
        if(f_qi)then
           qcld=qc(i,k,j)+qi(i,k,j)
        else
           qcld=qc(i,k,j)
        endif
        if(qcld.lt.-1..or.qcld.gt.1.)then
           write(6,'(a,g12.2,a,3i5)')'qcld=',qcld,' for i,k,j=',i,k,j
           call wrf_error_fatal("1")
        endif
        if(qcld.gt.1.e-20)then
           icldfra(k)=cldfra(i,k,j)*qi(i,k,j)/qcld
        else
           icldfra(k)=0.
        endif
       
       ! get cloud fraction, check for minimum
       icldfra(k)=max(icldfra(k),mincld)
     enddo
!    write(6,*)'inside cmaq_icenucleat here2'
    
!---------------------------------------------------------------
!diagnose sub-grid vertical velocity
!---------------------------------------------------------------
     do k=kts,kte
        cs(k)=rho(i,k,j) ! air density (kg/m3)
        dz(k)=dz8w(i,k,j) ! layer thickness (m)
        zn(k)=1./(cs(k)*dz(k))
        if(k>kts)then
           ekd(k)=kvh(i,k,j)
           ekd(k)=max(ekd(k),zkmin)
           ekd(k)=min(ekd(k),zkmax)
        else
           ekd(k)=0
        endif

        !diagnose subgrid vertical velocity from diffusivity
        wtke(k)=sq2pi*ekd(k)/dz(k)
!       wsubi(k)=w(i,k,j)+wtke(k)
        wsubi(k)=wtke(k)
        !set upper and lower bounds
        wsubi(k)=max(wmixmin,wsubi(k))
        wsubi(k)=min(wsubi(k),wmixmax)

     enddo
     zn(kte+1)          = 0.

!--------------------------------------------------------------------------
!Get humidity and saturation vapor pressures
     do k=kts,kte

! find temperature and saturation value for provisional t and q without
! condensation

      call vqsatd_water(temp(i,k,j),p_phy(i,k,j),esat,qsat,gammas) 
        relhum1d(k)=qv(i,k,j)/qsat


!-------------------------------------------
! not used
!      eslq(k)=polysvp1(temp(i,k,j),0)
!      esic(k)=polysvp1(temp(i,k,j),1)
!! hm fix, make sure when above freezing that esic=eslq, not active yet
!      if (temp(i,k,j).gt.tmelt) esic(k)=eslq(k)
!-------------------------------------------

! calculate nfice based on liquid and ice mmr (no rain and snow mmr available yet)
       nfice(k)=0.
       dumfice=qc(i,k,j)+qi(i,k,j)
       if (dumfice.gt.qsmall .and. qi(i,k,j).gt.qsmall) then
          nfice(k)=qi(i,k,j)/dumfice
       endif

     enddo ! k

!------------------------------------------------------------------
!  calculate num and mass mixing ratio for aerosol
!------------------------------------------------------------------

     do k=kts,kte 
      do n=1,ntype_aer
           ncomp(n)=ncomp_aer(n)
      do m=1,nsize_aer(n)
           lnum=numptr_aer(m,n,ai_phase)
           if(lnum>0)then
              scale0 = cs(k) 
              raercol(k,lnum)=chem(i,k,j,lnum)*scale0 ! #/kg to #/m3
            endif
           do l=1,ncomp(n)
              lmass=massptr_aer(l,m,n,ai_phase)
              scale0 = cs(k)*1.e-9 
              raercol(k,lmass)=chem(i,k,j,lmass)*scale0 ! ug/kg to kg/m3
           enddo
      enddo ! size
      enddo ! type
     enddo  ! k

!    write(6,*)'inside cmaq_icenucleat here3',i,j

!------------------------------------------------------------------
!ICE Nucleation following CAM
     do k=kts,kte 
     ! initialize
        nacon1(i,k,j)=0.0
        nacon2(i,k,j)=0.0
        nacon3(i,k,j)=0.0
        nacon4(i,k,j)=0.0
        rndst1(i,k,j)=1.e-7 ! 0.1 um to avoid divid by zero
        rndst2(i,k,j)=1.e-7
        rndst3(i,k,j)=1.e-7
        rndst4(i,k,j)=1.e-7

      if ((temp(i,k,j).lt.tmelt - 5.) .or. (thompson_indirect)) then

       !--------------------------------------------------------
       !Get medium diameter 
       pi = 4.*atan(1.0)
       do n=1,ntype_aer
       do m=1,nsize_aer(n)
         alogsig(m,n)=alog(sigmag_aer(i,k,j,m,n))
         !now use predicted size instead of prescribed size --czhao
!        dgnum_m(m,n) = (dpvolmean_aer(m,n)*1.e-2) * exp( -1.5*alogsig(m,n)*alogsig(m,n) ) !m
!        dpsrfmean_aer(m,n) = (dgnum_m(m,n)*1.e2) * exp(alogsig(m,n)*alogsig(m,n) ) !cm
         if (n==1.and.m==1) then ! Aitken mode... 
         dgnum_m(m,n) = dgnuc(i,k,j)  ! m 
         dpvolmean_aer(m,n) = (dgnum_m(m,n)*1.e2) * exp(1.5*alogsig(m,n)*alogsig(m,n) ) !cm 
         dpsrfmean_aer(m,n) = (dgnum_m(m,n)*1.e2) * exp(alogsig(m,n)*alogsig(m,n) ) !cm
         endif
!         if (n==1.and.m==2) then ! Accumulation mode... 
         if (n==2.and.m==1) then ! Accumulation mode for CMAQ aerosol... (Shaocai Yu)
         dgnum_m(m,n) = dgacc(i,k,j)  ! m 
         dpvolmean_aer(m,n) = (dgnum_m(m,n)*1.e2) * exp(1.5*alogsig(m,n)*alogsig(m,n) ) !cm 
         dpsrfmean_aer(m,n) = (dgnum_m(m,n)*1.e2) * exp(alogsig(m,n)*alogsig(m,n) ) !cm
         endif
!         if (n==2.and.m==1) then ! Coarse mode...
         if (n==3.and.m==1) then ! Coarse mode for CMAQ aerosol...(Shaocai Yu)	 
         dgnum_m(m,n) = dgcor(i,k,j)  ! m 
         dpvolmean_aer(m,n) = (dgnum_m(m,n)*1.e2) * exp(1.5*alogsig(m,n)*alogsig(m,n) ) !cm 
         dpsrfmean_aer(m,n) = (dgnum_m(m,n)*1.e2) * exp(alogsig(m,n)*alogsig(m,n) ) !cm
         endif

       enddo
       enddo
!    write(6,*)'inside cmaq_icenucleat here4',i,j,k,num_chem,raercol(k,1:num_chem)
    
        nacon1d(1:4)=0. ! initialize 
        rndst1d(1:4)=1.e-7 ! 0.1 um to avoid divid by zero
        num_wf1d=0 !R Atlas
        num_if1d=0 !R Atlas

       ! get ice nucleation rate
       call cmaq_nucleati_cam(wsubi(k),temp(i,k,j),p_phy(i,k,j),&
            relhum1d(k),icldfra(k),qc(i,k,j),nfice(k),cs(k),&
            raercol(k,1:num_chem),dpvolmean_aer,dpsrfmean_aer, &
            maxd_asize,maxd_atype,num_chem,ai_phase,dum2,nihf2,niimm2,&
            nidep2,nimey2, & 
       !for contact freezing in MG microphysics
            rndst1d,nacon1d, &
!-----------------------------------------------------------
!Added by R Atlas for indirect effect
            rho(i,k,j),num_wf1d,num_if1d )
!-----------------------------------------------------------
       
           qnice(k)=dum2
           nihf3d(i,k,j)=nihf2
           niimm3d(i,k,j)=niimm2
           nidep3d(i,k,j)=nidep2
           nimey3d(i,k,j)=nimey2
           
           nacon1(i,k,j)=nacon1d(1)
           nacon2(i,k,j)=nacon1d(2)
           nacon3(i,k,j)=nacon1d(3)
           nacon4(i,k,j)=nacon1d(4)
           rndst1(i,k,j)=rndst1d(1)
           rndst2(i,k,j)=rndst1d(2)
           rndst3(i,k,j)=rndst1d(3)
           rndst4(i,k,j)=rndst1d(4)
!-----------------------------------------------------------
!Added by R Atlas for indirect effect
           num_wf(i,k,j)=num_wf1d
           num_if(i,k,j)=num_if1d
!-----------------------------------------------------------
      endif ! temp
!---------------------------------
!  keep the old ice number (Shaocai Yu)
!------------------------------------------------
!       if (temp(i,k,j).gt.tmelt) then
!         qnice3d(i,k,j)= 0.0
!       else      
!         qnice3d(i,k,j)=qnice3d(i,k,j)+max(qnice(k),1.e-10)
!       endif
!-----------------------------------------------
!     Original code
!--------------------------------------------
       qnice3d(i,k,j)=max(qnice(k),1.e-10)
!-------------------------------------------------
       
     enddo  ! k

 
enddo OVERALL_MAIN_I_LOOP ! end of main loop over i
enddo OVERALL_MAIN_J_LOOP ! end of main loop over j
   
	write(*,*)'in cmaq_icenucleat22222'			 
  end subroutine cmaq_icenucleat 
!----------------------------------------------------------------------

!----------------------------------------------------------------------
  subroutine cmaq_nucleati_cam(wbar,tair,pmid,relhum,cldn,qc,nfice,&
             rhoair,chem1d,dpvolmean,dpsrfmean,maxd_asize,maxd_atype,&
             num_chem,ai_phase,nuci,onihf,oniimm,onidep,onimey, &
             !for contact freezing in MG microphysics
             rndst1d, nacon1d,  & 
!-----------------------------------------------------------
!Added by R Atlas for indirect effect
             dens, num_wf1d, num_if1d)
!---------------------------------------------------------------
! Purpose:
!  The parameterization of ice nucleation.
!
! Method: The current method is based on Liu & Penner (2005)
!  It related the ice nucleation with the aerosol number, temperature and the
!  updraft velocity. It includes homogeneous freezing of sulfate, immersion
!  freezing of soot, and Meyers et al. (1992) deposition nucleation
!
! Authors: Xiaohong Liu, 01/2005, modifications by A. Gettelman 2009-2010
!          modifications by Chun Zhao 01/2011 for WRF-Chem
!          modificated by Shaocai Yu 5/20/2011 for WRF-CMAQ
!----------------------------------------------------------------
 
  implicit none
!-----------------------------------------------------
! Input Arguments
!
  INTEGER, intent(in) ::  num_chem,maxd_atype,maxd_asize,ai_phase

  real,intent(in) :: wbar    ! subgrid vertical velocity for ice nucleation (m/s)
  real,intent(in) :: tair    ! air temperature
  real,intent(in) :: pmid    ! air pressure 
  real,intent(in) :: relhum  ! relative humidity
  real,intent(in) :: cldn    ! ice cloud fraction
  real,intent(in) :: qc      ! liquid water mixing ratio (kg/kg)
  real,intent(in) :: nfice   ! ice mass fraction
  real,intent(in) :: rhoair  ! air density (kg/m3)
  real,intent(in) :: chem1d(1:num_chem) ! aerosol number and mass mixing ratios (#/m3 or kg/m3) 
  real,intent(in) :: dpvolmean(maxd_asize, maxd_atype)  ! mean-volume diameter (cm) of mode
  real,intent(in) :: dpsrfmean(maxd_asize, maxd_atype)  ! mean-surface diameter (cm) of mode

!
! Output Arguments
!
  real, intent(out) :: nuci   ! ice number nucleated (#/kg)
  real, intent(out) :: onihf  ! nucleated number from homogeneous freezing of so4
  real, intent(out) :: oniimm ! nucleated number from immersion freezing
  real, intent(out) :: onidep ! nucleated number from deposition nucleation
  real, intent(out) :: onimey ! nucleated number from deposition nucleation  (meyers: mixed phase)

  REAL, intent(out), DIMENSION( 1:4 ) ::  nacon1d, rndst1d

!
! Local workspace
!
  real dpvolmean_m(maxd_asize, maxd_atype)  ! mean-volume diameter (m) of mode
  real dpsrfmean_m(maxd_asize, maxd_atype)  ! mean-surface diameter (m) of mode
  real so4_num           ! so4 aerosol number (#/cm^3)
  real soot_num          ! soot (hydrophilic) aerosol number (#/cm^3)
  real organic_num       ! organic (hydrophobic) aerosol number (#/cm^3)
  real dst_num                    ! total dust aerosol number (#/cm^3) 
  real soot_sfc, organic_sfc, dst_sfc   ! aerosol surface area (m2/cm^3)

  real nihf             ! nucleated number from homogeneous freezing of so4
  real niimm            ! nucleated number from immersion freezing
  real nidep            ! nucleated number from deposition nucleation
  real nimey            ! nucleated number from deposition nucleation (meyers)
  real ni,n1               ! nucleated number
  real tc,A,B,C,regm,RHw
  real subgrid
! Phillips et al. for mixed-phase cloud ice nucleation
  logical, parameter :: phi_in = .true.
  real  A_dust, n_dust,  &
            A_soot, n_soot,  &
            A_organic_phobic, num_organic
  real  num_ic_dust_imm, num_ic_soot_imm, num_ic_organic_imm
  real  esl,esi,deles
!-----------------------------------------------------------
!Added by R Atlas for thompson indirect effect
  real, intent(in) :: dens
  real, intent(out) :: num_wf1d, num_if1d
!-----------------------------------------------------------
  CHARACTER (LEN=200) :: message

    so4_num=0.0
    soot_num=0.0
    organic_num=0.0
    dst_num=0.0

    soot_sfc = 0.0
    organic_sfc = 0.0
    dst_sfc = 0.0
    
!    write(6,*)'inside cmaq_nucleati_cam  here1',chem1d

    ! get aerosol number and mass for modal aerosols 
    dpvolmean_m=dpvolmean*1.e-2
    dpsrfmean_m=dpsrfmean*1.e-2
    call get_modal_aer_cmaq(chem1d,num_chem,dpvolmean_m,dpsrfmean_m, & 
    so4_num,soot_num,organic_num,dst_num,soot_sfc,organic_sfc,dst_sfc, & 
    nacon1d,rndst1d,dens,num_wf1d,num_if1d) 

!    write(6,*)'inside cmaq_nucleati_cam  here2'

    ni=0.
    tc=tair-273.15

    ! initialize
    niimm=0.
    nidep=0.
    nihf=0.

    if(so4_num.ge.1.0e-10.or.dst_num.ge.1.0e-10.and.cldn.gt.0.) then
    !-----------------------------
    ! homogeneous nucelationi & heterogeneous immersion nucleation 
    ! in pure ice clouds
     subgrid = 1.2
     if ((relhum*polysvp1(tair,0)/polysvp1(tair,1)*subgrid).ge.1.2) then  ! supersaturation
     if ( ((tc.le.0.0).and.(tc.ge.-37.0).and.(qc.lt.1.e-12)).or.(tc.le.-37.0) ) then ! pure ice clouds
       A = -1.4938 * log(dst_num) + 12.884
       B = -10.41  * log(dst_num) - 67.69
       regm = A * log(wbar) + B
       if(tc.gt.regm) then    ! heterogeneous nucleation only
         if(tc.lt.-40. .and. wbar.gt.1.) then ! exclude T<-40 & W>1m/s from hetero. nucleation
           call hf(tc,wbar,relhum,subgrid,so4_num,nihf)
           niimm=0.
           nidep=0.
           n1=nihf
         else
           call hetero(tc,wbar,dst_num,niimm,nidep)
           nihf=0.
           n1=niimm+nidep
         endif
       elseif (tc.lt.regm-5.) then ! homogeneous nucleation only
         call hf(tc,wbar,relhum,subgrid,so4_num,nihf)
         niimm=0.
         nidep=0.
         n1=nihf
       else        ! transition between homogeneous and heterogeneous: interpolate in-between
         if(tc.lt.-40. .and. wbar.gt.1.) then ! exclude T<-40 & W>1m/s from hetero. nucleation
           call hf(tc,wbar,relhum,subgrid,so4_num,nihf)
           niimm=0.
           nidep=0.
           n1=nihf
         else

           call hf(regm-5.,wbar,relhum,subgrid,so4_num,nihf)
           call hetero(regm,wbar,dst_num,niimm,nidep)

           if(nihf.le.(niimm+nidep)) then
             n1=nihf
           else
             n1=(niimm+nidep)*((niimm+nidep)/nihf)**((tc-regm)/5.)
           endif
         endif
       endif

       ni=n1

     endif
     endif

    endif  ! cldn

! deposition/condensation nucleation in mixed clouds (-37<T<0C)

    if (phi_in) then  ! Phillips et al. 2008 use particles number and area with D>0.1 um

      n_dust = dst_num*1.e6/rhoair   ! #/kg
      A_dust = dst_sfc*1.e6/rhoair   ! m2/kg

      n_soot = soot_num*1.e6/rhoair   ! #/kg
      A_soot = soot_sfc*1.e6/rhoair   ! m2/kg

      num_organic = organic_num*1.e6/rhoair      ! #/kg
      A_organic_phobic = organic_sfc*1.e6/rhoair ! m2/kg

    endif

    if(tc.lt.0. .and. tc.gt.-37. .and. qc.gt.1.e-12) then

      if (phi_in) then      ! use Phillips et al

        call nucleatiphi(tair, pmid, n_dust, A_dust, n_soot, A_soot, & 
        num_organic, A_organic_phobic, num_ic_dust_imm, num_ic_soot_imm, &
        num_ic_organic_imm)

        nimey=(num_ic_dust_imm+num_ic_soot_imm+num_ic_organic_imm)*rhoair*1.e-6  ! cm^-3

      else                  ! use Meyers et al. 1992
        esl = polysvp1(tair,0)     ! over water in mixed clouds
        esi = polysvp1(tair,1)     ! over ice
        deles = (esl - esi)
        nimey=1.e-3*exp(12.96*deles/esi - 0.639)

      endif

    else
      nimey=0.
    endif

    nuci=ni+nimey
    if(nuci.gt.9999.0.or.nuci.lt.0.) then
        write(message,'(A)') 'Warning: incorrect ice nucleation number (nuci reset =0)' 
        CALL wrf_message( TRIM(message) )

!       write(message,' ') ni, tair, relhum, wbar, nihf, niimm, nidep,deles,esi,dst2_num,dst3_num,dst4_num
        nuci=0.
    endif

    nuci=nuci*1.e6/rhoair    ! change unit from #/cm3 to #/kg
    onimey=nimey*1.e6/rhoair
    onidep=nidep*1.e6/rhoair
    oniimm=niimm*1.e6/rhoair
    onihf=nihf*1.e6/rhoair

  return

  end subroutine cmaq_nucleati_cam
!----------------------------------------------------------------------

!----------------------------------------------------------------------
  subroutine get_modal_aer_cmaq(chem1d,num_chem,dpvolmean,dpsrfmean, &
    so4_num,soot_num,organic_num,dst_num,soot_sfc,organic_sfc,dst_sfc, & 
    nacon1d,rndst1d,&
!-----------------------------------------------------------
!Added by R Atlas for indirect effect
    dens,num_wf1d,num_if1d) 
!---------------------------------------------------------------
!          modificated by Shaocai Yu 5/20/2011 for WRF-CMAQ
!----------------------------------------------------------------

  USE module_data_cmaq !
  USE module_aerosols_cmaq_init_aercld !
  USE module_state_description, only:  param_first_scalar
 
  implicit none
!-----------------------------------------------------
! Input Arguments
!
  INTEGER, intent(in) ::  num_chem
  real,intent(in) :: chem1d(1:num_chem) ! aerosol number and mass mixing ratios (#/m3 or kg/m3) 
  real,intent(in) :: dpvolmean(maxd_asize, maxd_atype)  ! mean-volume diameter (m) of mode
  real,intent(in) :: dpsrfmean(maxd_asize, maxd_atype)  ! mean-surface diameter (m) of mode

!-----------------------------------------------------
! Output Arguments
!
  real, intent(inout) :: so4_num,soot_num,organic_num,dst_num, & 
                            soot_sfc,organic_sfc,dst_sfc 
  REAL, intent(out), DIMENSION( 1:4 ) ::  nacon1d, rndst1d

!-----------------------------------------------------
! Local Arguments
!
  integer isize,itype
  integer p1st 
  integer l
  !Aiken
  real  mass_so4i , mass_oci ,  mass_bci  , mass_dusti, mass_orgpai
	
  !Accumulation
  real  mass_so4j , mass_ocj ,  mass_bcj  , mass_dustj, mass_orgpaj, &
        mass_alkj,mass_xyl1j,mass_xyl2j,mass_xyl3j,mass_tol1j, &
        mass_tol2j,mass_tol3j,mass_bnz1j,mass_bnz2j,mass_bnz3j, &
        mass_trp1j,mass_trp2j,mass_iso1j,mass_iso2j,     &
        mass_sqtj,mass_ISO3j,mass_OLGAj,mass_OLGBj,mass_orgcj,mass_paj
	
  !Coarse      
  real  mass_so4k, mass_soil, &
        mass_acor,mass_soilj !R Atlas

  real  vol_so4i   , vol_oci   ,  vol_bci    , vol_dusti, &
        vol_so4j   , vol_ocj   ,  vol_bcj    , vol_dustj, &
        vol_so4k   , vol_soil  ,  &
        vol_acor   , vol_soilj !R Atlas

  real  sfc_so4i   , sfc_oci   ,  sfc_bci    , sfc_dusti, &
        sfc_so4j   , sfc_ocj   ,  sfc_bcj    , sfc_dustj, &
        sfc_so4k   , sfc_soil 

  real  num_so4i   , num_oci   ,  num_bci    , num_dusti, &
        num_so4j   , num_ocj   ,  num_bcj    , num_dustj, &
        num_so4k   , num_soil  ,  &
        num_acor   , num_soilj !R Atlas

  real  dpvol_i,dpvol_j,dpvol_c

!-----------------------------------------------------------
!Added by R Atlas for thompson indirect effect
  real, intent(in) :: dens
  real, intent(out) :: num_wf1d, num_if1d
!-----------------------------------------------------------

!    write(6,*)'inside get_modal_aer_cmaq  here1',num_chem
!    write(6,*) chem1d   
!    
!        p1st = param_first_scalar
!-------------------Shaocai
        p1st = 1
!--------------------------------
	
        ! Aitken mode...
        isize = 1 ; itype = 1 
	
	mass_so4i = chem1d(p_aso4ai)
        mass_orgpai = chem1d(p_aorgpaai)
        mass_bci = chem1d(p_aecai)	
	
        mass_oci = mass_orgpai
        mass_dusti = 0.0    ! no dust for aitken mode so far for CMAQ

        vol_so4i = mass_so4i/(dens_aso4_aer*1.e3)   ! m3/m3
        vol_oci = mass_oci/(dens_aorgc_aer*1.e3)
        vol_bci = mass_bci/(dens_aec_aer*1.e3)
        vol_dusti = mass_dusti/(dens_asoil_aer*1.e3)

        num_so4i = vol_so4i*f6dpi/(dpvolmean(isize,itype)**3)  ! #/m3
        num_oci = vol_oci*f6dpi/(dpvolmean(isize,itype)**3)
        num_bci = vol_bci*f6dpi/(dpvolmean(isize,itype)**3)
        num_dusti = vol_dusti*f6dpi/(dpvolmean(isize,itype)**3)

        sfc_so4i = num_so4i*pirs*(dpsrfmean(isize,itype)**2) ! m2/m3
        sfc_oci = num_oci*pirs*(dpsrfmean(isize,itype)**2)
        sfc_bci = num_bci*pirs*(dpsrfmean(isize,itype)**2)
        sfc_dusti = num_dusti*pirs*(dpsrfmean(isize,itype)**2)

        dpvol_i=dpvolmean(isize,itype)
	
!    write(6,*)'inside get_modal_aer_cmaq  here1a'
    
! Accumulation mode...
        isize = 1 ; itype = 2  
        mass_so4j= chem1d(p_aso4aj)
        mass_alkj= chem1d(p_aalkaj)
        mass_xyl1j= chem1d(p_axyl1aj)
        mass_xyl2j= chem1d(p_axyl2aj)
        mass_xyl3j= chem1d(p_axyl3aj)
        mass_tol1j= chem1d(p_atol1aj)
        mass_tol2j= chem1d(p_atol2aj)
        mass_tol3j= chem1d(p_atol3aj)
        mass_bnz1j= chem1d(p_abnz1aj)	
        mass_bnz2j= chem1d(p_abnz2aj)	
        mass_bnz3j= chem1d(p_abnz3aj)	
        mass_trp1j= chem1d(p_atrp1aj)	
        mass_trp2j= chem1d(p_atrp2aj)	
        mass_iso1j= chem1d(p_aiso1aj)	
        mass_iso2j= chem1d(p_aiso2aj)	
        mass_sqtj= chem1d(p_asqtaj)	
        mass_ISO3j= chem1d(p_aISO3aj)	
        mass_OLGAj= chem1d(p_aOLGAaj)	
        mass_OLGBj= chem1d(p_aOLGBaj)	
        mass_orgcj= chem1d(p_aorgcaj)	
        mass_orgpaj = chem1d(p_aorgpaaj)
        mass_bcj= chem1d(p_aecaj)
	
        mass_ocj=mass_alkj+mass_xyl1j+mass_xyl2j+mass_xyl3j+mass_tol1j   &
                 +mass_tol2j+mass_tol3j+mass_bnz1j+mass_bnz2j+mass_bnz3j &
                 +mass_trp1j+mass_trp2j+mass_iso1j+mass_iso2j+mass_sqtj  &
                 +mass_ISO3j+mass_OLGAj+mass_OLGBj+mass_orgcj+mass_orgpaj
        mass_dustj = 0.0    ! no dust for aitken mode so far for CMAQ
        mass_soilj= chem1d(p_asoilaj) !R Atlas

        vol_so4j = mass_so4j/(dens_aso4_aer*1.e3)  ! m3/m3
        vol_ocj = mass_ocj/(dens_aorgc_aer*1.e3)   ! m3/m3
        vol_bcj = mass_bcj/(dens_aec_aer*1.e3)   ! m3/m3
        vol_dustj = mass_dustj/(dens_asoil_aer*1.e3)  !m3/m3
        vol_soilj = mass_soilj/(dens_asoil_aer*1.e3)   !m3/m3 R Atlas

        num_so4j = vol_so4j*f6dpi/(dpvolmean(isize,itype)**3)  !#/m3
        num_ocj = vol_ocj*f6dpi/(dpvolmean(isize,itype)**3)
        num_bcj = vol_bcj*f6dpi/(dpvolmean(isize,itype)**3)
        num_dustj = vol_dustj*f6dpi/(dpvolmean(isize,itype)**3)
        num_soilj = vol_soilj*f6dpi/(dpvolmean(1,2)**3)  !#/m3 R Atlas

        sfc_so4j = num_so4j*pirs*(dpsrfmean(isize,itype)**2)   !m2/m3
        sfc_ocj = num_ocj*pirs*(dpsrfmean(isize,itype)**2)
        sfc_bcj = num_bcj*pirs*(dpsrfmean(isize,itype)**2)
        sfc_dustj = num_dustj*pirs*(dpsrfmean(isize,itype)**2)

        dpvol_j=dpvolmean(isize,itype)
	
!    write(6,*)'inside get_modal_aer_cmaq  here2'

! Coarse mode...
        isize = 1 ; itype = 3   ! after  march-2008 ordering

        mass_so4k= chem1d(p_aso4ak)
        mass_soil= chem1d(p_asoilak)
        mass_acor= chem1d(p_acorsak) !R Atlas

        vol_so4k = mass_so4k/(dens_aso4_aer*1.e3)  ! m3/m3
        vol_soil = mass_soil/(dens_asoil_aer*1.e3)   !m3/m3
        vol_acor = mass_acor/(dens_acors_aer*1.e3)   !m3/m3 R Atlas

        num_so4k = vol_so4k*f6dpi/(dpvolmean(isize,itype)**3)  !#/m3
        num_soil = vol_soil*f6dpi/(dpvolmean(isize,itype)**3)  !#/m3
        num_acor = vol_acor*f6dpi/(dpvolmean(isize,itype)**3)  !#/m3 R Atlas

        sfc_so4k = num_so4k*pirs*(dpsrfmean(isize,itype)**2)  !m2/m3
        sfc_soil = num_soil*pirs*(dpsrfmean(isize,itype)**2)  !m2/m3

        so4_num = (num_so4i + num_so4j+num_so4k)*1.e-6   ! #/cm3 
        organic_num = (num_oci + num_ocj)*1.e-6 
        soot_num = (num_bci + num_bcj)*1.e-6 
        dst_num = (num_dusti + num_dustj + num_soil)*1.e-6 

        organic_sfc = (sfc_oci + sfc_ocj)*1.e-6  ! m2/cm3 
        soot_sfc = (sfc_bci + sfc_bcj)*1.e-6 
        dst_sfc = (sfc_dusti + sfc_dustj + sfc_soil)*1.e-6

        dpvol_c=dpvolmean(isize,itype)

!-----------------------------------------------------------
!Added by R Atlas for indirect effect
!Convert from #/m^3 to #/kg
        num_if1d = (num_soil + num_acor + num_soilj)/dens
        num_wf1d = (chem1d(45) + chem1d(46) -num_bcj)/dens - num_if1d
!--------------------------------------------------------------------------------
! Contact freezing  (-40<T<-3 C) (Young, 1974) with hooks into simulated dust
! estimate rndst and nacon for 4 dust bins here to pass to MG microphysics
! rndst: radius of 4 dust bins for contact freezing
! nacon: number in 4 dust bins for contact freezing 
! the estimate here has large uncertainty, can be improved in future

       ! For modal aerosols: use size '2' for dust accumulation mode, and 
       ! size '3' for dust coarse mode
         nacon1d(2)=num_dustj   ! #/m3
         nacon1d(3)=num_soil

       !also redefine parameters based on size...
         rndst1d(2)=0.5*dpvol_j ! m
         rndst1d(3)=0.5*dpvol_c

  return

  end subroutine get_modal_aer_cmaq 
!----------------------------------------------------------------------
  subroutine hf(T,ww,RH,subgrid,Na,Ni)

      real, intent(in)  :: T, ww, RH, subgrid, Na
      real, intent(out) :: Ni 

      real    A1_fast,A21_fast,A22_fast,B1_fast,B21_fast,B22_fast
      real    A2_fast,B2_fast
      real    C1_fast,C2_fast,k1_fast,k2_fast
      real    A1_slow,A2_slow,B1_slow,B2_slow,B3_slow
      real    C1_slow,C2_slow,k1_slow,k2_slow
      real    regm
      real    A,B,C
      real    RHw

!---------------------------------------------------------------------
! parameters

      A1_fast  =0.0231
      A21_fast =-1.6387  !(T>-64 deg)
      A22_fast =-6.045   !(T<=-64 deg)
      B1_fast  =-0.008   
      B21_fast =-0.042   !(T>-64 deg)
      B22_fast =-0.112   !(T<=-64 deg)
      C1_fast  =0.0739   
      C2_fast  =1.2372

      A1_slow  =-0.3949
      A2_slow  =1.282
      B1_slow  =-0.0156
      B2_slow  =0.0111
      B3_slow  =0.0217
      C1_slow  =0.120
      C2_slow  =2.312

      Ni = 0.0

!----------------------------
!RHw parameters
      A = 6.0e-4*log(ww)+6.6e-3
      B = 6.0e-2*log(ww)+1.052
      C = 1.68  *log(ww)+129.35
      RHw=(A*T*T+B*T+C)*0.01

      !czhao turn off the RH criteria
!     if((T.le.-37.0) .and. ((RH*subgrid).ge.RHw)) then
      if(T.le.-37.0) then

        regm = 6.07*log(ww)-55.0

        if(T.ge.regm) then    ! fast-growth regime

          if(T.gt.-64.0) then
            A2_fast=A21_fast
            B2_fast=B21_fast
          else
            A2_fast=A22_fast
            B2_fast=B22_fast
          endif

          k1_fast = exp(A2_fast + B2_fast*T + C2_fast*log(ww))
          k2_fast = A1_fast+B1_fast*T+C1_fast*log(ww)

          Ni = k1_fast*Na**(k2_fast)
          Ni = min(Ni,Na)

        else       ! slow-growth regime

          k1_slow = exp(A2_slow + (B2_slow+B3_slow*log(ww))*T + C2_slow*log(ww))
          k2_slow = A1_slow+B1_slow*T+C1_slow*log(ww)

          Ni = k1_slow*Na**(k2_slow)
          Ni = min(Ni,Na)

        endif

      end if

      return
  end subroutine hf
!----------------------------------------------------------------------
!----------------------------------------------------------------------
  subroutine hetero(T,ww,Ns,Nis,Nid)

    real, intent(in)  :: T, ww, Ns
    real, intent(out) :: Nis, Nid

    real A11,A12,A21,A22,B11,B12,B21,B22
    real A,B,C

!---------------------------------------------------------------------
! parameters

      A11 = 0.0263
      A12 = -0.0185
      A21 = 2.758
      A22 = 1.3221
      B11 = -0.008
      B12 = -0.0468
      B21 = -0.2667
      B22 = -1.4588

!     ice from immersion nucleation (cm^-3)

      B = (A11+B11*log(Ns)) * log(ww) + (A12+B12*log(Ns))
      C =  A21+B21*log(Ns)

      Nis = exp(A22) * Ns**B22 * exp(B*T) * ww**C
      Nis = min(Nis,Ns)

      Nid = 0.0    ! don't include deposition nucleation for cirrus clouds when T<-37C

      return
  end subroutine hetero
!----------------------------------------------------------------------
  subroutine nucleatiphi(tc2, P_SAT, n_dust, A_dust, n_soot, A_soot, & 
  num_organic, A_organic_phobic, num_ic_dust_imm, num_ic_soot_imm, & 
  num_ic_organic_imm)

  real, intent(in)  :: tc2                               ! temperature in K
  real, intent(in)  :: P_SAT                             ! pressure in Pa
  real, intent(in)  :: n_dust, n_soot, num_organic       ! #/kg
  real, intent(in)  :: A_dust, A_soot, A_organic_phobic  ! m2/kg
  real, intent(out) :: num_ic_dust_imm, num_ic_soot_imm, num_ic_organic_imm  ! #/kg

!
! Local workspace
!
  INTEGER ::  IMM_COND_ON, SIZE_RESOLVED
  real :: BASE_DUST_OMEGA, BASE_SOOT_PHILIC_OMEGA, BASE_ORGANIC_PHOBIC_OMEGA, &
              FRACTION_DEPNUCL_WARM_DUST, FAC_CORRECT_RH,  &
              RHO_CFDC, ALPHA_DUST, ALPHA_SOOT, ALPHA_ORGANIC, PIE, EPS

  INTEGER ::  ij
  real :: n_iw, nin_a_nuc_dust, nin_a_nuc_soot, nin_a_nuc_organic
  real :: tc1, ESW, ESI, QSI, QSW, rho, qcb
  real :: Psi_dust, Psi_soot, Psi_organic
  real :: SS_iw, SS_i, tc_HM_degC, S_i_0, S_w_0, dep_frac, &
          H_frac_dust, H_frac_soot, H_frac_organic, tc2_dum
  real :: n_in, n_in_max, n_in_hat, n_in_tilde,  &
              n_in_dust, mu, D_mean_IN, sum_n,  &
              n_in_soot, n_in_organic, n_in_ultra, n_in_dust_ultra, n_in_soot_ultra, n_in_organic_ultra
  real :: CIHENC_dust, CIHENC_soot, CIHENC_organic

  PARAMETER(IMM_COND_ON = 1, SIZE_RESOLVED = 0)
  PARAMETER(BASE_DUST_OMEGA=2.0e-6,BASE_SOOT_PHILIC_OMEGA=2.7e-7,BASE_ORGANIC_PHOBIC_OMEGA=9.1e-7)
  PARAMETER(FRACTION_DEPNUCL_WARM_DUST=0.15,FAC_CORRECT_RH=2.,RHO_CFDC = 50000./(287.*228.15))
  PARAMETER(ALPHA_DUST=2./3., ALPHA_SOOT = 1./3. - 0.06, ALPHA_ORGANIC = 0.06, PIE=3.1415926, EPS = 0.622)

  n_iw = 0.
  nin_a_nuc_dust = 0.; nin_a_nuc_soot = 0.; nin_a_nuc_organic = 0.;
  num_ic_dust_imm = 0.; num_ic_soot_imm = 0.; num_ic_organic_imm = 0.;

  tc1 = 273.
  ESI = polysvp1(tc2,1)       ! Pa
  QSI = EPS*ESI/(P_SAT-ESI)
  ESW = polysvp1(tc2,0)       ! Pa
  QSW = EPS*ESW/(P_SAT-ESW)
  rho = P_SAT/(287.*tc2)

  if(IMM_COND_ON == 1) then
        qcb = QSW
  else
        qcb = QSI +  (QSW - QSI)*2./3.
        if(qcb > QSW*0.97) qcb = QSW*0.97
        if(qcb < QSI) qcb = QSI

  endif

                Psi_dust = A_dust/BASE_DUST_OMEGA
                Psi_soot = A_soot/BASE_SOOT_PHILIC_OMEGA
                Psi_organic = A_organic_phobic/BASE_ORGANIC_PHOBIC_OMEGA

                if(qcb >  QSI) then
                        if(tc2 < 273.15 .and. tc2 > 273.15 - 90. ) then
                                SS_iw = QSW/QSI - 1.

                                if(qcb < QSW) then
                                        SS_i = qcb/QSI - 1.
                                else
                                        SS_i = QSW/QSI - 1.
                                endif

                                tc_HM_degC = tc2 - 273.15

                                S_i_0 = 1. + 10.**(8.2584e-6*tc_HM_degC*tc_HM_degC*tc_HM_degC + 5.3938E-4*tc_HM_degC*tc_HM_degC &
                                         + 3.1656E-3*tc_HM_degC - 1.0261)

                                S_w_0 = 0.97

                                dep_frac = H_1_smooth(SS_i + 1., S_i_0,  S_i_0 + 0.1, 0.,1.)* &
                                H_1_smooth(-(tc2-273.15), 35., 40., FRACTION_DEPNUCL_WARM_DUST, 1.)/FAC_CORRECT_RH

                                H_frac_dust = dep_frac  + (1. - dep_frac)*H_1_smooth(qcb/QSW, S_w_0, 1., 0.,1.)
                                if(H_frac_dust > 1.) H_frac_dust = 1.

                                S_w_0 = 0.97
                                tc2_dum = tc2
                                if(tc2_dum > 239.) tc2_dum = 239.
                                if(tc2_dum < 198.) tc2_dum = 198.

                                S_i_0 = -3.118E-5*tc2_dum*tc2_dum + 1.085e-2*tc2_dum + 5.652e-1

                                dep_frac = H_1_smooth(SS_i + 1., S_i_0-0.1, S_i_0, 0.,1.)* H_1_smooth(-(tc2-273.15), 45., 50., 0.,1.) &
                                      /FAC_CORRECT_RH
                                H_frac_soot = dep_frac  + (1. - dep_frac)*H_1_smooth(qcb/QSW, S_w_0, 1., 0.,1.)
                                if(H_frac_soot > 1.) H_frac_soot = 1.

                                H_frac_organic = H_frac_soot

                                if(tc2 < 273.15 .and. tc2 >= 273.15 - 35.) then
                                        n_in = 1.E3* (exp(12.96*SS_i - 0.639)/RHO_CFDC) *0.0587*FAC_CORRECT_RH
                                        if( tc2 > 273.15 -5. .and. tc2 < 273.15 - 2. ) then
                                                n_in = n_in*H_1_smooth(-(tc2-273.15), 2., 5., 0., 1.)
                                        endif
                                        if(tc2 >= 273.15 - 2. ) n_in = 0.

                                        if(tc2 < 273.15 -25. ) then
                                                n_in_tilde = 1000.*(exp(0.1296*(SS_i*100.-10.))**0.3)*FAC_CORRECT_RH/RHO_CFDC  ! PNAS (2003) data from Paul
                                                n_in_hat = n_in

                                                if(tc2 >= 273.15 - 30.) n_in_max = 1.E3* (exp(12.96*SS_iw - 0.639)/RHO_CFDC) *0.0587*FAC_CORRECT_RH
                                                if(tc2 < 273.15 - 30.) n_in_max = 1000.*(exp(0.1296*(SS_iw*100.-10.))**0.3)*FAC_CORRECT_RH/RHO_CFDC

                                                if(n_in_hat > n_in_max) n_in_hat = n_in_max
                                                if(n_in_tilde > n_in_max) n_in_tilde = n_in_max

                                                n_in = n_in_hat * ((n_in_tilde/n_in_hat)**(H_1_smooth(-(tc2-273.15), 25., 35., 0., 1.)))

                                                if(n_in > n_in_max) n_in = n_in_max

                                        endif

                                        if(SIZE_RESOLVED == 0) then
                                                n_in_dust = Psi_dust*n_in*ALPHA_DUST*H_frac_dust
                                        endif

                                        if(SIZE_RESOLVED == 0) then
                                                n_in_soot = Psi_soot*n_in*ALPHA_SOOT*H_frac_soot
                                        endif

                                        if(SIZE_RESOLVED == 0) then
                                                n_in_organic = Psi_organic*n_in*ALPHA_ORGANIC*H_frac_organic
                                        endif
                                else
                                        n_in = 0.; n_in_ultra = 0.; n_in_dust = 0.;  n_in_soot  = 0.; n_in_organic = 0.;
                                endif

                                if(tc2 < 273.15 - 35.) then
                                        n_in_ultra = 1000.*(exp(0.1296*(SS_i*100.-10.))**0.3)*FAC_CORRECT_RH/RHO_CFDC  ! PNAS (2003) data from Paul

                                        if(SIZE_RESOLVED == 0) then
                                                n_in_dust_ultra = Psi_dust*n_in_ultra*ALPHA_DUST*H_frac_dust
                                        endif
                                        if(SIZE_RESOLVED == 0) then
                                                n_in_soot_ultra = Psi_soot*n_in_ultra*ALPHA_SOOT*H_frac_soot
                                        endif

                                        if(SIZE_RESOLVED == 0) then
                                                n_in_organic_ultra = Psi_organic*n_in_ultra*ALPHA_ORGANIC*H_frac_organic
                                        endif

                                else
                                        n_in_ultra = 0.; n_in_dust_ultra = 0.; n_in_soot_ultra = 0.;  n_in_organic_ultra = 0.;
                                endif

                                n_in_dust = n_in_dust + n_in_dust_ultra;
                                n_in_soot = n_in_soot + n_in_soot_ultra;
                                n_in_organic = n_in_organic + n_in_organic_ultra;

! PROBLEM:  how to ensure that the frozen fraction does not exceed 1 ?

                                if(n_in_dust + n_in_organic + n_in_soot > 0.) then

                                        CIHENC_dust = n_in_dust - nin_a_nuc_dust
                                        if(CIHENC_dust < 0.) CIHENC_dust = 0.

                                        CIHENC_soot = n_in_soot - nin_a_nuc_soot
                                        if(CIHENC_soot < 0.) CIHENC_soot = 0.

                                        CIHENC_organic = n_in_organic - nin_a_nuc_organic
                                        if(CIHENC_organic < 0.) CIHENC_organic = 0.

                                        n_iw =  n_iw + CIHENC_dust
                                        nin_a_nuc_dust = nin_a_nuc_dust + CIHENC_dust
                                        num_ic_dust_imm = num_ic_dust_imm + CIHENC_dust

                                        n_iw =  n_iw + CIHENC_soot
                                        nin_a_nuc_soot = nin_a_nuc_soot + CIHENC_soot
                                        num_ic_soot_imm = num_ic_soot_imm + CIHENC_soot

                                        n_iw =  n_iw + CIHENC_organic
                                        nin_a_nuc_organic = nin_a_nuc_organic + CIHENC_organic
                                        num_ic_organic_imm = num_ic_organic_imm + CIHENC_organic

                                endif

                endif
                endif

  return
  end subroutine nucleatiphi

real function H_1(X, X_1, X_2, Hlo)
real, intent(in) :: Hlo, X, X_1, X_2

if(X >= X_2) H_1 = 1
if(X <= X_1) H_1 = Hlo
if(X > X_1 .and. X < X_2) H_1 = (X - X_1)/(X_2 - X_1)

if( X_2 <= X_1) stop 91919

return
end function H_1

real function H_1_smooth(X, X_1, X_2, Hlo, Hhi)
real, intent(in) :: Hlo, Hhi, X, X_1, X_2
real :: a_0, a_1, a_2, a_3, A, B

if(X >= X_2) H_1_smooth = Hhi
if(X <= X_1) H_1_smooth = Hlo
if(X > X_1 .and. X < X_2) then
        A = 6.*(Hlo - Hhi)/(X_2**3. - X_1**3. + 3.*(X_2*X_1*X_1 - X_1*X_2*X_2) )
        a_3 = (A/3.)
        a_2 = -(A/2.)*(X_1 + X_2)
        a_1 = A*X_2*X_1
        B = Hlo + A*(X_1**3.)/6. - A*X_1*X_1*X_2/2.
        a_0 = B
        H_1_smooth = a_0 + a_1*X + a_2*X*X + a_3*X*X*X
endif

if( X_2 <= X_1) stop 91919

return
end function H_1_smooth

!----------------------------------------------------------------------
       subroutine vqsatd_water(t ,p ,esat ,qsat  ,gam)

       use module_model_constants, only: cp

!------------------------------Arguments--------------------------------
! Input arguments
!
   real, intent(in) :: t       ! temperature
   real, intent(in) :: p       ! pressure

!
! Output arguments
!
   real, intent(out) :: esat   ! saturation vapor pressure
   real, intent(out) :: qsat   ! saturation specific humidity
   real, intent(out) :: gam  ! (l/cp)*(d(qs)/dt)

!--------------------------Local Variables------------------------------
!
   integer i      ! index for vector calculations
!
   real omeps     ! 1. - 0.622
   real hlatv     ! Latent heat of fusion (J/Kg)
   real hltalt    ! appropriately modified hlat for T derivatives
   real hlatsb    ! hlat weighted in transition region
   real hlatvp    ! hlat modified for t changes above freezing
   real desdt     ! d(es)/dT
   real epsqs     ! Ratio of h2o to dry air molecular weights
   real tmelt     ! freezing T of fresh water (K)
   real rgasv     ! Gas constant for water vapor

   CHARACTER (LEN=200) :: message

!
!-----------------------------------------------------------------------
!
   rgasv = 6.02214e26*1.38065e-23/18.016
   tmelt = 273.15
   hlatv = 2.501e6
   epsqs = 18.016/28.966 
   omeps = 1.0 - epsqs

   esat = polysvp1(t,0)
!
! Saturation specific humidity
!
   qsat = epsqs*esat/(p - omeps*esat)
   if (qsat.le.0) then 
      write(message,'(A,1x,4(1x,e10.3))') 'czhao check qsat ',qsat,esat,p,t
      CALL wrf_message( TRIM(message) )
   endif

!
! The following check is to avoid the generation of negative
! values that can occur in the upper stratosphere and mesosphere
!
   qsat = min(1.0,qsat)
!
   if (qsat < 0.0) then
     qsat = 1.0
     esat = p
   end if
!
! No icephs or water to ice transition
!
! Account for change of hlatv with t above freezing where
! constant slope is given by -2369 j/(kg c) = cpv - cw
!
      hlatvp = hlatv - 2369.0*(t-tmelt)
      hlatsb = hlatv
      if (t < tmelt) then
         hltalt = hlatsb
      else
         hltalt = hlatvp
      end if
      desdt  = hltalt*esat/(rgasv*t*t)
      gam = hltalt*qsat*p*desdt/(cp*esat*(p - omeps*esat))
      if (qsat == 1.0) gam = 0.0
!
   return
!
        end subroutine vqsatd_water
      function polysvp1 (T,type)
!  Compute saturation vapor pressure by using
! function from Goff and Gatch (1946)

!  Polysvp returned in units of pa.
!  T is input in units of K.
!  type refers to saturation with respect to liquid (0) or ice (1)

      real dum
      real T,polysvp1
      integer type

! ice
      if (type.eq.1) then
! Goff Gatch equation (good down to -100 C)
         polysvp1 = 10.**(-9.09718*(273.16/t-1.)-3.56654* &
          log10(273.16/t)+0.876793*(1.-t/273.16)+ &
          log10(6.1071))*100.
      end if

! Goff Gatch equation, uncertain below -70 C

      if (type.eq.0) then
         polysvp1 = 10.**(-7.90298*(373.16/t-1.)+ &
             5.02808*log10(373.16/t)- &
             1.3816e-7*(10.**(11.344*(1.-t/373.16))-1.)+ &
             8.1328e-3*(10.**(-3.49149*(373.16/t-1.))-1.)+ &
             log10(1013.246))*100.
      endif


      end function polysvp1


END MODULE module_Cmixactivate
