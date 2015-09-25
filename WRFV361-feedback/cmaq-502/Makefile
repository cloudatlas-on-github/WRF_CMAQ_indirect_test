#   Makefile generated using program [bldmake]
#
#   Generation date [09/09/2014  09:12:43]
#   Configuration file [cfg.CCTM_v5.1_UT.bld]
#   Using GIT repository
#

MODEL = CCTM_v5.1_UT

FC = /usr/local/apps/intel/composer_xe_2013.2.146/bin/intel64/ifort
CC = /usr/bin/cc

f_FLAGS    = -fixed -132 -O3 -override-limits -fno-alias -mp1 -fp-model precise -I /home/wdx/lib_sol/x86_64/intel/ioapi_3.1/Linux2_x86_64ifort -I /home/wdx/lib_sol/x86_64/intel/mpich/include -I /home/wdx/lib_sol/x86_64/intel/ioapi_3.1/ioapi/fixed_src -I.
F_FLAGS    = $(f_FLAGS)
f90_FLAGS  = -free -O3 -fno-alias -mp1 -fp-model precise -I /home/wdx/lib_sol/x86_64/intel/ioapi_3.1/Linux2_x86_64ifort -I /home/wdx/lib_sol/x86_64/intel/mpich/include -I /home/wdx/lib_sol/x86_64/intel/ioapi_3.1/ioapi/fixed_src -I.
F90_FLAGS  = $(f90_FLAGS)
C_FLAGS    = -O2 -DFLDMN -I /home/wdx/lib_sol/x86_64/intel/mpich/include

LINKER     = $(FC)
LINK_FLAGS = -i-static

REPOROOT = /home/wdx/CCTM

CPP = $(FC)
CPP_FLAGS = \
  -Dparallel \
  -DSUBST_BARRIER=SE_BARRIER \
  -DSUBST_GLOBAL_MAX=SE_GLOBAL_MAX \
  -DSUBST_GLOBAL_MIN=SE_GLOBAL_MIN \
  -DSUBST_GLOBAL_MIN_DATA=SE_GLOBAL_MIN_DATA \
  -DSUBST_GLOBAL_TO_LOCAL_COORD=SE_GLOBAL_TO_LOCAL_COORD \
  -DSUBST_GLOBAL_SUM=SE_GLOBAL_SUM \
  -DSUBST_GLOBAL_LOGICAL=SE_GLOBAL_LOGICAL \
  -DSUBST_LOOP_INDEX=SE_LOOP_INDEX \
  -DSUBST_SUBGRID_INDEX=SE_SUBGRID_INDEX \
  -DSUBST_HI_LO_BND_PE=SE_HI_LO_BND_PE \
  -DSUBST_SUM_CHK=SE_SUM_CHK \
  -DSUBST_INIT_ARRAY=SE_INIT_ARRAY \
  -DSUBST_COMM=SE_COMM \
  -DSUBST_MY_REGION=SE_MY_REGION \
  -DSUBST_SLICE=SE_SLICE \
  -DSUBST_GATHER=SE_GATHER \
  -DSUBST_DATA_COPY=SE_DATA_COPY \
  -DSUBST_IN_SYN=SE_IN_SYN

IOAPI =  -L/home/wdx/lib_sol/x86_64/intel/ioapi_3.1/Linux2_x86_64ifort -lioapi
NETCDF =  -L/home/wdx/lib_sol/x86_64/intel/netcdf/lib -lnetcdf -lnetcdff
MPICH =  -L/home/wdx/lib_sol/x86_64/intel/mpich/lib -lmpich -lmpiif

LIBRARIES = $(IOAPI) $(NETCDF) $(MPICH)

BASE_INC = .
PA_INC = .

INCLUDES =  \
  -DSUBST_PE_COMM=\"./PE_COMM.EXT\" \
  -DSUBST_CONST=\"./CONST.EXT\" \
  -DSUBST_FILES_ID=\"./FILES_CTM.EXT\" \
  -DSUBST_EMISPRM=\"./EMISPRM.EXT\" \
  -DSUBST_PACTL_ID=\"./PA_CTL.EXT\" \
  -DSUBST_PACMN_ID=\"./PA_CMN.EXT\" \
  -DSUBST_PADAT_ID=\"./PA_DAT.EXT\" \
  -DSUBST_MPI=\"mpif.h\"

GLOBAL_MODULES = \
  se_bndy_copy_info_ext.o \
  se_pe_info_ext.o \
  se_comm_info_ext.o \
  se_data_send_module.o \
  se_domain_info_ext.o \
  se_internal_util_module.o \
  se_data_recv_module.o \
  se_disp_info_ext.o \
  se_subgrid_info_ext.o \
  se_global_sum_module.o \
  se_global_max_module.o \
  se_global_min_module.o \
  se_slice_module.o \
  se_init_module.o \
  se_bndy_copy_module.o \
  se_comm_module.o \
  se_reconfig_grid_info_ext.o \
  se_reconfig_grid_module.o \
  se_gather_module.o \
  se_data_copy_module.o \
  se_term_module.o \
  se_twoway_comm_module.o \
  se_util_module.o \
  se_modules.o \
  swap_sandia.o \
  pinterpb_mod.o \
  parutilio.o \
  twoway_data_module.o \
  UTILIO_DEFN.o \
  VGRD_DEFN.o \
  HGRD_DEFN.o \
  GRID_CONF.o \
  twoway_sd_time_series.o \
  PAGRD_DEFN.o \
  RXNS_DATA_MODULE.o \
  CGRID_SPCS.o \
  BIDI_MOD.o \
  LSM_MOD.o \
  DEPVVARS.o \
  MOSAIC_MOD.o \
  ABFLUX_MOD.o \
  AEROMET_DATA.o \
  BIOG_EMIS.o \
  BEIS_DEFN.o \
  AERO_DATA.o \
  LTNG_DEFN.o \
  UDTYPES.o \
  STK_EMIS.o \
  STK_PRMS.o \
  PTMAP.o \
  PTBILIN.o \
  SSEMIS.o \
  PTMET.o \
  PT3D_DEFN.o \
  SOA_DEFN.o \
  PCGRID_DEFN.o \
  STD_CONC.o \
  RXNS_FUNC_MODULE.o \
  DUST_EMIS.o \
  PRECURSOR_DATA.o \
  AERO_EMIS.o \
  AEROSOL_CHEMISTRY.o \
  DEPV_DEFN.o \
  EMIS_DEFN.o \
  WVEL_DEFN.o

STENEX = \
  swap_sandia_routines.o

PARIO = \
  alloc_data_mod.o \
  piomaps_mod.o \
  boundary.o \
  get_write_map.o \
  growbuf.o \
  gtndxhdv.o \
  interpol.o \
  pinterpb.o \
  pio_init.o \
  pio_re_init.o \
  pm3err.o \
  pm3exit.o \
  pm3warn.o \
  pshut3.o \
  ptrwrite3.o \
  pwrgrdd.o \
  pwrite3.o \
  readbndy.o \
  subdmap.o \
  wrsubmap.o

PAR = \
  distr_env.o \
  mpcomm_init.o

TWOWAY = \
  twoway_cgrid_aerosol_spc_map_module.o \
  twoway_header_data_module.o \
  twoway_met_param_module.o \
  twoway_util_module.o \
  twoway_aqprep.o \
  twoway_aqprep_util.o \
  twoway_feedback.o

DRIVER = \
  AVG_CONC.o \
  advstep.o \
  driver.o \
  hveloc.o \
  sciproc.o \
  wr_aconc.o \
  wr_cgrid.o \
  wr_conc.o

INIT = \
  diffmsg.o \
  flcheck.o \
  grdcheck.o \
  initscen.o \
  load_cgrid.o \
  opconc.o

COUPLE = \
  couple.o \
  couplewrf.o

HADV = \
  advbc_map.o \
  hadvyppm.o \
  hcontvel.o \
  hppm.o \
  rdbcon.o \
  x_ppm.o \
  x_yamo.o \
  y_ppm.o \
  y_yamo.o \
  zfdbc.o

VADV = \
  vppm.o \
  zadvppmwrf.o

HDIFF = \
  deform.o \
  hcdiff3d.o \
  hdiff.o \
  rho_j.o

VDIFF = \
  VDIFF_DIAG.o \
  VDIFF_MAP.o \
  VDIFF_MET.o \
  aero_sedv.o \
  conv_cgrid.o \
  eddyx.o \
  matrix1.o \
  opddep.o \
  opddep_fst.o \
  opddep_mos.o \
  rddepv.o \
  sedi.o \
  tri.o \
  vdiffacmx.o \
  vdiffproc.o

DEPV = \
  cgrid_depv.o \
  gas_depv_map.o \
  m3dry.o \
  opdepv_diag.o \
  opdepv_fst.o \
  opdepv_mos.o

EMIS = \
  cropcal.o \
  opemis.o \
  tfabove.o \
  tfbelow.o

BIOG = \
  beis3.o \
  checkmem.o \
  chkgrid.o \
  czangle.o \
  getfline.o \
  getparb.o \
  hrno.o \
  parsline.o \
  tmpbeis.o \
  wrdaymsg.o

PLRISE = \
  delta_zs.o \
  fire_plmris.o \
  openlayout.o \
  oppt3d_diag.o \
  plmris.o \
  plsprd.o \
  preplm.o \
  ungridb2.o \
  write3_distr.o

PHOT = \
  CSQY_DATA.o \
  PHOT_MOD.o \
  aero_photdata.o \
  o3totcol.o \
  opphot.o \
  phot.o

GAS = \
  hrdata_mod.o \
  hrdriver.o \
  hrg1.o \
  hrg2.o \
  hrg3.o \
  hrg4.o \
  hrinit.o \
  hrprodloss.o \
  hrrates.o \
  hrsolver.o

AERO = \
  aero_depv.o \
  aero_driver.o \
  aero_subs.o \
  coags.o \
  getpar.o \
  isocom.o \
  isofwd.o \
  isorev.o \
  opdiam.o \
  opvis.o

CLOUD = \
  AQ_DATA.o \
  acmcld.o \
  aqchem.o \
  aq_map.o \
  cldproc_acm.o \
  convcld_acm.o \
  getalpha.o \
  hlconst.o \
  indexn.o \
  rescld.o \
  scavwdep.o

PROCAN = \
  pa_init.o \
  pa_irr_ctl.o \
  pa_irr.o \
  pa_mkhdr.o \
  pa_output.o \
  pa_update.o

UTIL = \
  cksummer.o \
  findex.o \
  get_envlist.o \
  lstepf.o \
  setup_logdev.o \
  subhdomain.o \
  subhfile.o

OBJS = \
  $(GLOBAL_MODULES) \
  $(STENEX) \
  $(PARIO) \
  $(PAR) \
  $(TWOWAY) \
  $(DRIVER) \
  $(INIT) \
  $(COUPLE) \
  $(HADV) \
  $(VADV) \
  $(HDIFF) \
  $(VDIFF) \
  $(DEPV) \
  $(EMIS) \
  $(BIOG) \
  $(PLRISE) \
  $(PHOT) \
  $(GAS) \
  $(AERO) \
  $(CLOUD) \
  $(PROCAN) \
  $(UTIL)

.SUFFIXES: .F .f .c .F90 .f90

$(MODEL): $(OBJS)
	$(LINKER) $(LINK_FLAGS) $(OBJS) $(LIBRARIES) -o $@

.F.o:
	$(FC) -c $(F_FLAGS) $(CPP_FLAGS) $(INCLUDES) $<

.f.o:
	$(FC) -c $(f_FLAGS) $<

.F90.o:
	$(FC) -c $(F90_FLAGS) $(CPP_FLAGS) $(INCLUDES) $<

.f90.o:
	$(FC) -c $(f90_FLAGS) $<

.c.o:
	$(CC) -c $(C_FLAGS) $<

clean:
	rm -f $(OBJS) *.mod


# dependencies

se_bndy_copy_module.o:	se_bndy_copy_info_ext.f se_internal_util_module.f \
		se_data_send_module.f se_data_recv_module.f se_pe_info_ext.f \
		se_domain_info_ext.f
se_comm_module.o:	se_comm_info_ext.f se_internal_util_module.f \
		se_data_send_module.f se_data_recv_module.f se_pe_info_ext.f \
		se_domain_info_ext.f
se_data_copy_module.o:	se_subgrid_info_ext.f se_pe_info_ext.f \
		se_data_send_module.f se_data_recv_module.f se_domain_info_ext.f
se_data_recv_module.o:	se_pe_info_ext.f
se_data_send_module.o:	se_pe_info_ext.f
se_gather_module.o:	se_pe_info_ext.f se_slice_module.f se_global_sum_module.f
se_global_max_module.o:	se_pe_info_ext.f
se_global_min_module.o:	se_pe_info_ext.f
se_global_sum_module.o:	se_pe_info_ext.f
se_init_module.o:	se_pe_info_ext.f se_domain_info_ext.f se_comm_info_ext.f \
		se_disp_info_ext.f se_bndy_copy_info_ext.f se_internal_util_module.f
se_internal_util_module.o:	se_pe_info_ext.f
se_modules.o:	se_init_module.F se_term_module.f se_util_module.f \
		se_comm_module.f se_slice_module.f se_data_copy_module.f se_gather_module.f \
		se_reconfig_grid_module.f se_bndy_copy_module.f se_global_max_module.f \
		se_global_min_module.f se_global_sum_module.f se_twoway_comm_module.f
se_reconfig_grid_module.o:	se_reconfig_grid_info_ext.f se_pe_info_ext.f \
		se_data_send_module.f se_data_recv_module.f se_domain_info_ext.f
se_slice_module.o:	se_pe_info_ext.f
se_term_module.o:	se_pe_info_ext.f se_domain_info_ext.f se_reconfig_grid_info_ext.f \
		se_subgrid_info_ext.f
se_twoway_comm_module.o:	se_comm_info_ext.f
se_util_module.o:	se_comm_info_ext.f se_domain_info_ext.f se_pe_info_ext.f \
		se_subgrid_info_ext.f se_internal_util_module.f
swap_sandia_routines.o:	swap_sandia.f
get_write_map.o:	piomaps_mod.f
growbuf.o:	pinterpb_mod.f
gtndxhdv.o:	pinterpb_mod.f
parutilio.o:	pinterpb_mod.f
pinterpb.o:	pinterpb_mod.f
pio_init.o:	piomaps_mod.f
pio_re_init.o:	piomaps_mod.f
ptrwrite3.o:	piomaps_mod.f
pwrgrdd.o:	piomaps_mod.f alloc_data_mod.f
pwrite3.o:	piomaps_mod.f
readbndy.o:	piomaps_mod.f pinterpb_mod.f
mpcomm_init.o:	UTILIO_DEFN.F
twoway_aqprep.o:	twoway_util_module.F90 twoway_header_data_module.F90 \
		twoway_met_param_module.F90 twoway_data_module.F90 HGRD_DEFN.F \
		se_comm_info_ext.f UTILIO_DEFN.F ./CONST.EXT
twoway_aqprep_util.o:	./CONST.EXT
twoway_cgrid_aerosol_spc_map_module.o:	AERO_DATA.F
twoway_feedback.o:	twoway_header_data_module.F90 twoway_met_param_module.F90 \
		twoway_data_module.F90 twoway_util_module.F90 twoway_cgrid_aerosol_spc_map_module.F90 \
		UTILIO_DEFN.F HGRD_DEFN.F AERO_DATA.F CGRID_SPCS.F
twoway_sd_time_series.o:	HGRD_DEFN.F UTILIO_DEFN.F ./FILES_CTM.EXT
twoway_util_module.o:	twoway_header_data_module.F90 UTILIO_DEFN.F
advstep.o:	GRID_CONF.F UTILIO_DEFN.F se_modules.f ./CONST.EXT \
		./FILES_CTM.EXT
AVG_CONC.o:	VGRD_DEFN.F CGRID_SPCS.F UTILIO_DEFN.F
driver.o:	PCGRID_DEFN.F RXNS_DATA_MODULE.F90 CGRID_SPCS.F STD_CONC.F \
		AVG_CONC.F WVEL_DEFN.F PAGRD_DEFN.F UTILIO_DEFN.F twoway_sd_time_series.F90 \
		se_modules.f swap_sandia.f ./FILES_CTM.EXT ./PA_CTL.EXT
hveloc.o:	GRID_CONF.F UTILIO_DEFN.F se_modules.f ./PE_COMM.EXT \
		./CONST.EXT ./FILES_CTM.EXT
sciproc.o:	UTILIO_DEFN.F ./PA_CTL.EXT
STD_CONC.o:	HGRD_DEFN.F VGRD_DEFN.F CGRID_SPCS.F UTILIO_DEFN.F
wr_aconc.o:	GRID_CONF.F AVG_CONC.F UTILIO_DEFN.F se_modules.f \
		./FILES_CTM.EXT
wr_cgrid.o:	GRID_CONF.F CGRID_SPCS.F UTILIO_DEFN.F se_modules.f \
		./FILES_CTM.EXT
wr_conc.o:	GRID_CONF.F STD_CONC.F WVEL_DEFN.F UTILIO_DEFN.F ./FILES_CTM.EXT
WVEL_DEFN.o:	GRID_CONF.F UTILIO_DEFN.F
GRID_CONF.o:	HGRD_DEFN.F VGRD_DEFN.F
HGRD_DEFN.o:	UTILIO_DEFN.F
PAGRD_DEFN.o:	GRID_CONF.F se_modules.f ./PA_CTL.EXT ./PA_CMN.EXT
PCGRID_DEFN.o:	GRID_CONF.F CGRID_SPCS.F UTILIO_DEFN.F
VGRD_DEFN.o:	UTILIO_DEFN.F
diffmsg.o:	UTILIO_DEFN.F
flcheck.o:	UTILIO_DEFN.F ./FILES_CTM.EXT
grdcheck.o:	GRID_CONF.F UTILIO_DEFN.F ./FILES_CTM.EXT
initscen.o:	HGRD_DEFN.F CGRID_SPCS.F UTILIO_DEFN.F se_modules.f \
		./FILES_CTM.EXT
load_cgrid.o:	GRID_CONF.F CGRID_SPCS.F UTILIO_DEFN.F ./CONST.EXT \
		./FILES_CTM.EXT
opconc.o:	GRID_CONF.F WVEL_DEFN.F STD_CONC.F UTILIO_DEFN.F se_modules.f \
		./FILES_CTM.EXT
couple.o:	GRID_CONF.F CGRID_SPCS.F UTILIO_DEFN.F ./FILES_CTM.EXT
couplewrf.o:	GRID_CONF.F CGRID_SPCS.F UTILIO_DEFN.F ./FILES_CTM.EXT
advbc_map.o:	CGRID_SPCS.F UTILIO_DEFN.F ./FILES_CTM.EXT
hadvyppm.o:	GRID_CONF.F CGRID_SPCS.F UTILIO_DEFN.F se_modules.f \
		HGRD_DEFN.F ./PE_COMM.EXT ./FILES_CTM.EXT
hcontvel.o:	GRID_CONF.F UTILIO_DEFN.F se_modules.f ./PE_COMM.EXT \
		./FILES_CTM.EXT
hppm.o:	HGRD_DEFN.F UTILIO_DEFN.F se_modules.f
rdbcon.o:	GRID_CONF.F CGRID_SPCS.F UTILIO_DEFN.F ./CONST.EXT ./FILES_CTM.EXT
x_ppm.o:	HGRD_DEFN.F CGRID_SPCS.F UTILIO_DEFN.F se_modules.f swap_sandia.f \
		./PE_COMM.EXT ./CONST.EXT
x_yamo.o:	HGRD_DEFN.F CGRID_SPCS.F UTILIO_DEFN.F se_modules.f \
		swap_sandia.f ./PE_COMM.EXT ./CONST.EXT
y_ppm.o:	HGRD_DEFN.F CGRID_SPCS.F UTILIO_DEFN.F se_modules.f swap_sandia.f \
		./PE_COMM.EXT ./CONST.EXT
y_yamo.o:	HGRD_DEFN.F CGRID_SPCS.F UTILIO_DEFN.F se_modules.f \
		swap_sandia.f ./PE_COMM.EXT ./CONST.EXT
vppm.o:	CGRID_SPCS.F UTILIO_DEFN.F
zadvppmwrf.o:	GRID_CONF.F CGRID_SPCS.F WVEL_DEFN.F se_modules.f \
		UTILIO_DEFN.F ./FILES_CTM.EXT
deform.o:	GRID_CONF.F UTILIO_DEFN.F se_modules.f swap_sandia.f \
		./PE_COMM.EXT ./CONST.EXT ./FILES_CTM.EXT
hcdiff3d.o:	GRID_CONF.F UTILIO_DEFN.F se_modules.f ./PE_COMM.EXT \
		./CONST.EXT ./FILES_CTM.EXT
hdiff.o:	GRID_CONF.F CGRID_SPCS.F UTILIO_DEFN.F se_modules.f swap_sandia.f \
		./PE_COMM.EXT ./CONST.EXT
rho_j.o:	GRID_CONF.F UTILIO_DEFN.F ./FILES_CTM.EXT
aero_sedv.o:	GRID_CONF.F CGRID_SPCS.F UTILIO_DEFN.F VDIFF_MET.F \
		AERO_DATA.F SOA_DEFN.F AEROMET_DATA.F
conv_cgrid.o:	GRID_CONF.F CGRID_SPCS.F UTILIO_DEFN.F ./CONST.EXT \
		./FILES_CTM.EXT
eddyx.o:	GRID_CONF.F VDIFF_MET.F UTILIO_DEFN.F se_modules.f ./CONST.EXT
matrix1.o:	VGRD_DEFN.F CGRID_SPCS.F UTILIO_DEFN.F
opddep.o:	GRID_CONF.F CGRID_SPCS.F UTILIO_DEFN.F ./FILES_CTM.EXT
opddep_fst.o:	GRID_CONF.F CGRID_SPCS.F UTILIO_DEFN.F LSM_MOD.F \
		./FILES_CTM.EXT
opddep_mos.o:	GRID_CONF.F CGRID_SPCS.F UTILIO_DEFN.F LSM_MOD.F \
		./FILES_CTM.EXT
rddepv.o:	HGRD_DEFN.F RXNS_DATA_MODULE.F90 CGRID_SPCS.F UTILIO_DEFN.F \
		./FILES_CTM.EXT
sedi.o:	CGRID_SPCS.F GRID_CONF.F VDIFF_MET.F VDIFF_MAP.F UTILIO_DEFN.F \
		VDIFF_DIAG.F AERO_DATA.F SOA_DEFN.F AEROMET_DATA.F ./CONST.EXT
tri.o:	VGRD_DEFN.F CGRID_SPCS.F
vdiffacmx.o:	CGRID_SPCS.F GRID_CONF.F EMIS_DEFN.F DEPV_DEFN.F \
		VDIFF_MET.F VDIFF_MAP.F UTILIO_DEFN.F BIDI_MOD.F LSM_MOD.F VDIFF_DIAG.F \
		./CONST.EXT ./FILES_CTM.EXT ./PA_CTL.EXT
VDIFF_DIAG.o:	UTILIO_DEFN.F GRID_CONF.F CGRID_SPCS.F ./FILES_CTM.EXT
VDIFF_MAP.o:	CGRID_SPCS.F AERO_EMIS.F UTILIO_DEFN.F ./EMISPRM.EXT
VDIFF_MET.o:	GRID_CONF.F UTILIO_DEFN.F se_modules.f ./PE_COMM.EXT \
		./FILES_CTM.EXT
vdiffproc.o:	CGRID_SPCS.F GRID_CONF.F EMIS_DEFN.F DEPV_DEFN.F \
		VDIFF_MET.F VDIFF_MAP.F UTILIO_DEFN.F BIDI_MOD.F LSM_MOD.F VDIFF_DIAG.F \
		se_modules.f ./CONST.EXT ./FILES_CTM.EXT ./PA_CTL.EXT
ABFLUX_MOD.o:	HGRD_DEFN.F UTILIO_DEFN.F BIDI_MOD.F MOSAIC_MOD.F \
		LSM_MOD.F ./CONST.EXT ./FILES_CTM.EXT
BIDI_MOD.o:	GRID_CONF.F CGRID_SPCS.F UTILIO_DEFN.F ./FILES_CTM.EXT
cgrid_depv.o:	CGRID_SPCS.F HGRD_DEFN.F UTILIO_DEFN.F
DEPV_DEFN.o:	HGRD_DEFN.F RXNS_DATA_MODULE.F90 CGRID_SPCS.F DEPVVARS.F \
		UTILIO_DEFN.F MOSAIC_MOD.F LSM_MOD.F BIDI_MOD.F VGRD_DEFN.F \
		se_modules.f ABFLUX_MOD.F GRID_CONF.F AERO_DATA.F SOA_DEFN.F \
		AEROMET_DATA.F ./CONST.EXT ./FILES_CTM.EXT
DEPVVARS.o:	HGRD_DEFN.F UTILIO_DEFN.F
gas_depv_map.o:	CGRID_SPCS.F DEPVVARS.F UTILIO_DEFN.F
LSM_MOD.o:	HGRD_DEFN.F UTILIO_DEFN.F twoway_data_module.F90 ./FILES_CTM.EXT
m3dry.o:	HGRD_DEFN.F DEPVVARS.F VGRD_DEFN.F UTILIO_DEFN.F se_modules.f \
		ABFLUX_MOD.F LSM_MOD.F MOSAIC_MOD.F BIDI_MOD.F twoway_data_module.F90 \
		./PE_COMM.EXT ./CONST.EXT ./FILES_CTM.EXT
MOSAIC_MOD.o:	HGRD_DEFN.F LSM_MOD.F UTILIO_DEFN.F ./FILES_CTM.EXT
opdepv_diag.o:	GRID_CONF.F DEPVVARS.F UTILIO_DEFN.F ./FILES_CTM.EXT
opdepv_fst.o:	GRID_CONF.F DEPVVARS.F UTILIO_DEFN.F LSM_MOD.F ./FILES_CTM.EXT
opdepv_mos.o:	GRID_CONF.F DEPVVARS.F UTILIO_DEFN.F LSM_MOD.F ./FILES_CTM.EXT
BEIS_DEFN.o:	HGRD_DEFN.F CGRID_SPCS.F BIOG_EMIS.F UTILIO_DEFN.F \
		./CONST.EXT
BIOG_EMIS.o:	UTILIO_DEFN.F
cropcal.o:	HGRD_DEFN.F UTILIO_DEFN.F
DUST_EMIS.o:	AERO_DATA.F HGRD_DEFN.F UTILIO_DEFN.F GRID_CONF.F \
		./CONST.EXT ./FILES_CTM.EXT
EMIS_DEFN.o:	GRID_CONF.F CGRID_SPCS.F BEIS_DEFN.F LTNG_DEFN.F \
		PT3D_DEFN.F UTILIO_DEFN.F AERO_EMIS.F AERO_DATA.F DEPV_DEFN.F \
		BIOG_EMIS.F PTMAP.F ./CONST.EXT ./FILES_CTM.EXT
LTNG_DEFN.o:	GRID_CONF.F CGRID_SPCS.F UTILIO_DEFN.F ./CONST.EXT \
		./FILES_CTM.EXT
opemis.o:	VGRD_DEFN.F CGRID_SPCS.F UTILIO_DEFN.F ./FILES_CTM.EXT
PT3D_DEFN.o:	UDTYPES.F GRID_CONF.F CGRID_SPCS.F STK_PRMS.F STK_EMIS.F \
		PTMET.F PTBILIN.F PTMAP.F UTILIO_DEFN.F VGRD_DEFN.F ./CONST.EXT \
		./FILES_CTM.EXT
PTBILIN.o:	UDTYPES.F HGRD_DEFN.F VGRD_DEFN.F STK_PRMS.F
PTMAP.o:	UDTYPES.F CGRID_SPCS.F STK_EMIS.F AERO_DATA.F UTILIO_DEFN.F
PTMET.o:	UDTYPES.F STK_PRMS.F PTBILIN.F UTILIO_DEFN.F VGRD_DEFN.F \
		./FILES_CTM.EXT
SSEMIS.o:	HGRD_DEFN.F AERO_DATA.F UTILIO_DEFN.F GRID_CONF.F AEROMET_DATA.F \
		./FILES_CTM.EXT
STK_EMIS.o:	UDTYPES.F UTILIO_DEFN.F
STK_PRMS.o:	UDTYPES.F HGRD_DEFN.F UTILIO_DEFN.F
tfabove.o:	HGRD_DEFN.F
tfbelow.o:	HGRD_DEFN.F UTILIO_DEFN.F ./CONST.EXT ./FILES_CTM.EXT
beis3.o:	BIOG_EMIS.F ./CONST.EXT
chkgrid.o:	UTILIO_DEFN.F
czangle.o:	UTILIO_DEFN.F ./CONST.EXT
hrno.o:	HGRD_DEFN.F BIOG_EMIS.F UTILIO_DEFN.F
parsline.o:	UTILIO_DEFN.F
tmpbeis.o:	HGRD_DEFN.F UTILIO_DEFN.F BIOG_EMIS.F
wrdaymsg.o:	UTILIO_DEFN.F
fire_plmris.o:	./CONST.EXT
openlayout.o:	HGRD_DEFN.F UTILIO_DEFN.F
oppt3d_diag.o:	GRID_CONF.F UTILIO_DEFN.F
plmris.o:	./CONST.EXT
preplm.o:	UTILIO_DEFN.F
ungridb2.o:	HGRD_DEFN.F PTBILIN.F
CGRID_SPCS.o:	UTILIO_DEFN.F RXNS_DATA_MODULE.F90
aero_photdata.o:	RXNS_DATA_MODULE.F90 CGRID_SPCS.F UTILIO_DEFN.F \
		AERO_DATA.F SOA_DEFN.F AEROMET_DATA.F
o3totcol.o:	UTILIO_DEFN.F
opphot.o:	GRID_CONF.F UTILIO_DEFN.F RXNS_DATA_MODULE.F90 PHOT_MOD.F \
		./FILES_CTM.EXT
phot.o:	RXNS_DATA_MODULE.F90 CGRID_SPCS.F PCGRID_DEFN.F UTILIO_DEFN.F \
		LSM_MOD.F se_modules.f AERO_DATA.F PHOT_MOD.F SOA_DEFN.F AEROMET_DATA.F \
		./FILES_CTM.EXT
PHOT_MOD.o:	CSQY_DATA.F UTILIO_DEFN.F RXNS_DATA_MODULE.F90 ./CONST.EXT
hrdriver.o:	HGRD_DEFN.F VGRD_DEFN.F CGRID_SPCS.F UTILIO_DEFN.F \
		RXNS_DATA_MODULE.F90 AEROSOL_CHEMISTRY.F RXNS_FUNC_MODULE.F90 \
		hrdata_mod.F ./CONST.EXT ./FILES_CTM.EXT ./EMISPRM.EXT ./PA_CTL.EXT
hrg1.o:	hrdata_mod.F
hrg2.o:	hrdata_mod.F
hrg3.o:	hrdata_mod.F
hrg4.o:	hrdata_mod.F
hrinit.o:	UTILIO_DEFN.F hrdata_mod.F RXNS_DATA_MODULE.F90 RXNS_FUNC_MODULE.F90 \
		./CONST.EXT
hrprodloss.o:	hrdata_mod.F
hrrates.o:	hrdata_mod.F
hrsolver.o:	UTILIO_DEFN.F RXNS_DATA_MODULE.F90 hrdata_mod.F
RXNS_FUNC_MODULE.o:	RXNS_DATA_MODULE.F90 UTILIO_DEFN.F CGRID_SPCS.F
AERO_DATA.o:	CGRID_SPCS.F AEROMET_DATA.F UTILIO_DEFN.F ./FILES_CTM.EXT
aero_depv.o:	GRID_CONF.F RXNS_DATA_MODULE.F90 CGRID_SPCS.F UTILIO_DEFN.F \
		AERO_DATA.F SOA_DEFN.F AEROMET_DATA.F MOSAIC_MOD.F LSM_MOD.F \
		DEPV_DEFN.F ./FILES_CTM.EXT
aero_driver.o:	GRID_CONF.F RXNS_DATA_MODULE.F90 AERO_DATA.F PRECURSOR_DATA.F \
		SOA_DEFN.F AEROMET_DATA.F UTILIO_DEFN.F CGRID_SPCS.F AEROSOL_CHEMISTRY.F \
		./FILES_CTM.EXT ./PA_CTL.EXT
AERO_EMIS.o:	AERO_DATA.F GRID_CONF.F RXNS_DATA_MODULE.F90 CGRID_SPCS.F \
		SSEMIS.F DUST_EMIS.F PRECURSOR_DATA.F UTILIO_DEFN.F PT3D_DEFN.F \
		PTMAP.F ./CONST.EXT ./FILES_CTM.EXT
AEROMET_DATA.o:	./CONST.EXT
AEROSOL_CHEMISTRY.o:	UTILIO_DEFN.F GRID_CONF.F RXNS_DATA_MODULE.F90 \
		AERO_DATA.F PRECURSOR_DATA.F SOA_DEFN.F AEROMET_DATA.F ./CONST.EXT
aero_subs.o:	AERO_DATA.F SOA_DEFN.F AEROMET_DATA.F UTILIO_DEFN.F \
		PRECURSOR_DATA.F ./CONST.EXT
getpar.o:	AERO_DATA.F AEROMET_DATA.F
opdiam.o:	GRID_CONF.F UTILIO_DEFN.F ./FILES_CTM.EXT
opvis.o:	GRID_CONF.F UTILIO_DEFN.F ./FILES_CTM.EXT
PRECURSOR_DATA.o:	CGRID_SPCS.F UTILIO_DEFN.F AEROMET_DATA.F
SOA_DEFN.o:	RXNS_DATA_MODULE.F90 AERO_DATA.F AEROMET_DATA.F CGRID_SPCS.F \
		UTILIO_DEFN.F
aqchem.o:	RXNS_DATA_MODULE.F90 AQ_DATA.F UTILIO_DEFN.F ./CONST.EXT
AQ_DATA.o:	AERO_DATA.F CGRID_SPCS.F UTILIO_DEFN.F ./CONST.EXT
aq_map.o:	CGRID_SPCS.F AQ_DATA.F AERO_DATA.F UTILIO_DEFN.F
cldproc_acm.o:	GRID_CONF.F CGRID_SPCS.F UTILIO_DEFN.F ./CONST.EXT \
		./FILES_CTM.EXT
convcld_acm.o:	GRID_CONF.F CGRID_SPCS.F UTILIO_DEFN.F AQ_DATA.F \
		./CONST.EXT ./FILES_CTM.EXT
getalpha.o:	AERO_DATA.F ./CONST.EXT
hlconst.o:	UTILIO_DEFN.F
rescld.o:	GRID_CONF.F CGRID_SPCS.F UTILIO_DEFN.F ./CONST.EXT ./FILES_CTM.EXT
scavwdep.o:	CGRID_SPCS.F UTILIO_DEFN.F ./CONST.EXT
pa_init.o:	GRID_CONF.F PAGRD_DEFN.F UTILIO_DEFN.F se_modules.f \
		./FILES_CTM.EXT ./PA_CTL.EXT ./PA_CMN.EXT ./PA_DAT.EXT
pa_irr_ctl.o:	GRID_CONF.F CGRID_SPCS.F PAGRD_DEFN.F ./PA_CMN.EXT
pa_irr.o:	GRID_CONF.F RXNS_DATA_MODULE.F90 CGRID_SPCS.F UTILIO_DEFN.F \
		./PA_CMN.EXT
pa_mkhdr.o:	GRID_CONF.F CGRID_SPCS.F UTILIO_DEFN.F PAGRD_DEFN.F \
		./PA_CMN.EXT
pa_output.o:	GRID_CONF.F PAGRD_DEFN.F UTILIO_DEFN.F se_modules.f \
		./FILES_CTM.EXT ./PA_CTL.EXT ./PA_CMN.EXT
pa_update.o:	GRID_CONF.F CGRID_SPCS.F PAGRD_DEFN.F AERO_EMIS.F \
		UTILIO_DEFN.F se_modules.f ./CONST.EXT ./FILES_CTM.EXT ./EMISPRM.EXT \
		./PA_CTL.EXT ./PA_CMN.EXT
cksummer.o:	GRID_CONF.F CGRID_SPCS.F UTILIO_DEFN.F se_modules.f
lstepf.o:	UTILIO_DEFN.F ./FILES_CTM.EXT
setup_logdev.o:	UTILIO_DEFN.F
subhdomain.o:	UTILIO_DEFN.F
subhfile.o:	GRID_CONF.F UTILIO_DEFN.F ./FILES_CTM.EXT
UTILIO_DEFN.o:	parutilio.f
