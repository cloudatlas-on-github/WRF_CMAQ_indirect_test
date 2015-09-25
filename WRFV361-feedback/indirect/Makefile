#

LN      =       ln -sf
MAKE    =       make -i -r
RM      =       rm -f

MODULES = \
        module_data_cmaq.o \
        module_Cmixactivate.o \
        cmaq_mixactivate.o \

LIBTARGET    =  dyn_em
TARGETDIR    =  ./
$(LIBTARGET) :  $(MODULES)
		$(AR) $(ARFLAGS) ../main/$(LIBWRFLIB) $(MODULES) $(OBJS)

include ../configure.wrf

clean:
	@ echo 'use the clean script'

# DEPENDENCIES : only dependencies after this line (don't remove the word DEPENDENCIES)

module_Cmixactivate.o: module_data_cmaq.o \
		../dyn_em/aerosols_cmaq_init_aercld.o

cmaq_mixactivate.o : \
		../phys/module_Cmixactivate.o

# DO NOT DELETE
