#
# This is the make file for the user_lib library.
#
#				Gershon Elber, Aug 1990
#

include ../makeflag.sas

OBJS =  belts.o crvtranl.o crv_arng.o ddm_text.o dtr3d2im.o dtr3d3im.o \
	ff_krnl.o ffnt2bzr.o font3d.o fontlout.o \
	gcsetcvr.o gcvismap.o hrmt_crv.o \
	imgshd3d.o kinematc.o nc_tpath.o plyround.o \
	register.o scalimag.o rldmatch.o srfcrvtr.o srfpgeom.o \
        srf_mrch.o srf_cntr.o srf_ssi.o textwarp.o tv0jacob.o \
	user_ftl.o user_err.o userpick.o \
	usrcnvrt.o visible.o zcollide.o

all:	user.lib

user.lib: $(OBJS)
	rm -f user.lib
	oml user.lib a $(OBJS)

install: user.lib
	mv -f user.lib $(IRIT_LIB_DIR)

# DO NOT DELETE THIS LINE -- make depend depends on it.

srf_cntr.o: ../include/irit_sm.h ../include/cagd_lib.h ../include/miscattr.h
srf_cntr.o: ../include/misc_lib.h ../include/symb_lib.h ../include/iritprsr.h
srf_cntr.o: ../include/trim_lib.h ../include/triv_lib.h ../include/trng_lib.h
srf_cntr.o: ../include/mdl_lib.h ../include/mvar_lib.h ../include/ip_cnvrt.h
srf_cntr.o: ../include/attribut.h ../include/bool_lib.h ../include/allocate.h
srf_cntr.o: ../include/obj_dpnd.h ../include/geom_lib.h user_loc.h
srf_cntr.o: ../include/user_lib.h
srf_mrch.o: ../include/irit_sm.h ../include/cagd_lib.h ../include/miscattr.h
srf_mrch.o: ../include/misc_lib.h ../include/symb_lib.h ../include/iritprsr.h
srf_mrch.o: ../include/trim_lib.h ../include/triv_lib.h ../include/trng_lib.h
srf_mrch.o: ../include/mdl_lib.h ../include/mvar_lib.h ../include/allocate.h
srf_mrch.o: ../include/obj_dpnd.h ../include/geom_lib.h ../include/attribut.h
srf_mrch.o: user_loc.h ../include/user_lib.h
srfray.o: ../include/irit_sm.h ../include/geom_lib.h ../include/iritprsr.h
srfray.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
srfray.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
srfray.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
srfray.o: ../include/attribut.h user_loc.h ../include/user_lib.h
user_err.o: ../include/irit_sm.h user_loc.h ../include/user_lib.h
user_err.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
user_err.o: ../include/symb_lib.h ../include/iritprsr.h ../include/trim_lib.h
user_err.o: ../include/triv_lib.h ../include/trng_lib.h ../include/mdl_lib.h
user_err.o: ../include/mvar_lib.h
user_ftl.o: ../include/irit_sm.h user_loc.h ../include/user_lib.h
user_ftl.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
user_ftl.o: ../include/symb_lib.h ../include/iritprsr.h ../include/trim_lib.h
user_ftl.o: ../include/triv_lib.h ../include/trng_lib.h ../include/mdl_lib.h
user_ftl.o: ../include/mvar_lib.h
userpick.o: ../include/irit_sm.h user_loc.h ../include/user_lib.h
userpick.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
userpick.o: ../include/symb_lib.h ../include/iritprsr.h ../include/trim_lib.h
userpick.o: ../include/triv_lib.h ../include/trng_lib.h ../include/mdl_lib.h
userpick.o: ../include/mvar_lib.h ../include/geom_lib.h ../include/attribut.h
usrcnvrt.o: ../include/irit_sm.h ../include/iritprsr.h ../include/cagd_lib.h
usrcnvrt.o: ../include/miscattr.h ../include/misc_lib.h ../include/symb_lib.h
usrcnvrt.o: ../include/trim_lib.h ../include/triv_lib.h ../include/trng_lib.h
usrcnvrt.o: ../include/mdl_lib.h ../include/mvar_lib.h ../include/allocate.h
usrcnvrt.o: ../include/obj_dpnd.h user_loc.h ../include/user_lib.h
visible.o: ../include/irit_sm.h ../include/allocate.h ../include/iritprsr.h
visible.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
visible.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
visible.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
visible.o: ../include/obj_dpnd.h ../include/attribut.h ../include/user_lib.h
visible.o: ../include/geom_lib.h
