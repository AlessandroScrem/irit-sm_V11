#
# This is the make file for the cagd_lib library.
#
#				Gershon Elber, Aug 1990
#

#
# This library has a large table of I Choose K in cbzrtbl2.c.  If your
# compiler fails to compile this file, try to replace it by cbzrtbl1.c
# and change CAGD_MAX_BEZIER_CACHE_ORDER accordingly in cagd_lib.h
#

include ../makeflag.sas

OBJS =  afd_cube.o bez_clip.o blossom.o bsp2poly.o bsp_gen.o bsp_knot.o \
	bspboehm.o bspcoxdb.o bzr2poly.o bzr_gen.o bzr_intr.o \
	bzr_pwr.o cagd2ply.o cagd2pl2.o cagd_arc.o \
	cagd_aux.o cagd_cci.o cagd_cnc.o \
	cagd_dbg.o cagd_err.o cagd_ftl.o \
	cagd1gen.o cagd2gen.o cagdbbox.o cagdbsum.o cagdcmpt.o \
	cagdcmrg.o cagdcoer.o cagdcsrf.o cagdedit.o cagdextr.o cagdmesh.o \
	cagdoslo.o cagdprim.o cagdruld.o cagdsmrg.o cagdsrev.o \
	cagdswep.o cbsp_aux.o cbsp_fit.o cbsp_int.o cbspeval.o cbzr_aux.o \
	cbzr_tbl.o cbzr2tbl.o cbzreval.o cpwr_aux.o \
	crvmatch.o hermite.o mshplanr.o nrmleval.o \
	poly_err.o sbsp_aux.o sbsp_int.o sbspeval.o sbzr_aux.o sbzreval.o

all:	cagd.lib

cagd.lib: $(OBJS)
	rm -f cagd.lib
	oml cagd.lib a $(OBJS)

install: cagd.lib
	mv -f cagd.lib $(IRIT_LIB_DIR)

# DO NOT DELETE THIS LINE -- make depend depends on it.

afd_cube.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
afd_cube.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
afd_cube.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
afd_cube.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
blossom.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
blossom.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
blossom.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
blossom.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
bsp2poly.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
bsp2poly.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
bsp2poly.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
bsp2poly.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
bsp_gen.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
bsp_gen.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
bsp_gen.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
bsp_gen.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
bsp_knot.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
bsp_knot.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
bsp_knot.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
bsp_knot.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
bspboehm.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
bspboehm.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
bspboehm.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
bspboehm.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
bspcoxdb.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
bspcoxdb.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
bspcoxdb.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
bspcoxdb.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
bzr2poly.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
bzr2poly.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
bzr2poly.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
bzr2poly.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
bzr_gen.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
bzr_gen.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
bzr_gen.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
bzr_gen.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
bzr_pwr.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
bzr_pwr.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
bzr_pwr.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
bzr_pwr.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
cagd1gen.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
cagd1gen.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
cagd1gen.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
cagd1gen.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
cagd1gen.o: ../include/geom_lib.h ../include/attribut.h
cagd2gen.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
cagd2gen.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
cagd2gen.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
cagd2gen.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
cagd2gen.o: ../include/geom_lib.h ../include/attribut.h
cagd_arc.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
cagd_arc.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
cagd_arc.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
cagd_arc.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
cagd_aux.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
cagd_aux.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
cagd_aux.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
cagd_aux.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
cagd_cci.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
cagd_cci.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
cagd_cci.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
cagd_cci.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
cagd_cci.o: ../include/geom_lib.h ../include/attribut.h
cagd_dbg.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
cagd_dbg.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
cagd_dbg.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
cagd_dbg.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
cagd_err.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
cagd_err.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
cagd_err.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
cagd_err.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
cagd_ftl.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
cagd_ftl.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
cagd_ftl.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
cagd_ftl.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
cagdbbox.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
cagdbbox.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
cagdbbox.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
cagdbbox.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
cagdbsum.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
cagdbsum.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
cagdbsum.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
cagdbsum.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
cagdcmpt.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
cagdcmpt.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
cagdcmpt.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
cagdcmpt.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
cagdcmrg.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
cagdcmrg.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
cagdcmrg.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
cagdcmrg.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
cagdcmrg.o: ../include/geom_lib.h ../include/attribut.h
cagdcoer.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
cagdcoer.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
cagdcoer.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
cagdcoer.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
cagdcsrf.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
cagdcsrf.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
cagdcsrf.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
cagdcsrf.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
cagdedit.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
cagdedit.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
cagdedit.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
cagdedit.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
cagdextr.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
cagdextr.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
cagdextr.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
cagdextr.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
cagdmesh.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
cagdmesh.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
cagdmesh.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
cagdmesh.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
cagdoslo.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
cagdoslo.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
cagdoslo.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
cagdoslo.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
cagdprim.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
cagdprim.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
cagdprim.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
cagdprim.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
cagdruld.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
cagdruld.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
cagdruld.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
cagdruld.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
cagdsmrg.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
cagdsmrg.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
cagdsmrg.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
cagdsmrg.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
cagdsrev.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
cagdsrev.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
cagdsrev.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
cagdsrev.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
cagdsrev.o: ../include/geom_lib.h ../include/attribut.h
cagdswep.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
cagdswep.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
cagdswep.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
cagdswep.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
cbsp_aux.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
cbsp_aux.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
cbsp_aux.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
cbsp_aux.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
cbsp_int.o: ../include/geom_lib.h ../include/iritprsr.h ../include/irit_sm.h
cbsp_int.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
cbsp_int.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
cbsp_int.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
cbsp_int.o: ../include/attribut.h cagd_loc.h ../include/extra_fn.h
cbspeval.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
cbspeval.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
cbspeval.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
cbspeval.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
cbzr_aux.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
cbzr_aux.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
cbzr_aux.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
cbzr_aux.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
cbzr_tbl.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
cbzr_tbl.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
cbzr_tbl.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
cbzr_tbl.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
cbzreval.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
cbzreval.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
cbzreval.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
cbzreval.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
crvmatch.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
crvmatch.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
crvmatch.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
crvmatch.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
crvmatch.o: ../include/geom_lib.h ../include/attribut.h
hermite.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
hermite.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
hermite.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
hermite.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
mshplanr.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
mshplanr.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
mshplanr.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
mshplanr.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
nrmleval.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
nrmleval.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
nrmleval.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
nrmleval.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
poly_err.o: ../include/irit_sm.h cagd_loc.h ../include/iritprsr.h
poly_err.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
poly_err.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
poly_err.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
poly_err.o: ../include/geom_lib.h ../include/attribut.h
sbsp_aux.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
sbsp_aux.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
sbsp_aux.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
sbsp_aux.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
sbsp_int.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
sbsp_int.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
sbsp_int.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
sbsp_int.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
sbsp_int.o: ../include/extra_fn.h
sbspeval.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
sbspeval.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
sbspeval.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
sbspeval.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
sbzr_aux.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
sbzr_aux.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
sbzr_aux.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
sbzr_aux.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
sbzreval.o: cagd_loc.h ../include/iritprsr.h ../include/irit_sm.h
sbzreval.o: ../include/cagd_lib.h ../include/miscattr.h ../include/misc_lib.h
sbzreval.o: ../include/symb_lib.h ../include/trim_lib.h ../include/triv_lib.h
sbzreval.o: ../include/trng_lib.h ../include/mdl_lib.h ../include/mvar_lib.h
