#
# This is the make file for the mvar_lib library.
#
#                               Gershon Elber, May 1997
#

include ../makeflag.sas

MAT_OBJS = \
	mvbiscon.o \
	mvlowenv.o \
	mvsplmon.o \
	mvtrmbis.o \
	mvtrmpcr.o \
	mvvorcel.o \
	mvvorcrv.o

C2C_OBJS = \
	mv2cntct.o \
	mv2ctaux.o \
	mv2ctbvh.o

ZR_OBJS = \
	zret0d.o \
	zret1d.o \
	zret2d.o \
	zrmatlab.o \
	zrmv0d.o \
	zrmv1d.o \
	zrmv2dTJ.o \
	zrmv2dTp.o \
	zrmv2dTs.o \
	zrmvaux0.o \
	zrmvaux1.o \
	zrmvkant.o \
	zrsolver.o

OBJS =  contacts.o control.o crv_krnl.o flankmil.o hasdrf2d.o \
	hasdrf3d.o lndstsqr.o ms_circ.o ms_sphr.o mvaccess.o \
	mvar_aux.o mvar_dbg.o mvar_der.o mvar_det.o mvar_err.o \
	mvar_ftl.o mvar_gen.o mvar_int.o mvar_pll.o \
	mvar_ref.o mvar_rev.o mvar_sub.o \
	mvar_sym.o mvar_vec.o mvarcmpt.o mvarcoer.o mvardist.o mvaredit.o \
	mvareval.o mvarextr.o mvarintr.o mvarjimp.o \
	mvarmrph.o mvarpack.o mvarpck2.o mvarpole.o mvarproj.o \
	mvarprim.o mvarrais.o mvarsils.o mvarstpl.o mvartopo.o mvarzral.o \
	mv_crvtr.o mv_mat2d.o mvbisect.o mvbspsym.o mvbzrpwr.o \
	mvbzrsym.o mvcones.o mvtangnt.o mvtrivar.o offset2.o \
	offst2ni.o ray-trap.o raytrp3d.o selfintr.o skel2d.o

mvlowenv.o mvtrmbis.o mvvorcel.o mvvorcrv.o

all:	mvar.lib

mvar.lib: $(OBJS) $(MAT_OBJS) $(ZR_OBJS) $(C2C_OBJS)
	rm -f mvar.lib
	oml mvar.lib a $(OBJS) $(MAT_OBJS) $(ZR_OBJS) $(C2C_OBJS)

install: mvar.lib
	mv -f mvar.lib $(IRIT_LIB_DIR)


# DO NOT DELETE THIS LINE -- make depend depends on it.

