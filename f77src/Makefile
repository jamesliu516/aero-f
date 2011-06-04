LIB = libf77src.a

OBJS = $(LIB)(boundflux5.o) \
       $(LIB)(g_boundflux5.o) \
       $(LIB)(boundjac2.o) \
       $(LIB)(dsvdc.o) \
       $(LIB)(printmask.o) \
       $(LIB)(roeflux5.o) \
       $(LIB)(g_roeflux5.o) \
       $(LIB)(roeflux5prim.o) \
       $(LIB)(roeflux6.o) \
       $(LIB)(g_roeflux6.o) \
       $(LIB)(roejac2.o) \
       $(LIB)(roejac5.o) \
       $(LIB)(roejac6.o) \
       $(LIB)(roeflux5waterdissprec.o) \
       $(LIB)(roejac5waterdissprec.o) \
       $(LIB)(roejacwaterdissprim.o) \
       $(LIB)(roejac2water.o) \
       $(LIB)(eriemann_gw.o) \
       $(LIB)(eriemann_gg.o) \
       $(LIB)(eriemann_ww.o) \
       $(LIB)(torsionspring.o) \
       $(LIB)(ballvertex.o) \
       $(LIB)(genbcfluxgas.o) \
       $(LIB)(genbcfluxtait.o) \
       $(LIB)(genbcfluxjwl.o) \
       $(LIB)(roeflux5jwl.o) \
       $(LIB)(roejac6jwl.o) \
       $(LIB)(hlleflux.o) \
       $(LIB)(hllejac.o) \
       $(LIB)(hlleflux1.o) \
       $(LIB)(hllcflux.o)

ifeq ($(usescalapack), true)

  SOBJS = $(LIB)(thinsvd.o) \
          $(LIB)(lworksize.o)
else
  SBOJS =
endif

$(LIB): $(OBJS) $(SOBJS)
	echo $(LIB)

clean:
	rm -f $(LIB)
