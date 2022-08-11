# This file is part of the Astrometry.net suite.
# Licensed under a 3-clause BSD style license - see LICENSE

BASEDIR := ..
COMMON := $(BASEDIR)/util

CATALOGS := $(BASEDIR)/catalogs
CATS_INC := -I$(CATALOGS)
CATS_SLIB := $(CATALOGS)/libcatalogs.a
CATS_LIB :=

.PHONY: all
all:

# Detect GSL -- this minimum version was chosen to match the version in gsl-an.
# Earlier versions would probably work fine.
SYSTEM_GSL ?= $(shell (pkg-config --atleast-version=1.14 gsl && echo "yes") || echo "no")
# Make this variable visible to recursive "make" calls
export SYSTEM_GSL

include $(COMMON)/makefile.common
include $(COMMON)/makefile.anfiles
include $(COMMON)/makefile.cfitsio
# only for tweak-main.c
include $(COMMON)/makefile.png
include $(COMMON)/makefile.cairo
include $(COMMON)/makefile.jpeg
include $(COMMON)/makefile.netpbm
include $(COMMON)/makefile.cairo

ifneq ($(MAKECMDGOALS),clean)
    include $(COMMON)/makefile.os-features
endif

$(COMMON)/makefile.os-features:
	$(MAKE) -C $(COMMON) makefile.os-features

SO=$(SHAREDLIB_SUFFIX)

ENGINE_LIB := libastrometry.a
ENGINE_SO := libastrometry.$(SO)

LDFLAGS := $(LDFLAGS_DEF)

LDLIBS := $(LDLIBS_DEF)
LDLIBS += $(ANFILES_LIB)

SLIB := $(ENGINE_LIB)
SLIB += $(CATS_SLIB)
SLIB += $(ANFILES_SLIB)

CFLAGS += $(CFLAGS_DEF)
CFLAGS += $(CATS_INC)
CFLAGS += $(ANFILES_CFLAGS)
CFLAGS += $(ANFILES_INC)
# boilerplate.h
CFLAGS += -I$(ANUTILS_DIR)
# wcs-resample.h
CFLAGS += -I.
CFLAGS += $(CFITS_INC)

SHAREDLIBFLAGS := $(SHAREDLIBFLAGS_DEF)

ALL_WCSLIB_TARGETS :=

# inc:
# 	echo CFLAGS_DEF: $(CFLAGS_DEG)
# 	echo CATS_INC: $(CATS_INC)
# 	echo ANFILES_CFLAGS: $(ANFILES_CFLAGS)
# 	echo ANFILES_INC: $(ANFILES_INC)
# 	echo CFITS_INC: $(CFITS_INC)
# 	echo ZLIB_INC: $(ZLIB_INC)

CFITS_UTILS := tablist modhead fitscopy tabmerge liststruc \
		listhead imcopy imarith imstat
# fitstomatlab 

FITS_UTILS := $(CFITS_UTILS)

MISC_EXECS := project-usnob diffractionFlag_check hpgrid resort-xylist

resort-xylist: resort-xylist-main.o resort-xylist.o $(SLIB)

UTIL_OBJS := 

OTHER_OBJS := catalog.o codefile.o verify.o \
	solver.o solvedfile.o pnpoly.o tweak.o \
	quadcenters.o startree2rdls.o \
	solverutils.o engine-main.o engine.o tweak2.o

NOT_INSTALLED_PIPELINE := agreeable certifiable \
	hpquads codetree unpermute-quads unpermute-stars \
	printsolved mergesolved subwcs \
	augment-xylist merge-index index-to-table setsolved \
	uniformize-catalog local-index index-info control-program

PIPELINE := wcs-grab solve-field

MAIN_PROGS := image2xy new-wcs fits-guess-scale startree
# hpquads

SIMPLE_PROGS := wcs-grab get-wcs query-starkd
# hpowned

PROGS := astrometry-engine build-astrometry-index \
	$(MAIN_PROGS) $(SIMPLE_PROGS)

PROSPECTUS := quadidx codeprojections quadscales quadsperstar \
	quadcenters startree2rdls

ANLIBS := $(ANFILES_LIB) $(LIBKD_LIB) $(ANUTILS_LIB) $(GSL_LIB) $(QFITS_LIB)

INSTALL_LIB := $(ENGINE_LIB) $(ENGINE_SO)

ENGINE_OBJS := \
		engine.o solverutils.o onefield.o solver.o quad-utils.o \
		solvedfile.o tweak2.o \
		verify.o tweak.o

# These are required by solve-field and friends
ENGINE_OBJS += new-wcs.o fits-guess-scale.o cut-table.o \
	resort-xylist.o

BUILD_INDEX_OBJS := build-index.o uniformize-catalog.o startree.o hpquads.o \
	quad-builder.o quad-utils.o codefile.o codetree.o unpermute-stars.o \
	unpermute-quads.o merge-index.o 
ENGINE_OBJS += $(BUILD_INDEX_OBJS)

#augment-xylist.o

$(ENGINE_LIB): $(ENGINE_OBJS)
	-rm -f $@
	$(AR) rc $@ $(ENGINE_OBJS)
	$(RANLIB) $@

$(ENGINE_SO): $(ENGINE_OBJS) $(SLIB)
	$(CC) $(LDFLAGS) $(SHAREDLIBFLAGS) -o $@ $^ $(LDLIBS)

# old and miscellaneous executables that aren't part of the pipeline.
OLDEXECS := checkquads
OLDEXECS_OBJS := catalog.o verify.o $(UTIL_OBJS)

PIPELINE_MAIN_OBJ := $(addsuffix .o,$(PIPELINE))
PROSPECTUS_MAIN_OBJ := $(addsuffix .o,$(PROSPECTUS))
FITS_UTILS_MAIN_OBJ := $(addsuffix .o,$(FITS_UTILS))

INSTALL_EXECS := $(FITS_UTILS) fitsverify $(PIPELINE) $(PROGS)

INSTALL_H := allquads.h augment-xylist.h axyfile.h \
	engine.h onefield.h solverutils.h build-index.h catalog.h \
	codefile.h codetree.h fits-guess-scale.h hpquads.h \
	image2xy-files.h merge-index.h \
	new-wcs.h quad-builder.h quad-utils.h resort-xylist.h \
	solvedfile.h solver.h tweak.h uniformize-catalog.h \
	unpermute-quads.h unpermute-stars.h verify.h \
	tweak2.h

ALL_OBJ := $(UTIL_OBJS) $(KDTREE_OBJS) $(QFITS_OBJ) \
	$(PIPELINE_MAIN_OBJ) $(PROSPECTUS_MAIN_OBJ) $(FITS_UTILS_MAIN_OBJ) \
	$(OTHER_OBJS)
ALL_EXECS :=

all: $(QFITS_SLIB) $(LIBKD_LIB_FILE) \
	$(ANUTILS_LIB_FILE) $(ANFILES_LIB_FILE) \
	$(PROGS) $(PIPELINE) $(FITS_UTILS) fitsverify \
	$(ENGINE_LIB) $(ENGINE_SO) augment-xylist

$(MAIN_PROGS): %: %-main.o $(SLIB)
ALL_OBJ += $(addsuffix -main.o,$(MAIN_PROGS))

$(SIMPLE_PROGS): %: %.o $(SLIB)
ALL_OBJ += $(addsuffix .o,$(SIMPLE_PROGS))

.PHONY: pipeline
pipeline: $(PIPELINE)

.PHONY: prospectus_progs
prospectus_progs: $(PROSPECTUS)

.PHONY: fitsutils
fitsutils: $(FITS_UTILS)

simplexy: image2xy
.PHONY: simplexy

PY_INSTALL_DIR := $(PY_BASE_INSTALL_DIR)/solver
LINK_DIR := $(PY_BASE_LINK_DIR)/solver

PYTHON_EXECS := 
PYTHON_INSTALL := $(PYTHON_EXECS) __init__.py

install: $(INSTALL_EXECS) $(INSTALL_LIB)
	@echo Installing in directory '$(INSTALL_DIR)'
	$(MKDIR) '$(INSTALL_DIR)/bin'
	@for x in $(INSTALL_EXECS); do \
		echo $(CP) $$x '$(INSTALL_DIR)/bin'; \
		$(CP) $$x '$(INSTALL_DIR)/bin'; \
	done
	$(MKDIR) '$(ETC_INSTALL_DIR)'
	$(PYTHON) -c "import os; print(open('../etc/astrometry.cfg-dist').read().replace('DATA_INSTALL_DIR', '$(DATA_FINAL_DIR)'))" > '$(ETC_INSTALL_DIR)/astrometry.cfg'
	$(MKDIR) '$(INCLUDE_INSTALL_DIR)'
	@for x in $(INSTALL_H); do \
		echo $(CP) '$(INCLUDE_DIR)/'$$x '$(INCLUDE_INSTALL_DIR)'; \
		$(CP) '$(INCLUDE_DIR)/'$$x '$(INCLUDE_INSTALL_DIR)'; \
	done
	$(MKDIR) '$(LIB_INSTALL_DIR)'
	@for x in $(INSTALL_LIB); do \
		echo $(CP) $$x '$(LIB_INSTALL_DIR)'; \
		$(CP) $$x '$(LIB_INSTALL_DIR)'; \
	done
	$(MKDIR) '$(PY_INSTALL_DIR)'
	@for x in $(PYTHON_INSTALL); do \
		echo $(CP) $$x '$(PY_INSTALL_DIR)/'$$x; \
		$(CP) $$x '$(PY_INSTALL_DIR)/'$$x; \
	done
	@echo Making symlinks in directory '$(BIN_INSTALL_DIR)'
	$(MKDIR) '$(BIN_INSTALL_DIR)'
	@for x in $(PYTHON_EXECS); do \
		echo ln -f -s '$(LINK_DIR)/'$$x '$(BIN_INSTALL_DIR)/'$$x; \
		ln -f -s '$(LINK_DIR)/'$$x '$(BIN_INSTALL_DIR)/'$$x; \
	done

test-solver: test-solver.o solver_test.o $(SLIB)
test-solver-2: test-solver-2.o solver_test_2.o $(SLIB)

NODEP_OBJS += solver_test.o solver_test_2.o
ALL_OBJ += test-solver.o test-solver-2.o

CFLAGS_DEBUG = $(subst -DNDEBUG,,$(CFLAGS))

test-solver.o: test-solver.c
	$(CC) $(CPPFLAGS) $(CFLAGS_DEBUG) -o $@ -c $<
test-solver-2.o: test-solver-2.c
	$(CC) $(CPPFLAGS) $(CFLAGS_DEBUG) -o $@ -c $<

solver_test.o: solver.c
	$(CC) $(CPPFLAGS) $(CFLAGS_DEBUG) -DTESTING=1 -DTESTING_TRYALLCODES=1 -o $@ -c $<
solver_test_2.o: solver.c
	$(CC) $(CPPFLAGS) $(CFLAGS_DEBUG) -DTESTING=1 -DTESTING_TRYPERMUTATIONS=1 -o $@ -c $<


#######################################

# Add the basename of your test sources here...
ALL_TEST_FILES = test_solverutils \
	test_resort-xylist test_tweak test_multiindex2 test_predistort

#test_xscale -- requires a large index file...

#test_codefile -- takes a long time

$(ALL_TEST_FILES): $(SLIB)

ALL_TEST_EXTRA_OBJS :=
ALL_TEST_EXTRA_LDFLAGS := -lm
ALL_TEST_LIBS := $(SLIB)

# Add the dependencies here...
#test_multiindex2: test_multiindex2.o $(SLIB)

tests: $(ALL_TEST_FILES)
.PHONY: tests

include $(COMMON)/makefile.tests

ALL_OBJ += $(ALL_TEST_FILES_O) $(ALL_TEST_FILES_MAIN_O)
ALL_OBJ += $(ALL_TEST_EXTRA_OBJS) test.o

GENERATED_FILES += $(ALL_TEST_FILES_MAIN_C)

ALL_EXECS += $(ALL_TEST_FILES) test

# END OF KEIR'S SUPER SIMPLE TESTING FRAMEWORK
######################################################

demo_dsmooth: demo_dsmooth.o $(CAIRO_SLIB)
	$(CC) -o $@ $(LDFLAGS) $^ $(CAIRO_LIBS)

control-program: control-program.o $(SLIB)
ALL_OBJ += control-program.o

whynot: whynot.o $(PLOTSTUFF) $(CAIRO_SLIB) $(SLIB)
	$(CC) -o $@ $(LDFLAGS) $^ $(CAIRO_LIBS) $(LDLIBS)

quadidx: quadidx.o $(SLIB)

image2xy: image2xy-main.o image2xy-files.o $(CFITS_SLIB) $(SLIB)
	$(CC) -o $@ $(LDFLAGS) $^ $(CFITS_LIB) $(LDLIBS)

fit-wcs: fit-wcs-main.o $(SLIB)
	$(CC) -o $@ $(LDFLAGS) $^ $(LDLIBS)

ALL_OBJ += image2xy-files.o

hpgrid: hpgrid.o $(SLIB)

$(CFITS_UTILS) :: %: %.o $(CFITS_SLIB)
	$(CC) -o $@ $(CFLAGS) $(LDFLAGS) $^ $(CFITS_LIB) $(LDLIBS)

fitsverify: ftverify.c fvrf_data.c fvrf_file.c fvrf_head.c fvrf_key.c fvrf_misc.c $(CFITS_SLIB)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) -DSTANDALONE -trigraphs $(CFITS_INC) -o $@ $^ $(CFITS_LIB) -lm

$(OLDEXECS) :: %: %.o $(OLDEXECS_OBJS) $(SLIB)

build-astrometry-index: build-index-main.o $(SLIB) #$(BUILD_INDEX_OBJS) $(SLIB)
	$(CC) -o $@ $(LDFLAGS) $^ $(LDLIBS)
ALL_OBJ += build-index-main.o

allquads: allquads-main.o allquads.o $(SLIB)
	$(CC) -o $@ $(LDFLAGS) $^ $(LDLIBS)
ALL_OBJ += allquads-main.o allquads.o

unpermute-stars: unpermute-stars-main.o $(SLIB)
	$(CC) -o $@ $(LDFLAGS) $^ $(LDLIBS)
ALL_OBJ += unpermute-stars-main.o

astrometry-engine: engine-main.o $(SLIB)
	$(CC) -o $@ $(CFLAGS) $(LDFLAGS) $^ $(LDLIBS)

solve-field: solve-field.o augment-xylist.o image2xy-files.o $(SLIB) $(CFITS_SLIB)
	$(CC) -o $@ $(CFLAGS) $(LDFLAGS) $^ $(CFITS_LIB) $(LDLIBS)
ALL_OBJ += solve-field.o image2xy-files.o

augment-xylist: augment-xylist-main.o augment-xylist.o image2xy-files.o \
		$(SLIB) $(CFITS_SLIB)
	$(CC) -o $@ $(CFLAGS) $(LDFLAGS) $^ $(CFITS_LIB) $(LDLIBS)
ALL_OBJ += augment-xylist-main.o augment-xylist.o

PLOTDIR=$(BASEDIR)/plot/
PLOTSTUFF=$(addprefix $(PLOTDIR),plotstuff.o plotannotations.o plotfill.o plotgrid.o plothealpix.o plotimage.o plotindex.o plotmatch.o plotoutline.o plotradec.o plotxy.o)
tweak: tweak-main.o $(SLIB) $(CFITS_SLIB)
	$(CC) -o $@ $(CFLAGS) $(LDFLAGS) $^ $(CFITS_LIB) $(PLOTSTUFF) $(COMMON)/cairoutils.o $(CAIRO_LIB) $(PNG_LIB) $(JPEG_LIB) $(NETPBM_LIB) $(LDLIBS)
tweak-main.o: tweak-main.c
	$(CC) -c -o $@ $(CFLAGS) $(CAIRO_INC) tweak-main.c
ALL_OBJ += tweak-main.o

py: pysolver
.PHONY: py

pysolver: _solver$(PYTHON_SO_EXT)
.PHONY: pysolver

_solver$(PYTHON_SO_EXT): solver.i $(ENGINE_OBJS)
	LDFLAGS="$(LDFLAGS)" LDLIBS="$(LDLIBS)" \
	SLIB="$(SLIB)" \
	INC="$(ANFILES_INC)" \
	CFLAGS="$(ANFILES_CFLAGS)" \
	$(PYTHON) setup-solver.py build_ext -v --inplace --build-temp .

PYTHON_EXTRA_INSTALL := _solver$(PYTHON_SO_EXT)

install-extra: $(INSTALL_CAIRO_EXECS)
	@echo Installing in directory '$(INSTALL_DIR)'
	$(MAKE) $(PYTHON_EXTRA_INSTALL)
	$(MKDIR) '$(PY_INSTALL_DIR)'
	@for x in $(PYTHON_EXTRA_INSTALL); do \
		echo $(CP) $$x '$(PY_INSTALL_DIR)/'$$x; \
		$(CP) $$x '$(PY_INSTALL_DIR)/'$$x; \
	done

DEP_OBJ := $(ALL_OBJ)
DEP_PREREQS := $(QFITS_LIB)

include $(COMMON)/makefile.deps

.PHONY: clean

clean:
	rm -f $(EXECS) $(EXTRA_EXECS) $(SOLVER_EXECS) $(MISC_EXECS) $(PROGS) \
		$(PIPELINE) $(PROSPECTUS) $(DEPS) $(FITS_UTILS) $(ALL_OBJ) \
		$(NODEP_OBJS) fitsverify \
		$(ALL_EXECS) $(GENERATED_FILES) $(ALL_TESTS_CLEAN) \
		$(ENGINE_LIB) $(ENGINE_SO) *.o *~ *.dep deps
