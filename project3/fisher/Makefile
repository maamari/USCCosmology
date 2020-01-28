CFLAGS = -Wall -Wno-deprecated -I${GSLINCROOT}
OBJ = dnumrecipes.o dcomplex.o dcosmology2.o spline.o fisherCMB.o fisher.o fisherGAL.o neutrinos.o inflation.o fisherBAO.o fisherSN.o fisherVPAIR.o fisherGAL_CLUSTER.o #fnl_massf.o
#fisherGAL_CLUSTER.o fisherGAL_CLUSTER_fnl.o fnl_massf.o #fisherGAL_CLUSTER.o fisherGAL_CLUSTER_fnl.o #fnl_massf.o
CC = g++

GSLINCROOT = /usr/local/gsl/gcc4/64/1.16/include/gsl
GLIB = -L/usr/local/gsl/gcc4/64/1.16/lib/ -lgsl -lgslcblas

.SUFFIXES : .o .cc
.cc.o :
		$(CC) $(CFLAGS) -c $<
		# $(CC) -c $<
# all: driver_vpair.x
#
# driver_vpair.x : $(OBJ) driver_vpair.o
# 	$(CC) $(CFLAGS) -o $@ driver_vpair.o $(OBJ) $(GLIB) $(MISC)
# 	# $(CC) -o $@ driver_vpair.o $(OBJ) $(MISC)

all: driver.x

driver.x : $(OBJ) driver.o
	$(CC) $(CFLAGS) -o $@ driver.o $(OBJ) $(GLIB) $(MISC)
	# $(CC) -o $@ driver_vpair.o $(OBJ) $(MISC)

clean :
	rm -f *.o
