##=======================================================================
###############################################
# COPYRIGHT: See COPYRIGHT.txt                #
# 2015 by MCS, Argonne National Laboratory.   #
###############################################
##=======================================================================

##=======================================================================
##   PLEASE SET THESE VARIABLES BEFORE COMPILING
##=======================================================================

SZPATH		= /home/arham23/Software/SZ3/include/SZ3
HDF5PATH	= /home/arham23/Software/hdf5-1.12.1/hdf5

##=======================================================================
##   DIRECTORY TREE
##=======================================================================

LIB		= lib
OBJ		= obj
SRC		= src
INC		= include

##=======================================================================
##   COMPILERS
##=======================================================================

CC 		= gcc
MPICC 		= mpicc

##=======================================================================
##   FLAGS
##=======================================================================

SZFLAGS         = -I$(SZPATH)/api -L$(SZPATH) -Wl,-rpath,$(SZPATH)

HDF5FLAGS	= -I$(HDF5PATH)/include #$(HDF5PATH)/lib/libhdf5.a

##=======================================================================
##   TARGETS
##=======================================================================

OBJS		= $(OBJ)/H5Z_SZ3.o

SHARED		= libhdf5sz3.so
STATIC		= libhdf5sz3.a

all: 		$(LIB)/$(SHARED) $(LIB)/$(STATIC)

obj/H5Z_SZ3.o: src/H5Z_SZ3.cpp
	mkdir -p $(OBJ)
	$(CC) -c -fPIC -g -O3 -I./include $(HDF5FLAGS) $(SZFLAGS) $< -o $@

$(OBJ)/%.o:	$(SRC)/%.c
		@mkdir -p $(OBJ)
		$(CC) -c -fPIC -g -O3 -I./include $(HDF5FLAGS) $(SZFLAGS) $< -o $@
		
$(LIB)/$(SHARED):	$(OBJS)
		@mkdir -p $(LIB)
		$(CC) -O3 -shared -o $(LIB)/$(SHARED) $(OBJS) $(SZFLAGS) -L$(HDF5PATH)/lib -lc -lhdf5 -lzstd #-lzlib -lzstd

$(LIB)/$(STATIC):	$(OBJS)
		@mkdir -p $(LIB)
		$(RM) $@
		$(AR) -cvq $@ $(OBJS)
		
install: $(LIB)/$(SHARED) $(LIB)/$(STATIC)
		install -d $(SZPATH)/lib  $(SZPATH)/include
		install $(INC)/* $(SZPATH)/include/
		install $(LIB)/* $(SZPATH)/lib/

uninstall:
		$(RM) $(SZPATH)/$(LIB)/libhdf5sz3.* $(SZPATH)/$(INC)/H5Z_SZ3.h

clean:
		rm -f $(OBJ)/*.o $(LIB)/*

.PHONY:		$(SHARED) $(STATIC) install uninstall clean


