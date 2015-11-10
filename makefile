#	Makefile for c++ research programs
#
#   Default is to put all source files in Directory named 'source'
#    then compile all using "make all SDIR=source/"
#
#  on windows machine "make all SDIR=source/ APPEXT=.exe"
#
# the compiler to use.
CC=g++
# options passed to the compiler.
CFLAGS=-c -Wall
LDFLAGS=-Wall
# the file extension for executable .exe for windows
APPEXT=
# the source directory, blank if same as bin dir
SDIR=source/
# functions used by many programs
SOURCES=$(SDIR)xyz_analysis.cpp $(SDIR)xyz_manip.cpp $(SDIR)load_crdf.cpp
SOURCEF=$(SDIR)xyz_analysis     $(SDIR)xyz_manip     $(SDIR)load_crdf
SRC_OBJ=$(SOURCES:.cpp=.o)
SRC_HDR=$(SOURCES:.cpp=.h) $(SDIR)class_log.h
# programs that need the load_crdf.cpp
MLT_1=

# programs that need the load_crdf.cpp xyz_analysis.cpp xyz_manip.cpp
MLT_2=exciton_from_QM

# programs that stand alone
SNG_EXE=

.PHONY : set1 set2 set3 set4 all_exe VAR
all: $(SNG_EXE) $(SRC_OBJ) $(MLT_1) $(MLT_2)
	echo "`date`" > BUILD_DATE
set1: $(SNG_EXE)
set2: $(SDIR)$(SOURCES) $(MLT_1)
set3: $(SDIR)$(SOURCES) $(MLT_1) $(MLT_2)

VAR:
	@echo All stand alone programs $(SNG_EXE)
	@echo All funtion code $(SDIR)$(SOURCES)
	@echo All function code objects $(SDIR)$(SRC_OBJ)
	@echo All dependent code $(MLT_1) $(MLT_2)

.SECONDEXPANSION:
$(SNG_EXE): $(SDIR)$$@.cpp $(SDIR)class_log.h
	$(CC) $(LDFLAGS) $(SDIR)$@.cpp -o $@$(APPEXT)

$(SDIR)load_crd.o: $(SDIR)load_crd.cpp $(SDIR)load_crd.h $(SDIR)class_log.h
	$(CC) $(CFLAGS) $(@:.o=.cpp) -o $@

$(SDIR)xyz_analysis.o: $(SDIR)xyz_analysis.cpp $(SDIR)xyz_analysis.h $(SDIR)class_log.h
	$(CC) $(CFLAGS) $(@:.o=.cpp) -o $@

$(SDIR)xyz_manip.o: $(SDIR)xyz_manip.cpp $(SDIR)xyz_manip.h $(SDIR)class_log.h
	$(CC) $(CFLAGS) $(@:.o=.cpp) -o $@	

$(MLT_1): $(SDIR)$$@.cpp $(SDIR)load_crdf.o $(SDIR)load_crdf.h $(SDIR)class_log.h
	$(CC) $(CFLAGS) $(SDIR)$@.cpp -o $(SDIR)$@.o
	$(CC) $(LDFLAGS) $(SDIR)load_crdf.o $(SDIR)$@.o -o $@$(APPEXT)

$(MLT_2): $(SDIR)$$@.cpp $(SRC_OBJ) $(SDIR)class_log.h
	$(CC) $(CFLAGS) $(SDIR)$@.cpp -o $(SDIR)$@.o
	$(CC) $(LDFLAGS) $(SRC_OBJ) $(SDIR)$@.o -o $@$(APPEXT)

