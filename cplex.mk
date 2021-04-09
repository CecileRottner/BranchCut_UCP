##   @file   cplex.mk
##   @brief  Template makefile to compile a cplex project.
##   @author Nicolas Grebille
##
##   Adapted from the Makefile of cplex examples.
##   This makefile is written for (and tested with)
##   GNU Make, but I don't think it contains GNU
##   specific features.
##
##   You need to adapt the "local configuration" section
##   to your system. It is not "really" clean, but it
##   seems better than to use a separate file or to ask
##   the user to define the proper variables (this file
##   is very small after all).
##
##   The best way to use this template for simple projects
##   is to include this file at the top of your makefile
##   and then just use the main exported variables:
##
##    CPLEX_CXX_CFLAGS    C++ compiler flags
##    CPLEX_CXX_LDFLAGS   C++ linker flags        (like `-llib`)
##    CPLEX_CXX_LDDIRS    C++ linker directories  (like `-L$(PATH_TO_LIB)`)
##    CPLEX_C_CFLAGS      C compiler flags
##    CPLEX_C_LDFLAGS     C linker flags
##    CPLEX_C_LDDIRS      C linker directories
##
##   To ensure there will not be name clashed,
##   all cplex-related variables are prefixed with either
##   the CPLEX_ or CONCERT_ prefix.
##
##   For example, if you want to compile a project in C++ for
##   a release build, you would do:
##   @code 
##         # Change those default if necessary:
##         CPLEX_SYSTEM     := x86-64_sles10_4.1
##         CPLEX_LIBFORMAT  := static_pic
##         CPLEX_BASE_DIR   := /usr/local/cplex
##
##         include "cplex.mk"
##
##         CXXFLAGS   += $(CPLEX_CXX_CFLAGS) -DNDEBUG
##         LDFLAGS    += $(CPLEX_CXX_LDFLAGS)
##         LDDIRS     += $(CPLEX_CXX_LDDIRS)
##         CXX        := g++
##
##         prog: $(OBJ)
##               $(CXX) $(LDDIRS) -o $@ $^ $(LDFLAGS)
##         %.o: %.c
##         	$(CXX) $(CXXFLAGS) $<
##   @endcode
##
##   and then just write your makefile as usual.
##
##   There is a distinction between LDFLAGS (-lcplex, etc)
##   and LDDIRS (-L$(CPLEX_LIB_DIR)) because I consider it is
##   "cleaner" to put the LDDIRS flag at the beginning of the
##   linking command and the LDFLAGS at the end.
##   This is not strictly necessary, however the LDFLAGS must
##   be at the end of the line and the order must be kept if
##   you wish to change this.
##
##   @warning
##     The flag `-DNDEBUG` is _removed_ from the $(CPLEX_CXX_FLAGS)
##     in this file, because we think you may want to enable or
##     disable it yourself as needed: this is a generic template,
##     and the aim is not to prevent the user to use different flags
##     for debug and release builds.
##     Thus, you need to include it in your CFLAGS if you want it.
##
##   @remarks
##     The 'compiler selection' section of the original makefile
##     (for cplex version 12.5) recommends the compiler:
##     @code
##      CXX = g++ -O0
##      CC  = gcc -O0
##     @endcode
##     And the C++ flags have only -O0 for the optimization level.
##     I left the '-O' option in the C++ flags for now.
##   @todo
##     Check in the docs:
##     *  if there is a reason not to use -O2 (or even -O3) with cplex
##     *  if it is possible to use another compiler than gcc.
##

# ---------------------------------------------------------------------
#      Local configuration
# ---------------------------------------------------------------------

CPLEX_SYSTEM     ?= x86-64_linux
CPLEX_LIBFORMAT  ?= static_pic
CPLEX_BASE_DIR   ?= /logiciels/cplex/linux/x86/cplex-126/

#CPLEX_BASE_DIR   = /logiciels/ilog/cplex_studio12.4

# ---------------------------------------------------------------------

CPLEX_DIR            = $(CPLEX_BASE_DIR)/cplex
CONCERT_DIR          = $(CPLEX_BASE_DIR)/concert
CPLEX_LIB_DIR        = $(CPLEX_DIR)/lib/$(CPLEX_SYSTEM)/$(CPLEX_LIBFORMAT)
CONCERT_LIB_DIR      = $(CONCERT_DIR)/lib/$(CPLEX_SYSTEM)/$(CPLEX_LIBFORMAT)
CPLEX_INCLUDE_DIR    = $(CPLEX_DIR)/include
CONCERT_INCLUDE_DIR  = $(CONCERT_DIR)/include


#  Add -DNDEBUG to those if you need it.

CPLEX_CXX_CFLAGS     = -I$(CPLEX_INCLUDE_DIR) -I$(CONCERT_INCLUDE_DIR) \
                       -m64 -fPIC -fno-strict-aliasing \
                       -fexceptions -O -DIL_STD

CPLEX_C_CFLAGS       = -I$(CPLEX_INCLUDE_DIR) \
                       -m64 -fPIC -fno-strict-aliasing 

CPLEX_CXX_LDDIRS     = -L$(CPLEX_LIB_DIR) -L$(CONCERT_LIB_DIR)
CPLEX_CXX_LDFLAGS    = -lconcert -lilocplex -lcplex -lm -pthread
CPLEX_C_LDDIRS       = -L$(CPLEX_LIB_DIR)
CPLEX_C_LDFLAGS      = -lcplex -lm -pthread
