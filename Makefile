#   Makefile C9-1
#   This should work as a reasonable template for your work
#   on the project too.
#   Of course you can remove most of the comments once you
#   don't need them anymore, the file will be much more readable.


#  Global confituration

#  Automatically finds files in 'src' when looking for depandancies:

VPATH += src


#  CPLEX configuration

CPLEX_BASE_DIR := /logiciels/cplex/linux/x86/cplex-127/
#CPLEX_BASE_DIR := /usr/ensta/pack/ro/cplex/x86_64-linux/CPLEX_Studio_126


include cplex.mk

#  A list of 'special' targets:
#
#    clean     Remove all the temporary object files (.o)
#    mrproper  Clean everything, including the executable
#
#  You can execute those targets by typing 'make <target>'
#  (instead of simply 'make') at the command line.
#
#  The '.PHONY' target is special and internal to 'make'.
#  The purpose is that the actions for 'clean', 'doc'
#  and 'mrproper' are executed even if it happens to
#  be a file called 'clean', 'doc' or 'mrproper' in
#  the directory.
#  See the 'make' documentation for more information
#  (or google 'makefile PHONY').

.PHONY: clean mrproper doc

#  The list of the object files of the program
#  This list need to be updated if you add more files
#  to your program.
#
#  You can generate this list automatically from the
#  source file in 'src' with the 'wildcard' and
#  'patsubst' functions we saw in class (or, again,
#  google 'makefile wildcard' for example), although
#  I would advise you against it. (In particular
#  because it can be usefull to compile only a part of your
#  project when you have build errors.)

OBJ_LIST := InstanceUCP.o \
	    Separation.o \
	    ModeleUCP.o \
	    main.o \

#  The list of the path to the object files: since
#  all the object files should be in the 'build'
#  directory, we add the prefix 'build/' to all
#  the elements in OBJ_LIST.

OBJ      := $(addprefix build/, $(OBJ_LIST))

#   Name of the program (BIN stands for 'binary')

BIN      := mf

#   Compilation and link flags

CFLAGS   := -Wall -g
LDFLAGS  := $(CPLEX_CXX_LDFLAGS)

#   Aliases for functions.
#   This makes the Makefile a little bit more readable,
#   and allow to change e.g. the compiler easily.

compile  := g++ $(CFLAGS) $(CPLEX_CXX_CFLAGS) -c
link     := g++ $(CPLEX_CXX_LDDIRS)
remove   := rm -f

#   Compilation targets. As in class, except that
#   we now have separate directory: 'build' for
#   the object files and 'src' for the source files.

$(BIN):  $(OBJ)
	@$(link) $^ -o $@  $(LDFLAGS)

build/%.o: src/%.cc
	@$(compile) $< -o $@

#   The rules for the 'special' targets discussed
#   above (clean, mrproper, doc)

clean:
	@$(remove) $(OBJ)

mrproper: clean
	@$(remove) $(BIN)

#   Finally, the list of the depandancies of the projects.
#
#   For small projects, it is not difficult to do this by
#   hand: for each source file (.cpp), list all the headers
#   that it includes, directly or indirectly (through another
#   header).
#
#   Of course, this can be done automatically. The simplest
#   way with modern versions of gcc is to add the '-MP' and
#   '-MD' options to the compiler: this will create a file
#   with extension '.d' for each of your source file, that
#   you can then include (with the 'include' command) in
#   your makefile.
#   A more 'traditional' way (which will work with older
#   compilers) is described here:
#      http://mad-scientist.net/make/autodep.html
#   
#   But for this project, it is certainly less work for
#   you to just carefully write all your depandancies
#   by hand.


build/InstanceUCP.o   : InstanceUCP.h
build/ModeleUCP.o : InstanceUCP.h ModeleUCP.h
build/Separation.o : InstanceUCP.h Separation.h
build/main.o      : InstanceUCP.h ModeleUCP.h Separation.h 
