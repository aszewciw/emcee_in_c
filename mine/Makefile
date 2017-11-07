# Project subdirectories
SRCDIR = ./
INCDIR = ./
OBJDIR = ./
BINDIR = ./

# Compiler rules
CC=mpicc
CFLAGS=-cc=icc -I$(INCDIR) -Wall -xHost -O3 -vec_report2
# CFLAGS=-Wall -I$(INCDIR)

# Linker Flags and Libraries
LDFLAGS=-limf -lgsl -lgslcblas
LDLIBS=

# Name of final Executable
TARGET=run_mcmc

# inc files and object files
SOURCES := $(wildcard $(SRCDIR)/*.c)
INCLUDES := $(wildcard $(INCDIR)/*.h)
OBJECTS := $(SOURCES:$(SRCDIR)/%.c=$(OBJDIR)/%.o)

# Notes: $@ refers to the file name of the target of the rule
#        $< refers to the name of the first prerequisite
$(BINDIR)/$(TARGET): $(OBJECTS)
	$(CC) $(LDFLAGS) $(LDLIBS) $(OBJECTS) -o $@

# This is a poor way to do this, actually, because an object
# file depends on all header files. I suppose we don't expect
# the headers to change very often. It would be cleaner to build
# a list of header files that each source file depends on, especially
# since this will be a relatively small project.
$(OBJDIR)/%.o: $(SRCDIR)/%.c $(INCLUDES)
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: clean cleanall

clean:
	rm -f $(BINDIR)/$(TARGET)

cleanall: clean
	rm -f $(OBJECTS)