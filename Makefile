CXX = $(shell if [ `uname -s` = 'Darwin' ]; then echo '/opt/local/bin/g++-mp-4.5'; else echo `which g++`; fi)
CXXFLAGS = -W -Wall -Werror -pedantic-errors -std=c++0x -g
LDFLAGS += $(shell Magick++-config --ldflags)
LDFLAGS += $(shell Magick++-config --libs)
DESTDIR ?= ../bin/$(shell uname -s)/$(shell uname -m)

MAGICK_INCLUDE = $(shell Magick++-config --cppflags)

#Compiler flags
#if mode variable is empty, setting debug build mode
ifeq ($(mode),release)
	CXXFLAGS += -O2 -DNDEBUG
else
	mode = debug
	CXXFLAGS += -ggdb -DDEBUG -pg -fno-omit-frame-pointer
endif

OBJDIR := $(DESTDIR)
OBJS := $(addprefix $(OBJDIR)/,img.$(mode).o naive.$(mode).o border.$(mode).o matrix.$(mode).o img_2d.$(mode).o test.$(mode).o)

.PHONY: all

all: information $(OBJDIR)/tester

information:
ifneq ($(mode),release)
ifneq ($(mode),debug)
	@echo "Invalid build mode."
	@echo "Please use 'make mode=release' or 'make mode=debug'"
	@exit 1
endif
endif
	@echo "Building on "$(mode)" mode"
	@echo ".........................."

$(OBJDIR)/tester: $(OBJDIR)/tester.$(mode)

$(OBJDIR)/tester.$(mode): $(OBJS)
	${CXX} ${CXXFLAGS} ${MAGICK_INCLUDE} $^ ${LDFLAGS} -o $@
	cp "$(OBJDIR)/tester.$(mode)" "$(OBJDIR)/tester"

$(OBJDIR)/img_2d.$(mode).o: img_2d.cc img_2d.h img.h Makefile
	${CXX} -c ${CXXFLAGS} ${MAGICK_INCLUDE} -o $@ $<

$(OBJDIR)/mm.$(mode).o: mm.cc img.cc
	${CXX} -c ${CXXFLAGS} -o $@ $^

$(OBJDIR)/%.$(mode).o: %.cc %.h img.h Makefile
	${CXX} -c ${CXXFLAGS} -o $@ $<

$(OBJS): | $(OBJDIR)

$(OBJDIR):
	test ! -d "${DESTDIR}" && mkdir -p "${DESTDIR}" || sleep 0

.PHONY: help

help:
	@echo "usage: make [(mode=debug | mode=release) | help | clean]"
	@echo "	mode=debug:   generates binary with debug parameters"
	@echo "	mode=release: generates binary with optimizations"
	@echo "	help:         see this usage message"
	@echo "	clean:        removes output binaries"

.PHONY: clean

clean:
	rm -rf $(OBJS) $(OBJDIR)/tester*


