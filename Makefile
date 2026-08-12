# RapidMACS — standalone MACS3-compatible peak caller
# Produces:
#   lib/librapidmacs.a — static library
#   bin/rapidmacs      — standalone CLI
#
# htslib is referenced as a git submodule under third_party/htslib (pinned
# at tag 1.23). Set HTSLIB_DIR to point at an alternate build (e.g. an
# already-built copy in a sibling Chromap-suite checkout).

CXX ?= g++
AR ?= ar

# By default we link against system-installed htslib (libhts-dev provides
# kfunc.h and libhts.so). For a self-contained build, point HTSLIB_DIR at
# a built copy of third_party/htslib (or any other htslib install).
HTSLIB_DIR ?=

CXXFLAGS ?= -std=c++11 -Wall -Wextra -O3 -fopenmp -msse4.1
CPPFLAGS_LOCAL = -Iinclude
DEPFLAGS = -MMD -MP
LDFLAGS_LOCAL = -lhts -lm -lz -lpthread -fopenmp

ifneq ($(strip $(HTSLIB_DIR)),)
CPPFLAGS_LOCAL += -I$(HTSLIB_DIR)
LDFLAGS_LOCAL := -L$(HTSLIB_DIR) $(LDFLAGS_LOCAL)
endif

OBJDIR = build
LIBDIR = lib
BINDIR = bin

LIB_SOURCES = \
    src/bdgpeakcall.cc \
    src/binned_signal.cc \
    src/call_peaks.cc \
    src/frag_compact_store.cc \
    src/fragment_input.cc \
    src/fragments.cc \
    src/macs3_bed_callpeak.cc \
    src/macs3_frag_peak_pipeline.cc \
    src/macs3_frag_workspace.cc \
    src/macs3_frag_sweep.cc \
    src/macs3_pos_array.cc \
    src/peak_io.cc \
    src/signac_bed_projection.cc \
    src/stage_profile.cc \
    src/streaming.cc

LIB_OBJS = $(patsubst src/%.cc,$(OBJDIR)/%.o,$(LIB_SOURCES))

CLI_OBJ = $(OBJDIR)/cli_rapidmacs.o

LIB = $(LIBDIR)/librapidmacs.a
CLI = $(BINDIR)/rapidmacs

.PHONY: all clean rapidmacs librapidmacs \
        test test-unit \
        test-legacy-streaming-100k test-macs3-frag-100k \
        test-chipseq test-cutrun

all: rapidmacs

librapidmacs: $(LIB)

rapidmacs: $(CLI)

$(LIB): $(LIB_OBJS) | $(LIBDIR)
	$(AR) rcs $@ $^

$(CLI): $(CLI_OBJ) $(LIB) | $(BINDIR)
	$(CXX) $(CXXFLAGS) $(CLI_OBJ) $(LIB) -o $@ $(LDFLAGS_LOCAL)

$(CLI_OBJ): cli/rapidmacs.cc | $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS_LOCAL) $(DEPFLAGS) -c $< -o $@

$(OBJDIR)/%.o: src/%.cc | $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS_LOCAL) $(DEPFLAGS) -c $< -o $@

-include $(LIB_OBJS:.o=.d) $(CLI_OBJ:.o=.d)

$(OBJDIR) $(LIBDIR) $(BINDIR):
	mkdir -p $@

# 100K parity smoke for the LEGACY custom caller path (--out-prefix). Exercises
# the streaming legacy-caller path (uint16 bins, sliding-window lambda).
test-legacy-streaming-100k: $(CLI)
	tests/run_100k_legacy_caller_smoke.sh

# 100K parity smoke for the MACS3 FRAG path (--macs3-frag-narrowpeak +
# --macs3-frag-summits) — the path libchromap calls in production. Currently
# still runs in-memory.
test-macs3-frag-100k: $(CLI)
	tests/run_100k_macs3_frag_smoke.sh

test-chipseq: $(CLI)
	tests/run_chipseq_smoke.sh

test-cutrun: $(CLI)
	tests/run_cutrun_smoke.sh

test-keepdup: $(CLI)
	tests/run_keepdup_regression.sh

# Unit tests (compiled binaries under build/tests/, run sequentially).
test-unit: $(LIB)
	tests/run_unit_tests.sh

test: test-unit test-keepdup test-legacy-streaming-100k test-macs3-frag-100k \
      test-chipseq test-cutrun

clean:
	rm -rf $(OBJDIR) $(LIBDIR) $(BINDIR)
