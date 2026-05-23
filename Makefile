# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dan Alec Yamaguchi
#
# Licensed under GNU Affero General Public License v3.0.
# Commercial licenses available for enterprise use.
# Contact: danalec@gmail.com
# See LICENSE-COMMERCIAL for details.
#
# Makefile - RH Proof Build System (Paper Appendix A)
#
# Build all:     make all
# Build one:     make derive_k
# Clean:         make clean
# Run tests:     make test
#
# Requirements: GMP (optional, for derive_k_gmp)
#   Linux:   apt install libgmp-dev
#   macOS:   brew install gmp
#   Windows: Strawberry Perl provides GMP (C:/Strawberry/c/)

CC = gcc
CFLAGS = -O3 -Wall -Wextra -Wconversion -Wshadow -Isrc

SRCDIR  = src
OBJDIR  = obj
BINDIR  = bin

# OS detection via environment variable
ifdef OS
  DETECTED_OS := Windows
  GMP_INC ?= C:/Strawberry/c/include
  GMP_LIB ?= C:/Strawberry/c/lib
  CFLAGS += -I$(GMP_INC)
  LDFLAGS_GMP = -L$(GMP_LIB) -lgmp -lm
  LDFLAGS_STD = -lm
  EXE = .exe
else
  UNAME_S := $(shell uname -s 2>/dev/null || echo Unknown)
  ifeq ($(UNAME_S),Linux)
    DETECTED_OS := Linux
    LDFLAGS_GMP = $(shell pkg-config --libs gmp 2>/dev/null || echo -lgmp) -lm
    LDFLAGS_STD = -lm
    EXE =
  endif
  ifeq ($(UNAME_S),Darwin)
    DETECTED_OS := macOS
    GMP_INC ?= /opt/homebrew/include
    GMP_LIB ?= /opt/homebrew/lib
    CFLAGS += -I$(GMP_INC)
    LDFLAGS_GMP = -L$(GMP_LIB) -lgmp -lm
    LDFLAGS_STD = -lm
    EXE =
  endif
  ifneq ($(UNAME_S),Linux)
    ifneq ($(UNAME_S),Darwin)
      DETECTED_OS := Unix
      LDFLAGS_GMP = -lgmp -lm
      LDFLAGS_STD = -lm
      EXE =
    endif
  endif
endif

# Paper Appendix A: Core verification programs (12 standard + 1 GMP)
TARGETS_STD = \
	derive_k \
	derive_k2 \
	trace_verify \
	weyl_law_verify \
	heat_kernel_expansion \
	tauberian_argument \
	trace_error_bound \
	prove_epsilon_zero_closure \
	test_epsilon_paths \
	test_fejer_prime_sum \
	prove_path_a_determinant \
	prove_path_b_gaussian

TARGETS_GMP = derive_k_gmp

ALL_TARGETS = $(TARGETS_STD) $(TARGETS_GMP)

# Default: build all programs
all: dirs std gmp
	@echo ""
	@echo "Built $(words $(ALL_TARGETS)) programs in $(BINDIR)/"
	@echo "  Standard: $(words $(TARGETS_STD))"
	@echo "  GMP:      $(words $(TARGETS_GMP))"

dirs:
ifdef OS
	@-mkdir $(OBJDIR) 2>nul
	@-mkdir $(BINDIR) 2>nul
else
	@mkdir -p $(OBJDIR)
	@mkdir -p $(BINDIR)
endif

# Standard programs (no GMP)
std: dirs
	@echo "Building standard programs ($(words $(TARGETS_STD)))..."
	@for prog in $(TARGETS_STD); do \
		echo "  $$prog"; \
		$(CC) $(CFLAGS) $(SRCDIR)/$$prog.c -o $(BINDIR)/$$prog$(EXE) $(LDFLAGS_STD); \
	done

# GMP programs
gmp: dirs
	@echo "Building GMP-dependent programs ($(words $(TARGETS_GMP)))..."
	@for prog in $(TARGETS_GMP); do \
		echo "  $$prog"; \
		$(CC) $(CFLAGS) $(SRCDIR)/$$prog.c -o $(BINDIR)/$$prog$(EXE) $(LDFLAGS_GMP); \
	done

# Run verification suite
test: all
	@echo "Running verification suite..."
	@for prog in $(TARGETS_STD); do \
		echo ""; \
		echo "=== $$prog ==="; \
		$(BINDIR)/$$prog$(EXE) 2>&1 | head -20; \
	done

# Convenience targets
derive_k: dirs
	$(CC) $(CFLAGS) $(SRCDIR)/derive_k.c -o $(BINDIR)/derive_k$(EXE) $(LDFLAGS_STD)

derive_k2: dirs
	$(CC) $(CFLAGS) $(SRCDIR)/derive_k2.c -o $(BINDIR)/derive_k2$(EXE) $(LDFLAGS_STD)

trace: dirs
	$(CC) $(CFLAGS) $(SRCDIR)/trace_verify.c -o $(BINDIR)/trace_verify$(EXE) $(LDFLAGS_STD)

weyl: dirs
	$(CC) $(CFLAGS) $(SRCDIR)/weyl_law_verify.c -o $(BINDIR)/weyl_law_verify$(EXE) $(LDFLAGS_STD)

# Clean
clean:
ifdef OS
	@-del /q $(OBJDIR)\*.o 2>nul
	@-del /q $(BINDIR)\*$(EXE) 2>nul
else
	@rm -f $(OBJDIR)/*.o $(BINDIR)/*$(EXE) 2>/dev/null || true
endif
	@echo "Cleaned."

distclean:
ifdef OS
	@-rd /s /q $(OBJDIR) 2>nul
	@-rd /s /q $(BINDIR) 2>nul
else
	@rm -rf $(OBJDIR) $(BINDIR)
endif
	@echo "Removed $(OBJDIR)/ and $(BINDIR)/."

# Info
info:
	@echo "RH Proof Build System (Paper Appendix A)"
	@echo "========================================="
	@echo "OS:              $(DETECTED_OS)"
	@echo "Compiler:        $(CC)"
	@echo "CFLAGS:          $(CFLAGS)"
	@echo "GMP Link:        $(LDFLAGS_GMP)"
	@echo ""
	@echo "Programs: $(words $(ALL_TARGETS)) total"
	@echo "  Standard: $(TARGETS_STD)"
	@echo "  GMP:      $(TARGETS_GMP)"

.PHONY: all dirs std gmp test clean distclean info derive_k derive_k2 trace weyl