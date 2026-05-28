# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dan Alec Yamaguchi
#
# Licensed under GNU Affero General Public License v3.0.
# Commercial licenses available for enterprise use.
# Contact: danalec@gmail.com
# See LICENSE-COMMERCIAL for details.
#
# ===========================================================================
# Cross-Platform Build Support
# ===========================================================================
# Verified compilers: gcc (all platforms), clang (Linux/macOS)
#
# Platform        | Compiler     | OpenMP flags              | Notes
# ----------------|--------------|---------------------------|------------------
# Windows (MSYS2) | gcc          | -fopenmp                  | Requires MinGW-w64
# Windows (MSYS2) | clang        | -fopenmp                  | Requires libomp
# Linux           | gcc          | -fopenmp                  |
# Linux           | clang        | -fopenmp                  | Requires libomp
# macOS           | clang        | -Xpreprocessor -fopenmp   | Requires libomp (brew)
# macOS           | gcc (Homebrew)| -fopenmp                 |
#
# Portability:
#   - All math constants (M_PI, M_E, M_SQRT2) have portable fallbacks
#   - No POSIX-only headers in src/ (no unistd.h, sys/time.h, etc.)
#   - Windows UTF-8 console (SetConsoleOutputCP) guarded by #ifdef _WIN32
#   - OpenMP pragmas gracefully ignored when compiler lacks support
#   - MPFR/GMP required (arb.h) — install via package manager
# ===========================================================================

CC      ?= gcc
LDFLAGS ?= -lm

SRCDIR  = src
BINDIR  = bin

# Detect compiler type dynamically
IS_CLANG := $(shell $(CC) --version 2>/dev/null | grep -i -q "clang" && echo "yes" || echo "no")

ifeq ($(IS_CLANG),yes)
  LOOP_FLAGS :=
else
  LOOP_FLAGS := -fno-peel-loops -fno-unswitch-loops
endif

CFLAGS  ?= -O3 -Wall -Wextra -Wconversion -Wshadow -Werror \
	  -fno-strict-aliasing $(LOOP_FLAGS) -Isrc

UNAME_S := $(shell uname -s 2>/dev/null || echo "Windows")

# Set executable extension
ifneq (,$(filter %MINGW %MSYS %Windows,$(UNAME_S)))
  EXE := .exe
else
  EXE :=
endif

# Set OpenMP flags based on platform and compiler
ifeq ($(UNAME_S),Darwin)
  ifeq ($(IS_CLANG),yes)
    OPENMP := -Xpreprocessor -fopenmp
    LDFLAGS += -lomp
  else
    OPENMP := -fopenmp
  endif
else
  ifeq ($(IS_CLANG),yes)
    OPENMP := -fopenmp=libgomp
  else
    OPENMP := -fopenmp
  endif
endif

CFLAGS += $(OPENMP)

# Lemma I: Gram Jacobi construction, Weyl law, Sturm oscillation
PROOF_LEMMA_I = deboor_golub weyl_law_verify phase_shooting stronger_conditions

# Lemma II: Correction formula
PROOF_LEMMA_II = derive_k derive_k2 derive_k_gmp

# Theorem I: Central convergence, trace formula
PROOF_THEOREM_I = trace_verify trace_error_bound heat_kernel_expansion

# Theorem II: Spectral shift convergence, Birman-Krein + Guinand-Weil
PROOF_THEOREM_II = spectral_shift predict sincos_n25 krein_ssf_enhanced

# Theorem III: Spectral determinant identity
PROOF_THEOREM_III = det2_vs_xi det_uniform_bound det2_uniform_bound \
	hadamard_vs_analytic hadamard_extrapolation hadamard_10k \
	hadamard_terminal xi_hadamard_vs_mpmath spectral_rigidity

# Path A/B/C: Birman-Krein, Gaussian-Weil, epsilon closure
PROOF_PATHS = prove_path_a_determinant prove_path_b_gaussian prove_epsilon_zero_closure

# Level 3 barrier closure
PROOF_LEVEL3 = resolvent_trace

# Gap closures: P0 (alpha outliers), P1-3 (Gram!=zeros), k(p), L8
PROOF_GAPS = classify_primes gauge_invariance kp_verify derive_kp l8_verify \
	stationary_phase_test

# Infrastructure and verification
PROOF_INFRA = isospectral_flow phase234 trace_class \
	spectral_cauchy unnormalized_cauchy \
	mfunction_bridge m_function m_zeta_bridge rescaled_operator \
	spectral_chain spectral_chain_extended riemann_siegel_zeta \
	jacobian_analysis block_jacobi \
	analytic_detrend analytic_entry_formula det_xi_match \
	wkb_resolvent

PROOF_CORE = $(PROOF_LEMMA_I) $(PROOF_LEMMA_II) \
	$(PROOF_THEOREM_I) $(PROOF_THEOREM_II) $(PROOF_THEOREM_III) \
	$(PROOF_PATHS) $(PROOF_LEVEL3) $(PROOF_GAPS) $(PROOF_INFRA)

BINS = $(foreach p,$(PROOF_CORE),$(BINDIR)/$p$(EXE))

all: $(BINDIR) $(BINS)

$(BINDIR):
	@mkdir -p $(BINDIR)

$(BINDIR)/%$(EXE): $(SRCDIR)/%.c | $(BINDIR)
	@echo "  $*"
	@$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

test: all
	@for prog in $(PROOF_CORE); do \
		echo "=== $$prog ==="; \
		$(BINDIR)/$$prog$(EXE) 2>&1 | head -15; \
		echo ""; \
	done

clean:
	rm -rf $(BINDIR)

.PHONY: all test clean
