#!/usr/bin/env bash
# ==============================================================================
# Master Reproducibility Script for Rolling Cosine-Similarity Screening
# Evaluates approximate matching against Myers' bit-vector and Bitap
# ==============================================================================
set -e

# Move to repository root
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

echo "=== [1/3] Compiling Native C Benchmark Executables ==="
gcc -O3 src/run_precise_bench.c -lm -o run_precise_bench
gcc -O3 src/benchmark_suite.c -lm -o benchmark_suite
echo "C benchmark binaries built successfully."

echo ""
echo "=== [2/3] Running Approximate Matching Benchmark on E. coli & English ==="
./run_precise_bench

echo ""
echo "=== [3/3] Generating Publication Figures ==="
python3 scripts/generate_figures.py

echo ""
echo "=============================================================================="
echo "All benchmarks completed successfully and publication figures generated!"
echo "=============================================================================="
