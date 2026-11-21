# Rolling-Cosine-Screening

A C implementation and benchmarking suite for **Rolling Cosine-Similarity Screening** for exact and approximate string matching.

## Overview

Classical filter-and-verify string matching (e.g., Rabin--Karp) uses an algebraic rolling hash. While efficient, binary hashes cannot quantify partial similarity and therefore cannot extend to approximate string matching.

This repository implements a rolling filter over $q$-gram frequency vectors maintained in $O(1)$ amortized time per shift:
- **Exact Matching:** Eliminates all floating-point operations using the Cauchy--Schwarz identity $(A \cdot B)^2 = \|A\|^2 \|B\|^2$.
- **Approximate Matching:** Derives a lossless threshold from the $q$-gram lemma:
  $$A \cdot B \ge (m - q + 1) - k \cdot q$$
  guaranteeing zero false negatives (100% recall) under edit distance $k$.

## Repository Structure

```text
.
├── src/
│   ├── cosinescreen.c       # Fixed q=2 exact screening algorithm (SMART integration)
│   ├── cosinescreendqlin.c   # Dynamic q exact screening algorithm
│   ├── run_precise_bench.c  # Approximate benchmark (Ours vs. Myers Bit-Vector vs. Bitap)
│   └── benchmark_suite.c    # Extended experimental test harness
├── scripts/
│   ├── run_experiments.sh   # Master execution and verification script
│   └── generate_figures.py  # Publication figure plotting script
├── .gitignore
└── README.md
```

## Quick Start

### Prerequisites
- GCC (supporting C99 or later)
- POSIX-compliant environment (Linux / macOS)
- Python 3 with `matplotlib` and `numpy` (for plotting)

### Compiling and Running
To build the binaries and run the approximate string matching benchmarks:

```bash
# Clone the repository
git clone https://github.com/Sai0Nikhil/Rolling-Cosine-Screening.git
cd Rolling-Cosine-Screening

# Compile and run benchmarks
./scripts/run_experiments.sh
```

Or compile manually:
```bash
gcc -O3 src/run_precise_bench.c -lm -o run_precise_bench
./run_precise_bench
```

## Citation

If you use this codebase or benchmark suite in your research, please cite:

```bibtex
@article{nikhil2026cosinescreen,
  author    = {Sai Nikhil},
  title     = {Rolling Cosine-Similarity Screening for Exact and Approximate String Matching},
  year      = {2026}
}
```

## License

This project is licensed under the MIT License.
