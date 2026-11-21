import matplotlib.pyplot as plt
import numpy as np

# Set standard publication font styles
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'legend.fontsize': 9,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'figure.autolayout': True
})

# -------------------------------------------------------------
# Figure 1: Time vs Alphabet Size sigma
# -------------------------------------------------------------
sigmas = [2, 4, 8, 16, 32, 64, 128, 250]
bf_times = [6.66, 5.34, 3.21, 2.06, 1.49, 1.22, 1.25, 1.20]
kmp_times = [6.88, 5.52, 4.21, 3.31, 2.79, 2.71, 2.86, 2.62]
ours_times = [3.18, 3.07, 2.61, 2.67, 2.42, 2.24, 2.66, 2.92]
kr_times = [1.25, 1.38, 1.35, 1.24, 1.32, 1.31, 1.45, 1.27]
bndm_times = [1.40, 0.90, 0.66, 0.64, 0.63, 0.43, 0.39, 0.32]

fig, ax = plt.subplots(figsize=(6.5, 3.8))
ax.plot(sigmas, bf_times, 'b-o', label='Brute Force (BF)', linewidth=1.8, markersize=5)
ax.plot(sigmas, kmp_times, 'orange', linestyle='-', marker='s', label='Knuth-Morris-Pratt (KMP)', linewidth=1.8, markersize=5)
ax.plot(sigmas, ours_times, 'g-^', label='Ours (Fixed q=2)', linewidth=2.2, markersize=6)
ax.plot(sigmas, kr_times, 'gold', linestyle='--', marker='d', label='Rabin-Karp (KR)', linewidth=1.8, markersize=5)
ax.plot(sigmas, bndm_times, 'm-.', marker='p', label='BNDM (Bit-Parallel)', linewidth=1.8, markersize=5)

ax.set_xscale('log', base=2)
ax.set_xticks(sigmas)
ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
ax.set_xlabel(r'Alphabet Size $\sigma$ (log scale)')
ax.set_ylabel('Execution Time (ms)')
ax.set_title(r'Exact Matching Time vs. Alphabet Size $\sigma$ ($m=16$, $n=10^6$ bytes)')
ax.grid(True, which='both', linestyle=':', alpha=0.6)
ax.legend(loc='upper right', frameon=True, ncol=2)
fig.savefig('fig_alphasweep.pdf', dpi=300)
print('Saved fig_alphasweep.pdf')

# -------------------------------------------------------------
# Figure 2: Approximate Matching Runtime across m and k on DNA
# -------------------------------------------------------------
fig2, ax2 = plt.subplots(figsize=(6.5, 3.8))
labels = ['m=16, k=1', 'm=16, k=2', 'm=24, k=1', 'm=24, k=2', 'm=32, k=1', 'm=32, k=2', 'm=32, k=3']
ours_approx = [1.61, 3.13, 1.56, 1.56, 1.56, 1.55, 1.56]
myers_approx = [2.67, 2.61, 2.61, 2.60, 2.62, 2.60, 2.63]
bitap_approx = [1.35, 2.18, 1.32, 2.20, 1.32, 2.17, 2.43]

x = np.arange(len(labels))
width = 0.25

ax2.bar(x - width, ours_approx, width, label='Ours (Cosine Screen + Banded DP)', color='#2ca02c')
ax2.bar(x, myers_approx, width, label='Myers Bit-Vector (1999)', color='#1f77b4')
ax2.bar(x + width, bitap_approx, width, label='Bitap / Shift-Or (Wu-Manber)', color='#d62728')

ax2.set_ylabel('Execution Time (ms)')
ax2.set_title(r'Approximate String Matching Runtime on Real DNA ($E.\ coli$, $n=10^6$)')
ax2.set_xticks(x)
ax2.set_xticklabels(labels, rotation=25, ha='right')
ax2.legend(loc='upper right', frameon=True)
ax2.grid(True, axis='y', linestyle=':', alpha=0.6)
fig2.savefig('fig_approx_bench.pdf', dpi=300)
print('Saved fig_approx_bench.pdf')
