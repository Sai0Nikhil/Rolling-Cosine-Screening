#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

#define MAX_M 64
#define HASH_TABLE_SIZE 65536

static inline double get_time_ms(struct timespec start, struct timespec end) {
    return (end.tv_sec - start.tv_sec) * 1000.0 + (end.tv_nsec - start.tv_nsec) / 1000000.0;
}

static inline int edit_distance_le_k(const unsigned char *P, int m, const unsigned char *W, int len, int k) {
    int D[MAX_M + 1][2 * MAX_M + 1];
    for (int i = 0; i <= m; i++) D[i][0] = i;
    for (int j = 0; j <= len; j++) D[0][j] = 0;

    int min_last = m + 1;
    for (int j = 1; j <= len; j++) {
        for (int i = 1; i <= m; i++) {
            int cost = (P[i - 1] == W[j - 1]) ? 0 : 1;
            int del = D[i][j - 1] + 1;
            int ins = D[i - 1][j] + 1;
            int sub = D[i - 1][j - 1] + cost;
            int v = del < ins ? del : ins;
            if (sub < v) v = sub;
            D[i][j] = v;
        }
        if (D[m][j] < min_last) min_last = D[m][j];
    }
    return min_last <= k;
}

static inline int verify_candidate(const unsigned char *P, int m, const unsigned char *T, int n, int end_pos, int k) {
    int start = end_pos - m - k + 1;
    if (start < 0) start = 0;
    int len = (end_pos + k) - start + 1;
    if (start + len > n) len = n - start;
    return edit_distance_le_k(P, m, T + start, len, k);
}

int myers_search(const unsigned char *P, int m, const unsigned char *T, int n, int k) {
    if (m > 64) return 0;
    uint64_t B[256];
    for (int c = 0; c < 256; c++) B[c] = 0;
    for (int i = 0; i < m; i++) B[P[i]] |= (1ULL << i);

    uint64_t VP = ~0ULL;
    uint64_t VN = 0;
    int score = m;
    int matches = 0;

    for (int j = 0; j < n; j++) {
        uint64_t X = B[T[j]] | VN;
        uint64_t D0 = (((X & VP) + VP) ^ VP) | X;
        uint64_t HN = VP & D0;
        uint64_t HP = VN | ~(VP | D0);
        X = (HP << 1) | 1ULL;
        VN = X & D0;
        VP = (HN << 1) | ~(X | D0);
        if (HP & (1ULL << (m - 1))) score++;
        else if (HN & (1ULL << (m - 1))) score--;
        if (score <= k) matches++;
    }
    return matches;
}

int bitap_search(const unsigned char *P, int m, const unsigned char *T, int n, int k) {
    if (m > 64) return 0;
    uint64_t pattern_mask[256];
    for (int i = 0; i < 256; i++) pattern_mask[i] = ~0ULL;
    for (int i = 0; i < m; i++) pattern_mask[P[i]] &= ~(1ULL << i);

    uint64_t R[MAX_M + 1];
    for (int i = 0; i <= k; i++) R[i] = ~1ULL;

    int matches = 0;
    uint64_t match_bit = 1ULL << (m - 1);

    for (int j = 0; j < n; j++) {
        uint64_t old_R0 = R[0];
        R[0] = ((R[0] << 1) | 1ULL) | pattern_mask[T[j]];

        for (int d = 1; d <= k; d++) {
            uint64_t sub = (old_R0 << 1) | 1ULL;
            uint64_t ins = old_R0;
            uint64_t del = (R[d] << 1) | 1ULL;
            uint64_t mat = ((R[d] << 1) | 1ULL) | pattern_mask[T[j]];

            old_R0 = R[d];
            R[d] = sub & ins & del & mat;
        }
        if ((R[k] & match_bit) == 0) matches++;
    }
    return matches;
}

int choose_approx_q(int m, int k, int sigma) {
    int max_q = (m - 1) / (k + 1);
    if (max_q < 2) return 2;
    int q = max_q;
    if (sigma <= 4) {
        if (q > 5) q = 5;
    } else {
        if (q > 4) q = 4;
    }
    return q < 2 ? 2 : q;
}

int cosine_approx_search(const unsigned char *P, int m, const unsigned char *T, int n, int k, int q,
                         int *candidates_screened, int *verified_matches) {
    int min_shared = (m - q + 1) - k * q;
    if (min_shared < 1) min_shared = 1;

    static int patternCounts[HASH_TABLE_SIZE];
    static int windowCounts[HASH_TABLE_SIZE];

    uint64_t BASE = 131ULL;
    uint64_t power = 1;
    for (int i = 0; i < q - 1; i++) power *= BASE;

    long patternMagSq = 0;
    uint64_t ph = 0;
    for (int i = 0; i < q; i++) ph = ph * BASE + P[i];
    int pkey = ph % HASH_TABLE_SIZE;
    patternCounts[pkey]++;
    patternMagSq += 1;

    for (int i = q; i < m; i++) {
        ph = (ph - (uint64_t)P[i - q] * power) * BASE + P[i];
        pkey = ph % HASH_TABLE_SIZE;
        int old = patternCounts[pkey];
        patternMagSq += (2L * old + 1);
        patternCounts[pkey] = old + 1;
    }

    long dot = 0, windowMagSq = 0;
    uint64_t wh = 0;
    for (int i = 0; i < q; i++) wh = wh * BASE + T[i];
    int wkey = wh % HASH_TABLE_SIZE;
    dot += patternCounts[wkey];
    windowMagSq += 1;
    windowCounts[wkey]++;

    for (int i = q; i < m; i++) {
        wh = (wh - (uint64_t)T[i - q] * power) * BASE + T[i];
        wkey = wh % HASH_TABLE_SIZE;
        int oldW = windowCounts[wkey];
        dot += patternCounts[wkey];
        windowMagSq += (2L * oldW + 1);
        windowCounts[wkey] = oldW + 1;
    }

    uint64_t leavingHash = 0;
    for (int i = 0; i < q; i++) leavingHash = leavingHash * BASE + T[i];
    uint64_t enteringHash = wh;

    int screened = 0;
    int confirmed = 0;

    for (int j = 0; j <= n - m; j++) {
        if (dot >= min_shared) {
            screened++;
            if (verify_candidate(P, m, T, n, j + m - 1, k)) {
                confirmed++;
            }
        }

        if (j < n - m) {
            int lk = (int)(leavingHash % HASH_TABLE_SIZE);
            int oldL = windowCounts[lk];
            dot -= patternCounts[lk];
            windowMagSq -= (2L * oldL - 1);
            windowCounts[lk] = oldL - 1;

            enteringHash = (enteringHash - (uint64_t)T[j + m - q] * power) * BASE + T[j + m];
            int ek = (int)(enteringHash % HASH_TABLE_SIZE);
            int oldE = windowCounts[ek];
            dot += patternCounts[ek];
            windowMagSq += (2L * oldE + 1);
            windowCounts[ek] = oldE + 1;

            leavingHash = (leavingHash - (uint64_t)T[j] * power) * BASE + T[j + q];
        }
    }

    ph = 0;
    for (int i = 0; i < q; i++) ph = ph * BASE + P[i];
    patternCounts[ph % HASH_TABLE_SIZE] = 0;
    for (int i = q; i < m; i++) {
        ph = (ph - (uint64_t)P[i - q] * power) * BASE + P[i];
        patternCounts[ph % HASH_TABLE_SIZE] = 0;
    }

    wh = 0;
    int lastStart = n - m;
    for (int i = 0; i < q; i++) wh = wh * BASE + T[lastStart + i];
    windowCounts[wh % HASH_TABLE_SIZE] = 0;
    for (int i = q; i < m; i++) {
        wh = (wh - (uint64_t)T[lastStart + i - q] * power) * BASE + T[lastStart + i];
        windowCounts[wh % HASH_TABLE_SIZE] = 0;
    }

    *candidates_screened = screened;
    *verified_matches = confirmed;
    return confirmed;
}

int main() {
    const char *datasets[] = {"smart/data/genome/ecoli.txt", "smart/data/englishTexts/bible.txt"};
    const char *names[] = {"Genome (E. coli, DNA)", "English (Bible, text)"};
    int sigmas[] = {4, 64};

    for (int d = 0; d < 2; d++) {
        FILE *fp = fopen(datasets[d], "rb");
        if (!fp) continue;
        fseek(fp, 0, SEEK_END);
        long n_total = ftell(fp);
        fseek(fp, 0, SEEK_SET);

        int n = (n_total > 1000000) ? 1000000 : (int)n_total;
        unsigned char *T = (unsigned char *)malloc(n + 1);
        size_t rd = fread(T, 1, n, fp);
        fclose(fp);

        printf("\n=====================================================================================\n");
        printf("DATASET: %s (n = %d bytes)\n", names[d], n);
        printf("=====================================================================================\n");

        int test_m[] = {16, 24, 32};
        int test_k[] = {1, 2, 3};
        srand(42);

        for (int mi = 0; mi < 3; mi++) {
            int m = test_m[mi];
            for (int ki = 0; ki < 3; ki++) {
                int k = test_k[ki];
                if (k >= m / 3) continue;

                int q = choose_approx_q(m, k, sigmas[d]);
                int num_patterns = 30;

                double time_cos = 0, time_myers = 0, time_bitap = 0;
                int total_cand = 0, total_confirmed = 0;

                for (int p = 0; p < num_patterns; p++) {
                    int r = rand() % (n - m - 10);
                    unsigned char P[MAX_M + 1];
                    memcpy(P, T + r, m);

                    // Insert k errors (substitutions)
                    for (int e = 0; e < k; e++) {
                        int pos = (e * (m / k) + 2) % m;
                        P[pos] = (P[pos] == 'A') ? 'T' : 'A';
                    }

                    struct timespec ts0, ts1;

                    // 1. Cosine
                    int c_cand = 0, c_conf = 0;
                    clock_gettime(CLOCK_MONOTONIC, &ts0);
                    cosine_approx_search(P, m, T, n, k, q, &c_cand, &c_conf);
                    clock_gettime(CLOCK_MONOTONIC, &ts1);
                    time_cos += get_time_ms(ts0, ts1);
                    total_cand += c_cand;
                    total_confirmed += c_conf;

                    // 2. Myers
                    clock_gettime(CLOCK_MONOTONIC, &ts0);
                    volatile int my_res = myers_search(P, m, T, n, k); (void)my_res;
                    clock_gettime(CLOCK_MONOTONIC, &ts1);
                    time_myers += get_time_ms(ts0, ts1);

                    // 3. Bitap
                    clock_gettime(CLOCK_MONOTONIC, &ts0);
                    volatile int bi_res = bitap_search(P, m, T, n, k); (void)bi_res;
                    clock_gettime(CLOCK_MONOTONIC, &ts1);
                    time_bitap += get_time_ms(ts0, ts1);
                }

                double avg_cos = time_cos / num_patterns;
                double avg_myers = time_myers / num_patterns;
                double avg_bitap = time_bitap / num_patterns;

                double pruning = 100.0 * (1.0 - (double)total_cand / ((double)(n - m + 1) * num_patterns));
                double precision = (total_cand > 0) ? (100.0 * (double)total_confirmed / total_cand) : 100.0;

                printf("m=%2d, k=%d, q=%d | Pruning: %6.2f%% | Prec: %5.2f%% | Times: Ours=%5.2f ms, Myers=%5.2f ms, Bitap=%5.2f ms\n",
                       m, k, q, pruning, precision, avg_cos, avg_myers, avg_bitap);
            }
        }
        free(T);
    }
    return 0;
}
