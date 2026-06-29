/*
 * Generates golden reference data for JPL DE file header parsing.
 * Standalone parser — reads the binary format directly without linking
 * against the Swiss Ephemeris library. Mirrors the logic in swejpl.c
 * fsizer() and state() first-call init.
 *
 * Compile (from tests/c-gen/):
 *   cc -Wall -o gen_jpl_header gen_jpl_header.c -lm
 * Run (from tests/c-gen/):
 *   ./gen_jpl_header > ../golden-data/jpl_header.json
 *
 * Reads: ../../ephe/de441.eph
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

static void swap_bytes(void *x, int size) {
    char *p = (char *)x;
    char tmp;
    for (int i = 0; i < size / 2; i++) {
        tmp = p[i];
        p[i] = p[size - 1 - i];
        p[size - 1 - i] = tmp;
    }
}

static void reorder_doubles(double *arr, int n, int do_reorder) {
    if (!do_reorder) return;
    for (int i = 0; i < n; i++)
        swap_bytes(&arr[i], 8);
}

static void reorder_i32s(int32_t *arr, int n, int do_reorder) {
    if (!do_reorder) return;
    for (int i = 0; i < n; i++)
        swap_bytes(&arr[i], 4);
}

int main(void) {
    const char *path = "../../ephe/de441.eph";
    FILE *fp = fopen(path, "rb");
    if (!fp) {
        fprintf(stderr, "Cannot open %s\n", path);
        return 1;
    }

    /* Skip 252-byte title + 2400-byte constant names to reach offset 2652 */
    fseek(fp, 2652, SEEK_SET);

    /* Read ss[3] as raw bytes first, then detect byte order */
    double ss[3];
    fread(ss, sizeof(double), 3, fp);

    /* If ss[2] (segment length in days) is implausible, file needs byte-swap */
    int do_reorder = (ss[2] < 1.0 || ss[2] > 200.0) ? 1 : 0;
    reorder_doubles(ss, 3, do_reorder);

    if (ss[2] < 1.0 || ss[2] > 200.0) {
        fprintf(stderr, "Cannot detect byte order: ss[2]=%g\n", ss[2]);
        fclose(fp);
        return 1;
    }

    int32_t ncon;
    fread(&ncon, sizeof(int32_t), 1, fp);
    reorder_i32s(&ncon, 1, do_reorder);

    double au;
    fread(&au, sizeof(double), 1, fp);
    reorder_doubles(&au, 1, do_reorder);

    double emrat;
    fread(&emrat, sizeof(double), 1, fp);
    reorder_doubles(&emrat, 1, do_reorder);

    int32_t ipt[39];
    fread(ipt, sizeof(int32_t), 36, fp);
    reorder_i32s(ipt, 36, do_reorder);

    int32_t numde;
    fread(&numde, sizeof(int32_t), 1, fp);
    reorder_i32s(&numde, 1, do_reorder);

    /* lpt[3] -> ipt[36..38] */
    fread(&ipt[36], sizeof(int32_t), 3, fp);
    reorder_i32s(&ipt[36], 3, do_reorder);

    fclose(fp);

    /* Compute ksize from ipt[] (swejpl.c:275-291) */
    int32_t kmx = 0;
    int khi = 0;
    for (int i = 0; i < 13; i++) {
        if (ipt[i * 3] > kmx) {
            kmx = ipt[i * 3];
            khi = i + 1; /* 1-based */
        }
    }
    int nd = (khi == 12) ? 2 : 3;
    int32_t ksize = (ipt[khi*3-3] + nd * ipt[khi*3-2] * ipt[khi*3-1] - 1) * 2;
    if (ksize == 1546) ksize = 1652; /* DE102 padding */
    int ncoeffs = ksize / 2;

    /* Emit JSON */
    printf("{\n");
    printf("  \"ss\": [%.20e, %.20e, %.20e],\n", ss[0], ss[1], ss[2]);
    printf("  \"au\": %.20e,\n", au);
    printf("  \"emrat\": %.20e,\n", emrat);
    printf("  \"denum\": %d,\n", numde);
    printf("  \"ncon\": %d,\n", ncon);
    printf("  \"ipt\": [");
    for (int i = 0; i < 39; i++) {
        printf("%d", ipt[i]);
        if (i < 38) printf(", ");
    }
    printf("],\n");
    printf("  \"ksize\": %d,\n", ksize);
    printf("  \"ncoeffs\": %d\n", ncoeffs);
    printf("}\n");

    return 0;
}
