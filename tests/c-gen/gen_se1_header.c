/*
 * Generates golden reference data for SE1 file header parsing.
 * Standalone parser — reads the binary format directly without linking
 * against the Swiss Ephemeris library.
 *
 * Compile:
 *   cc -Wall -o gen_se1_header gen_se1_header.c -lm
 * Run:
 *   ./gen_se1_header > ../golden-data/se1_header.json
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define ENDIAN_TEST_VAL 0x616263
#define SEI_FLG_ELLIPSE 4
#define SE_MARS 4

static const char *basename_of(const char *path) {
    const char *p = strrchr(path, '/');
    return p ? p + 1 : path;
}

static void parse_se1(const char *path, const char *file_type_str, int is_last) {
    FILE *fp = fopen(path, "rb");
    if (!fp) {
        fprintf(stderr, "Cannot open %s\n", path);
        exit(1);
    }

    /* --- Text header: 3 lines for planet/moon --- */
    char line[512];
    int version = 0;
    for (int i = 0; i < 3; i++) {
        if (!fgets(line, sizeof(line), fp)) {
            fprintf(stderr, "Failed to read text header line %d\n", i + 1);
            exit(1);
        }
        if (i == 0) {
            char *p = line;
            while (*p && !(*p >= '0' && *p <= '9')) p++;
            version = atoi(p);
        }
    }

    /* --- Binary header --- */
    /* Byte order detection */
    uint32_t endian_raw;
    fread(&endian_raw, 4, 1, fp);
    int needs_swap = 0;
    const char *byte_order_str;
    if (endian_raw == ENDIAN_TEST_VAL) {
        needs_swap = 0;
        /* Host and file agree — detect which */
        unsigned char *c = (unsigned char *)&endian_raw;
        if (*c == 0x63)
            byte_order_str = "little";
        else
            byte_order_str = "big";
    } else {
        needs_swap = 1;
        unsigned char *c = (unsigned char *)&endian_raw;
        if (*c == 0x63)
            byte_order_str = "big";      /* first byte is LSB but file needs swap → file is BE */
        else
            byte_order_str = "little";
        /* Verify swapped value */
        uint32_t swapped = ((endian_raw >> 24) & 0xff)
                         | ((endian_raw >> 8) & 0xff00)
                         | ((endian_raw << 8) & 0xff0000)
                         | ((endian_raw << 24) & 0xff000000u);
        if (swapped != ENDIAN_TEST_VAL) {
            fprintf(stderr, "Invalid endian test value in %s\n", path);
            exit(1);
        }
    }
    /* For simplicity, assert no swap needed (all standard files are LE on LE hosts) */
    if (needs_swap) {
        fprintf(stderr, "Byte-swapped files not supported by this harness\n");
        exit(1);
    }

    int32_t file_len;
    fread(&file_len, 4, 1, fp);

    int32_t denum;
    fread(&denum, 4, 1, fp);

    double tfstart, tfend;
    fread(&tfstart, 8, 1, fp);
    fread(&tfend, 8, 1, fp);

    int16_t nplan_raw;
    fread(&nplan_raw, 2, 1, fp);
    int nbytes_ipl = 2;
    int nplan = nplan_raw;
    if (nplan > 256) {
        nbytes_ipl = 4;
        nplan %= 256;
    }

    int32_t ipl[50];
    for (int i = 0; i < nplan; i++) {
        if (nbytes_ipl == 4) {
            fread(&ipl[i], 4, 1, fp);
        } else {
            int16_t id16;
            fread(&id16, 2, 1, fp);
            ipl[i] = id16;
        }
    }

    /* CRC32 — skip */
    fseek(fp, 4, SEEK_CUR);

    /* Physical constants — skip 5 doubles */
    fseek(fp, 40, SEEK_CUR);

    /* --- Output --- */
    printf("  {\n");
    printf("    \"filename\": \"%s\",\n", basename_of(path));
    printf("    \"version\": %d,\n", version);
    printf("    \"file_type\": \"%s\",\n", file_type_str);
    printf("    \"tfstart\": %.20e,\n", tfstart);
    printf("    \"tfend\": %.20e,\n", tfend);
    printf("    \"denum\": %d,\n", denum);
    printf("    \"byte_order\": \"%s\",\n", byte_order_str);
    printf("    \"planets\": [\n");

    for (int k = 0; k < nplan; k++) {
        int32_t lndx0;
        fread(&lndx0, 4, 1, fp);

        uint8_t iflg_byte;
        fread(&iflg_byte, 1, 1, fp);
        uint32_t iflg = iflg_byte;

        uint8_t ncoe_byte;
        fread(&ncoe_byte, 1, 1, fp);
        int ncoe = ncoe_byte;

        int32_t rmax_raw;
        fread(&rmax_raw, 4, 1, fp);
        double rmax = rmax_raw / 1000.0;

        double orbital[10];
        fread(orbital, 8, 10, fp);

        double p_tfstart = orbital[0];
        double p_tfend   = orbital[1];
        double dseg      = orbital[2];
        double telem     = orbital[3];
        double prot      = orbital[4];
        double dprot     = orbital[5];
        double qrot      = orbital[6];
        double dqrot     = orbital[7];
        double peri      = orbital[8];
        double dperi     = orbital[9];

        int32_t nndx = (int32_t)((p_tfend - p_tfstart + 0.1) / dseg);

        printf("      {\n");
        printf("        \"body_id\": %d,\n", ipl[k]);
        printf("        \"iflg\": %u,\n", iflg);
        printf("        \"ncoe\": %d,\n", ncoe);
        printf("        \"neval\": %d,\n", ncoe);
        printf("        \"rmax\": %.20e,\n", rmax);
        printf("        \"dseg\": %.20e,\n", dseg);
        printf("        \"tfstart\": %.20e,\n", p_tfstart);
        printf("        \"tfend\": %.20e,\n", p_tfend);
        printf("        \"lndx0\": %d,\n", lndx0);
        printf("        \"nndx\": %d,\n", nndx);
        printf("        \"telem\": %.20e,\n", telem);
        printf("        \"prot\": %.20e,\n", prot);
        printf("        \"qrot\": %.20e,\n", qrot);
        printf("        \"dprot\": %.20e,\n", dprot);
        printf("        \"dqrot\": %.20e,\n", dqrot);
        printf("        \"peri\": %.20e,\n", peri);
        printf("        \"dperi\": %.20e\n", dperi);
        printf("      }%s\n", k < nplan - 1 ? "," : "");

        /* Skip refep if SEI_FLG_ELLIPSE */
        if (iflg & SEI_FLG_ELLIPSE) {
            fseek(fp, ncoe * 2 * 8, SEEK_CUR);
        }
    }

    printf("    ]\n");
    printf("  }%s\n", is_last ? "" : ",");

    fclose(fp);
}

int main(void) {
    printf("{\n");
    printf("  \"files\": [\n");
    parse_se1("../../../swisseph/ephe/sepl_18.se1", "planet", 0);
    parse_se1("../../../swisseph/ephe/semo_18.se1", "moon", 1);
    printf("  ]\n");
    printf("}\n");
    return 0;
}
