#ifndef YAK_COUNT_H
#define YAK_COUNT_H

#include "yak-priv.h"

#define YC_COUNT_CONTIG 1
#define YC_COUNT_READ 2
#define YC_COUNT_CONTIG_ALLOW_BF 4

yak_ch_v *yak_count(const char *fn, const yak_copt_t *copt, yak_opt_t *opt,
					yak_ch_t *h_asm, yak_ch_t *h_reads, yak_ch_t *h_reads_orphan, uint8_t mode);

#endif  // YAK_COUNT_H