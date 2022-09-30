#ifndef YAM_COMPARE_H
#define YAM_COMPARE_H
#include "yak-priv.h"
int main_katcompare_count_mode(const char *fn_asm, const char *fn_reads, yak_copt_t *copt);
int main_katcompare_validation_mode(const char *fn_asm_ht, const char *fn_reads, const char *fn_asm_fa, 
                                    yak_copt_t *copt, int read_low, int read_high, int asm_low, int asm_high);
int main_dist_via_sketch(const char *fn, const char *fn2, int fl_mask, int k, int n_hash);
#endif  // YAM_COMPARE_H