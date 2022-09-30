#ifndef YAK_HASHTABLE_H
#define YAK_HASHTABLE_H
#include "yak-priv.h"

typedef struct yak_ht_t yak_ht_t;
void radix_sort_u64t(uint64_t *beg, uint64_t *end);

yak_ch_t *yak_ch_restore(const char *fn);
int yak_ch_dump(const yak_ch_t *ch, const char *fn);
void yam_ht_destroy(yak_ht_t *h);

void yam_compare_count_2D(yak_ch_t *h_asm, yak_ch_t *h_reads, int **mat, int pre);
void yam_compare_zero_dump(yak_ch_t *h_asm, yak_ch_t *h_reads, int kmersize, int pre);
void yam_compare_count_2D_blue1mismatch(yak_ch_t *h_asm, yak_ch_t *h_reads, 
										int *buf, int pre, int kmersize);

yak_ht_t *yam_compare_get_kmers_of_interest(yak_ch_t *h_asm, yak_ch_t *h_reads, int pre,
											int read_low, int read_high, int asm_low, int asm_high);
int yam_compare_count_hits_given_seq_and_ht(yak_ht_t *h, char *seq, int seq_l, int32_t k);
int yam_compare_count_bed_given_seq_and_ht(FILE *fp, yak_ht_t *h, 
										  char *seq, int seq_l, 
										  char *seqname, int seqname_l, 
										  int32_t k, int l_gap);

vec_u64t yam_minhash_sketch(char *seq, int seq_l, int k, int n_hash);
double yam_minhash_mashdist_core(vec_u64t *h1, vec_u64t *h2, int k,  int n_hash, int *ret_share);

#endif  // YAK_HASHTABLE_H