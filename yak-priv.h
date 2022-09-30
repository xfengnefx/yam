#ifndef YAK_PRIV_H
#define YAK_PRIV_H

#include <stdint.h>
#include "kvec.h"

#define CALLOC(ptr, len) ((ptr) = (__typeof__(ptr))calloc((len), sizeof(*(ptr))))
#define MALLOC(ptr, len) ((ptr) = (__typeof__(ptr))malloc((len) * sizeof(*(ptr))))
#define REALLOC(ptr, len) ((ptr) = (__typeof__(ptr))realloc((ptr), (len) * sizeof(*(ptr))))
typedef kvec_t(uint32_t) vec_u32t;
typedef kvec_t(uint64_t) vec_u64t;
typedef kvec_t(vec_u64t) vec_u64v;
typedef kvec_t(char) vec_char;

#define YAK_MAX_KMER     31
#define YAK_COUNTER_BITS 10
#define YAK_N_COUNTS     (1<<YAK_COUNTER_BITS)
#define YAK_MAX_COUNT    ((1<<YAK_COUNTER_BITS)-1)

#define YAK_BLK_SHIFT  9 // 64 bytes, the size of a cache line
#define YAK_BLK_MASK   ((1<<(YAK_BLK_SHIFT)) - 1)

double Get_T(void);

extern unsigned char seq_nt4_table[256];

static inline uint64_t yak_hash64(uint64_t key, uint64_t mask) // invertible integer hash function
{
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}

static inline uint64_t yak_hash64_64(uint64_t key)
{
	key = ~key + (key << 21);
	key = key ^ key >> 24;
	key = (key + (key << 3)) + (key << 8);
	key = key ^ key >> 14;
	key = (key + (key << 2)) + (key << 4);
	key = key ^ key >> 28;
	key = key + (key << 31);
	return key;
}

static inline uint64_t yak_hash_long(uint64_t x[4])
{
	int j = x[1] < x[3]? 0 : 1;
	return yak_hash64_64(x[j<<1|0]) + yak_hash64_64(x[j<<1|1]);
}


// The inversion of hash_64(). Modified from <https://naml.us/blog/tag/invertible>
// https://gist.github.com/lh3/974ced188be2f90422cc
static inline uint64_t yak_hash64i(uint64_t key, uint64_t mask)
{
	uint64_t tmp;

	// Invert key = key + (key << 31)
	tmp = (key - (key << 31));
	key = (key - (tmp << 31)) & mask;

	// Invert key = key ^ (key >> 28)
	tmp = key ^ key >> 28;
	key = key ^ tmp >> 28;

	// Invert key *= 21
	key = (key * 14933078535860113213ull) & mask;

	// Invert key = key ^ (key >> 14)
	tmp = key ^ key >> 14;
	tmp = key ^ tmp >> 14;
	tmp = key ^ tmp >> 14;
	key = key ^ tmp >> 14;

	// Invert key *= 265
	key = (key * 15244667743933553977ull) & mask;

	// Invert key = key ^ (key >> 24)
	tmp = key ^ key >> 24;
	key = key ^ tmp >> 24;

	// Invert key = (~key) + (key << 21)
	tmp = ~key;
	tmp = ~(key - (tmp << 21));
	tmp = ~(key - (tmp << 21));
	key = ~(key - (tmp << 21)) & mask;

	return key;
}


typedef struct {
	int n_shift, n_hashes;
	uint8_t *b;
} yak_bf_t;  // hashtable (bloom filter)
typedef struct {
	struct yak_ht_t *h;
	yak_bf_t *b;
} yak_ch1_t;  // hashtable (companion)
typedef struct {
	int k, pre, n_hash, n_shift;
	uint64_t tot;
	yak_ch1_t *h;
} yak_ch_t;  // hashtable options
typedef struct {
	yak_ch_t **hh;
	int n;
}yak_ch_v;  // just a convinience for returning multiple hashtables. 
void yak_ch_destroy(yak_ch_t *h);
void yak_ch_destroy_bf(yak_ch_t *h);

typedef struct {
	int32_t bf_shift, bf_n_hash;
	int32_t k;
	int32_t pre;
	int32_t n_thread;
	int64_t chunk_size;
	char *output_name;
	yak_ch_t *h_ref;  //  optional, hastable from assemly 
	// special options for `compare`
	int is_dump_kmer_loc;
	int is_dump_zero_fasta;
	
} yak_copt_t;  // counter options


typedef struct {
	float read_hit_ratio;  // how much kmers shall a read have to be classified not orphan 
	
} yak_opt_t;  // global options

void yak_copt_init(yak_copt_t *o);
void yak_opt_init(yak_opt_t *o);

#endif
