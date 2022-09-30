#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "kthread.h"
#include "yak-priv.h"
#include "khashl.h" // hash table
#include "ksort.h"

#define YAK_LOAD_ALL       1
#define YAK_LOAD_TRIOBIN1  2
#define YAK_LOAD_TRIOBIN2  3

#define YAK_MAGIC "YAK\2"

#define yak_ch_eq(a, b) ((a)>>YAK_COUNTER_BITS == (b)>>YAK_COUNTER_BITS) // lower 8 bits for counts; higher bits for k-mer
#define yak_ch_hash(a) ((a)>>YAK_COUNTER_BITS)
KHASHL_SET_INIT(, yak_ht_t, yak_ht, uint64_t, yak_ch_hash, yak_ch_eq)
KRADIX_SORT_INIT(u64t, uint64_t, (uint64_t), 8)


yak_bf_t *yak_bf_init(int n_shift, int n_hashes)
{
	yak_bf_t *b;
	void *ptr = 0;
	if (n_shift + YAK_BLK_SHIFT > 64 || n_shift < YAK_BLK_SHIFT) return 0;
	b = calloc(1, sizeof(yak_bf_t));
	b->n_shift = n_shift;
	b->n_hashes = n_hashes;
	posix_memalign(&ptr, 1<<(YAK_BLK_SHIFT-3), 1ULL<<(n_shift-3));
	b->b = ptr;
	bzero(b->b, 1ULL<<(n_shift-3));
	return b;
}

void yak_bf_destroy(yak_bf_t *b)
{
	if (b == 0) return;
	free(b->b); free(b);
}

int yak_bf_insert(yak_bf_t *b, uint64_t hash)
{
	int x = b->n_shift - YAK_BLK_SHIFT;
	uint64_t y = hash & ((1ULL<<x) - 1);
	int h1 = hash >> x & YAK_BLK_MASK;
	int h2 = hash >> b->n_shift & YAK_BLK_MASK;
	uint8_t *p = &b->b[y<<(YAK_BLK_SHIFT-3)];
	int i, z = h1, cnt = 0;
	if ((h2&31) == 0) h2 = (h2 + 1) & YAK_BLK_MASK; // otherwise we may repeatedly use a few bits
	for (i = 0; i < b->n_hashes; z = (z + h2) & YAK_BLK_MASK) {
		uint8_t *q = &p[z>>3], u;
		u = 1<<(z&7);
		cnt += !!(*q & u);
		*q |= u;
		++i;
	}
	return cnt;
}


yak_ch_t *yak_ch_init(int k, int pre, int n_hash, int n_shift)
{
	yak_ch_t *h;
	int i;
	if (pre < YAK_COUNTER_BITS) return 0;
	CALLOC(h, 1);
	h->k = k, h->pre = pre;
	CALLOC(h->h, 1<<h->pre);
	for (i = 0; i < 1<<h->pre; ++i)
		h->h[i].h = yak_ht_init();
	if (n_hash > 0 && n_shift > h->pre) {
		h->n_hash = n_hash, h->n_shift = n_shift;
		for (i = 0; i < 1<<h->pre; ++i)
			h->h[i].b = yak_bf_init(h->n_shift - h->pre, h->n_hash);
	}
	return h;
}

void yak_ch_destroy_bf(yak_ch_t *h)
{
	int i;
	for (i = 0; i < 1<<h->pre; ++i) {
		if (h->h[i].b)
			yak_bf_destroy(h->h[i].b);
		h->h[i].b = 0;
	}
}

void yak_ch_destroy(yak_ch_t *h)
{
	int i;
	if (h == 0) return;
	yak_ch_destroy_bf(h);
	for (i = 0; i < 1<<h->pre; ++i)
		yak_ht_destroy(h->h[i].h);
	free(h->h); free(h);
}

int yak_ch_insert_list(yak_ch_t *h, int create_new, int n, const uint64_t *a)
{
	int j, mask = (1<<h->pre) - 1, n_ins = 0;
	yak_ch1_t *g;
	if (n == 0) return 0;
	g = &h->h[a[0]&mask];
	for (j = 0; j < n; ++j) {
		int ins = 1, absent;
		uint64_t x = a[j] >> h->pre;
		khint_t k;
		if ((a[j]&mask) != (a[0]&mask)) continue;
		if (create_new) {
			if (g->b)
				ins = (yak_bf_insert(g->b, x) == h->n_hash);
			if (ins) {
				k = yak_ht_put(g->h, x<<YAK_COUNTER_BITS, &absent);
				if (absent) ++n_ins;
				if ((kh_key(g->h, k)&YAK_MAX_COUNT) < YAK_MAX_COUNT)
					++kh_key(g->h, k);
			}
		} else {
			k = yak_ht_get(g->h, x<<YAK_COUNTER_BITS);
			if (k != kh_end(g->h) && (kh_key(g->h, k)&YAK_MAX_COUNT) < YAK_MAX_COUNT)
				++kh_key(g->h, k);
		}
	}
	return n_ins;
}

int yak_ch_get(const yak_ch_t *h, uint64_t x)
{
	int mask = (1<<h->pre) - 1;
	yak_ht_t *g = h->h[x&mask].h;
	khint_t k;
	k = yak_ht_get(g, x >> h->pre << YAK_COUNTER_BITS);
	return k == kh_end(g)? -1 : kh_key(g, k)&YAK_MAX_COUNT;
}

/*************************
 * Clear all counts to 0 *
 *************************/

static void worker_clear(void *data, long i, int tid) // callback for kt_for()
{
	yak_ch_t *h = (yak_ch_t*)data;
	yak_ht_t *g = h->h[i].h;
	khint_t k;
	uint64_t mask = ~1ULL >> YAK_COUNTER_BITS << YAK_COUNTER_BITS;
	for (k = 0; k < kh_end(g); ++k)
		if (kh_exist(g, k))
			kh_key(g, k) &= mask;
}

void yak_ch_clear(yak_ch_t *h, int n_thread)
{
	kt_for(n_thread, worker_clear, h, 1<<h->pre);
}

/*************
 * Histogram *
 *************/

typedef struct {
	uint64_t c[YAK_N_COUNTS];
} buf_cnt_t;

typedef struct {
	const yak_ch_t *h;
	buf_cnt_t *cnt;
} hist_aux_t;

static void worker_hist(void *data, long i, int tid) // callback for kt_for()
{
	hist_aux_t *a = (hist_aux_t*)data;
	uint64_t *cnt = a->cnt[tid].c;
	yak_ht_t *g = a->h->h[i].h;
	khint_t k;
	for (k = 0; k < kh_end(g); ++k)
		if (kh_exist(g, k))
			++cnt[kh_key(g, k)&YAK_MAX_COUNT];
}

void yak_ch_hist(const yak_ch_t *h, int64_t cnt[YAK_N_COUNTS], int n_thread)
{
	hist_aux_t a;
	int i, j;
	a.h = h;
	memset(cnt, 0, YAK_N_COUNTS * sizeof(uint64_t));
	CALLOC(a.cnt, n_thread);
	kt_for(n_thread, worker_hist, &a, 1<<h->pre);
	for (i = 0; i < YAK_N_COUNTS; ++i) cnt[i] = 0;
	for (j = 0; j < n_thread; ++j)
		for (i = 0; i < YAK_N_COUNTS; ++i)
			cnt[i] += a.cnt[j].c[i];
	free(a.cnt);
}

/**********
 * Shrink *
 **********/

typedef struct {
	int min, max;
	yak_ch_t *h;
} shrink_aux_t;

static void worker_shrink(void *data, long i, int tid) // callback for kt_for()
{
	shrink_aux_t *a = (shrink_aux_t*)data;
	yak_ch_t *h = a->h;
	yak_ht_t *g = h->h[i].h, *f;
	khint_t k;
	f = yak_ht_init();
	yak_ht_resize(f, kh_size(g));
	for (k = 0; k < kh_end(g); ++k) {
		if (kh_exist(g, k)) {
			int absent, c = kh_key(g, k) & YAK_MAX_COUNT;
			if (c >= a->min && c <= a->max)
				yak_ht_put(f, kh_key(g, k), &absent);
		}
	}
	yak_ht_destroy(g);
	h->h[i].h = f;
}

void yak_ch_shrink(yak_ch_t *h, int min, int max, int n_thread)
{
	int i;
	shrink_aux_t a;
	a.h = h, a.min = min, a.max = max;
	kt_for(n_thread, worker_shrink, &a, 1<<h->pre);
	for (i = 0, h->tot = 0; i < 1<<h->pre; ++i)
		h->tot += kh_size(h->h[i].h);
}

/*******
 * I/O *
 *******/

int yak_ch_dump(const yak_ch_t *ch, const char *fn)
{
	FILE *fp;
	uint32_t t[3];
	int i;
	if ((fp = strcmp(fn, "-")? fopen(fn, "wb") : stdout) == 0) return -1;
	fwrite(YAK_MAGIC, 1, 4, fp);
	t[0] = ch->k, t[1] = ch->pre, t[2] = YAK_COUNTER_BITS;
	fwrite(t, 4, 3, fp);
	for (i = 0; i < 1<<ch->pre; ++i) {
		yak_ht_t *h = ch->h[i].h;
		khint_t k;
		t[0] = kh_capacity(h), t[1] = kh_size(h);
		fwrite(t, 4, 2, fp);
		for (k = 0; k < kh_end(h); ++k)
			if (kh_exist(h, k))
				fwrite(&kh_key(h, k), 8, 1, fp);
	}
	fprintf(stderr, "[M::%s] dumpped the hash table to file '%s'.\n", __func__, fn);
	fclose(fp);
	return 0;
}

yak_ch_t *yak_ch_restore_core(yak_ch_t *ch0, const char *fn, int mode, ...)
{
	va_list ap;
	FILE *fp;
	uint32_t t[3];
	char magic[4];
	int i, j, absent, min_cnt = 0, mid_cnt = 0, mode_err = 0;
	uint64_t mask = (1ULL<<YAK_COUNTER_BITS) - 1, n_ins = 0, n_new = 0;
	yak_ch_t *ch;

	va_start(ap, mode);
	if (mode == YAK_LOAD_ALL) { // do nothing
	} else if (mode == YAK_LOAD_TRIOBIN1 || mode == YAK_LOAD_TRIOBIN2) {
		assert(YAK_COUNTER_BITS >= 4);
		min_cnt = va_arg(ap, int);
		mid_cnt = va_arg(ap, int);
		if (ch0 == 0 && mode == YAK_LOAD_TRIOBIN2)
			mode_err = 1;
	} else mode_err = 1;
	va_end(ap);
	if (mode_err) return 0;

	if ((fp = fopen(fn, "rb")) == 0) return 0;
	if (fread(magic, 1, 4, fp) != 4) return 0;
	if (strncmp(magic, YAK_MAGIC, 4) != 0) {
		fprintf(stderr, "ERROR: wrong file magic.\n");
		fclose(fp);
		return 0;
	}
	fread(t, 4, 3, fp);
	if (t[2] != YAK_COUNTER_BITS) {
		fprintf(stderr, "ERROR: saved counter bits: %d; compile-time counter bits: %d\n", t[2], YAK_COUNTER_BITS);
		fclose(fp);
		return 0;
	}

	ch = ch0 == 0? yak_ch_init(t[0], t[1], 0, 0) : ch0;
	assert((int)t[0] == ch->k && (int)t[1] == ch->pre);
	for (i = 0; i < 1<<ch->pre; ++i) {
		yak_ht_t *h = ch->h[i].h;
		fread(t, 4, 2, fp);
		if (ch0 == 0) yak_ht_resize(h, t[0]);
		for (j = 0; j < t[1]; ++j) {
			uint64_t key;
			fread(&key, 8, 1, fp);
			if (mode == YAK_LOAD_ALL) {
				++n_ins;
				yak_ht_put(h, key, &absent);
				if (absent) ++n_new;
			} else if (mode == YAK_LOAD_TRIOBIN1 || mode == YAK_LOAD_TRIOBIN2) {
				int cnt = key & mask, x, shift = mode == YAK_LOAD_TRIOBIN1? 0 : 2;
				if (cnt >= mid_cnt) x = 2<<shift;
				else if (cnt >= min_cnt) x = 1<<shift;
				else x = -1;
				if (x >= 0) {
					khint_t k;
					key = (key & ~mask) | x;
					++n_ins;
					k = yak_ht_put(h, key, &absent);
					if (absent) ++n_new;
					else kh_key(h, k) = kh_key(h, k) | x;
				}
			}
		}
	}
	fclose(fp);
	fprintf(stderr, "[M::%s] inserted %ld k-mers, of which %ld are new\n", __func__, (long)n_ins, (long)n_new);
	return ch;
}

yak_ch_t *yak_ch_restore(const char *fn)
{
	return yak_ch_restore_core(0, fn, YAK_LOAD_ALL);
}


/*
 other usages
*/

void yam_compare_count_2D(yak_ch_t *h_asm, yak_ch_t *h_reads, int **mat, int pre){
	khint_t k, k2;
	uint64_t mask = (1<<YAK_COUNTER_BITS)-1, xr, xa;
    int i, j, cnt_read, cnt_asm;
    yak_ht_t *h, *h2;

	for (i=0; i<(1<<pre)-1; i++){
		h = h_reads->h[i].h;
		for (k=0; k<kh_end(h); k++){
			if (kh_exist(h, k)){
				cnt_read = kh_key(h, k) & mask;
				xr = kh_key(h, k)>>YAK_COUNTER_BITS<<pre;

				// check if the kmer is present in assembly's hashtable
				h2 = h_asm->h[i].h;
				k2 = yak_ht_get(h2, xr);
				if (k2==kh_end(h2)){
					cnt_asm = 0;
				}else{
					cnt_asm = kh_key(h2, k2) & mask;
				}
				// update counter
				mat[cnt_read-1][cnt_asm]++;
				
			}
		}
    }
}

void yam_compare_zero_dump(yak_ch_t *h_asm, yak_ch_t *h_reads, int kmersize, int pre){
	khint_t k, k2;
	uint64_t mask = (1<<YAK_COUNTER_BITS)-1, xr, xa, kmer, hash;
    uint64_t mask_kmer = (1ULL<<kmersize*2) - 1;
	int i, j, cnt_read, cnt_asm, ii;
    yak_ht_t *h, *h2;
	int *indices = (int*)calloc(YAK_MAX_COUNT+1, sizeof(int));  assert(indices);
	char strmap[4]="ACGT";
	char *kmer_s = (char*)malloc(kmersize);

	for (i=0; i<(1<<pre)-1; i++){
		h = h_reads->h[i].h;
		for (k=0; k<kh_end(h); k++){
			if (kh_exist(h, k)){
				cnt_read = kh_key(h, k) & mask;
				xr = kh_key(h, k)>>YAK_COUNTER_BITS<<pre;

				// check if the kmer is present in assembly's hashtable
				h2 = h_asm->h[i].h;
				k2 = yak_ht_get(h2, xr);
				if (k2==kh_end(h2)){
					cnt_asm = 0;
				}else{
					cnt_asm = kh_key(h2, k2) & mask;
				}
				if (cnt_asm==0){
					fprintf(stdout, ">rcnt%d_%d\n", cnt_read, indices[cnt_read-1]++);
					hash = (kh_key(h, k)>>YAK_COUNTER_BITS<<pre) | (uint64_t)i;
					kmer = yak_hash64i(hash, mask_kmer);
					for (ii=0; ii<kmersize*2; ii+=2){
						kmer_s[kmersize-1-ii/2] = strmap[(kmer>>ii)&3];
					}
					fprintf(stdout, "%.*s\n", kmersize, kmer_s);
				}
			}
		}
    }
	free(indices);
	free(kmer_s);
}


// FUNC
//     For each kmer bucket, find the kmers that are not in the contigs,
//      but could be found if we allow one mismatch.
// PAR
//     length of `buf` is $max_count + 1 
// NOTE
//     too slow not using this; do zerodump and bwa aln instead.
void yam_compare_count_2D_blue1mismatch(yak_ch_t *h_asm, yak_ch_t *h_reads, 
										int *buf, int pre, int kmersize){
	khint_t k, k2;
	int pre2;
	uint64_t mask_pre = (1<<pre)-1;
	uint64_t mask = (1<<YAK_COUNTER_BITS)-1, xr, xa, kmer;
    int i, j, cnt_read, cnt_asm, ik;
    yak_ht_t *h, *h2;

	for (i=0; i<(1<<pre)-1; i++){
		h = h_reads->h[i].h;
		for (k=0; k<kh_end(h); k++){
			if (!kh_exist(h, k)) continue;
			cnt_read = kh_key(h, k) & mask;
			xr = kh_key(h, k)>>YAK_COUNTER_BITS<<pre;

			// check if the kmer is present in assembly's hashtable
			h2 = h_asm->h[i].h;
			k2 = yak_ht_get(h2, xr);
			if (k2==kh_end(h2)){  // not in the contigs, now try with one mismatch
				kmer = yak_hash64i(xr, mask) | (uint64_t)i;
				int max_outcome=0;
				for (ik=0; ik<kmersize; ik+=2){
					xr = xr^(1ULL<<ik);  // mutation 1
					pre2 = xr & mask_pre;
					h2 = h_asm->h[pre2].h;
					k2 = yak_ht_get(h2, xr>>pre<<YAK_COUNTER_BITS);
					if (k2!=kh_end(h2)){
						cnt_asm = kh_key(h2, k2) & YAK_MAX_COUNT;
						max_outcome = cnt_asm>max_outcome? cnt_asm : max_outcome; 
					}
					
					xr = xr^(1ULL<<(ik+1));  // mutation 2
					pre2 = xr & mask_pre;
					h2 = h_asm->h[pre2].h;
					k2 = yak_ht_get(h2, xr>>pre<<YAK_COUNTER_BITS);
					if (k2!=kh_end(h2)){
						cnt_asm = kh_key(h2, k2) & YAK_MAX_COUNT;
						max_outcome = cnt_asm>max_outcome? cnt_asm : max_outcome; 
					}

					xr = xr^(1ULL<<ik);  // mutaiton 3
					pre2 = xr & mask_pre;
					h2 = h_asm->h[pre2].h;
					k2 = yak_ht_get(h2, xr>>pre<<YAK_COUNTER_BITS);
					if (k2!=kh_end(h2)){
						cnt_asm = kh_key(h2, k2) & YAK_MAX_COUNT;
						max_outcome = cnt_asm>max_outcome? cnt_asm : max_outcome; 
					}
				}
				buf[max_outcome]++;
			}
		}
	}

}

yak_ht_t *yam_compare_get_kmers_of_interest(yak_ch_t *h_asm, yak_ch_t *h_reads, int pre,
											int read_low, int read_high, int asm_low, int asm_high){
    yak_ht_t *hret = yak_ht_init();
	double time = Get_T();
	khint_t k, k2;
	uint64_t mask = (1<<YAK_COUNTER_BITS)-1, x, xr, xa;
	int i, j, cnt_read, cnt_asm, absent;
	yak_ht_t *h, *h2;

	int sancheck_counter = 0;
	for (i=0; i<(1<<pre)-1; i++){
		h = h_reads->h[i].h;
		for (k=0; k!=kh_end(h); k++){
			if (kh_exist(h, k)){
				cnt_read = kh_key(h, k) & mask;
				if (cnt_read<read_low || cnt_read>=read_high) continue; 
				xr = kh_key(h, k)>>YAK_COUNTER_BITS<<pre;

				// check if the kmer is present in assembly's hashtable
				h2 = h_asm->h[i].h;
				k2 = yak_ht_get(h2, xr);
				if (k2==kh_end(h2)){
					cnt_asm = 0;
				}else{
					cnt_asm = kh_key(h2, k2) & mask;
				}
				if (cnt_asm<asm_low || cnt_asm>=asm_high) continue;

				// plug
				x = ((kh_key(h, k)>>YAK_COUNTER_BITS)<<pre) | ((uint64_t)i);
				yak_ht_put(hret, x, &absent);
				// assert(absent==1);  // must not have been put, since each kmer could only belong to one place. 
				sancheck_counter++;
			}
		}
	}
	fprintf(stderr, "[M::%s] collected %d kmers of interest, used %.2fs\n", 
					 __func__, sancheck_counter, Get_T()-time);

	return hret;
}

// FUNC
//     Simply count kmers of a seq and return a hashtable.
yak_ht_t *yam_count_seq(char *seq, int seq_l, int32_t k, int *n_distinct){
	yak_ht_t *h = yak_ht_init();
	int i, l, absent, tot=0;
	uint64_t x[2], hash, mask = (1ULL<<k*2) - 1, shift = (k - 1) * 2;
	khint_t key;
	for (i = l = 0, x[0] = x[1] = 0; i < seq_l; ++i) {
		int c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) { // not an "N" base
			x[0] = (x[0] << 2 | c) & mask;                  // forward strand
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
			if (++l >= k) { // we find a k-mer
				uint64_t y = x[0] < x[1]? x[0] : x[1];
				hash = yak_hash64(y, mask);
				yak_ht_put(h, hash, &absent);
				if (absent) tot++;
			}
		} else l = 0, x[0] = x[1] = 0; // if there is an "N", restart
	}
	*n_distinct = tot;
	return h;
}


// FUNC
//     Comparison mode2 core.
//     `h` contains only the desired kmers.
int yam_compare_count_hits_given_seq_and_ht(yak_ht_t *h, char *seq, int seq_l, int32_t k){
	int ret = 0, i, l;
	uint64_t x[2], hash, mask = (1ULL<<k*2) - 1, shift = (k - 1) * 2;
	khint_t key;
	for (i = l = 0, x[0] = x[1] = 0; i < seq_l; ++i) {
		int c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) { // not an "N" base
			x[0] = (x[0] << 2 | c) & mask;                  // forward strand
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
			if (++l >= k) { // we find a k-mer
				uint64_t y = x[0] < x[1]? x[0] : x[1];
				hash = yak_hash64(y, mask);
				key = yak_ht_get(h, hash);
				if (key!=kh_end(h))
					ret++;
			}
		} else l = 0, x[0] = x[1] = 0; // if there is an "N", restart
	}
	return ret;
}

// FUNC
//     Comparison mode2, write BED formatted lines for IGV eyeballing 
//      the kmer locations on the given contigs.
// PAR
//     `h` contains only the desired kmers.
//     `l_gap` for the threshold of merging two blocks
int yam_compare_count_bed_given_seq_and_ht(FILE *fp, yak_ht_t *h, 
										  char *seq, int seq_l, 
										  char *seqname, int seqname_l, 
										  int32_t k, int l_gap){
	int n_blocks = 0, n_bp=0;
	int prv_end=0, start=0, end=0;
	int i, l;
	uint64_t x[2], hash, mask = (1ULL<<k*2) - 1, shift = (k - 1) * 2;
	khint_t key;
	for (i = l = 0, x[0] = x[1] = 0; i < seq_l; ++i) {
		int c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) { // not an "N" base
			x[0] = (x[0] << 2 | c) & mask;                  // forward strand
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
			if (++l >= k) { // we find a k-mer
				uint64_t y = x[0] < x[1]? x[0] : x[1];
				hash = yak_hash64(y, mask);
				key = yak_ht_get(h, hash);
				if (key!=kh_end(h)){  // valid kmer
					if (prv_end+l_gap+k>=i){  // small gap
						prv_end = i;
						end = i;
					}else{  // large gap, flush and update
						fprintf(fp, "%.*s\t%d\t%d\n", seqname_l, seqname, start, end);
						n_bp += end-start;
						start = i-k;
						end = i;
						prv_end = i;
						n_blocks++;
					}
				}
			}
		} else l = 0, x[0] = x[1] = 0; // if there is an "N", restart
	}
	if (end>start){
		fprintf(fp, "%.*s\t%d\t%d\n", seqname_l, seqname, start, end);
		n_bp += end-start;
		n_blocks++;
	}
	fprintf(stderr, "[M::%s] wrote %d blocks, cover %d bp (out of %d).\n", __func__, n_blocks, n_bp, seq_l);
	return n_blocks;
}


void yam_ht_destroy(yak_ht_t *h){
	yak_ht_destroy(h);
}


// FUNC
//     Collect bottom minhash sketch of a sequence.
vec_u64t yam_minhash_sketch(char *seq, int seq_l, int kmersize, int n_hash){
	int h_size = 0, n=0, i;
	yak_ht_t *h = yam_count_seq(seq, seq_l, kmersize, &h_size);
	uint64_t *buf = (uint64_t*)malloc(sizeof(uint64_t)*h_size);
	vec_u64t ret; 
	kv_init(ret);
	kv_resize(uint64_t, ret, n_hash);

	khint_t k;
	uint64_t hash;
	for (k=0; k<kh_end(h); k++){
		if (kh_exist(h, k)){
			hash = kh_key(h, k);
			buf[n++] = hash;
		}
	}
	radix_sort_u64t(buf, buf+n);
	for (i=0; i<n_hash; i++){
		kv_push(uint64_t, ret, buf[i]);
	}
	
	free(buf);
	yak_ht_destroy(h);
	return ret;
}


// FUNC
//     mash distance: ANI estimation using Jaccard estimator. 
//     See Ondov et al (2016) equation 4.
double yam_minhash_mashdist_core(vec_u64t *h1, vec_u64t *h2, int k, int n_hash, int *ret_share){
	double jaccard, dist;
	int tot=0, share=0, i1=0, i2=0, n=h1->n, which=1; 
	assert(h1->n==h2->n);
	while (i1<n && i2<n){
		if (h1->a[i1]==h2->a[i2]){
			share++;
			i1++;
			i2++;
		}else if (h1->a[i1]>h2->a[i2]){
			i2++;
		}else{
			i1++;
		}
		tot++;
		if (tot==n_hash) break;
	}
	if (share==0) dist=1;
	else if (share==n_hash) dist=+0;
	else{
		jaccard = (double)share/n_hash;
		dist = -1.0/(double)k * log(2*jaccard/(1+jaccard));
	}
	if (ret_share) *ret_share = share;
	return dist;
}