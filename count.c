#include <stdio.h>
#include <stdint.h>
#include <zlib.h>
#include <string.h>
#include <assert.h>
#include "kthread.h" // multi-threading models: pipeline and multi-threaded for loop
#include "yak-priv.h"
#include "count.h"
#include "kseq.h" // FASTA/Q parser

KSEQ_INIT(gzFile, gzread)


unsigned char seq_nt4_table[256] = { // translate ACGT to 0123
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

typedef struct {
	int n, m;
	uint64_t n_ins;
	uint64_t *a;
} ch_buf_t;

typedef struct { // global data structure for kt_pipeline()
	const yak_copt_t *copt;
	const yak_opt_t *opt;
	int create_new;
	yak_ch_t *h;   // hashtable to be filled 
	int create_new2;
	yak_ch_t *h2;  // orphans' hashtable to be filled (optional)
	yak_ch_t *href;  // optional: the assembly's hashtable when counting reads
	uint8_t mode;  // bitflag
	kseq_t *ks;
	int n_tot;  // total read count
} pldat_t;

typedef struct { // data structure for each step in kt_pipeline()
	pldat_t *p;
	int n, m, sum_len, nk;
	int *len;
	char **seq;
	ch_buf_t *buf;
	ch_buf_t *buf2;
	vec_u64t v;  // mark orphan
} stepdat_t;

int yak_ch_get(const yak_ch_t *h, uint64_t x);
int yak_ch_insert_list(yak_ch_t *h, int create_new, int n, const uint64_t *a);
yak_ch_t *yak_ch_init(int k, int pre, int n_hash, int n_shift);


void yak_copt_init(yak_copt_t *o)
{
	memset(o, 0, sizeof(yak_copt_t));
	o->bf_shift = 0;
	o->bf_n_hash = 4;
	o->k = 31;
	o->pre = 10;
	o->n_thread = 32;
	o->chunk_size = 10000000;
	// special options for `compare`
	o->is_dump_kmer_loc = 0;
}


static inline void ch_insert_buf(ch_buf_t *buf, int p, uint64_t y) // insert a k-mer $y to a linear buffer
{
	int pre = y & ((1<<p) - 1);
	ch_buf_t *b = &buf[pre];
	if (b->n == b->m) {
		b->m = b->m < 8? 8 : b->m + (b->m>>1);
		REALLOC(b->a, b->m);
	}
	b->a[b->n++] = y;
}

static void count_seq_buf(stepdat_t *s, int i_read){  // ktfor
	int k = s->p->copt->k;
	int p = s->p->copt->pre;
	int len = s->len[i_read];
	char *seq = s->seq[i_read];
	vec_u64t *buf = &s->v;
	buf->n = 0;
	uint64_t hash;
	int i, l, key, hit_kmer=0, n_kmer=0, n_ins;

	uint64_t x[2], mask = (1ULL<<k*2) - 1, shift = (k - 1) * 2;
	for (i = l = 0, x[0] = x[1] = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) { // not an "N" base
			x[0] = (x[0] << 2 | c) & mask;                  // forward strand
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
			if (++l >= k) { // we find a k-mer
				uint64_t y = x[0] < x[1]? x[0] : x[1];
				hash = yak_hash64(y, mask);
				kv_push(uint64_t, (*buf), hash);
				// if hastable from assembly is present, check if kmer is in there
				if (s->p->mode & YC_COUNT_READ && hit_kmer<1000){
					n_kmer++;
					key = yak_ch_get(s->p->href, hash);
					if (key>=0) hit_kmer++;
				}
			}
		} else l = 0, x[0] = x[1] = 0; // if there is an "N", restart
	}

	// if (s->p->mode & YC_COUNT_CONTIG || (hit_kmer/ (float) n_kmer) > s->p->opt->read_hit_ratio){
	if (s->p->mode & YC_COUNT_CONTIG || hit_kmer>=1000){
		for (i=0; i<buf->n; i++){
			ch_insert_buf(s->buf, p, buf->a[i]);	
		}
	}else{
		for (i=0; i<buf->n; i++){
			ch_insert_buf(s->buf2, p, buf->a[i]);	
		}
	}

}


static void count_seq_buf_long(stepdat_t *s, int i_read){  // ktfor
	int k = s->p->copt->k;
	int p = s->p->copt->pre;
	int len = s->len[i_read];
	char *seq = s->seq[i_read];
	vec_u64t *buf = &s->v;
	buf->n = 0;
	uint64_t hash;
	int i, l, key, hit_kmer=0, n_kmer=0;

	uint64_t x[4], mask = (1ULL<<k) - 1, shift = k - 1;
	for (i = l = 0, x[0] = x[1] = x[2] = x[3] = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)seq[i_read]];
		if (c < 4) { // not an "N" base
			x[0] = (x[0] << 1 | (c&1))  & mask;
			x[1] = (x[1] << 1 | (c>>1)) & mask;
			x[2] = x[2] >> 1 | (uint64_t)(1 - (c&1))  << shift;
			x[3] = x[3] >> 1 | (uint64_t)(1 - (c>>1)) << shift;
			if (++l >= k){  // found a kmer
				hash = yak_hash_long(x);
				kv_push(uint64_t, (*buf), hash);
				// if hastable from assembly is present, check if kmer is in there
				if (s->p->mode & YC_COUNT_READ){
					n_kmer++;
					key = yak_ch_get(s->p->href, hash);
					if (key>=0) hit_kmer++;
				}
			}
		} else l = 0, x[0] = x[1] = x[2] = x[3] = 0; // if there is an "N", restart
	}
	if (s->p->mode & YC_COUNT_CONTIG || hit_kmer/ (float) n_kmer > s->p->opt->read_hit_ratio){
		for (i=0; i<buf->n; i++){
			ch_insert_buf(s->buf, p, buf->a[i]);	
		}
	}else{
		for (i=0; i<buf->n; i++){
			ch_insert_buf(s->buf2, p, buf->a[i]);	
		}
	}
}

static void worker_for(void *data, long i, int tid) // callback for kt_for()
{
	stepdat_t *s = (stepdat_t*)data;
	ch_buf_t *b = &s->buf[i];
	yak_ch_t *h = s->p->h;
	yak_ch_t *h2 = s->p->h2;
	b->n_ins += yak_ch_insert_list(h, s->p->create_new, b->n, b->a);
	if (s->p->mode & YC_COUNT_READ){
		ch_buf_t *b2 = &s->buf2[i];
		b2->n_ins += yak_ch_insert_list(h2, s->p->create_new, b2->n, b2->a);
	}
}


static void *worker_pipeline(void *data, int step, void *in) // callback for kt_pipeline()
{
	pldat_t *p = (pldat_t*)data;
	if (step == 0) { // step 1: read a block of sequences
		double tstart = Get_T();
		int ret;
		stepdat_t *s;
		CALLOC(s, 1);
		s->p = p;
		kv_init(s->v);
		s->v.a = (uint64_t*)realloc(s->v.a, sizeof(uint64_t)*1024);
		while ((ret = kseq_read(p->ks)) >= 0) {
			int l = p->ks->seq.l;
			if (l < p->copt->k) continue;
			if (s->n == s->m) {
				s->m = s->m < 16? 16 : s->m + (s->n>>1);
				REALLOC(s->len, s->m);
				REALLOC(s->seq, s->m);
			}
			MALLOC(s->seq[s->n], l);
			memcpy(s->seq[s->n], p->ks->seq.s, l);
			s->len[s->n++] = l;
			s->sum_len += l;
			s->nk += l - p->copt->k + 1;
			if (s->sum_len >= p->copt->chunk_size)
				break;
		}
		// fprintf(stderr, "step0 - %.2f\n", Get_T()-tstart);
		if (s->sum_len == 0) {
			if (s->p->mode & YC_COUNT_READ) {
				kv_destroy(s->v);
			}
			free(s);
		}
		else {
			p->n_tot += s->n;  // update read count
			return s;
		}
	} else if (step == 1) { // step 2: extract k-mers
		stepdat_t *s = (stepdat_t*)in;
		double tstart = Get_T();
		
		int i, n = 1<<p->copt->pre, m;
		m = (int)(s->nk * 1.2 / n) + 1;
		// (allocate kmer linear buffers)
		CALLOC(s->buf, n);
		if (s->p->mode & YC_COUNT_READ) CALLOC(s->buf2, n);
		for (i = 0; i < n; ++i) {
			s->buf[i].m = m;
			MALLOC(s->buf[i].a, m);
			if (s->p->mode & YC_COUNT_READ){
				s->buf2[i].m = m;
				MALLOC(s->buf2[i].a, m);
			}
		}
		for (i=0; i<s->n; i++){
			if (p->copt->k < 32){
				count_seq_buf(s, i);
			}else{
				count_seq_buf_long(s, i);
			}
			free(s->seq[i]);
		}
		// fprintf(stderr, "step1 - %.2f\n", Get_T()-tstart);
		free(s->seq); free(s->len);
		kv_destroy(s->v);
		return s;
	} else {  // step 3: insert to hash tables
		double tstart = Get_T();
		stepdat_t *s = (stepdat_t*)in;
		int i, n = 1<<p->copt->pre;
		uint64_t n_ins = 0, n_ins2 = 0;
		
		kt_for(p->copt->n_thread, worker_for, s, n);

		for (i = 0; i < n; ++i) {
			n_ins += s->buf[i].n_ins;
			free(s->buf[i].a);
			if (s->p->mode & YC_COUNT_READ){
				n_ins2 += s->buf2[i].n_ins;
				free(s->buf2[i].a);
			}
		}
		p->h->tot += n_ins;
		// fprintf(stderr, "step2 - %.2f\n", Get_T()-tstart);
		if (s->p->mode & YC_COUNT_READ) p->h2->tot += n_ins2;		
		free(s->buf);
		if (s->p->mode & YC_COUNT_READ) free(s->buf2);

		fprintf(stderr, "[M::%s] processed %d sequences, create new %d, is contig %d, is reads %d; "
				 "%ld distinct k-mers in the hash table, %ld in orphan h\n", __func__,
				 s->n, 
				 (int)!!(s->p->create_new),
				 (int)!!(s->p->mode & YC_COUNT_CONTIG), (int)!!(s->p->mode & YC_COUNT_READ),
				 (long)p->h->tot,
				p->h2? p->h2->tot : (long)-1);

		free(s);
	}
		
	return 0;
}



yak_ch_v *yak_count(const char *fn, const yak_copt_t *copt, yak_opt_t *opt,
					yak_ch_t *h_asm, yak_ch_t *h_reads, yak_ch_t *h_reads_orphan, uint8_t mode)
{
	// (support biparting read given markers; also report total number of reads)
	pldat_t pl;
	memset(&pl, 0, sizeof(pl));
	gzFile fp;
	if ((fp = gzopen(fn, "r")) == 0) {
		fprintf(stderr, "[E::%s] unable to open input file: %s\n", __func__, fn);
		exit(1);
	}

	pl.ks = kseq_init(fp);
	pl.copt = copt;
	pl.opt = opt;
	pl.n_tot = 0;
	pl.mode = 0;

	if (!h_asm){
		if (mode & YC_COUNT_CONTIG_ALLOW_BF) pl.h = yak_ch_init(copt->k, copt->pre, copt->bf_n_hash, copt->bf_shift);
		else pl.h = yak_ch_init(copt->k, copt->pre, copt->bf_n_hash, 0);  // override, disable bloom filter
		pl.create_new = 1;
		pl.mode |= YC_COUNT_CONTIG;
	}else{
		pl.mode = mode;
		if (mode & YC_COUNT_CONTIG){
			pl.h = h_asm;
			pl.create_new = 0;
		}else if (mode & YC_COUNT_READ){
			assert(h_asm);
			assert (h_asm->k == copt->k && h_asm->pre == copt->pre);
			pl.href = h_asm;
			if (!h_reads){
				pl.create_new = 1;
				pl.h = yak_ch_init(copt->k, copt->pre, copt->bf_n_hash, copt->bf_shift);
				pl.h2 = yak_ch_init(copt->k, copt->pre, copt->bf_n_hash, copt->bf_shift);
			}else{
				pl.create_new = 0;
				assert(h_reads);
				assert(h_reads_orphan);
				pl.h = h_reads;
				pl.h2 = h_reads_orphan;
				assert (h_reads->k == copt->k && h_reads->pre == copt->pre);
				assert (h_reads_orphan->k == copt->k && h_reads_orphan->pre == copt->pre);
			}
		}else{
			fprintf(stderr, "[E::%s] invalid bit flag\n", __func__);
		}
	}

	kt_pipeline(3, worker_pipeline, &pl, 3);
	kseq_destroy(pl.ks);
	gzclose(fp);

	yak_ch_v *ret = (yak_ch_v*)malloc(sizeof(yak_ch_v));
	ret->n = 2;
	ret->hh = (yak_ch_t**)malloc(ret->n * sizeof(yak_ch_t*));
	if (pl.mode & YC_COUNT_CONTIG){
		ret->hh[0] = pl.h;
		ret->n = 1;
	}else{
		ret->hh[0] = pl.h;
		ret->hh[1] = pl.h2;
	}
	return ret;
}