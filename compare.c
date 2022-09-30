#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include <assert.h>
#include "khashl.h"
#include "htab.h"
#include "yak-priv.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

// FUNC
//     Produces a `kat comp`-like result.
//     Write one file: prefix.mat (similar to kat's .mx file).
int main_katcompare_count_mode(const char *fn_asm, const char *fn_reads, yak_copt_t *copt){
    yak_ch_t *h_asm, *h_reads;
    
    // load from file
    h_asm = yak_ch_restore(fn_asm);
    h_reads = yak_ch_restore(fn_reads);
    
    // open files to write to
    char *fn = (char*)malloc(strlen(copt->output_name)+50);
    sprintf(fn, "%s.mat", copt->output_name);
    FILE *fp_mat = fopen(fn, "w"); assert(fp_mat);

    // init an empty matrix to accumulate counts
    int **mat = (int**)calloc(YAK_MAX_COUNT+1, sizeof(int*));
    int i, j;
    for (i=0; i<YAK_MAX_COUNT+1; i++){
        mat[i] = (int*)calloc(YAK_MAX_COUNT+1, sizeof(int));
    }

    // fill the matrix
    yam_compare_count_2D(h_asm, h_reads, mat, copt->pre);

    // write matrix
    for (i=0; i<YAK_MAX_COUNT+1; i++){  // iterate over kmer counts in reads
        for (j=0; j<YAK_MAX_COUNT+1; j++){  // iterate over kmer counts in the assembly
            fprintf(fp_mat, "%d ", mat[i][j]);  // separator is space, in order to be compatible with kat's mx
        }
        fprintf(fp_mat, "\n");
    }

    // dump non-contig kmers to fasta?
    if (copt->is_dump_zero_fasta){
        yam_compare_zero_dump(h_asm, h_reads, copt->k, copt->pre);
    }

    free(fn);
    fclose(fp_mat);
    for (i=0; i<YAK_MAX_COUNT+1; i++){
        free(mat[i]);
    }
    free(mat);
    yak_ch_destroy(h_asm);
    yak_ch_destroy(h_reads);
    return 0;
}


// FUNC
//     Print contig names of those who have kmers of interest.
// TODO
//     threading
int main_katcompare_validation_mode(const char *fn_asm_ht, const char *fn_reads, const char *fn_asm_fa, 
                                    yak_copt_t *copt, int read_low, int read_high, int asm_low, int asm_high){
    yak_ch_t *h_asm, *h_reads;
    FILE *fp;
    gzFile fg;
    yak_ht_t *h_ovlp;
    int ret;
    kseq_t *ks;

    // open or try to open files
    if ((fp = fopen(fn_asm_ht, "rb")) == 0) {
        fprintf(stderr, "[E::%s] can't open '%s'.\n", __func__, fn_asm_ht);
        exit(1);
    }
    fclose(fp);
    if ((fp = fopen(fn_reads, "rb")) == 0) {
        fprintf(stderr, "[E::%s] can't open '%s'.\n", __func__, fn_reads);
        exit(1);
    }
    fclose(fp);
    if ((fg=gzopen(fn_asm_fa, "r")) == 0 ){
        fprintf(stderr, "[E::%s] can't open '%s'.\n", __func__, fn_asm_fa);
    }
    ks = kseq_init(fg);

    // open the output file
    char *fn = (char*)malloc(strlen(copt->output_name)+50);
    sprintf(fn, "%s.val.tsv", copt->output_name);
    FILE *fp_out = fopen(fn, "w");
    FILE *fp_bed = 0;
    if (copt->is_dump_kmer_loc){
        sprintf(fn, "%s.val.bed", copt->output_name);
        fp_bed = fopen(fn, "w");
    }

    // load hashtables
    h_asm = yak_ch_restore(fn_asm_ht);
    h_reads = yak_ch_restore(fn_reads);

    // plug kmers-of-interest into a new hashtable
    h_ovlp = yam_compare_get_kmers_of_interest(h_asm, h_reads, copt->pre, 
                                                read_low, read_high, asm_low, asm_high);

    // go through the assembly fasta and collect kmer counts
    int total = 0, cnt, total_contigs=0, total_no_hit = 0;
    fprintf(fp_out, "#read_low,read_high,asm_low,asm_high\t%d,%d,%d,%d\n", read_low, read_high, asm_low, asm_high);
    while ((ret=kseq_read(ks))>=0){
        cnt = yam_compare_count_hits_given_seq_and_ht(h_ovlp, ks->seq.s, ks->seq.l, copt->k);
        fprintf(fp_out, "%.*s\t%d\t%d\n", (int)ks->name.l, ks->name.s, (int)ks->seq.l, cnt);

        if (fp_bed)
            yam_compare_count_bed_given_seq_and_ht(fp_bed, h_ovlp, 
                                                    ks->seq.s, ks->seq.l, ks->name.s, (int)ks->name.l,
                                                    copt->k, copt->k*5);

        total+=cnt;
        total_contigs++;
        if (cnt==0) total_no_hit++;
    }
    gzclose(fg);
    fprintf(stderr, "[M::%s] %d contigs, %d had no hit (%.2f%%); %d total hits.\n",
                      __func__, total_contigs, total_no_hit, (float)total_no_hit/total_contigs*100.0,
                      total);
    
    fclose(fp_out);
    if (fp_bed) fclose(fp_bed);
    yak_ch_destroy(h_asm);
    yak_ch_destroy(h_reads);
    yam_ht_destroy(h_ovlp);
    return 0;
}


#define get_seqname_l(ll,i) (ll).a[i+1]-(ll).a[i]
#define get_seqname_s(ll,i) (ll).a[i]
void sketch_one_file(kseq_t *ks, vec_u64v *sketches, vec_char *seqnames, vec_u32t *seqnamesl,
                    int k, int n_hash, int is_sketch_whole_file, char* name, int name_l){
    int stat, i;
    uint32_t tmpl;
    if (!is_sketch_whole_file){
        while ((stat=kseq_read(ks))>=0){
            kv_push(vec_u64t, *sketches, yam_minhash_sketch(ks->seq.s, ks->seq.l, k, n_hash));
            for (i=0; i<ks->name.l; i++){
                kv_push(char, *seqnames, ks->name.s[i]);
            }
            tmpl = seqnamesl->a[seqnamesl->n-1]+ks->name.l;
            kv_push(uint32_t, *seqnamesl, tmpl);
        }
    }else{
        fprintf(stderr, "[M::%s] sketching %.*s\n", __func__, name_l, name);
        vec_u64t dd;
        kv_init(dd);
        kv_resize(uint64_t, dd, 1024);
        while ((stat=kseq_read(ks))>=0){
            vec_u64t d = yam_minhash_sketch(ks->seq.s, ks->seq.l, k, n_hash);
            for (i=0; i<d.n; i++){
                kv_push(uint64_t, dd, d.a[i]);
            }
            kv_destroy(d);
        }
        radix_sort_u64t(dd.a, dd.a+dd.n);
        dd.n = n_hash;  // TODO: shrink the buffer
        kv_push(vec_u64t, *sketches, dd);
        for (i=0; i<name_l; i++){
            kv_push(char, *seqnames, name[i]);
        }
        tmpl = seqnamesl->a[seqnamesl->n-1]+name_l;
        kv_push(uint32_t, *seqnamesl, tmpl);
    }
}
int main_dist_via_sketch(const char *fn, const char *fn2, int fl_mask, int k, int n_hash){
    gzFile fg = gzopen(fn, "r");
    if (fg==0){
        fprintf(stderr, "[E::%s] can't open file '%s'\n", __func__, fn);
        return 1;
    }
    gzFile fg2 = gzopen(fn2, "r");
    if (fg2==0){
        fprintf(stderr, "[E::%s] can't open file '%s'\n", __func__, fn2);
        return 1;
    }

    kseq_t *ks;
    int stat, i, j;
    uint32_t tmpl;
    vec_u64v sketches, sketches2;
    vec_char seqnames, seqnames2, filename;
    vec_u32t seqnamesl, seqnamesl2;
    kv_init(filename); kv_resize(char, filename, 32);
    kv_init(sketches); kv_resize(vec_u64t, sketches, 32);
    kv_init(seqnames); kv_resize(char, seqnames, 64);
    kv_init(seqnamesl);  kv_resize(uint32_t, seqnamesl, 64);
    kv_init(sketches2); kv_resize(vec_u64t, sketches2, 32);
    kv_init(seqnames2); kv_resize(char, seqnames2, 64);
    kv_init(seqnamesl2);  kv_resize(uint32_t, seqnamesl2, 64);

    kv_push(uint32_t, seqnamesl, 0);
    kv_push(uint32_t, seqnamesl2, 0);

    // sketch all - file1
    if (fl_mask & 2){  // first file is a list of files, sketch each file
        char ch;
        while (gzeof(fg)==0){
            ch = gzgetc(fg);
            if (ch=='\n'){
                kv_push(char, filename, '\0');
                gzFile fg_sub = gzopen(filename.a, "r");
                if (!fg_sub){
                    fprintf(stderr, "[E::%s] can't open file: %s\n", __func__, filename.a);
                    filename.n = 0;
                    continue;
                }
                ks = kseq_init(fg_sub);
                sketch_one_file(ks, &sketches, &seqnames, &seqnamesl, k, n_hash, 1, filename.a, filename.n);
                gzclose(fg_sub);
                kseq_destroy(ks);
                filename.n = 0;
            }else{
                kv_push(char, filename, ch);
            }
        }
    }else{  // first file is a fasta/q, sketch each entry
        ks = kseq_init(fg);
        sketch_one_file(ks, &sketches, &seqnames, &seqnamesl, k, n_hash, 0, "", 0);
        kseq_destroy(ks);
    }
    gzclose(fg);

    // sketch all - file2
    if (fl_mask & 1){  // 2nd file is a list of files, sketch each file
        char ch;
        while (gzeof(fg2)==0){
            ch = gzgetc(fg2);
            if (ch=='\n'){
                kv_push(char, filename, '\0');
                gzFile fg_sub = gzopen(filename.a, "r");
                if (!fg_sub){
                    fprintf(stderr, "[E::%s] can't open file: %s\n", __func__, filename.a);
                    filename.n = 0;
                    continue;
                }
                ks = kseq_init(fg_sub);
                sketch_one_file(ks, &sketches2, &seqnames2, &seqnamesl2, k, n_hash, 1, filename.a, filename.n);
                gzclose(fg_sub);
                kseq_destroy(ks);
                filename.n = 0;
            }else{
                kv_push(char, filename, ch);
            }
        }
    }else{  // first file is a fasta/q, sketch each entry
        ks = kseq_init(fg2);
        sketch_one_file(ks, &sketches2, &seqnames2, &seqnamesl2, k, n_hash, 0, "", 0);
        kseq_destroy(ks);
    }
    gzclose(fg2);


    // pairwise comparison
    int share;
    for (i=0; i<sketches.n; i++){
        for (j=0; j<sketches2.n; j++){
            double dist = yam_minhash_mashdist_core(&sketches.a[i], &sketches2.a[j], k, n_hash, &share);
            fprintf(stdout, "%.*s\t%.*s\t%.6f\t%d/%d\n", 
                            (int)get_seqname_l(seqnamesl,i), &seqnames.a[get_seqname_s(seqnamesl, i)],
                            (int)get_seqname_l(seqnamesl2,j), &seqnames2.a[get_seqname_s(seqnamesl2, j)],
                            dist, share, n_hash
                    );
        }
    }
    fflush(stdout);
    
    for (i=0; i<sketches.n; i++){
        kv_destroy(sketches.a[i]);
    }
    for (i=0; i<sketches2.n; i++){
        kv_destroy(sketches2.a[i]);
    }
    kv_destroy(filename);
    kv_destroy(sketches);
    kv_destroy(seqnames);
    kv_destroy(seqnamesl);
    kv_destroy(sketches2);
    kv_destroy(seqnames2);
    kv_destroy(seqnamesl2);
    return 0;
}