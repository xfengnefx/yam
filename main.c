#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ketopt.h"
#include "yak-priv.h"
#include "count.h"
#include "compare.h"
#include "htab.h"

static inline int64_t mm_parse_num(const char *str)
{
	double x;
	char *p;
	x = strtod(str, &p);
	if (*p == 'G' || *p == 'g') x *= 1e9;
	else if (*p == 'M' || *p == 'm') x *= 1e6;
	else if (*p == 'K' || *p == 'k') x *= 1e3;
	return (int64_t)(x + .499);
}

void yak_opt_init(yak_opt_t *o){
	o->read_hit_ratio = 0.1;
}

// FUNC
//     Assembly completeness estimation by counting "orphan" kmer ratios.
//     An oprphan read is a read that have too many kmers not presented
//      in the assembly. All kmers belong to an orphaned read are
//      counted as orphaned kmers, regardless of individual kmer frequencies.
int main_completeness(int argc, char *argv[]){
    yak_copt_t copt;
    yak_opt_t opt;
    yak_ch_t *h_asm=0, *h_reads=0, *h_reads_orphan=0;
    yak_ch_v *hret;
    yak_copt_init(&copt);  // kmer counting options
    yak_opt_init(&opt);  // general options
    ketopt_t o = KETOPT_INIT;
    
    int c;
    while ((c = ketopt(&o, argc, argv, 1, "k:p:K:t:b:H:r:", 0)) >= 0) {
		if (c == 'k') copt.k = atoi(o.arg);
		else if (c == 'p') copt.pre = atoi(o.arg);
		else if (c == 'K') copt.chunk_size = mm_parse_num(o.arg);
		else if (c == 't') copt.n_thread = atoi(o.arg);
		else if (c == 'b') copt.bf_shift = atoi(o.arg);
		else if (c == 'H') copt.bf_n_hash = mm_parse_num(o.arg);
        else if (c == 'r') opt.read_hit_ratio = atof(o.arg);
	}
    if (argc - o.ind < 1) {
		fprintf(stderr, "Usage: yak completeness [options] <contig.fa> <reads.fa>\n");
        fprintf(stderr, "Options - general:\n");
        fprintf(stderr, "  -r FLOAT   orphan criteria\n");
		fprintf(stderr, "Options - kmer counting:\n");
		fprintf(stderr, "  -k INT     k-mer size [%d]\n", copt.k);
		fprintf(stderr, "  -p INT     prefix length [%d]\n", copt.pre);
		fprintf(stderr, "  -b INT     set Bloom filter size to 2**INT bits; 0 to disable [%d]\n", copt.bf_shift);
		fprintf(stderr, "  -H INT     use INT hash functions for Bloom filter [%d]\n", copt.bf_n_hash);
		fprintf(stderr, "  -t INT     number of worker threads [%d]\n", copt.n_thread);
		fprintf(stderr, "  -K INT     chunk size [100m]\n");
		fprintf(stderr, "Note: -b37 is recommended for human reads\n");
		return 1;
	}
    
    // collect contig kmers
    hret = yak_count(argv[o.ind], &copt, &opt, 0, 0, 0, (uint8_t)YC_COUNT_CONTIG);
    h_asm = hret->hh[0];
    free(hret->hh);
    free(hret);

    // collect read kmers
    hret = yak_count(argv[o.ind+1], &copt, &opt, h_asm, 0, 0, (uint8_t)YC_COUNT_READ);
    h_reads = hret->hh[0];
    h_reads_orphan = hret->hh[1];
    free(hret->hh);
    free(hret);

    
    yak_ch_destroy_bf(h_asm);
    yak_ch_destroy(h_asm);
    yak_ch_destroy_bf(h_reads);
    yak_ch_destroy(h_reads);
    yak_ch_destroy_bf(h_reads_orphan);
    yak_ch_destroy(h_reads_orphan);
    
    
    return 0;
}


int main_katcompare(int argc, char *argv[]){
    yak_copt_t copt;
    yak_copt_init(&copt);
    yak_opt_t opt;
    yak_opt_init(&opt);
    
    yak_ch_v *hret;
    yak_ch_t *h_asm, *h_reads;

    ketopt_t o = KETOPT_INIT;
    int c;
    while ((c = ketopt(&o, argc, argv, 1, "k:p:K:t:b:H:o:dz", 0)) >=0 ){
        if (c == 'k') copt.k = atoi(o.arg);
        else if (c == 'K') copt.chunk_size = mm_parse_num(o.arg);
		else if (c == 't') copt.n_thread = atoi(o.arg);
		else if (c == 'b') copt.bf_shift = atoi(o.arg);
		else if (c == 'H') copt.bf_n_hash = mm_parse_num(o.arg);
        else if (c == 'o') copt.output_name = o.arg;
        else if (c == 'd') copt.is_dump_kmer_loc = 1;
        else if (c == 'z') copt.is_dump_zero_fasta = 1;
    }
    if (argc - o.ind < 1) {
		fprintf(stderr, "Usage: yam compare [options] <asm.ht> <reads.ht> [1>zerodump.fa]"
                         "[asm.fa, read_low, read_high, asm_low, asm_high]\n");
        fprintf(stderr, " There are two modes: when only two files (or two streams) are supplied,\n");
        fprintf(stderr, "  yam does the kat-like read-to-contig comparison.\n");
        fprintf(stderr, "  If 1 file and 4 numbers are supplied, dump names of contigs that contain\n"
                          "kmers of certain read & assembly occurrences.\n\n");
        fprintf(stderr, "Options - general:\n");
        fprintf(stderr, "  -o STR     2D count matrix's prefix\n");
        fprintf(stderr, "  -d         (validation mode) also write a bed file showing all kmer hits (intervals merged)\n");
        fprintf(stderr, "  -z         write non-contig kmers to a fasta (to STDOUT)\n");
		fprintf(stderr, "Options - kmer counting:\n");
		fprintf(stderr, "  -k INT     k-mer size [%d]\n", copt.k);
		fprintf(stderr, "  -p INT     prefix length [%d]\n", copt.pre);
		fprintf(stderr, "  -b INT     set Bloom filter size to 2**INT bits; 0 to disable [%d]\n", copt.bf_shift);
		fprintf(stderr, "  -H INT     use INT hash functions for Bloom filter [%d]\n", copt.bf_n_hash);
		fprintf(stderr, "  -t INT     number of worker threads [%d]\n", copt.n_thread);
		fprintf(stderr, "  -K INT     chunk size [100m]\n");
		fprintf(stderr, "Note: -b37 is recommended for human reads\n");
		return 1;
	}
    // sanchecks
    if ( (argc-o.ind)>2 && (argc-o.ind)!=7){
        fprintf(stderr, "[E::%s] validation mode, but incorrect number of arguments. "
                         "You need 2 files and 4 numbers.\n", __func__);
        return 1;
    }

    if ((argc-o.ind)!=7){
        if (copt.output_name==0){
            fprintf(stderr, "[E::%s] must provide output prefix with -o.\n", __func__);
            return 1;
        }
        return main_katcompare_count_mode(argv[o.ind], argv[o.ind+1], &copt);        
    }else{
        return main_katcompare_validation_mode(argv[o.ind], argv[o.ind+1], argv[o.ind+2], &copt,
                                        atoi(argv[o.ind+3]), atoi(argv[o.ind+4]), atoi(argv[o.ind+5]), atoi(argv[o.ind+6]));
    }
}

int main_ht(int argc, char *argv[])
{
	yak_copt_t copt;
    yak_opt_t opt;
    yak_ch_t *h_asm=0, *h_reads=0, *h_reads_orphan=0;
    yak_ch_v *hret;
    yak_copt_init(&copt);  // kmer counting options
    yak_opt_init(&opt);  // general options
    ketopt_t o = KETOPT_INIT;

    yak_ch_t *h;
    char *fn_out;
    
    int c, ignore_assembly=0, ignore_reads=0;
    while ((c = ketopt(&o, argc, argv, 1, "k:p:K:t:b:H:12o:", 0)) >= 0) {
		if (c == 'k') copt.k = atoi(o.arg);
		else if (c == 'p') copt.pre = atoi(o.arg);
		else if (c == 'K') copt.chunk_size = mm_parse_num(o.arg);
		else if (c == 't') copt.n_thread = atoi(o.arg);
		else if (c == 'b') copt.bf_shift = atoi(o.arg);
		else if (c == 'H') copt.bf_n_hash = mm_parse_num(o.arg);
        else if (c == '1') ignore_assembly = 1;
        else if (c == '2') ignore_reads = 1;
        else if (c == 'o') copt.output_name = o.arg;
	}
    if (argc - o.ind < 1) {
		fprintf(stderr, "Usage: yak count -o prefix [options] <assembly.fa> <reads.fa>\n");
		fprintf(stderr, "Options - general:\n");
        fprintf(stderr, "  -o STR     output prefix\n");
        fprintf(stderr, "Options - to count only one file:\n");
        fprintf(stderr, "  -1         skip counting assembly\n");
        fprintf(stderr, "  -2         skip counting reads\n");
        fprintf(stderr, "Options - kmer counting:\n");
		fprintf(stderr, "  -k INT     k-mer size [%d]\n", copt.k);
		fprintf(stderr, "  -p INT     prefix length [%d]\n", copt.pre);
		fprintf(stderr, "  -b INT     set Bloom filter size to 2**INT bits; 0 to disable [%d]\n", copt.bf_shift);
		fprintf(stderr, "  -H INT     use INT hash functions for Bloom filter [%d]\n", copt.bf_n_hash);
		fprintf(stderr, "  -t INT     number of worker threads [%d]\n", copt.n_thread);
		fprintf(stderr, "  -K INT     chunk size [100m]\n");
		fprintf(stderr, "Note: -b37 is recommended for human reads\n");
		return 1;
	}
    if (argc-o.ind<2){
        fprintf(stderr, "[E::%s] must provide two files.\n", __func__);
        return 1;
    }
    if (!copt.output_name){
        fprintf(stderr, "[E::%s] must give output prefix.\n", __func__);
        return 1;
    }

    fn_out = (char*)malloc(strlen(copt.output_name)+25);

    // count the assembly
    if (!ignore_assembly){
        fprintf(stderr, "[M::%s] counting assembly '%s'..\n", __func__, argv[o.ind]);
        sprintf(fn_out, "%s.asm.hashtables", copt.output_name);
        h = yak_count(argv[o.ind], &copt, &opt, 0, 0, 0, YC_COUNT_CONTIG)->hh[0];
        yak_ch_dump(h, fn_out);
        yak_ch_destroy(h);
    }

    // count the reads
    if (!ignore_reads){
        fprintf(stderr, "[M::%s] counting reads '%s'..\n", __func__, argv[o.ind+1]);
        sprintf(fn_out, "%s.reads.hashtables", copt.output_name);
        h = yak_count(argv[o.ind+1], &copt, &opt, 0, 0, 0, YC_COUNT_CONTIG | YC_COUNT_CONTIG_ALLOW_BF)->hh[0];
        yak_ch_dump(h, fn_out);
        if (copt.bf_shift > 0) {
            yak_ch_destroy_bf(h);
        }
        yak_ch_destroy(h);
    }


    free(fn_out);
	return 0;
}

// mash dist-like; for testing
int main_dist(int argc, char *argv[]){
    int k=21, n_hash=1000, c;
    int read_files = 0;
    ketopt_t o = KETOPT_INIT;
    yak_ch_t *h;
    while ((c = ketopt(&o, argc, argv, 1, "k:s:f:", 0)) >= 0) {
		if (c == 'k') k = atoi(o.arg);
		else if (c == 's') n_hash = atoi(o.arg);
        else if (c == 'f') read_files = atoi(o.arg);  // one or both input is/are list of file names.  
	}
    if (argc - o.ind < 2) {
        fprintf(stderr, "Usage: yam dist <-k, -s, -f> in1.fa in2.fa 1>dist.tsv\n");   
        fprintf(stderr, "   -f can take 1, 2, or 3, corresponding to : in2 is list, in1 is list, both are lists.\n");
        fprintf(stderr, "      (^^^ if input is file, will sketch each fa/fq entry. If is file list, sketch per file.)\n");
        return 1;
    }
    return main_dist_via_sketch(argv[o.ind], argv[o.ind+1], read_files, k, n_hash);
}

int main(int argc, char *argv[]){
    int ret = 0;
    if (argc==1){
        fprintf(stderr, "Usage: yam <command> <arguments>\n");
        fprintf(stderr, "Command:\n");
        fprintf(stderr, "  ht              dump hashtable.\n");
        // fprintf(stderr, "  completeness    given contigs and reads, calculate completeness.\n");
        fprintf(stderr, "  compare         like kat comp.\n");
        // fprintf(stderr, "  dist            mash-dist pairwise for testing.\n");
        return 1;
    }
    if (strcmp(argv[1], "completeness")==0) ret = main_completeness(argc-1, argv+1);
    else if (strcmp(argv[1], "compare")==0) ret = main_katcompare(argc-1, argv+1);
    else if (strcmp(argv[1], "ht")==0) ret = main_ht(argc-1, argv+1);
    else if (strcmp(argv[1], "dist")==0) ret = main_dist(argc-1, argv+1);
    else{
        fprintf(stderr, "[E::%s] Unknown command\n", __func__);
        return 1;
    }

    return ret;

}
