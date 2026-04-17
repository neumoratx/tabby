/*  tabix.c -- Generic indexer for TAB-delimited genome position files.

    Copyright (C) 2009-2011 Broad Institute.
    Copyright (C) 2010-2012, 2014-2020, 2024 Genome Research Ltd.

    Author: Heng Li <lh3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <config.h>

#include <assert.h>
#include <stdio.h>
#include <stdint.h>
#include <limits.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <strings.h>
#include <getopt.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include "htslib/tbx.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "htslib/kseq.h"
#include "htslib/bgzf.h"
#include "htslib/hts.h"
#include "htslib/regidx.h"
#include "htslib/hts_defs.h"
#include "htslib/hts_log.h"
#include "htslib/thread_pool.h"
//#include "htslib/tbx.h"
#include "htslib/tbx_cw.h"

/* seqonly preset: index by chromosome name only, no coordinate columns.
 * Defined in tbx.c; TBX_SEQONLY must match the value there. */
#ifndef TBX_SEQONLY
#define TBX_SEQONLY (1 << 20)
#endif
extern const tbx_conf_t tbx_conf_seqonly;

//for easy coding
#define RELEASE_TPOOL(X) { hts_tpool *ptr = (hts_tpool*)(X); if (ptr) { hts_tpool_destroy(ptr); } }
#define bam_index_build3(fn, min_shift, nthreads) (sam_index_build3((fn), NULL, (min_shift), (nthreads)))

typedef struct
{
    char *regions_fname, *targets_fname;
    int print_header, header_only, cache_megs, download_index, separate_regs, threads;
}
args_t;

static void HTS_FORMAT(HTS_PRINTF_FMT, 1, 2) HTS_NORETURN
error(const char *format, ...)
{
    va_list ap;
    fflush(stdout);
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    fflush(stderr);
    exit(EXIT_FAILURE);
}

static void HTS_FORMAT(HTS_PRINTF_FMT, 1, 2) HTS_NORETURN
error_errno(const char *format, ...)
{
    va_list ap;
    int eno = errno;
    fflush(stdout);
    if (format) {
        va_start(ap, format);
        vfprintf(stderr, format, ap);
        va_end(ap);
    }
    if (eno) {
        fprintf(stderr, "%s%s\n", format ? ": " : "", strerror(eno));
    } else {
        fprintf(stderr, "\n");
    }
    fflush(stderr);
    exit(EXIT_FAILURE);
}


#define IS_GFF  (1<<0)
#define IS_BED  (1<<1)
#define IS_SAM  (1<<2)
#define IS_VCF  (1<<3)
#define IS_BCF  (1<<4)
#define IS_BAM  (1<<5)
#define IS_CRAM (1<<6)
#define IS_GAF  (1<<7)
#define IS_TXT  (IS_GFF|IS_BED|IS_SAM|IS_VCF)

int file_type(const char *fname)
{
    int l = strlen(fname);
    if (l>=7 && strcasecmp(fname+l-7, ".gff.gz") == 0) return IS_GFF;
    else if (l>=7 && strcasecmp(fname+l-7, ".bed.gz") == 0) return IS_BED;
    else if (l>=7 && strcasecmp(fname+l-7, ".sam.gz") == 0) return IS_SAM;
    else if (l>=7 && strcasecmp(fname+l-7, ".vcf.gz") == 0) return IS_VCF;
    else if (l>=4 && strcasecmp(fname+l-4, ".bcf") == 0) return IS_BCF;
    else if (l>=4 && strcasecmp(fname+l-4, ".bam") == 0) return IS_BAM;
    else if (l>=5 && strcasecmp(fname+l-5, ".cram") == 0) return IS_CRAM;
    else if (l>=7 && strcasecmp(fname+l-7, ".gaf.gz") == 0) return IS_GAF;

    htsFile *fp = hts_open(fname,"r");
    if (!fp) {
        if (errno == ENOEXEC) {
            // hts_open() uses this to report that it didn't understand the
            // file format.
            error("Couldn't understand format of \"%s\"\n", fname);
        } else {
            error_errno("Couldn't open \"%s\"", fname);
        }
    }
    enum htsExactFormat format = hts_get_format(fp)->format;
    hts_close(fp);
    if ( format == bcf ) return IS_BCF;
    if ( format == bam ) return IS_BAM;
    if ( format == cram ) return IS_CRAM;
    if ( format == vcf ) return IS_VCF;

    return 0;
}

static char **parse_regions(char *regions_fname, char **argv, int argc, int *nregs)
{
    kstring_t str = {0,0,0};
    int iseq = 0, ireg = 0;
    char **regs = NULL;
    *nregs = argc;

    if ( regions_fname )
    {
        // improve me: this is a too heavy machinery for parsing regions...

        regidx_t *idx = regidx_init(regions_fname, NULL, NULL, 0, NULL);
        if ( !idx ) {
            error_errno("Could not build region list for \"%s\"", regions_fname);
        }
        regitr_t *itr = regitr_init(idx);
        if ( !itr ) {
            error_errno("Could not initialize an iterator over \"%s\"",
                        regions_fname);
        }

        (*nregs) += regidx_nregs(idx);
        regs = (char**) malloc(sizeof(char*)*(*nregs));
        if (!regs) error_errno(NULL);

        int nseq;
        char **seqs = regidx_seq_names(idx, &nseq);
        for (iseq=0; iseq<nseq; iseq++)
        {
            if (regidx_overlap(idx, seqs[iseq], 0, HTS_POS_MAX, itr) < 0)
                error_errno("Failed to build overlapping regions list");

            while ( regitr_overlap(itr) )
            {
                str.l = 0;
                if (ksprintf(&str, "%s:%"PRIhts_pos"-%"PRIhts_pos, seqs[iseq], itr->beg+1, itr->end+1) < 0) {
                    error_errno(NULL);
                }
                regs[ireg] = strdup(str.s);
                if (!regs[ireg]) error_errno(NULL);
                ireg++;
            }
        }
        regidx_destroy(idx);
        regitr_destroy(itr);
    }
    free(str.s);

    if ( !ireg )
    {
        if ( argc )
        {
            regs = (char**) malloc(sizeof(char*)*argc);
            if (!regs) error_errno(NULL);
        }
        else
        {
            regs = (char**) malloc(sizeof(char*));
            if (!regs) error_errno(NULL);
            regs[0] = strdup(".");
            if (!regs[0]) error_errno(NULL);
            *nregs = 1;
        }
    }

    for (iseq=0; iseq<argc; iseq++, ireg++) {
        regs[ireg] = strdup(argv[iseq]);
        if (!regs[ireg]) error_errno(NULL);
    }
    return regs;
}
static int query_regions(args_t *args, tbx_conf_t *conf, char *fname, char **regs, int nregs)
{
    int i;
    htsThreadPool tpool = {NULL, 0};
    htsFile *fp = hts_open(fname,"r");
    if ( !fp ) error_errno("Could not open \"%s\"", fname);
    enum htsExactFormat format = hts_get_format(fp)->format;
    if (args->cache_megs)
        hts_set_cache_size(fp, args->cache_megs * 1048576);

    //set threads if needed, errors are logged and ignored
    if (args->threads >= 1) {
        if (!(tpool.pool = hts_tpool_init(args->threads))) {
            hts_log_info("Could not initialize thread pool!");
        }
        if (hts_set_thread_pool(fp, &tpool) < 0) {
            hts_log_info("Could not set thread pool!");
        }
    }

    regidx_t *reg_idx = NULL;
    if ( args->targets_fname )
    {
        reg_idx = regidx_init(args->targets_fname, NULL, NULL, 0, NULL);
        if (!reg_idx) {
            RELEASE_TPOOL(tpool.pool);
            error_errno("Could not build region list for \"%s\"",
                        args->targets_fname);
        }
    }

    if ( format == bcf )
    {
        htsFile *out = hts_open("-","w");
        if ( !out ) {
            RELEASE_TPOOL(tpool.pool);
            error_errno("Could not open stdout");
        }
        if (hts_set_thread_pool(out, &tpool) < 0) {
            hts_log_info("Could not set thread pool to output file!");
        }
        hts_idx_t *idx = bcf_index_load3(fname, NULL, args->download_index ? HTS_IDX_SAVE_REMOTE : 0);
        if ( !idx ) {
            RELEASE_TPOOL(tpool.pool);
            error_errno("Could not load .csi index of \"%s\"", fname);
        }

        bcf_hdr_t *hdr = bcf_hdr_read(fp);
        if ( !hdr ) {
            RELEASE_TPOOL(tpool.pool);
            error_errno("Could not read the header from \"%s\"", fname);
        }

        if ( args->print_header ) {
            if ( bcf_hdr_write(out,hdr)!=0 ) {
                RELEASE_TPOOL(tpool.pool);
                error_errno("Failed to write to stdout");
            }
        }
        if ( !args->header_only )
        {
            assert(regs != NULL);
            bcf1_t *rec = bcf_init();
            if (!rec) {
                RELEASE_TPOOL(tpool.pool);
                error_errno(NULL);
            }
            for (i=0; i<nregs; i++)
            {
                int ret, found = 0;
                hts_itr_t *itr = bcf_itr_querys(idx,hdr,regs[i]);
                if (!itr) continue;
                while ((ret = bcf_itr_next(fp, itr, rec)) >=0 )
                {
                    if ( reg_idx )
                    {
                        const char *chr = bcf_seqname(hdr,rec);
                        if (!chr) {
                            RELEASE_TPOOL(tpool.pool);
                            error("Bad BCF record in \"%s\" : "
                                  "Invalid CONTIG id %d\n",
                                  fname, rec->rid);
                        }
                        if ( !regidx_overlap(reg_idx,chr,rec->pos,rec->pos+rec->rlen-1, NULL) ) continue;
                    }
                    if (!found) {
                        if (args->separate_regs) printf("%c%s\n", conf->meta_char, regs[i]);
                        found = 1;
                    }
                    if ( bcf_write(out,hdr,rec)!=0 ) {
                        RELEASE_TPOOL(tpool.pool);
                        error_errno("Failed to write to stdout");
                    }
                }

                if (ret < -1) {
                    RELEASE_TPOOL(tpool.pool);
                    error_errno("Reading \"%s\" failed", fname);
                }
                bcf_itr_destroy(itr);
            }
            bcf_destroy(rec);
        }
        if ( hts_close(out) ) {
            RELEASE_TPOOL(tpool.pool);
            error_errno("hts_close returned non-zero status for stdout");
        }

        bcf_hdr_destroy(hdr);
        hts_idx_destroy(idx);
    }
    else if ( format==vcf || format==sam || format==bed || format==text_format || format==unknown_format )
    {
        tbx_t *tbx = tbx_index_load3(fname, NULL, args->download_index ? HTS_IDX_SAVE_REMOTE : 0);
        if ( !tbx ) {
            RELEASE_TPOOL(tpool.pool);
            error_errno("Could not load .tbi/.csi index of %s", fname);
        }
        kstring_t str = {0,0,0};
        if ( args->print_header )
        {
            int ret;
            while ((ret = hts_getline(fp, KS_SEP_LINE, &str)) >= 0)
            {
                if ( !str.l || str.s[0]!=tbx->conf.meta_char ) break;
                if (puts(str.s) < 0) {
                    RELEASE_TPOOL(tpool.pool);
                    error_errno("Error writing to stdout");
                }
            }
            if (ret < -1) {
                RELEASE_TPOOL(tpool.pool);
                error_errno("Reading \"%s\" failed", fname);
            }
        }
        if ( !args->header_only )
        {
            int nseq;
            const char **seq = NULL;
            if ( reg_idx ) {
                seq = tbx_seqnames(tbx, &nseq);
                if (!seq) {
                    RELEASE_TPOOL(tpool.pool);
                    error_errno("Failed to get sequence names list");
                }
            }
            /*
             * TBX_SEQONLY: records are indexed at [0,1) on each chromosome.
             * A bare chromosome name (no ':') must be expanded to "CHROM:1-1"
             * (1-based, so the query window [0,1) in 0-based is covered) so
             * that tbx_itr_querys performs a whole-chromosome fetch.
             */
            int is_seqonly = (tbx->conf.preset & TBX_SEQONLY) != 0;
            for (i=0; i<nregs; i++)
            {
                int ret, found = 0;
                char *region = regs[i];
                kstring_t expanded = {0, 0, NULL};
                if (is_seqonly && strchr(region, ':') == NULL) {
                    if (ksprintf(&expanded, "%s:1-1", region) < 0) {
                        RELEASE_TPOOL(tpool.pool);
                        error_errno(NULL);
                    }
                    region = expanded.s;
                }
                hts_itr_t *itr = tbx_itr_querys(tbx, region);
                free(expanded.s);
                if ( !itr ) continue;
                while ((ret = tbx_itr_next(fp, tbx, itr, &str)) >= 0)
                {
                    if ( reg_idx && !regidx_overlap(reg_idx,seq[itr->curr_tid],itr->curr_beg,itr->curr_end-1, NULL) ) continue;
                    //CW: TODO: for 2ndary filtering:
                    //want to get current offset [current block "ID"]: itr->off[itr->i] or is it itr->curr_off?
                    //then translate the current block offset into a block ID (starts at 0) that can be used to get block column min/max values
                    //in 2ndary index
                    if (!found) {
                        if (args->separate_regs) printf("%c%s\n", conf->meta_char, regs[i]);
                        found = 1;
                    }
                    if (puts(str.s) < 0) {
                        RELEASE_TPOOL(tpool.pool);
                        error_errno("Failed to write to stdout");
                    }
                }
                if (ret < -1) {
                    RELEASE_TPOOL(tpool.pool);
                    error_errno("Reading \"%s\" failed", fname);
                }
                tbx_itr_destroy(itr);
            }
            free(seq);
        }
        free(str.s);
        tbx_destroy(tbx);
    }
    else if ( format==bam ) {
        RELEASE_TPOOL(tpool.pool);
        error("Please use \"samtools view\" for querying BAM files.\n");
    }

    if ( reg_idx ) regidx_destroy(reg_idx);
    if ( hts_close(fp) ) {
        RELEASE_TPOOL(tpool.pool);
        error_errno("hts_close returned non-zero status: %s", fname);
    }

    for (i=0; i<nregs; i++) free(regs[i]);
    free(regs);
    RELEASE_TPOOL(tpool.pool);
    return 0;
}
static int query_chroms(char *fname, int download)
{
    const char **seq;
    int i, nseq, ftype = file_type(fname);
    if ( ftype & IS_TXT || !ftype )
    {
        tbx_t *tbx = tbx_index_load3(fname, NULL, download ? HTS_IDX_SAVE_REMOTE : 0);
        if ( !tbx ) error_errno("Could not load .tbi index of %s", fname);
        seq = tbx_seqnames(tbx, &nseq);
        if (!seq) error_errno("Couldn't get list of sequence names");
        for (i=0; i<nseq; i++) {
            if (printf("%s\n", seq[i]) < 0)
                error_errno("Couldn't write to stdout");
        }
        free(seq);
        tbx_destroy(tbx);
    }
    else if ( ftype==IS_BCF )
    {
        htsFile *fp = hts_open(fname,"r");
        if ( !fp ) error_errno("Could not open \"%s\"", fname);
        bcf_hdr_t *hdr = bcf_hdr_read(fp);
        if ( !hdr ) error_errno("Could not read the header: \"%s\"", fname);
        hts_close(fp);
        hts_idx_t *idx = bcf_index_load3(fname, NULL, download ? HTS_IDX_SAVE_REMOTE : 0);
        if ( !idx ) error_errno("Could not load .csi index of \"%s\"", fname);
        seq = bcf_index_seqnames(idx, hdr, &nseq);
        if (!seq) error_errno("Couldn't get list of sequence names");
        for (i=0; i<nseq; i++) {
            if (printf("%s\n", seq[i]) < 0)
                error_errno("Couldn't write to stdout");
        }
        free(seq);
        bcf_hdr_destroy(hdr);
        hts_idx_destroy(idx);
    }
    else if ( ftype==IS_BAM )   // todo: BAM
        error("BAM: todo\n");
    return 0;
}

/*
 * dump_blocks() -- scan a BGZF-compressed, tabix-indexed file and:
 *
 *  1. Print a TSV summary to stdout for every compressed block:
 *
 *       compressed_offset  seq(s)  coord_start  coord_end  n_lines  uncompressed_bytes
 *
 *     followed by per-column statistics (one row per column per block):
 *
 *       col_<N>  <type>  <min>  <max>
 *
 *  2. Write a binary secondary index file (<fname>.sidx) with the
 *     layout described below.
 *
 * ── Secondary index file format (.sidx) ─────────────────────────────────────
 *
 *  Field 0  magic[4]                  char[4]      "TBA\x02"
 *  Field 1  n_columns                 uint32_t     # of tab-separated columns
 *  Field 2a column_min_types[ncols]   char[ncols]  type code for each col's min
 *  Field 2b column_max_types[ncols]   char[ncols]  type code for each col's max
 *  Field 3  n_blocks                  uint64_t     # of BGZF blocks
 *  Field 4  block_record_length       uint64_t     byte length of one block record
 *  Field 5  offsets[n_blocks]         uint64_t[]   byte offset of each block record
 *                                                  within the block_records section
 *  Field 6  block_records             char[]       n_blocks * block_record_length bytes
 *
 * Type codes (stored as a single unsigned byte):
 *   0 = string / N/A   1 = char (boolean 0|1)   2 = int16_t   3 = uint16_t
 *   4 = int32_t         5 = uint32_t              6 = int64_t   7 = uint64_t
 *   8 = float           9 = double
 *
 * block_record_length:  for each column, the number of bytes reserved is
 *   sizeof(min_value) + sizeof(max_value).
 *   Exception: if BOTH the min-type and max-type for a column are 0 (string/NA)
 *   only 1 byte is reserved for that column (set to 0x00) to save space.
 *
 * "compressed_offset" is the byte position of the block in the .gz file
 * (upper 48 bits of the BGZF virtual offset).
 *
 * Column type inference (most-restrictive type that fits every value seen):
 *   boolean  – every value is "0" or "1" (single byte)
 *   int16    – fits in [-32768, 32767]
 *   uint16   – fits in [0, 65535]   (and not boolean)
 *   int32    – fits in [-2^31, 2^31-1]
 *   uint32   – fits in [0, 2^32-1]  (and not smaller)
 *   int64    – fits in [-2^63, 2^63-1]
 *   uint64   – fits in [0, 2^64-1]  (and not smaller signed)
 *   float    – parses as floating-point with <=7 significant decimal digits
 *   double   – parses as floating-point with >7 significant decimal digits
 *   string   – anything else
 */

/* ── column-statistics helpers ───────────────────────────────────────────── */

/* Ordered type ranks (lower = more restrictive). */
typedef enum {
    CT_BOOL   = 0,
    CT_INT16  = 1,
    CT_UINT16 = 2,
    CT_INT32  = 3,
    CT_UINT32 = 4,
    CT_INT64  = 5,
    CT_UINT64 = 6,
    CT_FLOAT  = 7,
    CT_DOUBLE = 8,
    CT_STRING = 9
} ColType;

static const char *coltype_name(ColType t) {
    switch (t) {
        case CT_BOOL:   return "boolean";
        case CT_INT16:  return "int16";
        case CT_UINT16: return "uint16";
        case CT_INT32:  return "int32";
        case CT_UINT32: return "uint32";
        case CT_INT64:  return "int64";
        case CT_UINT64: return "uint64";
        case CT_FLOAT:  return "float";
        case CT_DOUBLE: return "double";
        default:        return "string";
    }
}

/*
 * Map ColType to the binary type-code byte used in the .sidx format.
 *   CT_BOOL   → 1  (char/single-byte boolean)
 *   CT_INT16  → 2
 *   CT_UINT16 → 3
 *   CT_INT32  → 4
 *   CT_UINT32 → 5
 *   CT_INT64  → 6
 *   CT_UINT64 → 7
 *   CT_FLOAT  → 8
 *   CT_DOUBLE → 9
 *   CT_STRING → 0  (string / N/A)
 */
static uint8_t coltype_to_sidx(ColType t) {
    /* The ColType enum values 0..8 map to sidx codes 1..9 for all numeric
     * types; CT_STRING (9) maps to sidx code 0. */
    if (t == CT_STRING) return 0;
    return (uint8_t)(t + 1);   /* CT_BOOL=0 → 1, CT_DOUBLE=8 → 9 */
}

/* Return the byte width of a value whose sidx type code is tc. */
static size_t sidx_type_size(uint8_t tc) {
    switch (tc) {
        case 0: return 0;               /* string/NA – no value bytes (handled specially) */
        case 1: return sizeof(uint8_t);  /* char / boolean                                */
        case 2: return sizeof(int16_t);
        case 3: return sizeof(uint16_t);
        case 4: return sizeof(int32_t);
        case 5: return sizeof(uint32_t);
        case 6: return sizeof(int64_t);
        case 7: return sizeof(uint64_t);
        case 8: return sizeof(float);
        case 9: return sizeof(double);
        default: return 0;
    }
}

/* Per-column state for one block. */
typedef struct {
    ColType  type;       /* narrowest type consistent with all values seen    */
    /* numeric min/max stored as the widest possible representations          */
    int64_t  imin;       /* used when type <= CT_INT64                        */
    uint64_t umax;       /* used when type == CT_UINT64 or for unsigned bound */
    double   dmin;       /* used when type == CT_FLOAT or CT_DOUBLE           */
    double   dmax;
    int      seen;       /* 1 once the first value for this column is seen    */
} ColStat;

/* Count decimal significant digits (ignoring leading zeros, sign, decimal
 * point).  Returns the number of significant digits found, capped at 8.
 * We use this to decide float vs double: >7 sig-digits → double. */
static int count_sig_digits(const char *s) {
    int sig = 0, in_sig = 0;
    for (; *s && sig < 8; s++) {
        if (*s == '-' || *s == '+' || *s == '.' || *s == 'e' || *s == 'E')
            continue;
        if (*s >= '0' && *s <= '9') {
            if (*s != '0') in_sig = 1;
            if (in_sig) sig++;
        }
    }
    return sig;
}

/*
 * Classify a single field value and widen *stat if necessary.
 * On the very first call (stat->seen == 0) the type is set to the narrowest
 * possible; subsequent calls only widen.
 */
static void col_update(ColStat *stat, const char *field, size_t len) {
    char *end;

    /* ── Attempt integer parse ─────────────────────────────────────── */
    errno = 0;
    uint64_t uval = (uint64_t)strtoull(field, &end, 10);
    int is_uint = (errno == 0 && end == field + len);

    errno = 0;
    int64_t  ival = (int64_t)strtoll(field, &end, 10);
    int is_int  = (errno == 0 && end == field + len);

    /* ── Determine value's natural type ───────────────────────────── */
    ColType vtype;
    if (is_int || is_uint) {
        if (is_int && ival >= 0 && ival <= 1 && len == 1) {
            vtype = CT_BOOL;
        } else if (is_int && ival >= -32768 && ival <= 32767) {
            vtype = CT_INT16;
        } else if (is_uint && uval <= 65535) {
            vtype = CT_UINT16;
        } else if (is_int && ival >= (int64_t)INT32_MIN && ival <= (int64_t)INT32_MAX) {
            vtype = CT_INT32;
        } else if (is_uint && uval <= (uint64_t)UINT32_MAX) {
            vtype = CT_UINT32;
        } else if (is_int) {
            vtype = CT_INT64;
        } else {
            vtype = CT_UINT64;
        }
    } else {
        errno = 0;
        double dval = strtod(field, &end);
        (void)dval;
        if (errno == 0 && end == field + len && end != field) {
            int nsig = count_sig_digits(field);
            vtype = (nsig <= 7) ? CT_FLOAT : CT_DOUBLE;
        } else {
            vtype = CT_STRING;
        }
    }

    /* ── Widen accumulated type ───────────────────────────────────── */
    if (!stat->seen) {
        stat->type = vtype;
    } else {
        if (vtype != stat->type) {
            if (vtype <= CT_UINT64 && stat->type <= CT_UINT64) {
                if (vtype > stat->type) stat->type = vtype;
            } else if (vtype <= CT_UINT64 && stat->type >= CT_FLOAT) {
                stat->type = CT_DOUBLE;
            } else if (stat->type <= CT_UINT64 && vtype >= CT_FLOAT) {
                stat->type = CT_DOUBLE;
            } else if (vtype >= CT_FLOAT && vtype <= CT_DOUBLE &&
                       stat->type >= CT_FLOAT && stat->type <= CT_DOUBLE) {
                if (vtype > stat->type) stat->type = vtype;
            } else {
                stat->type = CT_STRING;
            }
        }
    }

    /* ── Update min / max ─────────────────────────────────────────── */
    if (stat->type == CT_STRING) {
        stat->seen = 1;
        return;
    }

    if (stat->type <= CT_INT64) {
        end = NULL; errno = 0;
        int64_t v = (int64_t)strtoll(field, &end, 10);
        if (errno || end == field) { stat->type = CT_STRING; stat->seen = 1; return; }
        if (!stat->seen || v < stat->imin) stat->imin = v;
        if (!stat->seen || (int64_t)stat->umax < v) stat->umax = (uint64_t)v;
    } else if (stat->type == CT_UINT32 || stat->type == CT_UINT64) {
        end = NULL; errno = 0;
        uint64_t v = (uint64_t)strtoull(field, &end, 10);
        if (errno || end == field) { stat->type = CT_STRING; stat->seen = 1; return; }
        if (!stat->seen || (int64_t)v < stat->imin) stat->imin = (int64_t)v;
        if (!stat->seen || v > stat->umax) stat->umax = v;
    } else {
        end = NULL; errno = 0;
        double v = strtod(field, &end);
        if (errno || end == field) { stat->type = CT_STRING; stat->seen = 1; return; }
        if (!stat->seen || v < stat->dmin) stat->dmin = v;
        if (!stat->seen || v > stat->dmax) stat->dmax = v;
    }

    stat->seen = 1;
}

/* Initial / reset value for a ColStat. */
static ColStat colstat_init(void) {
    ColStat s;
    memset(&s, 0, sizeof(s));
    return s;
}

/*
 * Serialise the min value of *cs into buf[] according to sidx type code tc.
 * Returns the number of bytes written (== sidx_type_size(tc)).
 * Caller must ensure buf has enough room.
 */
static size_t sidx_write_min(uint8_t tc, const ColStat *cs, uint8_t *buf) {
    switch (tc) {
        case 1: { uint8_t  v = (uint8_t)cs->imin;  memcpy(buf, &v, 1); return 1; }
        case 2: { int16_t  v = (int16_t)cs->imin;  memcpy(buf, &v, 2); return 2; }
        case 3: { uint16_t v = (uint16_t)(uint64_t)cs->imin; memcpy(buf, &v, 2); return 2; }
        case 4: { int32_t  v = (int32_t)cs->imin;  memcpy(buf, &v, 4); return 4; }
        case 5: { uint32_t v = (uint32_t)(uint64_t)cs->imin; memcpy(buf, &v, 4); return 4; }
        case 6: { int64_t  v = cs->imin;            memcpy(buf, &v, 8); return 8; }
        case 7: { uint64_t v = (uint64_t)cs->imin;  memcpy(buf, &v, 8); return 8; }
        case 8: { float    v = (float)cs->dmin;     memcpy(buf, &v, 4); return 4; }
        case 9: { double   v = cs->dmin;             memcpy(buf, &v, 8); return 8; }
        default: return 0;
    }
}

/*
 * Serialise the max value of *cs into buf[] according to sidx type code tc.
 * Returns the number of bytes written (== sidx_type_size(tc)).
 */
static size_t sidx_write_max(uint8_t tc, const ColStat *cs, uint8_t *buf) {
    switch (tc) {
        case 1: { uint8_t  v = (uint8_t)cs->umax;   memcpy(buf, &v, 1); return 1; }
        case 2: { int16_t  v = (int16_t)(int64_t)cs->umax;  memcpy(buf, &v, 2); return 2; }
        case 3: { uint16_t v = (uint16_t)cs->umax;  memcpy(buf, &v, 2); return 2; }
        case 4: { int32_t  v = (int32_t)(int64_t)cs->umax;  memcpy(buf, &v, 4); return 4; }
        case 5: { uint32_t v = (uint32_t)cs->umax;  memcpy(buf, &v, 4); return 4; }
        case 6: { int64_t  v = (int64_t)cs->umax;   memcpy(buf, &v, 8); return 8; }
        case 7: { uint64_t v = cs->umax;             memcpy(buf, &v, 8); return 8; }
        case 8: { float    v = (float)cs->dmax;      memcpy(buf, &v, 4); return 4; }
        case 9: { double   v = cs->dmax;             memcpy(buf, &v, 8); return 8; }
        default: return 0;
    }
}

/* ── Saved per-block data (accumulated in memory for the two-pass write) ─── */

/*
 * We cannot write the .sidx header until the whole file is scanned (we need
 * n_columns, n_blocks, and block_record_length up-front in the file).
 * Strategy: accumulate each block's serialised record into a growable byte
 * buffer (sidx_data / sidx_data_len) and also record its starting offset
 * in that buffer.  At the end, write the full .sidx in one pass.
 */

/* ── Main dump_blocks function ────────────────────────────────────────────── */

static int dump_blocks(const char *fname, int download)
{
    int ftype = file_type(fname);
    if (!(ftype & IS_TXT) && ftype != 0) {
        fprintf(stderr, "[tabix] --dump-blocks only supports text (TBI) indexed files\n");
        return 1;
    }

    /* Load the tabix index so we can resolve tid -> name. */
    tbx_t *tbx = tbx_index_load3(fname, NULL, download ? HTS_IDX_SAVE_REMOTE : 0);
    if (!tbx) {
        fprintf(stderr, "[tabix] Could not load .tbi index of %s\n", fname);
        return 1;
    }

    int nseq = 0;
    const char **seqnames = tbx_seqnames(tbx, &nseq);
    if (!seqnames) {
        fprintf(stderr, "[tabix] Could not get sequence names\n");
        tbx_destroy(tbx);
        return 1;
    }

    /* Open the compressed data file for scanning. */
    BGZF *fp = bgzf_open(fname, "r");
    if (!fp) {
        fprintf(stderr, "[tabix] Could not open %s\n", fname);
        free(seqnames);
        tbx_destroy(tbx);
        return 1;
    }

    /* -- per-block accumulators -- */
    int64_t  blk_compressed_off = 0;
    int64_t  blk_coord_beg      = -1;
    int64_t  blk_coord_end      = -1;
    int64_t  blk_nlines         = 0;
    int64_t  blk_uncomp_bytes   = 0;
    int64_t  last_utell = 0;
    int64_t  new_utell  = 0;

    int      blk_ntids          = 0;
    int      blk_tids_cap       = 16;
    int     *blk_tids           = malloc(blk_tids_cap * sizeof(int));
    if (!blk_tids) {
        fprintf(stderr, "[tabix] out of memory\n");
        bgzf_close(fp); free(seqnames); tbx_destroy(tbx);
        return 1;
    }

    /* -- per-column statistics (across the whole file, grown as needed) -- */
    int       ncols         = 0;
    int       ncols_cap     = 32;
    ColStat  *blk_col       = malloc(ncols_cap * sizeof(ColStat));
    ColType  *file_col_type = malloc(ncols_cap * sizeof(ColType));
    if (!blk_col || !file_col_type) {
        fprintf(stderr, "[tabix] out of memory\n");
        free(blk_tids); free(blk_col); free(file_col_type);
        bgzf_close(fp); free(seqnames); tbx_destroy(tbx);
        return 1;
    }

    /*
     * Secondary index in-memory accumulation.
     *
     * sidx_records      – raw bytes of all block records concatenated
     * sidx_records_len  – bytes used so far
     * sidx_records_cap  – allocated capacity
     * sidx_offsets      – offset of each block record within sidx_records[]
     * sidx_nblocks      – number of block records stored so far
     * sidx_nblocks_cap  – capacity of sidx_offsets[]
     *
     * block_record_length is not known until after the first block is flushed
     * (because ncols is discovered incrementally).  Once ncols stabilises we
     * compute it once and cache it in sidx_rec_len.  If ncols grows after
     * that point existing records are retroactively zero-extended.
     */
    uint8_t  *sidx_records     = NULL;
    size_t    sidx_records_len = 0;
    size_t    sidx_records_cap = 0;
    uint64_t *sidx_offsets     = NULL;
    uint64_t  sidx_nblocks     = 0;
    uint64_t  sidx_nblocks_cap = 0;

    /* Helper: grow column arrays to at least 'need' columns. */
#define ENSURE_COLS(need) do { \
    if ((need) > ncols_cap) { \
        int _newcap = ncols_cap; \
        while (_newcap < (need)) _newcap *= 2; \
        ColStat  *_bc  = realloc(blk_col,       _newcap * sizeof(ColStat));  \
        ColType  *_fct = realloc(file_col_type,  _newcap * sizeof(ColType)); \
        if (!_bc || !_fct) { \
            fprintf(stderr, "[tabix] out of memory\n"); \
            free(_bc ? _bc : blk_col); \
            free(_fct ? _fct : file_col_type); \
            blk_col = NULL; file_col_type = NULL; \
            goto done; \
        } \
        blk_col = _bc; file_col_type = _fct; \
        for (int _ci = ncols_cap; _ci < _newcap; _ci++) { \
            blk_col[_ci]       = colstat_init(); \
            file_col_type[_ci] = CT_BOOL; \
        } \
        ncols_cap = _newcap; \
    } \
} while (0)

    /* Helper: reset per-block column stats (keep ncols and file_col_type). */
#define RESET_BLK_COLS() do { \
    for (int _ci = 0; _ci < ncols; _ci++) \
        blk_col[_ci] = colstat_init(); \
} while (0)

    /* Print TSV header */
    printf("compressed_offset\tseq\tcoord_start\tcoord_end\tn_lines\tuncompressed_bytes\n");

    kstring_t str       = {0, 0, 0};
    char     *linebuf   = NULL;
    size_t    linebuf_cap = 0;
    int64_t   lineno    = 0;

    /*
     * FLUSH_BLOCK() – called at every block boundary and at end-of-file.
     *
     * Responsibilities:
     *   A. Print the human-readable TSV summary line + per-column stats.
     *   B. Merge this block's ColType into the file-level file_col_type[].
     *   C. Serialise the block record into sidx_records[] and record its
     *      starting offset in sidx_offsets[].
     *
     * The .sidx block record for one block has exactly block_record_length
     * bytes laid out as:
     *
     *   for each column ci:
     *     if (min_tc[ci] == 0 && max_tc[ci] == 0):
     *       1 zero byte
     *     else:
     *       min_value  (sidx_type_size(min_tc[ci]) bytes)
     *       max_value  (sidx_type_size(max_tc[ci]) bytes)
     *
     * where min_tc[ci] = coltype_to_sidx( blk_col[ci].type ) for the min
     * and max_tc[ci]   = coltype_to_sidx( blk_col[ci].type ) for the max.
     * (Currently min and max share the same inferred type for a column.)
     *
     * block_record_length is fixed for the whole file once ncols is known.
     * It is computed lazily on the first flush (or recomputed if ncols grows,
     * though in practice ncols stabilises on the first data block).
     */
#define FLUSH_BLOCK() do { \
    if (blk_nlines > 0) { \
        /* ── A: human-readable TSV output ── */ \
        if (blk_ntids == 0) { \
            printf("%"PRId64"\t.\t-1\t-1\t%"PRId64"\t%"PRId64"\n", \
                   blk_compressed_off, blk_nlines, blk_uncomp_bytes); \
        } else { \
            printf("%"PRId64"\t", blk_compressed_off); \
            for (int _i = 0; _i < blk_ntids; _i++) { \
                if (_i) putchar(','); \
                int _tid = blk_tids[_i]; \
                if (_tid >= 0 && _tid < nseq) \
                    fputs(seqnames[_tid], stdout); \
                else \
                    printf("tid%d", _tid); \
            } \
            printf("\t%"PRId64"\t%"PRId64"\t%"PRId64"\t%"PRId64"\n", \
                   blk_coord_beg, blk_coord_end, blk_nlines, blk_uncomp_bytes); \
        } \
        for (int _ci = 0; _ci < ncols; _ci++) { \
            ColStat *_cs = &blk_col[_ci]; \
            if (!_cs->seen) { \
                printf("\tcol_%d\t(no data)\t(no data)\n", _ci); \
            } else if (_cs->type == CT_STRING) { \
                printf("\tcol_%d\t%s\t(N/A)\t(N/A)\n", _ci, coltype_name(_cs->type)); \
            } else if (_cs->type <= CT_INT64) { \
                printf("\tcol_%d\t%s\t%"PRId64"\t%"PRId64"\n", \
                       _ci, coltype_name(_cs->type), _cs->imin, (int64_t)_cs->umax); \
            } else if (_cs->type == CT_UINT32 || _cs->type == CT_UINT64) { \
                printf("\tcol_%d\t%s\t%"PRIu64"\t%"PRIu64"\n", \
                       _ci, coltype_name(_cs->type), (uint64_t)_cs->imin, _cs->umax); \
            } else { \
                printf("\tcol_%d\t%s\t%g\t%g\n", \
                       _ci, coltype_name(_cs->type), _cs->dmin, _cs->dmax); \
            } \
            /* ── B: merge into file-level type ── */ \
            if (_cs->seen && file_col_type[_ci] != CT_STRING) { \
                ColType _bt = _cs->type, _ft = file_col_type[_ci]; \
                if (_bt != _ft) { \
                    if (_bt <= CT_UINT64 && _ft <= CT_UINT64) { \
                        if (_bt > _ft) file_col_type[_ci] = _bt; \
                    } else if ((_bt <= CT_UINT64 && _ft >= CT_FLOAT) || \
                               (_ft <= CT_UINT64 && _bt >= CT_FLOAT)) { \
                        file_col_type[_ci] = CT_DOUBLE; \
                    } else if (_bt >= CT_FLOAT && _bt <= CT_DOUBLE && \
                               _ft >= CT_FLOAT && _ft <= CT_DOUBLE) { \
                        if (_bt > _ft) file_col_type[_ci] = _bt; \
                    } else { \
                        file_col_type[_ci] = CT_STRING; \
                    } \
                } \
            } \
        } \
        /* ── C: serialise block record into sidx_records[] ── */ \
        /* \
         * Compute the record for this block into a local scratch buffer,    \
         * then append to sidx_records[].                                    \
         *                                                                   \
         * For each column ci:                                                \
         *   min_tc = max_tc = coltype_to_sidx(blk_col[ci].type)             \
         *   (the current implementation uses the same inferred ColType for  \
         *    both the minimum and maximum of a column within a block)        \
         *   if min_tc == 0 && max_tc == 0: write 1 zero byte               \
         *   else: write min_value bytes, then max_value bytes               \
         */ \
        { \
            /* Compute record length for this block snapshot.                \
             * If ncols is 0 (header-only block) we still write a 0-byte    \
             * record so the offset table stays aligned.                     */ \
            size_t _rlen = 0; \
            for (int _ci = 0; _ci < ncols; _ci++) { \
                uint8_t _mtc = coltype_to_sidx(blk_col[_ci].type); \
                if (_mtc == 0) { \
                    _rlen += 1; /* both NA → 1 sentinel byte */ \
                } else { \
                    _rlen += sidx_type_size(_mtc) * 2; /* min + max */ \
                } \
            } \
            /* Grow sidx_records[] if needed. */ \
            if (sidx_records_len + _rlen > sidx_records_cap) { \
                size_t _newcap = sidx_records_cap ? sidx_records_cap * 2 : 65536; \
                while (_newcap < sidx_records_len + _rlen) _newcap *= 2; \
                uint8_t *_nb = realloc(sidx_records, _newcap); \
                if (!_nb) { \
                    fprintf(stderr, "[tabix] out of memory (sidx records)\n"); \
                    goto done; \
                } \
                sidx_records = _nb; \
                sidx_records_cap = _newcap; \
            } \
            /* Grow sidx_offsets[] if needed. */ \
            if (sidx_nblocks + 1 > sidx_nblocks_cap) { \
                uint64_t _newcap = sidx_nblocks_cap ? sidx_nblocks_cap * 2 : 1024; \
                uint64_t *_no = realloc(sidx_offsets, (size_t)_newcap * sizeof(uint64_t)); \
                if (!_no) { \
                    fprintf(stderr, "[tabix] out of memory (sidx offsets)\n"); \
                    goto done; \
                } \
                sidx_offsets = _no; \
                sidx_nblocks_cap = _newcap; \
            } \
            /* Record this block's offset within the records section. */ \
            sidx_offsets[sidx_nblocks] = (uint64_t)sidx_records_len; \
            sidx_nblocks++; \
            /* Serialise each column's min/max. */ \
            uint8_t *_p = sidx_records + sidx_records_len; \
            for (int _ci = 0; _ci < ncols; _ci++) { \
                ColStat *_cs2 = &blk_col[_ci]; \
                uint8_t _tc = coltype_to_sidx(_cs2->type); \
                if (_tc == 0) { \
                    *_p++ = 0x00; /* sentinel for string/NA column */ \
                } else { \
                    _p += sidx_write_min(_tc, _cs2, _p); \
                    _p += sidx_write_max(_tc, _cs2, _p); \
                } \
            } \
            sidx_records_len += _rlen; \
        } \
    } \
} while (0)

    /* Initial virtual offset before the first read. */
    uint64_t voff_after = 0;
    blk_compressed_off = 0;

    int ret;
    while ((ret = bgzf_getline(fp, '\n', &str)) >= 0) {
        ++lineno;
        uint64_t voff_now = (uint64_t)bgzf_tell(fp);
        int64_t  cur_blk  = (int64_t)(voff_now >> 16);

        int crossed = (cur_blk != (int64_t)(voff_after >> 16));

        int is_data = 0;
        tbx_intv_t intv;
        if (str.l && str.s[0] != tbx->conf.meta_char && lineno > tbx->conf.line_skip) {
            if (tbx_parse1(&tbx->conf, str.l, str.s, &intv) == 0) {
                char save = *intv.se;
                *intv.se = '\0';
                int tid = tbx_name2id(tbx, intv.ss);
                *intv.se = save;

                if (crossed) {
                    new_utell = bgzf_utell(fp);
                    blk_uncomp_bytes = new_utell - last_utell;
                    last_utell = new_utell;
                    FLUSH_BLOCK();
                    blk_compressed_off = cur_blk;
                    blk_coord_beg    = -1;
                    blk_coord_end    = -1;
                    blk_nlines       = 0;
                    blk_ntids        = 0;
                    RESET_BLK_COLS();
                }

                blk_nlines++;
                is_data = 1;

                if (blk_coord_beg < 0 || intv.beg < blk_coord_beg)
                    blk_coord_beg = intv.beg;
                if (blk_coord_end < 0 || intv.end > blk_coord_end)
                    blk_coord_end = intv.end;

                if (tid >= 0) {
                    int found = 0;
                    for (int i = 0; i < blk_ntids; i++) {
                        if (blk_tids[i] == tid) { found = 1; break; }
                    }
                    if (!found) {
                        if (blk_ntids == blk_tids_cap) {
                            blk_tids_cap *= 2;
                            int *tmp = realloc(blk_tids, blk_tids_cap * sizeof(int));
                            if (tmp == NULL) {
                                fprintf(stderr, "[tabix] out of memory\n");
                                goto done;
                            }
                            blk_tids = tmp;
                        }
                        blk_tids[blk_ntids++] = tid;
                    }
                }

                /* ── Tokenise line and update per-column statistics ────── */
                if (str.l + 1 > linebuf_cap) {
                    size_t newcap = (str.l + 1) * 2;
                    char *nb = realloc(linebuf, newcap);
                    if (!nb) {
                        fprintf(stderr, "[tabix] out of memory\n");
                        goto done;
                    }
                    linebuf = nb;
                    linebuf_cap = newcap;
                }
                memcpy(linebuf, str.s, str.l + 1);

                int col_idx = 0;
                char *p = linebuf;
                char *line_end = linebuf + str.l;
                while (p <= line_end) {
                    char *tab = memchr(p, '\t', (size_t)(line_end - p));
                    char *field_end = tab ? tab : line_end;
                    *field_end = '\0';

                    ENSURE_COLS(col_idx + 1);
                    if (col_idx >= ncols) {
                        for (int ci = ncols; ci <= col_idx; ci++) {
                            blk_col[ci]       = colstat_init();
                            file_col_type[ci] = CT_BOOL;
                        }
                        ncols = col_idx + 1;
                    }

                    col_update(&blk_col[col_idx], p, (size_t)(field_end - p));

                    if (!tab) break;
                    p = tab + 1;
                    col_idx++;
                }
            }
        }

        if (!is_data && crossed) {
            blk_uncomp_bytes = bgzf_utell(fp) - last_utell;
            FLUSH_BLOCK();
            blk_compressed_off = cur_blk;
            blk_coord_beg    = -1;
            blk_coord_end    = -1;
            blk_nlines       = 1;
            blk_ntids        = 0;
            RESET_BLK_COLS();
            /* Claude added this mistakenly thinking the else if could be processed in
             * *addition* to the previous if clause, which it can't – but does not harm,
             * so leaving it here as a lesson (the !crossed is implied by the above
             * if statement's check of crossed). */
        } else if (!is_data && !crossed) {
            blk_nlines++;
        }

        voff_after = voff_now;
    }

    /* Flush the final block. */
    blk_uncomp_bytes = bgzf_utell(fp) - last_utell;
    FLUSH_BLOCK();

    /* ── Print whole-file column type summary ───────────────────────────── */
    if (ncols > 0) {
        printf("\n# File-level column type summary (%d column(s)):\n", ncols);
        for (int ci = 0; ci < ncols; ci++)
            printf("#   col_%d: %s\n", ci, coltype_name(file_col_type[ci]));
    }

    /* ── Write the binary secondary index (.sidx) ───────────────────────── */
    /*
     * Now that we know ncols, n_blocks, and every block's record bytes, we
     * can compute block_record_length and write the complete .sidx file.
     *
     * block_record_length: fixed across all blocks.  For each column:
     *   • if both min_type and max_type codes are 0 → 1 byte
     *   • else → sizeof(min_value) + sizeof(max_value)
     *
     * The file-level type in file_col_type[] is the widest type seen for
     * that column across all blocks; that is what we store in the header
     * arrays column_min_types[] and column_max_types[].
     * Each per-block record, however, uses the per-block inferred type (which
     * may be narrower); it is already serialised in sidx_records[].
     *
     * NOTE: because we serialise block records eagerly in FLUSH_BLOCK using
     * the per-block type, but the header stores the file-level (widest) type,
     * a reader must use the per-block record's own type codes to interpret
     * the payload.  The header types are a "worst-case" hint only.
     *
     * block_record_length is derived from the file-level types (the widest
     * possible per column) so every block record is guaranteed to fit.
     * Blocks that saw only string/NA for a column write 1 byte; all others
     * write sizeof(min)+sizeof(max) bytes for the file-level type of that
     * column.  This means the per-block record length is identical to the
     * file-level record length – the lengths are consistent.
     */
    {
        /* Build the column type-code arrays using file-level (widest) types. */
        uint8_t *col_min_tc = malloc(ncols > 0 ? (size_t)ncols : 1);
        uint8_t *col_max_tc = malloc(ncols > 0 ? (size_t)ncols : 1);
        if (!col_min_tc || !col_max_tc) {
            fprintf(stderr, "[tabix] out of memory (sidx header)\n");
            free(col_min_tc); free(col_max_tc);
            goto done;
        }
        uint64_t rec_len = 0;
        for (int ci = 0; ci < ncols; ci++) {
            uint8_t tc = coltype_to_sidx(file_col_type[ci]);
            col_min_tc[ci] = tc;
            col_max_tc[ci] = tc;
            if (tc == 0) {
                rec_len += 1;               /* both NA → 1 sentinel byte */
            } else {
                rec_len += sidx_type_size(tc) * 2; /* min + max */
            }
        }

        /* Build the output filename: <fname>.sidx */
        size_t fnlen = strlen(fname);
        char *sidx_fname = malloc(fnlen + 6);
        if (!sidx_fname) {
            fprintf(stderr, "[tabix] out of memory (sidx filename)\n");
            free(col_min_tc); free(col_max_tc);
            goto done;
        }
        memcpy(sidx_fname, fname, fnlen);
        memcpy(sidx_fname + fnlen, ".sidx", 6); /* includes NUL */

        FILE *sf = fopen(sidx_fname, "wb");
        if (!sf) {
            fprintf(stderr, "[tabix] Could not open %s for writing: %s\n",
                    sidx_fname, strerror(errno));
            free(col_min_tc); free(col_max_tc); free(sidx_fname);
            goto done;
        }

        /* Helper: write exactly n bytes or abort. */
#define SIDX_WRITE(ptr, n) do { \
    if (fwrite((ptr), 1, (n), sf) != (size_t)(n)) { \
        fprintf(stderr, "[tabix] Write error on %s: %s\n", sidx_fname, strerror(errno)); \
        fclose(sf); free(col_min_tc); free(col_max_tc); free(sidx_fname); \
        goto done; \
    } \
} while (0)

        /* Field 0: magic "TBA\x02" */
        const char magic[4] = {'T', 'B', 'A', '\x02'};
        SIDX_WRITE(magic, 4);

        /* Field 1: n_columns (uint32_t) */
        {
            uint32_t nc = (uint32_t)ncols;
            SIDX_WRITE(&nc, 4);
        }

        /* Field 2a: column_min_types[ncols] */
        if (ncols > 0) SIDX_WRITE(col_min_tc, (size_t)ncols);

        /* Field 2b: column_max_types[ncols] */
        if (ncols > 0) SIDX_WRITE(col_max_tc, (size_t)ncols);

        /* Field 3: n_blocks (uint64_t) */
        SIDX_WRITE(&sidx_nblocks, 8);

        /* Field 4: block_record_length (uint64_t) */
        SIDX_WRITE(&rec_len, 8);

        /* Field 5: block_record_uncompressed_byte_offsets (uint64_t[n_blocks]) */
        if (sidx_nblocks > 0)
            SIDX_WRITE(sidx_offsets, (size_t)(sidx_nblocks * sizeof(uint64_t)));

        /* Field 6: block_records */
        if (sidx_records_len > 0)
            SIDX_WRITE(sidx_records, sidx_records_len);

#undef SIDX_WRITE

        fclose(sf);
        fprintf(stderr, "[tabix] wrote secondary index: %s "
                "(%d col(s), %"PRIu64" block(s), %"PRIu64" bytes/record)\n",
                sidx_fname, ncols, sidx_nblocks, rec_len);

        free(col_min_tc);
        free(col_max_tc);
        free(sidx_fname);
    }

#undef FLUSH_BLOCK
#undef RESET_BLK_COLS
#undef ENSURE_COLS

done:
    free(sidx_records);
    free(sidx_offsets);
    free(linebuf);
    free(blk_col);
    free(file_col_type);
    free(str.s);
    free(blk_tids);
    free(seqnames);
    bgzf_close(fp);
    tbx_destroy(tbx);
    return 0;
}

int reheader_file(const char *fname, const char *header, int ftype, tbx_conf_t *conf, int threads)
{
    hts_tpool *tpool = NULL;
    if (threads >= 1) {
        if (!(tpool = hts_tpool_init(threads))) {
            hts_log_info("Could not initialize thread pool!");
        }
    }
    if ( ftype & IS_TXT || !ftype )
    {
        BGZF *fp = bgzf_open(fname,"r");
        if (!fp) {
            RELEASE_TPOOL(tpool);
            return -1;
        }
        if (bgzf_thread_pool(fp, tpool, 0) < 0) {
            hts_log_info("Could not set thread pool!");
        }
        if (bgzf_read_block(fp) != 0 || !fp->block_length ) {
            RELEASE_TPOOL(tpool);
            return -1;
        }

        char *buffer = fp->uncompressed_block;
        int skip_until = 0;

        // Skip the header: find out the position of the data block
        if ( buffer[0]==conf->meta_char )
        {
            skip_until = 1;
            while (1)
            {
                if ( buffer[skip_until]=='\n' )
                {
                    skip_until++;
                    if ( skip_until>=fp->block_length )
                    {
                        if ( bgzf_read_block(fp) != 0 || !fp->block_length ) {
                            RELEASE_TPOOL(tpool);
                            error("FIXME: No body in the file: %s\n", fname);
                        }
                        skip_until = 0;
                    }
                    // The header has finished
                    if ( buffer[skip_until]!=conf->meta_char ) break;
                }
                skip_until++;
                if ( skip_until>=fp->block_length )
                {
                    if (bgzf_read_block(fp) != 0 || !fp->block_length) {
                        RELEASE_TPOOL(tpool);
                        error("FIXME: No body in the file: %s\n", fname);
                    }
                    skip_until = 0;
                }
            }
        }

        // Output the new header
        FILE *hdr  = fopen(header,"r");
        if ( !hdr ) {
            RELEASE_TPOOL(tpool);
            error("%s: %s", header,strerror(errno));
        }
        const size_t page_size = 32768;
        char *buf = malloc(page_size);
        BGZF *bgzf_out = bgzf_open("-", "w");
        ssize_t nread;

        if (!buf) {
            RELEASE_TPOOL(tpool);
            error("%s\n", strerror(errno));
        }
        if (!bgzf_out) {
            RELEASE_TPOOL(tpool);
            error_errno("Couldn't open output stream");
        }
        if (bgzf_thread_pool(bgzf_out, tpool, 0) < 0) {
            hts_log_info("Could not set thread pool to output file!");
        }
        while ( (nread=fread(buf,1,page_size-1,hdr))>0 )
        {
            if ( nread<page_size-1 && buf[nread-1]!='\n' ) buf[nread++] = '\n';
            if (bgzf_write(bgzf_out, buf, nread) < 0) {
                RELEASE_TPOOL(tpool);
                error_errno("Write error %d", bgzf_out->errcode);
            }
        }
        if ( ferror(hdr) ) {
            RELEASE_TPOOL(tpool);
            error_errno("Failed to read \"%s\"", header);
        }
        if ( fclose(hdr) ) {
            RELEASE_TPOOL(tpool);
            error_errno("Closing \"%s\" failed", header);
        }

        // Output all remaining data read with the header block
        if ( fp->block_length - skip_until > 0 )
        {
            if (bgzf_write(bgzf_out, buffer+skip_until, fp->block_length-skip_until) < 0) {
                RELEASE_TPOOL(tpool);
                error_errno("Write error %d",fp->errcode);
            }
        }
        if (bgzf_flush(bgzf_out) < 0) {
            RELEASE_TPOOL(tpool);
            error_errno("Write error %d", bgzf_out->errcode);
        }

        while (1)
        {
            nread = bgzf_raw_read(fp, buf, page_size);
            if ( nread<=0 ) break;

            int count = bgzf_raw_write(bgzf_out, buf, nread);
            if (count != nread) {
                RELEASE_TPOOL(tpool);
                error_errno("Write failed, wrote %d instead of %d bytes", count,(int)nread);
            }
        }
        if (nread < 0) {
            RELEASE_TPOOL(tpool);
            error_errno("Error reading \"%s\"", fname);
        }
        if (bgzf_close(bgzf_out) < 0) {
            RELEASE_TPOOL(tpool);
            error_errno("Error %d closing output", bgzf_out->errcode);
        }
        if (bgzf_close(fp) < 0) {
            RELEASE_TPOOL(tpool);
            error_errno("Error %d closing \"%s\"", bgzf_out->errcode, fname);
        }
        free(buf);
    }
    else {
        RELEASE_TPOOL(tpool);
        error("todo: reheader BCF, BAM\n");  // BCF is difficult, records contain pointers to the header.
    }
    RELEASE_TPOOL(tpool);
    return 0;
}

static int usage(FILE *fp, int status)
{
    fprintf(fp, "\n");
    fprintf(fp, "Version: %s\n", hts_version());
    fprintf(fp, "Usage:   tabix [OPTIONS] [FILE] [REGION [...]]\n");
    fprintf(fp, "\n");
    fprintf(fp, "Indexing Options:\n");
    fprintf(fp, "   -0, --zero-based           coordinates are zero-based\n");
    fprintf(fp, "   -b, --begin INT            column number for region start [4]\n");
    fprintf(fp, "   -c, --comment CHAR         skip comment lines starting with CHAR [null]\n");
    fprintf(fp, "   -C, --csi                  generate CSI index for VCF (default is TBI)\n");
    fprintf(fp, "   -e, --end INT              column number for region end (if no end, set INT to -b) [5]\n");
    fprintf(fp, "   -f, --force                overwrite existing index without asking\n");
    fprintf(fp, "   -m, --min-shift INT        set minimal interval size for CSI indices to 2^INT [14]\n");
    fprintf(fp, "   -p, --preset STR           gff, bed, sam, vcf, gaf, seqonly\n");
    fprintf(fp, "   -s, --sequence INT         column number for sequence names (suppressed by -p) [1]\n");
    fprintf(fp, "   -S, --skip-lines INT       skip first INT lines [0]\n");
    fprintf(fp, "\n");
    fprintf(fp, "Querying and other options:\n");
    fprintf(fp, "   -h, --print-header         print also the header lines\n");
    fprintf(fp, "   -H, --only-header          print only the header lines\n");
    fprintf(fp, "   -l, --list-chroms          list chromosome names\n");
    fprintf(fp, "       --dump-blocks          print compressed block info (offset, seq(s), coord range, n_lines, uncompressed_bytes)\n");
    fprintf(fp, "   -r, --reheader FILE        replace the header with the content of FILE\n");
    fprintf(fp, "   -R, --regions FILE         restrict to regions listed in the file\n");
    fprintf(fp, "   -T, --targets FILE         similar to -R but streams rather than index-jumps\n");
    fprintf(fp, "   -D                         do not download the index file\n");
    fprintf(fp, "       --cache INT            set cache size to INT megabytes (0 disables) [10]\n");
    fprintf(fp, "       --separate-regions     separate the output by corresponding regions\n");
    fprintf(fp, "       --verbosity INT        set verbosity [3]\n");
    fprintf(fp, "   -@, --threads INT          number of additional threads to use [0]\n");
    fprintf(fp, "\n");
    fprintf(fp, "The 'seqonly' preset (-p seqonly) indexes tab-delimited files by\n");
    fprintf(fp, "sequence/chromosome name only. No start/end coordinate columns are\n");
    fprintf(fp, "needed. Use -s to specify the chromosome column (default: 1).\n");
    fprintf(fp, "Query with a bare name: tabix FILE chr1\n");
    fprintf(fp, "\n");
    return status;
}

int main(int argc, char *argv[])
{
    int c, detect = 1, min_shift = 0, is_force = 0, list_chroms = 0, do_csi = 0, dump_blocks_flag = 0;
    tbx_conf_t conf = tbx_conf_gff;
    char *reheader = NULL;
    args_t args;
    memset(&args,0,sizeof(args_t));
    args.cache_megs = 10;
    args.download_index = 1;
    int32_t new_line_skip = -1;

    static const struct option loptions[] =
    {
        {"help", no_argument, NULL, 2},
        {"regions", required_argument, NULL, 'R'},
        {"targets", required_argument, NULL, 'T'},
        {"csi", no_argument, NULL, 'C'},
        {"zero-based", no_argument, NULL, '0'},
        {"print-header", no_argument, NULL, 'h'},
        {"only-header", no_argument, NULL, 'H'},
        {"begin", required_argument, NULL, 'b'},
        {"comment", required_argument, NULL, 'c'},
        {"end", required_argument, NULL, 'e'},
        {"force", no_argument, NULL, 'f'},
        {"min-shift", required_argument, NULL, 'm'},
        {"preset", required_argument, NULL, 'p'},
        {"sequence", required_argument, NULL, 's'},
        {"skip-lines", required_argument, NULL, 'S'},
        {"list-chroms", no_argument, NULL, 'l'},
        {"reheader", required_argument, NULL, 'r'},
        {"version", no_argument, NULL, 1},
        {"verbosity", required_argument, NULL, 3},
        {"cache", required_argument, NULL, 4},
        {"separate-regions", no_argument, NULL, 5},
        {"threads", required_argument, NULL, '@'},
        {"dump-blocks", no_argument, NULL, 6},
        {NULL, 0, NULL, 0}
    };

    char *tmp;
    while ((c = getopt_long(argc, argv, "hH?0b:c:e:fm:p:s:S:lr:CR:T:D@:", loptions,NULL)) >= 0)
    {
        switch (c)
        {
            case 'R': args.regions_fname = optarg; break;
            case 'T': args.targets_fname = optarg; break;
            case 'C': do_csi = 1; break;
            case 'r': reheader = optarg; break;
            case 'h': args.print_header = 1; break;
            case 'H': args.print_header = 1; args.header_only = 1; break;
            case 'l': list_chroms = 1; break;
            case '0': conf.preset |= TBX_UCSC; detect = 0; break;
            case 'b':
                conf.bc = strtol(optarg,&tmp,10);
                if ( *tmp ) error("Could not parse argument: -b %s\n", optarg);
                detect = 0;
                break;
            case 'e':
                conf.ec = strtol(optarg,&tmp,10);
                if ( *tmp ) error("Could not parse argument: -e %s\n", optarg);
                detect = 0;
                break;
            case 'c': conf.meta_char = *optarg; detect = 0; break;
            case 'f': is_force = 1; break;
            case 'm':
                min_shift = strtol(optarg,&tmp,10);
                if ( *tmp ) error("Could not parse argument: -m %s\n", optarg);
                break;
            case 'p':
                detect = 0;
                if (strcmp(optarg, "gff") == 0) conf = tbx_conf_gff;
                else if (strcmp(optarg, "bed") == 0) conf = tbx_conf_bed;
                else if (strcmp(optarg, "sam") == 0) conf = tbx_conf_sam;
                else if (strcmp(optarg, "vcf") == 0) conf = tbx_conf_vcf;
                else if (strcmp(optarg, "gaf") == 0) conf = tbx_conf_gaf;
                else if (strcmp(optarg, "seqonly") == 0) conf = tbx_conf_seqonly;
                else if (strcmp(optarg, "bcf") == 0) detect = 1; // bcf is autodetected, preset is not needed
                else if (strcmp(optarg, "bam") == 0) detect = 1; // same as bcf
                else error("The preset string not recognised: '%s'\n", optarg);
                break;
            case 's':
                conf.sc = strtol(optarg,&tmp,10);
                if ( *tmp ) error("Could not parse argument: -s %s\n", optarg);
                detect = 0;
                if (conf.preset & TBX_SEQONLY) {
                    conf.bc = conf.ec = conf.sc;
                }
                break;
            case 'S':
                new_line_skip = strtol(optarg,&tmp,10);
                if ( *tmp ) error("Could not parse argument: -S %s\n", optarg);
                detect = 0;
                break;
            case 'D':
                args.download_index = 0;
                break;
            case 1:
                printf(
"tabix (htslib) %s\n"
"Copyright (C) 2025 Genome Research Ltd.\n", hts_version());
                return EXIT_SUCCESS;
            case 2:
                return usage(stdout, EXIT_SUCCESS);
            case 3: {
                int v = atoi(optarg);
                if (v < 0) v = 0;
                hts_set_log_level(v);
                break;
            }
            case 4:
                args.cache_megs = atoi(optarg);
                if (args.cache_megs < 0) {
                    args.cache_megs = 0;
                } else if (args.cache_megs >= INT_MAX / 1048576) {
                    args.cache_megs = INT_MAX / 1048576;
                }
                break;
            case 5:
                args.separate_regs = 1;
                break;
            case 6:
                dump_blocks_flag = 1;
                break;
            case '@':   //thread count
                args.threads = atoi(optarg);
                break;
            default: return usage(stderr, EXIT_FAILURE);
        }
    }

    if (new_line_skip >= 0)
        conf.line_skip = new_line_skip;

    if ( optind==argc ) return usage(stderr, EXIT_FAILURE);

    if ( list_chroms )
        return query_chroms(argv[optind], args.download_index);

    if ( dump_blocks_flag )
        return dump_blocks(argv[optind], args.download_index);

    char *fname = argv[optind];
    int ftype = file_type(fname);
    if ( detect )  // no preset given
    {
        if ( ftype==IS_GFF ) conf = tbx_conf_gff;
        else if ( ftype==IS_BED ) conf = tbx_conf_bed;
        else if ( ftype==IS_GAF ) conf = tbx_conf_gaf;
        else if ( ftype==IS_SAM ) conf = tbx_conf_sam;
        else if ( ftype==IS_VCF )
        {
            conf = tbx_conf_vcf;
            if ( !min_shift && do_csi ) min_shift = 14;
        }
        else if ( ftype==IS_BCF )
        {
            if ( !min_shift ) min_shift = 14;
        }
        else if ( ftype==IS_BAM )
        {
            if ( !min_shift ) min_shift = 14;
        }
    }
    if ( argc > optind+1 || args.header_only || args.regions_fname || args.targets_fname )
    {
        int nregs = 0;
        char **regs = NULL;
        if ( !args.header_only )
            regs = parse_regions(args.regions_fname, argv+optind+1, argc-optind-1, &nregs);
        return query_regions(&args, &conf, fname, regs, nregs);
    }
    if ( do_csi )
    {
        if ( !min_shift ) min_shift = 14;
        min_shift *= do_csi;  // positive for CSIv2, negative for CSIv1
    }
    if ( min_shift!=0 && !do_csi ) do_csi = 1;

    if ( reheader )
        return reheader_file(fname, reheader, ftype, &conf, args.threads);

    char *suffix = ".tbi";
    if ( do_csi ) suffix = ".csi";
    else if ( ftype==IS_BAM ) suffix = ".bai";
    else if ( ftype==IS_CRAM ) suffix = ".crai";

    char *idx_fname = calloc(strlen(fname) + 6, 1);
    if (!idx_fname) error("%s\n", strerror(errno));
    strcat(strcpy(idx_fname, fname), suffix);

    struct stat stat_tbi, stat_file;
    if ( !is_force && stat(idx_fname, &stat_tbi)==0 )
    {
        // Before complaining about existing index, check if the VCF file isn't
        // newer. This is a common source of errors, people tend not to notice
        // that tabix failed
        stat(fname, &stat_file);
        if ( stat_file.st_mtime <= stat_tbi.st_mtime )
            error("[tabix] the index file exists. Please use '-f' to overwrite.\n");
    }
    free(idx_fname);

    int ret;
    if ( ftype==IS_CRAM )
    {
        if ( bam_index_build3(fname, min_shift, args.threads)!=0 ) error("bam_index_build failed: %s\n", fname);
        return 0;
    }
    else if ( do_csi )
    {
        if ( ftype==IS_BCF )
        {
            if ( bcf_index_build3(fname, NULL, min_shift, args.threads)!=0 ) error("bcf_index_build failed: %s\n", fname);
            return 0;
        }
        if ( ftype==IS_BAM )
        {
            if ( bam_index_build3(fname, min_shift, args.threads)!=0 ) error("bam_index_build failed: %s\n", fname);
            return 0;
        }

        switch (ret = tbx_index_build3(fname, NULL, min_shift, args.threads, &conf))
        {
            case 0:
                return 0;
            case -2:
                error("[tabix] the compression of '%s' is not BGZF\n", fname);
            default:
                error("tbx_index_build3 failed: %s\n", fname);
        }
    }
    else    // TBI index
    {
        //CW:INDEXING
        switch (ret = tbx_index_build3(fname, NULL, min_shift, args.threads, &conf))
        {
            case 0:
                return 0;
            case -2:
                error("[tabix] the compression of '%s' is not BGZF\n", fname);
            default:
                error("tbx_index_build3 failed: %s\n", fname);
        }
    }

    return 0;
}
