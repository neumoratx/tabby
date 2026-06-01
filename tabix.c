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
#include <regex.h>
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
#include "htslib/hfile.h"
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

/* ── Filter expression (for -F / --filter) ───────────────────────────────── */
typedef enum {
    FILT_EQ = 0,   /* == */
    FILT_NE,       /* != */
    FILT_GE,       /* >= */
    FILT_LE,       /* <= */
    FILT_RE,       /* ~= (regex match)    */
    FILT_NRE,      /* !~ (regex no-match) */
} FilterOp;

typedef struct {
    int      col;         /* 0-based column index; -1 if col_name not yet resolved */
    char    *col_name;    /* non-NULL if column was specified by name (resolved at query time) */
    FilterOp op;          /* comparison operator                                */
    double   val;         /* numeric threshold (numeric filters only)           */
    char    *sval;        /* string/regex value (non-numeric filters; NULL if numeric) */
    regex_t  compiled_re; /* compiled regex (FILT_RE / FILT_NRE only)          */
    int      has_regex;   /* 1 if compiled_re is valid and must be freed        */
    int      is_or;       /* 1 if from -O (OR group); 0 if from -F (AND group)  */
} filter_expr_t;

typedef struct
{
    char *regions_fname, *targets_fname;
    int print_header, header_only, cache_megs, download_index, separate_regs, threads;
    filter_expr_t *filters;   /* array of parsed -F / -O expressions */
    int            nfilters;  /* total entries in filters[]           */
    int            n_or_filters; /* number of -O (OR-group) entries   */
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

/*
 * sidx_block_index_from_offset() -- map a BGZF compressed block offset to
 * the corresponding block index in the secondary index arrays built by
 * dump_blocks() (i.e. sidx_offsets[] / sidx_blk_offsets[]).
 *
 * Parameters:
 *   blk_off        - compressed file offset of the block to look up.
 *                    This is the upper 48 bits of a BGZF virtual offset,
 *                    i.e.  (uint64_t)itr->curr_off >> 16.
 *   sidx_blk_offsets - array of compressed block offsets, one per block,
 *                    in the order they were recorded by dump_blocks().
 *                    The array MUST be sorted in ascending order (which is
 *                    guaranteed because BGZF blocks are visited sequentially).
 *   n_blocks       - number of entries in sidx_blk_offsets[].
 *
 * Returns:
 *   The zero-based block index in [0, n_blocks) whose compressed offset
 *   equals blk_off, or -1 if blk_off is not found in the array.
 *
 * Complexity: O(log n) via binary search.
 *
 * Usage inside query_regions() at line 385:
 *
 *   uint64_t blk_off = itr->curr_off >> 16;
 *   int64_t  bidx    = sidx_block_index_from_offset(
 *                          blk_off, sidx_blk_offsets, sidx_nblocks);
 *   if (bidx >= 0) {
 *       // use sidx_offsets[bidx] to locate the block record in sidx_records[]
 *       // and read the per-column min/max for secondary filtering
 *   }
 */
static int64_t sidx_block_index_from_offset(uint64_t        blk_off,
                                             const uint64_t *sidx_blk_offsets,
                                             uint64_t        n_blocks)
{
    if (!sidx_blk_offsets || n_blocks == 0)
        return -1;

    int64_t lo = 0, hi = (int64_t)n_blocks - 1;
    while (lo <= hi) {
        int64_t mid = lo + (hi - lo) / 2;
        uint64_t mid_off = sidx_blk_offsets[mid];
        if (mid_off == blk_off)
            return mid;
        else if (mid_off < blk_off)
            lo = mid + 1;
        else
            hi = mid - 1;
    }
    return -1;   /* blk_off not found */
}

/*
 * sidx_block_index_from_offset_bounded() -- same contract as
 * sidx_block_index_from_offset(), but replaces the O(log n) binary search
 * with an O(1) index calculation followed by a short linear scan whose
 * length is bounded by the physics of BGZF block sizes.
 *
 * Algorithm
 * ─────────
 * BGZF blocks are written sequentially, so blk_comp_offsets[] is sorted
 * ascending and the gap between any two consecutive entries is at most
 * max_comp_block_size bytes.  That means:
 *
 *   blk_comp_offsets[i]  >=  blk_off - max_comp_block_size
 *
 * for all i where blk_comp_offsets[i] == blk_off, and therefore the
 * earliest array index that could possibly hold blk_off is:
 *
 *   floor_idx = blk_off / max_comp_block_size
 *               (treating max_comp_block_size as a stride)
 *
 * We clamp floor_idx to [0, n_blocks) and walk forward linearly from
 * there.  Because max_comp_block_size is the *maximum* compressed block
 * size, and a compressed BGZF block is at least 1 byte, the scan visits
 * at most max_comp_block_size + 1 array entries before it must either
 * find blk_off or overshoot it — in practice, typically 1–3 entries.
 *
 * Parameters:
 *   blk_off            - compressed file offset to look up
 *                        (upper 48 bits of a BGZF virtual offset,
 *                         i.e. (uint64_t)itr->curr_off >> 16).
 *   sidx_blk_offsets   - sorted-ascending array of compressed block offsets
 *                        as written by dump_blocks() / sidx_load().
 *   n_blocks           - number of entries in sidx_blk_offsets[].
 *   max_comp_block_size - largest compressed block size (bytes) across the
 *                        file; stored as sidx->max_comp_block_size.
 *                        Must be > 0; if 0 the function falls back to -1.
 *
 * Returns:
 *   Zero-based block index i such that sidx_blk_offsets[i] == blk_off,
 *   or -1 if blk_off is not present in the array.
 *
 * Complexity: O(1) amortised (bounded linear scan of length
 *   ≤ ceil(max_comp_block_size / min_compressed_block_size) + 1).
 */
static int64_t sidx_block_index_from_offset_bounded(
        uint64_t        blk_off,
        const uint64_t *sidx_blk_offsets,
        uint64_t        n_blocks,
        uint16_t        max_comp_block_size)
{
    if (!sidx_blk_offsets || n_blocks == 0 || max_comp_block_size == 0)
        return -1;

    /*
     * Compute the earliest index whose stored offset could equal blk_off.
     *
     * Any block at index i with sidx_blk_offsets[i] == blk_off satisfies:
     *
     *   blk_off == sidx_blk_offsets[i]
     *           >= sidx_blk_offsets[0] + i * 1          (minimum step = 1 byte)
     *
     * and more usefully the *previous* block satisfies:
     *
     *   sidx_blk_offsets[i-1] >= blk_off - max_comp_block_size
     *
     * So blk_off / max_comp_block_size is a lower bound on i (the target
     * index cannot be smaller than this, because every index below it
     * would require a previous block that spans more than
     * max_comp_block_size bytes — which is impossible by construction).
     */
    uint64_t floor_idx = blk_off / (uint64_t)max_comp_block_size;

    /* Clamp to valid range. */
    if (floor_idx >= n_blocks)
        floor_idx = n_blocks - 1;

    /*
     * Linear scan forward from floor_idx.
     *
     * We stop as soon as we find a match or the stored offset exceeds
     * blk_off (the array is sorted, so no later entry can match either).
     */
    for (uint64_t i = floor_idx; i < n_blocks; i++) {
        if (sidx_blk_offsets[i] == blk_off)
            return (int64_t)i;
        if (sidx_blk_offsets[i] > blk_off)
            break;   /* all subsequent entries are strictly larger */
    }

    return -1;   /* blk_off not found */
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

/*
 * tbx_parse_filter() -- parse a single -F/-O argument of the form
 *   <col><op><value>
 * into a filter_expr_t and append it to args->filters[].
 *
 * col    : non-negative integer (0-based column index) OR a column name from
 *          the file's header line (resolved to an index at query time).
 *          Column names must not contain any two-character operator substring.
 * op     : one of  ==  !=  >=  <=  ~=  !~
 * value  : numeric literal, exact string, or ERE regex pattern
 *
 * Returns 0 on success, -1 on parse error (message printed to stderr).
 */
static int tbx_parse_filter(args_t *args, const char *expr, int is_or)
{
    if (!expr || !*expr) {
        fprintf(stderr, "[tabby] -F: empty filter expression\n");
        return -1;
    }

    /* Parse the column spec: a non-negative integer OR a column name.
     * Column names are resolved to 0-based indexes at query time; they must
     * not contain any two-character operator string (==, !=, >=, <=, ~=, !~). */
    char *end;
    long col = strtol(expr, &end, 10);
    char *col_name = NULL;

    if (end == expr) {
        /* Not a numeric index — find the first two-char operator to locate
         * the end of the column name. */
        static const char * const two_char_ops[] =
            { "==", "!=", ">=", "<=", "~=", "!~", NULL };
        const char *op_pos = NULL;
        for (int i = 0; two_char_ops[i]; i++) {
            const char *p = strstr(expr, two_char_ops[i]);
            if (p && (!op_pos || p < op_pos))
                op_pos = p;
        }
        if (!op_pos) {
            fprintf(stderr, "[tabby] -F: no operator found in '%s' "
                            "(expected ==, !=, >=, <=, ~=, or !~)\n", expr);
            return -1;
        }
        size_t name_len = (size_t)(op_pos - expr);
        if (name_len == 0) {
            fprintf(stderr, "[tabby] -F: empty column spec in '%s'\n", expr);
            return -1;
        }
        col_name = malloc(name_len + 1);
        if (!col_name) {
            fprintf(stderr, "[tabby] -F: out of memory\n");
            return -1;
        }
        memcpy(col_name, expr, name_len);
        col_name[name_len] = '\0';
        col = -1;
        end = (char *)op_pos;
    } else if (col < 0) {
        fprintf(stderr, "[tabby] -F: expected non-negative column index at "
                        "start of '%s'\n", expr);
        return -1;
    }

    /* Parse the two-character operator. */
    FilterOp op;
    if      (end[0] == '=' && end[1] == '=') { op = FILT_EQ;  end += 2; }
    else if (end[0] == '!' && end[1] == '=') { op = FILT_NE;  end += 2; }
    else if (end[0] == '>' && end[1] == '=') { op = FILT_GE;  end += 2; }
    else if (end[0] == '<' && end[1] == '=') { op = FILT_LE;  end += 2; }
    else if (end[0] == '~' && end[1] == '=') { op = FILT_RE;  end += 2; }
    else if (end[0] == '!' && end[1] == '~') { op = FILT_NRE; end += 2; }
    else {
        free(col_name);
        fprintf(stderr, "[tabby] -F: unrecognised operator in '%s' "
                        "(expected ==, !=, >=, <=, ~=, or !~)\n", expr);
        return -1;
    }

    /* Parse the value: regex operators always take a string; others try numeric
     * first and fall back to string. */
    if (!*end) {
        free(col_name);
        fprintf(stderr, "[tabby] -F: missing value after operator in '%s'\n",
                expr);
        return -1;
    }

    int is_string = 0;
    double val = 0.0;

    if (op == FILT_RE || op == FILT_NRE) {
        is_string = 1;   /* regex operators always treat value as a pattern */
    } else {
        char *vend;
        errno = 0;
        val = strtod(end, &vend);
        is_string = (errno != 0 || vend == end || *vend != '\0');

        if (is_string && op != FILT_EQ && op != FILT_NE) {
            free(col_name);
            fprintf(stderr, "[tabby] -F: string filters only support ==, !=, ~=, and !~; "
                            "got '%s'\n", expr);
            return -1;
        }
    }

    /* Append to args->filters[]. */
    filter_expr_t *tmp = realloc(args->filters,
                                 (size_t)(args->nfilters + 1)
                                 * sizeof(filter_expr_t));
    if (!tmp) {
        free(col_name);
        fprintf(stderr, "[tabby] -F: out of memory\n");
        return -1;
    }
    args->filters = tmp;
    filter_expr_t *f = &args->filters[args->nfilters];
    f->col       = (int)col;
    f->col_name  = col_name;
    f->op        = op;
    f->val       = is_string ? 0.0 : val;
    f->sval      = NULL;
    f->has_regex = 0;
    f->is_or     = is_or;

    if (is_string) {
        f->sval = strdup(end);
        if (!f->sval) {
            free(f->col_name); f->col_name = NULL;
            fprintf(stderr, "[tabby] -F: out of memory\n");
            return -1;
        }
        if (op == FILT_RE || op == FILT_NRE) {
            int rc = regcomp(&f->compiled_re, end, REG_EXTENDED | REG_NOSUB);
            if (rc != 0) {
                char errbuf[256];
                regerror(rc, &f->compiled_re, errbuf, sizeof(errbuf));
                fprintf(stderr, "[tabby] -F: invalid regex '%s': %s\n",
                        end, errbuf);
                free(f->sval);
                f->sval = NULL;
                free(f->col_name); f->col_name = NULL;
                return -1;
            }
            f->has_regex = 1;
        }
    }
    args->nfilters++;
    if (is_or) args->n_or_filters++;
    return 0;
}

/*
 * tbx_apply_filters() -- test one decoded data line against all -F filter
 * expressions (per-field exact check only).
 *
 * Block/chunk-level pre-filtering is handled upstream by tbx_filter_chunks(),
 * which prunes iter->off[] before any I/O occurs.  By the time this function
 * is called the line has already been decompressed and decoded; we simply
 * tokenise it and compare each target field value against the operator+value.
 *
 * Parameters:
 *   args   - contains the filter list
 *   line   - raw tab-delimited text line (NUL-terminated)
 *   len    - strlen(line)
 *
 * Returns 1 if the line passes ALL filters, 0 if rejected by any filter.
 *
 * If the referenced column is absent from the line the line is kept
 * (conservative: no data != failed filter).  If the column value is
 * non-numeric an error is printed and the line is dropped.
 */
static int tbx_apply_filters(const args_t *args,
                              const char   *line,
                              size_t        len)
{
    if (!args->nfilters) return 1;   /* no filters → accept everything */

    /* Two-pass evaluation:
     *   Pass 1: AND filters (-F) — all must pass; fail fast on first failure.
     *   Pass 2: OR  filters (-O) — at least one must pass (skipped when none). */
    int or_passed = 0;

    for (int fi = 0; fi < args->nfilters; fi++) {
        const filter_expr_t *f = &args->filters[fi];

        /* Safety net: unresolved column name (resolve_filter_col_names should
         * have been called before we reach here; skip conservatively if not). */
        if (f->col < 0) continue;

        /* Chunk-level block pruning is handled upstream by tbx_filter_chunks().
         * Here we only perform the exact per-field check on the decoded line. */

        /* ── Exact per-field check: tokenise line to reach column f->col ── */
        const char *p    = line;
        const char *lend = line + len;
        int ci = 0;

        /* Advance past columns 0 .. f->col-1. */
        while (ci < f->col) {
            const char *tab = memchr(p, '\t', (size_t)(lend - p));
            if (!tab) break;   /* column not present in this line */
            p = tab + 1;
            ci++;
        }

        if (ci < f->col) {
            /* Column doesn't exist in this line; keep conservatively. */
            continue;
        }

        /* Find the end of this field. */
        const char *tab       = memchr(p, '\t', (size_t)(lend - p));
        const char *field_end = tab ? tab : lend;
        size_t      flen      = (size_t)(field_end - p);

        /* Apply the filter: string or numeric path. */
        int pass;
        if (f->sval) {
            if (f->has_regex) {
                /* Regex filter: copy field to a NUL-terminated buffer for regexec. */
                char *rbuf = NULL;
                char rbuf_stack[256];
                if (flen + 1 <= sizeof(rbuf_stack)) {
                    rbuf = rbuf_stack;
                } else {
                    rbuf = malloc(flen + 1);
                    if (!rbuf) {
                        fprintf(stderr,
                                "[tabby] -F: out of memory for regex match\n");
                        return 0;
                    }
                }
                memcpy(rbuf, p, flen);
                rbuf[flen] = '\0';
                int matched = (regexec(&f->compiled_re, rbuf, 0, NULL, 0) == 0);
                if (rbuf != rbuf_stack) free(rbuf);
                switch (f->op) {
                    case FILT_RE:  pass = matched;  break;
                    case FILT_NRE: pass = !matched; break;
                    default:       pass = 0;        break;
                }
            } else {
                /* Exact string filter: compare field bytes directly. */
                size_t slen = strlen(f->sval);
                int eq = (flen == slen && memcmp(p, f->sval, flen) == 0);
                switch (f->op) {
                    case FILT_EQ: pass = eq;  break;
                    case FILT_NE: pass = !eq; break;
                    default:      pass = 0;   break;
                }
            }
        } else {
            /* Numeric filter. */
            char field_buf[64];
            if (flen == 0 || flen >= sizeof(field_buf)) {
                fprintf(stderr,
                        "[tabby] -F: column %d value is empty or too long for "
                        "a numeric filter\n", f->col);
                return 0;
            }
            memcpy(field_buf, p, flen);
            field_buf[flen] = '\0';

            char *nend;
            errno = 0;
            double field_val = strtod(field_buf, &nend);
            if (errno || nend == field_buf || *nend != '\0') {
                fprintf(stderr,
                        "[tabby] -F: column %d value '%s' is non-numeric; "
                        "use == / != for exact string match or ~= / !~ for regex\n",
                        f->col, field_buf);
                return 0;
            }

            switch (f->op) {
                case FILT_EQ: pass = (field_val == f->val); break;
                case FILT_NE: pass = (field_val != f->val); break;
                case FILT_GE: pass = (field_val >= f->val); break;
                case FILT_LE: pass = (field_val <= f->val); break;
                default:      pass = 1; break;
            }
        }
        if (f->is_or) {
            if (pass) or_passed = 1;   /* at least one OR filter satisfied */
        } else {
            if (!pass) return 0;       /* AND filter failed → reject immediately */
        }
    }

    /* If any OR filters were registered, at least one must have passed. */
    if (args->n_or_filters > 0 && !or_passed) return 0;

    return 1;   /* all AND filters passed; OR condition satisfied */
}


/*
 * tbx_filter_chunks() -- prune the iterator's chunk list (iter->off[]) using
 * the secondary index BEFORE any BGZF blocks are seeked to or downloaded.
 *
 * Background
 * ──────────
 * After tbx_itr_querys() / hts_itr_query() returns, iter->off[0..n_off-1] is
 * a sorted, merged array of hts_pair64_max_t {u, v, max}.  Each entry is a
 * *chunk* — a contiguous byte range [u, v) in the compressed file that the
 * TBI primary index believes may contain records for the queried region.
 * hts_itr_next() will later seek to each chunk's start (u) and read forward
 * to v, decompressing every BGZF block in that range.
 *
 * For each chunk the BGZF block(s) it covers have compressed offsets in the
 * half-open interval  [ u>>16 , v>>16 ]  (upper 48 bits of the virtual
 * offset).  Using the secondary index we can inspect the per-column min/max
 * of every block in that range and, when *every* block in the chunk is
 * guaranteed to fail a filter, remove the chunk entirely from off[] before
 * any I/O happens.
 *
 * Algorithm
 * ─────────
 * For each chunk [u, v):
 *   1. Find the first sidx block whose compressed offset >= u>>16  (lower
 *      bound via binary search on blk_comp_offsets[]).
 *   2. Walk forward through sidx blocks while blk_comp_offsets[bi] <= v>>16.
 *   3. For each filter expression f, check whether ALL blocks in the chunk
 *      are impossible for f (block_max < f->val for >=, etc.).
 *      If so, the whole chunk cannot contain any row that passes f → drop it.
 *   4. If no filter can rule out the chunk, keep it.
 *
 * A chunk that spans zero known sidx blocks (e.g. it lies entirely between
 * recorded blocks) is kept conservatively.
 *
 * The compaction is done in-place: kept entries are shifted left, then
 * iter->n_off is updated.  No realloc is needed.
 *
 * Parameters:
 *   args  - contains the filter expression list
 *   sidx  - loaded secondary index (function is a no-op if NULL or no blocks)
 *   itr   - the iterator whose off[] will be pruned in-place
 */
static void tbx_filter_chunks(const args_t *args,
                               const sidx_t *sidx,
                               hts_itr_t    *itr)
{
    /* Nothing to do if there are no filters, no secondary index, or no chunks. */
    if (!args->nfilters || !sidx || !sidx->nblocks || !itr || itr->n_off == 0)
        return;

    int kept = 0;   /* write cursor for the in-place compaction */

    for (int ci = 0; ci < itr->n_off; ci++) {
        uint64_t chunk_u_blk = itr->off[ci].u >> 16;   /* first compressed block offset */
        uint64_t chunk_v_blk = itr->off[ci].v >> 16;   /* last  compressed block offset  */

        /* Binary search: first sidx block with offset >= chunk_u_blk. */
        int64_t lo = 0, hi = (int64_t)sidx->nblocks - 1, first_bi = -1;
        while (lo <= hi) {
            int64_t mid = lo + (hi - lo) / 2;
            if (sidx->blk_comp_offsets[mid] >= chunk_u_blk) {
                first_bi = mid;
                hi = mid - 1;
            } else {
                lo = mid + 1;
            }
        }

        /* If no sidx block starts at or after chunk_u_blk, keep conservatively. */
        if (first_bi < 0) {
            itr->off[kept++] = itr->off[ci];
            continue;
        }

        /* Count how many sidx blocks fall within this chunk. */
        int64_t last_bi = first_bi;
        while (last_bi + 1 < (int64_t)sidx->nblocks &&
               sidx->blk_comp_offsets[last_bi + 1] <= chunk_v_blk)
            last_bi++;

        /* If no sidx block falls within [chunk_u_blk, chunk_v_blk], keep. */
        if (sidx->blk_comp_offsets[first_bi] > chunk_v_blk) {
            itr->off[kept++] = itr->off[ci];
            continue;
        }

        /*
         * For each filter, check whether the filter is impossible for ALL
         * sidx blocks in [first_bi, last_bi].  If so, the chunk cannot
         * contain any matching row and can be dropped.
         *
         * "impossible for a block" uses the same block-level rules as
         * tbx_apply_filters():
         *   ==  impossible when val < blk_min  OR  val > blk_max
         *   !=  impossible when blk_min == blk_max == val
         *   >=  impossible when blk_max < val
         *   <=  impossible when blk_min > val
         *
         * We need ALL blocks to be impossible for the filter to drop the chunk.
         */
        int chunk_dropped = 0;
        for (int fi = 0; fi < args->nfilters && !chunk_dropped; fi++) {
            const filter_expr_t *f = &args->filters[fi];

            /* OR filters and string/regex filters cannot safely prune blocks — skip. */
            if (f->is_or || f->sval) continue;

            /* Only block-filter numeric columns that are within ncols. */
            if ((uint32_t)f->col >= sidx->ncols) continue;

            uint8_t tc_min = sidx->col_min_tc[f->col];
            uint8_t tc_max = sidx->col_max_tc[f->col];
            if (tc_min == 0 || tc_max == 0) continue;   /* string/NA column */

            /* Check every block in [first_bi, last_bi].
             * If even one block *might* satisfy the filter we stop early. */
            int all_blocks_impossible = 1;
            for (int64_t bi = first_bi; bi <= last_bi; bi++) {
                /* Walk the block record to column f->col. */
                const uint8_t *rec = sidx->records
                                   + sidx->blk_rec_offsets[bi];
                for (int col = 0; col < f->col; col++) {
                    uint8_t tcm = sidx->col_min_tc[col];
                    uint8_t tcM = sidx->col_max_tc[col];
                    if (tcm == 0 || tcM == 0) rec += 1;
                    else { rec += sidx_type_size(tcm); rec += sidx_type_size(tcM); }
                }

                /* Read blk_min and blk_max as double. */
                double blk_min, blk_max;
#define TBX_CHUNK_READ(tc, ptr, dst) do { switch (tc) { \
    case 1:{uint8_t  _v;memcpy(&_v,(ptr),1);(dst)=(double)_v;break;} \
    case 2:{int16_t  _v;memcpy(&_v,(ptr),2);(dst)=(double)_v;break;} \
    case 3:{uint16_t _v;memcpy(&_v,(ptr),2);(dst)=(double)_v;break;} \
    case 4:{int32_t  _v;memcpy(&_v,(ptr),4);(dst)=(double)_v;break;} \
    case 5:{uint32_t _v;memcpy(&_v,(ptr),4);(dst)=(double)_v;break;} \
    case 6:{int64_t  _v;memcpy(&_v,(ptr),8);(dst)=(double)_v;break;} \
    case 7:{uint64_t _v;memcpy(&_v,(ptr),8);(dst)=(double)_v;break;} \
    case 8:{float    _v;memcpy(&_v,(ptr),4);(dst)=(double)_v;break;} \
    case 9:{double   _v;memcpy(&_v,(ptr),8);(dst)=_v;        break;} \
    default:(dst)=0;break; \
}} while(0)
                TBX_CHUNK_READ(tc_min, rec, blk_min);
                rec += sidx_type_size(tc_min);
                TBX_CHUNK_READ(tc_max, rec, blk_max);
#undef TBX_CHUNK_READ

                int blk_impossible;
                switch (f->op) {
                    case FILT_EQ:
                        blk_impossible = (f->val < blk_min || f->val > blk_max); break;
                    case FILT_NE:
                        blk_impossible = (blk_min == f->val && blk_max == f->val); break;
                    case FILT_GE:
                        blk_impossible = (blk_max < f->val); break;
                    case FILT_LE:
                        blk_impossible = (blk_min > f->val); break;
                    default:
                        blk_impossible = 0; break;
                }

                if (!blk_impossible) {
                    all_blocks_impossible = 0;
                    break;   /* at least one block might satisfy f → keep chunk */
                }
            }

            if (all_blocks_impossible)
                chunk_dropped = 1;   /* every block in chunk fails this filter */
        }

        if (!chunk_dropped)
            itr->off[kept++] = itr->off[ci];
    }

    itr->n_off = kept;
    if (kept == 0)
        itr->finished = 1;   /* no chunks left → iterator is already done */
}

/*
 * State passed to sidx_block_filter_cb() and bgzf_decomp_skip_cb() via the
 * respective callback data pointers.  Stack-allocated once per region query.
 */
typedef struct {
    const sidx_t        *sidx;
    const filter_expr_t *filters;
    int                  nfilters;
    const hts_itr_t     *itr;    /* read-only; used by bgzf_decomp_skip_cb to
                                  * read itr->off[itr->i].v (the current chunk
                                  * end) and avoid skipping past it          */
} tbx_block_filter_state_t;

/*
 * sidx_block_filter_cb() -- hts_block_filter_func implementation for tabix.
 *
 * Called by hts_itr_next() at every BGZF block boundary (single-threaded,
 * non-gzip path only) with the compressed offset of the block about to be
 * loaded.  Uses the secondary index min/max statistics to decide whether
 * any record in that block can possibly satisfy the -F filters.
 *
 * Return value follows the hts_block_filter_func contract:
 *   UINT64_MAX  — block passes; read it normally.
 *   other       — block fails; returned value is the compressed offset of
 *                 the next sidx-recorded block to try (or UINT64_MAX-1 as
 *                 a sentinel meaning "no further blocks known, skip to
 *                 chunk end").
 *
 * Internally this is a two-step lookup:
 *   1. Binary-search sidx->blk_comp_offsets[] for block_address.
 *      If not found, the block was not recorded in the .sidx — accept it
 *      conservatively (we cannot rule it out).
 *   2. For each filter, apply the block-level impossibility test against
 *      the block's stored min/max.  If ANY filter is impossible for this
 *      block, the block is rejected and we return the next block's address.
 */
static uint64_t sidx_block_filter_cb(uint64_t block_address, void *data)
{
    const tbx_block_filter_state_t *st = (const tbx_block_filter_state_t *)data;
    const sidx_t *sidx = st->sidx;

    /* Binary-search for this block in blk_comp_offsets[]. */
    int64_t lo = 0, hi = (int64_t)sidx->nblocks - 1, bidx = -1;
    while (lo <= hi) {
        int64_t mid = lo + (hi - lo) / 2;
        uint64_t mid_off = sidx->blk_comp_offsets[mid];
        if (mid_off == block_address) { bidx = mid; break; }
        else if (mid_off < block_address) lo = mid + 1;
        else hi = mid - 1;
    }

    if (bidx < 0) {
        /* Block not in sidx — cannot prune, accept conservatively. */
        return UINT64_MAX;
    }

    /* Test each filter against this block's min/max. */
    for (int fi = 0; fi < st->nfilters; fi++) {
        const filter_expr_t *f = &st->filters[fi];

        if (f->is_or || f->sval) continue;  /* OR / string / regex: skip block pruning */
        if ((uint32_t)f->col >= sidx->ncols) continue;
        uint8_t tc_min = sidx->col_min_tc[f->col];
        uint8_t tc_max = sidx->col_max_tc[f->col];
        if (tc_min == 0 || tc_max == 0) continue;  /* string/NA column */

        /* Walk the block record to column f->col. */
        const uint8_t *rec = sidx->records + sidx->blk_rec_offsets[bidx];
        for (int ci = 0; ci < f->col; ci++) {
            uint8_t tcm = sidx->col_min_tc[ci];
            uint8_t tcM = sidx->col_max_tc[ci];
            if (tcm == 0 || tcM == 0) rec += 1;
            else { rec += sidx_type_size(tcm); rec += sidx_type_size(tcM); }
        }

        /* Read blk_min and blk_max as double. */
        double blk_min, blk_max;
#define SIDX_CB_READ(tc, ptr, dst) do { switch (tc) { \
    case 1:{uint8_t  _v;memcpy(&_v,(ptr),1);(dst)=(double)_v;break;} \
    case 2:{int16_t  _v;memcpy(&_v,(ptr),2);(dst)=(double)_v;break;} \
    case 3:{uint16_t _v;memcpy(&_v,(ptr),2);(dst)=(double)_v;break;} \
    case 4:{int32_t  _v;memcpy(&_v,(ptr),4);(dst)=(double)_v;break;} \
    case 5:{uint32_t _v;memcpy(&_v,(ptr),4);(dst)=(double)_v;break;} \
    case 6:{int64_t  _v;memcpy(&_v,(ptr),8);(dst)=(double)_v;break;} \
    case 7:{uint64_t _v;memcpy(&_v,(ptr),8);(dst)=(double)_v;break;} \
    case 8:{float    _v;memcpy(&_v,(ptr),4);(dst)=(double)_v;break;} \
    case 9:{double   _v;memcpy(&_v,(ptr),8);(dst)=_v;        break;} \
    default:(dst)=0;break; \
}} while(0)
        SIDX_CB_READ(tc_min, rec, blk_min);
        rec += sidx_type_size(tc_min);
        SIDX_CB_READ(tc_max, rec, blk_max);
#undef SIDX_CB_READ

        int impossible;
        switch (f->op) {
            case FILT_EQ: impossible = (f->val < blk_min || f->val > blk_max); break;
            case FILT_NE: impossible = (blk_min == f->val && blk_max == f->val); break;
            case FILT_GE: impossible = (blk_max < f->val); break;
            case FILT_LE: impossible = (blk_min > f->val); break;
            default:      impossible = 0; break;
        }

        if (impossible) {
            /*
             * Block rejected.  Return the next sidx-recorded block address
             * so hts_itr_next() can seek directly there without scanning.
             * If this was the last recorded block, return UINT64_MAX-1 as
             * a sentinel: hts_itr_next() will treat any value <= cur_blk or
             * > chunk_end as "skip to chunk end".
             */
            if (bidx + 1 < (int64_t)sidx->nblocks)
                return sidx->blk_comp_offsets[bidx + 1];
            else
                return UINT64_MAX - 1;  /* no next block known */
        }
    }

    return UINT64_MAX;  /* all filters pass for this block */
}

/*
 * bgzf_decomp_skip_cb() -- BGZF::decomp_skip_func implementation for tabix.
 *
 * This is the post-download, pre-decompression counterpart to
 * sidx_block_filter_cb().  It is installed on bgzfp->decomp_skip_func for
 * remote files, where bgzf_seek is expensive (HTTP round-trip) and therefore
 * intra-chunk seeking at the hts_itr_t level is disabled.  Instead, this
 * callback fires inside bgzf_read_block() after the compressed bytes are
 * already in memory but before inflate_block() runs, saving the decompression
 * CPU cost for blocks that cannot contain matching rows.
 *
 * Uses the same tbx_block_filter_state_t and sidx lookup as
 * sidx_block_filter_cb(), but has the signature required by decomp_skip_func:
 *   int(int64_t block_address, void *data)
 * returning 0 = decompress normally, non-zero = skip this block.
 *
 * Note: block_address is int64_t here (matching the BGZF struct field type)
 * vs uint64_t in sidx_block_filter_cb.  BGZF compressed offsets are always
 * non-negative so the cast is safe.
 */
static int bgzf_decomp_skip_cb(int64_t block_address, void *data)
{
    const tbx_block_filter_state_t *st = (const tbx_block_filter_state_t *)data;
    const sidx_t *sidx = st->sidx;

    /*
     * Chunk-boundary guard.
     *
     * bgzf_read_block() has no knowledge of iterator chunk boundaries.
     * If this block's compressed offset is at or beyond the current chunk
     * end (iter->off[iter->i].v >> 16), we must NOT skip it here — doing
     * so would cause bgzf_read_block() to loop and load the next block,
     * which lies outside the chunk hts_itr_next() is currently processing.
     * hts_itr_next() would then see curr_off past the chunk end on the
     * next bgzf_tell() read, but only after a spurious decompression of
     * a block we were never supposed to touch.
     *
     * Instead we return 0 (decompress normally) and let hts_itr_next()
     * handle the chunk transition in its own logic.  Note that itr->i may
     * be -1 before the first chunk is entered; guard against that too.
     */
    if (st->itr && st->itr->i >= 0) {
        uint64_t chunk_end_blk = st->itr->off[st->itr->i].v >> 16;
        if ((uint64_t)block_address >= chunk_end_blk)
            return 0;   /* at or past chunk end — do not skip */
    }

    /* Binary-search for this block in blk_comp_offsets[]. */
    int64_t lo = 0, hi = (int64_t)sidx->nblocks - 1, bidx = -1;
    while (lo <= hi) {
        int64_t mid = lo + (hi - lo) / 2;
        uint64_t mid_off = sidx->blk_comp_offsets[mid];
        if (mid_off == (uint64_t)block_address) { bidx = mid; break; }
        else if (mid_off < (uint64_t)block_address) lo = mid + 1;
        else hi = mid - 1;
    }

    if (bidx < 0)
        return 0;   /* block not in sidx — decompress conservatively */

    /* Test each filter against this block's min/max. */
    for (int fi = 0; fi < st->nfilters; fi++) {
        const filter_expr_t *f = &st->filters[fi];

        if (f->is_or || f->sval) continue;  /* OR / string / regex: skip block pruning */
        if ((uint32_t)f->col >= sidx->ncols) continue;
        uint8_t tc_min = sidx->col_min_tc[f->col];
        uint8_t tc_max = sidx->col_max_tc[f->col];
        if (tc_min == 0 || tc_max == 0) continue;

        const uint8_t *rec = sidx->records + sidx->blk_rec_offsets[bidx];
        for (int ci = 0; ci < f->col; ci++) {
            uint8_t tcm = sidx->col_min_tc[ci];
            uint8_t tcM = sidx->col_max_tc[ci];
            if (tcm == 0 || tcM == 0) rec += 1;
            else { rec += sidx_type_size(tcm); rec += sidx_type_size(tcM); }
        }

        double blk_min, blk_max;
#define BGZF_CB_READ(tc, ptr, dst) do { switch (tc) { \
    case 1:{uint8_t  _v;memcpy(&_v,(ptr),1);(dst)=(double)_v;break;} \
    case 2:{int16_t  _v;memcpy(&_v,(ptr),2);(dst)=(double)_v;break;} \
    case 3:{uint16_t _v;memcpy(&_v,(ptr),2);(dst)=(double)_v;break;} \
    case 4:{int32_t  _v;memcpy(&_v,(ptr),4);(dst)=(double)_v;break;} \
    case 5:{uint32_t _v;memcpy(&_v,(ptr),4);(dst)=(double)_v;break;} \
    case 6:{int64_t  _v;memcpy(&_v,(ptr),8);(dst)=(double)_v;break;} \
    case 7:{uint64_t _v;memcpy(&_v,(ptr),8);(dst)=(double)_v;break;} \
    case 8:{float    _v;memcpy(&_v,(ptr),4);(dst)=(double)_v;break;} \
    case 9:{double   _v;memcpy(&_v,(ptr),8);(dst)=_v;        break;} \
    default:(dst)=0;break; \
}} while(0)
        BGZF_CB_READ(tc_min, rec, blk_min);
        rec += sidx_type_size(tc_min);
        BGZF_CB_READ(tc_max, rec, blk_max);
#undef BGZF_CB_READ

        int impossible;
        switch (f->op) {
            case FILT_EQ: impossible = (f->val < blk_min || f->val > blk_max); break;
            case FILT_NE: impossible = (blk_min == f->val && blk_max == f->val); break;
            case FILT_GE: impossible = (blk_max < f->val); break;
            case FILT_LE: impossible = (blk_min > f->val); break;
            default:      impossible = 0; break;
        }
        if (impossible) return 1;   /* skip: don't decompress */
    }
    return 0;   /* all filters pass — decompress normally */
}

/*
 * resolve_filter_col_names() -- map any column-name filter specs to their
 * 0-based column indexes by reading the first line of the file as a header.
 *
 * The header line is expected to be the first line of the BGZF file.
 * A leading meta_char (e.g. '#') is stripped before splitting by tab.
 * This handles both comment-style headers (#COL1\tCOL2\t...) and plain
 * headers used with -S 1 (line_skip=1).
 *
 * Returns 0 on success (all names resolved or no named filters present),
 * -1 on error (message printed to stderr).
 */
static int resolve_filter_col_names(args_t *args, const char *fname,
                                    const tbx_conf_t *conf)
{
    int has_names = 0;
    for (int i = 0; i < args->nfilters; i++)
        if (args->filters[i].col_name) { has_names = 1; break; }
    if (!has_names) return 0;

    BGZF *fp = bgzf_open(fname, "r");
    if (!fp) {
        fprintf(stderr, "[tabby] cannot open '%s' to resolve column names: %s\n",
                fname, strerror(errno));
        return -1;
    }

    kstring_t line = {0, 0, NULL};
    int ret = bgzf_getline(fp, '\n', &line);
    bgzf_close(fp);

    if (ret < 0 || !line.l) {
        fprintf(stderr, "[tabby] cannot read header line from '%s'\n", fname);
        free(line.s);
        return -1;
    }

    /* Strip a leading meta_char (e.g. '#') if present. */
    char *hdr = line.s;
    if (conf->meta_char && (unsigned char)*hdr == (unsigned char)conf->meta_char)
        hdr++;

    int rc = 0;
    for (int fi = 0; fi < args->nfilters && rc == 0; fi++) {
        filter_expr_t *f = &args->filters[fi];
        if (!f->col_name) continue;

        char *p = hdr;
        int ci = 0, found = 0;
        while (*p && !found) {
            char *tab = strchr(p, '\t');
            size_t flen = tab ? (size_t)(tab - p) : strlen(p);
            if (flen == strlen(f->col_name) &&
                    memcmp(p, f->col_name, flen) == 0) {
                f->col = ci;
                found = 1;
            }
            if (!found) {
                if (!tab) break;
                p = tab + 1;
                ci++;
            }
        }
        if (!found) {
            fprintf(stderr, "[tabby] -F: column name '%s' not found in "
                            "header of '%s'\n", f->col_name, fname);
            rc = -1;
        }
    }

    free(line.s);
    return rc;
}

static int query_regions(args_t *args, tbx_conf_t *conf, char *fname, char **regs, int nregs,const sidx_t *sidx)
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

        /* Resolve any column names in -F/-O expressions to 0-based indexes
         * using the header line from the file. */
        if (resolve_filter_col_names(args, fname, &tbx->conf) != 0) {
            free(str.s);
            tbx_destroy(tbx);
            if (reg_idx) regidx_destroy(reg_idx);
            hts_close(fp);
            RELEASE_TPOOL(tpool.pool);
            return 1;
        }

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
            /*
             * Evaluate once before the per-region loop: both answers are
             * constant for the lifetime of this query call.
             *
             * file_is_remote — hisremote() is a pure scheme-prefix table
             *   lookup (no I/O), but calling it per-region would be
             *   needlessly repetitive.  Remote files must not use the
             *   intra-chunk block-seek optimisation because bgzf_seek on a
             *   remote hFILE triggers a new HTTP range request, which would
             *   be more expensive than reading the block.
             *
             * bgzfp — hts_get_bgzfp() is also O(1) and returns the same
             *   pointer every call; cache it here for the same reason.
             */
            int file_is_remote = hisremote(fname);
            BGZF *bgzfp = hts_get_bgzfp(fp);
            for (i=0; i<nregs; i++)
            {
                int ret, found = 0;
                char *region = regs[i];
                kstring_t expanded = {0, 0, NULL};
                if (is_seqonly && strchr(region, ':') == NULL && strcmp(region, ".") != 0) {
                    if (ksprintf(&expanded, "%s:1-1", region) < 0) {
                        RELEASE_TPOOL(tpool.pool);
                        error_errno(NULL);
                    }
                    region = expanded.s;
                }
                hts_itr_t *itr = tbx_itr_querys(tbx, region);
                free(expanded.s);
                if ( !itr ) continue;
                /* Phase 1 (chunk-level): prune entire BGZF chunks from the
                 * iterator's offset list using the secondary index min/max,
                 * before any blocks are seeked to or downloaded. */
                tbx_filter_chunks(args, sidx, itr);

                /*
                 * Phase 1b (intra-chunk block-level): attach a per-block
                 * filter callback to the iterator so hts_itr_next() can
                 * seek past individual BGZF blocks within kept chunks whose
                 * sidx min/max statistics guarantee no record can match.
                 *
                 * Activated only when:
                 *  • a .sidx is loaded and has at least one block
                 *  • at least one -F filter is active
                 *  • the BGZF stream is single-threaded and not gzip
                 *  • the file is local (not remote)
                 *
                 * The remote exclusion is critical: bgzf_seek on a remote
                 * hFILE issues a new HTTP range request per call.  Skipping
                 * individual blocks that way would create more round trips
                 * than reading them, making things strictly worse.  For
                 * remote files, chunk-level pruning (Phase 1) is the right
                 * granularity — one range request per kept chunk.
                 *
                 * The state struct is stack-allocated; its lifetime covers
                 * the tbx_itr_next loop below and nothing else.
                 */
                tbx_block_filter_state_t bfs = {0};
                if (sidx && sidx->nblocks > 0
                        && args->nfilters > 0
                        && bgzfp && !bgzfp->mt && !bgzfp->is_gzip) {
                    bfs.sidx     = sidx;
                    bfs.filters  = args->filters;
                    bfs.nfilters = args->nfilters;
                    bfs.itr      = itr;
                    if (!file_is_remote) {
                        /* Local file: bgzf_seek is a cheap OS seek.
                         * Use the hts_itr_t callback to skip entire blocks
                         * before they are even read from disk. */
                        itr->block_filter_func = sidx_block_filter_cb;
                        itr->block_filter_data = &bfs;
                    } else {
                        /* Remote file: bgzf_seek costs an HTTP round-trip,
                         * so we cannot skip blocks before download.
                         * Install the BGZF-level callback instead: it fires
                         * after the compressed bytes are in memory but before
                         * inflate_block() runs, saving the decompression CPU
                         * cost for blocks that cannot contain matching rows.
                         * Chunk-level pruning (Phase 1) already minimised
                         * the number of blocks downloaded. */
                        bgzfp->decomp_skip_func = bgzf_decomp_skip_cb;
                        bgzfp->decomp_skip_data = &bfs;
                    }
                }

                while ((ret = tbx_itr_next(fp, tbx, itr, &str)) >= 0)
                {
                    if ( reg_idx && !regidx_overlap(reg_idx,seq[itr->curr_tid],itr->curr_beg,itr->curr_end-1, NULL) ) continue;

                    /* Phase 2 (per-field exact check): chunk-level pruning was
                     * already done by tbx_filter_chunks() before any I/O. */
                    if (!tbx_apply_filters(args, str.s, str.l))
                        continue;

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
                /* Clear any BGZF-level decomp skip callback installed for
                 * this region: bfs is stack-allocated and goes out of scope
                 * at the end of this loop iteration. */
                if (bgzfp) {
                    bgzfp->decomp_skip_func = NULL;
                    bgzfp->decomp_skip_data = NULL;
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
 *  Field 0  magic[4]                      char[4]      "TBA\x02"
 *  Field 1  n_columns                     uint32_t     # of tab-separated columns
 *  Field 2a column_min_types[ncols]       char[ncols]  type code for each col's min
 *  Field 2b column_max_types[ncols]       char[ncols]  type code for each col's max
 *  Field 3  n_blocks                      uint64_t     # of BGZF blocks
 *  Field 3b max_comp_block_size           uint16_t     largest compressed block size (bytes)
 *                                                      across the entire .tbi index file
 *  Field 4  block_record_length           uint64_t     byte length of one block record
 *  Field 5  offsets[n_blocks]             uint64_t[]   byte offset of each block record
 *                                                      within the block_records section
 *  Field 6  block_records                 char[]       n_blocks * block_record_length bytes
 *  Field 7  compressed_offsets[n_blocks]  uint64_t[]   BGZF compressed file offset of each
 *                                                      block (upper 48 bits of virtual offset);
 *                                                      sorted ascending; used to map
 *                                                      itr->curr_off >> 16 → block index via
 *                                                      sidx_block_index_from_offset()
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

/* ── Main dump_blocks function ────────────────────────────────────────────── */



/*
 * sidx_free() -- release all memory owned by *sidx and zero the struct.
 * Safe to call on a zero-initialised sidx_t (all pointers are NULL).
 */
static void sidx_free(sidx_t *sidx)
{
    if (!sidx) return;
    free(sidx->col_min_tc);
    free(sidx->col_max_tc);
    free(sidx->blk_rec_offsets);
    free(sidx->records);
    free(sidx->blk_comp_offsets);
    memset(sidx, 0, sizeof(*sidx));
}

/*
 * sidx_load() -- load a .sidx file written by dump_blocks() into *out.
 *
 * Parameters:
 *   sidx_fname  - path to the .sidx file (e.g. "myfile.vcf.gz.sidx")
 *   out         - caller-allocated sidx_t to populate; must be zero-initialised
 *                 on entry (e.g. `sidx_t sidx; memset(&sidx,0,sizeof(sidx));`).
 *                 On success *out owns all allocated memory; call sidx_free(out)
 *                 when done.  On failure *out is left in a consistent, freeable
 *                 state (sidx_free() is safe to call).
 *
 * Returns:
 *    0  on success
 *   -1  on I/O error or format mismatch (a message is printed to stderr)
 *
 * The function validates:
 *   • the 4-byte magic "TBA\x02"
 *   • that n_blocks * rec_len bytes are present in Field 6
 *   • that Field 7 contains exactly n_blocks uint64_t values
 *
 * All integer fields in the file are stored in the native byte order used by
 * dump_blocks() (which writes them with fwrite on the same platform), so no
 * byte-swapping is performed here.
 */
static int sidx_load(const char *sidx_fname, sidx_t *out)
{
    FILE *sf = fopen(sidx_fname, "rb");
    if (!sf) {
        fprintf(stderr, "[tabby] sidx_load: cannot open '%s': %s\n",
                sidx_fname, strerror(errno));
        return -1;
    }

#define SIDX_READ(ptr, n) do { \
    if (fread((ptr), 1, (size_t)(n), sf) != (size_t)(n)) { \
        fprintf(stderr, "[tabby] sidx_load: read error in '%s': %s\n", \
                sidx_fname, feof(sf) ? "unexpected EOF" : strerror(errno)); \
        fclose(sf); \
        return -1; \
    } \
} while (0)

    /* Field 0: magic */
    char magic[4];
    SIDX_READ(magic, 4);
    if (magic[0] != 'T' || magic[1] != 'B' || magic[2] != 'A' || magic[3] != '\x02') {
        fprintf(stderr, "[tabby] sidx_load: bad magic in '%s' "
                "(expected 'TBA\\x02', got %02x %02x %02x %02x)\n",
                sidx_fname,
                (unsigned char)magic[0], (unsigned char)magic[1],
                (unsigned char)magic[2], (unsigned char)magic[3]);
        fclose(sf);
        return -1;
    }

    /* Field 1: n_columns */
    uint32_t ncols = 0;
    SIDX_READ(&ncols, 4);
    out->ncols = ncols;

    /* Field 2a: column_min_types[ncols] */
    if (ncols > 0) {
        out->col_min_tc = malloc((size_t)ncols);
        if (!out->col_min_tc) {
            fprintf(stderr, "[tabby] sidx_load: out of memory\n");
            fclose(sf);
            return -1;
        }
        SIDX_READ(out->col_min_tc, (size_t)ncols);
    }

    /* Field 2b: column_max_types[ncols] */
    if (ncols > 0) {
        out->col_max_tc = malloc((size_t)ncols);
        if (!out->col_max_tc) {
            fprintf(stderr, "[tabby] sidx_load: out of memory\n");
            fclose(sf);
            return -1;
        }
        SIDX_READ(out->col_max_tc, (size_t)ncols);
    }

    /* Field 3: n_blocks */
    SIDX_READ(&out->nblocks, 8);

    /* Field 3b: max_comp_block_size */
    SIDX_READ(&out->max_comp_block_size, 2);

    /* Field 4: block_record_length */
    SIDX_READ(&out->rec_len, 8);

    /* Field 5: blk_rec_offsets[nblocks] */
    if (out->nblocks > 0) {
        out->blk_rec_offsets = malloc((size_t)(out->nblocks * sizeof(uint64_t)));
        if (!out->blk_rec_offsets) {
            fprintf(stderr, "[tabby] sidx_load: out of memory\n");
            fclose(sf);
            return -1;
        }
        SIDX_READ(out->blk_rec_offsets, (size_t)(out->nblocks * sizeof(uint64_t)));
    }

    /* Field 6: block_records (nblocks * rec_len bytes) */
    {
        size_t total_records = (size_t)(out->nblocks * out->rec_len);
        if (total_records > 0) {
            out->records = malloc(total_records);
            if (!out->records) {
                fprintf(stderr, "[tabby] sidx_load: out of memory\n");
                fclose(sf);
                return -1;
            }
            SIDX_READ(out->records, total_records);
        }
    }

    /* Field 7: blk_comp_offsets[nblocks] */
    if (out->nblocks > 0) {
        out->blk_comp_offsets = malloc((size_t)(out->nblocks * sizeof(uint64_t)));
        if (!out->blk_comp_offsets) {
            fprintf(stderr, "[tabby] sidx_load: out of memory\n");
            fclose(sf);
            return -1;
        }
        SIDX_READ(out->blk_comp_offsets, (size_t)(out->nblocks * sizeof(uint64_t)));
    }

#undef SIDX_READ

    fclose(sf);

    fprintf(stderr, "[tabby] sidx_load: loaded '%s' "
            "(%"PRIu32" col(s), %"PRIu64" block(s), %"PRIu64" bytes/record, "
            "max_comp_block=%"PRIu16")\n",
            sidx_fname, out->ncols, out->nblocks, out->rec_len,
            out->max_comp_block_size);
    return 0;
}

/*
 * dump_blocks() uses TWO passes over the BGZF file to guarantee that every
 * block index record in the .sidx output has exactly the same byte length.
 *
 * PASS 1 – type discovery
 *   Scan every data line, tokenise it, and feed each field into a file-level
 *   ColStat array (file_col[]).  No per-block tracking, no output.
 *   After pass 1 we know:
 *     • ncols          – exact column count
 *     • file_col_type[]– the widest ColType seen for each column across the
 *                        entire file (derived from file_col[])
 *     • rec_len        – the fixed byte length of every block record, computed
 *                        once from file_col_type[] and never changed again
 *
 * PASS 2 – block statistics + output
 *   Re-open the file and scan again.  Per-block ColStat (blk_col[]) is reset
 *   at each block boundary.  When FLUSH_BLOCK() serialises a block record it
 *   uses file_col_type[ci] (the frozen file-level type) as the wire type for
 *   every column, so every record is exactly rec_len bytes.  The TSV summary
 *   and the .sidx binary file are both written during this pass.
 *
 * This eliminates the inconsistency that existed in the single-pass approach,
 * where blocks seen earlier in the file might have a narrower per-block type
 * than the file-level (widest) type that ends up in the header.
 */
static int dump_blocks(const char *fname, int download)
{
    int ftype = file_type(fname);
    if (!(ftype & IS_TXT) && ftype != 0) {
        fprintf(stderr, "[tabby] --dump-blocks only supports text (TBI) indexed files\n");
        return 1;
    }

    /* Load the tabix index so we can resolve tid -> name. */
    tbx_t *tbx = tbx_index_load3(fname, NULL, download ? HTS_IDX_SAVE_REMOTE : 0);
    if (!tbx) {
        fprintf(stderr, "[tabby] Could not load .tbi index of %s\n", fname);
        return 1;
    }

    int nseq = 0;
    const char **seqnames = tbx_seqnames(tbx, &nseq);
    if (!seqnames) {
        fprintf(stderr, "[tabby] Could not get sequence names\n");
        tbx_destroy(tbx);
        return 1;
    }

    /* ── Shared allocations (used in both passes) ── */
    int       ncols         = 0;
    int       ncols_cap     = 32;

    /* file_col[ci]  – file-level ColStat, accumulated during pass 1 only.
     * file_col_type[ci] – the widest ColType for column ci (derived from
     *                     file_col[] after pass 1; frozen before pass 2). */
    ColStat  *file_col      = malloc(ncols_cap * sizeof(ColStat));
    ColType  *file_col_type = malloc(ncols_cap * sizeof(ColType));

    /* blk_col[ci]   – per-block ColStat, reset at each block boundary in
     *                 pass 2.  Not used during pass 1. */
    ColStat  *blk_col       = malloc(ncols_cap * sizeof(ColStat));

    if (!file_col || !file_col_type || !blk_col) {
        fprintf(stderr, "[tabby] out of memory\n");
        free(file_col); free(file_col_type); free(blk_col);
        free(seqnames); tbx_destroy(tbx);
        return 1;
    }

    /* Initialize all pre-allocated entries; col_update() checks stat->seen to
     * distinguish first-occurrence from widening — garbage there corrupts types. */
    for (int ci = 0; ci < ncols_cap; ci++) {
        file_col[ci]      = colstat_init();
        file_col_type[ci] = CT_BOOL;
        blk_col[ci]       = colstat_init();
    }

    /* Mutable line buffer shared across both passes. */
    char     *linebuf     = NULL;
    size_t    linebuf_cap = 0;
    kstring_t str         = {0, 0, 0};

    /* Helper: grow all three column arrays to at least 'need' entries. */
#define ENSURE_COLS(need) do { \
    if ((need) > ncols_cap) { \
        int _newcap = ncols_cap; \
        while (_newcap < (need)) _newcap *= 2; \
        ColStat *_fc  = realloc(file_col,       _newcap * sizeof(ColStat)); \
        ColType *_fct = realloc(file_col_type,  _newcap * sizeof(ColType)); \
        ColStat *_bc  = realloc(blk_col,        _newcap * sizeof(ColStat)); \
        if (!_fc || !_fct || !_bc) { \
            fprintf(stderr, "[tabby] out of memory\n"); \
            free(_fc  ? _fc  : file_col);       \
            free(_fct ? _fct : file_col_type);  \
            free(_bc  ? _bc  : blk_col);        \
            file_col = NULL; file_col_type = NULL; blk_col = NULL; \
            goto done; \
        } \
        file_col = _fc; file_col_type = _fct; blk_col = _bc; \
        for (int _ci = ncols_cap; _ci < _newcap; _ci++) { \
            file_col[_ci]      = colstat_init(); \
            file_col_type[_ci] = CT_BOOL; /* will be widened on first use */ \
            blk_col[_ci]       = colstat_init(); \
        } \
        ncols_cap = _newcap; \
    } \
} while (0)

    /* ════════════════════════════════════════════════════════════════════════
     * PASS 1 – scan the whole file to determine the file-level (widest) type
     *           for every column.
     * ════════════════════════════════════════════════════════════════════════ */
    {
        BGZF *fp1 = bgzf_open(fname, "r");
        if (!fp1) {
            fprintf(stderr, "[tabby] Could not open %s (pass 1)\n", fname);
            goto done;
        }

        int64_t lineno1 = 0;
        int     ret1;
        while ((ret1 = bgzf_getline(fp1, '\n', &str)) >= 0) {
            ++lineno1;

            /* Skip header/comment lines and lines before line_skip. */
            if (!str.l || str.s[0] == tbx->conf.meta_char
                       || lineno1 <= tbx->conf.line_skip)
                continue;

            /* Only process lines that tbx_parse1 accepts as data records. */
            tbx_intv_t intv1;
            if (tbx_parse1(&tbx->conf, str.l, str.s, &intv1) != 0)
                continue;

            /* Tokenise into fields; update file_col[] for each column. */
            if (str.l + 1 > linebuf_cap) {
                size_t newcap = (str.l + 1) * 2;
                char *nb = realloc(linebuf, newcap);
                if (!nb) {
                    fprintf(stderr, "[tabby] out of memory\n");
                    bgzf_close(fp1);
                    goto done;
                }
                linebuf = nb;
                linebuf_cap = newcap;
            }
            memcpy(linebuf, str.s, str.l + 1);

            int   col_idx  = 0;
            char *p        = linebuf;
            char *line_end = linebuf + str.l;
            while (p <= line_end) {
                char *tab       = memchr(p, '\t', (size_t)(line_end - p));
                char *field_end = tab ? tab : line_end;
                *field_end = '\0';

                ENSURE_COLS(col_idx + 1);
                if (col_idx >= ncols)
                    ncols = col_idx + 1;

                col_update(&file_col[col_idx], p, (size_t)(field_end - p));

                if (!tab) break;
                p = tab + 1;
                col_idx++;
            }
        }
        bgzf_close(fp1);

        if (ret1 < -1) {
            fprintf(stderr, "[tabby] Read error during pass 1 of %s\n", fname);
            goto done;
        }
    }

    /* Derive file_col_type[] from file_col[] now that pass 1 is complete. */
    for (int ci = 0; ci < ncols; ci++)
        file_col_type[ci] = file_col[ci].seen ? file_col[ci].type : CT_STRING;

    /* Compute the fixed record length from the frozen file-level types. */
    uint64_t rec_len = 0;
    for (int ci = 0; ci < ncols; ci++) {
        uint8_t tc = coltype_to_sidx(file_col_type[ci]);
        if (tc == 0)
            rec_len += 1;
        else
            rec_len += sidx_type_size(tc) * 2;
    }

    /* ════════════════════════════════════════════════════════════════════════
     * PASS 2 – re-scan the file; emit TSV output and accumulate sidx records.
     * ════════════════════════════════════════════════════════════════════════ */

    /* Re-initialise per-block ColStat array for pass 2. */
    for (int ci = 0; ci < ncols; ci++)
        blk_col[ci] = colstat_init();

    /* -- per-block accumulators -- */
    int64_t  blk_compressed_off = 0;
    int64_t  blk_coord_beg      = -1;
    int64_t  blk_coord_end      = -1;
    int64_t  blk_nlines         = 0;
    int64_t  blk_uncomp_bytes   = 0;
    int64_t  last_utell         = 0;
    int64_t  new_utell          = 0;

    /* Largest compressed block size seen across the whole file (Field 3b). */
    uint16_t max_comp_block_size  = 0;
    /* Compressed file offset of the block *following* the one being flushed;
     * set by the caller just before each FLUSH_BLOCK() invocation so the
     * macro can compute this block's compressed size. */
    uint64_t flush_comp_block_end = 0;

    int      blk_ntids    = 0;
    int      blk_tids_cap = 16;
    int     *blk_tids     = malloc(blk_tids_cap * sizeof(int));
    if (!blk_tids) {
        fprintf(stderr, "[tabby] out of memory\n");
        goto done;
    }

    /*
     * Secondary index in-memory accumulation.
     * Every record is exactly rec_len bytes (guaranteed by using file_col_type[]
     * as the wire type in FLUSH_BLOCK, not the narrower per-block type).
     */
    uint8_t  *sidx_records     = NULL;
    size_t    sidx_records_len = 0;
    size_t    sidx_records_cap = 0;
    uint64_t *sidx_offsets     = NULL;   /* byte offset of each block record within sidx_records[] */
    uint64_t *sidx_blk_offsets = NULL;   /* compressed file offset (upper 48 bits of voff) per block */
    uint64_t  sidx_nblocks     = 0;
    uint64_t  sidx_nblocks_cap = 0;

    /* Helper: reset per-block column stats between blocks in pass 2. */
#define RESET_BLK_COLS() do { \
    for (int _ci = 0; _ci < ncols; _ci++) \
        blk_col[_ci] = colstat_init(); \
} while (0)

    /* Print TSV header. */
    printf("compressed_offset\tseq\tcoord_start\tcoord_end\tn_lines\tuncompressed_bytes\n");

    int64_t lineno = 0;

    /*
     * FLUSH_BLOCK() – called at every block boundary and at end-of-file.
     *
     * A. Print the human-readable TSV summary line + per-column stats.
     * B. Serialise the block record into sidx_records[] using the
     *    file-level type (file_col_type[ci]) as the wire type for each
     *    column, so every record is exactly rec_len bytes.
     *    Layout per column ci:
     *      if file_col_type sidx code == 0:  1 zero sentinel byte
     *      else:  min_value (sidx_type_size bytes) + max_value (same)
     *    Columns absent in this block (!seen) are written as zero bytes
     *    of the correct file-level width, keeping the record length fixed.
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
        } \
        /* ── B: serialise block record using frozen file-level types ── */ \
        { \
            /* Grow sidx_records[] if needed. */ \
            if (sidx_records_len + (size_t)rec_len > sidx_records_cap) { \
                size_t _newcap = sidx_records_cap ? sidx_records_cap * 2 : 65536; \
                while (_newcap < sidx_records_len + (size_t)rec_len) _newcap *= 2; \
                uint8_t *_nb = realloc(sidx_records, _newcap); \
                if (!_nb) { \
                    fprintf(stderr, "[tabby] out of memory (sidx records)\n"); \
                    goto done; \
                } \
                sidx_records = _nb; \
                sidx_records_cap = _newcap; \
            } \
            /* Grow sidx_offsets[] and sidx_blk_offsets[] if needed. */ \
            if (sidx_nblocks + 1 > sidx_nblocks_cap) { \
                uint64_t _newcap2 = sidx_nblocks_cap ? sidx_nblocks_cap * 2 : 1024; \
                uint64_t *_no = realloc(sidx_offsets, (size_t)_newcap2 * sizeof(uint64_t)); \
                uint64_t *_nb = realloc(sidx_blk_offsets, (size_t)_newcap2 * sizeof(uint64_t)); \
                if (!_no || !_nb) { \
                    fprintf(stderr, "[tabby] out of memory (sidx offsets)\n"); \
                    goto done; \
                } \
                sidx_offsets     = _no; \
                sidx_blk_offsets = _nb; \
                sidx_nblocks_cap = _newcap2; \
            } \
            sidx_offsets[sidx_nblocks]     = (uint64_t)sidx_records_len; \
            sidx_blk_offsets[sidx_nblocks] = (uint64_t)blk_compressed_off; \
            sidx_nblocks++; \
            { \
                uint64_t _csz = flush_comp_block_end - (uint64_t)blk_compressed_off; \
                if (_csz > max_comp_block_size) \
                    max_comp_block_size = (_csz > 0xFFFF) ? 0xFFFF : (uint16_t)_csz; \
            } \
            uint8_t *_p = sidx_records + sidx_records_len; \
            for (int _ci = 0; _ci < ncols; _ci++) { \
                uint8_t  _ftc = coltype_to_sidx(file_col_type[_ci]); \
                ColStat *_cs2 = &blk_col[_ci]; \
                if (_ftc == 0) { \
                    *_p++ = 0x00; /* string/NA: 1 sentinel zero byte */ \
                } else if (!_cs2->seen) { \
                    /* Column absent in this block: zero-fill to keep length fixed. */ \
                    size_t _vsz = sidx_type_size(_ftc); \
                    memset(_p, 0, _vsz * 2); \
                    _p += _vsz * 2; \
                } else { \
                    _p += sidx_write_min(_ftc, _cs2, _p); \
                    _p += sidx_write_max(_ftc, _cs2, _p); \
                } \
            } \
            sidx_records_len += (size_t)rec_len; \
        } \
    } \
} while (0)

    BGZF *fp = bgzf_open(fname, "r");
    if (!fp) {
        fprintf(stderr, "[tabby] Could not open %s (pass 2)\n", fname);
        goto done;
    }

    uint64_t voff_after = 0;
    blk_compressed_off  = 0;

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
                    new_utell        = bgzf_utell(fp);
                    blk_uncomp_bytes = new_utell - last_utell;
                    last_utell       = new_utell;
                    flush_comp_block_end = (uint64_t)cur_blk;
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
                            if (!tmp) {
                                fprintf(stderr, "[tabby] out of memory\n");
                                bgzf_close(fp);
                                goto done;
                            }
                            blk_tids = tmp;
                        }
                        blk_tids[blk_ntids++] = tid;
                    }
                }

                /* ── Tokenise line and update per-block column statistics ── */
                if (str.l + 1 > linebuf_cap) {
                    size_t newcap = (str.l + 1) * 2;
                    char *nb = realloc(linebuf, newcap);
                    if (!nb) {
                        fprintf(stderr, "[tabby] out of memory\n");
                        bgzf_close(fp);
                        goto done;
                    }
                    linebuf = nb;
                    linebuf_cap = newcap;
                }
                memcpy(linebuf, str.s, str.l + 1);

                int   col_idx  = 0;
                char *p2       = linebuf;
                char *line_end = linebuf + str.l;
                while (p2 <= line_end) {
                    char *tab       = memchr(p2, '\t', (size_t)(line_end - p2));
                    char *field_end = tab ? tab : line_end;
                    *field_end = '\0';

                    /* ncols is frozen after pass 1; ENSURE_COLS guards against
                     * malformed files that somehow have more columns in pass 2. */
                    ENSURE_COLS(col_idx + 1);
                    if (col_idx >= ncols) ncols = col_idx + 1;

                    col_update(&blk_col[col_idx], p2, (size_t)(field_end - p2));

                    if (!tab) break;
                    p2 = tab + 1;
                    col_idx++;
                }
            }
        }

        if (!is_data && crossed) {
            blk_uncomp_bytes = bgzf_utell(fp) - last_utell;
            flush_comp_block_end = (uint64_t)cur_blk;
            FLUSH_BLOCK();
            blk_compressed_off = cur_blk;
            blk_coord_beg    = -1;
            blk_coord_end    = -1;
            blk_nlines       = 1;   /* count this non-data line in the new block */
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
    flush_comp_block_end = (uint64_t)bgzf_tell(fp) >> 16;
    FLUSH_BLOCK();

    bgzf_close(fp);

    /* ── Print whole-file column type summary ───────────────────────────── */
    if (ncols > 0) {
        printf("\n# File-level column type summary (%d column(s)):\n", ncols);
        for (int ci = 0; ci < ncols; ci++)
            printf("#   col_%d: %s\n", ci, coltype_name(file_col_type[ci]));
    }

    /* ── Write the binary secondary index (.sidx) ───────────────────────── */
    /*
     * All fields are now available:
     *   n_columns            = ncols           (from pass 1)
     *   column_min/max_types = file_col_type[]  (frozen before pass 2)
     *   n_blocks             = sidx_nblocks    (accumulated in pass 2)
     *   block_record_length  = rec_len         (fixed; derived from file_col_type[])
     *   offsets[]            = sidx_offsets[]  (accumulated in pass 2)
     *   block_records        = sidx_records[]  (accumulated in pass 2)
     *
     * Every record in sidx_records[] is exactly rec_len bytes because
     * FLUSH_BLOCK() used file_col_type[] as the wire type throughout pass 2.
     */
    {
        uint8_t *col_min_tc = malloc(ncols > 0 ? (size_t)ncols : 1);
        uint8_t *col_max_tc = malloc(ncols > 0 ? (size_t)ncols : 1);
        if (!col_min_tc || !col_max_tc) {
            fprintf(stderr, "[tabby] out of memory (sidx header)\n");
            free(col_min_tc); free(col_max_tc);
            goto done;
        }
        for (int ci = 0; ci < ncols; ci++) {
            uint8_t tc = coltype_to_sidx(file_col_type[ci]);
            col_min_tc[ci] = tc;
            col_max_tc[ci] = tc;
        }

        size_t fnlen = strlen(fname);
        char *sidx_fname = malloc(fnlen + 6);
        if (!sidx_fname) {
            fprintf(stderr, "[tabby] out of memory (sidx filename)\n");
            free(col_min_tc); free(col_max_tc);
            goto done;
        }
        memcpy(sidx_fname, fname, fnlen);
        memcpy(sidx_fname + fnlen, ".sidx", 6); /* includes NUL */

        FILE *sf = fopen(sidx_fname, "wb");
        if (!sf) {
            fprintf(stderr, "[tabby] Could not open %s for writing: %s\n",
                    sidx_fname, strerror(errno));
            free(col_min_tc); free(col_max_tc); free(sidx_fname);
            goto done;
        }

#define SIDX_WRITE(ptr, n) do { \
    if (fwrite((ptr), 1, (n), sf) != (size_t)(n)) { \
        fprintf(stderr, "[tabby] Write error on %s: %s\n", sidx_fname, strerror(errno)); \
        fclose(sf); free(col_min_tc); free(col_max_tc); free(sidx_fname); \
        goto done; \
    } \
} while (0)

        /* Field 0: magic "TBA\x02" */
        const char magic[4] = {'T', 'B', 'A', '\x02'};
        SIDX_WRITE(magic, 4);

        /* Field 1: n_columns (uint32_t) */
        { uint32_t nc = (uint32_t)ncols; SIDX_WRITE(&nc, 4); }

        /* Field 2a: column_min_types[ncols] */
        if (ncols > 0) SIDX_WRITE(col_min_tc, (size_t)ncols);

        /* Field 2b: column_max_types[ncols] */
        if (ncols > 0) SIDX_WRITE(col_max_tc, (size_t)ncols);

        /* Field 3: n_blocks (uint64_t) */
        SIDX_WRITE(&sidx_nblocks, 8);

        /* Field 3b: max_comp_block_size (uint16_t) --
         * Largest compressed BGZF block size in bytes across the whole file.
         * Capped at UINT16_MAX (65535) if a block somehow exceeds that. */
        SIDX_WRITE(&max_comp_block_size, 2);

        /* Field 4: block_record_length (uint64_t) */
        SIDX_WRITE(&rec_len, 8);

        /* Field 5: block_record_uncompressed_byte_offsets (uint64_t[n_blocks]) */
        if (sidx_nblocks > 0)
            SIDX_WRITE(sidx_offsets, (size_t)(sidx_nblocks * sizeof(uint64_t)));

        /* Field 6: block_records (n_blocks * rec_len bytes) */
        if (sidx_records_len > 0)
            SIDX_WRITE(sidx_records, sidx_records_len);

        /* Field 7: compressed_offsets[n_blocks] (uint64_t[]) --
         * BGZF compressed file byte offset for each block (upper 48 bits of
         * the BGZF virtual offset, i.e. voff >> 16).  Written in the same
         * block order as Field 5 / Field 6.  Sorted ascending because BGZF
         * blocks are visited sequentially during dump_blocks().
         * Readers can binary-search this array with
         * sidx_block_index_from_offset() to map itr->curr_off >> 16 to a
         * block index, then use Field 5 to locate the block's column-stat
         * record within Field 6. */
        if (sidx_nblocks > 0)
            SIDX_WRITE(sidx_blk_offsets, (size_t)(sidx_nblocks * sizeof(uint64_t)));

#undef SIDX_WRITE

        fclose(sf);
        fprintf(stderr, "[tabby] wrote secondary index: %s "
                "(%d col(s), %"PRIu64" block(s), %"PRIu64" bytes/record, "
                "max_comp_block=%"PRIu16", compressed_offsets included)\n",
                sidx_fname, ncols, sidx_nblocks, rec_len, max_comp_block_size);

        free(col_min_tc);
        free(col_max_tc);
        free(sidx_fname);
    }

#undef FLUSH_BLOCK
#undef RESET_BLK_COLS
#undef ENSURE_COLS

done:
    if(sidx_records)
        free(sidx_records);
    if(sidx_offsets)
       free(sidx_offsets);
    if(sidx_blk_offsets)
        free(sidx_blk_offsets);
    if(blk_tids)
        free(blk_tids);
    if(linebuf)
        free(linebuf);
    if(str.s)
        free(str.s);
    if(blk_col)
        free(blk_col);
    if(file_col)
        free(file_col);
    if(file_col_type)
        free(file_col_type);
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
    fprintf(fp, "Usage:   tabby [OPTIONS] [FILE] [REGION [...]]\n");
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
    fprintf(fp, "   -F, --filter EXPR          filter records by column value (AND); EXPR format:\n");
    fprintf(fp, "                              <col><op><val>, col is 0-based index or column name.\n");
    fprintf(fp, "                              Column names are resolved from the file's first\n");
    fprintf(fp, "                              header line (leading # stripped if present).\n");
    fprintf(fp, "                              Numeric ops: ==  !=  >=  <=\n");
    fprintf(fp, "                              String ops:  ==  !=  (exact)  ~=  !~  (ERE regex)\n");
    fprintf(fp, "                              Repeatable. Examples: -F \"3<=12.5\" -F \"SCORE>=50\"\n");
    fprintf(fp, "   -O, --or-filter EXPR       like -F but expressions are OR'd; the AND group\n");
    fprintf(fp, "                              (-F) and OR group (-O) are combined with AND.\n");
    fprintf(fp, "       --verbosity INT        set verbosity [3]\n");
    fprintf(fp, "   -@, --threads INT          number of additional threads to use [0]\n");
    fprintf(fp, "\n");
    fprintf(fp, "The default preset is 'seqonly': indexes tab-delimited files by\n");
    fprintf(fp, "sequence/chromosome name only. No start/end coordinate columns are\n");
    fprintf(fp, "needed. Use -s to specify the chromosome column (default: 1).\n");
    fprintf(fp, "Query with a bare name: tabby FILE chr1\n");
    fprintf(fp, "Use -p gff/bed/sam/vcf to index coordinate-based files.\n");
    fprintf(fp, "\n");
    return status;
}

int main(int argc, char *argv[])
{
    if (argc == 2 && strcmp(argv[1], "cat") == 0) {
        printf("\n");
        printf("      /\\     /\\\n");
        printf("     (  \\___/  )\n");
        printf("     (  M   M  )\n");
        printf("     (  o   o  )\n");
        printf("      \\   ^   /    ~ meow ~\n");
        printf("       ) ~|~ (\n");
        printf("      /  |||  \\\n");
        printf("     (__|___|__)\n");
        printf("\n");
        printf("         tabby\n");
        printf("\n");
        return EXIT_SUCCESS;
    }

    int c, detect = 0, min_shift = 0, is_force = 0, list_chroms = 0, do_csi = 0, dump_blocks_flag = 0;
    tbx_conf_t conf = tbx_conf_seqonly;
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
        {"filter", required_argument, NULL, 'F'},
        {"or-filter", required_argument, NULL, 'O'},
        {NULL, 0, NULL, 0}
    };

    char *tmp;
    while ((c = getopt_long(argc, argv, "hH?0b:c:e:fm:p:s:S:lr:CR:T:D@:F:O:", loptions,NULL)) >= 0)
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
                conf.preset &= ~TBX_SEQONLY;   /* explicit coord column → not seqonly */
                detect = 0;
                break;
            case 'e':
                conf.ec = strtol(optarg,&tmp,10);
                if ( *tmp ) error("Could not parse argument: -e %s\n", optarg);
                conf.preset &= ~TBX_SEQONLY;   /* explicit coord column → not seqonly */
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
"Tabby (htslib) %s\n"
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
            case 'F':
                if (tbx_parse_filter(&args, optarg, 0) != 0)
                    return EXIT_FAILURE;
                break;
            case 'O':
                if (tbx_parse_filter(&args, optarg, 1) != 0)
                    return EXIT_FAILURE;
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
    if ( argc > optind+1 || args.header_only || args.regions_fname || args.targets_fname || args.nfilters > 0 || args.n_or_filters > 0 )
    {
        int nregs = 0;
        char **regs = NULL;
        if ( !args.header_only )
            regs = parse_regions(args.regions_fname, argv+optind+1, argc-optind-1, &nregs);

        /* Attempt to load the secondary index if present.
         * Failure is non-fatal: query_regions() handles sidx == NULL gracefully. */
        sidx_t sidx;
        memset(&sidx, 0, sizeof(sidx));
        sidx_t *sidx_ptr = NULL;
        {
            size_t fnlen = strlen(fname);
            char *sidx_fname = malloc(fnlen + 6);   /* ".sidx" + NUL */
            if (sidx_fname) {
                memcpy(sidx_fname, fname, fnlen);
                memcpy(sidx_fname + fnlen, ".sidx", 6);
                if (sidx_load(sidx_fname, &sidx) == 0)
                    sidx_ptr = &sidx;
                free(sidx_fname);
            }
        }

        int qret = query_regions(&args, &conf, fname, regs, nregs, sidx_ptr);
        sidx_free(&sidx);
        for (int i = 0; i < args.nfilters; i++) {
            if (args.filters[i].has_regex)
                regfree(&args.filters[i].compiled_re);
            free(args.filters[i].sval);
            free(args.filters[i].col_name);
        }
        free(args.filters);
        return qret;
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
            error("[tabby] the index file exists. Please use '-f' to overwrite.\n");
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
                error("[tabby] the compression of '%s' is not BGZF\n", fname);
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
                error("[tabby] the compression of '%s' is not BGZF\n", fname);
            default:
                error("tbx_index_build3 failed: %s\n", fname);
        }
    }

    return 0;
}
