#!/bin/bash
# test-sidx.sh -- Integration tests for the tabix secondary-index (sidx).
#
# Functions exercised:
#   dump_blocks()          -- creates .sidx and emits per-block TSV
#   sidx_load()            -- loads .sidx on query startup; graceful on failure
#   sidx_free()            -- called automatically after every query
#   tbx_filter_chunks()    -- chunk-level pruning with all four filter ops
#   sidx_block_filter_cb() -- block-level callback for local-file queries
#   tbx_apply_filters()    -- per-row exact check; string/regex/numeric/OOB column handling
#   tbx_parse_filter()     -- -F expression parsing; all error paths including
#                             string filters (==, !=) and regex filters (~=, !~)
#   coltype_to_sidx()      -- type-code mapping (exercised by dump_blocks)
#   sidx_write_min/max()   -- value serialisation (exercised by dump_blocks)
#   col_update()           -- type narrowing; boolean/int32/float/double paths
#
# Functions NOT exercised here (no CLI entry point or require remote HTTP):
#   sidx_block_index_from_offset()         -- static utility, not yet called
#   sidx_block_index_from_offset_bounded() -- static utility, not yet called
#   bgzf_decomp_skip_cb()                  -- remote-file path only
#
# Run from test/tabix/:
#   cd test/tabix && ./test-sidx.sh
#
# Requires tabix and bgzip compiled in ../../

set -u

TABIX="./tabby"
BGZIP="./bgzip"
TMP="${TMPDIR:-/tmp}/sidx_test_$$"

_np=0; _nf=0

pass() { echo "PASS : $1"; ((_np++)); }
fail() { echo "FAIL : $1"; ((_nf++)); }

# check <name> <cmd> [args...]  -- run command, expect exit 0
check() {
    local n="$1"; shift
    if "$@" >"${TMP}.out" 2>"${TMP}.err"; then
        pass "$n"
    else
        fail "$n (exit $?)"
        head -2 "${TMP}.err" | sed 's/^/       /' >&2
    fi
}

# check_fail <name> <cmd> [args...]  -- run command, expect non-zero exit
check_fail() {
    local n="$1"; shift
    if "$@" >"${TMP}.out" 2>"${TMP}.err"; then
        fail "$n (expected non-zero exit)"
    else
        pass "$n"
    fi
}

# check_out <name> <expected_file> <cmd> [args...]  -- compare stdout to file
check_out() {
    local n="$1" exp="$2"; shift 2
    "$@" >"${TMP}.actual" 2>"${TMP}.err"
    local ec=$?
    if [ "$ec" -ne 0 ]; then
        fail "$n (exit $ec)"
        head -2 "${TMP}.err" | sed 's/^/       /' >&2
    elif cmp -s "$exp" "${TMP}.actual"; then
        pass "$n"
    else
        local el al
        el=$(wc -l < "$exp"); al=$(wc -l < "${TMP}.actual")
        fail "$n (output mismatch: expected $el lines, got $al lines)"
        diff "$exp" "${TMP}.actual" | head -5 | sed 's/^/       /' >&2
    fi
}

cleanup() {
    rm -f sidx_data.tsv sidx_data.tsv.gz sidx_data.tsv.gz.tbi sidx_data.tsv.gz.sidx
    rm -f sidx_types.tsv sidx_types.tsv.gz sidx_types.tsv.gz.tbi sidx_types.tsv.gz.sidx
    rm -f sidx_same.tsv sidx_same.tsv.gz sidx_same.tsv.gz.tbi sidx_same.tsv.gz.sidx
    rm -f sidx_seqonly.tsv sidx_seqonly.tsv.gz sidx_seqonly.tsv.gz.tbi
    rm -f sidx_named.tsv sidx_named.tsv.gz sidx_named.tsv.gz.tbi sidx_named.tsv.gz.sidx
    rm -f sidx_coords.tsv sidx_coords.tsv.gz sidx_coords.tsv.gz.tbi
    rm -f sidx_notbi.tsv.gz
    rm -f sidx_dump.out sidx_types_dump.out sidx_same_dump.out
    rm -f sidx_exp_*.tsv
    rm -f sidx_exp_str_*.tsv
    rm -f sidx_exp_re_*.tsv
    rm -f "${TMP}".*
}

trap cleanup EXIT

die() { echo "INIT FAIL: $*" >&2; exit 1; }

[ -x "$TABIX" ] || die "tabix not found at $TABIX"
[ -x "$BGZIP" ] || die "bgzip not found at $BGZIP"

echo "Testing tabix secondary index (sidx)..."
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# SETUP: generate three test files
# ─────────────────────────────────────────────────────────────────────────────
echo "Generating test data..."

# ── sidx_data.tsv: main test file, 8000 records spanning 5-6 BGZF blocks ──
#
# 6 columns (0-based):
#   0: CHROM  string  (chr1 / chr2)
#   1: POS    int16   (1000–4999)
#   2: END    int16   (= POS, point records)
#   3: NAME   string  (row_<pos>)
#   4: SCORE  int16   (chr1: POS, range 1000–4999;
#                      chr2: POS+10000, range 11000–14999)
#   5: QUAL   float   (chr1: POS/100.0,       range 10.00–49.99;
#                      chr2: (POS-1000)/100.0, range  0.00–39.99)
#
# ~46 bytes/line × 8000 = ~368 KB → 5–6 BGZF blocks at default 64 KB block size.
# SCORE gap (chr1 max=4999, chr2 min=11000) enables block-level pruning tests.
awk 'BEGIN {
    for (i = 1000; i <= 4999; i++)
        printf "chr1\t%d\t%d\trow_%d\t%d\t%.2f\n", i, i, i, i, i/100.0
    for (i = 1000; i <= 4999; i++)
        printf "chr2\t%d\t%d\trow_%d\t%d\t%.2f\n", i, i, i, i+10000, (i-1000)/100.0
}' > sidx_data.tsv

$BGZIP -c sidx_data.tsv > sidx_data.tsv.gz  || die "bgzip sidx_data"
$TABIX -f -s 1 -b 2 -e 3 -c '#' sidx_data.tsv.gz || die "tabix index sidx_data"

# ── sidx_types.tsv: type-inference test, one BGZF block, 7 columns ──
#
# Covers boolean (CT_BOOL), int32 (CT_INT32), float (CT_FLOAT), double (CT_DOUBLE).
#   0: CHROM  string
#   1: POS    int16   (1, 2)
#   2: END    int16   (= POS)
#   3: BVAL   boolean (0 or 1 — single-char 0/1 → CT_BOOL)
#   4: BIG    int32   (100000 / 200000 — exceeds int16 range 32767)
#   5: FP     float   (3.14 / 2.72 — ≤7 significant digits → CT_FLOAT)
#   6: DP     double  (3.14159265358979 / 2.71828182845905 — >7 sig digits → CT_DOUBLE)
printf 'chr1\t1\t1\t0\t100000\t3.14\t3.14159265358979\n' >  sidx_types.tsv
printf 'chr1\t2\t2\t1\t200000\t2.72\t2.71828182845905\n' >> sidx_types.tsv

$BGZIP -c sidx_types.tsv > sidx_types.tsv.gz  || die "bgzip sidx_types"
$TABIX -f -s 1 -b 2 -e 3 -c '#' sidx_types.tsv.gz || die "tabix index sidx_types"

# ── sidx_same.tsv: all records have SCORE=42 — tests FILT_NE block pruning ──
#
# All three records land in one BGZF block.  dump_blocks records blk_min=42,
# blk_max=42.  A filter "4!=42" triggers the block-impossibility condition
# (blk_min == val && blk_max == val) in tbx_filter_chunks, dropping the whole
# chunk before any I/O.
printf 'chr1\t1\t1\tA\t42\t1.00\n' >  sidx_same.tsv
printf 'chr1\t2\t2\tB\t42\t2.00\n' >> sidx_same.tsv
printf 'chr1\t3\t3\tC\t42\t3.00\n' >> sidx_same.tsv

$BGZIP -c sidx_same.tsv > sidx_same.tsv.gz  || die "bgzip sidx_same"
$TABIX -f -s 1 -b 2 -e 3 -c '#' sidx_same.tsv.gz || die "tabix index sidx_same"

# ── Pre-generate expected output files from the uncompressed source ──
awk '$1 == "chr1"'                sidx_data.tsv > sidx_exp_chr1_all.tsv
awk '$5 == 2500'                  sidx_data.tsv > sidx_exp_score_eq2500.tsv
awk '$1=="chr1" && $5 <= 1010'    sidx_data.tsv > sidx_exp_score_le1010.tsv
awk '$1=="chr2" && $5 >= 14990'   sidx_data.tsv > sidx_exp_score_ge14990.tsv
awk '$1=="chr1" && $6 >= 49.0'    sidx_data.tsv > sidx_exp_qual_ge49.tsv
> sidx_exp_empty.tsv  # empty file: expected output when filter eliminates all

# ── Expected files for string filter tests ──
awk '$1=="chr1" && $2==2000'      sidx_data.tsv > sidx_exp_str_name_2000.tsv
printf 'chr1\t1\t1\tA\t42\t1.00\nchr1\t3\t3\tC\t42\t3.00\n' > sidx_exp_str_same_neqB.tsv
printf 'chr1\t2\t2\tB\t42\t2.00\n'                            > sidx_exp_str_same_eqB.tsv

# ── Expected files for regex filter tests ──
# ~= ^row_1[0-9][0-9][0-9]$ on chr1: rows 1000-1999 (NAME matches ^row_1...)
awk '$1=="chr1" && $2>=1000 && $2<=1999' sidx_data.tsv > sidx_exp_re_row1xxx.tsv
# !~ on NAME: chr1 rows whose NAME does NOT start with row_1
awk '$1=="chr1" && $2>=2000 && $2<=4999' sidx_data.tsv > sidx_exp_re_not_row1xxx.tsv

# ── sidx_seqonly.tsv: seqonly-indexed file (-p seqonly) ──
#
# 5 columns (0-based):
#   0: GENE   string (the indexed key — treated as the "chromosome")
#   1: CHROM  string
#   2: START  int
#   3: END    int
#   4: TYPE   string (protein_coding / lncRNA / rRNA)
#
# Mix of types so string/regex filters have meaningful pass/fail behaviour.
printf 'GeneA\tchr1\t1000\t2000\tprotein_coding\n' >  sidx_seqonly.tsv
printf 'GeneB\tchr1\t3000\t4000\tlncRNA\n'         >> sidx_seqonly.tsv
printf 'GeneC\tchr2\t5000\t6000\tprotein_coding\n' >> sidx_seqonly.tsv
printf 'GeneD\tchr2\t7000\t8000\trRNA\n'           >> sidx_seqonly.tsv
printf 'GeneE\tchr3\t9000\t9500\tprotein_coding\n' >> sidx_seqonly.tsv

$BGZIP -c sidx_seqonly.tsv > sidx_seqonly.tsv.gz || die "bgzip sidx_seqonly"
$TABIX -f -p seqonly sidx_seqonly.tsv.gz          || die "tabix index sidx_seqonly"

# ── sidx_named.tsv: column-name filter test file ──
#
# Has a '#' header line so column names can be used in -F/-O expressions.
# 6 columns (0-based):
#   0: CHROM   string
#   1: POS     int16
#   2: END     int16 (= POS, point records)
#   3: LABEL   string (A / B / C)
#   4: SCORE   int16  (42 / 100 / 42)
#   5: TYPE    string (protein_coding / lncRNA / rRNA)
printf '#CHROM\tPOS\tEND\tLABEL\tSCORE\tTYPE\n'     >  sidx_named.tsv
printf 'chr1\t1\t1\tA\t42\tprotein_coding\n'         >> sidx_named.tsv
printf 'chr1\t2\t2\tB\t100\tlncRNA\n'                >> sidx_named.tsv
printf 'chr1\t3\t3\tC\t42\trRNA\n'                   >> sidx_named.tsv

$BGZIP -c sidx_named.tsv > sidx_named.tsv.gz    || die "bgzip sidx_named"
$TABIX -f -s 1 -b 2 -e 3 -c '#' sidx_named.tsv.gz || die "tabix index sidx_named"

# ── Run dump_blocks once per file and capture TSV output ──
# (dump_blocks also creates the .sidx binary beside the .gz file)
$TABIX --dump-blocks sidx_data.tsv.gz  > sidx_dump.out       2>/dev/null \
    || die "dump_blocks sidx_data"
$TABIX --dump-blocks sidx_types.tsv.gz > sidx_types_dump.out 2>/dev/null \
    || die "dump_blocks sidx_types"
$TABIX --dump-blocks sidx_same.tsv.gz  > sidx_same_dump.out  2>/dev/null \
    || die "dump_blocks sidx_same"

echo ""

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: dump_blocks — .sidx creation and TSV format
# ─────────────────────────────────────────────────────────────────────────────
echo "--- Section 1: dump_blocks creates .sidx and emits TSV ---"

check ".sidx created for sidx_data (non-empty)" \
    test -s sidx_data.tsv.gz.sidx

# Verify binary magic "TBA\x02" at the start of the .sidx file.
printf 'TBA\002' > "${TMP}.magic_exp"
head -c 4 sidx_data.tsv.gz.sidx > "${TMP}.magic_got"
check ".sidx has correct magic bytes (TBA\\x02)" \
    cmp -s "${TMP}.magic_exp" "${TMP}.magic_got"

check "dump_blocks TSV: correct header line" \
    grep -q 'compressed_offset' sidx_dump.out

# File is ~368 KB; at 64 KB BGZF blocks we expect at least 4 block-summary rows.
# Block-summary rows start with a digit (the compressed offset).
check "dump_blocks TSV: >=4 BGZF block rows" \
    awk 'BEGIN{n=0}/^[0-9]/{n++}END{exit(n>=4?0:1)}' sidx_dump.out

# Per-column stat rows have a leading tab.
# CT_STRING columns show type "string" with "(N/A)" for both min and max.
check "dump_blocks TSV: col_0 (CHROM) annotated as string" \
    awk 'BEGIN{ok=0}/\tcol_0\tstring\t\(N\/A\)\t\(N\/A\)/{ok=1}END{exit(!ok)}' sidx_dump.out

# SCORE (col_4) has values 1000–14999; all fit in int16 (max 32767) → int16.
check "dump_blocks TSV: col_4 (SCORE) annotated as int16" \
    awk 'BEGIN{ok=0}/\tcol_4\tint16\t/{ok=1}END{exit(!ok)}' sidx_dump.out

# QUAL (col_5) is formatted with %.2f (≤7 significant digits) → float.
check "dump_blocks TSV: col_5 (QUAL) annotated as float" \
    awk 'BEGIN{ok=0}/\tcol_5\tfloat\t/{ok=1}END{exit(!ok)}' sidx_dump.out

check ".sidx created for sidx_types file" \
    test -s sidx_types.tsv.gz.sidx

# BVAL: single-char "0" or "1" → CT_BOOL → "boolean" in TSV
check "dump_blocks TSV: col_3 (BVAL) annotated as boolean" \
    awk 'BEGIN{ok=0}/\tcol_3\tboolean\t/{ok=1}END{exit(!ok)}' sidx_types_dump.out

# BIG: 100000 and 200000 exceed int16 max (32767) but fit in int32 → "int32"
check "dump_blocks TSV: col_4 (BIG) annotated as int32" \
    awk 'BEGIN{ok=0}/\tcol_4\tint32\t/{ok=1}END{exit(!ok)}' sidx_types_dump.out

# FP: 3.14 / 2.72 — strtod succeeds, count_sig_digits ≤ 7 → CT_FLOAT → "float"
check "dump_blocks TSV: col_5 (FP) annotated as float" \
    awk 'BEGIN{ok=0}/\tcol_5\tfloat\t/{ok=1}END{exit(!ok)}' sidx_types_dump.out

# DP: 3.14159265358979 — 15 significant digits → CT_DOUBLE → "double"
check "dump_blocks TSV: col_6 (DP) annotated as double" \
    awk 'BEGIN{ok=0}/\tcol_6\tdouble\t/{ok=1}END{exit(!ok)}' sidx_types_dump.out

# dump_blocks requires a .tbi index alongside the .gz; without one it returns 1.
$BGZIP -c sidx_data.tsv > sidx_notbi.tsv.gz || die "bgzip sidx_notbi"
check_fail "dump_blocks returns non-zero without .tbi" \
    "$TABIX" --dump-blocks sidx_notbi.tsv.gz

echo ""

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: query with -F
#   Exercises: sidx_load, tbx_filter_chunks, sidx_block_filter_cb,
#              tbx_apply_filters (row-level pass/fail).
# ─────────────────────────────────────────────────────────────────────────────
echo "--- Section 2: -F filter with .sidx (chunk + block pruning) ---"

# FILT_GE: value below all data — filter passes every row.
# Verifies sidx_load path (sidx is loaded, query proceeds normally).
check_out "FILT_GE passes all chr1 records (score>=1)" \
    sidx_exp_chr1_all.tsv \
    "$TABIX" -F "4>=1" sidx_data.tsv.gz chr1:1000-4999

# FILT_GE chunk pruning: chr1 SCORE max=4999 < 10000.
# All chr1 sidx blocks are impossible for >=10000 → tbx_filter_chunks drops
# every chunk before any BGZF block is read.
check_out "FILT_GE prunes all chr1 blocks (score>=10000)" \
    sidx_exp_empty.tsv \
    "$TABIX" -F "4>=10000" sidx_data.tsv.gz chr1:1000-4999

# FILT_LE chunk pruning: chr2 SCORE min=11000 > 4999.
# All chr2 sidx blocks are impossible for <=4999.
check_out "FILT_LE prunes all chr2 blocks (score<=4999)" \
    sidx_exp_empty.tsv \
    "$TABIX" -F "4<=4999" sidx_data.tsv.gz chr2:1000-4999

# FILT_EQ exact match: exactly one record has SCORE=2500 (chr1, pos=2500).
check_out "FILT_EQ selects exact record (score==2500)" \
    sidx_exp_score_eq2500.tsv \
    "$TABIX" -F "4==2500" sidx_data.tsv.gz chr1:1000-4999

# FILT_NE row-level: the only record in region chr1:4999-4999 has SCORE=4999.
# tbx_apply_filters sees 4999!=4999 = false → row is dropped → empty output.
check_out "FILT_NE drops exact match at position (score!=4999)" \
    sidx_exp_empty.tsv \
    "$TABIX" -F "4!=4999" sidx_data.tsv.gz chr1:4999-4999

# FILT_LE range: returns chr1 records with score<=1010 (positions 1000–1010).
check_out "FILT_LE range selection (score<=1010)" \
    sidx_exp_score_le1010.tsv \
    "$TABIX" -F "4<=1010" sidx_data.tsv.gz chr1:1000-4999

# FILT_GE range: returns chr2 records with score>=14990 (positions 4990–4999).
check_out "FILT_GE range selection (score>=14990)" \
    sidx_exp_score_ge14990.tsv \
    "$TABIX" -F "4>=14990" sidx_data.tsv.gz chr2:1000-4999

# Float column filter: chr1 records where qual>=49.0 (positions 4900–4999).
# Also prunes chr2 blocks (chr2 qual max=39.99 < 49.0).
check_out "FILT_GE on float column (qual>=49.0)" \
    sidx_exp_qual_ge49.tsv \
    "$TABIX" -F "5>=49.0" sidx_data.tsv.gz chr1:1000-4999

# FILT_NE block-level pruning: sidx_same.tsv.gz has all SCORE=42 in one block.
# dump_blocks recorded blk_min=42, blk_max=42.  The block-impossibility check
# for FILT_NE fires: (blk_min == val && blk_max == val) → drop chunk entirely.
check_out "FILT_NE block pruning: single-value SCORE block (score!=42)" \
    sidx_exp_empty.tsv \
    "$TABIX" -F "4!=42" sidx_same.tsv.gz chr1:1-3

echo ""

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: tbx_apply_filters edge cases
# ─────────────────────────────────────────────────────────────────────────────
echo "--- Section 3: tbx_apply_filters edge cases ---"

# Numeric filter applied to a string column (col_0 = CHROM):
#   Block level: col_min_tc == 0 (string) → skipped conservatively.
#   Row level: "chr1" is non-numeric → tbx_apply_filters prints an error
#              suggesting a string filter and returns 0 → all rows dropped.
check_out "numeric filter on string column drops all rows (col_0 >= 0)" \
    sidx_exp_empty.tsv \
    "$TABIX" -F "0>=0" sidx_data.tsv.gz chr1:1000-4999

# Out-of-bounds column (col_10 doesn't exist in the 6-column file):
#   Block level: f->col (10) >= sidx->ncols (6) → skipped conservatively.
#   Row level: tokeniser reaches end of line before column 10 → kept
#              conservatively (no data != failed filter).
check_out "OOB column filter keeps all rows conservatively (col_10 >= 0)" \
    sidx_exp_chr1_all.tsv \
    "$TABIX" -F "10>=0" sidx_data.tsv.gz chr1:1000-4999

echo ""

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: sidx_load error path — graceful degradation
# ─────────────────────────────────────────────────────────────────────────────
echo "--- Section 4: graceful degradation ---"

# No .sidx present: sidx_load can't open the file → sidx_ptr = NULL.
# Query falls back to row-level filtering only; results are identical.
rm -f sidx_data.tsv.gz.sidx
check_out "FILT_EQ correct without .sidx" \
    sidx_exp_score_eq2500.tsv \
    "$TABIX" -F "4==2500" sidx_data.tsv.gz chr1:1000-4999

# Corrupt .sidx (wrong magic "BAD!"): sidx_load reads magic, detects mismatch,
# returns -1 → sidx_ptr = NULL → graceful fallback; query still returns correct
# results via tbx_apply_filters.
printf 'BAD!' > sidx_data.tsv.gz.sidx
check_out "FILT_EQ correct with corrupt .sidx (bad magic)" \
    sidx_exp_score_eq2500.tsv \
    "$TABIX" -F "4==2500" sidx_data.tsv.gz chr1:1000-4999

echo ""

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: tbx_parse_filter error paths
# ─────────────────────────────────────────────────────────────────────────────
echo "--- Section 5: tbx_parse_filter error paths ---"

# Non-numeric column index: strtol("abc>=5") → end==expr → error message, exit 1.
check_fail "bad filter: non-numeric column index (abc>=5)" \
    "$TABIX" -F "abc>=5" sidx_data.tsv.gz chr1:1000-1001

# Unrecognised operator: after parsing col=4, "!!" is not in {==,!=,>=,<=} → error.
check_fail "bad filter: unrecognised operator (4!!5)" \
    "$TABIX" -F "4!!5" sidx_data.tsv.gz chr1:1000-1001

# Missing value: operator parsed, but *end == '\0' → error.
check_fail "bad filter: missing value after operator (4>=)" \
    "$TABIX" -F "4>=" sidx_data.tsv.gz chr1:1000-1001

# String value with ordered operator: strtod("abc") fails → is_string=1, but
# op is FILT_GE (>=) which is not permitted for strings → error.
check_fail "bad filter: string value with >= (4>=abc)" \
    "$TABIX" -F "4>=abc" sidx_data.tsv.gz chr1:1000-1001

# Same check for <= with a string value.
check_fail "bad filter: string value with <= (4<=word)" \
    "$TABIX" -F "4<=word" sidx_data.tsv.gz chr1:1000-1001

# Invalid regex pattern: regcomp returns error → exit 1.
check_fail "bad filter: invalid regex pattern (0~=[invalid)" \
    "$TABIX" -F "0~=[invalid" sidx_data.tsv.gz chr1:1000-1001

# Negative column index: col < 0 check in tbx_parse_filter → error.
check_fail "bad filter: negative column index (-1>=5)" \
    "$TABIX" -F "-1>=5" sidx_data.tsv.gz chr1:1000-1001

echo ""

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6: String filter support (== and != only)
#   Exercises: tbx_parse_filter string path, tbx_apply_filters string branch,
#              block-level functions skipping string filters (f->sval != NULL).
# ─────────────────────────────────────────────────────────────────────────────
echo "--- Section 6: string -F filters (== and !=) ---"

# Rebuild the sidx for sidx_data (removed in Section 4).
$TABIX --dump-blocks sidx_data.tsv.gz >/dev/null 2>/dev/null \
    || die "dump_blocks sidx_data for Section 6"

# FILT_EQ on CHROM (col 0): query is already restricted to chr1 by the TBI
# primary index, so 0==chr1 keeps every row in the result.
check_out "string FILT_EQ on CHROM passes all matching rows (0==chr1)" \
    sidx_exp_chr1_all.tsv \
    "$TABIX" -F "0==chr1" sidx_data.tsv.gz chr1:1000-4999

# FILT_NE on CHROM: within a chr1 query, 0!=chr1 rejects every row → empty.
check_out "string FILT_NE on CHROM rejects all rows in region (0!=chr1)" \
    sidx_exp_empty.tsv \
    "$TABIX" -F "0!=chr1" sidx_data.tsv.gz chr1:1000-4999

# FILT_EQ on NAME (col 3): only one record has NAME==row_2000.
check_out "string FILT_EQ on NAME selects single record (3==row_2000)" \
    sidx_exp_str_name_2000.tsv \
    "$TABIX" -F "3==row_2000" sidx_data.tsv.gz chr1:2000-2000

# sidx_same tests: three rows A/B/C, all SCORE=42.
# String FILT_EQ: keep only row B.
check_out "string FILT_EQ selects single row by name (3==B)" \
    sidx_exp_str_same_eqB.tsv \
    "$TABIX" -F "3==B" sidx_same.tsv.gz chr1:1-3

# String FILT_NE: reject row B, keep A and C.
check_out "string FILT_NE rejects named row, keeps others (3!=B)" \
    sidx_exp_str_same_neqB.tsv \
    "$TABIX" -F "3!=B" sidx_same.tsv.gz chr1:1-3

# Block-level functions must skip string filters (no numeric min/max in sidx).
# Combine a numeric filter that prunes blocks with a string filter that only
# acts at row level.  The numeric filter 4!=42 would drop the entire chunk
# (Section 2 test) but with 3==B added, the string filter must not interfere
# with the block pruning decision — the chunk is still dropped by 4!=42.
check_out "string filter does not block block-level pruning by numeric filter" \
    sidx_exp_empty.tsv \
    "$TABIX" -F "4!=42" -F "3==B" sidx_same.tsv.gz chr1:1-3

# Verify string filter alone does NOT trigger block-level pruning:
# 3==Z matches nothing, but the rows are checked at apply time (not pruned
# as a chunk) → result is empty, exit 0.
check_out "string FILT_EQ with no match returns empty output (3==Z)" \
    sidx_exp_empty.tsv \
    "$TABIX" -F "3==Z" sidx_same.tsv.gz chr1:1-3

echo ""

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 7: Regex filter support (~= and !~)
#   Exercises: tbx_parse_filter regex path (regcomp), tbx_apply_filters regex
#              branch (regexec), block-level functions skipping regex filters.
# ─────────────────────────────────────────────────────────────────────────────
echo "--- Section 7: regex -F filters (~= and !~) ---"

# Rebuild the sidx for sidx_data (may have been removed or corrupted in earlier sections).
$TABIX --dump-blocks sidx_data.tsv.gz >/dev/null 2>/dev/null \
    || die "dump_blocks sidx_data for Section 7"

# ~= anchored match: NAME (col 3) matches ^row_1[0-9][0-9][0-9]$ → rows 1000-1999.
check_out "regex FILT_RE anchored match selects rows 1000-1999 (3~=^row_1...)" \
    sidx_exp_re_row1xxx.tsv \
    "$TABIX" -F "3~=^row_1[0-9][0-9][0-9]\$" sidx_data.tsv.gz chr1:1000-4999

# !~ inverted match: NAME does NOT match ^row_1 → rows 2000-4999.
check_out "regex FILT_NRE inverted match selects rows 2000-4999 (3!~^row_1)" \
    sidx_exp_re_not_row1xxx.tsv \
    "$TABIX" -F "3!~^row_1" sidx_data.tsv.gz chr1:1000-4999

# ~= with alternation (ERE): CHROM matches chr1|chr2 — passes every row in the query.
check_out "regex FILT_RE ERE alternation passes all rows (0~=chr1|chr2)" \
    sidx_exp_chr1_all.tsv \
    "$TABIX" -F "0~=chr1|chr2" sidx_data.tsv.gz chr1:1000-4999

# ~= with pattern that matches nothing → empty output.
check_out "regex FILT_RE no-match pattern returns empty (3~=^nomatch)" \
    sidx_exp_empty.tsv \
    "$TABIX" -F "3~=^nomatch" sidx_data.tsv.gz chr1:1000-4999

# !~ with pattern that matches every row → empty output (every row excluded).
check_out "regex FILT_NRE universal match returns empty (3!~^row_)" \
    sidx_exp_empty.tsv \
    "$TABIX" -F "3!~^row_" sidx_data.tsv.gz chr1:1000-4999

# Regex filter must not interfere with numeric block-level pruning.
# 4!=42 drops the sidx_same chunk at block level; 3~=B is skipped at block level.
check_out "regex filter does not block block-level pruning by numeric filter" \
    sidx_exp_empty.tsv \
    "$TABIX" -F "4!=42" -F "3~=B" sidx_same.tsv.gz chr1:1-3

echo ""

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 8: Whole-file filter scan (no region argument)
#   When -F is given without any region, tabix enters the query path with the
#   "." wildcard, scanning the entire file and applying all filters.
# ─────────────────────────────────────────────────────────────────────────────
echo "--- Section 8: whole-file -F scan without region argument ---"

# Rebuild the sidx for sidx_data (may have been removed/corrupted earlier).
$TABIX --dump-blocks sidx_data.tsv.gz >/dev/null 2>/dev/null \
    || die "dump_blocks sidx_data for Section 8"

# Pre-build expected files that span both chromosomes.
awk '$5 == 2500' sidx_data.tsv > sidx_exp_whole_score_eq2500.tsv   # one chr1 row
awk '$5 >= 14990' sidx_data.tsv > sidx_exp_whole_score_ge14990.tsv # chr2 rows 4990-4999

# Numeric filter, no region: FILT_EQ on SCORE across the whole file.
# Block-level pruning (via sidx) should still fire for chr1 blocks.
check_out "whole-file FILT_EQ numeric filter (no region, score==2500)" \
    sidx_exp_whole_score_eq2500.tsv \
    "$TABIX" -F "4==2500" sidx_data.tsv.gz

# Numeric filter, no region: FILT_GE prunes chr1 blocks (max 4999 < 14990).
check_out "whole-file FILT_GE with block pruning (no region, score>=14990)" \
    sidx_exp_whole_score_ge14990.tsv \
    "$TABIX" -F "4>=14990" sidx_data.tsv.gz

# String filter, no region: exact match on CHROM returns only chr2 rows.
check_out "whole-file string FILT_EQ on CHROM (no region, 0==chr2)" \
    sidx_exp_whole_score_ge14990.tsv \
    "$TABIX" -F "0==chr2" -F "4>=14990" sidx_data.tsv.gz

# Regex filter, no region: NAME matches ^row_4999 — one row on each chromosome.
awk '$4=="row_4999"' sidx_data.tsv > sidx_exp_whole_re_row4999.tsv
check_out "whole-file regex filter (no region, 3~=^row_4999\$)" \
    sidx_exp_whole_re_row4999.tsv \
    "$TABIX" -F "3~=^row_4999\$" sidx_data.tsv.gz

echo ""

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 9: Whole-file -F scan on seqonly-indexed files
#
#   seqonly files index by gene/sequence name only (no position columns).
#   The "." all-records wildcard must NOT be expanded to ".:1-1" by the seqonly
#   expansion code — doing so looks for a literal sequence named "." and finds
#   nothing.  The fix: skip expansion when region == ".".
#
#   Also verifies that named-gene queries still work normally after the fix.
# ─────────────────────────────────────────────────────────────────────────────
echo "--- Section 9: whole-file -F scan on seqonly-indexed file ---"

# Named-gene query (sanity check — confirms the index works at all).
check_out "seqonly: named-gene query returns correct row (GeneC)" \
    <(grep '^GeneC' sidx_seqonly.tsv) \
    "$TABIX" sidx_seqonly.tsv.gz GeneC

# Whole-file string == filter: two genes are protein_coding.
grep 'protein_coding' sidx_seqonly.tsv > sidx_exp_seqonly_pc.tsv
check_out "seqonly: whole-file string FILT_EQ (4==protein_coding)" \
    sidx_exp_seqonly_pc.tsv \
    "$TABIX" -F "4==protein_coding" sidx_seqonly.tsv.gz

# Whole-file string != filter: excludes protein_coding, keeps lncRNA and rRNA.
grep -v 'protein_coding' sidx_seqonly.tsv > sidx_exp_seqonly_non_pc.tsv
check_out "seqonly: whole-file string FILT_NE (4!=protein_coding)" \
    sidx_exp_seqonly_non_pc.tsv \
    "$TABIX" -F "4!=protein_coding" sidx_seqonly.tsv.gz

# Whole-file regex filter: TYPE matches RNA (lncRNA and rRNA).
grep 'RNA' sidx_seqonly.tsv > sidx_exp_seqonly_rna.tsv
check_out "seqonly: whole-file regex FILT_RE (4~=RNA)" \
    sidx_exp_seqonly_rna.tsv \
    "$TABIX" -F "4~=RNA" sidx_seqonly.tsv.gz

# Whole-file regex !~ filter: TYPE does NOT match RNA → protein_coding rows.
check_out "seqonly: whole-file regex FILT_NRE (4!~RNA)" \
    sidx_exp_seqonly_pc.tsv \
    "$TABIX" -F "4!~RNA" sidx_seqonly.tsv.gz

# Filter that matches nothing → empty output, exit 0.
check_out "seqonly: whole-file filter with no match returns empty (4==miRNA)" \
    sidx_exp_empty.tsv \
    "$TABIX" -F "4==miRNA" sidx_seqonly.tsv.gz

echo ""

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 10: -O OR-filter support
#
#   -O filters are OR'd among themselves; the OR group and the -F AND group
#   are combined with AND:  row passes iff (all -F pass) AND (any -O passes).
#   OR filters must not participate in block-level chunk pruning.
# ─────────────────────────────────────────────────────────────────────────────
echo "--- Section 10: -O OR-filter expressions ---"

# Rebuild sidx for sidx_data (may have been invalidated by earlier sections).
$TABIX --dump-blocks sidx_data.tsv.gz >/dev/null 2>/dev/null \
    || die "dump_blocks sidx_data for Section 10"

# OR of two exact matches: rows with SCORE==2500 OR SCORE==3000.
awk '$5==2500 || $5==3000' sidx_data.tsv > sidx_exp_or_2500_3000.tsv
check_out "FILT_OR: two exact scores OR'd (4==2500 OR 4==3000)" \
    sidx_exp_or_2500_3000.tsv \
    "$TABIX" -O "4==2500" -O "4==3000" sidx_data.tsv.gz chr1:1000-4999

# OR combined with AND: chr2 only (AND) AND (SCORE==12500 OR SCORE==13000).
awk '$1=="chr2" && ($5==12500 || $5==13000)' sidx_data.tsv > sidx_exp_or_and.tsv
check_out "FILT_OR + FILT_AND: chr2 rows with score 12500 or 13000" \
    sidx_exp_or_and.tsv \
    "$TABIX" -F "0==chr2" -O "4==12500" -O "4==13000" sidx_data.tsv.gz chr2:1000-4999

# OR with regex: CHROM matches chr1 OR NAME matches row_2000.
# Within the chr1 query region, chr1 matches all rows; so all chr1 rows pass.
check_out "FILT_OR: regex OR exact string passes all chr1 rows" \
    sidx_exp_chr1_all.tsv \
    "$TABIX" -O "0~=chr1" -O "3==row_2000" sidx_data.tsv.gz chr1:1000-4999

# OR where neither branch matches → empty output.
check_out "FILT_OR: no match in either branch → empty output" \
    sidx_exp_empty.tsv \
    "$TABIX" -O "4==99999" -O "4==88888" sidx_data.tsv.gz chr1:1000-4999

# AND filter that fails overrides a passing OR filter → empty output.
# F: SCORE>=10000 (impossible on chr1, max=4999); O: NAME==row_1000 (passes).
# The AND group must pass first; since SCORE>=10000 fails, the row is rejected.
check_out "FILT_AND failure overrides passing FILT_OR" \
    sidx_exp_empty.tsv \
    "$TABIX" -F "4>=10000" -O "3==row_1000" sidx_data.tsv.gz chr1:1000-4999

# OR filters must not participate in block-level pruning.
#
# sidx_same has all SCORE=42.  -F "4!=42" drops the entire chunk at block level
# (blk_min==blk_max==val → impossible).  -O "4!=42" must NOT trigger the same
# pruning — doing so would drop the chunk before any row is examined, causing
# a second OR branch ("3==B") to never fire even though row B would pass it.
#
# Correct behaviour: chunk is kept → row B satisfies "3==B" → one row returned.
# Bug behaviour  : chunk pruned by "4!=42" OR branch → empty output.
check_out "FILT_OR: impossible OR branch does not prune chunk; passing branch fires" \
    sidx_exp_str_same_eqB.tsv \
    "$TABIX" -O "4!=42" -O "3==B" sidx_same.tsv.gz chr1:1-3

# Whole-file OR scan (no region), seqonly file: GeneA or GeneC — two rows.
awk '$1=="GeneA" || $1=="GeneC"' sidx_seqonly.tsv > sidx_exp_or_seqonly_AC.tsv
check_out "FILT_OR whole-file seqonly: two gene names OR'd (GeneA or GeneC)" \
    sidx_exp_or_seqonly_AC.tsv \
    "$TABIX" -O "0==GeneA" -O "0==GeneC" sidx_seqonly.tsv.gz

echo ""

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 11: Column-name filter support
#
#   -F / -O expressions can use a column name from the file's first header
#   line instead of a 0-based integer index.  The header line is the first
#   line of the BGZF file; a leading '#' is stripped before splitting by tab.
#   Names are resolved to 0-based column indexes at query time.
# ─────────────────────────────────────────────────────────────────────────────
echo "--- Section 11: column-name filters (-F with header line) ---"

# Expected output files (data rows only, no header).
grep -v '^#' sidx_named.tsv | awk '$4=="B"'        > sidx_exp_named_B.tsv
grep -v '^#' sidx_named.tsv | awk '$5>=50'          > sidx_exp_named_score50.tsv
grep -v '^#' sidx_named.tsv | grep 'RNA'            > sidx_exp_named_rna.tsv
grep -v '^#' sidx_named.tsv | awk '$6=="rRNA"'      > sidx_exp_named_rrna.tsv

# String filter by column name: LABEL==B → row B only.
check_out "column-name string FILT_EQ: LABEL==B" \
    sidx_exp_named_B.tsv \
    "$TABIX" -F "LABEL==B" sidx_named.tsv.gz chr1:1-3

# Numeric filter by column name: SCORE>=50 → row B (score=100).
check_out "column-name numeric FILT_GE: SCORE>=50" \
    sidx_exp_named_score50.tsv \
    "$TABIX" -F "SCORE>=50" sidx_named.tsv.gz chr1:1-3

# Regex filter by column name: TYPE~=RNA → rows B (lncRNA) and C (rRNA).
check_out "column-name regex FILT_RE: TYPE~=RNA" \
    sidx_exp_named_rna.tsv \
    "$TABIX" -F "TYPE~=RNA" sidx_named.tsv.gz chr1:1-3

# Numeric column index still works on a named-header file (backward compat).
check_out "column-index still works on named-header file (3==B)" \
    sidx_exp_named_B.tsv \
    "$TABIX" -F "3==B" sidx_named.tsv.gz chr1:1-3

# No match → empty output (not an error).
check_out "column-name filter no match returns empty (LABEL==Z)" \
    sidx_exp_empty.tsv \
    "$TABIX" -F "LABEL==Z" sidx_named.tsv.gz chr1:1-3

# Whole-file scan with column-name filter (no region argument).
check_out "column-name whole-file scan (TYPE==rRNA)" \
    sidx_exp_named_rrna.tsv \
    "$TABIX" -F "TYPE==rRNA" sidx_named.tsv.gz

# Unknown column name → non-zero exit.
check_fail "unknown column name exits with error (MISSINGCOL==X)" \
    "$TABIX" -F "MISSINGCOL==X" sidx_named.tsv.gz chr1:1-3

echo ""

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 12: Coordinate-based indexing (-s/-b/-e, no -p) with -F/-O filters
#
#   Regression guard for the "seqonly-is-default" change.  Providing -b or -e
#   must clear TBX_SEQONLY so that records are indexed at their real coordinates
#   and coordinate-range queries find them.  This section explicitly creates a
#   file, indexes it with -s -b -e (no -p flag), and verifies that -F and -O
#   filters work correctly alongside coordinate-range queries.
# ─────────────────────────────────────────────────────────────────────────────
echo "--- Section 12: coordinate-based (-s/-b/-e) indexing with -F/-O ---"

# Small coordinate-based test file (no comment lines, no header).
# 5 columns (0-based): CHROM  POS  END  NAME  SCORE
# chr1 rows: SCORE in {10, 50}; chr2 rows: SCORE in {75, 90}.
printf 'chr1\t100\t200\tA\t10\n' >  sidx_coords.tsv
printf 'chr1\t300\t400\tB\t50\n' >> sidx_coords.tsv
printf 'chr1\t500\t600\tC\t10\n' >> sidx_coords.tsv
printf 'chr2\t100\t200\tD\t75\n' >> sidx_coords.tsv
printf 'chr2\t300\t400\tE\t90\n' >> sidx_coords.tsv

$BGZIP -c sidx_coords.tsv > sidx_coords.tsv.gz  || die "bgzip sidx_coords"
# No -p flag: default preset is seqonly, but -b/-e clears TBX_SEQONLY so
# records are indexed at their actual POS/END coordinates.
$TABIX -f -s 1 -b 2 -e 3 sidx_coords.tsv.gz    || die "tabix index sidx_coords"

# ── Expected output files ──
awk '$4=="B"'                  sidx_coords.tsv > sidx_exp_coords_B.tsv
awk '$4=="A"'                  sidx_coords.tsv > sidx_exp_coords_A.tsv
awk '$4=="A" || $4=="C"'       sidx_coords.tsv > sidx_exp_coords_AC.tsv
awk '$1=="chr1" && $5>=50'     sidx_coords.tsv > sidx_exp_coords_chr1_ge50.tsv
awk '$1=="chr1" && $5==10'     sidx_coords.tsv > sidx_exp_coords_chr1_s10.tsv
awk '$5>=50'                   sidx_coords.tsv > sidx_exp_coords_ge50.tsv

# ── Sanity check: query that only resolves via real coordinates.
# In seqonly mode records sit at [0,1); chr1:300-400 (0-based [299,400)) does
# NOT overlap [0,1), so seqonly mode returns nothing.  Coordinate mode finds
# row B (POS=300, END=400) — confirming TBX_SEQONLY was cleared by -b/-e.
check_out "coord-indexed: coordinate-range chr1:300-400 finds row B" \
    sidx_exp_coords_B.tsv \
    "$TABIX" sidx_coords.tsv.gz chr1:300-400

# Numeric -F filter within a coordinate-range query.
# chr1:1-9999 returns all chr1 rows; SCORE>=50 keeps only row B.
check_out "coord-indexed: numeric FILT_GE with coord region (4>=50)" \
    sidx_exp_coords_chr1_ge50.tsv \
    "$TABIX" -F "4>=50" sidx_coords.tsv.gz chr1:1-9999

# String -F filter within a coordinate-range query.
check_out "coord-indexed: string FILT_EQ with coord region (3==A)" \
    sidx_exp_coords_A.tsv \
    "$TABIX" -F "3==A" sidx_coords.tsv.gz chr1:1-9999

# -O OR filter within a coordinate-range query.
check_out "coord-indexed: -O OR filter with coord region (3==A OR 3==C)" \
    sidx_exp_coords_AC.tsv \
    "$TABIX" -O "3==A" -O "3==C" sidx_coords.tsv.gz chr1:1-9999

# Combined -F AND + -O OR: CHROM==chr1 (AND) with SCORE==10 OR SCORE==90 (OR).
# chr1 rows: A(score=10) passes both; B(score=50) fails OR; C(score=10) passes.
check_out "coord-indexed: -F AND + -O OR with coord region" \
    sidx_exp_coords_chr1_s10.tsv \
    "$TABIX" -F "0==chr1" -O "4==10" -O "4==90" sidx_coords.tsv.gz chr1:1-9999

# Whole-file -F scan on coordinate-indexed file (no region argument).
check_out "coord-indexed: whole-file -F scan (4>=50)" \
    sidx_exp_coords_ge50.tsv \
    "$TABIX" -F "4>=50" sidx_coords.tsv.gz

echo ""

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 13: Cross-compatibility — .tbi from upstream (unmodified) tabix
#
#   A .tbi index created by any recent upstream tabix (using -s/-b/-e) must be
#   fully usable by our tabby with -F/-O filters.  Upstream tabix never sets
#   TBX_SEQONLY in the preset, so tabby should load the index, see no seqonly
#   flag, and apply coordinate-range queries + row-level filters correctly.
#   There will be no .sidx alongside the index, so tabby falls back to
#   row-level-only filtering — the results must still be correct.
#
#   This section runs only when an upstream-equivalent tabix binary is present
#   at ./tabix (the pre-rename binary left in the repo directory during
#   development).  In clean build environments it is silently skipped.
# ─────────────────────────────────────────────────────────────────────────────
echo "--- Section 13: cross-compat (.tbi from upstream tabix + tabby filters) ---"

UPSTREAM_TABIX="./tabix"
if [ ! -x "$UPSTREAM_TABIX" ]; then
    echo "  (skipped: $UPSTREAM_TABIX not present)"
else
    # Re-index sidx_coords.tsv.gz using the upstream tabix binary.
    # No .sidx is produced; tabby must fall back to row-level filtering.
    "$UPSTREAM_TABIX" -f -s 1 -b 2 -e 3 sidx_coords.tsv.gz \
        || die "upstream tabix indexing failed"

    # Sanity: coordinate-range query at chr1:300-400 must find row B.
    # In seqonly mode records sit at [0,1), so this range returns nothing —
    # a pass here proves TBX_SEQONLY is not set in the upstream-produced .tbi.
    check_out "cross-compat: upstream .tbi, coordinate-range query (no filter)" \
        sidx_exp_coords_B.tsv \
        "$TABIX" sidx_coords.tsv.gz chr1:300-400

    # Numeric -F filter on upstream-indexed file + coordinate region.
    check_out "cross-compat: upstream .tbi + -F FILT_GE (4>=50)" \
        sidx_exp_coords_chr1_ge50.tsv \
        "$TABIX" -F "4>=50" sidx_coords.tsv.gz chr1:1-9999

    # -O OR filter on upstream-indexed file + coordinate region.
    check_out "cross-compat: upstream .tbi + -O OR filter (3==A OR 3==C)" \
        sidx_exp_coords_AC.tsv \
        "$TABIX" -O "3==A" -O "3==C" sidx_coords.tsv.gz chr1:1-9999

    # Whole-file -F scan on upstream-indexed file (no region argument).
    check_out "cross-compat: upstream .tbi + whole-file -F scan (4>=50)" \
        sidx_exp_coords_ge50.tsv \
        "$TABIX" -F "4>=50" sidx_coords.tsv.gz
fi

echo ""

# ─────────────────────────────────────────────────────────────────────────────
# SUMMARY
# ─────────────────────────────────────────────────────────────────────────────
echo "Expected   passes:   $_np"
echo "Unexpected passes:   0"
echo "Expected   failures: 0"
echo "Unexpected failures: $_nf"
echo ""

[ "$_nf" -eq 0 ]
