#!/bin/bash
# test-sidx.sh -- Integration tests for the tabix secondary-index (sidx).
#
# Functions exercised:
#   dump_blocks()          -- creates .sidx and emits per-block TSV
#   sidx_load()            -- loads .sidx on query startup; graceful on failure
#   sidx_free()            -- called automatically after every query
#   tbx_filter_chunks()    -- chunk-level pruning with all four filter ops
#   sidx_block_filter_cb() -- block-level callback for local-file queries
#   tbx_apply_filters()    -- per-row exact check; string/OOB column handling
#   tbx_parse_filter()     -- -F expression parsing; all error paths
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

TABIX="../../tabix"
BGZIP="../../bgzip"
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
    rm -f sidx_notbi.tsv.gz
    rm -f sidx_dump.out sidx_types_dump.out sidx_same_dump.out
    rm -f sidx_exp_*.tsv
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

# String column (col_0 = CHROM):
#   Block level: col_min_tc == 0 (string) → skipped conservatively (block kept).
#   Row level: "chr1" is non-numeric → tbx_apply_filters prints error and
#              returns 0 → every row is dropped → empty output.
check_out "string column filter drops all rows (col_0 >= 0)" \
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

# Non-numeric value: strtod("abc") → vend==end → error.
check_fail "bad filter: non-numeric value (4>=abc)" \
    "$TABIX" -F "4>=abc" sidx_data.tsv.gz chr1:1000-1001

# Negative column index: col < 0 check in tbx_parse_filter → error.
check_fail "bad filter: negative column index (-1>=5)" \
    "$TABIX" -F "-1>=5" sidx_data.tsv.gz chr1:1000-1001

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
