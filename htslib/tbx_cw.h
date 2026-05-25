/// @file htslib/tbx_cw.h

#ifndef HTSLIB_TBX_CW_H
#define HTSLIB_TBX_CW_H

#ifdef __cplusplus
extern "C" {
#endif

/* ── Secondary index in-memory representation ───────────────────────────────
 *
 * sidx_t holds the complete contents of a .sidx file loaded from disk.
 * It maps 1-to-1 onto the binary layout written by dump_blocks():
 *
 *   Field 0  magic[4]                    – verified on load; not stored
 *   Field 1  n_columns        → ncols
 *   Field 2a column_min_types → col_min_tc[ncols]
 *   Field 2b column_max_types → col_max_tc[ncols]
 *   Field 3  n_blocks         → nblocks
 *   Field 3b max_comp_block_size → max_comp_block_size
 *   Field 4  block_record_length → rec_len
 *   Field 5  offsets[]        → blk_rec_offsets[nblocks]
 *   Field 6  block_records    → records[nblocks * rec_len bytes]
 *   Field 7  compressed_offsets[] → blk_comp_offsets[nblocks]
 *
 * To filter a candidate block during query_regions():
 *
 *   uint64_t blk_off = (uint64_t)itr->curr_off >> 16;
 *   int64_t  bidx    = sidx_block_index_from_offset(blk_off,
 *                          sidx->blk_comp_offsets, sidx->nblocks);
 *   if (bidx >= 0) {
 *       const uint8_t *rec = sidx->records + sidx->blk_rec_offsets[bidx];
 *       // Walk rec[] column by column using sidx->col_min_tc / col_max_tc
 *       // and sidx_type_size() to read per-column min/max values and decide
 *       // whether this block can possibly satisfy the query.
 *   }
 */
typedef struct {
    uint32_t  ncols;               /* number of tab-separated columns          */
    uint8_t  *col_min_tc;          /* sidx type code for each column's min     */
    uint8_t  *col_max_tc;          /* sidx type code for each column's max     */
    uint64_t  nblocks;             /* number of BGZF blocks recorded           */
    uint16_t  max_comp_block_size; /* largest compressed block size (bytes)    */
    uint64_t  rec_len;             /* byte length of one block record          */
    uint64_t *blk_rec_offsets;     /* blk_rec_offsets[i]: byte offset of block
                                    * i's record within records[]              */
    uint8_t  *records;             /* nblocks * rec_len bytes of block records */
    uint64_t *blk_comp_offsets;    /* blk_comp_offsets[i]: BGZF compressed file
                                    * offset of block i (upper 48 bits of the
                                    * virtual offset, i.e. voff >> 16); sorted
                                    * ascending; used by
                                    * sidx_block_index_from_offset()           */
} sidx_t;


typedef struct {
    int64_t beg, end;
    char *ss, *se;
    int tid;
} tbx_intv_t;

int tbx_parse1(const tbx_conf_t *conf, size_t len, char *line, tbx_intv_t *intv);
static int64_t sidx_block_index_from_offset(uint64_t blk_off, const uint64_t *sidx_blk_offsets, uint64_t n_blocks);
static size_t sidx_type_size(uint8_t tc);

#ifdef __cplusplus
}
#endif

#endif
