/// @file htslib/tbx_cw.h

#ifndef HTSLIB_TBX_CW_H
#define HTSLIB_TBX_CW_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int64_t beg, end;
    char *ss, *se;
    int tid;
} tbx_intv_t;

#ifdef __cplusplus
}
#endif

#endif
