/*
 * Copyright © 2019, Luca Barbato
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <stdlib.h>

#include "common/bitdepth.h"
#include "common/intops.h"

#include "src/cdef.h"
#include "src/cpu.h"

#include "src/mips/types.h"

#if BITDEPTH == 8
static inline v8i16 vconstrain(const v8i16 diff, const int16_t threshold,
                               const int damping)
{
    const v8i16 zero = __msa_fill_h(0);
    if (!threshold) return zero;
    const uint16_t shift = imax(0, damping - ulog2(threshold));
    const v8i16 abs_diff = __msa_add_a_h(diff, zero); /* Get abs */
    const v16u8 mask = __msa_clt_s_h(diff, zero);
    const v8i16 thr = __msa_fill_h(threshold);
    const v8i16 sub = __msa_subv_h(thr, __msa_sra_h(abs_diff, __msa_fill_h(shift)));
    const v8i16 max = __msa_max_s_h(zero, sub);
    const v8i16 min = __msa_min_s_h(abs_diff, max);
    const v8i16 neg = __msa_subv_h(zero, min);
    return __msa_bsel_v(min, neg, mask);
}

static inline void copy4xN(uint16_t *tmp, const ptrdiff_t tmp_stride,
                           const uint8_t *src, const ptrdiff_t src_stride,
                           const uint8_t (*left)[2], const uint8_t *const top,
                           const int w, const int h,
                           const enum CdefEdgeFlags edges)
{
    const v8u16 fill = __msa_fill_h((uint16_t)INT16_MAX);

    v8u16 l0;
    v8u16 l1;

    int y_start = -2, y_end = h + 2;

    // Copy top and bottom first
    if (!(edges & CDEF_HAVE_TOP)) {
        l0 = fill;
        l1 = fill;
        y_start = 0;
    } else {
        l0 = u8h_to_u16(__msa_ld_h((top + 0 * src_stride - 2), 0);
        l1 = u8h_to_u16(__msa_ld_h((top + 1 * src_stride - 2), 0);
    }

    __msa_st_h(l0, tmp - 2 * 8, 0);
    __msa_st_h(l1, tmp - 1 * 8, 0);

    if (!(edges & CDEF_HAVE_BOTTOM)) {
        l0 = fill;
        l1 = fill;
        y_end -= 2;
    } else {
        l0 = u8h_to_u16(__msa_ld_h((src - 2 + (h + 0) * src_stride), 0);
        l1 = u8h_to_u16(__msa_ld_h((src - 2 + (h + 1) * src_stride), 0);
    }

    __msa_st_h(l0, tmp + (h + 0) * 8, 0);
    __msa_st_h(l1, tmp + (h + 1) * 8, 0);

    for (int y = 0; y < h; y++) {
        v8u16 l = u8h_to_u16((0, src - 2 + y * src_stride), 0);
        __msa_st_h(l, tmp + y * 8, 0);
    }

    if (!(edges & CDEF_HAVE_LEFT)) {
        for (int y = y_start; y < y_end; y++) {
            tmp[y * 8] = INT16_MAX;
            tmp[1 + y * 8] = INT16_MAX;
        }
    } else {
        for (int y = 0; y < h; y++) {
            tmp[y * 8] = left[y][0];
            tmp[1 + y * 8] = left[y][1];
        }
    }
    if (!(edges & CDEF_HAVE_RIGHT)) {
        for (int y = y_start; y < y_end; y++) {
            tmp[- 2 + (y + 1) * 8] = INT16_MAX;
            tmp[- 1 + (y + 1) * 8] = INT16_MAX;
        }
    }
}

static inline void copy8xN(uint16_t *tmp, const ptrdiff_t tmp_stride,
                           const uint8_t *src, const ptrdiff_t src_stride,
                           const uint8_t (*left)[2], const uint8_t *const top,
                           const int w, const int h,
                           const enum CdefEdgeFlags edges)
{
    const v8u16 fill = __msa_fill_h((uint16_t)INT16_MAX);

    v8u16 l0h, l0l;
    v8u16 l1h, l1l;

    int y_start = -2, y_end = h + 2;

    // Copy top and bottom first
    if (!(edges & CDEF_HAVE_TOP)) {
        l0h = fill;
        l0l = fill;
        l1h = fill;
        l1l = fill;
        y_start = 0;
    } else {
        u8x16 l0 = __msa_ld_h(top + 0 * src_stride - 2, 0);
        u8x16 l1 = __msa_ld_h(top + 1 * src_stride - 2, 0);
        l0h = u8h_to_u16(l0);
        l0l = u8l_to_u16(l0);
        l1h = u8h_to_u16(l1);
        l1l = u8l_to_u16(l1);
    }

    __msa_st_h(l0h, tmp - 4 * 8, 0);
    __msa_st_h(l0l, tmp - 3 * 8, 0);
    __msa_st_h(l1h, tmp - 2 * 8, 0);
    __msa_st_h(l1l, tmp - 1 * 8, 0);

    if (!(edges & CDEF_HAVE_BOTTOM)) {
        l0h = fill;
        l0l = fill;
        l1h = fill;
        l1l = fill;
        y_end -= 2;
    } else {
        u8x16 l0 = __msa_ld_h(src - 2 + (h + 0) * src_stride, 0);
        u8x16 l1 = __msa_ld_h(src - 2 + (h + 1) * src_stride, 0);
        l0h = u8h_to_u16(l0);
        l0l = u8l_to_u16(l0);
        l1h = u8h_to_u16(l1);
        l1l = u8l_to_u16(l1);
    }

    __msa_st_h(l0h, tmp + (h + 0) * 16, 0);
    __msa_st_h(l0l, tmp + (h + 0) * 16 + 8, 0);
    __msa_st_h(l1h, tmp + (h + 1) * 16, 0);
    __msa_st_h(l1l, tmp + (h + 1) * 16 + 8, 0);

    for (int y = 0; y < h; y++) {
        u8x16 l = __msa_ld_h(src - 2 + y * src_stride, 0);
        v8u16 lh = u8h_to_u16(l);
        v8u16 ll = u8l_to_u16(l);
        __msa_st_h(lh, tmp + y * 16, 0);
        __msa_st_h(ll, tmp + 8 + y * 16, 0);
    }

    if (!(edges & CDEF_HAVE_LEFT)) {
        for (int y = y_start; y < y_end; y++) {
            tmp[y * 16] = INT16_MAX;
            tmp[1 + y * 16] = INT16_MAX;
        }
    } else {
        for (int y = 0; y < h; y++) {
            tmp[y * 16] = left[y][0];
            tmp[1 + y * 16] = left[y][1];
        }
    }
    if (!(edges & CDEF_HAVE_RIGHT)) {
        for (int y = y_start; y < y_end; y++) {
            tmp[- 6 + (y + 1) * 16] = INT16_MAX;
            tmp[- 5 + (y + 1) * 16] = INT16_MAX;
        }
    }
}

static inline v8i16 max_mask(v8i16 a, v8i16 b) {
    const v8i16 v8i16_INT16_MAX = __msa_fill_h((int16_t)INT16_MAX);

    const v16u8 mask = __msa_ceqi_h(a, v8i16_INT16_MAX);

    const v8i16 val = __msa_bsel_v(a, b, mask);

    return __msa_max_s_h(val, b);
}

#define LOAD_PIX(addr) \
    const v8i16 px = (v8i16)__msa_ld_h(addr, 0); \
    v8i16 max = px; \
    v8i16 min = px; \
    v8i16 sum = __msa_fill_h(0);

#define LOAD_PIX4(addr) \
    const v8i16 a = (v8i16)__msa_ld_h(addr, 0); \
    const v8i16 b = (v8i16)__msa_ld_h(addr + tmp_stride, 0); \
    const v8i16 px = __msa_ilvr_d(a, b, 0); \
    v8i16 max = px; \
    v8i16 min = px; \
    v8i16 sum = __msa_fill_h(0);

#define LOAD_DIR(p, addr, o0, o1) \
    const v8i16 p ## 0 = (v8i16)__msa_ld_h(addr + o0, 0); \
    const v8i16 p ## 1 = (v8i16)__msa_ld_h(addr - o0, 0); \
    const v8i16 p ## 2 = (v8i16)__msa_ld_h(addr + o1, 0); \
    const v8i16 p ## 3 = (v8i16)__msa_ld_h(addr - o1, 0);

#define LOAD_DIR4(p, addr, o0, o1) \
    LOAD_DIR(p ## a, addr, o0, o1) \
    LOAD_DIR(p ## b, addr + tmp_stride, o0, o1) \
    const v8i16 p ## 0 = __msa_ilvr_d(p ## a ## 0, p ## b ## 0); \
    const v8i16 p ## 1 = __msa_ilvr_d(p ## a ## 1, p ## b ## 1); \
    const v8i16 p ## 2 = __msa_ilvr_d(p ## a ## 2, p ## b ## 2); \
    const v8i16 p ## 3 = __msa_ilvr_d(p ## a ## 3, p ## b ## 3);

#define CONSTRAIN(p, strength) \
    const v8i16 p ## _d0 = __msa_subv_h(p ## 0, px); \
    const v8i16 p ## _d1 = __msa_subv_h(p ## 1, px); \
    const v8i16 p ## _d2 = __msa_subv_h(p ## 2, px); \
    const v8i16 p ## _d3 = __msa_subv_h(p ## 3, px); \
\
    v8i16 p ## _c0 = vconstrain(p ## _d0, strength, damping); \
    v8i16 p ## _c1 = vconstrain(p ## _d1, strength, damping); \
    v8i16 p ## _c2 = vconstrain(p ## _d2, strength, damping); \
    v8i16 p ## _c3 = vconstrain(p ## _d3, strength, damping);

#define MIN_MAX(p) \
    max = max_mask(p ## 0, max); \
    min = __msa_min_s_h(p ## 0, min); \
    max = max_mask(p ## 1, max); \
    min = __msa_min_s_h(p ## 1, min); \
    max = max_mask(p ## 2, max); \
    min = __msa_min_s_h(p ## 2, min); \
    max = max_mask(p ## 3, max); \
    min = __msa_min_s_h(p ## 3, min);

#define PRI_0(p) \
    p ## _c0 = __msa_addv_h(__msa_slli_h(p ## _c0, 1), __msa_sll_h(p ## _c0, __msa_fill_h(tap_even))); \
    p ## _c1 = __msa_addv_h(__msa_slli_h(p ## _c1, 1), __msa_sll_h(p ## _c1, __msa_fill_h(tap_even)));

#define PRI_1(p) \
    p ## _c2 = __msa_subv_h(__msa_slli_h(p ## _c2, 2), __msa_sll_h(p ## _c2, __msa_fill_h(tap_even))); \
    p ## _c3 = __msa_subv_h(__msa_slli_h(p ## _c3, 2)), __msa_sll_h((p ## _c3, __msa_fill_h(tap_even)));

#define SEC_0(p) \
    p ## _c0 = __msa_slli_h(p ## _c0, 1); \
    p ## _c1 = __msa_slli_h(p ## _c1, 1); \
    p ## _c2 = __msa_slli_h(p ## _c2, 1); \
    p ## _c3 = __msa_slli_h(p ## _c3, 1);

#define UPDATE_SUM(p) \
    const v8i16 p ## sum0 = __msa_addv_h(p ## _c0, p ## _c1); \
    const v8i16 p ## sum1 = __msa_addv_h(p ## _c2, p ## _c3); \
    sum = __msa_addv_h(sum, p ## sum0); \
    sum = __msa_addv_h(sum, p ## sum1);

static inline void
filter_4xN(pixel *dst, const ptrdiff_t dst_stride,
           const pixel (*left)[2], const pixel *const top,
           const int w, const int h, const int pri_strength,
           const int sec_strength, const int dir,
           const int damping, const enum CdefEdgeFlags edges,
           const ptrdiff_t tmp_stride, uint16_t *tmp)
{
    const int8_t cdef_directions[8 /* dir */][2 /* pass */] = {
        { -1 * tmp_stride + 1, -2 * tmp_stride + 2 },
        {  0 * tmp_stride + 1, -1 * tmp_stride + 2 },
        {  0 * tmp_stride + 1,  0 * tmp_stride + 2 },
        {  0 * tmp_stride + 1,  1 * tmp_stride + 2 },
        {  1 * tmp_stride + 1,  2 * tmp_stride + 2 },
        {  1 * tmp_stride + 0,  2 * tmp_stride + 1 },
        {  1 * tmp_stride + 0,  2 * tmp_stride + 0 },
        {  1 * tmp_stride + 0,  2 * tmp_stride - 1 }
    };
    const int bitdepth_min_8 = bitdepth_from_max(bitdepth_max) - 8;
    const uint16_t tap_even = !((pri_strength >> bitdepth_min_8) & 1);
    const int off1 = cdef_directions[dir][0];
    const int off1_1 = cdef_directions[dir][1];

    const int off2 = cdef_directions[(dir + 2) & 7][0];
    const int off3 = cdef_directions[(dir + 6) & 7][0];

    const int off2_1 = cdef_directions[(dir + 2) & 7][1];
    const int off3_1 = cdef_directions[(dir + 6) & 7][1];


    copy4xN(tmp - 2, tmp_stride, dst, dst_stride, left, top, w, h, edges);
    for (int y = 0; y < h / 2; y++) {
        LOAD_PIX4(tmp)

        // Primary pass
        LOAD_DIR4(p, tmp, off1, off1_1)

        CONSTRAIN(p, pri_strength)

        MIN_MAX(p)

        (p)
        PRI_1(p)

        UPDATE_SUM(p)

        // Secondary pass 1
        LOAD_DIR4(s, tmp, off2, off3)

        CONSTRAIN(s, sec_strength)

        MIN_MAX(s)

        SEC_0(s)

        UPDATE_SUM(s)

        // Secondary pass 2
        LOAD_DIR4(s2, tmp, off2_1, off3_1)

        CONSTRAIN(s2, sec_strength)

        MIN_MAX(s2)

        UPDATE_SUM(s2)

        // Store
        v8i16 bias = __msa_and_v(__msa_clt_s_h(sum, __msa_fill_h(0)), __msa_fill_h(1));
        bias = __msa_subv_h(__msa_fill_h(8), bias);
        v8i16 unclamped = __msa_addv_h(px, __msa_sra_h(__msa_addv_h(sum, bias), __msa_fill_h(4)));
        v8i16 vdst = __msa_max_s_h(__msa_min_s_h(unclamped, max), min);

        dst[0] = vdst[0];
        dst[1] = vdst[1];
        dst[2] = vdst[2];
        dst[3] = vdst[3];

        tmp += tmp_stride;
        dst += PXSTRIDE(dst_stride);
        dst[0] = vdst[4];
        dst[1] = vdst[5];
        dst[2] = vdst[6];
        dst[3] = vdst[7];

        tmp += tmp_stride;
        dst += PXSTRIDE(dst_stride);
    }
}

static inline void
filter_8xN(pixel *dst, const ptrdiff_t dst_stride,
           const pixel (*left)[2], const pixel *const top,
           const int w, const int h, const int pri_strength,
           const int sec_strength, const int dir,
           const int damping, const enum CdefEdgeFlags edges,
           const ptrdiff_t tmp_stride, uint16_t *tmp)
{
    const int8_t cdef_directions[8 /* dir */][2 /* pass */] = {
        { -1 * tmp_stride + 1, -2 * tmp_stride + 2 },
        {  0 * tmp_stride + 1, -1 * tmp_stride + 2 },
        {  0 * tmp_stride + 1,  0 * tmp_stride + 2 },
        {  0 * tmp_stride + 1,  1 * tmp_stride + 2 },
        {  1 * tmp_stride + 1,  2 * tmp_stride + 2 },
        {  1 * tmp_stride + 0,  2 * tmp_stride + 1 },
        {  1 * tmp_stride + 0,  2 * tmp_stride + 0 },
        {  1 * tmp_stride + 0,  2 * tmp_stride - 1 }
    };
    const int bitdepth_min_8 = bitdepth_from_max(bitdepth_max) - 8;


    const uint16_t tap_even = !((pri_strength >> bitdepth_min_8) & 1);
    const int off1 = cdef_directions[dir][0];
    const int off1_1 = cdef_directions[dir][1];

    const int off2 = cdef_directions[(dir + 2) & 7][0];
    const int off3 = cdef_directions[(dir + 6) & 7][0];

    const int off2_1 = cdef_directions[(dir + 2) & 7][1];
    const int off3_1 = cdef_directions[(dir + 6) & 7][1];

    copy8xN(tmp - 2, tmp_stride, dst, dst_stride, left, top, w, h, edges);

    for (int y = 0; y < h; y++) {
        LOAD_PIX(tmp)

        // Primary pass
        LOAD_DIR(p, tmp, off1, off1_1)

        CONSTRAIN(p, pri_strength)

        MIN_MAX(p)

        PRI_0(p)
        PRI_1(p)

        UPDATE_SUM(p)

        // Secondary pass 1
        LOAD_DIR(s, tmp, off2, off3)

        CONSTRAIN(s, sec_strength)

        MIN_MAX(s)

        SEC_0(s)

        UPDATE_SUM(s)

        // Secondary pass 2
        LOAD_DIR(s2, tmp, off2_1, off3_1)

        CONSTRAIN(s2, sec_strength)

        MIN_MAX(s2)

        UPDATE_SUM(s2)

        // Store
        v8i16 bias = __msa_and_v(__msa_clt_s_h(sum, __msa_fill_h(0)), __msa_fill_h(1));
        bias = __msa_subv_h(__msa_fill_h(8), bias);
        v8i16 unclamped = __msa_addv_h(px, __msa_sra_h(__msa_addv_h(sum, bias), __msa_fill_h(4)));
        v8i16 vdst = __msa_max_s_h(__msa_min_s_h(unclamped, max), min);

        dst[0] = vdst[0];
        dst[1] = vdst[1];
        dst[2] = vdst[2];
        dst[3] = vdst[3];
        dst[4] = vdst[4];
        dst[5] = vdst[5];
        dst[6] = vdst[6];
        dst[7] = vdst[7];

        tmp += tmp_stride;
        dst += PXSTRIDE(dst_stride);
    }

}


#define cdef_fn(w, h, tmp_stride) \
static void cdef_filter_##w##x##h##_msa(pixel *const dst, \
                                        const ptrdiff_t dst_stride, \
                                        const pixel (*left)[2], \
                                        const pixel *const top, \
                                        const int pri_strength, \
                                        const int sec_strength, \
                                        const int dir, \
                                        const int damping, \
                                        const enum CdefEdgeFlags edges) \
{ \
    ALIGN_STK_16(uint16_t, tmp_buf, 12 * tmp_stride,); \
    uint16_t *tmp = tmp_buf + 2 * tmp_stride + 2; \
    filter_##w##xN(dst, dst_stride, left, top, w, h, pri_strength, sec_strength, \
                   dir, damping, edges, tmp_stride, tmp); \
}

cdef_fn(4, 4, 8);
cdef_fn(4, 8, 8);
cdef_fn(8, 8, 16);
#endif

COLD void bitfn(dav1d_cdef_dsp_init_mips)(Dav1dCdefDSPContext *const c) {
    const unsigned flags = dav1d_get_cpu_flags();

    if (!(flags & DAV1D_MIPS_CPU_FLAG_MSA)) return;

#if BITDEPTH == 8
    // c->dir = dav1d_cdef_find_dir_msa;
    c->fb[0] = cdef_filter_8x8_msa;
    c->fb[1] = cdef_filter_4x8_msa;
    c->fb[2] = cdef_filter_4x4_msa;
#endif
}
