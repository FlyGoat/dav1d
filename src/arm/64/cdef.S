/*
 * Copyright © 2018, VideoLAN and dav1d authors
 * Copyright © 2019, Martin Storsjo
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

#include "src/arm/asm.S"
#include "util.S"
#include "cdef_tmpl.S"

.macro pad_top_bottom s1, s2, w, stride, rn, rw, ret
        tst             w6,  #1 // CDEF_HAVE_LEFT
        b.eq            2f
        // CDEF_HAVE_LEFT
        sub             \s1,  \s1,  #2
        sub             \s2,  \s2,  #2
        tst             w6,  #2 // CDEF_HAVE_RIGHT
        b.eq            1f
        // CDEF_HAVE_LEFT+CDEF_HAVE_RIGHT
        ldr             \rn\()0, [\s1]
        ldr             s1,      [\s1, #\w]
        ldr             \rn\()2, [\s2]
        ldr             s3,      [\s2, #\w]
        uxtl            v0.8h,   v0.8b
        uxtl            v1.8h,   v1.8b
        uxtl            v2.8h,   v2.8b
        uxtl            v3.8h,   v3.8b
        str             \rw\()0, [x0]
        str             d1,      [x0, #2*\w]
        add             x0,  x0,  #2*\stride
        str             \rw\()2, [x0]
        str             d3,      [x0, #2*\w]
.if \ret
        ret
.else
        add             x0,  x0,  #2*\stride
        b               3f
.endif

1:
        // CDEF_HAVE_LEFT+!CDEF_HAVE_RIGHT
        ldr             \rn\()0, [\s1]
        ldr             h1,      [\s1, #\w]
        ldr             \rn\()2, [\s2]
        ldr             h3,      [\s2, #\w]
        uxtl            v0.8h,   v0.8b
        uxtl            v1.8h,   v1.8b
        uxtl            v2.8h,   v2.8b
        uxtl            v3.8h,   v3.8b
        str             \rw\()0, [x0]
        str             s1,      [x0, #2*\w]
        str             s31,     [x0, #2*\w+4]
        add             x0,  x0,  #2*\stride
        str             \rw\()2, [x0]
        str             s3,      [x0, #2*\w]
        str             s31,     [x0, #2*\w+4]
.if \ret
        ret
.else
        add             x0,  x0,  #2*\stride
        b               3f
.endif

2:
        // !CDEF_HAVE_LEFT
        tst             w6,  #2 // CDEF_HAVE_RIGHT
        b.eq            1f
        // !CDEF_HAVE_LEFT+CDEF_HAVE_RIGHT
        ldr             \rn\()0, [\s1]
        ldr             h1,      [\s1, #\w]
        ldr             \rn\()2, [\s2]
        ldr             h3,      [\s2, #\w]
        uxtl            v0.8h,  v0.8b
        uxtl            v1.8h,  v1.8b
        uxtl            v2.8h,  v2.8b
        uxtl            v3.8h,  v3.8b
        str             s31, [x0]
        stur            \rw\()0, [x0, #4]
        str             s1,      [x0, #4+2*\w]
        add             x0,  x0,  #2*\stride
        str             s31, [x0]
        stur            \rw\()2, [x0, #4]
        str             s3,      [x0, #4+2*\w]
.if \ret
        ret
.else
        add             x0,  x0,  #2*\stride
        b               3f
.endif

1:
        // !CDEF_HAVE_LEFT+!CDEF_HAVE_RIGHT
        ldr             \rn\()0, [\s1]
        ldr             \rn\()1, [\s2]
        uxtl            v0.8h,  v0.8b
        uxtl            v1.8h,  v1.8b
        str             s31,     [x0]
        stur            \rw\()0, [x0, #4]
        str             s31,     [x0, #4+2*\w]
        add             x0,  x0,  #2*\stride
        str             s31,     [x0]
        stur            \rw\()1, [x0, #4]
        str             s31,     [x0, #4+2*\w]
.if \ret
        ret
.else
        add             x0,  x0,  #2*\stride
.endif
3:
.endm

.macro load_n_incr dst, src, incr, w
.if \w == 4
        ld1             {\dst\().s}[0], [\src], \incr
.else
        ld1             {\dst\().8b},   [\src], \incr
.endif
.endm

// void dav1d_cdef_paddingX_8bpc_neon(uint16_t *tmp, const pixel *src,
//                                    ptrdiff_t src_stride, const pixel (*left)[2],
//                                    const pixel *const top, int h,
//                                    enum CdefEdgeFlags edges);

.macro padding_func w, stride, rn, rw
function cdef_padding\w\()_8bpc_neon, export=1
        movi            v30.8h,  #0x80, lsl #8
        mov             v31.16b, v30.16b
        sub             x0,  x0,  #2*(2*\stride+2)
        tst             w6,  #4 // CDEF_HAVE_TOP
        b.ne            1f
        // !CDEF_HAVE_TOP
        st1             {v30.8h, v31.8h}, [x0], #32
.if \w == 8
        st1             {v30.8h, v31.8h}, [x0], #32
.endif
        b               3f
1:
        // CDEF_HAVE_TOP
        add             x9,  x4,  x2
        pad_top_bottom  x4,  x9, \w, \stride, \rn, \rw, 0

        // Middle section
3:
        tst             w6,  #1 // CDEF_HAVE_LEFT
        b.eq            2f
        // CDEF_HAVE_LEFT
        tst             w6,  #2 // CDEF_HAVE_RIGHT
        b.eq            1f
        // CDEF_HAVE_LEFT+CDEF_HAVE_RIGHT
0:
        ld1             {v0.h}[0], [x3], #2
        ldr             h2,      [x1, #\w]
        load_n_incr     v1,  x1,  x2,  \w
        subs            w5,  w5,  #1
        uxtl            v0.8h,  v0.8b
        uxtl            v1.8h,  v1.8b
        uxtl            v2.8h,  v2.8b
        str             s0,      [x0]
        stur            \rw\()1, [x0, #4]
        str             s2,      [x0, #4+2*\w]
        add             x0,  x0,  #2*\stride
        b.gt            0b
        b               3f
1:
        // CDEF_HAVE_LEFT+!CDEF_HAVE_RIGHT
        ld1             {v0.h}[0], [x3], #2
        load_n_incr     v1,  x1,  x2,  \w
        subs            w5,  w5,  #1
        uxtl            v0.8h,  v0.8b
        uxtl            v1.8h,  v1.8b
        str             s0,      [x0]
        stur            \rw\()1, [x0, #4]
        str             s31,     [x0, #4+2*\w]
        add             x0,  x0,  #2*\stride
        b.gt            1b
        b               3f
2:
        tst             w6,  #2 // CDEF_HAVE_RIGHT
        b.eq            1f
        // !CDEF_HAVE_LEFT+CDEF_HAVE_RIGHT
0:
        ldr             h1,      [x1, #\w]
        load_n_incr     v0,  x1,  x2,  \w
        subs            w5,  w5,  #1
        uxtl            v0.8h,  v0.8b
        uxtl            v1.8h,  v1.8b
        str             s31,     [x0]
        stur            \rw\()0, [x0, #4]
        str             s1,      [x0, #4+2*\w]
        add             x0,  x0,  #2*\stride
        b.gt            0b
        b               3f
1:
        // !CDEF_HAVE_LEFT+!CDEF_HAVE_RIGHT
        load_n_incr     v0,  x1,  x2,  \w
        subs            w5,  w5,  #1
        uxtl            v0.8h,  v0.8b
        str             s31,     [x0]
        stur            \rw\()0, [x0, #4]
        str             s31,     [x0, #4+2*\w]
        add             x0,  x0,  #2*\stride
        b.gt            1b

3:
        tst             w6,  #8 // CDEF_HAVE_BOTTOM
        b.ne            1f
        // !CDEF_HAVE_BOTTOM
        st1             {v30.8h, v31.8h}, [x0], #32
.if \w == 8
        st1             {v30.8h, v31.8h}, [x0], #32
.endif
        ret
1:
        // CDEF_HAVE_BOTTOM
        add             x9,  x1,  x2
        pad_top_bottom  x1,  x9, \w, \stride, \rn, \rw, 1
endfunc
.endm

padding_func 8, 16, d, q
padding_func 4, 8,  s, d

tables

filter 8, 8
filter 4, 8

find_dir 8
