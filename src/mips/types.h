/*
 * Copyright © 2019, VideoLAN and dav1d authors
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

#ifndef DAV1D_SRC_MIPS_TYPES_H
#define DAV1D_SRC_MIPS_TYPES_H

#include <msa.h>
#undef pixel

#define u8h_to_u16(v) ((u16x8) __msa_ilvev_b((v16u8) v, __msa_fill_b(0)))
#define u8l_to_u16(v) ((u16x8) __msa_ilvod_b((v16u8) v, __msa_fill_b(0)))
#define u16h_to_i32(v) ((v4i32) __msa_ilvev_h((u16x8) v, __msa_fill_h(0)))
//#define i16h_to_i32(v) ((v4i32) vec_unpackh((v8i16)v))
#define u16l_to_i32(v) ((v4i32) __msa_ilvod_h((u16x8) v, __msa_fill_h(0)))
//#define i16l_to_i32(v) ((v4i32) vec_unpackl((v8i16)v))

#endif /* DAV1D_SRC_PPC_TYPES_H */
