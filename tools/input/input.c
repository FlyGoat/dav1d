/*
 * Copyright © 2018, VideoLAN and dav1d authors
 * Copyright © 2018, Two Orioles, LLC
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

#include "config.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/attributes.h"
#include "common/intops.h"

#include "input/input.h"
#include "input/demuxer.h"

struct DemuxerContext {
    DemuxerPriv *data;
    const Demuxer *impl;
};

#define MAX_NUM_DEMUXERS 3
static const Demuxer *demuxers[MAX_NUM_DEMUXERS];
static int num_demuxers = 0;

#define register_demuxer(impl) { \
    extern const Demuxer impl; \
    assert(num_demuxers < MAX_NUM_DEMUXERS); \
    demuxers[num_demuxers++] = &impl; \
}

void init_demuxers(void) {
    register_demuxer(ivf_demuxer);
    register_demuxer(annexb_demuxer);
    register_demuxer(section5_demuxer);
}

int input_open(DemuxerContext **const c_out,
               const char *const name, const char *const filename,
               unsigned fps[2], unsigned *const num_frames, unsigned timebase[2])
{
    const Demuxer *impl;
    DemuxerContext *c;
    int res, i;

    if (name) {
        for (i = 0; i < num_demuxers; i++) {
            if (!strcmp(demuxers[i]->name, name)) {
                impl = demuxers[i];
                break;
            }
        }
        if (i == num_demuxers) {
            fprintf(stderr, "Failed to find demuxer named \"%s\"\n", name);
            return DAV1D_ERR(ENOPROTOOPT);
        }
    } else {
        int probe_sz = 0;
        for (i = 0; i < num_demuxers; i++)
            probe_sz = imax(probe_sz, demuxers[i]->probe_sz);
        uint8_t *const probe_data = malloc(probe_sz);
        if (!probe_data) {
            fprintf(stderr, "Failed to allocate memory\n");
            return DAV1D_ERR(ENOMEM);
        }
        FILE *f = fopen(filename, "rb");
        res = !!fread(probe_data, 1, probe_sz, f);
        fclose(f);
        if (!res) {
            free(probe_data);
            fprintf(stderr, "Failed to read probe data\n");
            return errno ? DAV1D_ERR(errno) : DAV1D_ERR(EIO);
        }

        for (i = 0; i < num_demuxers; i++) {
            if (demuxers[i]->probe(probe_data)) {
                impl = demuxers[i];
                break;
            }
        }
        free(probe_data);
        if (i == num_demuxers) {
            fprintf(stderr,
                    "Failed to probe demuxer for file %s\n",
                    filename);
            return DAV1D_ERR(ENOPROTOOPT);
        }
    }

    if (!(c = malloc(sizeof(DemuxerContext) + impl->priv_data_size))) {
        fprintf(stderr, "Failed to allocate memory\n");
        return DAV1D_ERR(ENOMEM);
    }
    memset(c, 0, sizeof(DemuxerContext) + impl->priv_data_size);
    c->impl = impl;
    c->data = (DemuxerPriv *) &c[1];
    if ((res = impl->open(c->data, filename, fps, num_frames, timebase)) < 0) {
        free(c);
        return res;
    }
    *c_out = c;

    return 0;
}

int input_read(DemuxerContext *const ctx, Dav1dData *const data) {
    return ctx->impl->read(ctx->data, data);
}

void input_close(DemuxerContext *const ctx) {
    ctx->impl->close(ctx->data);
    free(ctx);
}
