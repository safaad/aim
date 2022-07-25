/*
 *                             The MIT License
 *
 * Wavefront Alignments Algorithms
 * Copyright (c) 2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 * This file is part of Wavefront Alignments Algorithms.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * PROJECT: Wavefront Alignments Algorithms
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */
/* MIT License

Copyright (c) 2021 SAFARI Research Group at ETH Zürich

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE. */
#include "../common/common.h"
#include "wfa_backtracing.h"

#include "dpu_allocator_wram.h"
#include <barrier.h>

void edit_cigar_allocate(
    edit_cigar_t *edit_cigar,
    int pattern_length,
    int text_length,
    dpu_alloc_wram_t *allocator)
{
    edit_cigar->max_operations = pattern_length + text_length;
    edit_cigar->begin_offset = edit_cigar->max_operations - 1;
    edit_cigar->end_offset = edit_cigar->max_operations;
    edit_cigar->score = INT32_MIN;
}

void affine_wfa_reduce_wvs(wfa_component *wfa, awf_offset_t pattern_length, awf_offset_t text_length, int score)
{
    int min_wavefront_length = 10;
    int max_distance_threshold = 50;
    int alignment_k = AFFINE_WAVEFRONT_DIAGONAL(text_length, pattern_length);

    if (wfa == NULL || wfa->m_null)
        return;
    if ((wfa->khi - wfa->klo + 1) < min_wavefront_length)
        return;

    int min_distance = MAX(pattern_length, text_length);

    int klo = wfa->klo;
    int khi = wfa->khi;

    for (int k = klo; k <= khi; ++k)
    {
        awf_offset_t offset = wfa->mwavefront[k];
        int v = AFFINE_WAVEFRONT_V(k, offset);
        int h = AFFINE_WAVEFRONT_H(k, offset);
        int left_v = pattern_length - v;
        int left_h = text_length - h;
        int distance = MAX(left_v, left_h);
        min_distance = MIN(distance, min_distance);
    }

    // reduce m
    // reduce from bottom
    int top_limit = MIN(alignment_k - 1, khi);

    for (int k = wfa->klo; k < top_limit; ++k)
    {
        awf_offset_t offset = wfa->mwavefront[k];

        int v = AFFINE_WAVEFRONT_V(k, offset);
        int h = AFFINE_WAVEFRONT_H(k, offset);
        int left_v = pattern_length - v;
        int left_h = text_length - h;
        int distance = MAX(left_v, left_h);
        if ((distance - min_distance) <= max_distance_threshold)
            break;
        wfa->klo = wfa->klo + 1;
    }

    // reduce from top
    int bottom_limit = MAX(alignment_k + 1, wfa->klo);
    for (int k = khi; k > bottom_limit; --k)
    {
        awf_offset_t offset = wfa->mwavefront[k];

        int v = AFFINE_WAVEFRONT_V(k, offset);
        int h = AFFINE_WAVEFRONT_H(k, offset);
        int left_v = pattern_length - v;
        int left_h = text_length - h;
        int distance = MAX(left_v, left_h);

        if (distance - min_distance <= max_distance_threshold)
            break;
        wfa->khi = wfa->khi - 1;
    }

    if (wfa->klo > wfa->khi)
    {
        wfa->m_null = true;
        wfa->i_null = true;
        wfa->d_null = true;
        wfa->khi = khi;
        wfa->klo = klo;
        return;
    }
}

// insert new score
wfa_component *allocate_new_score(dpu_alloc_wram_t *allocator, int score, int lo, int hi, int kernel)
{

    int wv_len = hi - lo + 1;

    wfa_component *wfa_cmpnt = (wfa_component *)allocate_new(allocator, sizeof(wfa_component));
    awf_offset_t *offset_ptr = (awf_offset_t *)allocate_new(allocator, (wv_len * sizeof(awf_offset_t)));

    wfa_cmpnt->mwavefront = (awf_offset_t *)(offset_ptr - lo);
    if (kernel == 3 || kernel == 1)
    {

        awf_offset_t *offset_ptr = (awf_offset_t *)allocate_new(allocator, (wv_len * sizeof(awf_offset_t)));
        wfa_cmpnt->dwavefront = offset_ptr - lo;
        wfa_cmpnt->d_null = false;
    }
    else
    {
        wfa_cmpnt->dwavefront = NULL;
        wfa_cmpnt->d_null = true;
    }
    if (kernel == 3 || kernel == 2)
    {
        awf_offset_t *offset_ptr = (awf_offset_t *)allocate_new(allocator, (wv_len * sizeof(awf_offset_t)));
        wfa_cmpnt->iwavefront = offset_ptr - lo;
        wfa_cmpnt->i_null = false;
    }
    else
    {
        wfa_cmpnt->iwavefront = NULL;
        wfa_cmpnt->i_null = true;
    }
    wfa_cmpnt->m_null = false;

    wfa_cmpnt->klo = lo;
    wfa_cmpnt->khi = hi;
    wfa_cmpnt->lo_base = lo;
    wfa_cmpnt->hi_base = hi;

    return wfa_cmpnt;
}

// wavefront extend matching
void affine_wfa_extend(wfa_component *wfa, char *pattern, char *text, awf_offset_t pattern_len, awf_offset_t text_len, int score)
{
    if (wfa == NULL || wfa->m_null)
        return;

    int klo = wfa->klo;
    for (int k = klo; k <= wfa->khi; ++k)
    {
        int moffset = wfa->mwavefront[k];
        if (moffset < 0)
            continue;

        int v = moffset - k;
        int h = moffset;

        int count = 0;
        while ((v < pattern_len && h < text_len && v >= 0 && h >= 0) && pattern[v++] == text[h++])
        {
            ++count;
        }
        wfa->mwavefront[k] += count;
    }
}
// end reached
bool affine_wfa_end_reached(wfa_component *wfa, awf_offset_t pattern_len, awf_offset_t text_len, int score)
{

    if (wfa == NULL || wfa->m_null)
        return false;

    int alignment_k =
        AFFINE_WAVEFRONT_DIAGONAL(text_len, pattern_len);
    int alignment_offset =
        AFFINE_WAVEFRONT_OFFSET(text_len, pattern_len);

    if (wfa->klo <= alignment_k && wfa->khi >= alignment_k)
    {
        int offset = wfa->mwavefront[alignment_k];

        if (offset >= alignment_offset)
            return true;
    }

    return false;
}
void affine_wfa_compute_offsets(wfa_component *wfa, wfa_set wfa_set, int lo, int hi, int score, int kernel)
{

    for (int k = lo; k <= hi; ++k)
    {
        awf_offset_t ins = -10;
        if (!wfa_set.m_o_null || !wfa_set.i_e_null)
        {
            // Update I
            awf_offset_t ins_g = AFFINE_WAVEFRONT_COND_FETCH(wfa_set.m_o_null, wfa_set.m_o_lo, wfa_set.m_o_hi, k - 1, wfa_set.wfa_o->mwavefront[(k)-1]);
            awf_offset_t ins_i = AFFINE_WAVEFRONT_COND_FETCH(wfa_set.i_e_null, wfa_set.e_lo, wfa_set.e_hi, k - 1, wfa_set.wfa_e->iwavefront[(k)-1]);
            if (ins_g == AFFINE_WAVEFRONT_OFFSET_NULL && ins_i == AFFINE_WAVEFRONT_OFFSET_NULL)
                ins = AFFINE_WAVEFRONT_OFFSET_NULL;
            else
                ins = MAX(ins_g, ins_i) + 1;
            wfa->iwavefront[k] = ins;
        }
        awf_offset_t del = -10;
        if (!wfa_set.m_o_null || !wfa_set.d_e_null)
        {
            // Update D
            awf_offset_t del_g = AFFINE_WAVEFRONT_COND_FETCH(wfa_set.m_o_null, wfa_set.m_o_lo, wfa_set.m_o_hi, k + 1, wfa_set.wfa_o->mwavefront[(k) + 1]);
            awf_offset_t del_d = AFFINE_WAVEFRONT_COND_FETCH(wfa_set.d_e_null, wfa_set.e_lo, wfa_set.e_hi, k + 1, wfa_set.wfa_e->dwavefront[(k) + 1]);

            del = MAX(del_g, del_d);
            wfa->dwavefront[k] = del;
        }
        // Update M
        awf_offset_t sub = -10;
        if (!wfa_set.m_sub_null)
            sub = AFFINE_WAVEFRONT_COND_FETCH(wfa_set.m_sub_null, wfa_set.m_sub_lo, wfa_set.m_sub_hi, k, wfa_set.wfa_sub->mwavefront[k] + 1);

        awf_offset_t new = MAX(sub, ins);
        wfa->mwavefront[k] = MAX(del, new);
    }
}

void affine_wfa_compute_next(wfa_component **wavefronts, dpu_alloc_wram_t *alloc_obj, int score)
{
    wfa_set wfa_set;

    // get previous scores
    int mismatch_score = score - MISMATCH;
    int o_score = score - GAP_O - GAP_E;
    int e_score = score - GAP_E;

    // is null?
    wfa_set.m_sub_null = ((mismatch_score < 0) || (wavefronts[mismatch_score] == NULL) || (wavefronts[mismatch_score]->m_null));
    wfa_set.m_o_null = ((o_score < 0) || (wavefronts[o_score] == NULL) || (wavefronts[o_score]->m_null));
    wfa_set.i_e_null = ((e_score < 0) || (wavefronts[e_score] == NULL || wavefronts[e_score]->iwavefront == NULL || wavefronts[e_score]->i_null));
    wfa_set.d_e_null = ((e_score < 0) || (wavefronts[e_score] == NULL || wavefronts[e_score]->dwavefront == NULL || wavefronts[e_score]->d_null));

    wfa_set.i_out_null = wfa_set.m_o_null && wfa_set.i_e_null;
    wfa_set.d_out_null = wfa_set.m_o_null && wfa_set.d_e_null;

    // null mwavefront
    if (wfa_set.m_sub_null && (wfa_set.i_out_null && wfa_set.d_out_null))
    {
        wavefronts[score] = NULL;
        return;
    }

    // compute limits
    if (wfa_set.m_sub_null)
    {
        wfa_set.m_sub_lo = 1;
        wfa_set.m_sub_hi = -1;
    }
    else
    {
        wfa_set.m_sub_lo = wavefronts[mismatch_score]->klo;
        wfa_set.m_sub_hi = wavefronts[mismatch_score]->khi;
        wfa_set.wfa_sub = wavefronts[mismatch_score];
    }

    if (wfa_set.m_o_null)
    {
        wfa_set.m_o_lo = 1;
        wfa_set.m_o_hi = -1;
    }
    else
    {
        wfa_set.m_o_lo = wavefronts[o_score]->klo;
        wfa_set.m_o_hi = wavefronts[o_score]->khi;
        wfa_set.wfa_o = wavefronts[o_score];
    }

    if (wfa_set.i_e_null && wfa_set.d_e_null)
    {
        wfa_set.e_lo = 1;
        wfa_set.e_hi = -1;
    }
    else
    {
        wfa_set.e_lo = wavefronts[e_score]->klo;
        wfa_set.e_hi = wavefronts[e_score]->khi;
        wfa_set.wfa_e = wavefronts[e_score];
    }
    int lo = MIN(wfa_set.m_sub_lo, wfa_set.m_o_lo);
    lo = MIN(lo, wfa_set.e_lo) - 1;
    int hi = MAX(wfa_set.m_sub_hi, wfa_set.m_o_hi);
    hi = MAX(hi, wfa_set.e_hi) + 1;

    // Compute WFA
    int kernel = ((!wfa_set.i_out_null) << 1) | (!wfa_set.d_out_null);

    wavefronts[score] = allocate_new_score(alloc_obj, score, lo, hi, kernel);

    affine_wfa_compute_offsets(wavefronts[score], wfa_set, lo, hi, score, kernel);
}

void affine_wfa_compute(dpu_alloc_wram_t *dpu_alloc_wram, edit_cigar_t *cigar, char pattern[], char text[], int pattern_length, int text_length)
{

    wfa_component *wavefronts[MAX_SCORE + 1] = {NULL};

    wavefronts[0] = allocate_new_score(dpu_alloc_wram, 0, 0, 0, 0);
    wavefronts[0]->mwavefront[0] = 0;

    int score = 0;
    while (true)
    {

        affine_wfa_extend(wavefronts[score], pattern, text, pattern_length, text_length, score);

#ifdef REDUCE
        affine_wfa_reduce_wvs(wavefronts[score], pattern_length, text_length, score);
#endif
        if (affine_wfa_end_reached(wavefronts[score], pattern_length, text_length, score))
        {
#ifdef BACKTRACE
            affine_wavefronts_backtrace(wavefronts, cigar, pattern, pattern_length, text, text_length, score);
#endif
            cigar->score = score;
            return;
        }

        ++score;
        if (score > MAX_SCORE)
        {
            cigar->score = score;
#ifdef BACKTRACE
            affine_wavefronts_backtrace(wavefronts, cigar, pattern, pattern_length, text, text_length, score);
#endif
            return;
        }
        affine_wfa_compute_next(wavefronts, dpu_alloc_wram, score);
    }
}

int main()
{
    mem_reset();
    uint32_t tasklet_id = me();
    dpu_alloc_wram_t dpu_alloc_wram;

    // Load parameters
    uint32_t params_m = (uint32_t)DPU_MRAM_HEAP_POINTER;
    DPUParams params_w;
    mram_read((__mram_ptr void const *)params_m, &params_w, ROUND_UP_MULTIPLE_8(sizeof(DPUParams)));
    uint32_t nb_reads_per_dpu = params_w.dpuNumReads;

    if (nb_reads_per_dpu <= 0)
        return 0;

    int nb_reads_per_tasklets = (((nb_reads_per_dpu + NR_TASKLETS) / NR_TASKLETS));
    // Each tasklet allocates WRAM segment
    dpu_alloc_wram = init_dpu_alloc_wram(WRAM_SEGMENT);

    uint32_t dpuRequests_m = ((uint32_t)DPU_MRAM_HEAP_POINTER) + params_w.dpuRequests_m;
    uint32_t dpuResults_m = ((uint32_t)DPU_MRAM_HEAP_POINTER) + params_w.dpuResults_m;
    uint32_t dpuPatterns_m = ((uint32_t)DPU_MRAM_HEAP_POINTER) + params_w.dpuPatterns_m;
    uint32_t dpuTexts_m = ((uint32_t)DPU_MRAM_HEAP_POINTER) + params_w.dpuTexts_m;
    uint32_t dpuOperations_m = ((uint32_t)DPU_MRAM_HEAP_POINTER) + params_w.dpuOperations_m;

    for (int read_nb = 0; read_nb < nb_reads_per_tasklets; ++read_nb)
    {
        request_t *request_w = (request_t *)allocate_new(&dpu_alloc_wram, (sizeof(request_t)));
        result_t *result_w = (result_t *)allocate_new(&dpu_alloc_wram, (sizeof(result_t)));

        edit_cigar_t *cigar;
        cigar = (edit_cigar_t *)allocate_new(&dpu_alloc_wram, sizeof(edit_cigar_t));

        if (read_nb + tasklet_id * nb_reads_per_tasklets < nb_reads_per_dpu)
        {
            mram_read((__mram_ptr void const *)(dpuRequests_m + (read_nb + tasklet_id * nb_reads_per_tasklets) * (sizeof(request_t))), request_w, ROUND_UP_MULTIPLE_8(sizeof(request_t)));
            char *pattern = (char *)allocate_new(&dpu_alloc_wram, ROUND_UP_MULTIPLE_8(request_w->pattern_len));
            char *text = (char *)allocate_new(&dpu_alloc_wram, ROUND_UP_MULTIPLE_8(request_w->text_len));

            // DMA transfers must be less than 2048
            if (ROUND_UP_MULTIPLE_8(request_w->pattern_len) <= 2048)
            {
                mram_read((__mram_ptr void const *)(dpuPatterns_m + (read_nb + tasklet_id * nb_reads_per_tasklets) * (READ_SIZE)), pattern, ROUND_UP_MULTIPLE_8(request_w->pattern_len));
            }
            else
            {
                for (int segment_size = 0; segment_size <= ROUND_UP_MULTIPLE_8(request_w->pattern_len); segment_size += 2048)
                {
                    if (segment_size + 2048 <= ROUND_UP_MULTIPLE_8(request_w->pattern_len))
                    {
                        mram_read((__mram_ptr void const *)(dpuPatterns_m + (read_nb + tasklet_id * nb_reads_per_tasklets) * (READ_SIZE) + segment_size), &(pattern[segment_size]), 2048);
                    }
                    else
                    {
                        int size = ROUND_UP_MULTIPLE_8(request_w->pattern_len) - segment_size;
                        mram_read((__mram_ptr void const *)(dpuPatterns_m + (read_nb + tasklet_id * nb_reads_per_tasklets) * (READ_SIZE) + segment_size), &(pattern[segment_size]), ROUND_UP_MULTIPLE_8(size));
                    }
                }
            }

            if (ROUND_UP_MULTIPLE_8(request_w->text_len) <= 2048)
            {
                mram_read((__mram_ptr void const *)(dpuTexts_m + (read_nb + tasklet_id * nb_reads_per_tasklets) * (READ_SIZE)), text, ROUND_UP_MULTIPLE_8(request_w->text_len));
            }
            else
            {
                for (int segment_size = 0; segment_size <= ROUND_UP_MULTIPLE_8(request_w->text_len); segment_size += 2048)
                {
                    if (segment_size + 2048 <= ROUND_UP_MULTIPLE_8(request_w->text_len))
                    {
                        mram_read((__mram_ptr void const *)(dpuTexts_m + (read_nb + tasklet_id * nb_reads_per_tasklets) * (READ_SIZE) + segment_size), &(text[segment_size]), 2048);
                    }
                    else
                    {
                        int size = ROUND_UP_MULTIPLE_8(request_w->text_len) - segment_size;
                        mram_read((__mram_ptr void const *)(dpuTexts_m + (read_nb + tasklet_id * nb_reads_per_tasklets) * (READ_SIZE) + segment_size), &(text[segment_size]), ROUND_UP_MULTIPLE_8(size));
                    }
                }
            }
            edit_cigar_allocate(cigar, request_w->pattern_len, request_w->text_len, &dpu_alloc_wram);

#ifdef BACKTRACE
            cigar->operations = (char *)allocate_new(&dpu_alloc_wram, 2 * READ_SIZE);
            // initialize traceback operations segment
            memset(cigar->operations, 'M', 2 * READ_SIZE);
#endif
            affine_wfa_compute(&dpu_alloc_wram, cigar, pattern, text, request_w->pattern_len, request_w->text_len);

            result_w->idx = request_w->idx;

#ifdef BACKTRACE
            if (ROUND_UP_MULTIPLE_8(cigar->max_operations) <= 2048)
            {
                mram_write((cigar->operations), (__mram_ptr void *)(dpuOperations_m + (read_nb + tasklet_id * nb_reads_per_tasklets) * (2 * READ_SIZE)), ROUND_UP_MULTIPLE_8(cigar->max_operations));
            }
            else
            {
                for (int segment_size = 0; segment_size <= ROUND_UP_MULTIPLE_8(cigar->max_operations); segment_size += 2048)
                {
                    if (segment_size + 2048 <= ROUND_UP_MULTIPLE_8(cigar->max_operations))
                    {
                        mram_write(&(cigar->operations[segment_size]), (__mram_ptr void *)(dpuOperations_m + (read_nb + tasklet_id * nb_reads_per_tasklets) * (2 * READ_SIZE) + segment_size), 2048);
                    }
                    else
                    {
                        int size = ROUND_UP_MULTIPLE_8(cigar->max_operations) - segment_size;
                        mram_write(&(cigar->operations[segment_size]), (__mram_ptr void *)(dpuOperations_m + (read_nb + tasklet_id * nb_reads_per_tasklets) * (2 * READ_SIZE) + segment_size), ROUND_UP_MULTIPLE_8(size));
                    }
                }
            }
#endif
            result_w->score = cigar->score;
            result_w->max_operations = cigar->max_operations;
            result_w->begin_offset = cigar->begin_offset;
            result_w->end_offset = cigar->end_offset;

            mram_write(result_w, (__mram_ptr void *)(dpuResults_m + (read_nb + tasklet_id * nb_reads_per_tasklets) * (sizeof(result_t))), sizeof(result_t));
        }
        // reset WRAM segment for each tasklet after every alignment
        reset_dpu_alloc_wram(&dpu_alloc_wram);
    }
    return 0;
}