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
#include "dpu_allocator_mram.h"

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

    // //reduce from top
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
// change this function so it allocates using score % needed
// EDIT THIS FUNCTION!!!!
wfa_component *allocate_new_score(dpu_alloc_wram_t *allocator, int score, int lo, int hi, int kernel, uint32_t *mramIdx, dpu_alloc_mram_t *dpu_alloc_mram)
{
    int wv_len = hi - lo + 1;
    uint32_t cmpnt_size = 0;

// these two lines show that the size of the wfa_cmpnt is constant, but the size of the wfa cmpnt's offsets is not constant.
    wfa_component *wfa_cmpnt = (wfa_component *)allocate_new(allocator, sizeof(wfa_component));
    awf_offset_t *offset_ptr = (awf_offset_t *)allocate_new(allocator, (wv_len * sizeof(awf_offset_t)));

    wfa_cmpnt->mwavefront = (awf_offset_t *)(offset_ptr - lo);
    cmpnt_size += ROUND_UP_MULTIPLE_8(wv_len * sizeof(awf_offset_t));
    // kernel = 3 means I and D
    // kernel = 2 means I
    // kernel = 1 means D
    if (kernel == 3 || kernel == 1)
    {

        awf_offset_t *offset_ptr = (awf_offset_t *)allocate_new(allocator, (wv_len * sizeof(awf_offset_t)));

        wfa_cmpnt->dwavefront = offset_ptr - lo;
        wfa_cmpnt->d_null = false;
        cmpnt_size += ROUND_UP_MULTIPLE_8(wv_len * sizeof(awf_offset_t));
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
        cmpnt_size += ROUND_UP_MULTIPLE_8(wv_len * sizeof(awf_offset_t));
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
    cmpnt_size += ROUND_UP_MULTIPLE_8(sizeof(wfa_component));
    add_wfa_cmpnt_to_mram(mramIdx, cmpnt_size, dpu_alloc_mram);
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
    // Compute score wavefronts (core)
    for (int k = lo; k <= hi; ++k)
    {
        awf_offset_t ins = -10;
        if (!wfa_set.m_o_null || !wfa_set.i_e_null)
        {
            // Update I
            awf_offset_t ins_g = AFFINE_WAVEFRONT_COND_FETCH(wfa_set.m_o_null, wfa_set.m_o_lo, wfa_set.m_o_hi, k - 1, wfa_set.wfa_o_mwavefront[(k)-1]);
            awf_offset_t ins_i = AFFINE_WAVEFRONT_COND_FETCH(wfa_set.i_e_null, wfa_set.e_lo, wfa_set.e_hi, k - 1, wfa_set.wfa_e_iwavefront[(k)-1]);
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
            awf_offset_t del_g = AFFINE_WAVEFRONT_COND_FETCH(wfa_set.m_o_null, wfa_set.m_o_lo, wfa_set.m_o_hi, k + 1, wfa_set.wfa_o_mwavefront[(k) + 1]);
            awf_offset_t del_d = AFFINE_WAVEFRONT_COND_FETCH(wfa_set.d_e_null, wfa_set.e_lo, wfa_set.e_hi, k + 1, wfa_set.wfa_e_dwavefront[(k) + 1]);

            del = MAX(del_g, del_d);
            wfa->dwavefront[k] = del;
        }
        // Update M
        awf_offset_t sub = -10;
        if (!wfa_set.m_sub_null)
            sub = AFFINE_WAVEFRONT_COND_FETCH(wfa_set.m_sub_null, wfa_set.m_sub_lo, wfa_set.m_sub_hi, k, wfa_set.wfa_sub_mwavefront[k] + 1);

        awf_offset_t new = MAX(sub, ins);
        wfa->mwavefront[k] = MAX(del, new);
    }
}

wfa_component *affine_wfa_compute_next(int score, uint32_t *mramIdx, dpu_alloc_wram_t *alloc_obj, dpu_alloc_mram_t *dpu_alloc_mram)
{
    wfa_set wfa_set;

    // get previous scores
    int mismatch_score = score - MISMATCH;
    int o_score = score - GAP_O - GAP_E;
    int e_score = score - GAP_E;

// edit here so it computes differently and loads mwaefront from score % something
    wfa_component *wfa_mismatch = (mismatch_score < 0 || mramIdx[mismatch_score] == 0) ? NULL : load_mwavefront_cmpnt_from_mram(alloc_obj, mramIdx[mismatch_score]);
    wfa_component *wfa_o_score = (o_score < 0 || mramIdx[o_score] == 0) ? NULL : load_mwavefront_cmpnt_from_mram(alloc_obj, mramIdx[o_score]);
    wfa_component *wfa_e_score = (e_score < 0 || mramIdx[e_score] == 0) ? NULL : load_idwavefront_cmpnt_from_mram(alloc_obj, mramIdx[e_score]);

    // is null?
    wfa_set.m_sub_null = ((mismatch_score < 0) || (wfa_mismatch == NULL) || (wfa_mismatch->m_null));
    wfa_set.m_o_null = ((o_score < 0) || (wfa_o_score == NULL) || (wfa_o_score->m_null));
    wfa_set.i_e_null = ((e_score < 0) || (wfa_e_score == NULL || wfa_e_score->i_null || wfa_e_score->iwavefront == NULL));
    wfa_set.d_e_null = ((e_score < 0) || (wfa_e_score == NULL || wfa_e_score->d_null || wfa_e_score->dwavefront == NULL));

    wfa_set.i_out_null = wfa_set.m_o_null && wfa_set.i_e_null;
    wfa_set.d_out_null = wfa_set.m_o_null && wfa_set.d_e_null;

    // null mwavefront
    if (wfa_set.m_sub_null && (wfa_set.i_out_null && wfa_set.d_out_null))
    {
        //  if the wavefront is null store 0 in the mram idx
        mramIdx[score] = 0;
        return NULL;
    }

    // compute limits
    if (wfa_set.m_sub_null)
    {
        wfa_set.m_sub_lo = 1;
        wfa_set.m_sub_hi = -1;
    }
    else
    {
        wfa_set.m_sub_lo = wfa_mismatch->klo;
        wfa_set.m_sub_hi = wfa_mismatch->khi;
        wfa_set.wfa_sub_mwavefront = wfa_mismatch->mwavefront;
    }

    if (wfa_set.m_o_null)
    {
        wfa_set.m_o_lo = 1;
        wfa_set.m_o_hi = -1;
    }
    else
    {
        wfa_set.m_o_lo = wfa_o_score->klo;
        wfa_set.m_o_hi = wfa_o_score->khi;
        wfa_set.wfa_o_mwavefront = wfa_o_score->mwavefront;
    }

    if (wfa_set.i_e_null && wfa_set.d_e_null)
    {
        wfa_set.e_lo = 1;
        wfa_set.e_hi = -1;
    }
    else
    {
        wfa_set.e_lo = wfa_e_score->klo;
        wfa_set.e_hi = wfa_e_score->khi;
        wfa_set.wfa_e_iwavefront = wfa_e_score->iwavefront;
        wfa_set.wfa_e_dwavefront = wfa_e_score->dwavefront;
    }
    int lo = MIN(wfa_set.m_sub_lo, wfa_set.m_o_lo);
    lo = MIN(lo, wfa_set.e_lo) - 1;
    int hi = MAX(wfa_set.m_sub_hi, wfa_set.m_o_hi);
    hi = MAX(hi, wfa_set.e_hi) + 1;

    // Compute WF
    int kernel = ((!wfa_set.i_out_null) << 1) | (!wfa_set.d_out_null);

    wfa_component *wfa = allocate_new_score(alloc_obj, score, lo, hi, kernel, &mramIdx[score], dpu_alloc_mram);

    affine_wfa_compute_offsets(wfa, wfa_set, lo, hi, score, kernel);
    return wfa;
}

void affine_wfa_compute(dpu_alloc_wram_t *dpu_alloc_wram, edit_cigar_t *cigar, char *pattern, char *text, int pattern_length, int text_length, dpu_alloc_mram_t *dpu_alloc_mram)
{

    wfa_component *wfa_score;

    // MRAM base address for every WFA components
    uint32_t *wfa_mramIdx = (uint32_t *)allocate_new(dpu_alloc_wram, (MAX_SCORE + 1) * sizeof(uint32_t));

    wfa_score = allocate_new_score(dpu_alloc_wram, 0, 0, 0, 0, &(wfa_mramIdx[0]), dpu_alloc_mram);

    wfa_score->mwavefront[0] = 0;

    dpu_alloc_wram->WFA_PTR_WRAM = dpu_alloc_wram->CUR_PTR_WRAM;
    uint32_t mem_used_wram_old = dpu_alloc_wram->mem_used_wram;

    int score = 0;
    while (true)
    {

        affine_wfa_extend(wfa_score, pattern, text, pattern_length, text_length, score);

#ifdef REDUCE
        affine_wfa_reduce_wvs(wfa_score, pattern_length, text_length, score);
#endif

        if (affine_wfa_end_reached(wfa_score, pattern_length, text_length, score))
        {
#ifdef BACKTRACE
            store_wfa_cmpnt_to_mram(wfa_score, wfa_mramIdx[score]);
            dpu_alloc_wram->CUR_PTR_WRAM = dpu_alloc_wram->WFA_PTR_WRAM;
            dpu_alloc_wram->mem_used_wram = mem_used_wram_old;
            affine_wavefronts_backtrace(wfa_mramIdx, cigar, pattern, pattern_length, text, text_length, score, dpu_alloc_wram);
#endif
            cigar->score = score;
            return;
        }

        store_wfa_cmpnt_to_mram(wfa_score, wfa_mramIdx[score]);

        // reset wram after every iteration
        dpu_alloc_wram->CUR_PTR_WRAM = dpu_alloc_wram->WFA_PTR_WRAM;
        dpu_alloc_wram->mem_used_wram = mem_used_wram_old;

        ++score;
        if (score > MAX_SCORE)
        {
            cigar->score = score;
            return;
        }
        wfa_score = affine_wfa_compute_next(score, wfa_mramIdx, dpu_alloc_wram, dpu_alloc_mram);
    }
}

int main()
{
    mem_reset();
    uint32_t tasklet_id = me();

    dpu_alloc_wram_t dpu_alloc_wram;
    dpu_alloc_mram_t dpu_alloc_mram;

    // Load parameters
    uint32_t params_m = (uint32_t)DPU_MRAM_HEAP_POINTER;
    DPUParams params_w;
    mram_read((__mram_ptr void const *)params_m, &params_w, ROUND_UP_MULTIPLE_8(sizeof(DPUParams)));
    uint32_t nb_reads_per_dpu = params_w.dpuNumReads;

    if (nb_reads_per_dpu <= 0)
        return 0;

    int nb_reads_per_tasklets = (((nb_reads_per_dpu + NR_TASKLETS) / NR_TASKLETS));
    dpu_alloc_wram = init_dpu_alloc_wram(WRAM_SEGMENT);

    // Divide MRAM segments equally between tasklets
    dpu_alloc_mram.segment_size = ROUND_UP_MULTIPLE_8((64000000 - params_w.mramTotalAllocated) / NR_TASKLETS);
    dpu_alloc_mram.HEAD_PTR_MRAM = ROUND_UP_MULTIPLE_8(dpu_alloc_mram.segment_size * tasklet_id + params_w.mramTotalAllocated);
    dpu_alloc_mram.CUR_PTR_MRAM = dpu_alloc_mram.HEAD_PTR_MRAM;
    dpu_alloc_mram.mem_used_mram = 0;

    uint32_t dpuRequests_m = ((uint32_t)DPU_MRAM_HEAP_POINTER) + params_w.dpuRequests_m;
    uint32_t dpuResults_m = ((uint32_t)DPU_MRAM_HEAP_POINTER) + params_w.dpuResults_m;
    uint32_t dpuPatterns_m = ((uint32_t)DPU_MRAM_HEAP_POINTER) + params_w.dpuPatterns_m;
    uint32_t dpuTexts_m = ((uint32_t)DPU_MRAM_HEAP_POINTER) + params_w.dpuTexts_m;
#ifdef BACKTRACE
    uint32_t dpuOperations_m = ((uint32_t)DPU_MRAM_HEAP_POINTER) + params_w.dpuOperations_m;
#endif
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

            //  DMA transfers can't be of size grater than 2048
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
            // Initialize the operations memory
            memset(cigar->operations, 'M', 2 * READ_SIZE);
#endif

            affine_wfa_compute(&dpu_alloc_wram, cigar, pattern, text, request_w->pattern_len, request_w->text_len, &dpu_alloc_mram);

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

        // reset WRAM and MRAM segments after every read pair alignment
        reset_dpu_alloc_wram(&dpu_alloc_wram);
        dpu_alloc_mram.CUR_PTR_MRAM = dpu_alloc_mram.HEAD_PTR_MRAM;
        dpu_alloc_mram.mem_used_mram = 0;
    }
    return 0;
}