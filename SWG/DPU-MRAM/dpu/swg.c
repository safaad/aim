/*
 *                  The MIT License
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
#include "dpu_allocator_wram.h"
#include "dpu_allocator_mram.h"

void edit_cigar_allocate(
    edit_cigar_t *edit_cigar,
    int pattern_length,
    int text_length)
{
    edit_cigar->max_operations = pattern_length + text_length;
    edit_cigar->begin_offset = edit_cigar->max_operations - 1;
    edit_cigar->end_offset = edit_cigar->max_operations;
    edit_cigar->score = INT32_MIN;
}

void swg_traceback(int num_cols, int num_rows, edit_cigar_t *cigar, dpu_alloc_mram_t *dpu_alloc_mram, dp_cell_t *cell_cache, dp_cell_t *upper_cell_cache, dp_cell_t *diag_cell_cache, dp_cell_t *left_cell_cache)
{

    uint32_t matrix_offset = (uint32_t)DPU_MRAM_HEAP_POINTER + dpu_alloc_mram->CUR_PTR_MRAM;

    char *const operations = cigar->operations;
    int op_sentinel = cigar->end_offset - 1;
    int h, v;
    // Compute traceback
    h = num_cols - 1;
    v = num_rows - 1;
    swg_layer_type swg_layer = swg_M_layer;

    while (h > 0 && v > 0)
    {
        // Load needed cells from MAM
        mram_read((__mram_ptr void const *)(matrix_offset + (num_cols * h + v) * sizeof(dp_cell_t)), cell_cache, ROUND_UP_MULTIPLE_8(sizeof(dp_cell_t)));
        mram_read((__mram_ptr void const *)(matrix_offset + (num_cols * h + v - 1) * sizeof(dp_cell_t)), upper_cell_cache, ROUND_UP_MULTIPLE_8(sizeof(dp_cell_t)));
        mram_read((__mram_ptr void const *)(matrix_offset + (num_cols * (h - 1) + v - 1) * sizeof(dp_cell_t)), diag_cell_cache, ROUND_UP_MULTIPLE_8(sizeof(dp_cell_t)));
        mram_read((__mram_ptr void const *)(matrix_offset + (num_cols * (h - 1) + v) * sizeof(dp_cell_t)), left_cell_cache, ROUND_UP_MULTIPLE_8(sizeof(dp_cell_t)));

        switch (swg_layer)
        {
        case swg_D_layer:
            // Traceback D-matrix
            operations[op_sentinel--] = 'D';
            if (cell_cache->D == upper_cell_cache->M + GAP_O + GAP_E)
            {
                swg_layer = swg_M_layer;
            }
            --v;
            break;
        case swg_I_layer:
            // Traceback I-matrix
            operations[op_sentinel--] = 'I';
            if (cell_cache->I == left_cell_cache->M + GAP_O + GAP_E)
            {
                swg_layer = swg_M_layer;
            }
            --h;
            break;
        case swg_M_layer:
            // Traceback M-matrix
            if (cell_cache->M == cell_cache->D)
            {
                swg_layer = swg_D_layer;
            }
            else if (cell_cache->M == cell_cache->I)
            {
                swg_layer = swg_I_layer;
            }
            else if (cell_cache->M == diag_cell_cache->M + MATCH)
            {
                operations[op_sentinel--] = 'M';
                --h;
                --v;
            }
            else if (cell_cache->M == diag_cell_cache->M + MISMATCH)
            {
                operations[op_sentinel--] = 'X';
                --h;
                --v;
            }
            else
            {

                printf("SWG backtrace. No backtrace operation found");
                exit(1);
            }
            break;
        }
    }
    while (h > 0)
    {
        operations[op_sentinel--] = 'I';
        --h;
    }
    while (v > 0)
    {
        operations[op_sentinel--] = 'D';
        --v;
    }
    cigar->begin_offset = op_sentinel + 1;
}

void swg_compute(char *pattern, char *text, int pattern_length, int text_length, edit_cigar_t *cigar, dpu_alloc_mram_t *dpu_alloc_mram, dp_cell_t *cell_cache, dp_cell_t *upper_cell_cache, dp_cell_t *diag_cell_cache, dp_cell_t *left_cell_cache)
{
    int h, v;
    int num_rows = pattern_length + 1;
    int num_cols = text_length + 1;

    uint32_t matrix_offset = (uint32_t)DPU_MRAM_HEAP_POINTER + dpu_alloc_mram->CUR_PTR_MRAM;

    // Initialize DP
    cell_cache->D = MAX_SCORE;
    cell_cache->I = MAX_SCORE;
    cell_cache->M = 0;

    mram_write(cell_cache, (__mram_ptr void *)(matrix_offset), ROUND_UP_MULTIPLE_8(sizeof(dp_cell_t)));

    for (v = 1; v <= pattern_length; ++v)
    {
        // Init first column
        mram_read((__mram_ptr void const *)(matrix_offset + v * sizeof(dp_cell_t)), cell_cache, ROUND_UP_MULTIPLE_8(sizeof(dp_cell_t)));
        cell_cache->D = GAP_O + v * GAP_E;
        cell_cache->I = MAX_SCORE;
        cell_cache->M = cell_cache->D;
        mram_write(cell_cache, (__mram_ptr void *)(matrix_offset + v * sizeof(dp_cell_t)), ROUND_UP_MULTIPLE_8(sizeof(dp_cell_t)));
    }
    for (h = 1; h <= text_length; ++h)
    {
        // Init first row
        mram_read((__mram_ptr void const *)(matrix_offset + num_cols * h * sizeof(dp_cell_t)), cell_cache, ROUND_UP_MULTIPLE_8(sizeof(dp_cell_t)));
        cell_cache->D = MAX_SCORE;
        cell_cache->I = GAP_O + h * GAP_E;
        cell_cache->M = cell_cache->I;
        mram_write(cell_cache, (__mram_ptr void *)(matrix_offset + num_cols * h * sizeof(dp_cell_t)), ROUND_UP_MULTIPLE_8(sizeof(dp_cell_t)));
    }

    // Compute DP
    int score = 0;
    for (h = 1; h <= text_length; ++h)
    {
        for (v = 1; v <= pattern_length; ++v)
        {
            mram_read((__mram_ptr void const *)(matrix_offset + (num_cols * h + v - 1) * sizeof(dp_cell_t)), upper_cell_cache, ROUND_UP_MULTIPLE_8(sizeof(dp_cell_t)));
            mram_read((__mram_ptr void const *)(matrix_offset + (num_cols * (h - 1) + v - 1) * sizeof(dp_cell_t)), diag_cell_cache, ROUND_UP_MULTIPLE_8(sizeof(dp_cell_t)));
            mram_read((__mram_ptr void const *)(matrix_offset + (num_cols * (h - 1) + v) * sizeof(dp_cell_t)), left_cell_cache, ROUND_UP_MULTIPLE_8(sizeof(dp_cell_t)));

            // Update DP.D
            cell_size_t del_new = upper_cell_cache->M + GAP_O + GAP_E;
            cell_size_t del_ext = upper_cell_cache->D + GAP_E;
            cell_size_t del = MIN(del_new, del_ext);
            cell_cache->D = del;
            // Update DP.I
            cell_size_t ins_new = left_cell_cache->M + GAP_O + GAP_E;
            cell_size_t ins_ext = left_cell_cache->I + GAP_E;
            cell_size_t ins = MIN(ins_new, ins_ext);
            cell_cache->I = ins;
            // Update DP.M
            cell_size_t m_match = diag_cell_cache->M + ((pattern[v - 1] == text[h - 1]) ? MATCH : MISMATCH);
            cell_cache->M = MIN(m_match, MIN(ins, del));
            score = cell_cache->M;
            mram_write(cell_cache, (__mram_ptr void *)(matrix_offset + (num_cols * h + v) * sizeof(dp_cell_t)), ROUND_UP_MULTIPLE_8(sizeof(dp_cell_t)));
        }
    }

    cigar->score = score;
#ifdef BACKTRACE
    // Compute traceback
    swg_traceback(num_cols, num_rows, cigar, dpu_alloc_mram, cell_cache, upper_cell_cache, diag_cell_cache, left_cell_cache);
#endif
}

int main()
{
    mem_reset();
    uint32_t tasklet_id = me();

    dpu_alloc_mram_t dpu_alloc_mram;

    // Load parameters
    uint32_t params_m = (uint32_t)DPU_MRAM_HEAP_POINTER;
    DPUParams params_w;
    mram_read((__mram_ptr void const *)params_m, &params_w, ROUND_UP_MULTIPLE_8(sizeof(DPUParams)));
    uint32_t nb_reads_per_dpu = params_w.dpuNumReads;

    if (nb_reads_per_dpu <= 0)
        return 0;

    int nb_reads_per_tasklets = ROUND_UP_MULTIPLE_8(((nb_reads_per_dpu + NR_TASKLETS) / NR_TASKLETS));

    uint32_t dpuRequests_m = ((uint32_t)DPU_MRAM_HEAP_POINTER) + params_w.dpuRequests_m;
    uint32_t dpuResults_m = ((uint32_t)DPU_MRAM_HEAP_POINTER) + params_w.dpuResults_m;
    uint32_t dpuPatterns_m = ((uint32_t)DPU_MRAM_HEAP_POINTER) + params_w.dpuPatterns_m;
    uint32_t dpuTexts_m = ((uint32_t)DPU_MRAM_HEAP_POINTER) + params_w.dpuTexts_m;
#ifdef BACKTRACE
    uint32_t dpuOperations_m = ((uint32_t)DPU_MRAM_HEAP_POINTER) + params_w.dpuOperations_m;
#endif
    // Check if the size of the DP-table fits in the MRAM for all tasklets
    if (ROUND_UP_MULTIPLE_8(READ_SIZE * READ_SIZE * NR_TASKLETS * sizeof(dp_cell_t)) + params_w.mramTotalAllocated >= 64000000)
    {
        printf("Insufficient MRAM memory\n");
        exit(-1);
    }

    // Get the base address of the DP-table in the MRAM
    dpu_alloc_mram.HEAD_PTR_MRAM = ROUND_UP_MULTIPLE_8(READ_SIZE * READ_SIZE * tasklet_id * sizeof(dp_cell_t)) + params_w.mramTotalAllocated;
    dpu_alloc_mram.CUR_PTR_MRAM = dpu_alloc_mram.HEAD_PTR_MRAM;
    dpu_alloc_mram.mem_used_mram = 0;

    request_t *request_w = (request_t *)mem_alloc(ROUND_UP_MULTIPLE_8(sizeof(request_t)));
    result_t *result_w = (result_t *)mem_alloc(ROUND_UP_MULTIPLE_8(sizeof(result_t)));

    edit_cigar_t *cigar;
    cigar = (edit_cigar_t *)mem_alloc(ROUND_UP_MULTIPLE_8(sizeof(edit_cigar_t)));

    char *pattern = (char *)mem_alloc(ROUND_UP_MULTIPLE_8(READ_SIZE));
    char *text = (char *)mem_alloc(ROUND_UP_MULTIPLE_8(READ_SIZE));

    // We only need 4 cache cells in the WRAM
    dp_cell_t *cell_cache = (dp_cell_t *)mem_alloc(ROUND_UP_MULTIPLE_8(sizeof(dp_cell_t)));
    dp_cell_t *diag_cell_cache = (dp_cell_t *)mem_alloc(ROUND_UP_MULTIPLE_8(sizeof(dp_cell_t)));
    dp_cell_t *upper_cell_cache = (dp_cell_t *)mem_alloc(ROUND_UP_MULTIPLE_8(sizeof(dp_cell_t)));
    dp_cell_t *left_cell_cache = (dp_cell_t *)mem_alloc(ROUND_UP_MULTIPLE_8(sizeof(dp_cell_t)));
#ifdef BACKTRACE
    cigar->operations = (char *)mem_alloc(2 * READ_SIZE);
#endif

    for (int read_nb = 0; read_nb < nb_reads_per_tasklets; ++read_nb)
    {

        if (read_nb + tasklet_id * nb_reads_per_tasklets < nb_reads_per_dpu)
        {

            mram_read((__mram_ptr void const *)(dpuRequests_m + (read_nb + tasklet_id * nb_reads_per_tasklets) * (sizeof(request_t))), request_w, ROUND_UP_MULTIPLE_8(sizeof(request_t)));

            // DMA transfers can't be greater than 2048
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
            edit_cigar_allocate(cigar, request_w->pattern_len, request_w->text_len);

#ifdef BACKTRACE
            // Initialize the operations memory
            memset(cigar->operations, 'M', 2 * READ_SIZE);
#endif

            result_w->idx = request_w->idx;

            swg_compute(pattern, text, request_w->pattern_len, request_w->text_len, cigar, &dpu_alloc_mram, cell_cache, upper_cell_cache, diag_cell_cache, left_cell_cache);

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
    }

    return 0;
}