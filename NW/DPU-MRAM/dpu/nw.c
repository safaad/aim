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
 * DESCRIPTION: Dynam
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

#define CACHE_SIZE (ROUND_UP_MULTIPLE_8(sizeof(cell_type_t)))

void edit_cigar_print(
    edit_cigar_t *const edit_cigar)
{
    char last_op = edit_cigar->operations[edit_cigar->begin_offset];
    int last_op_length = 1;
    int i;
    for (i = edit_cigar->begin_offset + 1; i < edit_cigar->end_offset; ++i)
    {
        if (edit_cigar->operations[i] == last_op)
        {
            ++last_op_length;
        }
        else
        {
            printf("%d%c", last_op_length, last_op);
            last_op = edit_cigar->operations[i];
            last_op_length = 1;
        }
    }
    printf("%d%c\n", last_op_length, last_op);
}

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

void nw_traceback(int num_cols, int num_rows, edit_cigar_t *cigar, dpu_alloc_mram_t *dpu_alloc_mram,int tasklet_id, cell_type_t *cell_cache, cell_type_t *upper_cell_cache, cell_type_t *diag_cell_cache, cell_type_t *left_cell_cache)
{
    uint32_t matrix_offset = (uint32_t)DPU_MRAM_HEAP_POINTER + dpu_alloc_mram->CUR_PTR_MRAM;
    char *const operations = cigar->operations;
    int op_sentinel = cigar->end_offset - 1;
    int h, v;
    // Compute traceback
    h = num_cols - 1;
    v = num_rows - 1;

    while (h > 0 && v > 0)
    {
        int cell_offset = (num_cols * h + v) & (-4);
        int cell_index = (num_cols * h + v) & 3;

        int upper_cell_offset = (num_cols * h + v - 1) & (-4);
        int upper_cell_index = (num_cols * h + v - 1) & 3;

        int left_cell_offset = ((num_cols * (h - 1) + v)) & (-4);
        int left_cell_index = (num_cols * (h - 1) + v) & 3;

        int diag_cell_offset = ((num_cols * (h - 1) + v - 1)) & (-4);
        int diag_cell_index = (num_cols * (h - 1) + v - 1) & 3;

        mram_read((__mram_ptr void const *)(matrix_offset + upper_cell_offset*sizeof(cell_type_t)), upper_cell_cache, CACHE_SIZE);
        mram_read((__mram_ptr void const *)(matrix_offset + diag_cell_offset*sizeof(cell_type_t)), diag_cell_cache, CACHE_SIZE);
        mram_read((__mram_ptr void const *)(matrix_offset + left_cell_offset*sizeof(cell_type_t)), left_cell_cache, CACHE_SIZE);
        mram_read((__mram_ptr const void *)(matrix_offset + cell_offset*sizeof(cell_type_t)), cell_cache, CACHE_SIZE);

        if (cell_cache[cell_index] == upper_cell_cache[upper_cell_index] + GAP_D)
        {
            operations[op_sentinel--] = 'D';
            --v;
        }
        else if (cell_cache[cell_index] == left_cell_cache[left_cell_index] + GAP_I)
        {
            operations[op_sentinel--] = 'I';
            --h;
        }
        else
        {
            operations[op_sentinel--] =
                (cell_cache[cell_index] == diag_cell_cache[diag_cell_index] + MISMATCH) ? 'X' : 'M';
            --h;
            --v;
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

void nw_compute(char *pattern, char *text, int pattern_length, int text_length, edit_cigar_t *cigar,dpu_alloc_mram_t *dpu_alloc_mram, uint32_t tasklet_id, cell_type_t *cell_cache, cell_type_t *upper_cell_cache, cell_type_t *diag_cell_cache, cell_type_t *left_cell_cache)
{
    int h, v;
    int num_rows = pattern_length + 1;
    int num_cols = text_length + 1;

    // DP_table offset relative to each tasklet
    uint32_t matrix_offset = (uint32_t)DPU_MRAM_HEAP_POINTER + dpu_alloc_mram->CUR_PTR_MRAM;

    // Cell base address in the MRAM must be aligned to 8
    int cell_offset = (matrix_offset) & (-8);
    int cell_index = ((matrix_offset)-cell_offset);

    mram_read((__mram_ptr void const *)(matrix_offset), cell_cache, CACHE_SIZE);
    cell_type_t cell = 0;
    cell_cache[cell_index] = cell;
    mram_write(cell_cache,((__mram_ptr void *)(matrix_offset)), CACHE_SIZE);

    for (v = 1; v <= pattern_length; ++v)
    {
        // Init first column
        // Cell base address in the MRAM must be aligned to 8
        int cell_offset = (v) & (-4);
        int cell_index = v & 3;
        mram_read((__mram_ptr void const *)(matrix_offset + cell_offset*sizeof(cell_type_t)), cell_cache, CACHE_SIZE);
        cell = cell + GAP_D;
        cell_cache[cell_index] = cell;
        mram_write(cell_cache, (__mram_ptr void *)(matrix_offset + cell_offset*sizeof(cell_type_t)), CACHE_SIZE);
    }

    cell = 0;
    for (h = 1; h <= text_length; ++h)
    {
        // Init first row
        // Cell base address in the MRAM must be aligned to 8
        int cell_offset = (num_cols * h) & (-4);
        int cell_index = (num_cols * h) & 3;

        mram_read((__mram_ptr void const *)(matrix_offset + cell_offset*sizeof(cell_type_t)), cell_cache, CACHE_SIZE);
        cell = cell + GAP_I;
        cell_cache[cell_index] = cell;
        mram_write(cell_cache, (__mram_ptr void *)(matrix_offset + cell_offset*sizeof(cell_type_t)), CACHE_SIZE);
    }

    // Compute DP
    int score = 0;
    for (h = 1; h <= text_length; ++h)
    {
        for (v = 1; v <= pattern_length; ++v)
        {
            // Cell base address in the MRAM must be aligned to 8
            int cell_offset = (num_cols * h + v) & (-4);
            int cell_index = (num_cols * h + v) & 3;

            int upper_cell_offset = (num_cols * h + v - 1) & (-4);
            int upper_cell_index = (num_cols * h + v - 1) & 3;

            int left_cell_offset = ((num_cols * (h - 1) + v)) & (-4);
            int left_cell_index = (num_cols * (h - 1) + v) & 3;

            int diag_cell_offset = ((num_cols * (h - 1) + v - 1)) & (-4);
            int diag_cell_index = (num_cols * (h - 1) + v - 1) & 3;

            mram_read((__mram_ptr void const *)(matrix_offset + upper_cell_offset*sizeof(cell_type_t)), upper_cell_cache, CACHE_SIZE);
            mram_read((__mram_ptr void const *)(matrix_offset + diag_cell_offset*sizeof(cell_type_t)), diag_cell_cache, CACHE_SIZE);
            mram_read((__mram_ptr void const *)(matrix_offset + left_cell_offset*sizeof(cell_type_t)), left_cell_cache, CACHE_SIZE);
            mram_read((__mram_ptr const void *)(matrix_offset + cell_offset*sizeof(cell_type_t)), cell_cache, CACHE_SIZE);
            // Del
            cell_type_t del = (cell_type_t)upper_cell_cache[upper_cell_index] + GAP_D;
            // Ins
            cell_type_t ins = (cell_type_t)left_cell_cache[left_cell_index] + GAP_I;
            // Match
            cell_type_t m_match = (cell_type_t)diag_cell_cache[diag_cell_index] + ((pattern[v - 1] == text[h - 1]) ? 0 : MISMATCH);

            score = (cell_type_t)MIN(m_match, MIN(ins, del));

            cell_cache[cell_index] = score;

            mram_write(cell_cache, (__mram_ptr void *)(matrix_offset + cell_offset*sizeof(cell_type_t)), CACHE_SIZE);
        }
    }
    cigar->score = score;
#ifdef BACKTRACE
    // Compute traceback
    nw_traceback(num_cols, num_rows, cigar, dpu_alloc_mram, tasklet_id, cell_cache, upper_cell_cache, diag_cell_cache, left_cell_cache);
#endif
}

int main()
{
    mem_reset();
    uint32_t tasklet_id = me();

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

    dpu_alloc_mram_t dpu_alloc_mram;
    // Get the base address of the DP-table in the MRAM
    dpu_alloc_mram.HEAD_PTR_MRAM = ROUND_UP_MULTIPLE_8(READ_SIZE * READ_SIZE * tasklet_id * sizeof(cell_type_t)) + params_w.mramTotalAllocated;
    dpu_alloc_mram.CUR_PTR_MRAM = dpu_alloc_mram.HEAD_PTR_MRAM;
    dpu_alloc_mram.mem_used_mram = 0;

    request_t *request_w = (request_t *)mem_alloc(ROUND_UP_MULTIPLE_8(sizeof(request_t)));
    result_t *result_w = (result_t *)mem_alloc(ROUND_UP_MULTIPLE_8(sizeof(result_t)));

    edit_cigar_t *cigar;
    cigar = (edit_cigar_t *)mem_alloc(ROUND_UP_MULTIPLE_8(sizeof(edit_cigar_t)));

    char *pattern = (char *)mem_alloc(ROUND_UP_MULTIPLE_8(READ_SIZE));
    char *text = (char *)mem_alloc(ROUND_UP_MULTIPLE_8(READ_SIZE));

    // Only 4 cell caches are needed in the WRAM
    cell_type_t *cell_cache = (cell_type_t *)mem_alloc(CACHE_SIZE);
    cell_type_t *diag_cell_cache = (cell_type_t *)mem_alloc(CACHE_SIZE);
    cell_type_t *upper_cell_cache = (cell_type_t *)mem_alloc(CACHE_SIZE);
    cell_type_t *left_cell_cache = (cell_type_t *)mem_alloc(CACHE_SIZE);

#ifdef BACKTRACE
    cigar->operations = (char *)mem_alloc(2 * READ_SIZE);
    memset(cigar->operations, 'M', 2 * READ_SIZE);
#endif

    for (int read_nb = 0; read_nb < nb_reads_per_tasklets; ++read_nb)
    {
        if (read_nb + tasklet_id * nb_reads_per_tasklets < nb_reads_per_dpu)
        {

            mram_read((__mram_ptr void const *)(dpuRequests_m + (read_nb + tasklet_id * nb_reads_per_tasklets) * (sizeof(request_t))), request_w, ROUND_UP_MULTIPLE_8(sizeof(request_t)));

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

            nw_compute(pattern, text, request_w->pattern_len, request_w->text_len, cigar, &dpu_alloc_mram, tasklet_id, cell_cache, upper_cell_cache, diag_cell_cache, left_cell_cache);

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
    }
    return 0;
}