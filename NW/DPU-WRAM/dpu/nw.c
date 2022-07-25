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
    edit_cigar->score = 0;
}

void nw_traceback(int num_cols, int num_rows, edit_cigar_t *cigar, cell_type_t *dp_table)
{
    char *const operations = cigar->operations;
    int op_sentinel = cigar->end_offset - 1;
    int h, v;
    // Compute traceback
    h = num_cols - 1;
    v = num_rows - 1;

    while (h > 0 && v > 0)
    {
        if (dp_table[num_cols * h + v] == dp_table[num_cols * h + v - 1] + GAP_D)
        {
            operations[op_sentinel--] = 'D';
            --v;
        }
        else if (dp_table[num_cols * h + v] == dp_table[num_cols * (h - 1) + v] + GAP_I)
        {
            operations[op_sentinel--] = 'I';
            --h;
        }
        else
        {
            operations[op_sentinel--] =
                (dp_table[num_cols * h + v] == dp_table[num_cols * (h - 1) + v - 1] + MISMATCH) ? 'X' : 'M';
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

void nw_compute(char *pattern, char *text, int pattern_length, int text_length, edit_cigar_t *cigar, cell_type_t *dp_table)
{
    int h, v;
    int num_rows = pattern_length + 1;
    int num_cols = text_length + 1;

    int cell = 0;
    dp_table[0] = cell;

    for (v = 1; v <= pattern_length; ++v)
    {
        // Initialize first column
        cell = cell + GAP_D;
        dp_table[v] = cell;
    }
    cell = 0;
    for (h = 1; h <= text_length; ++h)
    {
        // Initialize first row
        cell = cell + GAP_I;
        dp_table[num_cols * h] = cell;
    }

    // Compute DP
    cell_type_t score = 0;
    for (h = 1; h <= text_length; ++h)
    {
        for (v = 1; v <= pattern_length; ++v)
        {
            // Del
            cell_type_t del = (cell_type_t)dp_table[(num_cols * (h) + v - 1)] + GAP_D;
            // Ins
            cell_type_t ins = (cell_type_t)dp_table[num_cols * (h - 1) + v] + GAP_I;
            // Match
            cell_type_t m_match = (cell_type_t)dp_table[(num_cols * (h - 1) + v - 1)] + ((pattern[v - 1] == text[h - 1]) ? 0 : MISMATCH);

            score = dp_table[num_cols * h + v] = (cell_type_t)MIN(m_match, MIN(ins, del));
        }
    }
    cigar->score = (int)score;
#ifdef BACKTRACE
    // Compute traceback
    nw_traceback(num_cols, num_rows, cigar, dp_table);
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
    uint32_t dpuOperations_m = ((uint32_t)DPU_MRAM_HEAP_POINTER) + params_w.dpuOperations_m;

    request_t *request_w = (request_t *)mem_alloc(ROUND_UP_MULTIPLE_8(sizeof(request_t)));
    result_t *result_w = (result_t *)mem_alloc(ROUND_UP_MULTIPLE_8(sizeof(result_t)));

    edit_cigar_t *cigar;
    cigar = (edit_cigar_t *)mem_alloc(ROUND_UP_MULTIPLE_8(sizeof(edit_cigar_t)));

    char *pattern = (char *)mem_alloc(ROUND_UP_MULTIPLE_8(READ_SIZE));
    char *text = (char *)mem_alloc(ROUND_UP_MULTIPLE_8(READ_SIZE));

    // Each taasklet has a DP-table stored in WRAM reused
    cell_type_t *dp_table = (cell_type_t *)mem_alloc(ROUND_UP_MULTIPLE_8(READ_SIZE * READ_SIZE * sizeof(cell_type_t)));

#ifdef BACKTRACE
    cigar->operations = (char *)mem_alloc(2 * READ_SIZE);
    // initialize traceback operations segment
    memset(cigar->operations, 'M', 2 * READ_SIZE);
#endif

    for (int read_nb = 0; read_nb < nb_reads_per_tasklets; ++read_nb)
    {
        if (read_nb + tasklet_id * nb_reads_per_tasklets < nb_reads_per_dpu)
        {

            mram_read((__mram_ptr void const *)(dpuRequests_m + (read_nb + tasklet_id * nb_reads_per_tasklets) * (sizeof(request_t))), request_w, ROUND_UP_MULTIPLE_8(sizeof(request_t)));

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
            edit_cigar_allocate(cigar, request_w->pattern_len, request_w->text_len);

            nw_compute(pattern, text, request_w->pattern_len, request_w->text_len, cigar, dp_table);

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