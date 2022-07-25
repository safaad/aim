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
#ifndef COMMON_H__
#define COMMON_H__

#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>

#define DPU_CAPACITY (64 << 20)

#ifndef MATCH
#define MATCH 0
#endif

#ifndef MISMATCH
#define MISMATCH 3
#endif

#ifndef GAP_I
#define GAP_I 4
#endif

#ifndef GAP_D
#define GAP_D 4
#endif

#ifndef MAX_SCORE
#define MAX_SCORE 300
#endif

#ifndef READ_SIZE
#define READ_SIZE 1120
#endif

#define NW_W16

#ifdef NW_W8
typedef int8_t cell_type_t;
#else
#ifdef NW_W16
typedef int16_t cell_type_t;
#else
typedef int cell_type_t;
#endif
#endif

#define MIN(a, b) (((a) <= (b)) ? (a) : (b))
#define MAX(a, b) (((a) >= (b)) ? (a) : (b))
#define ABS(a) (((a) >= 0) ? (a) : -1 * (a))

#define ROUND_UP_MULTIPLE_8(x) ((((x) + 7) / 8) * 8)

typedef struct
{
    int max_operations;
    char *operations;
    int begin_offset;
    int end_offset;
    int score;
} edit_cigar_t;

typedef struct request_t
{
    int pattern_len;
    int text_len;
    int padding; /* Padding to ensure the alignment of the struct */
    uint32_t idx;
} request_t;

typedef struct result_t
{
    int max_operations;
    int begin_offset;
    int end_offset;
    int score;
    uint64_t cycles;
    uint32_t idx;
} result_t;

typedef struct DPUParams
{
    uint32_t dpuNumReads;        /* Number of reads assigned to the DPU */
    uint32_t dpuRequests_m;      /* Base address of the requests in the MRAM */
    uint32_t dpuResults_m;       /* Base address of the results in the MRAM */
    uint32_t dpuPatterns_m;      /* Base address of the patterns sequences in the MRAM */
    uint32_t dpuTexts_m;         /* Base address of the text sequences in the MRAM */
    uint32_t dpuOperations_m;    /* Base address of the traceback operations in the MRAM */
    uint32_t mramTotalAllocated; /* Size of the MRAM memory allocated by the host */
    uint32_t padding;            /* Padding to ensure the alignment of the struct */
} DPUParams;

#endif
