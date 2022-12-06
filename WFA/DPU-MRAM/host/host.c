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
#define _GNU_SOURCE
#include <stdio.h>
#include "timer.h"
#include "common.h"
#include "mram-management.h"
#include <time.h>
#include <dpu.h>
#ifndef DPU_BINARY
#define DPU_BINARY "build/wfa_dpu"
#endif

#ifndef ENERGY
#define ENERGY 0
#endif
#if ENERGY
#include <dpu_probe.h>
#endif

void edit_cigar_print(
    edit_cigar_t *const edit_cigar, FILE *out)
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
            fprintf(out, "%d%c", last_op_length, last_op);
            last_op = edit_cigar->operations[i];
            last_op_length = 1;
        }
    }
    fprintf(out, "%d%c\n", last_op_length, last_op);
}

uint32_t get_reads(FILE *in, request_t *dpu_requests, char *dpu_patterns, char *dpu_texts, uint32_t nb_reads_per_dpu, int nb_sent_requests, int total_nb_reads)
{
    char *line1 = NULL, *line2 = NULL;
    int line1_length = 0, line2_length = 0;
    size_t line1_allocated = 0, line2_allocated = 0;
    char *pattern;
    char *text;
    int pattern_length;
    int text_length;

    uint32_t nb_reads = 0;
    for (nb_reads = 0; nb_reads < nb_reads_per_dpu; ++nb_reads)
    {
        line1_length = getline(&line1, &line1_allocated, in);
        if (line1_length == -1)
            break;

        line2_length = getline(&line2, &line2_allocated, in);
        if (line2_length == -1)
            break;

        pattern = line1 + 1;
        pattern_length = line1_length - 2;
        pattern[pattern_length] = '\0';
        text = line2 + 1;
        text_length = line2_length - 2;
        text[text_length] = '\0';

        if (text_length > READ_SIZE || pattern_length > READ_SIZE)
        {
            printf("READ LENGTH less than length of the input reads");
            exit(0);
        }
        if (pattern == NULL || text == NULL)
            exit(0);
        strcpy(&dpu_patterns[nb_reads * (READ_SIZE)], pattern);
        strcpy(&dpu_texts[nb_reads * (READ_SIZE)], text);
        dpu_requests[nb_reads].pattern_len = pattern_length;
        dpu_requests[nb_reads].text_len = text_length;

        dpu_requests[nb_reads].idx = nb_reads + nb_sent_requests;
    }
    return nb_reads;
}

int main(int argc, char *argv[])
{

    // Timing and profiling
    Timer timer;
    float loadTime = 0.0f, dpuTime = 0.0f, retrieveTime = 0.0f;
#if ENERGY
    struct dpu_probe_t probe;
    DPU_ASSERT(dpu_probe_init("energy_probe", &probe));
#endif

    struct dpu_set_t dpu_set, dpu;
    uint32_t nr_of_dpus;

    if (argc != 4)
    {
        printf("wrong number of arguments\n");
        exit(1);
    }

    char *in = argv[1];                      // input read pairs file
    char *out = argv[2];                     // output file
    uint32_t total_nb_reads = atoi(argv[3]); // total number of reads to align

    FILE *input_file = NULL;
    FILE *output_file = NULL;
    input_file = fopen(in, "r");
    output_file = fopen(out, "w");
    FILE *dpu_file = fopen("dpu-out", "w");
    if (input_file == NULL)
    {
        fprintf(stderr, "Input file '%s' couldn't be opened\n", in);
        exit(1);
    }
    if (output_file == NULL)
    {
        fprintf(stderr, "Output file '%s' couldn't be opened\n", out);
        exit(1);
    }
    if (total_nb_reads <= 0)
    {
        fprintf(stderr, "Invalid nb of reads\n");
        exit(1);
    }
    if (total_nb_reads <= NR_DPUS)
    {
        printf("Allocated DPUs more than needed\n");
        exit(1);
    }

    DPU_ASSERT(dpu_alloc(NR_DPUS, NULL, &dpu_set));
    DPU_ASSERT(dpu_load(dpu_set, DPU_BINARY, NULL));
    DPU_ASSERT(dpu_get_nr_dpus(dpu_set, &nr_of_dpus));
    printf("Allocated %d DPU(s)\n", nr_of_dpus);

    uint32_t nb_reads_per_dpu = ROUND_UP_MULTIPLE_8(((total_nb_reads) / nr_of_dpus));
    printf("NumReads per dpu = %u\n", nb_reads_per_dpu);
    int nb_sent_requests = 0;

    // Allocate Buffer
    struct DPUParams dpuParams[nr_of_dpus];
    request_t *dpu_requests[nr_of_dpus];
    char *dpu_patterns[nr_of_dpus];
    char *dpu_texts[nr_of_dpus];

    for (int dpu_idx = 0; dpu_idx < nr_of_dpus; ++dpu_idx)
    {
        dpu_requests[dpu_idx] = (request_t *)malloc(nb_reads_per_dpu * (sizeof(request_t)));
        dpu_patterns[dpu_idx] = (char *)malloc(nb_reads_per_dpu * (READ_SIZE));
        dpu_texts[dpu_idx] = (char *)malloc(nb_reads_per_dpu * (READ_SIZE));
        uint32_t nb_reads = get_reads(input_file, dpu_requests[dpu_idx], dpu_patterns[dpu_idx], dpu_texts[dpu_idx], nb_reads_per_dpu, nb_sent_requests, total_nb_reads);
        nb_sent_requests += nb_reads;
        dpuParams[dpu_idx].dpuNumReads = nb_reads;
    }

    startTimer(&timer);
    uint32_t each_dpu;
    uint32_t dpuParams_m = 0;

    DPU_FOREACH(dpu_set, dpu, each_dpu)
    {
        // Allocate needed MRAM memory to store Read Pairs Requests and Results
        struct mram_heap_allocator_t allocator;
        init_allocator(&allocator);
        dpuParams_m = mram_heap_alloc(&allocator, (sizeof(struct DPUParams)));
        uint32_t dpuRequests_m = mram_heap_alloc(&allocator, nb_reads_per_dpu * (sizeof(request_t)));
        uint32_t dpuResults_m = mram_heap_alloc(&allocator, (nb_reads_per_dpu * (sizeof(result_t))));
        uint32_t dpuPatterns_m = mram_heap_alloc(&allocator, nb_reads_per_dpu * (READ_SIZE));
        uint32_t dpuTexts_m = mram_heap_alloc(&allocator, nb_reads_per_dpu * (READ_SIZE));
#ifdef BACKTRACE
        uint32_t dpuOperations_m = mram_heap_alloc(&allocator, nb_reads_per_dpu * (2 * READ_SIZE));
#else
        uint32_t dpuOperations_m = 0;
#endif
        assert((sizeof(request_t)) % 8 == 0 && "Requests must be a multiple of 8 bytes!");
        assert((sizeof(result_t)) % 8 == 0 && "Results must be a multiple of 8 bytes!");
        assert((nb_reads_per_dpu * (READ_SIZE)) % 8 == 0 && "Input sequences must be a multiple of 8 bytes!");
        assert((sizeof(struct DPUParams)) % 8 == 0 && "DPUParams must be a multiple of 8 bytes!");

        dpuParams[each_dpu].dpuRequests_m = dpuRequests_m;
        dpuParams[each_dpu].dpuResults_m = dpuResults_m;
        dpuParams[each_dpu].dpuPatterns_m = dpuPatterns_m;
        dpuParams[each_dpu].dpuTexts_m = dpuTexts_m;
        dpuParams[each_dpu].dpuOperations_m = dpuOperations_m;
        dpuParams[each_dpu].mramTotalAllocated = ROUND_UP_MULTIPLE_8(allocator.totalAllocated);
    }

    // Parallel Transfers of the Input
    printf("Copying data to DPU\n");
    // Transfer DPU Params
    DPU_FOREACH(dpu_set, dpu, each_dpu)
    {
        DPU_ASSERT(dpu_prepare_xfer(dpu, (uint8_t *)&dpuParams[each_dpu]));
    }
    DPU_ASSERT(dpu_push_xfer(dpu_set, DPU_XFER_TO_DPU, DPU_MRAM_HEAP_POINTER_NAME, dpuParams_m, ROUND_UP_MULTIPLE_8(sizeof(struct DPUParams)), DPU_XFER_DEFAULT));
    // Transfer the Requests
    DPU_FOREACH(dpu_set, dpu, each_dpu)
    {
        DPU_ASSERT(dpu_prepare_xfer(dpu, (uint8_t *)dpu_requests[each_dpu]));
    }
    DPU_ASSERT(dpu_push_xfer(dpu_set, DPU_XFER_TO_DPU, DPU_MRAM_HEAP_POINTER_NAME, dpuParams[0].dpuRequests_m, nb_reads_per_dpu * ((sizeof(request_t))), DPU_XFER_DEFAULT));
    // Transfer Pattern sequences
    DPU_FOREACH(dpu_set, dpu, each_dpu)
    {
        DPU_ASSERT(dpu_prepare_xfer(dpu, (uint8_t *)dpu_patterns[each_dpu]));
    }
    DPU_ASSERT(dpu_push_xfer(dpu_set, DPU_XFER_TO_DPU, DPU_MRAM_HEAP_POINTER_NAME, dpuParams[0].dpuPatterns_m, nb_reads_per_dpu * (READ_SIZE), DPU_XFER_DEFAULT));
    // Transfer Text Sequences
    DPU_FOREACH(dpu_set, dpu, each_dpu)
    {
        DPU_ASSERT(dpu_prepare_xfer(dpu, (uint8_t *)dpu_texts[each_dpu]));
    }
    DPU_ASSERT(dpu_push_xfer(dpu_set, DPU_XFER_TO_DPU, DPU_MRAM_HEAP_POINTER_NAME, dpuParams[0].dpuTexts_m, nb_reads_per_dpu * (READ_SIZE), DPU_XFER_DEFAULT));

    stopTimer(&timer);
    loadTime += getElapsedTime(timer);
    printf("CPU-DPU: %f ms\n", loadTime * 1e3);

    for (int dpu_idx = 0; dpu_idx < nr_of_dpus; ++dpu_idx)
    {
        free(dpu_requests[dpu_idx]);
        free(dpu_patterns[dpu_idx]);
        free(dpu_texts[dpu_idx]);
    }
    fclose(input_file);

    // Run the DPU Kernel
    printf("Run program on DPU(s)\n");
    startTimer(&timer);
#if ENERGY
    DPU_ASSERT(dpu_probe_start(&probe));
#endif

    DPU_ASSERT(dpu_launch(dpu_set, DPU_SYNCHRONOUS));

#if ENERGY
    DPU_ASSERT(dpu_probe_stop(&probe));
    double energy;
    DPU_ASSERT(dpu_probe_get(&probe, DPU_ENERGY, DPU_AVERAGE, &energy));
    PRINT_INFO(p.verbosity >= 1, "    DPU Energy: %f J", energy);
#endif
    stopTimer(&timer);
    dpuTime += getElapsedTime(timer);
    printf("DPU Kernel: %f ms\n", dpuTime * 1e3);

    result_t *dpuResults[nr_of_dpus];
#ifdef BACKTRACE
    char *dpuOperations[nr_of_dpus];
#endif
    for (int dpu = 0; dpu < nr_of_dpus; ++dpu)
    {
        dpuResults[dpu] = (result_t *)malloc(nb_reads_per_dpu * (sizeof(result_t)));
#ifdef BACKTRACE
        dpuOperations[dpu] = (char *)malloc(nb_reads_per_dpu * (2 * READ_SIZE));
#endif
    }

    // DPU-CPU Transfers
    printf("Retrieve results\n");
    startTimer(&timer);
    DPU_FOREACH(dpu_set, dpu, each_dpu)
    {
        DPU_ASSERT(dpu_prepare_xfer(dpu, dpuResults[each_dpu]));
    }
    DPU_ASSERT(dpu_push_xfer(dpu_set, DPU_XFER_FROM_DPU, DPU_MRAM_HEAP_POINTER_NAME, dpuParams[0].dpuResults_m, nb_reads_per_dpu * ((sizeof(result_t))), DPU_XFER_DEFAULT));
#ifdef BACKTRACE
    DPU_FOREACH(dpu_set, dpu, each_dpu)
    {
        DPU_ASSERT(dpu_prepare_xfer(dpu, dpuOperations[each_dpu]));
    }
    DPU_ASSERT(dpu_push_xfer(dpu_set, DPU_XFER_FROM_DPU, DPU_MRAM_HEAP_POINTER_NAME, dpuParams[0].dpuOperations_m, nb_reads_per_dpu * (2 * READ_SIZE), DPU_XFER_DEFAULT));
#endif
    stopTimer(&timer);
    retrieveTime += getElapsedTime(timer);
    printf("DPU-CPU: %f ms\n", retrieveTime * 1e3);

    for (int dpu = 0; dpu < nr_of_dpus; ++dpu)
    {
        int i;
#ifdef BACKTRACE
        char *operations = dpuOperations[dpu];
#endif
        for (i = 0; i < dpuParams[dpu].dpuNumReads; ++i)
        {
            fprintf(output_file, "%d, %d, \n", dpuResults[dpu][i].idx, dpuResults[dpu][i].score);

#ifdef BACKTRACE
            edit_cigar_t cigar;
            cigar.score = dpuResults[dpu][i].score;
            cigar.max_operations = dpuResults[dpu][i].max_operations;
            cigar.begin_offset = dpuResults[dpu][i].begin_offset;
            cigar.end_offset = dpuResults[dpu][i].end_offset;
            cigar.operations = (char *)malloc(ROUND_UP_MULTIPLE_8(cigar.max_operations));
            strncpy(cigar.operations, &(operations[i * 2 * READ_SIZE]), ROUND_UP_MULTIPLE_8(cigar.max_operations));
            edit_cigar_print(&cigar, output_file);
#endif
        }
    }

    // DPU Logs
    uint32_t dpuIdx;
    DPU_FOREACH(dpu_set, dpu, dpuIdx)
    {
        fprintf(dpu_file, "DPU %u:", dpuIdx);
        DPU_ASSERT(dpu_log_read(dpu, dpu_file));
        ++dpuIdx;
    }

    // Free
    for (int dpu = 0; dpu < nr_of_dpus; ++dpu)
    {
        free(dpuResults[dpu]);
#ifdef BACKTRACE
        free(dpuOperations[dpu]);
#endif
    }
    DPU_ASSERT(dpu_free(dpu_set));

    fclose(dpu_file);
    fclose(output_file);
    return 0;
}
