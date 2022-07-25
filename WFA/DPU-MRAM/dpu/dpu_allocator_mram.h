#ifndef MRAM_ALLOCATOR_
#define MRAM_ALLOCATOR_

#include "common.h"
#include "dpu_allocator_wram.h"

typedef struct dpu_alloc_mram_t
{
    uint32_t segment_size;
    uint32_t HEAD_PTR_MRAM;
    uint32_t CUR_PTR_MRAM;
    uint32_t mem_used_mram;
} dpu_alloc_mram_t;

void add_wfa_cmpnt_to_mram(uint32_t *mramIdx, uint32_t size, dpu_alloc_mram_t *dpu_alloc_mram);

wfa_component *load_wfa_cmpnt_from_mram(dpu_alloc_wram_t *allocator, uint32_t mramIdx);
wfa_component *load_mwavefront_cmpnt_from_mram(dpu_alloc_wram_t *allocator, uint32_t mramIdx);
wfa_component *load_idwavefront_cmpnt_from_mram(dpu_alloc_wram_t *allocator, uint32_t mramIdx);

void store_wfa_cmpnt_to_mram(wfa_component *wfa, uint32_t mramIdx);
#endif