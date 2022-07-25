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

#endif