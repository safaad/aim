#ifndef CUSTOM_MEM_H_
#define CUSTOM_MEM_H_

#include <defs.h>
#include <mram.h>
#include <alloc.h>
typedef struct dpu_alloc_wram_t
{
    unsigned int segment_size;
    char *HEAD_PTR_WRAM;
    char *CUR_PTR_WRAM;
    uint32_t mem_used_wram;
} dpu_alloc_wram_t;

// allocate the first segment in the WRAM
dpu_alloc_wram_t init_dpu_alloc_wram(unsigned int segment_size);

char *allocate_new(dpu_alloc_wram_t *dpu_alloc_obj, unsigned int size);
void reset_dpu_alloc_wram(dpu_alloc_wram_t *dpu_alloc_obj);

#endif