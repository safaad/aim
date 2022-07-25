#include "dpu_allocator_wram.h"
#include <assert.h>
#include "common.h"

// allocate the first segment in the WRAM
dpu_alloc_wram_t init_dpu_alloc_wram(unsigned int segment_size)
{
    dpu_alloc_wram_t dpu_alloc_obj;
    dpu_alloc_obj.mem_used_wram = 0;
    // check
    dpu_alloc_obj.segment_size = ROUND_UP_MULTIPLE_8(segment_size);
    dpu_alloc_obj.HEAD_PTR_WRAM = (char *)mem_alloc(segment_size);
    dpu_alloc_obj.CUR_PTR_WRAM = dpu_alloc_obj.HEAD_PTR_WRAM;
    return dpu_alloc_obj;
}

char *allocate_new(dpu_alloc_wram_t *dpu_alloc_obj, unsigned int size)
{
    if (((ROUND_UP_MULTIPLE_8(size) + dpu_alloc_obj->mem_used_wram) >= dpu_alloc_obj->segment_size))
    {
        printf("out of memory\n");
        exit(1);
    }
    size = ROUND_UP_MULTIPLE_8(size);
    dpu_alloc_obj->mem_used_wram += size;
    char *allocated = (char *)dpu_alloc_obj->CUR_PTR_WRAM;
    dpu_alloc_obj->CUR_PTR_WRAM += size;

    return allocated;
}

void reset_dpu_alloc_wram(dpu_alloc_wram_t *dpu_alloc_obj)
{
    dpu_alloc_obj->mem_used_wram = 0;
    dpu_alloc_obj->CUR_PTR_WRAM = dpu_alloc_obj->HEAD_PTR_WRAM;
}
