#include "dpu_allocator_wram.h"

dpu_alloc_wram_t init_dpu_alloc_wram(unsigned int segment_size)
{
    dpu_alloc_wram_t dpu_alloc_obj;
    if (segment_size * NR_TASKLETS >= 62000)
    {
        printf("Out of WRAM memory\n");
        exit(1);
    }
    dpu_alloc_obj.mem_used_wram = 0;
    dpu_alloc_obj.segment_size = ROUND_UP_MULTIPLE_8(segment_size);
    dpu_alloc_obj.HEAD_PTR_WRAM = (char *)mem_alloc(segment_size);
    dpu_alloc_obj.CUR_PTR_WRAM = dpu_alloc_obj.HEAD_PTR_WRAM;
    return dpu_alloc_obj;
}

char *allocate_new(dpu_alloc_wram_t *dpu_alloc_obj, unsigned int size)
{
    if (size <= 0)
        return NULL;
    if (((ROUND_UP_MULTIPLE_8(size) + dpu_alloc_obj->mem_used_wram) >= dpu_alloc_obj->segment_size))
    {
        printf("Out of WRAM memory %d %d\n", size, (ROUND_UP_MULTIPLE_8(size)) + dpu_alloc_obj->mem_used_wram);
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
