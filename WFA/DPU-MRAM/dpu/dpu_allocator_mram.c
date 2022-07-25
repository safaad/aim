#include "dpu_allocator_mram.h"

void add_wfa_cmpnt_to_mram(uint32_t *mramIdx, uint32_t size, dpu_alloc_mram_t *dpu_alloc_mram)
{

    if (((ROUND_UP_MULTIPLE_8(size) + dpu_alloc_mram->mem_used_mram) >= dpu_alloc_mram->segment_size))
    {
        printf("Out of memory MRAM\n");
        exit(-1);
    }
    if (size <= 0)
    {
        *(mramIdx) = 0;
        return;
    }
    size = ROUND_UP_MULTIPLE_8(size);
    dpu_alloc_mram->mem_used_mram += size;
    *(mramIdx) = dpu_alloc_mram->CUR_PTR_MRAM;
    dpu_alloc_mram->CUR_PTR_MRAM += size;
}

wfa_component *load_wfa_cmpnt_from_mram(dpu_alloc_wram_t *allocator, uint32_t mramIdx)
{
    if (mramIdx == 0)
    {
        return NULL;
    }
    wfa_component *wfa = (wfa_component *)allocate_new(allocator, sizeof(wfa_component));
    uint32_t wfa_m = ((uint32_t)DPU_MRAM_HEAP_POINTER) + mramIdx;
    mram_read((__mram_ptr void const *)(wfa_m), wfa, ROUND_UP_MULTIPLE_8(sizeof(wfa_component)));
    int wv_len = wfa->hi_base - wfa->lo_base + 1;
    wfa_m += ROUND_UP_MULTIPLE_8(sizeof(wfa_component));
    awf_offset_t *moffset = (awf_offset_t *)allocate_new(allocator, wv_len * sizeof(awf_offset_t));
    if (ROUND_UP_MULTIPLE_8(wv_len * sizeof(awf_offset_t)) <= 2048)
    {
        mram_read((__mram_ptr void const *)(wfa_m), moffset, ROUND_UP_MULTIPLE_8(wv_len * sizeof(awf_offset_t)));
    }
    else
    {
        int wavefront_segment = 2048 / sizeof(awf_offset_t);
        for (int num_wavefronts = 0; num_wavefronts <= wv_len; num_wavefronts += wavefront_segment)
        {
            if (num_wavefronts + wavefront_segment <= wv_len)
            {
                mram_read((__mram_ptr void const *)(wfa_m + num_wavefronts * sizeof(awf_offset_t)), &moffset[num_wavefronts], (wavefront_segment * sizeof(awf_offset_t)));
            }
            else
            {
                int wv_remaining = wv_len - num_wavefronts;
                if (wv_remaining > 0)
                {
                    mram_read((__mram_ptr void const *)(wfa_m + num_wavefronts * sizeof(awf_offset_t)), &moffset[num_wavefronts], ROUND_UP_MULTIPLE_8(wv_remaining * sizeof(awf_offset_t)));
                }
            }
        }
    }

    wfa->mwavefront = (awf_offset_t *)(moffset - wfa->lo_base);
    wfa_m += ROUND_UP_MULTIPLE_8(wv_len * sizeof(awf_offset_t));
    if (wfa->i_null)
    {
        wfa->iwavefront = NULL;
    }
    else
    {
        awf_offset_t *ioffset = (awf_offset_t *)allocate_new(allocator, wv_len * sizeof(awf_offset_t));
        if (ROUND_UP_MULTIPLE_8(wv_len * sizeof(awf_offset_t)) <= 2048)
        {
            mram_read((__mram_ptr void const *)(wfa_m), ioffset, ROUND_UP_MULTIPLE_8(wv_len * sizeof(awf_offset_t)));
        }
        else
        {
            int wavefront_segment = 2048 / sizeof(awf_offset_t);
            for (int num_wavefronts = 0; num_wavefronts <= wv_len; num_wavefronts += wavefront_segment)
            {
                if (num_wavefronts + wavefront_segment <= wv_len)
                {
                    mram_read((__mram_ptr void const *)(wfa_m + num_wavefronts * sizeof(awf_offset_t)), &ioffset[num_wavefronts], (wavefront_segment * sizeof(awf_offset_t)));
                }
                else
                {
                    int wv_remaining = wv_len - num_wavefronts;
                    if (wv_remaining > 0)
                    {
                        mram_read((__mram_ptr void const *)(wfa_m + num_wavefronts * sizeof(awf_offset_t)), &ioffset[num_wavefronts], ROUND_UP_MULTIPLE_8(wv_remaining * sizeof(awf_offset_t)));
                    }
                }
            }
        }
        wfa->iwavefront = (awf_offset_t *)(ioffset - wfa->lo_base);
        wfa_m += ROUND_UP_MULTIPLE_8(wv_len * sizeof(awf_offset_t));
    }
    if (wfa->d_null)
    {
        wfa->dwavefront = NULL;
    }
    else
    {
        awf_offset_t *doffset = (awf_offset_t *)allocate_new(allocator, wv_len * sizeof(awf_offset_t));
        if (ROUND_UP_MULTIPLE_8(wv_len * sizeof(awf_offset_t)) <= 2048)
        {
            mram_read((__mram_ptr void const *)(wfa_m), doffset, ROUND_UP_MULTIPLE_8(wv_len * sizeof(awf_offset_t)));
        }
        else
        {
            int wavefront_segment = 2048 / sizeof(awf_offset_t);
            for (int num_wavefronts = 0; num_wavefronts <= wv_len; num_wavefronts += wavefront_segment)
            {
                if (num_wavefronts + wavefront_segment <= wv_len)
                {
                    mram_read((__mram_ptr void const *)(wfa_m + num_wavefronts * sizeof(awf_offset_t)), &doffset[num_wavefronts], (wavefront_segment * sizeof(awf_offset_t)));
                }
                else
                {
                    int wv_remaining = wv_len - num_wavefronts;
                    if (wv_remaining > 0)
                    {
                        mram_read((__mram_ptr void const *)(wfa_m + num_wavefronts * sizeof(awf_offset_t)), &doffset[num_wavefronts], ROUND_UP_MULTIPLE_8(wv_remaining * sizeof(awf_offset_t)));
                    }
                }
            }
        }
        wfa->dwavefront = (awf_offset_t *)(doffset - wfa->lo_base);
    }
    return wfa;
}

void store_wfa_cmpnt_to_mram(wfa_component *wfa, uint32_t mramIdx)
{
    if (wfa == NULL || mramIdx == 0)
    {
        return;
    }
    uint32_t wfa_m = ((uint32_t)DPU_MRAM_HEAP_POINTER) + mramIdx;
    mram_write(wfa, (__mram_ptr void *)(wfa_m), ROUND_UP_MULTIPLE_8(sizeof(wfa_component)));
    wfa_m += ROUND_UP_MULTIPLE_8(sizeof(wfa_component));

    int wv_len = wfa->hi_base - wfa->lo_base + 1;
    awf_offset_t *moffset = (awf_offset_t *)(wfa->mwavefront + wfa->lo_base);

    if (ROUND_UP_MULTIPLE_8(wv_len * sizeof(awf_offset_t)) <= 2048)
    {
        mram_write(moffset, (__mram_ptr void *)(wfa_m), ROUND_UP_MULTIPLE_8(wv_len * sizeof(awf_offset_t)));
    }
    else
    {
        int wavefront_segment = 2048 / sizeof(awf_offset_t);
        for (int num_wavefronts = 0; num_wavefronts <= wv_len; num_wavefronts += wavefront_segment)
        {
            if (num_wavefronts + wavefront_segment <= wv_len)
            {
                mram_write(&moffset[num_wavefronts], (__mram_ptr void *)(wfa_m + num_wavefronts * sizeof(awf_offset_t)), (wavefront_segment * sizeof(awf_offset_t)));
            }
            else
            {
                int wv_remaining = wv_len - num_wavefronts;
                if (wv_remaining > 0)
                {
                    mram_write(&moffset[num_wavefronts], (__mram_ptr void *)(wfa_m + num_wavefronts * sizeof(awf_offset_t)), ROUND_UP_MULTIPLE_8(wv_remaining * sizeof(awf_offset_t)));
                }
            }
        }
    }
    wfa_m += ROUND_UP_MULTIPLE_8(wv_len * sizeof(awf_offset_t));
    if (!wfa->i_null)
    {
        awf_offset_t *ioffset = (awf_offset_t *)(wfa->iwavefront + wfa->lo_base);
        if (ROUND_UP_MULTIPLE_8(wv_len * sizeof(awf_offset_t)) <= 2048)
        {
            mram_write(ioffset, (__mram_ptr void *)(wfa_m), ROUND_UP_MULTIPLE_8(wv_len * sizeof(awf_offset_t)));
        }
        else
        {
            int wavefront_segment = 2048 / sizeof(awf_offset_t);
            for (int num_wavefronts = 0; num_wavefronts <= wv_len; num_wavefronts += wavefront_segment)
            {
                if (num_wavefronts + wavefront_segment <= wv_len)
                {
                    mram_write(&ioffset[num_wavefronts], (__mram_ptr void *)(wfa_m + num_wavefronts * sizeof(awf_offset_t)), (wavefront_segment * sizeof(awf_offset_t)));
                }
                else
                {
                    int wv_remaining = wv_len - num_wavefronts;
                    if (wv_remaining > 0)
                    {
                        mram_write(&ioffset[num_wavefronts], (__mram_ptr void *)(wfa_m + num_wavefronts * sizeof(awf_offset_t)), ROUND_UP_MULTIPLE_8(wv_remaining * sizeof(awf_offset_t)));
                    }
                }
            }
        }
        wfa_m += ROUND_UP_MULTIPLE_8(wv_len * sizeof(awf_offset_t));
    }

    if (!wfa->d_null)
    {
        awf_offset_t *doffset = (awf_offset_t *)(wfa->dwavefront + wfa->lo_base);
        if (ROUND_UP_MULTIPLE_8(wv_len * sizeof(awf_offset_t)) <= 2048)
        {
            mram_write(doffset, (__mram_ptr void *)(wfa_m), ROUND_UP_MULTIPLE_8(wv_len * sizeof(awf_offset_t)));
        }
        else
        {
            int wavefront_segment = 2048 / sizeof(awf_offset_t);
            for (int num_wavefronts = 0; num_wavefronts <= wv_len; num_wavefronts += wavefront_segment)
            {
                if (num_wavefronts + wavefront_segment <= wv_len)
                {
                    mram_write(&doffset[num_wavefronts], (__mram_ptr void *)(wfa_m + num_wavefronts * sizeof(awf_offset_t)), (wavefront_segment * sizeof(awf_offset_t)));
                }
                else
                {
                    int wv_remaining = wv_len - num_wavefronts;
                    if (wv_remaining > 0)
                    {
                        mram_write(&doffset[num_wavefronts], (__mram_ptr void *)(wfa_m + num_wavefronts * sizeof(awf_offset_t)), ROUND_UP_MULTIPLE_8(wv_remaining * sizeof(awf_offset_t)));
                    }
                }
            }
        }
    }
}

wfa_component *load_mwavefront_cmpnt_from_mram(dpu_alloc_wram_t *allocator, uint32_t mramIdx)
{
    if (mramIdx == 0)
    {
        return NULL;
    }
    wfa_component *wfa = (wfa_component *)allocate_new(allocator, sizeof(wfa_component));
    uint32_t wfa_m = ((uint32_t)DPU_MRAM_HEAP_POINTER) + mramIdx;
    mram_read((__mram_ptr void const *)(wfa_m), wfa, ROUND_UP_MULTIPLE_8(sizeof(wfa_component)));
    int wv_len = wfa->hi_base - wfa->lo_base + 1;
    wfa_m += ROUND_UP_MULTIPLE_8(sizeof(wfa_component));
    awf_offset_t *moffset = (awf_offset_t *)allocate_new(allocator, wv_len * sizeof(awf_offset_t));
    if (ROUND_UP_MULTIPLE_8(wv_len * sizeof(awf_offset_t)) <= 2048)
    {
        mram_read((__mram_ptr void const *)(wfa_m), moffset, ROUND_UP_MULTIPLE_8(wv_len * sizeof(awf_offset_t)));
    }
    else
    {
        int wavefront_segment = 2048 / sizeof(awf_offset_t);
        for (int num_wavefronts = 0; num_wavefronts <= wv_len; num_wavefronts += wavefront_segment)
        {
            if (num_wavefronts + wavefront_segment <= wv_len)
            {
                mram_read((__mram_ptr void const *)(wfa_m + num_wavefronts * sizeof(awf_offset_t)), &moffset[num_wavefronts], (wavefront_segment * sizeof(awf_offset_t)));
            }
            else
            {
                int wv_remaining = wv_len - num_wavefronts;
                if (wv_remaining > 0)
                {
                    mram_read((__mram_ptr void const *)(wfa_m + num_wavefronts * sizeof(awf_offset_t)), &moffset[num_wavefronts], ROUND_UP_MULTIPLE_8(wv_remaining * sizeof(awf_offset_t)));
                }
            }
        }
    }
    wfa->mwavefront = (awf_offset_t *)(moffset - wfa->lo_base);
    wfa_m += ROUND_UP_MULTIPLE_8(wv_len * sizeof(awf_offset_t));

    wfa->iwavefront = NULL;

    wfa->dwavefront = NULL;

    return wfa;
}

wfa_component *load_idwavefront_cmpnt_from_mram(dpu_alloc_wram_t *allocator, uint32_t mramIdx)
{
    if (mramIdx == 0)
    {
        return NULL;
    }
    wfa_component *wfa = (wfa_component *)allocate_new(allocator, sizeof(wfa_component));
    uint32_t wfa_m = ((uint32_t)DPU_MRAM_HEAP_POINTER) + mramIdx;
    mram_read((__mram_ptr void const *)(wfa_m), wfa, ROUND_UP_MULTIPLE_8(sizeof(wfa_component)));
    int wv_len = wfa->hi_base - wfa->lo_base + 1;
    wfa_m += ROUND_UP_MULTIPLE_8(sizeof(wfa_component));
    wfa->mwavefront = NULL;
    wfa_m += ROUND_UP_MULTIPLE_8(wv_len * sizeof(awf_offset_t));
    if (wfa->i_null)
    {
        wfa->iwavefront = NULL;
    }
    else
    {
        awf_offset_t *ioffset = (awf_offset_t *)allocate_new(allocator, wv_len * sizeof(awf_offset_t));
        if (ROUND_UP_MULTIPLE_8(wv_len * sizeof(awf_offset_t)) <= 2048)
        {
            mram_read((__mram_ptr void const *)(wfa_m), ioffset, ROUND_UP_MULTIPLE_8(wv_len * sizeof(awf_offset_t)));
        }
        else
        {
            int wavefront_segment = 2048 / sizeof(awf_offset_t);
            for (int num_wavefronts = 0; num_wavefronts <= wv_len; num_wavefronts += wavefront_segment)
            {
                if (num_wavefronts + wavefront_segment <= wv_len)
                {
                    mram_read((__mram_ptr void const *)(wfa_m + num_wavefronts * sizeof(awf_offset_t)), &ioffset[num_wavefronts], (wavefront_segment * sizeof(awf_offset_t)));
                }
                else
                {
                    int wv_remaining = wv_len - num_wavefronts;
                    if (wv_remaining > 0)
                    {
                        mram_read((__mram_ptr void const *)(wfa_m + num_wavefronts * sizeof(awf_offset_t)), &ioffset[num_wavefronts], ROUND_UP_MULTIPLE_8(wv_remaining * sizeof(awf_offset_t)));
                    }
                }
            }
        }
        wfa->iwavefront = (awf_offset_t *)(ioffset - wfa->lo_base);
        wfa_m += ROUND_UP_MULTIPLE_8(wv_len * sizeof(awf_offset_t));
    }
    if (wfa->d_null)
    {
        wfa->dwavefront = NULL;
    }
    else
    {
        awf_offset_t *doffset = (awf_offset_t *)allocate_new(allocator, wv_len * sizeof(awf_offset_t));
        if (ROUND_UP_MULTIPLE_8(wv_len * sizeof(awf_offset_t)) <= 2048)
        {
            mram_read((__mram_ptr void const *)(wfa_m), doffset, ROUND_UP_MULTIPLE_8(wv_len * sizeof(awf_offset_t)));
        }
        else
        {
            int wavefront_segment = 2048 / sizeof(awf_offset_t);
            for (int num_wavefronts = 0; num_wavefronts <= wv_len; num_wavefronts += wavefront_segment)
            {
                if (num_wavefronts + wavefront_segment <= wv_len)
                {
                    mram_read((__mram_ptr void const *)(wfa_m + num_wavefronts * sizeof(awf_offset_t)), &doffset[num_wavefronts], (wavefront_segment * sizeof(awf_offset_t)));
                }
                else
                {
                    int wv_remaining = wv_len - num_wavefronts;
                    if (wv_remaining > 0)
                    {
                        mram_read((__mram_ptr void const *)(wfa_m + num_wavefronts * sizeof(awf_offset_t)), &doffset[num_wavefronts], ROUND_UP_MULTIPLE_8(wv_remaining * sizeof(awf_offset_t)));
                    }
                }
            }
        }
        wfa->dwavefront = (awf_offset_t *)(doffset - wfa->lo_base);
    }
    return wfa;
}