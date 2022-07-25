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

#include "wfa_backtracing.h"
#include "../common/common.h"
/*
 * Backtrace Detect Limits
 */
bool affine_wavefronts_valid_location(
    int k,
    awf_offset_t offset,
    int pattern_length,
    int text_length)
{
  // Locate offset (remember that backtrace is always +1 offset ahead)
  int v = AFFINE_WAVEFRONT_V(k, offset);
  int h = AFFINE_WAVEFRONT_H(k, offset);
  return (v > 0 && v <= pattern_length &&
          h > 0 && h <= text_length);
}
void affine_wavefronts_offset_add_trailing_gap(
    edit_cigar_t *edit_cigar,
    int k,
    int alignment_k)
{
  // Parameters
  char *operations = edit_cigar->operations;
  int op_sentinel = edit_cigar->begin_offset;
  // Add trailing gap
  int i;
  if (k < alignment_k)
  {
    for (i = k; i < alignment_k; ++i)
      operations[op_sentinel--] = 'I';
  }
  else if (k > alignment_k)
  {
    for (i = alignment_k; i < k; ++i)
      operations[op_sentinel--] = 'D';
  }
  edit_cigar->begin_offset = op_sentinel;
}
/*
 * Backtrace Paths Offsets
 */
awf_offset_t backtrace_wavefront_trace_deletion_open_offset(
    wfa_component **wavefronts,
    int score,
    int k,
    awf_offset_t offset)
{
  if (score < 0)
    return AFFINE_WAVEFRONT_OFFSET_NULL;
  wfa_component *wfa = wavefronts[score];
  if (wfa != NULL &&
      wfa->klo <= k + 1 &&
      k + 1 <= wfa->khi)
  {
    return wfa->mwavefront[k + 1];
  }
  else
  {
    return AFFINE_WAVEFRONT_OFFSET_NULL;
  }
}
awf_offset_t backtrace_wavefront_trace_deletion_extend_offset(
    wfa_component **wavefronts,
    int score,
    int k,
    awf_offset_t offset)
{
  if (score < 0)
    return AFFINE_WAVEFRONT_OFFSET_NULL;
  wfa_component *wfa = wavefronts[score];
  if (wfa != NULL && !wfa->d_null &&
      wfa->klo <= k + 1 &&
      k + 1 <= wfa->khi)
  {
    return wfa->dwavefront[k + 1];
  }
  else
  {
    return AFFINE_WAVEFRONT_OFFSET_NULL;
  }
}
awf_offset_t backtrace_wavefront_trace_insertion_open_offset(
    wfa_component **wavefronts,
    int score,
    int k,
    awf_offset_t offset)
{
  if (score < 0)
    return AFFINE_WAVEFRONT_OFFSET_NULL;
  wfa_component *wfa = wavefronts[score];
  if (wfa != NULL &&
      wfa->klo <= k - 1 &&
      k - 1 <= wfa->khi)
  {
    return wfa->mwavefront[k - 1] + 1;
  }
  else
  {
    return AFFINE_WAVEFRONT_OFFSET_NULL;
  }
}
awf_offset_t backtrace_wavefront_trace_insertion_extend_offset(
    wfa_component **wavefronts,
    int score,
    int k,
    awf_offset_t offset)
{
  if (score < 0)
    return AFFINE_WAVEFRONT_OFFSET_NULL;
  wfa_component *wfa = wavefronts[score];
  if (wfa != NULL && wfa->iwavefront != NULL &&
      wfa->klo <= k - 1 &&
      k - 1 <= wfa->khi)
  {
    return wfa->iwavefront[k - 1] + 1;
  }
  else
  {
    return AFFINE_WAVEFRONT_OFFSET_NULL;
  }
}
awf_offset_t backtrace_wavefront_trace_mismatch_offset(
    wfa_component **wavefronts,
    int score,
    int k,
    awf_offset_t offset)
{
  if (score < 0)
    return AFFINE_WAVEFRONT_OFFSET_NULL;
  wfa_component *wfa = wavefronts[score];
  if (wfa != NULL &&
      wfa->klo <= k &&
      k <= wfa->khi)
  {
    return wfa->mwavefront[k] + 1;
  }
  else
  {
    return AFFINE_WAVEFRONT_OFFSET_NULL;
  }
}

/*
 * Backtrace Operations
 */
void affine_wavefronts_backtrace_matches__check(
    char *pattern,
    char *text,
    int k,
    awf_offset_t offset,
    bool valid_location,
    int num_matches,
    edit_cigar_t *edit_cigar)
{
  int i;
  for (i = 0; i < num_matches; ++i)
  {
    // Set Match
    edit_cigar->operations[(edit_cigar->begin_offset)--] = 'M';
    // Update state
    --offset;
  }
}

void affine_wavefronts_backtrace_matches(
    edit_cigar_t *edit_cigar,
    int num_matches)
{
  // Set Matches
  int i;
  for (i = 0; i < num_matches; ++i)
  {
    edit_cigar->operations[(edit_cigar->begin_offset)--] = 'M';
  }
}
/*
 * Backtrace (single solution)
 */
void affine_wavefronts_backtrace(
    wfa_component **affine_wavefronts,
    edit_cigar_t *cigar,
    char *pattern,
    int pattern_length,
    char *text,
    int text_length,
    int alignment_score)
{

  // Parameters
  int alignment_k = AFFINE_WAVEFRONT_DIAGONAL(text_length, pattern_length);

  // Compute starting location
  int score = alignment_score;
  int k = alignment_k;
  awf_offset_t offset = affine_wavefronts[alignment_score]->mwavefront[k];
  bool valid_location = affine_wavefronts_valid_location(k, offset, pattern_length, text_length);
  // Trace the alignment back
  backtrace_wavefront_type backtrace_type = backtrace_wavefront_M;
  int v = AFFINE_WAVEFRONT_V(k, offset);
  int h = AFFINE_WAVEFRONT_H(k, offset);
  while (v > 0 && h > 0 && score > 0)
  {
    // Check location
    if (!valid_location)
    {
      valid_location = affine_wavefronts_valid_location(k, offset, pattern_length, text_length);
      if (valid_location)
      {
        affine_wavefronts_offset_add_trailing_gap(cigar, k, alignment_k);
      }
    }
    // Compute scores
    int gap_open_score = score - GAP_O - GAP_E;
    int gap_extend_score = score - GAP_E;
    int mismatch_score = score - MISMATCH;
    // Compute source offsets
    awf_offset_t del_ext = (backtrace_type == backtrace_wavefront_I) ? AFFINE_WAVEFRONT_OFFSET_NULL : backtrace_wavefront_trace_deletion_extend_offset(affine_wavefronts, gap_extend_score, k, offset);
    awf_offset_t del_open = (backtrace_type == backtrace_wavefront_I) ? AFFINE_WAVEFRONT_OFFSET_NULL : backtrace_wavefront_trace_deletion_open_offset(affine_wavefronts, gap_open_score, k, offset);
    awf_offset_t ins_ext = (backtrace_type == backtrace_wavefront_D) ? AFFINE_WAVEFRONT_OFFSET_NULL : backtrace_wavefront_trace_insertion_extend_offset(affine_wavefronts, gap_extend_score, k, offset);
    awf_offset_t ins_open = (backtrace_type == backtrace_wavefront_D) ? AFFINE_WAVEFRONT_OFFSET_NULL : backtrace_wavefront_trace_insertion_open_offset(affine_wavefronts, gap_open_score, k, offset);
    awf_offset_t misms = (backtrace_type != backtrace_wavefront_M) ? AFFINE_WAVEFRONT_OFFSET_NULL : backtrace_wavefront_trace_mismatch_offset(affine_wavefronts, mismatch_score, k, offset);
    // Compute maximum offset
    awf_offset_t max_del = MAX(del_ext, del_open);
    awf_offset_t max_ins = MAX(ins_ext, ins_open);
    awf_offset_t max_all = MAX(misms, MAX(max_ins, max_del));
    // Traceback Matches
    if (backtrace_type == backtrace_wavefront_M)
    {
      int num_matches = offset - max_all;
      affine_wavefronts_backtrace_matches__check(pattern, text, k, offset, valid_location, num_matches, cigar);
      offset = max_all;
      // Update coordinates
      v = AFFINE_WAVEFRONT_V(k, offset);
      h = AFFINE_WAVEFRONT_H(k, offset);
      if (v <= 0 || h <= 0)
        break;
    }
    // Traceback Operation
    if (max_all == del_ext)
    {
      // Add Deletion
      if (valid_location)
        cigar->operations[(cigar->begin_offset)--] = 'D';
      // Update state
      score = gap_extend_score;
      ++k;
      backtrace_type = backtrace_wavefront_D;
    }
    else if (max_all == del_open)
    {
      // Add Deletion
      if (valid_location)
        cigar->operations[(cigar->begin_offset)--] = 'D';
      // Update state
      score = gap_open_score;
      ++k;
      backtrace_type = backtrace_wavefront_M;
    }
    else if (max_all == ins_ext)
    {
      // Add Insertion
      if (valid_location)
        cigar->operations[(cigar->begin_offset)--] = 'I';
      // Update state
      score = gap_extend_score;
      --k;
      --offset;
      backtrace_type = backtrace_wavefront_I;
    }
    else if (max_all == ins_open)
    {
      // Add Insertion
      if (valid_location)
        cigar->operations[(cigar->begin_offset)--] = 'I';
      // Update state
      score = gap_open_score;
      --k;
      --offset;
      backtrace_type = backtrace_wavefront_M;
    }
    else if (max_all == misms)
    {
      // Add Mismatch
      if (valid_location)
        cigar->operations[(cigar->begin_offset)--] = 'X';
      // Update state
      score = mismatch_score;
      --offset;
    }
    else
    {
      printf("Backtrace error: No link found during backtrace\n");
      exit(1);
    }
    // Update coordinates
    v = AFFINE_WAVEFRONT_V(k, offset);
    h = AFFINE_WAVEFRONT_H(k, offset);
  }
  // Account for last operations
  if (score == 0)
  {
    // Account for last stroke of matches
    affine_wavefronts_backtrace_matches__check(pattern, text, k, offset, valid_location, offset, cigar);
  }
  else
  {
    // Account for last stroke of insertion/deletion
    while (v > 0)
    {
      cigar->operations[(cigar->begin_offset)--] = 'D';
      --v;
    };
    while (h > 0)
    {
      cigar->operations[(cigar->begin_offset)--] = 'I';
      --h;
    };
  }
  ++(cigar->begin_offset); // Set CIGAR length
}
