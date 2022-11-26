//
// Created by mhogan on 11/25/22.
//

#ifndef C_BOYER_MOORE_H
#define C_BOYER_MOORE_H


#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdlib.h>
#include <unistd.h>

#define ALPHABET_LEN 256
#define max(a, b) ((a < b) ? b : a)

// BAD CHARACTER RULE.
// delta1 table: delta1[c] contains the distance between the last
// character of pat and the rightmost occurrence of c in pat.
//
// If c does not occur in pat, then delta1[c] = patlen.
// If c is at string[i] and c != pat[patlen-1], we can safely shift i
//   over by delta1[c], which is the minimum distance needed to shift
//   pat forward to get string[i] lined up with some character in pat.
// c == pat[patlen-1] returning zero is only a concern for BMH, which
//   does not have delta2. BMH makes the value patlen in such a case.
//   We follow this choice instead of the original 0 because it skips
//   more. (correctness?)
//
// This algorithm runs in alphabet_len+patlen time.
void make_delta1(ptrdiff_t *delta1, uint8_t *pat, size_t patlen);


// true if the suffix of word starting from word[pos] is a prefix
// of word
bool is_prefix(uint8_t *word, size_t wordlen, ptrdiff_t pos);


// length of the longest suffix of word ending on word[pos].
// suffix_length("dddbcabc", 8, 4) = 2
size_t suffix_length(uint8_t *word, size_t wordlen, ptrdiff_t pos);


// GOOD SUFFIX RULE.
// delta2 table: given a mismatch at pat[pos], we want to align
// with the next possible full match could be based on what we
// know about pat[pos+1] to pat[patlen-1].
//
// In case 1:
// pat[pos+1] to pat[patlen-1] does not occur elsewhere in pat,
// the next plausible match starts at or after the mismatch.
// If, within the substring pat[pos+1 .. patlen-1], lies a prefix
// of pat, the next plausible match is here (if there are multiple
// prefixes in the substring, pick the longest). Otherwise, the
// next plausible match starts past the character aligned with
// pat[patlen-1].
//
// In case 2:
// pat[pos+1] to pat[patlen-1] does occur elsewhere in pat. The
// mismatch tells us that we are not looking at the end of a match.
// We may, however, be looking at the middle of a match.
//
// The first loop, which takes care of case 1, is analogous to
// the KMP table, adapted for a 'backwards' scan order with the
// additional restriction that the substrings it considers as
// potential prefixes are all suffixes. In the worst case scenario
// pat consists of the same letter repeated, so every suffix is
// a prefix. This loop alone is not sufficient, however:
// Suppose that pat is "ABYXCDBYX", and text is ".....ABYXCDEYX".
// We will match X, Y, and find B != E. There is no prefix of pat
// in the suffix "YX", so the first loop tells us to skip forward
// by 9 characters.
// Although superficially similar to the KMP table, the KMP table
// relies on information about the beginning of the partial match
// that the BM algorithm does not have.
//
// The second loop addresses case 2. Since suffix_length may not be
// unique, we want to take the minimum value, which will tell us
// how far away the closest potential match is.
void make_delta2(ptrdiff_t *delta2, uint8_t *pat, size_t patlen);


// Returns pointer to first match.
// See also glibc memmem() (non-BM) and std::boyer_moore_searcher (first-match).
uint8_t* boyer_moore(uint8_t *string, size_t stringlen, uint8_t *pat, size_t patlen);

#endif //C_BOYER_MOORE_H
