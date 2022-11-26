// Boyer-Moore string-matching algorithm in C17.
//
// Reads in a pattern P
// and a text T
// and returns the first substring, if any,
// that matches the pattern P
// Algorithms taken from Wikipedia page on Boyer-Moore

#include <stdio.h>
#include <string.h>
#include "boyer-moore.h"

#define MAX_STRING_LENGTH 1024

int main() {
    uint8_t pat[MAX_STRING_LENGTH];
    uint8_t tex[MAX_STRING_LENGTH];

    memset(pat, '\0', sizeof(uint8_t) * MAX_STRING_LENGTH );
    memset(tex, '\0', sizeof(uint8_t) * MAX_STRING_LENGTH );

    printf("Enter pattern P:\n");
    scanf("%s", pat);
    size_t pat_len = strlen((char*) pat);
    printf("Pattern P (length %lu) entered: %s \n",
           pat_len,
           pat);
    printf("Enter text T:\n");
    scanf("%s", tex);
    size_t tex_len = strlen((char*)tex);
    printf("Text T (length %lu) entered: %s\n",
           tex_len,
           tex);

    printf("Running Boyer-Moore\n");
    uint8_t* first_match = NULL;
    first_match = boyer_moore(
            tex,
            tex_len,
            pat,
            pat_len
            );
    printf("First match found is %s", first_match);
    return 0;
}
