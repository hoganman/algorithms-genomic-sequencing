//
// Created by mhogan on 11/25/22.
//

#ifndef CPLUSPLUS_BOYER_MOORE_H
#define CPLUSPLUS_BOYER_MOORE_H

#include <string>

// Return the first offset of the matching pattern against reference text
unsigned long long boyer_moore(std::string &pat, std::string &tex);

#endif //CPLUSPLUS_BOYER_MOORE_H
