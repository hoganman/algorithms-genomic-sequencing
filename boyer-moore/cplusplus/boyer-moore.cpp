//
// Created by mhogan on 11/25/22.
//

#include "boyer-moore.h"
#include <iomanip>
#include <algorithm>
#include <functional>

unsigned long long boyer_moore(std::string &pat, std::string &tex){
    const auto it = std::search(
            tex.begin(),
            tex.end(),
            std::boyer_moore_searcher(
                    pat.begin(),
                    pat.end()
                    ));
    unsigned long long match = it - tex.begin();
    return match;
}
