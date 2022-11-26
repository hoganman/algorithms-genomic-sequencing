#include <iostream>

#include "boyer-moore.h"

int main() {
    std::string pat, tex;

    std::cout << "Enter pattern P:" << std::endl;
    std::cin >> pat;
    auto pat_len = pat.length();
    printf("Pattern P (length %lu) entered: %s \n",
           pat_len,
           pat.c_str());

    std::cout << "Enter text T:" << std::endl;
    std::cin >> tex;
    auto tex_len = tex.length();
    printf("Text T (length %lu) entered: %s\n",
           tex_len,
           tex.c_str());

    std::cout << "Running Boyer-Moore" << std::endl;
    auto first_match_offset = boyer_moore(
            pat,
            tex
            );
    std::string first_match = tex.substr(first_match_offset);
    printf("First match found is at %s", first_match.c_str());

    return 0;
}
