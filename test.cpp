#include <type_traits>
#include <algorithm>
#include <sstream>
#include <array>
#include <numeric>
#include <limits>
#include "tiny_test.hpp"
#include "matrix.h"

using testing::make_pretty_test;
using testing::TestGroup;


template<typename T>
concept IMatrix = requires(T matrix) {

};

TestGroup all_tests[] = {
    {"Arithmetic",
        make_pretty_test("static tests", [](auto& test) {
               
        })
    }
};


int main() {
    bool success = true;
    for (auto& group : all_tests) {
        success &= group.run();
    }
    return success ? 0 : 1;
}

