#include <iostream>
#include <chrono>

#include "skygazing_test.h"

int main() {
    const auto t1 = std::chrono::high_resolution_clock::now();
    Skygazing::Testing::runAllTests();
    const auto t2 = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double, std::milli> ms = t2 - t1;
    std::cout << std::fixed << "finished in " << ms.count() / 1000. << "s\n";
    return 0;
}

