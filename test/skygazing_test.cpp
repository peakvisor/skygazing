#include <iostream>
#include <chrono>

#include "skygazing_test.h"

int main(int argc, char **argv) {
    const auto t1 = std::chrono::high_resolution_clock::now();
    bool runAnalytics = argc > 1 && argv[1] == "analytics"s;
    Skygazing::Testing::runAllTests(runAnalytics);
    const auto t2 = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double, std::milli> ms = t2 - t1;
    std::cout << std::fixed << "finished in " << ms.count() / 1000. << "s\n";
    return 0;
}

