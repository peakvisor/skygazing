#ifndef SKYGAZING_TESTS_H
#define SKYGAZING_TESTS_H
#include <random>
#include <chrono>
#include <set>

#include "skygazing_sun.h"
#include "skygazing_moon.h"
#include "skygazing_test_analytics.h"
#include "skygazing_test_cases.h"

namespace Skygazing {

struct Testing {
    static void testTimeTransforms() {
        std::random_device true_random{};
        RandomDataGenerator dataGenerator(true_random());

        const int kIterations = 1000;
        for (int year = 2017; year < 2030; ++year) {
            for (int i = 0; i < kIterations; ++i) {
                auto date = dataGenerator.getRandomDate(year);
                Seconds seconds = tmToSeconds(&date);
                std::string dateString = secondsToDateString(seconds);
                assert(dateString == secondsToDateString(seconds));
                Julian terrestrialJulian = secondsToTerrestrialJulian(seconds);
                SKYGAZING_ASSERT_SMALL(seconds - terrestrialJulianToSeconds(terrestrialJulian), 1e-5);
            }
        }
    }

    static void testSunTimes(bool showStatistics = false) {
        StatisticsAggregator statistics;
        std::vector<SolarNoonTests::TestCase> transitTestCases(SolarNoonTests::sunTransits.begin(), SolarNoonTests::sunTransits.end());
        std::sort(transitTestCases.begin(), transitTestCases.end(), [](const auto &lhs, const auto &rhs) {
            return dateStringToSeconds(lhs.dateString) < dateStringToSeconds(rhs.dateString);
        });
        for (const auto &test : transitTestCases) {
            auto seconds = dateStringToSeconds(test.dateString, test.timezone);
            auto transit = Sky::getTransit<Sun>(seconds, test.coordinates);
            Seconds sunTransitDiff = (transit.value() - seconds);
            SKYGAZING_ASSERT_SMALL(sunTransitDiff, 60);

            if (showStatistics) {
                statistics.account("sunTransitDiff", sunTransitDiff);
            }
        }

        std::vector<SolarNoonTests::TestCase> sunriseTestCases(SolarNoonTests::sunrises.begin(), SolarNoonTests::sunrises.end());
        std::sort(sunriseTestCases.begin(), sunriseTestCases.end(), [](const auto &lhs, const auto &rhs) {
            return dateStringToSeconds(lhs.dateString) < dateStringToSeconds(rhs.dateString);
        });

        for (const auto &test : sunriseTestCases) {
            auto seconds = dateStringToSeconds(test.dateString, test.timezone);
            auto times = Sky::findTimes<Sun>(seconds, test.coordinates, radsFromDegrees(0.5));
            Seconds sunriseTimeDiff = (times.rise.value() - seconds);
            SKYGAZING_ASSERT_SMALL(sunriseTimeDiff, 90);

            if (showStatistics) {
                statistics.account("sunriseTimeDiff", sunriseTimeDiff);
            }
        }
    }

    static void testMoonTimes(bool showStatistics = false) {
        StatisticsAggregator statistics;
        std::vector<SolarNoonTests::TestCase> transitTestCases(SolarNoonTests::moonrises.begin(), SolarNoonTests::moonrises.end());
        std::sort(transitTestCases.begin(), transitTestCases.end(), [](const auto &lhs, const auto &rhs) {
            return dateStringToSeconds(lhs.dateString) < dateStringToSeconds(rhs.dateString);
        });
        Rads bodySize = -radsFromDegrees(1.15);
        for (const auto &test : transitTestCases) {
            auto seconds = dateStringToSeconds(test.dateString, test.timezone);
            auto times = Sky::findTimes<Moon>(seconds, test.coordinates, bodySize);
            Seconds moonriseTimeDiff = times.rise.value() - seconds;
            SKYGAZING_ASSERT_SMALL(moonriseTimeDiff, 60);

            if (showStatistics) {
                statistics.account("moonriseTimeDiff", moonriseTimeDiff);
            }
        }

        std::vector<SolarNoonTests::TestCase> msetTestCases(SolarNoonTests::moonsets.begin(), SolarNoonTests::moonsets.end());
        std::sort(msetTestCases.begin(), msetTestCases.end(), [](const auto &lhs, const auto &rhs) {
            return dateStringToSeconds(lhs.dateString) < dateStringToSeconds(rhs.dateString);
        });
        for (const auto &test : msetTestCases) {
            auto seconds = dateStringToSeconds(test.dateString, test.timezone);
            auto times = Sky::findTimes<Moon>(seconds, test.coordinates, bodySize);
            Seconds moonsetTimeDiff = times.set.value() - seconds;
            SKYGAZING_ASSERT_SMALL(moonsetTimeDiff, 90);

            if (showStatistics) {
                statistics.account("moonsetTimeDiff", moonsetTimeDiff);
            }
        }
    }

    static void runAllTests() {
        const auto startTime = std::chrono::high_resolution_clock::now();

        testTimeTransforms();
        const auto timeTransformsTime = std::chrono::high_resolution_clock::now();
        std::cout << std::fixed << "tested timeTransformsTime in "
                  << std::chrono::duration<double, std::milli>(timeTransformsTime - startTime).count() / 1000. << "s\n";
        testSunTimes();
        const auto sunTimesTime = std::chrono::high_resolution_clock::now();
        std::cout << std::fixed << "tested sunTimes in "
                  << std::chrono::duration<double, std::milli>(sunTimesTime - timeTransformsTime).count() / 1000. << "s\n";
        testMoonTimes();
        const auto moonTimesTime = std::chrono::high_resolution_clock::now();
        std::cout << std::fixed << "tested moonTimes in "
                  << std::chrono::duration<double, std::milli>(moonTimesTime - sunTimesTime).count() / 1000. << "s\n";
    }
};

}

#endif // SKYGAZING_TESTS_H
