#ifndef SKYGAZING_TESTS_H
#define SKYGAZING_TESTS_H

#undef NDEBUG
#include <cassert>
#include <random>
#include <chrono>
#include <set>

#include "../skygazing.h"
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
                std::string dateString = dateStringFromSeconds(seconds);
                assert(dateString == dateStringFromSeconds(seconds));
                Julian terrestrialJulian = terrestrialJulianFromSeconds(seconds);
                SKYGAZING_ASSERT_SMALL(seconds - secondsFromTerrestrialJulian(terrestrialJulian), 1e-5);
            }
        }
    }

    template <typename CelestialBody>
    static void testTransitExistence() {
        std::random_device true_random{};
        RandomDataGenerator dataGenerator(true_random());

        const int kIterations = 100;
        for (int year = 2017; year < 2030; ++year) {
            for (int i = 0; i < kIterations; ++i) {
                auto date = dataGenerator.getRandomDate(year);
                Seconds seconds = tmToSeconds(&date);
                auto coordinates = dataGenerator.getRandomCoordinates(0, 85);
                if (!Sky::getTransit<CelestialBody>(seconds, coordinates).has_value()) {
                    std::cout << dateStringFromSeconds(seconds) << std::endl;
                    std::cout << coordinates.lat << ", " << coordinates.lng << std::endl;
                    abort();
                }
            }
        }
    }

    template <typename Func, typename Cases>
    static void testTimes(const Func &func, const Cases &cases, const std::string &description, Seconds acceptableDiff,
                          bool showStatistics = false) {
        StatisticsAggregator statistics;
        for (const auto &test : cases) {
            auto seconds = secondsFromDateString(test.dateString, test.timezone);
            auto result = func(seconds, test.coordinates);
            assert(result.has_value());
            Seconds diff = (*result - seconds);
            SKYGAZING_ASSERT_SMALL(diff, acceptableDiff);

            if (showStatistics) {
                statistics.account(description + " diff ", diff);
            }
        }
    }



    static void runAllTests() {
        SKYGAZING_TIME_EXECUTION(testTimes(Sky::getRise<Sun>, TimesTestCases::sunrises, "sunrises", 120));
        SKYGAZING_TIME_EXECUTION(testTimes(Sky::getTransit<Sun>, TimesTestCases::sunTransits, "transits", 60));
        SKYGAZING_TIME_EXECUTION(testTimes(Sky::getSet<Sun>, TimesTestCases::sunsets, "sunsets", 90));
        SKYGAZING_TIME_EXECUTION(testTimes(Sky::getRise<Moon>, TimesTestCases::moonrises, "moonrises", 60));
        SKYGAZING_TIME_EXECUTION(testTimes(Sky::getSet<Moon>, TimesTestCases::moonsets, "moonsets", 90));

        // check that moon transit exist
        Seconds kNoPrecision = std::numeric_limits<Seconds>::max();
        SKYGAZING_TIME_EXECUTION(testTimes(Sky::getTransit<Moon>, TimesTestCases::moonrises, "moonrises", kNoPrecision));
        SKYGAZING_TIME_EXECUTION(testTimes(Sky::getTransit<Moon>, TimesTestCases::moonsets, "moonrises", kNoPrecision));
        SKYGAZING_TIME_EXECUTION(testTimes(Sky::getTransit<Moon>, TimesTestCases::sunTransits, "moonrises", kNoPrecision));
        SKYGAZING_TIME_EXECUTION(testTransitExistence<Moon>());
    }
};

}

#endif // SKYGAZING_TESTS_H
