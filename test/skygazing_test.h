#ifndef SKYGAZING_TESTS_H
#define SKYGAZING_TESTS_H

#undef NDEBUG
#include <cassert>
#include <random>
#include <chrono>
#include <set>

#include "../include/skygazing.h"
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
                Seconds seconds = tmToUTC(&date);
                std::string dateString = dateStringFromUTC(seconds);
                assert(dateString == dateStringFromUTC(seconds));
                Julian terrestrialJulian = ttFromUTC(seconds);
                SKYGAZING_ASSERT_SMALL(seconds - utcFromTT(terrestrialJulian), 1e-5);
            }
        }
    }

    template <typename CelestialBody>
    static void testTransitExistence(int iterationsPerYear = 1000) {
        std::random_device true_random{};
        RandomDataGenerator dataGenerator(true_random());

        for (int year = 2017; year < 2030; ++year) {
            for (int i = 0; i < iterationsPerYear; ++i) {
                auto date = dataGenerator.getRandomDate(year);
                Seconds seconds = tmToUTC(&date);
                auto coordinates = dataGenerator.getRandomCoordinates(0, 89);
                Sky::getTrajectoryCycleFromUTC<CelestialBody>(seconds, coordinates);
            }
        }
    }

    template <typename Func, typename Cases>
    static void testCelestialTrajectory(const Func &func, const Cases &cases, const std::string
        &description, Seconds acceptableDiff, bool showStatistics = false)
    {
        StatisticsAggregator statistics;
        for (const auto &test : cases) {
            auto seconds = utcFromDateString(test.dateString, test.timezone);
            auto secondsShifted = seconds;
            auto result = func(secondsShifted, test.coordinates);

            if (!result.has_value()) {
                std::stringstream s;
                s << description << " test at " << test.dateString << " has no value";
                throw std::logic_error(s.str());
            }
            Seconds timesDiff = (*result - seconds);

            if (showStatistics) {
                statistics.account(description + " timesDiff ", timesDiff);
                auto moonLat =
                    degreesFromRads(Sky::observe<Moon>(seconds, test.coordinates).horizontal.lat);
                statistics.account(description + " moon lat", moonLat);
            } else {
                SKYGAZING_ASSERT_SMALL(timesDiff, acceptableDiff);
            }
        }
    }

    template <typename CelestialBody>
    static void testCelestialTrajectory(
        const std::vector<CelestialTrajectoryTestData::RiseTransitSetCase> &cases,
        const std::string &description, Seconds acceptableDiff, bool showStatistics = false)
    {
        StatisticsAggregator statistics;
        std::mt19937 gen(std::random_device{}());
        constexpr auto kRandomizeWindow = 0.4 * CelestialBody::approxDayLength() * kSecondsInDay;
        std::uniform_real_distribution<Seconds> shiftDistr(-kRandomizeWindow, kRandomizeWindow);
        for (const auto &test : cases) {
            auto testRiseSeconds = utcFromDateString(test.riseDateString, test.timezone);
            auto testTransitSeconds = utcFromDateString(test.transitDateString, test.timezone);
            auto testSetSeconds = utcFromDateString(test.setDateString, test.timezone);
            for (int i = 0; i < 10; ++i) {
                auto testSeconds = testTransitSeconds + shiftDistr(gen);
                auto trajectoryCycle = Sky::getTrajectoryCycleFromUTC<CelestialBody>(
                        testSeconds, test.coordinates);

                Seconds riseDiff = (trajectoryCycle.getTodaysRiseUTC().value() - testRiseSeconds);
                Seconds transitDiff = (trajectoryCycle.getTransitUTC() - testTransitSeconds);
                Seconds setDiff = (trajectoryCycle.getTodaysSetUTC().value() - testSetSeconds);

                if (showStatistics) {
                    statistics.account(description + " riseDiff ", riseDiff);
                    statistics.account(description + " transitDiff ", transitDiff);
                    statistics.account(description + " setDiff ", setDiff);
                } else {
                    SKYGAZING_ASSERT_SMALL(riseDiff, acceptableDiff);
                    SKYGAZING_ASSERT_SMALL(transitDiff, acceptableDiff);
                    SKYGAZING_ASSERT_SMALL(setDiff, acceptableDiff);
                }
            }
        }
    }

    static void printClosestHorizonCrossings(
            const CelestialTrajectoryTestData::SingleTimeCase &testCase)
    {
        auto utc = utcFromDateString(testCase.dateString, testCase.timezone);
        auto firstSunriseAfterCurrentTransit =
            Sky::getClosestToCurrentTransitHorizonCrossingUTC<Sun>(
                utc, testCase.coordinates, true, true);
        std::cout << "firstSunriseAfterCurrentTransit: "
            << (firstSunriseAfterCurrentTransit
                ? dateStringFromUTC(*firstSunriseAfterCurrentTransit, testCase.timezone) : "null"s)
            << std::endl;
        auto firstSunsetBeforeCurrentTransit =
            Sky::getClosestToCurrentTransitHorizonCrossingUTC<Sun>(
                utc, testCase.coordinates, false, false);
        std::cout << "firstSunsetBeforeCurrentTransit: "
            << (firstSunsetBeforeCurrentTransit
                ? dateStringFromUTC(*firstSunsetBeforeCurrentTransit, testCase.timezone) : "null"s)
            << std::endl;
    }

    static void printOnePoint(const CelestialTrajectoryTestData::SingleTimeCase &testCase) {
        std::cout << "Printing one point, dateString: " << testCase.dateString
            << " tz: " << testCase.timezone << " observer: "
            << testCase.coordinates.lat << ", " << testCase.coordinates.lng << std::endl;
        UTC utc = utcFromDateString(testCase.dateString, testCase.timezone);
        auto cycle = Sky::getTrajectoryCycleFromUTC<Sun>(utc, testCase.coordinates);
        std::cout << "last sun nadir: "
                  << dateStringFromUTC(cycle.getBeginNadir(), testCase.timezone) << std::endl;
        std::cout << "sun transit: "
                  << dateStringFromUTC(cycle.getTransitUTC(), testCase.timezone) << std::endl;
        std::cout << "next sun nadir: "
                  << dateStringFromUTC(cycle.getEndNadir(), testCase.timezone) << std::endl;
        auto horizontalCoordinates =
            Sky::getHorizontalDegreesCoordinates<Sun>(utc, testCase.coordinates);
        std::cout << "altitude angle: " << horizontalCoordinates.lat
            << ", azimuth: " << horizontalCoordinates.lng << std::endl;
    }

    static void printAltitudeAngleSeries(
        const CelestialTrajectoryTestData::SingleTimeCase &testCase, Seconds freq, int count)
    {
        Seconds seconds = utcFromDateString(testCase.dateString, testCase.timezone);
        for (int i = -count / 2; i < count / 2; ++i) {
            auto s = seconds + freq * i;
            auto hz = Sky::getHorizontalDegreesCoordinates<Sun>(s, testCase.coordinates);
            std::cout << hz.lat << "\t" << dateStringFromUTC(s) << std::endl;
        }
    }

    template <typename CelestialBody>
    static void printTrasitShifts(const CelestialTrajectoryTestData::SingleTimeCase &testCase) {
        auto seconds = utcFromDateString(testCase.dateString, testCase.timezone);
        auto startTransit = *Sky::getTransitUTC<CelestialBody>(seconds, testCase.coordinates);
        auto secondsDay = kSecondsInDay * CelestialBody::approxDayLength();

        for (int i = 0; i < 366; ++i) {
            auto nextTransit = *Sky::getTransitUTC<CelestialBody>(seconds, testCase.coordinates);
            std::cout << (nextTransit - i * secondsDay - startTransit) << std::endl;
            seconds = nextTransit + secondsDay;
        }
    }

    template <typename CelestialBody>
    static void analyzeCelestialTrajectoryStatistics(int seed) {
        RandomDataGenerator dataGenerator(seed);
        StatisticsAggregator statistics;
        const int kIterations = 1000;
        for (int year = 2017; year < 2030; ++year) {
            for (int i = 0; i < kIterations; ++i)
            {
                auto date = dataGenerator.getRandomDate(year);
                auto seconds = timegm(&date);
                auto observer = dataGenerator.getRandomCoordinates(0, 89);
                try {
                    auto trajectoryCycle = Sky::getTrajectoryCycleFromUTC<CelestialBody>(
                            seconds, observer);
                    auto firstHalfHours = (trajectoryCycle.getTransitUTC()
                        - trajectoryCycle.getBeginNadir()) / 3600;
                    auto dayLengthHours = trajectoryCycle.getLengthDays() * kHoursInDay;

                    auto prevCycle = trajectoryCycle.getShiftedCycle(-1);
                    auto nadirDayHours = (trajectoryCycle.getTransitUTC()
                        - prevCycle.getTransitUTC()) / 3600.;

                    statistics.account(SKYGAZING_NAME_COMMA_VAR(dayLengthHours));
                    statistics.account(SKYGAZING_NAME_COMMA_VAR(nadirDayHours));
                    statistics.account(SKYGAZING_NAME_COMMA_VAR(firstHalfHours));
                } catch (std::logic_error &e) {
                    std::cout << e.what() << std::endl;
                    std::cout << dateStringFromUTC(seconds) << " " << seconds << " "
                        << observer.lat << " " << observer.lng << std::endl;
                }
            }
        }
    }

    static void runAnalytics() {
        printClosestHorizonCrossings({{88, -163.519}, "2017-10-25T23:22:11Z"s, 0});
        CelestialTrajectoryTestData::SingleTimeCase testCase{
            {0, -163.519}, "2017-10-25T23:22:11Z"s, 0};
        std::cout << "sun shifts:\n";
        printTrasitShifts<Sun>(testCase);
        std::cout << "moon shifts:\n";
        printTrasitShifts<Moon>(testCase);
        std::cout << "one point:\n";
        printOnePoint(testCase);
        std::cout << "altitude angle series:\n";
        printAltitudeAngleSeries(testCase, 600, 24 * 6);
        std::random_device trueRandom{};
        std::cout << "sun celestial trajectory statistics:\n";
        analyzeCelestialTrajectoryStatistics<Sun>(trueRandom());
        std::cout << "moon celestial trajectory statistics:\n";
        analyzeCelestialTrajectoryStatistics<Moon>(trueRandom());
    }

    static void runAllTests(bool runAnalytics = false) {
        if (runAnalytics) { Testing::runAnalytics(); }
        SKYGAZING_TIME_EXECUTION(testTimeTransforms());
        SKYGAZING_TIME_EXECUTION(
            testCelestialTrajectory<Sun>(CelestialTrajectoryTestData::sunCases, "sun", 90));
        SKYGAZING_TIME_EXECUTION(testTransitExistence<Sun>(100));
        SKYGAZING_TIME_EXECUTION(testTransitExistence<Moon>(10));
        SKYGAZING_TIME_EXECUTION(testCelestialTrajectory(Sky::getTodaysRiseUTC<Moon>,
            CelestialTrajectoryTestData::moonrises, "moonrises", 90));
        SKYGAZING_TIME_EXECUTION(testCelestialTrajectory(Sky::getTodaysSetUTC<Moon>,
            CelestialTrajectoryTestData::moonsets, "moonsets", 90));
    }
};

}

#endif // SKYGAZING_TESTS_H
