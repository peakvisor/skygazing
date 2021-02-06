#include <random>
#include <chrono>
#include <set>
#include <ctime>

#include "AA+.h"
#include "../include/skygazing.h"
#include "skygazing_test.h"
#include "skygazing_test_analytics.h"

using namespace Skygazing;
using std::literals::operator""s;

struct AAPlus {
    static constexpr double kAU = 149597870691;

    static Sky::CelestialObjectObservation observeSun(double julian, DegreesCoordinates coordinates,
            bool highPrecision)
    {
        Sky::CelestialObjectObservation sun{};
        const double terrestrialJulian = CAADynamicalTime::UTC2TT(julian);
        sun.tt = terrestrialJulian - kJulianYear2000;
        double SunLong = CAASun::ApparentEclipticLongitude(terrestrialJulian, highPrecision);
        double SunLat = CAASun::ApparentEclipticLatitude(terrestrialJulian, highPrecision);
        sun.ecliptic = {radsFromDegrees(SunLat), radsFromDegrees(SunLong)};
        CAA2DCoordinate aaEquatorial =
            CAACoordinateTransformation::Ecliptic2Equatorial(SunLong, SunLat,// 23.4397);
                CAANutation::TrueObliquityOfEcliptic(terrestrialJulian));
        sun.equatorial = {CAACoordinateTransformation::DegreesToRadians(aaEquatorial.Y),
            CAACoordinateTransformation::HoursToRadians(aaEquatorial.X)};
        auto rad = CAAEarth::RadiusVector(terrestrialJulian, highPrecision);
        sun.distance = rad * kAU;
        double Height = 0;
        const CAA2DCoordinate SunTopo = CAAParallax::Equatorial2Topocentric(aaEquatorial.X,
            aaEquatorial.Y, sun.distance / kAU, -coordinates.lng, coordinates.lat,
            Height, terrestrialJulian);
        double AST = CAASidereal::ApparentGreenwichSiderealTime(julian);
        double LongtitudeAsHourAngle =
            CAACoordinateTransformation::DegreesToHours(-coordinates.lng);
        double LocalHourAngle = AST - LongtitudeAsHourAngle - SunTopo.X;
        sun.topocentricHourAngle = CAACoordinateTransformation::HoursToRadians(LocalHourAngle);
        sun.topocentric = {radsFromDegrees(SunTopo.Y),
            CAACoordinateTransformation::HoursToRadians(SunTopo.X)};
        CAA2DCoordinate SunHorizontal = CAACoordinateTransformation::Equatorial2Horizontal(
            LocalHourAngle, SunTopo.Y, coordinates.lat);
        SunHorizontal.Y += CAARefraction::RefractionFromTrue(SunHorizontal.Y, 1010, 10);
        Coordinates horizontal{radsFromDegrees(SunHorizontal.Y), radsFromDegrees(SunHorizontal.X)};
        sun.horizontal = horizontal.normalize();
        return sun;
    }

    static std::pair<Sky::CelestialObjectObservation, Moon::MoonDetails> observeMoon(double julian,
            double latitude, double longitude)
    {
        Sky::CelestialObjectObservation moon{};
        Moon::MoonDetails moonDetails{};
        moon.observer = Coordinates{radsFromDegrees(latitude), radsFromDegrees(longitude)};
        longitude = -longitude;
        const double jd = CAADynamicalTime::UTC2TT(julian);
        moon.tt = jd - kJulianYear2000;
        double eclipticLng = CAAMoon::EclipticLongitude(jd);
        double eclipticLat = CAAMoon::EclipticLatitude(jd);
        moon.ecliptic = Coordinates{radsFromDegrees(eclipticLat), radsFromDegrees(eclipticLng)};
        auto Equatorial = CAACoordinateTransformation::Ecliptic2Equatorial(
            eclipticLng, eclipticLat, CAANutation::TrueObliquityOfEcliptic(jd));
        moon.equatorial = Coordinates{CAACoordinateTransformation::DegreesToRadians(Equatorial.Y),
            CAACoordinateTransformation::HoursToRadians(Equatorial.X)};
        auto moonRad = CAAMoon::RadiusVector(jd);
        moon.distance = moonRad * 1000;
//        MoonRad /= 149597870.691; //Convert KM to AU
        auto Height = 0;
        const CAA2DCoordinate topocentric = CAAParallax::Equatorial2Topocentric(
            Equatorial.X, Equatorial.Y, moonRad / 149597870.691, longitude, latitude, Height, jd);
        auto AST = CAASidereal::ApparentGreenwichSiderealTime(julian);
        auto longtitudeAsHourAngle = CAACoordinateTransformation::DegreesToHours(longitude);
        auto localHourAngle = AST - longtitudeAsHourAngle - topocentric.X;
        moon.topocentric = {radsFromDegrees(topocentric.Y),
            CAACoordinateTransformation::HoursToRadians(topocentric.X)};
        moon.topocentricHourAngle = CAACoordinateTransformation::HoursToRadians(localHourAngle);
        CAA2DCoordinate MoonHorizontal = CAACoordinateTransformation::Equatorial2Horizontal(
            localHourAngle, topocentric.Y, latitude);
        MoonHorizontal.Y += CAARefraction::RefractionFromTrue(MoonHorizontal.Y, 1013, 10);
        moon.horizontal =
            Coordinates{radsFromDegrees(MoonHorizontal.Y), radsFromDegrees(MoonHorizontal.X)};

        auto sun = observeSun(julian, {latitude, longitude}, true);
        double moon_alpha = degreesFromRads(moon.rightAscension()) / 15.;
        double moon_delta = degreesFromRads(moon.declination());
        double sun_alpha = degreesFromRads(sun.rightAscension()) / 15.;
        double sun_delta = degreesFromRads(sun.declination());

        moonDetails.geocentricElongation =
            Sky::getGeocentricElongation(sun.equatorial, moon.equatorial);

        double position_angle =
            CAAMoonIlluminatedFraction::PositionAngle(sun_alpha, sun_delta, moon_alpha, moon_delta);
        double phase_angle = CAAMoonIlluminatedFraction::PhaseAngle(
            moonDetails.geocentricElongation, 368410.0, 149971520.0);
        double illuminated_fraction = CAAMoonIlluminatedFraction::IlluminatedFraction(phase_angle);
        moonDetails.angle = radsFromDegrees(position_angle);
        moonDetails.phase = radsFromDegrees(phase_angle);
        moonDetails.fraction = radsFromDegrees(illuminated_fraction);
        moon.parallacticAngle = radsFromDegrees(CAAParallactic::ParallacticAngle(
            CAACoordinateTransformation::RadiansToHours(moon.hourAngle),
            degreesFromRads(moon.observer.lat), degreesFromRads(moon.declination())
        ));
        return {std::move(moon), moonDetails};
    }
};

void testSunObservation(int seed = 23, bool showStatistics = false) {
    RandomDataGenerator dataGenerator(seed);
    StatisticsAggregator statistics;
    const int kIterations = 1000;
    for (int year = 2017; year < 2030; ++year) {
        for (int i = 0; i < kIterations; ++i) {
            auto date = dataGenerator.getRandomDate(year);
            auto seconds = tmToUTC(&date);
            auto observer = dataGenerator.getRandomCoordinates();

            auto sgSun = Sky::observe<Sun>(seconds, observer, true);
            auto aaSun = AAPlus::observeSun(julianSinceEpochFromUTC(seconds), observer, false);
            auto sunEclipticDistance = degreesFromRads(haversineDistance(
                sgSun.ecliptic, aaSun.ecliptic));
            auto sunEclipticLngDiff = degreesFromRads(normalizeRads(
                sgSun.ecliptic.lng - aaSun.ecliptic.lng));
            auto sunEquatorialDistance = degreesFromRads(haversineDistance(
                sgSun.equatorial, aaSun.equatorial));
            auto sunTopicentricDistance = degreesFromRads(haversineDistance(
                sgSun.topocentric, aaSun.topocentric));
            auto sunHorizontalDistance = degreesFromRads(haversineDistance(
                sgSun.horizontal, aaSun.horizontal));
            auto sunDistanceDiff = sgSun.distance - aaSun.distance;

            if (showStatistics) {
                statistics.account("sunEclipticDistance", sunEclipticDistance);
                statistics.account("sunEclipticLngDiff", sunEclipticLngDiff);
                statistics.account("sunEquatorialDistance", sunEquatorialDistance);
                statistics.account(SKYGAZING_NAME_COMMA_VAR(sunTopicentricDistance));
                statistics.account("sunHorizontalDistance", sunHorizontalDistance);
                statistics.accountForRadsDiff("sunHourAngleDiff",
                    sgSun.topocentricHourAngle, aaSun.topocentricHourAngle);
                statistics.account(SKYGAZING_NAME_COMMA_VAR(sunDistanceDiff));
                statistics.account(SKYGAZING_NAME_COMMA_VAR(sgSun.distance));
            } else {
                SKYGAZING_ASSERT_SMALL(sunEclipticDistance, 0.1);
                SKYGAZING_ASSERT_SMALL(sunEquatorialDistance, 0.1);
                SKYGAZING_ASSERT_SMALL(sunHorizontalDistance, 0.1);
            }
        }
    }
}

void testMoonObservation(int seed, bool showStatistics = false) {
    RandomDataGenerator dataGenerator(seed);
    StatisticsAggregator statistics;
    const int kIterations = 1000;
    for (int year = 2017; year < 2030; ++year) {
        for (int i = 0; i < kIterations; ++i) {
            auto date = dataGenerator.getRandomDate(year);
            auto seconds = timegm(&date);
            auto observer = dataGenerator.getRandomCoordinates();

            auto sgMoon = Sky::observe<Moon>(seconds, observer);
            auto sgSun = Sky::observe<Sun>(seconds, observer);
            auto sgMoonDetails = Moon::getMoonDetails(sgMoon, sgSun);
            sgMoon.accountForRefraction();
            auto [aaMoon, aaMoonDetails] = AAPlus::observeMoon(
                julianSinceEpochFromUTC(seconds), observer.lat, observer.lng);

            auto moonEclipticDistance = degreesFromRads(haversineDistance(
                sgMoon.ecliptic, aaMoon.ecliptic));
            auto moonEquatorialDistance = degreesFromRads(haversineDistance(
                sgMoon.equatorial, aaMoon.equatorial));
            auto moonTopicentricDistance = degreesFromRads(haversineDistance(
                sgMoon.topocentric, aaMoon.topocentric));
            auto moonHorizontalDistance = degreesFromRads(haversineDistance(
                sgMoon.horizontal, aaMoon.horizontal));
            auto geocentricElongationDiff = degreesFromRads(normalizeRads(
                sgMoonDetails.geocentricElongation - aaMoonDetails.geocentricElongation));
            auto moonDistanceDiff = sgMoon.distance - aaMoon.distance;
            SKYGAZING_ASSERT_SMALL(moonEclipticDistance, 0.01);
            SKYGAZING_ASSERT_SMALL(moonEquatorialDistance, 0.01);
            SKYGAZING_ASSERT_SMALL(moonHorizontalDistance, 0.1);

            if (showStatistics) {
                statistics.account("moonEclipticDistance", moonEclipticDistance);
                statistics.account("moonEquatorialDistance", moonEquatorialDistance);
                statistics.account(SKYGAZING_NAME_COMMA_VAR(moonTopicentricDistance));
                statistics.account("moonHorizontalDistance", moonHorizontalDistance);
                statistics.accountForRadsDiff(
                    "moonHourAngleDiff", sgMoon.topocentricHourAngle, aaMoon.topocentricHourAngle);
                statistics.accountForRadsDiff("moonParallacticAngleDiff",
                    sgMoon.getParallacticAngle(), aaMoon.parallacticAngle.value());
                statistics.account("geocentricElongationDiff", geocentricElongationDiff);
                statistics.account(SKYGAZING_NAME_COMMA_VAR(moonDistanceDiff));
                statistics.account(SKYGAZING_NAME_COMMA_VAR(sgMoon.distance));
            }
        }
    }
}

int main(int argc, char **argv) {
    const auto startTime = std::chrono::high_resolution_clock::now();
    bool runAnalytics = argc > 1 && argv[1] == "analytics"s;
    SKYGAZING_TIME_EXECUTION(Testing::runAllTests(runAnalytics));

    std::random_device trueRandom;
    SKYGAZING_TIME_EXECUTION(testSunObservation(trueRandom()));
    SKYGAZING_TIME_EXECUTION(testMoonObservation(trueRandom()));

    const auto finishTime = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double, std::milli> ms = finishTime - startTime;
    std::cout << std::fixed << "total test time " << ms.count() / 1000. << "s\n";
}
