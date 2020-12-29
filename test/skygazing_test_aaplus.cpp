#include <random>
#include <chrono>
#include <set>
#include <ctime>

#include "AA+.h"
#include "../skygazing.h"
#include "skygazing_test.h"
#include "skygazing_test_analytics.h"

using namespace Skygazing;
using std::literals::operator""s;

struct AAPlus {
    static Sky::CelestialObjectObservation observeSun(double julian, DegreesCoordinates coordinates, bool highPrecision) {
        Sky::CelestialObjectObservation sun{};
        const double terrestrialJulian = CAADynamicalTime::UTC2TT(julian);
        sun.julian = terrestrialJulian;
        sun.ecliptic = Coordinates{radsFromDegrees(CAASun::ApparentEclipticLatitude(terrestrialJulian, highPrecision)),
                                   radsFromDegrees(CAASun::ApparentEclipticLongitude(terrestrialJulian, highPrecision))};
        double SunLong = CAASun::ApparentEclipticLongitude(terrestrialJulian, highPrecision);
        double SunLat = CAASun::ApparentEclipticLatitude(terrestrialJulian, highPrecision);
        CAA2DCoordinate aaEquatorial =
                CAACoordinateTransformation::Ecliptic2Equatorial(SunLong, SunLat,
                                                                 CAANutation::TrueObliquityOfEcliptic(terrestrialJulian));
        sun.equatorial = {CAACoordinateTransformation::DegreesToRadians(aaEquatorial.Y),
                          CAACoordinateTransformation::HoursToRadians(aaEquatorial.X)};
        sun.distance = CAAEarth::RadiusVector(terrestrialJulian, highPrecision);
        double Height = 0;
        const CAA2DCoordinate SunTopo = CAAParallax::Equatorial2Topocentric(aaEquatorial.X, aaEquatorial.Y, sun.distance,
                                                                            -coordinates.lng, coordinates.lat, Height, terrestrialJulian);
        double AST = CAASidereal::ApparentGreenwichSiderealTime(julian);
        double LongtitudeAsHourAngle = CAACoordinateTransformation::DegreesToHours(-coordinates.lng);
        double LocalHourAngle = AST - LongtitudeAsHourAngle - SunTopo.X;
        sun.hourAngle = CAACoordinateTransformation::HoursToRadians(LocalHourAngle);
        CAA2DCoordinate SunHorizontal = CAACoordinateTransformation::Equatorial2Horizontal(LocalHourAngle, SunTopo.Y, coordinates.lat);
        SunHorizontal.Y += CAARefraction::RefractionFromTrue(SunHorizontal.Y, 1010, 10);
        Coordinates horizontal{radsFromDegrees(SunHorizontal.Y), radsFromDegrees(SunHorizontal.X)};
        sun.horizontal = horizontal.normalize();
        return sun;
    }

    static std::pair<Sky::CelestialObjectObservation, Moon::MoonDetails> observeMoon(double julian, double latitude, double longitude) {
        Sky::CelestialObjectObservation moon{};
        Moon::MoonDetails moonDetails{};
        moon.observer = Coordinates{radsFromDegrees(latitude), radsFromDegrees(longitude)};
        longitude = -longitude;
        const double jd = CAADynamicalTime::UTC2TT(julian);
        moon.julian = jd - kJulianYear2000;
        double eclipticLng = CAAMoon::EclipticLongitude(jd);
        double eclipticLat = CAAMoon::EclipticLatitude(jd);
        moon.ecliptic = Coordinates{radsFromDegrees(eclipticLat), radsFromDegrees(eclipticLng)};
        auto Equatorial = CAACoordinateTransformation::Ecliptic2Equatorial(eclipticLng, eclipticLat, CAANutation::TrueObliquityOfEcliptic(jd));
        moon.equatorial = Coordinates{CAACoordinateTransformation::DegreesToRadians(Equatorial.Y),
                                      CAACoordinateTransformation::HoursToRadians(Equatorial.X)};
        double MoonRad = CAAMoon::RadiusVector(jd);
        MoonRad /= 149597870.691; //Convert KM to AU
        auto Height = 0;
        const CAA2DCoordinate topocentric = CAAParallax::Equatorial2Topocentric(Equatorial.X, Equatorial.Y, MoonRad, longitude, latitude, Height, jd);
        auto AST = CAASidereal::ApparentGreenwichSiderealTime(julian);
        auto longtitudeAsHourAngle = CAACoordinateTransformation::DegreesToHours(longitude);
        auto localHourAngle = AST - longtitudeAsHourAngle - topocentric.X;
        moon.hourAngle = CAACoordinateTransformation::HoursToRadians(localHourAngle);
        CAA2DCoordinate MoonHorizontal = CAACoordinateTransformation::Equatorial2Horizontal(localHourAngle, topocentric.Y, latitude);
        MoonHorizontal.Y += CAARefraction::RefractionFromTrue(MoonHorizontal.Y, 1013, 10);
        moon.horizontal = Coordinates{radsFromDegrees(MoonHorizontal.Y), radsFromDegrees(MoonHorizontal.X)};

        auto sun = observeSun(julian, {latitude, longitude}, true);
        double moon_alpha = degreesFromRads(moon.rightAscension()) / 15.;
        double moon_delta = degreesFromRads(moon.declination());
        double sun_alpha = degreesFromRads(sun.rightAscension()) / 15.;
        double sun_delta = degreesFromRads(sun.declination());

        moonDetails.geocentricElongation = Sky::getGeocentricElongation(sun.equatorial, moon.equatorial);

        double position_angle = CAAMoonIlluminatedFraction::PositionAngle(sun_alpha, sun_delta, moon_alpha, moon_delta);
        double phase_angle = CAAMoonIlluminatedFraction::PhaseAngle(moonDetails.geocentricElongation, 368410.0, 149971520.0);
        double illuminated_fraction = CAAMoonIlluminatedFraction::IlluminatedFraction(phase_angle);
        moonDetails.angle = radsFromDegrees(position_angle);
        moonDetails.phase = radsFromDegrees(phase_angle);
        moonDetails.fraction = radsFromDegrees(illuminated_fraction);
        moon.parallacticAngle = radsFromDegrees(
                CAAParallactic::ParallacticAngle(CAACoordinateTransformation::RadiansToHours(moon.hourAngle),
                                                 degreesFromRads(moon.observer.lat), degreesFromRads(moon.declination())));
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
            auto seconds = tmToSeconds(&date);
            auto observer = dataGenerator.getRandomCoordinates();

            auto sgSun = Sky::observe<Sun>(seconds, observer);
            auto aaSun = AAPlus::observeSun(secondsToJulianSinceJulianEpoch(seconds), observer, true);
            auto sunEclipticDistance = degreesFromRads(haversineDistance(sgSun.ecliptic, aaSun.ecliptic));
            auto sunEquatorialDistance = degreesFromRads(haversineDistance(sgSun.equatorial, aaSun.equatorial));
            auto sunHorizontalDistance = degreesFromRads(haversineDistance(sgSun.horizontal, aaSun.horizontal));
            SKYGAZING_ASSERT_SMALL(sunEclipticDistance, 0.5);
            SKYGAZING_ASSERT_SMALL(sunEquatorialDistance, 0.5);
            SKYGAZING_ASSERT_SMALL(sunHorizontalDistance, 0.5);

            if (showStatistics) {
                statistics.account("sunEclipticDistance", sunEclipticDistance);
                statistics.account("sunEquatorialDistance", sunEquatorialDistance);
                statistics.account("sunHorizontalDistance", sunHorizontalDistance);
                statistics.accountForRadsDiff("sunHourAngleDiff", sgSun.hourAngle, aaSun.hourAngle);
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
            auto [aaMoon, aaMoonDetails] = AAPlus::observeMoon(secondsToJulianSinceJulianEpoch(seconds), observer.lat, observer.lng);

            auto moonEclipticDistance = degreesFromRads(haversineDistance(sgMoon.ecliptic, aaMoon.ecliptic));
            auto moonEquatorialDistance = degreesFromRads(haversineDistance(sgMoon.equatorial, aaMoon.equatorial));
            auto moonHorizontalDistance = degreesFromRads(haversineDistance(sgMoon.horizontal, aaMoon.horizontal));
            auto geocentricElongationDiff = degreesFromRads(sgMoonDetails.geocentricElongation - aaMoonDetails.geocentricElongation);
            SKYGAZING_ASSERT_SMALL(moonEclipticDistance, 0.01);
            SKYGAZING_ASSERT_SMALL(moonEquatorialDistance, 0.01);
            SKYGAZING_ASSERT_SMALL(moonHorizontalDistance, 1.5);

            if (showStatistics) {
                statistics.account("moonEclipticDistance", moonEclipticDistance);
                statistics.account("moonEquatorialDistance", moonEquatorialDistance);
                statistics.account("moonHorizontalDistance", moonHorizontalDistance);
                statistics.accountForRadsDiff("moonHourAngleDiff", sgMoon.hourAngle, aaMoon.hourAngle);
                statistics.accountForRadsDiff("moonParallacticAngleDiff",
                                              sgMoon.getParallacticAngle(), aaMoon.parallacticAngle.value());
                statistics.account("geocentricElongationDiff", geocentricElongationDiff);
            }
        }
    }
}


int main() {
    const auto startTime = std::chrono::high_resolution_clock::now();
    SKYGAZING_TIME_EXECUTION(Testing::runAllTests());

    std::random_device trueRandom;
    SKYGAZING_TIME_EXECUTION(testMoonObservation(trueRandom()));
    SKYGAZING_TIME_EXECUTION(testSunObservation(trueRandom()));

    const auto finishTime = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double, std::milli> ms = finishTime - startTime;
    std::cout << std::fixed << "total test time " << ms.count() / 1000. << "s\n";
}
