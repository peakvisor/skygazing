// [AA] "Astronomical Algorithms" 2nd edition by Jean Meeus (Willmann-Bell, Richmond) 1998.
#ifndef SKYGAZING_OBSERVATION_H
#define SKYGAZING_OBSERVATION_H

#include "skygazing_math.h"
#include "skygazing_time.h"

namespace Skygazing {

struct Sky {
    static constexpr Rads kEarthObliquity = 23.4397_deg;
    static constexpr Julian kTimesFindingPrecision = 10. / kSecondsInDay;
    static constexpr Julian kTimesFindingStep = 300. / kSecondsInDay;
    static constexpr Julian kTransitHalfWindow = 13 / kHoursInDay;

    struct EclipticPosition {
        Coordinates coordinates;
        double distance;
    };

    struct CelestialObjectObservation {
        CelestialObjectObservation() = default;
        CelestialObjectObservation(CelestialObjectObservation &&other) = default;
        CelestialObjectObservation(Julian julian, const EclipticPosition &eclipticPosition,
                                   const DegreesCoordinates &observerCoordinates, bool addRefraction = true)
            : julian{julian}
            , observer{radsFromDegrees(observerCoordinates)}
            , ecliptic{eclipticPosition.coordinates}
            , equatorial{eclipticToEquatorial(ecliptic)}
            , hourAngle{getLocalHourAngle(julian, observer.lng, equatorial.lng)}
            , horizontal{getHorizontalCoordinatesFromHourAngle(hourAngle, observer.lat, equatorial.lat)}
            , distance(eclipticPosition.distance)
        {
            if (addRefraction) {
                accountForRefraction();
            }
        }

        Rads declination() const { return equatorial.lat; }
        Rads rightAscension() const { return equatorial.lng; }
        Rads altitudeAngle() const {
            return normalizeHugeRads(horizontal.lat);
        }
        Rads azimuth() const { return horizontal.lng; }
        Rads getParallacticAngle() {
            if (!parallacticAngle) {
                parallacticAngle = std::atan2(std::sin(hourAngle),
                    std::tan(observer.lat) * std::cos(declination()) - std::sin(declination()) * std::cos(hourAngle));
            }
            return *parallacticAngle;
        }

        void accountForRefraction() {
            horizontal.lat += getRefractionFromTrue(horizontal.lat);
        }

        Julian julian;
        Coordinates observer;
        Coordinates ecliptic;
        Coordinates equatorial;
        double hourAngle;
        Coordinates horizontal;
        double distance;
        std::optional<double> parallacticAngle;
    };

    template <typename CelestialObject>
    static CelestialObjectObservation observeInTerrestrialTime(Julian terrestrialJulian, const DegreesCoordinates &observer) {
        auto eclipticPosition = CelestialObject::getEclipticPosition(terrestrialJulian);
        return CelestialObjectObservation(terrestrialJulian, eclipticPosition, observer);
    }

    template <typename CelestialObject>
    static CelestialObjectObservation observe(Seconds seconds, const DegreesCoordinates &observer) {
        Julian terrestrialJulian = terrestrialJulianFromSeconds(seconds);
        return observeInTerrestrialTime<CelestialObject>(terrestrialJulian, observer);
    }

    template <typename CelestialObject>
    static DegreesCoordinates getHorizontalDegreesCoordinates(Seconds seconds, const DegreesCoordinates &observer) {
        auto observation = observe<CelestialObject>(seconds, observer);
        return degreesFromRads(observation.horizontal);
    }

    static Rads getRefractionFromTrue(Rads altitudeAngle) noexcept {
        altitudeAngle = normalizeHugeRads(altitudeAngle);
        constexpr double lowerBound = -1.9006387000003735_deg;

        if (altitudeAngle <= lowerBound) {
            altitudeAngle = lowerBound;
        }

        // formula 16.4 of [AA]
        // 1.02 / tan(h + 10.3 / (h + 5.11)) + 0.0019279 h in degrees, result in arc minutes -> converted to rad:
        constexpr Rads coeff1 = 10.3 * kRadsPerDegree * kRadsPerDegree;
        constexpr Rads coeff2 = 5.11_deg;
        constexpr Rads numeratorRads = 1.02_deg / 60;
        constexpr Rads correctionRads = 0.0019279_deg / 60;
        return numeratorRads / (tan(altitudeAngle + coeff1 / (altitudeAngle + coeff2))) + correctionRads;
    }

    static Rads correctRefractionForPressureAndTemperature(double refraction, double pressureMillibars, double temperature) {
        return refraction * (pressureMillibars / 1010 * 283 / (273 + temperature));
    }

    static Rads getGreenwichSiderealTime(Julian julian) {
//    return 280.16_deg + 360.98564736629_deg * julian; // [AA] fomula (12.4)
        return 280.16_deg + 360.9856235_deg * julian;
    }

    static Rads getLocalHourAngle(Julian julian, Rads lng, Rads rightAscension) {
        Rads result = getGreenwichSiderealTime(julian) + lng - rightAscension;
        return normalizeHugeRads(result);
    }

    static Coordinates eclipticToEquatorial(const Coordinates &ecliptic, Rads obliquity = kEarthObliquity) {
        auto declination = std::asin(sin(ecliptic.lat) * std::cos(obliquity)
                                     + std::cos(ecliptic.lat) * std::sin(obliquity) * std::sin(ecliptic.lng));
        auto rightAscension = std::atan2(std::sin(ecliptic.lng) * std::cos(obliquity) - std::tan(ecliptic.lat) * std::sin(obliquity),
                                         std::cos(ecliptic.lng));
        return {declination, rightAscension};
    }

    static Coordinates getHorizontalCoordinatesFromHourAngle(Rads hourAngle, Rads lat, Rads declination) {
        Rads altitudeAngle = std::asin(std::sin(lat) * std::sin(declination) + std::cos(lat) * std::cos(declination) * std::cos(hourAngle));
        Rads azimuth = std::atan2(std::sin(hourAngle), std::cos(hourAngle) * std::sin(lat) - std::tan(declination) * std::cos(lat));

        return {normalizeRads(altitudeAngle), normalizeRads(azimuth)};
    }

    static Coordinates getHorizontalCoordinates(Julian days, const Coordinates &topographic, const Coordinates &equatorial) {
        Rads hourAngle = getLocalHourAngle(days, topographic.lng, equatorial.lng);
        return getHorizontalCoordinatesFromHourAngle(hourAngle, topographic.lat, equatorial.lat);
    }

    static double getGeocentricElongation(const Coordinates &sunEquatorial, const Coordinates &objectEquatorial) {
        double cosGeocentricElongation =
                std::sin(sunEquatorial.lat) * std::sin(objectEquatorial.lat)
                + std::cos(sunEquatorial.lat) * std::cos(objectEquatorial.lat) * std::cos(sunEquatorial.lng - objectEquatorial.lng);
        return std::acos(clamp(cosGeocentricElongation, 1.));
    }

    struct NotableTimesInSeconds {
        std::optional<Seconds> rise;
        std::optional<Seconds> transit;
        std::optional<Seconds> set;
    };

    struct NotableTimesInTerrestrialTime {
        std::optional<Julian> rise;
        std::optional<Julian> transit;
        std::optional<Julian> set;

        NotableTimesInSeconds toSeconds() const {
            return {secondsFromTerrestrialJulian(rise), secondsFromTerrestrialJulian(transit), secondsFromTerrestrialJulian(set)};
        }
    };

    template <typename CelestialBody>
    static std::optional<Julian> getTransitInTerrestrialFromTerrestrial(Julian terrestrialJulian,
                                                                        const DegreesCoordinates &observer)
    {
        auto getAltitudeAngle = [&observer](Julian julian) {
            auto observation = Sky::observeInTerrestrialTime<CelestialBody>(julian, observer);
            return observation.altitudeAngle();
        };
        auto altitudeAngleAnalyzer = FunctionAnalyzer{std::move(getAltitudeAngle)};
        auto maxSample = altitudeAngleAnalyzer.findMax(terrestrialJulian - kTransitHalfWindow,
                                                       terrestrialJulian + kTransitHalfWindow,
                                                       kTimesFindingStep, kTimesFindingPrecision);
        if (maxSample) {
            return maxSample->arg;
        } else {
            return std::nullopt;
        }
    }

    template <typename CelestialBody>
    static std::optional<Seconds> getTransit(Seconds seconds, const DegreesCoordinates &observer) {
        auto transit = getTransitInTerrestrialFromTerrestrial<CelestialBody>(terrestrialJulianFromSeconds(seconds), observer);
        if (transit) {
            transit = secondsFromTerrestrialJulian(*transit);
        }
        return transit;
    }

    template <typename CelestialBody>
    static NotableTimesInTerrestrialTime getNotableTimesInTerrestrialTime(Seconds seconds, const DegreesCoordinates &observer) {
        Julian someJulianAtTheDate = terrestrialJulianFromSeconds(seconds);

        auto transit = getTransitInTerrestrialFromTerrestrial<CelestialBody>(someJulianAtTheDate, observer);
        if (transit.has_value()) {
            auto upperEdgeAltitudeAngleAnalyzer = FunctionAnalyzer{[&observer](Julian julian) {
                auto observation = Sky::observeInTerrestrialTime<CelestialBody>(julian, observer);
                return observation.altitudeAngle() + CelestialBody::meanAngularSize() / 2;
            }};
            auto rise = upperEdgeAltitudeAngleAnalyzer.findSignChangeToPositive(*transit - kTransitHalfWindow, *transit,
                                                                                kTimesFindingStep, kTimesFindingPrecision);
            auto set = upperEdgeAltitudeAngleAnalyzer.findSignChangeToPositive(*transit + kTransitHalfWindow, *transit,
                                                                               kTimesFindingStep, kTimesFindingPrecision);
            return {rise, transit, set};
        }
        return {};
    }

    template <typename CelestialBody>
    static NotableTimesInSeconds getTimes(Seconds seconds, const DegreesCoordinates &observer) {
        return getNotableTimesInTerrestrialTime<CelestialBody>(seconds, observer).toSeconds();
    }

    template <typename CelestialBody>
    static std::optional<Seconds> getSet(Seconds seconds, const DegreesCoordinates &observer) {
        return getTimes<CelestialBody>(seconds, observer).set;
    }

    template <typename CelestialBody>
    static std::optional<Seconds> getRise(Seconds seconds, const DegreesCoordinates &observer) {
        return getTimes<CelestialBody>(seconds, observer).rise;
    }
};

}

#endif //SKYGAZING_OBSERVATION_H
