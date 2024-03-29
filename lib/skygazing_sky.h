// [AA] "Astronomical Algorithms" 2nd edition by Jean Meeus (Willmann-Bell, Richmond) 1998.
#pragma once

#include "skygazing_math.h"
#include "skygazing_time.h"
#include <cassert>

namespace Skygazing {

struct Earth {
    static constexpr double equatorialRadius = 6378137;
    static constexpr double flattening = 1 / 298.257223563;
    static constexpr Rads obliquity = 23.4397_deg;
};

struct Sky {
    static constexpr double kAstronomicalUnit = 149597870691.;

    static constexpr Days kTimesFindingPrecision = 1. / kSecondsInDay;
    static constexpr Days kTimesFindingStep = 300. / kSecondsInDay;
    static constexpr int kCyclesToSearchForSetRise = 3;
    static constexpr Rads kRefractionLowerBound = -1.9006387000003735_deg;

    struct EclipticPosition {
        Coordinates coordinates;
        double distance;
    };

    struct EquatorialPosition {
        Coordinates coordinates;
        double distance;
    };

    struct CelestialObjectObservation {
        CelestialObjectObservation() = default;
        CelestialObjectObservation(CelestialObjectObservation &&other) = default;
        CelestialObjectObservation(TT tt, const EclipticPosition &eclipticPosition,
                const DegreesCoordinates &observerCoordinates, double observerAltitude = 0)
                : tt{tt}
                , observer{radsFromDegrees(observerCoordinates)}
                , ecliptic{eclipticPosition.coordinates}
                , equatorial{eclipticToEquatorial(ecliptic)}
                , hourAngle{getLocalHourAngle(tt, observer.lng, equatorial.lng)}
                , topocentric{getTopocentricCoordinates(equatorial, eclipticPosition.distance,
                    observer, observerAltitude, hourAngle)} // non-rotating earth ellispoid
                , topocentricHourAngle{getLocalHourAngle(tt, observer.lng, topocentric.lng)}
                , horizontal{getHorizontalCoordinatesFromHourAngle(
                    topocentricHourAngle,
                    observer.lat, topocentric.lat)}
                , geocentricDistance(eclipticPosition.distance)
                , distance(getDistanceFromObserver(horizontal.lat, geocentricDistance)) {}

        CelestialObjectObservation(TT tt, const EquatorialPosition &equatorialPosition,
            const DegreesCoordinates &observerCoordinates, double observerAltitude = 0)
            : tt{tt}
            , observer{radsFromDegrees(observerCoordinates)}
            , equatorial{equatorialPosition.coordinates}
            , hourAngle{getLocalHourAngle(tt, observer.lng, equatorial.lng)}
            , topocentric{getTopocentricCoordinates(equatorial, equatorialPosition.distance,
                observer, observerAltitude, hourAngle)} // non-rotating earth ellispoid
            , topocentricHourAngle{getLocalHourAngle(tt, observer.lng, topocentric.lng)}
            , horizontal{getHorizontalCoordinatesFromHourAngle(
                topocentricHourAngle,
                observer.lat, topocentric.lat)}
            , geocentricDistance(equatorialPosition.distance)
            , distance(getDistanceFromObserver(horizontal.lat, geocentricDistance)) {}

        Rads declination() const { return equatorial.lat; }
        Rads rightAscension() const { return equatorial.lng; }
        Rads altitudeAngle() const {
            return normalizeHugeRads(horizontal.lat);
        }
        Rads azimuth() const { return horizontal.lng; }
        Rads getParallacticAngle() {
            if (!parallacticAngle) {
                const auto &dec = declination();
                parallacticAngle = std::atan2(std::sin(hourAngle),
                    std::tan(observer.lat) * std::cos(dec) - std::sin(dec) * std::cos(hourAngle));
            }
            return *parallacticAngle;
        }

        void accountForRefraction() {
            horizontal.lat += getRefractionFromTrue(horizontal.lat);
        }

        TT tt;
        Coordinates observer;
        Coordinates ecliptic;
        Coordinates equatorial;
        double hourAngle;
        Coordinates topocentric;
        double topocentricHourAngle;
        Coordinates horizontal;
        double geocentricDistance;
        double distance; // from observer
        std::optional<double> parallacticAngle;
    };

    template <typename CelestialObject>
    static CelestialObjectObservation observeInTT(TT tt, const DegreesCoordinates &observer,
            bool addRefraction = true)
    {
        auto eclipticPosition = CelestialObject::getEclipticPosition(tt);
        auto observation = CelestialObjectObservation(tt, eclipticPosition, observer);
        if (addRefraction) { observation.accountForRefraction(); }
        return observation;
    }

    template <typename CelestialObject>
    static CelestialObjectObservation observe(UTC seconds,
            const DegreesCoordinates &observer, bool addRefraction = true)
    {
        TT tt = ttFromUTC(seconds);
        auto observation = observeInTT<CelestialObject>(tt, observer, addRefraction);
        return observation;
    }

    static CelestialObjectObservation observeEquatorialInTT(
            TT tt, const DegreesCoordinates &observer,
            const DegreesCoordinates &equatorialCoordinates,
            double distance,
            bool addRefraction = true)
    {
        EquatorialPosition position{radsFromDegrees(equatorialCoordinates), distance};
        auto observation = CelestialObjectObservation(tt, position, observer);
        if (addRefraction) { observation.accountForRefraction(); }
        return observation;
    }

    static CelestialObjectObservation observeEquatorial(
            UTC seconds, const DegreesCoordinates &observer,
            const DegreesCoordinates &equatorialCoordinates,
            double distance,
            bool addRefraction = true)
   {
        TT tt = ttFromUTC(seconds);
        auto observation = observeEquatorialInTT(tt, observer, equatorialCoordinates,
                                                 distance, addRefraction);
        return observation;
    }

    template <typename CelestialObject>
    static DegreesCoordinates getHorizontalDegreesCoordinates(Seconds seconds,
            const DegreesCoordinates &observer)
    {
        auto observation = observe<CelestialObject>(seconds, observer);
        return degreesFromRads(observation.horizontal);
    }

    static Rads getRefractionFromTrue(Rads altitudeAngle) noexcept {
        altitudeAngle = normalizeHugeRads(altitudeAngle);
        constexpr double lowerBound = kRefractionLowerBound;

        if (altitudeAngle <= lowerBound) {
            altitudeAngle = lowerBound;
        }

        // formula 16.4 of [AA]
        // 1.02 / tan(h + 10.3 / (h + 5.11)) + 0.0019279 h in degrees, result in arc minutes
        // -> converted to rad:
        constexpr Rads coeff1 = 10.3 * kRadsPerDegree * kRadsPerDegree;
        constexpr Rads coeff2 = 5.11_deg;
        constexpr Rads numerator = 1.02_deg / 60;
        constexpr Rads correction = 0.0019279_deg / 60;
        return numerator / (tan(altitudeAngle + coeff1 / (altitudeAngle + coeff2))) + correction;
    }

    static Rads correctRefraction(double refraction, double pressureMillibars, double temperature) {
        return refraction * (pressureMillibars / 1010 * 283 / (273 + temperature));
    }

    static Rads getGreenwichSiderealTime(TT tt) {
        return Horner(tt)(280.16_deg, 360.98564736629_deg); // [AA] formula 12.4
    }

    static Rads getLocalHourAngle(TT tt, Rads lng, Rads rightAscension) {
        Rads result = getGreenwichSiderealTime(tt) + lng - rightAscension;
        return normalizeHugeRads(result);
    }

    static Coordinates eclipticToEquatorial(const Coordinates &e, Rads obliquity = Earth::obliquity)
    {
        auto declination = std::asin(sin(e.lat) * std::cos(obliquity)
            + std::cos(e.lat) * std::sin(obliquity) * std::sin(e.lng));
        auto rightAscension = std::atan2(
            std::sin(e.lng) * std::cos(obliquity) - std::tan(e.lat) * std::sin(obliquity),
            std::cos(e.lng));
        return {declination, rightAscension};
    }

    static Coordinates getHorizontalCoordinatesFromHourAngle(Rads ha, Rads lat, Rads declination) {
        Rads altitudeAngle = std::asin(
            std::sin(lat) * std::sin(declination)
            + std::cos(lat) * std::cos(declination) * std::cos(ha));
        Rads azimuth = std::atan2(
            std::sin(ha),
            std::cos(ha) * std::sin(lat) - std::tan(declination) * std::cos(lat));
        return {normalizeRads(altitudeAngle), normalizeRads(azimuth)};
    }

    static Coordinates getTopocentricCoordinates(const Coordinates &equatorial, double distance,
            const Coordinates &observer, double observerAltitude, Rads hourAngle)
    { // [AA] chapter 40 + chapter 11
        auto u = std::atan((1 - Earth::flattening) * std::tan(observer.lat));
        auto polarAxisComponent = (1 - Earth::flattening) * std::sin(u)
            + (observerAltitude / Earth::equatorialRadius) * std::sin(observer.lat);
        auto equatorialComponent = std::cos(u)
            + (observerAltitude / Earth::equatorialRadius) * std::cos(observer.lat);

        auto sinParallax = 4.2634515103856459e-05 / (distance / kAstronomicalUnit);
        auto cosDec = std::cos(equatorial.lat);
        auto cosHA = std::cos(hourAngle);
        auto deltaLng = std::atan2(-equatorialComponent * sinParallax * std::sin(hourAngle),
            cosDec - equatorialComponent * sinParallax * cosHA);
        auto topoLat = std::atan2(
            (std::sin(equatorial.lat) - polarAxisComponent * sinParallax) * cos(deltaLng),
            std::cos(equatorial.lat) - equatorialComponent * sinParallax * cosHA);
        return {
            topoLat,
            equatorial.lng + deltaLng
        };
    }

    static double getGeocentricElongation(const Coordinates &sunEquatorial,
            const Coordinates &objectEquatorial)
    {
        auto dec = objectEquatorial.lat;
        auto ra = objectEquatorial.lng;
        auto cosGeocentricElongation = std::sin(sunEquatorial.lat) * std::sin(dec)
            + std::cos(sunEquatorial.lat) * std::cos(dec) * std::cos(sunEquatorial.lng - ra);
        return std::acos(clamp(cosGeocentricElongation, 1.));
    }

    static double getDistanceFromObserver(double altitudeAngle, double geocentricDistance) {
        auto cosine = std::cos(altitudeAngle + M_PI / 2);
        auto c = geocentricDistance;
        auto b = Earth::equatorialRadius;

        // inverse law of cosines
        return b * cosine + std::sqrt(b * b * (cosine * cosine - 1) + c * c);
    }

    template <typename CelestialBody, bool switchToFindingNadir = false>
    static std::optional<TT> getTransitTT(TT tt, const DegreesCoordinates &observer) {
        auto getAltitudeAngle = [&observer](TT tt) {
            auto observation = Sky::observeInTT<CelestialBody>(tt, observer, true);
            return observation.altitudeAngle();
        };
        auto altitudeAngleAnalyzer = FunctionAnalyzer{std::move(getAltitudeAngle)};
        auto maxSample =
            altitudeAngleAnalyzer.template gradientWalkToExtremum<!switchToFindingNadir>(
                tt, kMaxDaysInAYear / 2,
                kTimesFindingStep, kTimesFindingPrecision);
        if (!maxSample) {
            return std::nullopt;
        }
        return maxSample->arg;
    }

    template <typename CelestialBody>
    static std::optional<TT> nadirFromTransitTT(TT transit, const DegreesCoordinates &observer,
            bool nextNadir)
    {
        auto mult = nextNadir ? 1 : -1;
        return getTransitTT<CelestialBody, true>(
            transit + mult * CelestialBody::approxDayLength() / 2, observer);
    }

    template <typename CelestialBody>
    static std::optional<UTC> getTransitUTC(UTC utc, const DegreesCoordinates &observer) {
        return utcFromTT(getTransitTT<CelestialBody>(ttFromUTC(utc), observer));
    }

    template <typename CelestialBody>
    struct TrajectoryCycle {
        static constexpr double kUpperEdgeBodySizeCorrection = CelestialBody::meanAngularSize() / 2;

        DegreesCoordinates observer;
        std::optional<TT> beginNadir;
        std::optional<TT> transit{};
        std::optional<TT> endNadir{};

        TrajectoryCycle(TT tt, const DegreesCoordinates &observer) : observer(observer) {
            transit = getTransitTT<CelestialBody>(tt, observer);
            if (!transit.has_value()) { return; }
            beginNadir = nadirFromTransitTT<CelestialBody>(*transit, observer, false);
            endNadir = nadirFromTransitTT<CelestialBody>(*transit, observer, true);
        }

        TrajectoryCycle() noexcept = default;
        TrajectoryCycle(TrajectoryCycle &&other) noexcept = default;
        TrajectoryCycle(const TrajectoryCycle &other) = default;
        TrajectoryCycle &operator=(TrajectoryCycle &&other) noexcept = default;
        TrajectoryCycle(
                const DegreesCoordinates &observer, std::optional<TT> firstNadir,
                std::optional<TT> transit, std::optional<TT> secondNadir)
                : observer(observer) , beginNadir(firstNadir)
                , transit(transit), endNadir(secondNadir) {}

        TrajectoryCycle getShiftedCycle(int shiftCycles) const {
            if (!isValid()) { return {}; }
            auto shiftedTransit = getTransitTT<CelestialBody>(
                *transit + shiftCycles * CelestialBody::approxDayLength(), observer);
            auto shiftedFirstNadir = shiftCycles == 1 ? endNadir
                : nadirFromTransitTT<CelestialBody>(*transit, observer, false);
            auto shiftedSecondNadir = shiftCycles == -1 ? beginNadir
                : nadirFromTransitTT<CelestialBody>(*transit, observer, true);
            return {observer, shiftedFirstNadir, shiftedTransit, shiftedSecondNadir};
        }

        bool isValid() const {
            return transit.has_value() && beginNadir.has_value() && endNadir.has_value();
        }

        std::optional<TT> getTTFromAngle(Rads angle, bool firstHalfOfTheCycle,
                double bodySizeCorrection = 0, bool accountForRefraction = true) const
        {
            if (!isValid()) { return std::nullopt; }
            auto altitudeAngleAnalyzer = FunctionAnalyzer{[this, bodySizeCorrection, angle, accountForRefraction](TT tt) {
                auto observation = Sky::observeInTT<CelestialBody>(tt, observer, false);
                observation.horizontal.lat += bodySizeCorrection;
                if (accountForRefraction) {
                    observation.accountForRefraction();
                }
                return observation.altitudeAngle() - angle;
            }};

            return altitudeAngleAnalyzer.getSignChangeToPositive(
                firstHalfOfTheCycle ? *beginNadir : *endNadir,
                *transit, kTimesFindingPrecision, 100);
        }

        std::optional<TT> getClosestTTWithRadsAltitudeAngle(Rads angle,
                bool firstHalfOfTheCycle, int shiftUpTo = 0, double bodySizeCorrection = 0,
                bool accountForRefraction = true) const
        {
            shiftUpTo = std::abs(shiftUpTo);
            auto trajectory = *this;
            int step = firstHalfOfTheCycle ? -1 : 1;
            int sanityCheck = kMaxDaysInAYear;
            while (true) {
                if (!trajectory.isValid()) { return std::nullopt; }
                if (sanityCheck-- <= 0) { return std::nullopt; }
                auto tt = trajectory.getTTFromAngle(angle, firstHalfOfTheCycle, bodySizeCorrection, accountForRefraction);
                if (tt) { return *tt; }
                if (shiftUpTo == 0) { return std::nullopt; }
                --shiftUpTo;
                trajectory = trajectory.getShiftedCycle(step);
            }
            return std::nullopt;
        }

        inline std::optional<UTC> getClosestUTCWithDegreesAltitudeAngle(Degrees angle,
                bool beforeTransit, int shiftUpTo = 0,
                bool caculateForTheUpperEdge = true,
                bool accountForRefraction = true) const
        {
            return utcFromTT(getClosestTTWithRadsAltitudeAngle(
                radsFromDegrees(angle), beforeTransit, shiftUpTo,
                caculateForTheUpperEdge ? kUpperEdgeBodySizeCorrection : 0.,
                accountForRefraction
            ));
        }

        inline std::optional<UTC> getTodaysSetUTC() const {
            return getClosestUTCWithDegreesAltitudeAngle(0, false, kCyclesToSearchForSetRise);
        }

        inline std::optional<UTC> getTodaysRiseUTC() const {
            return getClosestUTCWithDegreesAltitudeAngle(0, true, kCyclesToSearchForSetRise);
        }

        // rise or set
        inline std::optional<UTC> getClosestHorizonCrossingUTC(
                bool rise, bool forwardInTime)
        {
            return getClosestUTCWithDegreesAltitudeAngle(
                0, rise, (forwardInTime ? 1 : -1) * kMaxDaysInAYear);
        }

        inline std::optional<UTC> getTransitUTC() const {
            return utcFromTT(transit);
        }

        inline std::optional<UTC> getBeginNadir() const {
            return utcFromTT(beginNadir);
        }

        inline std::optional<UTC> getEndNadir() const {
            return utcFromTT(endNadir);
        }

        inline Days getLengthDays() const {
            if (!isValid()) { return 0.; }
            return *endNadir - *beginNadir;
        }
    };

    template <typename CelestialBody>
    static TrajectoryCycle<CelestialBody> getTrajectoryCycleFromUTC(UTC utc,
            const DegreesCoordinates &observer)
    {
        return {ttFromUTC(utc), observer};
    }

    template <typename CelestialBody>
    static inline std::optional<UTC> getTodaysSetUTC(UTC utc, const DegreesCoordinates &observer) {
        return getTrajectoryCycleFromUTC<CelestialBody>(utc, observer).getTodaysSetUTC();
    }

    template <typename CelestialBody>
    static inline std::optional<UTC> getTodaysRiseUTC(UTC utc, const DegreesCoordinates &observer) {
        return getTrajectoryCycleFromUTC<CelestialBody>(utc, observer).getTodaysRiseUTC();
    }

    template <typename CelestialBody>
    static inline std::optional<UTC> getClosestToCurrentTransitHorizonCrossingUTC(
            UTC utc, const DegreesCoordinates &observer, bool rise, bool forwardInTime)
    {
        return getTrajectoryCycleFromUTC<CelestialBody>(utc, observer).getClosestHorizonCrossingUTC(
            rise, forwardInTime);
    }

    struct NotableTimes {
        std::optional<UTC> rise;
        std::optional<UTC> transit;
        std::optional<UTC> set;
    };

    template <typename CelestialBody>
    static inline NotableTimes getNotableTimes(UTC utc, const DegreesCoordinates &observer) {
        auto trajectoryCycle = getTrajectoryCycleFromUTC<CelestialBody>(utc, observer);
        if (!trajectoryCycle.isValid()) { return {}; }
        return {trajectoryCycle.getTodaysRiseUTC(),
                trajectoryCycle.getTransitUTC(),
                trajectoryCycle.getTodaysSetUTC()};
    }
};

} // namespace Skygazing
