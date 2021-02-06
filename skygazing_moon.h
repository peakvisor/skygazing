#ifndef SKYGAZING_MOON_H
#define SKYGAZING_MOON_H

#include "skygazing_math.h"
#include "skygazing_time.h"
#include "skygazing_sky.h"
#include "skygazing_sun.h"

#include <array>
#include <cmath>

namespace Skygazing {

struct Moon {
    static constexpr Rads meanAngularSize() { return 0.5_deg; }
    static constexpr Days approxDayLength() { return 1.0345; }

    struct MoonDetails {
        double geocentricElongation;
        double angle;
        double fraction;
        double phase;
        double crescentApparentRotation;
    };

    static MoonDetails getMoonDetails(const Sky::CelestialObjectObservation &moon,
            const Sky::CelestialObjectObservation &sun)
    {
        MoonDetails details{};

        double elongation = Sky::getGeocentricElongation(sun.equatorial, moon.equatorial);
        details.geocentricElongation = elongation;

        double inc = std::atan2(sun.distance * std::sin(elongation),
            moon.distance - sun.distance * std::cos(elongation));
        details.fraction = (1 + std::cos(inc)) / 2;

        auto [sunDec, sunRa] = sun.equatorial;
        auto [dec, ra] = moon.equatorial;
        details.angle = std::atan2(cos(sunDec) * std::sin(sunRa - ra),
            std::sin(sunDec) * std::cos(dec)
                - std::cos(sunDec) * std::sin(dec) * std::cos(sunRa - ra));

        details.phase = 0.5 + 0.5 * inc * (details.angle < 0 ? -1 : 1) / M_PI;
        details.crescentApparentRotation = getCrescentApparentRotation(moon, sun, details.phase);
        return details;
    }

    static inline Rads getCrescentApparentRotation(const Sky::CelestialObjectObservation &moon,
            const Sky::CelestialObjectObservation &sun, double moonPhase)
    {
        auto b = M_PI_2 - moon.horizontal.lat, c = M_PI_2 - sun.horizontal.lat;
        auto a = haversineDistance(moon.horizontal, sun.horizontal);
        auto denominator = std::sin(a) * std::sin(b);
        auto rotation = std::abs(denominator) < std::numeric_limits<double>::epsilon() ? 0.0
            : std::acos((std::cos(c) - std::cos(a) * std::cos(b)) / denominator);

        auto deltaLng = normalizeRads(sun.horizontal.lng - moon.horizontal.lng);
        if (deltaLng < 0) {
            rotation = -rotation;
        }

        if (moonPhase > 0.5) {
            rotation += M_PI;
        }

        return rotation;
    }

    struct BaseParameter {
        double meanElongation;
        double solarMeanAnomaly;
        double moonMeanAnomaly;
        double moonArgumentOfLatitude;

        constexpr size_t eccentricityPower() const {
            auto k = static_cast<int>(solarMeanAnomaly);
            return k > 0 ? +k : -k; // aka constexpr std::abs
        }

        template <typename Container>
        constexpr static bool validateEccentricityPowers(const Container &coeffs, size_t limit) {
            for (auto &coef : coeffs) {
                if (coef.eccentricityPower() >= limit) { return false; }
            }
            return true;
        }

        constexpr double dot(const BaseParameter &other) const {
            return meanElongation * other.meanElongation
                + solarMeanAnomaly * other.solarMeanAnomaly
                + moonMeanAnomaly * other.moonMeanAnomaly
                + moonArgumentOfLatitude * other.moonArgumentOfLatitude;
        }
    };

    // [AA] chapter 47
    static Sky::EclipticPosition getEclipticPosition(TT tt) {
        auto julianCenturies = toJulianCenturies(tt);
        auto h = Horner(julianCenturies);

        auto meanLng = normalizeHugeRads(h(218.3164477_deg, 481267.88123421_deg, -0.0015786_deg,
            1_deg / 538841, -1_deg / 65194000));
        auto base = baseParameterFromJulianCenturies(julianCenturies);
        auto A1 = h(119.75_deg, 131.849_deg);
        auto A2 = h(53.09_deg, 479264.29_deg);
        auto A3 = h(313.45_deg, 481266.484_deg);
        auto eccentricityMultiplier = h(1, -0.002516, -0.0000074);
        std::array ePowers{1.0, eccentricityMultiplier,
            eccentricityMultiplier * eccentricityMultiplier};
        static_assert(decltype(base)::validateEccentricityPowers(latCoeffs, std::size(ePowers)),
            "Unexpected lat coefficient");
        static_assert(decltype(base)::validateEccentricityPowers(lngCoeffs, std::size(ePowers)),
            "Unexpected lng coefficient");

        auto distanceDelta = 0.0;
        auto lngDelta = 3958 * std::sin(A1) + 1962 * std::sin(meanLng - base.moonArgumentOfLatitude)
            + 318 * std::sin(A2);
        auto lat = -2235 * std::sin(meanLng) + 382 * std::sin(A3)
            + 175 * std::sin(A1 - base.moonArgumentOfLatitude)
            + 175 * std::sin(A1 + base.moonArgumentOfLatitude)
            + 127 * std::sin(meanLng - base.moonMeanAnomaly)
            - 115 * std::sin(meanLng + base.moonMeanAnomaly);

        for (auto &coeff : lngCoeffs) {
            auto a = base.dot(coeff);
            lngDelta += coeff.sin * std::sin(a) * ePowers[coeff.eccentricityPower()];
            distanceDelta += coeff.cos * std::cos(a) * ePowers[coeff.eccentricityPower()];
        }
        for (auto &coeff : latCoeffs) {
            auto b = base.dot(coeff);
            lat += coeff.sin * std::sin(b) * ePowers[coeff.eccentricityPower()];
        }

        lngDelta *= 1e-6_deg;
        lat *= 1e-6_deg;
        return {
            {normalizeHugeRads(lat), normalizeHugeRads(meanLng + lngDelta)},
            385000560 + distanceDelta
        };
    }

    static BaseParameter baseParameterFromJulianCenturies(JulianCenturies julianCenturies) {
        auto h = Horner(julianCenturies);
        return {
            normalizeHugeRads(h(297.8501921_deg, 445267.1114034_deg, -0.0018819_deg, 1_deg / 545868,
                -1_deg / 113065000)),
            normalizeHugeRads(Sun::getSolarMeanAnomalyFromJulianCenturies(julianCenturies)),
            normalizeHugeRads(h(134.9633964_deg, 477198.8675055_deg, 0.0087414_deg, 1_deg / 69699,
                -1_deg / 14712000)),
            normalizeHugeRads(h(93.272095_deg, 483202.0175233_deg, -0.0036539_deg, -1_deg / 3526000,
                1_deg / 863310000))
        };
    }

    struct LatitudeCoeff : BaseParameter {
        double sin;
        double cos;
    };

    struct LongitudeCoeff : BaseParameter {
        double sin;
    };

    static constexpr std::array<LatitudeCoeff, 60> lngCoeffs{{
        {0, 0, 1, 0, 6288774, -20905355},
        {2, 0, -1, 0, 1274027, -3699111},
        {2, 0, 0, 0, 658314, -2955968},
        {0, 0, 2, 0, 213618, -569925},

        {0, 1, 0, 0, -185116, 48888},
        {0, 0, 0, 2, -114332, -3149},
        {2, 0, -2, 0, 58793, 246158},
        {2, -1, -1, 0, 57066, -152138},

        {2, 0, 1, 0, 53322, -170733},
        {2, -1, 0, 0, 45758, -204586},
        {0, 1, -1, 0, -40923, -129620},
        {1, 0, 0, 0, -34720, 108743},

        {0, 1, 1, 0, -30383, 104755},
        {2, 0, 0, -2, 15327, 10321},
        {0, 0, 1, 2, -12528, 0},
        {0, 0, 1, -2, 10980, 79661},

        {4, 0, -1, 0, 10675, -34782},
        {0, 0, 3, 0, 10034, -23210},
        {4, 0, -2, 0, 8548, -21636},
        {2, 1, -1, 0, -7888, 24208},

        {2, 1, 0, 0, -6766, 30824},
        {1, 0, -1, 0, -5163, -8379},
        {1, 1, 0, 0, 4987, -16675},
        {2, -1, 1, 0, 4036, -12831},

        {2, 0, 2, 0, 3994, -10445},
        {4, 0, 0, 0, 3861, -11650},
        {2, 0, -3, 0, 3665, 14403},
        {0, 1, -2, 0, -2689, -7003},

        {2, 0, -1, 2, -2602, 0},
        {2, -1, -2, 0, 2390, 10056},
        {1, 0, 1, 0, -2348, 6322},
        {2, -2, 0, 0, 2236, -9884},

        {0, 1, 2, 0, -2120, 5751},
        {0, 2, 0, 0, -2069, 0},
        {2, -2, -1, 0, 2048, -4950},
        {2, 0, 1, -2, -1773, 4130},

        {2, 0, 0, 2, -1595, 0},
        {4, -1, -1, 0, 1215, -3958},
        {0, 0, 2, 2, -1110, 0},
        {3, 0, -1, 0, -892, 3258},

        {2, 1, 1, 0, -810, 2616},
        {4, -1, -2, 0, 759, -1897},
        {0, 2, -1, 0, -713, -2117},
        {2, 2, -1, 0, -700, 2354},

        {2, 1, -2, 0, 691, 0},
        {2, -1, 0, -2, 596, 0},
        {4, 0, 1, 0, 549, -1423},
        {0, 0, 4, 0, 537, -1117},

        {4, -1, 0, 0, 520, -1571},
        {1, 0, -2, 0, -487, -1739},
        {2, 1, 0, -2, -399, 0},
        {0, 0, 2, -2, -381, -4421},

        {1, 1, 1, 0, 351, 0},
        {3, 0, -2, 0, -340, 0},
        {4, 0, -3, 0, 330, 0},
        {2, -1, 2, 0, 327, 0},

        {0, 2, 1, 0, -323, 1165},
        {1, 1, -1, 0, 299, 0},
        {2, 0, 3, 0, 294, 0},
        {2, 0, -1, -2, 0, 8752}
    }};

    static constexpr std::array<LongitudeCoeff, 60> latCoeffs{{
        {0, 0, 0, 1, 5128122},
        {0, 0, 1, 1, 280602},
        {0, 0, 1, -1, 277693},
        {2, 0, 0, -1, 173237},

        {2, 0, -1, 1, 55413},
        {2, 0, -1, -1, 46271},
        {2, 0, 0, 1, 32573},
        {0, 0, 2, 1, 17198},

        {2, 0, 1, -1, 9266},
        {0, 0, 2, -1, 8822},
        {2, -1, 0, -1, 8216},
        {2, 0, -2, -1, 4324},

        {2, 0, 1, 1, 4200},
        {2, 1, 0, -1, -3359},
        {2, -1, -1, 1, 2463},
        {2, -1, 0, 1, 2211},

        {2, -1, -1, -1, 2065},
        {0, 1, -1, -1, -1870},
        {4, 0, -1, -1, 1828},
        {0, 1, 0, 1, -1794},

        {0, 0, 0, 3, -1749},
        {0, 1, -1, 1, -1565},
        {1, 0, 0, 1, -1491},
        {0, 1, 1, 1, -1475},

        {0, 1, 1, -1, -1410},
        {0, 1, 0, -1, -1344},
        {1, 0, 0, -1, -1335},
        {0, 0, 3, 1, 1107},

        {4, 0, 0, -1, 1021},
        {4, 0, -1, 1, 833},

        {0, 0, 1, -3, 777},
        {4, 0, -2, 1, 671},
        {2, 0, 0, -3, 607},
        {2, 0, 2, -1, 596},

        {2, -1, 1, -1, 491},
        {2, 0, -2, 1, -451},
        {0, 0, 3, -1, 439},
        {2, 0, 2, 1, 422},

        {2, 0, -3, -1, 421},
        {2, 1, -1, 1, -366},
        {2, 1, 0, 1, -351},
        {4, 0, 0, 1, 331},

        {2, -1, 1, 1, 315},
        {2, -2, 0, -1, 302},
        {0, 0, 1, 3, -283},
        {2, 1, 1, -1, -229},

        {1, 1, 0, -1, 223},
        {1, 1, 0, 1, 223},
        {0, 1, -2, -1, -220},
        {2, 1, -1, -1, -220},

        {1, 0, 1, 1, -185},
        {2, -1, -2, -1, 181},
        {0, 1, 2, 1, -177},
        {4, 0, -2, -1, 176},

        {4, -1, -1, -1, 166},
        {1, 0, 1, -1, -164},
        {4, 0, 1, -1, 132},
        {1, 0, -1, -1, -119},

        {4, -1, 0, -1, 115},
        {2, -2, 0, 1, 107}
    }};
};

}
#endif // SKYGAZING_MOON_H
