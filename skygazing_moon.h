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
    struct MoonDetails {
        double geocentricElongation;
        double angle;
        double fraction;
        double phase;
        double crescentApparentRotation;
    };

    static MoonDetails getMoonDetails(const Sky::CelestialObjectObservation &moon, const Sky::CelestialObjectObservation &sun) {
        MoonDetails details;
        auto sunRa = sun.rightAscension();
        auto sunDec = sun.declination();
        auto moonRa = moon.rightAscension();
        auto moonDec = moon.declination();
        details.geocentricElongation = Sky::getGeocentricElongation(sun.equatorial, moon.equatorial);
        const auto &ge = details.geocentricElongation;
        double inc = std::atan2(sun.distance * std::sin(ge), moon.distance - sun.distance * std::cos(ge));
        details.angle = std::atan2(
                cos(sunDec) * std::sin(sunRa - moonRa),
                std::sin(sunDec) * std::cos(moonDec) - std::cos(sunDec) * std::sin(moonDec) * std::cos(sunRa - moonRa));

        details.fraction = (1 + std::cos(inc)) / 2;
        details.phase = 0.5 + 0.5 * inc * (details.angle < 0 ? -1 : 1) / M_PI;
        details.crescentApparentRotation = getCrescentApparentRotation(moon, sun, details.phase);
        return details;
    }

    static inline Rads getCrescentApparentRotation(const Sky::CelestialObjectObservation &moon,
                                                   const Sky::CelestialObjectObservation &sun,
                                                   double moonPhase) {
        auto b = M_PI_2 - moon.horizontal.lat, c = M_PI_2 - sun.horizontal.lat;
        auto a = haversineDistance(moon.horizontal, sun.horizontal);
        auto denominator = std::sin(a) * std::sin(b);
        auto rotation = std::abs(denominator) < std::numeric_limits<double>::epsilon()
                        ? 0.0
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

    // [AA] chapter 47
    static Sky::EclipticPosition getEclipticPosition(Julian julian) {
        double T = toJulianCenturies(julian);
        double L_ = normalizeHugeRads(Horner(T)(218.3164477_deg, 481267.88123421_deg, -0.0015786_deg, 1_deg/538841, -1_deg/65194000));
        double D = normalizeHugeRads(Horner(T)(297.8501921_deg, 445267.1114034_deg, -0.0018819_deg, 1_deg/545868, -1_deg/113065000));

        double M = Horner(T)(357.5291092_deg, 35999.0502909_deg, -0.0001535_deg, 1_deg / 24490000);
        double M_ = Horner(T)(134.9633964_deg, 477198.8675055_deg, 0.0087414_deg, 1_deg / 69699, -1_deg / 14712000);
        double F = Horner(T)(93.272095_deg, 483202.0175233_deg, -0.0036539_deg, -1_deg/3526000, 1_deg / 863310000);

        double A1 = 119.75_deg + 131.849_deg * T;
        double A2 = 53.09_deg + 479264.29_deg * T;
        double A3 = 313.45_deg + 481266.484_deg * T;
        double E = Horner(T)(1, -0.002516, -0.0000074);
        double E2 = E * E;
        double suml = 3958 * std::sin(A1) + 1962 * std::sin(L_-F) + 318 * std::sin(A2);
        double sumr = 0;
        double sumb = -2235 * std::sin(L_) + 382 * std::sin(A3) + 175 * std::sin(A1-F) +
                      175 * std::sin(A1 + F) + 127 * std::sin(L_ - M_) - 115 * std::sin(L_+M_);

        for (const auto &r : ta) {
            // 0:D, 1:M, 2:M_, 3:F, 4:suml, 5:sumr
            double a = D * r[0] + M * r[1] + M_ * r[2] + F * r[3];
            double sa = std::sin(a);
            double ca = std::cos(a);
            switch (int(r[1])) { // M
                case 0:
                    suml += r[4] * sa;
                    sumr += r[5] * ca;
                    break;
                case  1:
                case -1:
                    suml += r[4] * sa * E;
                    sumr += r[5] * ca * E;
                    break;
                case  2:
                case -2:
                    suml += r[4] * sa * E2;
                    sumr += r[5] * ca * E2;
                    break;
                default:
                    throw std::invalid_argument("incorrect M");
            }
        }

        for (const auto &r : tb) {
            // 0:D, 1:M, 2:M_, 3:F, 4:sumb
            double b = D * r[0] + M * r[1] + M_ * r[2] + F * r[3];
            double sb = std::sin(b);

            switch (int(r[1])) { // M
                case 0:
                    sumb += r[4] * sb;
                    break;
                case  1:
                case -1:
                    sumb += r[4] * sb * E;
                    break;
                case  2:
                case -2:
                    sumb += r[4] * sb * E2;
                    break;
                default:
                    throw std::invalid_argument("invalid M in tb");
            }
        }
        return {
                {normalizeHugeRads(sumb * 1e-6_deg),
                 normalizeHugeRads(L_ + suml * 1e-6_deg)},
                385000.56 + sumr * 1e-3
        };
    }

    static constexpr std::array<std::array<double,6>, 60> ta{{
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

    // 0:D, 1:M, 2:Mʹ, 3:F, 4:sumb
    static constexpr std::array<std::array<double, 5>, 60> tb{{
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
