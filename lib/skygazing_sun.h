#ifndef SKYGAZING_SUN_H
#define	SKYGAZING_SUN_H

#include "skygazing_math.h"
#include "skygazing_time.h"
#include "skygazing_sky.h"

#include <iostream>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <vector>
#include <map>
#include <cassert>
#include <array>

namespace Skygazing {

struct Sun {
    static constexpr Rads meanAngularSize() { return 0.5_deg; }
    static constexpr Days approxDayLength() { return 1; }

    static Sky::EclipticPosition getEclipticPosition(TT tt) { // [AA] chapter 25
        auto julianCenturies = toJulianCenturies(tt);
        auto h = Horner(julianCenturies);
        Rads meanLng = h(280.46646_deg, 36000.76983_deg, 0.0003032_deg);
        Rads meanAnomaly = getSolarMeanAnomaly(tt);
        Rads equationOfTheCenter =
            h(1.914602_deg, -0.004817_deg, -0.000014_deg)  * std::sin(meanAnomaly)
            + h(0.019993_deg, -0.000101_deg) * std::sin(2 * meanAnomaly)
            + 0.000289_deg * std::sin(3 * meanAnomaly);
        Rads approxLat = 0;
        Rads lng = meanLng + equationOfTheCenter;
        auto omega = h(125.04_deg, -1934.136);
        Rads apparentLng = lng - 0.00569_deg - 0.00478_deg * std::sin(omega);
        return {{approxLat, normalizeRads(apparentLng)}, getDistanceFromMeanAnomaly(meanAnomaly)};
    }

    static Rads getSolarMeanAnomalyFromJulianCenturies(JulianCenturies julianCentury) {
        // [AA] (25.3), approximation
        return Horner(julianCentury)(357.5291092_deg, 35999.0502909_deg);
    }

    static Rads getSolarMeanAnomaly(TT tt) {
        return getSolarMeanAnomalyFromJulianCenturies(toJulianCenturies(tt));
    }

    static double getDistanceFromMeanAnomaly(Rads meanAnomaly) {
        constexpr auto astronomicalUnit = 149597870700.0;
        constexpr auto eccentricity = 0.0167086;
        constexpr auto semimajor = 1.000001018;
        auto r = semimajor * (1 - eccentricity * eccentricity)
            / (1 + eccentricity * std::cos(meanAnomaly));
        return astronomicalUnit * r;
    }
};

} // namespace Skygazing

#endif // SKYGAZING_SUN_H
