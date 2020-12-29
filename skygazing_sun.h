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
    static Sky::EclipticPosition getEclipticPosition(Julian julian) {
        Rads meanAnomaly = getSolarMeanAnomaly(julian);
        Rads equationOfTheCenter = 1.9148_deg * std::sin(meanAnomaly)
                                   + 0.02_deg * std::sin(2 * meanAnomaly)
                                   + 0.0003_deg * std::sin(3 * meanAnomaly);
        constexpr Rads kEarthPerihelionLongitude = 102.9372_deg;
        constexpr Rads kSkygazingOverfitting = 0.003;
        Rads eclipticLongitude = meanAnomaly + equationOfTheCenter + kEarthPerihelionLongitude + M_PI + kSkygazingOverfitting;

        return {{0., normalizeRads(eclipticLongitude)}, getDistanceFromMeanAnomaly(meanAnomaly)};
    }

    static Rads getSolarMeanAnomalyFromJulianCenturies(Julian julianCentury) {
        return Horner(julianCentury)(357.5291092_deg, 35999.0502909_deg); // [AA] (25.3), approximation
    }

    static Rads getSolarMeanAnomaly(Julian julian) {
        return getSolarMeanAnomalyFromJulianCenturies(toJulianCenturies(julian));
    }

    static double getDistanceFromMeanAnomaly(Rads meanAnomaly) {
        constexpr auto astronomicalUnit = 149597870700.0;
        constexpr auto eccentricity = 0.0167086;
        constexpr auto semimajor = 1.000001018;
        auto r = semimajor * (1 - eccentricity * eccentricity) / (1 + eccentricity * std::cos(meanAnomaly));
        return astronomicalUnit * r;
    }

    static constexpr Rads meanAngularSize() {
        return 0.5_deg;
    }
};

} // namespace Skygazing

#endif // SKYGAZING_SUN_H