#ifndef SKYGAZING_TIME_H
#define SKYGAZING_TIME_H

#include <sstream>
#include <iomanip>

namespace Skygazing {

using OriginalJulian = double; // Julian Day Number https://en.wikipedia.org/wiki/Julian_day
using Julian = double; // Julian days since kJulianYear2000
constexpr Seconds kSecondsInHour = 3600;
constexpr Seconds kSecondsInDay = kSecondsInHour * 24;
constexpr OriginalJulian kJulianYear1970 = 2440588;
constexpr OriginalJulian kJulianYear2000 = 2451545;
constexpr Julian kJulianFrom1970To2000 = 2451545 - 2440588;
constexpr Julian kJulianCentury = 36525;

inline Seconds dateStringToSeconds(const std::string &dateString, int timezone = 0) {
    tm tmp{};
    std::istringstream ss(dateString);
    ss >> std::get_time(&tmp, "%Y-%m-%dT%H:%M:%SZ");
    return static_cast<Seconds>(timegm(&tmp)) - timezone * kSecondsInHour;
}

inline std::string secondsToDateString(Seconds secondsSince1970, int timezone = 0) {
    time_t t(static_cast<int>(std::floor(secondsSince1970)));
    t += timezone * kSecondsInHour;
    tm date = *gmtime(&t);
    static auto toFixedLengthString = [](const auto &value, size_t length) {
        auto withoutLeadingZeros = std::to_string(value);
        if (withoutLeadingZeros.length() > length) {
            throw std::invalid_argument("invalid value in toFixedLengthString: " + withoutLeadingZeros);
        }
        return std::string(length - withoutLeadingZeros.length(), '0') + withoutLeadingZeros;
    };
    return std::to_string(date.tm_year + 1900) + "-" + toFixedLengthString(date.tm_mon + 1, 2)
           + "-" + toFixedLengthString(date.tm_mday, 2) + "T"
           + toFixedLengthString(date.tm_hour, 2) + ":" + toFixedLengthString(date.tm_min, 2) + ":" + toFixedLengthString(date.tm_sec, 2) + "Z";
}

constexpr OriginalJulian secondsToJulianSinceJulianEpoch(Seconds time) noexcept {
    return time / kSecondsInDay - 0.5 + kJulianYear1970;
}

constexpr Julian secondsToJulian(Seconds seconds) noexcept {
    return seconds / kSecondsInDay - 0.5 - kJulianFrom1970To2000;
}

constexpr Seconds julianToSeconds(Julian julian) {
    return (julian + kJulianFrom1970To2000 + 0.5) * kSecondsInDay;
}

inline Seconds tmToSeconds(std::tm *date) {
    return timegm(date);
}

inline std::tm* secondsToTm(Seconds seconds) {
    auto time = static_cast<time_t>(std::round(seconds));
    return gmtime(&time);
}

inline Julian tmToJulian(std::tm *date) {
    return secondsToJulian(tmToSeconds(date));
}

inline Julian tmToJulianSinceJulianEpoch(tm *date) {
    time_t time = timegm(date);
    return secondsToJulianSinceJulianEpoch(time);
}

inline std::tm* julianToTm(Julian julian) {
    return secondsToTm(julianToSeconds(julian));
}

constexpr Julian toJulianCenturies(Julian julian) {
    return julian / kJulianCentury;
}

constexpr Julian estimate(Julian jd) { // works between 2005 and 2050
    double t = jd / 365.25;
    return  62.92 + 0.32217 * t + 0.005589 * t * t;
}

constexpr Julian toTerrestrialTime(Julian jd) {
    double deltaT = estimate(jd);
    return jd + deltaT / 86400;
}

constexpr Julian fromTerrestrialTime(Julian jd) {
    double deltaT = estimate(jd);
    return jd - deltaT / 86400;
}

constexpr Julian secondsToTerrestrialJulian(Seconds seconds) noexcept {
    return toTerrestrialTime(secondsToJulian(seconds));
}

constexpr Seconds terrestrialJulianToSeconds(Julian terrestrialJulian) {
    return julianToSeconds(fromTerrestrialTime(terrestrialJulian));
}

} // namespace Skygazing

#endif //SKYGAZING_TIME_H
