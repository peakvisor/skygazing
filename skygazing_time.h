#ifndef SKYGAZING_TIME_H
#define SKYGAZING_TIME_H

#include <sstream>
#include <iomanip>
#include <chrono>

namespace Skygazing {

using Days = double;
using OriginalJulian = Days; // Julian Day Number https://en.wikipedia.org/wiki/Julian_day
using Julian = Days; // Julian days since kJulianYear2000
using TT = Days; // Terrestrial Time
using Seconds = double;
using UTC = Seconds; // Universal Coordinated Time in seconds since 1970 = time_since_epoch()
using JulianCenturies = double;
constexpr Seconds kSecondsInHour = 3600;
constexpr double kHoursInDay = 24.;
constexpr Seconds kSecondsInDay = kSecondsInHour * kHoursInDay;
constexpr OriginalJulian kJulianYear1970 = 2440588;
constexpr OriginalJulian kJulianYear2000 = 2451545;
constexpr Days kDaysFrom1970To2000 = 2451545 - 2440588;
constexpr Days kJulianCentury = 36525;
constexpr Days kMaxDaysInAYear = 366;

UTC getCurrentUTC() {
    auto timePointNow = std::chrono::high_resolution_clock::now().time_since_epoch();
    auto durationNow = std::chrono::duration_cast<std::chrono::milliseconds>(timePointNow);
    constexpr double kPrecision = 1.e-3;
    return static_cast<double>(durationNow.count()) * kPrecision;
}

inline UTC utcFromDateString(const std::string &dateString, int timezone = 0) {
    tm tmp{};
    std::istringstream ss(dateString);
    ss >> std::get_time(&tmp, "%Y-%m-%dT%H:%M:%SZ");
    return static_cast<Seconds>(timegm(&tmp)) - timezone * kSecondsInHour;
}

inline std::string dateStringFromUTC(UTC utc, int timezone = 0) {
    time_t t(static_cast<int>(std::floor(utc)));
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

constexpr OriginalJulian julianSinceEpochFromUTC(UTC utc) noexcept {
    return utc / kSecondsInDay - 0.5 + kJulianYear1970;
}

constexpr Julian julianFromUTC(UTC utc) noexcept {
    return utc / kSecondsInDay - 0.5 - kDaysFrom1970To2000;
}

constexpr UTC utcFromJulian(Julian julian) {
    return (julian + kDaysFrom1970To2000 + 0.5) * kSecondsInDay;
}

inline UTC tmToUTC(std::tm *date) {
    return timegm(date);
}

inline std::tm* utcToTm(UTC utc) {
    auto time = static_cast<time_t>(std::round(utc));
    return gmtime(&time);
}

constexpr JulianCenturies toJulianCenturies(Julian julian) {
    return julian / kJulianCentury;
}

constexpr Seconds estimate(Julian jd) { // works between 2005 and 2050
    double t = jd / 365.25;
    return  62.92 + 0.32217 * t + 0.005589 * t * t;
}

constexpr TT ttFromJulian(Julian jd) {
    double deltaT = estimate(jd);
    return jd + deltaT / 86400;
}

constexpr Julian julianFromTT(TT tt) {
    double deltaT = estimate(tt);
    return tt - deltaT / 86400;
}

constexpr TT ttFromUTC(UTC utc) noexcept {
    return ttFromJulian(julianFromUTC(utc));
}

constexpr std::optional<TT> ttFromUTC(const std::optional<UTC> &utc) noexcept {
    return utc ? std::optional<Julian>(ttFromUTC(*utc)) : std::nullopt;
}

constexpr UTC utcFromTT(TT tt) noexcept {
    return utcFromJulian(julianFromTT(tt));
}

constexpr std::optional<UTC> utcFromTT(const std::optional<TT> &tt) noexcept {
    return tt ? std::optional<UTC>(utcFromTT(*tt)) : std::nullopt;
}

} // namespace Skygazing

#endif //SKYGAZING_TIME_H
