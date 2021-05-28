#ifndef SKYGAZING_TIME_H
#define SKYGAZING_TIME_H

#include <sstream>
#include <iomanip>
#include <chrono>
#include <optional>
#include <cmath>
#include <cstring>

namespace Skygazing {

using Days = double;
using OriginalJulian = Days; // Julian Day Number https://en.wikipedia.org/wiki/Julian_day
using Julian = Days; // Julian days since kJulianYear2000
using TT = Days; // Terrestrial Time
using Seconds = double;
using UTC = Seconds; // Universal Coordinated Time in seconds since 1970 = time_since_epoch()
using JulianCenturies = double;
constexpr Seconds kSecondsInMinute = 60;
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

template<bool useCache, bool useValidation, bool useManualUnwind = true>
struct TimestampFormatter {
    static constexpr int kUnitWidth = 2, kYearWidth = 4;
    static constexpr int kDayOnlyWidth = 10, kWithoutMillisecondsWidth = 20;
    static constexpr bool isDigit(char c) {
        return c >= '0' && c <= '9';
    }

    std::tm cachedDate = {};
    const char *cursor = nullptr;

    std::array<char, kDayOnlyWidth> cachedDayString = {};
    time_t cachedDayStart = {};
    bool errorFlag = false;

    struct DecimalComposer {
        int value = 0;

        void push(char digit) { value = value * 10 + (digit - '0'); }
        operator int() const { return value; }
    };

    void skip(std::string_view token) {
        if (useValidation && std::memcmp(cursor, token.begin(), token.size())) {
            errorFlag = true;
        }
        cursor += token.size();
    }

    char scanChar() {
        return *cursor++;
    }

    int scanInt(int length, std::string_view delimiter = "") {
        DecimalComposer parsed;
        while (length-- > 0) {
            auto c = scanChar();
            parsed.push(c);

            if (useValidation && !isDigit(c)) {
                errorFlag = true;
            }
        }
        skip(delimiter);
        return parsed;
    }

    std::optional<UTC> utcFromString(std::string_view string) {
        if (string.size() < kWithoutMillisecondsWidth) {
            return std::nullopt;
        }
        errorFlag = false;

        time_t time{};
        cursor = string.begin();
        if (useCache && !std::memcmp(cachedDayString.begin(), cursor, kDayOnlyWidth)) {
            cursor += kDayOnlyWidth + 1;
            time = cachedDayStart;
        } else {
            auto year = scanInt(kYearWidth, "-");
            auto month = scanInt(kUnitWidth, "-");
            auto day = scanInt(kUnitWidth, "T");

            cachedDate.tm_mday = day;
            cachedDate.tm_mon = month - 1;
            cachedDate.tm_year = year - 1900;
            time = timegm(&cachedDate);

            if (useCache) {
                cachedDayStart = time;
                std::memcpy(cachedDayString.begin(), string.begin(), kDayOnlyWidth);
            }
        }

        if (useValidation || !useManualUnwind) {
            auto hour = scanInt(kUnitWidth, ":");
            auto minute = scanInt(kUnitWidth, ":");
            auto second = scanInt(kUnitWidth, "");
            time += second + kSecondsInMinute * minute + kSecondsInHour * hour;
        } else {
            auto toDigit =
                [iter = cursor](int offset) { return static_cast<int>(iter[offset] - '0'); };
            time += kSecondsInHour * 10 * toDigit(0) + kSecondsInHour * toDigit(1)
                + kSecondsInMinute * 10 * toDigit(3) + kSecondsInMinute * toDigit(4)
                + 10 * toDigit(6) + toDigit(7);
            cursor += 8;
        }

        Seconds millisecondsPart = 0;
        if (string.size() > kWithoutMillisecondsWidth && scanChar() == '.') {
            Seconds milliUnit = 1;
            while (cursor != string.end()) {
                auto c = scanChar();
                if (!isDigit(c)) { break; }
                millisecondsPart += (milliUnit /= 10) * (c - '0');
            }
        }

        cursor = nullptr;
        if (useValidation && errorFlag) {
            return std::nullopt;
        }
        return static_cast<UTC>(time) + millisecondsPart;
    }

    static constexpr int unwindingPow10(int k) { return k > 0 ? 10 * unwindingPow10(k - 1) : 1; }

    template<int width = kUnitWidth>
    static void pushInt(char *&iter, int value, char c = 0) {
        for (auto unit = unwindingPow10(width - 1); unit; unit /= 10) {
            auto digit = value / unit;
            value -= digit * unit;
            *iter++ = static_cast<char>('0' + digit);
        }
        if (c) { *iter++ = c; }
    }

    std::string stringFromUTC(UTC utc) {
        auto wholeSeconds = static_cast<time_t>(std::floor(utc));
        auto date = *gmtime(&wholeSeconds);

        constexpr int kDotWidth = 1;
        constexpr int kMillisecondsWidth = 3;
        constexpr int kTerminatingZeroWidth = 1;
        constexpr int kResultWidth =
            kWithoutMillisecondsWidth + kDotWidth + kMillisecondsWidth + kTerminatingZeroWidth;
        std::array<char, kResultWidth> result{};
        auto iter = &result[0];

        auto year = date.tm_year + 1900;
        auto megaUnit = unwindingPow10(kYearWidth);
        auto yearOverflow = year / megaUnit;
        year -= yearOverflow * megaUnit;
        pushInt<kYearWidth>(iter, year, '-');
        pushInt(iter, date.tm_mon + 1, '-');
        pushInt(iter, date.tm_mday, 'T');
        pushInt(iter, date.tm_hour, ':');
        pushInt(iter, date.tm_min, ':');
        pushInt(iter, date.tm_sec);

        if (auto milliseconds = static_cast<int>((utc - wholeSeconds) * 1.0e3)) {
            *iter++ = '.';
            pushInt<kMillisecondsWidth>(iter, milliseconds);
        }
        *iter++ = 'Z';
        *iter++ = '\0';
        if (yearOverflow) {
            return std::to_string(yearOverflow) + &result[0];
        }
        return &result[0];
    }
};

inline std::optional<UTC> utcFromDateString(const std::string &dateString, int timezone = 0) {
    static TimestampFormatter<true, true, true> formatter{};
    auto localTime = formatter.utcFromString(dateString);
    if (!localTime) {
        std::cout << dateString << std::endl;
    }

    return localTime
        ? std::optional<UTC>(localTime.value() - timezone * kSecondsInHour)
        : std::optional<UTC>{};
}

inline std::string dateStringFromUTC(UTC utc, int timezone = 0) {
    static TimestampFormatter<true, true, true> formatter{};
    return formatter.stringFromUTC(utc + timezone * kSecondsInHour);
}

} // namespace Skygazing

#endif //SKYGAZING_TIME_H
