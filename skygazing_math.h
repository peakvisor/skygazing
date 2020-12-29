#ifndef SKYGAZING_MATH_H
#define SKYGAZING_MATH_H

#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include <utility>
#include <optional>

namespace Skygazing {

using Seconds = double; // seconds since 1970
using Rads = double; // Radians
using Degrees = double;
constexpr Rads kRadsPerDegree = M_PI / 180;

constexpr Degrees degreesFromRads(Rads rads) {
    return rads / kRadsPerDegree;
}

constexpr Rads radsFromDegrees(Degrees degrees) {
    return degrees * kRadsPerDegree;
}

constexpr double clamp(double value, double maxAbs = 1) { // too fast to check that maxAbs > 0
    return std::min(std::max(value, -maxAbs), +maxAbs);
}

template <typename Double>
constexpr Double normalizeAdding(Double value, Double min, Double turn) {
    while (value < min) {
        value += turn;
    }

    while (value >= min + turn) {
        value -= turn;
    }
    return value;
}

template <typename Double>
inline Double normalizeDividing(Double value, Double oneTurn) {
    auto turnsCount = std::floor(value / oneTurn);
    return value - turnsCount * oneTurn;
}

constexpr Rads normalizeRads(Rads rads) {
    return normalizeAdding(rads, -M_PI, 2 * M_PI);
}

inline Rads normalizeHugeRads(Rads rads) {
    return normalizeDividing(rads + M_PI, 2 * M_PI) - M_PI;
}

constexpr Degrees operator"" _deg(long double d) { return radsFromDegrees(d); }
constexpr Degrees operator"" _deg(unsigned long long d) { return radsFromDegrees(d); }

struct DegreesCoordinates {
    Degrees lat;
    Degrees lng;
};

struct Coordinates {
    Rads lat; // latitude
    Rads lng; // longitude

    [[nodiscard]] Coordinates normalize() const {
        return {normalizeHugeRads(lat),normalizeHugeRads(lng)};
    }
};

constexpr DegreesCoordinates degreesFromRads(const Coordinates &coordinates) {
    return {degreesFromRads(coordinates.lat), degreesFromRads(coordinates.lng)};
}

constexpr inline Coordinates radsFromDegrees(const DegreesCoordinates &degreesCoordinates) {
    return {radsFromDegrees(degreesCoordinates.lat), radsFromDegrees(degreesCoordinates.lng)};
}

inline auto haversineDistance(const Coordinates &from, const Coordinates &to) {
    auto latSin = std::sin((from.lat - to.lat) / 2);
    auto lngSin = std::sin((from.lng - to.lng) / 2);
    auto asinArg = std::sqrt(latSin * latSin + std::cos(from.lat) * std::cos(to.lat) * lngSin * lngSin);

    return 2 * std::asin(clamp(asinArg));
}

template <typename Number>
struct Horner {
    explicit constexpr Horner(Number x) : x(x) {}

    Number x;

    template <typename ...Coeffs>
    constexpr Number operator()(Number k0, Coeffs ...coeffs) const {
        return k0 + x * (*this)(coeffs...);
    }

    constexpr Number operator()(Number kn) const {
        return kn;
    }
};

template <typename Func>
struct FunctionAnalyzer {
    using Argument = double;
    using Value = double;

    struct Sample {
        Argument arg;
        Value val;
    };

    const Func func;

    explicit FunctionAnalyzer(Func &&func_) : func(std::move(func_)) {}

    Sample getSample(Argument arg) const {
        return {arg, func(arg)};
    }

    std::optional<Argument> findSignChangeToPositive(Argument fromArg, Argument toArg, Argument step,
                                                     Argument precision, int maxIterations = 100) const
    {
        if (func(fromArg) > 0) {
            return std::nullopt;
        }
        auto steps = static_cast<int>(std::ceil(std::abs((toArg - fromArg) / step)));
        step = (toArg - fromArg) / steps;

        Argument arg = fromArg;
        for (int i = 1; i < steps; ++i) {
            arg += step;
            if (func(arg) >= 0) {
                return getSignChangeToPositive(arg - step, arg, precision, maxIterations);
            }
        }
        return std::nullopt;
    }

    std::optional<Argument> getSignChangeToPositive(Argument leftArg, Argument rightArg, Argument precision,
                                                    int maxIterations) const
    {
        auto left = getSample(leftArg);
        auto right = getSample(rightArg);
        if (left.val > 0 || right.val < 0) {
            throw std::invalid_argument("potentially no root in getSignChangeToPositive");
        }
        for (int i = 0; i < maxIterations; ++i) {
            if (std::abs(right.arg - left.arg) < precision) {
                break;
            }
            auto mid = getSample((left.arg + right.arg) / 2);
            if (mid.val >= 0) {
                std::exchange(right, mid);
            } else {
                std::exchange(left, mid);
            }
        }
        return (left.arg + right.arg) / 2;
    }

    static bool potentialMax(Value left, Value mid, Value right) {
        return mid > left && mid > right;
    }

    std::optional<Sample> findMax(Argument fromArg, Argument toArg, Argument step, Argument precision,
                                  int maxIterations = 100) const
    {
        auto left = getSample(fromArg);
        auto mid = getSample(fromArg + step);
        auto right = getSample(mid.arg + step);

        while (mid.arg < toArg) {
            if (potentialMax(left.val, mid.val, right.val)) {
                return getMax(left, mid, right, precision, maxIterations);
            }
            left = std::exchange(mid, right);
            right = getSample(mid.arg + step);
        }
        return std::nullopt;
    }

    std::optional<Sample> getMax(Sample left, Sample mid, Sample right, Argument precision,
                                 int maxIterations = 100) const
    {
        for (int i = 0; i < maxIterations; ++i) {
            if (right.arg - left.arg <= precision) {
                break;
            }
            auto midLeft = getSample((mid.arg + left.arg) / 2);
            if (potentialMax(left.val, midLeft.val, mid.val)) {
                right = std::exchange(mid, midLeft);
                continue;
            }
            auto midRight = getSample((mid.arg + right.arg) / 2);
            if (potentialMax(mid.val, midRight.val, right.val)) {
                left = std::exchange(mid, midRight);
                continue;
            }
            left = midLeft;
            right = midRight;
        }
        return mid;
    }
};

} // namespace Skygazing

#endif //SKYGAZING_MATH_H
