#ifndef SKYGAZING_MATH_H
#define SKYGAZING_MATH_H

#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include <utility>
#include <optional>

namespace Skygazing {

using Rads = double; // Radians
using Degrees = double;
constexpr Rads kRadsPerDegree = M_PI / 180;

constexpr Degrees degreesFromRads(Rads rads) {
    return rads / kRadsPerDegree;
}

constexpr Rads radsFromDegrees(Degrees degrees) {
    return degrees * kRadsPerDegree;
}

constexpr double clamp(double value, double maxAbs = 1) {
    // too fast to check that maxAbs > 0
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
    auto asinArg = std::sqrt(latSin * latSin
        + std::cos(from.lat) * std::cos(to.lat) * lngSin * lngSin);

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

    struct SampleTriangle {
        std::array<Sample, 3> samples;
        inline Sample &left() { return samples[0]; }
        const inline Sample &left() const { return samples[0]; }
        const inline Sample &mid() const { return samples[1]; }
        const inline Sample &right() const { return samples[2]; }

        void shift(const Func &func, Argument step) {
            auto i = step > 0 ? 2 : 0;
            auto newArg = samples[i].arg + step;
            samples[2 - i] = std::exchange(samples[1], samples[i]);
            samples[i] = {newArg, func(newArg)};
        }
    };

    const Func func;

    explicit FunctionAnalyzer(Func &&func_) : func(std::move(func_)) {}

    constexpr Sample getSample(Argument arg) const {
        return {arg, func(arg)};
    }

    constexpr Sample getMidSample(const Sample &left, const Sample &right) const {
        return getSample((left.arg + right.arg) / 2);
    }

    constexpr SampleTriangle getSampleTriangle(Argument mid, Argument step) const {
        return {{getSample(mid - step), getSample(mid), getSample(mid +  step)}};
    }

    std::optional<Argument> getSignChangeToPositive(Sample left, Sample right,
            Argument precision, int iterations = 100) const
    {
        if (left.val > 0 || right.val < 0) { return std::nullopt; }
        while (iterations-- > 0 && std::abs(right.arg - left.arg) > precision) {
            auto mid = getMidSample(left, right);
            (mid.val > 0 ? right : left) = mid;
        }
        return (left.arg + right.arg) / 2;
    }

    std::optional<Argument> getSignChangeToPositive(Argument leftArg, Argument rightArg,
            Argument precision, int iterations = 100) const
    {
        auto left = getSample(leftArg);
        auto right = getSample(rightArg);
        return getSignChangeToPositive(left, right, precision, iterations);
    }

    template <bool isMax>
    static constexpr bool hasExtremumBetween(Value left, Value mid, Value right) {
        return isMax ? mid > std::max(left, right) : mid < std::min(left, right);
    }

    template <bool isMax>
    static constexpr bool hasExtremumBetween(const Sample &left, const Sample &mid,
            const Sample &right) {
        return hasExtremumBetween<isMax>(left.val, mid.val, right.val);
    }

    template <bool isMax>
    static constexpr bool hasExtremumIn(const SampleTriangle &t) {
        return hasExtremumBetween<isMax>(t.left().val, t.mid().val, t.right().val);
    }

    template <bool isMax>
    inline Sample getExtremum(SampleTriangle t, Argument precision,
            int iterations = 100) const
    {
        while (iterations-- > 0 && std::abs(t.right().arg - t.left().arg) > precision) {
            auto newSample = getMidSample(t.left(), t.mid());
            if (hasExtremumBetween<isMax>(t.left(), newSample, t.mid())) {
                t.samples[2] = std::exchange(t.samples[1], newSample);
            } else {
                t.samples[0] = std::exchange(t.samples[2], newSample);
            }
        }
        return t.mid();
    }

    template <bool toMax = true>
    std::optional<Sample> gradientWalkToExtremum(Argument fromArg, Argument maxShift,
            Argument step, Argument precision, int maxIterations = 100) const
    {
        auto triangle = getSampleTriangle(fromArg, step);
        bool shouldGoLeft = (hasExtremumIn<!toMax>(triangle)
            && getExtremum<!toMax>(triangle, precision, maxIterations).arg > fromArg)
            || ((triangle.right().val < triangle.left().val) == toMax);
        step = std::copysign(step, shouldGoLeft ? -1 : +1);

        maxShift = std::abs(maxShift);
        while (std::abs(triangle.mid().arg - fromArg) < maxShift) {
            if (hasExtremumIn<toMax>(triangle)) {
                return getExtremum<toMax>(triangle, precision, maxIterations);
            }
            triangle.shift(func, step);
        }
        return std::nullopt;
    }
};

} // namespace Skygazing

#endif //SKYGAZING_MATH_H
