#pragma once

#include "../lib/skygazing_math.h"
#include "../lib/skygazing_time.h"

#include <random>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <map>
#include <execution>

namespace Skygazing {

template <typename L, typename R>
void assertEqualityImpl(const L &lhs, const R &rhs, const std::string &lhsString,
        const std::string &rhsString, double precision)
{
    if (std::abs(lhs - rhs) >= precision) {
        std::cout << "FAIL: " << lhsString << " == " << lhs << " != " << rhs << " == " << rhsString
            << std::endl;
        abort();
    }
}

template <typename T>
void assertSmallImpl(const T &x, const std::string &xString, double precision) {
    if (std::abs(x) > precision) {
        std::cout << "FAIL: std::abs(" << xString << ") = " << std::abs(x) << " > " << precision
            << std::endl;
        abort();
    }
}

void outputExecutionTimeMessage(UTC executionStart, const std::string &taskName) {
    std::cout << taskName << " completed in "
        << std::setprecision(6) << getCurrentUTC() - executionStart << "s\n";
}

#define SKYGAZING_ASSERT_NEAR(L, R, P) assertEqualityImpl((L), (R), #L, #R, (P));
#define SKYGAZING_ASSERT_SMALL(X, P) Skygazing::assertSmallImpl((X), #X, (P));
#define SKYGAZING_TIME_EXECUTION(T) {auto t = Skygazing::getCurrentUTC(); (T); \
Skygazing::outputExecutionTimeMessage(t, #T);}
#define SKYGAZING_NAME_COMMA_VAR(X) #X, (X)

template <int decimalPlaces, typename Value>
Value round(Value w) {
    constexpr auto toPower = [](int base, int pow) {
        int result = 1;
        for (int i = 0; i < pow; ++i) {
            result *= base;
        }
        return result;
    };
    constexpr int shift = toPower(10, decimalPlaces);
    return std::round(w * shift) / shift;
}

struct RandomDataGenerator {
    explicit RandomDataGenerator(int seed) : generator(seed) {}

    std::tm getRandomDate(int year) {
        std::tm date{0, 0, 0, 1, 0, year - 1900};
        auto seconds = timegm(&date);
        constexpr int kMaxSecondsInAYear = kMaxDaysInAYear * kHoursInDay * kSecondsInDay;

        while (true) {
            time_t secondsOfTheYear =
                std::uniform_int_distribution<time_t>(0, kMaxSecondsInAYear)(generator);
            time_t newTime = seconds + secondsOfTheYear;
            std::tm newDate = *gmtime(&newTime);
            if (newDate.tm_year == date.tm_year) {
                std::exchange(date, newDate);
                break;
            }
        }
        return date;
    }

    DegreesCoordinates getRandomCoordinates(Degrees minAbsLat = 0, Degrees maxAbsLat = 90.) {
        Degrees halfSpan = maxAbsLat - minAbsLat;
        Degrees lat = std::uniform_real_distribution<Degrees>(-halfSpan, halfSpan)(generator);
        lat = lat < 0 ? lat + minAbsLat : lat - minAbsLat;
        return {lat, std::uniform_real_distribution<Degrees>(-180., 180.)(generator)};
    }

    std::mt19937 generator;
};

template <typename Value = double>
class StatisticsWatcher {
public:
    StatisticsWatcher(std::string name, std::vector<Value> &&quantiles)
            : name_(std::move(name)), quantiles_(std::move(quantiles)) {}

    void account(Value value) {
        samples_.push_back(value);
    }

    void dump() {
        std::cout << "-----------------------\n";
        std::cout << name_ << ", samples: " << samples_.size() << std::endl;
        auto count = static_cast<double>(samples_.size());
        const int width = 10;
        const std::string sep = "  ";
        const int precision = width - static_cast<int>(sep.size());
        constexpr int decimalPoints = 4;
        std::cout << std::setw(width) << "min";
        for (auto q : quantiles_) {
            std::cout << std::setw(width) << std::setprecision(width - 6) << q;
        }
        std::cout << std::setw(width) << "max\n";
        std::sort(samples_.begin(), samples_.end());
        std::cout << sep << std::setw(precision) << std::setprecision(precision - 2)
                  << round<decimalPoints>(samples_.front());
        for (auto q : quantiles_) {
            int index = static_cast<int>(q * count);
            std::cout << sep << std::setw(precision) << std::setprecision(precision - 2)
                      << round<decimalPoints>(samples_.at(index));
        }
        std::cout << sep << std::setw(precision) << round<decimalPoints>(samples_.back())
            << std::endl;
        Value squareSum = std::accumulate(
                samples_.begin(), samples_.end(), Value{}, [](Value sum, Value next)
        {
            return sum + next * next;
        });
        std::cout << "span: " << samples_.back() - samples_.front() << std::endl;
        std::cout << "rmsd: " << std::sqrt(squareSum / count) << std::endl;
    }

private:
    const std::string name_;
    std::vector<Value> samples_;
    std::vector<Value> quantiles_;
};

struct StatisticsAggregator {
    void account(const std::string &name, double value) {
        auto ins = watchers.try_emplace(name, name, std::vector<double>{0.01, 0.1, 0.5, 0.9, 0.99});
        ins.first->second.account(value);
        if (ins.second) {
            declarationOrder.push_back(name);
        }
    }

    void accountForRadsDiff(const std::string &name, double lhsRads, double rhsRads) {
        auto diff = normalizeHugeRads(std::abs(lhsRads - rhsRads));
        account(name + " diff", Skygazing::degreesFromRads(diff));
    }

    void accountForHaversineDistance(const std::string &name,
            const Coordinates &lhs, const Coordinates &rhs)
    {
        accountForRadsDiff(name , haversineDistance(lhs, rhs), 0.);
    }

    ~StatisticsAggregator() {
        for (const auto &name : declarationOrder) {
            auto &[_, watcher] = *watchers.find(name);
            watcher.dump();
            watchers.erase(name);
        }
    }

    std::unordered_map<std::string, StatisticsWatcher<double>> watchers;
    std::vector<std::string> declarationOrder;
};

} // namespace Skygazing
