
#ifndef AFLOAT_H
#define AFLOAT_H

#include <cstdint>
#include <cmath>
#include <limits>
#include <type_traits>
#include <cassert>

// ----------------------------
// Generic parametric aFloat
// ----------------------------
template<int SIGN_BITS, int EXP_BITS, int MANT_BITS>
struct aFloat {
    static_assert(SIGN_BITS <= 1, "Only 0 or 1 sign bit allowed");
    static_assert(SIGN_BITS + EXP_BITS + MANT_BITS <= 32,
                  "Total bits cannot exceed 32.");

    uint32_t sign : SIGN_BITS;
    uint32_t exp  : EXP_BITS;
    uint32_t mant : MANT_BITS;

    static constexpr double bias = (1 << EXP_BITS) / 2.0;
    static constexpr int totalBits = SIGN_BITS + EXP_BITS + MANT_BITS;
    using storage_t = typename std::conditional<
        (totalBits <= 8), uint8_t,
        typename std::conditional<(totalBits <= 16), uint16_t, uint32_t>::type
    >::type;

    // Decode to double
    double toDouble() const {
        double e = static_cast<double>(exp) - bias;
        double m = 1.0 + static_cast<double>(mant) / (1 << MANT_BITS);
        double value = m * std::pow(2.0, e);
        if constexpr (SIGN_BITS > 0 && sign) {
            value = -value;
        }
        return value;
    }

    // Quantize double -> nearest aFloat
    static aFloat doubleToFlt(double target) {
        aFloat best{};
        double minDiff = std::numeric_limits<double>::infinity();

        int maxSign = (SIGN_BITS > 0) ? 2 : 1;
        int maxExp  = 1 << EXP_BITS;
        int maxMant = 1 << MANT_BITS;

        for (int s = 0; s < maxSign; ++s) {
            double t = target;
            if constexpr (SIGN_BITS > 0) {
                if (s == 1) t = -target; // handle negative
            }

            for (int e = 0; e < maxExp; ++e) {
                double val_base = std::pow(2.0, e - bias);
                double ideal_m = t / val_base;
                int m_quant = static_cast<int>(std::round((ideal_m - 1.0) * maxMant));
                if (m_quant < 0) m_quant = 0;
                if (m_quant >= maxMant) m_quant = maxMant - 1;

                aFloat candidate{};
                candidate.sign = s;
                candidate.exp  = e;
                candidate.mant = m_quant;

                double val = candidate.toDouble();
                double diff = std::abs(val - target);

                if (diff < minDiff) {
                    minDiff = diff;
                    best = candidate;
                }
            }
        }

        return best;
    }

    // Quantize float -> nearest aFloat
    static aFloat floatToFlt(float target) {
        return doubleToFlt(static_cast<double>(target));
    }

    // Pack into minimal bytes
    storage_t pack() const {
        storage_t result = 0;
        result |= mant;
        result |= exp << MANT_BITS;
        if constexpr (SIGN_BITS > 0) {
            result |= sign << (MANT_BITS + EXP_BITS);
        }
        return result;
    }

    // Unpack from minimal bytes
    static aFloat unpack(storage_t data) {
        aFloat res{};
        res.mant = data & ((1 << MANT_BITS) - 1);
        res.exp  = (data >> MANT_BITS) & ((1 << EXP_BITS) - 1);
        if constexpr (SIGN_BITS > 0) {
            res.sign = (data >> (MANT_BITS + EXP_BITS)) & 1;
        }
        return res;
    }
};

// ----------------------------
// Predefined typedefs
// ----------------------------
using aFlt2  = aFloat<0, 1, 1>;
using aFlt3  = aFloat<1, 1, 1>;
using aFlt4  = aFloat<1, 2, 1>;
using aFlt5  = aFloat<1, 2, 2>;
using aFlt6  = aFloat<1, 3, 2>;
using aFlt7  = aFloat<1, 4, 2>;
using aFlt8  = aFloat<1, 4, 3>;
using aFlt9  = aFloat<1, 5, 3>;
using aFlt10 = aFloat<1, 6, 3>;
using aFlt11 = aFloat<1, 6, 4>;
using aFlt12 = aFloat<1, 7, 4>;
using aFlt13 = aFloat<1, 8, 4>;
using aFlt14 = aFloat<1, 8, 5>;
using aFlt15 = aFloat<1, 9, 5>;
using aFlt16 = aFloat<1, 10, 5>;
using aFlt17 = aFloat<1, 10, 6>;
using aFlt18 = aFloat<1, 11, 6>;
using aFlt19 = aFloat<1, 12, 6>;
using aFlt20 = aFloat<1, 12, 7>;
using aFlt21 = aFloat<1, 13, 7>;
using aFlt22 = aFloat<1, 14, 7>;
using aFlt23 = aFloat<1, 14, 8>;
using aFlt24 = aFloat<1, 15, 8>;
using aFlt25 = aFloat<1, 16, 8>;
using aFlt26 = aFloat<1, 16, 9>;
using aFlt27 = aFloat<1, 17, 9>;
using aFlt28 = aFloat<1, 18, 9>;
using aFlt29 = aFloat<1, 18, 10>;
using aFlt30 = aFloat<1, 19, 10>;
using aFlt31 = aFloat<1, 20, 10>;
using aFlt32 = aFloat<1, 20, 11>;

#endif // AFLOAT_H
