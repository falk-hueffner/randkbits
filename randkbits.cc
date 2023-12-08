/* randkbits -- generating random words with exactly k 1-bits
   Copyright (C) 2015  Falk Hüffner

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along
   with this program; if not, write to the Free Software Foundation, Inc.,
   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.  */

#include <cstdint>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <limits>
#include <utility>

#include <immintrin.h>

//#undef FAST_RNG
//#define FAST_RNG 1

#ifndef WORD
typedef uint32_t word;
#else
typedef WORD word;
#endif

constexpr auto BITS_PER_WORD = std::numeric_limits<word>::digits;

inline int popcount(word x) {
    static_assert(sizeof (word) <= sizeof (int) || sizeof (word) == sizeof (long) || sizeof (word) == sizeof (long long),
		  "cannot determine popcount intrinsic");
    if (sizeof (word) <= sizeof (int))
	return __builtin_popcount(x);
    else if (sizeof (word) == sizeof (long))
	return __builtin_popcountl(x);
    else
	return __builtin_popcountll(x);
}

inline int clz(word x) {
    static_assert(sizeof (word) <= sizeof (int) || sizeof (word) == sizeof (long) || sizeof (word) == sizeof (long long),
		  "cannot determine clz intrinsic");
    if (sizeof (word) <= sizeof (int))
	return __builtin_clz(x) - (std::numeric_limits<unsigned int>::digits - std::numeric_limits<word>::digits);
    else if (sizeof (word) == sizeof (long))
	return __builtin_clzl(x);
    else
	return __builtin_clzll(x);
}

#ifndef FAST_RNG
#include <random>
word randomWord() {
    static std::default_random_engine generator;
    return std::uniform_int_distribution<word>()(generator);
}

word randomRange(word min, word max) {
    static std::default_random_engine generator;
    return std::uniform_int_distribution<word>(min, max)(generator);
}
#else
static_assert(sizeof (uint64_t) >= sizeof (word), "word too large");
word randomWord() {
    // xorshift128+ (http://arxiv.org/abs/1404.0390)
    static uint64_t s[2] = { 0xf04375ee4da57183, 0x458090208cb72372 }; // from /dev/random
    uint64_t x = s[0];
    uint64_t const y = s[1];
    s[0] = y;
    x ^= x << 23;
    x ^= x >> 17;
    x ^= y ^ (y >> 26);
    s[1] = x;
    return x + y;
}

word randomRange(word min, word max) {
    word range = max - min;
    if (range == 0)
	return min;		// clz(0) is undefined
    word mask = word(~word(0)) >> clz(range);
    word ret;
    int i = 0;
    do {
	ret = randomWord() & mask;
	++i;
    } while (ret > range);
    return min + ret;
}
#endif

// set random bits
word randomKBits0(int k) {
    bool inv = false;
    if (k > BITS_PER_WORD / 2) {
	k = BITS_PER_WORD - k;
	inv = true;
    }
    word w = 0;
    while (k) {
	int i = randomWord() % BITS_PER_WORD;
	if ((w & (word(1) << i)) == 0) {
	    w |= word(1) << i;
	    --k;
	}
    }
    return inv ? ~w : w;
}

// Fisher-Yates shuffle
word randomKBits1(int k) {
    bool inv = false;
    if (k > BITS_PER_WORD / 2) {
	k = BITS_PER_WORD - k;
	inv = true;
    }
    int bits[BITS_PER_WORD];
    for (int i = 0; i < BITS_PER_WORD; ++i)
	bits[i] = i;
    for (int i = 0; i < k; ++i) {
	int j = randomRange(i, BITS_PER_WORD - 1);
	std::swap(bits[i], bits[j]);
    }
    word w = 0;
    for (int i = 0; i < k; ++i)
	w |= word(1) << bits[i];
    return inv ? ~w : w;
}

// in-place Fisher-Yates shuffle
inline word bitswap(word w, int i, int j) {
    word x = ((w >> i) ^ (w >> j)) & 1;
    w ^= (x << i) | (x << j);
    return w;
}
word randomKBits2(int k) {
    bool inv = false;
    if (k > BITS_PER_WORD / 2) {
	k = BITS_PER_WORD - k;
	inv = true;
    }
    word w = k == 0 ? 0 : (word(1) << k) - 1;
    for (int i = 0; i < k; ++i) {
	int j = randomRange(i, BITS_PER_WORD - 1);
	w = bitswap(w, i, j);
    }
    return inv ? ~w : w;
}

// bisection method by Stack Overflow user mic006
// http://stackoverflow.com/a/28243691/378360
word randomKBits3(int k) {
    word min = 0;
    word max = word(~word(0));
    int n = 0;
    while (n != k) {
	word x = randomWord();
	x = min | (x & max);
	n = popcount(x);
	if (n > k)
	    max = x;
	else
	    min = x;
    }
    return min;
}

// enumeration method
word binCoeff[BITS_PER_WORD + 1][BITS_PER_WORD + 1];
__attribute__((constructor))
void initBinCoeff() {
    for (int n = 0; n <= BITS_PER_WORD; ++n)
        for (int k = 0; k <= n; ++k)
            if (k == 0)
                binCoeff[n][k] = 1;
            else
                binCoeff[n][k] = binCoeff[n-1][k-1] + binCoeff[n-1][k];
}
inline word ithCombination(int n, int k, word i) {
    // i is zero-based
    word x = 0;
    word b = 1;
    while (k) {
        word c = binCoeff[n - 1][k - 1];
        if (i < c) {
            x |= b;
            --k;
        } else {
            i -= c;
        }
        --n;
        b <<= 1;
    }
    return x;
}
word randomKBits4(int k) {
    bool inv = false;
    if (k > BITS_PER_WORD / 2) {
        k = BITS_PER_WORD - k;
        inv = true;
    }
    word i = randomRange(0, binCoeff[BITS_PER_WORD][k] - 1);
    word w = ithCombination(BITS_PER_WORD, k, i);
    return inv ? ~w : w;
}

// https://gist.github.com/zwegner/616aeac9a49a7e854c0743f2d7094791
// variant of above with PDEP

uint8_t choose_idx[8][8];
uint8_t choose_len[8][8];
uint8_t choose_table[1024];

// Get the next value with the same number of bits
// from https://www.chessprogramming.org/Traversing_Subsets_of_a_Set#Snoobing_the_Universe
uint64_t snoob(uint64_t x) {
   uint64_t smallest, ripple, ones;
   smallest = x & -x;
   ripple = x + smallest;
   ones = x ^ ripple;
   ones = (ones >> 2) / smallest;
   return ripple | ones;
}
void init_tables() {
    // Initialize n-choose-k tables
    int offset = 0;
    for (int i = 0; i < 8; i++) {
        int max = 1 << i;
        choose_len[i][0] = 1;
        choose_idx[i][0] = offset;
        choose_table[offset++] = 0;
        for (int j = 1; j <= i; j++) {
            int x = (1 << (j + 0)) - 1;
            int c = 0;
            choose_idx[i][j] = offset;
            while (x < max) {
                choose_table[offset++] = x;
                x = snoob(x);
                c++;
            }
            choose_len[i][j] = c;
        }
    }
}
word randomKBits5(int k) {
    word min = 0;
    word max = ~(word)0;
    int n = 0, min_n = 0, max_n = 8 * sizeof(word);
    while (max_n - min_n > 7) {
        word x = randomWord();
        x = min | (x & max);
        n = popcount(x);
        if (n > k) {
            max = x;
            max_n = n;
        }
        else {
            min = x;
            min_n = n;
        }
    }
    // Fill in extra bits
    n = max_n - min_n;
    k = k - min_n;
    int offset = randomWord() % choose_len[n][k];
    word bits = choose_table[offset + choose_idx[n][k]];
    word extra = _pdep_u64(bits, min ^ max);
    return min | extra;
}

void time(word (*f)(int)) {
    auto start = std::clock();
    for (int k = 0; k <= BITS_PER_WORD; ++k) {
	auto kstart = std::clock();
	for (int i = 0; i < 1000000; ++i)
	    volatile  __attribute__((unused)) word sink = f(k);
	auto kend = std::clock();
	if (k == 1 || k == BITS_PER_WORD / 8 || k == BITS_PER_WORD / 4 || k == BITS_PER_WORD / 2)
	    std::cout << "  k = " << k << ": "
		      << std::fixed << std::setprecision(2) << double(kend - kstart) / CLOCKS_PER_SEC;
    }
    auto end = std::clock();
    std::cout << "  total: "
	      << std::setw(5) << std::fixed << std::setprecision(2) << double(end - start) / CLOCKS_PER_SEC << std::endl;
}

int main() {
    init_tables();

    std::cout << "set random bits ";
    time(randomKBits0);
    std::cout << "shuffle         ";
    time(randomKBits1);
    std::cout << "in-place shuffle";
    time(randomKBits2);
    std::cout << "bisection       ";
    time(randomKBits3);
    std::cout << "enumeration     ";
    time(randomKBits4);
    std::cout << "bisection pdep  ";
    time(randomKBits5);
    return 0;
}
