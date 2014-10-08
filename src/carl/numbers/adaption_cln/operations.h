/** 
 * @file   adaption_cln/operations.h
 * @ingroup cln
 * @author Gereon Kremer <gereon.kremer@cs.rwth-aachen.de>
 * @author Sebastian Junges
 * 
 * @warning This file should never be included directly but only via operations.h
 * 
 */

#pragma once

#include <cassert>
#include <limits>
#include <cmath>
#include "typetraits.h"
#include "boost/algorithm/string.hpp"

namespace carl {

/**
 * Extract the numerator from a fraction.
 * @param n Fraction.
 * @return Numerator.
 */
inline const cln::cl_I getNum(const cln::cl_RA& n) {
	return cln::numerator(n);
}

/**
 * Extract the denominator from a fraction.
 * @param n Fraction.
 * @return Denominator.
 */
inline const cln::cl_I getDenom(const cln::cl_RA& n) {
	return cln::denominator(n);
}

/**
 * Check if a number is integral.
 * As cln::cl_I are always integral, this method returns true.
 * @param An integer.
 * @return true.
 */
inline bool isInteger(const cln::cl_I&) {
	return true;
}

/**
 * Check if a fraction is integral.
 * @param n A fraction.
 * @return true.
 */
inline bool isInteger(const cln::cl_RA& n) {
	return getDenom(n) == (cln::cl_I)(1);
}

/**
 * Get the bit size of the representation of a integer.
 * @param n An integer.
 * @return Bit size of n.
 */
inline std::size_t bitsize(const cln::cl_I& n) {
	return cln::integer_length(n);
}
/**
 * Get the bit size of the representation of a fraction.
 * @param n A fraction.
 * @return Bit size of n.
 */
inline std::size_t bitsize(const cln::cl_RA& n) {
	return cln::integer_length(getNum(n)) + cln::integer_length(getDenom(n));
}

/**
 * Converts the given fraction to a double.
 * @param n A fraction.
 * @return Double.
 */
inline double toDouble(const cln::cl_RA& n) {
	return cln::double_approx(n);
}
/**
 * Converts the given integer to a double.
 * @param n An integer.
 * @return Double.
 */
inline double toDouble(const cln::cl_I& n) {
	return cln::double_approx(n);
}

template<typename Integer>
inline Integer toInt(const cln::cl_I& n);
template<typename Integer>
inline Integer toInt(const cln::cl_RA& n);

/**
 * Convert a cln integer to an int.
 * @param n An integer.
 * @return n as int.
 */
template<>
inline int toInt<int>(const cln::cl_I& n) {
    assert(n <= std::numeric_limits<int>::max());
    assert(n >= std::numeric_limits<int>::min());
    return cln::cl_I_to_int(n);
}

/**
 * Convert a cln integer to an unsigned.
 * @param n An integer.
 * @return n as unsigned.
 */
template<>
inline unsigned toInt<unsigned>(const cln::cl_I& n) {
    assert(n <= std::numeric_limits<unsigned>::max());
	assert(n >= std::numeric_limits<unsigned>::min());
    return cln::cl_I_to_uint(n);
}

/**
 * Convert a cln integer to a longint.
 * @param n An integer.
 * @return n as long int.
 */
template<>
inline long int toInt<long int>(const cln::cl_I& n) {
    assert(n <= std::numeric_limits<long int>::max());
    assert(n >= std::numeric_limits<long int>::min());
    return cln::cl_I_to_long(n);
}

/**
 * Convert a fraction to an integer.
 * This method assert, that the given fraction is an integer, i.e. that the denominator is one.
 * @param n A fraction.
 * @return An integer.
 */
template<>
inline cln::cl_I toInt<cln::cl_I>(const cln::cl_RA& n) {
	assert(isInteger(n));
	return getNum(n);
}

/**
 * Convert a fraction to an unsigned.
 * @param n A fraction.
 * @return n as unsigned.
 */
template<>
inline unsigned toInt<unsigned>(const cln::cl_RA& n) {
	return toInt<unsigned>(toInt<cln::cl_I>(n));
}

/**
 * Convert a cln fraction to a cln long float.
 * @param n A fraction.
 * @return n as cln::cl_LF.
 */
inline cln::cl_LF toLF(const cln::cl_RA& n) {
	return cln::cl_R_to_LF(n, std::max(cln::integer_length(cln::numerator(n)), cln::integer_length(cln::denominator(n))));
}

template<typename T>
inline T rationalize(double n);
template<typename T>
inline T rationalize(float n);

static const cln::cl_RA ONE_DIVIDED_BY_10_TO_THE_POWER_OF_23 = cln::cl_RA(1)/cln::expt(cln::cl_RA(10), 23);
static const cln::cl_RA ONE_DIVIDED_BY_10_TO_THE_POWER_OF_52 = cln::cl_RA(1)/cln::expt(cln::cl_RA(10), 52);

template<>
inline cln::cl_RA rationalize<cln::cl_RA>(double n) {
	switch (std::fpclassify(n)) {
		case FP_NORMAL: // normalized are fully supported
			return cln::rationalize(cln::cl_R(n));
		case FP_SUBNORMAL: { // subnormals result in underflows, hence the value of the double is 0.f, where f is the significand precision
				assert(sizeof(n) == 8);
				long int significandBits = reinterpret_cast<long int>(&n);
				significandBits = (significandBits << 12) >> 12;
				if( n < 0 )
					significandBits = -significandBits;
				return cln::cl_RA( significandBits ) * ONE_DIVIDED_BY_10_TO_THE_POWER_OF_52;
            }
		case FP_ZERO:
			return 0;
		case FP_NAN: // NaN and infinite are not supported
		case FP_INFINITE:
			assert(false);
			break;
	}
	return 0;
}

template<>
inline cln::cl_RA rationalize<cln::cl_RA>(float n) {
	switch (std::fpclassify(n)) {
		case FP_NORMAL: // normalized are fully supported
			return cln::rationalize(cln::cl_R(n));
		case FP_SUBNORMAL: { // subnormals result in underflows, hence the value of the double is 0.f, where f is the significand precision
				assert(sizeof(n) == 4);
				long int significandBits = reinterpret_cast<long int>(&n);
				significandBits = (significandBits << 9) >> 9;
				if( n < 0 )
					significandBits = -significandBits;
				return cln::cl_RA( significandBits ) * ONE_DIVIDED_BY_10_TO_THE_POWER_OF_23;
            }
		case FP_ZERO:
			return 0;
		case FP_NAN: // NaN and infinite are not supported
		case FP_INFINITE:
			assert(false);
			break;
	}
	return 0;
}

/**
 * Get absolute value of an integer.
 * @param n An integer.
 * @return \f$|n|\f$.
 */
inline cln::cl_I abs(const cln::cl_I& n) {
	return cln::abs(n);
}

/**
 * Get absolute value of a fraction.
 * @param n A fraction.
 * @return \f$|n|\f$.
 */
inline cln::cl_RA abs(const cln::cl_RA& n) {
	return cln::abs(n);
}

/**
 * Round down a fraction.
 * @param n A fraction.
 * @return \f$\lfloor n \rfloor\f$.
 */
inline cln::cl_I floor(const cln::cl_RA& n) {
	return cln::floor1(n);
}

/**
 * Round down an integer.
 * @param n An integer.
 * @return \f$\lfloor n \rfloor\f$.
 */
inline cln::cl_I floor(const cln::cl_I& n) {
	return n;
}

/**
 * Round up a fraction.
 * @param n A fraction.
 * @return \f$\lceil n \rceil\f$.
 */
inline cln::cl_I ceil(const cln::cl_RA& n) {
	return cln::ceiling1(n);
}

/**
 * Round up an integer.
 * @param n An integer.
 * @return \f$\lceil n \rceil\f$.
 */
inline cln::cl_I ceil(const cln::cl_I& n) {
	return n;
}

/**
 * Calculate the greatest common divisor of two integers.
 * @param a First argument.
 * @param b Second argument.
 * @return Gcd of a and b.
 */
inline cln::cl_I gcd(const cln::cl_I& a, const cln::cl_I& b) {
	return cln::gcd(a,b);
}

/**
 * Calculate the greatest common divisor of two integers.
 * Stores the result in the first argument.
 * @param a First argument.
 * @param b Second argument.
 * @return Updated a.
 */
inline cln::cl_I& gcd_assign(cln::cl_I& a, const cln::cl_I& b) {
    a = cln::gcd(a,b);
	return a;
}

/**
 * Calculate the greatest common divisor of two fractions.
 * Stores the result in the first argument.
 * Asserts that the arguments are integral.
 * @param a First argument.
 * @param b Second argument.
 * @return Updated a.
 */
inline cln::cl_RA& gcd_assign(cln::cl_RA& a, const cln::cl_RA& b) {
    assert( carl::isInteger( a ) );
    assert( carl::isInteger( b ) );
	a = cln::gcd(carl::getNum(a),carl::getNum(b));
	return a;
}

/**
 * Calculate the greatest common divisor of two fractions.
 * Asserts that the arguments are integral.
 * @param a First argument.
 * @param b Second argument.
 * @return Gcd of a and b.
 */
inline cln::cl_RA gcd(const cln::cl_RA& a, const cln::cl_RA& b) {
    assert( carl::isInteger( a ) );
    assert( carl::isInteger( b ) );
	return cln::gcd(carl::getNum(a),carl::getNum(b));
}

/**
 * Calculate the least common multiple of two integers.
 * @param a First argument.
 * @param b Second argument.
 * @return Lcm of a and b.
 */
inline cln::cl_I lcm(const cln::cl_I& a, const cln::cl_I& b) {
	return cln::lcm(a,b);
}

/**
 * Calculate the least common multiple of two fractions.
 * Asserts that the arguments are integral.
 * @param a First argument.
 * @param b Second argument.
 * @return Lcm of a and b.
 */
inline cln::cl_RA lcm(const cln::cl_RA& a, const cln::cl_RA& b) {
    assert( carl::isInteger( a ) );
    assert( carl::isInteger( b ) );
	return cln::lcm(carl::getNum(a),carl::getNum(b));
}

/**
 * Calculate the power of some fraction to some positive integer.
 * @param n Basis.
 * @param e Exponent.
 * @return \f$n^e\f$
 */
inline cln::cl_RA pow(const cln::cl_RA& n, unsigned e) {
	return cln::expt(n, (int)e);
}

/**
 * Calculate the square root of a fraction.
 * 
 * If we are able to find a an \f$x\f$ such that \f$x\f$ is the exact root of \f$a\f$, \f$(x,x)\f$ is returned.
 * If we can not find such a number (note that such a number might not even exist), \f$(x,y)\f$ is returned with \f$ x < \sqrt{a} < y \f$.
 * Note that we try to find bounds that are very close to the actual square root. If a small representation is more important than a small interval, sqrt_fast should be used.
 * @param a A fraction.
 * @return Interval containing the square root of a.
 */
inline std::pair<cln::cl_RA, cln::cl_RA> sqrt(const cln::cl_RA& a) {
    assert( a >= 0 );
    cln::cl_R root = cln::sqrt(toLF(a));
    cln::cl_RA rroot = cln::rationalize(root);
    if( rroot == root ) // the result of the sqrt operation is a rational and thus an exact solution -> return a point-Interval
    {
        return std::make_pair(rroot, rroot);
    }
    else // we need to find the second bound of the overapprox. - the first is given by the rationalized result.
    {
		// Check if root^2 > a
        if( cln::expt_pos(rroot,2) > a ) // we need to find the lower bound
        {
            cln::cl_R lower = cln::sqrt(toLF(a-rroot));
            cln::cl_RA rlower = cln::rationalize(lower);
            if( rlower == lower )
            {
                return std::make_pair(rlower, rroot);
            }
            else
            {
                cln::cl_I num = cln::numerator(rlower);
                cln::cl_I den = cln::denominator(rlower);
                --num;
                return std::make_pair( num/den, rroot );
            }
        }
        else // we need to find the upper bound
        {
            cln::cl_R upper = cln::sqrt(toLF(a+rroot));
            cln::cl_RA rupper = cln::rationalize(upper);
            if( rupper == upper )
            {
                return std::make_pair(rroot, rupper);
            }
            else
            {
                cln::cl_I num = cln::numerator(rupper);
                cln::cl_I den = cln::denominator(rupper);
                ++num;
                return std::make_pair(rroot, num/den );
            }
        }
    }
}

/**
 * Compute square root in a fast but less precise way.
 * Use cln::sqrt() to obtain an approximation. If the result is rational, i.e. the result is exact, use this result.
 * Otherwise use the nearest integers as bounds on the square root.
 * @param a Some number.
 * @return [x,x] if sqrt(a) = x is rational, otherwise [y,z] for y,z integer and y < sqrt(a) < z. 
 */
inline std::pair<cln::cl_RA, cln::cl_RA> sqrt_fast(const cln::cl_RA& a) {
	assert(a >= 0);
	cln::cl_R tmp = cln::sqrt(toLF(a));
	cln::cl_RA root = cln::rationalize(tmp);
	if(root == tmp) {
		// root is a cln::cl_RA
		return std::make_pair(root, root);
	} else {
		// root is a cln::cl_LF. In this case, it is not integer and we can assume that the surrounding integers contain the actual root.
		cln::cl_I lower = carl::floor(root);
		cln::cl_I upper = carl::ceil(root);
		assert(cln::expt_pos(lower,2) < a);
		assert(cln::expt_pos(upper,2) > a);
		return std::make_pair(lower, upper);
    }
}

/**
 * Calculate the remainder of the integer division.
 * @param a First argument.
 * @param b Second argument.
 * @return \f$a \% b\f$.
 */
inline cln::cl_I mod(const cln::cl_I& a, const cln::cl_I& b) {
	return cln::rem(a, b);
}

/**
 * Divide two fractions.
 * @param a First argument.
 * @param b Second argument.
 * @return \f$ a / b \f$.
 */
inline cln::cl_RA div(const cln::cl_RA& a, const cln::cl_RA& b) {
	return (a / b);
}

/**
 * Divide two integers.
 * Asserts that the remainder is zero.
 * @param a First argument.
 * @param b Second argument.
 * @return \f$ a / b \f$.
 */
inline cln::cl_I div(const cln::cl_I& a, const cln::cl_I& b) {
	assert(cln::mod(a, b) == 0);
	return cln::exquo(a, b);
}

/**
 * Divide two fractions.
 * Stores the result in the first argument.
 * @param a First argument.
 * @param b Second argument.
 * @return \f$ a / b \f$.
 */
inline cln::cl_RA& div_assign(cln::cl_RA& a, const cln::cl_RA& b) {
    a /= b;
	return a;
}

/**
 * Divide two integers.
 * Asserts that the remainder is zero.
 * Stores the result in the first argument.
 * @param a First argument.
 * @param b Second argument.
 * @return \f$ a / b \f$.
 */
inline cln::cl_I& div_assign(cln::cl_I& a, const cln::cl_I& b) {
	assert(cln::mod(a,b) == 0);
	a = cln::exquo(a, b);
    return a;
}

/**
 * Divide two fractions.
 * @param a First argument.
 * @param b Second argument.
 * @return \f$ a / b \f$.
 */
inline cln::cl_RA quotient(const cln::cl_RA& a, const cln::cl_RA& b)
{
	return a / b;
}
/**
 * Divide two integers.
 * Discards the remainder of the division.
 * @param a First argument.
 * @param b Second argument.
 * @return \f$ a / b \f$.
 */
inline cln::cl_I quotient(const cln::cl_I& a, const cln::cl_I& b)
{
	return cln::exquo(a - cln::rem(a, b), b);
}


/**
 * Calculate the remainder of the integer division.
 * @param a First argument.
 * @param b Second argument.
 * @return \f$a \% b\f$.
 */
inline cln::cl_I remainder(const cln::cl_I& a, const cln::cl_I& b) {
	return cln::rem(a, b);
}


/**
 * Divide two integers.
 * Discards the remainder of the division.
 * @param a First argument.
 * @param b Second argument.
 * @return \f$ a / b \f$.
 */
inline cln::cl_I operator/(const cln::cl_I& a, const cln::cl_I& b)
{
	return quotient(a,b);
}

template<typename T>
inline T rationalize(const std::string& n);

template<>
inline cln::cl_RA rationalize<cln::cl_RA>(const std::string& inputstring) {
	std::vector<std::string> strs;
    boost::split(strs, inputstring, boost::is_any_of("."));

    if(strs.size() > 2)
    {
        throw std::invalid_argument("More than one delimiter in the string.");
    }
    cln::cl_RA result(0);
    if(!strs.front().empty())
    {
        result += cln::cl_RA(strs.front().c_str());
    }
    if(strs.size() > 1)
    {
        result += (cln::cl_RA(strs.back().c_str())/carl::pow(cln::cl_I(10),static_cast<unsigned>(strs.back().size())));
    }
    return result;
}

}
