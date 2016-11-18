#pragma once

#include "Number.h"

#include "NumberMpz.h"
#include "NumberMpq.h"
#include "NumberClI.h"





namespace carl {



#ifdef USE_CLN_NUMBERS

	//TODO: where to put this?
	//NOTE: this can only be uncommented once adaption_cln/operations.h is not used anymore (defined twice)
	//static const cln::cl_RA ONE_DIVIDED_BY_10_TO_THE_POWER_OF_23 = cln::cl_RA(1)/cln::expt(cln::cl_RA(10), 23);
	//static const cln::cl_RA ONE_DIVIDED_BY_10_TO_THE_POWER_OF_52 = cln::cl_RA(1)/cln::expt(cln::cl_RA(10), 52);

	template<>
	class Number<cln::cl_RA> : public BaseNumber<cln::cl_RA,Number<cln::cl_RA>> {
	private:
	cln::cl_RA scaleByPowerOfTwo(const cln::cl_RA& a, int exp); //auxiliary function 
	public:

		Number(): BaseNumber() {}
		explicit Number(const cln::cl_RA& t): BaseNumber(t) {}
		explicit Number(cln::cl_RA&& t): BaseNumber(t) {}
		Number(const Number<cln::cl_RA>& n): BaseNumber(n) {}
		Number(Number<cln::cl_RA>&& n) noexcept : BaseNumber(n) {}
		explicit Number(int n) : BaseNumber(n) {}
		explicit Number(long long int n) { mData = cln::cl_RA(n); }
		explicit Number(unsigned long long int n) { mData = cln::cl_RA(n); }


		//The following constructors can maybe be grouped together in a Rational-superclass	
		//TODO: explicit or not?
		explicit Number(double d) { mData = cln::rationalize(d);
			/*std::ostringstream os;
			os << d;
 			std::istringstream istr(os.str()); 
			cln::cl_read_flags flags = {cln::syntax_float, cln::lsyntax_all, 10, cln::default_float_format, cln::default_float_format, false}; 
			mData = cln::cl_RA(cln::read_float(istr,flags)); */
		} 
		explicit Number(float f) { mData = cln::rationalize(f);
		/*	std::ostringstream os;
			os << f;
 			std::istringstream istr(os.str()); 
			cln::cl_read_flags flags = {cln::syntax_float, cln::lsyntax_all, 10, cln::default_float_format, cln::default_float_format, false}; 
			mData = cln::cl_RA(cln::read_float(istr,flags)); */
		 } 


		Number(const std::string& s); /* {
			cln::cl_read_flags flags = {cln::syntax_float, cln::lsyntax_all, 10, cln::default_float_format, cln::default_float_format, false}; 
			std::istringstream istr(s);
			mData = cln::cl_RA(cln::read_rational(istr,flags));
		} */

		//constructs a/b:
		//(this looks hacky.. seems to be the only really functioning way though: take the integers as strings, put the sign at the front and construct
		//cl_RA from the string "[-]a/b")
		Number(const Number<cln::cl_I>& a,const Number<cln::cl_I>& b) :
			 Number(((a.isNegative() xor b.isNegative()) ? "-" : "") + a.abs().toString()+"/"+b.abs().toString()) {}

	
		Number(const Number<cln::cl_I>& n) { mData = cln::cl_RA(n.getValue()); }
		//Number(const cln::cl_I& n) { mData = cln::cl_RA(n); }

		//TODO: Is there a better way? Maybe retrieve "pieces" that fit into the data type and add them together again
		Number(const Number<mpq_class>& n) : Number(cln::cl_RA(n.toString().c_str())) {} 
		Number(const Number<mpz_class>& n) : Number(cln::cl_RA(n.toString().c_str())) {} 


		Number<cln::cl_RA>& operator=(const Number<cln::cl_RA>& n) {
			this->mData = n.mData;
			return *this;
		}

		template<typename Other>
		Number<cln::cl_RA>& operator=(const Other& n) {
			this->mData = n;
			return *this;
		}

		Number<cln::cl_RA>& operator=(Number<cln::cl_RA>&& n) noexcept {
			this->mData = std::move(n.mData);
			return *this;
		}


		

		
		//could probably be changed such that it's the same as in all other classes
		inline bool isZero() {
			return zerop(mData);
		}


		//these 3 methods are the same as for mpq, mpz
		inline bool isOne() {
			return mData  == carl::constant_one<cln::cl_RA>().get();
		}


		inline bool isPositive() const {
			return mData > carl::constant_zero<cln::cl_RA>().get();
		}


		inline bool isNegative() const {
			return mData < carl::constant_zero<cln::cl_RA>().get();
		}

		/**
		 * Extract the numerator from a fraction.
		 * @return Numerator.
		 */
		inline Number<cln::cl_I> getNum() const {
			return Number<cln::cl_I> (cln::numerator(mData));
		}

		/**
		 * Extract the denominator from a fraction.
		 * @return Denominator.
		 */
		inline Number<cln::cl_I> getDenom() const {
			return Number<cln::cl_I> (cln::denominator(mData));
		}

		/**
		 * Check if a fraction is integral.
		 * @return true.
		 */
		inline bool isInteger() {
			return this->getDenom().isOne();
		}


		/**
		 * Get the bit size of the representation of a fraction.
		 * @param n A fraction.
		 * @return Bit size of n.
		 */
		inline std::size_t bitsize() {
			return cln::integer_length(cln::numerator(mData)) + cln::integer_length(cln::denominator(mData));
		}

		/**
		 * Converts the given fraction to a double.
		 * @param n A fraction.
		 * @return Double.
		 */
		inline double toDouble() {
			return cln::double_approx(mData);
		}



		template<typename Integer>
		inline Integer toInt();



		
		//TODO: Rationalize as constructors!!
		/*template<>
		cln::cl_RA rationalize<cln::cl_RA>(double n);

		template<>
		cln::cl_RA rationalize<cln::cl_RA>(float n);

		template<>
		inline cln::cl_RA rationalize<cln::cl_RA>(int n) {
		    return cln::cl_RA(n);
		}

		template<>
		inline cln::cl_RA rationalize<cln::cl_RA>(uint n) {
			return cln::cl_RA(n);
		}

		template<>
		inline cln::cl_RA rationalize<cln::cl_RA>(sint n) {
			return cln::cl_RA(n);
		}

		template<>
		cln::cl_RA rationalize<cln::cl_RA>(const std::string& inputstring); */



		/**
		 * Get absolute value of a fraction.
		 * @return \f$|n|\f$.
		 */
		inline Number<cln::cl_RA> abs() const {
			return Number(cln::abs(mData));
		}

		/**
		 * Round a fraction to next integer.
		 * @return The next integer.
		 */
		inline Number<cln::cl_I> round() {
			return Number<cln::cl_I> (cln::round1(mData));
		}



		/**
		 * Round down a fraction.
		 * @return \f$\lfloor n \rfloor\f$.
		 */
		inline Number<cln::cl_I> floor() {
			return Number<cln::cl_I> (cln::floor1(mData));
		}



		/**
		 * Round up a fraction.
		 * @return \f$\lceil n \rceil\f$.
		 */
		inline Number<cln::cl_I> ceil() {
			return Number<cln::cl_I> (cln::ceiling1(mData));
		}





		/**
		 * Calculate the greatest common divisor of two fractions.
		 * Asserts that the arguments are integral.
		 * @param b Second argument.
		 * @return Gcd of a and b.
		 */
		inline Number<cln::cl_RA> gcd(const Number<cln::cl_RA>& b) {
			return Number(cln::gcd(cln::numerator(mData),cln::numerator(b.mData)) / cln::lcm(cln::denominator(mData),cln::denominator(b.mData)));
		}


		/**
		 * Calculate the least common multiple of two fractions.
		 * Asserts that the arguments are integral.
		 * @param b Second argument.
		 * @return Lcm of a and b.
		 */
		inline Number<cln::cl_RA> lcm(const Number<cln::cl_RA>& b) {
		    assert( this->isInteger());
		    assert( b.isInteger() );
			return Number(cln::lcm(cln::numerator(mData),cln::numerator(b.mData)));
		}

		/**
		 * Calculate the power of some fraction to some positive integer.
		 * @param e Exponent.
		 * @return \f$n^e\f$
		 */
		inline Number<cln::cl_RA> pow(std::size_t e) {
			return Number(cln::expt(mData, int(e)));
		}

		inline Number<cln::cl_RA> log() {
			return Number(cln::rationalize(cln::realpart(cln::log(mData))));
		}

		//Note that with std::cos, std::sin the result is more precise than with cln::sin, cln::cos!!
		inline Number<cln::cl_RA> sin() {
			return Number(std::sin(toDouble()));
		}

		inline Number<cln::cl_RA> cos() {
			return Number(std::cos(toDouble()));
		}

		/**
		 * Calculate the square root of a fraction if possible.
		 *
		 * @param b A reference to the rational, in which the result is stored.
		 * @return true, if the number to calculate the square root for is a square;
		 *         false, otherwise.
		 */
		bool sqrt_exact(Number<cln::cl_RA>& b);

		Number<cln::cl_RA> sqrt();

		/**
		 * Calculate the square root of a fraction.
		 *
		 * If we are able to find a an \f$x\f$ such that \f$x\f$ is the exact root of \f$a\f$, \f$(x,x)\f$ is returned.
		 * If we can not find such a number (note that such a number might not even exist), \f$(x,y)\f$ is returned with \f$ x < \sqrt{a} < y \f$.
		 * Note that we try to find bounds that are very close to the actual square root. If a small representation is more important than a small interval, sqrt_fast should be used.
		 * @return Interval containing the square root of a.
		 */
		std::pair<Number<cln::cl_RA>, Number<cln::cl_RA>> sqrt_safe();

		/**
		 * Compute square root in a fast but less precise way.
		 * Use cln::sqrt() to obtain an approximation. If the result is rational, i.e. the result is exact, use this result.
		 * Otherwise use the nearest integers as bounds on the square root.
		 * @param a Some number.
		 * @return [x,x] if sqrt(a) = x is rational, otherwise [y,z] for y,z integer and y < sqrt(a) < z.
		 */
		std::pair<Number<cln::cl_RA>, Number<cln::cl_RA>> sqrt_fast();



		/**
		 * Divide two fractions.
		 * @param b Second argument.
		 * @return \f$ this / b \f$.
		 */
		inline Number<cln::cl_RA> div(const Number<cln::cl_RA>& b) {
			return Number(mData / b.mData);
		}



		/**
		 * Divide two fractions.
		 * @param b Second argument.
		 * @return \f$ this / b \f$.
		 */
		inline Number<cln::cl_RA> quotient(const Number<cln::cl_RA>& b)
		{
			return Number(mData / b.mData);
		}


		inline Number<cln::cl_RA> reciprocal() {
			return Number(cln::recip(mData));
		}

		std::string toString(bool _infix=true) const;



	};

		/**
		 * Convert a fraction to an integer.
		 * This method assert, that the given fraction is an integer, i.e. that the denominator is one.
		 * @param n A fraction.
		 * @return An integer.
		 */
		template<>
		inline Number<cln::cl_I> Number<cln::cl_RA>::toInt<Number<cln::cl_I>>() {
			assert(isInteger());
			return Number<cln::cl_I>(getNum());
		}
		//TODO: is this correct?!
		template<>
		inline sint Number<cln::cl_RA>::toInt<sint>() {
			return toInt<Number<cln::cl_I>>().toInt<sint>();
		}
		template<>
		inline uint Number<cln::cl_RA>::toInt<uint>() {
		    return toInt<Number<cln::cl_I>>().toInt<uint>();
		}


#endif

}

