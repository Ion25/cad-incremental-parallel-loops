#pragma once


#include <carl-common/util/variant_util.h>
#include <carl-formula/bitvector/BVValue.h>
#include <carl-arith/extended/MultivariateRoot.h>
#include <carl-arith/ran/ran.h>
#include <carl-arith/vs/SqrtEx.h>
#include <carl-formula/uninterpreted/SortValue.h>
#include <carl-formula/uninterpreted/UFModel.h>

#include <memory>
#include <variant>

namespace carl
{
	template<typename Rational, typename Poly>
	class ModelValue;
	template<typename Rational, typename Poly>
	class ModelSubstitution;
	template<typename Rational, typename Poly>
	class ModelMVRootSubstitution;
	template<typename Rational, typename Poly>
	using ModelSubstitutionPtr = std::unique_ptr<ModelSubstitution<Rational,Poly>>;
	template<typename Rational, typename Poly, typename Substitution, typename... Args>
	inline ModelValue<Rational,Poly> createSubstitution(Args&&... args);
	template<typename Rational, typename Poly, typename Substitution, typename... Args>
	inline ModelSubstitutionPtr<Rational,Poly> createSubstitutionPtr(Args&&... args);
	template<typename Rational, typename Poly>
	inline ModelValue<Rational,Poly> createSubstitution(const MultivariateRoot<Poly>& mr);

	/**
	 * This class represents infinity or minus infinity, depending on its flag positive.
	 * The default is minus infinity.
	 */
	struct InfinityValue { bool positive = false; };

	inline bool operator==(InfinityValue lhs, InfinityValue rhs) {
		return lhs.positive == rhs.positive;
	}
	
	inline std::ostream& operator<<(std::ostream& os, const InfinityValue& iv) {
		return os << (iv.positive ? "+" : "-") << "oo";
	}

    template<typename Poly>
    struct SymbolicInterval {
        std::optional<MultivariateRoot<Poly>> lower = std::nullopt;
        bool lower_strict = true;
        std::optional<MultivariateRoot<Poly>> upper = std::nullopt;
        bool upper_strict = true;
    };
	
    template<typename Poly>
	inline std::ostream& operator<<(std::ostream& os, const SymbolicInterval<Poly>& i) {
        os << (i.lower_strict ? "(" : "[");
        if (i.lower) os << i.lower;
        else os << "-oo";
        os << ", ";
        if (i.upper) os << i.upper;
        else os << "oo";
        os << (i.upper_strict ? ")" : "]");
		return os;
	}

    template<typename Poly>
    inline bool operator==(const SymbolicInterval<Poly>& lhs, const SymbolicInterval<Poly>& rhs) {
		return (
               lhs.lower_strict == rhs.lower_strict
            && lhs.lower == rhs.lower
            && lhs.upper_strict == rhs.upper_strict
            && lhs.upper == rhs.upper
        );
	}


    template<typename RAN>
    struct Infinitesimal {
        RAN value;
        bool is_infimum; // otherwise it is a supremum
    };

    template<typename RAN>
	inline std::ostream& operator<<(std::ostream& os, const Infinitesimal<RAN>& d) {
        os << "(" << (d.is_infimum ? "+ " : "- ") << d.value << " epsilon)";
		return os;
	}

    template<typename RAN>
    inline bool operator==(const Infinitesimal<RAN>& lhs, const Infinitesimal<RAN>& rhs) {
		return lhs.value == rhs.value && lhs.is_infimum == rhs.is_infimum;
	}

	
	/**
	 * Represent a sum type/variant over the different kinds of values that
	 * can be assigned to the different kinds of variables that exist in
	 * CARL and to use them in a more uniform way,
	 * e.g. a plain "bool", "infinity", a "carl::RealAlgebraicNumber",
	 * a (bitvector) "carl::BVValue" etc.
	 */
	template<typename Rational, typename Poly>
	class ModelValue {
		template<typename R, typename P>
		friend std::ostream& operator<<(std::ostream& os, const ModelValue<R,P>& mv);
		template<typename R, typename P>
		friend bool operator==(const ModelValue<R,P>& lhs, const ModelValue<R,P>& rhs);
		template<typename R, typename P>
		friend bool operator<(const ModelValue<R,P>& lhs, const ModelValue<R,P>& rhs);

        using RAN = Poly::RootType;
		/**
		 * Base type we are deriving from.
		 */
		using Super = std::variant<
                                bool,
                                Rational,
                                SqrtEx<Poly>,
                                RAN,
                                BVValue,
                                SortValue,
                                UFModel,
                                InfinityValue,
                                ModelSubstitutionPtr<Rational,Poly>,
                                SymbolicInterval<Poly>,
                                Infinitesimal<RAN>
                                >;
		
		Super mData;

		template<typename Variant>
		Super clone(const Variant& v) const {
			return std::visit(overloaded {
				[](const MultivariateRoot<Poly>& mr) {
					return Super(createSubstitutionPtr<Rational,Poly,ModelMVRootSubstitution<Rational,Poly>>(mr));
				},
				[](const ModelSubstitutionPtr<Rational,Poly>& subs) { return Super(subs->clone()); },
				[](const auto& t) { return Super(t); }
			}, v);
		}
		
	public:
		/**
		 * Default constructor.
		 */
		ModelValue() = default;

		ModelValue(const ModelValue& mv)
			: mData(clone(mv.mData))
		{}
		ModelValue(ModelValue&& mv) = default;

		/**
		 * Initialize the Assignment from some valid type of the underlying variant.
		 */
		template<typename T, typename T2 = typename std::enable_if<convertible_to_variant<T, Super>::value, T>::type>
		ModelValue(const T& _t): mData(_t) {}

		template<typename T, typename T2 = typename std::enable_if<convertible_to_variant<T, Super>::value, T>::type>
		ModelValue(T&& _t): mData(std::move(_t)) {}
		
		template<typename ...Args>
		ModelValue(const std::variant<Args...>& variant): mData(clone(variant)) {}
		
		ModelValue(const MultivariateRoot<Poly>& mr): mData(createSubstitution<Rational,Poly>(mr).asSubstitution()) {}

		ModelValue& operator=(const ModelValue& mv) {
			// Attention: Replace copy assignment with constructor + move assignment to avoid problems with self-assignment
			mData = clone(mv.mData);
			return *this;
		}
		ModelValue& operator=(ModelValue&& mv) = default;
		
		/**
		 * Assign some value to the underlying variant.
		 * @param t Some value.
		 * @return *this.
		 */
		template<typename T>
		ModelValue& operator=(const T& t) {
			// Attention: Replace copy assignment with constructor + move assignment to avoid problems with self-assignment
			mData = Super(t);
			return *this;
		}
		template<typename ...Args>
		ModelValue& operator=(const std::variant<Args...>& variant) {
			// Attention: Replace copy assignment with constructor + move assignment to avoid problems with self-assignment
			mData = clone(variant);
			return *this;
		}
		ModelValue& operator=(const MultivariateRoot<Poly>& mr) {
			mData = createSubstitution<Rational,Poly>(mr).asSubstitution();
			return *this;
		}
		
		template<typename F>
		auto visit(F&& f) const {
			return std::visit(f, mData);
		}

		bool isBool()      const { return std::holds_alternative<bool>(mData); }
		bool isRational()  const { return std::holds_alternative<Rational>(mData); }
		bool isSqrtEx()    const { return std::holds_alternative<SqrtEx<Poly>>(mData); }
		bool isRAN()       const { return std::holds_alternative<RAN>(mData); }
		bool isBVValue()   const { return std::holds_alternative<BVValue>(mData); }
		bool isSortValue() const { return std::holds_alternative<SortValue>(mData); }
		bool isUFModel()   const { return std::holds_alternative<UFModel>(mData); }
		bool isSubstitution() const { return std::holds_alternative<ModelSubstitutionPtr<Rational,Poly>>(mData); }
		bool isPlusInfinity() const {
			return std::holds_alternative<InfinityValue>(mData) && std::get<InfinityValue>(mData).positive;
		}
		bool isMinusInfinity() const {
			return std::holds_alternative<InfinityValue>(mData) && !std::get<InfinityValue>(mData).positive;
		}
        bool isSymbolicInterval() const { return std::holds_alternative<SymbolicInterval<Poly>>(mData); }
        bool isInfinitesimal() const { return std::holds_alternative<Infinitesimal<RAN>>(mData); }

		bool asBool() const { assert(isBool()); return std::get<bool>(mData); }
		const Rational&   asRational() const { assert(isRational()); return std::get<Rational>(mData); }
		const SqrtEx<Poly>& asSqrtEx() const { assert(isSqrtEx()); return std::get<SqrtEx<Poly>>(mData); }
		const RAN&             asRAN() const { assert(isRAN()); return std::get<RAN>(mData); }
		const BVValue&     asBVValue() const { assert(isBVValue()); return std::get<BVValue>(mData); }
		const SortValue& asSortValue() const { assert(isSortValue()); return std::get<SortValue>(mData); }
		const UFModel&     asUFModel() const { assert(isUFModel()); return std::get<UFModel>(mData); }
		      UFModel&     asUFModel()       { assert(isUFModel()); return std::get<UFModel>(mData); }
		const InfinityValue& asInfinity() const { assert(isPlusInfinity() || isMinusInfinity()); return std::get<InfinityValue>(mData); }
		const ModelSubstitutionPtr<Rational,Poly>& asSubstitution() const { assert(isSubstitution()); return std::get<ModelSubstitutionPtr<Rational,Poly>>(mData); }
		      ModelSubstitutionPtr<Rational,Poly>& asSubstitution()       { assert(isSubstitution()); return std::get<ModelSubstitutionPtr<Rational,Poly>>(mData); }
        const SymbolicInterval<Poly>& asSymbolicInterval() const { assert(isSymbolicInterval()); return std::get<SymbolicInterval<Poly>>(mData); }
        const Infinitesimal<RAN>& asInfinitesimal() const { assert(isInfinitesimal()); return std::get<Infinitesimal<RAN>>(mData); }
	};

	/**
	 * Check if two Assignments are equal.
	 * Two Assignments are considered equal, if both are either bool or not bool and their value is the same.
	 * 
	 * If both Assignments are not bools, the check may return false although they represent the same value.
	 * If both are numbers in different representations, this comparison is only done as a "best effort".
	 * 
	 * @param lhs First Assignment.
	 * @param rhs Second Assignment.
	 * @return lhs == rhs.
	 */
	template<typename Rational, typename Poly>
	bool operator==(const ModelValue<Rational,Poly>& lhs, const ModelValue<Rational,Poly>& rhs) {
		return lhs.mData == rhs.mData;
	}

	template<typename Rational, typename Poly>
	bool operator<(const ModelValue<Rational,Poly>& lhs, const ModelValue<Rational,Poly>& rhs) {
		return lhs.mData < rhs.mData;
	}
	
	template<typename R, typename P>
	inline std::ostream& operator<<(std::ostream& os, const ModelValue<R,P>& mv) {
        if (mv.isBool())      { return os << mv.asBool(); }
		if (mv.isRational())  { return os << mv.asRational(); }
		if (mv.isSqrtEx())    { return os << mv.asSqrtEx(); }
		if (mv.isRAN())       { return os << mv.asRAN(); }
		if (mv.isBVValue())   { return os << mv.asBVValue(); }
		if (mv.isSortValue()) { return os << mv.asSortValue(); }
		if (mv.isUFModel())   { return os << mv.asUFModel(); }
		if (mv.isSubstitution()) { return os << mv.asSubstitution(); }
		if (mv.isPlusInfinity() || mv.isMinusInfinity()) { return os << mv.asInfinity(); }
        if (mv.isSymbolicInterval()) { return os << mv.asSymbolicInterval(); }
        if (mv.isInfinitesimal()) { return os << mv.asInfinitesimal(); }
		return os << mv.mData;
	}
}
