#include "LPPolynomial.h"
#include "helper.h"
#include <poly/variable_list.h>

#include <carl-common/config.h>
#ifdef USE_LIBPOLY

namespace carl {

LPPolynomial::LPPolynomial(const LPPolynomial& rhs)
    : m_internal(lp_polynomial_new_copy(rhs.m_internal)), m_context(rhs.m_context) {
    assert(lp_polynomial_check_order(get_internal()));
}

LPPolynomial::LPPolynomial(LPPolynomial&& rhs)
    : m_internal(rhs.m_internal), m_context(std::move(rhs.m_context)) {
    rhs.m_internal = 0;
    assert(lp_polynomial_check_order(get_internal()));
}

LPPolynomial& LPPolynomial::operator=(const LPPolynomial& rhs) {
    m_internal = lp_polynomial_new_copy(rhs.m_internal);
    m_context = rhs.m_context;
    assert(lp_polynomial_check_order(get_internal()));
    return *this;
}

LPPolynomial& LPPolynomial::operator=(LPPolynomial&& rhs) {
    if (m_internal) lp_polynomial_delete(m_internal);
    m_internal = rhs.m_internal;
    m_context = std::move(rhs.m_context);
    rhs.m_internal = 0;
    assert(lp_polynomial_check_order(get_internal()));
    return *this;
}

LPPolynomial::LPPolynomial(const LPContext& context)
    : m_internal(lp_polynomial_new(context.lp_context())), m_context(context) {
    //lp_polynomial_set_external(get_internal());
    assert(lp_polynomial_check_order(get_internal()));
}

LPPolynomial::LPPolynomial(lp_polynomial_t* p, const LPContext& context)
    : m_internal(p), m_context(context) {
    //lp_polynomial_set_external(get_internal());
    assert(lp_polynomial_check_order(get_internal()));
    assert(context.lp_context() == lp_polynomial_get_context(get_internal()));
}

LPPolynomial::LPPolynomial(const LPContext& context, long val)
    : m_internal(lp_polynomial_alloc()), m_context(context) {
    lp_polynomial_construct_simple(get_internal(), context.lp_context(), mpz_class(val).get_mpz_t(), 0, 0);
    //lp_polynomial_set_external(get_internal());
    assert(lp_polynomial_check_order(get_internal()));
}

LPPolynomial::LPPolynomial(const LPContext& context, const mpz_class& val)
    : m_internal(lp_polynomial_alloc()), m_context(context) {
    lp_polynomial_construct_simple(get_internal(), context.lp_context(), val.get_mpz_t(), lp_variable_null, 0) ;
    //lp_polynomial_set_external(get_internal());
    assert(lp_polynomial_check_order(get_internal()));
}

LPPolynomial::LPPolynomial(const LPContext& context, const mpq_class& val) : LPPolynomial(context, carl::get_num(val)) {}

LPPolynomial::LPPolynomial(const LPContext& context, const Variable& var, const mpz_class& coeff, unsigned int degree)
    : m_internal(lp_polynomial_alloc()), m_context(context) {
    lp_polynomial_construct_simple(get_internal(), context.lp_context(), mpz_class(coeff).get_mpz_t(), context.lp_variable(var), degree);
    //lp_polynomial_set_external(get_internal());
    assert(lp_polynomial_check_order(get_internal()));
}

LPPolynomial::LPPolynomial(const LPContext& context, const Variable& var)
    : m_internal(lp_polynomial_alloc()), m_context(context) {
    lp_polynomial_construct_simple(get_internal(), context.lp_context(), mpz_class(1).get_mpz_t(), context.lp_variable(var), 1);
    //lp_polynomial_set_external(get_internal());
    assert(lp_polynomial_check_order(get_internal()));
}

LPPolynomial::LPPolynomial(const LPContext& context, const Variable& mainVar, const std::initializer_list<mpz_class>& coefficients)
    : LPPolynomial(context) {

    auto var = context.lp_variable(mainVar);
    auto pow = coefficients.size();

    for (const mpz_class& coeff : coefficients) {
        pow--;
        if (is_zero(coeff)) continue;
        lp_monomial_t t;
		lp_monomial_construct(context.lp_context(), &t);
		lp_monomial_set_coefficient(context.lp_context(), &t, mpz_class(coeff).get_mpz_t());
        if (pow>0)
            lp_monomial_push(&t, var, (unsigned int)pow);
        lp_polynomial_add_monomial(get_internal(), &t);
		lp_monomial_destruct(&t);
    }
    //lp_polynomial_set_external(get_internal());
}

LPPolynomial::LPPolynomial(const LPContext& context, const Variable& mainVar, const std::vector<mpz_class>& coefficients)
    : LPPolynomial(context) {

    auto var = context.lp_variable(mainVar);
    auto pow = coefficients.size();

    for (const mpz_class& coeff : coefficients) {
        pow--;
        if (is_zero(coeff)) continue;
        lp_monomial_t t;
		lp_monomial_construct(context.lp_context(), &t);
		lp_monomial_set_coefficient(context.lp_context(), &t, mpz_class(coeff).get_mpz_t());
        if (pow>0)
            lp_monomial_push(&t, var, (unsigned int)pow);
        lp_polynomial_add_monomial(get_internal(), &t);
		lp_monomial_destruct(&t);
    }
    //lp_polynomial_set_external(get_internal());
}

LPPolynomial::LPPolynomial(const LPContext& context, const Variable& mainVar, std::vector<mpz_class>&& coefficients)
    : LPPolynomial(context) {

    auto var = context.lp_variable(mainVar);
    auto pow = coefficients.size();

    for (const mpz_class& coeff : coefficients) {
        pow--;
        if (is_zero(coeff)) continue;
        lp_monomial_t t;
		lp_monomial_construct(context.lp_context(), &t);
		lp_monomial_set_coefficient(context.lp_context(), &t, mpz_class(coeff).get_mpz_t());
        if (pow>0)
            lp_monomial_push(&t, var, (unsigned int)pow);
        lp_polynomial_add_monomial(get_internal(), &t);
		lp_monomial_destruct(&t);
    }
    //lp_polynomial_set_external(get_internal());
}

LPPolynomial::LPPolynomial(const LPContext& context, const Variable& mainVar, const std::map<unsigned int, mpz_class>& coefficients)
    : LPPolynomial(context) {

    auto var = context.lp_variable(mainVar);

    for (const auto& coeff : coefficients) {
        if (is_zero(coeff.second)) continue;
        lp_monomial_t t;
		lp_monomial_construct(context.lp_context(), &t);
		lp_monomial_set_coefficient(context.lp_context(), &t, mpz_class(coeff.second).get_mpz_t());
        if (coeff.first>0)
            lp_monomial_push(&t, var, (unsigned int)coeff.first);
        lp_polynomial_add_monomial(get_internal(), &t);
		lp_monomial_destruct(&t);
    }
    //lp_polynomial_set_external(get_internal());
}

LPPolynomial::~LPPolynomial() {
    // std::cout << ">>>>>>>>>>> K : " << m_context.lp_context()->K << "\n";
    if (m_internal) lp_polynomial_delete(m_internal);
}

bool LPPolynomial::has(const Variable& var) const {
    lp_variable_list_t varList;
    lp_variable_list_construct(&varList);
    lp_polynomial_get_variables(get_internal(), &varList);
    auto lp_variable = context().lp_variable_opt(var);
    if (!lp_variable) return false;
    bool contains = lp_variable_list_contains(&varList, *lp_variable);
    lp_variable_list_destruct(&varList);
    return contains;
}

bool operator==(const LPPolynomial& lhs, const LPPolynomial& rhs) {
    return lp_polynomial_eq(lhs.get_internal(), rhs.get_internal());
}
bool operator==(const LPPolynomial& lhs, const mpz_class& rhs) {
    if (!is_number(lhs)) {
        return false;
    }
    return lhs.constant_part() == rhs;
}
bool operator==(const mpz_class& lhs, const LPPolynomial& rhs) {
    return rhs == lhs;
}

bool operator!=(const LPPolynomial& lhs, const LPPolynomial& rhs) {
    return !(lhs == rhs);
}
bool operator!=(const LPPolynomial& lhs, const mpz_class& rhs) {
    return !(lhs == rhs);
}
bool operator!=(const mpz_class& lhs, const LPPolynomial& rhs) {
    return !(lhs == rhs);
}

inline auto cmp_util(const LPPolynomial& lhs, const mpz_class& rhs) {
    lp_polynomial_t* tmp = poly_helper::construct_lp_poly(lp_polynomial_get_context(lhs.get_internal()), rhs);
    auto res = lp_polynomial_cmp(lhs.get_internal(), tmp);
    lp_polynomial_delete(tmp);
    return res;
}

bool operator<(const LPPolynomial& lhs, const LPPolynomial& rhs) {
    return lp_polynomial_cmp(lhs.get_internal(), rhs.get_internal()) < 0;
}
bool operator<(const LPPolynomial& lhs, const mpz_class& rhs) {
    return cmp_util(lhs,rhs) < 0;
}
bool operator<(const mpz_class& lhs, const LPPolynomial& rhs) {
    return cmp_util(rhs,lhs) > 0;
}

bool operator<=(const LPPolynomial& lhs, const LPPolynomial& rhs) {
    return lp_polynomial_cmp(lhs.get_internal(), rhs.get_internal()) <= 0;
}
bool operator<=(const LPPolynomial& lhs, const mpz_class& rhs) {
    return cmp_util(lhs,rhs) <= 0;
}
bool operator<=(const mpz_class& lhs, const LPPolynomial& rhs) {
    return cmp_util(rhs,lhs) >= 0;
}

bool operator>(const LPPolynomial& lhs, const LPPolynomial& rhs) {
    return lp_polynomial_cmp(lhs.get_internal(), rhs.get_internal()) > 0;
}
bool operator>(const LPPolynomial& lhs, const mpz_class& rhs) {
    return cmp_util(lhs,rhs) > 0;
}
bool operator>(const mpz_class& lhs, const LPPolynomial& rhs) {
    return cmp_util(rhs,lhs) < 0;
}

bool operator>=(const LPPolynomial& lhs, const LPPolynomial& rhs) {
    return lp_polynomial_cmp(lhs.get_internal(), rhs.get_internal()) >= 0;
}
bool operator>=(const LPPolynomial& lhs, const mpz_class& rhs) {
    return cmp_util(lhs,rhs) >= 0;
}
bool operator>=(const mpz_class& lhs, const LPPolynomial& rhs) {
    return cmp_util(rhs,lhs) <= 0;
}

LPPolynomial operator+(const LPPolynomial& lhs, const LPPolynomial& rhs) {
    assert(rhs.context() == lhs.context());
    assert(lp_polynomial_context_equal(lp_polynomial_get_context(lhs.get_internal()), lp_polynomial_get_context(rhs.get_internal())));
    LPPolynomial result(lhs.context());
    lp_polynomial_add(result.get_internal(), lhs.get_internal(), rhs.get_internal());
    return result;
}
LPPolynomial operator+(const LPPolynomial& lhs, const mpz_class& rhs) {
    return lhs + LPPolynomial(lhs.context(), rhs);
}
LPPolynomial operator+(const mpz_class& lhs, const LPPolynomial& rhs) {
    return rhs + lhs;
}

LPPolynomial operator-(const LPPolynomial& lhs, const LPPolynomial& rhs) {
    assert(rhs.context() == lhs.context());
    assert(lp_polynomial_context_equal(lp_polynomial_get_context(lhs.get_internal()), lp_polynomial_get_context(rhs.get_internal())));
    LPPolynomial result(lhs.context());
    lp_polynomial_sub(result.get_internal(), lhs.get_internal(), rhs.get_internal());
    return result;
}
LPPolynomial operator-(const LPPolynomial& lhs, const mpz_class& rhs) {
    return lhs - LPPolynomial(lhs.context(), rhs);
}
LPPolynomial operator-(const mpz_class& lhs, const LPPolynomial& rhs) {
    return LPPolynomial(rhs.context(), lhs) - rhs;
}

LPPolynomial operator*(const LPPolynomial& lhs, const LPPolynomial& rhs) {
    assert(lhs.context() == rhs.context());
    assert(lp_polynomial_context_equal(lp_polynomial_get_context(lhs.get_internal()), lp_polynomial_get_context(rhs.get_internal())));
    LPPolynomial result(lhs.context());
    lp_polynomial_mul(result.get_internal(), lhs.get_internal(), rhs.get_internal());
    return result;
}
LPPolynomial operator*(const LPPolynomial& lhs, const mpz_class& rhs) {
    return lhs * LPPolynomial(lhs.context(), rhs);
}
LPPolynomial operator*(const mpz_class& lhs, const LPPolynomial& rhs) {
    return rhs * lhs;
}

LPPolynomial& operator+=(LPPolynomial& lhs, const LPPolynomial& rhs) {
    assert(rhs.context() == lhs.context());
    assert(lp_polynomial_context_equal(lp_polynomial_get_context(lhs.get_internal()), lp_polynomial_get_context(rhs.get_internal())));
    lp_polynomial_add(lhs.get_internal(), lhs.get_internal(), rhs.get_internal());
    return lhs;
}
LPPolynomial& operator+=(LPPolynomial& lhs, const mpz_class& rhs) {
    lp_polynomial_t* tmp = poly_helper::construct_lp_poly(lp_polynomial_get_context(lhs.get_internal()), rhs);
    lp_polynomial_add(lhs.get_internal(), lhs.get_internal(), tmp);
    lp_polynomial_delete(tmp);
    return lhs;
}

LPPolynomial& operator-=(LPPolynomial& lhs, const LPPolynomial& rhs) {
    assert(rhs.context() == lhs.context());
    assert(lp_polynomial_context_equal(lp_polynomial_get_context(lhs.get_internal()), lp_polynomial_get_context(rhs.get_internal())));
    lp_polynomial_sub(lhs.get_internal(), lhs.get_internal(), rhs.get_internal());
    return lhs;
}
LPPolynomial& operator-=(LPPolynomial& lhs, const mpz_class& rhs) {
    lp_polynomial_t* tmp = poly_helper::construct_lp_poly(lp_polynomial_get_context(lhs.get_internal()), rhs);
    lp_polynomial_sub(lhs.get_internal(), lhs.get_internal(), tmp);
    lp_polynomial_delete(tmp);
    return lhs;
}

LPPolynomial& operator*=(LPPolynomial& lhs, const LPPolynomial& rhs) {
    assert(rhs.context() == lhs.context());
    assert(lp_polynomial_context_equal(lp_polynomial_get_context(lhs.get_internal()), lp_polynomial_get_context(rhs.get_internal())));
    lp_polynomial_mul(lhs.get_internal(), lhs.get_internal(), rhs.get_internal());
    return lhs;
}
LPPolynomial& operator*=(LPPolynomial& lhs, const mpz_class& rhs) {
    lp_polynomial_t* tmp = poly_helper::construct_lp_poly(lp_polynomial_get_context(lhs.get_internal()), rhs);
    lp_polynomial_mul(lhs.get_internal(), lhs.get_internal(), tmp);
    lp_polynomial_delete(tmp);
    return lhs;
}

mpz_class LPPolynomial::coprime_factor() const {
    // TODO: can this be done with content/primitive part?
    struct coprime_factor_travers {
        std::vector<mpz_class> coefficients; // coefficients of the polynomial
    };

    auto getCoeffs = [](const lp_polynomial_context_t* /*ctx*/,
                        lp_monomial_t* m,
                        void* d) {
        coprime_factor_travers& v = *static_cast<coprime_factor_travers*>(d);
        v.coefficients.push_back(*reinterpret_cast<mpz_class*>(&m->a));
    };

    // first get the coefficients of every monomial
    coprime_factor_travers travers;
    lp_polynomial_traverse(get_internal(), getCoeffs, &travers);

    if (travers.coefficients.size() == 0) {
        return 0;
    }
    mpz_class res = travers.coefficients[0];
    for (size_t i = 1; i < travers.coefficients.size(); i++) {
        res = gcd(res, travers.coefficients[i]);
    }

    return res;
}

LPPolynomial LPPolynomial::coprime_coefficients() const {
    mpz_class g = coprime_factor();
    if (g == 1) return *this;
    lp_polynomial_t* temp = poly_helper::construct_lp_poly(lp_polynomial_get_context(get_internal()), g);
    lp_polynomial_t* res = lp_polynomial_new(context().lp_context());
    lp_polynomial_div(res, get_internal(), temp);
    lp_polynomial_delete(temp);
    return LPPolynomial(res, context());
}

std::size_t LPPolynomial::total_degree() const {

    struct degree_travers {
        std::size_t degree = 0;
    };

    auto getDegree = [](const lp_polynomial_context_t* /*ctx*/,
                        lp_monomial_t* m,
                        void* d) {
        degree_travers& v = *static_cast<degree_travers*>(d);

        size_t current_degree = 0;
        // iterate over the number of variables and add up their degrees
        for (size_t i = 0; i < m->n; i++) {
            current_degree += m->p[i].d;
        }
        v.degree = std::max(v.degree, current_degree);
    };

    degree_travers travers;
    lp_polynomial_traverse(get_internal(), getDegree, &travers);

    return travers.degree;
}

std::size_t LPPolynomial::degree(Variable::Arg var) const {
    struct degree_travers {
        std::size_t degree = 0;
        lp_variable_t var; // the variable we are looking for
    };

    auto getDegree = [](const lp_polynomial_context_t* /*ctx*/,
                        lp_monomial_t* m,
                        void* d) {
        degree_travers& v = *static_cast<degree_travers*>(d);

        size_t current_degree = 0;
        // iterate over the number of variables and add up their degrees
        for (size_t i = 0; i < m->n; i++) {
            if (m->p[i].x == v.var) {
                current_degree = m->p[i].d;
                break;
            }
        }
        v.degree = std::max(v.degree, current_degree);
    };

    degree_travers travers;
    travers.var = context().lp_variable(var);
    lp_polynomial_traverse(get_internal(), getDegree, &travers);

    return travers.degree;
}

std::vector<std::size_t>  LPPolynomial::monomial_total_degrees() const {
    struct degree_travers {
        std::vector<std::size_t> degree;
    };

    auto getDegree = [](const lp_polynomial_context_t* /*ctx*/,
                        lp_monomial_t* m,
                        void* d) {
        degree_travers& v = *static_cast<degree_travers*>(d);

        size_t current_degree = 0;
        // iterate over the number of variables and add up their degrees
        for (size_t i = 0; i < m->n; i++) {
            current_degree += m->p[i].d;
        }
        v.degree.push_back(current_degree);
    };

    degree_travers travers;
    lp_polynomial_traverse(get_internal(), getDegree, &travers);

    return travers.degree;
}

std::vector<std::size_t>  LPPolynomial::monomial_degrees(Variable::Arg var) const {
    struct degree_travers {
        std::vector<std::size_t> degree;
        lp_variable_t var; // the variable we are looking for
    };

    auto getDegree = [](const lp_polynomial_context_t* /*ctx*/,
                        lp_monomial_t* m,
                        void* d) {
        degree_travers& v = *static_cast<degree_travers*>(d);

        size_t current_degree = 0;
        // iterate over the number of variables and add up their degrees
        for (size_t i = 0; i < m->n; i++) {
            if (m->p[i].x == v.var) {
                current_degree = m->p[i].d;
                break;
            }
        }
        v.degree.push_back(current_degree);
    };

    degree_travers travers;
    travers.var = context().lp_variable(var);
    lp_polynomial_traverse(get_internal(), getDegree, &travers);

    return travers.degree;
}

std::size_t  LPPolynomial::degree_all_variables() const {
    struct degree_travers {
        std::size_t degree = 0;
    };

    auto getDegree = [](const lp_polynomial_context_t* /*ctx*/,
                        lp_monomial_t* m,
                        void* d) {
        degree_travers& v = *static_cast<degree_travers*>(d);

        // iterate over the number of variables and add up their degrees
        for (size_t i = 0; i < m->n; i++) {
            if (v.degree < m->p[i].d) {
                v.degree = m->p[i].d;
            }
        }
    };

    degree_travers travers;
    lp_polynomial_traverse(get_internal(), getDegree, &travers);

    return travers.degree;
}

mpz_class LPPolynomial::unit_part() const {
    //As we can only have integer coefficients, they do not form a field
    //Thus the unit part is the sign of the leading coefficient, if it is not zero
    //Is the Poly is zero unit part is one
    if(is_zero(*this)) {
        return 1 ; 
    }
    return lp_polynomial_lc_sgn(get_internal()) ;
}

LPPolynomial LPPolynomial::normalized() const {
    auto res = coprime_coefficients();
    auto unit = res.unit_part() ; 
    assert(!is_zero(unit)) ;
    return res * unit ;  

    /*
    auto unit = unit_part() ; 
    assert(!is_zero(unit)) ;
    return (*this) * unit ;  
    */
}

LPPolynomial LPPolynomial::coeff(Variable::Arg var, std::size_t exp) const {
    struct coeff_travers {
        std::vector<lp_monomial_t> coeff;
        const lp_polynomial_context_t* ctx; // context for the newly created monomials
        lp_variable_t var;            // the variable we are looking for
        std::size_t exp;              // the exponent we are looking for
    };


    auto getCoeff = [](const lp_polynomial_context_t* /*ctx*/,
                       lp_monomial_t* m,
                       void* d) {
        coeff_travers& v = *static_cast<coeff_travers*>(d);

        bool found = false;
        // iterate each monomials, and add the ones to add
        for (size_t i = 0; i < m->n; i++) {
            if (m->p[i].x == v.var && m->p[i].d == v.exp) {
                found = true;
                break;
            }
        }

        if (!found) {
            return;
        }

        // Make a copy (without var^exp) of each monomial and add it to the vector
        lp_monomial_t new_monomial;
        lp_monomial_construct(v.ctx, &new_monomial);
        new_monomial.n = m->n -1  ; // copy the number of variables (-1 because we remove var^exp)
        new_monomial.capacity = m->capacity ; // copy the capacity
        lp_integer_assign(v.ctx->K, &new_monomial.a ,&m->a); // copy the coefficient
        new_monomial.p = (power_t*)realloc(new_monomial.p, sizeof(power_t) * new_monomial.capacity); //allocate the memory for the power array

        size_t current_counter = 0 ; 
        for(size_t i = 0; i < m->n; i++){
            if(m->p[i].x == v.var && m->p[i].d == v.exp){
                continue;
            }
            new_monomial.p[current_counter].x = m->p[i].x;
            new_monomial.p[current_counter].d = m->p[i].d;
            current_counter++;
        }
        assert(current_counter == new_monomial.n);
        v.coeff.push_back(std::move(new_monomial));
    };

    coeff_travers travers;
    travers.var = context().lp_variable(var);
    travers.exp = exp;
	travers.ctx = lp_polynomial_get_context(get_internal());
    lp_polynomial_traverse(get_internal(), getCoeff, &travers);

    LPPolynomial res(context());
    for (auto m : travers.coeff) {
        lp_polynomial_add_monomial(res.get_internal(), &m);
    }

    //free the memory allocated for the monomials
    for(auto& m : travers.coeff){
        lp_monomial_destruct(&m);
    }

    return res;
}

bool LPPolynomial::is_normal() const {
    return carl::is_one(this->unit_part()) ; 
    //return carl::is_one(carl::abs(this->unit_part())) ; 
}

std::ostream& operator<<(std::ostream& os, const LPPolynomial& p) {
    os << lp_polynomial_to_string(p.get_internal());
    return os;
}

void LPPolynomial::set_context(const LPContext& c) {
    for (auto& v : variables(*this)) assert(c.has(v));
    if (context() == c) return;

    bool reorder = !(c.is_extension_of(context()) || context().is_extension_of(c));
    m_context = c;
    lp_polynomial_set_context(get_internal(), m_context.lp_context());
    if (reorder) {
        lp_polynomial_ensure_order(get_internal());
    }
    assert(lp_polynomial_check_order(get_internal()));
}

} // namespace carl

#endif