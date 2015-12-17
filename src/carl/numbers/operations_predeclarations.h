#pragma once

#ifndef INCLUDED_FROM_NUMBERS_H
static_assert(false, "This file may only be included indirectly by numbers.h");
#endif

#include <string>
#include "../util/SFINAE.h"

namespace carl {

	template<typename T>
	inline T rationalize(double n);

	template<typename T>
	inline T rationalize(float n);

	template<typename T>
	inline T rationalize(int n);

	template<typename T>
	inline T rationalize(std::size_t n);

	template<typename T>
	inline T rationalize(const std::string& n);

	template<typename From, typename To, carl::DisableIf< std::is_same<From,To> > = dummy >
	inline To convert(const From&);

	template<typename From, typename To, carl::EnableIf< std::is_same<From,To> > = dummy >
	inline To convert(const From& n);

	template<typename Number>
	inline int toInt(const Number& n);
}
