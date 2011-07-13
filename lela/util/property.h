/* lela/util/property.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * C++-structure which mimics C#-style properties
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LELA_UTIL_PROPERTY_H
#define __LELA_UTIL_PROPERTY_H

#include <vector>

namespace LELA
{

/// Closure which dereferences an iterator
template <class Iterator>
struct SimpleAccessor
{
	typedef typename std::iterator_traits<Iterator>::value_type value_type;
	static value_type get (const Iterator &i) { return *i; }
	static void set (Iterator &i, const value_type &T) { *i = T; }
};

/// Closure which gets the second entry in a pair
template <class Iterator>
struct SecondEntryAccessor
{
	typedef typename std::iterator_traits<Iterator>::value_type::second_type value_type;
	static value_type get (const Iterator &i) { return i->second; }
	static void set (Iterator &i, const value_type &T) { i->second = T; }
};

/// This structure mimics a C#-style property.
/// 
/// It wraps an iterator and makes its value appear to be an ordinary
/// value.
///
/// @param Iterator Iterator-type to be wrapped
/// @param Accessor Closure with which to access value of iterator. Normally just dereferences, but may be a more complicated transformation.
template <class Iterator, class Accessor = SimpleAccessor<Iterator> >
struct Property
{
	Iterator _i;

	typedef typename Accessor::value_type value_type;

	Property () {}
	Property (Iterator i) : _i (i) {}

	template <class I2, class A2>
	Property (const Property<I2, A2> &p) : _i (p._i) {}

	Property &operator = (const value_type &v)
		{ Accessor::set (_i, v); return *this; }

	template <class I2, class A2>
	Property &operator = (const Property<I2, A2> &v)
		{ Accessor::set (_i, A2::get (v._i)); return *this; }

	Property &operator += (const value_type &v)
		{ Accessor::set (_i, Accessor::get (_i) + v); return *this; }

	template <class I2, class A2>
	Property &operator += (const Property<I2, A2> &v)
		{ Accessor::set (_i, Accessor::get (_i) + A2::get (v._i)); return *this; }

	Property &operator -= (const value_type &v)
		{ Accessor::set (_i, Accessor::get (_i) - v); return *this; }

	template <class I2, class A2>
	Property &operator -= (const Property<I2, A2> &v)
		{ Accessor::set (_i, Accessor::get (_i) - A2::get (v._i)); return *this; }

	Property &operator *= (const value_type &v)
		{ Accessor::set (_i, Accessor::get (_i) * v); return *this; }

	template <class I2, class A2>
	Property &operator *= (const Property<I2, A2> &v)
		{ Accessor::set (_i, Accessor::get (_i) * A2::get (v._i)); return *this; }

	Property &operator %= (const value_type &v)
		{ Accessor::set (_i, Accessor::get (_i) % v); return *this; }

	template <class I2, class A2>
	Property &operator %= (const Property<I2, A2> &v)
		{ Accessor::set (_i, Accessor::get (_i) % A2::get (v._i)); return *this; }

	Property &operator &= (const value_type &v)
		{ Accessor::set (_i, Accessor::get (_i) & v); return *this; }

	template <class I2, class A2>
	Property &operator &= (const Property<I2, A2> &v)
		{ Accessor::set (_i, Accessor::get (_i) & A2::get (v._i)); return *this; }

	Property &operator |= (const value_type &v)
		{ Accessor::set (_i, Accessor::get (_i) | v); return *this; }

	template <class I2, class A2>
	Property &operator |= (const Property<I2, A2> &v)
		{ Accessor::set (_i, Accessor::get (_i) | A2::get (v._i)); return *this; }

	Property &operator ^= (const value_type &v)
		{ Accessor::set (_i, Accessor::get (_i) ^ v); return *this; }

	template <class I2, class A2>
	Property &operator ^= (const Property<I2, A2> &v)
		{ Accessor::set (_i, Accessor::get (_i) ^ A2::get (v._i)); return *this; }

	operator value_type () const
		{ return Accessor::get (_i); }
};

/// Similar to Property above but shifts the value by a fixed
/// amount. Used in sparse subvectors.

template <class Iterator>
struct ShiftedProperty
{
	typedef typename std::iterator_traits<Iterator>::value_type value_type;

	Iterator _i;
	value_type _shift;

	ShiftedProperty () {}

	ShiftedProperty (Iterator i, value_type shift)
		: _i (i), _shift (shift) {}

	ShiftedProperty &operator = (const value_type &v)
		{ *_i = v + _shift; return *this; }

	operator value_type () const
		{ return *_i - _shift; }
};

/// Similar to Property above but shifts the value by a fixed amount
/// and, treating Iterator::value_type as a pair, works with the first
/// entry. It therefore only supports simple reading and writing. Used
/// in sparse subvectors.

template <class Iterator>
struct ShiftedPropertyFirstEntry
{
	typedef typename std::iterator_traits<Iterator>::value_type::first_type value_type;

	Iterator _i;
	value_type _shift;

	ShiftedPropertyFirstEntry () {}

	ShiftedPropertyFirstEntry (Iterator i, value_type shift)
		: _i (i), _shift (shift) {}

	ShiftedPropertyFirstEntry &operator = (const value_type &v)
		{ _i->first = v + _shift; return *this; }

	operator value_type () const
		{ return _i->first - _shift; }
};

} // namespace LELA

#endif // __LELA_UTIL_PROPERTY_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax