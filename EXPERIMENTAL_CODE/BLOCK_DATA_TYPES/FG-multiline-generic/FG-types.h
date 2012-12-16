/*
 * FG-types.h
 *
 *  Created on: 14 juin 2012
 *      Author: martani
 */

#ifndef FG_TYPES_H_
#define FG_TYPES_H_

#include <iostream>
#include <vector>
#include <assert.h>
#include "lela/vector/sparse.h"


using namespace LELA;

template <typename Element>
class MultiLineVector {
public:

	MultiLineVector() : _bloc_height (2) {}
	MultiLineVector(uint16 nb_lines) : _bloc_height (nb_lines) {}
	MultiLineVector( const MultiLineVector& source ):
					_bloc_height(source._bloc_height)
	{
			assert(_bloc_height > 0);
			IndexData._index_vector = std::vector<uint32> ();
			ValuesData._data = std::vector<Element> ();
	}

	inline void		reserve		(size_t sz)		{ IndexData._index_vector.reserve (sz); ValuesData._data.reserve (sz*_bloc_height); }
	inline bool		empty		()	const		{ return IndexData._index_vector.empty (); }

	inline size_t 	size 		()	const		{ return IndexData._index_vector.size (); }
	uint16			nb_lines	()	const		{ return _bloc_height; }

	inline Element	at			(uint16 line_index, size_t n)	const
	{
		if (line_index >= _bloc_height)
		{
			std::cerr << "std::out_of_range (line_index)=" << line_index << " [_bloc_height=" << _bloc_height << "]" <<std::endl;
			throw std::out_of_range ("line_index");
		}

		if (n >= size ())
		{
			std::cerr << "std::out_of_range (n)";
			throw std::out_of_range ("n");
		}
		else
			return ValuesData._data[line_index + _bloc_height * n];
	}

	inline Element	at_unchecked(uint16 line_index, size_t n)	const
	{
		return ValuesData._data[line_index + _bloc_height * n];
	}


	inline Element* at_unchecked_pointer(uint16 line_index, size_t n)
	{
		return &(ValuesData._data[line_index + _bloc_height * n]);
	}

	struct IndexData
	{
		typedef typename std::vector<uint32>::iterator iterator;

		iterator	begin ()		{ return _index_vector.begin (); }
		iterator	end   ()		{ return _index_vector.end (); }

		inline uint32		&operator []	(size_t i) { return _index_vector[i]; }
		inline const uint32		&operator []	(size_t i)	const { return _index_vector[i]; }

		inline void		push_back (uint32 i)	{ _index_vector.push_back (i); }

		std::vector<uint32>		_index_vector;

	} IndexData;

	struct ValuesData
	{
		typedef typename std::vector<Element>::iterator iterator;

		iterator	begin ()		{ return this->_data.begin (); }
		iterator	end   ()		{ return this->_data.end (); }

		inline Element		&operator []	(size_t i) { return _data[i]; }
		inline const Element		&operator []	(size_t i) const { return _data[i]; }

		inline void		push_back (Element e)	{ this->_data.push_back (e); }

		std::vector<Element>	_data;

	} ValuesData;

	void swap(MultiLineVector<Element>& other)
	{
		this->_bloc_height = other._bloc_height;
		std::swap(this->IndexData._index_vector, other.IndexData._index_vector);
		std::swap(this->ValuesData._data, other.ValuesData._data);
	}

private:
	uint16					_bloc_height;			//number of lines per bloc

};

class WrongParameter {};

template <typename Element>
class SparseMultilineMatrix {

public:
	typedef MultiLineVector<Element> Row;
	typedef const Row ConstRow;
	typedef std::vector<Row> Rep;

	SparseMultilineMatrix () : _m (0), _n (0) {}
	SparseMultilineMatrix (size_t n, size_t m) :
								_m (m),
								_n (n),
								_bloc_height (2)		//number of lines per row
	{
		if(_bloc_height<1)
			throw WrongParameter ();

		MultiLineVector<Element> tmp_copy (_bloc_height);

		if(n%_bloc_height == 0)						// (n/_bloc_height) lines of blocs
			_A = Rep (n/_bloc_height, tmp_copy);
		else
			_A = Rep (n/_bloc_height + 1, tmp_copy);
	}

	SparseMultilineMatrix (size_t n, size_t m, uint16 bloc_height) :
								_m (m),
								_n (n),
								_bloc_height (bloc_height)
	{
		if(_bloc_height<1)
			throw WrongParameter ();

		MultiLineVector<Element> tmp_copy (_bloc_height);

		if(n%_bloc_height == 0)						// (n/_bloc_height) lines of blocs
			_A = Rep (n/_bloc_height, tmp_copy);
		else
			_A = Rep (n/_bloc_height + 1, tmp_copy);
	}

	~SparseMultilineMatrix () {}

	size_t rowdim () const { return _n; }
	size_t coldim () const { return _m; }

	size_t multiline_rowdim () const
	{
		if(rowdim () % _bloc_height == 0)						// (n/_bloc_height) lines of blocs
			return rowdim () / _bloc_height;
		else
			return rowdim () / _bloc_height + 1;
	}

	uint16 nb_lines_per_bloc 	() const { return _bloc_height; }

	typedef typename Rep::iterator RowIterator;
	typedef typename Rep::const_iterator ConstRowIterator;

	RowIterator      rowBegin ()			{ return _A.begin (); }
	ConstRowIterator rowBegin ()	const	{ return _A.begin (); }
	ConstRowIterator rowEnd ()		const   { return _A.end (); }
	RowIterator      rowEnd ()				{ return _A.end (); }

	Row			&getRow			(size_t i) { return _A[i]; }
	Row			&operator []	(size_t i) { return _A[i]; }
	ConstRow	&operator []	(size_t i) const { return _A[i]; }

private:
	Rep		_A;
	size_t		_m;
	size_t		_n;

	uint16		_bloc_height;

	SparseMultilineMatrix (const SparseMultilineMatrix& other) {}
};























#endif /* FG_TYPES_H_ */
