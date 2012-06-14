/*
 * FG-types.h
 *
 *  Created on: 14 juin 2012
 *      Author: martani
 */

#ifndef FG_TYPES_H_
#define FG_TYPES_H_

#include <vector>
#include "lela/vector/sparse.h"


using namespace LELA;

template <typename Element>
class MultiLineVector {
public:

	MultiLineVector() : _bloc_height (2) {}
	MultiLineVector(uint16 nb_lines) : _bloc_height (nb_lines) {}

	inline void		reserve		(size_t sz)		{ IndexData._index_vector.reserve (sz); ValuesData._data.reserve (sz*_bloc_height); }
	inline bool		empty		()	const		{ return IndexData._index_vector.empty (); }

	size_t 			size 		()	const		{ return IndexData._index_vector.size (); }
	uint16			nb_lines	()	const		{ return _bloc_height; }

	inline Element	at			(uint16 line_index, size_t n)	const
	{
		if (line_index >= _bloc_height)
		{
			std::cerr << "std::out_of_range (line_index) " << line_index << std::endl;
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

	struct IndexData
	{
		typedef typename std::vector<uint32>::iterator iterator;

		iterator	begin ()		{ return _index_vector.begin (); }
		iterator	end   ()		{ return _index_vector.end (); }

		uint32		&operator []	(size_t i) { return _index_vector[i]; }
		const uint32		&operator []	(size_t i)	const { return _index_vector[i]; }

		void		push_back (uint32 i)	{ _index_vector.push_back (i); }

		std::vector<uint32>		_index_vector;

	} IndexData;

	struct ValuesData
	{
		typedef typename std::vector<Element>::iterator iterator;

		iterator	begin ()		{ return this->_data.begin (); }
		iterator	end   ()		{ return this->_data.end (); }

		uint32		&operator []	(size_t i) { return _data[i]; }
		uint32		&operator []	(size_t i) const { return _data[i]; }

		void		push_back (Element e)	{ this->_data.push_back (e); }

		std::vector<Element>	_data;

	} ValuesData;

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

		if(n%_bloc_height == 0)						// (n/_bloc_height) lines of blocs
			_A = Rep (n/_bloc_height);
		else
			_A = Rep (n/_bloc_height + 1);
	}

	SparseMultilineMatrix (size_t n, size_t m, uint16 bloc_height, uint16 bloc_width) :
								_m (m),
								_n (n),
								_bloc_height (bloc_height)
	{
		if(_bloc_height<1)
			throw WrongParameter ();

		if(n%_bloc_height == 0)						// (n/_bloc_height) lines of blocs
			_A (n/_bloc_height);
		else
			_A (n/_bloc_height + 1);
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
};























#endif /* FG_TYPES_H_ */
