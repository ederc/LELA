/*
 * types.h
 *
 *  Created on: 4 juil. 2012
 *      Author: martani
 */

#ifndef TYPES_FG_BLOC_H_
#define TYPES_FG_BLOC_H_

#include <iostream>
#include <vector>
#include "lela/vector/sparse.h"


using namespace LELA;

#ifndef DEFAULT_BLOC_HEIGHT
#define DEFAULT_BLOC_HEIGHT -1
#endif

#ifndef DEFAULT_BLOC_WIDTH
#define DEFAULT_BLOC_WIDTH -1
#endif

template <typename Element, typename Index = uint16>
class SparseBloc {

public:
	typedef Index IndexType;

	typedef SparseVector <Element, std::vector<Index>, std::vector<Element> > Row;
	//typedef SparseVector<Element, std::vector<Index> > Row;
	typedef const Row ConstRow;
	typedef std::vector<Row> Rep;

	typedef typename Rep::iterator RowIterator;
	typedef typename Rep::const_iterator ConstRowIterator;

	SparseBloc () : _bloc_height (DEFAULT_BLOC_HEIGHT), _bloc_width (DEFAULT_BLOC_WIDTH)
	{
		_A = Rep (_bloc_height);
	}

	SparseBloc (uint16 bloc_height, uint16 bloc_width) :
				_bloc_height (bloc_height),
				_bloc_width (bloc_width)
	{
		if(_bloc_height < 1)
			throw std::invalid_argument ("bloc_height");

		if(_bloc_width < 1)
			throw std::invalid_argument ("bloc_width");

		Row tmp;
		_A = Rep (_bloc_height);
	}

	SparseBloc (uint16 bloc_height) :
				_bloc_height (bloc_height),
				_bloc_width (DEFAULT_BLOC_WIDTH)
	{
		if(_bloc_height < 1)
			throw std::invalid_argument ("bloc_height");

		if(_bloc_width < 1)
			throw std::invalid_argument ("bloc_width");

		_A = Rep (_bloc_height);
	}


	SparseBloc(const SparseBloc<Element, Index>& other) :
				_bloc_height (other._bloc_height),
				_bloc_width (other._bloc_width),
				_A (other._A)
	{

	}

	void	init	(uint16 height, uint16 width)
	{
		if(height < 1)
			throw std::invalid_argument ("bloc_height");

		if(width < 1)
			throw std::invalid_argument ("bloc_width");

		_bloc_height = height;
		_bloc_width = width;

		_A = Rep (_bloc_height);

		/*for(uint32 i=0; i<_bloc_height; ++i)
			_A[i].reserve (_bloc_width);*/

		//std::cout <<" in init " << std::endl;
	}

	uint16	bloc_height	()	const	{ return _bloc_height; }
	uint16	bloc_width	()	const	{ return _bloc_width; }

	RowIterator      rowBegin ()			{ return _A.begin (); }
	ConstRowIterator rowBegin ()	const	{ return _A.begin (); }
	ConstRowIterator rowEnd ()		const   { return _A.end (); }
	RowIterator      rowEnd ()				{ return _A.end (); }

	Row			&getRow			(size_t i) { return _A[i]; }
	Row			&operator []	(size_t i) { return _A[i]; }
	ConstRow	&operator []	(size_t i) const { return _A[i]; }

	size_t	size	()	const	{ return _A.size(); }
	bool	empty	()	const	{ return _A.empty(); }

private:
	uint16					_bloc_height;			//number of lines per bloc
	uint16					_bloc_width;
	//size_t					_bloc_start_colum_index;		//the row of blocs should keep this info
	Rep						_A;
	bool					_initiliazed;
};


template <typename Element, typename Index = uint16>
class SparseBlocMatrix {

public:

	typedef std::vector<SparseBloc<Element, Index> > Row;
	typedef const Row ConstRow;
	typedef std::vector<Row> Rep;

	typedef SparseBloc<Element, Index> BlocType;

	typedef typename Rep::iterator RowIterator;
	typedef typename Rep::const_iterator ConstRowIterator;

	/**
	 * Enumeration that specifies the arrangement of blocs inside the matrix.
	 * Blocs are arranged by lines:
	 * 		- example: from Top to Down and from Left to Right: any left rows (resp columns) are added to the
	 * 		last row of blocs (resp right most blocs)
	 */
	enum BlocArrangement {
		ArrangementTopDown_LeftRight,
		ArrangementTopDown_RightLeft,
		ArrangementDownTop_LeftRight,
		ArrangementDownTop_RightLeft
	};

	BlocArrangement blocArrangement;

	SparseBlocMatrix() :
			_m(0),
			_n(0),
			_bloc_height(DEFAULT_BLOC_HEIGHT),
			_bloc_width(DEFAULT_BLOC_WIDTH) {}

	SparseBlocMatrix (const SparseBlocMatrix<Element, Index>& other) :
							blocArrangement (other.blocArrangement),
							_m (other._m),
							_n (other._n),
							_bloc_height (other._bloc_height),		//number of lines per bloc
							_bloc_width (other._bloc_width),
							_fill_with_empty_blocs (other._fill_with_empty_blocs),
							_A (other._A) {}

	/**
	 * @param fill_with_empty_blocs: fills the empty areas of the matrix with empty (not initiliazed blocs)
	 * this parameter is to be use especially in the case of read-write matrices. If the matrix is read only,
	 * this parameter should be false to preserve memory.
	 * When false, blocs are added only where needed.
	 */
	SparseBlocMatrix (size_t n, size_t m,
					  BlocArrangement blocs_arrangement = ArrangementTopDown_LeftRight,
					  bool fill_with_empty_blocs = false,
					  uint16 bloc_height=DEFAULT_BLOC_HEIGHT,
					  uint16 bloc_width=DEFAULT_BLOC_WIDTH) :
						  	  	blocArrangement (blocs_arrangement),
								_m (m),
								_n (n),
								_bloc_height (bloc_height),		//number of lines per bloc
								_bloc_width (bloc_width),
								_fill_with_empty_blocs (fill_with_empty_blocs)
	{
		if(_bloc_height < 1)
			throw std::invalid_argument ("bloc_height");

		if(_bloc_width < 1)
			throw std::invalid_argument ("bloc_width");


		if(n%_bloc_height == 0)
		{
			_A = Rep (n/_bloc_height);
			FirstBlocsColumIndexes.init(n/_bloc_height, 0);
			_row_block_dim = n/_bloc_height;
		}
		else
		{
			_A = Rep (n/_bloc_height + 1);
			FirstBlocsColumIndexes.init(n/_bloc_height + 1, 0);
			_row_block_dim = n/_bloc_height + 1;
		}
	}

	size_t rowdim () const { return _n; }
	size_t coldim () const { return _m; }
	size_t rowBlocDim () const { return _row_block_dim; }

	uint16 block_width ()	const	{ return _bloc_width; }
	uint16 block_height ()	const	{ return _bloc_height; }

	RowIterator      rowBegin ()			{ return _A.begin (); }
	ConstRowIterator rowBegin ()	const	{ return _A.begin (); }
	ConstRowIterator rowEnd ()		const   { return _A.end (); }
	RowIterator      rowEnd ()				{ return _A.end (); }

	Row			&getRow			(size_t i) { return _A[i]; }
	Row			&operator []	(size_t i) { return _A[i]; }
	ConstRow	&operator []	(size_t i) const { return _A[i]; }

	struct FirstBlocsColumIndexes {
	public:
		void	init	(size_t sz, uint32 val)		{ _blocs_start_colum_indexes = std::vector<uint32> (sz, val); }

		uint32		&operator [] (size_t i)			{ return _blocs_start_colum_indexes[i]; }
		const uint32		&operator [] (size_t i)	const	{ return _blocs_start_colum_indexes[i]; }

	private:
		std::vector<uint32>	_blocs_start_colum_indexes;
	} FirstBlocsColumIndexes;

	bool		isFilledWithEmptyBlocs()	const	{ return _fill_with_empty_blocs; }

private:
	Rep					_A;
	//std::vector<uint32>	_blocs_start_colum_indexes;
	size_t				_m;
	size_t				_n;

	size_t				_row_block_dim;

	uint16				_bloc_height;
	uint16				_bloc_width;

	bool				_fill_with_empty_blocs;
};


#endif /* TYPES_FG_BLOC_H_ */
