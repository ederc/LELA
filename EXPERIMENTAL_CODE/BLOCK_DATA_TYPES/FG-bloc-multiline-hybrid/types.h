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
#include <assert.h>

#include "Allocator.h"

#include "lela/vector/sparse.h"
//#include "matrix-utils.h"

using namespace LELA;

#ifndef DEFAULT_BLOC_HEIGHT
#define DEFAULT_BLOC_HEIGHT 0
#error "must define DEFAULT_BLOC_HEIGHT"
#endif

#ifndef DEFAULT_BLOC_WIDTH
#define DEFAULT_BLOC_WIDTH 0
#error "must define DEFAULT_BLOC_WIDTH"
#endif

//#ifndef HYBRID_REPRESENTATION_THRESHOLD
#define HYBRID_REPRESENTATION_THRESHOLD 0.5f
//#error "must define HYBRID_REPRESENTATION_THRESHOLD"
//#endif

#ifndef NB_ROWS_PER_MULTILINE
#define NB_ROWS_PER_MULTILINE 2
#endif


template <typename Element, typename Index = uint16>
class MultiLineVector {
public:

	MultiLineVector() : _bloc_height (NB_ROWS_PER_MULTILINE)
	{
		IndexData._index_vector = std::vector<Index, aligned_allocator<Index, 16> > ();
		ValuesData._data = std::vector<Element, aligned_allocator<Element, 16> > ();
	}

	MultiLineVector(uint16 nb_lines) : _bloc_height (nb_lines)
	{
		IndexData._index_vector = std::vector<Index, aligned_allocator<Index, 16> > ();
		ValuesData._data = std::vector<Element, aligned_allocator<Element, 16> > ();
	}

//	MultiLineVector( const MultiLineVector& source ):
//					_bloc_height (source._bloc_height)
//	{
//			IndexData._index_vector = source.IndexData._index_vector;
//			ValuesData._data = source.ValuesData._data;
//	}

	inline void		reserve		(size_t sz)		{ IndexData._index_vector.reserve (sz); ValuesData._data.reserve (sz*_bloc_height); }

	inline bool		empty		()	const		{ return ValuesData._data.empty (); }

	inline size_t size() const
	{
		return ValuesData._data.size() / _bloc_height;
	}
	uint16			nb_lines	()	const		{ return _bloc_height; }

	inline void clear()
	{
		IndexData._index_vector.clear();
		ValuesData._data.clear();
	}

	inline Element	at			(uint16 line_index, size_t n)	const
	{
		if (line_index >= _bloc_height)
		{
			std::cerr << "std::out_of_range (line_index) " << line_index << "[_bloc_height = " << _bloc_height << "]" <<std::endl;
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

//	inline void		at_unchecked2(size_t n, Element& val1, Element& val2)	const
//	{
//		register uint32 idx = _bloc_height * n;
//		val1 = ValuesData._data[idx];
//		val2 = ValuesData._data[idx+1];
//	}
//
//	inline void		at_unchecked2(size_t n, uint32& val1, uint32& val2)	const
//	{
//		//register uint32 idx = _bloc_height * n;
//		val1 = ValuesData._data[_bloc_height * n];
//		val2 = ValuesData._data[_bloc_height * n+1];
//	}

//	inline Element* at_unchecked_pointer(uint16 line_index, size_t n)
//	{
//		return &(ValuesData._data[line_index + _bloc_height * n]);
//	}

	inline void push_back(uint32 index, Element e1, Element e2)
	{
		IndexData.push_back(index);
		ValuesData.push_back(e1);
		ValuesData.push_back(e2);
	}

	struct IndexData
	{
		typedef typename std::vector<Index>::iterator iterator;

		iterator	begin ()		{ return _index_vector.begin (); }
		iterator	end   ()		{ return _index_vector.end (); }

		inline Index		&operator []	(size_t i) { return _index_vector[i]; }
		inline const Index	&operator []	(size_t i)	const { return _index_vector[i]; }

		inline void		push_back (Index i)	{ _index_vector.push_back (i); }

		std::vector<Index, aligned_allocator<Index, 16> >		_index_vector;
		
		Index* 	getStartingPointer ()	{ return &(_index_vector[0]); }
		const Index* 	getStartingPointer ()	const	{ return &(_index_vector[0]); }

		void* 	getStartingPointerAsVoid ()	{ return (void *)&(_index_vector[0]); }
		const void* 	getStartingPointerAsVoid ()	const	{ return (void *)&(_index_vector[0]); }

	} IndexData;

	struct ValuesData
	{
		typedef typename std::vector<Element>::iterator iterator;

		iterator	begin ()		{ return this->_data.begin (); }
		iterator	end   ()		{ return this->_data.end (); }

		inline Element		&operator []	(size_t i) { return _data[i]; }
		inline const Element	&operator []	(size_t i) const { return _data[i]; }

		inline void		push_back (Element e)	{ this->_data.push_back (e); }

		std::vector<Element, aligned_allocator<Element, 16> >	_data;
		
		Element*		getStartingPointer ()	{ return &(_data[0]); }
		const Element*		getStartingPointer ()	const	{ return &(_data[0])	; }

		void*		getStartingPointerAsVoid ()	{ return (void *)&(_data[0]); }
		const void*		getStartingPointerAsVoid ()	const	{ return (void *)&(_data[0])	; }

	} ValuesData;

	void swap(MultiLineVector<Element, Index>& other)
	{
		this->_bloc_height = other._bloc_height;
		std::swap(this->IndexData._index_vector, other.IndexData._index_vector);
		std::swap(this->ValuesData._data, other.ValuesData._data);
	}

	static inline float get_HYBRID_REPRESENTATION_THRESHOLD ()
	{
		return (float) HYBRID_REPRESENTATION_THRESHOLD;
	}

	inline bool is_sparse() const
	{
		return ValuesData._data.size() < DEFAULT_BLOC_WIDTH * 2;
	}

private:
	uint16					_bloc_height;			//number of lines per bloc

};




template <typename Element, typename Index = uint16>
class SparseMultilineBloc {

public:
	typedef Index IndexType;

	typedef MultiLineVector<Element, Index> Row;
	typedef const Row ConstRow;
	typedef std::vector<Row> Rep;

	typedef typename Rep::iterator RowIterator;
	typedef typename Rep::const_iterator ConstRowIterator;

	SparseMultilineBloc() :
			_bloc_height(DEFAULT_BLOC_HEIGHT / NB_ROWS_PER_MULTILINE),
			_bloc_width(DEFAULT_BLOC_WIDTH)
	{
		//_A = Rep(_bloc_height);
	}

	SparseMultilineBloc (uint16 bloc_height, uint16 bloc_width) :
				_bloc_height (bloc_height/NB_ROWS_PER_MULTILINE),
				_bloc_width (bloc_width)
	{
		if(_bloc_height < 1)
			throw std::invalid_argument ("bloc_height");

		if(_bloc_width < 1)
			throw std::invalid_argument ("bloc_width");

		_A = Rep (_bloc_height);
	}

	SparseMultilineBloc (uint16 bloc_height) :
				_bloc_height (bloc_height/NB_ROWS_PER_MULTILINE),
				_bloc_width (DEFAULT_BLOC_WIDTH)
	{
		if(_bloc_height < 1)
			throw std::invalid_argument ("bloc_height");

		if(_bloc_width < 1)
			throw std::invalid_argument ("bloc_width");

		_A = Rep (_bloc_height);
	}


	SparseMultilineBloc(const SparseMultilineBloc<Element, Index>& other) :
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

		_bloc_height = height/NB_ROWS_PER_MULTILINE;
		_bloc_width = width;

		_A = Rep (_bloc_height);
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

	static inline float get_HYBRID_REPRESENTATION_THRESHOLD ()
	{
		return (float) HYBRID_REPRESENTATION_THRESHOLD;
	}

private:
	uint16					_bloc_height;			//number of lines per bloc
	uint16					_bloc_width;
	Rep						_A;
};

template <typename BlocType_>
class SparseBlocMatrix {

public:
	bool acceptRowsHybrid;

	typedef std::vector<BlocType_> Row;
	typedef const Row ConstRow;
	typedef std::vector<Row> Rep;

	typedef BlocType_ BlocType;

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
			acceptRowsHybrid (true),
			_m(0),
			_n(0),
			_bloc_height(DEFAULT_BLOC_HEIGHT),
			_bloc_width(DEFAULT_BLOC_WIDTH) {}


	SparseBlocMatrix (const SparseBlocMatrix<BlocType>& other) :
							blocArrangement (other.blocArrangement),
							_m (other._m),
							_n (other._n),
							_bloc_height (other._bloc_height),		//number of lines per bloc
							_bloc_width (other._bloc_width),
							_fill_with_empty_blocs (other._fill_with_empty_blocs),
							_A (other._A),
							acceptRowsHybrid (true){}

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

		acceptRowsHybrid = true;
	}

	size_t rowdim () const { return _n; }
	size_t coldim () const { return _m; }
	size_t rowBlocDim () const { return _row_block_dim; }

	uint16 bloc_width ()	const	{ return _bloc_width; }
	uint16 bloc_height ()	const	{ return _bloc_height; }

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

	static inline float get_HYBRID_REPRESENTATION_THRESHOLD ()
	{
		return (float) HYBRID_REPRESENTATION_THRESHOLD;
	}

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
								_bloc_height (NB_ROWS_PER_MULTILINE)		//number of lines per row
	{
		if(_bloc_height<1)
			throw std::logic_error ("bloc_height");

		//MultiLineVector<Element> tmp_copy (_bloc_height);

		if(n%_bloc_height == 0)					// (n/_bloc_height) lines of blocs
			_A = Rep (n/_bloc_height);
		else
			_A = Rep (n/_bloc_height + 1);
	}

	SparseMultilineMatrix (size_t n, size_t m, uint16 bloc_height) :
								_m (m),
								_n (n),
								_bloc_height (bloc_height)
	{
		if(_bloc_height<1)
			throw std::logic_error ("bloc_height");

		//MultiLineVector<Element> tmp_copy (_bloc_height);

		if(n%_bloc_height == 0)						// (n/_bloc_height) lines of blocs
			_A = Rep (n/_bloc_height);
		else
			_A = Rep (n/_bloc_height + 1);
	}

	SparseMultilineMatrix (const SparseMultilineMatrix<Element>& other) :
						_m (other._m),
						_n (other._n),
						_bloc_height (other._bloc_height),
						_A (other._A) {}

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

#endif /* TYPES_FG_BLOC_H_ */
