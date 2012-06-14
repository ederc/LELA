#ifndef __BLOC_TYPES_H
#define __BLOC_TYPES_H

#include <vector>

#include "lela/vector/sparse.h"


using namespace LELA;

namespace BlocTypes {
	struct Sparse {};
	struct Dense {};
};

namespace FG_Types {
	
class UnknownVectorType {};
class WrongParameter {};

template <typename Element, typename Index = uint32>
class HybridVector {

public:
	typedef LELA::SparseVector<Element, std::vector<Index>, std::vector<Element> > _SparseVector;
	typedef std::vector<Element> _DenseVector;
	
	enum VectorTypeEnum {
		SparseVectorType,
		DenseVectorType
	};
	VectorTypeEnum VectorType;

	typedef size_t size_type;

	HybridVector ()
	{
		this->VectorType = SparseVectorType;
		initVectorStorage(this->VectorType);
	}

	HybridVector(VectorTypeEnum type)
	{
		this->VectorType = type;
		initVectorStorage(this->VectorType);
	}

	HybridVector(VectorTypeEnum type, size_t size)
	{
		this->VectorType = type;
		initVectorStorage(this->VectorType);

		if(this->VectorType == SparseVectorType)
			this->_sparse.reserve(size);
		else if(this->VectorType == SparseVectorType)
			this->_dense.reserve(size);
		else
			throw UnknownVectorType ();

	}

	struct SparseRep
	{
		typedef typename _SparseVector::iterator iterator;

		iterator	begin ()		{ return this->_sparse.begin (); }
		iterator	end   ()		{ return this->_sparse.end (); }
	} SparseRep;

	struct DenseRep
	{
		typedef typename _DenseVector::iterator iterator;

		iterator	begin ()		{ return this->_dense.begin (); }
		iterator	end   ()		{ return this->_dense.end (); }
	} DenseRep;

	inline size_type	size ()	const
	{
		switch (VectorType) {
			case SparseVectorType:
				return this->_sparse.size ();
				break;
			case DenseVectorType:
				return this->_dense.size ();
				break;
			default:
				throw UnknownVectorType ();
				break;
		}

		return 0;
	}

	inline size_type	nbNonZeroElements	()	const
	{
		switch (VectorType) {
			case SparseVectorType:
				return this->_sparse.size ();
				break;
			case DenseVectorType:
				{
					size_type sz = 0;
					for (uint32 i = 0; i < this->_dense.size (); ++i) {
						if(this->_dense[i] != 0)
							sz++;
					}

					return sz;
				}
				break;
			default:
				throw UnknownVectorType ();
				break;
		}

		return 0;
	}

	inline bool		empty	()	const
	{
		switch (VectorType) {
			case SparseVectorType:
				return this->_sparse.empty ();
				break;
			case DenseVectorType:
				return this->_dense.empty ();
				break;
			default:
				throw UnknownVectorType ();
				break;
		}
		return true;
	}

	inline bool operator ==(const _SparseVector& v) const {
		switch (VectorType) {
			case SparseVectorType:
				return this->_sparse == v;
				break;

			case DenseVectorType:
				return sparseDenseVectorEqual(v, this->_dense);
				break;

			default:
				throw UnknownVectorType ();
				break;
		}

		return false;
	}

	inline bool operator ==(const _DenseVector& v) const {
		switch (VectorType) {
			case SparseVectorType:
				return sparseDenseVectorEqual(this->_sparse, v);
				break;

			case DenseVectorType:
				return this->_dense == v;
				break;

			default:
				throw UnknownVectorType ();
				break;
		}

		return false;
	}

	inline bool operator == (const HybridVector<Element, Index>& v) const
	{
		switch (VectorType) {
			case SparseVectorType:
				if(v.VectorType == SparseVectorType)
					return this->_sparse == v._sparse;
				else if(v.VectorType == DenseVectorType)
					return sparseDenseVectorEqual(this->_sparse, v._dense);
				else
					throw UnknownVectorType ();

				break;
			case DenseVectorType:
				if(v.VectorType == SparseVectorType)
					return sparseDenseVectorEqual(v._sparse, this->_dense);
				else if(v.VectorType == DenseVectorType)
					return this->_dense == v._dense;
				else
					throw UnknownVectorType ();

				break;
			default:
				throw UnknownVectorType ();
				break;
		}

		return false;
	}


	void swap (HybridVector<Element, Index>& v)
	{
		switch (VectorType) {
			case SparseVectorType:
				if(v.VectorType == SparseVectorType)
				    this->_sparse.swap(v._sparse);
				else if(v.VectorType == DenseVectorType)
					{
						this->VectorType = DenseVectorType;
						this->_sparse.clear ();
						_SparseVector _tmp_sparse;
						this->_sparse.swap(_tmp_sparse);

						this->_dense.swap(v._dense);
					}
				else
					throw UnknownVectorType ();

				break;
			case DenseVectorType:
				if(v.VectorType == SparseVectorType)
					{
						this->VectorType = SparseVectorType;
						this->_dense.clear ();
						_DenseVector _tmp_dense;
						this->_dense.swap(_tmp_dense);

						this->_sparse.swap(v._sparse);

					}
				else if(v.VectorType == DenseVectorType)
					this->_dense.swap(v._dense);
				else
					throw UnknownVectorType ();

				break;
			default:
				throw UnknownVectorType ();
				break;
		}

	}

	inline size_type		capacity	() const
	{
		switch (VectorType) {
			case SparseVectorType:
				return this->_sparse.capacity ();
				break;
			case DenseVectorType:
				return this->_dense.capacity ();
				break;
			default:
				throw UnknownVectorType ();
				break;
		}
		return 0;
	}

	inline void				reserve		(size_type s)
	{
		switch (VectorType) {
			case SparseVectorType:
				this->_sparse.reserve (s);
				break;
			case DenseVectorType:
				this->_dense.reserve (s);
				break;
			default:
				throw UnknownVectorType ();
				break;
		}
	}

	void	free	()
	{
		if(this->VectorType == SparseVectorType)
		{
			_SparseVector tmp;
			this->_sparse.clear ();
			this->_sparse.swap (tmp);
		}
		else if(this->VectorType == SparseVectorType)
		{
			_DenseVector tmp;
			this->_dense.clear ();
			this->_dense.swap (tmp);
		}
		else
			throw UnknownVectorType ();
	}

	template <typename Ring>
	int	head (Ring& R, Element& a)	const
	{
		switch (VectorType) {
			case SparseVectorType:
				return head(R, a, this->_sparse);
				break;

			case DenseVectorType:
				return head(R, a, this->_dense);
				break;

			default:
				throw UnknownVectorType ();
				break;
		}

		return -1;
	}

	int	head ()	const
	{
		switch (VectorType) {
			case SparseVectorType:
				return head(this->_sparse);
				break;

			case DenseVectorType:
				return head(this->_dense);
				break;

			default:
				throw UnknownVectorType ();
				break;
		}

		return -1;
	}

private:
	_SparseVector _sparse;	// takes 24*2=48 bytes
	_DenseVector _dense;	// takes 24 bytes

	void initVectorStorage(VectorTypeEnum type)
	{
		_SparseVector _tmp_sparse;
		_DenseVector _tmp_dense;

		switch (type) {
			case SparseVectorType:
				_sparse.clear ();
				_sparse.swap(_tmp_sparse);

				_dense.clear ();
				_dense.swap(_tmp_dense);
				break;
			case DenseVectorType:
				_sparse.clear ();
				_sparse.swap(_tmp_sparse);

				_dense.clear ();
				_dense.swap(_tmp_dense);
				break;
			default:
				break;
		}
	}

	//Notice that vectors could have different nb of elements!
	static bool sparseDenseVectorEqual(_SparseVector sparse, _DenseVector dense)
	{
		typename _SparseVector::const_iterator it = sparse.begin ();
		uint32 i;

		for(i = 0; i < dense.size (); ++i)
		{
			if(dense[i] == 0)
			{
				if(it == sparse.end () || i < it->first)
					continue;

				if(i >= it->first)
					return false;
			}


			if(i == it->first && dense[i] == it->second)
			{
				++it;
				continue;
			}

			return false;
		}

		if(it == sparse.end ())
			return true;
		else
			return false;
	}

	template <typename Ring>
	static inline int head(Ring& R, Element& a, _SparseVector& v)
	{
		if(v.empty ())
			return -1;
		else
		{
			R.copy(a, v.front ().second);
			return v.front ().first;
		}
	}

	template <typename Ring>
	static inline int head(Ring& R, Element& a, _DenseVector& v)
	{
		for (int i = 0; i < v.size(); ++i) {
			if(v[i] != 0)
			{
				R.copy(a, v[i]);
				return i;
			}
		}

		return -1;
	}

	static inline int head(_SparseVector& v)
	{
		if(v.empty ())
			return -1;
		else
			return v.front ().first;
	}

	static inline int head(_DenseVector& v)
	{
		for (int i = 0; i < v.size(); ++i) {
			if(v[i] != 0)
				return i;
		}

		return -1;
	}

};

///A sparse bloc is simply a container vector of sparse vectors.
///It represents a line of sparse blocs in a matrix
template <typename Vector>
class Bloc {

public:
	
	typedef Vector Row;
	typedef const Row ConstRow;
	typedef std::vector<Row> Rep;

	/*Bloc () : _data (0), _height (0), _width (0) {}
	Bloc (size_t n, size_t m) : _data (n), _height (n), _width (m) {}*/

	Bloc () : _data (0) {}
	Bloc (size_t n) : _data (n) {}

	~Bloc () {}

	/*size_t rowdim () const { return _height; }
	size_t coldim () const { return _width; }*/

	typedef typename Rep::iterator RowIterator;
	typedef typename Rep::const_iterator ConstRowIterator;

	RowIterator      rowBegin ()			{ return _data.begin (); }
	ConstRowIterator rowBegin ()	const	{ return _data.begin (); }
	ConstRowIterator rowEnd ()		const	{ return _data.end (); }
	RowIterator      rowEnd ()				{ return _data.end (); }

	Row			&getRow			(size_t i)	{ return _data[i]; }
	Row			&operator []	(size_t i)	{ return _data[i]; }
	ConstRow	&operator []	(size_t i)	const	{ return _data[i]; }

	void	init (uint16 n)
	{
		_data = Rep (n);
		/*_height = n;
		_width = m;*/
	}

private:
	Rep		_data;
	/*size_t		_height;
	size_t		_width;*/
};


///A SparseBlocMatrix is represented as a vector of vectors of sparse blocs
///The blocs are arranged horizontally by lines (vector<SparseBloc>)
///Lines are stored from down to top, and blocs are stored whitin each line
///from right to left
template <typename Element, typename Index = uint32>
class SparseBlocMatrix {

public:
	typedef LELA::SparseVector<Element, std::vector<Index>, std::vector<Element> > Vector;
	typedef Bloc<Vector> SparseBloc;
	typedef std::vector<SparseBloc> RowOfBlocs;
	typedef const RowOfBlocs ConstRowOfBlocs;
	typedef std::vector<RowOfBlocs> Rep;

	SparseBlocMatrix () : _m (0), _n (0) {}
	SparseBlocMatrix (size_t n, size_t m) :
								_m (m),
								_n (n),
								_bloc_height (128),
								_bloc_width (128)
	{
		if(_bloc_height<1 || _bloc_width < 1)
			throw WrongParameter ();

		if(n%_bloc_height == 0)						// (n/_bloc_height) lines of blocs
			_A = Rep (n/_bloc_height);
		else
			_A = Rep (n/_bloc_height + 1);
	}

	SparseBlocMatrix (size_t n, size_t m, uint16 bloc_height, uint16 bloc_width) :
								_m (m),
								_n (n),
								_bloc_height (bloc_height),
								_bloc_width (bloc_width)
	{
		if(_bloc_height<1 || _bloc_width < 1)
			throw WrongParameter ();

		if(n%_bloc_height == 0)						// (n/_bloc_height) lines of blocs
			_A (n/_bloc_height);
		else
			_A (n/_bloc_height + 1);
	}

	~SparseBlocMatrix () {}

	size_t rowdim () const { return _n; }
	size_t coldim () const { return _m; }

	uint16 bloc_rowdim () const { return _bloc_height; }
	uint16 bloc_coldim () const { return _bloc_width; }

	typedef typename Rep::iterator RowIterator;
	typedef typename Rep::const_iterator ConstRowIterator;

	RowIterator      rowBegin ()			{ return _A.begin (); }
	ConstRowIterator rowBegin ()	const	{ return _A.begin (); }
	ConstRowIterator rowEnd ()		const   { return _A.end (); }
	RowIterator      rowEnd ()				{ return _A.end (); }

	RowOfBlocs			&getRow			(size_t i) { return _A[i]; }
	RowOfBlocs			&operator []	(size_t i) { return _A[i]; }
	ConstRowOfBlocs	&operator []	(size_t i) const { return _A[i]; }
	
private:
	Rep		_A;
	size_t		_m;
	size_t		_n;

	uint16		_bloc_height;
	uint16		_bloc_width;
};


///A HybridBlocMatrix is represented as a vector of vectors of hybrid blocs
///The blocs are arranged vertically in columns (vector<HybridBloc>)
///Columns are stored from left to right, and blocs within the columns are
///stored from down to top
template <typename Element, typename Index = uint32>
class HybridBlocMatrix {
	
public:
	typedef HybridVector<Element, Index> Vector;
	typedef Bloc<Vector> SparseBloc;
	typedef std::vector<SparseBloc> ColumnOfBlocs;
	typedef const ColumnOfBlocs ConstColumnOfBlocs;
	typedef std::vector<ColumnOfBlocs> Rep;


	HybridBlocMatrix () : _m (0), _n (0) {}
	HybridBlocMatrix (size_t n, size_t m) :
								_m (m),
								_n (n),
								_bloc_height (128),
								_bloc_width (128)
	{
		if(_bloc_height<1 || _bloc_width < 1)
			throw WrongParameter ();

		if(n%_bloc_height == 0)					// (n/_bloc_height) lines of blocs
			_A = Rep (n/_bloc_height);
		else
			_A = Rep (n/_bloc_height + 1);
	}

	HybridBlocMatrix (size_t n, size_t m, uint16 bloc_height, uint16 bloc_width) :
								_A (n),
								_m (m),
								_n (n),
								_bloc_height (bloc_height),
								_bloc_width (bloc_width)
	{
		if(_bloc_height<1 || _bloc_width < 1)
			throw WrongParameter ();

		if(n%_bloc_height == 0)					// (n/_bloc_height) lines of blocs
			_A = Rep (n/_bloc_height);
		else
			_A = Rep (n/_bloc_height + 1);
	}

	~HybridBlocMatrix () {}

	size_t rowdim () const { return _n; }
	size_t coldim () const { return _m; }

	uint16 bloc_rowdim () const { return _bloc_height; }
	uint16 bloc_coldim () const { return _bloc_width; }

	typedef typename Rep::iterator RowIterator;
	typedef typename Rep::const_iterator ConstRowIterator;

	RowIterator      rowBegin ()			{ return _A.begin (); }
	ConstRowIterator rowBegin ()	const	{ return _A.begin (); }
	ConstRowIterator rowEnd ()		const   { return _A.end (); }
	RowIterator      rowEnd ()				{ return _A.end (); }

	ColumnOfBlocs			&getRow			(size_t i) { return _A[i]; }
	ColumnOfBlocs			&operator []	(size_t i) { return _A[i]; }
	ConstColumnOfBlocs		&operator []	(size_t i) const { return _A[i]; }

private:
	Rep		_A;
	size_t		_m;
	size_t		_n;

	uint16		_bloc_height;
	uint16		_bloc_width;
};



















}


#endif	//__BLOC_TYPES_H
