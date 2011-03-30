#ifndef HJ_SPARSE_FORMAT_H_
#define HJ_SPARSE_FORMAT_H_

#include <zjucad/matrix/matrix.h>

#include <vector>
#include <map>

namespace hj { namespace sparse {

typedef zjucad::matrix::size_type size_type;
typedef zjucad::matrix::idx_type idx_type;
typedef zjucad::matrix::offset_type offset_type;

template <typename K, typename T>
class std_map : public std::map<K, T>
{
public:
	void reserve(size_t){}
};

template <typename K, typename T>
class map_by_vector
{
public:
	typedef std::vector<std::pair<K, T> > container;
	typedef typename container::iterator iterator;
	typedef typename container::const_iterator const_iterator;

	typedef std::pair<const K, T> value_type;

	void reserve(size_t nnz) {
		data_.reserve(nnz);
	}
	size_t size(void) const {return data_.size();}

	const_iterator begin(void) const { return data_.begin(); }
	const_iterator end(void) const { return data_.end(); }
	iterator begin(void) { return data_.begin(); }
	iterator end(void) { return data_.end(); }

	void clear(void) { data_.clear(); }

	container data_;
};

template <typename K, typename T>
class map_by_sorted_vector : public map_by_vector<K, T>
{
public:
	typedef typename map_by_vector<K, T>::iterator iterator;
	static inline bool less(const std::pair<K, T> &a, const std::pair<K, T> &b) {
		return a.first < b.first;
	}
	iterator lower_bound(K i) {
		return std::lower_bound(this->data_.begin(), this->data_.end(), std::make_pair(i, T()), less);
	}
	iterator insert(iterator itr, const std::pair<K, T> &v) {
		return this->data_.insert(itr, v);
	}
};


template <typename K, typename T>
class map_by_unsorted_vector : public map_by_sorted_vector<K, T>
{
public:
	typedef typename map_by_vector<K, T>::iterator iterator;

	iterator lower_bound(const K &i) {
		iterator rtn;
		for(rtn = this->data_.begin(); rtn != this->data_.end(); ++rtn)
			if(rtn->first == i)
				break;
		return rtn;
	}
	iterator insert(const iterator &itr, const std::pair<K, T> &v) {
		if(itr == this->data_.end()) {
			this->data_.push_back(v);
			return this->data_.end()-1;
		}
		return itr;
	}
};


//! @brief map_vec
template <typename VAL_TYPE = double,
		  typename INT_TYPE = idx_type,
		  template <typename K, typename T> class MAP_TYPE = std_map>
class map_vec
{
public:
	typedef VAL_TYPE value_type;
	typedef VAL_TYPE &reference;
	typedef INT_TYPE idx_type;
	typedef const VAL_TYPE & const_reference;
	typedef MAP_TYPE<INT_TYPE, VAL_TYPE> container;

	typedef typename container::const_iterator nz_const_iterator;
	typedef typename container::iterator nz_iterator;

	map_vec(void):len_(0){}
	map_vec(idx_type len):len_(len){}
	map_vec(idx_type len, idx_type nnz):len_(len) {
		reserve(nnz);
	}

	void reserve(idx_type nnz) {
		data_.reserve(nnz);
	}

	idx_type size(void) const {return len_;}

	//! @brief change the size and clear the data, i.e. make nnz to be zero
	void resize(idx_type len) { len_ = len; data_.clear(); }
	void resize(idx_type len, idx_type nnz) { resize(len_); data_.reserve(nnz); }

	inline const_reference operator[](idx_type i) const {
		return get(i);
	}

	class assign_proxy {
	public:
		assign_proxy(map_vec<value_type, idx_type, MAP_TYPE> &th, idx_type i)
			:this_(&th), i_(i)
			{}
		inline operator const_reference() const {
			return this_->get(i_);
		}
		template <typename T>
		void inline operator = (const T &v) {
			this_->get_or_create(i_) = v;
		}
		template <typename T>
		void inline operator += (const T &v) {
			this_->get_or_create(i_) += v;
		}
		template <typename T>
		void inline operator -= (const T &v) {
			this_->get_or_create(i_) -= v;
		}
		template <typename T>
		void inline operator *= (const T &v) {
			this_->get_or_create(i_) *= v;
		}
		template <typename T>
		void inline operator /= (const T &v) {
			this_->get_or_create(i_) *= v;
			return *this;
		}
	private:
		map_vec<value_type, idx_type, MAP_TYPE> * const this_;
		const idx_type i_;
	};

	inline assign_proxy operator[](idx_type i) {
		return assign_proxy(*this, i);
	}

	inline const_reference get(idx_type i) const {
		typename container::const_iterator pos = data_.find(i);
		if(pos != data_.end())
			return pos->second;
		return zero_;
	}
	reference get_or_create(idx_type i) {
		typename container::iterator pos = data_.lower_bound(i);
		if(pos == data_.end() || pos->first != i)
			pos = data_.insert(pos, std::make_pair(i, 0));
		return pos->second;
	}

	nz_const_iterator begin_nz(void) const { return data_.begin(); }
	nz_const_iterator end_nz(void) const { return data_.end(); }
	nz_iterator begin_nz(void) { return data_.begin(); }
	nz_iterator end_nz(void) { return data_.end(); }

	container data_;
	idx_type len_;
	static const value_type zero_;
};

template <typename VAL_TYPE, typename INT_TYPE,
		  template <typename K, typename T> class MAP_TYPE>
const typename map_vec<VAL_TYPE, INT_TYPE, MAP_TYPE>::value_type
map_vec<VAL_TYPE, INT_TYPE, MAP_TYPE>::zero_ = value_type();

template <typename VAL_TYPE, typename INT_TYPE,
		  template <typename K, typename T> class MAP_TYPE>
inline size_t nnz(const map_vec<VAL_TYPE, INT_TYPE, MAP_TYPE> &v) {
	return v.data_.size();
}

//! @brief csc
template <typename T = double, typename INT_TYPE = idx_type>
class csc	// sparse matrix csc format
{
public:
	csc():rows_(0){}
	csc(INT_TYPE rows, INT_TYPE cols, INT_TYPE nnz = 0) {
		resize(rows, cols, nnz);
	}
	void resize(INT_TYPE rows, INT_TYPE cols, INT_TYPE nnz = 0) {
		assert(rows >= 0 && cols >= 0 && nnz >= 0);
		rows_ = rows;
		ptr_.resize(cols+1);
		ptr_[0] = 0;
		idx_.resize(nnz);
		val_.resize(nnz);
	}
	size_type size(int dim) const {
		return (dim == 1)?rows_:(ptr_.size()-1);
	}
	zjucad::matrix::matrix<INT_TYPE> ptr_, idx_;
	zjucad::matrix::matrix<T> val_;
	INT_TYPE rows_;
};

template <typename T, typename INT_TYPE>
inline INT_TYPE nnz(const csc<T, INT_TYPE> &A) {
	return A.idx_.size();
}

//! @brief csc_by_vm: csc by vector<map_vec>

template <typename VAL_TYPE = double, typename INT_TYPE = idx_type,
		  template <typename K, typename T> class MAP_TYPE = map_by_unsorted_vector>
class csc_by_vm
{
public:
	typedef map_vec<VAL_TYPE, INT_TYPE, MAP_TYPE> map_vec_type;
	typedef std::vector<map_vec_type> container;
	
	typedef typename container::const_iterator col_const_iterator;
	typedef typename container::iterator col_iterator;

	csc_by_vm(void)
		:nrows_(0) {
	}

	csc_by_vm(INT_TYPE nrows, INT_TYPE ncols) {
		resize(nrows, ncols);
	}
	csc_by_vm(INT_TYPE nrows, INT_TYPE ncols, INT_TYPE est_nz_per_col) {
		resize(nrows, ncols, est_nz_per_col);
	}

	INT_TYPE size(void) const {
		return nrows_*data_.size();
	}
	INT_TYPE size(int dim) const {
		return (dim == 1)?nrows_:data_.size();
	}

	//! @brief change the size and clear the data, i.e. make nnz to be zero
	void resize(INT_TYPE nrow, INT_TYPE ncol) {
		nrows_ = nrow;
		data_.resize(ncol);
		for(size_t ci = 0; ci < data_.size(); ++ci)
			data_[ci].resize(nrow);
	}

	//! @ param est_nz_per_col : estimated nz per col
	void resize(INT_TYPE nrow, INT_TYPE ncol, INT_TYPE est_nz_per_col) {
		nrows_ = nrow;
		data_.resize(ncol);

		for(size_t ci = 0; ci < data_.size(); ++ci)
			data_[ci].resize(nrow, est_nz_per_col);
	}

	const map_vec_type &operator[](INT_TYPE i) const { return data_[i]; }
	map_vec_type &operator[](INT_TYPE i) { return data_[i]; }

	col_const_iterator begin_col(void) const { return data_.begin(); }
	col_const_iterator end_col(void) const { return data_.end(); }
	col_iterator begin_col(void)  { return data_.begin(); }
	col_iterator end_col(void)  { return data_.end(); }

	INT_TYPE nrows_;
	container data_;
};

template <typename VAL_TYPE, typename INT_TYPE,
		  template <typename K, typename T> class MAP_TYPE>
INT_TYPE nnz(const csc_by_vm<VAL_TYPE, INT_TYPE, MAP_TYPE> &A) {
	INT_TYPE nz = 0;
	for(size_t i = 0; i < A.data_.size(); ++i)
		nz += nnz(A.data_[i]);
	return nz;
}

}}

#endif
