#ifndef HJ_SPARSE2_H_
#define HJ_SPARSE2_H_

#include <cassert>


#include <zjucad/matrix/matrix.h>

namespace hj { namespace sparse2 {

//! typical column of a csc matrix is map, paired/unpaired +
//! sorted/unsorted array.
template <typename VAL_TYPE, typename INT_TYPE>
class csc
{
public:
	size_t size(void) const { return nrows_*(ptr_.size()-1); }
	size_t size(int dim) const { return (dim == 1)?nrows_:(ptr_.size()-1); }
	size_t nnz(void) const {
		return ptr_[ptr_.size()-1]-ptr_[0];
	}

	class col_vec {
	public:
		col_vec(csc *th, INT_TYPE col)
			: this_(th), col_(col) {
		}
		class nz_iterator {
		public:
			nz_iterator(csc *th, size_t nzi)
				:this_(th), nzi_(nzi) {
			}
			nz_iterator &operator ++(void) {
				++nzi_;
				return *this;
			}
			template <typename I1>
			bool operator != (const I1 &i) {
				assert(this_ == i.this_);
				return nzi_ != i.nzi_;
			}
			INT_TYPE idx(void) const { return this_->idx_[nzi_]; }
			VAL_TYPE &operator *(void) const { return this_->val_[nzi_]; }

		private:
			csc *this_;
			size_t nzi_;
		};

		nz_iterator begin_nz(void) { return nz_iterator(this_, this_->ptr_[col_]);}
		nz_iterator end_nz(void) { return nz_iterator(this_, this_->ptr_[col_+1]); }

		/*
		  insert(INT_TYPE idx, VAL_TYPE val);
		  remove(INT_TYPE idx, VAL_TYPE val);
		*/
		nz_iterator lower_bound(INT_TYPE idx) {// no sort
			nz_iterator i;
			for(i = begin_nz(); i != end_nz(); ++i)
				if(i.idx() == idx)
					return i;
			return end_nz();
		}
		// nz_iterator lower_bound(INT_TYPE idx) {// sort
		// 	INT_TYPE *pos = std::lower_bound(
		// 		&this_->idx_[0]+this_->ptr_[col],
		// 		&this_->idx_[0]+this_->ptr_[col+1],
		// 		idx);
		// 	if(pos != &this_->idx_[0]+this_->ptr_[col+1])
		// 		return nz_iterator(this_, pos);
		// 	return end();
		// }

		csc *this_;
		INT_TYPE col_;
	};

	class col_iterator {
	public:
		col_iterator(){}
		col_iterator(csc *th, size_t col)
			:cv_(th, col) {
		}
		col_vec &operator*(void) {
			return cv_;
		}
		col_vec *operator->(void) {
			return &cv_;
		}
		col_iterator operator ++(void) {
			++cv_.col_;
			return *this;
		}
		template <typename I1>
		bool operator != (const I1 &i) const {
			assert(this_ == i.this_);
			return cv_.col_ != i.cv_.col_;
		}
		col_vec cv_;
	};

	col_iterator begin_col(void) {
		return col_iterator(this, 0);
	}
	col_iterator end_col(void) {
		return col_iterator(this, ptr_.size()-1);
	}

	zjucad::matrix::matrix<INT_TYPE> ptr_, idx_;
	zjucad::matrix::matrix<VAL_TYPE> val_;
	INT_TYPE nrows_;
};

template <typename VAL_TYPE, typename INT_TYPE, typename COL_VEC_TYPE>
class vov
{
public:
};


}}

#endif
