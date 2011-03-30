#ifndef HJ_SPARSE_CONVERT_H_
#define HJ_SPARSE_CONVERT_H_

namespace hj { namespace sparse {

template <typename E, typename T, typename INT_TYPE>
csc<T, INT_TYPE> &convert(const zjucad::matrix::matrix_expression<E> &A, csc<T, INT_TYPE> &B, T eps = 0)
{
	using namespace zjucad::matrix;
	B.resize(A().size(1), A().size(2), zjucad::matrix::sum<size_type>(fabs(A) >= eps));
	INT_TYPE ci = 0, ri = 0, vi = 0, ai = 0;
	for(; ci < A().size(2); ++ci) {	// for each column
		B.ptr_[ci+1] = B.ptr_[ci];
		for(ri = 0; ri < A().size(1); ++ri, ++ai) {
			if(fabs(A()[ai]) < eps) continue;
			++B.ptr_[ci+1];
			B.idx_[vi] = ri;
			B.val_[vi]= A()[ai];
			++vi;
		}
	}
	assert(vi == nnz(B));
	return B;
}

template <typename E, typename T, typename INT_TYPE>
zjucad::matrix::matrix_expression<E> &
convert(const csc<T, INT_TYPE> &A, zjucad::matrix::matrix_expression<E> &B)
{
	using namespace zjucad::matrix;
	B() = zeros(A.size(1), A.size(2));
	INT_TYPE ci = 0;
	for(; ci < A.size(2); ++ci) {	// for each column
		B()(zjucad::matrix::colon(), ci)
			(A.idx_(zjucad::matrix::colon(A.ptr_[ci], A.ptr_[ci+1]-1)))
			= A.val_(zjucad::matrix::colon(A.ptr_[ci], A.ptr_[ci+1]-1));
	}
	return B;
}

template <typename T1, typename INT_TYPE1,
		  template <typename K, typename T> class MAP_TYPE,
		  typename T2, typename INT_TYPE2>
csc<T2, INT_TYPE2> &convert(const csc_by_vm<T1, INT_TYPE1, MAP_TYPE> &vmM,
							csc<T2, INT_TYPE2> &cscM)
{
	cscM.rows_ = vmM.size(1);
	const INT_TYPE2 cols = static_cast<INT_TYPE2>(vmM.size(2));
	cscM.ptr_.resize(cols+1);
	cscM.ptr_[0] = 0;
	INT_TYPE2 i, row_i;
	for(i = 0; i < cols; ++i) // for each column
		cscM.ptr_[i+1] = cscM.ptr_[i]+static_cast<INT_TYPE2>(nnz(vmM[i]));
	cscM.idx_.resize(cscM.ptr_[cols]);
	cscM.val_.resize(cscM.ptr_[cols]);

	typename map_vec<T1, INT_TYPE1, MAP_TYPE>::nz_const_iterator mi;
	for(i = 0; i < cols; ++i) {	// for each column
		const map_vec<T1, INT_TYPE1, MAP_TYPE> &col = vmM[i];
		for(row_i = cscM.ptr_[i], mi = col.begin_nz(); row_i < cscM.ptr_[i+1]; ++row_i, ++mi) {
			cscM.idx_[row_i] = mi->first;
			cscM.val_[row_i] = mi->second;
		}
	}
	return cscM;
}

}}

#endif
