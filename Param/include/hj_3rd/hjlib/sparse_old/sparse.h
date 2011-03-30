#ifndef HJ_SPARSE_H_
#define HJ_SPARSE_H_

#include "config.h"

#include <hj_3rd/zjucad/matrix/matrix.h>

#include <vector>
#include <map>

namespace hj { namespace sparse {

#ifdef __GNUG__
#define DEPRECATED __attribute__((deprecated))
#endif

#ifdef _MSC_VER
#define DEPRECATED __declspec(deprecated)
#endif

template <typename T>
class spm_csc	// sparse matrix csc format
{
public:
	spm_csc(){}
	spm_csc(int rows, int cols, int nnz = 0) {
		resize(rows, cols, nnz);
	}
	void resize(int rows, int cols, int nnz = 0) {
		assert(rows >= 0 && cols >= 0 && nnz >= 0);
		rows_ = rows;
		ptr_.resize(cols+1);
		ptr_[0] = 0;
		idx_.resize(nnz);
		val_.resize(nnz);
	}
	int size(int dim) const {
		return (dim == 1)?rows_:(ptr_.size()-1);
	}
	zjucad::matrix::matrix<int> ptr_, idx_;
	zjucad::matrix::matrix<T> val_;
	int rows_;
};

template <typename T>
inline int nnz(const spm_csc<T> &A) {
	return A.idx_.size();
}

/// convert
template <typename E, typename T>
spm_csc<T> &convert(const zjucad::matrix::matrix_expression<E> &A, spm_csc<T> &B, T eps = 0)
{
	B.resize(A().size(1), A().size(2), zjucad::matrix::sum<int>(fabs(A) >= eps));
	int ci = 0, ri = 0, pi = 0, vi = 0, ai = 0;
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

template <typename E, typename T>
zjucad::matrix::matrix_expression<E> &
convert(const spm_csc<T> &A, zjucad::matrix::matrix_expression<E> &B)
{
	B() = zeros(A.size(1), A.size(2));
	int ci = 0, ri = 0, pi = 0, vi = 0, ai = 0;
	for(; ci < A.size(2); ++ci) {	// for each column
		B()(zjucad::matrix::colon(), ci)
			(A.idx_(zjucad::matrix::colon(A.ptr_[ci], A.ptr_[ci+1]-1)))
			= A.val_(zjucad::matrix::colon(A.ptr_[ci], A.ptr_[ci+1]-1));
	}
	return B;
}

template <typename T1, typename T2>
spm_csc<T2> &convert(const std::vector<std::map<int, T1> > &vmM, spm_csc<T2> &cscM, int rows)
{
	cscM.rows_ = rows;
	const int cols = static_cast<int>(vmM.size());
	cscM.ptr_.resize(cols+1);
	cscM.ptr_[0] = 0;
	int i, row_i;
	for(i = 0; i < cols; ++i) // for each column
		cscM.ptr_[i+1] = cscM.ptr_[i]+static_cast<int>(vmM[i].size());
	cscM.idx_.resize(cscM.ptr_[cols]);
	cscM.val_.resize(cscM.ptr_[cols]);

	std::map<int, double>::const_iterator mi;
	for(i = 0; i < cols; ++i) {	// for each column
		const std::map<int, double> &col = vmM[i];
		for(row_i = cscM.ptr_[i], mi = col.begin(); row_i < cscM.ptr_[i+1]; ++row_i, ++mi) {
			cscM.idx_[row_i] = mi->first;
			cscM.val_[row_i] = mi->second;
		}
	}
	return cscM;
}

/// operation

// NOTICE: C = A+B
template <typename T1, typename T2, typename T3>
spm_csc<T3> &spm_mpm(const spm_csc<T1> &A, const spm_csc<T2> &B, spm_csc<T3> &C)
{
	assert(A.size(1) == B.size(1) && A.size(2) == B.size(2));
	const int rows = A.size(1), cols = A.size(2);
	std::vector<std::map<int, double> > mapC(cols);
	int i, j;
	for(i = 0; i < cols; ++i) {
		std::map<int, double> &col = mapC[i];
		for(j = A.ptr_[i]; j < A.ptr_[i+1]; ++j)
			col[A.idx_[j]] += A.val_[j];
		for(j = B.ptr_[i]; j < B.ptr_[i+1]; ++j)
			col[B.idx_[j]] += B.val_[j];
	}

	return convert(mapC, C, rows);
}

// NOTICE: C += A*B
template <typename T, typename M1, typename M2>
M2 &spm_mm(bool transA, const spm_csc<T> &A,
		   bool transB, const zjucad::matrix::matrix_expression<M1> &B,
		   bool transC, zjucad::matrix::matrix_expression<M2> &C)
{
	using namespace zjucad::matrix;
#define spm_mm_MULT(op) \
	for(int ci = 0; ci < A.size(2); ++ci) { for(; vi < A.ptr_[ci+1]; ++vi) { \
		op;\
	}}

	int vi = A.ptr_[0];
	if(transA) {
		if(transB) {
			if(transC) {
				spm_mm_MULT(C()(colon(), ci) += A.val_[vi]*B()(colon(), A.idx_[vi]));
			}
			else {
				spm_mm_MULT(C()(ci, colon()) += A.val_[vi]*trans(B()(colon(), A.idx_[vi])));
			}
		}
		else {
			if(transC) {
				spm_mm_MULT(C()(colon(), ci) += A.val_[vi]*trans(B()(A.idx_[vi], colon())));
			}
			else {
				spm_mm_MULT(C()(ci, colon()) += A.val_[vi]*B()(A.idx_[vi], colon()));
			}
		}
	}
	else {
		if(transB) {
			if(transC) {
				spm_mm_MULT(C()(colon(), A.idx_[vi]) += A.val_[vi]*B()(colon(), ci));
			}
			else {
				spm_mm_MULT(C()(A.idx_[vi], colon()) += A.val_[vi]*trans(B()(colon(), ci)));
			}
		}
		else {
			if(transC) {
				spm_mm_MULT(C()(colon(), A.idx_[vi]) += A.val_[vi]*trans(B()(ci, colon())));
			}
			else {
				spm_mm_MULT(C()(A.idx_[vi], colon()) += A.val_[vi]*B()(ci, colon()));
			}
		}
	}
#undef spm_mm_MULT
	return C();
}

template <typename T, typename V1, typename V2>
V2 &spm_mxv(const spm_csc<T> &A, const V1 &b, V2 &c)
{
	int vi = A.ptr_[0];
	for(int ci = 0; ci < A.size(2); ++ci) {
		for(; vi < A.ptr_[ci+1]; ++vi) {
			c[A.idx_[vi]] += A.val_[vi]*b[ci];
		}
	}
	return c;
}

template <typename T, typename V1, typename V2>
V2 &spm_mTxv(const spm_csc<T> &A, const V1 &b, V2 &c)
{
	int vi = A.ptr_[0];
	for(int ci = 0; ci < A.size(2); ++ci) {
		for(; vi < A.ptr_[ci+1]; ++vi) {
			c[ci] += A.val_[vi]*b[A.idx_[vi]];
		}
	}
	return c;
}

template <typename T1, typename T2>
hj::sparse::spm_csc<T2> &AAT(const hj::sparse::spm_csc<T1> &A, hj::sparse::spm_csc<T2> &AAT)
{
	const int m = A.size(1), n = A.size(2);
	std::vector<std::map<int, T2> > mapAAT(m);	// [column][row]
	int i, col_i = 0, row_i = 0;
	for(i = 0; i < n; ++i) {	// vector outer product: vvT(i, j) = vi*vj, for each v
		for(col_i = A.ptr_[i]; col_i < A.ptr_[i+1]; ++col_i) {	// for each nz column
			if(A.val_[col_i] == 0) continue;
			std::map<int, T2> &col = mapAAT[A.idx_[col_i]];
			for(row_i = A.ptr_[i]; row_i < A.ptr_[i+1]; ++row_i)	// scalar * sparse, for each nz row
				col[A.idx_[row_i]] += A.val_[col_i]*A.val_[row_i];
		}
	}

	return convert(mapAAT, AAT, m);
}

int inline get_idx_from_row_col(const zjucad::matrix::matrix<int> &ptr, const zjucad::matrix::matrix<int> &row_idx,
						 int row, int col)
{
	int rtn = -1;
	for(rtn = ptr[col]; rtn < ptr[col+1]; ++rtn)
		if(row_idx[rtn] == row) break;
	if(rtn == ptr[col+1])
		return -1;
	return rtn;
}

// assume that AAT have the same patten
template <typename T1, typename T2>
hj::sparse::spm_csc<T2> &fast_AAT(const hj::sparse::spm_csc<T1> &A, hj::sparse::spm_csc<T2> &AAT)
{
	const int m = A.size(1), n = A.size(2);
	fill(AAT.val_.begin(), AAT.val_.end(), 0);
	int i, col_i = 0, row_i = 0;
	for(i = 0; i < n; ++i) {	// vector outer product: vvT(i, j) = vi*vj, for each v
		for(col_i = A.ptr_[i]; col_i < A.ptr_[i+1]; ++col_i) {	// for each nz column
			if(A.val_[col_i] == 0) continue;
			for(row_i = A.ptr_[i]; row_i < A.ptr_[i+1]; ++row_i) {	// scalar * sparse, for each nz row
				//mapAAT[A.idx_[col_i]][A.idx_[row_i]] += A.val_[col_i]*A.val_[row_i];
				int pos = get_idx_from_row_col(AAT.ptr_, AAT.idx_, A.idx_[row_i], A.idx_[col_i]);
				assert(pos >= 0);
				AAT.val_[pos] += A.val_[col_i]*A.val_[row_i];
			}
		}
	}

	return AAT;
}

class HJ_SPARSE_API solver
{
public:
	virtual ~solver();
	/**
	   @param id "umfpack" or "cholmod", default 0 is current best
	   implementation
	*/
	DEPRECATED static solver *create(const spm_csc<double> &A, const char *id = 0);
	virtual bool set_patten(int rows, const zjucad::matrix::matrix<int> &ptr, const zjucad::matrix::matrix<int> &idx) = 0;
	virtual bool set_value (const zjucad::matrix::matrix<double> &val) = 0;
	virtual bool solve(const double *b, double *x, int nrhs = 1) = 0;
protected:
	virtual bool init(const hj::sparse::spm_csc<double> &A) = 0;
};

#undef DEPRECATED

}}

#endif
