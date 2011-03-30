#ifndef MATRIX_CONVERTER_2_H
#define MATRIX_CONVERTER_2_H


#include "MeshSparseMatrix.h"
#ifdef WIN32
#include <hj_3rd/hjlib/sparse_old/sparse.h>
#else
#include <hj_3rd/hjlib/sparse/sparse.h>
#endif
using namespace zjucad::matrix;

class CMatrixConverter
{
public:
	CMatrixConverter();
	~CMatrixConverter();

public:
    static void CSparseMatrix2hjCscMatrix(hj::sparse::csc<double>& hjcscMatrix, CMeshSparseMatrix& cMatrix);
	static void CSparseMatrix2hjCscMatrix(hj::sparse::csc<double,int>& hjcscMatrix, CMeshSparseMatrix& cMatrix);
};

#endif
