#include "MatrixConverter.h"

//////////////////////////////////////////////////////////////////////
// CMatrixConverter Construction/Destruction
//////////////////////////////////////////////////////////////////////
CMatrixConverter::CMatrixConverter()
{
}
CMatrixConverter::~CMatrixConverter()
{
}

//////////////////////////////////////////////////////////////////////
// CMatrixConverter public methods
//////////////////////////////////////////////////////////////////////


void CMatrixConverter::CSparseMatrix2hjCscMatrix(hj::sparse::csc<double>& hjcscMatrix, CMeshSparseMatrix& cMatrix)
{
	//
	cMatrix.MakeMatrixIndexLessSeq();

	int i;
	int m_Row = cMatrix.GetRowNum();
	int m_Col = cMatrix.GetColNum();

	int nnz = cMatrix.GetNZNum();
	hjcscMatrix.resize(m_Row, m_Col, nnz);

	//
	int index = 0;
	double value = 0;

	for(i = 0; i < m_Col; ++ i)
	{
		(hjcscMatrix.ptr_)[i] = index;

		size_t n = cMatrix.m_ColIndex[i].size();

		for(size_t j = 0; j < n; ++ j)
		{
			(hjcscMatrix.idx_)[index] = cMatrix.m_ColIndex[i][j];
			(hjcscMatrix.val_)[index ++] = cMatrix.m_ColData[i][j];
		}
	}
	hjcscMatrix.ptr_[i] = nnz;
}

void CMatrixConverter::CSparseMatrix2hjCscMatrix(hj::sparse::csc<double,int>& hjcscMatrix, CMeshSparseMatrix& cMatrix)
{
	int i;
	int m_Row = cMatrix.GetRowNum();
	int m_Col = cMatrix.GetColNum();

	int nnz = cMatrix.GetNZNum();
	hjcscMatrix.resize(m_Row, m_Col, nnz);

	//
	int index = 0;
	double value = 0;

	for(i = 0; i < m_Col; ++ i)
	{
		(hjcscMatrix.ptr_)[i] = index;

		//size_t n = cMatrix.m_ColIndex[i].size();

		std::vector<std::pair<int, double> > col_info;
		cMatrix.GetMatrixLesseqCol(i, col_info);
		for(size_t j = 0; j < col_info.size(); ++ j)
		{
			(hjcscMatrix.idx_)[index] = col_info[j].first;
			(hjcscMatrix.val_)[index ++] = col_info[j].second;
		}

		/*
		for(size_t j = 0; j < n; ++ j)
		{
		(hjcscMatrix.idx_)[index] = cMatrix.m_ColIndex[i][j];
		(hjcscMatrix.val_)[index ++] = cMatrix.m_ColData[i][j];
		}*/
	}
	hjcscMatrix.ptr_[i] = nnz;

}
