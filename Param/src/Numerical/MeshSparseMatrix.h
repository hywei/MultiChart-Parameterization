// MeshSparseMatrix.h: interface for the CMeshSparseMatrix class.
//
// This class is the implementation of a row-dominant sparse matrix.
// Primarily, the class is designed for sparse mesh structure.
//
//////////////////////////////////////////////////////////////////////

#ifndef MESHSPARSEMATRIX_2_H
#define MESHSPARSEMATRIX_2_H

#include "../Common/BasicDataType.h"

#include <vector>
#include <algorithm>

class CMeshSparseMatrix
{
public:
    int m_nRow;         // Number of rows
    int m_nCol;         // Number of columns
    PolyIndexArray m_ColIndex;  // The non-zero element column indices for each col
    PolyDataArray  m_ColData;   // The corresponding value of non-zero elements for each row
    PolyIndexArray m_RowIndex;  // The non-zero element column indices for each row
    PolyDataArray  m_RowData;   // The corresponding value of non-zero elements for each row
	
    double m_Epsilon;   // Zero epsilon
	bool m_IsIndexSeqed;

public:
	CMeshSparseMatrix();
	virtual ~CMeshSparseMatrix();
	void ClearData();

public:
	bool LoadSparseMatrix(std::string filename);
	bool SaveSparseMatrix(std::string filename);

public:
	void MakeCopy(CMeshSparseMatrix& other);
	void SetRowCol(int nRow, int nCol);     // Initialize, and allocate memory

	void GetElement(int row, int col, double& d);
	void SetElement(int row, int col, double d);
	void AddElement(int row, int col, double d);
	bool RemoveElement(int i, int j);

	inline int GetColNum() { return m_nCol; }
	inline int GetRowNum() { return m_nRow; }

    bool IsValidRowAndCol(int i, int j) { return (i >= 0 && i < m_nRow && j >= 0 && j < m_nCol); }
    bool IsValidRowIndex(int i) { return (i >= 0 && i < m_nRow); }
    bool IsValidColIndex(int col) {return (col >=0 && col < m_nCol);};

	bool Find(int i, int j);

    // Methods    
	bool CheckSymmetric();
	bool CheckZeroElement();
	void RemoveZeroElement();
    int  GetNZNum();

	void Transpose(CMeshSparseMatrix& mtxAT);
	void ATA(CMeshSparseMatrix& mtxATA);
    void Multiply(CMeshSparseMatrix& mtxOther, CMeshSparseMatrix& result, bool fill_row = true);  // matrix * matrix
	void Multiply(CMeshSparseMatrix& leftmtx, CMeshSparseMatrix& rightmtx, CMeshSparseMatrix& result);  // matrix * matrix
    bool MultiplyVector(DoubleArray& vec, DoubleArray& result, bool lefMulti = false);             // matrix * vector
	void Add(CMeshSparseMatrix& mtxOther);       // matrix + matrix
	void MultiplyDouble(double weight);                // matrix * double
	void MultiplyDouble4Row(std::vector<double> weight);       // matrix * vector<double>
	void InsertMatrix(int LURow, int LUCol, CMeshSparseMatrix& inmtx);
 
	//
	void GetMatrixLesseqCol(int col, std::vector< std::pair<int, double> >& colInfo);
	void MakeMatrixIndexLessSeq();
	double CalFrobeniusNorm();

	void MatrixLUTri2Full();
	double MaxAii();
	double AvgAii();

private:
	void AddNewElement(int i, int j, double d);
	void UpdateElement(int i, int j, double d, bool is_add=false);
    bool DeleteElement(int i, int j);
	void fill_row_info();
};

#endif // !defined(AFX_MESHSPARSEMATRIX_H__E867A797_11BF_4F9B_8709_6B43E5E4D8A7__INCLUDED_)
