// MeshSparseMatrix.cpp: implementation of the CMeshSparseMatrix class.
//
//////////////////////////////////////////////////////////////////////

#include "MeshSparseMatrix.h"
#include "../Common/Utility.h"
#include <assert.h>
#include <math.h>
#include <iostream>
#include <fstream>
using namespace std;
#pragma warning(disable: 4018 4267)

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
struct DataCompare
{
	bool operator() (const pair<int, double>& lhs, const pair<int, double>& rhs) const 
	{
		return lhs.first < rhs.first;
	}
};

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CMeshSparseMatrix::CMeshSparseMatrix()
{
    ClearData();
}

CMeshSparseMatrix::~CMeshSparseMatrix()
{
    ClearData();
}

void CMeshSparseMatrix::ClearData()
{
	m_nRow = 0;
	m_nCol = 0;
	m_Epsilon = LARGE_ZERO_EPSILON;	// Zero epsilon, value less than this threshold is considered to be zero
	m_IsIndexSeqed = false;
	
    Utility util;
    util.FreeVector(m_ColIndex);
    util.FreeVector(m_ColData);
	util.FreeVector(m_RowIndex);
    util.FreeVector(m_RowData);
}
void CMeshSparseMatrix::MakeCopy(CMeshSparseMatrix& other)
{
	other.m_IsIndexSeqed = this->m_IsIndexSeqed;
	other.m_nRow = this->m_nRow;
	other.m_nCol = this->m_nCol;
	other.m_Epsilon = this->m_Epsilon;
	other.m_RowData.assign(this->m_RowData.begin(), this->m_RowData.end());
	other.m_RowIndex.assign(this->m_RowIndex.begin(), this->m_RowIndex.end());
	other.m_ColData.assign(this->m_ColData.begin(), this->m_ColData.end());
	other.m_ColIndex.assign(this->m_ColIndex.begin(), this->m_ColIndex.end());
}
// Initialize, and allocate memory
void CMeshSparseMatrix::SetRowCol(int nRow, int nCol)
{
	ClearData();

    assert(nRow >= 0 && nCol >= 0);

	m_nRow = nRow;
	m_nCol = nCol;
	m_ColIndex.resize(nCol);
	m_ColData.resize(nCol);
	m_RowIndex.resize(m_nRow);
	m_RowData.resize(m_nRow);

	int i;
	for(i = 0; i < nCol; ++ i)
	{
        m_ColIndex[i].resize(0);

        m_ColData[i].resize(0);
	}

	for(i = 0; i < m_nRow; ++ i)
	{
        m_RowIndex[i].resize(0);

        m_RowData[i].resize(0);
	}
}

void CMeshSparseMatrix::GetElement(int row, int col, double& d)
{
    assert(IsValidRowAndCol(row, col));

	d = 0.0;
    size_t k, n = m_RowIndex[row].size();

	for(k = 0; k < n; ++ k)
	{
		if(m_RowIndex[row][k] == col)
		{
			d = m_RowData[row][k];
			return;
		}
	}

	k, n = m_ColIndex[col].size();

	for(k = 0; k < n; ++ k)
	{
		if(m_ColIndex[col][k] == row)
		{
			d = m_ColData[col][k];
			return;
		}
	}
}

// Return true:	Add a new element
// Return false: Update an existed element
void CMeshSparseMatrix::SetElement(int row, int col, double d)
{
	assert(IsValidRowAndCol(row, col));

	UpdateElement(row, col, d);
	
	/*
	if(Find(row, col))
	{
		UpdateElement(row, col, d);
		return false;
	}
	else if(IsValidColIndex(col))
	{
		AddNewElement(row, col, d);
		return true;
	}
	else
	{
		printf("Error: Invalid Row and Col Index when setting.");
		assert(false);
		return false;
	}*/
}
void CMeshSparseMatrix::AddElement(int row, int col, double d)
{
	assert(IsValidRowAndCol(row, col));

	UpdateElement(row, col, d, true);
}
bool CMeshSparseMatrix::RemoveElement(int row, int col)
{
	if(Find(row, col))
	{
		// remove this element from matrix here.
		DeleteElement(row, col);

		return true;
	}
	else
	{
		printf("Error: Invalid Row and Col Index when deleting.");
		assert(false);
		return false;
	}
}
// Return true: Find Element(i,j)
// Return false: Not Found
bool CMeshSparseMatrix::Find(int row, int col)
{
	if(!IsValidColIndex(col))
		return false;

	size_t nElement = m_ColIndex[col].size();

	for(size_t k = 0; k < nElement; ++ k)
	{
		if(m_ColIndex[col][k] == row)
			return true;
	}

	return false;
}

void CMeshSparseMatrix::AddNewElement(int row, int col, double d)
{
	m_ColIndex[col].push_back(row);
	m_ColData[col].push_back(d);
	m_RowIndex[row].push_back(col);
	m_RowData[row].push_back(d);
}

void CMeshSparseMatrix::UpdateElement(int row, int col, double d, bool is_add)
{
	size_t nElement = m_ColIndex[col].size();
    size_t k;
	bool is_found = false;
	for(k = 0; k < nElement; ++ k)
	{
		if(m_ColIndex[col][k] == row)
		{
			if (is_add)
			{
				m_ColData[col][k] += d;
			}
			else
			{
				m_ColData[col][k] = d;
			}

			is_found = true;
			break;
		}
	}
	if (!is_found)
	{
		m_ColIndex[col].push_back(row);
		m_ColData[col].push_back(d);
	}

	nElement = m_RowIndex[row].size();
	is_found = false;
	for(k = 0; k < nElement; ++ k)
	{
		if(m_RowIndex[row][k] == col)
		{
			if (is_add)
			{
				m_RowData[row][k] += d;
			}
			else
			{
				m_RowData[row][k] = d;
			}
			
			is_found = true;
			break;
		}
	}
	if (!is_found)
	{
		m_RowIndex[row].push_back(col);
		m_RowData[row].push_back(d);
	}
}
bool CMeshSparseMatrix::DeleteElement(int row, int col)
{
	size_t nElement = m_ColIndex[col].size();
    size_t k;

	// delete the element form column index and data.
	for(k = 0; k < nElement; ++ k)
	{
		if(m_ColIndex[col][k] == row)
		{
			// remove the row index from the column index array.
			IndexArray& colIndex = m_ColIndex[col];
			colIndex.erase(remove(colIndex.begin(), colIndex.end(), row), colIndex.end());

			// remove the row data from the column data array.
			DoubleArray data;
			data.resize(nElement - 1);
			int index = 0;

			for (size_t i = 0; i < nElement; i++)
			{
				if (i != k)
				{
					data[index++] = m_ColData[col][i];
				}
			}
			m_ColData[col].assign(data.begin(), data.end());

			break;
		}
	}

	// delete the element form row index and data.

	nElement = m_RowIndex[row].size();
	for(k = 0; k < nElement; ++ k)
	{
		if(m_RowIndex[row][k] == col)
		{
			// remove the row index from the column index array.
			IndexArray& rowIndex = m_RowIndex[col];
			rowIndex.erase(remove(rowIndex.begin(), rowIndex.end(), col), rowIndex.end());

			// remove the row data from the column data array.
			DoubleArray data;
			data.resize(nElement - 1);
			int index = 0;
			
			for (size_t i = 0; i < nElement; i++)
			{
				if (i != k)
				{
					data[index++] = m_RowData[row][i];
				}
			}
			m_RowData[row].assign(data.begin(), data.end());

			break;
		}
	}

	return true;
}
bool CMeshSparseMatrix::CheckSymmetric()
{
	int nNotSymEle = 0;
    double v1, v2;
    int i;
	size_t j, nElement;
	for(i = 0; i < m_nCol; ++ i)
	{
		nElement = m_ColIndex[i].size();
		for(j = 0; j < nElement; ++ j)
		{
			GetElement(i, m_ColIndex[i][j], v1);
			GetElement(m_ColIndex[i][j], i, v2);

			double eps = fabs(v1-v2); 

			if(eps > m_Epsilon)
			{
				++ nNotSymEle;
			}
		}
	}

	fprintf(stdout,"\n# Not Symmetric Elements = %d\n",nNotSymEle);
	
    return (nNotSymEle == 0);
}
bool CMeshSparseMatrix::CheckZeroElement()
{
	int i;
	size_t j, nElement;
	double value;
	int nZero = 0;
	
	for(i = 0; i < m_nCol; ++ i)
	{
		nElement = m_ColIndex[i].size();
		for(j = 0; j < nElement; ++ j)
		{
			GetElement(m_ColIndex[i][j], i, value);
			
			if (ALMOST_EQUAL_SMALL(value, 0.0))
			{
				nZero++;
			}
		}
	}

	fprintf(stdout,"\n# Zero Elements = %d\n", nZero);
	
    return (nZero == 0);
}
void CMeshSparseMatrix::RemoveZeroElement()
{
	int i;
	size_t j, nElement;
	double value;

	for(i = 0; i < m_nCol; ++ i)
	{
		nElement = m_ColIndex[i].size();
		for(j = 0; j < nElement; ++ j)
		{
			GetElement(m_ColIndex[i][j], i, value);

			if (ALMOST_EQUAL_SMALL(value, 0.0))
			{
				DeleteElement(m_ColIndex[i][j], i);
			}
		}
	}
}
int CMeshSparseMatrix::GetNZNum()
{ 
	int i;
	size_t nElement;
	int nzNum = 0;

	for(i = 0; i < m_nCol; ++ i)
	{
		nElement = m_ColIndex[i].size();
		nzNum += (int) nElement;
	}

	return nzNum;
}
void CMeshSparseMatrix::ATA(CMeshSparseMatrix& mtxATA)
{
	CMeshSparseMatrix mtx;
	CMeshSparseMatrix mtxAT;
	Transpose(mtxAT);
	//MakeCopy(mtx);

	mtxAT.Multiply(*this, mtxATA, false);
}
void CMeshSparseMatrix::Transpose(CMeshSparseMatrix& mtxAT)
{
	//mtxAT.SetRowCol(m_nCol,m_nRow);

	mtxAT.m_nRow = this->m_nCol;
	mtxAT.m_nCol = this->m_nRow;
	mtxAT.m_Epsilon = this->m_Epsilon;
	mtxAT.m_RowData.assign(this->m_ColData.begin(), this->m_ColData.end());
	mtxAT.m_RowIndex.assign(this->m_ColIndex.begin(), this->m_ColIndex.end());
	mtxAT.m_ColData.assign(this->m_RowData.begin(), this->m_RowData.end());
	mtxAT.m_ColIndex.assign(this->m_RowIndex.begin(), this->m_RowIndex.end());

	/*
	int i;
	size_t j, nElement;
	
    int index;
	double data;
	for(i = 0; i < m_nRow; ++ i)
	{
		nElement = m_RowIndex[i].size();
		for(j = 0; j < nElement; ++ j)
		{
			index = m_RowIndex[i][j];
			data = m_RowData[i][j];
			mtxAT.m_RowIndex[index].push_back(i);
			mtxAT.m_RowData[index].push_back(data);
		}
	}

	for (i = 0; i < m_nCol; i++)
	{
		nElement = m_ColIndex[i].size();
		for (j = 0; j < nElement; j++)
		{
			index = m_ColIndex[i][j];
			data = m_ColData[i][j];
			mtxAT.m_ColIndex[index].push_back(i);
			mtxAT.m_ColData[index].push_back(data);
		}
	}

	mtxAT.m_IsIndexSeqed = this->m_IsIndexSeqed;*/
}
void CMeshSparseMatrix::Multiply(CMeshSparseMatrix& other, CMeshSparseMatrix& result, bool fill_row)
{
	// 首先检查行列数是否符合要求
	assert (m_nCol == other.GetRowNum());
 
	// construct the object we are going to return
	result.SetRowCol(GetRowNum(), other.GetColNum());

    printf("Start AT*A. \n");

	//
	//other.CalColInfo();

	// 矩阵乘法，即
	//
	// [A][B][C]   [G][H]     [A*G + B*I + C*K][A*H + B*J + C*L]
	// [D][E][F] * [I][J] =   [D*G + E*I + F*K][D*H + E*J + F*L]
	//             [K][L]
	//
	for(int i = 0; i < m_nRow; ++ i)
	{
		size_t rElement = m_RowIndex[i].size();
		for(size_t j = 0; j < rElement; ++ j)
		{
			int index = m_RowIndex[i][j];
			double data1 = m_RowData[i][j];

			size_t cElement = other.m_RowIndex[index].size();

			for (size_t n = 0; n < cElement; n++)
			{
				int cindex = other.m_RowIndex[index][n];
				double data2 = other.m_RowData[index][n];
				double d = data1 * data2;

				bool is_found = false;
				for(size_t k = 0; k < result.m_ColIndex[cindex].size(); ++ k)
				{
					if(result.m_ColIndex[cindex][k] == i)
					{
						result.m_ColData[cindex][k] += d;
						is_found = true;
						break;
					}
				}
				if (!is_found)
				{
					result.m_ColIndex[cindex].push_back(i);
					result.m_ColData[cindex].push_back(d);
				}

				/*
				result.GetElement(i, cindex, d);
				d += data1 * data2;
				result.SetElement(i, cindex, d);
				*/
			}
		}
	}

	if (fill_row)
	{
		result.fill_row_info();
	}
	printf("AT*A is finished.\n");
}
void CMeshSparseMatrix::Multiply(CMeshSparseMatrix& leftmtx, CMeshSparseMatrix& rightmtx, CMeshSparseMatrix& result)
{
	// result = leftmtx * this * rightmtx.
	CMeshSparseMatrix mtx, tmpmtx;
	MakeCopy(mtx);
	leftmtx.Multiply(mtx, tmpmtx);
	tmpmtx.Multiply(rightmtx, result);
}
void CMeshSparseMatrix::Add(CMeshSparseMatrix& mtxOther)
{
	// 首先检查行列数是否符合要求
	assert (m_nCol == mtxOther.GetColNum() && m_nRow == mtxOther.GetRowNum());

	//printf("Start A + B. \n");

	for (int i = 0; i < m_nRow; i++)
	{
		size_t rElement = mtxOther.m_RowIndex[i].size();

		for(size_t j = 0; j < rElement; ++ j)
		{
			int index = mtxOther.m_RowIndex[i][j];
			double data = mtxOther.m_RowData[i][j];

			double d;
			GetElement(i, index, d);
			d += data;
			SetElement(i, index, d);
		}
	}

	//printf("A + B is finished.\n");
}
void CMeshSparseMatrix::MultiplyDouble(double weight)
{
	for (int i = 0; i < m_nRow; i++)
	{
		size_t rElement = m_RowData[i].size();

		for (size_t j = 0; j < rElement; j++)
		{
			m_RowData[i][j] *= weight;
		}
	}

	for (int i = 0; i < m_nCol; i++)
	{
		size_t rElement = m_ColData[i].size();

		for (size_t j = 0; j < rElement; j++)
		{
			m_ColData[i][j] *= weight;
		}
	}
}
void CMeshSparseMatrix::MultiplyDouble4Row(vector<double> weight)
{
	assert(weight.size() == GetRowNum());

	for (int i = 0; i < m_nRow; i++)
	{
		size_t rElement = m_RowData[i].size();

		for (size_t j = 0; j < rElement; j++)
		{
			m_RowData[i][j] *= weight[i];
		}
	}

	for (int i = 0; i < m_nCol; i++)
	{
		size_t rElement = m_ColData[i].size();

		for (size_t j = 0; j < rElement; j++)
		{
			m_ColData[i][j] *= weight[m_ColIndex[i][j]];
		}
	}
}
bool CMeshSparseMatrix::MultiplyVector(DoubleArray& vec, DoubleArray& result, bool lefMulti)
{
	// if right multiply vector.
	if (!lefMulti)
	{
		if(m_nCol != vec.size())
			return false;

		result.resize(m_nRow);
		int i;
		size_t j, n;
		for(i = 0; i < m_nRow; ++ i)
		{
			n = m_RowIndex[i].size();
			double& value = result[i];
			value = 0.0;
			for(j = 0; j < n; ++ j)
				value += (m_RowData[i][j] * vec[m_RowIndex[i][j]]);
		}
	}
	// if left multiply vector.
	else
	{
		if(m_nRow != vec.size())
			return false;

		result.resize(m_nCol);

		int i;
		size_t j, n;
		for(i = 0; i < m_nCol; ++ i)
		{
			n = m_ColIndex[i].size();
			double& value = result[i];
			value = 0.0;
			for(j = 0; j < n; ++ j)
				value += (m_ColData[i][j] * vec[m_ColIndex[i][j]]);
		}
	}

	return true;
}
void CMeshSparseMatrix::InsertMatrix(int LURow, int LUCol, CMeshSparseMatrix& inmtx)
{
	assert(GetRowNum() >= (inmtx.GetRowNum() + LURow) &&
		GetColNum() >= (inmtx.GetColNum() + LUCol));


	
	// copy row data and index here.
	for (int i = 0; i < inmtx.GetRowNum(); i++)
	{
		vector<int>& rowindex = inmtx.m_RowIndex[i];
		vector<double>& rowdata = inmtx.m_RowData[i];

		vector<int>& nowrowindex = this->m_RowIndex[LURow + i];
		vector<double>& nowrowdata = this->m_RowData[LURow + i];
		nowrowdata.insert(nowrowdata.end(), rowdata.begin(), rowdata.end());
		for (size_t j = 0; j < rowindex.size(); j++)
		{
			nowrowindex.push_back(rowindex[j] + LUCol);
		}
	}

	// copy col data and index here.
	for (int i = 0; i < inmtx.GetColNum(); i++)
	{
		vector<int>& colindex = inmtx.m_ColIndex[i];
		vector<double>& coldata = inmtx.m_ColData[i];

		vector<int>& nowcolindex = this->m_ColIndex[LUCol + i];
		vector<double>& nowcoldata = this->m_ColData[LUCol + i];
		nowcoldata.insert(nowcoldata.end(), coldata.begin(), coldata.end());
		for (size_t j = 0; j < colindex.size(); j++)
		{
			nowcolindex.push_back(colindex[j] + LURow);
		}
	}
}
void CMeshSparseMatrix::GetMatrixLesseqCol(int col, vector<pair<int, double> >& colInfo)
{
	vector<double>& colData = m_ColData[col];
	vector<int>& colIndex = m_ColIndex[col];

	for (size_t j = 0; j < colData.size(); j++)
	{
		colInfo.push_back(make_pair(colIndex[j], colData[j]));
	}

	sort(colInfo.begin(), colInfo.end(), DataCompare());
}
void CMeshSparseMatrix::MakeMatrixIndexLessSeq()
{
	if (m_IsIndexSeqed)
	{
		return;
	}

	// re-sort row data and index here.
	for (size_t i = 0; i < m_RowData.size(); i++)
	{
		vector<double>& rowData = m_RowData[i];
		vector<int>& rowIndex = m_RowIndex[i];

		vector<pair<int, double> > rowInfo;

		for (size_t j = 0; j < rowData.size(); j++)
		{
			rowInfo.push_back(make_pair(rowIndex[j], rowData[j]));
		}

		sort(rowInfo.begin(), rowInfo.end(), DataCompare());

		for (size_t j = 0; j < rowData.size(); j++)
		{
			rowIndex[j] = rowInfo[j].first;
			rowData[j] = rowInfo[j].second;
		}
	}

	// re-sort column data and index here.
	for (size_t i = 0; i < m_ColData.size(); i++)
	{
		vector<double>& colData = m_ColData[i];
		vector<int>& colIndex = m_ColIndex[i];

		vector<pair<int, double> > colInfo;

		for (size_t j = 0; j < colData.size(); j++)
		{
			colInfo.push_back(make_pair(colIndex[j], colData[j]));
		}

		sort(colInfo.begin(), colInfo.end(), DataCompare());

		for (size_t j = 0; j < colData.size(); j++)
		{
			colIndex[j] = colInfo[j].first;
			colData[j] = colInfo[j].second;
		}
	}

	this->m_IsIndexSeqed = true;
}
double CMeshSparseMatrix::CalFrobeniusNorm()
{
	double sum = 0;

	for (size_t i = 0; i < m_ColData.size(); i++)
	{
		vector<double>& colData = m_ColData[i];

		for (size_t j = 0; j < colData.size(); j++)
		{
			sum += colData[j] * colData[j];
		}
	}

	return sqrt(sum);
}
bool CMeshSparseMatrix::LoadSparseMatrix(string filename)
{
	ifstream file(filename.c_str());
	//ifstream file(filename.c_str(), ios::binary);

	if (!file)
	{
		return false;
	}

	/*
	
	int row, col;
	file.read((char *)&row, sizeof(int));
	file.read((char *)&col, sizeof(int));

	this->SetRowCol(row, col);

	for (size_t i = 0; i < col; i++)
	{
		int col_size;
		file.read((char *)&col_size, sizeof(int));

		vector<int>& col_index = m_ColIndex[i];//(col_size);
		vector<double>& col_data = m_ColData[i];//(col_size);
		col_index.resize(col_size);
		col_data.resize(col_size);
		
		file.read((char *)&col_index[0], sizeof(int) * col_index.size());
		file.read((char *)&col_data[0], sizeof(double) * col_data.size());

	
		//for (size_t j = 0; j < col_size; j++)
		//{
		//	this->SetElement(col_index[j], i, col_data[j]);
		//}
	}

	for (size_t i = 0; i < row; i++)
	{
		int row_size;
		file.read((char *)&row_size, sizeof(int));

		vector<int>& row_index = m_RowIndex[i];//(col_size);
		vector<double>& row_data = m_RowData[i];//(col_size);
		row_index.resize(row_size);
		row_data.resize(row_size);

		file.read((char *)&row_index[0], sizeof(int) * row_index.size());
		file.read((char *)&row_data[0], sizeof(double) * row_data.size());

	}*/


	
	int row, col;
	file >> row >> col;

	this->SetRowCol(row, col);

	//
	size_t col_size;
	double col_value;
	for (size_t i = 0; i < col; i++)
	{
		file >> col_size;
		vector<int> col_index(col_size);

		//
		for (size_t j = 0; j < col_size; j++)
		{
			file >> col_index[j];
		}

		//
		for (size_t j = 0; j < col_size; j++)
		{
			file >> col_value;
			this->SetElement(col_index[j], i, col_value);
		}
	}

	file.close();
	return true;
}
bool CMeshSparseMatrix::SaveSparseMatrix(string filename)
{
	ofstream file(filename.c_str());
	//ofstream file(filename.c_str(), ios::binary);

	if (!file)
	{
		return false;
	}

	/*
	file.write((const char *)&m_nRow, sizeof(int));
	file.write((const char *)&m_nCol, sizeof(int));

	for (size_t i = 0; i < m_nCol; i++)
	{
		int col_size = m_ColIndex[i].size();
		file.write((const char *)&col_size, sizeof(int));

		IndexArray& col_index = m_ColIndex[i];
		DoubleArray& col_data = m_ColData[i];

		file.write((const char *)&col_index[0], sizeof(int) * col_index.size());
		file.write((const char *)&col_data[0], sizeof(double) * col_data.size());
	}

	for (size_t i = 0; i < m_nRow; i++)
	{
		int row_size = m_RowIndex[i].size();
		file.write((const char *)&row_size, sizeof(int));

		IndexArray& row_index = m_RowIndex[i];
		DoubleArray& row_data = m_RowData[i];

		file.write((const char *)&row_index[0], sizeof(int) * row_index.size());
		file.write((const char *)&row_data[0], sizeof(double) * row_data.size());
	}*/


	file << m_nRow << ' ' << m_nCol << '\n';
	
    size_t line_num = 0;
	
	for (size_t i = 0; i < m_ColIndex.size(); i++)
	{
		file  << m_ColIndex[i].size() << '\n';

		IndexArray& col_index = m_ColIndex[i];
		DoubleArray& col_data = m_ColData[i];
		
		line_num = 0;
		for (size_t j = 0; j < col_index.size(); j++)
		{
			file << col_index[j] << ' ';
			line_num++;
			if (line_num == 10)
			{
				line_num = 0;
				file << '\n';
			}
		}
		if (line_num < 10)
		{
			file << '\n';
		}

		//
		line_num = 0;
		for (size_t j = 0; j < col_data.size(); j++)
		{
			file << col_data[j] << ' ';
			line_num++;
			if (line_num == 10)
			{
				line_num = 0;
				file << '\n';
			}
		}
		if (line_num < 10)
		{
			file << '\n';
		}
	}

	file.close();
	return true;
}
void CMeshSparseMatrix::fill_row_info()
{
	for (size_t i = 0; i < m_nCol; i++)
	{
		vector<int>& col_index = m_ColIndex[i];
		vector<double>& col_data = m_ColData[i];

		for (size_t j = 0; j < col_index.size(); j++)
		{
			m_RowIndex[col_index[j]].push_back(i);
			m_RowData[col_index[j]].push_back(col_data[j]);
		}
	}
}
void CMeshSparseMatrix::MatrixLUTri2Full()
{
	CMeshSparseMatrix temp;
	MakeCopy(temp);

	for (size_t i = 0; i < temp.m_RowIndex.size(); i++)
	{
		int row = (int) i;
		vector<int>& row_index = temp.m_RowIndex[i];
		vector<double>& row_data = temp.m_RowData[i];

		for (size_t j = 0; j < row_index.size(); j++)
		{
			// if not diagnoal
			if (row_index[j] != row)
			{
				SetElement(row_index[j], row, row_data[j]);
			}
		}
	}
}
double CMeshSparseMatrix::MaxAii()
{
	double max_val = -1e20;

	for (size_t i = 0; i < m_RowIndex.size(); i++)
	{
		int row = (int) i;
		vector<int>& row_index = m_RowIndex[i];
		vector<double>& row_data = m_RowData[i];

		for (size_t j = 0; j < row_index.size(); j++)
		{
			double fabs_val = fabs(row_data[j]);
			if (row_index[j] == row && fabs_val > max_val)
			{
				max_val = fabs_val;
			}
		}
	}

	return max_val;
}
double CMeshSparseMatrix::AvgAii()
{
	double sum = 0;
	for (size_t i = 0; i < m_RowIndex.size(); i++)
	{
		int row = (int) i;
		vector<int>& row_index = m_RowIndex[i];
		vector<double>& row_data = m_RowData[i];

		for (size_t j = 0; j < row_index.size(); j++)
		{
			double fabs_val = fabs(row_data[j]);
			if (row_index[j] == row)
			{
				sum += fabs_val;
			}
		}
	}

	return sum / GetColNum();
}