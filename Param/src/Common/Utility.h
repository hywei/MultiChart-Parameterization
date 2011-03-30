/* ================== Library Information ================== */
// [Name] 
// MeshLib Library
//
// [Developer]
// Xu Dong
// State Key Lab of CAD&CG, Zhejiang University
// 
// [Date]
// 2005-08-05
//
// [Goal]
// A general, flexible, versatile and easy-to-use mesh library for research purpose.
// Supporting arbitrary polygonal meshes as input, but with a 
// primary focus on triangle meshes.

/* ================== File Information ================== */
// [Name]
// Utility.h
//
// [Developer]
// Xu Dong
// State Key Lab of CAD&CG, Zhejiang University
// 
// [Date]
// 2005-08-05
//
// [Goal]
// Defining various utility functions



#include "BasicDataType.h"
#include "../Numerical/matrix.h"
#include "../Numerical/Rotation.h"
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <ctype.h>

#ifdef WIN32
#include <float.h> // for ::_isnan()
#else
#include <stdlib.h>
#endif

#pragma once

typedef enum
{
	IS_EDGE,  // Intersection is on the edge (anywhere on the edge but at the end-points)
	IS_POINT, // Intersection with one of the two end-points
	IS_NONE
} IntersectionType;

class Utility
{
public:
    // string utilities
	void ResolveFileName(std::string filename, std::string& file_path, 
		std::string& file_title, std::string& file_ext);
	void MakeLower(std::string& str);
	void MakeUpper(std::string& str);

    // Free various vectors
    template <class T>
    void FreeVector(T& arr)
    {
        T tmp;
        arr.clear();
        arr.swap(tmp);
    }

    // Flag utilities
    template <class T>
        void SetFlag(T& flag_ele, T flag) { flag_ele |= flag; }

    template <class T>
        void ClearFlag(T& flag_ele, T flag) { flag_ele &= ~flag; }

    template <class T>
        bool IsSetFlag(T& flag_ele, T flag) { return ((flag_ele&flag) == flag); }

    template <class T>
        void ToggleFlag(T& flag_ele, T flag)
    {
        if(IsSetFlag(flag_ele, flag))
            ClearFlag(flag_ele, flag);
        else
            SetFlag(flag_ele, flag);
    }

    // Debug output utility
    template <class T>  // Make sure class T is a vector
	void PrintArray(T& arr, std::string arrName)
    {
        int n = arr.size();
        std::cout << '\n';
        for(int i = 0; i < n; ++ i)
        {
            std::cout << arrName;
            printf("[%8d] = ", i);
            std::cout << arr[i] << '\n';
        }
    }
    // Mean value and variance calculation
    template <class T>  // Make sure class T is a vector of basic data type (int, float, double, etc... )
	void CalMeanAndStdVar(std::vector<T>& arr, T& mean, T& std_var)
    {
        int n = arr.size();
        if(n == 0)
            return;

        double Mean = 0.0, Var = 0.0;
        for(int i = 0; i < n; ++ i)
            Mean += arr[i];
        Mean /= (double)n;

        for(int i = 0; i < n; ++ i)
        {
            double tmp = arr[i]-Mean;
            Var += tmp*tmp;
        }
        Var /= (double)n;

        mean = (T)Mean;
        std_var  = (T)sqrt(Var);
    }

    // 
    template <class T>  // Make sure class T is a vector of basic data type (int, float, double, etc... )
	bool IsInVector(std::vector<T>& arr, T& value)
    {
		size_t n = arr.size();
        if(n == 0)
            return false;

        for(size_t i = 0; i < n; ++ i)
		{
			if (arr[i] == value)
			{
				return true;
			}
		}

		return false;
	}

	template <class T>  // Make sure class T is a vector of basic data type (int, float, double, etc... )
	int IndexInVector(std::vector<T>& arr, T value)
	{
		int n = (int) arr.size();
		if(n == 0)
			return -1;

		for(int i = 0; i < n; ++ i)
		{
			if (arr[i] == value)
			{
				return i;
			}
		}

		return -1;
	}
	template <class T>  // Make sure class T is a vector of basic data type (int, float, double, etc... )
	void ReverseVector(std::vector<T>& arr)
	{
      std::vector<T> tmp_arr(arr);

		int n = (int) arr.size();
		if(n == 0)
			return;

		for(int i = 0; i < n; ++ i)
		{
			arr[n- 1 - i] = tmp_arr[i]; 
		}
	}

    // Data set filtering
    template <class T>  // Make sure class T is a vector of basic data type (int, float, double, etc... )
	int DataClamp(std::vector<T>& arr, T mean, T std_var, double coef)
    {
        int n = arr.size();
        int nClamp = 0;
        for(int i = 0; i < n; ++ i)
        {
            T& data = arr[i];
            double tmp = (data-mean)/std_var;
            if(tmp > coef)
            {
                data = coef*std_var + mean;
                ++ nClamp;
            }
            else if(tmp < -coef)
            {
                data = -coef*std_var + mean;
                ++ nClamp;
            }
        }
        
        return nClamp;
    }

	double VectorMultiplyVector(std::vector<double>& vec1, std::vector<double>& vec2)
	{
		if (vec1.size() != vec2.size())
		{
			return -1;
		}

		double value = 0;
		for (size_t i = 0; i < vec1.size(); i++)
		{
			value += vec1[i] * vec2[i];
		}

		return value;
	}

	void VectorMultiplyNumber(std::vector<double>& vec, double number)
	{
		for (size_t i = 0; i < vec.size(); i++)
		{
			vec[i] *= number;
		}
	}

	double NoramlizeVector(std::vector<double>& vec)
	{
		double sum = 0;
		for (size_t i = 0; i < vec.size(); i++)
		{
			sum += vec[i] * vec[i];
		}

		sum = sqrt(sum);

		for (size_t i = 0; i < vec.size(); i++)
		{
			vec[i] /= sum;
		}
		return sum;
	}
	double CenterVector(std::vector<double>& vec)
	{
		double sum = 0;
		for (size_t i = 0; i < vec.size(); i++)
		{
			sum += vec[i];
		}

		sum /= vec.size();

		for (size_t i = 0; i < vec.size(); i++)
		{
			vec[i] -= sum;
		}

		return sum;
	}

	double VectorDistance(std::vector<double>& vec1, std::vector<double>& vec2)
	{
		//
		double dis = 0;

		std::vector<double> nvec1(vec1);
		std::vector<double> nvec2(vec2);
		NoramlizeVector(nvec1);
		NoramlizeVector(nvec2);

		for (size_t i = 0; i < vec1.size(); i++)
		{
			dis += (nvec1[i] - nvec2[i]) * (nvec1[i] - nvec2[i]);
		}

		return dis;
	}

	//////////////////////////////////////////////////////////////////////////

	bool is_nan(double x) 
	{
#ifdef WIN32
		return (_isnan(x) != 0) || (_finite(x) == 0) ;
#else
		return isnan(x) || !finite(x);
#endif
	}    

    //////////////////////////////////////////////////////////////////////////
	bool CalTransMatrix(Coord src_v[3], Coord tar_v[3], CMatrix& TransMtx, bool bNormalScale = false);
    bool CalRigidTransMatrix(Coord src_v[3], Coord tar_v[3], CMatrix& TransMtx);
	bool CalTransMatrixAndSplit(Coord src_v[3], Coord tar_v[3], CMatrix& TransMtx, QUATERNION& quat, double S[3][3]);
	bool MatrixSplit(CMatrix TransMtx, QUATERNION& quat, double S[3][3]);
	Coord VectorTransform(CMatrix& TransMtx, Coord v);

    //////////////////////////////////////////////////////////////////////////
    // Given normal and another vector, compute the local frame
    inline void CalLocalFrame(Coord vector, Coord normal, CoordArray& Axis)
    {
        Axis[1] = (cross(normal, vector)).unit();
        Axis[0] = cross(Axis[1], normal);
        Axis[2] = normal;
    }

    template<class T> // Make sure that T is a index array
    void Unique(T& arr)
    {
        sort(arr.begin(), arr.end());
        arr.erase(unique(arr.begin(), arr.end()), arr.end());
    }

	//////////////////////////////////////////////////////////////////////////
	bool SegmentIntersection(double s1x1, double s1y1, double s1x2, double s1y2,
		                     double s2x1, double s2y1, double s2x2, double s2y2,
							 double& r, double &s);
	bool SegmentIntersection(const Coord& A, const Coord& B, const Coord& C, const Coord& D, double& r, double& s);
	bool TEdgeIntersection(std::pair<Coord, Coord>& edge, std::pair<Coord, Coord>& radial, 
		std::pair<Coord, IntersectionType>& interSection);

};
