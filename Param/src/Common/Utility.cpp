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
// Utility.cpp
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

#include "Utility.h"
#include <assert.h>
//#include "complex.h"

void Utility::ResolveFileName(std::string filename, std::string& file_path, 
							  std::string& file_title, std::string& file_ext)
{
	size_t i;
    size_t lash_index = filename.rfind('/');
	size_t dot_index = filename.rfind('.');
	size_t length = filename.length();

	for(i = 0; i <= lash_index; ++ i)	
		file_path += filename[i];
	for(i = lash_index+1; i < dot_index; ++ i)
		file_title += filename[i];
	for(i = dot_index; i < length; ++ i)
		file_ext += filename[i];
}

void Utility::MakeLower(std::string& str)
{
    size_t length = str.length();
    for(size_t i = 0; i < length; ++ i)
    {
        char& ch = str[i];
        if(ch >= 'A' && ch <= 'Z')
            ch = ch - 'A' + 'a';
    }
}

void Utility::MakeUpper(std::string& str)
{
    size_t length = str.length();
    for(size_t i = 0; i < length; ++ i)
    {
        char& ch = str[i];
        if(ch >= 'a' && ch <= 'a')
            ch = ch - 'a' + 'A';
    }
}

// Calculate the transform matrix for the source triangle src_v[3] and the target triangle tar_v[3]
bool Utility::CalTransMatrix(Coord src_v[3], Coord tar_v[3], CMatrix& TransMtx, bool bNormalScale)
{
    int j, k;
    Coord src_fNormal, tar_fNormal, temp[3];
    Coord src_v4, tar_v4;
    CMatrix SrcMtx(3, 3), TarMtx(3, 3);
    TransMtx.SetRowCol(3, 3);
    
    // Setup Source Matrix
    temp[0] = src_v[1]-src_v[0];
    temp[1] = src_v[2]-src_v[0];
    temp[2] = cross(temp[0], temp[1]);
    double src_area = temp[2].abs();
    temp[2].normalize();
    for(j = 0; j < 3; ++ j)
    {
        for(k = 0; k < 3; ++ k)
            SrcMtx.SetElement(k, j, temp[j][k]);
    }
    
    // Setup Target Matrix
    temp[0] = tar_v[1]-tar_v[0];
    temp[1] = tar_v[2]-tar_v[0];
    temp[2] = cross(temp[0], temp[1]);
	double tar_area = temp[2].abs();
    temp[2].normalize();

	// Stretch according to area ratio
	double area_ratio = (bNormalScale) ? sqrt(tar_area/src_area) : 1.0;
	temp[2] *= area_ratio;

    for(j = 0; j < 3; ++ j)
    {
        for(k = 0; k < 3; ++ k)
            TarMtx.SetElement(k, j, temp[j][k]);
    }
    
    if(SrcMtx.InvertGaussJordan())
    {
        TransMtx = TarMtx*SrcMtx;
        return true;
    }
    else
    {
        fprintf(stdout, "\nSource Matrix Invert Error!\n");
        return false;
    }
}

bool Utility::CalRigidTransMatrix(Coord src_v[3], Coord tar_v[3], CMatrix& TransMtx)
{
	int j, k;
	Coord src_fNormal, tar_fNormal, temp[3];
	Coord src_v4, tar_v4;
	CMatrix SrcMtx(3, 3), TarMtx(3, 3);
	TransMtx.SetRowCol(3, 3);

	// Setup Source Matrix
	temp[0] = src_v[1]-src_v[0];
	temp[1] = src_v[2]-src_v[0];
	temp[2] = cross(temp[0], temp[1]);
	temp[2].normalize();
	for(j = 0; j < 3; ++ j)
	{
		for(k = 0; k < 3; ++ k)
			SrcMtx.SetElement(k, j, temp[j][k]);
	}

	// Setup Target Matrix
	temp[0] = tar_v[1]-tar_v[0];
	temp[1] = tar_v[2]-tar_v[0];
	temp[2] = cross(temp[0], temp[1]);
	temp[2].normalize();
	for(j = 0; j < 3; ++ j)
	{
		for(k = 0; k < 3; ++ k)
			TarMtx.SetElement(k, j, temp[j][k]);
	}

	CMatrix mtxU, mtxV;
	CMatrix mtxUV, mtxS;
	if(SrcMtx.InvertGaussJordan())
	{
		TransMtx = TarMtx*SrcMtx;
		if(TransMtx.SplitUV(mtxU, mtxV, 1.0e-8))	// SVD
		{
			double det_u = mtxU.GetElement(0,0)*mtxU.GetElement(1,1)*mtxU.GetElement(2,2) + 
				mtxU.GetElement(0,1)*mtxU.GetElement(1,2)*mtxU.GetElement(2,0) + 
				mtxU.GetElement(1,0)*mtxU.GetElement(2,1)*mtxU.GetElement(0,2) - 
				mtxU.GetElement(0,2)*mtxU.GetElement(1,1)*mtxU.GetElement(2,0) - 
				mtxU.GetElement(0,1)*mtxU.GetElement(1,0)*mtxU.GetElement(2,2) - 
				mtxU.GetElement(0,0)*mtxU.GetElement(2,1)*mtxU.GetElement(1,2);
			double det_v = mtxV.GetElement(0,0)*mtxV.GetElement(1,1)*mtxV.GetElement(2,2) + 
				mtxV.GetElement(0,1)*mtxV.GetElement(1,2)*mtxV.GetElement(2,0) + 
				mtxV.GetElement(1,0)*mtxV.GetElement(2,1)*mtxV.GetElement(0,2) - 
				mtxV.GetElement(0,2)*mtxV.GetElement(1,1)*mtxV.GetElement(2,0) - 
				mtxV.GetElement(0,1)*mtxV.GetElement(1,0)*mtxV.GetElement(2,2) - 
				mtxV.GetElement(0,0)*mtxV.GetElement(2,1)*mtxV.GetElement(1,2);
			
			TransMtx = mtxU*mtxV;
		}
		else
		{
			fprintf(stdout, "Transformation Matrix SVD Error!\n");
			return false;
		}
	}
	else
	{
		fprintf(stdout, "\nSource Matrix Invert Error!\n");
		return false;
	}
	
	return true;
}

bool Utility::CalTransMatrixAndSplit(Coord src_v[3], Coord tar_v[3], CMatrix& TransMtx, QUATERNION& quat, double S[3][3])
{
	int j, k;
	Coord src_fNormal, tar_fNormal, temp[3];
	Coord src_v4, tar_v4;
	CMatrix SrcMtx(3, 3), TarMtx(3, 3);
	TransMtx.SetRowCol(3, 3);

	// Setup Source Matrix
	temp[0] = src_v[1]-src_v[0];
	temp[1] = src_v[2]-src_v[0];
	temp[2] = cross(temp[0], temp[1]);
	temp[2].normalize();
	for(j = 0; j < 3; ++ j)
	{
		for(k = 0; k < 3; ++ k)
			SrcMtx.SetElement(k, j, temp[j][k]);
	}

	// Setup Target Matrix
	temp[0] = tar_v[1]-tar_v[0];
	temp[1] = tar_v[2]-tar_v[0];
	temp[2] = cross(temp[0], temp[1]);
	temp[2].normalize();
	for(j = 0; j < 3; ++ j)
	{
		for(k = 0; k < 3; ++ k)
			TarMtx.SetElement(k, j, temp[j][k]);
	}

	CMatrix mtxU, mtxV;
	CMatrix mtxUV, mtxS;
	if(SrcMtx.InvertGaussJordan())
	{
		TransMtx = TarMtx*SrcMtx;
		if(TransMtx.SplitUV(mtxU, mtxV, 1.0e-8))	// SVD
		{
			double det_u = mtxU.GetElement(0,0)*mtxU.GetElement(1,1)*mtxU.GetElement(2,2) + 
				mtxU.GetElement(0,1)*mtxU.GetElement(1,2)*mtxU.GetElement(2,0) + 
				mtxU.GetElement(1,0)*mtxU.GetElement(2,1)*mtxU.GetElement(0,2) - 
				mtxU.GetElement(0,2)*mtxU.GetElement(1,1)*mtxU.GetElement(2,0) - 
				mtxU.GetElement(0,1)*mtxU.GetElement(1,0)*mtxU.GetElement(2,2) - 
				mtxU.GetElement(0,0)*mtxU.GetElement(2,1)*mtxU.GetElement(1,2);
			double det_v = mtxV.GetElement(0,0)*mtxV.GetElement(1,1)*mtxV.GetElement(2,2) + 
				mtxV.GetElement(0,1)*mtxV.GetElement(1,2)*mtxV.GetElement(2,0) + 
				mtxV.GetElement(1,0)*mtxV.GetElement(2,1)*mtxV.GetElement(0,2) - 
				mtxV.GetElement(0,2)*mtxV.GetElement(1,1)*mtxV.GetElement(2,0) - 
				mtxV.GetElement(0,1)*mtxV.GetElement(1,0)*mtxV.GetElement(2,2) - 
				mtxV.GetElement(0,0)*mtxV.GetElement(2,1)*mtxV.GetElement(1,2);
			
			mtxUV = mtxU*mtxV;
			mtxS = mtxV.Transpose()*TransMtx*mtxV;
			
			// Polar decomposition			
			double RUV[3][3];
			for(j = 0; j < 3; ++ j)
			{
				for(k = 0; k < 3; ++ k)
				{
					RUV[j][k] = mtxUV.GetElement(j, k);
					S[j][k] = mtxS.GetElement(j, k);
				}
			}
			MatrixToQuaternion2(quat, RUV);
		}
		else
		{
			fprintf(stdout, "Transformation Matrix SVD Error!\n");
			return false;
		}
	}
	else
	{
		fprintf(stdout, "\nSource Matrix Invert Error!\n");
		return false;
	}
	
	return true;
}

bool Utility::MatrixSplit(CMatrix TransMtx, QUATERNION& quat, double S[3][3])
{
	CMatrix mtxU, mtxV;
	CMatrix mtxUV, mtxS;
	int j, k;
	double RUV[3][3];
	if(TransMtx.SplitUV(mtxU, mtxV, 1.0e-8))	// SVD
	{
		double det_u = mtxU.GetElement(0,0)*mtxU.GetElement(1,1)*mtxU.GetElement(2,2) + 
			mtxU.GetElement(0,1)*mtxU.GetElement(1,2)*mtxU.GetElement(2,0) + 
			mtxU.GetElement(1,0)*mtxU.GetElement(2,1)*mtxU.GetElement(0,2) - 
			mtxU.GetElement(0,2)*mtxU.GetElement(1,1)*mtxU.GetElement(2,0) - 
			mtxU.GetElement(0,1)*mtxU.GetElement(1,0)*mtxU.GetElement(2,2) - 
			mtxU.GetElement(0,0)*mtxU.GetElement(2,1)*mtxU.GetElement(1,2);
		double det_v = mtxV.GetElement(0,0)*mtxV.GetElement(1,1)*mtxV.GetElement(2,2) + 
			mtxV.GetElement(0,1)*mtxV.GetElement(1,2)*mtxV.GetElement(2,0) + 
			mtxV.GetElement(1,0)*mtxV.GetElement(2,1)*mtxV.GetElement(0,2) - 
			mtxV.GetElement(0,2)*mtxV.GetElement(1,1)*mtxV.GetElement(2,0) - 
			mtxV.GetElement(0,1)*mtxV.GetElement(1,0)*mtxV.GetElement(2,2) - 
			mtxV.GetElement(0,0)*mtxV.GetElement(2,1)*mtxV.GetElement(1,2);
		
		mtxUV = mtxU*mtxV;
		mtxS = mtxV.Transpose()*TransMtx*mtxV;
		
		// Polar decomposition
		for(j = 0; j < 3; ++ j)
		{
			for(k = 0; k < 3; ++ k)
			{
				RUV[j][k] = mtxUV.GetElement(j, k);
				S[j][k] = mtxS.GetElement(j, k);
			}
		}
		MatrixToQuaternion2(quat, RUV);
		return true;
	}
	else
	{
		fprintf(stdout, "Transformation Matrix SVD Error!\n");
		return false;
	}
}

Coord Utility::VectorTransform(CMatrix& TransMtx, Coord v)
{
	Coord nv(0.0, 0.0, 0.0);
	assert(TransMtx.GetNumRows() == 3 && TransMtx.GetNumColumns() == 3);
	for(int i = 0; i < 3; ++ i)
	{
		for(int j = 0; j < 3; ++ j)
			nv[i] += TransMtx.GetElement(i, j)*v[j];
	}

	return nv;
}
bool Utility::SegmentIntersection(double s1x1, double s1y1, double s1x2, double s1y2,
								  double s2x1, double s2y1, double s2x2, double s2y2, 
								  double& r, double &s)
{
	// we need to normalize the vector, then compute 
	// kross(d0, d1)
	double sqrEpsilon = 1e-6;
	double denominator = (s1x2 - s1x1) * (s2y2 - s2y1) - (s1y2 - s1y1) * (s2x2 - s2x1);
	double sqrdenominator = denominator * denominator;
	double sqrLen0 = (s1x2 - s1x1) * (s1x2 - s1x1) + (s1y2 - s1y1) * (s1y2 - s1y1);
	double sqrLen1 = (s2x2 - s2x1) * (s2x2 - s2x1) + (s2y2 - s2y1) * (s2y2 - s2y1);

	if (sqrdenominator > sqrEpsilon * sqrLen0 * sqrLen1)
	{

		r = ((s1y1 - s2y1) * (s2x2 - s2x1) - (s1x1 - s2x1) * (s2y2 - s2y1)) / denominator;
		s = ((s1y1 - s2y1) * (s1x2 - s1x1) - (s1x1 - s2x1) * (s1y2 - s1y1)) / denominator;

		return true;
	}
	else
	{
		// Parallel
		return false;
	}
}
bool Utility::SegmentIntersection(const Coord& A, const Coord& B, const Coord& C, const Coord& D, double& r, double& s)
{
	// Cross product
	Coord BA = B - A;
	Coord DC = D - C;
	Coord BAcDC = cross(BA, DC).unit();

	// Find best axis projection
	double x = fabs(dot(BAcDC, Coord(1, 0, 0)));
	double y = fabs(dot(BAcDC, Coord(0, 1, 0)));
	double z = fabs(dot(BAcDC, Coord(0, 0, 1)));

	//
	if(x > y)
	{
		if(x > z)
			return SegmentIntersection(A[1], A[2], B[1], B[2], C[1], C[2], D[1], D[2], r, s);
		else
			return SegmentIntersection(A[0], A[1], B[0], B[1], C[0], C[1], D[0], D[1], r, s);
	}

	if(y > z)
		return SegmentIntersection(A[0], A[2], B[0], B[2], C[0], C[2], D[0], D[2], r, s);

	return SegmentIntersection(A[0], A[1], B[0], B[1], C[0], C[1], D[0], D[1], r, s);
}
bool Utility::TEdgeIntersection(std::pair<Coord, Coord>& edge, std::pair<Coord, Coord>& radial, 
								std::pair<Coord, IntersectionType>& interSection)
{
	// Intersection position on line 1 : m_intersection (0.0) -> m_intersection + m_direction (1.0)
	//                                    (radial.first)                        (radial.second)
	double r;

	// Intersection position on line 2 : p1 (0.0) ->  p2 (1.0)
	//                              (edge.first)   (edge.second) 
	double s;

	interSection.first = Coord(0, 0, 0);
	interSection.second = IS_NONE;

	// check intersection here
	if(SegmentIntersection(radial.first, radial.first + radial.second, edge.first, edge.second, r, s))
	{
		// Segment intersects with v1
		if(s == 0 || (r > 0 && fabs(s) < 1e-6))
		{
			interSection.first = edge.first;
			interSection.second = IS_POINT;
			return true;
		}
		// Segment intersects with v2
		else if(s == 1 || (r > 0 && fabs(s - 1) < 1e-6))
		{
			interSection.first = edge.second;
			interSection.second = IS_POINT;
			return true;
		}
		// Segments don't intersect
		else if(r < -1e-6)
		{
			return false;
		}
		// Intersection ~ 0 : m_intersection point is on (v1; v2) line
		else if(r < 1e-6)
		{
			interSection.first = radial.first;
			interSection.second = IS_EDGE;
			return true;
		}
		// Intersection is on half-infinite line starting from m_intersection following m_direction
		else if(s > 0 && s < 1)
		{
			interSection.first = edge.first + (edge.second - edge.first) * s;
			interSection.second = IS_EDGE;
			return true;
		}
		// Intersection is on half-infinite line but outside [v1;v2] segment
		else
		{
			return false;
		}
	}

	return false;
}
