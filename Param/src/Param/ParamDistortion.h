#ifndef PARAMDISTORTION_H_
#define PARAMDISTORTION_H_

#include "../ModelMesh/MeshModel.h"
#include "../hj_3rd/include/zjucad/matrix/matrix.h"

#include "QuadParam.h"

#include <vector>
#include <boost/shared_ptr.hpp>

namespace PARAM
{
	class ParamDistortion
	{
	public:
		ParamDistortion(const PARAM::QuadParam& param);
		~ParamDistortion();

		zjucad::matrix::matrix<double> ComputeLocalJacobiMatrix(int fid);
		zjucad::matrix::matrix<double> ComputeLocalJacobiMatrix(const std::vector<Coord>& tri_face, 
			const std::vector<ParamCoord>& param_coord_array);		

		zjucad::matrix::matrix<double> ComputeTriJacobiMatrix(int fid);
		zjucad::matrix::matrix<double> ComputeTriJacobiMatrix(const std::vector<Coord>& tri_face,
			const std::vector<ParamCoord>& param_coord_array);
		zjucad::matrix::matrix<double> ComputeTriJacobiMatrix(const std::vector<Coord2D>& tri_face_2d,
			const std::vector<ParamCoord>& param_coord_array);

		//! From a parameter coordinate (u, v , chart_id), try to find its coordinate in surface
		//! this function only could be call after computed parameterization
		/*
		\param chart_param_coord the parameter coordinate ( u, v, chart_id)
		\param surface_coord the parameter coordinate corresponding surface coordinate
		*/
		bool ComputeSurfaceCoord(const ChartParamCoord& chart_param_coord, SurfaceCoord& surface_coord);

	private:
		void ComputeJacobiMatrix();

	public:
		static std::vector<Coord2D> ComputeLocal2DCoord(const std::vector<Coord>& tri_vtx_coord);
  	   
		/// method for 
		static Coord ComputeBarycentric(const std::vector<ParamCoord>& tri_vtx_array, ParamCoord vtx_coord);
		static Coord ComputePosWithBarycentric(const std::vector<Coord>& vtx_coord_array, Coord barycentric);

	    static double CalOutBarycentricError(const Coord& barycentric);

		static bool IsValidBarycentric(const Coord& barycenter)
		{ 
			for(size_t k=0; k<3; ++k) 
			{
				if(!(GreaterEqual(barycenter[k], 0, EPSILON) && LessEqual(barycenter[k], 1, EPSILON)))
					return false;
			}
			double by_sum = barycenter[0] + barycenter[1] + barycenter[2];
			if(fabs(by_sum - 1) < EPSILON) 
				return true;
			return false;
		}

	private:
		static bool GreaterEqual(double a, double b, double epsilon)
		{
			return (a > b) || ( fabs(a - b) < epsilon);
		}
		static bool LessEqual(double a, double b, double epsilon)
		{
			return (a < b) || ( fabs(a - b) < epsilon);
		}


	private:
		const PARAM::QuadParam& m_param;

		const boost::shared_ptr<MeshModel> p_mesh;
		const std::vector<PARAM::ParamCoord>& m_param_coord;		

		std::vector<zjucad::matrix::matrix<double> > m_jacobi_matrix; /// each element is a 3x3 = 9 array


	public:
		const static double EPSILON;

	};
}

#endif