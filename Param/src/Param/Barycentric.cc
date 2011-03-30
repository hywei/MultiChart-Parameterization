#include "Barycentric.h"
#include <hj_3rd/hjlib/math/blas_lapack.h>
#include <hj_3rd/zjucad/matrix/lapack.h>
#include <hj_3rd/zjucad/matrix/io.h>

namespace PARAM
{
	std::vector<Coord2D> ComputeTriangleLocal2DCoord(const std::vector<Coord>& tri_vtx_coord)
	{
		assert(tri_vtx_coord.size() == 3);
		std::vector<Coord2D> local_coord(3);
		local_coord[0] = Coord2D(0, 0);

		Coord vec_1 = Coord(tri_vtx_coord[1] - tri_vtx_coord[0]);
		Coord vec_2 = Coord(tri_vtx_coord[2] - tri_vtx_coord[0]);

		double edge_len_1 = vec_1.abs();
		double edge_len_2 = vec_2.abs();
		double agl = angle(vec_1, vec_2);

		local_coord[1] = Coord2D(edge_len_1, 0);

		double x_3 = edge_len_2*cos(agl);
		double y_3 = edge_len_2*sin(agl);

		local_coord[2] = Coord2D(x_3, y_3);

		return local_coord;
	}

	Barycentrc ComputeVertexBarycentric(const std::vector<ParamCoord>& node_param_coord_array,
		const ParamCoord& vtx_param_coord)
	{
		zjucad::matrix::matrix<double> coord_matrix(3, 3);
		for(size_t c=0; c<3; ++c)
		{
			coord_matrix(0, c) = node_param_coord_array[c].s_coord;
			coord_matrix(1, c) = node_param_coord_array[c].t_coord;
			coord_matrix(2, c) = 1.0;
		}
		inv(coord_matrix);
		zjucad::matrix::matrix<double> right_matrix(3, 1);
		right_matrix(0, 0) = vtx_param_coord.s_coord;
		right_matrix(1, 0) = vtx_param_coord.t_coord;
		right_matrix(2, 0) = 1.0;

		zjucad::matrix::matrix<double> result_matrix = coord_matrix*right_matrix;

		return Coord(result_matrix(0, 0), result_matrix(1, 0), result_matrix(2, 0));
	}

	Barycentrc ComputeVertexBarycentric(const std::vector<Coord>& face_vert_coord_array,
		const Coord& vtx_coord)
	{
		zjucad::matrix::matrix<double> coord_matrix(3, 3);
		for(size_t r=0; r<3; ++r)
		{
			for(size_t c=0; c<3; ++c) coord_matrix(r, c) = face_vert_coord_array[c][r];
		}
		inv(coord_matrix);
		zjucad::matrix::matrix<double> right_matrix(3, 1);
		for(size_t r=0; r<3; ++r) right_matrix(r, 0) = vtx_coord[r];

		zjucad::matrix::matrix<double> res_matrix = coord_matrix*right_matrix;
		return Coord(res_matrix(0, 0), res_matrix(1, 0), res_matrix(2, 0));
	}

	Coord ComputeVertexCoord(const std::vector<Coord>& node_coord_array, const Barycentrc& baryc_coord)
	{
		
		zjucad::matrix::matrix<double> pos_matrix(3, 3);
		for(size_t i=0; i<3; ++i)
			for(size_t j=0; j<3; ++j) pos_matrix(i, j) = node_coord_array[j][i];
		zjucad::matrix::matrix<double> barycentric_matrix(3, 1);
		for(size_t i=0; i<3; ++i) barycentric_matrix(i, 0) = baryc_coord[i];

		zjucad::matrix::matrix<double> res_matrix = pos_matrix * barycentric_matrix;

		return Coord(res_matrix(0, 0), res_matrix(1, 0), res_matrix(2, 0));
	}

	bool IsValidBarycentic(const Barycentrc& baryc_coord)
	{
		for(size_t k=0; k<3; ++k) 
		{
			if(!(GreaterEqual(baryc_coord[k], 0, LARGE_ZERO_EPSILON) && LessEqual(
				baryc_coord[k], 1, LARGE_ZERO_EPSILON)))
				return false;
		}
		double by_sum = baryc_coord[0] + baryc_coord[1] + baryc_coord[2];
		if(fabs(by_sum - 1) < LARGE_ZERO_EPSILON) 
			return true;
		return false;
	}

	double ComputeErrorOutValidBarycentric(const Barycentrc& baryc_coord)
	{
		double error_value = 0.0;
		for(size_t k=0; k<3; ++k)
		{
			if(baryc_coord[k] < 0) 
			{
				error_value += baryc_coord[k]*baryc_coord[k];
			}else if(baryc_coord[k] > 1)
			{
				error_value += (baryc_coord[k] -1) * (baryc_coord[k] - 1);
			}
		}

		return error_value;
	}
} 
