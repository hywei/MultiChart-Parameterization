#include "ParamDistortion.h"

#include "../hj_3rd/include/math/blas_lapack.h"
#include "../hj_3rd/include/zjucad/matrix/lapack.h"
#include "../hj_3rd/include/zjucad/matrix/io.h"
#include <limits>
using namespace std;

namespace PARAM
{
	const double ParamDistortion::EPSILON = 1e-8;

	ParamDistortion::ParamDistortion(const PARAM::QuadParam& param) 
		: m_param(param), 
		p_mesh(param.GetMeshModel()),
		m_param_coord(param.GetParamCoord())
	{
	}

	ParamDistortion::~ParamDistortion(){}

	void ParamDistortion::ComputeJacobiMatrix()
	{
		size_t face_num = p_mesh->m_Kernel.GetModelInfo().GetFaceNum();

		m_jacobi_matrix.clear();
		m_jacobi_matrix.resize(face_num);
		for(size_t k=0; k<face_num; ++k)
		{
			//m_jacobi_matrix[k] = ComputeLocalJacobiMatrix(k);
			m_jacobi_matrix[k] = ComputeTriJacobiMatrix(k);
		}
	}

	zjucad::matrix::matrix<double> ParamDistortion::ComputeLocalJacobiMatrix(int fid)
	{
		const CoordArray& vtx_coord_array = p_mesh->m_Kernel.GetVertexInfo().GetCoord();
		const PolyIndexArray& face_index_array = p_mesh->m_Kernel.GetFaceInfo().GetIndex();
		const IndexArray& face_index = face_index_array[fid];

		vector<ParamCoord> inChart_param_coord(3);
		int face_chart_id = m_param.GetFaceChartID(fid);
		for(size_t k=0; k<3; ++k)
		{
			int vtx_idx = face_index[k];
			int vtx_chart_id = m_param.GetVertexChartID(vtx_idx);
			if(vtx_chart_id != face_chart_id)
			{
				m_param.TransParamCoordBetweenTwoChart(vtx_chart_id, face_chart_id,
					m_param_coord[vtx_idx], inChart_param_coord[k]);
			}else
			{
				inChart_param_coord[k] = m_param_coord[vtx_idx];
			}
		}

		vector<Coord> vtx_coord(3);
		for(size_t k=0; k<3; ++k)
		{
			int vtx_idx = face_index[k];
			vtx_coord[k] = vtx_coord_array[vtx_idx];
		}

		return ComputeLocalJacobiMatrix(vtx_coord, inChart_param_coord);
	}

	zjucad::matrix::matrix<double> ParamDistortion::ComputeLocalJacobiMatrix(const std::vector<Coord>& vtx_coord, 
		const std::vector<ParamCoord>& param_coord_array)
	{
		
		zjucad::matrix::matrix<double> uv_matrix(2, 3);
		for(size_t c=0; c<3; ++c) 
		{
			uv_matrix(0, c) = param_coord_array[c].s_coord;
			uv_matrix(1, c) = param_coord_array[c].t_coord;
		}
	 		
		vector<Coord2D> tri_local_coord = ComputeLocal2DCoord(vtx_coord);
		zjucad::matrix::matrix<double> mid_matrix(3, 3);
		for(size_t c=0; c<3; ++c)
		{
			mid_matrix(0, c) = tri_local_coord[c][0];
			mid_matrix(1, c) = tri_local_coord[c][1];
			mid_matrix(2, c) = 1.0;
		}
		inv(mid_matrix);

		zjucad::matrix::matrix<double> right_matrix(3,2);
		right_matrix(0, 0) = 1.0; right_matrix(0, 1) = 0.0;
		right_matrix(1, 0) = 0.0; right_matrix(1, 1) = 1.0;
		right_matrix(2, 0) = 0.0; right_matrix(2, 1) = 0.0;

		zjucad::matrix::matrix<double> temp_matrix = uv_matrix * mid_matrix;
		zjucad::matrix::matrix<double> ret_matrix = temp_matrix*right_matrix;

		return ret_matrix;
	}

	zjucad::matrix::matrix<double> ParamDistortion::ComputeTriJacobiMatrix(int fid)
	{		
		const CoordArray& vtx_coord_array = p_mesh->m_Kernel.GetVertexInfo().GetCoord();
		const PolyIndexArray& face_index_array = p_mesh->m_Kernel.GetFaceInfo().GetIndex();
		const IndexArray& face_index = face_index_array[fid];

		vector<ParamCoord> inChart_param_coord(3);
		//int face_chart_id = m_param.GetFaceChartID(fid);
		int face_chart_id = m_param.GetVertexChartID(face_index[0]);
		for(size_t k=0; k<3; ++k)
		{
			int vtx_idx = face_index[k];
			int vtx_chart_id = m_param.GetVertexChartID(vtx_idx);
			if(vtx_chart_id != face_chart_id)
			{
				m_param.TransParamCoordBetweenTwoChart(vtx_chart_id, face_chart_id,
					m_param_coord[vtx_idx], inChart_param_coord[k]);
			}else
			{
				inChart_param_coord[k] = m_param_coord[vtx_idx];
			}
		}

		vector<Coord> vtx_coord(3);
		for(size_t k=0; k<3; ++k)
		{
			int vtx_idx = face_index[k];
			vtx_coord[k] = vtx_coord_array[vtx_idx];
		}
		return ComputeTriJacobiMatrix(vtx_coord, inChart_param_coord);
	}

	zjucad::matrix::matrix<double> ParamDistortion::ComputeTriJacobiMatrix(const std::vector<Coord>& tri_face,
		const std::vector<ParamCoord>& param_coord_array)
	{
		double u[3], v[3];
		for(size_t i=0; i<3; ++i) { u[i] = param_coord_array[i].s_coord; v[i] = param_coord_array[i].t_coord;}

		double area_2 = fabs( (u[1] - u[0]) * (v[2] - v[0]) - (u[2] - u[0]) * (v[1] - v[0]));
		double area_2_inv = 1.0/area_2;
		vector<Coord2D> local_coord = ComputeLocal2DCoord(tri_face);

// 		cout<<"local coord: "<<endl;
// 		for(int k=0; k<3; ++k)
// 		{
// 			cout<< local_coord[k][0] <<' ' <<local_coord[k][1]<<endl;
// 		}
// 		cout<<"param coord: "<<endl;
// 		for(int k=0; k<3; ++k)
// 		{
// 			cout<< param_coord_array[k].s_coord << " " << param_coord_array[k].t_coord<<endl;
// 		}
// 
// 		cout<< "area" << area_2 << endl;
// 		cout<< v[1] - v[2] << " " << v[2] - v[0] << " " << v[0] - v[1] << endl;
// 		cout<< u[2] - u[1] << " " << u[0] - u[2] << " " << u[1] - u[0] << endl;

		Coord2D partial_u = (local_coord[0]*(v[1] - v[2]) + local_coord[1]*(v[2] - v[0]) + local_coord[2]*(v[0] - v[1]))*area_2_inv;
		Coord2D partial_v = (local_coord[0]*(u[2] - u[1]) + local_coord[1]*(u[0] - u[2]) + local_coord[2]*(u[1] - u[0]))*area_2_inv;

		zjucad::matrix::matrix<double> j_mat(2, 2);
		j_mat(0, 0) = partial_u[0]; j_mat(0, 1) = partial_v[0];
		j_mat(1, 0) = partial_u[1]; j_mat(1, 1) = partial_v[1];


		return j_mat;
	}

	zjucad::matrix::matrix<double> ParamDistortion::ComputeTriJacobiMatrix(const std::vector<Coord2D>& local_coord,
		const std::vector<ParamCoord>& param_coord_array)
	{
//		cout<<"local coord: "<<endl;
//		for(int k=0; k<3; ++k)
//		{
//			cout<< local_coord[k][0] <<' ' <<local_coord[k][1]<<endl;
//		}
//		cout<<"param coord: "<<endl;
//		for(int k=0; k<3; ++k)
//		{
//			cout<< param_coord_array[k].s_coord << " " << param_coord_array[k].t_coord<<endl;
//		}
		double u[3], v[3];
		for(size_t i=0; i<3; ++i) { u[i] = param_coord_array[i].s_coord; v[i] = param_coord_array[i].t_coord;}

		double area_2 = fabs( (u[1] - u[0]) * (v[2] - v[0]) - (u[2] - u[0]) * (v[1] - v[0]));
		if(fabs(area_2) < EPSILON) 
		{
			printf("Warnning, degenerate triangle in paramaterzation domain!\n");
		}
		double area_2_inv = 1.0/area_2;

//		cout<< "area" << area_2 << endl;
//		cout<< v[1] - v[2] << " " << v[2] - v[0] << " " << v[0] - v[1] << endl;
//		cout<< u[2] - u[1] << " " << u[0] - u[2] << " " << u[1] - u[0] << endl;

		Coord2D partial_u = (local_coord[0]*(v[1] - v[2]) + local_coord[1]*(v[2] - v[0]) + local_coord[2]*(v[0] - v[1]))*area_2_inv;
		Coord2D partial_v = (local_coord[0]*(u[2] - u[1]) + local_coord[1]*(u[0] - u[2]) + local_coord[2]*(u[1] - u[0]))*area_2_inv;

		zjucad::matrix::matrix<double> j_mat(2, 2);
		j_mat(0, 0) = partial_u[0]; j_mat(0, 1) = partial_v[0];
		j_mat(1, 0) = partial_u[1]; j_mat(1, 1) = partial_v[1];
		
		return j_mat;
	}

	vector<Coord2D> ParamDistortion::ComputeLocal2DCoord(const vector<Coord>& tri_vtx_coord) 
	{
		assert(tri_vtx_coord.size() == 3);
		vector<Coord2D> local_coord(3);
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


	bool ParamDistortion::ComputeSurfaceCoord(const ChartParamCoord& chart_param_coord, SurfaceCoord& surface_coord)
	{
		const PolyIndexArray& face_index_array = p_mesh->m_Kernel.GetFaceInfo().GetIndex();
		int vert_num = p_mesh->m_Kernel.GetModelInfo().GetVertexNum();
		int face_num = p_mesh->m_Kernel.GetModelInfo().GetFaceNum();

		std::vector< int > valid_face_array;
		std::vector< Coord > valid_barycentric_array;

		double min_out_bycentric_error = numeric_limits<double>::infinity();

		int valid_barycentric_num = 0;
		for(int k=0; k<face_num;  ++k)
		{
			const IndexArray& face_index = face_index_array[k];
			vector<ParamCoord> param_coord_array(3);
			for(int i=0; i<3; ++i)
			{
				int vid = face_index[i];				
				ParamCoord param_coord = m_param_coord[vid];
				param_coord_array[i] = param_coord;
				int cur_chart_id = m_param.GetVertexChartID(vid);
				if(cur_chart_id != chart_param_coord.chart_id)
				{
					m_param.TransParamCoordBetweenTwoChart(cur_chart_id, chart_param_coord.chart_id, 
						param_coord, param_coord_array[i]);
				}

			}
			Coord barycenter_coord = ComputeBarycentric(param_coord_array,
				ParamCoord(chart_param_coord.s_coord, chart_param_coord.t_coord));
			if(IsValidBarycentric(barycenter_coord))
			{
				surface_coord = SurfaceCoord(k, barycenter_coord);
				valid_barycentric_num++;

				valid_face_array.push_back(k);
				valid_barycentric_array.push_back(barycenter_coord);
				
				return true;
			}else
			{
				if(valid_barycentric_num != 0) continue;
				double cur_out_bycentric_error = CalOutBarycentricError(barycenter_coord);
				if(cur_out_bycentric_error < min_out_bycentric_error)
				{
					min_out_bycentric_error = cur_out_bycentric_error;
					surface_coord = SurfaceCoord(k, barycenter_coord);
				}
			}
		}

		if(valid_barycentric_num > 1)
		{
			printf("There are %d valid barycentric in surface!\n", valid_barycentric_num);

			for(int k=0; k<valid_barycentric_num; ++k)
			{

				Coord pos_coord(0, 0, 0);
				const CoordArray& vtx_coord_array = p_mesh->m_Kernel.GetVertexInfo().GetCoord();
				const IndexArray& face_index = face_index_array[valid_face_array[k]];

				vector<Coord> vtx_coord(3);
				for(size_t i=0; i<3; ++i) vtx_coord[i] = vtx_coord_array[face_index[i]];
				Coord pos = ComputePosWithBarycentric(vtx_coord, valid_barycentric_array[k]);

				printf("Coordinate in surface 1 : %lf %lf %lf\n", pos[0], pos[1], pos[2]);
			}
		}
		if(valid_barycentric_num !=0) return true;
		
		return false;
	}

	Coord ParamDistortion::ComputeBarycentric(const vector<ParamCoord>& tri_vtx_array, ParamCoord vtx_param_coord)
	{
		zjucad::matrix::matrix<double> coord_matrix(3, 3);
		for(size_t c=0; c<3; ++c)
		{
			coord_matrix(0, c) = tri_vtx_array[c].s_coord;
			coord_matrix(1, c) = tri_vtx_array[c].t_coord;
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

	Coord ParamDistortion::ComputePosWithBarycentric(const std::vector<Coord>& vtx_coord_array, Coord barycentric)
	{
		assert(vtx_coord_array.size() == 3);
		zjucad::matrix::matrix<double> pos_matrix(3, 3);
		for(size_t i=0; i<3; ++i)
			for(size_t j=0; j<3; ++j) pos_matrix(i, j) = vtx_coord_array[j][i];
		zjucad::matrix::matrix<double> barycentric_matrix(3, 1);
		for(size_t i=0; i<3; ++i) barycentric_matrix(i, 0) = barycentric[i];

		zjucad::matrix::matrix<double> res_matrix = pos_matrix * barycentric_matrix;

		return Coord(res_matrix(0, 0), res_matrix(1, 0), res_matrix(2, 0));
	}

	double ParamDistortion::CalOutBarycentricError(const Coord& barycentric)
	{
		double error_value = 0.0;
		for(size_t k=0; k<3; ++k)
		{
			if(barycentric[k] < 0) 
			{
				error_value += barycentric[k]*barycentric[k];
			}else if(barycentric[k] > 1)
			{
				error_value += (barycentric[k] -1) * (barycentric[k] - 1);
			}
		}

		return error_value;
	}
}