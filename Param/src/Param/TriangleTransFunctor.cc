#include "TriangleTransFunctor.h"
#include "Parameterization.h"
#include "Barycentric.h"

#include <hj_3rd/hjlib/math/blas_lapack.h>
#include <hj_3rd/zjucad/matrix/lapack.h>


#include <queue>
#include <set>
#include <map>
using namespace std;

namespace PARAM
{
	TriTransFunctor::TriTransFunctor(boost::shared_ptr<MeshModel> _mesh) 
		: p_mesh(_mesh)
	{
		SetChartNeighbors();
	}

	TriTransFunctor::~TriTransFunctor(){}

	bool TriTransFunctor::TransLocalCoordBetweenCharts(int from_chart_id, const Coord2D& from_chart_coord, 
		int to_chart_id, Coord2D& to_chart_coord) const
	{
		vector<int> trans_list; 
		if(! GetTransChartsChain(from_chart_id, to_chart_id, trans_list))
		{
			return false;
		}

		//printf("Transition chain length: %d %d %d\n", from_chart_id, to_chart_id, trans_list.size());
		if(trans_list.size() == 1) to_chart_coord = from_chart_coord;
		else
		{			
			zjucad::matrix::matrix<double> trans_mat(zjucad::matrix::eye<double>(3));		
	   
			for(size_t k=1; k<trans_list.size(); ++k)
			{			   
				zjucad::matrix::matrix<double> temp_trans_mat_1 = 
					GetTransMatrixOfAdjCharts(trans_list[k-1], trans_list[k]) ;
				zjucad::matrix::matrix<double> temp_trans_mat_2 = trans_mat;
				trans_mat = temp_trans_mat_1 * temp_trans_mat_2;
			
			}

			zjucad::matrix::matrix<double> origin_coord (3, 1);
			origin_coord(0, 0) = from_chart_coord[0]; 
			origin_coord(1, 0) = from_chart_coord[1];
			origin_coord(2, 0) = 1.0;
  
			zjucad::matrix::matrix<double> new_coord = trans_mat * origin_coord;

			to_chart_coord = Coord2D(new_coord(0, 0), new_coord(1, 0));
		
		}		

		return true;
	}

	zjucad::matrix::matrix<double> TriTransFunctor::GetTransMatrixOfAdjCharts(int from_chart_id, int to_chart_id) const
	{
		const PolyIndexArray& face_list = p_mesh->m_Kernel.GetFaceInfo().GetIndex();
		const IndexArray& face_1 = face_list[from_chart_id];
		const IndexArray& face_2 = face_list[to_chart_id];

		/// find the common edge of these two faces
		pair<int, int> com_edge_idx_1, com_edge_idx_2;
		for(int k=0; k<3; ++k)
		{			
			bool flag = false;
			com_edge_idx_1 = make_pair(k, (k+1)%3);
			pair<int, int> com_edge_1 = make_pair(face_1[k], face_1[(k+1)%3]);
			for(int i=0; i<3; ++i)
			{
				com_edge_idx_2 = make_pair(i, (i+1)%3);
				pair<int, int> com_edge_2 = make_pair(face_2[i], face_2[(i+1)%3]);
				if(com_edge_1 == com_edge_2 || 
					(com_edge_1.first == com_edge_2.second && com_edge_1.second == com_edge_2.first)) { flag = true; break; }
			}
			if(flag) break;
		}

		/// default:  each face's first vertex is the origin, the frist edge is the local frame's x
		zjucad::matrix::matrix<double> trans_mat_1 = GetTransMatrixInOneChart(from_chart_id,
			make_pair(0, 1), make_pair(com_edge_idx_1.first, com_edge_idx_1.second));
		zjucad::matrix::matrix<double> trans_mat_2 = GetTransMatrixInOneChart(to_chart_id,
			make_pair(0, 1), make_pair(com_edge_idx_2.second, com_edge_idx_2.first));
		inv(trans_mat_2);

		return trans_mat_2 * trans_mat_1;

	}

	zjucad::matrix::matrix<double> TriTransFunctor::GetTransMatrixInOneChart(int chart_id, 
		std::pair<int, int> old_x_index, std::pair<int, int> new_x_index) const
	{
		zjucad::matrix::matrix<double> t_mat( zjucad::matrix::eye<double>(3));
		zjucad::matrix::matrix<double> r_mat( zjucad::matrix::eye<double>(3));

		/// Get this triangle's local 2d frame
		const CoordArray& vtx_coord_array = p_mesh->m_Kernel.GetVertexInfo().GetCoord();
		const IndexArray& face = (p_mesh->m_Kernel.GetFaceInfo().GetIndex())[chart_id];
		vector<Coord> tri_vtx_coord_array(3);
		for(int k=0; k<3; ++k) tri_vtx_coord_array[k] = vtx_coord_array[face[k]];
		vector<Coord2D> local_coord_array = ComputeTriangleLocal2DCoord(tri_vtx_coord_array);

		int old_origin_idx = old_x_index.first;
		int new_origin_idx = new_x_index.first;
		
		if(old_origin_idx != new_origin_idx)
		{
			double t_x = local_coord_array[old_origin_idx][0] - local_coord_array[new_origin_idx][0];
			double t_y = local_coord_array[old_origin_idx][1] - local_coord_array[new_origin_idx][1];
			t_mat = GetTranslationMatrix2D(t_x, t_y);
		}

		if(old_x_index != new_x_index)
		{
 			Coord2D old_x_vec = local_coord_array[old_x_index.second] - local_coord_array[old_x_index.first];
 			Coord2D new_x_vec = local_coord_array[new_x_index.second] - local_coord_array[new_x_index.first];

			Coord2D tmp_a_vec (old_x_vec), tmp_b_vec (new_x_vec);
			tmp_a_vec.normalize(); tmp_b_vec.normalize();
			double cross_v = tmp_a_vec[0] * tmp_b_vec[1] - tmp_a_vec[1] * tmp_b_vec[0];
			double dot_v = tmp_a_vec[0] * tmp_b_vec[0] + tmp_a_vec[1] * tmp_b_vec[1];


			double asin_angle = asin(cross_v);
			double acos_angle = acos(dot_v);

			double r_angle;
			if(fabs( cross_v ) < 1e-12) 
			{
				if(dot_v < 0) r_angle = PI;
				else r_angle = 0;
			}else if (fabs( dot_v) < 1e-12)
			{
				if(cross_v > 0) r_angle = PI/2;
				else r_angle = PI*3 /2;
			}else 
			{
				if(cross_v < 0 && dot_v > 0) r_angle = asin_angle + 2*PI; /// 4 
				else if( cross_v < 0 && dot_v < 0) r_angle = PI - asin_angle; /// 3
				else if( cross_v > 0 && dot_v < 0) r_angle = acos_angle; /// 2
				else if( cross_v > 0 && dot_v > 0) r_angle = asin_angle; /// 1

				if(r_angle < 0) r_angle += 2*PI;
			}
			//printf("angle is %lf\n", r_angle*180/PI);		   		
			
			r_mat = GetRotationMatrix2D(r_angle);
		}

		return r_mat*t_mat;

	}

	bool TriTransFunctor::GetTransChartsChain(int from_chart_id, int to_chart_id, std::vector<int>& trans_list) const
	{
		trans_list.clear();

		std::set<int> visited_face_set;
		std::map<int, int> prev_chart_map;

		queue<int> q;
		q.push(to_chart_id);
		visited_face_set.insert(to_chart_id);
		prev_chart_map[to_chart_id] = -1;


		bool flag = false;
		while(!q.empty())
		{
			int cur_chart_id = q.front(); q.pop();
			if(cur_chart_id == from_chart_id) { flag = true; break;}

			const vector<int>& adj_charts = m_neighbor_charts[cur_chart_id];
			for(size_t k=0; k<adj_charts.size(); ++k)
			{
				int chart_id = adj_charts[k];
				if(visited_face_set.find(chart_id) == visited_face_set.end())
				{
					q.push(chart_id);
					visited_face_set.insert(chart_id);
					prev_chart_map[chart_id] = cur_chart_id;
				}
			}
		}

		if(!flag) return false;
		
		trans_list.push_back(from_chart_id);
		int prev_chart = prev_chart_map[from_chart_id];

		while(prev_chart != -1)
		{
			trans_list.push_back(prev_chart);
			prev_chart = prev_chart_map[prev_chart];
		}

		return true;
	}


	zjucad::matrix::matrix<double> TriTransFunctor::GetTranslationMatrix2D(double t_x, double t_y) const
	{
		zjucad::matrix::matrix<double> t_mat( zjucad::matrix::eye<double>(3));
		t_mat(0, 2) = t_x;
		t_mat(1, 2) = t_y;

		return t_mat;
	}
	
	zjucad::matrix::matrix<double> TriTransFunctor::GetRotationMatrix2D(double rotation_angle) const 
	{
		zjucad::matrix::matrix<double> r_mat(3, 3);
		double cos_v = cos(rotation_angle), sin_v = sin(rotation_angle);
		r_mat(0, 0) = cos_v;  r_mat(0, 1) = sin_v; r_mat(0, 2) = 0;
		r_mat(1, 0) = -sin_v; r_mat(1, 1) = cos_v; r_mat(1, 2) = 0;
		r_mat(2, 0) = 0; r_mat(2, 1) = 0; r_mat(2, 2) = 1;
		
		return r_mat;
	}

	void TriTransFunctor::SetChartNeighbors()
	{
		const PolyIndexArray& face_list_array = p_mesh->m_Kernel.GetFaceInfo().GetIndex();
		int face_num = p_mesh->m_Kernel.GetModelInfo().GetFaceNum();

		m_neighbor_charts.clear();
		m_neighbor_charts.resize(face_num);
	    
		map< pair<int, int>, vector<int> > edge_adj_face_map;

		for(int fid=0; fid < face_num; ++fid)
		{
			const IndexArray& vertices = face_list_array[fid];
			for(int i=0; i<3; ++i)
			{
				int vid1 = vertices[i];
				int vid2 = vertices[(i+1)%3];
				pair<int, int> edge = (vid1 > vid2) ? make_pair(vid2 , vid1) : make_pair(vid1, vid2);
				edge_adj_face_map[edge].push_back(fid);
			}
		}

		for(int fid=0; fid < face_num; ++fid)
		{
			const IndexArray& vertices = face_list_array[fid];
			m_neighbor_charts[fid].clear(); 			
			for(int i=0; i<3; ++i)
			{
				int vid1 = vertices[i];
				int vid2 = vertices[(i+1)%3];
				pair<int, int> edge = (vid1 > vid2) ? make_pair(vid2 , vid1) : make_pair(vid1, vid2);
				const vector<int>& adj_faces = edge_adj_face_map[edge];
				if(adj_faces.size() == 2)
				{
					int _fid = (fid == adj_faces[0]) ? adj_faces[1] : adj_faces[0];
					m_neighbor_charts[fid].push_back(_fid);
				}
			}
		}
	}

	void TriTransFunctor::TestInnerChartTrans(int chart_id) const
	{
		printf("Test One Chart Inner Transition function, chart id : %d\n", chart_id);
		zjucad::matrix::matrix<double> trans_mat = GetTransMatrixInOneChart( chart_id, 
			make_pair(0, 1), make_pair(1, 2));

		zjucad::matrix::matrix<double> origin_coord(3, 1);  

		origin_coord(0, 0) = 0; origin_coord(1, 0) = 0; origin_coord(2, 0) = 1;
		zjucad::matrix::matrix<double> new_coord = trans_mat * origin_coord;
		printf("new coord: %lf %lf %lf\n", new_coord(0, 0), new_coord(1, 0), new_coord(2, 0));

		trans_mat = GetTransMatrixInOneChart( chart_id, make_pair(0, 1), make_pair(2, 0));
		new_coord = trans_mat * origin_coord;
		printf("new coord: %lf %lf %lf\n", new_coord(0, 0), new_coord(1, 0), new_coord(2, 0));

		zjucad::matrix::matrix<double> old_coord(3, 1);
		old_coord(0, 0) = 0.13333; old_coord(1, 0) = 0.13333; old_coord(2, 0) = 1;
		new_coord = trans_mat * old_coord;
		printf("new coord: %lf %lf %lf\n", new_coord(0, 0), new_coord(1, 0), new_coord(2, 0));

		trans_mat = GetTransMatrixInOneChart( chart_id, make_pair(0, 1), make_pair(2, 1));
		new_coord = trans_mat * origin_coord;
		printf("new coord: %lf %lf %lf\n", new_coord(0, 0), new_coord(1, 0), new_coord(2, 0));
		new_coord = trans_mat * old_coord;
		printf("new coord: %lf %lf %lf\n", new_coord(0, 0), new_coord(1, 0), new_coord(2, 0));

		trans_mat = GetTransMatrixInOneChart( chart_id, make_pair(0, 1), make_pair(0, 2));
		new_coord = trans_mat * origin_coord;
		printf("new coord: %lf %lf %lf\n", new_coord(0, 0), new_coord(1, 0), new_coord(2, 0));
		new_coord = trans_mat * old_coord;
		printf("new coord: %lf %lf %lf\n", new_coord(0, 0), new_coord(1, 0), new_coord(2, 0));

	}

	void TriTransFunctor::TestAdjChartsTrans(int from_chart_id, int to_chart_id) const
	{
		printf("Test Adjacent Charts Transition Function, chart id: %d %d\n", from_chart_id, to_chart_id);
		zjucad::matrix::matrix<double> trans_mat = GetTransMatrixOfAdjCharts(from_chart_id, to_chart_id);

		zjucad::matrix::matrix<double> origin_coord(3, 1);  
		origin_coord(0, 0) = 0; origin_coord(1, 0) = 0; origin_coord(2, 0) = 1;
		zjucad::matrix::matrix<double> new_coord = trans_mat * origin_coord;
		printf("new coord: %lf %lf %lf\n", new_coord(0, 0), new_coord(1, 0), new_coord(2, 0));

		zjucad::matrix::matrix<double> old_coord(3, 1);
		old_coord(0, 0) = 0.13333; old_coord(1, 0) = 0.13333; old_coord(2, 0) = 1;
		new_coord = trans_mat * old_coord;
		printf("new coord: %lf %lf %lf\n", new_coord(0, 0), new_coord(1, 0), new_coord(2, 0));

		old_coord(0, 0) = 0.13333; old_coord(1, 0) = 0; old_coord(2, 0) = 1.0;
		new_coord = trans_mat * old_coord;
		printf("new coord: %lf %lf %lf\n", new_coord(0, 0), new_coord(1, 0), new_coord(2, 0));
	}

	void TriTransFunctor::TestTwoChartsTrans(int from_chart_id, int to_chart_id) const
	{
		printf("Test Two Charts Transition Function, chart id: %d %d\n", from_chart_id, to_chart_id);
		Coord2D old_coord(0, 0), new_coord;
		TransLocalCoordBetweenCharts(from_chart_id, old_coord, to_chart_id, new_coord);
		printf("new coord: %lf %lf\n", new_coord[0], new_coord[1]);
		old_coord = Coord2D(1, 1);
		TransLocalCoordBetweenCharts(from_chart_id, old_coord, to_chart_id, new_coord);
		printf("new coord: %lf %lf\n", new_coord[0], new_coord[1]);

		old_coord = Coord2D(0.133333, 0.133333);
		TransLocalCoordBetweenCharts(from_chart_id, old_coord, to_chart_id, new_coord);
		printf("new coord: %lf %lf\n", new_coord[0], new_coord[1]);

		old_coord = Coord2D(0.2, 0.1);
		TransLocalCoordBetweenCharts(from_chart_id, old_coord, to_chart_id, new_coord);
		printf("new coord: %lf %lf\n", new_coord[0], new_coord[1]);
	}
}
