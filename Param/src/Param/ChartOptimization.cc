#include "ChartOptimization.h"
#include "Parameter.h"
#include "ParamPatch.h"
#include "ChartCreator.h"
#include "../ModelMesh/MeshModel.h"
#include "quad_dof.h"

#include <hj_3rd/zjucad/matrix/matrix.h>
#include <hj_3rd/zjucad/matrix/lapack.h>
#include <hj_3rd/zjucad/matrix/io.h>
#include <hj_3rd/hjlib/math/blas_lapack.h>

using namespace std;

namespace PARAM{
    ChartOptimizor::ChartOptimizor(const Parameter& parameter, const ParamPatch& patch): m_parameter(parameter), m_patch(patch) {}

    ChartOptimizor::~ChartOptimizor(){}

    void ChartOptimizor::Optimization()
    {
        boost::shared_ptr<MeshModel> p_mesh = m_parameter.GetMeshModel();
        const vector<int>& face_index_vec = m_patch.m_face_index_array;

        const PolyIndexArray& mesh_face = p_mesh->m_Kernel.GetFaceInfo().GetIndex();

        m_dof_vec.clear(); m_dof_vec.resize(5, 0.0);
        
        zjucad::matrix::matrix<double> x(5, 1);
        zjucad::matrix::matrix<double> E_J(5, 1);
        zjucad::matrix::matrix<double> E_H(5, 5);

        zjucad::matrix::matrix<double> delta_x(5,1);

        SetNewtonMethodInitValue(x);
        
        int iter_times = 3;
        for(int i=0; i<iter_times; ++i){
            double E_total = 0;
            for(size_t k=0; k<face_index_vec.size(); ++k){
                int fid = face_index_vec[k];
                zjucad::matrix::matrix cur_c_mat = CalTempMatrix_1(fid) * CalTempMatrix_2(fid);
                double area = CalFaceAreaInParamDomain(fid);

                zjucad::matrix::matrix<double> cur_E_J(5,1);
                zjucad::matrix::matrix<double> cur_E_H(5, 5);
                calc_E_J(&cur_E_J[0], &x[0], &cur_c_mat);
                calc_E_H(&cur_E_H[0], &x[0], &cur_c_mat);

                E_J += cur_E_J * area;
                E_H += cur_E_H * area;

                zjucad::matrix::matrix<double> cur_E(1,1);
                calc_E(&cur_E[0], &x[0], &cur_c_mat);
                E_total += cur_E(0,0)*area;
            }

            SolveNewtonMethodIterValue(E_J, E_H, delta_x);
            UpdateNewtonMethodValue(x, delta_x);
        }
        
    }

    void ChartOptimizor::SetNewtonMethodInitValue(zjucad::matrix::matrix<double>& x) const{
        //! x is a 5x1 matrix
        //! {(0, 0), (1, 0), (1, 1), (0, 1)} = {(0, 0), (a, 0), (b, c), (d, e)}
        x(0,0) = x(1, 0) = x(2, 0) = x(4, 0) = 1; x(3,0) = 0;
    }

    void ChartOptimizor::SolveNewtonMethodIterValue(const zjucad::matrix::matrix<double>& j_mat, const zjucad::matrix::matrix<double>& h_mat, zjucad::matrix::matrix<double>& _x) const{
        zjucad::matrix::matrix<double> inv_h(h_mat);
        inv(inv_h);
        
        _x = inv_h * ( j_mat * (-1.0));

        
    }
    
    void ChartOptimizor::UpdateNewtonMethodValue(zjucad::matrix::matrix<double>& x, const zjucad::matrix::matrix<double>& delta_x) const{
        //! x and delta_x are both 5x1 matrix;
        x += delta_x;
    }
    
    zjucad::matrix::matrix<double> ChartOptimizor::CalTempMatrix_1(int fid)
    {
        boost::shared_ptr<MeshModel> p_mesh = m_parameter.GetMeshModel();
		const CoordArray& vert_coord_array = p_mesh->m_Kernel.GetVertexInfo().GetCoord();
		const PolyIndexArray& face_list_array = p_mesh->m_Kernel.GetFaceInfo().GetIndex();
		const IndexArray& faces = face_list_array[fid];
        
		std::vector<Coord> tri_vert_coord_3d(3);
		for(size_t k=0; k<3; ++k) tri_vert_coord_3d[k] = vert_coord_array[faces[k]];
		std::vector<Coord2D> local_coord = ComputeTriangleLocal2DCoord(tri_vert_coord_3d);
        zjucad::matrix::matrix<double> tmp_1(3, 3);
        for(int i=0; i<3; ++i){
            tmp_1(0, i) = local_coord[i][0];
            tmp_1(1, i) = local_coord[i][1];
            tmp_1(2, i) = 1;
        }

        inv(tmp_1);
        zjucad::matrix::matrix<double> tmp_2(3,2) = zjucad::matrix::eye<double>(3, 2);
        tmp_2(0,0) = tmp_2(1,1) = 1;

        return tmp_1*tmp_2;
    }

    zjucad::matrix::matrix<double> ChartOptimizor::CalTempMatrix_2(int fid)
    {
        const vector<ParamCoord>& vert_pc_vec = m_parameter.GetVertexParamCoordArray();
        boost::shared_ptr<MeshModel> p_mesh = m_parameter.GetMeshModel();
		const PolyIndexArray& face_list_array = p_mesh->m_Kernel.GetFaceInfo().GetIndex();
		const IndexArray& faces = face_list_array[fid];

        zjucad::matrix::matrix<double> tmp_mat(3, 3);
        for(size_t k=0; k<3; ++k){
            int vid = faces[k];
            const ParamCoord& pc = vert_pc_vec[vid];
            tmp_mat(k, 0) = pc.s_coord(1.0 - pc.t_coord);
            tmp_mat(k, 1) = pc.t_coord(1.0 - pc.s_coord);
            tmp_mat(k, 2) = pc.s_coord * pc.t_coord;
        }
        return tmp_mat;
    }

    double ChartOptimizor::CalFaceAreaInParamDomain(int fid)
    {
const vector<ParamCoord>& vert_pc_vec = m_parameter.GetVertexParamCoordArray();
        boost::shared_ptr<MeshModel> p_mesh = m_parameter.GetMeshModel();
		const PolyIndexArray& face_list_array = p_mesh->m_Kernel.GetFaceInfo().GetIndex();
		const IndexArray& faces = face_list_array[fid];

        double u[3], v[3];
        for(size_t i=0; i<3; ++i) {
            int vid = faces[i];
            u[i] = vert_pc_vec[vid].s_coord;
            v[i] = vert_pc_vec[vid].t_coord;
        }
        
        return fabs((u[1] - u[0]) * (v[2] - v[0]) - (u[2] - u[0]) * (v[1] - v[0]));
    }
}
