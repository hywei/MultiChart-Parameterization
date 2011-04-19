#include "ChartOptimization.h"
#include "Parameter.h"
#include "ParamPatch.h"
#include "ChartCreator.h"
extern "C" {
#include "quad_dof.h"
}
#include "Barycentric.h"

#include "../ModelMesh/MeshModel.h"

#include <hj_3rd/hjlib/math/blas_lapack.h>
#include <hj_3rd/zjucad/matrix/lapack.h>
#include <hj_3rd/zjucad/matrix/io.h>


using namespace std;
using namespace zjucad::matrix;

void energy_func_grad(const alglib::real_1d_array& x, double& func, alglib::real_1d_array& grad, void* ptr)
{
    PARAM::ChartOptimizor* cop = reinterpret_cast<PARAM::ChartOptimizor*> (ptr);
    assert(cop != 0);
    cop->EnergyFuncGrad(x, func, grad, ptr);
}

namespace PARAM{
    ChartOptimizor::ChartOptimizor(const Parameter& parameter, const ParamPatch& patch): m_parameter(parameter), m_patch(patch) {}

    ChartOptimizor::~ChartOptimizor(){}

    void ChartOptimizor::Optimization()
    {
        m_dof_vec.clear();
        m_dof_vec.resize(5, 0.0);
        
        m_dof_vec = SolveWithLbfgs();
    }
    
    matrix<double> ChartOptimizor::CalInitialValue() const
    {
        matrix<double> x(5, 1);        
        const vector<ParamCoord>& corner_pc_vec = m_patch.m_conner_pc_array;
        
        x(0,0) = corner_pc_vec[1].s_coord;
        x(1,0) = corner_pc_vec[2].s_coord;
        x(2,0) = corner_pc_vec[2].t_coord;

        x(3,0) = corner_pc_vec[3].s_coord;
        x(4,0) = corner_pc_vec[3].t_coord;
        return x;
    }

    vector<double> ChartOptimizor::SolveWithLbfgs() const
    {
        alglib::real_1d_array init_value = "[1,1,1,0,1]";
        double epsg = 1e-8;
        double epsf = 0;
        double epsx = 0;
        alglib::ae_int_t max_iters = 0;
        
        alglib::minlbfgsstate state;
        alglib::minlbfgsreport report;

        alglib::minlbfgscreate(5, init_value, state); //! create lbfgs solver
        alglib::minlbfgssetcond(state, epsg, epsf, epsx, max_iters);
        
        //alglib::minlbfgsoptimize(state, EnergyFuncGrad);
        double func_value;
        alglib::real_1d_array grad_vec;
        const void* ptr = reinterpret_cast<const void*> (this);
        assert(ptr != 0);
        alglib::minlbfgsoptimize(state, energy_func_grad, NULL, const_cast<void*> (ptr));
        //        alglib::minlbfgsoptimize(state, energy_func_grad);
        alglib::minlbfgsresults(state, init_value, report);

        cout << "Lbfgs Result: " << endl;
        //cout << "\t" << int(report.terminationtype) << endl;
        cout << "\t" << init_value.tostring(5).c_str() << endl;

        vector<double> ret(5);
        for(size_t k=0; k<5; ++k){
            ret[k] = init_value[k];
        }

        return ret;
    }

    void ChartOptimizor::EnergyFuncGrad(const alglib::real_1d_array& x, double& func_value, alglib::real_1d_array& grad_vec, void* ptr)
    {
        boost::shared_ptr<MeshModel> p_mesh = m_parameter.GetMeshModel();
        const vector<int>& face_index_vec = m_patch.m_face_index_array;
        const PolyIndexArray& mesh_face = p_mesh->m_Kernel.GetFaceInfo().GetIndex();
        const DoubleArray& face_area_vec = p_mesh->m_Kernel.GetFaceInfo().GetFaceArea();

        func_value = 0;
        vector< matrix<double>  > right_mat(face_index_vec.size());
        matrix<double> E_J(5, 1);
        
        for(size_t k=0; k<face_index_vec.size(); ++k){
            right_mat[k] = CalRightMatrix(face_index_vec[k]);
        }
        
        for(size_t k=0; k<face_index_vec.size(); ++k){
            int fid = face_index_vec[k];
            double area = face_area_vec[fid];
            const matrix<double>& r_mat = right_mat[k];

            matrix<double> cur_e_j(5, 1);
            calc_e_j(&cur_e_j[0], &x[0], &r_mat[0]);                
            E_J += cur_e_j * area;
                
            matrix<double> cur_E(1,1);
            calc_e(&cur_E[0], &x[0], &r_mat[0]);
            func_value += cur_E(0,0)*area;
        }

        for(size_t k=0; k<5; ++k) grad_vec[k] = E_J(k, 0);
    }

    matrix<double> ChartOptimizor::CalRightMatrix(int fid) const
    {
        boost::shared_ptr<MeshModel> p_mesh = m_parameter.GetMeshModel();
		const CoordArray& vert_coord_array = p_mesh->m_Kernel.GetVertexInfo().GetCoord();
		const IndexArray& faces = (p_mesh->m_Kernel.GetFaceInfo().GetIndex())[fid];
        const vector<ParamCoord>& vert_pc_vec = m_parameter.GetVertexParamCoordArray(); 
		std::vector<Coord> tri_vert_coord_3d(3);
		for(size_t k=0; k<3; ++k) tri_vert_coord_3d[k] = vert_coord_array[faces[k]];
		vector<Coord2D> local_coord = ComputeTriangleLocal2DCoord(tri_vert_coord_3d);
        matrix<double> lp_mat(3, 3); //! local coordinate for this triangle
        for(int i=0; i<3; ++i){
            lp_mat(0, i) = local_coord[i][0];
            lp_mat(1, i) = local_coord[i][1];
            lp_mat(2, i) = 1;
        }
        inv(lp_mat);
        
        matrix<double> r_mat(3, 2);
        r_mat(0, 0) = 1; r_mat(1, 0) = 0; r_mat(2, 0) = 0;
        r_mat(0, 1) = 0; r_mat(1, 1) = 1; r_mat(2, 1) = 0;

        matrix<double> l_mat(3, 3);
        for(size_t k=0; k<3; ++k){
            int vid = faces[k];
            ParamCoord pc = vert_pc_vec[vid];
            int vert_cid = m_parameter.GetVertexChartID(vid);
            if(vert_cid != m_patch.patch_id){
                m_parameter.TransParamCoordBetweenCharts( vert_cid, m_patch.patch_id, vid, vert_pc_vec[vid], pc);
            }
            l_mat(0, k) = pc.s_coord * (1.0 - pc.t_coord);
            l_mat(1, k) = pc.s_coord * pc.t_coord;
            l_mat(2, k) = pc.t_coord * (1.0 - pc.s_coord);
        }

        matrix<double> tmp_mat = lp_mat*r_mat;
        return l_mat*tmp_mat;
    }
}
