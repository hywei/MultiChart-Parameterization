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
        CheckJacobian();        
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
        cout << "\t" << int(report.terminationtype) << endl;
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
            int fid = face_index_vec[k];
            right_mat[k] = CalTempMatrix_2(fid) * CalTempMatrix_1(fid);
        }
        
        for(size_t k=0; k<face_index_vec.size(); ++k){
            int fid = face_index_vec[k];
            double area = face_area_vec[fid];
            const matrix<double>& r_mat = right_mat[k];

            zjucad::matrix::matrix<double> cur_e_j(5, 1);
            calc_e_j(&cur_e_j[0], &x[0], &r_mat[0]);                
            E_J += cur_e_j * area;
                
            matrix<double> cur_E(1,1);
            calc_e(&cur_E[0], &x[0], &r_mat[0]);
            func_value += cur_E(0,0)*area;
        }

        for(size_t k=0; k<5; ++k){
            grad_vec[k] = E_J(k, 0);
        }
    }

        /****************** Solve With Newton Method ****************************/
    vector<double> ChartOptimizor::SolveWithNewtonMethod() const
    {
        boost::shared_ptr<MeshModel> p_mesh = m_parameter.GetMeshModel();
        const vector<int>& face_index_vec = m_patch.m_face_index_array;

        const PolyIndexArray& mesh_face = p_mesh->m_Kernel.GetFaceInfo().GetIndex();
        const DoubleArray& face_area_vec = p_mesh->m_Kernel.GetFaceInfo().GetFaceArea();

        zjucad::matrix::matrix<double> x(5, 1);
        zjucad::matrix::matrix<double> delta_x(5,1);

        x = m_init_value;
        if(x.size() != 5) x = CalInitialValue();

        int face_num = face_index_vec.size();
        vector< matrix<double>  > right_mat(face_num);
        
        for(size_t k=0; k<face_index_vec.size(); ++k){
            int fid = face_index_vec[k];
            matrix<double> tmp_mat_1 = CalTempMatrix_1(fid);
            matrix<double> tmp_mat_2 = CalTempMatrix_2(fid);
            right_mat[k] = tmp_mat_2 * tmp_mat_1;
        }
        
        int iter_times = 20;
        double prev_energy = -1.0;
        for(int it=0; it<iter_times; ++it){
            double E_total = 0;
            zjucad::matrix::matrix<double> E_J(5, 1);
            zjucad::matrix::matrix<double> E_H(5, 5);
            
            for(size_t k=0; k<face_index_vec.size(); ++k){
                int fid = face_index_vec[k];
                double area = face_area_vec[fid];
                const matrix<double>& r_mat = right_mat[k];

                zjucad::matrix::matrix<double> cur_e_j(5, 1);
                zjucad::matrix::matrix<double> cur_e_h(5, 5);
                calc_e_j(&cur_e_j[0], &x[0], &r_mat[0]);
                calc_e_h(&cur_e_h[0], &x[0], &r_mat[0]);
                
                E_J += cur_e_j * area;
                E_H += cur_e_h * area;                
                
                zjucad::matrix::matrix<double> cur_E(1,1);
                calc_e(&cur_E[0], &x[0], &r_mat[0]);
                E_total += cur_E(0,0)*area;
            }

            matrix<double> inv_h(E_H);
            inv(inv_h);
            delta_x = inv_h * ( E_J * (-1.0));
            x += delta_x;

            cout << "Gradient : " << dot(delta_x, E_J) << std::endl;
            cout << "Current x: ";
            for(int i=0; i<5; ++i) std::cout << x(i,0) <<" ";
            std::cout << std::endl;
            std::cout << "Delta x: ";
            for(int i=0; i<5; ++i) std::cout << delta_x(i,0) << " ";
            std::cout << std::endl;
            std::cout << "Total Energy: " << E_total << std::endl;
            std::cout << std::endl;
            
            if(fabs(prev_energy - E_total) < 1e-8 || it == iter_times -1 ){  
                std::cout <<"Gradient Vector:";
                for(size_t j=0; j<5; ++j) cout << E_J(j,0) << " ";
                std::cout << endl;
                break;
            }
            prev_energy = E_total;
        }

        vector<double> ret(5);
        for(size_t k=0; k<5; ++k) ret[k] = x(k, 0);
        return ret;
    }

    
    double ChartOptimizor::CalTriIsometricDistoritionEnergy(int fid) const
    {
        boost::shared_ptr<MeshModel> p_mesh = m_parameter.GetMeshModel();
		const CoordArray& vert_coord_array = p_mesh->m_Kernel.GetVertexInfo().GetCoord();
		const IndexArray& faces = (p_mesh->m_Kernel.GetFaceInfo().GetIndex())[fid];
        const vector<ParamCoord>& pc_vec = m_parameter.GetVertexParamCoordArray();
        zjucad::matrix::matrix<double> uv_mat(2, 3);
        for(size_t i=0; i<3; ++i){
            int vid = faces[i];
            uv_mat(0, i) = pc_vec[vid].s_coord;
            uv_mat(1, i) = pc_vec[vid].t_coord;
        }
        
		vector<Coord> tri_vert_coord_3d(3);
		for(size_t k=0; k<3; ++k) tri_vert_coord_3d[k] = vert_coord_array[faces[k]];
		vector<Coord2D> local_coord = ComputeTriangleLocal2DCoord(tri_vert_coord_3d);
        zjucad::matrix::matrix<double> tmp_1(3, 3);
        for(int i=0; i<3; ++i){
            tmp_1(0, i) = local_coord[i][0];
            tmp_1(1, i) = local_coord[i][1];
            tmp_1(2, i) = 1;
        }
        inv(tmp_1);
        zjucad::matrix::matrix<double> tmp_2(3, 2);
        tmp_2(0,1) = tmp_2(1, 0) = tmp_2(2, 0) = tmp_2(2, 1) = 0;
        tmp_2(0,0) = tmp_2(1,1) = 1;

        zjucad::matrix::matrix<double> tmp_mat = uv_mat*tmp_1;
        zjucad::matrix::matrix<double> j_mat = tmp_mat*tmp_2;

        zjucad::matrix::matrix<double> e_mat = trans(j_mat)*j_mat - (zjucad::matrix::eye<double>(2));

        return e_mat(0,0)*e_mat(0,0) + e_mat(0,1)*e_mat(0,1) + e_mat(1,0)*e_mat(1,0) + e_mat(1,1)*e_mat(1,1);
    }

    
    matrix<double> ChartOptimizor::CalTempMatrix_1(int fid) const
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

        return lp_mat*r_mat;
    }

    matrix<double> ChartOptimizor::CalTempMatrix_2(int fid) const
    {
        const vector<ParamCoord>& vert_pc_vec = m_parameter.GetVertexParamCoordArray();
        boost::shared_ptr<MeshModel> p_mesh = m_parameter.GetMeshModel();
		const IndexArray& faces = (p_mesh->m_Kernel.GetFaceInfo().GetIndex())[fid];
        matrix<double> tmp_mat(3, 3);
        for(size_t k=0; k<3; ++k){
            int vid = faces[k];
            const ParamCoord& pc = vert_pc_vec[vid];
            tmp_mat(0, k) = pc.s_coord * (1.0 - pc.t_coord);
            tmp_mat(1, k) = pc.s_coord * pc.t_coord;
            tmp_mat(2, k) = pc.t_coord * (1.0 - pc.s_coord);
        }
        return tmp_mat;
    }

    void ChartOptimizor::CheckJacobian() const
    {
        boost::shared_ptr<MeshModel> p_mesh = m_parameter.GetMeshModel();
        const vector<int>& face_index_vec = m_patch.m_face_index_array;
        const PolyIndexArray& mesh_face = p_mesh->m_Kernel.GetFaceInfo().GetIndex();
        const DoubleArray& face_area_vec = p_mesh->m_Kernel.GetFaceInfo().GetFaceArea();
        const vector<ParamCoord>& pc_vec = m_parameter.GetVertexParamCoordArray();
        for(size_t k=0; k<pc_vec.size(); ++k){
            cout << pc_vec[k].s_coord << " " << pc_vec[k].t_coord << endl;
        }
        
        std::vector<double> face_area_in_param_domain(mesh_face.size(), 0);
        for(size_t k=0; k<face_index_vec.size(); ++k){
            int fid = face_index_vec[k];
            
            double iso_distortion = CalTriIsometricDistoritionEnergy(fid);
            cout << "Distortion: " << iso_distortion << endl;
            cout << "Face Area: " << face_area_vec[fid] << endl;

            zjucad::matrix::matrix<double> j_mat_1 = CalJacMat(fid);
            zjucad::matrix::matrix<double> j_mat_2 = CalJacMat_B(fid);

            std::cout << "Jacobi 1: " << std::endl;
            std::cout << j_mat_1 << std::endl;
            std::cout << "Jacobi 2: " << std::endl;
            std::cout << j_mat_2 << std::endl;
        }
    }
    
    zjucad::matrix::matrix<double> ChartOptimizor::CalJacMat(int fid) const
    {
        boost::shared_ptr<MeshModel> p_mesh = m_parameter.GetMeshModel();
		const CoordArray& vert_coord_array = p_mesh->m_Kernel.GetVertexInfo().GetCoord();
		const IndexArray& faces = (p_mesh->m_Kernel.GetFaceInfo().GetIndex())[fid];


		vector<Coord> tri_vert_coord_3d(3);
		for(size_t k=0; k<3; ++k) tri_vert_coord_3d[k] = vert_coord_array[faces[k]];
		vector<Coord2D> local_coord = ComputeTriangleLocal2DCoord(tri_vert_coord_3d);
		double area_2 = 2*(p_mesh->m_Kernel.GetFaceInfo().GetFaceArea())[fid];
        
		/// get these three vertices's parameter coordinate
		int cid = m_parameter.GetFaceChartID(fid);
		vector<ParamCoord> vert_pc(3);
		for(size_t k=0; k<3; ++k){
			int cur_cid = m_parameter.GetVertexChartID(faces[k]);
			ParamCoord cur_pc = m_parameter.GetVertexParamCoord(faces[k]);
			vert_pc[k] = cur_pc;
			if(cur_cid != cid){
				m_parameter.TransParamCoordBetweenCharts(cur_cid, cid, faces[k],cur_pc, vert_pc[k]);
			}
		}

        // Algorithm : Sig2007 parameterization course, p40, equation(4.8)
        matrix<double> left_mat(2,2); //! matrix: ([0, -1], [1, 0])
        left_mat(0, 0) = 0; left_mat(0, 1) = -1;
        left_mat(1, 0) = 1; left_mat(1, 1) = 0;

        matrix<double> lc_mat(2, 3); //! local coordinate of this triangle
        for(size_t k=0; k<3; ++k){
            lc_mat(0, k) = local_coord[k][0]; // X
            lc_mat(1, k) = local_coord[k][1]; // Y
        }        

        matrix<double> mid_mat(3, 3); //! ([0, 1, -1], [-1, 0, 1], [1, -1, 0])
        mid_mat(0, 0) = 0; mid_mat(0, 1) = 1; mid_mat(0, 2) = -1;
        mid_mat(1, 0) = -1; mid_mat(1, 1) = 0; mid_mat(1, 2) = 1;
        mid_mat(2, 0) = 1; mid_mat(2, 1) = -1; mid_mat(2, 2) = 0;

        matrix<double> uv_mat(3, 2); //! parameter coordinate matrix
        for(size_t k=0; k<3; ++k){
            uv_mat(k, 0) = vert_pc[k].s_coord;
            uv_mat(k, 1) = vert_pc[k].t_coord;
        }
        
        matrix<double> tmp_l = left_mat * lc_mat;
        matrix<double> tmp_r = mid_mat * uv_mat;
        matrix<double> res = tmp_l * tmp_r;

        return res * (1.0/area_2);
    }
    matrix<double> ChartOptimizor::CalJacMat_B(int fid) const
    {
        matrix<double> r_mat = CalTempMatrix_1(fid);
        matrix<double> uv_mat(2, 3);

        boost::shared_ptr<MeshModel> p_mesh = m_parameter.GetMeshModel();
		const IndexArray& faces = (p_mesh->m_Kernel.GetFaceInfo().GetIndex())[fid];
		/// get these three vertices's parameter coordinate
		int chart_id = m_parameter.GetFaceChartID(fid);
		vector<ParamCoord> vert_param_corod(3);
		for(size_t k=0; k<3; ++k){
			int vid = faces[k];
			int cur_chart_id = m_parameter.GetVertexChartID(vid);
			ParamCoord cur_param_coord = m_parameter.GetVertexParamCoord(vid);
			vert_param_corod[k] = cur_param_coord;
			if(cur_chart_id != chart_id){
				m_parameter.TransParamCoordBetweenCharts(cur_chart_id, chart_id, vid,cur_param_coord, vert_param_corod[k]);
			}

            uv_mat(0, k) = cur_param_coord.s_coord;
            uv_mat(1, k) = cur_param_coord.t_coord;
		}

        return uv_mat * r_mat;
    }
}
