#include "LocalParameter.h"
#include "../Numerical/MeshSparseMatrix.h"

using namespace std;

namespace PARAM{
    LocalParameter::LocalParameter(boost::shared_ptr<MeshModel> _mesh): p_mesh(_mesh){}
    LocalParameter::~LocalParameter(){}

    void LocalParameter::ArcLengthParameter(double agl1, double agl2, const std::vector<size_t>& _path, std::vector<ParamCoord>& pc_vec)
    {
        if(_path.size() <=2) return;
        pc_vec.clear();
        pc_vec.resize(_path.size());
        pc_vec[0] = ParamCoord(cos(agl1), sin(agl1));

        const CoordArray& vCoord = p_mesh->m_Kernel.GetVertInfo().GetCoord();

        double sum_len = 0;
        for(size_t i=1; i<_path.size(); ++i){
            sum_len += (vCoord[_path[i]] - vCoord[_path[i-1]]).abs();
        }

        double cur_len = 0, lambda;
        for(size_t i=1; i<_path.size(); ++i){
            cur_len += (vCoord[_path[i]] - vCoord[_path[i-1]]).abs();
            lambda = cur_len / sum_len;

            double cur_agl = (1-lambda)*agl1 + lambda*agl2;
            pc_vec[i] = ParamCoord(cos(cur_agl), sin(cur_agl));
        }
        
    }

    void LocalParameter::ComputeConformalMap()
    {
        CMeshSparseMatrix lap_mat;
        SetLapMatrixCoefWithMeanValueCoord(p_mesh, lap_mat);

        map<size_t, int> var_index_mapping;
        map<size_t, int> bdvert_index_mapping;
        for(size_t k=0; k<m_inner_vert_vec.size(); ++k){
            var_index_mapping[m_inner_vert_vec[k]] = (int)k;
        }
        for(size_t k=0; k<m_bound_vert_vec.size(); ++k){
            var_index_mapping[m_bound_vert_vec[k]] = -1;
            bdvert_index_mapping[m_bound_vert_vec[k]] = k;
        }
        
        int vari_num = m_inner_vert_vec.size()*2;

        LinearSolver linear_solver(vari_num);

        linear_solver.begin_equation();

        for(size_t k=0; k<m_inner_vert_vec.size(); ++k){
            int vid = m_inner_vert_vec[k];
            const vector<int>& row_index = lap_mat.m_RowIndex[vid];
            const vector<double>& row_data = lap_mat.m_RowData[vid];

            for(int st=0; st<2; ++st){
                linear_solver.begin_row();
                double right_b = 0;
                for(size_t i=0; i<row_index.size(); ++i){
                    int col_vert = row_index[i];
                    assert(var_index_mapping.find(col_vert) != var_index_mapping.end());
                    int var_index = var_index_mapping[col_vert];
                    double lap_weight = row_data[i];

                    if(var_index == -1){ //! boundary vertex
                        int bdv_index = bdvert_index_mapping[col_vert];
                        double st_value = (st==0) ? m_bound_pc_vec[bdv_index].s_coord : m_bound_pc_vec[bdv_index].t_coord;
                        right_b -= lap_*st_value;                
                    }else{
                        linear_solver.add_coefficient(var_index*2 + st, lap_weight);
                    }
                }
                linear_solver.set_right_hand_side(right_b);
                linear_solver.end_row();
            }
        }

        linear_solver.end_equation();

        linear_solver.solve();

        for(size_t k=0; k<m_inner_vert_vec.size(); ++k){
            int vid = m_inner_vert_vec[k];
            int var_index = var_index_mapping[vid];
            m_inner_pc_vec[var_index].s_coord = linear_solver.variable(var_index*2).value();
            m_inner_pc_vec[var_index].t_coord = linear_solver.variable(var_index*2+1).value();
        }
    }
    
}
