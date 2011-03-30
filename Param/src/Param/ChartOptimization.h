#ifndef CHARTOPTIMIZATION_H_
#define CHARTOPTIMIZATION_H_

#include <boost/shared_ptr.hpp>

namespace PARAM{

    class Parameter;
    class ParamPatch;
    
    class ChartOptimizor
    {
    public:
        ChartOptimizor(const Parameter& parameter, const ParamPatch& patch);
        ~ChartOptimizor();

        void Optimization();

        const std::vector<double>& GetDofArray() const { return m_dof_vec; }
        

    private:
        void SetNewtonMethodInitValue(zjucad::matrix::matrix<double>& x) const;
        void SolveNewtonMethodIterValue(const zjucad::matrix::matrix<double>& j_mat, const zjucad::matrix::matrix<double>& h_mat, zjucad::matrix::matrix<double>& _x) const;
        void UpdateNewtonMethodValue(zjucad::matrix::matrix<double>& x, const zjucad::matrix::matrix<double>& delta_x) const;
        
        zjucad::matrix::matrix<double> CalTempMatrix_1(int fid) const;
        zjucad::matrix::matrix<double> CalTempMatrix_2(int fid) const;

        double CalFaceAreaInParamDomain(int fid) const;
        
    private:
        const Parameter& m_parameter;
        const ParamPatch& m_patch;
        
        std::vector<double> m_dof_vec;
    };
    
}

#endif
