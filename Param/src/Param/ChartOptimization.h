#ifndef CHARTOPTIMIZATION_H_
#define CHARTOPTIMIZATION_H_

#include <boost/shared_ptr.hpp>
#include <vector>

/// hjlib header
#include <hj_3rd/hjlib/math/blas_lapack.h>
#include <hj_3rd/zjucad/matrix/matrix.h>
#include <hj_3rd/zjucad/matrix/io.h>

/// alglib header
#include "../AlgLib/optimization.h"
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
        
        void SetInitialValue(const zjucad::matrix::matrix<double>& x){
            m_init_value = x;
        }
        
        void EnergyFuncGrad(const alglib::real_1d_array& x, double& func, alglib::real_1d_array& grad, void* ptr);
        
    private:
        zjucad::matrix::matrix<double> CalInitialValue() const;

        //! solve with simple newton method, return the optimization value
        std::vector<double> SolveWithNewtonMethod() const;
        //! solve with lbfgs by alglib
        std::vector<double> SolveWithLbfgs() const;
        
        zjucad::matrix::matrix<double> CalTempMatrix_1(int fid) const;
        zjucad::matrix::matrix<double> CalTempMatrix_2(int fid) const;
        double CalFaceAreaInParamDomain(int fid) const;

        double CalTriIsometricDistoritionEnergy(int fid) const;

        void CheckJacobian() const;
        zjucad::matrix::matrix<double> CalJacMat(int ) const;
        zjucad::matrix::matrix<double> CalJacMat_B(int ) const;
    private:
        const Parameter& m_parameter;
        const ParamPatch& m_patch;

        zjucad::matrix::matrix<double> m_init_value;
        std::vector<double> m_dof_vec;
    };
    
}

#endif
