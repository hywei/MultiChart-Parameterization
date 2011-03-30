#ifndef LOCALPARAMETER_H_
#define LOCALPARAMETER_H_

#include "Parameterization.h"

#include <vector>
#include <boost/shared_ptr.hpp>

namespace PARAM{

    class LocalParameter
    {
    pubilc:
        LocalParameter(boost::shared_ptr<MeshModel> _mesh);
        ~LocalParameter();

        void SetBoundaryVertexArray(const std::vector<size_t>& bound_vert_array){ m_bound_vert_vec = bound_vert_array;}
        void SetInnerVertexArray(const std::vector<size_t>& inner_vert_array){ m_inner_vert_vec = inner_vert_array;}
        void SetBoundaryVertexPCArray(const std::vector<ParamCoord>& bd_pc_array){ m_bound_pc_vec = bd_pc_array; }
        const std::vector<ParamCoord>& GetInnerVertPCArray() const { return m_inner_pc_vec; }

        void ComputeConformalMap();

    private:
        void ArcLengthParameter(double agl1, double agl2, const std::vector<size_t>& _path, std::vector<ParamCoord>& pc_vec);

    private:
        boost::shared_ptr<MeshModel> p_mesh;
        std::vector<size_t> m_bound_vert_vec; //! boundary vertices set
        std::vector<size_t> m_inner_vert_vec; //! inner vertices set

        std::vector<ParamCoord> m_bound_pc_vec; //! boundary vertices parameter coordinate
        std::vector<ParamCoord> m_inner_pc_vec; //! inner vertices parameter coordinate
        std::vector<size_t> m_bound_corner_vec;
        
    };
}

#endif
