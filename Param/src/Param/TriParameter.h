#ifndef TRIPARAMETER_H_
#define TRIPARAMETER_H_

#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>
#include "Parameterization.h"

class MeshModel;
class LinearSolver;
class CMeshSparseMatrix;

namespace PARAM
{
    class TriChartCreator;

    class TriParameter
    {
    public:
        TriParameter(boost::shared_ptr<MeshModel> _p_mesh);
        ~TriParameter();

        bool ComputeParamCoord();

    public:
        //! IO Function
        bool LoadTriPatchFile(const std::string& tri_patch_file);

        int GetVertexChartID(int vid) const { return m_vert_chart_array[vid]; }
        ParamCoord GetVertexParamCoord(int vid) const { return m_vert_param_coord_array[vid]; }

        const std::vector<int>& GetVertexChartArray() const { return m_vert_chart_array; }

        boost::shared_ptr<MeshModel> GetMeshModel() const { return p_mesh; }

    private:
        boost::shared_ptr<MeshModel> p_mesh;
        boost::shared_ptr<TriChartCreator> p_quad_chart_creator;

        std::vector<int> m_vert_chart_array;
        std::vector<int> m_face_chart_array;

        std::vector<ParamCoord> m_vert_param_coord_array;

        std::vector < std::vector<int> > m_chart_vertices_array;

        
    };
        
}

#endif //TRIPARAMETER_H_
