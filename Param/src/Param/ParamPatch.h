#ifndef PARAMPATCH_H_
#define PARAMPATCH_H_

#include <vector>

namespace PARAM
{
    class PatchConner
    {
    public:
        PatchConner(){}
        ~PatchConner() {}

    public:
        int m_conner_type; //! 0: minimum, 1: saddle, 2:maximum
        int m_mesh_index;
        std::vector<int> m_nb_conner_index_array; //! neighbor conners index
        std::vector<int> m_nb_edge_index_array; //! neighbor edges index
    };

    class PatchEdge
    {
    public:
        PatchEdge(){}
        ~PatchEdge(){}

    public:
        std::pair< int, int> m_conner_pair_index; //! two endpoint conners' index
        std::vector<int> m_mesh_path; //! the mesh vertices index
        std::vector<int> m_nb_patch_index_array; //! the neighbor patch index
        
    };

    class ParamPatch
    {
    public:
        ParamPatch() { m_valid_range = std::make_pair(0.0, 1.0); patch_id = -1;}
        ~ParamPatch(){}

        
        //! check the parameter coordinate is valid or not
        bool InValidRangle(const ParamCoord& param_coord) const
        {
            return InsideValidChartRegion(param_coord);
        }

        //! check this parameter coordinate inside polygen or not
        bool InsideValidChartRegion(const ParamCoord& param_coord) const
        {
            Coord2D p(param_coord.s_coord, param_coord.t_coord);
            std::vector<Coord2D> polygen;
            for(size_t k=0; k<m_conner_pc_array.size(); ++k){
                const ParamCoord& pc = m_conner_pc_array[k];
                polygen.push_back(Coord2D(pc.s_coord, pc.t_coord));
            }

            return InsidePolygen(p, polygen);
        }
        
    public:
        std::vector<int> m_conner_index_array; //! index in patch conner array
        std::vector<int> m_edge_index_array; //! index in patch edge array
        std::vector<int> m_face_index_array; //! index in mesh face array
        std::vector<int> m_nb_patch_index_array; //! neighbor patch index

        //! the corresponding chart
        std::vector<ParamCoord> m_conner_pc_array; //! each corner's parameter coordinate
        std::pair<double, double> m_valid_range;

        int patch_id;
    };
}

#endif //PARAMPATCH_H_
