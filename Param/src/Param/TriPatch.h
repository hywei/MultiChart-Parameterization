#ifndef TRIPATCH_H_
#define TRIPATCH_H_

#include <vector>

namespace PARAM
{
    class TriPatch
    {
    public:
        TriPatch(){}
        ~TriPatch(){}

    public:
        std::vector<int> m_conner_index_array;
        std::vector<int> m_edge_index_array;
        std::vector<int> m_face_index_array;
        std::vector<int> m_neighbor_patch_array;
    };
  
}


#endif //TRIPATCH_H_
