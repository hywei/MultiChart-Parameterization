#ifndef TRICHARTCREATOR_H_
#define TRICHARTCREATOR_H_

#include "Parameterization.h"
#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>

class MeshModel;

namespace PARAM
{
    class TriChartCreator
    {
    public:
        TriChartCreator(boost::shared_ptr<MeshModel> _p_mesh);
        ~TriChartCreator();

        bool LoadPatchFile(const std::string& patch_file);

    public:
        
    private:
        boost::shared_ptr<MeshModel> p_mesh;
    }
}

#endif //TRICHARTCREATOR_H_
