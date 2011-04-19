#ifndef CHARTCREATOR_H_
#define CHARTCREATOR_H_

#include "Parameterization.h"
#include "ParamPatch.h"

#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>

#include <map>

class MeshModel;

namespace PARAM
{
  class ChartCreator
  {
    public:
      ChartCreator(boost::shared_ptr<MeshModel> _p_mesh);
      ~ChartCreator();

      bool LoadPatchFile(const std::string& patch_file);

      bool FormParamCharts();


    public:
      //! Get Methods
      boost::shared_ptr<MeshModel> GetMeshModel() const
      {
        return p_mesh;
      }
      const std::vector<ParamPatch>& GetPatchArray() const
      {
        return m_patch_array;
      }
      const std::vector<PatchConner>& GetPatchConnerArray() const
      {
        return m_patch_conner_array;
      }
      std::vector<PatchConner>& GetPatchConnerArray()
      {
        return m_patch_conner_array;
      }
      const std::vector<PatchEdge>& GetPatchEdgeArray() const
      {
        return m_patch_edge_array;
      }
      int GetPatchNumber() const
      {
        return (int) m_patch_array.size();
      }

      std::vector<ParamPatch>& GetPatchArray()
      {
        return m_patch_array;
      }

      //! check two charts is ambiguity(have more than one common edges)?
      bool IsAmbiguityChartPair(int chart_id_1, int chart_id_2) const;
    private:

      //! set each patch's conners
      void SetPatchConners();
      void SetPatchConners(int patch_id);

      //! set each patch's neighbor info
      void SetPatchNeighbors();
      double GetPatchEdgeLength(int edge_idx);

      void SetQuadChartConnerParamCoord(int);
      void SetTriangleChartConnerParamCoord(int);

    private:
      boost::shared_ptr<MeshModel> p_mesh;

      std::vector<ParamPatch> m_patch_array;
      std::vector<PatchConner> m_patch_conner_array;
      std::vector<PatchEdge> m_patch_edge_array;

      std::vector<int> m_unre_edge_index_array;

  };
}

#endif //CHARTCREATOR_H_
