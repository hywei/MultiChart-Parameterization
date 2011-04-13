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

      void OptimizeAmbiguityPatchShape();

      //! optimize the patch shape, won't change the topology
      void OptimizePatchShape();


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

      void SetMeshEdgePatchEdgeMapping(std::map< std::pair<int, int>, std::vector<int> >& me_pe_mapping) const;
      bool FloodFillFaceForAllPatchs();
      bool FloodFillFaceAPatch(int init_fid, std::vector<bool>& face_visited_flag,
                               const std::map< std::pair<int, int>, std::vector<int> >& me_pe_mapping);


      //! patch shape optimization

      void SpliteDegradedPatch();
      //! this will splite a degraded patch to two ambiguity patchs
      void SpliteDegradedPatch(int patch_id);

      void ReConnectDegradePatch(int patch_id);
      //! check a patch is degraded or not
      bool IsDegradedPatch(const ParamPatch& patch) const;

      //! method for optimizing patch shape
      void OptimizeAmbiguityPatchPairShape(int patch_id_1, int patch_id_2);

      void FindPatchInnerFace(int patch_id);
      void FormPatchBoundary(int patch_id, std::vector<int>& boundary) const;

      void ValenceControl();
      int FindLongestPatchEdge(int conner_index) const;

      void RemovePatchEdgeToReduceValance(int rm_edge_idx, int high_valence_conner_idx);

      void FindShortestPathInPatch(int patch_id, int start_vid, int end_vid, std::vector<int>& path) const;
      void FindShortestPathInTwoAdjPatch(int patch_id1, int patch_id2,
                                         int start_vid, int end_vid, std::vector<int>& pach) const;


      void FormTriPatch();
      void SpliteQuadPatchToTriPatch(int patch_id);
      void SpliteDegradedPatchToTriPatch(int patch_id);
      void EdgeFlipForTriPacth(int patch_id1, int patch_id2);
      void ValenceControlForTriPatch();


      std::vector<int> GetCommonEdgesBetweenTwoPatch(const ParamPatch& patch_1, const ParamPatch& patch_2);
      std::vector<int> GetCommonConnersBetweenTwoPatch(const ParamPatch& patch_1, const ParamPatch& patch_2);
      int GetCommonConnerBetweenTwoEdge(const PatchEdge& edge_1, const PatchEdge& edge_2);

      void SetTriangleChartConnerParamCoord(int chart_id);
      void SetQuadChartConnerParamCoord(int chart_id);

      double GetPatchEdgeLength(int edge_idx);

      void GetTransListBetweenTwoChart(int chart_id1, int chart_id2, int edge_idx1, int edge_idx2);

    private:
      void AlginQuadPatchOrient();

    private:
      boost::shared_ptr<MeshModel> p_mesh;

      std::vector<ParamPatch> m_patch_array;
      std::vector<PatchConner> m_patch_conner_array;
      std::vector<PatchEdge> m_patch_edge_array;

      HalfEdge m_half_edge;
      std::vector<int> m_unre_edge_index_array;

  };
}

#endif //CHARTCREATOR_H_
