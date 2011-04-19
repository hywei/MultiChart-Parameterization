#include "ChartCreator.h"
#include "../ModelMesh/MeshModel.h"
#include <fstream>
#include <cmath>
#include <iostream>
#include <queue>
#include <set>

namespace PARAM
{
    ChartCreator::ChartCreator(boost::shared_ptr<MeshModel> _p_mesh):
        p_mesh(_p_mesh) {}

    ChartCreator::~ChartCreator(){}

    bool ChartCreator::LoadPatchFile(const std::string& patch_file)
    {
        std::ifstream fin(patch_file.c_str());
        if(fin.fail()){
            std::cerr << "@@@Eroor : Load Patch File Fail!" << std::endl;
        }

        //! clear data
        m_patch_conner_array.clear();
        m_patch_edge_array.clear();
        m_patch_array.clear();
        
        int patch_conner_num(0), patch_edge_num(0), patch_num(0);

        //! read patch conner info
        fin >> patch_conner_num;
        m_patch_conner_array.resize(patch_conner_num);

        for(int k=0; k<patch_conner_num; ++k){
            PatchConner& patch_conner = m_patch_conner_array[k];
            fin >> patch_conner.m_conner_type >> patch_conner.m_mesh_index;
            //! read neighbor info
            int neighbor_num;
            fin >> neighbor_num;
            patch_conner.m_nb_conner_index_array.resize(neighbor_num);
            patch_conner.m_nb_edge_index_array.resize(neighbor_num);
            for(int i=0; i<neighbor_num; ++i){
                fin >> patch_conner.m_nb_conner_index_array[i];
                fin >> patch_conner.m_nb_edge_index_array[i];
            }
            
        }    

        //! read patch edge info
        fin >> patch_edge_num;
        m_patch_edge_array.resize(patch_edge_num);

        for(int k=0; k<patch_edge_num; ++k){
            PatchEdge& patch_edge = m_patch_edge_array[k];
            fin >> patch_edge.m_conner_pair_index.first >> patch_edge.m_conner_pair_index.second;
            
            int neighbor_num;
            fin >> neighbor_num;
            patch_edge.m_nb_patch_index_array.resize(neighbor_num);
            for(int i=0; i<neighbor_num; ++i){
                fin >> patch_edge.m_nb_patch_index_array[i];
            }

            int path_vert_num;
            fin >> path_vert_num;
            patch_edge.m_mesh_path.resize(path_vert_num);
            for(int i=0; i<path_vert_num; ++i){
                fin >> patch_edge.m_mesh_path[i];
            }
        }

        //! read patch info
        fin >> patch_num;
        m_patch_array.resize(patch_num);

        for(int k=0; k<patch_num; ++k){
            ParamPatch& patch = m_patch_array[k];
            patch.patch_id = (int) k;
            //! read face info
            int face_num;
            fin >> face_num;
            patch.m_face_index_array.resize(face_num);
            for(int i=0; i<face_num; ++i){
                fin >> patch.m_face_index_array[i];
            }

            //! read edge info
            int edge_num;
            fin >> edge_num;
            patch.m_edge_index_array.resize(edge_num);
            for(int i=0; i<edge_num; ++i){
                fin >> patch.m_edge_index_array[i];
            }                
        }
        fin.close();
        return true;
    }

    bool ChartCreator::FormParamCharts()
    {
		SetPatchConners();
		SetPatchNeighbors();
        
        for(size_t k=0; k<m_patch_array.size(); ++k){
            ParamPatch& patch = m_patch_array[k];
            std::vector<ParamCoord>& corner_pc_vec = patch.m_conner_pc_array;
            corner_pc_vec.clear();
            corner_pc_vec.resize(patch.m_conner_index_array.size());
            //!TODO: this is just set for quad-patch
            corner_pc_vec[0] = ParamCoord(0, 0);
            corner_pc_vec[1] = ParamCoord(1, 0);
            corner_pc_vec[2] = ParamCoord(1, 1);
            corner_pc_vec[3] = ParamCoord(0, 1);
        }

        std::cout << "Form " << m_patch_array.size() <<" Charts!" << std::endl;
        return true;
    }    

	void ChartCreator::SetQuadChartConnerParamCoord(int chart_id)
	{
		const ParamPatch& patch = m_patch_array[chart_id];
		std::vector<ParamCoord>& conner_pc_vec = m_patch_array[chart_id].m_conner_pc_array;
		conner_pc_vec.clear(); conner_pc_vec.resize(4);
		conner_pc_vec[0] = ParamCoord(0, 0);
		conner_pc_vec[1] = ParamCoord(1, 0);
		conner_pc_vec[2] = ParamCoord(1, 1);
		conner_pc_vec[3] = ParamCoord(0, 1);
	}

	void ChartCreator::SetTriangleChartConnerParamCoord(int chart_id)
	{
		const ParamPatch& patch = m_patch_array[chart_id];
		std::vector<ParamCoord>& conner_pc_vec = m_patch_array[chart_id].m_conner_pc_array;
		int vid1 = m_patch_conner_array[patch.m_conner_index_array[0]].m_mesh_index;
		int vid2 = m_patch_conner_array[patch.m_conner_index_array[1]].m_mesh_index;
		int vid3 = m_patch_conner_array[patch.m_conner_index_array[2]].m_mesh_index;

		const CoordArray& vCoord = p_mesh->m_Kernel.GetVertexInfo().GetCoord();
		double edge_len1 = (vCoord[vid1] - vCoord[vid2]).abs();	
		double edge_len2 = (vCoord[vid3] - vCoord[vid2]).abs();
		double edge_len3 = (vCoord[vid1] - vCoord[vid3]).abs();

		double agl = angle(vCoord[vid2] - vCoord[vid1], vCoord[vid3] - vCoord[vid1]);
//		edge_len1 = GetPatchEdgeLength(patch.m_edge_index_array[0]);
//		edge_len2 = GetPatchEdgeLength(patch.m_edge_index_array[1]);
//		edge_len3 = GetPatchEdgeLength(patch.m_edge_index_array[2]);

//		double val = (edge_len1*edge_len1 + edge_len3*edge_len3 - edge_len2*edge_len2)/(2.0*edge_len1*edge_len3);
//		agl = acos(val);

		conner_pc_vec.resize(3);
		conner_pc_vec[0] = ParamCoord(0, 0);
//		conner_pc_vec[1] = ParamCoord(1, 0); 
		conner_pc_vec[1] = ParamCoord(edge_len1, 0);
//		conner_pc_vec[2] = ParamCoord(0.5, sqrt(3.0)/2); 
		conner_pc_vec[2] = ParamCoord(edge_len3*cos(agl), edge_len3*sin(agl));

	}

	double ChartCreator::GetPatchEdgeLength(int edge_idx)
	{
		const PatchEdge& edge = m_patch_edge_array[edge_idx];
		const std::vector<int>& mesh_path = edge.m_mesh_path;

		const CoordArray& vCoord = p_mesh->m_Kernel.GetVertexInfo().GetCoord();
		double len=0.0;
		for(size_t i=1; i<mesh_path.size(); ++i)
		{
			len += (vCoord[mesh_path[i]] - vCoord[mesh_path[i-1]]).abs();
		}
		return len;
	}

    void ChartCreator::SetPatchConners()
    {
        //! we assume thate a patch's edges have been sorted by ccw
        //! then we can get conners with the same order

        for(size_t k=0; k<m_patch_array.size(); ++k){
			SetPatchConners(k);
        }            
    }

	void ChartCreator::SetPatchConners(int patch_id)
	{
		ParamPatch& patch = m_patch_array[patch_id];
		size_t edge_num = patch.m_edge_index_array.size();
		patch.m_conner_index_array.clear();
		patch.m_conner_index_array.resize(edge_num);

		std::vector<int> un_decide_conner;
		for(size_t i=0; i<patch.m_edge_index_array.size(); ++i){
			int cur_idx = patch.m_edge_index_array[i];
			int nxt_idx = patch.m_edge_index_array[(i+1)%edge_num];
			const PatchEdge& cur_edge = m_patch_edge_array[cur_idx];
			const PatchEdge& nxt_edge = m_patch_edge_array[nxt_idx];

			//! find the current edge's corresponding conner
			std::vector<int> conners_1(2), conners_2(2);
			conners_1[0] = cur_edge.m_conner_pair_index.first;
			conners_1[1] = cur_edge.m_conner_pair_index.second;
			conners_2[0] = nxt_edge.m_conner_pair_index.first;
			conners_2[1] = nxt_edge.m_conner_pair_index.second;

			if(conners_1[0]>conners_1[1]) std::swap(conners_1[0], conners_1[1]);
			if(conners_2[0]>conners_2[1]) std::swap(conners_2[0], conners_2[1]);
			if(conners_1 == conners_2){
				std::cout << "Warning : these two adjcent edges' endponts are all the same!" << std::endl;
				un_decide_conner.push_back(i);
			}else{
				if(conners_1[0] == conners_2[0] || conners_1[0] == conners_2[1]) patch.m_conner_index_array[i] = conners_1[1];
				else patch.m_conner_index_array[i] = conners_1[0];
			}
		}

		if(un_decide_conner.size() == patch.m_conner_index_array.size()){
			std::cerr << "Error: can't set conner for this patch " << patch_id << std::endl;
		}else{
			if(un_decide_conner.size() !=0 )
				std::cout << "Warning: there are undecide conner for patch " << patch_id << std::endl;
			//! TODO:
			for(size_t k=0; k<un_decide_conner.size(); ++k){

			}
		}
	}

    void ChartCreator::SetPatchNeighbors()
    {
        for(size_t k=0; k<m_patch_array.size(); ++k){
            ParamPatch& patch = m_patch_array[k];
            patch.m_nb_patch_index_array.clear();
            for(size_t i=0; i<patch.m_edge_index_array.size(); ++i){
                int edge_idx = patch.m_edge_index_array[i];
                const PatchEdge& patch_edge = m_patch_edge_array[edge_idx];
                if(patch_edge.m_nb_patch_index_array.size() == 2){
                    int nb_patch_idx = (patch_edge.m_nb_patch_index_array[0] == k) ? 
						patch_edge.m_nb_patch_index_array[1] : patch_edge.m_nb_patch_index_array[0];
                    patch.m_nb_patch_index_array.push_back(nb_patch_idx);
                }
                    
            }
			
        }
        
    }

    bool ChartCreator::IsAmbiguityChartPair(int chart_id_1, int chart_id_2) const
    {
		const ParamPatch& patch_1 = m_patch_array[chart_id_1];
		const ParamPatch& patch_2 = m_patch_array[chart_id_2];

		/// find the common edge of these two charts
		const std::vector<int>& patch_edge_array_1 = patch_1.m_edge_index_array;
		const std::vector<int>& patch_edge_array_2 = patch_2.m_edge_index_array;

		std::vector<int> common_edges;
		for(size_t k=0; k<patch_edge_array_1.size(); ++k){
			if(find(patch_edge_array_2.begin(), patch_edge_array_2.end(), patch_edge_array_1[k]) 
				!= patch_edge_array_2.end()) common_edges.push_back(patch_edge_array_1[k]);
		}

		if(common_edges.size() >=2 ) return true;
		return false;
        
    }
}
