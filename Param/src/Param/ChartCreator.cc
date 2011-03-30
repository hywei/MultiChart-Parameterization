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
		m_half_edge.CreateHalfEdge(p_mesh);	   
		SetPatchConners();
		SetPatchNeighbors();
		
		
	//	FormTriPatch();

		//ValenceControl();

        //! this function should be call after call LoadPatchFile(const std::string&)
        
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
    	
	void ChartCreator::SpliteDegradedPatch()
	{
		size_t origin_patch_num = m_patch_array.size();
		for(size_t k=0; k<origin_patch_num; ++k){
			const ParamPatch& patch = m_patch_array[k];
			if(!IsDegradedPatch(patch)){
				SpliteQuadPatchToTriPatch(k);
				//SpliteDegradedPatch(k);
				//ReConnectDegradePatch(k);
			}
		}
	}

	bool ChartCreator::IsDegradedPatch(const ParamPatch& patch) const
	{
		//! a patch is degraded if its two edges are same
		for(size_t k=0; k<patch.m_edge_index_array.size(); ++k){
			for(size_t i=k+1; i<patch.m_edge_index_array.size(); ++i){
				if(patch.m_edge_index_array[i] == patch.m_edge_index_array[k]) return true;
			}
		}
		return false;
	}

	void ChartCreator::FormTriPatch()
	{
		SetPatchConners();
		SetPatchNeighbors();

		size_t origin_patch_num = m_patch_array.size();
		for(size_t k=0; k<origin_patch_num; ++k){
			const ParamPatch& patch = m_patch_array[k];
			if(!IsDegradedPatch(patch)){						
				SpliteQuadPatchToTriPatch(k);
			}
		}

		for(size_t k=0; k<origin_patch_num; ++k){
			const ParamPatch& patch = m_patch_array[k];
			if(IsDegradedPatch(patch)){
				SpliteDegradedPatchToTriPatch(k);
			}
		}

		SetPatchConners();
		SetPatchNeighbors();
	}

	void ChartCreator::SpliteQuadPatchToTriPatch(int patch_id)
	{
		std::cout << "Splite Quad Patch : " << patch_id << std::endl;
		ParamPatch& cur_patch = m_patch_array[patch_id];
		std::vector<int>& cur_conners = cur_patch.m_conner_index_array;
		std::vector<int>& cur_edges = cur_patch.m_edge_index_array;
		int val[4];
		for(size_t k=0; k<cur_conners.size(); ++k)
		{
			const PatchConner& conner = m_patch_conner_array[cur_conners[k]];
			val[k] = conner.m_nb_conner_index_array.size();
		}
		int val1 = std::max(val[0], val[2]), val2 = std::max(val[1], val[3]);
		if(val1 > val2){
			int tmp = cur_edges[0];
			for(size_t k=0; k<3; ++k) cur_edges[k] = cur_edges[k+1];
			cur_edges[3] = tmp;
			tmp = cur_conners[0];
			for(size_t k=0; k<3; ++k) cur_conners[k] = cur_conners[k+1];
			cur_conners[3] = tmp;
		}
		
		
		int add_edge_idx = (int)m_patch_edge_array.size();
		int add_patch_idx = (int)m_patch_array.size();
		PatchEdge add_edge;
		add_edge.m_conner_pair_index = std::make_pair(cur_conners[0], cur_conners[2]);
		add_edge.m_nb_patch_index_array.push_back(patch_id); add_edge.m_nb_patch_index_array.push_back(add_patch_idx);
		FindShortestPathInPatch(patch_id, m_patch_conner_array[cur_conners[0]].m_mesh_index, 
			m_patch_conner_array[cur_conners[2]].m_mesh_index, add_edge.m_mesh_path);
		
		ParamPatch add_patch;
	    add_patch.m_conner_index_array.push_back(cur_conners[0]); add_patch.m_edge_index_array.push_back(add_edge_idx);
		add_patch.m_conner_index_array.push_back(cur_conners[2]); add_patch.m_edge_index_array.push_back(cur_edges[2]);
		add_patch.m_conner_index_array.push_back(cur_conners[3]); add_patch.m_edge_index_array.push_back(cur_edges[3]);
		
		PatchConner& md_conner_1 = m_patch_conner_array[cur_conners[0]];
		std::vector<int>& nb_conners1 = md_conner_1.m_nb_conner_index_array;
		std::vector<int>& nb_edges1 = md_conner_1.m_nb_edge_index_array;
		int pos = find(nb_edges1.begin(), nb_edges1.end(), cur_edges[0]) - nb_edges1.begin();
		nb_edges1.insert(nb_edges1.begin() + pos+1, add_edge_idx);
		nb_conners1.insert(nb_conners1.begin()+pos+1, cur_conners[2]);

		PatchConner& md_conner_2 = m_patch_conner_array[cur_conners[2]];
		std::vector<int>& nb_conners2 = md_conner_2.m_nb_conner_index_array;
		std::vector<int>& nb_edges2 = md_conner_2.m_nb_edge_index_array;
		pos = find(nb_edges2.begin(), nb_edges2.end(), cur_edges[2]) - nb_edges2.begin();
		nb_edges2.insert(nb_edges2.begin() + pos +1, add_edge_idx);
		nb_conners2.insert(nb_conners2.begin() + pos+1, cur_conners[0]);


		PatchEdge& md_edge_1 = m_patch_edge_array[cur_edges[2]];
		if(md_edge_1.m_nb_patch_index_array[0] == patch_id) md_edge_1.m_nb_patch_index_array[0] = add_patch_idx;
		else md_edge_1.m_nb_patch_index_array[1] = add_patch_idx;
		PatchEdge& md_edge_2 = m_patch_edge_array[cur_edges[3]];
		if(md_edge_2.m_nb_patch_index_array[0] == patch_id) md_edge_2.m_nb_patch_index_array[0] = add_patch_idx;
		else md_edge_2.m_nb_patch_index_array[1] = add_patch_idx;

		cur_conners.erase(cur_conners.begin() + 3); 
		cur_edges.erase(cur_edges.begin() + 2, cur_edges.end()); 
		cur_edges.push_back(add_edge_idx);
		
		m_patch_edge_array.push_back(add_edge);
		m_patch_array.push_back(add_patch);

		FindPatchInnerFace(patch_id);
		FindPatchInnerFace(add_patch_idx);

	}

	void ChartCreator::SpliteDegradedPatchToTriPatch(int patch_id)
	{
		std::cout << "Splite Degraded Quad Patch : " << patch_id << std::endl;
		//! it will create one ambiguity patch pair
		ParamPatch& patch = m_patch_array[patch_id];
		const std::vector<int>& patch_edge_vec = patch.m_edge_index_array;
		int merged_edge_idx(-1),prev_edge_idx(-1), next_edge_idx(-1);
		for(size_t k=0; k<patch_edge_vec.size(); ++k){
			for(size_t i=k+1; i<patch_edge_vec.size(); ++i){
				if(patch_edge_vec[k] == patch_edge_vec[i]){
					merged_edge_idx = patch_edge_vec[i]; 
					prev_edge_idx = patch_edge_vec[(k+patch_edge_vec.size()-1)%patch_edge_vec.size()];
					next_edge_idx = patch_edge_vec[(i+1)%patch_edge_vec.size()];
					break; 
				}
			}
			if(merged_edge_idx != -1) break;
		}


		int merged_conner_idx(-1), single_conner_idx(-1), third_conner_idx(-1);
		const PatchEdge& merged_patch_edge = m_patch_edge_array[merged_edge_idx];
		int temp_conner_idx = merged_patch_edge.m_conner_pair_index.first;
		const PatchConner& temp_conner = m_patch_conner_array[temp_conner_idx];
		if(temp_conner.m_nb_conner_index_array.size() == 1) {
			merged_conner_idx = merged_patch_edge.m_conner_pair_index.second;
			single_conner_idx = merged_patch_edge.m_conner_pair_index.first;
		}else{
			merged_conner_idx = merged_patch_edge.m_conner_pair_index.first;
			single_conner_idx = merged_patch_edge.m_conner_pair_index.second;
		}

		const PatchEdge& prev_patch_edge = m_patch_edge_array[prev_edge_idx];
		if(prev_patch_edge.m_conner_pair_index.first == merged_conner_idx){
			third_conner_idx = prev_patch_edge.m_conner_pair_index.second;
		}else{
			third_conner_idx = prev_patch_edge.m_conner_pair_index.first;
		}

		const std::vector<int>& inner_path = m_patch_edge_array[merged_edge_idx].m_mesh_path;

		const std::vector<int>& face_vec = patch.m_face_index_array;		
		const PolyIndexArray& face_list_array = p_mesh->m_Kernel.GetFaceInfo().GetIndex();

		std::set< std::pair<int, int> > region_mesh_edge_set;
		for(size_t k=0; k<patch_edge_vec.size(); ++k){
			if(patch_edge_vec[k] == merged_edge_idx) continue;
			const std::vector<int>& path = m_patch_edge_array[patch_edge_vec[k]].m_mesh_path;
			for(size_t i=1; i<path.size(); ++i){
				region_mesh_edge_set.insert(MakeEdge(path[i-1], path[i]));
			}
		}

		for(size_t k=0; k<face_vec.size(); ++k){
			int fid = face_vec[k];
			const IndexArray& face = face_list_array[fid];
			for(int i=0; i<3; ++i){
				if(find(inner_path.begin(), inner_path.end(), face[(i+1)%3]) != inner_path.end() &&
					find(inner_path.begin(), inner_path.end(), face[i]) != inner_path.end()) continue;
				region_mesh_edge_set.insert(MakeEdge(face[i], face[(i+1)%3]));
			}
		}
		std::vector<int> add_path;
		FindShortestPathInRegion(p_mesh, m_patch_conner_array[single_conner_idx].m_mesh_index, 
			m_patch_conner_array[third_conner_idx].m_mesh_index, region_mesh_edge_set, add_path);

		int add_edge_idx = (int)m_patch_edge_array.size();
		int add_patch_idx = (int)m_patch_array.size();

		PatchEdge add_edge;
		add_edge.m_conner_pair_index = std::make_pair(single_conner_idx, third_conner_idx);
		add_edge.m_nb_patch_index_array.push_back(patch_id);
		add_edge.m_nb_patch_index_array.push_back(add_patch_idx);
		add_edge.m_mesh_path = add_path;

		ParamPatch add_patch;
		add_patch.m_conner_index_array.push_back(merged_conner_idx); add_patch.m_edge_index_array.push_back(next_edge_idx);
		add_patch.m_conner_index_array.push_back(third_conner_idx); add_patch.m_edge_index_array.push_back(add_edge_idx);
		add_patch.m_conner_index_array.push_back(single_conner_idx); add_patch.m_edge_index_array.push_back(merged_edge_idx);

		std::vector<int>& cur_conners = m_patch_array[patch_id].m_conner_index_array;
		std::vector<int>& cur_edges = m_patch_array[patch_id].m_edge_index_array;
		cur_conners.clear(); cur_edges.clear();
		cur_conners.push_back(merged_conner_idx); cur_edges.push_back(merged_edge_idx);
		cur_conners.push_back(single_conner_idx); cur_edges.push_back(add_edge_idx);
		cur_conners.push_back(third_conner_idx); cur_edges.push_back(prev_edge_idx);


		PatchConner& md_conner_1 = m_patch_conner_array[third_conner_idx];
		int pos = find(md_conner_1.m_nb_edge_index_array.begin(), md_conner_1.m_nb_edge_index_array.end(), prev_edge_idx)
			- md_conner_1.m_nb_edge_index_array.begin();
		md_conner_1.m_nb_edge_index_array.insert(md_conner_1.m_nb_edge_index_array.begin() + pos + 1, add_edge_idx);
		md_conner_1.m_nb_conner_index_array.insert(md_conner_1.m_nb_conner_index_array.begin() + pos + 1, single_conner_idx);

		PatchConner& md_conner_2 = m_patch_conner_array[single_conner_idx];
		md_conner_2.m_nb_conner_index_array.push_back(third_conner_idx);
		md_conner_2.m_nb_edge_index_array.push_back(add_edge_idx);


		PatchEdge& md_edge_1 = m_patch_edge_array[merged_edge_idx];
		md_edge_1.m_nb_patch_index_array.clear(); 
		md_edge_1.m_nb_patch_index_array.push_back(patch_id);
		md_edge_1.m_nb_patch_index_array.push_back(add_patch_idx);
		
		PatchEdge& md_edge_2 = m_patch_edge_array[next_edge_idx];
		if(md_edge_2.m_nb_patch_index_array[0] == patch_id) md_edge_2.m_nb_patch_index_array[0] = add_patch_idx;
		else md_edge_2.m_nb_patch_index_array[1] = add_patch_idx;

		m_patch_edge_array.push_back(add_edge);
		m_patch_array.push_back(add_patch);

		FindPatchInnerFace(patch_id);
		FindPatchInnerFace(add_patch_idx);

		
		{
			int patch_id_1 = m_patch_edge_array[prev_edge_idx].m_nb_patch_index_array[0];
			int patch_id_2 = m_patch_edge_array[prev_edge_idx].m_nb_patch_index_array[1];
			EdgeFlipForTriPacth(patch_id_1, patch_id_2);
			patch_id_1 = m_patch_edge_array[next_edge_idx].m_nb_patch_index_array[0];
			patch_id_2 = m_patch_edge_array[next_edge_idx].m_nb_patch_index_array[1];
			EdgeFlipForTriPacth(patch_id_1, patch_id_2);
		}

	}

	void ChartCreator::EdgeFlipForTriPacth(int patch_id1, int patch_id2)
	{
	    ParamPatch& patch_1 = m_patch_array[patch_id1];
		ParamPatch& patch_2 = m_patch_array[patch_id2];

		std::vector<int>& edge_vec_1 = patch_1.m_edge_index_array;
		std::vector<int>& conner_vec_1 = patch_1.m_conner_index_array;
		std::vector<int>& edge_vec_2 = patch_2.m_edge_index_array;
		std::vector<int>& conner_vec_2 = patch_2.m_conner_index_array;

		int conner_idx_1(-1), conner_idx_2(-1), conner_idx_3(-1), conner_idx_4(-1);
		int edge_idx_1(-1), edge_idx_2(-1), edge_idx_3(-1), edge_idx_4(-1), com_edge_idx(-1);

		std::vector<int> com_edges = GetCommonEdgesBetweenTwoPatch(patch_1, patch_2);
		if(com_edges.size() != 1) return;
		com_edge_idx = com_edges[0];
		std::vector<int> com_conners = GetCommonConnersBetweenTwoPatch(patch_1, patch_2);
		if(com_conners.size() != 2) return;
		conner_idx_1 = com_conners[0]; conner_idx_3 = com_conners[1];
		
		int ps = find(edge_vec_1.begin(), edge_vec_1.end(), com_edge_idx) - edge_vec_1.begin();
		edge_idx_1 = edge_vec_1[(ps+1)%3]; edge_idx_2 = edge_vec_1[(ps+2)%3];
		ps = find(edge_vec_2.begin(), edge_vec_2.end(), com_edge_idx) - edge_vec_2.begin();
		edge_idx_3 = edge_vec_2[(ps+1)%3]; edge_idx_4 = edge_vec_2[(ps+2)%3];

		PatchEdge& edge_1 = m_patch_edge_array[edge_idx_1];
		PatchEdge& edge_2 = m_patch_edge_array[edge_idx_2];
		PatchEdge& edge_3 = m_patch_edge_array[edge_idx_3];
		PatchEdge& edge_4 = m_patch_edge_array[edge_idx_4];
		PatchEdge& com_edge = m_patch_edge_array[com_edge_idx];

		conner_idx_1 = GetCommonConnerBetweenTwoEdge(edge_1, com_edge);
		conner_idx_2 = GetCommonConnerBetweenTwoEdge(edge_1, edge_2);
		conner_idx_3 = GetCommonConnerBetweenTwoEdge(edge_3, com_edge);
		conner_idx_4 = GetCommonConnerBetweenTwoEdge(edge_3, edge_4);
		PatchConner& conner_1 = m_patch_conner_array[conner_idx_1];
		PatchConner& conner_2 = m_patch_conner_array[conner_idx_2];
		PatchConner& conner_3 = m_patch_conner_array[conner_idx_3];
		PatchConner& conner_4 = m_patch_conner_array[conner_idx_4];

		com_edge.m_conner_pair_index = std::make_pair(conner_idx_2, conner_idx_4);
		FindShortestPathInTwoAdjPatch(patch_id1, patch_id2, conner_2.m_mesh_index, conner_4.m_mesh_index, com_edge.m_mesh_path);
		
		ps = find(conner_1.m_nb_edge_index_array.begin(), conner_1.m_nb_edge_index_array.end(), com_edge_idx)
			- conner_1.m_nb_edge_index_array.begin();
		conner_1.m_nb_edge_index_array.erase(conner_1.m_nb_edge_index_array.begin() + ps);
		conner_1.m_nb_conner_index_array.erase(conner_1.m_nb_conner_index_array.begin() + ps);

		ps = find(conner_3.m_nb_edge_index_array.begin(), conner_3.m_nb_edge_index_array.end(), com_edge_idx)
			- conner_3.m_nb_edge_index_array.begin();
		conner_3.m_nb_edge_index_array.erase(conner_3.m_nb_edge_index_array.begin() + ps);
		conner_3.m_nb_conner_index_array.erase(conner_3.m_nb_conner_index_array.begin() + ps);

		ps = find(conner_2.m_nb_edge_index_array.begin(), conner_2.m_nb_edge_index_array.end(), edge_idx_2)
			- conner_2.m_nb_edge_index_array.begin();
		conner_2.m_nb_edge_index_array.insert(conner_2.m_nb_edge_index_array.begin() + ps+1, com_edge_idx);
		conner_2.m_nb_conner_index_array.insert(conner_2.m_nb_conner_index_array.begin() + ps+1, conner_idx_4);

		ps = find(conner_4.m_nb_edge_index_array.begin(), conner_4.m_nb_edge_index_array.end(), edge_idx_3)
			- conner_4.m_nb_edge_index_array.begin();
		conner_4.m_nb_edge_index_array.insert(conner_4.m_nb_edge_index_array.begin() + ps , com_edge_idx);
		conner_4.m_nb_conner_index_array.insert(conner_4.m_nb_conner_index_array.begin() + ps , conner_idx_2);

		if(edge_2.m_nb_patch_index_array[0] == patch_id1) edge_2.m_nb_patch_index_array[0] = patch_id2;
		else edge_2.m_nb_patch_index_array[1] = patch_id2;
		if(edge_4.m_nb_patch_index_array[0] == patch_id2) edge_4.m_nb_patch_index_array[0] = patch_id1;
		else edge_4.m_nb_patch_index_array[1] = patch_id1;

		edge_vec_1[0] = edge_idx_1; edge_vec_1[1] = com_edge_idx; edge_vec_1[2] = edge_idx_4;
		conner_vec_1[0] = conner_idx_1; conner_vec_1[1] = conner_idx_2; conner_vec_1[2] = conner_idx_4;
		edge_vec_2[0] = edge_idx_3; edge_vec_2[1] = com_edge_idx; edge_vec_2[2] = edge_idx_2;
		conner_vec_2[0] = conner_idx_3; conner_vec_2[1] = conner_idx_4; conner_vec_2[2] = conner_idx_2;

		FindPatchInnerFace(patch_id1);
		FindPatchInnerFace(patch_id2);
	}

	void ChartCreator::SpliteDegradedPatch(int patch_id)
	{
		//! it will create one ambiguity patch pair
		ParamPatch& patch = m_patch_array[patch_id];
		const std::vector<int>& patch_edge_vec = patch.m_edge_index_array;
		int merged_edge_idx(-1),prev_edge_idx(-1), next_edge_idx(-1);
		for(size_t k=0; k<patch_edge_vec.size(); ++k){
			for(size_t i=k+1; i<patch_edge_vec.size(); ++i){
				if(patch_edge_vec[k] == patch_edge_vec[i]){
					merged_edge_idx = patch_edge_vec[i]; 
					prev_edge_idx = patch_edge_vec[(k+patch_edge_vec.size()-1)%patch_edge_vec.size()];
					next_edge_idx = patch_edge_vec[(i+1)%patch_edge_vec.size()];
					break; 
				}
			}
			if(merged_edge_idx != -1) break;
		}


		int merged_conner_idx(-1), single_conner_idx(-1), third_conner_idx(-1);
		const PatchEdge& merged_patch_edge = m_patch_edge_array[merged_edge_idx];
		int temp_conner_idx = merged_patch_edge.m_conner_pair_index.first;
		const PatchConner& temp_conner = m_patch_conner_array[temp_conner_idx];
		if(temp_conner.m_nb_conner_index_array.size() == 1) {
			merged_conner_idx = merged_patch_edge.m_conner_pair_index.second;
			single_conner_idx = merged_patch_edge.m_conner_pair_index.first;
		}else{
			merged_conner_idx = merged_patch_edge.m_conner_pair_index.first;
			single_conner_idx = merged_patch_edge.m_conner_pair_index.second;
		}

		const PatchEdge& prev_patch_edge = m_patch_edge_array[prev_edge_idx];
		if(prev_patch_edge.m_conner_pair_index.first == merged_conner_idx){
			third_conner_idx = prev_patch_edge.m_conner_pair_index.second;
		}else{
			third_conner_idx = prev_patch_edge.m_conner_pair_index.first;
		}

		const std::vector<int>& inner_path = m_patch_edge_array[merged_edge_idx].m_mesh_path;

		const std::vector<int>& face_vec = patch.m_face_index_array;		
		const PolyIndexArray& face_list_array = p_mesh->m_Kernel.GetFaceInfo().GetIndex();

		std::set< std::pair<int, int> > region_mesh_edge_set;
		for(size_t k=0; k<patch_edge_vec.size(); ++k){
			if(patch_edge_vec[k] == merged_edge_idx) continue;
			const std::vector<int>& path = m_patch_edge_array[patch_edge_vec[k]].m_mesh_path;
			for(size_t i=1; i<path.size(); ++i){
				region_mesh_edge_set.insert(MakeEdge(path[i-1], path[i]));
			}
		}

		for(size_t k=0; k<face_vec.size(); ++k){
			int fid = face_vec[k];
			const IndexArray& face = face_list_array[fid];
			for(int i=0; i<3; ++i){
				if(find(inner_path.begin(), inner_path.end(), face[(i+1)%3]) != inner_path.end() &&
					find(inner_path.begin(), inner_path.end(), face[i]) != inner_path.end()) continue;
				region_mesh_edge_set.insert(MakeEdge(face[i], face[(i+1)%3]));
			}
		}
		std::vector<int> add_path;
	    FindShortestPathInRegion(p_mesh, m_patch_conner_array[single_conner_idx].m_mesh_index, 
			m_patch_conner_array[third_conner_idx].m_mesh_index, region_mesh_edge_set, add_path);

		int add_vert_num = add_path.size();
		assert(add_vert_num > 3);

		int mid_index = add_vert_num/2;
		int add_conner_vert = add_path[mid_index];
		std::vector<int> add_path_1 (add_path.begin(), add_path.begin() + mid_index+1);
		std::vector<int> add_path_2 (add_path.begin()+mid_index, add_path.end());

		int add_conner_index = (int)m_patch_conner_array.size();
		int add_edge_index_1 = (int) m_patch_edge_array.size();
		int add_edge_index_2 = (int) m_patch_edge_array.size() + 1;
		int add_patch_index = (int) m_patch_array.size();

		/// add a new patch conner		
		PatchConner add_conner;
		add_conner.m_conner_type = 1; /// assume it is a saddle, others will be ok
		add_conner.m_mesh_index = add_conner_vert;
		add_conner.m_nb_conner_index_array.push_back(single_conner_idx);  /// two neighbor
		add_conner.m_nb_conner_index_array.push_back(third_conner_idx);
		add_conner.m_nb_edge_index_array.push_back(add_edge_index_1);
		add_conner.m_nb_edge_index_array.push_back(add_edge_index_2);
		m_patch_conner_array.push_back(add_conner);

		/// add two new patch edges
		PatchEdge add_edge_1, add_edge_2;
		add_edge_1.m_conner_pair_index = std::make_pair(single_conner_idx, add_conner_index);
		add_edge_1.m_mesh_path = add_path_1;
		add_edge_1.m_nb_patch_index_array.push_back(patch_id);
		add_edge_1.m_nb_patch_index_array.push_back(add_patch_index);

		add_edge_2.m_conner_pair_index = std::make_pair(add_conner_index, third_conner_idx);
		add_edge_2.m_mesh_path = add_path_2;
		add_edge_2.m_nb_patch_index_array.push_back(patch_id);
		add_edge_2.m_nb_patch_index_array.push_back(add_patch_index);
		
		m_patch_edge_array.push_back(add_edge_1);
		m_patch_edge_array.push_back(add_edge_2);

		/// add a new patch
		ParamPatch add_patch;
		add_patch.m_edge_index_array.push_back(prev_edge_idx);
		add_patch.m_edge_index_array.push_back(merged_edge_idx);
		add_patch.m_edge_index_array.push_back(add_edge_index_1);
		add_patch.m_edge_index_array.push_back(add_edge_index_2);		
		add_patch.m_conner_index_array.push_back(third_conner_idx);
		add_patch.m_conner_index_array.push_back(merged_conner_idx);
		add_patch.m_conner_index_array.push_back(single_conner_idx);
		add_patch.m_conner_index_array.push_back(add_conner_index);
		m_patch_array.push_back(add_patch);

		/// modify some other patch info
		PatchConner& single_conner = m_patch_conner_array[single_conner_idx];

		single_conner.m_nb_conner_index_array.push_back(add_conner_index);
		single_conner.m_nb_edge_index_array.push_back(add_edge_index_1);

		PatchConner& third_conner = m_patch_conner_array[third_conner_idx];
		std::vector<int>& third_conner_nb_conners = third_conner.m_nb_conner_index_array;
		std::vector<int>& third_conner_nb_edges = third_conner.m_nb_edge_index_array;;
		std::vector<int>::iterator next_edge_iter = find(third_conner_nb_edges.begin(), 
			third_conner_nb_edges.end(), next_edge_idx);
		assert(next_edge_iter != third_conner_nb_edges.end());
		int pos = next_edge_iter - third_conner_nb_edges.begin();
		third_conner_nb_conners.insert(third_conner_nb_conners.begin() + pos, add_conner_index);
		third_conner_nb_edges.insert(next_edge_iter, add_edge_index_2);


		PatchEdge& merged_edge = m_patch_edge_array[merged_edge_idx];
		merged_edge.m_nb_patch_index_array.clear();
		merged_edge.m_nb_patch_index_array.push_back(patch_id);
		merged_edge.m_nb_patch_index_array.push_back(add_patch_index);

		PatchEdge& prev_edge = m_patch_edge_array[prev_edge_idx];
		for(size_t k=0; k<prev_edge.m_nb_patch_index_array.size(); ++k){
			if(prev_edge.m_nb_patch_index_array[k] == patch_id){
				prev_edge.m_nb_patch_index_array[k] = add_patch_index;
			}
		}

		ParamPatch& cur_patch = m_patch_array[patch_id];
		cur_patch.m_edge_index_array.clear();
		cur_patch.m_edge_index_array.push_back(next_edge_idx);
		cur_patch.m_edge_index_array.push_back(add_edge_index_2);
		cur_patch.m_edge_index_array.push_back(add_edge_index_1);
		cur_patch.m_edge_index_array.push_back(merged_edge_idx);
		cur_patch.m_conner_index_array.clear();
		cur_patch.m_conner_index_array.push_back(merged_conner_idx);
		cur_patch.m_conner_index_array.push_back(third_conner_idx);
		cur_patch.m_conner_index_array.push_back(add_conner_index);
		cur_patch.m_conner_index_array.push_back(single_conner_idx);

		std::vector<int> patch_boundary;
		FormPatchBoundary(patch_id, patch_boundary);
		FindInnerFace(p_mesh, patch_boundary, cur_patch.m_face_index_array, m_half_edge);
		FormPatchBoundary(add_patch_index, patch_boundary);
		FindInnerFace(p_mesh, patch_boundary, m_patch_array[add_patch_index].m_face_index_array, m_half_edge);

		
	}

	void ChartCreator::ReConnectDegradePatch(int patch_id)
	{
		std::cout << "Reconnect degrade patch " << patch_id << std::endl;
		ParamPatch& patch = m_patch_array[patch_id];
		const std::vector<int>& patch_edge_vec = patch.m_edge_index_array;
		int merged_edge_idx(-1),prev_edge_idx(-1), next_edge_idx(-1);
		for(size_t k=0; k<patch_edge_vec.size(); ++k){
			for(size_t i=k+1; i<patch_edge_vec.size(); ++i){
				if(patch_edge_vec[k] == patch_edge_vec[i]){
					merged_edge_idx = patch_edge_vec[i]; 
					prev_edge_idx = patch_edge_vec[(k+patch_edge_vec.size()-1)%patch_edge_vec.size()];
					next_edge_idx = patch_edge_vec[(i+1)%patch_edge_vec.size()];
					break; 
				}
			}
			if(merged_edge_idx != -1) break;
		}


		int merged_conner_idx(-1), single_conner_idx(-1), third_conner_idx(-1);
		const PatchEdge& merged_patch_edge = m_patch_edge_array[merged_edge_idx];
		int temp_conner_idx = merged_patch_edge.m_conner_pair_index.first;
		const PatchConner& temp_conner = m_patch_conner_array[temp_conner_idx];
		if(temp_conner.m_nb_conner_index_array.size() == 1) {
			merged_conner_idx = merged_patch_edge.m_conner_pair_index.second;
			single_conner_idx = merged_patch_edge.m_conner_pair_index.first;
		}else{
			merged_conner_idx = merged_patch_edge.m_conner_pair_index.first;
			single_conner_idx = merged_patch_edge.m_conner_pair_index.second;
		}

		const PatchEdge& prev_patch_edge = m_patch_edge_array[prev_edge_idx];
		if(prev_patch_edge.m_conner_pair_index.first == merged_conner_idx){
			third_conner_idx = prev_patch_edge.m_conner_pair_index.second;
		}else{
			third_conner_idx = prev_patch_edge.m_conner_pair_index.first;
		}
		const PatchEdge& next_patch_edge = m_patch_edge_array[next_edge_idx];

		const PatchConner& merged_conner = m_patch_conner_array[merged_conner_idx];
		const PatchConner& single_conner = m_patch_conner_array[single_conner_idx];
		const PatchConner& third_conner = m_patch_conner_array[third_conner_idx];
		

		int pre_patch_idx(-1), nxt_patch_idx(-1);
		for(size_t k=0; k<prev_patch_edge.m_nb_patch_index_array.size(); ++k){
			if(prev_patch_edge.m_nb_patch_index_array[k] != patch_id) 
				pre_patch_idx = prev_patch_edge.m_nb_patch_index_array[k];
		}

		for(size_t k=0; k<next_patch_edge.m_nb_patch_index_array.size(); ++k){
			if(next_patch_edge.m_nb_patch_index_array[k] != patch_id)
				nxt_patch_idx = next_patch_edge.m_nb_patch_index_array[k];
		}
		assert(pre_patch_idx != -1 && nxt_patch_idx != -1);


		const ParamPatch& pre_patch = m_patch_array[pre_patch_idx];		
		const ParamPatch& nxt_patch = m_patch_array[nxt_patch_idx];		
		
		int sad1_conner_idx(-1), sad1_third_edge_idx(-1), sad2_conner_idx(-1), sad2_third_edge_idx(-1);
		int pos = find(third_conner.m_nb_edge_index_array.begin(), third_conner.m_nb_edge_index_array.end(), prev_edge_idx)
			- third_conner.m_nb_edge_index_array.begin();
		sad1_conner_idx = third_conner.m_nb_conner_index_array[(pos+third_conner.m_nb_conner_index_array.size()-1)%third_conner.m_nb_conner_index_array.size()];
		sad1_third_edge_idx = third_conner.m_nb_edge_index_array[(pos+third_conner.m_nb_edge_index_array.size()-1)%third_conner.m_nb_edge_index_array.size()];
		pos = find(third_conner.m_nb_edge_index_array.begin(), third_conner.m_nb_edge_index_array.end(), next_edge_idx)
			- third_conner.m_nb_edge_index_array.begin();
		sad2_conner_idx = third_conner.m_nb_conner_index_array[(pos+1)%third_conner.m_nb_conner_index_array.size()];
		sad2_third_edge_idx = third_conner.m_nb_edge_index_array[(pos+1)%third_conner.m_nb_edge_index_array.size()];

		if(sad1_conner_idx == sad2_conner_idx)
		{
			const PatchConner& sad_conner = m_patch_conner_array[sad1_conner_idx];
			int pos1 = find(sad_conner.m_nb_edge_index_array.begin(), sad_conner.m_nb_edge_index_array.end(), sad1_third_edge_idx)
				- sad_conner.m_nb_edge_index_array.begin();
			int pos2 = find(sad_conner.m_nb_edge_index_array.begin(), sad_conner.m_nb_edge_index_array.end(), sad2_third_edge_idx)
				- sad_conner.m_nb_edge_index_array.begin();

			if( (pos1 + 2) % sad_conner.m_nb_edge_index_array.size() == pos2)
			{
				int sad_nb_num = sad_conner.m_nb_conner_index_array.size();
				int single_conner_idx_2 = sad_conner.m_nb_conner_index_array[(pos1+1)%sad_nb_num];
				int merger_edge_idx_2 = sad_conner.m_nb_edge_index_array[(pos1+1)%sad_nb_num];				
				int sp_conner_idx = sad_conner.m_nb_conner_index_array[(pos1+sad_nb_num-1)%sad_nb_num];
				int sn_conner_idx = sad_conner.m_nb_conner_index_array[(pos2+1)%sad_nb_num];
				int sp_edge_idx = sad_conner.m_nb_edge_index_array[(pos1+sad_nb_num-1)%sad_nb_num];
				int sn_edge_idx = sad_conner.m_nb_edge_index_array[(pos2+1)%sad_nb_num];
				int mid_patch_id(-1);
				const PatchEdge& sp_edge = m_patch_edge_array[sad1_third_edge_idx];
				if(sp_edge.m_nb_patch_index_array[0] == pre_patch_idx) mid_patch_id = sp_edge.m_nb_patch_index_array[1];
				else mid_patch_id = sp_edge.m_nb_patch_index_array[0];

				/// modified 
				PatchConner& md_conner_1 = m_patch_conner_array[third_conner_idx];
				md_conner_1.m_nb_conner_index_array.clear(); md_conner_1.m_nb_conner_index_array.resize(4);
				md_conner_1.m_nb_edge_index_array.clear(); md_conner_1.m_nb_edge_index_array.resize(4);
				md_conner_1.m_nb_conner_index_array[0] = sp_conner_idx; md_conner_1.m_nb_edge_index_array[0] = sad1_third_edge_idx;
				md_conner_1.m_nb_conner_index_array[1] = single_conner_idx; md_conner_1.m_nb_edge_index_array[1] = prev_edge_idx;
				md_conner_1.m_nb_conner_index_array[2] = sn_conner_idx; md_conner_1.m_nb_edge_index_array[2] = next_edge_idx;
				md_conner_1.m_nb_conner_index_array[3] = single_conner_idx_2; md_conner_1.m_nb_edge_index_array[3] = sad2_third_edge_idx;

				PatchConner& md_conner_2 = m_patch_conner_array[single_conner_idx];
				md_conner_2.m_nb_conner_index_array.clear(); md_conner_2.m_nb_conner_index_array.resize(2);
				md_conner_2.m_nb_edge_index_array.clear(); md_conner_2.m_nb_edge_index_array.resize(2);
				md_conner_2.m_nb_conner_index_array[0] = third_conner_idx; md_conner_2.m_nb_edge_index_array[0] = prev_edge_idx;
				md_conner_2.m_nb_conner_index_array[1] = merged_conner_idx; md_conner_2.m_nb_edge_index_array[1] = merged_edge_idx;

				PatchConner& md_conner_3 = m_patch_conner_array[single_conner_idx_2];
				md_conner_3.m_nb_conner_index_array.clear(); md_conner_3.m_nb_conner_index_array.resize(2);
				md_conner_3.m_nb_edge_index_array.clear(); md_conner_3.m_nb_edge_index_array.resize(2);
				md_conner_3.m_nb_conner_index_array[0] = third_conner_idx; md_conner_3.m_nb_edge_index_array[0] = sad2_third_edge_idx;
				md_conner_3.m_nb_conner_index_array[1] = sad1_conner_idx; md_conner_3.m_nb_edge_index_array[1] = merger_edge_idx_2;

				PatchConner& md_conner_4 = m_patch_conner_array[merged_conner_idx];
				pos = find(md_conner_4.m_nb_edge_index_array.begin(), md_conner_4.m_nb_edge_index_array.end(), prev_edge_idx)
					- md_conner_4.m_nb_edge_index_array.begin();
				md_conner_4.m_nb_conner_index_array.erase(md_conner_4.m_nb_conner_index_array.begin() + pos);
				md_conner_4.m_nb_edge_index_array.erase(md_conner_4.m_nb_edge_index_array.begin() + pos);
				pos = find(md_conner_4.m_nb_edge_index_array.begin(), md_conner_4.m_nb_edge_index_array.end(), next_edge_idx)
					- md_conner_4.m_nb_edge_index_array.begin();
				md_conner_4.m_nb_conner_index_array.erase(md_conner_4.m_nb_conner_index_array.begin() + pos);
				md_conner_4.m_nb_edge_index_array.erase(md_conner_4.m_nb_edge_index_array.begin() + pos);

				PatchConner& md_conner_5 = m_patch_conner_array[sad1_conner_idx];
				pos = find(md_conner_5.m_nb_edge_index_array.begin(), md_conner_5.m_nb_edge_index_array.end(), sad1_third_edge_idx)
					- md_conner_5.m_nb_edge_index_array.begin();
				md_conner_5.m_nb_conner_index_array.erase(md_conner_5.m_nb_conner_index_array.begin() + pos);
				md_conner_5.m_nb_edge_index_array.erase(md_conner_5.m_nb_edge_index_array.begin() + pos);
				pos = find(md_conner_5.m_nb_edge_index_array.begin(), md_conner_5.m_nb_edge_index_array.end(), sad2_third_edge_idx)
					- md_conner_5.m_nb_edge_index_array.begin();
				md_conner_5.m_nb_conner_index_array.erase(md_conner_5.m_nb_conner_index_array.begin() + pos);
				md_conner_5.m_nb_edge_index_array.erase(md_conner_5.m_nb_edge_index_array.begin() + pos);

				PatchConner& md_conner_6 = m_patch_conner_array[sp_conner_idx];
				pos = find(md_conner_6.m_nb_edge_index_array.begin(), md_conner_6.m_nb_edge_index_array.end(), sp_edge_idx)
					- md_conner_6.m_nb_edge_index_array.begin();
				md_conner_6.m_nb_edge_index_array.insert(md_conner_6.m_nb_edge_index_array.begin() + pos, sad1_third_edge_idx);
				md_conner_6.m_nb_conner_index_array.insert(md_conner_6.m_nb_conner_index_array.begin() + pos, third_conner_idx);

				PatchConner& md_conner_7 = m_patch_conner_array[sn_conner_idx];
				pos = find(md_conner_7.m_nb_edge_index_array.begin(), md_conner_7.m_nb_edge_index_array.end(), sn_edge_idx)
					- md_conner_7.m_nb_edge_index_array.begin() + 1;
				md_conner_7.m_nb_edge_index_array.insert(md_conner_7.m_nb_edge_index_array.begin() + pos, next_edge_idx);
				md_conner_7.m_nb_conner_index_array.insert(md_conner_7.m_nb_conner_index_array.begin() + pos, third_conner_idx);


				PatchEdge& modified_edge1 = m_patch_edge_array[sad1_third_edge_idx];
				modified_edge1.m_conner_pair_index = std::make_pair(sp_conner_idx, third_conner_idx);
				SetPatchConners(pre_patch_idx);
				FindShortestPathInPatch(pre_patch_idx, m_patch_conner_array[sp_conner_idx].m_mesh_index, third_conner.m_mesh_index, modified_edge1.m_mesh_path);

				PatchEdge& modified_edge2 = m_patch_edge_array[sad2_third_edge_idx];
				modified_edge2.m_conner_pair_index = std::make_pair(single_conner_idx_2, third_conner_idx);
				modified_edge2.m_nb_patch_index_array[0] = mid_patch_id; 
				modified_edge2.m_nb_patch_index_array[1] = patch_id;
				SetPatchConners(mid_patch_id);
				FindShortestPathInPatch(mid_patch_id, m_patch_conner_array[single_conner_idx_2].m_mesh_index, third_conner.m_mesh_index, modified_edge2.m_mesh_path);

				PatchEdge& modified_edge3 = m_patch_edge_array[prev_edge_idx];
				modified_edge3.m_conner_pair_index = std::make_pair(single_conner_idx, third_conner_idx);
				modified_edge3.m_nb_patch_index_array[0] = pre_patch_idx;
				modified_edge3.m_nb_patch_index_array[1] = nxt_patch_idx;
				SetPatchConners(patch_id);
				FindShortestPathInPatch(patch_id, single_conner.m_mesh_index, third_conner.m_mesh_index, modified_edge3.m_mesh_path);

				PatchEdge& modified_edge4 = m_patch_edge_array[next_edge_idx];
				modified_edge4.m_conner_pair_index = std::make_pair(sn_conner_idx, third_conner_idx);
				SetPatchConners(nxt_patch_idx);
				FindShortestPathInPatch(nxt_patch_idx, m_patch_conner_array[sn_conner_idx].m_mesh_index, third_conner.m_mesh_index, modified_edge4.m_mesh_path);

				PatchEdge& modified_edge5 = m_patch_edge_array[merged_edge_idx];
				modified_edge5.m_nb_patch_index_array[0] = pre_patch_idx;
				modified_edge5.m_nb_patch_index_array[1] = nxt_patch_idx;

				PatchEdge& modified_edge6 = m_patch_edge_array[merger_edge_idx_2];
				modified_edge6.m_nb_patch_index_array[0] = mid_patch_id;
				modified_edge6.m_nb_patch_index_array[1] = patch_id;

				PatchEdge& modified_edge7 = m_patch_edge_array[sp_edge_idx];
				for(size_t k=0; k<modified_edge7.m_nb_patch_index_array.size(); ++k) 
					if(modified_edge7.m_nb_patch_index_array[k] == pre_patch_idx) modified_edge7.m_nb_patch_index_array[k] = mid_patch_id;

				PatchEdge& modified_edge8 = m_patch_edge_array[sn_edge_idx];
				for(size_t k=0; k<modified_edge8.m_nb_patch_index_array.size(); ++k)
					if(modified_edge8.m_nb_patch_index_array[k] == nxt_patch_idx) modified_edge8.m_nb_patch_index_array[k] = patch_id;


				ParamPatch& md_patch1 = m_patch_array[mid_patch_id];
				md_patch1.m_conner_index_array.clear(); md_patch1.m_conner_index_array.resize(4);
				md_patch1.m_edge_index_array.clear(); md_patch1.m_edge_index_array.resize(4);
				md_patch1.m_conner_index_array[0] = sp_conner_idx; md_patch1.m_edge_index_array[0] = sad1_third_edge_idx;
				md_patch1.m_conner_index_array[1] = third_conner_idx; md_patch1.m_edge_index_array[1] = sad2_third_edge_idx;
				md_patch1.m_conner_index_array[2] = single_conner_idx_2, md_patch1.m_edge_index_array[2] = merger_edge_idx_2;
				md_patch1.m_conner_index_array[3] = sad1_conner_idx; md_patch1.m_edge_index_array[3] = sp_edge_idx;
				FindPatchInnerFace(mid_patch_id);

				ParamPatch& md_patch2 = m_patch_array[patch_id];
				md_patch2.m_conner_index_array.clear(); md_patch2.m_conner_index_array.resize(4);
				md_patch2.m_edge_index_array.clear(); md_patch2.m_edge_index_array.resize(4);
				md_patch2.m_conner_index_array[0] = sn_conner_idx; md_patch2.m_edge_index_array[0] = sn_edge_idx;
				md_patch2.m_conner_index_array[1] = sad1_conner_idx; md_patch2.m_edge_index_array[1] = merger_edge_idx_2;
				md_patch2.m_conner_index_array[2] = single_conner_idx_2; md_patch2.m_edge_index_array[2] = sad2_third_edge_idx;
				md_patch2.m_conner_index_array[3] = third_conner_idx; md_patch2.m_edge_index_array[3] = next_edge_idx;
				FindPatchInnerFace(single_conner_idx);

				ParamPatch& md_patch3 = m_patch_array[pre_patch_idx];
				for(size_t k=0; k<md_patch3.m_edge_index_array.size(); ++k){
					if(md_patch3.m_edge_index_array[k] == prev_edge_idx){
						md_patch3.m_edge_index_array[k] = merged_edge_idx;
						md_patch3.m_conner_index_array[k] = merged_conner_idx;
					}else if(md_patch3.m_edge_index_array[k] == sad1_third_edge_idx){
						md_patch3.m_edge_index_array[k] = prev_edge_idx;
						md_patch3.m_conner_index_array[k] = single_conner_idx;
					}else if(md_patch3.m_edge_index_array[k] == sp_edge_idx){
						md_patch3.m_edge_index_array[k] = sad1_third_edge_idx;
						md_patch3.m_conner_index_array[k] = third_conner_idx;
					}
				}
				FindPatchInnerFace(pre_patch_idx);

				ParamPatch& md_patch4 = m_patch_array[nxt_patch_idx];
				for(size_t k=0; k<md_patch4.m_edge_index_array.size(); ++k){
					if(md_patch4.m_edge_index_array[k] == sn_edge_idx){
						md_patch4.m_edge_index_array[k] = next_edge_idx;
						md_patch4.m_conner_index_array[k] = sn_conner_idx;
					}else if(md_patch4.m_edge_index_array[k] == sad2_third_edge_idx){
						md_patch4.m_edge_index_array[k] = prev_edge_idx;
						md_patch4.m_conner_index_array[k] = third_conner_idx;
					}else if(md_patch4.m_edge_index_array[k] == next_edge_idx){
						md_patch4.m_edge_index_array[k] = merged_edge_idx;
						md_patch4.m_conner_index_array[k] = single_conner_idx;
					}
				}
				FindPatchInnerFace(nxt_patch_idx);
			}
			return;
		}

		const PatchConner& sad1_conner = m_patch_conner_array[sad1_conner_idx];
		const PatchConner& sad2_conner = m_patch_conner_array[sad2_conner_idx];

		SetPatchConners(patch_id);
		SetPatchConners(pre_patch_idx);
		SetPatchConners(nxt_patch_idx);

		/// modify patch /edge / conner
		PatchEdge& modified_edge_1 = m_patch_edge_array[prev_edge_idx]; /// modify the prev edge
		modified_edge_1.m_conner_pair_index = std::make_pair(sad1_conner_idx, single_conner_idx);
		FindShortestPathInTwoAdjPatch(pre_patch_idx, patch_id, sad1_conner.m_mesh_index, single_conner.m_mesh_index, modified_edge_1.m_mesh_path);
		PatchEdge& modified_edge_2 = m_patch_edge_array[next_edge_idx]; /// modify the next edge
		modified_edge_2.m_conner_pair_index = std::make_pair(sad2_conner_idx, single_conner_idx);
		FindShortestPathInTwoAdjPatch(nxt_patch_idx, patch_id, sad2_conner.m_mesh_index, single_conner.m_mesh_index, modified_edge_2.m_mesh_path);

		PatchEdge& modified_edge_3 = m_patch_edge_array[sad1_third_edge_idx];
		for(size_t k=0; k<modified_edge_3.m_nb_patch_index_array.size(); ++k){
			if(modified_edge_3.m_nb_patch_index_array[k] == pre_patch_idx) modified_edge_3.m_nb_patch_index_array[k] = patch_id;
		}
		PatchEdge& modified_edge_4 = m_patch_edge_array[sad2_third_edge_idx];
		for(size_t k=0; k<modified_edge_4.m_nb_patch_index_array.size(); ++k){
			if(modified_edge_4.m_nb_patch_index_array[k] == nxt_patch_idx) modified_edge_4.m_nb_patch_index_array[k] = patch_id;
		}
		PatchEdge& modified_edge_5 = m_patch_edge_array[merged_edge_idx];
		modified_edge_5.m_nb_patch_index_array.clear();
		modified_edge_5.m_nb_patch_index_array.push_back(pre_patch_idx); 
		modified_edge_5.m_nb_patch_index_array.push_back(nxt_patch_idx);

		////
		PatchConner& modified_conner_1 = m_patch_conner_array[third_conner_idx]; 
		pos = find(modified_conner_1.m_nb_edge_index_array.begin(), modified_conner_1.m_nb_edge_index_array.end(), 
			prev_edge_idx) - modified_conner_1.m_nb_edge_index_array.begin();
		modified_conner_1.m_nb_conner_index_array.erase(modified_conner_1.m_nb_conner_index_array.begin() + pos);
		modified_conner_1.m_nb_edge_index_array.erase(modified_conner_1.m_nb_edge_index_array.begin() + pos);
		pos = find(modified_conner_1.m_nb_edge_index_array.begin(), modified_conner_1.m_nb_edge_index_array.end(), 
			next_edge_idx) - modified_conner_1.m_nb_edge_index_array.begin();
		modified_conner_1.m_nb_conner_index_array.erase(modified_conner_1.m_nb_conner_index_array.begin() + pos);
		modified_conner_1.m_nb_edge_index_array.erase(modified_conner_1.m_nb_edge_index_array.begin() + pos);

		PatchConner& modified_conner_2 = m_patch_conner_array[single_conner_idx];
		modified_conner_2.m_nb_conner_index_array.push_back(sad2_conner_idx);
		modified_conner_2.m_nb_edge_index_array.push_back(next_edge_idx);
		modified_conner_2.m_nb_conner_index_array.push_back(sad1_conner_idx);
		modified_conner_2.m_nb_edge_index_array.push_back(prev_edge_idx);

		PatchConner& modified_conner_3 = m_patch_conner_array[merged_conner_idx];
		pos = find(modified_conner_3.m_nb_edge_index_array.begin(), modified_conner_3.m_nb_edge_index_array.end(), 
			prev_edge_idx) - modified_conner_3.m_nb_edge_index_array.begin();
		modified_conner_3.m_nb_conner_index_array.erase(modified_conner_3.m_nb_conner_index_array.begin() + pos);
		modified_conner_3.m_nb_edge_index_array.erase(modified_conner_3.m_nb_edge_index_array.begin() + pos);
		pos = find(modified_conner_3.m_nb_edge_index_array.begin(), modified_conner_3.m_nb_edge_index_array.end(), 
			next_edge_idx) - modified_conner_3.m_nb_edge_index_array.begin();
		modified_conner_3.m_nb_conner_index_array.erase(modified_conner_3.m_nb_conner_index_array.begin() + pos);
		modified_conner_3.m_nb_edge_index_array.erase(modified_conner_3.m_nb_edge_index_array.begin() + pos);

		PatchConner& modified_conner_4 = m_patch_conner_array[sad1_conner_idx];
		pos = find(modified_conner_4.m_nb_edge_index_array.begin(), modified_conner_4.m_nb_edge_index_array.end(),
			sad1_third_edge_idx) - modified_conner_4.m_nb_edge_index_array.begin()  ;
		modified_conner_4.m_nb_edge_index_array.insert(modified_conner_4.m_nb_edge_index_array.begin() + pos, prev_edge_idx);
		modified_conner_4.m_nb_conner_index_array.insert(modified_conner_4.m_nb_conner_index_array.begin() + pos, single_conner_idx);

		PatchConner& modified_conner_5 = m_patch_conner_array[sad2_conner_idx];
		pos = find(modified_conner_5.m_nb_edge_index_array.begin(), modified_conner_5.m_nb_edge_index_array.end(), 
			sad2_third_edge_idx) - modified_conner_5.m_nb_edge_index_array.begin() + 1;
		modified_conner_5.m_nb_edge_index_array.insert(modified_conner_5.m_nb_edge_index_array.begin() + pos, next_edge_idx);
		modified_conner_5.m_nb_conner_index_array.insert(modified_conner_5.m_nb_conner_index_array.begin() + pos, single_conner_idx);


		////
		ParamPatch& modified_patch_1 = m_patch_array[patch_id];
		modified_patch_1.m_conner_index_array.clear();
		modified_patch_1.m_conner_index_array.push_back(third_conner_idx);
		modified_patch_1.m_conner_index_array.push_back(sad1_conner_idx);
		modified_patch_1.m_conner_index_array.push_back(single_conner_idx);
		modified_patch_1.m_conner_index_array.push_back(sad2_conner_idx);
		modified_patch_1.m_edge_index_array.clear();
		modified_patch_1.m_edge_index_array.push_back(sad1_third_edge_idx);
		modified_patch_1.m_edge_index_array.push_back(prev_edge_idx);
		modified_patch_1.m_edge_index_array.push_back(next_edge_idx);
		modified_patch_1.m_edge_index_array.push_back(sad2_third_edge_idx);
		FindPatchInnerFace(patch_id);

		ParamPatch& modified_patch_2 = m_patch_array[pre_patch_idx];
		for(size_t k=0; k<modified_patch_2.m_edge_index_array.size(); ++k){
			if(modified_patch_2.m_edge_index_array[k] == sad1_third_edge_idx){
				modified_patch_2.m_edge_index_array[k] = prev_edge_idx;
				modified_patch_2.m_conner_index_array[k] = single_conner_idx;
			}else if(modified_patch_2.m_edge_index_array[k] == prev_edge_idx){
				modified_patch_2.m_edge_index_array[k] = merged_edge_idx;
			}
		}
		FindPatchInnerFace(pre_patch_idx);
		
		ParamPatch& modified_patch_3 = m_patch_array[nxt_patch_idx];
		for(size_t k=0; k<modified_patch_3.m_edge_index_array.size(); ++k){
			if(modified_patch_3.m_edge_index_array[k] == sad2_third_edge_idx){
				modified_patch_3.m_edge_index_array[k] = next_edge_idx;
			}else if(modified_patch_3.m_edge_index_array[k] == next_edge_idx){
				modified_patch_3.m_edge_index_array[k] = merged_edge_idx;
				modified_patch_3.m_conner_index_array[k] = single_conner_idx;
			}
		}		
		FindPatchInnerFace(nxt_patch_idx);		
	}

	void ChartCreator::OptimizeAmbiguityPatchShape()
	{
		std::set< std::pair<int, int> > oped_pair;
		for(size_t k=0; k<m_patch_array.size(); ++k){
			const ParamPatch& patch = m_patch_array[k];
			const std::vector<int>& nb_patch_vec = patch.m_nb_patch_index_array;
			for(size_t i=0; i<nb_patch_vec.size(); ++i){
				int nb_pid = nb_patch_vec[i];
				if(oped_pair.find(std::make_pair(k, nb_pid)) == oped_pair.end() &&
					oped_pair.find(std::make_pair(nb_pid, k)) == oped_pair.end())
				{
					OptimizeAmbiguityPatchPairShape(k, nb_pid);
					oped_pair.insert(std::make_pair(k, nb_pid));
				}
			}
		}
	}

    void ChartCreator::OptimizeAmbiguityPatchPairShape(int patch_id_1, int patch_id_2)
    {
        ParamPatch& param_patch_1 = m_patch_array[patch_id_1];
        ParamPatch& param_patch_2 = m_patch_array[patch_id_2];

		/// find the two common edges
		std::vector<int> common_edges;
		const std::vector<int>& patch_edge_vec_1 = param_patch_1.m_edge_index_array;
		const std::vector<int>& patch_edge_vec_2 = param_patch_2.m_edge_index_array;
		for(size_t k=0; k<patch_edge_vec_1.size(); ++k){
			if(find(patch_edge_vec_2.begin(), patch_edge_vec_2.end(), patch_edge_vec_1[k]) 
				!= patch_edge_vec_2.end()) common_edges.push_back(patch_edge_vec_1[k]);
		}
		if(common_edges.size() !=2 ) return;

		std::cout << "Optimize Ambiguity patch " << patch_id_1 << " " << patch_id_2 << std::endl;

		int com_edge_idx_1 = common_edges[0], com_edge_idx_2 = common_edges[1], com_edge_idx_3(-1);
		PatchEdge& patch_edge_1 = m_patch_edge_array[com_edge_idx_1];
		PatchEdge& patch_edge_2 = m_patch_edge_array[com_edge_idx_2];

		bool flag(false);
		if(m_patch_conner_array[patch_edge_1.m_conner_pair_index.first].m_conner_type == 1
			|| m_patch_conner_array[patch_edge_1.m_conner_pair_index.second].m_conner_type == 1){
				flag = true;
		}
		com_edge_idx_3 = flag ? com_edge_idx_1 : com_edge_idx_2;

		const std::vector<int>& inner_path = m_patch_edge_array[com_edge_idx_3].m_mesh_path;

		/// find the common region of these two patchs
		const std::vector<int>& face_vec_1 = param_patch_1.m_face_index_array;
		const std::vector<int>& face_vec_2 = param_patch_2.m_face_index_array;
		const PolyIndexArray& face_list_array = p_mesh->m_Kernel.GetFaceInfo().GetIndex();

		std::set< std::pair<int, int> > region_mesh_edge_set;
		for(size_t k=0; k<patch_edge_vec_1.size(); ++k){
			const std::vector<int>& path = m_patch_edge_array[patch_edge_vec_2[k]].m_mesh_path;
			for(size_t i=1; i<path.size(); ++i){
				region_mesh_edge_set.insert(MakeEdge(path[i-1], path[i]));
			}
		}
		for(size_t k=0; k<patch_edge_vec_2.size(); ++k){			
			const std::vector<int>& path = m_patch_edge_array[patch_edge_vec_2[k]].m_mesh_path;
			for(size_t i=1; i<path.size(); ++i){
				region_mesh_edge_set.insert(MakeEdge(path[i-1],path[i]));
			}
		}

		for(size_t k=0; k<face_vec_1.size(); ++k){
			int fid = face_vec_1[k];
			const IndexArray& face = face_list_array[fid];
			for(int i=0; i<3; ++i){
				if(find(inner_path.begin(), inner_path.end(), face[(i+1)%3]) != inner_path.end() &&
					find(inner_path.begin(), inner_path.end(), face[i]) != inner_path.end()) continue;
				region_mesh_edge_set.insert(MakeEdge(face[i], face[(i+1)%3]));
			}
		}
		for(size_t k=0; k<face_vec_2.size(); ++k){
			int fid = face_vec_2[k];
			const IndexArray& face = face_list_array[fid];
			for(int i=0; i<3; ++i){
				if(find(inner_path.begin(), inner_path.end(), face[(i+1)%3]) != inner_path.end() &&
					find(inner_path.begin(), inner_path.end(), face[i]) != inner_path.end()) continue;
				region_mesh_edge_set.insert(MakeEdge(face[i], face[(i+1)%3]));
			}
		}

		int vid1, vid2;
		if(!flag){
			vid1 = m_patch_conner_array[patch_edge_1.m_conner_pair_index.first].m_mesh_index;
			vid2 = m_patch_conner_array[patch_edge_1.m_conner_pair_index.second].m_mesh_index;
			FindShortestPathInRegion(p_mesh, vid1, vid2, region_mesh_edge_set, patch_edge_1.m_mesh_path);
		}else{
			vid1 = m_patch_conner_array[patch_edge_2.m_conner_pair_index.first].m_mesh_index;
			vid2 = m_patch_conner_array[patch_edge_2.m_conner_pair_index.second].m_mesh_index;
			FindShortestPathInRegion(p_mesh, vid1, vid2, region_mesh_edge_set, patch_edge_2.m_mesh_path);
		}
		FindPatchInnerFace(patch_id_1);
		FindPatchInnerFace(patch_id_2);
    }
    
	void ChartCreator::FindPatchInnerFace(int patch_id)
	{
		ParamPatch& patch = m_patch_array[patch_id];

		std::vector<int> patch_bounary;
		FormPatchBoundary(patch_id, patch_bounary);
		FindInnerFace(p_mesh, patch_bounary, patch.m_face_index_array, m_half_edge);
	}

	void ChartCreator::FormPatchBoundary(int patch_id, std::vector<int>& boundary) const
	{
		boundary.clear();
		const ParamPatch& patch = m_patch_array[patch_id];
		const std::vector<int>& conner_index_vec = patch.m_conner_index_array;
		const std::vector<int>& edge_index_vec = patch.m_edge_index_array;

		std::vector < std::vector<int> > path_vec(edge_index_vec.size());
		for(size_t k=0; k<edge_index_vec.size(); ++k){
			const PatchEdge& patch_edge = m_patch_edge_array[edge_index_vec[k]];
			const std::vector<int>& path = patch_edge.m_mesh_path;
			if(patch_edge.m_conner_pair_index.first == conner_index_vec[k]){				
				for(size_t i=0; i<path.size()-1; ++i){
					boundary.push_back(path[i]);
				}
			}else{
				for(int i=(int)path.size()-1; i>0; --i){
					boundary.push_back(path[i]);
				}
			}
		}
		boundary.push_back(boundary[0]);
		reverse(boundary.begin(), boundary.end());
	}

	void ChartCreator::ValenceControl()
	{	
		int prev_hightest_val = -1;
		m_unre_edge_index_array.clear();
		while(true)
		{
			int highest_valence(-1);
			int highest_valence_conner_index(-1);
			for(size_t k=0; k<m_patch_conner_array.size(); ++k)
			{
				int cur_conner_val = 0;
				std::vector<int>& conner_nb_edges = m_patch_conner_array[k].m_nb_edge_index_array;
				for(size_t i=0; i<conner_nb_edges.size(); ++i)
				{
					if(find(m_unre_edge_index_array.begin(), m_unre_edge_index_array.end(), conner_nb_edges[i])
						== m_unre_edge_index_array.end()) cur_conner_val ++;
				}
				if(cur_conner_val > highest_valence)
				{
					highest_valence = cur_conner_val;
					highest_valence_conner_index = k;
				}
			}
			std::cout << "Highest Valence : " << highest_valence << std::endl;
			
			if(highest_valence <= 6 || highest_valence == prev_hightest_val) break;
			prev_hightest_val = highest_valence;
			int edge_idx = FindLongestPatchEdge(highest_valence_conner_index);
			RemovePatchEdgeToReduceValance(edge_idx, highest_valence_conner_index);
		}
	}

	int ChartCreator::FindLongestPatchEdge(int conner_index) const
	{
		const PatchConner& patch_conner = m_patch_conner_array[conner_index];
		int longest_edge_index(-1), longest_edge_num(-1);
		for(size_t k=0; k<patch_conner.m_nb_edge_index_array.size(); ++k){
			int edge_idx = patch_conner.m_nb_edge_index_array[k];
			if(find(m_unre_edge_index_array.begin(), m_unre_edge_index_array.end(), edge_idx) 
				!= m_unre_edge_index_array.end()) continue;
			const PatchEdge& patch_edge = m_patch_edge_array[edge_idx];
			if(patch_edge.m_nb_patch_index_array.size() != 2){ 
				/// don't consider boundary edge
				continue;
			}
			int edge_num = (int)patch_edge.m_mesh_path.size();
			if(edge_num >longest_edge_num){
				longest_edge_num = edge_num;
				longest_edge_index = edge_idx;
			}
		}
		return longest_edge_index;
	}
	
	void ChartCreator::RemovePatchEdgeToReduceValance(int rm_edge_idx, int high_val_conner_idx)
	{
		const PatchConner& high_val_conner = m_patch_conner_array[high_val_conner_idx];
		const PatchEdge& rm_edge = m_patch_edge_array[rm_edge_idx];
		assert(rm_edge.m_nb_patch_index_array.size() == 2);
		
		/// the prev and next is order by the conners ccw and the rm edge
		int pre_edge_idx(-1), nxt_edge_idx(-1);
		int nb_edge_num = (int) high_val_conner.m_nb_edge_index_array.size();
		for(size_t k=0; k<high_val_conner.m_nb_edge_index_array.size(); ++k){
			if(high_val_conner.m_nb_edge_index_array[k] == rm_edge_idx){
				pre_edge_idx = high_val_conner.m_nb_edge_index_array[(k+nb_edge_num-1)%nb_edge_num];
				nxt_edge_idx = high_val_conner.m_nb_edge_index_array[(k+1)%nb_edge_num];
				break;
			}
		}
		assert(pre_edge_idx != -1 && nxt_edge_idx != -1);
		const PatchEdge& pre_edge = m_patch_edge_array[pre_edge_idx];
		const PatchEdge& nxt_edge = m_patch_edge_array[nxt_edge_idx];

		int nb_patch_idx_1(rm_edge.m_nb_patch_index_array[0]);
		int nb_patch_idx_2(rm_edge.m_nb_patch_index_array[1]);
		int pre_patch_idx(-1), nxt_patch_idx(-1);
		const ParamPatch& nb_patch_1 = m_patch_array[nb_patch_idx_1];
		const ParamPatch& nb_patch_2 = m_patch_array[nb_patch_idx_2];

		if(find(nb_patch_1.m_edge_index_array.begin(), nb_patch_1.m_edge_index_array.end(), pre_edge_idx)
			!= nb_patch_1.m_edge_index_array.end()){
				pre_patch_idx = nb_patch_idx_1;
				nxt_patch_idx = nb_patch_idx_2;
		}else{
			pre_patch_idx = nb_patch_idx_2;
			nxt_patch_idx = nb_patch_idx_1;
		}
		assert(pre_patch_idx != -1 && nxt_patch_idx != -1);
		const ParamPatch& pre_patch = m_patch_array[pre_patch_idx];
		const ParamPatch& nxt_patch = m_patch_array[nxt_patch_idx];

		int pre_conner_idx(-1), nxt_conner_idx(-1), mid_conner_idx(-1);
		pre_conner_idx = (pre_edge.m_conner_pair_index.first == high_val_conner_idx) ? 
			pre_edge.m_conner_pair_index.second : pre_edge.m_conner_pair_index.first;
		nxt_conner_idx = (nxt_edge.m_conner_pair_index.first == high_val_conner_idx) ?
			nxt_edge.m_conner_pair_index.second : nxt_edge.m_conner_pair_index.first;
		mid_conner_idx = (rm_edge.m_conner_pair_index.first == high_val_conner_idx) ?
			rm_edge.m_conner_pair_index.second : rm_edge.m_conner_pair_index.first;

		if(pre_conner_idx == nxt_conner_idx){
			m_unre_edge_index_array.push_back(rm_edge_idx);
			return;
		}

		int mid_vert_idx = rm_edge.m_mesh_path.size() / 2;
		int mid_vert = rm_edge.m_mesh_path[mid_vert_idx];

		/// we will add one conner, two edges and one patch
		int add_conner_idx = (int) m_patch_conner_array.size();
		int add_edge_idx_1 = (int) m_patch_edge_array.size();
		int add_edge_idx_2 = (int) m_patch_edge_array.size() + 1;
		int add_patch_idx = (int) m_patch_array.size();
		/// add a new conner
		PatchConner add_conner;
		if(high_val_conner.m_conner_type == 0) add_conner.m_conner_type = 2;
		else add_conner.m_conner_type = 1;
		add_conner.m_mesh_index = mid_vert;
		add_conner.m_nb_conner_index_array.push_back(pre_conner_idx);
		add_conner.m_nb_conner_index_array.push_back(mid_conner_idx);
		add_conner.m_nb_conner_index_array.push_back(nxt_conner_idx);
		add_conner.m_nb_edge_index_array.push_back(add_edge_idx_1);
		add_conner.m_nb_edge_index_array.push_back(rm_edge_idx);
		add_conner.m_nb_edge_index_array.push_back(add_edge_idx_2);
		m_patch_conner_array.push_back(add_conner);

		/// add two edges
		PatchEdge add_edge_1, add_edge_2;
		add_edge_1.m_conner_pair_index = std::make_pair(pre_conner_idx, add_conner_idx);
		add_edge_1.m_nb_patch_index_array.push_back(pre_patch_idx);
		add_edge_1.m_nb_patch_index_array.push_back(add_patch_idx);
		FindShortestPathInPatch(pre_patch_idx, m_patch_conner_array[pre_conner_idx].m_mesh_index,
			mid_vert, add_edge_1.m_mesh_path);

		add_edge_2.m_conner_pair_index = std::make_pair(nxt_conner_idx, add_conner_idx);
		add_edge_2.m_nb_patch_index_array.push_back(nxt_patch_idx);
		add_edge_2.m_nb_patch_index_array.push_back(add_patch_idx);
		FindShortestPathInPatch(nxt_patch_idx, m_patch_conner_array[nxt_conner_idx].m_mesh_index,
			mid_vert, add_edge_2.m_mesh_path);

		m_patch_edge_array.push_back(add_edge_1);
		m_patch_edge_array.push_back(add_edge_2);

		/// add a new patch
		ParamPatch add_patch;
		add_patch.m_conner_index_array.push_back(high_val_conner_idx);
		add_patch.m_edge_index_array.push_back(pre_edge_idx);
		add_patch.m_conner_index_array.push_back(pre_conner_idx);
		add_patch.m_edge_index_array.push_back(add_edge_idx_1);
		add_patch.m_conner_index_array.push_back(add_conner_idx);
		add_patch.m_edge_index_array.push_back(add_edge_idx_2);
		add_patch.m_conner_index_array.push_back(nxt_conner_idx);
		add_patch.m_edge_index_array.push_back(nxt_edge_idx);
		m_patch_array.push_back(add_patch);
		FindPatchInnerFace(add_patch_idx);

		/// modify other conner/edge/patch info
		/// modify high valance conner's neighbor
		std::vector<int>& high_val_conner_nb_conners = m_patch_conner_array[high_val_conner_idx].m_nb_conner_index_array;
		std::vector<int>& high_val_conner_nb_edges = m_patch_conner_array[high_val_conner_idx].m_nb_edge_index_array;
		std::vector<int>::iterator rm_edge_iter = find(high_val_conner_nb_edges.begin(), 
			high_val_conner_nb_edges.end(), rm_edge_idx);
		assert(rm_edge_iter != high_val_conner_nb_edges.end());
		int pos = rm_edge_iter - high_val_conner_nb_edges.begin();
		high_val_conner_nb_edges.erase(rm_edge_iter);
		high_val_conner_nb_conners.erase(high_val_conner_nb_conners.begin() + pos);

		/// modify middle conner's neighbor
		std::vector<int>& mid_conner_nb_conners = m_patch_conner_array[mid_conner_idx].m_nb_conner_index_array;
		std::vector<int>& mid_conner_nb_edges = m_patch_conner_array[mid_conner_idx].m_nb_edge_index_array;
		for(size_t k=0; k<mid_conner_nb_edges.size(); ++k){
			if(mid_conner_nb_edges[k] == rm_edge_idx){ mid_conner_nb_conners[k] = add_conner_idx; break;}
		}

		/// modify previous conner's neighbor
		std::vector<int>& pre_conner_nb_conners = m_patch_conner_array[pre_conner_idx].m_nb_conner_index_array;
		std::vector<int>& pre_conner_nb_edges = m_patch_conner_array[pre_conner_idx].m_nb_edge_index_array;		
		std::vector<int>::iterator pre_edge_iter = find(pre_conner_nb_edges.begin(), pre_conner_nb_edges.end(), pre_edge_idx);
		assert(pre_edge_iter != pre_conner_nb_edges.end());
		pos = pre_edge_iter - pre_conner_nb_edges.begin();
		pre_conner_nb_conners.insert(pre_conner_nb_conners.begin()+pos, add_conner_idx);
		pre_conner_nb_edges.insert(pre_edge_iter, add_edge_idx_1);


		/// modify next conner's neighbor
		std::vector<int>& nxt_conner_nb_conners = m_patch_conner_array[nxt_conner_idx].m_nb_conner_index_array;
		std::vector<int>& nxt_conner_nb_edgs = m_patch_conner_array[nxt_conner_idx].m_nb_edge_index_array;
		std::vector<int>::iterator nxt_edge_iter = find(nxt_conner_nb_edgs.begin(), nxt_conner_nb_edgs.end(), nxt_edge_idx);
		assert(nxt_edge_iter != nxt_conner_nb_edgs.end());
		pos = nxt_edge_iter - nxt_conner_nb_edgs.begin() + 1;
		nxt_conner_nb_conners.insert(nxt_conner_nb_conners.begin() + pos, add_conner_idx);
		nxt_conner_nb_edgs.insert(nxt_edge_iter+1, add_edge_idx_2);

		/// modify pervious/next/rm edge's neighbor
		PatchEdge& prev_edge = m_patch_edge_array[pre_edge_idx];
		PatchEdge& next_edge = m_patch_edge_array[nxt_edge_idx];
		PatchEdge& del_edge = m_patch_edge_array[rm_edge_idx];

		for(size_t k=0; k<prev_edge.m_nb_patch_index_array.size(); ++k){
			if(prev_edge.m_nb_patch_index_array[k] == pre_patch_idx){
				prev_edge.m_nb_patch_index_array[k] = add_patch_idx;
			}
		}
		for(size_t k=0; k<next_edge.m_nb_patch_index_array.size(); ++k){
			if(next_edge.m_nb_patch_index_array[k] == nxt_patch_idx){
				next_edge.m_nb_patch_index_array[k] = add_patch_idx;
			}
		}
		if(del_edge.m_conner_pair_index.first == high_val_conner_idx){
			del_edge.m_conner_pair_index.first = add_conner_idx;
			del_edge.m_mesh_path.erase(del_edge.m_mesh_path.begin(), del_edge.m_mesh_path.begin() + mid_vert_idx);
		}else{
			del_edge.m_conner_pair_index.second = add_conner_idx;
			del_edge.m_mesh_path.erase(del_edge.m_mesh_path.begin() + mid_vert_idx + 1, del_edge.m_mesh_path.end());
		}

		/// modify prev/next patch
		std::vector<int>& pre_patch_edges = m_patch_array[pre_patch_idx].m_edge_index_array;
		std::vector<int>& nxt_patch_edges = m_patch_array[nxt_patch_idx].m_edge_index_array;
		for(size_t k=0; k<pre_patch_edges.size(); ++k){
			if(pre_patch_edges[k] == pre_edge_idx){
				pre_patch_edges[k] = add_edge_idx_1; break;
			}
		}
		for(size_t k=0; k<nxt_patch_edges.size(); ++k){
			if(nxt_patch_edges[k] == nxt_edge_idx){
				nxt_patch_edges[k] = add_edge_idx_2; break;
			}
		}
		SetPatchConners(pre_patch_idx);
		SetPatchConners(nxt_patch_idx);
		FindPatchInnerFace(pre_patch_idx);
		FindPatchInnerFace(nxt_patch_idx);
	}

	void ChartCreator::FindShortestPathInPatch(int patch_id, int start_vid, int end_vid, std::vector<int>& path) const
	{
		const ParamPatch& patch = m_patch_array[patch_id];
		const std::vector<int>& patch_edge_vec = patch.m_edge_index_array;
		const std::vector<int>& face_vec = patch.m_face_index_array;

		std::set< std::pair<int, int> > region_mesh_edge_set;
		
		const PolyIndexArray& face_list_array = p_mesh->m_Kernel.GetFaceInfo().GetIndex();
		for(size_t k=0; k<face_vec.size(); ++k){
			int fid = face_vec[k];
			const IndexArray& face = face_list_array[fid];
			for(int i=0; i<3; ++i){
				region_mesh_edge_set.insert(MakeEdge(face[i], face[(i+1)%3]));
			}
		}		

		if(!FindShortestPathInRegion(p_mesh, start_vid, end_vid, region_mesh_edge_set, path)){					
			
			for(size_t k=0; k<patch_edge_vec.size(); ++k){
				const std::vector<int>& path = m_patch_edge_array[patch_edge_vec[k]].m_mesh_path;
				for(size_t i=1; i<path.size(); ++i){
					region_mesh_edge_set.insert(MakeEdge(path[i-1], path[i]));
				}
			}
			if(!FindShortestPathInRegion(p_mesh, start_vid, end_vid, region_mesh_edge_set, path))
			{
				std::cout << "Can't find the patch !" << std::endl;
			}
		}
	}
	
	void ChartCreator::FindShortestPathInTwoAdjPatch(int patch_id1, int patch_id2, 
		int start_vid, int end_vid, std::vector<int>& path) const
	{
		const ParamPatch& patch1 = m_patch_array[patch_id1];
		const ParamPatch& patch2 = m_patch_array[patch_id2];

		const std::vector<int>& patch_edge_vec1 = patch1.m_edge_index_array;
		const std::vector<int>& patch_edge_vec2 = patch2.m_edge_index_array;
		const std::vector<int>& face_vec1 = patch1.m_face_index_array;
		const std::vector<int>& face_vec2 = patch2.m_face_index_array;

		std::set< std::pair<int, int> > region_mesh_edge_set;
		

		const PolyIndexArray& face_list_array = p_mesh->m_Kernel.GetFaceInfo().GetIndex();
		for(size_t k=0; k<face_vec1.size(); ++k){
			int fid = face_vec1[k];
			const IndexArray& face = face_list_array[fid];
			for(int i=0; i<3; ++i){
				region_mesh_edge_set.insert(MakeEdge(face[i], face[(i+1)%3]));
			}
		}		
		for(size_t k=0; k<face_vec2.size(); ++k){
			int fid = face_vec2[k];
			const IndexArray& face = face_list_array[fid];
			for(int i=0; i<3; ++i){
				region_mesh_edge_set.insert(MakeEdge(face[i], face[(i+1)%3]));
			}
		}

		if(!FindShortestPathInRegion(p_mesh, start_vid, end_vid, region_mesh_edge_set, path))
		{
			for(size_t k=0; k<patch_edge_vec1.size(); ++k){
				const std::vector<int>& path = m_patch_edge_array[patch_edge_vec1[k]].m_mesh_path;
				for(size_t i=1; i<path.size(); ++i){
					region_mesh_edge_set.insert(MakeEdge(path[i-1], path[i]));
				}
			}
			for(size_t k=0; k<patch_edge_vec2.size(); ++k){
				const std::vector<int>& path = m_patch_edge_array[patch_edge_vec2[k]].m_mesh_path;
				for(size_t i=1; i<path.size(); ++i){
					region_mesh_edge_set.insert(MakeEdge(path[i-1], path[i]));
				}
			}
			if(!FindShortestPathInRegion(p_mesh, start_vid, end_vid, region_mesh_edge_set, path))
			{
				std::cout <<" Can't find the path !" << std::endl;
			}
		}


	}

	void ChartCreator::OptimizePatchShape()
	{
		for(size_t k=0; k<m_patch_edge_array.size(); ++k)
		{
			PatchEdge& patch_edge = m_patch_edge_array[k];
			if(patch_edge.m_nb_patch_index_array.size() !=2 ) continue;
			int patch_id1 = patch_edge.m_nb_patch_index_array[0];
			int patch_id2 = patch_edge.m_nb_patch_index_array[1];		
			int start_vid = patch_edge.m_mesh_path[0];
			int end_vid = patch_edge.m_mesh_path[patch_edge.m_mesh_path.size()-1];
			FindShortestPathInTwoAdjPatch(patch_id1, patch_id2, start_vid, end_vid, patch_edge.m_mesh_path);
			FindPatchInnerFace(patch_id1);
			FindPatchInnerFace(patch_id2);
		}
	}

	int ChartCreator::GetCommonConnerBetweenTwoEdge(const PatchEdge& edge_1, const PatchEdge& edge_2)
	{
		int idx1 = edge_1.m_conner_pair_index.first;
		int idx2 = edge_1.m_conner_pair_index.second;
		int idx3 = edge_2.m_conner_pair_index.first;
		int idx4 = edge_2.m_conner_pair_index.second;

		if(idx1 == idx3 || idx1 == idx4) return idx1;
		else if(idx2 == idx3 || idx2 == idx4) return idx2;
		return -1;	 
	}

	std::vector<int> ChartCreator::GetCommonEdgesBetweenTwoPatch(const ParamPatch& patch_1, const ParamPatch& patch_2)
	{
		const std::vector<int>& edge_vec_1 = patch_1.m_edge_index_array;
		const std::vector<int>& edge_vec_2 = patch_2.m_edge_index_array;

		std::vector<int> com_edges;
		for(size_t i=0; i<edge_vec_1.size(); ++i)
		{
			int edge_idx = edge_vec_1[i];
			if(find(edge_vec_2.begin(), edge_vec_2.end(), edge_idx) != edge_vec_2.end())
			{
				com_edges.push_back(edge_idx);
			}
		}
		return com_edges;
	}

	std::vector<int> ChartCreator::GetCommonConnersBetweenTwoPatch(const ParamPatch& patch_1, const ParamPatch& patch_2)
	{
		const std::vector<int>& conner_vec_1 = patch_1.m_conner_index_array;
		const std::vector<int>& conner_vec_2 = patch_2.m_conner_index_array;

		std::vector<int> com_conners;
		for(size_t i=0; i<conner_vec_1.size(); ++i)
		{
			int conner_idx = conner_vec_1[i];
			if(find(conner_vec_2.begin(), conner_vec_2.end(), conner_idx) != conner_vec_2.end())
			{
				com_conners.push_back(conner_idx);
			}
		}
		return com_conners;
	}

	void ChartCreator::GetTransListBetweenTwoChart(int chart_id1, int chart_id2, int edge_idx1, int edge_idx2)
	{
		const ParamPatch& patch_1 = m_patch_array[chart_id1];
		const ParamPatch& patch_2 = m_patch_array[chart_id2];

		const PatchEdge& edge_1 = m_patch_edge_array[edge_idx1];
		const PatchEdge& edge_2 = m_patch_edge_array[edge_idx2];

		int com_conner_idx(-1);
		if(edge_1.m_conner_pair_index.first == edge_2.m_conner_pair_index.first || 
			edge_1.m_conner_pair_index.first == edge_2.m_conner_pair_index.second){
				com_conner_idx = edge_1.m_conner_pair_index.first;
		}else if(edge_1.m_conner_pair_index.second == edge_2.m_conner_pair_index.first ||
			edge_1.m_conner_pair_index.second == edge_2.m_conner_pair_index.second){
				com_conner_idx = edge_1.m_conner_pair_index.second;
		}else{
			std::cerr <<"Error, these two edges havn't common conner." << std::endl;
			return;
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
