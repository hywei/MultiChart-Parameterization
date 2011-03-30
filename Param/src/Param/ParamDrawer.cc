#include "ParamDrawer.h"
#include "ParamPatch.h"
#include "Parameter.h"
#include "ChartCreator.h"
#include "../ModelMesh/MeshModel.h"
#include <vector>
#include <limits>

#ifdef WIN32
#include <GL/GLAux.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

using namespace std;
namespace PARAM
{
	ParamDrawer::ParamDrawer(const Parameter& _param) 
		: m_parameter(_param) {
			m_selected_patch_id = -1;
			m_selected_vert_coord = -1; 
			m_draw_patch_conner = m_draw_patch_edge = m_draw_patch_face = false;
			m_draw_out_range_vertices = false;
			m_draw_select_patch = false;	
			m_draw_selected_vertex = false;
			m_draw_flip_face = false;
			m_draw_uncorrespoinding = false;
			m_selected_surface_coord = SurfaceCoord(-1, 0, 0, 0);
	}
	ParamDrawer::~ParamDrawer(){}

	void ParamDrawer::Draw() const
	{
		glEnable(GL_LIGHTING);
		glEnable(GL_POLYGON_SMOOTH);
		glPolygonMode(GL_FRONT, GL_FILL);
		glEnable(GL_COLOR_MATERIAL);
		glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);

		//if(m_draw_mode & DRAWPATCHCONNER) DrawPatchConner();
		if(m_draw_patch_conner) DrawPatchConner();
		if(m_draw_patch_face) DrawPatchFace();
		if(m_draw_out_range_vertices) DrawOutRangeVertex();
		if(m_draw_select_patch) DrawSelectedPatch();
		if(m_draw_flip_face) DrawFlipedFace();
		if(m_draw_selected_vertex) DrawSelectedVert();
        // std::cout << "Select Vertex : " << m_selected_vert_id << " : "<<
        //     m_selected_vert_coord[0] << " " << m_selected_vert_coord[1] << ' '
        //           << m_selected_vert_coord[2] << std::endl;
		
		if(m_draw_uncorrespoinding) DrawUnCorrespondingVertex();	   

		// DrawUnSetFace();

		
        
		glDisable(GL_COLOR_MATERIAL);
		glDisable(GL_POLYGON_SMOOTH);

		glDisable(GL_LIGHTING);
		glDepthFunc(GL_LEQUAL);


		glEnable(GL_LINE_SMOOTH);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
		glLineWidth(3.0f);
	
		glPolygonOffset(1.0, 1.0);
		glEnable(GL_POLYGON_OFFSET_FILL);

		//if(m_draw_mode & DRAWPATCHEDGE) DrawPatchEdge();
		if(m_draw_patch_edge) DrawPatchEdge();

		//DrawBaseDomain();

		glDisable(GL_POLYGON_OFFSET_FILL);

		glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
		glDisable(GL_BLEND);
		glDisable(GL_LINE_SMOOTH);
		glDisable(GL_POLYGON_OFFSET_FILL);
		glEnable(GL_LIGHTING);
		glDepthFunc(GL_LESS);
		glEnable(GL_LIGHTING);
	
	
	}

	void ParamDrawer::DrawPatchConner() const
	{
		boost::shared_ptr<MeshModel> p_mesh = m_parameter.GetMeshModel();
		boost::shared_ptr<ChartCreator> p_chart_creator = m_parameter.GetChartCreator();

		if(p_mesh == NULL) return ;
		if(p_chart_creator == NULL) return;
		const CoordArray& vCoord = p_mesh->m_Kernel.GetVertexInfo().GetCoord();

		int colors[56][3] = 
		{
			{255, 0, 0}, {0, 255, 0}, {0, 0, 255}, {255, 255, 0}, 
			{ 255, 0 , 255}, {0, 255, 255}, {255, 255, 255}, {0, 0, 0},

			{255, 128, 128}, {0, 64, 128},  {255, 128, 192}, {128, 255, 128}, 
			{0, 255, 128}, {128, 255, 255}, {0, 128, 255}, {255, 128, 255}, 

			{255, 0, 0}, {255, 255, 0}, {128, 255, 0}, {0, 255, 64}, 
			{0, 255, 255}, {0, 128, 192}, {128, 128, 192}, {255, 0, 255}, 

			{128, 64, 64}, {255, 128, 64}, {0, 255, 0}, {0, 128, 128}, 
			{0, 64, 128}, {128, 128, 255}, {128, 0, 64}, {255, 0, 128}, 

			{128, 0, 0}, {0, 128, 0}, {0, 128, 64}, {255, 255, 128},
			{0, 0, 255}, {0, 0, 160}, {128, 0, 128}, {128, 0, 255}, 

			{64, 0, 0}, {128, 64, 0}, {0, 64, 0}, {0, 64, 64}, 
			{0, 0, 128}, {0, 0, 64}, {64, 0, 64}, {64, 0, 128}, 
		};


		const std::vector<PatchConner>& conner_array = p_chart_creator->GetPatchConnerArray();
		for(size_t k=0; k<conner_array.size(); ++k)
		{            
			int c = k%56;
			glColor3ub(colors[c][0], colors[c][1], colors[c][2]);
			const PatchConner& conner = conner_array[k];
			if(conner.m_conner_type == 2){
				glColor3ub(255, 0, 0);
			}else if(conner.m_conner_type == 1){
				glColor3ub(0, 255, 0);
			}else if(conner.m_conner_type == 0){
				glColor3ub(0, 0, 255);
			}
            glColor3ub(colors[c][0], colors[c][1], colors[c][2]);
			int mesh_idx = conner_array[k].m_mesh_index;
			DrawSphere(vCoord[mesh_idx]);
		}

	}

	void ParamDrawer::DrawPatchEdge() const
	{
		boost::shared_ptr<MeshModel> p_mesh = m_parameter.GetMeshModel();
		boost::shared_ptr<ChartCreator> p_chart_creator = m_parameter.GetChartCreator();

		if(p_mesh == NULL) return;
		if(p_chart_creator == NULL) return;
		const CoordArray& vCoord = p_mesh->m_Kernel.GetVertexInfo().GetCoord();
		const std::vector<PatchEdge>& patch_edge_array = p_chart_creator->GetPatchEdgeArray();

		int colors[56][3] = 
		{			  
			{255, 0, 0}, {0, 255, 0}, {0, 0, 255}, {255, 255, 0}, 
			{ 255, 0 , 255}, {0, 255, 255}, {255, 255, 255}, {0, 0, 0},

			{255, 128, 128}, {0, 64, 128},  {255, 128, 192}, {128, 255, 128}, 
			{0, 255, 128}, {128, 255, 255}, {0, 128, 255}, {255, 128, 255}, 

			{255, 0, 0}, {255, 255, 0}, {128, 255, 0}, {0, 255, 64}, 
			{0, 255, 255}, {0, 128, 192}, {128, 128, 192}, {255, 0, 255}, 

			{128, 64, 64}, {255, 128, 64}, {0, 255, 0}, {0, 128, 128}, 
			{0, 64, 128}, {128, 128, 255}, {128, 0, 64}, {255, 0, 128}, 

			{128, 0, 0}, {0, 128, 0}, {0, 128, 64}, {255, 255, 128},
			{0, 0, 255}, {0, 0, 160}, {128, 0, 128}, {128, 0, 255}, 

			{64, 0, 0}, {128, 64, 0}, {0, 64, 0}, {0, 64, 64}, 
			{0, 0, 128}, {0, 0, 64}, {64, 0, 64}, {64, 0, 128}, 
		};

		for(size_t k=0; k<patch_edge_array.size(); ++k)
		{
			int c = k%56;
			glColor3ub(colors[c][0], colors[c][1], colors[c][2]);
			const std::vector<int>& path = patch_edge_array[k].m_mesh_path;
			glBegin(GL_LINE_STRIP);
			for(size_t i=0; i<path.size(); ++i)
			{
				const Coord& vtxCoord = vCoord[path[i]];
				glVertex3d(vtxCoord[0], vtxCoord[1], vtxCoord[2]);
			}
			glEnd();
		}
	}

    void ParamDrawer::DrawPatchFace() const
    {
        boost::shared_ptr<MeshModel> p_mesh = m_parameter.GetMeshModel();
        boost::shared_ptr<ChartCreator> p_chart_creator = m_parameter.GetChartCreator();
        
        if(p_mesh == NULL) return;
        if(p_chart_creator == NULL) return;

        int colors[56][3] = 
		{			  
			{255, 0, 0}, {0, 255, 0}, {0, 0, 255}, {255, 255, 0}, 
			{ 255, 0 , 255}, {0, 255, 255}, {255, 255, 255}, {0, 0, 0},

			{255, 128, 128}, {0, 64, 128},  {255, 128, 192}, {128, 255, 128}, 
			{0, 255, 128}, {128, 255, 255}, {0, 128, 255}, {255, 128, 255}, 

			{255, 0, 0}, {255, 255, 0}, {128, 255, 0}, {0, 255, 64}, 
			{0, 255, 255}, {0, 128, 192}, {128, 128, 192}, {255, 0, 255}, 

			{128, 64, 64}, {255, 128, 64}, {0, 255, 0}, {0, 128, 128}, 
			{0, 64, 128}, {128, 128, 255}, {128, 0, 64}, {255, 0, 128}, 

			{128, 0, 0}, {0, 128, 0}, {0, 128, 64}, {255, 255, 128},
			{0, 0, 255}, {0, 0, 160}, {128, 0, 128}, {128, 0, 255}, 

			{64, 0, 0}, {128, 64, 0}, {0, 64, 0}, {0, 64, 64}, 
			{0, 0, 128}, {0, 0, 64}, {64, 0, 64}, {64, 0, 128}, 
		};

        
        const std::vector<ParamPatch>& patch_array = p_chart_creator->GetPatchArray();
        const PolyIndexArray& face_list_array = p_mesh->m_Kernel.GetFaceInfo().GetIndex();
        const NormalArray& vert_norm_array = p_mesh->m_Kernel.GetVertexInfo().GetNormal();
        const CoordArray& vert_coord_array = p_mesh->m_Kernel.GetVertexInfo().GetCoord();
		const NormalArray& face_norm_array = p_mesh->m_Kernel.GetFaceInfo().GetNormal();

        glBegin(GL_TRIANGLES);
        for(size_t k=0; k<patch_array.size(); ++k){           
            const ParamPatch& patch = patch_array[k];
            const std::vector<int>& faces = patch.m_face_index_array;
            glColor3ub(colors[k%56][0], colors[k%56][1], colors[k%56][2]);
            for(size_t i=0; i<faces.size(); ++i){
                int fid = faces[i];
                const IndexArray& vertices = face_list_array[fid];
                for(int j=0; j<3; ++j){
                    const Coord& v = vert_coord_array[vertices[j]];
                    const Normal& n = vert_norm_array[vertices[j]];
                    glNormal3d(n[0], n[1], n[2]);
                    glVertex3d(v[0], v[1], v[2]);
                }
            }
        }
        glEnd();
                    
    }

	void ParamDrawer::DrawOutRangeVertex() const
	{
		boost::shared_ptr<MeshModel> p_mesh = m_parameter.GetMeshModel();		
		if(p_mesh == NULL) return ;
		const std::vector<int>& out_range_vert_array = m_parameter.GetOutRangeVertArray();
		const CoordArray& vCoord = p_mesh->m_Kernel.GetVertexInfo().GetCoord();

		glColor3ub(255, 255, 0);
		for(size_t k=0; k<out_range_vert_array.size(); ++k)
		{
			int vid = out_range_vert_array[k];
			const Coord& vtxCoord = vCoord[vid];
			DrawSphere(vtxCoord, 0.5);
		}
	}

	void ParamDrawer::DrawSelectedPatch() const
	{
		//std::cout << "Selected Patch id " << m_selected_patch_id << std::endl;
		if(m_selected_patch_id == -1) return;
		boost::shared_ptr<ChartCreator> p_chart_creator = m_parameter.GetChartCreator();
		if(p_chart_creator == NULL) return ;
		boost::shared_ptr<MeshModel> p_mesh = m_parameter.GetMeshModel();	
		if(p_mesh == NULL) return;

		int colors[8][3] = 
		{			  
			{255, 0, 0}, {0, 255, 0}, {0, 0, 255}, {255, 255, 0}, 
			{ 255, 0 , 255}, {0, 255, 255}, {255, 255, 255}, {0, 0, 0}
		};

		const CoordArray& vCoord = p_mesh->m_Kernel.GetVertexInfo().GetCoord();
		const PolyIndexArray& face_list_array = p_mesh->m_Kernel.GetFaceInfo().GetIndex();

		const std::vector<ParamPatch>& patch_array = p_chart_creator->GetPatchArray();
		const std::vector<PatchConner>& patch_conners = p_chart_creator->GetPatchConnerArray();
		const std::vector<PatchEdge>& patch_edges = p_chart_creator->GetPatchEdgeArray();

		const ParamPatch& patch = patch_array[m_selected_patch_id];

		for(size_t k=0; k<patch.m_conner_index_array.size(); ++k)
		{
			int conner_idx = patch.m_conner_index_array[k];
			int vid = patch_conners[conner_idx].m_mesh_index;
			int c = k%8;
			glColor3ub(colors[c][0], colors[c][1], colors[c][2]);
			DrawSphere(vCoord[vid]);
		}

		glDisable(GL_LIGHTING);
		for(size_t k=0; k<patch.m_edge_index_array.size(); ++k)
		{
			int edge_idx = patch.m_edge_index_array[k];
			const std::vector<int>& path = patch_edges[edge_idx].m_mesh_path;
			int c= k%8;
			glColor3ub(colors[c][0], colors[c][1], colors[c][2]);
			glBegin(GL_LINE_STRIP);
			for(size_t i=0; i<path.size(); ++i)
			{
				const Coord& vtxCoord = vCoord[path[i]];
				glVertex3d(vtxCoord[0], vtxCoord[1], vtxCoord[2]);
			}
			glEnd();

		}
		

		glColor3ub(0, 128, 0);
		for(size_t k=0; k<patch.m_face_index_array.size(); ++k)
		{
			int fid = patch.m_face_index_array[k];
			const IndexArray& fIndex = face_list_array[fid];
			glBegin(GL_TRIANGLES);

			for(int j = 0; j < 3; ++ j)
			{
				int vID = fIndex[j];
				const Coord& v = vCoord[vID];
				glVertex3d(v[0], v[1], v[2]);
			}

			glEnd();
		}

		const std::vector<int>& patch_nb_vec = patch.m_nb_patch_index_array;
		for(size_t i=0; i<patch_nb_vec.size(); ++i)
		{
			int c = i%8;
			glColor3ub(colors[c][0], colors[c][1], colors[c][2]);
			const ParamPatch& nb_patch = patch_array[patch_nb_vec[i]];
			for(size_t k=0; k<nb_patch.m_face_index_array.size(); ++k)
			{
				int fid = nb_patch.m_face_index_array[k];
				const IndexArray& fIndex = face_list_array[fid];
				glBegin(GL_TRIANGLES);

				for(int j = 0; j < 3; ++ j)
				{
					int vID = fIndex[j];
					const Coord& v = vCoord[vID];
					glVertex3d(v[0], v[1], v[2]);
				}

				glEnd();
			}
		}

		glEnable(GL_LIGHTING);
	}

	void ParamDrawer::DrawUnCorrespondingVertex() const
	{
		boost::shared_ptr<MeshModel> p_mesh = m_parameter.GetMeshModel();
		if(p_mesh == NULL) return;

		const CoordArray& vCoord = p_mesh->m_Kernel.GetVertexInfo().GetCoord();
		glColor3ub(128, 128, 0);
		for(size_t k=0; k<m_uncorrespondnig_vert_array.size(); ++k)
		{
			int vid = m_uncorrespondnig_vert_array[k];
			const Coord& vtxCoord = vCoord[vid];
			DrawSphere(vtxCoord, 1.0);
		}
	}

	void ParamDrawer::DrawUnSetFace() const
	{
		boost::shared_ptr<MeshModel> p_mesh = m_parameter.GetMeshModel();		
		if(p_mesh == NULL) return ;
		const std::vector<int>& unset_face_array = m_parameter.GetUnSetFaceArray();
		const PolyIndexArray& face_list_array = p_mesh->m_Kernel.GetFaceInfo().GetIndex();
		const CoordArray& vCoord = p_mesh->m_Kernel.GetVertexInfo().GetCoord();

		glColor3ub(0, 0, 255);
		for(size_t k=0; k<unset_face_array.size(); ++k)
		{
			int fid = unset_face_array[k];
			
			const IndexArray& fIndex = face_list_array[fid];
			glBegin(GL_TRIANGLES);

			for(int j = 0; j < 3; ++ j)
			{
					int vID = fIndex[j];
					const Coord& v = vCoord[vID];
					glVertex3d(v[0], v[1], v[2]);
			}
			
			glEnd();
		}

	}

	void ParamDrawer::DrawFlipedFace() const
	{
		boost::shared_ptr<MeshModel> p_mesh = m_parameter.GetMeshModel();		
		if(p_mesh == NULL) return ;
		const std::vector<int>& fliped_face_array = m_parameter.GetFlipedFaceArray();
		const PolyIndexArray& face_list_array = p_mesh->m_Kernel.GetFaceInfo().GetIndex();
		const CoordArray& vCoord = p_mesh->m_Kernel.GetVertexInfo().GetCoord();

		glDisable(GL_LIGHTING);

		glColor3ub(255, 255, 0);
		for(size_t k=0; k<fliped_face_array.size(); ++k)
		{
			int fid = fliped_face_array[k];

			const IndexArray& fIndex = face_list_array[fid];
			glBegin(GL_TRIANGLES);

			for(int j = 0; j < 3; ++ j)
			{
				int vID = fIndex[j];
				const Coord& v = vCoord[vID];
				glVertex3d(v[0], v[1], v[2]);
			}

			glEnd();
		}

		glEnable(GL_LIGHTING);
	}

	void ParamDrawer::DrawSelectedVert() const
	{
		if(m_selected_vert_id == -1) return;
		glDepthMask( GL_FALSE );
		//glColor3ub(64, 0, 128);
		glColor3ub(255, 0, 0);
		DrawSphere(m_selected_vert_coord, 2.2);
		glDepthMask( GL_TRUE);

	}

	int ParamDrawer::FindSelectedVertId(const Coord& select_coord)
	{
		boost::shared_ptr<MeshModel> p_mesh = m_parameter.GetMeshModel();
		if(p_mesh == NULL) return -1;
		const CoordArray& vCoord = p_mesh->m_Kernel.GetVertexInfo().GetCoord();
		int vert_num = p_mesh->m_Kernel.GetModelInfo().GetVertexNum();

		double min_dist = numeric_limits<double>::infinity();
		
		int nearest_vert_id(-1);

		for(int k=0; k<vert_num; ++k)
		{
			const Coord& vtx_coord = vCoord[k];
			double cur_dist = (vtx_coord - select_coord).abs();
			if(min_dist > cur_dist)
			{
				min_dist = cur_dist;
				nearest_vert_id = k;
			}
		}

		return nearest_vert_id;
	}

	SurfaceCoord ParamDrawer::FindSelectedSurfaceCoord(const Coord& select_coord)
	{
		boost::shared_ptr<MeshModel> p_mesh = m_parameter.GetMeshModel();
		if(p_mesh == NULL) return SurfaceCoord(-1, 0, 0, 0);
		const CoordArray& vCoord = p_mesh->m_Kernel.GetVertexInfo().GetCoord();
		const PolyIndexArray& fIndex = p_mesh->m_Kernel.GetFaceInfo().GetIndex();
		const PolyIndexArray& vAdjFaces = p_mesh->m_Kernel.GetVertexInfo().GetAdjFaces();

		if(m_selected_vert_id == -1) return SurfaceCoord(-1, Barycentrc(0, 0, 0));
		const IndexArray& adj_faces = vAdjFaces[m_selected_vert_id];
		

		double min_baryc_corod_error = std::numeric_limits<double>::infinity();
		int min_error_fid;
		Barycentrc min_error_baryc;

		for(int k = 0; k<adj_faces.size(); ++k)
		{
			int fid = adj_faces[k];
			std::vector<Coord> vert_coord_vec(3);
			const IndexArray& face = fIndex[fid];
			for(int i=0; i<3; ++i)
			{
				vert_coord_vec[i] = vCoord[face[i]];
			}
			Barycentrc bc = ComputeVertexBarycentric(vert_coord_vec, select_coord);
			if(IsValidBarycentic(bc)){
				return SurfaceCoord(fid, bc);
			}else{
				double error = ComputeErrorOutValidBarycentric(bc);
				if(error < min_baryc_corod_error){
					min_baryc_corod_error = error;
					min_error_fid = fid;
					min_error_baryc = bc;
				}
			}
		}

		return SurfaceCoord(min_error_fid, min_error_baryc);
	}

	void ParamDrawer::SetSelectedVertCoord(const Coord& select_coord)
	{
		m_selected_vert_coord = select_coord;		
		m_selected_vert_id = FindSelectedVertId(m_selected_vert_coord);
		m_selected_surface_coord = FindSelectedSurfaceCoord(m_selected_vert_coord);
	}

	void ParamDrawer::SetSelectedVertCoord(const SurfaceCoord& select_coord)
	{
		boost::shared_ptr<MeshModel> p_mesh = m_parameter.GetMeshModel();
		if(p_mesh == NULL) return;

		m_selected_surface_coord = select_coord;

		const CoordArray& vCoord = p_mesh->m_Kernel.GetVertexInfo().GetCoord();
		const PolyIndexArray& face_list_array = p_mesh->m_Kernel.GetFaceInfo().GetIndex();

		const IndexArray& face_vert_array = face_list_array[select_coord.face_index];
		Coord pos_coord(0, 0, 0);
		for(int k=0; k<3; ++k)
		{
			int vid = face_vert_array[k];
			pos_coord += vCoord[vid] * (select_coord.barycentric[k]);
		}
		m_selected_vert_coord = pos_coord;

		m_selected_vert_id = FindSelectedVertId(m_selected_vert_coord);
	}

	void ParamDrawer::SetSelectedPatchCoord(const Coord& select_coord)
	{
		int select_vert_id = FindSelectedVertId(select_coord);
		m_selected_patch_id = (m_parameter.GetVertexPatchArray())[select_vert_id];

		const std::vector<ParamPatch>& patch_array = m_parameter.GetChartCreator()->GetPatchArray();
		const std::vector<int>& conner_array = patch_array[m_selected_patch_id].m_conner_index_array;
		const std::vector<int>& edge_array = patch_array[m_selected_patch_id].m_edge_index_array;
		const std::vector<int>& nb_pach_array = patch_array[m_selected_patch_id].m_nb_patch_index_array;

		std::cout << "Select Patch id : " << m_selected_patch_id << std::endl;
		std::cout << "Patch Conners : " ;
		for(size_t i=0; i<conner_array.size(); ++i){
			std::cout << conner_array[i] <<" ";
		}
		std::cout << std::endl;
		std::cout << "Patch Edges : ";
		for(size_t i=0; i<edge_array.size(); ++i){
			std::cout << edge_array[i] <<" ";
		}
		std::cout << std::endl;
		std::cout << "Neighbor Patchs : ";
		for(size_t i=0; i<nb_pach_array.size(); ++i){
			std::cout << nb_pach_array[i] <<" ";
		}
		std::cout << std::endl;
		
	}

	void ParamDrawer::DrawBaseDomain() const
	{
		boost::shared_ptr<MeshModel> p_mesh = m_parameter.GetMeshModel();		
		if(p_mesh == NULL) return ;
		boost::shared_ptr<ChartCreator> p_chart_creator = m_parameter.GetChartCreator();
		
		const CoordArray& vCoord = p_mesh->m_Kernel.GetVertexInfo().GetCoord();
		const std::vector<ParamPatch>& patch_vec = p_chart_creator->GetPatchArray();
		const std::vector<PatchConner>& conner_vec = p_chart_creator->GetPatchConnerArray();

		int colors[56][3] = 
		{			  
			{255, 0, 0}, {0, 255, 0}, {0, 0, 255}, {255, 255, 0}, 
			{ 255, 0 , 255}, {0, 255, 255}, {255, 255, 255}, {0, 0, 0},

			{255, 128, 128}, {0, 64, 128},  {255, 128, 192}, {128, 255, 128}, 
			{0, 255, 128}, {128, 255, 255}, {0, 128, 255}, {255, 128, 255}, 

			{255, 0, 0}, {255, 255, 0}, {128, 255, 0}, {0, 255, 64}, 
			{0, 255, 255}, {0, 128, 192}, {128, 128, 192}, {255, 0, 255}, 

			{128, 64, 64}, {255, 128, 64}, {0, 255, 0}, {0, 128, 128}, 
			{0, 64, 128}, {128, 128, 255}, {128, 0, 64}, {255, 0, 128}, 

			{128, 0, 0}, {0, 128, 0}, {0, 128, 64}, {255, 255, 128},
			{0, 0, 255}, {0, 0, 160}, {128, 0, 128}, {128, 0, 255}, 

			{64, 0, 0}, {128, 64, 0}, {0, 64, 0}, {0, 64, 64}, 
			{0, 0, 128}, {0, 0, 64}, {64, 0, 64}, {64, 0, 128}, 
		};
		for(size_t k=0; k<patch_vec.size(); ++k)
		{
			const std::vector<int>& conner_idx_vec = patch_vec[k].m_conner_index_array;
			glColor3ub(colors[k%56][0], colors[k%56][1], colors[k%56][2]);	

			glBegin(GL_TRIANGLES);

			const Coord& v1 = conner_vec[conner_idx_vec[0]].m_mesh_index;
			const Coord& v2 = conner_vec[conner_idx_vec[1]].m_mesh_index;

			Normal n = cross(v1, v2);
			n.unit();

			for(int j = 0; j < 3; ++ j)
			{
				int vID = conner_vec[conner_idx_vec[j]].m_mesh_index;
				const Coord& v = vCoord[vID];
				glNormal3d(n[0], n[1], n[2]);
				glVertex3d(v[0], v[1], v[2]);
			}

			glEnd();
		}
	}


	void ParamDrawer::DrawSphere(const Coord& center, double point_size) const
	{
		boost::shared_ptr<MeshModel> p_mesh = m_parameter.GetMeshModel();

		double bRadius;
		Coord c;
		p_mesh->m_Kernel.GetModelInfo().GetBoundingSphere(c, bRadius);
		GLUquadricObj* qobj = gluNewQuadric();
		gluQuadricDrawStyle(qobj, GLU_FILL);
		gluQuadricNormals(qobj, GLU_SMOOTH);
		glTranslatef((float)center[0], (float)center[1], (float)center[2]);

		gluSphere(qobj, point_size *0.02 * bRadius , 15, 10);
		glTranslatef((float)-center[0], (float)-center[1], (float)-center[2]);
	}
}
