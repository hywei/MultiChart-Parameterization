#ifndef PARAMDRAWER_H_
#define PARAMDRAWER_H_

#include "../Common/BasicDataType.h"
#include "Barycentric.h"

namespace PARAM
{
	class Parameter;
	class SurfaceCoord;

	class ParamDrawer
	{
	public:
		enum DrawMode
		{
			DRAWNOTHING = 0x00000000,

			//! draw patch info
			DRAWLAYOUT = 0x00000001,
			DRAWPATCHCONNER = 0x00000011, 
			DRAWPATCHEDGE = 0x000000021,
			DRAWPATCHFACE = 0x000000041,

			//! draw distortion info
			DRAWDISTORTION = 0x00000002,
			DRAWFACEHARMONICDISTORTION = 0x00000012,
			DRAWVERTEXHARMONICDISTORTION = 0x00000022,
			
			//! Draw texture
			DRAWFACETEXTURE = 0x00000004,

			//! Draw Corresponding
			DRAWCORRESPONGING = 0x00000008,
			DRAWSELECTION = 0x000000018
		};
	public:
		void SetUnCorrespondingVertArray(const std::vector<int>& uncorresponding_vert_array)
		{
			m_uncorrespondnig_vert_array = uncorresponding_vert_array;
		}

		void SetSelectedVertCoord(const Coord& select_coord);
		void SetSelectedVertCoord(const SurfaceCoord& select_coord);
			
		void SetSelectedPatchCoord(const Coord& select_coord);

		int FindSelectedVertId(const Coord& select_coord);	
		SurfaceCoord FindSelectedSurfaceCoord(const Coord& select_coord);
		int GetSelectedVertID() const { return m_selected_vert_id; }
		SurfaceCoord GetSelectedSurfaceCorod() const { return m_selected_surface_coord; }

	public:
		void SetDrawMode(DrawMode mode){ m_draw_mode = mode;}
		void SetDrawPatchConner(bool is_draw) { m_draw_patch_conner = is_draw; }
		void SetDrawPatchEdge(bool is_draw) { m_draw_patch_edge = is_draw; }
		void SetDrawPatchFace(bool is_draw) { m_draw_patch_face = is_draw; }
		void SetDrawOutRangeVertices(bool is_draw) { m_draw_out_range_vertices = is_draw; }
		void SetDrawSelectedPatch(bool is_draw) { m_draw_select_patch = is_draw; }
		void SetDrawFlipFace(bool is_draw) { m_draw_flip_face = is_draw; }
		void SetDrawSelectedVertex(bool is_draw) { m_draw_selected_vertex = is_draw; }
		void SetDrawUnCorresponding(bool is_draw) { m_draw_uncorrespoinding = is_draw; }

	public:
		ParamDrawer(const Parameter& quad_param);
		~ParamDrawer();

		void Draw() const;

	private:
		void DrawPatchConner() const;
		void DrawPatchEdge() const;
		void DrawPatchFace() const;

		void DrawBaseDomain() const;

		void DrawOutRangeVertex() const;
		void DrawUnSetFace() const;
		void DrawFlipedFace() const;

		void DrawFaceDistortion() const;
		void DrawFaceTexture() const;

		void DrawCorresponding() const;
		void DrawUnCorrespondingVertex() const;

		void DrawSelectedPatch() const;

		void DrawSelectedVert() const;

		void DrawSphere(const Coord& center, double point_size = 1.0) const;

	private:
		const Parameter& m_parameter;

		std::vector<int> m_uncorrespondnig_vert_array;

		int m_selected_vert_id;
		Coord m_selected_vert_coord;
		SurfaceCoord m_selected_surface_coord;

		DrawMode m_draw_mode;

		bool m_draw_patch_conner;
		bool m_draw_patch_edge;
		bool m_draw_patch_face;
		bool m_draw_out_range_vertices;
		bool m_draw_flip_face;
		bool m_draw_selected_vertex;
		bool m_draw_select_patch;
		bool m_draw_uncorrespoinding;

		int m_selected_patch_id;
	};
}

#endif // PARAMDRAWER_H_
