#ifndef PARAMETERIZATION_H_
#define PARAMETERIZATION_H_

#include "../Common/BasicDataType.h"
#include <boost/shared_ptr.hpp>
#include <set>
#include <map>

class MeshModel;
class CMeshSparseMatrix;

namespace PARAM
{
	class ParamCoord 
	{
	public:
		ParamCoord(){}
		ParamCoord(double s, double t): s_coord(s), t_coord(t){}
		~ParamCoord(){}
        

        bool operator == (const ParamCoord& rhs);
        
		double s_coord;
		double t_coord;
	};

	class ChartParamCoord 
	{
	public:
		ChartParamCoord(double u=0, double v=0, int _chart_id = -1) : param_coord(u, v),
			chart_id(_chart_id){}
		ChartParamCoord(const ParamCoord& _param_coord, int _chart_id = -1) : 
		chart_id(_chart_id), param_coord(_param_coord) {}
		~ChartParamCoord(){}

	public:
		int chart_id;		
		ParamCoord param_coord;
	};

	class SurfaceCoord
	{
	public:
		SurfaceCoord() : face_index(-1), barycentric(0, 0, 0){}
		SurfaceCoord(int fid, const Coord& _coord) : face_index(fid), barycentric(_coord){}
		SurfaceCoord(int fid, double b0, double b1, double b2) : face_index(fid), barycentric(b0, b1, b2) {}
		~SurfaceCoord(){}

		int face_index;
		Coord barycentric;
	};

	void SetCotCoef(const boost::shared_ptr<MeshModel> p_mesh, std::vector<Coord>& cot_coef_vec);
	void SetLapMatrixCoef(const boost::shared_ptr<MeshModel> p_mesh, CMeshSparseMatrix& lap_mat);

	void SetTanCoef(const boost::shared_ptr<MeshModel> p_mesh, std::vector<Coord>& tan_coef_vec);
	void SetLapMatrixCoefWithMeanValueCoord(const boost::shared_ptr<MeshModel> p_mesh, CMeshSparseMatrix& lap_mat);

    int InsidePolygen(const Coord2D& p, const std::vector<Coord2D>& polygen);

	/// functions for computing the distance of a point to a segment
	double xmult(double x1,double y1,double x2,double y2,double x0,double y0);		
	double area_triangle(double x1,double y1,double x2,double y2,double x3,double y3);
	double dis_ptoline(double x1,double y1,double x2,double y2,double ex,double ey);

	double DistanceToTriangle( const Coord2D &p, const Coord2D &a, const Coord2D &b, const Coord2D &c );
	
	
	void FaceValue2VtxColor(boost::shared_ptr<MeshModel> _mesh, std::vector<double>& face_value);


	inline std::pair<int, int> MakeEdge(int vid1, int vid2){ 
		return (vid1 < vid2) ? std::make_pair(vid1, vid2) : std::make_pair(vid2, vid1); }

	bool FindShortestPathInRegion(boost::shared_ptr<MeshModel> p_mesh, int start_vid, int end_vid, 
		const std::set< std::pair<int, int> >& region_edge_set, std::vector<int>& path);

	/// get nearest vertex on a path from another vertex
	double GetNearestVertexOnPath(boost::shared_ptr<MeshModel> p_mesh, int from_vert, 
		const std::vector<int>& path, int& nearest_vid);
	std::vector<int> GetMeshEdgeAdjFaces(boost::shared_ptr<MeshModel> p_mesh, int vid1, int vid2);

	class _HE_edge
	{
	public:
		_HE_edge() : pair(NULL), next(NULL) {}
		_HE_edge(size_t _fid,  size_t _vid) : fid(_fid), vid(_vid), 
			pair(NULL), next(NULL) {} 
	public:
		_HE_edge* pair;
		_HE_edge* next;

		size_t fid;
		size_t vid;
	};

	typedef _HE_edge HE_edge;

	class HalfEdge
	{
	private:
		std::map< std::pair<size_t, size_t>, size_t> edge_map;
		std::vector< _HE_edge* > hf_edge;
	public:
		HalfEdge();
		~HalfEdge();

		int CreateHalfEdge(boost::shared_ptr<MeshModel> mesh);

		_HE_edge* GetHalfEdge(std::pair<size_t, size_t> e);
		const _HE_edge* GetHalfEdge(std::pair<size_t, size_t> e) const;
	};


	int FindInnerFace(boost::shared_ptr<MeshModel> p_mesh, const std::vector<int>& boundary_path,  
		std::vector<int>& face_set, const HalfEdge& half_edge);
}

#endif // PARAMTERIZATION_H_
