#ifndef BARYCENTRIC_H_
#define BARYCENTRIC_H_

#include "../Common/BasicDataType.h"
#include "Parameterization.h"
/*
 * This header is all about barycentric in parameterization *
 * Each vertex in a triangle can be expressed by this triangle's 
 * three nodes and this vertex's barycentric.
 * A valid barycenteric coordinate is in range [0.0 - 1.0]
 */
namespace PARAM
{
	typedef Coord Barycentrc;

	//! compute a triangle's local 2d coord array
	std::vector<Coord2D> ComputeTriangleLocal2DCoord(const std::vector<Coord>& tri_vtx_coord);	
	//! compute a vertex's barycentric in a triangle by each nodes' parameter coordinate
	Barycentrc ComputeVertexBarycentric(const std::vector<ParamCoord>& node_param_coord_array, 
		const ParamCoord& vtx_param_coord);
	Barycentrc ComputeVertexBarycentric(const std::vector<Coord>& face_vert_coord_array, 
		const Coord& vert_coord);
	//! compute a vertex's coordinate with a triangle's nodes coordinate and its barycentric
	Coord ComputeVertexCoord(const std::vector<Coord>& node_coord_array, const Barycentrc& baryc_coord);
	//! check the barycentric is valid or not
	bool IsValidBarycentic(const Barycentrc& baryc_coord);


	//! compute the error out of valid barycentric
	double ComputeErrorOutValidBarycentric(const Barycentrc& baryc_coord);
}

#endif // BARYCENTRIC_H_