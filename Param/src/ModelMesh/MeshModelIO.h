/* ================== Library Information ================== */
// [Name] 
// MeshLib Library
//
// [Developer]
// Xu Dong
// State Key Lab of CAD&CG, Zhejiang University
// 
// [Date]
// 2005-08-05
//
// [Goal]
// A general, flexible, versatile and easy-to-use mesh library for research purpose.
// Supporting arbitrary polygonal meshes as input, but with a 
// primary focus on triangle meshes.

/* ================== File Information ================== */
// [Name]
// MeshModelIO.h
//
// [Developer]
// Xu Dong
// State Key Lab of CAD&CG, Zhejiang University
// 
// [Date]
// 2005-08-05
//
// [Goal]
// Defining the input/output (I/O) functions of the kernel components of the mesh model
// Supporting various 3D model files, including .tm, .obj, .off, etc

#ifndef MESHMODELIO_H_
#define MESHMODELIO_H_

#include <string>
#include "MeshModelKernel.h"
#include "../Common/Utility.h"
#include "../Common/TriangularMesh.h"
#include <boost/shared_ptr.hpp>


class MeshModelIO
{
private:
    MeshModelKernel* kernel;
    Utility util;

public:
    // Constructor
    MeshModelIO();

    // Destructor
    ~MeshModelIO();

    // Initializer
    void ClearData();
    void AttachKernel(MeshModelKernel* pKernel = NULL);

public:
    // General I/O functions
	bool LoadModel(const std::string& filename);
	bool StoreModel(const std::string& filename);

    // .tm file I/O functions
	bool OpenTmFile(const std::string& filename);
	bool SaveTmFile(const std::string& filename);

    // .ply2 file I/O functions
	bool OpenPly2File(const std::string& filename);
	bool SavePly2File(const std::string& filename);

	// .off file I/O functions
	bool OpenOffFile(const std::string& filename);
	bool SaveOffFile(const std::string& filename);
	
	//
	bool OpenObjFile(const std::string& filename);
	bool SaveObjFile(const std::string& filename);

	/// added by hywei
	/* \function LoadTriangularMesh
	 * a minimize dependency tri_mesh is made up by a
	 * coord array for each vertex and the face list
	*/
	bool LoadTriangularMesh(boost::shared_ptr<const tri_mesh_3d> p_mesh);
};

#endif
