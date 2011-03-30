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
// MeshModelBasicOp.cpp
//
// [Developer]
// Xu Dong
// State Key Lab of CAD&CG, Zhejiang University
// 
// [Date]
// 2005-08-05
//
// [Goal]
// Defining the basic operations of the kernel components of the mesh model
// Including predictions, queries, Euler operations, etc
//
// This class is responsible for the manipulations of the kernel components and states

#include "MeshModelBasicOp.h"
#include <stack>
#include <algorithm>
#include "../Numerical/Heap.h"
#include <cmath>
#include <cassert>
#include <fstream>
using namespace std;

// Constructor
MeshModelBasicOp::MeshModelBasicOp()
{
    kernel = NULL ;
    auxdata = NULL;
}

// Destructor
MeshModelBasicOp::~MeshModelBasicOp()
{

}

// Initializer
void MeshModelBasicOp::ClearData()
{
    
}

void MeshModelBasicOp::AttachKernel(MeshModelKernel* pKernel)
{
    assert(pKernel != NULL);
    kernel = pKernel;
}

void MeshModelBasicOp::AttachAuxData(MeshModelAuxData* pAuxData)
{
    assert(pAuxData != NULL);
    auxdata = pAuxData;
}


/* ================== Vertex information calculation ================== */

// Calculate the adjacent face/vertex information for each vertex
void MeshModelBasicOp::CalAdjacentInfo()
{
    VertexInfo& vInfo = kernel->GetVertexInfo();
    FaceInfo& fInfo = kernel->GetFaceInfo();

    PolyIndexArray& vAdjVertices = kernel->GetVertexInfo().GetAdjVertices();
    PolyIndexArray& vAdjFaces = vInfo.GetAdjFaces();
    
    util.FreeVector(vAdjFaces);
    util.FreeVector(vAdjVertices);

    // Calculate the adjacent faces for each vertex
    size_t i, j, n;
    size_t nVertex = vInfo.GetCoord().size();
    vAdjFaces.resize(nVertex);

    PolyIndexArray& fIndex = fInfo.GetIndex();
    size_t nFace = fIndex.size();
    for(i = 0; i < nFace; ++ i)
    {
        IndexArray& f = fIndex[i];
        n = f.size();
        for(j = 0; j < n; ++ j)
        {
            VertexID vID = f[j];    // vID is the jth vertex of ith Face
            vAdjFaces[vID].push_back((int) i);
        }
    }

    // Calculate the adjacent vertices for each vertex
    vAdjVertices.resize(nVertex);
    for(i = 0; i < nVertex; ++ i)
    {
        // Gathering the adjacent vertices of vertex i from adjacent faces
        IndexArray& adjVertices = vAdjVertices[i];
        IndexArray& adjFaces = vAdjFaces[i];
        n = adjFaces.size();
        for(j = 0; j < n; ++ j)
        {
            FaceID fID = vAdjFaces[i][j];
            IndexArray& face = fIndex[fID];
            // Find the position of vertex i in face fID
            size_t idx = distance(face.begin(), find(face.begin(), face.end(), i));
            assert(idx >= 0 && idx < face.size());

            // Add previous and next vertices of vertex i to adjacent-vertex array
            size_t m = face.size();
            adjVertices.push_back(face[(idx+1)%m]);
            adjVertices.push_back(face[(idx+m-1)%m]);
        }

        // Validate the 1-ring neighborhood of vertex i
        // Make sure NOT contain central vertex i
        adjVertices.erase(remove(adjVertices.begin(), adjVertices.end(), i), adjVertices.end());
        // Sort for the next function unique
        sort(adjVertices.begin(), adjVertices.end());
        // Make sure NO duplicated adjacent vertex, but does not resize the vector
        // Resize the vector so that only contain all necessary elements
        adjVertices.erase(unique(adjVertices.begin(), adjVertices.end()), adjVertices.end());

        // Debug
//        util.PrintArray(adjVertices, "Adjvertices");
    }

    // Debug
//    for(i = 0; i < nVertex; ++ i)
//    {
//        int n = vAdjFaces[i].size();
//        fprintf(stdout, "Vertex [%8d], #Adjacent Faces = %2d, ", i, n);
//        n = vAdjVertices[i].size();
//        fprintf(stdout, "#Adjacent vertices = %2d\n", n);
//    }
}

// Calculate normal vector of all vertices
void MeshModelBasicOp::CalVertexNormal()
{
    CoordArray& vCoord = kernel->GetVertexInfo().GetCoord();
    NormalArray& vNormal = kernel->GetVertexInfo().GetNormal();
    NormalArray& fNormal = kernel->GetFaceInfo().GetNormal();
    PolyIndexArray& vAdjFaces = kernel->GetVertexInfo().GetAdjFaces();

    size_t nVertex = vCoord.size();
    vNormal.resize(nVertex);
    size_t i, j, n;
    for(i = 0; i < nVertex; ++ i)
    {
        n = vAdjFaces[i].size();
        Normal& vn = vNormal[i];
        vn.setCoords(0.0, 0.0, 0.0);
        for(j = 0; j < n; ++ j)
        {
            FaceID fID = vAdjFaces[i][j];
            vn += fNormal[fID];
        }

        if(!vn.normalize())
            vn = COORD_AXIS_Z;
    }
}

// Calculate normal vector of selected vertices
void MeshModelBasicOp::CalVertexNormal(IntArray& arrIndex)
{
    CoordArray& vCoord = kernel->GetVertexInfo().GetCoord();
    NormalArray& vNormal = kernel->GetVertexInfo().GetNormal();
    NormalArray& fNormal = kernel->GetFaceInfo().GetNormal();
    PolyIndexArray& vAdjFaces = kernel->GetVertexInfo().GetAdjFaces();

    size_t i, j, m, n = arrIndex.size();
    for(i = 0; i < n; ++ i)
    {
        VertexID vID = arrIndex[i];
        Normal& vn = vNormal[vID];
        IndexArray& adjFaces = vAdjFaces[vID];
        m = adjFaces.size();
        vn.setCoords(0.0, 0.0, 0.0);
        for(j = 0; j < m; ++ j)
        {
            FaceID fID = adjFaces[j];
            vn += fNormal[fID];
        }

        if(!vn.normalize())
            vn = COORD_AXIS_Z;
    }
}



/* ================== Face information calculation ================== */

// Calculate normal vector of all faces
void MeshModelBasicOp::CalFaceNormal()
{
    CoordArray& vCoord = kernel->GetVertexInfo().GetCoord();
    PolyIndexArray& fIndex = kernel->GetFaceInfo().GetIndex();

    size_t nVertex = vCoord.size();
    size_t nFace = fIndex.size();

    NormalArray& fNormal = kernel->GetFaceInfo().GetNormal();
    fNormal.resize(nFace);
    size_t i, j, n;
    Coord v[3];
    for(i = 0; i < nFace; ++ i)
    {
        IndexArray& f = fIndex[i];
        n = f.size();
        for(j = 0; j < 3; ++ j)
            v[j] = vCoord[f[j]];
        Normal& fn = fNormal[i];
        fn = cross(v[1]-v[0], v[2]-v[0]);
        if(!fn.normalize())
            fn = COORD_AXIS_Z;
    }
}

// Calculate normal vector of selected faces
void MeshModelBasicOp::CalFaceNormal(IntArray& arrIndex)
{
    CoordArray& vCoord = kernel->GetVertexInfo().GetCoord();
    PolyIndexArray& fIndex = kernel->GetFaceInfo().GetIndex();
    NormalArray& fNormal = kernel->GetFaceInfo().GetNormal();

    Coord v[3];
    size_t i, j, n = arrIndex.size();
    for(i = 0; i < n; ++ i)
    {
        FaceID fID = arrIndex[i];
        IndexArray& f = fIndex[fID];
        for(j = 0; j < 3; ++ j)
            v[j] = vCoord[f[j]];
        Normal& fn = fNormal[fID];
        fn = cross(v[1]-v[0], v[2]-v[0]);
        if(!fn.normalize())
            fn = COORD_AXIS_Z;
    }
}

void MeshModelBasicOp::CalFaceBaryCenter()
{
	// cal the BaryCenter of each face.
	FaceInfo& fInfo = kernel->GetFaceInfo();
	PolyIndexArray& polyIndexArray = fInfo.GetIndex();
	VertexInfo& vInfo = kernel->GetVertexInfo();
	CoordArray& arrCoord = vInfo.GetCoord();
	size_t nFace = polyIndexArray.size();

	CoordArray& fBaryCenter = kernel->GetFaceInfo().GetBaryCenter();
	fBaryCenter.resize(nFace);

	for (size_t i = 0; i < polyIndexArray.size(); i++)
	{
		IndexArray indexArray = polyIndexArray[i];
        Coord& center = fBaryCenter[i];
		center = (arrCoord[indexArray[0]] + arrCoord[indexArray[1]] + arrCoord[indexArray[2]]) / 3;
	}
}
void MeshModelBasicOp::CalFaceArea()
{
	CoordArray& vCoord = kernel->GetVertexInfo().GetCoord();
	PolyIndexArray& fIndex = kernel->GetFaceInfo().GetIndex();
	DoubleArray& faceArea = kernel->GetFaceInfo().GetFaceArea();
	
	size_t nFace = fIndex.size();
	faceArea.resize(nFace);

    size_t i, j;
	Coord v[3], fNormal;
	double area;
	
	for (j = 0; j < nFace; j++)

	{
		IndexArray& face = fIndex[j];
		
		for(i = 0; i < 3; i ++)
		{
            v[i] = vCoord[face[i]];
		}
		
		fNormal = cross(v[1]-v[0], v[2]-v[0]);
		area = fNormal.abs();
		
		faceArea[j] = area / 2;
	}
}

/* ================== Edge information calculation, optional functions ================== */

// Create halfedge, optional function
void MeshModelBasicOp::CreateHalfEdge()
{
}

// Calculate the halfedge information, optional function
void MeshModelBasicOp::CalHalfEdgeInfo()
{
}
// Calculate the edge information, optional function
void MeshModelBasicOp::CalVertexEdgeInfo()
{
	PolyIndexArray& polyIndexArray = kernel->GetEdgeInfo().GetVertexIndex();
	CoordArray& vCoord = kernel->GetVertexInfo().GetCoord();
	PolyIndexArray& adjVerticesArray = kernel->GetVertexInfo().GetAdjVertices();
	
	size_t i;
	size_t nVertex = vCoord.size();
	vector<bool> Visited;
	Visited.resize(nVertex);
	fill(Visited.begin(), Visited.end(), false);

	//
	for (i = 0; i < nVertex; i++)
	{
		IndexArray& adjVertices = adjVerticesArray[i];

		for (size_t j = 0; j < adjVertices.size(); j++)
		{
			if (!Visited[adjVertices[j]])
			{
				IndexArray edge;
				edge.clear();
				edge.push_back((int) i);
				edge.push_back(adjVertices[j]);
			    polyIndexArray.push_back(edge);
			}
		}
        Visited[i] = true;
	}
}

void MeshModelBasicOp::CalFaceEdgeInfo()
{
	PolyIndexArray& vtxIdxArray = kernel->GetEdgeInfo().GetVertexIndex();
	PolyIndexArray& faceIdxArray = kernel->GetEdgeInfo().GetFaceIndex();

	faceIdxArray.clear();
	faceIdxArray.resize(vtxIdxArray.size());

	for (size_t i = 0; i < vtxIdxArray.size(); i++)
	{
		IndexArray& edge_vtx = vtxIdxArray[i];
		IndexArray& edge_face = faceIdxArray[i];

		int fid1, fid2;
		GetAdjacentFace(edge_vtx[0], edge_vtx[1], fid1, fid2);
		if (fid1 >= 0)
		{
			edge_face.push_back(fid1);
		}
		if (fid2 >= 0)
		{
			edge_face.push_back(fid2);
		}
	}
}

// Calculate the dihedral angle information, optional function
void MeshModelBasicOp::CalDihedralAngle()
{
	PolyIndexArray& polyIndexArray = kernel->GetEdgeInfo().GetVertexIndex();
	size_t nEdge = polyIndexArray.size();
	DoubleArray& angleArray = kernel->GetEdgeInfo().GetDihedralAngle();
	angleArray.resize(nEdge);
	NormalArray& noramlArray = kernel->GetFaceInfo().GetNormal();

	// calculate the dihedral angle here.
	for (size_t i = 0; i < nEdge; i++)
	{
		IndexArray& indexArray = polyIndexArray[i];
		FaceID faceID1, faceID2;
		GetAdjacentFace(indexArray[0], indexArray[1], faceID1, faceID2);

		// if not boundary edge
		if (faceID1 >= 0 && faceID2 >= 0)
		{
			angleArray[i] = angle(noramlArray[faceID1], noramlArray[faceID2]);
		}
		else
		{
			angleArray[i] = PI / 2.0;
		}
		//double a = angle(noramlArray[faceID1], noramlArray[faceID2]);
	}
}
/* ================== Model information calculation ================== */

// Bounding box calculation
void MeshModelBasicOp::CalBoundingBox()
{
    CoordArray& vCoord = kernel->GetVertexInfo().GetCoord();
    size_t nVertex = vCoord.size();

    size_t i;
	int j;
    Coord BoxMin, BoxMax, SphereCenter;
    double SphereRadius;
    for(j = 0; j < 3; ++ j)
        BoxMin[j] = BoxMax[j] = SphereCenter[j] = vCoord[0][j];

    // Calculate the bounding box and the center of the bounding sphere
    for(i = 1; i < nVertex; ++ i)
    {
        Coord& v = vCoord[i];
		for(j = 0; j < 3; ++ j)
		{
			if(v[j] < BoxMin[j])
				BoxMin[j] = v[j];
			if(v[j] > BoxMax[j])
				BoxMax[j] = v[j];
            SphereCenter[j] += v[j];
		}
    }
    SphereCenter /= (double)nVertex;

    // Calculate the radius of the bounding sphere
    SphereRadius = (vCoord[0]-SphereCenter).abs();
    for(i = 1; i < nVertex; ++ i)
    {
        double r = (vCoord[i]-SphereCenter).abs();
        if(r > SphereRadius)
            SphereRadius = r;
    }

    // Update corresponding model information
    ModelInfo& mInfo = kernel->GetModelInfo();
    mInfo.SetBoundingBox(BoxMin, BoxMax, BoxMax-BoxMin);
    mInfo.SetBoundingSphere(SphereCenter, SphereRadius);

	printf("Model BoundingSphere Center : (%lf %lf %lf), Radius : %lf\n",
		SphereCenter[0], SphereCenter[1], SphereCenter[2], SphereRadius);
}

// Update bounding box according to an array of updated vertex coords
void MeshModelBasicOp::CalBoundingBox(IntArray& arrIndex)
{
    CoordArray& vCoord = kernel->GetVertexInfo().GetCoord();
    size_t nVertex = vCoord.size();

    // Get current bounding box information
    Coord BoxMin, BoxMax, BoxDim, Center;
    double radius;
    ModelInfo& mInfo = kernel->GetModelInfo();
    mInfo.GetBoundingBox(BoxMin, BoxMax, BoxDim);
    mInfo.GetBoundingSphere(Center, radius);

    Center *= (double)nVertex;
    size_t i, j, n = arrIndex.size();
    for(i = 0; i < n; ++ i)
    {
        VertexID vID = arrIndex[i];
        Coord v = vCoord[vID];
        for(j = 0; j < 3; ++ j)
        {
            if(v[j] < BoxMin[j])
                BoxMin[j] = v[j];
            if(v[j] > BoxMax[j])
                BoxMax[j] = v[j];
            Center[j] += v[j];
        }
    }
    Center /= (double)(nVertex+n);

    for(i = 0; i < n; ++ i)
    {
        VertexID vID = arrIndex[i];
        Coord v = vCoord[vID];
        double r = (v-Center).abs();
        if(r > radius)
            radius = r;
    }

    // Update corresponding model information
    mInfo.SetBoundingBox(BoxMin, BoxMax, BoxMax-BoxMin);
    mInfo.SetBoundingSphere(Center, radius);
}

// Component calculation
void MeshModelBasicOp::CalComponentInfo()
{
    PolyIndexArray& vAdjVertices = kernel->GetVertexInfo().GetAdjVertices();

    size_t nVertex = vAdjVertices.size();
    size_t i, j, n;

    vector<bool> Visited;
    Visited.resize(nVertex);
    fill(Visited.begin(), Visited.end(), false);

    int nComponent = 0;
    for(i = 0; i < nVertex; ++ i)
    {
        if(Visited[i])
            continue;

        stack<VertexID> Stack;
        Stack.push((int) i);
        Visited[i] = true;
        while(Stack.size())
        {
            VertexID vID = Stack.top();
            Stack.pop();

            IndexArray& adjVertices = vAdjVertices[vID];
            n = adjVertices.size();
            for(j = 0; j < n; ++ j)
            {
                if(!Visited[adjVertices[j]])
                {
                    Stack.push(adjVertices[j]);
                    Visited[adjVertices[j]] = true;
                }
            }
        }
        ++ nComponent;
    }

    kernel->GetModelInfo().SetComponentNum(nComponent);
    printf("Number of Components = %4d\n", nComponent);
}



/* ================== Model analysis functions ================== */

// Return the number of adjacent faces for given edge vID_vID2
int MeshModelBasicOp::AdjFaceNum(VertexID vID, VertexID vID2)
{
    int fID;
	int nAdjFace = 0;
    PolyIndexArray& vAdjFaces = kernel->GetVertexInfo().GetAdjFaces();
    PolyIndexArray& fIndex = kernel->GetFaceInfo().GetIndex();

	size_t n = vAdjFaces[vID].size();
	for(size_t i = 0; i < n; ++ i)
	{
		fID = vAdjFaces[vID][i];
        IndexArray& face = fIndex[fID];
        if(find(face.begin(), face.end(), vID2) != face.end()) // find end point in current polygon
			nAdjFace ++;
	}
	return nAdjFace;
}

void MeshModelBasicOp::TopologyAnalysis()
{

    CoordArray& vCoord = kernel->GetVertexInfo().GetCoord();
    FlagArray& vFlag = kernel->GetVertexInfo().GetFlag();
    PolyIndexArray& vAdjFaces = kernel->GetVertexInfo().GetAdjFaces();
    PolyIndexArray& vAdjVertices = kernel->GetVertexInfo().GetAdjVertices();
    PolyIndexArray& fIndex = kernel->GetFaceInfo().GetIndex();
    FlagArray& fFlag = kernel->GetFaceInfo().GetFlag();

    size_t nVertex = vCoord.size();
    size_t nFace = fIndex.size();

    size_t i, j, n;
    vFlag.resize(nVertex);
    fFlag.resize(nFace);
    fill(vFlag.begin(),vFlag.end(), 0);
    fill(fFlag.begin(),fFlag.end(), 0);

	ofstream fout("bad_topogloy.txt");

    bool bManifoldModel = true;
    for(i = 0; i < nVertex; ++ i)
    {
        IndexArray& adjFaces = vAdjFaces[i];
        Flag& flag = vFlag[i];

        n = adjFaces.size();
        if(!n)  // Isolated vertex
        {
            util.SetFlag(flag, VERTEX_FLAG_ISOLATED);
            continue;
        }

        // Check topology for the 1-ring neighborhood of vertex i
        IndexArray& adjVertices = vAdjVertices[i];
        n = adjVertices.size();
        int nBdyFace = 0, nNonManifoldFace = 0;
        for(j = 0; j < n; ++ j)
        {
            switch(AdjFaceNum((int) i, adjVertices[j]))
            {
            case 1:     // Boundary
                ++ nBdyFace;
		//		auxdata->AddLine(vCoord[i], vCoord[adjVertices[j]], DARK_GREEN);
				fout << i+1 <<' '<<adjVertices[j]+1 <<endl;
                break;
            case 2:     // 2-Manifold
                break;
            default:    // Non-Manifold
                ++ nNonManifoldFace;
				auxdata->AddLine(vCoord[i], vCoord[adjVertices[j]], DARK_GREEN);
				fout << i+1 <<' '<<adjVertices[j]+1 <<endl;
            }
			
        }
        if(nNonManifoldFace || nBdyFace > 2)    // Non-manifold vertex
        {
            Coord v = vCoord[i];
            auxdata->AddPoint(v, DARK_RED);
			fout<< i+1 << endl;
            bManifoldModel = false;
            continue;
        }

// 		auxdata->AddPoint(vCoord[0], DARK_RED);
// 		auxdata->AddPoint(vCoord[15], DARK_GREEN);
// 		auxdata->AddPoint(vCoord[240], DARK_BLUE);
// 		auxdata->AddPoint(vCoord[255], DARK_GREY);
        
        // Manifold vertex
        util.SetFlag(flag, VERTEX_FLAG_MANIFOLD);
        if(nBdyFace == 2)   // Boundary vertex
            util.SetFlag(flag, VERTEX_FLAG_BOUNDARY);
    }
	fout.close();

    // Set face flag
    bool bTriMesh = true;
    bool bQuadMesh = true;
    for(i = 0; i < nFace; ++ i)
    {
        IndexArray& f = fIndex[i];
        n = f.size();

        if(n != 3)
            bTriMesh = false;
        else if(n != 4)
            bQuadMesh = false;

        int nBdyVtx = 0;
        bool bManifoldFace = true;
        for(j = 0; j < n; ++ j)
        {
            if(!util.IsSetFlag(vFlag[f[j]], VERTEX_FLAG_MANIFOLD))
            {
                bManifoldFace = false;
                break;
            }
            else if(util.IsSetFlag(vFlag[f[j]], VERTEX_FLAG_BOUNDARY))
                nBdyVtx ++;
        }
        Flag& flag = fFlag[i];
        if(bManifoldFace)
        {
            util.SetFlag(flag, FACE_FLAG_MANIFOLD);
            if(nBdyVtx > 1)
                util.SetFlag(flag, FACE_FLAG_BOUNDARY);
        }
    }

    // Set model flag
    Flag& mFlag = kernel->GetModelInfo().GetFlag();
    mFlag = 0;
    if(bManifoldModel)  // Manifold model
    {
        util.SetFlag(mFlag, MODEL_FLAG_MANIFOLD);
        printf("Manifold Model\n");
    }
    if(bTriMesh)
    {
        util.SetFlag(mFlag, MODEL_FLAG_TRIMESH);
        printf("Triangle Mesh\n");
    }
    else if(bQuadMesh)
    {
        util.SetFlag(mFlag, MODEL_FLAG_QUADMESH);
        printf("Quadangle Mesh\n");
    }
    else
    {
        util.SetFlag(mFlag, MODEL_FLAG_GENERALMESH);
        printf("General Mesh\n");
    }
}

        

/* ================== Manifold functions ================== */

// Make sure the 1-ring neighbors are CCW order
void MeshModelBasicOp::SortAdjacentInfo()
{
    CoordArray& vCoord = kernel->GetVertexInfo().GetCoord();
    FlagArray& vFlag = kernel->GetVertexInfo().GetFlag();
    PolyIndexArray& vAdjFaces = kernel->GetVertexInfo().GetAdjFaces();
    PolyIndexArray& vAdjVertices = kernel->GetVertexInfo().GetAdjVertices();
    PolyIndexArray& fIndex = kernel->GetFaceInfo().GetIndex();
    FlagArray& fFlag = kernel->GetFaceInfo().GetFlag();

    size_t nVertex = vCoord.size();
    size_t nFace = fIndex.size();
    
    size_t i, j, n, m;
	size_t k;
    for(i = 0; i < nVertex; ++ i)
    {
        Flag& flag = vFlag[i];
        if(!util.IsSetFlag(flag, VERTEX_FLAG_MANIFOLD))  // Non-manifold vertex
            continue;
        if(util.IsSetFlag(flag, VERTEX_FLAG_ISOLATED))  // Isolated vertex
            continue;

        // Only for manifold vertex
        IndexArray& adjFaces = vAdjFaces[i];
        n = adjFaces.size();
        assert(n > 0);

        // Find the start face
        FaceID start_fID = adjFaces[0];
        if(util.IsSetFlag(flag, VERTEX_FLAG_BOUNDARY))   // Boundary vertex
        {
            start_fID = -1;
            for(j = 0; j < n; ++ j)
            {
                FaceID fID = adjFaces[j];
                if(!util.IsSetFlag(fFlag[fID], FACE_FLAG_BOUNDARY))
                    continue;
                IndexArray& face = fIndex[fID];
                size_t idx = distance(face.begin(), find(face.begin(), face.end(), i));
                VertexID next_vID = face[(idx+1)%face.size()];
                if(util.IsSetFlag(vFlag[next_vID], VERTEX_FLAG_BOUNDARY))   // Find it
                {
                    if(AdjFaceNum((int) i, next_vID) == 1)    // Make sure it is a boundary edge (i, next_vID)
                    {
                        start_fID = fID;
                        break;
                    }
                }
            }
            assert(start_fID != -1);
        }

        // Iteratively find the next neighboring face
        IndexArray Sorted;
        Sorted.reserve(n);
        
        IndexArray Idx;     // Record the index of vertex i in each sorted adjacent face
        Idx.reserve(n);
        for(j = 0; j < n; ++ j)
        {
            IndexArray& face = fIndex[start_fID];
            m = face.size();
            size_t idx = distance(face.begin(), find(face.begin(), face.end(), i));
            Idx.push_back((int) idx);
            
            // Vertex i may have only 1 adjacent face
            if(j == n-1)    // Last one, Not need to find the next neighboring face
                continue;

            Sorted.push_back(start_fID);
            adjFaces.erase(remove(adjFaces.begin(), adjFaces.end(), start_fID), adjFaces.end());
        
            VertexID prev_vID = face[(idx+m-1)%m];
            for(k = 0; k < adjFaces.size(); ++ k)
            {
                FaceID fID = adjFaces[k];
                IndexArray& f = fIndex[fID];
                size_t idx = distance(f.begin(), find(f.begin(), f.end(), i));
                VertexID next_vID = f[(idx+1)%f.size()];
                if(next_vID == prev_vID)    // Next neighboring face
                {
                    start_fID = fID;
                    break;
                }
            }
        }
        assert(adjFaces.size() >= 1);
        Sorted.push_back(adjFaces[0]);
        adjFaces.clear();
        adjFaces = Sorted;

        // Set sorted adjacent vertices for vertex i
        IndexArray& adjVertices = vAdjVertices[i];
        adjVertices.clear();
        n = adjFaces.size();
        for(j = 0; j < n; ++ j)
        {
            IndexArray& face = fIndex[adjFaces[j]];
            int idx = Idx[j];
            m = face.size();
            VertexID vID = face[(idx+1)%m];
            adjVertices.push_back(vID);
        }
        if(util.IsSetFlag(vFlag[i], VERTEX_FLAG_BOUNDARY))  // Boundary vertex, add one more adjacent vertex
        {
            IndexArray& face = fIndex[adjFaces[n-1]];
            int idx = Idx[n-1];
            m = face.size();
            VertexID vID = face[(idx+m-1)%m];
            adjVertices.push_back(vID);
        }
    }
}

// Boundary calculation
void MeshModelBasicOp::CalBoundaryInfo()
{
    CoordArray& vCoord = kernel->GetVertexInfo().GetCoord();
    FlagArray& vFlag = kernel->GetVertexInfo().GetFlag();
    PolyIndexArray& vAdjVertices = kernel->GetVertexInfo().GetAdjVertices();
    PolyIndexArray& Boundaries = kernel->GetModelInfo().GetBoundary();

    util.FreeVector(Boundaries);

    // Extract boundaries
    int i;
    int nVertex = (int) vCoord.size();
    IndexArray Bdy;
    
    // Allocate a bitset with nVertex bits initialized to zero
    vector<bool> visited;
    visited.resize(nVertex);
    fill(visited.begin(), visited.end(), false);

    for(i = 0; i < nVertex; ++ i)
    {
        if(!util.IsSetFlag(vFlag[i], VERTEX_FLAG_BOUNDARY))
            continue;
        if(visited[i])
            continue;
        
        Bdy.clear();
        VertexID start_vID = i;
        VertexID curr_vID = start_vID;
        VertexID next_vID = -1;
        do 
        {
        	Bdy.push_back(curr_vID);
            IndexArray& adjVertices = vAdjVertices[curr_vID];
            next_vID = adjVertices[0];
            if(next_vID != start_vID)
                curr_vID = next_vID;
            else
                break;
        } while(1);
        Boundaries.push_back(Bdy);

        size_t n = Bdy.size();
        for(size_t j = 0; j < n; ++ j)
            visited[Bdy[j]] = true;
    }
    
    // Update corresponding model information
    kernel->GetModelInfo().SetBoundaryNum((int) Boundaries.size());
    printf("Number of Boundaries = %4d\n", Boundaries.size());
}

// Topology analisis and initialization
// Analyzing and initializing the topological structure of the mesh model
// Currently, VertexInfo::m_TexCoord and VertexInfo::m_Color and FaceInfo::m_Color are not initialized.
void MeshModelBasicOp::InitModel()
{
    if(kernel->GetModelInfo().GetType() != MODEL_TYPE_POLYGON_SOAP)
        return;
    
    printf("Analyze the mesh model...\n");

    // Calculating the adjacent information for each vertex - the basic topological data structure
    CalAdjacentInfo();

	// Bounding box and bounding sphere calculations
	CalBoundingBox();

	// Check and Reorder the triangle vertex to CCW


    // Normal calculations
    CalFaceNormal();
    CalVertexNormal();
	CalFaceBaryCenter();
	CalFaceArea();

    // Calculating number of connected components
    CalComponentInfo();

    // Analyzing the topology of the model, manifold or not
    TopologyAnalysis();
    
    if(kernel->GetModelInfo().IsManifold())
    {
        SortAdjacentInfo(); // Sorting the adjacent face/vertex information for each vertex
        CalBoundaryInfo();  // Extracting the boundaries of the model
    }

	// cal avg edge length here.
	kernel->GetModelInfo().SetAvgEdgeLength(GetAvgEdgeLength());

	// cal the dihedral angle here.
	CalVertexEdgeInfo();
	CalFaceEdgeInfo();
	CalDihedralAngle();

	// cal the vertex curature here.
	//CalVertexCurvature();

	//m_VertexFlag.resize(kernel->GetVertexInfo().GetCoord().size());

	/// 
	int vert_num = (int)kernel->GetVertexInfo().GetCoord().size();
	int face_num = (int)kernel->GetFaceInfo().GetIndex().size();

	kernel->GetModelInfo().SetVertexNum(vert_num);
	kernel->GetModelInfo().SetFaceNum(face_num);

    printf("\n");
}

// Shortest path between two given vertices
void MeshModelBasicOp::GetShortestPath(VertexID vStart, VertexID vEnd, IndexArray& Path)
{
    assert(IsValidVertexIndex(vStart) && IsValidVertexIndex(vEnd));

    if(vStart == vEnd)
    {
        Path.clear();
        Path.push_back(vStart);
        return;
    }

    PolyIndexArray& vAdjVertices = kernel->GetVertexInfo().GetAdjVertices();
    CoordArray& vCoord = kernel->GetVertexInfo().GetCoord();
    FlagArray& vFlag = kernel->GetVertexInfo().GetFlag();

    int nVertex = (int) vAdjVertices.size();
    int j, n;
    DoubleArray m_VtxDist;
    IndexArray  m_VtxParent;
    m_VtxDist.resize(nVertex);
    m_VtxParent.resize(nVertex);
    fill(m_VtxDist.begin(), m_VtxDist.end(), INFINITE_DISTANCE);
    fill(m_VtxParent.begin(), m_VtxParent.end(), -1);

	CHeap heap;

	m_VtxDist[vEnd] = 0.0;
	Node* pNode = new Node;
	pNode->v = 0.0;
	pNode->type = vEnd;
	heap.insert(pNode);
	
	// Gather all neighborhood vertices
	while(!heap.heapEmpty())
	{
		Node* pNode = heap.Remove();
		VertexID vID = pNode->type;
		double v_dist = m_VtxDist[vID];
		Coord v = vCoord[vID];
        IndexArray& adjVertices = vAdjVertices[vID];
        n = (int) adjVertices.size();
        for(j = 0; j < n; ++ j)
		{
			VertexID vtxID = adjVertices[j];
			double vtx_dist = m_VtxDist[vtxID];
			Coord vtx = vCoord[vtxID];
			double edge_length = (vtx-v).abs();
			if(v_dist+edge_length < vtx_dist)		// Update
			{
				m_VtxDist[vtxID] = v_dist+edge_length;
				m_VtxParent[vtxID] = vID;
				int index = heap.heapFind(vtxID);
				if(index)	// Already in heap
				{
					heap.a[index]->v = m_VtxDist[vtxID];
					heap.upheap(index);
				}
				else
				{
					Node* pNewNode = new Node;
					pNewNode->v = m_VtxDist[vtxID];
					pNewNode->type = vtxID;
					heap.insert(pNewNode);
				}
			}
		}
		delete pNode;
		if(vID == vStart)
			break;
	}
	
	while(!heap.heapEmpty())
	{
		Node* pNode = heap.Remove();
		delete pNode;
	}

	// Extract Vertex Path
	Path.clear();
	int curr_vID = vStart;
	while(m_VtxParent[curr_vID] != -1)
	{
		Path.push_back(curr_vID);
		curr_vID = m_VtxParent[curr_vID];
	}
    Path.push_back(vEnd);
}

// Flood fill from a seed to marked boundary
void MeshModelBasicOp::SurfaceFloodFillVertex(VertexID vID, IndexArray& FillVtx, IndexArray& FillFace, IndexArray& SelectedVtx)
{
    assert(IsValidVertexIndex(vID));

    PolyIndexArray& vAdjVertices = kernel->GetVertexInfo().GetAdjVertices();
    PolyIndexArray& vAdjFaces = kernel->GetVertexInfo().GetAdjFaces();
	PolyIndexArray& fIndex = kernel->GetFaceInfo().GetIndex();

    size_t nVertex = vAdjVertices.size();
    size_t nFace = fIndex.size();

    BoolArray VtxVisited, FaceVisited;
    VtxVisited.resize(nVertex);
    FaceVisited.resize(nFace);
    fill(VtxVisited.begin(), VtxVisited.end(), false);
    fill(FaceVisited.begin(), FaceVisited.end(), false);
    VtxVisited[vID] = true;

    // Mark boundary vertices
    size_t i, n = SelectedVtx.size();
    for(i = 0; i < n; ++ i)
        VtxVisited[SelectedVtx[i]] = true;

    // Gather filled vertices
    stack<VertexID> Stack;
    Stack.push(vID);

    while(!Stack.empty())
    {
        VertexID vID = Stack.top();
        Stack.pop();
        IndexArray& adjVertices = vAdjVertices[vID];
        n = adjVertices.size();
        for(i = 0; i < n; ++ i)
        {
            VertexID vtxID = adjVertices[i];
            if(!VtxVisited[vtxID])
            {
                Stack.push(vtxID);
                VtxVisited[vtxID] = true;
            }
        }

        FillVtx.push_back(vID);
    }

    // Gather filled faces
    n = FillVtx.size();
    for(i = 0; i < n; ++ i)
    {
        VertexID vID = FillVtx[i];
        IndexArray& adjFaces = vAdjFaces[vID];
        size_t j, m = adjFaces.size();
        for(j = 0; j < m; ++ j)
        {
            FaceID fID = adjFaces[j];
            if(!FaceVisited[fID])
            {
                FillFace.push_back(fID);
                FaceVisited[fID] = true;
            }
        }
    }
}

void MeshModelBasicOp::SurfaceFloodFillFace(FaceID fID, IndexArray& FillVtx, IndexArray& FillFace, IndexArray& SelectedFace)
{
    assert(IsValidFaceIndex(fID));

    PolyIndexArray& vAdjVertices = kernel->GetVertexInfo().GetAdjVertices();
    PolyIndexArray& vAdjFaces = kernel->GetVertexInfo().GetAdjFaces();
	PolyIndexArray& fIndex = kernel->GetFaceInfo().GetIndex();

    size_t nVertex = vAdjVertices.size();
    size_t nFace = fIndex.size();

    BoolArray VtxVisited, FaceVisited;
    VtxVisited.resize(nVertex);
    FaceVisited.resize(nFace);
    fill(VtxVisited.begin(), VtxVisited.end(), false);
    fill(FaceVisited.begin(), FaceVisited.end(), false);
    FaceVisited[fID] = true;

    // Mark boundary faces
    size_t i, n = SelectedFace.size();
    for(i = 0; i < n; ++ i)
        FaceVisited[SelectedFace[i]] = true;

    // Gather filled vertices
    stack<VertexID> Stack;
    Stack.push(fID);

    while(!Stack.empty())
    {
        FaceID fID = Stack.top();
        Stack.pop();
        IndexArray& face = fIndex[fID];
        n = face.size();
        for(i = 0; i < n; ++ i)
        {
            VertexID vID = face[i];
            IndexArray& adjFaces = vAdjFaces[vID];
            size_t j, m = adjFaces.size();
            for(j = 0; j < m; ++ j)
            {
                FaceID faceID = adjFaces[j];
                if(!FaceVisited[faceID])
                {
                    Stack.push(faceID);
                    FaceVisited[faceID] = true;
                }
            }
        }

        FillFace.push_back(fID);
    }

    // Gather filled faces
    n = FillFace.size();
    for(i = 0; i < n; ++ i)
    {
        FaceID fID = FillFace[i];
        IndexArray& face = fIndex[fID];
        size_t j, m = face.size();
        for(j = 0; j < m; ++ j)
        {
            VertexID vID = face[j];
            if(!VtxVisited[vID])
            {
                FillVtx.push_back(vID);
                VtxVisited[vID] = true;
            }
        }
    }
}

// Get the neighbors of a given vertex inside a circle with a given radius
// DO NOT clear out arrays
void MeshModelBasicOp::GetNeighborhood(VertexID vID, IndexArray& NeiVtx, IndexArray& NeiFace, double radius)
{
    assert(IsValidVertexIndex(vID));
    
    PolyIndexArray& vAdjVertices = kernel->GetVertexInfo().GetAdjVertices();
    PolyIndexArray& vAdjFaces = kernel->GetVertexInfo().GetAdjFaces();
    CoordArray& vCoord = kernel->GetVertexInfo().GetCoord();
	PolyIndexArray& fIndex = kernel->GetFaceInfo().GetIndex();
    size_t nVertex = vCoord.size();

    // Multiply the distance factor
	radius *= GetDistanceFactor();

    size_t i, j, n;
    DoubleArray m_VtxDist;
    IndexArray  m_VtxParent;
    m_VtxDist.resize(nVertex);
    m_VtxParent.resize(nVertex);
    fill(m_VtxDist.begin(), m_VtxDist.end(), INFINITE_DISTANCE);
    fill(m_VtxParent.begin(), m_VtxParent.end(), -1);

    CHeap heap;

	m_VtxDist[vID] = 0.0;
	Node* pNode = new Node;
	pNode->v = 0.0;
	pNode->type = vID;
	heap.insert(pNode);

	// Mark selected vertices
    BoolArray VtxVisited;
    VtxVisited.resize(nVertex);
    fill(VtxVisited.begin(), VtxVisited.end(), false);
    n = NeiVtx.size();
    for(i = 0; i < n; ++ i)
        VtxVisited[NeiVtx[i]] = true;

    // Mark selected faces
    size_t nFace = fIndex.size();
    BoolArray FaceVisited;
    FaceVisited.resize(nFace);
    fill(FaceVisited.begin(), FaceVisited.end(), false);
    n = NeiFace.size();
    for(i = 0; i < n; ++ i)
        FaceVisited[NeiFace[i]] = true;

	// Gather all neighboring vertices
	while(!heap.heapEmpty())
	{
		Node* pNode = heap.Remove();
		VertexID vID = pNode->type;
		double v_dist = m_VtxDist[vID];
		Coord v = vCoord[vID];
        IndexArray& adjVertices = vAdjVertices[vID];
        n = adjVertices.size();
        for(i = 0; i < n; ++ i)
		{
			VertexID vtxID = adjVertices[i];
			double vtx_dist = m_VtxDist[vtxID];
			Coord vtx = vCoord[vtxID];
			double edge_length = (vtx-v).abs();
			if(v_dist+edge_length < vtx_dist)		// Update
			{
				m_VtxDist[vtxID] = v_dist+edge_length;
				int index = heap.heapFind(vtxID);
				if(index)	// Already in heap
				{
					heap.a[index]->v = m_VtxDist[vtxID];
                    m_VtxParent[vtxID] = vID;
					heap.upheap(index);
				}
				else if(m_VtxDist[vtxID] < radius)
				{
					Node* pNewNode = new Node;
					pNewNode->v = m_VtxDist[vtxID];
					pNewNode->type = vtxID;
					heap.insert(pNewNode);
				}
			}
		}
		delete pNode;

        // Add to neighboring vertex array
        if(!VtxVisited[vID])
        {
            NeiVtx.push_back(vID);
            VtxVisited[vID] = true;
        }

        // Add to neighboring face array
        IndexArray& adjFaces = vAdjFaces[vID];
        n = adjFaces.size();
        for(i = 0; i < n; ++ i)
        {
            FaceID fID = adjFaces[i];
            if(FaceVisited[fID])
                continue;
            IndexArray& face = fIndex[fID];
            size_t m = face.size();
            for(j = 0; j < m; ++ j)
                if(!VtxVisited[face[j]])
                    break;
            if(j == m)    // All vertices of the given face are visited
            {
                NeiFace.push_back(fID);
                FaceVisited[fID] = true;
            }
        }
	}
}

// Get the nearest vertex in a given face for a give position
void MeshModelBasicOp::GetNearestVertex(FaceID fID, Coord pos, VertexID& vID)
{
    IndexArray& face = kernel->GetFaceInfo().GetIndex()[fID];
    CoordArray& vCoord = kernel->GetVertexInfo().GetCoord();

    size_t i, n = face.size();
    DoubleArray Dist;
    Dist.resize(n);
    for(i = 0; i < n; ++ i)
    {
        Coord& v = vCoord[face[i]];
        double dist = (pos-v).abs();
        Dist[i] = dist;
    }

    size_t index = distance(Dist.begin(), min_element(Dist.begin(), Dist.end()));
    assert(index >= 0 && index < n);
    vID = face[index];
}
void MeshModelBasicOp::GetNearestVertex(Coord pos, VertexID& vID)
{
	CoordArray& vCoord = kernel->GetVertexInfo().GetCoord();

    size_t i;
    double minDist = 1e20;
    
    for(i = 0; i < vCoord.size(); ++ i)
    {
        Coord& v = vCoord[i];
        double dist = (pos - v).abs();

		if (minDist > dist)
		{
			minDist = dist;
			vID = (int) i;
		}
    }
}

// Get the neighbors of a given face
void MeshModelBasicOp::GetFaceNeighborhood(FaceID fID, IndexArray& NeiFace)
{
    assert(IsValidFaceIndex(fID));

    IndexArray& f = kernel->GetFaceInfo().GetIndex()[fID];
    PolyIndexArray& vAdjFaces = kernel->GetVertexInfo().GetAdjFaces();
    FlagArray& vFlag = kernel->GetVertexInfo().GetFlag();

    size_t i, n = f.size();
    for(i = 0; i < n; ++ i)
    {
        VertexID vID = f[i];
        Flag flag = vFlag[vID];
        IndexArray& adjFaces = vAdjFaces[vID];
        size_t m = adjFaces.size();
        size_t idx = distance(adjFaces.begin(), find(adjFaces.begin(), adjFaces.end(), fID));
        
        if(!util.IsSetFlag(flag, VERTEX_FLAG_BOUNDARY))
        {
            FaceID faceID = adjFaces[(idx+1)%m];
            if(find(NeiFace.begin(), NeiFace.end(), faceID) == NeiFace.end())
                NeiFace.push_back(faceID);
        }
        else if(idx != m-1)    // Not the last adjacent face
        {
            FaceID faceID =  adjFaces[(idx+1)%m];
            if(find(NeiFace.begin(), NeiFace.end(), faceID) == NeiFace.end())
                NeiFace.push_back(faceID);
        }
    }
}

// Geodesic distance functions
double MeshModelBasicOp::DistanceFromSeeds(IndexArray& Seeds, DoubleArray& VtxDist, double distance)
{
	size_t nVtx = Seeds.size();
	if(!nVtx)
		return -1.0;
	
    // Multiply the distance factor
    if(distance != INFINITE_DISTANCE)
	    distance *= GetDistanceFactor();

    CoordArray& vCoord = kernel->GetVertexInfo().GetCoord();
    PolyIndexArray& vAdjVertices = kernel->GetVertexInfo().GetAdjVertices();

    size_t nVertex = vCoord.size();

	VtxDist.clear();
	VtxDist.resize(nVertex);
    fill(VtxDist.begin(), VtxDist.end(), INFINITE_DISTANCE);

	CHeap heap;
    size_t i, n;

    for(i = 0; i < nVtx; ++ i)
	{
		VertexID vID = Seeds[i];
		VtxDist[vID] = 0.0;
		Node* pNode = new Node;
		pNode->v = 0.0;
		pNode->type = vID;
		heap.insert(pNode);
	}
	
	// Mark selected vertices
    BoolArray VtxVisited;
    VtxVisited.resize(nVertex);
    fill(VtxVisited.begin(), VtxVisited.end(), false);
    n = Seeds.size();
    for(i = 0; i < n; ++ i)
        VtxVisited[Seeds[i]] = true;

    double VtxMaxDist = 0.0;
    // Gather all neighborhood vertices
    Coord v, vtx;
	while(!heap.heapEmpty())
	{
		Node* pNode = heap.Remove();
		VertexID vID = pNode->type;
		double v_dist = VtxDist[vID];
		v = vCoord[vID];
        
        IndexArray& adjVertices = vAdjVertices[vID];
		n = adjVertices.size();
        for(i = 0; i < n; ++ i)
		{
			VertexID vtxID = adjVertices[i];
			double vtx_dist = VtxDist[vtxID];
			vtx = vCoord[vtxID];
			double edge_length = (vtx-v).abs();
			if(v_dist+edge_length < vtx_dist)		// Update
			{
				VtxDist[vtxID] = v_dist+edge_length;
				int index = heap.heapFind(vtxID);
				if(index)	// Already in heap
				{
					heap.a[index]->v = VtxDist[vtxID];
					heap.upheap(index);
				}
				else if(VtxDist[vtxID] < distance)
				{
					Node* pNewNode = new Node;
					pNewNode->v = VtxDist[vtxID];
					pNewNode->type = vtxID;
					heap.insert(pNewNode);
				}
			}
		}
		VtxMaxDist = VtxDist[vID];
		delete pNode;
	}

    return VtxMaxDist;
}

double MeshModelBasicOp::DistanceFromSeeds2(IndexArray& Seeds, DoubleArray& VtxDist, DoubleArray& FaceDist, double distance)
{
    DistanceFromSeeds(Seeds, VtxDist, distance);

    PolyIndexArray& fIndex = kernel->GetFaceInfo().GetIndex();
    size_t nFace = fIndex.size();

	FaceDist.clear();
	FaceDist.resize(nFace);
    fill(FaceDist.begin(), FaceDist.end(), INFINITE_DISTANCE);

    double FaceMaxDist = 0.0;

    size_t i, j, n;
    for(i = 0; i < nFace; ++ i)
	{
		IndexArray& f = fIndex[i];
		double dist = 0.0;
		bool bValid = true;
        n = f.size();
		for(j = 0; j < n; ++ j)
		{
			if(VtxDist[f[j]] == INFINITE_DISTANCE)
			{
				bValid = false;
				break;
			}
			else
				dist += VtxDist[f[j]];
		}
		if(bValid)
		{
            double& dist = FaceDist[i];
			dist = dist/3.0;
			if(!i || dist > FaceMaxDist)
				FaceMaxDist = dist;
		}
	}

    return FaceMaxDist;
}

double MeshModelBasicOp::DistanceFromSeeds3(IndexArray& Seeds, DoubleArray& FaceDist, double distance /* = INFINITE_DISTANCE */)
{
    DoubleArray VtxDist;
    DistanceFromSeeds(Seeds, VtxDist, distance);

    PolyIndexArray& fIndex = kernel->GetFaceInfo().GetIndex();
    size_t nFace = fIndex.size();

	FaceDist.clear();
	FaceDist.resize(nFace);
    fill(FaceDist.begin(), FaceDist.end(), INFINITE_DISTANCE);

    double FaceMaxDist = 0.0;

    size_t i, j, n;
    for(i = 0; i < nFace; ++ i)
	{
		IndexArray& f = fIndex[i];
		double dist = 0.0;
		bool bValid = true;
        n = f.size();
		for(j = 0; j < n; ++ j)
		{
			if(VtxDist[f[j]] == INFINITE_DISTANCE)
			{
				bValid = false;
				break;
			}
			else
				dist += VtxDist[f[j]];
		}
		if(bValid)
		{
            double& facedist = FaceDist[i];
			facedist = dist/3.0;
			if(!i || facedist > FaceMaxDist)
				FaceMaxDist = facedist;
		}
	}

    return FaceMaxDist;
}

// Compute the neighboring vertices of given seeds (not including seeds) within the given threshold
// and their corresponding distance from seeds
double MeshModelBasicOp::DistanceFromSeeds4(IndexArray& Seeds, IndexArray& NeiVtx, DoubleArray& NeiVtxDist, double distance /* = INFINITE_DISTANCE */)
{
	size_t nVtx = Seeds.size();
	if(!nVtx)
		return -1.0;
	
    // Multiply the distance factor
    if(distance != INFINITE_DISTANCE)
	    distance *= GetDistanceFactor();

    CoordArray& vCoord = kernel->GetVertexInfo().GetCoord();
    PolyIndexArray& vAdjVertices = kernel->GetVertexInfo().GetAdjVertices();
    size_t nVertex = vCoord.size();

    NeiVtx.clear();
    NeiVtxDist.clear();

    DoubleArray VtxDist;
	VtxDist.resize(nVertex);
    fill(VtxDist.begin(), VtxDist.end(), INFINITE_DISTANCE);

	CHeap heap;
    size_t i, n;

    for(i = 0; i < nVtx; ++ i)
	{
		VertexID vID = Seeds[i];
		VtxDist[vID] = 0.0;
		Node* pNode = new Node;
		pNode->v = 0.0;
		pNode->type = vID;
		heap.insert(pNode);
	}
	
	// Mark selected vertices
    BoolArray VtxVisited;
    VtxVisited.resize(nVertex);
    fill(VtxVisited.begin(), VtxVisited.end(), false);
    n = Seeds.size();
    for(i = 0; i < n; ++ i)
        VtxVisited[Seeds[i]] = true;

    double VtxMaxDist = 0.0;
    // Gather all neighborhood vertices
    Coord v, vtx;
	while(!heap.heapEmpty())
	{
		Node* pNode = heap.Remove();
		VertexID vID = pNode->type;
		double v_dist = VtxDist[vID];
		v = vCoord[vID];
        
        IndexArray& adjVertices = vAdjVertices[vID];
		n = adjVertices.size();
        for(i = 0; i < n; ++ i)
		{
			VertexID vtxID = adjVertices[i];
			double vtx_dist = VtxDist[vtxID];
			vtx = vCoord[vtxID];
			double edge_length = (vtx-v).abs();
			if(v_dist+edge_length < vtx_dist)		// Update
			{
				VtxDist[vtxID] = v_dist+edge_length;
				int index = heap.heapFind(vtxID);
				if(index)	// Already in heap
				{
					heap.a[index]->v = VtxDist[vtxID];
					heap.upheap(index);
				}
				else if(VtxDist[vtxID] < distance)
				{
					Node* pNewNode = new Node;
					pNewNode->v = VtxDist[vtxID];
					pNewNode->type = vtxID;
					heap.insert(pNewNode);

                    // Add to neighborhood
                    NeiVtx.push_back(vtxID);
				}
			}
		}
		VtxMaxDist = VtxDist[vID];
		delete pNode;
	}

    // Output
    n = NeiVtx.size();
    NeiVtxDist.resize(n);
    for(i = 0; i < n; ++ i)
        NeiVtxDist[i] = VtxDist[NeiVtx[i]];

    return VtxMaxDist;
}

double MeshModelBasicOp::GetDistanceFactor()
{
    // Multiply the radius by the 5% of the radius of the bounding sphere
    Coord Certer;
    double Radius;
    kernel->GetModelInfo().GetBoundingSphere(Certer, Radius);
	return Radius*0.05;
}

// Get the average edge length
double MeshModelBasicOp::GetAvgEdgeLength()
{
    CoordArray& vCoord = kernel->GetVertexInfo().GetCoord();
    PolyIndexArray& vAdjVertices = kernel->GetVertexInfo().GetAdjVertices();
    size_t nVertex = vCoord.size();
    
    time_t tt1;
    time(&tt1);
    srand((unsigned)tt1);

    double avg_edge = 0.0;
    int nEdge = 0;
    // Calculate by full sampling
    for(size_t i = 0; i < nVertex; ++ i)
    {
        IndexArray& adjVertices = vAdjVertices[i];
        size_t j, m = adjVertices.size();
        for(j = 0; j < m; ++ j)
        {
            VertexID vID = adjVertices[j];
            if(vID > (int) i)
                continue;
            avg_edge += (vCoord[i]-vCoord[vID]).abs();
            ++ nEdge;
        }
    }

//    // Calculate by part sampling
//    int i, n = (nVertex/10 > 0) ? nVertex/10 : nVertex;
//    for(i = 0; i < n; ++ i)
//    {
//        VertexID vID = rand()%nVertex;
//        IndexArray& adjVertices = vAdjVertices[vID];
//        int j, m = adjVertices.size();
//        for(j = 0; j < m; ++ j)
//        {
//            VertexID vtxID = adjVertices[j];
//            avg_edge += (vCoord[vtxID]-vCoord[vID]).abs();
//            ++ nEdge;
//        }
//    }

    return avg_edge/(double)nEdge;
}

// Discrete differential operator
void MeshModelBasicOp::CalMeanCurvature(CoordArray& VtxMeanCurv)
{
    if(!kernel->GetModelInfo().IsTriManifold())
        return;

    CoordArray& vCoord = kernel->GetVertexInfo().GetCoord();
    PolyIndexArray& vAdjVertices = kernel->GetVertexInfo().GetAdjVertices();
    PolyIndexArray& vAdjFaces = kernel->GetVertexInfo().GetAdjFaces();
    PolyIndexArray& fIndex = kernel->GetFaceInfo().GetIndex();
    size_t nFace = fIndex.size();
    size_t nVertex = vCoord.size();
    
    // Calculate the three inner angles of each face
    vector<float> FaceAngle, FaceLength;
    FaceAngle.resize(nFace*3);
    FaceLength.resize(nFace*3);
    size_t i, n = nFace;
	size_t k, l;
    for(i = 0; i < n; ++ i)
    {
        IndexArray& f = fIndex[i];
        size_t j, m = f.size();
        size_t idx = i*3;
        for(j = 0; j < 3; ++ j)
        {
            VertexID vID1 = f[(j+m-1)%m];
            VertexID vID2 = f[j];
            FaceLength[idx+j] = (float)(vCoord[vID1]-vCoord[vID2]).abs();
        }
        for(j = 0; j < 2; ++ j)
        {
            VertexID vID1 = f[(j+m-1)%m];
            VertexID vID2 = f[j];
            VertexID vID3 = f[(j+1)%m];
            Coord e1 = (vCoord[vID1]-vCoord[vID2])/FaceLength[idx+j];
            Coord e2 = (vCoord[vID3]-vCoord[vID2])/FaceLength[idx+(j+1)%3];
            double cos_angle = dot(e1, e2);
            if(cos_angle > 1.0)
                cos_angle = 1.0;
            else if(cos_angle < -1.0)
                cos_angle = -1.0;
            FaceAngle[idx+j] = (float)acos(cos_angle);
        }
        float& angle = FaceAngle[idx+2];
        angle = (float)(PI-FaceAngle[idx]-FaceAngle[idx+1]);
        if(angle < 0.0f)
            angle = 0.0f;
    }

    // Calculate the cotangent of angles
    n = nFace*3;
    for(i = 0; i < n; ++ i)
        FaceAngle[i] = 1.0f/tan(FaceAngle[i]);

    // Calculate the mean curvature for each vertex
    n = nVertex;
    VtxMeanCurv.resize(nVertex);
    fill(VtxMeanCurv.begin(), VtxMeanCurv.end(), Coord(0.0, 0.0, 0.0));
    for(i = 0; i < n; ++ i)
    {
        IndexArray& adjFaces = vAdjFaces[i];

        float mixed_ring_area = 0.0;
        Coord& mean_curv = VtxMeanCurv[i];

        size_t j, m = adjFaces.size();
        for(j = 0; j < m; ++ j)
        {
            FaceID fID = adjFaces[j];
            int idx = fID*3;
            for(k = 0; k < 3; ++ k)
            {
                if(FaceAngle[idx+k] < 0.0f)   // Obtuse angle, ctg(angle) < 0.0f
                    break;
            }

            IndexArray& f = fIndex[fID];
            for(l = 0; l < 3; ++ l)
            {
                if(f[l] == i)
                    break;
            }

            // Accumulate the mixed 1-ring Voronoi area
            if(k == 3)          // Non-obtuse triangle
            {
                size_t idx1 = idx+l;
                size_t idx2 = idx+(l+1)%3;
                mixed_ring_area += FaceLength[idx1]*FaceLength[idx1]*FaceAngle[idx2];
                mixed_ring_area += FaceLength[idx2]*FaceLength[idx2]*FaceAngle[idx+(l+2)%3];
                mixed_ring_area /= 8.0f;
            }
            else
            {
                // Calculate the face area using the Helen formula
                float a = FaceLength[idx];
                float b = FaceLength[idx+1];
                float c = FaceLength[idx+2];
                float p = (a+b+c)/2.0f;
                
                float sqr_area = p*(p-a)*(p-b)*(p-c);
                if(sqr_area < 0.0)
                    sqr_area = 0.0;
                float f_area = sqrt(sqr_area);

                mixed_ring_area += (k == l) ? f_area/2.0f : f_area/4.0f;     // Obtuse angle is at vertex i
            }

            // Accumulate the mean curvature operator
            mean_curv += (vCoord[f[l]]-vCoord[f[(l+1)%3]])*FaceAngle[idx+(l+2)%3];
            mean_curv += (vCoord[f[l]]-vCoord[f[(l+2)%3]])*FaceAngle[idx+(l+1)%3];
        }
        mean_curv /= (2.0f*mixed_ring_area);
    }
}

// Get the boundary index and the loop index
bool MeshModelBasicOp::GetBdyIndex(VertexID vID, BdyID& bdyID, Index& idx)
{
    assert(IsValidVertexIndex(vID));

    bdyID = idx = -1;
    if(!util.IsSetFlag(kernel->GetVertexInfo().GetFlag()[vID], VERTEX_FLAG_BOUNDARY))   // Not a boundary vertex
        return false;
    if(kernel->GetModelInfo().GetBoundaryNum() == 0)    // Closed model
        return false;

    PolyIndexArray& Boundaries = kernel->GetModelInfo().GetBoundary();
    size_t i, n = Boundaries.size();
    for(i = 0; i < n; ++ i)
    {
        IndexArray& Bdy = Boundaries[i];
        IndexArray::iterator iter = find(Bdy.begin(), Bdy.end(), vID);
        if(iter == Bdy.end())   // Not find
            continue;

        // Find it
        bdyID = (int) i;
        idx = (int) distance(Bdy.begin(), iter);
        return true;
    }
    return false;
}

// Get the adjacent face(s) of the edge (vID1, vID2)
void MeshModelBasicOp::GetAdjacentFace(VertexID vID1, VertexID vID2, FaceID& fID1, FaceID& fID2)
{
    IndexArray& adjFaces = kernel->GetVertexInfo().GetAdjFaces()[vID1];
    PolyIndexArray& fIndex = kernel->GetFaceInfo().GetIndex();
    
    fID1 = fID2 = -1;
    size_t i, n = adjFaces.size();
    for(i = 0; i < n; ++ i)
    {
        FaceID fID = adjFaces[i];
        IndexArray& f = fIndex[fID];
        IndexArray::iterator iter = find(f.begin(), f.end(), vID1);
        if(iter == f.end())
            continue;
        size_t idx = distance(f.begin(), iter);
        size_t m = f.size();
        if(f[(idx+1)%m] == vID2)
            fID1 = fID;
        else if(f[(idx+m-1)%m] == vID2)
            fID2 = fID;
    }
}
void MeshModelBasicOp::GetAdjacentFace(VertexID vID1, VertexID vID2, vector<FaceID>& face_vec)
{
	int fid_1, fid_2;
	GetAdjacentFace(vID1, vID2, fid_1, fid_2);

	//
	face_vec.clear();
	if (fid_1 >= 0)
	{
		face_vec.push_back(fid_1);
	}
	if (fid_2 >= 0)
	{
		face_vec.push_back(fid_2);
	}
}
void MeshModelBasicOp::GetAdjacentFace(FaceID fID, vector<FaceID>& face_vec)
{
	PolyIndexArray& fIndex = kernel->GetFaceInfo().GetIndex();
	IndexArray& face = fIndex[fID];

	//
	face_vec.clear();
	for (size_t i = 0; i < face.size(); i++)
	{
		int oppfid;
		GetEdgeOppositeFace(face[i], face[(i+1)%3], fID, oppfid);
		if (oppfid >= 0)
		{
			face_vec.push_back(oppfid);
		}
	}
}
void MeshModelBasicOp::GetEdgeOppositeFace(VertexID vID1, VertexID vID2, FaceID& fID, FaceID& oppfID)
{
	FaceID fid1, fid2;
	GetAdjacentFace(vID1, vID2, fid1, fid2);

	oppfID = fid1 == fID ? fid2 : fid1;
}

// Get the adjacent vertices of the edge (vID1, vID2)
void MeshModelBasicOp::GetAdjacentVertex(VertexID vID1, VertexID vID2, vector<VertexID>& oppVID)
{
    if(!IsBoundaryVertex(vID1))
    {
        IndexArray& adjVertices = kernel->GetVertexInfo().GetAdjVertices()[vID1];
        size_t idx = distance(adjVertices.begin(), find(adjVertices.begin(), adjVertices.end(), vID2));
        size_t m = adjVertices.size();
        oppVID.push_back(adjVertices[(idx+1)%m]);
        oppVID.push_back(adjVertices[(idx+m-1)%m]);
    }
    else
    { 
        IndexArray& adjFaces = kernel->GetVertexInfo().GetAdjFaces()[vID1];
        PolyIndexArray& fIndex = kernel->GetFaceInfo().GetIndex();
        
        size_t i, n = adjFaces.size();
        for(i = 0; i < n; ++ i)
        {
            FaceID fID = adjFaces[i];
            IndexArray& f = fIndex[fID];
            IndexArray::iterator iter = find(f.begin(), f.end(), vID1);
            if(iter == f.end())
                continue;
            size_t idx = distance(f.begin(), iter);
            size_t m = f.size();
            if(f[(idx+1)%m] == vID2)
                oppVID.push_back(f[(idx+m-1)%m]);
            else if(f[(idx+m-1)%m] == vID2)
                oppVID.push_back(f[(idx+1)%m]);
        }
    }
}

void MeshModelBasicOp::CalBaryCentricCoord(FaceID fID, Coord& pos, Coord& bCenter)
{
	int i;
	PolyIndexArray& fIndex = kernel->GetFaceInfo().GetIndex();
	CoordArray& vCoord =  kernel->GetVertexInfo().GetCoord();
    
	Coord v[3], fNormal;
	IndexArray& face = fIndex[fID];
	
	for(i = 0; i < 3; i ++)
	{
		v[i] = vCoord[face[i]];
	}
	
	fNormal = cross(v[1]-v[0], v[2]-v[0]);
	double fArea = fNormal.abs();
	
	double area0 = cross(v[1] - pos, v[2] - pos).abs();
	double area1 = cross(v[0] - pos, v[2] - pos).abs();
	
	bCenter[0] = area0 / fArea;
	bCenter[1] = area1 / fArea;
	bCenter[2] = 1 - bCenter[0] - bCenter[1];
}
void MeshModelBasicOp::CalTrglGradient(vector<Coord>& trgP, vector<Coord>& trgG)
{

}
void MeshModelBasicOp::CalVertexCurvature()
{
	/*
	CMeshModelCurvature meshCurv;
	meshCurv.AttachMeshModelKernel(kernel);
	meshCurv.EstimatorVertexCurvatures();
	*/
}
double MeshModelBasicOp::GetBaryAdjFaceArea(VertexID vID)
{
	DoubleArray& faceAreas = kernel->GetFaceInfo().GetFaceArea();
	PolyIndexArray& vAdjFaceArray = kernel->GetVertexInfo().GetAdjFaces();
	IndexArray& vAdjFace = vAdjFaceArray[vID];

	double area = 0.0;
	for (size_t i = 0; i < vAdjFace.size(); i++)
	{
		area += faceAreas[vAdjFace[i]];
	}

	return area / 3.0;
}
double MeshModelBasicOp::GetEdgeAdjFaceArea(VertexID vID_1, VertexID vID_2)
{
	vector<int> face_vec(2);
	GetAdjacentFace(vID_1, vID_2, face_vec[0], face_vec[1]);
	DoubleArray& faceAreas = kernel->GetFaceInfo().GetFaceArea();
	
	double area = 0;
	int num = 0;
	for (size_t i = 0; i < 2; i++)
	{
		if (face_vec[i] >= 0)
		{
			area += faceAreas[face_vec[i]];
			num++;
		}
	}

	return area / num;
}
void MeshModelBasicOp::GetNeighborhoodVertex(int vID, size_t neighRingSize, bool onlyRing, vector<int>& neighVIDs)
{
	PolyIndexArray& vAdjIndexArray = kernel->GetVertexInfo().GetAdjVertices();
	m_VertexFlag.clear(); m_VertexFlag.resize(vAdjIndexArray.size(), false);
	m_VertexFlag[vID] = true;

	vector<int> vidArray, tmpArray;
	vidArray.push_back(vID);

	for (size_t i = 0; i < neighRingSize; i++)
	{
		for (size_t j = 0; j < vidArray.size(); j++)
		{
			int vid = vidArray[j];
			IntArray& adjVIndex = vAdjIndexArray[vid];

			for (size_t k = 0; k < adjVIndex.size(); k++)
			{
				int vvvid = adjVIndex[k];

				if (!m_VertexFlag[vvvid])
				{
					tmpArray.push_back(vvvid);
					m_VertexFlag[vvvid] = true;
				}
			}
		}

		if (onlyRing)
		{
			vidArray.assign(tmpArray.begin(), tmpArray.end());
		}
		else
		{
			vidArray.insert(vidArray.end(), tmpArray.begin(), tmpArray.end());
		}

		tmpArray.clear();
	}

	neighVIDs.assign(vidArray.begin(), vidArray.end());
}
void MeshModelBasicOp::AddNoise2Model(double amp_)
{
	/*
	CoordArray& vCoord = kernel->GetVertexInfo().GetCoord();
	CoordArray& vNormal = kernel->GetVertexInfo().GetNormal();

	//
	srand((unsigned) time(NULL));
	for (size_t i = 0; i < vCoord.size(); i++)
	{
		double r_ = rand();
		double val_ = r_ / RAND_MAX - 0.5;
		val_ *= amp_;

		vCoord[i] += vNormal[i] * val_;
	}
	*/
}
