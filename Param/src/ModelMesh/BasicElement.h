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
// BasicElement.h
//
// [Developer]
// Xu Dong
// State Key Lab of CAD&CG, Zhejiang University
// 
// [Date]
// 2005-08-05
//
// [Goal]
// Defining 2 kinds of basic elements and 3 kinds of auxiliary elements
// 2 kinds of basic elements -- Edge & Face
// 3 kinds of auxiliary elements -- Point & Line & Polygon



#include "../Common/Macro.h"
#include "../Common/BasicDataType.h"
#pragma once

#include <vector>


/* ================== Basic Element - Mesh Edge ================== */

class Edge 
{
private:
	VertexID v1, v2;

public:
	Edge(VertexID v1 = -1, VertexID v2 = -1) { 
		this->v1 = v1;
		this->v2 = v2; 
	}
	void set(VertexID v1, VertexID v2) { 
		this->v1 = v1;
		this->v2 = v2; 
	}
	VertexID getV1() {
		return v1;
	}
	VertexID getV2() {
		return v2;
	}
	bool operator==(const Edge &edge) {
		if (edge.v1==v1 && edge.v2==v2) return true;
		if (edge.v1==v2 && edge.v2==v1) return true;
		return false;
	}
	void operator=(const Edge &edge) {
		v1=edge.v1;
		v2=edge.v2;
	}
};

/* ================== Basic Element - Mesh Face ================== */
// This class is used for defining new faces. 
// The orientation of the vertices is: v1 v2 v3 .. vn
class Face 
{
private:
	std::vector<VertexID> v;

public:
    // Constructor
    Face()
    {
        v.clear();
    }

	Face(std::vector<VertexID>& fv) 
    {
        set(fv);
	}

    // Initializer
	void set(std::vector<VertexID> fv) 
    { 
        v.clear();
        v = fv;
	}

    // Method
	VertexID getV(size_t vNum) const 
    {
		if (vNum>=v.size() || vNum<0)
			return -1;
		return v[vNum];
	}

    size_t getVertexNum() const
    {
        return v.size();
    }

    // Operator
	void operator =(const Face &face) 
    {
        v.clear();
        size_t nv = face.getVertexNum();

		for (size_t i = 0; i < nv; ++ i) 
        {
			v[i] = face.getV(i);
		}
	}

	bool operator ==(const Face &face) 
    {
        size_t nv = face.getVertexNum();
        if(nv != getVertexNum())
            return false;

		for (size_t i = 0; i < nv; ++ i) 
        {
			if (face.getV(i) != getV(i)) 
                return false;
		}
		return true;
	}
};



/* ================== Auxiliary Element - Point ================== */

class Point
{
private:
	Coord v;

public:
	Point() {}
	Point(const Coord& coord)
	{
		for(int i = 0; i < 3; i ++)
			v[i] = coord[i];
	}
	void setCoord(const Coord& coord)
	{
		for(int i = 0; i < 3; i ++)
			v[i] = coord[i];
	}
	void setCoord(const Coord coord)
	{
		for(int i = 0; i < 3; i ++)
			v[i] = coord[i];
	}
	void setCoord(const Coord* _coord)
	{
		if(_coord == NULL)
			return;
		for(int i = 0; i < 3; i ++)
			v[i] = (*_coord)[i];
	}

	void getCoord(Coord& coord)
	{
		for(int i = 0; i < 3; i ++)
			coord[i] = v[i];
	}
	void getCoord(Coord* _coord)
	{
		if(_coord == NULL)
			return;
		for(int i = 0; i < 3; i ++)
			(*_coord)[i] = v[i];
	}
};



/* ================== Auxiliary Element - Line ================== */

class Line {

private:
	Coord coords[2];

public:
	
	Line() {} 

	Line(const Coord coord1, const Coord coord2) {
		for (int i=0; i<3; i++) {
			coords[0][i] = coord1[i];
			coords[1][i] = coord2[i];
		}
	}
	void setCoords(const Coord& coord1, const Coord& coord2) {
		for(int i = 0; i < 3; i ++)
		{
			coords[0][i] = coord1[i];
			coords[1][i] = coord2[i];
		}
	}
	void setCoords(const Coord coord1, const Coord coord2) {
		for(int i = 0; i < 3; i ++)
		{
			coords[0][i] = coord1[i];
			coords[1][i] = coord2[i];
		}
	}
	void setCoords(const Coord* _coords) { 
		if (_coords==NULL) return;
		for (int i=0; i<3; i++) {
			coords[0][i] = _coords[0][i];
			coords[1][i] = _coords[1][i];
		}
	}

	void getCoords(Coord& coord1, Coord& coord2) {
		for(int i = 0; i < 3; i ++)
		{
			coord1[i] = coords[0][i];
			coord2[i] = coords[1][i];
		}
	}
	void getCoords(Coord* coord1, Coord* coord2) {
		for(int i = 0; i < 3; i ++)
		{
			(*coord1)[i] = coords[0][i];
			(*coord2)[i] = coords[1][i];
		}
	}
	void getCoords(Coord* _coords)  { 
		if (_coords==NULL) return;
		for (int i=0; i<3; i++) {
			_coords[0][i] = coords[0][i];
			_coords[1][i] = coords[1][i];
		}
	}
};



/* ================== Auxiliary Element - Polygon ================== */

class Polygon
{
private:
	Coord coords[3];

public:
	
	Polygon() {} 

	Polygon(const Coord coord1, const Coord coord2, const Coord coord3) {
		for (int i=0; i<3; i++) {
			coords[0][i] = coord1[i];
			coords[1][i] = coord2[i];
			coords[2][i] = coord3[i];
		}
	}

	void setCoords(const Coord& coord1, const Coord& coord2, const Coord& coord3) {
		for(int i = 0; i < 3; i ++)
		{
			coords[0][i] = coord1[i];
			coords[1][i] = coord2[i];
			coords[2][i] = coord3[i];
		}
	}
	void setCoords(const Coord coord1, const Coord coord2, const Coord coord3) {
		for(int i = 0; i < 3; i ++)
		{
			coords[0][i] = coord1[i];
			coords[1][i] = coord2[i];
			coords[2][i] = coord3[i];
		}
	}
	void setCoords(const Coord* _coords) { 
		if (_coords==NULL) return;
		for (int i=0; i<3; i++) {
			coords[0][i] = _coords[0][i];
			coords[1][i] = _coords[1][i];
			coords[2][i] = _coords[2][i];
		}
	}

	void getCoords(Coord& coord1, Coord& coord2, Coord& coord3) {
		for(int i = 0; i < 3; i ++)
		{
			coord1[i] = coords[0][i];
			coord2[i] = coords[1][i];
			coord3[i] = coords[2][i];
		}
	}
	void getCoords(Coord* coord1, Coord* coord2, Coord* coord3) {
		for(int i = 0; i < 3; i ++)
		{
			(*coord1)[i] = coords[0][i];
			(*coord2)[i] = coords[1][i];
			(*coord3)[i] = coords[2][i];
		}
	}
	void getCoords(Coord* _coords)  { 
		if (_coords==NULL) return;
		for (int i=0; i<3; i++) {
			_coords[0][i] = coords[0][i];
			_coords[1][i] = coords[1][i];
			_coords[2][i] = coords[2][i];
		}
	}
};


