#ifndef TRANSFUNCTOR_H_
#define TRANSFUNCTOR_H_

#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>
#include <hj_3rd/zjucad/matrix/matrix.h>

namespace PARAM
{
    class ChartCreator;
	class ParamCoord;

    class TransFunctor
    {
    public:
        TransFunctor(boost::shared_ptr<ChartCreator> _p_chart_creator);
        ~TransFunctor();

        //! get the transition matrix between two charts
		zjucad::matrix::matrix<double> GetTransMatrix(int from_chart_id, int to_chart_id, int from_vid) const;
		zjucad::matrix::matrix<double> GetTransMatrix(int from_vid, int from_chart_id, int to_vid, int to_chart_id) const;

		//! get the transition matrix between two ambiguity charts, the ambiguity charts is that two charts share two or more edges
		zjucad::matrix::matrix<double> GetTransMatrixBetweenAmbiguityCharts(int from_vert, int from_chart_id, int to_vert, int to_chart_id) const;
		void TransParamCoordBetweenAmbiguityCharts(int from_vert, int from_chart_id, int to_chart_id, const ParamCoord& from_coord, ParamCoord& to_coord) const;

		void SetTransList(const std::map< std::pair<int, int>, std::vector<int> >& trans_list) { m_edge_trans_list = trans_list; }
        
	public:
		//! get the transition list between two charts
		bool GetTranslistBetweenTwoCharts(int from_chart_id, int to_chart_id, std::vector<int>& trans_list) const;

		//! get the edge transition lists between two vertices, if the two vertices in same chart, then the lists is empty, else find the shortest path between them 
		bool GetPatchEdgeTransLists(int from_vid, int to_vid, std::vector<int>& edge_trans_list) const;

		//! get the matrix between two adjacent charts
		zjucad::matrix::matrix<double> GetTransMatrixOfAdjCharts(int from_chart, int to_chart, int vid) const;
		//! get the transition matrix in one chart from a frame to another frame
		zjucad::matrix::matrix<double> GetTransMatrixInOneChart(int chart_id, std::pair<int,int> old_x_axis, std::pair<int,int> new_x_axis) const;

		//! choose which edge to transition for ambiguity charts 
		int ChooseEdgeForAmbiguityChart( const std::vector<int>& common_edge_index_array, int from_vert, int to_vert=-1) const;

        zjucad::matrix::matrix<double> GetScaleMatrixBetweenTwoCharts(int from_cid, int to_cid, int pedge_index1, int pedge_index2) const;
        
	private:
		boost::shared_ptr<ChartCreator> p_chart_creator;
		std::map< std::pair<int, int>, std::vector<int> > m_edge_trans_list;
    };
}

#endif
