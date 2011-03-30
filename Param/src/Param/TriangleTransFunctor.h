#ifndef TRIANGLETRANSFUNCTOR_H_
#define TRIANGLETRANSFUNCTOR_H_

#include "../ModelMesh/MeshModel.h"
#include <hj_3rd/zjucad/matrix/matrix.h>

namespace PARAM
{
	/*
	 */
	class TriTransFunctor
	{
	public:
		TriTransFunctor(boost::shared_ptr<MeshModel> _mesh);
		~TriTransFunctor();

	public:
		bool TransLocalCoordBetweenCharts(int from_chart_id, const Coord2D& from_chart_coord, int to_chart_id, Coord2D& to_chart_coord) const;

	//protected:
		bool GetTransChartsChain(int from_chart_id, int to_chart_id, std::vector<int>& ) const;
		zjucad::matrix::matrix<double> GetTransMatrix(int from_chart_id, int to_chart_id) const;

		//! Get the transition matrix between two adjacent charts
		zjucad::matrix::matrix<double> GetTransMatrixOfAdjCharts(int from_chart_id, int to_chart_id) const;
		//! Get the transition matrix in the same chart, just changed the origin and x-axis
		zjucad::matrix::matrix<double> GetTransMatrixInOneChart(int chart_id, 
			std::pair<int, int> old_x_index,  std::pair<int, int> new_x_index) const;

		//! some matrix operator
		zjucad::matrix::matrix<double> GetTranslationMatrix2D(double t_x, double t_y) const;
		zjucad::matrix::matrix<double> GetRotationMatrix2D(double rotation_angle) const;

		void SetChartNeighbors();

	private:
		boost::shared_ptr<MeshModel> p_mesh;
		std::vector< std::vector<int> > m_neighbor_charts;

	public:
		/// Test 
		void TestInnerChartTrans(int chart_id) const;
		void TestAdjChartsTrans(int from_chart_id, int to_chart_id) const;
		void TestTwoChartsTrans(int from_Chart_id, int to_chart_id) const;

	};
}

#endif
