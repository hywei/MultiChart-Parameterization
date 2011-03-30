#ifndef TRANSFUNCTION_H_
#define TRANSFUNCTION_H_

#include "../hj_3rd/include/zjucad/matrix/matrix.h"

namespace PARAM
{
	class TransFunc
	{
	public:
		TransFunc();
		TransFunc(const TransFunc&);
		~TransFunc();

	public:
		virtual zjucad::matrix::matrix<double> GetTransMatrix(int from_chart_id, int to_chart_id);

	protected:
		std::vector< std::vector<int> >  

	};
}

#endif