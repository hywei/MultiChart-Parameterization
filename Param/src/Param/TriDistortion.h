#ifndef TRIDISTORTION_H_
#define TRIDISTORTION_H_

#include <hj_3rd/zjucad/matrix/matrix.h>
#include <vector>

namespace PARAM
{
	class Parameter;

	class TriDistortion
	{
	public:
		TriDistortion(const Parameter& parameter);
		~TriDistortion();

		void ComputeDistortion();

	public:
		//! IO 
		const std::vector<double>& GetFaceHarmonicDistortion() const { return m_face_harmonic_distortion; }
		const std::vector<double>& GetFaceIsometricDistortion() const{ return m_face_isometric_distortion; }

		//! compute the jacobi matrix from surface to parameter domain
		zjucad::matrix::matrix<double> ComputeParamJacobiMatrix(int fid) const;


	private:
		double ComputeHarmonicDistortion(const zjucad::matrix::matrix<double>& jacobi_mat) const;
		double ComputeIsometricDistortion(const zjucad::matrix::matrix<double>& jacobi_mat) const;		

	private:
		const Parameter& m_parameter;
		std::vector<double> m_face_harmonic_distortion; //! each face's harmonic map distortion
		std::vector<double> m_face_isometric_distortion; //! each face's isometric map distortion

	};
}

#endif //TRIDISTORTION