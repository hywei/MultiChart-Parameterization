#ifndef PARAMETERCONTROL_H_
#define PARAMETERCONTROL_H_

#include "QGLViewer.h"

namespace Param{
    class Parameter;
}

class ParameterControl : public QWidget
{
	Q_OBJECT

public:
	ParameterControl(QGLViewer* _gl_viewer, QWidget* parent = 0);

private:
	QGroupBox* CreateSurfaceGroup(QWidget* parent = 0);
	QGroupBox* CreateTextureGroup(QWidget* parent = 0);
	QGroupBox* CreateVisualizationGroup(QWidget* parent = 0);
    QGroupBox* CreateChartOptimizationGroup(QWidget* parent=0);
    
	void CreateMainLayout();

private slots:
	void SetPatchConnerDisplay(bool );
	void SetPatchEdgeDisplay(bool );
	void SetPatchFaceDisplay(bool );
	void SetOutRangeVertDisplay(bool );
	void SetSelectedPatchDisplay(bool);
	void SetFlippedTriangleDisplay(bool);

    void ChartOptimization(); 
private:
	QGLViewer* m_gl_viewer;

	QGroupBox* m_surface_group;
	QGroupBox* m_texture_setting_group;
	QGroupBox* m_visualization_group;
    QGroupBox* m_chart_optimization_group;

    boost::shared_ptr<PARAM::Parameter> p_parameter;

    QLineEdit* m_chart_init_value_le;
};

#endif
