#include "ParameterControl.h"
#include "../Param/Parameter.h"
#include "../Param/ParamDrawer.h"
#include "../ModelMesh/MeshModel.h"
#include <fstream>

ParameterControl::ParameterControl(QGLViewer* _gl_viewer, QWidget* parent)
: QWidget(parent)
{
    m_gl_viewer = _gl_viewer;
    
	m_surface_group = CreateSurfaceGroup(this);
	m_texture_setting_group = CreateTextureGroup(this);
	m_visualization_group = CreateVisualizationGroup(this);
    
	CreateMainLayout();
}

QGroupBox* ParameterControl::CreateSurfaceGroup(QWidget* parent /* = 0 */)
{
	QGroupBox* surface_group = new QGroupBox(parent);
	surface_group->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
	surface_group->setTitle(tr("Surface 1 Setting"));
    
	/// child widgets
	QPushButton* load_surface_mesh = new QPushButton(surface_group);
	load_surface_mesh->setText(tr("Load Surface 1 Mesh"));
	QPushButton* load_surface_patch = new QPushButton(surface_group);
	load_surface_patch->setText(tr("Load Surface 1 Patch"));
// 	QPushButton* optimize_ambiguity_patch = new QPushButton(surface_1_group);
// 	optimize_ambiguity_patch->setText(tr("Optimize Ambiguity Patch"));
	QPushButton* parameter = new QPushButton(surface_group);
	parameter->setText(tr("Parameter"));

	/// layout
	QVBoxLayout* surface_layout = new QVBoxLayout(surface_group);
	surface_layout->addWidget(load_surface_mesh);
	surface_layout->addWidget(load_surface_patch);
//	surface_1_layout->addWidget(optimize_ambiguity_patch);
	surface_layout->addWidget(parameter);

	/// connections
	connect(load_surface_mesh, SIGNAL(clicked()), m_gl_viewer, SLOT(loadMeshModel()));
	connect(load_surface_patch, SIGNAL(clicked()), m_gl_viewer, SLOT(loadQuadFile()));
//	connect(optimize_ambiguity_patch, SIGNAL(clicked()), m_gl_viewer, SLOT(OptimizeAmbiguityPatch()));
	connect(parameter, SIGNAL(clicked()), m_gl_viewer, SLOT(SolveParameter()));

	return surface_group;
}

QGroupBox* ParameterControl::CreateTextureGroup(QWidget* parent /* = 0 */)
{
	QGroupBox* texture_group = new QGroupBox(parent);
	texture_group->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
	texture_group->setTitle(tr("Texture"));

	/// child widgets
	QRadioButton* square_texture = new QRadioButton(texture_group);
	QRadioButton* line_texture = new QRadioButton(texture_group);
	QRadioButton* boundary_texture = new QRadioButton(texture_group);

	square_texture->setText(tr("Square"));
	line_texture->setText(tr("Line"));
	boundary_texture->setText(tr("Boundary"));
	square_texture->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
	line_texture->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
	boundary_texture->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);

	QVBoxLayout* texture_layout = new QVBoxLayout(texture_group);
	texture_layout->addWidget(square_texture);
	texture_layout->addWidget(line_texture);
	texture_layout->addWidget(boundary_texture);

	/// connections
	connect(square_texture, SIGNAL(clicked()), m_gl_viewer, SLOT(CreateSquareTexture()));
	connect(line_texture, SIGNAL(clicked()), m_gl_viewer, SLOT(CreateLineTexture()));
	connect(boundary_texture, SIGNAL(clicked()), m_gl_viewer, SLOT(CreateBoundaryTexture()));

	square_texture->setChecked(true);

	return texture_group;
}

QGroupBox* ParameterControl::CreateVisualizationGroup(QWidget* parent /* = 0 */)
{
	QGroupBox* visual_group = new QGroupBox(parent);
	visual_group->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
	visual_group->setTitle(tr("Visulzation"));

	/// child widgets
	// QRadioButton* patch_layout = new QRadioButton(visual_group);
	QCheckBox* patch_conner = new QCheckBox(visual_group);
	QCheckBox* patch_edge = new QCheckBox(visual_group);
	QCheckBox* patch_face = new QCheckBox(visual_group);
	QCheckBox* outrange_vert = new QCheckBox(visual_group);
	QCheckBox* select_patch = new QCheckBox(visual_group);
	QCheckBox* flipped_triangle = new QCheckBox(visual_group);
    
	patch_conner->setText(tr("Patch Conner"));
	patch_edge ->setText(tr("Patch Edge"));
	patch_face ->setText(tr("Patch Face"));
	outrange_vert->setText(tr("Out Range Vertex"));
	select_patch->setText(tr("Selected Patch"));
	flipped_triangle ->setText("Flipped Triangle");

	patch_conner->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
	patch_edge->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
	patch_face->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
	outrange_vert->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
	select_patch->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
	flipped_triangle->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);

	QVBoxLayout* visual_layout = new QVBoxLayout(visual_group);
	visual_layout->addWidget(patch_conner);
	visual_layout->addWidget(patch_edge);
	visual_layout->addWidget(patch_face);
	visual_layout->addWidget(outrange_vert);
	visual_layout->addWidget(select_patch);
	visual_layout->addWidget(flipped_triangle);

	/// connects
	connect(patch_conner, SIGNAL(toggled(bool)), this, SLOT(SetPatchConnerDisplay(bool)));
    connect(patch_edge, SIGNAL(toggled(bool)), this, SLOT(SetPatchEdgeDisplay(bool)));
    connect(patch_face, SIGNAL(toggled(bool)), this, SLOT(SetPatchFaceDisplay(bool)));
	connect(outrange_vert, SIGNAL(toggled(bool)), this, SLOT(SetOutRangeVertDisplay(bool)));
	connect(select_patch, SIGNAL(toggled(bool)), this, SLOT(SetSelectedPatchDisplay(bool)));
	connect(flipped_triangle, SIGNAL(toggled(bool)), this, SLOT(SetFlippedTriangleDisplay(bool)));

	return visual_group;
}

void ParameterControl::SetPatchConnerDisplay(bool toggled)
{
	if(m_gl_viewer && m_gl_viewer->p_param_drawer) {
		m_gl_viewer->p_param_drawer->SetDrawPatchConner(toggled);
		m_gl_viewer->updateGL();
	}
}

void ParameterControl::SetPatchEdgeDisplay(bool toggled)
{
	if(m_gl_viewer && m_gl_viewer->p_param_drawer){
		m_gl_viewer->p_param_drawer->SetDrawPatchEdge(toggled);
		m_gl_viewer->updateGL();
	}
}

void ParameterControl::SetPatchFaceDisplay(bool toggled)
{
	if(m_gl_viewer && m_gl_viewer->p_param_drawer){
		m_gl_viewer->p_param_drawer->SetDrawPatchFace(toggled);
		m_gl_viewer->updateGL();
	}
}

void ParameterControl::SetOutRangeVertDisplay(bool toggled)
{
	if(m_gl_viewer && m_gl_viewer->p_param_drawer){
		m_gl_viewer->p_param_drawer->SetDrawOutRangeVertices(toggled);
		m_gl_viewer->updateGL();
	}
}

void ParameterControl::SetSelectedPatchDisplay(bool toggled)
{
	if(m_gl_viewer && m_gl_viewer->p_param_drawer){
		m_gl_viewer->p_param_drawer->SetDrawSelectedPatch(toggled);
		m_gl_viewer->updateGL();
	}
}

void ParameterControl::SetFlippedTriangleDisplay(bool toggled)
{
	if(m_gl_viewer && m_gl_viewer->p_param_drawer){
		m_gl_viewer->p_param_drawer->SetDrawFlipFace(toggled);
		m_gl_viewer->updateGL();
	}
}

void ParameterControl::CreateMainLayout()
{
   QGroupBox* main_group = new QGroupBox(this);
   main_group->setFixedWidth(200);
   main_group->setTitle(tr("Cross Parameterization"));

   QVBoxLayout* main_group_layout = new QVBoxLayout(main_group);
   main_group_layout->setMargin(3);
   main_group_layout->setSpacing(30);
   main_group_layout->addWidget(m_surface_group);
   main_group_layout->addWidget(m_texture_setting_group);
   main_group_layout->addWidget(m_visualization_group);
   main_group_layout->addStretch(1);

}
