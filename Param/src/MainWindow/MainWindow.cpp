#include "MainWindow.h"

#include <QtGui>
#include <sstream>
#include <iostream>

#include "ParameterControl.h"
#include "QGLViewer.h"

MainWindow::MainWindow(QWidget *parent):QMainWindow(parent)
{

	setMinimumSize(20, 20);
	resize(1500, 900); 

	glViewer = new QGLViewer(this);

	m_parameter_control = new ParameterControl(glViewer, this);	

	createActions();
	createMenus();
	createToolBars();

	QList <int> sizes;
	sizes<<1200<<300;
	QSplitter* mainSplitter = new QSplitter(Qt::Horizontal);
	mainSplitter->addWidget(glViewer);
	mainSplitter->addWidget(m_parameter_control);
    mainSplitter->setSizes(sizes);
	mainSplitter->setStretchFactor(0, 1);

	setCentralWidget(mainSplitter);

	setWindowTitle("Parameter");
}

void MainWindow::setStyle(QStyle* style)
{
	QWidget::setStyle(style);
	if(m_parameter_control != 0)
	{
		m_parameter_control->setStyle(style);

		QList<QWidget*> widgets = qFindChildren<QWidget*> (m_parameter_control);
		foreach(QWidget* w, widgets)
			w->setStyle(style);
	}
}

void MainWindow::mouseSpin()
{
	glViewer->mouseSpin();
	glViewer->setCursor(Qt::PointingHandCursor);
}

void MainWindow::mouseMove()
{
	glViewer->mouseMove();
	glViewer->setCursor(Qt::ArrowCursor);
}

void MainWindow::mouseZoom()
{
	glViewer->mouseZoom();
	glViewer->setCursor(Qt::SizeAllCursor);
}

void MainWindow::createActions()
{
	// file menu actions
	openModelAct = new QAction(QIcon("../images/open.png"), tr("&Open Model"), this);
	openModelAct->setShortcut(QKeySequence::New);
	openModelAct->setStatusTip(tr("open a mesh model"));
	connect(openModelAct, SIGNAL(triggered()), this, SLOT(openModel()));

	openTextureFileAct = new QAction(QIcon("../images/open.png"), tr("&Open Texture Image"), this);
	connect(openTextureFileAct, SIGNAL(triggered()), this, SLOT(openTextureImage()));

	openQuadFileAct = new QAction(QIcon("../images/open.png"), tr("&Open Quad File"), this);
	connect(openQuadFileAct, SIGNAL(triggered()), this, SLOT(openQuadFile()));

	saveModelAct = new QAction(QIcon("../images/save.png"), tr("&Save Model"), this);
	saveModelAct->setShortcut(QKeySequence::Save);
	saveModelAct->setStatusTip(tr("save this mesh model"));
	connect(saveModelAct, SIGNAL(triggered()), this, SLOT(saveModel()));

	saveAsBmpAct = new QAction(tr("Save as BMP"), this);
	connect(saveAsBmpAct, SIGNAL(triggered()), this, SLOT(saveAsBmp()));

	recentFileAct = new QAction(tr("Recent Files"), this);
	connect(recentFileAct, SIGNAL(triggered()), this, SLOT(recentFiles()));

	exitAct = new QAction(tr("E&xit"), this);
	exitAct->setShortcuts(QKeySequence::Quit);
	exitAct->setStatusTip(tr("Exit the application"));
	connect(exitAct, SIGNAL(triggered()), this, SLOT(close()));


	// toolbar menu actions
	toolBarAct = new QAction(tr("ToolBar"), this);
	connect(toolBarAct, SIGNAL(triggered()), this, SLOT(toobBar()));

	stateBarAct = new QAction(tr("StateBar"), this);
	connect(stateBarAct, SIGNAL(triggered()), this, SLOT(stateBar()));

	// quad menu actions
	
	// help menu actions
	aboutAct = new QAction(tr("About"), this);
	connect(aboutAct, SIGNAL(triggered()), this, SLOT(about()));



	// mouse actions
	mouseSpinAct = new QAction(QIcon("../images/rotate-left.png"), tr("Spin"), this);
	connect(mouseSpinAct, SIGNAL(triggered()), this, SLOT(mouseSpin()));

	mouseMoveAct = new QAction(QIcon("../images/move.png"), tr("Move"), this);
	connect(mouseMoveAct, SIGNAL(triggered()), this, SLOT(mouseMove()));

	mouseZoomAct = new QAction(QIcon("../images/zoom.png"), tr("Zoom"), this);
	connect(mouseZoomAct, SIGNAL(triggered()), this, SLOT(mouseZoom()));

}


void MainWindow::createMenus()
{
	// create file menu
	fileMenu = menuBar()->addMenu(tr("&File"));
	fileMenu->addAction(openModelAct);
	fileMenu->addAction(openTextureFileAct);
	fileMenu->addAction(openQuadFileAct);
	fileMenu->addAction(saveModelAct);
	fileMenu->addSeparator();
	fileMenu->addAction(saveAsBmpAct);
	fileMenu->addAction(recentFileAct);
	fileMenu->addSeparator();
	fileMenu->addAction(exitAct);

	// create view menu
	viewMenu = menuBar()->addMenu(tr("&View"));
	viewMenu->addAction(toolBarAct);
	viewMenu->addAction(stateBarAct);
    
	// create help menu
	helpMenu = menuBar()->addMenu(tr("&Help"));
	helpMenu->addAction(aboutAct);
}

void MainWindow::createToolBars()
{
	fileToolBar = addToolBar(tr("File"));
	fileToolBar->addAction(openModelAct);
	fileToolBar->addAction(saveModelAct);

	mouseActToolBar = addToolBar(tr("Mouse Action"));
	mouseActToolBar->addAction(mouseSpinAct);
	mouseActToolBar->addAction(mouseMoveAct);
	mouseActToolBar->addAction(mouseZoomAct);
}


void MainWindow::openModel()
{
	glViewer->loadMeshModel();
}

void MainWindow::openTextureImage()
{
	glViewer->loadTextureImage();
}
void MainWindow::openQuadFile()
{
	glViewer->loadQuadFile();
}
void MainWindow::saveModel()
{
	glViewer->saveMeshModel();
}
void MainWindow::saveAsBmp(){}
void MainWindow::recentFiles(){}
void MainWindow::toobBar(){}
void MainWindow::stateBar(){}
void MainWindow::about(){}
MainWindow::~MainWindow(){}
