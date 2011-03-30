#include "MainWindow.h"
#include <QMainWindow>
#include <QtOpenGL>
#include <iostream>


int main(int argc, char *argv[])
{    
	QApplication::setColorSpec( QApplication::CustomColor );
	QApplication app(argc,argv);

	if ( !QGLFormat::hasOpenGL() ) {
		QString msg = "System has no OpenGL support!";
		QMessageBox::critical( 0, QString("OpenGL"), msg + QString(argv[1]) );
		return -1;
	}
	// create widget
    MainWindow mainWin;

    mainWin.show();

	return app.exec();
}
