#include "mainWindow.h"
#include "MainLogic.h"
#include <QtGui/QApplication>
#include <iostream>
#include <fstream>
#include <itkFileOutputWindow.h>
#include <vtkFileOutputWindow.h>
#include <string>
#include <QDir>
#include <QFileInfo>

#ifdef _WINDOWS //for chdir
#include <direct.h>
#else
#include <unistd.h>
#endif

using namespace std;

int main( int argc, char **argv )
{
	QApplication app(argc, argv);
	app.setOrganizationName("CG Siegen");
    app.setOrganizationDomain("www.cg.informatik.uni-siegen.de");
    app.setApplicationName("Spine Analyzer");

	#ifdef _WINDOWS
	//redirect output to files
	ofstream fout("D:\\cout.txt");
	ofstream ferr("D:\\cerr.txt");
	cout.rdbuf(fout.rdbuf());
	cerr.rdbuf(ferr.rdbuf());
	#endif

	MainWindow mainForm;
	mainForm.show();

    //itk::FileOutputWindow::Pointer win=itk::FileOutputWindow::New();
    //win->SetFileName(mainForm.logic->tempDir+"itkOutput.txt");
    //win->SetInstance(win);
    //vtkFileOutputWindow *winVtk=vtkFileOutputWindow::New();
    //winVtk->SetFileName((mainForm.logic->tempDir+"vtkOutput.txt").c_str());
    //winVtk->SetInstance(winVtk);
    
    QDir qd;
    if (!qd.exists("cascade.xml")) //the working directory does not contain support files
    {
        QFileInfo ed(argv[0]);
        chdir(ed.absoluteDir().canonicalPath().toStdString().c_str());
    }
    
    if (argc==1)
        return app.exec();
    //else segment and quit
    mainForm.actionInteractive_inflation->setChecked(false); //do not show inflation steps
    if (argc>=3)
        mainForm.logic->tempDir=argv[2];
    try
    {
        if (argc>=4 && strcmp(argv[3],"fullAuto")==0)
            mainForm.logic->fullAuto=true;
        QFileInfo qf(argv[1]);
        string fn(argv[1]);
        int pos=fn.find_last_of('.');
        string ext=fn.substr(pos+1);
        
        if (qf.isDir())
            mainForm.logic->dicomOpen(fn);
        else if (ext=="init")
            mainForm.openInitialization(argv[1]);
        else
        {
            mainForm.logic->volumeOpen(fn);
            return app.exec();
        }   
    }
    catch (...)
    {
        return 1;
    }
    return 0;
}

#ifdef _WINDOWS
#pragma comment(linker, "/subsystem:\"windows\" /entry:\"mainCRTStartup\"") //disables cmd window in VS2008
#endif
