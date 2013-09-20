#include <QtGui/QApplication>
#include "mainWindow.h"
#include <string>
#include <fstream>
#include <iostream>
#include <QFile>
using namespace std;

const string labels[]=
    {"S5", "S4", "S3", "S2", "S1", "L5", "L4", "L3", "L2", "L1",
    "T12", "T11", "T10", "T9", "T8", "T7", "T6", "T5", "T4", "T3", "T2", "T1",
    "C7", "C6", "C5", "C4", "C3", "C2", "C1"};

extern Slices openHelper(string filename, QColor color, bool sagOnly=false);

string base, folderManual, folderAuto, filename;
float dsc1(MainWindow& mf, int index)
{
    string as=folderAuto+filename+'_'+labels[index]+".mha";
    QFile qf(as.c_str());
    if (!qf.exists())
        return 0;

    mf.slices=openHelper(base, Qt::white);
    mf.viewDirRadioButtonToggled(true);
	Slices seg=openHelper(folderManual+labels[index]+".mha", c1, true);
	mf.blend(seg);
	seg=openHelper(as, c2, true);
	mf.blend(seg);
    return mf.DiceSimilarityCoefficient();
}

int main( int argc, char **argv )
{
	QApplication app(argc, argv);
	MainWindow mainForm;
	mainForm.show();
    if (argc==7)
    {
        base=string(argv[1]); //original image
        folderManual=string(argv[2]);
        folderAuto=string(argv[3]);
        int slash=base.find_last_of('/');
        int dot=base.find_last_of('.');
        filename=base.substr(slash+1,dot-slash-1);
        int ind1=0, ind2=0;
        string v1(argv[4]); //start vertebra
        string v2(argv[5]); //end vertebra
        while (labels[ind1]!=v1)
            ind1++;
        while (labels[ind2]!=v2)
            ind2++;
        if (ind1>ind2)
            swap(ind1, ind2);

        ofstream f(folderAuto+argv[6]);
        f<<"Vertebra,DSC\n";
        for (int i=ind2; i>=ind1; i--)
        {
            f<<labels[i]<<','<<dsc1(mainForm,i)<<endl;
        }
        f.close();
    }
    else if (argc>1)
    {
	    mainForm.slices=openHelper(argv[1], Qt::white);
	    mainForm.setWindowTitle(QString("Segmentation Comparer - ")+QString::fromStdString(argv[1]));
	    mainForm.viewDirRadioButtonToggled(true);
        return app.exec();
    }
    else
        return app.exec();    
}

#ifdef _WINDOWS
#pragma comment(linker, "/subsystem:\"windows\" /entry:\"mainCRTStartup\"") //disables cmd window
#endif
