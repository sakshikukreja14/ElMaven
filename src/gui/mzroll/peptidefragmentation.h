#ifndef PEPTIDEFRAGMENTATION_H
#define PEPTIDEFRAGMENTATION_H

#include "stable.h"
#include "ui_peptidefragmentation.h"

class MainWindow;
class MassCutoff;
class Scan;

class PeptideFragmentationWidget: public QDockWidget,  public Ui_PaptideFragmentationWidget {
      Q_OBJECT

public:
      PeptideFragmentationWidget(MainWindow* mw);

protected:
        

public Q_SLOTS: 
 	  void setPeptideSequence(QString sequence);
	  void setCharge(float charge);
	  void setResolution(MassCutoff *massCutoff);
      void setScan(Scan*);

private Q_SLOTS:
	  void compute();
      void showTable();
      void focusPrecursorMz();
   
private:
      MainWindow* _mw;
      QString _sequence;
	  double _charge;
      MassCutoff *_massCutoff;
      Scan* _scan;
      
};

#endif
