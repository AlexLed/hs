#ifndef MAINWINDOW_H
#define MAINWINDOW_H
//#include <QtCore>
#include <QtGui>
#include <QMainWindow>
#include <QtWidgets>
#include <QtCharts>
#include "solver.h"
#include "datadef.h"
#include "graphicswidget.h"
#include "colormap.h"
#include "viewarray.h"
#include "glwidget.h"

using namespace QtCharts;

QT_BEGIN_NAMESPACE
class RenderArea;
class ViewArray;
class WingTab;
class GraphicTab;
class ColorMapTab;
class GLTab;
class WingView;
class ColorMapWidget;
class GLWidget;
void updateWing();
QT_END_NAMESPACE

class MainWindow : public QMainWindow {
   Q_OBJECT
public:
   QString type;
   QTime timeFunc, timeTotal;
   MainWindow();
private:
   int i, j, k, NCases;
   void createTabWidget();
   void createActions();
   void createMenus();
   void createToolBars();
   void createStatusBar();
   void createParameters();
   void fillForm();
   void createResults();
   void createDockWindows();
   void readSettings();
   void writeSettings();
   void loadFile(const QString &fileName);
   bool saveFile(const QString &fileName);
   void saveResults();
   void enableUI(bool is_enabled);
   QTimer *timer;
   QTranslator translator;
   QWidget *mainwidget, *parameters, *results, *colormap;
   QTabWidget *tabWidget;
   WingTab *wingTab;
   GraphicTab *graphicTab;
   ColorMapTab *colorMapTab;
   GLTab *glTab;
   QLineEdit *lineEdit_geometryTheta0d, *lineEdit_geometryTheta1d, *lineEdit_geometryBetad, *lineEdit_geometryH, *lineEdit_geometryL,
             *lineEdit_parametersGamma, *lineEdit_parametersPr, *lineEdit_parametersHw, *lineEdit_parametersOmega,
             *lineEdit_solving_relU, *lineEdit_solving_relP, *lineEdit_solving3, *lineEdit_solving4,
             *lineEdit_solving_AcU, *lineEdit_solving_AcPr,
             *lineEdit_grid_dx, *lineEdit_grid_NX,
             *lineEdit_grid_dy, *lineEdit_grid_NY,
             *lineEdit_grid_dz, *lineEdit_grid_NZ;
   QLabel *label_results_Edge1_1, *label_results_Edge1_2, *label_results_Edge1_3,
          *label_results_Edge2_1, *label_results_Edge2_2, *label_results_Edge2_3,
          *label_results_Nose_1, *label_results_Nose_2, *label_results_Nose_3,
          *label_results_Wing_1, *label_results_Wing_2, *label_results_Wing_3,
          *label_results_TotalTime1, *label_results_FinishTime0, *label_results_FinishTime1;
   QAction *openAct, *saveAct, *saveAsAct, *exitAct, *solveAct, *stopAct, *solveCycleAct, *selfTestAct,
           *viewArrayAct, *colormapAct, *wingAct, *graphAct, *langEnAct, *langRuAct, *helpAct, *aboutQtAct, *aboutAct;
   QMenu *viewMenu;
   QCheckBox *checkBox_graph, *checkBox_MultiGrid;
   QRadioButton *radioC1, *radioC2;
private slots:
   void open();
   bool save();
   bool saveAs();
   void about();
   void aboutQt();
   void changeLang(QAction *action) {
       QSettings settings("settings.ini", QSettings::IniFormat);
       if (action == langEnAct) settings.setValue("Language", "en");
       else if(action == langRuAct)  settings.setValue("Language", "ru");
   }
public slots:
   void getParameters();
   void readParameters(QString line);
   void solve();
   void solveCycle();
   void solveProblem(QString filename);
   void abort() { converged = true; }
   void selfTest();
   void changeStopAct(bool trigger) {stopAct->setDisabled(trigger);}
   void changeCoordinates(bool checked) { cartesian = !checked; updateWing();}
   void changeMethod(bool checked) { thomas_method = checked; }
   void changeOrderTheta(bool checked) { is_first_order_theta = checked; }
   void changeRel(bool checked) { is_mod_rel = checked; }
   void resetProgress();
   void updateProgress();
   void finishProgress();
   void updateWing();
   void showWing() { tabWidget->setCurrentIndex(0); }
   void showGraphic() { tabWidget->setCurrentIndex(1); }
   void showArray() { tabWidget->setCurrentIndex(2); }
   void showColorMap() { tabWidget->setCurrentIndex(3); }
   void changed_dx(QString str) { this->lineEdit_grid_NX->setText(QString::number(L/str.toDouble()+1)); updateWing(); }
   void changed_NX(QString str) { this->lineEdit_grid_dx->setText(QString::number(L/(str.toInt()-1))); updateWing(); }
   void changed_dy(QString str) { this->lineEdit_grid_NY->setText(QString::number(h/str.toDouble()+1)); updateWing(); }
   void changed_NY(QString str) { this->lineEdit_grid_dy->setText(QString::number(h/(str.toInt()-1))); updateWing(); }
   void changed_dz(QString str) { this->lineEdit_grid_NZ->setText(QString::number(2./str.toDouble()+1)); updateWing(); }
   void changed_NZ(QString str) { this->lineEdit_grid_dz->setText(QString::number(2./(str.toInt()-1))); updateWing(); }
   void changed_h(QString str) {
       this->lineEdit_grid_NY->setText(QString::number(str.toDouble()/lineEdit_grid_dy->text().toDouble()+1));
   }

protected:
   void closeEvent(QCloseEvent * /*event*/) {
       converged = true;
       writeSettings();
   }
};

class WingView : public QWidget
{
   Q_OBJECT
public:
   WingView(QWidget *parent = 0);
private:
   QPixmap pixmap;
public slots:
protected:
   void paintEvent(QPaintEvent *event);
   void resizeEvent(QResizeEvent * /*event*/) { this->update(); }
};

class WingTab : public QWidget
{
    Q_OBJECT
public:
    WingTab(QWidget *parent = 0);
    WingView *wingWidget;
    QCheckBox *gridCheckBox;
public slots:
    void drawGrid(bool checked) {draw_grid = checked; wingWidget->update(); }
};
class GraphicTab : public QWidget
{
    Q_OBJECT
public:
    GraphicTab(QWidget *parent = 0);
    QChart *chart;
    QChartView *chartView;
    QValueAxis *axisX;
    QValueAxis *axisY;
    RenderArea *renderarea;
    QSpinBox *spinBoxLX, *spinBoxLY, *spinBoxLZ;
    QDoubleSpinBox *spinBoxMinX, *spinBoxMaxX, *spinBoxMinY, *spinBoxMaxY;
    QComboBox *fileComboBox;
public slots:
    void Fx() {axisX->setTickCount(2);}
    void loadData();
    void changeLX(int i) { LX = i; renderarea->update(); }
    void changeLY(int i) { LY = i; renderarea->update(); }
    void changeLZ(int i) { LZ = i; renderarea->update(); }
    void changeMinX(double i) { renderarea->minX = i; renderarea->update(); }
    void changeMaxX(double i) { renderarea->maxX = i; renderarea->update(); }
    void changeMinY(double i) { renderarea->minY = i; renderarea->update(); }
    void changeMaxY(double i) { renderarea->maxY = i; renderarea->update(); }
};

class ArrayTab : public QWidget
{
    Q_OBJECT
public:
    ArrayTab(QWidget *parent = 0);
private slots:
    void changeLX(int i) { LX = i; }
};

class ColorMapTab : public QWidget
{
    Q_OBJECT
public:
    ColorMapTab(QWidget *parent = 0);
    ColorMapWidget *colorMap;
private:
    QCheckBox *chBox_flowlines, *chBox_vector_plot;
private slots:
    void changeFlowlines(bool checked) {
        colorMap->flowlines = checked;
        colorMap->update();
    }
    void drawTauU() {
        colorMap->setupPlot(TauUDistribution, "TauU", chBox_vector_plot->isChecked(), chBox_flowlines->isChecked());
        colorMap->update();	}
    void drawTauW() {
        colorMap->setupPlot(TauWDistribution, "TauW", chBox_vector_plot->isChecked(), chBox_flowlines->isChecked());
        colorMap->update();	}
    void drawTauH() {
        colorMap->setupPlot(TauHDistribution, "TauH", chBox_vector_plot->isChecked(), chBox_flowlines->isChecked());
        colorMap->update();	}
    void drawPressure() {
        colorMap->setupPlot(PressureDistribution, "Pressure", chBox_vector_plot->isChecked(), chBox_flowlines->isChecked());
        colorMap->update();	}
    void drawDelta() {
        colorMap->setupPlot(DeltaDistribution, "Delta", chBox_vector_plot->isChecked(), chBox_flowlines->isChecked());
        colorMap->update(); }
};

class GLTab : public QWidget
{
    Q_OBJECT
public:
    GLTab (QWidget *parent = 0);
    GLWidget *glWidget;
};

#endif
