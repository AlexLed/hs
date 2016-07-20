#ifndef RENDERAREA_H
#define RENDERAREA_H
#include <QtGui>
#include <QWidget>
#include <QToolButton>
#include "solver.h"

class RenderArea : public QWidget
{
    Q_OBJECT
public:
    RenderArea(QWidget *parent = 0);
    void setupAxis(double min_x, double max_x, double min_y, double max_y, double step_x, double step_y, QString label_x, QString label_y);

    QToolButton *zoomInButton, *zoomOutButton;
    //QAction *exportPlotAsPictureAct;
    double minX, maxX, minY, maxY;
    struct Curve {
        double ***data;
        QPainterPath plot;
        QString caption;
        QPen pen;
        bool show;
    };
    struct DataFromFile {
        double ***array;
        QStringList title, caption;
    };
    int Narr;
    Curve curve[10];
    DataFromFile data;
    int i, j, k, n, Nmax, NXd, NYd, NZd, Ndata, Ncalc;
private:
    double sx, sy, stepX, stepY;
    QString labelX, labelY;
    //double m1, m2, c1, c2, M1, M2, C1, C2;
    int Margin; //offset of plot from boundaries
    void paintAxis(QPainter &painter);
    void plotGraph(QPainter &painter);
    void createActions();
    //void contextMenuEvent(QContextMenuEvent *event);

private slots:
    //void exportPlotAsPicture();

public slots:
   void Fx() {
       setupAxis(0, L, 0, 1, 0.2, 0.2, "r", "F(x)");
       this->update();
   }
   void Fy() {
       setupAxis(-1, 2, 0, dy*jmax, 0.5, 1, "F(y)", "y");
       this->update();
   }
   void Fz() {
       setupAxis(-1, 1, -0.4, 1, 0.2, 0.2, tr("%1").arg(QChar(0x03b8)), "F(z)");
       this->update();
   }
   void PaintTau() {
       setupAxis(-1, 1, 0.0, 1, 0.2, 0.2, tr("%1").arg(QChar(0x03b8)), "Tau");
       this->update();
   }
   void PaintDelta() {
       setupAxis(-1, 1, 0.0, 2, 0.2, 0.5, tr("%1").arg(QChar(0x03b8)), "D");
       this->update();
   }
   void PaintP() {
       setupAxis(-1, 1, -0.5, 2, 0.2, 0.5, tr("%1").arg(QChar(0x03b8)), "P");
       this->update();
   }
   void PaintDisturbances() {
       setupAxis(-1, 1, -1, 1, 0.1, 0.1, tr("%1").arg(QChar(0x03b8)), "Disturbances");
       this->update();
   }
   void PaintData(int index) {
       Narr = index;
       setupAxis(-1, 1, 0, 1, 0.2, 0.2, tr("%1").arg(QChar(0x03b8)), data.caption.at(index));
       this->update();
   }

protected:
   void paintEvent(QPaintEvent *event);
   void resizeEvent(QResizeEvent *event);
//	void mousePressEvent(QMouseEvent *event);
//	void mouseMoveEvent(QMouseEvent *event);
//	void mouseReleaseEvent(QMouseEvent *event);
//	void keyPressEvent(QKeyEvent *event);
//	void wheelEvent(QWheelEvent *event);
};

#endif
