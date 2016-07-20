#ifndef COLORMAP_H
#define COLORMAP_H
#include <QWidget>
#include <QDebug>

class ColorMapWidget : public QWidget
{
    Q_OBJECT
public:
    bool vector_plot, flowlines;
    ColorMapWidget(QWidget *parent);
    double **array_for_drawing;
    QString caption;
private:
    int i, k, n;
    void drawVectorPlot(QPainter *painter);
    void drawFlowLines(QPainter *painter);
public slots:
    void setupPlot(double **array, QString string, bool vp, bool fl){
        array_for_drawing = array;
        caption = string;
        vector_plot = vp;
        flowlines = fl;
    }
protected:
    void paintEvent(QPaintEvent *event);
};

#endif
