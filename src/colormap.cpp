#include <QtGui>
#include "colormap.h"
#include "math.h"
#include "solver.h"

ColorMapWidget::ColorMapWidget(QWidget *parent)
    : QWidget(parent)
{
    setAutoFillBackground(false);
    flowlines = false;
    setupPlot(PressureDistribution, "Pressure", false, false);
}
void ColorMapWidget::paintEvent(QPaintEvent *event)
{
    QPainter painter;
    painter.begin(this);
    painter.fillRect(event->rect(), QBrush(QColor(255, 255, 255)));
    double a, t, min = 1000, max = 0;
    long int f;
    for(i = 0; i < N; i++)
        for(k = 0; k < N; k++)
        {
            t = array_for_drawing[i][k]; //for speed
            if ( t != 0 )
            {
                if( max < t ) max = t;
                if( min > t ) min = t;
            }
        }
    a = 1023/(max - min);
    painter.translate(5,30);
    for(int iy = 0; iy < N; iy++)
        for(int ix = 0; ix < N; ix++)
        {
            f = int(a*(array_for_drawing[ix][iy] - min));
            if(array_for_drawing[ix][iy] == 0) painter.setPen(Qt::white);
            else {
                if(f < 0 || f > 1023) painter.setPen(QColor(0, 0, 0));
                else if(f <= 255) painter.setPen(QColor(0, f, 255));
                else if(f <= 511) painter.setPen(QColor(0, 255, 511-f));
                else if(f <= 767) painter.setPen(QColor(f-512, 255, 0));
                else painter.setPen(QColor(255, 1023-f, 0));
            }
            painter.drawPoint(ix,iy);
        }
    //Gradient scale
    int rect_w = 20, rect_l = 150;
    QRect gradient_rect(N+15, 20, rect_w, rect_l);
    painter.setRenderHint(QPainter::Antialiasing);
    QLinearGradient gradient(gradient_rect.topLeft(), gradient_rect.bottomRight());
    gradient.setColorAt(1.00, QColor(0,	0, 255));
    gradient.setColorAt(0.75, QColor(0, 255, 255));
    gradient.setColorAt(0.50, QColor(0, 255, 0));
    gradient.setColorAt(0.25, QColor(255, 255, 0));
    gradient.setColorAt(0.00, QColor(255, 0, 0));
    painter.fillRect(gradient_rect, gradient);

    painter.setPen(QPen(Qt::black, 2));
    painter.setFont(QFont("Arial", 9));
    double nticks = 4;
    for(i = 0; i < nticks; i++)
    {
        QPoint line_begin, line_end;
        line_begin = QPoint(gradient_rect.right()+2, gradient_rect.bottom()-i*(gradient_rect.height()-3)/(nticks-1));
        line_end = line_begin + QPoint(1,0);
        painter.drawLine(line_begin, line_end);
        painter.drawText(QRect(line_end + QPoint(3, -8), line_end + QPoint(35 , 4)),
                         Qt::AlignLeft,
                         QString::number(max*i/(nticks-1)+min*(nticks-1-i)/(nticks-1),'f',2));
    }
    painter.setFont(QFont("Helvetica", 15/*, QFont::Bold*/));
    painter.drawText(QRect(0, -30, N, 30), Qt::AlignCenter, caption);
    painter.setPen(QPen(Qt::black, 0.1, Qt::SolidLine, Qt::SquareCap, Qt::MiterJoin));
    painter.setRenderHint(QPainter::Antialiasing, false);
    if (vector_plot) drawVectorPlot(&painter);
    if ( flowlines ) drawFlowLines(&painter);
    painter.end();
}
void ColorMapWidget::drawFlowLines(QPainter *painter) {
    double phi, x, x2, y, z_e;
    QPainterPath flowline_l, flowline_r;
    int numlines = 19;
    double space=5;
    if (!cartesian) painter->translate(N*0.5,N*0.4);
    for (int nl = 0; nl < numlines; nl++) {
        if(cartesian) {
            space = 30;
            x = N/2+Ze[nl*2]/z0*N;
            x2 = N/2-Ze[nl*2]/z0*N;
            y = X[nl*2]*N;
            flowline_l.moveTo(x2, y);
            flowline_r.moveTo(x, y);
            for (int j = 0; j < 10000 && !((y < 0 || y > N || x < 0 || x > N)); j++) {
                i = int(y*(imax-1)*2*z0/N);
                z_e = Ze[i+1]/2/z0*N;
                k = int((x-N*0.5)/z_e*(kmax-1)*0.5+(kmax-1)*0.5);
                x += W[i][1][k]/dy;
                y += U[i][1][k]/dy;
                k = int((x2-N*0.5)/z_e*(kmax-1)*0.5+(kmax-1)*0.5);
                x2 += W[i][1][k]/dy;
                y += U[i][1][k]/dy;
                flowline_l.lineTo(x2, y);
                flowline_r.lineTo(x, y);
            }
        }
        else {
            x = space*sin(-Theta0-Beta);
            y = space*cos(-Theta0-Beta);
            flowline_l.moveTo(x, y);
            for (k = 0; k < 10000; k++) {
                phi = atan2(x,y)/Theta0 + Beta/Theta0;
                n = (int)floor(phi/dz+0.5) + kmax/2;
                x += (U[0][1][n]*sin(Theta0*Z[n]-Beta)+W[0][1][n]*cos(Theta0*Z[n]-Beta))/dy/4;
                y += (U[0][1][n]*cos(Theta0*Z[n]-Beta)-W[0][1][n]*sin(Theta0*Z[n]-Beta))/dy/4;
                if((y > -0.4*N) && (y < 0.6*N) && (x > -0.5*N) && (x < 0.5*N))
                    flowline_l.lineTo(x, y);
                else flowline_l.moveTo(x, y);
            };

            x = space*sin(Theta0-Beta);
            y = space*cos(Theta0-Beta);
            flowline_r.moveTo(x, y);
            for (k = 0; k < 10000; k++) {
                phi = atan2(x,y)/Theta0 + Beta/Theta0;
                n = (int)floor(phi/dz+0.5) + kmax/2;
                x += (U[0][1][n]*sin(Theta0*Z[n]-Beta)+W[0][1][n]*cos(Theta0*Z[n]-Beta))/dy/4;
                y += (U[0][1][n]*cos(Theta0*Z[n]-Beta)-W[0][1][n]*sin(Theta0*Z[n]-Beta))/dy/4;
                if((y > -0.4*N) && (y < 0.6*N) && (x > -0.5*N) && (x < 0.5*N))
                    flowline_r.lineTo(x, y);
                else flowline_r.moveTo(x, y);
            }
            space += pow(double(i+0.00001)/numlines, 1)*100;
        }
        painter->drawPath(flowline_l);
        painter->drawPath(flowline_r);
    }
}
void ColorMapWidget::drawVectorPlot(QPainter *painter) {
    double phi, x, y, r0;
    int x1, y1, x2, y2;
    int K = 20;
    for (i = 0; i < K; i++)
        for (k = 0; k < K; k++)
        {
            x = i*N/K - N/2;
            y = k*N/K - N/2;
            phi = atan2(x,y)/Theta0 + Beta/Theta0;
            r0 = sqrt(x*x+y*y);
            if(fabs(phi) <= 1)
            {
                n = (int)floor(phi/dz+0.5) + kmax/2;
                x1 = i*N/K+5;
                y1 = k*N/K+30;
                x2 = 5+i*N/K + int(30/dy*(U[0][1][n]*sin(Theta0*Z[n]-Beta)+W[0][1][n]*cos(Theta0*Z[n]-Beta))* pow(r0/4+0.1, -0.25)/(pow(1.0001-Z[n]*Z[n], 0.15)));
                y2 = 30+k*N/K + int(30/dy*(U[0][1][n]*cos(Theta0*Z[n]-Beta)-W[0][1][n]*sin(Theta0*Z[n]-Beta))* pow(r0/4+0.1, -0.25)/(pow(1.0001-Z[n]*Z[n], 0.15)));
                painter->drawLine(x1, y1, x2, y2);
                //arrow
                double bac = atan2(y2-y1, x2-x1);
                painter->drawLine(x2, y2, int(x2-7*cos(bac+M_PI/12)), int(y2-7*sin(bac+M_PI/12)));
                painter->drawLine(x2, y2, int(x2-7*cos(bac-M_PI/12)), int(y2-7*sin(bac-M_PI/12)));
            }
        }
}
