#include <QPainter>
#include <QPaintEvent>
#include <QPen>
#include <QWidget>
#include <QtGui>
#include "graphicswidget.h"
#include "solver.h"
#include "datadef.h"

RenderArea::RenderArea(QWidget *parent)
   : QWidget(parent)
{
   setPalette(QPalette(QColor(255, 255, 255)));
   setAutoFillBackground(true);
   Margin = 30;
   Nmax = 10;
   labelY = "F(z)";
   curve[0].caption = "U";
   curve[0].pen = QPen(Qt::black, 1, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);
   curve[0].show = true;
   curve[1].caption = "W";
   curve[1].pen = QPen(Qt::blue, 1, Qt::DashLine, Qt::RoundCap, Qt::RoundJoin);
   curve[1].show = true;
   curve[2].caption = "H";
   curve[2].pen = QPen(Qt::red, 1, Qt::DashDotLine, Qt::RoundCap, Qt::RoundJoin);
   curve[2].show = true;
   curve[3].caption = "V";
   curve[3].pen = QPen(Qt::black, 1, Qt::DotLine, Qt::RoundCap, Qt::RoundJoin);
   curve[3].show = true;
   curve[4].caption = "TauU";
   curve[4].pen = QPen(Qt::magenta, 1, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);
   curve[4].show = false;
   curve[5].caption = "TauW";
   curve[5].pen = QPen(Qt::blue, 1, Qt::DashLine, Qt::RoundCap, Qt::RoundJoin);
   curve[5].show = false;
   curve[6].caption = "TauH";
   curve[6].pen = QPen(Qt::red, 1, Qt::DashDotLine, Qt::RoundCap, Qt::RoundJoin);
   curve[6].show = false;

//   createActions();
   setupAxis(-1, 1, -0.4, 1, 0.2, 0.5, tr("%1").arg(QChar(0x03b8)), "F(z)");

//	zoomInButton = new QToolButton(this);
//	zoomInButton->setIcon(QIcon(":/images/zoomin.png"));
//	zoomInButton->adjustSize();
//	//connect(zoomInButton, SIGNAL(clicked()), this, SLOT(zoomIn()));
//
//	zoomOutButton = new QToolButton(this);
//	zoomOutButton->setIcon(QIcon(":/images/zoomout.png"));
//	zoomOutButton->adjustSize();
//	//connect(zoomOutButton, SIGNAL(clicked()), this, SLOT(zoomOut()));
//
//	int x = width() - (zoomInButton->width()
//					   + zoomOutButton->width() + 10);
//	zoomInButton->move(x, 5);
//	zoomOutButton->move(x + zoomInButton->width() + 5, 5);
}
void RenderArea::paintEvent(QPaintEvent *event)
{
   QPainter painter(this);
   paintAxis(painter);
   plotGraph(painter);
}
void RenderArea::resizeEvent(QResizeEvent * /* event */)
{
//	int x = width() - (zoomInButton->width()
//					   + zoomOutButton->width() + 10);
//	zoomInButton->move(x, 5);
//	zoomOutButton->move(x + zoomInButton->width() + 5, 5);
   this->update();
}
//void RenderArea::createActions(){
//    exportPlotAsPictureAct = new QAction(QIcon(":/images/document-open.png"), tr("&Save graphic as a picture"), this);
//    exportPlotAsPictureAct->setShortcut(tr("Ctrl+G"));
//    exportPlotAsPictureAct->setToolTip(tr("Save graphic as a picture"));
//    connect(exportPlotAsPictureAct, SIGNAL(triggered()), this, SLOT(exportPlotAsPicture()));
//}
void RenderArea::setupAxis(double min_x, double max_x, double min_y, double max_y, double step_x, double step_y, QString label_x, QString label_y) {
    minX = min_x;
    maxX = max_x;
    minY = min_y;
    maxY = max_y;
    stepX = step_x;
    stepY = step_y;
    labelX = label_x;
    labelY = label_y;
}
void RenderArea::plotGraph(QPainter &painter)
{
   painter.save();
   painter.setRenderHint(QPainter::Antialiasing, true);
   painter.setFont(QFont("Times",13,QFont::Bold));
   QPen penU = QPen(Qt::black, 1, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin),
       penW = QPen(Qt::blue, 1, Qt::DashLine, Qt::RoundCap, Qt::RoundJoin),
       penH = QPen(Qt::red, 1, Qt::DashDotLine, Qt::RoundCap, Qt::RoundJoin),
       penV = QPen(Qt::black, 1, Qt::DotLine, Qt::RoundCap, Qt::RoundJoin);
   // Cartesian-Windows conversion coordinates
//	m1=(double)(rect().width()-2*d)/(maxX-minX);
//	c1=(double)d-minX*m1;
//	m2=(double)(2*d-rect().height())/(maxY-minY);
//	c2=(double)rect().height()-d-minY*m2;
   sx = (rect().width()-2*Margin)/(maxX-minX);
   sy = (2*Margin-rect().height())/(maxY-minY);
   painter.setClipRect(0, 0, rect().width(), rect().height());
   painter.translate(Margin-minX*sx, Margin-maxY*sy);

   int i,j,k,n;
   QPainterPath plot;
   curve[0].data = U;
   curve[1].data = W;
   curve[2].data = H;
   curve[3].data = V;

   if(labelY == "F(x)") {
      for(n = 0; n < 4/*Nmax*/; n++)
         if(curve[n].show)
         {
            plot.moveTo(minX*sx, curve[n].data[0][LY][LZ]*sy);
            for (i = 1; i < imax; i++)
               if(fabs(curve[n].data[i][LY][LZ]) < maxY*2) plot.lineTo((minX+i*dx)*sx, curve[n].data[i][LY][LZ]*sy);
            painter.setPen(curve[n].pen);
            painter.drawPath(plot);
            painter.drawText(int(minX*sx), int(curve[n].data[0][LY][LZ]*sy), curve[n].caption);
         }
      QPainterPath plotP;
      plotP.moveTo(minX*sx, P[0][LZ]*sy);
      for (i = 1; i < imax+1; i++)
         /*if(fabs(P[i][LZ]) < maxY*2)*/ plotP.lineTo((minX+i*dx)*sx, P[i][LZ]*sy);
      painter.setPen(curve[4].pen);
      painter.drawPath(plotP);
      painter.drawText(int(minX*sx), int(P[0][LZ]*sy), "P");
      QPainterPath plotD;
      plotD.moveTo(minX*sx, Delta[0][LZ]*sy);
      for (i = 1; i < imax; i++)
         /*if(fabs(Delta[i][LZ]) < maxY*2) */plotD.lineTo((minX+i*dx)*sx, Delta[i][LZ]*sy);
      painter.setPen(curve[5].pen);
      painter.drawPath(plotD);
      painter.drawText(int(minX*sx), int(Delta[0][LZ]*sy), "D");
  }
  else if(labelX == "F(y)") {
      for(i = 0; i < 4/*Nmax*/; i++)
         if(curve[i].show) {
            QPainterPath plot;
            plot.moveTo(curve[i].data[LX][0][LZ]*sx, minY*sy);
            for (j = 1; j < jmax-1; j++)
               if(fabs(curve[i].data[LX][j][LZ]) < maxX*2) plot.lineTo(curve[i].data[LX][j][LZ]*sx, (minY+j*dy)*sy);
            painter.setPen(curve[i].pen);
            painter.drawPath(plot);
            painter.drawText(int(curve[i].data[LX][j-1][LZ]*sx)-6, int(maxY*sy)+3, curve[i].caption);
        }
      QPainterPath plot;
      plot.moveTo(Mach[LX][0][LZ]*sx, minY*sy);
      for (j = 1; j < jmax-1; j++)
          plot.lineTo(Mach[LX][j][LZ]*sx, (minY+j*dy)*sy);
          //plot.lineTo((U[0][j][LZ]*sin(Theta0*Z[LZ]-Beta)+W[0][j][LZ]*cos(Theta0*Z[LZ]-Beta))*sx, (minY+j*dy)*sy);
      painter.setPen(curve[4].pen);
      painter.drawPath(plot);
  }
  else if(labelY == "F(z)")
   {
      for(i = 0; i < 4/*Nmax*/; i++)
         if(curve[i].show)
         {
            QPainterPath plot;
            painter.setPen(curve[i].pen);
            plot.moveTo(minX*sx, curve[i].data[LX][LY][0]*sy);
            for (k = 1; k < kmax; k++)
                if(fabs(curve[i].data[LX][LY][k]) < maxY*2) plot.lineTo((minX+k*dz)*sx, curve[i].data[LX][LY][k]*sy);
            painter.drawPath(plot);
            painter.drawText(int(minX*sx)-16, int(curve[i].data[LX][LY][0]*sy)+6, curve[i].caption);
         }
//	  j = LY+1;
//	  k = 0;
//	  QPainterPath plot;
//	  plot.moveTo(minX*sx, (W[0][j][k+1]-W[0][j][k-1]+W[0][j-1][k+1]-W[0][j-1][k-1])/dz*sy);
//	  for (k = 1; k < NZ; k++)
//		 plot.lineTo((minX+k*dz)*sx, (P[0][k+1]-P[0][k-1])/dz/*(W[0][j][k+1]-W[0][j][k-1])/dz*/*sy);
//	  painter.setPen(curve[5].pen);
//	  painter.drawPath(plot);
//	  painter.drawText(int(minX*sx)-29, int((W[LX][1][0]-W[LX][0][0])/dy*sy)+12, "2");

      if(curve[4].show){
         curve[4].plot.moveTo(minX*sx, (U[LX][1][0]-U[LX][0][0])/dy*sy);
         for (k = 1; k < kmax; k++)
            curve[4].plot.lineTo((minX+k*dz)*sx, (U[LX][1][k]-U[LX][0][k])/dy*sy);
         painter.setPen(curve[4].pen);
         painter.drawPath(curve[4].plot);
         painter.drawText(int(minX*sx)-29, int((U[LX][1][0]-U[LX][0][0])/dy*sy)+12, "TauU");}
      if(curve[5].show){
         curve[5].plot.moveTo(minX*sx, (W[LX][1][0]-W[LX][0][0])/dy*sy);
         for (k = 1; k < kmax; k++)
            curve[5].plot.lineTo((minX+k*dz)*sx, (W[LX][1][k]-W[LX][0][k])/dy*sy);
         painter.setPen(penW);
         painter.drawPath(curve[5].plot);
         painter.drawText(int(minX*sx)-29, int((W[LX][1][0]-W[LX][0][0])/dy*sy)+12, "TauW");}
      if(curve[6].show){
         curve[6].plot.moveTo(minX*sx, (H[LX][1][0]-H[LX][0][0])/dy*sy);
         for (k = 1; k < kmax; k++)
            curve[6].plot.lineTo((minX+k*dz)*sx, (H[LX][1][k]-H[LX][0][k])/dy*sy);
         painter.setPen(penH);
         painter.drawPath(curve[6].plot);
         painter.drawText(int(minX*sx)-29, int((H[LX][1][0]-H[LX][0][0])/dy*sy)+12, "TauH");}
   }
  else if(labelY == "Tau")
   {
      QPainterPath plot;
      plot.moveTo(minX*sx, Tau[LX][1]*sy);
      for (k = 2; k < kmax-1; k++)
         plot.lineTo((minX+k*dz)*sx, Tau[LX][k]*sy);
      curve[5].plot.moveTo(minX*sx, TauW[LX][1]*sy);
      for (k = 2; k < kmax-1; k++)
         curve[5].plot.lineTo((minX+k*dz)*sx, TauW[LX][k]*sy);
      curve[6].plot.moveTo(minX*sx, TauH[LX][1]*sy);
      for (k = 2; k < kmax-1; k++)
         curve[6].plot.lineTo((minX+k*dz)*sx, TauH[LX][k]*sy);


      painter.setPen(penU);
      painter.drawPath(plot);
//	  painter.setPen(penW);
//	  painter.drawPath(curve[5].plot);
//	  painter.setPen(penH);
//	  painter.drawPath(curve[6].plot);
   }
  else if(labelY == "P")
   {
      plot.moveTo(minX*sx, P[LX][0]*sy);
      for (k = 1; k < kmax; k++)
         plot.lineTo((minX+k*dz)*sx, P[LX][k]*sy);
      painter.setPen(penU);
      painter.drawPath(plot);

//	  QPainterPath plot;
//	  double m;
//	  plot.moveTo(minX*sx, P[LX][1]*1/sqrt(sin(Theta0*(1+Z[1])))*sy);
//	  for (k = 1; k < NZ-1; k++) {
//		  if (Z[k] < 0) m = 1/sqrt(sin(Theta0*(1+Z[k])))/sqrt(cos((1+Z[k])*Theta0));
//		  else m = 1/sqrt(1-Z[k])/sqrt(cos(Theta0)/cos(Z[k]*Theta0));
//		  plot.lineTo((minX+k*dz)*sx, P[LX][k]*m*sy);
//	  }
//	  //plot.lineTo((minX+k*dz)*sx, (P[LX][k]-P[LX][k-1])/dz*sy);
//	  painter.setPen(penW);
//	  painter.drawPath(plot);

      QPainterPath plot2;
      j = 10;
      k = 0;
      plot2.moveTo(minX*sx, (1.0-Z[k]*Z[k])/P[0][k]*W[0][j][k]*sy);
      for (k = 1; k < kmax; k++)
         plot2.lineTo((minX+k*dz)*sx, (Z[k]/P[0][k])*sy);
      painter.setPen(curve[4].pen);
      //painter.drawPath(plot2);
   }
    else if(labelY == "D")
   {
      plot.moveTo(minX*sx, Delta[LX][0]*sy);
      for (k = 1; k < kmax; k++)
         plot.lineTo((minX+k*dz)*sx, Delta[LX][k]*sy);
      painter.setPen(penU);
      painter.drawPath(plot);
   }
    else if(labelY == "Disturbances") {
        if (sx > sy)
            sx = -sy;
        else
            sy = -sx;
        plot.moveTo(0, DisturbancesDistribution[0][LZ]*sy); // minus is to change main direction from top to bottom
        for (int angle = 1; angle < 360; angle++)
            plot.lineTo(DisturbancesDistribution[angle][LZ]*sin(angle*Pi/180)*sx,
                        DisturbancesDistribution[angle][LZ]*cos(angle*Pi/180)*sy);
        painter.setPen(penU);
        painter.drawPath(plot);

        QPainterPath plot2;
        plot2.moveTo(minX*sx, DisturbancesDistribution[180][0]*sy);
        for (k = 1; k < kmax; k++)
             plot2.lineTo((minX+k*dz)*sx, DisturbancesDistribution[180][k]*sy);
//        plot2.moveTo(minX*sx, DisturbancesDistribution[270-int(Theta0d+Betad)][k]*sy);
//        if (Theta0d < 90) {
//            for (k = 1; k < NZ/2; k++)
//               plot2.lineTo((minX+k*dz)*sx, DisturbancesDistribution[270-int(Theta0d+Betad)][k]*sy);
//            for (k = NZ/2; k < NZ; k++)
//               plot2.lineTo((minX+k*dz)*sx, DisturbancesDistribution[90+int(Theta0d-Betad)][k]*sy);
//        }
//        else {
//            for (k = 1; (1-fabs(Z[k]))*Theta0d < 90; k++)
//               plot2.lineTo((minX+k*dz)*sx, DisturbancesDistribution[270-int(Theta0d+Betad)][k]*sy);
//            for (; (1-fabs(Z[k]))*Theta0d > 90; k++)
//                plot2.lineTo((minX+k*dz)*sx, DisturbancesDistribution[180-int(-Theta0d*Z[k]+Betad)][k]*sy);
//            for (; k < NZ; k++)
//               plot2.lineTo((minX+k*dz)*sx, DisturbancesDistribution[90+int(Theta0d-Betad)][k]*sy);
//        }
        painter.setPen(curve[4].pen);
        painter.drawPath(plot2);

        painter.setPen(Qt::black);
        if(Theta0d >= 90) painter.translate(rect().width()*0.4, -rect().height()*0.4);
        else painter.translate(rect().width()*0.2, rect().height()*0.3);
        painter.rotate(Betad);
        int height = rect().height()*0.20;
        int width = rect().width()*0.20;
        width = int(height*0.49);
        if(!cartesian) {
           painter.drawPie(-width, -width, width*2, width*2,
                    -(int)(Theta0d+90)*16, (int)(2*Theta0d)*16);
           painter.drawLine(0, 0, 0, width);
           painter.setPen(QPen(Qt::red, 1.0, Qt::SolidLine, Qt::FlatCap, Qt::MiterJoin));
           painter.drawLine(0, 0, int(width*sin(Z[LZ]*Theta0)), int(width*cos(Z[LZ]*Theta0)));
           painter.setPen(QPen(Qt::red, 6.0, Qt::SolidLine, Qt::FlatCap, Qt::MiterJoin));
           painter.drawPoint(int(width*sin(Z[LZ]*Theta0)), int(width*cos(Z[LZ]*Theta0)));
        }
        else if(cartesian) {
           int x, y, i;
           x = int(0.8*width);
           y = int(x/z0);
     //	  if (x > width/2) { x = width*3/4; y = int(x/tan(Theta0)); }
     //	  if (y > height/2) { y = height*2/4; x = int(y*tan(Theta0)); }
           for (i = 0; i < imax-1; i++) {
               painter.drawLine(Ze[i]*y, y*i*dx, Ze[i+1]*y, y*(i+1)*dx);
               painter.drawLine(-Ze[i]*y, y*i*dx, -Ze[i+1]*y, y*(i+1)*dx);
           }
           painter.drawLine(Ze[i]*y, y*i*dx, -Ze[i]*y, y*i*dx);
           painter.drawLine(0, 0, 0, y);
        }
   }
   else
   {
       for(n = 0; n < Ncalc; n++) {
           QPainterPath plot;
           plot.moveTo(minX*sx, data.array[n][Narr][0]*sy);
           for (k = 1; k < kmax; k++)
               plot.lineTo((minX+k*dz)*sx, data.array[n][Narr][k]*sy);
           painter.setPen(curve[n].pen);
           painter.drawPath(plot);
           painter.setFont(QFont("Times", 8));
           painter.drawText(int(maxX*sx)+1, int(data.array[n][Narr][kmax-1]*sy), data.title.at(n));
       }
   }
   painter.restore();
}
void RenderArea::paintAxis(QPainter &painter)
{
   painter.save();
   painter.setRenderHint(QPainter::Antialiasing, true);
   painter.setFont(QFont("Times", 12));
   painter.setPen(QPen(Qt::black, 2, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));

//	// Cartesian-Windows conversion coordinates
//	M1=(double)rect().width()/(maxX-minX);
//	C1=(double)-minX*M1;
//	M2=(double)-rect().height()/(maxY-minY);
//	C2=(double)rect().height()-minY*M2;
//	m1=(double)(rect().width()-2*d)/(maxX-minX);
//	c1=(double)d-minX*m1;
//	m2=(double)(2*d-rect().height())/(maxY-minY);
//	c2=(double)rect().height()-d-minY*m2;

   int arroww = 4, arrowl = 11; //width and length of arrows
   sx = (rect().width()-2*Margin)/(maxX-minX);
   sy = (2*Margin-rect().height())/(maxY-minY);
   painter.translate(Margin-minX*sx, Margin-maxY*sy);
   // Draw & label the x,y axis
   //axis У (directed from top to bottom)
   painter.drawLine(0, int(minY*sy+Margin), 0, int(maxY*sy-Margin));
   //arrow of axis У
   painter.drawLine(0, int(maxY*sy-Margin), 0+arroww, int(maxY*sy-Margin+arrowl));
   painter.drawLine(0, int(maxY*sy-Margin), 0-arroww, int(maxY*sy-Margin+arrowl));
   painter.drawLine(0, int(maxY*sy-Margin+4), 0+arroww, int(maxY*sy-Margin+arrowl));
   painter.drawLine(0, int(maxY*sy-Margin+4), 0-arroww, int(maxY*sy-Margin+arrowl));
   //axis Х (directed from left to right)
   painter.drawLine(int(minX*sx-Margin), 0, int(maxX*sx+Margin), 0);
   //arrow of axis Х
   painter.drawLine(int(maxX*sx+Margin), 0, int(maxX*sx+Margin-arrowl), 0+arroww);
   painter.drawLine(int(maxX*sx+Margin), 0, int(maxX*sx+Margin-arrowl), 0-arroww);
   painter.drawLine(int(maxX*sx+Margin-4), 0, int(maxX*sx+Margin-arrowl), 0+arroww);
   painter.drawLine(int(maxX*sx+Margin-4), 0, int(maxX*sx+Margin-arrowl), 0-arroww);
   //Labels of axis X
   painter.setFont(QFont("Arial", 11));
   if ( (maxX - minX)/stepX < 15 )
       for (double x = minX; x <= maxX; x += stepX) {
           painter.drawLine((int)(x*sx), 0, (int)(x*sx), 5);
           if (round(x*100) != 0) painter.drawText((int)(x*sx)-15, 5, 32, 20, Qt::AlignCenter, QString::number(x));
       }
   else
       for (double x = minX; x <= maxX; x += (maxX - minX)/10) {
           painter.drawLine((int)(x*sx), 0, (int)(x*sx), 5);
           if (round(x*10) != 0) painter.drawText((int)(x*sx)-15, 5, 36, 20, Qt::AlignCenter, QString::number(x));
       }
   //Labels of axis Y
   if ( (maxY - minY)/stepY < 15 )
       for (double y = minY; y <= maxY; y += stepY) {
           painter.drawLine(0, (int)(y*sy), 5, (int)(y*sy));
           if (round(y*100) != 0) painter.drawText(8, (int)(y*sy)-10, 26, 20, Qt::AlignVCenter&&Qt::AlignLeft, QString::number(y));
       }
   else
       for (double y = minY; y <= maxY; y += (maxY - minY)/10) {
           painter.drawLine(0, (int)(y*sy), 5, (int)(y*sy));
           if (round(y*10) != 0) painter.drawText(8, (int)(y*sy)-10, 36, 20, Qt::AlignVCenter&&Qt::AlignLeft, QString::number(y));
       }
   painter.drawText(2, 5, 12, 20, Qt::AlignCenter, QString::number(0));
   //Labels
   painter.setPen(QPen(Qt::black, 2, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
   painter.setFont(QFont("Times", 15));
   painter.drawText((int)(maxX*sx)+Margin-35, -27, 32, 20, Qt::AlignVCenter | Qt::AlignRight, labelX);
   painter.setFont(QFont("Times", 14));
   painter.drawText(7, (int)(maxY*sy)-Margin, 80, 24, Qt::AlignVCenter | Qt::AlignLeft, labelY);

   painter.restore();
}
//void RenderArea::contextMenuEvent(QContextMenuEvent *event)
//{
//    QMenu menu(parentWidget());
//    menu.addAction(exportPlotAsPictureAct);
//    menu.exec(event->globalPos());
//}
//void RenderArea::exportPlotAsPicture() {
//    QString fname = QFileDialog::getSaveFileName(this, tr("Export plot as image"), QString(), tr("Images (*.png *.xpm *.bmp *.jpg )"));
//    if (fname.isEmpty()) return;
//    QImage img(this->width(), this->height(), QImage::Format_RGB32);
//    QPainter p(&img);
//    QPalette pal = palette();
//    QBrush b(pal.color(QPalette::Window));
//    p.fillRect(0, 0, this->width(), this->height(), b);
//    paintAxis(p);
//    plotGraph(p);
//    //drawAll(&p, true, false); // no zoom indicators, please
//    img.save(fname);
//}
