#include <QtGui>
#include <QTableWidget>
#include <QLayout>
#include "viewarray.h"

ViewArray::ViewArray(QWidget *parent)
    : QWidget(parent)
{
    tableWidget = new QTableWidget(this);
    tableWidget->setRowCount(jmax);
    tableWidget->setColumnCount(kmax);
    QStringList horizontalHeaders;
    for (int k = 0; k < kmax+1; k++)
          horizontalHeaders << QString::number(Z[k]);
    tableWidget->setHorizontalHeaderLabels(horizontalHeaders);

    QVBoxLayout *mainLayout = new QVBoxLayout;
    mainLayout->addWidget(tableWidget);
    mainLayout->setMargin(0);
    setLayout(mainLayout);
}
void ViewArray::show(double ***array)
{
    QStringList verticalHeaders;
    for (int j = 0; j < kmax; j++)
          verticalHeaders << QString::number(dy*j);
    tableWidget->setVerticalHeaderLabels(verticalHeaders);

    for (int row = 0; row < jmax; row++)
        for (int column = 0; column < kmax; column++) {
            QTableWidgetItem *item = new QTableWidgetItem(QString::number(array[LX][row][column]));
            tableWidget->setItem(row,column,item);
        }
    tableWidget->resizeColumnsToContents();
    tableWidget->resizeRowsToContents();
}
void ViewArray::showOther()
{
    QStringList verticalHeaders;
    verticalHeaders << "P" << "Delta" << "dP/dz" << "Prel" << "dD/dP" << "J[k]";
    tableWidget->setVerticalHeaderLabels(verticalHeaders);
    for (int column = 0; column < kmax; column++){
        QTableWidgetItem *item = new QTableWidgetItem(QString::number(P[LX][column]));
        tableWidget->setItem(0,column,item);
        for (int row = 1; row < jmax; row++) tableWidget->setItem(row, column, NULL);
    }
    for (int column = 0; column < kmax; column++){
        QTableWidgetItem *item = new QTableWidgetItem(QString::number(Delta[LX][column]));
        tableWidget->setItem(1,column,item);
    }
    for (int column = 0; column < kmax; column++){
        QTableWidgetItem *item = new QTableWidgetItem(QString::number((P[LX][column]-P[LX][column-1])/dz));
        tableWidget->setItem(2,column,item);
    }
    for (int column = 0; column < kmax; column++){
        QTableWidgetItem *item = new QTableWidgetItem(QString::number(Prel[LX][column]));
        tableWidget->setItem(3,column,item);
    }
    for (int column = 0; column < kmax; column++){
        QTableWidgetItem *item = new QTableWidgetItem(QString::number((Delta[LX][column]-Delta[LX][column-1])/(P[LX][column]-P[LX][column-1])));
        tableWidget->setItem(4,column,item);
    }
    for (int column = 0; column < kmax; column++){
        QTableWidgetItem *item = new QTableWidgetItem(QString::number(J[column]));
        tableWidget->setItem(5,column,item);
    }
    tableWidget->resizeColumnsToContents();
    tableWidget->resizeRowsToContents();
}
