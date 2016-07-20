#ifndef VIEWARRAY_H
#define VIEWARRAY_H

#include <QWidget>
#include "solver.h"

class QTableWidget;

class ViewArray : public QWidget
{
	Q_OBJECT
public:
	ViewArray(QWidget *parent);
	QTableWidget *tableWidget;
private slots:
	void show(double ***array);
public slots:
	void showU() { show(U); }
	void showW() { show(W); }
	void showH() { show(H); }
	void showV() { show(V); }
	void showT() { show(T); }
	void showOther();
};

#endif
