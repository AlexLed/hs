#ifndef GLWIDGET_H
#define GLWIDGET_H
#include <QGLWidget>

class GLWidget : public QGLWidget
{
	Q_OBJECT

public:
	GLWidget(QWidget *parent = 0);
	~GLWidget();

	QString plotType;
	bool show_axis;

	int xRotation() const { return xRot; }
	int yRotation() const { return yRot; }
	int zRotation() const { return zRot; }

public slots:
	void show3dU() { plotType = "U"; updateGL(); }
	void show3dW() { plotType = "W"; updateGL(); }
	void show3dH() { plotType = "H"; updateGL(); }
	void show3dPressure() { plotType = "Pressure"; updateGL(); }
	void show3dDelta() { plotType = "Delta"; updateGL(); }
	void show3dBody() { plotType = "Body"; updateGL(); }
	void showAxis(bool checked) { show_axis = checked; updateGL(); }
	void setXRotation(int angle);
	void setYRotation(int angle);
	void setZRotation(int angle);

signals:
	void xRotationChanged(int angle);
	void yRotationChanged(int angle);
	void zRotationChanged(int angle);

protected:
	void initializeGL();
	void paintGL();
	void resizeGL(int width, int height);
	void mousePressEvent(QMouseEvent *event);
	void mouseMoveEvent(QMouseEvent *event);
	void wheelEvent(QWheelEvent *event);

private slots:

private:
	GLuint makeGear(const GLfloat *reflectance, GLdouble innerRadius,
					GLdouble outerRadius, GLdouble thickness,
					GLdouble toothSize, GLint toothCount);
	void drawGear(GLuint gear, GLdouble dx, GLdouble dy, GLdouble dz,
				  GLdouble angle);
	void normalizeAngle(int *angle);

	int xRot, yRot, zRot;
	int gear1Rot;
	float scale;

	QPoint lastPos;
};

#endif
