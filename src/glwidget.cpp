#include <QtGui>
#include <QtOpenGL>
#include <QGLWidget>
#include <QGL>
#include "glwidget.h"
#include "math.h"
#include "solver.h"

GLWidget::GLWidget(QWidget *parent)
     : QGLWidget(parent)
{
    xRot = 0;
    yRot = 0;
    zRot = 0;
    gear1Rot = 0;
    scale = 8;
    plotType = "void";
    show_axis = true;
}

GLWidget::~GLWidget()
{
    makeCurrent();
}

void GLWidget::setXRotation(int angle)
{
    normalizeAngle(&angle);
    if (angle != xRot) {
        xRot = angle;
        emit xRotationChanged(angle);
        updateGL();
    }
}

void GLWidget::setYRotation(int angle)
{
    normalizeAngle(&angle);
    if (angle != yRot) {
        yRot = angle;
        emit yRotationChanged(angle);
        updateGL();
    }
}

void GLWidget::setZRotation(int angle)
{
    normalizeAngle(&angle);
    if (angle != zRot) {
        zRot = angle;
        emit zRotationChanged(angle);
        updateGL();
    }
}

void GLWidget::initializeGL()
{
    //static const GLfloat lightPos[4] = { 5.0f, 5.0f, 10.0f, 1.0f };
    static const GLfloat reflectance1[4] = { 0.8f, 0.1f, 0.0f, 1.0f };
    //static const GLfloat reflectance2[4] = { 0.0f, 0.8f, 0.2f, 1.0f };
    //static const GLfloat reflectance3[4] = { 0.2f, 0.2f, 1.0f, 1.0f };

//     glLightfv(GL_LIGHT0, GL_POSITION, lightPos);
//     glEnable(GL_LIGHTING);
//     glEnable(GL_LIGHT0);
    //glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    //glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
}

void GLWidget::paintGL()
{
//	GLUquadricObj* QuadrObj;
//	QuadrObj=gluNewQuadric();
    static const GLfloat LightDiffuse[]= { 0.8f, 0.8f, 0.8f, 1.0f };
    //static const GLfloat dir[3] = {-1,-1,-1};
    //static const GLfloat color[4] = { 0.1f, 0.3f, 0.4f, 1.0f };
    static const GLfloat ambient0[]={0.4, 0.4, 0.4, 1.0};

    //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

//    glLightfv(GL_LIGHT0, GL_POSITION, lightPos1);
//    glLightfv(GL_LIGHT0, GL_AMBIENT, ambient0);
//   glLightfv(GL_LIGHT0, GL_DIFFUSE, LightDiffuse);
    //glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, dir);
//    glLightfv(GL_LIGHT1, GL_POSITION, lightPos2);

//    glEnable(GL_LIGHTING);
//    glEnable(GL_LIGHT0);
    //glEnable(GL_LIGHT1);

//    glShadeModel(GL_SMOOTH);
//    glEnable(GL_COLOR_MATERIAL);
//    glEnable(GL_AUTO_NORMAL);

//    glPushMatrix();
//    glRotated(xRot / 16.0, 1.0, 0.0, 0.0);
//    glRotated(yRot / 16.0, 0.0, 1.0, 0.0);
//    glRotated(zRot / 16.0, 0.0, 0.0, 1.0);

//    glScalef(scale, scale, scale);
//    glColor4f(1.0, 1.0, 0.0, 1.0);

    if(show_axis) {
//        glPushMatrix();
//        glTranslatef(-0.5, 0.0, 0.0);
//        glColor3f(1.0, 0.0, 0.0);
//		gluCylinder(QuadrObj, 0.005, 0.005, 0.6, 70, 70); //axis
//        glTranslatef(0.0, 0.0, 0.6);
//		gluCylinder(QuadrObj, 0.015, 0.0, 0.05, 70, 70); //arrow
        //renderText(0, -0.01, 0.1, 0.2, Qt::AlignCenter, "Z"/*, QFont("Arial", 10, QFont::Bold)*/);
//        glTranslatef(0.0, 0.0, -0.6);
//        glRotated(-90, 1.0, 0.0, 0.0);
//        glColor3f(0.0, 1.0, 0.0);
//		gluCylinder(QuadrObj, 0.005, 0.005, 0.6, 70, 70); //axis
//        glTranslatef(0.0, 0.0, 0.6);
//		gluCylinder(QuadrObj, 0.015, 0.0, 0.05, 70, 70); //arrow
        //renderText(-0.01, 0.01, 0.1, "Y", QFont("Arial", 10, QFont::Bold));
//        glTranslatef(0.0, 0.0, -0.6);
//        glRotated(90, 0.0, 1.0, 0.0);
//        glColor3f(0.0, 0.0, 1.0);
//		gluCylinder(QuadrObj, 0.005, 0.005, 1.3, 70, 70); //axis
//        glTranslatef(0.0, 0.0, 1.3);
//		gluCylinder(QuadrObj, 0.015, 0.0, 0.05, 70, 70); //arrow
        //renderText(0.01, -0.01, 0.1, "X", QFont("Arial", 10, QFont::Bold));
//        glTranslatef(0.0, 0.0, -1.3);
        //glTranslatef(-2.0, 0.0, 0.0);
//        glPopMatrix();
    }

//    glColor3f(1.0, 0.0, 0.0);
//    glTranslated(-0.5, 0.0, 0.0);

//	glPushMatrix();
//	glMaterialfv(GL_FRONT, GL_AMBIENT, color);
//	glRotated(-90, 0.0, 0.0, 1.0);
//	glRotated(-90, 0.0, 1.0, 0.0);
//	gluPartialDisk(QuadrObj, 0, 1, 100, 1, -Theta0d, 2*Theta0d);
//	glPopMatrix();

//    glColor3f(0.0, 0.0, 0.0);
    //glLineWidth(1);
/*
    if (plotType == "U") {
        for(int i = 0; i < imax; i++) {
            glBegin(GL_LINE_STRIP);
                for (int k = 0; k < kmax; k++) glVertex3f(X[i], U[i][30][k], Z[k]/2);
            glEnd();
        }
        for(int k = 0; k < kmax; k++) {
            glBegin(GL_LINE_STRIP);
                for (int i = 0; i < imax; i++) glVertex3f(X[i], U[i][30][k], Z[k]/2);
            glEnd();
        }
    }
    else if (plotType == "W") {
        for(int i = 0; i < imax; i++) {
            glBegin(GL_LINE_STRIP);
                for (int k = 0; k < kmax; k++) glVertex3f(X[i], W[i][30][k], Z[k]/2);
            glEnd();
        }
        for(int k = 0; k < kmax; k++) {
            glBegin(GL_LINE_STRIP);
                for (int i = 0; i < imax; i++) glVertex3f(X[i], W[i][30][k], Z[k]/2);
            glEnd();
        }
    }
    else if (plotType == "H") {
        for(int i = 0; i < imax; i++) {
            glBegin(GL_LINE_STRIP);
                for (int k = 0; k < kmax; k++) glVertex3f(X[i], H[i][30][k], Z[k]/2);
            glEnd();
        }
        for(int k = 0; k < kmax; k++) {
            glBegin(GL_LINE_STRIP);
                for (int i = 0; i < imax; i++) glVertex3f(X[i], H[i][30][k], Z[k]/2);
            glEnd();
        }
    }
    else if (plotType == "Pressure") {
        for(int i = 0; i < imax; i++) {
            glBegin(GL_LINE_STRIP);
                for (int k = 0; k < kmax; k++) glVertex3f(X[i], P[i][k]-0.5, Z[k]/2);
            glEnd();
        }
        for(int k = 0; k < kmax; k++) {
            glBegin(GL_LINE_STRIP);
                for (int i = 0; i < imax; i++) glVertex3f(X[i], P[i][k]-0.5, Z[k]/2);
            glEnd();
        }
    }
    else if(plotType == "Delta") {
        for(int i = 0; i < imax; i++) {
            glBegin(GL_LINE_STRIP);
                for (int k = 0; k < kmax; k++) glVertex3f(X[i], Delta[i][k], Z[k]/2);
            glEnd();
        }
        for(int k = 0; k < kmax; k++) {
            glBegin(GL_LINE_STRIP);
                for (int i = 0; i < imax; i++) glVertex3f(X[i], Delta[i][k], Z[k]/2);
            glEnd();
        }
//		for(int i = 0; i < NX-1; i++) {
//			glBegin(GL_QUAD_STRIP);
//				for (int k = 0; k < NZ; k++) {
//					glColor3f(1.0f*(1-Delta[0][k]*pow((1-Z[k]*Z[k]), 0.75)*pow(X[i], 0.75)),
//							  0.2f*(1-Delta[0][k]*pow((1-Z[k]*Z[k]), 0.75)*pow(X[i], 0.75)),
//							  1.0f*Delta[0][k]*pow((1-Z[k]*Z[k]), 0.75)*pow(X[i], 0.75));
//					glVertex3f(X[i]*cos(Theta0), Delta[0][k]*pow((1-Z[k]*Z[k]), 0.75)*pow(X[i], 0.75)/4, Z[k]*sin(Theta0)*i/(NX-1));
//					glColor3f(1.0f*(1-Delta[0][k]*pow((1-Z[k]*Z[k]), 0.75)*pow(X[i+1], 0.75)),
//							  0.2f*(1-Delta[0][k]*pow((1-Z[k]*Z[k]), 0.75)*pow(X[i+1], 0.75)),
//							  1.0f*Delta[0][k]*pow((1-Z[k]*Z[k]), 0.75)*pow(X[i+1], 0.75));
//					glVertex3f(X[i+1]*cos(Theta0), Delta[0][k]*pow((1-Z[k]*Z[k]), 0.75)*pow(X[i+1], 0.75)/4, Z[k]*sin(Theta0)*(i+1)/(NX-1));
//				}
//			glEnd();
//		}
    }
    else if(plotType == "Grid") {
        for(int i = 0; i < imax; i++) {
            glBegin(GL_LINE_STRIP);
            for (int k = 0; k < kmax; k++) glVertex3f(X[i], Body[i][k], Z[k]*Ze[i]);
            glEnd();
        }
        for(int k = 0; k < kmax; k++) {
            glBegin(GL_LINE_STRIP);
                for (int i = 0; i < imax; i++) glVertex3f(X[i], Body[i][k], Z[k]*Ze[i]);
            glEnd();
        }
    }
    else if(plotType == "Body") {
        for(int i = 0; i < imax; i++) {
            glBegin(GL_LINE_STRIP);
                for (int k = 0; k < kmax/2; k++) glVertex3f(X[i], DeltaE[i][k], Z[k]*Ze[i]);
                for (int k = kmax/2-1; k < kmax; k++) glVertex3f(X[i], Body[i][k], Z[k]*Ze[i]);
            glEnd();
        }
        for(int k = 0; k < kmax/2; k++) {
            glBegin(GL_LINE_STRIP);
                for (int i = 0; i < imax; i++) glVertex3f(X[i], DeltaE[i][k], Z[k]*Ze[i]);
            glEnd();
        }
        for(int k = kmax/2-1; k < kmax; k++) {
            glBegin(GL_LINE_STRIP);
                for (int i = 0; i < imax; i++) glVertex3f(X[i], Body[i][k], Z[k]*Ze[i]);
            glEnd();
        }
        for(int j = 0; j < jmax; j++) {
            glBegin(GL_LINE_STRIP);
                for (int i = 0; i < imax; i++) glVertex3f(X[i], DeltaE[i][kmax/2-1]*j*1.0/jmax, Z[kmax/2-1]*Ze[i]);
            glEnd();
        }
    }
    else if(plotType == "void") {
            glBegin(GL_TRIANGLE_FAN);
                glColor3f(0, 1, 0);
                glVertex3d(0, 0, 0);
                glColor3f(0, 0, 1);
                for (int k = 0; k < kmax; k++)
                    glVertex3d(cos(Theta0*Z[k]), (1 - Z[k]*Z[k])/4, sin(Theta0*Z[k]));
            glEnd();
    }
*/
    //gluQuadricDrawStyle (QuadrObj, GLU_LINE);
    //gluSphere(QuadrObj, 1, 70, 70);
//	gluDeleteQuadric(QuadrObj);

//     drawGear(gear1, -3.0, -2.0, 0.0, gear1Rot / 16.0);
//     drawGear(gear2, +3.1, -2.0, 0.0, -2.0 * (gear1Rot / 16.0) - 9.0);
//
//     glRotated(+90.0, 1.0, 0.0, 0.0);
//     drawGear(gear3, -3.1, -1.8, -2.2, +2.0 * (gear1Rot / 16.0) - 2.0);

//    glPopMatrix();
}

void GLWidget::resizeGL(int width, int height)
{
    int side = qMin(width, height);
    //glViewport((width - side) / 2, (height - side) / 2, side, side);
//    glMatrixMode(GL_PROJECTION);
//    glLoadIdentity();
//    glFrustum(-1.0, +1.0, -1.0, 1.0, 5.0, 60.0);
//    glMatrixMode(GL_MODELVIEW);
//    glLoadIdentity();
//    glTranslated(0.0, 0.0, -40.0);
}

void GLWidget::mousePressEvent(QMouseEvent *event)
{
    lastPos = event->pos();
}

void GLWidget::mouseMoveEvent(QMouseEvent *event)
{
    int dx = event->x() - lastPos.x();
    int dy = event->y() - lastPos.y();

    if (event->buttons() & Qt::LeftButton) {
        setXRotation(xRot + 8 * dy);
        setYRotation(yRot + 8 * dx);
    } else if (event->buttons() & Qt::RightButton) {
        setXRotation(xRot + 8 * dy);
        setZRotation(zRot + 8 * dx);
    }
    lastPos = event->pos();
}

void GLWidget::wheelEvent(QWheelEvent *e)
{
    e->delta() > 0 ? scale += scale*0.1f : scale -= scale*0.1f;
    updateGL();
}

GLuint GLWidget::makeGear(const GLfloat *reflectance, GLdouble innerRadius,
                          GLdouble outerRadius, GLdouble thickness,
                          GLdouble toothSize, GLint toothCount)
{
    const double Pi = 3.14159265358979323846;
//    GLuint list = glGenLists(1);
//    glNewList(list, GL_COMPILE);
//    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, reflectance);

    GLdouble r0 = innerRadius;
    GLdouble r1 = outerRadius - toothSize / 2.0;
    GLdouble r2 = outerRadius + toothSize / 2.0;
    GLdouble delta = (2.0 * Pi / toothCount) / 4.0;
    GLdouble z = thickness / 2.0;
    int i, j;

//    glShadeModel(GL_FLAT);
/*
    for (i = 0; i < 2; ++i) {
        GLdouble sign = (i == 0) ? +1.0 : -1.0;

        glNormal3d(0.0, 0.0, sign);

        glBegin(GL_QUAD_STRIP);
        for (j = 0; j <= toothCount; ++j) {
            GLdouble angle = 2.0 * Pi * j / toothCount;
            glVertex3d(r0 * cos(angle), r0 * sin(angle), sign * z);
            glVertex3d(r1 * cos(angle), r1 * sin(angle), sign * z);
            glVertex3d(r0 * cos(angle), r0 * sin(angle), sign * z);
            glVertex3d(r1 * cos(angle + 3 * delta), r1 * sin(angle + 3 * delta),
                       sign * z);
        }
        glEnd();

        glBegin(GL_QUADS);
        for (j = 0; j < toothCount; ++j) {
            GLdouble angle = 2.0 * Pi * j / toothCount;
            glVertex3d(r1 * cos(angle), r1 * sin(angle), sign * z);
            glVertex3d(r2 * cos(angle + delta), r2 * sin(angle + delta),
                       sign * z);
            glVertex3d(r2 * cos(angle + 2 * delta), r2 * sin(angle + 2 * delta),
                       sign * z);
            glVertex3d(r1 * cos(angle + 3 * delta), r1 * sin(angle + 3 * delta),
                       sign * z);
        }
        glEnd();
    }

    glBegin(GL_QUAD_STRIP);
    for (i = 0; i < toothCount; ++i) {
        for (j = 0; j < 2; ++j) {
            GLdouble angle = 2.0 * Pi * (i + (j / 2.0)) / toothCount;
            GLdouble s1 = r1;
            GLdouble s2 = r2;
            if (j == 1)
                qSwap(s1, s2);

            glNormal3d(cos(angle), sin(angle), 0.0);
            glVertex3d(s1 * cos(angle), s1 * sin(angle), +z);
            glVertex3d(s1 * cos(angle), s1 * sin(angle), -z);

            glNormal3d(s2 * sin(angle + delta) - s1 * sin(angle),
                       s1 * cos(angle) - s2 * cos(angle + delta), 0.0);
            glVertex3d(s2 * cos(angle + delta), s2 * sin(angle + delta), +z);
            glVertex3d(s2 * cos(angle + delta), s2 * sin(angle + delta), -z);
        }
    }
    glVertex3d(r1, 0.0, +z);
    glVertex3d(r1, 0.0, -z);
    glEnd();

    glShadeModel(GL_SMOOTH);

    glBegin(GL_QUAD_STRIP);
    for (i = 0; i <= toothCount; ++i) {
        GLdouble angle = i * 2.0 * Pi / toothCount;
        glNormal3d(-cos(angle), -sin(angle), 0.0);
        glVertex3d(r0 * cos(angle), r0 * sin(angle), +z);
        glVertex3d(r0 * cos(angle), r0 * sin(angle), -z);
    }
    glEnd();

    glEndList();
*/
    return 0; //list;
}

void GLWidget::drawGear(GLuint gear, GLdouble dx, GLdouble dy, GLdouble dz,
                        GLdouble angle)
{
//    glPushMatrix();
//    glTranslated(dx, dy, dz);
//    glRotated(angle, 0.0, 0.0, 1.0);
//    glCallList(gear);
//    glPopMatrix();
}

void GLWidget::normalizeAngle(int *angle)
{
    while (*angle < 0)
        *angle += 360 * 16;
    while (*angle > 360 * 16)
        *angle -= 360 * 16;
}

//void GLWidget::initializeGL()
//{
//	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
//	glEnable(GL_DEPTH_TEST);
//	glHint(GL_PERSPECTIVE_CORRECTION_HINZ, GL_NICEST);
//	glMatrixMode(GL_PROJECTION);
//	glLoadIdentity();
//	glOrtho(-.5, .5, .5, -.5, -1000, 1000);
//	glMatrixMode(GL_MODELVIEW);
//	glLoadIdentity();
//	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
//
//	//makeObject();
//}
//void GLWidget::resizeGL(int w, int h)
//{
//	 glViewport(0, 0, w, h);
//     glMatrixMode(GL_PROJECTION);
//     glLoadIdentity();
//     float aspect = w/(float)(h ? h : 1);
//     glFrustum(-aspect, aspect, -1, 1, 10, 100);
//     glTranslatef(-0.5f, -0.5f, -0.5f);
//     glTranslatef(0.0f, 0.0f, -15.0f);
//}
//void GLWidget::paintGL()
//{
//	double a, d, max, w = 1;
//	GLfloat x1, x2, y1, y2, color_r, color_g, color_b;
//	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//	glLoadIdentity();
//	glOrtho(0.0f, 400.0f, 400.0f, 0.0f, -1.0f, 1.0f);
//	for(i = 0; i < N; i++)
//		for(k = 0; k < N; k++)
//			if( max < array_for_drawing[i][k] ) max = array_for_drawing[i][k];
//	a = 1.0/max;
//	glBegin(GL_QUADS);
//	for(int i = 0; i < N; ++i)
//		for(int k = 0; k < N; ++k)
//		{
//			d = a*array_for_drawing[i][k];
//			if(d == 0) {color_r = color_g = color_b = 1.0f;}
//			else
//			{
//				if(d < 0.25) {color_r = 0.0f; color_g = 1.0f*d/0.25; color_b = 1.0f;}
//				else if(d < 0.5) {color_r = 0.0f; color_g = 1.0f; color_b = 1.0f*(2-d/0.25);}
//				else if(d < 0.75) {color_r = 1.0f*(d-0.5)/0.25; color_g = 1.0f; color_b = 0.0f;}
//				else {color_r = 1.0f; color_g = 1.0f*(1-(d-0.75)/0.25); color_b = 0.0f;}
//			}
//			glColor3f(color_r, color_g, color_b);
//			x1 = i*w; x2 = i*w+w;
//			y1 = k*w; y2 = k*w+w;
//			glBegin(GL_QUADS);
//				glVertex3i(x1, y1, 0);
//				//glColor3f(1.0f, 0.0f, 0.0f);
//				glVertex3i(x2, y1, 0);
//				//glColor3f(0.0f, 0.0f, 1.0f);
//				glVertex3i(x2, y2, 0);
//				//glColor3f(0.0f, 0.5f, 0.5f);
//				glVertex3i(x1, y2, 0);
//			glEnd();
//		}
//	glColor3f(0.0f, 1.0f, 0.0f);
//	//glMatrixMode(GL_PROJECTION);
//	//glMatrixMode(GL_MODELVIEW);
//	glTranslated(200,200,0);
//	glRotated(+90.0, 1.0, 0.0, 0.0);
//	glTranslated(-200,-200,0);
//	glBegin(GL_LINES);
//		glVertex3i(200, 200, 0);
//		glVertex3i(200, 100, 0);
//	glEnd();
//	glBegin(GL_LINES);
//		glVertex3i(200, 200, 0);
//		glVertex3i(200, 200, 1);
//	glEnd();
//
//	glEnd();
//	glLoadIdentity();
//	glColor3f(0.0f, 1.0f, 0.0f);
//	glBegin(GL_QUADS);
//		glVertex3i(0, 0, 0);
//		glColor3f(1.0f, 0.0f, 0.0f);
//		glVertex3i(10, 0, 0);
//		glColor3f(0.0f, 0.0f, 1.0f);
//		glVertex3i(0, 10, 0);
//		glColor3f(0.0f, 0.5f, 0.5f);
//		glVertex3i(0, 10, 0);
//	glEnd();
//
/*	glBegin(GL_TRIANGLES);// Drawing Using Triangles
//		glColor3f(1.0f,0.0f,0.0f);glVertex3f( 0.0f, 1.0f, 0.0f);				// Top
//		glColor3f(0.0f,1.0f,0.0f);glVertex3f(-1.0f,-1.0f, 0.0f);				// Bottom Left
//		glColor3f(0.0f,0.0f,1.0f);glVertex3f( 1.0f,-1.0f, 0.0f);				// Bottom Right
//	glEnd();	*/						// Finished Drawing The Triangle
//
//	static float rot = 0.0;
//	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//	glMatrixMode(GL_MODELVIEW);
//	glPushMatrix();
//	glEnable(GL_MULTISAMPLE);
//	glTranslatef(-0.25f, -0.10f, 0.0f);
//	glScalef(0.75f, 1.15f, 0.0f);
//	glRotatef(rot, 0.0f, 0.f, 1.0f);
//	glCallList(list);
//	glPopMatrix();
//
//	glPushMatrix();
//	glDisable(GL_MULTISAMPLE);
//	glTranslatef(0.25f, -0.10f, 0.0f);
//	glScalef(0.75f, 1.15f, 0.0f);
//	glRotatef(rot, 0.0f, 0.0f, 1.0f);
//	glCallList(list);
//	glPopMatrix();
//
//	rot += 0.2f;
//
//	qglColor(Qt::black);
//	renderText(-0.35, 0.4, 0.0, "Multisampling enabled");
//	renderText(0.15, 0.4, 0.0, "Multisampling disabled");
//}
