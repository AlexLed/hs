TEMPLATE = app
TARGET = hs
MOC_DIR = tmp
OBJECTS_DIR = tmp
DESTDIR = compiled
QT += widgets opengl charts

HEADERS += src/mainwindow.h \
	src/solver.h \
	src/viewarray.h \
	src/glwidget.h \
	src/colormap.h \
        src/graphicswidget.h

SOURCES += src/main.cpp \
	src/mainwindow.cpp \
	src/solver.cpp \
	src/viewarray.cpp \
	src/colormap.cpp \
	src/glwidget.cpp \
        src/graphicswidget.cpp
TRANSLATIONS = translations/hs_ru.ts
RESOURCES = hs.qrc
RC_FILE += hs.rc
