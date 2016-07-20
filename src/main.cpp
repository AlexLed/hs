#include <QApplication>

#include "mainwindow.h"

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    qApp->setApplicationName("HS");
    QApplication::setApplicationVersion("0.4.0");
    qApp->setOrganizationName("MIPT");

    MainWindow window;
    window.show();
    if(argc > 1)
      window.solveProblem(argv[1]);
    //window.solve();
    return app.exec();
}
