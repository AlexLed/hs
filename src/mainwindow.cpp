#include <QtGui>
#include <QTime>
#include <QTranslator>
#include "mainwindow.h"
#include <QString>

MainWindow::MainWindow() {
    setWindowTitle("HS");
    CreateArrays(false, imax, jmax, kmax, dx, dz, N);
    createActions();
    createMenus();
    createParameters();
    fillForm();
    createTabWidget();
    createToolBars();
    createStatusBar();
    createResults();
    createDockWindows();
    readSettings();
    getParameters();
    enableUI(true);
}
void MainWindow::createTabWidget() {
    tabWidget = new QTabWidget;
    wingTab = new WingTab;
    graphicTab = new GraphicTab;
    colorMapTab = new ColorMapTab;
    glTab = new GLTab;
    tabWidget->setTabShape(QTabWidget::Rounded);
    tabWidget->addTab(wingTab, QIcon(":/images/plane.png"), tr("Wing"));
    tabWidget->addTab(graphicTab, QIcon(":/images/graphics.png"), tr("Graphics"));
    tabWidget->addTab(new ArrayTab, QIcon(":/images/table.png"), tr("Arrays"));
    tabWidget->addTab(colorMapTab, tr("ColorMap"));
    tabWidget->setTabIcon(3, QIcon(":/images/color.png"));
    tabWidget->addTab(glTab, QIcon(":/images/3D.png"), tr("3D plot"));
    this->setCentralWidget(tabWidget);
}
void MainWindow::createActions() {
    //File menu
    openAct = new QAction(QIcon(":/images/document-open.png"), tr("&Open..."), this);
    openAct->setShortcut(tr("Ctrl+O"));
    openAct->setToolTip(tr("Open an existing file"));
    connect(openAct, SIGNAL(triggered()), this, SLOT(open()));

    saveAct = new QAction(QIcon(":/images/document-save.png"), tr("&Save"), this);
    saveAct->setShortcut(tr("Ctrl+S"));
    saveAct->setToolTip(tr("Save results to disc"));
    connect(saveAct, SIGNAL(triggered()), this, SLOT(save()));

    saveAsAct = new QAction(QIcon(":/images/document-save-as.png"), tr("&Save as..."), this);
    saveAsAct->setShortcut(tr("Ctrl+Shift+S"));
    saveAsAct->setToolTip(tr("Save results to disc"));
    connect(saveAsAct, SIGNAL(triggered()), this, SLOT(saveAs()));

    exitAct = new QAction(QIcon(":/images/exit.png"), tr("E&xit"), this);
    exitAct->setShortcut(tr("Ctrl+Q"));
    exitAct->setToolTip(tr("Exit the application"));
    connect(exitAct, SIGNAL(triggered()), qApp, SLOT(quit()));

    //Solve menu
    solveAct = new QAction(QIcon(":/images/go-next-view.png"),tr("So&lve"),this);
    solveAct->setShortcut(tr("F5"));
    solveAct->setToolTip(tr("Start solver"));
    connect(solveAct, SIGNAL(triggered()), this, SLOT(solve()));

    stopAct = new QAction(QIcon(":/images/process-stop.png"),tr("Stop"),this);
    stopAct->setShortcut(tr("F4"));
    stopAct->setToolTip(tr("Stop solver"));
    connect(stopAct, SIGNAL(triggered()), this, SLOT(abort()));

    solveCycleAct = new QAction(QIcon(":/images/cache.png"),tr("Parametric calculations"),this);
    //stopAct->setShortcut(tr("Ctrl+R"));
    solveCycleAct->setToolTip(tr("Parametric calculations"));
    connect(solveCycleAct, SIGNAL(triggered()), this, SLOT(solveCycle()));

    selfTestAct = new QAction(tr("Self Test"), this);
    selfTestAct->setToolTip(tr("Run self test"));
    connect(selfTestAct, SIGNAL(triggered()), this, SLOT(selfTest()));

    wingAct = new QAction(QIcon(":/images/plane.png"),tr("Show wing"),this);
    wingAct->setShortcut(tr("Ctrl+W"));
    wingAct->setToolTip(tr("Show wing geometry"));
    connect(wingAct, SIGNAL(triggered()), this, SLOT(showWing()));

    viewArrayAct = new QAction(QIcon(":/images/table.png"),tr("View arrays"),this);
    viewArrayAct->setShortcut(tr("Ctrl+V"));
    viewArrayAct->setToolTip(tr("View arrays"));
    connect(viewArrayAct, SIGNAL(triggered()), this, SLOT(showArray()));

    colormapAct = new QAction(QIcon(":/images/color.png"),tr("View ColorMap"),this);
    //colormapAct->setShortcut(tr("Ctrl+V"));
    colormapAct->setToolTip(tr("View colormap"));
    connect(colormapAct, SIGNAL(triggered()), this, SLOT(showColorMap()));

    graphAct = new QAction(QIcon(":/images/graphics.png"),tr("Plot graphics"),this);
    graphAct->setShortcut(tr("Ctrl+G"));
    graphAct->setToolTip(tr("Plot graphics"));
    connect(graphAct, SIGNAL(triggered()), this, SLOT(showGraphic()));

    //Preferences menu
    langEnAct = new QAction(QIcon(":/images/english.png"),tr("English"), this);
    langEnAct->setCheckable(true);

    langRuAct = new QAction(QIcon(":/images/russian.png"),tr("Russian"), this);
    langRuAct->setCheckable(true);

    QActionGroup *languageGroup = new QActionGroup(this);
    languageGroup->addAction(langEnAct);
    languageGroup->addAction(langRuAct);
    QSettings settings("settings.ini", QSettings::IniFormat);
    if(settings.value("Language", "en").toString() == "ru")
        langRuAct->setChecked(true);
    else
        langEnAct->setChecked(true);
    connect(languageGroup, SIGNAL(triggered(QAction *)), this, SLOT(changeLang(QAction *)));

    //Help menu
    helpAct = new QAction(QIcon(":/images/help-contents.png"), tr("&Help"), this);
    helpAct->setStatusTip(tr("Help"));

    aboutQtAct = new QAction(QIcon(":/images/qt-logo.png"), tr("About &Qt"), this);
    connect(aboutQtAct, SIGNAL(triggered()), this, SLOT(aboutQt()));

    aboutAct = new QAction(QIcon(":/images/help-about.png"), tr("&About"), this);
    connect(aboutAct, SIGNAL(triggered()), this, SLOT(about()));
}
void MainWindow::createMenus() {
    QMenu *fileMenu = menuBar()->addMenu(tr("&File"));
    fileMenu->addAction(openAct);
    fileMenu->addAction(saveAct);
    fileMenu->addAction(saveAsAct);
    fileMenu->addSeparator();
    fileMenu->addAction(exitAct);

    QMenu *solveMenu = menuBar()->addMenu(tr("&Solve"));
    solveMenu->addAction(solveAct);
    solveMenu->addAction(stopAct);
    solveMenu->addAction(solveCycleAct);
    solveMenu->addAction(selfTestAct);
    solveMenu->addSeparator();
    solveMenu->addAction(wingAct);
    solveMenu->addAction(viewArrayAct);
    solveMenu->addAction(colormapAct);
    solveMenu->addAction(graphAct);

    QMenu *preferencesMenu = menuBar()->addMenu(tr("&Preferences"));
    QMenu *languageMenu = preferencesMenu->addMenu(tr("&Language"));
    languageMenu->addAction(langEnAct);
    languageMenu->addAction(langRuAct);

    viewMenu = menuBar()->addMenu(tr("&View"));

    QMenu *helpMenu = menuBar()->addMenu(tr("&Help"));
    helpMenu->addAction(helpAct);
    helpMenu->addSeparator();
    helpMenu->addAction(aboutQtAct);
    helpMenu->addAction(aboutAct);
}
void MainWindow::createToolBars() {
    QToolBar *fileToolBar = addToolBar(tr("File"));
    fileToolBar->addAction(openAct);
    fileToolBar->addAction(saveAct);
    fileToolBar->setIconSize(QSize(20,20));
    fileToolBar->setMovable(false);

    QToolBar *solveToolBar = addToolBar(tr("Solve"));
    solveToolBar->addAction(solveAct);
    solveToolBar->addAction(stopAct);
    solveToolBar->addSeparator();
    solveToolBar->addAction(solveCycleAct);
    solveToolBar->setMovable(false);
    solveToolBar->setIconSize(QSize(20,20));
    solveAct->setEnabled(true);
    stopAct->setEnabled(false);
}
void MainWindow::createStatusBar() {
    QLabel *newLabel = new QLabel("Ready");
    statusBar()->addWidget(newLabel, 0);
}
void MainWindow::createParameters() {
    int lineedit_w = 45, lineedit_h = 18, label_w = 22, label_h = 18;
//Wing geometry
    QGroupBox *groupBox_geometry = new QGroupBox(tr("Geometry"));
    groupBox_geometry->setFixedSize(150,61);
    QLabel *label_geometryTheta0d = new QLabel(tr("%1 ").arg(QChar(0x03b8)),groupBox_geometry);
    label_geometryTheta0d->setGeometry(QRect(5, 17, label_w, label_h));
    label_geometryTheta0d->setAlignment(Qt::AlignRight);
    label_geometryTheta0d->setFont(QFont("Helvetica",11));
    QLabel *label_geometryTheta1d = new QLabel(QString("%1 ").arg(QChar(0x03a6)),groupBox_geometry);
    label_geometryTheta1d->setGeometry(QRect(55, 19, label_w, label_h));
    label_geometryTheta1d->setAlignment(Qt::AlignRight);
    label_geometryTheta1d->setFont(QFont("Times",11));
    QLabel *label_geometry2 = new QLabel(tr("%1 ").arg(QChar(0x03b2)),groupBox_geometry);
    label_geometry2->setGeometry(QRect(98, 19, label_w, label_h));
    label_geometry2->setAlignment(Qt::AlignRight);
    label_geometry2->setFont(QFont("Times",11));
    QLabel *label_geometry3 = new QLabel(tr("h "),groupBox_geometry);
    label_geometry3->setGeometry(QRect(5, 38, label_w, label_h));
    label_geometry3->setAlignment(Qt::AlignRight);
    label_geometry3->setFont(QFont("Times",11));
    QLabel *label_geometry4 = new QLabel(tr("L "),groupBox_geometry);
    label_geometry4->setGeometry(QRect(78, 38, label_w, label_h));
    label_geometry4->setAlignment(Qt::AlignRight);
    label_geometry4->setFont(QFont("Times",11));
    lineEdit_geometryTheta0d = new QLineEdit(groupBox_geometry);
    lineEdit_geometryTheta0d->setGeometry(QRect(27, 19, lineedit_w*0.66, lineedit_h));
    lineEdit_geometryTheta0d->setToolTip(tr("Angle between bisectrix and leading edge"));
    lineEdit_geometryTheta1d = new QLineEdit(groupBox_geometry);
    lineEdit_geometryTheta1d->setGeometry(QRect(74, 19, lineedit_w*0.66, lineedit_h));
    lineEdit_geometryTheta1d->setToolTip(tr("Angle between bisectrix and leading edge"));
    lineEdit_geometryBetad = new QLineEdit(groupBox_geometry);
    lineEdit_geometryBetad->setGeometry(QRect(118, 19, lineedit_w*0.66, lineedit_h));
    lineEdit_geometryBetad->setToolTip(tr("Angle between bisectrix and U"));
    lineEdit_geometryH = new QLineEdit(groupBox_geometry);
    lineEdit_geometryH->setGeometry(QRect(27, 38, lineedit_w, lineedit_h));
    lineEdit_geometryH->setToolTip(tr("Height of the grid"));
    connect(lineEdit_geometryH, SIGNAL(textEdited(QString)), this, SLOT(changed_h(QString)));
    lineEdit_geometryL = new QLineEdit(groupBox_geometry);
    lineEdit_geometryL->setGeometry(QRect(100, 38, lineedit_w, lineedit_h));
    lineEdit_geometryL->setToolTip(tr("Length of the wing"));

    connect(lineEdit_geometryTheta0d, SIGNAL(textEdited(QString)), this, SLOT(updateWing()));
    connect(lineEdit_geometryTheta1d, SIGNAL(textEdited(QString)), this, SLOT(updateWing()));
    connect(lineEdit_geometryBetad, SIGNAL(textEdited(QString)), this, SLOT(updateWing()));
    connect(lineEdit_geometryL, SIGNAL(textEdited(QString)), this, SLOT(updateWing()));

//Parameters
    QGroupBox *groupBox_parameters = new QGroupBox(tr("Parameters"));
    groupBox_parameters->setFixedSize(150,61);
    QLabel *label_parameters1 = new QLabel(tr("%1 ").arg(QChar(0x03b3)),groupBox_parameters);
    label_parameters1->setGeometry(QRect(5, 19, label_w, label_h));
    label_parameters1->setAlignment(Qt::AlignRight);label_parameters1->setFont(QFont("Times",12));
    QLabel *label_parameters2 = new QLabel(tr("Pr "),groupBox_parameters);
    label_parameters2->setGeometry(QRect(5, 38, label_w, label_h));
    label_parameters2->setAlignment(Qt::AlignRight);
    QLabel *label_parameters3 = new QLabel(tr("Hw "),groupBox_parameters);
    label_parameters3->setGeometry(QRect(78, 19, label_w, label_h));
    label_parameters3->setAlignment(Qt::AlignRight);
    QLabel *label_parametersOmega = new QLabel(tr("%1 ").arg(QChar(0x03c9)),groupBox_parameters);
    label_parametersOmega->setGeometry(QRect(78, 38, label_w, label_h));
    label_parametersOmega->setAlignment(Qt::AlignRight);label_parametersOmega->setFont(QFont("Times",11));

    lineEdit_parametersGamma = new QLineEdit(groupBox_parameters);
    lineEdit_parametersGamma->setGeometry(QRect(27, 19, lineedit_w, lineedit_h));
    lineEdit_parametersGamma->setToolTip(tr("Cp/Cv"));
    lineEdit_parametersPr = new QLineEdit(groupBox_parameters);
    lineEdit_parametersPr->setGeometry(QRect(27, 38, lineedit_w, lineedit_h));
    lineEdit_parametersPr->setToolTip(tr("Prandtl number"));
    lineEdit_parametersHw = new QLineEdit(groupBox_parameters);
    lineEdit_parametersHw->setGeometry(QRect(100, 19, lineedit_w, lineedit_h));
    lineEdit_parametersHw->setToolTip(tr("Enthalpy on the wall"));
    lineEdit_parametersOmega = new QLineEdit(groupBox_parameters);
    lineEdit_parametersOmega->setGeometry(QRect(100, 38, lineedit_w, lineedit_h));
    lineEdit_parametersOmega->setToolTip(tr("Dependence of viscosity on temperature: mu ~ T ^ Omega"));

//Solving
    QGroupBox *groupBox_solving = new QGroupBox(tr("Solving"));
    groupBox_solving->setFixedWidth(160);
    QVBoxLayout *solving_layout = new QVBoxLayout;
    QGridLayout *solving_glayout = new QGridLayout;
    solving_glayout->setContentsMargins(0,0,0,0);
    QLabel *label_solving1 = new QLabel("RU ");
    label_solving1->setGeometry(QRect(5, 19, label_w, label_h));
    label_solving1->setAlignment(Qt::AlignRight);
    QLabel *label_solving2 = new QLabel("RP ");
    label_solving2->setGeometry(QRect(5, 38, label_w, label_h));
    label_solving2->setAlignment(Qt::AlignRight);
    QLabel *label_solving3 = new QLabel("MxI ");
    label_solving3->setGeometry(QRect(78, 19, label_w, label_h));
    label_solving3->setAlignment(Qt::AlignRight);
    QLabel *label_solving4 = new QLabel("RP2 "/*"%1 ".arg(QChar(0x03b1))*/);
    label_solving4->setGeometry(QRect(78, 38, label_w, label_h));
    label_solving4->setAlignment(Qt::AlignRight);
    QLabel *label_solving_AcU = new QLabel("AcU ");
    label_solving_AcU->setGeometry(QRect(5, 57, label_w, label_h));
    label_solving_AcU->setAlignment(Qt::AlignRight);
    QLabel *label_solving_AcPr = new QLabel("AcP ");
    label_solving_AcPr->setGeometry(QRect(78, 57, label_w, label_h));
    label_solving_AcPr->setAlignment(Qt::AlignRight);

    lineEdit_solving_relU = new QLineEdit("", groupBox_solving);
    lineEdit_solving_relU->setGeometry(QRect(27, 19, lineedit_w, lineedit_h));
    lineEdit_solving_relU->setToolTip(tr("Relaxation coefficient for velocity and enthalpy"));
    lineEdit_solving_relP = new QLineEdit(groupBox_solving);
    lineEdit_solving_relP->setGeometry(QRect(27, 38, lineedit_w, lineedit_h));
    lineEdit_solving_relP->setToolTip(tr("Relaxation coefficient for pressure"));
    lineEdit_solving3 = new QLineEdit(groupBox_solving);
    lineEdit_solving3->setGeometry(QRect(100, 19, lineedit_w, lineedit_h));
    lineEdit_solving3->setToolTip(tr("Initial profile of pressure"));
    lineEdit_solving4 = new QLineEdit(groupBox_solving);
    lineEdit_solving4->setGeometry(QRect(100, 38, lineedit_w, lineedit_h));
    lineEdit_solving4->setToolTip(tr("Coefficient in relaxation"));
    lineEdit_solving_AcU = new QLineEdit(groupBox_solving);
    lineEdit_solving_AcU->setGeometry(QRect(27, 57, lineedit_w, lineedit_h));
    lineEdit_solving_AcU->setToolTip(tr("Accuracy for velocity"));
    lineEdit_solving_AcPr = new QLineEdit(groupBox_solving);
    lineEdit_solving_AcPr->setGeometry(QRect(100, 57, lineedit_w, lineedit_h));
    lineEdit_solving_AcPr->setToolTip(tr("Accuracy for pressure"));
    solving_glayout->addWidget(label_solving1,0,1,Qt::AlignRight);
    solving_glayout->addWidget(lineEdit_solving_relU,0,2);
    solving_glayout->addWidget(label_solving2,0,3,Qt::AlignRight);
    solving_glayout->addWidget(lineEdit_solving_relP,0,4);
    solving_glayout->addWidget(label_solving3,1,1,Qt::AlignRight);
    solving_glayout->addWidget(lineEdit_solving3,1,2);
    solving_glayout->addWidget(label_solving4,1,3,Qt::AlignRight);
    solving_glayout->addWidget(lineEdit_solving4,1,4);
    solving_glayout->addWidget(label_solving_AcU,2,1,Qt::AlignRight);
    solving_glayout->addWidget(lineEdit_solving_AcU,2,2);
    solving_glayout->addWidget(label_solving_AcPr,2,3,Qt::AlignRight);
    solving_glayout->addWidget(lineEdit_solving_AcPr,2,4);

//	QButtonGroup *radioTheta = new QButtonGroup;
//	radio1th = new QRadioButton(tr("First order in %1").arg(QChar(0x03b8)),groupBox_solving);
//	radio2th = new QRadioButton(tr("Second order in %1").arg(QChar(0x03b8)),groupBox_solving);
//	radioTheta->addButton(radio1th);
//	radioTheta->addButton(radio2th);
//	radio1th->setGeometry(QRect(7, 77, 140, 16));
//	radio2th->setGeometry(QRect(7, 92, 140, 16));
//	radio1th->setChecked(is_first_order_theta);
//	radio2th->setChecked(!is_first_order_theta);
//	connect(radio1th, SIGNAL(toggled(bool)), this, SLOT(changeOrderTheta(bool)));

    QButtonGroup *radioCoordinates = new QButtonGroup;
    radioC1 = new QRadioButton(tr("Cylindrical"));
    radioC2 = new QRadioButton(tr("Cartesian"));
    radioCoordinates->addButton(radioC1);
    radioCoordinates->addButton(radioC2);
    radioC1->setChecked(!cartesian);
    radioC2->setChecked(cartesian);
    connect(radioC1, SIGNAL(toggled(bool)), this, SLOT(changeCoordinates(bool)));

    QButtonGroup *radioMethod = new QButtonGroup;
    QRadioButton *radioM1 = new QRadioButton(tr("Progonka"));
    QRadioButton *radioM2 = new QRadioButton(tr("Pseudotransient"));
    radioMethod->addButton(radioM1);
    radioMethod->addButton(radioM2);
    radioM1->setChecked(thomas_method);
    radioM2->setChecked(!thomas_method);
    connect(radioM1, SIGNAL(toggled(bool)), this, SLOT(changeMethod(bool)));

    checkBox_graph = new QCheckBox(tr("Plot graphic"));

    QCheckBox *checkBox_MR = new QCheckBox(tr("Modified relaxation"));
    checkBox_MR->setChecked(is_mod_rel);
    connect(checkBox_MR, SIGNAL(toggled(bool)), this, SLOT(changeRel(bool)));

    checkBox_MultiGrid = new QCheckBox(tr("MultiGrid"));
    checkBox_MultiGrid->setChecked(false);

    solving_layout->addLayout(solving_glayout);
    solving_layout->addWidget(radioC1);
    solving_layout->addWidget(radioC2);
    solving_layout->addWidget(radioM1);
    solving_layout->addWidget(radioM2);
    solving_layout->addWidget(checkBox_graph);
    solving_layout->addWidget(checkBox_MR);
    solving_layout->addWidget(checkBox_MultiGrid);
    groupBox_solving->setLayout(solving_layout);

//Grid parameters
    QGroupBox *groupBox_grid = new QGroupBox(tr("Grid"));
    groupBox_grid->setFixedSize(150,82);
    QLabel *label_grid_dx = new QLabel(tr("dx "),groupBox_grid);
    label_grid_dx->setAlignment(Qt::AlignRight);
    QLabel *label_grid_NX = new QLabel(tr("NX "),groupBox_grid);
    label_grid_NX->setAlignment(Qt::AlignRight);
    QLabel *label_grid_dy = new QLabel(tr("dy "),groupBox_grid);
    label_grid_dy->setAlignment(Qt::AlignRight);
    QLabel *label_grid_NY = new QLabel(tr("NY "),groupBox_grid);
    label_grid_NY->setAlignment(Qt::AlignRight);
    QLabel *label_grid_dz = new QLabel(tr("dz "),groupBox_grid);
    label_grid_dz->setAlignment(Qt::AlignRight);
    QLabel *label_grid_NZ = new QLabel(tr("NZ "),groupBox_grid);
    label_grid_NZ->setAlignment(Qt::AlignRight);
    label_grid_dx->setGeometry(QRect(5, 18, label_w, label_h));
    label_grid_NX->setGeometry(QRect(78, 18, label_w, label_h));
    label_grid_dy->setGeometry(QRect(5, 40, label_w, label_h));
    label_grid_NY->setGeometry(QRect(78, 40, label_w, label_h));
    label_grid_dz->setGeometry(QRect(5, 62, label_w, label_h));
    label_grid_NZ->setGeometry(QRect(78, 62, label_w, label_h));

    lineEdit_grid_dx = new QLineEdit(groupBox_grid);
    lineEdit_grid_dx->setToolTip(tr("Grid spacing in x direction"));
    lineEdit_grid_NX = new QLineEdit(groupBox_grid);
    lineEdit_grid_NX->setToolTip(tr("Number of meshes in x direction"));
    lineEdit_grid_dy = new QLineEdit(groupBox_grid);
    lineEdit_grid_dy->setToolTip(tr("Grid spacing in vertical direction"));
    lineEdit_grid_NY = new QLineEdit(groupBox_grid);
    lineEdit_grid_NY->setToolTip(tr("Number of meshes in vertical direction"));
    lineEdit_grid_dz = new QLineEdit(groupBox_grid);
    lineEdit_grid_dz->setToolTip(tr("Grid spacing in z direction"));
    lineEdit_grid_NZ = new QLineEdit(groupBox_grid);
    lineEdit_grid_NZ->setToolTip(tr("Number of meshes in z direction"));
    lineEdit_grid_dx->setGeometry(QRect(27, 18, lineedit_w, lineedit_h));
    lineEdit_grid_NX->setGeometry(QRect(100, 18, lineedit_w, lineedit_h));
    lineEdit_grid_dy->setGeometry(QRect(27, 40, lineedit_w, lineedit_h));
    lineEdit_grid_NY->setGeometry(QRect(100, 40, lineedit_w, lineedit_h));
    lineEdit_grid_dz->setGeometry(QRect(27, 62, lineedit_w, lineedit_h));
    lineEdit_grid_NZ->setGeometry(QRect(100, 62, lineedit_w, lineedit_h));

    connect(lineEdit_grid_dx, SIGNAL(textEdited(QString)), this, SLOT(changed_dx(QString)));
    connect(lineEdit_grid_NX, SIGNAL(textEdited(QString)), this, SLOT(changed_NX(QString)));
    connect(lineEdit_grid_dy, SIGNAL(textEdited(QString)), this, SLOT(changed_dy(QString)));
    connect(lineEdit_grid_NY, SIGNAL(textEdited(QString)), this, SLOT(changed_NY(QString)));
    connect(lineEdit_grid_dz, SIGNAL(textEdited(QString)), this, SLOT(changed_dz(QString)));
    connect(lineEdit_grid_NZ, SIGNAL(textEdited(QString)), this, SLOT(changed_NZ(QString)));

    QVBoxLayout *parameters_layout = new QVBoxLayout;
    parameters_layout->addWidget(groupBox_geometry);
    parameters_layout->addWidget(groupBox_parameters);
    parameters_layout->addWidget(groupBox_solving);
    parameters_layout->addWidget(groupBox_grid);
    //parameters_layout->addWidget(groupBox_results);
    parameters_layout->addStretch(1);
    parameters_layout->setMargin(1);

    parameters = new QWidget();
    parameters->setLayout(parameters_layout);
}
void MainWindow::fillForm() {
    lineEdit_geometryTheta0d->setText(QString::number(Theta0d));
    lineEdit_geometryTheta1d->setText(QString::number(Theta1d));
    lineEdit_geometryBetad->setText(QString::number(Betad));
    lineEdit_geometryH->setText(QString::number(h));
    lineEdit_geometryL->setText(QString::number(L));
    lineEdit_parametersGamma->setText(QString::number(Gamma));
    lineEdit_parametersPr->setText(QString::number(Pr));
    lineEdit_parametersHw->setText(QString::number(Hw));
    lineEdit_parametersOmega->setText(QString::number(Vw));//setText(QString::number(Omega));
    lineEdit_solving_relU->setText(QString::number(relU));
    lineEdit_solving_relP->setText(QString::number(relP));
    lineEdit_solving3->setText(QString::number(MaxIterations));
    lineEdit_solving4->setText(QString::number(relP2/*a0*/));
    lineEdit_solving_AcU->setText(QString::number(AccuracyU));
    lineEdit_solving_AcPr->setText(QString::number(AccuracyPressure));
    lineEdit_grid_dx->setText(QString::number(dx));
    lineEdit_grid_NX->setText(QString::number(imax));
    lineEdit_grid_dy->setText(QString::number(dy));
    lineEdit_grid_NY->setText(QString::number(jmax));
    lineEdit_grid_dz->setText(QString::number(dz));
    lineEdit_grid_NZ->setText(QString::number(kmax));
}
void MainWindow::createResults() {
    QLabel *label_results_labelIterations = new QLabel(tr("Iters"));
    QLabel *label_results_labelTime = new QLabel(tr("Time"));
    QLabel *label_results_labelSpeed = new QLabel(tr("Speed"));
    QLabel *label_results_Edge1_0 = new QLabel(tr("Edge (-1):"));
    label_results_Edge1_0 -> setToolTip(tr("Number of iterations for edge -1"));
    label_results_Edge1_1 = new QLabel(tr("0"));
    label_results_Edge1_2 = new QLabel(tr("0 s"));
    label_results_Edge1_3 = new QLabel(tr("0 ips"));
    QLabel *label_results_Edge2_0 = new QLabel(tr("Edge (1):"));
    label_results_Edge2_0 -> setToolTip(tr("Number of iterations for edge 1"));
    label_results_Edge2_1 = new QLabel(tr("0"));
    label_results_Edge2_2 = new QLabel(tr("0 s"));
    label_results_Edge2_3 = new QLabel(tr("0 ips"));
    QLabel *label_results_Nose_0 = new QLabel(tr("Nose:"));
    label_results_Nose_0 -> setToolTip(tr("Number of iterations for nose"));
    label_results_Nose_1 = new QLabel(tr("0"));
    label_results_Nose_2 = new QLabel(tr("0 s"));
    label_results_Nose_3 = new QLabel(tr("0 ips"));
    QLabel *label_results_Wing_0 = new QLabel(tr("Wing:"));
    label_results_Wing_0 -> setToolTip(tr("Number of iterations for wing"));
    label_results_Wing_1 = new QLabel(tr("0"));
    label_results_Wing_2 = new QLabel(tr("0 s"));
    label_results_Wing_3 = new QLabel(tr("0 ips"));
    QLabel *label_results_TotalTime0 = new QLabel(tr("RunTime:"));
    label_results_TotalTime1 = new QLabel(tr("0 s"));
    label_results_FinishTime0 = new QLabel(tr(""));
    label_results_FinishTime1 = new QLabel(tr(""));

    QGridLayout *results_layout = new QGridLayout;
    results_layout->addWidget(label_results_labelIterations, 0, 1, Qt::AlignCenter);
    results_layout->addWidget(label_results_labelTime, 0, 2, Qt::AlignCenter);
    results_layout->addWidget(label_results_labelSpeed, 0, 3, Qt::AlignCenter);
    results_layout->addWidget(label_results_Edge1_0, 1, 0, Qt::AlignRight);
    results_layout->addWidget(label_results_Edge1_1, 1, 1, Qt::AlignCenter);
    results_layout->addWidget(label_results_Edge1_2, 1, 2, Qt::AlignCenter);
    results_layout->addWidget(label_results_Edge1_3, 1, 3, Qt::AlignCenter);
    results_layout->addWidget(label_results_Edge2_0, 2, 0, Qt::AlignRight);
    results_layout->addWidget(label_results_Edge2_1, 2, 1, Qt::AlignCenter);
    results_layout->addWidget(label_results_Edge2_2, 2, 2, Qt::AlignCenter);
    results_layout->addWidget(label_results_Edge2_3, 2, 3, Qt::AlignCenter);
    results_layout->addWidget(label_results_Nose_0, 3, 0, Qt::AlignRight);
    results_layout->addWidget(label_results_Nose_1, 3, 1, Qt::AlignCenter);
    results_layout->addWidget(label_results_Nose_2, 3, 2, Qt::AlignCenter);
    results_layout->addWidget(label_results_Nose_3, 3, 3, Qt::AlignCenter);
    results_layout->addWidget(label_results_Wing_0, 4, 0, Qt::AlignRight);
    results_layout->addWidget(label_results_Wing_1, 4, 1, Qt::AlignCenter);
    results_layout->addWidget(label_results_Wing_2, 4, 2, Qt::AlignCenter);
    results_layout->addWidget(label_results_Wing_3, 4, 3, Qt::AlignCenter);
    results_layout->addWidget(label_results_TotalTime0, 6, 0, Qt::AlignRight);
    results_layout->addWidget(label_results_TotalTime1, 6, 1, 1, 2, Qt::AlignLeft);
    results_layout->addWidget(label_results_FinishTime0, 7, 0, Qt::AlignRight);
    results_layout->addWidget(label_results_FinishTime1, 7, 1, 1, 2, Qt::AlignLeft);

    results_layout->setMargin(1);

    results = new QWidget;
    results->setLayout(results_layout);
}
void MainWindow::createDockWindows() {
    QMainWindow::setCorner(Qt::BottomLeftCorner, Qt::LeftDockWidgetArea);
    QDockWidget *dock = new QDockWidget(tr("Parameters"), this);
    dock->setFeatures(QDockWidget::DockWidgetMovable/*NoDockWidgetFeatures*/);
    dock->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);
    dock->setWidget(parameters);
    this->addDockWidget(Qt::LeftDockWidgetArea, dock);
    dock->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);
    dock->setMaximumWidth(180);
    viewMenu->addAction(dock->toggleViewAction());

    QDockWidget *dock_results = new QDockWidget(tr("Results"), this);
    dock_results->setFeatures(QDockWidget::DockWidgetMovable);
    //dock_results->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::BottomDockWidgetArea);
    dock_results->setWidget(results);
    this->addDockWidget(Qt::LeftDockWidgetArea, dock_results);
    dock_results->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Expanding);
    dock_results->setMaximumWidth(180);
    viewMenu->addAction(dock_results->toggleViewAction());

//	dock = new QDockWidget(tr("Output"), this);
//	outputConsole = new QTextEdit(dock);
//	dock->setWidget(outputConsole);
//	this->addDockWidget(Qt::BottomDockWidgetArea, dock);
//	viewMenu->addAction(dock->toggleViewAction());
}
void MainWindow::open() {
    QString fileName = QFileDialog::getOpenFileName(this);
    if (!fileName.isEmpty())
    loadFile(fileName);
}
bool MainWindow::save() {
    saveFile(QString("results/%1-%2,Hw=%3,Pr=%4,G=%5,%6x%7x%8.dat").arg(Theta0d).arg(Betad).arg(Hw).arg(Pr).arg(Gamma).arg(imax).arg(jmax).arg(kmax));
    //saveFile(QString("results/Th=%1_B=%2_dy=%3_dz=%4_G=%5_Pr=%6_Hw=%7.dat").arg(Theta0d).arg(Betad).arg(dy).arg(dz).arg(Gamma).arg(Pr).arg(Hw));
    return true;
}
bool MainWindow::saveAs() {
    QString fileName = QFileDialog::getSaveFileName(this);
    saveFile(fileName);
    return true;
}
void MainWindow::about() {
    QMessageBox::about(this, tr("About HS"),
      tr("<b>HS Program</b><br>Revision 99<br>"
      "Й Ledovskiy Aleksey, <b>MIPT</b>, 2008-2012<br>"
      "<i>svn</i>: <a href=\"http://hypers.googlecode.com/svn/trunk\">http://hypers.googlecode.com</a><br>"
      "Built on %1 at %2").arg(__DATE__).arg(__TIME__));
}
void MainWindow::aboutQt() {
    QMessageBox::aboutQt(this, tr("About QT"));
}
void MainWindow::loadFile(const QString &fileName) {
    QFile file(fileName);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
                QMessageBox::warning(this, "HS",
                tr("Cannot read file %1:\n%2.")
                .arg(fileName).arg(file.errorString()));
                return;
        }
    QApplication::setOverrideCursor(Qt::WaitCursor);
    QTextStream in(&file);
    QStringList list;
    Theta0 = in.readLine(0).split("=").at(1).toDouble();
    Theta0d = in.readLine(0).split("=").at(1).toDouble();
    Beta = in.readLine(0).split("=").at(1).toDouble();
    Betad = in.readLine(0).split("=").at(1).toDouble();
    Gamma = in.readLine(0).split("=").at(1).toDouble();
    Pr = in.readLine(0).split("=").at(1).toDouble();
    Hw = in.readLine(0).split("=").at(1).toDouble();
    Omega = in.readLine(0).split("=").at(1).toDouble();
    relU = in.readLine(0).split("=").at(1).toDouble();
    relP = in.readLine(0).split("=").at(1).toDouble();
    AccuracyU = in.readLine(0).split("=").at(1).toDouble();
    AccuracyPressure = in.readLine(0).split("=").at(1).toDouble();
    imax = in.readLine(0).split("=").at(1).toInt();
    dx = in.readLine(0).split("=").at(1).toDouble();
    jmax = in.readLine(0).split("=").at(1).toInt();
    dy = in.readLine(0).split("=").at(1).toDouble();
    kmax = in.readLine(0).split("=").at(1).toInt();
    dz = in.readLine(0).split("=").at(1).toDouble();

    DeleteArrays(false);
    CreateArrays(false, imax, jmax, kmax, dx, dz, N);

    in.readLine(0); //X[NX]
    list = in.readLine(0).split(";");
    for(i = 0; i < imax; i++)
      X[i] = list.at(i).toDouble();

    in.readLine(0); //Z[NZ]
    list = in.readLine(0).split(";");
    for(k = 0; k < kmax; k++)
      Z[k] = list.at(k).toDouble();

    in.readLine(0); //Ze[NX+1]
    list = in.readLine(0).split(";");
    for(i = 0; i < imax+1; i++)
      Ze[i] = list.at(i).toDouble();

    in.readLine(0); //Body[NX][NZ]
    for(i = 0; i < imax; i++) {
        list = in.readLine(0).split(";");
        for(k = 0; k < kmax; k++)
          Body[i][k] = list.at(k).toDouble();
      }

    in.readLine(0); //U[NX][NY][NZ]
    for(i = 0; i < imax; i++)
      for(j = 0; j < jmax; j++) {
          list = in.readLine(0).split(";");
          for(k = 0; k < kmax; k++)
            U[i][j][k] = list.at(k).toDouble();
        }

    in.readLine(0); //W[NX][NY][NZ]
    for(i = 0; i < imax; i++)
      for(j = 0; j < jmax; j++) {
          list = in.readLine(0).split(";");
          for(k = 0; k < kmax; k++)
            W[i][j][k] = list.at(k).toDouble();
        }
    in.readLine(0); //H[NX][NY][NZ]
    for(i = 0; i < imax; i++)
      for(j = 0; j < jmax; j++) {
          list = in.readLine(0).split(";");
          for(k = 0; k < kmax; k++)
            H[i][j][k] = list.at(k).toDouble();
        }
    in.readLine(0); //V[NX][NY][NZ]
    for(i = 0; i < imax; i++)
      for(j = 0; j < jmax; j++) {
          list = in.readLine(0).split(";");
          for(k = 0; k < kmax; k++)
            V[i][j][k] = list.at(k).toDouble();
        }

    in.readLine(0); //P[NX+1][NZ]
    for(i = 0; i < imax+1; i++) {
        list = in.readLine(0).split(";");
        for(k = 0; k < kmax; k++)
          P[i][k] = list.at(k).toDouble();
      }
    in.readLine(0); //Delta[NX][NZ]
    for(i = 0; i < imax; i++) {
        list = in.readLine(0).split(";");
        for(k = 0; k < kmax; k++)
          Delta[i][k] = list.at(k).toDouble();
      }
    file.close();

    fillForm();
    Translate(N);
    Analysis();

    updateWing();
    colorMapTab->colorMap->setupPlot(PressureDistribution, "Pressure", false, false);
    colorMapTab->update();
    colorMapTab->setEnabled(true);
    colormapAct->setEnabled(true);

    QApplication::restoreOverrideCursor();
    statusBar()->showMessage(tr("File loaded"), 3000);
}
bool MainWindow::saveFile(const QString &fileName) {
    QFile file(fileName);
    file.open(QFile::WriteOnly | QFile::Text);

    QTextStream out(&file);
    //out.setFieldAlignment(QTextStream::AlignLeft);
    //out.setRealNumberNotation(QTextStream::FixedNotation);
    QApplication::setOverrideCursor(Qt::WaitCursor);

    out << "Theta0=" + QString::number(Theta0)
    + "\nTheta0d=" + QString::number(Theta0d)
    + "\nBeta=" + QString::number(Beta)
    + "\nBetad=" + QString::number(Betad)
    + "\nGamma=" + QString::number(Gamma)
    + "\nPr=" + QString::number(Pr)
    + "\nHw=" + QString::number(Hw)
    + "\nOmega=" + QString::number(Omega)
    + "\nrelU=" + QString::number(relU)
    + "\nrelP=" + QString::number(relP)
    + "\nAccuracyU=" + QString::number(AccuracyU)
    + "\nAccuracyPressure=" + QString::number(AccuracyPressure)
    + "\nimax=" + QString::number(imax)
    + "\ndx=" + QString::number(dx)
    + "\njmax=" + QString::number(jmax)
    + "\ndy=" + QString::number(dy)
    + "\nkmax=" + QString::number(kmax)
    + "\ndz=" + QString::number(dz);

    out << "\nX[NX]\n";
    for (i = 0; i < imax; i++) out << QString("%1;").arg(X[i]);
    out << "\nZ[NZ]\n";
    for (k = 0; k < kmax; k++) out << QString("%1;").arg(Z[k]);
    out << "\nZe[NX+1]\n";
    for (i = 0; i < imax+1; i++) out << QString("%1;").arg(Ze[i]);

    out << "\nBody[NX][NZ]\n";
    for (i = 0; i < imax; i++) {
        for (k = 0; k < kmax; k++)
            out << QString("%1;").arg(Body[i][k]);
        out << "\n";
    }

    out << "U[NX][NY][NZ]\n";
    for(i = 0; (i < imax) && ((imax > 2) || (i < 1)) ; i++)
        for(j = 0; j < jmax; j++) {
            for(k = 0; k < kmax; k++)
                out << QString("%1;").arg(U[i][j][k]);
            out << "\n";
        }
    out << "W[NX][NY][NZ]\n";
    for(i = 0; (i < imax) && ((imax > 2) || (i < 1)) ; i++)
        for(j = 0; j < jmax; j++) {
            for(k = 0; k < kmax; k++)
                out << QString("%1;").arg(W[i][j][k]);
            out << "\n";
        }
    out << "H[NX][NY][NZ]\n";
    for(i = 0; (i < imax) && ((imax > 2) || (i < 1)) ; i++)
        for(j = 0; j < jmax; j++) {
            for(k = 0; k < kmax; k++)
                out << QString("%1;").arg(H[i][j][k]);
            out << "\n";
        }
    out << "V[NX][NY][NZ]\n";
    for(i = 0; (i < imax) && ((imax > 2) || (i < 1)) ; i++)
        for(j = 0; j < jmax; j++) {
            for(k = 0; k < kmax; k++)
                out << QString("%1;").arg(V[i][j][k]);
            out << "\n";
        }
    out << "P[NX+1][NZ]\n";
    for (i = 0; i < imax+1; i++) {
        for (k = 0; k < kmax; k++)
            out << QString("%1;").arg(P[i][k]);
        out << "\n";
    }
    out << "Delta[NX][NZ]\n";
    for (i = 0; i < imax; i++) {
        for (k = 0; k < kmax; k++)
            out << QString("%1;").arg(Delta[i][k]);
        out << "\n";
    }
//	out << "TauU[NX][NZ]\n";
//	for(i = 0; i < NX; i++){
//		for(k = 0; k < NZ; k++)
//			out << QString("%1\n").arg(Tau[i][k]);
//		out << "\n";
//	}
//	out << "TauH[NX][NZ]\n";
//	for(i = 0; i < NX; i++){
//		for(k = 0; k < NZ; k++)
//			out << QString("%1\n").arg(TauH[i][k]);
//		out << "\n";
//	}
//	out << "P0[NX][NZ]\n";
//	for(i = 0; i < NX; i++){
//		for(k = 0; k < NZ; k++)
//			out << QString("%1\n").arg(P0[i][k]);
//		out << "\n";
//	}
//	out << "DeltaE[NX][NZ]\n";
//	for(i = 0; i < NX; i++){
//		for(k = 0; k < NZ; k++)
//			out << QString("%1\n").arg(DeltaE[i][k]);
//		out << "\n";
//	}
//	out << "PressureDistribution[N][N]\n";
//	for(i = 0; i < N; i++){
//		for(k = 0; k < N; k++)
//			out << QString("%1;").arg(PressureDistribution[i][k]);
//		out << "\n";
//	}

//   out << "Theta" << "Pressure" << "DeltaE" << "Tau" << "TauH";
//   out << "\n";
//   for(k = 0; k < NZ; k++) {out << Z[k] << P0[0][k] << DeltaE[0][k] << Tau[0][k] << TauH[0][k]; out << "\n";}
//   out << "\n";

    QApplication::restoreOverrideCursor();
    file.close();
    statusBar()->showMessage(tr("File saved"), 2000);
    return true;
}
void MainWindow::saveResults() {
    QString params = QString("%1-%2,Hw=%3,Vw=%4,Pr=%5,G=%6")
        .arg(Theta0d).arg(Betad).arg(Hw).arg(Vw, 0, 'f', 1).arg(Pr).arg(Gamma)/*.arg(imax).arg(jmax).arg(kmax)*/;
    QDir dir;
    dir.mkpath("results/"+params);
    //dir.cd("results/"+params);
    QFile file("results/"+params+"/TauU.dat");
    file.open(QFile::WriteOnly | QFile::Text);
    QTextStream out(&file);
    out.setRealNumberNotation(QTextStream::FixedNotation);
    out << "#Theta TauU\n";
    for(k = 0; k < kmax; k++) {
        out << QString("%1 %2\n").arg(Z[k]).arg(Tau[0][k]);
    }
    file.close();

    file.setFileName("results/"+params+"/TauH.dat");
    file.open(QFile::WriteOnly | QFile::Text);
    out << "#Theta TauH\n";
    for(k = 0; k < kmax; k++) {
        out << QString("%1 %2\n").arg(Z[k]).arg(TauH[0][k]);
    }
    file.close();

    file.setFileName("results/"+params+"/Pressure0.dat");
    file.open(QFile::WriteOnly | QFile::Text);
    out << "#Theta Pressure0\n";
    for(k = 0; k < kmax; k++) {
        out << QString("%1 %2\n").arg(Z[k]).arg(P0[0][k]);
    }
    file.close();

    file.setFileName("results/"+params+"/Pressure.dat");
    file.open(QFile::WriteOnly | QFile::Text);
    out << "#Theta Pressure\n";
    for(k = 0; k < kmax; k++) {
        out << QString("%1 %2\n").arg(Z[k]).arg(P[0][k]);
    }
    file.close();

    file.setFileName("results/"+params+"/DeltaE.dat");
    file.open(QFile::WriteOnly | QFile::Text);
    out << "#Theta DeltaE\n";
    for(k = 0; k < kmax; k++) {
        out << QString("%1 %2\n").arg(Z[k]).arg(DeltaE[0][k]);
    }
    file.close();

    file.setFileName("results/"+params+"/Delta.dat");
    file.open(QFile::WriteOnly | QFile::Text);
    out << "#Theta Delta\n";
    for(k = 0; k < kmax; k++) {
      out << QString("%1 %2\n").arg(Z[k]).arg(Delta[0][k]);
    }
    file.close();

    file.setFileName("results/"+params+"/ProfileU.dat");
    file.open(QFile::WriteOnly | QFile::Text);
    out << "#U Y\n";
    //int k = int((1+Beta/Theta0)*0.5);
    for(j = 0; j < jmax; j++) {
        out << j*dy;
        for(k = 0; k < kmax; k++) {
            out << QString("%1 %2\n").arg(sqrt(U[0][j][k]*U[0][j][k]+W[0][j][k]*W[0][j][k])).arg(j*dy);
        }
        out << "\n";
    }
    file.close();
    file.setFileName("results/"+params+"/ProfileH.dat");
    file.open(QFile::WriteOnly | QFile::Text);
    out << "#H Y\n";
    for(j = 0; j < jmax; j++) {
        out << j*dy;
        for(k = 0; k < kmax; k++) {
            out << H[0][j][k];
        }
        out << "\n";
    }
    file.close();

    file.setFileName("results/"+params+"/ProfilesU.dat");
    file.open(QFile::WriteOnly | QFile::Text);
    out << "#Y U\n";
    for(j = 0; j < jmax; j++) {
        out << j*dy;
        for(k = 0; k < kmax; k += kmax/10) {
            out << " " << U[0][j][k]+Z[k];
        }
        out << "\n";
    }
    file.close();

    file.setFileName("results/"+params+"/ProfilesW.dat");
    file.open(QFile::WriteOnly | QFile::Text);
    out << "#Y U\n";
    for(j = 0; j < jmax; j++) {
        out << j*dy;
        for(k = 0; k < kmax; k += kmax/10) {
            out << " " << W[0][j][k]+Z[k];
        }
        out << "\n";
    }
    file.close();

    file.setFileName("results/"+params+"/diagram.dat");
    file.open(QFile::WriteOnly | QFile::Text);
    out << "#angle velocity\n";
    for (int angle = 0; angle <= 359; angle++) {
        out << angle-90; //-90 to rotate polar plot in gnuplot
        for (k = 0; k < kmax; k++)
            out << " " << DisturbancesDistribution[angle][k];
        out << "\n";
    }
    file.close();

    file.setFileName("results/"+params+"/a[180][k].dat");
    file.open(QFile::WriteOnly | QFile::Text);
    out << "#Theta a[180]\n";
    for (k = 0; k < kmax; k++) {
        out << Z[k] << " " << DisturbancesDistribution[180][k] << "\n";
    }
    file.close();

    file.setFileName("results/"+params+"/PressureDistribution.dat");
    file.open(QFile::WriteOnly | QFile::Text);
    for (int m = 0; m < N; m+=4)
         for (int l = 0; l < N; l+=4) {
             out << m << " " << l << " " << PressureDistribution[m][l] << "\n";
         }
    file.close();

    out.setRealNumberNotation(QTextStream::SmartNotation);
    params = QString("%1-%2,Hw=%3,Pr=%4,G=%5")
        .arg(Theta0d).arg(Betad).arg(Hw).arg(Pr).arg(Gamma)/*.arg(imax).arg(jmax).arg(kmax)*/;
    file.setFileName("results/Vw_influence_"+params+".dat");
    file.open(QFile::Append | QFile::Text);
    out << Vw << "\t" << DisturbancesDistribution[180][kmax/2+int(Betad/Theta0d*kmax/2)] << "\t" << Delta99[0][kmax/2+int(Betad/Theta0d*kmax/2)] << "\t" << QString("Parameters_%1-%2,Hw=%3,Vw=%4,Pr=%5,G=%6,imax=%7,jmax=%8,kmax=%9").arg(Theta0d).arg(Betad).arg(Hw).arg(Vw, 0, 'f', 1).arg(Pr).arg(Gamma).arg(imax).arg(jmax).arg(kmax) << "\n";
    file.close();
}
void MainWindow::getParameters() {
    Theta0d = lineEdit_geometryTheta0d->text().toDouble();
    Theta0 = Theta0d*Pi/180;
    Theta1d = lineEdit_geometryTheta1d->text().toDouble();
    Theta1 = Theta1d*Pi/180;
    Betad = lineEdit_geometryBetad->text().toDouble();
    Beta = Betad*Pi/180;
    h = lineEdit_geometryH->text().toDouble();
    L = lineEdit_geometryL->text().toDouble();
    Gamma = lineEdit_parametersGamma->text().toDouble();
    Pr = lineEdit_parametersPr->text().toDouble();
    Hw = lineEdit_parametersHw->text().toDouble();
    Vw = lineEdit_parametersOmega->text().toDouble();
    relU = lineEdit_solving_relU->text().toDouble();
    relP = lineEdit_solving_relP->text().toDouble();
    MaxIterations = lineEdit_solving3->text().toInt();
    relP2 = lineEdit_solving4->text().toDouble();
    AccuracyU = lineEdit_solving_AcU->text().toDouble();
    AccuracyPressure = lineEdit_solving_AcPr->text().toDouble();
    dx = lineEdit_grid_dx->text().toDouble();
    imax = lineEdit_grid_NX->text().toInt();
    dy = lineEdit_grid_dy->text().toDouble();
    jmax = lineEdit_grid_NY->text().toInt();
    dz = lineEdit_grid_dz->text().toDouble();
    kmax = lineEdit_grid_NZ->text().toInt();
    graphicTab->spinBoxLX->setRange(0, imax-1);
    graphicTab->spinBoxLY->setRange(0, jmax-1);
    graphicTab->spinBoxLZ->setRange(0, kmax-1);
}
void MainWindow::solve() {
    QString problem = "3D"; // Symmetry, SelfSimilar, 3D, Disturbances, EigenValues
    enableUI(false);
    resetProgress();
    DeleteArrays(checkBox_MultiGrid->isChecked());
    getParameters();
    CreateArrays(checkBox_MultiGrid->isChecked(), imax, jmax, kmax, dx, dz, N);
    SetBoundaryConditions(radioC2->isChecked(), lineEdit_parametersHw->text().toDouble());
    timeFunc.start();
    tabWidget->setCurrentIndex(1); //Show plots
    QList<double> paramslist;

    SetZe(imax, dx, Theta0);

    if (problem == "Symmetry") {
        /*Hw*/ //for (double x = 0.05; x < 1.0; x+=0.05) params << x;
        /*Pr*/ //params << 0.7 << 0.8 << 0.9 << 1.0;
        /*Gamma*/ //params << 1.1 << 1.15 << 1.2 << 1.25 << 1.3 << 1.35 << 1.4;
         paramslist << 0.1 << 0.9;
        foreach(double x, paramslist) {
            Hw = x;
            qDebug() << "Hw =" << x;
            Simm();
            for(j = 0; j < jmax; j++) {
                qDebug() << U[0][j][0];
            }
            for(j = 0; j < jmax; j++) {
                qDebug() << H[0][j][0];
            }
        }
    }
    else if (problem == "SelfSimilar") {
        for(i = 0; i < imax/L; i++) {
            Edge(Theta0, i, radioC2->isChecked(), -1);
            Edge(Theta0, i, radioC2->isChecked(), 1);
        }
        type = "Nose";
        Nos(radioC2->isChecked());
        Translate(N);
        Analysis();
    }
    else if (problem == "EigenValues") {
        SelfNumbers(-20, 1, 0.5);
        /*z0*/ //params << 2 << 2 << 2 << 2 << 1 << 0.75 << 0.5 << 0.25 << 0.1 << 0.01;
        /*Gamma*/ //params << 1.6 << 1.5 << 1.4 << 1.3 << 1.2 << 1.1 << 1.05;
        /*P2*/ paramslist << -50 << -40 << -30 << -20 << -15 << -10 << -9 << -8 << -7 << -6 << -5 << -3 << -2 << -1 << -0.5 << -0.2;
    //    foreach(double x, params) {
    //        P2 = x;
    //        qDebug() << x << SelfNumbers(P2, z0, Hw) << Delta[1][0] << Delta[0][0];
    //    }
    //    qDebug() << SelfNumbers(P2, z0, Hw) << Delta[1][0] << Delta[0][0];
    }
    else if (problem == "Disturbances") {
        //paramslist << -2;//<< -5 << -4 << -3 << -2 << -1 << -0.5 << -0.2 << -0.1 << 0 << 0.1 << 0.2 << 0.3 << 0.4 << 0.5 << 0.6 << 0.7 << 0.8 << 0.9 << 1.0;
        //paramslist << 0.05 << 0.1 << 0.2 << 0.3 << 0.4 << 0.5 << 0.6 << 0.7 << 0.8 << 0.9 << 0.95;
        paramslist << -1.0 << 0.0 << 0.2;
        for (Hw=0.1; Hw<=0.9; Hw+=0.4) {
        foreach(double x, paramslist) {
           Vw = x;
           SetBoundaryConditions(radioC2->isChecked(), Hw/*lineEdit_parametersHw->text().toDouble()*/);
           qDebug() << "Hw =" << Hw << ", Vw =" << Vw;
           for(i = 0; i < imax/L; i++) {
               Edge(Theta0, i, radioC2->isChecked(), -1);
               Edge(Theta0, i, radioC2->isChecked(), 1);
           }
           type = "Nose";
           Nos(radioC2->isChecked());
           Translate(N);
           Analysis();
           Disturbances(180, 180);
           //save results to file
           QString params = QString("%1-%2,Hw=%3,Vw=%4,Pr=%5,G=%6")
               .arg(Theta0d).arg(Betad).arg(Hw).arg(Vw, 0, 'f', 1).arg(Pr).arg(Gamma)/*.arg(imax).arg(jmax).arg(kmax)*/;
           QFile file("results/Disturbs_"+params+".dat");
           file.open(QFile::WriteOnly | QFile::Text);
           QTextStream out(&file);
           //out.setRealNumberNotation(QTextStream::FixedNotation);
           out.setFieldWidth(10);
           //out.setFieldAlignment(QTextStream::AlignLeft);
           out.setRealNumberPrecision(4);
           //file.setFileName("results/Pressure_"+params+".dat");
           out << "#Theta" << "Pressure" << "Delta" << "Mach" << "Tau" << "a[180]" << "Delta99" << "Subsonic" << "\n";
           for(k = 0; k < kmax; k++) {
               out << Z[k] << P[0][k] << Delta[0][k] << J[k] << Tau[0][k] << DisturbancesDistribution[180][k] << Delta99[0][k] << Delta99[0][k]/Delta[0][k] << "\n";
           }
           file.close();
           saveResults();
       }
    }
    }
    else if (problem == "3D") {
        for(i = 0; i < imax/L; i++) {
            Edge(Theta0, i, radioC2->isChecked(), -1);
            Edge(Theta0, i, radioC2->isChecked(), 1);
        }
        type = "Nose";
        Nos(radioC2->isChecked());
        type = "Wing";
        timeFunc.restart();
        if (imax > 2) FullWing(radioC2->isChecked(), checkBox_MultiGrid->isChecked());
        Translate(N);
        Analysis();
    }
    else {
        if(checkBox_MultiGrid->isChecked()) {
            DeleteArrays(checkBox_MultiGrid->isChecked());
            CreateArrays(checkBox_MultiGrid->isChecked(),
                         lineEdit_grid_NX->text().toInt(),
                         lineEdit_grid_NY->text().toInt(),
                         lineEdit_grid_NZ->text().toInt(),
                         lineEdit_grid_dx->text().toDouble(),
                         lineEdit_grid_dz->text().toDouble(), N);
            InterpolateArrays(lineEdit_grid_NX->text().toInt(),
                              lineEdit_grid_NY->text().toInt(),
                              lineEdit_grid_NZ->text().toInt(),
                              imax, jmax, kmax);
            getParameters();
            SetBoundaryConditions(radioC2->isChecked(), lineEdit_parametersHw->text().toDouble());
        }
        else {
            for(i = 0; i < imax/L; i++) Edge(Theta0, i, radioC2->isChecked(), -1);
            //if((Theta1 !=0) && (Theta1 != Theta0)) Edge(Theta1, NX/2, radioC2->isChecked(), -1);
            label_results_Edge1_1->setText(QString::number(NIterations));
            label_results_Edge1_2->setText(QString::number(timeFunc.elapsed()/1000.));
            timeFunc.restart();

            for(i = 0; i < imax/L; i++) Edge(Theta0, i, radioC2->isChecked(), 1);
            //if((Theta1 !=0) && (Theta1 != Theta0)) Edge(Theta1, NX/2, radioC2->isChecked(), 1);
            label_results_Edge2_1->setText(QString::number(NIterations));
            label_results_Edge2_2->setText(QString::number(timeFunc.elapsed()/1000.));

            type = "Nose";
            timeFunc.restart();
            //Nos(radioC2->isChecked());

            label_results_Nose_1->setText(QString::number(NIterations));
        }
        type = "Wing";
        timeFunc.restart();
        if (imax > 2) FullWing(radioC2->isChecked(), checkBox_MultiGrid->isChecked());
        //	if ((NX > 2) && (L >= 2.0)) FullWingWithTrail(radioC2->isChecked());
        Translate(N);
        Analysis();
        //Disturbances(0, 359); //(0, 359) for a full diagramm
    }
    //saveResults();
    colorMapTab->colorMap->setupPlot(TauUDistribution, "TauU", false, true);
    colorMapTab->update();
    graphicTab->renderarea->update();

    finishProgress();
    enableUI(true);
}
void MainWindow::solveCycle() {
    QTimer *timer = new QTimer;
    connect(timer, SIGNAL(timeout()), this, SLOT(updateProgress()));
    if(checkBox_graph->isChecked())
      connect(timer, SIGNAL(timeout()), graphicTab->renderarea, SLOT(update()));
    timer->start(2000);

    DeleteArrays(false);
    getParameters();
    CreateArrays(false, imax, jmax, kmax, dx, dz, N);
    //Plate();
    QList<double> params;
    params << -5 << -4 << -3 << -2 << -1 << -0.5 << -0.25 << -0.1 << 0 << 0.1 << 0.2 << 0.3 << 0.4 << 0.5 << 0.6 << 0.7 << 0.8 << 0.9 << 1.0;
//    params << -5 <<  -3 << -1 << -0.5 << -0.25 << 0 << 0.2 << 0.4 << 0.6 << 0.8 << 1.0;
    Theta0d = 45;
    Betad = 0;
    for (Hw = 0.1; Hw <= 0.9; Hw += 0.4) {
        qDebug() << QString("Theta=%1,Beta=%2,Hw=%3").arg(Theta0d).arg(Betad).arg(Hw);
        foreach(double x, params) {
            Vw = x;
            fillForm();
            solve();
            if(status == "bad") {
                status = "ok";
                break;
              }
          }
    }
    Theta0d = 105;
    Betad = 35;
    for (Hw = 0.1; Hw <= 0.9; Hw += 0.4) {
        qDebug() << QString("Theta=%1,Beta=%2,Hw=%3").arg(Theta0d).arg(Betad).arg(Hw);
        foreach(double x, params) {
            Vw = x;
            fillForm();
            solve();
            if(status == "bad") {
                status = "ok";
                break;
              }
          }
    }
    Theta0d = 90;
    Betad = 0;
    for (Hw = 0.1; Hw <= 0.9; Hw += 0.4) {
        qDebug() << QString("Theta=%1,Beta=%2,Hw=%3").arg(Theta0d).arg(Betad).arg(Hw);
        foreach(double x, params) {
            Vw = x;
            fillForm();
            solve();
            if(status == "bad") {
                status = "ok";
                break;
              }
          }
    }
}
void MainWindow::solveProblem(QString filename) {
    QFile file(filename);
    if(!file.open(QFile::ReadOnly | QFile::Text)) {
          QMessageBox::warning(this, "HS", tr("Cannot read file %1:\n%2.")
                                           .arg(filename)
                                           .arg(file.errorString()));
          return;
    }
    QTextStream in(&file);
    QString line;
    while(!in.atEnd()) {
        line = in.readLine(0);
        if(line.contains("="))
          readParameters(line);
        else if(line.startsWith("#Case"))
          break;
        else
          continue;
    }
    dx = L/(imax-1);
    dy = h/(jmax-1);
    dz = 2.0/(kmax-1);
    while(!in.atEnd()) {
        line = in.readLine(0);
        if(line.contains("="))
            readParameters(line);
        else if (line.startsWith("#Case")) {
            DeleteArrays(false);
            CreateArrays(false, imax, jmax, kmax, dx, dz, N);
            fillForm();
            updateWing();
            solve();
        }
    }
    file.close();
}
void MainWindow::readParameters(QString line) {
  if(line.split("=").at(0).toUtf8() == "Theta0d") {
    Theta0d = line.split("=").at(1).toDouble();
    Theta0 = Theta0d*Pi/180;
  }
  else if(line.split("=").at(0).toUtf8() == "Betad") {
    Betad = line.split("=").at(1).toDouble();
    Beta = Betad*Pi/180;
  }
  else if(line.split("=").at(0).toUtf8() == "Gamma")
    Gamma = line.split("=").at(1).toDouble();
  else if(line.split("=").at(0).toUtf8() == "Pr")
    Pr = line.split("=").at(1).toDouble();
  else if(line.split("=").at(0).toUtf8() == "Hw")
    Hw = line.split("=").at(1).toDouble();
  else if(line.split("=").at(0).toUtf8() == "Vw")
    Vw = line.split("=").at(1).toDouble();
  else if(line.split("=").at(0).toUtf8() == "Omega")
    Omega = line.split("=").at(1).toDouble();
  else if(line.split("=").at(0).toUtf8() == "relU")
    relU = line.split("=").at(1).toDouble();
  else if(line.split("=").at(0).toUtf8() == "relP")
    relP = line.split("=").at(1).toDouble();
  else if(line.split("=").at(0).toUtf8() == "AccuracyU")
    AccuracyU = line.split("=").at(1).toDouble();
  else if(line.split("=").at(0).toUtf8() == "AccuracyPressure")
    AccuracyPressure = line.split("=").at(1).toDouble();
  else if(line.split("=").at(0).toUtf8() == "h")
    h = line.split("=").at(1).toDouble();
  else if(line.split("=").at(0).toUtf8() == "L")
    L = line.split("=").at(1).toDouble();
  else if(line.split("=").at(0).toUtf8() == "dt")
    dt = line.split("=").at(1).toDouble();
  else if(line.split("=").at(0).toUtf8() == "imax")
    imax = line.split("=").at(1).toInt();
  else if(line.split("=").at(0).toUtf8() == "jmax")
    jmax = line.split("=").at(1).toInt();
  else if(line.split("=").at(0).toUtf8() == "kmax")
    kmax = line.split("=").at(1).toInt();
  else if(line.split("=").at(0).toUtf8() == "NCases")
    NCases = line.split("=").at(1).toInt();
  else
    QMessageBox::warning(this, "HS", tr("Cannot read parameter \"%1\" from file.")
                            .arg(line));
}
void MainWindow::selfTest() {
    setWindowTitle(tr("SelfTest running... HS"));
    enableUI(false);
    solveAct->setEnabled(false);
    stopAct->setEnabled(true);
    QMessageBox::information(this, tr("SelfTest"),SelfTest());
    graphicTab->renderarea->update();
    solveAct->setEnabled(true);
    stopAct->setEnabled(false);
    enableUI(true);
    setWindowTitle("HS");
}
void MainWindow::updateWing() {
    getParameters();
    SetZe(imax, dx, Theta0);
    wingTab->wingWidget->update();
}
void MainWindow::resetProgress() {
    label_results_Edge1_1->setText(tr("0"));
    label_results_Edge1_2->setText(tr("0 s"));
    label_results_Edge1_3->setText(tr("0 ips"));
    label_results_Edge2_1->setText(tr("0"));
    label_results_Edge2_2->setText(tr("0 s"));
    label_results_Edge2_3->setText(tr("0 ips"));
    label_results_Nose_1->setText(tr("0"));
    label_results_Nose_2->setText(tr("0 s"));
    label_results_Nose_3->setText(tr("0 ips"));
    label_results_Wing_1->setText(tr("0"));
    label_results_Wing_2->setText(tr("0 s"));
    label_results_Wing_3->setText(tr("0 ips"));
    label_results_TotalTime1->setText(tr("0 s"));
    label_results_FinishTime0->setText(tr(""));
    label_results_FinishTime1->setText(tr(""));

    timeTotal.start();
    timer = new QTimer;
    connect(timer, SIGNAL(timeout()), this, SLOT(updateProgress()));
    if(checkBox_graph->isChecked())
        connect(timer, SIGNAL(timeout()), graphicTab->renderarea, SLOT(update()));
    timer->start(2000);
}
void MainWindow::updateProgress() {
    double t = timeFunc.elapsed()/1000.;
    if(type == "Nose") {
        label_results_Nose_1->setText(tr("<b>%1</b>").arg(QString::number(NIterations)));
        label_results_Nose_2->setText(QString::number(t));
        label_results_Nose_3->setText(QString::number(NIterations/t, 'f', 1));
    }
    else {
        label_results_Wing_1->setText(tr("<b>%1</b>").arg(QString::number(NIterations)));
        label_results_Wing_2->setText(QString::number(t));
        label_results_Wing_3->setText(QString::number(NIterations/t, 'f', 1));
    }
    label_results_TotalTime1->setText(QString::number(timeTotal.elapsed()/1000.)+" s");
}
void MainWindow::finishProgress() {
    timer->stop();

    label_results_Wing_1->setText(QString::number(NIterations));
    label_results_TotalTime1->setText(QString::number(timeTotal.elapsed()/1000.)+" s");

    label_results_FinishTime0->setText(tr("Finished at:"));
    label_results_FinishTime1->setText(QDateTime::currentDateTime().toString("HH:mm:ss"));
}
void MainWindow::readSettings() {
    QSettings settings("settings.ini", QSettings::IniFormat);
    //this->resize(settings.value("MainWindow/size", QSize(800, 600)).toSize());
    //this->move(settings.value("MainWindow/pos", QPoint(100, 100)).toPoint());
    if(settings.value("MainWindow/isMaximized", 1).toBool()) this->showMaximized();
    else this->resize(QSize(1000, 800));
    if(settings.value("Language", "en").toString() == "ru") {
        translator.load(QString(":/translations/hs_ru"));
        qApp->installTranslator(&translator);
    }
    checkBox_graph->setChecked(settings.value("Solving/updatePlot",1).toBool());
}
void MainWindow::writeSettings() {
    QSettings settings("settings.ini", QSettings::IniFormat);
    settings.setValue("MainWindow/size", this->size());
    settings.setValue("MainWindow/pos", this->pos());
    settings.setValue("MainWindow/isMaximized", isMaximized());
    settings.setValue("Solving/updatePlot", checkBox_graph->isChecked());
    //settings.setValue("Language", "en");
}
void MainWindow::enableUI(bool is_enabled) {
    colorMapTab->setEnabled(is_enabled);
    colormapAct->setEnabled(is_enabled);
    saveAct->setEnabled(is_enabled);
    saveAsAct->setEnabled(is_enabled);
    solveAct->setEnabled(is_enabled);
    stopAct->setEnabled(!is_enabled);
    if(is_enabled) {
        setWindowTitle("HS");
        QApplication::restoreOverrideCursor();
    }
    else {
        setWindowTitle("Solving...");
        QApplication::setOverrideCursor(Qt::WaitCursor);
    }
}
WingTab::WingTab(QWidget *parent)
     : QWidget(parent)  {
    wingWidget = new WingView;
    QToolBar *wingToolBar = new QToolBar(tr("Wing"), this);
    wingToolBar->setIconSize(QSize(15,15));
    gridCheckBox = new QCheckBox("Draw grid", this);
    gridCheckBox->setChecked(true);
    connect(gridCheckBox, SIGNAL(clicked(bool)), SLOT(drawGrid(bool)));
    wingToolBar->addWidget(gridCheckBox);
    QVBoxLayout *layout = new QVBoxLayout;
    layout->addWidget(wingToolBar);
    layout->addWidget(wingWidget);
    layout->setMargin(2);
    setLayout(layout);
 }
GraphicTab::GraphicTab(QWidget *parent)
     : QWidget(parent) {

    QLineSeries *series1 = new QLineSeries();
    *series1 << QPointF(1, 1) << QPointF(2, 73) << QPointF(3, 268) << QPointF(4, 17) << QPointF(5, 4325) << QPointF(6, 723);

    chart = new QChart();
    chart->addSeries(series1);
    //chart->legend()->hide();
    chart->setTitle("Example title");

    axisX = new QValueAxis;
    axisX->setTitleText("Data point");
    axisX->setTickCount(6);
    axisX->setLabelFormat("%i");
    chart->addAxis(axisX, Qt::AlignBottom);
    series1->attachAxis(axisX);

    axisY = new QValueAxis;
    axisY->setLabelFormat("%g");
    axisY->setTitleText("Values");
    chart->addAxis(axisY, Qt::AlignLeft);
    series1->attachAxis(axisY);

    chartView = new QChartView(chart);
    chartView->setRenderHint(QPainter::Antialiasing);

    QToolBar *plotToolBar = new QToolBar(tr("Graphic"), this);
    plotToolBar->setIconSize(QSize(15,15));

    QFont font = QFont("Arial",9);
    QPushButton *fxButton = new QPushButton("f(x)");
    fxButton->setFont(font); fxButton->setFixedSize(30,20);
    connect(fxButton, SIGNAL(clicked()), this, SLOT(Fx()));
    QPushButton *fyButton = new QPushButton("f(y)");
    fyButton->setFont(font); fyButton->setFixedSize(30,20);
    connect(fyButton, SIGNAL(clicked()), this, SLOT(Fy()));
    QPushButton *fzButton = new QPushButton("f(z)");
    fzButton->setFont(font); fzButton->setFixedSize(30,20);
    connect(fzButton, SIGNAL(clicked()), this, SLOT(Fz()));
    QPushButton *tauButton = new QPushButton("Tau");
    tauButton->setFont(font); tauButton->setFixedSize(30,20);
    //connect(tauButton, SIGNAL(clicked()), renderarea, SLOT(PaintTau()));
    QPushButton *pButton = new QPushButton("P");
    pButton->setFont(font); pButton->setFixedSize(30,20);
    //connect(pButton, SIGNAL(clicked()), renderarea, SLOT(PaintP()));
    QPushButton *deltaButton = new QPushButton(tr("%1").arg(QChar(0x0394)));
    deltaButton->setFont(font); deltaButton->setFixedSize(30,20);
    //connect(deltaButton, SIGNAL(clicked()), renderarea, SLOT(PaintDelta()));
    QPushButton *disturbButton = new QPushButton(tr("Disturbs"));
    disturbButton->setFont(font); disturbButton->setFixedSize(60,20);
    //connect(disturbButton, SIGNAL(clicked()), renderarea, SLOT(PaintDisturbances()));
    QPushButton *fromFileButton = new QPushButton(tr("From file..."));
    fromFileButton->setFont(font); fromFileButton->setFixedSize(80,20);
    //connect(fromFileButton, SIGNAL(clicked()), this, SLOT(loadData()));
    fileComboBox = new QComboBox(this);
    //fileComboBox->addItems();
    //connect(fileComboBox, SIGNAL(currentIndexChanged(int)), renderarea, SLOT(PaintData(int)));

    QLabel *label_LX = new QLabel("  LX:");
    spinBoxLX = new QSpinBox(this);
    spinBoxLX->setRange(0, imax-1);
    spinBoxLX->setValue(LX);
    connect(spinBoxLX, SIGNAL(valueChanged(int)), this, SLOT(changeLX(int)));
    QLabel *label_LY = new QLabel("  LY:");
    spinBoxLY = new QSpinBox(this);
    spinBoxLY->setRange(0, jmax-1);
    spinBoxLY->setValue(LY);
    connect(spinBoxLY, SIGNAL(valueChanged(int)), this, SLOT(changeLY(int)));
    QLabel *label_LZ = new QLabel("  LZ:");
    spinBoxLZ = new QSpinBox(this);
    spinBoxLZ->setRange(0, kmax-1);
    spinBoxLZ->setValue(LZ);
    connect(spinBoxLZ, SIGNAL(valueChanged(int)), this, SLOT(changeLZ(int)));

    QLabel *labelMinX = new QLabel(" MinX:");
    spinBoxMinX = new QDoubleSpinBox(this);
    spinBoxMinX->setDecimals(2);
    spinBoxMinX->setSingleStep(0.5);
    spinBoxMinX->setRange(-1000, 0);
    spinBoxMinX->setValue(renderarea->minX);
    spinBoxMinX->setFixedWidth(50);
    connect(spinBoxMinX, SIGNAL(valueChanged(double)), this, SLOT(changeMinX(double)));
    QLabel *labelMaxX = new QLabel(" MaxX:");
    spinBoxMaxX = new QDoubleSpinBox(this);
    spinBoxMaxX->setDecimals(2);
    spinBoxMaxX->setSingleStep(0.5);
    spinBoxMaxX->setRange(0, 1000);
    spinBoxMaxX->setValue(renderarea->maxX);
    spinBoxMaxX->setFixedWidth(50);
    connect(spinBoxMaxX, SIGNAL(valueChanged(double)), this, SLOT(changeMaxX(double)));
    QLabel *labelMinY = new QLabel(" MinY:");
    spinBoxMinY = new QDoubleSpinBox(this);
    spinBoxMinY->setDecimals(2);
    spinBoxMinY->setSingleStep(0.5);
    spinBoxMinY->setRange(-100, 0);
    spinBoxMinY->setValue(renderarea->minY);
    spinBoxMinY->setFixedWidth(50);
    connect(spinBoxMinY, SIGNAL(valueChanged(double)), this, SLOT(changeMinY(double)));
    QLabel *labelMaxY = new QLabel(" MaxY:");
    spinBoxMaxY = new QDoubleSpinBox(this);
    spinBoxMaxY->setDecimals(2);
    spinBoxMaxY->setSingleStep(0.5);
    spinBoxMaxY->setRange(0, 1000);
    spinBoxMaxY->setValue(renderarea->maxY);
    spinBoxMaxY->setFixedWidth(50);
    connect(spinBoxMaxY, SIGNAL(valueChanged(double)), this, SLOT(changeMaxY(double)));

    plotToolBar->addWidget(fxButton);
    plotToolBar->addWidget(fyButton);
    plotToolBar->addWidget(fzButton);
    plotToolBar->addWidget(tauButton);
    plotToolBar->addWidget(pButton);
    plotToolBar->addWidget(deltaButton);
    plotToolBar->addWidget(disturbButton);
    plotToolBar->addWidget(fromFileButton);
    plotToolBar->addSeparator();
    plotToolBar->addWidget(fileComboBox);
    plotToolBar->addSeparator();
    plotToolBar->addWidget(label_LX);
    plotToolBar->addWidget(spinBoxLX);
    plotToolBar->addWidget(label_LY);
    plotToolBar->addWidget(spinBoxLY);
    plotToolBar->addWidget(label_LZ);
    plotToolBar->addWidget(spinBoxLZ);
    plotToolBar->addSeparator();
    plotToolBar->addSeparator();
    plotToolBar->addWidget(labelMinX);
    plotToolBar->addWidget(spinBoxMinX);
    plotToolBar->addWidget(labelMaxX);
    plotToolBar->addWidget(spinBoxMaxX);
    plotToolBar->addWidget(labelMinY);
    plotToolBar->addWidget(spinBoxMinY);
    plotToolBar->addWidget(labelMaxY);
    plotToolBar->addWidget(spinBoxMaxY);


    QVBoxLayout *layout = new QVBoxLayout;
    layout->addWidget(plotToolBar);
    layout->addWidget(chartView);
    layout->setMargin(0);
    setLayout(layout);
}
void GraphicTab::loadData() {
    QString fileName = QFileDialog::getOpenFileName(this);
    if (fileName.isEmpty()) return;
    QFile file(fileName);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, "HS",
                             tr("Cannot read file %1:\n%2.")
                             .arg(fileName)
                             .arg(file.errorString()));
        return;
    }
    QApplication::setOverrideCursor(Qt::WaitCursor);
    QTextStream in(&file);
    QStringList list;
    renderarea->Ncalc  = in.readLine(0).split("=").at(1).toInt();
    renderarea->Ndata  = in.readLine(0).split("=").at(1).toInt();
    renderarea->NXd = in.readLine(0).split("=").at(1).toInt();
    renderarea->NYd = in.readLine(0).split("=").at(1).toInt();
    renderarea->NZd = in.readLine(0).split("=").at(1).toInt();

    renderarea->data.array = new double**[renderarea->Ncalc];
    for (int n = 0; n < renderarea->Ncalc; n++) {
        renderarea->data.array[n] = new double*[renderarea->Ndata];
        for (int j = 0; j < renderarea->Ndata; j++)
            renderarea->data.array[n][j] = new double [renderarea->NZd];
    }
    for(int n = 0; n < renderarea->Ncalc; n++) {
        renderarea->data.title << in.readLine(0);
        for(int j = 0; j < renderarea->Ndata; j++) {
            list = in.readLine(0).split(";");
            if(renderarea->data.caption.count() < renderarea->Ndata)
                renderarea->data.caption << list.at(0);
            for(int k = 0; k < renderarea->NZd; k++)
                renderarea->data.array[n][j][k] = list.at(k+1).toDouble();
        }
    }
    qDebug() << "Omega=0.7";
    for(int k = 0; k < renderarea->NZd; k++)
        qDebug() << renderarea->data.array[0][3][k];
    qDebug() << "Omega=1.0";
    for(int k = 0; k < renderarea->NZd; k++)
        qDebug() << renderarea->data.array[1][3][k];
    QApplication::restoreOverrideCursor();
    renderarea->Narr = 0;
    this->fileComboBox->addItems(renderarea->data.caption);
}
ArrayTab::ArrayTab(QWidget *parent)
     : QWidget(parent) {
    ViewArray *arrayView = new ViewArray(this);

    QToolBar *toolBar = new QToolBar(tr("Array"), this);
    toolBar->setIconSize(QSize(15,15));

    QFont font = QFont("Times",12,QFont::Bold);
    QPushButton *buttonU = new QPushButton("&U");
    buttonU->setToolTip(tr("Show array U"));
    buttonU->setFont(font);	buttonU->setFixedSize(40,25);
    connect(buttonU, SIGNAL(clicked()), arrayView, SLOT(showU()));
    toolBar->addWidget(buttonU);
    QPushButton *buttonW = new QPushButton("&W");
    buttonW->setToolTip(tr("Show array W"));
    buttonW->setFont(font); buttonW->setFixedSize(40,25);
    connect(buttonW, SIGNAL(clicked()), arrayView, SLOT(showW()));
    toolBar->addWidget(buttonW);
    QPushButton *buttonH = new QPushButton("&H");
    buttonH->setToolTip(tr("Show array H"));
    buttonH->setFont(font); buttonH->setFixedSize(40,25);
    connect(buttonH, SIGNAL(clicked()), arrayView, SLOT(showH()));
    toolBar->addWidget(buttonH);
    QPushButton *buttonV = new QPushButton("&V");
    buttonV->setToolTip(tr("Show array V"));
    buttonV->setFont(font); buttonV->setFixedSize(40,25);
    connect(buttonV, SIGNAL(clicked()), arrayView, SLOT(showV()));
    toolBar->addWidget(buttonV);
    QPushButton *buttonT = new QPushButton("&T");
    buttonT->setToolTip(tr("Show array T"));
    buttonT->setFont(font); buttonT->setFixedSize(40,25);
    connect(buttonT, SIGNAL(clicked()), arrayView, SLOT(showT()));
    toolBar->addWidget(buttonT);
    QPushButton *buttonOther = new QPushButton(tr("&Other"));
    buttonOther->setToolTip(tr("Show other arrays"));
    buttonOther->setFont(font); buttonOther->setFixedSize(55,25);
    connect(buttonOther, SIGNAL(clicked()), arrayView, SLOT(showOther()));
    toolBar->addWidget(buttonOther);
    toolBar->addSeparator();
    QLabel *labelLX = new QLabel("  LX:");
    toolBar->addWidget(labelLX);
    QSpinBox *spinBoxLX = new QSpinBox;
    spinBoxLX = new QSpinBox(this);
    spinBoxLX->setRange(0, imax-1);
    spinBoxLX->setValue(LX);
    toolBar->addWidget(spinBoxLX);

    connect(spinBoxLX, SIGNAL(valueChanged(int)), this, SLOT(changeLX(int)));

    QVBoxLayout *layout = new QVBoxLayout;
    layout->addWidget(toolBar);
    layout->addWidget(arrayView);
    layout->setMargin(0);
    setLayout(layout);
 }
ColorMapTab::ColorMapTab(QWidget *parent)
     : QWidget(parent) {
    colorMap = new ColorMapWidget(this);

    QToolBar *colorMapToolBar = new QToolBar(tr("Graphic"), this);
    colorMapToolBar->setIconSize(QSize(15,15));

    QFont font = QFont("Times",12,QFont::Bold);

    QPushButton *TauUButton = new QPushButton("&TauU");
    TauUButton->setFixedSize(45,25); TauUButton->setFont(font);
    connect(TauUButton, SIGNAL(clicked()),this, SLOT(drawTauU()));
    QPushButton *TauWButton = new QPushButton("&TauW");
    TauWButton->setFixedSize(45,25); TauWButton->setFont(font);
    connect(TauWButton, SIGNAL(clicked()),this, SLOT(drawTauW()));
    QPushButton *TauHButton = new QPushButton("&TauH");
    TauHButton->setFixedSize(45,25); TauHButton->setFont(font);
    connect(TauHButton, SIGNAL(clicked()),this, SLOT(drawTauH()));
    QPushButton *PressureButton = new QPushButton("&P");
    PressureButton->setFixedSize(45,25); PressureButton->setFont(font);
    connect(PressureButton, SIGNAL(clicked()),this, SLOT(drawPressure()));
    QPushButton *DeltaButton = new QPushButton(tr("%1").arg(QChar(0x0394)));
    DeltaButton->setFixedSize(45,25); DeltaButton->setFont(font);
    connect(DeltaButton, SIGNAL(clicked()),this, SLOT(drawDelta()));
    chBox_flowlines = new QCheckBox(tr("Flowlines"));
    chBox_flowlines->setChecked(false);
    connect(chBox_flowlines, SIGNAL(toggled(bool)), this, SLOT(changeFlowlines(bool)));
    chBox_vector_plot = new QCheckBox(tr("Vector plot"));
    chBox_vector_plot->setChecked(false);

    colorMapToolBar->addWidget(TauUButton);
    //colorMapToolBar->addWidget(TauWButton);
    colorMapToolBar->addWidget(TauHButton);
    colorMapToolBar->addWidget(PressureButton);
    colorMapToolBar->addWidget(DeltaButton);
    colorMapToolBar->addSeparator();
    colorMapToolBar->addWidget(chBox_flowlines);
    colorMapToolBar->addWidget(chBox_vector_plot);

    QVBoxLayout *layout = new QVBoxLayout;
    layout->addWidget(colorMapToolBar);
    layout->addWidget(colorMap);
    layout->setMargin(1);
    setLayout(layout);
 }
GLTab::GLTab(QWidget *parent)
     : QWidget(parent) {
    glWidget = new GLWidget(this);

    QToolBar *glToolBar = new QToolBar(tr("3D"), this);
    glToolBar->setIconSize(QSize(15,15));

    QFont font = QFont("Times",12,QFont::Bold);

    QPushButton *button3dU = new QPushButton("&U");
    button3dU->setFixedSize(45,25); button3dU->setFont(font);
    connect(button3dU, SIGNAL(clicked()),glWidget, SLOT(show3dU()));

    QPushButton *button3dW = new QPushButton("&W");
    button3dW->setFixedSize(45,25); button3dW->setFont(font);
    connect(button3dW, SIGNAL(clicked()),glWidget, SLOT(show3dW()));

    QPushButton *button3dH = new QPushButton("&H");
    button3dH->setFixedSize(45,25); button3dH->setFont(font);
    connect(button3dH, SIGNAL(clicked()),glWidget, SLOT(show3dH()));

    QPushButton *button3dPressure = new QPushButton("&Pressure");
    button3dPressure->setFixedSize(65,25); button3dPressure->setFont(font);
    connect(button3dPressure, SIGNAL(clicked()),glWidget, SLOT(show3dPressure()));

    QPushButton *button3dDelta = new QPushButton("&Delta");
    button3dDelta->setFixedSize(45,25); button3dDelta->setFont(font);
    connect(button3dDelta, SIGNAL(clicked()),glWidget, SLOT(show3dDelta()));

    QPushButton *button3dBody = new QPushButton("&Body");
    button3dBody->setFixedSize(45,25); button3dBody->setFont(font);
    connect(button3dBody, SIGNAL(clicked()),glWidget, SLOT(show3dBody()));

    QCheckBox *checkBoxAxis = new QCheckBox(tr("Axis"));
    checkBoxAxis->setChecked(true);
    connect(checkBoxAxis, SIGNAL(toggled(bool)), glWidget, SLOT(showAxis(bool)));
    //connect(checkcheckBoxAxis, SIGNAL(released()), colorMap, SLOT(update()));

    glToolBar->addWidget(button3dU);
    glToolBar->addWidget(button3dW);
    glToolBar->addWidget(button3dH);
    glToolBar->addWidget(button3dPressure);
    glToolBar->addWidget(button3dDelta);
    glToolBar->addWidget(button3dBody);
    glToolBar->addWidget(checkBoxAxis);

    QVBoxLayout *layout = new QVBoxLayout;
    layout->addWidget(glToolBar);
    layout->addWidget(glWidget);
    layout->setMargin(0);
    setLayout(layout);
 }
WingView::WingView(QWidget *parent)
     : QWidget(parent) {
    //QPainter painter(this);
    //setPalette(QPalette(QColor(250, 250, 210)));
    setPalette(QPalette(QColor(255, 255, 255)));
//	setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    setAutoFillBackground(true);
}
void WingView::paintEvent(QPaintEvent * /* event */) {
   QPainter painter(this);
   painter.setRenderHint(QPainter::Antialiasing, true);
   int height, width;
   painter.save();
   //arrow for showing direction of U
   painter.translate(rect().width()/2,0);
   //painter.setBrush(QColor(167,192,220,220));
   painter.setBrush(Qt::white);
   //painter.setPen(QPen(QColor(70,81,92),0.2));
   painter.setPen(QPen(QColor(0,0,0),0.2));
   bool draw_arrow = true;
   if(draw_arrow)
   {
      QPainterPath arrow;
      arrow.moveTo(10, 0);
      arrow.lineTo(10, 40);
      arrow.lineTo(20, 40);
      arrow.lineTo(0, 55);
      arrow.lineTo(-20, 40);
      arrow.lineTo(-10, 40);
      arrow.lineTo(-10, 0);
      arrow.closeSubpath();
      painter.drawPath(arrow);
      QFont font=QFont("Times",14,QFont::Bold);
      painter.setFont(font);
      painter.setPen(Qt::black);
      painter.drawText(QRect(-5,20,10,40), Qt::AlignCenter, tr("U"));
   }
   //wing
   painter.setPen(Qt::black);
   if(Theta0d >= 90) painter.translate(0, rect().height()/2);
   else painter.translate(0, rect().height()/5);
   painter.rotate(Betad);
   height = rect().height();
   width = rect().width();
   width = int(height*0.49);
   if(!cartesian)  {
      painter.drawPie(-width, -width, width*2, width*2,
               -(int)(Theta0d+90)*16, (int)(2*Theta0d)*16);
      painter.drawLine(0, 0, 0, width);
      if (draw_grid) {
          painter.setPen(QPen(Qt::black, 0.3, Qt::SolidLine, Qt::FlatCap, Qt::MiterJoin));
          for(int k = 0; k < kmax; k++)
              painter.drawLine(0, 0, int(width*sin(Z[k]*Theta0)), int(width*cos(Z[k]*Theta0)));
          for(int i = 0; i < imax; i++)
              painter.drawArc(-width*i/(imax-1), -width*i/(imax-1),
                              2*width*i/(imax-1), 2*width*i/(imax-1),
                              -(int)(Theta0d+90)*16, (int)(2*Theta0d)*16);
      }
   }
   else if(cartesian)
   {
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
      if( L == 2.0) {
          painter.setPen(QPen(Qt::black, 1, Qt::DashLine, Qt::FlatCap, Qt::MiterJoin));
          painter.drawLine(-x,y,-x,2*y);
          painter.drawLine( x,y, x,2*y);
      }
      if (draw_grid) {
          painter.setPen(QPen(Qt::black, 0.3, Qt::SolidLine, Qt::FlatCap, Qt::MiterJoin));
          for(int k = 0; k < kmax; k++)
              for(int i = 0; i < imax-1; i++)
                  painter.drawLine(Ze[i]*y*Z[k], y*i*dx, Ze[i+1]*y*Z[k], y*(i+1)*dx);
          for(int i = 0; i < imax; i++)
              painter.drawLine(-Ze[i]*y ,y*X[i], Ze[i]*y, y*X[i]);
      }
   }
//   painter.drawArc(-30, -30, 60, 60, -(int)(Theta0d+90)*16, (int)(1*Theta0d)*16);
//   painter.rotate(-Betad);
//   if ((abs((int)Betad) > 0) && ((Theta0d < 80) ))
//   {
//	 painter.setPen(QPen(Qt::DotLine));
//	  painter.drawArc(-40, -40, 80, 80, -(int)(Betad+90)*16, (int)(1*Betad)*16);
//	  painter.drawLine(0, 0, 0, 80);
//   }
//   painter.setPen(Qt::black);
//   painter.setFont(QFont("Times",12/*,QFont::Bold*/));
//   painter.drawText(QRect(int(-21*sin(Theta0*0.5+Beta)-10), int(21*cos(Theta0*0.5+Beta))-10, 20, 20), Qt::AlignCenter, tr("%1").arg(QChar(0x03b8)));
//   if((Betad > 15) && (Theta0d < 80))
//	  painter.drawText(QRect(int(-31*sin(Beta*0.5)-10), int(31*cos(Beta*0.5)-10), 20, 20), Qt::AlignCenter, tr("%1").arg(QChar(0x03b2)));
   painter.restore();
   painter.setFont(QFont("Times",12 /*,QFont::Bold*/));
   painter.drawText(rect().width()-60, 0, 60, 25, Qt::AlignCenter, tr("%1 = %2%3").arg(QChar(0x03b8)).arg(QString::number(Theta0d)).arg(QChar(0x00B0)));
   painter.drawText(rect().width()-60, 25, 60, 25, Qt::AlignCenter, tr("%1 = %2%3").arg(QChar(0x03b2)).arg(QString::number(Betad)).arg(QChar(0x00B0)));
}
