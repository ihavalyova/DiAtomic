# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'diatomWVoSsf.ui'
##
## Created by: Qt User Interface Compiler version 5.15.6
##
## WARNING! All changes made in this file will be lost when recompiling UI file!
################################################################################

from PyQt5.QtCore import *  # type: ignore
from PyQt5.QtGui import *  # type: ignore
from PyQt5.QtWidgets import *  # type: ignore

from pyqtgraph import PlotWidget


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        if not MainWindow.objectName():
            MainWindow.setObjectName(u"MainWindow")
        MainWindow.resize(1068, 848)
        self.actionNew = QAction(MainWindow)
        self.actionNew.setObjectName(u"actionNew")
        self.centralwidget = QWidget(MainWindow)
        self.centralwidget.setObjectName(u"centralwidget")
        self.groupBox = QGroupBox(self.centralwidget)
        self.groupBox.setObjectName(u"groupBox")
        self.groupBox.setGeometry(QRect(10, 0, 531, 521))
        self.gridLayout = QGridLayout(self.groupBox)
        self.gridLayout.setObjectName(u"gridLayout")
        self.tabWidget_7 = QTabWidget(self.groupBox)
        self.tabWidget_7.setObjectName(u"tabWidget_7")
        sizePolicy = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.tabWidget_7.sizePolicy().hasHeightForWidth())
        self.tabWidget_7.setSizePolicy(sizePolicy)
        self.tab_16 = QWidget()
        self.tab_16.setObjectName(u"tab_16")
        self.widget = QWidget(self.tab_16)
        self.widget.setObjectName(u"widget")
        self.widget.setGeometry(QRect(10, 9, 491, 411))
        self.gridLayout_10 = QGridLayout(self.widget)
        self.gridLayout_10.setObjectName(u"gridLayout_10")
        self.gridLayout_10.setContentsMargins(0, 0, 0, 0)
        self.paramsFileLabel = QLabel(self.widget)
        self.paramsFileLabel.setObjectName(u"paramsFileLabel")

        self.gridLayout_10.addWidget(self.paramsFileLabel, 0, 0, 1, 1)

        self.paramsFileLineEdit = QLineEdit(self.widget)
        self.paramsFileLineEdit.setObjectName(u"paramsFileLineEdit")

        self.gridLayout_10.addWidget(self.paramsFileLineEdit, 1, 0, 1, 1)

        self.paramsFileToolButton = QToolButton(self.widget)
        self.paramsFileToolButton.setObjectName(u"paramsFileToolButton")

        self.gridLayout_10.addWidget(self.paramsFileToolButton, 1, 1, 1, 1)

        self.previewParamsPushButton = QPushButton(self.widget)
        self.previewParamsPushButton.setObjectName(u"previewParamsPushButton")

        self.gridLayout_10.addWidget(self.previewParamsPushButton, 1, 2, 1, 1)

        self.obsEnergyFileLabel = QLabel(self.widget)
        self.obsEnergyFileLabel.setObjectName(u"obsEnergyFileLabel")

        self.gridLayout_10.addWidget(self.obsEnergyFileLabel, 2, 0, 1, 1)

        self.obsEnergyFileLineEdit = QLineEdit(self.widget)
        self.obsEnergyFileLineEdit.setObjectName(u"obsEnergyFileLineEdit")

        self.gridLayout_10.addWidget(self.obsEnergyFileLineEdit, 3, 0, 1, 1)

        self.obsEnergyFileToolButton = QToolButton(self.widget)
        self.obsEnergyFileToolButton.setObjectName(u"obsEnergyFileToolButton")

        self.gridLayout_10.addWidget(self.obsEnergyFileToolButton, 3, 1, 1, 1)

        self.previewObsEnrPushButton = QPushButton(self.widget)
        self.previewObsEnrPushButton.setObjectName(u"previewObsEnrPushButton")

        self.gridLayout_10.addWidget(self.previewObsEnrPushButton, 3, 2, 1, 1)

        self.obsWavenFileLabel = QLabel(self.widget)
        self.obsWavenFileLabel.setObjectName(u"obsWavenFileLabel")

        self.gridLayout_10.addWidget(self.obsWavenFileLabel, 4, 0, 1, 1)

        self.obsWavenFileLineEdit = QLineEdit(self.widget)
        self.obsWavenFileLineEdit.setObjectName(u"obsWavenFileLineEdit")

        self.gridLayout_10.addWidget(self.obsWavenFileLineEdit, 5, 0, 1, 1)

        self.obsWavenFileToolButton = QToolButton(self.widget)
        self.obsWavenFileToolButton.setObjectName(u"obsWavenFileToolButton")

        self.gridLayout_10.addWidget(self.obsWavenFileToolButton, 5, 1, 1, 1)

        self.previewWavenPushButton = QPushButton(self.widget)
        self.previewWavenPushButton.setObjectName(u"previewWavenPushButton")

        self.gridLayout_10.addWidget(self.previewWavenPushButton, 5, 2, 1, 1)

        self.obsIntFileLabel = QLabel(self.widget)
        self.obsIntFileLabel.setObjectName(u"obsIntFileLabel")

        self.gridLayout_10.addWidget(self.obsIntFileLabel, 6, 0, 1, 1)

        self.obsIntFileLineEdit = QLineEdit(self.widget)
        self.obsIntFileLineEdit.setObjectName(u"obsIntFileLineEdit")

        self.gridLayout_10.addWidget(self.obsIntFileLineEdit, 7, 0, 1, 1)

        self.obsIntFileToolButton = QToolButton(self.widget)
        self.obsIntFileToolButton.setObjectName(u"obsIntFileToolButton")

        self.gridLayout_10.addWidget(self.obsIntFileToolButton, 7, 1, 1, 1)

        self.previewIntPushButton = QPushButton(self.widget)
        self.previewIntPushButton.setObjectName(u"previewIntPushButton")

        self.gridLayout_10.addWidget(self.previewIntPushButton, 7, 2, 1, 1)

        self.tabWidget_7.addTab(self.tab_16, "")
        self.tab_15 = QWidget()
        self.tab_15.setObjectName(u"tab_15")
        self.tabWidget = QTabWidget(self.tab_15)
        self.tabWidget.setObjectName(u"tabWidget")
        self.tabWidget.setEnabled(True)
        self.tabWidget.setGeometry(QRect(-4, -1, 511, 451))
        sizePolicy1 = QSizePolicy(QSizePolicy.Minimum, QSizePolicy.Minimum)
        sizePolicy1.setHorizontalStretch(0)
        sizePolicy1.setVerticalStretch(0)
        sizePolicy1.setHeightForWidth(self.tabWidget.sizePolicy().hasHeightForWidth())
        self.tabWidget.setSizePolicy(sizePolicy1)
        self.tab1 = QWidget()
        self.tab1.setObjectName(u"tab1")
        self.layoutWidget = QWidget(self.tab1)
        self.layoutWidget.setObjectName(u"layoutWidget")
        self.layoutWidget.setGeometry(QRect(261, 0, 241, 58))
        self.horizontalLayout_11 = QHBoxLayout(self.layoutWidget)
        self.horizontalLayout_11.setObjectName(u"horizontalLayout_11")
        self.horizontalLayout_11.setContentsMargins(0, 0, 0, 0)
        self.label_10 = QLabel(self.layoutWidget)
        self.label_10.setObjectName(u"label_10")

        self.horizontalLayout_11.addWidget(self.label_10)

        self.verticalLayout_2 = QVBoxLayout()
        self.verticalLayout_2.setObjectName(u"verticalLayout_2")
        self.eqidGridRadioButton = QRadioButton(self.layoutWidget)
        self.eqidGridRadioButton.setObjectName(u"eqidGridRadioButton")
        self.eqidGridRadioButton.setChecked(True)

        self.verticalLayout_2.addWidget(self.eqidGridRadioButton)

        self.noneqidGridRadioButton = QRadioButton(self.layoutWidget)
        self.noneqidGridRadioButton.setObjectName(u"noneqidGridRadioButton")

        self.verticalLayout_2.addWidget(self.noneqidGridRadioButton)


        self.horizontalLayout_11.addLayout(self.verticalLayout_2)

        self.layoutWidget1 = QWidget(self.tab1)
        self.layoutWidget1.setObjectName(u"layoutWidget1")
        self.layoutWidget1.setGeometry(QRect(260, 160, 241, 52))
        self.verticalLayout = QVBoxLayout(self.layoutWidget1)
        self.verticalLayout.setObjectName(u"verticalLayout")
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.solverLabel = QLabel(self.layoutWidget1)
        self.solverLabel.setObjectName(u"solverLabel")

        self.verticalLayout.addWidget(self.solverLabel)

        self.horizontalLayout_10 = QHBoxLayout()
        self.horizontalLayout_10.setObjectName(u"horizontalLayout_10")
        self.sincSolverRadioButton = QRadioButton(self.layoutWidget1)
        self.sincSolverRadioButton.setObjectName(u"sincSolverRadioButton")
        self.sincSolverRadioButton.setChecked(True)

        self.horizontalLayout_10.addWidget(self.sincSolverRadioButton)

        self.fourierSolverRadioButton = QRadioButton(self.layoutWidget1)
        self.fourierSolverRadioButton.setObjectName(u"fourierSolverRadioButton")

        self.horizontalLayout_10.addWidget(self.fourierSolverRadioButton)

        self.fd5SolverRadioButton = QRadioButton(self.layoutWidget1)
        self.fd5SolverRadioButton.setObjectName(u"fd5SolverRadioButton")

        self.horizontalLayout_10.addWidget(self.fd5SolverRadioButton)


        self.verticalLayout.addLayout(self.horizontalLayout_10)

        self.layoutWidget2 = QWidget(self.tab1)
        self.layoutWidget2.setObjectName(u"layoutWidget2")
        self.layoutWidget2.setGeometry(QRect(262, 62, 241, 92))
        self.gridLayout_3 = QGridLayout(self.layoutWidget2)
        self.gridLayout_3.setObjectName(u"gridLayout_3")
        self.gridLayout_3.setContentsMargins(0, 0, 0, 0)
        self.nGridPointsLabel = QLabel(self.layoutWidget2)
        self.nGridPointsLabel.setObjectName(u"nGridPointsLabel")

        self.gridLayout_3.addWidget(self.nGridPointsLabel, 0, 0, 1, 2)

        self.nGridPointsLineEdit = QLineEdit(self.layoutWidget2)
        self.nGridPointsLineEdit.setObjectName(u"nGridPointsLineEdit")

        self.gridLayout_3.addWidget(self.nGridPointsLineEdit, 0, 2, 1, 1)

        self.rminLabel = QLabel(self.layoutWidget2)
        self.rminLabel.setObjectName(u"rminLabel")

        self.gridLayout_3.addWidget(self.rminLabel, 1, 0, 1, 1)

        self.rminLineEdit = QLineEdit(self.layoutWidget2)
        self.rminLineEdit.setObjectName(u"rminLineEdit")

        self.gridLayout_3.addWidget(self.rminLineEdit, 1, 1, 1, 2)

        self.rmaxLabel = QLabel(self.layoutWidget2)
        self.rmaxLabel.setObjectName(u"rmaxLabel")

        self.gridLayout_3.addWidget(self.rmaxLabel, 2, 0, 1, 1)

        self.rmaxLineEdit = QLineEdit(self.layoutWidget2)
        self.rmaxLineEdit.setObjectName(u"rmaxLineEdit")

        self.gridLayout_3.addWidget(self.rmaxLineEdit, 2, 1, 1, 2)

        self.layoutWidget3 = QWidget(self.tab1)
        self.layoutWidget3.setObjectName(u"layoutWidget3")
        self.layoutWidget3.setGeometry(QRect(11, 0, 241, 211))
        self.verticalLayout_4 = QVBoxLayout(self.layoutWidget3)
        self.verticalLayout_4.setObjectName(u"verticalLayout_4")
        self.verticalLayout_4.setContentsMargins(0, 0, 0, 0)
        self.moleculeLabel = QLabel(self.layoutWidget3)
        self.moleculeLabel.setObjectName(u"moleculeLabel")

        self.verticalLayout_4.addWidget(self.moleculeLabel)

        self.moleculeLineEdit = QLineEdit(self.layoutWidget3)
        self.moleculeLineEdit.setObjectName(u"moleculeLineEdit")

        self.verticalLayout_4.addWidget(self.moleculeLineEdit)

        self.isotopesListWidget = QListWidget(self.layoutWidget3)
        self.isotopesListWidget.setObjectName(u"isotopesListWidget")

        self.verticalLayout_4.addWidget(self.isotopesListWidget)

        self.layoutWidget4 = QWidget(self.tab1)
        self.layoutWidget4.setObjectName(u"layoutWidget4")
        self.layoutWidget4.setGeometry(QRect(10, 220, 491, 191))
        self.verticalLayout_6 = QVBoxLayout(self.layoutWidget4)
        self.verticalLayout_6.setObjectName(u"verticalLayout_6")
        self.verticalLayout_6.setContentsMargins(0, 0, 0, 0)
        self.isotopesTextEdit = QPlainTextEdit(self.layoutWidget4)
        self.isotopesTextEdit.setObjectName(u"isotopesTextEdit")

        self.verticalLayout_6.addWidget(self.isotopesTextEdit)

        self.horizontalLayout_12 = QHBoxLayout()
        self.horizontalLayout_12.setObjectName(u"horizontalLayout_12")
        self.refJLabel = QLabel(self.layoutWidget4)
        self.refJLabel.setObjectName(u"refJLabel")

        self.horizontalLayout_12.addWidget(self.refJLabel)

        self.refJLineEdit = QLineEdit(self.layoutWidget4)
        self.refJLineEdit.setObjectName(u"refJLineEdit")

        self.horizontalLayout_12.addWidget(self.refJLineEdit)

        self.refEnergyLabel = QLabel(self.layoutWidget4)
        self.refEnergyLabel.setObjectName(u"refEnergyLabel")

        self.horizontalLayout_12.addWidget(self.refEnergyLabel)

        self.refEnergyLineEdit = QLineEdit(self.layoutWidget4)
        self.refEnergyLineEdit.setObjectName(u"refEnergyLineEdit")

        self.horizontalLayout_12.addWidget(self.refEnergyLineEdit)


        self.verticalLayout_6.addLayout(self.horizontalLayout_12)

        self.tabWidget.addTab(self.tab1, "")
        self.tab4 = QWidget()
        self.tab4.setObjectName(u"tab4")
        self.scrollArea = QScrollArea(self.tab4)
        self.scrollArea.setObjectName(u"scrollArea")
        self.scrollArea.setGeometry(QRect(10, 80, 491, 331))
        sizePolicy2 = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        sizePolicy2.setHorizontalStretch(0)
        sizePolicy2.setVerticalStretch(0)
        sizePolicy2.setHeightForWidth(self.scrollArea.sizePolicy().hasHeightForWidth())
        self.scrollArea.setSizePolicy(sizePolicy2)
        self.scrollArea.setMinimumSize(QSize(0, 0))
        self.scrollArea.setWidgetResizable(True)
        self.scrollAreaWidgetContents = QWidget()
        self.scrollAreaWidgetContents.setObjectName(u"scrollAreaWidgetContents")
        self.scrollAreaWidgetContents.setGeometry(QRect(0, 0, 489, 146))
        sizePolicy3 = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        sizePolicy3.setHorizontalStretch(0)
        sizePolicy3.setVerticalStretch(0)
        sizePolicy3.setHeightForWidth(self.scrollAreaWidgetContents.sizePolicy().hasHeightForWidth())
        self.scrollAreaWidgetContents.setSizePolicy(sizePolicy3)
        self.gridLayout_7 = QGridLayout(self.scrollAreaWidgetContents)
        self.gridLayout_7.setObjectName(u"gridLayout_7")
        self.gridLayout_4 = QGridLayout()
        self.gridLayout_4.setObjectName(u"gridLayout_4")
        self.gridLayout_4.setSizeConstraint(QLayout.SetDefaultConstraint)
        self.gridLayout_4.setVerticalSpacing(6)
        self.gridLayout_4.setContentsMargins(-1, -1, -1, 0)
        self.doubleSpinBox_4 = QDoubleSpinBox(self.scrollAreaWidgetContents)
        self.doubleSpinBox_4.setObjectName(u"doubleSpinBox_4")
        self.doubleSpinBox_4.setDecimals(1)

        self.gridLayout_4.addWidget(self.doubleSpinBox_4, 2, 1, 1, 1)

        self.label_8 = QLabel(self.scrollAreaWidgetContents)
        self.label_8.setObjectName(u"label_8")

        self.gridLayout_4.addWidget(self.label_8, 3, 3, 1, 2)

        self.label_4 = QLabel(self.scrollAreaWidgetContents)
        self.label_4.setObjectName(u"label_4")

        self.gridLayout_4.addWidget(self.label_4, 0, 3, 1, 2)

        self.label = QLabel(self.scrollAreaWidgetContents)
        self.label.setObjectName(u"label")
        self.label.setMinimumSize(QSize(0, 0))

        self.gridLayout_4.addWidget(self.label, 0, 0, 1, 1)

        self.checkBox = QCheckBox(self.scrollAreaWidgetContents)
        self.checkBox.setObjectName(u"checkBox")
        self.checkBox.setChecked(True)

        self.gridLayout_4.addWidget(self.checkBox, 0, 5, 1, 1)

        self.spinBox_3 = QSpinBox(self.scrollAreaWidgetContents)
        self.spinBox_3.setObjectName(u"spinBox_3")

        self.gridLayout_4.addWidget(self.spinBox_3, 1, 5, 1, 1)

        self.doubleSpinBox_2 = QDoubleSpinBox(self.scrollAreaWidgetContents)
        self.doubleSpinBox_2.setObjectName(u"doubleSpinBox_2")
        sizePolicy4 = QSizePolicy(QSizePolicy.Minimum, QSizePolicy.Fixed)
        sizePolicy4.setHorizontalStretch(0)
        sizePolicy4.setVerticalStretch(0)
        sizePolicy4.setHeightForWidth(self.doubleSpinBox_2.sizePolicy().hasHeightForWidth())
        self.doubleSpinBox_2.setSizePolicy(sizePolicy4)
        self.doubleSpinBox_2.setDecimals(1)

        self.gridLayout_4.addWidget(self.doubleSpinBox_2, 3, 5, 1, 1)

        self.doubleSpinBox_3 = QDoubleSpinBox(self.scrollAreaWidgetContents)
        self.doubleSpinBox_3.setObjectName(u"doubleSpinBox_3")
        self.doubleSpinBox_3.setDecimals(1)

        self.gridLayout_4.addWidget(self.doubleSpinBox_3, 1, 1, 1, 1)

        self.spinBox_4 = QSpinBox(self.scrollAreaWidgetContents)
        self.spinBox_4.setObjectName(u"spinBox_4")
        self.spinBox_4.setEnabled(True)
        self.spinBox_4.setReadOnly(True)
        self.spinBox_4.setMinimum(1)

        self.gridLayout_4.addWidget(self.spinBox_4, 0, 1, 1, 1)

        self.label_3 = QLabel(self.scrollAreaWidgetContents)
        self.label_3.setObjectName(u"label_3")

        self.gridLayout_4.addWidget(self.label_3, 2, 0, 1, 1)

        self.label_2 = QLabel(self.scrollAreaWidgetContents)
        self.label_2.setObjectName(u"label_2")

        self.gridLayout_4.addWidget(self.label_2, 1, 0, 1, 1)

        self.doubleSpinBox = QDoubleSpinBox(self.scrollAreaWidgetContents)
        self.doubleSpinBox.setObjectName(u"doubleSpinBox")
        self.doubleSpinBox.setDecimals(1)

        self.gridLayout_4.addWidget(self.doubleSpinBox, 2, 5, 1, 1)

        self.label_6 = QLabel(self.scrollAreaWidgetContents)
        self.label_6.setObjectName(u"label_6")

        self.gridLayout_4.addWidget(self.label_6, 1, 3, 1, 2)

        self.lineEdit_2 = QLineEdit(self.scrollAreaWidgetContents)
        self.lineEdit_2.setObjectName(u"lineEdit_2")

        self.gridLayout_4.addWidget(self.lineEdit_2, 3, 1, 1, 2)

        self.label_9 = QLabel(self.scrollAreaWidgetContents)
        self.label_9.setObjectName(u"label_9")

        self.gridLayout_4.addWidget(self.label_9, 3, 0, 1, 1)

        self.label_7 = QLabel(self.scrollAreaWidgetContents)
        self.label_7.setObjectName(u"label_7")

        self.gridLayout_4.addWidget(self.label_7, 2, 3, 1, 2)

        self.checkBox_2 = QCheckBox(self.scrollAreaWidgetContents)
        self.checkBox_2.setObjectName(u"checkBox_2")

        self.gridLayout_4.addWidget(self.checkBox_2, 0, 6, 1, 1)


        self.gridLayout_7.addLayout(self.gridLayout_4, 0, 0, 1, 1)

        self.scrollArea.setWidget(self.scrollAreaWidgetContents)
        self.layoutWidget5 = QWidget(self.tab4)
        self.layoutWidget5.setObjectName(u"layoutWidget5")
        self.layoutWidget5.setGeometry(QRect(11, 2, 491, 65))
        self.verticalLayout_3 = QVBoxLayout(self.layoutWidget5)
        self.verticalLayout_3.setObjectName(u"verticalLayout_3")
        self.verticalLayout_3.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout = QHBoxLayout()
        self.horizontalLayout.setObjectName(u"horizontalLayout")
        self.label_11 = QLabel(self.layoutWidget5)
        self.label_11.setObjectName(u"label_11")

        self.horizontalLayout.addWidget(self.label_11)

        self.spinBox_5 = QSpinBox(self.layoutWidget5)
        self.spinBox_5.setObjectName(u"spinBox_5")
        self.spinBox_5.setMinimum(1)

        self.horizontalLayout.addWidget(self.spinBox_5)


        self.verticalLayout_3.addLayout(self.horizontalLayout)

        self.horizontalLayout_2 = QHBoxLayout()
        self.horizontalLayout_2.setObjectName(u"horizontalLayout_2")
        self.addStatePushButton = QPushButton(self.layoutWidget5)
        self.addStatePushButton.setObjectName(u"addStatePushButton")

        self.horizontalLayout_2.addWidget(self.addStatePushButton)

        self.removeStatePushButton = QPushButton(self.layoutWidget5)
        self.removeStatePushButton.setObjectName(u"removeStatePushButton")

        self.horizontalLayout_2.addWidget(self.removeStatePushButton)


        self.verticalLayout_3.addLayout(self.horizontalLayout_2)

        self.tabWidget.addTab(self.tab4, "")
        self.tab5 = QWidget()
        self.tab5.setObjectName(u"tab5")
        self.scrollArea_2 = QScrollArea(self.tab5)
        self.scrollArea_2.setObjectName(u"scrollArea_2")
        self.scrollArea_2.setGeometry(QRect(9, 49, 491, 361))
        self.scrollArea_2.setWidgetResizable(True)
        self.scrollAreaWidgetContents_2 = QWidget()
        self.scrollAreaWidgetContents_2.setObjectName(u"scrollAreaWidgetContents_2")
        self.scrollAreaWidgetContents_2.setGeometry(QRect(0, 0, 489, 207))
        sizePolicy3.setHeightForWidth(self.scrollAreaWidgetContents_2.sizePolicy().hasHeightForWidth())
        self.scrollAreaWidgetContents_2.setSizePolicy(sizePolicy3)
        self.gridLayout_8 = QGridLayout(self.scrollAreaWidgetContents_2)
        self.gridLayout_8.setObjectName(u"gridLayout_8")
        self.gridLayout_6 = QGridLayout()
        self.gridLayout_6.setObjectName(u"gridLayout_6")
        self.pairStates2SpinBox1 = QSpinBox(self.scrollAreaWidgetContents_2)
        self.pairStates2SpinBox1.setObjectName(u"pairStates2SpinBox1")
        self.pairStates2SpinBox1.setMinimum(1)

        self.gridLayout_6.addWidget(self.pairStates2SpinBox1, 1, 2, 1, 1)

        self.labelParamsLabel1 = QLabel(self.scrollAreaWidgetContents_2)
        self.labelParamsLabel1.setObjectName(u"labelParamsLabel1")

        self.gridLayout_6.addWidget(self.labelParamsLabel1, 3, 0, 1, 1)

        self.pairStatesLabel1 = QLabel(self.scrollAreaWidgetContents_2)
        self.pairStatesLabel1.setObjectName(u"pairStatesLabel1")

        self.gridLayout_6.addWidget(self.pairStatesLabel1, 1, 0, 1, 1)

        self.modelTypeComboBox1 = QComboBox(self.scrollAreaWidgetContents_2)
        self.modelTypeComboBox1.setObjectName(u"modelTypeComboBox1")

        self.gridLayout_6.addWidget(self.modelTypeComboBox1, 2, 1, 1, 2)

        self.labelRotCorrLineEdit1 = QLineEdit(self.scrollAreaWidgetContents_2)
        self.labelRotCorrLineEdit1.setObjectName(u"labelRotCorrLineEdit1")

        self.gridLayout_6.addWidget(self.labelRotCorrLineEdit1, 4, 1, 1, 2)

        self.modelTypeLabel1 = QLabel(self.scrollAreaWidgetContents_2)
        self.modelTypeLabel1.setObjectName(u"modelTypeLabel1")

        self.gridLayout_6.addWidget(self.modelTypeLabel1, 2, 0, 1, 1)

        self.labelRadParamsLineEdit1 = QLineEdit(self.scrollAreaWidgetContents_2)
        self.labelRadParamsLineEdit1.setObjectName(u"labelRadParamsLineEdit1")

        self.gridLayout_6.addWidget(self.labelRadParamsLineEdit1, 3, 1, 1, 2)

        self.pairStates1SpinBox1 = QSpinBox(self.scrollAreaWidgetContents_2)
        self.pairStates1SpinBox1.setObjectName(u"pairStates1SpinBox1")
        self.pairStates1SpinBox1.setMinimum(1)

        self.gridLayout_6.addWidget(self.pairStates1SpinBox1, 1, 1, 1, 1)

        self.operatorTypeComboBox1 = QComboBox(self.scrollAreaWidgetContents_2)
        self.operatorTypeComboBox1.setObjectName(u"operatorTypeComboBox1")

        self.gridLayout_6.addWidget(self.operatorTypeComboBox1, 0, 1, 1, 2)

        self.rotCorrLabel1 = QLabel(self.scrollAreaWidgetContents_2)
        self.rotCorrLabel1.setObjectName(u"rotCorrLabel1")

        self.gridLayout_6.addWidget(self.rotCorrLabel1, 4, 0, 1, 1)

        self.operatorTypeLabel1 = QLabel(self.scrollAreaWidgetContents_2)
        self.operatorTypeLabel1.setObjectName(u"operatorTypeLabel1")

        self.gridLayout_6.addWidget(self.operatorTypeLabel1, 0, 0, 1, 1)

        self.showRadParamsPushButton1 = QPushButton(self.scrollAreaWidgetContents_2)
        self.showRadParamsPushButton1.setObjectName(u"showRadParamsPushButton1")

        self.gridLayout_6.addWidget(self.showRadParamsPushButton1, 5, 0, 1, 1)

        self.plotRadParamsPushButton1 = QPushButton(self.scrollAreaWidgetContents_2)
        self.plotRadParamsPushButton1.setObjectName(u"plotRadParamsPushButton1")

        self.gridLayout_6.addWidget(self.plotRadParamsPushButton1, 5, 1, 1, 2)


        self.gridLayout_8.addLayout(self.gridLayout_6, 0, 0, 1, 1)

        self.scrollArea_2.setWidget(self.scrollAreaWidgetContents_2)
        self.layoutWidget6 = QWidget(self.tab5)
        self.layoutWidget6.setObjectName(u"layoutWidget6")
        self.layoutWidget6.setGeometry(QRect(10, 10, 481, 28))
        self.horizontalLayout_3 = QHBoxLayout(self.layoutWidget6)
        self.horizontalLayout_3.setObjectName(u"horizontalLayout_3")
        self.horizontalLayout_3.setContentsMargins(0, 0, 0, 0)
        self.addOperatorPushButton = QPushButton(self.layoutWidget6)
        self.addOperatorPushButton.setObjectName(u"addOperatorPushButton")

        self.horizontalLayout_3.addWidget(self.addOperatorPushButton)

        self.removeOperatorPushButton = QPushButton(self.layoutWidget6)
        self.removeOperatorPushButton.setObjectName(u"removeOperatorPushButton")

        self.horizontalLayout_3.addWidget(self.removeOperatorPushButton)

        self.tabWidget.addTab(self.tab5, "")
        self.tab6 = QWidget()
        self.tab6.setObjectName(u"tab6")
        self.scrollArea_3 = QScrollArea(self.tab6)
        self.scrollArea_3.setObjectName(u"scrollArea_3")
        self.scrollArea_3.setGeometry(QRect(10, 0, 491, 371))
        self.scrollArea_3.setWidgetResizable(True)
        self.scrollAreaWidgetContents_3 = QWidget()
        self.scrollAreaWidgetContents_3.setObjectName(u"scrollAreaWidgetContents_3")
        self.scrollAreaWidgetContents_3.setGeometry(QRect(0, 0, 491, 236))
        sizePolicy3.setHeightForWidth(self.scrollAreaWidgetContents_3.sizePolicy().hasHeightForWidth())
        self.scrollAreaWidgetContents_3.setSizePolicy(sizePolicy3)
        self.gridLayout_9 = QGridLayout(self.scrollAreaWidgetContents_3)
        self.gridLayout_9.setObjectName(u"gridLayout_9")
        self.gridLayout_5 = QGridLayout()
        self.gridLayout_5.setObjectName(u"gridLayout_5")
        self.calcEnergyLevelsCheckBox = QCheckBox(self.scrollAreaWidgetContents_3)
        self.calcEnergyLevelsCheckBox.setObjectName(u"calcEnergyLevelsCheckBox")
        sizePolicy5 = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        sizePolicy5.setHorizontalStretch(0)
        sizePolicy5.setVerticalStretch(0)
        sizePolicy5.setHeightForWidth(self.calcEnergyLevelsCheckBox.sizePolicy().hasHeightForWidth())
        self.calcEnergyLevelsCheckBox.setSizePolicy(sizePolicy5)

        self.gridLayout_5.addWidget(self.calcEnergyLevelsCheckBox, 2, 2, 1, 2)

        self.label_13 = QLabel(self.scrollAreaWidgetContents_3)
        self.label_13.setObjectName(u"label_13")

        self.gridLayout_5.addWidget(self.label_13, 1, 0, 1, 1)

        self.energyRange2LineEdit = QLineEdit(self.scrollAreaWidgetContents_3)
        self.energyRange2LineEdit.setObjectName(u"energyRange2LineEdit")
        sizePolicy3.setHeightForWidth(self.energyRange2LineEdit.sizePolicy().hasHeightForWidth())
        self.energyRange2LineEdit.setSizePolicy(sizePolicy3)
        self.energyRange2LineEdit.setMinimumSize(QSize(145, 0))

        self.gridLayout_5.addWidget(self.energyRange2LineEdit, 0, 2, 1, 1)

        self.calcEigenenergyCheckBox = QCheckBox(self.scrollAreaWidgetContents_3)
        self.calcEigenenergyCheckBox.setObjectName(u"calcEigenenergyCheckBox")
        sizePolicy5.setHeightForWidth(self.calcEigenenergyCheckBox.sizePolicy().hasHeightForWidth())
        self.calcEigenenergyCheckBox.setSizePolicy(sizePolicy5)
        self.calcEigenenergyCheckBox.setChecked(True)

        self.gridLayout_5.addWidget(self.calcEigenenergyCheckBox, 2, 1, 1, 1)

        self.energyIndex1LineEdit = QLineEdit(self.scrollAreaWidgetContents_3)
        self.energyIndex1LineEdit.setObjectName(u"energyIndex1LineEdit")
        sizePolicy6 = QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        sizePolicy6.setHorizontalStretch(0)
        sizePolicy6.setVerticalStretch(0)
        sizePolicy6.setHeightForWidth(self.energyIndex1LineEdit.sizePolicy().hasHeightForWidth())
        self.energyIndex1LineEdit.setSizePolicy(sizePolicy6)
        self.energyIndex1LineEdit.setMinimumSize(QSize(145, 0))

        self.gridLayout_5.addWidget(self.energyIndex1LineEdit, 1, 1, 1, 1)

        self.identLevelsByEnrRadioButton = QRadioButton(self.scrollAreaWidgetContents_3)
        self.buttonGroup_2 = QButtonGroup(MainWindow)
        self.buttonGroup_2.setObjectName(u"buttonGroup_2")
        self.buttonGroup_2.addButton(self.identLevelsByEnrRadioButton)
        self.identLevelsByEnrRadioButton.setObjectName(u"identLevelsByEnrRadioButton")
        sizePolicy5.setHeightForWidth(self.identLevelsByEnrRadioButton.sizePolicy().hasHeightForWidth())
        self.identLevelsByEnrRadioButton.setSizePolicy(sizePolicy5)

        self.gridLayout_5.addWidget(self.identLevelsByEnrRadioButton, 3, 2, 1, 2)

        self.label_15 = QLabel(self.scrollAreaWidgetContents_3)
        self.label_15.setObjectName(u"label_15")

        self.gridLayout_5.addWidget(self.label_15, 3, 0, 1, 1)

        self.label_12 = QLabel(self.scrollAreaWidgetContents_3)
        self.label_12.setObjectName(u"label_12")

        self.gridLayout_5.addWidget(self.label_12, 0, 0, 1, 1)

        self.saveEigenvectrosCheckBox = QCheckBox(self.scrollAreaWidgetContents_3)
        self.saveEigenvectrosCheckBox.setObjectName(u"saveEigenvectrosCheckBox")

        self.gridLayout_5.addWidget(self.saveEigenvectrosCheckBox, 4, 0, 1, 1)

        self.plotNlevelsLineEdit = QLineEdit(self.scrollAreaWidgetContents_3)
        self.plotNlevelsLineEdit.setObjectName(u"plotNlevelsLineEdit")
        sizePolicy5.setHeightForWidth(self.plotNlevelsLineEdit.sizePolicy().hasHeightForWidth())
        self.plotNlevelsLineEdit.setSizePolicy(sizePolicy5)

        self.gridLayout_5.addWidget(self.plotNlevelsLineEdit, 4, 1, 1, 3)

        self.energyIndex2LineEdit = QLineEdit(self.scrollAreaWidgetContents_3)
        self.energyIndex2LineEdit.setObjectName(u"energyIndex2LineEdit")
        sizePolicy6.setHeightForWidth(self.energyIndex2LineEdit.sizePolicy().hasHeightForWidth())
        self.energyIndex2LineEdit.setSizePolicy(sizePolicy6)
        self.energyIndex2LineEdit.setMinimumSize(QSize(145, 0))

        self.gridLayout_5.addWidget(self.energyIndex2LineEdit, 1, 2, 1, 1)

        self.label_14 = QLabel(self.scrollAreaWidgetContents_3)
        self.label_14.setObjectName(u"label_14")

        self.gridLayout_5.addWidget(self.label_14, 2, 0, 1, 1)

        self.identLevelsByContrRadioButton = QRadioButton(self.scrollAreaWidgetContents_3)
        self.buttonGroup_2.addButton(self.identLevelsByContrRadioButton)
        self.identLevelsByContrRadioButton.setObjectName(u"identLevelsByContrRadioButton")
        sizePolicy5.setHeightForWidth(self.identLevelsByContrRadioButton.sizePolicy().hasHeightForWidth())
        self.identLevelsByContrRadioButton.setSizePolicy(sizePolicy5)
        self.identLevelsByContrRadioButton.setChecked(True)

        self.gridLayout_5.addWidget(self.identLevelsByContrRadioButton, 3, 1, 1, 1)

        self.energyRange1LineEdit = QLineEdit(self.scrollAreaWidgetContents_3)
        self.energyRange1LineEdit.setObjectName(u"energyRange1LineEdit")
        sizePolicy6.setHeightForWidth(self.energyRange1LineEdit.sizePolicy().hasHeightForWidth())
        self.energyRange1LineEdit.setSizePolicy(sizePolicy6)
        self.energyRange1LineEdit.setMinimumSize(QSize(145, 0))

        self.gridLayout_5.addWidget(self.energyRange1LineEdit, 0, 1, 1, 1)


        self.gridLayout_9.addLayout(self.gridLayout_5, 1, 0, 2, 1)

        self.horizontalLayout_5 = QHBoxLayout()
        self.horizontalLayout_5.setObjectName(u"horizontalLayout_5")
        self.label_5 = QLabel(self.scrollAreaWidgetContents_3)
        self.label_5.setObjectName(u"label_5")

        self.horizontalLayout_5.addWidget(self.label_5)

        self.horizontalLayout_4 = QHBoxLayout()
        self.horizontalLayout_4.setObjectName(u"horizontalLayout_4")
        self.radioButton = QRadioButton(self.scrollAreaWidgetContents_3)
        self.buttonGroup = QButtonGroup(MainWindow)
        self.buttonGroup.setObjectName(u"buttonGroup")
        self.buttonGroup.addButton(self.radioButton)
        self.radioButton.setObjectName(u"radioButton")
        self.radioButton.setChecked(True)

        self.horizontalLayout_4.addWidget(self.radioButton)

        self.radioButton_2 = QRadioButton(self.scrollAreaWidgetContents_3)
        self.buttonGroup.addButton(self.radioButton_2)
        self.radioButton_2.setObjectName(u"radioButton_2")

        self.horizontalLayout_4.addWidget(self.radioButton_2)


        self.horizontalLayout_5.addLayout(self.horizontalLayout_4)


        self.gridLayout_9.addLayout(self.horizontalLayout_5, 0, 0, 1, 1)

        self.previewOutputPushButton = QPushButton(self.scrollAreaWidgetContents_3)
        self.previewOutputPushButton.setObjectName(u"previewOutputPushButton")

        self.gridLayout_9.addWidget(self.previewOutputPushButton, 3, 0, 1, 1)

        self.scrollArea_3.setWidget(self.scrollAreaWidgetContents_3)
        self.solvePushButton = QPushButton(self.tab6)
        self.solvePushButton.setObjectName(u"solvePushButton")
        self.solvePushButton.setGeometry(QRect(370, 380, 131, 31))
        self.tabWidget.addTab(self.tab6, "")
        self.tabWidget_7.addTab(self.tab_15, "")
        self.tab_17 = QWidget()
        self.tab_17.setObjectName(u"tab_17")
        self.widget1 = QWidget(self.tab_17)
        self.widget1.setObjectName(u"widget1")
        self.widget1.setGeometry(QRect(10, 10, 481, 191))
        self.gridLayout_2 = QGridLayout(self.widget1)
        self.gridLayout_2.setObjectName(u"gridLayout_2")
        self.gridLayout_2.setContentsMargins(0, 0, 0, 0)
        self.label_16 = QLabel(self.widget1)
        self.label_16.setObjectName(u"label_16")

        self.gridLayout_2.addWidget(self.label_16, 0, 0, 1, 2)

        self.numberOfIterLineEdit = QLineEdit(self.widget1)
        self.numberOfIterLineEdit.setObjectName(u"numberOfIterLineEdit")

        self.gridLayout_2.addWidget(self.numberOfIterLineEdit, 3, 2, 1, 3)

        self.label_19 = QLabel(self.widget1)
        self.label_19.setObjectName(u"label_19")

        self.gridLayout_2.addWidget(self.label_19, 3, 0, 1, 2)

        self.toleranceLineEdit = QLineEdit(self.widget1)
        self.toleranceLineEdit.setObjectName(u"toleranceLineEdit")

        self.gridLayout_2.addWidget(self.toleranceLineEdit, 1, 2, 1, 3)

        self.label_18 = QLabel(self.widget1)
        self.label_18.setObjectName(u"label_18")

        self.gridLayout_2.addWidget(self.label_18, 1, 0, 1, 2)

        self.label_17 = QLabel(self.widget1)
        self.label_17.setObjectName(u"label_17")

        self.gridLayout_2.addWidget(self.label_17, 2, 0, 1, 2)

        self.numericalDerivRadioButton = QRadioButton(self.widget1)
        self.buttonGroup_3 = QButtonGroup(MainWindow)
        self.buttonGroup_3.setObjectName(u"buttonGroup_3")
        self.buttonGroup_3.addButton(self.numericalDerivRadioButton)
        self.numericalDerivRadioButton.setObjectName(u"numericalDerivRadioButton")
        self.numericalDerivRadioButton.setChecked(True)

        self.gridLayout_2.addWidget(self.numericalDerivRadioButton, 2, 2, 1, 1)

        self.analyticalDerivRadioButton = QRadioButton(self.widget1)
        self.buttonGroup_3.addButton(self.analyticalDerivRadioButton)
        self.analyticalDerivRadioButton.setObjectName(u"analyticalDerivRadioButton")

        self.gridLayout_2.addWidget(self.analyticalDerivRadioButton, 2, 3, 1, 2)

        self.svdOptMethodRadioButton = QRadioButton(self.widget1)
        self.buttonGroup_4 = QButtonGroup(MainWindow)
        self.buttonGroup_4.setObjectName(u"buttonGroup_4")
        self.buttonGroup_4.addButton(self.svdOptMethodRadioButton)
        self.svdOptMethodRadioButton.setObjectName(u"svdOptMethodRadioButton")
        self.svdOptMethodRadioButton.setChecked(True)

        self.gridLayout_2.addWidget(self.svdOptMethodRadioButton, 0, 2, 1, 1)

        self.minuitOptMethodRadioButton = QRadioButton(self.widget1)
        self.buttonGroup_4.addButton(self.minuitOptMethodRadioButton)
        self.minuitOptMethodRadioButton.setObjectName(u"minuitOptMethodRadioButton")

        self.gridLayout_2.addWidget(self.minuitOptMethodRadioButton, 0, 3, 1, 2)

        self.tabWidget_7.addTab(self.tab_17, "")
        self.tab_18 = QWidget()
        self.tab_18.setObjectName(u"tab_18")
        self.tabWidget_7.addTab(self.tab_18, "")

        self.gridLayout.addWidget(self.tabWidget_7, 0, 1, 1, 1)

        self.groupBox_3 = QGroupBox(self.centralwidget)
        self.groupBox_3.setObjectName(u"groupBox_3")
        self.groupBox_3.setGeometry(QRect(10, 519, 1051, 281))
        self.tabWidget_5 = QTabWidget(self.groupBox_3)
        self.tabWidget_5.setObjectName(u"tabWidget_5")
        self.tabWidget_5.setGeometry(QRect(0, 20, 1041, 280))
        self.tab_10 = QWidget()
        self.tab_10.setObjectName(u"tab_10")
        self.tabWidget_3 = QTabWidget(self.tab_10)
        self.tabWidget_3.setObjectName(u"tabWidget_3")
        self.tabWidget_3.setGeometry(QRect(-4, -1, 1051, 252))
        self.tab_3 = QWidget()
        self.tab_3.setObjectName(u"tab_3")
        self.tableWidget = QTableWidget(self.tab_3)
        self.tableWidget.setObjectName(u"tableWidget")
        self.tableWidget.setGeometry(QRect(0, -1, 1051, 201))
        self.tabWidget_3.addTab(self.tab_3, "")
        self.tab_4 = QWidget()
        self.tab_4.setObjectName(u"tab_4")
        self.tableWidget_2 = QTableWidget(self.tab_4)
        self.tableWidget_2.setObjectName(u"tableWidget_2")
        self.tableWidget_2.setGeometry(QRect(0, -1, 1041, 225))
        self.tabWidget_3.addTab(self.tab_4, "")
        self.tab_5 = QWidget()
        self.tab_5.setObjectName(u"tab_5")
        self.tableWidget_3 = QTableWidget(self.tab_5)
        self.tableWidget_3.setObjectName(u"tableWidget_3")
        self.tableWidget_3.setGeometry(QRect(0, -1, 1041, 225))
        self.tabWidget_3.addTab(self.tab_5, "")
        self.tab_6 = QWidget()
        self.tab_6.setObjectName(u"tab_6")
        self.tableWidget_4 = QTableWidget(self.tab_6)
        self.tableWidget_4.setObjectName(u"tableWidget_4")
        self.tableWidget_4.setGeometry(QRect(0, -1, 1041, 225))
        self.tabWidget_3.addTab(self.tab_6, "")
        self.tabWidget_5.addTab(self.tab_10, "")
        self.tab_11 = QWidget()
        self.tab_11.setObjectName(u"tab_11")
        self.tabWidget_6 = QTabWidget(self.tab_11)
        self.tabWidget_6.setObjectName(u"tabWidget_6")
        self.tabWidget_6.setGeometry(QRect(-4, -1, 1041, 252))
        self.tab_12 = QWidget()
        self.tab_12.setObjectName(u"tab_12")
        self.tableWidget_5 = QTableWidget(self.tab_12)
        self.tableWidget_5.setObjectName(u"tableWidget_5")
        self.tableWidget_5.setGeometry(QRect(0, -1, 1041, 225))
        self.tabWidget_6.addTab(self.tab_12, "")
        self.tab_13 = QWidget()
        self.tab_13.setObjectName(u"tab_13")
        self.tableWidget_6 = QTableWidget(self.tab_13)
        self.tableWidget_6.setObjectName(u"tableWidget_6")
        self.tableWidget_6.setGeometry(QRect(0, -1, 1041, 225))
        self.tabWidget_6.addTab(self.tab_13, "")
        self.tab_14 = QWidget()
        self.tab_14.setObjectName(u"tab_14")
        self.tableWidget_7 = QTableWidget(self.tab_14)
        self.tableWidget_7.setObjectName(u"tableWidget_7")
        self.tableWidget_7.setGeometry(QRect(0, -1, 1041, 225))
        self.tabWidget_6.addTab(self.tab_14, "")
        self.tabWidget_5.addTab(self.tab_11, "")
        self.groupBox_4 = QGroupBox(self.centralwidget)
        self.groupBox_4.setObjectName(u"groupBox_4")
        self.groupBox_4.setGeometry(QRect(550, 0, 511, 521))
        self.tabWidget_4 = QTabWidget(self.groupBox_4)
        self.tabWidget_4.setObjectName(u"tabWidget_4")
        self.tabWidget_4.setGeometry(QRect(0, 20, 511, 501))
        self.tab_9 = QWidget()
        self.tab_9.setObjectName(u"tab_9")
        self.graphWidget2 = PlotWidget(self.tab_9)
        self.graphWidget2.setObjectName(u"graphWidget2")
        self.graphWidget2.setGeometry(QRect(0, 0, 511, 471))
        self.tabWidget_4.addTab(self.tab_9, "")
        self.tab_7 = QWidget()
        self.tab_7.setObjectName(u"tab_7")
        self.graphWidget1 = PlotWidget(self.tab_7)
        self.graphWidget1.setObjectName(u"graphWidget1")
        self.graphWidget1.setGeometry(QRect(0, 0, 511, 471))
        self.graphWidget1.setAutoFillBackground(True)
        self.tabWidget_4.addTab(self.tab_7, "")
        self.tab_8 = QWidget()
        self.tab_8.setObjectName(u"tab_8")
        self.tabWidget_4.addTab(self.tab_8, "")
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QMenuBar(MainWindow)
        self.menubar.setObjectName(u"menubar")
        self.menubar.setGeometry(QRect(0, 0, 1068, 23))
        self.menuFile = QMenu(self.menubar)
        self.menuFile.setObjectName(u"menuFile")
        self.menuEdit = QMenu(self.menubar)
        self.menuEdit.setObjectName(u"menuEdit")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QStatusBar(MainWindow)
        self.statusbar.setObjectName(u"statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuEdit.menuAction())
        self.menuFile.addAction(self.actionNew)

        self.retranslateUi(MainWindow)

        self.tabWidget_7.setCurrentIndex(1)
        self.tabWidget.setCurrentIndex(0)
        self.tabWidget_5.setCurrentIndex(1)
        self.tabWidget_3.setCurrentIndex(3)
        self.tabWidget_6.setCurrentIndex(0)
        self.tabWidget_4.setCurrentIndex(0)


        QMetaObject.connectSlotsByName(MainWindow)
    # setupUi

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QCoreApplication.translate("MainWindow", u"MainWindow", None))
        self.actionNew.setText(QCoreApplication.translate("MainWindow", u"New", None))
        self.groupBox.setTitle(QCoreApplication.translate("MainWindow", u"Levels, Fit and Spectrum", None))
        self.paramsFileLabel.setText(QCoreApplication.translate("MainWindow", u"Parameters file:", None))
        self.paramsFileToolButton.setText(QCoreApplication.translate("MainWindow", u"...", None))
        self.previewParamsPushButton.setText(QCoreApplication.translate("MainWindow", u"Preview", None))
        self.obsEnergyFileLabel.setText(QCoreApplication.translate("MainWindow", u"Observed energies file:", None))
        self.obsEnergyFileToolButton.setText(QCoreApplication.translate("MainWindow", u"...", None))
        self.previewObsEnrPushButton.setText(QCoreApplication.translate("MainWindow", u"Preview", None))
        self.obsWavenFileLabel.setText(QCoreApplication.translate("MainWindow", u"Observed wavenumbers file:", None))
        self.obsWavenFileToolButton.setText(QCoreApplication.translate("MainWindow", u"...", None))
        self.previewWavenPushButton.setText(QCoreApplication.translate("MainWindow", u"Preview", None))
        self.obsIntFileLabel.setText(QCoreApplication.translate("MainWindow", u"Observed intensities file:", None))
        self.obsIntFileToolButton.setText(QCoreApplication.translate("MainWindow", u"...", None))
        self.previewIntPushButton.setText(QCoreApplication.translate("MainWindow", u"Preview", None))
        self.tabWidget_7.setTabText(self.tabWidget_7.indexOf(self.tab_16), QCoreApplication.translate("MainWindow", u"Datasets", None))
        self.label_10.setText(QCoreApplication.translate("MainWindow", u"Grid:", None))
        self.eqidGridRadioButton.setText(QCoreApplication.translate("MainWindow", u"Equidistant", None))
        self.noneqidGridRadioButton.setText(QCoreApplication.translate("MainWindow", u"Non-equidistant", None))
        self.solverLabel.setText(QCoreApplication.translate("MainWindow", u"Solver:", None))
        self.sincSolverRadioButton.setText(QCoreApplication.translate("MainWindow", u"Sinc", None))
        self.fourierSolverRadioButton.setText(QCoreApplication.translate("MainWindow", u"Fourier", None))
        self.fd5SolverRadioButton.setText(QCoreApplication.translate("MainWindow", u"FD5", None))
        self.nGridPointsLabel.setText(QCoreApplication.translate("MainWindow", u"Number of points:", None))
        self.rminLabel.setText(QCoreApplication.translate("MainWindow", u"Rmin:", None))
        self.rmaxLabel.setText(QCoreApplication.translate("MainWindow", u"Rmax:", None))
        self.moleculeLabel.setText(QCoreApplication.translate("MainWindow", u"Molecule and isotopologues:", None))
        self.moleculeLineEdit.setText("")
        self.refJLabel.setText(QCoreApplication.translate("MainWindow", u"Ref J:", None))
        self.refEnergyLabel.setText(QCoreApplication.translate("MainWindow", u"Ref Energy:", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab1), QCoreApplication.translate("MainWindow", u"Molecule and Grid", None))
        self.label_8.setText(QCoreApplication.translate("MainWindow", u"Sigma =", None))
        self.label_4.setText(QCoreApplication.translate("MainWindow", u"Symmetry:", None))
        self.label.setText(QCoreApplication.translate("MainWindow", u"State Label:", None))
        self.checkBox.setText(QCoreApplication.translate("MainWindow", u"e", None))
        self.label_3.setText(QCoreApplication.translate("MainWindow", u"J end =", None))
        self.label_2.setText(QCoreApplication.translate("MainWindow", u"J start =", None))
        self.label_6.setText(QCoreApplication.translate("MainWindow", u"Lambda =", None))
        self.label_9.setText(QCoreApplication.translate("MainWindow", u"J values:", None))
        self.label_7.setText(QCoreApplication.translate("MainWindow", u"Spin =", None))
        self.checkBox_2.setText(QCoreApplication.translate("MainWindow", u"f", None))
        self.label_11.setText(QCoreApplication.translate("MainWindow", u"Number of states to use =", None))
        self.addStatePushButton.setText(QCoreApplication.translate("MainWindow", u"Add State", None))
        self.removeStatePushButton.setText(QCoreApplication.translate("MainWindow", u"Remove State", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab4), QCoreApplication.translate("MainWindow", u"States", None))
        self.labelParamsLabel1.setText(QCoreApplication.translate("MainWindow", u"Label of Radial Parameters:", None))
        self.pairStatesLabel1.setText(QCoreApplication.translate("MainWindow", u"Pair States:", None))
        self.labelRotCorrLineEdit1.setPlaceholderText(QCoreApplication.translate("MainWindow", u"0.0", None))
        self.modelTypeLabel1.setText(QCoreApplication.translate("MainWindow", u"Model Type:", None))
        self.rotCorrLabel1.setText(QCoreApplication.translate("MainWindow", u"Rot Correction:", None))
        self.operatorTypeLabel1.setText(QCoreApplication.translate("MainWindow", u"Operator Type:", None))
        self.showRadParamsPushButton1.setText(QCoreApplication.translate("MainWindow", u"Show Radial Parameters", None))
        self.plotRadParamsPushButton1.setText(QCoreApplication.translate("MainWindow", u"Plot Radial Function", None))
        self.addOperatorPushButton.setText(QCoreApplication.translate("MainWindow", u"Add Operator", None))
        self.removeOperatorPushButton.setText(QCoreApplication.translate("MainWindow", u"Remove Operator", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab5), QCoreApplication.translate("MainWindow", u"Operators", None))
        self.calcEnergyLevelsCheckBox.setText(QCoreApplication.translate("MainWindow", u"Energy Levels", None))
        self.label_13.setText(QCoreApplication.translate("MainWindow", u"Energy index:", None))
        self.calcEigenenergyCheckBox.setText(QCoreApplication.translate("MainWindow", u"Eigenenergies", None))
        self.identLevelsByEnrRadioButton.setText(QCoreApplication.translate("MainWindow", u"Closest Energy", None))
        self.label_15.setText(QCoreApplication.translate("MainWindow", u"Identify Levels by:", None))
        self.label_12.setText(QCoreApplication.translate("MainWindow", u"Energy range:", None))
        self.saveEigenvectrosCheckBox.setText(QCoreApplication.translate("MainWindow", u"Save/Plot Eigvecs:", None))
        self.label_14.setText(QCoreApplication.translate("MainWindow", u"Calculate:", None))
        self.identLevelsByContrRadioButton.setText(QCoreApplication.translate("MainWindow", u"State Contr.", None))
        self.label_5.setText(QCoreApplication.translate("MainWindow", u"Diagonalization routine:", None))
        self.radioButton.setText(QCoreApplication.translate("MainWindow", u"eigh", None))
        self.radioButton_2.setText(QCoreApplication.translate("MainWindow", u"eigsh", None))
        self.previewOutputPushButton.setText(QCoreApplication.translate("MainWindow", u"Preview Output", None))
        self.solvePushButton.setText(QCoreApplication.translate("MainWindow", u"Solve", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab6), QCoreApplication.translate("MainWindow", u"Hamiltonian", None))
        self.tabWidget_7.setTabText(self.tabWidget_7.indexOf(self.tab_15), QCoreApplication.translate("MainWindow", u"Energy Levels", None))
        self.label_16.setText(QCoreApplication.translate("MainWindow", u"Optimization method:", None))
        self.numberOfIterLineEdit.setPlaceholderText(QCoreApplication.translate("MainWindow", u"1", None))
        self.label_19.setText(QCoreApplication.translate("MainWindow", u"Number of iterations:", None))
        self.toleranceLineEdit.setPlaceholderText(QCoreApplication.translate("MainWindow", u"0.1", None))
        self.label_18.setText(QCoreApplication.translate("MainWindow", u"Tolerance value:", None))
        self.label_17.setText(QCoreApplication.translate("MainWindow", u"Derivative:", None))
        self.numericalDerivRadioButton.setText(QCoreApplication.translate("MainWindow", u"Numerical", None))
        self.analyticalDerivRadioButton.setText(QCoreApplication.translate("MainWindow", u"Analytical", None))
        self.svdOptMethodRadioButton.setText(QCoreApplication.translate("MainWindow", u"SVD", None))
        self.minuitOptMethodRadioButton.setText(QCoreApplication.translate("MainWindow", u"Minuit", None))
        self.tabWidget_7.setTabText(self.tabWidget_7.indexOf(self.tab_17), QCoreApplication.translate("MainWindow", u"Fitting", None))
        self.tabWidget_7.setTabText(self.tabWidget_7.indexOf(self.tab_18), QCoreApplication.translate("MainWindow", u"Spectrum", None))
        self.groupBox_3.setTitle(QCoreApplication.translate("MainWindow", u"Data Viewer and Editor", None))
        self.tabWidget_3.setTabText(self.tabWidget_3.indexOf(self.tab_3), QCoreApplication.translate("MainWindow", u"Parameters", None))
        self.tabWidget_3.setTabText(self.tabWidget_3.indexOf(self.tab_4), QCoreApplication.translate("MainWindow", u"Obs. Energies", None))
        self.tabWidget_3.setTabText(self.tabWidget_3.indexOf(self.tab_5), QCoreApplication.translate("MainWindow", u"Obs. Waven", None))
        self.tabWidget_3.setTabText(self.tabWidget_3.indexOf(self.tab_6), QCoreApplication.translate("MainWindow", u"Obs. Intensity", None))
        self.tabWidget_5.setTabText(self.tabWidget_5.indexOf(self.tab_10), QCoreApplication.translate("MainWindow", u"Input Data", None))
        self.tabWidget_6.setTabText(self.tabWidget_6.indexOf(self.tab_12), QCoreApplication.translate("MainWindow", u"Eigenenergies", None))
        self.tabWidget_6.setTabText(self.tabWidget_6.indexOf(self.tab_13), QCoreApplication.translate("MainWindow", u"Energy Levels", None))
        self.tabWidget_6.setTabText(self.tabWidget_6.indexOf(self.tab_14), QCoreApplication.translate("MainWindow", u"Eigenvectors", None))
        self.tabWidget_5.setTabText(self.tabWidget_5.indexOf(self.tab_11), QCoreApplication.translate("MainWindow", u"Output Data", None))
        self.groupBox_4.setTitle(QCoreApplication.translate("MainWindow", u"Plots", None))
        self.tabWidget_4.setTabText(self.tabWidget_4.indexOf(self.tab_9), QCoreApplication.translate("MainWindow", u"Rad. func", None))
        self.tabWidget_4.setTabText(self.tabWidget_4.indexOf(self.tab_7), QCoreApplication.translate("MainWindow", u"Eigenvectors", None))
        self.tabWidget_4.setTabText(self.tabWidget_4.indexOf(self.tab_8), QCoreApplication.translate("MainWindow", u"All funcs", None))
        self.menuFile.setTitle(QCoreApplication.translate("MainWindow", u"File", None))
        self.menuEdit.setTitle(QCoreApplication.translate("MainWindow", u"Edit", None))
    # retranslateUi

