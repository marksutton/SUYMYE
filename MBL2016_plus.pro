#-------------------------------------------------
#
# Project created by QtCreator 2015-04-24T09:42:56
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = MBL2016
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    lineage.cpp \
    simulation.cpp \
    genus.cpp \
    qcustomplot.cpp

HEADERS  += mainwindow.h \
    lineage.h \
    simulation.h \
    genus.h \
    qcustomplot.h

FORMS    += mainwindow.ui
