TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt


QMAKE_CXXFLAGS += /wd"4819"


QMAKE_CXXFLAGS += /openmp

DEFINES += _CRT_SECURE_NO_WARNINGS __OPENMP


SOURCES += \
        SuperLU_MT.cpp \
        main.cpp

LIBS += \
    $$PWD/../SuperLU_MT_3.1/SuperLU_MT_3.1.lib \
    $$PWD/../OpenBLAS-v0.2.15-Win64-int32/lib/libopenblas.dll.a \

INCLUDEPATH += \
#    $$PWD \
    $$PWD/../armadillo/include \
    $$PWD/../SuperLU_MT_3.1/SRC \
    $$PWD/../OpenBLAS-v0.2.15-Win64-int32/include \

CONFIG += debug_and_release

HEADERS += \
    SuperLU_MT.h
