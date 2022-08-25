QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++11

QMAKE_CXXFLAGS += /wd"4819"

QMAKE_CXXFLAGS += /openmp

# The following define makes your compiler emit warnings if you use
# any Qt feature that has been marked deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

DEFINES += _CRT_SECURE_NO_WARNINGS __OPENMP

# You can also make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    main.cpp \
#    mainwindow.cpp \
    matsolver.cpp \
    serialmatsolver.cpp \
    temp3dfemcore.cpp \
    ../metis-5.1.0/programs/mpmetis.c \
    ../metis-5.1.0/programs/cmdline_gpmetis.c \
    ../metis-5.1.0/programs/io.c \
    tlmcore.cpp

HEADERS += \
    datatype.h \
#    mainwindow.h \
    matsolver.h \
    serialmatsolver.h \
    temp3dfemcore.h \
    ../metis-5.1.0/programs/mpmetis.h \
    tlmcore.h

FORMS +=

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../armadillo/examples/lib_win64/ -lblas_win64_MT
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../armadillo/examples/lib_win64/ -lblas_win64_MT
else:unix: LIBS += -L$$PWD/../armadillo/examples/lib_win64/ -lblas_win64_MT

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../armadillo/examples/lib_win64/ -llapack_win64_MT
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../armadillo/examples/lib_win64/ -llapack_win64_MT
else:unix: LIBS += -L$$PWD/../armadillo/examples/lib_win64/ -llapack_win64_MT

INCLUDEPATH += $$PWD/../armadillo/include \
               $$PWD/../metis-5.1.0/GKlib\
               $$PWD/../metis-5.1.0/include \
               $$PWD/../metis-5.1.0/programs \
               $$PWD/../SuperLU_5.2.1/SRC \
               $$PWD/../SuperLU_5.2.1/CBLAS \

DEPENDPATH += $$PWD/../armadillo/include

LIBS += $$PWD/../metis-5.1.0/metis.lib \
        $$PWD/../SuperLU_5.2.1/lib/SuperLU.lib \
        $$PWD/../SuperLU_5.2.1/lib/CBLAS.lib \

CONFIG +=console


# SUPERLU_MT
LIBS += \
#    $$PWD/../SuperLU_MT/LIB/SuperLU_MT_3.1.lib \
#    $$PWD/../SuperLU_MT/LIB/RELEASE/openblas.lib \

INCLUDEPATH += \
#    $$PWD/../SuperLU_MT/SRC \
#    D:/SuperLU_MT_3.1/OpenBLAS-v0.2.15-Win64-int32/include \

