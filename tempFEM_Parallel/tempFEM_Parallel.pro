TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        main.cpp \
        matsolver.cpp \
        parallelmatsolver.cpp \
        temp3dfemcore.cpp \
        tlmcore.cpp

QMAKE_CXXFLAGS += /wd"4819"

QMAKE_CXXFLAGS += /openmp

DEFINES += QT_DEPRECATED_WARNINGS

DEFINES += _CRT_SECURE_NO_WARNINGS __OPENMP


SOURCES += \
        main.cpp \
    ../metis-5.1.0/programs/mpmetis.c \
    ../metis-5.1.0/programs/cmdline_gpmetis.c \
    ../metis-5.1.0/programs/io.c \

# SUPERLU_MT位置
LIBS += \
    D:/SuperLU_MT/SuperLU_MT_3.1r.lib \
    D:/SuperLU_MT/openblas.lib  \

INCLUDEPATH += \
    D:/SuperLU_MT_3.1/SRC \

CONFIG += debug_and_release


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
#               $$PWD/../SuperLU_5.2.1/SRC \
#               $$PWD/../SuperLU_5.2.1/CBLAS \

DEPENDPATH += $$PWD/../armadillo/include \

LIBS += $$PWD/../metis.lib \

CONFIG += console

HEADERS += \
    datatype.h \
    matsolver.h \
    parallelmatsolver.h \
    temp3dfemcore.h \
    tlmcore.h \
    ../metis-5.1.0/programs/mpmetis.h \
