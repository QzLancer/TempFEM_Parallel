TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        main.cpp

INCLUDEPATH += $$PWD/../SuperLU_5.2.1/SRC \
               $$PWD/../SuperLU_5.2.1/CBLAS \

LIBS += $$PWD/../SuperLU_5.2.1/SuperLU/x64/Release/SuperLU.lib \
        $$PWD/../SuperLU_5.2.1/SuperLU/x64/Release/CBLAS.lib \
