TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp

LIBS += -L$$PWD/../PartsEngine -lPartsEngine
LIBS += -L"C:\Program Files\Microsoft MPI\Lib\i386" -lmsmpi

INCLUDEPATH += ../PartsEngine/
INCLUDEPATH += "C:\Program Files\Microsoft MPI\Inc"

OTHER_FILES += \
    README.md \
    .gitignore
