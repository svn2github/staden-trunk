all: sprun.exe

run.res.o: sprun.rc
	windres -o sprun.res.o sprun.rc 

sprun.exe: main.cpp run.res.o
	gcc -g main.cpp -o sprun.exe sprun.res.o
