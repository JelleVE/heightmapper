all: main.cpp
	g++ -o build/heightmapper main.cpp

clean:
	rm -f build/*