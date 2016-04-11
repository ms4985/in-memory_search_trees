CC = gcc

build: p2randomv2.o
	$(CC)  p2randomv2.o -o build

p2randomv2.o: p2randomv2.c
	$(CC)  -c p2randomv2.c
	

.PHONY: clean
clean:
	rm -rf *.o build

.PHONY: all
all: clean build
