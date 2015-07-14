all:global.o
	gcc global.c -o global.o
	./global.o Global_2DR2-2J02 

clean:
	rm -f global.o
