all:align.o
	g++ align.cpp -o align.o
	python data_script.py

edit:edit_iPARTS2_align.o
	g++ edit_iPARTS2_align.cpp -o edit_iPARTS2_align.o


clean:
	rm -f global.o



