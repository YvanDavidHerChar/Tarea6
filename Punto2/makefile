transendental.gif: laShit495.dat creandoElGif.py
	python creandoElGif.py
laShit495.dat: el2D
	./el2D
el2D: burgers2d.c
	cc burgers2d.c -lm -o el2D

clean: 
	rm -f *.dat *.png el2D *.gif