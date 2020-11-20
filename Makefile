all: allalign

allalign: allalign.c
	gcc -Wall -o allalign allalign.c -lm -O3 -mcmodel=medium

clean: 
	rm -fr *~ allalign