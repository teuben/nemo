#



test1:	p100.tab ptest1
	ptest1

clean:
	rm ptest1 p100.tab tmpxyz.tab

ptest1: ptest1.c
	$(CC) $(CFLAGS) -o ptest1 ptest1.c

p100.tab:
	mkplummer - 100 | snapprint - x,y,z > p100.tab



