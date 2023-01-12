src=main.cc
src+= matrix.cc

exe=matmul

obj=${src:.cc=.o}


matmul: ${obj}
	mpicxx -o $@ $^

%.o: %.cc
	mpicxx -c $^ -o $@
clean:
	rm -rf ${obj} ${exe}
