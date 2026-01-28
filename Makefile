
all:
	make -Cmain
	make -Cbuild 
	#make -Cplot plot 
	#make -Clatex pdf 

clean:
	make -Cbuild clean
	make -Cplot clean 
	make -Clatex clean