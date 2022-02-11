all: ProyectoFmain
ProyectoFmain: ProyectoFmain.x
	./$<
gprof: gprof.x
	./$< >input-profiling.txt
	gprof ./$< >gprof-report.txt
cachegrind: cachegrind.x
	valgrind --tool=cachegrind --cachegrind-out-file=out ./$< > valgrindcgexe.txt
	cg_annotate out $< >cachegrind-report.txt
memcheck: memcheck.x
	valgrind --tool=memcheck ./$< 

testcatch.x: testcatch.cpp ProyectoF.cpp ProyectoF.h
	source ${HOME}/repos/spack/share/spack/setup-env.sh;
	spack load catch2;
	g++ -I  ${CMAKE_PREFIX_PATH} $< ProyectoF.cpp -o $@ ;
	spack unload catch2
ProyectoFmain.x: ProyectoFmain.cpp ProyectoF.cpp ProyectoF.h
	g++ -fopenmp -g -O2 -fsanitize=thread -fsanitize=undefined  $< ProyectoF.cpp -o $@

gprof.x: ProyectoFmain.cpp ProyectoF.cpp ProyectoF.h
	g++ -Wall -pg -O2 $< ProyectoF.cpp ProyectoF.h -o $@

cachegrind.x: ProyectoFmain.cpp ProyectoF.cpp ProyectoF.h
	g++ -g -O2 $< ProyectoF.cpp ProyectoF.h -o $@

memcheck.x: ProyectoFmain.cpp ProyectoF.cpp ProyectoF.h
	g++ -g -O2 $< ProyectoF.cpp ProyectoF.h -o $@

fig.png: plot.gp Simulacion1.txt Simulacion21.txt Simulacion22.txt Simulacion23.txt Simulacion3.txt Simulacion4.txt Simulacion5.txt Simulacion611.txt Simulacion612.txt Simulacion613.txt Simulacion621.txt Simulacion622.txt Simulacion623.txt Simulacion631.txt Simulacion632.txt Simulacion633.txt ProyectoFmain.x
	gnuplot plot.gp
	gnuplot plot1.gp
  
