seed:=../seed
ising_in:=../ising
mdmc_in:=../mdmc

CXX:=g++ -std=c++11
LD:=g++ -std=c++11
CXXFLAGS:=-Wall -O3 
CPPFLAGS:=-I random -D 'SEED="${seed}"'
LDFLAGS:=

.PHONY: all ex1 ex2 ex3 ex4 ex5 ex6 ex7 ex8 ex9 ex10 \
	clean clean_ex1 clean_ex2 clean_ex3 clean_ex4 \
	clean_ex5 clean_ex6 clean_ex7 clean_ex8 clean_ex9 clean_ex10

all: ex1 ex2 ex3 ex4 ex5 ex6 ex7 ex8 ex9 ex10

clean: clean_ex1 clean_ex2 clean_ex3 clean_ex4 clean_ex5 \
	clean_ex6 clean_ex7 clean_ex8 clean_ex9 clean_ex10

#Random number generator
random/obj/random.o: random/random.cpp random/random.h
	mkdir -p random/obj
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@

#Molecular dynamics / Monte-carlo 
mdmc/obj/mdmc.o: mdmc/mdmc.cpp mdmc/mdmc.h
	mkdir -p mdmc/obj
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@

#Monte-Carlo Ising 1D
ising/obj/ising.o: ising/ising.cpp ising/ising.h
	mkdir -p ising/obj
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@

#Exercise 1
ex1: ex1/1_1 ex1/1_2 ex1/1_3

ex1/1_1: ex1/obj/1_1.o random/obj/random.o
	${LD} ${LDFLAGS} $^ -o $@

ex1/obj/1_1.o: ex1/src/1_1.cpp random/random.h
	mkdir -p ex1/obj
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@

ex1/1_2: ex1/obj/1_2.o random/obj/random.o
	${LD} ${LDFLAGS} $^ -o $@

ex1/obj/1_2.o: ex1/src/1_2.cpp random/random.h
	mkdir -p ex1/obj
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@

ex1/1_3: ex1/obj/1_3.o random/obj/random.o
	${LD} ${LDFLAGS} $^ -o $@

ex1/obj/1_3.o: ex1/src/1_3.cpp random/random.h
	mkdir -p ex1/obj
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@

clean_ex1:
	rm -f random/obj/*
	rm -f ex1/obj/*
	rm -f ex1/1_1 ex1/1_2 ex1/1_3

#Exercise 2
ex2: ex2/2_1 ex2/2_2

ex2/2_1: ex2/obj/2_1.o random/obj/random.o
	${LD} ${LDFLAGS} $^ -o $@

ex2/obj/2_1.o: ex2/src/2_1.cpp random/random.h
	mkdir -p ex2/obj
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@

ex2/2_2: ex2/obj/2_2.o random/obj/random.o
	${LD} ${LDFLAGS} $^ -o $@

ex2/obj/2_2.o: ex2/src/2_2.cpp random/random.h
	mkdir -p ex2/obj
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@

clean_ex2:
	rm -f ex2/obj/*
	rm -f ex2/2_1 ex2/2_2
	  
#Exercise 3
ex3: ex3/3_1

ex3/3_1: ex3/obj/3_1.o random/obj/random.o
	${LD} ${LDFLAGS} $^ -o $@

ex3/obj/3_1.o: ex3/src/3_1.cpp random/random.h
	mkdir -p ex3/obj
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@

clean_ex3:
	rm -f ex3/obj/*
	rm -f ex3/3_1

#Exercise 4
ex4: ex4/4_2

ex4/4_2: ex4/obj/4_2.o random/obj/random.o mdmc/obj/mdmc.o
	${LD} ${LDFLAGS} $^ -o $@

ex4/obj/4_2.o: ex4/src/4_2.cpp mdmc/mdmc.h
	mkdir -p ex4/obj
	${CXX} ${CPPFLAGS} -I mdmc -D 'MDMC_DIR="${mdmc_in}"' ${CXXFLAGS} -c -I mdmc $< -o $@

clean_ex4:
	rm -f ex4/obj/*
	rm -f ex4/4_2

#Exercise 5
ex5: ex5/5_1

ex5/5_1: ex5/obj/5_1.o random/obj/random.o
	${LD} ${LDFLAGS} $^ -o $@

ex5/obj/5_1.o: ex5/src/5_1.cpp random/random.h
	mkdir -p ex5/obj
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@

clean_ex5:
	rm -f ex5/obj/*
	rm -f ex5/5_1

#Exercise 6
ex6: ex6/6_1

ex6/6_1: ex6/obj/6_1.o random/obj/random.o ising/obj/ising.o
	${LD} ${LDFLAGS} $^ -o $@

ex6/obj/6_1.o: ex6/src/6_1.cpp ising/ising.h
	mkdir -p ex6/obj
	${CXX} ${CPPFLAGS} -I ising -D 'IN_DIR="${ising_in}"' ${CXXFLAGS} -c -I ising $< -o $@

clean_ex6:
	rm -f ex6/obj/*
	rm -f ex6/6_1

#Exercise 7
ex7: ex7/7_2 ex7/7_4

ex7/7_2: ex7/obj/7_2.o random/obj/random.o mdmc/obj/mdmc.o
	${LD} ${LDFLAGS} $^ -o $@

ex7/obj/7_2.o: ex7/src/7_2.cpp mdmc/mdmc.h
	mkdir -p ex7/obj
	${CXX} ${CPPFLAGS} -I mdmc -D 'MDMC_DIR="${mdmc_in}"' ${CXXFLAGS} -c -I mdmc $< -o $@

ex7/7_4: ex7/obj/7_4.o random/obj/random.o mdmc/obj/mdmc.o
	${LD} ${LDFLAGS} $^ -o $@

ex7/obj/7_4.o: ex7/src/7_4.cpp mdmc/mdmc.h
	mkdir -p ex7/obj
	${CXX} ${CPPFLAGS} -I mdmc -D 'MDMC_DIR="${mdmc_in}"' ${CXXFLAGS} -c -I mdmc $< -o $@

clean_ex7:
	rm -f ex7/obj/*
	rm -f ex7/7_2 ex7/7_4

#Exercise 8
ex8: ex8/8_1 ex8/8_2

ex8/8_1: ex8/obj/8_1.o random/obj/random.o
	${LD} ${LDFLAGS} $^ -o $@

ex8/obj/8_1.o: ex8/src/8_1.cpp random/random.h
	mkdir -p ex8/obj
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@

ex8/8_2: ex8/obj/8_2.o random/obj/random.o
	${LD} ${LDFLAGS} $^ -o $@

ex8/obj/8_2.o: ex8/src/8_2.cpp random/random.h
	mkdir -p ex8/obj
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@

clean_ex8:
	rm -f ex8/obj/*
	rm -f ex8/8_1 ex8/8_2

#Exercise 9
ex9: ex9/9_1

ex9/9_1: ex9/obj/9_1.o random/obj/random.o
	${LD} ${LDFLAGS} $^ -o $@

ex9/obj/9_1.o: ex9/src/9_1.cpp random/random.h
	mkdir -p ex9/obj
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@

clean_ex9:
	rm -f ex9/obj/*
	rm -f ex9/9_1

#Exercise 10
ex10: ex10/10_1

ex10/10_1: ex10/obj/10_1.o random/obj/random.o
	mpicxx -std=c++11 ${LDFLAGS} $^ -o $@

ex10/obj/10_1.o: ex10/src/10_1.cpp random/random.h
	mkdir -p ex10/obj
	mpicxx -std=c++11 ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@

clean_ex10:
	rm -f ex10/obj/*
	rm -f ex10/10_1
