CXX = g++ -std=c++11
CXXFLAG = -g -wall -std=c++11

assignment1: odd_even_parallel_sort.cpp odd_even_sequential_sort.cpp check_sorted.cpp test_data_generator.cpp
	mpic++ -std=c++11 odd_even_parallel_sort.cpp -o psort
	g++ -std=c++11 odd_even_sequential_sort.cpp -o ssort
	g++ check_sorted.cpp -o check
	g++ test_data_generator.cpp -o gen
clean:
	rm -f psort
	rm -f ssort
	rm -f check
	rm -r gen