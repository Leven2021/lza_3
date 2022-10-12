#include <cstdlib>
#include <fstream>
#include <iostream>
#include <chrono>


int main (int argc, char **argv) {

    int num_elements; // number of elements to be sorted
    num_elements = atoi(argv[1]); // convert command line argument to num_elements

    int elements[num_elements]; // store elements
    int sorted_elements[num_elements]; // store sorted elements

    std::ifstream input(argv[2]);
    int element;
    int i = 0;
    while (input >> element) {
        elements[i] = element;
        i++;
    }
    std::cout << "The input array is:" << std::endl;
    for (int j = 0; j < 20; j++) {
        std::cout << elements[j] << std::endl;
    }
    std::cout << "actual number of elements:" << i << std::endl;

    std::chrono::high_resolution_clock::time_point t1;
    std::chrono::high_resolution_clock::time_point t2;
    std::chrono::duration<double> time_span;
    t1 = std::chrono::high_resolution_clock::now(); // record time

    /* TODO BEGIN
        Implement sequential odd even transposition sort
        Code in this block is not a necessary. 
        Replace it with your own code.
    */
    while (true) {
        bool changed = false;
        // odd comparison
        for (int j = 0; j < num_elements - 1; j += 2) {
            if (elements[j] > elements[j+1]) {
                changed = true;
                int temp = elements[j];
                elements[j] = elements[j+1];
                elements[j+1] = temp;
            }
        }
        // even comparison
        for (int j = 1; j < num_elements - 1; j += 2) {
            if (elements[j] > elements[j+1]) {
                changed = true;
                int temp = elements[j];
                elements[j] = elements[j+1];
                elements[j+1] = temp;
            }
        }
        if (!changed) break;
    }
    for (int j = 0; j < num_elements; j++) {
        sorted_elements[j] = elements[j];
    }
    /* TODO END */

    t2 = std::chrono::high_resolution_clock::now();  
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    std::cout << "Student ID: " << "119010211" << std::endl; // replace it with your student id
    std::cout << "Name: " << "Liu Ziang" << std::endl; // replace it with your name
    std::cout << "Assignment 1" << std::endl;
    std::cout << "Sequential version" << std::endl;
    std::cout << "Run Time: " << time_span.count() << " seconds" << std::endl;
    std::cout << "Input Size: " << num_elements << std::endl;
    std::cout << "Process Number: " << 1 << std::endl;

    std::cout << "The output array is:" << std::endl;
    for (int j = 0; j < 20; j++) {
        std::cout << sorted_elements[j] << std::endl;
    }
    
    std::ofstream output(argv[2]+std::string(".seq.out"), std::ios_base::out);
    for (int i = 0; i < num_elements; i++) {
        output << sorted_elements[i] << std::endl;
    }
    
    return 0;
}


