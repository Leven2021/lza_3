#include <mpi.h>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <chrono>


int main (int argc, char **argv) {

    MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int num_elements; // number of elements to be sorted
    
    num_elements = atoi(argv[1]); // convert command line argument to num_elements

    int elements[num_elements]; // store elements
    int sorted_elements[num_elements]; // store sorted elements

    if (rank == 0) { // read inputs from file (master process)
        std::ifstream input(argv[2]);
        int element;
        int i = 0;
        while (input >> element) {
            elements[i] = element;
            i++;
        }
        std::cout << "actual number of elements:" << i << std::endl;
        std::cout << "The input array is:" << std::endl;
        for (int j = 0; j < 20; j++) {
            std::cout << elements[j] << std::endl;
        }
    }

    std::chrono::high_resolution_clock::time_point t1;
    std::chrono::high_resolution_clock::time_point t2;
    std::chrono::duration<double> time_span;
    if (rank == 0){ 
        t1 = std::chrono::high_resolution_clock::now(); // record time
    }

    /* TODO BEGIN
        Implement parallel odd even transposition sort
        Code in this block is not a necessary. 
        Replace it with your own code.
        Useful MPI documentation: https://rookiehpc.github.io/mpi/docs
    */

    int num_my_element = (num_elements % world_size == 0) ? num_elements / world_size : num_elements / world_size + 1; // number of elements allocated to each process
    int my_element[num_my_element]; // store elements of each process

    // construct two new arrays to store old elements plus some 99999999s (larger than any other), this is for indivisible case only!!!
    int num_maxs = (num_elements % world_size == 0) ? 1 - num_elements : world_size - (num_elements % world_size);
    int new_elements[num_elements + num_maxs];
    int new_sorted_elements[num_elements + num_maxs];

    // see if the number of elements is divisible by the number of processors
    if (num_elements % world_size == 0) {
        MPI_Scatter(elements, num_my_element, MPI_INT, my_element, num_my_element, MPI_INT, 0, MPI_COMM_WORLD); // distribute elements to each process
    } else {
        for (int j = 0; j < num_elements + num_maxs; j++) {
            if (j < num_elements) {
                new_elements[j] = elements[j];
            } else {
                new_elements[j] = 99999999;
            }
        }
        MPI_Scatter(new_elements, num_my_element, MPI_INT, my_element, num_my_element, MPI_INT, 0, MPI_COMM_WORLD); // distribute elements to each process
    }

    if (num_my_element % 2 == 0) { // every process has even elements
        while (true) {
            int changed = 0;
            // odd comparison
            for (int j = 0; j < num_my_element - 1; j += 2) {
                if (my_element[j] > my_element[j+1]) {
                    changed = 1;
                    int temp = my_element[j];
                    my_element[j] = my_element[j+1];
                    my_element[j+1] = temp;
                }
            }
            // even comparison, leftovers happen
            for (int j = 1; j < num_my_element - 1; j += 2) {
                if (my_element[j] > my_element[j+1]) {
                    changed = 1;
                    int temp = my_element[j];
                    my_element[j] = my_element[j+1];
                    my_element[j+1] = temp;
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
            // boundary comparison: send the left boundary element to the anterior process, compare, and then send something back
            if (rank != 0) { // send the left boundary element to the anterior process
                MPI_Send(my_element, 1, MPI_INT, rank - 1, 7, MPI_COMM_WORLD);
                int from_left;
                MPI_Recv(&from_left, 1, MPI_INT, rank - 1, 77, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                my_element[0] = from_left;
            }
            if (rank != world_size - 1) { // compare, and then send something back
                int from_right;
                MPI_Recv(&from_right, 1, MPI_INT, rank + 1, 7, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                int left = my_element[num_my_element - 1];
                if (left > from_right) {
                    changed = 1;
                    my_element[num_my_element - 1] = from_right;
                    MPI_Send(&left, 1, MPI_INT, rank + 1, 77, MPI_COMM_WORLD);
                } else {
                    MPI_Send(&from_right, 1, MPI_INT, rank + 1, 77, MPI_COMM_WORLD);
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
            // see if it has been sorted
            int changed_all[world_size];
            MPI_Allgather(&changed, 1, MPI_INT, changed_all, 1, MPI_INT, MPI_COMM_WORLD);
            bool sorted = true;
            for (int j = 0; j < world_size; j++) {
                if (changed_all[j] == 1) {
                    sorted = false;
                    break;
                }
            }
            if (sorted) break;
        }
    } else { // every process has odd elements
        while (true) {
            int changed = 0;
            // odd comparison, leftovers happen
            for (int j = 0; j < num_my_element - 1; j += 2) {
                if (my_element[j] > my_element[j+1]) {
                    changed = 1;
                    int temp = my_element[j];
                    my_element[j] = my_element[j+1];
                    my_element[j+1] = temp;
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
            // boundary comparison
            if (rank != 0) {
                MPI_Send(my_element, 1, MPI_INT, rank - 1, 7, MPI_COMM_WORLD);
                int from_left;
                MPI_Recv(&from_left, 1, MPI_INT, rank - 1, 77, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                my_element[0] = from_left;
            }
            if (rank != world_size - 1) {
                int from_right;
                MPI_Recv(&from_right, 1, MPI_INT, rank + 1, 7, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                int left = my_element[num_my_element - 1];
                if (left > from_right) {
                    changed = 1;
                    my_element[num_my_element - 1] = from_right;
                    MPI_Send(&left, 1, MPI_INT, rank + 1, 77, MPI_COMM_WORLD);
                } else {
                    MPI_Send(&from_right, 1, MPI_INT, rank + 1, 77, MPI_COMM_WORLD);
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
            // even comparison, leftovers also happen
            for (int j = 1; j < num_my_element - 1; j += 2) {
                if (my_element[j] > my_element[j+1]) {
                    changed = 1;
                    int temp = my_element[j];
                    my_element[j] = my_element[j+1];
                    my_element[j+1] = temp;
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
            // boundary comparison
            if (rank != 0) {
                MPI_Send(my_element, 1, MPI_INT, rank - 1, 7, MPI_COMM_WORLD);
                int from_left;
                MPI_Recv(&from_left, 1, MPI_INT, rank - 1, 77, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                my_element[0] = from_left;
            }
            if (rank != world_size - 1) {
                int from_right;
                MPI_Recv(&from_right, 1, MPI_INT, rank + 1, 7, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                int left = my_element[num_my_element - 1];
                if (left > from_right) {
                    changed = 1;
                    my_element[num_my_element - 1] = from_right;
                    MPI_Send(&left, 1, MPI_INT, rank + 1, 77, MPI_COMM_WORLD);
                } else {
                    MPI_Send(&from_right, 1, MPI_INT, rank + 1, 77, MPI_COMM_WORLD);
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
            // see if it has been sorted
            int changed_all[world_size];
            MPI_Allgather(&changed, 1, MPI_INT, changed_all, 1, MPI_INT, MPI_COMM_WORLD);
            bool sorted = true;
            for (int j = 0; j < world_size; j++) {
                if (changed_all[j] == 1) {
                    sorted = false;
                    break;
                }
            }
            if (sorted) break;
        }
    }

    if (num_elements % world_size == 0) {
        MPI_Gather(my_element, num_my_element, MPI_INT, sorted_elements, num_my_element, MPI_INT, 0, MPI_COMM_WORLD); // collect result from each process
    } else {
        MPI_Gather(my_element, num_my_element, MPI_INT, new_sorted_elements, num_my_element, MPI_INT, 0, MPI_COMM_WORLD); // collect result from each process
        for (int j = 0; j < num_elements; j++) {
            sorted_elements[j] = new_sorted_elements[j]; // here we are sure that 99999999s must be at the very end!
        }
    }
    
    /* TODO END */

    if (rank == 0){ // record time (only executed in master process)
        t2 = std::chrono::high_resolution_clock::now();  
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
        std::cout << "Student ID: " << "119010211" << std::endl; // replace it with your student id
        std::cout << "Name: " << "Liu Ziang" << std::endl; // replace it with your name
        std::cout << "Assignment 1" << std::endl;
        std::cout << "MPI version" << std::endl;
        std::cout << "Run Time: " << time_span.count() << " seconds" << std::endl;
        std::cout << "Input Size: " << num_elements << std::endl;
        std::cout << "Process Number: " << world_size << std::endl;

        std::cout << "The output array is:" << std::endl;
        for (int j = 0; j < 20; j++) {
            std::cout << sorted_elements[j] << std::endl;
        }
    }

    if (rank == 0){ // write result to file (only executed in master process)
        std::ofstream output(argv[2]+std::string(".parallel.out"), std::ios_base::out);
        for (int i = 0; i < num_elements; i++) {
            output << sorted_elements[i] << std::endl;
        }
    }
    
    MPI_Finalize();
    
    return 0;
}


