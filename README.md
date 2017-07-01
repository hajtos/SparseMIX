# SparseMIX
This repository holds the code of the method SpaseMIX presented in (link to add after publication).\
There are two separate implementations, one in c++ and one in python. If used as standalone programs, both take a file and output the cluster assignments.\

# c++
To run the program in c++ as a standalone, compile the program using the 'make' command, then run the program with:\
./algo \<filename\> \<no. of clusters\> \<T parameter\> [\<seed\>]\
The program will load the file, run the algorithm and print the result. The result will consist of space-seperated integers, one for each point in the dataset, representing the number of the cluster the point was assigned to(numbers start with 0).\
If you want to use the method inside another program, you can use the runAlgorithm(...) function declared in algorithm.h.\

# python
To run the program in python as a standalone, simply run:\
python main.py \<filename\> \<no. of clusters\> \<T parameter\> [\<seed\>]\
The result of the program will be the same as in the c++ implementation(see above)\
If you want to use the method inside another program, you can use the runAlgorithm(...) function declared in Algorithm.py.\

# Input format
The first line of the input file consists of 2 integers: the number of points(n) and the number of dimensions(d), seperated by a space\
Each of the next n lines corresponds to a single point of the dataset and contains space-seperated integer positions from 0 to d-1 on which this point has non-zero values.\
The input file is ended with an empty line.\
An example of an input file is located in example/mushroom.txt\
