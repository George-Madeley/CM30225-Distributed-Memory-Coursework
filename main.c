/**
 * @file main.c
 * @author gm768 (gm768@bath.ac.uk)
 * @brief Computes the average of four neighboring cells for each element in a
 *        2D-array. Repeats this process until the difference between the
 *        previous average and current average is less than a given level of
 *        g_precision.
 * 
 *        To run the program, the program requires the following arguments in
 *        the following order:
 *        1)  (unsigned int)  The size of the square 2D array,
 *        2)  (double)        The level of precision,
 *        4)  (unsigned int)  Number of tests to run (must be bigger than 0),
 *        5)  (unsigned int)  The type of test to run:
 *                  0: Run the program using he given arguments,
 *                  1: Test for each level of precision starting at 0.9 to the
 *                     specified level of precision decreasing logarithmically,
 *                  2: Test for each array size starting at 10 to the specified
 *                     array size increasing logarithmically.
 *        6)  (unsigned int)  Whether to print the input and output arrays:
 *                  0: Don't print the arrays.
 *                  1: Print the arrays.
 *       7)  (unsigned int)   Whether to run the sequential version of the
 *                            algorithm:
 * @version 0.1
 * @date 2023-01-09
 * 
 * @bug No known bugs
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include <time.h>


/**
 * @brief Runs a test with a given size, array type, number of tests, and
 *       precision.
 * 
 * @param size          The size of the array.
 * @param arr_type      The type of array to use.
 * @param num_of_tests  The number of tests to conduct.
 * @param can_print     Whether or not to print the arrays.
 * @param do_sequential Whether or not to run the sequential version of the algorithm.
 * @param precision     The precision to use.
 */
void run_static_test(
  unsigned int const size,
  unsigned int const arr_type,
  unsigned int const num_of_tests,
  unsigned int const can_print,
  unsigned int const do_sequential,
  double const precision
);

/**
 * @brief Runs a test with a given size, array type, number of tests, and
 *        precision. The precision is varied from 0.9 to the given level of
 *        precision. The average time taken for each precision is printed.
 * 
 * @param size          The size of the array.
 * @param arr_type      The type of array to use.
 * @param num_of_tests  The number of tests to conduct.
 * @param can_print     Whether or not to print the arrays.
 * @param do_sequential Whether or not to run the sequential version of the algorithm.
 * @param precision     The precision to use.
 */
void run_precision_tests(
  unsigned int const size,
  unsigned int const arr_type,
  unsigned int const num_of_tests,
  unsigned int const can_print,
  unsigned int const do_sequential,
  double const precision
);

/**
 * @brief Runs a test with a given size, array type, number of tests, and
 *       precision. The size of the array is varied from 10 to the given size.
 * 
 * @param size          The size of the array.
 * @param arr_type      The type of array to use.
 * @param num_of_tests  The number of tests to conduct.
 * @param can_print     Whether or not to print the arrays.
 * @param do_sequential Whether or not to run the sequential version of the algorithm.
 * @param precision     The precision to use.
 */
void run_arr_size_tests(
  unsigned int const size,
  unsigned int const arr_type,
  unsigned int const num_of_tests,
  unsigned int const can_print,
  unsigned int const do_sequential,
  double const precision
);

/**
 * @brief Runs a given number of tests with a given size, array type, number of
 *        tests, and precision. Each test is timed at the average time is
 *        calculated.
 * 
 * @param size                      The size of the array.
 * @param arr_type                  The type of array to use.
 * @param num_of_tests              The number of tests to conduct.
 * @param can_print                 Whether or not to print the arrays.
 * @param do_sequential             Whether or not to run the sequential version of the algorithm.
 * @param average_parallel_time     The average time taken for the parallel
 *                                  version of the algorithm.                     
 * @param average_sequential_time   The average time taken for the sequential
 *                                  version of the algorithm.
 * @param precision                 The precision to use.
 * 
 * @return int 0 if the arrays match, 1 otherwise. 
 */
int run_tests(
  unsigned int const size,
  unsigned int const arr_type,
  unsigned int const num_of_tests,
  unsigned int const can_print,
  unsigned int const do_sequential,
  double *average_parallel_time,
  double *average_sequential_time,
  double const precision
);

/**
 * @brief Runs a single test with a given size, array type, number of tests, and
 *        precision. The test is timed and the time taken is returned.
 * 
 * @param size            The size of the array.
 * @param arr_type        The type of array to use.
 * @param can_print       Whether or not to print the arrays.
 * @param do_sequential   Whether or not to run the sequential version of the algorithm.
 * @param parallel_time   The time taken for the parallel version of the algorithm.
 * @param sequential_time The time taken for the sequential version of the
 *                        algorithm.
 * @param precision       The precision to use.
 * 
 * @return int 0 if the arrays match, 1 otherwise. 
 */
int run_test(
  unsigned int const size,
  unsigned int const arr_type,
  unsigned int const can_print,
  unsigned int const do_sequential,
  double *parallel_time,
  double *sequential_time,
  double const precision
);

/**
 * @brief Computes the average of four neighboring cells for each element in a
 *        2D-array. Repeats this process until the difference between the
 *        previous average and current average is less than a given level of
 *        g_precision.
 * 
 * @param ptr_in_arr  A pointer to the input array.
 * @param ptr_out_arr A pointer to the output array.
 * @param precision   The precision to use.
 * @param can_print   Whether or not to print the arrays.
 * @param size        The size of the array.
 */
void compute_parallel(
  double *ptr_in_arr,
  double *ptr_out_arr,
  double const precision,
  unsigned int can_print,
  unsigned int const size
);

/**
 * @brief Computes the average of four neighboring cells for each element in a
 *        2D-array. Repeats this process until the difference between the
 *        previous average and current average is less than a given level of
 *        g_precision.
 * 
 * @param ptr_in_arr  A pointer to the input array.
 * @param ptr_out_arr A pointer to the output array.
 * @param size        The size of the array.
 * @param precision   The precision to use.
 */
void compute_sequentially(
  double *ptr_in_arr,
  double *ptr_out_arr,
  unsigned int const size,
  double const precision
);

/**
 * @brief Computes the average of four neighboring cells for each element in a
 *        2D-array.
 * 
 * @param ptr_in_arr  A pointer to the input array.
 * @param ptr_out_arr A pointer to the output array.
 * @param precision   The precision to use.
 * @param is_precise  A pointer to an array of flags indicating whether or not
 *                    the average for each element is precise.
 * @param num_of_rows The number of rows in the array.
 * @param size        The size of the array.
 */
void calculate_averages(
  double *ptr_in_arr,
  double *ptr_out_arr,
  double precision,
  unsigned int *is_precise,
  unsigned int const num_of_rows,
  unsigned int const size
);

/**
 * @brief Populates an array based on the given array type.
 * 
 * @param input_arr A pointer to the array to populate.
 * @param size      The size of the array.
 * @param type      The type of array to use.
 */
void populate_array(
  double *input_arr,
  unsigned int const size,
  unsigned int const type
);

/**
 * @brief Populates an array with random values.
 * 
 * @param input_arr A pointer to the array to populate.
 * @param size      The size of the array.
 */
void populate_array_type_0(
  double *input_arr,
  unsigned int const size
);

/**
 * @brief Populates the array by setting the first row and column to 1 and the
 *       rest to 0.
 * 
 * @param input_arr A pointer to the array to populate.
 * @param size      The size of the array.
 */
void populate_array_type_1(
  double *input_arr,
  unsigned int const size
);

/**
 * @brief Populates the given array by setting each elements value to its y-coordinate.
 * 
 * @param input_arr A pointer to the array to populate.
 * @param size      The size of the array.
 */
void populate_array_type_2(
  double *input_arr,
  unsigned int const size
);

/**
 * @brief Populates the array by setting the first row to 1s and the rest to 0s.
 * 
 * @param input_arr A pointer to the array to populate.
 * @param size      The size of the array.
 */
void populate_array_type_3(
  double *input_arr,
  unsigned int const size
);

/**
 * @brief Compares two arrays and returns 1 if they match, 0 otherwise.
 * 
 * @param arr_1       A pointer to the first array.
 * @param arr_2       A pointer to the second array.
 * @param size        The size of the arrays.
 * @param can_print   Whether or not to print the arrays.
 * 
 * @return int        1 if the arrays match, 0 otherwise.
 */
int do_arrays_match(
  double *arr_1,
  double *arr_2,
  unsigned int const size,
  unsigned int const can_print
);

/**
 * @brief Prints an array.
 * 
 * @param arr   A pointer to the array.
 * @param size  The size of the array.
 */
void print_array_uniform(
  double *arr, 
  unsigned int const size
);

/**
 * @brief Prints a 2D-array.
 * 
 * @param arr   A pointer to the array.
 * @param cols  The number of columns in the array.
 * @param rows  The number of rows in the array.
 */
void print_array(
  double *arr, 
  unsigned int const cols,
  unsigned int const rows
);





int main(int argc, char **argv)
{
  // Check for valid command line arguments.
  if (argc < 8) {
    printf("Invalid Command Line Arguments");
    exit(-1);
  }

  // Get command line arguments.
  unsigned int const size = (unsigned int)(abs(atoi(argv[1])));
  double const precision = (double)(fabs(atof(argv[2])));
  unsigned int const test_type = (unsigned int)(abs(atoi(argv[3]))) % 3;
  unsigned int const num_of_tests = (unsigned int)(abs(atoi(argv[4])));
  unsigned int const arr_type = (unsigned int)(abs(atoi(argv[5]))) % 4;
  unsigned int const can_print = (unsigned int)(abs(atoi(argv[6]))) % 2;
  unsigned int const do_sequential = (unsigned int)(abs(atoi(argv[7]))) % 2;

  if (size < 3) {
    printf("Invalid size. Size must be greater than 2.\n");
    exit(-1);
  }

  if (precision < 0.) {
    printf("Invalid precision. Precision must be greater than 0.\n");
    exit(-1);
  }
  // Initialize MPI and check for errors.
  int rc = MPI_Init(&argc, &argv);
  if (rc != MPI_SUCCESS) {
    printf("Error starting MPI program. Terminating.\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
    exit(-1);
  }

  // Run the appropriate test based on the test type.
  if (test_type == 0) {
    // This will run a single test a given number of times.
    run_static_test(size, arr_type, num_of_tests, can_print, do_sequential, precision);
  } else if (test_type == 1) {
    // This will run a single test with a range of precisions starting at 0.1
    // and ending at the given precision. The size will vary logarithmically.
    // A given number of tests will be run for each precision.
    run_precision_tests(size, arr_type, num_of_tests, can_print, do_sequential, precision);
  } else if (test_type == 2) {
    // This will run a single test with a range of sizes starting at 10 and
    // ending at the given size. The size will vary logarithmically. A given
    // number of tests will be run for each size.
    run_arr_size_tests(size, arr_type, num_of_tests, can_print, do_sequential, precision);
  }
  
  // Finalize MPI.
  MPI_Finalize();

  return 0;
}


void run_static_test(
  unsigned int const size,
  unsigned int const arr_type,
  unsigned int const num_of_tests,
  unsigned int const can_print,
  unsigned int const do_sequential,
  double const precision
) {

  // Get the core id so that only core 0 prints the results instead of multiple
  // cores printing the same results.
  int core_id;
  MPI_Comm_rank(MPI_COMM_WORLD, &core_id);

  // Core 0 prints the results of the test in a CSV like style. When the data
  // is written to the .out file, the data can be easily interpreted as a CSV
  // file by Excel.
  if (core_id == 0) {
    printf("Num of Cores,\tSize,\tArr Type,\tPrecision,\tPass/Fail,\tSequential Time (s),\tParallel Time (s)\n");
  }

  // Run the tests and get the average sequential and parallel times.
  double average_sequential_time = 0.;
  double average_parallel_time = 0.;

  int average_has_passed = run_tests(
    size,
    arr_type,
    num_of_tests,
    can_print,
    do_sequential,
    &average_parallel_time,
    &average_sequential_time,
    precision
  );

  // Print the results.
  int num_cores;
  MPI_Comm_size(MPI_COMM_WORLD, &num_cores);
  if (core_id == 0) {
    printf("%d,\t%d,\t%d,\t%f,\t", num_cores, size, arr_type, precision);
    if (average_has_passed == 0) {
      printf("FAILED,\t");
    } else {
      printf("PASSED,\t");
    }
    printf("%f,\t%f\n", average_sequential_time, average_parallel_time);
  }
}


void run_precision_tests(
  unsigned int const size,
  unsigned int const arr_type,
  unsigned int const num_of_tests,
  unsigned int const can_print,
  unsigned int const do_sequential,
  double const precision
) {
  // Get the core id so that only core 0 prints the results instead of multiple
  // cores printing the same results.
  int core_id;
  MPI_Comm_rank(MPI_COMM_WORLD, &core_id);

  // Core 0 prints the results of the test in a CSV like style. When the data
  // is written to the .out file, the data can be easily interpreted as a CSV
  // file by Excel.
  if (core_id == 0) {
    printf("Num of Cores,\tSize,\tArr Type,\tPrecision,\tPass/Fail,\tSequential Time (s),\tParallel Time (s)\n");
  }

  // Run the tests for a range of precisions.
  double max_exponent = fabs(log10(precision));
  for(double exponent = 1; exponent <= max_exponent; exponent++) {
    for (double i = 9; i > 0; i--) {
      // Calculate the precision for this test.
      double var_precision = i / pow(10., exponent);

      // If the precision is less than the given precision, then break.
      if (var_precision < precision) { break; }

      // Run the tests and get the average sequential and parallel times.
      double average_sequential_time = 0.;
      double average_parallel_time = 0.;

      int average_has_passed = run_tests(
        size,
        arr_type,
        num_of_tests,
        can_print,
        do_sequential,
        &average_parallel_time,
        &average_sequential_time,
        var_precision
      );

      // Print the results.
      int num_cores;
      MPI_Comm_size(MPI_COMM_WORLD, &num_cores);
      if (core_id == 0) {
        printf("%d,\t%d,\t%d,\t%f,\t", num_cores, size, arr_type, var_precision);
        if (average_has_passed == 0) {
          printf("FAILED,\t");
        } else {
          printf("PASSED,\t");
        }
        printf("%f,\t%f\n", average_sequential_time, average_parallel_time);
      }
    }
  }
}


void run_arr_size_tests(
  unsigned int const size,
  unsigned int const arr_type,
  unsigned int const num_of_tests,
  unsigned int const can_print,
  unsigned int const do_sequential,
  double const precision
) {
  // Get the core id so that only core 0 prints the results instead of multiple
  // cores printing the same results.
  int core_id;
  MPI_Comm_rank(MPI_COMM_WORLD, &core_id);

  // Core 0 prints the results of the test in a CSV like style. When the data
  // is written to the .out file, the data can be easily interpreted as a CSV
  // file by Excel.
  if (core_id == 0) {
    printf("Num of Cores,\tSize,\tArr Type,\tPrecision,\tPass/Fail,\tSequential Time (s),\tParallel Time (s)\n");
  }

  // Run the tests for a range of array sizes.
  double max_exponent = fabs(log10(size));
  for(double exponent = 1; exponent <= max_exponent; exponent++) {
    for (unsigned int i = 1; i < 10; i ++) {
      // Calculate the array size for this test.
      unsigned int array_size = i * (unsigned int)pow(10, exponent);

      // If the array size is greater than the size specified, break out of the
      // loop. This is so that the array size does not exceed the size specified.
      if (array_size > size) { break; }

      // Run the tests and get the average sequential and parallel times.
      double average_sequential_time = 0.;
      double average_parallel_time = 0.;

      int average_has_passed = run_tests(
        array_size,
        arr_type,
        num_of_tests,
        can_print,
        do_sequential,
        &average_parallel_time,
        &average_sequential_time,
        precision
      );

      // Print the results.
      int num_cores;
      MPI_Comm_size(MPI_COMM_WORLD, &num_cores);
      if (core_id == 0) {
        printf("%d,\t%d,\t%d,\t%f,\t", num_cores, array_size, arr_type, precision);
        if (average_has_passed == 0) {
          printf("FAILED,\t");
        } else {
          printf("PASSED,\t");
        }
        printf("%f,\t%f\n", average_sequential_time, average_parallel_time);
      }
    }
  }
}


int run_tests(
  unsigned int const size,
  unsigned int const arr_type,
  unsigned int const num_of_tests,
  unsigned int const can_print,
  unsigned int const do_sequential,
  double *average_parallel_time,
  double *average_sequential_time,
  double const precision
) {
  // Assume that the tests have passed. This will be overwritten if any of the
  // tests fail. If no tests fail, then the value will remain 1 hence the tests
  // have passed.
  int average_has_passed = 1;

  // Run the tests a given number of times to calculate an average runtime for
  // the sequential and parallel implementations.
  for(unsigned int test_num = 0; test_num < num_of_tests; test_num++) {
    // Run the test and get the sequential and parallel times.
    double sequential_time = 0.;
    double parallel_time = 0.;
    
    int has_passed = run_test(
      size,
      arr_type,
      can_print,
      do_sequential,
      &parallel_time,
      &sequential_time,
      precision
    );
    // If any of the tests fail, then the average has failed.
    average_has_passed = has_passed == 0 ? 0 : average_has_passed;
    // Add the sequential and parallel times to the average sequential and
    *average_sequential_time += sequential_time;
    *average_parallel_time += parallel_time;
  }
  // Divide the average sequential and parallel times by the number of tests to get
  // the average sequential and parallel times.
  *average_sequential_time /= num_of_tests;
  *average_parallel_time /= num_of_tests;
  return average_has_passed;
}


int run_test(
  unsigned int const size,
  unsigned int const arr_type,
  unsigned int const can_print,
  unsigned int const do_sequential,
  double *parallel_time,
  double *sequential_time,
  double const precision
) {

  // Get the core rank and the number of cores so we can get an more accurate
  // time for the parallel and sequential algorithms.
  int core_id;
  int num_cores;
  MPI_Comm_rank(MPI_COMM_WORLD, &core_id);
  MPI_Comm_size(MPI_COMM_WORLD, &num_cores);

  // If the number of cores is larger than the size of the array, the rows of
  // the array cannot be equally divided up amongst the cores. The algorithm
  // could use a proportion of the cores in use but this would reduce the
  // overall efficiency of the algorithm at that number of cores. Therefore,
  // return 0 stating the parallel and sequential arrays did not match
  // resulting in a failed test.
  if (num_cores > (int)size) {
    return 0;
  }

  // Create pointers to the arrays. Scatter and Gather will use this later to
  // scatter and gather data from the 0th core which contains the populated
  // data.
  double *ptr_arr = NULL;
  double *ptr_p_input_arr = NULL;
  double *ptr_p_output_arr = NULL;

  // Only core 0 populates the array incase array type 0 is used generating
  // an array filled with random numbers. In this case, each core would generate
  // a different array. Hence, only core 0 generates the array and then scatters
  // it to the other cores.

  // This also reduces the amount of memory taken up by the program in each
  // cores memory.
  if (core_id == 0) {
    ptr_arr = malloc((size * size) * sizeof(double));
    populate_array(ptr_arr, size, arr_type);
    ptr_p_input_arr = ptr_arr;
    ptr_p_output_arr = malloc((size * size) * sizeof(double));
  }
  

  // Run the parallel algorithm first as there will be less variance in the
  // start time for each core as less code is executed before the parallel
  // algorithm is run.
  double parallel_time_start = (double)clock();
  compute_parallel(ptr_p_input_arr, ptr_p_output_arr, precision, can_print, size);
  double parallel_time_end = (double)clock();
  
  // To get the true parallel time we need to find the start time of the
  // first core and the end time of the last core. This is because the
  // cores are not guaranteed to start or finish at the same time.
  double min_parallel_time_start;
  double max_parallel_time_end;
  MPI_Reduce(&parallel_time_start, &min_parallel_time_start, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&parallel_time_end, &max_parallel_time_end, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  // Only the first core will have the correct parallel time. So, it
  // will need to broadcast it to the other cores.
  if (core_id == 0) {
    *parallel_time = (double)(max_parallel_time_end - min_parallel_time_start) / CLOCKS_PER_SEC;
  }
  MPI_Bcast(parallel_time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // As the sequential algorithm is optional to reduce batch testing run time,
  // we assume that the tests have passed. This will be overwritten if the
  // sequential is run and the comparision of the sequential and parallel
  // arrays do not match. This does run the risk of not being informed when
  // the sequential and parallel arrays do not match however, rigorous testing
  // should prove that the parallel and sequential arrays match and now one
  // can batch test to calculate parallel runtime.
  int has_passed = 1;

  // If the sequential algorithm is to be run, then run it. However, only the
  // 0th core should run it as it holds the original arrays and broadcasting
  // square arrays of doubles bigger than 16,000 x 16,000 is not feasible.
  if (do_sequential && core_id == 0) {

    // Allocate memory for the input and output arrays for the sequential
    // algorithm. The arrays are copied from the original array to reduce
    // runtime.
    double *ptr_s_input_arr  = malloc((size * size) * sizeof(double));
    double *ptr_s_output_arr = malloc((size * size) * sizeof(double));
    memcpy(ptr_s_input_arr , ptr_arr, (size * size) * sizeof(double));
    memcpy(ptr_s_output_arr, ptr_arr, (size * size) * sizeof(double));

    // Run the sequential algorithm
    clock_t sequential_time_start = clock();
    compute_sequentially(ptr_s_input_arr, ptr_s_output_arr, size, precision);
    clock_t sequential_time_end = clock();

    //Calculates the sequential time
    *sequential_time = (double)(sequential_time_end - sequential_time_start) / CLOCKS_PER_SEC;

    // Broadcast the sequential time to the other cores.
    MPI_Bcast(sequential_time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Check if the sequential and parallel algorithms have produced the same
    // output. If has_passed is 1, then the arrays match and the parallel
    // output is assumed correct. Else, the parallel output is assumed
    // incorrect.
    has_passed = do_arrays_match(ptr_p_output_arr, ptr_s_output_arr, size, can_print);
    MPI_Bcast(&has_passed, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Prints the sequential arrays.
    if (can_print) {
      printf("Sequential Output:\n");
      print_array(ptr_s_output_arr, size, size);
    }

    // free the memory
    free(ptr_s_output_arr);
    free(ptr_s_input_arr);
  } else if (do_sequential){
    // If the sequential algorithm is to be run, but the core rank is not 0,
    // then it will need to receive the sequential time and whether the
    // sequential and parallel algorithms have produced the same output.
    MPI_Bcast(sequential_time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&has_passed, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }
  

  // Prints the  parallel arrays.
  if (core_id == 0 && can_print) {
    printf("Parallel Output:\n");
    print_array(ptr_p_output_arr, size, size);
  }

  // Free the memory
  if (core_id == 0) {
    free(ptr_p_output_arr);
    free(ptr_arr);
  }

  // Return whether the sequential and parallel algorithms have produced the
  // same output.
  return has_passed;
}


void compute_parallel(
  double *ptr_in_arr,
  double *ptr_out_arr,
  double precision,
  unsigned int can_print,
  unsigned int size
) {

  // Get the core rank and the number of cores so we can split the array
  // between the cores.
  int core_id;
  int num_cores;
  MPI_Comm_rank(MPI_COMM_WORLD, &core_id);
  MPI_Comm_size(MPI_COMM_WORLD, &num_cores);

  if (can_print) { printf("Core %d: Beginning Parallel Algorithm\n", core_id); }
  
  MPI_Status status;

  // Each row in the array is divided up amongst the cores. Each core has the
  // index of the row it starts and ends on. However, for arrays which cannot
  // be equally divided, the remaining rows will be distributed amongst the
  // first x number of cores; where x is the number of remaining rows.

  // Calculate the minimum number of rows that each core will be responsible
  // for.
  unsigned int min_num_of_rows = (unsigned int)floor(size / (unsigned int)num_cores);

  // Calculate the number of rows that will be left over.
  unsigned int remainder_rows = size % (unsigned int)num_cores;

  // Calculate the y-index for the row the core will start and end on. This
  // will be different depending on the number of remaining rows there are
  // and whether the core has been allocated an extra row.
  unsigned int start_row;
  unsigned int end_row;
  if ((unsigned int)core_id < remainder_rows) {
    start_row = (unsigned int)core_id * (min_num_of_rows + 1);
    end_row = start_row + min_num_of_rows;
  } else {
    start_row = (unsigned int)core_id * min_num_of_rows + remainder_rows;
    end_row = start_row + min_num_of_rows - 1;
  }
  
  // Calculate the number of rows this core will be responsible for. This is
  // different to min_num_of_rows as this will include the extra row if there
  // are remainder rows to calculate for.
  unsigned int num_of_rows = end_row + 1 - start_row;

  // Calculate the core id of the core above and below this core so the core
  // can send/receive data to/from these cores. The prev_core_id for the
  // first core will be set to the last core and the next_core_id for the
  // last core will be set to the first core. The first and last core wont
  // use the data from sent from the other core as they do not calculate
  // the averages of the top and bottom row respectfully. This is done to
  // improve the simplicity of the row sharing computation.
  int next_core_id = (core_id + 1) % num_cores;
  int prev_core_id = core_id != 0 ? core_id - 1 : num_cores - 1;

  // Calculate the number of rows that the core above and below this core
  unsigned int is_precise = 0;

  // Used for debugging
  if (can_print) {
    printf("Core %d: Start Row: %d, End Row: %d, Num of Rows: %d, next core: %d, prev core: %d\n", 
      core_id, start_row, end_row, num_of_rows, next_core_id, prev_core_id
    );
  }
  
  // Create an array to hold the allocated rows of the input array. The array
  // is padded by one row at the top and bottom. This is to store the row sent
  // to them by neighboring cores.
  double sub_arr[size * (num_of_rows + 2)];

  // Calculates the number of elements the core is responsible for.
  int num_of_ele = (int)(size * num_of_rows);

  // Before each core can share their rows, the need to receive the data from
  // the root core. However, as some cores have more rows than others, the root
  // core will need to use a variable scatter to scatter the rows amongst the
  // cores. To compute a variable scatter, the root core will require a list of
  // the number of elements each core is responsible for (called counts) and a
  // list of start indexes of where the rows can be round relative to the array
  // (called displs). Both lists should be in rank order.

  // Creates pointers for the counts and displs array.
  int *displs = NULL;
  int *counts = NULL;

  // The root core gathers the number of elements each core is responsible for.
  if (core_id == 0) {
    counts = (int *)malloc(sizeof(int) * (unsigned int)num_cores);
  }
  MPI_Gather(&num_of_ele, 1, MPI_INT, counts, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // The root core calculates where in the input array each cores elements can
  // be found.
  if (core_id == 0) {
    displs = (int *)malloc(sizeof(int) * (unsigned int)num_cores);

    // The first sub array starts at the beginning of the output array.
    displs[0] = 0;

    // The displacements of the remaining sub arrays are calculated by adding
    // the size of the previous sub array to the displacement of the previous
    // sub array.
    for (int i = 1; i < num_cores; i++) {
      displs[i] = displs[i - 1] + counts[i - 1];
    }
  }

  // The root core scatters the sub arrays to each core.
  MPI_Scatterv(
    ptr_in_arr, counts, displs, MPI_DOUBLE,
    &sub_arr[size], num_of_ele, MPI_DOUBLE, 
    0, MPI_COMM_WORLD
  );


  // Creates pointers to where the top and bottom rows can be found. This is to
  // send the top and bottom rows to the prev and next core respectfully.
  double *core_top_row = &sub_arr[size];;
  double *core_bot_row = &sub_arr[size * num_of_rows];

  if (can_print) { printf("Core %d: Beginning Iterations\n", core_id); }
  // Runs the algorithm until the required precision is met.
  while (is_precise == 0) {
    // Make the assumption that all values in the sub_arr have already
    // met the required precision.
    is_precise = 1;

    // If there is only one core, then there is no need to pass messages.
    // A loop is used to determine which core is sending data, which core
    // is receiving data, and which core is skipping to the next iteration.

    // The second and last core will receive a row from the first core whilst
    // the other cores will wait to receive. Once the second and last core have
    // received the row and the first sent it, they will become unblocked and
    // continue to the next iteration where they will send and receive the
    // row from the core above them.

    // Essentially, each core takes it in turn to send data whilst the other
    // cores wait to receive. Once all cores have received the data, they can
    // calculate the averages.

    // The First core will receive a row from the bottom core and vice versa
    // however, the cores will not use this row as the edge rows and columns
    // averages are no calculated. This is just a redundant send and receive.
    if (can_print) { printf("Core %d: Passing Rows\n", core_id); }

    if (num_cores > 1) {
      for (int core_num = 0; core_num < num_cores; core_num++) {

        int next_core_num = (core_num + 1) % num_cores;
        int prev_core_num = core_num != 0 ? core_num - 1 : num_cores - 1;

        if (core_id == prev_core_num) {
          // Wait to receive a row from the next core.
          MPI_Recv(
            &sub_arr[size * (num_of_rows + 1)], (int)size, MPI_DOUBLE,
            next_core_id, core_num, MPI_COMM_WORLD, &status
          );
        } else if (core_id == core_num) {
          // Send the top row to the previous core.
          MPI_Send(
            core_top_row, (int)size, MPI_DOUBLE,
            prev_core_id, core_num, MPI_COMM_WORLD
          );
          // Send the bottom row to the next core.
          MPI_Send(
            core_bot_row, (int)size, MPI_DOUBLE,
            next_core_id, core_num, MPI_COMM_WORLD
          );
        } else if (core_id == next_core_num) {
          // Wait to receive a row from the previous core.
          MPI_Recv(
            &sub_arr[0], (int)size, MPI_DOUBLE,
            prev_core_id, core_num, MPI_COMM_WORLD, &status
          );
        } else {
          continue;
        }
      }
    }

    if (can_print) { printf("Core %d: Calculating Averages\n", core_id); }

    // Calculate the averages of the sub array and store them in a new array.
    double avg_arr[size * num_of_rows];
    memcpy(&avg_arr, &sub_arr[size], sizeof(double) * size * num_of_rows);
    calculate_averages(
      sub_arr,
      avg_arr,
      precision,
      &is_precise,
      num_of_rows,
      size
    );
    
    
    // Copy the averages to the sub array.
    memcpy(&sub_arr[size], avg_arr, sizeof(double) * size * num_of_rows);

    // Used for Debugging.
    if (can_print) { printf("Core %d: Checking Precision\n", core_id); }

    // Some cores may have reached the required level of precision before
    // others. To ensure the cores that have reached the required level of
    // precision still communicate with their neighboring cores, the root
    // core must gather the is_precise values from each core, use a MIN
    // to check if all cores have reached the required level of precision
    // and then broadcast the result to all cores.
    unsigned int global_is_precise = 0;
    MPI_Reduce(&is_precise, &global_is_precise, 1, MPI_UNSIGNED, MPI_MIN, 0, MPI_COMM_WORLD);
    is_precise = global_is_precise == 1 ? 1 : 0;
    MPI_Bcast(&is_precise, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);


  }

  // Used for debugging
  if (can_print) { printf("Core %d: Gathering Sub Arrays\n", core_id); }

  // Every value within the cores sub array has reached the required level
  // of precision. Therefore, the sub array can be copied to the output array.
  // To do this, the root core must gather every sub array from each core.

  MPI_Gatherv(&sub_arr[size], num_of_ele, MPI_DOUBLE, ptr_out_arr, counts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // The root core frees up memory taken up by the gatherv and scatterv
  // requirements.
  if (core_id == 0) {
    free(displs);
    free(counts);
  }

  // used for debugging
  if (can_print) { printf("Core %d: Gathering Complete\n", core_id); }
}

void calculate_averages(
  double *ptr_in_arr,
  double *ptr_out_arr,
  double const precision,
  unsigned int *is_precise,
  unsigned int const num_of_rows,
  unsigned int const size
) {
  // Gets the core id and the number of cores so that the first and last
  // rows can be ignored by the first and last cores respectively.
  int core_id;
  int num_cores;
  MPI_Comm_rank(MPI_COMM_WORLD, &core_id);
  MPI_Comm_size(MPI_COMM_WORLD, &num_cores);

  // Calculates the average of the four surrounding values for each element
  // and stores the result in the output array.
  for (unsigned int idx = size; idx < (size * num_of_rows + size); idx++) {

    // Calculates the x and y coordinates of the current element So that
    // the first and last column and row can be ignored.
    unsigned int x = idx % size;
    unsigned int y = (unsigned int)floor(idx / size);


    // Ignores the first and last column
    if ((x == 0) || (x == (size - 1))) { continue; }
    // Ignores the first and last row
    if (y == 1 && core_id == 0) { continue; }
    if (y == num_of_rows && core_id == num_cores - 1) { continue; }

    // Calculates the average of the four surrounding values
    double accumulator = 0;
    double val_1 = ptr_in_arr[idx - 1];
    double val_2 = ptr_in_arr[idx + 1];
    double val_3 = ptr_in_arr[idx - size];
    double val_4 = ptr_in_arr[idx + size];
    
    accumulator += (val_1 + val_2 + val_3 + val_4);
    double average = accumulator / 4.;
    ptr_out_arr[idx - size] = average;


    // Computes the difference between the current value and the average
    double current_val = ptr_in_arr[idx];
    double difference = fabs(current_val - average);

    // Checks if the difference meets the required level of precision.
    // If it does not, then is_precise is set to 0. Else, it is set to
    // its previous value incase a previous difference did not meet the
    // required level of precision.
    *is_precise = difference >= precision ? 0 : *is_precise;
  }
}


void compute_sequentially(
  double *ptr_in_arr,
  double *ptr_out_arr,
  unsigned int const size,
  double const precision
) {
  // Loops until the required level of precision is met.
  int is_precise = 0;
  while (is_precise == 0) {
    // Sets is_precise to 1 so that it can be set to 0 if a difference
    // does not meet the required level of precision.
    is_precise = 1;

    // Loops through the array and calculates the average of the four
    // surrounding values for each element.
    for (unsigned int j = 0; j < size; j++) {
      for (unsigned int i = 0; i < size; i++) {

        // Calculates the index of the current element
        unsigned int index = (j * size) + i;

        // Ignores the first and last column and row
        if ((i == 0) || (i == (size - 1)) || (j == 0) || (j == (size - 1))) { continue; }

        // Calculates the average of the four surrounding values
        double accumulator = 0;
        double val_1 = ptr_in_arr[index - 1];
        double val_2 = ptr_in_arr[index + 1];
        double val_3 = ptr_in_arr[index - size];
        double val_4 = ptr_in_arr[index + size];
        accumulator += (val_1 + val_2 + val_3 + val_4);
        double average = accumulator / 4.;

        // Stores the average in the output array
        ptr_out_arr[index] = average;

        // Computes the difference between the current value and the average
        // and checks if the difference meets the required level of precision.
        double current_val = ptr_in_arr[index];
        double difference = fabs(current_val - average);
        if (difference >= precision) {
          is_precise = 0;
        }
      }
    }
    // Copies the output array to the input array so that the next iteration
    // can use the output array as the input array.
    memcpy(ptr_in_arr, ptr_out_arr, size * size * sizeof(double));
  }
}

void populate_array_type_0(
  double *input_arr,
  unsigned int const size
) {
  // Every element is set to a random value between 0 and 9 inclusive.

  // This is great to stress test the program.
  srand((unsigned int)time(NULL));
  for (unsigned int j = 0; j < size; j++) {
    for (unsigned int i = 0; i < size; i++) {
      unsigned int index = (j * size) + i;
      input_arr[index] = (double)(rand() % 10);
    }
  }
};

void populate_array_type_1(
  double *input_arr,
  unsigned int const size
) {
  // Sets the first row and column to 1 and the rest to 0.

  // This is just the example in the assignment description.
  for (unsigned int j = 0; j < size; j++) {
    for (unsigned int i = 0; i < size; i++) {
      unsigned int index = (j * size) + i;
      if ((j == 0) || (i == 0)) {
        input_arr[index] = 1.;
      } else {
        input_arr[index] = 0.;
      }
    }
  }
};

void populate_array_type_2(
  double *input_arr,
  unsigned int const size
) {
  // Sets each elements value to its y coordinate.

  // This is great to run only a single iteration as the average of each
  // elements neighbors will be the same every iteration. Therefore, the
  // difference between the current value and the average will always be 0
  // causing the program only to run one iteration regardless of the level
  // of precision.
  for (unsigned int j = 0; j < size; j++) {
    for (unsigned int i = 0; i < size; i++) {
      unsigned int index = (j * size) + i;
      input_arr[index] = (double)j;
    }
  }
};

void populate_array_type_3(
  double *input_arr,
  unsigned int const size
) {
  // Sets the top row to 1 and the rest to 0.

  // This is great to see how many iterations have occurred as the number of
  // iterations is proportional to the number of rows that are not 0 minus 1.
  for (unsigned int i = 0; i < size; i++) {
    input_arr[i] = 1.;
  }
};

void populate_array(
  double *input_arr,
  unsigned int const size,
  unsigned int const type
) {
  // Populates the array based on the given type.
  if (type == 0) {
    // Type 0 sets every element to a random value between 0 and 9 inclusive.
    populate_array_type_0(input_arr, size);
  } else if (type == 1) {
    // Type 1 sets the first row and column to 1 and the rest to 0.
    populate_array_type_1(input_arr, size);
  } else if (type == 2) {
    // Type 2 sets each elements value to its y coordinate.
    populate_array_type_2(input_arr, size);
  } else if (type == 3) {
    // Type 3 sets the top row to 1 and the rest to 0.
    populate_array_type_3(input_arr, size);
  }
};


int do_arrays_match(
  double *arr_1,
  double *arr_2,
  unsigned int const size,
  unsigned int const can_print
) {
  // Gets the core id so that only the root core prints instead of all cores
  // printing.
  int core_id;
  MPI_Comm_rank(MPI_COMM_WORLD, &core_id);
  if (can_print && (core_id == 0)) {
    printf("Checking if arrays match...\n");
  }

  // Assumes that the arrays match at the start. If they do not match, then
  // is_same is set to 0. Else, it is set to its previous value incase a
  // previous difference did not meet the required level of precision.
  int is_same = 1;

  // Loops through the arrays and checks if the values match.
  for (unsigned int j = 0; j < size; j++) {
    for (unsigned int i = 0; i < size; i++) {
      // Calculates the index of the current element
      unsigned int index = (j * size) + i;

      // Checks if the values match.
      if (arr_1[index] != arr_2[index]) {
        // Sets is_same to 0 if the values do not match.
        is_same = 0;
        // Prints the values if the arrays do not match for debugging.
        if (can_print && (core_id == 0)) {
          printf("NO MATCH y: %d x: %d. Expected %f, got %f.\n", j, i, arr_2[index], arr_1[index]);
        }
      }
    }
  }
  // Prints if the arrays match.
  return is_same;
}

void print_array(
  double *arr,
  unsigned int const cols,
  unsigned int const rows
) {
  // Prints the array row by row.
  for (unsigned int j = 0; j < rows; j++) {
    printf("[");
    for (unsigned int i = 0; i < cols; i++) {
      unsigned int index = (j * cols) + i;
      double value = arr[index];
      printf(" %f,", value);
    }
    printf("]\n");
  }
}


void print_array_uniform(
  double *arr,
  unsigned int const size
) {
  // Prints the array in a uniform format.
  print_array(arr, size, size);
};