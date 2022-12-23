/**
 * @file main.c
 * @author gm768 (gm768@bath.ac.uk)
 * @brief Computes the average of four neighboring cells for each element in a
 *        2D-array. Repeats this process until the difference between the
 *        previous average and current average is less than a given level of
 *        g_precision.
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
  unsigned int const size = (unsigned int)atoi(argv[1]);
  double const precision = (double)atof(argv[2]);
  unsigned int const test_type = (unsigned int)atoi(argv[3]);
  unsigned int const num_of_tests = (unsigned int)atoi(argv[4]);
  unsigned int const arr_type = (unsigned int)atoi(argv[5]);
  unsigned int const can_print = (unsigned int)atoi(argv[6]);
  unsigned int const do_sequential = (unsigned int)atoi(argv[7]);

  // Initialize MPI and check for errors.
  int rc = MPI_Init(&argc, &argv);
  if (rc != MPI_SUCCESS) {
    printf("Error starting MPI program. Terminating.\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
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

  // Print the header for the results.
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

  // Print the header for the results.
  if (core_id == 0) {
    printf("Num of Cores,\tSize,\tArr Type,\tPrecision,\tPass/Fail,\tSequential Time (s),\tParallel Time (s)\n");
  }

  // Run the tests for a range of precisions.
  double max_exponent = fabs(log10(precision));
  for(double exponent = 1; exponent <= max_exponent; exponent++) {
    for (double i = 9; i > 0; i--) {
      // Calculate the precision for this test.
      double var_precision = i / pow(10., exponent);

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

  // Print the header for the results.
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

  // Run the tests a given number of times.
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

  if (num_cores > (int)size) {
    return 0;
  }

  // Allocate memory for the array then populate it based on the array type
  double *ptr_arr = malloc((size * size) * sizeof(double));
  populate_array(ptr_arr, size, arr_type);
  
  // Allocate memory for the input and output arrays for the parallel algorithm.
  // The arrays are copied from the original array to reduce runtime.
  double *ptr_p_input_arr  = malloc((size * size) * sizeof(double));
  double *ptr_p_output_arr = malloc((size * size) * sizeof(double));
  memcpy(ptr_p_input_arr , ptr_arr, (size * size) * sizeof(double));
  memcpy(ptr_p_output_arr, ptr_arr, (size * size) * sizeof(double));

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
  // arrays do not match.
  int has_passed = 1;

  // If the sequential algorithm is to be run, then run it.
  if (do_sequential) {

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

    // To get a more precise sequential time we need to find the average
    // sequential time of all the cores.
    double core_sequential_time = (double)(sequential_time_end - sequential_time_start) / CLOCKS_PER_SEC;
    double summed_sequential_time;
    MPI_Reduce(&core_sequential_time, &summed_sequential_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    // Only the first core will have the correct sequential time. So, it
    // will need to broadcast it to the other cores.
    if (core_id == 0) {
      *sequential_time = summed_sequential_time / num_cores;
    }
    MPI_Bcast(sequential_time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Check if the sequential and parallel algorithms have produced the same
    // output.
    has_passed = do_arrays_match(ptr_p_output_arr, ptr_s_output_arr, size, can_print);

    // Prints the sequential arrays.
    if (core_id == 0 && can_print) {
      printf("Sequential Output:\n");
      print_array(ptr_s_output_arr, size, size);
    }

    // free the memory
    free(ptr_s_output_arr);
    free(ptr_s_input_arr);
  }

  // Prints the  parallel arrays.
  if (core_id == 0 && can_print) {
    printf("Parallel Output:\n");
    print_array(ptr_p_output_arr, size, size);
  }

  // Free the memory
  free(ptr_p_input_arr);
  free(ptr_p_output_arr);
  free(ptr_arr);

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

  // Calculate the minimum number of rows that each core will be responsible
  // for.
  unsigned int rows_per_core = (unsigned int)floor(size / (unsigned int)num_cores);
  // Calculate the number of rows that will be left over.
  unsigned int remainder_rows = size % (unsigned int)num_cores;
  // Calculate the start and end row for this core.
  unsigned int start_row = (unsigned int)core_id < remainder_rows ? (unsigned int)core_id * (rows_per_core + 1) : (unsigned int)core_id * rows_per_core + remainder_rows;
  unsigned int end_row = (unsigned int)core_id < remainder_rows ? start_row + rows_per_core : start_row + rows_per_core - 1;
  // Calculate the number of rows this core will be responsible for.
  unsigned int num_of_rows = end_row + 1 - start_row;
  // Calculate the core id of the core above and below this core.
  int next_core_id = (core_id + 1) % num_cores;
  int prev_core_id = core_id != 0 ? core_id - 1 : num_cores - 1;
  // Calculate the number of rows that the core above and below this core
  unsigned int is_precise = 0;

  if (can_print) { printf("Core %d: Start Row: %d, End Row: %d, Num of Rows: %d, next core: %d, prev core: %d\n", core_id, start_row, end_row, num_of_rows, next_core_id, prev_core_id); }
  // Create an array to hold the allocated rows of the input array.
  // This is padded by one row at the top and another at the bottom
  // to allow for the core to access the rows above and below it that
  // are allocated to other cores.
  double sub_arr[size * (num_of_rows + 2)];
  memcpy(&sub_arr[size], &(ptr_in_arr[size * start_row]), sizeof(double) * size * num_of_rows);

  // Creates two arrays to store the top and bottom rows of the core
  // that are too be sent to the adjacent cores.
  double core_top_row[size];
  double core_bot_row[size];

  if (can_print) { printf("Core %d: Beginning Iterations\n", core_id); }
  // Runs the algorithm until the required precision is met.
  while (is_precise == 0) {
    // Make the assumption that all values in the sub_arr have already
    // met the required precision.
    is_precise = 1;

    // Copy each cores top and bottom rows to a space in memory
    // ready to be sent to their neighboring cores.
    memcpy(core_top_row, &sub_arr[size], sizeof(double) * size);
    memcpy(core_bot_row, &sub_arr[size * num_of_rows], sizeof(double) * size);

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
          MPI_Recv(&sub_arr[size * (num_of_rows + 1)], (int)size, MPI_DOUBLE, next_core_id, core_num, MPI_COMM_WORLD, &status);
        } else if (core_id == core_num) {
          MPI_Send(core_top_row, (int)size, MPI_DOUBLE, prev_core_id, core_num, MPI_COMM_WORLD);
          MPI_Send(core_bot_row, (int)size, MPI_DOUBLE, next_core_id, core_num, MPI_COMM_WORLD);
        } else if (core_id == next_core_num) {
          MPI_Recv(&sub_arr[0], (int)size, MPI_DOUBLE, prev_core_id, core_num, MPI_COMM_WORLD, &status);
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

    // Some cores may have reached the required level of precision before
    // others. To ensure the cores that have reached the required level of
    // precision still communicate with their neighboring cores, the root
    // core must gather the is_precise values from each core, use a MIN
    // to check if all cores have reached the required level of precision
    // and then broadcast the result to all cores.
    if (can_print) { printf("Core %d: Checking Precision\n", core_id); }

    unsigned int global_is_precise = 0;
    MPI_Reduce(&is_precise, &global_is_precise, 1, MPI_UNSIGNED, MPI_MIN, 0, MPI_COMM_WORLD);
    is_precise = global_is_precise == 1 ? 1 : 0;
    MPI_Bcast(&is_precise, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);


  }
  // Every value within the cores sub array has reached the required level
  // of precision. Therefore, the sub array can be copied to the output array.
  // To do this, the root core must gather every sub array from each core.

  // The issue is, some sub arrays are larger than others. Therefore, the
  // root must gather the larger arrays first, then the smaller ones.

  // For this to be done, the root core gathers the larger arrays first.
  // Whilst doing so, another core gathers the smaller arrays. Once completed,
  // a send/receive is sent to the root core containing the gathered smaller arrays.

  if (can_print) { printf("Core %d: Gathering Sub Arrays\n", core_id); }



  int count = (int)(size * num_of_rows);

  

  if (remainder_rows == 0) {
    // If there are no remainder rows, then all cores have the same number of rows.
    // Therefore, the root core can gather all the sub arrays in one go.
    if (core_id == 0) {
      if (can_print) { printf("Core %d: Same Size, Gathering Receive.\n", core_id); }
      MPI_Gather(&sub_arr[size], count, MPI_DOUBLE, ptr_out_arr, count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    } else {
      if (can_print) { printf("Core %d: Same Size, Gathering Send.\n", core_id); }
      MPI_Gather(&sub_arr[size], count, MPI_DOUBLE, NULL, count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
  } else {
    // If there are remainder rows, then the root core must gather the larger arrays whilst,
    // another core gathers the smaller arrays. Once completed, a send/receive is sent to the
    // root core containing the gathered smaller arrays.
    if (core_id < (int)remainder_rows) {
      // Larger Arrays are gathered by the root core.
      if (core_id == 0) {
        if (can_print) { printf("Core %d: Large Size, Gathering Receive.\n", core_id); }
        MPI_Gather(&sub_arr[size], count, MPI_DOUBLE, ptr_out_arr, count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        int send_count = (int)((int)size * ((int)num_of_rows - 1) * (num_cores - (int)remainder_rows));

        // Receives the smaller arrays from the other core. It immediately stores them in the output array.
        // This completes the gathering of all the arrays.
        if (can_print) { printf("Core %d: ROOT Waiting for Receive.\n", core_id); }
        MPI_Recv(&(ptr_out_arr[count * (int)remainder_rows]), send_count, MPI_DOUBLE, (int)remainder_rows, 0, MPI_COMM_WORLD, &status);
      } else {
        if (can_print) { printf("Core %d: Large Size, Gathering Send.\n", core_id); }
        MPI_Gather(&sub_arr[size], count, MPI_DOUBLE, NULL, count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      }
    } else {
      // Smaller Arrays are gathered by the other core.
      if (core_id == (int)remainder_rows) {
        // Temporarily store the smaller arrays in a temporary array.
        if (can_print) { printf("Core %d: Small Size, Gathering Receive.\n", core_id); }
        double ptr_temp_out_arr[count * (num_cores - (int)remainder_rows)];
        MPI_Gather(&sub_arr[size], count, MPI_DOUBLE, ptr_temp_out_arr, count, MPI_DOUBLE, (int)remainder_rows, MPI_COMM_WORLD);
        // Sends the smaller arrays to the root core.
        if (can_print) { printf("Core %d: SmaLL Size, Waiting for Send.\n", core_id); }
        MPI_Send(&ptr_temp_out_arr, (int)(count * (num_cores - (int)remainder_rows)), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
      } else {
        if (can_print) { printf("Core %d: Small Size, Gathering Send.\n", core_id); }
        MPI_Gather(&sub_arr[size], count, MPI_DOUBLE, NULL, count, MPI_DOUBLE, (int)remainder_rows, MPI_COMM_WORLD);
      }
    }
  }

  if (can_print) { printf("Core %d: Gathering Complete\n", core_id); }
  // The root core broadcasts the output array to all cores to ensure all cores
  // have the same output array.
  MPI_Bcast(ptr_out_arr, (int)(size * size), MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (can_print) { printf("Core %d: Finished\n", core_id); }
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


void populate_array_type_1(
  double *input_arr,
  unsigned int const size
) {
  // Sets the first row and column to 1 and the rest to 0.
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
  if (type == 1) {
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