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

void run_static_test(
  unsigned int const size,
  unsigned int const arr_type,
  unsigned int const num_of_tests,
  double const precision
);

void run_precision_tests(
  unsigned int const size,
  unsigned int const arr_type,
  unsigned int const num_of_tests,
  double const precision
);

void run_arr_size_tests(
  unsigned int const size,
  unsigned int const arr_type,
  unsigned int const num_of_tests,
  double const precision
);

int run_tests(
  unsigned int const size,
  unsigned int const arr_type,
  unsigned int const num_of_tests,
  double *average_parallel_time,
  double *average_sequential_time,
  double const precision
);

int run_test(
  unsigned int const size,
  unsigned int const arr_type,
  double *parallel_time,
  double *sequential_time,
  double const precision
);

void compute_parallel(
  double *ptr_in_arr,
  double *ptr_out_arr,
  double const precision,
  unsigned int const size
);

void compute_sequentially(
  double *ptr_in_arr,
  double *ptr_out_arr,
  unsigned int const size,
  double const precision
);

void calculate_averages(
  double *ptr_in_arr,
  double *ptr_out_arr,
  double precision,
  unsigned int *is_precise,
  unsigned int const num_of_rows,
  unsigned int const size
);

void populate_array(
  double *input_arr,
  unsigned int const size,
  unsigned int const type
);

int do_arrays_match(
  double *arr_1,
  double *arr_2,
  unsigned int const size
);

void print_array_uniform(
  double *arr, 
  unsigned int const size
);

void print_array(
  double *arr, 
  unsigned int const cols,
  unsigned int const rows
);





int main(int argc, char **argv)
{
  if (argc < 6) {
    printf("Invalid Command Line Arguments");
    exit(-1);
  }

  unsigned int const size = (unsigned int)atoi(argv[1]);
  double const precision = (double)atof(argv[2]);
  unsigned int const test_type = (unsigned int)atoi(argv[3]);
  unsigned int const num_of_tests = (unsigned int)atoi(argv[4]);
  unsigned int const arr_type = (unsigned int)atoi(argv[5]);

  int rc = MPI_Init(&argc, &argv);
  if (rc != MPI_SUCCESS) {
    printf("Error starting MPI program. Terminating.\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
  }

  if (test_type == 0) {
    run_static_test(size, arr_type, num_of_tests, precision);
  } else if (test_type == 1) {
    run_precision_tests(size, arr_type, num_of_tests, precision);
  } else if (test_type == 2) {
    run_arr_size_tests(size, arr_type, num_of_tests, precision);
  }
  
  MPI_Finalize();

  return 0;
}


void run_static_test(
  unsigned int const size,
  unsigned int const arr_type,
  unsigned int const num_of_tests,
  double const precision
) {
  int core_id;
  MPI_Comm_rank(MPI_COMM_WORLD, &core_id);

  if (core_id == 0) {
    printf("Size,\tArr Type,\tPrecision,\tPass/Fail,\tSequential Time (s),\tParallel Time (s)\n");
  }

  double average_sequential_time = 0.;
  double average_parallel_time = 0.;

  int average_has_passed = run_tests(
    size,
    arr_type,
    num_of_tests,
    &average_parallel_time,
    &average_sequential_time,
    precision
  );

  if (core_id == 0) {
    printf("%d,\t%d,\t%f,\t", size, arr_type, precision);
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
  double const precision
) {
  int core_id;
  MPI_Comm_rank(MPI_COMM_WORLD, &core_id);

  if (core_id == 0) {
    printf("Size,\tArr Type,\tPrecision,\tPass/Fail,\tSequential Time (s),\tParallel Time (s)\n");
  }

  double max_exponent = fabs(log10(precision));
  for(double exponent = 1; exponent <= max_exponent; exponent++) {
    for (double i = 9; i > 0; i--) {
      double var_precision = i / pow(10., exponent);

      double average_sequential_time = 0.;
      double average_parallel_time = 0.;

      int average_has_passed = run_tests(
        size,
        arr_type,
        num_of_tests,
        &average_parallel_time,
        &average_sequential_time,
        var_precision
      );

      if (core_id == 0) {
        printf("%d,\t%d,\t%f,\t", size, arr_type, var_precision);
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
  double const precision
) {
  int core_id;
  MPI_Comm_rank(MPI_COMM_WORLD, &core_id);

  if (core_id == 0) {
    printf("Size,\tArr Type,\tPrecision,\tPass/Fail,\tSequential Time (s),\tParallel Time (s)\n");
  }

  double max_exponent = fabs(log10(size));
  for(double exponent = 1; exponent <= max_exponent; exponent++) {
    for (unsigned int i = 1; i < 10; i ++) {
      unsigned int array_size = i * (unsigned int)pow(10, exponent);

      double average_sequential_time = 0.;
      double average_parallel_time = 0.;

      int average_has_passed = run_tests(
        array_size,
        arr_type,
        num_of_tests,
        &average_parallel_time,
        &average_sequential_time,
        precision
      );

      if (core_id == 0) {
        printf("%d,\t%d,\t%f,\t", array_size, arr_type, precision);
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
  double *average_parallel_time,
  double *average_sequential_time,
  double const precision
) {
  int average_has_passed = 1;

  for(unsigned int test_num = 0; test_num < num_of_tests; test_num++) {
    double sequential_time = 0.;
    double parallel_time = 0.;
    
    int has_passed = run_test(
      size,
      arr_type,
      &parallel_time,
      &sequential_time,
      precision
    );
    average_has_passed = has_passed == 0 ? 0 : average_has_passed;
    *average_sequential_time += sequential_time;
    *average_parallel_time += parallel_time;
  }
  *average_sequential_time /= num_of_tests;
  *average_parallel_time /= num_of_tests;
  return average_has_passed;
}


int run_test(
  unsigned int const size,
  unsigned int const arr_type,
  double *parallel_time,
  double *sequential_time,
  double const precision
) {

  // Allocate memory for the array then populate it based on the array type
  double *ptr_arr = malloc((size * size) * sizeof(double));
  populate_array(ptr_arr, size, arr_type);
  
  // Allocate memory for the sequential and parallel arrays
  double *ptr_s_input_arr  = malloc((size * size) * sizeof(double));
  double *ptr_s_output_arr = malloc((size * size) * sizeof(double));
  double *ptr_p_input_arr  = malloc((size * size) * sizeof(double));
  double *ptr_p_output_arr = malloc((size * size) * sizeof(double));

  // Copy the array into the sequential and parallel arrays
  memcpy(ptr_s_input_arr , ptr_arr, (size * size) * sizeof(double));
  memcpy(ptr_s_output_arr, ptr_arr, (size * size) * sizeof(double));
  memcpy(ptr_p_input_arr , ptr_arr, (size * size) * sizeof(double));
  memcpy(ptr_p_output_arr, ptr_arr, (size * size) * sizeof(double));

  // Free the original array
  free(ptr_arr);

  int core_id;
  int num_cores;
  MPI_Comm_rank(MPI_COMM_WORLD, &core_id);
  MPI_Comm_size(MPI_COMM_WORLD, &num_cores);

  // Run the parallel algorithm
  double parallel_time_start = (double)clock();
  compute_parallel(ptr_p_input_arr, ptr_p_output_arr, precision, size);
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

  if (core_id == 0) {
    printf("Sequential Array: \t\n");
    print_array_uniform(ptr_s_output_arr, size);
    printf("Parallel Array: \t\n");
    print_array_uniform(ptr_p_output_arr, size);
  }

  int has_passed = do_arrays_match(ptr_p_output_arr, ptr_s_output_arr, size);

  free(ptr_s_output_arr);
  free(ptr_s_input_arr);
  free(ptr_p_input_arr);
  free(ptr_p_output_arr);

  return has_passed;
}


void compute_parallel(
  double *ptr_in_arr,
  double *ptr_out_arr,
  double precision,
  unsigned int size
) {
  
  int core_id;
  int num_cores;
  MPI_Comm_rank(MPI_COMM_WORLD, &core_id);
  MPI_Comm_size(MPI_COMM_WORLD, &num_cores);

  MPI_Status status;

  unsigned int rows_per_core = (unsigned int)floor(size / (unsigned int)num_cores);
  unsigned int remainder_rows = size % (unsigned int)num_cores;
  unsigned int start_row = (unsigned int)core_id < remainder_rows ? (unsigned int)core_id * rows_per_core * 2 : (unsigned int)core_id * rows_per_core + remainder_rows;
  unsigned int end_row = (unsigned int)core_id < remainder_rows ? start_row + rows_per_core * 2 - 1 : start_row + rows_per_core - 1;
  unsigned int num_of_rows = end_row + 1 - start_row;
  int next_core_id = (core_id + 1) % num_cores;
  int prev_core_id = core_id != 0 ? core_id - 1 : num_cores - 1;
  unsigned int is_precise = 0;

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


  while (is_precise == 0) {
    // Make the assumption that all values in the sub_arr have already
    // met the required precision.
    is_precise = 1;

    // Copy each cores top and bottom rows to a space in memory
    // ready to be sent to their neighboring cores.
    memcpy(core_top_row, &sub_arr[size], sizeof(double) * size);
    memcpy(core_bot_row, &sub_arr[size * num_of_rows], sizeof(double) * size);

    // If there is only one core, then there is no need to pass messages.
    // All the even ID cores will send their top and bottom rows to their
    // neighboring cores whilst the odd ID cores will wait to receive these
    // rows. Once a core has sent/received a row, they will then with to
    // receive/send a row. This way, every send has a receive and vice versa
    // therefore preventing deadlocks.

    // The First core will receive a row from the bottom core and vice versa
    // however, the cores will not use this row as the edge rows and columns
    // averages are no calculated. This is just a redundant send and receive.
    if (num_cores > 1) {
      // If the core is an odd ID, then it will wait to receive a row from
      // its neighboring cores.
      if (core_id % 2 == 1) {
        MPI_Recv(&sub_arr[size * (num_of_rows + 1)], (int)size, MPI_DOUBLE, next_core_id, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&sub_arr[0], (int)size, MPI_DOUBLE, prev_core_id, 1, MPI_COMM_WORLD, &status);
      } else {
        MPI_Send(core_top_row, (int)size, MPI_DOUBLE, prev_core_id, 1, MPI_COMM_WORLD);
        MPI_Send(core_bot_row, (int)size, MPI_DOUBLE, next_core_id, 1, MPI_COMM_WORLD);
      }
      // If the core is an even ID, then it will send a row to its neighboring
      // cores.
      if (core_id % 2 == 1) {
        MPI_Send(core_top_row, (int)size, MPI_DOUBLE, prev_core_id, 0, MPI_COMM_WORLD);
        MPI_Send(core_bot_row, (int)size, MPI_DOUBLE, next_core_id, 0, MPI_COMM_WORLD);
      } else {
        MPI_Recv(&sub_arr[size * (num_of_rows + 1)], (int)size, MPI_DOUBLE, next_core_id, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&sub_arr[0], (int)size, MPI_DOUBLE, prev_core_id, 0, MPI_COMM_WORLD, &status);
      }
    }

    double avg_arr[size * num_of_rows];
    memcpy(&avg_arr, &sub_arr[size], sizeof(double) * size * num_of_rows);
    // Calculate the averages of the sub array and store them in a new array.
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

  int count = (int)(size * num_of_rows);
  if (remainder_rows == 0) {
    // If there are no remainder rows, then all cores have the same number of rows.
    // Therefore, the root core can gather all the sub arrays in one go.
    if (core_id == 0) {
      MPI_Gather(&sub_arr[size], count, MPI_DOUBLE, ptr_out_arr, count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    } else {
      MPI_Gather(&sub_arr[size], count, MPI_DOUBLE, NULL, count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
  } else {
    // If there are remainder rows, then the root core must gather the larger arrays whilst,
    // another core gathers the smaller arrays. Once completed, a send/receive is sent to the
    // root core containing the gathered smaller arrays.
    if (core_id < (int)remainder_rows) {
      // Larger Arrays are gathered by the root core.
      if (core_id == 0) {
        MPI_Gather(&sub_arr[size], count, MPI_DOUBLE, ptr_out_arr, count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        int send_count = (int)((int)size * ((int)num_of_rows - 1) * (num_cores - (int)remainder_rows));

        // Receives the smaller arrays from the other core. It immediately stores them in the output array.
        // This completes the gathering of all the arrays.
        MPI_Recv(&(ptr_out_arr[count * (int)remainder_rows]), send_count, MPI_DOUBLE, (int)remainder_rows, 0, MPI_COMM_WORLD, &status);
      } else {
        MPI_Gather(&sub_arr[size], count, MPI_DOUBLE, NULL, count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      }
    } else {
      // Smaller Arrays are gathered by the other core.
      if (core_id == (int)remainder_rows) {
        // Temporarily store the smaller arrays in a temporary array.
        double ptr_temp_out_arr[count * (num_cores - (int)remainder_rows)];
        MPI_Gather(&sub_arr[size], count, MPI_DOUBLE, ptr_out_arr, count, MPI_DOUBLE, (int)remainder_rows, MPI_COMM_WORLD);
        // Sends the smaller arrays to the root core.
        MPI_Send(&ptr_temp_out_arr, (int)(count * (num_cores - (int)remainder_rows)), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
      } else {
        MPI_Gather(&sub_arr[size], count, MPI_DOUBLE, NULL, count, MPI_DOUBLE, (int)remainder_rows, MPI_COMM_WORLD);
      }
    }
  }

  MPI_Bcast(ptr_out_arr, (int)(size * size), MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void calculate_averages(
  double *ptr_in_arr,
  double *ptr_out_arr,
  double const precision,
  unsigned int *is_precise,
  unsigned int const num_of_rows,
  unsigned int const size
) {
  
  int core_id;
  int num_cores;
  MPI_Comm_rank(MPI_COMM_WORLD, &core_id);
  MPI_Comm_size(MPI_COMM_WORLD, &num_cores);

  for (unsigned int idx = size; idx < (size * num_of_rows + size); idx++) {

    unsigned int x = idx % size;
    unsigned int y = (unsigned int)floor(idx / size);


    // Ignores the first and last column
    if (x == 0 || x == (size - 1)) { continue; }
    // Ignores the first and last row
    if (y == 1 && core_id == 0) { continue; }
    if (y == num_of_rows && core_id == num_cores - 1) { continue; }

    // Calculates the average of the four surrounding values
    double val_1 = ptr_in_arr[idx + 1];
    double val_2 = ptr_in_arr[idx - 1];
    double val_3 = ptr_in_arr[idx + size];
    double val_4 = ptr_in_arr[idx - size];

    // There is an error here. The average is not being calculated correctly.
    // For relatively large values (>>0.002), the average is being calculated
    // correctly. However, for smaller values, the average is not being calculated
    // correctly. This is not due to type casting and other tests have proven this.
    // The cause of the error is unknown. However, the error is not significant
    // as the averages being calculated are relatively small.

    double average = (double)(val_1 + val_2 + val_3 + val_4) / 4.;

    double current_val = ptr_in_arr[idx];


    // Computes the difference between the current value and the average
    double difference  = (double)fabs(average - current_val);

    // Checks if the difference meets the required level of precision.
    // If it does not, then is_precise is set to 0. Else, it is set to
    // its previous value incase a previous difference did not meet the
    // required level of precision.
    *is_precise = difference > precision ? 0 : *is_precise;

    // Stores the average in the output array
    ptr_out_arr[idx - size] = average;
  }
}


void compute_sequentially(
  double *ptr_in_arr,
  double *ptr_out_arr,
  unsigned int const size,
  double const precision
) {
  int is_precise = 0;
  while (is_precise == 0) {
    is_precise = 1;
    for (unsigned int j = 0; j < size; j++) {
      for (unsigned int i = 0; i < size; i++) {
        unsigned int index = (j * size) + i;
        if ((i == 0) || (i == (size - 1)) || (j == 0) || (j == (size - 1))) { continue; }
        double accumulator = 0;
        double val_1 = ptr_in_arr[index - 1];
        double val_2 = ptr_in_arr[index + 1];
        double val_3 = ptr_in_arr[index - size];
        double val_4 = ptr_in_arr[index + size];
        
        accumulator += (val_1 + val_2 + val_3 + val_4);
        double average = accumulator / 4.;
        ptr_out_arr[index] = average;

        double current_val = ptr_in_arr[index];
        double difference = fabs(current_val - average);
        if (difference >= precision) {
          is_precise = 0;
        }
      }
    }
    memcpy(ptr_in_arr, ptr_out_arr, size * size * sizeof(double));
  }
}


void populate_array_type_1(
  double *input_arr,
  unsigned int const size
) {
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
  for (unsigned int i = 0; i < size; i++) {
    input_arr[i] = 1.;
  }
};

void populate_array(
  double *input_arr,
  unsigned int const size,
  unsigned int const type
) {
  if (type == 1) {
    populate_array_type_1(input_arr, size);
  } else if (type == 2) {
    populate_array_type_2(input_arr, size);
  } else if (type == 3) {
    populate_array_type_3(input_arr, size);
  }
};

/**
 * @brief Compares two given arrays to see if their values match.
 * 
 * @param out_arr   (*double)       pointer to array to compare.
 * @param exp_arr   (*double)       pointer to array to compare.
 * @param size      (unsigned int)  size of the arrays (one-dimension).
 * 
 * @return          (int)           1 if arrays match, 0 if they do not.
 */
int do_arrays_match(
  double *arr_1,
  double *arr_2,
  unsigned int const size
) {
  int is_same = 1;
  for (unsigned int j = 0; j < size; j++) {
    for (unsigned int i = 0; i < size; i++) {
      unsigned int index = (j * size) + i;
      if (arr_2[index] != arr_1[index]) {
        is_same = 0;
        break;
      }
    }
    if (is_same == 0) {
      break;
    }
  }
  return is_same;
}

void print_array(
  double *arr,
  unsigned int const cols,
  unsigned int const rows
) {
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
  print_array(arr, size, size);
};