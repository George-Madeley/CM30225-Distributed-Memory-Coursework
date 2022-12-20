/**
 * @file main.c
 * @author gm768 (gm768@bath.ac.uk)
 * @brief Computes the average of four neighboring cells for each element in a
 *        2D-array. Repeats this process until the difference between the
 *        previous average and current average is less than a given level of
 *        g_precision.
 * @version 0.1
 * @date 2023-01-06
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
  unsigned int arr_rows,
  unsigned int num_of_tests,
  unsigned int arr_type,
  float precision,
  int argc,
  char const *argv[]
);

void run_precision_tests(
  unsigned int arr_rows,
  unsigned int num_of_tests,
  unsigned int arr_type,
  float precision,
  int argc,
  char const *argv[]
);

void run_arr_size_tests(
  unsigned int arr_rows,
  unsigned int num_of_tests,
  unsigned int arr_type,
  float precision,
  int argc,
  char const *argv[]
);

int run_tests(
  unsigned int arr_rows,
  unsigned int arr_type,
  unsigned int num_of_tests,
  float precision,
  double *average_sequential_time,
  double *average_parallel_time,
  int argc,
  char const *argv[]
);

int run_test(
  unsigned int arr_rows,
  unsigned int arr_type,
  float precision,
  double *sequential_time,
  double *parallel_time,
  int argc,
  char const *argv[]
);

void compute_parallel(
  double *ptr_in_arr,
  double *ptr_out_arr,
  unsigned int size,
  float precision,
  int argc,
  char const *argv[]
);

void compute_sequentially(
  double *ptr_in_arr,
  double *ptr_out_arr,
  unsigned int size,
  double precision
);

double *calculate_averages(
  double *ptr_in_arr,
  int num_of_rows,
  int size,
  int core_id,
  int num_cores,
  int *is_precise,
  double precision
);

void populate_array(
  double *input_arr,
  unsigned int size
);

int is_expected_outcome(
  double *out_arr,
  double *exp_arr,
  unsigned int size
);

void print_array_uniform(
  double *arr, 
  unsigned int size
);

void print_array(
  double *arr, 
  unsigned int rows,
  unsigned int cols
);





int main(int argc, char const *argv[])
{
  /* code */
  unsigned int arr_rows;
  float precision;
  unsigned int test_type;
  unsigned int num_of_tests;
  unsigned int arr_type;
  if (argc >= 6) {
    arr_rows = atoi(argv[1]);
    precision = atof(argv[2]);
    test_type = atoi(argv[3]);
    num_of_tests = atoi(argv[4]);
    arr_type = atoi(argv[5]);
  } else {
    printf("Invalid Command Line Arguments");
    exit(-1);
  }

  if (test_type == 0) {
    run_static_test(arr_rows, num_of_tests, arr_type, precision, argc, argv);
  } else if (test_type == 1) {
    run_precision_tests(arr_rows, num_of_tests, arr_type, precision, argc, argv);
  } else if (test_type == 2) {
    run_arr_size_tests(arr_rows, num_of_tests, arr_type, precision, argc, argv);
  }


  return 0;
}

/**
 * @brief Runs a given number of test using the given array size and level of precision.
 * 
 * @param arr_rows      (unsigned int)  The number of rows for the array to have,
 * @param num_of_tests  (unsigned int)  The number of tests to conduct,
 * @param arr_type      (unsigned int)  The array type,
 * @param precision     (float)         The level of precision.
 */
void run_static_test(
  unsigned int arr_rows,
  unsigned int num_of_tests,
  unsigned int arr_type,
  float precision,
  int argc,
  char const *argv[]
) {
  printf("Pass/Fail,\tSequential Time (s),\tParallel Time (s)\n");
  double average_sequential_time = 0.;
  double average_parallel_time = 0.;

  int average_has_passed = run_tests(
    arr_rows,
    arr_type,
    num_of_tests,
    precision,
    &average_sequential_time,
    &average_parallel_time,
    argc,
    argv
  );

  if (average_has_passed == 0) {
    printf("FAILED,\t");
  } else {
    printf("PASSED,\t");
  }
  printf("%f,\t%f\n", average_sequential_time, average_parallel_time);
}

/**
 * @brief Runs a given number of tests for varying level of precisions. The level of precision varies logarithmically from 0.9 to the given level of precision.
 * 
 * @param arr_rows      (unsigned int)  The number of rows for the array to have,
 * @param num_of_tests  (unsigned int)  The number of tests to conduct,
 * @param arr_type      (unsigned int)  The array type,
 * @param precision     (float)         The level of precision.
 */
void run_precision_tests(
  unsigned int arr_rows,
  unsigned int num_of_tests,
  unsigned int arr_type,
  float precision,
  int argc,
  char const *argv[]
) {
  printf("Precision,\tPass/Fail,\tSequential Time (s),\tParallel Time (s)\n");
  double max_exponent = fabs(log10(precision));
  for(double exponent = 1; exponent <= max_exponent; exponent++) {
    for (double i = 9; i > 0; i--) {
      double precision = i / pow(10., exponent);

      double average_sequential_time = 0.;
      double average_parallel_time = 0.;

      int average_has_passed = run_tests(
        arr_rows,
        arr_type,
        num_of_tests,
        precision,
        &average_sequential_time,
        &average_parallel_time,
        argc,
        argv
      );

      printf("%f,\t", precision);
      if (average_has_passed == 0) {
        printf("FAILED,\t");
      } else {
        printf("PASSED,\t");
      }
      printf("%f,\t%f\n", average_sequential_time, average_parallel_time);
    }
  }
}

/**
 * @brief Runs a given number of tests for each array size. The array size varies logarithmically from 10 to the given array size.
 * 
 * @param arr_rows      (unsigned int)  The number of rows for the array to have,
 * @param num_of_tests  (unsigned int)  The number of tests to conduct,
 * @param arr_type      (unsigned int)  The array type,
 * @param precision     (float)         The level of precision.
 */
void run_arr_size_tests(
  unsigned int arr_rows,
  unsigned int num_of_tests,
  unsigned int arr_type,
  float precision,
  int argc,
  char const *argv[]
) {
  printf("Array Size,\tPass/Fail,\tSequential Time (s),\tParallel Time (s)\n");
  double max_exponent = fabs(log10(arr_rows));
  for(double exponent = 1; exponent <= max_exponent; exponent++) {
    for (unsigned int i = 1; i < 10; i ++) {
      unsigned int array_size = i * (unsigned int)pow(10, exponent);

      double average_sequential_time = 0.;
      double average_parallel_time = 0.;

      int average_has_passed = run_tests(
        arr_rows,
        arr_type,
        num_of_tests,
        precision,
        &average_sequential_time,
        &average_parallel_time,
        argc,
        argv
      );

      printf("%d,\t", array_size);
      if (average_has_passed == 0) {
        printf("FAILED,\t");
      } else {
        printf("PASSED,\t");
      }
      printf("%f,\t%f\n", average_sequential_time, average_parallel_time);
    }
  }
}

/**
 * @brief Runs a single test a given number of times and calculates the average tie it took for the sequential program and the parallel program to run respectfully.
 * 
 * @param arr_rows                (unsigned int)  The number of rows for the array to have,
 * @param arr_type                (unsigned int)  The array type,
 * @param num_of_tests            (unsigned int)  The number of tests to conduct,
 * @param precision               (float)         The level of precision,
 * @param average_sequential_time (double *)      Pointer to the average sequential time,
 * @param average_parallel_time   (double *)      Pointer to the average parallel time.
 * @return (int)  Return 1 if all tests pass (A test passes when the parallel and sequential outputs match)
 */
int run_tests(
  unsigned int arr_rows,
  unsigned int arr_type,
  unsigned int num_of_tests,
  float precision,
  double *average_sequential_time,
  double *average_parallel_time,
  int argc,
  char const *argv[]
) {
  unsigned int average_has_passed = 1;

  for(unsigned int test_num = 0; test_num < num_of_tests; test_num++) {
    double sequential_time = 0.;
    double parallel_time = 0.;
    
    int has_passed = run_test(
      arr_rows,
      arr_type,
      precision,
      &sequential_time,
      &parallel_time,
      argc,
      argv
    );
    int average_has_passed = has_passed == 0 ? 0 : average_has_passed;
    *average_sequential_time += sequential_time;
    *average_parallel_time += parallel_time;
  }
  *average_sequential_time /= num_of_tests;
  *average_parallel_time /= num_of_tests;
}

/**
 * @brief Runs a single test of both the sequential and parallel parts. Times each parts and stores them in the provided pointers.
 * 
 * @param arr_rows                (unsigned int)  The number of rows for the array to have,
 * @param arr_type                (unsigned int)  The array type,
 * @param precision               (float)         The level of precision,
 * @param sequential_time         (double *)      Pointer to the sequential time,
 * @param parallel_time           (double *)      Pointer to the parallel time.
 * 
 * @return (int)  Return 1 if test passed (Test passes when the parallel and sequential outputs match)
 */
int run_test(
  unsigned int arr_rows,
  unsigned int arr_type,
  float precision,
  double *sequential_time,
  double *parallel_time,
  int argc,
  char const *argv[]
) {
  double sequential_time = 0.;
  double parallel_time = 0.;

  double *ptr_arr = malloc((arr_rows * arr_rows) * sizeof(double));
  populate_array(ptr_arr, arr_rows);
  
  double *ptr_s_input_arr  = malloc((arr_rows * arr_rows) * sizeof(double));
  double *ptr_s_output_arr = calloc((arr_rows * arr_rows), sizeof(double));
  double *ptr_p_input_arr  = malloc((arr_rows * arr_rows) * sizeof(double));
  double *ptr_p_output_arr = calloc((arr_rows * arr_rows), sizeof(double));

  memcpy(ptr_s_input_arr , ptr_arr, (arr_rows * arr_rows) * sizeof(double));
  memcpy(ptr_p_input_arr , ptr_arr, (arr_rows * arr_rows) * sizeof(double));

  free(ptr_arr);

  clock_t sequential_time_start = clock();

  compute_sequentially(ptr_s_input_arr, ptr_s_output_arr, arr_rows, precision);

  clock_t sequential_time_end = clock();

  *sequential_time = (double)(sequential_time_end - sequential_time_start) / CLOCKS_PER_SEC;

  clock_t parallel_time_start = clock();

  compute_parallel(ptr_p_input_arr, ptr_p_output_arr, arr_rows, precision, argc, argv);

  clock_t parallel_time_end = clock();
  
  *parallel_time = (double)(parallel_time_end - parallel_time_start) / CLOCKS_PER_SEC;

  int has_passed = is_expected_outcome(ptr_p_output_arr, ptr_s_output_arr, arr_rows);

  printf("Sequential Array: \t\n");
  print_array_uniform(ptr_s_output_arr, arr_rows);
  printf("Parallel Array: \t\n");
  print_array_uniform(ptr_p_output_arr, arr_rows);

  free(ptr_s_output_arr);
  free(ptr_s_input_arr);
  free(ptr_p_input_arr);
  free(ptr_p_output_arr);

  return has_passed;
}

/**
 * @brief Populates a given array with the expected outcome values.
 * 
 * @param pptr_in_arr   (**double)       pointer to pointer of the input array,
 * @param pptr_out_arr  (**double)       pointer to pointer of the input array,
 * @param size          (unsigned int)   size of the array (one-dimension),
 * @param precision     (double)         Level of precision to reach.
 * 
 * @return (void)
 */
void compute_parallel(
  double *ptr_in_arr,
  double *ptr_out_arr,
  unsigned int size,
  float precision,
  int argc,
  char const *argv[]
) {
  unsigned int num_cores;
  int core_id;

  MPI_Status status;

  if(MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    printf("MPI_Init error\n");
  }

  MPI_Comm_size(MPI_COMM_WORLD, &num_cores);
  MPI_Comm_rank(MPI_COMM_WORLD, &core_id);

  unsigned int rows_per_core = (unsigned int)floor(size / num_cores);
  unsigned int remainder_rows = size % num_cores;

  int row_idx_start = core_id < remainder_rows ? core_id * rows_per_core * 2 : core_id * rows_per_core + remainder_rows;
  int row_idx_end = core_id < remainder_rows ? row_idx_start + rows_per_core * 2 - 1 : row_idx_start + rows_per_core - 1;
  int num_of_rows = row_idx_end + 1 - row_idx_start;
  int next_core_id = (core_id + 1) % num_cores;
  int prev_core_id = core_id != 0 ? core_id - 1 : num_cores - 1;
  int is_precise = 0;

  // Create an array to hold the allocated rows of the input array.
  // This is padded by one row at the top and another at the bottom
  // to allow for the core to access the rows above and below it that
  // are allocated to other cores.
  double sub_arr[size * (num_of_rows + 2)];
  memcpy(&sub_arr[size], &(ptr_in_arr[size * row_idx_start]), sizeof(double) * size * num_of_rows);

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
    memcpy(core_top_row, &sub_arr[size], sizeof(float) * size);
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
        MPI_Recv(&sub_arr[size * (num_of_rows + 1)], size, MPI_DOUBLE, next_core_id, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&sub_arr[0], size, MPI_DOUBLE, prev_core_id, 1, MPI_COMM_WORLD, &status);
      } else {
        MPI_Send(core_top_row, size, MPI_DOUBLE, prev_core_id, 1, MPI_COMM_WORLD);
        MPI_Send(core_bot_row, size, MPI_DOUBLE, next_core_id, 1, MPI_COMM_WORLD);
      }
      // If the core is an even ID, then it will send a row to its neighboring
      // cores.
      if (core_id % 2 == 1) {
        MPI_Send(core_top_row, size, MPI_DOUBLE, prev_core_id, 0, MPI_COMM_WORLD);
        MPI_Send(core_bot_row, size, MPI_DOUBLE, next_core_id, 0, MPI_COMM_WORLD);
      } else {
        MPI_Recv(&sub_arr[size * (num_of_rows + 1)], size, MPI_DOUBLE, next_core_id, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&sub_arr[0], size, MPI_DOUBLE, prev_core_id, 0, MPI_COMM_WORLD, &status);
      }
    }

    printf("Core %d Sub Array: \t\n", core_id);
    print_array(sub_arr, num_of_rows + 2, size);

    // Calculate the averages of the sub array and store them in a new array.
    double *avg_arr = calculate_averages(
      sub_arr,
      num_of_rows,
      size,
      core_id,
      num_cores,
      &is_precise,
      precision
    );

    // Copy the averages to the sub array.
    memcpy(&sub_arr[size], avg_arr, sizeof(double) * size * num_of_rows);

    free(avg_arr);
  }
  // Every value within the cores sub array has reached the required level
  // of precision. Therefore, the sub array can be copied to the output array.
  // To do this, the root core must gather every sub array from each core.

  // The issue is, some sub arrays are larger than others. Therefore, the
  // root must gather the larger arrays first, then the smaller ones.

  // For this to be done, the root core gathers the larger arrays first.
  // Whilst doing so, another core gathers the smaller arrays. Once completed,
  // a send/receive is sent to the root core containing the gathered smaller arrays.

  if (remainder_rows == 0) {
    // If there are no remainder rows, then all cores have the same number of rows.
    // Therefore, the root core can gather all the sub arrays in one go.
    if (core_id == 0) {
      MPI_Gather(&sub_arr[size], size * num_of_rows, MPI_DOUBLE, ptr_out_arr, size * num_of_rows, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    } else {
      MPI_Gather(&sub_arr[size], size * num_of_rows, MPI_DOUBLE, NULL, size * num_of_rows, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
  } else {
    // If there are remainder rows, then the root core must gather the larger arrays whilst,
    // another core gathers the smaller arrays. Once completed, a send/receive is sent to the
    // root core containing the gathered smaller arrays.
    if (core_id < remainder_rows) {
      // Larger Arrays are gathered by the root core.
      if (core_id == 0) {
        // Temporarily store the larger arrays in a temporary array.
        double ptr_temp_out_arr[size * num_of_rows * remainder_rows];
        MPI_Gather(&sub_arr[size], size * num_of_rows, MPI_DOUBLE, ptr_out_arr, size * num_of_rows, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        // Receives the smaller arrays from the other core. It immediately stores them in the output array.
        // This completes the gathering of all the arrays.
        MPI_Recv(&(ptr_out_arr[size * num_of_rows * remainder_rows]), size * (num_of_rows - 1) * (num_cores - remainder_rows), MPI_DOUBLE, remainder_rows, 0, MPI_COMM_WORLD, &status);
      } else {
        MPI_Gather(&sub_arr[size], size * num_of_rows, MPI_DOUBLE, NULL, size * num_of_rows, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      }
    } else {
      // Smaller Arrays are gathered by the other core.
      if (core_id == remainder_rows) {
        // Temporarily store the smaller arrays in a temporary array.
        double ptr_temp_out_arr[size * num_of_rows * (num_cores - remainder_rows)];
        MPI_Gather(&sub_arr[size], size * num_of_rows, MPI_DOUBLE, ptr_out_arr, size * num_of_rows, MPI_DOUBLE, remainder_rows, MPI_COMM_WORLD);
        // Sends the smaller arrays to the root core.
        MPI_Send(&ptr_temp_out_arr, size * num_of_rows * (num_cores - remainder_rows), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
      } else {
        MPI_Gather(&sub_arr[size], size * num_of_rows, MPI_DOUBLE, NULL, size * num_of_rows, MPI_DOUBLE, remainder_rows, MPI_COMM_WORLD);
      }
    }
  }

  // Now that the output array has been populated, MPI can finalize.

  MPI_Finalize();
}

double *calculate_averages(
  double *ptr_in_arr,
  int num_of_rows,
  int size,
  int core_id,
  int num_cores,
  int *is_precise,
  double precision
) {
  double *ptr_out_arr = calloc(size * num_of_rows, sizeof(double));

  for (int idx = size; idx < (size * num_of_rows + size); idx++) {

    int x = idx % size;
    int y = (int)floor(idx / size);

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

  return ptr_out_arr;
}

/**
 * @brief Populates a given array with the expected outcome values.
 * 
 * @param pptr_in_arr   (**double)       pointer to pointer of the input array,
 * @param pptr_out_arr  (**double)       pointer to pointer of the input array,
 * @param size          (unsigned int)   size of the array (one-dimension),
 * @param precision     (double)         Level of precision to reach.
 * 
 * @return (void)
 */
void compute_sequentially(
  double *ptr_in_arr,
  double *ptr_out_arr,
  unsigned int size,
  double precision
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
        double average = accumulator / 4;
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

/**
 * @brief Populates a given array with the initial starting values.
 * 
 * @param input_arr (*double)       pointer to a square 2D-array.
 * @param size      (unsigned int)  size of the array (one-dimension).
 */
void populate_array(
  double *input_arr,
  unsigned int size
) {
  for (unsigned int j = 0; j < size; j++) {
    for (unsigned int i = 0; i < size; i++) {
      unsigned int index = (j * size) + i;
      input_arr[index] = (double)j;
    }
  }
};

/**
 * @brief Populates a given array with the initial starting values.
 * 
 * @param input_arr (*double)       pointer to a square 2D-array.
 * @param size      (unsigned int)  size of the array (one-dimension).
 */
void populate_array(
  double *input_arr,
  unsigned int size
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

/**
 * @brief Compares two given arrays to see if their values match.
 * 
 * @param out_arr   (*double)       pointer to array to compare.
 * @param exp_arr   (*double)       pointer to array to compare.
 * @param size      (unsigned int)  size of the arrays (one-dimension).
 * 
 * @return          (int)           1 if arrays match, 0 if they do not.
 */
int is_expected_outcome(
  double *out_arr,
  double *exp_arr,
  unsigned int size
) {
  int is_same = 1;
  for (unsigned int j = 0; j < size; j++) {
    for (unsigned int i = 0; i < size; i++) {
      unsigned int index = (j * size) + i;
      if (out_arr[index] != exp_arr[index]) {
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
  unsigned int rows,
  unsigned int cols
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

/**
 * @brief Prints a given array.
 * 
 * @param arr   (*double) pointer to array.
 * @param size  (int)     size of the array (one-dimension).
 */
void print_array_uniform(
  double *arr,
  unsigned int size
) {
  print_array(arr, size, size);
};