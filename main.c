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
  float precision
);

void run_precision_tests(
  unsigned int arr_rows,
  unsigned int num_of_tests,
  unsigned int arr_type,
  float precision
);

void run_arr_size_tests(
  unsigned int arr_rows,
  unsigned int num_of_tests,
  unsigned int arr_type,
  float precision
);

int run_tests(
  unsigned int arr_rows,
  unsigned int arr_type,
  unsigned int num_of_tests,
  float precision,
  double *average_sequential_time,
  double *average_parallel_time
);

int run_test(
  unsigned int arr_rows,
  unsigned int arr_type,
  float precision,
  double *sequential_time,
  double *parallel_time
);

void compute_parallel(
  double **pptr_in_arr,
  double **pptr_out_arr,
  unsigned int size,
  float precision
);

void compute_sequentially(
  double **pptr_in_arr,
  double **pptr_out_arr,
  unsigned int size,
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

void print_array(
  double *arr, 
  unsigned int size
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
    run_static_test(arr_rows, num_of_tests, arr_type, precision);
  } else if (test_type == 1) {
    run_precision_tests(arr_rows, num_of_tests, arr_type, precision);
  } else if (test_type == 2) {
    run_arr_size_tests(arr_rows, num_of_tests, arr_type, precision);
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
  float precision
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
    &average_parallel_time
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
  float precision
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
        &average_parallel_time
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
  float precision
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
        &average_parallel_time
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
  double *average_parallel_time
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
      &parallel_time
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
  double *parallel_time
) {
  double sequential_time = 0.;
  double parallel_time = 0.;

  double *ptr_arr = malloc((arr_rows * arr_rows) * sizeof(double));
  populate_array(ptr_arr, arr_rows);
  
  double *ptr_s_input_arr  = malloc((arr_rows * arr_rows) * sizeof(double));
  double *ptr_s_output_arr = malloc((arr_rows * arr_rows) * sizeof(double));
  double *ptr_p_input_arr  = malloc((arr_rows * arr_rows) * sizeof(double));
  double *ptr_p_output_arr = malloc((arr_rows * arr_rows) * sizeof(double));

  memcpy(ptr_s_input_arr , ptr_arr, (arr_rows * arr_rows) * sizeof(double));
  memcpy(ptr_s_output_arr, ptr_arr, (arr_rows * arr_rows) * sizeof(double));
  memcpy(ptr_p_input_arr , ptr_arr, (arr_rows * arr_rows) * sizeof(double));
  memcpy(ptr_p_output_arr, ptr_arr, (arr_rows * arr_rows) * sizeof(double));
  
  double **pptr_s_input_arr  = &ptr_s_input_arr; 
  double **pptr_s_output_arr = &ptr_s_output_arr; 
  double **pptr_p_input_arr  = &ptr_p_input_arr; 
  double **pptr_p_output_arr = &ptr_p_output_arr;

  free(ptr_arr);

  clock_t sequential_time_start = clock();

  compute_sequentially(pptr_s_input_arr, pptr_s_output_arr, arr_rows, precision);
  ptr_s_input_arr = *pptr_s_input_arr;
  ptr_s_output_arr = *pptr_s_output_arr;

  clock_t sequential_time_end = clock();

  *sequential_time = (double)(sequential_time_end - sequential_time_start) / CLOCKS_PER_SEC;

  clock_t parallel_time_start = clock();

  compute_parallel(pptr_p_input_arr, pptr_p_output_arr, arr_rows, precision);
  ptr_p_input_arr = *pptr_p_input_arr;
  ptr_p_output_arr = *pptr_p_output_arr;

  clock_t parallel_time_end = clock();
  
  *parallel_time = (double)(parallel_time_end - parallel_time_start) / CLOCKS_PER_SEC;

  if (print_array == 1) {
    printf("Input Array:\n");
    print_array(ptr_p_input_arr, arr_rows);
    printf("\nOutput Array:\n");
    print_array(ptr_p_output_arr, arr_rows);  
    printf("\nExpected Array:\n");
    print_array(ptr_s_output_arr, arr_rows);
  }

  int has_passed = is_expected_outcome(ptr_p_output_arr, ptr_s_output_arr, arr_rows);

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
  double **pptr_in_arr,
  double **pptr_out_arr,
  unsigned int size,
  float precision
) {

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
  double **pptr_in_arr,
  double **pptr_out_arr,
  unsigned int size,
  double precision
) {
  double *in_arr = *pptr_in_arr;
  double *out_arr = *pptr_out_arr;
  int is_precise = 0;
  while (is_precise == 0) {
    is_precise = 1;
    for (unsigned int j = 0; j < size; j++) {
      for (unsigned int i = 0; i < size; i++) {
        unsigned int index = (j * size) + i;
        if ((i == 0) || (i == (size - 1)) || (j == 0) || (j == (size - 1))) { continue; }
        double accumulator = 0;
        double val_1 = in_arr[index - 1];
        double val_2 = in_arr[index + 1];
        double val_3 = in_arr[index - size];
        double val_4 = in_arr[index + size];
        accumulator += (val_1 + val_2 + val_3 + val_4);
        double average = accumulator / 4;
        out_arr[index] = average;

        double current_val = in_arr[index];
        double difference = fabs(current_val - average);
        if (difference >= precision) {
          is_precise = 0;
        }
      }
    }
    double *temp_arr = out_arr;
    out_arr = in_arr;
    in_arr = temp_arr;
  }
  *pptr_in_arr = out_arr;
  *pptr_out_arr = in_arr;
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
      if ((j == 0) || (i == 0)) {
        input_arr[index] = 1.;
      } else {
        input_arr[index] = 0.;
      }
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

/**
 * @brief Prints a given array.
 * 
 * @param arr   (*double) pointer to array.
 * @param size  (int)     size of the array (one-dimension).
 */
void print_array(
  double *arr,
  unsigned int size
) {
  for (unsigned int j = 0; j < size; j++) {
    printf("[");
    for (unsigned int i = 0; i < size; i++) {
      unsigned int index = (j * size) + i;
      double value = arr[index];
      printf(" %f,", value);
    }
    printf("]\n");
  }
};