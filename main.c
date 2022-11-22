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
#include <math.h>
#include <time.h>

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
);

/**
 * @brief Populates a given array with the initial starting values.
 * 
 * @param input_arr (*double) pointer to a square 2D-array.
 * @param size      (int)     size of the array (one-dimension).
 */
void populate_array(
  double *input_arr,
  unsigned int size
);

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
);

/**
 * @brief Prints a given array.
 * 
 * @param arr   (*double)       pointer to array.
 * @param size  (unsigned int)  size of the array (one-dimension).
 */
void print_array(
  double *arr, 
  unsigned int size
);





int main(int argc, char const *argv[])
{
  /* code */
  return 0;
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