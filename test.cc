// Copyright 2012 Alexandre Yukio Harano
//
// Licensed under the ImageMagick License (the "License"); you may not use
// this file except in compliance with the License.  You may obtain a copy
// of the License at
//
//   http://www.imagemagick.org/script/license.php
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
// WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  See the
// License for the specific language governing permissions and limitations
// under the License.
//
// Author: Alexandre Yukio Harano <ayharano AT ime DOT usp DOT br>

// This file contains the implementation of an algorithm tester for both
// erosion and dilation for 2D images.

#include <cmath>
#include <fstream>

#include "border.h"
#include "naive.h"
#include "matrix.h"

#include "img_2d.h"
#include "test.h"

namespace {

static imaging::BoundingBox *bb = NULL;
static int half_se_length = 0;
static imaging::Position *p = NULL;

bool MakeP() {
  char n = imaging::Dimension::number();
  if (n == 0) n = imaging::Dimension::Set(2);
  if (n != 2) return false;
  if (::p != NULL) return false;
  ::p = new imaging::Position();
  if (::p == NULL) return false;
  return (::p)->SetAsOrigin();
}

bool SetBB() {
  char n = imaging::Dimension::number();
  bool ok_so_far = true;
  if (n == 0) n = imaging::Dimension::Set(2);
  if (n != 2) return false;
  if (::bb != NULL) return false;
  if (::half_se_length < 1 || ::half_se_length > 4) return false;
  imaging::Position lower;
  imaging::Position upper;
  ok_so_far = lower.set_value(0, -1*(::half_se_length));
  if (!ok_so_far) return ok_so_far;
  ok_so_far = lower.set_value(1, -1*(::half_se_length));
  if (!ok_so_far) return ok_so_far;
  ok_so_far = upper.set_value(0, ::half_se_length);
  if (!ok_so_far) return ok_so_far;
  ok_so_far = upper.set_value(1, ::half_se_length);
  if (!ok_so_far) return ok_so_far;
  ::bb = new imaging::BoundingBox(lower, upper);
  return ::bb != NULL;
}

bool MakeSEWithOriginOnly(imaging::binary::StructuringElement **se_to_be) {
  bool ok_so_far = true;
  if (se_to_be == NULL) return false;
  if (*se_to_be != NULL) return false;
  if (::bb == NULL) {
    ok_so_far = SetBB();
    if (::bb == NULL) ok_so_far = false;
    if (!ok_so_far) return ok_so_far;
  }
  if (::p == NULL) {
    ok_so_far = MakeP();
    if (::p == NULL) ok_so_far = false;
    if (!ok_so_far) return ok_so_far;
  }
  *se_to_be = new imaging::binary::StructuringElement(*bb, true);
  if (*se_to_be == NULL) return false;
  ok_so_far = (::p)->SetAsOrigin();
  if (!ok_so_far) return ok_so_far;
  ok_so_far = (*se_to_be)->set_value(*(::p), true);
  return ok_so_far;
}

bool SetRandomSE(imaging::binary::StructuringElement **se_to_be) {
  imaging::binary::StructuringElement *current_se = NULL;
  bool ok_so_far = true;
  int x = 0;
  int y = 0;
  ok_so_far = MakeSEWithOriginOnly(se_to_be);
  if (*se_to_be == NULL) ok_so_far = false;
  if (!ok_so_far) return ok_so_far;
  current_se = *se_to_be;
  bool value = true;
  for (x = -(::half_se_length); ok_so_far && x <= ::half_se_length; ++x) {
    ok_so_far = (::p)->set_value(0, x);
    if (!ok_so_far) continue;
    for (y = -(::half_se_length); ok_so_far && y <= ::half_se_length; ++y) {
      if (x == 0 && y == 0) continue;
      ok_so_far = (::p)->set_value(1, y);
      if (!ok_so_far) continue;
      value = (rand()%2) == 1;
      ok_so_far = current_se->set_value(*(::p), value);
      if (!ok_so_far) continue;
    }
  }
  return ok_so_far;
}

} // namespace

int tester(
    const bool debug,
    const std::string &file_path,
    const std::string &counter_data_prefix,
    const int number_of_se,
    const bool be_verbose,
    const bool do_save,
    const bool image_info,
    const long seed,
    const int total_algorithms,
    const bool algorithms[]) {
  std::ostream &debug_output = std::cout;
  std::vector<imaging::binary::StructuringElement*> actual_se;
  std::vector< imaging::ImagePositionIndex >
      *algorithm_determinate_border_comparison_counter = NULL;
  std::vector< imaging::ImagePositionIndex >
      *algorithm_insert_new_candidate_comparison_counter = NULL;
  std::vector< imaging::ImagePositionIndex >
      *algorithm_insert_new_candidate_memory_access_counter = NULL;
  std::vector< imaging::ImagePositionIndex >
      *algorithm_number_of_elements_in_border = NULL;
  std::vector< imaging::ImagePositionIndex >
      *algorithm_remove_candidate_comparison_counter = NULL;
  std::vector< imaging::ImagePositionIndex >
      *algorithm_remove_candidate_memory_access_counter = NULL;
  std::vector<int> background;
  imaging::grayscale::Image *current_output = NULL;
  imaging::binary::StructuringElement *current_se = NULL;
  std::vector<std::string>::const_iterator current_string;
  imaging::binary::morphology::Transform *current_transform = NULL;
  std::vector<std::string> data_counter_variable_name;
  int delta = 0;
  std::vector< std::vector< imaging::ImagePositionIndex >* >
      determinate_border_comparison_counter;
  const std::vector<imaging::ImagePositionIndex> empty_data;
  double end = 0.;
  std::vector<int> foreground;
  int height = 0;
  int i = 0;
  int i_se = 0;
  imaging::binary::Image *image = NULL;
  imaging::binary::Image *image_d = NULL;
  imaging::binary::Image *image_e = NULL;
  std::vector< std::vector< imaging::ImagePositionIndex >* >
      insert_new_candidate_comparison_counter;
  std::vector< std::vector< imaging::ImagePositionIndex >* >
      insert_new_candidate_memory_access_counter;
  int input_background = 0;
  int input_foreground = 0;
  bool iterate_for_save_or_info = true;
  int iteration = 0;
  int m = 0;
  int n = 0;
  std::vector< std::vector< imaging::ImagePositionIndex >* >
      number_of_elements_in_border;
  bool ok_so_far = true;
  std::vector<imaging::grayscale::Image *> output;
  std::vector<std::ofstream*> output_text;
  FILE *output_counter_data_file = NULL;
  int output_position_value = 0;
  bool position_value = true;
  std::vector< std::vector< imaging::ImagePositionIndex >* >
      remove_candidate_comparison_counter;
  std::vector< std::vector< imaging::ImagePositionIndex >* >
      remove_candidate_memory_access_counter;
  int result = 0;
  double start = 0.;
  std::vector<double> times;
  imaging::grayscale::Image *unpadded_output = NULL;
  int x = 0;
  int y = 0;
  int width = 0;
  // Vector data initializer.
  for (i = 0; i < total_algorithms; ++i) {
    // Algorithm verifier data.
    background.push_back(0);
    foreground.push_back(0);
    output.push_back(NULL);
    output_text.push_back(NULL);
    times.push_back(0.);
    // Algorithm complexity analyzer data.
    determinate_border_comparison_counter.push_back(NULL);
    insert_new_candidate_comparison_counter.push_back(NULL);
    insert_new_candidate_memory_access_counter.push_back(NULL);
    remove_candidate_comparison_counter.push_back(NULL);
    remove_candidate_memory_access_counter.push_back(NULL);
    number_of_elements_in_border.push_back(NULL);
  }
  for (i_se = 0; ok_so_far && i_se < number_of_se; ++i_se) {
    current_se = NULL;
    ok_so_far = ::SetRandomSE(&current_se);
    if (current_se == NULL) ok_so_far = false;
    if (!ok_so_far) continue;
    actual_se.push_back(current_se);
    current_se = NULL;
  }
  ok_so_far = bidimensional::LoadBinaryImage(file_path, &image_d, &image_e);
  if (!ok_so_far || image_d == NULL) return -2;
  if (do_save) {
    ok_so_far = bidimensional::SaveBinaryImage(file_path+".input_d.png", *image_d);
    if (!ok_so_far || image_d == NULL) return -3;
    ok_so_far = bidimensional::SaveBinaryImage(file_path+".input_e.png", *image_e);
    if (!ok_so_far || image_e == NULL) return -3;
  }
  width = image_e->Length(0);
  if (width < 1) return -4;
  height = image_e->Length(1);
  if (height < 1) return -4;
  if (image_info) {
    for (x = 0; ok_so_far && x < width; ++x) {
      ok_so_far = (::p)->set_value(0, x);
      if (!ok_so_far) continue;
      for (y = 0; ok_so_far && y < height; ++y) {
        ok_so_far = (::p)->set_value(1, y);
        if (!ok_so_far) continue;
        ok_so_far = image_e->value(*(::p), &position_value);
        if (!ok_so_far) continue;
        if (position_value) {
          ++input_foreground;
        } else {
          ++input_background;
        }
      }
    }
    if (input_foreground+input_background != width*height)
      ok_so_far = false;
  }
  if (!ok_so_far) return ok_so_far;
  if (be_verbose) {
    printf("Image size: %d x %d\n", width, height);
    if (image_info) {
      printf("Input image foreground pixels: %d"
             "\tInput image background pixels: %d\n",
             input_foreground, input_background);
    }
  }
  if (seed == -1) {
    // After setting SEs, set a random seed for rand().
    srand(time(NULL));
  } else {
    srand(seed);
  }
  // Obtain resulting images using selected algorithms.
  for (delta = 0; ok_so_far && delta < 2; ++delta) {
    const bool true_for_erosion = (delta == 0);
    for (i = delta*total_algorithms/2;
        ok_so_far && i < (delta+1)*total_algorithms/2; ++i) {
      if (!algorithms[i]) continue;
      if (debug || be_verbose) {
        debug_output << "Algorithm: ";
        switch (i%(total_algorithms/2)) {
          case 0:  debug_output << "Naive";
                   break;
          case 1:  debug_output << "Border";
                   break;
          case 2:  debug_output << "Matrix";
                   break;
          default: debug_output << "ERROR";
                   break;
        }
        debug_output << " ";
        if (true_for_erosion) {
          debug_output << "erosion";
        } else {
          debug_output << "dilation";
        }
        debug_output << "\n";
      }
      algorithm_determinate_border_comparison_counter =
          new std::vector< imaging::ImagePositionIndex >();
      algorithm_insert_new_candidate_comparison_counter =
          new std::vector< imaging::ImagePositionIndex >();
      algorithm_insert_new_candidate_memory_access_counter =
          new std::vector< imaging::ImagePositionIndex >();
      algorithm_remove_candidate_comparison_counter =
          new std::vector< imaging::ImagePositionIndex >();
      algorithm_remove_candidate_memory_access_counter =
          new std::vector< imaging::ImagePositionIndex >();
      algorithm_number_of_elements_in_border =
          new std::vector< imaging::ImagePositionIndex >();
      if (algorithm_determinate_border_comparison_counter == NULL
          || algorithm_insert_new_candidate_comparison_counter == NULL
          || algorithm_insert_new_candidate_memory_access_counter == NULL
          || algorithm_remove_candidate_comparison_counter == NULL
          || algorithm_remove_candidate_memory_access_counter == NULL
          || algorithm_number_of_elements_in_border == NULL) {
        ok_so_far = false;
        continue;
      }
      algorithm_determinate_border_comparison_counter->clear();
      algorithm_insert_new_candidate_comparison_counter->clear();
      algorithm_insert_new_candidate_memory_access_counter->clear();
      algorithm_remove_candidate_comparison_counter->clear();
      algorithm_remove_candidate_memory_access_counter->clear();
      algorithm_number_of_elements_in_border->clear();
      if (current_transform != NULL) {
        ok_so_far = false;
        continue;
      }
      switch (i%(total_algorithms/2)) {
        case 0:  if (true_for_erosion) {
                   current_transform = new imaging::binary::morphology::NaiveErosion(debug, debug_output);
                   if (current_transform == NULL) {
                    ok_so_far = false;
                    continue;
                   }
                 } else {
                   current_transform = new imaging::binary::morphology::NaiveDilation(debug, debug_output);
                   if (current_transform == NULL) {
                    ok_so_far = false;
                    continue;
                   }                  
                 }
                 break;
        case 1:  if (true_for_erosion) {
                   current_transform = new imaging::binary::morphology::BorderErosion(debug, debug_output);
                   if (current_transform == NULL) {
                    ok_so_far = false;
                    continue;
                   }
                 } else {
                   current_transform = new imaging::binary::morphology::BorderDilation(debug, debug_output);
                   if (current_transform == NULL) {
                    ok_so_far = false;
                    continue;
                   }                  
                 }
                 break;
        case 2:  if (true_for_erosion) {
                   current_transform = new imaging::binary::morphology::MatrixErosion(debug, debug_output);
                   if (current_transform == NULL) {
                    ok_so_far = false;
                    continue;
                   }
                 } else {
                   current_transform = new imaging::binary::morphology::MatrixDilation(debug, debug_output);
                   if (current_transform == NULL) {
                    ok_so_far = false;
                    continue;
                   }                  
                 }
                 break;
        default:
                 break;
      }
      if (true_for_erosion) {
        image = image_e;
      } else {
        image = image_d;
      }
      ok_so_far = current_transform->Calculate(*image, actual_se,
          &current_output, algorithm_determinate_border_comparison_counter,
          algorithm_insert_new_candidate_comparison_counter,
          algorithm_insert_new_candidate_memory_access_counter,
          algorithm_remove_candidate_comparison_counter,
          algorithm_remove_candidate_memory_access_counter,
          algorithm_number_of_elements_in_border, &start, &end);
      if (current_output == NULL) ok_so_far = false;
      if (!ok_so_far) {
        result |= 1 << 1;
        continue;
      }
      image = NULL;
      unpadded_output = new imaging::grayscale::Image(imaging::Size(), 0);
      ok_so_far = current_output->UnpaddedImage(unpadded_output);
      output.at(i) = unpadded_output;
      delete current_output;
      current_output = NULL;
      unpadded_output = NULL;
      determinate_border_comparison_counter.at(i) =
          algorithm_determinate_border_comparison_counter;
      insert_new_candidate_comparison_counter.at(i) =
          algorithm_insert_new_candidate_comparison_counter;
      insert_new_candidate_memory_access_counter.at(i) =
          algorithm_insert_new_candidate_memory_access_counter;
      remove_candidate_comparison_counter.at(i) =
          algorithm_remove_candidate_comparison_counter;
      remove_candidate_memory_access_counter.at(i) =
          algorithm_remove_candidate_memory_access_counter;
      number_of_elements_in_border.at(i) =
          algorithm_number_of_elements_in_border;
      algorithm_determinate_border_comparison_counter = NULL;
      algorithm_insert_new_candidate_comparison_counter = NULL;
      algorithm_insert_new_candidate_memory_access_counter = NULL;
      algorithm_remove_candidate_comparison_counter = NULL;
      algorithm_remove_candidate_memory_access_counter = NULL;
      algorithm_number_of_elements_in_border = NULL;
      times.at(i) = end-start;
      delete current_transform;
      current_transform = NULL;
    }
  }
  // Iterate to save output images or to display info, if requested.
  iterate_for_save_or_info = (result == 0) && (do_save || image_info);
  for (i = 0;
      ok_so_far && i < total_algorithms && iterate_for_save_or_info;
      ++i) {
    const bool true_for_erosion = (i < total_algorithms/2);
    if (!algorithms[i]) continue;
    if (!ok_so_far) continue;
    if (do_save) {
      std::string suffix(".");
      if (true_for_erosion)
        suffix += "erosion";
      else
        suffix += "dilation";
      suffix += "_";
      switch (i%(total_algorithms/2)) {
        case 0:  suffix += "naive";
                 break;
        case 1:  suffix += "border";
                 break;
        case 2:  suffix += "matrix";
                 break;
        default: break;
      }
      ok_so_far = bidimensional::SaveGrayscaleImage(
          (file_path+suffix+".png").c_str(), *(output.at(i)));
      if (!ok_so_far) continue;
      if (output_text.at(i) != NULL) ok_so_far = false;
      if (!ok_so_far) continue;
      output_text.at(i) = new std::ofstream((file_path+suffix+".txt").c_str(),
          std::ofstream::out);
      if (output_text.at(i) == NULL) ok_so_far = false;
      if (!ok_so_far) continue;
      if ((output_text.at(i))->is_open()) ok_so_far = false;
      ok_so_far = (output.at(i))->Print(*(output_text.at(i)));
    }
    if (!ok_so_far) continue;
    if (image_info) {
      for (y = 0; ok_so_far && y < height; ++y) {
        ok_so_far = (::p)->set_value(1, y);
        if (!ok_so_far) continue;
        for (x = 0; ok_so_far && x < width; ++x) {
          ok_so_far = (::p)->set_value(0, x);
          if (!ok_so_far) continue;
          ok_so_far = (output.at(i))->value(*(::p), &output_position_value);
          if (!ok_so_far) continue;
          if (output_position_value > -1) {
            ++(foreground.at(i));
          } else {
            ++(background.at(i));
          }
        }
      }
      if (foreground.at(i)+background.at(i) != width*height)
        ok_so_far = false;
    }
  }
  // Compair obtained images.
  for (delta = 0; ok_so_far && delta < 2; ++delta) {
    const bool true_for_erosion = (delta == 0);
    const int last_for_comparison = (delta+1)*total_algorithms/2-1;
    for (m = delta*(total_algorithms/2); m < last_for_comparison; ++m) {
      if (!algorithms[m]) continue;
      for (n = m+1; n < last_for_comparison+1; ++n) {
        if (!algorithms[n]) continue;
        if (!(output.at(m))->Equals(*(output.at(n)))) {
          debug_output << "error: ";
          if (true_for_erosion) {
            debug_output << "eroded";
          } else {
            debug_output << "dilated";
          }
          debug_output << " images differ.\n";
          if (be_verbose) {
            debug_output << "algorithm " << m
                << " - foreground pixels: " << foreground.at(m)
                << "\tbackground pixels: " << background.at(m) << "\n";
            debug_output << "algorithm " << n
                << " - foreground pixels: " << foreground.at(n)
                << "\tbackground pixels: " << background.at(n) << "\n";
          }
          result |= 1 << 2;
          ok_so_far = false;
        } else {
          if (be_verbose) {
            printf("verbose message: ");
            if (true_for_erosion) {
              printf("eroded");
            } else {
              printf("dilated");
            }
            printf(" images equal.\n");
          }
        }
      }
    }
  }
  // Printing output images.
  if (debug || be_verbose) {
    for (i = 0; ok_so_far && i < total_algorithms; ++i) {
      if (!algorithms[i]) continue;
      if (output.at(i) == NULL) ok_so_far = false;
      if (!ok_so_far) continue;
      ok_so_far = (output.at(i))->Print(std::cout);
    }
  }
  // Print timing data to stdout.
  if (result == 0) {
    printf("%d;%d;%d;%d;",
        1+2*(::half_se_length), number_of_se, width, height);
    for (i = 0; ok_so_far && i < total_algorithms; ++i) {
      const bool true_for_erosion = (i < total_algorithms/2);
      if (be_verbose) {
        switch (i%(total_algorithms/2)) {
          case 0:  printf("Naive");
                   break;
          case 1:  printf("Border");
                   break;
          case 2:  printf("Matrix");
                   break;
          default: break;
        }
        printf(" ");
        if (true_for_erosion) {
          printf("erosion");
        } else {
          printf("dilation");
        }
        printf(": %.4e us.\n", times.at(i));
      } else {
        if (algorithms[i]) printf("%.4e", times.at(i));
        if (i < total_algorithms-1) printf(";");
      }
    }
    printf("\n");
  }
  // Prepare to export counter values to csv files.
  data_counter_variable_name.push_back(
      "determinate border comparison counter");
  data_counter_variable_name.push_back(
      "insert new candidate comparison counter");
  data_counter_variable_name.push_back(
      "insert new candidate memory access counter");
  data_counter_variable_name.push_back(
      "remove candidate comparison counter");
  data_counter_variable_name.push_back(
      "remove candidate memory access counter");
  data_counter_variable_name.push_back(
      "number of elements in border");
  for (delta = 0; ok_so_far && delta < 2; ++delta) {
    int current_operation_iterations = 0;
    std::string suffix(".");
    const bool true_for_erosion = (delta == 0);
    // Verifies the number of iterations.
    for (i = delta*total_algorithms/2; i < (delta+1)*total_algorithms/2; ++i) {
      if (determinate_border_comparison_counter.at(i) == NULL) continue;
      int current_algorithm_iterations =
          static_cast<int>(determinate_border_comparison_counter.at(i)->size());
      if (current_algorithm_iterations == 0) continue;
      if (current_operation_iterations == 0) {
        current_operation_iterations = current_algorithm_iterations;
        continue;
      }
      if (fabs(current_operation_iterations-current_algorithm_iterations) > 1.)
        debug_output << "warning: number of iterations differs by more than 1.\n";
      if (current_operation_iterations < current_algorithm_iterations)
        current_operation_iterations = current_algorithm_iterations;
    }
    if (current_operation_iterations == 0) continue;
    for (i = delta*total_algorithms/2; i < (delta+1)*total_algorithms/2; ++i) {
      if (!algorithms[i]) continue;
      // Prepare output file.
      suffix = ".";
      if (true_for_erosion)
        suffix += "erosion";
      else
        suffix += "dilation";
      suffix += "_";
      switch (i%(total_algorithms/2)) {
        case 0:  suffix += "naive";
                 break;
        case 1:  suffix += "border";
                 break;
        case 2:  suffix += "matrix";
                 break;
        default: break;
      }
      suffix += ".csv";
      output_counter_data_file =
          fopen((counter_data_prefix+suffix).c_str(), "w");
      if (output_counter_data_file == NULL) continue;
      fprintf(output_counter_data_file,"iteration");
      if (debug || be_verbose) {
        printf("\n\nCONTENT OF FILE: %s\n\n",
          (counter_data_prefix+suffix).c_str());
        printf("iteration");
      }
      for (current_string = data_counter_variable_name.begin();
          current_string != data_counter_variable_name.end();
          ++current_string) {
        fprintf(output_counter_data_file, "; %s", current_string->c_str());
        if (debug || be_verbose) printf("; %s", current_string->c_str());
      }
      fprintf(output_counter_data_file,"\n");
      if (debug || be_verbose) printf("\n");
      algorithm_determinate_border_comparison_counter =
          determinate_border_comparison_counter.at(i);
      algorithm_insert_new_candidate_comparison_counter =
          insert_new_candidate_comparison_counter.at(i);
      algorithm_insert_new_candidate_memory_access_counter =
          insert_new_candidate_memory_access_counter.at(i);
      algorithm_remove_candidate_comparison_counter =
          remove_candidate_comparison_counter.at(i);
      algorithm_remove_candidate_memory_access_counter =
          remove_candidate_memory_access_counter.at(i);
      algorithm_number_of_elements_in_border =
          number_of_elements_in_border.at(i);
      for (iteration = 0; iteration < current_operation_iterations;
          ++iteration) {
        fprintf(output_counter_data_file, "%d;", iteration);
        if (debug || be_verbose) printf("%d;", iteration);
        if ( iteration <
            static_cast<int>(
            algorithm_determinate_border_comparison_counter->size()) ) {
          fprintf(output_counter_data_file, "%d",
              algorithm_determinate_border_comparison_counter->at(iteration));
          if (debug || be_verbose)
            printf("%d",
                algorithm_determinate_border_comparison_counter->at(
                iteration));
        }
        fprintf(output_counter_data_file, ";");
        if (debug || be_verbose) printf(";");
        if ( iteration <
            static_cast<int>(
            algorithm_insert_new_candidate_comparison_counter->size()) ) {
          fprintf(output_counter_data_file, "%d",
              algorithm_insert_new_candidate_comparison_counter->at(
              iteration));
          if (debug || be_verbose)
            printf("%d",
                algorithm_insert_new_candidate_comparison_counter->at(
                iteration));
        }
        fprintf(output_counter_data_file, ";");
        if (debug || be_verbose) printf(";");
        if ( iteration <
            static_cast<int>(
            algorithm_insert_new_candidate_memory_access_counter->size()) ) {
          fprintf(output_counter_data_file, "%d",
              algorithm_insert_new_candidate_memory_access_counter->at(
              iteration));
          if (debug || be_verbose)
            printf("%d",
                algorithm_insert_new_candidate_memory_access_counter->at(
                iteration));
        }
        fprintf(output_counter_data_file, ";");
        if (debug || be_verbose) printf(";");
        if ( iteration <
            static_cast<int>(
            algorithm_remove_candidate_comparison_counter->size()) ) {
          fprintf(output_counter_data_file, "%d",
              algorithm_remove_candidate_comparison_counter->at(iteration));
          if (debug || be_verbose)
            printf("%d",
                algorithm_remove_candidate_comparison_counter->at(iteration));
        }
        fprintf(output_counter_data_file, ";");
        if (debug || be_verbose) printf(";");
        if ( iteration <
            static_cast<int>(
            algorithm_remove_candidate_memory_access_counter->size()) ) {
          fprintf(output_counter_data_file, "%d",
              algorithm_remove_candidate_memory_access_counter->at(iteration));
          if (debug || be_verbose)
            printf("%d",
                algorithm_remove_candidate_memory_access_counter->at(iteration)
                );
        }
        fprintf(output_counter_data_file, ";");
        if (debug || be_verbose) printf(";");
        if ( iteration <
            static_cast<int>(
            algorithm_number_of_elements_in_border->size()) ) {
          fprintf(output_counter_data_file, "%d",
              algorithm_number_of_elements_in_border->at(iteration));
          if (debug || be_verbose)
            printf("%d",
                algorithm_number_of_elements_in_border->at(iteration));
        }
        fprintf(output_counter_data_file, "\n");
        if (debug || be_verbose) printf("\n");
      }
      algorithm_determinate_border_comparison_counter = NULL;
      algorithm_insert_new_candidate_comparison_counter = NULL;
      algorithm_insert_new_candidate_memory_access_counter = NULL;
      algorithm_remove_candidate_comparison_counter = NULL;
      algorithm_remove_candidate_memory_access_counter = NULL;
      algorithm_number_of_elements_in_border = NULL;
      fclose(output_counter_data_file);
      output_counter_data_file = NULL;
      if (debug || be_verbose)
        printf("\nFile %s was generated!\n",
            (counter_data_prefix+suffix).c_str());
    }
  }
  // Cleaning up 'output' data.
  for (i = 0; i < total_algorithms; ++i) {
    if (do_save) {
      if (output_text.at(i) != NULL) {
        (output_text.at(i))->close();
        delete (output_text.at(i));
        output_text.at(i) = NULL;
      }
    }
    if (output.at(i) != NULL) {
      delete output.at(i);
      output.at(i) = NULL;
    }
    // Algorithm complexity analyzer data.
    if (determinate_border_comparison_counter.at(i) != NULL) {
      delete determinate_border_comparison_counter.at(i);
      determinate_border_comparison_counter.at(i) = NULL;
    }
    if (insert_new_candidate_comparison_counter.at(i) != NULL) {
      delete insert_new_candidate_comparison_counter.at(i);
      insert_new_candidate_comparison_counter.at(i) = NULL;
    }
    if (insert_new_candidate_memory_access_counter.at(i) != NULL) {
      delete insert_new_candidate_memory_access_counter.at(i);
      insert_new_candidate_memory_access_counter.at(i) = NULL;
    }
    if (remove_candidate_comparison_counter.at(i) != NULL) {
      delete remove_candidate_comparison_counter.at(i);
      remove_candidate_comparison_counter.at(i) = NULL;
    }
    if (remove_candidate_memory_access_counter.at(i) != NULL) {
      delete remove_candidate_memory_access_counter.at(i);
      remove_candidate_memory_access_counter.at(i) = NULL;
    }
    if (number_of_elements_in_border.at(i) != NULL) {
      delete number_of_elements_in_border.at(i);
      number_of_elements_in_border.at(i) = NULL;
    }
  }
  // Cleaning up 'image' data.
  if (image_e != NULL) {
    delete image_e;
    image_e = NULL;
  }
  if (image_d != NULL) {
    delete image_d;
    image_d = NULL;
  }
  if (image != NULL) {
    delete image;
    image = NULL;
  }
  // Cleaning up 'se' data.
  for (i_se = 0; ok_so_far && i_se < number_of_se; ++i_se) {
    current_se = actual_se[i_se];
    if (current_se != NULL) {
      actual_se[i_se] = NULL;
      delete current_se;
      current_se = NULL;
    }
  }
  // Cleaning up static data.
  if (::bb != NULL) {
    delete ::bb;
    ::bb = NULL;
  }
  if (::p != NULL) {
    delete ::p;
    ::p = NULL;
  }
  return result;
}

int main(int argc, const char* argv[]) {
  char usage_buffer[8192];
  sprintf(usage_buffer,
          "usage: '%s' [-i] [-r] [-s] [-v]"
          " image_file_path counter_file_prefix"
          " se_length number_of_se algorithms seed\n"
          "\tOptional:\n"
          "\t\t-i: image information\n"
          "\t\t-r: randomize SEs\n"
          "\t\t-s: save each image\n"
          "\t\t-v: print human readable messages\n"
          "\tRequired:\n"
          "\t\timage_file_path: path of a valid 2D image file\n"
          "\t\tcounter_file_prefix: path of the prefix which will be used to"
          " store counter data\n"
          "\t\tse_length: an odd number to difine structuring elements size"
          " as 'se_length'x'se_length', from 3 to 9\n"
          "\t\tnumber_of_se: number of structuring elements\n"
          "\t\talgorithms: select algorithms using an integer\n"
          "\t\t\t\tbit 0 : naive erosion\n"
          "\t\t\t\tbit 1 : border erosion\n"
          "\t\t\t\tbit 2 : matrix erosion\n"
          "\t\t\t\tbit 3 : naive dilation\n"
          "\t\t\t\tbit 4 : border dilation\n"
          "\t\t\t\tbit 5 : matrix dilation\n"
          "\t\tseed: select seed for random ordering of arrays, "
          "-1 to use time for seed\n",
          argv[0]);
  const int optional = 4;
  const int required = 6;
  const int total_algorithms = 6;
  if (optional < 0 || required < 1) return 1;
  if (argc < required+1 || argc > required+optional+1) {
    printf("%s", usage_buffer);
    return -1;
  }
  bool algorithms[total_algorithms];
  bool be_random = false;
  bool be_verbose = false;
  bool debug = false;
  bool do_save = false;
  const std::string file_path(argv[argc-required]);
  const std::string counter_data_prefix(argv[argc-required+1]);
  int i = 0;
  bool image_info = false;
  const std::string info("-i");
  int number_of_se = 0;
  const std::string random("-r");
  bool return_value = false;
  const std::string save("-s");
  int se_length = 0;
  long seed = -1;
  int selected_algorithms = -1;
  const std::string verbose("-v");
#ifdef DEBUG
  debug = true;
#endif
  for (i = 0; i < total_algorithms; ++i) algorithms[i] = false;
  if (argc > required+1) {
    if (info.compare(argv[1]) != 0
        && random.compare(argv[1]) != 0
        && save.compare(argv[1]) != 0
        && verbose.compare(argv[1]) != 0) {
      printf("%s", usage_buffer);
      return -1;
    }
    switch (argc-required-1) {
      case 1:         if (info.compare(argv[1]) == 0) {
                        image_info = true;
                        break;
                      }
                      if (random.compare(argv[1]) == 0) {
                        be_random = true;
                        break;
                      }
                      if (save.compare(argv[1]) == 0) {
                        do_save = true;
                        break;
                      }
                      if (verbose.compare(argv[1]) == 0) {
                        be_verbose = true;
                        break;
                      }
      case 2:         if (info.compare(argv[1]) == 0
                          && random.compare(argv[2]) == 0) {
                        be_random = true;
                        image_info = true;
                        break;
                      }
                      if (random.compare(argv[1]) == 0
                          && save.compare(argv[2]) == 0) {
                        be_random = true;
                        do_save = true;
                        break;
                      }
                      if (save.compare(argv[1]) == 0
                          && verbose.compare(argv[2]) == 0) {
                        do_save = true;
                        be_verbose = true;
                        break;
                      }
                      if (info.compare(argv[1]) == 0
                          && save.compare(argv[2]) == 0) {
                        do_save = true;
                        image_info = true;
                        break;
                      }
                      if (info.compare(argv[1]) == 0
                          && verbose.compare(argv[2]) == 0) {
                        image_info = true;
                        be_verbose = true;
                        break;
                      }
                      if (random.compare(argv[1]) == 0
                          && verbose.compare(argv[2]) == 0) {
                        be_random = true;
                        be_verbose = true;
                        break;
                      }
      case 3:         if (info.compare(argv[1]) == 0
                          && random.compare(argv[2]) == 0
                          && save.compare(argv[3]) == 0) {
                        image_info = true;
                        be_random = true;
                        do_save = true;
                        break;
                      }
                      if (random.compare(argv[1]) == 0
                          && save.compare(argv[2]) == 0
                          && verbose.compare(argv[3]) == 0) {
                        be_random = true;
                        do_save = true;
                        be_verbose = true;
                        break;
                      }
                      if (info.compare(argv[1]) == 0
                          && random.compare(argv[2]) == 0
                          && verbose.compare(argv[3]) == 0) {
                        image_info = true;
                        be_random = true;
                        be_verbose = true;
                        break;
                      }
                      if (info.compare(argv[1]) == 0
                          && save.compare(argv[2]) == 0
                          && verbose.compare(argv[3]) == 0) {
                        image_info = true;
                        do_save = true;
                        be_verbose = true;
                        break;
                      }
      case optional:  if (info.compare(argv[1]) == 0
                          && random.compare(argv[2]) == 0
                          && save.compare(argv[3]) == 0
                          && verbose.compare(argv[4]) == 0) {
                        image_info = true;
                        be_random = true;
                        do_save = true;
                        be_verbose = true;
                        break;
                      }
      default: printf("%s", usage_buffer);
               return -1;
    }
  }
  if (be_random) {
    srand(time(NULL));
  } else {
    srand(0);
  }
  se_length = atoi(argv[argc-4]);
  number_of_se = atoi(argv[argc-3]);
  selected_algorithms = atoi(argv[argc-2]);
  seed = atoi(argv[argc-1]);
  if (se_length < 2 || se_length > 9 || number_of_se < 1
      || selected_algorithms < 0
      || selected_algorithms > (1<<total_algorithms)-1
      || seed < -1) {
    printf("%s", usage_buffer);
    return -1;
  }
  if ((se_length % 2) != 1) {
    printf("%s", usage_buffer);
    return -1;
  }
  ::half_se_length = se_length/2;
  if (be_verbose) {
    printf("Image info? ");
    if (image_info) printf("Yes"); else printf("No");
    printf("\n");
    printf("Random SEs? ");
    if (be_random) printf("Yes"); else printf("No");
    printf("\n");
    printf("Save? ");
    if (do_save) printf("Yes"); else printf("No");
    printf("\n");
  }
  for (i = 0; i < total_algorithms; ++i) {
    if ((selected_algorithms & (1<<i)) == (1<<i)) {
      algorithms[i] = true;
    } else {
      algorithms[i] = false;
    }
  }
  return_value = tester(debug, file_path, counter_data_prefix, number_of_se,
        be_verbose, do_save, image_info, seed, total_algorithms, algorithms);
  return return_value;
}
