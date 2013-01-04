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

// This file contains the implementation of naive erosion and
// dilation transforms.

#include <cstdio>

#include "naive.h"

bool imaging::binary::morphology::Naive::DetectBorder(
      const bool true_for_erosion,
      const imaging::SEIndex current_se_index,
      const std::vector<imaging::ImagePositionIndex> &current_se_indexes) {
  imaging::ImagePositionIndex current = 0;
  const imaging::ImagePositionIndex current_cardinality =
          se_cardinality_.at(current_se_index);
  const std::vector< imaging::ImagePositionIndex > &current_se_vector =
      se_elements_.at(current_se_index);
  imaging::ImagePositionIndex i = 0;
  bool keep_pixel = true;
  imaging::ImagePositionIndex next = 0;
  bool ok_so_far = true;
  bool position_value = true;
  imaging::Position target;
  current = candidate_next_.at(imaging::HEADER);
  while (ok_so_far && current != imaging::HEADER) {
    keep_pixel = true;
    const imaging::Position &current_p = candidate_position_.at(current);
    next = candidate_next_.at(current);
    for (i = 0; ok_so_far && keep_pixel && i < current_cardinality; ++i) {
      // For this iteration, calculate border, which is
      //  - internal morphological gradient in case of erosion
      //  - external morphological gradient in case of dilation
      const imaging::ImagePositionIndex current_element_index =
          current_se_indexes.at(i);
      const imaging::ImagePositionIndex element_index =
          current_se_vector.at(current_element_index);
      const imaging::Position &delta = u_elements_.at(element_index);
      algorithm_determinate_border_comparison_counter_->at(se_iteration_) += 1;
      if (true_for_erosion) {
        ok_so_far = current_p.Sum(delta, &target);
      } else {
        ok_so_far = current_p.Subtract(delta, &target);
      }
      if (!ok_so_far) continue;
      if (Y_->IsPositionValid(target)) {
        ok_so_far = Y_->value(target, &position_value);
        if (!ok_so_far) continue;
        if (position_value == true_for_erosion) continue;
      } else {
        if (!true_for_erosion) continue;
      }
      keep_pixel = false;
    }
    if (!keep_pixel) {
      border_.at(border_counter_) = current;
      ++border_counter_;
    }
    current = next;
  }
  return ok_so_far;
}

bool imaging::binary::morphology::Naive::InitialCandidatePositionFound(
    const bool /*true_for_erosion*/,
    const imaging::binary::Image &/*image*/,
    const imaging::ImagePositionIndex &image_position,
    const imaging::Position &/*value*/) {
  return this->EnqueueCandidateNode(image_position);
}

bool imaging::binary::morphology::Naive::InsertNewCandidateFromBorder(
    const bool /*true_for_erosion*/,
    imaging::grayscale::Image **output_image) {
  if (output_image == NULL) return false;
  if (*output_image == NULL) return false;
  return true;
}

