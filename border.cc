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

// This file contains the implementation of erosion and dilation transforms,
// which uses image's border.

#include "border.h"

// imaging::binary::morphology::BorderDilation

bool imaging::binary::morphology::BorderDilation::DetectBorder(
      const imaging::SEIndex current_se_index,
      const std::vector<imaging::ImagePositionIndex> &current_se_indexes) {
  imaging::ImagePositionIndex current = 0;
  const imaging::ImagePositionIndex &current_cardinality =
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
      ok_so_far = current_p.Subtract(delta, &target);
      if (!ok_so_far) continue;
      if (Y_->IsPositionValid(target)) {
        ok_so_far = Y_->value(target, &position_value);
        if (!ok_so_far) continue;
        if (!position_value) continue;
      } else {
        continue;
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

bool imaging::binary::morphology::BorderDilation::InitialCandidatePositionFound(
    const imaging::binary::Image &image,
    const imaging::ImagePositionIndex &image_position,
    const imaging::Position &value) {
  bool candidate_found = false;
  std::vector<imaging::Position>::const_iterator delta;
  imaging::Position neighbor;
  bool ok_so_far = true;
  bool position_value = true;
  for (delta = u_elements_.begin();
      ok_so_far && !candidate_found && delta != u_elements_.end();
      ++delta) {
    algorithm_insert_new_candidate_comparison_counter_->at(se_iteration_) += 1;
    ok_so_far = value.Subtract(*delta, &neighbor);
    if (!ok_so_far) continue;
    if (image.IsPositionValid(neighbor)) {
      ok_so_far = Y_->value(neighbor, &position_value);
      if (!ok_so_far) continue;
      if (!position_value) continue;
    } else {
      continue;
    }
    candidate_found = true;
  }
  // Candidate found.
  if (candidate_found) return EnqueueCandidateNode(image_position);
  return ok_so_far;
}

bool imaging::binary::morphology::BorderDilation::InsertNewCandidateFromBorder(
    imaging::grayscale::Image **output_image) {
  imaging::ImagePositionIndex i = 0;
  imaging::SEIndex j = 0;
  int marked_value = 0;
  imaging::ImagePositionIndex node = 0;
  int one_less_than_the_minimum_value = -1;
  bool ok_so_far = true;
  bool position_value = true;
  imaging::Position target;
  for (i = 0; ok_so_far && i < border_counter_; ++i) {
    const imaging::Position &p = candidate_position_.at(border_.at(i));
    for (j = 0; ok_so_far && j < u_cardinality(); ++j) {
      const imaging::Position &delta = u_elements_.at(j);
      algorithm_insert_new_candidate_comparison_counter_->at(se_iteration_)
          += 1;
      ok_so_far = p.Sum(delta, &target);
      if (!ok_so_far) continue;
      if (!Y_->IsPositionValid(target)) continue;
      position_value = true;
      ok_so_far = Y_->value(target, &position_value);
      if (!ok_so_far) continue;
      if (position_value) continue;
      marked_value = 1;
      ok_so_far = (*output_image)->value(target, &marked_value);
      if (!ok_so_far) continue;
      if (marked_value > one_less_than_the_minimum_value) continue;
      ok_so_far = position(target, &node);
      if (!ok_so_far) continue;
      ok_so_far = EnqueueCandidateNode(node);
      if (!ok_so_far) continue;
      node = imaging::HEADER;
    }
  }
  return ok_so_far;
}

// imaging::binary::morphology::BorderErosion

bool imaging::binary::morphology::BorderErosion::DetectBorder(
      const imaging::SEIndex current_se_index,
      const std::vector<imaging::ImagePositionIndex> &current_se_indexes) {
  imaging::ImagePositionIndex current = 0;
  const imaging::ImagePositionIndex &current_cardinality =
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
      ok_so_far = current_p.Sum(delta, &target);
      if (!ok_so_far) continue;
      if (Y_->IsPositionValid(target)) {
        ok_so_far = Y_->value(target, &position_value);
        if (!ok_so_far) continue;
        if (position_value) continue;
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

bool imaging::binary::morphology::BorderErosion::InitialCandidatePositionFound(
    const imaging::binary::Image &image,
    const imaging::ImagePositionIndex &image_position,
    const imaging::Position &value) {
  bool candidate_found = false;
  std::vector<imaging::Position>::const_iterator delta;
  imaging::Position neighbor;
  bool ok_so_far = true;
  bool position_value = true;
  for (delta = u_elements_.begin();
      ok_so_far && !candidate_found && delta != u_elements_.end();
      ++delta) {
    algorithm_insert_new_candidate_comparison_counter_->at(se_iteration_) += 1;
    ok_so_far = value.Sum(*delta, &neighbor);
    if (!ok_so_far) continue;
    if (image.IsPositionValid(neighbor)) {
      ok_so_far = Y_->value(neighbor, &position_value);
      if (!ok_so_far) continue;
      if (position_value) continue;
    }
    candidate_found = true;
  }
  // Candidate found.
  if (candidate_found) return EnqueueCandidateNode(image_position);
  return ok_so_far;
}

bool imaging::binary::morphology::BorderErosion::InsertNewCandidateFromBorder(
    imaging::grayscale::Image **output_image) {
  imaging::ImagePositionIndex i = 0;
  imaging::SEIndex j = 0;
  int marked_value = 0;
  imaging::ImagePositionIndex node = 0;
  int one_less_than_the_minimum_value = 0;
  bool ok_so_far = true;
  bool position_value = true;
  imaging::Position target;
  for (i = 0; ok_so_far && i < border_counter_; ++i) {
    const imaging::Position &p = candidate_position_.at(border_.at(i));
    for (j = 0; ok_so_far && j < u_cardinality(); ++j) {
      const imaging::Position &delta = u_elements_.at(j);
      algorithm_insert_new_candidate_comparison_counter_->at(se_iteration_)
          += 1;
      ok_so_far = p.Subtract(delta, &target);
      if (!ok_so_far) continue;
      if (!Y_->IsPositionValid(target)) continue;
      position_value = false;
      ok_so_far = Y_->value(target, &position_value);
      if (!ok_so_far) continue;
      if (!position_value) continue;
      marked_value = 1;
      ok_so_far = (*output_image)->value(target, &marked_value);
      if (!ok_so_far) continue;
      if (marked_value > one_less_than_the_minimum_value) continue;
      ok_so_far = position(target, &node);
      if (!ok_so_far) continue;
      ok_so_far = EnqueueCandidateNode(node);
      if (!ok_so_far) continue;
      node = imaging::HEADER;
    }
  }
  return ok_so_far;
}
