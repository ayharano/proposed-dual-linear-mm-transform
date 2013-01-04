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

// This file contains the implementation of proposed erosion and dilation.

#include "proposed.h"

bool imaging::binary::morphology::Proposed::clear() {
  link_next_.clear();
  link_next_link_.clear();
  link_previous_.clear();
  candidate_next_link_.clear();
  return DualOperation::clear();
}

bool imaging::binary::morphology::Proposed::Debug() {
  DualOperation::Debug();
  return true;
}

bool imaging::binary::morphology::Proposed::CustomInitialize() {
  imaging::SEIndex i = 0;
  std::vector<imaging::ImagePositionIndex> default_values;
  std::vector<imaging::ImagePositionIndex> next_link_default;
  algorithm_insert_new_candidate_memory_access_counter_->at(se_iteration_)
      += 3;
  // Initialize doubly-linked lists' vectors.
  default_values.push_back(imaging::HEADER);
  next_link_default.push_back(u_cardinality());
  for (i = 0; i < u_cardinality(); ++i) {
    algorithm_insert_new_candidate_memory_access_counter_->at(se_iteration_)
        += 3;
    link_next_link_.push_back(next_link_default);
    link_next_.push_back(default_values);
    link_previous_.push_back(default_values);
  }
  // Initialize next link node's vectors.
  candidate_next_link_.push_back(u_cardinality());
  return true;
}

bool imaging::binary::morphology::Proposed::DetectBorder(
      const bool /*true_for_erosion*/,
      const imaging::SEIndex current_se_index,
      const std::vector<imaging::ImagePositionIndex> &current_se_indexes) {
  imaging::ImagePositionIndex first = 0;
  const imaging::ImagePositionIndex &current_cardinality =
          se_cardinality_.at(current_se_index);
  const std::vector< imaging::ImagePositionIndex > &current_se_vector =
      se_elements_.at(current_se_index);
  imaging::ImagePositionIndex i = 0;
  bool ok_so_far = true;
  for (i = 0; ok_so_far && i < current_cardinality; ++i) {
    // For this iteration, calculate border, which is
    //  - internal morphological gradient in case of erosion
    //  - external morphological gradient in case of dilation
    const imaging::ImagePositionIndex current_element_index =
        current_se_indexes.at(i);
    const imaging::ImagePositionIndex element_index =
        current_se_vector.at(current_element_index);
    const std::vector<imaging::ImagePositionIndex> &current_element_next =
        link_next_.at(element_index);
    first = current_element_next.at(imaging::HEADER);
    while (ok_so_far && first != imaging::HEADER) {
      ok_so_far = RemoveCandidateNode(first);
      if (!ok_so_far) continue;
      border_.at(border_counter_) = first;
      ++border_counter_;
      first = current_element_next.at(imaging::HEADER);
    }
  }
  return ok_so_far;
}

bool imaging::binary::morphology::Proposed::InitialCandidatePositionFound(
    const bool true_for_erosion,
    const imaging::binary::Image &image,
    const imaging::ImagePositionIndex &image_position,
    const imaging::Position &value) {
  imaging::SEIndex i = 0;
  imaging::Position neighbor;
  bool ok_so_far = true;
  bool position_value = true;
  // Include new vector positions related to the new position found.
  algorithm_insert_new_candidate_memory_access_counter_->at(se_iteration_)
      += 1;
  candidate_next_link_.push_back(u_cardinality());
  for (i = 0; i < u_cardinality(); ++i) {
    algorithm_insert_new_candidate_memory_access_counter_->at(se_iteration_)
        += 3;
    (link_next_link_.at(i)).push_back(u_cardinality());
    (link_next_.at(i)).push_back(imaging::HEADER);
    (link_previous_.at(i)).push_back(imaging::HEADER);
  }
  // Verify new link nodes.
  for (i = 0; ok_so_far && i < u_cardinality(); ++i) {
    const imaging::Position &delta = u_elements_.at(i);
    algorithm_insert_new_candidate_comparison_counter_->at(se_iteration_) += 1;
    if (true_for_erosion) {
      ok_so_far = value.Sum(delta, &neighbor);
    } else {
      ok_so_far = value.Subtract(delta, &neighbor);
    }
    if (!ok_so_far) continue;
    if (image.IsPositionValid(neighbor)) {
      ok_so_far = Y_->value(neighbor, &position_value);
      if (!ok_so_far) continue;
      if (position_value == true_for_erosion) continue;
    } else {
      if (!true_for_erosion) continue;
    }
    ok_so_far = LinkingProcedure(image_position, i);
  }
  return ok_so_far;
}

bool imaging::binary::morphology::Proposed::InsertNewCandidateFromBorder(
    const bool true_for_erosion,
    imaging::grayscale::Image **output_image) {
  imaging::ImagePositionIndex i = 0;
  imaging::SEIndex j = 0;
  int marked_value = 0;
  imaging::ImagePositionIndex node = 0;
  int one_less_than_the_minimum_value = true_for_erosion ? 0 : -1;
  bool ok_so_far = true;
  imaging::Position target;
  bool position_value = true;
  for (i = 0; ok_so_far && i < border_counter_; ++i) {
    // Update the sparse matrix by including new Candidate and Link nodes
    // at 'candidates'.
    const imaging::Position &p = candidate_position_.at(border_.at(i));
    for (j = 0; ok_so_far && j < u_cardinality(); ++j) {
      const imaging::Position &delta = u_elements_.at(j);
      algorithm_insert_new_candidate_comparison_counter_->at(se_iteration_)
          += 1;
      if (true_for_erosion) {
        ok_so_far = p.Subtract(delta, &target);
      } else {
        ok_so_far = p.Sum(delta, &target);
      }
      if (!ok_so_far) continue;
      if (!Y_->IsPositionValid(target)) continue;
      position_value = !true_for_erosion;
      ok_so_far = Y_->value(target, &position_value);
      if (!ok_so_far) continue;
      if (position_value != true_for_erosion) continue;
      marked_value = 1;
      ok_so_far = (*output_image)->value(target, &marked_value);
      if (!ok_so_far) continue;
      if (marked_value > one_less_than_the_minimum_value) continue;
      ok_so_far = position(target, &node);
      if (!ok_so_far) continue;
      ok_so_far = LinkingProcedure(node, j);
      if (!ok_so_far) continue;
      node = imaging::HEADER;
    }
  }
  return ok_so_far;
}

bool imaging::binary::morphology::Proposed::RemoveCandidateNode(
    const imaging::ImagePositionIndex &image_position) {
  if (image_position == imaging::HEADER) return false;
  bool ok_so_far = true;
  const imaging::ImagePositionIndex first_se =
      candidate_next_link_.at(image_position);
  ok_so_far = RemoveLinkNode(image_position, first_se);
  if (!ok_so_far) return ok_so_far;
  return DualOperation::RemoveCandidateNode(image_position);
}

inline bool imaging::binary::morphology::Proposed::LinkingProcedure(
    const imaging::ImagePositionIndex &image_position,
    const imaging::SEIndex &link_se_index) {
  const imaging::ImagePositionIndex previous =
      (link_previous_.at(link_se_index)).at(imaging::HEADER);
  const imaging::SEIndex &next_se = candidate_next_link_.at(image_position);
  algorithm_insert_new_candidate_memory_access_counter_->at(se_iteration_)
      += 6;
  // Related to linked node.
  (link_next_.at(link_se_index)).at(image_position) = imaging::HEADER;
  (link_previous_.at(link_se_index)).at(image_position) = previous;
  (link_next_.at(link_se_index)).at(previous) = image_position;
  (link_previous_.at(link_se_index)).at(imaging::HEADER) = image_position;
  // Related to candidate node.
  (link_next_link_.at(link_se_index)).at(image_position) = next_se;
  candidate_next_link_.at(image_position) = link_se_index;
  return EnqueueCandidateNode(image_position);
}

inline bool imaging::binary::morphology::Proposed::RemoveLinkNode(
    const imaging::ImagePositionIndex &image_position,
    const imaging::SEIndex &link_se_index) {
  if (link_se_index == u_cardinality()) return true;
  const imaging::ImagePositionIndex next =
      (link_next_.at(link_se_index)).at(image_position);
  const imaging::ImagePositionIndex previous =
      (link_previous_.at(link_se_index)).at(image_position);
  const imaging::ImagePositionIndex next_link =
      (link_next_link_.at(link_se_index)).at(image_position);
  algorithm_remove_candidate_memory_access_counter_->at(se_iteration_) += 2;
  (link_previous_.at(link_se_index)).at(next) = previous;
  (link_next_.at(link_se_index)).at(previous) = next;
  return RemoveLinkNode(image_position, next_link);
}
