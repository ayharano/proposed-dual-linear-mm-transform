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

// This file contains the declaration of proposed erosion and
// dilation transforms.

#ifndef PROPOSED_H_
#define PROPOSED_H_

#include "img.h"

namespace imaging {


namespace binary {


namespace morphology {


class Proposed : public DualOperation {
 public:
  Proposed(const bool debug, std::ostream &debug_output)
      : DualOperation(true, false, debug, debug_output) {}
  virtual ~Proposed() {}
 protected:
  virtual bool clear();
  virtual bool Debug();
  virtual bool CustomInitialize();
  virtual bool DetectBorder(
      const bool true_for_erosion,
      const imaging::SEIndex current_se_index,
      const std::vector<imaging::ImagePositionIndex> &current_se_indexes);
  virtual bool InitialCandidatePositionFound(
      const bool true_for_erosion,
      const imaging::binary::Image &image,
      const imaging::ImagePositionIndex &image_position,
      const imaging::Position &value);
  virtual bool InsertNewCandidateFromBorder(const bool true_for_erosion,
      imaging::grayscale::Image **output_image);
  virtual bool RemoveCandidateNode(
      const imaging::ImagePositionIndex &image_position);
 private:
  Proposed() : DualOperation(false, true, false, std::cout) {}
  inline bool LinkingProcedure(
      const imaging::ImagePositionIndex &image_position,
      const imaging::SEIndex &link_se_index);
  inline bool RemoveLinkNode(
      const imaging::ImagePositionIndex &image_position,
      const imaging::SEIndex &link_se_index);
  std::vector< std::vector<imaging::ImagePositionIndex> > link_next_;
  std::vector< std::vector<imaging::ImagePositionIndex> > link_next_link_;
  std::vector< std::vector<imaging::ImagePositionIndex> > link_previous_;
  std::vector<imaging::ImagePositionIndex> candidate_next_link_;
  DISALLOW_COPY_AND_ASSIGN(Proposed);
}; // imaging:::binary::morphology::Proposed


} // namespace imaging::binary::morphology


} // namespace imaging::binary


} // namespace imaging

#endif // PROPOSED_H_
