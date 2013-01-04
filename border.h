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

// This file contains the declaration of erosion and dilation transforms,
// which uses image's border.

#ifndef BORDER_H_
#define BORDER_H_

#include "img.h"

namespace imaging {


namespace binary {


namespace morphology {


class Border : public DualOperation {
 public:
  Border(const bool debug, std::ostream &debug_output)
      : DualOperation(true, true, debug, debug_output) {}
  virtual ~Border() {}
 protected:
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
 private:
  Border() : DualOperation(true, true, false, std::cout) {}
  DISALLOW_COPY_AND_ASSIGN(Border);
}; // imaging:::binary::morphology::Border


} // namespace imaging::binary::morphology


} // namespace imaging::binary


} // namespace imaging

#endif // BORDER_H_
