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

// This file contains the declaration of naive erosion and
// dilation transforms.

#ifndef NAIVE_H_
#define NAIVE_H_

#include "img.h"

namespace imaging {


namespace binary {


namespace morphology {


class NaiveDilation : public DilationTransform {
 public:
  NaiveDilation(const bool debug, std::ostream &debug_output)
      : DilationTransform(false, true, debug, debug_output) {}
  virtual ~NaiveDilation() {}
 protected:
  virtual bool DetectBorder(
      const imaging::SEIndex current_se_index,
      const std::vector<imaging::ImagePositionIndex> &current_se_indexes);
  virtual bool InitialCandidatePositionFound(
      const imaging::binary::Image &image,
      const imaging::ImagePositionIndex &image_position,
      const imaging::Position &value);
  virtual bool InsertNewCandidateFromBorder(
      imaging::grayscale::Image **output_image);
 private:
  NaiveDilation() : DilationTransform(true, false, false, std::cout) {}
  DISALLOW_COPY_AND_ASSIGN(NaiveDilation);
}; // imaging:::binary::morphology::NaiveDilation


class NaiveErosion : public ErosionTransform {
 public:
  NaiveErosion(const bool debug, std::ostream &debug_output)
      : ErosionTransform(false, true, debug, debug_output) {}
  virtual ~NaiveErosion() {}
 protected:
  virtual bool DetectBorder(
      const imaging::SEIndex current_se_index,
      const std::vector<imaging::ImagePositionIndex> &current_se_indexes);
  virtual bool InitialCandidatePositionFound(
      const imaging::binary::Image &image,
      const imaging::ImagePositionIndex &image_position,
      const imaging::Position &value);
  virtual bool InsertNewCandidateFromBorder(
      imaging::grayscale::Image **output_image);
 private:
  NaiveErosion() : ErosionTransform(true, false, false, std::cout) {}
  DISALLOW_COPY_AND_ASSIGN(NaiveErosion);
}; // imaging:::binary::morphology::NaiveErosion


} // namespace imaging::binary::morphology


} // namespace imaging::binary


} // namespace imaging

#endif // NAIVE_H_

