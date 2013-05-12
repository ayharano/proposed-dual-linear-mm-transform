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

// This file contains the declaration of matrix border erosion and
// dilation transforms.

#ifndef MATRIX_H_
#define MATRIX_H_

#include "img.h"

namespace imaging {


namespace binary {


namespace morphology {


class MatrixDilation : public DilationTransform {
 public:
  MatrixDilation(const bool debug, std::ostream &debug_output)
      : DilationTransform(true, false, debug, debug_output) {}
  virtual ~MatrixDilation() {}
 protected:
  virtual bool clear();
  virtual bool CustomInitialize();
  virtual bool DetectBorder(
      const imaging::SEIndex current_se_index,
      const std::vector<imaging::ImagePositionIndex> &current_se_indexes);
  virtual bool InitialCandidatePositionFound(
      const imaging::binary::Image &image,
      const imaging::ImagePositionIndex &image_position,
      const imaging::Position &value);
  virtual bool InsertNewCandidateFromBorder(
      imaging::grayscale::Image **output_image);
  virtual bool RemoveCandidateNode(
      const imaging::ImagePositionIndex &image_position);
 private:
  MatrixDilation() : DilationTransform(false, true, false, std::cout) {}
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

  DISALLOW_COPY_AND_ASSIGN(MatrixDilation);
}; // imaging:::binary::morphology::MatrixDilation


class MatrixErosion : public ErosionTransform {
 public:
  MatrixErosion(const bool debug, std::ostream &debug_output)
      : ErosionTransform(true, false, debug, debug_output) {}
  virtual ~MatrixErosion() {}
 protected:
  virtual bool clear();
  virtual bool CustomInitialize();
  virtual bool DetectBorder(
      const imaging::SEIndex current_se_index,
      const std::vector<imaging::ImagePositionIndex> &current_se_indexes);
  virtual bool InitialCandidatePositionFound(
      const imaging::binary::Image &image,
      const imaging::ImagePositionIndex &image_position,
      const imaging::Position &value);
  virtual bool InsertNewCandidateFromBorder(
      imaging::grayscale::Image **output_image);
  virtual bool RemoveCandidateNode(
      const imaging::ImagePositionIndex &image_position);
 private:
  MatrixErosion() : ErosionTransform(false, true, false, std::cout) {}
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

  DISALLOW_COPY_AND_ASSIGN(MatrixErosion);
}; // imaging:::binary::morphology::MatrixErosion


} // namespace imaging::binary::morphology


} // namespace imaging::binary


} // namespace imaging

#endif // MATRIX_H_
