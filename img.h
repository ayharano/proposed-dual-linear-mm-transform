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

// This file contains the declaration of basic imaging classes.

#ifndef IMG_H_
#define IMG_H_

#include <climits>
#include <iostream>
#include <vector>

#include "disallow_ca.h"

namespace imaging {


typedef uint16_t SEIndex; // Type definition for number of SEs


typedef uint32_t ImagePositionIndex; // Type definition for image elements index

const imaging::ImagePositionIndex HEADER = 0;

class Dimension {
 public:
  ~Dimension() {}
  static char number();
  static char Set(const char number);
 private:
  Dimension();
  static char n_;
  DISALLOW_COPY_AND_ASSIGN(Dimension);
}; // imaging::Dimension


class Position {
 public:
  Position();
  Position(const Position &other);
  ~Position();
  Position& operator= (const Position &other);
  bool CopyFrom(const Position &other);
  bool CopyOppositeOf(const Position &other);
  bool Equals(const Position &other) const;
  bool IsOrigin() const;
  bool SetAsOrigin();
  bool set_value(const char index, const long value);
  bool Subtract(const Position &other, Position *result) const;
  bool Sum(const Position &other, Position *result) const;
  bool value(const char index, long *value) const;
 private:
  bool PlusFactor(const Position &other, int factor, Position *result) const;
  std::vector<long> values_;
}; // imaging::Position


// Forward declaration: imaging::Size
class Size;


class BoundingBox {
 public:
  BoundingBox();
  BoundingBox(const BoundingBox &other);
  BoundingBox(const Position &lower, const Position &upper);
  virtual ~BoundingBox();
  virtual BoundingBox& operator= (const BoundingBox &other);
  imaging::ImagePositionIndex capacity() const;
  virtual bool CopyFrom(const BoundingBox &other);
  bool Equals(const BoundingBox &other) const;
  virtual bool Expand(const BoundingBox &other);
  virtual bool Expand(const Position &position);
  bool Intersection(const BoundingBox &other, bool *empty,
      BoundingBox *result) const;
  bool IsValid(const Position &position) const;
  long Length(const char index) const;
  const imaging::Position& lower() const;
  virtual bool ReflectByOrigin(BoundingBox *result) const;
  bool Union(const BoundingBox &other, BoundingBox *result) const;
  const imaging::Position& upper() const;
  const imaging::Size& size() const;
 protected:
  bool set_lower(const char index, const long value);
  bool set_lower(const imaging::Position &position);
  bool set_upper(const char index, const long value);
  bool set_upper(const imaging::Position &position);
 private:
  bool RecalculateSize() const;
  Position lower_;
  mutable bool modified_;
  mutable Size *size_;
  mutable int size_access_counter_;
  Position upper_;
}; // imaging::BoundingBox


class Size : public BoundingBox {
 public:
  Size();
  Size(Position size);
  Size(const Size &other);
  virtual ~Size() {}
  virtual Size& operator= (const Size &other);
  virtual bool CopyFrom(const BoundingBox &other);
  virtual bool Expand(const BoundingBox &other);
  virtual bool Expand(const Position &position);
  virtual bool ReflectByOrigin(BoundingBox **result) const;
}; // imaging::Size


class PositionIterator {
 public:
  PositionIterator(const BoundingBox &target);
  virtual ~PositionIterator() {}
  bool begin();
  bool IsFinished() const;
  bool iterate();
  const Position& value() const;
 private:
  PositionIterator();
  mutable Position current_;
  const BoundingBox &target_;
}; // imaging::PositionIterator


class NDimensionalMatrixInterface {
 public:
  NDimensionalMatrixInterface(const imaging::Size &size, const bool bit_matrix);
  NDimensionalMatrixInterface(const NDimensionalMatrixInterface &other);
  virtual ~NDimensionalMatrixInterface();
  virtual NDimensionalMatrixInterface& operator= (
      const NDimensionalMatrixInterface &other);
  virtual bool CopyFrom(const NDimensionalMatrixInterface &other);
  virtual bool Equals(const NDimensionalMatrixInterface &other) const;
  virtual long Length(const char index) const;
  virtual const imaging::Size& size() const;
 protected:
  NDimensionalMatrixInterface();
  virtual bool CalculateIndexes(const imaging::Position &position,
      imaging::ImagePositionIndex *block_index,
      imaging::ImagePositionIndex *bit_index) const;
 private:
  bool bit_matrix_;
  imaging::Size size_;
}; // imaging::NDimensionalMatrixInterface


namespace grayscale {


namespace _internal {

template< class T >
class NumericalMatrix : public imaging::NDimensionalMatrixInterface {
 public:
  NumericalMatrix(const imaging::Size &size, const T default_value);
  NumericalMatrix(const NumericalMatrix &other);
  ~NumericalMatrix();
  NumericalMatrix& operator= (const NumericalMatrix &other);
  bool clear();
  bool CopyFrom(const NumericalMatrix &other);
  bool Equals(const NumericalMatrix &other) const;
  bool set_value(const imaging::Position &position, const T value);
  bool value(const imaging::Position &position, T *value) const;
 private:
  NumericalMatrix();
  std::vector<T> array_;
}; // imaging::grayscale::_internal::NumericalMatrix


} // namespace imaging::grayscale::_internal


class Image {
 public:
  Image(const imaging::Size &size, const int default_value);
  Image(const Image &other);
  ~Image() {}
  Image& operator= (const Image &other);
  const imaging::BoundingBox& bounding_box() const;
  bool CopyFrom(const Image &other);
  bool Equals(const Image &other) const;
  bool Maximum(const Image &other, bool *empty, Image *result) const;
  bool Minimum(const Image &other, bool *empty, Image *result) const;
  bool IsPositionValid(const imaging::Position &position) const;
  long Length(const char index) const;
  bool Print(std::ostream &out) const;
  bool set_value(const imaging::Position &position, const int value);
  const imaging::Size& size() const;
  bool value(const imaging::Position &position, int *value) const;
 protected:
  Image();
 private:
  mutable imaging::Position int_matrix_position_;
  imaging::BoundingBox bounding_box_;
  imaging::grayscale::_internal::NumericalMatrix<int> data_;
}; // imaging::grayscale::Image


} // namespace imaging::grayscale

namespace binary {


namespace _internal {

typedef int BLOCK;
const long BLOCK_MAX = INT_MAX;

class BitMatrix : public imaging::NDimensionalMatrixInterface {
 public:
  BitMatrix(const imaging::Size &size, const bool empty);
  BitMatrix(const BitMatrix &other);
  ~BitMatrix();
  BitMatrix& operator= (const BitMatrix &other);
  bool CopyFrom(const BitMatrix &other);
  bool Equals(const BitMatrix &other) const;
  bool InvertValues();
  bool set_value(const imaging::Position &position, const bool value);
  bool value(const imaging::Position &position, bool *value) const;
 private:
  BitMatrix();
  bool SetSize(const imaging::Size &size, const bool empty);
  std::vector<BLOCK> array_;
  long blocks_;
}; // imaging::binary::_internal::BitMatrix


} // namespace imaging::binary::_internal


class StructuringElement {
 public:
  StructuringElement(const imaging::BoundingBox &bounding_box,
      const bool empty);
  StructuringElement(const StructuringElement &other);
  virtual ~StructuringElement() {}
  virtual StructuringElement& operator= (const StructuringElement &other);
  const imaging::BoundingBox& bounding_box() const;
  virtual bool CopyFrom(const StructuringElement &other);
  bool DelimitedComplement(StructuringElement *result) const;
  bool Equals(const StructuringElement &other) const;
  bool Intersection(const StructuringElement &other, bool *empty,
      StructuringElement *result) const;
  bool IsPositionValid(const imaging::Position &position) const;
  long Length(const char index) const;
  virtual bool ReflectByOrigin(StructuringElement *result) const;
  bool set_value(const imaging::Position &position, const bool value);
  bool SetMinus(const StructuringElement &to_be_subtracted,
      StructuringElement *result) const;
  const imaging::Size& size() const;
  virtual bool Union(const StructuringElement &other,
      StructuringElement *result) const;
  bool value(const imaging::Position &position, bool *value) const;
 protected:
  StructuringElement();
 private:
  mutable imaging::Position bit_matrix_position_;
  imaging::BoundingBox bounding_box_;
  imaging::binary::_internal::BitMatrix data_;
}; // imaging::binary::StructuringElement


class Image : public StructuringElement {
 public:
  Image(const imaging::Size &size, const bool empty);
  Image(const Image &other);
  virtual ~Image() {}
  virtual Image& operator= (const Image &other);
  virtual bool CopyFrom(const StructuringElement &other);
  virtual bool ReflectByOrigin(StructuringElement *result) const;
  virtual bool Union(const StructuringElement &other,
      StructuringElement *result) const;
 private:
  Image() {}
}; // imaging::binary::Image


namespace morphology {

class DualOperation {
 public:
  DualOperation(const bool use_candidate_matrix, const bool regular_removal,
      const bool debug, std::ostream &debug_output)
      : algorithm_determinate_border_comparison_counter_(NULL),
        algorithm_insert_new_candidate_comparison_counter_(NULL),
        algorithm_insert_new_candidate_memory_access_counter_(NULL),
        algorithm_remove_candidate_comparison_counter_(NULL),
        algorithm_remove_candidate_memory_access_counter_(NULL),
        algorithm_number_of_elements_in_border_(NULL),
        border_counter_(0), candidate_matrix_(NULL), debug_(debug),
        debug_output_(debug_output), se_iteration_(0), Y_(NULL),
        regular_removal_(regular_removal),
        use_candidate_matrix_(use_candidate_matrix) {}
  virtual ~DualOperation();
  bool DilationTransform(const imaging::binary::Image &image,
      const std::vector<imaging::binary::StructuringElement*> &se,
      imaging::grayscale::Image **output,
      std::vector<imaging::ImagePositionIndex>
          *algorithm_determinate_border_comparison_counter,
      std::vector<imaging::ImagePositionIndex>
          *algorithm_insert_new_candidate_comparison_counter,
      std::vector<imaging::ImagePositionIndex>
          *algorithm_insert_new_candidate_memory_access_counter,
      std::vector<imaging::ImagePositionIndex>
          *algorithm_remove_candidate_comparison_counter,
      std::vector<imaging::ImagePositionIndex>
          *algorithm_remove_candidate_memory_access_counter,
      std::vector<imaging::ImagePositionIndex>
          *algorithm_number_of_elements_in_border,
       double *start, double *end);
  bool ErosionTransform(const imaging::binary::Image &image,
      const std::vector<imaging::binary::StructuringElement*> &se,
      imaging::grayscale::Image **output,
      std::vector<imaging::ImagePositionIndex>
          *algorithm_determinate_border_comparison_counter,
      std::vector<imaging::ImagePositionIndex>
          *algorithm_insert_new_candidate_comparison_counter,
      std::vector<imaging::ImagePositionIndex>
          *algorithm_insert_new_candidate_memory_access_counter,
      std::vector<imaging::ImagePositionIndex>
          *algorithm_remove_candidate_comparison_counter,
      std::vector<imaging::ImagePositionIndex>
          *algorithm_remove_candidate_memory_access_counter,
      std::vector<imaging::ImagePositionIndex>
          *algorithm_number_of_elements_in_border,
      double *start, double *end);
 protected:
  virtual bool clear();
  virtual bool CustomInitialize();
  virtual bool Debug();
  virtual bool DetectBorder(
      const bool true_for_erosion,
      const imaging::SEIndex current_se_index,
      const std::vector<imaging::ImagePositionIndex> &current_se_indexes) = 0;
  bool EnqueueCandidateNode(const imaging::ImagePositionIndex &position);
  virtual bool InitialCandidatePositionFound(
      const bool true_for_erosion,
      const imaging::binary::Image &image,
      const imaging::ImagePositionIndex &image_position,
      const imaging::Position &value) = 0;
  virtual bool InsertNewCandidateFromBorder(const bool true_for_erosion,
      imaging::grayscale::Image **output_image) = 0;
  bool position(const imaging::Position &image_position,
      imaging::ImagePositionIndex *value) const;
  bool position(const imaging::ImagePositionIndex &position_index,
      imaging::Position *value) const;
  virtual bool RemoveCandidateNode(
      const imaging::ImagePositionIndex &image_position);
  imaging::SEIndex u_cardinality() const;

  std::vector<imaging::ImagePositionIndex>
      *algorithm_determinate_border_comparison_counter_;
  std::vector<imaging::ImagePositionIndex>
      *algorithm_insert_new_candidate_comparison_counter_;
  std::vector<imaging::ImagePositionIndex>
      *algorithm_insert_new_candidate_memory_access_counter_;
  std::vector<imaging::ImagePositionIndex>
      *algorithm_remove_candidate_comparison_counter_;
  std::vector<imaging::ImagePositionIndex>
      *algorithm_remove_candidate_memory_access_counter_;
  std::vector<imaging::ImagePositionIndex>
      *algorithm_number_of_elements_in_border_;
  std::vector<imaging::ImagePositionIndex> border_;
  imaging::ImagePositionIndex border_counter_;
  std::vector<bool> candidate_initialized_;
  imaging::grayscale::_internal::NumericalMatrix<imaging::ImagePositionIndex>*
      candidate_matrix_;
  std::vector<imaging::ImagePositionIndex> candidate_next_;
  std::vector<imaging::Position> candidate_position_;
  std::vector<imaging::ImagePositionIndex> candidate_previous_;
  bool debug_;
  std::ostream &debug_output_;
  std::vector<imaging::ImagePositionIndex> se_cardinality_;
  std::vector< std::vector< imaging::ImagePositionIndex > > se_elements_;
  imaging::ImagePositionIndex se_iteration_;
  std::vector<imaging::Position> u_elements_;
  imaging::binary::Image* Y_;
 private:
  DualOperation()
      : algorithm_determinate_border_comparison_counter_(NULL),
        algorithm_insert_new_candidate_comparison_counter_(NULL),
        algorithm_insert_new_candidate_memory_access_counter_(NULL),
        algorithm_remove_candidate_comparison_counter_(NULL),
        algorithm_remove_candidate_memory_access_counter_(NULL),
        algorithm_number_of_elements_in_border_(NULL),
        border_counter_(0), candidate_matrix_(NULL), debug_(false),
        debug_output_(std::cout), se_iteration_(0), Y_(NULL),
        regular_removal_(false), use_candidate_matrix_(false) {}
  bool ActualAlgorithm(
      const bool true_for_erosion,
      imaging::grayscale::Image **output_image);
  bool InitializeCandidateData(
      const bool true_for_erosion,
      const imaging::binary::Image &image);
  bool InitializeCounters();
  bool InitializeSEData(
      const std::vector< std::vector<imaging::Position> > &vectorized_se);
  bool Operation(
      const bool true_for_erosion,
      const imaging::binary::Image &image,
      const std::vector<imaging::binary::StructuringElement*> &se,
      imaging::grayscale::Image **output_image,
      std::vector<imaging::ImagePositionIndex>
          *algorithm_determinate_border_comparison_counter,
      std::vector<imaging::ImagePositionIndex>
          *algorithm_insert_new_candidate_comparison_counter,
      std::vector<imaging::ImagePositionIndex>
          *algorithm_insert_new_candidate_memory_access_counter,
      std::vector<imaging::ImagePositionIndex>
          *algorithm_remove_candidate_comparison_counter,
      std::vector<imaging::ImagePositionIndex>
          *algorithm_remove_candidate_memory_access_counter,
      std::vector<imaging::ImagePositionIndex>
          *algorithm_number_of_elements_in_border,
      double *start, double *end);

  bool regular_removal_;
  bool use_candidate_matrix_;
  DISALLOW_COPY_AND_ASSIGN(DualOperation);
}; // imaging::binary::morphology::DualOperation


} // namespace imaging::binary::morphology


} // namespace imaging::binary


} // namespace imaging

#endif // IMG_H_

