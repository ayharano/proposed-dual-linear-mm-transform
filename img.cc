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

// This file contains the implementation of basic imaging classes.

#include <cassert>
#include <iomanip>

#include <sys/time.h>

#include "img.h"
#include "shuffle-inl.h"

namespace {

// Useful static data.

static imaging::ImagePositionIndex NBITS = 0;
static struct timeval timer;

imaging::ImagePositionIndex calculate_number_of_bits(
    imaging::binary::_internal::BLOCK value) {
  return value == 0 ? 0 : 1 + calculate_number_of_bits(value/2);
}

imaging::ImagePositionIndex number_of_bits() {
  if (NBITS != 0) return NBITS;
  NBITS = calculate_number_of_bits(imaging::binary::_internal::BLOCK_MAX);
  return NBITS;
}

bool InitializeAlgorithmsOutputImage(
    const imaging::binary::Image &image,
    imaging::grayscale::Image **output) {
  if (output == NULL) return false;
  if (*output != NULL) return false;
  bool ok_so_far = true;
  bool position_value = false;
  // Initialize transitional output image.
  *output = new imaging::grayscale::Image(image.size(), -1);
  if (*output == NULL) ok_so_far = false;
  if (!ok_so_far) return ok_so_far;
  imaging::PositionIterator iterator(image.size());
  ok_so_far = iterator.begin();
  if (!ok_so_far) return ok_so_far;
  do {
    const imaging::Position &current = iterator.value();
    ok_so_far = image.value(current, &position_value);
    if (!ok_so_far) continue;
    if (position_value) {
      ok_so_far = (*output)->set_value(current, 0);
      if (!ok_so_far) continue;
    }
  } while (ok_so_far && iterator.iterate());
  if (!iterator.IsFinished()) ok_so_far = false;
  return ok_so_far;
}

bool VectorizeSEs(
    const std::vector<imaging::binary::StructuringElement*> &se,
    std::vector< std::vector<imaging::Position> > *vectorized_se) {
  imaging::ImagePositionIndex i_se = 0;
  bool ok_so_far = true;
  bool position_value = true;
  const imaging::ImagePositionIndex &size = se.size();
  if (vectorized_se == NULL) return false;
  for (i_se = 0; ok_so_far && i_se < size; ++i_se) {
    std::vector<imaging::Position> new_vectorized_se;
    const imaging::binary::StructuringElement* current_se = se.at(i_se);
    imaging::PositionIterator iterator(current_se->bounding_box());
    ok_so_far = iterator.begin();
    if (!ok_so_far) continue;
    do {
      const imaging::Position &current = iterator.value();
      if (current.IsOrigin()) {
        ok_so_far = current_se->value(current, &position_value);
        if (!ok_so_far) continue;
        if (!position_value) {
          return false;
        } else {
          continue;
        }
      }
      ok_so_far = current_se->value(current, &position_value);
      if (!ok_so_far) continue;
      if (!position_value) continue;
      new_vectorized_se.push_back(current);
    } while (ok_so_far && iterator.iterate());
    if (!iterator.IsFinished()) ok_so_far = false;
    if (!ok_so_far) continue;
    vectorized_se->push_back(new_vectorized_se);
  }
  if (size !=
      static_cast<imaging::ImagePositionIndex>(vectorized_se->size())) {
    ok_so_far = false;
  }
  return ok_so_far;
}

} // namespace

// imaging::Dimension

char imaging::Dimension::n_ = 0;

char imaging::Dimension::number() {
  if (imaging::Dimension::n_ < 1) return 0;
  return imaging::Dimension::n_;
}

char imaging::Dimension::Set(const char number) {
  if (number < 1) return 0;
  imaging::Dimension::n_ = number;
  return number;
}

imaging::Dimension::Dimension() {
  ; // empty
}

// imaging::Position

imaging::Position::Position() {
  assert(imaging::Dimension::number() > 0);
  SetAsOrigin();
}

imaging::Position::Position(const imaging::Position &other) {
  assert(imaging::Dimension::number() > 0);
  SetAsOrigin();
  CopyFrom(other);
}

imaging::Position::~Position() {
  ; // empty
}

imaging::Position& imaging::Position::operator= (
    const imaging::Position &other) {
  if (this != &other) CopyFrom(other);
  return (*this);
}

bool imaging::Position::CopyFrom(const imaging::Position &other) {
  if (this == &other) return true;
  values_ = other.values_;
  return true;
}

bool imaging::Position::CopyOppositeOf(const imaging::Position &other) {
  if (this == &other) return true;
  long i = 0;
  for (i = 0; i < imaging::Dimension::number(); ++i)
    values_.at(i) = -1 * other.values_.at(i);
  return true;
}

bool imaging::Position::Equals(const imaging::Position &other) const {
  if (this == &other) return true;
  long i = 0;
  for (i = 0; i < imaging::Dimension::number(); ++i)
    if (values_.at(i) != other.values_.at(i)) return false;
  return true;
}

bool imaging::Position::IsOrigin() const{
  long i = 0;
  for (i = 0; i < imaging::Dimension::number(); ++i)
    if (values_.at(i) != 0) return false;
  return true;
}

bool imaging::Position::SetAsOrigin() {
  int i = 0;
  const int N = imaging::Dimension::number();
  if (static_cast<int>(values_.size()) != N) {
    values_.clear();
    for (i = 0; i < N; ++i) values_.push_back(0);
  }
  for (i = 0; i < N; ++i) values_.at(i) = 0;
  return true;
}

bool imaging::Position::set_value(const char index, const long value) {
  const int Index = index;
  const int N = imaging::Dimension::number();
  if (Index < 0 || Index >= N) return false;
  values_.at(Index) = value;
  return true;
}

bool imaging::Position::Subtract(const imaging::Position &other,
    imaging::Position *result) const {
  return PlusFactor(other, -1, result);
}

bool imaging::Position::Sum(const imaging::Position &other,
    imaging::Position *result) const {
  return PlusFactor(other, 1, result);
}

bool imaging::Position::value(char index, long *value) const {
  const int Index = index;
  const int N = imaging::Dimension::number();
  if (Index < 0 || Index >= N) return false;
  *value = values_.at(Index);
  return true;
}

bool imaging::Position::PlusFactor(const imaging::Position &other, int factor,
    imaging::Position *result) const {
  if (result == NULL) return false;
  if (this == &other) return result->SetAsOrigin();
  long i = 0;
  for (i = 0; i < imaging::Dimension::number(); ++i) {
    (result->values_).at(i) = values_.at(i) + factor*(other.values_).at(i);
  }
  return true;
}


// imaging::BoundingBox

imaging::BoundingBox::BoundingBox() {
  modified_ = true;
  size_ = NULL;
  size_access_counter_ = 0;
  lower_.SetAsOrigin();
  upper_.SetAsOrigin();
}

imaging::BoundingBox::BoundingBox(const imaging::BoundingBox &other) {
  modified_ = true;
  size_ = NULL;
  size_access_counter_ = 0;
  lower_.CopyFrom(other.lower_);
  upper_.CopyFrom(other.upper_);
}

imaging::BoundingBox::BoundingBox(const imaging::Position &lower,
                                  const imaging::Position &upper) {
  char i = 0;
  long maximum = 0;
  long minimum = 0;
  const char n = imaging::Dimension::number();
  bool ok_so_far = true;
  long temp = 0;
  modified_ = true;
  size_ = NULL;
  size_access_counter_ = 0;
  lower_.SetAsOrigin();
  upper_.SetAsOrigin();
  assert(n > 0);
  for (i = 0; ok_so_far && i < n; ++i) {
    ok_so_far = lower.value(i, &minimum);
    if (!ok_so_far) continue;
    ok_so_far = upper.value(i, &maximum);
    if (!ok_so_far) continue;
    if (minimum > maximum) {
      temp = maximum;
      maximum = minimum;
      minimum = temp;
    }
    lower_.set_value(i, minimum);
    upper_.set_value(i, maximum);
  }
}

imaging::BoundingBox::~BoundingBox() {
  if (size_ != NULL) {
    delete size_; // Somewhere may still be using this data.
    size_ = NULL;
  }
}

imaging::BoundingBox& imaging::BoundingBox::operator= (
    const imaging::BoundingBox &other) {
  if (this != &other) {
    if (!Equals(other)) {
      modified_ = true;
      lower_.CopyFrom(other.lower_);
      upper_.CopyFrom(other.upper_);
    }
  }
  return *this;
}

imaging::ImagePositionIndex imaging::BoundingBox::capacity() const {
  char i = 0;
  imaging::ImagePositionIndex result = 1;
  imaging::ImagePositionIndex length = 0;
  const char n = imaging::Dimension::number();
  if (n < 1) return 0;
  for (i = 0; i < n; ++i) {
    length = Length(i);
    assert(length > 0);
    result = result*length;
  }
  assert(result > 0);
  return result;
}

bool imaging::BoundingBox::CopyFrom(const imaging::BoundingBox &other) {
  bool ok_so_far = true;
  if (this == &other) return true;
  if (Equals(other)) return true;
  modified_ = true;
  ok_so_far = lower_.CopyFrom(other.lower_) && upper_.CopyFrom(other.upper_);
  return ok_so_far;
}

bool imaging::BoundingBox::Equals(const imaging::BoundingBox &other) const {
  if (this == &other) return true;
  return lower_.Equals(other.lower_) && upper_.Equals(other.upper_);
}

bool imaging::BoundingBox::Expand(const imaging::BoundingBox &other) {
  char i = 0;
  long maximum_here = 0;
  long maximum_there = 0;
  long minimum_here = 0;
  long minimum_there = 0;
  const char n = imaging::Dimension::number();
  bool ok_so_far = true;
  if (n < 1) return false;
  for (i = 0; ok_so_far && i < n; ++i) {
    ok_so_far = lower_.value(i, &minimum_here);
    if (!ok_so_far) continue;
    ok_so_far = upper_.value(i, &maximum_here);
    if (!ok_so_far) continue;
    ok_so_far = other.lower_.value(i, &minimum_there);
    if (!ok_so_far) continue;
    ok_so_far = other.upper_.value(i, &maximum_there);
    if (!ok_so_far) continue;
    if (minimum_here > minimum_there) {
      ok_so_far = lower_.set_value(i, minimum_there);
      if (!ok_so_far) continue;
      modified_ = true;
    }
    if (maximum_here < maximum_there) {
      ok_so_far = upper_.set_value(i, maximum_there);
      if (!ok_so_far) continue;
      modified_ = true;
    }
  }
  return ok_so_far;
}

bool imaging::BoundingBox::Expand(const imaging::Position &position) {
  char i = 0;
  long maximum_here = 0;
  long maximum_there = 0;
  long minimum_here = 0;
  long minimum_there = 0;
  const char n = imaging::Dimension::number();
  bool ok_so_far = true;
  if (n < 1) return false;
  for (i = 0; ok_so_far && i < n; ++i) {
    ok_so_far = lower_.value(i, &minimum_here);
    if (!ok_so_far) continue;
    ok_so_far = upper_.value(i, &maximum_here);
    if (!ok_so_far) continue;
    ok_so_far = position.value(i, &minimum_there);
    if (!ok_so_far) continue;
    maximum_there = minimum_there;
    if (minimum_here > minimum_there) {
      ok_so_far = lower_.set_value(i, minimum_there);
      if (!ok_so_far) continue;
      modified_ = true;
    }
    if (maximum_here < maximum_there) {
      ok_so_far = upper_.set_value(i, maximum_there);
      if (!ok_so_far) continue;
      modified_ = true;
    }
  }
  return ok_so_far;
}

bool imaging::BoundingBox::Intersection(
    const imaging::BoundingBox &other,
    bool *empty,
    imaging::BoundingBox *result) const {
  if (result == NULL || empty == NULL) return false;
  char i = 0;
  bool intersect = true;
  Position intersection_lower;
  Position intersection_upper;
  long max_minimum = 0;
  long maximum_here = 0;
  long maximum_there = 0;
  long min_maximum = 0;
  long minimum_here = 0;
  long minimum_there = 0;
  const char n = imaging::Dimension::number();
  bool ok_so_far = true;
  if (n < 1) return false;
  for (i = 0; ok_so_far && intersect && i < n; ++i) {
    ok_so_far = lower_.value(i, &minimum_here);
    if (!ok_so_far) continue;
    ok_so_far = upper_.value(i, &maximum_here);
    if (!ok_so_far) continue;
    ok_so_far = other.lower_.value(i, &minimum_there);
    if (!ok_so_far) continue;
    ok_so_far = other.upper_.value(i, &maximum_there);
    if (!ok_so_far) continue;
    if (maximum_there < minimum_here || maximum_here < minimum_there) {
      intersect = false;
      continue;
    }
    max_minimum = minimum_here;
    if (max_minimum < minimum_there) max_minimum = minimum_there;
    min_maximum = maximum_here;
    if (min_maximum > maximum_there) min_maximum = maximum_there;
    ok_so_far = intersection_lower.set_value(i, max_minimum);
    if (!ok_so_far) continue;
    ok_so_far = intersection_upper.set_value(i, min_maximum);
  }
  if (!ok_so_far) return ok_so_far;
  if (!intersect) {
    *empty = true;
    return true;
  }
  ok_so_far = result->set_lower(intersection_lower);
  if (!ok_so_far) return ok_so_far;
  ok_so_far = result->set_upper(intersection_upper);
  if (!ok_so_far) return ok_so_far;
  *empty = false;
  return true;
}

bool imaging::BoundingBox::IsValid(const imaging::Position &position) const {
  char i = 0;
  long maximum = 0;
  long minimum = 0;
  const char n = imaging::Dimension::number();
  bool ok_so_far = true;
  long value = 0;
  if (n < 1) return false;
  for (i = 0; ok_so_far && i < n; ++i) {
    ok_so_far = lower_.value(i, &minimum);
    if (!ok_so_far) continue;
    ok_so_far = upper_.value(i, &maximum);
    if (!ok_so_far) continue;
    ok_so_far = position.value(i, &value);
    if (!ok_so_far) continue;
    if (value < minimum || value > maximum) return false;
  }
  return ok_so_far;
}

long imaging::BoundingBox::Length(const char index) const {
  long maximum = 0;
  long minimum = 0;
  bool ok_so_far = false;
  ok_so_far = lower_.value(index, &minimum);
  if (!ok_so_far) return -1;
  ok_so_far = upper_.value(index, &maximum);
  if (!ok_so_far) return -1;
  if (maximum < minimum) return -1;
  return maximum-minimum+1;
}

const imaging::Position& imaging::BoundingBox::lower() const {
  return lower_;
}

bool imaging::BoundingBox::ReflectByOrigin(imaging::BoundingBox *result)
    const {
  if (result == NULL) return false;
  char i = 0;
  long maximum = 0;
  long minimum = 0;
  const char n = imaging::Dimension::number();
  imaging::Position new_lower;
  imaging::Position new_upper;
  bool ok_so_far = true;
  if (n < 1) return false;
  for (i = 0; ok_so_far && i < n; ++i) {
    ok_so_far = lower_.value(i, &minimum);
    if (!ok_so_far) continue;
    ok_so_far = upper_.value(i, &maximum);
    if (!ok_so_far) continue;
    ok_so_far = new_lower.set_value(i, -maximum);
    if (!ok_so_far) continue;
    ok_so_far = new_upper.set_value(i, -minimum);
    if (!ok_so_far) continue;
  }
  if (!ok_so_far) return ok_so_far;
  ok_so_far = result->set_lower(new_lower);
  if (!ok_so_far) return ok_so_far;
  ok_so_far = result->set_upper(new_upper);
  return ok_so_far;
}

bool imaging::BoundingBox::Union(const imaging::BoundingBox &other,
                                 imaging::BoundingBox *result) const {
  if (result == NULL) return false;
  char i = 0;
  long maximum = 0;
  long maximum_here = 0;
  long maximum_there = 0;
  long minimum = 0;
  long minimum_here = 0;
  long minimum_there = 0;
  char n = 0;
  bool ok_so_far = true;
  Position union_lower;
  Position union_upper;
  if (n < 1) return false;
  for (i = 0; ok_so_far && i < n; ++i) {
    ok_so_far = lower_.value(i, &minimum_here);
    if (!ok_so_far) continue;
    ok_so_far = upper_.value(i, &maximum_here);
    if (!ok_so_far) continue;
    ok_so_far = other.lower_.value(i, &minimum_there);
    if (!ok_so_far) continue;
    ok_so_far = other.upper_.value(i, &maximum_there);
    if (!ok_so_far) continue;
    minimum = minimum_here;
    if (minimum > minimum_there) minimum = minimum_there;
    maximum = maximum_here;
    if (maximum < maximum_there) maximum = maximum_there;
    ok_so_far = union_lower.set_value(i, minimum);
    if (!ok_so_far) continue;
    ok_so_far = union_upper.set_value(i, maximum);
  }
  if (!ok_so_far) return ok_so_far;
  ok_so_far = result->set_lower(union_lower);
  if (!ok_so_far) return ok_so_far;
  ok_so_far = result->set_upper(union_upper);
  return ok_so_far;
}

const imaging::Position& imaging::BoundingBox::upper() const {
  return upper_;
}

const imaging::Size& imaging::BoundingBox::size() const {
  RecalculateSize();
  assert(size_ != NULL);
  ++size_access_counter_;
  return const_cast<const imaging::Size&>(*size_);
}

bool imaging::BoundingBox::set_lower(const char index, const long value) {
  long maximum = 0;
  long minimum = 0;
  const char n = imaging::Dimension::number();
  bool ok_so_far = true;
  if (index < 0 || index >= n) return false;
  ok_so_far = lower_.value(index, &minimum);
  if (minimum == value) return true;
  ok_so_far = upper_.value(index, &maximum);
  if (!ok_so_far) return ok_so_far;
  if (value > maximum) ok_so_far = upper_.set_value(index, value);
  if (!ok_so_far) return ok_so_far;
  ok_so_far = lower_.set_value(index, value);
  if (!ok_so_far) return ok_so_far;
  modified_ = true;
  return ok_so_far;
}

bool imaging::BoundingBox::set_lower(const imaging::Position &position) {
  if (lower_.Equals(position)) return true;
  long coordinate_value = 0;
  char i = 0;
  long maximum = 0;
  const char n = imaging::Dimension::number();
  bool ok_so_far = true;
  if (n < 1) return false;
  for (i = 0; ok_so_far && i < n; ++i) {
    ok_so_far = upper_.value(i, &maximum);
    if (!ok_so_far) continue;
    ok_so_far = position.value(i, &(coordinate_value));
    if (!ok_so_far) continue;
    if (coordinate_value > maximum) {
      ok_so_far = upper_.set_value(i, coordinate_value);
      if (!ok_so_far) continue;
    }
    ok_so_far = lower_.set_value(i, coordinate_value);
  }
  if (!ok_so_far) return ok_so_far;
  modified_ = true;
  return ok_so_far;
}

bool imaging::BoundingBox::set_upper(const char index, const long value) {
  long maximum = 0;
  long minimum = 0;
  const char n = imaging::Dimension::number();
  bool ok_so_far = true;
  if (index < 0 || index >= n) return false;
  ok_so_far = upper_.value(index, &maximum);
  if (maximum == value) return true;
  ok_so_far = lower_.value(index, &minimum);
  if (!ok_so_far) return ok_so_far;
  if (value < minimum) {
    ok_so_far = lower_.set_value(index, value);
    if (!ok_so_far) return ok_so_far;
  }
  ok_so_far = upper_.set_value(index, value);
  if (!ok_so_far) return ok_so_far;
  modified_ = true;
  return ok_so_far;
}

bool imaging::BoundingBox::set_upper(const imaging::Position &position) {
  if (upper_.Equals(position)) return true;
  long coordinate_value = 0;
  char i = 0;
  long minimum = 0;
  const char n = imaging::Dimension::number();
  bool ok_so_far = true;
  if (n < 1) return false;
  for (i = 0; ok_so_far && i < n; ++i) {
    ok_so_far = lower_.value(i, &minimum);
    if (!ok_so_far) continue;
    ok_so_far = position.value(i, &(coordinate_value));
    if (!ok_so_far) continue;
    if (coordinate_value < minimum) {
      ok_so_far = lower_.set_value(i, coordinate_value);
      if (!ok_so_far) continue;
    }
    ok_so_far = upper_.set_value(i, coordinate_value);
  }
  if (!ok_so_far) return ok_so_far;
  modified_ = true;
  return ok_so_far;
}

bool imaging::BoundingBox::RecalculateSize() const {
  if (!modified_) return true;
  if (size_access_counter_ != 0) return false;
  if (size_ != NULL) {
    if (size_access_counter_ == 0) {
      delete size_; // Somewhere may still be using this data.
    }
    size_ = NULL;
  }
  imaging::Position external_upper;
  char i = 0;
  long length = 0;
  const char n = imaging::Dimension::number();
  bool ok_so_far = true;
  if (n < 1) return false;
  for (i = 0; ok_so_far && i < n; ++i) {
    length = Length(i);
    if (length <= 0) {
      ok_so_far = false;
      continue;
    }
    ok_so_far = external_upper.set_value(i, length);
  }
  size_ = new imaging::Size(external_upper);
  modified_ = false;
  size_access_counter_ = 0;
  return size_ != NULL;
}

// imaging::Size

imaging::Size::Size() {
  ; // empty
}

imaging::Size::Size(imaging::Position size) {
  char i = 0;
  const char n = imaging::Dimension::number();
  imaging::Position one;
  bool ok_so_far = true;
  imaging::Position upper;
  long value = 0;
  assert(n > 0);
  for (i = 0; ok_so_far && i < n; ++i) {
    ok_so_far = one.set_value(i, 1);
    ok_so_far = size.value(i, &value);
    if (!ok_so_far) continue;
    if (value < 1) ok_so_far = false;
  }
  if (ok_so_far) {
    ok_so_far = size.Subtract(one, &upper);
    if (!ok_so_far) return;
    this->Expand(BoundingBox(Position(), upper));
  }
}

imaging::Size::Size(const imaging::Size &other) : imaging::BoundingBox(other) {
  ; // empty
}

imaging::Size& imaging::Size::operator= (const imaging::Size &other) {
  if (this != &other) imaging::BoundingBox::operator=(other);
  return *this;
}

bool imaging::Size::CopyFrom(const imaging::BoundingBox &other) {
  if (this == &other) return true;
  char i = 0;
  imaging::Position lower_position;
  long maximum_there = 0;
  long minimum_there = 0;
  const char n = imaging::Dimension::number();
  bool ok_so_far = true;
  const imaging::Position &other_lower = other.lower();
  const imaging::Position &other_upper = other.upper();
  imaging::Position upper_position;
  if (n < 1) return false;
  for (i = 0; ok_so_far && i < n; ++i) {
    ok_so_far = other_lower.value(i, &minimum_there);
    if (!ok_so_far) continue;
    ok_so_far = other_upper.value(i, &maximum_there);
    if (!ok_so_far) continue;
    if (minimum_there < 0 ) {
      ok_so_far = lower_position.set_value(i, 0);
    } else {
      ok_so_far = lower_position.set_value(i, minimum_there);
    }
    if (!ok_so_far) continue;
    if (maximum_there < 0 ) {
      ok_so_far = upper_position.set_value(i, 0);
    } else {
      ok_so_far = upper_position.set_value(i, maximum_there);
    }
    if (!ok_so_far) continue;
  }
  if (!ok_so_far) return ok_so_far;
  return set_lower(lower_position) && set_upper(upper_position);
}

bool imaging::Size::Expand(const imaging::BoundingBox &other) {
  char i = 0;
  long maximum_here = 0;
  long maximum_there = 0;
  long minimum_here = 0;
  long minimum_there = 0;
  const char n = imaging::Dimension::number();
  bool ok_so_far = true;
  const imaging::Position &other_lower = other.lower();
  const imaging::Position &other_upper = other.upper();
  const imaging::Position &this_lower = this->lower();
  const imaging::Position &this_upper = this->upper();
  if (n < 1) return false;
  for (i = 0; ok_so_far && i < n; ++i) {
    ok_so_far = this_lower.value(i, &minimum_here);
    if (!ok_so_far) continue;
    ok_so_far = this_upper.value(i, &maximum_here);
    if (!ok_so_far) continue;
    ok_so_far = other_lower.value(i, &minimum_there);
    if (!ok_so_far) continue;
    ok_so_far = other_upper.value(i, &maximum_there);
    if (!ok_so_far) continue;
    if (minimum_here > minimum_there) {
      if (minimum_there < 0) {
        ok_so_far = set_lower(i, 0);
      } else {
        ok_so_far = set_lower(i, minimum_there);
      }
    }
    if (!ok_so_far) continue;
    if (maximum_here < maximum_there) {
      if (maximum_there < 0) {
        ok_so_far = set_upper(i, 0);
      } else {
        ok_so_far = set_upper(i, maximum_there);
      }
    }
    if (!ok_so_far) continue;
  }
  return ok_so_far;
}

bool imaging::Size::Expand(const imaging::Position &position) {
  char i = 0;
  long maximum_here = 0;
  long maximum_there = 0;
  long minimum_here = 0;
  long minimum_there = 0;
  const char n = imaging::Dimension::number();
  bool ok_so_far = true;
  const imaging::Position &this_lower = this->lower();
  const imaging::Position &this_upper = this->upper();
  if (n < 1) return false;
  for (i = 0; ok_so_far && i < n; ++i) {
    ok_so_far = this_lower.value(i, &minimum_here);
    if (!ok_so_far) continue;
    ok_so_far = this_upper.value(i, &maximum_here);
    if (!ok_so_far) continue;
    ok_so_far = position.value(i, &minimum_there);
    if (!ok_so_far) continue;
    maximum_there = minimum_there;
    if (minimum_here > minimum_there) {
      if (minimum_there < 0) {
        ok_so_far = set_lower(i, 0);
      } else {
        ok_so_far = set_lower(i, minimum_there);
      }
    }
    if (!ok_so_far) continue;
    if (maximum_here < maximum_there) {
      if (maximum_there < 0) {
        ok_so_far = set_upper(i, 0);
      } else {
        ok_so_far = set_upper(i, maximum_there);
      }
    }
    if (!ok_so_far) continue;
  }
  return ok_so_far;
}

bool imaging::Size::ReflectByOrigin(imaging::BoundingBox **result)
    const {
  if (result == NULL) return false;
  if (*result != NULL) return false;
  *result = new imaging::Size(*this);
  if (*result == NULL) return false;
  return true;
}

// imaging::PositionIterator

imaging::PositionIterator::PositionIterator(const imaging::BoundingBox &target)
    : target_(target) {
  ; // empty
}

bool imaging::PositionIterator::begin() {
  return current_.CopyFrom(target_.lower());
}

bool imaging::PositionIterator::IsFinished() const {
  return current_.Equals(target_.upper());
}

bool imaging::PositionIterator::iterate() {
  long coordinate_value = 0;
  char i = 0;
  const imaging::Position &lower = target_.lower();
  long maximum = 0;
  long minimum = 0;
  bool next_calculated = false;
  bool ok_so_far = true;
  const imaging::Position &upper = target_.upper();
  // Calculate next position
  if (current_.Equals(upper)) return false;
  while (ok_so_far && !next_calculated) {
    ok_so_far = current_.value(i, &(coordinate_value));
    if (!ok_so_far) continue;
    ok_so_far = upper.value(i, &maximum);
    if (!ok_so_far) continue;
    ++coordinate_value;
    if (coordinate_value > maximum) {
      ok_so_far = lower.value(i, &minimum);
      if (!ok_so_far) continue;
      ok_so_far = current_.set_value(i, minimum);
      ++i;
    } else {
      ok_so_far = current_.set_value(i, coordinate_value);
      next_calculated = true;
    }
  }
  return ok_so_far;
}

const imaging::Position& imaging::PositionIterator::value() const {
  return current_;
}

imaging::PositionIterator::PositionIterator() : target_(BoundingBox()) {
  ; // empty
}

// imaging::NDimensionalMatrixInterface

imaging::NDimensionalMatrixInterface::NDimensionalMatrixInterface(
    const imaging::Size &size, const bool bit_matrix)
    : bit_matrix_(bit_matrix), size_(size) {
  ; // empty
}

imaging::NDimensionalMatrixInterface::NDimensionalMatrixInterface(
    const imaging::NDimensionalMatrixInterface &other)
    : bit_matrix_(other.bit_matrix_), size_(other.size_) {
  ; // empty
}

imaging::NDimensionalMatrixInterface::~NDimensionalMatrixInterface() {
  ; // empty
}

imaging::NDimensionalMatrixInterface&
    imaging::NDimensionalMatrixInterface::operator= (
    const imaging::NDimensionalMatrixInterface &other) {
  if (this != &other) size_.CopyFrom(other.size_);
  return *this;
}

bool imaging::NDimensionalMatrixInterface::CopyFrom(
    const imaging::NDimensionalMatrixInterface &other) {
  if (this == &other) return true;
  bool ok_so_far = true;
  ok_so_far = size_.CopyFrom(other.size_);
  return ok_so_far;
}

bool imaging::NDimensionalMatrixInterface::Equals(
    const imaging::NDimensionalMatrixInterface &other) const {
  if (this == &other) return true;
  if (!size_.Equals(other.size_)) return false;
  return true;
}

long imaging::NDimensionalMatrixInterface::Length(const char index) const {
  return size_.Length(index);
}

const imaging::Size& imaging::NDimensionalMatrixInterface::size() const {
  return size_;
}

imaging::NDimensionalMatrixInterface::NDimensionalMatrixInterface() {
  ; // empty
}

bool imaging::NDimensionalMatrixInterface::CalculateIndexes(
    const imaging::Position &position,
    imaging::ImagePositionIndex *block_index,
    imaging::ImagePositionIndex *bit_index) const {
  if (block_index == NULL) return false;
  if (bit_matrix_ && bit_index == NULL) return false;
  if (!size_.IsValid(position)) return false;
  imaging::ImagePositionIndex absolute_index = 0;
  imaging::ImagePositionIndex accumulated_delta = 1;
  long coordinate_value = 0;
  char i = 0;
  long length = 0;
  const char n = imaging::Dimension::number();
  bool ok_so_far = true;
  if (n < 1) return false;
  for (i = 0; ok_so_far && i < n; ++i) {
    ok_so_far = position.value(i, &coordinate_value);
    if (coordinate_value < 0) ok_so_far = false;
    if (!ok_so_far) continue;
    length = size_.Length(i);
    if (length < 1) ok_so_far = false;
    if (coordinate_value >= length) ok_so_far = false;
    if (!ok_so_far) continue;
    if (i == 0) {
      absolute_index = coordinate_value;
    } else {
      absolute_index += coordinate_value*accumulated_delta;
    }
    accumulated_delta = accumulated_delta*length;
  }
  if (bit_matrix_) {
    *block_index = absolute_index/(::number_of_bits());
    *bit_index = absolute_index%(::number_of_bits());
  } else {
    *block_index = absolute_index;
  }
  return true;
}

// imaging::grayscale::_internal::NumericalMatrix

template< class T >
imaging::grayscale::_internal::NumericalMatrix<T>::NumericalMatrix(
    const imaging::Size &size,
    const T default_value)
    : imaging::NDimensionalMatrixInterface(size, false) {
  array_.resize(size.capacity(), default_value);
}

template< class T >
imaging::grayscale::_internal::NumericalMatrix<T>::NumericalMatrix(
    const imaging::grayscale::_internal::NumericalMatrix<T> &other)
    : imaging::NDimensionalMatrixInterface(other) {
  array_ = other.array_;
}

template< class T >
imaging::grayscale::_internal::NumericalMatrix<T>::~NumericalMatrix() {
  clear();
}

template< class T >
imaging::grayscale::_internal::NumericalMatrix<T>&
    imaging::grayscale::_internal::NumericalMatrix<T>::operator= (
    const imaging::grayscale::_internal::NumericalMatrix<T> &other) {
  if (this != &other) {
    imaging::NDimensionalMatrixInterface::operator=(other);
    array_ = other.array_;
  }
  return *this;
}

template< class T >
bool imaging::grayscale::_internal::NumericalMatrix<T>::clear() {
  array_.clear();
  return true;
}

template< class T >
bool imaging::grayscale::_internal::NumericalMatrix<T>::CopyFrom(
    const imaging::grayscale::_internal::NumericalMatrix<T> &other) {
  if (this == &other) return true;
  if (!imaging::NDimensionalMatrixInterface::CopyFrom(other)) return false;
  array_ = other.array_;
  return true;
}

template< class T >
bool imaging::grayscale::_internal::NumericalMatrix<T>::Equals(
    const imaging::grayscale::_internal::NumericalMatrix<T> &other) const {
  if (this == &other) return true;
  if (!imaging::NDimensionalMatrixInterface::Equals(other)) return false;
  return (array_ == other.array_);
}

template< class T >
bool imaging::grayscale::_internal::NumericalMatrix<T>::set_value(
    const imaging::Position &position, const T value) {
  imaging::ImagePositionIndex index = 0;
  bool ok_so_far = true;
  ok_so_far = CalculateIndexes(position, &index, NULL);
  if (!ok_so_far) return ok_so_far;
  array_.at(index) = value;
  return ok_so_far;
}

template< class T >
bool imaging::grayscale::_internal::NumericalMatrix<T>::value(
    const imaging::Position &position, T *value) const {
  if (value == NULL) return false;
  imaging::ImagePositionIndex index = 0;
  bool ok_so_far = true;
  ok_so_far = CalculateIndexes(position, &index, NULL);
  if (!ok_so_far) return ok_so_far;
  *value = array_.at(index);
  return ok_so_far;
}

template< class T >
imaging::grayscale::_internal::NumericalMatrix<T>::NumericalMatrix()
    : array_(NULL) {
  ; // empty
}

// imaging::grayscale::Image

imaging::grayscale::Image::Image(
    const imaging::Size &size, const int default_value)
    : data_(size, default_value) {
  ; // empty
}

imaging::grayscale::Image::Image(const imaging::grayscale::Image &other)
    : data_(other.data_) {
  ; // empty
}

imaging::grayscale::Image& imaging::grayscale::Image::operator= (
    const imaging::grayscale::Image &other) {
  if (this != &other) data_.CopyFrom(other.data_);
  return *this;
}

const imaging::BoundingBox& imaging::grayscale::Image::bounding_box() const {
  return data_.size();
}

bool imaging::grayscale::Image::CopyFrom(
    const imaging::grayscale::Image &other) {
  if (this == &other) return true;
  return data_.CopyFrom(other.data_);
}

bool imaging::grayscale::Image::Equals(const imaging::grayscale::Image &other)
    const {
  if (this == &other) return true;
  return data_.Equals(other.data_);
}

bool imaging::grayscale::Image::Maximum(
    const Image &other,
    bool *empty,
    imaging::grayscale::Image *result) const {
  if (result == NULL || empty == NULL) return false;
  imaging::Size intersection_size;
  bool ok_so_far = true;
  int that_position_value = 0;
  int this_position_value = 0;
  ok_so_far = data_.size().Intersection(other.data_.size(), empty,
                                         &intersection_size);
  if (!ok_so_far) return ok_so_far;
  if (*empty) return true;
  imaging::grayscale::Image output(intersection_size, true);
  imaging::PositionIterator iterator(intersection_size);
  if (ok_so_far) ok_so_far = iterator.begin();
  if (!ok_so_far) return ok_so_far;
  do {
    ok_so_far = this->value(iterator.value(), &this_position_value);
    if (!ok_so_far) continue;
    ok_so_far = other.value(iterator.value(), &that_position_value);
    if (!ok_so_far) continue;
    if (this_position_value > that_position_value)
      ok_so_far = output.set_value(iterator.value(), this_position_value);
    else
      ok_so_far = output.set_value(iterator.value(), that_position_value);
    if (!ok_so_far) continue;
  } while (ok_so_far && iterator.iterate());
  if (!iterator.IsFinished()) ok_so_far = false;
  if (!ok_so_far) return ok_so_far;
  ok_so_far = result->CopyFrom(output);
  return ok_so_far;
}

bool imaging::grayscale::Image::Minimum(
    const Image &other,
    bool *empty,
    imaging::grayscale::Image *result) const {
  if (result == NULL || empty == NULL) return false;
  imaging::Size intersection_size;
  bool ok_so_far = true;
  int that_position_value = 0;
  int this_position_value = 0;
  ok_so_far = data_.size().Intersection(other.data_.size(), empty,
                                         &intersection_size);
  if (!ok_so_far) return ok_so_far;
  if (*empty) return true;
  imaging::grayscale::Image output(intersection_size, true);
  imaging::PositionIterator iterator(intersection_size);
  if (ok_so_far) ok_so_far = iterator.begin();
  if (!ok_so_far) return ok_so_far;
  do {
    ok_so_far = this->value(iterator.value(), &this_position_value);
    if (!ok_so_far) continue;
    ok_so_far = other.value(iterator.value(), &that_position_value);
    if (!ok_so_far) continue;
    if (this_position_value > that_position_value)
      ok_so_far = output.set_value(iterator.value(), that_position_value);
    else
      ok_so_far = output.set_value(iterator.value(), this_position_value);
    if (!ok_so_far) continue;
  } while (ok_so_far && iterator.iterate());
  if (!iterator.IsFinished()) ok_so_far = false;
  if (!ok_so_far) return ok_so_far;
  ok_so_far = result->CopyFrom(output);
  return ok_so_far;
}

bool imaging::grayscale::Image::IsPositionValid(
    const imaging::Position &position) const {
  return data_.size().IsValid(position);
}

long imaging::grayscale::Image::Length(const char index) const {
  return data_.size().Length(index);
}

bool imaging::grayscale::Image::Print(std::ostream &out) const {
  bool after_last = false;
  std::vector<long> c;
  bool carry = true;
  std::vector<long> delta;
  char i = 0;
  bool is_last = false;
  std::vector<long> last;
  const int N = imaging::Dimension::number();
  imaging::Position next;
  bool ok_so_far = true;
  imaging::Position p;
  bool print_newline = true;
  long value = 0;
  int value_p = 0;
  c.clear();
  c.resize(N, 0);
  delta.clear();
  delta.resize(N, 0);
  last.clear();
  last.resize(N, 0);
  out << "\nSize: ";
  for (i = 0; i < N; ++i) {
    value = this->Length(i);
    out << value;
    if (i != N-1) {
      out << "x";
    } else {
      out << "\n";
    }
    last.at(i) = value-1;
  }
  is_last = (c == last);
  while (!after_last && ok_so_far) {
    ok_so_far = this->value(p, &value_p);
    if (!ok_so_far) continue;
    out << std::setfill(' ') << std::setw (11) << value_p;
    carry = true;
    for (i = N-1; i >= 0 && carry && !is_last; --i) {
      if (c.at(i) + 1 == this->Length(i)) {
        delta.at(i) = -(c.at(i));
      } else {
        delta.at(i) = 1;
        carry = false;
      }
    }
    for (; i >= 0; --i) delta.at(i) = 0;
    for (i = 0; i < N && ok_so_far; ++i) {
      c.at(i) += delta.at(i);
      ok_so_far = p.set_value(i, c.at(i));
    }
    if (!ok_so_far) continue;
    if (!is_last) {
      is_last = (c == last);
    } else {
      after_last = true;
    }
    print_newline = true;
    for (i = N-1; i >= 0 && print_newline; --i) {
      if (delta.at(i) < 0 || after_last) {
        out << "\n";
      } else {
        print_newline = false;
      }
    }
  }
  return ok_so_far;
}

bool imaging::grayscale::Image::set_value(
    const imaging::Position &position, const int value) {
  return data_.set_value(position, value);
}

const imaging::Size& imaging::grayscale::Image::size() const {
  return data_.size();
}

bool imaging::grayscale::Image::value(
    const imaging::Position &position, int *value) const {
  if (value == NULL) return false;
  return data_.value(position, value);
}

imaging::grayscale::Image::Image() : data_(imaging::Size(), -1) {
  ; // empty
}

// imaging::binary::_internal::BitMatrix

imaging::binary::_internal::BitMatrix::BitMatrix(
    const imaging::Size &size,
    const bool empty)
    : imaging::NDimensionalMatrixInterface(size, true), blocks_(0) {
  SetSize(size, empty);
}

imaging::binary::_internal::BitMatrix::BitMatrix(
    const imaging::binary::_internal::BitMatrix &other)
    : imaging::NDimensionalMatrixInterface(other), array_(other.array_),
      blocks_(other.blocks_) {
  ; // empty
}

imaging::binary::_internal::BitMatrix::~BitMatrix() {
  array_.clear();
}

imaging::binary::_internal::BitMatrix&
    imaging::binary::_internal::BitMatrix::operator= (
    const imaging::binary::_internal::BitMatrix &other) {
  if (this != &other) {
    imaging::NDimensionalMatrixInterface::operator=(other);
    array_.clear();
    array_ = other.array_;
    blocks_ = other.blocks_;
  }
  return *this;
}

bool imaging::binary::_internal::BitMatrix::CopyFrom(
    const imaging::binary::_internal::BitMatrix &other) {
  if (this == &other) return true;
  if (!imaging::NDimensionalMatrixInterface::CopyFrom(other)) return false;
  array_.clear();
  array_ = other.array_;
  blocks_ = other.blocks_;
  return true;
}

bool imaging::binary::_internal::BitMatrix::Equals(
    const imaging::binary::_internal::BitMatrix &other) const {
  if (this == &other) return true;
  if (!imaging::NDimensionalMatrixInterface::Equals(other)) return false;
  long i = 0;
  if (blocks_ != other.blocks_) return false;
  for (i = 0; i < blocks_; ++i) {
    if (array_.at(i) != (other.array_).at(i)) return false;
  }
  return true;
}

bool imaging::binary::_internal::BitMatrix::InvertValues() {
  long i = 0;
  for (i = 0; i < blocks_; ++i) array_.at(i) = ~(array_.at(i));
  return true;
}

bool imaging::binary::_internal::BitMatrix::set_value(
    const imaging::Position &position, const bool value) {
  imaging::ImagePositionIndex bit_index = 0;
  imaging::binary::_internal::BLOCK bit_value = 0;
  imaging::binary::_internal::BLOCK block = 0;
  imaging::ImagePositionIndex block_index = 0;
  bool ok_so_far = true;
  ok_so_far = CalculateIndexes(position, &block_index, &bit_index);
  if (!ok_so_far) return ok_so_far;
  block = array_.at(block_index);
  bit_value = 1<<(bit_index);
  if (value) {
    block |= bit_value;
  } else {
    block &= ~bit_value;
  }
  array_.at(block_index) = block;
  return ok_so_far;
}

bool imaging::binary::_internal::BitMatrix::value(
    const imaging::Position &position, bool *value) const {
  if (value == NULL) return false;
  imaging::ImagePositionIndex bit_index = 0;
  imaging::binary::_internal::BLOCK bit_value = 0;
  imaging::binary::_internal::BLOCK block = 0;
  imaging::ImagePositionIndex block_index = 0;
  bool ok_so_far = true;
  ok_so_far = CalculateIndexes(position, &block_index, &bit_index);
  if (!ok_so_far) return ok_so_far;
  block = array_.at(block_index);
  bit_value = 1<<(bit_index);
  *value = (block & bit_value) == bit_value;
  return ok_so_far;
}

imaging::binary::_internal::BitMatrix::BitMatrix() {
  blocks_ = 0;
}

bool imaging::binary::_internal::BitMatrix::SetSize(
    const imaging::Size &size,
    const bool empty) {
  imaging::binary::_internal::BLOCK block_value = 0;
  long blocks = 0;
  long i = 0;
  const long items = size.capacity();
  blocks = items/(::number_of_bits());
  if (items%(::number_of_bits()) != 0) blocks++;
  if (!empty) block_value = ~(block_value);
  array_.clear();
  for (i = 0; i < blocks; ++i) array_.push_back(block_value);
  blocks_ = blocks;
  return true;
}

// imaging::binary::StructuringElement

imaging::binary::StructuringElement::StructuringElement(
    const imaging::BoundingBox &bounding_box, const bool empty)
    : bounding_box_(bounding_box), data_(bounding_box.size(), empty) {
  ; // empty
}

imaging::binary::StructuringElement::StructuringElement(
    const imaging::binary::StructuringElement &other)
    : bounding_box_(other.bounding_box_), data_(other.data_) {
  ; // empty
}

imaging::binary::StructuringElement&
    imaging::binary::StructuringElement::operator= (
    const imaging::binary::StructuringElement &other) {
  if (this != &other) {
    bounding_box_.CopyFrom(other.bounding_box_);
    data_.CopyFrom(other.data_);
  }
  return *this;
}

const imaging::BoundingBox& imaging::binary::StructuringElement::bounding_box()
    const {
  return bounding_box_;
}

bool imaging::binary::StructuringElement::CopyFrom(
    const imaging::binary::StructuringElement &other) {
  if (this == &other) return true;
  bool ok_so_far = true;
  ok_so_far = bounding_box_.CopyFrom(other.bounding_box_);
  if (!ok_so_far) return ok_so_far;
  ok_so_far = data_.CopyFrom(other.data_);
  return ok_so_far;
}

bool imaging::binary::StructuringElement::DelimitedComplement(
    imaging::binary::StructuringElement *result) const {
  if (result == NULL) return false;
  if (this == result) return false;
  bool ok_so_far = true;
  ok_so_far = result->CopyFrom(*this);
  if (!ok_so_far) return ok_so_far;
  ok_so_far = result->data_.InvertValues();
  return ok_so_far;
}

bool imaging::binary::StructuringElement::Equals(
    const imaging::binary::StructuringElement &other) const {
  if (this == &other) return true;
  if (!bounding_box_.Equals(other.bounding_box_)) return false;
  return data_.Equals(other.data_);
}

bool imaging::binary::StructuringElement::Intersection(
    const StructuringElement &other,
    bool *empty,
    imaging::binary::StructuringElement *result) const {
  if (result == NULL || empty == NULL) return false;
  imaging::BoundingBox intersection_bb;
  bool ok_so_far = true;
  bool position_value = true;
  ok_so_far = bounding_box_.Intersection(other.bounding_box_, empty,
                                         &intersection_bb);
  if (!ok_so_far) return ok_so_far;
  if (*empty) return true;
  imaging::binary::StructuringElement output(intersection_bb, true);
  imaging::PositionIterator iterator(intersection_bb);
  ok_so_far = iterator.begin();
  if (!ok_so_far) return ok_so_far;
  do {
    ok_so_far = this->value(iterator.value(), &position_value);
    if (!ok_so_far) continue;
    if (position_value) {
      ok_so_far = other.value(iterator.value(), &position_value);
      if (!ok_so_far) continue;
      if (position_value) {
        ok_so_far = output.set_value(iterator.value(), position_value);
        if (!ok_so_far) continue;
      }
    }
  } while (ok_so_far && iterator.iterate());
  if (!iterator.IsFinished()) ok_so_far = false;
  if (!ok_so_far) return ok_so_far;
  ok_so_far = result->CopyFrom(output);
  return ok_so_far;
}

bool imaging::binary::StructuringElement::IsPositionValid(
    const imaging::Position &position) const {
  return bounding_box_.IsValid(position);
}


long imaging::binary::StructuringElement::Length(const char index) const {
  return bounding_box_.Length(index);
}

bool imaging::binary::StructuringElement::ReflectByOrigin(
    imaging::binary::StructuringElement *result) const {
  if (result == NULL) return false;
  bool ok_so_far = true;
  bool position_value = true;
  imaging::BoundingBox reflected_bb;
  imaging::Position reflected_current;
  ok_so_far = bounding_box_.ReflectByOrigin(&reflected_bb);
  if (!ok_so_far) return ok_so_far;
  imaging::binary::StructuringElement output(reflected_bb, true);
  imaging::PositionIterator iterator(bounding_box_);
  ok_so_far = iterator.begin();
  if (!ok_so_far) return ok_so_far;
  do {
    ok_so_far = value(iterator.value(), &position_value);
    if (!ok_so_far) continue;
    if (position_value) {
      ok_so_far = reflected_current.CopyOppositeOf(iterator.value());
      if (!ok_so_far) continue;
      ok_so_far = output.set_value(reflected_current, true);
    }
  } while (ok_so_far && iterator.iterate());
  if (!iterator.IsFinished()) ok_so_far = false;
  if (!ok_so_far) return ok_so_far;
  ok_so_far = result->CopyFrom(output);
  return ok_so_far;
}

bool imaging::binary::StructuringElement::set_value(
    const imaging::Position &position, const bool value) {
  bool ok_so_far = true;
  ok_so_far = position.Subtract(bounding_box_.lower(), &bit_matrix_position_);
  if (!ok_so_far) return ok_so_far;
  return data_.set_value(bit_matrix_position_, value);
}

bool imaging::binary::StructuringElement::SetMinus(
    const imaging::binary::StructuringElement &to_be_subtracted,
    imaging::binary::StructuringElement *result) const {
  if (result == NULL) return false;
  bool empty = true;
  imaging::BoundingBox intersection_bb;
  bool ok_so_far = true;
  bool position_value = true;
  imaging::binary::StructuringElement output(*this);
  ok_so_far = bounding_box_.Intersection(to_be_subtracted.bounding_box_,
                                         &empty, &intersection_bb);
  if (!ok_so_far) return ok_so_far;
  if (empty) return true;
  imaging::PositionIterator iterator(intersection_bb);
  ok_so_far = iterator.begin();
  if (!ok_so_far) return ok_so_far;
  do {
    ok_so_far = to_be_subtracted.value(iterator.value(), &position_value);
    if (!ok_so_far) continue;
    if (position_value) {
      ok_so_far = output.set_value(iterator.value(), false);
      if (!ok_so_far) continue;
    }
  } while (ok_so_far && iterator.iterate());
  if (!iterator.IsFinished()) ok_so_far = false;
  if (!ok_so_far) return ok_so_far;
  ok_so_far = result->CopyFrom(output);
  return ok_so_far;
}

const imaging::Size& imaging::binary::StructuringElement::size() const {
  return bounding_box_.size();
}

bool imaging::binary::StructuringElement::Union(
    const imaging::binary::StructuringElement &other,
    imaging::binary::StructuringElement *result) const {
  if (result == NULL) return false;
  bool ok_so_far = true;
  bool position_value = true;
  imaging::BoundingBox union_bb;
  ok_so_far = bounding_box_.Union(other.bounding_box_, &union_bb);
  if (!ok_so_far) return ok_so_far;
  imaging::binary::StructuringElement output(union_bb, true);
  imaging::PositionIterator iterator(union_bb);
  ok_so_far = iterator.begin();
  if (!ok_so_far) return ok_so_far;
  do {
    ok_so_far = this->value(iterator.value(), &position_value);
    if (ok_so_far && position_value) {
      ok_so_far = output.set_value(iterator.value(), position_value);
      if (!ok_so_far) continue;
    } else {
      ok_so_far = other.value(iterator.value(), &position_value);
      if (ok_so_far && position_value) {
        ok_so_far = output.set_value(iterator.value(), position_value);
        if (!ok_so_far) continue;
      }
    }
    ok_so_far = true;
  } while (ok_so_far && iterator.iterate());
  if (!iterator.IsFinished()) ok_so_far = false;
  if (!ok_so_far) return ok_so_far;
  ok_so_far = result->CopyFrom(output);
  return ok_so_far;
}

bool imaging::binary::StructuringElement::value(
    const imaging::Position &position, bool *value) const {
  if (value == NULL) return false;
  bool ok_so_far = true;
  ok_so_far = position.Subtract(bounding_box_.lower(), &bit_matrix_position_);
  if (!ok_so_far) return ok_so_far;
  return data_.value(bit_matrix_position_, value);
}

imaging::binary::StructuringElement::StructuringElement()
    : bounding_box_(imaging::BoundingBox()),
      data_(imaging::BoundingBox().size(), true) {
  ; // empty
}

// imaging::binary::Image

imaging::binary::Image::Image(const imaging::Size &size, const bool empty)
    : imaging::binary::StructuringElement(size, empty) {
  ; // empty
}

imaging::binary::Image::Image(const imaging::binary::Image &other)
    : imaging::binary::StructuringElement(other) {
  ; // empty
}

imaging::binary::Image& imaging::binary::Image::Image::operator= (
    const imaging::binary::Image &other) {
  if (this != &other) imaging::binary::StructuringElement::operator=(other);
  return *this;
}

bool imaging::binary::Image::CopyFrom(
    const imaging::binary::StructuringElement &other) {
  if (this == &other) return true;
  imaging::Size limited_size;
  bool ok_so_far = true;
  bool position_value = true;
  ok_so_far = limited_size.CopyFrom(other.bounding_box());
  if (!ok_so_far) return ok_so_far;
  imaging::binary::Image to_be_copied(limited_size, true);
  imaging::PositionIterator iterator(limited_size);
  ok_so_far = iterator.begin();
  if (!ok_so_far) return ok_so_far;
  do {
    ok_so_far = other.value(iterator.value(), &position_value);
    if (!ok_so_far) continue;
    if (position_value) {
      ok_so_far = to_be_copied.set_value(iterator.value(), position_value);
    }
  } while (ok_so_far && iterator.iterate());
  if (!iterator.IsFinished()) ok_so_far = false;
  if (!ok_so_far) return ok_so_far;
  ok_so_far = imaging::binary::StructuringElement::CopyFrom(to_be_copied);
  return ok_so_far;
}

bool imaging::binary::Image::ReflectByOrigin(
    imaging::binary::StructuringElement *result) const {
  if (result == NULL) return false;
  long coordinate_value = 0;
  char i = 0;
  long maximum = 0;
  const char n = imaging::Dimension::number();
  bool ok_so_far = true;
  bool position_value = true;
  imaging::Position reflected_current;
  const imaging::Size &size = this->size();
  const imaging::Position &upper = size.upper();
  imaging::binary::Image output(size, true);
  if (n < 1) return false;
  imaging::PositionIterator iterator(size);
  ok_so_far = iterator.begin();
  if (!ok_so_far) return ok_so_far;
  do {
    ok_so_far = value(iterator.value(), &position_value);
    if (!ok_so_far) continue;
    if (position_value) {
      const imaging::Position &current = iterator.value();
      for (i = 0; ok_so_far && i < n; ++i) {
        ok_so_far = current.value(i, &coordinate_value);
        if (!ok_so_far) continue;
        ok_so_far = upper.value(i, &maximum);
        if (!ok_so_far) continue;
        ok_so_far =
            reflected_current.set_value(i, maximum-coordinate_value);
      }
      if (!ok_so_far) continue;
      ok_so_far = output.set_value(reflected_current, true);
    }
  } while (ok_so_far && iterator.iterate());
  if (!ok_so_far) return ok_so_far;
  ok_so_far = result->CopyFrom(output);
  return ok_so_far;
}

bool imaging::binary::Image::Union(
    const imaging::binary::StructuringElement &other,
    imaging::binary::StructuringElement *result) const {
  if (result == NULL) return false;
  bool ok_so_far = true;
  bool position_value = true;
  imaging::BoundingBox union_bb;
  imaging::Size union_size;
  const imaging::BoundingBox &bounding_box = this->bounding_box();
  const imaging::BoundingBox &other_bounding_box = other.bounding_box();
  ok_so_far = bounding_box.Union(other_bounding_box, &union_bb);
  if (!ok_so_far) return ok_so_far;
  ok_so_far = union_size.CopyFrom(union_bb);
  if (!ok_so_far) return ok_so_far;
  imaging::binary::StructuringElement output(union_size, true);
  imaging::PositionIterator iterator(union_size);
  ok_so_far = iterator.begin();
  if (!ok_so_far) return ok_so_far;
  do {
    ok_so_far = this->value(iterator.value(), &position_value);
    if (ok_so_far && position_value) {
      ok_so_far = output.set_value(iterator.value(), position_value);
      if (!ok_so_far) continue;
    } else {
      ok_so_far = other.value(iterator.value(), &position_value);
      if (ok_so_far && position_value) {
        ok_so_far = output.set_value(iterator.value(), position_value);
        if (!ok_so_far) continue;
      }
    }
    ok_so_far = true;
  } while (ok_so_far && iterator.iterate());
  if (!iterator.IsFinished()) ok_so_far = false;
  if (!ok_so_far) return ok_so_far;
  ok_so_far = result->CopyFrom(output);
  return ok_so_far;
}

// imaging::binary::morphology::DualOperation

imaging::binary::morphology::DualOperation::~DualOperation() {
  DualOperation::clear();
}

bool imaging::binary::morphology::DualOperation::DilationTransform(
    const imaging::binary::Image &image,
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
    double *start,
    double *end) {
  return Operation(false, image, se, output,
      algorithm_determinate_border_comparison_counter,
      algorithm_insert_new_candidate_comparison_counter,
      algorithm_insert_new_candidate_memory_access_counter,
      algorithm_remove_candidate_comparison_counter,
      algorithm_remove_candidate_memory_access_counter,
      algorithm_number_of_elements_in_border,
      start, end);
}

bool imaging::binary::morphology::DualOperation::ErosionTransform(
    const imaging::binary::Image &image,
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
    double *start,
    double *end) {
  return Operation(true, image, se, output,
      algorithm_determinate_border_comparison_counter,
      algorithm_insert_new_candidate_comparison_counter,
      algorithm_insert_new_candidate_memory_access_counter,
      algorithm_remove_candidate_comparison_counter,
      algorithm_remove_candidate_memory_access_counter,
      algorithm_number_of_elements_in_border,
      start, end);
}

bool imaging::binary::morphology::DualOperation::clear() {
  // Clear the instance data.
  algorithm_determinate_border_comparison_counter_ = NULL;
  algorithm_insert_new_candidate_comparison_counter_ = NULL;
  algorithm_insert_new_candidate_memory_access_counter_ = NULL;
  algorithm_remove_candidate_comparison_counter_ = NULL;
  algorithm_remove_candidate_memory_access_counter_ = NULL;
  border_.clear();
  border_counter_ = 0;
  candidate_initialized_.clear();
  candidate_next_.clear();
  candidate_position_.clear();
  candidate_previous_.clear();
  se_cardinality_.clear();
  se_elements_.clear();
  se_iteration_ = 0;
  u_elements_.clear();
  if (candidate_matrix_ != NULL) {
    delete candidate_matrix_;
    candidate_matrix_ = NULL;
  }
  if (Y_ != NULL) {
    delete Y_;
    Y_ = NULL;
  }
  return true;
}

bool imaging::binary::morphology::DualOperation::CustomInitialize() {
  return true;
}

bool imaging::binary::morphology::DualOperation::Debug() {
  imaging::ImagePositionIndex counter = 0;
  imaging::ImagePositionIndex current = 0;
  imaging::ImagePositionIndex current_size = 0;
  imaging::ImagePositionIndex i = 0;
  if (!debug_) {
    debug_output_ << "ERROR: debug not set.\n";
    return false;
  }
  if (algorithm_determinate_border_comparison_counter_ == NULL) {
    debug_output_ << "\t\tERROR:"
        << " algorithm_determinate_border_comparison_counter_ is NULL.\n";
    return false;
  } else {
    debug_output_ << "\t\talgorithm_determinate_border_comparison_counter_"
        << " not NULL.\n";
  }
  current_size = static_cast<imaging::ImagePositionIndex>(
      algorithm_determinate_border_comparison_counter_->size());
  if (se_iteration_+1 != current_size) {
    debug_output_ << "\t\tERROR:"
        << " size of algorithm_determinate_border_comparison_counter_ ("
        << current_size << ") differs from se_iteration_+1 ("
        << se_iteration_+1 << ").\n";
    return false;
  } else {
    debug_output_ << "\t\t"
        << "size of algorithm_determinate_border_comparison_counter_ ("
        << current_size << ") matches se_iteration_+1 ("
        << se_iteration_+1 << ").\n";
  }
  if (algorithm_insert_new_candidate_comparison_counter_ == NULL) {
    debug_output_ << "\t\tERROR:"
        << " algorithm_insert_new_candidate_comparison_counter_ is NULL.\n";
    return false;
  } else {
    debug_output_ << "\t\talgorithm_insert_new_candidate_comparison_counter_"
        << " not NULL.\n";
  }
  current_size = static_cast<imaging::ImagePositionIndex>(
      algorithm_insert_new_candidate_comparison_counter_->size());
  if (se_iteration_+1 != current_size) {
    debug_output_ << "\t\tERROR:"
        << " size of algorithm_insert_new_candidate_comparison_counter_ ("
        << current_size << ") differs from se_iteration_+1 ("
        << se_iteration_+1 << ").\n";
    return false;
  } else {
    debug_output_ << "\t\t"
        << "size of algorithm_insert_new_candidate_comparison_counter_ ("
        << current_size << ") matches se_iteration_+1 ("
        << se_iteration_+1 << ").\n";
  }
  if (algorithm_insert_new_candidate_memory_access_counter_ == NULL) {
    debug_output_ << "\t\tERROR:"
        << " algorithm_insert_new_candidate_memory_access_counter_ is NULL.\n";
    return false;
  } else {
    debug_output_ << "\t\talgorithm_insert_new_candidate_memory_access_counter_"
        << " not NULL.\n";
  }
  current_size = static_cast<imaging::ImagePositionIndex>(
      algorithm_insert_new_candidate_memory_access_counter_->size());
  if (se_iteration_+1 != current_size) {
    debug_output_ << "\t\tERROR:"
        << " size of algorithm_insert_new_candidate_memory_access_counter_ ("
        << current_size << ") differs from se_iteration_+1 ("
        << se_iteration_+1 << ").\n";
    return false;
  } else {
    debug_output_ << "\t\t"
        << "size of algorithm_insert_new_candidate_memory_access_counter_ ("
        << current_size << ") matches se_iteration_+1 ("
        << se_iteration_+1 << ").\n";
  }
  if (algorithm_remove_candidate_comparison_counter_ == NULL) {
    debug_output_ << "\t\tERROR:"
        << " algorithm_remove_candidate_comparison_counter_ is NULL.\n";
    return false;
  } else {
    debug_output_ << "\t\talgorithm_remove_candidate_comparison_counter_"
        << " not NULL.\n";
  }
  current_size = static_cast<imaging::ImagePositionIndex>(
      algorithm_remove_candidate_comparison_counter_->size());
  if (se_iteration_+1 != current_size) {
    debug_output_ << "\t\tERROR:"
        << " size of algorithm_remove_candidate_comparison_counter_ ("
        << current_size << ") differs from se_iteration_+1 ("
        << se_iteration_+1 << ").\n";
    return false;
  } else {
    debug_output_ << "\t\t"
        << "size of algorithm_remove_candidate_comparison_counter_ ("
        << current_size << ") matches se_iteration_+1 ("
        << se_iteration_+1 << ").\n";
  }
  if (algorithm_remove_candidate_memory_access_counter_ == NULL) {
    debug_output_ << "\t\tERROR:"
        << " algorithm_remove_candidate_memory_access_counter_ is NULL.\n";
    return false;
  } else {
    debug_output_ << "\t\talgorithm_remove_candidate_memory_access_counter_"
        << " not NULL.\n";
  }
  current_size = static_cast<imaging::ImagePositionIndex>(
      algorithm_remove_candidate_memory_access_counter_->size());
  if (se_iteration_+1 != current_size) {
    debug_output_ << "\t\tERROR:"
        << " size of algorithm_remove_candidate_memory_access_counter_ ("
        << current_size << ") differs from se_iteration_+1 ("
        << se_iteration_+1 << ").\n";
    return false;
  } else {
    debug_output_ << "\t\t"
        << "size of algorithm_remove_candidate_memory_access_counter_ ("
        << current_size << ") matches se_iteration_+1 ("
        << se_iteration_+1 << ").\n";
  }
  if (use_candidate_matrix_ && candidate_matrix_ == NULL) {
    debug_output_ << "\t\tERROR:"
        << " algorithm uses candidate matrix but it is not set.\n";
    return false;
  }
  current = candidate_next_.at(imaging::HEADER);
  counter = 0;
  while (current != imaging::HEADER) {
    ++counter;
    current = candidate_next_.at(current);
  }
  for (i = 0; i < border_counter_; ++i) {
    if (border_.at(i) == imaging::HEADER) {
      debug_output_ << "\t\tERROR: border must not include header index.\n";
      return false;
    }
    if (!candidate_initialized_.at(border_.at(i))) {
      debug_output_ << "\t\tERROR:"
          << " an uninitialized candidate was found in border.\n";
      return false;
    }
  }
  debug_output_ << "\t\tse_iteration_: " << se_iteration_ << "\n";
  debug_output_ << "\t\tcandidate nodes: " << counter << "\n";
  debug_output_ << "\t\tborder nodes: " << border_counter_ << "\n";
  return true;
}

bool imaging::binary::morphology::DualOperation::EnqueueCandidateNode(
    const imaging::ImagePositionIndex &image_position) {
  const imaging::ImagePositionIndex previous = candidate_previous_.at(
      imaging::HEADER);
  if (image_position == imaging::HEADER) return false;
  if (candidate_initialized_.at(image_position)) {
    return true; // already enqueued
  }
  algorithm_insert_new_candidate_memory_access_counter_->at(se_iteration_)
      += 5;
  candidate_previous_.at(imaging::HEADER) = image_position;
  candidate_next_.at(previous) = image_position;
  candidate_previous_.at(image_position) = previous;
  candidate_next_.at(image_position) = imaging::HEADER;
  candidate_initialized_.at(image_position) = true;
  return true;
}

bool imaging::binary::morphology::DualOperation::position(
    const imaging::Position &image_position,
    imaging::ImagePositionIndex *value) const {
  if (!use_candidate_matrix_ || value == NULL) return false;
  return candidate_matrix_->value(image_position, value);
}

bool imaging::binary::morphology::DualOperation::position(
    const imaging::ImagePositionIndex &position_index,
    imaging::Position *value) const {
  if (value == NULL) return false;
  return value->CopyFrom(candidate_position_.at(position_index));
}

bool imaging::binary::morphology::DualOperation::RemoveCandidateNode(
    const imaging::ImagePositionIndex &image_position) {
  if (image_position == imaging::HEADER) return false;
  if (!candidate_initialized_.at(image_position)) return false;
  const imaging::ImagePositionIndex next = candidate_next_.at(image_position);
  const imaging::ImagePositionIndex previous =
      candidate_previous_.at(image_position);
  algorithm_remove_candidate_memory_access_counter_->at(se_iteration_) += 4;
  candidate_next_.at(previous) = next;
  candidate_previous_.at(next) = previous;
  candidate_next_.at(image_position) = image_position;
  candidate_previous_.at(image_position) = image_position;
  return true;
}

imaging::SEIndex imaging::binary::morphology::DualOperation::u_cardinality()
    const {
  return static_cast<imaging::SEIndex>(u_elements_.size());
}


bool imaging::binary::morphology::DualOperation::ActualAlgorithm(
    const bool true_for_erosion,
    imaging::grayscale::Image **output_image) {
  imaging::ImagePositionIndex current = 0;
  imaging::SEIndex current_se = 0;
  imaging::SEIndex current_se_index = 0;
  imaging::ImagePositionIndex current_se_index_element = 0;
  imaging::ImagePositionIndex i = 0;
  const imaging::ImagePositionIndex initial_counter_value = 0;
  imaging::SEIndex not_done = 0;
  bool ok_so_far = true;
  const imaging::SEIndex number_of_se = se_elements_.size();
  std::vector< imaging::SEIndex > se_index;
  std::vector< std::vector<imaging::ImagePositionIndex> > shuffled_indexes;
  // Setting up SE element index shuffler's data.
  for (current_se = 0; ok_so_far && current_se < number_of_se; ++current_se) {
    std::vector<imaging::ImagePositionIndex> current_se_indexes;
    const imaging::ImagePositionIndex current_cardinality =
          se_cardinality_.at(current_se);
    for (current_se_index = 0; current_se_index < current_cardinality;
        ++current_se_index) {
      current_se_indexes.push_back(current_se_index);
    }
    shuffled_indexes.push_back(current_se_indexes);
    se_index.push_back(current_se);
  }
  // Transform itself.
  se_iteration_ = 0;
  not_done = 0;
  current = candidate_next_.at(imaging::HEADER);
  while (ok_so_far && current != imaging::HEADER && not_done < number_of_se) {
    // Insert current iteration counter data value.
    algorithm_determinate_border_comparison_counter_->push_back(
        initial_counter_value);
    algorithm_insert_new_candidate_comparison_counter_->push_back(
        initial_counter_value);
    algorithm_insert_new_candidate_memory_access_counter_->push_back(
        initial_counter_value);
    algorithm_remove_candidate_comparison_counter_->push_back(
        initial_counter_value);
    algorithm_remove_candidate_memory_access_counter_->push_back(
        initial_counter_value);
    algorithm_number_of_elements_in_border_->push_back(
        initial_counter_value);
    // Set up data with current values.
    ++se_iteration_;
    ok_so_far = shuffle(&se_index); // shuffle SE sequence
    not_done = 0;
    if (debug_) {
      debug_output_ << "\tBefore running iteration " << se_iteration_ << ":\n";
      debug_output_ << "\n\tSE sequence: [ ";
      for (current_se = 0; current_se < number_of_se; ++current_se) {
        debug_output_ << se_index.at(current_se);
        if (current_se != number_of_se-1) debug_output_ << ",";
        debug_output_ << " ";
      }
      debug_output_ << "]\n\n";
      this->Debug();
      debug_output_ << "\n";
    }
    for (current_se = 0;
        ok_so_far && current != imaging::HEADER && current_se < number_of_se;
        ++current_se) {
      border_counter_ = 0;
      current_se_index = se_index.at(current_se);
      border_counter_ = 0;
      // Determinate which candidates belongs to current border.
      ok_so_far = shuffle(&(shuffled_indexes.at(current_se_index)));
      if (!ok_so_far) continue;
      if (debug_) {
        const imaging::ImagePositionIndex se_cardinality
            = se_cardinality_.at(current_se_index);
        const std::vector<imaging::ImagePositionIndex> &se_elements
            = shuffled_indexes.at(current_se_index);
        debug_output_ << "\tSE current_se_index " << current_se_index
            << " [" << current_se+1 << "/" << number_of_se
            << "] elements sequence: [ ";
        for (current_se_index_element = 0;
            current_se_index_element < se_cardinality;
            ++current_se_index_element) {
          debug_output_ << se_elements.at(current_se_index_element);
          if (current_se_index_element != se_cardinality-1)
            debug_output_ << ",";
          debug_output_ << " ";
        }
        debug_output_ << "]\n\n";
      }
      ok_so_far = DetectBorder(true_for_erosion, current_se_index,
          shuffled_indexes.at(current_se_index));
      if (!ok_so_far) continue;
      if (debug_) {
        debug_output_ << "\tAfter determinating border [" << current_se+1
            << "/" << number_of_se << "]:\n";
        this->Debug();
        debug_output_ << "\n";
      }
      // Remove border pixels.
      for (i = 0; ok_so_far && i < border_counter_; ++i) {
        const imaging::Position &p = candidate_position_.at(border_.at(i));
        ok_so_far = Y_->set_value(p, !true_for_erosion);
        if (!ok_so_far) continue;
        if (regular_removal_) {
          ok_so_far = RemoveCandidateNode(border_.at(i));
          if (!ok_so_far) continue;
        }
        algorithm_remove_candidate_memory_access_counter_->at(se_iteration_)
            += 1;
        ok_so_far = (*output_image)->set_value(p, se_iteration_);
        if (!ok_so_far) continue;
      }
      if (!ok_so_far) continue;
      // Insert new candidate pixels into candidate queue.
      ok_so_far = this->InsertNewCandidateFromBorder(true_for_erosion,
          output_image);
      if (!ok_so_far) continue;
      // Verifies if operation was possible for the original SE.
      if (border_counter_ == 0) {
        ++not_done;
      } else {
        not_done = 0;
      }
      algorithm_number_of_elements_in_border_->at(se_iteration_)
          += border_counter_;
      border_counter_ = 0; // setting zero for debug info
      if (debug_) {
        debug_output_ << "\tAfter removing border and inserting new candidates ["
            << current_se+1 << "/" << number_of_se << "]:\n";
        this->Debug();
        debug_output_ << "\n";
      }
      current = candidate_next_.at(imaging::HEADER);
    }
  }
  return ok_so_far;
}

bool imaging::binary::morphology::DualOperation::InitializeCandidateData(
    const bool true_for_erosion,
    const imaging::binary::Image &image) {
  if (Y_ == NULL) return false;
  bool ok_so_far = true;
  imaging::ImagePositionIndex position_counter = 0;
  bool position_value = true;
  // Initialize candidate nodes vectors.
  border_.push_back(imaging::HEADER);
  candidate_initialized_.push_back(false);
  candidate_next_.push_back(imaging::HEADER);
  candidate_position_.push_back(imaging::Position());
  candidate_previous_.push_back(imaging::HEADER);
  // Put each candidate pixel of the image into a list.
  imaging::PositionIterator iterator(image.size());
  ok_so_far = iterator.begin();
  if (!ok_so_far) return ok_so_far;
  do {
    const imaging::Position &current = iterator.value();
    algorithm_insert_new_candidate_comparison_counter_->at(se_iteration_) += 1;
    algorithm_insert_new_candidate_memory_access_counter_->at(se_iteration_)
        += 5;
    ok_so_far = image.value(current, &position_value);
    if (!ok_so_far) continue;
    if (position_value != true_for_erosion) continue;
    ++position_counter;
    border_.push_back(imaging::HEADER);
    candidate_initialized_.push_back(false);
    candidate_next_.push_back(position_counter);
    candidate_position_.push_back(current);
    candidate_previous_.push_back(position_counter);
    if (use_candidate_matrix_) {
      algorithm_insert_new_candidate_memory_access_counter_->at(se_iteration_)
          += 1;
      ok_so_far = candidate_matrix_->set_value(current, position_counter);
      if (!ok_so_far) continue;
    }
    ok_so_far = this->InitialCandidatePositionFound(
        true_for_erosion, image, position_counter, current);
  } while (ok_so_far && iterator.iterate());
  if (!iterator.IsFinished()) ok_so_far = false;
  if (ok_so_far && debug_) {
    debug_output_ << "\tAfter initialize candidate data:\n";
    ok_so_far = this->Debug();
    debug_output_ << "\n";
  }
  return ok_so_far;
}

bool imaging::binary::morphology::DualOperation::InitializeCounters() {
  const imaging::ImagePositionIndex initial_counter_value = 0;
  algorithm_determinate_border_comparison_counter_->push_back(
      initial_counter_value);
  algorithm_insert_new_candidate_comparison_counter_->push_back(
      initial_counter_value);
  algorithm_insert_new_candidate_memory_access_counter_->push_back(
      initial_counter_value);
  algorithm_remove_candidate_comparison_counter_->push_back(
      initial_counter_value);
  algorithm_remove_candidate_memory_access_counter_->push_back(
      initial_counter_value);
  algorithm_number_of_elements_in_border_->push_back(
      initial_counter_value);
  return true;
}


bool imaging::binary::morphology::DualOperation::InitializeSEData(
    const std::vector< std::vector<imaging::Position> > &vectorized_se) {
  imaging::ImagePositionIndex i = 0;
  imaging::SEIndex i_se = 0;
  imaging::ImagePositionIndex j_se = 0;
  const int number_of_se = vectorized_se.size();
  bool ok_so_far = true;
  bool position_value = true;
  std::vector<imaging::Position> q;
  imaging::ImagePositionIndex u_cardinality = 0;
  imaging::BoundingBox union_bb;
  // Calculate the bounding box of the union of the SE.
  for (i_se = 0; ok_so_far && i_se < number_of_se; ++i_se) {
    const std::vector<imaging::Position> &current_se = vectorized_se.at(i_se);
    const long current_cardinality = current_se.size();
    for (j_se = 0; ok_so_far && j_se < current_cardinality; ++j_se) {
      const imaging::Position &current_position = current_se.at(j_se);
      ok_so_far = union_bb.Expand(current_position);
    }
  }
  if (!ok_so_far) return ok_so_far;
  imaging::binary::StructuringElement U(union_bb, true);
  se_cardinality_.clear();
  se_cardinality_.resize(number_of_se, -1);
  // Iterate SEs, counting the cardinality of the union of SEs.
  for (i_se = 0; ok_so_far && i_se < number_of_se; ++i_se) {
    const std::vector<imaging::Position> &current_se = vectorized_se.at(i_se);
    const imaging::ImagePositionIndex current_cardinality = current_se.size();
    se_cardinality_.at(i_se) = current_cardinality;
    for (j_se = 0; ok_so_far && j_se < current_cardinality; ++j_se) {
      const imaging::Position &current_position = current_se.at(j_se);
      ok_so_far = U.value(current_position, &position_value);
      if (!ok_so_far) continue;
      if (position_value) continue;
      ++u_cardinality;
      ok_so_far = U.set_value(current_position, true);
      q.push_back(current_position);
    }
  }
  if (!ok_so_far) return ok_so_far;
  u_elements_.clear();
  u_elements_.resize(u_cardinality);
  if (u_cardinality != static_cast<imaging::ImagePositionIndex>(q.size())) {
    ok_so_far = false;
  }
  if (!ok_so_far) return ok_so_far;
  // Put the position of each foreground pixel of the union of SEs
  // into an index of 'u_elements'.
  for (i = 0; i < u_cardinality; ++i) u_elements_.at(i) = q.at(i);
  se_elements_.clear();
  se_elements_.resize(number_of_se);
  for (i_se = 0; ok_so_far && i_se < number_of_se; ++i_se) {
    // Put the index of foreground pixels of each SE, according to 'u_elements'
    // into a vector in 'se_elements'.
    imaging::BoundingBox bb;
    const imaging::ImagePositionIndex &current_cardinality =
        se_cardinality_.at(i_se);
    imaging::SEIndex j = 0;
    imaging::SEIndex m = 0;
    std::vector< imaging::ImagePositionIndex > se_indexes =
        se_elements_.at(i_se);
    // Temporarily creates a SE for direct position access.
    for (j_se = 0; ok_so_far && j_se < current_cardinality; ++j_se) {
      const std::vector<imaging::Position> &current_se = vectorized_se.at(i_se);
      const imaging::Position &current_position = current_se.at(j_se);
      ok_so_far = bb.Expand(current_position);
    }
    if (!ok_so_far) continue;
    imaging::binary::StructuringElement actual_se(bb, true);
    for (j_se = 0; ok_so_far && j_se < current_cardinality; ++j_se) {
      const std::vector<imaging::Position> &current_se = vectorized_se.at(i_se);
      const imaging::Position &current_position = current_se.at(j_se);
      ok_so_far = actual_se.set_value(current_position, true);
    }
    if (!ok_so_far) continue;
    se_indexes.clear();
    se_indexes.resize(current_cardinality, -1);
    while (ok_so_far && m < current_cardinality && j < u_cardinality) {
      const imaging::Position &current = u_elements_.at(j);
      if (actual_se.IsPositionValid(current)) {
        ok_so_far = actual_se.value(current, &position_value);
        if (!ok_so_far) continue;
        if (position_value) {
          se_indexes.at(m) = j;
          ++m;
        }
      }
      ++j;
    }
    if (m != current_cardinality) ok_so_far = false;
    se_elements_.at(i_se) = se_indexes;
  }
  return ok_so_far;
}

bool imaging::binary::morphology::DualOperation::Operation(
    const bool true_for_erosion,
    const imaging::binary::Image &image,
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
    double *start,
    double *end) {
  if (output == NULL
      || algorithm_determinate_border_comparison_counter == NULL
      || algorithm_insert_new_candidate_comparison_counter == NULL
      || algorithm_insert_new_candidate_memory_access_counter == NULL
      || algorithm_remove_candidate_comparison_counter == NULL
      || algorithm_remove_candidate_memory_access_counter == NULL
      || algorithm_number_of_elements_in_border == NULL
      || start == NULL || end == NULL)
    return false;
  if (*output != NULL) return false;
  const size_t empty = 0;
  bool ok_so_far = true;
  std::vector< std::vector<imaging::Position> > vectorized_se;
  // Clear the instance data.
  this->clear();
  // Set current algorithm counter data.
  if (algorithm_determinate_border_comparison_counter->size() != empty
      || algorithm_insert_new_candidate_comparison_counter->size() != empty
      || algorithm_insert_new_candidate_memory_access_counter->size() != empty
      || algorithm_remove_candidate_comparison_counter->size() != empty
      || algorithm_remove_candidate_memory_access_counter->size() != empty
      || algorithm_number_of_elements_in_border->size() != empty)
    return false;
  algorithm_determinate_border_comparison_counter_ =
      algorithm_determinate_border_comparison_counter;
  algorithm_insert_new_candidate_comparison_counter_ =
      algorithm_insert_new_candidate_comparison_counter;
  algorithm_insert_new_candidate_memory_access_counter_ =
      algorithm_insert_new_candidate_memory_access_counter;
  algorithm_remove_candidate_comparison_counter_ =
      algorithm_remove_candidate_comparison_counter;
  algorithm_remove_candidate_memory_access_counter_ =
      algorithm_remove_candidate_memory_access_counter;
  algorithm_number_of_elements_in_border_ =
      algorithm_number_of_elements_in_border;
  // Initialize candidate matrix.
  if (use_candidate_matrix_) {
    using namespace imaging;
    using namespace imaging::grayscale::_internal;
    candidate_matrix_ =
        new NumericalMatrix<ImagePositionIndex>(image.size(), imaging::HEADER);
    if (candidate_matrix_ == NULL) return false;
  }
  // Initialize temporary image.
  Y_ = new imaging::binary::Image(image);
  if (Y_ == NULL) return false;
  // Initialize transitional output image.
  ok_so_far = ::InitializeAlgorithmsOutputImage(image, output);
  if (!ok_so_far) return ok_so_far;
  // Vectorize SEs.
  ok_so_far = ::VectorizeSEs(se, &vectorized_se);
  if (!ok_so_far) return ok_so_far;
  // Start!
  gettimeofday(&(::timer), NULL);
  *start = (::timer).tv_sec*1000000.+(::timer).tv_usec;
  // Initialize SE data.
  ok_so_far = InitializeSEData(vectorized_se);
  if (!ok_so_far) return ok_so_far;
  // Initialize algorithm memory and comparison counters.
  ok_so_far = InitializeCounters();
  if (!ok_so_far) return ok_so_far;
  // Initialize according to chosen algorithm.
  ok_so_far = this->CustomInitialize();
  if (!ok_so_far) return ok_so_far;
  // Initialize candidate data.
  ok_so_far = InitializeCandidateData(true_for_erosion, image);
  if (!ok_so_far) return ok_so_far;
  // Actual algorithm!
  border_.resize(candidate_position_.size(), imaging::HEADER);
  ok_so_far = this->ActualAlgorithm(true_for_erosion, output);
  if (!ok_so_far) return ok_so_far;
  if (debug_) {
    debug_output_ << "\tAfter actual algorithm:\n";
    ok_so_far = this->Debug();
    if (!ok_so_far) return ok_so_far;
    debug_output_ << "\n";
  }
  // Finally!
  gettimeofday(&(::timer), NULL);
  *end = (::timer).tv_sec*1000000.+(::timer).tv_usec;
  // Clear the instance data.
  this->clear();
  return ok_so_far;
}
