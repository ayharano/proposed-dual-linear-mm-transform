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

// This file contains the implementation of 2D image file handlers.

#include <Magick++.h>

#include "img_2d.h"

bool bidimensional::LoadBinaryImage(const std::string &file_path,
                                    imaging::binary::Image **image) {
  if (image == NULL) return false;
  if (*image != NULL) return false;
  char n = imaging::Dimension::number();
  if (n == 0) n = imaging::Dimension::Set(2);
  if (n != 2) return false;
  bool ok_so_far = true;
  bool position_value = true;
  unsigned long u_x = 0;
  unsigned long u_y = 0;
  imaging::Position external_upper;
  imaging::Position p;
  long padding = 0;
  Magick::ColorMono value;
  long x = 0;
  long y = 0;
  Magick::Image original_image(file_path);
  // Set image size.
  const long real_width = static_cast<long>(original_image.columns());
  const long real_height = static_cast<long>(original_image.rows());
  if (real_width < real_height) {
    padding = 3*real_height/2;
  } else {
    padding = 3*real_width/2;
  }
  ok_so_far = external_upper.set_value(0, real_width+2*padding);
  if (!ok_so_far) return ok_so_far;
  ok_so_far = external_upper.set_value(1, real_height+2*padding);
  if (!ok_so_far) return ok_so_far;
  ok_so_far = external_upper.value(0, &x);
  ok_so_far = external_upper.value(1, &y);
  imaging::Size size(external_upper, padding);
  imaging::binary::Image output(size, true);
  imaging::PositionIterator iterator(size);
  ok_so_far = iterator.begin();
  if (!ok_so_far) return ok_so_far;
  // Copy image data.
  do {
    const imaging::Position &p = iterator.value();
    ok_so_far = p.value(0, &x);
    if (!ok_so_far) continue;
    if (x < padding || x > real_width+padding) continue;
    ok_so_far = p.value(1, &y);
    if (!ok_so_far) continue;
    if (y < padding || y > real_height+padding) continue;
    u_x = static_cast<unsigned long>(x-padding);
    u_y = static_cast<unsigned long>(y-padding);
    value = original_image.pixelColor(u_x, u_y);
    position_value = value.mono();
    ok_so_far = output.set_value(p, position_value);
    if (!ok_so_far) continue;
  } while (ok_so_far && iterator.iterate());
  if (!iterator.IsFinished()) ok_so_far = false;
  if (!ok_so_far) return ok_so_far;
  *image = new imaging::binary::Image(output);
  if (*image == NULL) return false;
  return ok_so_far;
}

bool bidimensional::SaveBinaryImage(const std::string &file_path,
                                    const imaging::binary::Image &image) {
  return SaveBinaryImage(file_path, 1, image);
}

bool bidimensional::SaveBinaryImage(const std::string &file_path,
                                    const int pixel_size,
                                    const imaging::binary::Image &image) {
  char n = imaging::Dimension::number();
  if (n != 2) return false;
  if (pixel_size < 1) return false;
  const imaging::Size &size = image.size();
  Magick::ColorMono color;
  long d_x = 0;
  long d_y = 0;
  bool ok_so_far = true;
  bool position_value = true;
  unsigned long u_x = 0;
  unsigned long u_y = 0;
  long x = 0;
  long y = 0;
  // Get image size.
  const long width = size.Length(0);
  if (width < 1) return false;
  const long height = size.Length(1);
  if (height < 1) return false;
  const long padding = size.padding();
  if (padding < 0) return false;
  if (2*padding >= width) return false;
  if (2*padding >= height) return false;
  Magick::Geometry magick_size(
      static_cast<size_t>(pixel_size*(width-2*padding)),
      static_cast<size_t>(pixel_size*(height-2*padding)));
  Magick::Image output(magick_size, Magick::ColorMono(false));
  imaging::PositionIterator iterator(size);
  ok_so_far = iterator.begin();
  if (!ok_so_far) return ok_so_far;
  // Copy image data.
  do {
    const imaging::Position &p = iterator.value();
    ok_so_far = p.value(0, &x);
    if (!ok_so_far) continue;
    if (x < padding || x >= width+padding) continue;
    ok_so_far = p.value(1, &y);
    if (!ok_so_far) continue;
    if (y < padding || y >= height+padding) continue;
    position_value = false;
    ok_so_far = image.value(p, &position_value);
    if (!ok_so_far) continue;
    if (!position_value) continue;
    color.mono(position_value);
    for (d_x = 0; d_x < pixel_size; ++d_x) {
      u_x = static_cast<unsigned long>(pixel_size*(x-padding)+d_x);
      for (d_y = 0; d_y < pixel_size; ++d_y) {
        u_y = static_cast<unsigned long>(pixel_size*(y-padding)+d_y);
        output.pixelColor(u_x, u_y, color);
      }
    }
  } while (ok_so_far && iterator.iterate());
  if (!iterator.IsFinished()) ok_so_far = false;
  if (!ok_so_far) return ok_so_far;
  output.write(file_path);
  return ok_so_far;
}

bool bidimensional::SaveGrayscaleImage(const std::string &file_path,
                                       const imaging::grayscale::Image &image) {
  return SaveGrayscaleImage(file_path, 1, image);
}

bool bidimensional::SaveGrayscaleImage(const std::string &file_path,
                                       const int pixel_size,
                                       const imaging::grayscale::Image &image) {
  char n = imaging::Dimension::number();
  if (n != 2) return false;
  if (pixel_size < 1) return false;
  const imaging::Size &size = image.size();
  Magick::ColorGray color;
  long d_x = 0;
  long d_y = 0;
  const int levels = 256;
  bool ok_so_far = true;
  int position_value = true;
  unsigned long u_x = 0;
  unsigned long u_y = 0;
  long x = 0;
  long y = 0;
  // Get image size.
  const long width = size.Length(0);
  if (width < 1) return false;
  const long height = size.Length(1);
  if (height < 1) return false;
  const long padding = size.padding();
  if (padding < 0) return false;
  if (2*padding >= width) return false;
  if (2*padding >= height) return false;
  Magick::Geometry magick_size(
      static_cast<size_t>(pixel_size*(width-2*padding)),
      static_cast<size_t>(pixel_size*(height-2*padding)));
  Magick::Image output(magick_size, Magick::ColorGray(0.));
  imaging::PositionIterator iterator(size);
  ok_so_far = iterator.begin();
  if (!ok_so_far) return ok_so_far;
  // Copy image data.
  do {
    const imaging::Position &p = iterator.value();
    ok_so_far = p.value(0, &x);
    if (!ok_so_far) continue;
    if (x < padding || x >= width+padding) continue;
    ok_so_far = p.value(1, &y);
    if (!ok_so_far) continue;
    if (y < padding || y >= height+padding) continue;
    position_value = 0;
    ok_so_far = image.value(p, &position_value);
    if (!ok_so_far) continue;
    ++position_value; // adjustment of pixel level
    color.shade((1.*position_value)/(1.*(levels-1)));
    for (d_x = 0; d_x < pixel_size; ++d_x) {
      u_x = static_cast<unsigned long>(pixel_size*(x-padding)+d_x);
      for (d_y = 0; d_y < pixel_size; ++d_y) {
        u_y = static_cast<unsigned long>(pixel_size*(y-padding)+d_y);
        output.pixelColor(u_x, u_y, color);
      }
    }
  } while (ok_so_far && iterator.iterate());
  if (!iterator.IsFinished()) ok_so_far = false;
  if (!ok_so_far) return ok_so_far;
  output.quantizeColorSpace( Magick::GRAYColorspace );
  output.quantizeColors( levels );
  output.quantize( );
  output.write(file_path);
  return ok_so_far;
}
