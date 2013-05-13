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

// This file contains the declaration of 2D image file handlers.

#ifndef IMG_2D_H_
#define IMG_2D_H_

#include <string>

#include "img.h"

namespace bidimensional {

  bool LoadBinaryImage(const std::string &file_path,
                       imaging::binary::Image **image_d,
                       imaging::binary::Image **image_e);

  bool SaveBinaryImage(const std::string &file_path,
                       const imaging::binary::Image &image);

  bool SaveBinaryImage(const std::string &file_path, const int pixel_size,
                       const imaging::binary::Image &image);

  bool SaveGrayscaleImage(const std::string &file_path,
                          const imaging::grayscale::Image &image);

  bool SaveGrayscaleImage(const std::string &file_path, const int pixel_size,
                          const imaging::grayscale::Image &image);

} // namespace bidimensional

#endif // IMG_2D_H_
