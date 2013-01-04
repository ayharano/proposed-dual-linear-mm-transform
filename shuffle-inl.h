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

// This file contains the implementation of inline function described
// as Algorithm P.

#ifndef SHUFFLE_INL_H_
#define SHUFFLE_INL_H_

#include <cmath>
#include <cstdlib>
#include <vector>

template< class T >
inline bool shuffle(std::vector<T> *data) {
  if (data == NULL) return false;
  long j = static_cast<long>(data->size())-1;
  long k = 0;
  T tmp;
  double U = 0.;
  while (j > 0) {
    U = (1.*rand())/(1.*RAND_MAX);
    k = static_cast<long>(floor(j*U));
    tmp = data->at(j);
    data->at(j) = data->at(k);
    data->at(k) = tmp;
    --j;
  }
  return true;
}

#endif // SHUFFLE_INL_H_
