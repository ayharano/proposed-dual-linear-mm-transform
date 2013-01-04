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

// This file contains the definition of an algorithm tester for both
// erosion and dilation for 2D images.

#ifndef TEST_H_
#define TEST_H_

int tester(const std::string &file_path, const std::vector<short> &se,
           const int times, const bool be_verbose);

int main(int argc, const char* argv[]);

#endif // TEST_H_
