#pragma once

#include <iostream>
#include <ctime>
#include <cstdlib>
#include <string>
#include <string_view>
#include <bitset>
#include <chrono>
#include <vector>
#include <assert.h>
#include <fstream>
#include <array>
#include <map>
#include <mutex>
#include <atomic>

#include "embedding.h"

#include <boost/program_options.hpp>
#include <boost/algorithm/string/case_conv.hpp>

#include "ksw2.h"

#include "tbb/tbb.h"
#include "tbb/task_scheduler_init.h"
#include "tbb/parallel_sort.h"

#include "type.h"
#include "read.h"
#include "const.h"

//extern "C"{
//#include "gap_affine/affine_wavefront_align.h"
//}
