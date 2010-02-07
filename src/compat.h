#ifndef COMPAT_H
#define COMPAT_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_UNORDERED_MAP
#include <unordered_map>
using std::unordered_map;
#else
#ifdef HAVE_TR1_UNORDERED_MAP
#include <tr1/unordered_map>
using std::tr1::unordered_map;
#endif
#endif

#endif
