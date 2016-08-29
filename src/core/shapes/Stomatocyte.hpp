/*
  Copyright (C) 2010,2011,2012,2013,2014 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
  Max-Planck-Institute for Polymer Research, Theory Group
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

#ifndef __STOMATOCYTE_HPP
#define __STOMATOCYTE_HPP

#include "Shape.hpp"

namespace Shapes {
  struct Stomatocyte : public Shape {
    virtual const std::string name() const { return std::string("Stomatocyte"); }
    int calculate_dist(const double *ppos, double *dist, double *vec);

    /** Stomatocyte position. */
    double position_x;
    double position_y;
    double position_z;

    /** Stomatocyte orientation. */
    double orientation_x;
    double orientation_y;
    double orientation_z;

    /** Stomatocyte dimensions. */
    double outer_radius;
    double inner_radius;
    double layer_width;

    /** Inside/Outside (+1 outside -1 inside interaction direction)*/
    double direction;
  };
};

#endif
