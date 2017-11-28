from __future__ import print_function, absolute_import
cdef extern from "particle_data.hpp":
    int init_type_map(int type) except +
    int find_particle_type(int type, int * id)  except +
    int number_of_particles_with_type(int type, int * number)  except +
