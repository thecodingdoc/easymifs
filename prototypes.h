/******************************************************************************
 ******************************************************************************
 ** protypes.h (part of the easymifs program)                                **
 ** Author:  Dario Ghersi                                                    **
 ** Version: 20160418                                                        **
 ******************************************************************************
 ******************************************************************************/

#include "structures.h"

#define VERSION "20160418"

double euclidean_distance(Point *, Atom *);
void set_parameters(int, char **, Parameters *);
Atom_type *read_atom_types(Parameters);
Nb_table *read_nb_table(Parameters);
Atom *extract_atoms(Parameters, Nb_table *);
Box calculate_box(Atom *, Parameters);
Box set_box(Parameters);
Point *inizialize_point_list(Box, Parameters);
void calculate_MIFs(Point *, Atom *, Nb_table *, Atom_type *, Parameters);
void compress(FILE *, double *, unsigned int);
void create_dx(Parameters, Point *, Box);
void free_memory_atom_types_list(Atom_type *);
void free_memory_nb_table(Nb_table *);
void free_memory_atoms(Atom *);
void free_memory_atom_types(Atom_type *);
void free_memory_points(Point *);
