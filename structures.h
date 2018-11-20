/******************************************************************************
 ******************************************************************************
 ** structures.h (part of the easymifs program)                              **
 ** Author:  Dario Ghersi                                                    **
 ** Version: 011210                                                          **
 ******************************************************************************
 ******************************************************************************/
#include "constants.h"

typedef struct point { /* linked list of the grid points */
  double x; /* cartesian coordinates of the point */
  double y;
  double z;
  double energy[MAX_NUM_PROBES]; /* interaction energies (one for each probe) */
  struct point *next; /* pointer to the next point */
} Point;

/* -------------------------------------------------------------------------- */

typedef struct box {
  double center[3]; /* center of the box */
  double resolution; /* resolution of the grid */
  unsigned int npoints[3]; /* number of points in x, y, z */
  double lower_corner[3]; /* coordinates of the lower corner of the box */
} Box;

/* -------------------------------------------------------------------------- */

typedef struct atom { /* linked list of the atoms in the protein */
  double x; /* cartesian coordinates of the atom */
  double y;
  double z;
  double partial_charge;
  double *c6; /* the parameters for the Lennard Jones term */
  double *c12;
  char type[TYPE_LEN + 1]; /* must match the force field table */
  struct atom *next; /* pointer to the next atom */
} Atom;

/* -------------------------------------------------------------------------- */

typedef struct atom_type { /* linked list of the atom types in the ff */
  char type[TYPE_LEN + 1];
  bool has_charge; /* boolean variable to check if charge is defined */
  double charge;
  struct atom_type *next; /* pointer to the next atom_type */
} Atom_type;

/* -------------------------------------------------------------------------- */

typedef struct nb_table { /* table of nonbonded interactions */
  char type1[TYPE_LEN + 1];
  char type2[TYPE_LEN + 1];
  double c6; /* Lennard-Jones C6 */
  double c12; /* Lennard-Jones C12 */
  struct nb_table *next; /* pointer to the next nb_table */
} Nb_table;

/* -------------------------------------------------------------------------- */

typedef struct parameters { /* user-defined parameters */
  char filename[MAX_STR_LEN]; /* .easymif input file */
  char atom_types_filename[MAX_STR_LEN]; /* atom types file name */
  char nb_table[MAX_STR_LEN]; /* the filename with the nonbonded parameters */
  char probes[MAX_NUM_PROBES][MAX_STR_LEN]; /* vector of probes */
  unsigned int num_probes; /* number of probes */
  double resolution; /* grid spacing */
  double clearance; /* clearance from the box sides */
  double center[3]; /*  center of the box */
  unsigned int npoints[3]; /* number of points in x, y, z (must be odd) */
  bool set_box; /* if true then the user-defined parameters for the box are used */
  bool compressed; /* if true then write compressed dx files */
} Parameters;

/* -------------------------------------------------------------------------- */

typedef struct hash {
  char value;
  int first;
  int next;
} Hash;
