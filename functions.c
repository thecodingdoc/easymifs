/******************************************************************************
 ******************************************************************************
 ** functions.c (part of the "easymifs" program)                             **
 **                                                                          **
 ** Author:  Dario Ghersi                                                    **
 ** Version: 011210                                                          **
 ** Goal:    contains the core functions of easyMIFs                         **
 ******************************************************************************
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "structures.h"

/******************************************************************************
 * FUNCTIONS                                                                  *
 ******************************************************************************/

char *basename(char *s)
{
  /* return */
  char *b;
  char p[1000];
  char *temp;

  b = malloc(sizeof(char) * MAX_STR_LEN);
  temp = strtok(s, ".");
  strncpy(b, temp, MAX_STR_LEN);
  p[0] = '\0';
  while(temp = strtok(NULL, ".")) {
    strcat(b, p);
    strcpy(p, ".");
    strcat(p, temp);
  }
  return b;
}

/* -------------------------------------------------------------------------- */

double euclidean_distance(Point *p, Atom *a)
{
  /* return the squared euclidean distance between a point and an atom */

  double ed = 0.0; /* the squared euclidean distance */

  ed = (p->x - a->x) * (p->x - a->x) + (p->y - a->y) * (p->y - a->y) +
       (p->z - a->z) * (p->z - a->z);
  return(ed);
}

/* -------------------------------------------------------------------------- */

void set_parameters(int argc, char **argv, Parameters *params)
{
  /* get the command line parameters and set them in the global variables */
  int return_value;
  char *probe;

  /* set default parameters */
  params->resolution = 1.0;
  params->clearance = 5.0;
  params->npoints[X] = params->npoints[Y] = params->npoints[Z] = 0;
  params->center[X] = params->center[Y] = params->center[Z] = 0.0;
  params->set_box = FALSE;
  params->compressed = FALSE;
  params->num_probes = 0;

  /* process the ARGV list */
  while ((argc > 1) && (argv[1][0] == '-')) {
    switch (argv[1][1]) {
    case 'f': /* filename */
      strtok(argv[1], "=");
      strncpy(params->filename, (char *)strtok(NULL, " "), MAX_STR_LEN);
      break;

    case 'p': /* probe */
      strtok(argv[1], "=");
      while ((probe = (char *) strtok(NULL, ","))) {
        strncpy(params->probes[params->num_probes], probe, MAX_STR_LEN);
        params->num_probes++;
      }
      break;

    case 'r': /* resolution */
      strtok(argv[1], "=");
      return_value = sscanf((char *)strtok(NULL, " "), "%lf", &params->resolution);
      if (return_value == 0) {
        fprintf(stderr, "%s option not valid\nAborting...\n", argv[1]);
        fprintf(stderr, "Usage: %s\n", USAGE);
        exit(1);
      }
      break;

    case 'n': /* number of points in X, Y, Z */
      strtok(argv[1], "=");
      return_value = sscanf((char *)strtok(NULL, ","), "%u", &params->npoints[X]);
      if (return_value == 0) {
        fprintf(stderr, "%s option not valid\nAborting...\n", argv[1]);
        fprintf(stderr, "Usage: %s\n", USAGE);
        exit(1);
      }
      return_value = sscanf((char *)strtok(NULL, ","), "%u", &params->npoints[Y]);
      if (return_value == 0) {
        fprintf(stderr, "%s option not valid\nAborting...\n", argv[1]);
        fprintf(stderr, "Usage: %s\n", USAGE);
        exit(1);
      }
      return_value = sscanf((char *)strtok(NULL, ","), "%u", &params->npoints[Z]);
      if (return_value == 0) {
        fprintf(stderr, "%s option not valid\nAborting...\n", argv[1]);
        fprintf(stderr, "Usage: %s\n", USAGE);
        exit(1);
      }
      break;

    case 'c': /* center of the grid */
      strtok(argv[1], "=");
      return_value = sscanf((char *)strtok(NULL, ","), "%lf", &params->center[X]);
      if (return_value == 0) {
        fprintf(stderr, "%s option not valid\nAborting...\n", argv[1]);
        fprintf(stderr, "Usage: %s\n", USAGE);
        exit(1);
      }
      return_value = sscanf((char *)strtok(NULL, ","), "%lf", &params->center[Y]);
      if (return_value == 0) {
        fprintf(stderr, "%s option not valid\nAborting...\n", argv[1]);
        fprintf(stderr, "Usage: %s\n", USAGE);
        exit(1);
      }
      return_value = sscanf((char *)strtok(NULL, ","), "%lf", &params->center[Z]);
      if (return_value == 0) {
        fprintf(stderr, "%s option not valid\nAborting...\n", argv[1]);
        fprintf(stderr, "Usage: %s\n", USAGE);
        exit(1);
      }
      break;

    case 'z': /* use compressed maps */
      params->compressed = TRUE;
      break;
    }
    ++argv;
    --argc;
  }

  /* sanity check (very limited!) */
  if (params->num_probes == 0) {
    fprintf(stderr, "Please enter at least one probe\nAborting...\n");
    fprintf(stderr, "Usage: %s\n", USAGE);
    exit(1);
  }
  if (params->npoints[X] != 0 && params->npoints[Y] != 0 && \
      params->npoints[Z] != 0)
    params->set_box = TRUE;
}

/* -------------------------------------------------------------------------- */

Atom_type *read_atom_types(Parameters params)
{
  /* read an atom_types file and fill in the atom_types structure with
     atom symbol names and charges where available */

  FILE *infile;
  char line[MAX_STR_LEN], *temp;
  Atom_type *first_atom_type = NULL, *current_atom_type;
  
  /* open the input file */
  if (!(infile = fopen(params.atom_types_filename, "r"))) {
    fprintf(stderr, "Cannot open %s\nAborting...\n", 
            params.atom_types_filename);
    exit(1);
  }

  /* parse the infile and fill in the atom_type list */
  while (fgets(line, MAX_STR_LEN, infile) != NULL) {
    current_atom_type = malloc(sizeof(Atom_type));
    temp = (char *)strtok(line, "\t ");
    strncpy(current_atom_type->type, temp, TYPE_LEN);
    strtok(NULL, "\t "); /* discard the atom type description */
    temp = (char *)strtok(NULL, "\t "); /* charge */
    if (temp[strlen(temp) - 1] == '\n')
      temp[strlen(temp) - 1] = '\0';
    if (strcmp(temp, "NA") != 0) {
      sscanf(temp, "%lf", &current_atom_type->charge);
      current_atom_type->has_charge = TRUE;
    }
    else
      current_atom_type->has_charge = FALSE;
    
    current_atom_type->next = first_atom_type;
    first_atom_type = current_atom_type;
  }
  return(first_atom_type);
}

/* -------------------------------------------------------------------------- */

Nb_table *read_nb_table(Parameters params)
{
  /* read in the non-bonded interactions parameters table and store it in
     a linked list  */

  Nb_table *nb_table = NULL, *current_nb;
  FILE *infile;
  char atom_type1[MAX_STR_LEN], atom_type2[MAX_STR_LEN], *temp;
  char line[MAX_STR_LEN];
  double c6, c12; /* lennard-jones coefficients */

  /* open the table file */
  infile = fopen(params.nb_table, "r");
  if (infile == NULL) {
    fprintf(stderr, "Error: unable to open nonbonded table file\n");
    exit(1);
  }

  /* process the table file */
  while (fgets(line, MAX_STR_LEN, infile) != NULL) {
    /* get the atom type 1, atom type 2, c6 and c12 parameters */
    temp = (char *)strtok(line, "\t ");
    strncpy(atom_type1, temp, 20);
    temp = (char *)strtok(NULL, "\t ");
    strcpy(atom_type2, temp);
    temp = (char *)strtok(NULL, "\t ");
    sscanf(temp, "%le", &c6);
    temp = (char *)strtok(NULL, "\t ");
    sscanf(temp, "%le", &c12);

    /* store the pair */
    current_nb = malloc(sizeof(Nb_table));
    strcpy(current_nb->type1, atom_type1);
    strcpy(current_nb->type2, atom_type2);
    current_nb->c6 = c6;
    current_nb->c12 = c12;
    current_nb->next = nb_table;
    nb_table = current_nb;
  }
  fclose(infile);
  return(nb_table);
}

/* -------------------------------------------------------------------------- */

Atom *extract_atoms(Parameters params, Nb_table *nb_list)
{
  /* read an .easymif file and stores the atomic coordinates, the atom types,
     the partial charges and the corresponding L6 and L12 parameters */

  char line[MAX_STR_LEN], temp[MAX_STR_LEN];
  char type[TYPE_LEN + 1]; /* atom type */
  Atom *first_atom = NULL, *current_atom;
  Nb_table *current_nb;
  FILE *infile;
  unsigned int n, i;

  /* open the input file */
  if (!(infile = fopen(params.filename, "r"))) {
    fprintf(stderr, "Cannot open %s\nAborting...\n", params.filename);
    exit(1);
  }

  /* parse the infile and fill in the atom_list */
  while(fgets(line, MAX_STR_LEN, infile) != NULL) {
    if (strncmp(line, "ATOM", 4) == 0 || strncmp(line, "HETATM", 6) == 0) {
      /* allocate memory for the atom */
      current_atom = malloc(sizeof(Atom));

      /* allocate memory for the Lennard-Jones term */
      current_atom->c6 = malloc(sizeof(double) * params.num_probes);
      current_atom->c12 = malloc(sizeof(double) * params.num_probes);

      /* extract the atom type removing blanks */
      strncpy(type, line + PDB_ATOM_SHIFT, 4);
      n = 0;
      while ( n <=4  ) {
        if (type[n] == ' ') {
    	  break;
        }
        n++;
      }
      strncpy(current_atom->type, type, n);
      current_atom->type[n] = '\0';

      /* extract the coordinates */
      strncpy(temp, line + PDB_X_SHIFT, 8);
      temp[8] = '\0';
      sscanf(temp, "%lf", &current_atom->x);
      strncpy(temp, line + PDB_Y_SHIFT, 8);
      temp[8] = '\0';
      sscanf(temp, "%lf", &current_atom->y);
      strncpy(temp, line + PDB_Z_SHIFT, 8);
      temp[8] = '\0';
      sscanf(temp, "%lf", &current_atom->z);

      
      /* extract the partial charge */
      strncpy(temp, line + PDB_OCCUPANCY, 7);
      temp[7] = '\0';
      sscanf(temp, "%lf", &current_atom->partial_charge);

      /* assign the L6 and L12 parameters */
      for (i = 0; i < params.num_probes; i++) { /* for each probe */
        current_nb = nb_list;
        while (current_nb != NULL) {
          if (((strcmp(current_nb->type1, params.probes[i]) == 0) && 
               (strcmp(current_nb->type2, current_atom->type) == 0)) ||
             (strcmp(current_nb->type2, params.probes[i]) == 0 && 
               (strcmp(current_nb->type1, current_atom->type) == 0))) {
            current_atom->c6[i] = current_nb->c6;
            current_atom->c12[i] = current_nb->c12;
            break;
          }
          current_nb = current_nb->next;
        }
      }

      /* allocate the next element of the list */      
      current_atom->next = first_atom;
      first_atom = current_atom;
    }
  }

  fclose(infile);

  return(first_atom); /* return the first element to the linked list */
}

/* -------------------------------------------------------------------------- */

Box calculate_box(Atom *atom_list, Parameters params)
{
  /* compute the center, number of points in x, y, z and the lower corner of
     the enclosing box, given an atom list and a set of user-defined
     parameters (resolution, clearance) */

  Box box;
  Atom *current_atom;
  unsigned int total_num_atoms = 1;
  double lower_coords[3], upper_coords[3]; /* the extreme coords of the atoms */

  /* set the resolution */
  box.resolution = params.resolution;

  /* find the center and extreme points of the atom coordinates */
  current_atom = atom_list;
  box.center[X] = current_atom->x;
  box.center[Y] = current_atom->y;
  box.center[Z] = current_atom->z;
  lower_coords[X] = current_atom->x;
  lower_coords[Y] = current_atom->y;
  lower_coords[Z] = current_atom->z;
  upper_coords[X] = current_atom->x;
  upper_coords[Y] = current_atom->y;
  upper_coords[Z] = current_atom->z;

  current_atom = current_atom->next;
  while (current_atom != NULL) {

    /* add the current coordinates to the center */
    box.center[X] += current_atom->x;
    box.center[Y] += current_atom->y;
    box.center[Z] += current_atom->z;

    /* check if coordinates are current extreme */
    if (lower_coords[X] > current_atom->x)
      lower_coords[X] = current_atom->x;
    if (lower_coords[Y] > current_atom->y)
      lower_coords[Y] = current_atom->y;

    if (lower_coords[Z] > current_atom->z)
      lower_coords[Z] = current_atom->z;

    if (upper_coords[X] < current_atom->x)
      upper_coords[X] = current_atom->x;
    if (upper_coords[Y] < current_atom->y)
      upper_coords[Y] = current_atom->y;
    if (upper_coords[Z] < current_atom->z)
      upper_coords[Z] = current_atom->z;

    total_num_atoms++;
    current_atom = current_atom->next;
  }

  /* calculate the average (center of the molecule = center of the box) */ 
  box.center[X] /= total_num_atoms;
  box.center[Y] /= total_num_atoms;
  box.center[Z] /= total_num_atoms;

  /* calculate the number of points in x, y, z, making sure it's an odd num */
  box.npoints[X] = ceil((2.0 * params.clearance + upper_coords[X] - 
                         lower_coords[X]) / box.resolution);
  if ((box.npoints[X] % 2) == 0)
    box.npoints[X]++; 
  box.npoints[Y] = ceil((2.0 * params.clearance + upper_coords[Y] - 
                         lower_coords[Y]) / box.resolution);
  if ((box.npoints[Y] % 2) == 0)
    box.npoints[Y]++;
  box.npoints[Z] = ceil((2.0 * params.clearance + upper_coords[Z] - 
                         lower_coords[Z]) / box.resolution);
  if ((box.npoints[Z] % 2) == 0)
    box.npoints[Z]++;

  /* calculate the lower corner of the box */
  box.lower_corner[X] = box.center[X] - box.resolution * (box.npoints[X] - 1) / 2.0;
  box.lower_corner[Y] = box.center[Y] - box.resolution * (box.npoints[Y] - 1) / 2.0;
  box.lower_corner[Z] = box.center[Z] - box.resolution * (box.npoints[Z] - 1) / 2.0;

  return(box);
}

/* -------------------------------------------------------------------------- */

Box set_box(Parameters params)
{
  /* set the user-defined parameters in the box structure (this function is
     to be used in alternative to calculate_box, whenever the user wants to
     focus on a specific region on the protein) */

  Box box;

  /* user-defined center of the box */
  box.center[X] = params.center[X];
  box.center[Y] = params.center[Y];
  box.center[Z] = params.center[Z];

  /* user-defined resolution of the box */
  box.resolution = params.resolution;
  
  /* user-defined number of points in x, y, z */
  box.npoints[X] = params.npoints[X];
  box.npoints[Y] = params.npoints[Y];
  box.npoints[Z] = params.npoints[Z];

  /* make sure the number of points is an odd  number */
  if ((box.npoints[X] % 2) == 0)
    box.npoints[X]++;
  if ((box.npoints[Y] % 2) == 0)
    box.npoints[Y]++;
  if ((box.npoints[Z] % 2) == 0)
    box.npoints[Z]++;
  
  /* compute the lower_corner of the box */
  box.lower_corner[X] = box.center[X] - box.resolution * (box.npoints[X] - 1) / 2.0;
  box.lower_corner[Y] = box.center[Y] - box.resolution * (box.npoints[Y] - 1) / 2.0;
  box.lower_corner[Z] = box.center[Z] - box.resolution * (box.npoints[Z] - 1) / 2.0;
  return(box);
}

/* -------------------------------------------------------------------------- */

Point *inizialize_point_list(Box box, Parameters params)
{
  /* create the linked list of points and initialize them */

  Point *first_point = NULL, *current_point;
  unsigned int x, y, z;
  unsigned int i;
  
  for (x = 0; x < box.npoints[X]; x++)
    for (y = 0; y < box.npoints[Y]; y++)
      for (z = 0; z < box.npoints[Z]; z++) {
        current_point = malloc(sizeof(Point));
        for (i = 0; i < params.num_probes; i++) /* inizialize the energy vals */
          current_point->energy[i] = 0.0;
        
        /* initialize the coordinates */
        current_point->x = box.lower_corner[X] + x * box.resolution;
        current_point->y = box.lower_corner[Y] + y * box.resolution;
        current_point->z = box.lower_corner[Z] + z * box.resolution;
        
        /* allocate the next element of the list */      
        current_point->next = first_point;
        first_point = current_point;
      }
  return(first_point);
}

/* -------------------------------------------------------------------------- */

void calculate_MIFs(Point *point_list, Atom *atom_list, Nb_table *nb_list,
                    Atom_type *atom_type_list, Parameters params)
{
  /* carry out the actual calculations for each point and probe */

  Point *current_point;
  Atom *current_atom;
  Atom_type *current_atom_type;
  double rij_squared, rij_Angstrom, rij, rij_6, rij_12; /* pairwise distances */
  double probe_charge = 0.0; /* charge of the current probe */
  double Velectr; 
  unsigned int i;

  for (i = 0; i < params.num_probes; i++) { /* for each probe */
    /* get the charge of the current probe */
    current_atom_type = atom_type_list;
    while (current_atom_type != NULL) {
      if (strcmp(current_atom_type->type, params.probes[i]) == 0) {
        probe_charge = current_atom_type->charge;
        break;
      }
      current_atom_type = current_atom_type->next;
    }

    current_point = point_list;
    while (current_point != NULL) { /* for each point in the box */
      current_atom = atom_list;
      while (current_atom != NULL) { /* for each atom */
        /* calculate the distances */
        rij_squared = euclidean_distance(current_point, current_atom);
        rij_Angstrom = sqrt(rij_squared);
        rij_squared *= 0.01; /* convert from sq Angstrom to sq nm */
        rij = sqrt(rij_squared);
        rij_6 = pow3(rij_squared);
        rij_12 = pow6(rij_squared);

        /* move to the next atom if distance > cutoff */
        if (rij > DISTANCE_CUTOFF) {
          current_atom = current_atom->next; /* go to the next atom */
          continue;
	}
	
        /* Lennard-Jones part */
	current_point->energy[i] += (current_atom->c12[i] / rij_12 - 
                                     current_atom->c6[i] / rij_6);

        /* electrostatics part with reaction field Cui and Mezei 
           distance-dependent dielectric */
        Velectr = fcoul * current_atom->partial_charge * probe_charge / rij;
        if (rij_Angstrom < 1.32)
          Velectr /= 8.0;
        else
          Velectr /= A + B / (1 + kappa * exp(-lambda * B * rij_Angstrom));
	
        current_point->energy[i] += Velectr;
	current_atom = current_atom->next; /* go to the next atom */
      } 
      current_point = current_point->next; /* go to the next point */
    }
  }

}

/* -------------------------------------------------------------------------- */

void compress(FILE *outfile, double *energy_values, unsigned int num_values)
{
  /* implementation of the compression algorithm                              *
   * LZW (Lempel-Ziv-Welch) Welch, IEEE Computer 1984                         
   * with O(1) dictionary search                                              */

  const char dictionary[] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '.', '-', ' '};
  const unsigned int length_dictionary = 13;
  const unsigned int size_dict = 1 << BITS;
  unsigned int i, current_index = 57, string_pos = 0;
  unsigned int current_num_value = num_values - 1;
  const unsigned int num_bytes = BITS / 8;
  int prefix, index;
  Hash *hash_table;
  char string[MAX_STR_LEN], byte;

  /* initialize the hash table (with characters from 0 to 9, point and space) */
  hash_table = malloc(sizeof(Hash) * size_dict);
  for (i = 0; i < length_dictionary; i++) {
    hash_table[(int) dictionary[i]].value = dictionary[i];
    hash_table[(int) dictionary[i]].first = -1;
    hash_table[(int) dictionary[i]].next = -1;
  }

  /* carry out the compression */
  sprintf(string, "%.3f ", energy_values[num_values]);
  prefix = string[0];
  while (current_num_value >= -1 || string[string_pos] != '\0') {
    /* check if floating point number has been fully read */
    string_pos++;
    if (string[string_pos] == '\0') {
      current_num_value--;
      if (current_num_value == -1)
        break;
      string_pos = 0;
      sprintf(string, "%.3f ", energy_values[current_num_value]);
    }
    byte = string[string_pos];

    /* search the current string in the hash table */
    index = hash_table[prefix].first;
    while (index != -1) {
      if (byte == hash_table[index].value) {
        prefix = index;
        break;
      }
      else {
        index = hash_table[index].next;
        }
    }

    /* string not found - add to the hash table */
    if (index == -1) {
      current_index++;
      /* output the prefix code */
      fwrite(&prefix, num_bytes, 1, outfile);
      if (current_index < size_dict) {
        hash_table[current_index].value = byte;
        hash_table[current_index].next = hash_table[prefix].first;
        hash_table[current_index].first = -1;
        hash_table[prefix].first = current_index;
      }
      prefix = byte;
    }
  }
  /* write the last byte */
  fwrite(&index, num_bytes, 1, outfile);

  /* clean up the memory */
  free(hash_table);
}

/* -------------------------------------------------------------------------- */

void create_dx(Parameters params, Point *point_list, Box box)
{
  /* print the maps on dx file(s), assuming the point are already arranged
     according to dx convention (z fast, y medium, x slow) */

  char outfilename[MAX_STR_LEN], *temp, open_mode[MAX_STR_LEN];
  Point *current_point;
  FILE *outfile;
  unsigned int i, total_num_points, counter;
  double *energy_values; /* to store the energy values */

  total_num_points = box.npoints[X] * box.npoints[Y] * box.npoints[Z];
  energy_values = malloc(sizeof(double) * total_num_points);
  for (i = 0; i <  params.num_probes; i++) { /* for every probe */
    /* open the output file */
  
    temp = basename(params.filename);
    strncpy(outfilename, temp, MAX_STR_LEN);
    free(temp);
    strcat(outfilename, "_");
    strcat(outfilename, params.probes[i]);
    if (params.compressed) {
      strcat(outfilename, ".cmp");
    }
    else {
      strcat(outfilename, ".dx");
    }

    /* set the file opening mode */
    if (params.compressed) {
      strncpy(open_mode, "wb", MAX_STR_LEN);
    }
    else {
      strncpy(open_mode, "w", MAX_STR_LEN);
    }
    outfile = fopen(outfilename, open_mode);
    if (outfile == NULL) {
      fprintf(stderr, "Error: unable to write %s file\n", outfilename);
      return;
    }

    /* write the header */
    fprintf(outfile, "# easymifs output\n#\n#\n#\n");
    fprintf(outfile, "object 1 class gridpositions counts %d %d %d\n", 
            box.npoints[X], box.npoints[Y], box.npoints[Z]);
    fprintf(outfile, "origin %.3f %.3f %.3f\n", box.lower_corner[X], 
            box.lower_corner[Y], box.lower_corner[Z]);
    fprintf(outfile, "delta %.3f 0 0\ndelta 0 %.3f 0\ndelta 0 0 %.3f\n",
            box.resolution, box.resolution, box.resolution);
    fprintf(outfile, "object 2 class gridconnections counts %d %d %d\n", 
            box.npoints[X], box.npoints[Y], box.npoints[Z]);
    fprintf(outfile, "object 3 class array type double rank 0 items %d ", 
            total_num_points);
    fprintf(outfile, "data follows\n");

    /* store the energy values for the current probe */
    current_point = point_list;
    counter = 0;
    while (current_point != NULL) {
      if (current_point->energy[i] > ENERGY_CUTOFF)
        energy_values[counter] = ENERGY_CUTOFF;
      else
        energy_values[counter] = current_point->energy[i];
      current_point = current_point->next;

      counter++;
    }

    if (params.compressed) { /* write compressed maps */
      compress(outfile, energy_values, total_num_points);
    }
    else { /* write the points to file in uncompressed form */
      while(counter > 0) {
        counter--;
        fprintf(outfile, "%.3f\n", energy_values[counter]);
      }
    }
    fclose(outfile);
  }
  free(energy_values);
}

/* -------------------------------------------------------------------------- */

void free_memory_atom_types(Atom_type *atom_type_list)
{
  /* free the memory used by the atom_type_list */

  Atom_type *next_atom_type;
  while (atom_type_list != NULL) {
    next_atom_type = atom_type_list->next;
    free(atom_type_list);
    atom_type_list = next_atom_type;
  }
}

/* -------------------------------------------------------------------------- */

void free_memory_nb_table(Nb_table *nb_list)
{
  /* free the memory used by the atom_type_list */

  Nb_table *next_nb_list;

  while (nb_list != NULL) {
    next_nb_list = nb_list->next;
    free(nb_list);
    nb_list = next_nb_list;
  }
}

/* -------------------------------------------------------------------------- */

void free_memory_atoms(Atom *atom_list)
{
  /* free the memory used by the atom_list */

  Atom *next_atom;
  while (atom_list != NULL) {
    next_atom = atom_list->next;
    free(atom_list->c6);
    free(atom_list->c12);
    free(atom_list);
    atom_list = next_atom;
  }
}

/* -------------------------------------------------------------------------- */

void free_memory_points(Point *point_list)
{
  /* free the memory used by the point_list */

  Point *next_point;
  while (point_list != NULL) {
    next_point = point_list->next;
    free(point_list);
    point_list = next_point;
  }
}
