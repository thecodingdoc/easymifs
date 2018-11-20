/******************************************************************************
 ******************************************************************************
 ** easymifs.c                                                               **
 ** Author:  Dario Ghersi                                                    **
 ** Version: 011210                                                          **
 ** Goal:    The program aims to become a flexible platform to calculate MIFs**
 **          (Molecular Interaction Fields) using a variety of Force Fields  **
 **          Currently the program supports the GROMOS force field           **
 **                                                                          **
 ** Revisions:                                                               **
 ** 011210   Added the option for writing compressed maps                    ** 
 ******************************************************************************
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "prototypes.h"

/******************************************************************************
 * MAIN PROGRAM                                                               *
 ******************************************************************************/

int main(int argc, char **argv) 
{
  Box box; /* the box surrounding the molecule where the MIFs calculations are
            performed */
  Atom_type *atom_type_list; /* the linked list of the atom types */
  Atom *atom_list; /* the linked list of the atoms making up the molecule */
  Nb_table *nb_list; /* the linked list of non-bonded parameters */
  Point *point_list; /* the linked list of the MIF map(s) */
  Parameters params; /* used defined parameters */
  Point *current_point;

  /* print a welcome message */
  fprintf(stdout, "\n***********************************************\n");
  fprintf(stdout, "* Welcome to easyMIFs, version: %s        *\n", VERSION);
  fprintf(stdout, "***********************************************\n\n");

  /* process the command-line parameters */
  set_parameters(argc, argv, &params);
  strcpy(params.atom_types_filename, "atom_types.txt");
  strcpy(params.nb_table, "ffG43b1nb.params");

  /* read in the atom types table */
  atom_type_list = read_atom_types(params);
  
  /* read the non-bonded table information and store the values */
  nb_list = read_nb_table(params);

  /* extract the atom information from the .easymifs file */
  atom_list = extract_atoms(params, nb_list);

  if (params.set_box) /* use user-defined parameters */
    box = set_box(params);
  else /* calculate the corresponding box */
    box = calculate_box(atom_list, params);
  
  /* initialize the point_list */
  point_list = inizialize_point_list(box, params);
  current_point = point_list;

  /* calculate the MIF(s) */
  calculate_MIFs(point_list, atom_list, nb_list, atom_type_list, params);

  /* print the dx file */
  create_dx(params, point_list, box);

  /* free the memory */
  free_memory_atom_types(atom_type_list);
  free_memory_nb_table(nb_list);
  free_memory_atoms(atom_list);
  free_memory_points(point_list);

  return 0;
}
