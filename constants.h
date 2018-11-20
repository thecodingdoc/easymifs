/******************************************************************************
 ******************************************************************************
 ** constants.h  (part of the easymifs program)                              **
 ** Author:  Dario Ghersi                                                    **
 ** Version: 20160418                                                        **
 ******************************************************************************
 ******************************************************************************/

#define USAGE "easyMIFs -f=INFILE -p=PROBE [-r=RESOLUTION -c=X,Y,Z -n=DX,DY,DZ -z]"

#define MAX_STR_LEN 1000
#define MAX_NUM_PROBES 10
#define TYPE_LEN 5 /* lenght of the atom type field */

#define ENERGY_CUTOFF 5.0 /* the capping value for energy */
#define DISTANCE_CUTOFF 1.4 /* in nm, not Angstrom! */
#define X 0
#define Y 1
#define Z 2
#define BITS 16 /* size of the code for compression (byte multiple) */

typedef enum { FALSE, TRUE } bool;

/*****************************************************************************
 * GROMOS FF CONSTANTS                                                       *
 *****************************************************************************/
#define fcoul  138.935485 /* 1 / 4_pi_eps0 */
#define epsilon_r 1.0
#define epsilon_rf 78.0
#define r_c 0.9
#define r_c_cube 0.729

/*****************************************************************************
 * DISTANCE-DEPENDENT DIELECTRIC FUNCTION TAKEN FROM the paper:              * 
 * Cui M et al, Protein Engineering, Design & Selection, 2008                *
 *****************************************************************************/
#define lambda 0.018733345
#define kappa 213.5782
#define e0 78.4
#define A 6.02944
#define B 72.37056 /* e0 - A */

/*****************************************************************************
 * PDB FORMAT CONSTANTS                                                      *
 *****************************************************************************/
#define PDB_ATOM_SHIFT 12
#define PDB_X_SHIFT 30
#define PDB_Y_SHIFT 38
#define PDB_Z_SHIFT 46
#define PDB_OCCUPANCY 54

/*****************************************************************************
 * MACROS                                                                    *
 *****************************************************************************/
#define pow3(x) ((x)*(x)*(x))
#define pow6(x) ((x)*(x)*(x)*(x)*(x)*(x))
