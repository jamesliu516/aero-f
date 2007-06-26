#ifndef _BC_DEF_H_
#define _BC_DEF_H_

#define BC_SYMMETRY 6
#define BC_OUTLET_MOVING -5
#define BC_OUTLET_FIXED 5
#define BC_INLET_MOVING -4
#define BC_INLET_FIXED 4
#define BC_ADIABATIC_WALL_MOVING -3
#define BC_ADIABATIC_WALL_FIXED 3
#define BC_SLIP_WALL_MOVING -2
#define BC_SLIP_WALL_FIXED 2
#define BC_ISOTHERMAL_WALL_MOVING -1
#define BC_ISOTHERMAL_WALL_FIXED 1
#define BC_INTERNAL 0

#define BC_MIN_CODE -5
#define BC_MAX_CODE 6

#define BC_FREE        0
#define BC_FIXED       1
#define BC_MATCHED     2
#define BC_CONSTRAINED 3 
#endif
