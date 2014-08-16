//
//  deformerConst.h
//

#pragma once

#ifndef ProbeDeformer_deformerConst_h
#define ProbeDeformer_deformerConst_h

// parametrisation mode
#define BM_SRL 0
#define BM_SES 1
#define BM_LOG3 3
#define BM_LOG4 4
#define BM_QSL 5
#define BM_AFF 10
#define BM_OFF -1

// weight mode
#define WM_INV_DISTANCE 0
#define WM_CUTOFF_DISTANCE 1
#define WM_DRAW 2
#define WM_HARMONIC 3
#define WM_HARMONIC_NEIBOUR 4

// constraint mode
#define CONSTRAINT_NEIBOUR 0
#define CONSTRAINT_CLOSEST 1

// stiffness mode
#define SM_NONE 0
#define SM_PAINT 1
#define SM_LEARN 2

// visualisation mode
#define VM_OFF 0
#define VM_ENERGY 1
#define VM_EFFECT 2
#define VM_CONSTRAINT 3
#define VM_STIFFNESS 4

// error codes
#define NUMPRB_AND_ATTR_DIFFERENT 2
#define INCOMPATIBLE_MESH 3
#define ERROR_ATTR 4


#endif
