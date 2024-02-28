#ifndef STDLIB_LINALG_STRIDE_C_H_INCLUDED
#define STDLIB_LINALG_STRIDE_C_H_INCLUDED

#include <ISO_Fortran_binding.h>
#include <stdlib.h>

// Definitions matching the interface in stdlib_linalg_stride
typedef void*    CFI_ADDRESS;
typedef int8_t   CFI_RANK;
typedef int      CFI_FLAG;
typedef int16_t  CFI_TYPE;
typedef intptr_t CFI_SIZE;
typedef int8_t   CFI_ATTR;

#define MAX_FORTRAN_RANK 15

// Compiler independent type flags
enum CFI_TYPES {
   INTEGER = 1,
   LOGICAL,
   REAL,
   COMPLEX,
   CHARACTER,
   STRUCT,
   CPTR,
   CFUNPTR,
   OTHER,
   SIGNED_CHAR,
   SHORT,
   INT,
   LONG,
   LONG_LONG,
   SIZE_T,
   INT8_T,
   INT16_T,
   INT32_T,
   INT64_T,
   INT_LEAST8_T,
   INT_LEAST16_T,
   INT_LEAST32_T,
   INT_LEAST64_T,
   INT_FAST8_T,
   INT_FAST16_T,
   INT_FAST32_T,
   INT_FAST64_T,
   INTMAX_T,
   INTPTR_T,
   PTRDIFF_T,
   BOOL,
   FLOAT,
   DOUBLE,
   FLOAT_COMPLEX,
   DOUBLE_COMPLEX
};

// C Descriptor: dimension information
typedef struct array_dimension {
    CFI_SIZE lower_bound;
    CFI_SIZE extent;
    CFI_SIZE stride_bytes;
} array_dimension;

// Global C Array Descriptor
typedef struct array_descriptor {
    CFI_ADDRESS base_address;
    CFI_SIZE elem_bytes;
    CFI_FLAG version;
    CFI_RANK rank;
    CFI_TYPE type;
    CFI_ATTR attribute;
    array_dimension dim[MAX_FORTRAN_RANK];
} array_descriptor;


#endif // STDLIB_LINALG_STRIDE_C_H_INCLUDED
