// Get stride/contiguous attribute from an array using ISO_Fortran_binding.h
#include <ISO_Fortran_binding.h>
#include <stdio.h>
#include <stdlib.h>
#include "stdlib_linalg_stride_c.h"

// Initialize empty CFI dim info
array_dimension new_CFI_dim_info()
{
    array_dimension CFI;
    CFI.lower_bound = 0;
    CFI.extent = 0;
    CFI.stride_bytes = 0;
    return CFI;
}

// Initialize empty CFI descriptor
array_descriptor new_array_descriptor()
{
    array_descriptor CFI;
    CFI.base_address = NULL;
    CFI.elem_bytes   = 0;
    CFI.version      = 0;
    CFI.rank         = 0;
    CFI.type         = 0;
    CFI.attribute    = 0;
    for (CFI_RANK i=0; i<MAX_FORTRAN_RANK; i++)
    {
        CFI.dim[i] = new_CFI_dim_info();
    }
    return CFI;
}

// Read whole CFI descriptor of a variable; pass it back to Fortran
array_descriptor CFI_to_Fortran(const CFI_cdesc_t * descr)
{
   // Initialize descriptor
   array_descriptor fdesc = new_array_descriptor();

   if (descr) {

        fdesc.base_address = (CFI_ADDRESS) descr->base_addr;
        fdesc.elem_bytes = (CFI_SIZE) descr->elem_len;
        fdesc.version = (CFI_FLAG) descr->version;
        fdesc.rank = (CFI_RANK) descr->rank;
        fdesc.type = (CFI_TYPE) descr->type;
        fdesc.attribute = (CFI_ATTR) descr->attribute;

        for (CFI_RANK rank=0; rank<descr->rank; rank++)
        {
           fdesc.dim[rank].lower_bound  = descr->dim[rank].lower_bound;
           fdesc.dim[rank].extent       = descr->dim[rank].extent;
           fdesc.dim[rank].stride_bytes = descr->dim[rank].sm;
        }

   }

   return fdesc;
}


