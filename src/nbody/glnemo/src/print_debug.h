// ============================================================================
//
//
// THis file allow to activate debug print statement on screen
//


#define ALL_DEBUG 0  // set to 1 if you want to debug all files

#if (ALL_DEBUG || LOCAL_DEBUG)
#  define PRINT_D if (1)
#else
#  define PRINT_D if (0)
#endif

// ============================================================================
