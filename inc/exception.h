#ifndef	NEMO_EXCEP
#define NEMO_EXCEP

#include <stdinc.h>

// marks point of return from exception handler
extern int BeginBlock (void);

// turns off exception handling for the block begun
// with the last BeginBlock; does cleanup of memory and streams
extern int EndBlock (void);

// Marks allocated memory as resource to be freed when
// error occurs
extern int RegisterPointer (void * ptr);

// Unmarks resource - for example, when it is freed
extern int UnRegisterPointer (void * ptr);

extern int RegisterStream (stream fptr);
extern int UnRegisterStream (stream fptr);

// Raises exception when irretrievable error occurs
// Transfers control to exception handling mechanism
extern void RaiseException (int errNumber);

// Restores context after exception handling is done
extern void RestoreUserContext (void);

// extern void FreeAllNemoResources (void);

#endif /* NEMO_EXCEP */
