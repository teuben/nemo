#ifndef	NEMO_EXCEP
#define NEMO_EXCEP

extern int BeginBlock (void);
extern int EndBlock (void);

extern int RegisterPointer (void * ptr);
extern int UnRegisterPointer (void * ptr);

extern int RegisterStream (stream fptr);
extern int UnRegisterStream (stream fptr);

extern void RaiseException (int errNumber);
extern void RestoreUserContext (void);

// extern void FreeAllNemoResources (void);

#endif /* NEMO_EXCEP */
