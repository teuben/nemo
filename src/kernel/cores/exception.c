/* exceptional.c -- will store source for 
 * exception handlers eventually
 * For now, stores code being tested/developed
 * Dec 10 2001 - NSA
 */

#include <stdinc.h>
#include <errno.h>
#include <setjmp.h>
#include <signal.h>
#include <assert.h>

#include "exception.h"


#define	InBlock()	((inNemo) > 0)

#define SUCCESS		0


typedef struct memBlkElem {
    void*		ptr;
    struct memBlkElem*	next;
} MemElem;

typedef struct fileElem {
    stream		str;
    struct fileElem*	next;
} StreamElem;

typedef struct excepBlkElem {
//   char*			blkName;
    sigjmp_buf			context;
    MemElem*			headMem;
    StreamElem*			headStream;
//    MessageFunction		dosomething;
//	(void *)(*printer_handler)(...);
    struct excepBlkElem*	next;

} BlockRes;

static BlockRes		topLevelRes = { {}, NULL, NULL, NULL};
static BlockRes*	firstRes = &topLevelRes;

static int		inNemo = 0; // shows level of nesting for Nemo
				    // routines; 0 if we are not in a 
				    // Nemo routine
static BlockRes*
PushBlock (BlockRes* blkResources) // function will do more if full scale
{ 				   // stacking is implemented
    inNemo++;
    assert (inNemo > 0);
    return blkResources;
}


static BlockRes*
GetCurrentBlock (void)
{
    return firstRes;
}



static BlockRes*
RegisterBlock (void)
{
    if (!InBlock()) {
	firstRes->headMem =  NULL; 
	firstRes->headStream = NULL;
	firstRes->next = NULL;
    }

    return PushBlock (firstRes);
}


// Used by caller to indicate beginning of execution in Nemo space
// If called in nested fashion, has no effect; in other words, you're
// either in Nemo space or not; no nested Nemo "levels"
int
BeginBlock (void)
{
    BlockRes*	blkResources = NULL;

    blkResources = RegisterBlock();
    assert (blkResources == firstRes);
    if (sigsetjmp (blkResources->context, 1) != 0)  {
	return -1;
    }
	
    return 0;
}

// Free all memory allocated in this block
static void
FreeAllBlkMem (MemElem* blkMemPtr)
{
    MemElem*	memResPtr = blkMemPtr;
    MemElem*	prevPtr = NULL;

    for (prevPtr = memResPtr;
		memResPtr != NULL && memResPtr->next != NULL;
		prevPtr = memResPtr, memResPtr = memResPtr->next, free(prevPtr))
	    if (memResPtr->ptr != NULL)  {
		free (memResPtr->ptr);
		memResPtr->ptr = NULL;
	    }
}


// Close all file streams opened in this block
static void
FreeAllBlkStreams (StreamElem* blkStreamPtr)
{
    StreamElem*	fileResPtr = blkStreamPtr;
    StreamElem*	prevPtr = NULL;

    for (prevPtr = fileResPtr;
		fileResPtr != NULL && fileResPtr->next != NULL;
		prevPtr = fileResPtr, fileResPtr = fileResPtr->next,
			free(prevPtr))
	    if (fileResPtr->str != NULL)  {
		fclose (fileResPtr->str);
		fileResPtr->str = NULL;
	    }
}


static int
FreeResBlock (BlockRes* resPtr)
{
    FreeAllBlkMem (resPtr->headMem);
    FreeAllBlkStreams (resPtr->headStream);
    // No need to free sigjmp_buf as it is a const
    inNemo = (InBlock()) ? --inNemo : inNemo;
    assert (inNemo >= 0);
}


static void
UnwindBlock (BlockRes* blkResources)
{
    FreeResBlock (blkResources);
    blkResources->headMem = NULL;
    blkResources->headStream = NULL;
    assert (blkResources->next == NULL);
}


// Called when a "clean" return is desired by the caller; assume that all 
// open files and allocated memory are still required, so dont free them.
// Instead transfer them to the next higher block if there is one. 
// QQQ What do we do if there isnt one ? QQQ AAA There should always be one !
static void
UnRegisterBlock (void)
{
    BlockRes*		blkResources = NULL;

    if (! InBlock())
	return;

    blkResources = GetCurrentBlock();
    inNemo = (InBlock()) ? --inNemo : inNemo;
    assert (inNemo >= 0);
}


static MemElem*
PushPointer (MemElem* head, void* ptr)
{
    MemElem*	memElem = NULL;

    if ((memElem = (MemElem *) calloc (1, (sizeof(MemElem)))) == NULL)
	return NULL;

    memElem->ptr = ptr;
    memElem->next = head;
    head = memElem;    // attach memElem at head of list

    return head;
}



int
RegisterPointer (void * ptr)
{
    BlockRes*	blkRes = NULL;

    if (! InBlock())
	return SUCCESS;

    blkRes = GetCurrentBlock();

    // Todo - error check here
    blkRes->headMem = PushPointer (blkRes->headMem, ptr);

    return SUCCESS;
}


static StreamElem*
PushStream (StreamElem* head, stream str)
{
    StreamElem*	streamElem = NULL;

    if ((streamElem = (StreamElem *)calloc(1, (sizeof(StreamElem)))) == NULL)
	return NULL;

    streamElem->str = str;
    streamElem->next = head;
    head = streamElem;    // attach streamElem at head of list

    return head;
}


int
RegisterStream (stream fptr)
{
    BlockRes*	blkRes = NULL;

    if (! InBlock())
	return SUCCESS;

    blkRes = GetCurrentBlock();
    blkRes->headStream = PushStream (blkRes->headStream, fptr);

    return SUCCESS;
}

static MemElem*
PopMem (MemElem* head, void* ptr)
{
    MemElem*	prevPtr = NULL;
    MemElem*	temp = NULL;

    for (temp = prevPtr = head; temp != NULL && temp->ptr != ptr;
	    prevPtr = temp, temp = temp->next)
	;

    if (temp == NULL)
	return NULL;

    prevPtr->next = temp->next;

    free (temp);
    return (prevPtr == head) ? prevPtr : head;
}


// Routine is used when memory resource is freed 
// normally by user
int
UnRegisterPointer (void * ptr)
{
    BlockRes*	blkRes = NULL;

    if (ptr == NULL || ! InBlock())
	return SUCCESS;

    blkRes = GetCurrentBlock();

    blkRes->headMem = PopMem (blkRes->headMem, ptr);
    return SUCCESS;
}

static StreamElem*
PopStream (StreamElem* head, stream fptr)
{
    StreamElem*	prevPtr = NULL;
    StreamElem*	temp    = NULL;

    for (temp = prevPtr = head; temp != NULL && temp->str != fptr;
	    prevPtr = temp->next, temp = temp->next)
	;

    if (temp == NULL)
	return NULL;

    prevPtr->next = temp->next;

    free (temp);
    return (prevPtr == head) ? prevPtr : head;
}

// Routine is used when file resource is freed 
// normally by user
int
UnRegisterStream (stream fptr)
{
    BlockRes*	blkRes = NULL;

    if (fptr == NULL || ! InBlock())
	return SUCCESS;

    blkRes = GetCurrentBlock();

    blkRes->headStream = PopStream (blkRes->headStream, fptr);
    return SUCCESS;
}


int
EndBlock (void)
{
    BlockRes*	blkRes = NULL;

    inNemo--;
    assert (inNemo >= 0);
    if (InBlock())
	return SUCCESS;

    blkRes = GetCurrentBlock();
    UnwindBlock (blkRes);
    return SUCCESS;
}



static int	nemoError	= 0;

static sigset_t	sigMask;

static struct sigaction		oldAction;
static struct sigaction		nemoErrorAction;


static void
ExceptionHandler (int signal)
{
    sigjmp_buf	context;
    BlockRes*	blkRes = NULL;

    // We must be in Nemo space
    assert (InBlock());

    blkRes = GetCurrentBlock();
    assert (blkRes == firstRes);

    UnwindBlock (blkRes);
    assert (firstRes->next == NULL);

    siglongjmp (blkRes->context, -1);
    assert (FALSE);		// shouldnt be here
}


void
RaiseException (int errNumber)
{
    // We must be in Nemo space
    assert (InBlock());

    nemoError = errNumber;	// save system error; restore for processing
    				// by Nemo user
#if 1
    nemoErrorAction.sa_handler = ExceptionHandler;
    sigfillset (&sigMask);
    nemoErrorAction.sa_mask = sigMask;
    nemoErrorAction.sa_flags = SA_ONESHOT;
    nemoErrorAction.sa_restorer = NULL;
    if (sigaction (SIGUSR1, &nemoErrorAction, &oldAction) == -1)
	assert (FALSE);
    raise (SIGUSR1);

    // never returns here
    assert (FALSE);
#else
	error("RaiseException: can't get here....");
#endif

}


void
RestoreUserContext(void)
{
    sigaction (SIGUSR1, &oldAction, NULL);
    errno = nemoError;
}


