/* STRINGS.H: a quirck for SYS5 type systems
 *
 *	This header file must only be added to $NEMO/inc
 *	when the local host does only support <string.h>,
 *	i.e. the SYS5 convention, plus adds str(r)chr to
 *	the function list by a #define
 *
 */

#include <string.h>

#define index  strchr
#define rindex strrchr

