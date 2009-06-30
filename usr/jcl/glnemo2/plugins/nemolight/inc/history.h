/*
 *      history.h:
 */

#ifndef _history_h
#define _history_h
#ifndef  HistoryTag
#define HistoryTag "History"            /* used to tag history in input     */
#endif

#ifdef __cplusplus
extern "C" {
#endif	
int    get_history   (stream);
int    put_history   (stream);
int    app_history   (string);
void   reset_history (void);
void   set_headline  (string);
string ask_headline  (void);
string *ask_history  (void);
#ifdef __cplusplus
}
#endif	

#endif /* _history_h */

