#ifndef _LICENCE_H
#define _LICENCE_H

#define LICENCE_FULL    'f'
#define LICENCE_DEMO    'd'
#define LICENCE_VIEWER  'v'

#define LICENCE_UNIX    'u'
#define LICENCE_WINDOWS 'w'

extern int check_licence(void);
extern int get_licence_type(void);
extern int get_licence_os(void);
extern char *get_licence_id(void);
extern int get_licence_expire(void);
extern int get_licence_users(void);
extern int valid_seq(char *seq, int len);
extern void viewer_mode(void);

#endif /* _LICENCE_H */
