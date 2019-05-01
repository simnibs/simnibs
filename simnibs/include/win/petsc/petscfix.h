#if !defined(INCLUDED_PETSCFIX_H)
#define INCLUDED_PETSCFIX_H

typedef int gid_t;
typedef int uid_t;
typedef int int32_t;
typedef int mode_t;
typedef int pid_t;
#if defined(__cplusplus)
extern "C" {
#include <stddef.h>
int getdomainname(char *, size_t);
double drand48(void);
void   srand48(long int);
}
#else
#include <stddef.h>
int getdomainname(char *, size_t);
double drand48(void);
void   srand48(long int);
#endif
#endif
