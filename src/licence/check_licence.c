#include "licence.h"

int get_licence_type(void) { return 'f'; }
int get_licence_os(void)   { return 'u'; }
char *get_licence_id(void) { return "no-licence-needed"; }
int get_licence_users(void) { return 0; }
int get_licence_expire(void) { return 0; }
int check_licence(void)    { return 0; }
void viewer_mode(void) {}

