/* 
  Program Name: mess
  File name: mess.c
  Purpose: put user information to the screen for the text version
         of osp--maybe error messages or any kind of user information
	 or question.

  Last Update: Tuesday 13 August 1991

  Copyright 1991: LaDeana Hillier and Philip Green

  Change Log:

  Modified to work with xdap output window
*/

/* --- includes --- */
#include <stdio.h>
#include <stdarg.h>

#include "text_output.h"

/* ---- Exports ---- */
void messagef(char *format, ...)
{
    va_list args;
    va_start (args,format);
    vfprintf (stdout, format, args);
    va_end(args);
    UpdateTextOutput();
}


