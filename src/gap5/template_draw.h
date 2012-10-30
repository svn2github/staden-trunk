/*
 * template_draw.h - fast line drawing for the template display
 *
 * Andrew Whitwham, January 2010
 * Wellcome Trust Sanger Institute
 *
 */

#ifndef _TEMPLATE_DRAW_H_
#define _TEMPLATE_DRAW_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <inttypes.h>

#include <tcl.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>

typedef struct {
    void *buf;  // image buffer
    int height; // pixel (Raster) size
    int width;	// pixel (Raster) size
    Display *dis;
    int screen;
    int depth;
    void *map_col;     // colour array
    int colour_count;  
    int colour_size;
    int single_col;    // special colours below
    int span_col[10];
    int inconsistent_col;
    int fwd_col;
    int rev_col;
    int fwd_col3;
    int rev_col3;
    int black;
    XImage *img;
} image_t;

image_t *initialise_image(Display *dis);
int add_colour(image_t *image, unsigned int red, unsigned int green, unsigned int blue);
int create_image_buffer(image_t *image, int width, int height, int bg_colour);
int draw_line(image_t *image, int x1, int x2, int y, int colour);
void create_image_from_buffer(image_t *image);
void image_remove(image_t *image);
void image_destroy(image_t *image);

#endif

