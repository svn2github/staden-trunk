/* Copyright Genome Research Limited (GRL). All rights reserved */

#ifndef _FINISH_FILTER_H
#define _FINISH_FILTER_H

#include "finish.h"

#define ctrl(a) ((a) & 0x1f)
#define FILTER_LOWCOMPLEX	'#'

#if 0
#define FILTER_POLYA		ctrl('A')
#define FILTER_POLYC		ctrl('C')
#define FILTER_POLYG		ctrl('G')
#define FILTER_POLYT		ctrl('T')
#define FILTER_POLYK		ctrl('K') /* GT */
#define FILTER_POLYM		ctrl('M') /* AC */
#define FILTER_POLYR		ctrl('R') /* AG */
#define FILTER_POLYS		ctrl('S') /* CG */
#define FILTER_POLYW		ctrl('W') /* AT */
#define FILTER_POLYY		ctrl('Y') /* CT */
#endif

#define FILTER_POLYA		'0'
#define FILTER_POLYC		'1'
#define FILTER_POLYG		'2'
#define FILTER_POLYT		'3'
#define FILTER_POLYK		'4' /* GT */
#define FILTER_POLYM		'5' /* AC */
#define FILTER_POLYR		'6' /* AG */
#define FILTER_POLYS		'7' /* CG */
#define FILTER_POLYW		'8' /* AT */
#define FILTER_POLYY		'9' /* CT */

#define is_filtered(c) ((c) == '#' || ((c) >=  FILTER_POLYA && (c) <= FILTER_POLYY))

void finish_filter(finish_t *fin, char *seq, int len);

#endif /*_FINISH_FILTER_H*/
