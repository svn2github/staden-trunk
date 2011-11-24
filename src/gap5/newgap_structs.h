#ifndef _NEWGAP_STRUCTS_H
#define _NEWGAP_STRUCTS_H
#include "tg_gio.h"

typedef struct {
    GapIO *io;
} io_arg;

typedef struct {
    int handle;
} handle_arg;

typedef struct {
    GapIO *io;
    tg_rec contig;
} contig_arg;

typedef struct {
    GapIO *io;
    char *inlist;
} list2_arg;

typedef struct {
    char *db_name;
    char *version;
    char *access_str;
    int create;
} open_db_arg;

typedef struct {
    GapIO *io;
    char *version;
    int collect;
} copy_db_arg;

typedef struct {
    GapIO *io;
    char *contig;
    int max_ht;
    int tick_dist;
    char *win_name;
    int line_width;
    int template_id;
} display_item_arg;

typedef struct {
     GapIO *io;
     int id;
     char *win_ruler;
     int disp_ruler;
     int line_width;
     char *colour;
     int offset;
     int tick_height;
     int tick_width;
     char *tick_colour;
     int text_offset;
     int tag_offset;
     int tag_width;
} ruler_arg;

typedef struct {
    GapIO *io;
    int id;
    int ticks;
} r_ticks_arg;

typedef struct {
    GapIO *io;
    int id;
    char *contigs;
    int x;
} update_order_arg;


typedef struct {
    GapIO *io;
    char *contigs;
    int ordered;
} show_relationships_arg;

typedef struct {
    GapIO *io;
    char *contigs;
    int avg_len;
} long_gels_arg;

typedef struct {
    GapIO *io;
    char *contigs;
    int avg_len;
} taq_terms_arg;

typedef struct {
    GapIO *io;
    int idir;
    int minmat;
    char *inlist;
    char *outfile;
    char *tag_list;
} find_repeats_arg;

typedef struct {
    GapIO *io;
    char *inlist;
    char *mode;
    int end_size;
    int min_map_qual;
    int min_freq;
    char *libraries;
} readpair_arg;

typedef struct {
    GapIO *io;
    char *list;
    int move;
    int remove_holes;
    int duplicate_tags;
} dis_reading_arg;

typedef struct {
    GapIO *io;
    tg_rec contig;
    int pos;
} contig_pos_arg;

typedef struct {
    GapIO *io;
    tg_rec contig1;
    tg_rec contig2;
    int pos1;
    int pos2;
} contig_pos2_arg;

typedef struct {
    GapIO *io;
    char *mask;
    char *mode;
    int min_overlap;
    float max_mis;
    int word_len;
    float max_prob;
    int min_match;
    int band;
    int win_size;
    int dash;
    int min_conf;
    int use_conf;
    int use_hidden;
    int max_display;
    int fast_mode;
    float filter_words;
    char *tag_list;
    char *inlist;
} fij_arg;

typedef struct {
    int handle;
    char *inlist;
    int disp_mode;
    int min_mat;
    int max_pad;
    float max_mis;
    int align;
    int joins;
    int fail_mode;
    int win_size;
    int dash;
    char *tag_list;
    int ignore_prev;
    int min_ovr;
} aa_arg;

typedef struct {
    GapIO *io;
    char *list;
    char *dir;
    int format;
} extract_arg;

typedef struct {
    int handle;
    char *list;
} pre_ass_arg;

typedef struct {
    GapIO *io;
    char *contig;
} delete_contig_arg;
 
typedef struct {
    GapIO *io;
    int anno;
} ann_addr_arg;

typedef struct {
    GapIO *io;
    char *list;
    int maxmis;
    float maxperc;
} dstrand_arg;

typedef struct {
    GapIO *io;
    char *inlist;
    int search_from;
    int search_to;
    int num_primers;
    int primer_start;
    char *primer_defs;
} primer_arg;

typedef struct {
    GapIO *io;
    char *inlist;
    char *type;
    char *mask;
    int win_size;
    int dash;
    int format;
    int gel_anno;
    int truncate;
    int gel_notes;
    char *out_file;
    char *tag_list;
    int nopads;
    int min_conf;
    int use_conf;
    int name_format;
} consensus_arg;
 
typedef struct {
    GapIO *io;
    char *contig;
    char *frame;
    char *quality;
    int cursor_wd;
    char *cursor_fill;
} disp_quality_arg;

typedef struct {
    GapIO *io;
    char *contig;
    char *frame;
    char *quality;
    int id;
} quality_arg;

typedef struct {
    GapIO *io;
    int id;
    int t_num;
} template_read_arg;

typedef struct {
    GapIO *io;
    char *contig;
    char *frame;
    char *tag_list;
    char *win_template;
    char *win_ruler;
    int line_width;
    int line_bold;
    int cursor_wd;
    char *cursor_fill;
} disp_template_arg;

typedef struct {
    GapIO *io;
    int id;
    int recalc;
} template_arg;

typedef struct {
    GapIO *io;
    int id;
    int readings;
    int ruler;
    char *tag_list;
} r_tags_arg;

typedef struct {
    GapIO *io;
    char *file;
    int unpadded;
} enter_tags_arg;

typedef struct {
    GapIO *io;
    int id;
} cs_tags_arg;

typedef struct {
    char *filename;
} read_enz_arg;

typedef struct {
    GapIO *io;
    char *frame;
    char *names;
    char *plot;
    int strand;
    char *contigs;
    int tick_ht;
    int tick_wd;
    char *tick_fill;
    int cursor_wd;
    char *cursor_fill;
    int yoffset;
} stop_codon_arg;

typedef struct {
    GapIO *io;
    int id;
    tg_rec contig;
    int strand;
    int update;
} refresh_codon_arg;

typedef struct {
    GapIO *io;
    char *inlist;
    int display;
    float mism;
    int align;
    int enter_failures;
    int ignore_vec;
} ass_direct_arg;

typedef struct {
    GapIO *io;
    char *inlist;
    int win_size;
    int ignore_N;
    float max_mismatch;
} check_ass_arg;

typedef struct {
    GapIO *io;
    char *type;
} anno_list_arg;

typedef struct {
    GapIO *io;
    char *annos;
} delete_anno_list_arg;

typedef struct {
    GapIO *io;
    int min_size;
    int max_size;
    float max_perc;
    int from;
    int to;
    char *vectors;
    char *contigs;
    char *primer_arg;
} find_probes_arg;

typedef struct {
    GapIO *io;
    char *inlist;
    float mis_match;
    char *tag_list;
    char *seq;
    int   consensus_only;
    int   cutoffs;
    char *file;
} oligo_arg;

typedef struct {
    GapIO *io;
    char *tag_list;
    int unpadded;
} add_tags_arg;

typedef struct {
    GapIO *io;
    int order;
} ord2num_arg;

typedef struct {
    GapIO *io;
    int id;
    int r_id;
    float amount;
    int x1;
    int y1;
    int x2;
    int y2;
    char *scroll;
} zoom_arg;

typedef struct {
    GapIO *io;
    int id;
    char *xscroll;
    char *yscroll;
} scroll_arg;

typedef struct {
    GapIO *io;
    int id;
} resize_arg;

typedef struct {
    GapIO *io;
    int id;
    tg_rec contig;
    int cx;
} l_cursor_arg;

typedef struct {
    GapIO *io;
    int id;
    char *window;
} delete_arg;

typedef struct {
    GapIO *io;
    char *window;
    char *frame;
    int tick_ht;
    int tick_wd;
    char *tick_fill;
    int tag_wd;
    int tag_offset;
    int cursor_wd;
    char *cursor_fill;
} display_cs_arg;

typedef struct {
    GapIO *io;
    int id;
    char *window;
    char *v_window;
} display_cp_arg;

typedef struct {
    GapIO *io;
    int id;
    int cx;
    char *contig;
} t_order_arg;

typedef struct {
    GapIO *io;
    char *contig;
} duplicates_arg;

typedef struct {
      GapIO *io;
      char *contig;
      int quality;
} qclip_arg;

typedef struct {
    GapIO *io;
    char *inlist;
    int tag;
} dclip_arg;

typedef struct {
    GapIO *io;
    char *inlist;
} reads_in_contig_arg;

typedef struct {
    GapIO *io;
    char *inlist;
    int summary;
} list_conf_arg;

typedef struct {
    GapIO *io;
    char *contigs;
    char *tag_list;
} find_tags_arg;

typedef struct {
    GapIO *io;
    char *contigs;
    char *frame;
    char *r_win;
    int cursor_wd;
    char *cursor_fill;
} consistency_arg;

typedef struct {
    GapIO *io;
    int id;
    char *frame;
    char *conf_win;
    char *r_win;
    int two_alleles;
} confidence_arg;

typedef struct {
    GapIO *io;
    int id;
    char *frame;
    char *conf_win;
    char *r_win;
    int strand;
} reading_cov_arg;

typedef struct {
    GapIO *io;
    int id;
} cons_world_arg;

typedef struct {
  double min;
  double max;
  int num_ticks;
} ticks_arg;

typedef struct {
    GapIO *io;
    int id;
    char *window;
} del_cons_ruler_arg;

typedef struct {
    GapIO *io;
    int id;
    char *frame;
    char *win;
    int strand;
    int problems;
} strand_arg;

typedef struct {
    GapIO *io;
    char *inlist;
    int band;
} shuffle_arg;

#ifdef USE_BIOLIMS
typedef struct {
  GapIO *io;
  char  *assembly;
} biolims_io_arg;
#endif

typedef struct {
    GapIO *io;
    char *inlist;
    float score;
    int by_consensus;
} abreak_arg;

#endif /* _NEWGAP_STRUCTS_H */
