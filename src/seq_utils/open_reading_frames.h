#ifndef _OPEN_READING_FRAMES_H_
#define _OPEN_READING_FRAMES_H_

int write_open_frames_f ( FILE *fp, char *dna, int num_bases, 
			  int user_start, int user_end, int min_open );
int write_screen_open_frames_f ( char *dna, int num_bases, int user_start,
				 int user_end, int min_open );
void write_open_frames_f_ft ( FILE *fp, char *dna, int num_bases, 
			      int user_start, int user_end, int min_open );
void write_screen_open_frames_f_ft (char *dna, int num_bases, int user_start,
				    int user_end, int min_open );
int write_open_frames_r ( FILE *fp, char *dna, int num_bases, 
			  int user_start,int user_end, int min_open );
int write_screen_open_frames_r ( char *dna, int num_bases, int user_start,
				 int user_end, int min_open );
void write_open_frames_r_ft ( FILE *fp, char *dna, int num_bases, 
			      int user_start, int user_end, int min_open );
void write_screen_open_frames_r_ft (char *dna, int num_bases, int user_start,
				    int user_end, int min_open );
#endif
