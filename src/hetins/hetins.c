#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <unistd.h>
#include "misc.h"

#include "Read.h"
#include "traceType.h"

#define FULL_TEST 2
#define TEST 1

typedef struct HETINS_PARAMS_ {
  int    window;
  int    mode;
  double worst_envelope;
  double worst_half_signal;
  double good_signal_grad_indel;
}HETINS_PARAMS;



void calc_peak_trough_values(Read *r, double *envelope, double *trough_buf) {
  int i, j, k, start_pos, end_pos, m;
  double t, b;

  m = r->NBases - 2;
  for (i = 1; i < m; i++) {
    k = ((r->basePos)[i] - (r->basePos[i-1])) /2;
    start_pos = (r->basePos)[i] - k;
    end_pos   = (r->basePos)[i] + k;
    for (t=-999.0,j=start_pos; j<end_pos; j++) {
      t = MAX(t, envelope[j]);
    }
    k = ((r->basePos)[i] - (r->basePos[i-1])) /4;
    start_pos = (r->basePos)[i-1] + k;
    end_pos   = (r->basePos)[i] + 3*k;
    for (b=9999.0,j=start_pos; j<end_pos; j++) {
      b = MIN(b, envelope[j]);
    }
    if(t>0.0) {
      trough_buf[i] = b / t;
    }
    /*printf("i %d %f %f %f\n",i,b,t,trough_buf[i]);*/
  }
}


int get_envelope(TRACE *tx, TRACE *ty, TRACE *tz, TRACE *tw, double *envelope,
int len_trace)
{
    double x,y,z,w;
    int i;

    for(i=0;i<len_trace;i++) {
      x = tx[i];
      y = ty[i];
      z = tz[i];
      w = tw[i];
      envelope[i] = MAX(MAX(x,y),MAX(z,w));
    }
    return 0;
}

double get_area(TRACE *trace, int startp, int endp)
{
    int i;
    double sum=1.0e-10;

    for (i=startp; i<endp; i++)
	sum += trace[i];

    return(sum);
}

double total_area(TRACE *tx, TRACE *ty, TRACE *tz, TRACE *tw, int stp, int endp)
{
    double x,y,z,w,u;
    x = get_area(tx,stp,endp);
    y = get_area(ty,stp,endp);
    z = get_area(tz,stp,endp);
    w = get_area(tw,stp,endp);
    /*u = (w + x + y + z)/4.0;*/
    u = (w + x + y + z);
    /*printf("%f %f %f %f %f \n",w,x,y,z,u);*/
    return u;
}

/** return sum_area_for_smallest_two_peaks / sum_area_for_largest_two_peaks
 */
double max_area2(TRACE *tx, TRACE *ty, TRACE *tz, TRACE *tw, int stp, int endp)
{
    double x,y,z,w,u,v;
    x = get_area(tx,stp,endp);
    y = get_area(ty,stp,endp);
    z = get_area(tz,stp,endp);
    w = get_area(tw,stp,endp);

    if ((w >= x) && (w >= y) && (w >= z)) {
      if ((x >= y) && (x >= z)) {
	u = w + x;
	v = y + z;
      }
      else {
	if ((y >= x) && (y >= z)) {
	  u = w + y;
	  v = x + z;
	}
	else {
	  if ((z >= x) && (z >= y)) {
	    u = w + z;
	    v = x + y;
	  }
	}
      }
      return (v/u);
    }
    if ((x >= w) && (x >= y) && (x >= z)) {
      if ((w >= y) && (w >= z)) {
	u = x + w;
	v = y + z;
      }
      else {
	if ((y >= w) && (y >= z)) {
	  u = x + y;
	  v = w + z;
	}
	else {
	  if ((z >= w) && (z >= y)) {
	    u = x + z;
	    v = y + w;
	  }
	}
      }
      return (v/u);
    }
    if ((y >= w) && (y >= x) && (y >= z)) {
      if ((x >= w) && (x >= z)) {
	u = y + x;
	v = w + z;
      }
      else {
	if ((w >= x) && (w >= z)) {
	  u = y + w;
	  v = x + z;
	}
	else {
	  if ((z >= x) && (z >= w)) {
	    u = y + z;
	    v = x + w;
	  }
	}
      }
      return (v/u);
    }
    /* must be z */
    if ((x >= y) && (x >= w)) {
      u = z + x;
      v = y + w;
    }
    else {
      if ((y >= x) && (y >= w)) {
	u = z + y;
	v = w + x;
      }
      else {
	if ((w >= x) && (w >= y)) {
	  u = z + w;
	  v = x + y;
	}
      }
    }
    return (v/u);
}

/**
 * return maximum peak area
 */
double max_peak_area(TRACE *tx, TRACE *ty, TRACE *tz, TRACE *tw, int stp, int endp)
{
    double x,y,z,w;
    x = get_area(tx,stp,endp);
    y = get_area(ty,stp,endp);
    z = get_area(tz,stp,endp);
    w = get_area(tw,stp,endp);
    /*printf("%f %f %f %f\n",w,x,y,z);*/

    if ((w >= x) && (w >= y) && (w >= z)) {
      return w;
    }
    if ((x >= w) && (x >= y) && (x >= z)) {
      return x;
    }
    if ((y >= w) && (y >= x) && (y >= z)) {
      return y;
    }
    /* must be z */
    return z;
}

/**
 * smooth using window
 * and make the bits at the edge = first window
 */

void smoothe(double i[], int l, int w) {
  int r, f, m, wo2;
  double s, *a;

  if (NULL == (a = (double *)xcalloc(l, sizeof(double)))) {
    return;
  }

  wo2 = w / 2;

  for(r=0;r<l;r++) a[r] = i[r];
  for(r=0,s=0.0;r<w;r++) s += i[r];

  a[wo2] = s;

  for(r=0,f=w,m=wo2+1;f<l-1;r++,f++,m++) {
    /*printf("r %d f %d m %d Am %f Ir %f If %f\n",r,f,m,a[m-1],i[r],i[f]);*/
    a[m] = a[m-1] - i[r] + i[f];
  }

  for(r=0;r<wo2;r++) a[r] = a[wo2];
  for(r=l-wo2-1;r<l;r++) a[r] = a[l-wo2-2];
  for(r=0;r<l-1;r++) a[r] = a[r] / w;
  a[l-1] = a[l-2];
  for(r=0;r<l;r++) i[r] = a[r];
  xfree(a);
}

/**
 * calc gradient using window
 */

void grad(double i[], double a[], int l, int w) {
  int r, f, m, wo2;

  wo2 = w / 2;

  a[0] = a[l-1] = a[l-2] = 0.0;
  for(r=0,f=w-1,m=wo2;f<l-1;r++,f++,m++) {
    /*printf("r %d f %d m %d Am %f Ir %f If %f\n",r,f,m,a[m-1],i[r],i[f]);*/
    a[m] = (i[f] - i[r])/MAX(1.0,i[f]);
    /*
    printf("%d %f %f\n",m,i[m],a[m]);
    */
  }

}

void get_good_signal(double *signal, double *max_peak, double *good_signal,
		     int good_bases) {
  int i;
  for (i = 0; i < good_bases-1; i++) {
    if (signal[i] > 0.0) {
      good_signal[i] = max_peak[i] / signal[i];
    }
    else {
      good_signal[i] = 0.5;
    }
  }
}

int find_heterozygous_indel(double *half_signal, double *envelope_signal,
			    double *good_signal_grad, double *good_signal_grad_grad,
			    int good_bases, HETINS_PARAMS params) {

  double min_grad= 99.0;
  int min_grad_pos = 0, k;
  int i, end_good = good_bases;
  double worst_envelope = params.worst_envelope;
  double worst_half_signal = params.worst_half_signal;
  double good_signal_grad_indel = params.good_signal_grad_indel;

  /* find where the ratio of two_lowest_peaks/two_highest_peaks gets too high
   * this means 2 worst peaks ~= 2 best ie even for indel data it is crap
   * or where the envelope has small peak_to_trough height
   */
  for (i = 0; i < good_bases-1; i++) {
    if((half_signal[i] > worst_half_signal)||(envelope_signal[i] > worst_envelope)) {
      end_good = i;
      break;
    }
  }

  for (i = 0; i < end_good-1; i++) {
    min_grad = MIN(min_grad,good_signal_grad[i]);
    if(good_signal_grad[i] < good_signal_grad_indel) {
      /* if the gradient of ratio(max_peak/signal) is sufficiently -ve
       * look for it bottoming out within next 5 samples
       * ie in indel region we expect the area of the highest peak
       * divided by the total signal to decline
       */
      for(k=MAX(1,i-1);k<MIN(i+4,good_bases-2);k++) {
	if((good_signal_grad_grad[k-1] < 0.0) && (good_signal_grad_grad[k] > 0.0)) {
	  min_grad_pos = i;
	  goto got_indel;
	}
      }
    }
  }

 got_indel:
  /*
  if(min_grad < good_signal_grad_indel) {
    printf("indel at %d\n",min_grad_pos);
  }
  else {
    printf("indel at 9999\n");
  }
  return 0;
  */
  return min_grad_pos;
}



int heterozygous_indels(Read *r, HETINS_PARAMS params) {
    int i, ret = -1, mode;
    int start_pos,end_pos,good_bases;
    int win_len;
    double *signal, *good_signal, *good_signal_grad, *good_signal_grad_grad;
    double *envelope, *envelope_signal, *half_signal, *max_peak;

    win_len = params.window;
    mode    = params.mode;
    good_bases = r->NBases - 1;

    signal = good_signal = good_signal_grad = good_signal_grad_grad = NULL;
    envelope = envelope_signal = half_signal = max_peak = NULL;
    if (NULL == (signal = (double *)xcalloc(good_bases, sizeof(double))))goto bail_out;
    if (NULL == (max_peak = (double *)xcalloc(good_bases, sizeof(double))))goto bail_out;
    if (NULL == (good_signal = (double *)xcalloc(good_bases, sizeof(double))))goto bail_out;
    if (NULL == (good_signal_grad = (double *)xcalloc(good_bases, sizeof(double))))goto bail_out;
    if (NULL == (good_signal_grad_grad = (double *)xcalloc(good_bases, sizeof(double))))goto bail_out;
    if (NULL == (envelope = (double *)xcalloc(r->NPoints, sizeof(double))))goto bail_out;
    if (NULL == (envelope_signal = (double *)xcalloc(good_bases, sizeof(double))))goto bail_out;
    if (NULL == (half_signal = (double *)xcalloc(good_bases, sizeof(double))))goto bail_out;
    for (i = 1; i < good_bases-1; i++) {
      start_pos = (r->basePos)[i]-((r->basePos)[i] - (r->basePos[i-1])) /2;
      end_pos   = (r->basePos)[i]+((r->basePos)[i+1] - (r->basePos[i])) /2;
      signal[i] = total_area(r->traceC, r->traceG, r->traceA,
			     r->traceT, start_pos, end_pos);
      max_peak[i] = max_peak_area(r->traceC, r->traceG, r->traceA,
				  r->traceT, start_pos, end_pos);
      half_signal[i] = max_area2(r->traceC, r->traceG, r->traceA,
				 r->traceT, start_pos, end_pos);
    }

    if(get_envelope(r->traceC, r->traceG, r->traceA,
		    r->traceT, envelope, r->NPoints)) goto bail_out;
    calc_peak_trough_values(r, envelope, envelope_signal);
    smoothe(signal,good_bases-1,win_len);
    smoothe(max_peak,good_bases-1,win_len);
    smoothe(envelope_signal,good_bases-1,win_len);
    smoothe(half_signal,good_bases-1,win_len);
    get_good_signal(signal, max_peak, good_signal, good_bases);
    grad(good_signal,good_signal_grad,good_bases-1,win_len);
    grad(good_signal_grad,good_signal_grad_grad,good_bases-1,10);
    if(mode == FULL_TEST) {
      printf("\n");
      for (i = 0; i < good_bases-1; i++) {
	printf("%d %f %f %f %f %f %f %f\n",
	       i,signal[i],max_peak[i],good_signal[i],
	       good_signal_grad[i],envelope_signal[i],
	       good_signal_grad_grad[i],half_signal[i]);
      }
      ret = 0;
    }
    else {
      ret = find_heterozygous_indel(half_signal, envelope_signal,
				    good_signal_grad, good_signal_grad_grad,
				    good_bases, params);
    }
 bail_out:
    xfree(good_signal_grad_grad);
    xfree(good_signal_grad);
    xfree(good_signal);
    xfree(half_signal);
    xfree(signal);
    xfree(max_peak);
    xfree(envelope_signal);
    xfree(envelope);
    return ret;
}

void usage(HETINS_PARAMS params) {

    fprintf(stderr,
	    "Usage: hetins [options] file_name\n"
	    "Where options are:\n"
	    "    [-w window_length (%d)]           [-e worst_envelope (%f)]\n"
	    "    [-h worst_half_signal (%f)]  [-g good_signal_grad_indel (%f)]\n"
	    "    [-t debug only]       file_name\n",
	    params.window, params.worst_envelope, params.worst_half_signal,
	    params.good_signal_grad_indel);
    exit(1);
}

/* 8/1/99 johnt - must explictly import globals from DLLs with Visual C++ */
#ifdef _MSC_VER
#  define DLL_IMPORT __declspec(dllimport)
#else
#  define DLL_IMPORT
#endif

#define BUFSIZE 1024

int main(int argc, char **argv) {
  char buffer[BUFSIZE];
  char *fn;
  int c;
  HETINS_PARAMS params;
  Read *read = NULL;
  Exp_info *exp_file = NULL;
  int file_type, ret;
  int qr;
  extern DLL_IMPORT char *optarg;
  extern DLL_IMPORT int optind;

  params.window = 101;
  params.worst_envelope = 0.5;
  params.worst_half_signal     = 0.15;
  params.good_signal_grad_indel = -0.146;
  params.mode = 0;


  while ((c = getopt(argc, argv, "w:e:h:g:tT")) != -1) {
    switch (c) {
    case 'T':
      params.mode = FULL_TEST;
      break;
    case 't':
      params.mode = TEST;
      break;
    case 'w':
      params.window = atoi(optarg);
      break;
    case 'e':
      params.worst_envelope = atof(optarg);
      break;
    case 'h':
      params.worst_half_signal = atof(optarg);
      break;
    case 'g':
      params.good_signal_grad_indel = atof(optarg);
      break;
    default:
      usage(params);
    }
  }
  if (optind == argc) usage(params);

  fn = argv[optind];
  file_type = determine_trace_type(fn);

  if ((file_type == TT_PLN) || (file_type == TT_UNK)) {
    fprintf(stderr,"Input file not EXP or trace\n");
    goto bail_out;
  }

  if (file_type == TT_EXP) {

    if(NULL ==(exp_file = exp_read_info(fn))) {
      fprintf(stderr, "Couldn't read reading file %s\n",fn);
      goto bail_out;
    }
    /* get the QR value for later */

    if (exp_Nentries(exp_file, EFLT_QR))
	qr = atoi(exp_get_entry ( exp_file, EFLT_QR ));
    else
	qr = 0;

    /* Extract LN record, gives us input trace file name */
    if( exp_get_str(exp_file,EFLT_LN,buffer,BUFSIZE) ) {
      fprintf( stderr, "Unable to read LN record from experiment file %s,\n", fn);
      goto bail_out;
    }


    /* Open input trace */

    if(NULL ==(read = read_reading( buffer, TT_ANY ))) {
      fprintf( stderr, "Unable to read trace for experiment file %s,\n", fn);
      goto bail_out;
    }
  }
  else {
    if(NULL ==(read = read_reading( fn, TT_ANY ))) {
      fprintf( stderr, "Unable to read trace file %s,\n", fn);
      goto bail_out;
    }
    params.mode = TEST;
  }

  ret = heterozygous_indels(read, params);

  if (params.mode == TEST ) {
    printf("%s %d\n",fn,ret);
  }
  else if (params.mode != FULL_TEST ) {

    fprintf(stdout,"%d\n",ret);
    if (ret && (file_type == TT_EXP) && (params.mode != 1)) {
      /* write out a tag */
      sprintf(buffer, "HETI = %d..%d\n %d %5.3f %5.3f %6.3f",
	      ret,read->NBases-1,params.window,params.worst_envelope,
	      params.worst_half_signal,params.good_signal_grad_indel);
      if (exp_put_str(exp_file, EFLT_TG, buffer, strlen(buffer))) {
	fprintf( stderr, "Unable to write to experiment file %s,\n", fn);
	goto bail_out;
      }
      /* Only shift QR to the right, never to the left */
      /* if qr left of the indel we won't see it in gap4! */
      if( ret > qr ) {
	sprintf(buffer, "%d", ret+1);
	if (exp_put_str(exp_file, EFLT_QR, buffer, strlen(buffer))) {
	  fprintf( stderr, "Unable to write to experiment file %s,\n", fn);
	  goto bail_out;
	}
      }
    }
  }
  read_deallocate(read);
  exp_destroy_info ( exp_file );
  return 0;

  bail_out:
  read_deallocate(read);
  exp_destroy_info ( exp_file );
  return -1;
}

