extern struct vector
{
  float x,y,z;
} *pos, *vel;
extern int swap;

extern void swap_Nbyte(char *data,int n,int m);
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream);
extern int find_block(FILE *fd,char *label);
extern int read_gadget_float(float *data,char *label,FILE *fd);
extern int read_gadget_int(int *data,char *label,FILE *fd);
extern int read_gadget_float3(float *data,char *label,FILE *fd);
extern int read_gadget_head(int *npart,double *massarr,double *time,double *redshift,FILE *fd);
/*
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream);
extern int write_gadget_float(float *data,char *label,FILE *fd);
extern int write_gadget_int(int *data,char *label,FILE *fd);
extern int write_gadget_float3(float *data,char *label,FILE *fd);
extern int write_gadget_head(int *npart,double *massarr,double *time,double *redshift,FILE *fd);
*/
void make_histogram(float *x, int Nsamples, float *xbin, int *F, int Nbins);
