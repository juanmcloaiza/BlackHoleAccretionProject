#include <stdio.h>

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "readgadget.h"
#define int4bytes int
#define GAMMA 7.0/5.0
/*--------- comment/uncomment to remove/enable DEBUG outputs ------------------*/

#define MY_DEBUG


/*-----------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------*/
/*---------------------------- Low Level Routines -----------------------------*/
/*-----------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------*/

int4bytes blksize,swap=0;
#define SKIPr  {my_fread(&blksize,sizeof(int),1,fd); swap_Nbyte((char*)&blksize,1,4);}
#define SKIPw  {my_fwrite(&blksize,sizeof(int),1,fp);}

/*---------------------- Basic routine to read data from a file ---------------*/
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
	size_t nread;

	if((nread = fread(ptr, size, nmemb, stream)) != nmemb)
	{
		printf("I/O error (fread) !\n");
		exit(3);
	}
	return nread;
}

/*---------------------- Basic routine to write data to a file ---------------*/
size_t my_fwrite(const void *ptr, size_t size, size_t nmemb, FILE * stream)
{
	size_t nwrite;

	if((nwrite = fwrite(ptr, size, nmemb, stream)) != nmemb)
	{
		printf("I/O error (fwrite) !\n");
		exit(3);
	}
	return nwrite;
}

/*-----------------------------------------------------------------------------*/
/*---------------------- Routine to swap ENDIAN -------------------------------*/
/*-------- char *data:    Pointer to the data ---------------------------------*/
/*-------- int n:         Number of elements to swap --------------------------*/
/*-------- int m:         Size of single element to swap ----------------------*/
/*--------                int,float = 4 ---------------------------------------*/
/*--------                double    = 8 ---------------------------------------*/
/*-----------------------------------------------------------------------------*/
void swap_Nbyte(char *data,int n,int m)
{
	int i,j;
	char old_data[16];

	if(swap>0)
	{
		for(j=0;j<n;j++)
		{
			memcpy(&old_data[0],&data[j*m],m);
			for(i=0;i<m;i++)
			{
				data[j*m+i]=old_data[m-i-1];
			}
		}
	}
}

/*-----------------------------------------------------------------------------*/
/*---------------------- Routine find a block in a snapfile -------------------*/
/*-------- FILE *fd:      File handle -----------------------------------------*/
/*-------- char *label:   4 byte identifyer for block -------------------------*/
/*-------- returns length of block found, -------------------------------------*/
/*-------- the file fd points to starting point of block ----------------------*/
/*-----------------------------------------------------------------------------*/
int find_block(FILE *fd,char *label)
{
	int4bytes blocksize=0;
	char blocklabel[5]={"    "};

	rewind(fd);

	while(!feof(fd) && blocksize == 0)
	{
		SKIPr;
		if(blksize == 134217728)
		{
#ifdef MY_DEBUG
			printf("Enable ENDIAN swapping !\n");
#endif
			swap=1-swap;
			swap_Nbyte((char*)&blksize,1,4);
		}
		if(blksize != 8)
		{
			printf("incorrect format (blksize=%d)!\n",blksize);
			exit(1);
		}
		else
		{
			my_fread(blocklabel, 4*sizeof(char), 1, fd);
			my_fread(&blocksize, sizeof(int4bytes), 1, fd);
			swap_Nbyte((char*)&blocksize,1,4);
#ifdef MY_DEBUG
			printf("Found Block <%s> with %d bytes\n",blocklabel,blocksize);
#endif
			SKIPr;
			if(strcmp(label,blocklabel)!=0)
			{ 
				fseek(fd,blocksize,1);
				blocksize=0;
			}
		}
	}
	return(blocksize-8);
}

/*-----------------------------------------------------------------------------*/
/*----------------- Routine cover for a block in a snapfile -------------------*/
/*-------- FILE *fd:      File handle -----------------------------------------*/
/*-------- char *label:   4 byte identifyer for block -------------------------*/
/*-------- returns length of block found, -------------------------------------*/
/*-------- the file fd points to starting point of block ----------------------*/
/*-----------------------------------------------------------------------------*/
void cover_block(FILE *fp,char *label, int blocksize)
{
	int totblock;
	blksize = 8;
	totblock = blocksize + blksize;
	printf("Writing cover for <%s> of %d bytes\n", label, blocksize);
	SKIPw;
	my_fwrite(label, 4*sizeof(char), 1, fp);
	my_fwrite(&totblock, sizeof(int4bytes), 1, fp);
	SKIPw;
}

/*-----------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------*/
/*---------------------------- High Level Routines ----------------------------*/
/*-----------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------*/
/*---------------------- Routine to read the header information ---------------*/
/*-------- int *npart:    List of Particle numbers for spezies [0..5] ---------*/
/*-------- int *massarr:  List of masses for spezies [0..5] -------------------*/
/*-------- int *time:     Time of snapshot ------------------------------------*/
/*-------- int *massarr:  Redshift of snapshot --------------------------------*/
/*-------- FILE *fd:      File handle -----------------------------------------*/
/*-------- returns number of read bytes ---------------------------------------*/
/*-----------------------------------------------------------------------------*/
int read_gadget_head(int *npart,double *massarr,double *time,double *redshift,FILE *fd)
{
	int blocksize,dummysize;
	int i;

	blocksize = find_block(fd,"HEAD");
	if(blocksize <= 0)
	{
		printf("Block <%s> not fond !\n","HEAD");
		exit(5);
	}
	else
	{
		dummysize=blocksize - 6 * sizeof(int) - 8 * sizeof(double);
		SKIPr;
		my_fread(npart,6*sizeof(int), 1, fd);        swap_Nbyte((char*)npart,6,4);
		for(i=0;i<6;i++) printf("npart[%d] = %d\n",i,npart[i]);
		my_fread(massarr,6*sizeof(double), 1, fd);   swap_Nbyte((char*)massarr,6,8);
		for(i=0;i<6;i++) printf("massarr[%d] = %f\n",i,massarr[i]);
		my_fread(time,sizeof(double), 1, fd);        swap_Nbyte((char*)time,1,8);
		my_fread(redshift,sizeof(double), 1, fd);    swap_Nbyte((char*)redshift,1,8);
		fseek(fd,dummysize,1);
		SKIPr;
	}
	return(blocksize);
}

/*-----------------------------------------------------------------------------*/
/*--------------------- Routine to write the header information ---------------*/
/*-------- int *npart:    List of Particle numbers for spezies [0..5] ---------*/
/*-------- int *massarr:  List of masses for spezies [0..5] -------------------*/
/*-------- int *time:     Time of snapshot ------------------------------------*/
/*-------- int *massarr:  Redshift of snapshot --------------------------------*/
/*-------- FILE *fd:      File handle -----------------------------------------*/
/*-------- returns number of read bytes ---------------------------------------*/
/*-----------------------------------------------------------------------------*/
int write_gadget_head(int *npart,double *massarr,double *time,double *redshift,FILE *fp)
{
	int blocksize,dummysize;
	char blocklabel[5]={"HEAD"};
	printf("writing header\n");
	blocksize = 256;
	if(blocksize <= 0)
	{
		printf("Block <%s> not written !\n","HEAD");
		exit(5);
	}
	else
	{
		dummysize=blocksize - 6 * sizeof(int) - 8 * sizeof(double);

		cover_block(fp,blocklabel,blocksize);
		blksize = blocksize; 
		SKIPw;
		my_fwrite(npart,6*sizeof(int), 1, fp);  
		my_fwrite(massarr,6*sizeof(double), 1, fp); 
		my_fwrite(time,sizeof(double), 1, fp);     
		my_fwrite(redshift,sizeof(double), 1, fp);
		fseek(fp,dummysize,1);
		SKIPw;
	}
	printf("done writing header\n");
	return(blocksize);
}

/*-----------------------------------------------------------------------------*/
/*---------------------- Routine to read a 1D float array ---------------------*/
/*-------- int *data:     Pointer where the data are stored to ----------------*/
/*-------- char *label:   Identifyer for the datafield to read ----------------*/
/*-------- FILE *fd:      File handle -----------------------------------------*/
/*-------- returns length of dataarray ----------------------------------------*/
/*-----------------------------------------------------------------------------*/
int read_gadget_float(float *data,char *label,FILE *fd)
{
	int blocksize;

	blocksize = find_block(fd,label);
	if(blocksize <= 0)
	{
		printf("Block <%s> not fond !\n",label);
		exit(5);
	}
	else
	{
#ifdef MY_DEBUG
		printf("Reading %d bytes of data from <%s>...\n",blocksize,label);
#endif
		SKIPr;
		my_fread(data,blocksize, 1, fd);
		swap_Nbyte((char*)data,blocksize/sizeof(float),4);
		SKIPr;
	}
	return(blocksize/sizeof(float));
}

/*-----------------------------------------------------------------------------*/
/*--------------------- Routine to write a 1D float array ---------------------*/
/*-------- int *data:     Pointer where the data are stored to ----------------*/
/*-------- char *label:   Identifyer for the datafield to read ----------------*/
/*-------- FILE *fd:      File handle -----------------------------------------*/
/*-------- returns length of dataarray ----------------------------------------*/
/*-----------------------------------------------------------------------------*/
int write_gadget_float(int npart, float *data,char *label,FILE *fp)
{
	int blocksize;
	char blocklabel[5];
	strcpy(blocklabel,label);


	blocksize = npart * sizeof(float);
	if(blocksize <= 0)
	{
		printf("Block <%s> couldn't be written!\n",label);
		exit(5);
	}
	else
	{
#ifdef MY_DEBUG
		printf("Writing %d bytes of data to <%s>...\n",blocksize,label);
#endif
		cover_block(fp,label,blocksize);
		blksize = blocksize; 
		SKIPw;
		my_fwrite(data,blocksize, 1, fp);
		SKIPw;
	}
	return(blocksize/sizeof(float));
}

/*-----------------------------------------------------------------------------*/
/*---------------------- Routine to read a 1D int array -----------------------*/
/*-------- int *data:     Pointer where the data are stored to ----------------*/
/*-------- char *label:   Identifyer for the datafield to read ----------------*/
/*-------- FILE *fd:      File handle -----------------------------------------*/
/*-------- returns length of dataarray ----------------------------------------*/
/*-----------------------------------------------------------------------------*/
int read_gadget_int(int *data,char *label,FILE *fd)
{
	int blocksize;

	blocksize = find_block(fd,label);
	if(blocksize <= 0)
	{
		printf("Block <%s> not fond !\n",label);
		exit(5);
	}
	else
	{
#ifdef MY_DEBUG
		printf("Reading %d bytes of data from <%s>...\n",blocksize,label);
#endif
		SKIPr;
		my_fread(data,blocksize, 1, fd);
		swap_Nbyte((char*)data,blocksize/sizeof(int),4);
		SKIPr;
	}
	return(blocksize/sizeof(int));
}

/*-----------------------------------------------------------------------------*/
/*---------------------- Routine to write a 1D int array -----------------------*/
/*-------- int *data:     Pointer where the data are stored to ----------------*/
/*-------- char *label:   Identifyer for the datafield to write ----------------*/
/*-------- FILE *fd:      File handle -----------------------------------------*/
/*-------- returns length of dataarray ----------------------------------------*/
/*-----------------------------------------------------------------------------*/
int write_gadget_int(int npart, int *data,char *label,FILE *fp)
{
	int blocksize;
	char blocklabel[5];
	strcpy(blocklabel,label);

	blocksize = npart*sizeof(int);
	if(blocksize <= 0)
	{
		printf("Block <%s> couldn't be written!\n",label);
		exit(5);
	}
	else
	{
#ifdef MY_DEBUG
		printf("Writing %d bytes of data to <%s>...\n",blocksize,label);
#endif 
		cover_block(fp,label,blocksize);
		blksize = blocksize; 
		SKIPw;
		my_fwrite(data,blocksize, 1, fp);
		SKIPw;
	}
	return(blocksize/sizeof(int));
}

/*-----------------------------------------------------------------------------*/
/*---------------------- Routine to read a 3D float array ---------------------*/
/*-------- int *data:     Pointer where the data are stored to ----------------*/
/*-------- char *label:   Identifyer for the datafield to read ----------------*/
/*-------- FILE *fd:      File handle -----------------------------------------*/
/*-------- returns length of dataarray ----------------------------------------*/
/*-----------------------------------------------------------------------------*/
int read_gadget_float3(float *data,char *label,FILE *fd)
{
	int blocksize;

	blocksize = find_block(fd,label);
	if(blocksize <= 0)
	{
		printf("Block <%s> not fond !\n",label);
		exit(5);
	}
	else
	{
#ifdef MY_DEBUG
		printf("Reding %d bytes of data from <%s>...\n",blocksize,label);
#endif
		SKIPr;
		my_fread(data,blocksize, 1, fd);
		swap_Nbyte((char*)data,blocksize/sizeof(float),4);
		SKIPr;
	}
	return(blocksize/sizeof(float)/3);
}

/*-----------------------------------------------------------------------------*/
/*---------------------- Routine to write a 3D float array ---------------------*/
/*-------- int *data:     Pointer where the data are stored to ----------------*/
/*-------- char *label:   Identifyer for the datafield to write ----------------*/
/*-------- FILE *fd:      File handle -----------------------------------------*/
/*-------- returns length of dataarray ----------------------------------------*/
/*-----------------------------------------------------------------------------*/
int write_gadget_float3(int npart, float *data,char *label,FILE *fp)
{
	int blocksize;
	char blocklabel[5];
	strcpy(blocklabel,label);

	blocksize = 3*npart * sizeof(float);
	if(blocksize <= 0)
	{
		printf("Block <%s> couldn't be written!\n",label);
		exit(5);
	}
	else
	{
#ifdef MY_DEBUG
		printf("Writing %d bytes of data to <%s>...\n",blocksize,label);
#endif
		cover_block(fp,label,blocksize);
		blksize = blocksize; 
		SKIPw;
		my_fwrite(data,blocksize, 1, fp);
		SKIPw;
	}
	return(blocksize/sizeof(float)/3);
}

/*-----------------------------------------------------------------------------*/
/*----------- Routine to make a histogram from some array ---------------------*/
/*-------- float *x:      Pointer where the data are stored to ----------------*/
/*-------- int Nsamples:  Number of data samples ------------------------------*/
/*-------- float *xbin: bin array (calculated here) ---------------------------*/
/*-------- int *F:  frequency array (calculated here) -------------------------*/
/*---------int Nbins: desired number of bins ----------------------------------*/
/*-----------------------------------------------------------------------------*/
void make_histogram(float *x, int Nsamples, float *xbin, int *F, int Nbins){
	int i,j, count;
	float xmin = 1e10;
	float xmax = -1e10;
	float dL;
	/*Check xmin and xmax*/
	for( i = 0 ; i < Nsamples ; i++ ){
		if (x[i] < xmin) xmin = x[i];
		if (x[i] > xmax) xmax = x[i];
	}
	dL = (xmax - xmin)/(float)Nbins;
	for( j = 0; j < Nbins ; j++ ){
		xbin[j] = xmin + j * dL;
	}
	for( j = 0 ; j < Nbins ; j++ ){
		count = 0;
		for( i = 0 ; i < Nsamples ; i++ ){
			if ( ( x[i] > xbin[j] ) && ( x[i] <= xbin[j] + dL) ){
				count = count + 1;
			}
			if(j == 0 && x[i] == xmin) count = count + 1;
		}
		F[j] = count;
	} 
}

/*-----------------------------------------------------------------------------*/
/*---- Routine to make a cumulative distribution from some array --------------*/
/*-------- float *x:      Pointer where the data are stored to ----------------*/
/*-------- int Nsamples:  Number of data samples ------------------------------*/
/*-------- float *xbin: bin array (calculated here) ---------------------------*/
/*-------- int *F:  cumulative frequency array (calculated here) --------------*/
/*-------- int Nbins: desired number of bins ----------------------------------*/
/*-----------------------------------------------------------------------------*/
void make_cumulative_histogram(float *x, int Nsamples, float *xbin, int *F, int Nbins){
	int i,j, count;
	float xmin = 1e10;
	float xmax = -1e10;
	float dL;
	/*Check xmin and xmax*/
	for( i = 0 ; i < Nsamples ; i++ ){
		if (x[i] < xmin) xmin = x[i];
		if (x[i] > xmax) xmax = x[i];
	}
	dL = (xmax - xmin)/(float)Nbins;
	for( j = 0; j < Nbins ; j++ ){
		xbin[j] = xmin + j * dL;
	}
	count = 0;
	for( j = 0 ; j < Nbins ; j++ ){
		for( i = 0 ; i < Nsamples ; i++ ){
			if ( ( x[i] > xbin[j] ) && ( x[i] <= xbin[j] + dL) ){
				count = count + 1;
			}
			if(j == 0 && x[i] == xmin) count = count + 1;
		}
		F[j] = count;
	} 
}




/*-----------------------------------------------------------------------------*/
/*-------- Routine to make a distribution of p as function of x ---------------*/
/*-------- float *x:         Pointer to data of the independent variable -------*/
/*-------- float *p:         Pointer to data of the dependent variable ---------*/
/*-------- float *xbin:     Pointer to which the x_i bins are stored ----------*/
/*-------- float *pbin:     Pointer to which the p_i value is stored ----------*/
/*-------- int Nsamples:      Number of data samples ----------------------------*/
/*-------- float Nbins:     Number of desired bins ----------------------------*/
/*-----------------------------------------------------------------------------*/
void get_distribution(float *x, float *p, float *xbin, float *pbin, int Nsamples, int Nbins){
	int i, j, count;
	float xmin = 1e10;
	float xmax = -1e10;
	float dL, p_sum;
	/*Check xmin and xmax*/
	for( i = 0 ; i < Nsamples ; i++ ){
		if (x[i] < xmin) xmin = x[i];
		if (x[i] > xmax) xmax = x[i];
	}
	printf("Warning! overriding xmax and xmin\n");
	xmin = 0.0; xmax = 1.0;
	printf("xmin = %f, xmax = %f\n", xmin, xmax);
	dL = (xmax - xmin)/(float)Nbins;
	/*Binning of x...*/
	for( j = 0; j < Nbins ; j++ ){
		xbin[j] = xmin + j * dL;
	}

	for( j = 0 ; j < Nbins ; j++ ){
		p_sum = 0;
		count = 0;
		for( i = 0 ; i < Nsamples ; i++ ){
			//printf("%f\n",p[i]);
			if ( ( x[i] > xbin[j] ) && ( x[i] <= xbin[j] + dL) ){
				p_sum = p_sum + p[i];
				count  = count + 1;
				//printf("%f\n",p_sum);
			}
		}
		pbin[j] = p_sum/(float)count;
	} 
}

float sign(float x){
	if ( x < 0 ) return -1;
	if ( x > 0 ) return 1;
	return 0;
}


/*-----------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------*/
/*---------------------------- Small Example HowToUse -------------------------*/
/*-----------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------*/

int main(int argc, char **argv)
{
	FILE *fd = 0;
	FILE *fp = 0;
	char filein[265];
	char fileout[265];
	int i,n,ntot;
	int Nbins = 1000;
	int npart[6];
	double masstab[6],redshift,time;
	float (*pos)[3], (*vel)[3];
	int *id;
	float *r_, *l_, (*lxyz)[3];
	int *Nr, *Nl;
	float *rbin, *lbin;
	
	//Abort if arguments are not what they should:
	if(argc != 3){
		printf("Usage: %s <snapshot_input> <filename_output>\n",argv[0]);
		exit(4);
	}
	strcpy( filein,  argv[1]);
	strcpy( fileout, argv[2]);
	//Abort if not able to open files:
	if(!(fd = fopen(filein,"r"))){
		printf("Cannot open file <%s> !\n",filein);
		exit(2);
	}
	if(!(fp = fopen(fileout,"w+"))){
		printf("Cannot open file <%s> !\n",fileout);
		exit(2);
	}
	/*----------- RED HEADER TO GET GLOBAL PROPERTIES -------------*/
	n = read_gadget_head(npart,masstab,&time,&redshift,fd);

	ntot=0;
	for(i=0;i<6;i++){
		printf("PartSpezies %d, anz=%d, masstab=%f\n",i,npart[i],masstab[i]);
		ntot += npart[i];
	}
	printf("Time of snapshot=%f, z=%f, ntot=%d\n",time,redshift,ntot);

	/*---------- ALLOCATE MEMORY ---------------------------------*/
	pos = malloc(3*ntot*sizeof(float));
	vel = malloc(3*ntot*sizeof(float));
	id = malloc(ntot*sizeof(int));
	
	lxyz = malloc(3*ntot*sizeof(float));
	r_ = malloc(ntot*sizeof(float));
	l_ = malloc(ntot*sizeof(float));

	rbin = malloc(Nbins * sizeof(float));
	lbin =  malloc(Nbins * sizeof(float));
	Nr =  malloc(Nbins * sizeof(int));
	Nl =  malloc(Nbins * sizeof(int));

	/*---------- READ DATA BLOCKS --------------------------------*/
	n = read_gadget_float3((float*)pos,"POS ",fd);
	n = read_gadget_float3((float*)vel,"VEL ",fd);
	n = read_gadget_int(id,"ID  ",fd);

	/*---------- DO STUFF ----------------------------------------*/
	
	/*Compute radius and angular momentum*/
	for(i=0;i<npart[0];i++){
		lxyz[i][0] = (pos[i][1] * vel[i][2] - pos[i][2] * vel[i][1]);
		lxyz[i][1] = (pos[i][2] * vel[i][0] - pos[i][0] * vel[i][2]);
    	lxyz[i][2] = (pos[i][0] * vel[i][1] - pos[i][1] * vel[i][0]);
		r_[i] = sqrt( pow(pos[i][0],2) + pow(pos[i][1],2) + pow(pos[i][2],2) ); 
		l_[i] = sqrt( pow(lxyz[i][0],2) + pow(pos[i][1],2) + pow(pos[i][2],2) ) * sign(lxyz[i][2]); 
	}

	/*Make histograms*/
	printf("making histogram for r...\n");
	make_cumulative_histogram(r_, ntot, rbin, Nr, Nbins);
	printf("making histogram for l...\n");
	make_cumulative_histogram(l_, ntot, lbin, Nl, Nbins);
	printf("writing output file...\n");
	fprintf(fp,"#rbin[i] Nr[i] lbin[i] Nl[i]\n");
	for(i = 0 ; i < Nbins ; i++){
		fprintf(fp," %f %f %f %f \n",rbin[i],(float)Nr[i]/(float)ntot,lbin[i],(float)Nl[i]/(float)ntot);
	}

	/*Free bins*/

	/*---------- CLOSE FILES AND FREE DATA ------------------------*/
	printf("Closing files...\n");
	fclose(fd);
	fclose(fp);

	printf("Freing memory...\n");
	free(pos);
	free(vel);
	free(id);
	
	free(lxyz);
	free(r_);
	free(l_);

	free(rbin);
	free(Nr);
	free(lbin); 
	free(Nl);


	printf("Output file writen to %s\n", fileout);
	return 0;
} 
