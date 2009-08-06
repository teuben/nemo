/* 
	ConvertSnap.cpp   Andrew Davis
		This file takes a Gadget Snapshot of type 1 and puts it out as a snapshot
		of type 2 -- including the 4 character header to each block of data.
		
		WARNING: I've assumed you only have particles of type 0 (gas) and 1 (halo).  
		WARNING: I've also assumend you have not turned on the output potential, accellerations, entropy
				 or timesteps in the makefile
 
		


*/


#include<iostream>
#include<cstdio>
#include<cstdlib>

#define SKIP fread(&dummy, sizeof(dummy), 1, fp);

#define PUT fwrite(&dummy, sizeof(dummy), 1, fout);

struct io_header {
	long npart[6];
	double mass[6];
	double time;
	double redshift;
	int flag_sfr;
	int flag_feedback;
	long npartTotal[6];
	int flag_cooling;
	int num_files;
	double BoxSize;
	double Omega0;
	double OmegaL;
	double H0; // in km/s/Mpc
	char fill[256-6*4-6*8-2*8-2*4-6*4-2*4-4*8]; //fill to 256 bytes
};

void show(io_header &in);

using namespace std;

int main(int argc, char * argv[]) {
	
	FILE *fp, *fout;
	string InFile,OutFile;
	unsigned int junk,ct,id,size;
	int dummy,i,j, blksize, next_blksize;
	long ntot, ngas, ndm, nmass;
	io_header header1;
	float tmp;
	bool individual_mass = false;
	
	if (argc != 3) {
		cout << "Wrong number of inputs!" << endl;
		cout << "Usage: ./ConvertSnap INFILE OUTFILE\n";
		exit(1);
	}
	InFile=argv[1];
	OutFile=argv[2];
	
	
   	fp=fopen(InFile.c_str(),"r");
	fout=fopen(OutFile.c_str(),"wb");
   	if (fp == NULL) {cout<<"Bad readin of file!\nGood Bye\n";
      	exit(2);
   	}
	
	
	blksize=sizeof(char)*4+sizeof(int);
	next_blksize=sizeof(header1)+2*sizeof(int);
	fwrite(&blksize,sizeof(blksize),1,fout);
	fwrite("HEAD", 4*sizeof(char), 1, fout);
	fwrite(&next_blksize,sizeof(next_blksize),1,fout);
	fwrite(&blksize,sizeof(blksize),1,fout);
   	
	SKIP
	fread(&header1.npart,sizeof(header1.npart),1,fp);
   	fread(&header1.mass,sizeof(header1.mass),1,fp);
   	fread(&header1.time,sizeof(header1.time),1,fp);
   	fread(&header1.redshift,sizeof(header1.redshift),1,fp);
   	fread(&header1.flag_sfr,sizeof(header1.flag_sfr),1,fp);
   	fread(&header1.flag_feedback,sizeof(header1.flag_feedback),1,fp);
   	fread(&header1.npartTotal,sizeof(header1.npartTotal),1,fp);
   	fread(&header1.flag_cooling,sizeof(header1.flag_cooling),1,fp);
   	fread(&header1.num_files,sizeof(header1.num_files),1,fp);
   	fread(&header1.BoxSize,sizeof(header1.BoxSize),1,fp);
   	fread(&header1.Omega0,sizeof(header1.Omega0),1,fp);
   	fread(&header1.OmegaL,sizeof(header1.OmegaL),1,fp);
   	fread(&header1.H0,sizeof(header1.H0),1,fp);
   	fread(&header1.fill,sizeof(header1.fill),1,fp);
	
	cout << "Header size = " << sizeof(header1) << " junk: " << junk << endl;
	
   	SKIP
	
	PUT		
	fwrite(&header1,sizeof(header1),1,fout);
	PUT
	
	
	
	show(header1);
   	cout << endl;
   	
	ntot = header1.npartTotal[0]+header1.npartTotal[1];
	ngas = header1.npartTotal[0];
	ndm = header1.npartTotal[1];
	
	
   	//reading positions:
	// This reads in the positions into the tmp (a float) variable three times, 
	///  and lets you decide where to put it.
	
	blksize=sizeof(char)*4+sizeof(int);
	next_blksize=sizeof(float)*3*ntot+2*sizeof(int);
	fwrite(&blksize,sizeof(blksize),1,fout);
	fwrite("POS ", 4*sizeof(char), 1, fout);
	fwrite(&next_blksize,sizeof(next_blksize),1,fout);
	fwrite(&blksize,sizeof(blksize),1,fout);
   	
	
	
	SKIP
	PUT
		for (j=0;j<ntot;j++) {
			for (int k=0;k<3;k++) {
				fread(&tmp,sizeof(tmp),1,fp);
				fwrite(&tmp,sizeof(tmp),1,fout);
			}
		}
	SKIP	
	PUT
			
   	//reading velocities:
	// See the same note about how this reads it above in positions.	
		
		
	blksize=sizeof(char)*4+sizeof(int);
	next_blksize=sizeof(float)*3*ntot+2*sizeof(int);
	fwrite(&blksize,sizeof(blksize),1,fout);
	fwrite("VEL ", 4*sizeof(char), 1, fout);
	fwrite(&next_blksize,sizeof(next_blksize),1,fout);
	fwrite(&blksize,sizeof(blksize),1,fout);	
	
	SKIP	
	PUT
		for (j=0;j<ntot;j++) {
			for (int k=0;k<3;k++) {
				fread(&tmp,sizeof(tmp),1,fp);
				fwrite(&tmp,sizeof(tmp),1,fout);
			}
		}
	
   	SKIP	
	PUT
			
   	//read id's:
		
	blksize=sizeof(char)*4+sizeof(int);
	next_blksize=sizeof(int)*ntot+2*sizeof(int);
	fwrite(&blksize,sizeof(blksize),1,fout);
	fwrite("ID ", 4*sizeof(char), 1, fout);
	fwrite(&next_blksize,sizeof(next_blksize),1,fout);
	fwrite(&blksize,sizeof(blksize),1,fout);

	SKIP
	PUT
		for (j=0;j<ntot;j++) {
			fread(&id,sizeof(id),1,fp);
			fwrite(&id,sizeof(id),1,fout);
		}
	
   	SKIP
	PUT;
	
	nmass = 0;
	for (i=0;i<2;i++) if (header1.npart[i] != 0 && header1.mass[i] == 0) {
			individual_mass = true;
			nmass+=header1.npart[i];
	}
	
		
		if (individual_mass) {
			blksize=sizeof(char)*4+sizeof(int);
			next_blksize=sizeof(float)*nmass+2*sizeof(int);
			fwrite(&blksize,sizeof(blksize),1,fout);
			fwrite("MASS", 4*sizeof(char), 1, fout);
			fwrite(&next_blksize,sizeof(next_blksize),1,fout);
			fwrite(&blksize,sizeof(blksize),1,fout);
			
			SKIP
			PUT;
			for (j=0;j<nmass;j++) {
				fread(&tmp,sizeof(tmp),1,fp);
				fwrite(&tmp,sizeof(tmp),1,fout);
			}
			SKIP			
			PUT;
		}
		
		if (ngas > 0 ) {
			blksize=sizeof(char)*4+sizeof(int);
			next_blksize=sizeof(float)*ngas+2*sizeof(int);
			fwrite(&blksize,sizeof(blksize),1,fout);
			fwrite("U   ", 4*sizeof(char), 1, fout);
			fwrite(&next_blksize,sizeof(next_blksize),1,fout);
			fwrite(&blksize,sizeof(blksize),1,fout);
			
			SKIP			
			PUT;
			for (j=0;j<ngas;j++) {
				fread(&tmp,sizeof(tmp),1,fp);
				fwrite(&tmp,sizeof(tmp),1,fout);
			}
			SKIP			
			PUT;

			blksize=sizeof(char)*4+sizeof(int);
			next_blksize=sizeof(float)*ngas+2*sizeof(int);
			fwrite(&blksize,sizeof(blksize),1,fout);
			fwrite("RHO ", 4*sizeof(char), 1, fout);
			fwrite(&next_blksize,sizeof(next_blksize),1,fout);
			fwrite(&blksize,sizeof(blksize),1,fout);
			
			SKIP			
			PUT;
			for (j=0;j<ngas;j++) {
				fread(&tmp,sizeof(tmp),1,fp);
				fwrite(&tmp,sizeof(tmp),1,fout);
			}
			SKIP			
			PUT;
			
			blksize=sizeof(char)*4+sizeof(int);
			next_blksize=sizeof(float)*ngas+2*sizeof(int);
			fwrite(&blksize,sizeof(blksize),1,fout);
			fwrite("HSML", 4*sizeof(char), 1, fout);
			fwrite(&next_blksize,sizeof(next_blksize),1,fout);
			fwrite(&blksize,sizeof(blksize),1,fout);
			
			SKIP			
			PUT;
			for (j=0;j<ngas;j++) {
				fread(&tmp,sizeof(tmp),1,fp);
				fwrite(&tmp,sizeof(tmp),1,fout);
			}
			SKIP			
			PUT;
		}
		
		return 0;
	
}

void show(io_header &in) {
	int i=0;
	cout <<"\nShowing current Header: \n";
	for (i=0;i<6;i++)
		cout << in.npart[i] << " ";
	cout << endl;
	for (i=0;i<6;i++) 
		cout << in.mass[i] << " ";
	cout << endl;
	cout << in.time << endl;
	cout << in.redshift << endl;
	cout << in.flag_sfr << endl;
	cout << in.flag_feedback << endl;
	for (i=0;i<6;i++) 
		cout << in.npartTotal[i] << " ";
	cout << endl;
	cout << in.flag_cooling << endl;
	cout << in.num_files << endl;
	cout << in.BoxSize << endl;
	cout << in.Omega0 << endl;
	cout << in.OmegaL << endl;
	cout << in.H0 << endl;
	
	return;
}
