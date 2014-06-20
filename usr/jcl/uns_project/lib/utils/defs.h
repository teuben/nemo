#define BLOCKSIZE 2048

/*
struct particle {
float pos[3];
float vel[3];
int  id;
};

*/
    
struct particle {
float mass;
float pos[3];
float vel[3];
int  id;
};


struct posit {
  short int pos[3];
};

struct veloc {
 short int vel[3];
};
