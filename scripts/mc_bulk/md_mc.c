/* Grand-canonical monomer-dimer model on an L x L periodic square lattice.
 * Weight: z^{#monomers}. Moves:
 *   link move on random (i, dir):
 *     dimer (i,j)            -> remove, two monomers   acc min(1, z^2)
 *     monomers at i and j     -> add dimer              acc min(1, z^-2)
 *     monomer at i, j in (j,k)-> slide: (i,j) dimer, monomer to k   acc 1
 *     (and mirror with i<->j)
 *   plaquette flip of two parallel dimers.
 * Measurements: monomer density n, wall-wall charge correlator Cq(x),
 *   wall-wall density correlator Cd(x) (raw), charge structure factor Sq(k_m).
 * Output per block (for jackknife).
 * usage: md_mc L z therm_sweeps meas_sweeps nblocks seed
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

static uint64_t s[4];
static inline uint64_t rotl(const uint64_t x, int k){return (x<<k)|(x>>(64-k));}
static inline uint64_t nxt(void){
  const uint64_t r = rotl(s[1]*5,7)*9; const uint64_t t = s[1]<<17;
  s[2]^=s[0]; s[3]^=s[1]; s[1]^=s[2]; s[0]^=s[3]; s[2]^=t; s[3]=rotl(s[3],45); return r;}
static inline double urand(void){return (nxt()>>11)*0x1.0p-53;}
static void seed_rng(uint64_t sd){ for(int i=0;i<4;i++){ sd += 0x9e3779b97f4a7c15ULL; uint64_t z=sd;
  z=(z^(z>>30))*0xbf58476d1ce4e5b9ULL; z=(z^(z>>27))*0x94d049bb133111ebULL; s[i]=z^(z>>31);} }

int L, N; int *p; int *nb; int nmono;

static inline int nbr(int i,int d){ return nb[4*i+d]; }

void link_move(double z2, double iz2){
  int i = nxt()%N; int d = nxt()%4; int j = nbr(i,d);
  int pi=p[i], pj=p[j];
  if(pi==j){ if(z2>=1.0 || urand()<z2){ p[i]=-1; p[j]=-1; nmono+=2; } return; }
  if(pi==-1 && pj==-1){ if(iz2>=1.0 || urand()<iz2){ p[i]=j; p[j]=i; nmono-=2; } return; }
  if(pi==-1 && pj>=0){ int k=pj; p[k]=-1; p[i]=j; p[j]=i; return; }
  if(pj==-1 && pi>=0){ int k=pi; p[k]=-1; p[j]=i; p[i]=j; return; }
}
void plaq_move(void){
  int i = nxt()%N; int a=nbr(i,0), b=nbr(i,1), c=nbr(a,1); /* i a / b c */
  if(p[i]==a && p[b]==c){ p[i]=b; p[b]=i; p[a]=c; p[c]=a; }
  else if(p[i]==b && p[a]==c){ p[i]=a; p[a]=i; p[b]=c; p[c]=b; }
}

int main(int argc,char**argv){
  if(argc<7){fprintf(stderr,"usage\n");return 1;}
  L=atoi(argv[1]); double z=atof(argv[2]); long therm=atol(argv[3]); long meas=atol(argv[4]);
  int nblk=atoi(argv[5]); seed_rng(atoll(argv[6]));
  N=L*L; p=malloc(N*sizeof(int)); nb=malloc(4*N*sizeof(int));
  for(int y=0;y<L;y++)for(int x=0;x<L;x++){int i=x+L*y;
    nb[4*i+0]=((x+1)%L)+L*y; nb[4*i+1]=x+L*((y+1)%L);
    nb[4*i+2]=((x-1+L)%L)+L*y; nb[4*i+3]=x+L*((y-1+L)%L);}
  for(int y=0;y<L;y++)for(int x=0;x<L;x+=2){int i=x+L*y; p[i]=i+1; p[i+1]=i;}
  nmono=0;
  double z2=z*z, iz2=1.0/z2;
  int H=L/2; int KM=12;
  double *rq=malloc(L*sizeof(double)), *cq=malloc(L*sizeof(double));
  double *rd=malloc(L*sizeof(double)), *cd=malloc(L*sizeof(double));
  double *aCq=calloc(H+1,sizeof(double)), *aCd=calloc(H+1,sizeof(double));
  double *aSq=calloc(KM+1,sizeof(double)), *aSd=calloc(KM+1,sizeof(double));
  double an=0, aD=0; long cnt=0;
  for(long sw=0; sw<therm; sw++) for(int t=0;t<N;t++){ if(urand()<0.8) link_move(z2,iz2); else plaq_move(); }
  long per = meas/nblk;
  printf("# L %d z %g therm %ld meas %ld nblk %d\n",L,z,therm,meas,nblk);
  for(int b=0;b<nblk;b++){
    for(int k=0;k<=H;k++){aCq[k]=0;aCd[k]=0;} for(int k=0;k<=KM;k++){aSq[k]=0;aSd[k]=0;} an=0;aD=0;cnt=0;
    for(long sw=0; sw<per; sw++){
      for(int t=0;t<N;t++){ if(urand()<0.8) link_move(z2,iz2); else plaq_move(); }
      /* measure: rows (sum over y at fixed x) and columns */
      for(int dirn=0; dirn<2; dirn++){
        for(int a=0;a<L;a++){rq[a]=0;rd[a]=0;}
        for(int y=0;y<L;y++)for(int x=0;x<L;x++){int i=x+L*y; if(p[i]<0){
            int a = dirn? y : x; int e=((x+y)&1)?-1:1; rq[a]+=e; rd[a]+=1; }}
        for(int r=0;r<=H;r++){double sq=0,sd=0; for(int a=0;a<L;a++){int a2=(a+r)%L; sq+=rq[a]*rq[a2]; sd+=rd[a]*rd[a2];}
          aCq[r]+=sq/L; aCd[r]+=sd/L;}
        for(int m=1;m<=KM;m++){double cr=0,ci=0,dr=0,di=0; double k=2*M_PI*m/L;
          for(int a=0;a<L;a++){cr+=rq[a]*cos(k*a); ci+=rq[a]*sin(k*a); dr+=rd[a]*cos(k*a); di+=rd[a]*sin(k*a);}
          aSq[m]+=(cr*cr+ci*ci)/N; aSd[m]+=(dr*dr+di*di)/N;}
        double D=0; for(int a=0;a<L;a++) D+=rd[a]; aD+=D/L;
        cnt++;
      }
      an += (double)nmono/N;
    }
    printf("B %d n %.10g D %.10g\n", b, an/per, aD/cnt);
    printf("Cq"); for(int r=0;r<=H;r++) printf(" %.10g",aCq[r]/cnt); printf("\n");
    printf("Cd"); for(int r=0;r<=H;r++) printf(" %.10g",aCd[r]/cnt); printf("\n");
    printf("Sq"); for(int m=1;m<=KM;m++) printf(" %.10g",aSq[m]/cnt); printf("\n");
    printf("Sd"); for(int m=1;m<=KM;m++) printf(" %.10g",aSd[m]/cnt); printf("\n");
    fflush(stdout);
  }
  return 0;
}
