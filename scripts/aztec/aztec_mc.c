/* Grand-canonical monomer-dimer Monte Carlo on the Aztec diamond of order N.
 * Weight z^{#monomers}. Cells (x,y) in [0,2N)^2 with |x-N+1/2|+|y-N+1/2| <= N.
 * Moves (detailed balance, symmetric proposals):
 *   link move: random cell i, random direction d; j = neighbour (reject if outside)
 *     i,j paired          -> break into two monomers     acc min(1,z^2)
 *     i,j both monomers   -> join                        acc min(1,z^-2)
 *     i monomer, j paired with k -> slide: (i,j) paired, k monomer   acc 1   (and i<->j)
 *   plaquette flip: random cell as lower-left of a 2x2 plaquette inside the region.
 * Measurements every `every` sweeps:
 *   Q_top = sum over monomers of colour sign * triangle weight (top triangle, as az_analyze.py)
 *   n_tot, per-cell histogram of state {mono, +x, +y, -x, -y} -> frozenness phi.
 * Output: block averages of Q_top, n_tot, and final per-cell histogram (binary-free text).
 * usage: aztec_mc N z therm_sweeps meas_sweeps nblk seed every outprefix
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

static uint64_t s[4];
static inline uint64_t rotl(const uint64_t x, int k){return (x<<k)|(x>>(64-k));}
static inline uint64_t nxt(void){ const uint64_t r = rotl(s[1]*5,7)*9; const uint64_t t = s[1]<<17;
  s[2]^=s[0]; s[3]^=s[1]; s[1]^=s[2]; s[0]^=s[3]; s[2]^=t; s[3]=rotl(s[3],45); return r;}
static inline double urand(void){return (nxt()>>11)*0x1.0p-53;}
static void seed_rng(uint64_t sd){ for(int i=0;i<4;i++){ sd += 0x9e3779b97f4a7c15ULL; uint64_t z=sd;
  z=(z^(z>>30))*0xbf58476d1ce4e5b9ULL; z=(z^(z>>27))*0x94d049bb133111ebULL; s[i]=z^(z>>31);} }

int N, M, NC; int *inreg, *p, *nb, *cells; double *wtop; int *sgn;
int nmono; double Qtop;

int main(int argc, char **argv){
  if(argc<9){fprintf(stderr,"usage\n");return 1;}
  N=atoi(argv[1]); double z=atof(argv[2]); long therm=atol(argv[3]); long meas=atol(argv[4]);
  int nblk=atoi(argv[5]); seed_rng(atoll(argv[6])); int every=atoi(argv[7]); const char *pre=argv[8];
  M=2*N; inreg=calloc(M*M,sizeof(int)); p=malloc(M*M*sizeof(int)); nb=malloc(4*M*M*sizeof(int));
  wtop=calloc(M*M,sizeof(double)); sgn=calloc(M*M,sizeof(int)); cells=malloc(M*M*sizeof(int)); NC=0;
  for(int y=0;y<M;y++)for(int x=0;x<M;x++){ double u=fabs(x-N+0.5)+fabs(y-N+0.5); if(u<=N+1e-9){inreg[x+M*y]=1; cells[NC++]=x+M*y;} }
  int dx[4]={1,0,-1,0}, dy[4]={0,1,0,-1};
  for(int y=0;y<M;y++)for(int x=0;x<M;x++){int i=x+M*y; for(int d=0;d<4;d++){int X=x+dx[d],Y=y+dy[d];
      nb[4*i+d]=(X>=0&&X<M&&Y>=0&&Y<M&&inreg[X+M*Y])?X+M*Y:-1;}
    if(inreg[i]){ int k=x+y-(N-1), b=x-y, uu=k-N; double w=0;
      if(abs(uu)>abs(b)) w=(uu<0)?1:0; else if(abs(b)>abs(uu)) w=0; else w=(uu<0)?0.5:0;
      wtop[i]=w; sgn[i]=(k%2==0)?1:-1; } }
  /* initial: horizontal dimers in every row (row widths are even) */
  for(int y=0;y<M;y++){ int x0=-1; for(int x=0;x<M;x++) if(inreg[x+M*y]){x0=x;break;}
    if(x0<0) continue; for(int x=x0; x<M && inreg[x+M*y]; x+=2){ int i=x+M*y; p[i]=i+1; p[i+1]=i; } }
  nmono=0; Qtop=0;
  double z2=z*z, iz2=(z>0)?1.0/z2:1e300;
  long *hist=calloc(5*(size_t)M*M,sizeof(long)); long nmeas=0;
  long nsteps=(long)NC;
  for(long sw=0; sw<therm+meas; sw++){
    for(long t=0;t<nsteps;t++){
      if(urand()<0.75){
        int i=cells[nxt()%NC]; int d=nxt()&3; int j=nb[4*i+d]; if(j<0) continue;
        int pi=p[i], pj=p[j];
        if(pi==j){ if(z2>=1.0||urand()<z2){ p[i]=-1;p[j]=-1;nmono+=2; Qtop+=sgn[i]*wtop[i]+sgn[j]*wtop[j]; } }
        else if(pi==-1&&pj==-1){ if(iz2>=1.0||urand()<iz2){ p[i]=j;p[j]=i;nmono-=2; Qtop-=sgn[i]*wtop[i]+sgn[j]*wtop[j]; } }
        else if(pi==-1&&pj>=0){ int k=pj; p[k]=-1; p[i]=j; p[j]=i; Qtop+=sgn[k]*wtop[k]-sgn[i]*wtop[i]; }
        else if(pj==-1&&pi>=0){ int k=pi; p[k]=-1; p[j]=i; p[i]=j; Qtop+=sgn[k]*wtop[k]-sgn[j]*wtop[j]; }
      } else {
        int i=cells[nxt()%NC]; int a=nb[4*i+0], b=nb[4*i+1]; if(a<0||b<0) continue; int c=nb[4*a+1]; if(c<0) continue;
        if(p[i]==a&&p[b]==c){ p[i]=b;p[b]=i;p[a]=c;p[c]=a; }
        else if(p[i]==b&&p[a]==c){ p[i]=a;p[a]=i;p[b]=c;p[c]=b; }
      }
    }
    if(sw>=therm && (sw-therm)%every==0){
      for(int q=0;q<NC;q++){ int i=cells[q]; int st=0; if(p[i]>=0){ for(int d=0;d<4;d++) if(nb[4*i+d]==p[i]) {st=1+d;break;} } hist[5*i+st]++; }
      nmeas++;
    }
    static double aQ=0, an=0; static long cnt=0;
    if(sw>=therm){ aQ+=Qtop; an+=nmono; cnt++;
      long per=meas/nblk; if(cnt==per){ printf("B Q %.8f n %.8f\n", aQ/cnt, an/cnt/NC); fflush(stdout); aQ=0; an=0; cnt=0; } }
  }
  char fn[512]; snprintf(fn,sizeof fn,"%s_hist.txt",pre); FILE *f=fopen(fn,"w");
  fprintf(f,"# N %d z %g nmeas %ld\n",N,z,nmeas);
  for(int q=0;q<NC;q++){int i=cells[q]; fprintf(f,"%d %d %ld %ld %ld %ld %ld\n",i%M,i/M,hist[5*i],hist[5*i+1],hist[5*i+2],hist[5*i+3],hist[5*i+4]);}
  fclose(f); return 0;
}
