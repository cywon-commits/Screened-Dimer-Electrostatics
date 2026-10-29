/* Parallel-tempering (replica exchange in fugacity z) grand-canonical monomer-dimer MC
 * on the Aztec diamond of order N.  Replica r has fugacity z_r (ladder read from the command line).
 * Within a replica: local moves of aztec_mc.c (link: break/join/slide, plaquette flip).
 * Every sweep: attempt swaps between neighbouring replicas r, r+1 with
 *     P = min(1, (z_r / z_{r+1})^(n_{r+1} - n_r))      (n = number of monomers)
 * Measurements per replica: block averages of Q_top and monomer density; per-cell state histograms.
 * usage: aztec_pt N therm meas nblk seed every outprefix z_1 z_2 ... z_R   (z ascending)
 * Output: <outprefix>_r<k>.txt (blocks), <outprefix>_r<k>_hist.txt, swap acceptance on stdout.
 * Checkpoint: <outprefix>.ckpt is written every AZTEC_CKPT_SEC seconds (default 300) at a block boundary;
 * rerunning the same command resumes from it (bitwise identical to an uninterrupted run).
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

static uint64_t s[4];
static inline uint64_t rotl(const uint64_t x, int k){return (x<<k)|(x>>(64-k));}
static inline uint64_t nxt(void){ const uint64_t r = rotl(s[1]*5,7)*9; const uint64_t t = s[1]<<17;
  s[2]^=s[0]; s[3]^=s[1]; s[1]^=s[2]; s[0]^=s[3]; s[2]^=t; s[3]=rotl(s[3],45); return r;}
static inline double urand(void){return (nxt()>>11)*0x1.0p-53;}
static void seed_rng(uint64_t sd){ for(int i=0;i<4;i++){ sd += 0x9e3779b97f4a7c15ULL; uint64_t z=sd;
  z=(z^(z>>30))*0xbf58476d1ce4e5b9ULL; z=(z^(z>>27))*0x94d049bb133111ebULL; s[i]=z^(z>>31);} }

typedef struct { int *p; int nmono; double Q; double z, z2, iz2; } Rep;
int N, M, NC; int *inreg, *nb, *cells, *sgn; double *wtop;

static void sweep(Rep *R){
  int *p=R->p;
  for(long t=0;t<NC;t++){
    if(urand()<0.75){
      int i=cells[nxt()%NC]; int d=nxt()&3; int j=nb[4*i+d]; if(j<0) continue;
      int pi=p[i], pj=p[j];
      if(pi==j){ if(R->z2>=1.0||urand()<R->z2){ p[i]=-1;p[j]=-1;R->nmono+=2; R->Q+=sgn[i]*wtop[i]+sgn[j]*wtop[j]; } }
      else if(pi==-1&&pj==-1){ if(R->iz2>=1.0||urand()<R->iz2){ p[i]=j;p[j]=i;R->nmono-=2; R->Q-=sgn[i]*wtop[i]+sgn[j]*wtop[j]; } }
      else if(pi==-1&&pj>=0){ int k=pj; p[k]=-1; p[i]=j; p[j]=i; R->Q+=sgn[k]*wtop[k]-sgn[i]*wtop[i]; }
      else if(pj==-1&&pi>=0){ int k=pi; p[k]=-1; p[j]=i; p[i]=j; R->Q+=sgn[k]*wtop[k]-sgn[j]*wtop[j]; }
    } else {
      int i=cells[nxt()%NC]; int a=nb[4*i+0], b=nb[4*i+1]; if(a<0||b<0) continue; int c=nb[4*a+1]; if(c<0) continue;
      if(p[i]==a&&p[b]==c){ p[i]=b;p[b]=i;p[a]=c;p[c]=a; }
      else if(p[i]==b&&p[a]==c){ p[i]=a;p[a]=i;p[b]=c;p[c]=b; }
    }
  }
}

int main(int argc, char **argv){
  if(argc<9){fprintf(stderr,"usage: aztec_pt N therm meas nblk seed every outprefix z1 z2 ...\n");return 1;}
  N=atoi(argv[1]); long therm=atol(argv[2]); long meas=atol(argv[3]); int nblk=atoi(argv[4]);
  seed_rng(atoll(argv[5])); int every=atoi(argv[6]); const char *pre=argv[7];
  int NR=argc-8; double *zl=malloc(NR*sizeof(double)); for(int r=0;r<NR;r++) zl[r]=atof(argv[8+r]);
  M=2*N; inreg=calloc(M*M,sizeof(int)); nb=malloc(4*M*M*sizeof(int)); cells=malloc(M*M*sizeof(int));
  wtop=calloc(M*M,sizeof(double)); sgn=calloc(M*M,sizeof(int)); NC=0;
  for(int y=0;y<M;y++)for(int x=0;x<M;x++){ double u=fabs(x-N+0.5)+fabs(y-N+0.5); if(u<=N+1e-9){inreg[x+M*y]=1; cells[NC++]=x+M*y;} }
  int dx[4]={1,0,-1,0}, dy[4]={0,1,0,-1};
  for(int y=0;y<M;y++)for(int x=0;x<M;x++){int i=x+M*y; for(int d=0;d<4;d++){int X=x+dx[d],Y=y+dy[d];
      nb[4*i+d]=(X>=0&&X<M&&Y>=0&&Y<M&&inreg[X+M*Y])?X+M*Y:-1;}
    if(inreg[i]){ int k=x+y-(N-1), b=x-y, uu=k-N; double w;
      if(abs(uu)>abs(b)) w=(uu<0)?1:0; else if(abs(b)>abs(uu)) w=0; else w=(uu<0)?0.5:0;
      wtop[i]=w; sgn[i]=(k%2==0)?1:-1; } }
  /* replicas: slot r always holds fugacity zl[r]; configurations are swapped */
  Rep *R=calloc(NR,sizeof(Rep));
  for(int r=0;r<NR;r++){ R[r].p=malloc(M*M*sizeof(int)); R[r].nmono=0; R[r].Q=0;
    R[r].z=zl[r]; R[r].z2=zl[r]*zl[r]; R[r].iz2=(zl[r]>0)?1.0/R[r].z2:1e300;
    for(int y=0;y<M;y++){ int x0=-1; for(int x=0;x<M;x++) if(inreg[x+M*y]){x0=x;break;}
      if(x0<0) continue; for(int x=x0;x<M&&inreg[x+M*y];x+=2){int i=x+M*y; R[r].p[i]=i+1; R[r].p[i+1]=i;} } }
  long **hist=malloc(NR*sizeof(long*)); for(int r=0;r<NR;r++) hist[r]=calloc(5*(size_t)M*M,sizeof(long));
  double *aQ=calloc(NR,sizeof(double)), *an=calloc(NR,sizeof(double)); long cnt=0, per=meas/nblk;
  long *sw_try=calloc(NR,sizeof(long)), *sw_acc=calloc(NR,sizeof(long));
  FILE **fo=malloc(NR*sizeof(FILE*)); char fn[512], ck[512], cktmp[520];
  snprintf(ck,sizeof ck,"%s.ckpt",pre); snprintf(cktmp,sizeof cktmp,"%s.tmp",ck);
  double ckint=getenv("AZTEC_CKPT_SEC")?atof(getenv("AZTEC_CKPT_SEC")):300; time_t tck=time(NULL);
  long sw0=0; FILE *fc=fopen(ck,"rb");
  if(fc){ /* resume: header must match this run exactly */
    int hN,hNR; long hth,hme,hper; int ok=fread(&hN,sizeof hN,1,fc)==1&&fread(&hNR,sizeof hNR,1,fc)==1&&fread(&hth,sizeof hth,1,fc)==1
      &&fread(&hme,sizeof hme,1,fc)==1&&fread(&hper,sizeof hper,1,fc)==1&&hN==N&&hNR==NR&&hth==therm&&hme==meas&&hper==per;
    double *hz=malloc(NR*sizeof(double)); ok=ok&&fread(hz,sizeof(double),NR,fc)==(size_t)NR&&memcmp(hz,zl,NR*sizeof(double))==0; free(hz);
    ok=ok&&fread(&sw0,sizeof sw0,1,fc)==1&&fread(s,sizeof s,1,fc)==1;
    for(int r=0;ok&&r<NR;r++) ok=fread(R[r].p,sizeof(int),M*M,fc)==(size_t)(M*M)&&fread(&R[r].nmono,sizeof(int),1,fc)==1
      &&fread(&R[r].Q,sizeof(double),1,fc)==1&&fread(hist[r],sizeof(long),5*(size_t)M*M,fc)==5*(size_t)M*M;
    ok=ok&&fread(sw_try,sizeof(long),NR,fc)==(size_t)NR&&fread(sw_acc,sizeof(long),NR,fc)==(size_t)NR;
    fclose(fc);
    if(!ok){fprintf(stderr,"checkpoint %s does not match this run; remove it to start over\n",ck);return 1;}
    long nbw=sw0>therm?(sw0-therm)/per:0; /* keep header + blocks written before the checkpoint */
    for(int r=0;r<NR;r++){ snprintf(fn,sizeof fn,"%s_r%d.txt",pre,r); FILE *f=fopen(fn,"r"); char line[256], *keep=calloc(1,(nbw+2)*256); long nl=0;
      while(f&&nl<=nbw&&fgets(line,sizeof line,f)){strcat(keep,line);nl++;} if(f)fclose(f);
      if(nl!=nbw+1){fprintf(stderr,"%s has %ld lines, checkpoint expects %ld\n",fn,nl,nbw+1);return 1;}
      fo[r]=fopen(fn,"w"); fputs(keep,fo[r]); fflush(fo[r]); free(keep); }
    fprintf(stderr,"resumed from %s at sweep %ld\n",ck,sw0);
  } else
    for(int r=0;r<NR;r++){ snprintf(fn,sizeof fn,"%s_r%d.txt",pre,r); fo[r]=fopen(fn,"w"); fprintf(fo[r],"# N %d z %g\n",N,zl[r]); fflush(fo[r]); }
  for(long sw=sw0; sw<therm+meas; sw++){
    for(int r=0;r<NR;r++) sweep(&R[r]);
    for(int r=(int)(sw&1); r+1<NR; r+=2){
      double lnP=(R[r+1].nmono-R[r].nmono)*log(zl[r]/zl[r+1]);
      sw_try[r]++;
      if(lnP>=0||urand()<exp(lnP)){
        int *tp=R[r].p; R[r].p=R[r+1].p; R[r+1].p=tp;
        int tn=R[r].nmono; R[r].nmono=R[r+1].nmono; R[r+1].nmono=tn;
        double tq=R[r].Q; R[r].Q=R[r+1].Q; R[r+1].Q=tq; sw_acc[r]++; }
    }
    if(sw>=therm){
      for(int r=0;r<NR;r++){ aQ[r]+=R[r].Q; an[r]+=R[r].nmono; }
      cnt++;
      if(cnt==per){ for(int r=0;r<NR;r++){ fprintf(fo[r],"B Q %.8f n %.8f\n",aQ[r]/cnt,an[r]/cnt/NC); fflush(fo[r]); aQ[r]=0; an[r]=0; } cnt=0; }
      if((sw-therm)%every==0)
        for(int r=0;r<NR;r++) for(int q=0;q<NC;q++){ int i=cells[q]; int st=0; int *p=R[r].p;
          if(p[i]>=0){ for(int d=0;d<4;d++) if(nb[4*i+d]==p[i]){st=1+d;break;} } hist[r][5*i+st]++; }
    }
    if(cnt==0&&sw+1<therm+meas&&difftime(time(NULL),tck)>=ckint){ /* state after sweep sw, at a block boundary */
      long nx=sw+1; FILE *f=fopen(cktmp,"wb"); fwrite(&N,sizeof N,1,f); fwrite(&NR,sizeof NR,1,f); fwrite(&therm,sizeof therm,1,f);
      fwrite(&meas,sizeof meas,1,f); fwrite(&per,sizeof per,1,f); fwrite(zl,sizeof(double),NR,f); fwrite(&nx,sizeof nx,1,f); fwrite(s,sizeof s,1,f);
      for(int r=0;r<NR;r++){ fwrite(R[r].p,sizeof(int),M*M,f); fwrite(&R[r].nmono,sizeof(int),1,f); fwrite(&R[r].Q,sizeof(double),1,f);
        fwrite(hist[r],sizeof(long),5*(size_t)M*M,f); }
      fwrite(sw_try,sizeof(long),NR,f); fwrite(sw_acc,sizeof(long),NR,f);
      if(fflush(f)==0&&fsync(fileno(f))==0&&fclose(f)==0) rename(cktmp,ck); else fprintf(stderr,"checkpoint write failed\n");
      tck=time(NULL); }
  }
  remove(ck);
  printf("swap acceptance:"); for(int r=0;r+1<NR;r++) printf(" %.3f",sw_try[r]?(double)sw_acc[r]/sw_try[r]:0); printf("\n");
  for(int r=0;r<NR;r++){ fclose(fo[r]); snprintf(fn,sizeof fn,"%s_r%d_hist.txt",pre,r); FILE *f=fopen(fn,"w");
    fprintf(f,"# N %d z %g\n",N,zl[r]);
    for(int q=0;q<NC;q++){int i=cells[q]; long *h=hist[r]+5*i; fprintf(f,"%d %d %ld %ld %ld %ld %ld\n",i%M,i/M,h[0],h[1],h[2],h[3],h[4]);}
    fclose(f); }
  return 0;
}
