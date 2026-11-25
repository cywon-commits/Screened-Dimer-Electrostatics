// Metropolis sampling of mu_n(x) ~ prod_{i<j} |sin(pi (x_i-x_j))| * prod_i g_i on a ring of length 1 (g_i = nearest-neighbour gaps).
// Output: kappa_n = <(1/n) sum_i 1/g_i> / n, with block errors.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
static unsigned long long s[2];
static inline double ur(void){ unsigned long long a=s[0],b=s[1]; s[0]=b; a^=a<<23; s[1]=a^b^(a>>17)^(b>>26); return ((s[1]+b)>>11)*(1.0/9007199254740992.0); }
int main(int argc,char**argv){
  int n=atoi(argv[1]); long sw=atol(argv[2]); s[0]=88172645463325252ULL+atoi(argv[3]); s[1]=0x9E3779B97F4A7C15ULL;
  double *x=malloc(n*sizeof(double)); for(int i=0;i<n;i++) x[i]=(i+0.5)/n;
  int nb=20; double bs[20]={0}; long per=sw/nb;
  for(long t=-sw/10;t<sw;t++){
    for(int m=0;m<n;m++){
      int i=(int)(ur()*n); int il=(i+n-1)%n, ir=(i+1)%n, ill=(i+n-2)%n, irr=(i+2)%n;
      double lo=x[il]-(i==0), hi=x[ir]+(i==n-1);          // unwrapped neighbours
      double y=lo+(hi-lo)*ur(); double xo=x[i];
      double dl=0; for(int j=0;j<n;j++) if(j!=i) dl+=log(fabs(sin(M_PI*(y-x[j]))/sin(M_PI*(xo-x[j]))));
      dl+=log((y-lo)*(hi-y)/((xo-lo)*(hi-xo)));
      (void)ill;(void)irr;
      if(dl>=0||ur()<exp(dl)) x[i]=y;
    }
    if(t>=0){ double a=0; for(int i=0;i<n;i++){ double g=x[(i+1)%n]-x[i]+(i==n-1); a+=1.0/g; } bs[t/per<nb?t/per:nb-1]+=a/n/n/per; }
  }
  double m=0,v=0; for(int b=0;b<nb;b++) m+=bs[b]/nb; for(int b=0;b<nb;b++) v+=(bs[b]-m)*(bs[b]-m);
  printf("n=%d kappa=%.5f +- %.5f\n",n,m,sqrt(v/nb/(nb-1))); return 0; }
