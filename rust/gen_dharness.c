#include <stdio.h>
#include <stdlib.h>
#include <string.h>
extern double D__align( double **, char **, char **, double *, double *, int, int, int, int, double *, char *, char *, char *, char *, int *, int, int *, int, int );
extern void constants( int nseq, char **seq );
extern int njob, nlenmax, penalty, penalty_ex, penalty_shift, nalphabets;
extern double **n_dis_consweight_multi;
extern char *newgapstr;
extern int outgap, legacygapcost, trywarp, specificityconsideration;
extern char dorp;
extern double penalty_shift_factor;
extern int ppenalty, ppenalty_dist, ppenalty_OP, ppenalty_ex, ppenalty_EX, poffset, pamN, kimuraR;
extern int nblosum, scoremtx, fmodel, nwildcard;
#define NOTSPEC 100009
static const char *B(double v){ static char b[32]; unsigned long long u; memcpy(&u,&v,8); snprintf(b,32,"%016llx",u); return b; }
int main( int argc, char **argv )
{
    int n1 = atoi(argv[1]); int a = 2;
    char **s1 = calloc(n1+1, sizeof(char*)); double *e1 = calloc(n1, sizeof(double));
    for (int i=0;i<n1;i++) s1[i]=argv[a++];
    for (int i=0;i<n1;i++) e1[i]=atof(argv[a++]);
    int n2 = atoi(argv[a++]);
    char **s2 = calloc(n2+1, sizeof(char*)); double *e2 = calloc(n2, sizeof(double));
    for (int i=0;i<n2;i++) s2[i]=argv[a++];
    for (int i=0;i<n2;i++) e2[i]=atof(argv[a++]);
    int headgp = atoi(argv[a++]), tailgp = atoi(argv[a++]);
    int l1=strlen(s1[0]), l2=strlen(s2[0]);
    int alloclen = (l1+l2)*2 + 100;
    njob = n1+n2; nlenmax = l1>l2?l1:l2;
    // build the scoring matrix the same way disttbfast does
    char **all = calloc(njob+1, sizeof(char*));
    for(int i=0;i<n1;i++) all[i]=s1[i];
    for(int i=0;i<n2;i++) all[n1+i]=s2[i];
    dorp = 'p';
    ppenalty = -1530; ppenalty_dist = NOTSPEC; ppenalty_OP = NOTSPEC;
    ppenalty_ex = 0; ppenalty_EX = NOTSPEC; poffset = 0;
    pamN = NOTSPEC; kimuraR = NOTSPEC;
    nblosum = 62; scoremtx = 1; fmodel = 0; nwildcard = 0;
    penalty_shift_factor = 100.0;
    constants( njob, all );
    outgap = 1; legacygapcost = 0; trywarp = 0;
    // grow the rows to alloclen, as the caller in treebase does
    char **g1 = calloc(n1+1,sizeof(char*)), **g2 = calloc(n2+1,sizeof(char*));
    for(int i=0;i<n1;i++){ g1[i]=calloc(alloclen+10,1); strcpy(g1[i], s1[i]); }
    for(int i=0;i<n2;i++){ g2[i]=calloc(alloclen+10,1); strcpy(g2[i], s2[i]); }
    double dumdb = 0.0;
    double sc = D__align( n_dis_consweight_multi, g1, g2, e1, e2, n1, n2, alloclen,
                          0, &dumdb, NULL, NULL, NULL, NULL, NULL, 0, NULL, headgp, tailgp );
    printf("score %s\n", B(sc));
    for(int i=0;i<n1;i++) printf("1 %d %s\n", i, g1[i]);
    for(int i=0;i<n2;i++) printf("2 %d %s\n", i, g2[i]);
    return 0;
}
