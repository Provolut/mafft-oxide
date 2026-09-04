#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static const char *B(double v){ static char b[4][32]; static int t=0; t=(t+1)&3; unsigned long long u; memcpy(&u,&v,8); snprintf(b[t],32,"%016llx",u); return b[t]; }
typedef struct _gaplenvec { int relend; int idatnext; int idatprev; int npat; int len; double freq; } Gaplen;
extern void gaplencount( int n, int l, Gaplen **mtx, char **seq, double *eff );
extern double calcpfacnoidatend( Gaplen **g1, Gaplen **g2, int i, int j, char *s1, char *s2, int disp );
extern double calcpfac_gap_noidatend( Gaplen **g1, Gaplen **g2, int newgaplen, int i, int j, char *s1, char *s2, int disp );
extern double calcpfac_gapex_noidatend( Gaplen **g1, Gaplen **g2, int i, int j, int newgaplen, char *s1, char *s2, int disp );
int main( int argc, char **argv )
{
    // argv: n1 seqs... eff... n2 seqs... eff...
    int n1 = atoi(argv[1]);
    char **s1 = calloc(n1, sizeof(char*));
    double *e1 = calloc(n1, sizeof(double));
    int a = 2;
    for (int i=0;i<n1;i++) s1[i]=argv[a++];
    for (int i=0;i<n1;i++) e1[i]=atof(argv[a++]);
    int n2 = atoi(argv[a++]);
    char **s2 = calloc(n2, sizeof(char*));
    double *e2 = calloc(n2, sizeof(double));
    for (int i=0;i<n2;i++) s2[i]=argv[a++];
    for (int i=0;i<n2;i++) e2[i]=atof(argv[a++]);
    int l1 = strlen(s1[0]), l2 = strlen(s2[0]);
    Gaplen **g1 = calloc(l1+2, sizeof(Gaplen*));
    Gaplen **g2 = calloc(l2+2, sizeof(Gaplen*));
    gaplencount(n1, l1, g1, s1, e1);
    gaplencount(n2, l2, g2, s2, e2);
    for (int j=0;j<=l1;j++) {
        if (!g1[j]) { printf("1 %d NULL\n", j); continue; }
        for (int k=0; g1[j][k].idatnext!=-1; k++)
            printf("1 %d %d len=%d relend=%d freq=%s next=%d prev=%d npat=%d\n",
                j,k,g1[j][k].len,g1[j][k].relend,B(g1[j][k].freq),g1[j][k].idatnext,g1[j][k].idatprev,g1[j][k].npat);
    }
    for (int j=0;j<=l2;j++) {
        if (!g2[j]) { printf("2 %d NULL\n", j); continue; }
        for (int k=0; g2[j][k].idatnext!=-1; k++)
            printf("2 %d %d len=%d relend=%d freq=%s next=%d prev=%d npat=%d\n",
                j,k,g2[j][k].len,g2[j][k].relend,B(g2[j][k].freq),g2[j][k].idatnext,g2[j][k].idatprev,g2[j][k].npat);
    }
    for (int i=0;i<=l1;i++) for (int j=0;j<=l2;j++)
        printf("pf %d %d %s gap=%s gapex=%s\n", i, j,
            B(calcpfacnoidatend(g1,g2,i,j,s1[0],s2[0],0)),
            B(calcpfac_gap_noidatend(g1,g2,2,i,j,s1[0],s2[0],0)),
            B(calcpfac_gapex_noidatend(g1,g2,i,j,2,s1[0],s2[0],0)));
    return 0;
}
