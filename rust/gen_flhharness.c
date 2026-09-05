/* Stage one `Falign_localhom()` call in the real C and print its result.
 *
 * Compiled against a built `core/` by verify_flh.sh and by gen_flhgolden.py;
 * `examples/flhref.rs` takes the same argument vector and runs the port. Every
 * number is printed as a raw f64 bit pattern, so a formatting difference
 * cannot hide a value difference.
 *
 * The configuration staged here is the only one `dvtditr` reaches:
 * kobetsubunkatsu=1, alg='A', constraint=2, scoringmatrices=NULL. See
 * `GOAL-C-FREE.md`.
 *
 * argv: n1 <n1 rows> <n1 effs> <n1 orinums>
 *       n2 <n2 rows> <n2 effs> <n2 orinums>
 *       outgap legacygapcost fastathreshold fftWinSize fftThreshold
 *       nlh <nlh x (i j start1 end1 start2 end2 opt importance)>
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mltaln.h"

extern double Falign_localhom( int **, double ***, double **, char **, char **, double *, double *,
                               double **, double **, int, int, int, int, double *, int *, int *,
                               int *, int, int * );
extern void part_imp_match_init_strict( double *, int, int, int, int, char **, char **, double *,
                                        double *, double *, double *, LocalHom ***, char *, int,
                                        int *, int * );
extern void constants( int nseq, char **seq );
extern double part_imp_match_out_sc( int i1, int j1 );

#define NOTSPEC 100009

static const char *bits( double v )
{
    static char b[32];
    unsigned long long u;
    memcpy( &u, &v, 8 );
    snprintf( b, 32, "%016llx", u );
    return b;
}

int main( int argc, char **argv )
{
    int a = 1;
    int n1 = atoi( argv[a++] );
    char **s1 = calloc( n1 + 1, sizeof( char * ) );
    double *e1 = calloc( n1, sizeof( double ) );
    int *o1 = calloc( n1, sizeof( int ) );
    for( int i = 0; i < n1; i++ ) s1[i] = argv[a++];
    for( int i = 0; i < n1; i++ ) e1[i] = atof( argv[a++] );
    for( int i = 0; i < n1; i++ ) o1[i] = atoi( argv[a++] );
    int n2 = atoi( argv[a++] );
    char **s2 = calloc( n2 + 1, sizeof( char * ) );
    double *e2 = calloc( n2, sizeof( double ) );
    int *o2 = calloc( n2, sizeof( int ) );
    for( int i = 0; i < n2; i++ ) s2[i] = argv[a++];
    for( int i = 0; i < n2; i++ ) e2[i] = atof( argv[a++] );
    for( int i = 0; i < n2; i++ ) o2[i] = atoi( argv[a++] );
    int outgapv = atoi( argv[a++] );
    int legacyv = atoi( argv[a++] );
    double fthr = atof( argv[a++] );
    int fftwin = atoi( argv[a++] );
    int fftthr = atoi( argv[a++] );
    int nlh = atoi( argv[a++] );

    int lgth = strlen( s1[0] );
    int alloclen = lgth * 4 + 100;

    njob = n1 + n2;
    nlenmax = lgth;
    dorp = 'p';
    ppenalty = -1530; ppenalty_dist = NOTSPEC; ppenalty_OP = NOTSPEC;
    ppenalty_ex = 0; ppenalty_EX = NOTSPEC; poffset = 0;
    pamN = NOTSPEC; kimuraR = NOTSPEC;
    nblosum = 62; scoremtx = 1; fmodel = 0; nwildcard = 0;
    penalty_shift_factor = 100.0;
    fftscore = 1;
    fftWinSize = fftwin;   /* what a traced `mafft --auto` run reports: 20 */
    fftThreshold = fftthr; /* `dvtditr -z 50` */
    fftNoAnchStop = 0;
    fftRepeatStop = 0;
    kobetsubunkatsu = 1; /* the only configuration dvtditr reaches */
    alg = 'A';
    fastathreshold = fthr;
    newgapstr = "-";

    char **all = calloc( njob + 1, sizeof( char * ) );
    for( int i = 0; i < n1; i++ ) all[i] = s1[i];
    for( int i = 0; i < n2; i++ ) all[n1 + i] = s2[i];
    constants( njob, all );
    outgap = outgapv;
    legacygapcost = legacyv;
    trywarp = 0;

    /* grow the rows to alloclen, as tditeration's localcopy does */
    char **g1 = calloc( n1 + 1, sizeof( char * ) );
    char **g2 = calloc( n2 + 1, sizeof( char * ) );
    for( int i = 0; i < n1; i++ ) { g1[i] = calloc( alloclen + 10, 1 ); strcpy( g1[i], s1[i] ); }
    for( int i = 0; i < n2; i++ ) { g2[i] = calloc( alloclen + 10, 1 ); strcpy( g2[i], s2[i] ); }

    double *ekozo1 = calloc( n1, sizeof( double ) );
    double *ekozo2 = calloc( n2, sizeof( double ) );

    LocalHom ***lh = calloc( n1, sizeof( LocalHom ** ) );
    for( int i = 0; i < n1; i++ ) lh[i] = calloc( n2, sizeof( LocalHom * ) );
    for( int k = 0; k < nlh; k++ )
    {
        int i = atoi( argv[a++] ), j = atoi( argv[a++] );
        LocalHom *node = calloc( 1, sizeof( LocalHom ) );
        node->start1 = atoi( argv[a++] );
        node->end1 = atoi( argv[a++] );
        node->start2 = atoi( argv[a++] );
        node->end2 = atoi( argv[a++] );
        node->opt = atof( argv[a++] );
        node->importance = atof( argv[a++] );
        node->rimportance = node->importance;
        node->overlapaa = node->end1 - node->start1 + 1;
        node->extended = -1;
        node->korh = 'h';
        node->nokori = 0;
        node->next = NULL;
        node->last = node;
        if( lh[i][j] == NULL ) lh[i][j] = node;
        else { LocalHom *p = lh[i][j]; while( p->next ) p = p->next; p->next = node; }
    }

    part_imp_match_init_strict( NULL, n1, n2, lgth, lgth, g1, g2, e1, e2, ekozo1, ekozo2,
                                lh, NULL, 1, o1, o2 );

    if( getenv( "FLH_DUMP_IMP" ) )
    {
        double oimp = 0.0;
        for( int i = lgth - 1; i >= 0; i-- ) oimp += part_imp_match_out_sc( i, i );
        printf( "oimp %s\n", bits( oimp ) );
    }

    int *gapmap1 = calloc( alloclen + 10, sizeof( int ) );
    int *gapmap2 = calloc( alloclen + 10, sizeof( int ) );
    double totalimp = 0.0;
    double sc = Falign_localhom( NULL, NULL, n_dis_consweight_multi, g1, g2, e1, e2, NULL, NULL,
                                 n1, n2, alloclen, 2, &totalimp, gapmap1, gapmap2, NULL, 0, NULL );

    printf( "score %s\n", bits( sc ) );
    printf( "imp %s\n", bits( totalimp ) );
    for( int i = 0; i < n1; i++ ) printf( "1 %d %s\n", i, g1[i] );
    for( int i = 0; i < n2; i++ ) printf( "2 %d %s\n", i, g2[i] );
    return 0;
}
