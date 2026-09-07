/* Stage a guide tree in the real C and print every branch weight, distance and
 * path weight `tditeration` reads from it.
 *
 * Compiled against a built `core/` by verify_td.sh; `examples/tdref.rs` takes
 * the same argument vector and runs the port. Every number is printed as a raw
 * f64 bit pattern, so a formatting difference cannot hide a value difference.
 *
 * argv: nseq  then nseq-1 x ( members0 members1 len0 len1 ), where a member
 *       list is comma-separated and this harness appends the C's -1 terminator.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mltaln.h"

extern void treeCnv( Node *stopol, int locnseq, int ***topol, double **len, double **bw );
extern void calcBranchWeight( double **bw, int locnseq, Node *stopol, int ***topol, double **len );
extern void distFromABranch( int nseq, double *result, Node *stopol, int ***topol, double **len, int step, int LorR );
extern void weightFromABranch( int nseq, double *result, Node *stopol, int ***topol, int step, int LorR );
extern void OneClusterAndTheOther_fast( int locnjob, int *memlist1, int *memlist2, int *s1, int *s2, char *pair, int ***topol, int step, int branch, double **smalldistmtx, double **distmtx, double *distontree );

static const char *bits( double v )
{
    static char b[32];
    unsigned long long u;
    memcpy( &u, &v, 8 );
    snprintf( b, 32, "%016llx", u );
    return b;
}

/* "0,3,7" -> { 0, 3, 7, -1 } */
static int *parse_members( const char *s, int nseq )
{
    int *v = calloc( nseq + 1, sizeof( int ) );
    int n = 0;
    const char *p = s;
    while( *p )
    {
        v[n++] = atoi( p );
        while( *p && *p != ',' ) p++;
        if( *p == ',' ) p++;
    }
    v[n] = -1;
    return v;
}

int main( int argc, char **argv )
{
    int a = 1;
    int nseq = atoi( argv[a++] );
    njob = nseq;

    int ***topol = calloc( nseq, sizeof( int ** ) );
    double **len = calloc( nseq, sizeof( double * ) );
    double **bw = calloc( nseq, sizeof( double * ) );
    for( int i = 0; i < nseq; i++ )
    {
        topol[i] = calloc( 2, sizeof( int * ) );
        len[i] = calloc( 2, sizeof( double ) );
        bw[i] = calloc( 2, sizeof( double ) );
    }
    for( int i = 0; i < nseq - 1; i++ )
    {
        topol[i][0] = parse_members( argv[a++], nseq );
        topol[i][1] = parse_members( argv[a++], nseq );
        len[i][0] = atof( argv[a++] );
        len[i][1] = atof( argv[a++] );
    }

    Node *stopol = calloc( nseq * 2, sizeof( Node ) );
    treeCnv( stopol, nseq, topol, len, bw );
    calcBranchWeight( bw, nseq, stopol, topol, len );

    for( int i = 0; i < nseq - 1; i++ )
        for( int j = 0; j < 2; j++ )
            printf( "bw %d %d %s\n", i, j, bits( bw[i][j] ) );

    /* `len` is clamped in place by treeCnv, so print it back: a port that
     * skipped checkMinusLength would still match every weight on a tree whose
     * branches are all above MINLEN, and differ here. */
    for( int i = 0; i < nseq - 1; i++ )
        for( int j = 0; j < 2; j++ )
            printf( "len %d %d %s\n", i, j, bits( len[i][j] ) );

    double *result = calloc( nseq, sizeof( double ) );
    for( int i = 0; i < nseq - 1; i++ )
    {
        for( int j = 0; j < 2; j++ )
        {
            distFromABranch( nseq, result, stopol, topol, len, i, j );
            for( int k = 0; k < nseq; k++ )
                printf( "dist %d %d %d %s\n", i, j, k, bits( result[k] ) );
            weightFromABranch( nseq, result, stopol, topol, i, j );
            for( int k = 0; k < nseq; k++ )
                printf( "wt %d %d %d %s\n", i, j, k, bits( result[k] ) );
        }
    }
    /* The group split every refinement step starts from. */
    int *memlist1 = calloc( nseq + 1, sizeof( int ) );
    int *memlist2 = calloc( nseq + 1, sizeof( int ) );
    char *pairbuf = calloc( nseq, sizeof( char ) );
    for( int i = 0; i < nseq - 1; i++ )
    {
        for( int j = 0; j < 2; j++ )
        {
            int s1 = -1, s2 = -1;
            OneClusterAndTheOther_fast( nseq, memlist1, memlist2, &s1, &s2, pairbuf, topol, i, j,
                                        NULL, NULL, NULL );
            printf( "split %d %d s %d %d\n", i, j, s1, s2 );
            for( int k = 0; memlist1[k] != -1; k++ ) printf( "split %d %d m1 %d %d\n", i, j, k, memlist1[k] );
            for( int k = 0; memlist2[k] != -1; k++ ) printf( "split %d %d m2 %d %d\n", i, j, k, memlist2[k] );
        }
    }
    return 0;
}
