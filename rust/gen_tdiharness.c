/* Stage one whole run of `TreeDependentIteration` in the real C.
 *
 * Compiled against a built `core/` by verify_tdi.sh; `examples/tdiref.rs`
 * takes the same argument vector and runs the port. The refinement loop is a
 * driver, not a function with a small answer, so what is compared is what it
 * produces: every row of the refined alignment, and the outcome code.
 *
 * The configuration staged is the one `scripts/mafft` runs and the only one
 * the port implements: utree=1, weight=4, score_check=1, nkozo=0, alg='A',
 * use_fft=1, sequential.
 *
 * argv: nseq niter cut constraint fastathreshold
 *       nseq rows
 *       nseq-1 x ( members0 members1 len0 len1 )
 *       nseq-1 x ( skip0 skip1 )
 *       nlh  then nlh x ( i j start1 end1 start2 end2 opt )
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mltaln.h"

extern int TreeDependentIteration( int locnjob, char **name, int nlen[M], char **aseq, char **bseq,
                                   int ***topol, double **len, double **distmtx, int **skipthisbranch,
                                   int alloclen, LocalHom **localhomtable, RNApair ***singlerna,
                                   int nkozo, char *kozoarivec, int ntarget, int *targetmap,
                                   int *targetmapr );
extern void constants( int nseq, char **seq );

#define NOTSPEC 100009

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
    (void)argc;
    int a = 1;
    int nseq = atoi( argv[a++] );
    int niterv = atoi( argv[a++] );
    double cutv = atof( argv[a++] );
    int constraintv = atoi( argv[a++] );
    double fthr = atof( argv[a++] );

    char **seq = calloc( nseq + 1, sizeof( char * ) );
    for( int i = 0; i < nseq; i++ ) seq[i] = argv[a++];
    int lgth = strlen( seq[0] );
    /* `dvtditr` sizes it this way (dvtditr.c: `alloclen = nlenmax * 9`), and the
     * refinement loop can genuinely grow an alignment: sized smaller, the C
     * aborts with LENGTH OVER on cases the real program handles. */
    int alloclen = lgth * 9;

    njob = nseq;
    nlenmax = lgth;
    dorp = 'p';
    ppenalty = -1530; ppenalty_dist = NOTSPEC; ppenalty_OP = NOTSPEC;
    ppenalty_ex = 0; ppenalty_EX = NOTSPEC; poffset = 0;
    pamN = NOTSPEC; kimuraR = NOTSPEC;
    nblosum = 62; scoremtx = 1; fmodel = 0; nwildcard = 0;
    penalty_shift_factor = 100.0;
    fftscore = 1;
    fftWinSize = NOTSPEC;
    fftThreshold = NOTSPEC;
    fftNoAnchStop = 0;
    fftRepeatStop = 0;
    kobetsubunkatsu = 1;
    alg = 'A';
    use_fft = 1;
    utree = 1;
    weight = 4;
    score_check = 1;
    cooling = 1;
    devide = 0;
    scoreout = 0;
    checkC = 0;
    rnakozo = 0;
    rnaprediction = 'n';
    RNAscoremtx = 'n';
    specificityconsideration = 0.0;
    parallelizationstrategy = BAATARI2;
    nthread = 0;
    niter = niterv;
    cut = cutv;
    minimumweight = 0.00001;
    constraint = constraintv;
    fastathreshold = fthr;
    newgapstr = "-";
    outgap = 1;
    legacygapcost = 0;
    trywarp = 0;
    treemethod = 'X';

    constants( njob, seq );

    /* The loop writes progress to `trace` and the current best to `pre`; both
     * are dvtditr's, and neither is compared. */
    trap_g = fopen( "/dev/null", "w" );

    int ***topol = calloc( nseq, sizeof( int ** ) );
    double **len = calloc( nseq, sizeof( double * ) );
    int **skipthisbranch = AllocateIntMtx( nseq, 2 );
    for( int i = 0; i < nseq; i++ )
    {
        topol[i] = calloc( 2, sizeof( int * ) );
        len[i] = calloc( 2, sizeof( double ) );
        topol[i][0] = calloc( nseq + 1, sizeof( int ) );
        topol[i][1] = calloc( nseq + 1, sizeof( int ) );
        topol[i][0][0] = topol[i][1][0] = -1;
        skipthisbranch[i][0] = skipthisbranch[i][1] = 0;
    }
    for( int i = 0; i < nseq - 1; i++ )
    {
        free( topol[i][0] ); free( topol[i][1] );
        topol[i][0] = parse_members( argv[a++], nseq );
        topol[i][1] = parse_members( argv[a++], nseq );
        len[i][0] = atof( argv[a++] );
        len[i][1] = atof( argv[a++] );
    }
    for( int i = 0; i < nseq - 1; i++ )
    {
        skipthisbranch[i][0] = atoi( argv[a++] );
        skipthisbranch[i][1] = atoi( argv[a++] );
    }

    LocalHom **localhomtable = NULL;
    int nlh = atoi( argv[a++] );
    if( constraint )
    {
        localhomtable = (LocalHom **)calloc( njob, sizeof( LocalHom * ) );
        for( int i = 0; i < njob; i++ )
        {
            localhomtable[i] = (LocalHom *)calloc( njob - i, sizeof( LocalHom ) );
            for( int j = 0; j < njob - i; j++ )
            {
                localhomtable[i][j].start1 = -1;
                localhomtable[i][j].end1 = -1;
                localhomtable[i][j].start2 = -1;
                localhomtable[i][j].end2 = -1;
                localhomtable[i][j].opt = -1.0;
                localhomtable[i][j].next = NULL;
                localhomtable[i][j].nokori = 0;
                localhomtable[i][j].extended = -1;
                localhomtable[i][j].korh = 'h';
            }
        }
        for( int k = 0; k < nlh; k++ )
        {
            int i = atoi( argv[a++] ), j = atoi( argv[a++] );
            LocalHom *slot = localhomtable[i] + ( j - i );
            LocalHom *node = slot;
            if( slot->opt != -1.0 )
            {
                /* Append, as readlocalhomtable2_half does for a second segment
                 * between the same pair. */
                while( slot->next ) slot = slot->next;
                node = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
                slot->next = node;
                slot->nokori++;
                node->next = NULL;
                node->nokori = 0;
            }
            node->start1 = atoi( argv[a++] );
            node->end1 = atoi( argv[a++] );
            node->start2 = atoi( argv[a++] );
            node->end2 = atoi( argv[a++] );
            node->opt = atof( argv[a++] );
            node->overlapaa = node->end1 - node->start1 + 1;
            node->extended = -1;
            node->korh = 'h';
        }
    }
    else
    {
        for( int k = 0; k < nlh; k++ ) a += 7;
    }

    char **aseq = AllocateCharMtx( nseq, alloclen + 10 );
    char **bseq = AllocateCharMtx( nseq, alloclen + 10 );
    char **name = AllocateCharMtx( nseq, B );
    int *nlenarr = calloc( nseq, sizeof( int ) );
    for( int i = 0; i < nseq; i++ )
    {
        strcpy( aseq[i], seq[i] );
        snprintf( name[i], B, "=seq%d", i );
        nlenarr[i] = strlen( seq[i] );
    }

    int value = TreeDependentIteration( nseq, name, nlenarr, aseq, bseq, topol, len, NULL,
                                        skipthisbranch, alloclen, localhomtable, NULL,
                                        0, NULL, njob, NULL, NULL );

    printf( "outcome %d\n", value );
    for( int i = 0; i < nseq; i++ ) printf( "%d %s\n", i, bseq[i] );
    return 0;
}
