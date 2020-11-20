#include <stdio.h>  
#include <stdlib.h>
#include <unistd.h>  
#include <math.h>
#include <getopt.h>
#include <string.h>
#include <limits.h>

#define MAX(a, b) ( (a > b) ? a : b )
#define MIN(a, b) ( (a < b) ? a : b )
#define MAX_SEQ_LEN 8192

                                //   A   C   G   T     (for match/mismatch)
const int SCORING_MATRIX[4][4] = { {  2, -3, -3, -3},   // A
                                   { -3,  2, -3, -3},   // C
                                   { -3, -3,  2, -3},   // G
                                   { -3, -3, -3,  2} }; // T

/*const int SCORING_MATRIX[4][4] = {{  1, -1, -1, -1},   // A
                                   { -1,  1, -1, -1},   // C
                                   { -1, -1,  1, -1},   // G
                                   { -1, -1, -1,  1} }; // T*/
                                   
int seq1[MAX_SEQ_LEN]; // space bound
int seq2[MAX_SEQ_LEN];
int alignment_matrix[MAX_SEQ_LEN+1][MAX_SEQ_LEN+1];

int getMax3(int a, int b, int c)
{
    return ( MAX( MAX(a,b),c));
}

int getMax4(int a, int b, int c, int d)
{
    return ( MAX( MAX(a,b),MAX(c,d)));
}

char decode( int a) 
{
    switch(a)
    {
        case 0: return 'A'; 
        case 1: return 'C'; 
        case 2: return 'G'; 
        case 3: return 'T'; 
        case -1: return '-';
    }
}

int encode( char c) 
{
    switch(c)
    {
        case 'A': return 0; 
        case 'C': return 1; 
        case 'G': return 2; 
        case 'T': return 3; 
    }
}

void generateOutput(int search_type, int seq1len, int seq2len, int starti, int* seq1algn, int* seq2algn, int max)
{
    FILE* output_file;
    switch(search_type)
    {
        case 0: 
        {
            output_file = fopen("global-naiveGap.aln", "w"); 
            break;
        }
        case 1: 
        {
            output_file = fopen("global-affineGap.aln", "w"); 
            break;
        }
        case 2: 
        {
            output_file = fopen("local-naiveGap.aln", "w"); 
            break;
        }
        case 3: 
        {
            output_file = fopen("local-affineGap.aln", "w"); 
            break;
        } 
    }
    if (!output_file) 
    {
        printf("ERROR: File could not open.");
        return 1;
    }
    if(search_type == 2 || search_type == 3) fprintf(output_file, "Score = %d", max);
    else if(search_type == 0 || search_type == 1) fprintf(output_file, "Score = %d", alignment_matrix[seq1len][seq2len]);
    for(int i = starti; i <= seq1len + seq2len; i+=60)
    {
        fprintf(output_file, "\nmy_first_sequence\t");
        for(int j = i; j < i + 60 && j <= seq1len + seq2len; j++)
        {
            if( seq1algn[j] == -1) fprintf(output_file, "-");
            else fprintf(output_file, "%c", decode(seq1algn[j]));
        }
        fprintf(output_file, "\nanother_sequence\t");
        for(int j = i; j < i + 60 && j <= seq1len + seq2len; j++)
        {
            if( seq2algn[j] == -1) fprintf(output_file, "-");
            else fprintf(output_file, "%c", decode(seq2algn[j]));
        }
    }
    fclose(output_file);
}

int getStarti(int search_type, int seq1len, int seq2len, int* seq1algn, int* seq2algn)
{
    int starti = 1; 
    for (int i = seq2len + seq1len; i >= 1; i--) 
    { 
        if ((search_type == 1 && seq1algn[i] == -1 && seq2algn[i] == -1) ||
            (search_type == 3 && (seq1algn[i] == -1 || seq2algn[i] == -1)) ||
            (search_type == 0 && seq1algn[i] == -1 && seq2algn[i] == -1) ||
            (search_type == 2 && (seq1algn[i] == -1 || seq2algn[i] == -1)) ) 
        {
            starti = i + 1;
            break;
        }
    }
    return starti;
}

void fillWithGaps(int y, int x, int seq1pos, int seq2pos, int*seq1algn, int*seq2algn)
{
    // fill with gaps
    while(y > 0)
    {
        if(seq2pos > 0)
        {
            seq2pos--;
            seq2algn[y] = seq2[seq2pos];
        }
        else 
        {
            seq2algn[y] = -1; 
        }
        y--;
    }
    while(x > 0)
    {
        if(seq1pos > 0)
        {
            seq1pos--;
            seq1algn[x] = seq1[seq1pos]; 
        }
        else 
        {
            seq1algn[x] = -1; 
        }
        x--;
    }
}

void affineAlignmentMethod(int search_type, int seq1len, int seq2len, int gapopen_p, int gapext_p)
{
    // init matrices
    int e[seq1len+1][seq2len+1];
    int f[seq1len+1][seq2len+1];
    int g[seq1len+1][seq2len+1];
    int traceback_table[seq1len+1][seq2len+1]; // 1: diagonal, 2: vertical, 3: horizontal 
    memset(e, 0, sizeof e);
    memset(f, 0, sizeof f);
    memset(g, 0, sizeof g);
    memset(alignment_matrix, 0, sizeof alignment_matrix);
    memset(traceback_table, 0, sizeof traceback_table);

    e[0][0] = INT_MIN;
    f[0][0] = INT_MIN;
    g[0][0] = INT_MIN;
    alignment_matrix[0][0] = 0;

    int maxi = 0;
    int maxj = 0;
    int max = 0;
    
    for(int i = 1; i <= seq1len; i++)
    {
        e[i][0] = gapopen_p + i * gapext_p;
        f[i][0] = 0;
        alignment_matrix[i][0] = gapopen_p + i * gapext_p;
    }
    
    for(int i = 1; i <= seq2len; i++)
    {
        e[0][i] = 0;
        f[0][i] = gapopen_p + i * gapext_p;
        alignment_matrix[0][i] = gapopen_p + i * gapext_p;
    }

    for( int i = 1; i <= seq1len; i++)
    {
        for( int j = 1; j <= seq2len; j++)
        {
            e[i][j] = MAX( e[i][j-1] + gapext_p,
                           alignment_matrix[i][j-1] + gapopen_p + gapext_p);
            
            f[i][j] = MAX( f[i-1][j] + gapext_p,
                           alignment_matrix[i-1][j] + gapopen_p + gapext_p);

            // match
            if( seq1[i] == seq2[j])
            {
                g[i][j] = alignment_matrix[i-1][j-1] + SCORING_MATRIX[seq1[i-1]][seq2[j-1]];
                traceback_table[i][j] = 1;
            }
            // non match
            else
            {
                g[i][j] = alignment_matrix[i-1][j-1] + SCORING_MATRIX[seq1[i-1]][seq2[j-1]];
            }

            if (search_type == 1) 
            {
                alignment_matrix[i][j] = getMax3( g[i][j], e[i][j], f[i][j]);
            }
            else if (search_type == 3)
            {
                alignment_matrix[i][j] = getMax4( 0, g[i][j], e[i][j], f[i][j]);
            }
                
            if(alignment_matrix[i][j] == g[i][j]) traceback_table[i][j] = 1;
            else if(alignment_matrix[i][j] == f[i][j]) traceback_table[i][j] = 2;
            else if(alignment_matrix[i][j] == e[i][j]) traceback_table[i][j] = 3;
            
            if( search_type == 3 && alignment_matrix[i][j] > max)
            {
                maxi = i;
                maxj = j;
                max = alignment_matrix[i][j];
            }

             
        }
    }

    // TRACEBACK (uses a traceback table)
    int seq1pos = seq1len;
    int seq2pos = seq2len;
    int seq1algn[seq1len + seq2len + 1];
    int seq2algn[seq1len + seq2len + 1];
    int x = seq1len + seq2len;
    int y = seq1len + seq2len;

    if( search_type == 3)
    {
        //printf("maxi: %d, maxj: %d, max: %d\n",maxi,maxj,max);
        seq1pos = maxi;
        seq2pos = maxj;
    }

    while( (search_type == 1 && !(seq1pos == 0 || seq2pos == 0)) || (search_type == 3 && alignment_matrix[seq1pos][seq2pos] != 0))
    {
        // traceback diagonal
        if( traceback_table[seq1pos][seq2pos] == 1)
        {
            seq1algn[x] = seq1[seq1pos-1];
            seq2algn[y] = seq2[seq2pos-1];
            x -= 1;
            y -= 1;
            seq1pos -= 1;
            seq2pos -= 1;
        }
        // traceback vertical
        else if( traceback_table[seq1pos][seq2pos] == 2)
        {
            seq1algn[x] = seq1[seq1pos-1];
            seq2algn[y] = -1;
            x -= 1;
            y -= 1;
            seq1pos -= 1;
        }
        // traceback horizontal
        else if( traceback_table[seq1pos][seq2pos] == 3)
        {
            seq1algn[x] = -1;
            seq2algn[y] = seq2[seq2pos-1];
            x -= 1;
            y -= 1;
            seq2pos -= 1;
        }
    }

    // fill with gaps
    fillWithGaps( y, x, seq1pos, seq2pos,&seq1algn, &seq2algn);

    // find the start point where alignment starts
    int starti = getStarti(search_type, seq1len, seq2len, &seq1algn, &seq2algn);

    // output alignment result
    generateOutput(search_type, seq1len, seq2len, starti, seq1algn, seq2algn, max);
}

void naiveAlignmentMethod(int search_type, int seq1len, int seq2len, int gapopen_p)
{
    memset(alignment_matrix, 0, sizeof alignment_matrix);
    int maxi = 0;
    int maxj = 0;
    int max = 0;

    // init matrix's first row and col with gap penalties
    // base cases
    for( int i = 0; i <= seq1len; i++)
    {
        if( search_type == 0) alignment_matrix[i][0] = i * gapopen_p;
        else if( search_type == 2) alignment_matrix[i][0] = 0;
    }
    for( int i = 0; i <= seq2len; i++)
    {
        if( search_type == 0) alignment_matrix[0][i] = i * gapopen_p;
        else if( search_type == 2) alignment_matrix[0][i] = 0;
    }
    // calculate maximum score using DP
    for( int i = 1; i <= seq1len; i++)
    {
        for( int j = 1; j <= seq2len; j++)
        {
            // matching
            if( seq1[i-1] == seq2[j-1])
            {
                alignment_matrix[i][j] = alignment_matrix[i-1][j-1] + SCORING_MATRIX[seq1[i-1]][seq2[j-1]];
                
                if( search_type == 2)
                {
                    if( alignment_matrix[i][j] > max)
                    {
                        maxi = i;
                        maxj = j;
                        max = alignment_matrix[i][j];
                    }
                }
            }
            // not matching
            else
            {
                // GLOBAL ALIGNMENT
                if( search_type == 0)
                {
                    alignment_matrix[i][j] = getMax3( alignment_matrix[i-1][j-1] + SCORING_MATRIX[seq1[i-1]][seq2[j-1]], // mismatch
                                                    alignment_matrix[i-1][j] + gapopen_p, // deletion
                                                    alignment_matrix[i][j-1] + gapopen_p); // insertion
                }
                // LOCAL ALIGNMENT
                else if( search_type == 2)
                {
                    alignment_matrix[i][j] = getMax4( 0,
                                                    alignment_matrix[i-1][j-1] + SCORING_MATRIX[seq1[i-1]][seq2[j-1]], // mismatch
                                                    alignment_matrix[i-1][j] + gapopen_p, // deletion
                                                    alignment_matrix[i][j-1] + gapopen_p); // insertion
                    if( alignment_matrix[i][j] > max)
                    {
                        maxi = i;
                        maxj = j;
                        max = alignment_matrix[i][j];
                    }
                }
            }
        }
    }

    // TRACEBACK (does not use a traceback table)

    int seq1pos = seq1len;
    int seq2pos = seq2len;
    int seq1algn[seq1len + seq2len + 1];
    int seq2algn[seq1len + seq2len + 1];
    int x = seq1len + seq2len;
    int y = seq1len + seq2len;

    if( search_type == 2)
    {
        //printf("maxi: %d, maxj: %d, max: %d\n",maxi,maxj,max);
        seq1pos = maxi;
        seq2pos = maxj;
    }
    /*
    for(int i = 0; i <= seq1len; i++)
    {
        for(int j = 0; j<=seq2len; j++)
        {
            printf("%d ",alignment_matrix[i][j]);
        }
        printf("\n");
        
    }
    printf("------------------");
    */
    while( (search_type == 0 && !(seq1pos == 0 || seq2pos == 0)) || (search_type == 2 && alignment_matrix[seq1pos][seq2pos] != 0))
    {
        // match
        if( seq1[seq1pos-1] == seq2[seq2pos-1])
        {
            seq1algn[x] = seq1[seq1pos-1];
            seq2algn[y] = seq2[seq2pos-1];
            x -= 1;
            y -= 1;
            seq1pos -= 1;
            seq2pos -= 1;
        }
        // mismatch
        else if( alignment_matrix[seq1pos - 1][seq2pos - 1] + SCORING_MATRIX[seq1[seq1pos-1]][seq2[seq2pos-1]] == alignment_matrix[seq1pos][seq2pos])
        {
            seq1algn[x] = seq1[seq1pos-1];
            seq2algn[y] = seq2[seq2pos-1];
            x -= 1;
            y -= 1;
            seq1pos -= 1;
            seq2pos -= 1;
        }
        // deletion
        else if( alignment_matrix[seq1pos - 1][seq2pos] + gapopen_p == alignment_matrix[seq1pos][seq2pos])
        {
            seq1algn[x] = seq1[seq1pos-1];
            seq2algn[y] = -1;
            x -= 1;
            y -= 1;
            seq1pos -= 1;
        }
        // insertion
        else if( alignment_matrix[seq1pos][seq2pos - 1] + gapopen_p == alignment_matrix[seq1pos][seq2pos])
        {
            seq1algn[x] = -1;
            seq2algn[y] = seq2[seq2pos-1];
            x -= 1;
            y -= 1;
            seq2pos -= 1;
        }
    }
    
    // fill with gaps
    fillWithGaps( y, x, seq1pos, seq2pos,&seq1algn, &seq2algn);
    // find the start point where alignment starts
    int starti = getStarti(search_type, seq1len, seq2len, &seq1algn, &seq2algn);
    // output alignment result
    generateOutput(search_type, seq1len, seq2len, starti, seq1algn, seq2algn, max);
}

int main(int argc, char **argv) {
    
    // read inputs
    int option;
	FILE* input_fp;
	char* input_fn;
    int gapopen_p = 0;
    int gapext_p = 0;
    int search_type = -1;

    static struct option long_options[] =
    {
        {"mode", required_argument, NULL, 'mode'},
        {"input", required_argument, NULL, 'input'},
        {"gapopen", required_argument, NULL, 'gapopen'},
        {"gapext", required_argument, NULL, 'gapext'},
        {NULL, 0, NULL, 0}
    };

	while(( option = getopt_long( argc, argv, "mode:input:gapopen:gapext:", long_options, NULL)) != -1)
	{
		switch(option)
		{   
			case 'mode':
                if (strcmp(optarg, "global") == 0) {
                    search_type = 0;
                }
                else if (strcmp(optarg, "aglobal") == 0)  {
                    search_type = 1;
                }
                else if (strcmp(optarg, "local") == 0)  {
                    search_type = 2;
                }
                else if (strcmp(optarg, "alocal") == 0)  {
                    search_type = 3;
                }
            	break;
			case 'input':
				input_fn = optarg;
				input_fp = fopen(optarg,"r");
            	break;
            case 'gapopen':
                gapopen_p = atoi(optarg);
                break;
            case 'gapext':
                gapext_p = atoi(optarg);
                break;
		}
	}
    //printf("gapext: %d gapopen: %d mode: %d input: %s\n",gapext_p,gapopen_p,search_type,input_fn);

    // read first sequence
    char c;
    int seq1len = 0;
    int seq2len = 0;
    int flag = 0;

	int i = 0;
	while(c != EOF)
    {   
        while((c = getc(input_fp)) != '\n' && c != EOF)
        {
            if( c == '>')
            {
                while((c = getc(input_fp)) != '\n' && c != EOF);
                flag++;
            }
            else
            {
                if (flag == 2 && c != '\n') 
                {
                    seq2[seq2len] = encode(c);
                    seq2len++;
                }
                else if (flag == 1 && c != '\n') 
                {
                    seq1[seq1len] = encode(c);
                    seq1len++;
                }   
            }
        }
    }
    fclose(input_fp);

    // the alignment method for naive global and local alignment options
    if( search_type == 0 || search_type == 2)
    {
        naiveAlignmentMethod(search_type, seq1len, seq2len, gapopen_p);
    }
    else if( search_type == 1 || search_type == 3)
    {
        affineAlignmentMethod(search_type, seq1len, seq2len, gapopen_p, gapext_p);
    }

    printf("file generated.\n");
    return 0;
}