/*  synzero v0.1: script that generates sequences from PWM */
/*  usage: ./synzero [PWM filename] [number of sequences]  */
/*  written Monday 22 July 2024 by J. Taipale              */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

/* NORMALIZED PWM STRUCTURE */
short int max_width_of_pwm = 30;
short int contacts = 0;

struct normalized_pwm {char *name; char *seed; short int width; long int max_counts; double *information_content; short int *original_position; double *position_score; long int *total_counts_for_column; double **fraction; short int negative_values_allowed;};
short int normalized_pwm_init (struct normalized_pwm *i, char *name, short int width, double initial_value)
{
short int maximum_width = max_width_of_pwm;
short int counter;
short int counter2;
(*i).negative_values_allowed = 0;
(*i).name = malloc(100);
strcpy ((*i).name, name);
(*i).seed = malloc(1000);
strcpy ((*i).seed, "UNKNOWN");
(*i).width = width;
(*i).max_counts = initial_value;
(*i).fraction = malloc(sizeof(double *) * (5 + contacts * 12) + 5);
(*i).information_content = malloc(sizeof(double) * maximum_width + 5);
(*i).position_score = malloc(sizeof(double) * maximum_width + 5);
(*i).original_position = malloc(sizeof(short int) * maximum_width + 5);
(*i).total_counts_for_column = malloc(sizeof(long int) * maximum_width + 5);

for (counter = 0; counter < 5 + contacts * 12; counter++)
{
(*i).fraction[counter] = malloc(sizeof(double) * maximum_width + 5);
for (counter2 = 0; counter2 < maximum_width; counter2++) (*i).fraction[counter][counter2] = initial_value;
}
for (counter2 = 0; counter2 < maximum_width; counter2++)
{
(*i).information_content[counter2] = 0;
(*i).position_score[counter2] = 0;
(*i).original_position[counter2] = counter2;
(*i).total_counts_for_column[counter2] = 0;
}
return(0);
}
short int normalized_pwm_free (struct normalized_pwm *i)
{
short int counter;
free((*i).name);
free((*i).information_content);
free((*i).position_score);
free((*i).total_counts_for_column);
for (counter = 0; counter < 5; counter++) free((*i).fraction[counter]);
free((*i).fraction);
return(0);
}


/* SUBROUTINE THAT RENORMALIZES NORMALIZED PWM (ROWS IN EACH COLUMN ADD TO 1) */
short int Normalize_pwm (struct normalized_pwm *n)
{
    short int counter;
    short int position;
    double total_nucleotides = 0;
    double normalized_value = 0;
    // printf("\nWidth: %d", (*n).width);
    for (position = 0; position < (*n).width; position++)
    {
        for (counter = 0, total_nucleotides = 0; counter < 4; counter++)
        {
        if ((*n).fraction[counter][position] > 0) total_nucleotides += (*n).fraction[counter][position];
        else if ((*n).negative_values_allowed == 1) total_nucleotides += -(*n).fraction[counter][position];
        }
        
        for (counter = 0; counter < 4; counter++)
        {
        normalized_value = ((double) (*n).fraction[counter][position]) / total_nucleotides;
        if ((normalized_value < 0) && ((*n).negative_values_allowed == 0)) normalized_value = 0;
        (*n).fraction[counter][position] = normalized_value;
        }
    }
    return (0);
}


/* SUBROUTINE THAT LOADS A PWM AND NORMALIZES IT */
short int Load_pwm (struct normalized_pwm *p, char *filename, short int normalize)
{
long int counter;
char text1;
short int line = 0;
short int pwm_position = 0;
char *current_string;
current_string = malloc(200);
FILE *pwmfile;
if ((pwmfile = fopen(filename, "r")) == (void *)0) {printf("\nNo File: %s", filename); exit (2);}
for(line = 0; line <= 3 + contacts * 12;)
{
    for(counter = 0; counter < 30; counter++)
    {
        text1 = getc(pwmfile);
        if (text1 == EOF || text1 == '\n' || text1 == '\t')
        {
        current_string[counter] = '\0';
        if (counter > 0 && (current_string[0] == '0' || current_string[0] == '1' || current_string[0] == '2' || current_string[0] == '3' || current_string[0] == '4' || current_string[0] == '5' || current_string[0] == '6' || current_string[0] == '7' || current_string[0] == '8' || current_string[0] == '9' || current_string[0] == ' ' || current_string[0] == '-'))
        {
        (*p).fraction[line][pwm_position] = atof(current_string);
        //printf("\n%f", (*p).fraction[line][pwm_position]);
        pwm_position++;
        }
        if (text1 == '\n' || text1 == EOF) {(*p).width = pwm_position; line++; pwm_position = 0;}
        break;
        }
        current_string[counter]=text1;
        /* printf ("%c", text1); */
    }
}
free (current_string);
if (normalize == 1) Normalize_pwm(p);
// for (line = 0; line < 4 + contacts * 12; line++) {printf("\n"); for (pwm_position = 0; pwm_position < (*p).width; pwm_position++) printf("\t%f", (*p).fraction[line][pwm_position]);}
if (text1 == EOF && line != 3) return(1);
else return (0);
}

/* ########################## MAIN PROGRAM ############################# */
/* ########################## MAIN PROGRAM ############################# */
/* ########################## MAIN PROGRAM ############################# */
/* ########################## MAIN PROGRAM ############################# */
/* ########################## MAIN PROGRAM ############################# */
int main (int argc, char *argv[])

{
    
    srand(time(NULL));
    
    char *DNA = "ACGTN";
    char *dnalc = "acgtn";
    short int Nlength;
    
    char *backgroundPWM_name;
    backgroundPWM_name = malloc(1000);
    strcpy(backgroundPWM_name, argv[1]);
    
    long int number_of_generated_sequences;
    number_of_generated_sequences = atoi(argv[2]);
    
    /* BACKGROUND PWM STRUCTURE */
    struct normalized_pwm qb;
    normalized_pwm_init(&qb, "empty", Nlength * 2, 0);

    Load_pwm (&qb, backgroundPWM_name, 1);
    strcpy(qb.name, backgroundPWM_name);
    Normalize_pwm(&qb);
    Nlength = qb.width+1;
    
    signed short int use_background_position = 0;
    long int current_sequence_position;
    double current_random_number;
    double cutoff;
    short int base;
    short int revcomp;
    short int temp_position;
    
    long int start_positions[Nlength];
    for (temp_position = 0; temp_position < Nlength; temp_position++) start_positions[temp_position] = 0;
    
    for (; number_of_generated_sequences > 1; number_of_generated_sequences--, use_background_position = 0)
    {
        printf("\n");
        for (current_sequence_position = 0; current_sequence_position < Nlength-1; current_sequence_position++)
    {
        current_random_number = (double) rand () / (double) RAND_MAX;
           for(cutoff = 0, base = 0; base < 4; base++)
            {
                cutoff += qb.fraction[base][use_background_position];
                if (current_random_number < cutoff) break;
            }
            use_background_position++;
            if (use_background_position >= qb.width) use_background_position = 0;
            printf("%c",DNA[base]);


    }
    }
    
    printf("\n");
    
    //for(temp_position=0; temp_position < Nlength; temp_position++) printf ("%d: %ld\n", temp_position, start_positions[temp_position]);
    //printf("\n");
}
