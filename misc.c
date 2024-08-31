#include <stdio.h>
#include <stdlib.h>
/*
 *	Global variables for rstart & uni
 */

float uni_u[98];	/* Was U(97) in Fortran version -- too lazy to fix */
float uni_c, uni_cd, uni_cm;
int uni_ui, uni_uj;

float uni(void)
{
	float luni;			/* local variable for uni */

	luni = uni_u[uni_ui] - uni_u[uni_uj];
	if (luni < 0.0)
		luni += 1.0;
	uni_u[uni_ui] = luni;
	if (--uni_ui == 0)
		uni_ui = 97;
	if (--uni_uj == 0)
		uni_uj = 97;
	if ((uni_c -= uni_cd) < 0.0)
		uni_c += uni_cm;
	if ((luni -= uni_c) < 0.0)
		luni += 1.0;
	return (float) luni;
}

void rstart(int i, int j, int k, int l)
{
	int ii, jj, m;
	float s, t;

	for (ii = 1; ii <= 97; ii++) {
		s = 0.0;
		t = 0.5;
		for (jj = 1; jj <= 24; jj++) {
			m = ((i*j % 179) * k) % 179;
			i = j;
			j = k;
			k = m;
			l = (53*l+1) % 169;
			if (l*m % 64 >= 32)
				s += t;
			t *= 0.5;
		}
		uni_u[ii] = s;
	}
	uni_c  = 362436.0   / 16777216.0;
	uni_cd = 7654321.0  / 16777216.0;
	uni_cm = 16777213.0 / 16777216.0;
	uni_ui = 97;	/*  There is a bug in the original Fortran version */
	uni_uj = 33;	/*  of UNI -- i and j should be SAVEd in UNI()     */
}


/* ~rinit: this takes a single integer in the range
		0 <= ijkl <= 900 000 000
	and produces the four smaller integers needed for rstart. It is
 *	based on the ideas contained in the RMARIN subroutine in
 *		F. James, "A Review of Pseudorandom Number Generators",
 *			Comp. Phys. Commun. Oct 1990, p.340
 *	To reduce the modifications to the existing code, rinit now
 *	takes the role of a preprocessor for rstart.
 *
 *	This is useful for the parallel version of the code as James
 *	states that any integer ijkl will produce a statistically
 *	independent sequence of random numbers.
 *
 *     Very funny. If that statement was worth anything he would have provided
 *     a proof to go with it. spb 12/12/90 
 */

void rinit(int ijkl)
{
	int i, j, k, l, ij, kl;

	/* check ijkl is within range */
	if( (ijkl < 0) || (ijkl > 900000000) )
		{
		printf("rinit: ijkl = %d -- out of range\n\n", ijkl);
		exit(3);
               	}

/*        printf("rinit: seed_ijkl = %d\n", ijkl); */

	/* decompose the long integer into the the equivalent four
	 * integers for rstart. This should be a 1-1 mapping
 	 *	ijkl <--> (i, j, k, l)
	 * though not quite all of the possible sets of (i, j, k, l)
	 * can be produced.
	 */

	ij = ijkl/30082;
	kl = ijkl - (30082 * ij);

	i = ((ij/177) % 177) + 2;
	j = (ij % 177) + 2;
	k = ((kl/169) % 178) + 1;
	l = kl % 169;

	if( (i <= 0) || (i > 178) )
		{
		printf("rinit: i = %d -- out of range\n\n", i);
		exit(3);
               	}

	if( (j <= 0) || (j > 178) )
		{
		printf("rinit: j = %d -- out of range\n\n", j);
		exit(3);
               	}

	if( (k <= 0) || (k > 178) )
		{
		printf("rinit: k = %d -- out of range\n\n", k);
		exit(3);
               	}

	if( (l < 0) || (l > 168) )
		{
		printf("rinit: l = %d -- out of range\n\n", l);
		exit(3);
               	}

	if (i == 1 && j == 1 && k == 1)
		{
                printf("rinit: 1 1 1 not allowed for 1st 3 seeds\n\n");
		exit(4);
                }

/*        printf("rinit: initialising RNG via rstart(%d, %d, %d, %d)\n",
				i, j, k, l); */

        rstart(i, j, k, l);

}




/*
 *  Function to write a percolation map in black and white Portable
 *  Bit Map (PBM) format.
 *
 *  Note that this version expects the map array to have been
 *  dynamically allocated, e.g. using the arralloc() routine:
 *
 *  int **cell;
 *  cell = (int **) arralloc(sizeof(int), 2, L, L);
 *  ...
 *  writecelldynamic("cell.pbm", cell, L);
 */

void writecelldynamic(char *cellfile, int **cell, int l)
{
  FILE *fp;

  int i, j, npix, col;
  static int pixperline = 32; // the PGM format limits to 70 characters per line

  /*
   *  Write the file
   */

  printf("writecelldynamic: opening file <%s>\n", cellfile);

  fp = fopen(cellfile, "w");

  printf("writecelldynamic: writing data ...\n");

  /*
   *  Start with the PBM header
   */

  fprintf(fp, "P1\n");
  fprintf(fp, "# Written by writecelldynamic\n");
  fprintf(fp, "%d %d\n", l, l) ;

  /*
   *  Now write the cells to file so that cell[0][0] is in the
   *  bottom-left-hand corner and cell[l-1][l-1] is in the
   *  top-right-hand corner
   */

  npix = 0;

  for (j=l-1; j >= 0; j--)
    {
      for (i=0; i < l; i++)
	{
	  npix++;

          // Strangely, PBM files have 1 for black and 0 for white
          
          col = 1;
          if (cell[i][j] == 1) col = 0;

	  // Make sure lines wrap after "npix" pixels

	  if (npix == 1)
	    {
	      fprintf(fp, "%1d", col);
	    }
	  else if (npix < pixperline)
	    {
	      fprintf(fp, " %1d", col);
	    }
	  else
	    {
	      fprintf(fp, " %1d\n", col);
	      npix = 0;
	    }
	}
    }

  if (npix != 0) fprintf(fp, "\n");

  printf("writecelldynamic: ... done\n");

  fclose(fp);
  printf("writecelldynamic: file closed\n");
}