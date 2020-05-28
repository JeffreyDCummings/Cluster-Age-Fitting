#include <math.h>
#include <stdio.h>
#include <string.h>

//This version is the final version and corrects for cases where the triangulation has larger than
//90 degrees relative to the isochrone and for the case of where the minimum actually occurs
//between the central and the more distant of the 3 minimization points.  For example, if the 
//step size between one is much larger than the other.

int main()
{
    static double Vobs, BVobs, Diff[50][50], DiffW[50][50], Difftemp[3], Sen[50];
    double Diffstore, Difftot, DifftotW, BVdiff, Sep, Bm[3]={0.0, 0.0, 0.0};
    double comp1, comp2, Vm[3]={0.0, 0.0, 0.0}, BVm, Vdiff;
    int n, x=0, x2, c, m, m2, ID[50];
    char line[700], file[30];
    FILE *in1, *in2, *in3, *out;
    
//This versions takes MSTO photometry where reddening/extinction corrections are already applied
    puts ("Enter name of the corrected MSTO photometry file:");
    scanf ("%s", file);
    in1 = fopen (file, "r");

    puts ("Enter name of input MIST isochrones file list:");
    scanf ("%s", file);
    in2 = fopen (file, "r");

    puts ("Enter name of output minimization file:");
    scanf ("%s", file);
    out = fopen (file, "w");

    while (fgets(line, 50, in1) != NULL) {	
		sscanf (line, "%d %lf %lf", &ID[x], &Vobs, &BVobs);
		m=0;
    while (fgets(line, 50, in2) != NULL) {	
		sscanf (line, "%s", file);	
		in3 = fopen (file, "r");
		n=1;
		c=0;
    	while (fgets(line, 700, in3) != NULL) {		
			if (n >= 14) {
				Bm[0] = Bm[1];
				Bm[1] = Bm[2];
				Vm[0] = Vm[1];
				Vm[1] = Vm[2];
				sscanf (line, "%*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %lf %lf", &Bm[2], &Vm[2]);
				BVm = Bm[2] - Vm[2];
				BVdiff = BVm - BVobs;
				Vdiff = Vm[2] - Vobs - 0.06;
				Difftemp[0] = Difftemp[1];
				Difftemp[1] = Difftemp[2];
				Difftemp[2] = pow(BVdiff,2) + pow(Vdiff,2);
				if (Difftemp[1] <= Difftemp[2] && Difftemp[1] <= Difftemp[0] && BVm < 0.1 && Vm[2] < 2) {
					Sep = pow(Vm[2] - Vm[1], 2) + pow(Bm[2] - Vm[2] - Bm[1] + Vm[1], 2);
					comp1 = Difftemp[1] + Sep - Difftemp[2];
					comp2 = Difftemp[2] + Sep - Difftemp[1];
					if (comp1 < 0 && (Diff[x][m] > Difftemp[1] || c == 0))
						Diff[x][m] = Difftemp[1];
					if (comp2 < 0 && (Diff[x][m] > Difftemp[2] || c == 0))
						Diff[x][m] = Difftemp[2];
					if (comp1 >= 0 && comp2 >= 0) {
						Diffstore = Difftemp[2] - pow((Difftemp[2] - Difftemp[1] + Sep) / (2*sqrt(Sep)), 2);
					    if (Diff[x][m] > Diffstore || c == 0)
					    	Diff[x][m] = Diffstore; }
					c++;
					Sep = pow(Vm[0] - Vm[1], 2) + pow(Bm[0] - Vm[0] - Bm[1] + Vm[1], 2);
					comp1 = Difftemp[1] + Sep - Difftemp[0];
					comp2 = Difftemp[0] + Sep - Difftemp[1];
					if (comp1 < 0 && Diff[x][m] > Difftemp[1])
						Diff[x][m] = Difftemp[1];
					if (comp2 < 0 && Diff[x][m] > Difftemp[0])
						Diff[x][m] = Difftemp[0];
					if (comp1 >= 0 && comp2 >= 0) {
						Diffstore = Difftemp[0] - pow((Difftemp[0] - Difftemp[1] + Sep) / (2*sqrt(Sep)), 2);
					    if (Diff[x][m] > Diffstore)
					    	Diff[x][m] = Diffstore; }
				}
		    }
		 n++;		     
		 }
		 m++;
	 }
	 x++;
	 rewind(in2);
	 }
     fprintf(out,"Isochrone               %8d", ID[0]);
     x--;
	 for (x2 = 1; x2 < x; x2++)
		 fprintf(out, "%8d ", ID[x2]);
	 fprintf(out,"%8d    TOTAL\n", ID[x]);
	 m2=0;
     while (fgets(line, 50, in2) != NULL) {	
		sscanf (line, "%s", file);	
	    fprintf(out,"%s  ", file);
	    Difftot = 0;
		for (x2=0; x2 < x; x2++) {
			fprintf(out,"%8lf ", Diff[x2][m2]);
			Difftot = Difftot+Diff[x2][m2]; 
			if (m2 > 0)
				Sen[x2] = Sen[x2] + pow(sqrt(Diff[x2][m2]) - sqrt(Diff[x2][m2-1]), 2);
		}
		fprintf(out, "%8lf %8lf\n", Diff[x][m2], Difftot);			
		if (m2 > 0)
			Sen[x] = Sen[x] + pow(sqrt(Diff[x][m2]) - sqrt(Diff[x][m2-1]), 2);
		m2++;
	 }
	 rewind(in2);
	 fprintf(out, "The Following Values Are Quadratically Weighted By A Star's Age Sensitivity\n");
	 m2 = 0;
	 while (fgets(line, 50, in2) != NULL) {	
		sscanf (line, "%s", file);	
	    fprintf(out, "%s  ", file);
	    DifftotW = 0;
		for (x2 = 0; x2 < x; x2++) {
			DiffW[x2][m2] = Diff[x2][m2] * Sen[x2] * 1000000;
			fprintf(out, "%8lf ", DiffW[x2][m2]);
			DifftotW = DifftotW + DiffW[x2][m2]; }
		DiffW[x][m2] = Diff[x][m2] * Sen[x] * 1000000;
		DifftotW = DifftotW + DiffW[x][m2]; 
		fprintf(out, "%8lf %8lf\n", DiffW[x][m2], DifftotW);
		m2++;
	 }

	 fclose (in1);
	 fclose (in2);
	 fclose (in3);
	 fclose (out);
}
