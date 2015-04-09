#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "geodist.h"

static const double DEG2RAD = 0.01745329251994329576923690768;
static const double EARTH_RADIUS_METERS = 6372797.56;

double ArcInRadians(Position from, Position to) {
    double latArc  = (from.lat - to.lat) * DEG2RAD;
    double lonArc = (from.lon - to.lon) * DEG2RAD;
    double latH = sin(latArc * 0.5);
    latH *= latH;
    double lontitudeH = sin(lonArc * 0.5);
    lontitudeH *= lontitudeH;
    double tmp = cos(from.lat*DEG2RAD) * cos(to.lat*DEG2RAD);
    return 2.0 * asin(sqrt(latH + tmp*lontitudeH));
}

double DistanceInMeters(Position from, Position to) {
    return EARTH_RADIUS_METERS*ArcInRadians(from, to);
}

/* int main(int argc, char *argv[])  */
/* {  */
/*   int done = 0, n, numprocs, i=0, j, rc;  */
/*   int number=0, start=0; */

/*   Position *points = malloc(10000 * sizeof *points); */
/*   int *myid = malloc(10000 * sizeof(int)); */
/*   double distance; */

/*   FILE *fp; */
/*   fp = fopen("latlng.txt", "r"); */
/*   int index = 0; */
/*   while(fscanf(fp, "%d %lf %lf", &myid[index], &points[index].lon, &points[index].lat) == 3) */
/*       index++; */
/*   fclose(fp);  */
   
/*   Position from, to; */

/*   FILE *fs; */
/*   fs = fopen("dist.txt", "a+"); */

/*   for (i=0; i <= (index-2); i++) { */
/*     from=points[i]; */
/*     for (j=i+1; j<=(index-1); j++) { */
/*         to=points[j]; */
/*         distance=DistanceInMeters(from,to); */
/* 	fprintf(fs, "%d %d %lf %lf %lf %lf %lf\n", myid[i], myid[j], */
/* 		from.lon, from.lat, to.lon, to.lat, distance);  */
/* 	/\* printf("from %d to %d: (%lf,%lf) to (%lf,%lf): %lf meters\n", myid[i], */
/* 	   myid[j], from.lon, from.lat, to.lon, to.lat, distance); *\/ */
/*     } */
/*   }  */
/*   fclose(fs); */

/* }  */
