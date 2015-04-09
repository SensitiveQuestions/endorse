typedef struct {
    double lon;
    double lat;
} Position;

double ArcInRadians (Position from, Position to);

double DistanceInMeters (Position from, Position to);

int main (int argc, char *argv[]);
