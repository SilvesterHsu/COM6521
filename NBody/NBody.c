#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <omp.h>

#include "NBody.h"
#include "NBodyVisualiser.h"

#define USER_NAME "elt18sx"		//replace with your username

#define FILE_CACHE_SIZE 255

void print_help();
void step(void);

// Arguments  N, D, M and I.
struct argument{
	unsigned int n;
	unsigned int d;
	enum MODE m;
	unsigned int iter;
	char* input_file;
	boolean visualisation;
};

struct point {
	float x;
	float y;
};

void raise_error(char* error_message,boolean print_msg,int exit_type) {
	printf("%s\n",error_message);
	if (print_msg)
		print_help();
	if (exit_type != 0)
		exit(exit_type);
}

struct argument load_args(int argc, char* argv[]) {
	/*------------------------------------------------------
	Check the validity of input parameters from command line

	Return:
		An argument struct that stores the arguments. 

	Raises:
		exit(1): encounter with unexpected input.
	--------------------------------------------------------*/

	struct argument args = {200,8,CPU,1,NULL,TRUE};

	// Part 1: check the number of argc
	if (argc != 4 && argc != 6 && argc != 8)
		raise_error("Eorror: Incomplete parameters.", TRUE, 1);

	else{
		// Part 2: check arguments N ,D and Mode validation
		if (!atoi(argv[1]) || !atoi(argv[2]))
			raise_error("Error: N and D should be a number.", TRUE, 1);
		else {
			args.n = atoi(argv[1]); // N
			args.d = atoi(argv[2]); // D
		}
		if (strcmp(argv[3], "CPU") != 0 && strcmp(argv[3], "OPENMP") != 0)
			raise_error("Error: mode should be CPU or OPENMP", TRUE, 1);
		else
			args.m = (strcmp(argv[3], "CPU") == 0) ? CPU : OPENMP; // Mode

		// Part 3: check options
		if (argc == 6 || argc == 8) {
			if (strcmp(argv[4], "-i") == 0) {
				args.iter = atoi(argv[5]);
				args.visualisation = FALSE;
				if (argc == 8 && strcmp(argv[6], "-f") == 0)
					args.input_file = argv[7];
			}
			else if (strcmp(argv[4], "-f") == 0) {
				args.input_file = argv[5];
				if (argc == 8 && strcmp(argv[6], "-i") == 0) {
					args.iter = atoi(argv[7]);
					args.visualisation = FALSE;
				}
			}
		}
	}
	return args;
}

void read_file(struct argument args, struct nbody* bodies) {
	FILE* f = fopen(args.input_file,"r");
	if (f == NULL)
		raise_error("Error: file doesn't exist.", FALSE, 1);
	char buff[FILE_CACHE_SIZE];
	int file_index = 0;
	while (fgets(buff, FILE_CACHE_SIZE, (FILE*)f)) //lines
		if (buff[0] != '#') {
			//printf("%s", buff);
			char* token;
			char delims[] = ",";
			char* tokenremain = buff;
			float* body_member = &bodies[file_index++];
			int i = 0;

			for (token = strtok(buff, delims); i < 5; token = strtok(NULL, delims)) {
				//printf("%s\n", token);
				if (token!=NULL && strlen(token)>1)
					*(body_member + i++) = atof(token);
				else {
					if (i == 0 || i == 1) *(body_member + i++) = (float) rand()/0x8000;
					else if (i == 2 || i == 3) *(body_member + i++) = 0.0;
					else *(body_member + i++) = 1.0/args.n;
				}
			}
		}
	fclose(f);
	if (file_index != args.n)
		raise_error("Error: the data in file is incompatible with argument N.",FALSE,1);
}

void generate_data(struct argument args, struct nbody* bodies) {

	for (unsigned int file_index = 0; file_index < args.n;) {
		float* body_member = &bodies[file_index++];
		for (int i = 0; i < 5;) {
			if (i == 0 || i == 1) *(body_member + i++) = (float)rand() / 0x8000;
			else if (i == 2 || i == 3) *(body_member + i++) = 0.0;
			else *(body_member + i++) = 1.0 / args.n;
		}
	}

}

struct point calculate_single_body_acceleration(struct nbody* bodies,int body_index, struct argument args) {
	const float G_const = G;
	struct point acceleration = { 0,0 };
	struct nbody* target_bodies = bodies + body_index;
	for (unsigned int i = 0; i < args.n; i++) {
		struct nbody* external_body = bodies + i;
		if (i != body_index) {
			float x_diff = external_body->x - target_bodies->x;
			float y_diff = external_body->y - target_bodies->y;
			float r = sqrt((double)x_diff * x_diff + (double)y_diff * y_diff);
			float temp = G_const * external_body->m / (float)pow(((double)r + (double)SOFTENING), 3.0 / 2);
			acceleration.x += temp * x_diff;
			acceleration.y += temp * y_diff;
		}
	}
	return acceleration;
}

void compute_volocity(struct nbody* bodies, float time_step, struct argument args) {
	for (unsigned int i = 0; i < args.n; i++) {
		struct point acceleration = calculate_single_body_acceleration(bodies,i,args);
		(bodies + i)->vx += acceleration.x * time_step;
		(bodies + i)->vy += acceleration.y * time_step;
	}
}

void update_location(struct nbody* bodies, float time_step, struct argument args) {
	for (unsigned int i = 0; i < args.n; i++) {
		(bodies + i)->x += (bodies + i)->vx * time_step;
		(bodies + i)->y += (bodies + i)->vy * time_step;
	}
}

void update_heat_map(float* heat_map, struct nbody* bodies, struct argument args) {
	float grid_length = 1.0 / args.d;
	for (unsigned int i = 0; i < args.d * args.d; i++)
		*(heat_map + i) = 0.0;
	for (unsigned int i = 0; i < args.n; i++) {
		struct point body_location = { (bodies+i)->x, (bodies+i)->y };
		if (!(body_location.x < 0 || body_location.x>1 || body_location.y < 0 || body_location.y>1)) {
			int row = (int)(body_location.x / grid_length);
			int line = (int)(body_location.y / grid_length);
			*(heat_map+line*args.d+row) += 1.0;
		}
	}
	for (unsigned int i = 0; i < args.d * args.d; i++)
		*(heat_map + i) = *(heat_map + i) / args.n * args.d;
	printf("");
}

struct argument args;
struct nbody* bodies;
float* heat_map;

int main(int argc, char *argv[]) {
	//Processes the command line arguments
		//argc in the count of the command arguments
		//argv is an array (of length argc) of the arguments. The first argument is always the executable name (including path)
	args = load_args(argc,argv);

	//Allocate any heap memory
	bodies = (struct nbody*) malloc(sizeof(struct nbody) * args.n);
	heat_map = (float*)malloc(sizeof(float) * args.d * args.d);
	
	//Depending on program arguments, either read initial data from file or generate random data.
	if (args.input_file != NULL)
		read_file(args, bodies);
	else
		generate_data(args, bodies);

	//Depending on program arguments, either configure and start the visualiser or perform a fixed number of simulation steps (then output the timing results).
	//args.visualisation = TRUE;
	if (args.visualisation == TRUE) {
		initViewer(args.n, args.d, args.m, &step);
		setNBodyPositions(bodies);
		setHistogramData(heat_map);
		startVisualisationLoop();
	}
	else {
		clock_t tic = clock();
		//Sleep(123456);
		step();
		clock_t toc = clock();
		int seconds = (toc - tic) / CLOCKS_PER_SEC;
		int milliseconds = (toc - tic - seconds * CLOCKS_PER_SEC);
		printf("Execution time %d seconds %d milliseconds\n", seconds, milliseconds);
	}

	free(bodies);
	free(heat_map);
	return 0;
}

void step(void)
{
	float time_step = dt;
	for (unsigned int i = 0; i < args.iter; i++) {
		compute_volocity(bodies, time_step, args);
		update_location(bodies, time_step, args);
		if (args.visualisation == TRUE)
			update_heat_map(heat_map, bodies, args);
	}
}


void print_help(){
	printf("nbody_%s N D M [-i I] [-i input_file]\n", USER_NAME);

	printf("where:\n");
	printf("\tN                Is the number of bodies to simulate.\n");
	printf("\tD                Is the integer dimension of the activity grid. The Grid has D*D locations.\n");
	printf("\tM                Is the operation mode, either  'CPU' or 'OPENMP'\n");
	printf("\t[-i I]           Optionally specifies the number of simulation iterations 'I' to perform. Specifying no value will use visualisation mode. \n");
	printf("\t[-f input_file]  Optionally specifies an input file with an initial N bodies of data. If not specified random data will be created.\n");
}
