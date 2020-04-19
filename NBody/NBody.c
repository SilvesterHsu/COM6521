#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <omp.h>

#include "NBody.h"
#include "NBodyVisualiser.h"

#define USER_NAME "elt18sx"		//replace with your username

#define FILE_CACHE_SIZE 255		//cache used to store one line content while reading file

void print_help();
void step(void);

struct argument{
	// Arguments to store essential parameters, such
	// as N, D, M and I, etc.

	unsigned int n;
	unsigned int d;
	enum MODE m;
	unsigned int iter;
	char* input_file;
	boolean visualisation;
};

struct point {
	// An Structure to store float data in both x and y
	// axis. Usually, it is used to store the acceleration

	float x;
	float y;
};

void raise_error(char* error_message,boolean print_msg,int exit_type) {
	/*------------------------------------------------------
	When encounters an error, the scheduled output will be
	displayed and the program will exit.

	Args:
		error_message: Message to output
		print_msg: A flag of whether to output the message
		exit_type: Exit type, for example, 1, 2 or 3

	Return:
		void
	------------------------------------------------------*/

	printf("%s\n",error_message);
	if (print_msg)
		print_help();
	if (exit_type != 0)
		exit(exit_type);
}

struct argument load_args(int argc, char* argv[]) {
	/*------------------------------------------------------
	Check the validity of input parameters from command line

	Args:
		argc: The number of input parameters, and should be
			always equal or larger than 1.
		argv: An array which stores the parameters as string.
			In addition, `argv[0]` is always the PATH of the
			program.
		args: A structure which is used to store a set of
			parameters.

	Return:
		args: An argument struct that stores the arguments.

	Raises:
		exit(1): encounter with unexpected input.
	--------------------------------------------------------*/

	struct argument args = {200,8,CPU,1,NULL,TRUE};

	// Part 1: check the number of argc
	if (argc == 1)
		return args;
	else if (argc != 4 && argc != 6 && argc != 8)
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
	/*------------------------------------------------------
	Read the data from file according to the arguments

	Args:
		args: A structure which is used to store a set of
			parameters.
		f: A pointer to read file.
		buff: An array work as cache to store a line of content
		file_index: The index of data line, and comments are
			already ignored.
		token, tokenremain, delims: char pointer used for
			separate function.
		body_member: A single body data.

	Return:
		void

	Raises:
		exit(1): if file doesn't exist.
		exit(1): if the number of data in file is not equal
			to N.
	--------------------------------------------------------*/

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
	/*------------------------------------------------------
	Generate random data for all bodies.

	Args:
		args: A structure which is used to store a set of
				parameters.
		file_index: The index of data line.
		body_member: A single body data.

	Return:
		void
	--------------------------------------------------------*/

	if (args.m == CPU) {
		for (unsigned int file_index = 0; file_index < args.n;) {
			float* body_member = &bodies[file_index++];
			for (int i = 0; i < 5;) {
				if (i == 0 || i == 1) *(body_member + i++) = (float)rand() / 0x8000;
				else if (i == 2 || i == 3) *(body_member + i++) = 0.0;
				else *(body_member + i++) = 1.0 / args.n;
			}
		}
	}
	else if (args.m == OPENMP) {
		int file_index;
		#pragma omp parallel
		{
			srand(omp_get_thread_num() * 1000);
			#pragma omp for schedule(dynamic,1)
				for (file_index = 0; file_index < args.n;file_index++) {
					float* body_member = &bodies[file_index];
					for (int i = 0; i < 5;) {
						if (i == 0 || i == 1) *(body_member + i++) = (float)rand() / 0x8000;
						else if (i == 2 || i == 3) *(body_member + i++) = 0.0;
						else *(body_member + i++) = 1.0 / args.n;
					}
				}
		}
	}
}

struct point calculate_single_body_acceleration(struct nbody* bodies,int body_index, struct argument args) {
	/*------------------------------------------------------
	For a single body, calculate the total acceleration
	generated from other bodies.

	Args:
		args: A structure which is used to store a set of
				parameters.
		target_bodies: A single body data.
		SOFTENING_square: pre-calculate the square of SOFTING
			to accelerate the program speed.

	Return:
		acceleration: An structure with only 2 members which 
			stores the accelerate in both x and y axis.
	--------------------------------------------------------*/
	boolean inner_loop = FALSE;
	const float G_const = G;
	double SOFTENING_square = (double)SOFTENING * SOFTENING;
	struct point acceleration = { 0,0 };
	struct nbody* target_bodies = bodies + body_index;
	//double tic = omp_get_wtime();
	if (args.m == CPU || args.m == OPENMP && !inner_loop) {
		for (unsigned int i = 0; i < args.n; i++) {
			struct nbody* external_body = bodies + i;
			if (i != body_index) {
				float x_diff = external_body->x - target_bodies->x;
				float y_diff = external_body->y - target_bodies->y;
				//float r = sqrt((double)x_diff * x_diff + (double)y_diff * y_diff);
				double r = (double)x_diff * x_diff + (double)y_diff * y_diff;
				float temp = G_const * external_body->m / (float)(sqrt((r + SOFTENING_square))*(r + SOFTENING_square));
				//float temp = G_const * external_body->m / (float)pow(((double)r + SOFTENING_square), 3.0 / 2);
				acceleration.x += temp * x_diff;
				acceleration.y += temp * y_diff;
			}
		}
	}
	// Code for inner loop
	else if (args.m == OPENMP && inner_loop)
	{
		int i;
		int thread_num = omp_get_num_threads();
		struct point *local_acceleration = (struct point *) malloc (sizeof(struct point) * thread_num);
		#pragma omp parallel for
		// set local acceleration to zero
		for (i = 0; i < thread_num; i++)
			local_acceleration[i] = acceleration;
		#pragma omp barrier
		#pragma omp parallel for firstprivate(body_index) //schedule(dynamic,1)
			for (i = 0; i < args.n; i++) {
				struct nbody* external_body = bodies + i;
				if (i != body_index) {
					float x_diff = external_body->x - target_bodies->x;
					float y_diff = external_body->y - target_bodies->y;
					float r = sqrt((double)x_diff * x_diff + (double)y_diff * y_diff);
					float temp = G_const * external_body->m / (float)pow(((double)r + (double)SOFTENING), 3.0 / 2);
					local_acceleration[omp_get_thread_num()].x += temp * x_diff;
					local_acceleration[omp_get_thread_num()].y += temp * y_diff;
				}
			}
		#pragma barrier
		#pragma omp master
			for (i = 0; i < thread_num; i++) {
				acceleration.x += local_acceleration[i].x;
				acceleration.y += local_acceleration[i].y;
			}
	}
	//double t = omp_get_wtime() - tic;
	return acceleration;
}

void compute_volocity(struct nbody* bodies, float time_step, struct argument args) {
	/*------------------------------------------------------
	Calculate the volocity for bodies according to their
	accelerate.

	Args:
		args: A structure which is used to store a set of
				parameters.
		time_step: The time refers to dt.
		acceleration: An structure with only 2 members which
			stores the accelerate in both x and y axis.

	Return:
		void
	--------------------------------------------------------*/

	//double tic = omp_get_wtime();
	if (args.m == CPU) {
		for (unsigned int i = 0; i < args.n; i++) {
			struct point acceleration = calculate_single_body_acceleration(bodies,i,args);
			(bodies + i)->vx += acceleration.x * time_step;
			(bodies + i)->vy += acceleration.y * time_step;
		}
	}
	else if (args.m == OPENMP) {
		//omp_set_nested(1);
		int i;
		#pragma omp parallel for schedule(dynamic,2)
			for (i = 0; i < args.n; i++) {
				struct point acceleration = calculate_single_body_acceleration(bodies, i, args);
				(bodies + i)->vx += acceleration.x * time_step;
				(bodies + i)->vy += acceleration.y * time_step;
			}
	}
	//double t = omp_get_wtime() - tic;
	//printf("");
}

void update_location(struct nbody* bodies, float time_step, struct argument args) {
	/*------------------------------------------------------
	Calculate the new location for bodies according to their
	present location and speed.

	Args:
		args: A structure which is used to store a set of
				parameters.
		time_step: The time refers to dt.

	Return:
		void
	--------------------------------------------------------*/

	//double tic = omp_get_wtime();
	if (args.m == CPU) {
		for (unsigned int i = 0; i < args.n; i++) {
			(bodies + i)->x += (bodies + i)->vx * time_step;
			(bodies + i)->y += (bodies + i)->vy * time_step;
		}
	}
	else if (args.m == OPENMP) {
		int i;
		#pragma omp parallel for schedule(dynamic,2)
			for (i = 0; i < args.n; i++) {
				(bodies + i)->x += (bodies + i)->vx * time_step;
				(bodies + i)->y += (bodies + i)->vy * time_step;
			}
	}
	//double t = omp_get_wtime() - tic;
	//printf("");
}

void update_heat_map(float* heat_map, struct nbody* bodies, struct argument args) {
	/*------------------------------------------------------
	Calculate the heat map based on present bodies.

	Args:
		args: A structure which is used to store a set of
				parameters.

	Return:
		void
	--------------------------------------------------------*/

	float grid_length = 1.0 / args.d;
	if (args.m == CPU) {
		// Initial heat map
		for (unsigned int i = 0; i < args.d * args.d; i++)
			*(heat_map + i) = 0.0;
		// Iterate over all data points
		for (unsigned int i = 0; i < args.n; i++) {
			struct point body_location = { (bodies+i)->x, (bodies+i)->y };
			if (!(body_location.x < 0 || body_location.x>1 || body_location.y < 0 || body_location.y>1)) {
				int row = (int)(body_location.x / grid_length);
				int line = (int)(body_location.y / grid_length);
				*(heat_map+line*args.d+row) += 1.0;
			}
		}
		// Normalize heat
		for (unsigned int i = 0; i < args.d * args.d; i++)
			*(heat_map + i) = *(heat_map + i) / args.n * args.d;
	}
	else if (args.m == OPENMP) {
		// Initial heat map
		int i;
		#pragma omp parallel for schedule(dynamic,2)
		for (i = 0; i < args.d * args.d; i++)
			*(heat_map + i) = 0.0;
		// Iterate over all data points
		#pragma omp parallel for
		for (i = 0; i < args.n; i++) {
			struct point body_location = { (bodies + i)->x, (bodies + i)->y };
			if (!(body_location.x < 0 || body_location.x>1 || body_location.y < 0 || body_location.y>1)) {
				int row = (int)(body_location.x / grid_length);
				int line = (int)(body_location.y / grid_length);
				//#pragma omp atomic
				#pragma omp critical
				{*(heat_map + line * args.d + row) += 1.0;}
			}
		}
		// Normalize heat
		#pragma omp parallel for schedule(dynamic,2)
		for (i = 0; i < args.d * args.d; i++)
			*(heat_map + i) = *(heat_map + i) / args.n * args.d;
	}
}

struct argument args;
struct nbody* bodies;
float* heat_map;

int main(int argc, char *argv[]) {
	//Processes the command line arguments
		//argc in the count of the command arguments
		//argv is an array (of length argc) of the arguments. The first argument is always the executable name (including path)
	args = load_args(argc,argv);
	if (args.m == OPENMP)
		omp_set_num_threads(omp_get_max_threads());
		//omp_set_num_threads(10);

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
	/*------------------------------------------------------
	A single step to update all bodies.
	1. compute the volocity for all bodies.
	2. update the location for all bodies accoreding to their
	present speed(volocity).

	Args:
		args: A structure which is used to store a set of
				parameters.
		time_step: The time refers to dt.

	Return:
		void
	--------------------------------------------------------*/

	float time_step = dt;
	for (unsigned int i = 0; i < args.iter; i++) {
		compute_volocity(bodies, time_step, args);
		#pragma omp barrier
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
