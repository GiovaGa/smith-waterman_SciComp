
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "benchmark.h"
#include "sequence.h"

int smith_waterman(const int N, const int M, const char* restrict A, const char* restrict B, const struct scores_t*scores, int* restrict H);
int smith_waterman_parallel(const int N, const int M, const char* restrict A, const char* restrict B, const struct scores_t*scores, int* restrict H);
int smith_waterman_quadratic(const int N, const int M, const char* restrict A, const char* restrict B, const struct scores_t*scores, int* restrict H);
long smith_waterman_flops(int, int);
long smith_waterman_flops_quadratic(int, int);
int smith_waterman_quadratic_opt(const int N, const int M, const char* restrict A, const char* restrict B, const struct scores_t*scores, int* restrict H);
long smith_waterman_flops_quadratic_opt(int, int);
int smith_waterman_quadratic_parallel(const int N, const int M, const char* restrict A, const char* restrict B, const struct scores_t*scores, int* restrict H);
long smith_waterman_flops_quadratic_parallel(int, int);

void program_usage(const char* program_name)
{
	printf(
		"Usage: %s [sequence1] [sequence2] [scores]\n"
		"\twhere [sequence1] and [sequence2] are either:\n"
		"\t\t-f <file> - read the sequence from the file <file>\n"
		"\t\t-r <length> - generate a random sequence of length <length>\n"
		"\t\t-p - read characters from stdin.\n"
        "\tand [scores] contains ALL these flags (if none are provided, defaults will be used)\n"
		"\t\t-m <match score> - Positive integrer for matching nucleotides.\n"
		"\t\t-x <mismatch score> - Negative integrer for mismatching nucleotides.\n"
		"\t\t-o <gap opening score> - Negative integrer for opening a gap.\n"
		"\t\t-o <gap extension score> - Negative integrer for extending a gap.\n",
		program_name);
}

void bad_usage(const char* program_name)
{
	printf("Bad usage.\n");
	program_usage(program_name);
	exit(-1);
}


struct sequence_input
{
	enum class { FROM_FILE, RANDOM, FROM_STDIN } type;
	union
	{
		const char* filename;
		int random_length;
	};
};

/*
* Parses the command line arguments and returns the offset to the next argument.
* Exits the program on failure.
*/
int parse_input(int argc, char* argv[], int argc_offset, struct sequence_input* seq_in)
{
	if (argc == argc_offset)
	{
		seq_in->type = FROM_STDIN;
		return argc_offset;
	}
	if (strcmp(argv[argc_offset], "-f") == 0)
	{
		if (argc < argc_offset + 1) bad_usage(argv[0]);
		seq_in->type = FROM_FILE;
		seq_in->filename = argv[argc_offset + 1];
		return argc_offset + 2;
	}
	else if (strcmp(argv[argc_offset], "-r") == 0)
	{
		if (argc < argc_offset + 1) bad_usage(argv[0]);
		int len = atoi(argv[argc_offset + 1]);
		if (len <= 0)
		{
			fprintf(stderr, "Invalid length\n");
			exit(-1);
		}
		seq_in->type = RANDOM;
		seq_in->random_length = len;
		return argc_offset + 2;
	}
	else if (strcmp(argv[argc_offset], "-p") == 0)
	{
		seq_in->type = FROM_STDIN;
		return argc_offset + 1;
	}
	bad_usage(argv[0]);
    return -1; // 
}

struct sequence_t get_sequence(struct sequence_input seq_in)
{
	switch (seq_in.type)
	{
		case FROM_FILE:
			return read_sequence_from_file(seq_in.filename);
		case RANDOM:
			return get_random_sequence(seq_in.random_length);
		case FROM_STDIN:
			printf("Enter the sequence:\n");
			return read_sequence_from_stdin();
		default:
			fprintf(stderr, "Invalid sequence type\n");
			exit(-1);
	}
}

int parse_scores(int argc, char *argv[], int argc_offset, struct scores_t *scores) {
	for (int i=0; i<argc; i++) {
	}
    while (argc_offset < argc) {
        if (strcmp(argv[argc_offset], "-m") == 0) {
            if (argc < argc_offset + 1) bad_usage(argv[0]);
            scores->match = atoi(argv[argc_offset + 1]);
            argc_offset += 2;
        }

        else if (strcmp(argv[argc_offset], "-x") == 0) {
            if (argc < argc_offset + 1) bad_usage(argv[0]);
            scores->mismatch = atoi(argv[argc_offset + 1]);
            argc_offset += 2;
        }

        else if (strcmp(argv[argc_offset], "-o") == 0) {
            if (argc < argc_offset + 1) bad_usage(argv[0]);
            scores->gap_opening = atoi(argv[argc_offset + 1]);
            argc_offset += 2;
        }


        else if (strcmp(argv[argc_offset], "-e") == 0) {
            if (argc < argc_offset + 1) bad_usage(argv[0]);
            scores->gap_extension = atoi(argv[argc_offset + 1]);
            argc_offset += 2;
        }
    }
	return argc_offset;
}


int main(int argc, char* argv[])
{
	if (argc == 2 && strcmp(argv[1], "-h") == 0)
	{
		program_usage(argv[0]);
		return 0;
	}
	struct sequence_t A, B;

	struct sequence_input A_input, B_input;
	int argc_offset = 1;

	argc_offset = parse_input(argc, argv, argc_offset, &A_input);
	argc_offset = parse_input(argc, argv, argc_offset, &B_input);

	A = get_sequence(A_input);
	B = get_sequence(B_input);

		struct scores_t score = {4, -2, -3, -1};
    //there are still flags to be consumed
    if (argc > argc_offset) {
        if (argc < argc_offset + 8) {
            fprintf(stderr, "Missing %d arguments\n", argc_offset + 8 - argc);
            bad_usage(argv[0]);
            exit(-1);
        }
        argc_offset = parse_scores(argc, argv, argc_offset, &score);
    }
	if (argc != argc_offset)
		bad_usage(argv[0]);


	//printf("\n\nA (%d) =\n%s\n\nB (%d) =\n%s\n\n", A.length, A.data, B.length, B.data);
	printf("SCORES: -m %d, -x %d, -o %d, -e %d\n", score.match, score.mismatch, score.gap_opening, score.gap_extension);


	struct sw_implementation implementations[] = {
		//{ "Smith-Waterman", smith_waterman, smith_waterman_flops },
		//{ "Smith-Waterman Par", smith_waterman_parallel, smith_waterman_flops },
		{ "Smith-Waterman Quad", smith_waterman_quadratic, smith_waterman_flops_quadratic },
		{"Sw-Quad optimized", smith_waterman_quadratic_opt, smith_waterman_flops_quadratic_opt},
		{"Sw-Quad parallel", smith_waterman_quadratic_parallel, smith_waterman_flops_quadratic_parallel}
	};

	benchmark(implementations, sizeof(implementations) / sizeof(implementations[0]), A, B, &score);

	deallocate_sequence(A);
	deallocate_sequence(B);

	return 0;
}
