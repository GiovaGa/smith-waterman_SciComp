#include "sequence.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>
#include <stdbool.h>

/*
* Allocates memory for a sequence of the specified length.
* Adds null terminator at the end.
* Exits the program on failure.
*/
static char* allocate_sequence(int length)
{
	if (length <= 0 || length == INT_MAX)
	{
		fprintf(stderr, "Error: %d is not a valid sequence length\n", length);
		exit(-1);
	}
	char* data = (char*)malloc((length + 1) * sizeof(char));
	if (!data)
	{
		fprintf(stderr, "Cannot allocate memory for sequence\n");
		exit(-1);
	}
	data[length] = '\0';
	return data;
}


void deallocate_sequence(struct sequence_t sequence)
{
	free((char*)sequence.data);
}

/*
* Validates the sequence and converts it to upper case.
* Returns true if the sequence is valid, false otherwise.
*/
static bool validate_sequence(char* data, int length)
{
	for (int i = 0; i < length; i++)
	{
		char c = data[i];
		if (c == 'a' || c == 'g' || c == 'c' || c == 't')
		{
			data[i] = toupper(c);
			continue;
		}
		if (c != 'A' && c != 'G' && c != 'C' && c != 'T')
			return false;
	}
	return true;
}


struct sequence_t get_random_sequence(int length)
{
	char* data = allocate_sequence(length);
	for (int i = 0; i < length; i++)
	{
		switch (rand() % 4)
		{
		case 0: data[i] = 'A'; break;
		case 1: data[i] = 'G'; break;
		case 2: data[i] = 'C'; break;
		case 3: data[i] = 'T'; break;
		}
	}
	struct sequence_t sequence = { data, length };
	return sequence;
}

struct sequence_t read_sequence_from_file(const char* filename)
{
	FILE* file = fopen(filename, "r");
	if (!file)
	{
		perror("Cannot open file");
		exit(-1);
	}
	fseek(file, 0, SEEK_END);
	long file_length = ftell(file);
	rewind(file);
	if (file_length <= 0)
	{
		fprintf(stderr, "File is empty\n");
		exit(-1);
	}
	if (file_length > INT_MAX)
	{
		fprintf(stderr, "File is too big\n");
		exit(-1);
	}
	if (fgetc(file) == '>')
	{
		char c;
		do
			c = fgetc(file);
		while (c != '\n' && c != EOF);
	}
	else
		rewind(file);
	file_length -= ftell(file);
	char* data = allocate_sequence((int)file_length);
	int chars_read = fread(data, sizeof(char), file_length, file);
	if (chars_read != file_length)
	{
		fprintf(stderr, "Cannot read the entire file\n");
		exit(-1);
	}
	fclose(file);
	if (!validate_sequence(data, (int)file_length))
	{
		fprintf(stderr, "Invalid sequence\n");
		exit(-1);
	}
	struct sequence_t sequence = { data, (int)file_length };
	return sequence;
}


struct sequence_t read_sequence_from_stdin(void)
{
	int buffer_size = 128;
	char* data = NULL;
	int length = 0;
	while (true)
	{
		// here seq.data is NULL or has length characters + 1 for the null terminator
		data = (char*)realloc(data, buffer_size * sizeof(char));
		if (!data)
		{
			fprintf(stderr, "Cannot allocate memory for sequence\n");
			exit(-1);
		}
		if (!fgets(data + length, buffer_size - length, stdin))
		{
			fprintf(stderr, "Cannot read sequence\n");
			exit(-1);
		}
		int chars_read = strlen(data + length);
		length += chars_read;
		if (data[length - 1] == '\n' || feof(stdin))
		{
			length--;
			data[length] = '\0';
			break;
		}
		buffer_size *= 2;
	}
	if (!validate_sequence(data, length))
	{
		fprintf(stderr, "Invalid sequence\n");
		exit(-1);
	}
	struct sequence_t sequence = { data, length };
	return sequence;
}
