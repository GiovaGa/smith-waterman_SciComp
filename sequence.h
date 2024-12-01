#ifndef SEQUENCE_H
#define SEQUENCE_H

/*
* Structure representing a sequence of characters.
* The sequence is guaranteed to contain only uppercase A, C, G, T and be null-terminated.
* 'length' is the number of characters in the sequence, excluding the null terminator.
*/
struct sequence_t
{
	const char* restrict data;
	int length;
};

/*
* Deallocates memory for a sequence.
*/
void deallocate_sequence(struct sequence_t sequence);

/*
* Generates a random sequence of the specified length.
* Exits the program on failure.
* The caller is responsible for deallocating the sequence.
*/
struct sequence_t get_random_sequence(int length);

/*
* Reads a sequence from a file.
* Converts the sequence to uppercase.
* Exits the program on failure.
* The caller is responsible for deallocating the sequence.
*/
struct sequence_t read_sequence_from_file(const char* filename);


/*
* Reads a sequence from standard input.
* Converts the sequence to uppercase.
* Exits the program on failure.
* The caller is responsible for deallocating the sequence.
*/
struct sequence_t read_sequence_from_stdin(void);



#endif