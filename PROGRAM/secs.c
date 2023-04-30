#include <stdio.h>
#include <time.h>
#include <unistd.h>

/* #include <resource.h> */

long int secs_ ()
{
	long int	seconds;

	time ( &seconds );
	return seconds;
}

long int mygetpid_ ()
{
	long int	id;

	id = getpid ( );
	return id;
}
