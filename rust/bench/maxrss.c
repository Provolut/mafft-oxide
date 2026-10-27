/* Exact peak RSS of one command, from wait4's ru_maxrss.
 *
 * Why a C helper and not a script: ru_maxrss of a child includes whatever the
 * parent's address space cost before exec, because fork copies it. Measuring
 * from Python put a floor of 8188 KB under every case -- the interpreter's own
 * RSS -- and five cases read exactly that, which is how the mistake showed.
 * This parent is a few hundred KB, so the floor is out of the way.
 *
 * Polling /proc/<pid>/status for VmHWM, which bench/profile.sh does, has the
 * opposite failure: a case that exits in a millisecond can be gone before the
 * first read, and the poller reports start-up or nothing at all.
 *
 *   maxrss <dir> <cmd> [args...]   -> prints peak RSS in KB, exits with cmd's code
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/wait.h>
#include <sys/resource.h>

int main(int argc, char **argv)
{
	if (argc < 3) { fprintf(stderr, "usage: %s <dir> <cmd> [args...]\n", argv[0]); return 2; }
	pid_t pid = fork();
	if (pid < 0) { perror("fork"); return 2; }
	if (pid == 0) {
		if (chdir(argv[1]) != 0) { perror("chdir"); _exit(127); }
		int fd = open("/dev/null", O_WRONLY);
		if (fd >= 0) { dup2(fd, 1); dup2(fd, 2); }
		execv(argv[2], &argv[2]);
		_exit(127);
	}
	int status; struct rusage ru;
	if (wait4(pid, &status, 0, &ru) < 0) { perror("wait4"); return 2; }
	printf("%ld\n", ru.ru_maxrss);
	return WIFEXITED(status) ? WEXITSTATUS(status) : 1;
}
