#include <HsFFI.h>
#include <stddef.h>

void hs_start(void)
{
    int argc = 3;
    char *argv[] = {
        "cvpy",
        "+RTS",
        "-N",
        "-RTS",
        NULL
    };
    char **pargv = argv;

    hs_init(&argc, &pargv);
}

void hs_stop(void)
{
    hs_exit();
}

