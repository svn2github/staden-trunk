#include "IO.h"

int main(int argc, char **argv) {
    GapIO *io;
    int status;

    io = open_db("thrash", "0", &status, 1, 0);
    printf("io=%p, status=%d\n", io, status);

    thrash_0(io->server);
    thrash_1(io->server);
    thrash_2(io->server);
    thrash(io->server);
    close_db(io);
}
