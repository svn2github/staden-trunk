int main() {
    int i;
    for (i = 0; i < 10000; i++) {
	putchar("ACGT"[random() % 4]);
    }
}
