#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

void x25519(void *dst, const void *src, const void *scalar);

static size_t
hextobin(void *dst, size_t max_len, const char *src)
{
	uint8_t *buf;
	size_t u;
	unsigned acc;
	int z;

	buf = dst;
	u = 0;
	acc = 0;
	z = 0;
	while (*src) {
		int c;

		c = *src ++;
		if (c >= '0' && c <= '9') {
			c -= '0';
		} else if (c >= 'A' && c <= 'F') {
			c -= 'A' - 10;
		} else if (c >= 'a' && c <= 'f') {
			c -= 'a' - 10;
		} else if (c == ' ' || c == ':') {
			continue;
		} else {
			fprintf(stderr, "Invalid hex character: %c\n", c);
			exit(EXIT_FAILURE);
		}
		if (z) {
			if (u >= max_len) {
				fprintf(stderr, "hextobin: overflow\n");
				exit(EXIT_FAILURE);
			}
			buf[u ++] = (uint8_t)((acc << 4) + (unsigned)c);
		} else {
			acc = (unsigned)c;
		}
		z = !z;
	}
	if (z) {
		fprintf(stderr, "hextobin: half final byte\n");
		exit(EXIT_FAILURE);
	}
	return u;
}

#define HEXTOBIN_LEN(dst, len, src)   do { \
		size_t p_hextobin_t = (len); \
		if (hextobin((dst), p_hextobin_t, (src)) != p_hextobin_t) { \
			fprintf(stderr, "Unexpected string length\n"); \
			exit(EXIT_FAILURE); \
		} \
	} while (0)

#define HEXTOBIN(dst, src)   HEXTOBIN_LEN(dst, sizeof (dst), src)

static void
check_equals(const void *a1, const void *a2, size_t len, const char *msg)
{
	const uint8_t *b1, *b2;
	size_t u;

	if (memcmp(a1, a2, len) == 0) {
		return;
	}
	fprintf(stderr, "ERR: %s\n", msg);
	b1 = a1;
	b2 = a2;
	fprintf(stderr, "a1 = ");
	for (u = 0; u < len; u ++) {
		fprintf(stderr, "%02x", b1[u]);
	}
	fprintf(stderr, "\n");
	fprintf(stderr, "a2 = ");
	for (u = 0; u < len; u ++) {
		fprintf(stderr, "%02x", b2[u]);
	}
	fprintf(stderr, "\n");
	exit(EXIT_FAILURE);
}

/*
 * Each triplet is: input scalar, input u coordinate, output u coordinate.
 */
static const char *const KAT_X25519[] = {
	"a546e36bf0527c9d3b16154b82465edd62144c0ac1fc5a18506a2244ba449ac4",
	"e6db6867583030db3594c1a424b15f7c726624ec26b3353b10a903a6d0ab1c4c",
	"c3da55379de9c6908e94ea4df28d084f32eccf03491c71f754b4075577a28552",

	"4b66e9d4d1b4673c5ad22691957d6af5c11b6421e0ea01d42ca4169e7918ba0d",
	"e5210f12786811d3f4b7959d0538ae2c31dbe7106fc03c3efc4cd549c715a493",
	"95cbde9476e8907d7aade45cb4b873f88b595a68799fa152e6f8f7647aac7957",

	NULL
};

void
test_x25519(void)
{
	const char *const *s;
	uint8_t scalar[32], src[32], dst[32], ref[32];
	int i;

	printf("Test x25519: ");
	fflush(stdout);

	s = KAT_X25519;
	while (*s != NULL) {
		HEXTOBIN(scalar, *s ++);
		HEXTOBIN(src, *s ++);
		HEXTOBIN(ref, *s ++);

		x25519(dst, src, scalar);
		check_equals(dst, ref, 32, "KAT");

		printf(".");
		fflush(stdout);
	}

	printf(" ");
	fflush(stdout);

	memset(src, 0, sizeof src);
	src[0] = 0x09;
	memcpy(scalar, src, 32);
	for (i = 1; i <= 1000; i ++) {
		uint8_t tmp[32];

		x25519(tmp, src, scalar);
		memcpy(src, scalar, 32);
		memcpy(scalar, tmp, 32);

		if (i == 1 || i == 1000) {
			HEXTOBIN(ref, i == 1
   ? "422c8e7a6227d7bca1350b3e2bb7279f7897b87bb6854b783c60e80311ae3079"
   : "684cf59ba83309552800ef566f2f4d3c1c3887c49360e3875f2eb94d99532c51");
			check_equals(scalar, ref, 32, "KAT MC");
		}

		if (i % 50 == 0) {
			printf(".");
			fflush(stdout);
		}
	}

	printf(" done.\n");
	fflush(stdout);
}

#if 0
void test_gf_add(uint32_t *d, const uint32_t *a, const uint32_t *b);
int32_t test_gf_legendre(uint32_t *x);

void
test_legendre(void)
{
	int i;
	uint32_t t[8];

	printf("Test Legendre: ");
	fflush(stdout);

	/* 1 is a QR; 2 is not. */
	memset(t, 0, sizeof t);
	t[0] = 1;
	for (i = 0; i < 100000; i ++) {
		int qr, qr2;

		qr = test_gf_legendre(t);
		if ((i & 1) == 0) {
			qr2 = 1;
		} else {
			qr2 = -1;
		}
		if (qr != qr2) {
			fprintf(stderr, "wrong QR status\n");
			printf("0x%08X%08X%08X%08X%08X%08X%08X%08X\n",
				t[7], t[6], t[5], t[4],
				t[3], t[2], t[1], t[0]);
			printf("(got %d, expected %d)\n", qr, qr2);
			exit(EXIT_FAILURE);
		}
		test_gf_add(t, t, t);

		if (i % 1000 == 0) {
			printf(".");
			fflush(stdout);
		}
	}

	printf(" done.\n");
	fflush(stdout);
}
#endif

int
main(void)
{
	test_x25519();
	return 0;
}
