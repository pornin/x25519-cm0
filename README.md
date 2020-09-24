# X25519 Implementation for ARM Cortex-M0/M0+

This code implements the X25519 function (Montgomery ladder
multiplication in Curve25519) as specified by [RFC
7748](https://tools.ietf.org/html/rfc7748), in ARMv6-M assembly. It is
suitable for running on ARM Cortex-M0 and M0+ CPUs. It is fully
constant-time, without any assumption on the behaviour of the RAM system
(e.g. presence of cache or of wait states in cases of DMA conflict).
The core function is executed in **3229950 clock cycles** (measured on an
Atmel SAM D20 Xplained Pro board, clocked at 8 MHz).

This could should also run on larger 32-bit ARM systems, e.g. ARM
Cortex M3 and M4 (but substantially faster implementations exist
for these larger CPUs).

# Usage

The implementation and test code can be compiled with a simple `make`;
this assumes that the compiler is available in the `PATH` under the
name `arm-linux-gcc`, and is able to link against an appropriate
C library.

On non-ARM hosts, [Buildroot](https://buildroot.org/) can be used to
setup a cross-compilation toolchain, with compiler, linker, and libc,
for a target Linux system. Then, [QEMU](https://www.qemu.org/) user
mode emulation can run the binary (`qemu-arm ./test_x25519`).

However, the whole point of this implementation is to be integrated in
larger code bases. It consists of [a single assembly source
file](src/x25519-cm0.S) (not counting the tests).

# Integration Notes

X25519 is implemented as specified in RFC 7748. Notably:

  - In the input 'u' coordinate, the last bit (most significant bit of
    the last byte) is ignored.

  - All possible 255-bit patterns are supported, including values in
    the 2^255-19 to 2^255-1 range.

  - "Clamping" is automatically applied on the scalar:

      - bit 255 is forced to 0
      - bit 254 is forced to 1
      - bits 0, 1 and 2 are forced to 0

  - All operations are strictly "constant-time". In particular, the
    memory access pattern does not depend on secret values. There is
    no reliance on RAM accesses to be free of timing-based side
    channels (the cumulative cost of the constant-time "conditional
    swap" operations is about 46080 cycles).

The implementation has a single dependency on the external `memcpy()`
function.

The `x25519()` function is the external API; it is callable from C code,
since it respects the standard ABI. It takes three parameters:

```
void x25519(void *dst, const void *src, const void *scalar);
```

All three buffers have size 32 bytes, exactly. `dst` points to the
output buffer that receives the result. `src` is the source point
(encoded 'u' coordinate of a curve point, in RFC 7748 terminology).
`scalar` is the scalar, i.e. the integer by which the source point
is to be multiplied (unsigned little-endian encoding). There are
no alignment constraints on the inputs and outputs (`memcpy()` is
internally used to copy data in and out). Input and output buffers
may freely overlap.

None of the code contained therein modifies or even reads the r9
register. This register may be reserved for all purposes by the ABI
(in some systems using this ABI, r9 may be reserved at all times and
not even usable for local temporary storage; this was the case, for
instance, in iOS up to version 2).

# Internal Implementation Notes

Internal functions do not follow the normal ABI, in that they do not
preserve any register. The caller must save whichever registers it
wishes to conserve. This leads to performance gains because most
internal calls don't have many values to save, and most register saving
actions on the callee side would be lost. However, global functions,
which are potentially callable from the C side, must save the registers
mandated by the ABI (namely registers r4 to r11).

Field elements are internally represented as sequences of eight 32-bit
words, in little-endian order. All 256-bit patterns are allowed; i.e.
values up to and including 2^256-1 are properly handled. The
`gf_normalize_inner()` function performs full reduction and ensures the
output is in the 0..p-1 range.

Implementation uses the macro `MQ` to designate the value 19. All
functions can in fact be used with other moduli 2^255-t for odd values
of t in the 1..32767 range; the `MQ` macro just has to be adjusted
accordingly. The `gf_inv_inner()` and `gf_legendre_inner()` functions have
some extra requirements:

  - Both functions expect the modulus to be a prime integer (i.e. that
    we work in a finite field).

  - `gf_inv_inner()` requires the `INVT510` macro to evaluate to the
    correct constant value of 1/2^510 modulo p. If `MQ` is changed, then
    this value must be adjusted accordingly.

All other functions work just as well with any modulus as long as `MQ`
is odd and in the supported range.

(If we assumed that `MQ` is in the 1..127 range, then a few cycles could
be saved here and there by replacing some 2-cycle `ldr` opcodes with
1-cycle `movs` opcodes. This might save up to 2000 cycles in total, i.e.
less than 0.1% of the total cost. This optimization has not been
implemented here.)

### Modular Inversions

The `gf_inv_inner()` function implements inversion in the field using the
optimized binary GCD algorithm described in:
[https://eprint.iacr.org/2020/972](https://eprint.iacr.org/2020/972)

On the ARM Cortex-M0+, this function is much faster than the
traditional method using Fermat's Little Theorem: the latter requires
254 squarings and about 11 extra multiplications, for a total cost of
at least 270000 cycles, while `gf_inv_inner()` completes in 54793 cycles
only. `gf_inv_inner()` is fully constant-time, like the rest of the code.

`gf_inv_inner()` requires the modulus to be prime only because it assumes
that the GCD will be 1 unless the source value is zero. The reliance on
a prime modulus could be removed by instead performing a multiplication
at the end to verify that an inverse has truly been obtained; the
overhead would be about 1500 to 2000 cycles.

The constant precomputed value 1/2^510 modulo p would be removed by
replacing that multiplication with two Montgomery reductions. This would
also imply an overhead of a few thousand cycles, and need some extra code
for Montgomery reduction.

### Legendre Symbol

The `gf_legendre_inner()` function computes the Legendre symbol for a field
element. This is not actually needed for X25519, and is included here
only because it could be helpful in other operations adjacent to
X25519, e.g. the use of the Elligator2 map for encoding/hashing values
into curve points in a constant-time way. The Legendre symbol of x is:

  - 1 if x is a non-zero quadratic residue in the field;

  - -1 if x is not a quadratic residue in the field;

  - 0 if x is zero.

The traditional method is again Fermat's Little Theorem: for a prime
p, the Legendre symbol of x is equal to x^((p-1)/2) mod p. This would
again require about 270000 cycles for p = 2^255-19.

The algorithm implemented here is roughly the same as the binary GCD
used for inversion. It internally computes the GCD of x and p with the
exact same steps (hence, it always converges with the same number of
iterations); it does not keep track of the Bezout coefficients, since
these are not needed for a Legendre symbol; however, it follows value
updates to compute the symbol. What is actually computed is the
[Kronecker symbol](https://en.wikipedia.org/wiki/Kronecker_symbol)
(x|p), with the following properties:

  - (x|n) is equal to the Legendre symbol of x modulo n when n is a
    nonnegative odd prime.

  - (x|n) == (y|n) if x == y mod n and either n > 0, or x and y have the
    same sign.

  - If n and m are not both negative, then (n|m) == (m|n), unless both n
    == 3 mod 4 and m == 3 mod 4, in which case (n|m) == -(m|n). (This is
    the law of quadratic reciprocity.)

  - (2|n) == 1 if n == 1 or 7 mod 8, or -1 if n == 3 or 5 mod 8.

In the course of the binary GCD algorithm, we work over two values a
and b, such that they both converge toward 0 and 1. b is always odd.
Each iteration consists in three successive steps:

 1. If a and b are odd and a < b, then a and b are exchanged.

 2. If a is odd, then a is replaced with a-b.

 3. a <- a/2

When adapted to the Legendre symbol computation, we use the same steps,
but also maintain the expected Kronecker symbol in a variable j which
is initially 1, and is negated when approriate:

  - Step 1 exercises the law of quadratic reciprocity; j is negated if
    both a and b are equal to 3 modulo 4 at the time of the swap.

  - Step 2 does not change the Kronecker symbol; a critical observation
    here is that throughout the optimized binary GCD algorithm, it can
    never happen that a and b are both negative.

  - Step 3 negates j if and only if b == 3 or 5 mod 8 at that point.

These updates to j only need to look at the low bits of a and b (up to
three bits) and is thus largely compatible with the intermediate values
maintained by the optimized binary GCD in its inner loop. This implies
a relatively low overhead for the inner loop iterations. Combined with
the savings obtained by not keeping track of the Bezout coefficients,
we finally achieve the Legendre symbol computation in 43726 cycles, i.e.
even faster than inversions. This implementation is fully constant-time.

`gf_legendre_inner()` requires the modulus p to be prime only because it
assumes the GCD of x and p to be 1 as long as x != 0. The implementation
could be modified to support a non-prime modulus (in which case this
would compute the Jacobi symbol) with a slight overhead (the final
iterations may not be specialized out, and we would need an extra
comparison of b with 1 to check for a non-invertible input case).
