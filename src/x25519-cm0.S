@ =======================================================================
@ This is an X25519 implementation in ARMv6-M assembly, suitable for
@ ARM Cortex M0 and M0+. It will also run on later 32-bit ARM CPUs,
@ e.g. ARM Cortex M3 and M4 (but substantially faster implementations
@ exist for these larger CPUs).
@
@ This source file uses from C preprocessor directives and macros, and
@ thus should be compiled with a preprocessor pass. With GCC and Clang,
@ this is activated by default when the file name ends with '.S' with
@ an uppercase 'S'.
@
@ X25519 is implemented as specified in RFC 7748. Notably:
@  - In the input 'u' coordinate, the last bit (most significant bit of
@    the last byte) is ignored.
@  - All possible 255-bit patterns are supported, including values in
@    the 2^255-19 to 2^255-1 range.
@  - "Clamping" is automatically applied on the scalar:
@      - bit 255 is forced to 0
@      - bit 254 is forced to 1
@      - bits 0, 1 and 2 are forced to 0
@  - All operations are strictly "constant-time". In particular, the
@    memory access pattern does not depend on secret values. There is
@    no reliance on RAM accesses to be free of timing-based side
@    channels (the cumulative cost of the constant-time "conditional
@    swap" operations is about 46080 cycles).
@
@ The implementation has a single dependency on the external memcpy()
@ function.
@
@ The x25519() function is the external API; it is callable from C code,
@ since it respects the standard ABI.
@ None of the code contained therein modifies or even reads the r9
@ register. This register may be reserved for all purposes by the ABI
@ (in some systems using this ABI, r9 may be reserved at all times and
@ not even usable for local temporary storage; this was the case, for
@ instance, in iOS up to version 2).
@
@ -----------------------------------------------------------------------
@
@ IMPLEMENTATION NOTES:
@ =====================
@
@ Field elements are internally represented as sequences of eight 32-bit
@ words, in little-endian order. All 256-bit patterns are allowed; i.e.
@ values up to and including 2^256-1 are properly handled. The
@ gf_normalize_inner() function performs full reduction and ensures the
@ output is in the 0..p-1 range.
@
@ Implementation uses the macro MQ to designate the value 19. All
@ functions can in fact be used with other moduli 2^255-t for odd values
@ of t in the 1..32767 range; the MQ macro just has to be adjusted
@ accordingly. The gf_inv_inner() and gf_legendre_inner() functions have
@ some extra requirements:
@  - Both functions expect the modulus to be a prime integer (i.e. that
@    we work in a finite field).
@  - gf_inv_inner() requires the INVT510 macro to evaluate to the correct
@    constant value of 1/2^510 modulo p. If MQ is changed, then this
@    value must be adjusted accordingly.
@ All other functions work just as well with any modulus as long as MQ
@ is odd and in the supported range.
@
@ (If we assumed that MQ is in the 1..127 range, then a few cycles could
@ be saved here and there by replacing some 2-cycle 'ldr' opcodes with
@ 1-cycle 'movs' opcodes. This might save up to 2000 cycles in total,
@ i.e. less than 0.1% of the total cost. This optimization has not been
@ implemented here.)
@
@
@ The gf_inv_inner() function implements inversion in the field using the
@ optimized binary GCD algorithm described in:
@    https://eprint.iacr.org/2020/972
@ On the ARM Cortex-M0+, this function is much faster than the
@ traditional method using Fermat's Little Theorem: the latter requires
@ 254 squarings and about 11 extra multiplications, for a total cost of
@ at least 270000 cycles, while gf_inv_inner() completes in 54793 cycles
@ only. gf_inv_inner() is fully constant-time, like the rest of the code.
@
@ gf_inv_inner() requires the modulus to be prime only because it assumes
@ that the GCD will be 1 unless the source value is zero. The reliance on
@ a prime modulus could be removed by instead performing a multiplication
@ at the end to verify that an inverse has truly been obtained; the
@ overhead would be about 1500 to 2000 cycles.
@ The constant precomputed value 1/2^510 modulo p would be removed by
@ replacing that multiplication with two Montgomery reductions. This would
@ also imply an overhead of a few thousand cycles, and need some extra code
@ for Montgomery reduction.
@
@
@ The gf_legendre_inner() function computes the Legendre symbol for a field
@ element. This is not actually needed for X25519, and is included here
@ only because it could be helpful in other operations adjacent to
@ X25519, e.g. the use of the Elligator2 map for encoding/hashing values
@ into curve points in a constant-time way. The Legendre symbol of x is:
@   1   if x is a non-zero quadratic residue in the field
@  -1   if x is not a quadratic residue in the field
@   0   if x is zero
@ The traditional method is again Fermat's Little Theorem: for a prime
@ p, the Legendre symbol of x is equal to x^((p-1)/2) mod p. This would
@ again require about 270000 cycles for p = 2^255-19.
@
@ The algorithm implemented here is roughly the same as the binary GCD
@ used for inversion. It internally computes the GCD of x and p with the
@ exact same steps (hence, it always converges with the same number of
@ iterations); it does not keep track of the Bezout coefficients, since
@ these are not needed for a Legendre symbol; however, it follows value
@ updates to compute the symbol. What is actually computed is the
@ Kronecker symbol (x|p), with the following properties:
@
@   (x|n) is equal to the Legendre symbol of x modulo n when n is a
@   nonnegative odd prime.
@
@   (x|n) == (y|n) if x == y mod n and either n > 0, or x and y have
@   the same sign.
@
@   If n and m are not both negative, then (n|m) == (m|n), unless
@   both n == 3 mod 4 and m == 3 mod 4, in which case (n|m) == -(m|n).
@   (This is the law of quadratic reciprocity.)
@
@   (2|n) == 1 if n == 1 or 7 mod 8, or -1 if n == 3 or 5 mod 8.
@
@ In the course of the binary GCD algorithm, we work over two values a
@ and b, such that they both converge toward 0 and 1. b is always odd.
@ Each iteration consists in three successive steps:
@
@  1. If a and b are odd and a < b, then a and b are exchanged.
@  2. If a is odd, then a is replaced with a-b.
@  3. a <- a/2
@
@ When adapted to the Legendre symbol computation, we use the same steps,
@ but also maintain the expected Kronecker symbol in a variable j which
@ is initially 1, and is negated when approriate:
@
@   - Step 1 exercises the law of quadratic reciprocity; j is negated if
@     both a and b are equal to 3 modulo 4 at the time of the swap.
@
@   - Step 2 does not change the Kronecker symbol; a critical observation
@     here is that throughout the optimized binary GCD algorithm, it can
@     never happen that a and b are both negative.
@
@   - Step 3 negates j if and only if b == 3 or 5 mod 8 at that point.
@
@ These updates to j only need to look at the low bits of a and b (up to
@ three bits) and is thus largely compatible with the intermediate values
@ maintained by the optimized binary GCD in its inner loop. This implies
@ a relatively low overhead for the inner loop iterations. Combined with
@ the savings obtained by not keeping track of the Bezout coefficients,
@ we finally achieve the Legendre symbol computation in 43726 cycles, i.e.
@ even faster than inversions. This implementation is fully constant-time.
@
@ gf_legendre_inner() requires the modulus p to be prime only because it
@ assumes the GCD of x and p to be 1 as long as x != 0. The implementation
@ could be modified to support a non-prime modulus (in which case this
@ would compute the Jacobi symbol) with a slight overhead (the final
@ iterations may not be specialized out, and we would need an extra
@ comparison of b with 1 to check for a non-invertible input case).
@
@ =======================================================================

	.syntax	unified
	.cpu	cortex-m0
	.file	"x25519-cm0.S"
	.text

@ =======================================================================
@ All 'inner' functions are local to this file and may use non-conformant
@ ABIs (about preservation and modification of registers). All global
@ functions (callable from the outside) conform to the ABI.
@ =======================================================================

@ We work in finite field 2^255-MQ, with 1 <= MQ <= 32767.
#define MQ         19

@ Constant 1/2^510 mod p, in little-endian order.
#define INVT510    0xE72E181B, 0x4A75B7AA, 0x209ED8FF, 0x9E237502, 0x2595A0F9, 0x8F3F1D13, 0x5242A8C6, 0x093805AC

@ =======================================================================
@ void gf_add_inner(gf *d, const gf *a, const gf *b)
@
@ Registers r0 to r8 and r10 to r12 are modified.
@
@ Cost: 71
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	gf_add_inner, %function
gf_add_inner:
	@ Load and perform main addition.
	ldm	r1!, { r3, r4 }
	ldm	r2!, { r5, r6 }
	adds	r3, r5
	adcs	r4, r6
	mov	r8, r3
	mov	r10, r4
	ldm	r1!, { r3, r4 }
	ldm	r2!, { r5, r6 }
	adcs	r3, r5
	adcs	r4, r6
	mov	r11, r3
	mov	r12, r4
	ldm	r1!, { r3, r4 }
	ldm	r2!, { r5, r6 }
	adcs	r3, r5
	adcs	r4, r6
	ldm	r1!, { r5, r6 }
	ldm	r2, { r1, r2 }
	adcs	r5, r1
	adcs	r6, r2

	@ Sum is now in r8:r10:r11:r12:r3:r4:r5:r6, with carry in C flag.
	@ We recover the carry and combine it with the top bit of the top
	@ word to get the result modulo 2^255 (in r7). We also truncate
	@ the top word.
	sbcs	r1, r1
	adds	r1, #1
	lsls	r1, r1, #1
	lsrs	r7, r6, #31
	adds	r7, r1
	lsls	r6, r6, #1
	lsrs	r6, r6, #1

	@ Using 2^255 = MQ mod p, we perform the reduction. Since we
	@ cleared the top bit of the top word, this cannot create an
	@ extra carry.
	ldr	r1, const_gf_add_mq
	muls	r7, r1
	mov	r1, r8
	mov	r2, r10
	adds	r1, r7
	eors	r7, r7
	adcs	r2, r7
	stm	r0!, { r1, r2 }
	mov	r1, r11
	mov	r2, r12
	adcs	r1, r7
	adcs	r2, r7
	adcs	r3, r7
	adcs	r4, r7
	adcs	r5, r7
	adcs	r6, r7
	stm	r0!, { r1, r2, r3, r4, r5, r6 }

	bx	lr
	.align	2
const_gf_add_mq:
	.long	MQ
	.size	gf_add_inner, .-gf_add_inner

@ =======================================================================
@ void gf_sub_inner(gf *d, const gf *a, const gf *b)
@
@ Registers r0 r1 to r8 and r10 to r12 are modified.
@
@ Cost: 72
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	gf_sub_inner, %function
gf_sub_inner:
	@ Load and perform main subtraction.
	ldm	r1!, { r3, r4 }
	ldm	r2!, { r5, r6 }
	subs	r3, r5
	sbcs	r4, r6
	mov	r8, r3
	mov	r10, r4
	ldm	r1!, { r3, r4 }
	ldm	r2!, { r5, r6 }
	sbcs	r3, r5
	sbcs	r4, r6
	mov	r11, r3
	mov	r12, r4
	ldm	r1!, { r3, r4 }
	ldm	r2!, { r5, r6 }
	sbcs	r3, r5
	sbcs	r4, r6
	ldm	r1!, { r5, r6 }
	ldm	r2, { r1, r2 }
	sbcs	r5, r1
	sbcs	r6, r2

	@ Difference is now in r8:r10:r11:r12:r3:r4:r5:r6, with borrow in
	@ C flag.
	@
	@ Suppose that we recover the borrow (of value -1 or 0) and combine
	@ it with the top bit of the value into a value x (with -2 <= x <= 1);
	@ then, we truncate the remaining value to 255 bits (the top bit
	@ of the top word, currently in r6, is set to 0).
	@ We then need to add x*MQ to the trunacted value (which is in
	@ the 0..2^255-1 range). To get a value that fits in 256 bits, we
	@ also add an extra p = 2^255 - MQ; i.e. in total we add
	@ 2^255 + (x-1)*MQ. Since x-1 <= 0, the result cannot exceed 2^256-1.
	@ This is equivalent to adding 2^255 and then adding (x-1)*MQ. The
	@ addition of 2^255 to the truncated value is equivalent to forcing
	@ its top bit to 1.

	@ Get the borrow and top bit into value x. We in fact recover -(x-1)
	@ into r7.
	sbcs	r1, r1
	lsls	r1, r1, #1
	asrs	r7, r6, #31
	subs	r7, r1
	adds	r7, #1

	@ Force top bit of top word to 1.
	movs	r2, #1
	lsls	r2, r2, #31
	orrs	r6, r2

	@ Add (x-1)*MQ, i.e. subtract -(x-1)*MQ. We have -(x-1) in r7,
	@ and it is nonnegative.
	ldr	r1, const_gf_sub_mq
	muls	r7, r1
	mov	r1, r8
	mov	r2, r10
	subs	r1, r7
	eors	r7, r7
	sbcs	r2, r7
	stm	r0!, { r1, r2 }
	mov	r1, r11
	mov	r2, r12
	sbcs	r1, r7
	sbcs	r2, r7
	sbcs	r3, r7
	sbcs	r4, r7
	sbcs	r5, r7
	sbcs	r6, r7
	stm	r0!, { r1, r2, r3, r4, r5, r6 }

	bx	lr
	.align	2
const_gf_sub_mq:
	.long	MQ
	.size	gf_sub_inner, .-gf_sub_inner

@ =======================================================================

@ Multiply rl (32 bits) by rb (32 bits) and store result in rl:rh.
@ rb is consumed. On output, CF = 0.
@ ru and rt are two scratch registers.
@ Cost: 17
.macro MUL32x32  rl, rh, rb, rt, ru
	@ hi(a)*hi(b) -> rh
	lsrs	\rt, \rb, #16
	lsrs	\rh, \rl, #16
	muls	\rh, \rt

	@ lo(a)*lo(b) -> rl
	@ lo(a)*hi(b) -> rt
	@ hi(a)*lo(b) -> ru
	uxth	\rb, \rb
	lsrs	\ru, \rl, #16
	uxth	\rl, \rl
	muls	\rt, \rl
	muls	\ru, \rb
	muls	\rl, \rb

	@ Add rt and ru at the right place. We can use rb as scratch.
	lsrs	\rb, \rt, #16
	lsls	\rt, \rt, #16
	adds	\rl, \rt
	adcs	\rh, \rb
	lsrs	\rb, \ru, #16
	lsls	\ru, \ru, #16
	adds	\rl, \ru
	adcs	\rh, \rb
.endm

@ Perform 64x64->128 multiplication.
@   ra0:ra1  first operand
@   rb0:rb1  second operand
@ Result is written in ra0:rb0:ra1:rb1
@ rt1, rt2, rt3, rt4, rh5 and rh6 are scratch registers; rh5 and rh6 are
@ high registers (r8+).
@ Cost: 30 + 3*cost(MUL32x32) = 81
.macro MUL64x64  ra0, ra1, rb0, rb1, rt1, rt2, rt3, rt4, rh5, rh6
	@ Save rb0 and rb1 into high registers.
	mov	\rh5, \rb0
	mov	\rh6, \rb1

	@ rt1 <- |a0 - a1|, sign in rt3 (-1 for _positive_, 0 otherwise)
	subs	\rt1, \ra1, \ra0
	sbcs	\rt3, \rt3
	eors	\rt1, \rt3
	subs	\rt1, \rt3

	@ rt4 <- |b0 - b1|, sign in ra0 (-1 for negative, 0 otherwise)
	subs	\rt4, \rb0, \rb1
	sbcs	\rb1, \rb1
	eors	\rt4, \rb1
	subs	\rt4, \rb1

	@ rt3 <- -sign(a0-a1)*sign(b0-b1) (-1 for same sign, 0 otherwise)
	eors	\rt3, \rb1

	@ rb0 and rb1 are now scratch registers.

	@ rt1:rt2:rt3 <- xee = -sign(a0-a1)*sign(b0-b1) XOR |a0-a1|*|b0-b1|
	MUL32x32  \rt1, \rt2, \rt4, \rb0, \rb1
	eors	\rt1, \rt3
	eors	\rt2, \rt3

	@ a0*b0 -> ra0:rb0. Also save lo(xee) into rh5.
	mov	\rt4, \rh5
	mov	\rh5, \rt1
	MUL32x32  \ra0, \rb0, \rt4, \rt1, \rb1

	@ a1*b1 -> ra1:rb1. Also save hi(xee) into rh6.
	mov	\rt4, \rh6
	mov	\rh6, \rt2
	MUL32x32  \ra1, \rb1, \rt4, \rt1, \rt2

	@ We have a0*b0 and a1*b1 in place. We must still:
	@  - add 2^16*a0*b0
	@  - add 2^16*a1*b1
	@  - add ee (with sign)
	@ Adding ee is really adding xee - sign(ee), with:
	@    xee = rh5:rh6:rt3
	@    sign(ee) = rt3:rt3:rt3
	@ Thus, we need to add rh5:rh6:rt3, and subtract rt3:rt3:rt3.
	@ Since rt3 = 0 or -1, subtraction of rt3:rt3:rt3 really is
	@ addition of 0 or 1, which we can smuggle as an extra carry
	@ when adding xee.

	@ Recover xee into rt1:rt2; sign is still in rt3.
	mov	\rt1, \rh5
	mov	\rt2, \rh6

	@ Set rt4 to 0.
	eors	\rt4, \rt4

	@ Add a0*b0 into rt1:rt2:rt3. We also put the extra carry there.
	asrs	\rt3, \rt3, #1
	adcs	\rt1, \ra0
	adcs	\rt2, \rb0
	adcs	\rt3, \rt4

	@ Add a1*b1 into r1:rt2:rt3.
	adds	\rt1, \ra1
	adcs	\rt2, \rb1
	adcs	\rt3, \rt4

	@ Add rt1:rt2:rt3 into the result.
	adds	\rb0, \rt1
	adcs	\ra1, \rt2
	adcs	\rb1, \rt3
.endm

@ Input:
@   id    output buffer index (in stack)
@   alt   if zero, then:
@            r1 = pointer to first integer (a)
@            r2 = pointer to second integer (b)
@         else:
@            first integer (a) is at stack index iin
@            b[0]..b[3] is in r0:r1:r2:r3
@            second integer (b) is not on the stack
@         alt must be 0 or 1
@   iin   index for inputs, if alt != 0 (ignored if alt == 0)
@
@ Output: result goes to index itmp (8 words). If alt != 0, then the
@ eight output words are modified, but only the first two contain the
@ correct values; the remaining 6 words of output are returned in r1..r6.
@
@ All registers except r11 are consumed (including operands).
@ IMPORTANT: destination MUST NOT overlap with either source.
@ Cost:
@   alt = 0: 114 + 3*cost(MUL64x64) = 357
@   alt = 1: 103 + 3*cost(MUL64x64) = 346
.macro MUL128x128  id, alt, iin
	@ If alt == 0, then we save input pointers in r12 and r14.
	@ If alt != 0, then we save b0 on the stack and b1 (b[2]:b[3])
	@ in r12:r14; in that case, stack indexes are skewed until we
	@ pop b0 again.

	@ a0*b0 -> dest
	.if (\alt) != 0
	push	{ r0, r1 }
	mov	r12, r2
	mov	r14, r3
	movs	r2, r1
	add	r1, sp, #(4 * ((\iin) + 2))
	ldm	r1, { r1, r3 }
	MUL64x64  r0, r2, r1, r3, r4, r5, r6, r7, r8, r10
	add	r6, sp, #(4 * ((\id) + 2))
	stm	r6!, { r0, r1, r2, r3 }
	.else
	ldm	r1!, { r4, r6 }
	mov	r12, r1
	ldm	r2!, { r5, r7 }
	mov	r14, r2
	MUL64x64  r4, r6, r5, r7, r0, r1, r2, r3, r8, r10
	add	r0, sp, #(4 * (\id))
	stm	r0!, { r4, r5, r6, r7 }
	.endif

	@ a1*b1 -> dest + 4
	.if (\alt) != 0
	add	r0, sp, #(4 * ((\iin) + 2 + 2))
	ldm	r0, { r0, r2 }
	mov	r1, r12
	mov	r3, r14
	.else
	mov	r0, r12
	ldm	r0, { r0, r2 }
	mov	r1, r14
	ldm	r1, { r1, r3 }
	.endif
	MUL64x64  r0, r2, r1, r3, r4, r5, r6, r7, r8, r10
	.if (\alt) != 0
	add	r4, sp, #(4 * ((\id) + 4 + 2))
	.else
	add	r4, sp, #(4 * ((\id) + 4))
	.endif
	stm	r4!, { r0, r1, r2, r3 }

	@ |a1-a0| -> r2:r4, sign in r6
	.if (\alt) != 0
	add	r0, sp, #(4 * ((\iin) + 2))
	.else
	mov	r0, r12
	subs	r0, #8
	.endif
	ldm	r0, { r0, r1, r2, r4 }
	subs	r2, r0
	sbcs	r4, r1
	sbcs	r6, r6
	eors	r2, r6
	eors	r4, r6
	subs	r2, r6
	sbcs	r4, r6

	@ |b0-b1| -> r1:r3, sign in r0
	.if (\alt) != 0
	pop	{ r1, r3 }   @ This removes the skew on stack indexes
	mov	r5, r12
	mov	r7, r14
	.else
	mov	r1, r14
	subs	r1, #8
	ldm	r1, { r1, r3, r5, r7 }
	.endif
	subs	r1, r5
	sbcs	r3, r7
	sbcs	r0, r0
	eors	r1, r0
	eors	r3, r0
	subs	r1, r0
	sbcs	r3, r0

	@ -sign(a0-a1)*sign(b0-b1) -> r12 (-1 if same sign, 0 otherwise)
	eors	r6, r0
	mov	r12, r6

	@ -sign(a0-a1)*sign(b0-b1) XOR |a0-a1|*|b0-b1| -> r1:r2:r4:r5
	@ Also retrieve -sign(a0-a1)*sign(b0-b1) into r5.
	MUL64x64  r1, r3, r2, r4, r0, r5, r6, r7, r8, r10
	mov	r5, r12
	eors	r1, r5
	eors	r2, r5
	eors	r3, r5
	eors	r4, r5

	@ Final assembly. The destination array already contains a0*b0
	@ and a1*b1. We must add into (shifted by two words):
	@  - a0*b0
	@  - a1*b1
	@  - xee:r12
	@ and also subtract r12:r12:r12:r12:r12. Since r12 = 0 or -1, this
	@ can be done by adding 0 or 1, which is only an initial carry
	@ in one of the additions above.

	@ We will accumulate the intermediate result in r1:r2:r3:r4:r5.

	@ Add a0*b0 into the accumulator; also perform subtraction of
	@ r12:r12:r12:r12:r12 here.
	add	r0, sp, #(4 * (\id))
	ldm	r0!, { r6, r7 }
	asrs	r5, r5, #1
	adcs	r1, r6
	adcs	r2, r7
	ldm	r0!, { r6, r7 }
	adcs	r3, r6
	adcs	r4, r7
	eors	r7, r7
	adcs	r5, r7

	@ Add a1*b1 into the accumulator.
	ldm	r0!, { r6, r7 }
	adds	r1, r6
	adcs	r2, r7
	ldm	r0!, { r6, r7 }
	adcs	r3, r6
	adcs	r4, r7
	eors	r7, r7
	adcs	r5, r7

	@ Add the accumulator into the result, at the right position.
	@ At that point, the destination buffer contains:
	@   a0*b0 + 2^128*a1*b1
	@ The mathematical result is (a0+2^64*a1)*(b0+2^64*b1); therefore,
	@ the accumulator contains:
	@   a0*b1 + a1*b0
	@ which is necessarily positive. The sign extension of that value
	@ into 6 words is then done by appending a word of value zero.

	subs	r0, #24
	ldm	r0!, { r6, r7 }
	adds	r1, r6
	adcs	r2, r7
	ldm	r0!, { r6, r7 }
	adcs	r3, r6
	adcs	r4, r7
	ldm	r0!, { r6, r7 }
	adcs	r5, r6
	eors	r6, r6
	adcs	r6, r7

	@ We write out the final accumulator only in non-alt mode; otherwise,
	@ the caller will use the registers directly.
	.if (\alt) == 0
	subs	r0, #24
	stm	r0!, { r1, r2, r3, r4, r5, r6 }
	.endif
.endm

@ Input:
@   id     output buffer index (in stack)
@   r1     pointer to first integer (a)
@   r2     pointer to second integer (b)
@   itmp   index in stack for temporary area (12 words) (counted in words)
@ IMPORTANT: destination MUST NOT overlap with either source.
@ Cost: 222 + 2*cost(MUL128x128[alt=0]) + cost(MUL128x128[alt=1]) = 1282
.macro MUL256x256  id, itmp
	@ Save b pointer to stack. This means that \id and \itmp have to be
	@ adjusted whenever used (since we pushed 1 extra word), until we
	@ pop the value again.
	push	{ r2 }
	mov	r11, r1

	@ a0*b0 -> dest
	MUL128x128  ((\id) + 1), 0, 0

	@ a1*b1 -> dest + 8
	mov	r1, r11
	ldr	r2, [sp]
	adds	r1, #16
	adds	r2, #16
	MUL128x128  ((\id) + 9), 0, 0

	@ |a1-a0| -> tmp, save sign into r11.
	mov	r1, r11
	ldm	r1, { r0, r1, r2, r3, r4, r5, r6, r7 }
	subs	r4, r0
	sbcs	r5, r1
	sbcs	r6, r2
	sbcs	r7, r3
	sbcs	r0, r0
	eors	r4, r0
	eors	r5, r0
	eors	r6, r0
	eors	r7, r0
	subs	r4, r0
	sbcs	r5, r0
	sbcs	r6, r0
	sbcs	r7, r0
	mov	r11, r0
	add	r1, sp, #(4 * ((\itmp) + 1))
	stm	r1!, { r4, r5, r6, r7 }

	@ |b0-b1| -> r0..r3, and combine sign with the one in stack slot.
	pop	{ r2 }   @ This removes the skew on stack indexes
	ldm	r2, { r0, r1, r2, r3, r4, r5, r6, r7 }
	subs	r0, r4
	sbcs	r1, r5
	sbcs	r2, r6
	sbcs	r3, r7
	sbcs	r7, r7
	eors	r0, r7
	eors	r1, r7
	eors	r2, r7
	eors	r3, r7
	subs	r0, r7
	sbcs	r1, r7
	sbcs	r2, r7
	sbcs	r3, r7
	mov	r6, r11
	eors	r7, r6
	mov	r11, r7
	@ We don't store |b0-b1| on the stack, MUL128x128 can reuse it
	@ from the registers themselves.

	@ |a0-a1|*|b0-b1| -> tmp + 4
	MUL128x128  ((\itmp) + 4), 1, (\itmp)

	@ Load |a0-a1|*|b0-b1| and XOR with sign word.
	@ Words 2..7 of the value are already in r1..r6.
	mov	r7, r11
	eors	r3, r7
	eors	r4, r7
	eors	r5, r7
	eors	r6, r7
	mov	r10, r3
	mov	r11, r4
	mov	r12, r5
	mov	r14, r6
	add	r4, sp, #(4 * ((\itmp) + 4))
	ldm	r4!, { r0, r3 }
	eors	r0, r7
	eors	r1, r7
	eors	r2, r7
	eors	r3, r7

	@ Accumulator: r0:r3:r1:r2:r10:r11:r12:r14:r7

	@ Add a0*b0 to accumulator. Also handle extra carry (for the
	@ subtraction of the sign, which is 0 or -1).
	add	r4, sp, #(4 * (\id))
	asrs	r7, r7, #1
	mov	r8, r7
	ldm	r4!, { r5, r6, r7 }
	adcs	r0, r5
	adcs	r3, r6
	adcs	r1, r7
	mov	r5, r10
	mov	r6, r11
	mov	r7, r12
	mov	r10, r0
	mov	r11, r3
	mov	r12, r1
	@ Accumulator: r10:r11:r12:r2:r5:r6:r7:r14:r8
	ldm	r4!, { r0, r1, r3 }
	adcs	r2, r0
	adcs	r5, r1
	adcs	r6, r3
	mov	r1, r14
	mov	r0, r8
	mov	r14, r2
	@ Accumulator: r10:r11:r12:r14:r5:r6:r7:r1:r0
	ldm	r4!, { r2, r3 }
	adcs	r7, r2
	adcs	r1, r3
	@ eors does not change the carry flag.
	eors	r2, r2
	adcs	r0, r2

	@ Accumulator: r10:r11:r12:r14:r5:r6:r7:r1:r0

	@ Reorganize to get the low accumulator words in low registers.
	mov	r8, r0
	mov	r0, r10
	mov	r10, r5
	mov	r2, r12
	mov	r3, r14
	mov	r14, r1
	mov	r1, r11
	mov	r11, r6
	mov	r12, r7

	@ Accumulator: r0:r1:r2:r3:r10:r11:r12:r14:r8

	@ Add a1*b1 to accumulator.
	add	r4, sp, #(4 * ((\id) + 8))
	ldm	r4!, { r5, r6, r7 }
	adds	r0, r5
	adcs	r1, r6
	adcs	r2, r7
	mov	r5, r10
	mov	r6, r11
	mov	r7, r12
	mov	r10, r0
	mov	r11, r1
	mov	r12, r2
	@ Accumulator: r10:r11:r12:r3:r5:r6:r7:r14:r8
	ldm	r4!, { r0, r1, r2 }
	adcs	r3, r0
	adcs	r5, r1
	adcs	r6, r2
	mov	r1, r14
	mov	r0, r8
	mov	r14, r3
	@ Accumulator: r10:r11:r12:r14:r5:r6:r7:r1:r0
	ldm	r4!, { r2, r3 }
	adcs	r7, r2
	adcs	r1, r3
	eors	r2, r2
	adcs	r0, r2

	@ Accumulator: r10:r11:r12:r14:r5:r6:r7:r1:r0

	@ The accumulator is necessarily nonnegative. We (virtually)
	@ extend it with three zeros, and add it to the destination
	@ buffer at the right place.

	@ Reorganization to get the low limbs into low registers.
	mov	r8, r0
	mov	r0, r10
	mov	r10, r5
	mov	r2, r12
	mov	r3, r14
	mov	r14, r1
	mov	r1, r11
	mov	r11, r6
	mov	r12, r7

	@ Accumulator: r0:r1:r2:r3:r10:r11:r12:r14:r8

	add	r4, sp, #(4 * ((\id) + 4))
	ldr	r5, [r4]
	ldr	r6, [r4, #4]
	adds	r0, r5
	adcs	r1, r6
	ldr	r5, [r4, #8]
	ldr	r6, [r4, #12]
	adcs	r2, r5
	adcs	r3, r6
	stm	r4!, { r0, r1, r2, r3 }
	mov	r0, r10
	mov	r1, r11
	mov	r2, r12
	mov	r3, r14
	ldr	r5, [r4]
	ldr	r6, [r4, #4]
	adcs	r0, r5
	adcs	r1, r6
	ldr	r5, [r4, #8]
	ldr	r6, [r4, #12]
	adcs	r2, r5
	adcs	r3, r6
	stm	r4!, { r0, r1, r2, r3 }
	ldm	r4!, { r0, r1, r2, r3 }
	mov	r7, r8
	adcs	r0, r7
	sbcs	r5, r5
	adds	r5, #1
	eors	r6, r6
	adds	r1, r5
	adcs	r2, r6
	adcs	r3, r6
	subs	r4, #16
	stm	r4!, { r0, r1, r2, r3 }
.endm

@ Input:
@   rx      input value to square (consumed)
@ Output:
@   rl:rh   square of the input
@ rx is consumed. rt is scratch.
@ Cost: 10
.macro SQR32  rl, rh, rx, rt
	uxth	\rl, \rx
	lsrs	\rx, \rx, #16
	movs	\rh, \rx
	muls	\rx, \rl
	muls	\rl, \rl
	muls	\rh, \rh
	lsls	\rt, \rx, #17
	lsrs	\rx, \rx, #15
	adds	\rl, \rt
	adcs	\rh, \rx
.endm

@ Input:
@   rx1:rx2    value to square (64-bit) (consumed)
@ Output:
@   rd0..rd3   square of the input value
@ rx1 and rx2 are consumed. rt1 and rt2 are scratch. All registers must
@ be different.
@ Cost: 14 + 3*cost(SQR32) = 44
.macro SQR64  rd0, rd1, rd2, rd3, rx1, rx2, rt1, rt2
	@ |a0-a1| -> rd0
	subs	\rd0, \rx1, \rx2
	sbcs	\rd1, \rd1
	eors	\rd0, \rd1
	subs	\rd0, \rd1

	@ a1^2 -> rd2:rd3
	SQR32  \rd2, \rd3, \rx2, \rd1

	@ (a0-a1)^2 -> rt1:rx2
	SQR32  \rt1, \rx2, \rd0, \rd1

	@ a0^2 -> rd0:rd1
	SQR32  \rd0, \rd1, \rx1, \rt2

	@ Subtract a0^2 and a1^2 from (a0-a1)^2 (in rt1:rx2:rt2)
	eors	\rx1, \rx1
	subs	\rt1, \rd0
	sbcs	\rx2, \rd1
	sbcs	\rt2, \rt2
	subs	\rt1, \rd2
	sbcs	\rx2, \rd3
	sbcs	\rt2, \rx1

	@ rt1:rx2:rt2 now contains -2*a0*a1, which we subtract from
	@ the result (with a one-word shift).
	subs	\rd1, \rt1
	sbcs	\rd2, \rx2
	sbcs	\rd3, \rt2
.endm

@ Input:
@   r1  pointer to value to square (128-bit)
@ Output:
@   Eight stack words, at index id
@ Cost: 89 + 3*cost(SQR64) = 221
@ Input:
@   id    output buffer index (in stack)
@   alt   if zero, then:
@            r1 = pointer to source integer (a)
@         else:
@            source integer (a) is in r0:r1:r2:r3
@         alt must be 0 or 1
@
@ Output: result goes to index itmp (8 words).
@ If alt is non-zero, then all output words are modified but only the
@ first two contain the correct result; the remaining 6 words are
@ returned in r2..r7.
@
@ All registers except r8 are consumed.
@ IMPORTANT: destination MUST NOT overlap with the source.
@ Cost:
@   alt = 0: 89 + 3*cost(SQR64) = 221
@   alt = 1: 76 + 3*cost(SQR64) = 208
.macro SQR128  id, alt
	.if (\alt) != 0
	@ Save source in high registers.
	mov	r10, r0
	mov	r11, r1
	mov	r12, r2
	mov	r14, r3
	@ a0^2 -> dst
	SQR64  r4, r5, r6, r7, r0, r1, r2, r3
	add	r0, sp, #(4 * (\id))
	stm	r0!, { r4, r5, r6, r7 }

	@ a1^2 -> dst + 4
	mov	r0, r12
	mov	r1, r14
	SQR64  r4, r5, r6, r7, r0, r1, r2, r3
	add	r0, sp, #(4 * ((\id) + 4))
	stm	r0!, { r4, r5, r6, r7 }

	@ Load source into r4..r7
	mov	r4, r10
	mov	r5, r11
	mov	r6, r12
	mov	r7, r14
	.else
	@ a0^2 -> dst
	ldm	r1!, { r2, r3 }
	mov	r10, r1
	SQR64  r4, r5, r6, r7, r2, r3, r0, r1
	add	r0, sp, #(4 * (\id))
	stm	r0!, { r4, r5, r6, r7 }

	@ a1^2 -> dst + 4
	mov	r1, r10
	ldm	r1, { r0, r1 }
	SQR64  r4, r5, r6, r7, r0, r1, r2, r3
	add	r0, sp, #(4 * ((\id) + 4))
	stm	r0!, { r4, r5, r6, r7 }

	@ Load source into r0..r3
	mov	r0, r10
	subs	r0, #8
	ldm	r0!, { r4, r5, r6, r7 }
	.endif

	@ Compute (a0-a1)^2 -> r0..r3
	subs	r4, r6
	sbcs	r5, r7
	sbcs	r3, r3
	eors	r4, r3
	eors	r5, r3
	subs	r4, r3
	sbcs	r5, r3
	SQR64  r0, r1, r2, r3, r4, r5, r6, r7

	@ Subtract a0^2 from accumulator.
	add	r4, sp, #(4 * (\id))
	ldm	r4, { r4, r5, r6, r7 }
	subs	r0, r4
	sbcs	r1, r5
	sbcs	r2, r6
	sbcs	r3, r7
	sbcs	r6, r6

	@ Accumulator is: r0:r1:r2:r3:r6

	@ Subtract a1^2 from accumulator.
	add	r4, sp, #(4 * ((\id) + 4))
	ldm	r4!, { r5, r7 }
	subs	r0, r5
	sbcs	r1, r7
	ldm	r4!, { r5, r7 }
	sbcs	r2, r5
	sbcs	r3, r7
	eors	r4, r4
	sbcs	r6, r4

	@ Sign-extend accumulator and move some words to high registers.
	mov	r10, r1
	mov	r11, r2
	mov	r12, r3
	mov	r14, r6
	asrs	r1, r6, #31

	@ Accumulator is: r8:r10:r11:r12:r14:r7

	@ Subtract accumulator from result.
	add	r2, sp, #(4 * ((\id) + 2))
	ldm	r2, { r2, r3, r4, r5, r6, r7 }
	subs	r2, r0
	mov	r0, r10
	sbcs	r3, r0
	mov	r0, r11
	sbcs	r4, r0
	mov	r0, r12
	sbcs	r5, r0
	mov	r0, r14
	sbcs	r6, r0
	sbcs	r7, r1

	.if (\alt) == 0
	@ Write out the result (low two words are already there).
	add	r1, sp, #(4 * ((\id) + 2))
	stm	r1!, { r2, r3, r4, r5, r6, r7 }
	.endif
.endm

@ Input:
@   r1     pointer to value to square (128-bit)
@   id     index for output buffer (in stack)
@   itmp   index to free stack slots (8 slots)
@ Output:
@   16 stack words, at index id
@ Cost: 165 + 2*cost(SQR128[alt=0]) + 1*cost(SQR128[alt=1]) = 815
.macro SQR256  id, itmp
	@ Save destination pointer in r8 (which is not modified by SQR128).
	mov	r8, r1

	@ a0^2 -> dst
	SQR128  ((\id) + 0), 0

	@ a1^2 -> dst + 8
	mov	r1, r8
	adds	r1, #16
	SQR128  ((\id) + 8), 0

	@ (a0-a1)^2 -> tmp
	mov	r1, r8
	ldm	r1, { r0, r1, r2, r3, r4, r5, r6, r7 }
	subs	r0, r4
	sbcs	r1, r5
	sbcs	r2, r6
	sbcs	r3, r7
	sbcs	r4, r4
	eors	r0, r4
	eors	r1, r4
	eors	r2, r4
	eors	r3, r4
	subs	r0, r4
	sbcs	r1, r4
	sbcs	r2, r4
	sbcs	r3, r4
	SQR128  (\itmp), 1

	@ Load accumulator into r0:r1:r2:r3:r10:r11:r12:r14
	add	r0, sp, #(4 * (\itmp))
	ldm	r0, { r0, r1 }
	mov	r10, r4
	mov	r11, r5
	mov	r12, r6
	mov	r14, r7

	@ Subtract a0^2 from accumulator.
	add	r4, sp, #(4 * (\id))
	ldm	r4, { r4, r5, r6, r7 }
	subs	r0, r4
	sbcs	r1, r5
	sbcs	r2, r6
	sbcs	r3, r7
	mov	r4, r10
	mov	r5, r11
	mov	r6, r12
	mov	r7, r14
	mov	r10, r0
	mov	r11, r1
	mov	r12, r2
	mov	r14, r3
	add	r0, sp, #(4 * ((\id) + 4))
	ldm	r0, { r0, r1, r2, r3 }
	sbcs	r4, r0
	sbcs	r5, r1
	sbcs	r6, r2
	sbcs	r7, r3
	sbcs	r0, r0
	mov	r8, r0

	@ Accumulator is: r10:r11:r12:r14:r4:r5:r6:r7:r8

	@ Subtract a1^2 from accumulator.
	mov	r0, r10
	mov	r1, r11
	mov	r2, r12
	mov	r3, r14
	mov	r10, r4
	mov	r11, r5
	mov	r12, r6
	mov	r14, r7
	add	r4, sp, #(4 * ((\id) + 8))
	ldm	r4, { r4, r5, r6, r7 }
	subs	r0, r4
	sbcs	r1, r5
	sbcs	r2, r6
	sbcs	r3, r7
	mov	r4, r10
	mov	r5, r11
	mov	r6, r12
	mov	r7, r14
	mov	r10, r0
	mov	r11, r1
	mov	r12, r2
	mov	r14, r3
	add	r0, sp, #(4 * ((\id) + 12))
	ldm	r0, { r0, r1, r2, r3 }
	sbcs	r4, r0
	sbcs	r5, r1
	sbcs	r6, r2
	sbcs	r7, r3
	sbcs	r0, r0
	add	r8, r0

	@ Accumulator is: r10:r11:r12:r14:r4:r5:r6:r7:r8

	@ Subtract accumulator from result.
	mov	r0, r10
	mov	r1, r11
	mov	r2, r12
	mov	r3, r14
	mov	r10, r4
	mov	r11, r5
	mov	r12, r6
	mov	r14, r7
	add	r4, sp, #(4 * ((\id) + 4))
	ldm	r4, { r4, r5, r6, r7 }
	subs	r4, r0
	sbcs	r5, r1
	sbcs	r6, r2
	sbcs	r7, r3
	add	r0, sp, #(4 * ((\id) + 4))
	stm	r0!, { r4, r5, r6, r7 }
	mov	r0, r10
	mov	r1, r11
	mov	r2, r12
	mov	r3, r14
	add	r4, sp, #(4 * ((\id) + 8))
	ldm	r4, { r4, r5, r6, r7 }
	sbcs	r4, r0
	sbcs	r5, r1
	sbcs	r6, r2
	sbcs	r7, r3
	add	r0, sp, #(4 * ((\id) + 8))
	stm	r0!, { r4, r5, r6, r7 }
	@ For the remaining four words, we need to sign-extend the upper
	@ word of the accumulator. This can be done with rev and sxtb
	@ because that upper word is small in absolute value, so its
	@ top byte is either 0x00 or 0xFF; rev and sxtb do not modify the 
	@ flags.
	mov	r0, r8
	rev	r1, r0
	sxtb	r1, r1
	add	r4, sp, #(4 * ((\id) + 12))
	ldm	r4, { r4, r5, r6, r7 }
	sbcs	r4, r0
	sbcs	r5, r1
	sbcs	r6, r1
	sbcs	r7, r1
	add	r0, sp, #(4 * ((\id) + 12))
	stm	r0!, { r4, r5, r6, r7 }
.endm

@ Input:
@   rx = input word to multiply (32-bit)
@   rm = 2*MQ (unmodified)
@ This function computes rx*rm, and outputs the result as rx + ru:rc
@ i.e. rx, ru and rc are modified so that the result can be obtained by
@ performing the addition above. Moreover, the output value of rc fits on
@ 16 bits.
@ Cost: 6
.macro MUL32x16_START  rx, rc, ru, rm
	lsrs	\ru, \rx, #16
	uxth	\rx, \rx
	muls	\ru, \rm
	muls	\rx, \rm
	lsrs	\rc, \ru, #16
	lsls	\ru, \ru, #16
.endm

@ Input:
@   rx = input word to multiply (32-bit)
@   rc = carry (16-bit)
@   rm = 2*MQ (unmodified)
@ This function computes rx*rm + rc, and outputs the result as rx + ru:rc
@ i.e. rx, ru and rc are modified so that the result can be obtained by
@ performing the addition above. Moreover, the output value of rc fits on
@ 16 bits.
@ rt is a scratch register.
@ Cost: 7
.macro MUL32x16_CC_1  rx, rc, ru, rm, rt
	lsrs	\rt, \rx, #16
	uxth	\rx, \rx
	muls	\rt, \rm
	muls	\rx, \rm
	lsls	\ru, \rt, #16
	orrs	\ru, \rc
	lsrs	\rc, \rt, #16
.endm

@ Input:
@   rx = input word to multiply (32-bit)
@   rc = carry (16-bit)
@   rm = 2*MQ (unmodified)
@ This function computes rx*rm + rc, and outputs the result as rx + ru:rc
@ i.e. rx, ru and rc are modified so that the result can be obtained by
@ performing the addition above. Moreover, the output value of rc fits on
@ 16 bits.
@ This macro differs from MUL32x16_CC_1 in the following:
@  - on output, ru anc rc are exchanged (carry is in ru, not rc)
@  - no extra scratch register is needed
@  - cost is higher by 2 cycles
@ Cost: 9
.macro MUL32x16_CC_2  rx, rc, ru, rm
	lsrs	\ru, \rx, #16
	uxth	\rx, \rx
	muls	\ru, \rm
	muls	\rx, \rm
	rev	\ru, \ru
	rev16	\ru, \ru
	eors	\rc, \ru
	uxth	\ru, \ru
	eors	\rc, \ru
.endm

@ iin   input buffer index (on stack)
@ mq2   label for constant 2*MQ
@ Output is in r0:r1:r2:r3:r4:r10:r11:r7
@ Cost: 159
.macro REDUCE  iin, mq2
	@ Set r5 to zero (for carry propagation).
	eors	r5, r5

	@ Set r7 to multiplier (2*MQ).
	@ (We could gain one cycle by using
	ldr	r7, \mq2

	@ Multiply high word by 2*MQ.
	add	r0, sp, #(4 * ((\iin) + 8))
	ldm	r0, { r0, r1 }
	MUL32x16_START  r0, r6, r2, r7
	MUL32x16_CC_1   r1, r6, r3, r7, r4
	adds	r0, r2
	adcs	r1, r3
	adcs	r6, r5
	mov	r10, r0
	mov	r11, r1
	add	r0, sp, #(4 * ((\iin) + 10))
	ldm	r0, { r0, r1 }
	MUL32x16_CC_1   r0, r6, r2, r7, r4
	MUL32x16_CC_1   r1, r6, r3, r7, r4
	adds	r0, r2
	adcs	r1, r3
	adcs	r6, r5
	mov	r12, r0
	mov	r14, r1
	add	r0, sp, #(4 * ((\iin) + 12))
	ldm	r0, { r0, r1 }
	MUL32x16_CC_1   r0, r6, r2, r7, r4
	MUL32x16_CC_1   r1, r6, r3, r7, r4
	adds	r0, r2
	adcs	r1, r3
	adcs	r6, r5
	add	r2, sp, #(4 * ((\iin) + 14))
	ldm	r2, { r2, r3 }
	MUL32x16_CC_1   r2, r6, r4, r7, r5
	MUL32x16_CC_2   r3, r6, r5, r7
	adds	r2, r4
	adcs	r3, r6
	sbcs	r4, r4
	adds	r4, #1
	adds	r5, r4

	@ We have the value in r10:r11:r12:r14:r0:r1:r2:r3:r5.
	@ We now add the low word into this value.

	@ Reorganization.
	mov	r8, r5
	mov	r4, r10
	mov	r5, r11
	mov	r6, r12
	mov	r7, r14
	mov	r10, r0
	mov	r11, r1
	mov	r12, r2
	mov	r14, r3

	@ Accumulator is now r4:r5:r6:r7:r10:r11:r12:r14:r8
	add	r0, sp, #(4 * (\iin))
	ldm	r0, { r0, r1, r2, r3 }
	adds	r4, r0
	adcs	r5, r1
	adcs	r6, r2
	adcs	r7, r3
	mov	r0, r10
	mov	r1, r11
	mov	r2, r12
	mov	r3, r14
	mov	r10, r4
	mov	r11, r5
	mov	r12, r6
	mov	r14, r7
	@ Accumulator is now r10:r11:r12:r14:r0:r1:r2:r3:r8
	add	r4, sp, #(4 * ((\iin) + 4))
	ldm	r4, { r4, r5, r6, r7 }
	adcs	r4, r0
	adcs	r5, r1
	adcs	r6, r2
	adcs	r7, r3
	sbcs	r0, r0
	add	r0, r8
	adds	r0, #1
	@ Accumulator is now r10:r11:r12:r14:r4:r5:r6:r7:r0

	@ Second reduction round: multiply high word by 2*MQ, and add it.
	mov	r8, r7
	@ Accumulator is now r10:r11:r12:r14:r4:r5:r6:r8:r0
	eors	r7, r7
	ldr	r3, \mq2
	muls	r3, r0
	mov	r0, r10
	mov	r1, r11
	mov	r2, r12
	adds	r0, r3
	adcs	r1, r7
	adcs	r2, r7
	mov	r3, r14
	adcs	r3, r7
	adcs	r4, r7
	adcs	r5, r7
	adcs	r6, r7
	mov	r10, r5
	mov	r11, r6
	mov	r5, r8
	adcs	r7, r5
	@ Accumulator is now r0:r1:r2:r3:r4:r10:r11:r7
	@ There may still be a carry, which we then must wrap again; but
	@ it won't trigger any extra carry.
	sbcs	r5, r5
	ldr	r6, \mq2
	bics	r6, r5
	adds	r0, r6
.endm

@ =======================================================================
@ void gf_mul_inner(gf *d, const gf *a, const gf *b)
@
@ All registers are modified.
@
@ Cost: 23 + cost(MUL256x256) + cost(REDUCE) = 1464
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	gf_mul_inner, %function
gf_mul_inner:
	sub	sp, #112
	push	{ r0, lr }

	MUL256x256  2, 18

	@ Reduction modulo p.
	REDUCE  2, const_gf_mul_2mq

	@ Value is in r0:r1:r2:r3:r4:r10:r11:r7.
	pop	{ r5 }
	stm	r5!, { r0, r1, r2, r3, r4 }
	mov	r0, r10
	mov	r1, r11
	stm	r5!, { r0, r1, r7 }

	pop	{ r0 }
	add	sp, #112
	bx	r0
	.align	2
const_gf_mul_2mq:
	.long	2 * MQ
	.size	gf_mul_inner, .-gf_mul_inner

@ =======================================================================
@ void gf_sqr_inner(gf *d, const gf *a)
@
@ All registers are modified.
@
@ Cost: 23 + cost(SQR256) + cost(REDUCE) = 997
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	gf_sqr_inner, %function
gf_sqr_inner:
	sub	sp, #96
	push	{ r0, lr }

	SQR256  2, 18

	@ Reduction modulo p.
	REDUCE  2, const_gf_sqr_2mq

	@ Value is in r0:r1:r2:r3:r4:r10:r11:r7.
	pop	{ r5 }
	stm	r5!, { r0, r1, r2, r3, r4 }
	mov	r0, r10
	mov	r1, r11
	stm	r5!, { r0, r1, r7 }

	pop	{ r0 }
	add	sp, #96
	bx	r0
	.align	2
const_gf_sqr_2mq:
	.long	2 * MQ
	.size	gf_sqr_inner, .-gf_sqr_inner

@ =======================================================================
@ void gf_normalize_inner(gf *d, const gf *a)
@
@ Registers r0 to r8 and r10 to r12 are modified.
@
@ Cost: 78
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	gf_normalize_inner, %function
gf_normalize_inner:
	@ Load complete source.
	mov	r8, r0
	ldm	r1, { r0, r1, r2, r3, r4, r5, r6, r7 }

	@ If the top bit is set, propagate it (with 2^255 = MQ mod p)
	mov	r10, r6
	lsls	r6, r7, #1
	lsrs	r6, r6, #1
	mov	r11, r6
	asrs	r7, r7, #31
	ldr	r6, const_gf_normalize_mq
	ands	r6, r7
	eors	r7, r7
	adds	r0, r6
	adcs	r1, r7
	adcs	r2, r7
	adcs	r3, r7
	adcs	r4, r7
	adcs	r5, r7
	mov	r6, r10
	adcs	r6, r7
	adcs	r7, r7
	add	r11, r7

	@ Value is now at most 2^255+MQ-1, in r0:r1:r2:r3:r4:r5:r6:r11
	@ (thus, max value for r11 is 0x80000000).
	@ We subtract p.
	ldr	r7, const_gf_normalize_minusmq
	subs	r0, r7
	movs	r7, #0xFF
	sxtb	r7, r7
	sbcs	r1, r7
	sbcs	r2, r7
	sbcs	r3, r7
	sbcs	r4, r7
	sbcs	r5, r7
	sbcs	r6, r7
	sbcs	r7, r7
	add	r7, r11
	@ Top word: we only processed the borrow but have not subtracted
	@ 0x7FFFFFFF from it; we do it now.
	mov	r10, r6
	ldr	r6, const_gf_normalize_two31m1
	subs	r7, r6

	@ The value of r11 was at most 0x80000000; thus, after processing
	@ the borrow, value of r7 was in -1..+2^31-1. Subtracting the
	@ top word of p (2^31-1) leads to a value in -2^31..0, i.e.
	@ entirely within the range of signed values. The top bit of r7
	@ thus gives us the sign of the value. If it is negative, then
	@ we must add back p.

	@ Value is currently in: r0:r1:r2:r3:r4:r5:r10:r7
	mov	r11, r7
	asrs	r7, r7, #31
	lsrs	r6, r7, #1
	adds	r6, #1       @ 0x80000000 if adding, 1 if not adding.
	add	r11, r6
	ldr	r6, const_gf_normalize_minusmq
	ands	r6, r7

	@ Now: if adding, then r6 = -MQ and r7 = -1; if not adding, then
	@ r6 = r7 = 0.
	adds	r0, r6
	adcs	r1, r7
	adcs	r2, r7
	adcs	r3, r7
	adcs	r4, r7
	adcs	r5, r7
	mov	r6, r10
	adcs	r6, r7
	sbcs	r7, r7
	add	r11, r7
	@ Top word handling:
	@ 'sbcs r7, r7' sets r7 to c-1, where c is the carry at that point.
	@ We want to add c+0x7FFFFFFF to the initial top word if adding,
	@ or 0 if not. However, if not adding, then c is zero at this point
	@ (we just added zero to all other words); thus, we can
	@ systematically add (c-1)+1 to the top word, and 0x7FFFFFFF
	@ conditionally. This is equivalent to adding c-1 here (and we got
	@ c-1 in r7), and either 1 or 0x80000000 conditionally, which we
	@ did earlier.

	@ Result is now normalized, in r0:r1:r2:r3:r4:r5:r6:r11.
	mov	r7, r8
	stm	r7!, { r0, r1, r2, r3, r4, r5, r6 }
	mov	r0, r11
	str	r0, [r7]
	bx	lr
	.align	2
const_gf_normalize_mq:
	.long	MQ
const_gf_normalize_minusmq:
	.long	-MQ
const_gf_normalize_two31m1:
	.long	0x7FFFFFFF
	.size	gf_normalize_inner, .-gf_normalize_inner

@ =======================================================================
@ uint32_t gf_iszero(const gf *a)
@
@ Registers r0 to r7 are modified (r0 is the returned value).
@
@ Cost: 55
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	gf_iszero_inner, %function
gf_iszero_inner:
	@ There are three possible representations for 0: 0, p and 2*p.
	@ We accumulate three masks in r1, r6 and r7; at the end:
	@  - if value is 0, then r1 = 0; otherwise, r1 != 0
	@  - if value is p, then r6 = -1; otherwise, r6 != -1
	@  - if value is 2*p, then r7 = -1; otherwise, r7 != -1
	ldm	r0!, { r1, r2, r3, r4 }
	ldr	r5, const_gf_iszero_mq
	adds	r6, r1, r5
	adds	r5, r5
	adds	r7, r1, r5
	mvns	r6, r6
	mvns	r7, r7
	orrs	r1, r2
	orrs	r1, r3
	orrs	r1, r4
	ands	r6, r2
	ands	r6, r3
	ands	r6, r4
	ands	r7, r2
	ands	r7, r3
	ands	r7, r4
	ldm	r0, { r0, r2, r3, r4 }
	orrs	r1, r0
	orrs	r1, r2
	orrs	r1, r3
	orrs	r1, r4
	ands	r6, r0
	ands	r6, r2
	ands	r6, r3
	ands	r7, r0
	ands	r7, r2
	ands	r7, r3
	ands	r7, r4
	movs	r0, #1
	lsls	r0, r0, #31
	eors	r4, r0
	ands	r6, r4

	@ We do a bitwise NOT of r6 and r7.
	mvns	r6, r6
	mvns	r7, r7
	@ Now, result should be 1 if and only if one of r1, r6 and r7 is
	@ equal to zero.
	rsbs	r0, r1, #0
	orrs	r1, r0
	rsbs	r0, r6, #0
	orrs	r6, r0
	rsbs	r0, r7, #0
	orrs	r7, r0
	ands	r1, r6
	ands	r1, r7
	asrs	r1, r1, #31
	adds	r0, r1, #1
	bx	lr
	.align	2
const_gf_iszero_mq:
	.long	MQ
const_gf_iszero_minusmq:
	.long	-MQ
const_gf_iszero_two31m1:
	.long	0x7FFFFFFF
	.size	gf_iszero_inner, .-gf_iszero_inner

@ =======================================================================

@ If the top n bits of r2 are zero, then shift r2 to the left by n bits,
@ and add n to r6; otherwise, do not modify r2 or r6.
@ r3, r4 and r5 are scratch; other registers are unmodified.
@ Value n is provided as its base-2 logarithm logn, which must be in the
@ 0..4 range.
@ Cost: 9 if logn != 0, 7 if logn == 0
.macro CLZ_HI  logn
	.if (\logn) != 0
	lsrs	r3, r2, #(32 - (1 << \logn))
	subs	r3, #1
	asrs	r3, r3, #31
	lsls	r4, r2, #(1 << \logn)
	subs	r4, r2
	ands	r4, r3
	adds	r2, r4
	lsls	r3, r3, #(\logn)
	subs	r6, r3
	.else
	asrs	r3, r2, #31
	mvns	r3, r3
	lsls	r4, r2, #1
	subs	r4, r2
	ands	r4, r3
	adds	r2, r4
	subs	r6, r3
	.endif
.endm

@ =======================================================================

@ Start of gathering of the top two words of a and b.
@   r0   pointer to a (updated)
@   r1   pointer to b (updated)
@   r7   remaining index (updated)
@ Top words are accumulated in r10 and r12; second-to-top words go to r8
@ and r11. This macro initializes r8..r12 and processes the lowest two
@ words of a and b.
@ Output:
@   r6   -1 on index match, 0 otherwise
@ Cost: 23
.macro APPROX_START
	ldm	r0!, { r2, r4 }
	ldm	r1!, { r3, r5 }
	@ Second word is r2/r3 iff r7 <= 1
	subs	r7, #2
	sbcs	r6, r6
	ands	r2, r6
	ands	r3, r6
	mov	r8, r2
	mov	r11, r3
	movs	r2, r4
	ands	r2, r6
	mov	r10, r2
	ands	r6, r5
	mov	r12, r6
	subs	r7, #1
	sbcs	r6, r6
	ands	r4, r6
	ands	r5, r6
	add	r8, r4
	add	r11, r5
.endm

@ Continuation of gathering of the top two words of a and b.
@   r0   pointer to a (updated)
@   r1   pointer to b (updated)
@   r7   remaining index (updated)
@ One word of a and b is processed by this macro.
@ Output:
@   r6   -1 on index match, 0 otherwise
@ Cost: 15
.macro APPROX_CONT
	ldm	r0!, { r2 }
	ldm	r1!, { r3 }
	movs	r4, r2
	ands	r4, r6
	add	r10, r4
	ands	r6, r3
	add	r12, r6
	subs	r7, #1
	sbcs	r6, r6
	ands	r2, r6
	ands	r3, r6
	add	r8, r2
	add	r11, r3
.endm

@ =======================================================================
@ (uint32_t,uint32_t) approximate_ab(const gf *a, const gf *b)
@
@ This routine approximates the two provided 255-bit nonnegative integers
@ into two 32-bit values (b must be non-zero). If n = max(len(a), len(b)),
@ then, on output:
@
@  - If n <= 32, then r0 contains a[0] and r1 contains b[0]
@  - Otherwise:
@      r0 contains (a mod 2^15) + 2^15 * floor(a / 2^15)
@      r1 contains (b mod 2^15) + 2^15 * floor(b / 2^15)
@
@ All registers are consumed.
@
@ Cost: 251
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	approximate_ab, %function
approximate_ab:
	@ Find index of top non-zero word in c = a OR b.
	@ We already know that b[0] != 0.
	ldm	r0!, { r2, r3, r4 }
	ldm	r1!, { r2, r5, r6 }
	orrs	r3, r5
	orrs	r4, r6
	cmp	r3, #1
	sbcs	r7, r7
	cmp	r4, #1
	sbcs	r6, r6
	ands	r7, r6
	adds	r7, r6
	ldm	r0!, { r3, r4 }
	ldm	r1!, { r5, r6 }
	orrs	r3, r5
	orrs	r4, r6
	cmp	r3, #1
	sbcs	r6, r6
	ands	r7, r6
	adds	r7, r6
	cmp	r4, #1
	sbcs	r6, r6
	ands	r7, r6
	adds	r7, r6
	ldm	r0!, { r3, r4 }
	ldm	r1!, { r5, r6 }
	orrs	r3, r5
	orrs	r4, r6
	cmp	r3, #1
	sbcs	r6, r6
	ands	r7, r6
	adds	r7, r6
	cmp	r4, #1
	sbcs	r6, r6
	ands	r7, r6
	adds	r7, r6
	ldm	r0!, { r3 }
	ldm	r1!, { r5 }
	orrs	r3, r5
	cmp	r3, #1
	sbcs	r6, r6
	ands	r7, r6
	adds	r7, r6

	@ If j is the max index such that c[j] != 0, then r7 contains j-7.
	@ We translate it into the 0..7 range.
	adds	r7, #7

	@ We now scan the words again to extract the corresponding top
	@ two words of a and b (in r8:r10 and r11:r12, respectively).
	subs	r0, #32
	subs	r1, #32
	APPROX_START
	APPROX_CONT
	APPROX_CONT
	APPROX_CONT
	APPROX_CONT
	APPROX_CONT
	APPROX_CONT

	@ We have the two top words; we OR them to count the required
	@ shift.
	mov	r2, r10
	mov	r3, r12
	orrs	r2, r3
	eors	r6, r6
	CLZ_HI  4
	CLZ_HI  3
	CLZ_HI  2
	CLZ_HI  1
	CLZ_HI  0

	@ Shift count is in r6 (0 to 31 range). If r2 != 0, then c did
	@ not fit in a single word, and the count is correct; we
	@ assemble the output words. If r2 == 0, then c was 32 bits or
	@ shorter; this will be fixed afterwards.
	cmp	r2, #1
	sbcs	r7, r7   @ -1 if r2 == 0, 0 otherwise
	subs	r0, #32
	subs	r1, #32
	ldr	r2, [r0]
	ldr	r3, [r1]
	mov	r0, r10
	mov	r1, r12
	lsls	r0, r6
	lsls	r1, r6
	movs	r5, #32
	subs	r6, r5, r6
	mov	r4, r8
	mov	r5, r11
	lsrs	r4, r6
	lsrs	r5, r6
	orrs	r0, r4
	orrs	r1, r5

	@ We have the shifted top 32 bits of a and b in r0 and r1,
	@ respectively. We also have a[0] and b[0] in r2 and r3. We
	@ replace the low 15 bits of r0 and r1 with those of r2 and
	@ r3. However, if r7 == -1, then the two inputs could fit in
	@ 32 bits each, and in that case we want to replace r0 and r1
	@ with r2 and r3, respectively.
	ldr	r6, const_approximate_ab_mask15
	orrs	r6, r7
	bics	r0, r6
	bics	r1, r6
	ands	r2, r6
	ands	r3, r6
	orrs	r0, r2
	orrs	r1, r3

	bx	lr
	.align	2
const_approximate_ab_mask15:
	.long	0x7FFF
	.size	approximate_ab, .-approximate_ab

@ =======================================================================

@ One step of lin2.
@ Input:
@   r0    pointer to a and b (unmodified)
@   r4    current carry word for a (ignored if i == 0)
@   r5    current carry word for b (ignored if i == 0)
@   r8    update factor f0 (unmodified)
@   r10   update factor g0 (unmodified)
@   r11   update factor f1 (unmodified)
@   r12   update factor g1 (unmodified)
@ Output:
@   r4    new carry word for a
@   r5    new carry word for b
@ r0, r1 and r8+ are unmodified. r4 and r5 are updated. r2, r3, r6 and r7
@ are scratch.
@ Cost: 20 if i == 0, 22 otherwise
.macro LIN2_STEP  i
	@ Load the next 16-bit limbs xa and xb.
	ldrh	r2, [r0, #(2 * (\i))]
	ldrh	r3, [r0, #(2 * (\i) + 32)]

	@ Multiply with f0 and g0, and add (into r6). We know that
	@ |f0|+|g0| <= 2^15; therefore:
	@    |xa*f0 + xb*g0| <= 2^15*(2^16-1) < 2^31
	@ Thus, this value fits in 32 bits (with its sign bit).
	mov	r6, r8
	mov	r7, r10
	muls	r6, r2
	muls	r7, r3
	adds	r6, r7

	@ Similarly compute |xa*f1 + xb*g1| into r7 (xa is consumed).
	mov	r7, r11
	muls	r2, r7
	mov	r7, r12
	muls	r7, r3
	adds	r7, r2

	@ Add carries from lower limbs, and store.
	@ Suppose that the input carry for a is c such that:
	@    -2^15 + 1 <= c <= 2^15 - 1
	@ Then, we know that:
	@    -2^15 <= f0 + g0 <= 2^15
	@    0 <= xa, xb <= 2^16-1
	@ Therefore:
	@    -(2^16-1)*2^15 <= xa*f0 + xb*g0 <= (2^16-1)*2^15
	@ and:
	@    -2^31 + 2^15 <= xa*f0 + xb*g0 + c <= 2^31 - 1
	@ Thus, addition of the carry does not overflow. Moreover,
	@ after right-shifting by 16 bits, the new carry will also
	@ fulfill the inequation above (i.e. |c| <= 2^15-1).
	.if (\i) != 0
	adds	r6, r4
	adds	r7, r5
	.endif
	asrs	r4, r6, #16
	asrs	r5, r7, #16

	@ Write the two produced 16-bit limbs.
	strh	r6, [r0, #(2 * (\i))]
	strh	r7, [r0, #(2 * (\i) + 32)]
.endm

@ =======================================================================
@ (int32_t, int32_t) lin2(gf *ab, int32_t f0, int32_t g0,
@                         int32_t f1, int32_t g1)
@
@ This routine replaces a and b with, respectively, a*f0+b*g0 and
@ a*f1+b*g1. f0, g0, f1 and g1 are signed integers, in the -2^15+1..+2^15
@ range. The results are also signed; the two extra words are returned. Values
@ a and b are assumed to be consecutive in RAM.
@
@ ABI:
@   r0 = pointer to ab
@   r8 = f0
@   r10 = g0
@   r11 = f1
@   r12 = g1
@ High registers (r8+) are preserved, except r14 (used for the call).
@
@ The extra words for a and b are returned in r4 and r5, respectively.
@ r0 is unmodified. r1 is unmodified.
@
@ Cost: 352
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	lin2, %function
lin2:
	LIN2_STEP  0
	LIN2_STEP  1
	LIN2_STEP  2
	LIN2_STEP  3
	LIN2_STEP  4
	LIN2_STEP  5
	LIN2_STEP  6
	LIN2_STEP  7
	LIN2_STEP  8
	LIN2_STEP  9
	LIN2_STEP  10
	LIN2_STEP  11
	LIN2_STEP  12
	LIN2_STEP  13
	LIN2_STEP  14
	LIN2_STEP  15

	bx	lr
	.size	lin2, .-lin2

@ =======================================================================
@ uint32_t s256_lin_div15_abs(gf *ab, int32_t f0, int32_t g0,
@                             int32_t f1, int32_t g1)
@
@ This routine replaces a and b with, respectively, (a*f0+b*g0)/2^15 and
@ (a*f1+b*g1)/2^15. The divisions are assumed to be exact (i.e. the 15
@ low bits are dropped. The two 255-bit integers a and b are consecutive
@ in RAM, starting at address 'ab'.
@
@ f0, g0, f1 and g1 are signed integers, in the -2^15+1..+2^15 range.
@ Moreover, it is known that f0+g0 and f1+g1 are also both in the
@ -2^15+1..+2^15 range.
@
@ If the resulting value of a is negative, then it is negated. Similarly,
@ if b is negative, then it is negated. Returned value (in r0) is:
@   0  both values were positive, no negation applied
@   1  a was negative, negation applied on a
@   2  b was negative, negation applied on b
@   3  both values were negative, two negations applied
@
@ ABI:
@   r0 = pointer to ab
@   r2 = f0
@   r3 = g0
@   r4 = f1
@   r5 = g1
@ No register is preserved.
@
@ Cost: 163 + cost(lin2) = 515
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	s256_lin_div15_abs, %function
s256_lin_div15_abs:
	push	{ lr }

	@ Move update factors to high registers (this leaves more room
	@ for computations, and does not imply much overhead since we
	@ have to copy them for each multiplication).
	mov	r8, r2
	mov	r10, r3
	mov	r11, r4
	mov	r12, r5

	@ Perform linear combination.
	bl	lin2

	@ We have done all 16 iterations. r0 now points to b. We
	@ must apply the right-shift by 15 bits (division by 2^15)
	@ and also negate values if they are negative.

	@ Save extra word of b in r8.
	mov	r8, r5

	@ Process a. Extra word is in r4.

	@ We keep the current write pointer in r11.
	mov	r11, r0
	@ Read a0, a1, a2, a3.
	ldm	r0, { r0, r1, r2, r3 }
	@ Conditionally negate a, and keep the borrow in r5.
	asrs	r6, r4, #31
	eors	r0, r6
	eors	r1, r6
	eors	r2, r6
	eors	r3, r6
	subs	r0, r6
	sbcs	r1, r6
	sbcs	r2, r6
	sbcs	r3, r6
	sbcs	r5, r5
	@ Right-shift words by 15 bits. High bits of b3 are delayed (in r7).
	lsrs	r0, r0, #15
	lsls	r7, r1, #17
	orrs	r0, r7
	lsrs	r1, r1, #15
	lsls	r7, r2, #17
	orrs	r1, r7
	lsrs	r2, r2, #15
	lsls	r7, r3, #17
	orrs	r2, r7
	lsrs	r7, r3, #15
	@ Write the new a0, a1 and a2.
	mov	r6, r11
	stm	r6!, { r0, r1, r2 }
	mov	r11, r6
	@ Read a4, a5, a6 and a7. Word a8 is the extra word (already in r4).
	adds	r6, #4
	ldm	r6!, { r0, r1, r2, r3 }
	@ Finish conditional negation; also, save the original sign in r10.
	asrs	r6, r4, #31
	mov	r10, r6
	eors	r0, r6
	eors	r1, r6
	eors	r2, r6
	eors	r3, r6
	eors	r4, r6
	rsbs	r5, r5, #0
	sbcs	r0, r6
	sbcs	r1, r6
	sbcs	r2, r6
	sbcs	r3, r6
	sbcs	r4, r6
	@ Right-shift words by 15 bits. We have some delayed bits in r7.
	lsrs	r6, r0, #15
	lsls	r0, r0, #17
	orrs	r0, r7
	lsrs	r7, r1, #15
	lsls	r1, r1, #17
	orrs	r1, r6
	lsrs	r6, r2, #15
	lsls	r2, r2, #17
	orrs	r2, r7
	lsrs	r7, r3, #15
	lsls	r3, r3, #17
	orrs	r3, r6
	lsls	r5, r4, #17
	orrs	r5, r7
	@ Write the new b3, b4, b5, b6 and b7.
	mov	r7, r11
	stm	r7!, { r0, r1, r2, r3, r5 }

	@ Process b. Same code as above, extra word in r4.

	mov	r4, r8
	mov	r11, r7
	@ Read a0, a1, a2, a3.
	ldm	r7!, { r0, r1, r2, r3 }
	@ Conditionally negate b, and keep the borrow in r5.
	asrs	r6, r4, #31
	eors	r0, r6
	eors	r1, r6
	eors	r2, r6
	eors	r3, r6
	subs	r0, r6
	sbcs	r1, r6
	sbcs	r2, r6
	sbcs	r3, r6
	sbcs	r5, r5
	@ Right-shift words by 15 bits. High bits of a3 are delayed (in r7).
	lsrs	r0, r0, #15
	lsls	r7, r1, #17
	orrs	r0, r7
	lsrs	r1, r1, #15
	lsls	r7, r2, #17
	orrs	r1, r7
	lsrs	r2, r2, #15
	lsls	r7, r3, #17
	orrs	r2, r7
	lsrs	r7, r3, #15
	@ Write the new a0, a1 and a2.
	mov	r6, r11
	stm	r6!, { r0, r1, r2 }
	mov	r11, r6
	@ Read a4, a5, a6 and a7. Word a8 is the extra word (already in r4).
	adds	r6, #4
	ldm	r6!, { r0, r1, r2, r3 }
	@ Finish conditional negation; also, save the original sign in r8.
	asrs	r6, r4, #31
	mov	r8, r6
	eors	r0, r6
	eors	r1, r6
	eors	r2, r6
	eors	r3, r6
	eors	r4, r6
	rsbs	r5, r5, #0
	sbcs	r0, r6
	sbcs	r1, r6
	sbcs	r2, r6
	sbcs	r3, r6
	sbcs	r4, r6
	@ Right-shift words by 15 bits. We have some delayed bits in r7.
	lsrs	r6, r0, #15
	lsls	r0, r0, #17
	orrs	r0, r7
	lsrs	r7, r1, #15
	lsls	r1, r1, #17
	orrs	r1, r6
	lsrs	r6, r2, #15
	lsls	r2, r2, #17
	orrs	r2, r7
	lsrs	r7, r3, #15
	lsls	r3, r3, #17
	orrs	r3, r6
	lsls	r5, r4, #17
	orrs	r5, r7
	@ Write the new a3, a4, a5, a6 and a7.
	mov	r7, r11
	stm	r7!, { r0, r1, r2, r3, r5 }

	@ Return value computation. The original signs of a and b are
	@ in r10 and r8, respectively. Each is either -1 or 0.
	@ We just compute -(r10+2*r8).
	mov	r0, r8
	lsls	r0, r0, #1
	add	r0, r10
	rsbs	r0, r0, #0

	pop	{ pc }
	.size	s256_lin_div15_abs, .-s256_lin_div15_abs

@ =======================================================================
@ void gf_lin2(gf *ab, int32_t f0, int32_t g0, int32_t f1, int32_t g1)
@
@ This routine replaces a and b with, respectively, a*f0+b*g0 and
@ a*f1+b*g1. Computations are done modulo p.
@
@ f0, g0, f1 and g1 are signed integers, in the -2^15..+2^15 range.
@ Moreover, it is known that f0+g0 and f1+g1 are also both in the
@ -2^15..+2^15 range.
@
@ ABI:
@   r0 = pointer to ab
@   r2 = f0
@   r3 = g0
@   r4 = f1
@   r5 = g1
@ No register is preserved.
@
@ Cost: 106 + cost(lin2) = 458
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	gf_lin2_inner, %function
gf_lin2_inner:
	push	{ lr }

	@ Move update factors to high registers (this leaves more room
	@ for computations, and does not imply much overhead since we
	@ have to copy them for each multiplication).
	mov	r8, r2
	mov	r10, r3
	mov	r11, r4
	mov	r12, r5

	@ Perform linear combination.
	bl	lin2

	@ We have two signed results, each with an extra carry (in r4 and
	@ r5, respectively). We apply reduction modulo 2^255-MQ. Each
	@ extra carry is signed, and less than 2^15 in absolute value.

	@ Reduction of the first value.
	ldr	r7, const_gf_lin2_2mq
	muls	r7, r4
	asrs	r6, r4, #31
	mov	r11, r0
	ldm	r0!, { r1, r2, r3, r4 }
	adds	r1, r7
	adcs	r2, r6
	adcs	r3, r6
	adcs	r4, r6
	mov	r0, r11
	stm	r0!, { r1, r2, r3, r4 }
	mov	r11, r0
	ldm	r0!, { r1, r2, r3, r4 }
	adcs	r1, r6
	adcs	r2, r6
	adcs	r3, r6
	adcs	r4, r6
	mov	r0, r11
	stm	r0!, { r1, r2, r3, r4 }
	eors	r1, r1
	adcs	r6, r1
	@ Now r6 contains the extra word of the result.
	@ If the original value was negative, then r6 can only be -1 or 0,
	@ and in the former case, value in a is at least 2^256-2^31.
	@ If the original value was positive, then r6 can only be 0 or 1,
	@ and in the latter case, value in a is at most 2^31. In both
	@ cases, the second reduction will not trigger a carry/borrow
	@ beyond the first limb.
	ldr	r7, const_gf_lin2_2mq
	muls	r6, r7
	subs	r0, #32
	ldr	r1, [r0]
	adds	r1, r6
	str	r1, [r0]

	@ Reduction of the second value.
	adds	r0, #32
	muls	r7, r5
	asrs	r6, r5, #31
	mov	r11, r0
	ldm	r0!, { r1, r2, r3, r4 }
	adds	r1, r7
	adcs	r2, r6
	adcs	r3, r6
	adcs	r4, r6
	mov	r0, r11
	stm	r0!, { r1, r2, r3, r4 }
	mov	r11, r0
	ldm	r0!, { r1, r2, r3, r4 }
	adcs	r1, r6
	adcs	r2, r6
	adcs	r3, r6
	adcs	r4, r6
	mov	r0, r11
	stm	r0!, { r1, r2, r3, r4 }
	eors	r1, r1
	adcs	r6, r1
	ldr	r7, const_gf_lin2_2mq
	muls	r6, r7
	subs	r0, #32
	ldr	r1, [r0]
	adds	r1, r6
	str	r1, [r0]

	pop	{ pc }
	.align	2
const_gf_lin2_2mq:
	.long	2 * MQ
	.size	gf_lin2_inner, .-gf_lin2_inner

@ =======================================================================

@ One step of binary GCD.
@   r0 = xa
@   r1 = xb
@   r2 = fg0
@   r3 = fg1
@ r0, r1, r2 and r3 are updated. r4, r5 and r6 are scratch.
@ Cost: 20
.macro BINGCD_STEP
	@ r4 <- -1 if xa is odd, 0 otherwise
	lsls	r4, r0, #31
	asrs	r4, r4, #31

	@ r5 <- -1 if xa < xb, 0 otherwise
	cmp	r0, r1
	sbcs	r5, r5

	@ Perform conditional swap
	ands	r5, r4
	subs	r6, r0, r1
	ands	r6, r5
	adds	r1, r6
	subs	r0, r6
	subs	r6, r2, r3
	ands	r6, r5
	adds	r3, r6
	subs	r2, r6

	@ Subtract xb from xa if xa is odd
	movs	r5, r1
	ands	r5, r4
	subs	r0, r5
	ands	r4, r3
	subs	r2, r4

	@ Divide xa by 2; multiply f1 and g1 by 2
	lsrs	r0, r0, #1
	lsls	r3, r3, #1
.endm

@ One step of binary Kronecker symbol.
@ This macro works only for the first 13 iterations, because it requires
@ exact knowledge of the three least significant bits of a and b.
@   r0 = xa
@   r1 = xb
@   r2 = fg0
@   r3 = fg1
@   r8 = current symbol value (in least significant bit)
@ r0, r1, r2, r3 and r8 are updated. r4, r5 and r6 are scratch.
@ Cost: 28
.macro BINKRONECKER_STEP
	@ r4 <- -1 if xa is odd, 0 otherwise
	lsls	r4, r0, #31
	asrs	r4, r4, #31

	@ r5 <- -1 if xa and xb must be swapped, 0 otherwise
	cmp	r0, r1
	sbcs	r5, r5
	ands	r5, r4

	@ If swapping, then Kronecker symbol is flipped if and only
	@ if both values are equal to 3 mod 4. This relies on the
	@ following facts:
	@  - if swapping, then both a and b are odd
	@  - it never happens that a < 0 and b < 0 at the same time
	movs	r6, r0
	ands	r6, r1
	lsrs	r6, r6, #1
	ands	r6, r5
	add	r8, r6

	@ Perform conditional swap
	subs	r6, r0, r1
	ands	r6, r5
	adds	r1, r6
	subs	r0, r6
	subs	r6, r2, r3
	ands	r6, r5
	adds	r3, r6
	subs	r2, r6

	@ Subtract xb from xa if xa is odd.
	@ This cannot change the Kronecker symbol, because:
	@  - if b > 0, then (a|b) = (a+k*b|b) for any k \in Z
	@  - if b < 0, then a > 0 and a-b > 0, thus (a|b) = (a-b|b)
	movs	r5, r1
	ands	r5, r4
	subs	r0, r5
	ands	r4, r3
	subs	r2, r4

	@ Divide xa by 2; multiply f1 and g1 by 2
	lsrs	r0, r0, #1
	lsls	r3, r3, #1

	@ Kronecker symbol update: flip if and only if b = 3 or 5 mod 8
	@ at this point.
	adds	r4, r1, #2
	lsrs	r4, r4, #2
	add	r8, r4
.endm

@ One of the two final steps for binary Kronecker symbol. The _updated_
@ low words of a and b are also provided in registers r10 and r11 (for
@ access to the relevant low bits).
@   r0 = xa
@   r1 = xb
@   r2 = fg0
@   r3 = fg1
@   r8 = current symbol value (in least significant bit)
@   r10 = a0
@   r11 = b0
@ r0, r1, r2, r3 and r8 are updated. r4, r5, r6 and r7 are scratch.
@ Cost: 43
.macro BINKRONECKER_STEP_LAST
	@ r4 <- -1 if xa is odd, 0 otherwise
	lsls	r4, r0, #31
	asrs	r4, r4, #31

	@ r5 <- -1 if xa and xb must be swapped, 0 otherwise
	cmp	r0, r1
	sbcs	r5, r5
	ands	r5, r4

	@ If swapping, then Kronecker symbol is flipped if and only
	@ if both values are equal to 3 mod 4. We use the two
	@ values in r10 and r11.
	mov	r6, r10
	mov	r7, r11
	ands	r6, r7
	lsrs	r6, r6, #1
	ands	r6, r5
	add	r8, r6

	@ Perform conditional swap
	subs	r6, r0, r1
	ands	r6, r5
	adds	r1, r6
	subs	r0, r6
	subs	r6, r2, r3
	ands	r6, r5
	adds	r3, r6
	subs	r2, r6
	@ Also on values in r10 and r11
	mov	r6, r10
	mov	r7, r11
	subs	r6, r7
	ands	r6, r5
	rsbs	r7, r6, #0
	add	r10, r7
	add	r11, r6

	@ Subtract xb from xa if xa is odd
	movs	r5, r1
	ands	r5, r4
	subs	r0, r5
	movs	r5, r3
	ands	r5, r4
	subs	r2, r5

	@ Divide xa by 2; multiply f1 and g1 by 2
	lsrs	r0, r0, #1
	lsls	r3, r3, #1

	@ Also do conditional subtraction, and division by 2, on
	@ values in r10 and r11.
	mov	r6, r10
	mov	r7, r11
	ands	r4, r7
	subs	r6, r4
	lsrs	r6, r6, #1
	mov	r10, r6

	@ Kronecker symbol update: flip if and only if b = 3 or 5 mod 8
	@ at this point.
	adds	r4, r7, #2
	lsrs	r4, r4, #2
	add	r8, r4
.endm

@ =======================================================================
@ bingcd15()
@
@ Inputs:
@   r0 = xa
@   r1 = xb
@
@ Outputs:
@   r0 = new xa
@   r1 = new xb
@   r2 = f0
@   r3 = f1
@   r4 = g0
@   r5 = g1
@
@ r6 is scratch. r14 is consumed. All other registers are unmodified.
@
@ Cost: 317
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	bingcd15, %function
bingcd15:
	movs	r2, #1
	lsls	r3, r2, #16
	BINGCD_STEP
	BINGCD_STEP
	BINGCD_STEP
	BINGCD_STEP
	BINGCD_STEP
	BINGCD_STEP
	BINGCD_STEP
	BINGCD_STEP
	BINGCD_STEP
	BINGCD_STEP
	BINGCD_STEP
	BINGCD_STEP
	BINGCD_STEP
	BINGCD_STEP
	BINGCD_STEP

	@ Unpack the update factors.
	ldr	r6, const_bingcd15_unpack
	adds	r2, r6
	adds	r3, r6
	lsrs	r5, r3, #16
	uxth	r4, r3
	lsrs	r3, r2, #16
	uxth	r2, r2
	uxth	r6, r6
	subs	r5, r6
	subs	r4, r6
	subs	r3, r6
	subs	r2, r6

	bx	lr
	.align	2
const_bingcd15_unpack:
	.long	0x7FFF7FFF
	.size	bingcd15, .-bingcd15

@ =======================================================================
@ binkronecker15()
@
@ Inputs:
@   r0 = xa
@   r1 = xb
@   r10 = a0
@   r11 = b0
@
@ Outputs:
@   r0 = new xa
@   r1 = new xb
@   r2 = f0
@   r3 = f1
@   r4 = g0
@   r5 = g1
@   r6 = Kronecker symbol update (lsb)
@
@ r7, r8, r10, r11 and r14 are scratch. r12 is unmodified.
@
@ Cost: 49 + 13*cost(BINKRONECKER_STEP) + 2*cost(BINKRONECKER_STEP_LAST)
@       = 499
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	binkronecker15, %function
binkronecker15:
	movs	r2, #1
	lsls	r3, r2, #16
	eors	r4, r4
	mov	r8, r4
	BINKRONECKER_STEP
	BINKRONECKER_STEP
	BINKRONECKER_STEP
	BINKRONECKER_STEP
	BINKRONECKER_STEP
	BINKRONECKER_STEP
	BINKRONECKER_STEP
	BINKRONECKER_STEP
	BINKRONECKER_STEP
	BINKRONECKER_STEP
	BINKRONECKER_STEP
	BINKRONECKER_STEP
	BINKRONECKER_STEP

	@ For the two last steps, we need to compute a few more low bits for
	@ a and b. We have the two original a0 and b0 in r10 and r11,
	@ respectively, so we can apply the update factors on them.

	@ Unpack f0 and g0
	ldr	r7, const_binkronecker15_unpack
	adds	r4, r2, r7
	lsrs	r5, r4, #16
	uxth	r4, r4
	uxth	r7, r7
	subs	r4, r7
	subs	r5, r7
	@ Apply f0 and g0; also, keep original a0 in r6.
	mov	r6, r10
	muls	r4, r6
	mov	r7, r11
	muls	r5, r7
	adds	r4, r5
	lsrs	r4, r4, #13
	mov	r10, r4
	@ Unpack f1 and g1
	ldr	r7, const_binkronecker15_unpack
	adds	r4, r3, r7
	lsrs	r5, r4, #16
	uxth	r4, r4
	uxth	r7, r7
	subs	r4, r7
	subs	r5, r7
	@ Apply f1 and g1
	muls	r4, r6
	mov	r7, r11
	muls	r5, r7
	adds	r4, r5
	lsrs	r4, r4, #13
	mov	r11, r4

	BINKRONECKER_STEP_LAST
	BINKRONECKER_STEP_LAST

	@ Unpack the update factors.
	ldr	r6, const_binkronecker15_unpack
	adds	r2, r6
	adds	r3, r6
	lsrs	r5, r3, #16
	uxth	r4, r3
	lsrs	r3, r2, #16
	uxth	r2, r2
	uxth	r6, r6
	subs	r5, r6
	subs	r4, r6
	subs	r3, r6
	subs	r2, r6

	@ Return the Kronecker symbol updates.
	mov	r6, r8

	bx	lr
	.align	2
const_binkronecker15_unpack:
	.long	0x7FFF7FFF
	.size	binkronecker15, .-binkronecker15

@ =======================================================================
@ uint32_t gf_inv(gf *d, const gf *y)
@
@ ABI: all registers are consumed.
@
@ Cost: 1568 + cost(gf_normalize_inner) + 33*cost(approximate_ab)
@       + 34*cost(bingcd15) + 33*cost(s256_lin_div15_abs)
@       + 34*cost(gf_lin2_inner) + cost(gf_mul_inner) + cost(gf_iszero_inner)
@       = 54793
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	gf_inv_inner, %function
gf_inv_inner:
	push	{ r0, lr }
	sub	sp, #136

	@ Initial values:
	@ name  index  value
	@   a      0     y (normalized)
	@   b      8     p (modulus)
	@   u     16     1
	@   v     24     0

	@ b <- 2^255-MQ
	add	r0, sp, #32
	ldr	r2, const_gf_inv_minusmq
	asrs	r3, r2, #31
	movs	r4, r3
	movs	r5, r3
	stm	r0!, { r2, r3, r4, r5 }
	movs	r2, r3
	lsrs	r5, r5, #1
	stm	r0!, { r2, r3, r4, r5 }

	@ u <- 1
	movs	r2, #1
	eors	r3, r3
	eors	r4, r4
	eors	r5, r5
	stm	r0!, { r2, r3, r4, r5 }
	eors	r2, r2
	stm	r0!, { r2, r3, r4, r5 }

	@ v <- 0
	stm	r0!, { r2, r3, r4, r5 }
	stm	r0!, { r2, r3, r4, r5 }

	@ a <- y (normalized)
	mov	r0, sp
	bl	gf_normalize_inner

	@ Perform 33*15 = 495 iterations. Loop counter is kept on
	@ the stack.
	movs	r0, #33
	str	r0, [sp, #128]
.Lgf_inv_inner:
	@ Get approximations of a and b in r0 and r1
	mov	r0, sp
	add	r1, sp, #32
	bl	approximate_ab

	@ Run 15 inner iterations of the binary GCD
	bl	bingcd15

	@ Apply update factors on a and b. Take care to save the
	@ update factors.
	mov	r0, sp
	push	{ r2, r3, r4, r5 }
	bl	s256_lin_div15_abs

	@ Restore the update factors. If negation was applied on
	@ a and/or b, then propagate it to the relevant update
	@ factors.
	pop	{ r2, r3, r4, r5 }
	lsls	r6, r0, #31
	asrs	r6, r6, #31
	eors	r2, r6
	subs	r2, r6
	eors	r3, r6
	subs	r3, r6
	lsls	r6, r0, #30
	asrs	r6, r6, #31
	eors	r4, r6
	subs	r4, r6
	eors	r5, r6
	subs	r5, r6

	@ Apply update factors to u and v.
	add	r0, sp, #64
	bl	gf_lin2_inner

	@ Loop control.
	ldr	r7, [sp, #128]
	subs	r7, #1
	str	r7, [sp, #128]
	bne	.Lgf_inv_inner

	@ We have finished the main loop. 495 inner iterations have
	@ been performed. Values are now at most 15 bits in length;
	@ we can perform the final 15 iterations without any
	@ approximation.
	ldr	r0, [sp]
	ldr	r1, [sp, #32]
	bl	bingcd15

	@ TODO: we only need to update v here, not u.
	add	r0, sp, #64
	bl	gf_lin2_inner

	@ We have 510 delayed halvings to perform.
	ldr	r0, [sp, #136]
	add	r1, sp, #96
	adr	r2, const_gf_inv_invt510
	bl	gf_mul_inner

	@ API asks for a returned value of 1 on success, 0 on error.
	@ Algorithm above always succeeds, except when the input is not
	@ invertible modulo p, i.e. when it is zero (since p is prime).
	@ In the latter case, we got an output of zero, which is what
	@ we want; we can thus compare the output with zero to obtain
	@ the proper status to return.
	ldr	r0, [sp, #136]
	bl	gf_iszero_inner
	movs	r1, #1
	eors	r0, r1

	add	sp, #140
	pop	{ pc }
	.align	2
const_gf_inv_minusmq:
	.long	-MQ
const_gf_inv_invt510:
	.long	INVT510
	.size	gf_inv_inner, .-gf_inv_inner

@ =======================================================================
@ int32_t gf_legendre(const gf *y)
@
@ ABI: all registers are consumed.
@
@ Cost: 1349 + cost(gf_normalize_inner)
@       + 33*cost(approximate_ab) + 34*cost(binkronecker15)
@       + 33*cost(s256_lin_div15_abs) + cost(gf_iszero_inner)
@       = 43726
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	gf_legendre_inner, %function
gf_legendre_inner:
	push	{ r0, lr }
	sub	sp, #72

	@ Initial values:
	@ name  index  value
	@   a      0     y (normalized)
	@   b      8     p (modulus)

	@ a <- y (normalized)
	movs	r1, r0
	mov	r0, sp
	bl	gf_normalize_inner

	@ b <- 2^255-MQ
	add	r0, sp, #32
	ldr	r2, const_gf_legendre_minusmq
	asrs	r3, r2, #31
	movs	r4, r3
	movs	r5, r3
	stm	r0!, { r2, r3, r4, r5 }
	movs	r2, r3
	lsrs	r5, r5, #1
	stm	r0!, { r2, r3, r4, r5 }

	@ Initialize Kronecker symbol accumulator.
	eors	r0, r0
	str	r0, [sp, #68]

	@ Perform 33*15 = 495 iterations. Loop counter is kept on
	@ the stack.
	movs	r0, #33
	str	r0, [sp, #64]
.Lgf_legendre_inner:
	@ Get approximations of a and b in r0 and r1
	mov	r0, sp
	add	r1, sp, #32
	bl	approximate_ab

	@ Run 15 inner iterations of the binary Kronecker symbol
	ldr	r2, [sp]
	mov	r10, r2
	ldr	r3, [sp, #32]
	mov	r11, r3
	bl	binkronecker15
	ldr	r7, [sp, #68]
	eors	r7, r6
	str	r7, [sp, #68]

	@ Apply update factors on a and b. Take care to save the
	@ update factors.
	mov	r0, sp
	bl	s256_lin_div15_abs

	@ We may had to negate a or b (never both, though). Negating
	@ b does not change the result (because, in that case, a > 0).
	@ However, negating a flips the symbol if and only if b = 3 mod 4.
	ldr	r1, [sp, #32]
	lsrs	r1, r1, #1
	ands	r0, r1
	ldr	r2, [sp, #68]
	eors	r2, r0
	str	r2, [sp, #68]

	@ Loop control.
	ldr	r7, [sp, #64]
	subs	r7, #1
	str	r7, [sp, #64]
	bne	.Lgf_legendre_inner

	@ We have finished the main loop. 495 inner iterations have
	@ been performed. Values are now at most 15 bits in length;
	@ we can perform the final 15 iterations without any
	@ approximation.
	ldr	r0, [sp]
	ldr	r1, [sp, #32]
	mov	r10, r0
	mov	r11, r1
	bl	binkronecker15

	@ Since the values were exact, there will be no negation on
	@ a and b.
	ldr	r1, [sp, #68]
	eors	r1, r6
	movs	r0, #1
	ands	r1, r0
	lsls	r1, r1, #1
	subs	r0, r1
	mov	r8, r0

	@ Value is correct, unless the source operand was 0, in which
	@ case we need to return 0 and not the value computed above
	@ (as per the API).
	ldr	r0, [sp, #72]
	bl	gf_iszero_inner
	subs	r1, r0, #1
	mov	r0, r8
	ands	r0, r1

	add	sp, #76
	pop	{ pc }
	.align	2
const_gf_legendre_minusmq:
	.long	-MQ
	.size	gf_legendre_inner, .-gf_legendre_inner

@ =======================================================================
@ void cswap(gf *a, gf *b, uint32_t cn)
@
@ Contents of a and b are swapped if cn = -1; they are not modified if
@ cn = 0. Value of cn MUST be 0 or -1.
@
@ Registers r0 to r7 are modified (except r2). High registers are untouched.
@
@ Cost: 90
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	cswap, %function
cswap:
	ldm	r0!, { r3, r4 }
	ldm	r1!, { r5, r6 }
	subs	r7, r3, r5
	ands	r7, r2
	subs	r3, r7
	adds	r5, r7
	subs	r7, r4, r6
	ands	r7, r2
	subs	r4, r7
	adds	r6, r7
	subs	r0, #8
	subs	r1, #8
	stm	r0!, { r3, r4 }
	stm	r1!, { r5, r6 }
	ldm	r0!, { r3, r4 }
	ldm	r1!, { r5, r6 }
	subs	r7, r3, r5
	ands	r7, r2
	subs	r3, r7
	adds	r5, r7
	subs	r7, r4, r6
	ands	r7, r2
	subs	r4, r7
	adds	r6, r7
	subs	r0, #8
	subs	r1, #8
	stm	r0!, { r3, r4 }
	stm	r1!, { r5, r6 }
	ldm	r0!, { r3, r4 }
	ldm	r1!, { r5, r6 }
	subs	r7, r3, r5
	ands	r7, r2
	subs	r3, r7
	adds	r5, r7
	subs	r7, r4, r6
	ands	r7, r2
	subs	r4, r7
	adds	r6, r7
	subs	r0, #8
	subs	r1, #8
	stm	r0!, { r3, r4 }
	stm	r1!, { r5, r6 }
	ldm	r0!, { r3, r4 }
	ldm	r1!, { r5, r6 }
	subs	r7, r3, r5
	ands	r7, r2
	subs	r3, r7
	adds	r5, r7
	subs	r7, r4, r6
	ands	r7, r2
	subs	r4, r7
	adds	r6, r7
	subs	r0, #8
	subs	r1, #8
	stm	r0!, { r3, r4 }
	stm	r1!, { r5, r6 }
	bx	lr
	.size	cswap, .-cswap

@ =======================================================================

@ One step of gf_mul_a24_inner().
@   r1   pointer to next input word (updated)
@   r0   pointer to next output word (updated)
@   r6   carry word
@   r8   constant a24 = 121665
@ Registers r2, r3, r4, r5 and r7 are scratch.
@ Cost: 25
.macro MULA24_STEP
	ldm	r1!, { r2 }
	mov	r7, r8
	MUL32x32  r2, r3, r7, r4, r5
	adds	r2, r6
	eors	r6, r6
	adcs	r6, r3
	stm	r0!, { r2 }
.endm

@ =======================================================================
@ void gf_mul_a24(gf *d, const gf *a)
@
@ Multiply a value by the small constant a24 = 121665.
@
@ Registers r0 to r11 are modified. r12 is untouched.
@
@ Cost: 249
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	gf_mul_a24_inner, %function
gf_mul_a24_inner:
	ldr	r7, const_mul_a24_a24
	mov	r8, r7
	eors	r6, r6

	MULA24_STEP
	MULA24_STEP
	MULA24_STEP
	MULA24_STEP
	MULA24_STEP
	MULA24_STEP
	MULA24_STEP
	MULA24_STEP

	@ Extra word is in r6; it cannot be larger than 121664; we
	@ can multiply it by 38, and the result will still fit in 32 bits.
	movs	r7, #38
	muls	r7, r6

	@ Do first round of reduction.
	subs	r0, #32
	ldm	r0!, { r1, r2, r3, r4, r5, r6 }
	adds	r1, r7
	eors	r7, r7
	adcs	r2, r7
	adcs	r3, r7
	adcs	r4, r7
	adcs	r5, r7
	adcs	r6, r7
	mov	r8, r5
	mov	r10, r6
	ldm	r0!, { r5, r6 }
	adcs	r5, r7
	adcs	r6, r7
	adcs	r7, r7

	@ Value is in r1:r2:r3:r4:r8:r10:r5:r6, with an extra bit (0 or 1)
	@ in r7. If r7 is distinct from 0, then we must multiply it by 38
	@ and add it to r1; however, that addition will not trigger a
	@ carry.
	mov	r11, r5
	movs	r5, #38
	muls	r7, r5
	adds	r1, r7

	@ Write final result.
	subs	r0, #32
	stm	r0!, { r1, r2, r3, r4 }
	mov	r3, r8
	mov	r4, r10
	mov	r5, r11
	stm	r0!, { r3, r4, r5, r6 }

	bx	lr
	.align	2
const_mul_a24_a24:
	.long	121665
	.size	gf_mul_a24_inner, .-gf_mul_a24_inner

@ =======================================================================
@ void x25519(void *dst, const void *src, const void *scalar)
@
@ This function follows the external ABI and is callable from C code.
@ =======================================================================

	.align	1
	.global	x25519
	.thumb
	.thumb_func
	.type	x25519, %function
x25519:
	push	{ r4, r5, r6, r7, lr }
	mov	r4, r8
	mov	r5, r10
	mov	r6, r11
	push	{ r4, r5, r6 }

	@ Stack layout:
	@   k     0     decoded scalar (clamped)
	@   x1    32    (x2, z2, x3 and z3 must be consecutive, in that order)
	@   x2    64
	@   z2    96
	@   x3    128
	@   z3    160
	@   A     192
	@   AA    224
	@   B     256
	@   BB    288
	@   C           (alias on x2)
	@   D           (alias on z2)
	@   CB          (alias on x2)
	@   DA          (alias on z2)
	@   E           (alias on x2)
	@   t     320   loop counter (254 to 0)
	@   dst   328   destination address
	@   src   332   source point address
	@   swap        "swap" flag (0 or -1) (same address as src)
	push	{ r0, r1 }
	sub	sp, #328

#define x1     32
#define x2     64
#define z2     96
#define x3     128
#define z3     160
#define tA     192
#define tAA    224
#define tB     256
#define tBB    288
#define tC     x2
#define tD     z2
#define tCB    x2
#define tDA    z2
#define tE     x2
#define dst    328
#define src    332
#define swap   src
#define tcc    320

	@ Copy scalar into stack buffer, then clamp it.
	mov	r0, sp
	movs	r1, r2
	movs	r2, #32
	bl	memcpy
	ldr	r0, [sp]
	lsrs	r0, r0, #3
	lsls	r0, r0, #3
	str	r0, [sp]
	ldr	r0, [sp, #28]
	lsls	r0, r0, #1
	lsrs	r0, r0, #1
	movs	r1, #1
	lsls	r1, r1, #30
	orrs	r0, r1
	str	r0, [sp, #28]

	@ Decode source value, ignoring the top bit of the last byte.
	add	r0, sp, #x1
	ldr	r1, [sp, #src]
	movs	r2, #32
	bl	memcpy
	ldr	r0, [sp, #(x1 + 28)]
	lsls	r0, r0, #1
	lsrs	r0, r0, #1
	str	r0, [sp, #(x1 + 28)]

	@ Init:
	@   x2 = 1
	@   z2 = 0
	@   x3 = x1
	@   z3 = 1
	@   swap = 0
	movs	r1, #1
	eors	r2, r2
	eors	r3, r3
	eors	r4, r4
	eors	r5, r5
	add	r0, sp, #x2
	stm	r0!, { r1, r2, r3, r4 }
	stm	r0!, { r2, r3, r4, r5 }
	stm	r0!, { r2, r3, r4, r5 }
	stm	r0!, { r2, r3, r4, r5 }
	adds	r0, #32
	stm	r0!, { r1, r2, r3, r4 }
	stm	r0!, { r2, r3, r4, r5 }
	subs	r0, #64
	add	r1, sp, #x1
	ldm	r1!, { r2, r3, r4, r5 }
	stm	r0!, { r2, r3, r4, r5 }
	ldm	r1!, { r2, r3, r4, r5 }
	stm	r0!, { r2, r3, r4, r5 }
	eors	r0, r0
	str	r0, [sp, #swap]

	@ Main loop: 255 iterations.
	@ Each loop iteration expects the counter in register r2.
	movs	r2, #254
.Lx25519_loop:
	str	r2, [sp, #tcc]

	@ Get next scalar bit (kt) into r0 (into 0/-1 convention).
	lsrs	r1, r2, #3
	mov	r0, sp
	ldrb	r0, [r0, r1]
	movs	r1, #7
	ands	r1, r2
	lsrs	r0, r1
	lsls	r0, r0, #31
	asrs	r0, r0, #31

	@ Do conditional swaps and update the 'swap' flag.
	ldr	r2, [sp, #swap]
	eors	r2, r0
	str	r0, [sp, #swap]
	add	r0, sp, #x2
	add	r1, sp, #x3
	bl	cswap
	add	r0, sp, #z2
	add	r1, sp, #z3
	bl	cswap

	@ A <- x2 + z2
	add	r0, sp, #tA
	add	r1, sp, #x2
	add	r2, sp, #z2
	bl	gf_add_inner

	@ AA <- A^2
	add	r0, sp, #tAA
	add	r1, sp, #tA
	bl	gf_sqr_inner

	@ B <- x2 - z2
	add	r0, sp, #tB
	add	r1, sp, #x2
	add	r2, sp, #z2
	bl	gf_sub_inner

	@ BB <- B^2
	add	r0, sp, #tBB
	add	r1, sp, #tB
	bl	gf_sqr_inner

	@ C <- x3 + z3
	add	r0, sp, #tC
	add	r1, sp, #x3
	add	r2, sp, #z3
	bl	gf_add_inner

	@ D <- x3 - z3
	add	r0, sp, #tD
	add	r1, sp, #x3
	add	r2, sp, #z3
	bl	gf_sub_inner

	@ CB <- C * B
	add	r0, sp, #tCB
	add	r1, sp, #tC
	add	r2, sp, #tB
	bl	gf_mul_inner

	@ DA <- D * A
	add	r0, sp, #tDA
	add	r1, sp, #tD
	add	r2, sp, #tA
	bl	gf_mul_inner

	@ x_3 <- (DA + CB)^2
	add	r0, sp, #x3
	add	r1, sp, #tDA
	add	r2, sp, #tCB
	bl	gf_add_inner
	add	r0, sp, #x3
	add	r1, sp, #x3
	bl	gf_sqr_inner

	@ z3 <- x_1 * (DA - CB)^2
	add	r0, sp, #z3
	add	r1, sp, #tDA
	add	r2, sp, #tCB
	bl	gf_sub_inner
	add	r0, sp, #z3
	add	r1, sp, #z3
	bl	gf_sqr_inner
	add	r0, sp, #z3
	add	r1, sp, #z3
	add	r2, sp, #x1
	bl	gf_mul_inner

	@ E <- AA - BB
	add	r0, sp, #tE
	add	r1, sp, #tAA
	add	r2, sp, #tBB
	bl	gf_sub_inner

	@ z2 <- E * (AA + a24 * E)
	add	r0, sp, #z2
	add	r1, sp, #tE
	bl	gf_mul_a24_inner
	add	r0, sp, #z2
	add	r1, sp, #z2
	add	r2, sp, #tAA
	bl	gf_add_inner
	add	r0, sp, #z2
	add	r1, sp, #z2
	add	r2, sp, #tE
	bl	gf_mul_inner

	@ x2 <- AA * BB
	add	r0, sp, #x2
	add	r1, sp, #tAA
	add	r2, sp, #tBB
	bl	gf_mul_inner

	@ Loop control
	ldr	r2, [sp, #tcc]
	subs	r2, #1
	bcs	.Lx25519_loop

	@ Main loop is finished. We still have a pair of conditional
	@ swaps to perform.
	ldr	r2, [sp, #swap]
	add	r0, sp, #x2
	add	r1, sp, #x3
	bl	cswap
	add	r0, sp, #z2
	add	r1, sp, #z3
	bl	cswap

	@ Result is x2/z2; if z2 == 0 then result is zero. The gf_inv_inner()
	@ function will invert zero into zero, which is correct.
	add	r0, sp, #z2
	add	r1, sp, #z2
	bl	gf_inv_inner
	add	r0, sp, #x2
	add	r1, sp, #x2
	add	r2, sp, #z2
	bl	gf_mul_inner

	@ Encode the result, with proper normalization into the 0..p-1 range.
	@ Encoding must be into a stack buffer, to avoid alignment issues;
	@ memcpy() is then used to write into the output buffer.
	add	r0, sp, #x1
	add	r1, sp, #x2
	bl	gf_normalize_inner
	ldr	r0, [sp, #dst]
	add	r1, sp, #x1
	movs	r2, #32
	bl	memcpy

	add	sp, #336

	pop	{ r4, r5, r6 }
	mov	r8, r4
	mov	r10, r5
	mov	r11, r6
	pop	{ r4, r5, r6, r7, pc }

#undef x1
#undef x2
#undef z2
#undef x3
#undef z3
#undef tA
#undef tAA
#undef tB
#undef tBB
#undef tC
#undef tD
#undef tCB
#undef tDA
#undef tE
#undef dst
#undef src
#undef swap
#undef tcc

	.size	x25519, .-x25519

@ =======================================================================
@ Public wrappers used only for tests; commented out by default.
@
@@ =======================================================================
@@ void test_gf_add(gf *c, const gf *a, const gf *b)
@@ =======================================================================
@
@	.align	1
@	.global	test_gf_add
@	.thumb
@	.thumb_func
@	.type	test_gf_add, %function
@test_gf_add:
@	push	{ r4, r5, r6, r7, lr }
@	mov	r4, r8
@	mov	r5, r10
@	mov	r6, r11
@	push	{ r4, r5, r6 }
@	bl	gf_add_inner
@	pop	{ r4, r5, r6 }
@	mov	r8, r4
@	mov	r10, r5
@	mov	r11, r6
@	pop	{ r4, r5, r6, r7, pc }
@	.size	test_gf_add, .-test_gf_add
@
@@ =======================================================================
@@ int32_t test_gf_legendre(const gf *a)
@@ =======================================================================
@
@	.align	1
@	.global	test_gf_legendre
@	.thumb
@	.thumb_func
@	.type	test_gf_legendre, %function
@test_gf_legendre:
@	push	{ r4, r5, r6, r7, lr }
@	mov	r4, r8
@	mov	r5, r10
@	mov	r6, r11
@	push	{ r4, r5, r6 }
@	bl	gf_legendre_inner
@	pop	{ r4, r5, r6 }
@	mov	r8, r4
@	mov	r10, r5
@	mov	r11, r6
@	pop	{ r4, r5, r6, r7, pc }
@	.size	test_gf_legendre, .-test_gf_legendre
