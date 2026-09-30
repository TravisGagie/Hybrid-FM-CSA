# psie: an FM-index / CSA hybrid for repetitive texts over large alphabets

This is the prototype accompanying the note *Mixing FM-indexes and CSAs: backward search over an order-1 rank
encoding* (Travis Gagie). It counts occurrences of patterns in a sequence of integers.

## The idea in brief

Each character of the text T is replaced by its rank by frequency among the characters that follow the previous
character, giving an encoding E over a small, skewed alphabet. A pattern P is encoded the same way and searched for
in an FM-type index of E. The encoding of P only determines P given its first character, so the search ends with one
CSA-style step on Psi_E, the array of unencoded characters that precede the encoded suffixes. Counting is exact.
See the note for details, experiments and caveats.

## Building

Requires a C++17 compiler and make. The code uses BMI2 instructions when compiled with `-march=native` on a CPU that
has them, and portable fallbacks otherwise.

    make

## Usage

Inputs are binary files of little-endian unsigned integers, 4 bytes each by default.

    ./psie build <input> <index> [width] [huff|ef|explicit]
    ./psie count <index> <patterns.txt>
    ./psie bench <input> [width] [queries] [m ...]

- `build` builds the hybrid index and writes it to `<index>`. `width` is the number of bytes per integer (1, 2, 4
  or 8). The last argument chooses the structure for BWT(E): `huff` (default, run-length Huffman-shaped wavelet
  tree), `ef` (Elias-Fano run lists) or `explicit` (uncompressed run lists, fastest, largest).
- `count` reads one pattern per line, as whitespace-separated decimal integers, and prints the number of
  occurrences of each.
- `bench` builds eight indexes in memory: the hybrid with each of four structures on BWT(E), and an index of T with
  each of the same four structures on BWT(T) (an RLCSA, a compressed RLCSA, an FM-index and an RLFM-index). It
  checks that they agree, and prints their sizes and counting times for random substrings of the input of each
  length `m`.

Synthetic test data like that in the note can be generated with

    python3 gen.py <out.bin> <sigma> <base_len> <copies> <mutation_rate>

for example `python3 gen.py test.bin 1000 1000000 50 0.01` (requires numpy).

## Files

- `psie.cpp`: encoding, index construction, counting, and the `build`, `count` and `bench` modes.
- `runrank.hpp`: explicit per-character run lists (an uncompressed RLCSA).
- `efrank.hpp`: Elias-Fano run lists (a compressed RLCSA, following Brown, Gagie, Manzini, Navarro and Sciortino,
  "Faster run-length compressed suffix arrays").
- `huffwt.hpp`: Huffman-shaped wavelet trees, plain and run-length (an RLFM-index).
- `gen.py`: generator for synthetic repetitive data.
- `libsais/`: Ilya Grebnov's libsais, used to build suffix arrays, included under its Apache 2.0 licence.

## Limitations

- Counting only; locating is not implemented.
- Inputs must have fewer than 2^31 integers and at most 2^32 distinct values.
- This is research code and has not been tuned; in particular the RLFM-index is slower than a careful
  implementation would be.

## Licence

GNU General Public License v3.0 or later; see `LICENSE`. The bundled libsais is under the Apache License 2.0; see
`libsais/LICENSE`.
