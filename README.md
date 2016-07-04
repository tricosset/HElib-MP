# HElib-MP

HElib-MP is an adaptation of [HElib] [1] to handle multiprecision plaintext moduli.

At its present state, this library provides only basic SHE features
(key generation, encryption, decryption, addition and multiplication).

  [1]: https://github.com/shaih/HElib         "HElib"

## Tests

Three tests are provided for the library

- **src/Test_SHE**: this program benchmarks the encryption, homomorphic multiplication and decryption.
- **src/Test_RSA**: this program computes homomorphically a modular exponentiation with a 2048-bit plaintext modulus.
- **src/Test_ECC**: this program computes homomorphically an elliptic curve point multiplication over the NIST P-256 curve.

## License

GPLv2
