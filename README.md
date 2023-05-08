# babyjubjub-ark [![Test](https://github.com/arnaucube/babyjubjub-ark/workflows/Test/badge.svg)](https://github.com/arnaucube/babyjubjub-ark/actions?query=workflow%3ATest)

> **Note**: this repo is a fork from https://github.com/arnaucube/babyjubjub-rs , porting it to arkworks [ff](https://github.com/arkworks-rs/algebra/tree/master/ff).

BabyJubJub elliptic curve implementation in Rust. A twisted edwards curve embedded in the curve of BN128/BN256.

BabyJubJub curve explanation: https://medium.com/zokrates/efficient-ecc-in-zksnarks-using-zokrates-bd9ae37b8186

Uses:
- Poseidon hash function https://github.com/arnaucube/poseidon-ark

Compatible with the BabyJubJub implementations in:
- Go, from https://github.com/iden3/go-iden3-crypto
- circom & javascript, from https://github.com/iden3/circomlib

## Warning
Doing this in my free time, **do not use in production**.

### References
- BabyJubJub curve explanation: https://medium.com/zokrates/efficient-ecc-in-zksnarks-using-zokrates-bd9ae37b8186
	- C++ & Explanation https://github.com/barryWhiteHat/baby_jubjub
		- C++ https://github.com/barryWhiteHat/baby_jubjub_ecc
	- Javascript & Circom: https://github.com/iden3/circomlib
	- Go https://github.com/iden3/go-iden3-crypto
- JubJub curve explanation: https://z.cash/technology/jubjub/
	- Rust: https://github.com/zkcrypto/jubjub
	- Python: https://github.com/daira/jubjub
