# ConSTRain
**Co**py **n**umber-guided **STR** **a**llele **in**ference

## Installation
### Building from source
In order to build ConSTRain from the source code, you will need to [install the Rust programming language](https://www.rust-lang.org/tools/install).
Then, you can clone this repository and build ConSTRain as follows:

```bash
git clone https://github.com/acg-team/ConSTRain.git
cd ConSTRain/ConSTRain
cargo build --release --bin ConSTRain
```

The binary will be created at `./target/release/ConSTRain`.

### Known issues
#### libclang is not found
Building `hts-sys` fails with the following error: 
```bash
Unable to find libclang: "couldn't find any valid shared libraries matching: ['libclang.so', 'libclang-*.so', 'libclang.so.*' 'libclang-*.so.*'], set the LIBCLANG_PATH environment variable to a path where one of these files can be found (invalid: [])"
```

This is caused because `Bindgen` requires `libclang`. The problem can be solved by installing `clang` (see instructions [here](https://rust-lang.github.io/rust-bindgen/requirements.html)) or, if it is already installed, providing the location of the directory containing `libclang.so` through the `LIBCLANG_PATH` environment variable:
```bash
LIBCLANG_PATH=/path/to/clang/lib cargo build --release --bin ConSTRain
```

#### System header files are not found
Building `hts-sys` fails with an error like this:
```bash
./htslib/htslib/hts.h:31:10: fatal error: 'stddef.h' file not found
```

Caused because some system header file (e.g., `stddef.h`) is not found. Address this by providing the location of the directory containing system header files through the `C_INCLUDE_PATH` environment variable:
```bash
C_INCLUDE_PATH=/path/to/include cargo build --release --bin ConSTRain
```
