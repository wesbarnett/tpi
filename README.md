# Test particle insertion

Personal code that does test particle insertion for trajectories saved with
GROMACS. Uses some SIMD intrinsics to make it fast. Requires libgmxcpp 5.0+ with
AVX instructions enabled.

To install do:

    make
    make PREFIX=/usr/local install

Change `PREFIX` to your liking. The above installs the binary to
`/usr/local/bin/tpi`.

An example configuration input file is found in `dat`.

Was this software helpful? Send me a BTC tip: `1PZziQoUJfhMKZC8gXQZtS5ebHWMba3Geb`
