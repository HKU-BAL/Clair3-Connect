#!/usr/bin/env bash

export LIBTORCH="${LIBTORCH:-/path/to/libtorch}"
export LD_LIBRARY_PATH="$LIBTORCH/lib:${LD_LIBRARY_PATH:-}"
export LIBTORCH_BYPASS_VERSION_CHECK=1

export CC=/usr/bin/gcc
export CXX=/usr/bin/g++
export CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_LINKER=/usr/bin/gcc

export PROTOC="${PROTOC:-$(command -v protoc 2>/dev/null)}"
[ -n "$PROTOC" ] || echo "  WARNING: protoc not found — install conda 'libprotobuf' or system 'protobuf-compiler'."

export PYTORCH_TENSOREXPR=0

echo "Clair3 Connect env ready:"
echo "  LIBTORCH=$LIBTORCH"
echo "  PROTOC=$PROTOC"
[ -e "$LIBTORCH/lib/libc10.so" ] || echo "  WARNING: \$LIBTORCH/lib/libc10.so not found — set LIBTORCH correctly (server only)."
case "$(ldd --version 2>/dev/null | head -1)" in
  *2.2[0-8]*|*2.1*) echo "  WARNING: system glibc may be < 2.29; LibTorch 2.5 needs >= 2.29.";;
esac
