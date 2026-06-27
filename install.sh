#!/usr/bin/env bash

set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
RUST="$HERE/rust"
TORCH_VER="2.5.0"
WHAT="${1:-all}"
CUDA_SPEC="${2:-auto}"

say() { printf '\033[1;36m[install]\033[0m %s\n' "$*" >&2; }
die() { printf '\033[1;31m[install] ERROR:\033[0m %s\n' "$*" >&2; exit 1; }

if ! command -v cargo >/dev/null 2>&1; then
  say "Rust not found — installing rustup (non-interactive)…"
  curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
  source "$HOME/.cargo/env"
fi
command -v cargo >/dev/null 2>&1 || die "cargo still not on PATH; run: source \$HOME/.cargo/env"

IN_CONDA="${CONDA_PREFIX:-}"
need=()
command -v cc >/dev/null 2>&1 || command -v gcc >/dev/null 2>&1 || need+=("a C compiler (gcc)")
command -v protoc >/dev/null 2>&1 || need+=("protoc (protobuf compiler)")
if [ "${#need[@]}" -gt 0 ]; then
  echo
  say "Missing build dependencies: ${need[*]}"
  echo "  NO sudo?  Use conda/mamba (user-space, no root needed):"
  echo "      mamba env create -f $HERE/environment.yml   # or: conda env create -f ..."
  echo "      conda activate clair3-connect && ./install.sh ${WHAT}"
  echo "      (already have the env? add the missing bits: mamba install -c conda-forge libprotobuf)"
  echo "  Have sudo?"
  if   command -v apt-get >/dev/null; then
    echo "      sudo apt-get install -y build-essential protobuf-compiler zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev"
  elif command -v dnf >/dev/null; then
    echo "      sudo dnf install -y gcc gcc-c++ protobuf-compiler zlib-devel bzip2-devel xz-devel libcurl-devel"
  else
    echo "      install: gcc, protobuf-compiler, zlib, bzip2, xz/lzma, libcurl (dev headers)"
  fi
  die "install the deps above (system OR conda), then re-run ./install.sh"
fi
[ -n "$IN_CONDA" ] && say "Using conda toolchain at $IN_CONDA"

export PROTOC="${PROTOC:-$(command -v protoc)}"

build_client() {
  say "Building the CLIENT (no GPU / no LibTorch required)…"
  ( cd "$RUST" && cargo build --release -p clair3-connect-client )
  say "CLIENT ready -> $RUST/target/release/clair3-connect"
}

pick_cuda() {

  case "$CUDA_SPEC" in
    cpu|cu118|cu121|cu124) echo "$CUDA_SPEC"; return;;
  esac
  if ! command -v nvidia-smi >/dev/null 2>&1; then echo cpu; return; fi

  local v; v="$(nvidia-smi 2>/dev/null | grep -oE 'CUDA Version: [0-9]+\.[0-9]+' | grep -oE '[0-9]+\.[0-9]+' | head -1)"
  local maj="${v%%.*}" min="${v##*.}"; maj="${maj:-0}"; min="${min:-0}"
  if   [ "$maj" -ge 12 ] && [ "$min" -ge 4 ]; then echo cu124
  elif [ "$maj" -ge 12 ] && [ "$min" -ge 1 ]; then echo cu121
  elif [ "$maj" -ge 12 ];                      then echo cu121
  elif [ "$maj" -ge 11 ] && [ "$min" -ge 8 ]; then echo cu118
  else echo cpu; fi
}

fetch_libtorch() {
  local flavour="$1" dest="$HERE/libtorch"
  if [ -e "$dest/lib/libtorch.so" ]; then say "LibTorch already present at $dest"; echo "$dest"; return; fi
  local url="https://download.pytorch.org/libtorch/${flavour}/libtorch-cxx11-abi-shared-with-deps-${TORCH_VER}%2B${flavour}.zip"
  local zip="$HERE/libtorch-${flavour}.zip"
  say "Downloading LibTorch ${TORCH_VER} (${flavour})… (~2 GB for CUDA builds)"
  command -v curl >/dev/null && curl -L --fail -o "$zip" "$url" || die "download failed: $url"
  say "Unzipping…"; ( cd "$HERE" && unzip -q -o "$zip" ) || die "unzip failed (need 'unzip')"
  rm -f "$zip"
  [ -e "$dest/lib/libtorch.so" ] || die "LibTorch not found after unzip at $dest"
  echo "$dest"
}

build_server() {
  local flavour; flavour="$(pick_cuda)"
  say "Server CUDA flavour: $flavour  (override: ./install.sh server cpu|cu118|cu121|cu124)"
  local LT; LT="$(fetch_libtorch "$flavour")"
  export LIBTORCH="$LT"
  export LD_LIBRARY_PATH="$LT/lib:${LD_LIBRARY_PATH:-}"
  export LIBTORCH_BYPASS_VERSION_CHECK=1

  if [ -x /usr/bin/gcc ]; then
    export CC=/usr/bin/gcc CXX=/usr/bin/g++ CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_LINKER=/usr/bin/gcc
    say "Linking with host gcc (/usr/bin/gcc)."
  elif [ -n "${CONDA_PREFIX:-}" ]; then
    say "WARNING: no /usr/bin/gcc — linking with the conda compiler. If you hit"
    say "  'undefined reference: log@GLIBC_2.29', the host glibc is < 2.29 and the"
    say "  server can't be built/run here; use a host with glibc >= 2.29."
  fi
  say "Building the SERVER…"
  ( cd "$RUST" && cargo build --release -p clair3-connect-server )

  ln -sfn "$LT/lib" "$RUST/target/lib"
  say "SERVER ready -> $RUST/target/release/clair3-connect-server"
}

mkdir -p "$HERE/bin"
case "$WHAT" in
  client) build_client ;;
  server) build_server ;;
  all)    build_client; build_server ;;
  *) die "usage: ./install.sh [client|server|all] [cpu|cu118|cu121|cu124|auto]";;
esac

[ -e "$RUST/target/release/clair3-connect" ]        && ln -sfn "$RUST/target/release/clair3-connect"        "$HERE/bin/clair3-connect"
[ -e "$RUST/target/release/clair3-connect-server" ] && ln -sfn "$RUST/target/release/clair3-connect-server" "$HERE/bin/clair3-connect-server"
echo
say "Done. Binaries are in:  $HERE/bin/"
say "Add to PATH:  export PATH=\"$HERE/bin:\$PATH\""
say "Next: see README.md for run modes. The client needs no env; the server"
say "      finds its LibTorch automatically via a baked rpath."
