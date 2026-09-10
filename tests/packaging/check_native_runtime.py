"""Check a Linux shared CUDA/cuQuantum install without modifying it.

Example (after building the ordinary packaging/consumer fixture):
  python3 tests/packaging/check_native_runtime.py --build-dir build \
    --library /opt/QuEST/lib/libQuEST.so --consumer consumer-build/consumer_c \
    --consumer consumer-build/consumer_cpp
"""

import argparse
import os
from pathlib import Path
import re
import subprocess


def run(*command):
    # Also remove less common loader controls, including LD_RUN_PATH and LD_AUDIT.
    environment = {key: value for key, value in os.environ.items()
                   if not key.startswith(("LD_", "DYLD_"))
                   and key not in ("LIBPATH", "SHLIB_PATH")}
    result = subprocess.run(command, env=environment, text=True, check=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print(result.stdout, end="", flush=True)
    return result.stdout


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def dynamic(path):
    print(f"Checking ELF dynamic section: {path}", flush=True)
    output = run("readelf", "-d", str(path))
    needed = re.findall(r"\(NEEDED\).*\[([^]]+)\]", output)
    sonames = re.findall(r"\(SONAME\).*\[([^]]+)\]", output)
    search = re.findall(r"\((?:RPATH|RUNPATH)\).*\[([^]]*)\]", output)
    paths = [entry for value in search for entry in value.split(":")]
    for entry in paths:
        require(entry, f"Empty loader search directory in {path}")
        require(not is_stub(Path(entry)), f"Stub loader search directory in {path}: {entry}")
    return needed, sonames, paths


def is_stub(path):
    return any(part.lower() in ("stub", "stubs")
               for part in (*path.parts, *path.resolve().parts))


def loaded(path, expected):
    print(f"Checking loader resolution with loader variables unset: {path}", flush=True)
    output = run("ldd", str(path))
    require("not found" not in output, f"Unresolved runtime dependency in {path}")
    resolved = dict(re.findall(r"^\s*(\S+)\s+=>\s+(/.*?)\s+\(", output, re.MULTILINE))
    for soname, location in resolved.items():
        require(not is_stub(Path(location)), f"Loader resolved {soname} to a stub: {location}")
    for soname, library in expected.items():
        require(soname in resolved, f"{path} did not load {soname}")
        require(Path(resolved[soname]).resolve() == library,
                f"{path} loaded {soname} from {resolved[soname]}, expected {library}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--library", type=Path, required=True)
    parser.add_argument("--consumer", type=Path, action="append", default=[])
    parser.add_argument("--require-shared-cudart", action="store_true")
    args = parser.parse_args()
    cache = dict(re.findall(r"^([^/#:\r\n][^:\r\n]*):[^=\r\n]+=(.*)$",
                            (args.build_dir / "CMakeCache.txt").read_text(), re.MULTILINE))
    library = args.library.resolve(strict=True)
    needed, sonames, search = dynamic(library)
    require(len(sonames) == 1, f"Expected one shared QuEST SONAME in {library}")
    search_dirs = {Path(entry).resolve() for entry in search if Path(entry).is_absolute()}
    expected = {}
    # CUDA's default runtime can be static. Verify shared cudart when requested,
    # and always verify the real cuStateVec and CUDA directories in QuEST's ELF.
    for key, required in (("CUQUANTUM_cuStateVec_LIBRARY", True),
                          ("CUDA_cudart_LIBRARY", args.require_shared_cudart),
                          ("CUDA_cublas_LIBRARY", False),
                          ("CUDA_cublasLt_LIBRARY", False)):
        require(key in cache, f"Missing SDK discovery result {key}")
        dependency = Path(cache[key]).resolve(strict=True)
        require(not is_stub(dependency), f"SDK discovery selected a stub: {dependency}")
        require(dependency.parent in search_dirs,
                f"QuEST RPATH/RUNPATH lacks real SDK directory {dependency.parent}")
        _, dependency_sonames, _ = dynamic(dependency)
        require(len(dependency_sonames) == 1, f"Expected shared SDK library: {dependency}")
        soname = dependency_sonames[0]
        if required:
            require(soname in needed, f"QuEST has no NEEDED entry for {soname}")
        if soname in needed:
            expected[soname] = dependency
    loaded(library, expected)
    for consumer in args.consumer:
        consumer = consumer.resolve(strict=True)
        consumer_needed, _, _ = dynamic(consumer)
        require(sonames[0] in consumer_needed, f"{consumer} does not link shared QuEST")
        loaded(consumer, {**expected, sonames[0]: library})
        run(str(consumer))
    print("Native CUDA/cuQuantum runtime check passed.", flush=True)


if __name__ == "__main__":
    main()
