"""Exercise QuEST's packaging policy without its library or external dependencies.
Run: python3 tests/packaging/cpack/test_cpack.py
"""
import pathlib
import subprocess
import shutil
import sys
import tempfile
import unittest

REPO = pathlib.Path(__file__).resolve().parents[3]


class Packaging(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="quest cpack ")
        self.addCleanup(self.tmp.cleanup)
        self.source = pathlib.Path(self.tmp.name) / "source"
        self.build = pathlib.Path(self.tmp.name) / "build"
        self.source.mkdir()
        (self.source / "payload").write_text("QuEST fixture\n")
        (self.source / "AUTHORS.txt").write_text("Contact: fixture@example.com\n")
        (self.source / "LICENCE.txt").write_text("MIT License\n")
        (self.source / "CMakeLists.txt").write_text(f'''
cmake_minimum_required(VERSION 3.28)
project(QuEST VERSION 4.3.0 LANGUAGES CXX)
include(GNUInstallDirs)
set(QUEST_ENABLE_INSTALL ON)
set(QUEST_ENABLE_PACKAGING ON)
set(QUEST_FLOAT_PRECISION 2)
option(QUEST_BUILT_SHARED "" ON)
install(FILES payload DESTINATION include COMPONENT Development)
if(QUEST_BUILT_SHARED)
 install(FILES payload DESTINATION lib COMPONENT Runtime)
endif()
install(FILES payload DESTINATION foreign COMPONENT ForeignDependency)
include("{REPO.as_posix()}/cmake/QuESTPackaging.cmake")
''')

    def run_command(self, *args, success=True):
        result = subprocess.run(args, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if success:
            self.assertEqual(result.returncode, 0, result.stdout)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout)
        return result.stdout

    def configure(self, *args):
        self.run_command("cmake", "-S", str(self.source), "-B", str(self.build),
                         "-DCMAKE_BUILD_TYPE=Release", *args)

    def policy(self, generator, *settings, success=True):
        script = self.build / "inspect.cmake"
        script.write_text(f'include("{self.build.as_posix()}/CPackConfig.cmake")\n'
                          f'set(CPACK_GENERATOR {generator})\n' + "\n".join(settings) + '''
include("${CPACK_PROJECT_CONFIG_FILE}")
file(WRITE "${CMAKE_CURRENT_LIST_DIR}/policy.txt" "${CPACK_COMPONENTS_ALL}\n${CPACK_DEBIAN_DEVELOPMENT_PACKAGE_DEPENDS}\n${CPACK_RPM_DEVELOPMENT_PACKAGE_REQUIRES}\n${CPACK_PACKAGE_FILE_NAME}\n")
''')
        return self.run_command("cmake", "-P", str(script), success=success)

    def test_complete_archives_only_contain_quest_components(self):
        # A parent project must not turn dependency installation into bundled content.
        self.configure("-DCPACK_MONOLITHIC_INSTALL=ON")
        for generator, extension in [("TGZ", "tar.gz"), ("ZIP", "zip")]:
            self.run_command("cpack", "--config", str(self.build / "CPackConfig.cmake"),
                             "-G", generator, "-B", str(self.build / "packages"))
            archive, = (self.build / "packages").glob(f"*.{extension}")
            self.assertIn("Release-shared-fp2-cpu", archive.name)
            listing = self.run_command("cmake", "-E", "tar", "tf", str(archive))
            archive_root = archive.name[:-(len(extension) + 1)]
            self.assertIn(f"{archive_root}/include/payload", listing)
            self.assertIn(f"{archive_root}/lib/payload", listing)
            self.assertNotIn("foreign", listing)

    @unittest.skipUnless(shutil.which("dpkg-deb") and shutil.which("dpkg-shlibdeps"),
                         "Debian tools are needed for native component scanning")
    def test_deb_scanner_resolves_library_in_sibling_runtime_component(self):
        (self.source / "runtime.cpp").write_text("int runtime_function() { return 0; }\n")
        (self.source / "example.cpp").write_text(
            "extern int runtime_function(); int main() { return runtime_function(); }\n")
        cmake_file = self.source / "CMakeLists.txt"
        content = cmake_file.read_text().replace(f'include("{REPO.as_posix()}/cmake/QuESTPackaging.cmake")',
            f'''add_library(runtime SHARED runtime.cpp)
set_target_properties(runtime PROPERTIES SOVERSION 4)
add_executable(example example.cpp)
target_link_libraries(example PRIVATE runtime)
set_target_properties(example PROPERTIES INSTALL_RPATH "$ORIGIN/../lib")
install(TARGETS runtime LIBRARY DESTINATION lib COMPONENT Runtime NAMELINK_COMPONENT Development)
install(TARGETS example RUNTIME DESTINATION bin COMPONENT Examples)
set(QUEST_HAVE_INSTALLABLE_EXAMPLES ON)
include("{REPO.as_posix()}/cmake/QuESTPackaging.cmake")''')
        cmake_file.write_text(content)
        self.configure("-DCMAKE_INSTALL_PREFIX=/usr", "-DCMAKE_INSTALL_LIBDIR=lib",
                       "-DQUEST_NATIVE_PACKAGE_PROFILE=ubuntu24.04")
        self.run_command("cmake", "--build", str(self.build), "--parallel", "2")
        self.run_command("cpack", "--config", str(self.build / "CPackConfig.cmake"),
                         "-G", "DEB", "-B", str(self.build / "packages"))
        examples, = (self.build / "packages").glob("quest-examples_*.deb")
        metadata = self.run_command("dpkg-deb", "--field", str(examples), "Depends")
        self.assertIn("libquest4 (= 4.3.0-1)", metadata)

    def test_explicit_archive_prefix_is_preserved(self):
        self.configure("-DCPACK_PACKAGING_INSTALL_PREFIX=/custom-prefix")
        self.run_command("cpack", "--config", str(self.build / "CPackConfig.cmake"),
                         "-G", "TGZ", "-B", str(self.build / "packages"))
        archive, = (self.build / "packages").glob("*.tar.gz")
        listing = self.run_command("cmake", "-E", "tar", "tf", str(archive))
        self.assertIn("/custom-prefix/include/payload", listing)

    @unittest.skipUnless(sys.platform.startswith("linux"), "Native package policy requires Linux")
    def test_native_shared_and_static_dependencies(self):
        for shared in ["ON", "OFF"]:
            self.configure(f"-DQUEST_BUILT_SHARED={shared}", "-DQUEST_NATIVE_PACKAGE_PROFILE=ubuntu24.04",
                           "-DCMAKE_INSTALL_PREFIX=/usr", "-DQUEST_ENABLE_OMP=ON", "-DQUEST_ENABLE_NUMA=ON")
            self.policy("DEB")
            content = (self.build / "policy.txt").read_text()
            self.assertIn("g++", content)
            self.assertEqual("libquest4 (= 4.3.0-1)" in content, shared == "ON")
            self.assertEqual("libnuma-dev" in content, shared == "OFF")
            self.policy("RPM", 'set(CPACK_QUEST_NATIVE_PROFILE fedora44)')
            content = (self.build / "policy.txt").read_text()
            self.assertEqual("quest = 4.3.0-1" in content, shared == "ON")

    @unittest.skipUnless(sys.platform.startswith("linux"), "Native package policy requires Linux")
    def test_native_validation_is_deferred_and_vendor_metadata_required(self):
        self.configure("-DQUEST_ENABLE_CUDA=ON", "-DQUEST_NATIVE_PACKAGE_PROFILE=ubuntu24.04",
                       "-DCMAKE_INSTALL_PREFIX=/usr")
        self.policy("TGZ")
        error = self.policy("DEB", success=False)
        self.assertIn("CPACK_DEBIAN_RUNTIME_PACKAGE_DEPENDS", error)
        self.policy("DEB", 'set(CPACK_DEBIAN_DEVELOPMENT_PACKAGE_DEPENDS "g++, cuda-toolkit")',
                    'set(CPACK_DEBIAN_RUNTIME_PACKAGE_DEPENDS "cuda-cudart")')

    @unittest.skipUnless(sys.platform.startswith("linux"), "Native package policy requires Linux")
    def test_overrides_and_unknown_profiles(self):
        self.configure("-DCPACK_PACKAGE_FILE_NAME=custom", "-DQUEST_NATIVE_PACKAGE_PROFILE=unknown",
                       "-DCMAKE_INSTALL_PREFIX=/usr")
        self.policy("TGZ")
        self.assertIn("custom", (self.build / "policy.txt").read_text())
        self.assertIn("QUEST_NATIVE_PACKAGE_PROFILE", self.policy("DEB", success=False))

    def test_cpack_time_filename_override_and_configuration(self):
        self.configure()
        self.policy("TGZ", 'set(CPACK_BUILD_CONFIG Debug)')
        self.assertIn("Debug-shared", (self.build / "policy.txt").read_text())
        self.policy("TGZ", 'set(CPACK_PACKAGE_FILE_NAME explicit-at-package-time)')
        self.assertIn("explicit-at-package-time", (self.build / "policy.txt").read_text())

    @unittest.skipUnless(sys.platform.startswith("linux"), "Native package policy requires Linux")
    def test_exact_runtime_constraint_respects_native_version_overrides(self):
        self.configure("-DCMAKE_INSTALL_PREFIX=/usr", "-DQUEST_NATIVE_PACKAGE_PROFILE=ubuntu24.04")
        self.policy("DEB", 'set(CPACK_DEBIAN_PACKAGE_VERSION 4.3.1)',
                    'set(CPACK_DEBIAN_PACKAGE_RELEASE 7)', 'set(CPACK_DEBIAN_PACKAGE_EPOCH 2)',
                    'set(CPACK_DEBIAN_RUNTIME_PACKAGE_NAME alternate-runtime)',
                    'set(CPACK_DEBIAN_DEVELOPMENT_PACKAGE_DEPENDS custom-dependency)')
        self.assertIn("custom-dependency, alternate-runtime (= 2:4.3.1-7)",
                      (self.build / "policy.txt").read_text())
        self.policy("RPM", 'set(CPACK_QUEST_NATIVE_PROFILE fedora44)',
                    'set(CPACK_RPM_PACKAGE_VERSION 4.3.1)', 'set(CPACK_RPM_PACKAGE_RELEASE 7)',
                    'set(CPACK_RPM_PACKAGE_EPOCH 2)')
        self.assertIn("quest = 2:4.3.1-7", (self.build / "policy.txt").read_text())

    @unittest.skipUnless(sys.platform.startswith("linux"), "Native package policy requires Linux")
    def test_native_rejects_wrong_prefix_and_requires_non_gcc_metadata(self):
        self.configure("-DQUEST_NATIVE_PACKAGE_PROFILE=ubuntu24.04")
        self.assertIn("CMAKE_INSTALL_PREFIX=/usr", self.policy("DEB", success=False))
        self.configure("-DCMAKE_INSTALL_PREFIX=/usr")
        error = self.policy("DEB", 'set(CPACK_QUEST_COMPILER_ID Clang)', success=False)
        self.assertIn("requires explicit", error)

    @unittest.skipUnless(sys.platform.startswith("linux") and
                         (shutil.which("rpm") or shutil.which("dpkg-query")),
                         "Native package ownership tools required")
    def test_stock_profiles_reject_unowned_mpi_artifacts(self):
        generator = "RPM" if shutil.which("rpm") else "DEB"
        profile = "fedora44" if generator == "RPM" else "ubuntu24.04"
        self.configure("-DCMAKE_INSTALL_PREFIX=/usr", "-DQUEST_ENABLE_MPI=ON",
                       f"-DQUEST_NATIVE_PACKAGE_PROFILE={profile}",
                       "-DMPI_CXX_LIBRARIES=/unowned-sdk/libmpich.so")
        error = self.policy(generator, success=False)
        self.assertIn("distro OpenMPI", error)
        self.assertIn("QUEST_NATIVE_PACKAGE_PROFILE=custom", error)
        prefix = "CPACK_RPM" if generator == "RPM" else "CPACK_DEBIAN"
        suffix = "PACKAGE_REQUIRES" if generator == "RPM" else "PACKAGE_DEPENDS"
        self.policy(generator, 'set(CPACK_QUEST_NATIVE_PROFILE custom)',
                    f'set({prefix}_RUNTIME_{suffix} custom-mpi-runtime)',
                    f'set({prefix}_DEVELOPMENT_{suffix} custom-mpi-development)')

    def test_explicit_subproject_packaging_uses_quest_source_and_install_tree(self):
        parent = pathlib.Path(self.tmp.name) / "parent"
        parent.mkdir()
        (parent / "parent-only").write_text("parent payload")
        (parent / "CMakeLists.txt").write_text(f'''cmake_minimum_required(VERSION 3.28)
project(Parent LANGUAGES CXX)
install(FILES parent-only DESTINATION parent COMPONENT Development)
add_subdirectory("{self.source.as_posix()}" quest)
''')
        self.run_command("cmake", "-S", str(parent), "-B", str(self.build))
        for config, directory in [("CPackConfig.cmake", "binary"),
                                  ("CPackSourceConfig.cmake", "source-archive")]:
            self.run_command("cpack", "--config", str(self.build / "quest" / config), "-G", "TGZ",
                             "-B", str(self.build / directory))
            archive, = (self.build / directory).glob("*.tar.gz")
            listing = self.run_command("cmake", "-E", "tar", "tf", str(archive))
            self.assertIn("payload", listing)
            self.assertNotIn("parent-only", listing)

    def test_source_archive_survives_build_named_checkout_parent(self):
        parent = self.source.parent / "build-source-parent"
        parent.mkdir()
        self.source = self.source.rename(parent / "QuEST")
        self.configure()
        self.run_command("cpack", "--config", str(self.build / "CPackSourceConfig.cmake"),
                         "-G", "TGZ", "-B", str(self.build / "source-packages"))
        archive, = (self.build / "source-packages").glob("*.tar.gz")
        listing = self.run_command("cmake", "-E", "tar", "tf", str(archive))
        self.assertIn("QuEST-4.3.0-Source/CMakeLists.txt", listing)
        self.assertIn("QuEST-4.3.0-Source/payload", listing)

    def test_source_archives_exclude_arbitrary_build_trees(self):
        build_dir = self.source / "strangely named compilation"
        build_dir.mkdir()
        (build_dir / "CMakeCache.txt").write_text("cache")
        (build_dir / "junk").write_text("do not distribute")
        (self.source / "CMakeUserPresets.json").write_text("{}")
        self.configure()
        self.run_command("cpack", "--config", str(self.build / "CPackSourceConfig.cmake"),
                         "-G", "TGZ;ZIP", "-B", str(self.build / "source-packages"))
        for archive in (self.build / "source-packages").glob("QuEST-4.3.0-Source.*"):
            listing = self.run_command("cmake", "-E", "tar", "tf", str(archive))
            self.assertIn("QuEST-4.3.0-Source/CMakeLists.txt", listing)
            self.assertNotIn("strangely named", listing)
            self.assertNotIn("CMakeUserPresets", listing)
        self.assertEqual(len(list((self.build / "source-packages").glob("QuEST-4.3.0-Source.*"))), 2)


if __name__ == "__main__":
    unittest.main()
