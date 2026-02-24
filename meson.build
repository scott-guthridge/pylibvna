project('libvna', 'c', 'cython',
  version: '0.2.2',
  license: 'GPL-3',
  default_options: [
    'buildtype=debugoptimized',
    'b_ndebug=if-release',
    'c_std=c99',
  ],
)

py_mod = import('python')
py3 = py_mod.find_installation(pure: false)
dep_py = py3.dependency()
dep_vna = dependency('libvna', version: '>=0.3.9')

srcdir = meson.project_source_root()
testdir = join_paths(srcdir, 'tests')
builddir = meson.current_build_dir()

incdir_numpy = meson.get_external_property('numpy-include-dir', 'none')
if incdir_numpy == 'none'
  incdir_numpy = run_command(py3, ['-c',
'''
import numpy as np
import os

try:
  incdir = os.path.relpath(np.get_include())
except Exception:
  incdir = np.get_include()
print(incdir)
'''],
    check: true).stdout().strip()
endif
inc_np = include_directories(incdir_numpy)
dep_np = declare_dependency(include_directories: inc_np)


# Create test environment that includes build directory in Python path
test_env = environment()
test_env.prepend('PYTHONPATH', meson.current_build_dir())

test(
    'test_cal',
    py3,
    args: [join_paths(testdir, 'test_cal.py')],
    workdir: builddir,
    env: test_env
)
test(
    'test_conv',
    py3,
    args: [join_paths(testdir, 'test_conv.py')],
    workdir: builddir,
    env: test_env
)
test(
    'test_data',
    py3,
    args: [join_paths(testdir, 'test_data.py')],
    workdir: builddir,
    env: test_env
)

subdir('libvna')
