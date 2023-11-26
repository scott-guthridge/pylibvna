project('libvna', 'cython',
    version: '0.0.5-7')

py_mod = import('python')
py3 = py_mod.find_installation(pure: false)
dep_py = py3.dependency()
dep_vna = dependency('libvna', version: '>=0.3.9')

srcdir = meson.source_root()
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


test(
    'test_conv',
    py3,
    args: [join_paths(testdir, 'test_conv.py')],
    workdir: builddir
)
test(
    'test_data',
    py3,
    args: [join_paths(testdir, 'test_data.py')],
    workdir: builddir
)

subdir('libvna')
