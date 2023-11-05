project('libvna', 'cython',
    version: '0.0.4')

py_mod = import('python')
py3 = py_mod.find_installation(pure: false)
dep_py = py3.dependency()
dep_vna = dependency('libvna', version: '>=0.3.9')

srcdir = meson.source_root()
testdir = join_paths(srcdir, 'tests')
builddir = meson.current_build_dir()

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
