project('libvna', 'cython',
    version: '0.0.3')

py_mod = import('python')
py3 = py_mod.find_installation(pure: false)
dep_py = py3.dependency()
dep_vna = dependency('libvna', version: '>=0.3.9')

subdir('libvna')
