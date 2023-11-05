py3.install_sources(
  ['__init__.py'],
  subdir: 'libvna'
)

py3.extension_module(
    'cal',
    ['cal.pyx'],
    dependencies: [dep_py, dep_vna],
    subdir: 'libvna',
    install: true
)
py3.extension_module(
    'conv',
    ['conv.pyx'],
    dependencies: [dep_py, dep_vna],
    subdir: 'libvna',
    install: true
)
py3.extension_module(
    'data',
    ['data.pyx'],
    dependencies: [dep_py, dep_vna],
    subdir: 'libvna',
    install: true
)

