py3.install_sources(
  ['__init__.py'],
  subdir: 'libvna'
)

py3.extension_module(
    'libvna.cal',
    ['cal.pyx'],
    dependencies: dep_py,
    subdir: 'libvna',
    install: true
)
py3.extension_module(
    'libvna.conv',
    ['conv.pyx'],
    dependencies: dep_py,
    subdir: 'libvna',
    install: true
)
py3.extension_module(
    'libvna.data',
    ['data.pyx'],
    dependencies: dep_py,
    subdir: 'libvna',
    install: true
)

