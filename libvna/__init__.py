import importlib

__all__ = [ "cal", "conv", "data" ]

# If libvna only is imported, lazily load the submodules as they're
# accessed.
def __getattr__(name):
    if name in __all__:
        module = importlib.import_module(f"{__name__}.{name}")
        globals()[name] = module
        return module
    raise AttributeError(f"module {__name__} has no attribute {name}")
