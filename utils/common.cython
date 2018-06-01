from distutils.extension import Extension


DEFAULT_EXTENSION_KWARGS = {
    "extra_compile_args": ["-w"]
}


def pyx_extension(**kwargs):
    for key, value in DEFAULT_EXTENSION_KWARGS.items():
        if key not in kwargs:
            kwargs[key] = value

    return Extension(**kwargs)


def make_ext(modname, pyxfilename):
    return pyx_extension(name=modname, sources=[pyxfilename])
