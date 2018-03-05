from cffi import FFI
ffi = FFI()


ffi.set_source("_axsym", open("axsym.cpp").read(), source_extension='.cpp')
ffi.cdef(open("axsym.h").read().replace('extern "C" ',""))

if __name__ == "__main__":
    ffi.compile()