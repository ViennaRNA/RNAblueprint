# This file was automatically generated by SWIG (http://www.swig.org).
# Version 3.0.7
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.





from sys import version_info
if version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_RNAdesign', [dirname(__file__)])
        except ImportError:
            import _RNAdesign
            return _RNAdesign
        if fp is not None:
            try:
                _mod = imp.load_module('_RNAdesign', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _RNAdesign = swig_import_helper()
    del swig_import_helper
else:
    import _RNAdesign
del version_info
try:
    _swig_property = property
except NameError:
    pass  # Python < 2.2 doesn't have 'property'.


def _swig_setattr_nondynamic(self, class_type, name, value, static=1):
    if (name == "thisown"):
        return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name, None)
    if method:
        return method(self, value)
    if (not static):
        if _newclass:
            object.__setattr__(self, name, value)
        else:
            self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)


def _swig_setattr(self, class_type, name, value):
    return _swig_setattr_nondynamic(self, class_type, name, value, 0)


def _swig_getattr_nondynamic(self, class_type, name, static=1):
    if (name == "thisown"):
        return self.this.own()
    method = class_type.__swig_getmethods__.get(name, None)
    if method:
        return method(self)
    if (not static):
        return object.__getattr__(self, name)
    else:
        raise AttributeError(name)

def _swig_getattr(self, class_type, name):
    return _swig_getattr_nondynamic(self, class_type, name, 0)


def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object:
        pass
    _newclass = 0


class SwigPyIterator(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, SwigPyIterator, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, SwigPyIterator, name)

    def __init__(self, *args, **kwargs):
        raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_destroy__ = _RNAdesign.delete_SwigPyIterator
    __del__ = lambda self: None

    def value(self):
        return _RNAdesign.SwigPyIterator_value(self)

    def incr(self, n=1):
        return _RNAdesign.SwigPyIterator_incr(self, n)

    def decr(self, n=1):
        return _RNAdesign.SwigPyIterator_decr(self, n)

    def distance(self, x):
        return _RNAdesign.SwigPyIterator_distance(self, x)

    def equal(self, x):
        return _RNAdesign.SwigPyIterator_equal(self, x)

    def copy(self):
        return _RNAdesign.SwigPyIterator_copy(self)

    def next(self):
        return _RNAdesign.SwigPyIterator_next(self)

    def __next__(self):
        return _RNAdesign.SwigPyIterator___next__(self)

    def previous(self):
        return _RNAdesign.SwigPyIterator_previous(self)

    def advance(self, n):
        return _RNAdesign.SwigPyIterator_advance(self, n)

    def __eq__(self, x):
        return _RNAdesign.SwigPyIterator___eq__(self, x)

    def __ne__(self, x):
        return _RNAdesign.SwigPyIterator___ne__(self, x)

    def __iadd__(self, n):
        return _RNAdesign.SwigPyIterator___iadd__(self, n)

    def __isub__(self, n):
        return _RNAdesign.SwigPyIterator___isub__(self, n)

    def __add__(self, n):
        return _RNAdesign.SwigPyIterator___add__(self, n)

    def __sub__(self, *args):
        return _RNAdesign.SwigPyIterator___sub__(self, *args)
    def __iter__(self):
        return self
SwigPyIterator_swigregister = _RNAdesign.SwigPyIterator_swigregister
SwigPyIterator_swigregister(SwigPyIterator)

class StringVector(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, StringVector, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, StringVector, name)
    __repr__ = _swig_repr

    def iterator(self):
        return _RNAdesign.StringVector_iterator(self)
    def __iter__(self):
        return self.iterator()

    def __nonzero__(self):
        return _RNAdesign.StringVector___nonzero__(self)

    def __bool__(self):
        return _RNAdesign.StringVector___bool__(self)

    def __len__(self):
        return _RNAdesign.StringVector___len__(self)

    def pop(self):
        return _RNAdesign.StringVector_pop(self)

    def __getslice__(self, i, j):
        return _RNAdesign.StringVector___getslice__(self, i, j)

    def __setslice__(self, *args):
        return _RNAdesign.StringVector___setslice__(self, *args)

    def __delslice__(self, i, j):
        return _RNAdesign.StringVector___delslice__(self, i, j)

    def __delitem__(self, *args):
        return _RNAdesign.StringVector___delitem__(self, *args)

    def __getitem__(self, *args):
        return _RNAdesign.StringVector___getitem__(self, *args)

    def __setitem__(self, *args):
        return _RNAdesign.StringVector___setitem__(self, *args)

    def append(self, x):
        return _RNAdesign.StringVector_append(self, x)

    def empty(self):
        return _RNAdesign.StringVector_empty(self)

    def size(self):
        return _RNAdesign.StringVector_size(self)

    def clear(self):
        return _RNAdesign.StringVector_clear(self)

    def swap(self, v):
        return _RNAdesign.StringVector_swap(self, v)

    def get_allocator(self):
        return _RNAdesign.StringVector_get_allocator(self)

    def begin(self):
        return _RNAdesign.StringVector_begin(self)

    def end(self):
        return _RNAdesign.StringVector_end(self)

    def rbegin(self):
        return _RNAdesign.StringVector_rbegin(self)

    def rend(self):
        return _RNAdesign.StringVector_rend(self)

    def pop_back(self):
        return _RNAdesign.StringVector_pop_back(self)

    def erase(self, *args):
        return _RNAdesign.StringVector_erase(self, *args)

    def __init__(self, *args):
        this = _RNAdesign.new_StringVector(*args)
        try:
            self.this.append(this)
        except:
            self.this = this

    def push_back(self, x):
        return _RNAdesign.StringVector_push_back(self, x)

    def front(self):
        return _RNAdesign.StringVector_front(self)

    def back(self):
        return _RNAdesign.StringVector_back(self)

    def assign(self, n, x):
        return _RNAdesign.StringVector_assign(self, n, x)

    def resize(self, *args):
        return _RNAdesign.StringVector_resize(self, *args)

    def insert(self, *args):
        return _RNAdesign.StringVector_insert(self, *args)

    def reserve(self, n):
        return _RNAdesign.StringVector_reserve(self, n)

    def capacity(self):
        return _RNAdesign.StringVector_capacity(self)
    __swig_destroy__ = _RNAdesign.delete_StringVector
    __del__ = lambda self: None
StringVector_swigregister = _RNAdesign.StringVector_swigregister
StringVector_swigregister(StringVector)

class IntVector(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, IntVector, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, IntVector, name)
    __repr__ = _swig_repr

    def iterator(self):
        return _RNAdesign.IntVector_iterator(self)
    def __iter__(self):
        return self.iterator()

    def __nonzero__(self):
        return _RNAdesign.IntVector___nonzero__(self)

    def __bool__(self):
        return _RNAdesign.IntVector___bool__(self)

    def __len__(self):
        return _RNAdesign.IntVector___len__(self)

    def pop(self):
        return _RNAdesign.IntVector_pop(self)

    def __getslice__(self, i, j):
        return _RNAdesign.IntVector___getslice__(self, i, j)

    def __setslice__(self, *args):
        return _RNAdesign.IntVector___setslice__(self, *args)

    def __delslice__(self, i, j):
        return _RNAdesign.IntVector___delslice__(self, i, j)

    def __delitem__(self, *args):
        return _RNAdesign.IntVector___delitem__(self, *args)

    def __getitem__(self, *args):
        return _RNAdesign.IntVector___getitem__(self, *args)

    def __setitem__(self, *args):
        return _RNAdesign.IntVector___setitem__(self, *args)

    def append(self, x):
        return _RNAdesign.IntVector_append(self, x)

    def empty(self):
        return _RNAdesign.IntVector_empty(self)

    def size(self):
        return _RNAdesign.IntVector_size(self)

    def clear(self):
        return _RNAdesign.IntVector_clear(self)

    def swap(self, v):
        return _RNAdesign.IntVector_swap(self, v)

    def get_allocator(self):
        return _RNAdesign.IntVector_get_allocator(self)

    def begin(self):
        return _RNAdesign.IntVector_begin(self)

    def end(self):
        return _RNAdesign.IntVector_end(self)

    def rbegin(self):
        return _RNAdesign.IntVector_rbegin(self)

    def rend(self):
        return _RNAdesign.IntVector_rend(self)

    def pop_back(self):
        return _RNAdesign.IntVector_pop_back(self)

    def erase(self, *args):
        return _RNAdesign.IntVector_erase(self, *args)

    def __init__(self, *args):
        this = _RNAdesign.new_IntVector(*args)
        try:
            self.this.append(this)
        except:
            self.this = this

    def push_back(self, x):
        return _RNAdesign.IntVector_push_back(self, x)

    def front(self):
        return _RNAdesign.IntVector_front(self)

    def back(self):
        return _RNAdesign.IntVector_back(self)

    def assign(self, n, x):
        return _RNAdesign.IntVector_assign(self, n, x)

    def resize(self, *args):
        return _RNAdesign.IntVector_resize(self, *args)

    def insert(self, *args):
        return _RNAdesign.IntVector_insert(self, *args)

    def reserve(self, n):
        return _RNAdesign.IntVector_reserve(self, n)

    def capacity(self):
        return _RNAdesign.IntVector_capacity(self)
    __swig_destroy__ = _RNAdesign.delete_IntVector
    __del__ = lambda self: None
IntVector_swigregister = _RNAdesign.IntVector_swigregister
IntVector_swigregister(IntVector)


def initialize_library(debug):
    return _RNAdesign.initialize_library(debug)
initialize_library = _RNAdesign.initialize_library
class DependencyGraphMT(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, DependencyGraphMT, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, DependencyGraphMT, name)
    __repr__ = _swig_repr

    def __init__(self, *args):
        this = _RNAdesign.new_DependencyGraphMT(*args)
        try:
            self.this.append(this)
        except:
            self.this = this
    __swig_destroy__ = _RNAdesign.delete_DependencyGraphMT
    __del__ = lambda self: None

    def get_sequence(self):
        return _RNAdesign.DependencyGraphMT_get_sequence(self)

    def set_sequence(self, *args):
        return _RNAdesign.DependencyGraphMT_set_sequence(self, *args)

    def mutate_local(self, *args):
        return _RNAdesign.DependencyGraphMT_mutate_local(self, *args)

    def mutate_global(self, *args):
        return _RNAdesign.DependencyGraphMT_mutate_global(self, *args)

    def mutate(self, *args):
        return _RNAdesign.DependencyGraphMT_mutate(self, *args)

    def number_of_sequences(self, *args):
        return _RNAdesign.DependencyGraphMT_number_of_sequences(self, *args)

    def number_of_connected_components(self):
        return _RNAdesign.DependencyGraphMT_number_of_connected_components(self)

    def component_vertices(self, connected_component_ID):
        return _RNAdesign.DependencyGraphMT_component_vertices(self, connected_component_ID)

    def special_vertices(self):
        return _RNAdesign.DependencyGraphMT_special_vertices(self)
DependencyGraphMT_swigregister = _RNAdesign.DependencyGraphMT_swigregister
DependencyGraphMT_swigregister(DependencyGraphMT)

# This file is compatible with both classic and new-style classes.


