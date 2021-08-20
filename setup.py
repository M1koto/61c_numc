from distutils.core import setup, Extension
import sysconfig


def main():
    CFLAGS = ['-g', '-Wall', '-std=c99', '-fopenmp', '-mavx', '-mfma', '-pthread', '-O3']
    LDFLAGS = ['-fopenmp']
    # Use the setup function we imported and set up the modules.
    # You may find this reference helpful: https://docs.python.org/3.6/extending/building.html
    module1 = Extension('numc',
                        sources=['numc.c', 'matrix.c'],
                        extra_link_args=LDFLAGS,
                        extra_compile_args=CFLAGS,)
    setup(name='numc',
          version='1.0',
          description='CS61c Project 4 numc module',
          ext_modules=[module1])

    # alternatively...
    #sysconfig.get_config_vars()['CFLAGS'] = CFLAGS
    #sysconfig.get_config_var()['LDFLAGS'] = LDFLAGS


if __name__ == "__main__":
    main()
