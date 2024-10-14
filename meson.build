project('edipy2', 'fortran')

python = import('python').find_installation(pure: false)

scifor_dep= dependency('scifor', required:true)
edipack_dep= dependency('edipack2', required:true)

fortran_src = ['src/ED_INPUT_VARS.f90',
               'src/python/c_bindings/edi2py.f90']
python_src  = ['src/python/edipy2/edi2py.py',
               'src/python/edipy2/func_read_input.py',
               'src/python/edipy2/func_aux_funx.py',
               'src/python/edipy2/func_bath.py',
               'src/python/edipy2/func_main.py',
               'src/python/edipy2/func_io.py',
               'src/python/edipy2/func_bath_fit.py',
               'src/python/edipy2/__init__.py']

library('edi2py', 
                    fortran_src,
                    fortran_args: ['-ffree-line-length-none', '-cpp', '-D_MPI'],
                    dependencies: [scifor_dep,edipack_dep],
                    install: true,
                    install_dir: python.get_install_dir() / 'edipy2'
)

python.install_sources(
    python_src,
    subdir: 'edipy2'
)
