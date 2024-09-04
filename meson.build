project('edipy2', 'fortran')

python = import('python').find_installation(pure: false)

scifor_dep= dependency('scifor', required:true)
edipack_dep= dependency('edipack2', required:true)

fortran_src = ['src/ED_INPUT_VARS.f90',
               'src/edi2py/edi2py.f90']
python_src  = ['python/edi2py.py',
               'python/func_read_input.py',
               'python/func_aux_funx.py',
               'python/func_bath.py',
               'python/func_main.py',
               'python/func_io.py',
               'python/func_bath_fit.py',
               'python/__init__.py']

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
