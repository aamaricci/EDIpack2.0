project('edipy2', 'fortran')

python = import('python').find_installation(pure: false)

scifor_dep= dependency('scifor', required:true)
edipack_dep= dependency('edipack2', required:true)

fortran_src = ['src/ED_INPUT_VARS.f90',
               'src/edi2py/edi2py.f90']
python_src  = ['edipy2/edi2py.py',
               'edipy2/func_read_input.py',
               'edipy2/func_aux_funx.py',
               'edipy2/func_bath.py',
               'edipy2/func_main.py',
               'edipy2/func_io.py',
               'edipy2/func_bath_fit.py',
               'edipy2/__init__.py']

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
