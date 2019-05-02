make clean
sphinx-apidoc -o .  ../simnibs\
    ../simnibs/GUI/ \
    ../simnibs/cython_code/
make html
firefox _build/html/index.html
rm simnibs*.rst
rm modules.rst
