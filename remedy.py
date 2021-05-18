from os import listdir, chdir, path, rename
from shutil import rmtree

items = listdir('.')

for item in items:
    if path.isdir(item):
        chdir(item)
        if not path.isfile('coord.png'):
            rename(f'{item}.mol2', f'../{item}.mol2')
            chdir('../')
            rmtree(item)
            print(f'{item} restored!')
        else: chdir('../')