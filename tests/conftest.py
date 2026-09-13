"""Execute each restored chapter once for numerical checks (browser tested separately)."""
from pathlib import Path
import os,runpy,sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pytest
ROOT=Path(__file__).resolve().parents[1]

@pytest.fixture(scope='session')
def restored_chapter():
    cache={}
    def load(name):
        if name not in cache:
            old=Path.cwd();os.chdir(ROOT/'book');sys.path.insert(0,str(ROOT/'book'))
            import _support
            original_lab,original_diagram=_support.lab,_support.diagram
            _support.lab=lambda *args,**kwargs:None
            _support.diagram=lambda *args,**kwargs:None
            try:cache[name]=runpy.run_path(str(ROOT/'book'/name))
            finally:
                _support.lab,_support.diagram=original_lab,original_diagram
                plt.close('all');os.chdir(old)
        return cache[name]
    return load
