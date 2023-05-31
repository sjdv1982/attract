#!/usr/bin/env python

from distutils.core import setup
setup(name='spydersilk',
  version='0.03',
  description='Framework and universal editor to define, validate, convert and edit parametric data',
  
  author='Sjoerd de Vries',
  author_email='sjdv1982@gmail.com',
  url='http://www.spyderware.nl/',
  packages=['spyder', 'spyder.gtkform','spyder.qtform','spyder.formtools',
            'spyder.stalkwalk', 'spyder.stalkwalk.markup', 'spyder.stalkwalk.hamlet', 'spyder.stalkwalk.html', 'spyder.stalkwalk.webdemo'],
  package_data = {'spyder': 
    ['WISDOM',
    'tempdir/*.*', 
    'modules/*.*', 
    'modules/*/*.*', 
    'modules/*/*/*.*', 
    'modules/*/*/*/*.*',
    'modules/atom/tools/src/Makefile']
  },
  scripts = ['spyder/silk/silk'],
  long_description = """Spydersilk is a Python framework for parametric data modelling. Spyder models are well-defined hierarchical data structures that can be easily constructed and manipulated in code. Spyder supports powerful data validation and conversion that can be controlled with simple Python rules. Silk is a universal editor for Spyder models, with Qt, HTML and IPython integration.
""",
  classifiers =[
    'License :: OSI Approved :: BSD License',
    'Programming Language :: Python',
    'Development Status :: 3 - Alpha',
    'Operating System :: OS Independent',
    'Topic :: Software Development :: Code Generators',
    'Topic :: Scientific/Engineering :: Visualization',
    'Topic :: Internet :: WWW/HTTP',
    'Topic :: Software Development :: User Interfaces',
    'Topic :: Software Development :: Testing',
    'Topic :: Software Development :: Widget Sets',
  ],
)   
