from setuptools import setup
import os
import re


def read_version():
    path = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'dimred/__init__.py')
    with open(path, 'r') as fh:
        return re.search(r'__version__\s?=\s?[\'"](.+)[\'"]', fh.read()).group(1)


def readme():
    with open('README.md') as f:
        return f.read()


setup(name='dimred',
      version=read_version(),
      description='',
      long_description=readme(),
      long_description_content_type='text/markdown',
      author='Ebony Watson',
      author_email='ebonyrwatsion@gmail.com',
      url='https://github.com/Ebony-Watson/scProximitE',
      license='GPL3',
      project_urls={
          "Bug Tracker": "https://github.com/Ebony-Watson/scProximitE/issues",
          "Documentation": "https://github.com/Ebony-Watson/scProximitE",
          "Source Code": "https://github.com/Ebony-Watson/scProximitE",
      },
      classifiers=[
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Natural Language :: English',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7',
          'Programming Language :: Python :: 3.8',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
      ],
      keywords='util',
      packages=['dimred'],
      entry_points={
          'console_scripts': [
              'dimred = dimred.__main__:main'
          ]
      },
      install_requires=['pandas', 'numpy', 'tensorflow', 'stats', 'seaborn', 'matplotlib', 'sklearn', 'jupyterlab',
                        'scanpy==1.8.2', 'sciviso', 'scikit_posthocs', 'torch', 'genieclust', 'natsort',
                        'otscomics @ http://github.com/SergeySatskiy/cdm-pythonparser/archive/v2.0.1.tar.gz'],
      python_requires='>=3.6',
      data_files=[("", ["LICENSE"])]
      )
