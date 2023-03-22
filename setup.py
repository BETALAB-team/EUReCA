from setuptools import setup, find_packages

classifiers = [
  'Development Status :: 1 - Planning',
  'Intended Audience :: Education',
  'Operating System :: OS Independent',
  'License :: OSI Approved :: MIT License',
  'Programming Language :: Python :: 3.8'
]

setup(
   name='eureca-ubem',
   version='0.1.0',
   author='betalab group UNIPD',
   author_email='enrico.prataviera@unipd.it',
   packages=find_packages(),
   scripts=[],
   url='https://github.com/BETALAB-team/EUReCA',
   license='LICENSE',
   description='Package to simulate district heating cooling systems for eureca',
   long_description=open('README.md').read(),
   install_requires=[
       "numpy",
       "pandas",
       "matplotlib",
       "geopandas",
       "pytest",
       "cjio",
       "shapely",
       "eureca-building",
   ],
)