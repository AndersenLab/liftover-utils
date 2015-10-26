from setuptools import setup
import glob

setup(name='liftover.py',
      version='0.1',
      packages=['liftover'],
      description='C. elegans liftover utility',
      url='https://github.com/AndersenLab/liftover-utils',
      author='Daniel Cook',
      author_email='danielecook@gmail.com',
      license='MIT',
      entry_points="""
      [console_scripts]
      liftover = liftover.liftover:main
      """,
      install_requires=["docopt"],
      data_files=[('CHROMOSOME_DIFFERENCES', glob.glob("data/CHROMOSOME_DIFFERENCES/sequence*")),
      			   'remap_gff_between_releases.pl'],
      zip_safe=False)