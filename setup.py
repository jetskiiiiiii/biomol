from setuptools import find_packages, setup

setup(
    name='biomol',
    packages=find_packages(include=['biomol']),
    version='0.1.0',
    description='A collection of functions related to molecular biology.',
    author='Treverrio Primandaru',
    install_requires=[],
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    test_suite='tests',

)