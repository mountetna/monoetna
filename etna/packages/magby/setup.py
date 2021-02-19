from setuptools import setup, find_packages

setup(
    name='magby',
    version='0.0.1',
    packages=find_packages(),
    url='https://github.com/mountetna/monoetna/etna/packages/magby',
    license='MIT',
    author='Anton Gvaihir Ogorodnikov',
    author_email='anton.ogorodnikov@ucsf.edu',
    description='Python client for magma DB',
    python_requires='>=3.6',
	install_requires=['requests',
                      'pandas',
                      'betamax'],

)
