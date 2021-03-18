from setuptools import setup, find_packages



setup(
    name='geoline',
    version='0.0.1',
    packages=find_packages(),
    url='https://github.com/mountetna/monoetna/geoline',
    license='MIT',
    author='Anton Gvaihir Ogorodnikov',
    author_email='anton.ogorodnikov@ucsf.edu',
    description='Interactive submissions to GEO using data in DSCo Lab Data Library',
    python_requires='>=3.6',
	install_requires=['betamax',
                      'betamax_serializers',
                      'pandas'],

)