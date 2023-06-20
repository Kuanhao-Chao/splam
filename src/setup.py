from setuptools import setup, Extension


# Define the extension module
# extension_module = Extension('splam_extract', sources=['your_module.c'], libraries=['build/splam_extract'])

setup(
	name="splam",
	version="0.1.0",
	author="Kuan-Hao Chao",
	author_email="kh.chao@cs.jhu.edu",
	description="Splice junction scoring tool",
	url="https://github.com/Kuanhao-Chao/SPLAM",
	install_requires=['torch>=1.12.0', 'pandas>=2.0.0', 'pybedtools>=0.9.0'],
	python_requires='>=3.6',
	packages=['splam'],
    package_data={'splam_extract': ['build/splam_extract.cpython-38-darwin.so']},
    # ext_modules=[extension_module],
	entry_points={'console_scripts': ['splam = splam.splam:main'], },
)