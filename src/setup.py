import setuptools


setuptools.setup(
	name="splam",
	version="0.1.0",
	author="Kuan-Hao Chao",
	author_email="kh.chao@cs.jhu.edu",
	description="Splice junction scoring tool",
	url="https://github.com/Kuanhao-Chao/SPLAM",
	install_requires=['torch>=1.12.0', 'pandas>=2.0.0', 'pybedtools>=0.9.0'],
	python_requires='>=3.6',
	packages=['splam'],
	entry_points={'console_scripts': ['splam = splam.run_splam:main'], },
)
