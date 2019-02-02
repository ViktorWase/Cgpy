import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="cgpy",
    version="0.0.2",
    author="Viktor Wase",
    author_email="viktorwase@gmail.com",
    description="A Python library for Cartesian genetic programming.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="http://wase.io/",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)