from setuptools import Extension, setup


setup(
    ext_modules=[
        Extension(
            "fastop._native",
            sources=["src/fastop/_native.c"],
        ),
    ],
)
