# Contributing

When contributing to this repository, please first discuss the change you wish to make via issue,
email, or any other method with the maintainers of this repository. This will make life easier for everyone.

## Report Issues

Please use the [issue tracker](https://gitlab.com/vibes-developers/vibes/-/issues) to report issues. Please try to answer these questions:

- Has this issue been discussed before? Please have a quick look at the existing issues. If not:
- What is the issue? What is the expected behavior?
- Is the problem reproducible? Please provide a _minimal_ example.


## Contribute Code via Merge Request

In order to contribute code to `FHI-vibes`, please follow the usual steps for [preparing and creating a merge request](https://docs.gitlab.com/ee/user/project/merge_requests/creating_merge_requests.html). A few remarks regarding our guidelines for code and code style:

- We use [black](https://black.readthedocs.io/en/stable/) with default settings and [isort](https://pycqa.github.io/isort/) for formatting the code. The settings for `isort` are included in `setup.cfg`.
- Please _document_ and _test_ your changes. Tests are found in `vibes/tests` and written with [pytest](https://docs.pytest.org/en/stable/).
- Please use [google-type docstrings](https://google.github.io/styleguide/pyguide.html) for your functions. Optionally you can use type hints, but we currently don't enforce this.
- We loosely keep track of code coverage, please try not to decrease coverage when contributing new code.
