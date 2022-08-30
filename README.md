[![halem status](https://github.com/TUDelft-CITG/halem/actions/workflows/action.yml/badge.svg)](https://github.com/TUDelft-CITG/halem/actions/workflows/action.yml)
[![halem coverage](https://TUDelft-CITG.github.io/halem/coverage/coverage.svg)](https://TUDelft-CITG.github.io/halem/coverage/)

# HALEM

**H**ydrodynamic **A**lgorithm for **L**ogistic **E**nhancement **M**odule: route optimisation in dynamic currents.

Documentation can be found [here](https://halem.readthedocs.io).

## Features

This package contains route optimization for given currents. The following features are taken into account in this version:

* Spatial varying currents
* Temporal changing currents
* Variable shipping velocity
* minimal water depth
* Squad

Does not take into account:

* Inertial behavior of the ship

Different routes that can be optimized are:

* Shortest route (halem.HALEM_time)
* Fastest route (halem.HALEM_space)
* Cheapest route route (halem.HALEM_cost)
* Cleanest route (halem.HALEM_co2)

## The Python package

The Python package is setup using `PyScaffold`:

```bash
pip install pyscaffold
putup halem
```

Make sure not to run the `putup` command in a Git repository clone as it will raise an exception. If you already have a Git repo to store your package, run the `putup` command outside the Git repository clone and copy the generated contents into that directory.

`PyScaffold` provides you with a full Python package, including source, test and doc directories. You can create a development environment using:

```
cd halem
pip install -e .[testing]
```

## The Test environment

The `.github/workflows/action.yml`, `pytest-full.ini` and parts of the `setup.cfg` files form your test environment. Your test environment is built around the `pytest` framework. We use **fixtures** to configure our tests, please read about them. Also, there are basically two test environments:

- `local`:

  Your local test environment that is initiated by running `pytest`. It is meant to be fast, but not exhaustive. It only runs the tests you wrote yourself.

- `full`:

  Your full test environment that is used on Github to test your code more thoroughtly, including code style and consistency checks. It is initiated by running `pytest -c pytest-full.ini`. It is meant to be exhaustive, but not fast.
  
The full test environment enables various extensions that ensure that your code is correct, consistent and well readible:

- `black`

  Code formatter. Your code should be formatted according to Black's rules. If you are working in VSCode, the formatting will be done automatically for you.

- `isort`

  Import sorter. Your imports should be at the top of your module, in one particular order. If you are working in VSCode, the sorting will be done automatically for you.

- `flake8`

  Code conciseness checker. Your code should be concise. We don't want variables to be defined, while not being used. We don't want modules be imported, while not being used. We don't want code that can never be reached. Et cetera. If you are working in VSCode, you will be warned automatically if your code fails the conciseness checks.

- `mypy`

  Code consistency checker. Your code should be consistent. If you pass an integer to a function that expects a string, something is wrong. If you return a dictionairy, while a float was expected, something is wrong. If you are working in VSCode, you will be warned automatically if your code fails the consitency checks.

- `coverage`

  Code coverage checker. Your code should be tested. If your test do not use all of your code, users of your package might get unexpected results. Your code coverage percentage tells you how much of your code is actually used in your tests. It should at least be 70%, but preferably higher. The code coverage tools shows you how much of your code is tested and what code requires more tests.

Github runs the full test environment on several combinations of operating systems and Python versions. By default all combinations of the latest versions of the operating systems Windows, MacOS and Ubuntu, and Python versions 3.8, 3.9 and 3.10 are tested. Github also publishes teh coverage results and documentaion to your Github project page: https://TUDelft-CITG.github.io/halem/docs/ and https://TUDelft-CITG.github.io/halem/coverage if enabled on Github (`Settings` -> `Pages`). See `.github/workflows/action.yml` for details.

## The Dev environment

The `.vscode` directory forms your Dev environment. It contains configuration files that makes Visual Studio Code:

- Install your package in development mode.
- Constantly performs above-mentioned checks and balances for easy, consistent code development.
- Allows you to quickly test parts of your code.
- Allows you to quickly debug parts of your code.

It is recommended to create a virtual Python environment for each Python package you develop. Make sure that the paths in your `.vscode/settings.json` file point to the right interpreters. To enable automated tasks in Visual Studio Code choose `Ctrl+Shift+P` -> `Tasks: Manage Automated Tasks in Folder` -> `Allow Automated Tasks in Folder`.

## Python package installation

Once your code is on Github, you can tag it with a semantic version. A semantic version had the form `v1.0.0`, where the first number indicates a major revision, the second a minor revision, and the third a patch. The use of increments is a but subjective, but a useful rule of thumb is:

- `major`:

  Breaking changes, not compatible with previous versions. Users should update their code and environments that depend on your package.

- `minor`:
  
  No breaking changes, but new features are added. Users should review the changes to see if they can optimize their usage.

- `patch`:

  Small changes and bug fixes. Users should be able the update their code without review.

Once your code has a semantic version tagged to it, you can install it using `pip` (Linux/Docker syntax):

```bash
pip install halem
```

Or on Windows/Anaconda:

```bash
pip install halem
```

## The Dockerfile (optional)

The `Dockerfile` and `docker-compose.*.yml` files allow you to run your code in a containerized environment. A container is like a virtual machine. It provides a workspace for you to work in, completely separate from your own laptop. The advantage of containers is that you, your colleagues and the test servers are all working in exactly the same environment. The disadvantage is that Docker requires a lot of resources from your laptop, and, if you don't know how to run Docker on bare WSL, you need a paid license.

Using containers locally is not required. But if you are used to work with Docker it might be convenient.

### The Dev container

The `.devcontainer` directory provides a Visual Studio Code Dev environment. It leverages your Docker container environment and your Test environment in a full fletched automated Dev environment that is:

- Guaranteed to be identical as on the Github test servers, and your colleagues laptop.
- Constantly performs above-mentioned checks and balances for easy, consistent code development.
- Allows you to quickly test parts of your code.
- Allows you to quickly debug parts of your code.
