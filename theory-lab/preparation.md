# Preparation

For the theory part of this LO1 project, we will use the programming language Python. I ask you to have a working installation of Python before LO1 starts. The less time we need to spend installing stuff, the more time we can spend on chemistry.

I recommend using Miniconda as Python distribution and Visual Studio Code for writing and running code. You can find a short manual below.

## Miniconda 🐍

Python's success is due to many free packages being available. For example, Numpy makes mathematics fast and the Atomic Simulation Environment (ASE) has many tools for chemistry simulations. Sometimes different packages are incompatible with each other. It can then be useful to have different Python installations ("environments") that are independent of one another. Conda is a popular tool for managing Python installations.

To get Conda you can install Anaconda or Miniconda. Anaconda contains some unnecessary extras whereas Miniconda is a more lightweight version. If you already have an Anaconda installation you don't need to install Miniconda.

To install Miniconda:

1. Download the Miniconda installer from the official website: https://www.anaconda.com/download/success.
2. Run the downloaded installer with the default options (you don't need to check or uncheck any checkboxes).

If the installation fails on Mac you could try the Anaconda installer instead (same link as before, but get the 'Distribution' rather than the 'Miniconda' installer).

You can check if your installation works by searching for the 'Anaconda prompt' (Windows logo key + S, then type `prompt` for example). On Mac, you can just open the Terminal.

Once in the prompt or terminal, you can type `conda list` and it should list some  packages installed in your 'base' (default) environment, including a `python` package.

If you want, you can already type `pip install mace-torch` to install a package that we will need. The dependencies such as `torch` should be installed automatically.

## Visual Studio Code ⌨️
VS Code is arguably the most popular software for writing code. If you already have an editor like Spyder installed, you can use it, but I highly recommend VS Code. It has powerful features like interactive notebooks and AI coding suggestions, and it's easier to help with errors if everyone in the project uses the same editor.

To install VS Code:

1. Download the installer from the official website: https://code.visualstudio.com/
2. Run the installer and follow the steps with default settings.

If you want, you can already get used to VS Code by following a tutorial like [this YouTube video](https://youtu.be/6i3e-j3wSf0?si=HsgdzThwTTn1JvZb). (You can skip their installation instructions if you already installed Ana/Miniconda).

Some tips for VS Code and coding in general:

* Create a separate folder for your project, and open this folder in VS Code (`File -> Open Folder...`). This way you don't lose track of where you saved your project files.

* When there's a white dot next to the filename on top of the screen, your file contains unsaved changes.

* Save quickly with Ctrl + S.