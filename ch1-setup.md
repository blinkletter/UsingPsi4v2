# 1. Building Your System

First we must build a system of ***tools*** for computational chemistry. You can **acquire** the software and **install** it on your computer using the following **instructions**. Instructions are given for both *Windows* and *MacOS* computers.

```{attention} Stop! You don't have to do any of this
:class: dropdown
You do not need to install anything at all. **All you really need** is the web browser that you are using to read this book. Every notebook with *Psi4* exercises can be run in your **web browser** using [Google Colab](https://colab.research.google.com). This is useful if your primary computing device is an iPad, a phone, or a laptop computer. The rest of this setup page is for those enthusiasts who want to use *Psi4* on their own **personal computer**. 
```

We have some intuition about **chemical structure and energy**. We know that butane is at its lowest energy in the *anti* conformation and that there is a small **energy barrier** to bond rotation due to steric clash. We know that the amide bond the nitrogen is planar, rather than tetrahedral, and has a higher energy barrier to rotation due to its **electronic structure**. We can demonstrate these effects and **explore consequences of changes in structure** on our desktop using a computer and software that can calculate the energy of molecules. There are software packages that cost many thousands of dollars (per year!) and there are options that are **free to all** or free for academic use.

## The Software

I have lived in the *Unix* world my whole academic life and am largely unfamiliar with *Windows* systems. However, the tools below all work on *Windows*, *MacOS* and *Linux*, and the websites have information on installing them for each system. Where the steps differ, I have given separate instructions for *Windows* and *MacOS*.



### Psi4

***Psi4*** {cite:p}`psi4` ([Website](https://psicode.org)) is an open source *ab initio* computation package that performs **quantum mechanical calculations** for investigating molecular structure and energy. 

There are others options, such as ***GAMESS*** {cite:p}`gamess1`, {cite:p}`gamess2`, ***NWChem*** {cite:p}`nwchem` and ***ORCA*** {cite:p}`orca`. *Psi4* is different from these in that it is installed as a ***Python* module**. We control it by writing *Python* code, and the results come back to us as *Python* variables that we can process, tabulate and plot using all the other **tools** in the *Python* ecosystem.

### Python and Anaconda

*Python* is a programming language that has become the common language of scientific computing. [***Anaconda***](https://www.anaconda.com) is a **distribution** of *Python* that comes with a package manager called ***conda***. *Conda* will download and install *Psi4* and all the other packages it depends on. It also lets us keep the *Psi4* tools in their own separate **environment** so that they do not interfere with anything else on your computer. *Anaconda* can be obtained [here](https://www.anaconda.com/download).


### Jupyter

*Jupyter* notebooks let us mix **text, code and results** in a single document. We will write our *Psi4* calculations in notebooks, run them one cell at a time, and see the results appear right below the code. *Jupyter* is installed with *conda* along with *Psi4*.


### Avogadro

***Avogadro*** {cite:p}`avogadro` ([Website](https://avogadro.cc)) is a modern, open source application for **building molecules** and **visualizing** the results of calculations. 

You can use *Psi4* and *Avogadro* together. You can build larger molecules using Avogadro and export the input file to use in *Psi4*. You can then run *Psi4* from ***Jupyter* notebooks** using *Visual Studio Code*. After that you can load in various result files to visualize the optimized structure, molecular orbitals and surfaces.

In this note-based web book, all structures will be entered directly from the notebooks and molecular visualization will be created in the notebook itself using *Python* and *py3Dmol* library. I will note when the input file was built with the help of *Avogadro* or when I am displaying images created using *Avogadro*. 

### py3Dmol

[*py3Dmol*](https://3dmol.csb.pitt.edu) is a *Python* package that draws **interactive 3-D molecules** right inside a *Jupyter* notebook. You can rotate and zoom the structures with your mouse. We will use it to quickly look at our molecules before and after a calculation without leaving the notebook. *py3Dmol* is installed using *conda*, like *Psi4*. 


### Other Required Tools

We will need a text editor. Microsoft provides a **free cross-platform text editor** designed for writing code called *[Visual Studio Code](https://code.visualstudio.com/)*. It can also open and run *Jupyter* notebooks. I will be using *Visual Studio Code* in my demonstrations as it works on all our computers.

You will need a **Terminal Emulator**. We will be typing a few commands to install the software. On a *MacOS* system there is an application called *Terminal* that will allow us to access the real computer that lives inside your computer. On *Windows*, *Anaconda* provides an application called *Anaconda Prompt* for this purpose. *Visual Studio Code* also has a built-in terminal emulator.

## Installing *Anaconda*

*Anaconda* provides *Python* and the *conda* package manager. We will use *conda* to install everything else.

### Obtain and Install

Go to the [*Anaconda* download page](https://www.anaconda.com/download) and download the **installer** for your computer. Be sure to choose the correct one for your computer. Follow the instructions for your operating system below.

```{admonition} MacOS
:class: dropdown
Download the **graphical installer** (`.pkg`) for *MacOS*. If you have a newer Mac with an Apple processor (M1, M2, M3...) choose the *Apple Silicon* version. Older Macs with Intel processors need the *Intel* version. You can check by selecting `Apple menu` {math}`\rightarrow` `About This Mac`.

Double-click the downloaded file and follow the instructions. Accept the defaults. When it is done, open the *Terminal* application (it is in `/Applications/Utilities`). You should see `(base)` at the start of the command prompt. This means that *conda* is working and the base environment is active.
```

```{admonition} Windows
:class: dropdown
Download the **graphical installer** (`.exe`) for *Windows*. Double-click the downloaded file and follow the instructions. Install for **"Just Me"** and accept the default location. It is best **not** to check the box that adds *Anaconda* to your `PATH`.

When it is done, open the *Start* menu and find the application called ***Anaconda Prompt***. This is the terminal we will use on *Windows*. You should see `(base)` at the start of the command prompt. This means that *conda* is working and the base environment is active.
```

You can test the installation by typing the following at the command prompt in the terminal emulator.

```bash
conda --version
python --version
```

```{note}
:class: dropdown
If you would rather not install the whole *Anaconda* distribution, which is large and comes with many packages that we will not use, you can install ***Miniconda*** or ***Miniforge*** instead. These are small installers that provide only *Python* and *conda*. Every *conda* command below works the same way.
```

## Installing *Psi4* and *Jupyter*

*Psi4* is distributed as a *conda* package on the ***conda-forge*** channel. We will create a new **environment** called `psi4` that contains *Psi4*, *Jupyter* and the other *Python* packages we will use for analyzing and plotting our data. The same commands work in the *MacOS Terminal* and the *Windows Anaconda Prompt*.

### Create an Environment

Type the following command. It will take a few minutes while *conda* works out all the packages that are needed and downloads them. Answer `y` when it asks if you want to proceed.

```bash
conda create -n psi4 -c conda-forge psi4 jupyter ipykernel \
                        numpy pandas matplotlib py3dmol
```

| Package       | Description                                                       |
| :------------ | :---------------------------------------------------------------- |
| *psi4*        | The *Psi4* quantum chemistry package                              |
| *jupyter*     | The *Jupyter* notebook system                                     |
| *ipykernel*   | Lets *Visual Studio Code* run notebooks using this environment    |
| *numpy*       | Numerical arrays and math                                         |
| *pandas*      | Tables of data                                                    |
| *matplotlib*  | Plotting                                                          |
| *py3dmol*     | Interactive 3-D molecule viewer for notebooks                     |

### Activate the Environment

Every time you open a new terminal you must **activate** the environment before using *Psi4* from the command line.

```bash
conda activate psi4
```

The command prompt will change from `(base)` to `(psi4)`. To leave the environment, type `conda deactivate`.

### Test the Installation

With the environment active, start *Python* and ask *Psi4* to calculate something small.

```bash
python -c "import psi4; print(psi4.__version__)"
```

You should see the version number of *Psi4* printed. 

## Installing *Visual Studio Code*

*Visual Studio Code* is a free text editor from Microsoft that is designed for writing code. With a couple of **extensions** it becomes a very nice place to run *Jupyter* notebooks.

### Obtain and Install

Download *Visual Studio Code* from the [project website](https://code.visualstudio.com/). Follow the instructions for your operating system below.

```{admonition} MacOS
:class: dropdown
Download the *MacOS* version. It arrives as a `.zip` file that expands to the *Visual Studio Code* application. Drag the application into your `/Applications` folder. 
```

```{admonition} Windows
:class: dropdown
Download the *Windows* **User Installer**. Double-click the downloaded file and follow the instructions. Accept the defaults.
```

### Install the Extensions

Open *Visual Studio Code* and click on the **Extensions** icon in the bar on the left side of the window (it looks like four squares). Search for and install the following two extensions, both published by Microsoft.

| Extension | Description                                              |
| :-------- | :------------------------------------------------------- |
| Python    | Support for the *Python* language                        |
| Jupyter   | Support for opening and running *Jupyter* notebooks      |


*Visual Studio Code* has a **built-in terminal** that can be opened using `View` {math}`\rightarrow` `Terminal`. On *MacOS* it will behave like the *Terminal* application. On *Windows* you may have to type `conda activate psi4` in this terminal before *conda* commands will work. If they don't, use the *Anaconda Prompt* instead.

## Installing Avogadro

*Avogadro* is an open source application for **building molecules** and can export input files and interpret output files for many computational chemistry packages. It has many interesting features to explore. One very useful tool is the auto-optimization feature that will minimize the energy of your molecule using **molecular mechanics**. This will produce a structure that is closer to ground state and therefore faster to optimize in *Psi4*. It can also display **molecular orbitals** and **vibrations** from the files that *Psi4* writes.

### Obtain and Install

You can [download *Avogadro*](https://avogadro.cc) from the project's website. Follow the instructions for your operating system below.

```{admonition} MacOS
:class: dropdown
Download the *MacOS* `.dmg` file. Be sure to choose the correct one for your processor (*Apple Silicon* or *Intel*). Open the `.dmg` file and drag *Avogadro* into your `/Applications` folder.

The first time you open it, *MacOS* may complain that it cannot verify the developer. If so, right-click on the application icon, choose `Open` from the menu, and then click `Open` again in the warning box. You will only need to do this once.
```

```{admonition} Windows
:class: dropdown
Download the *Windows* installer (`.exe`). Double-click the downloaded file and follow the instructions. Accept the defaults. *Windows* may warn you about an unrecognized application. If so, click `More info` and then `Run anyway`.
```

Install the software and that's it. Open it and **explore** the menus.

### Using Avogadro

An [instruction manual](https://avogadro.cc) is available at the same website. We will explore *Avogadro* as we perform our experiments. There is no substitute for trial and error.

## Summary of System Setup

The following software was installed.

| Software              | Where it lives                    | How it is used                         |
| :-------------------- | :-------------------------------- | :------------------------------------- |
| *Anaconda*            | Your home directory               | Provides *Python* and *conda*          |
| *Psi4*, *Jupyter* etc. | The `psi4` *conda* environment   | Activate with `conda activate psi4`    |
| *Visual Studio Code*  | Applications folder      | Edit and run notebooks (kernel `psi4`) |
| *Avogadro*            | Applications folder      | Build and visualize molecules          |

You can **create your own setup** however you like. Just remember to select the `psi4` environment whenever you run a notebook.

## Next Steps

Next you should **test** that everything works. Run the *VSCode* and *Avogadro* programs. Then follow the tests in the next section to **confirm** that *Psi4* is installed. The next page will be the first executable *Jupyter* notebook in this web book. there you will run **your first test calculation**.


```{attention} AI Statement
:class: dropdown
TBA
```
