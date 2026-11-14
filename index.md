# Fun With *Psi4*

## Preface

This resource is intended to help organic chemistry students set up the [***Psi4*** computational chemistry software](https://psicode.org) {cite:p}`psi4` on their own computer and perform a series of small calculations. It supports my course in **physical organic chemistry**, CHEM 4410. We will be using computational chemistry to demonstrate some important concepts in **organic chemistry**, and we won't be thinking too much about how it works. There is another course available to students that explores computational chemistry in detail. CHEM 4320 will take you much farther if you catch the bug here. I hear that CHEM 4320 is all math, but CHEM 4410 is all fun. Let us have some **fun with *Psi4***.

## How to Use this Book

Navigate to the chapter of your choice using the **left sidebar** or the hamburger menu in the upper left corner. Most of the pages in this book are [*Jupyter* notebooks](https://jupyter.org) {cite:p}`granger2021jupyter` that you can opened in [*Google Colab*](https://colab.research.google.com) {cite:p}`googlecolab2024` to run the code. You can also **follow the instructions** in this book to set up your own personal computer to run *Python* and *Psi4*. **Have fun**!

This book is a *MyST* book and is produced using markdown files and *Jupyter* notebooks that are converted to a web book format using the [*MyST* parser](https://mystmd.org) {cite:p}`mystmd_cli`.

%```{figure} images/SN2_2_copy.gif
%```

## The Origin of this Book

Many years ago **I had a dream** of incorporating the topic of **computational chemistry** into my course in physical organic chemistry. I spent some time over a few summers banging my head against a wall as I **learned** how to use the powerful *GAMESS* package. I built a web book on [using *GAMESS*](https://blinkletter.github.io/UsingGAMESS) to perform various exercises in computational chemistry. The goal was **never achieved**. Students had to install *GAMESS* on their own computers and this turned into a neverending problem of tech support. 

A few years ago **I once again revisited the idea**. Over the years I had become familiar with *Python* and chemical research tools that it could control. The *Psi4* package can be run as a *Python* module. I can now present computational chemistry from **within interactive *Jupyter* notebooks** running *Python*. These notebooks can be presented as a web book and linked through to *Google Colab* so that students can run them **without needing to install** anything on their own computer.

### An AI Heel Turn

So, I finally paid for a subscription to an **LLM system** and tasked it with **converting** my previous notebooks to *Python* versions running *Psi4* instead of describing how to manage *GAMESS* via *Unix* commands. It did an **amazing job**. That was version one of my [web book for using Psi4](https://blinkletter.github.io/UsingPsi4). I also tasked it with **rewriting** the original book to use GAMESS but with *Python* code controlling the system from within a notebook. I was stunned to see it produce an excellent **re-interpretation** of the original book. That became [Using *GAMESS*, v2](https://blinkletter.github.io/UsingGAMESSv2).

I **learned an enormous amount** about *Python*, *Psi4* and *GAMESS*. The LLM was able to read the documentation and the source code of both packages. I was able to **solve issues** that had prevented me from compiling *GAMESS* on newer MacOS computers and **delve deeper** into the options available in *Psi4*. The code produced by the LLM was very sophisticated and taught me new ways to approach using *Python* as an interface to both these systems. However, both of those books might be **too sophisticated** for the absolute beginner. I decided to use the two revised web books as source material for a **new book** describing how to use *Psi4* in my course.

### My Redemption Arc

Welcome to **Using Psi4, 2{sup}`nd` ed.**. After reading and editing the two converted books, and rereading the original, I have observed lots of room for improvement. I am using the original writing from the first book, adapting the text written by the LLM that describes the *Python* code and **copying and pasting** much of the *Python* code that was produced by the LLM when it converted jobs meant to run in *GAMESS* in *Unix* to jobs that run in *Python* notebooks controlling *Psi4*. I tried to **simplify** the *Python* so that a student who is new to computer programming can follow along.

In the end, **this book** will perhaps be a 50/50 effort between what I **steal** from my original book and what I use from the two books that were adapted from that original source material by the LLM system. I will ask the LLM to analyze the new book **compared** to the original and the two conversions to create an "AI statement" for each chapter that will describe **how much** came from my **original work**, how much is from the LLM work in the converted books, and how much is my new **original writing** in this book.

```{attention} AI Statement
:class: dropdown
This page is my own **original writing**. The preface is **adapted** from the front pages of my earlier web books, [*Fun with GAMESS*](https:/blinkletter.github.io/UsingGAMESS) and the [first edition of *Fun with Psi4*](https:/blinkletter.github.io/UsingPsi4). The rest was written **new** for this edition. No text on this page was generated by an LLM. Claude (Anthropic) with Claude Code was used to **compare** this page with the earlier versions.