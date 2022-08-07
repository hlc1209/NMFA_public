# NMFA

Version 11/25/2021

This is the code repository for the paper "Deciphering the generating rules and functionalities of complex networks" 

https://www.nature.com/articles/s41598-021-02203-4

## File Structure
```
├── Tutorial.ipynb          // Quick Start Guide
├── NFD.py                  // Core Function Library
└── README.md
```

## Requirements

- Python 3.6 or higher (If you want to use Networkit, use python 3.7.x)
- Networkx 
- Networkit (only for Linux and MacOS according to its official guide)
- tqdm
- scipy


## Reference of Core Functions
***
### node_dimension_single(G, node, weight=None, fig=True)
Calculate the NFD value of the target node
#### input
G: (Networkx Graph) - Input graph

weight: (String) - Edge attribute of weight.

node: Label of the target node

#### output

[1] (double) - NFD value

[2] (double) - r value of linear regression

***
### nfd(G, Q = [q/100 for q in range(-2000,2001,10)])
Calculate tau of unweighted network
#### input
G: (Networkx Graph) - Input graph

Q: (list) - Distortion exponents

#### output

[1] (numpy.array) - tau_list

***
### nfd_nk(G, Q = [q/100 for q in range(-2000,2001,10)])
Calculate tau of unweighted network using networkit package. Much faster but only available on Linux and MacOS
#### input
G: (Networkx Graph) - Input graph

Q: (list) - Distortion exponents

#### output

[1] (numpy.array) - tau_list

***
### nfd_multiprocess(G, Q = [q/100 for q in range(-2000,2001,10)])
Calculate tau of unweighted network using python built-in Multiprocessing. Much faster but only available on Linux and MacOS
#### input
G: (Networkx Graph) - Input graph

Q: (list) - Distortion exponents

#### output

[1] (numpy.array) - tau_list

***
### wnfd(G, Q = [q/100 for q in range(-2000,2001,10)], fdigi=3)
Calculate tau of weighted network
#### input
G: (Networkx Graph) - Input graph

Q: (list) - Distortion exponents

fdigi: (int) - Precision of path length

#### output

[1] (numpy.array) - tau_list

***
### wnfd_nk(G, Q = [q/100 for q in range(-2000,2001,10)], fdigi=3)
Calculate tau of weighted network using networkit package. Much faster but only available on Linux and MacOS
#### input
G: (Networkx Graph) - Input graph

Q: (list) - Distortion exponents

fdigi: (int) - Precision of path length

#### output

[1] (numpy.array) - tau_list

***
### nspectrum(tau_list,q_list=[q/100 for q in range(-2000,2001,10)],name=None,linewidth=5)
Calculate $\alpha$, $f(\alpha)$ and asymmetry, and draw multifractal spectrum
#### input
tau_list: (list) - tau

q_list: (list) - Distortion exponents

name: (String) - Name of the network (for legend)

linewidth: (int) - line width in figure

#### output

[1] (numpy.array) - $\alpha$

[2] (numpy.array) - $f(\alpha)$

[3] (double) - asymmetry

***
### ndimension(tau_list,q_list=[q/100 for q in range(-2000,2001,10)],name=None,linewidth=5)
Calculate $D(q)$, and draw generalized fractal dimension
#### input
tau_list: (list) - tau

q_list: (list) - Distortion exponents

name: (String) - Name of the network (for legend)

linewidth: (int) - line width in figure

#### output

[1] (numpy.array) - $D(q)$

[2] (numpy.array) - Distortion exponents

***
### distance(ndim_list)
Calculate structure distance
#### input
ndim_list: (list) - $D(q)$

#### output

[1] (numpy.array) - N by N matrix of structure distance

***
### nheat(tau_list, q_list=[q/100 for q in range(-2000,2001,10)])
Calculate specific heat
#### input
tau_list: (list) - tau

q_list: (list) - Distortion exponents

#### output

[1] (numpy.array) - Specific heat

[2] (numpy.array) - Distortion exponents