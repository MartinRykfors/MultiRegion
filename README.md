# MultiRegion
## A cache efficient solver for max-flow/min-cut problems on grid graphs.

This was my 2012 master's thesis project. My thesis publication can be read at [https://lup.lub.lu.se/student-papers/search/publication/3233492](https://lup.lub.lu.se/student-papers/search/publication/3233492).

Its purpose is to solve max-flow problems on flat 3D grids that are complete along the last dimension. This allows for more complicated image segmentations using a single graph cut. For inctance, it allows splitting an image into three segments (two objects and background) optimally. The technique is described in the article "Globally Optimal Segmentation of Multi Region Objects" by A.Delong and Y.Boykov, [http://www.psi.toronto.edu/~andrew/](http://www.psi.toronto.edu/~andrew/).

It implements a slightly modified version of the algorithm by A. Golberg et al. [research.microsoft.com/pubs/150437/ibfs-proc.pdf](http://research.microsoft.com/pubs/150437/ibfs-proc.pdf).

Currently there are a variety of small examples showing how the MRGraph class is used.

Example 3 and the benchmark requires the GridCut headers that can be downloaded from [www.gridcut.com](http://www.gridcut.com).

This project was developed before I had a good understanding on how to structure an open-source project. I am currently revisiting it to make it easier to use for everyone. Work in progress is to replace the examples with unit tests, and to provide an easier way of building it.

## License

MIT
