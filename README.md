# Ground Filtering

**--INFORMATION**

The _Tin Refinement_ algorithm is from my collaborator Yitong Xia(xiayitong0630@gmail.com), thank her for her efforts and dedication.

_CSF_ algorithm is based on the cloth simulation with verlet integration.

The idea of _CSF_ is from this paper : http://www.mdpi.com/2072-4292/8/6/501/htm

The basic architectures and ideas are inspired by this article (its related code is open source).

**--HOW TO USE**

It's a cmake project, in principle, it can run on all platforms if everything goes fine. Use the source code to build the project and run it.

**--EXAMPLE**

**Input data set**

![image](https://user-images.githubusercontent.com/72781910/149506682-5fb2b4d1-6caa-4480-9a9e-f4f17e281dca.png)

**Output data set**

![image](https://user-images.githubusercontent.com/72781910/149506841-ed5a9585-05f7-43a5-ac06-83760ad73c32.png)

The original data set is classified as **_ground points_** and **_non-ground points_**.

**--INPUT**

The input files are in **_data_** folder.

**--OUTPUT**

The output files are in **_data_** folder. The names of generated files are in: **_data/params.json_**.

**--ADJUSTMENT**

Modify input file and parameters: enter into **_data/params.json_**.






