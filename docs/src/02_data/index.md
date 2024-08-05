# [Data formatting and outputs](@id data)
## [Data and annotation formatting](@id data)

KinBiont can operate directly on data files or inside the julia notebook.
When are in a julia notebook the  format of single time series that want to be analyzed is a 2 x n_time_points Matrix of FLoat64, e.g.,


```
 0.0        2.0       4.0       6.0       8.0        10.0       10.0       12.0       14.0       16.0       18.0       20.0       22.0      24.0      26.0       28.0       30.0       32.0       34.0       36.0       …  
 0.0912154  0.107956  0.105468  0.101727  0.0931484   0.106318   0.103697   0.139821   0.173598   0.204888   0.251052   0.289018   0.31298   0.33752   0.359356   0.370861   0.376347   0.383732   0.398496   0.384511 …  

```
The first row should be time and the second the quantity to be fitted (e.g., Optical Density or CFU)

Instead, three APIs call direclty the files: the user must input  the paths to  a .csv data file and a .csv annotation to the functions of KinBiont.jl
; In these cases KinBiont expect for data a matrix where the first row are the names of the wells and the columns the numerical value of the measurements. Note that the first one will be used as time:

```
Time,  A1,     A2,      A3, 
0.0,   0.09,   0.09,    0.087,
1.0,   0.08,   0.011,   0.012,
2.0,   0.011,  0.18,    0.1,
3.0,   0.012,  0.32,    0.22,
4.0,   0.008,  0.41,    0.122,
```
KinBiont expect a "," as separator between columns

The annotation file instead should be a two columns .csv file where the number of rows correspond to the number of wells, note that the name of the well should be the same between the data.csv and annotation.csv:

```
A1, b
A2, X
A3, unique_ID

```
as unique_ID the user can insert anything but consider that if two wells has the same ID the will be considered replicates. 'b' indicates that the well should be cosidered a blank and 'X' that the well should be discarded from any analysis


See the folders  XXXXX for some examples. 


If a OD calibration curved is provided it should have the following format XXXXXXX


## [Outputs of KinBiont]
