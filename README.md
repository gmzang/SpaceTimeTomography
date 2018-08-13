----------------------------
Space-time Tomography for Continuously Deforming Objects
----------------------------
`This is the open source repository for paper in SIGGRAPH 2018`:

[**Space-time Tomography for Continuously Deforming Objects**](https://vccimaging.org/Publications/Zang2018Space-time/)

[Guangming Zang](https://vccimaging.org/People/zangg/), [Ramzi Idoughi](https://vccimaging.org/People/idoughr/), [Ran Tao](https://cohmas.kaust.edu.sa/Pages/Ran%20Tao.aspx/), [Gilles Lubineau](https://cohmas.kaust.edu.sa/Pages/Gilles%20Lubineau.aspx), [Peter Wonka](http://peterwonka.net/) and [Wolfgang Heidrich](http://vccimaging.org/People/heidriw/)

[King Abdullah University of Science and Technology (KAUST)](https://www.kaust.edu.sa/en)



## Abstract:
X-ray computed tomography is a valuable tool for analyzing objects with
interesting internal structure or complex geometries that are not accessible
with optical means. Unfortunately, tomographic reconstruction of complex
shapes requires a multitude (often hundreds or thousands) of projections
from different viewpoints, which can only be achieved with mechanical
motion. This significantly limits the ability to use x-ray tomography for
either objects that undergo uncontrolled shape change at the time scale of a
scan, or else for analyzing dynamic phenomena, where the motion itself is
under investigation



## Usage
The code is tested in Visual Studio 2015 and 2018 on 64 bits Windows 7/ Windows 10.

To run the ST tomography framwork, you need to install  [argtable], [openmp], and [eigen] libraries

For command usage, you can simpily use  `STTomography --help`  to find more detail:
```
  -s, --Sigma=<double>                         Sigma in volume density reconstruction (default as 0.2)
  -t, --Tau=<double>                           Lambda in volume density reconstruction (default as 0.2)
  -l, --Lambda=<double>                        TV prior weight in volume density reconstruction
  -k, --TemporalPrior=<double>                 Temporal prior weight in volume density reconstruction
  -e, --Bright constancy prior=<double>        Bright constancy weight in volume density reconstruction
  --huber=<double>                             Trade-off weight for huber norm epsilon
  -a, --Alpha=<double>                         Stepsize for SART algorithm  (default as 0.3)
  -o <output>                                  Output 3d image (default as "hello.tif")
  -i, --projsvolume=<file>                     Input Dataset file(.mha.tif.tiff 2d bmp and other raw file are supported)
  -f, --XYZ=<int>                              Takes an integer value (defaults to 9)
  -u, --imgUV=<int>                            w and h of projection image default:512
  -c, --voxelspacing=<double>                  Voxel spacing (default as 1)
  -d, --detector spacing=<double>              Detector spacing (default as 1)
  -g, --sid=<double>                           Source to Iso-Object Distance  (default as 600mm)
  -j, --sdd=<double>                           Source to Detector Distance (default as 1000mm)
  --oxyz=<double>                              Offset x y z of the center of volume (default as 0)
  --sdg=<double>                               Start degree for each proj sequence
  -v, --nframes=<int>                          Number of frames
  -r, --nframes=<int>                          Number of rounds for each frames
  -b, --AlgorithmIter=<int>                    Number of algorithm iterations, default as 20
  -p, --SartIter=<int>                         Number of SART nested iterations, default as 1
  --prior={STV,SAD,ATV}                        Specify the prior you are applying from {STV,SAD,TV} , default as TV
  --bp={Voxelbased,Raybased}                   Specify the backprojection mode {Voxelbased,Raybased}, default as voxelbased        
  -h, --help                                   Print this help and exit
  --version                                    Print version information and exit
```



## License and citation
This research is released under the [CC BY-NC 3.0 US license](https://creativecommons.org/licenses/by-nc/3.0/us/). We encourage an attitude of reproducible work for academic-only purpose. Please kindly cite our work in your publications if it helps your research:

```

@article{Zang2018STTomography,
 author = {Zang, Guangming and Idoughi, Ramzi and Tao, Ran and Lubineau, Gilles and Wonka, Peter and Heidrich, Wolfgang},
 title = {Space-time Tomography for Continuously Deforming Objects},
 journal = {ACM Trans. Graph.},
 issue_date = {August 2018},
 volume = {37},
 number = {4},
 month = jul,
 year = {2018},
 issn = {0730-0301},
 pages = {100:1--100:14},
 articleno = {100},
 numpages = {14},
 url = {http://doi.acm.org/10.1145/3197517.3201298},
 doi = {10.1145/3197517.3201298},
 acmid = {3201298},
 publisher = {ACM},
 address = {New York, NY, USA},
 keywords = {4D reconstruction, X-ray computed tomography, optimization},
} 

```



## Contact: 
If you find any bug or if you have any suggestion or comment, please contact: 

Guangming : guangming.zang@kaust.edu.sa

Copyright: [Visual computing center], CEMSE, KAUST


[openmp]: <http://openmp.org/wp/>
[eigen]: <http://eigen.tuxfamily.org/index.php?title=Main_Page>
[argtable]: <http://argtable.sourceforge.net/>
[previous work]: <https://github.com/vccimaging/TRex-astra-1.7.1beta>
[Visual computing center]: <https://vcc.kaust.edu.sa/Pages/Home.aspx>

Enjoy!




