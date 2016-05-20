# Prior-less Compressible Structure from Motion
## Information
This is the release of the Prior-less Compressible Structure from Motion, described in the paper http://ci2cv.net/media/papers/Chen-CVPR-2016.pdf.
Project webpage: http://www.cs.cmu.edu/~chenk/projects/cvpr\_2016.html
Code git: git@github.com:kongchen1992/compressible-sfm.git
For questions concerning the code please contact Chen Kong at <chenk@cs.cmu.edu> and Simon Lucey <slucey@cs.cmu.edu>.
The code was developed when Chen Kong and Simon Lucey were at Carnegie Mellon University.

## Cite
If you use this code in your research, please cite:
@article{KongCVPR16, 
  title = {Prior-Less Compressible Structure from Motion},
  author = {Kong, Chen and Lucey, Simon},
  booktitle={Computer Vision and Pattern Recognition (CVPR)},
  year = {2016},
  organization={IEEE}
}

## Folders
`CamRecovery/` contains the necessary code for recovering camera motions and 3D structures.

`Exps/` contains the necessary code for various experiments mentioned in the paper.

`FacMethod/` contains the necessary code for block sparse dictionary learning including block KSVD, block OMP, and block FOCUSS.

`Viz/` contains the necessary code for visualization.
