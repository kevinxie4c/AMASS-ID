## Dependencies

[Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page)

[nlohmann::json](https://github.com/nlohmann/json)

[DART](https://dartsim.github.io/)


## Build

The default building process uses python library in the current environment. See `CMakeLists.txt` if you want to use system python. To build the binaries, do
```
cmake .
make
```
This will generate `inverse.out`, `CharacterViewer.out` in the current working directory.

## Usage

To perform inverse dynamics on the pose file `data/0005_Walking001_poses.npz` (with the default parameters), do
```
perl inverse.pl data/0005_Walking001_poses.npz
```
The result joint torques and the contact forces will be saved in `forces.txt` and `contact_forces.txt` under the output directory (`output/` by default).
You are specify the parameters in the command line options:
```
usage: inverse.pl [options] pose_file

options:
    --char_file=string
    --frame_time=double
    --cutoff_freq=double
    --ground_offset=double
    --outdir=string

```
where `char_file` is the JSON file defining the character skeleton and properties (mass, moment of inertia, etc.), `frame_time` is the duration per frame, `cutoff_freq` is the cutoff frequency of the lowpass filter, `ground_offset` is the height of the ground, and `outdir` is the directory for storing the result.

For example, to perform inverse dynamics on `data/0005_Walking001_poses.npz` with cutoff frequency of 5 Hz, do
```
perl --cutoff_freq=2 inverse.pl data/0005_Walking001_poses.npz
```
The script actually executes `find-contact-node.py`
```
usage: find-contact-node.py filename.npz [outdir]
```
and `inverse.out`
```

options:
-j --char_file=string
-f --frame_time=double
-c --cutoff_freq=double
-g --ground_offset=double
-o --outdir=string
```
You can also run them seperately.

To view the result contact force, run `CharacterViewer.out`. For example,
```
./CharacterViewer.out data/0005_2FeetJump001_poses.npz
```
You can also speify the character definition with option `-j, --char_file`, and the directory where the program read the result contact forces with option `-o, --oudir`:
```
usage: ./CharacterViewer.out [options] pose_file

options:
-j --char_file=string
-o --outdir=string
```
For example,
```
./CharacterViewer.out -j data/character.json -o output
```
