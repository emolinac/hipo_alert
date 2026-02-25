# hipo (emolinac fork)

This is a software to produce TTree with branches **exclusively** related ALERT's AHDC.  

The original software can be found [here]{https://github.com/gavalian/hipo}.

## Installing the package

Package has a dependency on LZ4 compression library, which is a submodule.
Use the following command to clone the repo and its dependency:

```
 git clone --recurse-submodules git@github.com:emolinac/hipo_alert.git
``` 

Load the respective modules in the ifarm:
```
module use /scigroup/cvmfs/hallb/clas12/sw/modulefiles
module load clas12
```
Then compile with:
```
cd ./hipo_alert/
make

cd examples/root/
make converter_alertbanks
```

## Usage

To convert a hipo file, run the following command:
```
./converter_alertbanks.exe /absolute-path-to-hipo-file/hipo-file.hipo 0
```
To convert at large scale use the script `perform_conversion.sh` which:
1. It will dump the hipo files names into a txt file.
2. The files listed in the txt file are copied into a temporary folder of your choice
3. The conversion happens and then the file in the temporary folder is deleted.