# S3DRepoConvert
Stanford 3D Scanning Repository Config file parser and ply exporter


### Prerequisites
- Meshlab (currently you need meshlab server to triangulate the range scans)
- clone glm header library into the ext folder
- a c++17 capable compiler

### Usage
'''
./Stanford3dRepoConvert -i <parent folder where config file is> -o <folder where to store the results>
'''
Subfolders will be iterated recursively - all found config files will be used


## Compatibility
Currently only tested on linux with gcc 9.2

## Authors
* **Johannes Jakob**

## License
This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments
* Hat tip to: *https://github.com/ddiakopoulos/tinyply* for tinyply implementation!
