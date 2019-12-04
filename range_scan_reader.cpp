#include "range_scan_reader.h"

#include <iostream>
#include <vector>
#include <cstring>
#include <fstream>
#include <algorithm>
#include <memory>

#include "tinyply.h"

#define GLM_FORCE_SWIZZLE
#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtx/quaternion.hpp>
#include <glm/vector_relational.hpp>



namespace stanford_repo { namespace range_scan {


struct vertex {
    vertex() : p(0,0,0), n(0,0,0), scale(0) {}
    vertex(const float x, const float y, const float z,
           const float nx, const float ny, const float nz,
           const float s)
        : p(x,y,z), n(nx,ny,nz), scale(s) {}


    inline vertex& transform(const glm::mat4x4& trafo, const glm::mat4x4& inv_trafo)
    {
        const glm::vec4 tmp = trafo * glm::vec4(p, 1.f);
        p = tmp.xyz() / tmp.w;

        const glm::vec4 tmpn = inv_trafo * glm::vec4(n, 1.f);
        n = glm::normalize( tmpn.xyz() );

        return *this;
    }

    glm::vec3 p;
    glm::vec3 n;
    float scale;
};

struct face {
    int v0, v1, v2;
};

struct mesh
{
    std::vector<vertex> vertices;
    std::vector<face> indices;
};



struct camera {
    glm::vec3 pos = glm::vec3(0.f, 0.f, 0.f);
    glm::quat rot = glm::quat(1.f, 0.f, 0.f, 0.f);

    inline friend std::istream &operator>>(std::istream &iss, camera &cam) {
        std::string tmp;
        iss >> tmp;
        iss >> cam.pos.x >> cam.pos.y >> cam.pos.z  >> cam.rot.w >> cam.rot.x >> cam.rot.y >> cam.rot.z;
        return iss;
    }

    void print() const {
        std::cout << "camera: "
                  << "pos( " << pos.x << ", " << pos.y << ", " << pos.z << " ) "
                  << "quat( " << rot.x << ", " << rot.y << ", " << rot.z << ", " << rot.w << " )"
                  << std::endl;
    }
};


struct partial_scan {
    std::string component;
    std::string filename;
    glm::vec3 t   = glm::vec3(0.f, 0.f, 0.f);
    glm::quat rot = glm::quat(1.f, 0.f, 0.f, 0.f);

    inline friend std::istream &operator>>(std::istream &iss, partial_scan &scan) {
        iss >> scan.component;
        iss >> scan.filename;
        iss >> scan.t.x >> scan.t.y >> scan.t.z >> scan.rot.w >> scan.rot.x >> scan.rot.y >> scan.rot.z;
        return iss;
    }

    inline glm::mat4x4 getTransform() const {
        return (glm::translate(t)*glm::mat4_cast(rot));
    }

    void print() const {
        std::cout << "partial_scan: "
                  << component << filename << "\t"
                  << "trans( " << t.x << ", " << t.y << ", " << t.z << " ) "
                  << "quat( " << rot.x << ", " << rot.y << ", " << rot.z << ", " << rot.w << " )"
                  << std::endl;
    }
};


struct conf_file {
    camera cam;
    partial_scan compound;
    std::vector<partial_scan> single_scans;

    void print() const {
        std::cout << "conf_file:\n";
        cam.print();
        compound.print();
        for(const auto& s : single_scans) s.print();
        std::cout << std::endl << std::flush;
    }
};



conf_file parse_conf_file( const std::filesystem::path & conf_file)
{
    struct conf_file conf;

    std::fstream ss(conf_file, std::ios::in | std::ios::binary);
    if (ss.fail()) throw std::runtime_error("failed to open " + conf_file.string());
    std::string line;

    {
        // first extract camera info:
        if( !std::getline(ss, line) ) throw std::runtime_error("invalid conf file");
        std::istringstream iss(line);
        iss >> conf.cam;

        // extract compound
        if( !std::getline(ss, line) ) throw std::runtime_error("invalid conf file");
        iss = std::istringstream(line);
        iss >> conf.compound;
    }

    // extract single scans:
    while (std::getline(ss, line))
    {
        std::istringstream iss(line);
        partial_scan scan;
        if (!(iss >> scan)) { break; }

        conf.single_scans.push_back(scan);
    }

    ss.close();

    return conf;
}



std::shared_ptr<mesh> read_ply_file( const std::filesystem::path& file )
{
    using namespace tinyply;

    try
    {
        std::ifstream ss(file, std::ios::binary);
        if (ss.fail()) throw std::runtime_error("failed to open " + file.string());

        PlyFile file;
        file.parse_header(ss);

        std::shared_ptr<PlyData> plyvertices, plyfaces;

        try { plyvertices = file.request_properties_from_element("vertex", { "x", "y", "z", "nx", "ny", "nz" }); }
        catch (const std::exception &) { throw std::invalid_argument("file does not contain vertices"); }

        try { plyfaces = file.request_properties_from_element("face", { "vertex_indices" }, 3); }
        catch (const std::exception &) { throw std::invalid_argument("file does not face indices"); }

        file.read(ss);

        auto result_mesh = std::make_shared<mesh>();
        result_mesh->vertices = std::vector<vertex>(plyvertices->count);
        result_mesh->indices  = std::vector<face>(plyfaces->count);

        const float* const vptr = reinterpret_cast<float*>(plyvertices->buffer.get());
        const int* const iptr   = reinterpret_cast<int*>(plyfaces->buffer.get());

        for (auto i(0u); i<plyvertices->count; ++i)
            result_mesh->vertices[i] = {vptr[6*i], vptr[6*i+1], vptr[6*i+2], vptr[6*i+3], vptr[6*i+4], vptr[6*i+5], 0.f};

        for (auto i(0u); i<plyfaces->count; ++i)
            result_mesh->indices[i] = {iptr[3*i], iptr[3*i+1], iptr[3*i+2]};

        return result_mesh;
    }
    catch (const std::exception & e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        return std::shared_ptr<mesh>();
    }
}



bool write_ply_file( const std::filesystem::path& file,
                     std::vector<vertex>& vertices )
{
    using namespace tinyply;

    try
    {
        std::filebuf fb;
        fb.open(file, std::ios::out | std::ios::binary);
        std::ostream outstream(&fb);
        if(outstream.fail()) throw std::runtime_error("failed to open " + file.string());

        PlyFile file;

        file.add_properties_to_element("vertex", { "x", "y", "z", "nx", "ny", "nz", "scale" },
            Type::FLOAT32, vertices.size(), reinterpret_cast<uint8_t*>(vertices.data()), Type::INVALID, 0);

        file.write(outstream, true);

        return true;
    }
    catch (const std::exception & e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        return false;
    }
}


inline void transform_points( std::vector<vertex>& vertices, const glm::mat4& trafo )
{
    const auto inv_trafo = glm::transpose(glm::inverse(trafo));

    for (auto& v : vertices)
        v.transform(trafo, inv_trafo);
}


bool transform( const std::filesystem::path& conf_file,
                const std::filesystem::path& output )
{
    const std::filesystem::path folder = conf_file.parent_path();
    const struct conf_file      conf   = parse_conf_file(conf_file);
    std::vector<vertex>         compound;


        for (auto & scan : conf.single_scans)
        {
            const std::filesystem::path input_file{ folder.string() + "/" + scan.filename };
            const std::filesystem::path mesh_file{ "./tmp.ply" };

            // execute meshlab server to compute vertex normals and triangulate range scan
            const auto cmd = std::string("meshlabserver -i " + input_file.string() + " -o " + mesh_file.string() + " -m vn");
            if (system(cmd.c_str()) != 0) { printf("ERROR: no meshlabserver available!"); return false; }

            // read mesh, transform points, compute scale from edge lengths and insert points into commpound
            std::shared_ptr<mesh> part = read_ply_file(mesh_file);
            transform_points(part->vertices,  scan.getTransform());
            compound.insert(compound.end(), part->vertices.begin(), part->vertices.end());

            //write_ply_file(input_file.string() + ".normals.ply", part->vertices);
        }

    write_ply_file(output, compound);
}




} // end namespace range_scan
} // end namespace stanford_repo











#if 0

// OTHER SOLUTIONS
boost::math::quaternion<double> translation, quaternionRotation;

//Get Transformation
translationVec = glm::dvec4(lineData[2].toDouble(), lineData[3].toDouble(), lineData[4].toDouble(),0.0);
quaternionRotation = boost::math::quaternion<double>(lineData[8].toDouble(),lineData[5].toDouble(),lineData[6].toDouble(),lineData[7].toDouble());

//calculate the unit quaternion
double magnitude = std::sqrt(
      quaternionRotation.R_component_1()*quaternionRotation.R_component_1()+
      quaternionRotation.R_component_2()*quaternionRotation.R_component_2()+
      quaternionRotation.R_component_3()*quaternionRotation.R_component_3()+
      quaternionRotation.R_component_4()*quaternionRotation.R_component_4());
quaternionRotation /= magnitude;
rotationMat = this->quat_to_mat(quaternionRotation);

//do some file related stuff
//...

//for each line: read the point data and transform it and store the point in a data array
pointData[j].x = stringPointData[0].toDouble();
pointData[j].y = stringPointData[1].toDouble();
pointData[j].z = stringPointData[2].toDouble();
//transform the curren point
glm::dvec4 curPoint =  glm::dvec4(pointData[j].x,pointData[j].y,pointData[j].z,1.0);
//first rotation
curPoint = rotationMat*curPoint;
//then translation
curPoint += translationVec;
//store the data in a data array
pointData[j].x = curPoint.x;
pointData[j].y = curPoint.y;
pointData[j].z = curPoint.z;





// OTHER SOLUTIONS
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265
#endif

class LineInput {
  public:
    LineInput(const std::string& filename) {
        F_ = fopen(filename.c_str(), "r" ) ;
        ok_ = (F_ != 0) ;
    }
    ~LineInput() {
        if(F_ != 0) {
            fclose(F_); F_ = 0 ;
        }
    }
    bool OK() const { return ok_ ; }
    bool eof() const { return feof(F_) ; }
    bool get_line() {
        line_[0] = '\0' ;
        // Skip the empty lines
        while(!isprint(line_[0])) {
            if(fgets(line_, MAX_LINE_LEN, F_) == 0) {
                return false ;
            }
        }
        // If the line ends with a backslash, append
        // the next line to the current line.
        bool check_multiline = true ;
        int total_length = MAX_LINE_LEN ;
        char* ptr = line_ ;
        while(check_multiline) {
            int L = strlen(ptr) ;
            total_length -= L ;
            ptr = ptr + L - 2;
            if(*ptr == '\\' && total_length > 0) {
                *ptr = ' ' ;
                ptr++ ;
                fgets(ptr, total_length, F_) ;
            } else {
                check_multiline = false ;
            }
        }
        if(total_length < 0) {
            std::cerr
                << "MultiLine longer than "
                << MAX_LINE_LEN << " bytes" << std::endl ;
        }
        return true ;
    }
    int nb_fields() const { return field_.size() ;  }
    char* field(int i) { return field_[i] ;   }
    int field_as_int(int i) {
        int result ;
        ok_ = ok_ && (sscanf(field(i), "%d", &result) == 1) ;
        return result ;
    }
    double field_as_double(int i) {
        double result ;
        ok_ = ok_ && (sscanf(field(i), "%lf", &result) == 1) ;
        return result ;
    }
    bool field_matches(int i, const char* s) {
        return !strcmp(field(i), s) ;
    }
    void get_fields(const char* separators=" \t\r\n") {
        field_.resize(0) ;
        char* tok = strtok(line_,separators) ;
        while(tok != 0) {
            field_.push_back(tok) ;
            tok = strtok(0,separators) ;
        }
    }

  private:
    enum { MAX_LINE_LEN = 65535 } ;
    FILE* F_ ;
    char line_[MAX_LINE_LEN] ;
    std::vector<char*> field_ ;
    bool ok_ ;
} ;

std::string to_string(int x, int mindigits) {
    char buff[100] ;
    sprintf(buff, "%03d", x) ;
    return std::string(buff) ;
}

double M[4][4] ;

void transform(double* xyz) {
    double xyzw[4] ;
    for(unsigned int c=0; c<4; c++) {
        xyzw[c] = M[3][c] ;
    }
    for(unsigned int j=0; j<4; j++) {
        for(unsigned int i=0; i<3; i++) {
            xyzw[j] += M[i][j] * xyz[i] ;
        }
    }
    for(unsigned int c=0; c<3; c++) {
        xyz[c] = xyzw[c] / xyzw[3] ;
    }
}

bool read_frames_file(int no) {
    std::string filename = "scan" + to_string(no,3) + ".frames" ;
    std::cerr << "Reading frames from:" << filename << std::endl ;
    LineInput in(filename) ;
    if(!in.OK()) {
       std::cerr << " ... not found" << std::endl ;
       return false ;
    }
    while(!in.eof() && in.get_line()) {
        in.get_fields() ;
        if(in.nb_fields() == 17) {
            int f = 0 ;
            for(unsigned int i=0; i<4; i++) {
                for(unsigned int j=0; j<4; j++) {
                    M[i][j] = in.field_as_double(f) ; f++ ;
                }
            }
        }
    }
    return true ;
}

bool read_pose_file(int no) {
    std::string filename = "scan" + to_string(no,3) + ".pose" ;
    std::cerr << "Reading pose from:" << filename << std::endl ;
    LineInput in(filename) ;
    if(!in.OK()) {
       std::cerr << " ... not found" << std::endl ;
       return false ;
    }
    double xyz[3] ;
    double euler[3] ;
    in.get_line() ;
    in.get_fields() ;
    xyz[0] = in.field_as_double(0) ;
    xyz[1] = in.field_as_double(1) ;
    xyz[2] = in.field_as_double(2) ;
    in.get_line() ;
    in.get_fields() ;
    euler[0] = in.field_as_double(0) * M_PI / 180.0 ;
    euler[1] = in.field_as_double(1) * M_PI / 180.0 ;
    euler[2] = in.field_as_double(2) * M_PI / 180.0 ;

   double sx = sin(euler[0]);
   double cx = cos(euler[0]);
   double sy = sin(euler[1]);
   double cy = cos(euler[1]);
   double sz = sin(euler[2]);
   double cz = cos(euler[2]);

   M[0][0] = cy*cz;
   M[0][1] = sx*sy*cz + cx*sz;
   M[0][2] = -cx*sy*cz + sx*sz;
   M[0][3] = 0.0;
   M[1][0] = -cy*sz;
   M[1][1] = -sx*sy*sz + cx*cz;
   M[1][2] = cx*sy*sz + sx*cz;
   M[1][3] = 0.0;
   M[2][0] = sy;
   M[2][1] = -sx*cy;
   M[2][2] = cx*cy;
   M[2][3] = 0.0;
   M[3][0] = xyz[0];
   M[3][1] = xyz[1];
   M[3][2] = xyz[2];
   M[3][3] = 1.0;
   return true ;
}

void setup_transform_from_translation_and_quaternion(
    double Tx, double Ty, double Tz,
    double Qx, double Qy, double Qz, double Qw
) {
    /* for unit q, just set s = 2 or set xs = Qx + Qx, etc. */

    double s = 2.0 / (Qx*Qx + Qy*Qy + Qz*Qz + Qw*Qw);

    double xs = Qx * s;
    double ys = Qy * s;
    double zs = Qz * s;

    double wx = Qw * xs;
    double wy = Qw * ys;
    double wz = Qw * zs;

    double xx = Qx * xs;
    double xy = Qx * ys;
    double xz = Qx * zs;

    double yy = Qy * ys;
    double yz = Qy * zs;
    double zz = Qz * zs;

    M[0][0] = 1.0 - (yy + zz);
    M[0][1] = xy - wz;
    M[0][2] = xz + wy;
    M[0][3] = 0.0;

    M[1][0] = xy + wz;
    M[1][1] = 1 - (xx + zz);
    M[1][2] = yz - wx;
    M[1][3] = 0.0;

    M[2][0] = xz - wy;
    M[2][1] = yz + wx;
    M[2][2] = 1 - (xx + yy);
    M[2][3] = 0.0;

    M[3][0] = Tx;
    M[3][1] = Ty;
    M[3][2] = Tz;
    M[3][3] = 1.0;
}

bool read_points_file(int no) {
    std::string filename = "scan" + to_string(no,3) + ".3d" ;
    std::cerr << "Reading points from:" << filename << std::endl ;
    LineInput in(filename) ;
    if(!in.OK()) {
       std::cerr << " ... not found" << std::endl ;
       return false ;
    }
    while(!in.eof() && in.get_line()) {
        in.get_fields() ;
        double xyz[3] ;
        if(in.nb_fields() >= 3) {
            for(unsigned int c=0; c<3; c++) {
                xyz[c] = in.field_as_double(c) ;
            }
            transform(xyz) ;
            printf("%f %f %f\n",xyz[0],xyz[1],xyz[2]) ;
        }
    }
    return true ;
}

.0530127 0.138516 0.0990356 0.908911 -0.0569874 0.154429 0.383126

/* only works for ASCII PLY files */
void read_ply_file(char* filename) {
    std::cerr << "Reading points from:" << filename << std::endl;
    LineInput in(filename) ;
    if(!in.OK()) {
        std::cerr << filename << ": could not open" << std::endl ;
        return;
    }
    bool reading_vertices = false;
    int nb_vertices = 0 ;
    int nb_read_vertices = 0 ;
    while(!in.eof() && in.get_line()) {
        in.get_fields();
        if(reading_vertices) {
            double xyz[3] ;
            for(unsigned int c=0; c<3; c++) {
                xyz[c] = in.field_as_double(c) ;
            }
            transform(xyz) ;
            printf("%f %f %f\n",xyz[0],xyz[1],xyz[2]) ;
            ++nb_read_vertices;
            if(nb_read_vertices == nb_vertices) {
                return;
            }
        } else if(
            in.field_matches(0,"element") &&
            in.field_matches(1,"vertex")
        ) {
            nb_vertices = in.field_as_int(2);
        } else if(in.field_matches(0,"end_header")) {
            reading_vertices = true;
        }
    }
}

/* For Stanford scanning repository */
void read_conf_file(char* filename) {
    LineInput in(filename) ;
    if(!in.OK()) {
        std::cerr << filename << ": could not open" << std::endl ;
        return;
    }
    while(!in.eof() && in.get_line()) {
        in.get_fields();
        if(in.nb_fields() == 0) { continue ; }
        if(in.field_matches(0,"bmesh")) {
            char* filename = in.field(1);
            // Translation vector
            double Tx = in.field_as_double(2);
            double Ty = in.field_as_double(3);
            double Tz = in.field_as_double(4);
            /// Quaternion
            double Qx = in.field_as_double(5);
            double Qy = in.field_as_double(6);
            double Qz = in.field_as_double(7);
            double Qw = in.field_as_double(8);
            setup_transform_from_translation_and_quaternion(Tx,Ty,Tz,Qx,Qy,Qz,Qw);
            read_ply_file(filename);
        }
    }
 }


int main(int argc, char** argv) {
    if(argc != 2) { return -1 ; }
    if(strstr(argv[1],".conf")) {
        read_conf_file(argv[1]);
    } else {
        int max_i = atoi(argv[1]) ;
        for(int i=0; i<=max_i; i++) {
            if(!read_frames_file(i)) {
                read_pose_file(i) ;
            }
            read_points_file(i) ;
        }
    }
    return 0 ;
}




#endif

