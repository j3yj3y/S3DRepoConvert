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
    unsigned int v0, v1, v2;
};

struct mesh
{
    std::vector<vertex> vertices;
    std::vector<face> indices;
};



struct camera {
    glm::vec3 pos = glm::vec3(0.f, 0.f, 0.f);
    glm::quat rot = glm::quat_identity<glm::quat::value_type, glm::highp>();

    inline friend std::istream &operator>>(std::istream &iss, camera &cam) {
        iss >> cam.pos.x >> cam.pos.y >> cam.pos.z
            >> cam.rot.x >> cam.rot.y >> cam.rot.z >> cam.rot.w;
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
    std::string filename;
    glm::vec3 t   = glm::vec3(0.f, 0.f, 0.f);
    glm::quat rot = glm::quat_identity<glm::quat::value_type, glm::highp>();

    inline friend std::istream &operator>>(std::istream &iss, partial_scan &scan) {
        iss >> scan.filename;
        iss >> scan.t.x >> scan.t.y >> scan.t.z
            >> scan.rot.x >> scan.rot.y >> scan.rot.z >> scan.rot.w;
        return iss;
    }

    inline glm::mat4x4 getTransform() const {
        return ( glm::translate(t)*glm::mat4_cast(glm::inverse(rot)) );
    }

    void print() const {
        std::cout << filename << "\t"
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
        printf("compound: ");
        compound.print();
        for(const auto& s : single_scans) {
            printf("scan: ");
            s.print();
        }
        std::cout << std::endl << std::flush;
    }
};



conf_file parse_conf_file( const std::filesystem::path & conf_file)
{
    struct conf_file conf;

    std::fstream ss(conf_file, std::ios::in | std::ios::binary);
    if (ss.fail()) throw std::runtime_error("failed to open " + conf_file.string());
    std::string line;

    while (std::getline(ss, line))
    {
        std::istringstream iss(line);
        std::string type;
        if (!(iss >> type)) { break; }

        if (type == "camera")
            iss >> conf.cam;

        if (type == "mesh")
            iss >> conf.compound;

        if (type == "bmesh") {
            partial_scan scan;
            iss >> scan;
            conf.single_scans.push_back(scan);
        }
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
            result_mesh->indices[i] = {static_cast<unsigned int>(iptr[3*i]),
                                       static_cast<unsigned int>(iptr[3*i+1]),
                                       static_cast<unsigned int>(iptr[3*i+2])};

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

        file.write(outstream, false);

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


void compute_scale_values( std::shared_ptr<mesh>& m )
{
    struct edge { unsigned int begin, end; };
    std::vector<edge> edges(m->indices.size()*3);

    // create edges
    for (auto i(0u); i<m->indices.size(); ++i) {
        const auto f = m->indices[i];
        edges[3*i+0] = ((f.v0 < f.v1) ? edge{f.v0, f.v1} : edge{f.v1, f.v0});
        edges[3*i+1] = ((f.v1 < f.v2) ? edge{f.v1, f.v2} : edge{f.v2, f.v1});
        edges[3*i+2] = ((f.v0 < f.v2) ? edge{f.v0, f.v2} : edge{f.v2, f.v0});
    }

    // sort edges
    std::stable_sort(edges.begin(), edges.end(), [](const edge& a, const edge& b){return (a.end < b.end);});
    std::stable_sort(edges.begin(), edges.end(), [](const edge& a, const edge& b){return (a.begin < b.begin);});

    // remove duplicates
    std::vector<edge> sorted_edges;
    auto isDuplicate = [](const edge& a, const edge& b){ return ((a.end==b.end) && (a.begin==b.begin)); };
    for (auto i(0u); i<edges.size(); ++i)
        if (i==0 || !isDuplicate(edges[i-1], edges[i]))
            sorted_edges.push_back(edges[i]);
    edges.clear();

    // accumulate edge-lengths
    struct vertexscale { double s; int adjacency; };
    std::vector<vertexscale> scales(m->vertices.size(), {0., 0});
    for (const auto& e : sorted_edges) {
        const auto v0 = m->vertices[e.begin];
        const auto v1 = m->vertices[e.end];
        const auto length = glm::length( v1.p-v0.p );

        scales[e.begin].s           += static_cast<double>(length);
        scales[e.begin].adjacency   += 1;

        scales[e.end].s         += static_cast<double>(length);
        scales[e.end].adjacency += 1;
    }

    // compute avg length and store scale value in mesh vertex
    for (auto i(0u); i<scales.size(); ++i)
        m->vertices[i].scale = static_cast<float>(scales[i].s / double(scales[i].adjacency));
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
        transform_points(part->vertices, scan.getTransform());
            compute_scale_values(part);
        compound.insert(compound.end(), part->vertices.begin(), part->vertices.end());
    }

    // transform points as given in compound
    //transform_points(compound, conf.compound.getTransform());

    write_ply_file(output, compound);
    return true;
}





#if 0
const float s = 2.f / (rot.x*rot.x + rot.y*rot.y + rot.z*rot.z + rot.w*rot.w);

const float xs = rot.x * s;
const float ys = rot.y * s;
const float zs = rot.z * s;

const float wx = rot.w * xs;
const float wy = rot.w * ys;
const float wz = rot.w * zs;

const float xx = rot.x * xs;
const float xy = rot.x * ys;
const float xz = rot.x * zs;

const float yy = rot.y * ys;
const float yz = rot.y * zs;
const float zz = rot.z * zs;

glm::mat4 trafo;
trafo[0][0] = 1.f - (yy + zz);
trafo[0][1] = xy - wz;
trafo[0][2] = xz + wy;
trafo[0][3] = 0.f;

trafo[1][0] = xy + wz;
trafo[1][1] = 1.f - (xx + zz);
trafo[1][2] = yz - wx;
trafo[1][3] = 0.f;

trafo[2][0] = xz - wy;
trafo[2][1] = yz + wx;
trafo[2][2] = 1.f - (xx + yy);
trafo[2][3] = 0.f;

trafo[3][0] = t.x;
trafo[3][1] = t.y;
trafo[3][2] = t.z;
trafo[3][3] = 1.f;

return trafo;
#endif





} // end namespace range_scan
} // end namespace stanford_repo
