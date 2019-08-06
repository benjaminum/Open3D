#pragma GCC push_options
#pragma GCC optimize ("O0")
// ----------------------------------------------------------------------------
// -                        Open3D: www.open3d.org                            -
// ----------------------------------------------------------------------------
// The MIT License (MIT)
//
// Copyright (c) 2019 www.open3d.org
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.
// ----------------------------------------------------------------------------

#include "Open3D/Geometry/TetraMesh.h"
#include "Open3D/Geometry/TriangleMesh.h"
//#include "Open3D/Geometry/IntersectionTest.h"
//#include "Open3D/Geometry/KDTreeFlann.h"
#include "Open3D/Geometry/PointCloud.h"
#include "Open3D/Geometry/Qhull.h"
#include "Open3D/Geometry/earcut.hpp"

#include <Eigen/Dense>
#include <numeric>
#include <set>
#include <queue>
#include <random>
#include <tuple>
#include <atomic>

#include "Open3D/Utility/Console.h"
#include "Open3D/Utility/DebugHelper.h"
#include <iostream> // TODO rm

namespace open3d {
namespace geometry {

TetraMesh &TetraMesh::Clear() {
    vertices_.clear();
    tetras_.clear();
    return *this;
}

bool TetraMesh::IsEmpty() const { return !HasVertices(); }

Eigen::Vector3d TetraMesh::GetMinBound() const {
    if (!HasVertices()) {
        return Eigen::Vector3d(0.0, 0.0, 0.0);
    }
    return std::accumulate(
            vertices_.begin(), vertices_.end(), vertices_[0],
            [](const Eigen::Vector3d &a, const Eigen::Vector3d &b) {
                return a.array().min(b.array()).matrix();
            });
}

Eigen::Vector3d TetraMesh::GetMaxBound() const {
    if (!HasVertices()) {
        return Eigen::Vector3d(0.0, 0.0, 0.0);
    }
    return std::accumulate(
            vertices_.begin(), vertices_.end(), vertices_[0],
            [](const Eigen::Vector3d &a, const Eigen::Vector3d &b) {
                return a.array().max(b.array()).matrix();
            });
}

TetraMesh &TetraMesh::Transform(const Eigen::Matrix4d &transformation) {
    for (auto &vertex : vertices_) {
        Eigen::Vector4d new_point =
                transformation *
                Eigen::Vector4d(vertex(0), vertex(1), vertex(2), 1.0);
        vertex = new_point.block<3, 1>(0, 0);
    }
    return *this;
}

TetraMesh &TetraMesh::Translate(const Eigen::Vector3d &translation) {
    for (auto &vertex : vertices_) {
        vertex += translation;
    }
    return *this;
}

TetraMesh &TetraMesh::Scale(const double scale, bool center) {
    Eigen::Vector3d vertex_center(0, 0, 0);
    if (center && !vertices_.empty()) {
        vertex_center = std::accumulate(vertices_.begin(), vertices_.end(),
                                        vertex_center);
        vertex_center /= double(vertices_.size());
    }
    for (auto &vertex : vertices_) {
        vertex = (vertex - vertex_center) * scale + vertex_center;
    }
    return *this;
}

TetraMesh &TetraMesh::Rotate(const Eigen::Vector3d &rotation,
                                   bool center,
                                   RotationType type) {
    Eigen::Vector3d vertex_center(0, 0, 0);
    if (center && !vertices_.empty()) {
        vertex_center = std::accumulate(vertices_.begin(), vertices_.end(),
                                        vertex_center);
        vertex_center /= double(vertices_.size());
    }
    const Eigen::Matrix3d R = GetRotationMatrix(rotation, type);
    for (auto &vertex : vertices_) {
        vertex = R * (vertex - vertex_center) + vertex_center;
    }
    return *this;
}

TetraMesh &TetraMesh::operator+=(const TetraMesh &mesh) {
    if (mesh.IsEmpty()) return (*this);
    size_t old_vert_num = vertices_.size();
    size_t add_vert_num = mesh.vertices_.size();
    size_t new_vert_num = old_vert_num + add_vert_num;
    size_t old_tetra_num = tetras_.size();
    size_t add_tetra_num = mesh.tetras_.size();
    vertices_.resize(new_vert_num);
    for (size_t i = 0; i < add_vert_num; i++)
        vertices_[old_vert_num + i] = mesh.vertices_[i];

    tetras_.resize(tetras_.size() + mesh.tetras_.size());
    Eigen::Vector4i64 index_shift = Eigen::Vector4i64::Constant((int64_t)old_vert_num);
    for (size_t i = 0; i < add_tetra_num; i++) {
        tetras_[old_tetra_num + i] = mesh.tetras_[i] + index_shift;
    }
    return (*this);
}

TetraMesh TetraMesh::operator+(const TetraMesh &mesh) const {
    return (TetraMesh(*this) += mesh);
}


TetraMesh &TetraMesh::RemoveDuplicatedVertices() {
    typedef decltype(tetras_)::value_type::Scalar Index;
    typedef std::tuple<double, double, double> Coordinate3;
    std::unordered_map<Coordinate3, size_t,
                       utility::hash_tuple::hash<Coordinate3>>
            point_to_old_index;
    std::vector<Index> index_old_to_new(vertices_.size());
    size_t old_vertex_num = vertices_.size();
    size_t k = 0;                                  // new index
    for (size_t i = 0; i < old_vertex_num; i++) {  // old index
        Coordinate3 coord = std::make_tuple(vertices_[i](0), vertices_[i](1),
                                            vertices_[i](2));
        if (point_to_old_index.find(coord) == point_to_old_index.end()) {
            point_to_old_index[coord] = i;
            vertices_[k] = vertices_[i];
            index_old_to_new[i] = (Index)k;
            k++;
        } else {
            index_old_to_new[i] = index_old_to_new[point_to_old_index[coord]];
        }
    }
    vertices_.resize(k);
    if (k < old_vertex_num) {
        for (auto &tetra : tetras_) {
            tetra(0) = index_old_to_new[tetra(0)];
            tetra(1) = index_old_to_new[tetra(1)];
            tetra(2) = index_old_to_new[tetra(2)];
            tetra(3) = index_old_to_new[tetra(3)];
        }
    }
    utility::LogDebug(
            "[RemoveDuplicatedVertices] {:d} vertices have been removed.\n",
            (int)(old_vertex_num - k));

    return *this;
}

TetraMesh &TetraMesh::RemoveDuplicatedTetras() {
    typedef decltype(tetras_)::value_type::Scalar Index;
    typedef std::tuple<Index, Index, Index, Index> Index4;
    std::unordered_map<Index4, size_t, utility::hash_tuple::hash<Index4>>
            tetra_to_old_index;
    size_t old_tetra_num = tetras_.size();
    size_t k = 0;
    for (size_t i = 0; i < old_tetra_num; i++) {
        Index4 index;
        std::array<Index,4> t{ 
          tetras_[i](0), tetras_[i](1), tetras_[i](2), tetras_[i](3) };

        // We sort the indices to find dubplicates, because tetra (0-1-2-3)
        // and tetra (2-0-3-1) are the same.
        std::sort(t.begin(), t.end());
        index = std::make_tuple(t[0], t[1], t[2], t[3]);

        if (tetra_to_old_index.find(index) == tetra_to_old_index.end()) {
            tetra_to_old_index[index] = i;
            tetras_[k] = tetras_[i];
            k++;
        }
    }
    tetras_.resize(k);
    utility::LogDebug(
            "[RemoveDuplicatedTetras] {:d} tetras have been removed.\n",
            (int)(old_tetra_num - k));

    return *this;
}

TetraMesh &TetraMesh::RemoveUnreferencedVertices() {
    typedef decltype(tetras_)::value_type::Scalar Index;
    std::vector<bool> vertex_has_reference(vertices_.size(), false);
    for (const auto &tetra : tetras_) {
        vertex_has_reference[tetra(0)] = true;
        vertex_has_reference[tetra(1)] = true;
        vertex_has_reference[tetra(2)] = true;
        vertex_has_reference[tetra(3)] = true;
    }
    std::vector<Index> index_old_to_new(vertices_.size());
    size_t old_vertex_num = vertices_.size();
    size_t k = 0;                                  // new index
    for (size_t i = 0; i < old_vertex_num; i++) {  // old index
        if (vertex_has_reference[i]) {
            vertices_[k] = vertices_[i];
            index_old_to_new[i] = (Index)k;
            k++;
        } else {
            index_old_to_new[i] = -1;
        }
    }
    vertices_.resize(k);
    if (k < old_vertex_num) {
        for (auto &tetra : tetras_) {
            tetra(0) = index_old_to_new[tetra(0)];
            tetra(1) = index_old_to_new[tetra(1)];
            tetra(2) = index_old_to_new[tetra(2)];
            tetra(3) = index_old_to_new[tetra(3)];
        }
    }
    utility::LogDebug(
            "[RemoveUnreferencedVertices] {:d} vertices have been removed.\n",
            (int)(old_vertex_num - k));

    return *this;
}

TetraMesh &TetraMesh::RemoveDegenerateTetras() {
    size_t old_tetra_num = tetras_.size();
    size_t k = 0;
    for (size_t i = 0; i < old_tetra_num; i++) {
        const auto &tetra = tetras_[i];
        if (tetra(0) != tetra(1) && 
            tetra(0) != tetra(2) &&
            tetra(0) != tetra(3) &&
            tetra(1) != tetra(2) &&
            tetra(1) != tetra(3) &&
            tetra(2) != tetra(3)) {
            tetras_[k] = tetras_[i];
            k++;
        }
    }
    tetras_.resize(k);
    utility::LogDebug(
            "[RemoveDegenerateTetras] {:d} tetras have been "
            "removed.\n",
            (int)(old_tetra_num - k));
    return *this;
}


std::shared_ptr<TriangleMesh> TetraMesh::ExtractTriangleMesh(const std::vector<double>& values, double level) {
    typedef decltype(tetras_)::value_type::Scalar Index;
    static_assert(std::is_signed<Index>(), "Index type must be signed");
    typedef std::tuple<Index,Index> Index2;
    typedef std::tuple<Index,Index,Index> Index3;

    struct Face {
        Face(): tetra_count(0) { }
        Index3 vertices;
        std::array<size_t,2> tetras;
        uint8_t tetra_count;

        void add_tetra(size_t t){
          if( tetra_count < 2 )
            tetras[tetra_count] = t;
          ++tetra_count;
        }

        bool operator<(const Face& other) const
        {
          return vertices < other.vertices;
        }
        bool operator==(const Face& other) const
        {
          return vertices == other.vertices;
        }
    };

    struct Edge {
        Edge(): adjacent_tetras_count(0) {}
        int adjacent_tetras_count;
        // unique tetra indices
        std::vector<Face> faces;

        void add_or_update_Face( const Face& face )
        {
          auto it = std::find(faces.begin(),faces.end(),face);
          if( it != faces.end() )
            it->add_tetra(face.tetras[0]); // face already exists -> add tetra idx
          else
            faces.push_back(face);
        }
    };


    //sendPolyData("test", 0, vertices_, {3,4}, {0,1,2});


    auto surface_intersection_test = []( double v0, double v1, double level ){
      return (v0 < level && v1 >= level) || (v0 >= level && v1 < level);
    };

    auto make_sorted_tuple2 = []( Index a, Index b ){
      if ( a < b ) return std::make_tuple(a,b);
      else return std::make_tuple(b,a);
    };

    auto make_sorted_tuple3 = []( Index a, Index b, Index c ){
      std::array<Index,3> arr({a,b,c});
      if( arr[0] > arr[1] ) std::swap(arr[0], arr[1]);
      if( arr[0] > arr[2] ) std::swap(arr[0], arr[2]);
      if( arr[1] > arr[2] ) std::swap(arr[1], arr[2]);
      return std::make_tuple(arr[0], arr[1], arr[2]);
    };

    auto create_poly = []( std::set<Face>& faces ){
      std::vector<size_t> poly;
      if( faces.size() < 3 )
        return poly;
      poly.push_back(faces.begin()->tetras[0]);
      size_t next_tetra = faces.begin()->tetras[1];
      faces.erase(faces.begin());

      bool success = true;
      while( faces.size() && success )
      {
        success = false;
        for( const Face& f : faces )
        {
          if( f.tetras[0] == next_tetra || f.tetras[1] == next_tetra )
          {
            // all entries in poly must be unique
            if( std::find(poly.begin(), poly.end(), next_tetra) != poly.end() )
              return std::vector<size_t>();
            poly.push_back(next_tetra);
            next_tetra = f.tetras[0] == next_tetra ? f.tetras[1] : f.tetras[0];
            faces.erase(f);
            success = true;
            break;
          }
        }
      } 
      // faces should be empty if the poly is valid
      if( faces.size() )
        return std::vector<size_t>();
      // the polygon should be closed
      if( poly[0] != next_tetra )
        return std::vector<size_t>();

      return poly;
    };
    
    std::unordered_map<Index2, Edge, utility::hash_tuple::hash<Index2>> intersecting_edges;


    const int tetra_edges[][2] = {{0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3}};
#undef _OPENMP
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for( size_t tetra_i = 0; tetra_i < tetras_.size(); ++tetra_i )
    {
      const auto& tetra = tetras_[tetra_i];
      for( int tet_edge_i = 0; tet_edge_i < 6; ++tet_edge_i )
      {
        Index edge_vert1 = tetra[tetra_edges[tet_edge_i][0]];
        Index edge_vert2 = tetra[tetra_edges[tet_edge_i][1]];
        double vert_value1 = values[edge_vert1];
        double vert_value2 = values[edge_vert2];
        if( surface_intersection_test(vert_value1, vert_value2, level) )
        {
          Index2 index = make_sorted_tuple2(edge_vert1, edge_vert2);
          
          std::array<Face,2> tmp_faces;
          int faces_count = 0;
          for( int i = 0; i < 4; ++i )
          {
            if( i != tetra_edges[tet_edge_i][0] && i != tetra_edges[tet_edge_i][1] )
            {
              // add potential new face
              Face& f = tmp_faces[faces_count];
              f.vertices = make_sorted_tuple3(tetra[i], edge_vert1, edge_vert2);
              f.add_tetra(tetra_i);
              ++faces_count;
            }
          }
#ifdef _OPENMP
#pragma omp critical (intersecting_edges)
#endif
          {
            Edge& e = intersecting_edges[index];
            ++e.adjacent_tetras_count;
            for( int i = 0; i < faces_count; ++i )
              e.add_or_update_Face(tmp_faces[i]);
          } // critical
        }
      }
    }

    // maps the tetrahedron index to the vertex index of the resulting triangle mesh
    std::unordered_map<size_t, size_t> tetra_vertex_index_map;

    // count the number of triangles we have to create
    std::atomic<size_t> triangle_count(0);

    // count the number of vertices we need to generate and do some checks.
    // The checks are not necessary if the tetra mesh is sane.
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for( size_t i = 0; i < intersecting_edges.size(); ++i )
    {
      auto iter = intersecting_edges.begin();
      advance(iter, i);

      Edge& edge = iter->second;

      if( edge.faces.size() < 3 )
      {
        edge.faces.clear();
        edge.adjacent_tetras_count = -1; // flag as ignore
        continue;
      }

      // check if all faces have two tetras and if the ids are ok
      bool face_tetras_ok = true;
      for( const Face& f : edge.faces )
      {
        if( f.tetra_count != 2 )
        {
          face_tetras_ok = false;
          break;
        }
        if( f.tetras[0] == f.tetras[1] )
        {
          face_tetras_ok = false;
          break;
        }
      }
      if( !face_tetras_ok )
      {
        edge.faces.clear();
        edge.adjacent_tetras_count = -1; // flag as ignore
        continue;
      }

      // try to create a valid polygon

      std::set<Face> tmp_faces_set(edge.faces.begin(), edge.faces.end());
      std::vector<size_t> poly = create_poly(tmp_faces_set);
      if( poly.empty() )
      {
        edge.faces.clear();
        edge.adjacent_tetras_count = -1; // flag as ignore
        continue;
      }

      if( poly.size() == 3 )
      {
        ++triangle_count;
      }
      else if( poly.size() == 4 ) // quad
      {
        triangle_count += 2;
      }
      else if( poly.size() >= 5 ) // poly
      {
        // simple polygons (no holes, no self intersections)
        // have N-2 triangles
        triangle_count += poly.size()-2;
      }
      else // this should not happen
      {
        std::cerr << poly.size() << " this should not happen: poly with less than 3 verts\n";
        edge.faces.clear();
        edge.adjacent_tetras_count = -1; // flag as ignore
        continue;
      }


#ifdef _OPENMP
#pragma omp critical(tetra_vertex_index_map)
#endif
      {
        for( size_t tetra_idx : poly )
        {
          if( !tetra_vertex_index_map.count(tetra_idx) )
          {
            size_t index = tetra_vertex_index_map.size();
            tetra_vertex_index_map[tetra_idx] = index;
          }
        }
      } // critical

    }

    const size_t total_vertex_count = tetra_vertex_index_map.size();
    std::vector<int> vertex_write_count(total_vertex_count, 0);


    auto triangle_mesh = std::make_shared<TriangleMesh>();
    triangle_mesh->vertices_.resize(total_vertex_count, Eigen::Vector3d(-1,-1,-1));
    triangle_mesh->triangles_.resize(triangle_count, Eigen::Vector3i(-1,-1,-1));

    std::atomic<size_t> current_triangle_index(0);

    auto compute_edge_vertex = [&values, level, this](Index idx1, Index idx2) {
      double v1 = values[idx1];
      double v2 = values[idx2];
      
      double t = (level-v2)/(v1-v2);
      if( std::isnan(t) || t < 0 || t > 1 )
      {
        std::cerr << "compute_edge_vertex!\n";
        t = 0.5;
      }
      Eigen::Vector3d intersection = t*vertices_[idx1] + (1-t)*vertices_[idx2];
      
      return intersection;
    };

    auto compute_tetra_vertex = [&]( size_t idx ){

      const auto& tetra = tetras_[idx];
      
      Eigen::Vector3d vert(0,0,0);
      double sum = 0;

      for( int i = 0; i < 6; ++i )
      {
        Index edge_vert1 = tetra[tetra_edges[i][0]];
        Index edge_vert2 = tetra[tetra_edges[i][1]];
        double vert_value1 = values[edge_vert1];
        double vert_value2 = values[edge_vert2];
        if( surface_intersection_test(vert_value1, vert_value2, level) )
        {
          vert += compute_edge_vertex(edge_vert1, edge_vert2);
          sum += 1;
        }
      }
      vert /= sum;
      return vert;
    };

    auto compute_triangle_normal = [](
        const Eigen::Vector3d& a, 
        const Eigen::Vector3d& b, 
        const Eigen::Vector3d& c
        ){
      Eigen::Vector3d ab(b-a);
      Eigen::Vector3d ac(c-a);
      return ab.cross(ac);
    };
    
    auto compute_points_eigenvectors = [](const std::vector<Eigen::Vector3d>& points)
    {
      Eigen::Matrix3d covariance;
      Eigen::Matrix<double, 9, 1> cumulants;
      cumulants.setZero();
      for (size_t i = 0; i < points.size(); i++) {
          const Eigen::Vector3d &point = points[i];
          cumulants(0) += point(0);
          cumulants(1) += point(1);
          cumulants(2) += point(2);
          cumulants(3) += point(0) * point(0);
          cumulants(4) += point(0) * point(1);
          cumulants(5) += point(0) * point(2);
          cumulants(6) += point(1) * point(1);
          cumulants(7) += point(1) * point(2);
          cumulants(8) += point(2) * point(2);
      }
      cumulants /= (double)points.size();
      covariance(0, 0) = cumulants(3) - cumulants(0) * cumulants(0);
      covariance(1, 1) = cumulants(6) - cumulants(1) * cumulants(1);
      covariance(2, 2) = cumulants(8) - cumulants(2) * cumulants(2);
      covariance(0, 1) = cumulants(4) - cumulants(0) * cumulants(1);
      covariance(1, 0) = covariance(0, 1);
      covariance(0, 2) = cumulants(5) - cumulants(0) * cumulants(2);
      covariance(2, 0) = covariance(0, 2);
      covariance(1, 2) = cumulants(7) - cumulants(1) * cumulants(2);
      covariance(2, 1) = covariance(1, 2);

      Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver;
      solver.compute(covariance, Eigen::ComputeEigenvectors);
      // store normal to intermediate var for compiler
      Eigen::Matrix3d ev = solver.eigenvectors();
      return ev;
    };
    
    auto create_basis_vectors = [](const Eigen::Vector3d& n) {
      // adapted from Duff et al. "Building an Orthonormal Basis, Revisited".
      double sign = std::copysignf(1.0f, n.z());
      const double a = -1.0f / (sign + n.z());
      const double b = n.x() * n.y() * a;
      Eigen::Vector3d b1 = Eigen::Vector3d(1.0f + sign * n.x() * n.x() * a, sign * b, -sign * n.x());
      Eigen::Vector3d b2 = Eigen::Vector3d(b, sign + n.y() * n.y() * a, -n.y());
      return std::make_tuple(b1,b2);
    };


#ifdef _OPENMP
#pragma omp parallel for
#endif
    for( size_t i = 0; i < intersecting_edges.size(); ++i )
    {
      auto iter = intersecting_edges.begin();
      advance(iter, i);

      const Edge& edge = iter->second;
      // skip edges for which we cannot create a polygon
      if( edge.adjacent_tetras_count == -1 )
        continue;

      Index edge_vert1, edge_vert2;
      std::tie(edge_vert1, edge_vert2) = iter->first;

      // make edge_vert1 be the vertex that is smaller than level (inside)
      if( values[edge_vert1] > values[edge_vert2] )
        std::swap(edge_vert1, edge_vert2);


      // create a valid polygon
      std::set<Face> tmp_faces_set(edge.faces.begin(), edge.faces.end());
      std::vector<size_t> poly = create_poly(tmp_faces_set);
      std::vector<Eigen::Vector3d> verts(poly.size());

      // compute the vertices within the tetrahedra
      for( size_t i = 0; i < poly.size(); ++i )
      {
        verts[i] = compute_tetra_vertex(poly[i]);
      }
#ifdef _OPENMP
#pragma omp critical(vertices)
#endif
      {
        for( size_t i = 0; i < poly.size(); ++i )
        {
          triangle_mesh->vertices_[tetra_vertex_index_map[poly[i]]] = verts[i];
          vertex_write_count.at(tetra_vertex_index_map[poly[i]]) += 1;
        }
      } // critical


      Eigen::Vector3d edge_dir(vertices_[edge_vert2]-vertices_[edge_vert1]);


      if( poly.size() == 3 )
      {
        // check triangle orientation
        double dot = edge_dir.dot(compute_triangle_normal(verts[0], verts[1], verts[2]));
        if( dot < 0 )
        {
          std::reverse(verts.begin(), verts.end());
          std::reverse(poly.begin(), poly.end());
        }
        size_t triangle_index = current_triangle_index.fetch_add(1);
        Eigen::Vector3i tri(tetra_vertex_index_map[poly[0]], tetra_vertex_index_map[poly[1]], tetra_vertex_index_map[poly[2]]);;
        triangle_mesh->triangles_[triangle_index] = tri;

        {
          std::vector<Eigen::Vector3d> tmp(verts); 
          std::vector<int> tris;
          tris.insert(tris.end(), {0,1,2});
          //sendPolyData("tris", triangle_index, tmp, {}, tris);
        }

      }
      else if( poly.size() == 4 )
      {
        // check triangle orientation
        double dot = edge_dir.dot(compute_triangle_normal(verts[0], verts[1], verts[2]));
        dot += edge_dir.dot(compute_triangle_normal(verts[2], verts[3], verts[0]));
        if( dot < 0 )
        {
          std::reverse(verts.begin(), verts.end());
          std::reverse(poly.begin(), poly.end());
        }

        size_t triangle_index = current_triangle_index.fetch_add(2);
        // select shorter diaginal
        if( (verts[0]-verts[2]).squaredNorm() < (verts[1]-verts[3]).squaredNorm() )
        {
          triangle_mesh->triangles_[triangle_index] << tetra_vertex_index_map[poly[0]], tetra_vertex_index_map[poly[1]], tetra_vertex_index_map[poly[2]];
          triangle_mesh->triangles_[triangle_index+1] << tetra_vertex_index_map[poly[2]], tetra_vertex_index_map[poly[3]], tetra_vertex_index_map[poly[0]];
        }
        else
        {
          triangle_mesh->triangles_[triangle_index] << tetra_vertex_index_map[poly[0]], tetra_vertex_index_map[poly[1]], tetra_vertex_index_map[poly[3]];
          triangle_mesh->triangles_[triangle_index+1] << tetra_vertex_index_map[poly[1]], tetra_vertex_index_map[poly[2]], tetra_vertex_index_map[poly[3]];
        }

        {
          std::vector<Eigen::Vector3d> tmp(verts); 
          std::vector<int> tris;
          if( (verts[0]-verts[2]).squaredNorm() < (verts[1]-verts[3]).squaredNorm() )
          {
            tris.insert(tris.end(), {0,1,2});
            tris.insert(tris.end(), {2,3,0});
          }
          else
          {
            tris.insert(tris.end(), {0,1,3});
            tris.insert(tris.end(), {1,2,3});
          }
          //sendPolyData("quads", triangle_index, tmp, {}, tris);
        }

      }
      else if( poly.size() >= 5 )
      {
        // project to a plane and triangulate the 2d poly with earcutting
        bool success = false;
        int attempt = 0;
        std::vector<std::vector<std::array<double,2>>> poly2d(1);
        poly2d[0].resize(poly.size());
        while( attempt < 3 && !success )
        {
          // usually the first ev is the normal but in degenerate cases
          // the other eigenvectors provide a better basis for projecting
          // the polygon
          Eigen::Vector3d poly_normal = compute_points_eigenvectors(verts).col(attempt);
          Eigen::Vector3d b1, b2;
          std::tie(b1,b2) = create_basis_vectors(poly_normal);

          // define a simple polygon without holes
          for( size_t i = 0; i < poly.size(); ++i )
          {
            poly2d[0][i] = {b1.dot(verts[i]-verts[0]), b2.dot(verts[i]-verts[0])};
          }
          std::vector<int> indices = mapbox::earcut<int>(poly2d);
          const int new_triangles = indices.size()/3;

          if( indices.size() % 3 == 0 && new_triangles == int(poly.size()-2) )
          {
            success = true;
            double dot = 0;
            for( int i = 0; i < new_triangles; ++i )
            {
              dot += edge_dir.dot(compute_triangle_normal(verts[indices[3*i+0]],verts[indices[3*i+1]],verts[indices[3*i+2]]));
            }
            if( dot < 0 )
            {
              for( int i = 0; i < new_triangles; ++i )
              {
                std::swap(indices[i*3+1], indices[i*3+2]);
              }
            }
            size_t triangle_index = current_triangle_index.fetch_add(new_triangles);
            for( int i = 0; i < new_triangles; ++i )
            {
              triangle_mesh->triangles_[triangle_index+i] << tetra_vertex_index_map[poly[indices[3*i+0]]], tetra_vertex_index_map[poly[indices[3*i+1]]], tetra_vertex_index_map[poly[indices[3*i+2]]];
            }





        //if( indices.size() != 3*(poly.size()-2) )
        //{
          //{
            //std::vector<Eigen::Vector3d> tmp;
            //tmp.push_back(verts[0]); 
            //tmp.push_back(verts[0]+poly_normal); 

            //std::vector<int> lines = {0,1};
            //sendPolyData("poly_normal", triangle_index, tmp, lines, {});
          //}

          //{
            //std::vector<Eigen::Vector3d> tmp(verts); 

            //std::vector<int> lines;
            //for( size_t i = 0; i < poly.size(); ++i )
            //{
              //lines.insert(lines.end(), {int(i), int((i+1)%poly.size())});
            //}
            //sendPolyData("poly", triangle_index, tmp, lines, {});
          //}

          //{
            //std::vector<Eigen::Vector3d> tmp; 
            //for( const auto& v2 : poly2d[0] ){
              //tmp.push_back(Eigen::Vector3d(v2[0], v2[1],0));
            //}

            //std::vector<int> lines;
            //for( size_t i = 0; i < tmp.size(); ++i )
            //{
              //lines.insert(lines.end(), {int(i), int((i+1)%poly.size())});
            //}
            //sendPolyData("poly2d", triangle_index, tmp, lines, indices);
          //}
        //}






          }
          ++attempt;
        }

        if( !success )
        {
            std::cerr << " this should not happen: cannot triangulate polygon\n";
        }





      }


    }

    triangle_mesh->triangles_.resize(current_triangle_index);
    return triangle_mesh;
}


}  // namespace geometry
}  // namespace open3d
