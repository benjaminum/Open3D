// ----------------------------------------------------------------------------
// -                        Open3D: www.open3d.org                            -
// ----------------------------------------------------------------------------
// The MIT License (MIT)
//
// Copyright (c) 2018 www.open3d.org
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
#include "Open3D/Geometry/PointCloud.h"
#include "TestUtility/UnitTest.h"

using namespace Eigen;
using namespace open3d;
using namespace std;
using namespace unit_test;

// ----------------------------------------------------------------------------
//
// ----------------------------------------------------------------------------
TEST(TetraMesh, Constructor) {
    geometry::TetraMesh tm;

    // inherited from Geometry2D
    EXPECT_EQ(geometry::Geometry::GeometryType::TetraMesh,
              tm.GetGeometryType());
    EXPECT_EQ(3, tm.Dimension());

    // public member variables
    EXPECT_EQ(0u, tm.vertices_.size());
    EXPECT_EQ(0u, tm.tetras_.size());

    // public members
    EXPECT_TRUE(tm.IsEmpty());

    ExpectEQ(Zero3d, tm.GetMinBound());
    ExpectEQ(Zero3d, tm.GetMaxBound());

    EXPECT_FALSE(tm.HasVertices());
    EXPECT_FALSE(tm.HasTetras());
}

// ----------------------------------------------------------------------------
//
// ----------------------------------------------------------------------------
TEST(TetraMesh, DISABLED_MemberData) { unit_test::NotImplemented(); }

// ----------------------------------------------------------------------------
//
// ----------------------------------------------------------------------------
TEST(TetraMesh, Clear) {
    int size = 100;

    Vector3d dmin(0.0, 0.0, 0.0);
    Vector3d dmax(1000.0, 1000.0, 1000.0);

    Vector4i imin(0, 0, 0, 0);
    Vector4i imax(size - 1, size - 1, size - 1, size - 1);

    geometry::TetraMesh tm;

    tm.vertices_.resize(size);
    tm.tetras_.resize(size);

    Rand(tm.vertices_, dmin, dmax, 0);
    Rand(tm.tetras_, imin, imax, 0);

    EXPECT_FALSE(tm.IsEmpty());

    ExpectEQ(Vector3d(19.607843, 0.0, 0.0), tm.GetMinBound());
    ExpectEQ(Vector3d(996.078431, 996.078431, 996.078431), tm.GetMaxBound());

    EXPECT_TRUE(tm.HasVertices());
    EXPECT_TRUE(tm.HasTetras());

    tm.Clear();

    // public members
    EXPECT_TRUE(tm.IsEmpty());

    ExpectEQ(Zero3d, tm.GetMinBound());
    ExpectEQ(Zero3d, tm.GetMaxBound());

    EXPECT_FALSE(tm.HasVertices());
    EXPECT_FALSE(tm.HasTetras());
}

// ----------------------------------------------------------------------------
//
// ----------------------------------------------------------------------------
TEST(TetraMesh, IsEmpty) {
    int size = 100;

    geometry::TetraMesh tm;

    EXPECT_TRUE(tm.IsEmpty());

    tm.vertices_.resize(size);

    EXPECT_FALSE(tm.IsEmpty());
}

// ----------------------------------------------------------------------------
//
// ----------------------------------------------------------------------------
TEST(TetraMesh, GetMinBound) {
    int size = 100;

    Vector3d dmin(0.0, 0.0, 0.0);
    Vector3d dmax(1000.0, 1000.0, 1000.0);

    geometry::TetraMesh tm;

    tm.vertices_.resize(size);
    Rand(tm.vertices_, dmin, dmax, 0);

    ExpectEQ(Vector3d(19.607843, 0.0, 0.0), tm.GetMinBound());
}

// ----------------------------------------------------------------------------
//
// ----------------------------------------------------------------------------
TEST(TetraMesh, GetMaxBound) {
    int size = 100;

    Vector3d dmin(0.0, 0.0, 0.0);
    Vector3d dmax(1000.0, 1000.0, 1000.0);

    geometry::TetraMesh tm;

    tm.vertices_.resize(size);
    Rand(tm.vertices_, dmin, dmax, 0);

    ExpectEQ(Vector3d(996.078431, 996.078431, 996.078431), tm.GetMaxBound());
}

// ----------------------------------------------------------------------------
//
// ----------------------------------------------------------------------------
TEST(TetraMesh, Transform) {
    vector<Vector3d> ref_vertices = {{396.870588, 1201.976471, 880.472941},
                                     {320.792157, 1081.976471, 829.139608},
                                     {269.027451, 818.447059, 406.786667},
                                     {338.831373, 1001.192157, 614.237647},
                                     {423.537255, 1153.349020, 483.727843},
                                     {432.949020, 1338.447059, 964.512157},
                                     {140.007843, 444.721569, 189.296471},
                                     {292.164706, 763.152941, 317.178824},
                                     {134.517647, 407.858824, 192.002353},
                                     {274.909804, 802.368627, 218.747451}};

    int size = 10;

    Vector3d dmin(0.0, 0.0, 0.0);
    Vector3d dmax(1000.0, 1000.0, 1000.0);

    geometry::TetraMesh tm;

    tm.vertices_.resize(size);

    Rand(tm.vertices_, dmin, dmax, 0);

    Matrix4d transformation;
    transformation << 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90,
            0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16;

    tm.Transform(transformation);

    ExpectEQ(ref_vertices, tm.vertices_);
}

// ----------------------------------------------------------------------------
//
// ----------------------------------------------------------------------------
TEST(TetraMesh, OperatorAppend) {
    size_t size = 100;

    Vector3d dmin(0.0, 0.0, 0.0);
    Vector3d dmax(1000.0, 1000.0, 1000.0);

    Vector4i imin(0, 0, 0, 0);
    Vector4i imax(size - 1, size - 1, size - 1, size - 1);

    geometry::TetraMesh tm0;
    geometry::TetraMesh tm1;

    tm0.vertices_.resize(size);
    tm0.tetras_.resize(size);

    tm1.vertices_.resize(size);
    tm1.tetras_.resize(size);

    Rand(tm0.vertices_, dmin, dmax, 0);
    Rand(tm0.tetras_, imin, imax, 0);

    Rand(tm1.vertices_, dmin, dmax, 0);
    Rand(tm1.tetras_, imin, imax, 0);

    geometry::TetraMesh tm(tm0);
    tm += tm1;

    EXPECT_EQ(2 * size, tm.vertices_.size());
    for (size_t i = 0; i < size; i++) {
        ExpectEQ(tm0.vertices_[i], tm.vertices_[i + 0]);
        ExpectEQ(tm1.vertices_[i], tm.vertices_[i + size]);
    }

    EXPECT_EQ(2 * size, tm.tetras_.size());
    for (size_t i = 0; i < size; i++) {
        ExpectEQ(tm0.tetras_[i], tm.tetras_[i + 0]);
        ExpectEQ(Vector4i(tm1.tetras_[i](0, 0) + size,
                          tm1.tetras_[i](1, 0) + size,
                          tm1.tetras_[i](2, 0) + size,
                          tm1.tetras_[i](3, 0) + size),
                 tm.tetras_[i + size]);
    }
}

// ----------------------------------------------------------------------------
//
// ----------------------------------------------------------------------------
TEST(TetraMesh, OperatorADD) {
    size_t size = 100;

    Vector3d dmin(0.0, 0.0, 0.0);
    Vector3d dmax(1000.0, 1000.0, 1000.0);

    Vector4i imin(0, 0, 0, 0);
    Vector4i imax(size - 1, size - 1, size - 1, size - 1);

    geometry::TetraMesh tm0;
    geometry::TetraMesh tm1;

    tm0.vertices_.resize(size);
    tm0.tetras_.resize(size);

    tm1.vertices_.resize(size);
    tm1.tetras_.resize(size);

    Rand(tm0.vertices_, dmin, dmax, 0);
    Rand(tm0.tetras_, imin, imax, 0);

    Rand(tm1.vertices_, dmin, dmax, 0);
    Rand(tm1.tetras_, imin, imax, 0);

    geometry::TetraMesh tm = tm0 + tm1;

    EXPECT_EQ(2 * size, tm.vertices_.size());
    for (size_t i = 0; i < size; i++) {
        ExpectEQ(tm0.vertices_[i], tm.vertices_[i + 0]);
        ExpectEQ(tm1.vertices_[i], tm.vertices_[i + size]);
    }

    EXPECT_EQ(2 * size, tm.tetras_.size());
    for (size_t i = 0; i < size; i++) {
        ExpectEQ(tm0.tetras_[i], tm.tetras_[i + 0]);
        ExpectEQ(Vector4i(tm1.tetras_[i](0, 0) + size,
                          tm1.tetras_[i](1, 0) + size,
                          tm1.tetras_[i](2, 0) + size,
                          tm1.tetras_[i](3, 0) + size),
                 tm.tetras_[i + size]);
    }
}

// ----------------------------------------------------------------------------
//
// ----------------------------------------------------------------------------
TEST(TetraMesh, Purge) {
    vector<Vector3d> ref_vertices = {{839.215686, 392.156863, 780.392157},
                                     {796.078431, 909.803922, 196.078431},
                                     {333.333333, 764.705882, 274.509804},
                                     {552.941176, 474.509804, 627.450980},
                                     {364.705882, 509.803922, 949.019608},
                                     {913.725490, 635.294118, 713.725490},
                                     {141.176471, 603.921569, 15.686275},
                                     {239.215686, 133.333333, 803.921569},
                                     {152.941176, 400.000000, 129.411765},
                                     {105.882353, 996.078431, 215.686275},
                                     {509.803922, 835.294118, 611.764706},
                                     {294.117647, 635.294118, 521.568627},
                                     {490.196078, 972.549020, 290.196078},
                                     {768.627451, 525.490196, 768.627451},
                                     {400.000000, 890.196078, 282.352941},
                                     {349.019608, 803.921569, 917.647059},
                                     {66.666667, 949.019608, 525.490196},
                                     {82.352941, 192.156863, 662.745098},
                                     {890.196078, 345.098039, 62.745098},
                                     {19.607843, 454.901961, 62.745098},
                                     {235.294118, 968.627451, 901.960784},
                                     {847.058824, 262.745098, 537.254902},
                                     {372.549020, 756.862745, 509.803922},
                                     {666.666667, 529.411765, 39.215686}};

    vector<Vector3i> ref_triangles = {
            {20, 9, 18},  {19, 21, 4}, {8, 18, 6}, {13, 11, 15}, {8, 12, 22},
            {21, 15, 17}, {3, 14, 0},  {5, 3, 19}, {2, 23, 5},   {12, 20, 14},
            {7, 15, 12},  {11, 23, 6}, {9, 21, 6}, {8, 19, 22},  {1, 22, 12},
            {1, 4, 15},   {21, 8, 1},  {0, 10, 1}, {5, 23, 21},  {20, 6, 12},
            {8, 18, 12},  {16, 12, 0}};

    vector<Vector3d> ref_triangle_normals = {
            {839.215686, 392.156863, 780.392157},
            {796.078431, 909.803922, 196.078431},
            {333.333333, 764.705882, 274.509804},
            {552.941176, 474.509804, 627.450980},
            {364.705882, 509.803922, 949.019608},
            {913.725490, 635.294118, 713.725490},
            {141.176471, 603.921569, 15.686275},
            {239.215686, 133.333333, 803.921569},
            {105.882353, 996.078431, 215.686275},
            {509.803922, 835.294118, 611.764706},
            {294.117647, 635.294118, 521.568627},
            {490.196078, 972.549020, 290.196078},
            {400.000000, 890.196078, 282.352941},
            {349.019608, 803.921569, 917.647059},
            {66.666667, 949.019608, 525.490196},
            {82.352941, 192.156863, 662.745098},
            {890.196078, 345.098039, 62.745098},
            {19.607843, 454.901961, 62.745098},
            {235.294118, 968.627451, 901.960784},
            {847.058824, 262.745098, 537.254902},
            {372.549020, 756.862745, 509.803922},
            {666.666667, 529.411765, 39.215686}};

    int size = 25;

    Vector3d dmin(0.0, 0.0, 0.0);
    Vector3d dmax(1000.0, 1000.0, 1000.0);

    Vector4i imin(0, 0, 0, 0);
    Vector4i imax(size - 1, size - 1, size - 1, size - 1);

    geometry::TetraMesh tm0;
    geometry::TetraMesh tm1;

    tm0.vertices_.resize(size);
    tm0.tetras_.resize(size);

    tm1.vertices_.resize(size);
    tm1.tetras_.resize(size);

    Rand(tm0.vertices_, dmin, dmax, 0);
    Rand(tm0.tetras_, imin, imax, 0);

    Rand(tm1.vertices_, dmin, dmax, 0);
    Rand(tm1.tetras_, imin, imax, 0);

    geometry::TetraMesh tm = tm0 + tm1;

    tm.RemoveDuplicatedVertices();
    tm.RemoveDuplicatedTetras();
    tm.RemoveUnreferencedVertices();
    tm.RemoveDegenerateTetras();

    ExpectEQ(ref_vertices, tm.vertices_);
    // ExpectEQ(ref_triangles, tm.tetras_);
}

// ----------------------------------------------------------------------------
//
// ----------------------------------------------------------------------------
TEST(TetraMesh, HasVertices) {
    int size = 100;

    geometry::TetraMesh tm;

    EXPECT_FALSE(tm.HasVertices());

    tm.vertices_.resize(size);

    EXPECT_TRUE(tm.HasVertices());
}

// ----------------------------------------------------------------------------
//
// ----------------------------------------------------------------------------
TEST(TetraMesh, HasTetras) {
    int size = 100;

    geometry::TetraMesh tm;

    EXPECT_FALSE(tm.HasTetras());

    tm.vertices_.resize(size);
    tm.tetras_.resize(size);

    EXPECT_TRUE(tm.HasTetras());
}
