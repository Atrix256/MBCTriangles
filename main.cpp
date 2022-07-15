#define _CRT_SECURE_NO_WARNINGS

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include <random>
#include <vector>
#include <array>
#include <direct.h>

static const int c_numPoints = 500;
#define DETERMINISTIC() false


typedef std::array<float, 2> Vec2;
typedef std::array<float, 3> Vec3;

static const float c_pi = 3.14159265359f;

Vec2 Calculate3rdEquilateralTrianglePoint(const Vec2& A, const Vec2& B);

static const Vec2 PointA{ 0.0f, 0.0f };
static const Vec2 PointB{ 1.0f, 0.0f };
static const Vec2 PointC = Calculate3rdEquilateralTrianglePoint(PointA, PointB);

// https://stackoverflow.com/questions/2861904/how-to-find-coordinates-of-a-2d-equilateral-triangle-in-c
Vec2 Calculate3rdEquilateralTrianglePoint(const Vec2& A, const Vec2& B)
{
    float s60 = sin(-60.0f * c_pi / 180.0f);
    float c60 = cos(-60.0f * c_pi / 180.0f);

    return Vec2
    {
      c60 * (A[0] - B[0]) - s60 * (A[1] - B[1]) + B[0],
      s60 * (A[0] - B[0]) + c60 * (A[1] - B[1]) + B[1]
    };
}

template <typename T, size_t N>
std::array<T, N> operator+(const std::array<T, N>& A, const std::array<T, N>& B)
{
    std::array<T, N> ret;
    for (size_t i = 0; i < N; ++i)
        ret[i] = A[i] + B[i];
    return ret;
}

template <typename T, size_t N>
std::array<T, N> operator-(const std::array<T, N>& A, const std::array<T, N>& B)
{
    std::array<T, N> ret;
    for (size_t i = 0; i < N; ++i)
        ret[i] = A[i] - B[i];
    return ret;
}

template <typename T, size_t N>
std::array<T, N> operator*(const std::array<T, N>& A, const T& B)
{
    std::array<T, N> ret;
    for (size_t i = 0; i < N; ++i)
        ret[i] = A[i] * B;
    return ret;
}

template <typename T, size_t N>
std::array<T, N> operator/(const std::array<T, N>& A, const T& B)
{
    std::array<T, N> ret;
    for (size_t i = 0; i < N; ++i)
        ret[i] = A[i] / B;
    return ret;
}

template <typename T, size_t N>
T Dot(const std::array<T, N>& A, const std::array<T, N>& B)
{
    T ret = T(0);
    for (size_t i = 0; i < N; ++i)
        ret += A[i] * B[i];
    return ret;
}

template <typename T, size_t N>
T Length(const std::array<T, N>& A)
{
    return sqrtf(Dot(A, A));
}

float SmoothStep(float edge0, float edge1, float x)
{
    if (x < edge0)
        return 0;

    if (x >= edge1)
        return 1;

    // Scale/bias into [0..1] range
    x = (x - edge0) / (edge1 - edge0);

    return x * x * (3 - 2 * x);
}

float Lerp(float A, float B, float t)
{
    return A * (1.0f - t) + B * t;
}

float Clamp(float value, float themin, float themax)
{
    if (value <= themin)
        return themin;
    if (value >= themax)
        return themax;
    return value;
}

// TODO: may not need this, dot, or the - operator
// Compute barycentric coordinates (u, v, w) for
// point p with respect to triangle (a, b, c)
Vec3 Barycentric(const Vec2& p, const Vec2& a, const Vec2& b, const Vec2& c)
{
    Vec2 v0 = b - a;
    Vec2 v1 = c - a;
    Vec2 v2 = p - a;
    float d00 = Dot(v0, v0);
    float d01 = Dot(v0, v1);
    float d11 = Dot(v1, v1);
    float d20 = Dot(v2, v0);
    float d21 = Dot(v2, v1);
    float denom = d00 * d11 - d01 * d01;

    Vec3 ret;
    ret[1] = (d11 * d20 - d01 * d21) / denom;
    ret[2] = (d00 * d21 - d01 * d20) / denom;
    ret[0] = 1.0f - ret[1] - ret[2];
    return ret;
}

// TODO: maybe don't need this
float Signed2DTriArea(const Vec2& A, const Vec2& B, const Vec2& C)
{
    return (A[0] - C[0]) * (B[1] - C[1]) - (A[1] - C[1]) * (B[0] - C[0]);
}

bool PointInTriangle(const Vec2& P, const Vec2& A, const Vec2& B, const Vec2& C)
{
    Vec3 bary = Barycentric(P, A, B, C);

    return
        bary[0] >= 0.0f && bary[0] <= 1.0f &&
        bary[1] >= 0.0f && bary[1] <= 1.0f &&
        bary[2] >= 0.0f && bary[2] <= 1.0f;
}

bool BaryPointInTriangle(const Vec3& bary)
{
    return
        bary[0] >= 0.0f && bary[0] <= 1.0f &&
        bary[1] >= 0.0f && bary[1] <= 1.0f &&
        bary[2] >= 0.0f && bary[2] <= 1.0f;
}

Vec2 BarycentricToCartesian(const Vec3& bary, const Vec2& A, const Vec2& B, const Vec2& C)
{
    return A * bary[0] + B * bary[1] + C * bary[2];
}

Vec2 ReflectPoint(const Vec2& P, const Vec2& A, const Vec2& B)
{
    Vec2 AB = B - A;
    Vec2 ABNorm = AB / Length(AB);

    float t = Dot(P - A, ABNorm);
    Vec2 linePoint = A + ABNorm * t;

    Vec2 PToLine = (P - linePoint);
    return P - PToLine * 2.0f;
}

float MinDistanceBetweenBaryPoints(const Vec3& candidatePos, const Vec3& point)
{
    Vec2 A = BarycentricToCartesian(candidatePos, PointA, PointB, PointC);
    Vec2 B = BarycentricToCartesian(point, PointA, PointB, PointC);

    // direct distance
    float distance = Length(B - A);

#if 0
    // reflect one of the points across each edge and keep the minimum distance found
    distance = std::min(distance, Length(ReflectPoint(B, PointA, PointB) - A));
    distance = std::min(distance, Length(ReflectPoint(B, PointB, PointC) - A));
    distance = std::min(distance, Length(ReflectPoint(B, PointC, PointA) - A));

    // Also check the candidate point against it's reflected self
    distance = std::min(distance, Length(ReflectPoint(A, PointA, PointB) - A));
    distance = std::min(distance, Length(ReflectPoint(A, PointB, PointC) - A));
    distance = std::min(distance, Length(ReflectPoint(A, PointC, PointA) - A));

#endif
    // TODO: this doesn't handle other edges being against this edge. (rotations?) need to do that.

    return distance;
}

void DrawPoints(const std::vector<Vec3>& points, int pointCount, const char* fileName)
{
    float minx = std::min(PointA[0], std::min(PointB[0], PointC[0]));
    float maxx = std::max(PointA[0], std::max(PointB[0], PointC[0]));
    float miny = std::min(PointA[1], std::min(PointB[1], PointC[1]));
    float maxy = std::max(PointA[1], std::max(PointB[1], PointC[1]));

    float offsety = 0.5f - (maxy - miny) / 2.0f;

    static const int c_imgSize = 512;

    std::vector<unsigned char> pixels(c_imgSize * c_imgSize * 3, 255);

    unsigned char* pixel = pixels.data();
    for (int iy = 0; iy < c_imgSize; ++iy)
    {
        for (int ix = 0; ix < c_imgSize; ++ix)
        {
            Vec2 uv = { float(ix) / float(c_imgSize), 1.0f - (offsety + float(iy) / float(c_imgSize)) };
            uv = uv * 2.0f;
            uv[0] = uv[0] - 0.5f;
            uv[1] = uv[1] - 0.5f;

            Vec3 bary = Barycentric(uv, PointA, PointB, PointC);

            Vec3 color = Vec3{ 1.0f, 1.0f, 1.0f };

            // handle edges
            if (BaryPointInTriangle(bary))
            {
                float offset = -0.01f;

                float dist = std::min(std::min(std::abs(bary[0]+ offset), std::abs(bary[1]+ offset)), std::abs(bary[2]+ offset));
                float shade = SmoothStep(0.0025f, 0.01f, dist);

                color[0] = Lerp(0.1f, 1.0f, shade);
                color[1] = Lerp(0.1f, 1.0f, shade);
                color[2] = Lerp(0.1f, 1.0f, shade);
            }

            // handle points
            for (size_t pointIndex = 0; pointIndex < pointCount; ++pointIndex)
            {
                const Vec3& point = points[pointIndex];

                Vec2 pointCar = BarycentricToCartesian(point, PointA, PointB, PointC);

                float dist = Length(pointCar - uv);
                float shade = 1.0f - SmoothStep(0.0025f, 0.01f, dist);

                color[0] = Lerp(color[0], 0.0f, shade);
                color[1] = Lerp(color[1], 0.0f, shade);
                color[2] = Lerp(color[2], 1.0f, shade);
            }

            pixel[0] = (unsigned char)Clamp(color[0] * 255.0f, 0.0f, 255.0f);
            pixel[1] = (unsigned char)Clamp(color[1] * 255.0f, 0.0f, 255.0f);
            pixel[2] = (unsigned char)Clamp(color[2] * 255.0f, 0.0f, 255.0f);

            pixel += 3;
        }
    }

    stbi_write_png(fileName, c_imgSize, c_imgSize, 3, pixels.data(), 0);
}

void DrawPointsRepeating(const std::vector<Vec3>& points, int pointCount, const char* fileName)
{
    float minx = std::min(PointA[0], std::min(PointB[0], PointC[0]));
    float maxx = std::max(PointA[0], std::max(PointB[0], PointC[0]));
    float miny = std::min(PointA[1], std::min(PointB[1], PointC[1]));
    float maxy = std::max(PointA[1], std::max(PointB[1], PointC[1]));

    float offsety = 0.5f - (maxy - miny) / 2.0f;

    static const int c_imgSize = 512;

    std::vector<unsigned char> pixels(c_imgSize * c_imgSize * 3, 255);

    unsigned char* pixel = pixels.data();
    for (int iy = 0; iy < c_imgSize; ++iy)
    {
        for (int ix = 0; ix < c_imgSize; ++ix)
        {
            Vec2 uv = { float(ix) / float(c_imgSize), 1.0f - (offsety + float(iy) / float(c_imgSize)) };
            uv = uv * 2.0f;
            uv[0] = uv[0] - 0.5f;
            uv[1] = uv[1] - 0.5f;

            Vec3 bary = Barycentric(uv, PointA, PointB, PointC);

            Vec3 color = Vec3{ 1.0f, 1.0f, 1.0f };

            // handle points
            for (size_t pointIndex = 0; pointIndex < pointCount; ++pointIndex)
            {
                const Vec3& point = points[pointIndex];

                Vec2 pointCar = BarycentricToCartesian(point, PointA, PointB, PointC);

                {
                    float dist = Length(pointCar - uv);
                    float shade = 1.0f - SmoothStep(0.0025f, 0.01f, dist);

                    color[0] = Lerp(color[0], 0.0f, shade);
                    color[1] = Lerp(color[1], 0.0f, shade);
                    color[2] = Lerp(color[2], 1.0f, shade);
                }

                {
                    float dist = Length(ReflectPoint(pointCar, PointA, PointB) - uv);
                    float shade = 1.0f - SmoothStep(0.0025f, 0.01f, dist);

                    color[0] = Lerp(color[0], 0.0f, shade);
                    color[1] = Lerp(color[1], 0.0f, shade);
                    color[2] = Lerp(color[2], 1.0f, shade);
                }

                {
                    float dist = Length(ReflectPoint(pointCar, PointB, PointC) - uv);
                    float shade = 1.0f - SmoothStep(0.0025f, 0.01f, dist);

                    color[0] = Lerp(color[0], 0.0f, shade);
                    color[1] = Lerp(color[1], 0.0f, shade);
                    color[2] = Lerp(color[2], 1.0f, shade);
                    {
                        float dist = Length(ReflectPoint(pointCar, PointC, PointA) - uv);
                        float shade = 1.0f - SmoothStep(0.0025f, 0.01f, dist);

                        color[0] = Lerp(color[0], 0.0f, shade);
                        color[1] = Lerp(color[1], 0.0f, shade);
                        color[2] = Lerp(color[2], 1.0f, shade);
                    }
                }
            }

            pixel[0] = (unsigned char)Clamp(color[0] * 255.0f, 0.0f, 255.0f);
            pixel[1] = (unsigned char)Clamp(color[1] * 255.0f, 0.0f, 255.0f);
            pixel[2] = (unsigned char)Clamp(color[2] * 255.0f, 0.0f, 255.0f);

            pixel += 3;
        }
    }

    stbi_write_png(fileName, c_imgSize, c_imgSize, 3, pixels.data(), 0);
}

int main(int argc, char** argv)
{
    _mkdir("out");

    std::vector<Vec3> points;

#if DETERMINISTIC()
    std::mt19937 rng;
#else
    std::random_device rd;
    std::mt19937 rng(rd());
#endif
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);

    for (int pointIndex = 0; pointIndex < c_numPoints; ++pointIndex)
    {
        // Generate multiple candidates for the next point and keep the best one
        int numCandidates = (int)points.size() + 1;
        Vec3 bestCandidatePos;
        float bestCandidateScore = -FLT_MAX;
        for (int candidateIndex = 0; candidateIndex < numCandidates; ++candidateIndex)
        {
            // Find a point in the triangle through rejection sampling
            Vec3 candidatePos = Vec3{ dist(rng), dist(rng), 0.0f };
            candidatePos[2] = 1.0f - candidatePos[0] - candidatePos[1];
            while (!BaryPointInTriangle(candidatePos))
            {
                candidatePos = Vec3{ dist(rng), dist(rng), 0.0f };
                candidatePos[2] = 1.0f - candidatePos[0] - candidatePos[1];
            }

            // calculate the distance to the nearest point
            float candidateScore = FLT_MAX;
            for (const Vec3& point : points)
                candidateScore = std::min(candidateScore, MinDistanceBetweenBaryPoints(candidatePos, point));

            // keep the candidate that has the largest distance to the nearest point
            if (candidateScore > bestCandidateScore)
            {
                bestCandidateScore = candidateScore;
                bestCandidatePos = candidatePos;
            }
        }

        // Add the point
        points.push_back(bestCandidatePos);
    }

    DrawPoints(points, int(float(points.size() * 0.1f)), "out/points_10.png");
    DrawPoints(points, int(float(points.size() * 0.25f)), "out/points_25.png");
    DrawPoints(points, int(float(points.size() * 0.50f)), "out/points_50.png");
    DrawPoints(points, int(float(points.size() * 0.75f)), "out/points_75.png");
    DrawPoints(points, int(points.size()), "out/points_100.png");

    DrawPointsRepeating(points, int(float(points.size() * 0.1f)), "out/points_r_10.png");
    DrawPointsRepeating(points, int(float(points.size() * 0.25f)), "out/points_r_25.png");
    DrawPointsRepeating(points, int(float(points.size() * 0.50f)), "out/points_r_50.png");
    DrawPointsRepeating(points, int(float(points.size() * 0.75f)), "out/points_r_75.png");
    DrawPointsRepeating(points, int(points.size()), "out/points_r_100.png");

    // write out file
    {
        FILE* file = nullptr;
        fopen_s(&file, "out/points.txt", "w+t");

        fprintf(file, "Barycentric coordinates of blue noise points on an equilateral triangle\n\n");

        for (Vec3& point : points)
            fprintf(file, "    {%f, %f, %f},\n", point[0], point[1], point[2]);

        fclose(file);
    }

    return 0;
}

/*
TODO:
- draw out the points to an image. different numbers of points. Could draw a non equilateral triangle too.
- draw out the points tiled randomly
- generate candidates using distances at each of the 3 edges, 3 orientations, and flipped/not

*/