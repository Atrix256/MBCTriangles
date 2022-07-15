#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <random>
#include <vector>
#include <array>

static const int c_numPoints = 100;
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
std::array<T, N> operator-(const std::array<T, N>& A, const std::array<T, N>& B)
{
    std::array<T, N> ret;
    for (size_t i = 0; i < N; ++i)
        ret[i] = A[i] - B[i];
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

int main(int argc, char** argv)
{
    std::vector<Vec2> points;

#if DETERMINISTIC()
    std::mt19937 rng;
#else
    std::random_device rd;
    std::mt19937 rng(rd());
#endif
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);

    for (int pointIndex = 0; pointIndex < c_numPoints; ++pointIndex)
    {
        int numCandidates = (int)points.size() + 1;

        Vec2 bestCandidatePos;
        float bestCandidateScore = -FLT_MAX;

        for (int candidateIndex = 0; candidateIndex < numCandidates; ++candidateIndex)
        {
            // Find a point in the triangle through rejection sampling
            Vec2 candidatePos = Vec2{ dist(rng), dist(rng) };
            while (!PointInTriangle(candidatePos, PointA, PointB, PointC))
                candidatePos = Vec2{ dist(rng), dist(rng) };

            float candidateScore = FLT_MAX;
            for (const Vec2& point : points)
            {
                float distance = DistanceBetweenPoints();
            }

            // TODO: calculate score

            if (candidateScore > bestCandidateScore)
            {
                bestCandidateScore = candidateScore;
                bestCandidatePos = candidatePos;
            }
        }

        // Add the point
        points.push_back(bestCandidatePos);
    }

    return 0;
}

// TODO: switch to working in barycentric coordinates, not actual coordinates. generate u, v, calculate w.

/*
TODO:
- draw out the points to an image. different numbers of points
- draw out the points tiled randomly
- generate candidates using distances at each of the 3 edges, 3 orientations, and flipped/not

- TODO: convert points to barycentric coordinates

*/