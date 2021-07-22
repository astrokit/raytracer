


#ifndef INCLUDED_TASK1B
#define INCLUDED_TASK1B

#pragma once

#include <cstddef>
#include <array>

#include <framework/math/vector.h>
#include <framework/math/matrix.h>
#include <framework/image.h>


class Scene;

struct Cube;
struct Plane;

struct HitPoint
{
	float3 position;
	float3 normal;
	float3 k_d;
	float3 k_s;
	float m;
};

struct Camera
{
	float w_s;
	float f;
	float3 eye;
	float3 lookat;
	float3 up;
};

struct Pointlight
{
	float3 position;
	float3 color;
};

struct Spotlight
{
  float3 position;
  float3 color;
  float3 direction;
  float angle;
};

float4x4 getRotationMatrix(const float3& rotation);
float4x4 getTranslationMatrix(const float3& position);
float4x4 getModelMatrix(const float4x4& translation_matrix, const float4x4& rotation_matrix);

bool intersectsRayCube(const float3& p, const float3& d, const Cube* cubes, std::size_t num_cubes, float t_min,
	float t_max);
const Cube* findClosestHitCubes(const float3& p, const float3& d, const Cube* cubes, std::size_t num_cubes,
	float& t);
bool intersectRayPlane(const float3& p, const float3& d, const float4& plane, float& t);
bool intersectsRayPlane(const float3& p, const float3& d, const Plane* planes, std::size_t num_planes, float t_min, float t_max);
const Plane* findClosestHitPlanes(const float3& p, const float3& d, const Plane* planes, std::size_t num_planes, float& t);
bool intersectRayCube(const float3& p, const float3& d, const Cube* cube, float& t);


float3 shade(const float3& p, const float3& d, const HitPoint& hit, const Scene& scene, const Pointlight* lights, std::size_t num_lights, const Spotlight* spotlights, std::size_t num_spotlights);

void render(image2D<float3>& framebuffer, int left, int top, int right, int bottom, const Scene& scene, const Camera& camera, const Pointlight* lights, std::size_t num_lights, const Spotlight* spotlights, std::size_t num_spotlights, const float3& background_color, int max_bounces);

#endif  // INCLUDED_TASK1B
