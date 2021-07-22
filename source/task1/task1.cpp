

#include <limits>
#include "Scene.h"
#include "task1.h"

constexpr float epsilon = 0.0001f;

#define FULL_VERSION

float4x4 getRotationMatrix(const float3& rotation)
{
    // TODO return the rotation matrix according to the 
    // x, y, and z angles given by rotation.

  float sin_x = std::sin(rotation.x);
  float cos_x = std::cos(rotation.x);
  float sin_y = std::sin(rotation.y);
  float cos_y = std::cos(rotation.y);
  float sin_z = std::sin(rotation.z);
  float cos_z = std::cos(rotation.z);

  float4x4 R_x = math::float4x4(1.0f, 0.0f, 0.0f, 0.0f,
                                0.0f, cos_x, -sin_x, 0.0f,
                                0.0f, sin_x, cos_x, 0.0f,
                                0.0f, 0.0f, 0.0f, 1.0f
  );

  float4x4 R_y = math::float4x4(cos_y, 0, sin_y, 0.0f,
                                0.0f, 1.0f, 0.0f, 0.0f,
                                -sin_y, 0.0f, cos_y, 0.0f,
                                0.0f, 0.0f, 0.0f, 1.0f
  );

  float4x4 R_z = math::float4x4(cos_z, -sin_z, 0.0f, 0.0f,
                                sin_z,  cos_z, 0.0f, 0.0f,
                                0.0f, 0.0f, 1.0f, 0.0f,
                                0.0f, 0.0f, 0.0f, 1.0f
  );

  return R_z * R_y * R_x;
}

float4x4 getTranslationMatrix(const float3& position)
{
    // TODO return the translation matrix according to the x, y, and z coordinates
    // of the given position.

  float4x4 translation_matrix = math::float4x4(1, 0, 0, position.x,
                                               0, 1, 0, position.y,
                                               0, 0, 1, position.z,
                                               0, 0, 0, 1);

  return translation_matrix;
}

float4x4 getModelMatrix(const float4x4& translation_matrix, const float4x4& rotation_matrix)
{
    // TODO return the model matrix

  return translation_matrix * rotation_matrix;
}

bool intersectRayPlane(const float3& p, const float3& d, const float4& plane, float& t)
{
    // TODO implement the intersection test between a ray and a plane.
    //
    // A plane is defined by the plane normal n and the offset along the normal w.
    // The plane is given as parameter plane where (plane.x, plane.y, plane.z) represent the plane normal
    // and plane.w is the offset w.
    //
    // If there is no intersection (Hint: or one we do not care about), return false.
    // Otherwise, compute and set the parameter t such that p + t * d yields the intersection point and return true.

  float3 normal = { plane.x, plane.y, plane.z };

  float check = dot(d, normal);
  if (check != 0) // -> ray cuts plane exactly at one point
  {
    t = (plane.w - dot(p, normal)) / check;

    return (t >= 0); // 0 = eye -> intersection only useful if not behind the eye
  }

  return false;
}


bool intersectRayCube(const float3& p, const float3& d, 
  const Cube* cube, float& t)
{
    // TODO Implement the intersection test between a ray and a cube.
    //
    // The parameter cube is a pointer to a struct Cube that holds 
    // neccessary information about the position, orientation and 
    // scaling of the cube (see Scene.h).
    //
    // 1. Use the inverse model and rotation matrices to transfrom 
    //    the ray from world space into the cube's local space (object space).
    // 2. Compute b_min and b_max of the cube using the cube's scale.
    // 3. Compute the intersection between the local ray and the
    //    axis aligned bounding box (AABB) given by b_min and b_max.
    // 4. If there is a valid intersection, set the parameter t (distance to intersection)
    //    and return true, otherwise return false.
    //
    // Note: As there is no scaling involved when transforming the ray into the cube's
    //       local space, the distance to the intersection (t) is the same in world and object space.

  // 1. Use the inverse model and rotation matrices to transfrom
  //    the ray from world space into the cube's local space (object space).

  float4 p_local = cube->inverse_model_matrix * float4(p, 1.0f);
  float4 d_local = cube->inverse_rotation_matrix * float4(d, 1.0f);

  // 2. Compute b_min and b_max of the cube using the cube's scale.

  float3 b_min = -cube->scale * 0.5f;
  float3 b_max = cube->scale * 0.5f;

  // 3. Compute the intersection between the local ray and the
  //    axis aligned bounding box (AABB) given by b_min and b_max.

  float tx1 = (b_min.x - p_local.x) / d_local.x;
  float tx2 = (b_max.x - p_local.x) / d_local.x;

  float ty1 = (b_min.y - p_local.y) / d_local.y;
  float ty2 = (b_max.y - p_local.y) / d_local.y;

  float tz1 = (b_min.z - p_local.z) / d_local.z;
  float tz2 = (b_max.z - p_local.z) / d_local.z;

  float t_near = std::max(std::max(std::min(tx1, tx2), std::min(ty1, ty2)), std::min(tz1, tz2));
  float t_far = std::min(std::min(std::max(tx1, tx2), std::max(ty1, ty2)), std::max(tz1, tz2));

  // 4. If there is a valid intersection, set the parameter t (distance to intersection)
  //    and return true, otherwise return false.

  if (t_near <= t_far)
  {
    t = t_near;
    return true;
  }

  return false;
}

bool intersectsRayPlane(const float3& p, const float3& d, const Plane* planes, std::size_t num_planes, float t_min, float t_max)
{
    // TODO: implement intersection test between a ray and a set of planes.
    // This method only has to detect whether there is an intersection with ANY
    // plane along the given subset of the ray. The ray is given by its start
    // point p and direction d. 
    // A plane is defined by the plane normal n and the offset along the normal w.
    // Each plane in planes contains a float4 member p where the plane normal n is 
    // (p.x, p.y, p.z) and w is p.w. 
    // If an intersection is found that falls on a point on the ray
    // between t_min and t_max, return true. Otherwise, return false.

  float t_hit;
  for (size_t i = 0; i < num_planes; i++)
  {
    if (intersectRayPlane(p, d, planes[i].p, t_hit) && t_hit >= t_min && t_hit <= t_max)
      return true;
  }

  return false;
}

const Plane* findClosestHitPlanes(const float3& p, const float3& d, const Plane* planes, std::size_t num_planes, float& t)
{
    // TODO: implement intersection test between a ray and a set of planes.
    // This function should find the CLOSEST intersection with a plane along
    // the ray. The ray is given by its start point p and direction d.
    // A plane is defined by its distance to the origin w and the plane normal n.
    // Each plane in planes contains a float4 member p where the plane normal n is 
    // (p.x, p.y, p.z) and w is p.w.
    // If an intersection is found, set t to the ray parameter and
    // return a pointer to the hit plane.
    // If no intersection is found, return nullptr.

  const Plane* hit = nullptr;
  float t_hit;

  for (size_t i = 0; i < num_planes; i++)
  {
    // Intersection is closer than previous if t_hit < current t
    if (intersectRayPlane(p, d, planes[i].p, t_hit) && t_hit < t)
    {
      t = t_hit;
      hit = &planes[i];
    }
  }

  return hit;
}

bool intersectsRayCube(const float3& p, const float3& d, const Cube* cubes, std::size_t num_cubes, float t_min,
    float t_max)
{
    // TODO: implement intersection test between a ray and a set of cubes.
    // This method only has to detect whether there is an intersection with ANY
    // cube along the given subset of the ray. The ray is given by its start
    // point p and direction d.
    // If an intersection is found that falls on a point on the ray
    // between t_min and t_max, return true. Otherwise, return false.

  float t_hit;
  for (size_t i = 0; i < num_cubes; i++)
  {
    if (intersectRayCube(p, d, &cubes[i], t_hit) && t_hit >= t_min && t_hit <= t_max)
      return true;
  }

  return false;
}

const Cube* findClosestHitCubes(const float3& p, const float3& d, const Cube* cubes, std::size_t num_cubes,
    float& t)
{
    // TODO: implement intersection test between a ray and a set of cubes.
    // This function should find the CLOSEST intersection with a cube along
    // the ray. The ray is given by its start point p and direction d.
    // If an intersection is found, set t to the ray parameter and
    // return a pointer to the hit cube.
    // If no intersection is found, return nullptr.

  const Cube* hit = nullptr;
  float t_hit;

  for (size_t i = 0; i < num_cubes; i++)
  {
    // Intersection is closer than previous if t_hit < current t
    if (intersectRayCube(p, d, &cubes[i], t_hit) && t_hit < t)
    {
      t = t_hit;
      hit = &cubes[i];
    }
  }

  return hit;
}


float3 shade(const float3& p, const float3& d, const HitPoint& hit,
             const Scene& scene, const Pointlight* lights,
             std::size_t num_lights, const Spotlight* spotlights,
             std::size_t num_spotlights) {
  // TODO: implement phong shading for point and spot ligths.
  // hit represents the surface point to be shaded.
  // hit.position, hit.normal, and hit.k_d and hit.k_s contain the position,
  // surface normal, and diffuse and specular reflection coefficients,
  // hit.m the specular power.
  // lights is a pointer to the first element of an array of num_lights
  // point light sources to consider.
  // Each light contains a member to give its position and color.
  // spotlights is a pointer to the first element of an array of num_spotlights
  // spot light sources to consider.
  // Each light contains members to give its position, color, direction and falloff angle.
  // Return the shaded color.

  // To implement shadows, use scene.intersectsRay(p, d, t_min, t_max) to test
  // whether a ray given by start point p and direction d intersects any
  // object on the section between t_min and t_max.

  /*
   * p = illuminated point | hit.position
   * n = surface normal | hit.normal
   * l = direction of light source
   * r = direction where light source is reflected to
   * v = part of light reflection the eye sees
   * theta = angle between n and l
   * alpha_r = angle between v and r
   *
   * k_d = diffuse reflection coefficient -> reflected to r | hit.k_d
   * c_d = diffuse part -> equally spread into all directions
   *     = k_d * max(cos(theta), 0)
   *
   * k_s = specular reflection coefficient | hit.k_s
   * c_s = specular part -> reflected into direction of v
   *     = (cos(theta) > 0) ? k_s * pow(max(cos(alpha_r), 0), m) : 0
   * m = sharpness level of specular reflection | hit.m
   */

  float3 color = { 0.0f, 0.0f, 0.0f };
  float3 shifted_hit_position = hit.position + hit.normal * epsilon;

  for (size_t i = 0; i < num_lights; i++)
  {
    // If ray starting at shifted_hit_position towards light is intersected by some object before it hits the position
    // of the light -> point is in shadow
    float3 shadow_ray_direction = lights[i].position - shifted_hit_position;
    if (!std::isinf(hit.m) && scene.intersectsRay(shifted_hit_position,
          normalize(shadow_ray_direction), 0, length(shadow_ray_direction)))
      continue;

    float3 l = normalize(lights[i].position - hit.position);

    float theta_cos = dot(hit.normal, l); // = cos(theta)
    float3 c_d = hit.k_d * std::max(theta_cos, 0.0f); // Lambert Shading

    float3 r = normalize(2.0f * hit.normal * theta_cos - l);
    float alpha_r_cos = dot(-d, r); // = cos(alpha_r)

    float3 c_s = { 0.0f, 0.0f, 0.0f };
    if (theta_cos > 0.0f)
      c_s = hit.k_s * std::pow(std::max(alpha_r_cos, 0.0f), hit.m);

    float3 cc = c_d + c_s;
    float3 light_color = lights[i].color;

    color += { cc.x * light_color.x, cc.y * light_color.y, cc.z * light_color.z }; // Hadamard multiplication
  }

  for (size_t i = 0; i < num_spotlights; i++)
  {
    // If ray starting at shifted_hit_position towards light is intersected by some object before it hits the position
    // of the light -> point is in shadow
    float3 shadow_ray_direction = spotlights[i].position - shifted_hit_position;
    if (!std::isinf(hit.m) && scene.intersectsRay(shifted_hit_position,
          normalize(shadow_ray_direction), 0, length(shadow_ray_direction)))
      continue;

    float3 l = normalize(spotlights[i].position - hit.position);

    float theta_cos = dot(hit.normal, l);
    float3 c_d = hit.k_d * std::max(theta_cos, 0.0f); // Lambert Shading

    float3 r = normalize(2.0f * hit.normal * theta_cos - l);
    float alpha_r_cos = dot(-d, r); // = cos(alpha_r)

    float3 c_s = { 0.0f, 0.0f, 0.0f };
    if (theta_cos > 0.0f)
      c_s = hit.k_s * std::pow(std::max(alpha_r_cos, 0.0f), hit.m);

    // Check if object gets some light at all and if yes, calculate attenuation
    float s = 0.0f;
    float theta = spotlights[i].angle;
    float alpha_cos = dot(-l, spotlights[i].direction);
    float alpha = std::acos(alpha_cos);
    if (alpha <= theta / 2.0f)
    {
      //s = 1.0f;
      // -> Less light around corners
      s = std::clamp((alpha_cos - std::cos(theta / 2.0f)) /
                     (1.0f - std::cos(theta / 2.0f)), 0.0f, 1.0f);
    }

    float3 cc = (c_d + c_s) * s;
    float3 light_color = spotlights[i].color;

    color += { cc.x * light_color.x, cc.y * light_color.y, cc.z * light_color.z }; // Hadamard multiplication
  }

  return color;
}

void render(image2D<float3>& framebuffer, int left, int top, int right,
            int bottom, const Scene& scene, const Camera& camera,
            const Pointlight* lights, std::size_t num_lights,
            const Spotlight* spotlights, std::size_t num_spotlights,
            const float3& background_color, int max_bounces) {
  // TODO: implement raytracing, render the given portion of the framebuffer.
  // left, top, right, and bottom specify the part of the image to compute
  // (inclusive left, top and exclusive right and bottom).
  // camera.eye, camera.lookat, and camera.up specify the position and
  // orientation of the camera, camera.w_s the width of the image plane,
  // and camera.f the focal length to use.
  // Use scene.findClosestHit(p, d) to find the closest point where the ray
  // hits an object.
  // The method returns an std::optional<HitPoint>.
  // If an object is hit, call the function shade to compute the color value
  // of the HitPoint illuminated by the given array of lights.
  // If the ray does not hit an object, use background_color.

  // BONUS: extend your implementation to recursive ray tracing.
  // max_bounces specifies the maximum number of secondary rays to trace.

  float r_x = (float) width(framebuffer);
  float r_y = (float) height(framebuffer);
  float h_s = r_y / r_x * camera.w_s; // height of image plane

  // Camera vectors u, v, w
  float3 w = normalize(camera.eye - camera.lookat); // = z
  float3 u = normalize(cross(camera.up, w)); // = x
  float3 v = cross(w, u); // = y

  float3 image_center = camera.eye - camera.f * w;
  float3 top_left_pixel = image_center + (camera.w_s / 2.0f) * (-u) + (h_s / 2.0f) * v;

  for (int y = top; y < bottom; ++y)
  {
    for (int x = left; x < right; ++x)
    {

      /*
       * r(t) = p + t * d
       *
       * p = start position of ray
       * d = ray direction
       * t = ray parameter > 0 -> distance to image -> || d ||
       *
       * camera.eye = position of eye point
       * camera.lookat = lookat point
       * camera.f = focal width
       * camera.w_s = width of image plane
       * camera.up = up vector
       */

      float3 p = top_left_pixel + (x + 0.5f) / r_x * camera.w_s * u + (y + 0.5f) / r_y * h_s * (-v);
      float3 d = normalize(p - camera.eye);

      std::optional<HitPoint> hit = scene.findClosestHit(p, d);

      /*
       * Recursive raytracing to handle mirror surfaces
       * (actually done iteratively because I was not sure if we are allowed to change the method signature)
       *
       * 1. Send reflection ray if hit is a mirror
       * 2. Multiply current hit->k_s with total hit_k_s value that is afterwards weighed with the final color
       * 3. Repeat until no more bounces left or ray does not hit anything
       */
      float3 hit_k_s_total = { 1.0f, 1.0f, 1.0f }; // default specular reflection coefficient
      int remaining_bounces = max_bounces;

      while (hit && std::isinf(hit->m) && remaining_bounces-- > 0)
      {
        // r(t) = (p + n * epsilon) + t * d
        p = hit->position + hit->normal * epsilon;
        d = normalize(-2.0f * hit->normal * dot(hit->normal, d) + d); // -> like r with phong shading

        // Multiply hit->k_s recursively to the total color
        hit_k_s_total = { hit_k_s_total.x * hit->k_s.x, hit_k_s_total.y * hit->k_s.y, hit_k_s_total.z * hit->k_s.z };

        hit = scene.findClosestHit(p, d);
      }

      float3 color = background_color;

      // Color shading if the ray hits something and it is not a mirror
      if (hit && (!std::isinf(hit->m) || remaining_bounces == 0))
        color = shade(p, d, hit.value(), scene, lights, num_lights, spotlights, num_spotlights);
      // else: set background color

      // Weigh final color with total specular reflection coefficient
      // Only needed if there was a hit to a mirror. Otherwise, a multiplication with { 1, 1, 1 } is done -> no effect
      color = { hit_k_s_total.x * color.x, hit_k_s_total.y * color.y, hit_k_s_total.z * color.z };

      framebuffer(x, y) = color;
    }
  }
}
