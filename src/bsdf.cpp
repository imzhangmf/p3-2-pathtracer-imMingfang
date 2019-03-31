#include "bsdf.h"

#include <iostream>
#include <algorithm>
#include <utility>

using std::min;
using std::max;
using std::swap;

namespace CGL {

void make_coord_space(Matrix3x3& o2w, const Vector3D& n) {
  Vector3D z = Vector3D(n.x, n.y, n.z);
  Vector3D h = z;
  if (fabs(h.x) <= fabs(h.y) && fabs(h.x) <= fabs(h.z)) h.x = 1.0;
  else if (fabs(h.y) <= fabs(h.x) && fabs(h.y) <= fabs(h.z)) h.y = 1.0;
  else h.z = 1.0;

  z.normalize();
  Vector3D y = cross(h, z);
  y.normalize();
  Vector3D x = cross(z, y);
  x.normalize();

  o2w[0] = x;
  o2w[1] = y;
  o2w[2] = z;
}

// Mirror BSDF //

Spectrum MirrorBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum MirrorBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  // TODO: 1.2
  // Using BSDF::reflect(), implement sample_f for a mirror surface

  reflect(wo, wi);
  *pdf = 1.0;
  return reflectance / abs_cos_theta(*wi);
}

// Microfacet BSDF //

double MicrofacetBSDF::G(const Vector3D& wo, const Vector3D& wi) {
    return 1.0 / (1.0 + Lambda(wi) + Lambda(wo));
}

double MicrofacetBSDF::D(const Vector3D& h) {
  // TODO: 2.2
  // Compute Beckmann normal distribution function (NDF) here.
  // You will need the roughness alpha.

  // alpha 是平面的粗糙程度，再输入wi wo得到平分线向量h，得到h方向上微观镜面的比例

  double theta = acos(h.z);
  return exp(-pow(tan(theta)/alpha, 2.0))/(PI*pow(alpha, 2.0)*pow(h.z, 4.0));
}

Spectrum MicrofacetBSDF::F(const Vector3D& wi) {
  // TODO: 2.3
  // Compute Fresnel term for reflection on dielectric-conductor interface.
  // You will need both eta and etaK, both of which are Spectrum.

  Spectrum cosSquare = Spectrum(wi.z*wi.z, wi.z*wi.z, wi.z*wi.z);
  Spectrum Rs = (eta*eta + k*k - 2.0*eta*wi.z + cosSquare)/(eta*eta + k*k + 2.0*eta*wi.z + cosSquare);
  Spectrum Rp = ((eta*eta + k*k)*cosSquare - 2.0*eta*wi.z + Spectrum(1.0, 1.0, 1.0))/((eta*eta + k*k)*cosSquare + 2.0*eta*wi.z + Spectrum(1.0, 1.0, 1.0));
  return (Rs + Rp)/2.0;
}

Spectrum MicrofacetBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  // TODO: 2.1
  // Implement microfacet model here

  if (wi.z <= 0 || wo.z <= 0) return Spectrum();

  Vector3D h = (wo+wi)/(wo+wi).norm();
  Spectrum result = (F(wi)*G(wo, wi)*D(h))/(4*wo.z*wi.z);

  return result;
}

Spectrum MicrofacetBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  // TODO: 2.4
  // *Importance* sample Beckmann normal distribution function (NDF) here.
  // Note: You should fill in the sampled direction *wi and the corresponding *pdf,
  //       and return the sampled BRDF value.

  Vector2D random = sampler.get_sample();
  double r1 = random.x;
  double r2 = random.y;

  double theta = atan(sqrt(-pow(alpha, 2.0)*log(1.0 - r1)));
  double phi = 2.0*PI*r2;
  Vector3D h = Vector3D(sin(theta)*sin(phi), sin(theta)*cos(phi), cos(theta));
  *wi = 2*h*dot(h, wo) - wo;

  if (wi->z <= 0) {
    *pdf = 0;
    return Spectrum();
  }

  double pTheta = 2.0*sin(theta)/(pow(alpha, 2.0)*pow(cos(theta), 3.0))/exp(pow(tan(theta)/alpha, 2.0));
  double pPhi  = 1/(2.0*PI);
  double pw_h = (pTheta * pPhi)/sin(theta);
  double pw_wi = pw_h/(4.0*dot(*wi, h));
  *pdf = pw_wi;
  return MicrofacetBSDF::f(wo, *wi);
}

// Refraction BSDF //

Spectrum RefractionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum RefractionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  return Spectrum();
}

// Glass BSDF //

Spectrum GlassBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum GlassBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

  // TODO: 1.4
  // Compute Fresnel coefficient and either reflect or refract based on it.

  if (!refract(wo, wi, ior)) {
    reflect(wo, wi);
    *pdf = 1.0;
    return reflectance / abs_cos_theta(*wi);
  } 
  else {
    double R0 = (1.0 - ior)*(1.0 - ior)/(1.0 + ior)/(1.0 + ior);
    double R = R0 + (1.0 - R0)*pow(1.0 - abs_cos_theta(wo), 5.0);
    if (coin_flip(R)) {
      reflect(wo, wi);
      *pdf = R;
      return R * reflectance / abs_cos_theta(*wi);
    } 
    else {
      refract(wo, wi, ior);
      *pdf = 1.0 - R;
      if (wo.z >= 0)
        return (1.0 - R) * transmittance / abs_cos_theta(*wi) / (1.0/ior/ior);
      else
        return (1.0 - R) * transmittance / abs_cos_theta(*wi) / (ior*ior);
    }
  }
}

void BSDF::reflect(const Vector3D& wo, Vector3D* wi) {

  // TODO: 1.1
  // Implement reflection of wo about normal (0,0,1) and store result in wi.
  *wi = Vector3D(-wo.x, -wo.y, wo.z);

}

bool BSDF::refract(const Vector3D& wo, Vector3D* wi, float ior) {

  // TODO: 1.3
  // Use Snell's Law to refract wo surface and store result ray in wi.
  // Return false if refraction does not occur due to total internal reflection
  // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
  // ray entering the surface through vacuum.

  double nl, sign;
  if (wo.z >= 0) {
    nl = 1.0/ior;
    sign = -1.0;
  } else {
    nl = ior;
    sign = 1.0;
  }

  double dis = 1.0 - nl*nl*(1.0 - wo.z*wo.z);
  if (dis < 0) return false;
  else {
    *wi = Vector3D(-nl*wo.x, -nl*wo.y, sign*sqrt(dis));
    return true;
  }
}

// Emission BSDF //

Spectrum EmissionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum EmissionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  *pdf = 1.0 / PI;
  *wi  = sampler.get_sample(pdf);
  return Spectrum();
}

} // namespace CGL
