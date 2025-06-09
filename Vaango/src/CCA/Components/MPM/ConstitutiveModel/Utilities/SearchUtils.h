#ifndef VAANGO_MPM_CONSTITUTIVE_MODEL_SEARCH_UTIL_H
#define VAANGO_MPM_CONSTITUTIVE_MODEL_SEARCH_UTIL_H

#include <CCA/Components/MPM/ConstitutiveModel/Utilities/nanoflann.hpp>
#include <memory>

//#define VAANGO_UTIL_USE_POLYLINE_POINTERS

namespace Vaango {

namespace Util {

/**
 * A point cloud for kd-tree seraches
 */
#ifdef VAANGO_UTIL_USE_POLYLINE_POINTERS
  struct PolylinePointCloud
  {
    const std::vector<Uintah::Point>* pts;

    PolylinePointCloud(const std::vector<Uintah::Point>* polyline)
    {
      pts = polyline;
    }

    inline size_t kdtree_get_point_count() const { return pts->size(); }

    inline double kdtree_get_pt(const size_t idx, const size_t dim) const
    {
      if (dim == 0) return (*pts)[idx].x();
      else if (dim == 1) return (*pts)[idx].y();
      else return (*pts)[idx].z();
    }

    template <class BBOX>
    bool kdtree_get_bbox(BBOX& bb) const { return false; }
  };
#else
  struct PolylinePointCloud
  {
    std::vector<Uintah::Point> pts;

    PolylinePointCloud(const std::vector<Uintah::Point>& polyline)
    {
      pts = polyline;
    }

    inline size_t kdtree_get_point_count() const { return pts.size(); }

    inline double kdtree_get_pt(const size_t idx, const size_t dim) const
    {
      if (dim == 0) return pts[idx].x();
      else if (dim == 1) return pts[idx].y();
      else return pts[idx].z();
    }

    template <class BBOX>
    bool kdtree_get_bbox([[maybe_unused]] BBOX& bb) const { return false; }
  };
#endif

using PolylineKDTree  = nanoflann::KDTreeSingleIndexAdaptor<
                          nanoflann::L2_Simple_Adaptor<double, PolylinePointCloud>,
                          PolylinePointCloud, 2 /* dim */>;
using PolylineKDTreeP = std::shared_ptr<PolylineKDTree>;


} // end namespace Vaango::Util
} // end namespace Vaango

#endif // VAANGO_MPM_CONSTITUTIVE_MODEL_SEARCH_UTIL_H
