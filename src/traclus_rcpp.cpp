#include <Rcpp.h>
static const double ZERO_THRESHOLD = 1e-15;
using namespace Rcpp;

// perpendicular distance
static double d_perp(double six, double siy, double eix, double eiy,
                     double sjx, double sjy, double ejx, double ejy) {
  double dx   = eix - six, dy = eiy - siy;
  double len2 = dx*dx + dy*dy;

  if (len2 < ZERO_THRESHOLD) {
    double ddx = sjx-six, ddy = sjy-siy;
    return std::sqrt(ddx*ddx + ddy*ddy);
  }

  double us  = ((sjx-six)*dx + (sjy-siy)*dy) / len2;
  double ue  = ((ejx-six)*dx + (ejy-siy)*dy) / len2;
  double psx = six+us*dx, psy = siy+us*dy;
  double pex = six+ue*dx, pey = siy+ue*dy;

  double l1sq = (sjx-psx)*(sjx-psx) + (sjy-psy)*(sjy-psy);
  double l2sq = (ejx-pex)*(ejx-pex) + (ejy-pey)*(ejy-pey);
  double den  = std::sqrt(l1sq) + std::sqrt(l2sq);
  if (den < ZERO_THRESHOLD) return 0.0;
  return (l1sq + l2sq) / den;
}

// parallel distance
static double d_par(double six, double siy, double eix, double eiy,
                    double sjx, double sjy, double ejx, double ejy) {
  double dx   = eix-six, dy = eiy-siy;
  double len2 = dx*dx + dy*dy;
  if (len2 < ZERO_THRESHOLD) return 0.0;

  double us  = ((sjx-six)*dx + (sjy-siy)*dy) / len2;
  double ue  = ((ejx-six)*dx + (ejy-siy)*dy) / len2;
  double psx = six+us*dx, psy = siy+us*dy;
  double pex = six+ue*dx, pey = siy+ue*dy;

  auto dist2 = [](double ax,double ay,double bx,double by){
    return std::sqrt((ax-bx)*(ax-bx)+(ay-by)*(ay-by));
  };

  double l1 = std::min(dist2(psx,psy,six,siy), dist2(psx,psy,eix,eiy));
  double l2 = std::min(dist2(pex,pey,six,siy), dist2(pex,pey,eix,eiy));
  return std::min(l1, l2);
}

// angle distance
static double d_angle(double six, double siy, double eix, double eiy,
                      double sjx, double sjy, double ejx, double ejy) {
  double lix = eix-six, liy = eiy-siy;
  double ljx = ejx-sjx, ljy = ejy-sjy;
  double li_len = std::sqrt(lix*lix+liy*liy);
  double lj_len = std::sqrt(ljx*ljx+ljy*ljy);
  if (li_len < ZERO_THRESHOLD || lj_len < ZERO_THRESHOLD) return 0.0;

  double cos_t = (lix*ljx+liy*ljy)/(li_len*lj_len);
  cos_t = std::max(-1.0, std::min(1.0, cos_t));
  double theta = std::acos(cos_t);
  return (theta < M_PI/2.0) ? lj_len*std::sin(theta) : lj_len;
}
// (weighted) total distance
static double traclus_dist(double six, double siy, double eix, double eiy,
                           double sjx, double sjy, double ejx, double ejy,
                           double w_perp, double w_par, double w_angle) {
  // Li has to be the longest distance
  double len_i = (eix-six)*(eix-six) + (eiy-siy)*(eiy-siy);
  double len_j = (ejx-sjx)*(ejx-sjx) + (ejy-sjy)*(ejy-sjy);
  if (len_j > len_i) {
    std::swap(six, sjx); std::swap(siy, sjy);
    std::swap(eix, ejx); std::swap(eiy, ejy);
  }
  return w_perp  * d_perp (six,siy,eix,eiy,sjx,sjy,ejx,ejy)
    + w_par   * d_par  (six,siy,eix,eiy,sjx,sjy,ejx,ejy)
    + w_angle * d_angle(six,siy,eix,eiy,sjx,sjy,ejx,ejy);
}

// Compute epsilon-neighbourhoods for all line segments
//
// Returns a sparse neighbourhood list. Each element of the returned list contains the 1-based
// indices of all segments within eps of segment i, including i itself.
 // [[Rcpp::export]]
 List compute_neighbourhoods(
     NumericVector sx, NumericVector sy,
     NumericVector ex, NumericVector ey,
     double eps,
     double w_perp  = 1.0,
     double w_par   = 1.0,
     double w_angle = 1.0) {

   int n = sx.size();
   // Pre-allocate neighbour lists — each segment is at least its own neighbour
   std::vector<std::vector<int>> nbrs(n);
   for (int i = 0; i < n; i++) nbrs[i].push_back(i + 1);  // 1-based

   for (int i = 0; i < n - 1; i++) {
     for (int j = i + 1; j < n; j++) {
       double d = traclus_dist(sx[i],sy[i],ex[i],ey[i],
                               sx[j],sy[j],ex[j],ey[j],
                                                   w_perp, w_par, w_angle);
       if (d <= eps) {
         nbrs[i].push_back(j + 1);  // 1-based for R
         nbrs[j].push_back(i + 1);
       }
     }
   }

   List result(n);
   for (int i = 0; i < n; i++) {
     result[i] = IntegerVector(nbrs[i].begin(), nbrs[i].end());
   }
   return result;
 }
