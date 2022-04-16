//
//  Academic License - for use in teaching, academic research, and meeting
//  course requirements at degree granting institutions only.  Not for
//  government, commercial, or other organizational use.
//
//  c_butterworth.cpp
//
//  Code generation for function 'c_butterworth'
//


// Include files
#include "c_butterworth.h"
#include "rt_nonfinite.h"
#include <cmath>
#include <cstring>

// Variable Definitions
static bool isInitialized_c_butterworth = false;

// Function Declarations
static void b_sqrt(creal_T *x);
static double b_xnrm2(int n, const double x[3]);
static void b_xzlartg(const creal_T f, const creal_T g, double *cs, creal_T *sn);
static int eml_dlahqr(double h[36]);
static void filter(const double b[7], const double a[7], const coder::array<
                   double, 1U> &x, coder::array<double, 1U> &y);
static void flip(coder::array<double, 1U> &x);
static double rt_hypotd_snf(double u0, double u1);
static void xdlanv2(double *a, double *b, double *c, double *d, double *rt1r,
                    double *rt1i, double *rt2r, double *rt2i, double *cs, double
                    *sn);
static void xgehrd(double a[36]);
static void xgetrf(double A[36], int ipiv[6], int *info);
static double xnrm2(int n, const double x[36], int ix0);
static void xzgeev(const double A[36], int *info, creal_T alpha1[6], creal_T
                   beta1[6]);
static void xzhgeqz(const creal_T A[36], int ilo, int ihi, int *info, creal_T
                    alpha1[6], creal_T beta1[6]);
static void xzlartg(const creal_T f, const creal_T g, double *cs, creal_T *sn,
                    creal_T *r);

// Function Definitions
static void b_sqrt(creal_T *x)
{
  double xr;
  double xi;
  double absxi;
  double absxr;
  xr = x->re;
  xi = x->im;
  if (xi == 0.0) {
    if (xr < 0.0) {
      absxi = 0.0;
      xr = std::sqrt(-xr);
    } else {
      absxi = std::sqrt(xr);
      xr = 0.0;
    }
  } else if (xr == 0.0) {
    if (xi < 0.0) {
      absxi = std::sqrt(-xi / 2.0);
      xr = -absxi;
    } else {
      absxi = std::sqrt(xi / 2.0);
      xr = absxi;
    }
  } else if (std::isnan(xr)) {
    absxi = xr;
  } else if (std::isnan(xi)) {
    absxi = xi;
    xr = xi;
  } else if (std::isinf(xi)) {
    absxi = std::abs(xi);
    xr = xi;
  } else if (std::isinf(xr)) {
    if (xr < 0.0) {
      absxi = 0.0;
      xr = xi * -xr;
    } else {
      absxi = xr;
      xr = 0.0;
    }
  } else {
    absxr = std::abs(xr);
    absxi = std::abs(xi);
    if ((absxr > 4.4942328371557893E+307) || (absxi > 4.4942328371557893E+307))
    {
      absxr *= 0.5;
      absxi = rt_hypotd_snf(absxr, absxi * 0.5);
      if (absxi > absxr) {
        absxi = std::sqrt(absxi) * std::sqrt(absxr / absxi + 1.0);
      } else {
        absxi = std::sqrt(absxi) * 1.4142135623730951;
      }
    } else {
      absxi = std::sqrt((rt_hypotd_snf(absxr, absxi) + absxr) * 0.5);
    }

    if (xr > 0.0) {
      xr = 0.5 * (xi / absxi);
    } else {
      if (xi < 0.0) {
        xr = -absxi;
      } else {
        xr = absxi;
      }

      absxi = 0.5 * (xi / xr);
    }
  }

  x->re = absxi;
  x->im = xr;
}

static double b_xnrm2(int n, const double x[3])
{
  double y;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = std::abs(x[1]);
    } else {
      double scale;
      double absxk;
      double t;
      scale = 3.3121686421112381E-170;
      absxk = std::abs(x[1]);
      if (absxk > 3.3121686421112381E-170) {
        y = 1.0;
        scale = absxk;
      } else {
        t = absxk / 3.3121686421112381E-170;
        y = t * t;
      }

      absxk = std::abs(x[2]);
      if (absxk > scale) {
        t = scale / absxk;
        y = y * t * t + 1.0;
        scale = absxk;
      } else {
        t = absxk / scale;
        y += t * t;
      }

      y = scale * std::sqrt(y);
    }
  }

  return y;
}

static void b_xzlartg(const creal_T f, const creal_T g, double *cs, creal_T *sn)
{
  double scale_tmp;
  double f2;
  double scale;
  double fs_re;
  double fs_im;
  double gs_re;
  double gs_im;
  bool guard1 = false;
  double g2;
  scale_tmp = std::abs(f.re);
  f2 = std::abs(f.im);
  if (f2 > scale_tmp) {
    scale_tmp = f2;
  }

  f2 = std::abs(g.re);
  scale = std::abs(g.im);
  if (scale > f2) {
    f2 = scale;
  }

  scale = scale_tmp;
  if (f2 > scale_tmp) {
    scale = f2;
  }

  fs_re = f.re;
  fs_im = f.im;
  gs_re = g.re;
  gs_im = g.im;
  guard1 = false;
  if (scale >= 7.4428285367870146E+137) {
    do {
      fs_re *= 1.3435752215134178E-138;
      fs_im *= 1.3435752215134178E-138;
      gs_re *= 1.3435752215134178E-138;
      gs_im *= 1.3435752215134178E-138;
      scale *= 1.3435752215134178E-138;
    } while (!(scale < 7.4428285367870146E+137));

    guard1 = true;
  } else if (scale <= 1.3435752215134178E-138) {
    if ((g.re == 0.0) && (g.im == 0.0)) {
      *cs = 1.0;
      sn->re = 0.0;
      sn->im = 0.0;
    } else {
      do {
        fs_re *= 7.4428285367870146E+137;
        fs_im *= 7.4428285367870146E+137;
        gs_re *= 7.4428285367870146E+137;
        gs_im *= 7.4428285367870146E+137;
        scale *= 7.4428285367870146E+137;
      } while (!(scale > 1.3435752215134178E-138));

      guard1 = true;
    }
  } else {
    guard1 = true;
  }

  if (guard1) {
    f2 = fs_re * fs_re + fs_im * fs_im;
    g2 = gs_re * gs_re + gs_im * gs_im;
    scale = g2;
    if (1.0 > g2) {
      scale = 1.0;
    }

    if (f2 <= scale * 2.0041683600089728E-292) {
      if ((f.re == 0.0) && (f.im == 0.0)) {
        *cs = 0.0;
        g2 = rt_hypotd_snf(gs_re, gs_im);
        sn->re = gs_re / g2;
        sn->im = -gs_im / g2;
      } else {
        double g2s;
        g2s = std::sqrt(g2);
        *cs = rt_hypotd_snf(fs_re, fs_im) / g2s;
        if (scale_tmp > 1.0) {
          g2 = rt_hypotd_snf(f.re, f.im);
          fs_re = f.re / g2;
          fs_im = f.im / g2;
        } else {
          f2 = 7.4428285367870146E+137 * f.re;
          scale = 7.4428285367870146E+137 * f.im;
          g2 = rt_hypotd_snf(f2, scale);
          fs_re = f2 / g2;
          fs_im = scale / g2;
        }

        gs_re /= g2s;
        gs_im = -gs_im / g2s;
        sn->re = fs_re * gs_re - fs_im * gs_im;
        sn->im = fs_re * gs_im + fs_im * gs_re;
      }
    } else {
      scale = std::sqrt(g2 / f2 + 1.0);
      *cs = 1.0 / scale;
      g2 += f2;
      fs_re = scale * fs_re / g2;
      fs_im = scale * fs_im / g2;
      sn->re = fs_re * gs_re - fs_im * -gs_im;
      sn->im = fs_re * -gs_im + fs_im * gs_re;
    }
  }
}

static int eml_dlahqr(double h[36])
{
  int info;
  double v[3];
  int i;
  bool exitg1;
  double s;
  double ba;
  int knt;
  double d;
  double tst;
  double bb;
  double ab;
  double aa;
  double h22;
  double rt1r;
  int iy;
  info = 0;
  v[0] = 0.0;
  h[2] = 0.0;
  h[3] = 0.0;
  v[1] = 0.0;
  h[9] = 0.0;
  h[10] = 0.0;
  v[2] = 0.0;
  h[16] = 0.0;
  h[17] = 0.0;
  h[23] = 0.0;
  i = 5;
  exitg1 = false;
  while ((!exitg1) && (i + 1 >= 1)) {
    int L;
    bool goto150;
    int its;
    bool exitg2;
    int k;
    int b_i;
    int hoffset;
    int nr;
    L = 1;
    goto150 = false;
    its = 0;
    exitg2 = false;
    while ((!exitg2) && (its < 301)) {
      bool exitg3;
      k = i;
      exitg3 = false;
      while ((!exitg3) && (k + 1 > L)) {
        b_i = k + 6 * (k - 1);
        ba = std::abs(h[b_i]);
        if (ba <= 6.0125050800269183E-292) {
          exitg3 = true;
        } else {
          knt = k + 6 * k;
          bb = std::abs(h[knt]);
          hoffset = b_i - 1;
          tst = std::abs(h[hoffset]) + bb;
          if (tst == 0.0) {
            if (k - 1 >= 1) {
              tst = std::abs(h[(k + 6 * (k - 2)) - 1]);
            }

            if (k + 2 <= 6) {
              tst += std::abs(h[knt + 1]);
            }
          }

          if (ba <= 2.2204460492503131E-16 * tst) {
            tst = std::abs(h[knt - 1]);
            if (ba > tst) {
              ab = ba;
              ba = tst;
            } else {
              ab = tst;
            }

            tst = std::abs(h[hoffset] - h[knt]);
            if (bb > tst) {
              aa = bb;
              bb = tst;
            } else {
              aa = tst;
            }

            s = aa + ab;
            if (ba * (ab / s) <= std::fmax(6.0125050800269183E-292,
                 2.2204460492503131E-16 * (bb * (aa / s)))) {
              exitg3 = true;
            } else {
              k--;
            }
          } else {
            k--;
          }
        }
      }

      L = k + 1;
      if (k + 1 > 1) {
        h[k + 6 * (k - 1)] = 0.0;
      }

      if (k + 1 >= i) {
        goto150 = true;
        exitg2 = true;
      } else {
        int m;
        if (its == 10) {
          hoffset = k + 6 * k;
          s = std::abs(h[hoffset + 1]) + std::abs(h[(k + 6 * (k + 1)) + 2]);
          tst = 0.75 * s + h[hoffset];
          aa = -0.4375 * s;
          ab = s;
          h22 = tst;
        } else if (its == 20) {
          s = std::abs(h[i + 6 * (i - 1)]) + std::abs(h[(i + 6 * (i - 2)) - 1]);
          tst = 0.75 * s + h[i + 6 * i];
          aa = -0.4375 * s;
          ab = s;
          h22 = tst;
        } else {
          hoffset = i + 6 * (i - 1);
          tst = h[hoffset - 1];
          ab = h[hoffset];
          aa = h[(i + 6 * i) - 1];
          h22 = h[i + 6 * i];
        }

        s = ((std::abs(tst) + std::abs(aa)) + std::abs(ab)) + std::abs(h22);
        if (s == 0.0) {
          rt1r = 0.0;
          tst = 0.0;
          bb = 0.0;
          ab = 0.0;
        } else {
          tst /= s;
          ab /= s;
          aa /= s;
          h22 /= s;
          ba = (tst + h22) / 2.0;
          tst = (tst - ba) * (h22 - ba) - aa * ab;
          ab = std::sqrt(std::abs(tst));
          if (tst >= 0.0) {
            rt1r = ba * s;
            bb = rt1r;
            tst = ab * s;
            ab = -tst;
          } else {
            rt1r = ba + ab;
            bb = ba - ab;
            if (std::abs(rt1r - h22) <= std::abs(bb - h22)) {
              rt1r *= s;
              bb = rt1r;
            } else {
              bb *= s;
              rt1r = bb;
            }

            tst = 0.0;
            ab = 0.0;
          }
        }

        m = i - 1;
        exitg3 = false;
        while ((!exitg3) && (m >= k + 1)) {
          hoffset = m + 6 * (m - 1);
          knt = hoffset - 1;
          aa = h[knt] - bb;
          s = (std::abs(aa) + std::abs(ab)) + std::abs(h[hoffset]);
          ba = h[hoffset] / s;
          hoffset = m + 6 * m;
          v[0] = (ba * h[hoffset - 1] + (h[knt] - rt1r) * (aa / s)) - tst * (ab /
            s);
          v[1] = ba * (((h[knt] + h[hoffset]) - rt1r) - bb);
          v[2] = ba * h[hoffset + 1];
          s = (std::abs(v[0]) + std::abs(v[1])) + std::abs(v[2]);
          v[0] /= s;
          v[1] /= s;
          v[2] /= s;
          if (m == k + 1) {
            exitg3 = true;
          } else {
            b_i = m + 6 * (m - 2);
            if (std::abs(h[b_i - 1]) * (std::abs(v[1]) + std::abs(v[2])) <=
                2.2204460492503131E-16 * std::abs(v[0]) * ((std::abs(h[b_i - 2])
                  + std::abs(h[knt])) + std::abs(h[hoffset]))) {
              exitg3 = true;
            } else {
              m--;
            }
          }
        }

        for (int b_k = m; b_k <= i; b_k++) {
          int j;
          nr = (i - b_k) + 2;
          if (3 < nr) {
            nr = 3;
          }

          if (b_k > m) {
            hoffset = (b_k + 6 * (b_k - 2)) - 1;
            for (j = 0; j < nr; j++) {
              v[j] = h[j + hoffset];
            }
          }

          ab = v[0];
          ba = 0.0;
          if (nr > 0) {
            tst = b_xnrm2(nr - 1, v);
            if (tst != 0.0) {
              aa = rt_hypotd_snf(v[0], tst);
              if (v[0] >= 0.0) {
                aa = -aa;
              }

              if (std::abs(aa) < 1.0020841800044864E-292) {
                knt = -1;
                do {
                  knt++;
                  for (iy = 2; iy <= nr; iy++) {
                    v[iy - 1] *= 9.9792015476736E+291;
                  }

                  aa *= 9.9792015476736E+291;
                  ab *= 9.9792015476736E+291;
                } while (!(std::abs(aa) >= 1.0020841800044864E-292));

                aa = rt_hypotd_snf(ab, b_xnrm2(nr - 1, v));
                if (ab >= 0.0) {
                  aa = -aa;
                }

                ba = (aa - ab) / aa;
                tst = 1.0 / (ab - aa);
                for (iy = 2; iy <= nr; iy++) {
                  v[iy - 1] *= tst;
                }

                for (iy = 0; iy <= knt; iy++) {
                  aa *= 1.0020841800044864E-292;
                }

                ab = aa;
              } else {
                ba = (aa - v[0]) / aa;
                tst = 1.0 / (v[0] - aa);
                for (iy = 2; iy <= nr; iy++) {
                  v[iy - 1] *= tst;
                }

                ab = aa;
              }
            }
          }

          v[0] = ab;
          if (b_k > m) {
            h[(b_k + 6 * (b_k - 2)) - 1] = ab;
            b_i = b_k + 6 * (b_k - 2);
            h[b_i] = 0.0;
            if (b_k < i) {
              h[b_i + 1] = 0.0;
            }
          } else {
            if (m > k + 1) {
              h[(b_k + 6 * (b_k - 2)) - 1] *= 1.0 - ba;
            }
          }

          s = v[1];
          tst = ba * v[1];
          if (nr == 3) {
            d = v[2];
            aa = ba * v[2];
            for (j = b_k; j < 7; j++) {
              iy = b_k + 6 * (j - 1);
              hoffset = iy - 1;
              knt = iy + 1;
              ab = (h[hoffset] + s * h[iy]) + d * h[knt];
              h[hoffset] -= ab * ba;
              h[iy] -= ab * tst;
              h[knt] -= ab * aa;
            }

            if (b_k + 3 < i + 1) {
              b_i = b_k + 2;
            } else {
              b_i = i;
            }

            for (j = 0; j <= b_i; j++) {
              iy = j + 6 * (b_k - 1);
              hoffset = j + 6 * b_k;
              knt = j + 6 * (b_k + 1);
              ab = (h[iy] + s * h[hoffset]) + d * h[knt];
              h[iy] -= ab * ba;
              h[hoffset] -= ab * tst;
              h[knt] -= ab * aa;
            }
          } else {
            if (nr == 2) {
              for (j = b_k; j < 7; j++) {
                iy = b_k + 6 * (j - 1);
                hoffset = iy - 1;
                ab = h[hoffset] + s * h[iy];
                h[hoffset] -= ab * ba;
                h[iy] -= ab * tst;
              }

              for (j = 0; j <= i; j++) {
                iy = j + 6 * (b_k - 1);
                hoffset = j + 6 * b_k;
                ab = h[iy] + s * h[hoffset];
                h[iy] -= ab * ba;
                h[hoffset] -= ab * tst;
              }
            }
          }
        }

        its++;
      }
    }

    if (!goto150) {
      info = i + 1;
      exitg1 = true;
    } else {
      if ((L != i + 1) && (L == i)) {
        b_i = i + 6 * i;
        hoffset = b_i - 1;
        s = h[hoffset];
        nr = 6 * (i - 1);
        knt = i + nr;
        d = h[knt];
        tst = h[b_i];
        xdlanv2(&h[(i + 6 * (i - 1)) - 1], &s, &d, &tst, &ab, &aa, &ba, &bb,
                &h22, &rt1r);
        h[hoffset] = s;
        h[knt] = d;
        h[b_i] = tst;
        if (6 > i + 1) {
          hoffset = 4 - i;
          iy = i + (i + 1) * 6;
          knt = iy - 1;
          for (k = 0; k <= hoffset; k++) {
            tst = h22 * h[knt] + rt1r * h[iy];
            h[iy] = h22 * h[iy] - rt1r * h[knt];
            h[knt] = tst;
            iy += 6;
            knt += 6;
          }
        }

        if (i - 1 >= 1) {
          iy = i * 6;
          for (k = 0; k <= i - 2; k++) {
            tst = h22 * h[nr] + rt1r * h[iy];
            h[iy] = h22 * h[iy] - rt1r * h[nr];
            h[nr] = tst;
            iy++;
            nr++;
          }
        }
      }

      i = L - 2;
    }
  }

  return info;
}

static void filter(const double b[7], const double a[7], const coder::array<
                   double, 1U> &x, coder::array<double, 1U> &y)
{
  int nx;
  int loop_ub;
  int naxpy;
  nx = x.size(0) - 1;
  loop_ub = x.size(0);
  y.set_size(x.size(0));
  for (naxpy = 0; naxpy < loop_ub; naxpy++) {
    y[naxpy] = 0.0;
  }

  for (int k = 0; k <= nx; k++) {
    int j;
    int y_tmp;
    double as;
    loop_ub = nx - k;
    naxpy = loop_ub + 1;
    if (naxpy >= 7) {
      naxpy = 7;
    }

    for (j = 0; j < naxpy; j++) {
      y_tmp = k + j;
      y[y_tmp] = y[y_tmp] + x[k] * b[j];
    }

    if (loop_ub < 6) {
      naxpy = loop_ub;
    } else {
      naxpy = 6;
    }

    as = -y[k];
    for (j = 0; j < naxpy; j++) {
      y_tmp = (k + j) + 1;
      y[y_tmp] = y[y_tmp] + as * a[j + 1];
    }
  }
}

static void flip(coder::array<double, 1U> &x)
{
  int dim;
  int vstride;
  dim = 2;
  if (x.size(0) != 1) {
    dim = 1;
  }

  if (x.size(0) != 0) {
    if (dim <= 1) {
      vstride = x.size(0);
    } else {
      vstride = 1;
    }

    if (vstride > 1) {
      int k;
      int nd2;
      int i;
      vstride = 1;
      for (k = 0; k <= dim - 2; k++) {
        vstride *= x.size(0);
      }

      if (dim <= 1) {
        dim = x.size(0) - 1;
      } else {
        dim = 0;
      }

      nd2 = (dim + 1) >> 1;
      i = vstride - 1;
      for (int b_i = 0; b_i <= i; b_i++) {
        for (k = 0; k < nd2; k++) {
          int tmp_tmp;
          double tmp;
          tmp_tmp = b_i + k * vstride;
          tmp = x[tmp_tmp];
          x[tmp_tmp] = x[b_i + (dim - k) * vstride];
          x[b_i + (dim - k) * vstride] = tmp;
        }
      }
    }
  }
}

static double rt_hypotd_snf(double u0, double u1)
{
  double y;
  double a;
  a = std::abs(u0);
  y = std::abs(u1);
  if (a < y) {
    a /= y;
    y *= std::sqrt(a * a + 1.0);
  } else if (a > y) {
    y /= a;
    y = a * std::sqrt(y * y + 1.0);
  } else {
    if (!std::isnan(y)) {
      y = a * 1.4142135623730951;
    }
  }

  return y;
}

static void xdlanv2(double *a, double *b, double *c, double *d, double *rt1r,
                    double *rt1i, double *rt2r, double *rt2i, double *cs, double
                    *sn)
{
  if (*c == 0.0) {
    *cs = 1.0;
    *sn = 0.0;
  } else if (*b == 0.0) {
    double bcmax;
    *cs = 0.0;
    *sn = 1.0;
    bcmax = *d;
    *d = *a;
    *a = bcmax;
    *b = -*c;
    *c = 0.0;
  } else {
    double tau;
    tau = *a - *d;
    if ((tau == 0.0) && ((*b < 0.0) != (*c < 0.0))) {
      *cs = 1.0;
      *sn = 0.0;
    } else {
      double p;
      double bcmax;
      double bcmis;
      double scale;
      int b_b;
      int b_c;
      double z;
      p = 0.5 * tau;
      bcmis = std::abs(*b);
      scale = std::abs(*c);
      bcmax = std::fmax(bcmis, scale);
      if (!(*b < 0.0)) {
        b_b = 1;
      } else {
        b_b = -1;
      }

      if (!(*c < 0.0)) {
        b_c = 1;
      } else {
        b_c = -1;
      }

      bcmis = std::fmin(bcmis, scale) * static_cast<double>(b_b) * static_cast<
        double>(b_c);
      scale = std::fmax(std::abs(p), bcmax);
      z = p / scale * p + bcmax / scale * bcmis;
      if (z >= 8.8817841970012523E-16) {
        *a = std::sqrt(scale) * std::sqrt(z);
        if (p < 0.0) {
          *a = -*a;
        }

        z = p + *a;
        *a = *d + z;
        *d -= bcmax / z * bcmis;
        tau = rt_hypotd_snf(*c, z);
        *cs = z / tau;
        *sn = *c / tau;
        *b -= *c;
        *c = 0.0;
      } else {
        bcmis = *b + *c;
        tau = rt_hypotd_snf(bcmis, tau);
        *cs = std::sqrt(0.5 * (std::abs(bcmis) / tau + 1.0));
        if (!(bcmis < 0.0)) {
          b_b = 1;
        } else {
          b_b = -1;
        }

        *sn = -(p / (tau * *cs)) * static_cast<double>(b_b);
        bcmax = *a * *cs + *b * *sn;
        scale = -*a * *sn + *b * *cs;
        z = *c * *cs + *d * *sn;
        bcmis = -*c * *sn + *d * *cs;
        *b = scale * *cs + bcmis * *sn;
        *c = -bcmax * *sn + z * *cs;
        bcmax = 0.5 * ((bcmax * *cs + z * *sn) + (-scale * *sn + bcmis * *cs));
        *a = bcmax;
        *d = bcmax;
        if (*c != 0.0) {
          if (*b != 0.0) {
            if ((*b < 0.0) == (*c < 0.0)) {
              bcmis = std::sqrt(std::abs(*b));
              z = std::sqrt(std::abs(*c));
              *a = bcmis * z;
              if (!(*c < 0.0)) {
                p = *a;
              } else {
                p = -*a;
              }

              tau = 1.0 / std::sqrt(std::abs(*b + *c));
              *a = bcmax + p;
              *d = bcmax - p;
              *b -= *c;
              *c = 0.0;
              scale = bcmis * tau;
              bcmis = z * tau;
              bcmax = *cs * scale - *sn * bcmis;
              *sn = *cs * bcmis + *sn * scale;
              *cs = bcmax;
            }
          } else {
            *b = -*c;
            *c = 0.0;
            bcmax = *cs;
            *cs = -*sn;
            *sn = bcmax;
          }
        }
      }
    }
  }

  *rt1r = *a;
  *rt2r = *d;
  if (*c == 0.0) {
    *rt1i = 0.0;
    *rt2i = 0.0;
  } else {
    *rt1i = std::sqrt(std::abs(*b)) * std::sqrt(std::abs(*c));
    *rt2i = -*rt1i;
  }
}

static void xgehrd(double a[36])
{
  int i;
  double work[6];
  int b_i;
  double alpha1;
  double tau[5];
  double temp;
  double beta1;
  int k;
  int jA;
  for (i = 0; i < 6; i++) {
    work[i] = 0.0;
  }

  for (i = 0; i < 5; i++) {
    int in;
    int alpha1_tmp;
    int c_i;
    int coltop;
    int iv0_tmp;
    int knt;
    int lastv;
    int lastc;
    bool exitg2;
    int ix;
    int exitg1;
    int d_i;
    b_i = i * 6 + 2;
    in = (i + 1) * 6;
    alpha1_tmp = (i + 6 * i) + 1;
    alpha1 = a[alpha1_tmp];
    if (i + 3 < 6) {
      c_i = i + 1;
    } else {
      c_i = 4;
    }

    coltop = c_i + b_i;
    tau[i] = 0.0;
    temp = xnrm2(4 - i, a, coltop);
    if (temp != 0.0) {
      beta1 = rt_hypotd_snf(a[alpha1_tmp], temp);
      if (a[alpha1_tmp] >= 0.0) {
        beta1 = -beta1;
      }

      if (std::abs(beta1) < 1.0020841800044864E-292) {
        knt = -1;
        c_i = (coltop - i) + 3;
        do {
          knt++;
          for (k = coltop; k <= c_i; k++) {
            a[k - 1] *= 9.9792015476736E+291;
          }

          beta1 *= 9.9792015476736E+291;
          alpha1 *= 9.9792015476736E+291;
        } while (!(std::abs(beta1) >= 1.0020841800044864E-292));

        beta1 = rt_hypotd_snf(alpha1, xnrm2(4 - i, a, coltop));
        if (alpha1 >= 0.0) {
          beta1 = -beta1;
        }

        tau[i] = (beta1 - alpha1) / beta1;
        temp = 1.0 / (alpha1 - beta1);
        c_i = (coltop - i) + 3;
        for (k = coltop; k <= c_i; k++) {
          a[k - 1] *= temp;
        }

        for (k = 0; k <= knt; k++) {
          beta1 *= 1.0020841800044864E-292;
        }

        alpha1 = beta1;
      } else {
        tau[i] = (beta1 - a[alpha1_tmp]) / beta1;
        temp = 1.0 / (a[alpha1_tmp] - beta1);
        c_i = (coltop - i) + 3;
        for (k = coltop; k <= c_i; k++) {
          a[k - 1] *= temp;
        }

        alpha1 = beta1;
      }
    }

    a[alpha1_tmp] = 1.0;
    iv0_tmp = (i + b_i) - 1;
    coltop = in + 1;
    if (tau[i] != 0.0) {
      lastv = 4 - i;
      b_i = (iv0_tmp - i) + 4;
      while ((lastv + 1 > 0) && (a[b_i] == 0.0)) {
        lastv--;
        b_i--;
      }

      lastc = 6;
      exitg2 = false;
      while ((!exitg2) && (lastc > 0)) {
        knt = in + lastc;
        k = knt;
        do {
          exitg1 = 0;
          if (k <= knt + lastv * 6) {
            if (a[k - 1] != 0.0) {
              exitg1 = 1;
            } else {
              k += 6;
            }
          } else {
            lastc--;
            exitg1 = 2;
          }
        } while (exitg1 == 0);

        if (exitg1 == 1) {
          exitg2 = true;
        }
      }
    } else {
      lastv = -1;
      lastc = 0;
    }

    if (lastv + 1 > 0) {
      if (lastc != 0) {
        if (0 <= lastc - 1) {
          std::memset(&work[0], 0, lastc * sizeof(double));
        }

        ix = iv0_tmp;
        c_i = (in + 6 * lastv) + 1;
        for (b_i = coltop; b_i <= c_i; b_i += 6) {
          knt = 0;
          d_i = (b_i + lastc) - 1;
          for (k = b_i; k <= d_i; k++) {
            work[knt] += a[k - 1] * a[ix];
            knt++;
          }

          ix++;
        }
      }

      if (!(-tau[i] == 0.0)) {
        jA = in;
        knt = iv0_tmp;
        for (b_i = 0; b_i <= lastv; b_i++) {
          if (a[knt] != 0.0) {
            temp = a[knt] * -tau[i];
            ix = 0;
            c_i = jA + 1;
            d_i = lastc + jA;
            for (coltop = c_i; coltop <= d_i; coltop++) {
              a[coltop - 1] += work[ix] * temp;
              ix++;
            }
          }

          knt++;
          jA += 6;
        }
      }
    }

    jA = (i + in) + 2;
    if (tau[i] != 0.0) {
      lastv = 5 - i;
      b_i = (iv0_tmp - i) + 4;
      while ((lastv > 0) && (a[b_i] == 0.0)) {
        lastv--;
        b_i--;
      }

      lastc = 4 - i;
      exitg2 = false;
      while ((!exitg2) && (lastc + 1 > 0)) {
        coltop = jA + lastc * 6;
        k = coltop;
        do {
          exitg1 = 0;
          if (k <= (coltop + lastv) - 1) {
            if (a[k - 1] != 0.0) {
              exitg1 = 1;
            } else {
              k++;
            }
          } else {
            lastc--;
            exitg1 = 2;
          }
        } while (exitg1 == 0);

        if (exitg1 == 1) {
          exitg2 = true;
        }
      }
    } else {
      lastv = 0;
      lastc = -1;
    }

    if (lastv > 0) {
      if (lastc + 1 != 0) {
        if (0 <= lastc) {
          std::memset(&work[0], 0, (lastc + 1) * sizeof(double));
        }

        knt = 0;
        c_i = jA + 6 * lastc;
        for (b_i = jA; b_i <= c_i; b_i += 6) {
          ix = iv0_tmp;
          temp = 0.0;
          d_i = (b_i + lastv) - 1;
          for (k = b_i; k <= d_i; k++) {
            temp += a[k - 1] * a[ix];
            ix++;
          }

          work[knt] += temp;
          knt++;
        }
      }

      if (!(-tau[i] == 0.0)) {
        knt = 0;
        for (b_i = 0; b_i <= lastc; b_i++) {
          if (work[knt] != 0.0) {
            temp = work[knt] * -tau[i];
            ix = iv0_tmp;
            c_i = lastv + jA;
            for (coltop = jA; coltop < c_i; coltop++) {
              a[coltop - 1] += a[ix] * temp;
              ix++;
            }
          }

          knt++;
          jA += 6;
        }
      }
    }

    a[alpha1_tmp] = alpha1;
  }
}

static void xgetrf(double A[36], int ipiv[6], int *info)
{
  int i;
  int iy;
  int jA;
  int ix;
  for (i = 0; i < 6; i++) {
    ipiv[i] = i + 1;
  }

  *info = 0;
  for (int j = 0; j < 5; j++) {
    int mmj_tmp;
    int b;
    int jj;
    int jp1j;
    double smax;
    int k;
    mmj_tmp = 4 - j;
    b = j * 7;
    jj = j * 7;
    jp1j = b + 2;
    iy = 6 - j;
    jA = 0;
    ix = b;
    smax = std::abs(A[jj]);
    for (k = 2; k <= iy; k++) {
      double s;
      ix++;
      s = std::abs(A[ix]);
      if (s > smax) {
        jA = k - 1;
        smax = s;
      }
    }

    if (A[jj + jA] != 0.0) {
      if (jA != 0) {
        iy = j + jA;
        ipiv[j] = iy + 1;
        ix = j;
        for (k = 0; k < 6; k++) {
          smax = A[ix];
          A[ix] = A[iy];
          A[iy] = smax;
          ix += 6;
          iy += 6;
        }
      }

      i = (jj - j) + 6;
      for (iy = jp1j; iy <= i; iy++) {
        A[iy - 1] /= A[jj];
      }
    } else {
      *info = j + 1;
    }

    iy = b + 6;
    jA = jj;
    for (b = 0; b <= mmj_tmp; b++) {
      smax = A[iy];
      if (A[iy] != 0.0) {
        ix = jj + 1;
        i = jA + 8;
        k = (jA - j) + 12;
        for (jp1j = i; jp1j <= k; jp1j++) {
          A[jp1j - 1] += A[ix] * -smax;
          ix++;
        }
      }

      iy += 6;
      jA += 6;
    }
  }

  if ((*info == 0) && (!(A[35] != 0.0))) {
    *info = 6;
  }
}

static double xnrm2(int n, const double x[36], int ix0)
{
  double y;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = std::abs(x[ix0 - 1]);
    } else {
      double scale;
      int kend;
      scale = 3.3121686421112381E-170;
      kend = (ix0 + n) - 1;
      for (int k = ix0; k <= kend; k++) {
        double absxk;
        absxk = std::abs(x[k - 1]);
        if (absxk > scale) {
          double t;
          t = scale / absxk;
          y = y * t * t + 1.0;
          scale = absxk;
        } else {
          double t;
          t = absxk / scale;
          y += t * t;
        }
      }

      y = scale * std::sqrt(y);
    }
  }

  return y;
}

static void xzgeev(const double A[36], int *info, creal_T alpha1[6], creal_T
                   beta1[6])
{
  int ii;
  creal_T At[36];
  double anrm;
  int jcol;
  bool exitg1;
  double absxk;
  creal_T atmp;
  for (ii = 0; ii < 36; ii++) {
    At[ii].re = A[ii];
    At[ii].im = 0.0;
  }

  *info = 0;
  anrm = 0.0;
  jcol = 0;
  exitg1 = false;
  while ((!exitg1) && (jcol < 36)) {
    absxk = rt_hypotd_snf(At[jcol].re, At[jcol].im);
    if (std::isnan(absxk)) {
      anrm = rtNaN;
      exitg1 = true;
    } else {
      if (absxk > anrm) {
        anrm = absxk;
      }

      jcol++;
    }
  }

  if (std::isinf(anrm) || std::isnan(anrm)) {
    for (int i = 0; i < 6; i++) {
      alpha1[i].re = rtNaN;
      alpha1[i].im = 0.0;
      beta1[i].re = rtNaN;
      beta1[i].im = 0.0;
    }
  } else {
    bool ilascl;
    int i;
    double anrmto;
    bool guard1 = false;
    int ilo;
    double ctoc;
    int ihi;
    bool notdone;
    int exitg3;
    double stemp_im;
    int j;
    double cto1;
    double a;
    int nzcount;
    bool exitg4;
    int At_tmp;
    ilascl = false;
    anrmto = anrm;
    guard1 = false;
    if ((anrm > 0.0) && (anrm < 6.7178761075670888E-139)) {
      anrmto = 6.7178761075670888E-139;
      ilascl = true;
      guard1 = true;
    } else {
      if (anrm > 1.4885657073574029E+138) {
        anrmto = 1.4885657073574029E+138;
        ilascl = true;
        guard1 = true;
      }
    }

    if (guard1) {
      absxk = anrm;
      ctoc = anrmto;
      notdone = true;
      while (notdone) {
        stemp_im = absxk * 2.0041683600089728E-292;
        cto1 = ctoc / 4.9896007738368E+291;
        if ((stemp_im > ctoc) && (ctoc != 0.0)) {
          a = 2.0041683600089728E-292;
          absxk = stemp_im;
        } else if (cto1 > absxk) {
          a = 4.9896007738368E+291;
          ctoc = cto1;
        } else {
          a = ctoc / absxk;
          notdone = false;
        }

        for (ii = 0; ii < 36; ii++) {
          At[ii].re *= a;
          At[ii].im *= a;
        }
      }
    }

    ilo = 1;
    ihi = 6;
    do {
      exitg3 = 0;
      i = 0;
      j = 0;
      notdone = false;
      ii = ihi;
      exitg1 = false;
      while ((!exitg1) && (ii > 0)) {
        nzcount = 0;
        i = ii;
        j = ihi;
        jcol = 0;
        exitg4 = false;
        while ((!exitg4) && (jcol <= ihi - 1)) {
          At_tmp = (ii + 6 * jcol) - 1;
          if ((At[At_tmp].re != 0.0) || (At[At_tmp].im != 0.0) || (ii == jcol +
               1)) {
            if (nzcount == 0) {
              j = jcol + 1;
              nzcount = 1;
              jcol++;
            } else {
              nzcount = 2;
              exitg4 = true;
            }
          } else {
            jcol++;
          }
        }

        if (nzcount < 2) {
          notdone = true;
          exitg1 = true;
        } else {
          ii--;
        }
      }

      if (!notdone) {
        exitg3 = 2;
      } else {
        if (i != ihi) {
          for (jcol = 0; jcol < 6; jcol++) {
            ii = (i + 6 * jcol) - 1;
            atmp = At[ii];
            At_tmp = (ihi + 6 * jcol) - 1;
            At[ii] = At[At_tmp];
            At[At_tmp] = atmp;
          }
        }

        if (j != ihi) {
          for (jcol = 0; jcol < ihi; jcol++) {
            ii = jcol + 6 * (j - 1);
            atmp = At[ii];
            At_tmp = jcol + 6 * (ihi - 1);
            At[ii] = At[At_tmp];
            At[At_tmp] = atmp;
          }
        }

        ihi--;
        if (ihi == 1) {
          exitg3 = 1;
        }
      }
    } while (exitg3 == 0);

    if (exitg3 != 1) {
      int exitg2;
      do {
        exitg2 = 0;
        i = 0;
        j = 0;
        notdone = false;
        jcol = ilo;
        exitg1 = false;
        while ((!exitg1) && (jcol <= ihi)) {
          nzcount = 0;
          i = ihi;
          j = jcol;
          ii = ilo;
          exitg4 = false;
          while ((!exitg4) && (ii <= ihi)) {
            At_tmp = (ii + 6 * (jcol - 1)) - 1;
            if ((At[At_tmp].re != 0.0) || (At[At_tmp].im != 0.0) || (ii == jcol))
            {
              if (nzcount == 0) {
                i = ii;
                nzcount = 1;
                ii++;
              } else {
                nzcount = 2;
                exitg4 = true;
              }
            } else {
              ii++;
            }
          }

          if (nzcount < 2) {
            notdone = true;
            exitg1 = true;
          } else {
            jcol++;
          }
        }

        if (!notdone) {
          exitg2 = 1;
        } else {
          if (i != ilo) {
            for (jcol = ilo; jcol < 7; jcol++) {
              ii = 6 * (jcol - 1);
              nzcount = (i + ii) - 1;
              atmp = At[nzcount];
              At_tmp = (ilo + ii) - 1;
              At[nzcount] = At[At_tmp];
              At[At_tmp] = atmp;
            }
          }

          if (j != ilo) {
            for (jcol = 0; jcol < ihi; jcol++) {
              ii = jcol + 6 * (j - 1);
              atmp = At[ii];
              At_tmp = jcol + 6 * (ilo - 1);
              At[ii] = At[At_tmp];
              At[At_tmp] = atmp;
            }
          }

          ilo++;
          if (ilo == ihi) {
            exitg2 = 1;
          }
        }
      } while (exitg2 == 0);
    }

    if (ihi >= ilo + 2) {
      for (jcol = ilo - 1; jcol + 1 < ihi - 1; jcol++) {
        int jcolp1;
        jcolp1 = jcol + 2;
        for (int jrow = ihi - 1; jrow + 1 > jcol + 2; jrow--) {
          At_tmp = jrow + 6 * jcol;
          xzlartg(At[At_tmp - 1], At[At_tmp], &absxk, &atmp, &At[(jrow + 6 *
                   jcol) - 1]);
          At[At_tmp].re = 0.0;
          At[At_tmp].im = 0.0;
          for (j = jcolp1; j < 7; j++) {
            ii = jrow + 6 * (j - 1);
            nzcount = ii - 1;
            ctoc = absxk * At[nzcount].re + (atmp.re * At[ii].re - atmp.im *
              At[ii].im);
            stemp_im = absxk * At[nzcount].im + (atmp.re * At[ii].im + atmp.im *
              At[ii].re);
            cto1 = At[nzcount].re;
            At[ii].re = absxk * At[ii].re - (atmp.re * At[nzcount].re + atmp.im *
              At[nzcount].im);
            At[ii].im = absxk * At[ii].im - (atmp.re * At[nzcount].im - atmp.im *
              cto1);
            At[nzcount].re = ctoc;
            At[nzcount].im = stemp_im;
          }

          atmp.re = -atmp.re;
          atmp.im = -atmp.im;
          for (i = 1; i <= ihi; i++) {
            ii = (i + 6 * (jrow - 1)) - 1;
            nzcount = (i + 6 * jrow) - 1;
            ctoc = absxk * At[nzcount].re + (atmp.re * At[ii].re - atmp.im *
              At[ii].im);
            stemp_im = absxk * At[nzcount].im + (atmp.re * At[ii].im + atmp.im *
              At[ii].re);
            cto1 = At[nzcount].re;
            At[ii].re = absxk * At[ii].re - (atmp.re * At[nzcount].re + atmp.im *
              At[nzcount].im);
            At[ii].im = absxk * At[ii].im - (atmp.re * At[nzcount].im - atmp.im *
              cto1);
            At[nzcount].re = ctoc;
            At[nzcount].im = stemp_im;
          }
        }
      }
    }

    xzhgeqz(At, ilo, ihi, info, alpha1, beta1);
    if ((*info == 0) && ilascl) {
      notdone = true;
      while (notdone) {
        stemp_im = anrmto * 2.0041683600089728E-292;
        cto1 = anrm / 4.9896007738368E+291;
        if ((stemp_im > anrm) && (anrm != 0.0)) {
          a = 2.0041683600089728E-292;
          anrmto = stemp_im;
        } else if (cto1 > anrmto) {
          a = 4.9896007738368E+291;
          anrm = cto1;
        } else {
          a = anrm / anrmto;
          notdone = false;
        }

        for (ii = 0; ii < 6; ii++) {
          alpha1[ii].re *= a;
          alpha1[ii].im *= a;
        }
      }
    }
  }
}

static void xzhgeqz(const creal_T A[36], int ilo, int ihi, int *info, creal_T
                    alpha1[6], creal_T beta1[6])
{
  creal_T b_A[36];
  int i;
  double eshift_re;
  double eshift_im;
  creal_T ctemp;
  double anorm;
  double scale;
  double b_atol;
  bool firstNonZero;
  int j;
  int ctemp_tmp;
  double ascale;
  int jp1;
  double temp1;
  bool guard1 = false;
  bool guard2 = false;
  double temp2;
  int ilast;
  creal_T shift;
  creal_T b_ascale;
  std::memcpy(&b_A[0], &A[0], 36U * sizeof(creal_T));
  *info = 0;
  for (i = 0; i < 6; i++) {
    alpha1[i].re = 0.0;
    alpha1[i].im = 0.0;
    beta1[i].re = 1.0;
    beta1[i].im = 0.0;
  }

  eshift_re = 0.0;
  eshift_im = 0.0;
  ctemp.re = 0.0;
  ctemp.im = 0.0;
  anorm = 0.0;
  if (ilo <= ihi) {
    scale = 0.0;
    anorm = 0.0;
    firstNonZero = true;
    for (j = ilo; j <= ihi; j++) {
      ctemp_tmp = j + 1;
      if (ihi < j + 1) {
        ctemp_tmp = ihi;
      }

      for (i = ilo; i <= ctemp_tmp; i++) {
        jp1 = (i + 6 * (j - 1)) - 1;
        if (A[jp1].re != 0.0) {
          temp1 = std::abs(A[jp1].re);
          if (firstNonZero) {
            anorm = 1.0;
            scale = temp1;
            firstNonZero = false;
          } else if (scale < temp1) {
            temp2 = scale / temp1;
            anorm = anorm * temp2 * temp2 + 1.0;
            scale = temp1;
          } else {
            temp2 = temp1 / scale;
            anorm += temp2 * temp2;
          }
        }

        if (A[jp1].im != 0.0) {
          temp1 = std::abs(A[jp1].im);
          if (firstNonZero) {
            anorm = 1.0;
            scale = temp1;
            firstNonZero = false;
          } else if (scale < temp1) {
            temp2 = scale / temp1;
            anorm = anorm * temp2 * temp2 + 1.0;
            scale = temp1;
          } else {
            temp2 = temp1 / scale;
            anorm += temp2 * temp2;
          }
        }
      }
    }

    anorm = scale * std::sqrt(anorm);
  }

  scale = 2.2204460492503131E-16 * anorm;
  b_atol = 2.2250738585072014E-308;
  if (scale > 2.2250738585072014E-308) {
    b_atol = scale;
  }

  scale = 2.2250738585072014E-308;
  if (anorm > 2.2250738585072014E-308) {
    scale = anorm;
  }

  ascale = 1.0 / scale;
  firstNonZero = true;
  ctemp_tmp = ihi + 1;
  for (j = ctemp_tmp; j < 7; j++) {
    alpha1[j - 1] = A[(j + 6 * (j - 1)) - 1];
  }

  guard1 = false;
  guard2 = false;
  if (ihi >= ilo) {
    int ifirst;
    int istart;
    int ilastm1;
    int ilastm;
    int iiter;
    bool goto60;
    bool goto70;
    bool goto90;
    int jiter;
    ifirst = ilo;
    istart = ilo;
    ilast = ihi - 1;
    ilastm1 = ihi - 2;
    ilastm = ihi;
    iiter = 0;
    goto60 = false;
    goto70 = false;
    goto90 = false;
    jiter = 0;
    int exitg1;
    do {
      exitg1 = 0;
      if (jiter <= 30 * ((ihi - ilo) + 1) - 1) {
        bool b_guard1 = false;
        bool exitg2;
        b_guard1 = false;
        if (ilast + 1 == ilo) {
          goto60 = true;
          b_guard1 = true;
        } else {
          ctemp_tmp = ilast + 6 * ilastm1;
          if (std::abs(b_A[ctemp_tmp].re) + std::abs(b_A[ctemp_tmp].im) <=
              b_atol) {
            b_A[ctemp_tmp].re = 0.0;
            b_A[ctemp_tmp].im = 0.0;
            goto60 = true;
            b_guard1 = true;
          } else {
            bool guard3 = false;
            j = ilastm1;
            guard3 = false;
            exitg2 = false;
            while ((!exitg2) && (j + 1 >= ilo)) {
              if (j + 1 == ilo) {
                guard3 = true;
                exitg2 = true;
              } else {
                ctemp_tmp = j + 6 * (j - 1);
                if (std::abs(b_A[ctemp_tmp].re) + std::abs(b_A[ctemp_tmp].im) <=
                    b_atol) {
                  b_A[ctemp_tmp].re = 0.0;
                  b_A[ctemp_tmp].im = 0.0;
                  guard3 = true;
                  exitg2 = true;
                } else {
                  j--;
                  guard3 = false;
                }
              }
            }

            if (guard3) {
              ifirst = j + 1;
              goto70 = true;
            }

            if (goto70) {
              b_guard1 = true;
            } else {
              for (i = 0; i < 6; i++) {
                alpha1[i].re = rtNaN;
                alpha1[i].im = 0.0;
                beta1[i].re = rtNaN;
                beta1[i].im = 0.0;
              }

              *info = 1;
              exitg1 = 1;
            }
          }
        }

        if (b_guard1) {
          if (goto60) {
            goto60 = false;
            alpha1[ilast] = b_A[ilast + 6 * ilast];
            ilast = ilastm1;
            ilastm1--;
            if (ilast + 1 < ilo) {
              firstNonZero = false;
              guard2 = true;
              exitg1 = 1;
            } else {
              iiter = 0;
              eshift_re = 0.0;
              eshift_im = 0.0;
              ilastm = ilast + 1;
              jiter++;
            }
          } else {
            if (goto70) {
              double ad22_re;
              double ad22_im;
              goto70 = false;
              iiter++;
              if (iiter - iiter / 10 * 10 != 0) {
                double ascale_re;
                double ascale_im;
                double t1_re;
                double t1_im;
                double t1_im_tmp;
                jp1 = ilastm1 + 6 * ilastm1;
                anorm = ascale * b_A[jp1].re;
                scale = ascale * b_A[jp1].im;
                if (scale == 0.0) {
                  shift.re = anorm / 0.40824829046386307;
                  shift.im = 0.0;
                } else if (anorm == 0.0) {
                  shift.re = 0.0;
                  shift.im = scale / 0.40824829046386307;
                } else {
                  shift.re = anorm / 0.40824829046386307;
                  shift.im = scale / 0.40824829046386307;
                }

                jp1 = ilast + 6 * ilast;
                anorm = ascale * b_A[jp1].re;
                scale = ascale * b_A[jp1].im;
                if (scale == 0.0) {
                  ad22_re = anorm / 0.40824829046386307;
                  ad22_im = 0.0;
                } else if (anorm == 0.0) {
                  ad22_re = 0.0;
                  ad22_im = scale / 0.40824829046386307;
                } else {
                  ad22_re = anorm / 0.40824829046386307;
                  ad22_im = scale / 0.40824829046386307;
                }

                t1_re = 0.5 * (shift.re + ad22_re);
                t1_im = 0.5 * (shift.im + ad22_im);
                t1_im_tmp = t1_re * t1_im;
                jp1 = ilastm1 + 6 * ilast;
                anorm = ascale * b_A[jp1].re;
                scale = ascale * b_A[jp1].im;
                if (scale == 0.0) {
                  ascale_re = anorm / 0.40824829046386307;
                  ascale_im = 0.0;
                } else if (anorm == 0.0) {
                  ascale_re = 0.0;
                  ascale_im = scale / 0.40824829046386307;
                } else {
                  ascale_re = anorm / 0.40824829046386307;
                  ascale_im = scale / 0.40824829046386307;
                }

                jp1 = ilast + 6 * ilastm1;
                anorm = ascale * b_A[jp1].re;
                scale = ascale * b_A[jp1].im;
                if (scale == 0.0) {
                  temp2 = anorm / 0.40824829046386307;
                  anorm = 0.0;
                } else if (anorm == 0.0) {
                  temp2 = 0.0;
                  anorm = scale / 0.40824829046386307;
                } else {
                  temp2 = anorm / 0.40824829046386307;
                  anorm = scale / 0.40824829046386307;
                }

                scale = shift.re * ad22_re - shift.im * ad22_im;
                temp1 = shift.re * ad22_im + shift.im * ad22_re;
                shift.re = ((t1_re * t1_re - t1_im * t1_im) + (ascale_re * temp2
                  - ascale_im * anorm)) - scale;
                shift.im = ((t1_im_tmp + t1_im_tmp) + (ascale_re * anorm +
                  ascale_im * temp2)) - temp1;
                b_sqrt(&shift);
                if ((t1_re - ad22_re) * shift.re + (t1_im - ad22_im) * shift.im <=
                    0.0) {
                  shift.re += t1_re;
                  shift.im += t1_im;
                } else {
                  shift.re = t1_re - shift.re;
                  shift.im = t1_im - shift.im;
                }
              } else {
                double ascale_re;
                double ascale_im;
                jp1 = ilast + 6 * ilastm1;
                anorm = ascale * b_A[jp1].re;
                scale = ascale * b_A[jp1].im;
                if (scale == 0.0) {
                  ascale_re = anorm / 0.40824829046386307;
                  ascale_im = 0.0;
                } else if (anorm == 0.0) {
                  ascale_re = 0.0;
                  ascale_im = scale / 0.40824829046386307;
                } else {
                  ascale_re = anorm / 0.40824829046386307;
                  ascale_im = scale / 0.40824829046386307;
                }

                eshift_re += ascale_re;
                eshift_im += ascale_im;
                shift.re = eshift_re;
                shift.im = eshift_im;
              }

              j = ilastm1;
              jp1 = ilastm1 + 1;
              exitg2 = false;
              while ((!exitg2) && (j + 1 > ifirst)) {
                istart = j + 1;
                ctemp_tmp = j + 6 * j;
                ctemp.re = ascale * b_A[ctemp_tmp].re - shift.re *
                  0.40824829046386307;
                ctemp.im = ascale * b_A[ctemp_tmp].im - shift.im *
                  0.40824829046386307;
                anorm = std::abs(ctemp.re) + std::abs(ctemp.im);
                jp1 += 6 * j;
                temp2 = ascale * (std::abs(b_A[jp1].re) + std::abs(b_A[jp1].im));
                scale = anorm;
                if (temp2 > anorm) {
                  scale = temp2;
                }

                if ((scale < 1.0) && (scale != 0.0)) {
                  anorm /= scale;
                  temp2 /= scale;
                }

                ctemp_tmp = j + 6 * (j - 1);
                if ((std::abs(b_A[ctemp_tmp].re) + std::abs(b_A[ctemp_tmp].im)) *
                    temp2 <= anorm * b_atol) {
                  goto90 = true;
                  exitg2 = true;
                } else {
                  jp1 = j;
                  j--;
                }
              }

              if (!goto90) {
                istart = ifirst;
                ctemp_tmp = (ifirst + 6 * (ifirst - 1)) - 1;
                ctemp.re = ascale * b_A[ctemp_tmp].re - shift.re *
                  0.40824829046386307;
                ctemp.im = ascale * b_A[ctemp_tmp].im - shift.im *
                  0.40824829046386307;
              }

              goto90 = false;
              jp1 = istart + 6 * (istart - 1);
              b_ascale.re = ascale * b_A[jp1].re;
              b_ascale.im = ascale * b_A[jp1].im;
              b_xzlartg(ctemp, b_ascale, &anorm, &shift);
              j = istart;
              jp1 = istart - 2;
              while (j < ilast + 1) {
                int ad22_re_tmp;
                if (j > istart) {
                  ctemp_tmp = j + 6 * jp1;
                  xzlartg(b_A[ctemp_tmp - 1], b_A[ctemp_tmp], &anorm, &shift,
                          &b_A[(j + 6 * jp1) - 1]);
                  b_A[ctemp_tmp].re = 0.0;
                  b_A[ctemp_tmp].im = 0.0;
                }

                for (jp1 = j; jp1 <= ilastm; jp1++) {
                  ctemp_tmp = j + 6 * (jp1 - 1);
                  ad22_re_tmp = ctemp_tmp - 1;
                  ad22_re = anorm * b_A[ad22_re_tmp].re + (shift.re *
                    b_A[ctemp_tmp].re - shift.im * b_A[ctemp_tmp].im);
                  ad22_im = anorm * b_A[ad22_re_tmp].im + (shift.re *
                    b_A[ctemp_tmp].im + shift.im * b_A[ctemp_tmp].re);
                  scale = b_A[ad22_re_tmp].re;
                  b_A[ctemp_tmp].re = anorm * b_A[ctemp_tmp].re - (shift.re *
                    b_A[ad22_re_tmp].re + shift.im * b_A[ad22_re_tmp].im);
                  b_A[ctemp_tmp].im = anorm * b_A[ctemp_tmp].im - (shift.re *
                    b_A[ad22_re_tmp].im - shift.im * scale);
                  b_A[ad22_re_tmp].re = ad22_re;
                  b_A[ad22_re_tmp].im = ad22_im;
                }

                shift.re = -shift.re;
                shift.im = -shift.im;
                jp1 = j;
                if (ilast + 1 < j + 2) {
                  jp1 = ilast - 1;
                }

                for (i = ifirst; i <= jp1 + 2; i++) {
                  ctemp_tmp = (i + 6 * (j - 1)) - 1;
                  ad22_re_tmp = (i + 6 * j) - 1;
                  ad22_re = anorm * b_A[ad22_re_tmp].re + (shift.re *
                    b_A[ctemp_tmp].re - shift.im * b_A[ctemp_tmp].im);
                  ad22_im = anorm * b_A[ad22_re_tmp].im + (shift.re *
                    b_A[ctemp_tmp].im + shift.im * b_A[ctemp_tmp].re);
                  scale = b_A[ad22_re_tmp].re;
                  b_A[ctemp_tmp].re = anorm * b_A[ctemp_tmp].re - (shift.re *
                    b_A[ad22_re_tmp].re + shift.im * b_A[ad22_re_tmp].im);
                  b_A[ctemp_tmp].im = anorm * b_A[ctemp_tmp].im - (shift.re *
                    b_A[ad22_re_tmp].im - shift.im * scale);
                  b_A[ad22_re_tmp].re = ad22_re;
                  b_A[ad22_re_tmp].im = ad22_im;
                }

                jp1 = j - 1;
                j++;
              }
            }

            jiter++;
          }
        }
      } else {
        guard2 = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  } else {
    guard1 = true;
  }

  if (guard2) {
    if (firstNonZero) {
      *info = ilast + 1;
      for (jp1 = 0; jp1 <= ilast; jp1++) {
        alpha1[jp1].re = rtNaN;
        alpha1[jp1].im = 0.0;
        beta1[jp1].re = rtNaN;
        beta1[jp1].im = 0.0;
      }
    } else {
      guard1 = true;
    }
  }

  if (guard1) {
    for (j = 0; j <= ilo - 2; j++) {
      alpha1[j] = b_A[j + 6 * j];
    }
  }
}

static void xzlartg(const creal_T f, const creal_T g, double *cs, creal_T *sn,
                    creal_T *r)
{
  double scale_tmp;
  double f2s;
  double scale;
  double fs_re;
  double fs_im;
  double gs_re;
  double gs_im;
  int count;
  int rescaledir;
  bool guard1 = false;
  double f2;
  scale_tmp = std::abs(f.re);
  f2s = std::abs(f.im);
  if (f2s > scale_tmp) {
    scale_tmp = f2s;
  }

  f2s = std::abs(g.re);
  scale = std::abs(g.im);
  if (scale > f2s) {
    f2s = scale;
  }

  scale = scale_tmp;
  if (f2s > scale_tmp) {
    scale = f2s;
  }

  fs_re = f.re;
  fs_im = f.im;
  gs_re = g.re;
  gs_im = g.im;
  count = -1;
  rescaledir = 0;
  guard1 = false;
  if (scale >= 7.4428285367870146E+137) {
    do {
      count++;
      fs_re *= 1.3435752215134178E-138;
      fs_im *= 1.3435752215134178E-138;
      gs_re *= 1.3435752215134178E-138;
      gs_im *= 1.3435752215134178E-138;
      scale *= 1.3435752215134178E-138;
    } while (!(scale < 7.4428285367870146E+137));

    rescaledir = 1;
    guard1 = true;
  } else if (scale <= 1.3435752215134178E-138) {
    if ((g.re == 0.0) && (g.im == 0.0)) {
      *cs = 1.0;
      sn->re = 0.0;
      sn->im = 0.0;
      *r = f;
    } else {
      do {
        count++;
        fs_re *= 7.4428285367870146E+137;
        fs_im *= 7.4428285367870146E+137;
        gs_re *= 7.4428285367870146E+137;
        gs_im *= 7.4428285367870146E+137;
        scale *= 7.4428285367870146E+137;
      } while (!(scale > 1.3435752215134178E-138));

      rescaledir = -1;
      guard1 = true;
    }
  } else {
    guard1 = true;
  }

  if (guard1) {
    double g2;
    f2 = fs_re * fs_re + fs_im * fs_im;
    g2 = gs_re * gs_re + gs_im * gs_im;
    scale = g2;
    if (1.0 > g2) {
      scale = 1.0;
    }

    if (f2 <= scale * 2.0041683600089728E-292) {
      if ((f.re == 0.0) && (f.im == 0.0)) {
        *cs = 0.0;
        r->re = rt_hypotd_snf(g.re, g.im);
        r->im = 0.0;
        f2 = rt_hypotd_snf(gs_re, gs_im);
        sn->re = gs_re / f2;
        sn->im = -gs_im / f2;
      } else {
        g2 = std::sqrt(g2);
        *cs = rt_hypotd_snf(fs_re, fs_im) / g2;
        if (scale_tmp > 1.0) {
          f2 = rt_hypotd_snf(f.re, f.im);
          fs_re = f.re / f2;
          fs_im = f.im / f2;
        } else {
          scale = 7.4428285367870146E+137 * f.re;
          f2s = 7.4428285367870146E+137 * f.im;
          f2 = rt_hypotd_snf(scale, f2s);
          fs_re = scale / f2;
          fs_im = f2s / f2;
        }

        gs_re /= g2;
        gs_im = -gs_im / g2;
        sn->re = fs_re * gs_re - fs_im * gs_im;
        sn->im = fs_re * gs_im + fs_im * gs_re;
        r->re = *cs * f.re + (sn->re * g.re - sn->im * g.im);
        r->im = *cs * f.im + (sn->re * g.im + sn->im * g.re);
      }
    } else {
      f2s = std::sqrt(g2 / f2 + 1.0);
      r->re = f2s * fs_re;
      r->im = f2s * fs_im;
      *cs = 1.0 / f2s;
      f2 += g2;
      f2s = r->re / f2;
      scale = r->im / f2;
      sn->re = f2s * gs_re - scale * -gs_im;
      sn->im = f2s * -gs_im + scale * gs_re;
      if (rescaledir > 0) {
        for (rescaledir = 0; rescaledir <= count; rescaledir++) {
          r->re *= 7.4428285367870146E+137;
          r->im *= 7.4428285367870146E+137;
        }
      } else {
        if (rescaledir < 0) {
          for (rescaledir = 0; rescaledir <= count; rescaledir++) {
            r->re *= 1.3435752215134178E-138;
            r->im *= 1.3435752215134178E-138;
          }
        }
      }
    }
  }
}

void c_butterworth(const coder::array<double, 1U> &x, double fc, double fs,
                   coder::array<double, 1U> &y)
{
  double wo;
  double a[36];
  int i;
  int b_i;
  double c[6];
  double t1_tmp[36];
  creal_T b_c[3];
  static const creal_T dcv[6] = { { -0.96592582628906831,// re
      -0.25881904510252057             // im
    }, { -0.96592582628906831,         // re
      0.25881904510252057              // im
    }, { -0.70710678118654746,         // re
      -0.70710678118654757             // im
    }, { -0.70710678118654746,         // re
      0.70710678118654757              // im
    }, { -0.25881904510252063,         // re
      -0.96592582628906831             // im
    }, { -0.25881904510252063,         // re
      0.96592582628906831              // im
    } };

  int k;
  double re;
  double wn;
  double d;
  int ipiv[6];
  int jBcol;
  double T[36];
  int j;
  double b1_idx_1;
  double t_idx_3;
  double b_a[4];
  bool p;
  int kAcol;
  int T_tmp;
  int a_tmp;
  creal_T alpha1[6];
  creal_T c_c[7];
  double a1[4];
  creal_T beta1[6];
  double c_a[7];
  double b[7];
  coder::array<double, 1U> r;
  static const signed char d_a[7] = { 1, 6, 15, 20, 15, 6, 1 };

  if (!isInitialized_c_butterworth) {
    c_butterworth_initialize();
  }

  wo = 4.0 * std::tan(3.1415926535897931 * (fc / (fs / 2.0)) / 2.0);
  std::memset(&a[0], 0, 36U * sizeof(double));
  for (i = 0; i < 6; i++) {
    c[i] = 0.0;
  }

  for (b_i = 0; b_i <= 5; b_i += 2) {
    b_c[1].re = -dcv[b_i].re - -dcv[b_i].im * 0.0;
    re = dcv[b_i + 1].re;
    wn = dcv[b_i + 1].im;
    b_c[2].re = -re * b_c[1].re - -wn * -dcv[b_i].im;
    d = b_c[1].re;
    for (k = 2; k >= 2; k--) {
      d -= re - wn * 0.0;
    }

    wn = std::sqrt(rt_hypotd_snf(dcv[b_i].re, dcv[b_i].im) * rt_hypotd_snf(re,
      wn));
    b1_idx_1 = 1.0 / wn;
    t_idx_3 = 1.0 / wn;
    wn = (1.0 - -d * 0.0) / t_idx_3;
    b_a[1] = wn;
    b_a[0] = -d - wn * 0.0;
    wn = (0.0 - -b_c[2].re * 0.0) / t_idx_3;
    b_a[3] = wn;
    b_a[2] = -b_c[2].re - wn * 0.0;
    for (i = 0; i < 2; i++) {
      d = b_a[i + 2];
      a1[i] = b_a[i] + d * 0.0;
      a1[i + 2] = b_a[i] * 0.0 + d * b1_idx_1;
    }

    if (b_i - 1 == -1) {
      a[0] = a1[0];
      a[1] = a1[1];
      c[0] = 0.0;
      a[6] = a1[2];
      a[7] = a1[3];
      c[1] = t_idx_3;
    } else {
      jBcol = b_i + 6 * b_i;
      kAcol = b_i + 6 * (b_i + 1);
      for (i = 0; i < b_i; i++) {
        a_tmp = b_i + 6 * i;
        a[a_tmp] = c[i];
        a[a_tmp + 1] = 0.0 * c[i];
        c[i] *= 0.0;
      }

      a[jBcol] = a1[0];
      a[jBcol + 1] = a1[1];
      a[kAcol] = a1[2];
      a[kAcol + 1] = a1[3];
      c[b_i] = 0.0;
      c[b_i + 1] = t_idx_3;
    }
  }

  std::memset(&t1_tmp[0], 0, 36U * sizeof(double));
  for (k = 0; k < 6; k++) {
    t1_tmp[k + 6 * k] = 1.0;
  }

  for (i = 0; i < 36; i++) {
    d = wo * a[i] * 0.5 / 2.0;
    a[i] = d;
    T[i] = t1_tmp[i] + d;
    t1_tmp[i] -= d;
  }

  xgetrf(t1_tmp, ipiv, &jBcol);
  for (b_i = 0; b_i < 5; b_i++) {
    if (ipiv[b_i] != b_i + 1) {
      for (j = 0; j < 6; j++) {
        jBcol = b_i + 6 * j;
        wn = T[jBcol];
        T_tmp = (ipiv[b_i] + 6 * j) - 1;
        T[jBcol] = T[T_tmp];
        T[T_tmp] = wn;
      }
    }
  }

  for (j = 0; j < 6; j++) {
    jBcol = 6 * j;
    for (k = 0; k < 6; k++) {
      kAcol = 6 * k;
      i = k + jBcol;
      if (T[i] != 0.0) {
        a_tmp = k + 2;
        for (b_i = a_tmp; b_i < 7; b_i++) {
          T_tmp = (b_i + jBcol) - 1;
          T[T_tmp] -= T[i] * t1_tmp[(b_i + kAcol) - 1];
        }
      }
    }
  }

  for (j = 0; j < 6; j++) {
    jBcol = 6 * j;
    for (k = 5; k >= 0; k--) {
      kAcol = 6 * k;
      i = k + jBcol;
      if (T[i] != 0.0) {
        T[i] /= t1_tmp[k + kAcol];
        for (b_i = 0; b_i < k; b_i++) {
          T_tmp = b_i + jBcol;
          T[T_tmp] -= T[i] * t1_tmp[b_i + kAcol];
        }
      }
    }
  }

  p = true;
  for (k = 0; k < 36; k++) {
    if ((!p) || (std::isinf(T[k]) || std::isnan(T[k]))) {
      p = false;
    }
  }

  if (!p) {
    for (b_i = 0; b_i < 6; b_i++) {
      alpha1[b_i].re = rtNaN;
      alpha1[b_i].im = 0.0;
    }
  } else {
    bool exitg2;
    p = true;
    j = 0;
    exitg2 = false;
    while ((!exitg2) && (j < 6)) {
      int exitg1;
      b_i = 0;
      do {
        exitg1 = 0;
        if (b_i <= j) {
          if (!(T[b_i + 6 * j] == T[j + 6 * b_i])) {
            p = false;
            exitg1 = 1;
          } else {
            b_i++;
          }
        } else {
          j++;
          exitg1 = 2;
        }
      } while (exitg1 == 0);

      if (exitg1 == 1) {
        exitg2 = true;
      }
    }

    if (p) {
      p = true;
      for (k = 0; k < 36; k++) {
        if ((!p) || (std::isinf(T[k]) || std::isnan(T[k]))) {
          p = false;
        }
      }

      if (!p) {
        for (i = 0; i < 36; i++) {
          T[i] = rtNaN;
        }

        jBcol = 2;
        for (j = 0; j < 5; j++) {
          if (jBcol <= 6) {
            std::memset(&T[(j * 6 + jBcol) + -1], 0, (7 - jBcol) * sizeof(double));
          }

          jBcol++;
        }
      } else {
        xgehrd(T);
        eml_dlahqr(T);
        jBcol = 4;
        for (j = 0; j < 3; j++) {
          if (jBcol <= 6) {
            std::memset(&T[(j * 6 + jBcol) + -1], 0, (7 - jBcol) * sizeof(double));
          }

          jBcol++;
        }
      }

      for (k = 0; k < 6; k++) {
        alpha1[k].re = T[k + 6 * k];
        alpha1[k].im = 0.0;
      }
    } else {
      xzgeev(T, &jBcol, alpha1, beta1);
      for (i = 0; i < 6; i++) {
        if (beta1[i].im == 0.0) {
          if (alpha1[i].im == 0.0) {
            re = alpha1[i].re / beta1[i].re;
            wn = 0.0;
          } else if (alpha1[i].re == 0.0) {
            re = 0.0;
            wn = alpha1[i].im / beta1[i].re;
          } else {
            re = alpha1[i].re / beta1[i].re;
            wn = alpha1[i].im / beta1[i].re;
          }
        } else if (beta1[i].re == 0.0) {
          if (alpha1[i].re == 0.0) {
            re = alpha1[i].im / beta1[i].im;
            wn = 0.0;
          } else if (alpha1[i].im == 0.0) {
            re = 0.0;
            wn = -(alpha1[i].re / beta1[i].im);
          } else {
            re = alpha1[i].im / beta1[i].im;
            wn = -(alpha1[i].re / beta1[i].im);
          }
        } else {
          b1_idx_1 = std::abs(beta1[i].re);
          wn = std::abs(beta1[i].im);
          if (b1_idx_1 > wn) {
            wn = beta1[i].im / beta1[i].re;
            t_idx_3 = beta1[i].re + wn * beta1[i].im;
            re = (alpha1[i].re + wn * alpha1[i].im) / t_idx_3;
            wn = (alpha1[i].im - wn * alpha1[i].re) / t_idx_3;
          } else if (wn == b1_idx_1) {
            if (beta1[i].re > 0.0) {
              t_idx_3 = 0.5;
            } else {
              t_idx_3 = -0.5;
            }

            if (beta1[i].im > 0.0) {
              wn = 0.5;
            } else {
              wn = -0.5;
            }

            re = (alpha1[i].re * t_idx_3 + alpha1[i].im * wn) / b1_idx_1;
            wn = (alpha1[i].im * t_idx_3 - alpha1[i].re * wn) / b1_idx_1;
          } else {
            wn = beta1[i].re / beta1[i].im;
            t_idx_3 = beta1[i].im + wn * beta1[i].re;
            re = (wn * alpha1[i].re + alpha1[i].im) / t_idx_3;
            wn = (wn * alpha1[i].im - alpha1[i].re) / t_idx_3;
          }
        }

        alpha1[i].re = re;
        alpha1[i].im = wn;
      }
    }
  }

  c_c[0].re = 1.0;
  c_c[0].im = 0.0;
  for (j = 0; j < 6; j++) {
    d = c_c[j].re;
    c_c[j + 1].re = -alpha1[j].re * c_c[j].re - -alpha1[j].im * c_c[j].im;
    c_c[j + 1].im = -alpha1[j].re * c_c[j].im + -alpha1[j].im * d;
    for (k = j + 1; k >= 2; k--) {
      wn = alpha1[j].re * c_c[k - 2].im + alpha1[j].im * c_c[k - 2].re;
      c_c[k - 1].re -= alpha1[j].re * c_c[k - 2].re - alpha1[j].im * c_c[k - 2].
        im;
      c_c[k - 1].im -= wn;
    }
  }

  t_idx_3 = 0.0;
  b1_idx_1 = 0.0;
  for (i = 0; i < 7; i++) {
    c_a[i] = c_c[i].re;
    t_idx_3 += c_c[i].re;
    b1_idx_1 += 0.0 * c_c[i].re;
  }

  for (i = 0; i < 7; i++) {
    wn = static_cast<double>(d_a[i]) * t_idx_3;
    if (static_cast<double>(d_a[i]) * b1_idx_1 == 0.0) {
      wn /= 64.0;
    } else if (wn == 0.0) {
      wn = 0.0;
    } else {
      wn /= 64.0;
    }

    b[i] = wn;
  }

  filter(b, c_a, x, r);
  flip(r);
  filter(b, c_a, r, y);
  flip(y);
}

void c_butterworth_initialize()
{
  rt_InitInfAndNaN();
  isInitialized_c_butterworth = true;
}

void c_butterworth_terminate()
{
  // (no terminate code required)
  isInitialized_c_butterworth = false;
}

// End of code generation (c_butterworth.cpp)
