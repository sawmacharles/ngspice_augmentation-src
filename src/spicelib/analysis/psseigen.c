/**********
card #458 second pass (plan T4.1-T4.3), SUSPECTED CAUSE 1: replace the
power-iteration spectral-radius estimate (which cannot resolve a
complex-conjugate dominant multiplier pair -- exactly what a slow
oscillation produces) with a full eigenvalue decomposition of the small
(n<=~10) monodromy matrix: real Hessenberg reduction by Householder
reflections, then the classical real double-shift implicit QR algorithm
("Francis QR", the method underlying EISPACK's hqr / LAPACK's dhseqr;
implemented here from the standard mathematical description -- see e.g.
Golub & Van Loan, "Matrix Computations", ch.7 -- as stdlib-only C, no
external library dependency, since n is always small here: 1 on R2, 7 on
R3/R4).

VALIDATED (this card) against an INDEPENDENT root-finder (Faddeev-LeVerrier
characteristic polynomial + Durand-Kerner simultaneous root iteration,
implemented separately in Python, tools/pss_ref.py:eig_reference()) on six
cases before being wired into pssshoot.c: a 2x2 pure-rotation matrix
(known exact eigenvalues +-i), a 3x3 diagonal matrix, a 3x3 companion
matrix with known real roots 1/2/3, a second 2x2 rotation (known exact
eigenvalues 0.6+-0.8i), a 4x4 asymmetric matrix with two complex-conjugate
pairs, and a 7x7 matrix (matching the R3/R4 state count) with a mix of
three complex pairs and one real eigenvalue. All six matched the
independent method to the digits printed (double precision). See the card
#458 report for the transcript.
**********/

#include "ngspice/ngspice.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "pssshoot.h"

#define PSS_SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

/* Householder Hessenberg reduction, 1-based (n+1)-sized arrays (row/col 0
 * unused) -- standard algorithm (Golub & Van Loan Alg. 7.4.2 structure).
 * Eigenvectors are not needed, so no transformation is accumulated. */
static void
pss_hessenberg_reduce(double **a, int n)
{
    int i, j, m;
    double scale, h, f, g;
    double *tmp = TMALLOC(double, n + 1);

    for (m = 2; m <= n - 1; m++) {
        scale = 0.0;
        for (i = m; i <= n; i++)
            scale += fabs(a[i][m-1]);
        if (scale != 0.0) {
            h = 0.0;
            for (i = n; i >= m; i--) {
                tmp[i] = a[i][m-1] / scale;
                h += tmp[i] * tmp[i];
            }
            g = (tmp[m] >= 0.0) ? -sqrt(h) : sqrt(h);
            h -= tmp[m] * g;
            tmp[m] -= g;
            for (j = m; j <= n; j++) {
                f = 0.0;
                for (i = n; i >= m; i--)
                    f += tmp[i] * a[i][j];
                f /= h;
                for (i = m; i <= n; i++)
                    a[i][j] -= f * tmp[i];
            }
            for (i = 1; i <= n; i++) {
                f = 0.0;
                for (j = n; j >= m; j--)
                    f += tmp[j] * a[i][j];
                f /= h;
                for (j = m; j <= n; j++)
                    a[i][j] -= f * tmp[j];
            }
            a[m][m-1] = scale * g;
            for (i = m + 1; i <= n; i++)
                a[i][m-1] = 0.0;
        }
    }
    tfree(tmp);
}

/* Real double-shift implicit QR with deflation ("hqr"), 1-based. Destroys
 * a[][]. Returns 0 on success, -1 if a sub-block failed to converge within
 * the iteration budget (should not happen for the well-scaled, physically
 * bounded monodromy matrices this analysis builds; reported honestly to
 * the caller rather than silently returning a wrong answer). */
static int
pss_hqr(double **a, int n, double *wr, double *wi)
{
    int nn, m, l, k, j, its, i, mmin;
    double z = 0, y = 0, x = 0, w = 0, s = 0, r = 0, q = 0, p = 0, anorm = 0.0, t;

    for (i = 1; i <= n; i++)
        for (j = (i - 1 > 1 ? i - 1 : 1); j <= n; j++)
            anorm += fabs(a[i][j]);

    nn = n;
    t = 0.0;
    while (nn >= 1) {
        its = 0;
        do {
            for (l = nn; l >= 2; l--) {
                s = fabs(a[l-1][l-1]) + fabs(a[l][l]);
                if (s == 0.0) s = anorm;
                if (fabs(a[l][l-1]) <= 1e-14 * s) break;
            }
            x = a[nn][nn];
            if (l == nn) {
                wr[nn] = x + t;
                wi[nn] = 0.0;
                nn--;
            } else {
                y = a[nn-1][nn-1];
                w = a[nn][nn-1] * a[nn-1][nn];
                if (l == nn - 1) {
                    p = 0.5 * (y - x);
                    q = p * p + w;
                    z = sqrt(fabs(q));
                    x += t;
                    if (q >= 0.0) {
                        z = p + PSS_SIGN(z, p);
                        wr[nn-1] = wr[nn] = x + z;
                        if (z != 0.0) wr[nn] = x - w / z;
                        wi[nn-1] = wi[nn] = 0.0;
                    } else {
                        wr[nn-1] = wr[nn] = x + p;
                        wi[nn-1] = -(wi[nn] = z);
                    }
                    nn -= 2;
                } else {
                    if (its == 60)
                        return -1;
                    if (its == 10 || its == 20 || its == 30 || its == 40 || its == 50) {
                        t += x;
                        for (i = 1; i <= nn; i++) a[i][i] -= x;
                        s = fabs(a[nn][nn-1]) + fabs(a[nn-1][nn-2]);
                        y = x = 0.75 * s;
                        w = -0.4375 * s * s;
                    }
                    its++;
                    for (m = nn - 2; m >= l; m--) {
                        z = a[m][m];
                        r = x - z;
                        s = y - z;
                        p = (r * s - w) / a[m+1][m] + a[m][m+1];
                        q = a[m+1][m+1] - z - r - s;
                        r = a[m+2][m+1];
                        s = fabs(p) + fabs(q) + fabs(r);
                        p /= s; q /= s; r /= s;
                        if (m == l) break;
                        {
                            double u = fabs(a[m][m-1]) * (fabs(q) + fabs(r));
                            double v = fabs(p) * (fabs(a[m-1][m-1]) + fabs(z) + fabs(a[m+1][m+1]));
                            if (u <= 1e-14 * v) break;
                        }
                    }
                    for (i = m + 2; i <= nn; i++) {
                        a[i][i-2] = 0.0;
                        if (i != m + 2) a[i][i-3] = 0.0;
                    }
                    for (k = m; k <= nn - 1; k++) {
                        if (k != m) {
                            p = a[k][k-1];
                            q = a[k+1][k-1];
                            r = (k != nn - 1) ? a[k+2][k-1] : 0.0;
                            x = fabs(p) + fabs(q) + fabs(r);
                            if (x != 0.0) { p /= x; q /= x; r /= x; }
                        }
                        s = PSS_SIGN(sqrt(p*p + q*q + r*r), p);
                        if (s != 0.0) {
                            if (k == m) {
                                if (l != m) a[k][k-1] = -a[k][k-1];
                            } else {
                                a[k][k-1] = -s * x;
                            }
                            p += s;
                            x = p / s; y = q / s; z = r / s;
                            q /= p; r /= p;
                            for (j = k; j <= nn; j++) {
                                double pp = a[k][j] + q * a[k+1][j];
                                if (k != nn - 1) {
                                    pp += r * a[k+2][j];
                                    a[k+2][j] -= pp * z;
                                }
                                a[k+1][j] -= pp * y;
                                a[k][j] -= pp * x;
                            }
                            mmin = (nn < k + 3) ? nn : k + 3;
                            for (i = l; i <= mmin; i++) {
                                double pp = x * a[i][k] + y * a[i][k+1];
                                if (k != nn - 1) {
                                    pp += z * a[i][k+2];
                                    a[i][k+2] -= pp * r;
                                }
                                a[i][k+1] -= pp * q;
                                a[i][k] -= pp;
                            }
                        }
                    }
                }
            }
        } while (l < nn - 1);
    }
    return 0;
}

/* Public entry: eigenvalues of the n*n row-major real matrix Ain. eig[]
 * (caller-allocated, n entries) receives the eigenvalues, complex
 * conjugate pairs adjacent. Returns 0 on success, -1 if a sub-block did
 * not converge (caller must treat this as "no trustworthy Floquet report",
 * not silently use a partial/garbage result). */
int
pss_eigenvalues(const double *Ain, int n, PSScomplex *eig)
{
    double **a;
    double *wr, *wi;
    int i, j, err;

    if (n <= 0)
        return -1;

    a  = TMALLOC(double *, n + 1);
    wr = TMALLOC(double, n + 1);
    wi = TMALLOC(double, n + 1);
    if (!a || !wr || !wi) {
        if (a) tfree(a);
        if (wr) tfree(wr);
        if (wi) tfree(wi);
        return -1;
    }
    for (i = 1; i <= n; i++) {
        a[i] = TMALLOC(double, n + 1);
        for (j = 1; j <= n; j++)
            a[i][j] = Ain[(i - 1) * n + (j - 1)];
    }

    if (n == 1) {
        eig[0].re = a[1][1];
        eig[0].im = 0.0;
        err = 0;
    } else {
        pss_hessenberg_reduce(a, n);
        err = pss_hqr(a, n, wr, wi);
        for (i = 0; i < n; i++) {
            eig[i].re = wr[i + 1];
            eig[i].im = wi[i + 1];
        }
    }

    for (i = 1; i <= n; i++) tfree(a[i]);
    tfree(a);
    tfree(wr);
    tfree(wi);
    return err;
}

/* Classify the Floquet multiplier set: spectral radius (max |lambda_i|)
 * plus which multiplier is dominant and why. See pssshoot.h for the
 * PSSfloquet fields. SUSPECTED CAUSE 1's specific failure mode (a slow
 * oscillation from an unresolved complex-conjugate dominant pair) is
 * exactly what pss_eigenvalues() resolves that power iteration could not:
 * report the dominant pair's angle, not just its magnitude. */
void
pss_classify_floquet(const PSScomplex *eig, int n, PSSfloquet *out)
{
    int i, dom = 0;
    double domabs;

    memset(out, 0, sizeof *out);
    if (n <= 0) {
        out->valid = 0;
        return;
    }
    out->valid = 1;
    domabs = hypot(eig[0].re, eig[0].im);
    for (i = 1; i < n; i++) {
        double a = hypot(eig[i].re, eig[i].im);
        if (a > domabs) { domabs = a; dom = i; }
    }
    out->dominant = eig[dom];
    out->spectral_radius = domabs;
    out->unstable = (domabs >= 1.0);

    /* Period-doubling: a real eigenvalue near -1 (angle ~180deg), whether
     * or not it is the single largest-magnitude one -- report it if
     * present among multipliers close to the unit circle, since a real
     * lambda<=-1 elsewhere in the spectrum but still >=1 in magnitude is
     * exactly the period-2 subharmonic signature (peak-current-mode
     * without slope compensation, refs/README.md's R4 D~0.55 case). */
    out->period_doubling = 0;
    out->oscillation = 0;
    for (i = 0; i < n; i++) {
        double a = hypot(eig[i].re, eig[i].im);
        if (a < 0.999 * domabs - 1e-12)
            continue;   /* only consider multipliers at/near the dominant magnitude */
        if (fabs(eig[i].im) < 1e-6 * (a > 1e-12 ? a : 1.0) && eig[i].re < 0.0) {
            out->period_doubling = 1;
        } else if (fabs(eig[i].im) > 1e-6 * (a > 1e-12 ? a : 1.0)) {
            out->oscillation = 1;
        }
    }
}
