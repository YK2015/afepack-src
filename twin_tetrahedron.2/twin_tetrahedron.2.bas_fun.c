/*************************************************************************
 *  twin_tetrahedron.2.bas_fun.c
 *
 *  Tutorial implementation of a continuous piecewise P2 basis on the
 *  twin-tetrahedron macro element.
 */

double get_volume(const double * v0,
                  const double * v1,
                  const double * v2,
                  const double * v3)
{
  return ((v1[0] - v0[0])*(v2[1] - v0[1])*(v3[2] - v0[2])
          + (v1[1] - v0[1])*(v2[2] - v0[2])*(v3[0] - v0[0])
          + (v1[2] - v0[2])*(v2[0] - v0[0])*(v3[1] - v0[1])
          - (v1[0] - v0[0])*(v2[2] - v0[2])*(v3[1] - v0[1])
          - (v1[1] - v0[1])*(v2[0] - v0[0])*(v3[2] - v0[2])
          - (v1[2] - v0[2])*(v2[1] - v0[1])*(v3[0] - v0[0]));
}

#define ZERO 1.0e-12
#define co_det(v, m, n) (\
    ((m%2==0)?-1.:1.) * (\
      (v[(m+2)%4][(n+1)%3] - v[(m+1)%4][(n+1)%3]) \
    * (v[(m+3)%4][(n+2)%3] - v[(m+1)%4][(n+2)%3]) \
    - (v[(m+2)%4][(n+2)%3] - v[(m+1)%4][(n+2)%3]) \
    * (v[(m+3)%4][(n+1)%3] - v[(m+1)%4][(n+1)%3]) \
  ) \
  )

static void zero_vec(double * val)
{
  val[0] = 0.0;
  val[1] = 0.0;
  val[2] = 0.0;
}

static void linear_value(int idx, const double * p, const double ** v, double * out)
{
  double volume, sign;

  if (idx == 0) {
    volume = get_volume(v[0], v[1], v[3], v[4]);
    out[0] = get_volume(p, v[1], v[3], v[4])/volume;
    return;
  }

  if (idx == 1) {
    volume = get_volume(v[0], v[1], v[2], v[4]);
    out[0] = get_volume(v[0], p, v[2], v[4])/volume;
    if (out[0] < 0.0) out[0] = 0.0;
    return;
  }

  if (idx == 2) {
    sign = get_volume(v[0], p, v[2], v[4]);
    volume = 0.5 * get_volume(v[0], v[1], v[3], v[4]);
    if (sign > 0.0) {
      out[0] = get_volume(v[0], v[1], p, v[4])/volume;
    } else {
      out[0] = get_volume(v[0], p, v[3], v[4])/volume;
    }
    return;
  }

  if (idx == 3) {
    volume = get_volume(v[0], v[2], v[3], v[4]);
    out[0] = get_volume(v[0], v[2], p, v[4])/volume;
    if (out[0] < 0.0) out[0] = 0.0;
    return;
  }

  volume = get_volume(v[0], v[1], v[3], v[4]);
  out[0] = get_volume(v[0], v[1], v[3], p)/volume;
}

static void linear_gradient(int idx, const double * p, const double ** v, double * out)
{
  const double * v0[4] = {v[0], v[1], v[3], v[4]};
  const double * v1[4] = {v[0], v[1], v[2], v[4]};
  const double * v2[4] = {v[0], v[2], v[3], v[4]};
  double sign, volume;

  if (idx == 0) {
    volume = get_volume(v[0], v[1], v[3], v[4]);
    out[0] = co_det(v0, 0, 0)/volume;
    out[1] = co_det(v0, 0, 1)/volume;
    out[2] = co_det(v0, 0, 2)/volume;
    return;
  }

  if (idx == 1) {
    sign = get_volume(v[0], p, v[2], v[4]);
    volume = get_volume(v[0], v[1], v[2], v[4]);
    if (sign > ZERO) {
      out[0] = co_det(v1, 1, 0)/volume;
      out[1] = co_det(v1, 1, 1)/volume;
      out[2] = co_det(v1, 1, 2)/volume;
    } else if (sign < -ZERO) {
      zero_vec(out);
    } else {
      out[0] = 0.5*co_det(v1, 1, 0)/volume;
      out[1] = 0.5*co_det(v1, 1, 1)/volume;
      out[2] = 0.5*co_det(v1, 1, 2)/volume;
    }
    return;
  }

  if (idx == 2) {
    sign = get_volume(v[0], p, v[2], v[4]);
    volume = 0.5*get_volume(v[0], v[1], v[3], v[4]);
    if (sign > ZERO) {
      out[0] = co_det(v1, 2, 0)/volume;
      out[1] = co_det(v1, 2, 1)/volume;
      out[2] = co_det(v1, 2, 2)/volume;
    } else if (sign < -ZERO) {
      out[0] = co_det(v2, 1, 0)/volume;
      out[1] = co_det(v2, 1, 1)/volume;
      out[2] = co_det(v2, 1, 2)/volume;
    } else {
      out[0] = 0.5*(co_det(v1, 2, 0) + co_det(v2, 1, 0))/volume;
      out[1] = 0.5*(co_det(v1, 2, 1) + co_det(v2, 1, 1))/volume;
      out[2] = 0.5*(co_det(v1, 2, 2) + co_det(v2, 1, 2))/volume;
    }
    return;
  }

  if (idx == 3) {
    sign = get_volume(v[0], v[2], p, v[4]);
    volume = get_volume(v[0], v[1], v[2], v[4]);
    if (sign > ZERO) {
      out[0] = co_det(v2, 2, 0)/volume;
      out[1] = co_det(v2, 2, 1)/volume;
      out[2] = co_det(v2, 2, 2)/volume;
    } else if (sign < -ZERO) {
      zero_vec(out);
    } else {
      out[0] = 0.5*co_det(v2, 2, 0)/volume;
      out[1] = 0.5*co_det(v2, 2, 1)/volume;
      out[2] = 0.5*co_det(v2, 2, 2)/volume;
    }
    return;
  }

  volume = get_volume(v[0], v[1], v[3], v[4]);
  out[0] = co_det(v0, 3, 0)/volume;
  out[1] = co_det(v0, 3, 1)/volume;
  out[2] = co_det(v0, 3, 2)/volume;
}

static void vertex_quad_value(int idx, const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lam;
  linear_value(idx, p, v, &lam);
  val[0] = lam*(2.0*lam - 1.0);
}

static void vertex_quad_gradient(int idx, const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lam, grad[3];
  linear_value(idx, p, v, &lam);
  linear_gradient(idx, p, v, grad);
  val[0] = (4.0*lam - 1.0)*grad[0];
  val[1] = (4.0*lam - 1.0)*grad[1];
  val[2] = (4.0*lam - 1.0)*grad[2];
}

static void edge_quad_value(int ia, int ib, const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double la, lb;
  linear_value(ia, p, v, &la);
  linear_value(ib, p, v, &lb);
  val[0] = 4.0*la*lb;
}

static void edge_quad_gradient(int ia, int ib, const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double la, lb, ga[3], gb[3];
  linear_value(ia, p, v, &la);
  linear_value(ib, p, v, &lb);
  linear_gradient(ia, p, v, ga);
  linear_gradient(ib, p, v, gb);
  val[0] = 4.0*(la*gb[0] + lb*ga[0]);
  val[1] = 4.0*(la*gb[1] + lb*ga[1]);
  val[2] = 4.0*(la*gb[2] + lb*ga[2]);
}

void phi_1(const double * p, const double ** v, void * value) { vertex_quad_value(0, p, v, value); }
void phi_2(const double * p, const double ** v, void * value) { vertex_quad_value(1, p, v, value); }
void phi_3(const double * p, const double ** v, void * value) { vertex_quad_value(2, p, v, value); }
void phi_4(const double * p, const double ** v, void * value) { vertex_quad_value(3, p, v, value); }
void phi_5(const double * p, const double ** v, void * value) { vertex_quad_value(4, p, v, value); }

void gradient_phi_1(const double * p, const double ** v, void * value) { vertex_quad_gradient(0, p, v, value); }
void gradient_phi_2(const double * p, const double ** v, void * value) { vertex_quad_gradient(1, p, v, value); }
void gradient_phi_3(const double * p, const double ** v, void * value) { vertex_quad_gradient(2, p, v, value); }
void gradient_phi_4(const double * p, const double ** v, void * value) { vertex_quad_gradient(3, p, v, value); }
void gradient_phi_5(const double * p, const double ** v, void * value) { vertex_quad_gradient(4, p, v, value); }

void psi_1(const double * p, const double ** v, void * value) { edge_quad_value(0, 1, p, v, value); }
void psi_2(const double * p, const double ** v, void * value) { edge_quad_value(0, 3, p, v, value); }
void psi_3(const double * p, const double ** v, void * value) { edge_quad_value(0, 4, p, v, value); }
void psi_4(const double * p, const double ** v, void * value) { edge_quad_value(3, 4, p, v, value); }
void psi_5(const double * p, const double ** v, void * value) { edge_quad_value(1, 4, p, v, value); }
void psi_6(const double * p, const double ** v, void * value) { edge_quad_value(1, 2, p, v, value); }
void psi_7(const double * p, const double ** v, void * value) { edge_quad_value(2, 3, p, v, value); }

void gradient_psi_1(const double * p, const double ** v, void * value) { edge_quad_gradient(0, 1, p, v, value); }
void gradient_psi_2(const double * p, const double ** v, void * value) { edge_quad_gradient(0, 3, p, v, value); }
void gradient_psi_3(const double * p, const double ** v, void * value) { edge_quad_gradient(0, 4, p, v, value); }
void gradient_psi_4(const double * p, const double ** v, void * value) { edge_quad_gradient(3, 4, p, v, value); }
void gradient_psi_5(const double * p, const double ** v, void * value) { edge_quad_gradient(1, 4, p, v, value); }
void gradient_psi_6(const double * p, const double ** v, void * value) { edge_quad_gradient(1, 2, p, v, value); }
void gradient_psi_7(const double * p, const double ** v, void * value) { edge_quad_gradient(2, 3, p, v, value); }

void chi_1(const double * p, const double ** v, void * value) { edge_quad_value(0, 2, p, v, value); }
void chi_2(const double * p, const double ** v, void * value) { edge_quad_value(2, 4, p, v, value); }

void gradient_chi_1(const double * p, const double ** v, void * value) { edge_quad_gradient(0, 2, p, v, value); }
void gradient_chi_2(const double * p, const double ** v, void * value) { edge_quad_gradient(2, 4, p, v, value); }

/*
 *  end of file
 *************************************************************************/
