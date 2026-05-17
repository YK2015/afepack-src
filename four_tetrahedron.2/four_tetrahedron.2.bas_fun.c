/*************************************************************************
 *  four_tetrahedron.2.bas_fun.c
 *
 *  Tutorial implementation of a continuous piecewise P2 basis on the
 *  four-tetrahedron macro element. The design reuses the piecewise linear
 *  nodal basis from four_tetrahedron.1 and lifts it through standard
 *  Lagrange product formulas.
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
  double volume, d45, d56, d64;

  if (idx == 0) {
    volume = get_volume(v[0], v[1], v[2], v[3]);
    out[0] = get_volume(p, v[1], v[2], v[3])/volume;
    return;
  }

  if (idx == 1) {
    volume = get_volume(v[0], v[1], v[6], v[5]);
    out[0] = get_volume(v[0], p, v[6], v[5])/volume;
    if (out[0] < 0.0) out[0] = 0.0;
    return;
  }

  if (idx == 2) {
    volume = get_volume(v[0], v[2], v[4], v[6]);
    out[0] = get_volume(v[0], p, v[4], v[6])/volume;
    if (out[0] < 0.0) out[0] = 0.0;
    return;
  }

  if (idx == 3) {
    volume = get_volume(v[0], v[3], v[5], v[4]);
    out[0] = get_volume(v[0], p, v[5], v[4])/volume;
    if (out[0] < 0.0) out[0] = 0.0;
    return;
  }

  volume = 0.25 * get_volume(v[0], v[1], v[2], v[3]);
  d45 = get_volume(v[0], p, v[4], v[5]);
  d56 = get_volume(v[0], p, v[5], v[6]);
  d64 = get_volume(v[0], p, v[6], v[4]);

  if (idx == 4) {
    if (d56 <= 0.0) {
      out[0] = 0.0;
      return;
    }
    if (d45 <= 0.0) {
      out[0] = get_volume(v[0], p, v[3], v[5])/volume;
      return;
    }
    if (d64 <= 0.0) {
      out[0] = get_volume(v[0], p, v[6], v[2])/volume;
      return;
    }
    out[0] = d56/volume;
    return;
  }

  if (idx == 5) {
    if (d64 <= 0.0) {
      out[0] = 0.0;
      return;
    }
    if (d45 <= 0.0) {
      out[0] = get_volume(v[0], p, v[4], v[3])/volume;
      return;
    }
    if (d56 <= 0.0) {
      out[0] = get_volume(v[0], p, v[1], v[6])/volume;
      return;
    }
    out[0] = d64/volume;
    return;
  }

  if (d45 <= 0.0) {
    out[0] = 0.0;
    return;
  }
  if (d56 <= 0.0) {
    out[0] = get_volume(v[0], p, v[5], v[1])/volume;
    return;
  }
  if (d64 <= 0.0) {
    out[0] = get_volume(v[0], p, v[2], v[4])/volume;
    return;
  }
  out[0] = d45/volume;
}

static void linear_gradient(int idx, const double * p, const double ** v, double * out)
{
  double volume, sign, d45, d56, d64, count;
  const double * v1[4] = {v[0], v[1], v[6], v[5]};
  const double * v2[4] = {v[0], v[2], v[4], v[6]};
  const double * v3[4] = {v[0], v[3], v[5], v[4]};
  const double * v4[4] = {v[0], v[4], v[5], v[6]};

  if (idx == 0) {
    volume = get_volume(v[0], v[1], v[2], v[3]);
    out[0] = co_det(v, 0, 0)/volume;
    out[1] = co_det(v, 0, 1)/volume;
    out[2] = co_det(v, 0, 2)/volume;
    return;
  }

  if (idx == 1) {
    volume = get_volume(v[0], v[1], v[6], v[5]);
    sign = get_volume(v[0], p, v[6], v[5])/volume;
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
    volume = get_volume(v[0], v[2], v[4], v[6]);
    sign = get_volume(v[0], p, v[4], v[6])/volume;
    if (sign > ZERO) {
      out[0] = co_det(v2, 1, 0)/volume;
      out[1] = co_det(v2, 1, 1)/volume;
      out[2] = co_det(v2, 1, 2)/volume;
    } else if (sign < -ZERO) {
      zero_vec(out);
    } else {
      out[0] = 0.5*co_det(v2, 1, 0)/volume;
      out[1] = 0.5*co_det(v2, 1, 1)/volume;
      out[2] = 0.5*co_det(v2, 1, 2)/volume;
    }
    return;
  }

  if (idx == 3) {
    volume = get_volume(v[0], v[3], v[5], v[4]);
    sign = get_volume(v[0], p, v[5], v[4])/volume;
    if (sign > ZERO) {
      out[0] = co_det(v3, 1, 0)/volume;
      out[1] = co_det(v3, 1, 1)/volume;
      out[2] = co_det(v3, 1, 2)/volume;
    } else if (sign < -ZERO) {
      zero_vec(out);
    } else {
      out[0] = 0.5*co_det(v3, 1, 0)/volume;
      out[1] = 0.5*co_det(v3, 1, 1)/volume;
      out[2] = 0.5*co_det(v3, 1, 2)/volume;
    }
    return;
  }

  volume = get_volume(v4[0], v4[1], v4[2], v4[3]);
  d45 = get_volume(v[0], p, v[4], v[5])/volume;
  d56 = get_volume(v[0], p, v[5], v[6])/volume;
  d64 = get_volume(v[0], p, v[6], v[4])/volume;
  count = 0.0;
  zero_vec(out);

  if (idx == 4) {
    if (d45 <= ZERO) {
      out[0] += co_det(v3, 3, 0)/volume;
      out[1] += co_det(v3, 3, 1)/volume;
      out[2] += co_det(v3, 3, 2)/volume;
      count += 1.0;
    }
    if (d56 <= ZERO) {
      count += 1.0;
    }
    if (d64 <= ZERO) {
      out[0] += co_det(v2, 2, 0)/volume;
      out[1] += co_det(v2, 2, 1)/volume;
      out[2] += co_det(v2, 2, 2)/volume;
      count += 1.0;
    }
    if (d45 >= -ZERO && d56 >= -ZERO && d64 >= -ZERO) {
      out[0] += co_det(v4, 1, 0)/volume;
      out[1] += co_det(v4, 1, 1)/volume;
      out[2] += co_det(v4, 1, 2)/volume;
      count += 1.0;
    }
  } else if (idx == 5) {
    if (d45 <= ZERO) {
      out[0] += co_det(v3, 2, 0)/volume;
      out[1] += co_det(v3, 2, 1)/volume;
      out[2] += co_det(v3, 2, 2)/volume;
      count += 1.0;
    }
    if (d56 <= ZERO) {
      out[0] += co_det(v1, 3, 0)/volume;
      out[1] += co_det(v1, 3, 1)/volume;
      out[2] += co_det(v1, 3, 2)/volume;
      count += 1.0;
    }
    if (d64 <= ZERO) {
      count += 1.0;
    }
    if (d45 >= -ZERO && d56 >= -ZERO && d64 >= -ZERO) {
      out[0] += co_det(v4, 2, 0)/volume;
      out[1] += co_det(v4, 2, 1)/volume;
      out[2] += co_det(v4, 2, 2)/volume;
      count += 1.0;
    }
  } else {
    if (d45 <= ZERO) {
      count += 1.0;
    }
    if (d56 <= ZERO) {
      out[0] += co_det(v1, 2, 0)/volume;
      out[1] += co_det(v1, 2, 1)/volume;
      out[2] += co_det(v1, 2, 2)/volume;
      count += 1.0;
    }
    if (d64 <= ZERO) {
      out[0] += co_det(v2, 3, 0)/volume;
      out[1] += co_det(v2, 3, 1)/volume;
      out[2] += co_det(v2, 3, 2)/volume;
      count += 1.0;
    }
    if (d45 >= -ZERO && d56 >= -ZERO && d64 >= -ZERO) {
      out[0] += co_det(v4, 3, 0)/volume;
      out[1] += co_det(v4, 3, 1)/volume;
      out[2] += co_det(v4, 3, 2)/volume;
      count += 1.0;
    }
  }

  out[0] /= count;
  out[1] /= count;
  out[2] /= count;
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
void phi_6(const double * p, const double ** v, void * value) { vertex_quad_value(5, p, v, value); }
void phi_7(const double * p, const double ** v, void * value) { vertex_quad_value(6, p, v, value); }

void gradient_phi_1(const double * p, const double ** v, void * value) { vertex_quad_gradient(0, p, v, value); }
void gradient_phi_2(const double * p, const double ** v, void * value) { vertex_quad_gradient(1, p, v, value); }
void gradient_phi_3(const double * p, const double ** v, void * value) { vertex_quad_gradient(2, p, v, value); }
void gradient_phi_4(const double * p, const double ** v, void * value) { vertex_quad_gradient(3, p, v, value); }
void gradient_phi_5(const double * p, const double ** v, void * value) { vertex_quad_gradient(4, p, v, value); }
void gradient_phi_6(const double * p, const double ** v, void * value) { vertex_quad_gradient(5, p, v, value); }
void gradient_phi_7(const double * p, const double ** v, void * value) { vertex_quad_gradient(6, p, v, value); }

void psi_1(const double * p, const double ** v, void * value) { edge_quad_value(0, 1, p, v, value); }
void psi_2(const double * p, const double ** v, void * value) { edge_quad_value(0, 2, p, v, value); }
void psi_3(const double * p, const double ** v, void * value) { edge_quad_value(0, 3, p, v, value); }
void psi_4(const double * p, const double ** v, void * value) { edge_quad_value(2, 4, p, v, value); }
void psi_5(const double * p, const double ** v, void * value) { edge_quad_value(3, 4, p, v, value); }
void psi_6(const double * p, const double ** v, void * value) { edge_quad_value(3, 5, p, v, value); }
void psi_7(const double * p, const double ** v, void * value) { edge_quad_value(1, 5, p, v, value); }
void psi_8(const double * p, const double ** v, void * value) { edge_quad_value(1, 6, p, v, value); }
void psi_9(const double * p, const double ** v, void * value) { edge_quad_value(2, 6, p, v, value); }
void psi_10(const double * p, const double ** v, void * value) { edge_quad_value(4, 6, p, v, value); }
void psi_11(const double * p, const double ** v, void * value) { edge_quad_value(4, 5, p, v, value); }
void psi_12(const double * p, const double ** v, void * value) { edge_quad_value(5, 6, p, v, value); }

void gradient_psi_1(const double * p, const double ** v, void * value) { edge_quad_gradient(0, 1, p, v, value); }
void gradient_psi_2(const double * p, const double ** v, void * value) { edge_quad_gradient(0, 2, p, v, value); }
void gradient_psi_3(const double * p, const double ** v, void * value) { edge_quad_gradient(0, 3, p, v, value); }
void gradient_psi_4(const double * p, const double ** v, void * value) { edge_quad_gradient(2, 4, p, v, value); }
void gradient_psi_5(const double * p, const double ** v, void * value) { edge_quad_gradient(3, 4, p, v, value); }
void gradient_psi_6(const double * p, const double ** v, void * value) { edge_quad_gradient(3, 5, p, v, value); }
void gradient_psi_7(const double * p, const double ** v, void * value) { edge_quad_gradient(1, 5, p, v, value); }
void gradient_psi_8(const double * p, const double ** v, void * value) { edge_quad_gradient(1, 6, p, v, value); }
void gradient_psi_9(const double * p, const double ** v, void * value) { edge_quad_gradient(2, 6, p, v, value); }
void gradient_psi_10(const double * p, const double ** v, void * value) { edge_quad_gradient(4, 6, p, v, value); }
void gradient_psi_11(const double * p, const double ** v, void * value) { edge_quad_gradient(4, 5, p, v, value); }
void gradient_psi_12(const double * p, const double ** v, void * value) { edge_quad_gradient(5, 6, p, v, value); }

void chi_1(const double * p, const double ** v, void * value) { edge_quad_value(0, 4, p, v, value); }
void chi_2(const double * p, const double ** v, void * value) { edge_quad_value(0, 5, p, v, value); }
void chi_3(const double * p, const double ** v, void * value) { edge_quad_value(0, 6, p, v, value); }

void gradient_chi_1(const double * p, const double ** v, void * value) { edge_quad_gradient(0, 4, p, v, value); }
void gradient_chi_2(const double * p, const double ** v, void * value) { edge_quad_gradient(0, 5, p, v, value); }
void gradient_chi_3(const double * p, const double ** v, void * value) { edge_quad_gradient(0, 6, p, v, value); }

/*
 *  end of file
 *************************************************************************/
