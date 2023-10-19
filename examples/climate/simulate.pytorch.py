import torch

def simulate(mesh, operators, hodge=None):
    if hodge is None:
        hodge = GeometricHodge()

    # Extracting operators and matrices
    M_d0, d0 = default_dec_matrix_generate(mesh, "d0", hodge)
    M_star_1, star_1 = default_dec_matrix_generate(mesh, "star_1", hodge)
    M_dual_d1, dual_d1 = default_dec_matrix_generate(mesh, "dual_d1", hodge)
    M_star_0_inv, star_0_inv = default_dec_matrix_generate(mesh, "star_0_inv", hodge)
    sharp = operators(mesh, "sharp")
    mag = operators(mesh, "mag")
    avg01 = operators(mesh, "avg01")
    power = operators(mesh, "power")

    # Initializing tensors
    dynamics_dot_1 = torch.empty(nparts(mesh, "E"))
    dynamics_dot_6 = torch.empty(nparts(mesh, "E"))
    dot_8_1 = torch.empty(nparts(mesh, "E"))
    dot_8_2 = torch.empty(nparts(mesh, "V"))
    dynamics_h_hat = torch.empty(nparts(mesh, "V"))

    def f(du, u, p, t):
        h = findnode(u, "h").values
        n = p.n
        stress_A = findnode(u, "stress_A").values
        stress_rho = p.stress_rho
        stress_g = p.stress_g
        one = torch.tensor(1.0)
        two = torch.tensor(2.0)

        dynamics_dot_1 = torch.matmul(M_d0, h)
        dynamics_dot_6 = torch.matmul(M_d0, h)
        dynamics_dot_5 = sharp(dynamics_dot_6)
        dynamics_dot_4 = mag(dynamics_dot_5)
        dynamics_dot_7 = n - one
        dynamics_dot_3 = torch.pow(dynamics_dot_4, dynamics_dot_7)
        stress_dot_3 = stress_rho * stress_g
        stress_dot_2 = torch.pow(stress_dot_3, n)
        dynamics_sum_1 = n + two
        stress_sum_1 = n + two
        dynamics_dot_2 = avg01(dynamics_dot_3)
        dynamics_dot_9 = torch.pow(h, dynamics_sum_1)
        stress_dot_1 = two / stress_sum_1
        stress_mult_1 = stress_dot_1 * stress_A
        gamma = stress_mult_1 * stress_dot_2
        dynamics_dot_8 = avg01(dynamics_dot_9)
        dynamics_mult_3 = gamma * dynamics_dot_1
        dynamics_mult_1 = dynamics_mult_3 * dynamics_dot_2
        dynamics_mult_2 = dynamics_mult_1 * dynamics_dot_8
        dot_8_1 = torch.matmul(M_star_1, dynamics_mult_2)
        dot_8_2 = torch.matmul(M_dual_d1, dot_8_1)
        dynamics_h_hat = torch.matmul(M_star_0_inv, dot_8_2)

        findnode(du, "h").values = dynamics_h_hat
        return du
    return f










