import numpy as np

from .components import make_2d_wavepacket, make_gaussian_2d
from .point_sampling import modified_poisson_disk_sampling


def make_synthetic_gravity_wave_composite(
    N=512, r_sigma=5, lw0=20.0, r0=100, return_envelope=False
):
    """
    On a 2D unit-grid of `(N,N)` points create a amplitude composite of
    randomly oriented synthetic gravity waves with each gravity wave packet
    having a 2D Gaussian envelope with characteristic size `r0` (and spaced by
    `2*r0` to aviod overlap) and carrier wave wavelength uniformly sampled in
    `[0, lw0]` and envelope size aspect ratio normally distributed (mean=0.5,
    sigma=0.5).
    """
    x, y = np.meshgrid(np.arange(N), np.arange(N), indexing="xy")

    n_pts, pts, radii = modified_poisson_disk_sampling(N=N, r0=r0, r_sigma=r_sigma)

    phi = np.zeros((N, N))
    if return_envelope:
        phi_envelope = np.zeros((N, N))

    for pt, r in zip(pts, radii):
        theta = 180 * np.random.uniform()
        lw = lw0 * np.random.normal()
        a = 0.5 + 0.5 * np.random.uniform()

        lx = r / 2
        ly = a * r / 2
        x0, y0 = pt

        if return_envelope:
            fn_envelope = make_gaussian_2d(
                x0=x0, y0=y0, theta=theta, sigma_x=lx / 2, sigma_y=ly / 2
            )
            phi_envelope += fn_envelope(x, y)

        fn_wavepacket = make_2d_wavepacket(
            x0=x0, y0=y0, lx=lx, ly=ly, theta=theta, lw=lw
        )
        phi += fn_wavepacket(x, y)

    if return_envelope:
        return phi, phi_envelope
    else:
        return phi
