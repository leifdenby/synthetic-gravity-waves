import numpy as np

import synthetic_gravity_waves as sgw


def test_point_sampling():
    r0 = 100
    r_sigma = 1.0
    (
        n_pts,
        pts,
        radii,
    ) = sgw.modified_poisson_disk_sampling(N=800, r0=r0, r_sigma=r_sigma)

    np.testing.assert_allclose(r0, np.mean(radii), rtol=r_sigma / r0)
    np.testing.assert_allclose(r_sigma, np.std(radii), rtol=0.2)


def test_wave_packet():
    fn_wavepacket = sgw.make_2d_wavepacket(
        x0=150, y0=200, lx=100, ly=30, theta=0.0, lw=20
    )
    N = 100
    x, y = np.meshgrid(np.arange(N), np.arange(N), indexing="ij")
    phi = fn_wavepacket(x, y)
    assert phi.shape == (N, N)


def test_gravity_wave_composite():
    N = 256
    phi, phi_envelope = sgw.make_synthetic_gravity_wave_composite(
        return_envelope=True, N=N
    )
    assert phi.shape == phi_envelope.shape == (N, N)


def test_version():
    assert sgw.__version__ is not None
