import numpy as np


def make_gaussian_2d(x0=0, y0=0, theta=0, sigma_x=10, sigma_y=10):
    """
    Return a function of `(x, y)` that evaluates the value of a 2D Gaussian
    distribution centered on `(x0, y0)` and std.div. in x- and
    y-directions `sigma_x` and `sigma_y` rotated by angle `theta` in degrees
    """
    # x_center and y_center will be the center of the gaussian, theta will be the rotation angle
    # sigma_x and sigma_y will be the stdevs in the x and y axis before rotation
    # x_size and y_size give the size of the frame
    sx = sigma_x
    sy = sigma_y

    theta = theta * np.pi / 180

    def fn(x, y):
        a = np.cos(theta) * x - np.sin(theta) * y
        b = np.sin(theta) * x + np.cos(theta) * y
        a0 = np.cos(theta) * x0 - np.sin(theta) * y0
        b0 = np.sin(theta) * x0 + np.cos(theta) * y0

        return np.exp(
            -(((a - a0) ** 2) / (2 * (sx**2)) + ((b - b0) ** 2) / (2 * (sy**2)))
        )

    return fn


def make_carrier_wave_2d(theta, lw):
    """
    Return function of `(x,y)` which evaluates the amplitude of a plane-parallel sinosoidal with wavelength `lw` and orientation angle `theta` (in degrees)
    """
    theta = theta * np.pi / 180

    def fn(x, y):
        x_ = np.cos(theta) * x - np.sin(theta) * y
        return np.cos(2 * np.pi / lw * x_)

    return fn


def make_2d_wavepacket(x0, y0, lx, ly, theta, lw):
    """
    Return a function of `(x,y)` which produces the amplitude of a plane-parallel sinosoidal wave with wavelength `lw` and orientation angle `theta` modulated by a 2D Gaussian envelope centered on `(x0, y0)` of characteristic widths `(lx, ly)`
    oriented with width `lx` along the direction angle `theta`.

    In summary:
        (x0, y0): position of wavepacket
        (lx, ly): x- and y- length-scale of wavepacket
        theta: orientation of wavepacket and sinosoidal carrier wave [deg]
        lw: length-scale of gravity "carrier wave"
    """
    fn_envelope = make_gaussian_2d(
        x0=x0, y0=y0, theta=theta, sigma_x=lx / 2, sigma_y=ly / 2
    )
    fn_carrier = make_carrier_wave_2d(theta=theta, lw=lw)

    def fn(x, y):
        return fn_envelope(x, y) * fn_carrier(x, y)

    return fn
