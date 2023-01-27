import numpy as np


def get_grid_coordinates(coords):
    return np.floor(coords).astype("int")


def modified_poisson_disk_sampling(N=100, r0=10, r_sigma=1, k=50, radiusType="default"):
    """
    Implementation of the Poisson Disk Sampling algorithm, but modified so that for each
    point added the radius between points is sampled from a Gaussian distribution

    :param N: grid size in number of pixels (assumed square domain)
    :param r0: length-scale for distance between points
    :param r_sigma: std div for distance between points
    :param k: Number of iterations to use to place a new point before giving up
    :return: nParticle: Number of particles in the sampling.
             particleCoordinates: 2d array containing the coordinates of the created particles.
             radii: radii for the sampled particles

    based on https://gitlab.com/abittner/poissondisksampling/-/blob/master/poissonDiskSampling/bridsonVariableRadius.py
    """

    def _gen_radius():
        return np.random.normal(r0, r_sigma)

    # Set-up background grid
    gridHeight = gridWidth = N
    grid = np.zeros((gridHeight, gridWidth))

    # Pick initial (active) point
    coords = (np.random.random() * gridHeight, np.random.random() * gridWidth)
    idx = get_grid_coordinates(coords)
    nParticle = 1
    grid[idx[0], idx[1]] = nParticle

    # Initialise active queue
    # Appending to list is much quicker than to numpy array, if you do it very often
    queue = [coords]
    # List containing the exact positions of the final particles
    particleCoordinates = [coords]
    radii = [_gen_radius()]
    activeRadii = [radii[0]]

    # Continue iteration while there is still points in active list
    while queue:

        # Pick random element in active queue
        idx = np.random.randint(len(queue))
        r_active = activeRadii[idx]
        activeCoords = queue[idx]

        success = False
        for _ in range(k):

            # Pick the angle to the sample particle and determine its coordinates
            angle = 2 * np.pi * np.random.random()
            newCoords = np.zeros(2)
            newCoords[0] = activeCoords[0] + 2 * r_active * np.sin(angle)
            newCoords[1] = activeCoords[1] + 2 * r_active * np.cos(angle)

            # Prevent that the new particle is outside of the grid
            if not (0 <= newCoords[1] <= gridWidth and 0 <= newCoords[0] <= gridHeight):
                continue

            # Check that particle is not too close to other particle
            newGridCoords = get_grid_coordinates((newCoords[1], newCoords[0]))

            radiusThere = _gen_radius()
            gridRangeX = (
                np.max([newGridCoords[0] - radiusThere, 0]).astype("int"),
                np.min([newGridCoords[0] + radiusThere + 1, gridWidth]).astype("int"),
            )
            gridRangeY = (
                np.max([newGridCoords[1] - radiusThere, 0]).astype("int"),
                np.min([newGridCoords[1] + radiusThere + 1, gridHeight]).astype("int"),
            )

            searchGrid = grid[
                slice(gridRangeY[0], gridRangeY[1]), slice(gridRangeX[0], gridRangeX[1])
            ]
            conflicts = np.where(searchGrid > 0)

            if len(conflicts[0]) == 0 and len(conflicts[1]) == 0:
                # No conflicts detected. Create a new particle at this position!
                queue.append(newCoords)
                activeRadii.append(radiusThere)
                radii.append(radiusThere)
                particleCoordinates.append(newCoords)
                nParticle += 1
                grid[newGridCoords[1], newGridCoords[0]] = nParticle
                success = True

            else:
                # There is a conflict. Do NOT create a new particle at this position!
                continue

        if not success:
            # No new particle could be associated to the currently active particle.
            # Remove current particle from the active queue!
            del queue[idx]
            del activeRadii[idx]

    return (nParticle, np.array(particleCoordinates), np.array(radii))
