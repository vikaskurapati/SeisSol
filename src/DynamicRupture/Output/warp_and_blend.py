import numpy as np
import matplotlib.pyplot as plt


def equidistant_points(N, corners):
    M = (N+2) * (N+1) // 2
    cartesian_coords = np.zeros((2, M))
    barycentric_coords = np.zeros((3, M))
    k = 0
    for i in range(N+1):
        for j in range(N-i+1):
            barycentric_coords[0, k] = i / N
            barycentric_coords[1, k] = j / N
            barycentric_coords[2, k] = 1 - barycentric_coords[0, k] - barycentric_coords[1, k]
            cartesian_coords[:,k] = np.dot(corners, barycentric_coords[:,k])
            k += 1
    return cartesian_coords, barycentric_coords

def warp_and_blend(barycentric_coords, corners):
    w = lambda x: -np.sin(0.5*x*np.pi) / (1 - x**2)
    result = np.zeros((3, 2, barycentric_coords.shape[1]))

    def w_a_b(barycentric_coords, corners, i, j):
        warp = w(barycentric_coords[i,:] - barycentric_coords[j,:])
        blend = 4 * barycentric_coords[i,:] * barycentric_coords[j,:]
        vector = corners[:,i] - corners[:,j]
        vector = vector / np.linalg.norm(vector)
        return np.multiply(blend, np.outer(warp, vector).T)

    result[0,:,:] = w_a_b(barycentric_coords, corners, 2, 1)
    result[1,:,:] = w_a_b(barycentric_coords, corners, 0, 2)
    result[2,:,:] = w_a_b(barycentric_coords, corners, 1, 0)

    return result

corners = np.array([[0, 1, 0.1], [0.1, 0.2, 1]])
points, barycentric_coords = equidistant_points(10, corners)
warp = warp_and_blend(barycentric_coords, corners)
fig, axes = plt.subplots(2, 4)
axes[0, 0].scatter(points[0,:], points[1,:])
axes[0, 0].plot([corners[0,0], corners[0,1], corners[0,2], corners[0,0]], [corners[1,0], corners[1,1], corners[1,2], corners[1,0]])
axes[0, 1].scatter(points[0,:], points[1,:], c=barycentric_coords[0,:])
axes[0, 1].set_title("l_0")
axes[0, 2].scatter(points[0,:], points[1,:], c=barycentric_coords[1,:])
axes[0, 2].set_title("l_1")
axes[0, 3].scatter(points[0,:], points[1,:], c=barycentric_coords[2,:])
axes[0, 3].set_title("l_2")

for k in range(3):
    axes[1, k+1].quiver(points[0,:], points[1,:], warp[k, 0,:], warp[k, 1,:])
axes[1, 0].quiver(points[0,:], points[1,:], np.sum(warp[:, 0,:], axis=0), np.sum(warp[:,1,:], axis=0))

for ax in axes.flatten():
    ax.set_aspect(1)
plt.show()

