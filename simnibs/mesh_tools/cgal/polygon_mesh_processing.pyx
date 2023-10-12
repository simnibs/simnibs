# from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.pair cimport pair
# import cython
from typing import Union
import numpy as np
import numpy.typing as npt
cimport numpy as np


cdef extern from "polygon_mesh_processing_src.cpp" nogil:
    vector[vector[int]] pmp_self_intersections(
        vector[vector[float]] vertices,
        vector[vector[int]] faces,
    )

    pair[vector[int], vector[int]] pmp_connected_components(
        vector[vector[float]] vertices,
        vector[vector[int]] faces,
        vector[int] contrained_faces,
    )

    # pair[vector[int], vector[int]] pmp_volume_connected_components(
    #     vector[vector[int]] faces,
    #     bool do_orientation_tests,
    #     bool do_self_intersection_tests,
    # )

    vector[vector[float]] pmp_smooth_shape(
        vector[vector[float]] vertices,
        vector[vector[int]] faces,
        vector[int] constrained_vertices,
        double time,
        int niter
    )

    # vector[vector[float]] pmp_angle_and_area_smoothing(
    #     vector[vector[float]] vertices,
    #     vector[vector[int]] faces,
    #     vector[int] constrained_vertices,
    #     int niter,
    #     bool use_angle_smoothing,
    #     bool use_area_smoothing,
    #     bool use_delaunay_flips,
    #     bool use_safety_constraints
    # )

    # vector[vector[float]] pmp_fair(
    #     vector[vector[float]] vertices,
    #     vector[vector[int]] faces,
    #     vector[int] indices,
    # )

    pair[vector[vector[float]], vector[vector[int]]] pmp_isotropic_remeshing(
        vector[vector[float]] vertices,
        vector[vector[int]] faces,
        double target_edge_length,
        int n_iterations,
    )

    # pair[vector[vector[float]], vector[vector[int]]] pmp_corefine_and_union(
    #     vector[vector[float]] vertices1,
    #     vector[vector[int]] faces1,
    #     vector[vector[float]] vertices2,
    #     vector[vector[int]] faces2,
    # )


def self_intersections(vertices: npt.ArrayLike, faces: npt.ArrayLike) -> npt.NDArray:
    """Compute the intersecting pairs of triangles in a surface mesh.

    Parameters
    ----------
    vertices : npt.ArrayLike

    faces : npt.ArrayLike

    Returns
    -------
    intersecting_pairs : npt

    """
    cdef np.ndarray[float, ndim=2] cpp_v = np.ascontiguousarray(vertices, dtype=np.float32)
    cdef np.ndarray[int, ndim=2] cpp_f = np.ascontiguousarray(faces, dtype=np.int32)

    cdef vector[vector[int]] intersecting_pairs # list of lists

    intersecting_pairs = pmp_self_intersections(cpp_v, cpp_f)

    return np.array(intersecting_pairs, dtype=int)

def connected_components(
        vertices: npt.ArrayLike,
        faces: npt.ArrayLike,
        constrained_faces: Union[npt.ArrayLike, None] = None
    ) -> tuple[npt.NDArray, npt.NDArray]:
    """Label connected components on a surface (graph).

    Parameters
    ----------
    vertices: npt.ArrayLike
        Vertices of the surfaces, shape = (N, 3).
    faces: npt.ArrayLike
        Faces of the surface, shape = (M, 3).
    constrained_faces: Union[npt.ArrayLike, None]
        The relevant function in CGAL takes a set of *edges* to constrain,
        however, this interface allows specifying *faces* instead such that
        only the *outer* edges of these faces are used as constraints. By
        default, no faces (edges) are constrained.

    Returns
    -------
    component_label : npt.NDArray
        The label associated with each face.
    component_size : npt.NDArray
        The size associated with each label.
    """
    cdef np.ndarray[float, ndim=2] cpp_v = np.ascontiguousarray(vertices, dtype=np.float32)
    cdef np.ndarray[int, ndim=2] cpp_f = np.ascontiguousarray(faces, dtype=np.int32)
    cdef np.ndarray[int] cpp_constrained_faces = np.ascontiguousarray(constrained_faces or [], dtype=np.int32)
    cdef pair[vector[int], vector[int]] out

    out = pmp_connected_components(cpp_v, cpp_f, cpp_constrained_faces)
    component_label = np.array(out.first, dtype=int)
    component_size = np.array(out.second, dtype=int)

    return component_label, component_size


# def volume_connected_components(faces, do_orientation_tests = False, do_self_intersection_tests = False):
#     cdef np.ndarray[int, ndim=2] cpp_f = np.ascontiguousarray(faces, dtype=np.int32)
#     cdef bool cpp_do_orientation_tests = do_orientation_tests
#     cdef bool cpp_do_self_intersection_tests = do_self_intersection_tests

#     cdef pair[vector[int], vector[int]] out

#     out = pmp_volume_connected_components(
#         cpp_f,
#         cpp_do_orientation_tests,
#         cpp_do_self_intersection_tests
#     )
#     component_label = np.array(out.first, dtype=int)
#     component_size = np.array(out.second, dtype=int)

#     return component_label, component_size


def smooth_shape(
        vertices: npt.ArrayLike,
        faces: npt.ArrayLike,
        constrained_vertices: Union[npt.ArrayLike, None] = None,
        time: float = 0.01,
        niter: int = 10
    ) -> npt.NDArray:
    """Shape smoothing using mean curvature flow.

    Parameters
    ----------
    time : float
        Determines the step size in the smoothing procedure (higher values
        means more aggressive smoothing).


    Returns
    -------
    The smoothed vertices.

    References
    ----------
    https://doc.cgal.org/latest/Polygon_mesh_processing/
    https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__meshing__grp.html#ga57fa999abe8dc557003482444df2a189
    """
    cdef np.ndarray[float, ndim=2] cpp_v = np.ascontiguousarray(vertices, dtype=np.float32)
    cdef np.ndarray[int, ndim=2] cpp_f = np.ascontiguousarray(faces, dtype=np.int32)
    cdef np.ndarray[int] cpp_constrained_vertices = np.ascontiguousarray(constrained_vertices or [], dtype=np.int32)
    cdef vector[vector[float]] v

    v = pmp_smooth_shape(cpp_v, cpp_f, cpp_constrained_vertices, time, niter)

    return np.array(v, dtype=float)


# def angle_and_area_smoothing(
#         vertices: npt.ArrayLike,
#         faces: npt.ArrayLike,
#         constrained_vertices: Union[npt.ArrayLike, None] = None,
#         niter: int = 10,
#         use_angle_smoothing: bool = True,
#         use_area_smoothing: bool = True,
#         use_delaunay_flips: bool = True,
#         use_safety_constraints: bool = False
#     ):
#     """Vertex smoothing preserving shape.

#     Parameters
#     ----------
#     vertices: npt.ArrayLike
#     faces: npt.ArrayLike
#     constrained_vertices : Union[npt.ArrayLike, None]
#     niter: int
#     use_angle_smoothing: bool = True
#     use_area_smoothing: bool = True
#     use_delaunay_flips: bool = True
#     use_safety_constraints: bool = False

#     Returns
#     -------


#     """
#     cdef np.ndarray[float, ndim=2] cpp_v = np.ascontiguousarray(vertices, dtype=np.float32)
#     cdef np.ndarray[int, ndim=2] cpp_f = np.ascontiguousarray(faces, dtype=np.int32)
#     cdef np.ndarray[int] cpp_constrained_vertices = np.ascontiguousarray(constrained_vertices or [], dtype=np.int32)
#     cdef vector[vector[float]] v

#     v = pmp_angle_and_area_smoothing(
#         cpp_v,
#         cpp_f,
#         cpp_constrained_vertices,
#         niter,
#         use_angle_smoothing,
#         use_area_smoothing,
#         use_delaunay_flips,
#         use_safety_constraints
#     )
#     return np.array(v, dtype=float)


# def fair(vertices, faces, vertex_indices):
#     """Mesh fairing."""
#     cdef np.ndarray[float, ndim=2] cpp_v = np.ascontiguousarray(vertices, dtype=np.float32)
#     cdef np.ndarray[int, ndim=2] cpp_f = np.ascontiguousarray(faces, dtype=np.int32)
#     cdef np.ndarray[int] cpp_vi = np.ascontiguousarray(vertex_indices, dtype=np.int32)
#     cdef vector[vector[float]] v

#     v = pmp_fair(cpp_v, cpp_f, cpp_vi)

#     return np.array(v, dtype=float)


def isotropic_remeshing(
        vertices: npt.ArrayLike,
        faces: npt.ArrayLike,
        target_edge_length: double,
        n_iterations: int = 1,
    ):
    """Isotropic surface remeshing. Remeshing is achieved by a combination of
    edge splits/flips/collapses, tangential relaxation, and projection back
    onto the original surface.

    Parameters
    ----------
    vertices: npt.ArrayLike
    faces: npt.ArrayLike
    target_edge_length: double
        The target edge length for the isotropic remesher. This defines the
        resolution of the resulting surface.
    n_iterations: int
        Number of iterations of the above-mentioned atomic operations.

    Returns
    -------
    v : npt.NDArray
        The new vertices.
    f : npt.NDArray
        The new faces.

    References
    ----------

    https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__meshing__grp.html#gaa5cc92275df27f0baab2472ecbc4ea3f

    """
    cdef np.ndarray[float, ndim=2] cpp_v = np.ascontiguousarray(vertices, dtype=np.float32)
    cdef np.ndarray[int, ndim=2] cpp_f = np.ascontiguousarray(faces, dtype=np.int32)
    cdef pair[vector[vector[float]], vector[vector[int]]] out

    out = pmp_isotropic_remeshing(
        cpp_v, cpp_f, target_edge_length, n_iterations
    )
    v = np.array(out.first, dtype=float)
    f = np.array(out.second, dtype=int)

    return v, f


# def corefine_and_union(vertices1, faces1, vertices2, faces2):
#     cdef np.ndarray[float, ndim=2] cpp_v1 = np.ascontiguousarray(vertices1, dtype=np.float32)
#     cdef np.ndarray[int, ndim=2] cpp_f1 = np.ascontiguousarray(faces1, dtype=np.int32)
#     cdef np.ndarray[float, ndim=2] cpp_v2 = np.ascontiguousarray(vertices2, dtype=np.float32)
#     cdef np.ndarray[int, ndim=2] cpp_f2 = np.ascontiguousarray(faces2, dtype=np.int32)
#     cdef pair[vector[vector[float]], vector[vector[int]]] out

#     out = pmp_corefine_and_union(cpp_v1, cpp_f1, cpp_v2, cpp_f2)
#     v = np.array(out.first, dtype=float)
#     f = np.array(out.second, dtype=int)

#     return v, f
