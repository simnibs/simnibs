# -*- coding: utf-8 -*-
"""
Created on Sat Jan 21 15:33:52 2023

@author: khm
"""
import numpy as np
import json
from simnibs.simulation import coil_numpy
from scipy.spatial.transform import Rotation


def read_tcd(fn):
    with open(fn, 'r') as fid:
        coil = json.loads(fid.read())

    for i, cElm in enumerate(coil['elementList']):
        coil['elementList'][i]['points'] = np.array(cElm['points'])
        coil['elementList'][i]['values'] = np.array(cElm['values'])
        for k in coil['coilModel'][i].keys():
            try:
                coil['coilModel'][i][k] = np.array(coil['coilModel'][i][k])
            except:
                pass

    for j, deform in enumerate(coil['deformationList']):
        for k in deform.keys():
            if isinstance(deform[k], list):
                coil['deformationList'][j][k] = np.array(deform[k])

    return coil


def write_tcd(coildict, fname):
    coil = coildict.copy()
    for i, cElm in enumerate(coil['elementList']):
        if isinstance(coil['elementList'][i]['points'],np.ndarray):
            coil['elementList'][i]['points'] = cElm['points'].tolist()
        if isinstance(coil['elementList'][i]['values'],np.ndarray):
            coil['elementList'][i]['values'] = cElm['values'].tolist()
        for k in coil['coilModel'][i].keys():
            try:
                coil['coilModel'][i][k] = coil['coilModel'][i][k].tolist()
            except:
                pass

    for j, deform in enumerate(coil['deformationList']):
        for k in deform.keys():
            if isinstance(deform[k], np.ndarray):
                coil['deformationList'][j][k] = deform[k].tolist()

    json_tcd = json.dumps(coil)
    with open(fname, 'w') as fid:
        fid.write(json_tcd)


def rot2p(p1, p2, a, unit='deg', affine=None):
    p1 = np.asanyarray(p1, dtype='float64')
    p2 = np.asanyarray(p2, dtype='float64')
    if affine is not None:
        p1 = affine @ p1
        p2 = affine @ p2
    p = p2 - p1
    p /= np.sqrt(np.sum(p**2))
    T = np.identity(4)
    T[:3, 3] = p1
    iT = np.identity(4)
    iT[:3, 3] = -p1
    R = np.identity(4)
    R[:3, :3] = Rotation.from_rotvec(p * a, degrees=unit == 'deg').as_matrix()
    Q = T @ R @ iT
    return Q


def get_transform(deform, parameter, affine=None):
    if deform['type'][:3] == 'rot':
        if deform['type'] == 'rot2p':
            R = rot2p(deform['point1'], deform['point2'],
                      parameter, affine=None)
        else:
            R = np.identity(4)
            R[:3, :3] = Rotation.from_euler(deform['type'][3:], parameter,
                                            degrees=deform['unit'] == 'deg'
                                            ).as_matrix()
            if not deform['affine'] is None:
                iA = np.linalg.inv(deform['affine'])
                R = iA @ R @ deform['affine']
    else:
        R = np.identity(4)
        i = ['x', 'y', 'z'].index(deform['type'])
        R[i, 3] = parameter
    return R


def transform_points(points, R):
    return points @ R[:3, :3].T + R[None, :3, 3]


def transform_vectors(vectors, R):
    return vectors @ R[:3, :3].T


def combine_transforms(deforms, parameters, idx=None):
    affine = None
    for i, p in enumerate(parameters):
        if idx is None:
            j = i
        else:
            j = idx[i]
        R = get_transform(deforms[j], p, affine=affine)
        if affine is None:
            affine = R
        else:
            affine = R @ affine
    return affine


def get_all_transforms(coil, parameters):
    transforms = []
    for j, cElm in enumerate(coil['elementList']):
        transforms.append(combine_transforms(coil['deformationList'],
                                             parameters[cElm['deformations']],
                                             idx=cElm['deformations']))
    return transforms


def get_model_points(coil, parameters,
                     fieldnames=['points', 'minDistancePoints'],
                     affine=None, collapse=True):
    
    if parameters is None:
        parameters = np.zeros(len(coil['deformationList']))
    parameters = np.asanyarray(parameters)
    points = []
    for k, f in enumerate(fieldnames):
        points.append([])
    for j, cElm in enumerate(coil['elementList']):
        R = combine_transforms(coil['deformationList'],
                               parameters[cElm['deformations']],
                               idx=cElm['deformations'])
        if affine is not None:
            R = affine @ R
        for k, f in enumerate(fieldnames):
            points[k].append(transform_points(coil['coilModel'][j][f], R))

    for k, f in enumerate(fieldnames):
        if collapse:
            points[k] = np.concatenate(points[k], axis=0)
    return points


def get_params(coil):
    out = []
    for deform in coil['deformationList']:
        out.append(f'{deform["type"]}')
    return out


def transform_coil(coil, parameters,
                   fieldnames=['points', 'minDistancePoints']):
    parameters = np.asanyarray(parameters)
    for j, cElm in enumerate(coil['elementList']):
        R = combine_transforms(coil['deformationList'],
                               parameters[cElm['deformations']],
                               idx=cElm['deformations'])
        cElm['points'] = transform_points(cElm['points'], R)
        if 'values' in cElm:
            cElm['values'] = transform_vectors(cElm['values'], R)
        if 'coilModel' in coil and len(coil['coilModel']) > j:
            for k, f in enumerate(fieldnames):
                if f in coil['coilModel'][j]:
                    coil['coilModel'][j][f] = transform_points(
                        coil['coilModel'][j][f], R)


def get_Afield(coil, positions, parameters=None, affine=None, dIdt=1.0, eps=1e-3):
    # positions in SI units (meters)
    # coil wire positions are converted to meters *1e-3
    points = []
    values = []
    if affine is None:
        affine = np.identity(4)
    for j, cElm in enumerate(coil['elementList']):
        if parameters is not None:
            R = affine @ combine_transforms(coil['deformationList'],
                                            parameters[cElm['deformations']],
                                            idx=cElm['deformations'])
        else:
            R = affine
            points.append(transform_points(cElm['points'], R))
            if 'values' in cElm:
                values.append(transform_vectors(cElm['values'], R))

    types = [cElm['type'] for cElm in coil['elementList']]
    if np.all(np.array(types) == types[0]):
        points = [np.concatenate(points, axis=0)]
        if len(values) > 0:
            values = [np.concatenate(values, axis=0)]

    types = [cElm['type'] for cElm in coil['elementList']]
    if np.all(np.array(types) == types[0]):
        points = [np.concatenate(points, axis=0)]
    A = np.zeros((positions.T.shape))
    k = 0
    for i, p in enumerate(points):
        # dipole case
        if coil['elementList'][i]['type'] == 1:
            A += dIdt * coil_numpy.A_from_dipoles(values[k], p*1e-3,
                                                  positions.T, eps=eps)
            k += 1
        # line integral with directions
        elif coil['elementList'][i]['type'] == 2:
            A += coil_numpy.A_biot_savart_path_fmm(p.T*1e-3,
                                                   positions, values[k].T*1e-3,
                                                   eps=eps, I=dIdt).T
            k += 1
        # line integral without directions
        elif coil['elementList'][i]['type'] == 3:
            A += coil_numpy.A_biot_savart_path_fmm(p.T*1e-3, positions,
                                                   eps=eps, I=dIdt).T
    return A


def get_costs(coil, msh, affine):
    tree = msh.get_AABBTree()
    trans = lambda x: x@affine[:3, :3].T + affine[:3, 3][None]
    subcost = lambda x: [tree.any_point_inside(trans(x[0])),
                         np.mean(np.sqrt(tree.min_sqdist(trans(x[1]))))]
    if 'intersectPoints' in coil['coilModel'][0]:
        cost = lambda x: subcost(
            get_model_points(coil, x, fieldnames=['intersectPoints',
                                                  'minDistancePoints']))
    else:
        cost = lambda x: subcost(
            get_model_points(coil, x, fieldnames=['points',
                                                  'minDistancePoints']))
    return cost

# def get_dist(coil, tree):
#     cost = lambda x:  np.sum(tree.min_sqdist(get_model_points(coil, x, fieldnames=['minDistancePoints'])[0]))
#     return cost

# def get_inside(coil, tree):
#     cost = lambda x: tree.any_point_inside(get_model_points(coil, x, fieldnames=['points'])[0])
#     return cost


def optimize_position(coil, skin, affine, M=None):
    cost = get_costs(coil, skin, affine)

    bounds = []
    for deform in coil['deformationList']:
        bounds.append(deform['range'])

    if M is None:
        p = np.array([d['initial'] for d in coil['deformationList']])

    constr, best_f = cost(p)
    while constr:
        p[-1] -= 2
        constr, best_f = cost(p)

    while not constr:
        p[-1] += 1
        constr, best_f = cost(p)
        if constr:
            p[-1] -= 2
            break

    M = np.abs(cost(p)[1])
    print(M)
    best_x = p

    def costf_x0(x, x0):
        fm = cost(x0 + x)
        f = M * fm[0] + fm[1]
        if not fm[0]:
            nonlocal best_f
            if f < best_f:
                if np.all([xi >= bounds[i][0] and
                           xi <= bounds[i][1] for i, xi in enumerate(x0 + x)]):
                    nonlocal best_x
                    best_x = x0 + x
                    best_f = f
        return f
    costf = lambda x: costf_x0(x, p)
    best_f = costf(p)

    def printsol(*xk):
        xk = np.round(np.array(p+xk).ravel(), 1)
        nonlocal best_f
        print(f'Testing parameter setting: {xk}, best cost: {best_f:.2f}')

    import scipy.optimize as opt
    res = opt.direct(costf, bounds=bounds, callback=printsol)
    xi = best_x.copy()
    xi[-1] -= 1
    print(best_x)
    costf = lambda x: costf_x0(x, 0)
    res2 = opt.minimize(costf, x0=xi, callback=printsol,
                        method='L-BFGS-B',
                        options={'eps': 0.001, 'maxls': 100}, bounds=bounds)

    print(cost(res2.x))
    print(cost(best_x))

    res2.x = best_x
    res2.fun = best_f
    print(res2)

    # refine univariate
    for i in range(len(res2.x)):
        cost1 = lambda xx: costf(np.concatenate(
            (res2.x[:i], [xx], res2.x[i+1:]), axis=None))
        res = opt.minimize(cost1, x0=res2.x[i], method='L-BFGS-B',
                           options={'eps': 0.001, 'maxls': 100},
                           bounds=[bounds[i]])
        print(best_x)
    print(cost(best_x))
    res2.x = best_x
    return res2


def get_coil_msh(coil, affine, parameters=None, filename=None):
    import simnibs.mesh_tools
    points = get_model_points(coil, parameters, fieldnames=['points'],
                              affine=affine, collapse=False)[0]
    coilmsh = simnibs.mesh_tools.mesh_io.Msh()
    for i, cm in enumerate(coil['coilModel']):
        coilmsh = coilmsh.join_mesh(
            simnibs.mesh_tools.mesh_io.Msh(
                nodes=simnibs.mesh_tools.mesh_io.Nodes(points[i]),
                elements=simnibs.mesh_tools.mesh_io.Elements(
                triangles=cm['cells'] + 1)))
    if filename is not None:
        simnibs.mesh_tools.mesh_io.write_stl(coilmsh, filename)
    return coilmsh


def vizCoil(coil, affine, parameters, msh, filename=None):
    import pyvista as pv
    points = get_model_points(coil, parameters, fieldnames=['points'],
                              affine=affine, collapse=False)[0]
    cells = [cm['cells'] for cm in coil['coilModel']]
    # cells = np.concatenate(cells,axis=0)
    pcoil = []
    for i, cell in enumerate(cells):
        f = np.empty((cells[i].shape[0], 4), dtype='uint32')
        f[:, 1:] = cells[i].astype('uint32', copy=False)
        f[:, 0] = 3

        pcoil.append(pv.PolyData(points[i], f))

    f = msh.elm.node_number_list[:, [3, 0, 1, 2]].astype(
        'uint32', copy=False, order='C') - 1
    f[:, 0] = 3
    pmsh = pv.PolyData(msh.nodes.node_coord, f)

    if filename is None:
        p = pv.Plotter()
    else:
        p = pv.Plotter(off_screen=True)
    p.set_background('w')
    p.add_mesh(pmsh, color='gray')
    for pc in pcoil:
        p.add_mesh(pc, color='blue')
    if filename is None:
        p.show()
    else:
        p.show(screenshot=filename)


def aniCoil(coil, affine, parameters, msh):
    import pyvista as pv
    import time
    points = get_model_points(coil, parameters[0], fieldnames=['points'],
                              affine=affine, collapse=False)[0]
    cells = [cm['cells'] for cm in coil['coilModel']]
    # cells = np.concatenate(cells,axis=0)
    pcoil = []
    for i, cell in enumerate(cells):
        f = np.empty((cells[i].shape[0], 4), dtype='uint32')
        f[:, 1:] = cells[i].astype('uint32', copy=False)
        f[:, 0] = 3

        pcoil.append(pv.PolyData(points[i], f))

    f = msh.elm.node_number_list[:, [3, 0, 1, 2]].astype(
        'uint32', copy=False, order='C') - 1
    f[:, 0] = 3
    pmsh = pv.PolyData(msh.nodes.node_coord, f)
    actors = []
    p = pv.Plotter()
    p.set_background('w')
    p.add_mesh(pmsh, color='gray')
    for j, pc in enumerate(pcoil):
        actors.append(p.add_mesh(pc, color='blue', name=f'c{j}'))
    running = True

    def quit_func():
        nonlocal running
        running = False

    p.add_key_event('q', quit_func)
    p.show(interactive_update=True)
    i = 0

    while True:
        i += 1
        if i >= len(parameters):
            i = 0
        points = get_model_points(coil, parameters[i], fieldnames=['points'],
                                  affine=affine, collapse=False)[0]
        for j, pp in enumerate(points):
            # p.update_coordinates(pp, pcoil[j],render=False)
            pcoil[j].points = pp
            # p.add_mesh(pcoil[j], color='blue',name=f'c{j}')
        # p.render()
        p.update(10, True)

        print(parameters[i])
        time.sleep(.05)
        if not running:
            break
    p.close()
