import sys
import os

import pytest
import numpy as np

sys.path.insert(0, os.path.abspath('..'))
import micromorphic

data = []

# Test the linear elastic material model
model_name = "LinearElasticity"
time = np.array([ 10, 2.7 ])

fparams = np.array([ 2, 1.7, 1.8, 5, 2.8, .76, .15, 9.8, 5.4, 11, 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 2, .76, 5.4])

current_grad_u = np.array([[ -1.07901185, -1.09656192, -0.04629144 ],\
                           [ -0.77749189, -1.27877771, -0.82648234 ],\
                           [  0.66484637, -0.05552567, -1.65125738 ]])

current_phi = np.array([-1.40391532, -0.42715691,  0.75393369,\
                         0.2849511 , -2.06484257, -0.52190902,\
                         1.07238446,  0.19155907, -0.39704566 ])

current_grad_phi = np.array([[ 0.14940184, 0.12460812, 0.31971128 ],\
                             [ 0.67550862, 0.61095383, 0.87972732 ],\
                             [ 0.30872424, 0.32158187, 0.25480281 ],\
                             [ 0.45570006, 0.69090695, 0.72388584 ],\
                             [ 0.14880964, 0.67520596, 0.15106516 ],\
                             [ 0.77810545, 0.07641724, 0.09367471 ],\
                             [ 0.15905979, 0.0651695 , 0.52150417 ],\
                             [ 0.91873444, 0.5622355 , 0.50199447 ],\
                             [ 0.26729942, 0.89858519, 0.09043229 ]])

previous_grad_u = np.array([[ 0, 0, 0],\
                            [ 0, 0, 0],\
                            [ 0, 0, 0]])

previous_phi = np.array([ 0, 0, 0,\
                          0, 0, 0,\
                          0, 0, 0])

previous_grad_phi = np.array([[ 0, 0, 0],\
                              [ 0, 0, 0],\
                              [ 0, 0, 0],\
                              [ 0, 0, 0],\
                              [ 0, 0, 0],\
                              [ 0, 0, 0],\
                              [ 0, 0, 0],\
                              [ 0, 0, 0],\
                              [ 0, 0, 0]]);


SDVS = np.zeros(3);
current_ADD_DOF, previous_ADD_DOF = np.zeros(3), np.zeros(3)

current_ADD_grad_DOF, previous_ADD_grad_DOF = np.zeros((3, 3)), np.zeros((3, 3))

PK2_answer = np.array([-26.78976487,  91.99831835, 135.04096376,
                       -63.68792655, 149.68226149, 186.67587146,
                       -42.54105342, 125.2317492 , 150.55767059])

SIGMA_answer = np.array([-47.59920949,  20.84881327,  93.02392773,
                          20.84881327, 302.43209139, 311.0104045 ,
                          93.02392773, 311.0104045 , 312.60512922])

M_answer = np.array([-50.37283054, -23.25778149, -37.92963077, -19.16962188,
                     -32.97279228, -14.89104497, -33.4026237 , -15.47947779,
                     -40.31460994, -16.29637436, -36.63942799, -18.22777296,
                     -39.33546661, -86.69472439, -59.29150146, -15.76480164,
                     -55.42039768, -35.09720118, -28.94394503, -17.96726082,
                     -45.09734176, -16.46568416, -50.79898863, -39.19129183,
                     -47.46372724, -42.98201472, -45.57864883])

answers = {'PK2':PK2_answer, 'SIGMA':SIGMA_answer, 'M':M_answer}

data += [(model_name, time, fparams, current_grad_u, current_phi, current_grad_phi, previous_grad_u, previous_phi, previous_grad_phi, SDVS,\
          current_ADD_DOF, previous_ADD_DOF, current_ADD_grad_DOF, previous_ADD_grad_DOF, answers)]

# Test the micromorphic elasto plastic model
model_name = "LinearElasticityDruckerPragerPlasticity"
time = np.array([ 10., 2.5 ])

fparams = np.array([2, 2.4e2, 1.5e1,
                    2, 1.4e2, 2.0e1,
                    2, 2.0e0, 2.7e1,
                    2, 0.56, 0.2,
                    2, 0.15,-0.2,
                    2, 0.82, 0.1,
                    2, 0.70, 0.3,
                    2, 0.40,-0.3,
                    2, 0.52, 0.4,
                    2, 696.47, 65.84,
                    5, -7.69, -51.92, 38.61, -27.31, 5.13,
                    11, 1.85, -0.19, -1.08, -1.57, 2.29, -0.61, 5.97, -2.02, 2.38, -0.32, -3.25,
                    2, -51.92, 5.13,
                    0.4, 0.3, 0.35, 1e-8, 1e-8,
                   ])

current_grad_u = np.array([[0.200, 0.100, 0.000 ],
                           [0.100, 0.001, 0.000 ],
                           [0.000, 0.000, 0.000 ]])

previous_grad_u = np.array([[0, 0, 0],
                            [0, 0, 0],
                            [0, 0, 0]])

current_phi = np.array([0.100, 0.000, 0.000,
                        0.000, 0.000, 0.000,
                        0.000, 0.000, 0.000])

previous_phi = np.array([ 0, 0, 0, 0, 0, 0, 0, 0, 0 ])

current_grad_phi = np.array([[ 0.13890017, -0.3598602 , -0.08048856 ],
                             [-0.18572739,  0.06847269,  0.22931628 ],
                             [-0.01829735, -0.48731265, -0.25277529 ],
                             [ 0.26626212,  0.4844646 , -0.31965177 ],
                             [ 0.49197846,  0.19051656, -0.0365349  ],
                             [-0.06607774, -0.33526875, -0.15803078 ],
                             [ 0.09738707, -0.49482218, -0.39584868 ],
                             [-0.45599864,  0.08585038, -0.09432794 ],
                             [ 0.23055539,  0.07564162,  0.24051469 ]])

previous_grad_phi = np.array([[0, 0, 0],
                              [0, 0, 0],
                              [0, 0, 0],
                              [0, 0, 0],
                              [0, 0, 0],
                              [0, 0, 0],
                              [0, 0, 0],
                              [0, 0, 0],
                              [0, 0, 0]])

SDVS = np.zeros([55])

PK2_answer = np.array([ 172.484,   15.3785,   -0.917177,
                         13.4848, 142.823,    -0.0214307,
                         -1.7635,   1.77719, 141.069 ])

SIGMA_answer = np.array([ 176.916,   15.8646,   -2.83731,
                           15.8646, 144.538,     1.85836,
                           -2.83731,  1.85836, 142.013 ])

M_answer = np.array([ 0.598283, -0.512218,  0.620664,    3.22636,   1.16682,
                      1.20593,   0.562825, -2.52317,     1.62616,  -2.61391,
                     -0.60994,  -1.02147,   0.668187,    0.49348,  -0.23916,
                     -2.77419,   0.760483,  1.71784,    -0.499389,  2.62828,
                     -0.761044,  1.23369,  -0.00778206, -2.25643,  -0.729551,
                      0.743204,  0.910521 ])

SDVS_answer = np.array([ -1.79592e-24,  0.0243222,    0.0822384,    0.0430345,   0.0435752,
                         -8.96006e-25,  0.00852191,   0.0465339,    0.0243507,   0.0246566,
                          0.00742998,   0.00500421,  -0.000296486,  0.00498757, -0.00260492,
                          0.000284355, -0.000367318,  0.000222511, -0.0015603,   0.00863313,
                          0.00537105,  -0.000347686,  0.00643802,  -0.00298667,  0.000297105,
                         -0.000422398,  0.000293946, -0.00181986,   0.0385522,  -0.0244021,
                         -0.00912035,   0.0105171,   -0.000328615,  0.0222843,   0.0185626,
                         -0.0118234,   -0.00785555,   0.0451085,    0.031607,    0.0212748,
                          0.0116981,    0.0161821,    0.00126031,   0.0147688,   0.00691805,
                         -0.0241431,    0.00942608,  -0.0350366,    0.0221571,  -0.0249697,
                          0.00935849,   0.0214931,    0.0169609,    0.0177352,   0.0203838 ])

answers = {'PK2':PK2_answer, 'SIGMA':SIGMA_answer, 'M':M_answer, 'SDVS':SDVS_answer}

data += [(model_name, time, fparams, current_grad_u, current_phi, current_grad_phi, previous_grad_u, previous_phi, previous_grad_phi, SDVS,\
          current_ADD_DOF, previous_ADD_DOF, current_ADD_grad_DOF, previous_ADD_grad_DOF, answers)]

@pytest.mark.parametrize("model_name, time, fparams, current_grad_u, current_phi, current_grad_phi, previous_grad_u, previous_phi, previous_grad_phi, SDVS, current_ADD_DOF, previous_ADD_DOF, current_ADD_grad_DOF, previous_ADD_grad_DOF, answers", data)
def test_evaluate_model(model_name, time, fparams, current_grad_u, current_phi, current_grad_phi, previous_grad_u, previous_phi, previous_grad_phi, SDVS, current_ADD_DOF, previous_ADD_DOF, current_ADD_grad_DOF, previous_ADD_grad_DOF, answers):
    """
    Test the evaluation of the model

    :param str model_name: The name of the model to evaluate
    :param np.ndarray time: The time
    :param np.ndarray fparams: The model parameters
    :param np.ndarray current_grad_u: The current gradient of the macro displacement
    :param np.ndarray current_phi: The current micro displacement
    :param np.ndarray current_grad_phi: The current gradient of the micro displacement
    :param np.ndarray previous_grad_u: The previous gradient of the macro displacement
    :param np.ndarray previous_phi: The previous micro displacement
    :param np.ndarray previous_grad_phi: The previous gradient of the micro displacement
    :param np.ndarray SDVS: The solution dependent state variables
    :param np.ndarray current_ADD_DOF: The current values additional degrees of freedom required for the model
    :param np.ndarray previous_ADD_DOF: The previous values of the additional degrees of freedom required for the model
    :param np.ndarray current_ADD_grad_DOF: The current values of the gradients of the additional degrees of freedom required for the model
    :param np.ndarray previous_ADD_grad_DOF: The previous values of the gradients of the additional degrees of freedom required for the model
    :param dict answers: The answer dictionary
    """

    keys = ['errorCode', 'PK2', 'SIGMA', 'M', 'SDVS',\
            'DPK2Dgrad_u', 'DPK2Dphi', 'DPK2Dgrad_phi',\
            'DSIGMADgrad_u', 'DSIGMADphi', 'DSIGMADgrad_phi',\
            'DMDgrad_u', 'DMDphi', 'DMDgrad_phi',\
            'ADD_TERMS', 'ADD_JACOBIANS', 'output_message']

    values = micromorphic.evaluate_model(model_name,\
                                         time, fparams,\
                                         current_grad_u, current_phi, current_grad_phi,\
                                         previous_grad_u, previous_phi, previous_grad_phi,\
                                         SDVS,\
                                         current_ADD_DOF, current_ADD_grad_DOF,\
                                         previous_ADD_DOF, previous_ADD_grad_DOF)

    assert len(keys) == len(values)

    results = dict(zip(keys, values))

    if (results['errorCode'] != 0):
        print(results['output_message'].decode('utf-8'))
        assert results['errorCode'] == 0

    for key in answers:

        assert np.allclose(answers[key], results[key])
