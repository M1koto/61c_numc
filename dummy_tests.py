from utils import *
from unittest import TestCase


def test_small_add():
    # TODO: YOUR CODE HERE
    dp_mat1, nc_mat1 = rand_dp_nc_matrix(2, 2, seed=0)
    # print('mat1')
    dp_mat2, nc_mat2 = rand_dp_nc_matrix(2, 2, seed=1)
    # print('mat2')
    is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
    self.assertTrue(is_correct)
    print_speedup(speed_up)


test_small_add()
