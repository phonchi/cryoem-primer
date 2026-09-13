"""Numerical checks for the recovered step/Snorlax/page/coins examples."""
import numpy as np
from scipy.signal import convolve2d

def test_original_step_averages(restored_chapter):
    c=restored_chapter('02_filter_segment.py')
    assert len(c['noisy_signal'])==100
    assert len(c['smooth_signal3'])==98 and len(c['smooth_signal11'])==90
    np.testing.assert_allclose(c['smooth_signal3'],np.convolve(c['noisy_signal'],np.ones(3)/3,'valid'))
    np.testing.assert_allclose(c['smooth_signal11'],np.convolve(c['noisy_signal'],np.ones(11)/11,'valid'))

def test_original_bright_square_and_hog(restored_chapter):
    c=restored_chapter('02_filter_segment.py')
    a=c['bright_square'];assert a.shape==(7,7)
    out=convolve2d(a,np.ones((3,3))/9,mode='same')
    assert np.isclose(out[3,3],1) and np.isclose(out[2,2],4/9)
    # Original 128x128 Snorlax: 16x16 single-cell blocks, eight bins.
    assert c['fd'].size==16*16*8

def test_watershed_original_labels_are_not_coin_instances(restored_chapter):
    c=restored_chapter('02_filter_segment.py')
    np.testing.assert_array_equal(np.unique(c['segmentation_coins']),[1,2])
    areas=[p.area for p in c['properties']]
    assert areas==[77442,38910]
    assert sum(areas)==c['coins'].size
