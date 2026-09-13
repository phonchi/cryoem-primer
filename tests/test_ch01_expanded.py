"""Checks against the original chapter-1 examples and geometry parameters."""
import numpy as np
from skimage.transform import EuclideanTransform,SimilarityTransform,AffineTransform,rotate,warp

def test_original_assets_channels_and_crop(restored_chapter):
    c=restored_chapter('01_image_basics.py')
    assert c['venusaurf'].shape==(80,80,4)
    np.testing.assert_array_equal(c['venusaurf'][40,75],[0,0,0,255])
    np.testing.assert_array_equal(c['venusaurf'][40,76],[0,255,0,255])
    np.testing.assert_array_equal(c['top_left'],c['charizardy'][:100,:100])
    assert c['image_rescaled'].shape==c['image_resized'].shape==c['image_downscaled'].shape==(128,128)

def test_original_control_points_rectify_text(restored_chapter):
    c=restored_chapter('01_image_basics.py')
    np.testing.assert_allclose(c['tform3'](c['src']),c['dst'],atol=1e-10)
    assert c['warped'].shape==(50,300)
    np.testing.assert_allclose(c['warped'],warp(c['text'],c['tform3'],output_shape=(50,300)))

def test_original_rotation_similarity_and_shear_matrices():
    angle=np.pi/12;R=np.array([[np.cos(angle),-np.sin(angle),30],[np.sin(angle),np.cos(angle),-20],[0,0,1]])
    np.testing.assert_allclose(EuclideanTransform(rotation=angle,translation=(30,-20)).params,R)
    S=R.copy();S[:2,:2]*=.5
    np.testing.assert_allclose(SimilarityTransform(scale=.5,rotation=angle,translation=(30,-20)).params,S)
    np.testing.assert_allclose(AffineTransform(shear=np.pi/6).params,[[1,-1/np.sqrt(3),0],[0,1,0],[0,0,1]])

def test_center_rotation_positive_direction():
    a=np.zeros((33,33));a[16,24]=1
    result=rotate(a,90,order=0,preserve_range=True)
    np.testing.assert_array_equal(np.argwhere(result>.5),[[8,16]])
