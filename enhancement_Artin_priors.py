import nibabel as nib
# import Image
import matplotlib.pylab as plt
from PIL import ImageEnhance , Image , ImageFilter
import numpy as np
import nifti
import os
from glob import glob

# EnhMethod = 'Sharpness' #'Contrast' # Sharpness   +'_int16DivideMultipliedBy7'  # Contrast
def enhancing(im , scaleEnhance):
    # im = im.astype(dtype=np.int16)
    im = Image.fromarray(im)
    im = im.convert('L')

    if scaleEnhance[0] != 1:
        im2 = ImageEnhance.Sharpness(im)
        im = im2.enhance(scaleEnhance[0])

    if scaleEnhance[1] != 1:
        im2 = ImageEnhance.Contrast(im)
        im = im2.enhance(scaleEnhance[1])

    return im


Directory = '/media/data1/artin/thomas/priors/'
# subDirs = os.listdir(Directory)

subDirsFull = glob(Directory+'/*/')
# subDirs = subDirsFull[0]
scaleEnhance = [[1,2],[1,3],[4,1],[6,1],[4,2],[4,3]]

for subDirs in subDirsFull:print subDirs
    print subDirs
    for name in ['WMnMPRAGE_bias_corr' , 'WMnMPRAGEdeformed']:

        im = nib.load(subDirs + name + '.nii.gz')
        imD = im.get_data()
        MaxValue = imD.max()
        # imD = imD.astype(float)*256/imD.max()
        imD = imD*256/imD.max()

        sz = imD.shape

        for s in scaleEnhance:
            print s
            imEnhanced = np.zeros(sz)
            for i in range(sz[2]):
                imEnhanced[:,:,i] = enhancing(imD[:,:,i] , s)

            imEnhanced = imEnhanced/256*MaxValue
            imEnhanced_nifti = nib.Nifti1Image(imEnhanced , im.affine , im.header)

            string = subDirs + name + '_' + 'Sharpness_' + str(s[0]) + '_Contrast_' + str(s[1]) + '.nii.gz'
            print string

            nib.save(imEnhanced_nifti,string)

            # sliceNum = 40
            # fig , axes = plt.subplots(1,2 , figsize=(10,5))
            # axes[0].imshow(imEnhanced_nifti.get_data()[:,:,sliceNum],cmap='gray',aspect='auto')
            # axes[1].imshow(imD[:,:,sliceNum],cmap='gray',aspect='auto')
            # plt.show()

Directory = '/media/data1/artin/thomas/'
scaleEnhance = [[1,2],[1,3],[4,1],[6,1],[4,2],[4,3]]
for name in ['templ_93x187x68']:

    im = nib.load(Directory + name + '.nii.gz')
    imD = im.get_data()
    MaxValue = imD.max()
    # imD = imD.astype(float)*256/imD.max()
    imD = imD*256/imD.max()

    sz = imD.shape

    for s in scaleEnhance:
        print s
        imEnhanced = np.zeros(sz)
        for i in range(sz[2]):
            imEnhanced[:,:,i] = enhancing(imD[:,:,i] , s)

        imEnhanced = imEnhanced/256*MaxValue
        imEnhanced_nifti = nib.Nifti1Image(imEnhanced , im.affine , im.header)

        string = Directory + name + '_' + 'Sharpness_' + str(s[0]) + '_Contrast_' + str(s[1]) + '.nii.gz'
        print string

        nib.save(imEnhanced_nifti,string)
