import nibabel as nib
# import Image
import matplotlib.pylab as plt
from PIL import ImageEnhance , Image , ImageFilter
import numpy as np
import nifti
import os

EnhMethod = 'Sharpness' #'Contrast' # Sharpness   +'_int16DivideMultipliedBy7'  # Contrast
def enhancing(im , scaleEnhance):
    # im = im.astype(dtype=np.int16)
    im = Image.fromarray(im)
    im = im.convert('L')
    # print 'EnhMethod:' + im.EnhMethod +'end'

    # im2 = ImageEnhance.Contrast(im)
    if EnhMethod == 'Sharpness':
        im2 = ImageEnhance.Sharpness(im)
    elif EnhMethod == 'Contrast':
        im2 = ImageEnhance.Contrast(im)

    im = im2.enhance(scaleEnhance)

    return im


for InputIm in ['template' , 'testFiles']:

    print InputIm

    if InputIm == 'template':
        Directory = '/media/data1/artin/thomas/'
        name = 'templ_93x187x68.nii.gz'

    elif InputIm == 'testFiles':
        Directory = '/media/data1/artin/vimp2_TestA/'
        name = 'WMnMPRAGEdeformed.nii.gz'  # WMnMPRAGEdeformed origtemplate  wmnpre wmnpre_Contrast_int16DivideMultipliedBy70_5.nii.gz

    im = nib.load(Directory+name)
    imD = im.get_data()
    MaxValue = imD.max()
    # imD = imD.astype(float)*256/imD.max()
    imD = imD*256/imD.max()

    scaleEnhance = [0.5 , 1.5 , 2.1 , 4.1 , 6.1 , 8.3 , 10.2]
    sz = imD.shape

    for EnhMethod in ['Sharpness' , 'Contrast']:
        print EnhMethod

        for s in range(len(scaleEnhance)):

            Divider = 1
            imEnhanced = np.zeros(sz)
            for i in range(sz[2]):
                imEnhanced[:,:,i] = enhancing(imD[:,:,i]/Divider , scaleEnhance[s])
                imEnhanced[:,:,i] = Divider*imEnhanced[:,:,i]

            imEnhanced = imEnhanced/256*MaxValue
            # imEnhanced = np.int16(imEnhanced)

            imEnhanced_nifti = nib.Nifti1Image(imEnhanced , im.affine , im.header)

            if InputIm == 'testFiles':

                if EnhMethod == 'Sharpness':
                    Directory2 = Directory + 'Test' + str(s+1) + '/'
                elif EnhMethod == 'Contrast':
                    Directory2 = Directory + 'Test' + str(s+1+len(scaleEnhance)) + '/'

                try:
                    os.stat(Directory2)
                except:
                    os.mkdir(Directory2)

                try:
                    os.stat(Directory2 + 'temp/')
                except:
                    os.mkdir(Directory2 + 'temp/')

                string = Directory2 + name

            elif InputIm == 'template':
                string = Directory + name

            string = string.split('.nii.gz')[0] + '_' + EnhMethod+str(scaleEnhance[s]) + '.nii.gz'
            print string

            nib.save(imEnhanced_nifti,string)


            # sliceNum = 40
            # fig , axes = plt.subplots(1,2 , figsize=(10,5))
            # axes[0].imshow(imEnhanced_nifti.get_data()[:,:,sliceNum],cmap='gray',aspect='auto')
            # axes[1].imshow(imD[:,:,sliceNum],cmap='gray',aspect='auto')
            # plt.show()
