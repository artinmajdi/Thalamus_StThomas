import os
import sys
import require
from parallel import command


def check_run(fname, func, *args, **kwargs):
    """
    Checks if fname exists first before execuring func.
    """
    if os.path.exists(fname):
        print 'Skipped, using %s' % fname
        return fname
    else:
        return func(*args, **kwargs)


"""
segment.py
"""
def check_warps(warp_path):
    """
    Checks if the necessary ANTS warps exist.
    """
    warp_file = os.path.join(f'{warp_path}InverseWarp.nii.gz')
    affine_file = os.path.join(f'{warp_path}Affine.txt')
    return bool(os.path.exists(warp_file) and os.path.exists(affine_file))


def sanitize_input(input_image, output_image, command=os.system):
    """
    Standardizes the input to neurological coordinates and can flip to segment
    right thalamus.
    """
    if sys.platform in ['linux2', 'darwin']:
        command(f'fslreorient2std {input_image} {output_image}')
    return output_image


def flip_lr(input_image, output_image, command=os.system):
    if sys.platform in ['linux2', 'darwin']:
        command(f'fslswapdim {input_image} -x y z {output_image}')
    # else:
    #     command('fsl5.0-fslswapdim %s -x y z %s' % (input_image, output_image))
    return output_image


"""
warp_to_all_via_crop
"""
def copy_header(reference, target, output, switches='1 1 1', echo=False):
    """
    Copies NIFTI header information from reference to target resulting in output.
    Usage:  CopyImageHeaderInformation refimage.ext imagetocopyrefimageinfoto.ext imageout.ext   boolcopydirection  boolcopyorigin boolcopyspacing  {bool-Image2-IsTensor}
    """
    cmd = f'CopyImageHeaderInformation {reference} {target} {output} {switches}'
    command(cmd, echo=echo)
    return output, cmd


"""
warp_to_all
"""
def ants_compose_a_to_b(a_transform_prefix, b_path, b_transform_prefix, output, **exec_options):
    """
    Compose a to b via an intermediate template space
    """
    a_affine = f'{a_transform_prefix}Affine.txt'
    a_warp = f'{a_transform_prefix}Warp.nii.gz'

    b_affine = f'-i {b_transform_prefix}Affine.txt'
    b_warp = f'{b_transform_prefix}InverseWarp.nii.gz'
    cmd = f'ComposeMultiTransform 3 {output} {b_affine} {b_warp} {a_warp} {a_affine} -R {b_path}'
    command(cmd, **exec_options)
    return output, cmd


def ants_apply_only_warp(template, input_image, input_warp, output_image, switches='', **exec_options):
    cmd = f'WarpImageMultiTransform 3 {input_image} {output_image} {input_warp} -R {template} {switches}'
    command(cmd, **exec_options)
    return output_image, cmd


def sanitize_label_image(input_image, output_image, **exec_options):
    # Slicer puts data in LR convention for some reason
    cmd = f'fslswapdim {input_image} LR PA IS {output_image}; '
    # Sometimes values are > 1
    cmd += f'fslmaths {output_image} -bin {output_image}'
    command(cmd, **exec_options)
    return output_image, cmd


"""
cv.py
"""
def create_atlas(label, path, subjects, target, output_atlas, echo=False):
    # Create 4D atlas for a label previously registered to a target subject
    label_paths = [
        os.path.join(path, subj, target, f'{label}.nii.gz')
        for subj in subjects
    ]
    cmd = f"fslmerge -t {output_atlas} {' '.join(label_paths)}"
    # if echo:
    #     print cmd
    # else:
    #     os.system(cmd)
    command(cmd, echo=echo)
    return output_atlas, cmd



"""
cv_registration_method
"""
def crop_by_mask(input_image, output_image, mask, label=1, padding=0):
    return f'ExtractRegionFromImageByMask 3 {input_image} {output_image} {mask} {label} {padding}'


def ants_label_fusions(output_prefix, labels, images=None):
    """
    Returns commands for various ANTS label fusion schemes.
    For correlation voting, the last element of images should be the target image to compare priors against.
    """
    l = ' '.join(labels)
    # Maximum
    output = f'{output_prefix}_maximum.nii.gz'
    cmd = f'AverageImages 3 {output} 0 {l};'
    cmd += f'ThresholdImage 3 {output} {output} 0.01 1000'
    cmds = [cmd]
    outputs = [output]
    # Majority
    output = f'{output_prefix}_majority.nii.gz'
    cmd = f'ImageMath 3 {output} MajorityVoting {l}'
    cmds.append(cmd)
    outputs.append(output)
    # STAPLE
    output = f'{output_prefix}_staple.nii.gz'
    output_probability = f'{output_prefix}_staple0001.nii.gz'
    confidence = 0.5
    cmd = f'ImageMath 3 {output} STAPLE {confidence} {l};'
    # Threshold at 0.5 even though this is known to be loose (Cardoso STEPS 2013), that's okay for this purpose
    cmd += f'ThresholdImage 3 {output_probability} {output} 0.5 1000'
    cmds.append(cmd)
    outputs.append(output)
    if images is not None:
        # Correlation Vote
        output = f'{output_prefix}_correlation.nii.gz'
        template = images.pop()
        cmd = f"ImageMath 3 {output} CorrelationVoting {template} {' '.join(images)} {l}"
        cmds.append(cmd)
        outputs.append(output)
    return outputs, cmds


def crop_prior_using_transform(output, crop, mask, padding, prior, affine, prior_padding=0, includes=[], output_mask=None):
    """
    Takes a cropped box with the mask and padding amount that defined it and brings it to another space via an (inverse) affine transform.
    - includes is a list of label masks in the prior space that must be included in the cropped_prior, assumes label masks are 0/1
    """
    import random

    if output_mask is None:
        output_mask = output
    # Make all 1s
    # CreateImage imageDimension referenceImage outputImage constant [random?]
    # Need to use a temporary file here beause uncrop.py overwrites with the empty output image at the first step
    # TODO use python temp file
    ones = f'{output}_DELETEME_{str(random.random())[2:]}.nii.gz'
    cmds = [
        f'CreateImage 3 {crop} {ones} 1',
        f"{os.path.join(sys.path[0], 'uncrop.py')} {ones} {output_mask} {mask} {padding}",
        f'WarpImageMultiTransform 3 {output_mask} {output_mask} -R {prior} --use-NN -i {affine}',
    ]
    # Incorporate must-include label masks, assumes label masks are 0/1, otherwise overadd will be problematic
    cmds.extend(
        f'ImageMath 3 {output_mask} overadd {output_mask} {include}'
        for include in includes
    )
    cmds.extend(
        (
            crop_by_mask(prior, output, output_mask, padding=prior_padding),
            f'rm {ones}',
        )
    )
    return '; '.join(cmds)


"""
cv_picsl
"""
def label_fusion_steps(input_image, image_atlas, label_atlas, output_label, sigma, X, mrf=0., echo=False):
    # Perform steps label fusion.  Parameter naming comes from Cardoso et al. 2013
    # verbose and only consider non-consensus voxels.
    cmd = f'seg_LabFusion -v -unc -in {label_atlas} -STEPS {sigma} {X} {input_image} {image_atlas} -out {output_label}'
    if 0 < mrf <= 5:
        cmd += ' -MRF_beta %g' % mrf
    # command(cmd, echo=echo)
    if not echo:
        command(cmd)
    return output_label, cmd


def label_fusion_picsl(input_image, atlas_images, atlas_labels, output_label, rp=[2, 2, 2], rs=[3, 3, 3], alpha=0.1, beta=2, **exec_options):
    """
    H Wang. Multi-Atlas Sementation with Joint Label Fusion. 2013.
    Joint Label Fusion:
    usage:
     jointfusion dim mod [options] output_image
    required options:
      dim                             Image dimension (2 or 3)
      mod                             Number of modalities or features
      -g atlas1_mod1.nii atlas1_mod2.nii ...atlasN_mod1.nii atlasN_mod2.nii ...
                                      Warped atlas images
      -tg target_mod1.nii ... target_modN.nii
                                      Target image(s)
      -l label1.nii ... labelN.nii    Warped atlas segmentation
      -m <method> [parameters]        Select voting method. Options: Joint (Joint Label Fusion)
                                      May be followed by optional parameters in brackets, e.g., -m Joint[0.1,2].
                                      See below for parameters
    other options:
      -rp radius                      Patch radius for similarity measures, scalar or vector (AxBxC)
                                      Default: 2x2x2
      -rs radius                      Local search radius.
                                      Default: 3x3x3
      -x label image.nii              Specify an exclusion region for the given label.
      -p filenamePattern              Save the posterior maps (probability that each voxel belongs to each label) as images.
                                      The number of images saved equals the number of labels.
                                      The filename pattern must be in C printf format, e.g. posterior%04d.nii.gz
    Parameters for -m Joint option:
      alpha                           Regularization term added to matrix Mx for inverse
                                      Default: 0.1
      beta                            Exponent for mapping intensity difference to joint error
                                      Default: 2
    """
    dim = 3
    mod = 1
    g = ' '.join(atlas_images)
    tg = input_image
    l = ' '.join(atlas_labels)
    m = 'Joint[%g,%g]' % (alpha, beta)
    rp = '%dx%dx%d' % tuple(rp)
    rs = '%dx%dx%d' % tuple(rs)
    cmd = f'jointfusion {dim} {mod} -g {g} -tg {tg} -l {l} -m {m} -rp {rp} -rs {rs} {output_label}'
    command(cmd, **exec_options)
    # if not echo:
    #     command(cmd)
    return output_label, cmd


def label_fusion_picsl_ants(input_image, atlas_images, atlas_labels, output_label, rp=[2, 2, 2], rs=[3, 3, 3], alpha=0.1, beta=2, mask='', **exec_options):
    """
    COMMAND:
         antsJointFusion
              antsJointFusion is an image fusion algorithm developed by Hongzhi Wang and Paul
              Yushkevich which won segmentation challenges at MICCAI 2012 and MICCAI 2013. The
              original label fusion framework was extended to accommodate intensities by Brian
              Avants. This implementation is based on Paul's original ITK-style implementation
              and Brian's ANTsR implementation. References include 1) H. Wang, J. W. Suh, S.
              Das, J. Pluta, C. Craige, P. Yushkevich, Multi-atlas segmentation with joint
              label fusion IEEE Trans. on Pattern Analysis and Machine Intelligence, 35(3),
              611-623, 2013. and 2) H. Wang and P. A. Yushkevich, Multi-atlas segmentation
              with joint label fusion and corrective learning--an open source implementation,
              Front. Neuroinform., 2013.
    OPTIONS:
         -d, --image-dimensionality 2/3/4
              This option forces the image to be treated as a specified-dimensional image. If
              not specified, the program tries to infer the dimensionality from the input
              image.
         -t, --target-image targetImage
                            [targetImageModality0,targetImageModality1,...,targetImageModalityN]
              The target image (or multimodal target images) assumed to be aligned to a common
              image domain.
         -g, --atlas-image atlasImage
                           [atlasImageModality0,atlasImageModality1,...,atlasImageModalityN]
              The atlas image (or multimodal atlas images) assumed to be aligned to a common
              image domain.
         -l, --atlas-segmentation atlasSegmentation
              The atlas segmentation images. For performing label fusion the number of
              specified segmentations should be identical to the number of atlas image sets.
         -a, --alpha 0.1
              Regularization term added to matrix Mx for calculating the inverse. Default =
              0.1
         -b, --beta 2.0
              Exponent for mapping intensity difference to the joint error. Default = 2.0
         -r, --retain-label-posterior-images (0)/1
              Retain label posterior probability images. Requires atlas segmentations to be
              specified. Default = false
         -f, --retain-atlas-voting-images (0)/1
              Retain atlas voting images. Default = false
         -c, --constrain-nonnegative (0)/1
              Constrain solution to non-negative weights.
         -p, --patch-radius 2
                            2x2x2
              Patch radius for similarity measures. Default = 2x2x2
         -m, --patch-metric (PC)/MSQ
              Metric to be used in determining the most similar neighborhood patch. Options
              include Pearson's correlation (PC) and mean squares (MSQ). Default = PC (Pearson
              correlation).
         -s, --search-radius 3
                             3x3x3
                             searchRadiusMap.nii.gz
              Search radius for similarity measures. Default = 3x3x3. One can also specify an
              image where the value at the voxel specifies the isotropic search radius at that
              voxel.
         -e, --exclusion-image label[exclusionImage]
              Specify an exclusion region for the given label.
         -x, --mask-image maskImageFilename
              If a mask image is specified, fusion is only performed in the mask region.
         -o, --output labelFusionImage
                      intensityFusionImageFileNameFormat
                      [labelFusionImage,intensityFusionImageFileNameFormat,<labelPosteriorProbabilityImageFileNameFormat>,<atlasVotingWeightImageFileNameFormat>]
              The output is the intensity and/or label fusion image. Additional optional
              outputs include the label posterior probability images and the atlas voting
              weight images.
         --version
              Get version information.
         -v, --verbose (0)/1
              Verbose output.
         -h
              Print the help menu (short version).
         --help
              Print the help menu.
    """
    dim = 3
    g = ' '.join(atlas_images)
    tg = input_image
    l = ' '.join(atlas_labels)
    rp = 'x'.join(['%d' % el for el in rp])
    rs = 'x'.join(['%d' % el for el in rs])
    if mask:
        mask = f'-x {mask}'
    cmd = 'antsJointFusion -d %s -g %s -t %s -l %s -a %g -b %g -p %s -s %s %s -o %s' % (dim, g, tg, l, alpha, beta, rp, rs, mask, output_label)
    command(cmd, **exec_options)
    return output_label, cmd
