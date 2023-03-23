#!/usr/bin/env python
"""
Segment subject for selected thalamic nuclei using whole-brain registration via a template and PICSL label fusion.
"""

import os
import sys
import argparse
import tempfile
import time
import libraries.parallel as parallel
from shutil import rmtree
from functools import partial
from datetime import timedelta
from libraries.imgtools_v0 import check_run, check_warps, sanitize_input, flip_lr, label_fusion_picsl_ants, label_fusion_picsl, ants_compose_a_to_b, ants_apply_only_warp, ants_WarpImageMultiTransform, crop_by_mask, label_fusion_majority
from libraries.ants_nonlinear_v0 import ants_nonlinear_registration, bias_correct, ants_linear_registration
from THOMAS_constants_v0 import image_name, orig_template, template_93, mask_93, this_path, prior_path, subjects, roi, roi_choices, optimal
import nibabel
import numpy as np
import nibabel as nib
from nilearn import image


parser = argparse.ArgumentParser(description='Shortened Template and THalamus for Optimal Multi Atlas Segmentation (ST THOMAS) for a given WMnMPRAGE image using Majority Voting or Joint Label Fusion and the Tourdias atlas.')
parser.add_argument('input_image', help='input WMnMPRAGE NiFTI image, may need to be in LR PA IS format')
#parser.add_argument('output_path', help='the output file for single ROI or directory for multiple ROIs')
parser.add_argument(
    'roi_names',
    metavar='roi_names',
    choices=roi_choices,
    nargs='+',
    help=f"a space separated list of one or more ROIs.  Valid targets are: {', '.join(roi_choices)}",
)
parser.add_argument('-a', '--algorithm', type=str, required=True, help='version of THOMAS: v0 or v1 or v2')
parser.add_argument('-w', '--warp', metavar='path', help='looks for {path}InverseWarp.nii.gz and {path}Affine.txt instead of basing it off input_image.')
parser.add_argument('-F', '--forcereg', action='store_true', help='force ANTS registration to WMnMPRAGE mean brain template. The --warp argument can be then used to specify the output path.')
parser.add_argument('-p', '--processes', nargs='?', default=None, const=None, type=int, help='number of parallel processes to use.  If unspecified, automatically set to number of CPUs.')
parser.add_argument('-v', '--verbose', action='store_true', help='verbose mode')
parser.add_argument('-d', '--debug', action='store_true', help='debug mode, interactive prompts')
parser.add_argument('-R', '--right', action='store_true', help='segment right thalamus')
parser.add_argument('-M', '--majorityvoting', action='store_true', help='use majority voting for joint fusion')
parser.add_argument('--jointfusion', action='store_true', help='use older jointfusion instead of antsJointFusion')
parser.add_argument('--tempdir', help='temporary directory to store registered atlases.  This will not be deleted as usual.')
parser.add_argument('--mask', help='custom mask if 93x187x68 mask size is not wanted')
parser.add_argument('--template', help='custom template if 93x187x68 size is not wanted')
parser.add_argument('--Downsample', help='Downsample the input (linear part) for increasing the speed')
parser.add_argument('--BiosCorrectName', help='WMnMPRAGE_bias_corr')
parser.add_argument('--MaskInput', help='MaskInput')
parser.add_argument('--output_path', help='specify a different output_path for the output file for single ROI or directory for multiple ROIs')
# TODO handle single roi, single output file case
# TODO fix verbose and debug
# TODO go back to shell=False for command to suppress output and then fix sanitize labels 4

# print '---------------------------------------'
# if args.BiosCorrectName is not None:
#     BiosCorrectName = args.BiosCorrectName
# else:
#     BiosCorrectName = 'WMnMPRAGE_bias_corr'
#
# print BiosCorrectName
# print '---------------------------------------'


def warp_atlas_subject(subject, path, labels, input_image , input_transform_prefix, output_path, BiosCorrectNameV, exec_options={}):
    """
    Warp a training set subject's labels to input_image.
    """

    if args.BiosCorrectName is not None:
        image_name = BiosCorrectNameV +'.nii.gz'

    a_transform_prefix = os.path.join(path, subject + '/WMnMPRAGE')
    output_path = os.path.join(output_path, subject)
    try:
        os.mkdir(output_path)
    except OSError:
        # Exists
        pass
    combined_warp = os.path.join(output_path, 'Warp.nii.gz')
    if not os.path.exists(combined_warp):
        check_run(
            combined_warp,
            ants_compose_a_to_b,
            a_transform_prefix,
            b_path=input_image,
            b_transform_prefix=input_transform_prefix,
            output=combined_warp,
            **exec_options
        )
    output_labels = {}
    # OPT parallelize, or merge parallelism with subject level
    for label in labels:
        label_fname = os.path.join(path, subject, 'sanitized_rois', label + '.nii.gz')
        warped_label = os.path.join(output_path, label + '.nii.gz')
        switches = '--use-NN'
        check_run(
            warped_label,
            ants_apply_only_warp,
            template=input_image,
            input_image=label_fname,
            input_warp=combined_warp,
            output_image=warped_label,
            switches=switches,
            **exec_options
        )
        output_labels[label] = warped_label
    # Warp anatomical WMnMPRAGE_bias_corr too
    # TODO merge this into previous for loop to be DRY?
    output_labels[BiosCorrectNameV] = output_image = os.path.join(output_path, image_name)
    if not os.path.exists(output_labels[BiosCorrectNameV]):
        print output_labels[BiosCorrectNameV]
        check_run(
            output_image,
            ants_apply_only_warp,
            template=input_image,
            input_image=os.path.join(path, subject, image_name),
            input_warp=combined_warp,
            output_image=output_image,
            switches='--use-BSpline',
            **exec_options
        )
    return output_labels


def conservative_mask(input_masks, output_path, dilation=0, fill=False):
    """
    Estimates a conservative maximum mask given a list of input masks.
    - for dilation > 0 and fill=True, each side is padded by dilation instead
    - fill will fill the bounding box of the mask producing a cube
    """
    # Maximum label fusion
    # Taken from cv_registration_method.ants_label_fusions
    # cmd = 'AverageImages 3 %s 0 %s' % (output_path, ' '.join(input_masks))
    # cmd = 'ThresholdImage 3 %s %s 0.01 1000' % (output_path, output_path)
    cmd = f"c3d {' '.join(input_masks)} -accum -max -endaccum -binarize -o {output_path}"
    # cmd = 'fslmaths %s -bin %s' % (' -add '.join(input_masks), output_path)
    if sys.platform in ['linux2', 'darwin']:
        parallel_command(cmd)
    if fill:
        # get bounding box
        bbox = map(int, os.popen(f'fslstats {output_path} -w').read().strip().split())
        padding = (-dilation, 2 * dilation) * 3  # min index, size change for 3 spatial dimensions
        if dilation > 0:
            # edit bounding box
            for i, inc in enumerate(padding):  # ignore time dimensions
                bbox[i] += inc
        roi = ' '.join(map(str, bbox))
        # fill bounding box
        cmd = f'fslmaths {output_path} -add 1 -bin -roi {roi} {output_path}'
    elif dilation > 0:
        kernel = '%dx%dx%dvox' % (dilation, dilation, dilation)
        cmd = f'c3d {output_path} -dilate 1 {kernel} -o {output_path}'
    if sys.platform in ['linux2', 'darwin']:
        parallel_command(cmd)
    return output_path


def get_bounding_box(A):
    B = np.argwhere(A)
    start, stop = B.min(0), B.max(0) + 1
    return zip(start, stop)


def split_roi(roi, axis, split_axis):
    """
    Pages through the roi along axis and cuts each slice in half in the split_axis dimension.
    axis=None cuts the 3D bounding box in half.
    """
    first = np.zeros_like(roi)
    second = np.zeros_like(roi)
    def split_halves(roi, first, second, sl, axis):
        N = len(roi.shape)
        idx = [slice(sl, sl+1) if el is axis else slice(None) for el in xrange(N)]
        try:
            box = get_bounding_box(roi[idx])
        except ValueError:
            # No ROI in this slice
            return
        try:
            box[axis] = tuple(el+sl for el in box[axis])
        except TypeError:
            # Occurs for axis=None case
            pass
        first_idx = [slice(a, a + (b-a)/2) if i is split_axis else slice(a, b) for i, (a, b) in enumerate(box)]
        second_idx = [slice(a + (b-a)/2, b) if i is split_axis else slice(a, b) for i, (a, b) in enumerate(box)]
        # Try pasting in the original ROI to the half boxes
        try:
            first[first_idx] = roi[first_idx]
        except ValueError:
            # exception if half box is 0 along one dimension
            pass
        try:
            second[second_idx] = roi[second_idx]
        except ValueError:
            pass
    if axis is None:
        split_halves(roi, first, second, 0, axis)
    else:
        for sl in xrange(roi.shape[axis]):
            split_halves(roi, first, second, sl, axis)
    return first, second

def main(args, temp_path, pool):


    print orig_template

    # orig_template = orig_template.split('orig_template.nii.gz') + 'origtemplate_Contrast0.5.nii.gz'
# """
    input_image = orig_input_image = args.input_image



    # assigning default value of mask
    mask = mask_93

    #setting up output path
    if args.output_path:
        output_path = args.output_path
    else:
        output_path = os.path.dirname(orig_input_image)

    #setting up the ROIs
    if roi['param_all'] in args.roi_names:
        labels = list(roi['label_names'])
    else:
        roi_dict = dict(zip(roi['param_names'], roi['label_names']))
        labels = [roi_dict[el] for el in args.roi_names]

    #setting up the template
    if args.algorithm == "v2":
        if args.template is not None and args.mask is not None:
            template = args.template
            mask = args.mask
            print "Custom template and mask"
        elif args.template is not None and args.mask is None:
            sys.exit("!!!!!!! Both template and mask need to be specified simultaneously and they need to be of the same size !!!!!!!")
        elif args.template is None and args.mask is not None:
            sys.exit("!!!!!!! Both template and mask need to be specified simultaneously and they need to be of the same size !!!!!!!")
        else:
            template = template_93
            mask = mask_93
            print "Algorithm is v2"
    elif args.algorithm == "v1":
        sys.exit("!!!!!!! v1 algorithm not yet implemented !!!!!!!")
    elif args.algorithm == "v0":
        template = orig_template
        print "Template is origtemplate.nii.gz"
    else:
        sys.exit("!!!!!!! Algorithm incorrectly specified !!!!!!!")

    # print 'Template being used is'
    # print os.path.abspath(template)

    print('-------------*********output_path************------------')
    print(output_path)

    # TODO prevent both jointfusion and majority voting being set
	# if args.jointfusion is None:
		# print "args.jointfusion has been set (value is %s)" % args.jointfusion
		# if args.majorityvoting is None:
			# print "args.majorityvoting has been set (value is %s)" % args.majorityvoting
			# sys.exit("!!!!!!! Only one label fusion can be selected at any time (default is antsJointFusion) !!!!!!!")

    if args.warp:
        warp_path = args.warp
    else:
        # TODO remove this as the default behavior, instead do ANTS?
        head, tail = os.path.split(input_image)
        tail = tail.replace('.nii', '').replace('.gz', '') #split('.', 1)[0]
        warp_path = os.path.join(temp_path, tail)

    t = time.time()
    t1 = time.time()
    if args.algorithm == "v2":
        # Crop the input
        # Affine registering input to template

        print "1.   Linear Registration of input and template \n"

        print 'orig_template: %s' % orig_template
        print 'orig_input_image: %s' % orig_input_image

        print '---------------------------------------'
        if args.BiosCorrectName is not None:
            BiosCorrectName = args.BiosCorrectName
        else:
            BiosCorrectName = 'WMnMPRAGE_bias_corr'

        print BiosCorrectName
        print '---------------------------------------'

        if args.MaskInput is None:
            if args.Downsample is not None:

                t1 = time.time()
                SamplingRate = int(args.Downsample)

                # SamplingRate = 4
                print '--------------------------DownSample_orig_template----------------------------------------'
                imm = nib.load(orig_template)

                NewAffine = imm.affine.copy()
                for i in range(0,3):
                    NewAffine[i,i] = imm.affine[i,i]*SamplingRate

                NewSize = [imm.shape[0]/SamplingRate , imm.shape[1]/SamplingRate , imm.shape[2]/SamplingRate]

                img_DS = image.resample_img(imm, target_affine=NewAffine,target_shape=NewSize )

                orig_template_DS = orig_template.split('.nii.gz')[0]+'_DS.nii.gz'
                nib.save(img_DS,orig_template_DS)
                # img_DS_US = image.resample_img(img_DS, target_affine=img.affine,target_shape=img.shape )

                print '-------------------------DownSample_orig_input_image-------------------------------------'
                print orig_input_image
                imm = nib.load(orig_input_image)
                NewAffine = imm.affine.copy()
                for i in range(0,3):
                    NewAffine[i,i] = imm.affine[i,i]*SamplingRate

                NewSize = [imm.shape[0]/SamplingRate , imm.shape[1]/SamplingRate , imm.shape[2]/SamplingRate]

                img_DS = image.resample_img(imm, target_affine=NewAffine,target_shape=NewSize )

                orig_input_image_DS = orig_input_image.split('.nii.gz')[0]+'_DS.nii.gz'
                nib.save(img_DS,orig_input_image_DS)


                print '0.   -----Down Sampling Time Elapsed: %s  \n \n' % timedelta(seconds=time.time()-t1)
                print '0.   -----Down Sampling Time Elapsed: %s  \n \n' % timedelta(seconds=time.time()-t)

                t1 = time.time()
                ants_linear_registration(orig_template_DS, orig_input_image_DS)

            else:

                t1 = time.time()
                ants_linear_registration(orig_template, orig_input_image)

            print '1.   ----- ants_linear_registration Time Elapsed : %s  \n \n' % timedelta(seconds=time.time()-t1)
            print '1.   ----- ants_linear_registration Time Elapsed : %s  \n \n' % timedelta(seconds=time.time()-t)

            mask_input = os.path.join(os.path.dirname(orig_input_image), 'mask_inp.nii.gz')

            print "2.   Transform mask from template space to input space \n"
            t1 = time.time()
            # Transform mask from template space to input space
            ants_WarpImageMultiTransform(mask, mask_input, orig_input_image)
            print '2.   ----- Transform mask from template space to input space Time Elapsed: %s  \n \n' % timedelta(seconds=time.time()-t1)
            print '2.   ----- Transform mask from template space to input space Time Elapsed: %s  \n \n' % timedelta(seconds=time.time()-t)
        else:
            mask_input = os.path.join(os.path.dirname(orig_input_image), 'mask_inp.nii.gz')

        file_name = os.path.basename(orig_input_image)
        index_of_dot = file_name.index('.')
        file_name_without_extension = file_name[:index_of_dot]
        input_image = os.path.join(os.path.dirname(orig_input_image), 'crop_'+file_name_without_extension+'.nii.gz')


        # Cropping input using this mask
        print "3.   Cropping input using this mask \n"
        t1 = time.time()
        parallel_command(crop_by_mask(orig_input_image, input_image, mask_input))
        print '3.   ----- Cropping input using this mask Time Elapsed: %s  \n \n' % timedelta(seconds=time.time()-t1)
        print '3.   ----- Cropping input using this mask Time Elapsed: %s  \n \n' % timedelta(seconds=time.time()-t)

    print('-------------*********output_path************------------')
    print(output_path)


    # FSL automatically converts .nii to .nii.gz
    sanitized_image = os.path.join(temp_path, os.path.basename(input_image) + ('.gz' if input_image.endswith('.nii') else ''))
    if not os.path.exists(sanitized_image):
        input_image = sanitize_input(input_image, sanitized_image, parallel_command)
        if args.right:

            print "4.   lipping along L-R \n"
            t1 = time.time()
            flip_lr(input_image, input_image, parallel_command)
            print '4.   ----- lipping along L-R Time Elapsed: %s  \n \n' % timedelta(seconds=time.time()-t1)
            print '4.   ----- lipping along L-R Time Elapsed: %s  \n \n' % timedelta(seconds=time.time()-t)

        print "5.   Correcting bias \n"
        t1 = time.time()
        bias_correct(input_image, input_image, **exec_options)
        print '5.   ----- Correcting bias Time Elapsed: %s  \n \n' % timedelta(seconds=time.time()-t1)
        print '5.   ----- Correcting bias Time Elapsed: %s  \n \n' % timedelta(seconds=time.time()-t)


    else:
        print 'Skipped, using %s' % sanitized_image
        input_image = sanitized_image


    print "5b.   ants_nonlinear_registration \n"
    if args.forcereg or not check_warps(warp_path):
        if args.warp:
            print 'Saving output as %s' % warp_path
        else:
            warp_path = os.path.join(temp_path, tail)
            print 'Saving output to temporary path.'
        ants_nonlinear_registration(template, input_image, warp_path, **exec_options)
    else:
        print 'Skipped, using %sInverseWarp.nii.gz and %sAffine.txt' % (warp_path, warp_path)

    print('-------------*********output_path************------------')
    print(output_path)

    print '5b.   ----- ants_nonlinear_registration Time Elapsed: %s  \n \n' % timedelta(seconds=time.time()-t1)
    print '5b.   ----- ants_nonlinear_registration Time Elapsed: %s  \n \n' % timedelta(seconds=time.time()-t)

    # generating the warped output
    print "6.   Warping prior labels and images \n \n"
    t1 = time.time()
    registered = os.path.join(temp_path, 'registered.nii.gz')
    cmd = 'WarpImageMultiTransform 3 %s %s -R %s %sWarp.nii.gz %sAffine.txt' % (input_image, registered, template, warp_path, warp_path)
    parallel_command(cmd)

    print '6.   --- Warping prior labels and images Elapsed: %s  \n \n' % timedelta(seconds=time.time()-t1)
    print '6.   --- Warping prior labels and images Elapsed: %s  \n \n' % timedelta(seconds=time.time()-t)

    t1 = time.time()
    # TODO should probably use output from warp_atlas_subject instead of hard coding paths in create_atlas
    # TODO make this more parallel
    warped_labels = pool.map(partial(
        warp_atlas_subject,
        path=prior_path,
        # TODO cleanup this hack to always have whole thalamus so can estimate mask
        labels=set(labels + ['1-THALAMUS']),
        input_image=input_image,
        input_transform_prefix=warp_path,
        output_path=temp_path,
        BiosCorrectNameV=BiosCorrectName,
        exec_options=exec_options,
    ), subjects)
    warped_labels = {label: {subj: d[label] for subj, d in zip(subjects, warped_labels)} for label in warped_labels[0]}
    # # print '--- Forming subject-registered atlases. --- Elapsed: %s' % timedelta(seconds=time.time()-t)
    # atlases = pool.map(partial(create_atlas, path=temp_path, subjects=subjects, target='', echo=exec_options['echo']),
    # [{'label': label, 'output_atlas': os.path.join(temp_path, label+'_atlas.nii.gz')} for label in warped_labels])
    # atlases = dict(zip(warped_labels, zip(*atlases)[0]))
    # atlas_image = atlases[BiosCorrectName]
    atlas_images = warped_labels[BiosCorrectName].values()

    # FIXME use whole-brain template registration optimized parameters instead, these are from crop pipeline
    optimal_picsl = optimal['PICSL']
    # for k, v in warped_labels.iteritems():
    #     print k, v
    # for label in labels:
    #     print optimal_picsl[label]
    if args.jointfusion:
        pool.map(partial(label_fusion_picsl, input_image, atlas_images),
                 [dict(
                     atlas_labels=warped_labels[label].values(),
                     output_label=os.path.join(temp_path, label + '.nii.gz'),
                     rp=optimal_picsl[label]['rp'],
                     rs=optimal_picsl[label]['rs'],
                     beta=optimal_picsl[label]['beta'],
                     **exec_options
                 ) for label in labels])
    elif args.majorityvoting:
        pool.map(partial(label_fusion_majority),
                 [dict(
                     atlas_labels=warped_labels[label].values(),
                     output_label=os.path.join(temp_path, label + '.nii.gz'),
                     rp=optimal_picsl[label]['rp'],
                     rs=optimal_picsl[label]['rs'],
                     beta=optimal_picsl[label]['beta'],
                     **exec_options
                 ) for label in labels])
    else:
        # Estimate mask to restrict computation
        mask = os.path.join(temp_path, 'mask.nii.gz')
        check_run(
            mask,
            conservative_mask,
            warped_labels['1-THALAMUS'].values(),
            mask,
            dilation=10,
        )
        pool.map(partial(label_fusion_picsl_ants, input_image, atlas_images),
                 [dict(
                     atlas_labels=warped_labels[label].values(),
                     output_label=os.path.join(temp_path, label + '.nii.gz'),
                     rp=optimal_picsl[label]['rp'],
                     rs=optimal_picsl[label]['rs'],
                     beta=optimal_picsl[label]['beta'],
                     mask=mask,
                     **exec_options
                 ) for label in labels])

    print '7. --- Performing label fusion. Elapsed: %s  \n \n' % timedelta(seconds=time.time() - t1)
    print '7. --- Performing label fusion. Elapsed: %s  \n \n' % timedelta(seconds=time.time() - t)

    # STEPS
    # pool_small.map(partial(label_fusion, input_image=input_image, image_atlas=atlases[BiosCorrectName], echo=exec_options['echo']),
    #     [{
    #         'label_atlas': atlases[label],
    #         'output_label': os.path.join(output_path, label+'.nii.gz'),
    #         'sigma': optimal_steps[label]['steps_sigma'],
    #         'X': optimal_steps[label]['steps_X'],
    #         'mrf': optimal_steps[label]['steps_mrf'],
    #     } for label in labels]
    # )
    # for label in labels:
    #     print {
    #         'label': label,
    #         'sigma': optimal_steps[label]['steps_sigma'],
    #         'X': optimal_steps[label]['steps_X'],
    #         'mrf': optimal_steps[label]['steps_mrf'],
    #     }
    #     partial_fusion = partial(label_fusion, input_image=input_image, image_atlas=atlases[BiosCorrectName], echo=exec_options['echo'])
    #     label_fusion_args = {
    #         'label_atlas': atlases[label],
    #         'output_label': os.path.join(output_path, label+'.nii.gz'),
    #         'sigma': optimal_steps[label]['steps_sigma'],
    #         'X': optimal_steps[label]['steps_X'],
    #         'mrf': optimal_steps[label]['steps_mrf'],
    #     }
    #     partial_fusion(**label_fusion_args)

    files = [(os.path.join(temp_path, label + '.nii.gz'), os.path.join(output_path, label + '.nii.gz')) for label in labels]
    if args.right:
        pool.map(flip_lr, files)
        files = [(os.path.join(output_path, label + '.nii.gz'), os.path.join(output_path, label + '.nii.gz')) for label in labels]

    t1 = time.time()
    # Resort output to original ordering
    pool.map(parallel_command,
        ['%s %s %s %s' % (os.path.join(this_path, 'swapdimlike.py'), in_file, orig_input_image, out_file) for in_file, out_file in files])

    print '8. --- Resort output to original ordering. Elapsed: %s  \n \n' % timedelta(seconds=time.time() - t1)
    print '8. --- Resort output to original ordering. Elapsed: %s  \n \n' % timedelta(seconds=time.time() - t)

    # get the vlp file path for splitting
    print('-------------*********output_path************------------')
    print(output_path)
    vlp_file = os.path.join(output_path, '6-VLP.nii.gz')
    print('-------------*********************------------')
    print('vlp_file:   ' + vlp_file)


    # Re-orient to standard space - LR PA IS format
    san_vlp_file = os.path.join(output_path, 'san_6-VLP.nii.gz')

    t1 = time.time()
    input_image1 = sanitize_input(vlp_file, san_vlp_file, parallel_command)
    print '9. --- sanitize_input. Elapsed: %s  \n \n' % timedelta(seconds=time.time() - t1)
    print '9. --- sanitize_input. Elapsed: %s  \n \n' % timedelta(seconds=time.time() - t)

    # get the sanitized vlp for processing
    input_nii = nibabel.load(input_image1)
    data = input_nii.get_data()
    hdr = input_nii.get_header()
    affine = input_nii.get_affine()

    t1 = time.time()
    # Coronal axis for RL PA IS orientation
    vlps = split_roi(data, None, 2)
    for fname, sub_vlp in zip(['6_VLPv.nii.gz', '6_VLPd.nii.gz'], vlps):
        output_nii = nibabel.Nifti1Image(sub_vlp, affine, hdr)
        output_nii.to_filename(os.path.join(os.path.dirname(out_file), fname))
    print '10. --- Coronal axis for RL PA IS orientation. Elapsed: %s  \n \n' % timedelta(seconds=time.time() - t1)
    print '10. --- Coronal axis for RL PA IS orientation. Elapsed: %s  \n \n' % timedelta(seconds=time.time() - t)

    print '10.--- Finished --- Elapsed: %s' % timedelta(seconds=time.time() - t)
# """

if __name__ == '__main__':
    args = parser.parse_args()
    # print args
    # exec_options.update({'debug': args.debug, 'verbose': args.verbose})
    exec_options = {'echo': False, 'suppress': True}
    if args.verbose:
        exec_options['verbose'] = True
    if args.debug:
        print 'Debugging mode forces serial execution.'
        # exec_options['echo'] = True
        args.processes = 1
    parallel_command = partial(parallel.command, **exec_options)
    pool = parallel.BetterPool(args.processes)
    print 'Running with %d processes.' % pool._processes
    # TODO don't hard code this number of processors
    # pool_small = parallel.BetterPool(4)
    # TODO Add path of script to command()
    # os.environ['PATH'] += os.pathsep + os.path.abspath(os.path.dirname(sys.argv[0]))
    if args.tempdir:
        temp_path = args.tempdir
        if not os.path.exists(temp_path):
            print 'Making %s' % os.path.abspath(temp_path)
            os.makedirs(temp_path)
    else:
        temp_path = tempfile.mkdtemp(dir=os.path.dirname(args.output_path))
    try:
        main(args, temp_path, pool)
    finally:
        pool.close()
        # Clean up temp folders
        if not args.debug and not args.tempdir:
            try:
                rmtree(temp_path)
            except OSError as exc:
                if exc.errno != 2:  # Code 2 - no such file or directory
                    raise
