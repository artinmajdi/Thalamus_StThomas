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
from THOMAS_constants_v0 import image_name, orig_template, template_61, template_93, mask_61, mask_93, this_path, prior_path, subjects, roi, roi_choices, optimal


def warp_atlas_subject(subject, path, labels, input_image, input_transform_prefix, output_path, exec_options={}):
    """
    Warp a training set subject's labels to input_image.
    """
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
    output_labels['WMnMPRAGE_bias_corr'] = output_image = os.path.join(output_path, image_name)
    if not os.path.exists(output_labels['WMnMPRAGE_bias_corr']):
        print output_labels['WMnMPRAGE_bias_corr']
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
parser.add_argument('--output_path', help='specify a different output_path for the output file for single ROI or directory for multiple ROIs')
# TODO handle single roi, single output file case
# TODO fix verbose and debug
# TODO go back to shell=False for command to suppress output and then fix sanitize labels


def main(args, temp_path, pool):
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

    if args.algorithm == "v2":
        # Crop the input
        # Affine registering input to template
        ants_linear_registration(orig_template, orig_input_image)
        print "Completed a quick linear registration of input and full template"
        mask_input = os.path.join(os.path.dirname(orig_input_image), 'mask_inp.nii.gz')
        # Transform mask from template space to input space
        ants_WarpImageMultiTransform(mask, mask_input, orig_input_image)
        print "Completed transforming the mask from template space to input space"
        file_name = os.path.basename(orig_input_image)
        index_of_dot = file_name.index('.')
        file_name_without_extension = file_name[:index_of_dot]
        input_image = os.path.join(os.path.dirname(orig_input_image), 'crop_'+file_name_without_extension+'.nii.gz')
        # Cropping input using this mask
        parallel_command(crop_by_mask(orig_input_image, input_image, mask_input))
        print "Completed cropping the input"


    # FSL automatically converts .nii to .nii.gz
    sanitized_image = os.path.join(temp_path, os.path.basename(input_image) + ('.gz' if input_image.endswith('.nii') else ''))
    print '--- Reorienting image. --- Elapsed: %s' % timedelta(seconds=time.time()-t)
    if not os.path.exists(sanitized_image):
        input_image = sanitize_input(input_image, sanitized_image, parallel_command)
        if args.right:
            print '--- Flipping along L-R. --- Elapsed: %s' % timedelta(seconds=time.time()-t)
            flip_lr(input_image, input_image, parallel_command)
        print '--- Correcting bias. --- Elapsed: %s' % timedelta(seconds=time.time()-t)
        bias_correct(input_image, input_image, **exec_options)
    else:
        print 'Skipped, using %s' % sanitized_image
        input_image = sanitized_image


    print '--- Registering to mean brain template. --- Elapsed: %s' % timedelta(seconds=time.time()-t)
    if args.forcereg or not check_warps(warp_path):
        if args.warp:
            print 'Saving output as %s' % warp_path
        else:
            warp_path = os.path.join(temp_path, tail)
            print 'Saving output to temporary path.'
        ants_nonlinear_registration(template, input_image, warp_path, **exec_options)
    else:
        print 'Skipped, using %sInverseWarp.nii.gz and %sAffine.txt' % (warp_path, warp_path)

    # generating the warped output
    registered = os.path.join(temp_path, 'registered.nii.gz')
    cmd = 'WarpImageMultiTransform 3 %s %s -R %s %sWarp.nii.gz %sAffine.txt' % (input_image, registered, template, warp_path, warp_path)
    parallel_command(cmd)

    print '--- Warping prior labels and images. --- Elapsed: %s' % timedelta(seconds=time.time() - t)
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
        exec_options=exec_options,
    ), subjects)
    warped_labels = {label: {subj: d[label] for subj, d in zip(subjects, warped_labels)} for label in warped_labels[0]}
    # # print '--- Forming subject-registered atlases. --- Elapsed: %s' % timedelta(seconds=time.time()-t)
    # atlases = pool.map(partial(create_atlas, path=temp_path, subjects=subjects, target='', echo=exec_options['echo']),
    # [{'label': label, 'output_atlas': os.path.join(temp_path, label+'_atlas.nii.gz')} for label in warped_labels])
    # atlases = dict(zip(warped_labels, zip(*atlases)[0]))
    # atlas_image = atlases['WMnMPRAGE_bias_corr']
    atlas_images = warped_labels['WMnMPRAGE_bias_corr'].values()

    print '--- Performing label fusion. --- Elapsed: %s' % timedelta(seconds=time.time() - t)
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
    # STEPS
    # pool_small.map(partial(label_fusion, input_image=input_image, image_atlas=atlases['WMnMPRAGE_bias_corr'], echo=exec_options['echo']),
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
    #     partial_fusion = partial(label_fusion, input_image=input_image, image_atlas=atlases['WMnMPRAGE_bias_corr'], echo=exec_options['echo'])
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
    # Resort output to original ordering
    pool.map(parallel_command,
        ['%s %s %s %s' % (os.path.join(this_path, 'swapdimlike.py'), in_file, orig_input_image, out_file) for in_file, out_file in files])
    print '--- Finished --- Elapsed: %s' % timedelta(seconds=time.time() - t)


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
