#!/bin/bash

set -x

ROOT=/home/larsonke # this is where you cloned the git repository into
WDIR=$ROOT/Code/Synthetic-Atrophy-For-Longtidunal-Cortical-Surface-Analyses
mkdir -p $WDIR

# External tools you will need to run this example:
c3d=$ROOT/bin/c3d # (https://sourceforge.net/p/c3d/git/ci/master/tree/)
slicer=$ROOT/Libraries/Slicer-4.10.1-linux-amd64/Slicer # (https://download.slicer.org/)
greedy=$ROOT/bin/greedy # (https://sites.google.com/view/greedyreg/)
imagemath=$ROOT/bin/ImageMath # (https://github.com/NIRALUser/niral_utilities/tree/master/ImageMath)
holedetection=$ROOT/bin/LOGB-Tools/HoleDetection # (https://www.insight-journal.org/browse/publication/43)


# Synthetic-Atrophy tools (you will have to manually add the paths for the executables once you compile them all)
addFloatArray=$ROOT/bin/VTK-Tools/parcellateSingleSurface
distanceBetweenSurfaces=$ROOT/Code/Synthetic-Atrophy-dev/Distance-Between-Surfaces-build/distanceBetweenSurfaces
fillholes=$ROOT/bin/LOGB-Tools/FillHoles
getBlurMask=$ROOT/Code/Synthetic-Atrophy-dev/Get-Blur-Mask-build/getBlurMask
getTrimRegion=$ROOT/Code/Synthetic-Atrophy-dev/Get-Trim-Region-build/getTrimRegion
maskField=$ROOT/Code/Synthetic-Atrophy-dev/Mask-Field-build/maskField
padImage=$ROOT/Code/Synthetic-Atrophy-dev/Pad-Image-build/padImage
padVectorImage=$ROOT/Code/Synthetic-Atrophy-dev/Pad-Vector-Image-build/padVectorImage
parcellateAtrophyRegion=$ROOT/Code/Synthetic-Atrophy-dev/Parcellate-Atrophy-Region-build/parcellateAtrophyRegion
performErosion=$ROOT/Code/Synthetic-Atrophy-dev/Perform-Erosion-build/performErosion
readMeshArray=$ROOT/bin/VTK-Tools/readMeshArray
trimFloatImage=$ROOT/Code/Synthetic-Atrophy-dev/Trim-Float-Image-build/trimFloatImage
trimImage=$ROOT/Code/Synthetic-Atrophy-dev/Trim-Image-build/trimImage
warpSurfaces=$ROOT/Code/Synthetic-Atrophy-dev/Warp-Surfaces-build/warpSurfs
VTK2VTP=$ROOT/bin/Switch-File-Types/VTK2VTP
VTP2VTK=$ROOT/bin/Switch-File-Types/VTP2VTK




# Crops an image and upsamples to 400%
function highRes_crop(){
    original=$1
    rescaleFactor=$2
    cropRegion_file=$3
    outbase=$4

    trim=${outbase}_trim.nii.gz
    cropRegion_info=$(sed -n 1p $cropRegion_file)
    $trimImage $original $trim $cropRegion_info

    cropRegion_origin=$(sed -n 2p $cropRegion_file)

    trim_highRes=${outbase}_trim_highRes.nii.gz
    c3d -interpolation NearestNeighbor $trim -resample ${rescaleFactor}% -holefill 1 0 -o $trim_highRes
    cropRegion_origin=$(sed -n 2p $cropRegion_file)
}


# Get the high resolution cropped data
function Get-HighRes-Mask(){
    SUBJECT=$1
    label=$2
    RAW_DIR=$DATA_ROOT/Raw/$PROJECT_NAME/$SUBJECT/Extras
    WDIR=$DATA_ROOT/Processed/$PROJECT_NAME/$SUBJECT/$label
    MASK_DIR=$DATA_ROOT/Processed/$PROJECT_NAME/$SUBJECT/$label/HighRes-Mask
    ROI_DIR=$DATA_ROOT/Processed/$PROJECT_NAME/$SUBJECT/$label/ROI
    mkdir -p $MASK_DIR $ROI_DIR

    crop_region_info_file=$MASK_DIR/crop_region_info.txt
    if [ -f $crop_region_info_file ] ; then
	rm $crop_region_info_file
    fi
    
    aparc_aseg=$RAW_DIR/aparc+aseg_native.nii.gz
    label_mask=$ROI_DIR/label_original.nii.gz
    c3d $aparc_aseg -threshold $label $label 1 0 -holefill 1 0 -o $label_mask
    
    brain_mask=$RAW_DIR/skullstrip_mask.nii.gz
    brain_std=$RAW_DIR/brain_std.nii.gz
    
    total_mask=$MASK_DIR/total_mask.nii.gz
    $getBlurMask $label_mask $brain_std $brain_mask $total_mask 0.1
    trim_info=$($getTrimRegion $total_mask 4)
    echo $trim_info >> $crop_region_info_file
    
    total_mask_trim=$MASK_DIR/total_mask_trim.nii.gz
    total_mask_trim_highRes=$MASK_DIR/total_mask_trim_highRes.nii.gz
    $trimImage $total_mask $total_mask_trim $trim_info
    ref_info=$(c3d $total_mask_trim -info)
    x=$(echo $ref_info | cut -d "[" -f 3 | cut -d " " -f 1)
    y=$(echo $ref_info | cut -d "[" -f 3 | cut -d " " -f 2)
    z=$(echo $ref_info | cut -d "[" -f 3 | cut -d " " -f 3 | cut -d "]" -f 1)
    echo $x $y $z >> $crop_region_info_file
 
    c3d -interpolation NearestNeighbor $total_mask_trim -resample 400% -holefill 1 0 -o $total_mask_trim_highRes

    # High res + crop
    highRes_crop $brain_mask 400 $crop_region_info_file $MASK_DIR/brain_mask
    highRes_crop $aparc_aseg 400 $crop_region_info_file $MASK_DIR/aparc+aseg_native
    highRes_crop $label_mask 400 $crop_region_info_file $ROI_DIR/label_original
    
    c3d $MASK_DIR/aparc+aseg_native_trim.nii.gz -threshold 1 inf 1 0 -sdt -o $MASK_DIR/brain_std_trim.nii.gz
    c3d $MASK_DIR/aparc+aseg_native_trim_highRes.nii.gz -threshold 1 inf 1 0 -sdt -o $MASK_DIR/brain_std_trim_highRes.nii.gz    

    hem_id=$(echo $label | cut -b 1)
    if [ $hem_id == 1 ] ; then
        wm_val=2
    elif [ $hem_id == 2 ] ; then
        wm_val=41
    fi
    c3d $MASK_DIR/aparc+aseg_native_trim.nii.gz -threshold $wm_val $wm_val 1 0 -o $MASK_DIR/wm_trim.nii.gz
    c3d $MASK_DIR/aparc+aseg_native_trim_highRes.nii.gz -threshold $wm_val $wm_val 1 0 -o $MASK_DIR/wm_trim_highRes.nii.gz

    c3d $MASK_DIR/aparc+aseg_native_trim_highRes.nii.gz -threshold 1 inf 1 0 -o $ROI_DIR/label_original_joint_trim_highRes.nii.gz
    
    $performErosion $total_mask_trim $MASK_DIR/wm_trim.nii.gz $MASK_DIR/total_mask_eroded_trim.nii.gz $blurKernel 0
    c3d $MASK_DIR/total_mask_eroded_trim.nii.gz $ROI_DIR/label_original_trim.nii.gz -add -threshold 1 inf 1 0 -o $MASK_DIR/total_mask_eroded_trim.nii.gz

}





function Remove-Handles(){
    dir=${1?}
    mkdir -p $dir

    image0=${2?}
    mesh0=${3?}
    outbase=${4?}

    cp $image0 $dir/image_0.nii.gz
    cp $mesh0 $dir/mesh_0.vtk

    flip=0

    for i in `seq 1 10`
    do
      echo $i
      prev=`expr $i - 1`
      holes=$dir/holes_$i.vtk
      mesh=$dir/mesh_${prev}.vtk
      newmesh=$dir/mesh_$i
      image=$dir/image_$prev.nii.gz
      newimage=$dir/image_$i.nii.gz
      finalimage=${outbase}holeFill.nii.gz
      finalmesh=${outbase}holeFill

      flip=`expr 1 - $flip`

      $holedetection 0 $mesh $holes

      holecount=`grep LINES $holes | wc -l `
      if [ $holecount -eq 0 ]
      then
         $imagemath $image -outfile $finalimage -conComp 1
         $slicer --launch ModelMaker $finalimage --generateAll --smooth 10 --decimate 0 --pointnormals --name $finalmesh >> /dev/null
         mv ${finalmesh}_1.vtk $finalmesh.vtk
         break
      else
        echo "holes found"
        $fillholes $holes $image $newimage $flip
	$slicer --launch ModelMaker $newimage --generateAll --smooth 0 --decimate 0 --pointnormals --name $newmesh >> /dev/null
        mv ${newmesh}_1.vtk $newmesh.vtk
      fi
    done
    rm *mrml
    rm -r $dir
}



function Get-Mesh(){
    inputImage=$1
    outputMesh=$2
    DIR=$3

    mkdir -p $DIR/Dump $DIR/Handle-Removal
    modelMakerMesh_base=$DIR/Dump/mesh
    removeHandles_base=$DIR/Dump/mesh_noHandles_

    $slicer --launch ModelMaker $inputImage --generateAll --smooth 10 --decimate 0 --pointnormals --name $modelMakerMesh_base >> /dev/null
    Remove-Handles $DIR/Handle-Removal $inputImage ${modelMakerMesh_base}_1.vtk $removeHandles_base
    $VTK2VTP ${removeHandles_base}holeFill.vtk $outputMesh

    rm -r $DIR/Dump
}


function Change-To-Working-Space(){
    input=$1
    output=$2

    $slicer --launch OrientScalarVolume --orientation LPI $input $output >> /dev/null
    $c3d $output -origin-voxel 0x0x0vox -o $output
}





function main {
    ######################################
    ##    Set up all the data to run    ##
    ######################################    

    # Set up directories
    DATA_DIR=$WDIR/Data
    INPUT_DIR=$DATA_DIR/Input
    ROI_DIR=$DATA_DIR/ROI
    MASK_DIR=$DATA_DIR/HighRes-Mask
    TRANSFORM_DIR=$DATA_DIR/Transforms
    SURFACE_DIR=$DATA_DIR/Surfaces
    RESULTS_DIR=$DATA_DIR/Results
    mkdir -p $ROI_DIR $MASK_DIR $TRANSFORM_DIR $SURFACE_DIR $RESULTS_DIR

    # Provided data
    T1=$INPUT_DIR/T1.nii.gz
    FLAIR=$INPUT_DIR/FLAIR.nii.gz
    parcMap=$INPUT_DIR/inputParcellation.nii.gz
    ribbon=$INPUT_DIR/corticalRibbon.nii.gz
    skullStripMask=$INPUT_DIR/skullStripMask.nii.gz

    # Get remaining inputs from provided data
    wm=$INPUT_DIR/wm.nii.gz
    WM_lh_mask=$SURFACE_DIR/WM_lh.nii.gz
    WM_rh_mask=$SURFACE_DIR/WM_rh.nii.gz
    fullBrain_mask=$DATA_DIR/fullBrainMask.nii.gz
    GM_lh_mask=$SURFACE_DIR/GM_lh.nii.gz
    GM_rh_mask=$SURFACE_DIR/GM_rh.nii.gz
    sdt=$INPUT_DIR/brainSDT.nii.gz
    
    if false ; then
    $c3d $ribbon -as ribbon -thresh 2 2 1 0 -as LWM -o $WM_lh_mask \
	 -push ribbon -thresh 41 41 1 0 -as RWM -o $WM_rh_mask \
	 -push LWM -push RWM -add -o $wm \
	 -push ribbon -thresh 2 3 1 0 -as Lfull -o $GM_lh_mask \
	 -push ribbon -thresh 41 42 1 0 -as Rfull -o $GM_rh_mask \
	 -push Lfull -push Rfull -add -as M -o $fullBrain_mask \
	 -push M -sdt -o $sdt
    fi

    # Convert to working space
    parcMapLPI=$INPUT_DIR/inputParcellation_LPI.nii.gz
    skullStripMaskLPI=$INPUT_DIR/skullStripMask_LPI.nii.gz
    wmLPI=$INPUT_DIR/wm_LPI.nii.gz
    sdtLPI=$INPUT_DIR/brainSDT_LPI.nii.gz
    
    if false ; then
    Change-To-Working-Space $parcMap $parcMapLPI
    Change-To-Working-Space $skullStripMask $skullStripMaskLPI
    Change-To-Working-Space $wm $wmLPI
    Change-To-Working-Space $sdt $sdtLPI
    fi

    
    
    ##############################################################
    ##    Get high-resolution images cropped to selected ROI    ##
    ##############################################################

    # Get label image
    label=1030 # label for atrophy induction (this is the ROI we used as an example in our publication)
    labelMask_original=$ROI_DIR/label_original.nii.gz
    if false ; then
    $c3d $parcMapLPI -threshold $label $label 1 0 -holefill 1 0 -o $labelMask_original
    fi
    # Blur mask (mask used for blurring/smoothing the atrophy deformation field)
    blurMask=$MASK_DIR/blurMask.nii.gz
    if false ; then
    $getBlurMask $labelMask_original $sdtLPI $skullStripMaskLPI $blurMask 0.1 
    fi
    
    # Get cropping data from blur mask (we are going to crop all the images to the blur mask + a 4 voxel buffer on all sides) 
    crop_region_info_file=$MASK_DIR/crop_region_info.txt
    trim_info=$($getTrimRegion $blurMask 4)
    if false ; then
    echo $trim_info >> $crop_region_info_file
    fi
    blurMask_trim=$MASK_DIR/blurMask_trim.nii.gz
    $trimImage $blurMask $blurMask_trim $trim_info
    ref_info=$($c3d $blurMask_trim -info)
    x=$(echo $ref_info | cut -d "[" -f 3 | cut -d " " -f 1)
    y=$(echo $ref_info | cut -d "[" -f 3 | cut -d " " -f 2)
    z=$(echo $ref_info | cut -d "[" -f 3 | cut -d " " -f 3 | cut -d "]" -f 1)
    echo $x $y $z >> $crop_region_info_file


    # Upsample and crop image data
    blurMask_trim_highRes=$MASK_DIR/blurMask_trim_highRes.nii.gz
    if false ; then
    c3d -interpolation NearestNeighbor $blurMask_trim -resample 400% -holefill 1 0 -o $blurMask_trim_highRes

    highRes_crop $skullStripMaskLPI 400 $crop_region_info_file $MASK_DIR/skullStripMask
    highRes_crop $parcMapLPI 400 $crop_region_info_file $MASK_DIR/inputParcellation
    highRes_crop $labelMask_original 400 $crop_region_info_file $ROI_DIR/label_original
    highRes_crop $wmLPI 400 $crop_region_info_file $MASK_DIR/wm
    highRes_crop $sdtLPI 400 $crop_region_info_file $MASK_DIR/sdt
    fi



    #########################################
    ##    Create synthetic atrophy data    ##
    #########################################

    # Perform erosion
    mkdir -p $ROI_DIR/Erosion-Iterations
    labelMask_atrophy=$ROI_DIR/label_atrophy_trim_highRes.nii.gz
    nIterations=4 # number of iterations for erosion via binary morphology operations

    if false ; then
    n=1
    while [ $n -le $nIterations ] ; do
	if [ $n == 1 ] ; then
	    prev=$ROI_DIR/label_original_trim_highRes.nii.gz
	else
	    let n_prev=$n-1
	    prev=$ROI_DIR/Erosion-Iterations/${n_prev}.nii.gz
	fi
	next=$ROI_DIR/Erosion-Iterations/${n}.nii.gz
	$performErosion $prev $MASK_DIR/wm_trim_highRes.nii.gz $next 1 0
	let n=$n+1
    done
    

    cp $next $labelMask_atrophy
    #rm $ROI_DIR/Erosion-Iterations
    fi


    # Get full brain mask for atrophied timepoint
    joint_labelMask_original=$ROI_DIR/label_original_joint_trim_highRes.nii.gz
    c3d $MASK_DIR/inputParcellation_trim_highRes.nii.gz -as P -threshold 1 inf 1 0 -o $joint_labelMask_original
    
    joint_labelMask_atrophy=$ROI_DIR/label_atrophy_joint_trim_highRes.nii.gz
    let lower_threshold=$label-1
    let upper_threshold=$label+1
    if false ; then
    c3d $MASK_DIR/inputParcellation_trim_highRes.nii.gz -as P -threshold 1 $lower_threshold 1 0 -popas LM \
	-push P -threshold $upper_threshold inf 1 0 -popas UM \
	-push LM -push UM -add $labelMask_atrophy -add -o $joint_labelMask_atrophy
    fi


    # Register original --> atrophied
    field_outbase=$TRANSFORM_DIR/field
    field_inverse_outbase=$TRANSFORM_DIR/field_inverse
    if false ; then
    $greedy -d 3 -m MSQ -i $joint_labelMask_atrophy $joint_labelMask_original -o ${field_outbase}_full_trim_highRes.nii.gz \
	    -oinv ${field_inverse_outbase}_full_trim_highRes.nii.gz \
	    -mm $joint_labelMask_original -n 100x100x50x10 >> /dev/null
    $c3d -mcs ${field_outbase}_full_trim_highRes.nii.gz -foreach -resample 25% -origin 0x0x0mm -endfor -omc ${field_outbase}_full_trim.nii.gz
    $c3d -mcs ${field_inverse_outbase}_full_trim_highRes.nii.gz -foreach -resample 25% -origin 0x0x0mm -endfor -omc ${field_inverse_outbase}_full_trim.nii.gz
    fi

    # Mask/blur fields
    pad_info=$(sed -n 1p $crop_region_info_file)
    mkdir -p $TRANSFORM_DIR/Dump

    if false ; then
    c3d $MASK_DIR/sdt_trim.nii.gz -origin 0x0x0mm -o $TRANSFORM_DIR/Dump/brainSDT_trim.nii.gz
    c3d $ROI_DIR/label_original_trim.nii.gz -origin 0x0x0mm -o $TRANSFORM_DIR/Dump/label_original_trim.nii.gz
    c3d $MASK_DIR/skullStripMask_trim.nii.gz -origin 0x0x0mm -o $TRANSFORM_DIR/Dump/skullStripMask_trim.nii.gz
    c3d $MASK_DIR/wm_trim.nii.gz -origin 0x0x0mm -o $TRANSFORM_DIR/Dump/wm_trim.nii.gz

    $maskField $TRANSFORM_DIR/Dump/brainSDT_trim.nii.gz $TRANSFORM_DIR/Dump/label_original_trim.nii.gz $TRANSFORM_DIR/Dump/skullStripMask_trim.nii.gz ${field_outbase}_full_trim.nii.gz \
	       $TRANSFORM_DIR/Dump/wm_trim.nii.gz 0.4 $TRANSFORM_DIR/Dump/points.vtp ${field_outbase}_maskBlur_trim.nii.gz
    $padVectorImage $MASK_DIR/blurMask.nii.gz ${field_outbase}_maskBlur_trim.nii.gz ${field_outbase}_maskBlur_pad.nii.gz $pad_info
    c3d -mcs ${field_outbase}_maskBlur_pad.nii.gz -foreach -origin 0x0x0vox -endfor -omc ${field_outbase}_maskBlur_pad.nii.gz
    
    $maskField $TRANSFORM_DIR/Dump/brainSDT_trim.nii.gz $TRANSFORM_DIR/Dump/label_original_trim.nii.gz $TRANSFORM_DIR/Dump/skullStripMask_trim.nii.gz ${field_inverse_outbase}_full_trim.nii.gz \
	       $TRANSFORM_DIR/Dump/wm_trim.nii.gz 0.4 $TRANSFORM_DIR/Dump/points.vtp ${field_inverse_outbase}_maskBlur_trim.nii.gz
    $padVectorImage $MASK_DIR/blurMask.nii.gz ${field_inverse_outbase}_maskBlur_trim.nii.gz ${field_inverse_outbase}_maskBlur_pad.nii.gz $pad_info
        c3d -mcs ${field_inverse_outbase}_maskBlur_pad.nii.gz -foreach -origin 0x0x0vox -endfor -omc ${field_inverse_outbase}_maskBlur_pad.nii.gz

    rm -r $TRANSFORM_DIR/Dump
    fi

    
    
    #################################
    ##    Create synthetic data    ##
    #################################
    
    # Apply warp to T1 image
    cp $T1 $RESULTS_DIR/T1.nii.gz
    T1_atrophy=$RESULTS_DIR/T1_atrophy.nii.gz

    if false ; then
    $greedy -d 3 -rf $T1 -rm $T1 $T1_atrophy -r ${field_outbase}_maskBlur_pad.nii.gz
    fi
    

    # Register FLAIR to T1 image
    FLAIR_T1space=$RESULTS_DIR/FLAIR_T1space.nii.gz
    FLAIR_atrophy=$RESULTS_DIR/FLAIR_T1space_atrophy.nii.gz
    moments=$DATA_DIR/FLAIR-Registration/moments.mat
    rigid=$DATA_DIR/FLAIR-Registration/rigid.mat
    mkdir -p $DATA_DIR/FLAIR-Registration

    if false ; then
    $greedy -d 3 -i $T1 $FLAIR -o $moments -moments 1 >> /dev/null
    $greedy -d 3 -i $T1 $FLAIR -o $rigid -a -dof 6 -n 100x50x10 -m 4x4x4 -ia $moments >> /dev/null
    $greedy -d 3 -r $rigid -rf $T1 -rm $FLAIR $FLAIR_T1space >> /dev/null
    $greedy -d 3 -rf $T1 -rm $FLAIR_T1space $FLAIR_atrophy -r ${field_outbase}_maskBlur_pad.nii.gz

    fi
    
    #rm -r $$DATA_DIR/FLAIR-Registration

    

    ####################################
    ##    Measure thickness change    ##
    ####################################

    # Make label map for ROI
    fullBrainMask_parc=$DATA_DIR/parcellation_ROI.nii.gz
    if false ; then
    c3d $fullBrain_mask $labelMask_original -as ROI -scale -1 -add -popas M \
	-push ROI -dilate 1 4x4x4vox -push M -times -push M -add -push ROI -scale 3 -add -o $fullBrainMask_parc
    fi
    
    # Get meshes of original (unatrophied) data
    WM_lh_mesh=$SURFACE_DIR/WM_lh.vtp
    WM_rh_mesh=$SURFACE_DIR/WM_rh.vtp
    GM_lh_mesh=$SURFACE_DIR/GM_lh.vtp
    GM_rh_mesh=$SURFACE_DIR/GM_rh.vtp

    if false ; then
	Get-Mesh $WM_lh_mask $WM_lh_mesh $SURFACE_DIR
	Get-Mesh $WM_rh_mask $WM_rh_mesh $SURFACE_DIR
	Get-Mesh $GM_lh_mask $GM_lh_mesh $SURFACE_DIR
	Get-Mesh $GM_rh_mask $GM_rh_mesh $SURFACE_DIR
    fi

    
    # Parcellate meshes
    WM_lh_mesh_original=$SURFACE_DIR/WM_lh_original.vtp
    WM_rh_mesh_original=$SURFACE_DIR/WM_rh_original.vtp
    GM_lh_mesh_original=$SURFACE_DIR/GM_lh_original.vtp
    GM_rh_mesh_original=$SURFACE_DIR/GM_rh_original.vtp

    WM_lh_mesh_parc=$RESULTS_DIR/WM_lh_atrophyLabelParcellation.txt
    WM_rh_mesh_parc=$RESULTS_DIR/WM_rh_atrophyLabelParcellation.txt
    GM_lh_mesh_parc=$RESULTS_DIR/GM_lh_atrophyLabelParcellation.txt
    GM_rh_mesh_parc=$RESULTS_DIR/GM_rh_atrophyLabelParcellation.txt

    if false ; then
    $parcellateAtrophyRegion $WM_lh_mesh $GM_lh_mesh $fullBrainMask_parc $WM_lh_mesh_original $GM_lh_mesh_original 0
    $parcellateAtrophyRegion $WM_rh_mesh $GM_rh_mesh $fullBrainMask_parc $WM_rh_mesh_original $GM_rh_mesh_original 0
    fi

    
    # Apply warp to surfaces
    WM_lh_mesh_atrophy=$SURFACE_DIR/WM_lh_atrophy.vtp
    WM_rh_mesh_atrophy=$SURFACE_DIR/WM_rh_atrophy.vtp
    GM_lh_mesh_atrophy=$SURFACE_DIR/GM_lh_atrophy.vtp
    GM_rh_mesh_atrophy=$SURFACE_DIR/GM_rh_atrophy.vtp

    if false ; then
    $warpSurfaces $WM_lh_mesh_original ${field_inverse_outbase}_maskBlur_pad.nii.gz $MASK_DIR/blurMask.nii.gz $WM_lh_mesh_atrophy 0 1
    $warpSurfaces $WM_rh_mesh_original ${field_inverse_outbase}_maskBlur_pad.nii.gz $MASK_DIR/blurMask.nii.gz $WM_rh_mesh_atrophy 0 1
    $warpSurfaces $GM_lh_mesh_original ${field_inverse_outbase}_maskBlur_pad.nii.gz $MASK_DIR/blurMask.nii.gz $GM_lh_mesh_atrophy 0 0
    $warpSurfaces $GM_rh_mesh_original ${field_inverse_outbase}_maskBlur_pad.nii.gz $MASK_DIR/blurMask.nii.gz $GM_rh_mesh_atrophy 0 0
    fi


    # Get distance between surfaces
    WM_lh_mesh_original_surfaceDistance=$SURFACE_DIR/WM_lh_original_surfaceDistance.vtp
    WM_rh_mesh_original_surfaceDistance=$SURFACE_DIR/WM_rh_original_surfaceDistance.vtp
    GM_lh_mesh_original_surfaceDistance=$SURFACE_DIR/GM_lh_original_surfaceDistance.vtp
    GM_rh_mesh_original_surfaceDistance=$SURFACE_DIR/GM_rh_original_surfaceDistance.vtp

    WM_lh_mesh_atrophy_surfaceDistance=$SURFACE_DIR/WM_lh_atrophy_surfaceDistance.vtp
    WM_rh_mesh_atrophy_surfaceDistance=$SURFACE_DIR/WM_rh_atrophy_surfaceDistance.vtp
    GM_lh_mesh_atrophy_surfaceDistance=$SURFACE_DIR/GM_lh_atrophy_surfaceDistance.vtp
    GM_rh_mesh_atrophy_surfaceDistance=$SURFACE_DIR/GM_rh_atrophy_surfaceDistance.vtp

    $distanceBetweenSurfaces $WM_lh_mesh_original $WM_lh_mesh_atrophy $WM_lh_mesh_original_surfaceDistance
    $distanceBetweenSurfaces $WM_lh_mesh_atrophy $WM_lh_mesh_original $WM_lh_mesh_atrophy_surfaceDistance
    $distanceBetweenSurfaces $WM_rh_mesh_original $WM_rh_mesh_atrophy $WM_rh_mesh_original_surfaceDistance
    $distanceBetweenSurfaces $WM_rh_mesh_atrophy $WM_rh_mesh_original $WM_rh_mesh_atrophy_surfaceDistance
    $distanceBetweenSurfaces $GM_lh_mesh_original $GM_lh_mesh_atrophy $GM_lh_mesh_original_surfaceDistance
    $distanceBetweenSurfaces $GM_lh_mesh_atrophy $GM_lh_mesh_original $GM_lh_mesh_atrophy_surfaceDistance
    $distanceBetweenSurfaces $GM_rh_mesh_original $GM_rh_mesh_atrophy $GM_rh_mesh_original_surfaceDistance
    $distanceBetweenSurfaces $GM_rh_mesh_atrophy $GM_rh_mesh_original $GM_rh_mesh_atrophy_surfaceDistance

    
    

    # Write everything out into a text file
    WM_lh_mesh_parc_txt=$RESULTS_DIR/WM_lh_atrophyLabelParcellation.txt
    WM_rh_mesh_parc_txt=$RESULTS_DIR/WM_rh_atrophyLabelParcellation.txt
    GM_lh_mesh_parc_txt=$RESULTS_DIR/GM_lh_atrophyLabelParcellation.txt
    GM_rh_mesh_parc_txt=$RESULTS_DIR/GM_rh_atrophyLabelParcellation.txt

    $readMeshArray $WM_lh_mesh_original $WM_lh_mesh_parc_txt "Atrophy-Region"
    $readMeshArray $WM_rh_mesh_original $WM_rh_mesh_parc_txt "Atrophy-Region"
    $readMeshArray $GM_lh_mesh_original $GM_lh_mesh_parc_txt "Atrophy-Region"
    $readMeshArray $GM_rh_mesh_original $GM_rh_mesh_parc_txt "Atrophy-Region"
    
    WM_lh_mesh_original_surfaceDistance_txt=$RESULTS_DIR/WM_lh_surfaceDistance_original2atrophy.txt
    WM_rh_mesh_original_surfaceDistance_txt=$RESULTS_DIR/WM_rh_surfaceDistance_original2atrophy.txt
    WM_rh_mesh_atrophy_surfaceDistance_txt=$RESULTS_DIR/WM_lh_surfaceDistance_atrophy2original.txt
    WM_rh_mesh_atrophy_surfaceDistance_txt=$RESULTS_DIR/WM_rh_surfaceDistance_atrophy2original.txt
    GM_lh_mesh_original_surfaceDistance_txt=$RESULTS_DIR/GM_lh_surfaceDistance_original2atrophy.txt
    GM_rh_mesh_original_surfaceDistance_txt=$RESULTS_DIR/GM_rh_surfaceDistance_original2atrophy.txt
    GM_rh_mesh_atrophy_surfaceDistance_txt=$RESULTS_DIR/GM_lh_surfaceDistance_atrophy2original.txt
    GM_rh_mesh_atrophy_surfaceDistance_txt=$RESULTS_DIR/GM_rh_surfaceDistance_atrophy2original.txt

    $readMeshArray $WM_lh_mesh_original_surfaceDistance $WM_lh_mesh_original_surfaceDistance_txt "Surface-Distance"
    $readMeshArray $WM_rh_mesh_original_surfaceDistance $WM_rh_mesh_original_surfaceDistance_txt "Surface-Distance"
    $readMeshArray $WM_lh_mesh_atrophy_surfaceDistance $WM_lh_mesh_atrophy_surfaceDistance_txt "Surface-Distance"
    $readMeshArray $WM_rh_mesh_atrophy_surfaceDistance $WM_rh_mesh_atrophy_surfaceDistance_txt "Surface-Distance"
    $readMeshArray $GM_lh_mesh_atrophy_surfaceDistance $GM_lh_mesh_original_surfaceDistance_txt "Surface-Distance"
    $readMeshArray $GM_rh_mesh_original_surfaceDistance $GM_rh_mesh_original_surfaceDistance_txt "Surface-Distance"
    $readMeshArray $GM_lh_mesh_atrophy_surfaceDistance $GM_lh_mesh_atrophy_surfaceDistance_txt "Surface-Distance"
    $readMeshArray $GM_rh_mesh_atrophy_surfaceDistance $GM_rh_mesh_atrophy_surfaceDistance_txt "Surface-Distance"
}



if [[ $1 ]]; then
  command=$1
  echo $1
  shift
  $command $@
else
    main
fi

