#!/bin/bash

# Specify SLURM data
#SBATCH --mail-user=larsonke@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=0-20:00:00
#SBATCH --mem=20G
#SBATCH --output=out_SA


RUN_TYPE=accre
#RUN_TYPE=debug


ROOT=/home/larsonke
DATA_ROOT=/data/h_oguz_lab/larsonke
PROJECT_NAME=SA_resolutionTesting_Kirby

greedy=$ROOT/bin/greedy
slicer=/home/oguzi/bin/Slicer
holedetection=$ROOT/bin/LOGB-Tools/HoleDetection
fillholes=$ROOT/bin/LOGB-Tools/FillHoles
imagemath=$ROOT/bin/ImageMath
performErosion=$ROOT/Code/Synthetic-Atrophy/tools/performErosion
maskField=$ROOT/Code/Synthetic-Atrophy/tools/maskField
maskFieldOld=$ROOT/Code/Synthetic-Atrophy/tools/maskField_old
VTP2VTK=$ROOT/bin/VTK-Tools/VTP2VTK
VTK2VTP=$ROOT/bin/VTK-Tools/VTK2VTP
copyScalars=$ROOT/bin/VTK-Tools/copyScalars
warpSurfaces_GM=$ROOT/Code/Synthetic-Atrophy/tools/warpSurfs_noCheck
warpSurfaces_WM=$ROOT/Code/Synthetic-Atrophy/tools/warpSurfs_checkMask
parcellateAtrophyRegion=$ROOT/Code/Synthetic-Atrophy/tools/parcellateAtrophyRegion
FSThickness=$ROOT/bin/Cortical-Thickness-Measurement/FSThickness
readThickness=$ROOT/bin/Cortical-Thickness-Measurement/readThickness
getBlurMask=$ROOT/Code/Synthetic-Atrophy/tools/getBlurMask
getTrimRegion=$ROOT/Code/Synthetic-Atrophy/tools/getTrimRegion
trimImage=$ROOT/Code/Synthetic-Atrophy/tools/trimImage
padImage=$ROOT/Code/Synthetic-Atrophy/tools/padImage
padVectorImage=$ROOT/Code/Synthetic-Atrophy/tools/padVectorImage
trimFloatImage=$ROOT/Code/Synthetic-Atrophy/tools/trimFloatImage
angle_normal2warp=$ROOT/Code/Synthetic-Atrophy/tools/normal2warp
addFloatArray=$ROOT/bin/VTK-Tools/parcellateSingleSurface
readMeshArray=$ROOT/bin/VTK-Tools/readMeshArray
parcellateROImesh=$ROOT/Code/Synthetic-Atrophy/tools/parcellateROIMesh
distanceBetweenSurfaces=$ROOT/Code/Synthetic-Atrophy/tools/distanceBetweenSurfaces



INITIAL_DATASET=Kirby
SUBJECT_FILE=$DATA_ROOT/Raw/$INITIAL_DATASET/all-ids.txt
SUBJECT_IDS=( $( awk -F " " '{print $1}' $SUBJECT_FILE ) )
NSUBJECTS=${#SUBJECT_IDS[@]}
NSUBJECTS=$(($NSUBJECTS-1))

LABEL_FILE=$ROOT/Code/Data-Processing/$PROJECT_NAME/labels_DKTatlas_rh.txt
LABELS_DKT=( $( awk -F " " '{print $1}' $LABEL_FILE ) )



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
    c3d -interpolation NearestNeighbor $trim -resample 400% -holefill 1 0 -o $trim_highRes
    cropRegion_origin=$(sed -n 2p $cropRegion_file)
}



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



function Get-True-Thickness-Original(){
    SUBJECT=$1
    label=$2

    RAW_DIR=$DATA_ROOT/Raw/$PROJECT_NAME/$SUBJECT
    WDIR=$DATA_ROOT/Processed/$PROJECT_NAME/$SUBJECT
    THICKNESS_DIR=$DATA_ROOT/Results/$PROJECT_NAME/$SUBJECT/$label


    # Make label maps
    mkdir -p $WDIR/$label/dump
    label_original=$WDIR/$label/ROI/label_original.nii.gz
    ribbon=$RAW_DIR/Extras/ribbon_native.nii.gz
    aparc_aseg=$RAW_DIR/Extras/aparc+aseg_native.nii.gz
    wholeBrainImage_noRegion=$WDIR/$label/dump/wholeBrain_noRegion.nii.gz
    wholeBrainImage_original_parc=$WDIR/$label/wholeBrain_original_withLabel.nii.gz   
    
    if false ; then
	c3d $aparc_aseg -threshold 1 inf 1 0 $label_original -scale -1 -add -o $wholeBrainImage_noRegion
        c3d $label_original -dilate 1 4x4x4vox $wholeBrainImage_noRegion -times -o $WDIR/$label/dump/label_surround.nii.gz
        c3d $wholeBrainImage_noRegion $WDIR/$label/dump/label_surround.nii.gz -add $label_original -scale 3 -add -o $wholeBrainImage_original_parc
    fi
    rm -r  $WDIR/$label/dump

    # Apply label map
    GMmesh_original_lh=$WDIR/Original/GM_lh_0.vtp
    GMmesh_original_rh=$WDIR/Original/GM_rh_0.vtp
    WMmesh_original_lh=$WDIR/Original/WM_lh_0.vtp
    WMmesh_original_rh=$WDIR/Original/WM_rh_0.vtp

    PARC_DIR=$WDIR/$label/Surfaces
    mkdir -p $PARC_DIR
    GMmesh_original_parc_lh=$PARC_DIR/GM_lh_0_parc.vtp
    WMmesh_original_parc_lh=$PARC_DIR/WM_lh_0_parc.vtp
    GMmesh_original_parc_rh=$PARC_DIR/GM_rh_0_parc.vtp
    WMmesh_original_parc_rh=$PARC_DIR/WM_rh_0_parc.vtp

    if false ; then
        $parcellateAtrophyRegion $WMmesh_original_lh $GMmesh_original_lh $wholeBrainImage_original_parc $WMmesh_original_parc_lh $GMmesh_original_parc_lh 0
	$parcellateAtrophyRegion $WMmesh_original_rh $GMmesh_original_rh $wholeBrainImage_original_parc $WMmesh_original_parc_rh $GMmesh_original_parc_rh 0
    fi
    

    atrophyRegion_lh_pial_txt=$THICKNESS_DIR/lh_pial_atrophyRegion.txt
    atrophyRegion_lh_white_txt=$THICKNESS_DIR/lh_white_atrophyRegion.txt
    atrophyRegion_rh_pial_txt=$THICKNESS_DIR/rh_pial_atrophyRegion.txt
    atrophyRegion_rh_white_txt=$THICKNESS_DIR/rh_white_atrophyRegion.txt

    if true ; then
        if [ ! -f $atrophyRegion_lh_pial_txt ] ; then
	    $readMeshArray $GMmesh_original_parc_lh $atrophyRegion_lh_pial_txt Atrophy-Region
	fi
	if [ ! -f $atrophyRegion_lh_white_txt ] ; then
            $readMeshArray $WMmesh_original_parc_lh $atrophyRegion_lh_white_txt Atrophy-Region
	fi
	if [ ! -f $atrophyRegion_rh_pial_txt ] ; then
            $readMeshArray $GMmesh_original_parc_rh $atrophyRegion_rh_pial_txt Atrophy-Region
	fi
	if [ ! -f $atrophyRegion_rh_white_txt ] ; then
            $readMeshArray $WMmesh_original_parc_rh $atrophyRegion_rh_white_txt Atrophy-Region
	fi
    fi



    # Calculate thickness
    mkdir -p $THICKNESS_DIR/surf
    thickness_original_lh=$THICKNESS_DIR/surf/lh_thickness-0.vtp 
    thickness_original_rh=$THICKNESS_DIR/surf/rh_thickness-0.vtp

    if false ; then
        $FSThickness $WMmesh_original_parc_lh $GMmesh_original_parc_lh $thickness_original_lh
        $readThickness $thickness_original_lh $THICKNESS_DIR/lh_0 Atrophy-Region
        $FSThickness $WMmesh_original_parc_rh $GMmesh_original_parc_rh $thickness_original_rh
        $readThickness $thickness_original_rh $THICKNESS_DIR/rh_0 Atrophy-Region
    fi


}




function Get-True-Thickness-In-Region(){
    SUBJECT=$1
    label=$2
    kernel=$3
    
    RAW_DIR=$DATA_ROOT/Raw/$PROJECT_NAME/$SUBJECT
    if [ $kernel -lt 10 ] ; then
	K_ID=0$kernel
    else
	K_ID=$kernel
    fi
    WDIR=$DATA_ROOT/Processed/$PROJECT_NAME/$SUBJECT
    THICKNESS_DIR=$DATA_ROOT/Results/$PROJECT_NAME/$SUBJECT/$label/${K_ID}/Surface-Distance

    # Perform erosion on ROI
    ROI_DIR=$WDIR/$label/ROI
    MASK_DIR=$WDIR/$label/HighRes-Mask
    aparc_aseg_trim_highRes=$MASK_DIR/aparc+aseg_native_trim_highRes.nii.gz
    wm_trim_highRes=$MASK_DIR/wm_trim_highRes.nii.gz
    atrophyROI_trim_highRes=$ROI_DIR/label_atrophy_${K_ID}_trim_highRes.nii.gz
    
    cropRegion_file=$MASK_DIR/crop_region_info.txt
    pad_info=$(sed -n 1p $cropRegion_file)
    cropRegion_origin=$(sed -n 2p $cropRegion_file)

    if false ; then
	if [ $kernel == 1 ] ; then
	    originalROI_trim_highRes=$ROI_DIR/label_original_trim_highRes.nii.gz
	    c3d $originalROI_trim_highRes -holefill 1 0 -o $originalROI_trim_highRes
	    c3d $aparc_aseg_trim_highRes -threshold 1 inf 1 0 -o $ROI_DIR/label_original_joint_trim_highRes.nii.gz
	elif [ $kernel -ge 2 ] && [ $kernel -le 9 ] ; then
            let kernel_last=$kernel-1
            originalROI_trim_highRes=$ROI_DIR/label_atrophy_0${kernel_last}_trim_highRes.nii.gz
	elif [ $kernel == 10 ] ; then
	    originalROI_trim_highRes=$ROI_DIR/label_atrophy_09_trim_highRes.nii.gz
	else
	    let kernel_last=$kernel-1 
	    originalROI_trim_highRes=$ROI_DIR/label_atrophy_${kernel_last}_trim_highRes.nii.gz
	fi
    fi

    joint_label_original=$ROI_DIR/label_original_joint_trim_highRes.nii.gz
    joint_label_atrophy=$ROI_DIR/label_atrophy_joint_${K_ID}_trim_highRes.nii.gz

    if false ; then
	$performErosion $originalROI_trim_highRes $wm_trim_highRes $atrophyROI_trim_highRes 1 0

	let lower_threshold=$label-1
	let upper_threshold=$label+1
	c3d $aparc_aseg_trim_highRes -threshold 1 $lower_threshold 1 0 -popas LM $aparc_aseg_trim_highRes -threshold $upper_threshold inf 1 0 -popas UM -push LM -push UM -add $atrophyROI_trim_highRes -add -o $joint_label_atrophy
    fi

    
    # Get the transform
    TRANSFORM_DIR=$WDIR/$label/Transforms/$K_ID
    mkdir -p $TRANSFORM_DIR $TRANSFORM_DIR/Dump
    originalROI_trim_highRes=$ROI_DIR/label_original_trim_highRes.nii.gz
    field_outbase=$TRANSFORM_DIR/field_${K_ID}
    field_inverse_outbase=$TRANSFORM_DIR/field_inverse_${K_ID}

    if false ; then
	if [ ! -f ${field_inverse_outbase}_full_trim_highRes.nii.gz ] ; then
	    $greedy -d 3 -m MSQ -i $joint_label_atrophy $joint_label_original -o ${field_outbase}_full_trim_highRes.nii.gz -oinv ${field_inverse_outbase}_full_trim_highRes.nii.gz -mm $joint_label_original -n 100x100x50x10 >> /dev/null
	fi
	
	if false ; then
	    if [ ! -f ${field_outbase}_maskBlur_trim.nii.gz ] ; then
		c3d $MASK_DIR/brain_std_trim.nii.gz -origin 0x0x0mm -o $TRANSFORM_DIR/Dump/brain_std_trim.nii.gz
		c3d $ROI_DIR/label_original_trim.nii.gz -origin 0x0x0mm -o $TRANSFORM_DIR/Dump/label_original_trim.nii.gz
		c3d $MASK_DIR/brain_mask_trim.nii.gz -origin 0x0x0mm -o $TRANSFORM_DIR/Dump/brain_mask_trim.nii.gz
		c3d $MASK_DIR/wm_trim.nii.gz -origin 0x0x0mm -o $TRANSFORM_DIR/Dump/wm_trim.nii.gz
		c3d -mcs ${field_inverse}_full_trim_highRes.nii.gz -foreach -resample 25% -origin 0x0x0mm -endfor -omc ${field_inverse}_full_trim.nii.gz
		
		$maskField $TRANSFORM_DIR/Dump/brain_std_trim.nii.gz $TRANSFORM_DIR/Dump/label_original_trim.nii.gz $TRANSFORM_DIR/Dump/brain_mask_trim.nii.gz ${field_outbase}_full_trim.nii.gz $TRANSFORM_DIR/Dump/wm_trim.nii.gz 0.2 $TRANSFORM_DIR/Dump/points.vtp ${field_outbase}_maskBlur_trim.nii.gz
		
		$padVectorImage $MASK_DIR/total_mask.nii.gz ${field_outbase}_maskBlur_trim.nii.gz ${field_outbase}_maskBlur_pad.nii.gz $pad_info
		c3d -mcs ${field_outbase}_maskBlur_pad.nii.gz -foreach -origin 0x0x0vox -endfor -omc ${field_outbase}_maskBlur_pad.nii.gz
	    fi
	fi

	if true ; then
	    if [ ! -f ${field_inverse_outbase}_maskBlur_trim.nii.gz ] ; then
		c3d $MASK_DIR/brain_std_trim.nii.gz -origin 0x0x0mm -o $TRANSFORM_DIR/Dump/brain_std_trim.nii.gz
		c3d $ROI_DIR/label_original_trim.nii.gz -origin 0x0x0mm -o $TRANSFORM_DIR/Dump/label_original_trim.nii.gz
		c3d $MASK_DIR/brain_mask_trim.nii.gz -origin 0x0x0mm -o $TRANSFORM_DIR/Dump/brain_mask_trim.nii.gz
		c3d $MASK_DIR/wm_trim.nii.gz -origin 0x0x0mm -o $TRANSFORM_DIR/Dump/wm_trim.nii.gz
		c3d -mcs ${field_inverse_outbase}_full_trim_highRes.nii.gz -foreach -resample 25% -origin 0x0x0mm -endfor -omc ${field_inverse_outbase}_full_trim.nii.gz
		
		$maskField $TRANSFORM_DIR/Dump/brain_std_trim.nii.gz $TRANSFORM_DIR/Dump/label_original_trim.nii.gz $TRANSFORM_DIR/Dump/brain_mask_trim.nii.gz ${field_inverse_outbase}_full_trim.nii.gz $TRANSFORM_DIR/Dump/wm_trim.nii.gz 0.2 $TRANSFORM_DIR/Dump/points.vtp ${field_inverse_outbase}_maskBlur_trim.nii.gz

		$padVectorImage $MASK_DIR/total_mask.nii.gz ${field_inverse_outbase}_maskBlur_trim.nii.gz ${field_inverse_outbase}_maskBlur_pad.nii.gz $pad_info
		c3d -mcs ${field_inverse_outbase}_maskBlur_pad.nii.gz -foreach -origin 0x0x0vox -endfor -omc ${field_inverse_outbase}_maskBlur_pad.nii.gz
	    fi
	fi
	rm -r $TRANSFORM_DIR/Dump
    fi


    # Apply warp to images to get dataset	
    mkdir -p $RAW_DIR/$label
    label_original_pad=$ROI_DIR/label_original_pad.nii.gz
    label_atrophy_pad=$ROI_DIR/label_atrophy_${K_ID}_pad.nii.gz
    MPRAGE_original=$RAW_DIR/$SUBJECT-0-MPRAGE.nii.gz
    FLAIR_original=$RAW_DIR/$SUBJECT-0-FLAIR.nii.gz
    MPRAGE_atrophy=$RAW_DIR/$label/$SUBJECT-1-$K_ID-MPRAGE.nii.gz
    FLAIR_atrophy=$RAW_DIR/$label/$SUBJECT-1-$K_ID-FLAIR.nii.gz

    if false ; then
        $greedy -d 3 -rf $MPRAGE_original -rm $MPRAGE_original $MPRAGE_atrophy -r ${field_outbase}_maskBlur_pad.nii.gz >> /dev/null
	$slicer --launch OrientScalarVolume --orientation LAI $MPRAGE_original $MPRAGE_original >> /dev/null
	$slicer --launch OrientScalarVolume --orientation LAI $MPRAGE_atrophy $MPRAGE_atrophy >> /dev/null

	$greedy -d 3 -rf $MPRAGE_original -rm $FLAIR_original $FLAIR_atrophy -r ${field_outbase}_maskBlur_pad.nii.gz >> /dev/null
	$slicer --launch OrientScalarVolume --orientation LAI $FLAIR_original $FLAIR_original >> /dev/null
	$slicer --launch OrientScalarVolume --orientation LAI $FLAIR_atrophy $FLAIR_atrophy >> /dev/null
    fi
    
    

    # Get final surfaces for thickness
    PARC_DIR=$WDIR/$label/Surfaces
    label_original_mesh=$PARC_DIR/label_original_mesh.vtp
    
    GMmesh_original_parc_lh=$PARC_DIR/GM_lh_0_parc.vtp
    GMmesh_original_parc_rh=$PARC_DIR/GM_rh_0_parc.vtp
    WMmesh_original_parc_lh=$PARC_DIR/WM_lh_0_parc.vtp
    WMmesh_original_parc_rh=$PARC_DIR/WM_rh_0_parc.vtp
    GMmesh_atrophy_parc_lh=$PARC_DIR/GM_lh_1_${K_ID}_parc.vtp
    WMmesh_atrophy_parc_lh=$PARC_DIR/WM_lh_1_${K_ID}_parc.vtp
    GMmesh_atrophy_parc_rh=$PARC_DIR/GM_rh_1_${K_ID}_parc.vtp
    WMmesh_atrophy_parc_rh=$PARC_DIR/WM_rh_1_${K_ID}_parc.vtp

    if true ; then
	if [ -f ${field_inverse_outbase}_maskBlur_pad.nii.gz ] ; then
	    if [ ! -f $GMmesh_atrophy_parc_lh ] ; then
		$warpSurfaces_GM $GMmesh_original_parc_lh ${field_inverse_outbase}_maskBlur_pad.nii.gz $MASK_DIR/total_mask.nii.gz $GMmesh_atrophy_parc_lh 0
	    fi
	    if [ ! -f $GMmesh_atrophy_parc_rh ] ; then
		$warpSurfaces_GM $GMmesh_original_parc_rh ${field_inverse_outbase}_maskBlur_pad.nii.gz $MASK_DIR/total_mask.nii.gz $GMmesh_atrophy_parc_rh 0
	    fi
	    if [ ! -f $WMmesh_atrophy_parc_lh ] ; then
		$warpSurfaces_WM $WMmesh_original_parc_lh ${field_inverse_outbase}_maskBlur_pad.nii.gz $MASK_DIR/total_mask.nii.gz $WMmesh_atrophy_parc_lh 0
            fi
            if [ ! -f $WMmesh_atrophy_parc_rh ] ; then
		$warpSurfaces_WM $WMmesh_original_parc_rh ${field_inverse_outbase}_maskBlur_pad.nii.gz $MASK_DIR/total_mask.nii.gz $WMmesh_atrophy_parc_rh 0
            fi
	fi
    fi
    
    
    
    # Extract field components at each surface
    mkdir -p $THICKNESS_DIR
    array_name="Surface-Distance"
    
    surfWithTrueChange_original_lh_pial=$THICKNESS_DIR/lh_pial_trueChange.vtp
    surfWithTrueChange_original_lh_white=$THICKNESS_DIR/lh_white_trueChange.vtp
    surfWithTrueChange_original_rh_pial=$THICKNESS_DIR/rh_pial_trueChange.vtp
    surfWithTrueChange_original_rh_white=$THICKNESS_DIR/rh_white_trueChange.vtp

    surfWithTrueChange_atrophy_lh_pial=$THICKNESS_DIR/lh_pial_trueChange1.vtp
    surfWithTrueChange_atrophy_lh_white=$THICKNESS_DIR/lh_white_trueChange1.vtp
    surfWithTrueChange_atrophy_rh_pial=$THICKNESS_DIR/rh_pial_trueChange1.vtp
    surfWithTrueChange_atrophy_rh_white=$THICKNESS_DIR/rh_white_trueChange1.vtp

    surfWithTrueChange_original_lh_pial_txt=$THICKNESS_DIR/lh_pial_trueChange.txt
    surfWithTrueChange_original_lh_white_txt=$THICKNESS_DIR/lh_white_trueChange.txt
    surfWithTrueChange_original_rh_pial_txt=$THICKNESS_DIR/rh_pial_trueChange.txt
    surfWithTrueChange_original_rh_white_txt=$THICKNESS_DIR/rh_white_trueChange.txt

    surfWithTrueChange_atrophy_lh_pial_txt=$THICKNESS_DIR/lh_pial_trueChange1.txt
    surfWithTrueChange_atrophy_lh_white_txt=$THICKNESS_DIR/lh_white_trueChange1.txt
    surfWithTrueChange_atrophy_rh_pial_txt=$THICKNESS_DIR/rh_pial_trueChange1.txt
    surfWithTrueChange_atrophy_rh_white_txt=$THICKNESS_DIR/rh_white_trueChange1.txt

    # Measure distance between surfaces
    if true ; then
	#$distanceBetweenSurfaces $GMmesh_original_parc_lh $GMmesh_atrophy_parc_lh $surfWithTrueChange_original_lh_pial
	$distanceBetweenSurfaces $GMmesh_atrophy_parc_lh $GMmesh_original_parc_lh $surfWithTrueChange_atrophy_lh_pial
	
	#$distanceBetweenSurfaces $GMmesh_original_parc_rh $GMmesh_atrophy_parc_rh $surfWithTrueChange_original_rh_pial
	$distanceBetweenSurfaces $GMmesh_atrophy_parc_rh $GMmesh_original_parc_rh $surfWithTrueChange_atrophy_rh_pial

	#$distanceBetweenSurfaces $WMmesh_original_parc_lh $WMmesh_atrophy_parc_lh $surfWithTrueChange_original_lh_white
	$distanceBetweenSurfaces $WMmesh_atrophy_parc_lh $WMmesh_original_parc_lh $surfWithTrueChange_atrophy_lh_white

	#$distanceBetweenSurfaces $WMmesh_original_parc_rh $WMmesh_atrophy_parc_rh $surfWithTrueChange_original_rh_white
	$distanceBetweenSurfaces $WMmesh_atrophy_parc_rh $WMmesh_original_parc_rh $surfWithTrueChange_atrophy_rh_white
    fi

    if true ; then
	#$readMeshArray $surfWithTrueChange_original_lh_pial $surfWithTrueChange_original_lh_pial_txt $array_name
	$readMeshArray $surfWithTrueChange_atrophy_lh_pial $surfWithTrueChange_atrophy_lh_pial_txt $array_name
	
	#$readMeshArray $surfWithTrueChange_original_lh_white $surfWithTrueChange_original_lh_white_txt $array_name
	$readMeshArray $surfWithTrueChange_atrophy_lh_white $surfWithTrueChange_atrophy_lh_white_txt $array_name

	#$readMeshArray $surfWithTrueChange_original_rh_pial $surfWithTrueChange_original_rh_pial_txt $array_name
	$readMeshArray $surfWithTrueChange_atrophy_rh_pial $surfWithTrueChange_atrophy_rh_pial_txt $array_name
	
	#$readMeshArray $surfWithTrueChange_original_rh_white $surfWithTrueChange_original_rh_white_txt $array_name
	$readMeshArray $surfWithTrueChange_atrophy_rh_white $surfWithTrueChange_atrophy_rh_white_txt $array_name
    fi

    if true ; then
	rm $THICKNESS_DIR/*.vtp
    fi
}




function Copy-Initial-Files(){

    SA_SUBJECT=$1
    INITIAL_SUBJECT=$2

    
    FS_DIR=$DATA_ROOT/Processed/$INITIAL_DATASET/FreeSurfer/$INITIAL_SUBJECT/$INITIAL_SUBJECT-01
    RAW_DIR=$DATA_ROOT/Raw/$PROJECT_NAME/$SA_SUBJECT
    mkdir -p $RAW_DIR $RAW_DIR/Extras
    
    # Copy aseg file
    ribbon_FS=$FS_DIR/mri/ribbon_native.nii.gz
    ribbon=$RAW_DIR/Extras/ribbon_native.nii.gz
    labels_FS=$FS_DIR/mri/aparc+aseg_native.nii.gz
    aparc_aseg=$RAW_DIR/Extras/aparc+aseg_native.nii.gz
    brain_mask_FS=$FS_DIR/mri/brainmask_native.nii.gz
    brain_mask=$RAW_DIR/Extras/skullstrip_mask.nii.gz
    wm_FS=$FS_DIR/mri/wm_native.nii.gz
    wm=$RAW_DIR/Extras/wm.nii.gz
    brain_std=$RAW_DIR/Extras/brain_std.nii.gz

    if false ; then
	$slicer --launch OrientScalarVolume --orientation LPI $ribbon_FS $ribbon >> /dev/null
	c3d $ribbon -origin-voxel 0x0x0vox -o $ribbon 
	$slicer --launch OrientScalarVolume --orientation LPI $labels_FS $aparc_aseg >> /dev/null
	c3d $aparc_aseg -origin-voxel 0x0x0vox -o $aparc_aseg
	$slicer --launch OrientScalarVolume --orientation LPI $brain_mask_FS $brain_mask >> /dev/null
	c3d $brain_mask -thresh 1 inf 1 0 -origin-voxel 0x0x0vox -o $brain_mask
	$slicer --launch OrientScalarVolume --orientation LPI $wm_FS $wm >> /dev/null
	c3d $wm -thresh 1 inf 1 0 -origin-voxel 0x0x0vox -o $wm
        c3d $aparc_aseg -threshold 1 inf 1 0 -sdt -o $brain_std
    fi

 
    # Register FLAIR to T1
    MPRAGE_root=$DATA_ROOT/Raw/Kirby/$INITIAL_SUBJECT/$INITIAL_SUBJECT-01-MPRAGE.nii.gz
    FLAIR_root=$DATA_ROOT/Raw/Kirby/$INITIAL_SUBJECT/$INITIAL_SUBJECT-01-FLAIR.nii.gz 

    moments=$RAW_DIR/Extras/FLAIR_to_MPRAGE_0_moments.mat
    rigid=$RAW_DIR/Extras/FLAIR_to_MPRAGE_0_rigid.mat
    FLAIR_registered=$RAW_DIR/Extras/$SA_SUBJECT-0-FLAIR_T1space.nii.gz

    if false ; then
	$greedy -d 3 -i $MPRAGE_root $FLAIR_root -o $moments -moments 1 >> /dev/null
	$greedy -d 3 -i $MPRAGE_root $FLAIR_root -o $rigid -a -dof 6 -n 100x50x10 -m NCC 4x4x4 -ia $moments >> /dev/null
	$greedy -d 3 -r $rigid -rf $MPRAGE_root -rm $FLAIR_root $FLAIR_registered >> /dev/null
    fi

    
    # Re-orient and re-sample so I can get the masks
    MPRAGE=$RAW_DIR/$SA_SUBJECT-0-MPRAGE.nii.gz
    FLAIR=$RAW_DIR/$SA_SUBJECT-0-FLAIR.nii.gz
    if false ; then
	$slicer --launch OrientScalarVolume --orientation LPI $MPRAGE_root $MPRAGE >> /dev/null
	c3d $MPRAGE -origin-voxel 0x0x0vox -o $MPRAGE
	$slicer --launch OrientScalarVolume --orientation LPI $FLAIR_registered $FLAIR >> /dev/null
	c3d $FLAIR -origin-voxel 0x0x0vox -o $FLAIR
    fi

    # Get the original meshes
    WDIR=$DATA_ROOT/Processed/$PROJECT_NAME/$SUBJECT/Original
    mkdir -p $WDIR
    WMmask_original_lh=$WDIR/WM_lh_0.nii.gz
    WMmask_original_rh=$WDIR/WM_rh_0.nii.gz    
    wholeBrainImage_original_lh=$WDIR/GM_lh_0.nii.gz
    wholeBrainImage_original_rh=$WDIR/GM_rh_0.nii.gz

    if false ; then
	c3d $ribbon -threshold 2 3 1 0 -as LM $aparc_aseg -threshold 1 inf 1 0 -times -o $wholeBrainImage_original_lh -pop -push LM $wm -times -o $WMmask_original_lh
        c3d $ribbon -threshold 41 42 1 0 -as RM $aparc_aseg -threshold 1 inf 1 0 -times -o $wholeBrainImage_original_rh -pop -push RM $wm -times -o $WMmask_original_rh	
    fi

    # Get final surfaces for thickness
    GMmesh_original_lh=$WDIR/GM_lh_0.vtp
    GMmesh_original_rh=$WDIR/GM_rh_0.vtp
    WMmesh_original_lh=$WDIR/WM_lh_0.vtp
    WMmesh_original_rh=$WDIR/WM_rh_0.vtp

    if false ; then
	Get-Mesh $wholeBrainImage_original_lh $GMmesh_original_lh $WDIR
        Get-Mesh $wholeBrainImage_original_rh $GMmesh_original_rh $WDIR
        Get-Mesh $WMmask_original_lh $WMmesh_original_lh $WDIR
        Get-Mesh $WMmask_original_rh $WMmesh_original_rh $WDIR
    fi

}





function run-SA-subject(){
    if [ $RUN_TYPE == "accre" ] ; then
        S=$SLURM_ARRAY_TASK_ID
    elif [ $RUN_TYPE == "debug" ] ; then
	S=$1
    fi

    if [ $S -lt 10 ] ; then
	SUBJECT=0$S
    else
        SUBJECT=$S
    fi
    blurKernel=1
    
    KIRBY_SUBJECT=${SUBJECT_IDS[$S]}
    #Copy-Initial-Files $SUBJECT $KIRBY_SUBJECT
    NLABELS=${#LABELS_DKT[@]}

    #n=$(ls $DATA_ROOT/Results/$PROJECT_NAME/$SUBJECT -1 | wc -l)
    #let n=$n-1
    n=0
    if true ; then
	while [ $n -lt $NLABELS ] ; do
	    LABEL=${LABELS_DKT[$n]}
	    echo "%~%~%~%~%~%~%~%" $SUBJECT $LABEL "%~%~%~%~%~%~%~%"
	    
	    #Get-HighRes-Mask $SUBJECT $LABEL
	    Get-True-Thickness-Original $SUBJECT $LABEL $kernel
	    
	    # Find true change
	    kernel=1
	    maxKernel=12
	    if true ; then
		while [ $kernel -le $maxKernel ] ; do 
		    echo Get-True-Thickness-In-Region $SUBJECT $LABEL $kernel
		    if [ $kernel -lt 10 ] ; then
			K_ID=0$kernel
		    else
			K_ID=$kernel
		    fi
		    fileRef=$DATA_ROOT/Results/$PROJECT_NAME/$SUBJECT/$LABEL/$K_ID/Surface-Distance/rh_white_trueChange1.txt
		    if [ -f $fileRef ] ; then
			nLines=$(cat $fileRef | wc -l)
			if [ $nLines -le 0 ] ; then
			    Get-True-Thickness-In-Region $SUBJECT $LABEL $kernel
			fi
		    else
			Get-True-Thickness-In-Region $SUBJECT $LABEL $kernel
		    fi
		    let kernel=$kernel+1

		done
	    fi
	    let n=n+1
	done
    fi
}



function main {
    NSUBJECTS=${#SUBJECT_IDS[@]}
    let NSUBJECTS1=$NSUBJECTS-1
    if [ $RUN_TYPE == "accre" ] ; then
	sbatch --array=0-$NSUBJECTS --output=$ROOT/Code/Data-Processing/$PROJECT_NAME/out-%a-SA_v6_rh $0 run-SA-subject
    elif [ $RUN_TYPE == "debug" ] ; then
	set -x
	ii=0
	while [ $ii -lt $NSUBJECTS ] ; do
            run-SA-subject $ii
            let ii=$ii+1
	    #exit
        done
    fi
}



if [[ $1 ]]; then
  command=$1
  echo $1
  shift
  $command $@
else
    main
fi

