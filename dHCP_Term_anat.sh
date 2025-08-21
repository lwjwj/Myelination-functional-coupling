###### Subject list and path
subjlist='/data/users/wliu/demo_dHCP_Analysis/TermList_myelin.txt'
for s in $(cat $subjlist)
do
if [[ $s == CC* ]] 
then
subj=$s
continue
else
sess=$s
fi

subjid=$subj
session=$sess
WorkbenchBinary='/opt/workbench/exe_rh_linux64/wb_command'  #path of Workbench
StandardSpherePath=/data/users/wliu/demo_dHCP_Analysis/wb_creat_xk
subPath=/data/users/wliu/dHCP_org_out
OutPath=/data/users/wliu/demo_dHCP_Analysis/
resolution=5
Population='TermAsymmetry'
mkdir -p ${OutPath}/1.5smooth/sub-${subjid}/ses-${session}  


for hemi in left right ; do

###### Surface downsampling
echo "surface downsampling"
for anatomical in midthickness pial; do 
SphereReg=/data/users/wliu/code/CorticalAsymmetry/Templates/dHCP/week-40_hemi-${hemi}_space-dhcpSym_dens-32k_sphere.surf.gii
TemplateSphere=${StandardSpherePath}/sphere.${resolution}k.${hemi}.surf.gii
NativeSurface=${subPath}/new-dHCPAsymmetry-struc-MSMself/${Population}/sub-${subjid}/ses-${session}/sub-${subjid}_ses-${session}_hemi-${hemi}_${anatomical}.MSMStrain.surf.gii
OutSurface=${OutPath}/sub-${subjid}/ses-${session}/sub-${subjid}_ses-${session}_hemi-${hemi}_${anatomical}.${resolution}k.surf.gii
${WorkbenchBinary} -surface-resample ${NativeSurface} ${SphereReg} ${TemplateSphere} BARYCENTRIC ${OutSurface} 
done

echo "Midthickness surface smoothing"
anatomical=midthickness
MidthicknessNosmoothed=${OutPath}/sub-${subjid}/ses-${session}/sub-${subjid}_ses-${session}_hemi-${hemi}_${anatomical}.${resolution}k.surf.gii
MidthicknessSmoothed=${OutPath}/sub-${subjid}/ses-${session}/sub-${subjid}_ses-${session}_hemi-${hemi}_${anatomical}.smoothed.${resolution}k.surf.gii
${WorkbenchBinary} -surface-smoothing ${MidthicknessNosmoothed} 0.75 10 ${MidthicknessSmoothed}

###### Myelination downsampling
AnalysisDirectory=/data/users/wliu/dHCP_org_out
DataDirectory=/data/users/wliu/dHCP_org
Metric=myelinmap
echo "myelinmap form native to 40w_32k"
SphereReg=${AnalysisDirectory}/new-surface_transforms/sub-${subjid}/ses-${session}/sub-${subjid}_ses-${session}_hemi-${hemi}_from-native_to-dhcpSym40_dens-32k_mode-sphere.reg40.surf.gii
TemplateSphere=/data/users/wliu/code/CorticalAsymmetry/Templates/dHCP/week-40_hemi-${hemi}_space-dhcpSym_dens-32k_sphere.surf.gii
NativeMidthickness=${DataDirectory}/sub-${subjid}/ses-${session}/anat/sub-${subjid}_ses-${session}_hemi-${hemi}_midthickness.surf.gii
TemplateMidthickness=${AnalysisDirectory}/new-dHCPAsymmetry-struc-MSMself/${Population}/sub-${subjid}/ses-${session}/sub-${subjid}_ses-${session}_hemi-${hemi}_midthickness.MSMStrain.surf.gii
NativeMetric=${DataDirectory}/sub-${subjid}/ses-${session}/anat/sub-${subjid}_ses-${session}_hemi-${hemi}_${Metric}.shape.gii
OutMetric=${AnalysisDirectory}/new-dHCPAsymmetry-struc-MSMself/${Population}/sub-${subjid}/ses-${session}/sub-${subjid}_ses-${session}_hemi-${hemi}_${Metric}.MSMStrain.shape.gii
${WorkbenchBinary} -metric-resample ${NativeMetric} ${SphereReg} ${TemplateSphere} ADAP_BARY_AREA ${OutMetric} -area-surfs ${NativeMidthickness} ${TemplateMidthickness} 

echo "myelinmap from 40w_32k to 5k" 
SphereReg=/data/users/wliu/code/CorticalAsymmetry/Templates/dHCP/week-40_hemi-${hemi}_space-dhcpSym_dens-32k_sphere.surf.gii
TemplateSphere=${StandardSpherePath}/sphere.${resolution}k.${hemi}.surf.gii
NativeMidthickness=${subPath}/new-dHCPAsymmetry-struc-MSMself/${Population}/sub-${subjid}/ses-${session}/sub-${subjid}_ses-${session}_hemi-${hemi}_midthickness.MSMStrain.surf.gii
TemplateMidthickness=${OutPath}/sub-${subjid}/ses-${session}/sub-${subjid}_ses-${session}_hemi-${hemi}_midthickness.${resolution}k.surf.gii
NativeMetric=${subPath}/new-dHCPAsymmetry-struc-MSMself/${Population}/sub-${subjid}/ses-${session}/sub-${subjid}_ses-${session}_hemi-${hemi}_${Metric}.MSMStrain.shape.gii
OutMetric=${OutPath}/sub-${subjid}/ses-${session}/sub-${subjid}_ses-${session}_hemi-${hemi}_${Metric}.${resolution}k.shape.gii
${WorkbenchBinary} -metric-resample ${NativeMetric} ${SphereReg} ${TemplateSphere} ADAP_BARY_AREA ${OutMetric} -area-surfs ${NativeMidthickness} ${TemplateMidthickness} 


###### Smooth the left and right myelinmap with corresponding midthickness surface
echo "myelination metric smoothing"
MidthicknessSmoothed=${OutPath}/sub-${subjid}/ses-${session}/sub-${subjid}_ses-${session}_hemi-${hemi}_midthickness.smoothed.5k.surf.gii
MetricIn=/data/users/wliu/demo_dHCP_Analysis/dHCP_Analysis_out/sub-${subjid}/ses-${session}/sub-${subjid}_ses-${session}_hemi-${hemi}_${Metric}.5k.shape.gii
SmoothingKernel=1.5
MetricOut=${OutPath}/sub-${subjid}/ses-${session}/sub-${subjid}_ses-${session}_hemi-${hemi}_${Metric}.5k.s15.shape.gii
${WorkbenchBinary} -metric-smoothing ${MidthicknessSmoothed} ${MetricIn} ${SmoothingKernel} ${MetricOut} 

done

###### Connect left myelinmap and right myelinmap
echo "myelinmap L and R"
MetricLeft=${OutPath}/sub-${subjid}/ses-${session}/sub-${subjid}_ses-${session}_hemi-left_${Metric}.5k.s15.shape.gii
MetricRight=${OutPath}/sub-${subjid}/ses-${session}/sub-${subjid}_ses-${session}_hemi-right_${Metric}.5k.s15.shape.gii
MetricLR=${OutPath}/1.5smooth/term_myelin_364/sub-${subjid}_ses-${session}_hemi-LR_${Metric}.hemi_5k.dtseries.nii
${WorkbenchBinary} -cifti-create-dense-timeseries ${MetricLR} -left-metric ${MetricLeft} -right-metric ${MetricRight}

###### Cross-sub merge
MetricMerge=${OutPath}/1.5smooth/term_myelin_364/term364_hemi-LR_${Metric}.hemi_5k.dtseries.nii
if [ -f ${MetricMerge} ]; then  
    ${WorkbenchBinary} -cifti-merge ${MetricMerge} -cifti ${MetricMerge} -cifti ${OutPath}/1.5smooth/term_myelin_364/sub-${subjid}_ses-${session}_hemi-LR_${Metric}.hemi_5k.dtseries.nii
else
    echo "term440 not exists! ${subjid} ${session}" ###Be sure to delete the previous 'term440_metric_asymmetry.shape.gii' before running again###
    cp ${OutPath}/1.5smooth/term_myelin_364/sub-${subjid}_ses-${session}_hemi-LR_${Metric}.hemi_5k.dtseries.nii  ${MetricMerge}
fi

if [ "${subjid}" == "CC01236XX16" ] && [ "${session}" == "155830" ]; then
    echo "Calculating ${Metric} cross-sub correlation"
    ${WorkbenchBinary} -cifti-correlation ${MetricMerge} ${OutPath}/1.5smooth/term_myelin_364/term364_hemi-LR_${Metric}.hemi_5k.dconn.nii
fi 

done
