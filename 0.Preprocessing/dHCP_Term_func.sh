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
VolumePath=/data/data/dHCP/dHCP/rel3data_fMRI/rel3_dhcp_fmri_pipeline
NativePath=/data/users/wliu/dHCP_org
OutPath=/data/users/wliu/demo_dHCP_Analysis
TemplatesDirectory=/data/users/wliu/code/CorticalAsymmetry/Templates/dHCP
AnalysisDirectory=/data/users/wliu/dHCP_org_out
StandardSpherePath=/data/users/wliu/demo_dHCP_Analysis/wb_creat_xk
Population=TermAsymmetry
RegName=MSMStrain
mkdir -p ${OutPath}/1.5smooth/sub-${subjid}/ses-${session}  


###### Volume-mapping-surface
echo "volume-mapping-surface"
for hemi in left right ; do
InVolume=${VolumePath}/sub-${subjid}/ses-${session}/func/bandpass_nomask2_gsr_f2a215mm.nii.gz
MidthicknessSurface=${NativePath}/sub-${subjid}/ses-${session}/anat/sub-${subjid}_ses-${session}_hemi-${hemi}_midthickness.surf.gii
OutBOLD=${OutPath}/sub-${subjid}/ses-${session}/sub-${subjid}_ses-${session}_hemi-${hemi}_BOLD_map-midthickness.native.func.gii
${WorkbenchBinary} -volume-to-surface-mapping  ${InVolume} ${MidthicknessSurface} ${OutBOLD} -trilinear

###### Registration and Resampling to 40week32k
echo "native to 40week32k"
TemplateSphere=${TemplatesDirectory}/week-40_hemi-${hemi}_space-dhcpSym_dens-32k_sphere.surf.gii
SphereReg=${AnalysisDirectory}/new-surface_transforms/sub-${subjid}/ses-${session}/sub-${subjid}_ses-${session}_hemi-${hemi}_from-native_to-dhcpSym40_dens-32k_mode-sphere.reg40.surf.gii
NativeMidthickness=${NativePath}/sub-${subjid}/ses-${session}/anat/sub-${subjid}_ses-${session}_hemi-${hemi}_midthickness.surf.gii
TemplateMidthickness=${AnalysisDirectory}/new-dHCPAsymmetry-struc-MSMself/${Population}/sub-${subjid}/ses-${session}/sub-${subjid}_ses-${session}_hemi-${hemi}_midthickness.${RegName}.surf.gii
NativeMetric=${OutPath}/sub-${subjid}/ses-${session}/sub-${subjid}_ses-${session}_hemi-${hemi}_BOLD_map-midthickness.native.func.gii
OutMetric=${OutPath}/sub-${subjid}/ses-${session}/sub-${subjid}_ses-${session}_hemi-${hemi}_BOLD_map-midthickness.${RegName}.func.gii
${WorkbenchBinary} -metric-resample ${NativeMetric} ${SphereReg} ${TemplateSphere} ADAP_BARY_AREA ${OutMetric} -area-surfs ${NativeMidthickness} ${TemplateMidthickness} 

###### Downsampling to 5k
echo "downsampling to 5k"
SphereReg=${TemplatesDirectory}/week-40_hemi-${hemi}_space-dhcpSym_dens-32k_sphere.surf.gii
TemplateSphere=${StandardSpherePath}/sphere.5k.${hemi}.surf.gii
NativeSurface=${AnalysisDirectory}/new-dHCPAsymmetry-struc-MSMself/${Population}/sub-${subjid}/ses-${session}/sub-${subjid}_ses-${session}_hemi-${hemi}_midthickness.${RegName}.surf.gii
OutSurface=${OutPath}/sub-${subjid}/ses-${session}/sub-${subjid}_ses-${session}_hemi-${hemi}_midthickness.5k.surf.gii
${WorkbenchBinary} -surface-resample ${NativeSurface} ${SphereReg} ${TemplateSphere} BARYCENTRIC ${OutSurface}

SphereReg=${TemplatesDirectory}/week-40_hemi-${hemi}_space-dhcpSym_dens-32k_sphere.surf.gii
TemplateSphere=${StandardSpherePath}/sphere.5k.${hemi}.surf.gii
NativeMidthickness=${AnalysisDirectory}/new-dHCPAsymmetry-struc-MSMself/${Population}/sub-${subjid}/ses-${session}/sub-${subjid}_ses-${session}_hemi-${hemi}_midthickness.${RegName}.surf.gii
TemplateMidthickness=${OutPath}/sub-${subjid}/ses-${session}/sub-${subjid}_ses-${session}_hemi-${hemi}_midthickness.5k.surf.gii
NativeMetric=${OutPath}/sub-${subjid}/ses-${session}/sub-${subjid}_ses-${session}_hemi-${hemi}_BOLD_map-midthickness.${RegName}.func.gii
OutMetric=${OutPath}/sub-${subjid}/ses-${session}/sub-${subjid}_ses-${session}_hemi-${hemi}_BOLD_map-midthickness.5k.func.gii
${WorkbenchBinary} -metric-resample ${NativeMetric} ${SphereReg} ${TemplateSphere} ADAP_BARY_AREA ${OutMetric} -area-surfs ${NativeMidthickness} ${TemplateMidthickness} 

###### Surface smoothing
echo "smoothing"
Midthickness=${OutPath}/sub-${subjid}/ses-${session}/sub-${subjid}_ses-${session}_hemi-${hemi}_midthickness.5k.surf.gii
MidthicknessSmoothed=${OutPath}/sub-${subjid}/ses-${session}/sub-${subjid}_ses-${session}_hemi-${hemi}_midthickness.smoothed.5k.surf.gii
${WorkbenchBinary} -surface-smoothing ${Midthickness} 0.75 10 ${MidthicknessSmoothed}

###### Metric smoothing
MidthicknessSmoothed=${OutPath}/sub-${subjid}/ses-${session}/sub-${subjid}_ses-${session}_hemi-${hemi}_midthickness.smoothed.5k.surf.gii
MetricIn=${OutPath}/sub-${subjid}/ses-${session}/sub-${subjid}_ses-${session}_hemi-${hemi}_BOLD_map-midthickness.5k.func.gii
SmoothingKernel=1.5
MetricOut=${OutPath}/sub-${subjid}/ses-${session}/sub-${subjid}_ses-${session}_hemi-${hemi}_BOLD_map-midthickness.5k.s15.func.gii
${WorkbenchBinary} -metric-smoothing ${MidthicknessSmoothed} ${MetricIn} ${SmoothingKernel} ${MetricOut} 
done

###### Connect left and right
echo "BOLD_map L and R"
BOLDLeft=${OutPath}/sub-${subjid}/ses-${session}/sub-${subjid}_ses-${session}_hemi-left_BOLD_map-midthickness.5k.s15.func.gii
BOLDRight=${OutPath}/sub-${subjid}/ses-${session}/sub-${subjid}_ses-${session}_hemi-right_BOLD_map-midthickness.5k.s15.func.gii
BOLDLR=${OutPath}/1.5smooth/sub-${subjid}/ses-${session}/sub-${subjid}_ses-${session}_hemi-LR_BOLD_map-midthickness.hemi_5k.dtseries.nii
${WorkbenchBinary} -cifti-create-dense-timeseries ${BOLDLR} -left-metric ${BOLDLeft} -right-metric ${BOLDRight}

###### Pearson's correlation to creat functional connectome
BOLDLR=${OutPath}/1.5smooth/sub-${subjid}/ses-${session}/sub-${subjid}_ses-${session}_hemi-LR_BOLD_map-midthickness.hemi_5k.dtseries.nii
CorrLR=${OutPath}/1.5smooth/sub-${subjid}/ses-${session}/sub-${subjid}_ses-${session}_hemi-LR_BOLD_correlation.hemi_5k.dconn.nii
${WorkbenchBinary} -cifti-correlation ${BOLDLR} ${CorrLR}


