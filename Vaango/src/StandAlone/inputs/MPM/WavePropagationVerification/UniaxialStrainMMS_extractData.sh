#------------------------------------------------------------------------------
../selectpart -mat 0 -box -0.15 0 0 0.15 0.0005 0.0005 -timesteplow 0 -timestephigh 0 UniaxialStrain_MMS.uda.000 > UniaxialStrainMMS_part.dat
../extractV -m 0 -p UniaxialStrainMMS_part.dat -uda UniaxialStrain_MMS.uda.000 -o UniaxialStrainMMS_pvel.dat -timefiles
../extractPvec -m 0 -partvar p.externalforce -p UniaxialStrainMMS_part.dat -uda UniaxialStrain_MMS.uda.000 -o UniaxialStrainMMS_pfext.dat -timefiles
../extractPvec -m 0 -partvar p.displacement -p UniaxialStrainMMS_part.dat -uda UniaxialStrain_MMS.uda.000 -o UniaxialStrainMMS_pdisp.dat -timefiles
../extractPmat -m 0 -partvar p.deformationGradient -p UniaxialStrainMMS_part.dat -uda UniaxialStrain_MMS.uda.000 -o UniaxialStrainMMS_pF.dat -timefiles
../extractPmat -m 0 -partvar p.stress -p UniaxialStrainMMS_part.dat -uda UniaxialStrain_MMS.uda.000 -o UniaxialStrainMMS_pstress.dat -timefiles

#------------------------------------------------------------------------------
#../selectpart -mat 0 -box 0 0 0 0.1 0.0005 0.0005 -timesteplow 0 -timestephigh 0 UniaxialStrain_MMS_pvelBC.uda.000 > UniaxialStrainMMS_part.dat
#../extractV -m 0 -p UniaxialStrainMMS_part.dat -uda UniaxialStrain_MMS_pvelBC.uda.000 -o UniaxialStrainMMS_pvel.dat -timefiles
#../extractPvec -m 0 -partvar p.externalforce -p UniaxialStrainMMS_part.dat -uda UniaxialStrain_MMS_pvelBC.uda.000 -o UniaxialStrainMMS_pfext.dat -timefiles
#../extractPvec -m 0 -partvar p.displacement -p UniaxialStrainMMS_part.dat -uda UniaxialStrain_MMS_pvelBC.uda.000 -o UniaxialStrainMMS_pdisp.dat -timefiles
#../extractPmat -m 0 -partvar p.deformationGradient -p UniaxialStrainMMS_part.dat -uda UniaxialStrain_MMS_pvelBC.uda.000 -o UniaxialStrainMMS_pF.dat -timefiles
#../extractPmat -m 0 -partvar p.stress -p UniaxialStrainMMS_part.dat -uda UniaxialStrain_MMS_pvelBC.uda.000 -o UniaxialStrainMMS_pstress.dat -timefiles
#------------------------------------------------------------------------------
#../selectpart -mat 0 -box 0 0 0 0.1 0.0005 0.0005 -timesteplow 0 -timestephigh 0 UniaxialStrainHomoLin_MMS.uda.000 > UniaxialStrainMMSHomoLin_part.dat
#../extractV -m 0 -p UniaxialStrainMMSHomoLin_part.dat -uda UniaxialStrainHomoLin_MMS.uda.000 -o UniaxialStrainMMSHomoLin_pvel.dat -timefiles
#../extractPvec -m 0 -partvar p.externalforce -p UniaxialStrainMMSHomoLin_part.dat -uda UniaxialStrainHomoLin_MMS.uda.000 -o UniaxialStrainMMSHomoLin_pfext.dat -timefiles
#../extractPvec -m 0 -partvar p.displacement -p UniaxialStrainMMSHomoLin_part.dat -uda UniaxialStrainHomoLin_MMS.uda.000 -o UniaxialStrainMMSHomoLin_pdisp.dat -timefiles
#../extractPmat -m 0 -partvar p.deformationGradient -p UniaxialStrainMMSHomoLin_part.dat -uda UniaxialStrainHomoLin_MMS.uda.000 -o UniaxialStrainMMSHomoLin_pF.dat -timefiles
#../extractPmat -m 0 -partvar p.stress -p UniaxialStrainMMSHomoLin_part.dat -uda UniaxialStrainHomoLin_MMS.uda.000 -o UniaxialStrainMMSHomoLin_pstress.dat -timefiles
#------------------------------------------------------------------------------
#../selectpart -mat 0 -box 0 0 0 0.1 0.0005 0.0005 -timesteplow 0 -timestephigh 0 UniaxialStrainHomoQuad_MMS.uda.000 > UniaxialStrainMMSHomoQuad_part.dat
#../extractV -m 0 -p UniaxialStrainMMSHomoQuad_part.dat -uda UniaxialStrainHomoQuad_MMS.uda.000 -o UniaxialStrainMMSHomoQuad_pvel.dat -timefiles
#../extractPvec -m 0 -partvar p.externalforce -p UniaxialStrainMMSHomoQuad_part.dat -uda UniaxialStrainHomoQuad_MMS.uda.000 -o UniaxialStrainMMSHomoQuad_pfext.dat -timefiles
#../extractPvec -m 0 -partvar p.displacement -p UniaxialStrainMMSHomoQuad_part.dat -uda UniaxialStrainHomoQuad_MMS.uda.000 -o UniaxialStrainMMSHomoQuad_pdisp.dat -timefiles
#../extractPmat -m 0 -partvar p.deformationGradient -p UniaxialStrainMMSHomoQuad_part.dat -uda UniaxialStrainHomoQuad_MMS.uda.000 -o UniaxialStrainMMSHomoQuad_pF.dat -timefiles
#../extractPmat -m 0 -partvar p.stress -p UniaxialStrainMMSHomoQuad_part.dat -uda UniaxialStrainHomoQuad_MMS.uda.000 -o UniaxialStrainMMSHomoQuad_pstress.dat -timefiles
#------------------------------------------------------------------------------
