1. Before 2015, Sofia downloaded and process all the RIL sequence.
  1 Raw sequence and links to RIL libraries were saved in FC133_RIL_39 like directories
FC133_RIL_39
FC153_RIL_13_24
FC153_RIL_1_12
FC153_RIL_25-26_29-31_33-38
FC153_RIL_42-49_51-53
FC193_HEG4_NB_25_99_40_60
FC193_RIL_12_13_16_43_45_47
FC197_RIL_115_87_93_94_96_134_136_137_139_141_151_152
FC197_RIL_120_121_240_241_243_256_178_190_191_204_210_226
FC197_RIL_123_166_231_256_268_277
FC197_RIL_50_62_55_57_56_58_63_66_68_69_70
FC197_RIL_54_59_64_65_78_80_81_102_132_213_221_222
FC197_RIL_61_74_67_73_74_77_76_112_113_116_118_119
FC197_RIL_75_140_179_234_259_A160
FC205_RIL_92_98_104_131_169_205_207_209_217_219_227_228
FC251_RIL_97_100_135_145_150_162_167_168_187_190_206_208_211_216-219_223_243_250_252_253_273_278
FC271_RIL_78_97_167_177_228_230_232256_265
FC279_RIL242
cornell_RILs_0813
cornell_RILs_1213
  2 Format of file name
Orinial file name:
flowcell205_lane2_pair1_ACAGTG.fastq
flowcell205_lane2_pair1_ACAGTG.fastq
Converted file name:
RIL104_0_ATGTCA_FC205L2_p1.fq
RIL104_0_ATGTCA_FC205L2_p1.fq
  3 Raw sequence directory for analysis was "Illumina"
  4 Genotype results were in "genotypes/MSU_r7.corrected"
  5 Processing scripts

2. After 2015, Jinfeng and Lulu corrected some errors in RIL ID and created the final version of RIL sequence and genotyping results.
  1 Process undetermined sequence in "FC_undetermined"
Found flowcell251_lane4_pair1_GTGAAA.fastq and flowcell251_lane4_pair2_GTGAAA.fastq in FC251, which should RIL167
  2 Correct RIL ID as verified by Jinfeng and Lulu. Put sequence in Illumina_correct and genotype in genotypes_correct
    1 We put all raw sequence in new directory FC_raw.
Merge genotype_sample_illuminaID.txt file from all UCR flowcell together.
Merge cornell.rename.txt file from all cornal flowcell together.
These two file need to manual edited to make sure all the names are corrected, using # to note why and when we made the decision. 
In this dir, we write a script to link the raw fastq to right RIL with correct ID.
We need one file have library which had bad quality/heterozygous we need to archive these: achieve.txt
    2 Do genotype for these switched ID, or changed ID or FC251, then link other unchanged file from Sofia's results.


